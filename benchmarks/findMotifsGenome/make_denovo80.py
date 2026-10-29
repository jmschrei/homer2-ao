"""HOMER v5.1 references for the de novo check compare.py runs (run once).

	python3 make_denovo80.py --cases <dir>

with HOMER_FMG_REF and HOMER_BENCH_GENOME set (fmg.py). Inputs: 80
subsamples of the homer2 denovo input findMotifsGenome.pl writes for the
objective case (100, 200, 300 or 600 targets with 4x background; seeds
101-148 at -len 10/12, 201-232 at -len 14/16, -S 10 or 25), made by a script
that is not part of this directory. <dir> holds tri.log, one line
"seed <n> len <L> S <S> ..." per subsample, and cases/s<seed>/seq.tsv and
group.adj. They are copied gzipped into <reference>/denovo80/inputs/. Each is
run with HOMER v5.1's homer2 (<reference>/v51/homer) at -p 1 -T 12 -mis 2
-cache 500 in the harness environment (fmg.py homer_env); the exit status and
sha256 of motifs.out and stdout go to <reference>/denovo80/expected.json, the
motif files beside it.
"""

import os
import re
import sys
import gzip
import json
import shutil
import hashlib
import tempfile
import subprocess
import argparse
import concurrent.futures

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fmg

CASES_DIR = None  # --cases
OUT = os.path.join(fmg.REF, 'denovo80')
WORKERS = 16


def cases():
	rows = []
	with open(os.path.join(CASES_DIR, 'tri.log')) as handle:
		for line in handle:
			m = re.match(r'seed (\d+) len (\d+) S (\d+) ', line)
			if m:
				rows.append({'name': 's{}'.format(m.group(1)), 'len': int(m.group(2)),
					'S': int(m.group(3))})
	assert len(rows) == 80, len(rows)
	return rows


def run(home, inputs, case, p, extra, scratch):
	d = tempfile.mkdtemp(dir=scratch, prefix='{}.p{}.'.format(case['name'], p))
	argv = [os.path.join(home, 'bin', 'homer2'), 'denovo', '-s',
		os.path.join(inputs, case['name'] + '.seq.tsv'), '-g',
		os.path.join(inputs, case['name'] + '.group.adj'), '-S', str(case['S']),
		'-mis', '2', '-cache', '500', '-p', str(p), '-T', '12', '-len',
		str(case['len']), '-o', 'motifs.out'] + list(extra)
	with open(os.path.join(d, '.stdout'), 'wb') as out:
		code = subprocess.call(argv, cwd=d, stdout=out, stderr=subprocess.DEVNULL,
			stdin=subprocess.DEVNULL, env=fmg.homer_env(home, p, d))
	got = {'exit': code}
	for key, name in (('motifs', 'motifs.out'), ('stdout', '.stdout')):
		path = os.path.join(d, name)
		got[key] = fmg.digest(open(path, 'rb').read()) if os.path.exists(path) else None
	return got, d


def main():
	if os.path.exists(OUT):
		sys.exit('{} exists: remove it first'.format(OUT))
	rows = cases()
	os.makedirs(os.path.join(OUT, 'inputs'))
	os.makedirs(os.path.join(OUT, 'motifs'))
	for case in rows:
		for kind in ('seq.tsv', 'group.adj'):
			with open(os.path.join(CASES_DIR, 'cases', case['name'], kind), 'rb') as src, \
					gzip.open(os.path.join(OUT, 'inputs', '{}.{}.gz'.format(case['name'],
						kind)), 'wb') as dst:
				shutil.copyfileobj(src, dst)
	scratch = tempfile.mkdtemp(dir=fmg.SCRATCH, prefix='denovo80-ref-')
	try:
		inputs = os.path.join(scratch, 'in')
		os.makedirs(inputs)
		for case in rows:
			for kind in ('seq.tsv', 'group.adj'):
				shutil.copy(os.path.join(CASES_DIR, 'cases', case['name'], kind),
					os.path.join(inputs, '{}.{}'.format(case['name'], kind)))
		with concurrent.futures.ThreadPoolExecutor(WORKERS) as pool:
			futures = [pool.submit(run, fmg.V51, inputs, case, 1, (), scratch)
				for case in rows]
			for case, future in zip(rows, futures):
				got, d = future.result()
				case.update(got)
				shutil.copy(os.path.join(d, 'motifs.out'), os.path.join(OUT, 'motifs',
					'{}.L{}.S{}.motifs'.format(case['name'], case['len'], case['S'])))
	finally:
		shutil.rmtree(scratch, ignore_errors=True)
	with open(os.path.join(OUT, 'expected.json'), 'w') as handle:
		json.dump(rows, handle, indent=1)
	bad = [c['name'] for c in rows if c['exit'] != 0]
	print('{} cases, {} non-zero exits {}'.format(len(rows), len(bad), bad))


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--cases', required=True)
	CASES_DIR = os.path.abspath(parser.parse_args().cases)
	main()

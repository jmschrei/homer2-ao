"""Build the HOMER v5.1 references the checks compare against (run once).

	python3 make_refs.py inputs --atac <ENCFF333TAT .narrowPeak.gz> --ctcf <ENCFF901CBP .narrowPeak.gz>
	python3 make_refs.py build <HOMER v5.1 tree>   # built into <reference>/v51/homer
	python3 make_refs.py preparse        # v5.1's cold preparse -> <reference>/preparsed (read-only)
	python3 make_refs.py case <name>     # v5.1 at -p 1 with homer2 denovo -T of the case's -p

with HOMER_FMG_REF and HOMER_BENCH_GENOME set (fmg.py). Every reference is
HOMER v5.1 built from the release. Each case's output directory is kept under
<reference>/cases/<case>/ (motifFindingParameters.txt normalized), with its
manifest.json of sha256 per file, stdout and exit status.
"""

import os
import sys
import gzip
import json
import time
import shutil
import argparse
import hashlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fmg

# (file, source, prefix, first rank, last rank), ranks 1-based by signalValue
# (column 7) descending, ties in file order; source is --atac or --ctcf.
PEAK_SETS = [
	('peaks10k.bed', 'atac', 'atac', 1, 10000),
	('peaks2k.bed', 'atac', 'atac', 10001, 12000),
	('ctcf10k.bed', 'ctcf', 'ctcf', 1, 10000),
	('peaks300.bed', 'atac', 'atac', 12001, 12300),
	('peaks500.bed', 'atac', 'atac', 10001, 10500),
]


def ranked(path):
	rows = []
	with gzip.open(path, 'rt') as handle:
		for line in handle:
			cols = line.rstrip('\n').split('\t')
			rows.append((cols[0], int(cols[1]), float(cols[6]), int(cols[9])))
	return sorted(rows, key=lambda row: -row[2])


def inputs(sources):
	os.makedirs(fmg.INPUTS, exist_ok=True)
	cache = {}
	for name, source, prefix, first, last in PEAK_SETS:
		source = sources[source]
		if source not in cache:
			cache[source] = ranked(source)
		rows = cache[source][first - 1:last]
		assert len(rows) == last - first + 1, (name, len(rows))
		lines = []
		for rank, (chrom, start, signal, summit) in enumerate(rows, first):
			left = max(0, start + summit - 100)
			lines.append('{}\t{}\t{}\t{}{}\t{:g}\t+\n'.format(chrom, left, left + 200,
				prefix, rank, signal))
		path = os.path.join(fmg.INPUTS, name)
		with open(path, 'w') as handle:
			handle.writelines(lines)
		with open(path, 'rb') as handle:
			digest = hashlib.sha256(handle.read()).hexdigest()
		print('{} {} peaks, ranks {}-{}, sha256 {}'.format(name, len(lines), first,
			last, digest))


def build(src):
	home = fmg.V51
	if os.path.exists(home):
		sys.exit('{} exists: remove it first'.format(home))
	fmg.copy_tree(src, home)
	seconds = fmg.configure(home, os.path.join(os.path.dirname(home), 'build.log'))
	print('HOMER v5.1 built in {} in {:.1f} s'.format(home, seconds))


def save(case, result, extra):
	folder = os.path.join(fmg.CASES_REF, case)
	if os.path.exists(folder):
		shutil.rmtree(folder)
	os.makedirs(folder)
	meta = {key: result[key] for key in ('code', 'seconds', 'cpu_seconds', 'argv')}
	meta.update(extra)
	meta['date'] = time.strftime('%Y-%m-%d %H:%M:%S')
	meta['steps'] = result.get('steps')
	with open(os.path.join(folder, 'meta.json'), 'w') as handle:
		json.dump(meta, handle, indent=1)
	with open(os.path.join(folder, 'manifest.json'), 'w') as handle:
		json.dump(result['manifest'], handle, indent=1, sort_keys=True)
	return folder


def preparse():
	if os.path.exists(fmg.PREPARSED):
		sys.exit('{} exists: remove it first'.format(fmg.PREPARSED))
	keep = os.path.join(fmg.SCRATCH, 'fmg-ref-preparse-{}'.format(os.getpid()))
	result = fmg.run_preparse(fmg.V51, keep=keep)
	if result['code'] != 0:
		sys.exit('HOMER v5.1 preparse exited {}; see {}'.format(result['code'], keep))
	folder = save('preparse', result, {})
	for name in ('.stdout', '.stderr'):
		shutil.copy(os.path.join(keep, name), os.path.join(folder, name[1:] + '.txt'))
	shutil.copytree(os.path.join(keep, 'preparsed'), fmg.PREPARSED)
	for name in os.listdir(fmg.PREPARSED):
		os.chmod(os.path.join(fmg.PREPARSED, name), 0o444)
	os.chmod(fmg.PREPARSED, 0o555)
	shutil.rmtree(keep)
	print('preparse: {:.1f} s, {} files: {}'.format(result['seconds'],
		len(result['manifest']) - 2, sorted(result['manifest'])))


def case(name):
	spec = fmg.CASES[name]
	keep = os.path.join(fmg.SCRATCH, 'fmg-ref-{}-{}'.format(name, os.getpid()))
	result = fmg.run_findmotifs(fmg.V51, name, p=1, trials=spec['T'], keep=keep,
		label='ref-' + name)
	folder = save(name, result, {'p': 1, 'T': spec['T'], 'case': spec})
	with open(os.path.join(keep, '.steps.tsv')) as handle:
		denovo = [line for line in handle if '\tdenovo ' in line]
	assert denovo and all(line.rstrip('\n').endswith(' -T {}'.format(spec['T']))
		for line in denovo), denovo
	for item in ('.stdout', '.stderr', '.steps.tsv'):
		shutil.copy(os.path.join(keep, item), os.path.join(folder, item[1:]))
	shutil.copytree(os.path.join(keep, 'out'), os.path.join(folder, 'out'))
	params = os.path.join(folder, 'out', 'motifFindingParameters.txt')
	if os.path.exists(params):
		with open(params, 'rb') as handle:
			data = fmg.normalize(params, handle.read())
		with open(params, 'wb') as handle:
			handle.write(data)
	shutil.rmtree(keep)
	print('{}: exit {}, {:.1f} s, {} files, denovo calls {}'.format(name,
		result['code'], result['seconds'], len(result['manifest']) - 2,
		len(denovo)))


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('what', choices=['inputs', 'build', 'preparse', 'case'])
	parser.add_argument('name', nargs='?', help='the case, or for build the HOMER v5.1 tree')
	parser.add_argument('--atac')
	parser.add_argument('--ctcf')
	args = parser.parse_args()
	if args.what == 'inputs':
		inputs({'atac': args.atac, 'ctcf': args.ctcf})
	elif args.what == 'build':
		build(os.path.abspath(args.name))
	elif args.what == 'preparse':
		preparse()
	else:
		case(args.name)

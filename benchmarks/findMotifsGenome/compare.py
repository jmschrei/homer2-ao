"""HOMER's other entry points on small inputs, against HOMER v5.1's outputs.

	python3 compare.py <src> <reference>
	HOMER_FMG_WRITE_REFERENCE=1 python3 compare.py <v5.1 tree> <reference>   # cache v5.1's

with HOMER_FMG_REF and HOMER_BENCH_GENOME set (fmg.py). <src> is a HOMER tree;
it is built in place once (fmg.py build_inplace). <reference> is the
directory the outputs are cached in (<reference dir>/checks/outputs). Each
case runs in a fresh /dev/shm directory with the benchmark's fixed
environment (Perl hash seed 0, srand(1) at every perl start, PATH = the
tree's bin and the system's, LC_ALL=C), on the small inputs make_inputs.py
writes (hg38 chr21+chr22 as a custom genome, ATAC peaks there, FASTA and
seq/group tables, a motif subset, edge-case peak files). Their paths are part
of some outputs, so run with the inputs directory the cache was made with.
Every file the case leaves in its directory, stdout and the exit status are
hashed (sha256); a case differs when any of them does, or when a file is on
one side only. stderr is not compared.

bench.py checks the timed call itself against v5.1 on every run; these cases
cover the other callers of the code a change may touch: findMotifsGenome.pl's
other options and odd inputs, findMotifs.pl, preparseGenome.pl, every homer2
subcommand, homerTools, mergePeaks, annotatePeaks.pl, scanMotifGenomeWide.pl
and compareMotifs.pl. Cases with homer2 denovo run at -p 1 (v5.1's -p > 1
denovo race); the others at -p 1 or 4 where v5.1 is deterministic.

Exit 0 when every case matches, 1 when any differs, 2 when <reference> holds
no outputs. One line per case; the last line is the summary.
"""

import os
import sys
import gzip
import json
import time
import shutil
import tempfile
import subprocess
import concurrent.futures

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fmg

I = fmg.CHECK_INPUTS
MINI = os.path.join(I, 'mini.fa')
WORKERS = 8
TIMEOUT = 900
SECONDS = {}

# The de novo check: homer2 denovo on 80 subsamples of the objective's de
# novo input at -p 2 -T 12, and 4 of them that are sensitive to how the
# oligo screen re-checks oligos after each mask, at -p 1, 4 and 16 with -T 12
# -Tcpus 1, each against HOMER v5.1 at -p 1 -T 12 (<reference dir>/denovo80,
# make_denovo80.py): exit status, the motif file and stdout. Skipped in write
# mode (its reference is v5.1's own run) and on a tree whose homer2 sources are
# v5.1's.
DENOVO = os.path.join(fmg.REF, 'denovo80')
RESCREEN = ('s106', 's133', 's148', 's219')
HOMER2_SOURCES = ('cpp/homer2.cpp', 'cpp/Motif2.cpp', 'cpp/Motif2.h',
	'cpp/Hashtable.cpp', 'cpp/Hashtable.h', 'cpp/statistics.cpp', 'cpp/statistics.h',
	'cpp/Makefile')
DENOVO_WORKERS = 12


def homer2_is_v51(src):
	for rel in HOMER2_SOURCES:
		try:
			with open(os.path.join(src, rel), 'rb') as a, open(os.path.join(fmg.V51,
					rel), 'rb') as b:
				if a.read() != b.read():
					return False
		except OSError:
			return False
	return True


def denovo_run(home, inputs, case, p, extra):
	d = tempfile.mkdtemp(dir=inputs, prefix='{}.p{}.'.format(case['name'], p))
	argv = [os.path.join(home, 'bin', 'homer2'), 'denovo', '-s',
		os.path.join(inputs, case['name'] + '.seq.tsv'), '-g',
		os.path.join(inputs, case['name'] + '.group.adj'), '-S', str(case['S']),
		'-mis', '2', '-cache', '500', '-p', str(p), '-T', '12', '-len',
		str(case['len']), '-o', 'motifs.out'] + list(extra)
	with open(os.path.join(d, '.stdout'), 'wb') as out:
		try:
			code = subprocess.call(argv, cwd=d, stdout=out, stderr=subprocess.DEVNULL,
				stdin=subprocess.DEVNULL, env=fmg.homer_env(home, p, d), timeout=TIMEOUT)
		except subprocess.TimeoutExpired:
			code = 'timeout'
	got = [code]
	for name in ('motifs.out', '.stdout'):
		path = os.path.join(d, name)
		got.append(fmg.digest(open(path, 'rb').read()) if os.path.exists(path) else None)
	shutil.rmtree(d, ignore_errors=True)
	return got


def denovo_check(home):
	"""{group: (runs, [differing run names], seconds)}."""

	with open(os.path.join(DENOVO, 'expected.json')) as handle:
		cases = json.load(handle)
	byname = {case['name']: case for case in cases}
	jobs = [('denovo80_p2', case, 2, ()) for case in cases] + [
		('rescreen_v51', byname[name], p, ('-Tcpus', '1')) for name in RESCREEN
		for p in (1, 4, 16)]
	start = time.time()
	inputs = tempfile.mkdtemp(dir=fmg.SCRATCH, prefix='fmgcmp-denovo-')
	try:
		for case in cases:
			for kind in ('seq.tsv', 'group.adj'):
				with gzip.open(os.path.join(DENOVO, 'inputs', '{}.{}.gz'.format(
						case['name'], kind)), 'rb') as src, open(os.path.join(inputs,
						'{}.{}'.format(case['name'], kind)), 'wb') as dst:
					shutil.copyfileobj(src, dst)
		with concurrent.futures.ThreadPoolExecutor(DENOVO_WORKERS) as pool:
			futures = [pool.submit(denovo_run, home, inputs, case, p, extra)
				for _, case, p, extra in jobs]
			out = {}
			for (group, case, p, _), future in zip(jobs, futures):
				runs, bad, _ = out.setdefault(group, (0, [], 0))
				if future.result() != [case['exit'], case['motifs'], case['stdout']]:
					bad.append('{} -len {} -S {} -p {}'.format(case['name'], case['len'],
						case['S'], p))
				out[group] = (runs + 1, bad, 0)
	finally:
		shutil.rmtree(inputs, ignore_errors=True)
	seconds = time.time() - start
	return {group: (runs, bad, seconds) for group, (runs, bad, _) in out.items()}


def inp(name):
	return os.path.join(I, name)


def fmg_call(peaks, extra, p=1, size='200', N='2000'):
	return ['findMotifsGenome.pl', inp(peaks), MINI, 'out', '-size', size, '-N', N,
		'-p', str(p), '-preparsedDir', 'pre'] + extra


# (name, -p for OMP_NUM_THREADS, argv, stage the size-200 preparse of mini.fa)
CASES = [
	('fmg_default_p1', 1, fmg_call('peaks.bed', ['-S', '5', '-len', '8,10']), True),
	('fmg_known_p4', 4, fmg_call('peaks.bed', ['-nomotif'], p=4), True),
	('fmg_size_given', 1, fmg_call('peaks2.bed', ['-S', '3', '-len', '8',
		'-noknown'], size='given'), True),
	('fmg_size100_cold', 1, fmg_call('peaks.bed', ['-nomotif'], size='100'), True),
	('fmg_mask_cold', 1, fmg_call('peaks.bed', ['-nomotif', '-mask']), True),
	('fmg_custom_bg', 1, fmg_call('peaks.bed', ['-bg', inp('bg.bed'), '-S', '3',
		'-len', '8']), True),
	('fmg_find', 1, fmg_call('peaks.bed', ['-find', inp('ctcf.motif')]), True),
	('fmg_mknown_hypergeo', 4, fmg_call('peaks.bed', ['-mknown',
		inp('known_sub.motifs'), '-nomotif', '-h', '-nlen', '2'], p=4), True),
	('fmg_options_mix', 1, fmg_call('peaks2.bed', ['-norevopp', '-len', '9', '-S',
		'3', '-mis', '1', '-nlen', '0', '-mknown', inp('known_sub.motifs'), '-mcheck',
		inp('known_sub.motifs'), '-dumpFasta', '-basic', '-nofacts']), True),
	('fmg_mset_bits_p2', 2, fmg_call('peaks2.bed', ['-mset', 'vertebrates',
		'-bits', '-nomotif'], p=2), True),
	('fmg_empty', 1, fmg_call('empty.bed', ['-S', '2', '-len', '8'], N='100'), True),
	('fmg_one_peak', 1, fmg_call('one.bed', ['-S', '2', '-len', '8'], N='100'), True),
	('fmg_edge_peaks', 1, fmg_call('edge.bed', ['-nomotif'], N='500'), True),
	('findMotifs_fasta', 1, ['findMotifs.pl', inp('t200.fa'), 'fasta', 'fout',
		'-fasta', inp('b1k.fa'), '-len', '8', '-S', '3', '-p', '1'], False),
	('homer2_denovo_fasta', 1, ['homer2', 'denovo', '-i', inp('t200.fa'), '-b',
		inp('b1k.fa'), '-len', '8', '-S', '3', '-p', '1', '-o', 'denovo.motifs'], False),
	('homer2_denovo_sg', 1, ['homer2', 'denovo', '-s', inp('seq.tsv'), '-g',
		inp('group.tsv'), '-len', '10', '-S', '4', '-p', '1', '-stat', 'hypergeo',
		'-o', 'denovo.motifs'], False),
	('homer2_known_fasta_p4', 4, ['homer2', 'known', '-i', inp('t200.fa'), '-b',
		inp('b1k.fa'), '-m', inp('known_sub.motifs'), '-p', '4', '-mout',
		'known.motifs'], False),
	('homer2_known_sg', 1, ['homer2', 'known', '-s', inp('seq.tsv'), '-g',
		inp('group.tsv'), '-m', inp('known_sub.motifs'), '-p', '1', '-stat',
		'hypergeo', '-nlen', '2', '-nout', 'norm.txt', '-o', 'known.txt'], False),
	('homer2_find_p1', 1, ['homer2', 'find', '-i', inp('t200.fa'), '-m',
		inp('known_sub.motifs'), '-p', '1'], False),
	('homer2_find_mscore', 1, ['homer2', 'find', '-i', inp('t200.fa'), '-m',
		inp('known_sub.motifs'), '-mscore', '-p', '1'], False),
	('homer2_mask', 1, ['homer2', 'mask', '-s', inp('seq.tsv'), '-m',
		inp('known_sub.motifs'), '-p', '1'], False),
	('homer2_norm', 1, ['homer2', 'norm', '-s', inp('seq.tsv'), '-g',
		inp('group.tsv'), '-nlen', '3', '-nout', 'norm.txt'], False),
	('homer2_bg_seeded', 1, ['homer2', 'bg', '-p', inp('peaks.bed'), '-g', MINI,
		'-N', '1000', '-NN', '1000000', '-seed', '7', '-o', 'bgout'], False),
	('homerTools_extract', 1, ['homerTools', 'extract', inp('edge.bed'), MINI],
		False),
	('homerTools_extract_fa', 1, ['homerTools', 'extract', inp('peaks.bed'), MINI,
		'-fa'], False),
	('homerTools_freq', 1, ['homerTools', 'freq', inp('t200.fa'), '-gc', 'gc.txt'],
		False),
	('mergePeaks_cobound', 1, ['mergePeaks', inp('peaks.bed'), inp('peaks2.bed'),
		inp('bg.bed'), '-d', '200', '-cobound', '1', '-prefix', 'cb', '-venn',
		'venn.txt'], False),
	('mergePeaks_given', 1, ['mergePeaks', inp('peaks.bed'), inp('peaks2.bed'),
		'-d', 'given'], False),
	('annotatePeaks_motif', 1, ['annotatePeaks.pl', inp('peaks.bed'), MINI, '-m',
		inp('ctcf.motif'), '-size', '200', '-noann', '-nogene'], False),
	('scanMotifGenomeWide', 1, ['scanMotifGenomeWide.pl', inp('ctcf.motif'), MINI,
		'-bed', '-p', '1'], False),
	('compareMotifs_direct', 1, ['compareMotifs.pl', inp('known_sub.motifs'), 'cmp/',
		'-known', inp('ctcf.motif'), '-cpu', '1'], False),
]


def preparse_case():
	return ('preparse_mini', 1, ['preparseGenome.pl', MINI, '-size', '200',
		'-preparsedDir', 'pre'], False)


def run_case(home, case, staged):
	name, threads, argv, stage = case
	scratch = tempfile.mkdtemp(dir=fmg.SCRATCH, prefix='fmgcmp-{}-'.format(name))
	try:
		if stage:
			shutil.copytree(staged, os.path.join(scratch, 'pre'))
		env = fmg.homer_env(home, threads, scratch)
		start = time.time()
		with open(os.path.join(scratch, '.stdout'), 'wb') as out, open(
				os.path.join(scratch, '.stderr'), 'wb') as err:
			try:
				code = subprocess.run(argv, cwd=scratch, env=env, stdout=out,
					stderr=err, stdin=subprocess.DEVNULL, timeout=TIMEOUT).returncode
			except subprocess.TimeoutExpired:
				code = 'timeout'
		SECONDS[name] = time.time() - start
		files = {}
		for root, dirs, names in os.walk(scratch):
			dirs[:] = [d for d in dirs if not d.startswith('.')]
			for item in names:
				if item.startswith('.'):
					continue
				path = os.path.join(root, item)
				rel = os.path.relpath(path, scratch)
				files[rel] = fmg.file_digest(path, rel)
		with open(os.path.join(scratch, '.stdout'), 'rb') as handle:
			files['<stdout>'] = fmg.digest(handle.read())
		files['<exit>'] = str(code)
		keep = None
		if name == 'preparse_mini':
			keep = tempfile.mkdtemp(dir=fmg.SCRATCH, prefix='fmgcmp-pre-')
			shutil.rmtree(keep)
			shutil.copytree(os.path.join(scratch, 'pre'), keep)
		with open(os.path.join(scratch, '.stderr'), 'rb') as handle:
			tail = handle.read()[-600:].decode(errors='replace')
		return files, keep, tail
	finally:
		shutil.rmtree(scratch, ignore_errors=True)


def main():
	src = os.path.realpath(sys.argv[1])
	reference = os.path.realpath(sys.argv[2])
	write = os.environ.get('HOMER_FMG_WRITE_REFERENCE') == '1'
	if not write and not os.path.exists(os.path.join(reference, 'preparse_mini.json')):
		print('no outputs in {}'.format(reference))
		sys.exit(2)
	start = time.time()
	fmg.build_inplace(src)

	results = {}
	denovo, skipped = None, None
	if write:
		skipped = 'not cached: compared with v5.1\'s references in {}'.format(DENOVO)
	elif homer2_is_v51(src):
		skipped = 'skipped: homer2\'s sources are v5.1\'s, the reference\'s own code'
	side = concurrent.futures.ThreadPoolExecutor(1)
	if skipped is None:
		denovo = side.submit(denovo_check, src)
	first = preparse_case()
	results[first[0]], staged, _ = run_case(src, first, None)
	try:
		with concurrent.futures.ThreadPoolExecutor(WORKERS) as pool:
			futures = {case[0]: pool.submit(run_case, src, case, staged)
				for case in CASES}
			tails = {}
			for name, future in futures.items():
				results[name], _, tails[name] = future.result()
	finally:
		shutil.rmtree(staged, ignore_errors=True)
	denovo = denovo.result() if denovo is not None else None
	side.shutdown()

	if write:
		os.makedirs(reference, exist_ok=True)
		for name, files in results.items():
			with open(os.path.join(reference, name + '.json'), 'w') as handle:
				json.dump(files, handle, indent=1, sort_keys=True)
			print('{:24s} cached ({} files, exit {}, {:.0f} s)'.format(name,
				len(files) - 2, files['<exit>'], SECONDS.get(name, 0)))
		print('{:24s} {}'.format('denovo80_p2', skipped))
		print('{:24s} {}'.format('rescreen_v51', skipped))
		print('{} inputs cached in {:.0f} s'.format(len(results), time.time() - start))
		return 0

	differ = []
	for name, files in results.items():
		path = os.path.join(reference, name + '.json')
		if not os.path.exists(path):
			print('{:24s} DIFFERS: no cached v5.1 output'.format(name))
			differ.append(name)
			continue
		with open(path) as handle:
			want = json.load(handle)
		keys = fmg.mismatches(files, want)
		if keys:
			differ.append(name)
			print('{:24s} DIFFERS: {} of {} keys: {}'.format(name, len(keys),
				len(set(files) | set(want)), ', '.join(keys[:6])))
			if '<exit>' in keys and name in tails:
				print('    stderr tail: {}'.format(' | '.join(tails[name].strip()
					.splitlines()[-3:])))
		else:
			print('{:24s} identical ({} files, exit {}, {:.0f} s)'.format(name,
				len(files) - 2, files['<exit>'], SECONDS.get(name, 0)))
	total = len(results)
	for group in ('denovo80_p2', 'rescreen_v51'):
		if denovo is None:
			print('{:24s} {}'.format(group, skipped))
			continue
		total += 1
		runs, bad, seconds = denovo[group]
		if bad:
			differ.append(group)
			print('{:24s} DIFFERS from v5.1 -p 1 -T 12: {} of {} runs: {}'.format(group,
				len(bad), runs, ', '.join(bad)))
		else:
			print('{:24s} identical to v5.1 -p 1 -T 12 ({} runs; de novo check {:.0f} s)'
				.format(group, runs, seconds))
	if differ:
		print('{} of {} inputs differ: {}'.format(len(differ), total, ', '.join(differ)))
		return 1
	print('{} inputs identical ({:.0f} s)'.format(total, time.time() - start))
	return 0


if __name__ == '__main__':
	sys.exit(main())

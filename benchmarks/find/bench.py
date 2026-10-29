"""One run of the homer2 find benchmark: build, warm up, time, check.

Objective: wall seconds of one `homer2 find -p 4` call on all 269,800 K562
ATAC-seq peaks of ENCODE ENCFF333TAT (200 bp around each summit, hg38) x the
472 known vertebrate motifs, after one identical warmup call. Every homer2
call runs in its own fresh directory on /dev/shm.

Every reference comes from HOMER v5.1's homer2 (make_refs.py,
make_small_refs.py), so a missing reference is an error. The checks, each a
count of mismatches (0 when every output matches):

- n_find_mismatch_vs_v51, n_mscore_mismatch_vs_v51,
  n_mask_mismatch_vs_v51: the sorted lines of each output against v5.1's,
  summed over every case of one mode;
- n_order_mismatch_vs_v51_p1: every output, unsorted, against v5.1's at
  -p 1 for the same arguments, whatever -p the case runs at.

For big270k, output byte-identical to v5.1 -p 1 implies the sorted check, so
the sorted comparison (13 M lines) runs only when the bytes differ.

	HOMER_FIND_REF=<reference> python3 bench.py --tree <HOMER tree> --seed 0 \\
		--out <dir> --cache <dir>

homer2 is built from <HOMER tree>/cpp with this directory's Makefile.
HOMER_FIND_INPUTS, if set, is the directory holding inputs/ (default: the
reference).
"""

import os
import re
import sys
import json
import time
import fcntl
import gzip
import hashlib
import argparse
import tempfile
import subprocess
import collections

import numpy

parser = argparse.ArgumentParser()
parser.add_argument('--tree', required=True)
parser.add_argument('--seed', type=int, required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--cache', required=True)
args = parser.parse_args()

REF = os.environ['HOMER_FIND_REF']
INPUTS = os.path.join(os.environ.get('HOMER_FIND_INPUTS', REF), 'inputs')
OUTPUTS = os.path.join(REF, 'outputs')
ORDER = os.path.join(REF, 'order_p1')
SRC = os.path.abspath(args.tree)
MAKEFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Makefile')

os.makedirs(args.out, exist_ok=True)


def log(*items):
	print(*items, flush=True)


def inp(name):
	return os.path.join(INPUTS, name)


# --- the portable-build rule, checked mechanically -----------------------------

NONPORTABLE = re.compile(r'-march=native|-mtune=native|-march=(?!x86-64\b)\S+'
	r'|-m(?:avx\w*|sse4\w*|fma\w*|bmi\w*|popcnt|f16c|lzcnt)\b|-ffast-math|-Ofast'
	r'|-funsafe-math\S*|-ffinite-math\S*|-fno-signed-zeros|-fassociative-math'
	r'|-freciprocal-math|-fno-trapping-math')
PRAGMA = re.compile(r'#\s*pragma\s+GCC\s+(?:target|optimize)\b.*'
	r'|__attribute__\s*\(\(\s*(?:target|optimize)\s*\(')


def nonportable_build():
	"""Count build or source lines that make the binary CPU-specific or relax
	floating-point semantics. target_clones (runtime dispatch) is allowed."""

	hits = []
	# What the build reads: this directory's Makefile and the top level of
	# <tree>/cpp. HOMER's own Makefiles build nothing here.
	for root, files in [(SRC, [MAKEFILE]), (os.path.join(SRC, 'cpp'),
			sorted(os.listdir(os.path.join(SRC, 'cpp'))))]:
		for name in files:
			path = os.path.join(root, name)
			is_build = path == MAKEFILE
			is_source = name.endswith(('.cpp', '.h', '.hpp', '.cc'))
			if not (is_build or is_source):
				continue
			with open(path, errors='replace') as handle:
				for number, line in enumerate(handle, 1):
					if (is_build and NONPORTABLE.search(line)) or (
							is_source and PRAGMA.search(line)):
						hits.append('{}:{}: {}'.format(os.path.relpath(path, SRC),
							number, line.strip()))
	for hit in hits:
		log('nonportable:', hit)
	return len(hits)


# --- build ----------------------------------------------------------------------

def build():
	"""Build homer2 once per commit into {cache}; later runs reuse it."""

	os.makedirs(args.cache, exist_ok=True)
	binary = os.path.join(args.cache, 'build', 'homer2')
	with open(os.path.join(args.cache, '.build.lock'), 'w') as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		if not os.path.exists(binary):
			start = time.perf_counter()
			result = subprocess.run(['make', '-f', MAKEFILE, '-C', SRC, '-j', '4',
				'OUT=' + os.path.join(args.cache, 'build')],
				capture_output=True, text=True)
			sys.stdout.write(result.stdout[-4000:])
			sys.stdout.write(result.stderr[-4000:])
			if result.returncode != 0 or not os.path.exists(binary):
				sys.exit('build failed')
			log('built {} in {:.1f} s'.format(binary, time.perf_counter() - start))
	return binary


# --- one homer2 call ------------------------------------------------------------

def homer(binary, mode, fasta, motifs, threads, extra=(), peak_rss=False):
	"""Run one homer2 call in a fresh /dev/shm directory.

	Returns (wall seconds, peak RSS in MB or None, user+sys CPU seconds,
	output bytes). With peak_rss the call runs under /usr/bin/time -f %M
	(os.wait4's ru_maxrss counts the forking Python parent), so that call's
	seconds include time's own start-up and are not used for timing.
	"""

	with tempfile.TemporaryDirectory(dir='/dev/shm', prefix='homerbench-') as tmp:
		output = os.path.join(tmp, 'out.txt')
		command = [binary, mode, '-i', fasta, '-m', motifs, '-p', str(threads),
			'-o', output] + list(extra)
		if peak_rss:
			command = ['/usr/bin/time', '-f', '%M', '-o', os.path.join(tmp, 'rss.txt')] + command

		with open(os.path.join(tmp, 'stderr.txt'), 'wb') as err:
			start = time.perf_counter()
			process = subprocess.Popen(command, cwd=tmp, stdin=subprocess.DEVNULL,
				stdout=subprocess.DEVNULL, stderr=err)
			_, status, usage = os.wait4(process.pid, 0)
			seconds = time.perf_counter() - start
			process.returncode = os.waitstatus_to_exitcode(status)

		if process.returncode != 0:
			with open(os.path.join(tmp, 'stderr.txt'), 'rb') as handle:
				sys.stdout.write(handle.read()[-3000:].decode(errors='replace'))
			sys.exit('homer2 exited {}: {}'.format(process.returncode,
				' '.join(command)))

		with open(output, 'rb') as handle:
			data = handle.read()
		rss = None
		if peak_rss:
			with open(os.path.join(tmp, 'rss.txt')) as handle:
				rss = int(handle.read().split()[-1]) / 1024

	return seconds, rss, usage.ru_utime + usage.ru_stime, data


# --- comparison against HOMER v5.1 ----------------------------------------------

def compare(case, data):
	"""Lines of `data`, sorted, against HOMER v5.1's for `case`.

	Returns the size of the multiset symmetric difference: 0 when the sorted
	outputs are byte-identical.
	"""

	lines = data.split(b'\n')
	if lines and lines[-1] == b'':
		lines.pop()
	lines.sort()
	blob = b'\n'.join(lines) + b'\n'
	digest = hashlib.sha256(blob).hexdigest()

	path = os.path.join(OUTPUTS, case + '.sorted.txt.gz')
	if not os.path.exists(path):
		sys.exit('no HOMER v5.1 reference for case {} at {}'.format(case, path))

	with open(path + '.sha256') as handle:
		if handle.read().strip() == digest:
			return 0

	with gzip.open(path, 'rb') as handle:
		expected = handle.read().split(b'\n')
	if expected and expected[-1] == b'':
		expected.pop()

	got, want = collections.Counter(lines), collections.Counter(expected)
	n = sum(((got - want) + (want - got)).values())
	log('MISMATCH {}: {} lines differ ({} here, {} in v5.1)'.format(
		case, n, len(lines), len(expected)))
	for line in list((got - want).elements())[:3]:
		log('  only here:     ', line.decode(errors='replace'))
	for line in list((want - got).elements())[:3]:
		log('  only v5.1:     ', line.decode(errors='replace'))
	return n


def compare_order(case, data):
	"""`data`, unsorted, against HOMER v5.1's output for `case` at -p 1.

	Returns the number of line positions that differ plus the difference in
	line count: 0 when the outputs are byte-identical.
	"""

	path = os.path.join(ORDER, case + '.txt.gz')
	if not os.path.exists(path):
		sys.exit('no HOMER v5.1 -p 1 order reference for case {} at {}'.format(case, path))
	with open(path + '.sha256') as handle:
		if handle.read().strip() == hashlib.sha256(data).hexdigest():
			return 0

	with gzip.open(path, 'rb') as handle:
		expected = handle.read()
	got, want = data.split(b'\n'), expected.split(b'\n')
	n = abs(len(got) - len(want)) + sum(a != b for a, b in zip(got, want))
	first = next((i for i, (a, b) in enumerate(zip(got, want)) if a != b),
		min(len(got), len(want)))
	log('ORDER MISMATCH {}: {} line positions differ ({} lines here, {} in v5.1 '
		'-p 1); first at line {}'.format(case, n, len(got), len(want), first + 1))
	if first < len(got):
		log('  here:       ', got[first][:200].decode(errors='replace'))
	if first < len(want):
		log('  v5.1 -p 1:  ', want[first][:200].decode(errors='replace'))
	return n


# --- the seed's extra input -----------------------------------------------------

def read_fasta(path):
	records, name, chunks = [], None, []
	with open(path) as handle:
		for line in handle:
			if line.startswith('>'):
				if name is not None:
					records.append((name, ''.join(chunks)))
				name, chunks = line[1:].rstrip('\n'), []
			else:
				chunks.append(line.rstrip('\n'))
	if name is not None:
		records.append((name, ''.join(chunks)))
	return records


def read_motifs(path):
	motifs = []
	with open(path) as handle:
		for line in handle:
			if line.startswith('>'):
				motifs.append([line])
			elif motifs and line.strip():
				motifs[-1].append(line)
	return motifs


def seed_inputs(directory, seed):
	"""300 peaks and 60 of the 1,049 known motifs drawn by the seed, so each run also checks a
	different input. Deterministic in the seed alone."""

	rng = numpy.random.RandomState(seed)
	peaks = read_fasta(inp('rest.fa'))
	motifs = read_motifs(inp('known1049.motifs'))
	pick = numpy.sort(rng.choice(len(peaks), 300, replace=False))
	chosen = numpy.sort(rng.choice(len(motifs), 60, replace=False))

	fasta = os.path.join(directory, 'seed_peaks.fa')
	with open(fasta, 'w') as handle:
		for i in pick:
			handle.write('>{}\n{}\n'.format(*peaks[i]))
	motif_path = os.path.join(directory, 'seed.motifs')
	with open(motif_path, 'w') as handle:
		for i in chosen:
			handle.writelines(motifs[i])
	return fasta, motif_path


# --- the cases ------------------------------------------------------------------

K472 = inp('known472.motifs')
FIND, MSCORE, MASK = ('n_find_mismatch_vs_v51', 'n_mscore_mismatch_vs_v51',
	'n_mask_mismatch_vs_v51')

# The objective's input; its warmup and timed calls are both checked.
OBJECTIVE = ('big270k_find_p4', 'find', inp('big270k.fa'), K472, 4, ())

# Further inputs, each timed once and checked: (measurement, case, mode,
# fasta, motifs, threads, extra arguments, sorted check). A case with no
# sorted check has no sorted reference (v5.1's output at that -p varies)
# and is checked only against v5.1 -p 1.
CASES = [
	('seconds_bench5k', 'bench5k_find_p4', 'find', inp('bench5k.fa'), K472, 4, (), FIND),
	('seconds_1thread', 'bench1k_find_p1', 'find', inp('bench1k.fa'), K472, 1, (), FIND),
	('seconds_16threads', 'bench5k_find_p16', 'find', inp('bench5k.fa'), K472, 16, (), FIND),
	('seconds_many_motifs', 'many_motifs_find_p4', 'find', inp('many_motifs200.fa'),
		inp('known1049.motifs'), 4, (), FIND),
	('seconds_long', 'long_find_p4', 'find', inp('long.fa'), inp('known50.motifs'), 4, (),
		FIND),
	('seconds_degenerate', 'degenerate_find_p4', 'find', inp('degenerate.fa'),
		inp('degenerate.motifs'), 4, (), FIND),
	('seconds_options', 'options_strand_plus', 'find', inp('options500.fa'), K472, 4,
		('-strand', '+'), FIND),
	('seconds_options', 'options_strand_minus', 'find', inp('options500.fa'), K472, 4,
		('-strand', '-'), FIND),
	('seconds_options', 'options_offset_0', 'find', inp('options500.fa'), K472, 4,
		('-offset', '0'), FIND),
	('seconds_mscore', 'bench1k_mscore_p4', 'find', inp('bench1k.fa'), K472, 4,
		('-mscore',), MSCORE),
	(None, 'degenerate_mscore_p4', 'find', inp('degenerate.fa'),
		inp('degenerate.motifs'), 4, ('-mscore',), MSCORE),
	('seconds_mask', 'mask500_mask_p1', 'mask', inp('mask500.fa'), K472, 1, (), MASK),
	(None, 'degenerate_mask_p1', 'mask', inp('degenerate.fa'),
		inp('degenerate.motifs'), 1, (), MASK),
	(None, 'mask500_mask_p4', 'mask', inp('mask500.fa'), K472, 4, (), None),
]


def seed_cases(seed, fasta, motifs):
	"""The seed's extra input, in all three modes, plus mask at -p 4."""

	return [(None, 'seed{}_{}'.format(seed, case), mode, fasta, motifs, threads, extra, sorted_check)
		for case, mode, threads, extra, sorted_check in [
			('find_p4', 'find', 4, (), FIND),
			('mscore_p4', 'find', 4, ('-mscore',), MSCORE),
			('mask_p1', 'mask', 1, (), MASK),
			('mask_p4', 'mask', 4, (), None)]]


# --- the run --------------------------------------------------------------------

binary = build()
checks = {FIND: 0, MSCORE: 0, MASK: 0, 'nonportable_build': nonportable_build(),
	'n_order_mismatch_vs_v51_p1': 0}
measurements = {}


def check(case, data, sorted_check):
	order = compare_order(case, data)
	checks['n_order_mismatch_vs_v51_p1'] += order
	if sorted_check is not None and not (order == 0 and case == OBJECTIVE[0]):
		# bytes identical to v5.1 -p 1 imply identical sorted lines; the
		# sorted check of the big output runs only when they differ
		checks[sorted_check] += compare(case, data)


# The objective: warm up (under /usr/bin/time, for peak RSS), then one timed call.
case, mode, fasta, motifs, threads, extra = OBJECTIVE
_, rss, _, warm = homer(binary, mode, fasta, motifs, threads, extra, peak_rss=True)
check(case, warm, FIND)
del warm

seconds, _, cpu, data = homer(binary, mode, fasta, motifs, threads, extra)
check(case, data, FIND)
measurements['peak_rss_mb'] = rss
measurements['cpu_seconds'] = cpu
measurements['n_sites'] = data.count(b'\n')
log('objective: {:.3f} s, {:.0f} MB peak RSS (warmup), {:.1f} CPU s, {} sites'.format(
	seconds, rss, cpu, measurements['n_sites']))
del data

for measure, case, mode, fasta, motifs, threads, extra, sorted_check in CASES:
	t, _, _, data = homer(binary, mode, fasta, motifs, threads, extra)
	check(case, data, sorted_check)
	if measure is not None:
		measurements[measure] = measurements.get(measure, 0.0) + t
	log('{:28s} {:7.3f} s'.format(case, t))

with tempfile.TemporaryDirectory(dir='/dev/shm', prefix='homerseed-') as tmp:
	for _, case, mode, fasta, motifs, threads, extra, sorted_check in seed_cases(args.seed,
			*seed_inputs(tmp, args.seed)):
		_, _, _, data = homer(binary, mode, fasta, motifs, threads, extra)
		check(case, data, sorted_check)

log('checks:', json.dumps(checks))

with open(os.path.join(args.out, 'metrics.json'), 'w') as handle:
	json.dump({'objective': seconds, 'checks': checks, 'measurements': measurements},
		handle, indent=1)

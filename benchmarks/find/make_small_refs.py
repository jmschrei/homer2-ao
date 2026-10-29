"""HOMER v5.1 references for the homer2 find benchmark's small cases.

	HOMER_FIND_REF=<reference> python3 make_small_refs.py --homer2 <HOMER v5.1 homer2> \\
		[--seeds 0-9] [--jobs 4]

For every fixed case of bench.py and every seed's extra input, writes
<reference>/order_p1/<case>.txt.gz (+ .sha256), v5.1's unsorted output at
-p 1, and for each case with a sorted check <reference>/outputs/<case>.sorted.txt.gz
(+ .sha256), its sorted lines at the case's own -p. Existing references are
kept. The objective's (big270k) come from make_refs.py. The cases below are
bench.py's and must stay equal to them.
"""

import os
import sys
import time
import gzip
import hashlib
import argparse
import tempfile
import subprocess
import concurrent.futures

import numpy

parser = argparse.ArgumentParser()
parser.add_argument('--homer2', required=True)
parser.add_argument('--seeds', default='0-9', help='seeds whose extra input gets references, as A-B')
parser.add_argument('--jobs', type=int, default=4, help='homer2 calls run at once')
args = parser.parse_args()

REF = os.environ['HOMER_FIND_REF']
INPUTS = os.path.join(os.environ.get('HOMER_FIND_INPUTS', REF), 'inputs')
OUTPUTS = os.path.join(REF, 'outputs')
ORDER = os.path.join(REF, 'order_p1')


def log(*items):
	print(*items, flush=True)


def inp(name):
	return os.path.join(INPUTS, name)


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


# --- the references -------------------------------------------------------------

def write(path, blob):
	partial = '{}.tmp{}'.format(path, os.getpid())
	with gzip.open(partial, 'wb', compresslevel=3) as handle:
		handle.write(blob)
	with open(path + '.sha256.tmp{}'.format(os.getpid()), 'w') as handle:
		handle.write(hashlib.sha256(blob).hexdigest() + '\n')
	os.replace(path + '.sha256.tmp{}'.format(os.getpid()), path + '.sha256')
	os.replace(partial, path)


def one(item):
	"""The -p 1 order reference and, with a sorted check, the sorted reference
	at the case's own -p, each unless it exists."""

	case, mode, fasta, motifs, threads, extra, sorted_check = item
	done = []
	path = os.path.join(ORDER, case + '.txt.gz')
	if os.path.exists(path):
		done.append('order kept')
	else:
		_, _, _, data = homer(args.homer2, mode, fasta, motifs, 1, extra)
		write(path, data)
		done.append('order {} lines'.format(data.count(b'\n')))
	path = os.path.join(OUTPUTS, case + '.sorted.txt.gz')
	if sorted_check is None:
		pass
	elif os.path.exists(path):
		done.append('sorted kept')
	else:
		_, _, _, data = homer(args.homer2, mode, fasta, motifs, threads, extra)
		lines = data.split(b'\n')
		if lines and lines[-1] == b'':
			lines.pop()
		lines.sort()
		write(path, b'\n'.join(lines) + b'\n')
		done.append('sorted {} lines at -p {}'.format(len(lines), threads))
	return case, ', '.join(done)


os.makedirs(ORDER, exist_ok=True)
os.makedirs(OUTPUTS, exist_ok=True)
first, last = [int(s) for s in args.seeds.split('-')]
start = time.perf_counter()
with tempfile.TemporaryDirectory(dir='/dev/shm', prefix='homerrefs-') as tmp:
	todo = [c[1:] for c in CASES]
	for seed in range(first, last + 1):
		directory = os.path.join(tmp, str(seed))
		os.makedirs(directory)
		todo += [c[1:] for c in seed_cases(seed, *seed_inputs(directory, seed))]
	with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
		for case, what in pool.map(one, todo):
			log('{:28s} {}'.format(case, what))
log('{} cases in {:.0f} s'.format(len(todo), time.perf_counter() - start))

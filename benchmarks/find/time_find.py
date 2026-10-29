"""Time one homer2 call from a HOMER tree, for paired timing of two trees.

Builds `<src>`'s homer2 with this directory's Makefile into
`<src>/.paired-cache/<hash>/` (the hash covers that Makefile and cpp/), runs
one call in a fresh /dev/shm directory and prints the wall seconds last. With
--order the output's bytes must equal HOMER v5.1's at -p 1 (which implies its
sorted lines match); without it, its sorted lines must equal v5.1's. Exits
non-zero on a wrong output, so a paired timer stops. Reads the references in
$HOMER_FIND_REF.

	python3 time_find.py --src <tree> --order [--case big270k_find_p4]

Cases: big270k_find_p4 (the objective), idr_find_p4 (172,732 peaks) and the
small fixed cases (bench5k_find_p4, ..., all_find_p4). The all_ and idr_
cases read validation_order_p1/ and validation_outputs/ references, which
nothing in this directory writes.
"""

import os
import sys
import time
import fcntl
import hashlib
import argparse
import tempfile
import subprocess

REF = os.environ['HOMER_FIND_REF']
MAKEFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Makefile')
INPUTS = os.path.join(REF, 'inputs')

CASES = {
	'big270k_find_p4': ('find', 'big270k.fa', 'known472.motifs', 4, ()),
	'idr_find_p4': ('find', 'idr172k.fa', 'known472.motifs', 4, ()),
	'bench5k_find_p4': ('find', 'bench5k.fa', 'known472.motifs', 4, ()),
	'bench1k_find_p1': ('find', 'bench1k.fa', 'known472.motifs', 1, ()),
	'bench5k_find_p16': ('find', 'bench5k.fa', 'known472.motifs', 16, ()),
	'many_motifs_find_p4': ('find', 'many_motifs200.fa', 'known1049.motifs', 4, ()),
	'long_find_p4': ('find', 'long.fa', 'known50.motifs', 4, ()),
	'degenerate_find_p4': ('find', 'degenerate.fa', 'degenerate.motifs', 4, ()),
	'bench1k_mscore_p4': ('find', 'bench1k.fa', 'known472.motifs', 4, ('-mscore',)),
	'mask500_mask_p1': ('mask', 'mask500.fa', 'known472.motifs', 1, ()),
	'all_find_p4': ('find', 'peaks_all.fa', 'known472.motifs', 4, ()),
}

parser = argparse.ArgumentParser()
parser.add_argument('--src', required=True)
parser.add_argument('--case', default='big270k_find_p4', choices=sorted(CASES))
parser.add_argument('--order', action='store_true',
	help='also require the unsorted output to equal HOMER v5.1 -p 1')
args = parser.parse_args()

src = os.path.abspath(args.src)
digest = hashlib.sha256()
for root in ['Makefile'] + sorted('cpp/' + n for n in os.listdir(os.path.join(src, 'cpp'))):
	path = MAKEFILE if root == 'Makefile' else os.path.join(src, root)
	if os.path.isfile(path):
		digest.update(root.encode())
		with open(path, 'rb') as handle:
			digest.update(handle.read())
build = os.path.join(src, '.paired-cache', digest.hexdigest()[:16])
binary = os.path.join(build, 'homer2')

os.makedirs(os.path.join(src, '.paired-cache'), exist_ok=True)
with open(os.path.join(src, '.paired-cache', '.lock'), 'w') as lock:
	fcntl.flock(lock, fcntl.LOCK_EX)
	if not os.path.exists(binary):
		result = subprocess.run(['make', '-f', MAKEFILE, '-C', src, '-j', '4', 'OUT=' + build],
			capture_output=True, text=True)
		if result.returncode != 0 or not os.path.exists(binary):
			sys.stderr.write(result.stdout[-3000:] + result.stderr[-3000:])
			sys.exit('build failed')

mode, fasta, motifs, threads, extra = CASES[args.case]
with tempfile.TemporaryDirectory(dir='/dev/shm', prefix='homertime-') as tmp:
	output = os.path.join(tmp, 'out.txt')
	command = [binary, mode, '-i', os.path.join(INPUTS, fasta), '-m',
		os.path.join(INPUTS, motifs), '-p', str(threads), '-o', output] + list(extra)
	start = time.perf_counter()
	result = subprocess.run(command, cwd=tmp, stdin=subprocess.DEVNULL,
		stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	seconds = time.perf_counter() - start
	if result.returncode != 0:
		sys.exit('homer2 exited {}'.format(result.returncode))
	with open(output, 'rb') as handle:
		data = handle.read()
		lines = data.split(b'\n')

validation = args.case.startswith(('all_', 'idr_'))
if args.order:
	order_ref = os.path.join(REF, 'validation_order_p1' if validation else 'order_p1',
		args.case + '.txt.gz.sha256')
	with open(order_ref) as handle:
		if handle.read().strip() != hashlib.sha256(data).hexdigest():
			sys.exit('output order differs from HOMER v5.1 -p 1 on {}'.format(args.case))
	# bytes identical to v5.1 -p 1: the sorted lines are too
	print('{} {:.4f}'.format(args.case, seconds))
	sys.exit(0)

if lines and lines[-1] == b'':
	lines.pop()
lines.sort()
got = hashlib.sha256(b'\n'.join(lines) + b'\n').hexdigest()
with open(os.path.join(REF, 'validation_outputs' if validation else 'outputs',
		args.case + '.sorted.txt.gz.sha256')) as handle:
	if handle.read().strip() != got:
		sys.exit('output differs from HOMER v5.1 on {}'.format(args.case))

print('{} {:.4f}'.format(args.case, seconds))

"""HOMER v5.1 homer2 references for a large find case, built by motif groups.

At -p 1 HOMER v5.1's homer2 writes find's lines motif by motif in file order, so
running it on consecutive groups of the motif file and concatenating the
outputs in group order gives its whole-file -p 1 output; --verify checks
that against an existing whole-file reference before it is trusted.

	python3 make_refs.py --homer2 <HOMER v5.1 homer2> --fasta <fa> --motifs <m> \\
		--case <name> --ref <reference dir> --groups 48 --jobs 24 [--verify <sha256 file>]

Writes <ref>/order_p1/<case>.txt.gz (+ .sha256), the unsorted -p 1 output,
and <ref>/outputs/<case>.sorted.txt.gz (+ .sha256), its lines sorted.
"""

import os
import gzip
import hashlib
import argparse
import tempfile
import subprocess
import concurrent.futures

parser = argparse.ArgumentParser()
parser.add_argument('--homer2', required=True)
parser.add_argument('--fasta', required=True)
parser.add_argument('--motifs', required=True)
parser.add_argument('--case', required=True)
parser.add_argument('--ref', required=True)
parser.add_argument('--groups', type=int, default=48)
parser.add_argument('--jobs', type=int, default=24)
parser.add_argument('--verify', default=None, help='sha256 file of a whole-file v5.1 -p 1 output; check only, write nothing')
parser.add_argument('--order-only', action='store_true')
args = parser.parse_args()

motifs = []
with open(args.motifs) as handle:
	for line in handle:
		if line.startswith('>'):
			motifs.append([line])
		elif motifs:
			motifs[-1].append(line)
k = min(args.groups, len(motifs))
bounds = [round(i * len(motifs) / k) for i in range(k + 1)]


def one(g):
	with tempfile.TemporaryDirectory(dir='/dev/shm', prefix='v51ref-') as tmp:
		mpath = os.path.join(tmp, 'group.motifs')
		with open(mpath, 'w') as handle:
			for m in motifs[bounds[g]:bounds[g+1]]:
				handle.writelines(m)
		out = os.path.join(tmp, 'out.txt')
		subprocess.run([args.homer2, 'find', '-i', os.path.abspath(args.fasta), '-m', mpath,
			'-p', '1', '-o', out], cwd=tmp, stdin=subprocess.DEVNULL,
			stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
		with open(out, 'rb') as handle:
			return handle.read()


with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
	parts = list(pool.map(one, range(k)))
data = b''.join(parts)
digest = hashlib.sha256(data).hexdigest()
print('{}: {} motifs in {} groups, {} lines, sha256 {}'.format(args.case, len(motifs), k,
	data.count(b'\n'), digest))

if args.verify is not None:
	with open(args.verify) as handle:
		want = handle.read().strip()
	print('verify:', 'identical to the whole-file reference' if want == digest else 'DIFFERS from ' + want)
	raise SystemExit(0 if want == digest else 1)

os.makedirs(os.path.join(args.ref, 'order_p1'), exist_ok=True)
path = os.path.join(args.ref, 'order_p1', args.case + '.txt.gz')
with gzip.open(path + '.tmp', 'wb', compresslevel=3) as handle:
	handle.write(data)
with open(path + '.sha256', 'w') as handle:
	handle.write(digest + '\n')
os.replace(path + '.tmp', path)
print('wrote', path)
if args.order_only:
	raise SystemExit(0)

lines = data.split(b'\n')
del data, parts
if lines and lines[-1] == b'':
	lines.pop()
lines.sort()
blob = b'\n'.join(lines) + b'\n'
os.makedirs(os.path.join(args.ref, 'outputs'), exist_ok=True)
path = os.path.join(args.ref, 'outputs', args.case + '.sorted.txt.gz')
with gzip.open(path + '.tmp', 'wb', compresslevel=3) as handle:
	handle.write(blob)
with open(path + '.sha256', 'w') as handle:
	handle.write(hashlib.sha256(blob).hexdigest() + '\n')
os.replace(path + '.tmp', path)
print('wrote', path)

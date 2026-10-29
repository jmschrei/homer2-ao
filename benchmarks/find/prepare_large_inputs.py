"""Build the large inputs of the homer2 find benchmark: 200 bp around each
summit, in file order, as samtools returns them (soft-masked, wrapped at 60).

	python3 prepare_large_inputs.py --out <reference>/inputs --genome <hg38.fa> \\
		--peaks <ENCFF333TAT .narrowPeak.gz> --idr-peaks <ENCFF855PCP .narrowPeak.gz>

big270k.fa: all 269,800 ENCODE ENCFF333TAT K562 pseudoreplicated peaks (the
objective). idr172k.fa: all 172,732 ENCFF855PCP K562 IDR-conservative peaks
(external validation). The small inputs come from prepare_small_inputs.py.
"""

import os
import gzip
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--out', required=True)
parser.add_argument('--genome', required=True)
parser.add_argument('--peaks', required=True)
parser.add_argument('--idr-peaks', required=True)
args = parser.parse_args()

HG38 = args.genome
SETS = [
	('big270k.fa', args.peaks, 269800),
	('idr172k.fa', args.idr_peaks, 172732),
]

for name, peaks, expected in SETS:
	regions = []
	with gzip.open(peaks, 'rt') as handle:
		for line in handle:
			cols = line.split('\t')
			start = max(0, int(cols[1]) + int(cols[9]) - 100)
			regions.append('{}:{}-{}'.format(cols[0], start + 1, start + 200))
	assert len(regions) == expected, (name, len(regions))
	path = os.path.join(args.out, '_regions.txt')
	with open(path, 'w') as handle:
		handle.write(''.join(r + '\n' for r in regions))
	with open(os.path.join(args.out, name), 'w') as handle:
		subprocess.run(['samtools', 'faidx', HG38, '-r', path, '-n', '60'],
			stdout=handle, check=True)
	os.remove(path)
	# every record 200 bp
	lengths, n = [], 0
	with open(os.path.join(args.out, name)) as handle:
		cur = None
		for line in handle:
			if line.startswith('>'):
				if cur is not None:
					lengths.append(cur)
				cur = 0
				n += 1
			else:
				cur += len(line.rstrip('\n'))
		lengths.append(cur)
	assert n == expected and all(x == 200 for x in lengths), (name, n, set(lengths))
	print('wrote', name, n, 'records')

"""Build the small inputs the homer2 find benchmark reads, once, from fixed seeds.

Writes FASTA and motif files into <out>/ (<reference>/inputs):

	python3 prepare_small_inputs.py --out <reference>/inputs --genome <hg38.fa> \\
		--peaks <medium K562 peaks .narrowPeak.gz> --homer <HOMER v5.1 tree>

--genome needs a samtools index (.fai). --peaks is every 10th line of ENCODE
ENCFF333TAT, coordinate-sorted (see ../README.md). --homer supplies
data/knownTFs.
"""

import os
import gzip
import argparse
import subprocess

import numpy

parser = argparse.ArgumentParser()
parser.add_argument('--out', required=True)
parser.add_argument('--genome', required=True)
parser.add_argument('--peaks', required=True)
parser.add_argument('--homer', required=True)
args = parser.parse_args()
os.makedirs(args.out, exist_ok=True)

HG38 = args.genome
PEAKS = args.peaks
KNOWN = os.path.join(args.homer, 'data', 'knownTFs')


def fetch(regions):
	"""Sequences for 'chr:start-end' (1-based, inclusive) regions, as samtools
	returns them: soft-masked, wrapped at 60."""

	path = os.path.join(args.out, '_regions.txt')
	with open(path, 'w') as handle:
		handle.write(''.join(r + '\n' for r in regions))

	text = subprocess.run(['samtools', 'faidx', HG38, '-r', path, '-n', '60'],
		capture_output=True, text=True, check=True).stdout
	os.remove(path)

	records = []
	for block in text.split('>')[1:]:
		lines = block.split('\n')
		records.append((lines[0], ''.join(lines[1:])))
	return records


def write_fasta(name, records, width=60):
	with open(os.path.join(args.out, name), 'w') as handle:
		for header, seq in records:
			handle.write('>{}\n'.format(header))
			for i in range(0, len(seq), width):
				handle.write(seq[i:i+width] + '\n')


def read_motifs(path):
	"""HOMER motif file -> list of (header line, [row lines])."""

	motifs = []
	with open(path) as handle:
		for line in handle:
			if line.startswith('>'):
				motifs.append([line, []])
			elif motifs and line.strip():
				motifs[-1][1].append(line)
	return motifs


def write_motifs(name, motifs):
	with open(os.path.join(args.out, name), 'w') as handle:
		for header, rows in motifs:
			handle.write(header)
			handle.writelines(rows)


# Peaks: 200 bp around each summit, in file order.
regions = []
with gzip.open(PEAKS, 'rt') as handle:
	for line in handle:
		cols = line.split('\t')
		start = max(0, int(cols[1]) + int(cols[9]) - 100)
		regions.append('{}:{}-{}'.format(cols[0], start + 1, start + 200))

peaks = fetch(regions)
assert len(peaks) == 26980, len(peaks)
assert all(len(s) == 200 for _, s in peaks)

order = numpy.random.RandomState(2026).permutation(len(peaks))
bench_idx = numpy.sort(order[:5000])
rest_idx = numpy.sort(order[5000:])

bench = [peaks[i] for i in bench_idx]
rest = [peaks[i] for i in rest_idx]

write_fasta('peaks_all.fa', peaks)                 # external validation
write_fasta('bench5k.fa', bench)                   # the objective
write_fasta('bench1k.fa', bench[:1000])            # seconds_1thread, mscore check
write_fasta('mask500.fa', bench[:500])             # mask check, -p 1
write_fasta('rest.fa', rest)                       # the seed's extra input draws from here
write_fasta('options500.fa', rest[:500])           # -strand / -offset
write_fasta('many_motifs200.fa', rest[500:700])
write_fasta('known_target2k.fa', rest[1000:3000])  # external validation: homer2 known
write_fasta('known_bg2k.fa', rest[3000:5000])

# Motifs.
known = read_motifs(os.path.join(KNOWN, 'vertebrates', 'known.motifs'))
every = read_motifs(os.path.join(KNOWN, 'all.motifs'))
assert len(known) == 472 and len(every) == 4668, (len(known), len(every))

write_motifs('known472.motifs', known)
write_motifs('all4668.motifs', every)

# HOMER's curated library for all organisms. all.motifs is not used for
# timing: most of its JASPAR entries are count matrices with threshold 0,
# which report a site at every window.
known_all = read_motifs(os.path.join(KNOWN, 'known.motifs'))
assert len(known_all) == 1049, len(known_all)
write_motifs('known1049.motifs', known_all)

rs = numpy.random.RandomState(7)
pick = numpy.sort(rs.choice(len(known), 50, replace=False))
write_motifs('known50.motifs', [known[i] for i in pick])

# Long, messy genomic chunks: soft-masked, with N runs, each under the
# 1,000,000-character buffer homer2 copies a FASTA record into.
long = fetch(['chr21:8000001-8900000', 'chr17:21600001-22500000'])
write_fasta('long.fa', long)

# Degenerate sequences.
rs = numpy.random.RandomState(11)
def random_seq(n, alphabet='ACGT'):
	return ''.join(rs.choice(list(alphabet), n))

degenerate = [
	('one_bp', 'A'),
	('five_bp', 'ACGTA'),
	('exact_width_10', random_seq(10)),
	('all_N', 'N' * 300),
	('all_lower', random_seq(300, 'acgt')),
	('mixed_case', random_seq(300, 'ACGTacgt')),
	('iupac', random_seq(400, 'ACGTRYKMSWBDHVN')),
	('n_islands', random_seq(80) + 'N' * 7 + random_seq(15) + 'N' + random_seq(120)),
	('dup', random_seq(200)),
	('dup', random_seq(200)),
	('dup', random_seq(200)),
	('name with spaces', random_seq(250)),
	('poly_A', 'A' * 300),
	('dinuc_repeat', 'CA' * 200),
	('gc_rich', random_seq(300, 'GGGCCCGCGA')),
	('one_line_50k', random_seq(50000)),
]
with open(os.path.join(args.out, 'degenerate.fa'), 'w') as handle:
	for name, seq in degenerate:
		handle.write('>{}\n'.format(name))
		if name == 'one_line_50k':
			handle.write(seq + '\n')
		elif name == 'mixed_case':
			for i in range(0, len(seq), 60):
				handle.write(seq[i:i+60] + '\r\n')
		else:
			for i in range(0, len(seq), 60):
				handle.write(seq[i:i+60] + '\n')
	handle.write('>empty_record\n')
	handle.write('>last_after_empty\n' + random_seq(120) + '\n')

# Degenerate motifs: the narrowest and widest known ones, plus synthetic
# matrices for ties, unreachable thresholds and a single column.
widths = [len(rows) for _, rows in every]
narrow = [every[i] for i in numpy.argsort(widths, kind='stable')[:3]]
wide = [every[i] for i in numpy.argsort(widths, kind='stable')[-3:]]
rs = numpy.random.RandomState(13)
some = [known[i] for i in numpy.sort(rs.choice(len(known), 20, replace=False))]

def matrix(name, rows, threshold):
	header = '>{}\t{}\t{}\t-10.0\t0\tT:0(0.0%),B:0(0.0%),P:1e0\n'.format(
		'N' * len(rows), name, threshold)
	return [header, ['\t'.join('{:.3f}'.format(v) for v in row) + '\n'
		for row in rows]]

synthetic = [
	matrix('uniform_tie_at_zero', [[0.25] * 4] * 8, 0.0),
	matrix('one_hot_TGACTCA', [[float(b == c) for c in 'ACGT']
		for b in 'TGACTCA'], 5.0),
	matrix('unreachable', [[0.7, 0.1, 0.1, 0.1]] * 6, 100.0),
	matrix('single_column', [[0.1, 0.6, 0.2, 0.1]], 0.5),
	matrix('negative_threshold', [[0.4, 0.1, 0.1, 0.4]] * 5, -3.0),
]
write_motifs('degenerate.motifs', narrow + wide + some + synthetic)

print('wrote', sorted(os.listdir(args.out)))

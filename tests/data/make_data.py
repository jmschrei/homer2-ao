#!/usr/bin/env python3
"""Write the synthetic inputs the test suite reads from tests/data/.

	python3 tests/data/make_data.py

Running it again writes the same bytes: everything comes from one seeded
random.Random, and nothing is read from outside this file.

genome.fa
	Three chromosomes (chr1 120,000 bp, chr2 90,000 bp, chr3 60,000 bp) of
	random sequence whose GC fraction changes every 3,000 bp (0.35 to 0.65),
	60 bases per line. Twenty-four stretches of 600-1,500 bp are lowercase
	(soft-masked), chr2:40,001-42,000 is N, and ten peaks are lowercase over
	their whole 200 bp.
peaks.txt
	200 peaks of 200 bp in HOMER's peak format (1-based start and end,
	both inclusive), IDs peak001-peak200 in genome order, every third peak on
	the - strand. A header line starts with '#'.
peaks.bed
	The same peaks as BED6 (0-based start, end exclusive).
targets.fa
	The 200 bp of each peak, + strand of the genome, case kept, named by ID.
background.fa
	1,000 genome windows of 200 bp (bg0001-bg1000) that contain no N and lie
	at least 200 bp from every peak, case kept.
planted.tsv
	Every planted site: peak, motif, chromosome, 0-based start, end, strand,
	and whether it is lowercase.
motifs/
	Motif files in HOMER's format, written out below with their thresholds.

Three sites are planted, each as its exact consensus, in the orientation the
peak's strand column does not decide (alternating + and - by plant order):

	siteA  TAGCGTACCGTA (12 bp)          120 peaks, centred on the peak centre
	siteB  GGATCCTTAC (10 bp)            80 peaks, 40 bp right of the centre
	REST   GGAGCTGTCCATGGTGCTGA (20 bp)  60 peaks, 60 bp left of the centre

The REST site is the consensus of HOMER's REST-NRSF known motif; siteA and
siteB are made up and match no motif in HOMER's library. Ten of the siteA
peaks are lowercase. Random sequence can add chance matches to any motif, so
tests count matches by scanning the files, and planted.tsv is a lower bound.
"""

import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261028

CHROMOSOMES = [('chr1', 120000), ('chr2', 90000), ('chr3', 60000)]
N_RUN = ('chr2', 40000, 42000)   # 0-based, end exclusive
PEAK_WIDTH = 200
PEAKS_PER_CHROMOSOME = {'chr1': 90, 'chr2': 65, 'chr3': 45}

SITES = {
	'siteA': 'TAGCGTACCGTA',
	'siteB': 'GGATCCTTAC',
	'REST': 'GGAGCTGTCCATGGTGCTGA',
}
# Offset of each site's first base from the window's first base (0-based).
SITE_OFFSET = {'siteA': 94, 'siteB': 140, 'REST': 30}

COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def revcomp(seq):
	return seq.translate(COMPLEMENT)[::-1]


def probability_rows(consensus, top):
	"""One row per base: `top` for the consensus base, the rest split evenly."""

	low = (1.0 - top) / 3.0
	rows = []
	for base in consensus:
		rows.append([top if b == base else low for b in 'ACGT'])
	return rows


def motif_text(consensus, name, threshold, rows):
	lines = ['>{}\t{}\t{:.6f}'.format(consensus, name, threshold)]
	for row in rows:
		lines.append('\t'.join('{:.3f}'.format(value) for value in row))
	return '\n'.join(lines) + '\n'


MOTIFS = {
	# 0.91 / 0.03: an exact match scores 15.50, one mismatch 12.09,
	# two 8.68, so threshold 10.0 accepts at most one mismatch.
	'siteA': ('TAGCGTACCGTA', 'siteA', 10.0, probability_rows('TAGCGTACCGTA', 0.91)),
	# 0.85 / 0.05: exact 12.24, one mismatch 9.40, two 6.57; threshold 9.0.
	'siteB': ('GGATCCTTAC', 'siteB', 9.0, probability_rows('GGATCCTTAC', 0.85)),
	# A palindrome, so every match is found on both strands; 0.94 / 0.02:
	# exact 7.95, one mismatch 4.10; threshold 6.0 accepts exact matches only.
	'ebox': ('CACGTG', 'ebox', 6.0, probability_rows('CACGTG', 0.94)),
	# Two motifs that were not planted.
	'decoy1': ('CTTTAAGGCA', 'decoy1', 8.0, probability_rows('CTTTAAGGCA', 0.70)),
	'decoy2': ('ATATGCGCAT', 'decoy2', 8.5, probability_rows('ATATGCGCAT', 0.75)),
}


def make_genome(rng):
	genome = {}
	for chrom, length in CHROMOSOMES:
		bases = []
		for block in range(0, length, 3000):
			gc = rng.uniform(0.35, 0.65)
			for _ in range(min(3000, length - block)):
				if rng.random() < gc:
					bases.append(rng.choice('GC'))
				else:
					bases.append(rng.choice('AT'))
		genome[chrom] = bases
	chrom, start, end = N_RUN
	for i in range(start, end):
		genome[chrom][i] = 'N'
	return genome


def place_peaks(rng):
	"""Window starts (0-based) per chromosome, at least 500 bp apart, 1,000 bp
	from the chromosome ends and 500 bp from the N run."""

	peaks = []
	for chrom, length in CHROMOSOMES:
		starts = []
		while len(starts) < PEAKS_PER_CHROMOSOME[chrom]:
			start = rng.randrange(1000, length - 1000 - PEAK_WIDTH)
			if chrom == N_RUN[0] and N_RUN[1] - 500 - PEAK_WIDTH < start < N_RUN[2] + 500:
				continue
			if any(abs(start - other) < 500 for other in starts):
				continue
			starts.append(start)
		peaks += [(chrom, start) for start in sorted(starts)]
	return peaks


def mask_stretches(rng, genome, peak_windows):
	"""Lowercase 24 stretches that touch no peak window."""

	stretches = []
	while len(stretches) < 24:
		chrom, length = rng.choice(CHROMOSOMES)
		size = rng.randrange(600, 1501)
		start = rng.randrange(0, length - size)
		end = start + size
		if any(c == chrom and start < s + PEAK_WIDTH + 100 and s - 100 < end
				for c, s in peak_windows):
			continue
		if chrom == N_RUN[0] and start < N_RUN[2] and N_RUN[1] < end:
			continue
		stretches.append((chrom, start, end))
	for chrom, start, end in stretches:
		for i in range(start, end):
			genome[chrom][i] = genome[chrom][i].lower()


def write_fasta(path, records, width=60):
	with open(path, 'w') as out:
		for name, seq in records:
			out.write('>' + name + '\n')
			for i in range(0, len(seq), width):
				out.write(seq[i:i + width] + '\n')


def main():
	rng = random.Random(SEED)
	genome = make_genome(rng)
	windows = place_peaks(rng)
	mask_stretches(rng, genome, windows)

	ids = ['peak{:03d}'.format(i + 1) for i in range(len(windows))]
	strands = ['-' if i % 3 == 2 else '+' for i in range(len(windows))]

	order = list(range(len(windows)))
	rng.shuffle(order)
	carriers = {
		'siteA': order[0:120],
		'siteB': order[80:160],
		'REST': order[140:200],
	}
	lowercase_peaks = set(order[0:10])

	planted = []
	for motif in ('siteA', 'siteB', 'REST'):
		for k, index in enumerate(carriers[motif]):
			chrom, w0 = windows[index]
			strand = '+' if k % 2 == 0 else '-'
			site = SITES[motif] if strand == '+' else revcomp(SITES[motif])
			start = w0 + SITE_OFFSET[motif]
			for j, base in enumerate(site):
				genome[chrom][start + j] = base
			planted.append([ids[index], motif, chrom, start, start + len(site), strand,
				'yes' if index in lowercase_peaks else 'no'])

	for index in lowercase_peaks:
		chrom, w0 = windows[index]
		for i in range(w0, w0 + PEAK_WIDTH):
			genome[chrom][i] = genome[chrom][i].lower()

	sequences = {chrom: ''.join(bases) for chrom, bases in genome.items()}
	write_fasta(os.path.join(HERE, 'genome.fa'), [(c, sequences[c]) for c, _ in CHROMOSOMES])

	with open(os.path.join(HERE, 'peaks.txt'), 'w') as out:
		out.write('#PeakID\tchr\tstart\tend\tstrand\tNormalized Tag Count\n')
		for i, (chrom, w0) in enumerate(windows):
			out.write('{}\t{}\t{}\t{}\t{}\t{:.1f}\n'.format(ids[i], chrom, w0 + 1,
				w0 + PEAK_WIDTH, strands[i], 10.0 + (i * 7) % 50))
	with open(os.path.join(HERE, 'peaks.bed'), 'w') as out:
		for i, (chrom, w0) in enumerate(windows):
			out.write('{}\t{}\t{}\t{}\t{:.1f}\t{}\n'.format(chrom, w0, w0 + PEAK_WIDTH,
				ids[i], 10.0 + (i * 7) % 50, strands[i]))

	write_fasta(os.path.join(HERE, 'targets.fa'), [(ids[i],
		sequences[chrom][w0:w0 + PEAK_WIDTH]) for i, (chrom, w0) in enumerate(windows)])

	background = []
	taken = []
	while len(background) < 1000:
		chrom, length = rng.choice(CHROMOSOMES)
		w0 = rng.randrange(0, length - PEAK_WIDTH)
		seq = sequences[chrom][w0:w0 + PEAK_WIDTH]
		if 'N' in seq:
			continue
		if any(c == chrom and abs(w0 - s) < PEAK_WIDTH + 200 for c, s in windows):
			continue
		if (chrom, w0) in taken:
			continue
		taken.append((chrom, w0))
		background.append(('bg{:04d}'.format(len(background) + 1), seq))
	write_fasta(os.path.join(HERE, 'background.fa'), background)

	with open(os.path.join(HERE, 'planted.tsv'), 'w') as out:
		out.write('peak\tmotif\tchr\tstart\tend\tstrand\tlowercase\n')
		for row in sorted(planted, key=lambda r: (r[1], r[0])):
			out.write('\t'.join(str(v) for v in row) + '\n')

	os.makedirs(os.path.join(HERE, 'motifs'), exist_ok=True)
	for key, (consensus, name, threshold, rows) in MOTIFS.items():
		with open(os.path.join(HERE, 'motifs', key + '.motif'), 'w') as out:
			out.write(motif_text(consensus, name, threshold, rows))
	with open(os.path.join(HERE, 'motifs', 'known_small.motifs'), 'w') as out:
		for key in ('siteA', 'siteB', 'ebox', 'decoy1', 'decoy2'):
			out.write(motif_text(*MOTIFS[key]))
	with open(os.path.join(HERE, 'motifs', 'siteA_siteB.motifs'), 'w') as out:
		for key in ('siteA', 'siteB'):
			out.write(motif_text(*MOTIFS[key]))


if __name__ == '__main__':
	main()

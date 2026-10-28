"""Python references for the nucleotide statistics HOMER's C++ programs
print: the position frequency table of `homerTools freq` (also makeTagDirectory's
tagFreq.txt), the per-sequence CpG/GC line of `homerTools freq -gc`, and the
GC-content histograms makeTagDirectory writes with -checkGC.

Every value is a ratio of whole-number counts, divided once in double
precision as the C++ code does, so formatting it with '%f' gives the same
text the program prints."""

import math

BASES = 'ACGT'
INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}

FREQ_HEADER = ['Offset', 'A', 'C', 'G', 'T', 'A/T', 'C/G', 'A/G', 'C/T'] + \
	[a + b for a in BASES for b in BASES] + ['N', 'NN']


def freq_table(records, length, offset=0):
	"""Rows of the position frequency table, as lists of strings, for
	`records` [(sequence, weight)]: one row per position 0..length-1, its
	label the position plus `offset`. A position counts the weight of every
	sequence holding A, C, G or T there (either case); a dinucleotide counts
	at its first base, when both bases are A, C, G or T. A/T, C/G, A/G and C/T
	are the fractions of dinucleotides made of those two bases only. N and NN
	are the base and dinucleotide totals."""

	mono = [[0.0] * 4 for _ in range(length)]
	di = [[0.0] * 16 for _ in range(length)]
	totals = [0.0] * length
	ditotals = [0.0] * length
	for seq, weight in records:
		for i in range(min(len(seq), length)):
			a = INDEX.get(seq[i])
			if a is None:
				continue
			mono[i][a] += weight
			totals[i] += weight
			if i + 1 < len(seq):
				b = INDEX.get(seq[i + 1])
				if b is not None:
					di[i][a * 4 + b] += weight
					ditotals[i] += weight
	rows = []
	for i in range(length):
		row = [str(i + offset)]
		t, dt, d = totals[i], ditotals[i], di[i]
		row += ['0.0' if t < 0.001 else '%f' % (mono[i][j] / t) for j in range(4)]
		if dt < 0.001:
			row += ['0.0'] * 4
		else:
			row += ['%f' % ((d[0] + d[3] + d[12] + d[15]) / dt),
				'%f' % ((d[5] + d[6] + d[9] + d[10]) / dt),
				'%f' % ((d[0] + d[2] + d[8] + d[10]) / dt),
				'%f' % ((d[5] + d[7] + d[13] + d[15]) / dt)]
		row += ['0.0' if dt < 0.001 else '%f' % (d[j] / dt) for j in range(16)]
		row += ['%.1f' % t, '%.1f' % dt]
		rows.append(row)
	return rows


def gc_line(name, seq):
	"""The `freq -gc` line of a sequence: name, CpG, GC, AG, AC and length.
	Only A, C, G and T (either case) count: GC, AG and AC are fractions of
	those bases, CpG is the number of CG dinucleotides over the number of
	dinucleotides made of two such bases, and length is the number of such
	bases. None when the sequence has none."""

	upper = seq.upper()
	n = sum(1 for b in upper if b in BASES)
	if n == 0:
		return None
	pairs = [upper[i:i + 2] for i in range(len(upper) - 1)]
	nn = sum(1 for p in pairs if p[0] in BASES and p[1] in BASES)
	cpg = sum(1 for p in pairs if p == 'CG')
	count = lambda chars: sum(1 for b in upper if b in chars)
	return [name, '%f' % (cpg / max(nn, 0.01)), '%f' % (count('CG') / n),
		'%f' % (count('AG') / n), '%f' % (count('AC') / n), '%d' % n]


def gc_bins(fragment_length):
	"""(bin width, number of bins) of a GC histogram over windows of
	`fragment_length`, as NucleotideFreq::initMatrix sizes it."""

	inc = 1.0 / fragment_length
	limit = int((1.0 / inc) + 4.0)
	z, total = 0, 0.0
	while z < limit:
		total += inc
		if total >= 1.0:
			break
		z += 1
	return inc, z + 2


def gc_bin(gc, total, inc):
	return int(math.floor((gc / total) / inc + 0.01))


def gc_histogram_rows(counts, inc):
	"""The rows of tagGCcontent.txt / genomeGCcontent.txt: bin start (3
	decimals), total weight and that weight as a density."""

	total = sum(counts)
	return [['%.3f' % (inc * i), '%.1f' % c, '%f' % (c / total / inc)]
		for i, c in enumerate(counts)]

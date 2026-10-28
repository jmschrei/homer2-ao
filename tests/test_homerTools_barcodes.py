"""homerTools barcodes and homerTools truseq: splitting FASTQ reads by the
barcode at the start of each read, and counting the TruSeq index at the end
of each read name.

The tests write small FASTQ files with known barcode counts into their own
directory (both commands write beside the input: <input>.freq.txt, and for
barcodes <input>.<barcode> per barcode) and check the files against those
counts: a barcode's frequency is its read count over all counted reads,
printed with six decimals.
"""

import os

import pytest

import homer_helpers as hh

INSERT = 'TTGACCAGTCAGGT'


def barcode_reads():
	"""(name, sequence, quality): 32 AAAA, 18 CCCC, 9 GGGG and 1 TTTT reads,
	interleaved; qualities vary inside the insert only."""

	plan = ['AAAA'] * 32 + ['CCCC'] * 18 + ['GGGG'] * 9 + ['TTTT']
	order = sorted(range(len(plan)), key=lambda i: (i * 37) % len(plan))
	reads = []
	for k, i in enumerate(order):
		insert = INSERT[k % 5:] + INSERT[:k % 5]
		qual = 'IIII' + ''.join(chr(53 + (k + j) % 20) for j in range(len(insert)))
		reads.append(('read{:02d} extra'.format(k), plan[i] + insert, qual))
	return reads


def write_fastq(path, reads):
	with open(str(path), 'w') as out:
		for name, seq, qual in reads:
			out.write('@{}\n{}\n+\n{}\n'.format(name, seq, qual))
	return str(path)


def read_fastq(path):
	lines = open(str(path)).read().splitlines()
	assert len(lines) % 4 == 0
	return [(lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


def read_freq(path):
	"""(total from the header, [(barcode, count, frequency text)])."""

	lines = [line.split('\t') for line in open(str(path)).read().splitlines()]
	header = lines[0]
	assert header[0] == 'Barcode' and header[2] == 'Frequency'
	assert header[1].startswith('Total Reads(of ') and header[1].endswith(')')
	total = int(header[1][len('Total Reads(of '):-1])
	return total, [(row[0], int(row[1]), row[2]) for row in lines[1:]]


def expected_freq(counts):
	total = sum(counts.values())
	return total, {(bc, n, '%f' % (n / total)) for bc, n in counts.items()}


def check_freq(path, counts):
	"""The table holds every barcode once, most frequent first."""

	total, rows = read_freq(path)
	assert (total, set(rows)) == expected_freq(counts)
	assert len(rows) == len(counts)
	assert [r[1] for r in rows] == sorted((r[1] for r in rows), reverse=True)


def split_files(directory, name):
	prefix = name + '.'
	return sorted(f[len(prefix):] for f in os.listdir(str(directory))
		if f.startswith(prefix) and f != name + '.freq.txt')


COUNTS = {'AAAA': 32, 'CCCC': 18, 'GGGG': 9, 'TTTT': 1}


def test_barcodes_split_reads(tmp_path):
	"""Each read goes, without its first 4 bases and their qualities, to
	<input>.<barcode>, in input order; a barcode under 2% of the reads (TTTT,
	1 of 60) gets no file."""

	reads = barcode_reads()
	write_fastq(tmp_path / 'lib.fq', reads)
	result = hh.run(['homerTools', 'barcodes', 4, 'lib.fq'], tmp_path)
	assert result.stdout == ''
	assert split_files(tmp_path, 'lib.fq') == ['AAAA', 'CCCC', 'GGGG']
	for barcode in ('AAAA', 'CCCC', 'GGGG'):
		want = [(n, s[4:], q[4:]) for n, s, q in reads if s.startswith(barcode)]
		assert read_fastq(tmp_path / ('lib.fq.' + barcode)) == want
	check_freq(tmp_path / 'lib.fq.freq.txt', COUNTS)


@pytest.mark.parametrize('minimum, kept', [
	('0', ['AAAA', 'CCCC', 'GGGG', 'TTTT']),
	('0.2', ['AAAA', 'CCCC']),
	('0.5', ['AAAA']),
])
def test_barcodes_min(tmp_path, minimum, kept):
	"""-min <#> keeps the files of barcodes with at least that fraction of
	the reads (32/60, 18/60, 9/60, 1/60); the frequency table lists every
	barcode."""

	write_fastq(tmp_path / 'lib.fq', barcode_reads())
	hh.run(['homerTools', 'barcodes', 4, '-min', minimum, 'lib.fq'], tmp_path)
	assert split_files(tmp_path, 'lib.fq') == kept
	check_freq(tmp_path / 'lib.fq.freq.txt', COUNTS)


@pytest.mark.parametrize('length', [1, 6])
def test_barcode_length(tmp_path, length):
	reads = barcode_reads()
	write_fastq(tmp_path / 'lib.fq', reads)
	hh.run(['homerTools', 'barcodes', length, '-min', 0, 'lib.fq'], tmp_path)
	counts = {}
	for _, seq, _ in reads:
		counts[seq[:length]] = counts.get(seq[:length], 0) + 1
	assert split_files(tmp_path, 'lib.fq') == sorted(counts)
	check_freq(tmp_path / 'lib.fq.freq.txt', counts)
	for barcode in counts:
		want = [(n, s[length:], q[length:]) for n, s, q in reads if s[:length] == barcode]
		assert read_fastq(tmp_path / ('lib.fq.' + barcode)) == want


def quality_reads():
	"""Barcode qualities: IIII (Phred 40), II5I (one base at Phred 20) and
	##II (Phred 2); the insert is always low quality."""

	reads = []
	for k, (barcode, qual) in enumerate([('ACGT', 'IIII'), ('ACGT', 'II5I'), ('ACGT', '##II'),
			('TGCA', 'IIII'), ('TGCA', 'II5I'), ('TGCA', 'IIII'), ('GGCC', '##II')]):
		reads.append(('q{}'.format(k), barcode + INSERT, qual + '#' * len(INSERT)))
	return reads


@pytest.mark.parametrize('options, kept', [
	([], ['q0', 'q1', 'q2', 'q3', 'q4', 'q5', 'q6']),
	(['-qual', 20, '-qualBase', '!'], ['q0', 'q1', 'q3', 'q4', 'q5']),
	(['-qual', 21, '-qualBase', '!'], ['q0', 'q3', 'q5']),
	(['-qual', 0], ['q0', 'q3', 'q5']),
	(['-qual', 0, '-qualBase', '!'], ['q0', 'q1', 'q2', 'q3', 'q4', 'q5', 'q6']),
])
def test_barcodes_quality(tmp_path, options, kept):
	"""-qual <#> drops reads whose lowest barcode quality, the character
	minus -qualBase (default B), is below #; dropped reads are not counted.
	Qualities after the barcode do not matter."""

	reads = quality_reads()
	write_fastq(tmp_path / 'lib.fq', reads)
	hh.run(['homerTools', 'barcodes', 4, '-min', 0] + options + ['lib.fq'], tmp_path)
	got = []
	counts = {}
	for barcode in split_files(tmp_path, 'lib.fq'):
		for name, seq, qual in read_fastq(tmp_path / ('lib.fq.' + barcode)):
			got.append(name)
			counts[barcode] = counts.get(barcode, 0) + 1
	assert sorted(got) == kept
	check_freq(tmp_path / 'lib.fq.freq.txt', counts)


def test_barcodes_several_files(tmp_path):
	write_fastq(tmp_path / 'one.fq', barcode_reads()[:20])
	write_fastq(tmp_path / 'two.fq', barcode_reads()[20:])
	hh.run(['homerTools', 'barcodes', 4, '-min', 0, 'one.fq', 'two.fq'], tmp_path)
	for name, reads in (('one.fq', barcode_reads()[:20]), ('two.fq', barcode_reads()[20:])):
		counts = {}
		for _, seq, _ in reads:
			counts[seq[:4]] = counts.get(seq[:4], 0) + 1
		assert split_files(tmp_path, name) == sorted(counts)
		check_freq(tmp_path / (name + '.freq.txt'), counts)


TRUSEQ_COUNTS = {'ACGTAC': 12, 'TTGCAA': 5, 'GGCCTA': 2, 'CAGATC': 1}


def write_truseq(path):
	"""Illumina read names ending in the index (1:N:0:<index>)."""

	plan = []
	for barcode, n in sorted(TRUSEQ_COUNTS.items()):
		plan += [barcode] * n
	plan = [plan[(i * 7) % len(plan)] for i in range(len(plan))]
	reads = [('M01234:12:000000000-ABCDE:1:1101:{}:{} 1:N:0:{}'.format(1000 + i, 2000 + 3 * i,
		barcode), INSERT, 'I' * len(INSERT)) for i, barcode in enumerate(plan)]
	assert sorted(set(plan)) == sorted(TRUSEQ_COUNTS)
	return write_fastq(path, reads)


def test_truseq_counts_indexes(tmp_path):
	"""truseq counts the index after the last ':' of each read name into
	<input>.freq.txt and writes no other file."""

	write_truseq(tmp_path / 'idx.fq')
	result = hh.run(['homerTools', 'truseq', 'idx.fq'], tmp_path)
	assert result.stdout == ''
	assert sorted(os.listdir(str(tmp_path))) == ['idx.fq', 'idx.fq.freq.txt']
	check_freq(tmp_path / 'idx.fq.freq.txt', TRUSEQ_COUNTS)


def test_truseq_quality_base_alone_keeps_the_table(tmp_path):
	"""-qualBase only sets the character -qual counts from, so on its own it
	leaves the table as it is."""

	write_truseq(tmp_path / 'idx.fq')
	hh.run(['homerTools', 'truseq', '-qualBase', '!', 'idx.fq'], tmp_path)
	assert sorted(os.listdir(str(tmp_path))) == ['idx.fq', 'idx.fq.freq.txt']
	check_freq(tmp_path / 'idx.fq.freq.txt', TRUSEQ_COUNTS)


def test_truseq_several_files(tmp_path):
	write_truseq(tmp_path / 'a.fq')
	write_truseq(tmp_path / 'b.fq')
	hh.run(['homerTools', 'truseq', 'a.fq', 'b.fq'], tmp_path)
	for name in ('a.fq', 'b.fq'):
		check_freq(tmp_path / (name + '.freq.txt'), TRUSEQ_COUNTS)

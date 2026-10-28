"""homerTools freq: position-dependent nucleotide and dinucleotide
frequencies of a set of sequences, and with -gc the CpG/GC content of each.

Expected tables come from freq_reference.freq_table() and gc_line(), which
count bases in Python: A, C, G and T in either case count, anything else is
skipped; a dinucleotide counts at the position of its first base. The
sequences are seeded random text of several lengths with lowercase stretches
and Ns, written by the tests as tab-separated text, FASTA and FASTQ.
"""

import random

import pytest

import freq_reference as fr
import homer_helpers as hh


def make_records():
	rng = random.Random(7)
	records = []
	for i in range(40):
		length = 60 if i == 0 else rng.choice([30, 45, 60, 60, 75])
		seq = [rng.choice('ACGT' if rng.random() < 0.8 else 'ACGTGC') for _ in range(length)]
		if i % 5 == 1:
			start = rng.randrange(0, length - 10)
			seq[start:start + 10] = [b.lower() for b in seq[start:start + 10]]
		if i % 7 == 2:
			seq[rng.randrange(length)] = 'N'
			seq[rng.randrange(length)] = 'R'
		records.append(('seq{:02d}'.format(i), ''.join(seq)))
	return records


RECORDS = make_records()


def write_fastq(path, records):
	with open(str(path), 'w') as out:
		for name, seq in records:
			out.write('@{}\n{}\n+\n{}\n'.format(name, seq, 'I' * len(seq)))
	return str(path)


@pytest.fixture(scope='module')
def inputs(tmp_path_factory):
	directory = tmp_path_factory.mktemp('freq')
	return {
		'tsv': hh.write_seq_tsv(str(directory / 'seqs.tsv'), RECORDS),
		'fasta': hh.write_fasta(str(directory / 'seqs.fa'), RECORDS, width=25),
		'fastq': write_fastq(directory / 'seqs.fq', RECORDS),
	}


def table(text):
	lines = [line.split('\t') for line in text.splitlines()]
	assert lines[0] == fr.FREQ_HEADER
	return lines[1:]


def expected(length=None, offset=0, records=RECORDS):
	if length is None:
		length = len(records[0][1])
	return fr.freq_table([(seq, 1.0) for _, seq in records], length, offset)


def test_tsv_frequencies(tmp_path, inputs):
	"""Without -maxlen the table has one row per position of the first
	sequence (60), longer sequences counting only up to there."""

	result = hh.run(['homerTools', 'freq', inputs['tsv']], tmp_path)
	got = table(result.stdout)
	assert len(got) == 60
	assert got == expected()
	assert max(len(s) for _, s in RECORDS) > 60


@pytest.mark.parametrize('fmt', ['fasta', 'fastq'])
def test_fasta_and_fastq_are_detected(tmp_path, inputs, fmt):
	result = hh.run(['homerTools', 'freq', inputs[fmt]], tmp_path)
	assert table(result.stdout) == expected()


@pytest.mark.parametrize('fmt', ['tsv', 'fasta', 'fastq'])
def test_format_option(tmp_path, inputs, fmt):
	result = hh.run(['homerTools', 'freq', '-format', fmt, inputs[fmt]], tmp_path)
	assert table(result.stdout) == expected()


@pytest.mark.parametrize('offset', [0, 1, -50])
def test_offset(tmp_path, inputs, offset):
	"""-offset <#> is added to the position label; the frequencies do not
	change."""

	result = hh.run(['homerTools', 'freq', inputs['tsv'], '-offset', offset], tmp_path)
	got = table(result.stdout)
	assert got == expected(offset=offset)
	assert [int(row[0]) for row in got] == list(range(offset, offset + 60))


@pytest.mark.parametrize('maxlen', [1, 50, 75])
def test_maxlen(tmp_path, inputs, maxlen):
	"""-maxlen <#> sets the number of positions: shorter than the first
	sequence cuts the table, longer counts the longer sequences further
	(positions past a sequence's end get nothing from it)."""

	result = hh.run(['homerTools', 'freq', inputs['fasta'], '-maxlen', maxlen], tmp_path)
	got = table(result.stdout)
	assert len(got) == maxlen
	assert got == expected(length=maxlen)


def test_output_file(tmp_path, inputs):
	out = tmp_path / 'freq.txt'
	result = hh.run(['homerTools', 'freq', inputs['tsv'], '-o', out], tmp_path)
	assert result.stdout == ''
	stdout = hh.run(['homerTools', 'freq', inputs['tsv']], tmp_path).stdout
	assert out.read_text() == stdout


@pytest.mark.parametrize('fmt', ['tsv', 'fasta', 'fastq'])
def test_gc_file(tmp_path, inputs, fmt):
	"""-gc <file> writes, per sequence in input order, name, CpG (CG over
	all dinucleotides of two A/C/G/T bases), the GC, AG and AC fractions of
	the A/C/G/T bases and their number; the frequency table still goes to
	stdout."""

	out = tmp_path / 'gc.txt'
	result = hh.run(['homerTools', 'freq', inputs[fmt], '-gc', out], tmp_path)
	assert table(result.stdout) == expected()
	got = [line.split('\t') for line in out.read_text().splitlines()]
	assert got == [fr.gc_line(name, seq) for name, seq in RECORDS]


def test_gc_values_by_hand(tmp_path):
	"""CGCG: 3 dinucleotides, 2 of them CG; all 4 bases G or C, 2 of them A
	or G (the Gs) and 2 A or C. AANNT: CpG 0, GC 0, AG 2/3, AC 2/3, length 3."""

	seqs = hh.write_seq_tsv(str(tmp_path / 's.tsv'), [('a', 'CGCG'), ('b', 'AANNT'),
		('c', 'acgt')])
	out = tmp_path / 'gc.txt'
	hh.run(['homerTools', 'freq', seqs, '-gc', out], tmp_path)
	assert out.read_text().splitlines() == [
		'a\t0.666667\t1.000000\t0.500000\t0.500000\t4',
		'b\t0.000000\t0.000000\t0.666667\t0.666667\t3',
		'c\t0.333333\t0.500000\t0.500000\t0.500000\t4',
	]


def test_frequencies_by_hand(tmp_path):
	"""Two sequences, AC and AG: position 0 is all A with dinucleotides AC and
	AG (A/G column: AG only, half); position 1 has no dinucleotide."""

	seqs = hh.write_seq_tsv(str(tmp_path / 's.tsv'), [('x', 'AC'), ('y', 'AG')])
	got = table(hh.run(['homerTools', 'freq', seqs], tmp_path).stdout)
	half, zero, one = '0.500000', '0.000000', '1.000000'
	assert got[0] == ['0', one, zero, zero, zero, zero, zero, half, zero,
		zero, half, half] + [zero] * 13 + ['2.0', '2.0']
	assert got[1] == ['1', zero, half, half, zero] + ['0.0'] * 20 + ['2.0', '0.0']

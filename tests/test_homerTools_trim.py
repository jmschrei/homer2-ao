"""homerTools trim: trimming FASTQ reads by fixed lengths, adapter
sequences and base quality.

Each test writes a few reads into its own directory (trim writes its output
beside the input: <input>.trimmed, or <input>.<suffix>, and <input>.lengths)
and checks the trimmed reads, worked out by hand from the usage text. Read
bodies are made of C and T only, so the adapters used here (which start with
A or G) can only match where a test puts them. A partial adapter at a
read's 3' end is either shorter than -minMatchLength or longer than it by at
least one base, so the expectation does not depend on whether the read's
last base is compared. Quality cuts are placed where the first window
checked (at -qstart) decides them, or where no window falls below -q.
"""

import os

import pytest

import homer_helpers as hh

ADAPTER = 'AGATCGGAAGAGC'   # 13 bp, so -minMatchLength defaults to 6
BODY = 'CTTCCTCTTTCCTCCTTCTCTTCCCTTCTCTCCTTTCTCC'

GOOD = 'I'   # Phred 40 with offset 33
BAD = '#'    # Phred 2


def body(n, start=0):
	return (BODY * 3)[start:start + n]


def write_fastq(path, reads):
	"""reads: (name, sequence) or (name, sequence, quality); quality
	defaults to all GOOD."""

	with open(str(path), 'w') as out:
		for read in reads:
			name, seq = read[0], read[1]
			qual = read[2] if len(read) > 2 else GOOD * len(seq)
			assert len(qual) == len(seq)
			out.write('@{}\n{}\n+\n{}\n'.format(name, seq, qual))
	return str(path)


def read_fastq(path):
	lines = open(str(path)).read().splitlines()
	assert len(lines) % 4 == 0, lines
	reads = []
	for i in range(0, len(lines), 4):
		assert lines[i].startswith('@') and lines[i + 2] == '+', lines[i:i + 4]
		assert len(lines[i + 1]) == len(lines[i + 3]), lines[i:i + 4]
		reads.append((lines[i][1:], lines[i + 1], lines[i + 3]))
	return reads


def trim(tmp_path, reads, *options, name='reads.fq'):
	"""Run trim in tmp_path on `reads` written to tmp_path/name, the file
	named relative to the working directory; return the trimmed reads as
	(name, sequence, quality)."""

	path = write_fastq(tmp_path / name, reads)
	result = hh.run(['homerTools', 'trim'] + list(options) + [name], tmp_path)
	assert result.stdout == ''
	return read_fastq(path + '.trimmed')


def seqs(reads):
	return [(r[0], r[1]) for r in reads]


def test_no_options_copies_reads(tmp_path):
	reads = [('r1', body(30), GOOD * 30), ('r2', body(12, 3), GOOD * 6 + BAD * 6)]
	assert trim(tmp_path, reads) == reads


def test_fixed_3prime_length(tmp_path):
	"""-3 <#> removes # bases from the 3' end of the sequence and quality; a
	read no longer than # is left empty and dropped (-min 1)."""

	reads = [('long', body(30), GOOD * 25 + BAD * 5), ('short', body(9)),
		('exact', body(5)), ('tiny', body(3))]
	got = trim(tmp_path, reads, '-3', 5)
	assert got == [('long', body(25), GOOD * 25), ('short', body(4), GOOD * 4)]


def test_fixed_5prime_length(tmp_path):
	reads = [('long', body(30), BAD * 4 + GOOD * 26), ('short', body(9)), ('exact', body(4))]
	got = trim(tmp_path, reads, '-5', 4)
	assert got == [('long', body(26, 4), GOOD * 26), ('short', body(5, 4), GOOD * 5)]


def test_3prime_adapter(tmp_path):
	"""-3 <sequence> cuts each read at the first copy of the adapter, and at
	a partial copy running off the 3' end when it has at least
	-minMatchLength bases (default: half the adapter, 6); reads without one
	are left alone."""

	reads = [
		('inside', body(15) + ADAPTER + body(7, 15)),
		('at_start', ADAPTER + body(10)),
		('partial8', body(20) + ADAPTER[:8]),
		('partial3', body(20) + ADAPTER[:3]),
		('none', body(25)),
	]
	got = seqs(trim(tmp_path, reads, '-3', ADAPTER))
	assert got == [('inside', body(15)), ('partial8', body(20)),
		('partial3', body(20) + ADAPTER[:3]), ('none', body(25))]


def test_adapter_trimming_keeps_quality_in_step(tmp_path):
	qual = ''.join(chr(33 + 20 + (i % 15)) for i in range(35))
	reads = [('inside', body(15) + ADAPTER + body(7, 15), qual)]
	assert trim(tmp_path, reads, '-3', ADAPTER) == [('inside', body(15), qual[:15])]


@pytest.mark.parametrize('length, trimmed', [
	(2, ['partial8', 'partial3']),
	(6, ['partial8']),
	(10, []),
])
def test_min_match_length(tmp_path, length, trimmed):
	"""-minMatchLength <#> is the shortest partial adapter at the 3' end
	that is cut; a full copy inside the read is always cut."""

	reads = [
		('inside', body(15) + ADAPTER + body(7, 15)),
		('partial8', body(20) + ADAPTER[:8]),
		('partial3', body(20) + ADAPTER[:3]),
	]
	got = dict(seqs(trim(tmp_path, reads, '-3', ADAPTER, '-minMatchLength', length)))
	assert got['inside'] == body(15)
	for name in ('partial8', 'partial3'):
		assert (got[name] == body(20)) == (name in trimmed), name


def with_mismatches(seq, positions):
	swap = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}
	seq = list(seq)
	for p in positions:
		seq[p] = swap[seq[p]]
	return ''.join(seq)


@pytest.mark.parametrize('mis', [0, 1, 2])
def test_mismatches(tmp_path, mis):
	"""-mis <#> allows up to # mismatches between the adapter and the read."""

	reads = [
		('exact', body(15) + ADAPTER + body(7)),
		('one', body(15) + with_mismatches(ADAPTER, [5]) + body(7)),
		('two', body(15) + with_mismatches(ADAPTER, [3, 9]) + body(7)),
	]
	got = dict(seqs(trim(tmp_path, reads, '-3', ADAPTER, '-mis', mis)))
	for name, needed in (('exact', 0), ('one', 1), ('two', 2)):
		assert (got[name] == body(15)) == (mis >= needed), name


def test_match_start(tmp_path):
	"""-matchStart <#> ignores adapter copies that start before position #,
	as in the usage's Hi-C example (-3 GATC -mis 0 -matchStart 20 -min 20)."""

	read = body(10) + 'GATC' + body(12, 10) + 'GATC' + body(8, 22)
	reads = [('hic', read), ('early', body(10) + 'GATC' + body(20))]
	assert seqs(trim(tmp_path, reads, '-3', 'GATC', '-mis', 0)) == \
		[('hic', body(10)), ('early', body(10))]
	got = seqs(trim(tmp_path, reads, '-3', 'GATC', '-mis', 0, '-matchStart', 20, '-min', 20))
	assert got == [('hic', read[:26]), ('early', body(10) + 'GATC' + body(20))]


def test_5prime_adapter(tmp_path):
	"""-5 <sequence> removes everything up to and including the adapter, and
	a partial copy at the 5' end of at least half the adapter's length."""

	adapter = 'GATCGGA'
	reads = [
		('whole', adapter + body(25)),
		('inside', body(6) + adapter + body(20)),
		('partial', adapter[3:] + body(20)),
		('none', body(24)),
	]
	got = seqs(trim(tmp_path, reads, '-5', adapter))
	assert got == [('whole', body(25)), ('inside', body(20)), ('partial', body(20)),
		('none', body(24))]


def test_5prime_adapter_keeps_quality_in_step(tmp_path):
	qual = ''.join(chr(33 + 20 + i) for i in range(32))
	reads = [('whole', 'GATCGGA' + body(25), qual)]
	assert trim(tmp_path, reads, '-5', 'GATCGGA') == [('whole', body(25), qual[7:])]


def test_n_matches_anything_unless_explicit(tmp_path):
	"""An N in the read or in the adapter matches any base; with -explicitN
	it matches only N."""

	reads = [
		('n_in_read', body(15) + ADAPTER[:4] + 'N' + ADAPTER[5:] + body(6)),
		('plain', body(15) + ADAPTER + body(6)),
	]
	got = dict(seqs(trim(tmp_path, reads, '-3', ADAPTER)))
	assert got == {'n_in_read': body(15), 'plain': body(15)}
	got = dict(seqs(trim(tmp_path, reads, '-3', ADAPTER, '-explicitN')))
	assert got == {'n_in_read': reads[0][1], 'plain': body(15)}
	n_adapter = ADAPTER[:6] + 'N' + ADAPTER[7:]
	got = dict(seqs(trim(tmp_path, reads, '-3', n_adapter)))
	assert got == {'n_in_read': body(15), 'plain': body(15)}
	got = dict(seqs(trim(tmp_path, reads, '-3', n_adapter, '-explicitN')))
	assert got == {'n_in_read': reads[0][1], 'plain': reads[1][1]}


def test_fixed_length(tmp_path):
	"""-len <#> keeps the first # bases of longer reads."""

	reads = [('long', body(30), GOOD * 15 + BAD * 15), ('short', body(10)), ('exact', body(15))]
	got = trim(tmp_path, reads, '-len', 15)
	assert got == [('long', body(15), GOOD * 15), ('short', body(10), GOOD * 10),
		('exact', body(15), GOOD * 15)]


def test_min_length(tmp_path):
	"""-min <#> drops reads shorter than # after trimming."""

	reads = [('a', body(20)), ('b', body(12)), ('c', body(11)),
		('d', body(10) + ADAPTER + body(5)), ('e', body(14) + ADAPTER)]
	assert [r[0] for r in trim(tmp_path, reads, '-min', 12)] == ['a', 'b', 'd', 'e']
	got = seqs(trim(tmp_path, reads, '-3', ADAPTER, '-min', 12))
	assert got == [('a', body(20)), ('b', body(12)), ('e', body(14))]


def test_max_length(tmp_path):
	"""-max <#> drops reads longer than #."""

	reads = [('a', body(20)), ('b', body(30)), ('c', body(10))]
	assert [r[0] for r in trim(tmp_path, reads, '-max', 25)] == ['a', 'c']
	assert [r[0] for r in trim(tmp_path, reads, '-max', 25, '-3', 10)] == ['a', 'b']


def quality_reads():
	read = body(30)
	return [
		('drop_at_10', read, GOOD * 10 + BAD * 20),
		('all_good', read, GOOD * 30),
		('dip_at_10', read, GOOD * 10 + BAD * 2 + GOOD * 18),
	]


@pytest.mark.parametrize('options, lengths', [
	# the 5 bp window starting at position 10 (the default -qstart) averages
	# (40 * 3 + 2 * 2) / 5 = 24 over the two-base dip
	(['-q', 20], [10, 30, 30]),
	(['-q', 20, '-qwindow', 1], [10, 30, 10]),
	(['-q', 20, '-qstart', 12], [12, 30, 30]),
	(['-q', 40], [10, 30, 10]),
	(['-q', 0], [30, 30, 30]),
])
def test_quality(tmp_path, options, lengths):
	"""-q <#> cuts a read at the first position, from -qstart on, where the
	average Phred quality (offset 33, integer division) of the -qwindow bases
	starting there is below #; the quality string is cut with it. Quality 40
	is not below -q 40."""

	got = trim(tmp_path, quality_reads(), *options)
	want = [(name, seq[:n], qual[:n]) for (name, seq, qual), n in zip(quality_reads(), lengths)]
	assert got == want


def test_quality_start(tmp_path):
	"""A read whose quality drops at base 12 is cut at 10 by default, where
	the window 10-14 averages (40 * 2 + 2 * 3) / 5 = 17, and at 12 with
	-qstart 12."""

	reads = [('drop_at_12', body(30), GOOD * 12 + BAD * 18)]
	assert seqs(trim(tmp_path, reads, '-q', 20)) == [('drop_at_12', body(10))]
	assert seqs(trim(tmp_path, reads, '-q', 20, '-qstart', 12)) == [('drop_at_12', body(12))]


def test_quality_by_default_is_off(tmp_path):
	assert trim(tmp_path, quality_reads()) == quality_reads()


def test_suffix(tmp_path):
	path = write_fastq(tmp_path / 'r.fq', [('a', body(20))])
	hh.run(['homerTools', 'trim', '-3', 5, '-suffix', 'cut', 'r.fq'], tmp_path)
	assert read_fastq(path + '.cut') == [('a', body(15), GOOD * 15)]
	assert not os.path.exists(path + '.trimmed')


def test_length_distribution(tmp_path):
	"""<input>.lengths lists, for each length from 1 to the longest kept
	read, the number of kept reads of that length and their percentage of all
	reads. Row 0 (reads dropped) is not checked."""

	reads = [('a', body(20)), ('b', body(20)), ('c', body(17)), ('d', body(10) + ADAPTER),
		('e', body(4))]
	path = write_fastq(tmp_path / 'r.fq', reads)
	hh.run(['homerTools', 'trim', '-3', ADAPTER, '-min', 5, 'r.fq'], tmp_path)
	lines = [line.split('\t') for line in open(path + '.lengths').read().splitlines()]
	assert lines[0] == ['Length', '# reads', 'Fraction']
	table = {int(row[0]): row[1:] for row in lines[1:] if int(row[0]) > 0}
	assert sorted(table) == list(range(1, 21))
	counts = {10: 1, 17: 1, 20: 2}
	for length, row in table.items():
		n = counts.get(length, 0)
		assert row == [str(n), '%f%%' % (n / 5.0 * 100)], (length, row)


def test_unique_names(tmp_path):
	"""-uniq names the reads Read1, Read2, ... in input order."""

	reads = [('x y z', body(20)), ('x y z', body(21)), ('other', body(22))]
	got = trim(tmp_path, reads, '-uniq')
	assert seqs(got) == [('Read1', body(20)), ('Read2', body(21)), ('Read3', body(22))]


def test_several_files(tmp_path):
	one = write_fastq(tmp_path / 'one.fq', [('a', body(20))])
	two = write_fastq(tmp_path / 'two.fq', [('b', body(24)), ('c', body(12))])
	hh.run(['homerTools', 'trim', '-len', 10, 'one.fq', 'two.fq'], tmp_path)
	assert seqs(read_fastq(one + '.trimmed')) == [('a', body(10))]
	assert seqs(read_fastq(two + '.trimmed')) == [('b', body(10)), ('c', body(10))]


MATES = [
	('p1', body(30), body(28, 2)),
	('p2', body(15) + ADAPTER + body(5), body(33, 1)),
	('p3', body(30), body(8) + ADAPTER + body(12)),
	('p4', body(26), body(26, 4)),
]


def write_mates(directory, first, second):
	write_fastq(directory / first, [(n, a) for n, a, _ in MATES])
	write_fastq(directory / second, [(n, b) for n, _, b in MATES])
	return str(directory / first), str(directory / second)


@pytest.mark.parametrize('first, second', [
	('lib_R1_001.fastq', 'lib_R2_001.fastq'),
	('lib_R1.fq', 'lib_R2.fq'),
	('lib_1.fastq', 'lib_2.fastq'),
])
def test_paired_end_mate_detection(tmp_path, first, second):
	"""-pe finds the mate file (_R1_ / _R1. / _1.f renamed to _R2_ / _R2. /
	_2.f) and trims both mates; a pair is kept only when both mates are at
	least -min long, in both outputs."""

	r1, r2 = write_mates(tmp_path, first, second)
	hh.run(['homerTools', 'trim', '-pe', '-3', ADAPTER, '-min', 12, first], tmp_path)
	got1, got2 = read_fastq(r1 + '.trimmed'), read_fastq(r2 + '.trimmed')
	assert seqs(got1) == [('p1', body(30)), ('p2', body(15)), ('p4', body(26))]
	assert seqs(got2) == [('p1', body(28, 2)), ('p2', body(33, 1)), ('p4', body(26, 4))]


def test_paired_end_files_named_together(tmp_path):
	"""<R1>,<R2> on the command line names the mate file directly."""

	r1, r2 = write_mates(tmp_path, 'left.fq', 'right.fq')
	hh.run(['homerTools', 'trim', '-len', 20, '-uniq', 'left.fq,right.fq'], tmp_path)
	got1, got2 = read_fastq(r1 + '.trimmed'), read_fastq(r2 + '.trimmed')
	assert seqs(got1) == [('Read{}'.format(i + 1), a[:20]) for i, (_, a, _) in enumerate(MATES)]
	assert seqs(got2) == [('Read{}'.format(i + 1), b[:20]) for i, (_, _, b) in enumerate(MATES)]


def test_single_end_without_mate_file(tmp_path):
	"""-pe with no matching mate file trims the file as single-end."""

	path = write_fastq(tmp_path / 'solo_R1_001.fq', [('a', body(20))])
	hh.run(['homerTools', 'trim', '-pe', '-len', 10, 'solo_R1_001.fq'], tmp_path)
	assert seqs(read_fastq(path + '.trimmed')) == [('a', body(10))]
	assert sorted(os.listdir(str(tmp_path))) == ['solo_R1_001.fq', 'solo_R1_001.fq.lengths',
		'solo_R1_001.fq.trimmed']

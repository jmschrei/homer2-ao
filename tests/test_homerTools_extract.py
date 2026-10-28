"""homerTools extract: the genome sequence under each peak, and
`homerTools extract stats`, the length, N% and GC% of each chromosome.

Expected sequences are Python slices of tests/data/genome.fa: a HOMER peak
file's start and end are 1-based and both inclusive, a BED file's start is
0-based and its end exclusive, a - strand peak gives the reverse complement,
lowercase bases are printed in upper case (or as N with -mask) and anything
other than A, C, G or T becomes N. The statistics are counted in Python from
the same file.
"""

import os
import shutil

import pytest

import homer_helpers as hh

EDGE_PEAKS = [
	# name, chromosome, start, end (1-based, inclusive), strand
	('first_base', 'chr1', 1, 50, '+'),
	('last_base', 'chr3', 59951, 60000, '+'),
	('last_base_minus', 'chr3', 59901, 60000, '-'),
	('overlap_a', 'chr1', 5000, 5300, '+'),
	('overlap_b', 'chr1', 5100, 5150, '-'),
	('overlap_c', 'chr1', 5200, 5600, '+'),
	('wide', 'chr1', 4000, 9000, '+'),
	('into_n_run', 'chr2', 39901, 40100, '+'),
	('same_plus', 'chr2', 1000, 1100, '+'),
	('same_minus', 'chr2', 1000, 1100, '-'),
	('single_base', 'chr2', 777, 777, '-'),
]


def expected_sequence(genome, chrom, start, end, strand, mask=False):
	seq = genome[chrom][start - 1:end]
	if mask:
		seq = ''.join('N' if b.islower() else b for b in seq)
	seq = hh.clean_sequence(seq)
	return hh.revcomp(seq) if strand == '-' else seq


def tsv_records(text):
	"""{name: sequence}, failing on a repeated name."""

	records = {}
	for line in text.splitlines():
		name, seq = line.split('\t')
		assert name not in records, name
		records[name] = seq
	return records


def fasta_records(path):
	records = {}
	for name, seq in hh.read_fasta(path):
		assert name not in records, name
		records[name] = seq
	return records


@pytest.fixture(scope='module')
def all_peaks(genome, tmp_path_factory):
	"""tests/data/peaks.txt plus EDGE_PEAKS, as (path, {name: expected
	sequence}, {name: expected sequence with -mask})."""

	directory = tmp_path_factory.mktemp('extract_peaks')
	rows = [(p['id'], p['chr'], p['start'], p['end'], p['strand']) for p in hh.read_peaks()]
	rows += EDGE_PEAKS
	path = directory / 'peaks.txt'
	with open(str(path), 'w') as out:
		out.write('#PeakID\tchr\tstart\tend\tstrand\n')
		for row in rows:
			out.write('{}\t{}\t{}\t{}\t{}\n'.format(*row))
	plain = {r[0]: expected_sequence(genome, *r[1:]) for r in rows}
	masked = {r[0]: expected_sequence(genome, *r[1:], mask=True) for r in rows}
	return str(path), plain, masked


@pytest.fixture(scope='module')
def genome_layouts(tmp_path_factory):
	"""The genome as one FASTA file, as a directory of per-chromosome .fa
	files (50 bases per line, a description after the name), as a directory
	of per-chromosome .fa.masked files (80 per line), and as a directory
	holding genome.fa."""

	directory = tmp_path_factory.mktemp('extract_genomes')
	records = hh.read_fasta(hh.GENOME)
	layouts = {'file': hh.GENOME}
	for key, suffix, width in (('chroms', '.fa', 50), ('masked', '.fa.masked', 80)):
		os.makedirs(str(directory / key))
		for name, seq in records:
			hh.write_fasta(str(directory / key / (name + suffix)),
				[(name + ' assembled chromosome', seq)], width=width)
		layouts[key] = str(directory / key)
	os.makedirs(str(directory / 'whole'))
	shutil.copy(hh.GENOME, str(directory / 'whole' / 'genome.fa'))
	layouts['whole'] = str(directory / 'whole')
	return layouts


def test_peak_file(tmp_path, all_peaks):
	"""Every peak of a HOMER peak file is printed once, tab-separated, with
	its sequence: + strand peaks as the genome reads, - strand peaks reverse
	complemented, all upper case; peaks that overlap, sit at a chromosome's
	first or last base, cover one base or run into Ns are cut the same way."""

	path, plain, _ = all_peaks
	result = hh.run(['homerTools', 'extract', path, hh.GENOME], tmp_path)
	assert tsv_records(result.stdout) == plain


def test_bed_file_uses_zero_based_starts(tmp_path, genome, peaks):
	"""tests/data/peaks.bed holds the same 200 peaks as peaks.txt, its starts
	one less; both give the same sequences, the BED file's - strand column
	included."""

	result = hh.run(['homerTools', 'extract', hh.PEAKS_BED, hh.GENOME], tmp_path)
	want = {p['id']: hh.peak_sequence(p, genome).upper() for p in peaks}
	assert sum(p['strand'] == '-' for p in peaks) > 50
	assert tsv_records(result.stdout) == want


def test_bed_without_name_or_strand(tmp_path, genome):
	"""A three-column BED file gives + strand sequences, one per line."""

	regions = [('chr1', 99, 150), ('chr2', 2000, 2060), ('chr3', 0, 25)]
	bed = tmp_path / 'three.bed'
	bed.write_text(''.join('{}\t{}\t{}\n'.format(*r) for r in regions))
	result = hh.run(['homerTools', 'extract', bed, hh.GENOME], tmp_path)
	got = sorted(tsv_records(result.stdout).values())
	assert got == sorted(genome[c][s:e].upper() for c, s, e in regions)


def test_bed_strand_column(tmp_path, genome):
	"""In BED6 the sixth column is the strand."""

	bed = tmp_path / 'six.bed'
	bed.write_text('chr1\t1000\t1040\tfwd\t0\t+\nchr1\t1000\t1040\trev\t0\t-\n')
	got = tsv_records(hh.run(['homerTools', 'extract', bed, hh.GENOME], tmp_path).stdout)
	assert got == {'fwd': genome['chr1'][1000:1040].upper(),
		'rev': hh.revcomp(genome['chr1'][1000:1040].upper())}


def test_fasta_output(tmp_path, all_peaks):
	"""-fa prints the same records as FASTA."""

	path, plain, _ = all_peaks
	result = hh.run(['homerTools', 'extract', path, hh.GENOME, '-fa'], tmp_path)
	out = tmp_path / 'out.fa'
	out.write_text(result.stdout)
	assert result.stdout.startswith('>')
	assert fasta_records(str(out)) == plain


def test_mask_turns_lowercase_into_n(tmp_path, all_peaks):
	"""-mask prints soft-masked (lowercase) bases as N; the planted
	lowercase peaks become all N."""

	path, plain, masked = all_peaks
	result = hh.run(['homerTools', 'extract', path, hh.GENOME, '-mask'], tmp_path)
	got = tsv_records(result.stdout)
	assert got == masked
	lowercase = {p['peak'] for p in hh.read_planted() if p['lowercase'] == 'yes'}
	assert len(lowercase) == 10
	for name in lowercase:
		assert set(got[name]) == {'N'} and set(plain[name]) != {'N'}


def test_mask_and_fasta_together(tmp_path, all_peaks):
	path, _, masked = all_peaks
	result = hh.run(['homerTools', 'extract', path, hh.GENOME, '-mask', '-fa'], tmp_path)
	out = tmp_path / 'out.fa'
	out.write_text(result.stdout)
	assert fasta_records(str(out)) == masked


@pytest.mark.parametrize('layout', ['chroms', 'masked', 'whole'])
def test_genome_directory(tmp_path, all_peaks, genome_layouts, layout):
	"""A directory of per-chromosome <chr>.fa or <chr>.fa.masked files, or one
	holding genome.fa, gives the same sequences as the single FASTA file,
	whatever the line width."""

	path, plain, masked = all_peaks
	directory = genome_layouts[layout]
	got = tsv_records(hh.run(['homerTools', 'extract', path, directory], tmp_path).stdout)
	assert got == plain
	got = tsv_records(hh.run(['homerTools', 'extract', path, directory, '-mask'],
		tmp_path).stdout)
	assert got == masked


@pytest.mark.parametrize('layout', ['file', 'chroms'])
def test_peaks_on_missing_chromosomes_are_skipped(tmp_path, genome, genome_layouts, layout):
	peaks = tmp_path / 'peaks.txt'
	peaks.write_text('kept\tchr1\t101\t130\t+\ngone\tchrUn\t101\t130\t+\n')
	got = tsv_records(hh.run(['homerTools', 'extract', peaks, genome_layouts[layout]],
		tmp_path).stdout)
	assert got == {'kept': genome['chr1'][100:130].upper()}


def test_other_characters_become_n(tmp_path):
	"""IUPAC codes and other letters in the genome are printed as N; the
	reverse complement of N is N."""

	fasta = hh.write_fasta(str(tmp_path / 'odd.fa'), [('chrA', 'ACGTRYKMacgtnNSW' * 3)])
	peaks = tmp_path / 'peaks.txt'
	peaks.write_text('p\tchrA\t1\t16\t+\nm\tchrA\t17\t32\t-\n')
	got = tsv_records(hh.run(['homerTools', 'extract', peaks, fasta], tmp_path).stdout)
	assert got == {'p': 'ACGTNNNNACGTNNNN', 'm': hh.revcomp('ACGTNNNNACGTNNNN')}


def stats_rows(genome):
	"""{name: [length, N%, GC%]} per chromosome and for 'genome': N and GC
	(either case) counted over all bases, N included."""

	rows = {}
	totals = [0, 0, 0]
	for name, seq in genome.items():
		upper = seq.upper()
		counts = [len(seq), upper.count('N'), upper.count('G') + upper.count('C')]
		totals = [a + b for a, b in zip(totals, counts)]
		rows[name] = counts
	rows['genome'] = totals

	def text(length, n, gc):
		return [str(length), '%.2f%%' % (n / length * 100.0), '%.2f%%' % (gc / length * 100.0)]

	return {name: text(*counts) for name, counts in rows.items()}


@pytest.mark.parametrize('layout', ['file', 'chroms', 'masked', 'whole'])
def test_stats(tmp_path, genome, genome_layouts, layout):
	"""`extract stats` prints a header, one row per chromosome (any order for a
	directory of files) and a last 'genome' row with the totals."""

	result = hh.run(['homerTools', 'extract', 'stats', genome_layouts[layout]], tmp_path)
	lines = [line.split('\t') for line in result.stdout.splitlines()]
	assert lines[0] == ['ChrName', 'Length', 'N%', 'GC%']
	assert lines[-1][0] == 'genome'
	got = {row[0]: row[1:] for row in lines[1:]}
	assert len(got) == len(lines) - 1
	want = stats_rows(genome)
	assert got == want
	assert want['chr2'][1] == '2.22%'
	if layout == 'file':
		assert [row[0] for row in lines[1:-1]] == ['chr1', 'chr2', 'chr3']


def test_stats_of_several_records_in_one_file(tmp_path):
	"""Each record of a multi-record file is its own row; lowercase bases
	count like uppercase ones."""

	records = [('chrA', 'ACGTNNnngc'), ('chrB', 'aaaa'), ('chrC', 'GGGGCCCCAT')]
	fasta = hh.write_fasta(str(tmp_path / 'small.fa'), records, width=3)
	result = hh.run(['homerTools', 'extract', 'stats', fasta], tmp_path)
	assert result.stdout.splitlines() == [
		'ChrName\tLength\tN%\tGC%',
		'chrA\t10\t40.00%\t40.00%',
		'chrB\t4\t0.00%\t0.00%',
		'chrC\t10\t0.00%\t80.00%',
		'genome\t24\t16.67%\t50.00%',
	]

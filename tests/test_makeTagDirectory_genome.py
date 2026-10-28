"""makeTagDirectory with -genome: the GC content and nucleotide frequencies
around tags (-checkGC, -freqStart, -freqEnd), the 5' oligo enrichment of tags
(-oligoStart, -oligoEnd) and the normalization of tag counts by it
(-normFixedOligo).

The reads are BED records the module draws from a seeded random generator
on tests/data/genome.fa's chromosomes, on both strands, some repeated, some
next to the N run. Expected files are computed in Python from the genome:

	a tag sits at the read's 5' end (BED start + 1 on +, BED end on -) and
	holds the number of reads of one length there (reads of other lengths
	at the same place are separate tags);
	with fragment length L (fixed with -fragLength 100), a tag's window runs
	from freqStart to L + freqEnd bp downstream of it (upstream for - strand
	tags, reverse complemented); a window off the chromosome is skipped;
	its GC content is taken over the L + 1 bases from the tag on (A, C, G and
	T only), binned at 1/L (rounded down after adding 0.01) and weighted by
	the tag's reads (tagGCcontent.txt); genomeGCcontent.txt bins the
	(L + 1)-base windows of the chromosomes that have tags, the same way;
	tagFreq.txt and tagFreqUniq.txt are the position frequency tables of the
	windows, weighted by reads or once per tag;
	the oligo of a tag is the genome from oligoStart to oligoEnd around it
	(reverse complemented on -), counted once per tag and skipped
	when it holds an N; the genome count is the windows of that length
	without N; enrichment is reads / (genome count * total reads / total
	genome windows), each raised to at least 0.5.

HOMER counts each chromosome's genome windows up to, but not including, the
one that ends at the chromosome's last base. Whether that window should be
left out is not documented and looks like an off-by-one, so the genome GC,
oligo enrichment and -normFixedOligo tests accept the expectation computed
either with or without those last windows.
"""

import os
import random
import struct

import pytest

import freq_reference as fr
import homer_helpers as hh

FRAG = 100
EDGE = 400   # reads stay this far from chromosome ends


def make_reads(genome):
	"""(chromosome, 0-based BED start, end, strand) for 700 reads, 40 of them
	repeated, 10 around each end of the N run on chr2."""

	rng = random.Random(4)
	reads = []
	names = sorted(genome)
	for _ in range(640):
		chrom = rng.choice(names)
		length = rng.choice([30, 36, 50])
		start = rng.randrange(EDGE, len(genome[chrom]) - EDGE)
		reads.append((chrom, start, start + length, rng.choice('+-')))
	for edge in (40000, 42000):
		for _ in range(10):
			start = edge + rng.randrange(-150, 150)
			reads.append(('chr2', start, start + 36, rng.choice('+-')))
	reads += [reads[rng.randrange(len(reads))] for _ in range(40)]
	return reads


def write_bed(path, reads):
	with open(str(path), 'w') as out:
		for i, (chrom, start, end, strand) in enumerate(reads):
			out.write('{}\t{}\t{}\tread{}\t0\t{}\n'.format(chrom, start, end, i, strand))
	return str(path)


def tag_counts(reads):
	"""{(chromosome, position, strand, length): reads}"""

	tags = {}
	for chrom, start, end, strand in reads:
		key = (chrom, start + 1 if strand == '+' else end, strand, end - start)
		tags[key] = tags.get(key, 0) + 1
	return tags


class Reference(object):
	"""Expected file contents for the module's reads, from the genome in
	upper case."""

	def __init__(self, genome, reads):
		self.genome = {c: s.upper() for c, s in genome.items()}
		self.tags = tag_counts(reads)
		self.chroms = sorted({k[0] for k in self.tags})

	def window(self, chrom, start, end, strand):
		"""The genome from start to end (1-based, inclusive), '' unless the
		window starts at base 1 or later and ends before the last base."""

		seq = self.genome[chrom]
		if start < 1 or end >= len(seq):
			return ''
		s = seq[start - 1:end]
		return hh.revcomp(s) if strand == '-' else s

	def fragment(self, chrom, p, strand, freq_start=-50, freq_end=50):
		if strand == '+':
			return self.window(chrom, p + freq_start, p + FRAG + freq_end, strand)
		return self.window(chrom, p - FRAG - freq_end, p - freq_start, strand)

	def oligo(self, chrom, p, strand, start, end):
		if strand == '+':
			return self.window(chrom, p + start, p + end, strand)
		return self.window(chrom, p - end, p - start, strand)

	def tag_gc(self, freq_start=-50, freq_end=50):
		inc, nbins = fr.gc_bins(FRAG)
		counts = [0.0] * nbins
		for (chrom, p, strand, _), v in sorted(self.tags.items()):
			w = self.fragment(chrom, p, strand, freq_start, freq_end)
			part = w[-freq_start:-freq_start + FRAG + 1]
			n = sum(1 for b in part if b in 'ACGT')
			if w and n:
				counts[fr.gc_bin(sum(1 for b in part if b in 'GC'), n, inc)] += v
		return fr.gc_histogram_rows(counts, inc)

	def genome_gc(self, last_window):
		"""genomeGCcontent.txt rows, with or without each chromosome's
		last window."""

		inc, nbins = fr.gc_bins(FRAG)
		counts = [0.0] * nbins
		size = FRAG + 1
		for chrom in self.chroms:
			seq = self.genome[chrom]
			valid = [1 if b in 'ACGT' else 0 for b in seq]
			gc = [1 if b in 'GC' else 0 for b in seq]
			n, g = sum(valid[:size]), sum(gc[:size])
			stop = len(seq) - size + (1 if last_window else 0)
			for i in range(stop):
				if n > 0:
					counts[fr.gc_bin(g, n, inc)] += 1
				if i + size < len(seq):
					n += valid[i + size] - valid[i]
					g += gc[i + size] - gc[i]
		return fr.gc_histogram_rows(counts, inc)

	def freq(self, unique, freq_start=-50, freq_end=50):
		records = []
		for (chrom, p, strand, _), v in sorted(self.tags.items()):
			w = self.fragment(chrom, p, strand, freq_start, freq_end)
			if w:
				records.append((w, 1.0 if unique else float(v)))
		return fr.freq_table(records, freq_end + FRAG - freq_start, freq_start)

	def enrichment(self, start, end, last_window):
		"""{oligo: (reads, genome count, enrichment)} as numbers, the genome
		counted with or without each chromosome's last window."""

		length = end - start + 1
		reads = {}
		for chrom, p, strand, _ in sorted(self.tags):
			o = self.oligo(chrom, p, strand, start, end)
			if o and 'N' not in o:
				reads[o] = reads.get(o, 0) + 1
		background = {}
		for chrom in self.chroms:
			seq = self.genome[chrom]
			for i in range(len(seq) - length + (1 if last_window else 0)):
				o = seq[i:i + length]
				if 'N' not in o:
					background[o] = background.get(o, 0) + 1
		ratio = sum(reads.values()) / float(sum(background.values()))
		table = {}
		for o in set(reads) | set(background):
			nt, nb = reads.get(o, 0), background.get(o, 0)
			table[o] = (nt, nb, max(nt, 0.5) / max(nb * ratio, 0.5))
		return table


@pytest.fixture(scope='module')
def ref(genome):
	return Reference(genome, make_reads(genome))


@pytest.fixture(scope='module')
def runs(genome, tmp_path_factory):
	"""Every makeTagDirectory call of the module, run once: {key: tag
	directory}."""

	directory = tmp_path_factory.mktemp('mtd_genome')
	bed = write_bed(directory / 'reads.bed', make_reads(genome))
	records = hh.read_fasta(hh.GENOME)
	os.makedirs(str(directory / 'chroms'))
	for name, seq in records:
		hh.write_fasta(str(directory / 'chroms' / (name + '.fa')), [(name, seq)], width=70)
	os.makedirs(str(directory / 'whole'))
	hh.write_fasta(str(directory / 'whole' / 'genome.fa'), records, width=80)
	fixed = ['-fragLength', FRAG]
	calls = {
		'gc': ['-genome', hh.GENOME, '-checkGC'],
		'freq_range': ['-genome', hh.GENOME, '-checkGC', '-freqStart', -20, '-freqEnd', 10],
		'oligo6': ['-genome', hh.GENOME, '-oligoStart', -2, '-oligoEnd', 3],
		'oligo1': ['-genome', hh.GENOME, '-oligoStart', 0, '-oligoEnd', 0],
		'oligo13': ['-genome', hh.GENOME, '-oligoStart', -6, '-oligoEnd', 6],
		'chroms': ['-genome', 'chroms', '-checkGC', '-oligoStart', -2, '-oligoEnd', 3],
		'whole': ['-genome', 'whole', '-checkGC', '-oligoStart', -2, '-oligoEnd', 3],
		'plain': [],
		'norm_default': ['-genome', hh.GENOME, '-oligoStart', 0, '-oligoEnd', 1,
			'-normFixedOligo', 'default'],
	}
	hh.run_many({key: (['makeTagDirectory', key, bed] + options + fixed, directory)
		for key, options in calls.items()})
	table = directory / 'norm_default' / 'oligoEnrichment0to1.txt'
	hh.run(['makeTagDirectory', 'norm_file', bed, '-genome', hh.GENOME, '-oligoStart', 0,
		'-oligoEnd', 1, '-normFixedOligo', table] + fixed, directory)
	return {key: directory / key for key in list(calls) + ['norm_file']}


def rows(path):
	return [line.split('\t') for line in open(str(path)).read().splitlines()]


def test_tag_gc_content(runs, ref):
	got = rows(runs['gc'] / 'tagGCcontent.txt')
	assert got[0] == ['GC%', 'Total', 'Normalized Fraction(PDF)']
	assert got[1:] == ref.tag_gc()
	assert sum(float(r[1]) for r in got[1:]) > 600


def test_genome_gc_content(runs, ref):
	got = rows(runs['gc'] / 'genomeGCcontent.txt')
	assert got[0] == ['GC%', 'Total', 'Normalized Fraction(PDF)']
	assert got[1:] in (ref.genome_gc(False), ref.genome_gc(True))


@pytest.mark.parametrize('name, unique', [('tagFreq.txt', False), ('tagFreqUniq.txt', True)])
def test_tag_frequencies(runs, ref, name, unique):
	"""Frequencies from 50 bp before the tag to 50 bp past the fragment
	(-freqStart -50, -freqEnd 50 by default)."""

	got = rows(runs['gc'] / name)
	assert got[0] == fr.FREQ_HEADER
	assert got[1:] == ref.freq(unique)
	assert got[1][0] == '-50' and len(got) - 1 == FRAG + 100


def test_freq_start_and_end(runs, ref):
	"""-freqStart -20 -freqEnd 10 moves the frequency window; the GC content
	is still taken over the fragment."""

	got = rows(runs['freq_range'] / 'tagFreq.txt')
	assert got[1:] == ref.freq(False, -20, 10)
	assert got[1][0] == '-20' and got[-1][0] == str(FRAG + 10 - 1)
	assert rows(runs['freq_range'] / 'tagGCcontent.txt')[1:] == ref.tag_gc(-20, 10)


def test_no_gc_files_without_genome_options(runs):
	files = set(os.listdir(str(runs['plain'])))
	assert 'tagInfo.txt' in files
	assert not files & {'tagGCcontent.txt', 'genomeGCcontent.txt', 'tagFreq.txt',
		'tagFreqUniq.txt'}
	assert not [f for f in files if f.startswith('oligoEnrichment')]


def check_enrichment(path, tables, length):
	"""The file matches one of `tables`."""

	lines = rows(path)
	start, end = path.name[len('oligoEnrichment'):-len('.txt')].split('to')
	assert lines[0] == ['Oligo', 'Read Count({} to {})'.format(start, end), 'Genome Count',
		'Enrichment']
	got = {r[0]: tuple(r[1:]) for r in lines[1:]}
	assert len(got) == len(lines) - 1
	assert got in [{o: ('%.1f' % nt, '%d' % nb, '%.4f' % e) for o, (nt, nb, e) in table.items()}
		for table in tables]
	assert all(len(o) == length for o in got)
	values = [float(r[3]) for r in lines[1:]]
	assert values == sorted(values, reverse=True)


@pytest.mark.parametrize('key, start, end', [
	('oligo6', -2, 3),
	('oligo1', 0, 0),
	('oligo13', -6, 6),
])
def test_oligo_enrichment(runs, ref, key, start, end):
	"""-oligoStart/-oligoEnd write oligoEnrichment<start>to<end>.txt: the
	oligos of that length seen at a tag or in the genome, strands kept apart,
	most enriched first."""

	path = runs[key] / 'oligoEnrichment{}to{}.txt'.format(start, end)
	tables = [ref.enrichment(start, end, last) for last in (False, True)]
	check_enrichment(path, tables, end - start + 1)
	if end - start + 1 <= 6:
		assert len(tables[0]) == 4 ** (end - start + 1)


@pytest.mark.parametrize('layout', ['chroms', 'whole'])
def test_genome_directory(runs, layout):
	"""A directory of per-chromosome FASTA files, or one holding genome.fa,
	gives the same files as the single FASTA file."""

	for name in ('tagGCcontent.txt', 'genomeGCcontent.txt', 'tagFreq.txt', 'tagFreqUniq.txt'):
		assert rows(runs[layout] / name) == rows(runs['gc'] / name), name
	got = rows(runs[layout] / 'oligoEnrichment-2to3.txt')
	want = rows(runs['oligo6'] / 'oligoEnrichment-2to3.txt')
	assert got[0] == want[0] and sorted(got[1:]) == sorted(want[1:])


def read_tags(directory):
	tags = {}
	for name in os.listdir(str(directory)):
		if name.endswith('.tags.tsv'):
			for line in open(os.path.join(str(directory), name)):
				cols = line.rstrip('\n').split('\t')
				key = (cols[1], int(cols[2]), '+' if cols[3] == '0' else '-', int(cols[5]))
				assert key not in tags
				tags[key] = cols[4]
	return tags


def float32(x):
	return struct.unpack('f', struct.pack('f', x))[0]


def test_tags_without_normalization(runs, ref):
	assert read_tags(runs['gc']) == {k: '%.1f' % v for k, v in ref.tags.items()}


@pytest.mark.parametrize('key', ['norm_default', 'norm_file'])
def test_norm_fixed_oligo(runs, ref, key):
	"""-normFixedOligo default multiplies each tag's reads by 1 / the
	enrichment of its oligo, kept between 0.25 and 2 (the default
	-minNormRatio and -maxNormRatio); a tag whose oligo holds an N is left
	alone. An oligoEnrichment file from an earlier run with the same reads
	gives the same values."""

	wants = []
	for last in (False, True):
		table = ref.enrichment(0, 1, last)
		want = {}
		for tag, v in ref.tags.items():
			o = ref.oligo(tag[0], tag[1], tag[2], 0, 1)
			value = float(v)
			if o and 'N' not in o:
				ratio = table[o][2]
				ratio = 1.0 / ratio if ratio > 0.0000001 else ratio
				value = float32(v * min(max(ratio, 0.25), 2.0))
			want[tag] = '%.1f' % value
		wants.append(want)
	got = read_tags(runs[key])
	assert got in wants
	assert len({v for v in got.values()}) > 3

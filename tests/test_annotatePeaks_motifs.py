"""annotatePeaks.pl <peaks> <genome FASTA> -noann -nogene: the motif options.

Gene annotation needs organism data HOMER does not ship, so every call skips
it with -noann -nogene. Each call runs once, in its own directory, in a
module-scoped batch; the tests read its stdout and output files.

Expected values come from a Python model of what the script does, written
from bin/annotatePeaks.pl, bin/adjustPeakFile.pl and homerTools:

Windows. Without -size (and without -hist) a peak is used as given. -size N
	puts N + 1 bases around the peak's centre (adjustPeakFile.pl): the centre
	is floor((start + end) / 2) on the + strand and that rounded up on the -
	strand; a + strand window starts floor(N / 2) before the centre, a - strand
	window ends floor(N / 2) after it. -size a,b first moves the peak by
	floor((a + b) / 2) along its strand and then does the same with
	N = b - a + 1, which gives b - a + 2 bases, one more than the usage's
	"from a to b" (the tests accept either; see
	test_windows_and_sites_for_each_size). The output's Start and End columns
	are the window's. The
	window's sequence is reverse-complemented for a - strand peak, and with
	-mask its lowercase bases are N.
Positions. homer2 find scans each window with -offset -floor(N / 2) (a for
	-size a,b, 0 for a peak used as given), so positions count from the centre
	in the peak's orientation. homer2 reports a - strand site at its last base;
	annotatePeaks.pl subtracts length - 1, so a reported position is the
	site's leftmost base in the window. Sites come from homer_helpers.scan(),
	which lists them in homer2's order: + strand from the left, then - strand
	from the right.
Motif columns. One per motif name: the sites sorted by position, each
	"pos(site,strand,0.00)" with the window's (+ strand) text of the site,
	comma-separated.
CpG% and GC% (homerTools freq -gc). GC is (C + G) / (A + C + G + T) over the
	window, case ignored and other characters not counted; CpG is the number
	of adjacent CG pairs over the number of adjacent pairs that are both
	A, C, G or T. Both printed with six decimals.
"""

import math
import os

import pytest

import homer_helpers as hh

MOTIF_COLUMN = ' Distance From Peak(sequence,strand,conservation)'
SCORE_COLUMN = ' Best Motif log-odds Score'


# --- the model ---------------------------------------------------------------------

def window(peak, size):
	"""(start, end, offset): the window's 1-based inclusive coordinates and
	the position its first base is reported at. `size` is 'given', a number
	or an (a, b) tuple."""

	s, e = peak['start'], peak['end']
	if size == 'given':
		return s, e, 0
	if isinstance(size, tuple):
		a, b = size
		n, move, offset = b - a + 1, (a + b) // 2, a
	else:
		n, move, offset = size, 0, -(size // 2)
	if peak['strand'] == '-':
		s, e = s - move, e - move
		end = (s + e + 1) // 2 + n // 2
		return end - n, end, offset
	s, e = s + move, e + move
	start = (s + e) // 2 - n // 2
	return start, start + n, offset


def window_sequence(genome, peak, start, end, mask=False):
	seq = genome[peak['chr']][start - 1:end]
	if mask:
		seq = ''.join('N' if base.islower() else base for base in seq)
	if peak['strand'] == '-':
		seq = hh.revcomp(seq)
	return seq


def peak_sites(genome, peak, motif, size, strand='both', mask=False, w=None):
	"""The sites annotatePeaks.pl records in one peak, in homer2's order, as
	(position, site, strand, score); `w` replaces window(peak, size)."""

	start, end, offset = w or window(peak, size)
	seq = window_sequence(genome, peak, start, end, mask)
	sites = []
	for pos, site, site_strand, score in hh.scan(seq, motif, strand, offset):
		if site_strand == '-':
			pos -= len(motif) - 1
		sites.append((pos, site, site_strand, score))
	return sites


def column_text(sites):
	return ','.join('{}({},{},0.00)'.format(pos, site, strand)
		for pos, site, strand, _ in sorted(sites, key=lambda s: s[0]))


def parse_column(text):
	"""[(position, site, strand)] from a motif column, in column order."""

	if text == '':
		return []
	sites = []
	for item in text.split('),'):
		pos, rest = item.rstrip(')').split('(')
		site, strand, cons = rest.split(',')
		assert cons == '0.00', text
		sites.append((int(pos), site, strand))
	return sites


def assert_column(text, sites):
	"""The column holds exactly `sites`, sorted by position (the order of
	equal positions is homer2's)."""

	got = parse_column(text)
	assert [g[0] for g in got] == sorted(g[0] for g in got), text
	assert sorted(got) == sorted(s[:3] for s in sites), text


def cpg_gc(seq):
	seq = seq.upper()
	bases = sum(1 for b in seq if b in 'ACGT')
	gc = sum(1 for b in seq if b in 'CG')
	pairs = cpg = 0
	for a, b in zip(seq, seq[1:]):
		if a in 'ACGT' and b in 'ACGT':
			pairs += 1
			cpg += a == 'C' and b == 'G'
	return '{:.6f}'.format(cpg / max(pairs, 0.01)), '{:.6f}'.format(gc / bases)


def remove_close(sites, distance, opposite, keep_first):
	"""annotatePeaks.pl's -rm / -rmrevopp, step for step: after each site
	is added, every still-valid site i invalidates each other site j within
	`distance` bp (positions differing by less than it) on the same strand
	(or, with `opposite`, the other strand), and itself too unless
	`keep_first`."""

	valid = []
	for k in range(len(sites)):
		valid.append(True)
		for i in range(k + 1):
			if not valid[i]:
				continue
			p1, d1 = sites[i][0], sites[i][2]
			for j in range(k + 1):
				if j == i:
					continue
				p2, d2 = sites[j][0], sites[j][2]
				if abs(p1 - p2) < distance and (d1 != d2) == opposite:
					if not keep_first:
						valid[i] = False
					valid[j] = False
	return [site for site, ok in zip(sites, valid) if ok]


def motif(name):
	return hh.read_motifs(os.path.join(hh.MOTIFS, name + '.motif'))[0]


def peak_scores():
	return {row[0]: row[5] for row in hh.read_tsv(hh.PEAKS) if not row[0].startswith('#')}


# --- the runs ----------------------------------------------------------------------

def custom_motifs(directory):
	"""Motif files the tests write: the 3-mer ACG (many sites, often close
	together) and REST-NRSF from HOMER's library."""

	library = hh.read_motifs(os.path.join(hh.HOMER_HOME, 'data', 'knownTFs',
		'vertebrates', 'known.motifs'))
	rest = [m for m in library if m.name.startswith('REST-NRSF')]
	assert len(rest) == 1
	paths = {
		'acg': hh.write_motifs(str(directory / 'acg.motif'),
			[hh.motif_from_consensus('ACG', 'acg', threshold=3.5)]),
		'rest': hh.write_motifs(str(directory / 'rest.motif'), rest),
	}
	for name in ('siteA', 'siteB', 'ebox', 'siteA_siteB', 'known_small'):
		ext = '.motif' if os.path.exists(os.path.join(hh.MOTIFS, name + '.motif')) else '.motifs'
		paths[name] = os.path.join(hh.MOTIFS, name + ext)
	return paths


def calls(m):
	"""key -> (peak file, options)."""

	a, ab = m['siteA'], m['siteA_siteB']
	peaks = hh.PEAKS
	return {
		'default': (peaks, ['-m', ab]),
		'given': (peaks, ['-m', ab, '-size', 'given']),
		'size200': (peaks, ['-m', ab, '-size', 200]),
		'size1000': (peaks, ['-m', ab, '-size', 1000]),
		'range': (peaks, ['-m', ab, '-size', '-100,50']),
		'files': (peaks, ['-m', m['siteA'], m['siteB'], m['ebox'], m['rest'], '-size', 200]),
		'nmotifs': (peaks, ['-m', m['ebox'], m['siteA'], '-size', 200, '-nmotifs']),
		'mdist': (peaks, ['-m', m['ebox'], m['siteA'], '-size', 200, '-mdist']),
		'mscore': (peaks, ['-m', ab, '-size', 200, '-mscore']),
		'norevopp': (peaks, ['-m', ab, '-size', 200, '-norevopp']),
		'norevopp_mscore': (peaks, ['-m', ab, '-size', 200, '-norevopp', '-mscore']),
		'strand+': (peaks, ['-m', ab, '-size', 200, '-strand', '+']),
		'strand-': (peaks, ['-m', ab, '-size', 200, '-strand', '-']),
		'mbed': (peaks, ['-m', ab, '-size', 200, '-mbed', 'sites.bed', '-mfasta', 'sites.fa']),
		'mpeak': (peaks, ['-m', ab, '-size', 200, '-mpeak', 'sites.txt']),
		'fm': (peaks, ['-m', m['acg'], '-fm', m['ebox'], '-size', 200]),
		'rm': (peaks, ['-m', m['acg'], '-size', 200, '-rm', 10]),
		'rm_opposite': (peaks, ['-m', m['acg'], '-size', 200, '-rm', -10]),
		'rmrevopp': (peaks, ['-m', m['acg'], '-size', 200, '-rmrevopp', 10]),
		'rmrevopp_nmotifs': (peaks, ['-m', m['ebox'], '-size', 200, '-rmrevopp', 1, '-nmotifs']),
		'matrix': (peaks, ['-m', m['known_small'], '-size', 200, '-matrix', 'mx']),
		'matrix_min': (peaks, ['-m', ab, '-size', 200, '-matrix', 'mx', '-matrixMinDist', 46]),
		'matrix_max': (peaks, ['-m', ab, '-size', 200, '-matrix', 'mx', '-matrixMaxDist', 45]),
		'matrix_min0': (peaks, ['-m', ab, '-size', 200, '-matrix', 'mx', '-matrixMinDist', 0]),
		'mlogic': (peaks, ['-m', ab, '-size', 200, '-mlogic', 'logic.txt']),
		'cpu2': (peaks, ['-m', ab, '-size', 200, '-cpu', 2]),
		'mask': (peaks, ['-m', a, '-size', 200, '-mask']),
		'center': (peaks, ['-center', a, '-size', 200]),
		'center_mirror': (peaks, ['-center', a, '-size', 200, '-mirror']),
		'center_multi': (peaks, ['-center', a, '-size', 200, '-multi']),
		'center_norevopp': (peaks, ['-center', a, '-size', 200, '-norevopp']),
		'hist': (peaks, ['-m', a, '-size', 200, '-hist', 10]),
		'hist_di': (peaks, ['-m', a, '-size', 200, '-hist', 10, '-di']),
		'hist_given': (peaks, ['-m', a, '-size', 'given', '-hist', 10]),
		'hist_nuc': (peaks, ['-size', 200, '-hist', 20, '-nuc']),
		'ghist': (peaks, ['-m', ab, '-size', 200, '-hist', 10, '-ghist']),
		'cpg': (peaks, ['-size', 200, '-CpG']),
		'nuc': (peaks, ['-size', 200, '-nuc']),
		'di': (peaks, ['-size', 200, '-di']),
		'bed': (hh.PEAKS_BED, ['-m', ab, '-size', 200]),
		'homer1': (peaks, ['-m', ab, '-size', 200, '-homer1']),
		'homer2': (peaks, ['-m', ab, '-size', 200, '-homer2']),
	}


class Run(object):

	def __init__(self, result, cwd):
		self.result = result
		self.cwd = cwd

	def ok(self):
		assert self.result.returncode == 0, self.result.stderr[-2000:]
		return self

	def table(self):
		"""(header, {peak ID: columns}) of the annotation table."""

		self.ok()
		lines = self.result.stdout.rstrip('\n').split('\n')
		header = lines[0].split('\t')
		rows = {}
		for line in lines[1:]:
			cols = line.split('\t')
			assert len(cols) == len(header), line
			assert cols[0] not in rows
			rows[cols[0]] = cols
		return header, rows

	def lines(self):
		self.ok()
		return [line.split('\t') for line in self.result.stdout.splitlines()]

	def path(self, name):
		self.ok()
		return os.path.join(str(self.cwd), name)


@pytest.fixture(scope='module')
def runs(tmp_path_factory):
	root = tmp_path_factory.mktemp('annotatePeaks')
	paths = custom_motifs(root)
	todo = {}
	for key, (peak_file, options) in calls(paths).items():
		args = ['annotatePeaks.pl', peak_file, hh.GENOME, '-noann', '-nogene'] + options
		todo[key] = (args, root / key.replace('+', 'plus').replace('-', 'minus'))
	results = hh.run_many(todo, check=False)
	runs = {key: Run(results[key], todo[key][1]) for key in todo}
	runs['paths'] = paths
	return runs


def column(header, name, kind=MOTIF_COLUMN):
	return header.index(name + kind)


# --- motif columns ---------------------------------------------------------------------

SIZES = {'default': 'given', 'given': 'given', 'size200': 200, 'size1000': 1000,
	'range': (-100, 50)}
LENGTHS = {'default': 200, 'given': 200, 'size200': 201, 'size1000': 1001}


@pytest.mark.parametrize('key', sorted(SIZES))
def test_windows_and_sites_for_each_size(runs, genome, peaks, key):
	"""-size N, -size a,b and -size given (the default without -hist): the
	Start and End columns are the window (N + 1 bases, or the 200 bp peak)
	and each motif column lists the window's sites, positions relative to
	its centre (or first base). The usage says -size a,b covers a to b from
	the centre, b - a + 1 bases, but annotatePeaks.pl gives b - a + 2 (to
	b + 1); either is accepted, with the window's sites, positions counted
	from a in both."""

	header, rows = runs[key].table()
	size = SIZES[key]
	motifs = hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA_siteB.motifs'))
	assert sorted(rows) == sorted(p['id'] for p in peaks)
	lengths = set()
	for peak in peaks:
		row = rows[peak['id']]
		start, end, offset = window(peak, size)
		windows = [(start, end, offset)]
		if key == 'range':
			windows.append((start, end - 1, offset) if peak['strand'] == '+'
				else (start + 1, end, offset))
		matching = [w for w in windows if row[1:5] == [peak['chr'], str(w[0]), str(w[1]),
			peak['strand']]]
		assert len(matching) == 1, (peak['id'], row[1:5], windows)
		w = matching[0]
		lengths.add(w[1] - w[0] + 1)
		for m in motifs:
			want = peak_sites(genome, peak, m, size, w=w)
			assert row[column(header, m.name)] == column_text(want), (peak['id'], m.name)
	assert len(lengths) == 1
	if key == 'range':
		assert lengths <= {151, 152}
	else:
		assert lengths == {LENGTHS[key]}


def test_explicit_given_is_the_default(runs):
	_, default = runs['default'].table()
	_, given = runs['given'].table()
	assert {k: v[1:] for k, v in default.items()} == {k: v[1:] for k, v in given.items()}


def test_planted_sites_are_at_the_centre(runs):
	"""siteA was planted 94 bases into each 200 bp peak, so with -size 200 a
	planted copy sits at -5 in either orientation, and every one is found."""

	header, rows = runs['size200'].table()
	col = column(header, 'siteA')
	planted = [p for p in hh.read_planted() if p['motif'] == 'siteA']
	assert len(planted) == 120
	for p in planted:
		positions = [site[0] for site in parse_column(rows[p['peak']][col])]
		assert -5 in positions, p


def test_cpg_and_gc_columns(runs, genome, peaks):
	header, rows = runs['size200'].table()
	assert header[7:9] == ['CpG%', 'GC%']
	for peak in peaks:
		start, end, _ = window(peak, 200)
		want = cpg_gc(window_sequence(genome, peak, start, end))
		assert tuple(rows[peak['id']][7:9]) == want, peak['id']


@pytest.mark.parametrize('flag', ['cpg', 'nuc', 'di'])
def test_sequence_options_add_cpg_and_gc(runs, genome, peaks, flag):
	"""-CpG, and -nuc and -di outside histogram mode, extract the sequence and
	add the CpG% and GC% columns, and nothing else."""

	header, rows = runs[flag].table()
	assert header[7:] == ['CpG%', 'GC%']
	for peak in peaks:
		start, end, _ = window(peak, 200)
		assert tuple(rows[peak['id']][7:9]) == cpg_gc(window_sequence(genome, peak, start,
			end)), peak['id']


def test_several_motif_files(runs, genome, peaks):
	"""-m with four files, one of them a motif from HOMER's library: a column
	per motif, and the 60 planted REST sites all found."""

	header, rows = runs['files'].table()
	paths = runs['paths']
	motifs = [hh.read_motifs(paths[name])[0] for name in ('siteA', 'siteB', 'ebox', 'rest')]
	assert len(header) == 9 + 4
	for peak in peaks:
		for m in motifs:
			assert_column(rows[peak['id']][column(header, m.name)],
				peak_sites(genome, peak, m, 200))
	rest = column(header, motifs[3].name)
	planted = {p['peak'] for p in hh.read_planted() if p['motif'] == 'REST'}
	assert len(planted) == 60
	assert all(rows[p][rest] != '' for p in planted)


def test_bed_peak_file_gives_the_same_table(runs):
	"""The same peaks as a BED file (0-based starts) give the same windows,
	sites and CpG/GC values."""

	header, txt = runs['size200'].table()
	_, bed = runs['bed'].table()
	assert sorted(bed) == sorted(txt)
	for peak_id, row in txt.items():
		assert bed[peak_id][1:5] == row[1:5]
		assert bed[peak_id][7:] == row[7:]


def test_threads_give_the_same_table(runs):
	_, one = runs['size200'].table()
	_, two = runs['cpu2'].table()
	assert {k: v[1:] for k, v in two.items()} == {k: v[1:] for k, v in one.items()}


@pytest.mark.parametrize('key', ['homer1', 'homer2'])
def test_homer1_and_homer2_find_the_same_sites(runs, genome, peaks, key):
	"""-homer1 scans with the original homer program, -homer2 (the default)
	with homer2: both give the reference sites."""

	header, rows = runs[key].table()
	for m in hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA_siteB.motifs')):
		col = column(header, m.name)
		for peak in peaks:
			assert_column(rows[peak['id']][col], peak_sites(genome, peak, m, 200))


def test_nmotifs_counts_sites(runs, genome, peaks):
	header, rows = runs['nmotifs'].table()
	for name in ('ebox', 'siteA'):
		col = column(header, name)
		for peak in peaks:
			want = len(peak_sites(genome, peak, motif(name), 200))
			assert rows[peak['id']][col] == str(want), (peak['id'], name)


def test_mdist_reports_the_closest_site(runs, genome, peaks):
	"""-mdist: the position nearest 0, the leftmost of two equally near;
	empty without a site."""

	header, rows = runs['mdist'].table()
	for name in ('ebox', 'siteA'):
		col = column(header, name)
		for peak in peaks:
			positions = sorted(s[0] for s in peak_sites(genome, peak, motif(name), 200))
			want = ''
			if positions:
				want = str(min(positions, key=lambda p: (abs(p), p)))
			assert rows[peak['id']][col] == want, (peak['id'], name)


@pytest.mark.parametrize('key,strand', [('mscore', 'both'), ('norevopp_mscore', '+')])
def test_mscore_reports_the_best_score(runs, genome, peaks, key, strand):
	"""-mscore: the best log-odds score of any window in the peak, threshold
	or not; with -norevopp only + strand windows count."""

	header, rows = runs[key].table()
	for m in hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA_siteB.motifs')):
		col = column(header, m.name, SCORE_COLUMN)
		for peak in peaks:
			start, end, _ = window(peak, 200)
			best = hh.best_score(window_sequence(genome, peak, start, end), m, strand)
			assert float(rows[peak['id']][col]) == pytest.approx(best, abs=1e-5), peak['id']


def test_norevopp_scans_only_the_peak_strand(runs, genome, peaks):
	header, rows = runs['norevopp'].table()
	for m in hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA_siteB.motifs')):
		col = column(header, m.name)
		for peak in peaks:
			want = peak_sites(genome, peak, m, 200, strand='+')
			assert rows[peak['id']][col] == column_text(want), peak['id']


@pytest.mark.parametrize('strand', ['+', '-'])
def test_strand_keeps_sites_on_one_strand(runs, genome, peaks, strand):
	"""-strand + or - keeps only the sites on that strand relative to the
	peak."""

	header, rows = runs['strand' + strand].table()
	for m in hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA_siteB.motifs')):
		col = column(header, m.name)
		for peak in peaks:
			want = [s for s in peak_sites(genome, peak, m, 200) if s[2] == strand]
			assert rows[peak['id']][col] == column_text(want), peak['id']


def test_mask_hides_lowercase_sequence(runs, genome, peaks):
	"""-mask extracts lowercase bases as N: the ten all-lowercase peaks lose
	their planted siteA, and CpG/GC count only the uppercase bases."""

	header, rows = runs['mask'].table()
	col = column(header, 'siteA')
	m = motif('siteA')
	for peak in peaks:
		row = rows[peak['id']]
		start, end, _ = window(peak, 200)
		assert row[col] == column_text(peak_sites(genome, peak, m, 200, mask=True)), peak['id']
		seq = window_sequence(genome, peak, start, end, mask=True)
		assert tuple(row[7:9]) == cpg_gc(seq), peak['id']
	lowercase = [p['peak'] for p in hh.read_planted() if p['motif'] == 'siteA'
		and p['lowercase'] == 'yes']
	assert len(lowercase) == 10
	assert all(rows[p][col] == '' for p in lowercase)


# --- filtering and removing sites -------------------------------------------------------

def test_fm_removes_sites_overlapping_a_filter_motif(runs, genome, peaks):
	"""-fm motifs are scanned first and their sites mask the -m sites: an
	ACG site that shares a base with an E-box site (CACGTG holds ACG) is
	dropped, one with at least a base between it and every E-box site is
	kept (sites directly next to one are not checked). The filter motif's
	own column stays empty."""

	header, rows = runs['fm'].table()
	acg = hh.read_motifs(runs['paths']['acg'])[0]
	ebox = motif('ebox')
	col = column(header, 'acg')
	dropped = kept = 0
	for peak in peaks:
		row = rows[peak['id']]
		assert row[column(header, 'ebox')] == ''
		filters = [(q, q + len(ebox) - 1) for q, _, _, _ in peak_sites(genome, peak, ebox, 200)]
		got = {site[:3] for site in parse_column(row[col])}
		for site in peak_sites(genome, peak, acg, 200):
			first, last = site[0], site[0] + len(acg) - 1
			if any(first <= f_last and f_first <= last for f_first, f_last in filters):
				assert site[:3] not in got, (peak['id'], site)
				dropped += 1
			elif all(first > f_last + 1 or last < f_first - 1 for f_first, f_last in filters):
				assert site[:3] in got, (peak['id'], site)
				kept += 1
		assert got <= {site[:3] for site in peak_sites(genome, peak, acg, 200)}
	assert dropped > 0 and kept > 0


@pytest.mark.parametrize('key,distance,opposite,keep_first', [
	('rm', 10, False, False),
	('rm_opposite', 10, True, False),
	('rmrevopp', 10, True, True),
])
def test_close_sites_are_removed(runs, genome, peaks, key, distance, opposite, keep_first):
	"""-rm <#> drops every site that has another site of the same motif and
	strand less than # bp away (a negative # means the other strand);
	-rmrevopp <#> keeps the site found first and drops later other-strand
	sites within # bp of it. Modelled step by step in remove_close()."""

	header, rows = runs[key].table()
	m = hh.read_motifs(runs['paths']['acg'])[0]
	col = column(header, 'acg')
	removed = 0
	for peak in peaks:
		sites = peak_sites(genome, peak, m, 200)
		want = remove_close(sites, distance, opposite, keep_first)
		removed += len(sites) - len(want)
		assert_column(rows[peak['id']][col], want)
	assert removed > 0


def test_rmrevopp_counts_a_palindrome_once(runs, genome, peaks):
	"""CACGTG is found on both strands at the same position; -rmrevopp 1
	keeps one of each pair, so -nmotifs counts each site once."""

	header, rows = runs['rmrevopp_nmotifs'].table()
	col = column(header, 'ebox')
	total = 0
	for peak in peaks:
		positions = {s[0] for s in peak_sites(genome, peak, motif('ebox'), 200)}
		total += len(positions)
		assert rows[peak['id']][col] == str(len(positions)), peak['id']
	assert total > 0


# --- site files ---------------------------------------------------------------------------

def genome_interval(peak, w, m, site):
	"""The site's 0-based, end-exclusive genome interval and the strand the
	motif reads on: a + strand peak's site starts (position - offset) bases
	into the window, a - strand peak's ends that many bases before the
	window's end, and the strand flips for a - strand peak."""

	start, end, offset = w
	index = site[0] - offset
	if peak['strand'] == '+':
		bed_start = start - 1 + index
		strand = site[2]
	else:
		bed_start = end - index - len(m)
		strand = '-' if site[2] == '+' else '+'
	return bed_start, bed_start + len(m), strand


def site_records(genome, peaks, size=200):
	"""(peak, motif, site, genome interval) for every siteA and siteB site."""

	records = []
	for m in hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA_siteB.motifs')):
		for peak in peaks:
			w = window(peak, size)
			for site in peak_sites(genome, peak, m, size):
				records.append((peak, m, site, genome_interval(peak, w, m, site)))
	return records


def assert_scores(got, want):
	assert sorted(got) == sorted(want)
	for key, score in got.items():
		assert score == pytest.approx(want[key], abs=1e-5), key


def test_mbed_writes_genome_coordinates(runs, genome, peaks):
	"""-mbed: a track line, then chr, start, end, motif, score and strand for
	every site; each interval holds the site, read on the reported strand."""

	with open(runs['mbed'].path('sites.bed')) as handle:
		lines = handle.read().splitlines()
	assert lines[0].startswith('track name=')
	rows = [line.split('\t') for line in lines[1:]]
	got = {(r[0], int(r[1]), int(r[2]), r[3], r[5]): float(r[4]) for r in rows}
	assert len(got) == len(rows)
	assert_scores(got, {(peak['chr'], start, end, m.name, strand): site[3]
		for peak, m, site, (start, end, strand) in site_records(genome, peaks)})
	motifs = {m.name: m for m in hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA_siteB.motifs'))}
	for chrom, start, end, name, strand in got:
		text = genome[chrom][start:end].upper()
		if strand == '-':
			text = hh.revcomp(text)
		assert hh.scan(text, motifs[name], '+'), (chrom, start, end, name)


def test_mpeak_writes_a_peak_file(runs, genome, peaks):
	"""-mpeak: motif, chr, start (1-based), end, strand, score and the site
	text, without a track line."""

	rows = hh.read_tsv(runs['mpeak'].path('sites.txt'))
	got = {(r[0], r[1], int(r[2]), int(r[3]), r[4], r[6]): float(r[5]) for r in rows}
	assert len(got) == len(rows)
	assert_scores(got, {(m.name, peak['chr'], start + 1, end, strand, site[1]): site[3]
		for peak, m, site, (start, end, strand) in site_records(genome, peaks)})
	for r in rows:
		assert genome[r[1]][int(r[2]) - 1:int(r[3])].upper() in (r[6], hh.revcomp(r[6]))


def test_mfasta_writes_each_site_in_motif_orientation(runs, genome, peaks):
	"""-mfasta: one record per site named "motif score peak position strand",
	the site's sequence read in the motif's orientation."""

	records = hh.read_fasta(runs['mbed'].path('sites.fa'))
	got = {}
	for name, seq in records:
		motif_name, score, peak_id, pos, strand = name.split(' ')
		got[(motif_name, peak_id, int(pos), strand, seq)] = float(score)
	assert len(got) == len(records)
	assert_scores(got, {(m.name, peak['id'], site[0], site[2],
		site[1] if site[2] == '+' else hh.revcomp(site[1])): site[3]
		for peak, m, site, _ in site_records(genome, peaks)})


# --- co-occurrence and orientation ----------------------------------------------------------

def midpoints(genome, peaks, motifs, size=200):
	"""{motif: {peak: [site midpoints]}} for peaks with a site; a midpoint is
	position + floor(length / 2)."""

	out = {}
	for m in motifs:
		out[m.name] = {}
		for peak in peaks:
			sites = peak_sites(genome, peak, m, size)
			if sites:
				out[m.name][peak['id']] = [s[0] + len(m) // 2 for s in sites]
	return out


def read_matrix(path):
	rows = hh.read_tsv(path)
	names = [c.rsplit(' (', 1)[0] for c in rows[0][1:]]
	values = {}
	for row in rows[1:]:
		for name, value in zip(names, row[1:]):
			values[(row[0].rsplit(' (', 1)[0], name)] = float(value)
	return names, values


@pytest.mark.parametrize('key,file,min_dist,max_dist', [
	('matrix', 'known_small.motifs', 4, 1e10),
	('matrix_min', 'siteA_siteB.motifs', 46, 1e10),
	('matrix_max', 'siteA_siteB.motifs', 4, 45),
	('matrix_min0', 'siteA_siteB.motifs', 0, 1e10),
])
def test_matrix_counts_co_occurring_peaks(runs, genome, peaks, key, file, min_dist, max_dist):
	"""-matrix: the count matrix holds, for each pair of motifs, the peaks
	with a site of each whose midpoints are -matrixMinDist (default 4, and at
	least 4 for a motif with itself) to -matrixMaxDist apart; a midpoint is
	position + floor(length / 2), so the planted siteA-siteB pairs are 45 bp
	apart, inside -matrixMaxDist 45 and outside -matrixMinDist 46. Off the diagonal the ratio is that
	count over the expected count n1 * n2 / peaks (at least 0.5), and the log
	p-value is ln P(X >= count) for X ~ Binomial(peaks, expected / peaks), or
	-ln P(X <= count) when the ratio is below 1."""

	motifs = hh.read_motifs(os.path.join(hh.MOTIFS, file))
	mids = midpoints(genome, peaks, motifs)
	total = len(peaks)
	names, counts = read_matrix(runs[key].path('mx.count.matrix.txt'))
	_, ratios = read_matrix(runs[key].path('mx.ratio.matrix.txt'))
	_, logps = read_matrix(runs[key].path('mx.logPvalue.matrix.txt'))
	assert names == [m.name for m in motifs]
	stats = {r[0]: r[1:] for r in hh.read_tsv(runs[key].path('mx.stats.txt'))[1:]}
	for a in names:
		for b in names:
			low = max(min_dist, 4) if a == b else min_dist
			n = sum(1 for peak in mids[a] if peak in mids[b] and any(
				low <= abs(x - y) <= max_dist for x in mids[a][peak] for y in mids[b][peak]))
			assert counts[(a, b)] == n, (a, b)
			if a == b:
				continue
			n1, n2 = len(mids[a]), len(mids[b])
			expected = n1 * n2 / total
			ratio = n / max(expected, 0.5)
			r = min(max(max(expected, 0.5) / total, 1 / total), (total - 1) / total)
			if ratio >= 1:
				logp = hh.log_binomial_tail(n, total, r)
			else:
				logp = -hh.log_binomial_tail(total - n, total, 1 - r)
			assert ratios[(a, b)] == pytest.approx(ratio, rel=1e-9), (a, b)
			assert logps[(a, b)] == pytest.approx(logp, rel=1e-4, abs=1e-6), (a, b)
			pair = stats.get(a + ':' + b)
			if pair is not None:
				assert pair[0] == str(total)
				assert sorted(int(v) for v in pair[1:3]) == sorted([n1, n2])
				assert int(pair[3]) == n
				assert float(pair[4]) == pytest.approx(expected, rel=1e-9)
	pair = counts[('siteA', 'siteB')]
	if key == 'matrix_min':
		assert pair < 40
	else:
		assert pair >= 40


def revopp_code(code):
	return code[::-1].swapcase()


def test_mlogic_reports_motif_orders(runs, genome, peaks):
	"""-mlogic: each motif gets a letter (A, B, ... on the + strand,
	lowercase on the -), each peak the letters of its sites in position
	order, written in whichever of its two orientations sorts first; the
	table counts the peaks per code exactly and the peaks whose code holds
	it (either way round) as a substring."""

	motifs = hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA_siteB.motifs'))
	codes = {}
	for peak in peaks:
		hits = []
		for k, m in enumerate(motifs):
			for pos, _, strand, _ in peak_sites(genome, peak, m, 200):
				hits.append((pos, chr(65 + k) if strand == '+' else chr(97 + k)))
		assert len({h[0] for h in hits}) == len(hits), 'positions should not tie'
		code = ''.join(c for _, c in sorted(hits))
		codes[peak['id']] = min(code, revopp_code(code))
	with open(runs['mlogic'].path('logic.txt')) as handle:
		sections = handle.read().split('\n\n\n')
	assert sections[0].split('\n') == ['Motif Name\t+ strand code\t- strand code',
		'siteA\tA\ta', 'siteB\tB\tb']
	table = [line.split('\t') for line in sections[1].split('\n')]
	assert table[0] == ['CRM codes', '# of Peaks', '% of Peaks', '# of Peaks (exact)',
		'% of Peaks (exact)']
	got = {row[0]: row[1:] for row in table[1:]}
	want = {}
	for code in set(codes.values()):
		within = sum(1 for c in codes.values() if code in c or revopp_code(code) in c)
		exact = sum(1 for c in codes.values() if c == code)
		want[code or 'None'] = [str(within), '{:.3f}'.format(within / 200), str(exact),
			'{:.3f}'.format(exact / 200)]
	assert got == want
	assert [int(row[1]) for row in table[1:]] == sorted((int(row[1]) for row in table[1:]),
		reverse=True)
	per_peak = [line.split('\t') for line in sections[2].strip('\n').split('\n')]
	assert per_peak[0] == ['PeakID', 'chr', 'start', 'end', 'strand', 'CRM code']
	by_id = {p['id']: p for p in peaks}
	for row in per_peak[1:]:
		peak = by_id[row[0]]
		start, end, _ = window(peak, 200)
		assert row[1:] == [peak['chr'], str(start), str(end), peak['strand'], codes[row[0]]]
	assert len(per_peak) == 201


# --- centring ----------------------------------------------------------------------

def centre_lines(genome, peaks, mirror=False, multi=False, strand='both'):
	"""What -center prints: homer2 find positions as reported (a - strand
	site at its last base), negated with -mirror; each peak's site nearest 0
	(the first found of equally near ones) moves the window by that many bases
	along the peak's strand and turns the peak to the - strand when the site
	is on the - strand relative to it; peaks without a site within the size
	are dropped. -multi prints a line for every site instead."""

	m = motif('siteA')
	scores = peak_scores()
	lines = []
	for peak in peaks:
		start, end, offset = window(peak, 200)
		seq = window_sequence(genome, peak, start, end)
		best = None
		for pos, _, site_strand, score in hh.scan(seq, m, strand, offset):
			pos = -pos if mirror else pos
			if peak['strand'] == '+':
				line = (start + pos, end + pos, 1 if site_strand == '-' else 0)
			else:
				line = (start - pos, end - pos, 0 if site_strand == '-' else 1)
			line = (peak['id'], peak['chr']) + line + (scores[peak['id']], score)
			if multi:
				lines.append(line)
			elif best is None or abs(pos) < abs(best[0]):
				best = (pos, line)
		if best is not None and abs(best[0]) <= 200:
			lines.append(best[1])
	return sorted(lines)


def parse_centre(run):
	return sorted((r[0], r[1], int(r[2]), int(r[3]), int(r[4]), r[5], float(r[6]))
		for r in run.lines())


def assert_centre(got, want):
	assert [g[:6] for g in got] == [w[:6] for w in want]
	for g, w in zip(got, want):
		assert g[6] == pytest.approx(w[6], abs=1e-5)


@pytest.mark.parametrize('key,options', [
	('center', {}),
	('center_mirror', {'mirror': True}),
	('center_multi', {'multi': True}),
	('center_norevopp', {'strand': '+'}),
])
def test_center_moves_peaks_onto_the_motif(runs, genome, peaks, key, options):
	"""-center <motif> prints only the recentred peaks (ID, chr, start, end,
	strand 0/1, peak score, motif score); -mirror flips positions, -multi
	prints every site and -norevopp scans the peak's strand only."""

	got = parse_centre(runs[key])
	want = centre_lines(genome, peaks, **options)
	assert_centre(got, want)
	if key == 'center':
		assert len(got) == 120


def test_centred_windows_start_on_the_site(runs, genome):
	"""After -center, a 201 bp window around the new centre has the site's
	first base in its middle: siteA in the motif's orientation starts at
	the new centre."""

	m = motif('siteA')
	for peak_id, chrom, start, end, strand, _, _ in parse_centre(runs['center']):
		assert end - start == 200
		centre = (start + end) // 2 if strand == 0 else (start + end + 1) // 2
		if strand == 0:
			text = genome[chrom][centre - 1:centre - 1 + len(m)]
		else:
			text = hh.revcomp(genome[chrom][centre - len(m):centre])
		assert hh.scan(text, m, '+'), peak_id


# --- histograms ----------------------------------------------------------------------------

def frequency_table(genome, peaks, size, offset):
	"""homerTools freq over the windows, as -hist reads it: per position
	(offset + index, the first `size` bases), each column's value as printed
	(six decimals)."""

	seqs = []
	for peak in peaks:
		start, end, _ = window(peak, size)
		seqs.append(window_sequence(genome, peak, start, end).upper())
	pairs = [a + b for a in 'ACGT' for b in 'ACGT']
	names = ['A', 'C', 'G', 'T', 'A/T', 'C/G', 'A/G', 'C/T'] + pairs
	groups = {'A/T': ('AA', 'AT', 'TA', 'TT'), 'C/G': ('CC', 'CG', 'GC', 'GG'),
		'A/G': ('AA', 'AG', 'GA', 'GG'), 'C/T': ('CC', 'CT', 'TC', 'TT')}
	table = {}
	for i in range(size):
		mono = [s[i] for s in seqs if s[i] in 'ACGT']
		di = [s[i:i + 2] for s in seqs if s[i] in 'ACGT' and s[i + 1] in 'ACGT']
		row = {b: float('{:.6f}'.format(mono.count(b) / len(mono))) for b in 'ACGT'}
		for pair in pairs:
			row[pair] = float('{:.6f}'.format(di.count(pair) / len(di)))
		for name, members in groups.items():
			row[name] = float('{:.6f}'.format(sum(di.count(p) for p in members) / len(di)))
		table[offset + i] = row
	return names, table


def read_histogram(run):
	rows = run.lines()
	return rows[0][1:], {int(r[0]): [float(v) for v in r[1:]] for r in rows[1:]}


def bin_of(value, width):
	return int(math.floor(value / width + 0.5)) * width


def motif_histogram(genome, peaks, m, width):
	"""Site counts per bin over all peaks, divided by bin width x peaks: a
	site's 5' end (first base on the + strand, last on the -) goes to the
	+ or - sites column, its middle (5' end plus or minus floor(length / 2))
	to the total column."""

	total, plus, minus = {}, {}, {}
	for peak in peaks:
		for pos, _, strand, _ in peak_sites(genome, peak, m, 200):
			half = len(m) // 2
			if strand == '+':
				plus[bin_of(pos, width)] = plus.get(bin_of(pos, width), 0) + 1
				middle = pos + half
			else:
				pos += len(m) - 1
				minus[bin_of(pos, width)] = minus.get(bin_of(pos, width), 0) + 1
				middle = pos - half
			total[bin_of(middle, width)] = total.get(bin_of(middle, width), 0) + 1
	scale = width * len(peaks)
	return [{k: v / scale for k, v in h.items()} for h in (total, plus, minus)]


def frequency_bins(table, name, width):
	sums, counts = {}, {}
	for pos in sorted(table):
		b = bin_of(pos, width)
		sums[b] = sums.get(b, 0.0) + table[pos][name]
		counts[b] = counts.get(b, 0) + 1
	return {b: sums[b] / counts[b] for b in sums}


@pytest.mark.parametrize('key,di', [('hist', False), ('hist_di', True)])
def test_hist_motif_and_nucleotide_columns(runs, genome, peaks, key, di):
	"""-hist 10 with -m: rows -100 to 100 by 10; per motif the total, + and -
	site densities (model in motif_histogram()), then the mean frequency of
	each base over the positions in a bin (from homerTools freq over the
	first 200 bases of the windows, positions from -100); -di adds the
	dinucleotide columns."""

	names, rows = read_histogram(runs[key])
	m = motif('siteA')
	freq_names, table = frequency_table(genome, peaks, 200, -100)
	freq_names = freq_names if di else freq_names[:4]
	assert names == ['siteA total sites', 'siteA + sites', 'siteA - sites'] + [
		name + ' frequency' for name in freq_names]
	assert sorted(rows) == list(range(-100, 101, 10))
	hists = motif_histogram(genome, peaks, m, 10)
	freqs = [frequency_bins(table, name, 10) for name in freq_names]
	for b, values in rows.items():
		for got, hist in zip(values[:3], hists):
			assert got == pytest.approx(hist.get(b, 0.0), abs=1e-12), b
		for got, freq in zip(values[3:], freqs):
			assert got == pytest.approx(freq[b], abs=1e-9), b
	assert rows[0][0] > 0.0


def test_hist_given_size_bins_by_tenths_of_the_peak(runs, genome, peaks):
	"""-hist 10 with -size given: ten rows labelled 0.1 to 1, one per tenth
	of the peak (a position p of a peak with end - start = W falls in tenth
	floor(10 p / W)); a site adds 10 / W to the tenth its 5' end or middle
	falls in, and the sums are divided by the number of peaks."""

	rows = runs['hist_given'].lines()
	assert rows[0][1:4] == ['siteA total sites', 'siteA + sites', 'siteA - sites']
	assert [float(r[0]) for r in rows[1:]] == pytest.approx([(i + 1) / 10.0 for i in range(10)])
	m = motif('siteA')
	hists = [{}, {}, {}]
	for peak in peaks:
		width = peak['end'] - peak['start']
		share = 10.0 / width
		for pos, _, strand, _ in peak_sites(genome, peak, m, 'given'):
			if strand == '+':
				five, middle, h = pos, pos + len(m) // 2, hists[1]
			else:
				five = pos + len(m) - 1
				middle, h = five - len(m) // 2, hists[2]
			for value, target in ((five, h), (middle, hists[0])):
				if 0 <= value <= width:
					b = int(math.floor(value / width * 10))
					target[b] = target.get(b, 0.0) + share / len(peaks)
	for i, row in enumerate(rows[1:]):
		# The frequency columns are skipped: in v5.1 they are 0 after the first row.
		for got, hist in zip(row[1:4], hists):
			assert float(got) == pytest.approx(hist.get(i, 0.0), abs=1e-12), i
	assert sum(float(r[1]) for r in rows[1:]) > 0


def test_hist_nuc_without_motifs(runs, genome, peaks):
	"""-hist 20 -nuc alone: only the four base-frequency columns, rows -100
	to 100 by 20."""

	names, rows = read_histogram(runs['hist_nuc'])
	assert names == ['A frequency', 'C frequency', 'G frequency', 'T frequency']
	assert sorted(rows) == list(range(-100, 101, 20))
	_, table = frequency_table(genome, peaks, 200, -100)
	for k, base in enumerate('ACGT'):
		freq = frequency_bins(table, base, 20)
		for b, values in rows.items():
			assert values[k] == pytest.approx(freq[b], abs=1e-9), (b, base)


def test_ghist_writes_one_profile_per_peak(runs, genome, peaks):
	"""-ghist: a row per peak with, for each motif, the number of site
	middles in each 10 bp bin from -100 to 100."""

	rows = runs['ghist'].lines()
	bins = list(range(-100, 101, 10))
	assert rows[0] == ['Gene'] + [str(b) for b in bins] * 2
	got = {r[0]: [int(v) for v in r[1:]] for r in rows[1:]}
	assert sorted(got) == sorted(p['id'] for p in peaks)
	motifs = hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA_siteB.motifs'))
	for peak in peaks:
		want = []
		for m in motifs:
			counts = motif_histogram(genome, [peak], m, 10)[0]
			want += [int(round(counts.get(b, 0.0) * 10)) for b in bins]
		assert got[peak['id']] == want, peak['id']

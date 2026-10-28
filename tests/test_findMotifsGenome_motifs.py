"""findMotifsGenome.pl options that change how motifs are scored, found and
reported, with a custom genome: -norevopp, -h, -flip, -float, -rand,
-maskMotif, -find, -mset, -mknown, -mcheck, -S, -mis, -oligo, -e, -olen,
-cache, -minlp, -quickMask, -noknown, -nomotif, -rna, -fdr, -basic, -bits,
-nofacts and -nocheck. Each is checked against a run without it or against
a reference computation.

All runs start in one batch (module fixture `runs`). Expected values come from
the inputs and from the run's own sequence and group files (-keepFiles):
sites from homer_helpers.scan, known and de novo motif counts from that scan
weighted by the group file, p-values from the binomial or hypergeometric
tail, empirical FDRs as Statistics.pm defines them, logo heights from the
motif's information content. Options passed through to homer2 denovo are
checked by their documented effect where one is visible, or against homer2
denovo called directly on the same sequence and group files.
"""

import math
import os

import pytest

import homer_helpers as hh
import findmotifsgenome_helpers as fh

KS = fh.KNOWN_SMALL
SITE_A = os.path.join(hh.MOTIFS, 'siteA.motif')
SITE_B = os.path.join(hh.MOTIFS, 'siteB.motif')
KNOWN = ['-size', 200, '-nomotif', '-mknown', KS, '-N', 1000]
KEEP = ['-dumpFasta', '-keepFiles']
DENOVO = ['-size', 200, '-len', 8, '-S', 2, '-N', 400, '-mknown', KS, '-mcheck', KS]
SEVEN = ['-size', 200, '-len', 8, '-S', 7, '-N', 400, '-noknown', '-nocheck']


@pytest.fixture(scope='module')
def runs(tmp_path_factory, peaks, preparsed):
	tmp = tmp_path_factory.mktemp('fmg_motifs')
	few = str(tmp / 'few.txt')
	with open(few, 'w') as handle:
		for p in peaks[:20]:
			handle.write('{id}\t{chr}\t{start}\t{end}\t{strand}\n'.format(**p))
	pp = {'preparsed': preparsed}
	calls = {
		'base': (KNOWN + KEEP, pp),
		'norevopp': (KNOWN + KEEP + ['-norevopp'], pp),
		'hypergeometric': (KNOWN + KEEP + ['-h'], pp),
		'float': (KNOWN + KEEP + ['-float'], pp),
		'flip': (KNOWN + KEEP + ['-flip'], pp),
		'rand': (KNOWN + KEEP + ['-rand'], pp),
		'mask_motif': (KNOWN + KEEP + ['-maskMotif', SITE_A], pp),
		'mask_motifs': (KNOWN + KEEP + ['-maskMotif', SITE_A, SITE_B], pp),
		'find': (['-size', 200, '-find', os.path.join(hh.MOTIFS, 'siteA_siteB.motifs')], pp),
		'find_given': (['-size', 'given', '-find', SITE_A, '-norevopp', '-preparse'],
			{'preparsed': str(tmp / 'pp199')}),
		'mset': (['-size', 200, '-len', 8, '-S', 2, '-N', 400, '-mset', 'worms'], pp),
		'mset_auto': (['-size', 200, '-nomotif', '-N', 40], dict(pp, peaks=few)),
		'dn': (DENOVO + ['-keepFiles', '-oligo'], pp),
		'dn_S1': (['-size', 200, '-len', 8, '-S', 1, '-N', 400, '-mknown', KS, '-mcheck', KS], pp),
		'dn_mis0': (DENOVO + ['-oligo', '-mis', 0, '-nocheck'], pp),
		'dn_e': (DENOVO + ['-e', 1], pp),
		'dn_olen': (DENOVO + ['-olen', 2], pp),
		'dn_cache': (DENOVO + ['-cache', 1], pp),
		'dn_minlp_low': (SEVEN + ['-minlp', -100000], pp),
		'dn_minlp_high': (SEVEN + ['-minlp', 1], pp),
		'dn_quickmask': (SEVEN + ['-minlp', -100000, '-quickMask'], pp),
		'dn_norevopp': (DENOVO + ['-norevopp', '-keepFiles'], pp),
		'dn_rna': (['-rna'] + DENOVO, pp),
		'dn_noknown': (DENOVO + ['-noknown'], pp),
		'dn_fdr': (DENOVO + ['-fdr', 2], pp),
		'dn_basic': (DENOVO + ['-basic'], pp),
		'dn_bits': (DENOVO + ['-bits'], pp),
		'dn_nofacts': (DENOVO + ['-nofacts'], pp),
		'dn_nocheck': (DENOVO + ['-nocheck'], pp),
	}
	return fh.run_calls(tmp, calls)


@pytest.fixture(scope='module')
def library():
	return hh.read_motifs(KS)


# --- strands, statistics ------------------------------------------------------------------

def test_norevopp_counts_forward_strand_sites_only(runs, library):
	"""-norevopp: known motif counts take + strand sites only, and homer2 norm
	counts oligos on the given strand without pooling reverse complements."""

	fh.check_known_run(runs, 'norevopp', library, strand='+')
	rows, base = fh.known_rows(runs, 'norevopp'), fh.known_rows(runs, 'base')
	assert rows['siteA']['targets'] < base['siteA']['targets']
	header, rows = hh.read_tsv(fh.out(runs, 'norevopp', 'seq.autonorm.tsv'), header=True)
	targets = [(s, 1.0) for s in fh.out_fasta(runs, 'norevopp', 'target.fa').values()]
	texts = fh.autonorm_target_texts(targets, (1, 2, 3), both_strands=False)
	for row in rows:
		assert row[1] in texts.get(row[0], {'0.00%'}), row


def test_hypergeometric_p_values(runs, library):
	"""-h: p-values from the hypergeometric tail instead of the binomial."""

	fh.check_known_run(runs, 'hypergeometric', library, hypergeometric=True)


def test_flip_swaps_targets_and_background(runs, library):
	"""-flip swaps the labels after normalization, so each motif's target
	and background counts trade places."""

	base, flip = fh.known_rows(runs, 'base'), fh.known_rows(runs, 'flip')
	for motif in library:
		assert flip[motif.name]['targets'] == base[motif.name]['background']
		assert flip[motif.name]['background'] == base[motif.name]['targets']


def test_float_adjusts_known_motif_thresholds(runs, library):
	"""-float: each enriched motif's file in knownResults/ carries the
	threshold homer2 chose, its counts are those of a scan at that threshold,
	and no p-value is worse than at the library threshold."""

	seqs = fh.read_seq_tsv(fh.out(runs, 'float', 'seq.tsv'))
	groups = fh.read_groups(fh.out(runs, 'float', 'group.adj'))
	rows = hh.read_known_results(fh.out(runs, 'float', 'knownResults.txt'))
	given = {m.name: m.threshold for m in library}
	changed = 0
	for rank, row in enumerate([r for r in rows if math.exp(r['logp']) <= 0.01], 1):
		motif = hh.read_motifs(fh.out(runs, 'float', 'knownResults',
			'known{}.motif'.format(rank)))[0]
		changed += abs(motif.threshold - given[motif.name]) > 1e-3
		hits, background_hits = fh.weighted_hits(seqs, groups, motif)
		assert row['targets'] == pytest.approx(hits, abs=0.051)
		assert row['background'] == pytest.approx(background_hits, abs=0.051)
	assert changed >= 1
	base = fh.known_rows(runs, 'base')
	for row in rows:
		assert row['logp'] <= base[row['name']]['logp'] * (1 - 1e-3) + 1e-6


def test_rand_shuffles_the_labels(runs, library, peaks):
	"""-rand shuffles target and background labels: as many targets, mostly
	other sequences, and siteA no longer concentrated in them (about 120 x 200
	/ 1,000 = 24 expected)."""

	groups = fh.read_groups(fh.out(runs, 'rand', 'group.adj'))
	targets = {name for name, (g, _) in groups.items() if g == 1}
	assert len(targets) == len(peaks)
	assert len(targets & {p['id'] for p in peaks}) < 100
	fh.check_known_run(runs, 'rand', library)
	assert fh.known_rows(runs, 'rand')['siteA']['targets'] < 60


def test_mask_motif_removes_its_sites_before_finding(runs, library):
	"""-maskMotif runs homer2 mask on the sequences: every site becomes N, so
	the masked motifs match nothing."""

	base = fh.read_seq_tsv(fh.out(runs, 'base', 'seq.tsv'))
	for key, paths in (('mask_motif', [SITE_A]), ('mask_motifs', [SITE_A, SITE_B])):
		motifs = [hh.read_motifs(path)[0] for path in paths]
		want = {}
		for name, seq in base.items():
			for motif in motifs:
				seq = hh.mask_sequence(seq, motif)
			want[name] = seq
		assert fh.read_seq_tsv(fh.out(runs, key, 'seq.tsv')) == want
		rows = fh.known_rows(runs, key)
		for motif in motifs:
			assert (rows[motif.name]['targets'], rows[motif.name]['background']) == (0, 0)
		fh.check_known_run(runs, key, library)


# --- -find, motif sets ----------------------------------------------------------------------

def find_rows(text):
	lines = text.splitlines()
	assert lines[0] == 'PositionID\tOffset\tSequence\tMotif Name\tStrand\tMotifScore'
	return sorted((c[0], int(c[1]), c[2], c[3], c[4], round(float(c[5]), 4))
		for c in (line.split('\t') for line in lines[1:]))


@pytest.mark.parametrize('key', ['find', 'find_given'])
def test_find_reports_every_site(runs, peaks, genome, key):
	"""-find prints homer2 find's sites in each window, positions relative to
	-floor(size / 2) (0 for -size given), and writes no enrichment files."""

	if key == 'find':
		path, strand, offset = os.path.join(hh.MOTIFS, 'siteA_siteB.motifs'), 'both', -100
		seqs = {p['id']: fh.window_sequence(p, 200, genome=genome) for p in peaks}
	else:
		path, strand, offset = SITE_A, '+', 0
		seqs = {p['id']: hh.peak_sequence(p, genome) for p in peaks}
	want = sorted((name, pos, site, m.name, s, round(score, 4)) for m in hh.read_motifs(path)
		for name, seq in seqs.items() for pos, site, s, score in hh.scan(seq, m, strand, offset))
	assert find_rows(runs[key].stdout) == want
	assert os.listdir(fh.out(runs, key)) == ['motifFindingParameters.txt']


def test_mset_chooses_the_motif_libraries(runs):
	"""-mset worms: known enrichment over data/knownTFs/worms/known.motifs,
	de novo matches from worms/all.motifs."""

	directory = os.path.join(hh.HOMER_HOME, 'data', 'knownTFs', 'worms')
	library = hh.read_motifs(os.path.join(directory, 'known.motifs'))
	assert sorted(fh.known_rows(runs, 'mset')) == sorted(m.name for m in library)
	names = {m.name for m in hh.read_motifs(os.path.join(directory, 'all.motifs'))}
	found = hh.read_motifs(fh.out(runs, 'mset', 'nonRedundant.motifs'))
	assert found
	for motif in found:
		assert motif.name.split('BestGuess:')[1].rsplit('(', 1)[0] in names


def test_custom_genome_defaults_to_the_all_motif_set(runs):
	library = hh.read_motifs(os.path.join(hh.HOMER_HOME, 'data', 'knownTFs', 'all',
		'known.motifs'))
	rows = hh.read_known_results(fh.out(runs, 'mset_auto', 'knownResults.txt'))
	assert sorted(r['name'] for r in rows) == sorted(m.name for m in library)


# --- de novo options --------------------------------------------------------------------------

def motif_text(results, key, name='homerMotifs.motifs8'):
	return open(fh.out(results, key, name)).read()


def motifs(results, key, name='homerMotifs.motifs8'):
	return hh.read_motifs(fh.out(results, key, name))


@pytest.mark.parametrize('key,strand', [('dn', 'both'), ('dn_norevopp', '+')])
def test_de_novo_motif_statistics(runs, key, strand):
	"""Each de novo motif's T: and B: counts are the weighted sequences with a
	site of its matrix at its threshold (+ strand only with -norevopp), and its
	ln p-value follows from them."""

	seqs = fh.read_seq_tsv(fh.out(runs, key, 'seq.tsv'))
	groups = fh.read_groups(fh.out(runs, key, 'group.adj'))
	total, background_total = fh.group_totals(seqs, groups)
	found = motifs(runs, key)
	assert found
	for motif in found:
		hits, background_hits = fh.weighted_hits(seqs, groups, motif, strand)
		stats = motif.header[2]
		assert stats.startswith('T:{:.1f}('.format(hits)), (stats, hits)
		assert ',B:{:.1f}('.format(background_hits) in stats, (stats, background_hits)
		logp = fh.log_pvalue(hits, background_hits, total, background_total)
		assert float(motif.header[0]) == pytest.approx(logp, rel=1e-3, abs=1e-3)


def test_S_sets_the_number_of_motifs(runs):
	"""-S 1 finds one motif, the first of those -S 2 finds."""

	one, two = motifs(runs, 'dn_S1'), motifs(runs, 'dn')
	assert (len(one), len(two)) == (1, 2)
	assert one[0].text() == two[0].text()


def test_mis_and_oligo(runs):
	"""-oligo writes oligos.<len>.txt, the oligo table of homer2's global
	search, whose header names the -mis mismatches allowed."""

	tables = {}
	for key, mis in (('dn', 2), ('dn_mis0', 0)):
		header, rows = hh.read_tsv(fh.out(runs, key, 'oligos.8.txt'), header=True)
		assert header[3] == 'Enrichment(ln p, up to {} mismatches)'.format(mis)
		assert rows and all(len(r[0]) == 8 and set(r[0]) <= set('ACGT') for r in rows)
		tables[key] = rows
	assert tables['dn'] != tables['dn_mis0']


def direct_denovo(tmp_path, results, extra):
	"""homer2 denovo on the 'dn' run's sequence and group files with the
	options findMotifsGenome.pl passes, plus `extra`."""

	out = tmp_path / 'direct.motifs'
	hh.run(['homer2', 'denovo', '-s', fh.out(results, 'dn', 'seq.tsv'), '-g',
		fh.out(results, 'dn', 'group.adj'), '-S', 2, '-stat', 'binomial', '-mis', 2, '-cache',
		500, '-minlp', -10, '-p', 1] + extra + ['-o', out, '-len', 8], tmp_path)
	return out.read_text()


@pytest.mark.parametrize('key,extra', [('dn_e', ['-e', 1]), ('dn_olen', ['-olen', 2,
	'-omax', 160])])
def test_options_passed_to_homer2_denovo(runs, tmp_path, key, extra):
	"""-e and -olen (with -omax set to -nmax) reach homer2 denovo and change
	the motifs found."""

	assert direct_denovo(tmp_path, runs, extra) == motif_text(runs, key)
	assert motif_text(runs, key) != motif_text(runs, 'dn')


@pytest.mark.parametrize('key', ['dn_cache', 'dn_noknown'])
def test_options_that_leave_de_novo_motifs_unchanged(runs, key):
	"""-cache sizes the statistics cache and -noknown skips another step; the
	motifs are the same."""

	assert motif_text(runs, key) == motif_text(runs, 'dn')


def test_minlp_and_quick_mask(runs):
	"""After five motifs, homer2 stops at seeds with ln p above -minlp:
	-minlp -100000 stops at five of -S 7, -minlp 1 never stops. -quickMask
	skips the full masking that this stop belongs to."""

	assert len(motifs(runs, 'dn_minlp_low')) == 5
	assert len(motifs(runs, 'dn_minlp_high')) == 7
	assert len(motifs(runs, 'dn_quickmask')) == 7


def test_noknown_and_nomotif_skip_their_analysis(runs):
	"""-noknown writes no known motif results, -nomotif no de novo ones."""

	listing = set(os.listdir(fh.out(runs, 'dn_noknown')))
	assert not listing & {'knownResults.txt', 'knownResults.html', 'knownResults'}
	assert 'homerResults.html' in listing
	listing = set(os.listdir(fh.out(runs, 'base')))
	assert not [n for n in listing if n.startswith('homerMotifs') or n.startswith('homerResults')]


def test_rna_mode(runs):
	"""-rna searches the forward strand only, skips known enrichment and
	draws logos with U."""

	assert motif_text(runs, 'dn_rna') == motif_text(runs, 'dn_norevopp')
	assert not os.path.exists(fh.out(runs, 'dn_rna', 'knownResults.txt'))
	logo = open(fh.out(runs, 'dn_rna', 'homerResults', 'motif1.logo.svg')).read()
	assert '>U</text>' in logo and '>T</text>' not in logo


def test_fdr_adds_an_empirical_fdr(runs):
	"""-fdr 2 reruns homer2 denovo on two label shufflings
	(randomizations/) and appends FDR:<value> computed from their p-values."""

	real = motifs(runs, 'dn_fdr')
	assert real
	random = []
	directory = fh.out(runs, 'dn_fdr', 'randomizations')
	assert sorted(os.listdir(directory)) == ['homerMotifs.r0.motifs8', 'homerMotifs.r1.motifs8']
	for name in os.listdir(directory):
		random += [float(m.header[0]) for m in hh.read_motifs(os.path.join(directory, name))]
	want = fh.empirical_fdr([float(m.header[0]) for m in real], random, 2)
	for motif, fdr in zip(real, want):
		assert float(motif.header[2].split(',FDR:')[1]) == pytest.approx(fdr, rel=1e-9)


def test_basic_skips_the_known_motif_comparison(runs):
	"""-basic: no best matches and no similar-motif pages."""

	directory = fh.out(runs, 'dn_basic')
	names = os.listdir(os.path.join(directory, 'homerResults'))
	assert not [n for n in names if '.similar' in n]
	assert [n for n in names if n.endswith('.info.html')]
	for motif in hh.read_motifs(os.path.join(directory, 'nonRedundant.motifs')):
		assert 'BestGuess' not in motif.name


def logo_column_heights(path):
	"""Summed letter heights per column of an SVG logo (columns are 25
	units apart)."""

	heights = {}
	for part in open(path).read().split('transform="matrix(')[1:]:
		values = [float(v) for v in part.split(')')[0].split(',')]
		column = int((values[4] + 5) // 25)
		heights[column] = heights.get(column, 0.0) + values[3]
	return [heights[c] for c in sorted(heights)]


def test_bits_scales_logos_by_information_content(runs, library):
	"""-bits: a column's letters fill (2 - H) / 2 of the height, H its
	entropy in bits (probabilities raised to 0.001); without it, all of it."""

	first = hh.read_known_results(fh.out(runs, 'dn', 'knownResults.txt'))[0]['name']
	motif = [m for m in library if m.name == first][0]
	want = []
	for row in motif.rows:
		entropy = -sum(max(p, 0.001) * math.log(max(p, 0.001), 2) for p in row)
		want.append((2 - entropy) / 2)
	for key, expected in (('dn', [1.0] * len(motif)), ('dn_bits', want)):
		heights = logo_column_heights(fh.out(runs, key, 'knownResults', 'known1.logo.svg'))
		assert heights == pytest.approx(expected, abs=0.01), key


def test_nofacts_leaves_the_closing_paragraph_empty(runs):
	"""-nofacts: the closing paragraph of each page, where a fact would go,
	is empty."""

	directory = fh.out(runs, 'dn_nofacts')
	pages = [os.path.join(directory, 'homerResults.html')] + [os.path.join(directory,
		'homerResults', n) for n in os.listdir(os.path.join(directory, 'homerResults'))
		if n.endswith('.html')]
	for page in pages:
		text = open(page).read()
		assert '<P></P>' in text and text.count('<P>') == 1, page


def test_nocheck_skips_motif_comparison(runs):
	"""-nocheck skips compareMotifs.pl: de novo and known results only."""

	listing = set(os.listdir(fh.out(runs, 'dn_nocheck')))
	assert 'homerMotifs.motifs8' in listing and 'knownResults.txt' in listing
	assert not listing & {'homerResults', 'homerResults.html', 'homerMotifs.all.motifs',
		'nonRedundant.motifs'}


def test_mcheck_names_the_best_matches(runs):
	"""-mcheck <file>: the de novo motifs' best matches come from the file;
	the first motif's is siteA."""

	names = {m.name for m in hh.read_motifs(KS)}
	kept = motifs(runs, 'dn', 'nonRedundant.motifs')
	assert kept[0].name.split('BestGuess:')[1].startswith('siteA(')
	for motif in kept:
		assert motif.name.split('BestGuess:')[1].rsplit('(', 1)[0] in names

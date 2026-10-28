"""findMotifsGenome.pl with a custom genome: the files a run writes and the
numbers in them.

Two runs feed most tests: one with de novo motif finding (-len 8,10 -S 3) and
one with known motif enrichment only that keeps its intermediate files
(-dumpFasta -keepFiles). Expected values come from the inputs: the windows
are cut from tests/data/genome.fa the way resizePosFile.pl places them
(size + 1 bases), motif sites from homer_helpers.scan, known motif counts from
that scan weighted by the group file, p-values from the binomial tail, q-values
by Benjamini-Hochberg, and planted.tsv says where siteA lies.
"""

import math
import os

import pytest

import homer_helpers as hh
import findmotifsgenome_helpers as fh

DENOVO = ['-size', 200, '-len', '8,10', '-S', 3, '-N', 1000, '-mknown', fh.KNOWN_SMALL,
	'-mcheck', fh.KNOWN_SMALL]


@pytest.fixture(scope='module')
def runs(tmp_path_factory, preparsed):
	tmp = tmp_path_factory.mktemp('fmg_outputs')
	return fh.run_calls(tmp, {
		'denovo': (DENOVO, {'preparsed': preparsed}),
		'known': (['-size', 200, '-nomotif', '-mknown', fh.KNOWN_SMALL, '-dumpFasta',
			'-keepFiles'], {'preparsed': preparsed}),
	})


@pytest.fixture(scope='module')
def threaded(tmp_path_factory, preparsed):
	"""The de novo run at -p 4, alone, so that no more than four threads run."""

	tmp = tmp_path_factory.mktemp('fmg_threads')
	return fh.run_calls(tmp, {'p4': (DENOVO + ['-p', 4], {'preparsed': preparsed})}, jobs=1)


@pytest.fixture(scope='module')
def library():
	return hh.read_motifs(fh.KNOWN_SMALL)


@pytest.fixture(scope='module')
def windows(genome, peaks):
	return {p['id']: fh.window_sequence(p, 200, genome=genome) for p in peaks}


def motif_files(directory):
	return sorted(name for name in os.listdir(directory) if name.startswith('homerMotifs.motifs'))


# --- files -------------------------------------------------------------------------

def test_a_run_writes_the_documented_files_and_removes_temporary_ones(runs):
	"""Sequence, group and position files, homer2's .tmp.motifs and
	mergePeaks' coBoundBy files are deleted at the end unless -keepFiles."""

	assert sorted(os.listdir(fh.out(runs, 'denovo'))) == sorted([
		'homerMotifs.all.motifs', 'homerMotifs.motifs10', 'homerMotifs.motifs8',
		'homerResults', 'homerResults.html', 'knownResults', 'knownResults.html',
		'knownResults.txt', 'motifFindingParameters.txt', 'nonRedundant.motifs',
		'seq.autonorm.tsv'])
	assert os.listdir(runs['denovo'].cwd) == ['out']


def test_keep_files_and_dump_fasta_keep_the_intermediate_files(runs):
	listing = set(os.listdir(fh.out(runs, 'known')))
	assert {'seq.tsv', 'group.adj', 'target.fa', 'background.fa', 'target200.pos',
		'target200.seq'} <= listing
	assert not listing & {'homerMotifs.all.motifs', 'homerResults', 'homerResults.html'}


def test_motif_finding_parameters_record_the_command(runs):
	text = open(fh.out(runs, 'denovo', 'motifFindingParameters.txt')).read()
	assert text.splitlines()[0] == 'cmd = ' + ' '.join(runs['denovo'].argv[1:])


# --- de novo motifs -----------------------------------------------------------------------

def test_de_novo_motif_files_hold_motifs_of_each_length(runs):
	"""-len 8,10 writes homerMotifs.motifs8 and .motifs10; -S 3 gives three
	motifs of each, named <rank>-<consensus> and ordered by p-value."""

	directory = fh.out(runs, 'denovo')
	assert motif_files(directory) == ['homerMotifs.motifs10', 'homerMotifs.motifs8']
	for length in (8, 10):
		motifs = hh.read_motifs(os.path.join(directory, 'homerMotifs.motifs{}'.format(length)))
		assert len(motifs) == 3
		logps = [float(m.header[0]) for m in motifs]
		assert logps == sorted(logps)
		for rank, motif in enumerate(motifs, 1):
			assert len(motif) == length
			assert len(motif.consensus) == length
			assert motif.name == '{}-{}'.format(rank, motif.consensus)
			for row in motif.rows:
				assert sum(row) == pytest.approx(1.0, abs=0.01)


def test_all_motifs_file_concatenates_the_length_files(runs):
	"""homerMotifs.all.motifs is cat homerMotifs.motifs*, in the shell's
	sorted order."""

	directory = fh.out(runs, 'denovo')
	joined = ''.join(open(os.path.join(directory, name)).read() for name in motif_files(directory))
	assert open(os.path.join(directory, 'homerMotifs.all.motifs')).read() == joined


@pytest.mark.parametrize('length', [8, 10])
def test_de_novo_recovers_siteA(runs, genome, length):
	"""The best motif of each length has a site over at least 110 of the 120
	planted siteA copies."""

	directory = fh.out(runs, 'denovo')
	best = hh.read_motifs(os.path.join(directory, 'homerMotifs.motifs{}'.format(length)))[0]
	assert len(fh.peaks_with_overlapping_site(best, 'siteA', genome=genome)) >= 110


def test_non_redundant_motifs_account_for_every_de_novo_motif(runs):
	"""compareMotifs.pl keeps each motif either in nonRedundant.motifs or as a
	homerResults/motif<i>.similar<j>.motif of the motif it merged into."""

	directory = fh.out(runs, 'denovo')
	everything = hh.read_motifs(os.path.join(directory, 'homerMotifs.all.motifs'))
	kept = hh.read_motifs(os.path.join(directory, 'nonRedundant.motifs'))
	results = os.path.join(directory, 'homerResults')
	similar = [hh.read_motifs(os.path.join(results, name))[0] for name in sorted(os.listdir(results))
		if '.similar' in name and name.endswith('.motif')]
	assert sorted(m.consensus for m in kept + similar) == sorted(m.consensus for m in everything)
	assert len(kept) >= 1


def test_homer_results_pages_name_each_motif(runs):
	"""For each non-redundant motif i: homerResults/motif<i>.motif holds it,
	motif<i>RV.motif its reverse complement, logos and an information page
	exist, and homerResults.html links the motif and names its best match
	in the -mcheck library."""

	directory = fh.out(runs, 'denovo')
	kept = hh.read_motifs(os.path.join(directory, 'nonRedundant.motifs'))
	html = open(os.path.join(directory, 'homerResults.html')).read()
	names = {m.name for m in hh.read_motifs(fh.KNOWN_SMALL)}
	for i, motif in enumerate(kept, 1):
		results = os.path.join(directory, 'homerResults')
		page = hh.read_motifs(os.path.join(results, 'motif{}.motif'.format(i)))[0]
		assert page.consensus == motif.consensus
		assert page.rows == motif.rows
		reverse = hh.read_motifs(os.path.join(results, 'motif{}RV.motif'.format(i)))[0]
		assert reverse.rows == [row[::-1] for row in motif.rows[::-1]]
		for suffix in ('.logo.svg', 'RV.logo.svg', '.info.html'):
			assert os.path.exists(os.path.join(results, 'motif{}{}'.format(i, suffix)))
		assert 'homerResults/motif{}.motif'.format(i) in html
		guess = motif.name.split('BestGuess:')[1].rsplit('(', 1)[0]
		assert guess in names
		assert guess in html
	assert kept[0].name.split('BestGuess:')[1].startswith('siteA(')


# --- known motifs ---------------------------------------------------------------------------

def test_target_fasta_holds_the_peak_windows(runs, windows):
	"""-dumpFasta's target.fa: the size + 1 bases around each peak centre,
	reverse-complemented for - strand peaks, in upper case."""

	got = fh.read_fasta_upper(fh.out(runs, 'known', 'target.fa'))
	assert got == {name: seq.upper() for name, seq in windows.items()}
	assert {len(seq) for seq in got.values()} == {201}


def test_known_target_counts_equal_a_reference_scan(runs, windows, library):
	rows = {r['name']: r for r in hh.read_known_results(fh.out(runs, 'known', 'knownResults.txt'))}
	for motif in library:
		want = sum(1 for seq in windows.values() if fh.has_site(seq, motif))
		assert rows[motif.name]['targets'] == want, motif.name
	assert rows['siteA']['targets'] == 120
	assert rows['siteA']['target_total'] == 200


def test_known_background_counts_follow_the_weights(runs, library):
	"""Background counts are the summed group.adj weights of the
	background.fa sequences with a site."""

	background = fh.read_fasta_upper(fh.out(runs, 'known', 'background.fa'))
	groups = fh.read_groups(fh.out(runs, 'known', 'group.adj'))
	assert set(background) == {name for name, (g, _) in groups.items() if g == 0}
	weights = {groups[name][1] for name in background}
	assert len(weights) > 10, 'the GC-matched background should carry unequal weights'
	rows = {r['name']: r for r in hh.read_known_results(fh.out(runs, 'known', 'knownResults.txt'))}
	for motif in library:
		want = sum(groups[name][1] for name, seq in background.items() if fh.has_site(seq, motif))
		assert rows[motif.name]['background'] == pytest.approx(want, abs=0.051), motif.name


def test_sequence_file_holds_targets_and_background(runs):
	seqs = fh.read_seq_tsv(fh.out(runs, 'known', 'seq.tsv'))
	target = fh.read_fasta_upper(fh.out(runs, 'known', 'target.fa'))
	background = fh.read_fasta_upper(fh.out(runs, 'known', 'background.fa'))
	assert seqs == dict(target, **background)


def test_known_results_columns_follow_from_the_counts(runs, library):
	"""Percentages, binomial p-values, Benjamini q-values and the ranking,
	all from the weighted reference counts."""

	directory = fh.out(runs, 'known')
	fh.check_known_results(os.path.join(directory, 'knownResults.txt'),
		fh.read_seq_tsv(os.path.join(directory, 'seq.tsv')),
		fh.read_groups(os.path.join(directory, 'group.adj')), library)


def test_known_results_html_and_motif_files(runs, library):
	"""Each motif with p <= 0.01 gets a row naming it in knownResults.html,
	and knownResults/known<rank>.motif holding the library matrix and
	threshold, with a logo; the others get neither."""

	directory = fh.out(runs, 'known')
	rows = hh.read_known_results(os.path.join(directory, 'knownResults.txt'))
	html = open(os.path.join(directory, 'knownResults.html')).read()
	by_name = {m.name: m for m in library}
	significant = [row for row in rows if math.exp(row['logp']) <= 0.01]
	assert [row['name'] for row in significant] == ['siteA', 'siteB']
	for rank, row in enumerate(significant, 1):
		path = os.path.join(directory, 'knownResults', 'known{}.motif'.format(rank))
		motif = hh.read_motifs(path)[0]
		assert motif.name == row['name']
		assert motif.rows == by_name[row['name']].rows
		assert motif.threshold == pytest.approx(by_name[row['name']].threshold)
		assert os.path.exists(path.replace('.motif', '.logo.svg'))
		assert '<TD>{}</TD>'.format(row['name']) in html
		assert 'known{}.motif'.format(rank) in html
	extra = len(significant) + 1
	assert not os.path.exists(os.path.join(directory, 'knownResults', 'known{}.motif'.format(extra)))
	for row in rows[len(significant):]:
		assert '<TD>{}</TD>'.format(row['name']) not in html


def test_autonorm_reports_the_target_oligo_frequencies(runs):
	"""seq.autonorm.tsv lists every 1-, 2- and 3-mer (-nlen 3) with the
	targets' frequency over both strands, each sequence's last k-mer counted
	or not."""

	directory = fh.out(runs, 'known')
	header, rows = hh.read_tsv(os.path.join(directory, 'seq.autonorm.tsv'), header=True)
	assert header[0] == 'Oligo'
	assert sorted(len(row[0]) for row in rows) == [1] * 4 + [2] * 16 + [3] * 64
	targets = [(seq, 1.0) for seq in fh.read_fasta_upper(os.path.join(directory,
		'target.fa')).values()]
	texts = fh.autonorm_target_texts(targets, (1, 2, 3))
	for row in rows:
		assert row[1] in texts.get(row[0], {'0.00%'}), row


# --- threads ----------------------------------------------------------------------------------

def test_threads_give_the_same_results(runs, threaded, genome):
	"""-p 4 changes no known motif's counts. P-values computed in threads go
	through homer2's shared statistics cache, and de novo motif finding at
	-p 4 can differ in detail from -p 1 (its trial count follows -p), so only
	the known counts and percentages, the files, the motif counts and the
	recovery of siteA are compared."""

	one, four = fh.out(runs, 'denovo'), fh.out(threaded, 'p4')
	columns = ('targets', 'target_percent', 'background', 'background_percent')
	counts = [{name: tuple(row[c] for c in columns) for name, row in
		fh.known_rows(results, key).items()} for results, key in ((runs, 'denovo'),
		(threaded, 'p4'))]
	assert counts[0] == counts[1]
	assert sorted(os.listdir(four)) == sorted(os.listdir(one))
	for name in motif_files(one):
		a = hh.read_motifs(os.path.join(one, name))
		b = hh.read_motifs(os.path.join(four, name))
		assert len(a) == len(b)
		assert len(fh.peaks_with_overlapping_site(b[0], 'siteA', genome=genome)) >= 110

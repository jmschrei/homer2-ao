"""homer2 denovo: de novo motif discovery in target sequences against
background sequences.

The inputs are tests/data/targets.fa (200 peaks; siteA TAGCGTACCGTA is
planted in 120 of them, siteB in 80, REST in 60) and background.fa (1,000
genome windows with no planted site), as one sequence/group file pair or as
the two FASTA files. Every run happens once, in a module fixture, at -p 1
unless a test is about -p, with -len 8 or less so that the module is fast.

Expected values are derived independently:

- a motif's statistics: T and B in its header are the numbers of target and
  background sequences with at least one window scoring at least the written
  threshold, counted with a Python scan of the written matrix; the log
  p-value is the binomial tail ln P(X >= T), X ~ Binomial(200, B / 1000)
  (B / 1000 raised to 1 / 1000), or the hypergeometric tail with
  -stat hypergeo;
- a threshold that HOMER chooses on whole sequences is the best-scoring
  sequence-level threshold, found by a Python search over each sequence's
  best window score;
- recovery of a planted site: the motif's consensus (IUPAC) aligns inside the
  site or its reverse complement with at most one mismatch, and its sites
  cover the planted copies;
- -oligos counts and -oout composition are counted in Python from the input
  sequences;
- options documented as defaults, or as affecting only memory, give the same
  file as the run without them.
"""

import collections
import itertools
import os
import sys

import pytest

import homer_helpers as hh


SITE_A = 'TAGCGTACCGTA'
SITE_B = 'GGATCCTTAC'
N_TARGETS = 200
N_BACKGROUND = 1000

IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG',
	'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
	'N': 'ACGT'}
IUPAC_CODE = {frozenset(bases): code for code, bases in IUPAC.items()}

# key -> arguments after the inputs. Every run writes its motifs to
# out.motifs except 'stdout', which writes them to stdout.
LEN8 = ['-len', 8]
RUNS = {
	'default': LEN8 + ['-S', 3],
	'stdout': LEN8 + ['-S', 3],
	'fasta': LEN8 + ['-S', 3],
	'full_mask': LEN8 + ['-S', 3, '-fullMask'],
	'strand_both': LEN8 + ['-S', 3, '-strand', 'both'],
	'binomial': LEN8 + ['-S', 3, '-stat', 'binomial'],
	'mis2': LEN8 + ['-S', 3, '-mis', 2],
	'trials10': LEN8 + ['-S', 3, '-T', 10],
	'max_back_half': LEN8 + ['-S', 3, '-maxBack', 0.5],
	'minlp_default': LEN8 + ['-S', 3, '-minlp', -10],
	'cache0': LEN8 + ['-S', 3, '-cache', 0],
	'cache1': LEN8 + ['-S', 3, '-cache', 1],
	'len6': ['-len', 6, '-S', 2],
	'len10': ['-len', 10, '-S', 1, '-mis', 1],
	'one_motif': LEN8 + ['-S', 1],
	'hypergeo': LEN8 + ['-S', 3, '-stat', 'hypergeo'],
	'hypergeo_cache0': LEN8 + ['-S', 3, '-stat', 'hypergeo', '-cache', 0],
	'plus': LEN8 + ['-S', 3, '-strand', '+', '-oligos', 'oligos.txt'],
	'minus': LEN8 + ['-S', 3, '-strand', '-'],
	'oligos': LEN8 + ['-S', 1, '-oligos', 'oligos.txt'],
	'oligos_mis0': LEN8 + ['-S', 1, '-mis', 0, '-oligos', 'oligos.txt'],
	'olen': LEN8 + ['-S', 1, '-olen', 2, '-oout', 'norm.txt', '-oligos', 'oligos.txt'],
	'omax1': LEN8 + ['-S', 1, '-olen', 2, '-omax', 1, '-oout', 'norm.txt',
		'-oligos', 'oligos.txt'],
	'threads2': LEN8 + ['-S', 2, '-p', 2],
	'threads4': LEN8 + ['-S', 2, '-p', 4],
	'quick': LEN8 + ['-S', 3, '-quickMask'],
	'quick_nozoops': LEN8 + ['-S', 3, '-quickMask', '-nozoops'],
	'quick_blen0': LEN8 + ['-S', 3, '-quickMask', '-blen', 0],
	'quick_blen1': LEN8 + ['-S', 3, '-quickMask', '-blen', 1],
	'quick_blen3': LEN8 + ['-S', 3, '-quickMask', '-blen', 3],
	'expect_strict': LEN8 + ['-S', 3, '-e', 0.001],
	'expect_loose': LEN8 + ['-S', 3, '-e', 0.1],
	'trials1': LEN8 + ['-S', 3, '-T', 1],
	'max_back0': LEN8 + ['-S', 3, '-maxBack', 0],
	'minlp_lenient': LEN8 + ['-S', 8, '-minlp', 0],
	'minlp_strict': LEN8 + ['-S', 8, '-minlp', -1000],
	'tmp_default': ['-len', 6, '-S', 1],
	'tmp_option': ['-len', 6, '-S', 1, '-tmp', os.path.join('sub', 'results.tmp')],
}
SENTINEL = 'not written by homer2\n'


@pytest.fixture(scope='module')
def records():
	return hh.read_fasta(hh.TARGETS), hh.read_fasta(hh.BACKGROUND)


@pytest.fixture(scope='module')
def runs(tmp_path_factory, records):
	"""Every run of the module, made at once: key -> (cwd, CompletedProcess)."""

	targets, background = records
	top = tmp_path_factory.mktemp('denovo')
	seq = hh.write_seq_tsv(str(top / 'seq.tsv'), targets + background)
	group = hh.write_group(str(top / 'group.txt'),
		[(name, 1, 1) for name, _ in targets] + [(name, 0, 1) for name, _ in background])
	calls = {}
	for key, extra in RUNS.items():
		cwd = top / key
		os.makedirs(str(cwd))
		if key == 'fasta':
			inputs = ['-i', hh.TARGETS, '-b', hh.BACKGROUND]
		else:
			inputs = ['-s', seq, '-g', group]
		output = [] if key == 'stdout' else ['-o', 'out.motifs']
		calls[key] = (['homer2', 'denovo'] + inputs + extra + output, cwd)
	for key in ('tmp_default', 'tmp_option'):
		(top / key / '.tmp.motifs').write_text(SENTINEL)
	os.makedirs(str(top / 'tmp_option' / 'sub'))
	(top / 'tmp_option' / 'sub' / 'results.tmp').write_text(SENTINEL)
	threaded = {key: call for key, call in calls.items() if '-p' in RUNS[key]}
	results = hh.run_many({key: call for key, call in calls.items() if key not in threaded})
	results.update(hh.run_many(threaded, jobs=1))
	return {key: (top / key, results[key]) for key in calls}


def motif_text(runs, key):
	cwd, result = runs[key]
	if key == 'stdout':
		return result.stdout
	return (cwd / 'out.motifs').read_text()


def motifs_of(runs, key):
	cwd, result = runs[key]
	if key == 'stdout':
		path = cwd / 'stdout.motifs'
		path.write_text(result.stdout)
		return hh.read_motifs(str(path))
	return hh.read_motifs(str(cwd / 'out.motifs'))


def run_strand(key):
	extra = RUNS[key]
	return extra[extra.index('-strand') + 1] if '-strand' in extra else 'both'


def run_stat(key):
	return 'hypergeo' if 'hypergeo' in RUNS[key] else 'binomial'


# --- reference computations ----------------------------------------------------------

def log_p(stat, t, b):
	"""HOMER's enrichment of t of 200 targets against b of 1,000 background
	sequences: the binomial tail with the background rate raised to 1/1000, or
	the hypergeometric tail."""

	if stat == 'hypergeo':
		return hh.log_hypergeometric_tail(t, t + b, N_TARGETS, N_TARGETS + N_BACKGROUND)
	rate = min(max(b / float(N_BACKGROUND), 1.0 / N_BACKGROUND),
		(N_BACKGROUND - 1.0) / N_BACKGROUND)
	return hh.log_binomial_tail(t, N_TARGETS, rate)


_BEST = {}


def best_scores(records, motif, strand):
	"""Best window score of every target, then of every background sequence
	(None for a sequence with no scorable window)."""

	key = (motif.text(), strand)
	if key not in _BEST:
		targets, background = records
		_BEST[key] = ([hh.best_score(s, motif, strand) for _, s in targets],
			[hh.best_score(s, motif, strand) for _, s in background])
	return _BEST[key]


def counts_at_threshold(records, motif, strand):
	targets, background = best_scores(records, motif, strand)
	t = sum(1 for s in targets if s is not None and s >= motif.threshold)
	b = sum(1 for s in background if s is not None and s >= motif.threshold)
	return t, b


def optimal_threshold(records, motif, strand, stat):
	"""(log p, threshold score, targets, background) of the best sequence-level
	threshold: sequences are taken in order of their best window score, each
	distinct score is a candidate, and the search stops once more than half
	of the background is included."""

	targets, background = best_scores(records, motif, strand)
	scored = sorted([(s, 1) for s in targets if s is not None]
		+ [(s, 0) for s in background if s is not None], reverse=True)
	best = None
	t = b = 0
	for i in range(len(scored) - 1):
		score, group = scored[i]
		t += group
		b += 1 - group
		if score < scored[i + 1][0] + 1e-6:
			continue
		lp = log_p(stat, t, b)
		if best is None or lp < best[0]:
			best = (lp, score, t, b)
		if b > 0.5 * N_BACKGROUND:
			break
	return best


def header_stats(motif):
	"""(log p, T, T%, B, B%, P exponent) from a motif header."""

	stats = motif.header[2]
	parts = dict(part.split(':', 1) for part in stats.split(','))
	t, t_pct = parts['T'].rstrip('%)').split('(')
	b, b_pct = parts['B'].rstrip('%)').split('(')
	return (float(motif.header[0]), float(t), float(t_pct), float(b), float(b_pct),
		int(parts['P'][2:]))


def site_mismatches(consensus, site):
	"""Fewest mismatches of the IUPAC consensus placed inside the site or its
	reverse complement (a consensus longer than the site may overhang it by
	the difference)."""

	best = len(consensus)
	for s in (site, hh.revcomp(site)):
		if len(consensus) <= len(s):
			placements = [(0, s[k:k + len(consensus)])
				for k in range(len(s) - len(consensus) + 1)]
		else:
			placements = [(k, s) for k in range(len(consensus) - len(s) + 1)]
		for k, part in placements:
			cons = consensus[k:k + len(part)]
			best = min(best, sum(1 for c, base in zip(cons, part) if base not in IUPAC[c]))
	return best


@pytest.fixture(scope='module')
def planted_spans(records):
	"""peak -> (start, end) of its planted siteA copy, 0-based in targets.fa."""

	peaks = {p['id']: p for p in hh.read_peaks()}
	spans = {}
	for p in hh.read_planted():
		if p['motif'] == 'siteA':
			start = int(p['start']) - (peaks[p['peak']]['start'] - 1)
			spans[p['peak']] = (start, start + len(SITE_A))
	assert len(spans) == 120
	return spans


def covered_sites(records, motif, spans, strand='both'):
	"""How many planted siteA copies a site of the motif overlaps."""

	targets = dict(records[0])
	width = len(motif)
	covered = 0
	for peak, (start, end) in spans.items():
		for pos, _, site_strand, _ in hh.scan(targets[peak], motif, strand, offset=0):
			first = pos if site_strand == '+' else pos - width + 1
			if first < end and start < first + width:
				covered += 1
				break
	return covered


def assert_recovers_site_a(records, motif, spans, strand='both'):
	assert site_mismatches(motif.consensus, SITE_A) <= 1, motif.consensus
	assert covered_sites(records, motif, spans, strand) >= 108, motif.consensus


def assert_header_matches_scan(records, motif, strand, stat):
	lp, t, t_pct, b, b_pct, exponent = header_stats(motif)
	want_t, want_b = counts_at_threshold(records, motif, strand)
	assert (t, b) == (want_t, want_b), motif.consensus
	assert t_pct == pytest.approx(100.0 * t / N_TARGETS, abs=0.006)
	assert b_pct == pytest.approx(100.0 * b / N_BACKGROUND, abs=0.006)
	assert lp == pytest.approx(log_p(stat, want_t, want_b), abs=1e-3), motif.consensus
	assert exponent == int(lp / 2.3026)


def consensus_of(rows):
	"""The IUPAC symbol of the bases within 0.15 of each position's largest
	probability; None for a position too close to that cut to decide from
	three-decimal values."""

	symbols = []
	for row in rows:
		top = max(row)
		gaps = [top - p for p in row]
		if any(abs(g - 0.15) < 0.002 for g in gaps):
			symbols.append(None)
			continue
		symbols.append(IUPAC_CODE[frozenset(b for b, g in zip('ACGT', gaps) if g < 0.15)])
	return symbols


# --- the motif file ------------------------------------------------------------------------

@pytest.mark.parametrize('key', sorted(k for k in RUNS if k != 'stdout'))
def test_motif_file_is_well_formed(runs, key):
	"""Each motif has -len rows of four probabilities summing to 1, a
	consensus that is the IUPAC code of each row, the name <rank>-<consensus>,
	and motifs come in order of significance, at most -S of them."""

	extra = RUNS[key]
	length = int(extra[extra.index('-len') + 1])
	most = int(extra[extra.index('-S') + 1])
	motifs = motifs_of(runs, key)
	assert 1 <= len(motifs) <= most
	for rank, m in enumerate(motifs, 1):
		assert len(m) == length and len(m.consensus) == length
		assert m.name == '{}-{}'.format(rank, m.consensus)
		for row in m.rows:
			assert len(row) == 4 and all(0.0 <= p <= 1.0 for p in row)
			assert sum(row) == pytest.approx(1.0, abs=0.003)
		for symbol, c in zip(consensus_of(m.rows), m.consensus):
			assert symbol is None or symbol == c, (m.consensus, m.rows)
	logps = [float(m.header[0]) for m in motifs]
	assert logps == sorted(logps)


@pytest.mark.parametrize('key', ['default', 'len6', 'len10', 'hypergeo', 'plus', 'minus',
	'quick', 'quick_nozoops', 'threads4', 'minlp_lenient', 'trials1', 'olen'])
def test_header_statistics_match_a_scan(runs, records, key):
	"""T and B are the target and background sequences with a site at the
	written threshold on the searched strands, and the log p-value is the
	statistic's tail probability of those counts."""

	for m in motifs_of(runs, key):
		assert_header_matches_scan(records, m, run_strand(key), run_stat(key))


@pytest.mark.parametrize('key', ['default', 'one_motif', 'len6', 'len10', 'hypergeo',
	'plus', 'minus', 'quick', 'threads2', 'threads4', 'expect_strict', 'expect_loose',
	'trials1', 'max_back0', 'fasta'])
def test_top_motif_is_the_planted_site(runs, records, planted_spans, key):
	"""The most significant motif is siteA, planted in 120 of 200 targets and
	in none of the background sequences."""

	top = motifs_of(runs, key)[0]
	assert_recovers_site_a(records, top, planted_spans, run_strand(key))


def test_second_motif_is_site_b(runs):
	motifs = motifs_of(runs, 'default')
	assert site_mismatches(motifs[1].consensus, SITE_B) <= 1, motifs[1].consensus


@pytest.mark.parametrize('key,strand,stat', [('default', 'both', 'binomial'),
	('hypergeo', 'both', 'hypergeo'), ('plus', '+', 'binomial')])
def test_first_motif_threshold_is_optimal(runs, records, key, strand, stat):
	"""The first motif found is scored on the unmasked sequences, so its
	threshold is the sequence-level threshold with the best enrichment, less
	1e-6."""

	top = motifs_of(runs, key)[0]
	lp, score, t, b = optimal_threshold(records, top, strand, stat)
	assert top.threshold == pytest.approx(score - 1e-6, abs=1e-4)
	assert float(top.header[0]) == pytest.approx(lp, abs=1e-3)


def test_pinned_default_motifs(runs):
	"""pins the current output: the consensus, length and target and
	background counts of each motif of -len 8 -S 3 at -p 1."""

	got = [(m.consensus, len(m), header_stats(m)[1], header_stats(m)[3])
		for m in motifs_of(runs, 'default')]
	# Two candidates for the third motif tie on every statistic, and qsort
	# picks between them; macOS's qsort picks the other one, in v5.1 too.
	third = 'CATGGTGC' if sys.platform == 'darwin' else 'GGAGCTGT'
	assert got == [('TASCGTAC', 8, 121.0, 6.0), ('TAAGGATC', 8, 80.0, 3.0),
		(third, 8, 61.0, 2.0)]


# --- inputs and outputs -----------------------------------------------------------------

def test_fasta_inputs_give_the_same_motifs(runs):
	"""-i/-b reads the two FASTA files as targets and background, the same
	as a sequence file with a group file that marks them 1 and 0."""

	assert motif_text(runs, 'fasta') == motif_text(runs, 'default')


def test_fasta_inputs_leave_no_files(runs):
	"""The sequence and group files made from -i/-b and the temporary results
	file are removed at exit."""

	cwd, _ = runs['fasta']
	assert sorted(os.listdir(str(cwd))) == ['out.motifs']


def test_stdout_without_o(runs):
	cwd, result = runs['default']
	assert result.stdout == ''
	assert runs['stdout'][1].stdout == motif_text(runs, 'default')


def test_tmp_file_is_removed_and_default_left_alone(runs):
	"""-tmp names the temporary results file: a file at that path is
	replaced and removed at exit, and .tmp.motifs in the working directory is
	not touched. Without -tmp, .tmp.motifs is the one replaced and removed."""

	default_cwd, _ = runs['tmp_default']
	assert not (default_cwd / '.tmp.motifs').exists()
	option_cwd, _ = runs['tmp_option']
	assert (option_cwd / '.tmp.motifs').read_text() == SENTINEL
	assert not (option_cwd / 'sub' / 'results.tmp').exists()
	assert motif_text(runs, 'tmp_option') == motif_text(runs, 'tmp_default')


@pytest.mark.parametrize('key', ['full_mask', 'strand_both', 'binomial', 'mis2', 'trials10',
	'max_back_half', 'minlp_default'])
def test_documented_defaults(runs, key):
	"""-fullMask, -strand both, -stat binomial, -mis 2, -T 10, -maxBack 0.5
	and -minlp -10 are the documented defaults."""

	assert motif_text(runs, key) == motif_text(runs, 'default')


@pytest.mark.parametrize('key,same_as', [('cache0', 'default'), ('cache1', 'default'),
	('hypergeo_cache0', 'hypergeo')])
def test_cache_size_changes_only_memory(runs, key, same_as):
	assert motif_text(runs, key) == motif_text(runs, same_as)


# --- options ----------------------------------------------------------------------------

def test_one_motif(runs):
	assert len(motifs_of(runs, 'one_motif')) == 1
	assert len(motifs_of(runs, 'default')) == 3


def test_motif_length(runs):
	for key, length in (('len6', 6), ('default', 8), ('len10', 10)):
		assert {len(m) for m in motifs_of(runs, key)} == {length}


def test_minus_strand_motif_is_the_reverse_complement(runs):
	"""-strand - searches the reverse complement of every sequence, so its
	top motif is the -strand + top motif reverse-complemented, with the same
	threshold and statistics."""

	plus = motifs_of(runs, 'plus')[0]
	minus = motifs_of(runs, 'minus')[0]
	assert minus.consensus == hh.revcomp(plus.consensus)
	assert minus.rows == [row[::-1] for row in plus.rows[::-1]]
	assert minus.threshold == plus.threshold
	assert minus.header[:3] == plus.header[:3]


def test_hypergeo_differs_from_binomial(runs):
	assert motif_text(runs, 'hypergeo') != motif_text(runs, 'default')


def test_threads_find_the_planted_sites(runs, records, planted_spans):
	"""At -p 2 and 4 the two planted sites are still the two motifs; their
	statistics agree with a scan (the run itself may vary slightly with
	thread timing)."""

	for key in ('threads2', 'threads4'):
		motifs = motifs_of(runs, key)
		assert len(motifs) == 2
		assert_recovers_site_a(records, motifs[0], planted_spans)
		assert site_mismatches(motifs[1].consensus, SITE_B) <= 1
		for m in motifs:
			assert_header_matches_scan(records, m, 'both', 'binomial')


def test_quick_mask_differs_from_full_mask(runs):
	assert motif_text(runs, 'quick') != motif_text(runs, 'default')


def test_quick_mask_rescores_thresholds_on_whole_sequences(runs, records):
	"""With -quickMask the final zoops scoring chooses every motif's threshold
	again, on the whole sequences."""

	for m in motifs_of(runs, 'quick'):
		lp, score, t, b = optimal_threshold(records, m, 'both', 'binomial')
		assert m.threshold == pytest.approx(score - 1e-6, abs=1e-4), m.consensus
		assert float(m.header[0]) == pytest.approx(lp, abs=1e-3)


def test_nozoops_keeps_the_oligo_thresholds(runs, records):
	"""-nozoops skips the final zoops scoring: the matrices are those
	of the -quickMask run, but at least one threshold is not the
	sequence-level optimum, and the statistics still describe the written
	threshold."""

	zoops = motifs_of(runs, 'quick')
	nozoops = motifs_of(runs, 'quick_nozoops')
	assert sorted(m.rows for m in nozoops) == sorted(m.rows for m in zoops)
	off = 0
	for m in nozoops:
		lp, score, t, b = optimal_threshold(records, m, 'both', 'binomial')
		if abs(m.threshold - (score - 1e-6)) > 1e-4:
			off += 1
		assert_header_matches_scan(records, m, 'both', 'binomial')
	assert off >= 1


def test_blen_one_is_the_default_at_length_8(runs):
	assert motif_text(runs, 'quick_blen1') == motif_text(runs, 'quick')


def test_blen_skips_seeds_next_to_found_motifs(runs, records, planted_spans):
	"""-blen is how far a seed may be shifted against an oligo already bound
	by a motif and still be skipped as redundant: with -blen 0 more of the
	motifs are shifted copies of siteA than with -blen 3."""

	def site_a_motifs(key):
		return sum(1 for m in motifs_of(runs, key)
			if covered_sites(records, m, planted_spans) >= 108)

	assert site_a_motifs('quick_blen0') > site_a_motifs('quick_blen3') >= 1


@pytest.mark.parametrize('key,other', [('expect_strict', 'expect_loose'),
	('trials1', 'default'), ('max_back0', 'default')])
def test_search_settings_change_the_motifs(runs, key, other):
	"""-e, -T and -maxBack change how matrices are refined: a strict
	setting gives other motifs than the default or a loose one (each still
	finds siteA, tested above)."""

	assert motif_text(runs, key) != motif_text(runs, other)


def test_minlp_stops_at_unpromising_seeds(runs):
	"""Seeds less significant than -minlp are not refined into motifs once
	the first motifs are found: -minlp -1000 gives fewer motifs than -minlp 0
	at the same -S."""

	lenient = motifs_of(runs, 'minlp_lenient')
	strict = motifs_of(runs, 'minlp_strict')
	assert len(lenient) == 8
	assert 1 <= len(strict) < len(lenient)


# --- -oligos, -olen, -oout, -omax ------------------------------------------------------------

def window_counts(records, length):
	"""(target, background) Counters of every length-mer window holding only
	A, C, G and T, case ignored."""

	out = []
	for recs in records:
		c = collections.Counter()
		for _, seq in recs:
			seq = seq.upper()
			for i in range(len(seq) - length + 1):
				w = seq[i:i + length]
				if all(b in 'ACGT' for b in w):
					c[w] += 1
		out.append(c)
	return out


def read_oligos(runs, key):
	cwd, _ = runs[key]
	header, rows = hh.read_tsv(str(cwd / 'oligos.txt'), header=True)
	return header, [(r[0], float(r[1]), float(r[2]), float(r[3])) for r in rows]


@pytest.fixture(scope='module')
def counts8(records):
	return window_counts(records, 8)


@pytest.mark.parametrize('key,strand,mis', [('oligos', 'both', 2), ('oligos_mis0', 'both', 0),
	('plus', '+', 2)])
def test_oligo_table_counts(runs, counts8, key, strand, mis):
	"""-oligos lists every distinct 8-mer once with its occurrences in the
	targets and background; on both strands an oligo and its reverse
	complement are one entry. Rows are in order of enrichment."""

	header, rows = read_oligos(runs, key)
	assert header[3] == 'Enrichment(ln p, up to {} mismatches)'.format(mis)
	target, background = counts8

	def count(c, oligo):
		if strand == '+':
			return c[oligo]
		rc = hh.revcomp(oligo)
		return c[oligo] + (c[rc] if rc != oligo else 0)

	seen = set()
	for oligo, t, b, _ in rows:
		assert (t, b) == (count(target, oligo), count(background, oligo)), oligo
		key_ = oligo if strand == '+' else min(oligo, hh.revcomp(oligo))
		assert key_ not in seen
		seen.add(key_)
	every = set(target) | set(background)
	if strand != '+':
		every = {min(o, hh.revcomp(o)) for o in every}
	assert seen == every
	assert sum(r[1] for r in rows) == N_TARGETS * 193
	assert sum(r[2] for r in rows) == N_BACKGROUND * 193
	values = [r[3] for r in rows]
	assert values == sorted(values)


def zoops_tables(n_targets, n_background):
	"""HOMER's approximate number of sequences holding n occurrences:
	a(1) = 1, a(n) = a(n - 1) + (N - a(n - 1)) / N."""

	size = int(5 * (n_targets + n_background))
	tables = []
	for total in (n_targets, n_background):
		a = [0.0, 1.0]
		for _ in range(2, size):
			a.append(a[-1] + (total - a[-1]) / total)
		tables.append(a)
	return tables


def test_oligo_enrichment_without_mismatches(runs):
	"""At -mis 0 an oligo's value is the binomial tail of its occurrences,
	each count first turned into an approximate number of sequences, against
	1,000 background sequences (the background is 1,000 sequences of the
	targets' length)."""

	_, rows = read_oligos(runs, 'oligos_mis0')
	ta, ba = zoops_tables(float(N_TARGETS), float(N_BACKGROUND))
	checked = 0
	for oligo, t, b, value in rows[:60]:
		k = int(ta[int(t)])
		rate = min(max(int(ba[int(b)]) / float(N_BACKGROUND), 1.0 / N_BACKGROUND),
			(N_BACKGROUND - 1.0) / N_BACKGROUND)
		assert value == pytest.approx(hh.log_binomial_tail(k, N_TARGETS, rate), abs=2e-3), oligo
		checked += 1
	assert checked == 60


def test_mismatches_only_make_oligos_more_significant(runs):
	"""With -mis 2 an oligo keeps the best value over 0, 1 and 2 mismatches,
	so no value is higher than at -mis 0, and some are lower."""

	_, zero = read_oligos(runs, 'oligos_mis0')
	_, two = read_oligos(runs, 'oligos')
	zero = {r[0]: r for r in zero}
	two = {r[0]: r for r in two}
	assert set(zero) == set(two)
	lower = 0
	for oligo, row in two.items():
		assert row[1:3] == zero[oligo][1:3]
		assert row[3] <= zero[oligo][3] + 1e-9, oligo
		if row[3] < zero[oligo][3] - 1e-3:
			lower += 1
	assert lower > 0


def test_olen_normalizes_background_counts(runs):
	"""-olen reweights background oligos to match the targets' lower-order
	composition: target counts are unchanged, background counts are not."""

	_, plain = read_oligos(runs, 'oligos')
	_, normalized = read_oligos(runs, 'olen')
	plain = {r[0]: r for r in plain}
	normalized = {r[0]: r for r in normalized}
	assert set(plain) == set(normalized)
	assert all(normalized[o][1] == plain[o][1] for o in plain)
	assert any(abs(normalized[o][2] - plain[o][2]) > 0.05 for o in plain)


def test_omax_limits_normalization_rounds(runs):
	"""-omax caps the number of reweighting rounds: with -omax 1 the
	background counts differ from those of the run with the default -omax 160,
	and the target counts are the same."""

	_, many = read_oligos(runs, 'olen')
	_, one = read_oligos(runs, 'omax1')
	many = {r[0]: r for r in many}
	assert any(abs(r[2] - many[r[0]][2]) > 0.05 for r in one)
	assert all(r[1] == many[r[0]][1] for r in one)


@pytest.fixture(scope='module')
def word_shares(records):
	"""word -> (target %, background %) of every 1-mer and 2-mer in the
	8-mer windows: each sub-word of a window counts 1 / (9 - k), and a word's
	share is averaged with its reverse complement's."""

	shares = {}
	for k in (1, 2):
		levels = []
		for recs in records:
			level = collections.Counter()
			for _, seq in recs:
				seq = seq.upper()
				n = len(seq)
				for p in range(n - k + 1):
					windows = min(p, n - 8) - max(0, p + k - 8) + 1
					level[seq[p:p + k]] += windows / float(9 - k)
			levels.append(level)
		for word in (''.join(w) for w in itertools.product('ACGT', repeat=k)):
			rc = hh.revcomp(word)
			shares[word] = tuple(100.0 * (lv[word] + lv[rc]) / 2.0 / sum(lv.values())
				for lv in levels)
	return shares


@pytest.mark.parametrize('key', ['olen', 'omax1'])
def test_oout_table(runs, word_shares, key):
	"""-oout lists each 1-mer and 2-mer with its share of the targets' and the
	background's 8-mer windows and the factor target share / background
	share. The 2-mer background shares are taken after the 1-mer reweighting,
	so only their target shares are compared."""

	cwd, _ = runs[key]
	header, rows = hh.read_tsv(str(cwd / 'norm.txt'), header=True)
	assert header == ['Oligo', 'TargetCounts', 'BackgroundCounts', 'NormalizationFactor']
	got = {}
	for row in rows:
		assert row[:2] == ['', '']
		got[row[2]] = (float(row[3].rstrip('%')), float(row[4].rstrip('%')), float(row[5]))
	assert sorted(got) == sorted(word_shares)
	for word, (target, background) in word_shares.items():
		assert got[word][0] == pytest.approx(target, abs=0.0051), word
		if len(word) == 1:
			assert got[word][1] == pytest.approx(background, abs=0.0051), word
			assert got[word][2] == pytest.approx(target / background, abs=6e-4), word

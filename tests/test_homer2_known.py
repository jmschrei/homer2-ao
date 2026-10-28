"""homer2 known: enrichment of known motifs in target over background
sequences.

Expected values come from a Python reference. A sequence counts for a motif
when its best window (homer_helpers.best_score(), the same scoring
homer_helpers.scan() uses) reaches the motif's threshold, so that it holds at
least one site; the target and background counts are the sums of the
group-file weights of those sequences. The log p-value is the one
EnrichmentScorer computes: for the default binomial statistic,
ln P(X >= t) for X ~ Binomial(n, b / N), where t and b are the target and
background counts and n and N the target and background totals, each
truncated to an integer, with b / N kept within [1/N, (N-1)/N]; for
-stat hypergeo, ln P(X >= t) for t of t + b draws from n + N sequences of
which n are targets. -opt is checked against a Python walk over the
sequences' best scores, and against a scan at the threshold it writes.
"""

import math
import os

import pytest

import homer_helpers as hh

LOGP_TOLERANCE = 2e-3


# --- reading the output ---------------------------------------------------------------

def parse_table(text):
	"""The enrichment table as (totals, rows, names in output order): totals
	is (targets, background) from the header, rows maps a motif name to its
	columns."""

	lines = [line.split('\t') for line in text.splitlines() if line]
	header = lines[0]
	assert header[:5] == ['Motif Name', 'Consensus', 'p-value', 'Log p-value',
		'Log2 enrichment ratio'], header
	totals = []
	for column in (header[5], header[7]):
		assert '(total=' in column, column
		totals.append(int(column[column.index('(total=') + 7:column.rindex(')')]))
	rows, names = {}, []
	for cols in lines[1:]:
		assert len(cols) == 9, cols
		assert cols[0] not in rows
		names.append(cols[0])
		rows[cols[0]] = {
			'consensus': cols[1], 'pvalue': cols[2], 'logp': float(cols[3]),
			'ratio': float(cols[4]), 'targets': int(cols[5]),
			'target_percent': float(cols[6].rstrip('%')), 'background': int(cols[7]),
			'background_percent': float(cols[8].rstrip('%')),
		}
	return tuple(totals), rows, names


def parse_stats(field):
	"""'T:120.0(60.00%),B:0.0(0.00%),P:1e-302' or 'Tpos:-0.5,...,Multiplicity:1.00'
	as a dict of floats (the P: value as text)."""

	stats = {}
	for part in field.split(','):
		key, value = part.split(':', 1)
		if '(' in value:
			count, percent = value.rstrip('%)').split('(')
			stats[key] = float(count)
			stats[key + '%'] = float(percent)
		elif key == 'P':
			stats[key] = value
		else:
			stats[key] = float(value)
	return stats


def parse_mout(path):
	"""-mout's motifs by name: the Motif, the log p-value column and the
	statistics after it."""

	motifs = {}
	for motif in hh.read_motifs(str(path)):
		assert motif.header[1] == '0', motif.header
		entry = {'motif': motif, 'logp': float(motif.header[0])}
		entry.update(parse_stats(motif.header[2]))
		if len(motif.header) > 3:
			entry.update(parse_stats(motif.header[3]))
		motifs[motif.name] = entry
	return motifs


# --- the reference ------------------------------------------------------------------------

def log_pvalue(t, b, n_targets, n_background, stat='binomial'):
	t, b = int(t), int(b)
	n, N = int(n_targets), int(n_background)
	if stat == 'hypergeo':
		return hh.log_hypergeometric_tail(t, t + b, n, n + N)
	ratio = min(max(b / N, 1.0 / N), (N - 1.0) / N)
	return hh.log_binomial_tail(t, n, ratio)


def log2_ratio(t, b, n_targets, n_background):
	tfreq = max(t / n_targets, 1.0 / n_targets)
	bfreq = max(b / n_background, 1.0 / n_background)
	return math.log(tfreq / bfreq, 2)


class Reference(object):
	"""Best window scores of each motif in each sequence, kept per strand,
	and the counts they give."""

	def __init__(self, records, motifs):
		self.records = records
		self.motifs = {m.name: m for m in motifs}
		self._best = {}

	def best(self, name, strand='both'):
		if name not in self._best:
			motif = self.motifs[name]
			self._best[name] = [(hh.best_score(seq, motif, '+'), hh.best_score(seq, motif, '-'))
				for _, seq in self.records]
		scores = []
		for plus, minus in self._best[name]:
			pick = {'+': [plus], '-': [minus], 'both': [plus, minus]}[strand]
			pick = [s for s in pick if s is not None]
			scores.append(max(pick) if pick else None)
		return scores

	def counts(self, name, groups, strand='both', threshold=None):
		"""(target, background) weighted counts; groups maps a sequence name
		to (group, weight) and leaves out sequences the group file omits."""

		if threshold is None:
			threshold = self.motifs[name].threshold
		t = b = 0.0
		for (seq_name, _), score in zip(self.records, self.best(name, strand)):
			if seq_name not in groups or score is None or score < threshold:
				continue
			group, weight = groups[seq_name]
			if group == 1:
				t += weight
			else:
				b += weight
		return t, b


def totals(groups):
	n = sum(w for g, w in groups.values() if g == 1)
	N = sum(w for g, w in groups.values() if g == 0)
	return n, N


def assert_rows_match_reference(table, reference, groups, strand='both', stat='binomial'):
	(total_t, total_b), rows, names = table
	n, N = totals(groups)
	assert (total_t, total_b) == (int(n), int(N))
	assert sorted(names) == sorted(reference.motifs)
	logps = [rows[name]['logp'] for name in names]
	assert logps == sorted(logps), 'rows are sorted by log p-value'
	for name, row in rows.items():
		t, b = reference.counts(name, groups, strand)
		logp = log_pvalue(t, b, n, N, stat)
		assert row['consensus'] == reference.motifs[name].consensus
		assert (row['targets'], row['background']) == (int(t), int(b)), name
		assert row['target_percent'] == pytest.approx(100.0 * t / n, abs=0.0051), name
		assert row['background_percent'] == pytest.approx(100.0 * b / N, abs=0.0051), name
		assert row['logp'] == pytest.approx(logp, abs=LOGP_TOLERANCE), name
		assert row['pvalue'] == '1e{}'.format(int(row['logp'] / 2.3026)), name
		assert row['ratio'] == pytest.approx(log2_ratio(t, b, n, N), abs=0.0011), name


def optimize_threshold(scores, groups_in_order, n_targets, n_background, max_back=0.5):
	"""The walk -opt makes: sequences by best score, highest first; at each
	distinct score (scores within 1e-6 of the next one are merged into it) the
	sequences at or above it give a log p-value; the lowest wins, and the
	walk stops after the first cut holding more than max_back of the
	background. Returns (threshold, targets, background, logp)."""

	order = sorted(((-1e10 if s is None else s), g, w)
		for s, (g, w) in zip(scores, groups_in_order))[::-1]
	best = (None, 0.0, 0.0, 1.0)
	t = b = 0.0
	limit = int(n_background) * max_back
	for i in range(len(order) - 1):
		score, group, weight = order[i]
		if group == 1:
			t += weight
		else:
			b += weight
		if score < order[i + 1][0] + 1e-6:
			continue
		logp = log_pvalue(t, b, n_targets, n_background)
		if logp < best[3]:
			best = (score - 1e-6, t, b, logp)
		if b > limit:
			break
	return best


def site_statistics(records, groups, motif, threshold=None):
	"""The Tpos/Tstd/Bpos/Bstd/StrandBias/Multiplicity statistics of the
	sites at `threshold` (default: the motif's own): site centres relative
	to the middle of the sequence, their weighted mean and spread in targets
	and in background, log2 of + over - target sites, and target sites per
	target sequence with a site, counting sites whose centres lie within half
	a motif of a later site once."""

	if threshold is not None:
		motif = hh.Motif(motif.consensus, motif.name, threshold, motif.rows)
	half = len(motif) // 2
	sums = {1: [0.0, 0.0], 0: [0.0, 0.0]}
	centres = {1: [], 0: []}
	plus = minus = 0.0
	unique = with_sites = 0.0
	for name, seq in records:
		group, weight = groups[name]
		sites = hh.scan(seq, motif)
		if not sites:
			continue
		cs = [pos + half if strand == '+' else pos - half for pos, _, strand, _ in sites]
		for j, (c, site) in enumerate(zip(cs, sites)):
			sums[group][0] += weight * c
			sums[group][1] += weight
			centres[group].append((c, weight))
			if group == 1:
				if site[2] == '+':
					plus += weight
				else:
					minus += weight
				if not any(-half < c - c2 < half for c2 in cs[j + 1:]):
					unique += weight
		if group == 1:
			with_sites += weight
	stats = {}
	for group, key in ((1, 'T'), (0, 'B')):
		total = sums[group][1]
		mean = sums[group][0] / total if total else 0.0
		var = sum((w * (c - mean)) ** 2 for c, w in centres[group]) / total if total else 0.0
		stats[key + 'pos'] = mean
		stats[key + 'std'] = math.sqrt(var)
	if minus < 1e-5 and plus > 1e-5:
		stats['StrandBias'] = 10.0
	elif minus > 1e-5 and plus < 1e-5:
		stats['StrandBias'] = -10.0
	elif minus < 1e-5 and plus < 1e-5:
		stats['StrandBias'] = 0.0
	else:
		stats['StrandBias'] = math.log(plus / minus, 2)
	stats['Multiplicity'] = unique / with_sites if with_sites > 0.01 else 1.0
	return stats


IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'AC': 'M', 'AG': 'R', 'AT': 'W',
	'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B',
	'ACGT': 'N'}


def iupac_consensus(rows):
	"""One IUPAC letter per position for the bases whose probability is
	within 0.15 of the position's highest."""

	letters = []
	for row in rows:
		top = max(row)
		letters.append(IUPAC[''.join(base for base, p in zip('ACGT', row) if top - p < 0.15)])
	return ''.join(letters)


# --- fixtures -------------------------------------------------------------------------------

LIBRARY = os.path.join(hh.HOMER_HOME, 'data', 'knownTFs', 'vertebrates', 'known.motifs')


@pytest.fixture(scope='module')
def inputs(tmp_path_factory):
	"""Targets and background as one sequence file, two group files (unit
	weights, and other weights with 100 background sequences left out), and two
	motif files: the five test motifs with REST-NRSF from HOMER's library,
	and 40 library motifs."""

	d = tmp_path_factory.mktemp('known_inputs')
	targets = hh.read_fasta(hh.TARGETS)
	background = hh.read_fasta(hh.BACKGROUND)
	records = targets + background
	unit = {name: (1, 1.0) for name, _ in targets}
	unit.update({name: (0, 1.0) for name, _ in background})
	weighted = {name: (1, 2.0 if i % 2 == 0 else 0.5) for i, (name, _) in enumerate(targets)}
	weighted.update({name: (0, 1.5 if i < 500 else 0.5)
		for i, (name, _) in enumerate(background[:900])})
	library = hh.read_motifs(LIBRARY)
	rest = [m for m in library if m.name.startswith('REST-NRSF')]
	assert len(rest) == 1
	motifs = hh.read_motifs(os.path.join(hh.MOTIFS, 'known_small.motifs')) + rest
	paths = {
		'seqs': hh.write_seq_tsv(str(d / 'seqs.tsv'), records),
		'groups': hh.write_group(str(d / 'groups.txt'),
			[(name, unit[name][0], unit[name][1]) for name, _ in records]),
		'weighted': hh.write_group(str(d / 'weighted.txt'),
			[(name,) + weighted[name] for name, _ in records if name in weighted]),
		'motifs': hh.write_motifs(str(d / 'test.motifs'), motifs),
		'library': hh.write_motifs(str(d / 'library.motifs'), library[:40]),
	}
	return {'dir': d, 'paths': paths, 'records': records, 'unit': unit,
		'weighted': weighted, 'reference': Reference(records, motifs)}


def known_args(paths, *extra, group='groups', motifs='motifs'):
	return ['homer2', 'known', '-s', paths['seqs'], '-g', paths[group], '-m',
		paths[motifs]] + list(extra)


@pytest.fixture(scope='module')
def runs(inputs):
	"""Every homer2 known call the module checks, run once. Calls with -p 4
	run one at a time."""

	d, p = inputs['dir'], inputs['paths']

	def mout(key):
		return str(d / key / 'out.motifs')

	calls = {
		'plain': known_args(p, '-mout', mout('plain')),
		'fasta': ['homer2', 'known', '-i', hh.TARGETS, '-b', hh.BACKGROUND, '-m', p['motifs']],
		'outfile': known_args(p, '-o', str(d / 'outfile' / 'table.txt')),
		'binomial': known_args(p, '-stat', 'binomial'),
		'hypergeo': known_args(p, '-stat', 'hypergeo'),
		'cache0': known_args(p, '-cache', 0),
		'cache1': known_args(p, '-cache', 1),
		'cache0_hypergeo': known_args(p, '-stat', 'hypergeo', '-cache', 0),
		'plus': known_args(p, '-strand', '+'),
		'minus': known_args(p, '-strand', '-'),
		'both': known_args(p, '-strand', 'both'),
		'weighted': known_args(p, '-mout', mout('weighted'), group='weighted'),
		'opt': known_args(p, '-opt', '-mout', mout('opt')),
		'opt_maxback': known_args(p, '-opt', '-maxBack', 0.01, '-mout', mout('opt_maxback')),
		'nlen': known_args(p, '-nlen', 2, '-nout', str(d / 'nlen' / 'nout.txt'),
			'-mout', mout('nlen')),
		'norm': ['homer2', 'norm', '-s', p['seqs'], '-g', p['groups'], '-nlen', 2,
			'-nout', str(d / 'norm' / 'nout.txt'), '-o', str(d / 'norm' / 'groups.txt')],
		'library_p1': known_args(p, '-p', 1, motifs='library'),
	}
	threaded = {
		'library_p4': known_args(p, '-p', 4, motifs='library'),
	}
	results = hh.run_many({key: (args, d / key) for key, args in calls.items()})
	results.update(hh.run_many({key: (args, d / key) for key, args in threaded.items()}, jobs=1))
	results['mout'] = {key: mout(key) for key in list(calls) + list(threaded)}
	return results


def table(runs, key):
	return parse_table(runs[key].stdout)


# --- tests -------------------------------------------------------------------------------------

def test_counts_and_pvalues_match_the_reference(inputs, runs):
	assert_rows_match_reference(table(runs, 'plain'), inputs['reference'], inputs['unit'])


def test_planted_sites_are_enriched(runs):
	"""siteA, siteB and REST-NRSF are planted in 120, 80 and 60 of the 200
	targets and nowhere in the background; they lead the table."""

	_, rows, names = table(runs, 'plain')
	assert set(names[:3]) == {'siteA', 'siteB', 'REST-NRSF(Zf)/Jurkat-NRSF-ChIP-Seq/Homer'}
	assert rows['siteA']['targets'] == 120
	assert rows['siteB']['targets'] >= 80
	assert rows['REST-NRSF(Zf)/Jurkat-NRSF-ChIP-Seq/Homer']['targets'] >= 60
	assert all(rows[name]['logp'] < -100 for name in names[:3])


def test_fasta_input_matches_sequence_and_group_files(inputs, runs):
	"""-i and -b give the same table as the same sequences in a sequence and
	a group file, and leave no temporary files behind."""

	assert runs['fasta'].stdout == runs['plain'].stdout
	assert os.listdir(str(inputs['dir'] / 'fasta')) == []


def test_output_file(inputs, runs):
	assert runs['outfile'].stdout == ''
	with open(str(inputs['dir'] / 'outfile' / 'table.txt')) as handle:
		assert handle.read() == runs['plain'].stdout


def test_binomial_is_the_default(runs):
	assert runs['binomial'].stdout == runs['plain'].stdout


def test_hypergeometric(inputs, runs):
	assert_rows_match_reference(table(runs, 'hypergeo'), inputs['reference'], inputs['unit'],
		stat='hypergeo')
	assert table(runs, 'hypergeo')[1]['siteB']['logp'] != table(runs, 'plain')[1]['siteB']['logp']


@pytest.mark.parametrize('key,same_as', [('cache0', 'plain'), ('cache1', 'plain'),
	('cache0_hypergeo', 'hypergeo')])
def test_cache_size_does_not_change_results(runs, key, same_as):
	"""-cache only sets how many p-values are kept for reuse."""

	assert runs[key].stdout == runs[same_as].stdout


@pytest.mark.parametrize('strand', ['+', '-'])
def test_strand(inputs, runs, strand):
	key = 'plus' if strand == '+' else 'minus'
	assert_rows_match_reference(table(runs, key), inputs['reference'], inputs['unit'], strand)
	plain = table(runs, 'plain')[1]
	assert table(runs, key)[1]['siteA']['targets'] < plain['siteA']['targets']


def test_both_strands_is_the_default(runs):
	assert runs['both'].stdout == runs['plain'].stdout


def test_group_file_weights(inputs, runs):
	"""Counts are sums of the weights in the group file, totals the sums of
	all target and all background weights, and sequences the group file
	leaves out are not counted."""

	groups = inputs['weighted']
	assert totals(groups) == (250.0, 950.0)
	reference = inputs['reference']
	assert_rows_match_reference(table(runs, 'weighted'), reference, groups)
	mout = parse_mout(runs['mout']['weighted'])
	for name, entry in mout.items():
		t, b = reference.counts(name, groups)
		assert (entry['T'], entry['B']) == pytest.approx((t, b), abs=0.051), name
	assert table(runs, 'weighted')[1]['siteA']['targets'] != 120


def test_mout_holds_each_motif_with_its_statistics(inputs, runs):
	"""-mout writes every motif with its threshold, the log p-value, and
	"T:<targets>(<%>),B:<background>(<%>),P:1e<exponent>" in the header; the
	matrix is the input's (rows renormalized, three decimals), and the
	consensus is made again from it (see iupac_consensus)."""

	_, rows, names = table(runs, 'plain')
	assert iupac_consensus(inputs['reference'].motifs['siteA'].rows) == 'TAGCGTACCGTA'
	written = hh.read_motifs(runs['mout']['plain'])
	assert [m.name for m in written] == names
	mout = parse_mout(runs['mout']['plain'])
	for name, original in inputs['reference'].motifs.items():
		entry, row = mout[name], rows[name]
		motif = entry['motif']
		assert motif.consensus == iupac_consensus(motif.rows)
		assert motif.threshold == pytest.approx(original.threshold, abs=1e-6)
		assert entry['logp'] == pytest.approx(row['logp'], abs=6e-4)
		assert entry['P'] == row['pvalue']
		t, b = inputs['reference'].counts(name, inputs['unit'])
		assert (entry['T'], entry['B']) == (t, b)
		assert entry['T%'] == pytest.approx(100.0 * t / 200, abs=0.0051)
		assert entry['B%'] == pytest.approx(100.0 * b / 1000, abs=0.0051)
		assert len(motif.rows) == len(original.rows)
		for got, want in zip(motif.rows, original.rows):
			total = sum(max(p, hh.PROBABILITY_MINIMUM) for p in want)
			assert got == pytest.approx([max(p, hh.PROBABILITY_MINIMUM) / total for p in want],
				abs=6e-4)


@pytest.mark.parametrize('key', ['opt', 'opt_maxback'])
def test_opt_counts_match_a_scan_at_the_written_threshold(inputs, runs, key):
	"""With -opt, the counts and p-value reported for each motif are those of
	a scan at the threshold -mout writes for it."""

	reference = inputs['reference']
	_, rows, _ = table(runs, key)
	mout = parse_mout(runs['mout'][key])
	for name, row in rows.items():
		threshold = mout[name]['motif'].threshold
		t, b = reference.counts(name, inputs['unit'], threshold=threshold)
		assert (row['targets'], row['background']) == (int(t), int(b)), name
		assert (mout[name]['T'], mout[name]['B']) == (t, b), name
		assert row['logp'] == pytest.approx(log_pvalue(t, b, 200, 1000), abs=LOGP_TOLERANCE)


@pytest.mark.parametrize('key,max_back', [('opt', 0.5), ('opt_maxback', 0.01)])
def test_opt_picks_the_most_enriched_threshold(inputs, runs, key, max_back):
	"""-opt moves each threshold to the best-scoring cut among the
	sequences' best window scores, searching only until more than -maxBack
	of the background is above the cut (see optimize_threshold)."""

	reference = inputs['reference']
	groups = [inputs['unit'][name] for name, _ in inputs['records']]
	_, rows, _ = table(runs, key)
	mout = parse_mout(runs['mout'][key])
	for name, row in rows.items():
		threshold, t, b, logp = optimize_threshold(reference.best(name), groups, 200, 1000,
			max_back)
		assert mout[name]['motif'].threshold == pytest.approx(threshold, abs=2e-6), name
		assert (mout[name]['T'], mout[name]['B']) == (t, b), name
		assert row['logp'] == pytest.approx(logp, abs=LOGP_TOLERANCE), name


def test_opt_improves_a_motif_with_a_loose_threshold(runs):
	"""siteB's own threshold admits one mismatch, which chance matches in the
	background reach; -opt raises it to exact matches only."""

	plain, opt = table(runs, 'plain')[1]['siteB'], table(runs, 'opt')[1]['siteB']
	assert plain['background'] > 0
	assert opt['background'] == 0
	assert opt['logp'] < plain['logp']


def test_max_back_changes_the_search(runs):
	opt = parse_mout(runs['mout']['opt'])
	limited = parse_mout(runs['mout']['opt_maxback'])
	assert any(opt[name]['motif'].threshold != limited[name]['motif'].threshold for name in opt)


def statistics_mismatches(entry, want):
	keys = [key for key in ('Tpos', 'Tstd', 'Bpos', 'Bstd', 'StrandBias')
		if abs(entry[key] - want[key]) > 0.051]
	if abs(entry['Multiplicity'] - want['Multiplicity']) > 0.0051:
		keys.append('Multiplicity')
	return keys


@pytest.mark.parametrize('name', ['siteA', 'siteB', 'ebox', 'decoy1'])
def test_opt_mout_site_statistics(inputs, runs, name):
	"""With -opt, -mout adds where the motif's sites lie (see
	site_statistics). HOMER computes these statistics from the sites at the
	motif's own threshold, before it optimizes the threshold, while the
	counts beside them are those at the optimized threshold, which looks
	like a defect; the test accepts statistics of the sites at either
	threshold. For siteA, siteB and decoy1 the two differ."""

	entry = parse_mout(runs['mout']['opt'])[name]
	motif = inputs['reference'].motifs[name]
	versions = {
		'own threshold': site_statistics(inputs['records'], inputs['unit'], motif),
		'optimized threshold': site_statistics(inputs['records'], inputs['unit'], motif,
			entry['motif'].threshold),
	}
	found = {label: statistics_mismatches(entry, want) for label, want in versions.items()}
	assert any(not keys for keys in found.values()), found


def test_threads_give_the_same_counts(runs):
	"""At -p 4 every motif has the counts and percentages it has at -p 1.
	Only these are compared: the threads share one cache of p-values, so the
	p-values, and the row order that follows them, are left out."""

	def counts(key):
		totals, rows, _ = table(runs, key)
		return totals, {name: (row['targets'], row['target_percent'], row['background'],
			row['background_percent']) for name, row in rows.items()}

	one, four = counts('library_p1'), counts('library_p4')
	assert len(one[1]) == 40
	assert four == one


def read_weights(path):
	return {row[0]: (int(row[1]), float(row[2])) for row in hh.read_tsv(str(path))}


def test_nlen_weights_background_as_homer2_norm_does(inputs, runs):
	"""-nlen reweights the background before counting: the counts are sums
	of the weights homer2 norm writes for the same option."""

	key = 'nlen'
	groups = read_weights(inputs['dir'] / 'norm' / 'groups.txt')
	assert groups != inputs['unit']
	reference = inputs['reference']
	(total_t, total_b), rows, _ = table(runs, key)
	n, N = totals(groups)
	assert (total_t, total_b) == (int(n), int(N))
	mout = parse_mout(runs['mout'][key])
	for name, row in rows.items():
		t, b = reference.counts(name, groups)
		assert mout[name]['T'] == t
		assert mout[name]['B'] == pytest.approx(b, abs=0.051)
		assert row['background_percent'] == pytest.approx(100.0 * b / N, abs=0.0051)
		assert row['logp'] == pytest.approx(log_pvalue(t, b, n, N), abs=LOGP_TOLERANCE)


def test_nout_writes_the_normalization_table(inputs, runs):
	"""-nout writes the same table homer2 norm writes for the same input."""

	with open(str(inputs['dir'] / 'nlen' / 'nout.txt')) as known, \
			open(str(inputs['dir'] / 'norm' / 'nout.txt')) as norm:
		text = known.read()
		assert text == norm.read()
	lines = text.splitlines()
	assert lines[0] == 'Oligo\tTargetCounts\tBackgroundCounts\tNormalizationFactor'
	assert len(lines) == 1 + 4 + 16

"""homer2 norm: reweight background sequences so that their short k-mer
frequencies match the targets', and write the new group file.

The expectations come from normalize() below, a Python reference of the
normalization as SequenceArray::normalizeSequence in cpp/Motif2.cpp does it.
For each k from 1 to -nlen, every sequence spreads its weight evenly over its
k-mer windows that hold no N; with both strands a k-mer's level is the mean
of its own and its reverse complement's. A k-mer's factor is its target
fraction over its background fraction (each fraction at least 0.5 / total),
and each background sequence's weight moves toward its weight times the mean
factor of its windows, in a damped step. With -neutral, targets and
background are both moved toward equal frequencies (4^-k per k-mer) instead.
Rounds repeat while the summed squared deviation of the factors from 1
falls, up to the -nmax limit.

Two details of HOMER's code look like defects, so no test depends on them.
It leaves out each sequence's last k-mer window, and after up to -nmax rounds
it runs one more, so -nmax N allows N + 1 rounds where the usage promises at
most N. The reference runs each way (with and without the last window; when
the rounds reach -nmax, -nmax rounds and a final one, or -nmax rounds in
all), and a test that compares weights or the -nout table with it passes
when HOMER matches any of these versions.

Small inputs are checked against the reference weight by weight. Larger
inputs are checked for what can be derived without running the whole
reference: target weights stay, background k-mer frequencies move toward the
targets', and the -nout table, which describes the first round.
"""

import os

import pytest

import homer_helpers as hh

INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
MINIMUM_WEIGHT = 0.001


# --- the reference ----------------------------------------------------------------------

def kmer_codes(seq, k, last_window=True):
	"""The code of each k-mer window, first base least significant, or None
	for a window holding N; without `last_window`, the window that starts at
	the last k-th base is left out."""

	codes = []
	for j in range(len(seq) - k + (1 if last_window else 0)):
		window = seq[j:j + k]
		if 'N' in window:
			codes.append(None)
		else:
			codes.append(sum(INDEX[base] * 4 ** pos for pos, base in enumerate(window)))
	return codes


def kmer_text(code, k):
	return ''.join('ACGT'[code // 4 ** j % 4] for j in range(k))


def revcomp_code(code, k):
	return sum((3 - (code // 4 ** j % 4)) * 4 ** (k - 1 - j) for j in range(k))


def kmer_levels(seqs, k, strand='both'):
	"""(target levels, background levels, target total, background total),
	levels indexed by k-mer code: each sequence adds its weight divided by
	its number of counted windows to the level of each window's k-mer."""

	total = 4 ** k
	levels = {1: [0.0] * total, 0: [0.0] * total}
	for s in seqs:
		codes = s['codes'][k]
		s['good'] = sum(1 for c in codes if c is not None)
		for c in codes:
			if c is not None:
				levels[1 if s['group'] == 1 else 0][c] += s['weight'] / s['good']
	sums = {}
	for group in (1, 0):
		sums[group] = 0.0
		for value in levels[group]:
			sums[group] += value
	if strand == 'both':
		for group in (1, 0):
			raw = levels[group]
			levels[group] = [(raw[i] + raw[revcomp_code(i, k)]) / 2.0 for i in range(total)]
	return levels[1], levels[0], sums[1], sums[0]


def normalization_round(seqs, nlen, strand, neutral, table=None):
	"""One round over k = 1 .. nlen; updates the weights in place, appends
	the -nout rows to `table` when given, and returns the round's error."""

	error = 0.0
	for k in range(1, nlen + 1):
		total = 4 ** k
		target, background, sum_t, sum_b = kmer_levels(seqs, k, strand)
		factors = {1: [1.0] * total, 0: [1.0] * total}
		for i in range(total):
			t = max(target[i], 0.5 / sum_t)
			b = max(background[i], 0.5 / sum_b)
			if not neutral:
				factor = sum_b * t / sum_t / b
				factors[0][i] = factor
				error += (factor - 1.0) * (factor - 1.0) / total
				row = (factor,)
			else:
				factors[1][i] = sum_t / t * 4.0 ** -k
				factors[0][i] = sum_b / b * 4.0 ** -k
				error += (factors[1][i] - 1.0) * (factors[1][i] - 1.0) / total
				error += (factors[0][i] - 1.0) * (factors[0][i] - 1.0) / total
				row = (factors[1][i], factors[0][i])
			if table is not None:
				table.append((kmer_text(i, k), 100.0 * t / sum_t, 100.0 * b / sum_b) + row)
		for s in seqs:
			if s['group'] == 1 and not neutral:
				continue
			factor = factors[1 if s['group'] == 1 else 0]
			score = 0.0
			for c in s['codes'][k]:
				if c is not None:
					score += factor[c]
			if s['good'] > 0.5:
				score /= s['good']
			new = min(max(score * s['weight'], MINIMUM_WEIGHT), 1.0 / MINIMUM_WEIGHT)
			penalty = (new if new >= 1 else 1.0 / new) ** 2
			delta = new - s['weight']
			if penalty > 1 and ((delta > 0 and new > 1) or (delta < 0 and new < 1)):
				s['weight'] += delta / penalty
			else:
				s['weight'] += delta
	return error


def prepare(rows, nlen, last_window):
	# HOMER holds the sequences last read first; summing in that order keeps
	# the rounding the same.
	return [{'name': name, 'group': group, 'weight': weight,
		'codes': {k: kmer_codes(hh.clean_sequence(seq), k, last_window)
			for k in range(1, nlen + 1)}}
		for name, seq, group, weight in reversed(rows)]


def normalize(rows, nlen, strand='both', nmax=160, neutral=False, last_window=True,
		extra_round=True):
	"""rows: (name, sequence, group, weight). Rounds run while the error
	falls, the round that stops it included, up to `nmax` of them, and then
	one final round; without `extra_round` the final round counts toward
	`nmax`. Returns ({name: (group, weight)}, the number of rounds before the
	final one)."""

	seqs = prepare(rows, nlen, last_window)
	rounds = 0
	if nlen > 0:
		last = 1e100
		for _ in range(nmax if extra_round else nmax - 1):
			error = normalization_round(seqs, nlen, strand, neutral)
			rounds += 1
			if error >= last:
				break
			last = error
		if extra_round or nmax > 0:
			normalization_round(seqs, nlen, strand, neutral)
	return {s['name']: (s['group'], s['weight']) for s in seqs}, rounds


def first_round_table(rows, nlen, strand='both', neutral=False, last_window=True):
	table = []
	normalization_round(prepare(rows, nlen, last_window), nlen, strand, neutral, table)
	return table


def weight_versions(rows, **options):
	"""The reference's weights with and without each sequence's last k-mer
	window and, when the rounds reach -nmax, with -nmax rounds and a final
	one or -nmax rounds in all."""

	versions = {}
	for last_window in (True, False):
		label = 'with last window' if last_window else 'without last window'
		weights, rounds = normalize(rows, last_window=last_window, **options)
		versions[label + ', -nmax rounds and a final one'] = weights
		if rounds == options.get('nmax', 160):
			versions[label + ', -nmax rounds in all'] = normalize(rows, last_window=last_window,
				extra_round=False, **options)[0]
	return versions


def both_versions(function, *args, **kwargs):
	"""The reference's result counting each sequence's last k-mer window and
	not counting it."""

	return {'with last window': function(*args, last_window=True, **kwargs),
		'without last window': function(*args, last_window=False, **kwargs)}


def kmer_fractions(rows, weights, k):
	"""Weighted fraction of each k-mer (with its reverse complement) in the
	targets and in the background, every window of a sequence counting
	equally."""

	seqs = [{'group': weights[name][0], 'weight': weights[name][1],
		'codes': {k: kmer_codes(hh.clean_sequence(seq), k)}} for name, seq, _, _ in rows]
	target, background, sum_t, sum_b = kmer_levels(seqs, k)
	return [v / sum_t for v in target], [v / sum_b for v in background]


def distance(a, b):
	return sum(abs(x - y) for x, y in zip(a, b))


# --- reading the output ---------------------------------------------------------------------

def read_groups(text):
	groups = {}
	for line in text.splitlines():
		cols = line.split('\t')
		assert len(cols) == 3, line
		assert cols[0] not in groups
		groups[cols[0]] = (int(cols[1]), float(cols[2]))
	return groups


def read_nout(path):
	rows = hh.read_tsv(str(path))
	assert rows[0] == ['Oligo', 'TargetCounts', 'BackgroundCounts', 'NormalizationFactor']
	return [(r[0], float(r[1].rstrip('%')), float(r[2].rstrip('%')))
		+ tuple(float(f) for f in r[3].split(',')) for r in rows[1:]]


def assert_same_weights(got, want, tolerance=2e-6):
	assert sorted(got) == sorted(want)
	for name, (group, weight) in want.items():
		assert got[name][0] == group, name
		assert got[name][1] == pytest.approx(weight, abs=tolerance), name


def weight_mismatches(got, want, tolerance=2e-6):
	"""Names whose group or weight differs, or that only one side has."""

	names = sorted(set(got) ^ set(want))
	for name in sorted(set(got) & set(want)):
		if got[name][0] != want[name][0] or abs(got[name][1] - want[name][1]) > tolerance:
			names.append(name)
	return names


def table_mismatches(got, want):
	"""K-mers whose -nout row differs beyond the rounding of the output."""

	if [row[0] for row in got] != [row[0] for row in want]:
		return ['k-mer column']
	return [g[0] for g, w in zip(got, want)
		if any(abs(a - b) > 0.0051 for a, b in zip(g[1:3], w[1:3]))
		or len(g) != len(w) or any(abs(a - b) > 0.00051 for a, b in zip(g[3:], w[3:]))]


def assert_matches_a_version(got, versions, mismatches):
	found = {label: mismatches(got, want) for label, want in versions.items()}
	assert any(not differ for differ in found.values()), \
		{label: differ[:5] for label, differ in found.items()}


# --- fixtures --------------------------------------------------------------------------------

@pytest.fixture(scope='module')
def inputs(tmp_path_factory):
	"""`full`: the 200 targets and 1,000 background sequences with weight 1.
	`skewed`: the 200 targets and the 500 background sequences poorest in G
	and C, with weight 1. `small`: 40 targets and 80 background sequences,
	made-up background sequences (one with N, one all G and C, one all N, one
	of 3 bases), some weights other than 1, and one sequence the group file
	leaves out."""

	d = tmp_path_factory.mktemp('norm_inputs')
	targets = hh.read_fasta(hh.TARGETS)
	background = hh.read_fasta(hh.BACKGROUND)
	full = [(n, s, 1, 1.0) for n, s in targets] + [(n, s, 0, 1.0) for n, s in background]
	by_gc = sorted(background, key=lambda r: (sum(r[1].upper().count(b) for b in 'GC'), r[0]))
	skewed = [(n, s, 1, 1.0) for n, s in targets] + [(n, s, 0, 1.0) for n, s in by_gc[:500]]
	small = [(n, s, 1, 0.5 if i % 5 == 0 else 1.0) for i, (n, s) in enumerate(targets[:40])]
	small += [(n, s, 0, 2.0 if i % 7 == 0 else 1.0) for i, (n, s) in enumerate(background[:80])]
	small += [('with_n', 'ACGTNNACGTAGGCTANCGATCGGGCCCATAT' * 3, 0, 1.0),
		('gc_only', 'GCGCGGCCGCGCGGGCCC' * 5, 0, 1.0), ('all_n', 'N' * 50, 0, 1.0),
		('short', 'ACG', 0, 1.0)]
	files = {}
	for key, rows in (('full', full), ('skewed', skewed), ('small', small)):
		extra = [('not_in_group_file', 'ACGT' * 30)] if key == 'small' else []
		files[key] = {
			'seqs': hh.write_seq_tsv(str(d / (key + '.tsv')),
				[(n, s) for n, s, _, _ in rows] + extra),
			'groups': hh.write_group(str(d / (key + '.groups')),
				[(n, g, w) for n, _, g, w in rows]),
		}
	return {'dir': d, 'files': files, 'full': full, 'skewed': skewed, 'small': small}


SMALL_RUNS = {
	'nlen1': dict(nlen=1),
	'nlen2': dict(nlen=2),
	'nlen3': dict(nlen=3),
	'plus': dict(nlen=2, strand='+'),
	'minus': dict(nlen=2, strand='-'),
	'nmax1': dict(nlen=2, nmax=1),
	'nmax3': dict(nlen=2, nmax=3),
	'neutral': dict(nlen=2, neutral=True),
	'nlen5': dict(nlen=5, nmax=3),
	'nlen5_threads': dict(nlen=5, nmax=3, threads=4),
}


def norm_args(files, nlen=None, strand=None, nmax=None, neutral=False, threads=None):
	args = ['homer2', 'norm', '-s', files['seqs'], '-g', files['groups']]
	if threads is not None:
		args += ['-p', threads]
	if nlen is not None:
		args += ['-nlen', nlen]
	if strand is not None:
		args += ['-strand', strand]
	if nmax is not None:
		args += ['-nmax', nmax]
	if neutral:
		args.append('-neutral')
	return args


@pytest.fixture(scope='module')
def runs(inputs):
	"""Every homer2 norm call the module checks, run once. Calls with -p 4
	run one at a time."""

	d, files = inputs['dir'], inputs['files']
	full, small = files['full'], files['small']
	calls = {'small_' + key: norm_args(small, **options) for key, options in SMALL_RUNS.items()
		if 'threads' not in options}
	calls['small_neutral'] += ['-nout', str(d / 'small_neutral' / 'nout.txt')]
	calls.update({
		'small_default': norm_args(small),
		'small_nlen0': norm_args(small, nlen=0),
		'full_default': norm_args(full),
		'full_nlen2': norm_args(full, nlen=2),
		'skewed_nlen2': norm_args(files['skewed'], nlen=2) + ['-nout',
			str(d / 'skewed_nlen2' / 'nout.txt')],
		'full_fasta': ['homer2', 'norm', '-i', hh.TARGETS, '-b', hh.BACKGROUND, '-nlen', 2],
		'full_outfile': norm_args(full, nlen=2) + ['-o', str(d / 'full_outfile' / 'groups.txt')],
		'full_both': norm_args(full, nlen=2, strand='both'),
		'full_nlen3_p1': norm_args(full, nlen=3, threads=1),
		'full_nlen5_p1': norm_args(full, nlen=5, threads=1),
	})
	threaded = {'small_' + key: norm_args(small, **options) for key, options in SMALL_RUNS.items()
		if 'threads' in options}
	threaded.update({
		'full_nlen3_p4': norm_args(full, nlen=3, threads=4),
		'full_nlen5_p4': norm_args(full, nlen=5, threads=4),
	})
	results = hh.run_many({key: (args, d / key) for key, args in calls.items()})
	results.update(hh.run_many({key: (args, d / key) for key, args in threaded.items()}, jobs=1))
	return results


def weights_of(runs, key):
	return read_groups(runs[key].stdout)


# --- tests -------------------------------------------------------------------------------------

@pytest.mark.parametrize('key', ['small_default', 'small_nlen0'])
def test_without_nlen_weights_are_unchanged(inputs, runs, key):
	"""-nlen defaults to 0, which normalizes nothing: the output is the group
	file read, for the sequences it names."""

	want = {name: (group, weight) for name, _, group, weight in inputs['small']}
	assert_same_weights(weights_of(runs, key), want, tolerance=0)


def test_full_inputs_without_nlen(inputs, runs):
	want = {name: (group, weight) for name, _, group, weight in inputs['full']}
	assert_same_weights(weights_of(runs, 'full_default'), want, tolerance=0)


@pytest.mark.parametrize('key', sorted(SMALL_RUNS))
def test_weights_match_the_reference(inputs, runs, key):
	"""The weights match one version of the reference (weight_versions):
	counting each sequence's last k-mer window or not and, when the rounds
	reach -nmax (they do in every run here, the default 160 included),
	running -nmax rounds and a final one or -nmax rounds in all."""

	options = {k: v for k, v in SMALL_RUNS[key].items() if k != 'threads'}
	got = weights_of(runs, 'small_' + key)
	assert_matches_a_version(got, weight_versions(inputs['small'], **options), weight_mismatches)
	assert 'not_in_group_file' not in got


@pytest.mark.parametrize('key', ['full', 'skewed'])
def test_target_weights_are_unchanged(inputs, runs, key):
	got = weights_of(runs, key + '_nlen2')
	targets = [name for name, _, group, _ in inputs[key] if group == 1]
	background = [name for name, _, group, _ in inputs[key] if group == 0]
	assert all(got[name] == (1, 1.0) for name in targets)
	assert all(got[name][0] == 0 for name in background)
	assert sum(1 for name in background if got[name][1] != 1.0) > 0.9 * len(background)


@pytest.mark.parametrize('k', [1, 2])
def test_background_kmer_frequencies_move_toward_the_targets(inputs, runs, k):
	"""With a background poorer in G and C than the targets, the
	background's weighted 1-mer and 2-mer frequencies (both strands) end
	closer to the targets' than they start: the summed absolute difference
	falls by more than a third."""

	rows = inputs['skewed']
	before = {name: (group, weight) for name, _, group, weight in rows}
	after = weights_of(runs, 'skewed_nlen2')
	target, background = kmer_fractions(rows, before, k)
	start = distance(target, background)
	target_after, background_after = kmer_fractions(rows, after, k)
	assert target_after == pytest.approx(target, abs=1e-12)
	assert start > 0.1
	assert distance(target, background_after) < start * 2 / 3


def test_nout_table(inputs, runs):
	"""-nout writes the first round's table: each k-mer up to -nlen, its
	target and background percentages, and its factor, as the reference
	computes them counting each sequence's last k-mer window or not."""

	got = read_nout(inputs['dir'] / 'skewed_nlen2' / 'nout.txt')
	assert len(got) == 4 + 16
	assert [row[0] for row in got[:6]] == ['A', 'C', 'G', 'T', 'AA', 'CA']
	assert_matches_a_version(got, both_versions(first_round_table, inputs['skewed'], 2),
		table_mismatches)


def test_nout_table_with_neutral(inputs, runs):
	"""With -neutral the factor column holds the target and the background
	factor, comma-separated; values as in test_nout_table."""

	got = read_nout(inputs['dir'] / 'small_neutral' / 'nout.txt')
	assert all(len(row) == 5 for row in got)
	assert_matches_a_version(got, both_versions(first_round_table, inputs['small'], 2,
		neutral=True), table_mismatches)


def test_neutral_moves_both_groups_toward_equal_frequencies(inputs, runs):
	rows = inputs['small']
	before = {name: (group, weight) for name, _, group, weight in rows}
	after = weights_of(runs, 'small_neutral')
	assert any(after[name][1] != weight for name, _, group, weight in rows if group == 1)
	equal = [0.25] * 4
	for start, end in zip(kmer_fractions(rows, before, 1), kmer_fractions(rows, after, 1)):
		assert distance(end, equal) < distance(start, equal)


def test_minus_strand_gives_the_plus_strand_weights(runs):
	"""Only -strand both pairs each k-mer with its reverse complement. The
	- strand holds the reverse complement of each + strand k-mer at the same
	frequency, so -strand - gives the weights of -strand +."""

	assert runs['small_minus'].stdout == runs['small_plus'].stdout
	assert weights_of(runs, 'small_plus') != weights_of(runs, 'small_nlen2')


def test_both_strands_is_the_default(runs):
	assert runs['full_both'].stdout == runs['full_nlen2'].stdout


def test_nmax_limits_the_rounds(runs):
	assert weights_of(runs, 'small_nmax1') != weights_of(runs, 'small_nmax3')
	assert weights_of(runs, 'small_nmax3') != weights_of(runs, 'small_nlen2')


def test_fasta_input_matches_sequence_and_group_files(inputs, runs):
	"""With -i and -b the output is still a group file, naming the FASTA
	records, and no temporary files are left behind."""

	assert_same_weights(weights_of(runs, 'full_fasta'), weights_of(runs, 'full_nlen2'),
		tolerance=0)
	assert os.listdir(str(inputs['dir'] / 'full_fasta')) == []


def test_output_file(inputs, runs):
	assert runs['full_outfile'].stdout == ''
	with open(str(inputs['dir'] / 'full_outfile' / 'groups.txt')) as handle:
		assert read_groups(handle.read()) == weights_of(runs, 'full_nlen2')


@pytest.mark.parametrize('nlen', [3, 5])
def test_threads_do_not_change_the_weights(runs, nlen):
	one, four = 'full_nlen{}_p1'.format(nlen), 'full_nlen{}_p4'.format(nlen)
	assert runs[four].stdout == runs[one].stdout
	assert weights_of(runs, one) != weights_of(runs, 'full_nlen2')

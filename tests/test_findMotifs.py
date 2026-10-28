"""findMotifs.pl in FASTA mode: known motif enrichment, de novo motif finding
and motif scanning on tests/data/targets.fa against tests/data/background.fa
(or against a scrambled copy of the targets), and the options that work
without a promoter set.

Expected values are derived independently:

- known motif counts come from homer_helpers.scan() over the targets and the
  background; with weights, a background count is the sum of the weights in
  the group file the run kept (-keepFiles), and p-values are the
  hypergeometric (default) or binomial (-b) tails of those counts;
- GC (or CpG) weights are recomputed from the sequences: the bins of
  freq2group.pl and the weights of assignGeneWeights.pl;
- options findMotifs.pl hands to homer2 (de novo and normalization options)
  are checked by running homer2 with those options on the files the run kept
  and comparing outputs, and by their effect on the output;
- -find output is compared with the reference scan, -dumpFasta, -min/-max and
  -chopify output with the input sequences cut by hand.

Every de novo run is small (-len 8, -S 1 or 2, one thread) and compares its
motifs against a small library (-mcheck) or not at all, except the -rna and
-mset runs, which use HOMER's RNA and worm libraries.
"""

import math
import os
import re
import shutil

import pytest

import homer_helpers as hh

KNOWN = os.path.join(hh.MOTIFS, 'known_small.motifs')
SITE_A = os.path.join(hh.MOTIFS, 'siteA.motif')
SITE_B = os.path.join(hh.MOTIFS, 'siteB.motif')
SITE_AB = os.path.join(hh.MOTIFS, 'siteA_siteB.motifs')
KNOWN_TFS = os.path.join(hh.HOMER_HOME, 'data', 'knownTFs')

BG = ['-fasta', hh.BACKGROUND]
DENOVO = ['-len', 8, '-S', 2, '-noknown', '-nocheck', '-keepFiles']
# -e as findMotifs.pl's usage gives its default
DENOVO_E = DENOVO + ['-e', 0.01]

RUNS = {
	'known': BG + ['-nomotif', '-mknown', KNOWN, '-keepFiles'],
	'known_fastaBg': ['-fastaBg', hh.BACKGROUND, '-nomotif', '-mknown', KNOWN],
	'known_p2': BG + ['-nomotif', '-mknown', KNOWN, '-p', 2],
	'known_cpu': BG + ['-nomotif', '-mknown', KNOWN, '-cpu', 2],
	'known_homer2': BG + ['-nomotif', '-mknown', KNOWN, '-homer2'],
	'noweight': BG + ['-nomotif', '-noweight', '-mknown', KNOWN, '-keepFiles'],
	'binomial': BG + ['-nomotif', '-noweight', '-b', '-mknown', KNOWN],
	'norevopp': BG + ['-nomotif', '-noweight', '-norevopp', '-mknown', KNOWN],
	'scrambled': ['-nomotif', '-noweight', '-mknown', KNOWN],
	'float': BG + ['-nomotif', '-float', '-mknown', KNOWN],
	'mask': BG + ['-nomotif', '-noweight', '-maskMotif', SITE_A, SITE_B, '-mknown', KNOWN,
		'-keepFiles'],
	'rand': BG + ['-nomotif', '-noweight', '-rand', '-mknown', KNOWN, '-keepFiles'],
	'nlen0': BG + ['-nomotif', '-noknown', '-nlen', 0, '-keepFiles'],
	'cpg': BG + ['-nomotif', '-noknown', '-nlen', 0, '-cpg', '-keepFiles'],
	'CpG': BG + ['-nomotif', '-noknown', '-nlen', 0, '-CpG', '-keepFiles'],
	'nlen2': BG + ['-nomotif', '-noknown', '-nlen', 2, '-keepFiles'],
	'nmax5': BG + ['-nomotif', '-noknown', '-nmax', 5, '-keepFiles'],
	'neutral': BG + ['-nomotif', '-noknown', '-neutral', '-keepFiles'],
	'find': ['-find', SITE_AB],
	'find_norevopp': ['-find', SITE_AB, '-norevopp'],
	'denovo': BG + ['-len', 8, '-S', 2, '-noknown', '-mcheck', KNOWN, '-nofacts', '-keepFiles'],
	'basic': BG + ['-len', 8, '-S', 2, '-noknown', '-mcheck', KNOWN, '-basic'],
	'bits': BG + ['-len', 8, '-S', 2, '-noknown', '-mcheck', KNOWN, '-bits'],
	'lens': BG + ['-len', '8,10', '-S', 1, '-noknown', '-mcheck', KNOWN],
	'denovo_e': BG + DENOVO_E,
	'mis0': BG + DENOVO_E + ['-mis', 0],
	'e': BG + DENOVO + ['-e', 0.0001],
	'quickMask': BG + DENOVO_E + ['-quickMask'],
	'minlp': BG + DENOVO_E + ['-minlp', -50],
	'cache': BG + DENOVO_E + ['-cache', 50],
	'denovo_norevopp': BG + DENOVO_E + ['-norevopp'],
	'denovo_binomial': BG + DENOVO_E + ['-b'],
	'fdr': BG + ['-len', 8, '-S', 2, '-noknown', '-nocheck', '-fdr', 2],
	'homer1': BG + ['-len', 8, '-S', 2, '-noknown', '-nocheck', '-homer1', '-keepFiles'],
	'mset': BG + ['-len', 8, '-S', 1, '-mset', 'worms'],
	'rna': BG + ['-len', 8, '-S', 1, '-rna', '-nofacts', '-keepFiles'],
	'rna_known': BG + ['-nomotif', '-rna', '-known', '-noweight'],
	'dump': ['-fasta', 'BACKGROUND_200', '-nomotif', '-noknown', '-dumpFasta'],
	'dump_cut': ['-fasta', 'BACKGROUND_200', '-nomotif', '-noknown', '-dumpFasta', '-min', 50,
		'-max', 250, '-chopify'],
}

# Runs with two threads, made one at a time.
THREADED = ('known_p2', 'known_cpu')

# The de novo option runs and the homer2 denovo options each one stands for.
DENOVO_OPTIONS = {
	'mis0': ['-mis', 0],
	'e': ['-e', 0.0001],
	'quickMask': ['-quickMask'],
	'minlp': ['-minlp', -50],
	'cache': ['-cache', 50],
	'denovo_norevopp': ['-strand', '+'],
	'denovo_binomial': ['-stat', 'binomial'],
}


def varied_targets():
	"""Targets of different lengths for -min/-max/-chopify: twenty of 100 bp,
	one each of 50 and 250 bp (kept: the limits are inclusive), and five each
	of 40 and 300 bp (removed)."""

	targets = hh.read_fasta(hh.TARGETS)
	records = [('t{:02d}'.format(i), seq[:100]) for i, (_, seq) in enumerate(targets[:20])]
	records.append(('t_min', targets[20][1][:50]))
	records.append(('t_max', targets[21][1] + targets[22][1][:50]))
	records += [('short{}'.format(i), seq[:40]) for i, (_, seq) in enumerate(targets[30:35])]
	records += [('long{}'.format(i), seq + targets[40 + i][1][:100])
		for i, (_, seq) in enumerate(targets[35:40])]
	return records


@pytest.fixture(scope='module')
def inputs(tmp_path_factory):
	d = tmp_path_factory.mktemp('fm_inputs')
	return {
		'varied': hh.write_fasta(str(d / 'varied.fa'), varied_targets()),
		# the first 200 background sequences, to keep -dumpFasta runs small
		'background_200': hh.write_fasta(str(d / 'background_200.fa'),
			hh.read_fasta(hh.BACKGROUND)[:200]),
	}


@pytest.fixture(scope='module')
def runs(tmp_path_factory, inputs):
	"""Every findMotifs.pl call of this module, run once: key -> the output
	directory (inside the directory the call ran in)."""

	base = tmp_path_factory.mktemp('fm_runs')
	calls = {}
	for key, options in RUNS.items():
		targets = inputs['varied'] if key == 'dump_cut' else hh.TARGETS
		options = [inputs['background_200'] if o == 'BACKGROUND_200' else o for o in options]
		calls[key] = (['findMotifs.pl', targets, 'fasta', 'out'] + options, base / key)
	results = hh.run_many({k: v for k, v in calls.items() if k not in THREADED})
	results.update(hh.run_many({k: calls[k] for k in THREADED}, jobs=1))
	for key, result in results.items():
		with open(str(base / key / 'stdout.txt'), 'w') as out:
			out.write(result.stdout)
	return {key: base / key / 'out' for key in results}


def homer2_denovo(out, options, cwd, length=8, motifs=2):
	"""homer2 denovo as findMotifs.pl runs it, on the sequences and weights
	the run kept, with `options` added."""

	return ['homer2', 'denovo', '-s', out / 'seq.tsv', '-g', out / 'group.adj', '-S', motifs,
		'-stat', 'hypergeo', '-mis', 1, '-cache', 500, '-minlp', -10, '-p', 1, '-len', length] + \
		options + ['-o', cwd / 'direct.motifs']


@pytest.fixture(scope='module')
def direct(tmp_path_factory, runs):
	"""homer2 run directly on the files findMotifs.pl kept, with the options
	findMotifs.pl should have passed, and a -reuse run on an edited copy of
	a kept output directory: key -> the directory the call ran in."""

	base = tmp_path_factory.mktemp('fm_direct')
	calls = {}
	for key, options in DENOVO_OPTIONS.items():
		cwd = base / key
		cwd.mkdir()
		if key != 'e':
			options = ['-e', 0.01] + options
		if key == 'denovo_binomial':
			args = homer2_denovo(runs[key], ['-e', 0.01], cwd)
			args[args.index('hypergeo')] = 'binomial'
		else:
			args = homer2_denovo(runs[key], options, cwd)
		calls[key] = (args, cwd)
	cwd = base / 'rna'
	cwd.mkdir()
	calls['rna'] = (homer2_denovo(runs['rna'], ['-strand', '+'], cwd, motifs=1), cwd)
	gc = runs['nlen0']
	for key, options in (('nlen3', ['-nlen', 3, '-nmax', 160]),
			('nlen2', ['-nlen', 2, '-nmax', 160]), ('nmax5', ['-nlen', 3, '-nmax', 5]),
			('neutral', ['-nlen', 3, '-nmax', 160, '-neutral'])):
		cwd = base / key
		cwd.mkdir()
		calls[key] = (['homer2', 'norm', '-g', gc / 'group.adj', '-s', gc / 'seq.tsv'] + options +
			['-nout', cwd / 'autonorm.tsv'], cwd)

	# -reuse: the kept files of the de novo run, every weight set to 1
	reuse = base / 'reuse'
	shutil.copytree(str(runs['denovo']), str(reuse / 'out'))
	rows = hh.read_tsv(str(reuse / 'out' / 'group.adj'))
	hh.write_group(str(reuse / 'out' / 'group.adj'), [(r[0], r[1], '1.000000') for r in rows])
	for name in os.listdir(str(reuse / 'out')):
		if name.startswith('homerMotifs') or name.startswith('homerResults') or \
				name == 'nonRedundant.motifs':
			path = reuse / 'out' / name
			if path.is_dir():
				shutil.rmtree(str(path))
			else:
				path.unlink()
	(base / 'reuse_direct').mkdir()
	calls['reuse_direct'] = (homer2_denovo(reuse / 'out', [], base / 'reuse_direct'),
		base / 'reuse_direct')
	results = hh.run_many(calls)
	for key in ('nlen3', 'nlen2', 'nmax5', 'neutral'):
		with open(str(base / key / 'group.adj'), 'w') as out:
			out.write(results[key].stdout)
	results['reuse'] = hh.run(['findMotifs.pl', hh.TARGETS, 'fasta', 'out'] + BG +
		['-len', 8, '-S', 2, '-noknown', '-nocheck', '-reuse', '-keepFiles'], reuse)
	return {key: base / key for key in results}


# --- reference values ---------------------------------------------------------------------

_SITES = {}


def has_sites(label, records, motif, strand='both'):
	"""{name: whether the reference scan finds a site}, cached per sequence
	set, motif and strand."""

	key = (label, motif.name, motif.threshold, tuple(map(tuple, motif.rows)), strand)
	if key not in _SITES:
		_SITES[key] = {name: hh.has_site(seq, motif, strand) for name, seq in records}
	return _SITES[key]


def read_group(path):
	"""{name: (group, weight)}; weight 1 when the file has no weight column."""

	return {row[0]: (int(row[1]), float(row[2]) if len(row) > 2 else 1.0)
		for row in hh.read_tsv(str(path))}


def read_seq_tsv(path):
	return [(row[0], row[1]) for row in hh.read_tsv(str(path))]


def known_by_name(out):
	return {row['name']: row for row in hh.read_known_results(str(out / 'knownResults.txt'))}


def counts(records, group, motif, label, strand='both'):
	"""(target weight with a site, background weight with a site) over the
	sequences in the group file."""

	found = has_sites(label, records, motif, strand)
	t = b = 0.0
	for name, _ in records:
		if name in group and found[name]:
			if group[name][0] == 1:
				t += group[name][1]
			else:
				b += group[name][1]
	return t, b


def input_records():
	"""The sequences findMotifs.pl analyses: the targets and the background,
	upper case, background names prefixed with "BG"."""

	records = [(n, hh.clean_sequence(s)) for n, s in hh.read_fasta(hh.TARGETS)]
	records += [('BG' + n, hh.clean_sequence(s)) for n, s in hh.read_fasta(hh.BACKGROUND)]
	return records


def benjamini(pvalues):
	"""Statistics::benjaminiFDR."""

	order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
	fdr = [1.0] * len(pvalues)
	last, last_index = -1e10, -1
	for i, k in enumerate(order):
		if i < len(order) - 1 and pvalues[order[i + 1]] == pvalues[k]:
			continue
		value = min(1.0, max(last, pvalues[k] * len(pvalues) / (i + 1)))
		for j in range(i, last_index, -1):
			fdr[order[j]] = value
		last, last_index = value, i
	return fdr


GC_BINS = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8]
CPG_BINS = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16]


def gc_weights(records, cpg=False):
	"""The group file assignGeneWeights.pl writes: sequences binned by GC
	fraction (CpG dinucleotide fraction with -cpg) as homerTools freq prints
	it (six decimals), bins without targets or without background dropped,
	targets weighted 1 and each bin's background weighted so that it holds
	the bin's share of the targets."""

	bins = CPG_BINS if cpg else GC_BINS
	members = {}
	for name, seq in records:
		if cpg:
			value = seq.count('CG') / max(len(seq) - 1, 0.01)
		else:
			value = sum(1 for b in seq if b in 'CG') / len(seq)
		value = float('%f' % value)
		index = next((i for i, edge in enumerate(bins) if value < edge), len(bins))
		members.setdefault(index, []).append(name)
	kept = {}
	for index, names in members.items():
		p = sum(1 for n in names if not n.startswith('BG'))
		if p >= 1 and len(names) - p >= 1:
			kept[index] = (p, len(names))
	total = sum(a for _, a in kept.values())
	positives = sum(p for p, _ in kept.values())
	weights = {}
	for index, (p, a) in kept.items():
		w = p / positives * (total - positives) / (a - p)
		for n in members[index]:
			weights[n] = ('1', '1.00000') if not n.startswith('BG') else ('0', '%.5f' % w)
	return weights


# --- known motif enrichment -----------------------------------------------------------------

def test_known_counts_without_weights(runs):
	"""-noweight: each motif's counts are the numbers of targets and of
	background sequences the reference scan finds a site in; the log p-value
	is the hypergeometric tail, the q-value Benjamini-Hochberg's, the p-value
	text 1e<int(logp / 2.3026)>, and rows are sorted by p-value."""

	records = input_records()
	group = read_group(runs['noweight'] / 'group.adj')
	assert sorted(g for g, _ in group.values()) == [0] * 1000 + [1] * 200
	rows = hh.read_known_results(str(runs['noweight'] / 'knownResults.txt'))
	motifs = hh.read_motifs(KNOWN)
	assert sorted(r['name'] for r in rows) == sorted(m.name for m in motifs)
	got = {r['name']: r for r in rows}
	logps = {}
	for m in motifs:
		t, b = counts(records, group, m, 'input')
		row = got[m.name]
		assert (row['targets'], row['background']) == (t, b), m.name
		assert row['target_total'] == 200
		assert row['target_percent'] == pytest.approx(100 * t / 200, abs=0.005)
		assert row['background_percent'] == pytest.approx(100 * b / 1000, abs=0.005)
		logps[m.name] = hh.log_hypergeometric_tail(int(t), int(t + b), 200, 1200)
		assert row['logp'] == pytest.approx(logps[m.name], rel=1e-3, abs=1e-3), m.name
		assert row['pvalue'] == '1e{}'.format(int(logps[m.name] / 2.3026))
	q = benjamini([math.exp(logps[m.name]) for m in motifs])
	for m, value in zip(motifs, q):
		assert got[m.name]['qvalue'] == pytest.approx(value, abs=1.5e-4)
	assert [r['logp'] for r in rows] == sorted(r['logp'] for r in rows)
	assert got['siteA']['targets'] == 120


def test_known_counts_with_weights(runs):
	"""By default background sequences are weighted (GC bins, then oligo
	normalization): a background count is the sum of the weights, in the
	group file the run kept, of the background sequences with a site, and
	the p-value uses the counts and totals truncated to integers."""

	out = runs['known']
	records = read_seq_tsv(out / 'seq.tsv')
	assert records == input_records()
	group = read_group(out / 'group.adj')
	assert {w for g, w in group.values() if g == 1} == {1.0}
	assert len({w for g, w in group.values() if g == 0}) > 10
	background_total = sum(w for g, w in group.values() if g == 0)
	got = known_by_name(out)
	for m in hh.read_motifs(KNOWN):
		t, b = counts(records, group, m, 'input')
		row = got[m.name]
		assert row['targets'] == t
		assert row['background'] == pytest.approx(b, abs=0.051), m.name
		assert row['background_percent'] == pytest.approx(100 * b / background_total, abs=0.006)
		want = hh.log_hypergeometric_tail(int(t), int(t) + int(b), 200,
			200 + int(background_total))
		assert row['logp'] == pytest.approx(want, rel=1e-3, abs=1e-3), m.name


def test_known_motif_files_for_significant_motifs(runs):
	"""knownResults/known<N>.motif is written for each motif with p-value at
	most 0.01, in the order of knownResults.txt, with its statistics."""

	out = runs['known']
	rows = hh.read_known_results(str(out / 'knownResults.txt'))
	significant = [r for r in rows if math.exp(r['logp']) <= 0.01]
	assert [r['name'] for r in significant] == ['siteA', 'siteB']
	known = {m.name: m for m in hh.read_motifs(KNOWN)}
	for n, row in enumerate(significant, 1):
		m = hh.read_motifs(str(out / 'knownResults' / 'known{}.motif'.format(n)))[0]
		assert m.name == row['name']
		assert m.rows == known[m.name].rows
		assert m.threshold == known[m.name].threshold
		assert m.header[2].startswith('T:{:.1f}('.format(row['targets']))
		assert (out / 'knownResults' / 'known{}.logo.svg'.format(n)).exists()
	assert not (out / 'knownResults' / 'known{}.motif'.format(len(significant) + 1)).exists()
	assert (out / 'knownResults.html').exists()


@pytest.mark.parametrize('key', ['known_fastaBg', 'known_homer2'])
def test_known_results_do_not_depend_on_alias(runs, key):
	"""-fastaBg is -fasta and -homer2 is the default."""

	a = (runs['known'] / 'knownResults.txt').read_text()
	assert (runs[key] / 'knownResults.txt').read_text() == a


@pytest.mark.parametrize('key', ['known_p2', 'known_cpu'])
def test_known_counts_do_not_depend_on_threads(runs, key):
	"""-p 2 (or its alias -cpu 2) finds the same motifs in the same
	sequences as one thread."""

	columns = ('name', 'targets', 'target_percent', 'background', 'background_percent')
	one = hh.read_known_results(str(runs['known'] / 'knownResults.txt'))
	two = hh.read_known_results(str(runs[key] / 'knownResults.txt'))
	assert sorted(tuple(r[c] for c in columns) for r in two) == \
		sorted(tuple(r[c] for c in columns) for r in one)


def test_binomial_statistic(runs):
	"""-b: log p-value is the binomial tail of the target count with the
	background fraction (at least 1 / background total) as probability."""

	records = input_records()
	group = read_group(runs['noweight'] / 'group.adj')
	got = known_by_name(runs['binomial'])
	for m in hh.read_motifs(KNOWN):
		t, b = counts(records, group, m, 'input')
		assert (got[m.name]['targets'], got[m.name]['background']) == (t, b)
		want = hh.log_binomial_tail(int(t), 200, max(b / 1000.0, 1 / 1000.0)) if t else 0.0
		assert got[m.name]['logp'] == pytest.approx(want, rel=1e-3, abs=1e-3), m.name
	default = known_by_name(runs['noweight'])
	assert got['siteA']['logp'] < default['siteA']['logp'] - 100


def test_norevopp_counts_the_forward_strand_only(runs):
	records = input_records()
	group = read_group(runs['noweight'] / 'group.adj')
	got = known_by_name(runs['norevopp'])
	both = known_by_name(runs['noweight'])
	for m in hh.read_motifs(KNOWN):
		t, b = counts(records, group, m, 'input', strand='+')
		assert (got[m.name]['targets'], got[m.name]['background']) == (t, b), m.name
	assert got['siteA']['targets'] < both['siteA']['targets']


def test_scrambled_background_without_fasta(runs):
	"""Without a background FASTA file the background is scrambleBg.fasta:
	five sequences per target, each as long as its target and drawn from its
	letters; known counts use it as the background."""

	out = runs['scrambled']
	targets = hh.read_fasta(hh.TARGETS)
	scrambled = hh.read_fasta(str(out / 'scrambleBg.fasta'))
	assert [name for name, _ in scrambled] == [str(i) for i in range(1, 1001)]
	for i, (_, seq) in enumerate(scrambled):
		source = targets[i // 5][1]
		assert len(seq) == len(source)
		assert set(seq) <= set(source)
	assert len({seq for _, seq in scrambled}) == 1000
	records = [(n, hh.clean_sequence(s)) for n, s in targets]
	records += [('BG' + n, hh.clean_sequence(s)) for n, s in scrambled]
	group = {n: (0 if n.startswith('BG') else 1, 1.0) for n, _ in records}
	got = known_by_name(out)
	for m in hh.read_motifs(KNOWN):
		t, b = counts(records, group, m, 'scrambled')
		assert (got[m.name]['targets'], got[m.name]['background']) == (t, b), m.name


def test_float_optimizes_known_thresholds(runs):
	"""-float lets the known motif thresholds move to improve each p-value:
	no motif gets a worse p-value, the written motifs carry the new
	thresholds, and the counts are those of a scan at that threshold."""

	floated = known_by_name(runs['float'])
	default = known_by_name(runs['known'])
	for name, row in floated.items():
		assert row['logp'] <= default[name]['logp'] + 1e-9
	assert floated['siteB']['logp'] < default['siteB']['logp'] - 10
	records = read_seq_tsv(runs['known'] / 'seq.tsv')
	group = read_group(runs['known'] / 'group.adj')
	known = {m.name: m for m in hh.read_motifs(KNOWN)}
	files = sorted((runs['float'] / 'knownResults').glob('known*.motif'))
	assert files
	changed = 0
	for path in files:
		m = hh.read_motifs(str(path))[0]
		changed += m.threshold != known[m.name].threshold
		t, b = counts(records, group, m, 'input')
		assert floated[m.name]['targets'] == t
		assert floated[m.name]['background'] == pytest.approx(b, abs=0.051)
	assert changed


def test_mask_motif_masks_before_counting(runs):
	"""-maskMotif <file> <file>: every site of the given motifs becomes N in
	the analysed sequences (as homer2 mask scans), then the known motifs are
	counted on the masked sequences."""

	out = runs['mask']
	a, b = hh.read_motifs(SITE_A)[0], hh.read_motifs(SITE_B)[0]
	want = [(n, hh.mask_sequence(hh.mask_sequence(s, a), b)) for n, s in input_records()]
	assert sorted(read_seq_tsv(out / 'seq.tsv')) == sorted(want)
	group = read_group(out / 'group.adj')
	got = known_by_name(out)
	for m in hh.read_motifs(KNOWN):
		t, bg = counts(want, group, m, 'masked')
		assert (got[m.name]['targets'], got[m.name]['background']) == (t, bg), m.name
	assert got['siteA']['targets'] == got['siteB']['targets'] == 0


def test_rand_shuffles_labels(runs):
	"""-rand: the same sequences with 200 target labels given at random;
	counts follow the new labels."""

	group = read_group(runs['rand'] / 'group.adj')
	records = input_records()
	assert sorted(group) == sorted(n for n, _ in records)
	targets = {n for n, (g, _) in group.items() if g == 1}
	assert len(targets) == 200
	assert targets != {n for n, _ in records if not n.startswith('BG')}
	got = known_by_name(runs['rand'])
	for m in hh.read_motifs(KNOWN):
		t, b = counts(records, group, m, 'input')
		assert (got[m.name]['targets'], got[m.name]['background']) == (t, b), m.name


def test_mset_selects_the_library(runs):
	"""-mset worms: known enrichment uses data/knownTFs/worms/known.motifs
	and de novo motifs are compared with worms/all.motifs."""

	out = runs['mset']
	known = hh.read_motifs(os.path.join(KNOWN_TFS, 'worms', 'known.motifs'))
	rows = hh.read_known_results(str(out / 'knownResults.txt'))
	assert sorted(r['name'] for r in rows) == sorted(m.name for m in known)
	library = {m.name for m in hh.read_motifs(os.path.join(KNOWN_TFS, 'worms', 'all.motifs'))}
	guesses = re.findall(r',BestGuess:(.*)\([-\d.]+\)\t', (out / 'nonRedundant.motifs').read_text())
	assert guesses and set(guesses) <= library


def test_rna_known_uses_the_rna_library_on_one_strand(runs):
	"""-rna turns known enrichment off (-known turns it back on) and searches
	the given strand only; the known motifs are data/knownTFs/known.rna.motifs."""

	got = hh.read_known_results(str(runs['rna_known'] / 'knownResults.txt'))
	library = hh.read_motifs(os.path.join(KNOWN_TFS, 'known.rna.motifs'))
	assert sorted(r['name'] for r in got) == sorted(m.name for m in library)
	records = input_records()
	group = {n: (0 if n.startswith('BG') else 1, 1.0) for n, _ in records}
	got = {r['name']: r for r in got}
	for m in library[:3]:
		t, b = counts(records, group, m, 'input', strand='+')
		assert (got[m.name]['targets'], got[m.name]['background']) == (t, b), m.name
	assert not (runs['rna'] / 'knownResults.txt').exists()


# --- weights and normalization -----------------------------------------------------------

@pytest.mark.parametrize('key,cpg', [('nlen0', False), ('cpg', True), ('CpG', True)])
def test_gc_weights(runs, key, cpg):
	"""-nlen 0 skips oligo normalization, so the group file holds the GC
	weights (CpG weights with -cpg or -CpG) and no seq.autonorm.tsv is
	written."""

	out = runs[key]
	records = read_seq_tsv(out / 'seq.tsv')
	got = {row[0]: (row[1], row[2]) for row in hh.read_tsv(str(out / 'group.adj'))}
	assert got == gc_weights(records, cpg)
	assert not (out / 'seq.autonorm.tsv').exists()
	assert len(got) < len(records)


@pytest.mark.parametrize('key,norm', [('known', 'nlen3'), ('nlen2', 'nlen2'), ('nmax5', 'nmax5'),
	('neutral', 'neutral')])
def test_normalization_options_reach_homer2_norm(runs, direct, key, norm):
	"""-nlen, -nmax and -neutral: the weights equal those of homer2 norm run
	with the same options on the GC-weighted sequences, and the oligo table
	lists oligos of length 1 to -nlen."""

	got = hh.read_tsv(str(runs[key] / 'group.adj'))
	want = hh.read_tsv(str(direct[norm] / 'group.adj'))
	assert sorted(got) == sorted(want)
	table = hh.read_tsv(str(runs[key] / 'seq.autonorm.tsv'))[1:]
	nlen = 2 if key == 'nlen2' else 3
	assert sorted({len(row[0]) for row in table}) == list(range(1, nlen + 1))
	assert (runs[key] / 'seq.autonorm.tsv').read_text() == \
		(direct[norm] / 'autonorm.tsv').read_text()


def test_normalization_options_change_weights(runs):
	weights = {key: sorted(hh.read_tsv(str(runs[key] / 'group.adj')))
		for key in ('known', 'nlen2', 'nmax5', 'neutral', 'nlen0')}
	for key in ('nlen2', 'nmax5', 'neutral', 'nlen0'):
		assert weights[key] != weights['known'], key
	neutral = read_group(runs['neutral'] / 'group.adj')
	assert {w for g, w in neutral.values() if g == 1} != {1.0}


def test_keep_files(runs):
	"""Without -keepFiles the sequence and group files are removed."""

	for name in ('seq.tsv', 'group.adj'):
		assert (runs['known'] / name).exists()
		assert not (runs['known_fastaBg'] / name).exists()


# --- de novo motif finding ---------------------------------------------------------------

def test_de_novo_recovers_planted_site_a(runs):
	"""The best motif (-len 8) is part of siteA: its sites are in exactly
	the 120 peaks siteA was planted in, it best matches siteA in the -mcheck
	library, and its T: and B: counts are those of the reference scan."""

	out = runs['denovo']
	motifs = hh.read_motifs(str(out / 'homerMotifs.motifs8'))
	assert len(motifs) == 2 and all(len(m) == 8 for m in motifs)
	records = read_seq_tsv(out / 'seq.tsv')
	group = read_group(out / 'group.adj')
	planted = {p['peak'] for p in hh.read_planted() if p['motif'] == 'siteA'}
	found = has_sites('input', records, motifs[0])
	assert {n for n, _ in records if not n.startswith('BG') and found[n]} == planted
	for m in motifs:
		t, b = counts(records, group, m, 'input')
		stats = re.match(r'T:([\d.]+)\([\d.]+%\),B:([\d.]+)\(', m.header[2])
		assert float(stats.group(1)) == t
		assert float(stats.group(2)) == pytest.approx(b, abs=0.051)
	text = (out / 'nonRedundant.motifs').read_text()
	assert re.findall(r',BestGuess:([^(]*)\(', text)[0] == 'siteA'


def test_de_novo_output_files(runs):
	"""homerMotifs.all.motifs gathers every length's motifs and
	homerResults.html shows the motifs left after comparing them."""

	out = runs['denovo']
	assert (out / 'homerMotifs.all.motifs').read_text() == \
		(out / 'homerMotifs.motifs8').read_text()
	page = (out / 'homerResults.html').read_text()
	kept = len(hh.read_motifs(str(out / 'nonRedundant.motifs')))
	for n in range(1, kept + 1):
		assert 'homerResults/motif{}.info.html'.format(n) in page
		assert (out / 'homerResults' / 'motif{}.motif'.format(n)).exists()
		assert (out / 'homerResults' / 'motif{}.similar.html'.format(n)).exists()
	for path in [out / 'homerResults.html'] + sorted((out / 'homerResults').glob('*.html')):
		assert re.findall(r'<P>(.*?)</P>', path.read_text(), re.S) == ['']


def test_lengths_and_motif_count(runs):
	"""-len 8,10 writes one motif file per length, -S 1 one motif in each."""

	out = runs['lens']
	for length in (8, 10):
		motifs = hh.read_motifs(str(out / 'homerMotifs.motifs{}'.format(length)))
		assert [len(m) for m in motifs] == [length]
	assert (out / 'homerMotifs.all.motifs').read_text() == \
		(out / 'homerMotifs.motifs10').read_text() + (out / 'homerMotifs.motifs8').read_text()


def test_nocheck_and_noknown(runs):
	out = runs['mis0']
	assert (out / 'homerMotifs.motifs8').exists()
	for name in ('homerResults.html', 'nonRedundant.motifs', 'homerMotifs.all.motifs',
			'knownResults.txt'):
		assert not (out / name).exists()


def test_nomotif(runs):
	out = runs['known']
	assert (out / 'knownResults.txt').exists()
	assert not list(out.glob('homerMotifs*'))
	assert not (out / 'homerResults.html').exists()


@pytest.mark.parametrize('key', sorted(DENOVO_OPTIONS))
def test_de_novo_options_reach_homer2(runs, direct, key):
	"""-mis, -e, -quickMask, -minlp, -cache, -norevopp (-strand +) and -b
	(-stat binomial): the motifs equal homer2 denovo's with that option; both
	sides get -e 0.01 (the default findMotifs.pl's usage gives) except in the
	-e case."""

	assert (runs[key] / 'homerMotifs.motifs8').read_text() == \
		(direct[key] / 'direct.motifs').read_text()


@pytest.mark.parametrize('key', ['mis0', 'e', 'quickMask', 'denovo_norevopp', 'denovo_binomial'])
def test_de_novo_options_change_the_motifs(runs, key):
	"""Each option changes the motifs of a run that differs from it only in
	that option (-e 0.01 for the -e 0.0001 run)."""

	assert (runs[key] / 'homerMotifs.motifs8').read_text() != \
		(runs['denovo_e'] / 'homerMotifs.motifs8').read_text()


def test_fdr_adds_empirical_fdr(runs):
	"""-fdr 2: denovo runs twice more on randomized labels (kept under
	randomizations/) and each motif's statistics gain FDR:<#>, the empirical
	FDR of its log p-value among the randomized runs' (Statistics::empiricalFDR2)."""

	out = runs['fdr']
	real = hh.read_motifs(str(out / 'homerMotifs.motifs8'))
	random_logp = []
	for i in range(2):
		path = out / 'randomizations' / 'homerMotifs.r{}.motifs8'.format(i)
		random_logp += [float(m.header[0]) for m in hh.read_motifs(str(path))]
	random_logp.sort()
	logp = [float(m.header[0]) for m in real]
	order = sorted(range(len(logp)), key=lambda i: logp[i])
	want, last, index = {}, -1e10, 0
	for rank, i in enumerate(order):
		while index < len(random_logp) and random_logp[index] < logp[i]:
			index += 1
		value = min(1.0, max(last, (index + 1) / 2.0 / (rank + 1)))
		want[i], last = value, value
	for i, m in enumerate(real):
		assert m.header[2].endswith(',FDR:' + '%.15g' % want[i])
	plain = hh.read_motifs(str(runs['denovo'] / 'homerMotifs.motifs8'))
	assert [m.rows for m in real] == [m.rows for m in plain]


def test_reuse_keeps_the_old_group_file(direct):
	"""-reuse runs de novo finding on the sequence and group files left in
	the output directory (here with every weight edited to 1) instead of
	making new ones."""

	got = (direct['reuse'] / 'out' / 'homerMotifs.motifs8').read_text()
	assert got == (direct['reuse_direct'] / 'direct.motifs').read_text()
	weights = {row[2] for row in hh.read_tsv(str(direct['reuse'] / 'out' / 'group.adj'))}
	assert weights == {'1.000000'}


def test_homer1_uses_the_original_homer(runs):
	"""-homer1: motifs from the original homer, whose statistics column
	starts with the total and target weights."""

	out = runs['homer1']
	group = read_group(out / 'group.adj')
	total = '%.1f' % sum(w for _, w in group.values())
	motifs = hh.read_motifs(str(out / 'homerMotifs.motifs8'))
	assert motifs and all(len(m) == 8 for m in motifs)
	for m in motifs:
		fields = m.header[2].split(',')
		assert fields[:2] == [total, '200.0']
		assert not m.header[2].startswith('T:')


def letter_heights(svg):
	found = re.findall(r'transform="matrix\([^,]+,0,0,([^,]+),[^)]+\)">(\w)</text>', svg)
	return [[(letter, float(h)) for h, letter in found[i:i + 4]] for i in range(0, len(found), 4)]


def expected_heights(rows, bits=False, letters='ACGT'):
	columns = []
	for row in rows:
		if bits:
			h = -sum(max(v, 0.001) * math.log(max(v, 0.001)) / math.log(2) for v in row)
			row = [v * (2 - h) / 2 for v in row]
		order = sorted(range(4), key=lambda j: row[j])
		columns.append([(letters[j], row[j]) for j in order])
	return columns


def assert_heights(got, want):
	assert [[x[0] for x in c] for c in got] == [[x[0] for x in c] for c in want]
	assert [x[1] for c in got for x in c] == pytest.approx([x[1] for c in want for x in c],
		rel=1e-6, abs=1e-9)


@pytest.mark.parametrize('key,bits', [('denovo', False), ('bits', True)])
def test_logos(runs, key, bits):
	"""Logos show each base's probability; -bits scales them by the
	column's information content, p * (2 - H) / 2 with H in bits."""

	out = runs[key]
	motifs = hh.read_motifs(str(out / 'nonRedundant.motifs'))
	assert motifs
	for n, m in enumerate(motifs, 1):
		svg = (out / 'homerResults' / 'motif{}.logo.svg'.format(n)).read_text()
		assert_heights(letter_heights(svg), expected_heights(m.rows, bits))


def test_basic_skips_known_comparison(runs):
	"""-basic: the same non-redundant motifs, without best guesses or
	similar-motif pages."""

	out = runs['basic']
	text = (out / 'nonRedundant.motifs').read_text()
	assert 'BestGuess' not in text
	checked = (runs['denovo'] / 'nonRedundant.motifs').read_text()
	assert text == re.sub(r',BestGuess:[^\t]*', '', checked)
	assert not list((out / 'homerResults').glob('*.similar.html'))


def test_rna_mode(runs, direct):
	"""-rna: de novo on the given strand only (homer2 denovo -strand +, and
	the counts are those of a + strand scan), logos with U, compared with
	the RNA library data/knownTFs/all.rna.motifs, no known enrichment."""

	out = runs['rna']
	assert (out / 'homerMotifs.motifs8').read_text() == \
		(direct['rna'] / 'direct.motifs').read_text()
	motif = hh.read_motifs(str(out / 'homerMotifs.motifs8'))[0]
	records = read_seq_tsv(out / 'seq.tsv')
	group = read_group(out / 'group.adj')
	t, b = counts(records, group, motif, 'input', strand='+')
	stats = re.match(r'T:([\d.]+)\([\d.]+%\),B:([\d.]+)\(', motif.header[2])
	assert float(stats.group(1)) == t
	assert float(stats.group(2)) == pytest.approx(b, abs=0.051)
	library = {m.name for m in hh.read_motifs(os.path.join(KNOWN_TFS, 'all.rna.motifs'))}
	guesses = re.findall(r',BestGuess:(.*)\([-\d.]+\)\t', (out / 'nonRedundant.motifs').read_text())
	assert guesses and set(guesses) <= library
	svg = (out / 'homerResults' / 'motif1.logo.svg').read_text()
	assert_heights(letter_heights(svg), expected_heights(motif.rows, letters='ACGU'))
	assert not (out / 'knownResults.txt').exists()


# --- scanning and sequence files -----------------------------------------------------------

def find_rows(stdout):
	lines = stdout.splitlines()
	assert lines[0].split('\t') == ['FASTA ID', 'Offset', 'Sequence', 'Motif Name', 'Strand',
		'MotifScore']
	rows = []
	for line in lines[1:]:
		cols = line.split('\t')
		rows.append((cols[0], int(cols[1]), cols[2], cols[3], cols[4], float(cols[5])))
	return sorted(rows)


@pytest.mark.parametrize('key,strand', [('find', 'both'), ('find_norevopp', '+')])
def test_find_prints_every_site(runs, key, strand):
	"""-find <motif file> prints a header and every site of the motifs in the
	targets (positions relative to the sequence middle), with -norevopp on
	the + strand only, and runs nothing else."""

	stdout = (runs[key].parent / 'stdout.txt').read_text()
	got = find_rows(stdout)
	want = []
	for m in hh.read_motifs(SITE_AB):
		for name, seq in hh.read_fasta(hh.TARGETS):
			for pos, site, site_strand, score in hh.scan(seq, m, strand):
				want.append((name, pos, site, m.name, site_strand, score))
	want.sort()
	assert [r[:5] for r in got] == [w[:5] for w in want]
	assert [r[5] for r in got] == pytest.approx([w[5] for w in want], abs=1e-6)
	assert not (runs[key] / 'knownResults.txt').exists()
	assert not list(runs[key].glob('homerMotifs*'))


def test_dump_fasta(runs, inputs):
	"""-dumpFasta writes the analysed targets and background (upper case,
	background names prefixed with "BG") as target.fa and background.fa."""

	out = runs['dump']
	want = {n: hh.clean_sequence(s) for n, s in hh.read_fasta(hh.TARGETS)}
	assert dict(hh.read_fasta(str(out / 'target.fa'))) == want
	want = {'BG' + n: hh.clean_sequence(s) for n, s in hh.read_fasta(inputs['background_200'])}
	assert dict(hh.read_fasta(str(out / 'background.fa'))) == want


def test_min_max_and_chopify(runs, inputs):
	"""-min 50 -max 250 keep targets of 50 to 250 bp; -chopify cuts each
	background sequence into pieces as long as the kept targets' mean length
	(rounded down), named <name>-chopped<#>, as many as the length divided by
	that mean, rounded."""

	out = runs['dump_cut']
	kept = {n: hh.clean_sequence(s) for n, s in hh.read_fasta(inputs['varied'])
		if 50 <= len(s) <= 250}
	assert dict(hh.read_fasta(str(out / 'target.fa'))) == kept
	assert {'t_min', 't_max'} <= set(kept) and len(kept) == 22
	mean = sum(len(s) for s in kept.values()) // len(kept)
	want = {}
	for name, seq in hh.read_fasta(inputs['background_200']):
		i = 0
		while i < len(seq) / mean - 0.5:
			want['BG{}-chopped{}'.format(name, i + 1)] = hh.clean_sequence(seq[i * mean:(i + 1) * mean])
			i += 1
	assert dict(hh.read_fasta(str(out / 'background.fa'))) == want
	assert len(want) > 200

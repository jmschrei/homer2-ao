"""homer2 bg (alias homer2 background): selecting background sequences that
match the GC content and k-mer composition of a set of targets.

The targets are the 200 peaks of tests/data/peaks.txt (-p, with -g
genome.fa) or tests/data/targets.fa (-i); the background comes from the
genome, from explicit positions (-bg), from genome ranges (-bgr), from a FASTA
file (-b) or from a model of the targets (-model). Every run passes -seed and
happens once, in a module fixture; most screen a sample of the genome
(-NN 20000) to stay fast.

Expected values are derived independently:

- every selected sequence equals the genome at its reported coordinates
  (reverse-complemented on the - strand, soft-masked bases as N with -mask),
  extracted in Python;
- no selected window holds a larger fraction of N than -maxFractionN, none
  overlaps a target or another selected window unless the options allow it,
  and when overlaps are allowed the count exceeds what non-overlapping
  windows could reach (a count computed from the genome or region lengths);
- the GC bins: thresholds are the GC fractions at every 1/nbins quantile of
  the targets, computed in Python; every sequence's reported GC and bin
  follow from its own sequence; every bin holds the same number of selected
  sequences, and the selection's GC distribution passes a two-sample
  Kolmogorov-Smirnov test against the targets' at alpha = 0.001;
- the k-mer frequency files and the profile scores of -pscore are computed in
  Python from the sequences.
"""

import bisect
import collections
import math
import os
import random

import pytest

import homer_helpers as hh


SEED = ['-seed', 1]
SCREEN = ['-NN', 20000]
OUTPUTS = ['bg.sequences.fasta', 'bg.stats.txt', 'bg.positions.bed', 'target.sequences.fasta',
	'target.stats.txt', 'target.positions.bed', 'bg.freq.txt', 'target.freq.txt', 'info.txt']
N_RUN = ('chr2', 40000, 42000)
NRUN_PEAKS = {'nrun25': 50, 'nrun60': 120, 'nrun90': 180}
REGIONS = [('chr1', 20000, 23000), ('chr2', 60000, 63000), ('chr3', 36000, 39000)]
WINDOW = 200


def write_inputs(top):
	"""Inputs written for this module: peaks overlapping the N run, explicit
	background positions, background ranges and targets of two lengths."""

	paths = {}
	peaks = hh.read_peaks()
	genome = hh.read_genome()
	path = top / 'nrun_peaks.txt'
	with open(str(path), 'w') as out:
		for p in peaks[:150]:
			out.write('{}\t{}\t{}\t{}\t{}\n'.format(p['id'], p['chr'], p['start'], p['end'],
				p['strand']))
		# a 200 bp peak ending k bases into the N run holds k Ns
		for name, k in sorted(NRUN_PEAKS.items()):
			start0 = N_RUN[1] - (WINDOW - k)
			out.write('{}\t{}\t{}\t{}\t+\n'.format(name, N_RUN[0], start0 + 1, start0 + WINDOW))
	paths['nrun_peaks'] = str(path)

	rng = random.Random(7)
	rows = []
	while len(rows) < 400:
		chrom = rng.choice(sorted(genome))
		start = rng.randrange(0, len(genome[chrom]) - WINDOW)
		window = genome[chrom][start:start + WINDOW]
		if 'N' in window or any(p['chr'] == chrom and start < p['end'] and p['start'] - 1 <
				start + WINDOW for p in peaks):
			continue
		rows.append((chrom, start, start + WINDOW, 'pos{}'.format(len(rows)), 0, rng.choice('+-')))
	path = top / 'bg_positions.bed'
	with open(str(path), 'w') as out:
		for row in rows:
			out.write('\t'.join(str(v) for v in row) + '\n')
	paths['bg_positions'] = str(path)

	path = top / 'regions.bed'
	with open(str(path), 'w') as out:
		for i, (chrom, start, end) in enumerate(REGIONS):
			out.write('{}\t{}\t{}\tregion{}\t0\t+\n'.format(chrom, start, end, i))
	paths['regions'] = str(path)

	targets = hh.read_fasta(hh.TARGETS)
	mixed = [(name, seq[:180] if i % 2 else seq) for i, (name, seq) in enumerate(targets)]
	paths['mixed'] = hh.write_fasta(str(top / 'mixed.fa'), mixed)
	return paths


def run_table(paths):
	peaks = ['-p', hh.PEAKS, '-g', hh.GENOME]
	big = peaks + ['-N', 3000, '-NN', 60000] + SEED
	return {
		'peaks': peaks + ['-N', 300, '-gs'] + SCREEN + SEED,
		'seed_again': peaks + ['-N', 300, '-gs'] + SCREEN + SEED,
		'alias': peaks + ['-N', 300, '-gs'] + SCREEN + SEED,
		'seed2': peaks + ['-N', 300] + SCREEN + ['-seed', 2],
		'peaks_bed': ['-p', hh.PEAKS_BED, '-g', hh.GENOME, '-N', 300, '-gs'] + SCREEN + SEED,
		'default_prefix': peaks + ['-N', 30] + SCREEN + SEED,
		'fasta_targets': ['-i', hh.TARGETS, '-g', hh.GENOME, '-N', 300] + SCREEN + SEED,
		'mixed_size': ['-i', paths['mixed'], '-g', hh.GENOME, '-N', 150] + SCREEN + SEED,
		'size100': peaks + ['-N', 150, '-size', 100] + SCREEN + SEED,
		'screen_small': peaks + ['-N', 3000, '-NN', 1000] + SEED,
		'big': big,
		'big_allow_bg': big + ['-allowBgOverlaps'],
		'big_strand': big + ['-strand'],
		'big_overlap_late': big + ['-overlapIteration', 25],
		'overlap_first': peaks + ['-N', 300, '-overlapIteration', 1] + SCREEN + SEED,
		'allow_target': peaks + ['-N', 300, '-allowTargetOverlaps'] + SCREEN + SEED,
		'mask': peaks + ['-N', 300, '-mask', '-maxFractionN', 0] + SCREEN + SEED,
		'no_n': peaks + ['-N', 300, '-maxFractionN', 0] + SCREEN + SEED,
		'nrun_default': ['-p', paths['nrun_peaks'], '-g', hh.GENOME, '-N', 150] + SCREEN + SEED,
		'nrun_02': ['-p', paths['nrun_peaks'], '-g', hh.GENOME, '-N', 150, '-maxFractionN', 0.2]
			+ SCREEN + SEED,
		'nrun_095': ['-p', paths['nrun_peaks'], '-g', hh.GENOME, '-N', 150, '-maxFractionN', 0.95]
			+ SCREEN + SEED,
		'nbins5': peaks + ['-N', 300, '-nbins', 5] + SCREEN + SEED,
		'nsub5': peaks + ['-N', 300, '-nsubBins', 5] + SCREEN + SEED,
		'pkmer2': peaks + ['-N', 300, '-pkmer', 2] + SCREEN + SEED,
		'ikmer3': peaks + ['-N', 300, '-ikmer', 3] + SCREEN + SEED,
		'include_ns': peaks + ['-N', 300, '-includeNs'] + SCREEN + SEED,
		'exclude_ns': peaks + ['-N', 300, '-excludeNs'] + SCREEN + SEED,
		'decay': peaks + ['-N', 300, '-decayRate', 0.1] + SCREEN + SEED,
		'model': ['-i', hh.TARGETS, '-model', '-N', 300] + SEED,
		'model_exclude_ns': ['-i', hh.TARGETS, '-model', '-N', 300, '-excludeNs'] + SEED,
		'bg_fasta': ['-i', hh.TARGETS, '-b', hh.BACKGROUND, '-N', 300, '-pscore', 'pscore.txt']
			+ SEED,
		'bg_positions': ['-i', hh.TARGETS, '-g', hh.GENOME, '-bg', paths['bg_positions'],
			'-N', 300] + SEED,
		'bg_regions': peaks + ['-bgr', paths['regions'], '-N', 300, '-nbins', 1] + SEED,
		'bg_regions_allow': peaks + ['-bgr', paths['regions'], '-N', 300, '-nbins', 1,
			'-allowBgOverlaps'] + SEED,
	}


@pytest.fixture(scope='module')
def inputs(tmp_path_factory):
	return write_inputs(tmp_path_factory.mktemp('bg_inputs'))


@pytest.fixture(scope='module')
def table(inputs):
	return run_table(inputs)


@pytest.fixture(scope='module')
def runs(tmp_path_factory, table):
	"""key -> (directory, output prefix) of every run of the module."""

	top = tmp_path_factory.mktemp('bg')
	calls = {}
	for key, args in table.items():
		command = 'background' if key == 'alias' else 'bg'
		output = [] if key == 'default_prefix' else ['-o', 'sel']
		calls[key] = (['homer2', command] + args + output, top / key)
	hh.run_many(calls)
	return {key: (top / key, 'out' if key == 'default_prefix' else 'sel') for key in calls}


def out_path(runs, key, name):
	cwd, prefix = runs[key]
	return str(cwd / '{}.{}'.format(prefix, name))


def read_output(runs, key, name):
	with open(out_path(runs, key, name)) as handle:
		return handle.read()


def read_selection(runs, key, which='bg'):
	"""Rows of <which>.positions.bed joined with the sequences and stats of the
	same name: dicts with chr, start, end, strand, seq, gc, bin, pscore."""

	bed = hh.read_tsv(out_path(runs, key, which + '.positions.bed')) \
		if os.path.getsize(out_path(runs, key, which + '.positions.bed')) else []
	fasta = hh.read_fasta(out_path(runs, key, which + '.sequences.fasta'))
	_, stats = hh.read_tsv(out_path(runs, key, which + '.stats.txt'), header=True)
	assert len(bed) == len(fasta) == len(stats)
	rows = []
	for b, (name, seq), st in zip(bed, fasta, stats):
		assert name.split(' ')[0] == b[3] == st[0]
		assert b[5] == st[3]
		rows.append({'chr': b[0], 'start': int(b[1]), 'end': int(b[2]), 'strand': b[5],
			'seq': seq, 'pscore': float(st[5]), 'gc': float(st[6]), 'bin': int(st[7])})
	return rows


# --- reference computations ---------------------------------------------------------------

def gc_fraction(seq):
	"""GC over the bases that are not N, as HOMER computes it."""

	bases = sum(1 for b in seq if b != 'N')
	return sum(1 for b in seq if b in 'GC') / float(bases)


def genome_sequence(chrom, start, end, strand, mask=False):
	seq = hh.read_genome()[chrom][start:end]
	if mask:
		seq = ''.join('N' if b.islower() else b for b in seq)
	seq = seq.upper()
	return hh.revcomp(seq) if strand == '-' else seq


def peak_intervals(path=hh.PEAKS):
	return [(p['chr'], p['start'] - 1, p['end'], p['strand']) for p in hh.read_peaks(path)]


def overlapping_pairs(intervals):
	"""Pairs of (chr, start, end, strand) intervals that share a base."""

	by_chr = collections.defaultdict(list)
	for iv in intervals:
		by_chr[iv[0]].append(iv)
	pairs = []
	for ivs in by_chr.values():
		ivs.sort(key=lambda iv: iv[1])
		for i in range(len(ivs)):
			for j in range(i + 1, len(ivs)):
				if ivs[j][1] >= ivs[i][2]:
					break
				pairs.append((ivs[i], ivs[j]))
	return pairs


def overlaps_with(rows, intervals):
	"""(row, interval) pairs where a selected window shares a base with an
	interval."""

	by_chr = collections.defaultdict(list)
	for iv in intervals:
		by_chr[iv[0]].append(iv)
	hits = []
	for row in rows:
		for iv in by_chr[row['chr']]:
			if row['start'] < iv[2] and iv[1] < row['end']:
				hits.append((row, iv))
	return hits


def capacity(lengths):
	"""Most 200 bp windows without a shared base that fit in sequences of
	these lengths (a window may start anywhere before the last 200 bases)."""

	return sum((n - WINDOW - 1) // WINDOW + 1 for n in lengths)


def bin_thresholds(target_gcs, nbins):
	"""-GC at every ceil(n / nbins)-th target in order of -GC; the last bin is
	open."""

	stats = sorted(-g for g in target_gcs)
	per = int(math.ceil(len(stats) / float(nbins)))
	thresholds = [stats[min((i + 1) * per, len(stats) - 1)] for i in range(nbins)]
	thresholds[-1] = float('inf')
	return thresholds


def bin_of(gc, thresholds):
	for i, t in enumerate(thresholds):
		if -gc <= t:
			return i
	raise AssertionError(gc)


def ks_statistic(a, b):
	a, b = sorted(a), sorted(b)
	return max(abs(bisect.bisect_right(a, v) / float(len(a)) - bisect.bisect_right(b, v) /
		float(len(b))) for v in set(a) | set(b))


def kmer_names(k, alphabet):
	"""Column order of the frequency files: the first base of a k-mer is its
	least significant digit."""

	n = len(alphabet)
	return [''.join(alphabet[(i // n ** d) % n] for d in range(k)) for i in range(n ** k)]


def kmer_frequencies(seqs, k, alphabet, positional, length):
	"""Frequencies of each k-mer over all windows (or per start position),
	skipping windows with a letter outside the alphabet."""

	names = kmer_names(k, alphabet)
	rows = length if positional else 1
	counts = [collections.Counter() for _ in range(rows)]
	for seq in seqs:
		for j in range(len(seq) - k + 1):
			word = seq[j:j + k]
			if all(b in alphabet for b in word):
				counts[j if positional else 0][word] += 1
	freqs = []
	for c in counts:
		total = float(sum(c.values()))
		freqs.append([c[w] / total if total else 0.0 for w in names])
	return names, freqs


def read_freq_file(path):
	header, rows = hh.read_tsv(path, header=True)
	return header, [(r[0], [float(v) for v in r[1:]]) for r in rows]


# --- invariants of every selection from a genome ----------------------------------------------

def check_genome_selection(runs, key, mask=False, max_n=0.5, size=WINDOW):
	rows = read_selection(runs, key)
	assert rows
	for row in rows:
		assert row['end'] - row['start'] == size == len(row['seq'])
		assert row['seq'] == genome_sequence(row['chr'], row['start'], row['end'],
			row['strand'], mask), row
		assert row['seq'].count('N') <= max_n * size
		assert row['gc'] == pytest.approx(gc_fraction(row['seq']), abs=1e-6)
	return rows


GENOME_RUNS = ['peaks', 'seed2', 'peaks_bed', 'default_prefix', 'fasta_targets', 'screen_small',
	'big', 'big_allow_bg', 'big_strand', 'big_overlap_late', 'overlap_first', 'allow_target',
	'nrun_default', 'nbins5', 'nsub5', 'pkmer2', 'ikmer3', 'include_ns', 'exclude_ns', 'decay',
	'bg_positions', 'bg_regions', 'bg_regions_allow']


@pytest.mark.parametrize('key', GENOME_RUNS)
def test_selected_sequences_are_the_genome(runs, key):
	"""Each selected sequence is the genome at its reported 0-based start and
	end, reverse-complemented on the - strand, with at most half of it N (the
	default -maxFractionN) and the reported GC fraction."""

	check_genome_selection(runs, key)


@pytest.mark.parametrize('key', ['peaks', 'seed2', 'peaks_bed', 'default_prefix', 'screen_small',
	'big', 'big_allow_bg', 'big_overlap_late', 'overlap_first', 'mask', 'no_n', 'nbins5', 'nsub5',
	'pkmer2', 'ikmer3', 'include_ns', 'exclude_ns', 'decay', 'bg_regions', 'bg_regions_allow'])
def test_no_selected_window_overlaps_a_target(runs, key):
	"""Windows overlapping the -p targets are not selected unless
	-allowTargetOverlaps is given."""

	assert overlaps_with(read_selection(runs, key), peak_intervals()) == []


@pytest.mark.parametrize('key', ['peaks', 'seed2', 'default_prefix', 'fasta_targets',
	'screen_small', 'big', 'overlap_first', 'allow_target', 'mask', 'no_n', 'nbins5', 'nsub5',
	'pkmer2', 'ikmer3', 'decay', 'bg_regions'])
def test_selected_windows_do_not_overlap(runs, key):
	rows = read_selection(runs, key)
	assert overlapping_pairs([(r['chr'], r['start'], r['end'], r['strand']) for r in rows]) == []


@pytest.mark.parametrize('key,nbins', [('peaks', 3), ('seed2', 3), ('fasta_targets', 3),
	('nbins5', 5), ('nsub5', 10), ('include_ns', 3), ('decay', 3), ('bg_fasta', 3),
	('bg_positions', 3), ('mask', 3)])
def test_gc_bins(runs, key, nbins):
	"""Targets and selected sequences are binned by GC at the targets'
	1/nbins quantiles; each bin holds the same number of selected sequences.
	Without -nbins or -nsubBins, 200 targets give 3 bins; setting -nsubBins
	alone keeps the documented -nbins default of 10."""

	targets = read_selection(runs, key, 'target')
	thresholds = bin_thresholds([gc_fraction(t['seq']) for t in targets], nbins)
	for t in targets:
		assert t['bin'] == bin_of(gc_fraction(t['seq']), thresholds)
	rows = read_selection(runs, key)
	for row in rows:
		assert row['bin'] == bin_of(gc_fraction(row['seq']), thresholds), row
	per_bin = collections.Counter(row['bin'] for row in rows)
	assert sorted(per_bin) == list(range(nbins))
	assert len(set(per_bin.values())) == 1


@pytest.mark.parametrize('key', ['peaks', 'fasta_targets', 'mask', 'no_n', 'nbins5', 'ikmer3',
	'decay', 'bg_fasta', 'bg_positions', 'big_strand'])
def test_gc_distribution_tracks_the_targets(runs, key):
	"""The selection is drawn bin by bin to match the targets' GC, so a
	two-sample Kolmogorov-Smirnov test does not tell the two GC distributions
	apart at alpha = 0.001 (critical D = 1.949 * sqrt((n + m) / (n m)))."""

	targets = [gc_fraction(t['seq']) for t in read_selection(runs, key, 'target')]
	selected = [gc_fraction(r['seq']) for r in read_selection(runs, key)]
	n, m = len(targets), len(selected)
	assert ks_statistic(targets, selected) < 1.949 * math.sqrt((n + m) / float(n * m))


# --- targets ---------------------------------------------------------------------------------

@pytest.mark.parametrize('key', ['peaks', 'peaks_bed'])
def test_targets_from_positions(runs, key):
	"""-p reads a HOMER peak file (1-based) or a BED file (0-based); either
	way the targets are the genome under each peak, on the peak's strand."""

	targets = read_selection(runs, key, 'target')
	got = sorted((t['chr'], t['start'], t['end'], t['strand']) for t in targets)
	assert got == sorted(peak_intervals())
	for t in targets:
		assert t['seq'] == genome_sequence(t['chr'], t['start'], t['end'], t['strand'])


def test_bed_and_peak_files_select_the_same_background(runs):
	for name in ('bg.positions.bed', 'bg.sequences.fasta', 'target.positions.bed'):
		assert read_output(runs, 'peaks', name) == read_output(runs, 'peaks_bed', name), name


def test_targets_from_fasta(runs):
	"""-i targets are the FASTA records as given, upper-cased."""

	targets = read_selection(runs, 'fasta_targets', 'target')
	records = dict(hh.read_fasta(hh.TARGETS))
	assert sorted(t['chr'] for t in targets) == sorted(records)
	for t in targets:
		assert (t['start'], t['end'], t['strand']) == (0, len(records[t['chr']]), '+')
		assert t['seq'] == records[t['chr']].upper()


@pytest.mark.parametrize('key,max_n,kept', [('nrun_default', 0.5, ['nrun25']),
	('nrun_02', 0.2, []), ('nrun_095', 0.95, ['nrun25', 'nrun60', 'nrun90'])])
def test_max_fraction_n_filters_targets(runs, key, max_n, kept):
	"""Peaks with 25%, 60% and 90% N are kept only when their fraction of N is
	at most -maxFractionN; no selected background window holds more N than
	that either."""

	targets = read_selection(runs, key, 'target')
	nrun = sorted(name for name, k in NRUN_PEAKS.items() if any(
		t['chr'] == N_RUN[0] and t['start'] == N_RUN[1] - (WINDOW - k) for t in targets))
	assert nrun == kept
	assert len(targets) == 150 + len(kept)
	check_genome_selection(runs, key, max_n=max_n)


def test_max_fraction_n_zero_selects_no_n(runs):
	rows = check_genome_selection(runs, 'no_n', max_n=0.0)
	assert all('N' not in r['seq'] for r in rows)


def test_mask_turns_soft_masked_bases_into_n(runs):
	"""-mask reads lowercase genome bases as N: with -maxFractionN 0 the
	targets are the peaks without a lowercase base, and no selected window
	holds a soft-masked base."""

	genome = hh.read_genome()
	unmasked = [iv for iv in peak_intervals() if not any(b.islower()
		for b in genome[iv[0]][iv[1]:iv[2]])]
	targets = read_selection(runs, 'mask', 'target')
	assert sorted((t['chr'], t['start'], t['end'], t['strand']) for t in targets) == \
		sorted(unmasked)
	assert len(unmasked) < len(peak_intervals())
	for t in targets:
		assert t['seq'] == genome_sequence(t['chr'], t['start'], t['end'], t['strand'], mask=True)
	rows = check_genome_selection(runs, 'mask', mask=True, max_n=0.0)
	for r in rows:
		assert not any(b.islower() for b in genome[r['chr']][r['start']:r['end']])


def test_without_mask_soft_masked_bases_are_used(runs):
	"""Without -mask, lowercase genome bases are ordinary sequence: the ten
	lowercase peaks are targets, upper-cased."""

	targets = read_selection(runs, 'peaks', 'target')
	genome = hh.read_genome()
	lowercase = [t for t in targets if genome[t['chr']][t['start']:t['end']].islower()]
	assert len(targets) == 200 and len(lowercase) == 10


# --- counts and sizes -----------------------------------------------------------------------

@pytest.mark.parametrize('key,n', [('peaks', 300), ('default_prefix', 30), ('size100', 150),
	('nbins5', 300), ('mixed_size', 150)])
def test_n_sequences_are_selected(runs, key, n):
	"""-N sequences are selected when the genome has room for them (each -N
	here is a multiple of the number of bins)."""

	assert len(read_selection(runs, key)) == n


def option(args, name, default):
	return args[args.index(name) + 1] if name in args else default


@pytest.mark.parametrize('key', sorted(run_table({'nrun_peaks': 'x', 'bg_positions': 'x',
	'regions': 'x', 'mixed': 'x'})))
def test_at_most_n_sequences(runs, table, key):
	"""No run selects more than -N sequences, and every GC bin contributes
	the same number (each -N here is a multiple of the number of bins)."""

	rows = read_selection(runs, key)
	assert len(rows) <= int(option(table[key], '-N', 100000))
	if '-model' not in table[key]:
		per_bin = collections.Counter(r['bin'] for r in rows)
		assert len(set(per_bin.values())) <= 1


def test_size_sets_the_window_length(runs):
	rows = check_genome_selection(runs, 'size100', size=100)
	assert len(rows) == 150


def test_default_size_is_the_mean_target_length(runs):
	"""Targets of 200 and 180 bp in turn average 190 bp, the default window
	length."""

	rows = read_selection(runs, 'mixed_size')
	lengths = [len(seq) for _, seq in hh.read_fasta(out_path(runs, 'mixed_size',
		'target.sequences.fasta'))]
	assert sorted(set(lengths)) == [180, 200]
	assert int(sum(lengths) / float(len(lengths))) == 190
	for row in rows:
		assert len(row['seq']) == 190 == row['end'] - row['start']
		assert row['seq'] == genome_sequence(row['chr'], row['start'], row['end'], row['strand'])


def test_screen_size_limits_the_candidates(runs):
	"""-NN 1000 keeps about 1,000 of the genome's windows (both strands) as
	candidates, so at most 1,200 (more than six standard deviations above
	1,000) can be selected, fewer than from a 60,000-window screen."""

	small = read_selection(runs, 'screen_small')
	assert 0 < len(small) <= 1200
	assert len(small) < len(read_selection(runs, 'big'))


def test_non_overlapping_selection_fits_the_genome(runs):
	"""Asked for 3,000 sequences, a non-overlapping selection is limited by
	how many separate 200 bp windows the genome holds."""

	rows = read_selection(runs, 'big')
	assert len(rows) <= capacity(len(s) for s in hh.read_genome().values())


@pytest.mark.parametrize('key', ['big_allow_bg', 'big_overlap_late'])
def test_overlaps_allowed(runs, key):
	"""With -allowBgOverlaps, or with -overlapIteration past the last of the
	20 selection rounds, all 3,000 sequences are selected: more than fit without
	overlap, so selected windows overlap."""

	rows = read_selection(runs, key)
	assert len(rows) == 3000 > capacity(len(s) for s in hh.read_genome().values())
	assert overlapping_pairs([(r['chr'], r['start'], r['end'], r['strand']) for r in rows])


def test_strand_allows_overlaps_on_opposite_strands(runs):
	"""-strand lets windows overlap other windows and targets on the other
	strand, never on the same strand."""

	rows = read_selection(runs, 'big_strand')
	windows = [(r['chr'], r['start'], r['end'], r['strand']) for r in rows]
	pairs = overlapping_pairs(windows)
	target_hits = overlaps_with(rows, peak_intervals())
	assert all(a[3] != b[3] for a, b in pairs)
	assert all(row['strand'] != iv[3] for row, iv in target_hits)
	assert pairs or target_hits
	assert len(rows) > len(read_selection(runs, 'big'))


def test_allow_target_overlaps(runs):
	"""Peaks and their flanks cover about a third of the genome's window
	starts, so 300 windows chosen without the target mask overlap peaks."""

	assert overlaps_with(read_selection(runs, 'allow_target'), peak_intervals())


# --- background sources ----------------------------------------------------------------------

def test_background_positions_are_taken_as_given(runs, inputs):
	"""-bg: every selected window is one of the given positions, on its given
	strand."""

	given = {(r[0], int(r[1]), int(r[2]), r[5]) for r in hh.read_tsv(inputs['bg_positions'])}
	rows = read_selection(runs, 'bg_positions')
	assert 0 < len(rows) <= 300 and len(rows) % 3 == 0
	for r in rows:
		assert (r['chr'], r['start'], r['end'], r['strand']) in given


def test_background_ranges(runs):
	"""-bgr: selected windows lie inside the ranges; without overlaps no more
	fit than separate windows in the ranges allow, and with -allowBgOverlaps
	all 300 are selected."""

	cap = capacity(end - start for _, start, end in REGIONS)
	for key in ('bg_regions', 'bg_regions_allow'):
		for r in read_selection(runs, key):
			assert any(r['chr'] == c and s <= r['start'] and r['end'] <= e for c, s, e in REGIONS)
	assert 0 < len(read_selection(runs, 'bg_regions')) <= cap
	assert len(read_selection(runs, 'bg_regions_allow')) == 300 > cap


def test_background_fasta(runs):
	"""-b: selected sequences are records of the FASTA file, whole and
	upper-cased."""

	records = dict(hh.read_fasta(hh.BACKGROUND))
	rows = read_selection(runs, 'bg_fasta')
	assert len(rows) == 300
	for r in rows:
		assert (r['start'], r['end'], r['strand']) == (0, len(records[r['chr']]), '+')
		assert r['seq'] == records[r['chr']].upper()


def test_profile_scores(runs):
	"""-pscore writes a table of every target, then a table of every
	candidate, with its profile score: the sum over its bases of (background -
	target) base frequency within its GC bin, frequencies taken over all
	targets and all candidates of the bin. The stats file reports the same
	scores for the selection."""

	cwd, _ = runs['bg_fasta']
	header = ['SeqID', 'chr', 'position', 'strand', 'weight', 'ProfileScore', 'GC', 'bin']
	lines = hh.read_tsv(str(cwd / 'pscore.txt'))
	assert lines[0] == header and lines[201] == header
	rows = lines[1:201] + lines[202:]
	targets = [(n, s.upper()) for n, s in hh.read_fasta(hh.TARGETS)]
	background = [(n, s.upper()) for n, s in hh.read_fasta(hh.BACKGROUND)]
	thresholds = bin_thresholds([gc_fraction(s) for _, s in targets], 3)
	bins = {}
	for group, recs in (('t', targets), ('b', background)):
		for name, seq in recs:
			bins[(group, name)] = bin_of(gc_fraction(seq), thresholds)
	freqs = {}
	for group, recs in (('t', targets), ('b', background)):
		for b in range(3):
			_, f = kmer_frequencies([s for n, s in recs if bins[(group, n)] == b], 1, 'ACGT',
				False, 1)
			freqs[(group, b)] = dict(zip('ACGT', f[0]))

	def score(seq, b):
		return sum(freqs[('b', b)][base] - freqs[('t', b)][base] for base in seq)

	assert len(rows) == len(targets) + len(background)
	expected = {}
	for group, recs in (('t', targets), ('b', background)):
		for name, seq in recs:
			expected[(group, name)] = (score(seq, bins[(group, name)]), gc_fraction(seq),
				bins[(group, name)])
	got = {}
	for i, row in enumerate(rows):
		group = 't' if i < len(targets) else 'b'
		got[(group, row[1])] = (float(row[5]), float(row[6]), int(row[7]))
	assert sorted(got) == sorted(expected)
	for k, (ps, gc, b) in expected.items():
		assert got[k][0] == pytest.approx(ps, abs=1e-5), k
		assert got[k][1] == pytest.approx(gc, abs=1e-6), k
		assert got[k][2] == b, k
	for r in read_selection(runs, 'bg_fasta'):
		assert r['pscore'] == pytest.approx(expected[('b', r['chr'])][0], abs=1e-5)


def test_model_generates_sequences(runs):
	"""-model writes -N generated sequences of the target length, named
	Random, made of the targets' letters; the first, middle and last third come
	from the high, middle and low GC bins, and each third's mean GC is within
	0.02 of its target bin's (a 100-sequence mean of 200 bp sequences varies
	by about 0.004)."""

	rows = read_selection(runs, 'model')
	assert len(rows) == 300
	for r in rows:
		assert r['chr'] == 'Random'
		assert len(r['seq']) == 200 == r['end'] - r['start']
		assert set(r['seq']) <= set('ACGT')
	targets = [gc_fraction(t['seq']) for t in read_selection(runs, 'model', 'target')]
	thresholds = bin_thresholds(targets, 3)
	for b in range(3):
		target_mean = sum(g for g in targets if bin_of(g, thresholds) == b) / \
			float(sum(1 for g in targets if bin_of(g, thresholds) == b))
		block = rows[100 * b:100 * (b + 1)]
		mean = sum(gc_fraction(r['seq']) for r in block) / 100.0
		assert mean == pytest.approx(target_mean, abs=0.02), b


# --- frequency, info and homer-style outputs ---------------------------------------------------

@pytest.mark.parametrize('key,k,alphabet,positional', [('peaks', 1, 'ACGT', False),
	('exclude_ns', 1, 'ACGT', False), ('include_ns', 1, 'ACGTN', False),
	('ikmer3', 3, 'ACGT', False), ('pkmer2', 2, 'ACGT', True), ('model', 1, 'ACGTN', False),
	('model_exclude_ns', 1, 'ACGT', False), ('bg_fasta', 1, 'ACGT', False)])
def test_frequency_files(runs, key, k, alphabet, positional):
	"""bg.freq.txt and target.freq.txt hold the k-mer frequencies of the
	selected background and of the targets, overall (-ikmer) or per start
	position (-pkmer); N is a letter only with -includeNs, or with -model
	unless -excludeNs; info.txt reports the count and the summed absolute
	frequency difference."""

	total_diff = 0.0
	tables = {}
	for which in ('bg', 'target'):
		seqs = [s for _, s in hh.read_fasta(out_path(runs, key, which + '.sequences.fasta'))]
		length = len(read_selection(runs, key)[0]['seq'])
		names, want = kmer_frequencies(seqs, k, alphabet, positional, length)
		header, rows = read_freq_file(out_path(runs, key, which + '.freq.txt'))
		assert header == ['Position'] + names
		if positional:
			assert [r[0] for r in rows] == [str(i) for i in range(length)]
		else:
			assert [r[0] for r in rows] == ['PositionIndependent']
		for (_, got), exp in zip(rows, want):
			assert got == pytest.approx(exp, abs=1e-6)
		tables[which] = [v for _, row in rows for v in row]
	total_diff = sum(abs(a - b) for a, b in zip(tables['bg'], tables['target']))
	info = read_output(runs, key, 'info.txt').splitlines()
	assert info[0].startswith('cmd = ')
	assert info[1] == 'total bg sequences selected = {}'.format(len(read_selection(runs, key)))
	assert info[2].startswith('frequency difference = ')
	assert float(info[2].split(' = ')[1]) == pytest.approx(total_diff, abs=len(tables['bg']) * 2e-6)


def test_homer_style_output(runs):
	"""-gs also writes a group file (targets 1, background 0, weight 1) and a
	sequence file with the same sequences as the FASTA outputs."""

	group = hh.read_tsv(out_path(runs, 'peaks', 'group.tsv'))
	seqs = hh.read_tsv(out_path(runs, 'peaks', 'seq.tsv'))
	assert [g[0] for g in group] == [r[0] for r in seqs]
	assert len(set(g[0] for g in group)) == len(group)
	assert all(float(g[2]) == 1.0 for g in group)

	def site(name):
		# names are <chr>_<position>_<strand>_<serial number>
		return name.split(' ')[0].rsplit('_', 1)[0]

	want = {}
	for label, which in (('1', 'target'), ('0', 'bg')):
		for name, seq in hh.read_fasta(out_path(runs, 'peaks', which + '.sequences.fasta')):
			want[site(name)] = (label, seq)
	got = {site(g[0]): (g[1], r[1]) for g, r in zip(group, seqs)}
	assert got == want
	assert not os.path.exists(out_path(runs, 'seed2', 'group.tsv'))


def test_homer_style_output_feeds_denovo(runs, tmp_path):
	"""The -gs files are homer2 denovo input: the planted siteA comes out as
	the top motif (consensus within one mismatch of the site or its reverse
	complement)."""

	iupac = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG',
		'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
		'N': 'ACGT'}
	hh.run(['homer2', 'denovo', '-s', out_path(runs, 'peaks', 'seq.tsv'), '-g',
		out_path(runs, 'peaks', 'group.tsv'), '-len', 8, '-S', 1, '-o', 'out.motifs'], tmp_path)
	top = hh.read_motifs(str(tmp_path / 'out.motifs'))[0]
	site = 'TAGCGTACCGTA'
	best = min(sum(1 for c, b in zip(top.consensus, s[i:i + 8]) if b not in iupac[c])
		for s in (site, hh.revcomp(site)) for i in range(len(site) - 8 + 1))
	assert best <= 1, top.consensus


# --- seed, prefix and alias ------------------------------------------------------------------

def test_same_seed_same_output(runs):
	for name in OUTPUTS + ['group.tsv', 'seq.tsv']:
		assert read_output(runs, 'peaks', name) == read_output(runs, 'seed_again', name), name


def test_other_seed_other_selection(runs):
	a = {(r['chr'], r['start'], r['strand']) for r in read_selection(runs, 'peaks')}
	b = {(r['chr'], r['start'], r['strand']) for r in read_selection(runs, 'seed2')}
	assert a != b


def test_background_is_an_alias_of_bg(runs):
	"""Both names give the same files, apart from info.txt, which records the
	command line."""

	for name in OUTPUTS[:-1] + ['group.tsv', 'seq.tsv']:
		assert read_output(runs, 'peaks', name) == read_output(runs, 'alias', name), name


def test_output_prefix(runs):
	"""Without -o every file is named out.<kind>; -o sets the prefix."""

	cwd, _ = runs['default_prefix']
	assert sorted(os.listdir(str(cwd))) == sorted('out.' + name for name in OUTPUTS)
	cwd, _ = runs['seed2']
	assert sorted(os.listdir(str(cwd))) == sorted('sel.' + name for name in OUTPUTS)


@pytest.mark.parametrize('key', ['decay', 'overlap_first', 'nsub5'])
def test_selection_settings_change_the_selection(runs, key):
	"""-decayRate, -overlapIteration and -nsubBins change which windows are
	picked with the same seed (the invariants above hold for each)."""

	a = {(r['chr'], r['start'], r['strand']) for r in read_selection(runs, 'peaks')}
	b = {(r['chr'], r['start'], r['strand']) for r in read_selection(runs, key)}
	assert a != b

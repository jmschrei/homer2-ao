"""findMotifsGenome.pl options that choose and weight the sequences, with a
custom genome: -size, -mask, -maxN, -bg, -chopify, -local, -ref, -N, -gc,
-cpg, -noweight, -nlen, -nmax, -neutral, -redundant, -useOldBg, -preparse,
-preparsedDir, -useNewBg with its options, BED input and a genome directory.
Each is checked against a run without it or against a reference computation.

All runs start in one batch (module fixture `runs`). Expected values come from
the inputs: windows cut from tests/data/genome.fa as resizePosFile.pl places
them, sequences dropped as removePoorSeq.pl and assignGeneWeights.pl drop them
(more than -maxN N; GC bins lacking targets or background), GC and CpG bins
with freq2group.pl's edges, oligo frequencies counted as homer2 norm counts
them, sites from homer_helpers.scan. Options passed through to homer2 bg are
checked by their visible effect on its output files and by the command line
findMotifsGenome.pl records in motifFindingParameters.txt; homer2 bg draws its
background with a clock seed, so those tests assert only what holds for any
draw.
"""

import os
import shutil

import pytest

import homer_helpers as hh
import findmotifsgenome_helpers as fh

KS = fh.KNOWN_SMALL
KNOWN = ['-size', 200, '-nomotif', '-mknown', KS, '-N', 1000]
KEEP = ['-dumpFasta', '-keepFiles']
NEWBG = ['-size', 200, '-nomotif', '-mknown', KS, '-useNewBg', '-N', 1000, '-NN', 20000] + KEEP
PREPARSED_FILES = ['cgbins', 'cgfreq', 'gcbins', 'pos', 'seq']


def write_peaks(path, rows):
	with open(str(path), 'w') as handle:
		for row in rows:
			handle.write('\t'.join(str(v) for v in row) + '\n')
	return str(path)


def gaps(genome, peaks, margin=500):
	"""Stretches at least `margin` bp from every peak, the chromosome ends and
	the N run."""

	result = []
	for chrom in sorted(genome):
		spans = sorted((p['start'], p['end']) for p in peaks if p['chr'] == chrom)
		edges = [(0, 1000)] + spans + [(len(genome[chrom]) - 1000, len(genome[chrom]))]
		for (_, left), (right, _) in zip(edges, edges[1:]):
			lo, hi = left + margin, right - margin
			if 'N' in genome[chrom][lo - margin:hi + margin].upper():
				continue
			if hi - lo >= 200:
				result.append((chrom, lo, hi))
	return result


@pytest.fixture(scope='module')
def inputs(tmp_path_factory, genome, peaks, preparsed):
	tmp = tmp_path_factory.mktemp('fmg_inputs')
	free = gaps(genome, peaks)
	paths = {
		'bg': write_peaks(tmp / 'bg.txt', [('bg{:03d}'.format(i), c, lo + 1, lo + 200, '+')
			for i, (c, lo, hi) in enumerate(free, 1)] + [('ov_' + p['id'], p['chr'],
			p['start'], p['end'], p['strand']) for p in peaks[:5]]),
		'big': write_peaks(tmp / 'big.txt', [('big{:03d}'.format(i), c, lo + 1, lo + 1000, '+')
			for i, (c, lo, hi) in enumerate([g for g in free if g[2] - g[1] >= 1000][:30], 1)]),
		# Windows ending 20, 100 and 160 bases into chr2's N run (40,001-42,000).
		'npeaks': write_peaks(tmp / 'npeaks.txt', [('N{}'.format(k), 'chr2', 40000 + k - 199,
			40000 + k, '+') for k in (20, 100, 160)] + [(p['id'], p['chr'], p['start'],
			p['end'], p['strand']) for p in peaks[:40]]),
		'ids': write_peaks(tmp / 'ids.txt', [(name,) for name in sorted({row['peak']
			for row in hh.read_planted() if row['motif'] == 'siteA'})]),
	}
	# A genome path ending in r would be read as the masked form of the name.
	chroms = tmp / 'chroms'
	chroms.mkdir()
	for name in sorted(genome):
		hh.write_fasta(str(chroms / (name + '.fa')), [(name, genome[name])])
	paths['chroms'] = str(chroms)
	hh.run_many({
		'pp199': (['preparseGenome.pl', hh.GENOME, '-size', 199, '-preparsedDir',
			tmp / 'pp199'], tmp / 'pp199'),
		'ppmask': (['preparseGenome.pl', hh.GENOME, '-mask', '-size', 200, '-preparsedDir',
			tmp / 'ppmask'], tmp / 'ppmask'),
	})
	stale = tmp / 'stale'
	shutil.copytree(str(preparsed), str(stale))
	for ext in ('pos', 'seq'):
		path = stale / 'genome.fa.200.{}'.format(ext)
		path.write_text(''.join(path.read_text().splitlines(True)[:10]))
	for name in ('pp199', 'ppmask', 'stale'):
		paths[name] = str(tmp / name)
	for name in ('pp150', 'ppdir', 'missing/sub', 'newbg'):
		paths[name] = str(tmp / name)
	return paths


@pytest.fixture(scope='module')
def runs(tmp_path_factory, inputs, preparsed):
	pp = {'preparsed': preparsed}
	newbg = {'preparsed': inputs['newbg']}
	calls = {
		'base': (KNOWN + KEEP, pp),
		'size_given': (['-size', 'given', '-nomotif', '-mknown', KS, '-N', 1000, '-dumpFasta'],
			{'preparsed': inputs['pp199']}),
		'size_offsets': (['-size', '-100,50', '-nomotif', '-mknown', KS, '-N', 1000,
			'-dumpFasta', '-preparse'], {'preparsed': inputs['pp150']}),
		'mask': (KNOWN + ['-mask', '-dumpFasta'], {'preparsed': inputs['ppmask']}),
		'mask_maxN': (KNOWN + ['-mask', '-maxN', 1, '-dumpFasta'], {'preparsed': inputs['ppmask']}),
		'maxN': (KNOWN + ['-dumpFasta'], dict(pp, peaks=inputs['npeaks'])),
		'maxN_03': (KNOWN + ['-maxN', 0.3, '-dumpFasta'], dict(pp, peaks=inputs['npeaks'])),
		'bg': (KNOWN + KEEP + ['-bg', inputs['bg']], pp),
		'chopify': (KNOWN + KEEP + ['-bg', inputs['big'], '-chopify'], pp),
		'local': (KNOWN + KEEP + ['-local', 1], pp),
		'ref': (KNOWN + KEEP + ['-ref', hh.PEAKS], dict(pp, peaks=inputs['ids'])),
		'N_low': (['-size', 200, '-nomotif', '-mknown', KS, '-N', 100, '-dumpFasta'], pp),
		'nlen2': (KNOWN + KEEP + ['-nlen', 2], pp),
		'nlen0': (['-size', 200, '-nomotif', '-mknown', KS, '-nlen', 0] + KEEP, pp),
		'nmax1': (KNOWN + KEEP + ['-nmax', 1], pp),
		'neutral': (KNOWN + KEEP + ['-neutral'], pp),
		'cpg': (KNOWN + KEEP + ['-cpg'], pp),
		'noweight': (KNOWN + KEEP + ['-noweight', '-nlen', 0], pp),
		'redundant': (KNOWN + ['-redundant', 2], pp),
		'old_bg': (KNOWN + ['-useOldBg'], pp),
		'gc': (KNOWN + ['-gc'], pp),
		'bed': (KNOWN, dict(pp, peaks=hh.PEAKS_BED)),
		'genome_dir': (KNOWN + ['-dumpFasta', '-preparse'],
			{'genome': inputs['chroms'], 'preparsed': inputs['ppdir']}),
		'preparse': (KNOWN + ['-preparse'], {'preparsed': inputs['stale']}),
		'preparsed_missing': (KNOWN, {'preparsed': inputs['missing/sub']}),
		'newbg': (NEWBG, newbg),
		'newbg_model': (NEWBG + ['-modelBg'], newbg),
		'newbg_model_genome': (NEWBG + ['-modelBg', '-genomeBg'], newbg),
		'newbg_bins': (NEWBG + ['-numBins', 3], newbg),
		'newbg_ikmer': (NEWBG + ['-ikmer', 3], newbg),
		'newbg_pkmer': (NEWBG + ['-pkmer', 2], newbg),
		'newbg_target_overlaps': (NEWBG + ['-allowTargetOverlaps'], newbg),
		'newbg_bg_overlaps': (NEWBG + ['-allowBgOverlaps'], newbg),
		'newbg_genome_dir': (NEWBG, dict(newbg, genome=inputs['chroms'])),
	}
	return fh.run_calls(tmp_path_factory.mktemp('fmg_options'), calls)


@pytest.fixture(scope='module')
def library():
	return hh.read_motifs(KS)


@pytest.fixture(scope='module')
def windows(genome, peaks):
	return {p['id']: fh.window_sequence(p, 200, genome=genome).upper() for p in peaks}


def bin_counts(seqs, bin_of=fh.gc_bin):
	counts = {}
	for seq in seqs:
		counts[bin_of(seq)] = counts.get(bin_of(seq), 0) + 1
	return counts


# --- -size, -mask, -maxN ----------------------------------------------------------------

def test_size_given_uses_the_regions_as_they_are(runs, inputs, peaks, genome, library):
	"""-size given extracts each peak's own bases and takes the background
	size from the mean width, end - start = 199, so the preparsed 199 bp files
	are used and no other size is made."""

	targets = fh.out_fasta(runs, 'size_given', 'target.fa')
	assert targets == {p['id']: hh.peak_sequence(p, genome).upper() for p in peaks}
	assert sorted(os.listdir(inputs['pp199'])) == ['genome.fa.199.' + e for e in PREPARSED_FILES]
	rows = fh.known_rows(runs, 'size_given')
	for motif in library:
		assert rows[motif.name]['targets'] == sum(1 for s in targets.values()
			if fh.has_site(s, motif))


def test_size_with_two_offsets_is_relative_to_the_peak_centre(runs, inputs, peaks, genome):
	"""-size -100,50 takes -100 to +50 around the centre along the peak's
	strand: a 150 bp window moved by floor((-100 + 50) / 2) = -25. -preparse
	makes the 150 bp background files in the empty -preparsedDir."""

	want = {p['id']: fh.window_sequence(p, 150, -25, genome).upper() for p in peaks}
	assert fh.out_fasta(runs, 'size_offsets', 'target.fa') == want
	assert sorted(os.listdir(inputs['pp150'])) == ['genome.fa.150.' + e for e in PREPARSED_FILES]


def test_mask_turns_lowercase_into_N_and_drops_masked_peaks(runs, inputs, peaks, genome):
	"""-mask reads soft-masked bases as N; the ten all-lowercase siteA peaks
	are then over 70% N and dropped, leaving 110 siteA targets of 190. The
	background comes from the 'r' (masked) preparsed files."""

	want = {}
	for p in peaks:
		seq = fh.masked(fh.window_sequence(p, 200, genome=genome))
		if fh.kept_by_max_n(seq):
			want[p['id']] = seq
	assert len(want) == 190
	assert fh.out_fasta(runs, 'mask', 'target.fa') == want
	rows = fh.known_rows(runs, 'mask')
	assert (rows['siteA']['targets'], rows['siteA']['target_total']) == (110, 190)
	assert sorted(os.listdir(inputs['ppmask'])) == ['genome.far.200.' + e for e in PREPARSED_FILES]


def test_mask_with_max_N_1_keeps_masked_peaks_until_gc_binning(runs, peaks, genome):
	"""With -maxN 1 the masked peaks stay in target.fa, all but one base N and
	without their siteA site; with one base left their GC fraction is 0 or 1,
	bins that hold no background, so assignGeneWeights.pl leaves them out of
	the analysis."""

	lowercase = {row['peak'] for row in hh.read_planted() if row['lowercase'] == 'yes'}
	assert len(lowercase) == 10
	targets = fh.out_fasta(runs, 'mask_maxN', 'target.fa')
	assert targets == {p['id']: fh.masked(fh.window_sequence(p, 200, genome=genome))
		for p in peaks}
	site = hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA.motif'))[0]
	for name in lowercase:
		assert targets[name].count('N') == 200 and not fh.has_site(targets[name], site)
	background_bins = set(bin_counts(fh.out_fasta(runs, 'mask_maxN', 'background.fa').values()))
	analysed = [name for name, seq in targets.items() if fh.gc_bin(seq) in background_bins]
	assert not lowercase & set(analysed)
	rows = fh.known_rows(runs, 'mask_maxN')
	assert (rows['siteA']['targets'], rows['siteA']['target_total']) == (110, len(analysed))


@pytest.mark.parametrize('key,max_n', [('maxN', 0.7), ('maxN_03', 0.3)])
def test_max_N_drops_sequences_with_too_many_N(runs, inputs, genome, key, max_n):
	"""Windows with 20, 100 and 160 of their 201 bases in the N run: those
	over the -maxN fraction (default 0.7) of N are dropped."""

	want = [p['id'] for p in hh.read_peaks(inputs['npeaks'])
		if fh.kept_by_max_n(fh.window_sequence(p, 200, genome=genome), max_n)]
	assert sorted(fh.out_fasta(runs, key, 'target.fa')) == sorted(want)
	assert len(want) == 40 + (2 if max_n == 0.7 else 1)


# --- background choice ------------------------------------------------------------------

def test_background_count_per_gc_bin(runs):
	"""The background matches the targets' GC bins in proportion. The usage
	calls -N the number of background sequences, and the script's own message
	the total including targets, so either reading is accepted: with 200
	targets, -N 1000 gives 5 or 4 background sequences per target in each bin,
	and -N 100, raised to twice the targets, 2 or 1."""

	for key, n in (('base', 1000), ('N_low', 400)):
		targets = bin_counts(fh.out_fasta(runs, key, 'target.fa').values())
		background = bin_counts(fh.out_fasta(runs, key, 'background.fa').values())
		total = sum(targets.values())
		assert background in [{b: ratio * c for b, c in targets.items()}
			for ratio in (n // total, (n - total) // total)], key


def test_gc_bins_and_cpg_bins(runs):
	"""Targets are binned by GC fraction by default (-gc) and by CpG ratio
	with -cpg, with freq2group.pl's edges."""

	for key, bin_of in (('base', fh.gc_bin), ('cpg', fh.cpg_bin)):
		bins = {row[0]: int(row[1]) for row in hh.read_tsv(fh.out(runs, key, 'target200.cgbins'))}
		targets = fh.out_fasta(runs, key, 'target.fa')
		assert bins == {name: bin_of(seq) for name, seq in targets.items()}, key


def test_gc_weights_without_oligo_normalization(runs):
	"""With -nlen 0 the weights are assignGeneWeights.pl's: equal within a GC
	bin, and each bin's background weight is its share of the targets."""

	assert not os.path.exists(fh.out(runs, 'nlen0', 'seq.autonorm.tsv'))
	groups = fh.read_groups(fh.out(runs, 'nlen0', 'group.adj'))
	targets = bin_counts(fh.out_fasta(runs, 'nlen0', 'target.fa').values())
	weights = {}
	for name, seq in fh.out_fasta(runs, 'nlen0', 'background.fa').items():
		weights.setdefault(fh.gc_bin(seq), []).append(groups[name][1])
	total = sum(sum(w) for w in weights.values())
	assert sorted(weights) == sorted(targets)
	for b, w in weights.items():
		assert max(w) - min(w) < 1e-9
		assert sum(w) / total == pytest.approx(targets[b] / 200.0, abs=1e-4)


def test_noweight_leaves_background_unweighted(runs, library):
	"""-noweight (no GC correction) with -nlen 0: the group file carries no
	weights, so known motif background counts are plain counts. -N 1000 gives
	1,000 or 800 background sequences (see test_background_count_per_gc_bin)."""

	groups = [row for row in hh.read_tsv(fh.out(runs, 'noweight', 'group.adj'))]
	assert {len(row) for row in groups} == {2}
	assert len(fh.out_fasta(runs, 'noweight', 'background.fa')) in (1000, 800)
	rows = fh.check_known_run(runs, 'noweight', library)
	assert all(row['background'] == int(row['background']) for row in rows)


def test_background_peak_file(runs, inputs, windows, genome, library):
	"""-bg: background windows from the given peaks, prefixed BG; those
	overlapping a target are removed (the five copies of target peaks), and
	GC bins without targets or background are dropped."""

	background = {'BG' + p['id']: fh.window_sequence(p, 200, genome=genome).upper()
		for p in hh.read_peaks(inputs['bg']) if not p['id'].startswith('ov_')}
	targets, background = fh.gc_matched(windows, background)
	assert fh.out_fasta(runs, 'bg', 'target.fa') == targets
	assert fh.out_fasta(runs, 'bg', 'background.fa') == background
	fh.check_known_run(runs, 'bg', library)


def test_chopify_cuts_background_regions_to_the_target_size(runs, inputs, windows, genome):
	"""-chopify cuts each background region into pieces of the mean target
	window size (200): <id>-chopped1, -chopped2, ... starting at the region's
	start, as many as size / 200 rounds to."""

	background = {}
	for p in hh.read_peaks(inputs['big']):
		count = 0
		while count < (p['end'] - p['start']) / 200.0 - 0.5:
			start = p['start'] + count * 200
			piece = {'chr': p['chr'], 'start': start, 'end': start + 200, 'strand': '+'}
			count += 1
			background['BG{}-chopped{}'.format(p['id'], count)] = \
				fh.window_sequence(piece, 200, genome=genome).upper()
	assert len(background) == 5 * len(hh.read_peaks(inputs['big'])) >= 100
	targets, background = fh.gc_matched(windows, background)
	assert fh.out_fasta(runs, 'chopify', 'target.fa') == targets
	assert fh.out_fasta(runs, 'chopify', 'background.fa') == background


def test_local_background_flanks_each_peak(runs, peaks, windows, genome, library):
	"""-local 1: the windows 200 bp to either side of each target window."""

	background = {}
	for p in peaks:
		start, end = fh.window(p, 200)
		for name, shift in (('-local+1', 200), ('-local-1', -200)):
			background[p['id'] + name] = fh.region_sequence(p['chr'], start + shift,
				end + shift, p['strand'], genome).upper()
	targets, background = fh.gc_matched(windows, background)
	assert fh.out_fasta(runs, 'local', 'target.fa') == targets
	assert fh.out_fasta(runs, 'local', 'background.fa') == background
	fh.check_known_run(runs, 'local', library)


def test_reference_peaks_split_into_targets_and_background(runs, inputs, windows, library):
	"""-ref <peak file> with a list of IDs as the first argument: the listed
	peaks are the targets and the rest of the file the background."""

	ids = {row[0] for row in hh.read_tsv(inputs['ids'])}
	targets, background = fh.gc_matched({n: s for n, s in windows.items() if n in ids},
		{n: s for n, s in windows.items() if n not in ids})
	assert fh.out_fasta(runs, 'ref', 'target.fa') == targets
	assert fh.out_fasta(runs, 'ref', 'background.fa') == background
	fh.check_known_run(runs, 'ref', library)


@pytest.mark.parametrize('key', ['redundant', 'old_bg', 'gc', 'bed'])
def test_runs_equal_to_the_default(runs, key):
	"""-redundant 2 (keep every sequence), -useOldBg and -gc (documented
	defaults) and the same peaks as BED give the default results."""

	assert fh.same_file(runs, key, 'base', 'knownResults.txt')


def test_genome_as_a_directory_of_fasta_files(runs, inputs):
	"""A directory of per-chromosome FASTA files gives the same targets; its
	preparsed files are named after the directory."""

	assert fh.out_fasta(runs, 'genome_dir', 'target.fa') == fh.out_fasta(runs, 'base', 'target.fa')
	a, b = fh.known_rows(runs, 'genome_dir'), fh.known_rows(runs, 'base')
	assert {n: r['targets'] for n, r in a.items()} == {n: r['targets'] for n, r in b.items()}
	assert sorted(os.listdir(inputs['ppdir'])) == ['chroms.200.' + e for e in PREPARSED_FILES]


@pytest.mark.parametrize('key,directory', [('preparse', 'stale'),
	('preparsed_missing', 'missing/sub')])
def test_preparsed_files_are_made_when_forced_or_missing(runs, inputs, preparsed, key, directory):
	"""-preparse rewrites damaged preparsed files, and a missing -preparsedDir
	is created and filled; either way the files equal a fresh preparse."""

	assert runs[key].returncode == 0
	for ext in PREPARSED_FILES:
		name = 'genome.fa.200.' + ext
		assert open(os.path.join(inputs[directory], name)).read() == \
			open(os.path.join(str(preparsed), name)).read(), name


# --- oligo normalization ------------------------------------------------------------------

def test_nlen_sets_the_oligo_lengths_normalized(runs):
	header, rows = hh.read_tsv(fh.out(runs, 'nlen2', 'seq.autonorm.tsv'), header=True)
	assert sorted(len(row[0]) for row in rows) == [1] * 4 + [2] * 16
	targets = [(s, 1.0) for s in fh.out_fasta(runs, 'nlen2', 'target.fa').values()]
	texts = fh.autonorm_target_texts(targets, (1, 2))
	for row in rows:
		assert row[1] in texts[row[0]], row


def weighted_sets(results, key):
	seqs = fh.read_seq_tsv(fh.out(results, key, 'seq.tsv'))
	groups = fh.read_groups(fh.out(results, key, 'group.adj'))
	return ([(seqs[n], w) for n, (g, w) in groups.items() if g == 1],
		[(seqs[n], w) for n, (g, w) in groups.items() if g == 0])


def test_more_normalization_matches_oligos_more_closely(runs):
	"""The weighted background's 3-mer content is closer to the targets'
	after the default -nlen 3 (up to -nmax 160 rounds) than after -nmax 1,
	and -nlen 2 matches 2-mers better than a run without normalization
	(-noweight -nlen 0)."""

	distance = {key: [fh.oligo_distance(*weighted_sets(runs, key), k=k) for k in (2, 3)]
		for key in ('base', 'nmax1', 'nlen2', 'noweight')}
	assert distance['base'][1] * 1.2 < distance['nmax1'][1]
	assert distance['nlen2'][0] * 1.2 < distance['noweight'][0]


def test_neutral_weights_toward_equal_base_frequencies(runs, library):
	targets, background = weighted_sets(runs, 'neutral')
	assert {w for _, w in targets} != {1.0}
	for seqs in (targets, background):
		for base, freq in fh.oligo_levels(seqs, 1).items():
			assert freq == pytest.approx(0.25, abs=0.002), base
	plain = fh.oligo_levels(weighted_sets(runs, 'base')[0], 1)
	assert max(abs(f - 0.25) for f in plain.values()) > 0.004
	fh.check_known_run(runs, 'neutral', library)


# --- the new background selection (homer2 bg) ------------------------------------------------

def bed(path):
	return [(row[0], int(row[1]), int(row[2])) for row in hh.read_tsv(path)]


def overlaps(a, b):
	return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]


def bg_command(results, key):
	return open(fh.out(results, key, 'motifFindingParameters.txt')).read().splitlines()[-1]


def in_genome(seqs, genome):
	text = [g.upper() for g in genome.values()]
	return sum(1 for s in seqs if any(s in t or hh.revcomp(s) in t for t in text))


def test_new_background_comes_from_the_genome(runs, inputs, windows, genome, library):
	"""-useNewBg runs homer2 bg with -N and -NN instead of reading preparsed
	files: background windows are genome sequence that overlaps neither a
	target nor another background window."""

	assert not os.path.exists(inputs['newbg'])
	assert '-NN 20000 -N 1000' in bg_command(runs, 'newbg')
	targets = fh.out_fasta(runs, 'newbg', 'homerBg.target.sequences.fasta')
	assert sorted(targets.values()) == sorted(windows.values())
	background = fh.out_fasta(runs, 'newbg', 'homerBg.bg.sequences.fasta')
	assert in_genome(background.values(), genome) == len(background) > 0
	target_bed = bed(fh.out(runs, 'newbg', 'homerBg.target.positions.bed'))
	bg_bed = sorted(bed(fh.out(runs, 'newbg', 'homerBg.bg.positions.bed')))
	assert not any(overlaps(b, t) for b in bg_bed for t in target_bed)
	assert not any(overlaps(a, b) for a, b in zip(bg_bed, bg_bed[1:]))
	fh.check_known_run(runs, 'newbg', library)


def test_model_background_is_generated(runs, genome):
	"""-modelBg generates sequences instead of taking them from the genome;
	a later -genomeBg switches back."""

	assert ' -model ' in bg_command(runs, 'newbg_model')
	background = fh.out_fasta(runs, 'newbg_model', 'homerBg.bg.sequences.fasta').values()
	assert len(background) > 0 and in_genome(background, genome) == 0
	assert '-model' not in bg_command(runs, 'newbg_model_genome')
	background = fh.out_fasta(runs, 'newbg_model_genome', 'homerBg.bg.sequences.fasta').values()
	assert in_genome(background, genome) == len(background) > 0


def test_num_bins_sets_the_gc_bins(runs):
	assert '-nbins 3 ' in bg_command(runs, 'newbg_bins')
	for key, n in (('newbg', 10), ('newbg_bins', 3)):
		header, rows = hh.read_tsv(fh.out(runs, key, 'homerBg.target.stats.txt'), header=True)
		assert {int(row[header.index('bin')]) for row in rows} == set(range(n)), key


@pytest.mark.parametrize('key,k,positional', [('newbg', 2, False), ('newbg_ikmer', 3, False),
	('newbg_pkmer', 2, True)])
def test_kmer_matching_options(runs, key, k, positional):
	"""-ikmer <k> matches k-mer content over the whole sequence (one row of
	k-mer frequencies), -pkmer <k> per position (one row per position)."""

	header, rows = hh.read_tsv(fh.out(runs, key, 'homerBg.target.freq.txt'), header=True)
	assert len(header) == 1 + 4 ** k and {len(mer) for mer in header[1:]} == {k}
	labels = [row[0] for row in rows]
	if positional:
		assert labels == [str(i) for i in range(len(labels))] and len(labels) > 150
	else:
		assert labels == ['PositionIndependent']
	assert ('-pkmer' if positional else '-ikmer') + ' {}'.format(k) in bg_command(runs, key)


@pytest.mark.parametrize('key,within', [('newbg_target_overlaps', False),
	('newbg_bg_overlaps', True)])
def test_overlap_options(runs, key, within):
	"""-allowTargetOverlaps lets background windows overlap targets,
	-allowBgOverlaps lets them overlap one another; with 1,000 windows in
	this small genome some do."""

	assert ('-allowBgOverlaps' if within else '-allowTargetOverlaps') in bg_command(runs, key)
	bg_bed = sorted(bed(fh.out(runs, key, 'homerBg.bg.positions.bed')))
	if within:
		assert any(overlaps(a, b) for a, b in zip(bg_bed, bg_bed[1:]))
	else:
		target_bed = bed(fh.out(runs, key, 'homerBg.target.positions.bed'))
		assert any(overlaps(b, t) for b in bg_bed for t in target_bed)


def test_new_background_with_a_genome_directory(runs, windows, genome):
	"""With -useNewBg a genome directory is joined into one FASTA file in the
	output directory for homer2 bg."""

	command = bg_command(runs, 'newbg_genome_dir')
	joined = command.split('-g "')[1].split('"')[0]
	path = os.path.join(runs['newbg_genome_dir'].cwd, joined)
	assert {name: seq for name, seq in hh.read_fasta(path)} == genome
	targets = fh.out_fasta(runs, 'newbg_genome_dir', 'homerBg.target.sequences.fasta')
	assert sorted(targets.values()) == sorted(windows.values())

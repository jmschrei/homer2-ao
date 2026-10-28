"""Reference computations for the findMotifsGenome.pl tests.

Each function restates a convention of a script or program findMotifsGenome.pl
runs (resizePosFile.pl, removePoorSeq.pl, homerTools freq, freq2group.pl,
assignGeneWeights.pl, homer2 known and norm, findKnownMotifs.pl,
Statistics.pm), so that expected values are computed from the inputs and not
read back from HOMER's own output.
"""

import math
import os

import pytest

import homer_helpers as hh

KNOWN_SMALL = os.path.join(hh.MOTIFS, 'known_small.motifs')

# freq2group.pl's bin edges for GC (column 2 of homerTools freq -gc) and CpG
# (column 1).
GC_EDGES = (0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8)
CPG_EDGES = (0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16)


def command(options, peaks=hh.PEAKS, genome=hh.GENOME, preparsed=None):
	"""findMotifsGenome.pl's arguments; results go to `out` under the call's
	working directory."""

	args = ['findMotifsGenome.pl', peaks, genome, 'out']
	if preparsed is not None:
		args += ['-preparsedDir', preparsed]
	return [str(arg) for arg in args + list(options)]


def out(results, key, *parts):
	"""Path of a file in the output directory of run `key`, after checking
	that the run exited 0."""

	result = results[key]
	assert result.returncode == 0, '{} exited {}\n{}'.format(key, result.returncode,
		result.stderr[-2000:])
	return os.path.join(result.cwd, 'out', *parts)


def run_calls(tmp, calls, jobs=None):
	"""Run {key: (options, keyword arguments of command())} at once, each in
	its own directory under `tmp`; returns key -> CompletedProcess with
	`argv` and `cwd` attributes. A failed run fails only the tests that read
	it."""

	prepared = {}
	for key, (options, kwargs) in calls.items():
		prepared[key] = (command(options, **kwargs), str(tmp / key))
	results = hh.run_many(prepared, check=False, jobs=jobs)
	for key, result in results.items():
		result.argv, result.cwd = prepared[key]
	return results


def out_fasta(results, key, name):
	"""A FASTA file of run `key`'s output directory, upper-cased, by name."""

	return read_fasta_upper(out(results, key, name))


def known_rows(results, key):
	"""knownResults.txt of run `key` by motif name."""

	return {row['name']: row for row in hh.read_known_results(out(results, key,
		'knownResults.txt'))}


def check_known_run(results, key, library, **kwargs):
	"""knownResults.txt of a -keepFiles run against its seq.tsv and group.adj
	(check_known_results)."""

	return check_known_results(out(results, key, 'knownResults.txt'),
		read_seq_tsv(out(results, key, 'seq.tsv')),
		read_groups(out(results, key, 'group.adj')), library, **kwargs)


def same_file(results, key_a, key_b, name):
	return open(out(results, key_a, name)).read() == open(out(results, key_b, name)).read()


# --- regions and sequences --------------------------------------------------------

def window(peak, size, move=0):
	"""(start, end), 1-based and inclusive, of the region extracted for a peak
	(resizePosFile.pl): the centre (start + end) / 2 moves by `move` along the
	peak's strand, and the window runs from floor(centre - size/2) to
	floor(centre + size/2), size + 1 bases."""

	centre = (peak['start'] + peak['end']) / 2.0
	centre += -move if peak['strand'] == '-' else move
	return int(math.floor(centre - size / 2.0)), int(math.floor(centre + size / 2.0))


def region_sequence(chrom, start, end, strand='+', genome=None):
	genome = genome or hh.read_genome()
	seq = genome[chrom][start - 1:end]
	return hh.revcomp(seq) if strand == '-' else seq


def window_sequence(peak, size, move=0, genome=None):
	start, end = window(peak, size, move)
	return region_sequence(peak['chr'], start, end, peak['strand'], genome)


def masked(seq):
	"""Soft-masked (lowercase) bases as N, the rest upper case."""

	return ''.join('N' if base.islower() else base for base in seq).upper()


def kept_by_max_n(seq, max_n=0.7):
	"""removePoorSeq.pl drops a sequence when more than `max_n` of it is N."""

	return 1.0 - (len(seq) - seq.upper().count('N')) / float(len(seq)) <= max_n


def read_fasta_upper(path):
	return {name: seq.upper() for name, seq in hh.read_fasta(path)}


def read_seq_tsv(path):
	return {row[0]: row[1] for row in hh.read_tsv(path)}


def read_groups(path):
	"""Group file as {name: (group, weight)}; a missing weight is 1, as
	homer2 reads it."""

	groups = {}
	for row in hh.read_tsv(path):
		groups[row[0]] = (int(row[1]), float(row[2]) if len(row) > 2 else 1.0)
	return groups


# --- sequence content ---------------------------------------------------------------

def gc_cpg(seq):
	"""(CpG, GC) as homerTools freq -gc reports them: CG dinucleotides per
	dinucleotide of A/C/G/T, and G + C per A/C/G/T base."""

	seq = seq.upper()
	bases = [base for base in seq if base in 'ACGT']
	pairs = [seq[i:i + 2] for i in range(len(seq) - 1)
		if seq[i] in 'ACGT' and seq[i + 1] in 'ACGT']
	cpg = pairs.count('CG') / max(float(len(pairs)), 0.01)
	gc = sum(1 for base in bases if base in 'GC') / float(len(bases))
	return cpg, gc


def content_bin(value, edges):
	"""freq2group.pl: the index of the first edge above the value (as printed,
	six decimals), or len(edges) when none is."""

	value = float('{:f}'.format(value))
	for i, edge in enumerate(edges):
		if value < edge:
			return i
	return len(edges)


def gc_bin(seq):
	return content_bin(gc_cpg(seq)[1], GC_EDGES)


def cpg_bin(seq):
	return content_bin(gc_cpg(seq)[0], CPG_EDGES)


def gc_matched(targets, background, bin_of=gc_bin):
	"""The targets and background sequences assignGeneWeights.pl keeps: those
	in content bins holding at least one target and one background
	sequence."""

	target_bins = {name: bin_of(seq) for name, seq in targets.items()}
	background_bins = {name: bin_of(seq) for name, seq in background.items()}
	shared = set(target_bins.values()) & set(background_bins.values())
	return ({name: seq for name, seq in targets.items() if target_bins[name] in shared},
		{name: seq for name, seq in background.items() if background_bins[name] in shared})


def oligo_levels(seqs, k, both_strands=True, last=True):
	"""Weighted oligo frequencies as homer2 norm defines them: each sequence
	adds weight / (number of its k-mers without N) for each such k-mer; with
	both strands a k-mer's level is the mean of its own and its reverse
	complement's. With `last` False the sequence's last k-mer is left out, as
	homer2 norm's loop (positions 0 to length - k - 1) does."""

	levels = {}
	for seq, weight in seqs:
		end = len(seq) - k + (1 if last else 0)
		mers = [seq[j:j + k] for j in range(end) if 'N' not in seq[j:j + k]]
		for mer in mers:
			levels[mer] = levels.get(mer, 0.0) + weight / len(mers)
	total = sum(levels.values())
	result = {}
	for mer in levels:
		level = levels[mer]
		if both_strands:
			level = (level + levels.get(hh.revcomp(mer), 0.0)) / 2.0
		result[mer] = level / total
	return result


def autonorm_target_texts(seqs, lengths, both_strands=True):
	"""{oligo: accepted TargetCounts texts} for seq.autonorm.tsv: the
	percentage with two decimals, with or without each sequence's last
	k-mer counted."""

	texts = {}
	for k in lengths:
		for last in (True, False):
			for mer, level in oligo_levels(seqs, k, both_strands, last).items():
				texts.setdefault(mer, set()).add('{:.2f}%'.format(100 * level))
	return texts


def oligo_distance(seqs_a, seqs_b, k):
	"""Sum over k-mers of the absolute difference in weighted frequency, both
	strands counted."""

	a = oligo_levels(seqs_a, k)
	b = oligo_levels(seqs_b, k)
	return sum(abs(a.get(mer, 0.0) - b.get(mer, 0.0)) for mer in set(a) | set(b))


# --- motif enrichment -----------------------------------------------------------------

_WORDS = {}
_STRAND_WORDS = {}


def site_words(motif, limit=200000):
	"""Every + strand word whose score, summed position by position as
	homer_helpers.scan sums it, reaches the motif's threshold; None when there
	are more than `limit`. A - strand site is a word whose reverse complement
	is in the set."""

	key = (tuple(tuple(row) for row in motif.rows), motif.threshold)
	if key not in _WORDS:
		lo = motif.log_odds()
		width = len(lo)
		rest = [0.0] * (width + 1)
		for j in range(width - 1, -1, -1):
			rest[j] = rest[j + 1] + max(lo[j])
		words = []
		stack = [('', 0.0)]
		while stack and len(words) <= limit:
			prefix, score = stack.pop()
			j = len(prefix)
			if j == width:
				if score >= motif.threshold:
					words.append(prefix)
				continue
			for code, base in enumerate('ACGT'):
				value = score + lo[j][code]
				if value + rest[j + 1] >= motif.threshold - 1e-9:
					stack.append((prefix + base, value))
		_WORDS[key] = set(words) if len(words) <= limit else None
	return _WORDS[key]


def has_site(seq, motif, strand='both'):
	"""Whether homer_helpers.scan would report a site of `motif` in `seq`."""

	plus = site_words(motif)
	if plus is None:
		return bool(hh.scan(seq, motif, strand))
	key = (id(plus), strand)
	if key not in _STRAND_WORDS:
		words = set()
		if strand in ('both', '+'):
			words |= plus
		if strand in ('both', '-'):
			words |= {hh.revcomp(word) for word in plus}
		_STRAND_WORDS[key] = words
	words = _STRAND_WORDS[key]
	seq = hh.clean_sequence(seq)
	width = len(motif)
	return any(seq[i:i + width] in words for i in range(len(seq) - width + 1))


def weighted_hits(seqs, groups, motif, strand='both'):
	"""(target, background): summed group-file weights of the sequences with
	at least one site, as homer2 known counts them."""

	target = background = 0.0
	for name, (group, weight) in groups.items():
		if name in seqs and has_site(seqs[name], motif, strand):
			if group == 1:
				target += weight
			else:
				background += weight
	return target, background


def group_totals(seqs, groups):
	target = sum(w for name, (g, w) in groups.items() if g == 1 and name in seqs)
	background = sum(w for name, (g, w) in groups.items() if g == 0 and name in seqs)
	return target, background


def log_pvalue(hits, background_hits, total, background_total, hypergeometric=False):
	"""homer2's enrichment ln p-value. Counts and totals are truncated to
	integers. Binomial: ln P(X >= hits), X ~ Binomial(total, r), r =
	background_hits / background_total held within [1/N, (N-1)/N].
	Hypergeometric: ln P(X >= hits) for hits + background_hits draws from
	total + background_total with total successes."""

	t, b = int(hits), int(background_hits)
	n, m = int(total), int(background_total)
	if hypergeometric:
		return hh.log_hypergeometric_tail(t, t + b, n, n + m)
	ratio = min(max(b / float(m), 1.0 / m), (m - 1.0) / m)
	return hh.log_binomial_tail(t, n, ratio)


def benjamini(pvalues):
	"""Statistics.pm benjaminiFDR: p * n / rank over sorted p-values, ties
	sharing the rank of their last member, made non-decreasing and capped
	at 1."""

	n = len(pvalues)
	order = sorted(range(n), key=lambda i: pvalues[i])
	qvalues = [1.0] * n
	last, last_index = -1e10, -1
	for k in range(n):
		if k < n - 1 and pvalues[order[k + 1]] == pvalues[order[k]]:
			continue
		q = min(max(pvalues[order[k]] * n / (k + 1), last), 1.0)
		for j in range(k, last_index, -1):
			qvalues[order[j]] = q
		last, last_index = q, k
	return qvalues


def pvalue_texts(logp):
	"""The P-value column: 1e<ln p / 2.3026, truncated>; both neighbours of a
	value this close to a boundary are accepted."""

	return {'1e{}'.format(int(logp * f / 2.3026)) for f in (1 - 1e-6, 1 + 1e-6)}


def check_known_results(path, seqs, groups, motifs, strand='both', hypergeometric=False):
	"""Every column of knownResults.txt against the reference: counts from a
	scan of `seqs` weighted by `groups`, percentages of the weighted totals,
	the p-value from the counts, Benjamini q-values over all rows, rows ranked
	by p-value."""

	rows = hh.read_known_results(path)
	assert sorted(row['name'] for row in rows) == sorted(m.name for m in motifs)
	total, background_total = group_totals(seqs, groups)
	by_name = {row['name']: row for row in rows}
	for motif in motifs:
		row = by_name[motif.name]
		hits, background_hits = weighted_hits(seqs, groups, motif, strand)
		assert row['targets'] == pytest.approx(hits, abs=0.051), motif.name
		assert row['background'] == pytest.approx(background_hits, abs=0.051), motif.name
		assert row['target_percent'] == pytest.approx(100.0 * hits / total, abs=0.0051)
		assert row['background_percent'] == pytest.approx(
			100.0 * background_hits / background_total, abs=0.0051)
		logp = log_pvalue(hits, background_hits, total, background_total, hypergeometric)
		assert row['logp'] == pytest.approx(logp, rel=1e-3, abs=1e-3), motif.name
		assert row['pvalue'] in pvalue_texts(logp), motif.name
	logps = [row['logp'] for row in rows]
	assert logps == sorted(logps)
	qvalues = benjamini([math.exp(lp) for lp in logps])
	for row, q in zip(rows, qvalues):
		assert row['qvalue'] == pytest.approx(q, rel=5e-3, abs=2e-4), row['name']
	return rows


def empirical_fdr(real, random, randomizations):
	"""Statistics.pm empiricalFDR2: for the i-th best real ln p-value (rank r
	from 1), (number of random values below it + 1) / randomizations / r,
	made non-decreasing and capped at 1; returned in the order of `real`."""

	order = sorted(range(len(real)), key=lambda i: real[i])
	random = sorted(random)
	fdr = [1.0] * len(real)
	index, last = 0, -1e10
	for rank, i in enumerate(order):
		while index < len(random) and random[index] < real[i]:
			index += 1
		value = min(max((index + 1.0) / randomizations / (rank + 1), last), 1.0)
		fdr[i] = last = value
	return fdr


# --- planted sites ------------------------------------------------------------------------

def peaks_with_overlapping_site(motif, motif_name='siteA', size=200, genome=None):
	"""IDs of peaks carrying a planted `motif_name` site that a site of
	`motif` in the peak's window overlaps, mapped back to the genome."""

	genome = genome or hh.read_genome()
	peaks = {p['id']: p for p in hh.read_peaks()}
	found = set()
	for plant in hh.read_planted():
		if plant['motif'] != motif_name:
			continue
		peak = peaks[plant['peak']]
		start, end = window(peak, size)
		seq = region_sequence(peak['chr'], start, end, peak['strand'], genome)
		width = len(motif)
		for pos, _, site_strand, _ in hh.scan(seq, motif, offset=0):
			first = pos if site_strand == '+' else pos - width + 1
			last = first + width - 1
			if peak['strand'] == '+':
				lo, hi = start - 1 + first, start - 1 + last
			else:
				lo, hi = end - 1 - last, end - 1 - first
			if lo < int(plant['end']) and int(plant['start']) <= hi:
				found.add(peak['id'])
				break
	return found

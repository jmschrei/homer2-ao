"""scanMotifGenomeWide.pl <motif> <genome FASTA>: every site of a motif in a
genome, in genome coordinates.

The script uppercases each chromosome (with -mask it turns lowercase bases
into N instead), cuts it into overlapping pieces of about 10 kb, runs
homer2 find -offset 0 on each piece and converts the sites to 1-based,
inclusive genome coordinates, dropping the copies a site has in two pieces.
Expected values come from homer_helpers.scan() over each whole chromosome,
which gives the same sites, converted the same way: a + strand site starts at
its first base, a - strand site (reported by homer2 at its last base) starts
length - 1 bases before. The script then sorts a chromosome's sites by start
and, without -keepAll, drops a site when a site of the same motif whose middle
(floor((start + end) / 2)) is at most half the motif's length away scores
higher, or scores at least as high and comes later; reduce() below follows
its loop. Each call runs once in a module-scoped batch.
"""

import math
import os

import pytest

import homer_helpers as hh


def motif_path(name):
	ext = '.motif' if os.path.exists(os.path.join(hh.MOTIFS, name + '.motif')) else '.motifs'
	return os.path.join(hh.MOTIFS, name + ext)


def genome_sites(genome, motifs, mask=False):
	"""{chromosome: [site]} with each site a dict: name, chr, start, end,
	strand, score (as printed, six decimals), seq (the genome's + strand
	text, uppercase) and mid."""

	out = {}
	for chrom, seq in genome.items():
		if mask:
			seq = ''.join('N' if base.islower() else base for base in seq)
		sites = []
		for m in motifs:
			for pos, site, strand, score in hh.scan(seq, m, offset=0):
				start = pos + 1 if strand == '+' else pos - len(m) + 2
				end = start + len(m) - 1
				sites.append({'name': m.name, 'chr': chrom, 'start': start, 'end': end,
					'strand': strand, 'score': float('{:.6f}'.format(score)), 'seq': site,
					'mid': (start + end) // 2})
		out[chrom] = sites
	return out


def five_prime(site):
	"""-5p: 201 bp centred on the site's 5' end (its start on the + strand,
	its end on the -)."""

	centre = site['start'] if site['strand'] == '+' else site['end']
	return dict(site, start=centre - 100, end=centre + 100)


def reduce(sites):
	"""The sites kept without -keepAll, from the start-sorted list: walking
	back from site i while the middles are at most length / 2 apart, a
	higher-scoring site of the same motif drops it; walking forward, one
	scoring at least as high does."""

	kept = []
	for i, site in enumerate(sites):
		half = len(site['seq']) / 2.0
		bad = False
		for j in range(i - 1, -1, -1):
			if site['mid'] - sites[j]['mid'] > half:
				break
			if sites[j]['name'] == site['name'] and site['score'] < sites[j]['score']:
				bad = True
				break
		if not bad:
			for j in range(i + 1, len(sites)):
				if sites[j]['mid'] - site['mid'] > half:
					break
				if sites[j]['name'] == site['name'] and site['score'] <= sites[j]['score']:
					bad = True
					break
		if not bad:
			kept.append(site)
	return kept


def report(genome, motifs, mask=False, keep_all=False, p5=False):
	"""The sites the script prints, in its order, as dicts with an `id`
	(1, 2, ... per motif in output order)."""

	out = []
	ids = {}
	for chrom, sites in genome_sites(genome, motifs, mask).items():
		if p5:
			sites = [five_prime(s) for s in sites]
		starts = [s['start'] for s in sites]
		assert len(set(starts)) == len(starts), 'no two sites should share a start'
		sites = sorted(sites, key=lambda s: s['start'])
		for site in (sites if keep_all else reduce(sites)):
			if site['start'] < 0 or site['start'] >= site['end']:
				continue
			ids[site['name']] = ids.get(site['name'], 0) + 1
			out.append(dict(site, id=ids[site['name']]))
	return out


def parse(text):
	rows = []
	for line in text.splitlines():
		cols = line.split('\t')
		assert len(cols) == 7, line
		name, ident = cols[0].rsplit('-', 1)
		rows.append({'name': name, 'id': int(ident), 'chr': cols[1], 'start': int(cols[2]),
			'end': int(cols[3]), 'strand': cols[4], 'score': float(cols[5]), 'seq': cols[6]})
	return rows


KEYS = ('name', 'id', 'chr', 'start', 'end', 'strand', 'seq')


def assert_report(got, want):
	assert [[g[k] for k in KEYS] for g in got] == [[w[k] for k in KEYS] for w in want]
	for g, w in zip(got, want):
		assert g['score'] == pytest.approx(w['score'], abs=1e-5), g


@pytest.fixture(scope='module')
def runs(tmp_path_factory, genome):
	"""key -> CompletedProcess. The genome is also written as a directory of
	one FASTA file per chromosome (50 bases per line, next to a file the
	script must ignore), and siteA is copied under the name "site A/test"."""

	root = tmp_path_factory.mktemp('scanMotifGenomeWide')
	directory = root / 'genome'
	directory.mkdir()
	for chrom, seq in genome.items():
		hh.write_fasta(str(directory / (chrom + '.fa')), [(chrom, seq)], width=50)
	(directory / 'notes.txt').write_text('>chrX\nACGT\n')
	site_a = hh.read_motifs(motif_path('siteA'))[0]
	renamed = hh.write_motifs(str(root / 'renamed.motif'),
		[hh.Motif(site_a.consensus, 'site A/test', site_a.threshold, site_a.rows)])
	a = motif_path('siteA')
	calls = {
		'default': [a, hh.GENOME],
		'two': [motif_path('siteA_siteB'), hh.GENOME],
		'keepall': [a, hh.GENOME, '-keepAll'],
		'keepall_known': [motif_path('known_small'), hh.GENOME, '-keepAll'],
		'bed': [renamed, hh.GENOME, '-bed'],
		'bed_int': [renamed, hh.GENOME, '-bed', '-int'],
		'p5': [a, hh.GENOME, '-5p'],
		'mask': [a, hh.GENOME, '-mask'],
		'threads': [a, hh.GENOME, '-p', 2],
		'homer2': [a, hh.GENOME, '-homer2'],
		'homer1': [a, hh.GENOME, '-homer1'],
		'directory': [a, str(directory)],
	}
	todo = {key: (['scanMotifGenomeWide.pl'] + args, root / key) for key, args in calls.items()}
	results = hh.run_many(todo, check=False)
	for key, result in results.items():
		result.ok = result.returncode == 0
	return results


def output(runs, key):
	result = runs[key]
	assert result.ok, result.stderr[-2000:]
	return result.stdout


@pytest.fixture(scope='module')
def site_a():
	return hh.read_motifs(motif_path('siteA'))


def test_default_output(runs, genome, site_a):
	"""One line per site: motif-<n>, chr, start, end (1-based, inclusive),
	strand, score and the genome's + strand text, chromosomes in file order
	and sites by start."""

	got = parse(output(runs, 'default'))
	assert_report(got, report(genome, site_a))


def test_reported_intervals_hold_the_site(runs, genome):
	"""Each interval's genome sequence is the reported text; read on the
	reported strand it scores the reported score."""

	motifs = {m.name: m for m in hh.read_motifs(motif_path('siteA_siteB'))}
	rows = parse(output(runs, 'two'))
	assert {r['name'] for r in rows} == set(motifs)
	for r in rows:
		text = genome[r['chr']][r['start'] - 1:r['end']].upper()
		assert text == r['seq'], r
		if r['strand'] == '-':
			text = hh.revcomp(text)
		hits = hh.scan(text, motifs[r['name']], '+')
		assert len(hits) == 1 and hits[0][3] == pytest.approx(r['score'], abs=1e-5), r


def test_planted_sites_are_all_reported(runs):
	"""Each of the 120 planted siteA copies, lowercase ones included, on the
	strand it was planted on."""

	rows = parse(output(runs, 'default'))
	found = {(r['chr'], r['start'], r['end'], r['strand']) for r in rows}
	planted = [p for p in hh.read_planted() if p['motif'] == 'siteA']
	assert len(planted) == 120
	for p in planted:
		assert (p['chr'], int(p['start']) + 1, int(p['end']), p['strand']) in found, p


def test_several_motifs_are_reduced_per_motif(runs, genome):
	"""Two motifs in one file: each keeps its own numbering, and a site is
	only ever dropped for a better site of the same motif."""

	got = parse(output(runs, 'two'))
	assert_report(got, report(genome, hh.read_motifs(motif_path('siteA_siteB'))))


def test_keep_all_reports_every_site(runs, genome, site_a):
	got = parse(output(runs, 'keepall'))
	assert_report(got, report(genome, site_a, keep_all=True))


def test_overlapping_sites_are_reduced_to_the_best(runs, site_a):
	"""Without -keepAll a site is left out only when a site of the same
	motif within half its length (comparing middles) scores at least as
	high, and the planted copies (exact matches) all survive."""

	every = parse(output(runs, 'keepall'))
	kept = parse(output(runs, 'default'))
	key = lambda r: (r['chr'], r['start'], r['strand'])
	kept_keys = {key(r) for r in kept}
	assert kept_keys < {key(r) for r in every}
	half = len(site_a[0]) / 2.0
	for r in every:
		if key(r) in kept_keys:
			continue
		mid = (r['start'] + r['end']) // 2
		assert any(k['chr'] == r['chr'] and abs((k['start'] + k['end']) // 2 - mid) <= half
			and k['score'] >= r['score'] for k in kept), r


def test_keep_all_with_a_palindrome_and_several_motifs(runs, genome):
	"""-keepAll with five motifs, the palindromic E-box among them (found on
	both strands at the same place): every site, numbered 1, 2, ... per
	motif."""

	motifs = hh.read_motifs(motif_path('known_small'))
	got = parse(output(runs, 'keepall_known'))
	want = []
	for sites in genome_sites(genome, motifs).values():
		want += sites
	fields = ('name', 'chr', 'start', 'end', 'strand', 'seq')
	assert sorted(tuple(g[k] for k in fields) for g in got) == sorted(
		tuple(w[k] for k in fields) for w in want)
	assert {m.name for m in motifs} >= {g['name'] for g in got} >= {'siteA', 'siteB', 'ebox'}
	for name in {g['name'] for g in got}:
		assert [g['id'] for g in got if g['name'] == name] == list(range(1,
			sum(1 for g in got if g['name'] == name) + 1))
	ebox = [g for g in got if g['name'] == 'ebox']
	assert {(g['chr'], g['start']) for g in ebox if g['strand'] == '+'} == {
		(g['chr'], g['start']) for g in ebox if g['strand'] == '-'}


def bed_rows(text):
	return [line.split('\t') for line in text.splitlines()]


def test_bed_output(runs, genome, site_a):
	"""-bed: chr, start (0-based), end, the motif name cut at the first "/"
	with whitespace made "_", the score, strand; no other lines."""

	want = report(genome, site_a)
	got = bed_rows(output(runs, 'bed'))
	assert [r[:4] + r[5:] for r in got] == [[w['chr'], str(w['start'] - 1), str(w['end']),
		'site_A', w['strand']] for w in want]
	for r, w in zip(got, want):
		assert float(r[4]) == pytest.approx(w['score'], abs=1e-5)


def test_bed_int_rounds_scores(runs, genome, site_a):
	"""-int: scores rounded to the nearest integer, kept within 0-1000."""

	want = report(genome, site_a)
	got = bed_rows(output(runs, 'bed_int'))
	assert [r[4] for r in got] == [str(min(1000, max(0, int(math.floor(w['score'] + 0.5)))))
		for w in want]
	assert [r[:4] for r in got] == [r[:4] for r in bed_rows(output(runs, 'bed'))]


def test_5p_centres_on_the_motif_start(runs, genome, site_a):
	"""-5p: 201 bp from 100 before to 100 after the site's 5' end; sites are
	sorted by that start before overlapping ones are reduced."""

	got = parse(output(runs, 'p5'))
	assert_report(got, report(genome, site_a, p5=True))
	plain = parse(output(runs, 'default'))
	assert len(got) == len(plain)
	ends = {(r['chr'], r['start'] if r['strand'] == '+' else r['end']) for r in plain}
	assert {(r['chr'], r['start'] + 100) for r in got} == ends


def test_mask_skips_lowercase_sequence(runs, genome, site_a):
	"""-mask: lowercase bases are N, so the ten planted copies in lowercase
	peaks disappear and the other 110 remain."""

	got = parse(output(runs, 'mask'))
	assert_report(got, report(genome, site_a, mask=True))
	found = {(r['chr'], r['start']) for r in got}
	planted = [p for p in hh.read_planted() if p['motif'] == 'siteA']
	for p in planted:
		assert ((p['chr'], int(p['start']) + 1) in found) == (p['lowercase'] == 'no'), p


def test_threads_give_the_same_sites(runs):
	assert sorted(output(runs, 'threads').splitlines()) == sorted(
		output(runs, 'default').splitlines())


def test_homer2_is_the_default(runs):
	assert output(runs, 'homer2') == output(runs, 'default')


def test_homer1_finds_the_same_sites(runs, genome, site_a):
	"""-homer1 scans with the original homer program."""

	assert_report(parse(output(runs, 'homer1')), report(genome, site_a))


def test_genome_directory(runs):
	"""A directory of per-chromosome FASTA files (the files matching *fa*)
	gives the same output as one file."""

	assert output(runs, 'directory') == output(runs, 'default')

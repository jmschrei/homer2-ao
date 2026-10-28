"""homer2 find: every site of each motif in a set of sequences.

Expected sites come from homer_helpers.scan(), a Python scan with the same
scoring: ln(p / 0.25) summed over the site, a site when the sum reaches the
motif's threshold, positions relative to -offset (default: the middle of the
sequence), and windows holding anything but A, C, G or T skipped.
"""

import os

import pytest

import homer_helpers as hh


def find_lines(text):
	"""homer2 find output as a sorted list of (name, position, site, motif,
	strand, score) tuples; the order of lines is not part of the output's
	meaning and changes with -p."""

	rows = []
	for line in text.splitlines():
		if not line.strip():
			continue
		cols = line.split('\t')
		assert len(cols) == 6, line
		rows.append((cols[0], int(cols[1]), cols[2], cols[3], cols[4], float(cols[5])))
	return sorted(rows)


def expected_lines(records, motifs, strand='both', offset=None):
	rows = []
	for motif in motifs:
		for name, seq in records:
			for pos, site, site_strand, score in hh.scan(seq, motif, strand, offset):
				rows.append((name, pos, site, motif.name, site_strand, score))
	return sorted(rows)


def assert_same_sites(got, want):
	assert [row[:5] for row in got] == [row[:5] for row in want]
	for g, w in zip(got, want):
		assert g[5] == pytest.approx(w[5], abs=1e-6), (g, w)


@pytest.fixture(scope='module')
def targets():
	return hh.read_fasta(hh.TARGETS)


@pytest.fixture(scope='module')
def targets_tsv(tmp_path_factory, targets):
	path = tmp_path_factory.mktemp('find') / 'targets.tsv'
	return hh.write_seq_tsv(str(path), targets)


@pytest.mark.parametrize('motif', ['siteA', 'siteB', 'ebox'])
def test_tsv_input_reports_every_site(tmp_path, targets, targets_tsv, motif):
	path = os.path.join(hh.MOTIFS, motif + '.motif')
	result = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path], tmp_path)
	want = expected_lines(targets, hh.read_motifs(path))
	assert want, 'the inputs should hold sites of ' + motif
	assert_same_sites(find_lines(result.stdout), want)


def test_planted_sites_are_all_found(tmp_path, targets_tsv):
	"""Each of the 120 planted siteA copies is reported, lowercase ones
	included, on the strand it was planted on."""

	result = hh.run(['homer2', 'find', '-s', targets_tsv, '-m',
		os.path.join(hh.MOTIFS, 'siteA.motif')], tmp_path)
	found = {(row[0], row[4]) for row in find_lines(result.stdout)}
	planted = [p for p in hh.read_planted() if p['motif'] == 'siteA']
	assert len(planted) == 120
	for p in planted:
		assert (p['peak'], p['strand']) in found
	lowercase = [p for p in planted if p['lowercase'] == 'yes']
	assert len(lowercase) == 10


def test_fasta_input_matches_tsv_input(tmp_path, targets, targets_tsv):
	path = os.path.join(hh.MOTIFS, 'siteA_siteB.motifs')
	tsv = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path], tmp_path)
	fasta = hh.run(['homer2', 'find', '-i', hh.TARGETS, '-m', path], tmp_path)
	assert find_lines(fasta.stdout) == find_lines(tsv.stdout)
	assert_same_sites(find_lines(fasta.stdout), expected_lines(targets, hh.read_motifs(path)))


def test_several_motifs_in_one_file(tmp_path, targets, targets_tsv):
	path = os.path.join(hh.MOTIFS, 'known_small.motifs')
	result = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path], tmp_path)
	motifs = hh.read_motifs(path)
	assert len(motifs) == 5
	got = find_lines(result.stdout)
	assert_same_sites(got, expected_lines(targets, motifs))
	assert {row[3] for row in got} >= {'siteA', 'siteB', 'ebox'}


@pytest.mark.parametrize('offset', [0, 50, -25])
def test_offset(tmp_path, targets, targets_tsv, offset):
	"""-offset <#> sets where positions are counted from: a site's position
	is its 0-based index plus the offset."""

	path = os.path.join(hh.MOTIFS, 'siteA.motif')
	result = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path,
		'-offset', offset], tmp_path)
	assert_same_sites(find_lines(result.stdout),
		expected_lines(targets, hh.read_motifs(path), offset=offset))


def test_default_offset_is_the_middle(tmp_path, targets_tsv):
	"""Without -offset, positions are relative to the middle: the 200 bp
	targets' positions are the -offset 0 positions minus 100."""

	path = os.path.join(hh.MOTIFS, 'siteB.motif')
	default = find_lines(hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path],
		tmp_path).stdout)
	zero = find_lines(hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path,
		'-offset', 0], tmp_path).stdout)
	assert default == sorted((r[0], r[1] - 100) + r[2:] for r in zero)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_strand(tmp_path, targets, targets_tsv, strand):
	path = os.path.join(hh.MOTIFS, 'siteA.motif')
	result = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path,
		'-strand', strand], tmp_path)
	got = find_lines(result.stdout)
	assert {row[4] for row in got} == {strand}
	assert_same_sites(got, expected_lines(targets, hh.read_motifs(path), strand=strand))


def test_both_strands_of_a_palindrome(tmp_path, targets_tsv):
	"""CACGTG reads the same on both strands, so each exact match is
	reported twice: + at its first base and - at its last."""

	path = os.path.join(hh.MOTIFS, 'ebox.motif')
	got = find_lines(hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path,
		'-offset', 0], tmp_path).stdout)
	plus = {(r[0], r[1]) for r in got if r[4] == '+'}
	minus = {(r[0], r[1] - 5) for r in got if r[4] == '-'}
	assert plus and plus == minus


def test_output_file(tmp_path, targets_tsv):
	path = os.path.join(hh.MOTIFS, 'siteA.motif')
	out = tmp_path / 'sites.txt'
	result = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path, '-o', out], tmp_path)
	assert result.stdout == ''
	stdout = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path], tmp_path).stdout
	assert find_lines(out.read_text()) == find_lines(stdout)


@pytest.mark.parametrize('threads', [2, 4])
def test_threads_give_the_same_sites(tmp_path, targets_tsv, threads):
	path = os.path.join(hh.MOTIFS, 'known_small.motifs')
	one = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path, '-p', 1], tmp_path)
	many = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path, '-p', threads], tmp_path)
	assert find_lines(many.stdout) == find_lines(one.stdout)


def mscore_lines(text):
	rows = {}
	for line in text.splitlines():
		cols = line.split('\t')
		assert len(cols) == 6, line
		assert cols[0] not in rows
		rows[cols[0]] = (int(cols[1]), cols[2], cols[3], cols[4], float(cols[5]))
	return rows


@pytest.mark.parametrize('motif', ['siteA', 'ebox'])
def test_mscore_reports_the_best_window(tmp_path, targets, targets_tsv, motif):
	"""-mscore prints one line per sequence: the best-scoring window's
	position and strand, the motif's consensus and name, and the score,
	whether or not it reaches the threshold."""

	path = os.path.join(hh.MOTIFS, motif + '.motif')
	m = hh.read_motifs(path)[0]
	got = mscore_lines(hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path,
		'-mscore'], tmp_path).stdout)
	assert sorted(got) == sorted(name for name, _ in targets)
	for name, seq in targets:
		score, pos, strand = hh.best_site(seq, m)
		g = got[name]
		assert g[1:4] == (m.consensus, m.name, strand), name
		assert g[0] == pos, name
		assert g[4] == pytest.approx(score, abs=1e-6), name


def test_mscore_with_threads(tmp_path, targets_tsv):
	path = os.path.join(hh.MOTIFS, 'siteB.motif')
	one = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path, '-mscore'], tmp_path)
	four = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path, '-mscore', '-p', 4],
		tmp_path)
	assert mscore_lines(four.stdout) == mscore_lines(one.stdout)


def test_sites_at_both_ends_and_case(tmp_path):
	"""Windows at the first and last base are scanned, and lowercase bases
	score as uppercase ones."""

	site = 'TAGCGTACCGTA'
	records = [
		('first', site + 'T' * 30),
		('last', 'T' * 30 + hh.revcomp(site)),
		('lower', 'T' * 10 + site.lower() + 'T' * 10),
		('whole', site),
		('short', site[:-1]),
	]
	seqs = hh.write_seq_tsv(str(tmp_path / 'ends.tsv'), records)
	path = os.path.join(hh.MOTIFS, 'siteA.motif')
	got = find_lines(hh.run(['homer2', 'find', '-s', seqs, '-m', path, '-offset', 0],
		tmp_path).stdout)
	assert [(r[0], r[1], r[4]) for r in got] == [('first', 0, '+'), ('last', 41, '-'),
		('lower', 10, '+'), ('whole', 0, '+')]
	assert_same_sites(got, expected_lines(records, hh.read_motifs(path), offset=0))


def test_windows_with_other_characters_are_skipped(tmp_path):
	"""N and any other character outside ACGT make a window unscorable."""

	site = 'TAGCGTACCGTA'
	records = [
		('n_inside', 'TTTT' + site[:5] + 'N' + site[6:] + 'TTTT'),
		('iupac_inside', 'TTTT' + site[:5] + 'R' + site[6:] + 'TTTT'),
		('n_outside', 'NNNN' + site + 'NNNN'),
	]
	seqs = hh.write_seq_tsv(str(tmp_path / 'n.tsv'), records)
	got = find_lines(hh.run(['homer2', 'find', '-s', seqs, '-m',
		os.path.join(hh.MOTIFS, 'siteA.motif'), '-offset', 0], tmp_path).stdout)
	assert [(r[0], r[1]) for r in got] == [('n_outside', 4)]


def test_crlf_and_wrapped_fasta(tmp_path, targets):
	"""The same sequences give the same sites whether the FASTA file wraps
	lines at 60 or 17 bases, or ends lines with CR LF."""

	path = os.path.join(hh.MOTIFS, 'siteA_siteB.motifs')
	sixty = hh.write_fasta(str(tmp_path / 'sixty.fa'), targets[:60])
	short = hh.write_fasta(str(tmp_path / 'short.fa'), targets[:60], width=17)
	crlf = tmp_path / 'crlf.fa'
	crlf.write_bytes(open(sixty, 'rb').read().replace(b'\n', b'\r\n'))
	results = [find_lines(hh.run(['homer2', 'find', '-i', f, '-m', path], tmp_path).stdout)
		for f in (sixty, short, str(crlf))]
	assert results[0] == results[1] == results[2]
	assert_same_sites(results[0], expected_lines(targets[:60], hh.read_motifs(path)))


def test_duplicate_fasta_names(tmp_path):
	"""pins the current output: a repeated FASTA name gets "-2", "-3", ...
	appended in the order the records appear."""

	site = 'TAGCGTACCGTA'
	fasta = hh.write_fasta(str(tmp_path / 'dup.fa'), [('x', 'TT' + site + 'TT'),
		('x', 'TTTT' + site), ('x', site + 'TTTT'), ('y', 'TTTT' + site)])
	got = find_lines(hh.run(['homer2', 'find', '-i', fasta, '-m',
		os.path.join(hh.MOTIFS, 'siteA.motif'), '-offset', 0], tmp_path).stdout)
	assert [(r[0], r[1]) for r in got] == [('x', 2), ('x-2', 4), ('x-3', 0), ('y', 4)]


def test_library_motif_finds_planted_rest_sites(tmp_path, targets, targets_tsv):
	"""REST-NRSF from HOMER's own library (data/knownTFs) finds the 60
	planted copies of its consensus, and nothing the reference scan does not
	find."""

	library = hh.read_motifs(os.path.join(hh.HOMER_HOME, 'data', 'knownTFs',
		'vertebrates', 'known.motifs'))
	rest = [m for m in library if m.name.startswith('REST-NRSF')]
	assert len(rest) == 1
	path = hh.write_motifs(str(tmp_path / 'rest.motif'), rest)
	got = find_lines(hh.run(['homer2', 'find', '-s', targets_tsv, '-m', path], tmp_path).stdout)
	assert_same_sites(got, expected_lines(targets, rest))
	planted = {p['peak'] for p in hh.read_planted() if p['motif'] == 'REST'}
	assert planted <= {row[0] for row in got}


def test_header_columns_after_the_threshold_are_ignored(tmp_path, targets_tsv):
	"""A motif file as HOMER writes it, with log p-value and statistics after
	the threshold, finds the same sites as one with three header columns."""

	m = hh.read_motifs(os.path.join(hh.MOTIFS, 'siteA.motif'))[0]
	plain = hh.write_motifs(str(tmp_path / 'plain.motif'), [m])
	m.header = ['-123.456', '0', 'T:120.0(60.00%),B:0.0(0.00%),P:1e-53']
	rich = hh.write_motifs(str(tmp_path / 'rich.motif'), [m])
	a = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', plain], tmp_path)
	b = hh.run(['homer2', 'find', '-s', targets_tsv, '-m', rich], tmp_path)
	assert find_lines(a.stdout) == find_lines(b.stdout) != []


def test_probabilities_below_the_minimum_are_raised(tmp_path):
	"""A probability of 0 scores as 0.001: one mismatch at a zero costs
	ln(0.001 / 0.25), so a threshold just below that site's score keeps it."""

	consensus = 'ACGTACGT'
	rows = [[1.0 if b == base else 0.0 for b in 'ACGT'] for base in consensus]
	exact = 8 * hh.math.log(1.0 / 0.25)
	one_off = 7 * hh.math.log(1.0 / 0.25) + hh.math.log(0.001 / 0.25)
	motif = hh.Motif(consensus, 'zeros', one_off - 0.01, rows)
	path = hh.write_motifs(str(tmp_path / 'zeros.motif'), [motif])
	seqs = hh.write_seq_tsv(str(tmp_path / 's.tsv'), [('exact', 'GG' + consensus + 'GG'),
		('one_off', 'GG' + 'ACGTTCGT' + 'GG')])
	got = find_lines(hh.run(['homer2', 'find', '-s', seqs, '-m', path, '-offset', 0,
		'-strand', '+'], tmp_path).stdout)
	assert [(r[0], r[1]) for r in got] == [('exact', 2), ('one_off', 2)]
	assert got[0][5] == pytest.approx(exact, abs=1e-5)
	assert got[1][5] == pytest.approx(one_off, abs=1e-5)

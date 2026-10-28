"""homer2 mask: every sequence, with each motif's sites replaced by N.

Expected sequences come from homer_helpers.mask_sequence(), applied for one
motif after another in the order of the motif file: + strand windows from the
left, then - strand windows from the right, each site masked as soon as it is
found, so a window that overlaps a masked site holds an N and is not scored.
Sequences are written upper case, with any character other than A, C, G or T
as N. The output is a tab-separated sequence file whose line order is not
part of its meaning, so it is compared as a mapping from name to sequence.
"""

import itertools
import os

import pytest

import homer_helpers as hh

SITE_A = 'TAGCGTACCGTA'


def masked(text):
	rows = {}
	for line in text.splitlines():
		cols = line.split('\t')
		assert len(cols) == 2, line
		assert cols[0] not in rows
		rows[cols[0]] = cols[1]
	return rows


def expected(records, motifs, strand='both'):
	result = {}
	for name, seq in records:
		for motif in motifs:
			seq = hh.mask_sequence(seq, motif, strand)
		result[name] = hh.clean_sequence(seq)
	return result


def site_spans(seq, motif):
	spans = []
	for pos, _, strand, _ in hh.scan(seq, motif, offset=0):
		start = pos if strand == '+' else pos - len(motif) + 1
		spans.append((start, start + len(motif)))
	return spans


def sites_overlap(records, motifs):
	"""Whether a site of one motif overlaps a site of another in any
	sequence."""

	for _, seq in records:
		spans = [site_spans(seq, motif) for motif in motifs]
		for a, b in itertools.combinations(spans, 2):
			if any(x[0] < y[1] and y[0] < x[1] for x in a for y in b):
				return True
	return False


def motif_file(name):
	return os.path.join(hh.MOTIFS, name)


@pytest.fixture(scope='module')
def targets():
	return hh.read_fasta(hh.TARGETS)


@pytest.fixture(scope='module')
def targets_tsv(tmp_path_factory, targets):
	path = tmp_path_factory.mktemp('mask') / 'targets.tsv'
	return hh.write_seq_tsv(str(path), targets)


@pytest.fixture(scope='module')
def rest_motif(tmp_path_factory):
	library = hh.read_motifs(os.path.join(hh.HOMER_HOME, 'data', 'knownTFs', 'vertebrates',
		'known.motifs'))
	rest = [m for m in library if m.name.startswith('REST-NRSF')]
	assert len(rest) == 1
	return rest[0]


@pytest.mark.parametrize('motif', ['siteA', 'siteB', 'ebox'])
def test_one_motif(tmp_path, targets, targets_tsv, motif):
	path = motif_file(motif + '.motif')
	got = masked(hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', path], tmp_path).stdout)
	want = expected(targets, hh.read_motifs(path))
	assert got == want
	assert sum(1 for name, seq in targets if want[name] != hh.clean_sequence(seq)) > 0


def test_planted_sites_are_masked(tmp_path, targets, targets_tsv):
	"""Every planted siteA copy, lowercase ones included, becomes 12 Ns, and
	nothing else in a sequence changes besides case."""

	got = masked(hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', motif_file('siteA.motif')],
		tmp_path).stdout)
	sequences = dict(targets)
	peak_starts = {peak['id']: peak['start'] for peak in hh.read_peaks()}
	for p in hh.read_planted():
		if p['motif'] != 'siteA':
			continue
		start = int(p['start']) - peak_starts[p['peak']] + 1
		assert got[p['peak']][start:start + 12] == 'N' * 12, p
		original = hh.clean_sequence(sequences[p['peak']])
		assert original[start:start + 12] in (SITE_A, hh.revcomp(SITE_A))
	for name, seq in targets:
		original = hh.clean_sequence(seq)
		assert all(a == b or a == 'N' for a, b in zip(got[name], original))
		assert len(got[name]) == len(original)


def test_fasta_input_matches_tsv_input(tmp_path, targets_tsv):
	"""-i gives the same sequences as -s and leaves no temporary files
	behind."""

	path = motif_file('siteA_siteB.motifs')
	tsv = hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', path], tmp_path)
	fasta = hh.run(['homer2', 'mask', '-i', hh.TARGETS, '-m', path], tmp_path)
	assert masked(fasta.stdout) == masked(tsv.stdout)
	assert os.listdir(str(tmp_path)) == []


@pytest.mark.parametrize('motifs', ['siteA_siteB.motifs', 'known_small.motifs'])
def test_several_motifs_are_masked_in_file_order(tmp_path, targets, targets_tsv, motifs):
	"""known_small.motifs holds siteB before ebox, and some siteB and ebox
	sites overlap, so the order in which motifs are masked matters here."""

	path = motif_file(motifs)
	got = masked(hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', path], tmp_path).stdout)
	assert got == expected(targets, hh.read_motifs(path))


def test_overlapping_sites_of_two_motifs(tmp_path):
	"""The first motif in the file masks its site, so the second motif's
	site, which overlaps it, keeps the bases it does not share; listing the
	motifs the other way round masks the second motif's site instead."""

	site_a = hh.read_motifs(motif_file('siteA.motif'))[0]
	other = hh.motif_from_consensus('ACCGTAGGCA', 'other')
	seq = 'TTTT' + SITE_A + 'GGCATTTT'
	seqs = hh.write_seq_tsv(str(tmp_path / 's.tsv'), [('s', seq)])
	ab = hh.write_motifs(str(tmp_path / 'ab.motifs'), [site_a, other])
	ba = hh.write_motifs(str(tmp_path / 'ba.motifs'), [other, site_a])
	got_ab = masked(hh.run(['homer2', 'mask', '-s', seqs, '-m', ab], tmp_path).stdout)['s']
	got_ba = masked(hh.run(['homer2', 'mask', '-s', seqs, '-m', ba], tmp_path).stdout)['s']
	assert got_ab == 'TTTT' + 'N' * 12 + 'GGCATTTT'
	assert got_ba == 'TTTT' + 'TAGCGT' + 'N' * 10 + 'TTTT'
	assert got_ab == expected([('s', seq)], [site_a, other])['s']
	assert got_ba == expected([('s', seq)], [other, site_a])['s']


def test_library_motif(tmp_path, targets, targets_tsv, rest_motif):
	path = hh.write_motifs(str(tmp_path / 'rest.motif'), [rest_motif])
	got = masked(hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', path], tmp_path).stdout)
	assert got == expected(targets, [rest_motif])
	planted = {p['peak'] for p in hh.read_planted() if p['motif'] == 'REST'}
	assert all('N' * 20 in got[name] for name in planted)


def test_output_file(tmp_path, targets_tsv):
	path = motif_file('siteA.motif')
	out = tmp_path / 'masked.tsv'
	result = hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', path, '-o', out], tmp_path)
	assert result.stdout == ''
	stdout = hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', path], tmp_path).stdout
	assert masked(out.read_text()) == masked(stdout)


@pytest.mark.parametrize('strand', ['+', '-', 'both'])
def test_strand(tmp_path, targets, targets_tsv, strand):
	path = motif_file('siteA_siteB.motifs')
	got = masked(hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', path, '-strand', strand],
		tmp_path).stdout)
	assert got == expected(targets, hh.read_motifs(path), strand)


def test_one_strand_masks_fewer_sites(tmp_path, targets_tsv):
	"""Planted sites alternate between strands, so masking one strand
	leaves sites that masking both removes."""

	path = motif_file('siteA.motif')
	runs = {strand: masked(hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', path,
		'-strand', strand], tmp_path).stdout) for strand in ('+', '-', 'both')}
	for strand in ('+', '-'):
		count = sum(seq.count('N') for seq in runs[strand].values())
		assert 0 < count < sum(seq.count('N') for seq in runs['both'].values())


@pytest.mark.parametrize('threads', [2, 4])
def test_threads_give_the_same_sequences(tmp_path, targets, targets_tsv, rest_motif, threads):
	"""With motifs whose sites never overlap, the result does not depend on
	which motif is masked first, so any number of threads gives the same
	sequences."""

	motifs = [m for m in hh.read_motifs(motif_file('known_small.motifs')) if m.name != 'ebox']
	motifs.append(rest_motif)
	assert not sites_overlap(targets, motifs)
	path = hh.write_motifs(str(tmp_path / 'set.motifs'), motifs)
	one = hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', path, '-p', 1], tmp_path)
	many = hh.run(['homer2', 'mask', '-s', targets_tsv, '-m', path, '-p', threads], tmp_path)
	assert masked(many.stdout) == masked(one.stdout) == expected(targets, motifs)


def test_case_other_characters_and_ends(tmp_path):
	"""Sites at the first and last base and lowercase sites are masked;
	lowercase is written upper case and other characters as N; a palindrome
	is masked once, on the + strand."""

	records = [
		('first', SITE_A + 'CCCC'),
		('last', 'CCCC' + hh.revcomp(SITE_A)),
		('lower', 'cc' + SITE_A.lower() + 'cc'),
		('iupac', 'CCRYCC' + SITE_A),
		('ebox', 'TTCACGTGTT'),
		('none', 'acgtACGT'),
	]
	seqs = hh.write_seq_tsv(str(tmp_path / 's.tsv'), records)
	path = hh.write_motifs(str(tmp_path / 'm.motifs'),
		hh.read_motifs(motif_file('siteA.motif')) + hh.read_motifs(motif_file('ebox.motif')))
	got = masked(hh.run(['homer2', 'mask', '-s', seqs, '-m', path], tmp_path).stdout)
	assert got == {
		'first': 'N' * 12 + 'CCCC',
		'last': 'CCCC' + 'N' * 12,
		'lower': 'CC' + 'N' * 12 + 'CC',
		'iupac': 'CCNNCC' + 'N' * 12,
		'ebox': 'TT' + 'N' * 6 + 'TT',
		'none': 'ACGTACGT',
	}
	assert got == expected(records, hh.read_motifs(path))

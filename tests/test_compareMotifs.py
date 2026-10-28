"""compareMotifs.pl: similarity scores between motifs, the reduction of
near-duplicate motifs, matches against a library of known motifs, the
filters on motif statistics, and the HTML pages and motif files it writes.

Expected values come from a Python implementation of compareMotifs.pl's
scoring in this module (read from bin/compareMotifs.pl and Statistics.pm):

- a motif file is read as compareMotifs.pl reads it: each row of a motif
  whose sixth header column does not start with "T:" is divided by its sum
  (so count matrices become frequencies) and its consensus is rebuilt from
  the rows (the most likely base, N when it is below 0.4);
- two motifs are compared at every offset that leaves at least one column
  overlapping, on the given strand and on the reverse complement of the
  second motif; unaligned columns are padded with 0.25 rows so that both
  matrices have the same length, and the score of an alignment is the
  Pearson correlation of all values of the two padded matrices (the default
  -stat correlation). The best alignment wins, ties keep the first found
  (offsets from the first motif's last column down to minus the second
  motif's length plus one, the forward strand winning ties with the reverse);
- -stat absDifference sums p1 * p2 - 1/16 over the padded columns and bases,
  -stat freqError sums 1 - (squared error / the error expected by chance)
  per column; both are divided by the aligned length when reducing and by
  the motif's length when matching;
- reduction sorts motifs by their log p-value (stable), keeps each motif
  not yet merged and merges every later motif scoring above -reduceThresh
  into it;
- known matches are sorted by score (stable), the first being the best guess.

The reference reproduces the order of floating point operations, so scores
are compared as the text compareMotifs.pl prints.
"""

import math
import os
import random
import re

import pytest

import homer_helpers as hh

UNIFORM = [0.25, 0.25, 0.25, 0.25]


# --- reading motifs as compareMotifs.pl does ------------------------------------

def perl_number(x):
	"""Perl's default text for a number."""

	return '%.15g' % x


def consensus(rows):
	code = ''
	for row in rows:
		best, base = 0.0, 'A'
		for j, v in enumerate(row):
			if v > best:
				best, base = v, 'ACGT'[j]
		code += base if best >= 0.4 else 'N'
	return code


def parse_motifs(text):
	"""Motifs as dicts: header fields (text as written; trailing empty
	fields dropped as Perl's split does), 'raw' rows as text, 'rows' as the
	numbers compareMotifs.pl compares, 'cons' and 'homer2' (a sixth header
	column starting with "T:")."""

	motifs = []
	for line in text.splitlines():
		cols = line.replace('\r', '').split('\t')
		while cols and cols[-1] == '':
			cols.pop()
		if cols and cols[0].startswith('>'):
			header = cols + [None] * (6 - len(cols))
			motifs.append({'cons': cols[0][1:], 'name': header[1] or '', 'v': header[2] or '0',
				'logp': header[3] or '1', 'gapinfo': header[4] or '0', 'str': header[5] or '',
				'raw': []})
		else:
			motifs[-1]['raw'].append(cols)
	for m in motifs:
		values = [[float(v) for v in row] for row in m['raw']]
		m['homer2'] = m['str'].startswith('T:')
		if m['homer2']:
			m['rows'] = values
			m['text_rows'] = m['raw']
		else:
			m['cons'] = consensus(values)
			rows = []
			for row in values:
				total = 0.0
				for v in row:
					total += v
				rows.append([v / total for v in row] if total >= 1e-100 else list(UNIFORM))
			m['rows'] = rows
			m['text_rows'] = [[perl_number(v) for v in row] for row in rows]
			if m['name'] == '':
				m['name'] = m['cons']
	return motifs


def read_perl_motifs(path):
	with open(str(path)) as handle:
		return parse_motifs(handle.read())


def homer2_stats(m):
	"""(targets, target %, background, background %) of a "T:" statistics
	column."""

	found = re.match(r'T:([^(]*)\(([^%]*)%\),B:([^(]*)\(([^%]*)%\)', m['str'])
	return tuple(float(v) for v in found.groups())


def motif_file_text(motifs, best_guess=None):
	"""The motif file compareMotifs.pl writes: header columns consensus,
	name (with ",BestGuess:..." when given), threshold, log p-value, gap
	info and statistics, then the rows."""

	lines = []
	for m in motifs:
		name = m['name']
		if best_guess and m['name'] in best_guess:
			name += ',BestGuess:' + best_guess[m['name']]
		lines.append('>' + '\t'.join([m['cons'], name, m['v'], m['logp'], m['gapinfo'],
			m['str']]))
		lines += ['\t'.join(row) for row in m['text_rows']]
	return '\n'.join(lines) + '\n'


# --- the similarity score ------------------------------------------------------------

def correlation(x, y):
	"""Statistics::correlation: None where the variance is zero ("NA")."""

	n = len(x)
	xy = xs = ys = x2 = y2 = 0.0
	for a, b in zip(x, y):
		xy += a * b
		xs += a
		ys += b
		x2 += a * a
		y2 += b * b
	numerator = xy - xs * ys / n
	denominator = (x2 - xs * xs / n) * (y2 - ys * ys / n)
	if denominator <= 0:
		return None
	return numerator / math.sqrt(denominator)


def column_score(mm1, mm2, stat):
	if stat == 'correlation':
		return correlation([v for row in mm1 for v in row], [v for row in mm2 for v in row])
	score = 0.0
	if stat == 'absDifference':
		for r1, r2 in zip(mm1, mm2):
			for j in range(4):
				score += r1[j] * r2[j] - 0.0625
		return score
	for r1, r2 in zip(mm1, mm2):
		current = expected = 0.0
		for j in range(4):
			d = r1[j] - r2[j]
			current -= d * d
			for k in range(4):
				d2 = r1[j] - r2[k]
				expected -= d2 * d2 / 4
		score += -1 * (current - expected) / expected
	return score


def similarity(m1, m2, stat='correlation', norevopp=False):
	"""(score, offset, reverse, aligned length) of the best alignment of two
	row lists; offset is where m2's first column sits relative to m1's."""

	len1, len2 = len(m1), len(m2)
	rv2 = [row[::-1] for row in m2[::-1]]
	best = (-1e10, 0, 0, 0)
	for off in range(len1 - 1, -len2, -1):
		mm1 = [UNIFORM] * max(0, -off) + m1 + [UNIFORM] * max(0, off + len2 - len1)
		mm2 = [UNIFORM] * max(0, off) + m2
		mm2 += [UNIFORM] * (len(mm1) - len(mm2))
		mm2r = [UNIFORM] * max(0, off) + rv2
		mm2r += [UNIFORM] * (len(mm1) - len(mm2r))
		forward = column_score(mm1, mm2, stat) or 0.0
		reverse = column_score(mm1, mm2r, stat) or 0.0
		if not norevopp and reverse > forward:
			if reverse > best[0]:
				best = (reverse, off, 1, len(mm1))
		elif forward > best[0]:
			best = (forward, off, 0, len(mm1))
	return best


def reduction_score(m1, m2, stat='correlation', norevopp=False):
	score, _, _, length = similarity(m1['rows'], m2['rows'], stat, norevopp)
	return score if stat == 'correlation' else score / length


def reduce_motifs(motifs, threshold, stat='correlation', norevopp=False):
	"""[(kept motif, [(merged motif, score), ...]), ...] in output order."""

	order = sorted(motifs, key=lambda m: float(m['logp']))
	merged = set()
	kept = []
	for i, m in enumerate(order):
		if i in merged:
			continue
		similar = []
		for j in range(i + 1, len(order)):
			if j in merged:
				continue
			score = reduction_score(m, order[j], stat, norevopp)
			if score > threshold:
				merged.add(j)
				similar.append((order[j], score))
		kept.append((m, similar))
	return kept


def known_matches(motif, library, stat='correlation'):
	"""[(score, library motif, offset, reverse), ...] best first."""

	results = []
	for known in library:
		score, offset, reverse, _ = similarity(motif['rows'], known['rows'], stat)
		if stat != 'correlation':
			score /= len(motif['rows'])
		results.append((score, known, offset, reverse))
	results.sort(key=lambda r: -r[0])
	return results


def revopp_consensus(code):
	pairs = dict(zip('ACGTRYMKBVDH', 'TGCAYRKMVBHD'))
	return ''.join(pairs.get(c, c) for c in reversed(code))


def alignment(motif, known, offset, reverse):
	"""The two alignment lines of the info page: '-' where a motif has no
	column."""

	if reverse:
		if known['homer2']:
			second = revopp_consensus(known['cons'])
		else:
			second = consensus([row[::-1] for row in known['rows'][::-1]])
	else:
		second = known['cons']
	a1 = '-' * max(0, -offset) + motif['cons']
	a2 = '-' * max(0, offset) + second
	width = max(len(a1), len(a2))
	return a1.ljust(width, '-'), a2.ljust(width, '-')


# --- reading compareMotifs.pl output ------------------------------------------------------

def read_text(path):
	with open(str(path)) as handle:
		return handle.read()


def read_matrix(path):
	rows = hh.read_tsv(str(path))
	names = rows[0][1:]
	return names, {(row[0], names[j]): value for row in rows[1:]
		for j, value in enumerate(row[1:])}


INFO_MATCH = re.compile(r'<H4>(.*?)</H4>.*?Match Rank:</TD><TD>(\d+)</TD>.*?'
	r'Score:</TD><TD>([^<]*)</TD.*?Offset:</TD><TD>([^<]*)</TD.*?'
	r'Orientation:</TD><TD>(forward|reverse) strand</TD>.*?<FONT face="courier">([^<]*)'
	r'</FONT><BR><FONT face="courier">([^<]*)</FONT>', re.S)


def info_matches(path):
	"""(name, rank, score, offset, orientation, line 1, line 2) of each known
	match listed on a motif's info page."""

	return [m.groups() for m in INFO_MATCH.finditer(read_text(path))]


def main_page_best_guesses(path):
	return re.findall(r'<TD>([^<]*)<BR/><A target="_blank" HREF="homerResults/motif\d+\.info',
		read_text(path))


def similar_scores(path):
	text = read_text(path)
	text = text[text.index('Similar de novo motifs found'):]
	return re.findall(r'<TR><TD>(\d+)</TD><TD>([^<]*)</TD><TD><svg', text)


def svg_letters(svg):
	"""[[(letter, height), ...] per column] in drawing order."""

	found = re.findall(r'transform="matrix\(([^,]+),0,0,([^,]+),([^,]+),([^)]+)\)">(\w)</text>',
		svg)
	columns = [found[i:i + 4] for i in range(0, len(found), 4)]
	return [[(f[4], float(f[1])) for f in column] for column in columns]


def expected_letters(rows, bits=False, rna=False):
	alphabet = 'ACGU' if rna else 'ACGT'
	columns = []
	for row in rows:
		if bits:
			h = 0.0
			for v in row:
				v = max(v, 0.001)
				h += v * math.log(v) / math.log(2)
			h *= -1
			row = [v * (2 - (h + 0)) / 2 for v in row]
		order = sorted(range(4), key=lambda j: row[j])
		columns.append([(alphabet[j], row[j]) for j in order])
	return columns


def assert_letters(got, want):
	assert len(got) == len(want)
	for g, w in zip(got, want):
		assert [x[0] for x in g] == [x[0] for x in w]
		assert [x[1] for x in g] == pytest.approx([x[1] for x in w], rel=1e-9, abs=1e-12)


# --- inputs --------------------------------------------------------------------------------

def from_consensus(code, name, top=0.91):
	m = hh.motif_from_consensus(code, name, top=top)
	return m.text()


def reverse_complement_text(code, name, top=0.91):
	m = hh.motif_from_consensus(code, name, top=top)
	rows = [row[::-1] for row in m.rows[::-1]]
	return hh.Motif(hh.revcomp(code), name, m.threshold, rows).text()


COUNTS = [[10 if b == c else 1 for b in 'ACGT'] for c in 'TGACGTCA']

HAND = ''.join([
	from_consensus('TGACGTCATC', 'tga'),
	from_consensus('TGACGTCATC', 'tga_copy'),
	reverse_complement_text('TGACGTCATC', 'tga_rc'),
	from_consensus('GTCATCGGAA', 'tga_shift4'),
	from_consensus('AAAAAAAA', 'polyA', top=0.7),
	from_consensus('CCGGTTAA', 'ccgg'),
	'>TGACGTCA\tcounts\t6.0\n' + ''.join('\t'.join(str(v) for v in row) + '\n' for row in COUNTS),
	'>TGACGTCA\tcounts_freq\t6.0\n' + ''.join('\t'.join(repr(v / 13.0) for v in row) + '\n'
		for row in COUNTS),
	from_consensus('GCGGAAGT', 'ets', top=0.8),
	from_consensus('CACGTGAC', 'ebox_like', top=0.85),
])

# homer2 de novo output: motifs with "T:" statistics, listed in an order other
# than by log p-value; h_dup copies h_top's rows with a weaker p-value.
HOMER2 = '\n'.join([
	'>TAGCGTAC\th_dup\t9.0\t-50.0\t0\tT:110.0(55.00%),B:20.0(2.00%),P:1e-21',
	'0.030\t0.030\t0.030\t0.910', '0.910\t0.030\t0.030\t0.030', '0.030\t0.030\t0.910\t0.030',
	'0.030\t0.910\t0.030\t0.030', '0.030\t0.030\t0.910\t0.030', '0.030\t0.030\t0.030\t0.910',
	'0.910\t0.030\t0.030\t0.030', '0.030\t0.910\t0.030\t0.030',
	'>TAGCGTAC\th_top\t9.0\t-100.0\t0\tT:150.0(75.00%),B:30.0(3.00%),P:1e-43\t'
	'Tpos:100.1,Tstd:20.0,Bpos:90.0,Bstd:50.0,StrandBias:0.1,Multiplicity:1.05',
	'0.030\t0.030\t0.030\t0.910', '0.910\t0.030\t0.030\t0.030', '0.030\t0.030\t0.910\t0.030',
	'0.030\t0.910\t0.030\t0.030', '0.030\t0.030\t0.910\t0.030', '0.030\t0.030\t0.030\t0.910',
	'0.910\t0.030\t0.030\t0.030', '0.030\t0.910\t0.030\t0.030',
	'>GGATCCTT\th_mid\t8.0\t-30.0\t0\tT:60.0(30.00%),B:100.0(10.00%),P:1e-13',
	'0.050\t0.050\t0.850\t0.050', '0.050\t0.050\t0.850\t0.050', '0.850\t0.050\t0.050\t0.050',
	'0.050\t0.050\t0.050\t0.850', '0.050\t0.850\t0.050\t0.050', '0.050\t0.850\t0.050\t0.050',
	'0.050\t0.050\t0.050\t0.850', '0.050\t0.050\t0.050\t0.850',
	'>CTTTAAGG\th_rare\t8.0\t-12.0\t0\tT:20.0(10.00%),B:5.0(0.50%),P:1e-5',
	'0.100\t0.700\t0.100\t0.100', '0.100\t0.100\t0.100\t0.700', '0.100\t0.100\t0.100\t0.700',
	'0.100\t0.100\t0.100\t0.700', '0.700\t0.100\t0.100\t0.100', '0.700\t0.100\t0.100\t0.100',
	'0.100\t0.100\t0.700\t0.100', '0.100\t0.100\t0.700\t0.100',
	'>NNSNNNSN\th_flat\t5.0\t-25.0\t0\tT:90.0(45.00%),B:300.0(30.00%),P:1e-11',
	'0.250\t0.300\t0.250\t0.200', '0.300\t0.200\t0.250\t0.250', '0.200\t0.300\t0.300\t0.200',
	'0.250\t0.250\t0.200\t0.300', '0.300\t0.250\t0.250\t0.200', '0.200\t0.300\t0.300\t0.200',
	'0.250\t0.250\t0.300\t0.200', '0.200\t0.250\t0.300\t0.250',
]) + '\n'


def random_motifs(count=40, seed=11):
	"""`count` motifs of 6-14 columns: random rows with three decimals, every
	fifth one a noisy copy of the motif before it, every seventh a count
	matrix."""

	rng = random.Random(seed)
	texts, rows_by_index = [], []
	for i in range(count):
		if i % 5 == 4:
			rows = [[max(0.001, v + rng.uniform(-0.05, 0.05)) for v in row]
				for row in rows_by_index[-1]]
		else:
			rows = []
			for _ in range(rng.randint(6, 14)):
				weights = [rng.random() ** 3 for _ in range(4)]
				rows.append([w / sum(weights) for w in weights])
		rows_by_index.append(rows)
		if i % 7 == 6:
			lines = ['\t'.join(str(int(round(v * 50))) for v in row) for row in rows]
		else:
			lines = ['\t'.join('%.3f' % v for v in row) for row in rows]
		texts.append('>{}\tr{:02d}\t8.0\n'.format(consensus(rows), i) + '\n'.join(lines) + '\n')
	return ''.join(texts)


LIBRARY_QUERY = 'REST-NRSF(Zf)/Jurkat-NRSF-ChIP-Seq/Homer'


def library_entry(path, name):
	"""The lines of one motif of a library file, copied as they are."""

	lines, inside = [], False
	with open(path) as handle:
		for line in handle:
			if line.startswith('>'):
				inside = line.split('\t')[1] == name
			if inside:
				lines.append(line)
	return ''.join(lines)


DEFAULT_LIBRARY = os.path.join(hh.HOMER_HOME, 'data', 'knownTFs', 'all', 'all.motifs')
RNA_LIBRARY = os.path.join(hh.HOMER_HOME, 'data', 'knownTFs', 'all.rna.motifs')


@pytest.fixture(scope='module')
def inputs(tmp_path_factory):
	d = tmp_path_factory.mktemp('cmp_inputs')
	files = {
		'hand': HAND,
		'homer2': HOMER2,
		'many': random_motifs(),
		'library_query': library_entry(DEFAULT_LIBRARY, LIBRARY_QUERY),
		'one': from_consensus('TGACGTCATC', 'tga'),
	}
	paths = {}
	for key, text in files.items():
		paths[key] = str(d / (key + '.motifs'))
		with open(paths[key], 'w') as out:
			out.write(text)
	# compareMotifs.pl -rna replaces a -known file whose name contains
	# "all.motifs", so the small library gets a name that does not.
	paths['known'] = str(d / 'known.motifs')
	with open(paths['known'], 'w') as out:
		out.write(read_text(os.path.join(hh.MOTIFS, 'known_small.motifs')))
	return paths


CALLS = {
	'matrix': ('hand', ['-matrix', 'm.txt']),
	'matrix_norevopp': ('hand', ['-matrix', 'm.txt', '-norevopp']),
	'matrix_absDifference': ('hand', ['-matrix', 'm.txt', '-stat', 'absDifference']),
	'matrix_freqError': ('hand', ['-matrix', 'm.txt', '-stat', 'freqError']),
	'matrix_many': ('many', ['-matrix', 'm.txt']),
	'html': ('hand', ['-known', 'KNOWN']),
	'html_nofacts': ('hand', ['-known', 'KNOWN', '-nofacts']),
	'basic': ('hand', ['-known', 'KNOWN', '-basic']),
	'bits': ('hand', ['-known', 'KNOWN', '-bits']),
	'rna': ('hand', ['-known', 'KNOWN', '-rna']),
	'top2': ('hand', ['-known', 'KNOWN', '-matchThresh', 'T2']),
	'above03': ('hand', ['-known', 'KNOWN', '-matchThresh', 0.3]),
	'absDifference': ('hand', ['-known', 'KNOWN', '-stat', 'absDifference']),
	'reduce_03': ('hand', ['-reduceMotifs', 'r.motifs', '-reduceThresh', 0.3]),
	'reduce_09': ('hand', ['-reduceMotifs', 'r.motifs', '-reduceThresh', 0.9]),
	'reduce_default': ('hand', ['-reduceMotifs', 'r.motifs']),
	'reduce_norevopp': ('hand', ['-reduceMotifs', 'r.motifs', '-reduceThresh', 0.9,
		'-norevopp']),
	'reduce_freqError': ('hand', ['-reduceMotifs', 'r.motifs', '-stat', 'freqError',
		'-reduceThresh', 0.5]),
	'many_reduce': ('many', ['-reduceMotifs', 'r.motifs', '-reduceThresh', 0.5]),
	'many_cpu1': ('many', ['-known', 'KNOWN', '-reduceThresh', 0.5]),
	'many_cpu4': ('many', ['-known', 'KNOWN', '-reduceThresh', 0.5, '-cpu', 4]),
	'homer2': ('homer2', ['-known', 'KNOWN']),
	'filter_none': ('homer2', ['-reduceMotifs', 'r.motifs']),
	'filter_pvalue': ('homer2', ['-reduceMotifs', 'r.motifs', '-pvalue', 1e-20]),
	'filter_fold': ('homer2', ['-reduceMotifs', 'r.motifs', '-F', 5]),
	'filter_background': ('homer2', ['-reduceMotifs', 'r.motifs', '-b', 0.01]),
	'filter_targets': ('homer2', ['-reduceMotifs', 'r.motifs', '-minT', 50]),
	'filter_all': ('homer2', ['-reduceMotifs', 'r.motifs', '-pvalue', 1e-200]),
	'dbview': ('hand', ['-known', 'KNOWN', '-dbview', '-dbTable', 'db.txt']),
	'library': ('library_query', []),
	'rna_library': ('one', ['-rna']),
}


@pytest.fixture(scope='module')
def runs(tmp_path_factory, inputs):
	"""Every compareMotifs.pl call of this module, run once: key -> the
	directory it ran in (the output directory is 'out' inside it)."""

	base = tmp_path_factory.mktemp('cmp_runs')
	calls = {}
	for key, (motifs, options) in CALLS.items():
		cwd = base / key
		cwd.mkdir()
		if key == 'dbview':
			# readDBTable: column 2 is the motif name, column 12 a gene name
			with open(str(cwd / 'db.txt'), 'w') as out:
				out.write('\t'.join(['x', 'tga'] + [''] * 9 + ['GENEA']) + '\n')
				out.write('\t'.join(['x', 'tga'] + [''] * 9 + ['GENEB']) + '\n')
		options = [inputs['known'] if o == 'KNOWN' else o for o in options]
		calls[key] = (['compareMotifs.pl', inputs[motifs], 'out'] + options, cwd)
	# the -cpu 4 run is made on its own
	results = hh.run_many({k: v for k, v in calls.items() if k != 'many_cpu4'})
	results.update(hh.run_many({'many_cpu4': calls['many_cpu4']}, jobs=1))
	return {key: base / key for key in results}


@pytest.fixture(scope='module')
def hand(inputs):
	return read_perl_motifs(inputs['hand'])


@pytest.fixture(scope='module')
def known(inputs):
	return read_perl_motifs(inputs['known'])


def by_name(motifs):
	return {m['name']: m for m in motifs}


# --- the similarity matrix ------------------------------------------------------------------

def expected_matrix(motifs, stat='correlation', norevopp=False):
	"""-matrix: motifs sorted by log p-value, 1 on the diagonal, and the
	score of motif i against motif j (i < j) in both cells."""

	order = sorted(motifs, key=lambda m: float(m['logp']))
	names = [m['name'] for m in order]
	cells = {}
	for i, a in enumerate(order):
		cells[(a['name'], a['name'])] = '1'
		for j in range(i + 1, len(order)):
			b = order[j]
			text = perl_number(reduction_score(a, b, stat, norevopp))
			cells[(a['name'], b['name'])] = cells[(b['name'], a['name'])] = text
	return names, cells


@pytest.mark.parametrize('key,stat,norevopp', [
	('matrix', 'correlation', False),
	('matrix_norevopp', 'correlation', True),
	('matrix_absDifference', 'absDifference', False),
	('matrix_freqError', 'freqError', False),
])
def test_matrix_scores_match_the_reference(runs, hand, key, stat, norevopp):
	names, cells = read_matrix(runs[key] / 'm.txt')
	assert (names, cells) == expected_matrix(hand, stat, norevopp)


def test_matrix_writes_no_html(runs):
	"""-matrix computes the scores, writes them and quits."""

	assert (runs['matrix'] / 'm.txt').exists()
	assert not (runs['matrix'] / 'out').exists()


def test_identical_and_reverse_complement_motifs_score_one(runs):
	_, cells = read_matrix(runs['matrix'] / 'm.txt')
	assert float(cells[('tga', 'tga_copy')]) == pytest.approx(1.0, abs=1e-12)
	assert float(cells[('tga', 'tga_rc')]) == pytest.approx(1.0, abs=1e-12)
	assert float(cells[('tga_copy', 'tga_rc')]) == pytest.approx(1.0, abs=1e-12)


def test_norevopp_compares_the_given_strand_only(runs, hand):
	"""The reverse complement no longer matches perfectly; its score is the
	best forward-strand alignment."""

	_, cells = read_matrix(runs['matrix_norevopp'] / 'm.txt')
	m = by_name(hand)
	assert float(cells[('tga', 'tga_rc')]) < 0.9
	assert cells[('tga', 'tga_rc')] == perl_number(similarity(m['tga']['rows'],
		m['tga_rc']['rows'], norevopp=True)[0])
	assert float(cells[('tga', 'tga_copy')]) == pytest.approx(1.0, abs=1e-12)


def test_partial_overlap_and_unrelated_motifs(runs, hand):
	"""A motif sharing six of its ten columns scores between unrelated and
	identical; the reference places the overlap at offset 4 on the forward
	strand."""

	_, cells = read_matrix(runs['matrix'] / 'm.txt')
	m = by_name(hand)
	score, offset, reverse, length = similarity(m['tga']['rows'], m['tga_shift4']['rows'])
	assert (offset, reverse, length) == (4, 0, 14)
	assert 0.5 < float(cells[('tga', 'tga_shift4')]) < 0.7
	for other in ('polyA', 'ccgg'):
		assert float(cells[('tga', other)]) < 0.3


def test_count_matrix_scores_as_its_frequencies(runs):
	"""Rows that do not sum to 1 are divided by their sum first: the count
	matrix and the same matrix as frequencies score alike against every
	motif."""

	names, cells = read_matrix(runs['matrix'] / 'm.txt')
	assert float(cells[('counts', 'counts_freq')]) == pytest.approx(1.0, abs=1e-12)
	for name in names:
		if name not in ('counts', 'counts_freq'):
			assert float(cells[('counts', name)]) == pytest.approx(
				float(cells[('counts_freq', name)]), abs=1e-12)


def test_many_motifs_matrix(runs, inputs):
	"""Forty motifs of different lengths, count matrices among them: all
	780 pairs as the reference scores them."""

	names, cells = read_matrix(runs['matrix_many'] / 'm.txt')
	want = expected_matrix(read_perl_motifs(inputs['many']))
	assert len(names) == 40
	assert (names, cells) == want


# --- reduction ------------------------------------------------------------------------------

@pytest.mark.parametrize('key,threshold,stat,norevopp', [
	('reduce_default', 0.6, 'correlation', False),
	('reduce_03', 0.3, 'correlation', False),
	('reduce_09', 0.9, 'correlation', False),
	('reduce_norevopp', 0.9, 'correlation', True),
	('reduce_freqError', 0.5, 'freqError', False),
])
def test_reduce_motifs_file(runs, hand, key, threshold, stat, norevopp):
	"""-reduceMotifs <file> writes the motifs left after merging those
	scoring above -reduceThresh (default 0.6), as compareMotifs.pl prints
	them (rows divided by their sums), and quits."""

	kept = reduce_motifs(hand, threshold, stat, norevopp)
	assert read_text(runs[key] / 'r.motifs') == motif_file_text([m for m, _ in kept])
	assert not (runs[key] / 'out').exists()


def test_reduce_thresholds_merge_more_when_lower(runs):
	counts = {key: len(read_perl_motifs(runs[key] / 'r.motifs'))
		for key in ('reduce_03', 'reduce_default', 'reduce_09')}
	assert counts['reduce_03'] < counts['reduce_default'] < counts['reduce_09']
	names = {m['name'] for m in read_perl_motifs(runs['reduce_09'] / 'r.motifs')}
	assert 'tga_copy' not in names and 'tga_rc' not in names
	names = {m['name'] for m in read_perl_motifs(runs['reduce_norevopp'] / 'r.motifs')}
	assert 'tga_copy' not in names and 'tga_rc' in names


def test_many_motifs_reduce_merges_noisy_copies(runs, inputs):
	motifs = read_perl_motifs(inputs['many'])
	kept = reduce_motifs(motifs, 0.5)
	assert read_text(runs['many_reduce'] / 'r.motifs') == motif_file_text([m for m, _ in kept])
	merged = {s['name'] for _, similar in kept for s, _ in similar}
	assert {'r04', 'r09', 'r14'} <= merged


def test_non_redundant_motifs_and_similar_pages(runs, hand, known):
	"""The main page lists the motifs left after reduction; each one's
	similar page lists the motifs merged into it with their scores (three
	decimals), and nonRedundant.motifs holds them with their best guess."""

	out = runs['html'] / 'out'
	kept = reduce_motifs(hand, 0.6)
	guesses = {}
	for n, (m, similar) in enumerate(kept, 1):
		best = known_matches(m, known)[0]
		guesses[m['name']] = '{}({:.3f})'.format(best[1]['name'], best[0])
		page = out / 'homerResults' / 'motif{}.similar.html'.format(n)
		assert similar_scores(page) == [(str(k), '%.3f' % s) for k, (_, s) in
			enumerate(similar, 1)]
		for k, (s, _) in enumerate(similar, 1):
			text = read_text(out / 'homerResults' / 'motif{}.similar{}.motif'.format(n, k))
			assert text == motif_file_text([s])
		assert read_text(out / 'homerResults' / 'motif{}.motif'.format(n)) == \
			motif_file_text([m], guesses)
	assert not (out / 'homerResults' / 'motif{}.motif'.format(len(kept) + 1)).exists()
	assert main_page_best_guesses(out / 'homerResults.html') == [guesses[m['name']]
		for m, _ in kept]
	assert read_text(out / 'nonRedundant.motifs') == motif_file_text([m for m, _ in kept],
		guesses)


def test_reverse_opposite_motif_file(runs, hand):
	"""motif<N>RV.motif holds the rows in reverse order with A/T and C/G
	swapped."""

	out = runs['html'] / 'out' / 'homerResults'
	m = reduce_motifs(hand, 0.6)[0][0]
	rv = read_perl_motifs(out / 'motif1RV.motif')[0]
	assert rv['raw'] == [[perl_number(v) for v in row[::-1]] for row in m['rows'][::-1]]


# --- known motifs ---------------------------------------------------------------------------

def expected_info_matches(motif, library, top=10, threshold=None, stat='correlation'):
	rows = []
	for rank, (score, k, offset, reverse) in enumerate(known_matches(motif, library, stat), 1):
		if threshold is None and rank > top:
			break
		if threshold is not None and score < threshold:
			break
		a1, a2 = alignment(motif, k, offset, reverse)
		rows.append((k['name'], str(rank), '%.2f' % score, str(offset),
			'reverse' if reverse else 'forward', a1, a2))
	return rows


@pytest.mark.parametrize('key,top,threshold,stat', [
	('html', 10, None, 'correlation'),
	('top2', 2, None, 'correlation'),
	('above03', None, 0.3, 'correlation'),
	('absDifference', 10, None, 'absDifference'),
])
def test_known_matches_on_info_pages(runs, hand, known, key, top, threshold, stat):
	"""Each motif's info page lists the known motifs best first: the top 10
	by default, the top # with -matchThresh T#, or all scoring at least # with
	-matchThresh #; with score (two decimals), offset, strand and alignment."""

	out = runs[key] / 'out' / 'homerResults'
	kept = reduce_motifs(hand, 0.6, stat)
	for n, (m, _) in enumerate(kept, 1):
		got = info_matches(out / 'motif{}.info.html'.format(n))
		assert got == expected_info_matches(m, known, top, threshold, stat), m['name']


def test_library_motif_best_match_is_itself(runs):
	"""A motif copied from the default library (data/knownTFs/all/all.motifs)
	matches itself best, with score 1, and the info page lists the
	reference's top 10."""

	out = runs['library'] / 'out'
	library = read_perl_motifs(DEFAULT_LIBRARY)
	query = [m for m in library if m['name'] == LIBRARY_QUERY]
	assert len(query) == 1
	want = expected_info_matches(query[0], library)
	assert want[0][:5] == (LIBRARY_QUERY, '1', '1.00', '0', 'forward')
	assert info_matches(out / 'homerResults' / 'motif1.info.html') == want
	assert main_page_best_guesses(out / 'homerResults.html') == [LIBRARY_QUERY + '(1.000)']


def test_rna_uses_the_rna_library_and_u(runs, inputs):
	"""-rna compares to data/knownTFs/all.rna.motifs when no -known is given,
	and logos show U for T."""

	out = runs['rna_library'] / 'out'
	motif = read_perl_motifs(inputs['one'])[0]
	best = known_matches(motif, read_perl_motifs(RNA_LIBRARY))[0]
	assert main_page_best_guesses(out / 'homerResults.html') == [
		'{}({:.3f})'.format(best[1]['name'], best[0])]
	svg = read_text(out / 'homerResults' / 'motif1.logo.svg')
	assert_letters(svg_letters(svg), expected_letters(motif['rows'], rna=True))


def test_rna_keeps_a_given_library(runs, hand, known):
	out = runs['rna'] / 'out'
	kept = reduce_motifs(hand, 0.6)
	got = info_matches(out / 'homerResults' / 'motif1.info.html')
	assert got == expected_info_matches(kept[0][0], known)
	assert '>T</text>' not in read_text(out / 'homerResults.html')
	assert '>U</text>' in read_text(out / 'homerResults.html')


def test_basic_skips_known_matches_and_similar_pages(runs, hand):
	out = runs['basic'] / 'out'
	kept = reduce_motifs(hand, 0.6)
	assert read_text(out / 'nonRedundant.motifs') == motif_file_text([m for m, _ in kept])
	for n in range(1, len(kept) + 1):
		assert not (out / 'homerResults' / 'motif{}.similar.html'.format(n)).exists()
		info = read_text(out / 'homerResults' / 'motif{}.info.html'.format(n))
		assert 'Matches to Known Motifs' not in info
	assert 'Similar Motifs Found' not in read_text(out / 'homerResults.html')


# --- logos, facts, threads -------------------------------------------------------------------

def test_logos_show_probabilities(runs, hand):
	"""Each column's letters are drawn smallest first, each as tall as its
	probability."""

	out = runs['html'] / 'out' / 'homerResults'
	for n, (m, _) in enumerate(reduce_motifs(hand, 0.6), 1):
		svg = read_text(out / 'motif{}.logo.svg'.format(n))
		assert_letters(svg_letters(svg), expected_letters(m['rows']))


def test_bits_scales_logos_by_information(runs, hand):
	"""-bits: heights are p * (2 - H) / 2, H the column's entropy in bits."""

	out = runs['bits'] / 'out' / 'homerResults'
	for n, (m, _) in enumerate(reduce_motifs(hand, 0.6), 1):
		svg = read_text(out / 'motif{}.logo.svg'.format(n))
		assert_letters(svg_letters(svg), expected_letters(m['rows'], bits=True))


def test_nofacts_leaves_fact_paragraphs_empty(runs):
	out = runs['html_nofacts'] / 'out'
	pages = [out / 'homerResults.html'] + sorted((out / 'homerResults').glob('*.html'))
	for page in pages:
		assert re.findall(r'<P>(.*?)</P>', read_text(page), re.S) == ['']


def output_files(directory):
	files = {}
	for root, _, names in os.walk(str(directory)):
		for name in names:
			path = os.path.join(root, name)
			files[os.path.relpath(path, str(directory))] = read_text(path)
	return files


def test_cpu_gives_the_same_output(runs):
	"""-cpu 4 writes the same files as -cpu 1."""

	a = output_files(runs['many_cpu1'] / 'out')
	b = output_files(runs['many_cpu4'] / 'out')
	assert len(a) > 100
	assert sorted(a) == sorted(b)
	for name in a:
		assert a[name] == b[name], name


# --- homer2 statistics and filters ---------------------------------------------------------

def test_homer2_statistics_on_pages(runs, inputs):
	"""Motifs with "T:" statistics: the totals are the counts divided by
	their percentages, rounded, and kept motifs are sorted by p-value, with
	the more significant of two duplicates kept."""

	out = runs['homer2'] / 'out'
	page = read_text(out / 'homerResults.html')
	assert 'Total target sequences = 200<BR/>' in page
	assert 'Total background sequences = 1000<BR/>' in page
	kept = reduce_motifs(read_perl_motifs(inputs['homer2']), 0.6)
	assert [m['name'] for m, _ in kept][:2] == ['h_top', 'h_mid']
	assert [s['name'] for s, _ in kept[0][1]] == ['h_dup']
	text = read_text(out / 'nonRedundant.motifs')
	assert [line.split('\t')[1].split(',')[0] for line in text.splitlines()
		if line.startswith('>')] == [m['name'] for m, _ in kept]
	info = read_text(out / 'homerResults' / 'motif1.info.html')
	assert '<TR><TD>Number of Target Sequences with motif</TD><TD>150.0</TD></TR>' in info
	assert '<TR><TD>Average Position of motif in Targets</TD><TD>100.1 +/- 20.0bp</TD></TR>' in info


def keep_names(motifs, pvalue=1e10, fold=-1e10, background=-1, min_targets=0):
	"""Names of motifs passing compareMotifs.pl's filters: log p-value below
	ln(-pvalue), fold enrichment (target % / background %) above -F,
	background % at least 100 * -b and target count at least -minT."""

	names = []
	for m in motifs:
		targets, tp, back, bp = homer2_stats(m)
		ratio = tp / bp if bp > 0 else 10000
		if float(m['logp']) < math.log(pvalue) and ratio > fold and bp >= background * 100 \
				and targets >= min_targets:
			names.append(m['name'])
	return names


@pytest.mark.parametrize('key,kwargs', [
	('filter_none', {}),
	('filter_pvalue', {'pvalue': 1e-20}),
	('filter_fold', {'fold': 5}),
	('filter_background', {'background': 0.01}),
	('filter_targets', {'min_targets': 50}),
])
def test_filters(runs, inputs, key, kwargs):
	motifs = read_perl_motifs(inputs['homer2'])
	names = keep_names(motifs, **kwargs)
	kept = reduce_motifs([m for m in motifs if m['name'] in names], 0.6)
	got = read_text(runs[key] / 'r.motifs')
	assert got == motif_file_text([m for m, _ in kept])
	if kwargs:
		assert len(names) < len(motifs)


def test_filtering_out_every_motif_writes_nothing(runs):
	assert not (runs['filter_all'] / 'r.motifs').exists()
	assert not (runs['filter_all'] / 'out').exists()


def test_homer2_rows_are_kept_as_written(runs, inputs):
	"""Rows of "T:" motifs are not rescaled: -reduceMotifs prints them as
	they were read."""

	text = read_text(runs['filter_none'] / 'r.motifs')
	assert '0.030\t0.030\t0.030\t0.910\n' in text


# --- internal options --------------------------------------------------------------------------

def test_dbview_keeps_every_motif(runs, hand):
	"""-dbview (internal) skips the reduction; -dbTable copies a motif's info
	page to <gene>.html for each gene the table gives it."""

	out = runs['dbview'] / 'out'
	page = read_text(out / 'homerResults.html')
	assert '<H1>Homer Motifs</H1>' in page
	assert len(re.findall(r'HREF="homerResults/motif\d+\.info\.html"', page)) == len(hand)
	first = read_text(out / 'homerResults' / 'motif1.info.html')
	assert 'Information for tga (Motif 1)' in first
	for gene in ('GENEA', 'GENEB'):
		assert read_text(out / 'homerResults' / (gene + '.html')) == first

"""homerTools cluster: hierarchical clustering of a distance matrix (-d) or
of the rows of a data matrix (-i), written as a TreeView .gtr/.cdt pair;
clusters cut from a .gtr file at a threshold (-gtr, -thresh); and scores of
each region against lists of regions (-c).

The matrices are small and clustered by hand: groups whose members are much
closer to each other than to anything else, so the merge order and the leaf
order are clear. The .gtr file lists merges in the order they happen, as
NODE<k>X <left> <right> <similarity>, leaves named GENE<row>X (0-based input
row), the similarity being minus the distance (or the correlation with -i).
Where a merge joins two groups, its distance depends on the linkage rule,
which the usage does not give; those values come from a Python
average-linkage reference (the average of all member-to-member distances),
read from Clustering.cpp.
"""

import math

import pytest

import homer_helpers as hh

FLT_MAX = 3.4028234663852886e+38

NAMES = ['a', 'b', 'c', 'd', 'e']
# a-b 1 and c-d 2 apart, the two pairs 10 apart, e 20 from everything
DIST = [
	[0, 1, 10, 10, 20],
	[1, 0, 10, 10, 20],
	[10, 10, 0, 2, 20],
	[10, 10, 2, 0, 20],
	[20, 20, 20, 20, 0],
]
TREE = [
	('NODE1X', 'GENE0X', 'GENE1X', -1.0),
	('NODE2X', 'GENE2X', 'GENE3X', -2.0),
	('NODE3X', 'NODE1X', 'NODE2X', -10.0),
	('NODE4X', 'NODE3X', 'GENE4X', -20.0),
]


def write_matrix(path, names, rows, header=None, ann=()):
	"""A tab-separated matrix: a header line, then per row its name, any
	extra annotation columns, and the values."""

	with open(str(path), 'w') as out:
		out.write('\t'.join(header or (['name'] + [a for a, _ in ann] + list(names))) + '\n')
		for i, (name, row) in enumerate(zip(names, rows)):
			out.write('\t'.join([name] + [v[i] for _, v in ann] + [repr(x) for x in row]) + '\n')
	return str(path)


def read_gtr(path):
	rows = []
	for line in open(str(path)).read().splitlines():
		cols = line.split('\t')
		assert len(cols) == 4, line
		rows.append((cols[0], cols[1], cols[2], float(cols[3])))
	return rows


def assert_same_tree(got, want):
	assert [r[:3] for r in got] == [r[:3] for r in want]
	for g, w in zip(got, want):
		assert g[3] == pytest.approx(w[3], abs=2e-6), (g, w)


def read_cdt(path):
	"""(header, [(gid, name, values)]) after checking the fixed columns and
	the EWEIGHT line."""

	lines = [line.split('\t') for line in open(str(path)).read().splitlines()]
	assert lines[0][:4] == ['GID', 'NA', 'NAME', 'GWEIGHT']
	n = len(lines[0]) - 4
	assert lines[1] == ['EWEIGHT', '', '', ''] + ['1.0'] * n
	rows = []
	for cols in lines[2:]:
		assert len(cols) == n + 4 and cols[1] == cols[2] and cols[3] == '1.0', cols
		rows.append((cols[0], cols[1], [float(v) for v in cols[4:]]))
	return lines[0][4:], rows


def average_linkage(dist):
	"""The .gtr rows HOMER's clustering gives for distance matrix `dist`:
	join the closest pair (first in row order on ties), the joined group's
	distance to every other being the size-weighted average of its two
	parts' distances, a merge never closer than the merges below it. Also
	returns the leaf order (left branch first; a group's left branch is the
	part holding the lower row)."""

	n = len(dist)
	m = [[float(x) for x in row] for row in dist]
	mask = [0] * n
	nodes = []
	for step in range(n - 1):
		best, left, right = FLT_MAX, 0, 0
		for i in range(n):
			if mask[i] == 1:
				continue
			for j in range(i + 1, n):
				if mask[j] != 1 and m[i][j] < best:
					best, left, right = m[i][j], i, j
		distance, weights, refs = best, [], []
		for k in (left, right):
			if mask[k] < 0:
				node = -mask[k] - 1
				weights.append(nodes[node][3])
				distance = max(distance, nodes[node][2])
				refs.append(mask[k])
			else:
				weights.append(1.0)
				refs.append(k)
		weight = weights[0] + weights[1]
		new = [(m[left][j] * weights[0] + m[right][j] * weights[1]) / weight for j in range(n)]
		mask[left] = -(step + 1)
		mask[right] = 1
		for j in range(n):
			m[left][j] = new[j]
			m[j][left] = new[j]
		nodes.append((refs[0], refs[1], distance, weight))

	def label(ref):
		return 'NODE{}X'.format(-ref) if ref < 0 else 'GENE{}X'.format(ref)

	rows, last = [], FLT_MAX
	for k, (left, right, distance, _) in enumerate(nodes):
		value = min(-distance, last)
		rows.append(('NODE{}X'.format(k + 1), label(left), label(right), value))
		last = value

	order = []

	def walk(ref):
		if ref >= 0:
			order.append(ref)
		else:
			walk(nodes[-ref - 1][0])
			walk(nodes[-ref - 1][1])

	walk(-(n - 1))
	return rows, order


def pearson(a, b):
	n = len(a)
	x1 = y1 = x2 = y2 = xy = 0.0
	for x, y in zip(a, b):
		x1 += x
		y1 += y
		x2 += x * x
		y2 += y * y
		xy += x * y
	den = (n * x2 - x1 * x1) * (n * y2 - y1 * y1)
	if den <= 0:
		return -2.0
	return (n * xy - x1 * y1) / math.sqrt(den)


@pytest.fixture(scope='module')
def dist_run(tmp_path_factory):
	"""-d on DIST with the default prefix, run once for several tests."""

	directory = tmp_path_factory.mktemp('cluster_d')
	write_matrix(directory / 'dist.txt', NAMES, DIST)
	result = hh.run(['homerTools', 'cluster', '-d', 'dist.txt'], directory)
	return directory, result


def test_distance_matrix_tree(dist_run):
	"""The pairs merge first (a-b at 1, c-d at 2), then the two pairs (all
	four cross distances are 10), then e (20 from all)."""

	directory, result = dist_run
	assert result.stdout == ''
	assert_same_tree(read_gtr(directory / 'out.gtr'), TREE)
	rows, order = average_linkage(DIST)
	assert_same_tree(rows, TREE)
	assert order == [0, 1, 2, 3, 4]


def test_distance_matrix_cdt(dist_run):
	"""out.cdt holds the input matrix with rows and columns in tree order."""

	directory, _ = dist_run
	header, rows = read_cdt(directory / 'out.cdt')
	assert header == NAMES
	assert [(r[0], r[1]) for r in rows] == [('GENE{}X'.format(i), n) for i, n in enumerate(NAMES)]
	assert [r[2] for r in rows] == [[float(x) for x in row] for row in DIST]


def shuffled():
	"""DIST with rows and columns in the order e, c, a, d, b."""

	perm = [4, 2, 0, 3, 1]
	return [NAMES[p] for p in perm], [[DIST[p][q] for q in perm] for p in perm], perm


def test_leaf_order_follows_the_tree(tmp_path):
	"""With the rows shuffled the groups still come out together in the
	.cdt; GENE ids are the shuffled input rows."""

	names, dist, perm = shuffled()
	write_matrix(tmp_path / 'dist.txt', names, dist)
	hh.run(['homerTools', 'cluster', '-d', 'dist.txt', '-o', 'shuf'], tmp_path)
	rows, order = average_linkage(dist)
	assert_same_tree(read_gtr(tmp_path / 'shuf.gtr'), rows)
	assert [r[3] for r in rows] == [-1.0, -2.0, -10.0, -20.0]
	header, cdt = read_cdt(tmp_path / 'shuf.cdt')
	assert header == [names[i] for i in order]
	assert [r[1] for r in cdt] == [names[i] for i in order]
	assert [r[1] for r in cdt] in (['e', 'c', 'd', 'a', 'b'], ['e', 'a', 'b', 'c', 'd'],
		['c', 'd', 'a', 'b', 'e'], ['a', 'b', 'c', 'd', 'e'])
	assert [r[2] for r in cdt] == [[float(dist[i][j]) for j in order] for i in order]


def test_output_prefix(tmp_path):
	write_matrix(tmp_path / 'dist.txt', NAMES, DIST)
	(tmp_path / 'sub').mkdir()
	hh.run(['homerTools', 'cluster', '-d', 'dist.txt', '-o', 'sub/my'], tmp_path)
	assert sorted(p.name for p in (tmp_path / 'sub').iterdir()) == ['my.cdt', 'my.gtr']
	assert not (tmp_path / 'out.gtr').exists()
	assert_same_tree(read_gtr(tmp_path / 'sub' / 'my.gtr'), TREE)


def test_annotation_columns(tmp_path):
	"""-annCols 2: a second annotation column before the values is skipped."""

	write_matrix(tmp_path / 'dist.txt', NAMES, DIST,
		ann=[('description', ['x{}'.format(i) for i in range(5)])])
	hh.run(['homerTools', 'cluster', '-d', 'dist.txt', '-annCols', 2], tmp_path)
	assert_same_tree(read_gtr(tmp_path / 'out.gtr'), TREE)
	header, rows = read_cdt(tmp_path / 'out.cdt')
	assert [r[1] for r in rows] == NAMES
	assert [r[2] for r in rows] == [[float(x) for x in row] for row in DIST]


def test_reverse_for_similarities(tmp_path):
	"""-r treats high values as close: a similarity matrix 30 - DIST gives
	the same tree, the .gtr listing the similarities."""

	sim = [[30 - x for x in row] for row in DIST]
	write_matrix(tmp_path / 'sim.txt', NAMES, sim)
	hh.run(['homerTools', 'cluster', '-d', 'sim.txt', '-r'], tmp_path)
	assert_same_tree(read_gtr(tmp_path / 'out.gtr'),
		[r[:3] + (30 + r[3],) for r in TREE])


def test_average_linkage(tmp_path):
	"""pins the current output: the distance between two groups is the
	average of their members' distances (UPGMA). Here a-b join at 1; c is 2
	from a but 8 from b (average 5), while c-d are 4 apart, so c-d join
	second, then everything at the averages."""

	names = ['a', 'b', 'c', 'd']
	dist = [
		[0, 1, 2, 9],
		[1, 0, 8, 7],
		[2, 8, 0, 4],
		[9, 7, 4, 0],
	]
	write_matrix(tmp_path / 'dist.txt', names, dist)
	hh.run(['homerTools', 'cluster', '-d', 'dist.txt'], tmp_path)
	rows, order = average_linkage(dist)
	assert rows == [('NODE1X', 'GENE0X', 'GENE1X', -1.0), ('NODE2X', 'GENE2X', 'GENE3X', -4.0),
		('NODE3X', 'NODE1X', 'NODE2X', -6.5)]
	assert_same_tree(read_gtr(tmp_path / 'out.gtr'), rows)


DATA_NAMES = ['up1', 'up2', 'down1', 'down2', 'flat', 'up3']
DATA = [
	[1.0, 2.0, 3.0, 4.0, 5.0],
	[2.0, 4.1, 6.0, 7.9, 10.2],
	[9.0, 7.0, 5.5, 3.0, 1.0],
	[4.0, 3.0, 2.0, 1.5, 0.0],
	[3.0, 2.5, 3.5, 2.0, 3.0],
	[0.5, 1.0, 1.4, 2.1, 2.4],
]
EXPS = ['t0', 't1', 't2', 't3', 't4']


@pytest.fixture(scope='module')
def data_file(tmp_path_factory):
	directory = tmp_path_factory.mktemp('cluster_i')
	return write_matrix(directory / 'data.txt', DATA_NAMES, DATA, header=['gene'] + EXPS)


@pytest.mark.parametrize('options', [[], ['-dist', 'correlation'], ['-dist', 'corr']])
def test_correlation_of_data_rows(tmp_path, data_file, options):
	"""-i clusters rows by Pearson correlation (the default): the rising
	rows, then the falling rows, join first; the .gtr lists correlations
	and the .cdt the data in tree order, columns as in the input."""

	hh.run(['homerTools', 'cluster', '-i', data_file] + options, tmp_path)
	dist = [[-1.0 if i == j else -pearson(a, b) for j, b in enumerate(DATA)]
		for i, a in enumerate(DATA)]
	rows, order = average_linkage(dist)
	got = read_gtr(tmp_path / 'out.gtr')
	assert_same_tree(got, rows)
	merged = [set(r[1:3]) for r in got[:3]]
	assert {'GENE0X', 'GENE1X'} in merged or {'GENE0X', 'GENE5X'} in merged
	assert {'GENE2X', 'GENE3X'} in merged
	assert all(r[3] > 0.95 for r in got[:3])
	header, cdt = read_cdt(tmp_path / 'out.cdt')
	assert header == EXPS
	assert [r[1] for r in cdt] == [DATA_NAMES[i] for i in order]
	assert [r[2] for r in cdt] == [DATA[i] for i in order]


def test_euclidean_distance_of_data_rows(tmp_path, data_file):
	hh.run(['homerTools', 'cluster', '-i', data_file, '-dist', 'euclidean'], tmp_path)
	dist = [[math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b))) for b in DATA] for a in DATA]
	rows, order = average_linkage(dist)
	assert_same_tree(read_gtr(tmp_path / 'out.gtr'), rows)
	assert [r[3] for r in rows] == sorted((r[3] for r in rows), reverse=True)
	_, cdt = read_cdt(tmp_path / 'out.cdt')
	assert [r[1] for r in cdt] == [DATA_NAMES[i] for i in order]


@pytest.mark.parametrize('sub', [5, 6, 7000])
def test_sub_not_below_the_row_count(tmp_path, data_file, sub):
	"""-sub <#> at 90% of the rows or more (5 of 6 rounds to the 90% mark)
	clusters all rows: the output equals the run without -sub."""

	hh.run(['homerTools', 'cluster', '-i', data_file, '-sub', sub, '-o', 'sub'], tmp_path)
	hh.run(['homerTools', 'cluster', '-i', data_file, '-o', 'all'], tmp_path)
	for suffix in ('.gtr', '.cdt'):
		assert (tmp_path / ('sub' + suffix)).read_text() == (tmp_path / ('all' + suffix)).read_text()


def clusters_text(text):
	return [tuple(line.split('\t')) for line in text.splitlines()]


@pytest.mark.parametrize('thresh, want', [
	(None, ['-1', '-1', '-1', '-1', '-1']),
	(1.5, ['1', '1', '-1', '-1', '-1']),
	(5, ['1', '1', '2', '2', '-1']),
	(15, ['1', '1', '1', '1', '-1']),
	(25, ['1', '1', '1', '1', '1']),
])
def test_clusters_from_gtr(dist_run, tmp_path, thresh, want):
	"""-gtr <file> -thresh <#> prints each leaf (in .cdt order, read from the
	.cdt beside the .gtr) with a cluster number: leaves under a merge closer
	than # share a number, counted from 1 in tree order; others get -1. The
	default threshold is 0."""

	directory, _ = dist_run
	options = [] if thresh is None else ['-thresh', thresh]
	result = hh.run(['homerTools', 'cluster', '-gtr', directory / 'out.gtr'] + options, tmp_path)
	assert clusters_text(result.stdout) == list(zip(NAMES, want))


def test_clusters_from_correlation_gtr(tmp_path, data_file):
	"""For a correlation tree the threshold is on minus the correlation:
	-thresh -0.95 groups leaves joined at a correlation above 0.95."""

	hh.run(['homerTools', 'cluster', '-i', data_file], tmp_path)
	got = clusters_text(hh.run(['homerTools', 'cluster', '-gtr', 'out.gtr', '-thresh', -0.95],
		tmp_path).stdout)
	groups = {}
	for name, number in got:
		groups.setdefault(number, set()).add(name)
	assert sorted(map(sorted, groups.values())) == [['down1', 'down2'], ['flat'],
		['up1', 'up2', 'up3']]
	assert dict(got)['flat'] == '-1'


def write_lists(directory):
	(directory / 'pairs.txt').write_text('a\tfirst\nb\nnot_a_region\n')
	(directory / 'far.txt').write_text('e\n')
	(directory / 'mid.txt').write_text('c\textra\td\n')


def expected_scores(names, dist, lists, sign=1):
	rows = []
	for i, name in enumerate(names):
		row = [name]
		for members in lists:
			idx = [names.index(m) for m in members]
			row.append(sign * sum(dist[i][k] for k in idx) / len(idx))
		rows.append(row)
	return rows


def read_scores(text):
	lines = [line.split('\t') for line in text.splitlines()]
	return lines[0], [[cols[0]] + [float(v) for v in cols[1:]] for cols in lines[1:]]


@pytest.mark.parametrize('reverse', [False, True])
def test_score_regions_against_lists(tmp_path, reverse):
	"""-c <file> ... prints, for each region of the matrix in input order, its
	average distance to the regions named in the first column of each file
	(names not in the matrix are skipped); with -r the distances are
	negated."""

	write_matrix(tmp_path / 'dist.txt', NAMES, DIST)
	write_lists(tmp_path)
	files = ['pairs.txt', 'far.txt', 'mid.txt']
	result = hh.run(['homerTools', 'cluster', '-d', 'dist.txt', '-c'] + files +
		(['-r'] if reverse else []), tmp_path)
	header, rows = read_scores(result.stdout)
	assert header == ['Region Name'] + files
	want = expected_scores(NAMES, DIST, [['a', 'b'], ['e'], ['c']], -1 if reverse else 1)
	assert rows == want
	assert not (tmp_path / 'out.gtr').exists()


def test_score_regions_from_gtr(dist_run, tmp_path):
	"""-gtr with -c scores against the matrix stored in the .cdt."""

	directory, _ = dist_run
	write_lists(tmp_path)
	result = hh.run(['homerTools', 'cluster', '-gtr', directory / 'out.gtr', '-c', 'pairs.txt',
		'mid.txt'], tmp_path)
	header, rows = read_scores(result.stdout)
	assert header == ['Region Name', 'pairs.txt', 'mid.txt']
	assert rows == expected_scores(NAMES, DIST, [['a', 'b'], ['c']])

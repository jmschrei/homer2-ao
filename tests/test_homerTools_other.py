"""homerTools decontaminate, homerTools special uniqmap, and the homerTools
command dispatcher.

decontaminate: the tag directories are made with makeTagDirectory from BED
reads the tests write. With -frac f, each contaminant tag position (value v)
removes v * f * (contaminated total / contaminant total) tags from the
nearest same-strand tags of the contaminated directory, nearest first
(downstream first on a tie), no farther than -size bp; the expectations for
the hand-made directory are worked out by hand. Without -frac the fraction is
estimated: regions of -size bp around contaminant positions with at least
-min tags are compared (log2 of 1 + tags per 10 million in each directory),
the differences binned by 0.1, and the fraction read from that histogram; the
reference below repeats that procedure from the region counts, which the test
data make plain (stacks of reads at single positions).

special uniqmap: the expected bit files are packed in Python from the input
positions.
"""

import math
import os
import shutil
import struct

import pytest

import homer_helpers as hh


def write_bed(path, reads):
	"""reads: (chromosome, 1-based tag position, strand, count); a + strand
	read starts at the position, a - strand read ends there."""

	with open(str(path), 'w') as out:
		for chrom, pos, strand, count in reads:
			start = pos - 1 if strand == '+' else pos - 30
			for _ in range(count):
				out.write('{}\t{}\t{}\tr\t0\t{}\n'.format(chrom, start, start + 30, strand))
	return str(path)


def read_tags(directory):
	"""{(chromosome, position, strand): value text} over every *.tags.tsv."""

	tags = {}
	for name in sorted(os.listdir(str(directory))):
		if not name.endswith('.tags.tsv'):
			continue
		for line in open(os.path.join(str(directory), name)):
			cols = line.rstrip('\n').split('\t')
			key = (cols[1], int(cols[2]), '+' if cols[3] == '0' else '-')
			assert key not in tags
			tags[key] = cols[4]
	return tags


# Contaminated: 10 tags on chr1 near the contaminant tags, 90 filler tags on
# chr2. Contaminant: 7 tags on chr1, 3 on chr3. Total 100 and 10, so with
# -frac 0.1 each contaminant read removes one tag.
EXP_READS = [
	('chr1', 10000, '+', 1), ('chr1', 10010, '+', 1),   # contaminant at 10004
	('chr1', 19990, '+', 1), ('chr1', 20010, '+', 1),   # 20000: a tie
	('chr1', 30002, '-', 1), ('chr1', 30005, '+', 1),   # 30000: wrong strand nearer
	('chr1', 40100, '+', 1),                            # 40000: 100 bp away
	('chr1', 50000, '+', 2), ('chr1', 50003, '+', 1),   # 50001, 3 reads
] + [('chr2', 1000 * (i + 1), '+', 1) for i in range(90)]
CTL_READS = [
	('chr1', 10004, '+', 1), ('chr1', 20000, '+', 1), ('chr1', 30000, '+', 1),
	('chr1', 40000, '+', 1), ('chr1', 50001, '+', 3),
	('chr3', 5000, '+', 3),
]


def changed(tags, values):
	out = dict(tags)
	for (pos, strand), value in values.items():
		key = ('chr1', pos, strand)
		assert key in out
		out[key] = value
	return out


@pytest.fixture(scope='module')
def tagdirs(tmp_path_factory):
	"""The hand-made contaminated (exp) and contaminant (ctl) directories,
	and one in-place decontamination of a copy of exp, made once here because
	it counts down 5 s before rewriting."""

	directory = tmp_path_factory.mktemp('decontaminate')
	write_bed(directory / 'exp.bed', EXP_READS)
	write_bed(directory / 'ctl.bed', CTL_READS)
	hh.run_many({
		'exp': (['makeTagDirectory', 'exp', 'exp.bed'], directory),
		'ctl': (['makeTagDirectory', 'ctl', 'ctl.bed'], directory),
	})
	shutil.copytree(str(directory / 'exp'), str(directory / 'inplace'))
	hh.run(['homerTools', 'decontaminate', 'inplace', 'ctl', '-frac', 0.1], directory)
	return directory


def test_tag_directories_as_written(tagdirs):
	exp = read_tags(tagdirs / 'exp')
	assert len(exp) == 99 and exp[('chr1', 50000, '+')] == '2.0'
	assert len(read_tags(tagdirs / 'ctl')) == 6


def test_frac_removes_nearest_tags(tagdirs, tmp_path):
	"""-frac 0.1 -o: 10004 removes the nearer 10000; 20000 sits 10 bp from
	both neighbours and removes the downstream one; 30000 skips the nearer
	- strand tag; 40000 reaches 40100 (within the default 250 bp); the
	3-read stack at 50001 empties 50000 (2 reads) and then 50003. Removed
	tags stay with value 0; the input directory is untouched."""

	exp = read_tags(tagdirs / 'exp')
	hh.run(['homerTools', 'decontaminate', tagdirs / 'exp', tagdirs / 'ctl', '-frac', 0.1,
		'-o', 'out'], tmp_path)
	got = read_tags(tmp_path / 'out')
	assert got == changed(exp, {(10000, '+'): '0.0', (20010, '+'): '0.0', (30005, '+'): '0.0',
		(40100, '+'): '0.0', (50000, '+'): '0.0', (50003, '+'): '0.0'})
	assert read_tags(tagdirs / 'exp') == exp
	info = (tmp_path / 'out' / 'tagInfo.txt').read_text().splitlines()[1].split('\t')
	assert float(info[2]) == 93.0


def test_in_place_without_o(tagdirs, tmp_path):
	"""Without -o the contaminated directory itself is rewritten."""

	hh.run(['homerTools', 'decontaminate', tagdirs / 'exp', tagdirs / 'ctl', '-frac', 0.1,
		'-o', 'out'], tmp_path)
	assert read_tags(tagdirs / 'inplace') == read_tags(tmp_path / 'out')
	assert read_tags(tagdirs / 'inplace') != read_tags(tagdirs / 'exp')


def test_frac_scales_the_removal(tagdirs, tmp_path):
	"""-frac 0.05 removes half a tag per contaminant read: 1.5 from the
	stack at 50000 (2 reads) leaves 0.5 and spares 50003."""

	exp = read_tags(tagdirs / 'exp')
	hh.run(['homerTools', 'decontaminate', tagdirs / 'exp', tagdirs / 'ctl', '-frac', 0.05,
		'-o', 'out'], tmp_path)
	assert read_tags(tmp_path / 'out') == changed(exp, {(10000, '+'): '0.5',
		(20010, '+'): '0.5', (30005, '+'): '0.5', (40100, '+'): '0.5', (50000, '+'): '0.5'})


def test_frac_zero_removes_nothing(tagdirs, tmp_path):
	hh.run(['homerTools', 'decontaminate', tagdirs / 'exp', tagdirs / 'ctl', '-frac', 0,
		'-o', 'out'], tmp_path)
	assert read_tags(tmp_path / 'out') == read_tags(tagdirs / 'exp')


@pytest.mark.parametrize('size, values', [
	(50, {(10000, '+'): '0.0', (20010, '+'): '0.0', (30005, '+'): '0.0', (50000, '+'): '0.0',
		(50003, '+'): '0.0'}),
	(3, {(50000, '+'): '0.0', (50003, '+'): '0.0'}),
])
def test_size_limits_the_distance(tagdirs, tmp_path, size, values):
	"""-size 50 spares 40100 (100 bp away); -size 3 spares everything but
	the stack at 50000 and 50003 (1 and 2 bp away)."""

	hh.run(['homerTools', 'decontaminate', tagdirs / 'exp', tagdirs / 'ctl', '-frac', 0.1,
		'-size', size, '-o', 'out'], tmp_path)
	assert read_tags(tmp_path / 'out') == changed(read_tags(tagdirs / 'exp'), values)


# Estimation: 30 contaminant stacks of 40 reads, 10 kb apart on chr1, and a
# stack of 10; the contaminated directory has 25 (first 20) or 60 (last 10)
# reads at each stack position, 5 more 100 bp downstream, and filler on chr2.
STACKS = [10000 * (k + 1) for k in range(30)]
EST_CTL = [('chr1', p, '+', 40) for p in STACKS] + [('chr1', 315000, '+', 10)]
EST_EXP = [('chr1', p, '+', 25 if k < 20 else 60) for k, p in enumerate(STACKS)] + \
	[('chr1', p + 100, '+', 5) for p in STACKS] + \
	[('chr2', 1000 + 60 * i, '+', 1) for i in range(750)]


@pytest.fixture(scope='module')
def estimate_dirs(tmp_path_factory):
	directory = tmp_path_factory.mktemp('estimate')
	write_bed(directory / 'exp.bed', EST_EXP)
	write_bed(directory / 'ctl.bed', EST_CTL)
	hh.run_many({
		'exp': (['makeTagDirectory', 'exp', 'exp.bed'], directory),
		'ctl': (['makeTagDirectory', 'ctl', 'ctl.bed'], directory),
	})
	return directory


def totals():
	return float(sum(r[3] for r in EST_EXP)), float(sum(r[3] for r in EST_CTL))


def region_counts(size, minimum):
	"""(contaminated, contaminant) tag counts of the regions found around
	contaminant stacks of at least `minimum` reads: the stack itself and, when
	the region reaches 100 bp downstream (-size above 200), the 5 reads
	there."""

	counts = []
	for k, p in enumerate(STACKS):
		if 40 < minimum:
			continue
		tc = (25 if k < 20 else 60) + (5 if size > 200 else 0)
		counts.append((tc, 40))
	return counts


def log_values(tc, ic):
	t, i = totals()
	return (math.log(1.0 + tc / t * 1e7) / math.log(2.0),
		math.log(1.0 + ic / i * 1e7) / math.log(2.0))


def estimate(counts, minimum):
	"""The histogram rows and the fraction estimateContamination derives:
	log differences (contaminant - contaminated) of regions with more than
	`minimum` contaminated tags, in 0.1 bins from -20 to 20 (truncated),
	smoothed over 3 bins; rows from the first non-empty bin; the estimate is
	1 / 2^d for the first bin d, scanning down from the top, whose smoothed
	count falls below the bin above once more than 10 regions are passed."""

	inc, half = 0.1, 200
	limit = 2 * half + 1
	hist = [0.0] * limit
	for tc, ic in counts:
		lt, li = log_values(tc, ic)
		if tc > minimum:
			index = min(max(int((li - lt) / inc) + half, 0), limit - 1)
			hist[index] += 1.0
	smooth = [0.0] * limit
	rows, started = [], False
	for i in range(limit):
		if hist[i] < 0.01 and not started:
			continue
		started = True
		near = [hist[i]] + ([hist[i - 1]] if i > 0 else []) + ([hist[i + 1]] if i + 1 < limit else [])
		smooth[i] = sum(near) / len(near)
		value = (i - half) * inc
		rows.append(['%.2f' % value, '%.3f' % min(1.0, 1.0 / 2.0 ** value), '%.1f' % hist[i],
			'%.1f' % smooth[i]])
	fraction, last, total = 0.0, 0.0, 0.0
	for i in range(limit - 1, -1, -1):
		total += smooth[i]
		if smooth[i] < last and total > 10.0:
			fraction = 1.0 / 2.0 ** ((i - half) * inc)
			break
		last = smooth[i]
	return rows, fraction


def float32(x):
	return struct.unpack('f', struct.pack('f', x))[0]


def remove(exp, ctl, fraction, size, norm):
	"""The removal pass on {(chr, pos, strand): value} dictionaries, float32
	values as in the tag files' Tag records."""

	out = dict(exp)
	for chrom in sorted({k[0] for k in exp}):
		tags = sorted((k[1], k[2]) for k in exp if k[0] == chrom)
		values = [float32(exp[(chrom, p, s)]) for p, s in tags]
		t = 0
		for p, s in sorted((k[1], k[2]) for k in ctl if k[0] == chrom):
			bad = float32(ctl[(chrom, p, s)] * norm * fraction)
			if bad <= 0.000001:
				continue
			while p > tags[t][0] and t < len(tags) - 1:
				t += 1
			down, up = t - 1, t
			while bad > 0.000001:
				dist_down = p - tags[down][0] if 0 <= down < len(tags) else size * 2
				dist_up = tags[up][0] - p if 0 <= up < len(tags) else size * 2
				k = up if dist_up <= dist_down else down
				if min(dist_up, dist_down) > size:
					break
				if tags[k][1] == s:
					if values[k] > bad:
						values[k] = float32(values[k] - bad)
						bad = 0.0
					else:
						bad = float32(bad - values[k])
						values[k] = 0.0
				if k == up:
					up += 1
				else:
					down -= 1
		for (q, s), v in zip(tags, values):
			out[(chrom, q, s)] = v
	return {k: '%.1f' % v for k, v in out.items()}


def scatter_values(path):
	lines = [line.split('\t') for line in open(str(path)).read().splitlines()]
	assert lines[0] == ['RegionID', 'chr', 'start', 'end', 'strand', 'Contaminated (log2)',
		'Contaminant (log2)']
	return sorted((float(r[5]), float(r[6])) for r in lines[1:]), lines[1:]


def assert_same_values(got, want):
	assert len(got) == len(want)
	for g, w in zip(got, want):
		assert g == pytest.approx(w, abs=2e-6)


@pytest.mark.parametrize('options, size, minimum', [
	([], 250, 20.0),
	(['-size', 50], 50, 20.0),
	(['-min', 35], 250, 35.0),
	(['-min', 45], 250, 45.0),
])
def test_estimate_only(estimate_dirs, tmp_path, options, size, minimum):
	"""-estimateOnly writes contaminationScatterPlot.txt (one row per region,
	each region containing its stack) and contaminationHistogram.txt into the
	output directory and leaves the tags alone. -size 50 leaves the reads
	100 bp downstream out of each region; -min 45 finds no region (stacks
	hold 40); -min 35 keeps the regions but bins only those with more than
	35 contaminated tags."""

	d = estimate_dirs
	hh.run(['homerTools', 'decontaminate', d / 'exp', d / 'ctl', '-estimateOnly', '-o', 'out']
		+ options, tmp_path)
	counts = region_counts(size, minimum)
	values, rows = scatter_values(tmp_path / 'out' / 'contaminationScatterPlot.txt')
	assert_same_values(values, sorted(log_values(tc, ic) for tc, ic in counts))
	stacks = sorted(STACKS)
	starts = sorted((int(r[2]), int(r[3])) for r in rows)
	for (start, end), p in zip(starts, stacks if counts else []):
		assert start <= p <= end and end - start <= size
	lines = [line.split('\t') for line in
		(tmp_path / 'out' / 'contaminationHistogram.txt').read_text().splitlines()]
	assert lines[0] == ['Log Difference (contaminant/contaminated)', 'Contamination Estimate',
		'Count(min %f)' % minimum, 'Smoothed']
	assert lines[1:] == estimate(counts, minimum)[0]
	assert read_tags(tmp_path / 'out') == read_tags(d / 'exp')


@pytest.mark.parametrize('options', [[], ['-frac', -1]])
def test_estimated_fraction_is_removed(estimate_dirs, tmp_path, options):
	"""Without -frac (or with a negative one) the estimated fraction is
	removed: here the estimate comes from the 20 stacks with 25 reads,
	emptying those (and the 5 reads beside them) and leaving part of each
	60-read stack."""

	d = estimate_dirs
	hh.run(['homerTools', 'decontaminate', d / 'exp', d / 'ctl', '-o', 'out'] + options, tmp_path)
	t, i = totals()
	_, fraction = estimate(region_counts(250, 20.0), 20.0)
	assert 0.0 < fraction < 1.0
	exp, ctl = read_tags(d / 'exp'), read_tags(d / 'ctl')
	want = remove({k: float(v) for k, v in exp.items()}, {k: float(v) for k, v in ctl.items()},
		fraction, 250, t / i)
	got = read_tags(tmp_path / 'out')
	assert got == want
	assert got[('chr1', STACKS[0], '+')] == '0.0'
	assert 0.0 < float(got[('chr1', STACKS[-1], '+')]) < 60.0


def test_uniqmap(tmp_path):
	"""special uniqmap reads chromosome (cut at the first '.'), position and
	strand (0 +, 1 -) per line, and writes <chr>.p.uniqmap and
	<chr>.n.uniqmap: the size (a multiple of 8 covering the last set
	position) as a 4-byte integer, then one bit per position, and
	uniqMapStats.txt with the number of distinct positions set and twice the
	size per chromosome and for the genome."""

	positions = {
		'chr1': [(0, 0), (7, 0), (8, 1), (100, 0), (100, 1), (7, 0)],
		'chrX': [(5, 1), (20, 0)],
	}
	with open(str(tmp_path / 'map.txt'), 'w') as out:
		for pos, strand in positions['chr1'][:4]:
			out.write('chr1\t{}\t{}\n'.format(pos, strand))
		out.write('chrX.part2\t5\t1\nshort\tline\n')
		for pos, strand in positions['chr1'][4:]:
			out.write('chr1.alt\t{}\t{}\textra\n'.format(pos, strand))
		out.write('chrX\t20\t0\n')
	result = hh.run(['homerTools', 'special', 'uniqmap', 'umap', 'map.txt'], tmp_path)
	assert result.stdout == ''
	stats = {}
	for chrom, entries in positions.items():
		top = max(pos for pos, _ in entries) // 8
		size = (top + 1) * 8
		for strand, letter in ((0, 'p'), (1, 'n')):
			bits = bytearray(top + 2)
			for pos, s in entries:
				if s == strand:
					bits[pos // 8] |= 1 << (pos % 8)
			data = (tmp_path / 'umap' / '{}.{}.uniqmap'.format(chrom, letter)).read_bytes()
			assert data == struct.pack('<I', size) + bytes(bits)
		stats[chrom] = [str(len(set(entries))), str(size)]
	lines = [line.split('\t') for line in (tmp_path / 'umap' / 'uniqMapStats.txt').read_text()
		.splitlines()]
	assert lines[0] == ['#Chr', 'Mappable bp (both strands)', 'Total bp (both strand)']
	assert {r[0]: r[1:] for r in lines[1:-1]} == stats
	assert lines[-1] == ['genome', '7', '128']
	assert sorted(os.listdir(str(tmp_path / 'umap'))) == ['chr1.n.uniqmap', 'chr1.p.uniqmap',
		'chrX.n.uniqmap', 'chrX.p.uniqmap', 'uniqMapStats.txt']


COMMANDS = ['barcodes', 'truseq', 'trim', 'freq', 'extract', 'decontaminate', 'cluster',
	'special']


def test_usage_lists_the_commands(tmp_path):
	result = hh.run(['homerTools'], tmp_path)
	assert result.stdout == ''
	for command in COMMANDS:
		assert command in result.stderr
	assert os.listdir(str(tmp_path)) == []


@pytest.mark.parametrize('args', [['--help'], ['nonsense']])
def test_help_and_unknown_command_print_the_usage(tmp_path, args):
	usage = hh.run(['homerTools'], tmp_path).stderr
	result = hh.run(['homerTools'] + args, tmp_path)
	assert result.stdout == ''
	assert result.stderr.endswith(usage)


@pytest.mark.parametrize('args', [[c] for c in COMMANDS] + [['special', 'uniqmap'],
	['special', 'uniqmap', 'umap']])
def test_command_without_arguments_prints_its_usage(tmp_path, args):
	"""Each command run without its arguments prints its own usage (not the
	general one), exits 0 and writes nothing."""

	usage = hh.run(['homerTools'], tmp_path).stderr
	result = hh.run(['homerTools'] + args, tmp_path)
	assert result.stdout == ''
	assert result.stderr.strip() and not result.stderr.endswith(usage)
	assert os.listdir(str(tmp_path)) == []

"""Shared helpers for HOMER's tests: where HOMER is, how to run its programs,
the synthetic inputs, and small reference implementations that expected
values are computed with.

Environment variables:

	HOMER_HOME      a HOMER install configured with `perl configureHomer.pl -make`
	                (default: the directory above tests/)
	HOMER_TEST_BIN  where the tests find the programs they call directly
	                (default: $HOMER_HOME/bin); HOMER's own scripts always run the
	                programs in $HOMER_HOME/bin
	HOMER_TEST_JOBS how many program calls a batch runs at once (default: the
	                number of CPUs, at most 8)
"""

import concurrent.futures
import math
import os
import subprocess

TESTS = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(TESTS, 'data')
MOTIFS = os.path.join(DATA, 'motifs')
PERL_DIR = os.path.join(TESTS, 'perl')

HOMER_HOME = os.path.abspath(os.environ.get('HOMER_HOME') or os.path.dirname(TESTS))
HOMER_BIN = os.path.join(HOMER_HOME, 'bin')
TEST_BIN = os.path.abspath(os.environ.get('HOMER_TEST_BIN') or HOMER_BIN)
JOBS = max(1, int(os.environ.get('HOMER_TEST_JOBS') or min(8, os.cpu_count() or 1)))

GENOME = os.path.join(DATA, 'genome.fa')
PEAKS = os.path.join(DATA, 'peaks.txt')
PEAKS_BED = os.path.join(DATA, 'peaks.bed')
TARGETS = os.path.join(DATA, 'targets.fa')
BACKGROUND = os.path.join(DATA, 'background.fa')
PLANTED = os.path.join(DATA, 'planted.tsv')

# Motif.cpp and Motif2.cpp raise every probability below this before taking
# the log (MOTIF_PROBABILITY_MINIMUM).
PROBABILITY_MINIMUM = 0.001

COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


# --- running HOMER -------------------------------------------------------------

def homer_env(tmpdir, extra=None):
	"""The environment every call gets. Perl's hash order and rand() are fixed:
	PERL_HASH_SEED and PERL_PERTURB_KEYS for the order, and HomerTestSrand.pm,
	loaded into every perl through PERL5OPT, for srand(1)."""

	env = {key: value for key, value in os.environ.items()
		if not key.startswith('PERL')}
	env.update({
		'PATH': HOMER_BIN + os.pathsep + os.environ.get('PATH', '/usr/bin:/bin'),
		'PERL_HASH_SEED': '0',
		'PERL_PERTURB_KEYS': '0',
		'PERL5LIB': PERL_DIR,
		'PERL5OPT': '-MHomerTestSrand',
		'LC_ALL': 'C',
		'TMPDIR': str(tmpdir),
	})
	env.pop('LANG', None)
	env.pop('LANGUAGE', None)
	if extra:
		env.update(extra)
	return env


def program(name):
	return os.path.join(TEST_BIN, name)


def run(args, cwd, check=True, timeout=600, env=None, stdin=None):
	"""Run HOMER program args[0] (a name in $HOMER_TEST_BIN) in `cwd` and
	return the CompletedProcess, stdout and stderr as text. With `check`, a
	non-zero exit fails the test and shows the end of stderr."""

	args = [str(a) for a in args]
	result = subprocess.run([program(args[0])] + args[1:], cwd=str(cwd),
		env=homer_env(cwd, env), stdin=stdin if stdin is not None else subprocess.DEVNULL,
		stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
	result.stdout = result.stdout.decode('utf-8', 'replace')
	result.stderr = result.stderr.decode('utf-8', 'replace')
	if check and result.returncode != 0:
		raise AssertionError('{} exited {}\n{}'.format(' '.join(args),
			result.returncode, result.stderr[-3000:]))
	return result


def run_many(calls, check=True, timeout=900, jobs=None):
	"""Run several calls at once, `jobs` at a time. `calls` maps a key to a
	(args, cwd) or (args, cwd, env) tuple; returns key -> CompletedProcess.
	Each call needs its own cwd."""

	def one(item):
		key, call = item
		args, cwd = call[0], call[1]
		env = call[2] if len(call) > 2 else None
		os.makedirs(str(cwd), exist_ok=True)
		return key, run(args, cwd, check=False, timeout=timeout, env=env)

	results = {}
	with concurrent.futures.ThreadPoolExecutor(max_workers=jobs or JOBS) as pool:
		for key, result in pool.map(one, sorted(calls.items())):
			results[key] = result
	if check:
		failed = ['{}: {} exited {}\n{}'.format(key, ' '.join(str(a) for a in
			calls[key][0]), r.returncode, r.stderr[-1500:])
			for key, r in sorted(results.items()) if r.returncode != 0]
		if failed:
			raise AssertionError('\n'.join(failed))
	return results


# --- files ---------------------------------------------------------------------

def revcomp(seq):
	return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
	"""[(name, sequence)] in file order; the name is the whole header line
	after '>'."""

	records = []
	name, parts = None, []
	with open(path) as handle:
		for line in handle:
			line = line.rstrip('\r\n')
			if line.startswith('>'):
				if name is not None:
					records.append((name, ''.join(parts)))
				name, parts = line[1:], []
			elif name is not None:
				parts.append(line.strip())
	if name is not None:
		records.append((name, ''.join(parts)))
	return records


def write_fasta(path, records, width=60):
	with open(path, 'w') as out:
		for name, seq in records:
			out.write('>{}\n'.format(name))
			for i in range(0, len(seq), width):
				out.write(seq[i:i + width] + '\n')
	return path


def write_seq_tsv(path, records):
	with open(path, 'w') as out:
		for name, seq in records:
			out.write('{}\t{}\n'.format(name, seq))
	return path


def write_group(path, rows):
	"""rows: (name, group, weight); group 1 is target, 0 background."""

	with open(path, 'w') as out:
		for name, group, weight in rows:
			out.write('{}\t{}\t{}\n'.format(name, group, weight))
	return path


def read_tsv(path, header=False):
	rows = []
	with open(path) as handle:
		for line in handle:
			line = line.rstrip('\r\n')
			if line:
				rows.append(line.split('\t'))
	if header:
		return rows[0], rows[1:]
	return rows


_GENOME_CACHE = {}


def read_genome(path=GENOME):
	"""{chromosome: sequence} with case kept; the name ends at the first
	whitespace, as HOMER reads FASTA headers."""

	if path not in _GENOME_CACHE:
		_GENOME_CACHE[path] = {name.split()[0]: seq for name, seq in read_fasta(path)}
	return _GENOME_CACHE[path]


def read_peaks(path=PEAKS):
	"""HOMER peak file rows as dicts: id, chr, start, end (1-based, both
	inclusive) and strand ('+' or '-')."""

	peaks = []
	for row in read_tsv(path):
		if row[0].startswith('#'):
			continue
		strand = row[4] if len(row) > 4 else '+'
		strand = '-' if strand in ('-', '1') else '+'
		peaks.append({'id': row[0], 'chr': row[1], 'start': int(row[2]),
			'end': int(row[3]), 'strand': strand})
	return peaks


def read_planted(path=PLANTED):
	header, rows = read_tsv(path, header=True)
	return [dict(zip(header, row)) for row in rows]


def resize(peak, size):
	"""(start, end), 1-based and inclusive, of the `size` window HOMER's
	resizePosFile.pl puts around a peak: centre c = (start + end) / 2, start
	floor(c - size/2), end floor(c + size/2), so the window holds size + 1
	bases."""

	centre = (peak['start'] + peak['end']) / 2.0
	return int(math.floor(centre - size / 2.0)), int(math.floor(centre + size / 2.0))


def peak_sequence(peak, genome=None, size=None, strand=True):
	"""The genome under a peak (1-based start and end, both inclusive), or
	under its resize() window; reverse-complemented for a - strand peak when
	`strand`."""

	genome = genome or read_genome()
	start, end = (peak['start'], peak['end']) if size is None else resize(peak, size)
	seq = genome[peak['chr']][start - 1:end]
	if strand and peak['strand'] == '-':
		seq = revcomp(seq)
	return seq


# --- motifs ----------------------------------------------------------------------

class Motif(object):
	"""A motif as HOMER's motif file holds it: the header's consensus, name and
	log-odds threshold, the rest of the header, and one row of A, C, G, T
	probabilities per position."""

	def __init__(self, consensus, name, threshold, rows, header=None):
		self.consensus = consensus
		self.name = name
		self.threshold = threshold
		self.rows = rows
		self.header = header or []

	def __len__(self):
		return len(self.rows)

	def log_odds(self):
		"""ln(p / 0.25) per position and base, probabilities raised to
		PROBABILITY_MINIMUM first: the score a site gets is the sum over its
		positions."""

		return [[math.log(max(p, PROBABILITY_MINIMUM) / 0.25) for p in row]
			for row in self.rows]

	def text(self):
		header = ['>' + self.consensus, self.name, '{:.6f}'.format(self.threshold)]
		header += self.header
		lines = ['\t'.join(header)]
		lines += ['\t'.join('{:.3f}'.format(p) for p in row) for row in self.rows]
		return '\n'.join(lines) + '\n'


def read_motifs(path):
	"""Every motif in a HOMER motif file, in order. A header line starts with
	'>' and holds consensus, name and threshold (tab-separated); each line
	after it holds the four probabilities of one position."""

	motifs = []
	with open(path) as handle:
		for line in handle:
			line = line.rstrip('\r\n')
			if not line.strip():
				continue
			cols = line.split('\t')
			if cols[0].startswith('>'):
				motifs.append(Motif(cols[0][1:], cols[1], float(cols[2]), [], cols[3:]))
			elif motifs:
				motifs[-1].rows.append([float(v) for v in cols[:4]])
	return motifs


def write_motifs(path, motifs):
	with open(path, 'w') as out:
		for motif in motifs:
			out.write(motif.text())
	return path


def motif_from_consensus(consensus, name, top=0.91, threshold=None):
	low = (1.0 - top) / 3.0
	rows = [[top if b == base else low for b in 'ACGT'] for base in consensus]
	if threshold is None:
		threshold = 0.8 * len(consensus) * math.log(top / 0.25)
	return Motif(consensus, name, threshold, rows)


def clean_sequence(seq):
	"""Upper case, and N for anything that is not A, C, G or T."""

	seq = seq.upper()
	return ''.join(base if base in 'ACGT' else 'N' for base in seq)


BASE_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def scan(seq, motif, strand='both', offset=None):
	"""Every site in `seq` scoring at least the motif's threshold, as
	(position, site, strand, score). A + strand site's position is its first
	base and a - strand site's its last, both 0-based, plus `offset`
	(default: minus half the sequence length, rounded towards zero, so that
	positions are relative to the middle); `site` is the + strand text; a
	window holding anything but A, C, G or T is not scored."""

	seq = clean_sequence(seq)
	if offset is None:
		offset = -(len(seq) // 2)
	lo = motif.log_odds()
	width = len(lo)
	sites = []
	if strand in ('both', '+'):
		for i in range(len(seq) - width + 1):
			window = seq[i:i + width]
			if 'N' in window:
				continue
			score = 0.0
			for j in range(width):
				score += lo[j][BASE_INDEX[window[j]]]
			if score >= motif.threshold:
				sites.append((i + offset, window, '+', score))
	if strand in ('both', '-'):
		for i in range(len(seq) - 1, width - 2, -1):
			window = seq[i - width + 1:i + 1]
			if 'N' in window:
				continue
			score = 0.0
			for j in range(width):
				score += lo[j][3 - BASE_INDEX[seq[i - j]]]
			if score >= motif.threshold:
				sites.append((i + offset, window, '-', score))
	return sites


def best_score(seq, motif, strand='both'):
	"""The highest score of any window, or None when no window can be
	scored."""

	seq = clean_sequence(seq)
	lo = motif.log_odds()
	width = len(lo)
	best = None
	for i in range(len(seq) - width + 1):
		window = seq[i:i + width]
		if 'N' in window:
			continue
		if strand in ('both', '+'):
			score = 0.0
			for j in range(width):
				score += lo[j][BASE_INDEX[window[j]]]
			best = score if best is None else max(best, score)
		if strand in ('both', '-'):
			score = 0.0
			for j in range(width):
				score += lo[j][3 - BASE_INDEX[window[width - 1 - j]]]
			best = score if best is None else max(best, score)
	return best


def best_site(seq, motif, strand='both', offset=None):
	"""(score, position, strand) of the best window, the first of equal
	scores in scan order (+ strand from the left, then - strand from the
	right), positions as scan() gives them; None when no window can be
	scored."""

	seq = clean_sequence(seq)
	if offset is None:
		offset = -(len(seq) // 2)
	lo = motif.log_odds()
	width = len(lo)
	best = None
	if strand in ('both', '+'):
		for i in range(len(seq) - width + 1):
			window = seq[i:i + width]
			if 'N' in window:
				continue
			score = 0.0
			for j in range(width):
				score += lo[j][BASE_INDEX[window[j]]]
			if best is None or score > best[0]:
				best = (score, i + offset, '+')
	if strand in ('both', '-'):
		for i in range(len(seq) - 1, width - 2, -1):
			window = seq[i - width + 1:i + 1]
			if 'N' in window:
				continue
			score = 0.0
			for j in range(width):
				score += lo[j][3 - BASE_INDEX[seq[i - j]]]
			if best is None or score > best[0]:
				best = (score, i + offset, '-')
	return best


def mask_sequence(seq, motif, strand='both'):
	"""`seq` with every site replaced by N, the way homer2 mask scans: + strand
	windows from the left, then - strand windows from the right, each site
	masked as soon as it is found, so a later window that overlaps it holds an
	N and is not scored."""

	s = list(clean_sequence(seq))
	lo = motif.log_odds()
	width = len(lo)
	if strand in ('both', '+'):
		for i in range(len(s) - width + 1):
			window = s[i:i + width]
			if 'N' in window:
				continue
			score = 0.0
			for j in range(width):
				score += lo[j][BASE_INDEX[window[j]]]
			if score >= motif.threshold:
				s[i:i + width] = ['N'] * width
	if strand in ('both', '-'):
		for i in range(len(s) - 1, width - 2, -1):
			window = s[i - width + 1:i + 1]
			if 'N' in window:
				continue
			score = 0.0
			for j in range(width):
				score += lo[j][3 - BASE_INDEX[s[i - j]]]
			if score >= motif.threshold:
				s[i - width + 1:i + 1] = ['N'] * width
	return ''.join(s)


def has_site(seq, motif, strand='both'):
	return bool(scan(seq, motif, strand))


# --- statistics --------------------------------------------------------------------

def log_binomial_tail(k, n, p):
	"""ln P(X >= k) for X ~ Binomial(n, p), summed in log space."""

	if k <= 0:
		return 0.0
	if p <= 0.0:
		return -math.inf
	if p >= 1.0:
		return 0.0
	terms = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
		+ i * math.log(p) + (n - i) * math.log1p(-p) for i in range(k, n + 1)]
	top = max(terms)
	return top + math.log(sum(math.exp(t - top) for t in terms))


def log_hypergeometric_tail(k, n, K, N):
	"""ln P(X >= k) for X ~ Hypergeometric(population N, K successes in it,
	n draws)."""

	def log_choose(a, b):
		return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)

	if k <= 0:
		return 0.0
	terms = [log_choose(K, i) + log_choose(N - K, n - i) - log_choose(N, n)
		for i in range(k, min(n, K) + 1) if n - i <= N - K]
	if not terms:
		return -math.inf
	top = max(terms)
	return top + math.log(sum(math.exp(t - top) for t in terms))


# --- HOMER output files ----------------------------------------------------------------

def read_known_results(path):
	"""knownResults.txt (findMotifsGenome.pl, findMotifs.pl) or homer2 known
	output as a list of dicts, one per motif in file order: name, consensus,
	pvalue (the text), logp, qvalue, targets, target_total, target_percent,
	background, background_total, background_percent. The totals come from the
	"(of N)" in the header."""

	header, rows = read_tsv(path, header=True)
	totals = []
	for column in (header[5], header[7]):
		inside = column[column.rfind('(of ') + 4:column.rfind(')')] if '(of ' in column else ''
		totals.append(float(inside) if inside else None)
	results = []
	for row in rows:
		results.append({
			'name': row[0], 'consensus': row[1], 'pvalue': row[2],
			'logp': float(row[3]), 'qvalue': float(row[4]),
			'targets': float(row[5]), 'target_total': totals[0],
			'target_percent': float(row[6].rstrip('%')),
			'background': float(row[7]), 'background_total': totals[1],
			'background_percent': float(row[8].rstrip('%')),
		})
	return results

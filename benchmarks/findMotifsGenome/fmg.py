"""The findMotifsGenome.pl benchmark's harness: build, environment, one case, manifests.

bench.py, time_one.py, make_refs.py, make_denovo80.py, make_inputs.py and
compare.py import this file from this directory.

Paths come from the environment:

- HOMER_FMG_REF: the reference directory (inputs/, preparsed/, cases/,
  v51/homer, the HOMER v5.1 build, denovo80/).
- HOMER_BENCH_GENOME: the hg38 FASTA.
- HOMER_FMG_CHECK_INPUTS: compare.py's small inputs (default:
  $HOMER_FMG_REF/checks/inputs).

findMotifsGenome.pl writes its command line into motifFindingParameters.txt,
so the genome and preparsed paths are part of every output a reference
holds: run with the paths the references were made with.

Every HOMER call runs with the same fixed randomness and environment, for
HOMER v5.1 and for every tree compared with it:

- PERL_HASH_SEED=0, PERL_PERTURB_KEYS=0: Perl's hash order, which decides the
  order of rand() draws in preparseGenome.pl.
- PERL5OPT=-MBenchSrand: srand(SRAND) at every perl start (rand() picks the
  temporary file names, the preparse positions and the background subsample).
- PATH: the step wrappers, then the build's bin/, then the system directories
  only (no other HOMER on PATH); LC_ALL=C; TMPDIR and cwd a fresh directory on
  /dev/shm; OMP_NUM_THREADS = the case's -p.
"""

import os
import re
import sys
import json
import time
import fcntl
import shutil
import hashlib
import tempfile
import subprocess

REF = os.environ['HOMER_FMG_REF']
HARNESS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'harness')
INPUTS = os.path.join(REF, 'inputs')
CHECK_INPUTS = os.environ.get('HOMER_FMG_CHECK_INPUTS') or os.path.join(REF, 'checks',
	'inputs')
PREPARSED = os.path.join(REF, 'preparsed')
CASES_REF = os.path.join(REF, 'cases')
V51 = os.path.join(REF, 'v51', 'homer')

GENOME = os.environ['HOMER_BENCH_GENOME']
SCRATCH = '/dev/shm'
SRAND = 1
BUILD_JOBS = 8

# case: peaks file, -N, the -p it runs at, and the trial count homer2 denovo
# uses at that -p (v5.1 rounds -T 10 to a multiple of -p), which the v5.1
# -p 1 reference is given through the homer2 wrapper.
CASES = {
	'objective': {'peaks': 'peaks10k.bed', 'N': 50000, 'p': 4, 'T': 12},
	# p1_2k: a larger -p 1 case, which bench.py does not run
	'p1': {'peaks': 'peaks500.bed', 'N': 2500, 'p': 1, 'T': 10},
	'p1_2k': {'peaks': 'peaks2k.bed', 'N': 10000, 'p': 1, 'T': 10},
	'p16': {'peaks': 'peaks2k.bed', 'N': 10000, 'p': 16, 'T': 16},
	'ctcf': {'peaks': 'ctcf10k.bed', 'N': 50000, 'p': 4, 'T': 12},
	# determinism tests only
	'tiny': {'peaks': 'peaks300.bed', 'N': 2000, 'p': 1, 'T': 10},
	'tiny4': {'peaks': 'peaks300.bed', 'N': 2000, 'p': 4, 'T': 12},
}

# Names in a tree that are not HOMER: never copied into a build.
NOT_HOMER = {'.git', 'runs', 'cache', 'probe', 'check', 'profile', 'unit-tests',
	'__pycache__', 'bench.py', 'validate.py', 'result.json', 'curves.json',
	'report.txt', 'meta.json', 'README.md', 'notes.md', 'check.json',
	'cold.json', 'config.txt', 'reference', 'log.txt', 'metrics.json'}

# Programs findMotifsGenome.pl starts by name whose wall time is measured.
WRAPPED = ('homer2', 'findKnownMotifs.pl', 'compareMotifs.pl')


def log(*items):
	print(*items, flush=True)


# --- the portable-build rule --------------------------------------------------

NONPORTABLE = re.compile(r'-march=native|-mtune=native|-march=(?!x86-64\b)\S+'
	r'|-mcpu=native|-m(?:avx\w*|sse4\w*|sse3|ssse3|fma\w*|bmi\w*|popcnt|f16c|lzcnt)\b'
	r'|-ffast-math|-Ofast|-funsafe-math\S*|-ffinite-math\S*|-fno-signed-zeros'
	r'|-fassociative-math|-freciprocal-math|-fno-trapping-math')
PRAGMA = re.compile(r'#\s*pragma\s+(?:GCC|clang)\s+(?:target|optimize)\b.*'
	r'|__attribute__\s*\(\(\s*(?:target|optimize)\s*\(')


def nonportable_build(src):
	"""Build or source lines that make a binary CPU-specific or relax
	floating-point semantics; target_clones (runtime dispatch) is allowed."""

	hits = []
	for root, dirs, files in os.walk(src):
		dirs[:] = [d for d in dirs if d not in NOT_HOMER]
		for name in files:
			path = os.path.join(root, name)
			is_build = (name.startswith(('Makefile', 'makefile', 'GNUmakefile',
				'CMakeLists')) or name.endswith(('.sh', '.mk', '.cmake'))
				or name == 'configureHomer.pl')
			is_source = name.endswith(('.cpp', '.h', '.hpp', '.cc', '.c', '.hh'))
			if not (is_build or is_source):
				continue
			with open(path, errors='replace') as handle:
				for number, line in enumerate(handle, 1):
					# a comment in a Makefile, script or configureHomer.pl sets nothing
					code = line.split('#', 1)[0] if is_build else line
					if (is_build and NONPORTABLE.search(code)) or (
							is_source and PRAGMA.search(line)):
						hits.append('{}:{}: {}'.format(os.path.relpath(path, src),
							number, line.strip()[:200]))
	return hits


# --- build ----------------------------------------------------------------------

def configure(home, logfile):
	"""`configureHomer.pl -make` in `home`: sets the scripts' `use lib` and
	$homeDir to `home` and compiles cpp/ (make clean, make) with
	MAKEFLAGS=-j BUILD_JOBS. Raises when a core binary is missing."""

	shutil.copy(os.path.join(HARNESS, 'config.txt'), os.path.join(home,
		'config.txt'))
	env = dict(os.environ, MAKEFLAGS='-j{}'.format(BUILD_JOBS), LC_ALL='C')
	for key in [k for k in env if k.startswith('PERL')]:
		del env[key]
	start = time.perf_counter()
	with open(logfile, 'w') as out:
		code = subprocess.run(['perl', os.path.join(home, 'configureHomer.pl'),
			'-make'], cwd=home, stdout=out, stderr=subprocess.STDOUT,
			env=env).returncode
	missing = [name for name in ('homer2', 'homerTools', 'mergePeaks')
		if not os.access(os.path.join(home, 'bin', name), os.X_OK)]
	if code != 0 or missing:
		with open(logfile, errors='replace') as handle:
			sys.stdout.write(handle.read()[-4000:])
		raise RuntimeError('build in {} failed (exit {}, missing {})'.format(home,
			code, missing))
	return time.perf_counter() - start


def copy_tree(src, dest):
	"""HOMER's files from a tree, without the runner's and the harness's own
	files."""

	os.makedirs(dest)
	for name in os.listdir(src):
		if name in NOT_HOMER:
			continue
		path = os.path.join(src, name)
		if os.path.isdir(path) and not os.path.islink(path):
			shutil.copytree(path, os.path.join(dest, name), symlinks=True,
				ignore=shutil.ignore_patterns('__pycache__', '*.o'))
		else:
			shutil.copy2(path, os.path.join(dest, name), follow_symlinks=False)


def build_cached(src, cache):
	"""Build `src` once per commit into `<cache>/homer`; later runs of the same
	commit reuse it. Returns the install directory."""

	os.makedirs(cache, exist_ok=True)
	home = os.path.join(cache, 'homer')
	stamp = os.path.join(home, '.bench-built')
	with open(os.path.join(cache, '.build.lock'), 'w') as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		if not os.path.exists(stamp):
			if os.path.exists(home):
				shutil.rmtree(home)
			copy_tree(src, home)
			seconds = configure(home, os.path.join(cache, 'build.log'))
			with open(stamp, 'w') as handle:
				handle.write('{:.1f}\n'.format(seconds))
			log('built {} in {:.1f} s'.format(home, seconds))
	return home


def build_inplace(src):
	"""Build a tree in place, once."""

	stamp = os.path.join(src, '.bench-built')
	with open(os.path.join(src, '.build.lock'), 'w') as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		if not os.path.exists(stamp):
			seconds = configure(src, os.path.join(src, '.build.log'))
			with open(stamp, 'w') as handle:
				handle.write('{:.1f}\n'.format(seconds))
	return src


# --- environment and wrappers ---------------------------------------------------

def homer_env(home, threads, scratch, wrappers=None):
	env = {key: value for key, value in os.environ.items()
		if not key.startswith(('PERL', 'OMP_'))}
	path = ([wrappers] if wrappers else []) + [os.path.join(home, 'bin'),
		'/usr/local/bin', '/usr/bin', '/bin']
	env.update({
		'PATH': ':'.join(path),
		'PERL_HASH_SEED': '0',
		'PERL_PERTURB_KEYS': '0',
		'PERL5LIB': os.path.join(HARNESS, 'perl'),
		'PERL5OPT': '-MBenchSrand',
		'OMP_NUM_THREADS': str(threads),
		'LC_ALL': 'C',
		'TMPDIR': scratch,
	})
	env.pop('LANG', None)
	env.pop('LANGUAGE', None)
	return env


WRAPPER = r'''#!/bin/bash
s=$EPOCHREALTIME
{extra}{env}"{real}" "$@"
rc=$?
printf '%s\t%s\t%s\t%s\t%s\n' "{name}" "$s" "$EPOCHREALTIME" "$rc" "$*" >> "{steplog}"
exit $rc
'''


COMPILED = ('homer2', 'homerTools', 'mergePeaks')


def make_wrappers(directory, home, steplog, denovo_trials=None, profile=None):
	"""Wrappers on PATH around the programs whose wall time is measured.
	With `denovo_trials`, `homer2 denovo` also gets `-T <n>` (the v5.1 -p 1
	reference's trial count). With `profile` = (sampler.so, directory), the
	compiled programs run with the SIGPROF sampler preloaded, one file each."""

	os.makedirs(directory, exist_ok=True)
	for name in WRAPPED + (COMPILED[1:] if profile else ()):
		extra = env = ''
		if name == 'homer2' and denovo_trials is not None:
			extra = 'if [ "$1" = denovo ]; then set -- "$@" -T {}; fi\n'.format(
				int(denovo_trials))
		if profile and name in COMPILED:
			env = 'SAMPLER_OUT="{}/{}.$$.txt" LD_PRELOAD="{}" '.format(profile[1], name,
				profile[0])
		path = os.path.join(directory, name)
		with open(path, 'w') as handle:
			handle.write(WRAPPER.format(extra=extra, env=env, real=os.path.join(home,
				'bin', name), name=name, steplog=steplog))
		os.chmod(path, 0o755)
	return directory


def steps(steplog, t0):
	"""Per-step wall seconds from the wrappers' log."""

	out = {'seconds_bg': 0.0, 'seconds_norm': 0.0, 'seconds_known': 0.0,
		'seconds_denovo8': 0.0, 'seconds_denovo10': 0.0, 'seconds_denovo12': 0.0,
		'seconds_compare': 0.0}
	if not os.path.exists(steplog):
		return out, []
	rows = []
	with open(steplog) as handle:
		for line in handle:
			parts = line.rstrip('\n').split('\t', 4)
			if len(parts) == 5:
				rows.append((parts[0], float(parts[1]), float(parts[2]),
					parts[3], parts[4]))
	starts = [row[1] for row in rows if row[0] in WRAPPED]
	if starts:
		out['seconds_bg'] = max(0.0, min(starts) - t0)
	for name, start, end, _, argv in rows:
		seconds = end - start
		words = argv.split()
		if name == 'findKnownMotifs.pl':
			out['seconds_known'] += seconds
		elif name == 'compareMotifs.pl':
			out['seconds_compare'] += seconds
		elif name == 'homer2' and words[:1] == ['norm']:
			out['seconds_norm'] += seconds
		elif name == 'homer2' and words[:1] == ['denovo'] and '-len' in words:
			key = 'seconds_denovo' + words[words.index('-len') + 1]
			if key in out:
				out[key] += seconds
	return out, rows


# --- manifests -----------------------------------------------------------------

P_ECHO = re.compile(rb'(?<=\s)-p \d+(?=\s|$)')


def normalize(relpath, data):
	"""The -p value echoed in motifFindingParameters.txt is not an output."""

	if os.path.basename(relpath) == 'motifFindingParameters.txt':
		return P_ECHO.sub(b'-p <P>', data)
	return data


def digest(data):
	return hashlib.sha256(data).hexdigest()


def file_digest(path, relpath):
	if os.path.basename(relpath) == 'motifFindingParameters.txt':
		with open(path, 'rb') as handle:
			return digest(normalize(relpath, handle.read()))
	sha = hashlib.sha256()
	with open(path, 'rb') as handle:
		for block in iter(lambda: handle.read(1 << 22), b''):
			sha.update(block)
	return sha.hexdigest()


def manifest(directory, stdout_path, code):
	"""sha256 of every file under `directory` (recursively, normalized), of
	stdout, and the exit status."""

	out = {}
	if directory is not None and os.path.isdir(directory):
		for root, _, files in os.walk(directory):
			for name in files:
				path = os.path.join(root, name)
				relpath = os.path.relpath(path, directory)
				if os.path.islink(path):
					out[relpath] = 'symlink:' + os.readlink(path)
				else:
					out[relpath] = file_digest(path, relpath)
	with open(stdout_path, 'rb') as handle:
		out['<stdout>'] = digest(handle.read())
	out['<exit>'] = str(code)
	return out


def mismatches(got, want):
	"""Keys whose hash differs, or that one side lacks."""

	return sorted(key for key in set(got) | set(want) if got.get(key) != want.get(key))


def load_reference(case):
	path = os.path.join(CASES_REF, case, 'manifest.json')
	if not os.path.exists(path):
		raise FileNotFoundError('no HOMER v5.1 reference for case {} at {}: run '
			'make_refs.py'.format(case, path))
	with open(path) as handle:
		return json.load(handle)


def first_difference(path_a, path_b, limit=200):
	"""The first differing line of two text files, for the log."""

	try:
		with open(path_a, 'rb') as a, open(path_b, 'rb') as b:
			for number, (x, y) in enumerate(zip(a, b), 1):
				if x != y:
					return 'line {}: here {!r} | v5.1 {!r}'.format(number,
						x[:limit], y[:limit])
		return 'one file is a prefix of the other'
	except OSError as error:
		return str(error)


# --- one case ------------------------------------------------------------------

def run_command(argv, cwd, env, stdout_path, stderr_path, rss_path=None):
	"""Run argv; returns (exit status, wall seconds, user+sys CPU seconds of
	the whole tree, peak RSS in MB of its largest process or None)."""

	if rss_path is not None:
		argv = ['/usr/bin/time', '-f', '%M', '-o', rss_path] + argv
	with open(stdout_path, 'wb') as out, open(stderr_path, 'wb') as err:
		start = time.perf_counter()
		process = subprocess.Popen(argv, cwd=cwd, env=env, stdin=subprocess.DEVNULL,
			stdout=out, stderr=err, start_new_session=True)
		_, status, usage = os.wait4(process.pid, 0)
		seconds = time.perf_counter() - start
	process.returncode = os.waitstatus_to_exitcode(status)
	rss = None
	if rss_path is not None and os.path.exists(rss_path):
		with open(rss_path) as handle:
			words = handle.read().split()
		if words and words[-1].isdigit():
			rss = int(words[-1]) / 1024
	return process.returncode, seconds, usage.ru_utime + usage.ru_stime, rss


def run_findmotifs(home, case, p=None, trials=None, rss=False, keep=None,
		label=None, profile=None):
	"""One findMotifsGenome.pl call on `case`'s peaks in a fresh directory on
	/dev/shm, laid out the same for every run (peaks.bed, out/, the reference
	preparsed directory), so the echoed command line matches. `p` overrides
	the case's -p; `trials` gives homer2 denovo -T (the v5.1 reference).
	With `keep`, the run directory is moved there; otherwise removed.

	Returns a dict: code, seconds, cpu_seconds, peak_rss_mb, manifest, steps,
	and dir (kept directory or None)."""

	spec = CASES[case]
	threads = spec['p'] if p is None else p
	scratch = tempfile.mkdtemp(dir=SCRATCH, prefix='fmg-{}-'.format(label or case))
	try:
		shutil.copy(os.path.join(INPUTS, spec['peaks']), os.path.join(scratch,
			'peaks.bed'))
		steplog = os.path.join(scratch, '.steps.tsv')
		wrappers = make_wrappers(os.path.join(scratch, '.wrap'), home, steplog,
			trials, profile)
		env = homer_env(home, threads, scratch, wrappers)
		argv = ['findMotifsGenome.pl', 'peaks.bed', GENOME, 'out', '-size', '200',
			'-N', str(spec['N']), '-p', str(threads), '-preparsedDir', PREPARSED]
		t0 = time.time()
		code, seconds, cpu, peak = run_command(argv, scratch, env,
			os.path.join(scratch, '.stdout'), os.path.join(scratch, '.stderr'),
			os.path.join(scratch, '.rss') if rss else None)
		result = {'code': code, 'seconds': seconds, 'cpu_seconds': cpu,
			'peak_rss_mb': peak, 'argv': argv, 'threads': threads,
			'manifest': manifest(os.path.join(scratch, 'out'), os.path.join(
				scratch, '.stdout'), code)}
		result['steps'], result['step_rows'] = steps(steplog, t0)
		result['dir'] = None
		if keep is not None:
			shutil.move(scratch, keep)
			result['dir'] = keep
		return result
	finally:
		if os.path.exists(scratch):
			shutil.rmtree(scratch, ignore_errors=True)


def run_preparse(home, keep=None, label='preparse'):
	"""A cold preparseGenome.pl for the genome at -size 200 into an empty
	directory, as findMotifsGenome.pl calls it when the files are missing."""

	scratch = tempfile.mkdtemp(dir=SCRATCH, prefix='fmg-{}-'.format(label))
	try:
		target = os.path.join(scratch, 'preparsed')
		env = homer_env(home, 4, scratch)
		argv = ['preparseGenome.pl', GENOME, '-size', '200', '-preparsedDir', target]
		code, seconds, cpu, _ = run_command(argv, scratch, env,
			os.path.join(scratch, '.stdout'), os.path.join(scratch, '.stderr'))
		result = {'code': code, 'seconds': seconds, 'cpu_seconds': cpu,
			'argv': argv, 'manifest': manifest(target, os.path.join(scratch,
				'.stdout'), code), 'dir': None}
		if keep is not None:
			shutil.move(scratch, keep)
			result['dir'] = keep
		return result
	finally:
		if os.path.exists(scratch):
			shutil.rmtree(scratch, ignore_errors=True)


def warm_files():
	"""Read the genome and the preparsed files once, so the timed call finds
	them in the page cache (not timed)."""

	start = time.perf_counter()
	paths = [GENOME] + [os.path.join(PREPARSED, name) for name in
		sorted(os.listdir(PREPARSED))]
	for path in paths:
		with open(path, 'rb') as handle:
			while handle.read(1 << 24):
				pass
	return time.perf_counter() - start


def report_mismatch(case, got, want, here_dir=None, reference_dir=None):
	"""Log each differing key; returns the count."""

	keys = mismatches(got, want)
	for key in keys[:20]:
		detail = ''
		if here_dir and reference_dir and not key.startswith('<') and key in got \
				and key in want:
			detail = ' ' + first_difference(os.path.join(here_dir, key),
				os.path.join(reference_dir, key))
		log('MISMATCH {}: {} ({}){}'.format(case, key, 'missing here' if key not in
			got else 'missing in v5.1' if key not in want else 'differs', detail))
	if len(keys) > 20:
		log('MISMATCH {}: ... {} more'.format(case, len(keys) - 20))
	return len(keys)


if __name__ == '__main__':
	# python3 fmg.py build <tree>  (in place, once; as bench.py builds its
	# copy: configureHomer.pl -make, MAKEFLAGS=-j8)
	if len(sys.argv) == 3 and sys.argv[1] == 'build':
		build_inplace(os.path.realpath(sys.argv[2]))
		print('built', os.path.realpath(sys.argv[2]))
	else:
		sys.exit('usage: fmg.py build <tree>')

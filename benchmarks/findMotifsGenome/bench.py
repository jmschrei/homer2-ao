"""One run of the findMotifsGenome.pl benchmark: build, warm, time, check.

Objective: wall seconds of one call with HOMER's own interface

	findMotifsGenome.pl peaks.bed <hg38.fa> out -size 200 -N 50000 -p 4 \
		-preparsedDir <reference>/preparsed

on the 10,000 highest-signal K562 ATAC-seq peaks of ENCODE ENCFF333TAT (200 bp
at the summit), process start to exit with every child, the preparsed hg38
background files already present (read-only) and the genome and those files
read once beforehand (not timed). Then, each timed once and checked: 500 peaks
(ranks 10,001-10,500) with -N 2500 at -p 1, 2,000 peaks (ranks
10,001-12,000) with -N 10000 at -p 16, and a cold preparseGenome.pl of hg38
at -size 200 into an empty directory.

Check n_file_mismatch_vs_v51: every file each call leaves (the output
directory recursively; the preparsed directory), stdout and the exit status,
against HOMER v5.1 run once per case at -p 1 with homer2 denovo given the
case's trial count (make_refs.py). The -p value echoed in
motifFindingParameters.txt is normalized. The fixed environment (Perl hash
seed 0, srand(1) at every perl start, PATH, LC_ALL=C, a fresh /dev/shm
directory per call) is fmg.py's, the same for v5.1 and every tree.

The seed varies nothing: every input and the Perl random stream are fixed, so
runs of one tree are repeats and their spread is machine noise.

	HOMER_FMG_REF=<reference> HOMER_BENCH_GENOME=<hg38.fa> python3 bench.py \
		--tree <HOMER tree> --seed 0 --out <dir> --cache <dir>
"""

import os
import sys
import json
import shutil
import tarfile
import argparse
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fmg

parser = argparse.ArgumentParser()
parser.add_argument('--tree', required=True)
parser.add_argument('--seed', type=int, required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--cache', required=True)
args = parser.parse_args()

SRC = os.path.abspath(args.tree)
OUT = os.path.abspath(args.out)
os.makedirs(OUT, exist_ok=True)
log = fmg.log

# Every reference first: a missing one is an error, never written here.
references = {case: fmg.load_reference(case) for case in ('objective', 'p1', 'p16',
	'preparse')}

hits = fmg.nonportable_build(SRC)
for hit in hits:
	log('nonportable:', hit)
checks = {'n_file_mismatch_vs_v51': 0, 'nonportable_build': len(hits)}
measurements = {}
cases = {}

home = fmg.build_cached(SRC, os.path.abspath(args.cache))
log('warm read of the genome and preparsed files: {:.1f} s'.format(fmg.warm_files()))


def check(case, result):
	kept = result.get('dir')
	sub, ref = ('preparsed', fmg.PREPARSED) if case == 'preparse' else ('out',
		os.path.join(fmg.CASES_REF, case, 'out'))
	n = fmg.report_mismatch(case, result['manifest'], references[case],
		os.path.join(kept, sub) if kept else None, ref)
	checks['n_file_mismatch_vs_v51'] += n
	cases[case] = {'seconds': result['seconds'], 'cpu_seconds': result['cpu_seconds'],
		'exit': result['code'], 'files': len(result['manifest']) - 2, 'mismatches': n}
	log('{:10s} {:9.2f} s  exit {}  {} files  {} mismatched'.format(case,
		result['seconds'], result['code'], len(result['manifest']) - 2, n))
	return n


def keep_dir(case):
	return tempfile.mkdtemp(dir=fmg.SCRATCH, prefix='fmg-keep-{}-'.format(case))


def save_logs(case, kept, outputs=False):
	"""stderr, the step log and, for the objective, its output directory."""

	for name in ('.stderr', '.stdout', '.steps.tsv'):
		path = os.path.join(kept, name)
		if os.path.exists(path):
			shutil.copy(path, os.path.join(OUT, '{}.{}'.format(case, name[1:])))
	if outputs and os.path.isdir(os.path.join(kept, 'out')):
		with tarfile.open(os.path.join(OUT, case + '_out.tar.gz'), 'w:gz') as tar:
			tar.add(os.path.join(kept, 'out'), arcname='out')


# The objective, under /usr/bin/time for peak RSS, the step wrappers on PATH.
kept = keep_dir('objective')
os.rmdir(kept)
result = fmg.run_findmotifs(home, 'objective', rss=True, keep=kept)
try:
	check('objective', result)
	save_logs('objective', kept, outputs=True)
finally:
	shutil.rmtree(kept, ignore_errors=True)
if result['code'] != 0:
	sys.exit('the objective call exited {}; see objective.stderr'.format(result['code']))
seconds = result['seconds']
measurements['peak_rss_mb'] = result['peak_rss_mb']
measurements['cpu_seconds'] = result['cpu_seconds']
measurements.update(result['steps'])

for case in ('p1', 'p16'):
	kept = keep_dir(case)
	os.rmdir(kept)
	result = fmg.run_findmotifs(home, case, keep=kept)
	try:
		check(case, result)
		save_logs(case, kept)
	finally:
		shutil.rmtree(kept, ignore_errors=True)
	measurements['seconds_' + case] = result['seconds']

kept = keep_dir('preparse')
os.rmdir(kept)
result = fmg.run_preparse(home, keep=kept)
try:
	check('preparse', result)
	save_logs('preparse', kept)
finally:
	shutil.rmtree(kept, ignore_errors=True)
measurements['seconds_preparse'] = result['seconds']

log('objective {:.2f} s; checks {}'.format(seconds, json.dumps(checks)))
log('measurements {}'.format(json.dumps({k: round(v, 2) if isinstance(v, float) else v
	for k, v in measurements.items()})))

with open(os.path.join(OUT, 'cases.json'), 'w') as handle:
	json.dump(cases, handle, indent=1)
with open(os.path.join(OUT, 'metrics.json'), 'w') as handle:
	json.dump({'objective': seconds, 'checks': checks, 'measurements': measurements,
		'note': 'objective_out.tar.gz holds the objective call\'s output directory'},
		handle, indent=1)

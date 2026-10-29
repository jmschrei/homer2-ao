"""Time one benchmark call on a HOMER tree, for paired timing of two trees.

	HOMER_FMG_REF=<reference> HOMER_BENCH_GENOME=<hg38.fa> python3 time_one.py --src <tree>

Builds the tree in place once (fmg.py build_inplace), reads the genome and
the preparsed files once, then times one findMotifsGenome.pl call laid out
exactly as bench.py's objective (10,000 ATAC peaks, -N 50000, -p 4). With
OMP_NUM_THREADS=16 in the environment it times the -p 16 case instead (2,000
peaks, -N 10000). Every output file, stdout and the exit status must match
HOMER v5.1's -p 1 reference, or it exits 1. The seconds are the last number
printed.
"""

import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fmg

parser = argparse.ArgumentParser()
parser.add_argument('--src', required=True)
args = parser.parse_args()

case = 'p16' if os.environ.get('OMP_NUM_THREADS') == '16' else 'objective'
reference = fmg.load_reference(case)
home = fmg.build_inplace(os.path.realpath(args.src))
fmg.warm_files()
result = fmg.run_findmotifs(home, case, label='time-' + case)
n = fmg.report_mismatch(case, result['manifest'], reference)
if n:
	print('{} files differ from HOMER v5.1 -p 1'.format(n))
	sys.exit(1)
print('{} -p {}: {:.3f}'.format(case, result['threads'], result['seconds']))

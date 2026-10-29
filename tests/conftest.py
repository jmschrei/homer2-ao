"""Fixtures shared by HOMER's tests. See homer_helpers.py for the
environment variables the suite reads.

Behaviour that contradicts HOMER's usage text, or is otherwise wrong, in
v5.1 is left untested, so that fixing it will not break the suite. No test
asserts any of these cases; where a test runs an affected option to check
something else, it does not check the affected output:

homer2
- find -mscore prints an unset position and strand for a sequence with no
  scorable window (all N).
- A tab in a FASTA name is read as the start of the sequence; 999 tabs
  overflow a fixed column array.
- known -mscore segfaults.
- A FASTA record of 1,000,000 bases or more overflows a fixed buffer in the
  conversion -i and -b use.
- bg with all-N target sequences segfaults.
- denovo -opt frees uninitialized pointers, so its outcome depends on what
  the heap held.
- denovo accepts a negative -len or -p.
- At -p above 1, the shared p-value cache hands threads entries that are not
  filled yet, so denovo's oligo-screen statistics vary from run to run (the
  known tests at -p 4 compare counts only, which the cache cannot change).
- An -i or -b FASTA file that is missing or empty ends the run with exit
  status 0 and leaves the temporary .seq and .group files behind.
- known -siteReduce reports every motif with p-value 1 (log p 0), and
  Motif::scoreOverlap shifts the first motif's site instead of the second's.
- known -opt writes the optimized threshold with site statistics (Tpos, Tstd,
  Bpos, Bstd, StrandBias, Multiplicity) computed at the input threshold.
- norm and known: -nmax N runs N + 1 rounds, -nmax 0 still normalizes once
  and writes no -nout file, and the k-mer count leaves out the last window of
  every sequence (the norm tests accept either count, and either round count).
- bg: -maxIterations does nothing; -maxFractionN is not applied to -i targets
  or -b background; -p targets with -b or -bg background read an
  uninitialized array; a -N that is not a multiple of the number of GC bins
  selects up to one extra sequence per bin.
- denovo: -nozoops does nothing with the default -fullMask; -omax 0 writes no
  -oout file; the -e default depends on -len although the usage gives 0.005.

homerTools
- trim compares the adapter without the read's last base; -q averages a
  window shifted by one base and records the kept length minus one;
  -stats <file> writes nothing; -lenSuffix drops the suffix for single-end
  files; -5 # with -3 # gives sequence and quality strings of different
  lengths; the .lengths table counts a phantom read (and -min 0 writes an
  empty record), its fractions can add up to more than 100%, and with -pe both
  mates' tables go to the first file.
- barcodes and truseq ignore -freq <file>; truseq documents -qual but reads it
  and its value as input files.
- cluster rejects -dist euclidian, the usage's spelling, and -sub never
  finishes when it has to sample.
- freq segfaults on a FASTA or FASTQ file with one record.
- extract mishandles peaks that start before base 1 or end past the end of
  the chromosome.

makeTagDirectory
- -genome drops a trailing "r" from a custom genome path (findMotifsGenome.pl
  does the same), -normFixedOligo reports the old plus the new total in
  tagInfo.txt, and -checkGC's genome GC and oligo tables leave out the window
  that ends at each chromosome's last base (the tests accept either count).

findMotifsGenome.pl, findMotifs.pl and compareMotifs.pl
- findMotifsGenome.pl -tss and -cgtss never write the group file, and the run
  reports no enrichment with exit status 0.
- findMotifsGenome.pl with an output directory it cannot use regenerates the
  files in -preparsedDir and exits 0.
- findMotifsGenome.pl -pc and -depth report that they took effect but change
  nothing.
- knownResults.txt and knownResults.html give a background total of 1 when
  the first motif has no background hits.
- findMotifsGenome.pl -N counts targets and background together, while the
  usage calls it the number of background sequences.
- findMotifsGenome.pl ignores -mask with -useNewBg, leaves a coBoundBy1 file in
  the output directory with -bg or -ref, and exits 0 on errors (a missing
  genome, -size 50,-100, an unknown option).
- findMotifs.pl -olen passes -o to homer2 denovo, -basic and -nofacts
  overwrite each other, -flip uses an undefined file name, and -find leaves a
  .seq2.tmp file in the output directory.
- compareMotifs.pl documents -siteReduce and -reducedMotifs but rejects both,
  -rna replaces any -known file whose path contains "all.motifs", and its
  information content uses ln(p)/2 instead of log2(p), so -info removes
  almost nothing.

annotatePeaks.pl
- -size a,b takes b - a + 2 bases instead of b - a + 1, and with an odd span
  shifts every motif position by one; -hist with
  -size given leaves the nucleotide columns empty after the first row; the
  -matrix statistics file swaps its "Peaks with Motif1" and "Motif2" columns;
  -fm also drops sites that only touch a filter site.
"""

import os
import re

import pytest

import homer_helpers as hh


@pytest.fixture(scope='session', autouse=True)
def homer_install():
	"""Stop the session with one message when HOMER_HOME is not a configured,
	built install, instead of failing every test."""

	missing = [name for name in ('homer2', 'homerTools', 'mergePeaks')
		if not os.access(os.path.join(hh.HOMER_BIN, name), os.X_OK)]
	script = os.path.join(hh.HOMER_BIN, 'findMotifsGenome.pl')
	configured = False
	if os.path.exists(script):
		with open(script) as handle:
			head = handle.read(2000)
		found = re.search(r'my \$homeDir = "([^"]*)"', head)
		configured = bool(found) and os.path.realpath(found.group(1)) == \
			os.path.realpath(hh.HOMER_HOME)
		# the pip package's scripts find their directory at run time
		configured = configured or 'my $homeDir = Cwd::abs_path(' in head
	if missing or not configured:
		pytest.exit('HOMER_HOME={} is not a built HOMER install: run '
			'"perl configureHomer.pl -make" there, or point HOMER_HOME at one '
			'(missing programs: {}; scripts configured for this directory: {})'.format(
			hh.HOMER_HOME, ', '.join(missing) or 'none', 'yes' if configured else 'no'),
			returncode=4)
	return hh.HOMER_HOME


@pytest.fixture(scope='session')
def genome():
	return hh.read_genome()


@pytest.fixture(scope='session')
def peaks():
	return hh.read_peaks()


@pytest.fixture(scope='session')
def preparsed(tmp_path_factory):
	"""A directory of the genome's preparsed background for -size 200, made
	once per session; tests only read it."""

	directory = tmp_path_factory.mktemp('preparsed')
	hh.run(['preparseGenome.pl', hh.GENOME, '-size', 200, '-preparsedDir',
		directory], cwd=directory)
	return directory

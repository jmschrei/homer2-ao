# Validation

How homer2-ao's output was compared with HOMER v5.1's. Each comparison ran both builds with the
same inputs, Perl's hash seed and `rand()` fixed, and counted a case as identical only when
every file the run left, stdout and the exit status matched byte for byte. The differences that
remain are listed in `docs/DIFFERENCES.md`; the v5.1 defects these runs found are in
`docs/V5.1-DEFECTS.md`.

## Unit tests

`tests/` holds 741 pytest tests of the motif programs: `findMotifsGenome.pl`, `findMotifs.pl`,
`compareMotifs.pl`, every `homer2` command, `homerTools`, the motif options of
`annotatePeaks.pl`, `scanMotifGenomeWide.pl` and `makeTagDirectory -genome`. Expected values are
computed from the inputs by the tests rather than recorded from HOMER's output. The suite passes
on v5.1 (the `main` branch before the merge) and on homer2-ao, built with GCC and with Clang.

```
perl configureHomer.pl -make
touch config.txt
python3 -m pytest tests
```

On 2026-10-28, on fresh clones with Python 3.13.5 and pytest 8.3.4, eight program calls at a
time: `main` before the merge (v5.1) built with g++ 13.3, 741 passed in 98 s; `optimize` built
with g++ 13.3, 741 passed in 66 s; `optimize` built with Clang 23.1.4 and the system's
libstdc++, 741 passed in 71 s. The same day the pip package, installed with `uv pip install`
into fresh environments from its wheel and from its sdist (built there with g++ 13.3), passed
741 of 741 each (`HOMER_HOME=$(homer2-ao --home)`, 64 s). A small `findMotifsGenome.pl` run
(the tests' 200 peaks, `-len 8,10 -S 5 -p 4`) through the package's command and through the
in-place build gave the same 102 output files, stdout and stderr.

The manylinux wheel is compiled in the `manylinux_2_28` image (AlmaLinux 8, GCC 14.2.1), so it
was checked again on 2026-10-28:

- Installed on the machine above, it passed 741 of 741 tests.
- Its `homer2 find` on the 269,800-peak benchmark (`-p 4`) wrote output identical to v5.1's at
  `-p 1`: 13,233,229 lines with the same SHA-256.
- The three `findMotifsGenome.pl` benchmark cases (`-p 4`, `-p 1` and `-p 16`) matched v5.1 in
  every one of their 1,015, 588 and 798 outputs (files, stdout and exit status).
- Installed in AlmaLinux 8 itself (glibc 2.28, Perl 5.26.3), it passed 740 of 741 tests. The
  other, `test_options_passed_to_homer2_denovo[dn_e]`, checks that `-e 1` changes the de novo
  motifs on the test data. With Perl 5.26.3 the run with `-e 1` finds the same motifs as the
  run without it, and v5.1 built in the same environment fails that test the same way.
- Through the package's commands, `configureHomer.pl -install sacCer3` installed the yeast
  genome, and `annotatePeaks.pl` and `findMotifsGenome.pl` then ran with it by name. After the
  package was reinstalled, the genome was still listed in `config.txt` and `annotatePeaks.pl`
  gave the same annotations; uninstalling removed the programs and kept `config.txt` and the
  genome.

## Every commit

Each prefix of the `optimize` branch (v5.1 plus the commits up to that one) was built and run
through the `findMotifsGenome.pl` benchmark's three cases (`-p 4`, `-p 1` and `-p 16`) and a
cold `preparseGenome.pl`, comparing every file against v5.1 at `-p 1`, and through a check of 34
inputs: 32 cases of HOMER's other entry points, `homer2 denovo` on 80 small subsamples, and 12
runs that exercise the de novo rescreen. Every prefix gave identical output on all of them.

## homer2 on a collection of public test inputs (2026-10-12)

The homer2 find changes, applied alone to v5.1, against v5.1 on 7,163 cases covering all seven
`homer2` commands: 412 real peak sets from about 20 assays and 7 assemblies, 457 malformed or
unusual FASTA files, and 165 motif files including malformed ones and other formats. Of 10,333
runs, 10,154 were identical. The other 179:

- 105 were `homer2 denovo -p 4`, where v5.1's multithreaded oligo screen varies from run to run;
  every variant homer2-ao gave is one v5.1 also gives when run repeatedly.
- 52 crash identically in both, but v5.1 leaves its temporary `.seq` and `.group` files behind.
- 15 read gzip data as FASTA, which overflows a fixed column array in both builds; what follows
  is undefined and varies between runs of v5.1 too.
- 7 are one 4.6 Mb single-record FASTA file, which v5.1 aborts on (fixed buffer) and homer2-ao
  reads.

The 285 runs where v5.1 ends with one of its own error messages and the 894 where it exits 0 with
empty output matched byte for byte.

## findMotifsGenome.pl and the programs it runs (2026-10-27)

A development version, which differs from the released code only in comments and in the two
error-path fixes described below, against v5.1 on 1,988 cases: every tiny and small
interval file of 10 assemblies through `findMotifsGenome.pl` (724), every option value on one
hg38 file (198), `-p` 2 to 16, every malformed peak file (459), invalid arguments,
`compareMotifs.pl` with `compareMotifMatrices`, `homer2 denovo`, `known` and `norm` with their
options, threads, malformed and invalid input, and a cold preparse per genome and size.

- 1,962 cases gave identical files, stdout and exit status; 546 of these also had identical
  stderr. The stderr differences of the rest are of the kinds in `docs/DIFFERENCES.md`, or v5.1's
  own counters varying at `-p` above 1.
- 26 cases differ in more than stderr:
  - 9 are undecided: v5.1 was stopped at the time limit (600 s for malformed input, 7,200 s
    otherwise) while homer2-ao finished, and every file v5.1 had finished equals homer2-ao's.
  - 4 give invalid values to the new `-Tcpus` option.
  - 8 are v5.1's own run-to-run variation: 6 with temporary files named from the clock,
    `homer2 denovo -rand`, and one `-seqlogo` case.
  - 2 are a missing FASTA file, where v5.1 leaves temporary files behind and homer2-ao does not.
  - 3 were changes in homer2-ao: a bare `-p` at the end of the command line made `homer2 norm`
    segfault, `homer2 known -mscore` exits 0, and `homer2 denovo -len -1` crashes with another
    signal. Inside one of the undecided cases, `-nlen 15` at `-p 1` also aborted in
    `homer2 norm`. The released code fixes the two `norm` cases (its buffer is sized in `size_t`, and
    `-p` is passed only when it is a positive integer), and both now give v5.1's output; the
    other two are described in `docs/DIFFERENCES.md`.

`findMotifsGenome.pl` on the 10,000 strongest K562 CTCF ChIP-seq peaks (IDR conservative) gave
the same 1,095 files as v5.1 at `-p 1` (same development version).

## A larger comparison (in progress)

A larger comparison of the released code with v5.1 covers the commands the comparison above did
not run (`findMotifs.pl`, `homerTools` called directly, `makeTagDirectory`, `annotatePeaks.pl`,
`homer2 find`, `mask` and `bg`, and the batch scripts), compares the released code with the
development version tested above, and runs AddressSanitizer, UndefinedBehaviorSanitizer and
ThreadSanitizer builds of both under Clang: 33,401 cases. Its results will be added here.

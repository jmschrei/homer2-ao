# How homer2-ao differs from HOMER v5.1

homer2-ao is meant to give the same output files as HOMER v5.1, faster. On the benchmarks and
the validation runs in `docs/VALIDATION.md`, every output file, stdout and the exit status match
v5.1's, except for the cases this page lists and the cases `docs/VALIDATION.md` reports as
undecided. This page lists every difference known; each commit message on the `optimize` branch
gives the same information for its own change.

## New programs and options

- `homer2 denovo -Tcpus <#>` rounds the trial count `-T` for `<#>` CPUs instead of `-p`. v5.1
  rounds `-T` up to a multiple of `-p` (`-T 10` becomes 12 at `-p 4`), so the motifs depend on
  `-p` through that rounding. Without `-Tcpus` the rounding is v5.1's. `findMotifsGenome.pl`
  passes `-Tcpus <-p>` to every `homer2 denovo` it runs.
- `compareMotifMatrices` (new, `cpp/compareMotifMatrices.cpp`) scores motif pairs as
  `compareMotifs.pl` does, with the same arithmetic. `compareMotifs.pl` calls it and falls back
  to its own Perl code when the program is missing.
- `homerTools bgbins`, `randbg` and `filterids` (new) do the background-selection steps of
  `findMotifsGenome.pl` with a preparsed background, which v5.1 runs as `mergePeaks -cobound`,
  `filterListBy.pl`, `makeBinaryFile.pl` and `randRemoveBackground.pl`. The script keeps the old
  commands when overlap checking is off and for paths Perl's two-argument `open` would read
  differently. It falls back to `mergePeaks` when `bgbins` declines its input, and to
  `randRemoveBackground.pl` when `randbg` does.

## Output

- **Line order at `-p` above 1.** `homer2 find`, `find -mscore` and `mask` write at any `-p` the
  output v5.1 writes at `-p 1`. v5.1's own line order at `-p` above 1 depends on thread timing.
- **Enrichment score cache.** v5.1 publishes a new cache row before filling it, so at `-p` above
  1 a thread of `homer2 known` or `homer2 denovo` can read an unfilled value as a log p-value.
  homer2-ao fills each row before publishing it.
- **Long FASTA records.** `homer2 known`, `find`, `mask` and `norm` read records of 1,000,000
  bytes or more in full; v5.1 overflows a fixed buffer on them. `homer2 denovo -i` still reads
  FASTA the v5.1 way.
- **Missing FASTA input.** `homer2 known`, `find`, `mask` and `norm` no longer leave temporary
  `.seq` and `.group` files behind when an `-i` or `-b` file is missing.
- **`homer2 known -mscore`**, which is not an option of `known`: v5.1 segfaults; homer2-ao exits
  0 and reports every motif in 100% of targets and background. Neither result is meaningful.
- **`homer2 denovo -len -1`** crashes in both: v5.1 aborts (exit status 134), homer2-ao segfaults
  (139).
- **Tied sort keys in `homer2 denovo`.** homer2-ao sorts oligos with `std::stable_sort`, which
  gave glibc `qsort`'s order in all 2,748 sorts of a 10,000-peak run. Other sorts, such as the
  final motif list, still use the C library's `qsort`, whose order for ties differs on macOS: in
  one test two candidate motifs tie on every statistic, and v5.1 and homer2-ao both report the
  same one on macOS and the other on Linux. Output on macOS was not compared with v5.1 beyond the
  test suite.

## stderr

The output files are the same in each of these cases; only stderr differs.

- `mergePeaks`' co-bound progress lines are not printed when `homerTools bgbins` does that step.
  When `bgbins` declines its input (no targets, NUL bytes, over-long names), `mergePeaks` runs
  later and its background messages come in a different order.
- `filterListBy.pl`'s `readline() on closed filehandle` warning is gone: `homerTools filterids`
  reads a missing file as empty.
- With `-p 4` or more and two or more de novo lengths, `findMotifsGenome.pl` runs the longest
  length beside the rest, so step messages and the children's progress come in another order.
- The `homer2 denovo` screen's progress bar can be one `=` shorter; v5.1 visits one index past
  the end.
- `homer2 denovo`'s usage text has a `-Tcpus` line, and its `-T` adjustment message names the
  CPU count used for rounding.
- Perl warnings for a malformed `-cpu`, `-p`, `-stat` or `-fdr` value come earlier, because the
  scripts read these values earlier.
- `-tss` and `-cgtss`, which fail in both, print extra `rm` lines.

## Threads and memory

- `findMotifs.pl` and `findMotifsGenome.pl` pass `-p` to `homer2 norm`, which now uses that many
  threads. `findMotifsGenome.pl` passes it only when the value is a positive integer.
- With `-p 4` or more, de novo finding on, no `-fdr` or `-oligos`, and two or more distinct
  lengths, `findMotifsGenome.pl` forks `homer2 denovo` for the longest length on
  `-p - ceil(-p/4)` threads and runs known motifs and the other lengths on the remaining
  `ceil(-p/4)`, so no more than `-p` threads run at once.
- Peak memory falls for `homer2 find` (5,579 to 544 MB on 269,800 peaks at `-p 4`) and is about
  the same for `findMotifsGenome.pl` (1,977 and 1,882 MB on the benchmark). In one test of
  `homer2 norm` and `known` at `-nlen 15 -p 1` on a small sequence set, v5.1 peaked at 58 GB and
  homer2-ao at 50 GB.

## The pip package

`pip install` and `uv pip install` configure a copy of the tree while building the wheel, in
place of `configureHomer.pl -make`, so that the installed package works wherever it lands:

- Each script configureHomer.pl configures gets the same first four lines, except that lines 3
  and 4 (`use lib` and `$homeDir`) take the script's own directory at run time instead of an
  absolute path. Nine more scripts carry the same two lines with the original author's path
  (`createIGVhtml.pl`, `filterTADsAndCPs.pl`, `findcsRNATSR.pl`, `flipPC1toMatch.pl`,
  `homer2jaspar.pl`, `HomerSeqBias.pm`, `makeHiCWashUfile.pl`, `selectRepeatBg.pl`,
  `subtractBedGraphsDirectory.pl`); configureHomer.pl leaves them as they are, and the package
  rewrites them too.
- `cpp/SeqTag.cpp`, whose first lines configureHomer.pl rewrites with the install path, finds the
  installation from the running program's path (`/proc/self/exe` on Linux,
  `_NSGetExecutablePath` on macOS), or from `HOMER2_AO_HOME`, which the package's commands set.
- The compiled programs have their debug information stripped. On Linux other than x86-64 they
  are compiled with `-fsigned-char`, because HOMER's code assumes a signed `char`
  (`docs/V5.1-DEFECTS.md`, entry 63); on macOS with `-std=gnu++17`. The x86-64 Linux build uses
  the Makefile's flags unchanged.
- The package ships the `config.txt` that `configureHomer.pl -install` writes for a fresh tree
  as `config.default.txt`, and its commands copy it to `config.txt` on first use. pip and uv
  therefore do not own `config.txt`, so an upgrade keeps it and the genomes it lists.
- The package leaves out `update/`, `bin/old/`, `cpp/backup/` and the Doughnut documentation PDF.

## Outside the optimization

- `bin/HomerSVGLogo.pm` calls `plotXY` as a function where v5.1 wrote `$plotXY(...)`, which Perl
  5.44 does not compile (`docs/V5.1-DEFECTS.md`, entry 62). Nothing calls that function, so no
  output changes.

## Building

- `cpp/Motif2.cpp` includes `<string_view>`, and the Makefile sets no `-std`, so the compiler's
  default must be C++17 or later (GCC 11 and later). It uses `__builtin_prefetch`,
  `__builtin_popcount` and `vector_size`, which GCC and Clang provide.
- `compareMotifMatrices.cpp` is compiled with `-ffp-contract=off`, so its scores match Perl's
  arithmetic bit for bit.
- The build uses no CPU-specific or fast-math flags. Continuous integration builds it and runs
  the tests on x86-64 and ARM Linux (GCC 13) and on macOS on Apple silicon (Apple clang 21). The
  byte-for-byte comparisons with v5.1 in `docs/VALIDATION.md` were run on x86-64 Linux.

# homer2-ao

> **If you use this software, you MUST cite the HOMER paper:**
>
> Heinz S, Benner C, Spann N, Bertolino E, Lin YC, Laslo P, Cheng JX, Murre C, Singh H, Glass CK.
> Simple combinations of lineage-determining transcription factors prime cis-regulatory elements
> required for macrophage and B cell identities. *Molecular Cell* 38(4):576-589 (2010).
> https://doi.org/10.1016/j.molcel.2010.05.004
>
> HOMER and its documentation: http://homer.ucsd.edu/homer/

homer2-ao is HOMER v5.1 with faster motif scanning and motif finding. `homer2` (`find`, `known`,
`denovo`, `norm`, `mask`), `findMotifsGenome.pl` and `compareMotifs.pl` give the same output files
as HOMER v5.1, byte for byte, in a small fraction of the time. Every HOMER program keeps its name
and options, so existing commands and pipelines run unchanged. This optimization was found with
auto-optimize, a tool developed by the Programmable Genomics Laboratory at UMass Chan Medical
School.

| Command (K562 ATAC-seq peaks, hg38, `-p 4`) | HOMER v5.1 | homer2-ao |
|---|---:|---:|
| `homer2 find`, 269,800 peaks x 472 known motifs | 603.6 s, 5.6 GB | 1.1 s, 0.5 GB |
| `findMotifsGenome.pl`, 10,000 peaks, `-N 50000` | 769 s | 17 s |

How these were measured, with the other cases and the time per step, is in
`docs/PERFORMANCE.md`.

## Installation

```
uv pip install homer2-ao        # or: pip install homer2-ao
uv tool install homer2-ao       # the commands on PATH, in an environment of their own
```

HOMER's scripts need `perl`. On x86-64 Linux with glibc 2.27 or later (most distributions
released since 2018), the package comes with HOMER's programs already compiled. On other Linux
systems it compiles them during installation (macOS has not been tested), which needs `make` and
a C++ compiler whose default language standard is C++17 or later (GCC 11 or later, or Clang).
Every HOMER program and script is installed as a command of the same name. From a clone of this
repository, `uv pip install .` compiles and installs it the same way.

The repository is also a complete HOMER tree, and can be built the way HOMER is, in place:

```
perl configureHomer.pl -make
```

## Usage

The commands and options are HOMER v5.1's; HOMER's documentation covers them. For example:

```
findMotifsGenome.pl peaks.bed hg38.fa motifs_out -size 200 -p 8
homer2 find -i peaks.fa -m known.motifs -p 8 -o sites.txt
```

A genome can be given as a FASTA file, as above, or installed by name into the package's HOMER
directory with `configureHomer.pl -install hg38`. `homer2-ao` prints the version and that
directory, and `homer2-ao <program> [args]` runs a HOMER program by name. From Python,
`homer2_ao.home()` and `homer2_ao.bin_dir()` give the same locations.

## What changed

Nine commits on top of v5.1, merged into `main`:

1. `homer2 denovo -Tcpus`: round the trial count for a given CPU count.
2. `compareMotifMatrices`, a C++ program that scores motif pairs as `compareMotifs.pl` does.
3. `compareMotifs.pl` calls it.
4. `homerTools bgbins`, `randbg` and `filterids` for `findMotifsGenome.pl`'s background
   selection.
5. `homer2 find`: an indexed motif scan, output in `-p 1` order at any `-p`, FASTA read in
   memory.
6. `homer2 known`: ZOOPS counts through the index, a window filter, a thread-safe score cache.
7. `homer2 norm`: autonormalization codes found once, on `-p` threads.
8. `homer2 denovo`: a faster oligo table, mismatch screen, scoring and mask.
9. `findMotifsGenome.pl`: the longest de novo length runs beside the rest.

Each commit message says what changes, why the output stays the same, and which functions it
changes. Every way homer2-ao's behaviour differs from v5.1's (new options, crash behaviour on bad
input, stderr, build requirements) is in `docs/DIFFERENCES.md`. HOMER's other programs are
unchanged except where they share this code.

## Tests

`tests/` holds 741 pytest tests of the motif programs, which pass on v5.1 and on homer2-ao, built
in place or installed as a package:

```
HOMER_HOME=$(homer2-ao --home) python -m pytest tests     # an installed package
perl configureHomer.pl -make && touch config.txt && python -m pytest tests     # in place
```

`HOMER_HOME`, `HOMER_TEST_BIN` and `HOMER_TEST_JOBS` choose the install, the programs and the
number of parallel program calls (`tests/homer_helpers.py`). How homer2-ao's output was compared
with v5.1's beyond these tests is in `docs/VALIDATION.md`.

## Benchmarks

`benchmarks/` holds the two benchmarks behind `docs/PERFORMANCE.md`, `homer2 find` and
`findMotifsGenome.pl`, with the scripts that build their inputs and v5.1 references; see
`benchmarks/README.md`.

## Defects in v5.1

The tests and comparisons found 61 defects in HOMER v5.1 itself, listed with reproductions in
`docs/V5.1-DEFECTS.md`. homer2-ao fixes 3 of them and behaves differently in 4 (also in
`docs/DIFFERENCES.md`); 41 are unchanged, most in code homer2-ao does not touch, and 13 were not
checked.

## License

HOMER is written by Christopher Benner and distributed under the GNU General Public License,
version 3 or later. homer2-ao is a modified version of HOMER v5.1, changed in October 2026, under
the same license (`COPYING`). The changes are the commits on the `optimize` branch, each dated.
The package carries HOMER's C++ source alongside the compiled programs.

The `main` branch holds the v5.1 release (tag `v5.1`, without its 15 prebuilt programs), the
tests, the merge of `optimize`, the packaging (`pyproject.toml`, `hatch_build.py`, `homer2_ao/`)
and these documents.

# Performance

Every number here is a wall-clock time measured on one machine: an AMD EPYC 9575F, 2 x 64 cores
with SMT (256 logical CPUs), shared with other jobs, so the load from those jobs is given where
it was recorded. Builds are `configureHomer.pl -make` with g++ 13.3 and the release's `-O3`. The
benchmarks themselves are in `benchmarks/`; each compares every output byte for byte with HOMER
v5.1 run at `-p 1` before a time counts.

## homer2 find

`homer2 find -i <peaks.fa> -m <motifs> -p N -o out.txt` on K562 ATAC-seq peaks (ENCODE
ENCFF333TAT; ENCFF855PCP for the 172,732-peak set), 200 bp around each summit on hg38, against
HOMER's 472 known vertebrate motifs. Whole-process wall time including input and output, with
peak RSS from `/usr/bin/time`, measured 2026-10-12 with the homer2 find changes applied alone to
v5.1. The v5.1 calls ran side by side with one another (medians of 3 at 5,000 and 26,980 peaks,
single runs at the two larger sizes). The changed homer2 ran alone, one call at a time, median
of 5 after a warm-up. Every output of the changed homer2 had the same SHA-256 as v5.1's at
`-p 1`.

| Peaks | v5.1 `-p 1` | v5.1 `-p 4` | homer2-ao `-p 1` | homer2-ao `-p 4` | Peak RSS at `-p 4`, v5.1 / homer2-ao |
|---:|---:|---:|---:|---:|---:|
| 5,000 | 47.1 s | 11.76 s | 0.072 s | 0.051 s | 95 / 28 MB |
| 26,980 | 240.6 s | 62.2 s | 0.409 s | 0.121 s | 536 / 84 MB |
| 172,732 | 1,609.9 s | 396.7 s | 2.27 s | 0.571 s | 3,533 / 350 MB |
| 269,800 | 2,428.5 s | 603.6 s | 3.81 s | 1.105 s | 5,579 / 544 MB |

The later commits change code that homer2 find shares (the window filter, the score cache, the
known ZOOPS path). On 2026-10-28 the released code was timed against the homer2 find changes alone
with paired runs (24 pairs per case, alternating which build ran first, sign test):

| Case | Final minus find-only, median | Final faster in | p |
|---|---:|---:|---:|
| 269,800 peaks, `-p 4` | -0.046 s | 17 of 24 | 0.064 |
| 5,000 peaks, `-p 4` | -0.0004 s | 13 of 24 | 0.84 |

Neither difference is resolved, so the table above stands for the released code. In one checked
benchmark run each on the same day, the 269,800-peak call took 0.915 s, 3.49 CPU s and 537 MB with the final
code, and 1.044 s, 3.96 CPU s and 519 MB with the find changes alone; both gave output identical
to v5.1 at `-p 1` for find, `-mscore` and mask on every case.

## findMotifsGenome.pl

`findMotifsGenome.pl peaks.bed hg38.fa out -size 200 -N 50000 -p 4 -preparsedDir <dir>` on the
10,000 highest-signal peaks of ENCFF333TAT, 200 bp at the summit, with the preparsed background
already present. HOMER's defaults apply: `-len 8,10,12`, `-S 25` and the 1,049 known motifs of
`-mset auto`. Every output file, stdout and the exit status are compared with v5.1 run at
`-p 1`, with `homer2 denovo` given the trial count the case's `-p` implies (`-T 12` at `-p 4`).

Paired runs, nine pairs per case, alternating which build ran first (2026-10-25):

| Case | v5.1, median | homer2-ao, median | homer2-ao faster in | p |
|---|---:|---:|---:|---:|
| 10,000 peaks, `-N 50000`, `-p 4` | 769.19 s | 16.94 s | 9 of 9 | 0.0039 |
| 2,000 peaks, `-N 10000`, `-p 16` | 161.91 s | 8.19 s | 9 of 9 | 0.0039 |

These pairs were measured on a development version that differs from the released code only in
comments and in two error-path fixes (`docs/VALIDATION.md`). The released code took 18.21 s on the
first case (median of 3 runs) against 17.64 s for that version, within the 1.14 s spread between
repeated runs. Built again from the release zip and measured at a higher machine load (61 to 112),
v5.1 took 851.79 s and that version 20.754 s (median of 3 each).

Per step, for the same development version (medians of 3 runs):

| Step | v5.1 | homer2-ao |
|---|---:|---:|
| Background selection | 27.21 s | 3.64 s |
| homer2 norm (autonormalization) | 13.16 s | 1.248 s |
| Known motifs | 220.7 s | 1.069 s |
| homer2 denovo, `-len 8` | 42.33 s | 2.484 s |
| homer2 denovo, `-len 10` | 119.7 s | 7.618 s |
| homer2 denovo, `-len 12` | 314.4 s | 12.09 s |
| compareMotifs.pl | 22.43 s | 0.6462 s |
| Cold `preparseGenome.pl` of hg38 (not changed) | 37.69 s | 38.28 s |
| CPU time of the whole call | 2,456 s | 41.11 s |
| Peak RSS | 1,977 MB | 1,882 MB |

homer2-ao runs known motifs, `-len 8` and `-len 10` in one lane beside `-len 12`, so those three
rows are not step-for-step comparable with v5.1's sequential steps.

### By commit

The benchmark run once on each prefix of the commits (v5.1 plus every commit up to and including
that one), one call at a time under a 16-core cap, 2026-10-28, at loads of 30 to 176. These are
single runs and show where the time goes, not paired comparisons. In these runs the known-motif
step took 236.7 s before the find commit, 21.0 s after it (`homer2 known` scans with the new
tables) and 0.39 s after the next commit (ZOOPS counts through the index).

| Commit | `-p 4` case |
|---|---:|
| homer2 denovo: add -Tcpus | 867.3 s |
| compareMotifMatrices | 852.2 s |
| compareMotifs.pl calls compareMotifMatrices | 855.0 s |
| homerTools: bgbins, randbg and filterids | 782.3 s |
| homer2 find: indexed motif scan | 477.5 s |
| homer2 known: index ZOOPS counts, window filter, score cache | 416.9 s |
| homer2 norm: codes once, on -p threads | 359.1 s |
| homer2 denovo: oligo table, mismatch screen, scoring, mask | 24.0 s |
| findMotifsGenome.pl: longest de novo length beside the rest | 19.3 s |

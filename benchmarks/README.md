# Benchmarks

Two benchmarks time a HOMER tree and check that its output is byte-identical to HOMER v5.1's.
Each needs its own reference directory: HOMER v5.1's outputs on fixed inputs, made once with
the scripts here. Results of these benchmarks are in `docs/PERFORMANCE.md`.

| Directory | Times | Checks against HOMER v5.1 |
|---|---|---|
| `find/` | `homer2 find -p 4` on 269,800 K562 ATAC-seq peaks (200 bp) x 472 known motifs | sorted lines and -p 1 line order of `find`, `find -mscore` and `mask` on 15 fixed cases and 4 on a per-seed input |
| `findMotifsGenome/` | `findMotifsGenome.pl` on the 10,000 strongest of those peaks against 50,000 background sequences, `-p 4` | every output file, stdout and exit status of four cases (`bench.py`), and of 32 small cases of HOMER's other programs plus two groups of `homer2 denovo` runs (`compare.py`) |

## Requirements

Linux with `/dev/shm`, g++ and make, Perl, `samtools` on `PATH`, and Python 3. numpy is needed by
`find/prepare_small_inputs.py`, `find/bench.py`, `find/make_small_refs.py` and
`findMotifsGenome/make_inputs.py`.

## Inputs

| Input | Source | Used as |
|---|---|---|
| hg38 | UCSC `hg38.fa.gz`, https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz, uncompressed and indexed with `samtools faidx` (455 sequences, soft-masked, 50 bp lines; md5 of the uncompressed file used: `b2aee9f885accc00531e59c4736bee63`) | genome for every input |
| ENCFF333TAT | K562 ATAC-seq, pooled pseudoreplicated peaks, https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz, unmodified (269,800 lines; md5 of the decompressed text `c6b01945e6e5dd27d7903326f5e0dce3`) | both benchmarks |
| ENCFF855PCP | K562 ATAC-seq, conservative IDR peaks, https://www.encodeproject.org/files/ENCFF855PCP/@@download/ENCFF855PCP.bed.gz, unmodified (172,732 lines; `92fc14129cbc446eac8dfbf8462fa76e`) | `find/`: idr172k.fa |
| medium peaks | every 10th line of ENCFF333TAT, coordinate-sorted: `zcat ENCFF333TAT.bed.gz \| awk 'NR % 10 == 0' \| LC_ALL=C sort -k1,1 -k2,2n \| gzip > medium.narrowPeak.gz` (26,980 lines; `94a4d3f818907c7a2d5528927a1cf30a`) | `find/`: the small inputs |
| ENCFF901CBP | K562 CTCF ChIP-seq, conservative IDR peaks, https://www.encodeproject.org/files/ENCFF901CBP/@@download/ENCFF901CBP.bed.gz, sorted with `LC_ALL=C sort -k1,1 -k2,2n -k3,3n` (56,067 lines; `d8afbe84058157e06a25c1c95881465a`) | `findMotifsGenome/`: ctcf10k.bed (used only by the `ctcf` case, which `bench.py` does not run) |
| known motifs | HOMER's own `data/knownTFs` and `motifs/ctcf.motif` | both benchmarks |

The md5 sums are of the files as these benchmarks' references were made. The find inputs keep
the peak file's line order, and the findMotifsGenome peak sets are ranked by signalValue with
ties in file order, so a file in a different line order gives different inputs.

## Environment variables

| Variable | Read by | Meaning |
|---|---|---|
| `HOMER_FIND_REF` | `find/bench.py`, `find/time_find.py`, `find/make_small_refs.py` | the find reference directory: `inputs/`, `outputs/` (sorted), `order_p1/` (-p 1 order) |
| `HOMER_FIND_INPUTS` | the same, optional | a directory holding `inputs/`, when it is not the reference directory |
| `HOMER_FMG_REF` | every script in `findMotifsGenome/` | the findMotifsGenome reference directory: `inputs/`, `preparsed/`, `cases/`, `v51/homer` (the HOMER v5.1 build), `denovo80/`, `checks/` |
| `HOMER_BENCH_GENOME` | every script in `findMotifsGenome/` | the hg38 FASTA |
| `HOMER_FMG_CHECK_INPUTS` | `findMotifsGenome/make_inputs.py`, `compare.py`, optional | `compare.py`'s small inputs (default `$HOMER_FMG_REF/checks/inputs`) |
| `HOMER_FMG_WRITE_REFERENCE` | `findMotifsGenome/compare.py` | `1` caches the given tree's outputs as the reference instead of comparing |

`findMotifsGenome.pl` writes its own command line into `motifFindingParameters.txt`, so the
genome path, the preparsed directory and `compare.py`'s input paths are part of the reference
outputs. Run each benchmark with the same paths its references were made with.

## homer2 find

`find/bench.py` builds homer2 from `<tree>/cpp` with `find/Makefile` (homer2 only, v5.1's
flags), times one `-p 4` call on big270k after an identical warmup, then runs the small cases
and the seed's extra input (300 peaks and 60 of 1,049 motifs drawn by `--seed`) and checks every
output. It writes `metrics.json` (`objective`, `checks`, `measurements`) to `--out`; each check
counts mismatches against v5.1 and is 0 when every output matches. References
exist for seeds 0-9 unless `make_small_refs.py` is given another range.

```
export HOMER_FIND_REF=<find reference>
python3 benchmarks/find/prepare_small_inputs.py --out $HOMER_FIND_REF/inputs \
	--genome hg38.fa --peaks medium.narrowPeak.gz --homer <HOMER v5.1 tree>
python3 benchmarks/find/prepare_large_inputs.py --out $HOMER_FIND_REF/inputs \
	--genome hg38.fa --peaks ENCFF333TAT.bed.gz --idr-peaks ENCFF855PCP.bed.gz
make -f benchmarks/find/Makefile -C <HOMER v5.1 tree> OUT=<v5.1 build>
python3 benchmarks/find/make_refs.py --homer2 <v5.1 build>/homer2 \
	--fasta $HOMER_FIND_REF/inputs/big270k.fa --motifs $HOMER_FIND_REF/inputs/known472.motifs \
	--case big270k_find_p4 --ref $HOMER_FIND_REF --groups 48 --jobs 24
python3 benchmarks/find/make_small_refs.py --homer2 <v5.1 build>/homer2 --seeds 0-9 --jobs 4

python3 benchmarks/find/bench.py --tree <HOMER tree> --seed 0 --out <dir> --cache <dir>
python3 benchmarks/find/time_find.py --src <HOMER tree> --order    # one timed call, seconds last
```

`make_refs.py` relies on v5.1 writing `find`'s lines motif by motif at `-p 1`: it runs the
motif file in consecutive groups in parallel and concatenates them. `--verify <sha256 file>`
checks the result against a whole-file `-p 1` output first. `time_find.py` builds into
`<tree>/.paired-cache/` and prints one call's seconds last, for an A/B timer that alternates
two trees; its `idr_` and `all_` cases need `validation_order_p1/` and `validation_outputs/`
references that no script here writes.

## findMotifsGenome.pl

`findMotifsGenome/bench.py` copies the tree into `--cache`, builds it with
`configureHomer.pl -make`, reads the genome and preparsed files once, then times the objective
call and checks it and three more cases: 500 peaks at `-p 1`, 2,000 peaks at `-p 16`, and a
cold `preparseGenome.pl`. Every call runs with Perl's hash order and `rand()` fixed
(`harness/perl/BenchSrand.pm`), so outputs repeat exactly; the references give `homer2 denovo`
the trial count the case's `-p` gives (v5.1 rounds `-T 10` up to a multiple of `-p`).
`compare.py` builds the tree in place and runs 32 small cases of HOMER's other programs plus
the de novo check (80 subsamples, and four of them at `-p 1`, 4 and 16), on chr21 and chr22.

```
export HOMER_FMG_REF=<findMotifsGenome reference> HOMER_BENCH_GENOME=<hg38.fa>
python3 benchmarks/findMotifsGenome/make_refs.py inputs \
	--atac ENCFF333TAT.bed.gz --ctcf ENCFF901CBP.sorted.narrowPeak.gz
python3 benchmarks/findMotifsGenome/make_refs.py build <HOMER v5.1 tree>
python3 benchmarks/findMotifsGenome/make_refs.py preparse
for case in objective p1 p16; do python3 benchmarks/findMotifsGenome/make_refs.py case $case; done
python3 benchmarks/findMotifsGenome/make_inputs.py --atac ENCFF333TAT.bed.gz
HOMER_FMG_WRITE_REFERENCE=1 python3 benchmarks/findMotifsGenome/compare.py \
	<copy of the HOMER v5.1 tree> $HOMER_FMG_REF/checks/outputs
python3 benchmarks/findMotifsGenome/make_denovo80.py --cases <subsample directory>

python3 benchmarks/findMotifsGenome/bench.py --tree <HOMER tree> --seed 0 --out <dir> --cache <dir>
python3 benchmarks/findMotifsGenome/compare.py <HOMER tree> $HOMER_FMG_REF/checks/outputs
python3 benchmarks/findMotifsGenome/time_one.py --src <HOMER tree>    # one timed call, seconds last
```

`compare.py` and `time_one.py` build the tree in place, and `configureHomer.pl -make` rewrites
the paths in HOMER's scripts, so give them an exported copy (`git archive`), not a working
checkout. The de novo check's 80 subsamples of the objective's `homer2 denovo` input were made
by a script that is not in this directory; `make_denovo80.py --cases` takes a directory with
their `tri.log` and `cases/s<seed>/seq.tsv` and `group.adj`, and `compare.py` stops at the
de novo check when `$HOMER_FMG_REF/denovo80` is missing.

## Run cost

Measured on a shared 256-core server, so approximate:

- find references: v5.1 `homer2 find -p 1` on all 269,800 peaks took 2,428.5 s as one call;
  `make_refs.py` spreads it over `--jobs` processes. `make_small_refs.py --seeds 0-1 --jobs 4`
  (22 cases) took 67 s. The reference directory takes about 700 MB.
- findMotifsGenome references at `-p 1`: objective 1,832 s, p1 318 s, p16 773 s, cold
  preparse 42 s. The preparsed directory takes 1.2 GB.
- One run on the optimized tree: find `bench.py` 20 s, findMotifsGenome `bench.py` 90 s and
  `compare.py` 60 s, under an 8-core cap. On a tree running at v5.1's speed, one
  findMotifsGenome `bench.py` run took 1,415 s and `compare.py` 304 s.

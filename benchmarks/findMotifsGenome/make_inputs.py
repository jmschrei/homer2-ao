"""Build compare.py's small inputs once, from public data.

	HOMER_FMG_REF=<reference> HOMER_BENCH_GENOME=<hg38.fa> python3 make_inputs.py \
		--atac <ENCFF333TAT .narrowPeak.gz>

Writes into $HOMER_FMG_CHECK_INPUTS (default <reference>/checks/inputs) and
reads HOMER v5.1's motif files from <reference>/v51/homer (make_refs.py
build first). mini.fa: hg38 chr21 and chr22 as they are in the genome
(soft-masked, 50 bp lines). Peaks: ENCODE ENCFF333TAT K562 ATAC-seq on
chr21/chr22, ranked by signalValue, 200 bp at the summit. Background: 1,000
positions drawn with a fixed numpy seed from non-N stretches. Motifs: every
9th of HOMER's 472 vertebrate known motifs, and HOMER's ctcf.motif. Edge
cases: an empty peak file, one peak, and peaks at the chromosome ends, in an
all-N region, on a chromosome the genome lacks, and with a repeated name.
"""

import os
import sys
import gzip
import argparse

import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fmg

OUT = fmg.CHECK_INPUTS
HG38 = fmg.GENOME
ATAC = None  # --atac
V51 = fmg.V51
CHROMS = ('chr21', 'chr22')


def read_chrom(name):
	with open(HG38 + '.fai') as handle:
		for line in handle:
			cols = line.split('\t')
			if cols[0] == name:
				length, offset, width, stride = map(int, cols[1:5])
				break
	rows = (length + width - 1) // width
	with open(HG38, 'rb') as handle:
		handle.seek(offset)
		data = handle.read(rows * stride)
	seq = data.replace(b'\n', b'')[:length]
	assert len(seq) == length
	return seq


def fasta(path, records):
	with open(path, 'w') as handle:
		for name, seq in records:
			handle.write('>{}\n{}\n'.format(name, seq))


def main():
	os.makedirs(OUT, exist_ok=True)
	genome = {name: read_chrom(name) for name in CHROMS}
	with open(os.path.join(OUT, 'mini.fa'), 'wb') as handle:
		for name in CHROMS:
			handle.write(b'>' + name.encode() + b'\n')
			seq = genome[name]
			for i in range(0, len(seq), 50):
				handle.write(seq[i:i + 50] + b'\n')

	peaks = []
	with gzip.open(ATAC, 'rt') as handle:
		for line in handle:
			cols = line.split('\t')
			if cols[0] in CHROMS:
				peaks.append((cols[0], int(cols[1]) + int(cols[9]), float(cols[6])))
	peaks.sort(key=lambda row: -row[2])

	def bed(rows, prefix, first=1):
		return ['{}\t{}\t{}\t{}{}\t0\t+\n'.format(chrom, summit - 100, summit + 100,
			prefix, rank) for rank, (chrom, summit, _) in enumerate(rows, first)]

	top = bed(peaks[:200], 'p')
	with open(os.path.join(OUT, 'peaks.bed'), 'w') as handle:
		handle.writelines(top)
	with open(os.path.join(OUT, 'peaks2.bed'), 'w') as handle:
		handle.writelines(bed(peaks[150:350], 'q', 151))

	rng = numpy.random.RandomState(20261024)
	background = []
	while len(background) < 1000:
		chrom = CHROMS[rng.randint(2)]
		start = int(rng.randint(0, len(genome[chrom]) - 200))
		seq = genome[chrom][start:start + 200]
		if seq.upper().count(b'N') == 0:
			background.append((chrom, start + 100, 0.0))
	with open(os.path.join(OUT, 'bg.bed'), 'w') as handle:
		handle.writelines(bed(background, 'b'))

	def seqs(rows, prefix):
		return [('{}{}'.format(prefix, rank), genome[chrom][summit - 100:summit + 100]
			.decode()) for rank, (chrom, summit, _) in enumerate(rows, 1)]

	targets, bgseqs = seqs(peaks[:200], 'p'), seqs(background, 'b')
	fasta(os.path.join(OUT, 't200.fa'), targets)
	fasta(os.path.join(OUT, 'b1k.fa'), bgseqs)
	with open(os.path.join(OUT, 'seq.tsv'), 'w') as handle:
		for name, seq in targets + bgseqs:
			handle.write('{}\t{}\n'.format(name, seq))
	with open(os.path.join(OUT, 'group.tsv'), 'w') as handle:
		for name, _ in targets:
			handle.write('{}\t1\t1\n'.format(name))
		for name, _ in bgseqs:
			handle.write('{}\t0\t1\n'.format(name))

	motifs, current = [], None
	with open(os.path.join(V51, 'data', 'knownTFs', 'vertebrates',
			'known.motifs')) as handle:
		for line in handle:
			if line.startswith('>'):
				current = [line]
				motifs.append(current)
			elif current is not None:
				current.append(line)
	with open(os.path.join(OUT, 'known_sub.motifs'), 'w') as handle:
		for motif in motifs[::9]:
			handle.writelines(motif)
	with open(os.path.join(V51, 'motifs', 'ctcf.motif')) as src, open(
			os.path.join(OUT, 'ctcf.motif'), 'w') as dst:
		dst.write(src.read())

	open(os.path.join(OUT, 'empty.bed'), 'w').close()
	with open(os.path.join(OUT, 'one.bed'), 'w') as handle:
		handle.write(top[0])
	n21, n22 = len(genome['chr21']), len(genome['chr22'])
	with open(os.path.join(OUT, 'edge.bed'), 'w') as handle:
		handle.writelines(top[:40])
		handle.write('chr21\t0\t200\tedge_start\t0\t+\n')
		handle.write('chr21\t{}\t{}\tedge_end\t0\t-\n'.format(n21 - 150, n21 + 50))
		handle.write('chr22\t{}\t{}\tedge_end22\t0\t+\n'.format(n22 - 200, n22))
		handle.write('chr21\t1000000\t1000200\tall_n\t0\t+\n')
		handle.write('chrUn_missing\t5000\t5200\tno_chrom\t0\t+\n')
		handle.write(top[50].replace('\tp51\t', '\tp1\t'))
	assert genome['chr21'][1000000:1000200].upper().count(b'N') == 200
	print('inputs written to', OUT)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--atac', required=True)
	ATAC = parser.parse_args().atac
	main()

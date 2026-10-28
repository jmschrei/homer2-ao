#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";


my $maxCPUs = 50;
my $motifFile = '';
my $date = '000000';
my $outDir = "/homer_data/www/html/homer/data/motifs/";
my $filePrefix = "homer.KnownMotifs";

sub printCMD {
	print STDERR "\n\tusage: updateMotifScanningTracks.pl [options] -date <YYMMDD> -genomes <genome1> [genome2] ...\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-date <YYMMDD> (identifier to apend to result filenames)\n";
	print STDERR "\t\t-genomes <genome1> [genome2] ... (genomes to scan)\n";
	print STDERR "\t\t-m <motifFile> (motif file to scan)\n";
	print STDERR "\t\t-p <#> (number of CPUs to use, def: $maxCPUs)\n";
	print STDERR "\t\t-dir <output directory> (def: $outDir)\n";
	print STDERR "\t\t-prefix <output file prefix> (def: $filePrefix)\n";
	print STDERR "\n\tNote that this script will create temporary files in your current working directory.\n";
	print STDERR "\tFor each genome a gzipped bed file of motif positions and a bigBed file will be produced.\n\n";
	exit;
}


if (@ARGV < 1) {
	printCMD();
}
my %genomes = ();
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-date') {
		$date = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-genomes' || $ARGV[$i] eq '-genome') {
		$i++;
		while ($i < @ARGV && $ARGV[$i] !~ /\-/) {
			$genomes{$ARGV[$i]} = 1;
			$i++;
		}	
		if ($ARGV[$i] =~ /\-/) {
			$i--;
			next;
		}
	} elsif ($ARGV[$i] eq '-p') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-m') {
		$motifFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-dir') {
		$outDir = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-prefix') {
		$filePrefix = $ARGV[++$i];
	} else {
		printCMD();
	}
}

if ($motifFile eq '') {
	print STDERR "!!! Error: Motif file (-m) required!!!\n";
	exit;
}
my @genomes = sort {$a cmp $b} keys %genomes;
if (@genomes < 1) {
	print STDERR "!!! Error: No genomes specified!!!\n";
	exit;
}

print STDERR "\n\tScanning script options:\n";
print STDERR "\t\tdate: $date\n";
print STDERR "\t\tCPUs: $maxCPUs\n";
print STDERR "\t\tmotif file: $motifFile\n";
print STDERR "\t\toutput directory: $outDir\n";
print STDERR "\t\toutput file prefix: $filePrefix\n";
print STDERR "\t\tGenomes to scan:\n";
foreach(@genomes) {
	print STDERR "\t\t\t$_\n";
}
print STDERR "\n\tSleeping for 5 seconds before starting...\n";
`sleep 5`;

my $rand = rand();
my $tmpBed = $rand . ".bed";
my $tmpBedSorted = $rand . ".sorted.bed";
my $tmpBedSortedZip = $rand . ".sorted.bed.gz";
my $tmpBigBed = $rand . ".bigBed";

foreach(@genomes) {
	my $genome = $_;

	my $genomenor = $genome;
	$genomenor =~ s/r$//;
	my $sizes = $homeDir . "/data/genomes/" . $genomenor . "/chrom.sizes";
	unless (-e $sizes) {
		print STDERR "!!! Error - could not find chrom.sizes file ($sizes)\n";
		next;
	}
	my $fname = $filePrefix . "." . $genome . "." . $date . ".bigBed";
	my $out1 = $outDir . "/" . $filePrefix . "." . $genome . "." . $date . ".bigBed";
	my $out2 = $out1 . ".track.txt";
	my $out3 = $outDir . "/" . $filePrefix . "." . $genome . "." . $date . ".bed.gz";
	print STDERR "\n\tAnalyzing genome: $genome\n";
	print STDERR "\t\tchromSizes file: $sizes\n";
	print STDERR "\t\tbb Filename  : $fname\n";
	print STDERR "\t\tBigBed output: $out1\n";
	print STDERR "\t\tBigBed tracks: $out2\n";
	print STDERR "\t\tgzipped bed:   $out3\n";

	`scanMotifGenomeWide.pl "$motifFile" $genome -bed -int -keepAll -p $maxCPUs > "$tmpBed"`;
	`sort -k1,1 -k2,2n "$tmpBed" > "$tmpBedSorted"`;
	`bedToBigBed "$tmpBedSorted" "$sizes" "$tmpBigBed"`;
	`mv "$tmpBigBed" "$out1"`;

	open OUT, ">$out2";
	print OUT "track type=bigBed name=\"HOMER Known Motifs $genome ($date)\" description=\"HOMER Known Motifs $genome ($date)\" bigDataUrl=http://homer.ucsd.edu/homer/data/motifs/$fname visibility=3\n";
	close OUT;

	`gzip "$tmpBedSorted"`;
	`mv "$tmpBedSortedZip" "$out3"`;
}


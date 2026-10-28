#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";

use HomerConfig;
use POSIX;

my $config = HomerConfig::loadConfigFile();


my $maxN = 0.7;
my $rmZeroFlag = 0;
my $mask = '';
my $size = '';

sub printCMD {
	print STDERR "\n\tNew usage: createScoredFasta.pl <peak/BED> <genome> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-size <#> or -size <#,#> (resize peak/BED file, default: $size)\n";
	print STDERR "\t\t-mask (use masked version of genome)\n";
	print STDERR "\t\t-maxN <#> (maximum fraction of sequence that can be masked, def: $maxN)\n";
	print STDERR "\t\t-removeZeros (remove sequences with score equal to '0')\n";

	print STDERR "\n";
	exit;
}

my $peakFile = '';
my $genome = '';

if (@ARGV < 2) {
	printCMD();
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-mask') {
		$mask = " -mask ";
	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-i') {
		$peakFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-removeZeros') {
		$rmZeroFlag = 1;
	} elsif ($ARGV[$i] eq '-maxN') {
		$maxN = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-size') {
		$size = "-size " . $ARGV[++$i];
	} elsif ($i == 0) {
		$peakFile = $ARGV[$i];
	} elsif ($i == 1) {
		$genome = $ARGV[$i];
	} else {
		printCMD();
	}
}

if ($peakFile eq '') {
	print STDERR "!!! peak file must be specified (argument 1)\n";
	exit;
}
if ($genome eq '' || $genome eq 'none') {
	print STDERR "!!! genome must be specified (argument 2)\n";
	exit;
}
if (-e $genome) {
	print STDERR "\tUsing genome found in $genome\n";
} else {
	if ($genome =~ s/r$//) {
		$mask = "-mask";
	}
	if (exists($config->{'GENOMES'}->{$genome})) {
		$genome = $config->{'GENOMES'}->{$genome}->{'directory'};
	} else {
		print STDERR "!!! Error - genome must be a valid FASTA file, directory of chromosome FASTA files, or configured HOMER genome\n";
		exit;
	}
}

my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";

`adjustPeakFile.pl "$peakFile" $size > "$tmpFile"`;

my $seqLen = $size;
if ($size =~ /\,/) {
	my @a = split /\,/, $size;
	$seqLen = $a[1]-$a[0];
}
print STDERR "\tSequence length will be $seqLen\n";


my $totalSize = 0;
my $totalN = 0;
my $total=0;
my $zeroCount=0;
my %data = ();
open IN, $tmpFile;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line = split /\t/;
	next if (scalar(@line) < 6);
	my $id = $line[0];
	my $v = $line[5];
	$total++;
	if ($rmZeroFlag) {
		if ($v eq '0' || $v eq 'na' || $v eq 'NA') {
			$zeroCount++;
			next;
		}
	}
	my $len = $line[3]-$line[2]+1;
	$totalSize += $len;
	$totalN+=1.0;
	$data{$id} = $v;
}
close IN;

my $avgSize = $totalSize/$totalN;
my $lenTarget = floor($avgSize+0.5);
print STDERR "\n\tTotal scored regions: $total\n";
print STDERR "\tTotal with zero-scores removed: $zeroCount\n\n";
print STDERR "\tAverage length of regions: $avgSize\n";

if (0) {
open OUT, ">covars.txt";
open IN, $ARGV[0];
my $z =0;
while (<IN>) {
	$z++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($z == 1 && $line[1] =~ /tart/) {
		print OUT "peakID";
		for (my $i=6;$i<@line;$i++)  {
			print OUT "\t$line[$i]";
		}
		print OUT "\n";
		next;
	}
	next if (/^#/);
	$data{$line[3]} = $line[4];
	print OUT "$line[3]";
	for (my $i=6;$i<@line;$i++) {
		print OUT "\t$line[$i]";
	}
	print OUT "\n";
}
close IN;
}

`homerTools extract "$tmpFile" "$genome" $mask > "$tmpFile2"`;
$total=0;
my $removedN=0;
my $removed1=0;
my $removedNotRightLen = 0;
open IN, $tmpFile2;
while (<IN> ){
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $id = $line[0];
	my $seq = $line[1];
	if (!exists($data{$id})) {
		print STDERR "!!! Could not find $id\n";
		next;
	}
	my $s = $data{$id};
	my $mseq = $seq;
	$mseq =~ s/N//g;
	my $l1 = length($seq);
	my $l2 = length($mseq);
	my $numN = $l1-$l2;
	if ($numN < 0) {
		print STDERR "!!! There is a problem...\n";
	}
	$total++;

	if ($l1 != $lenTarget) {
		$removedNotRightLen++;
		next;
	}
	if ($l1 < 1) {
		$removed1++;
		next;
	}
	#print STDERR "$id\t$numN\t$l1\t$seq\n";
	if (($numN/$l1) > $maxN) {
		$removedN++;
		next;
	}

	print ">$id $s\n";
	print "$seq\n";
}
`rm "$tmpFile" "$tmpFile2"`;	
print STDERR "\n\tTotal Sequences: $total\n";
my $x3 = floor($removedNotRightLen/$total*100.0);
print STDERR "\tTotal removed (not right length): $removedNotRightLen ($x3 %)\n";
my $x1 = floor($removed1/$total*100.0);
print STDERR "\tTotal removed (too short): $removed1 ($x1 %)\n";
my $x2 = floor($removedN/$total*100.0);
print STDERR "\tTotal removed (too manyN): $removedN ($x2 %)\n";

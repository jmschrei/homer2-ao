#!/usr/bin/env perl
use warnings;

#
## Copyright 2009 - 2019 Christopher Benner <cbenner@ucsd.edu>
##
## This file is part of HOMER
##
## HOMER is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## HOMER is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.


sub printCMD {
	print STDERR "\n\tfilterTagDirectory.pl <tagDirectory> <peak/BED file> [options] > output.tags\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-strand (filter reads from the same strand only)\n";
	print STDERR "\t\t-tagAdjust <#> (adjust reads #bp before assigning overlap of 5' position)\n";
	print STDERR "\t\t-keep (keep tags in peaks, discard others. def: discard peak tags)\n";
	print STDERR "\t\t-update (process tags, then replace tag directory)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 2) {
	printCMD();
}
my $tagDir = $ARGV[0];
my $peakFile = $ARGV[1];
my $strandFlag = 0;
my $tagAdjust = 0;
my $keepFlag = 0;
my $updateFlag = 0;
for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-strand') {
		$strandFlag = 1;
	} elsif ($ARGV[$i] eq '-tagAdjust') {
		$tagAdjust = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-keep') {
		$keepFlag = 1;
	} elsif ($ARGV[$i] eq '-update') {
		$updateFlag = 1;
	} else {
		printCMD();
	}
}

my $rand = rand();
my $tmpFile = $rand . ".tmp";
$tmpFile2 = $rand . ".2.tmp";
`bed2pos.pl "$peakFile" -check -unique > "$tmpFile"`;
`checkPeakFile.pl "$tmpFile"`;
`cleanUpPeakFile.pl "$tmpFile" > "$tmpFile2"`;
`mv "$tmpFile2" "$tmpFile"`;

my %peaks = ();
open IN, $tmpFile;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line = split /\t/;
	my $id = $line[0];
	my $chr = $line[1];
	my $s = $line[2];
	my $e = $line[3];
	my $d = $line[4];
	$d = '+' if ($d eq '0');
	$d = '-' if ($d eq '1');
	if (!exists($peaks{$chr})) {
		my @a = ();
		$peaks{$chr} = \@a;
	}
	my $p = {id=>$id,s=>$s,e=>$e,d=>$d};
	push(@{$peaks{$chr}}, $p);
}
close IN;
`rm "$tmpFile"`;

if ($updateFlag) {
	open OUT, ">$tmpFile";
}

my $adjustBuffer = 2*abs($tagAdjust);
my $totalKept = 0;
my $total =0;
my $totalKeptV = 0;
my $totalV =0;

$ls = `ls "$tagDir"/*.tags.tsv`;
my @files = split /\n/,$ls;
foreach(@files) {
	my $tagFile = $_;
	my $chr= $tagFile;
	$chr =~ s/^.*\///;
	$chr =~ s/\.tags.tsv$//;
	print STDERR "\tProcessing file $tagFile ($chr)\n";

	my $curChr = '';
	my $curPos = -1;
	my $curIndex = 0;
	my @peaks = ();

	open IN, $tagFile;
	while (<IN>) {
		my $og = $_;
		chomp;
		s/\r//g;
		$total++;
		my @line = split /\t/;
		my $c = $line[1];
		my $p = $line[2];
		my $d = $line[3];
		my $v = $line[4];
		$totalV+=$v;
		if ($d == 0) {
			$p += $tagAdjust;
			$d = '+';
		} else {
			$p -= $tagAdjust;
			$d = '-';
		}

		if ($c ne $curChr) {
			print STDERR "\tAnalyzing $c\n";
			$curChr = $c;
			@peaks = sort {$a->{'s'} <=> $b->{'s'}} values @{$peaks{$c}};
			$curIndex = 0;
			$curPos = -1;
		}
		my $overlap = 0;
		for (my $i=$curIndex;$i<@peaks;$i++) {
			last if ($peaks[$i]->{'s'} > $p);
			if ($peaks[$i]->{'e'} < $p-$adjustBuffer && $i==$curIndex) {
				$curIndex = $i+1;
				next;
			}
			if ($peaks[$i]->{'s'} <= $p && $peaks[$i]->{'e'} >= $p) {
				if ($strandFlag) {
					if ($peaks[$i]->{'d'} eq $d) {
						$overlap=1;
						last;
					}
				} else {
					$overlap=1;
					last;
				}
			}
		}
		if ($overlap && $keepFlag) {
			$totalKept++;
			$totalKeptV+=$v;
			if ($updateFlag) {
				print OUT $og;
			} else {
				print $og;
			}
		}
		if ($overlap==0 && $keepFlag==0) {
			$totalKept++;
			$totalKeptV+=$v;
			if ($updateFlag) {
				print OUT $og;
			} else {
				print $og;
			}
		}
	}
	close IN;
}
print STDERR "\n\tPositions:\n";
print STDERR "\tTotal = $total ($totalV)\n";
print STDERR "\tKept  = $totalKept ($totalKeptV)\n";

if ($updateFlag) {
	close OUT;
	`makeTagDirectory "$tagDir" -t "$tmpFile"`;
	`rm "$tmpFile"`;
}

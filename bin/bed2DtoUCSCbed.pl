#!/usr/bin/env perl
use warnings;


# Copyright 2009 - 2018 Christopher Benner <cbenner@ucsd.edu>
#
# This file is part of HOMER
#
# HOMER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HOMER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


use POSIX;

sub printCMD {
	print STDERR "\n\tusage: bed2DtoUCSCbed.pl <2D bed file> [-peaks|-interactions|-ucsc|-ends]\n";
	print STDERR "\n\t\t-ucsc (outputs BED12 format for visualizing 2D BED file in the UCSC browser)\n";
	print STDERR "\n\t\t-ucscInteract (outputs UCSC interact format to visualize loops)\n";
	print STDERR "\t\t-peaks (outputs HOMER-style peak file spanning interaction region)\n";
	print STDERR "\t\t-interactions (outputs HOMER-style interaction file for use with annotateInteractions.pl)\n";
	print STDERR "\t\t-ends (outputs HOMER-style peak files representing the anchor points of each interaction)\n";
	print STDERR "\t\t-trackName <name> (default: filename)\n";
	print STDERR "\t\t-1 (first column is an ID, chr starts in column 2)\n";
	print STDERR "\t\t-scoreCol <#> (column [starting with 1] harboring interaction strength)\n";
	print STDERR "\t\t-cisOnly (only output cis-interactions)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 1) {
	printCMD();
}

my $trackName = $ARGV[0];
my $flag = '';
my $scoreCol = '';
my $colOneFlag = 0;
my $cisOnly = 0;
for (my $i=1;$i<@ARGV;$i++) { 
	if ($ARGV[$i] eq '-peaks') {
		$flag = 'peaks';
	} elsif ($ARGV[$i] eq '-interactions') {
		$flag = 'interactions';
	} elsif ($ARGV[$i] eq '-ucscInteract') {
		$flag = 'ucscInteract';
	} elsif ($ARGV[$i] eq '-ends') {
		$flag = 'ends';
	} elsif ($ARGV[$i] eq '-ucsc') {
		$flag = '';
	} elsif ($ARGV[$i] eq '-cisOnly') {
		$cisOnly = 1;
	} elsif ($ARGV[$i] eq '-scoreCol') {
		$scoreCol = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-trackName') {
		$trackName = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-1') {
		$colOneFlag = 1;
		print STDERR "\tWill treat 1st column as an ID\n";
	} else {
		print STDERR "!!! \"$ARGV[$i]\" not recognized!!\n";
		printCMD();
	}
}

my $score = 1;
my $color = "150,150,150";
if ($flag eq '') {
	print "track name=\"$trackName\" visibility=4\n";
} elsif ($flag eq 'ucscInteract') {
	if ($scoreCol eq '') {
		print STDERR "!!! Error - need to specify the -scoreCol when using -ucscInteract\n";
		exit;
	}

	print "track name=\"$trackName\" type=interact useScore=on maxHeightPixels=200:100:50 visibility=full\n";
} elsif ($flag eq 'interactions') {
	print "#InteractionID\tPeakID(1)\tchr(1)\tstart(1)\tend(1)\tstrand(1)\tTotal Reads(1)\tPeakID(2)\tchr(2)\tstart(2)\tend(2)\tstrand(2)\tTotal Reads(2)\tDistance\tInteraction Reads\tExpected Reads\tZ-score\tLogP\tFDR\n";
}

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	s/\"//g;
	next if (/^\s*\#/);

	my @line = split /\t/;
	if (@line < 6) {
		next;
	}
	my $id = $line[0];
	if ($colOneFlag) {
	} else {
		$id = $line[0] . ":" . $line[1] . "-" . $line[5];
	}
	my $c1 = $line[$colOneFlag+0];
	my $s1 = $line[$colOneFlag+1];
	my $e1 = $line[$colOneFlag+2];
	my $c2 = $line[$colOneFlag+3];
	my $s2 = $line[$colOneFlag+4];
	my $e2 = $line[$colOneFlag+5];
	my $s = $s1;
	$s = $s2 if ($s2 < $s1);
	my $e = $e2;
	$e = $e1 if ($e2 < $e1);
	my $dist = $e-$s;
	if ($c1 ne $c2) {
		$dist = -1;
		if ($cisOnly) {
			next;
		}
	}

	if ($flag eq 'interactions') {
		my $totalReads = 0;
		my $expectReads = 0;
		my $zscore = 0;
		my $logp = 1;
		my $fdr = 1;
		$totalReads = $line[$colOneFlag+7] if (@line > 7+$colOneFlag);
		
		print "$c1:$s-$e\t$c1:$s1-$e1\t$c1\t$s1\t$e1\t+\t1\t$c1:$s2-$e2\t$c2\t$s2\t$e2\t+\t1\t$dist";
		print "\t$totalReads\t$expectReads\t$zscore\t$logp\t$fdr\n";
		next;
	}
	if ($flag eq 'peaks') {
		next if ($c1 ne $c2);
		my $v = 1;
		$v = $line[$colOneFlag+7] if (@line > $colOneFlag+7);
		print "$c1:$s-$e\t$c1\t$s\t$e\t+\t$v";
		if (@line > $colOneFlag+8) {
			my $v2 = $line[$colOneFlag+8];
			print "\t$v2";
		}
		my $size = $e-$s;
		print "\t$size";
		print "\n";
		next;
	}
	if ($flag eq 'ends') {
		my $v = 1;
		$v = $line[$colOneFlag+7] if (@line > $colOneFlag+7);
		if ($c1 eq $c2 && $s1 eq $s2 && $e1 eq $e2) {
			my $size = $e1-$s1;
			my $halfSize = floor($size/2);
			my $ss = $s-$halfSize;
			my $ee = $s+$halfSize;
			print "$c1:$s-$e-5p\t$c1\t$ss\t$ee\t+\t$v";
			if (@line > $colOneFlag+8) {
				my $v2 = $line[$colOneFlag+8];
				print "\t$v2";
			}
			print "\t$size\n";
			$ss = $e-$halfSize;
			$ee = $e+$halfSize;
			print "$c1:$s-$e-3p\t$c1\t$ss\t$ee\t-\t$v";
			if (@line > $colOneFlag+8) {
				my $v2 = $line[$colOneFlag+8];
				print "\t$v2";
			}
			print "\t$size\n";

		} else {
			print "$c1:$s-$e-5p\t$c1\t$s1\t$e1\t+\t$v";
			if (@line > $colOneFlag+8) {
				my $v2 = $line[$colOneFlag+8];
				print "\t$v2";
			}
			my $size = $e-$s;
			print "\t$size";
			print "\n";
			print "$c1:$s-$e-3p\t$c2\t$s2\t$e2\t-\t$v";
			if (@line > $colOneFlag+8) {
				my $v2 = $line[$colOneFlag+8];
				print "\t$v2";
			}
			$size = $e-$s;
			print "\t$size";
			print "\n";
		}
		next;
	}

	if ($flag eq 'ucscInteract') {
		my $v = '.';
		if ($scoreCol > -1) {
			$v = $line[$scoreCol-1];
			$v *= 10;
		} else {
			$v = 1000;
		}
		if ($c1 eq $c2) {
			my $s = $s1;
			$s = $s2 if ($s2 < $s);
			my $e = $e2;
			$e = $e1 if ($e1 > $e2);
			#chrom  chromStart  chromEnd  name  score  value  exp  color  sourceChrom  sourceStart  sourceEnd  sourceName  sourceStrand  targetChrom  targetStart  targetEnd  targetName  targetStrand
			print "$c1\t$s\t$e\t.\t$v\t$v\t.\t0\t$c1\t$s1\t$e1\t.\t.\t$c2\t$s2\t$e2\t.\t.\n";
		}
		next;
	}

	if ($c1 eq $c2 && $s1 eq $s2 && $e1 eq $e2) {
		#TAD style
		my $dist = $e1-$s1-50;
		print "$c1\t$s1\t$e1\tTAD-$c1:$s1-$e1\t$score\t+\t$s1\t$e1\t$color\t2\t50,50\t0,$dist\n";
		next;
	} 
	if ($c1 ne $c2) {
		#interchromosomal interaction
		my $dist = $e1-$s1-100;
		print "$c1\t$s1\t$e1\tInterchromosomal-Interaction->$c2:$s2-$e2\t$score\t+\t$s1\t$e1\t$color\t1\t$dist,\t0,\n";
		$dist = $e2-$s2-100;
		print "$c2\t$s2\t$e2\tInterchromosomal-Interaction->$c1:$s1-$e1\t$score\t+\t$s2\t$e2\t$color\t1\t$dist,\t0,\n";
		next;
	}
	if ($s1 < $s2) {
		my $dist1 = $e1-$s1;
		my $dist2 = $e2-$s2;
		my $dist = $e2-$s1-$dist2;
		print "$c1\t$s1\t$e2\tInteraction-$c1:$s1-$e2\t$score\t+\t$s1\t$e2\t$color\t2\t$dist1,$dist2,\t0,$dist,\n";
	} else {
		my $dist1 = $e1-$s1;
		my $dist2 = $e2-$s2;
		my $dist = $e1-$s2-$dist1;
		print "$c1\t$s2\t$e1\tInteraction-$c1:$s2-$e1\t$score\t+\t$s2\t$e1\t$color\t2\t$dist2,$dist1,\t0,$dist,\n";
	}

}
close IN;

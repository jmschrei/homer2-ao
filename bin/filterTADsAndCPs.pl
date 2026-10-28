#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";


# Copyright 2009 - 2017 Christopher Benner <cbenner@ucsd.edu>
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
#use HomerConfig;
#use Statistics;

my $prefix = "out";
my $bed2Dtad = '';
my $bed2Dint = '';
my $genome = "none";
my $res = 15000;
my $badRegionsFile = "";
my $tadCovFactor = 0.25;
my $intCovFactor = 0.25;
my $minIntScore = 0;
my $maxCPUs = 1;
my $minTADir = 0;
my %badChr = ();
#$badChr{'chrX'} = 1;
#$badChr{'chrY'} = 1;
#$badChr{'chrM'} = 1;

my $tadScoreFile = '';
my $intScoreFile = '';

my $tadFilter = 1;
my @tagDirs = ();

sub printCMD {
	print STDERR "\n\tfilterTADsAndInts.pl [options] -cp <CP 2D bed> -tad <TAD 2D bed> -d <tagDir> -genome <genome> -o <outPrefix>\n";
	print STDERR "\n\tParamters: (no strict order needed)\n";
	print STDERR "\t\t-o <prefix for output files> (def: $prefix)\n";
	print STDERR "\t\t-cp <contact point 2D bed> (input contact point file, 2D BED format\n";
	print STDERR "\t\t-tad <TAD 2D bed> (input TAD file, 2D BED format\n";
	print STDERR "\t\t-d <Hi-C tag dirctory> [tagdir2] ... (Used to calculate experiment coverage)\n";
	print STDERR "\t\t\t-cpu <#> (number of processes to use to count reads in mult. tagdirs at once)]n";
	print STDERR "\t\t-genome <genome version> (genome version, def: none)\n";
	print STDERR "\t**\t-p <peak/BED file of duplication/bad/blacklist regions> (def: gaps/segDup from UCSC)\n";
	print STDERR "\t\t-badChr <chr> [chr] ... (chromosomes to skip e.g. chrY, def: none)\n";
	print STDERR "\n\tThresholds:\n";
	print STDERR "\t\t-res <#> (resolution for interaction/region overlap, def: $res)\n";
	print STDERR "\t\t-minTADir <#> (minimum TAD inclusion ratio, def: $minTADir)\n";
	print STDERR "\t\t-minTADcov <#> (minimum TAD coverage ratio relative to average, def: $tadCovFactor)\n";
	print STDERR "\t\t-minCPscore <#> (minimum int score, def: $minIntScore)\n";
	print STDERR "\t\t-minCPcov <#> (minimum int coverage ratio relative to average, def: $intCovFactor)\n";
	print STDERR "\n\tco-filter score files:\n";
	print STDERR "\t\t-TADscore <TAD score file>\n";
	print STDERR "\t\t-CPscore <CP score file>\n";
	print STDERR "\n\n\t** Notes on -p <peak/BED bad regions>:\n";
	print STDERR "\t\tRecommend to use genome assembly gaps and known duplicated regions. To get these if you're using a UCSC\n";
	print STDERR "\t\tannotated genome, try something like the following to create BED file with these regions (for hg38):\n";
	print STDERR "\t\t\twget -O gap.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz\n";
	print STDERR "\t\t\twget -O dups.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz\n";
	print STDERR "\t\t\tzcat gap.txt.gz dups.txt.gz | cut -f2-4 > badRegions.bed\n";
	print STDERR "\t\t\t(then use \"badRegions.bed\" with the -p option\n";
	print STDERR "\t\tNote that depending on the genome, UCSC might have different file names...\n";
	print STDERR "\n";
	exit;
}
if (@ARGV < 2) {
	printCMD();
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-cp' || $ARGV[$i] eq '-ints') {
		$bed2Dint = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-tad' || $ARGV[$i] eq '-tads') {
		$bed2Dtad = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-prefix' || $ARGV[$i] eq '-o') {
		$prefix = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-res' || $ARGV[$i] eq '-superRes' || $ARGV[$i] eq '-window') {
		$res = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-p') {
		$badRegionsFile = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-minTADir') {
		$minTADir = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-minTADcov') {
		$tadCovFactor = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-minCPscore') {
		$minIntScore = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-minCPcov') {
		$intCovFactor = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-CPscore') {
		$intScoreFile = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-TADscore') {
		$tadScoreFile = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
 	} elsif ($ARGV[$i] eq '-d') {
		$i++;
		while ($i<@ARGV) {
			if ($ARGV[$i] =~ /^\-/) {
				$i--;
				last;
			}
			push(@tagDirs,$ARGV[$i]);
			$i++;
		}
 	} elsif ($ARGV[$i] eq '-badChr') {
		%badChr = ();
		$i++;
		while ($i<@ARGV) {
			if ($ARGV[$i] =~ /^\-/) {
				$i--;
				last;
			}
			$badChr{$ARGV[$i]}=1;
			$i++;
		}
	} else {
		print STDERR "!!! What is \"$ARGV[$i]\"?\n";
		printCMD();
	}
}


if ($bed2Dtad eq '' || $bed2Dint eq '') {
	print STDERR "TAD and int files are required\n";
	exit;
}


my $rand = rand();
#$rand = "0.264615196109848";
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $halfRes = floor($res/2);
my $tads = read2Dbed($bed2Dtad,"tad");
my $ints = read2Dbed($bed2Dint,"int");

open OUT, ">$tmpFile";
foreach(keys %$tads) {
	my $id = $_;
	my $tad = $tads->{$_};
	my $s = $tad->{'s'}-floor($tad->{'l'}/2);
	my $e = $tad->{'s'};
	print OUT "$id-upstream\t$tad->{'c'}\t$s\t$e\t+\n";
	$s = $tad->{'s'};
	$e = $tad->{'e'};
	print OUT "$id-tad\t$tad->{'c'}\t$s\t$e\t+\n";
	$s = $tad->{'e'};
	$e = $tad->{'e'}+floor($tad->{'l'}/2);
	print OUT "$id-downstream\t$tad->{'c'}\t$s\t$e\t+\n";
	$s = $tad->{'s'};
	$e = $tad->{'e'};
	my $ss = $s-$halfRes;
	my $ee = $s+$halfRes;
	print OUT "$id-upBoundary\t$tad->{'c'}\t$ss\t$ee\t+\n";
	$ss = $e-$halfRes;
	$ee = $e+$halfRes;
	print OUT "$id-downBoundary\t$tad->{'c'}\t$ss\t$ee\t-\n";
}

foreach(keys %$ints) {
	my $id = $_;
	my $int = $ints->{$id};
	my $s = $int->{'s'};
	my $e = $int->{'e'};
	my $ss = $s-$halfRes;
	my $ee = $s+$halfRes;
	print OUT "$id-upAnchor\t$int->{'c'}\t$ss\t$ee\t+\n";
	$ss = $e-$halfRes;
	$ee = $e+$halfRes;
	print OUT "$id-downAnchor\t$int->{'c'}\t$ss\t$ee\t-\n";
}

close OUT;

my $tagDirStr = '';
my $numTagDirs = scalar(@tagDirs);
foreach(@tagDirs) {
	$tagDirStr .= " \"$_\" ";
}
if ($tagDirStr ne '') {
	$tagDirStr = " -d " . $tagDirStr;
}
my $badRegionsStr = "";
if ($badRegionsFile ne '') {
	$badRegionsStr = " -p \"$badRegionsFile\" ";
}
	

`annotatePeaks.pl $tmpFile $genome -pcount $badRegionsStr $tagDirStr -fpkm -cpu $maxCPUs > $tmpFile2`;


my %badRegionOverlap = ();
my @tadAvgCoverage = ();
my @tadNCoverage = ();
my @intAvgCoverage = ();
my @intNCoverage = ();
for (my $i=0;$i<$numTagDirs;$i++) {
	push(@tadAvgCoverage,0);
	push(@tadNCoverage,0);
	push(@intAvgCoverage,0);
	push(@intNCoverage,0);
}

my $covColID = -1;
my $badColID = -1;
if ($numTagDirs > 0) {
	$covColID = 19;
	if ($badRegionsFile ne '') {
		$badColID = $covColID + $numTagDirs;
	}
} elsif ($badRegionsFile ne '') {
	$badColID = 19;
}




open IN, $tmpFile2;

while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $type = 0;
	my $id = '';
	if ($line[0] =~ /^(.+)-upstream/) {
		$type = 'tad-up-fpkm';
		$id = $1;
 	} elsif ($line[0] =~ /^(.+)-downstream/) {
		$type = 'tad-down-fpkm';
		$id = $1;
 	} elsif ($line[0] =~ /^(.+)-tad/) {
		$type = 'tad-tad-fpkm';
		$id = $1;
 	} elsif ($line[0] =~ /^(.+)-upBoundary/) {
		$type = 'tad-boundary';
		$id = $1;
 	} elsif ($line[0] =~ /^(.+)-downBoundary/) {
		$type = 'tad-boundary';
		$id = $1;
 	} elsif ($line[0] =~ /^(.+)-upAnchor/) {
		$type = 'int-up-fpkm';
		$id = $1;
 	} elsif ($line[0] =~ /^(.+)-downAnchor/) {
		$type = 'int-dn-fpkm';
		$id = $1;
	} 
	next if ($id eq '' || $type eq '');

	if ($badColID > 0 && @line > $badColID && $line[$badColID] > 0) {
		if ($type =~ /^int/) {
			$ints->{$id}->{'bad'} = 1;
		} elsif ($type eq 'tad-boundary') {
			$tads->{$id}->{'bad'} = 1;
		}
	}
	next if ($type eq 'tad-boundary');
	for (my $i=0;$i<$numTagDirs;$i++) {
		my $coverage = $line[$covColID+$i];

		if ($type =~ /^tad/) {
			next if (!exists($tads->{$id}));
			if (!exists($tads->{$id}->{$type})) {
				my @a = ();
				$tads->{$id}->{$type} = \@a;
			}
			$tads->{$id}->{$type}->[$i] = $coverage;
			if ($type eq 'tad-tad-fpkm') {
				$tadAvgCoverage[$i] += $coverage;
				$tadNCoverage[$i] += 1.0;
			}
		} elsif ($type =~ /^int/) {
			next if (!exists($ints->{$id}));
			if (!exists($ints->{$id}->{$type})) {
				my @a = ();
				$ints->{$id}->{$type} = \@a;
			}
			$ints->{$id}->{$type}->[$i] = $coverage;
			$intAvgCoverage[$i] += $coverage;
			$intNCoverage[$i] +=1.0;
		}
	}
}
close IN;

for (my $i=0;$i<$numTagDirs;$i++) {
	if ($tadNCoverage[$i] > 0) {
		$tadAvgCoverage[$i] /= $tadNCoverage[$i];
	}
	if ($intNCoverage[$i] > 0) {
		$intAvgCoverage[$i] /= $intNCoverage[$i];
	}
	print STDERR "\t[$tagDirs[$i]] Average TAD fpkm coverage = $tadAvgCoverage[$i] ($tadNCoverage[$i] regions total)\n";
	print STDERR "\t[$tagDirs[$i]] Average CP anchor fpkm coverage = $intAvgCoverage[$i] ($intNCoverage[$i] regions total)\n";
}

#`rm $tmpFile $tmpFile2`;

open OUT, ">$prefix.$bed2Dtad";

my %scores = ();

if ($tadScoreFile ne '') {
	open OUT2, ">$prefix.$tadScoreFile";
	open IN, $tadScoreFile;
	my $c = 0;
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my $og = $_;
		if ($c==1) {
			print OUT2 "$og\n";
			next;
		}
		my @line = split /\t/;
		my $x = $line[0];
		if ($x =~ /(.*?)\:(\d+)\-(\d+)/) {
			$x = $1 . ":" . ($2-1) . "-" . $3;
		}
		my $id = "tad-" . $x;
		$scores{$id} = $og;
	}
	close IN;
}


my $total = 0;	
my $numBadChr = 0;
my $numBadRegion = 0;
my $numBadCoverage = 0;
my $numBadIR = 0;
my $numGood=0;
foreach(keys %$tads) {
	$total++;
	my $id= $_;
	my $c = $tads->{$id}->{'c'};
	if ($tads->{$id}->{'ir'} ne 'NA') {
		if ($tads->{$id}->{'ir'} < $minTADir) {
			$numBadIR++;
			next;
		}
	}
	if (exists($badChr{$c})) {
		$numBadChr++;
		next;
	}
	if ($tads->{$id}->{'bad'}) {
		$numBadRegion++;
		next;
	}
	my $bcov = 0;
	for (my $i=0;$i<$numTagDirs;$i++) {
		if (exists($tads->{$id}->{'tad-tad-fpkm'}) && exists($tads->{$id}->{'tad-down-fpkm'}) 
										&& exists($tads->{$id}->{'tad-up-fpkm'})) {
			my $tc = $tads->{$id}->{'tad-tad-fpkm'}->[$i];
			my $uc = $tads->{$id}->{'tad-up-fpkm'}->[$i];
			my $dc = $tads->{$id}->{'tad-down-fpkm'}->[$i];
			if ($tc/$tadAvgCoverage[$i] < $tadCovFactor) {
				$bcov =1;
				last;
			}
			if ($uc/$tc < $tadCovFactor) {
				$bcov =1;
				last;
			}
			if ($dc/$tc < $tadCovFactor) {
				$bcov =1;
				last;
			}
		} else {
			$bcov =1;
			last;
		}
	}
	if ($bcov == 1) {
		$numBadCoverage++;
		next;
	}
	my $og = $tads->{$id}->{'og'};
	print OUT "$og\n";
	if ($tadScoreFile ne '') {
		if (exists($scores{$id})) {
			my $og = $scores{$id};
			print OUT2 "$og\n";
		} else {
				#print STDERR "Can't find $id\n";
		}
	}
	$numGood++;
}
close OUT;
if ($tadScoreFile ne '') {
	close OUT2;
}

print STDERR "\n\tTAD Stats:\n";
my $ratio = sprintf("%.3lf", 100*$numGood/$total);
print STDERR "\t\tGood = $numGood/$total ($ratio%)\n";
$ratio = sprintf("%.3lf", 100*$numBadIR/$total);
print STDERR "\t\tBad IR = $numBadIR/$total ($ratio%)\n";
$ratio = sprintf("%.3lf", 100*$numBadChr/$total);
print STDERR "\t\tBad chr = $numBadChr/$total ($ratio%)\n";
$ratio = sprintf("%.3lf", 100*$numBadRegion/$total);
print STDERR "\t\tBad region = $numBadRegion/$total ($ratio%) *not enforced for TADs\n";
$ratio = sprintf("%.3lf", 100*$numBadCoverage/$total);
print STDERR "\t\tBad coverage = $numBadCoverage/$total ($ratio%)\n";


open OUT, ">$prefix.$bed2Dint";

%scores = ();
if ($intScoreFile ne '') {
	open OUT2, ">$prefix.$intScoreFile";
	open IN, $intScoreFile;
	my $c = 0;
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my $og = $_;
		if ($c==1) {
			print OUT2 "$og\n";
			next;
		}
		my @line = split /\t/;
		my $id = "int-" . $line[0];
		$scores{$id} = $og;
	}
	close IN;
}


$total = 0;	
$numBadChr = 0;
$numBadRegion = 0;
$numBadCoverage = 0;
$numBadScore = 0;
$numGood=0;
foreach(keys %$ints) {
	$total++;
	my $id= $_;
	my $c = $ints->{$id}->{'c'};
	if ($ints->{$id}->{'ir'} ne 'NA') {
		if ($ints->{$id}->{'ir'} < $minIntScore) {
			$numBadScore++;
			next;
		}
	}
	if (exists($badChr{$c})) {
		$numBadChr++;
		next;
	}
	if ($ints->{$id}->{'bad'}) {
		$numBadRegion++;
		next;
	}
	my $bcov = 0;
	for (my $i=0;$i<$numTagDirs;$i++) {
		if (exists($ints->{$id}->{'int-up-fpkm'}) && exists($ints->{$id}->{'int-dn-fpkm'})) {
			my $up = $ints->{$id}->{'int-up-fpkm'}->[$i];
			my $dn = $ints->{$id}->{'int-dn-fpkm'}->[$i];
			if ($up/$intAvgCoverage[$i] < $intCovFactor) {
				$bcov=1;
				last;
			}
			if ($dn/$intAvgCoverage[$i] < $intCovFactor) {
				$bcov=1;
				last;
			}
		} else {
			$bcov=1;
			last;
		}
	}
	if ($bcov == 1) {
		$numBadCoverage++;
		next;
	}
	my $og = $ints->{$id}->{'og'};
	print OUT "$og\n";
	if ($intScoreFile ne '') {
		if (exists($scores{$id})) {
			my $og = $scores{$id};
			print OUT2 "$og\n";
		}
	}
	$numGood++;
}
close OUT;
if ($intScoreFile ne '') {
	close OUT2;
}

print STDERR "\n\tCP Stats:\n";
$ratio = sprintf("%.3lf", 100*$numGood/$total);
print STDERR "\t\tGood = $numGood/$total ($ratio%)\n";
$ratio = sprintf("%.3lf", 100*$numBadScore/$total);
print STDERR "\t\tBad IR = $numBadScore/$total ($ratio%)\n";
$ratio = sprintf("%.3lf", 100*$numBadChr/$total);
print STDERR "\t\tBad chr = $numBadChr/$total ($ratio%)\n";
$ratio = sprintf("%.3lf", 100*$numBadRegion/$total);
print STDERR "\t\tBad region = $numBadRegion/$total ($ratio%)\n";
$ratio = sprintf("%.3lf", 100*$numBadCoverage/$total);
print STDERR "\t\tBad coverage = $numBadCoverage/$total ($ratio%)\n";
print STDERR "\n";

exit;




sub read2Dbed {
	my ($file,$type) = @_;
	open IN, $file or die "Could not open $file\n";
	my @otherdata = ();
	my %regions = ();
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^\#/);
		my $og = $_;
		my @line = split /\t/;
		my $chr = $line[0];
		my $start = $line[1];
		my $end = $line[5];
		my $id = "$type-$chr:$start-$end";
		my $ir = "NA";
		for (my $i=6;$i<@line;$i++) {
			push(@otherdata, $line[$i]);
			if ($i==7) {
				$ir = $line[$i];
			}
		}
		$regions{$id} = {s=>$start,e=>$end,og=>$og,c=>$chr,l=>$end-$start,bad=>0,ir=>$ir};
	}
	close IN;
	return \%regions;
}


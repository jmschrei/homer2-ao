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
# GNU General Public License for more details./bin/bash: use: command not found


use POSIX;

my $minDist = 25000;
my $maxDist = 2000000;
my $insDist = 150000;
my $diDist = 1000000;
my $minTADsize = 45000;
my $minTADscore = 1.75;
my $sizeTADoptima = '';
my $findTADoptions = "";
my $noOverlapOption = " -noOverlap ";
my $res = 3000;
my $superRes = 15000;
my $minIntReads = 10;
my $stdFilter = 8;
my $minFilter = 0.05;
my $prefix = 'out';
my $sizeCP = '';
my $sizeCPlocal = '';
my $minCPdist = '';
my $minCPreads = '';
my $foldCPlocal = 1.5;
my $foldCPglobal = 2.0;
my $poissonCP = 1e-5;
my $rawFlag = 0;
my $normTotal = '';
my $balanceFlag = " ";
#my $balanceFlag = " -balance ";

my $findTADexe = "findTADsFromRelativeMatrix";
my $analyzeHiCexe = "analyzeHiC";

my $maxCPUs = 1;

print STDERR "\tPlease use findTADsAndLoops.pl instead...\n";
exit;

sub printCMD {
	print STDERR "\n\tfindTADsAndCPs.pl <cmd> [options]\n";	
	print STDERR "\n\tUsage Modes:\n";
	print STDERR "\n\t\tfindTADsAndCPs.pl find <tagDir> [options]\n";
	print STDERR "\t\t\tFinds TADs and CPs (contact points/anchors) in the given Hi-C tag directory. Output files\n";
	print STDERR "\t\t\twill have the same prefix, and include 2D bed files specifying the locations of the features.\n";
	print STDERR "\n\t\tfindTADsAndCPs.pl score -tad <tad.2D.bed> -cp <cp.2D.bed> -d <tagDir1> [tagDir2] ... [options]\n";
	print STDERR "\t\t\tScores TADs and CPs cross several Hi-C tag directories so that they can be compared across\n";
	print STDERR "\t\t\tdifferent replicates and experimental conditions.\n";
	print STDERR "\n\tCommon Options:\n";
	print STDERR "\t\t-prefix <output prefix> (default: $prefix)\n";
	print STDERR "\t\t-cpu <#> (max CPUs to use)\n";
	print STDERR "\t\t-res <#> (resolution, default: $res)\n";
	print STDERR "\t\t-window <#> (overlapping window size, i.e. superRes, default: $superRes)\n";
	#print STDERR "\t\t-balance | -nobalance (analyzed balanced contact matrices, default: uses balanced matrices)\n";
	print STDERR "\t\t-minDist <#> (minimum interaction distance to consider, default: $minDist)\n";
	print STDERR "\t\t-maxDist <#> (maximum interaction distance to consider, default: $maxDist)\n";

	print STDERR "\n\tOptions specific to 'find' mode:\n";
	print STDERR "\t\tTAD options:\n";
	print STDERR "\t\t  -minTADsize <#> (default: $minTADsize)\n";
	print STDERR "\t\t  -minTADscore <#> (minimum inclusion score to call TAD, default: $minTADscore)\n";
	print STDERR "\t\t  -keepOverlappingTADs (by default, program removes overlapping TAD calls)\n";
	print STDERR "\t\tCP options:\n";
	print STDERR "\t\t  -minCPdist <#> (minimum length of contact points to find, default: 3x -window/-superRes size)\n";
	print STDERR "\t\t  -sizeCP <#> (size of region used to calculate CP score, default: -window/-superRes size)\n";
	print STDERR "\t\t  -sizeCPlocal <#> (size of region used to calculate local bg, default: 3x -window/-superRes size)\n";
	print STDERR "\t\t  -foldCPlocal <#> (fold threshold for local CP enrichment, default: $foldCPlocal)\n";
	print STDERR "\t\t  -foldCPglobal <#> (fold threshold for global CP vs. interactions of same distance, default: $foldCPglobal)\n";
	print STDERR "\t\t  -poissonCP <#> (poisson p-value threshold for minimum interactions per CP, default: 1e-5)\n";
	print STDERR "\t\t  -minCPreads <#> (minimum number of adjusted reads for contact point, default: use -poissonCP)\n";
	#print STDERR "\t\t  -estFDR <#> (set CP threshold based on estimated FDR, default: 0.01)\n";
	print STDERR "\t\tOther options:\n";
	print STDERR "\t\t  -insDist <#> (distance across which to calculate insulation ratio, default: $insDist)\n";
	print STDERR "\t\t  -diDist <#> (distance across which to calculate directionality index , default: $diDist)\n";
	print STDERR "\n\tOptions specific to 'score' mode:\n";
	print STDERR "\t\t-d <HiC TagDir1> [HiC TagDir2] ...\n";
	print STDERR "\t\t-cp <contact point 2D BED file> (will score CP interaction rates from this file)\n";
	print STDERR "\t\t-tad <tad 2D BED file> (will score TAD inclusion ratios from this file)\n";
	print STDERR "\t\tCP score reporting: (default: Hi-C depth normalized counts [use for use with getDiffExpression.pl...])\n";
	print STDERR "\t\t  -raw (report raw interaction counts, will impact TAD inclusion ratio scoring too)\n";
	print STDERR "\t\t  -normTotal <#> (normalize based on this depth of sequencing for Hi-C experiments)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 2) {
	printCMD();
}

my $tagDir='';
my $tadfile = '';
my @tagDirs = ();
my $tagDirStr = '';
my $int2dBed = '';
my $tad2dBed = '';
my $startIndex = 1;
my $mode = $ARGV[0];
if ($mode eq 'find') {
	$tagDir = $ARGV[1];
	push(@tagDirs, $tagDir);
	$startIndex++;
} elsif ($mode eq 'score') {
} else {
	print STDERR "!!!Error - first argument must be either \"find\" or \"score\"!!!\n";
	printCMD();
}
for (my $i=$startIndex;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-res') {
		$res = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-superRes' || $ARGV[$i] eq '-window' || $ARGV[$i] eq '-win') {
		$superRes = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-prefix') {
		$prefix = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minTADsize') {
		$minTADsize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minTADscore') {
		$minTADscore = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-keepOverlappingTADs') {
		$noOverlapOption = "";
	} elsif ($ARGV[$i] eq '-diDist') {
		$diDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-balance') {
		$balanceFlag = " -balance ";
	} elsif ($ARGV[$i] eq '-nobalance') {
		$balanceFlag = " ";
	} elsif ($ARGV[$i] eq '-insDist') {
		$insDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxDist') {
		$maxDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minDist') {
		$minDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minCPreads') {
		$minCPreads = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-poissonCP') {
		$poissonCP = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minCPdist') {
		$minCPdist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cp') {
		$int2dBed = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-tad') {
		$tad2dBed = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-sizeCP') {
		$sizeCP = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-sizeCPlocal') {
		$sizeCPlocal = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-foldCPlocal') {
		$foldCPlocal = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-foldCPglobal') {
	  	$foldCPglobal = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-sizeTADoptima') {
	  	$sizeTADoptima = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-raw') {
		$rawFlag = 1;
	} elsif ($ARGV[$i] eq '-normTotal') {
		$normTotal = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-d') {
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@tagDirs, $ARGV[$i]);
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail);
		$i--;
	} else {
		print STDERR "!!! Error: \"$ARGV[$i]\" not recognized!!!\n";
		printCMD();
	}
}

if ($tagDir eq '') {
	if (scalar(@tagDirs) > 0) {
		$tagDir = $tagDirs[0];
	}
	if ($tagDir eq '') {
		print STDERR "Somethings wrong - no tag directories specified!!\n";
		exit;
	}
}
my $chrfiles = `ls -1 "$tagDir"/*.tags.tsv`;
my @chrs = split /\n/,$chrfiles;
foreach(@chrs) {
	$_ =~ s/^.*\///;
	$_ =~ s/\.tags\.tsv//;
	#print STDERR "$_\n";
}

my $rand = rand();

my @resultFiles =();

$sizeCP = $superRes if ($sizeCP eq '');
$sizeCPlocal = $superRes*3 if ($sizeCPlocal eq '');
$minCPdist = $superRes*3 if ($minCPdist eq '');
$sizeTADoptima = $superRes*3 if ($sizeTADoptima eq '');


if ($poissonCP eq '' && $mode eq 'find') {
	my $line = `egrep "^gsizeEstimate" $tagDir/tagInfo.txt`;
	if (my @line = split /\=/,$line) {
		my $size = $line[1];
		my $estTests = ($size / $sizeCPlocal);
		#my $estTests = ($size / $sizeCPlocal) * (($maxDist-$minDist)/$sizeCPlocal);
		$poissonCP = $estFDR / $estTests;
		print STDERR "\tEstimated poisson p-value threshold = $poissonCP (FDR=$estFDR)\n";
	} else {
		print STDERR "\t!!Error, couldn't determine genome size from tagdir ($tagDir)\n";
		print STDERR "\tSetting -poissonCP to 1e-7\n";
		$poissonCP = 1e-5;
	}
} 


#$minCPreads = 10 if ($minCPreads eq '');

print STDERR "\tWill analyze chrs: @chrs\n";

#for testing
#@chrs = ("chr17");

my $cpus = 0;
#first make all of the matricies
print STDERR "\tCreating relative matrix files...\n";
foreach(@chrs) {
	my $chr = $_;
	my $prefix = "$rand.$chr";

	for (my $z=0;$z<@tagDirs;$z++) {
		my $td = $tagDirs[$z];

		my $rmatrixFile = "$rand.$chr.$z.rmatrix.txt";

		my $pid = fork();
		$cpus++;
		if ($pid == 0) {
			#child process
			my $opt = " -res $res ";
			$opt .= " -superRes $superRes ";
			$opt .= " -maxDist $maxDist ";
			$opt .= " -chr $chr ";
			$opt .= " -std $stdFilter ";
			$opt .= " -min $minFilter ";
			$opt .= " -relative ";
			if ($rawFlag) {
				$opt .= " -raw ";
			} else {
				$opt .= " -simpleNorm ";
				$opt .= " -minExpect 1.0 ";
			}
			$opt .= " -normTotal $normTotal " if ($normTotal  ne '');
			$opt .= $balanceFlag;

			`$analyzeHiCexe "$td" $opt > "$rmatrixFile"`;
			exit(0);
		}
		if ($cpus >= $maxCPUs) {
			my $id = wait();
			$cpus--;
		}
	}
}
my $id = 0;
while ($id >=0) {
	$id = wait();
}

print STDERR "\tProcessing matricies into TADs/CPs...\n";
$cpus = 0;
foreach(@chrs) {
	my $chr = $_;
	my $prefix = "$rand.$chr";
	push(@resultFiles, $prefix);
	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		#child process
		if ($mode eq 'find') {

			my $rmatrix = "$rand.$chr.0.rmatrix.txt";

			my $opt = $findTADoptions;
			$opt .= " -prefix $prefix ";
			$opt .= " -minDist $minDist ";
			$opt .= " -diDist $diDist ";
			$opt .= " -insDist $insDist ";

			$opt .= " -minTADscore $minTADscore ";
			$opt .= " -minTADsize $minTADsize ";
			$opt .= " -sizeTADoptima $sizeTADoptima ";
			$opt .= $noOverlapOption;

			$opt .= " -sizeCP $sizeCP ";
			$opt .= " -sizeCPlocal $sizeCPlocal ";
			$opt .= " -minCPreads $minCPreads " if ($minCPreads ne '');
			$opt .= " -minCPdist $minCPdist ";
			$opt .= " -foldCPlocal $foldCPlocal ";
			$opt .= " -foldCPglobal $foldCPglobal ";
			$opt .= " -poissonCP $poissonCP ";

			print STDERR "`$findTADexe find $rmatrix $opt`;\n";
			`$findTADexe find $rmatrix $opt`;

			#`$findTADexe find $rmatrix -minTADsize $minTADsize -prefix $prefix $findTADoptions -sizeCP $sizeCP -sizeCPlocal $sizeCPlocal -minCPreads $minCPreads -diDist $diDist -insDist $insDist -minCPdist $minCPdist -minDist $minDist -minTADscore $minTADscore -foldCPlocal $foldCPlocal -foldCPglobal $foldCPglobal -sizeTADoptima $sizeTADoptima`;
			`rm "$rmatrix"`;

		} elsif ($mode eq 'score') {

			my $tagDirStr = "";
			my $z = 0;
			for (my $z=0;$z<@tagDirs;$z++) {
				my $td = $tagDirs[$z];
				$tagDirStr .= " \"$rand.$chr.$z.rmatrix.txt\"";
			}
			my $tadfileChr = $rand . ".tadfile.$chr";
			my $intfileChr = $rand . ".cpfile.$chr";

			my $opt = '';
			if ($int2dBed ne '') {
				my $match = $chr . '\t\d+\t\d+\t' . $chr . '\t';
				`grep -P "$match" $int2dBed > $intfileChr`;
				$opt .= " -cp $intfileChr";
			}
			if ($tad2dBed ne '') {
				my $match = $chr . '\t\d+\t\d+\t' . $chr . '\t';
				`grep -P "$match" $tad2dBed > $tadfileChr`;
				$opt .= " -tad $tadfileChr";
			}
			$opt .= " -prefix $prefix ";
			$opt .= " -matrix $tagDirStr ";
			$opt .= " -minDist $minDist ";
			$opt .= " -sizeCP $sizeCP ";
			$opt .= " -sizeCPlocal $sizeCPlocal ";

			`$findTADexe score $opt`;

			`rm $tagDirStr`;
			`rm $intfileChr` if ($int2dBed ne '');
			`rm $tadfileChr` if ($tad2dBed ne '');
		}
		exit(0);
	}
	if ($cpus >= $maxCPUs) {
		my $id = wait();
		$cpus--;
	}
}
$id = 0;
while ($id >=0) {
	$id = wait();
}


if ($mode eq 'score') {
	open OUT, ">$prefix.tad.scores.txt";
	print OUT "#TAD name\tchr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tInclustionRatio(IR)\tLength";
	foreach(@tagDirs) {
		print OUT "\t$_ IR";
	}
	print OUT "\n";
	for (my $i=0;$i<@resultFiles;$i++) {
		my $scoreFile = "$resultFiles[$i].tad.scores.txt";
		open IN, "$scoreFile" or print STDERR "\t!!! Can't open $scoreFile";
		my $c = 0;
		while (<IN>) {
			$c++;
			next if ($c < 2);
			print OUT $_;
		}
		close IN;
		`rm $scoreFile`;
	}
	close OUT;

	if ($int2dBed ne '') {
		open OUT, ">$prefix.cp.scores.txt";
		print OUT "#CP name\tchr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tScore\tlocalRatio";
		foreach(@tagDirs) {
			print OUT "\t$_ Score";
		}
		print OUT "\n";
		for (my $i=0;$i<@resultFiles;$i++) {
			my $scoreFile = "$resultFiles[$i].cp.scores.txt";
			open IN, "$scoreFile" or print STDERR "\t!!! Can't open $scoreFile";
			my $c = 0;
			while (<IN>) {
				$c++;
				next if ($c < 2);
				print OUT $_;
			}
			close IN;
			`rm $scoreFile`;
		}
		close OUT;
	}
	exit;
}
open BED2D, ">$prefix.tad.2D.bed";
print BED2D "#TAD chr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tInclusionRatio\tInclusionRatio\n";
open BED2DINT, ">$prefix.cp.2D.bed";
print BED2DINT "#CP chr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tCPscore\tLocalFold\n";
open BEDGRAPH, ">$prefix.DI.bedGraph";
print BEDGRAPH "track name=\"$prefix DI ($tagDir)\" type=bedGraph visibility=2\n";
open INSBEDGRAPH, ">$prefix.Insulation.bedGraph";
print INSBEDGRAPH "track name=\"$prefix Insulation Ratio ($tagDir)\" type=bedGraph visibility=2\n";
#open TAD, ">$prefix.tad.txt";
#print TAD "name\tchr\tstart\tend\tdirection\tInclusionRatio\tLength\n";
my $c = 0;
my $totalTADs=0;
my $totalInts=0;

for (my $z=0;$z<@resultFiles;$z++) {
	my $f = $resultFiles[$z];

	open IN, "$f.DI.bedGraph";
	$c=0;
	while (<IN>) {
		$c++;
		next if ($c==1);
		print BEDGRAPH "$_";
	}
	close IN;
	`rm $f.DI.bedGraph`;

	open IN, "$f.Insulation.bedGraph";
	$c=0;
	while (<IN>) {
		$c++;
		next if ($c==1);
		print INSBEDGRAPH "$_";
	}
	close IN;
	`rm $f.Insulation.bedGraph`;

	open IN, "$f.tad.2D.bed";
	$c=0;
	while (<IN>) {
		$c++;
		#next if ($c==1);
		print BED2D "$_";
		$totalTADs++;
	}
	close IN;
	`rm $f.tad.2D.bed`;

	open IN, "$f.cp.2D.bed";
	$c=0;
	while (<IN>) {
		$c++;
		#next if ($c==1);
		print BED2DINT "$_";
		$totalInts++;
	}
	close IN;
	`rm $f.cp.2D.bed`;
}
close BED2D;
close BED2DINT;
#close TAD;
close BEDGRAPH;
close INSBEDGRAPH;
`bed2DtoUCSCbed.pl "$prefix.tad.2D.bed" > "$prefix.ucsc.bed"`;
`bed2DtoUCSCbed.pl "$prefix.cp.2D.bed" >> "$prefix.ucsc.bed"`;
`cat "$prefix.DI.bedGraph" "$prefix.Insulation.bedGraph" >> "$prefix.ucsc.bed"`;
`bed2DtoUCSCbed.pl "$prefix.tad.2D.bed" -peaks > "$prefix.tad.peaks.txt"`;
`bed2DtoUCSCbed.pl "$prefix.cp.2D.bed" -interactions > "$prefix.interactions.txt"`;
print STDERR "\n\tTotal TADs: $totalTADs\n";
print STDERR "\n\tTotal CPs : $totalInts\n";
print STDERR "\n";

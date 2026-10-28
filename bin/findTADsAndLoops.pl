#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";


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
# GNU General Public License for more details./bin/bash: use: command not found


use POSIX;

my $minDist = 25000;
my $maxDist = 2000000;
my $insDist = 150000;
my $diDist = 1000000;
my $minTADsize = 45000;
my $minTADscore = 1.5;
my $sizeTADoptima = '';
my $findTADoptions = "";
my $noOverlapOption = " -noOverlap ";
my $res = 3000;
my $superRes = 15000;
my $minIntReads = 10;
my $stdFilter = 8;
my $minFilter = 0.05;
my $prefix = '';
my $anchorSize = '';
my $anchorLocalBgSize = '';
my $minLoopDist = '';
my $minLoopReads = '';
my $foldLoopLocal = 1.5;
my $foldLoopGlobal = 2.0;
my $poissonLoopLocal = 1e-2;
my $poissonLoopGlobal = 1e-5;
my $genome = '';
my $rawFlag = 0;
my $normTotal = '';
my $badChr = "";
my $badRegionsFile = '';
my $minCoverageTAD = 0.25;
my $minCoverageAnchor = 0.25;
my $balanceFlag = " ";
my $filterFlag = 1;
my $keepUnfilteredFlag = 0;
my $overrideFlag = '';
my @chrs = ();
my $keepFlag = 0;

my $findTADexe = "findTADsAndLoopsFromRelMatrix";
my $analyzeHiCexe = "analyzeHiC";

my $maxCPUs = 1;

sub printCMD {
	print STDERR "\n\tfindTADsAndLoops.pl <cmd> [options]\n";	
	print STDERR "\n\tUsage Modes:\n";
	print STDERR "\n\t\tfindTADsAndLoops.pl find <tagDir> [options]\n";
	print STDERR "\t\t\tFinds TADs and Loops in the given Hi-C tag directory. Output files\n";
	print STDERR "\t\t\twill have the same prefix, and include 2D bed files specifying the locations of the features.\n";
	print STDERR "\n\t\tfindTADsAndLoops.pl score -tad <tad.2D.bed> -loop <loop.2D.bed> -d <tagDir1> [tagDir2] ... [options]\n";
	print STDERR "\t\t\tScores TADs and Loops cross several Hi-C tag directories so that they can be compared across\n";
	print STDERR "\t\t\tdifferent replicates and experimental conditions.\n";
	print STDERR "\t\tNote: When comparing experiments, it is recommended to find TADs/Loops on individual experiments first,\n";
	print STDERR "\t\t\tthen merge the features with merge2Dbed.pl, then 'score' the commbined features across multiple tagDirs\n";
	print STDERR "\n\tCommon Options:\n";
	print STDERR "\t\t-cpu <#> (max CPUs to use)\n";
	print STDERR "\t\t-res <#> (resolution, default: $res)\n";
	print STDERR "\t\t-window <#> (overlapping window size, i.e. superRes, default: $superRes)\n";
	print STDERR "\t\t-balance | -nobalance (balance contact maps, default: -nobalance)\n";
	print STDERR "\t\t-minDist <#> (minimum interaction distance to consider, default: $minDist)\n";
	print STDERR "\t\t-maxDist <#> (maximum interaction distance to consider, default: $maxDist)\n";

	print STDERR "\n\tOptions specific to 'find' mode:\n";
	print STDERR "\t\t-o <output prefix> (by default will place output files within the Tag Directory)\n";
	print STDERR "\t\tTAD options:\n";
	print STDERR "\t\t  -minTADsize <#> (default: $minTADsize)\n";
	print STDERR "\t\t  -minTADscore <#> (minimum inclusion score to call TAD, default: $minTADscore)\n";
	print STDERR "\t\t  -keepOverlappingTADs (by default, program removes overlapping TAD calls)\n";
	print STDERR "\t\tLoop options:\n";
	print STDERR "\t\t  -foldLoopLocalBg <#> (fold threshold for local loop enrichment, default: $foldLoopLocal)\n";
	print STDERR "\t\t  -foldLoopGlobalBg <#> (fold threshold over avg. interactions of similar distance, default: $foldLoopGlobal)\n";
	print STDERR "\t\t  -poissonLoopLocalBg <#> (poisson p-value threshold for differential interactions per loop, default: $poissonLoopLocal)\n";
	print STDERR "\t\t  -poissonLoopGlobalBg <#> (poisson p-value threshold for minimum interactions per loop, default: $poissonLoopGlobal)\n";
	print STDERR "\t\t  Advanced options (probably best to avoid changing):\n";
	print STDERR "\t\t    -minLoopDist <#> (minimum length of contact points to find, default: 3x window size)\n";
	print STDERR "\t\t    -anchorSize <#> (size of region used to calculate loop score, default: window size)\n";
	print STDERR "\t\t    -anchorLocalBgSize <#> (size of region used to calculate local bg, default: 5x window size)\n";
	print STDERR "\t\t    -minLoopReads <#> (minimum number of adjusted reads for contact point, default: use -poissonLoopGlobalBg)\n";
	print STDERR "\t\tArtifact filtering options: (runs filterTADsAndLoops.pl)\n";
	print STDERR "\t\t  -skipFiltering (skip post-feature finding filtering step)\n";
	print STDERR "\t\t  -keepUnfiltered (keep unfiltered version of TADs/Loops)\n";
	print STDERR "\t\t  -p <peak/BED file> (File of genome duplications/bad/blacklisted regions to exclude from Loop/TAD endpoints)\n";
	print STDERR "\t\t  -genome <genome version> (i.e. hg19, hg38)\n";
	print STDERR "\t\t  -badChr <chr> [chr] ... (chromosomes to skip/exclude)\n";
	print STDERR "\t\t  -minCoverageTAD <#> (minimum ratio of read coverage in TAD relative to genome average, def: $minCoverageTAD)\n";
	print STDERR "\t\t  -minCoverageAnchor <#> (minimum ratio of read coverage at anchors relative to average, def: $minCoverageAnchor)\n";
	#print STDERR "\t\t\t-chr <chr> (only analyze this chromosome)\n";
	#print STDERR "\t\t  -estFDR <#> (set loop threshold based on estimated FDR, default: 0.01)\n";
	print STDERR "\t\tOther options:\n";
	print STDERR "\t\t  -insDist <#> (distance across which to calculate insulation ratio, default: $insDist)\n";
	print STDERR "\t\t  -diDist <#> (distance across which to calculate directionality index , default: $diDist)\n";
	print STDERR "\t\t  -override (if res/maxDist are set to produce huge matrices and you have the RAM, you may need to include this)\n";
	print STDERR "\n\tOptions specific to 'score' mode:\n";
	print STDERR "\t\t-o <output prefix> (default: \"score\")\n";
	print STDERR "\t\t-d <HiC TagDir1> [HiC TagDir2] ...\n";
	print STDERR "\t\t-loop <loop 2D BED file> (will score loop interaction rates from this file)\n";
	print STDERR "\t\t-tad <tad 2D BED file> (will score TAD inclusion ratios from this file)\n";
	print STDERR "\t\tLoop score reporting: (default: Hi-C depth normalized counts [use for use with getDiffExpression.pl...])\n";
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
	} elsif ($ARGV[$i] eq '-prefix' || $ARGV[$i] eq '-o') {
		$prefix = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minTADsize') {
		$minTADsize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minTADscore') {
		$minTADscore = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-balance') {
		$balanceFlag = " -balance ";
	} elsif ($ARGV[$i] eq '-nobalance') {
		$balanceFlag = " ";
	} elsif ($ARGV[$i] eq '-keepOverlappingTADs') {
		$noOverlapOption = "";
	} elsif ($ARGV[$i] eq '-diDist') {
		$diDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-chr') {
		push(@chrs, $ARGV[++$i]);
	} elsif ($ARGV[$i] eq '-keep') {
		$keepFlag = 1;
	} elsif ($ARGV[$i] eq '-insDist') {
		$insDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxDist') {
		$maxDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minDist') {
		$minDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-override') {
		$overrideFlag = " -override ";
	} elsif ($ARGV[$i] eq '-minLoopReads') {
		$minLoopReads = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-poissonLoopLocalBg') {
		$poissonLoopLocal = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-poissonLoopGlobalBg') {
		$poissonLoopGlobal = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minLoopDist') {
		$minLoopDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cp' || $ARGV[$i] eq '-loop' || $ARGV[$i] eq '-loops') {
		$int2dBed = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-tad' || $ARGV[$i] eq '-tads' || $ARGV[$i] eq '-TAD' || $ARGV[$i] eq '-TADs') {
		$tad2dBed = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-anchorSize') {
		$anchorSize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-anchorLocalBgSize') {
		$anchorLocalBgSize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-foldLoopLocalBg') {
		$foldLoopLocal = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-foldLoopGlobalBg') {
	  	$foldLoopGlobal = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-keepUnfiltered') {
	  	$keepUnfilteredFlag = 1;
	} elsif ($ARGV[$i] eq '-sizeTADoptima') {
	  	$sizeTADoptima = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-raw') {
		$rawFlag = 1;
	} elsif ($ARGV[$i] eq '-normTotal') {
		$normTotal = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-skipFiltering') {
		$filterFlag = 0;
	} elsif ($ARGV[$i] eq '-p') {
		$badRegionsFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minCoverageTAD') {
		$minCoverageTAD = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minCoverageAnchor') {
		$minCoverageAnchor = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
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
	} elsif ($ARGV[$i] eq '-badChr') {
		$badChr = " -badChr ";
		$i++;
		while ($i<@ARGV) {
			if ($ARGV[$i] =~ /^\-/) {
				$i--;
				last;
			}
			$badChr .= " $_ ";
			$i++;
		}
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

if ($prefix eq '' && $mode eq 'find') {
	my $n = $tagDir;
	$n =~ s/\/+$//;
	$n =~ s/^.*\///;
	$n =~ s/\s/\-/g;
	$prefix = $tagDir  . "/" . $n;
} elsif ($prefix eq '' && $mode eq 'score') {
	$prefix = "score";
}
print STDERR "\tOutput files will use prefix (-o): \"$prefix\"\n";
	


if (scalar(@chrs) < 1) {
	my $chrfiles = `ls -1 "$tagDir"/*.tags.tsv`;
	@chrs = split /\n/,$chrfiles;
	foreach(@chrs) {
		$_ =~ s/^.*\///;
		$_ =~ s/\.tags\.tsv//;
				#print STDERR "$_\n";
	}
	#@chrs = ();
	#push(@chrs,'chr21');
}

if ($int2dBed ne '' || $tad2dBed ne '') {
	my %c = ();
	if ($int2dBed ne '') {
		my $c = getChrsFrom2Dbed($int2dBed);
		foreach(keys %$c) {
			$c{$_}=1;
		}
	}
	if ($tad2dBed ne '') {
		my $c = getChrsFrom2Dbed($tad2dBed);
		foreach(keys %$c) {
			$c{$_}=1;
		}
	}
	my @newchrs = ();
	foreach(@chrs) {
		if (exists($c{$_})) {
			push(@newchrs, $_);
		}
	}
	@chrs = @newchrs;
}

my $rand = rand();

my @resultFiles =();

$anchorSize = $superRes if ($anchorSize eq '');
$anchorLocalBgSize = $superRes*5 if ($anchorLocalBgSize eq '');
$minLoopDist = $superRes*3 if ($minLoopDist eq '');
$sizeTADoptima = $superRes*3 if ($sizeTADoptima eq '');


if ($poissonLoopGlobal eq '' && $mode eq 'find') {
	my $line = `egrep "^gsizeEstimate" $tagDir/tagInfo.txt`;
	if (my @line = split /\=/,$line) {
		my $size = $line[1];
		my $estTests = ($size / $anchorLocalBgSize);
		#my $estTests = ($size / $anchorLocalBgSize) * (($maxDist-$minDist)/$anchorLocalBgSize);
		$poissonLoopGlobal = $estFDR / $estTests;
		print STDERR "\tEstimated poisson p-value threshold = $poissonLoopGlobal (FDR=$estFDR)\n";
	} else {
		print STDERR "\t!!Error, couldn't determine genome size from tagdir ($tagDir)\n";
		print STDERR "\tSetting -poissonLoopGlobal to 1e-5 [default]\n";
		$poissonLoopGlobal = 1e-5;
	}
} 



print STDERR "\tWill analyze chrs: @chrs\n";

#for testing
#@chrs = ("chr22");

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
			$opt .= " -window $superRes ";
			$opt .= " -maxDist $maxDist ";
			$opt .= " -chr $chr ";
			$opt .= " -std $stdFilter ";
			$opt .= " -min $minFilter ";
			$opt .= " -relative ";
			$opt .= " -normArea given ";
			if ($rawFlag) {
				$opt .= " -raw ";
			} else {
				$opt .= " -simpleNorm ";
				$opt .= " -minExpect 1.0 ";
			}
			$opt .= $balanceFlag;
			$opt .= $overrideFlag;
			if ($normTotal  ne '') {
				$opt .= " -normTotal $normTotal ";
			} else {
				$opt .= " -normTotal given ";
			}

			if ($keepFlag) {
				print STDERR "$analyzeHiCexe $td $opt > $rmatrixFile\n";
			}
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

print STDERR "\tProcessing matricies into TADs/Loops...\n";
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

			$opt .= " -anchorSize $anchorSize ";
			$opt .= " -anchorLocalBgSize $anchorLocalBgSize ";
			$opt .= " -minLoopReads $minLoopReads " if ($minLoopReads ne '');
			$opt .= " -minLoopDist $minLoopDist ";
			$opt .= " -foldLocalBg $foldLoopLocal ";
			$opt .= " -foldGlobalBg $foldLoopGlobal ";
			$opt .= " -poissonLocalBg $poissonLoopLocal ";
			$opt .= " -poissonGlobalBg $poissonLoopGlobal ";

			#print STDERR "`$findTADexe find $rmatrix $opt`;\n";
			`$findTADexe find $rmatrix $opt`;
			if ($keepFlag) {
				print STDERR "$findTADexe find $rmatrix $opt\n";
			} else {
				`rm "$rmatrix"`;
			}

		} elsif ($mode eq 'score') {

			my $tagDirStr = "";
			my $z = 0;
			for (my $z=0;$z<@tagDirs;$z++) {
				my $td = $tagDirs[$z];
				$tagDirStr .= " \"$rand.$chr.$z.rmatrix.txt\"";
			}
			my $tadfileChr = $rand . ".tadfile.$chr";
			my $intfileChr = $rand . ".loopfile.$chr";

			my $opt = '';
			if ($int2dBed ne '') {
				my $match = $chr . '\t\d+\t\d+\t' . $chr . '\t';
				`grep -P "$match" $int2dBed > $intfileChr`;
				$opt .= " -loop $intfileChr";
			}
			if ($tad2dBed ne '') {
				my $match = $chr . '\t\d+\t\d+\t' . $chr . '\t';
				`grep -P "$match" $tad2dBed > $tadfileChr`;
				$opt .= " -tad $tadfileChr";
			}
			$opt .= " -prefix $prefix ";
			$opt .= " -matrix $tagDirStr ";
			$opt .= " -minDist $minDist ";
			$opt .= " -anchorSize $anchorSize ";
			$opt .= " -anchorLocalBgSize $anchorLocalBgSize ";

			`$findTADexe score $opt`;
			if ($keepFlag) {
				print STDERR "$findTADexe score $opt\n";
			}

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

	my $totalLoops = 0;
	my $totalTADs = 0;

	if ($tad2dBed ne '') {
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
				$totalTADs++;
			}
			close IN;
			`rm $scoreFile`;
		}
		close OUT;
	}

	if ($int2dBed ne '') {
		open OUT, ">$prefix.loop.scores.txt";
		print OUT "#Loop name\tchr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tScore\tlocalRatio";
		foreach(@tagDirs) {
			print OUT "\t$_ Score";
		}
		print OUT "\n";
		for (my $i=0;$i<@resultFiles;$i++) {
			my $scoreFile = "$resultFiles[$i].loop.scores.txt";
			open IN, "$scoreFile" or print STDERR "\t!!! Can't open $scoreFile";
			my $c = 0;
			while (<IN>) {
				$c++;
				next if ($c < 2);
				print OUT $_;
				$totalLoops++;
			}
			close IN;
			`rm $scoreFile`;
		}
		close OUT;
	}

	print STDERR "\n\tFinished scoring TADs ($totalTADs total) and Loops ($totalLoops total)\n\n";
	exit;
}

my $prePrefix = $prefix . "_2";
open BED2D, ">$prePrefix.tad.2D.bed";
print BED2D "#TAD chr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tInclusionRatio\tInclusionRatio\n";
open BED2DINT, ">$prePrefix.loop.2D.bed";
print BED2DINT "#Loop chr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tLoopScore\tLocalFold\n";
open BEDGRAPH, ">$prefix.DI.bedGraph";
print BEDGRAPH "track name=\"$prefix DI ($tagDir)\" type=bedGraph visibility=2\n";
open INSBEDGRAPH, ">$prefix.Insulation.bedGraph";
print INSBEDGRAPH "track name=\"$prefix Insulation Ratio ($tagDir)\" type=bedGraph visibility=2\n";
#open TAD, ">$prefix.tad.txt";
#print TAD "name\tchr\tstart\tend\tdirection\tInclusionRatio\tLength\n";
my $c = 0;
my $totalTADs=0;
my $totalLoops=0;

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

	open IN, "$f.loop.2D.bed";
	$c=0;
	while (<IN>) {
		$c++;
		#next if ($c==1);
		print BED2DINT "$_";
		$totalLoops++;
	}
	close IN;
	`rm $f.loop.2D.bed`;
}
close BED2D;
close BED2DINT;
#close TAD;
close BEDGRAPH;
close INSBEDGRAPH;

#filtering interactions
if ($filterFlag) {
	print STDERR "\n\tFiltering TADs/Loops...\n\n";
	my $filterOpt = " -minCoverageTAD $minCoverageTAD ";
	$filterOpt .= " -minCoverageAnchor $minCoverageAnchor ";
	$genome = "none" if ($genome eq '');
	$filterOpt .= " -genome $genome " if ($genome ne '');
	$filterOpt .= " -p \"$badRegionsFile\" " if ($badRegionsFile ne '');
	$filterOpt .= $badChr;
	`filterTADsAndLoops.pl -loop "$prePrefix.loop.2D.bed" -tad "$prePrefix.tad.2D.bed" -o $prefix -d "$tagDir" $filterOpt`;
	#print STDERR "`filterTADsAndLoops.pl -loop $prePrefix.loop.2D.bed -tad $prePrefix.tad.2D.bed -o $prefix -d $tagDir $filterOpt`;\n";
	if ($keepUnfilteredFlag) {
		`mv "$prePrefix.loop.2D.bed" "$prefix.unfiltered.loop.2D.bed"`;
		`mv "$prePrefix.tad.2D.bed" "$prefix.unfiltered.tad.2D.bed"`;
	} else {
		`rm "$prePrefix.loop.2D.bed" "$prePrefix.tad.2D.bed"`;
	}
	$totalTADs = `wc -l "$prefix.tad.2D.bed"`;
	chomp $totalTADs;
	$totalTADs =~ s/^\s*//;
	$totalTADs =~ s/\s+.*//;
	$totalLoops = `wc -l "$prefix.loop.2D.bed"`;
	chomp $totalLoops;
	$totalLoops =~ s/^\s*//;
	$totalLoops =~ s/\s+.*//;
} else {
	`mv "$prePrefix.loop.2D.bed" "$prefix.loop.2D.bed"`;
	`mv "$prePrefix.tad.2D.bed" "$prefix.tad.2D.bed"`;
}



`bed2DtoUCSCbed.pl "$prefix.tad.2D.bed" > "$prefix.ucsc.bed"`;
`bed2DtoUCSCbed.pl "$prefix.loop.2D.bed" >> "$prefix.ucsc.bed"`;
`cat "$prefix.DI.bedGraph" "$prefix.Insulation.bedGraph" >> "$prefix.ucsc.bed"`;
`bed2DtoUCSCbed.pl "$prefix.tad.2D.bed" -peaks > "$prefix.tad.peaks.txt"`;
`bed2DtoUCSCbed.pl "$prefix.loop.2D.bed" -interactions > "$prefix.interactions.txt"`;
print STDERR "\n\n\tTotal TADs: $totalTADs\n";
print STDERR "\n\tTotal Loops : $totalLoops\n";
print STDERR "\n";


sub getChrsFrom2Dbed {
	my ($file) = @_;
	my %c = ();
	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^\#/);
		my @line = split /\t/;
		$c{$line[0]} = 1;
		$c{$line[4]} = 1;
	}
	close IN;
	return \%c;
}

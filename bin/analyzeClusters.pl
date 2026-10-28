#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";


# Copyright 2009-2021 Christopher Benner <cbenner@ucsd.edu>
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



use HomerConfig;

my $maxClusterCount = 10000;

my $gofile= "biological_process.txt";

my $stopFlag = 0;
my $inputFile = '';
my $gtrFile = '';
my $maxCPUs = 1;
my $maxpCPUs = 1;
my $annCols = 1;
my $peakFlag = 0;
my $genome = '';
my $org = '';
my $directory = '';
my $thresh = -0.7;
my $membershipFile = "";
my $noDenovo = '-nomotif';
my $size = 100;
my $S = 0;
my $lens = '8,10,12';
my $topResults = 20;
my $minGenes = 500;
my $minExp = -1e21;
my $minAbsExp = '';
my $minDiff = -1e21;
my $quantileFlag = 0;
my $centerFlag = 0;
my $maxGOtermGenes = 500;
my $logxformFlag = 0;
my $pseudo = 0;
my $bgPeakFile = '';
my $lastSampleInput = 0;
my $recastAsOligoFile = 0;
my $norevopp = 0;
my $peakCmpStr = "";
my @peakCmp=();
my $peakCmpDist = 400;
my $tssGOfile = '';
my $peakTSSGODist = 400;
my $peakTSSGOStrand = 'same';
my $includeNoBg = 0;
my $minThresh = -0.9;
my $maxThresh = 0.5;
my $incThresh = 0.05;

my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";
my $tmpFile4 = $rand . ".4.tmp";
my $tmpFile5 = $rand . ".5.tmp";

sub printCMD {
	print STDERR "\n\tusage: analyzeClusters.pl [options] -i <data file>\n";

	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-o <OutputDirectory> (output directory)\n";
	print STDERR "\t\t-i <data file> (data file to cluster/analyze, tab-delimited)\n";
	print STDERR "\t\t\t-peaks (annotatePeaks output file)\n";
	print STDERR "\t\t\t-repeats (analyzeRepeats.pl output file)\n";
	print STDERR "\t\t\t-annCols (number of annotation columns if not peak or expression file, def:1)\n";
	print STDERR "\n\t\tClustering options:\n";
	print STDERR "\t\t\t-gtr <gtr file> (will use this and skip clustering, i.e. directory/out.gtr)\n";
	print STDERR "\t\t\t-membership <file> (membership file assigning clusters)\n";
	print STDERR "\t\t\t-thresh <#> (threshold used for clusters (-1 to 1), default: $thresh)\n";
	print STDERR "\t\t\t-min <#> (minimum number of entries per cluster, def: $minGenes)\n";
	print STDERR "\n\t\tData value normalization/filtering:\n";
	print STDERR "\t\t\tOrder: pseudo,log, quantile, filtering, centering\n";
	print STDERR "\t\t\t-pseudo <#> (add pseudo count to data, applied before transformations etc.)\n";
	print STDERR "\t\t\t-log (log transform the data)\n";
	print STDERR "\t\t\t-quantile (quantile normalize matrix before clustering)\n";
	print STDERR "\t\t\t(note, filtering happens after quantile/log transformations)\n";
	print STDERR "\t\t\t-minExp <#> (minimum value of the row max to include in analysis)\n";
	#print STDERR "\t\t\t\t-minExpAll <#> (exlude if any of the samples have exp below this value)\n";
	print STDERR "\t\t\t-minAbsExp <#> (minimum value of the row max or abs(min) to include in analysis)\n";
	print STDERR "\t\t\t-minDiff <#> (minimum difference between min and max value in each row to include)\n";
	print STDERR "\t\t\t-center (center rows by average)\n";
	print STDERR "\t\t\t-seedClusters <#> (default: $maxClusterCount)\n";
	print STDERR "\t\t\t-lastSampleInput (last sample is control - remove if this value is higher than row max)\n";
	print STDERR "\n\t\tOther analysis options:\n";
	print STDERR "\t\t\t-cpu <#> (number of concurrent clusters to analyze, def: $maxCPUs)\n";
	print STDERR "\t\t\t-p <#> (number of CPUs to use per cluster, def: $maxpCPUs)\n";
	print STDERR "\t\t\t-genome <version> (Genome version to use for peak files)\n";
	print STDERR "\t\t\t-org <organism> (Organism to use for gene based analysis)\n";
	print STDERR "\t\t\t-size <#> (size of region for motif finding, def: 100 [peaks] or -200,50 [genes])\n";
	print STDERR "\t\t\t-S <#> (number of de novo motifs to find, default: none - known enrichment only)\n";
	print STDERR "\t\t\t-len <#,#,...> (lengths of motifs to look for)\n";
	print STDERR "\t\t\t-bg <bg peak file for motifs> (i.e. use input file here to get differential cluster motifs)\n";
	print STDERR "\t\t\t\t-includeNoBg (rum motifs with and without background file)\n";
	print STDERR "\t\t\t-maxGOgenes <#> (only GO categories with less than this number of genes considered, $maxGOtermGenes)\n";
	print STDERR "\t\t\t-peakCmp <peakFile1> [peakfile2]... (compare cluster to peaks, report top one in summary)\n";
	print STDERR "\t\t\t-tssGO <annotated TSS file> (only perform GO with overlapping TSS regions)\n";
	print STDERR "\t\t\t-oligos <# oligo length> (length of oligos to use, i.e. 8)\n";
	print STDERR "\t\t\t\t-norevopp (don't collapse reverseOpposites)\n";
	print STDERR "\t\t\t-stop (Stops analysis after clustering so that one can change the threshold/other paramters)\n";
	print STDERR "\n";
	exit;

}
if (@ARGV < 1) {
	printCMD();
}
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-i') {
		$inputFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minExp') {
		$minExp = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minAbsExp') {
		$minAbsExp = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minDiff') {
		$minDiff = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-stop') {
		$stopFlag = 1;
	} elsif ($ARGV[$i] eq '-len') {
		$lens = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-quantile') {
		$quantileFlag = 1;
	} elsif ($ARGV[$i] eq '-includeNoBg') {
		$includeNoBg = 1;
	} elsif ($ARGV[$i] eq '-center') {
		$centerFlag = 1;
	} elsif ($ARGV[$i] eq '-peakCmp') {
		$i++;
		while ($i<@ARGV && $ARGV[$i] !~ /^\-/) {
			push(@peakCmp, $ARGV[$i]);
			$peakCmpStr .= " " . $ARGV[$i];
			$i++;
		}
		if ($i<@ARGV && $ARGV[$i] =~ /^\-/) {
			$i--;
		}
		next;
	} elsif ($ARGV[$i] eq '-o') {
		$directory = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-gtr') {
		$gtrFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-membership') {
		$membershipFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-p') {
		$maxpCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-seedClusters') {
		$maxClusterCount = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxGOgenes') {
		$maxGOtermGenes = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-size') {
		$size = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$minGenes = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-log') {
		$logxformFlag = 1;
	} elsif ($ARGV[$i] eq '-tssGO') {
		$tssGOfile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bg') {
		$bgPeakFile = " -bg \"" . $ARGV[++$i] . "\" ";
	} elsif ($ARGV[$i] eq '-pseudo') {
		$pseudo = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-S') {
		$noDenovo = '';
		$S = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-annCols') {
		$annCols = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-peak' || $ARGV[$i] eq '-peaks') {
		$annCols = 19;
		$peakFlag = 1;
	} elsif ($ARGV[$i] eq '-repeats') {
		$annCols = 8;
		$peakFlag = 0;
	} elsif ($ARGV[$i] eq '-thresh') {
		$thresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-org') {
		$org = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-lastSampleInput') {
		$lastSampleInput = 1;
	} elsif ($ARGV[$i] eq '-oligos') {
		$recastAsOligoFile = $ARGV[++$i];
		print STDERR "\tWill recast regions as oligos\n";
	} elsif ($ARGV[$i] eq '-norevopp') {
		$norevopp = 1;
	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
		$peakFlag = 1;
	} else {
		print STDERR "!!! What is \"$ARGV[$i]\"?\n";
		printCMD();
	}
}


if ($directory eq '') {
	print STDERR "!!! -o <outputdirectory> is required!!!\n";
	exit;
}
if ($inputFile eq '' && $gtrFile eq '') {
	print STDERR "!!! -i <inputdata> or -gtr <gtr file>  is required!!!\n";
	exit;
}
if ($org eq '' && $genome eq '') {
	print STDERR "!!! -org <organism> or -genome <version> is required!!!\n";
	exit;
}
if ($peakFlag && ($inputFile eq '' || $genome eq '')) {
	print STDERR "!!! If analyzing regions/peaks you need an input file (-i) and genome (-genome)!!!\n";
	exit;
}
if ($org eq '' && $genome =~ /^hg/) {
	$org = 'human';
}
if ($org eq '' && $genome =~ /^mm/) {
	$org = 'mouse';
}
if ($tssGOfile ne '' && $org eq '') {
	print STDERR "!!! Error - '-org' required when using a tssGO file\n";
	exit;
}

`mkdir -p "$directory"`;
my $normFile = "$directory/normalized.txt";
my $oligoFile = "$directory/oligoFile.txt";

my $logFile = "$directory/log.txt";
open LOG, ">$logFile";
print LOG "inputFile = $inputFile\n";
print LOG "maxClusterCount = $maxClusterCount\n";
print LOG "gofile =  $gofile\n";
print LOG "gtrFile = $gtrFile\n";
print LOG "maxCPUs = $maxCPUs\n";
print LOG "maxpCPUs = $maxpCPUs\n";
print LOG "annCols = $annCols\n";
print LOG "peakFlag = $peakFlag\n";
print LOG "genome = $genome\n";
print LOG "org = $org\n";
print LOG "directory = $directory\n";
print LOG "thresh = $thresh\n";
print LOG "membershipFile = $membershipFile\n";
print LOG "noDevno = $noDenovo\n";
print LOG "size = $size\n";
print LOG "S = $S\n";
print LOG "lens = $lens\n";
print LOG "topResults = $topResults\n";
print LOG "minGenes = $minGenes\n";
print LOG "minExp = $minExp\n";
print LOG "minAbsExp = $minAbsExp\n";
print LOG "minDiff = $minDiff\n";
print LOG "quantileFlag = $quantileFlag\n";
print LOG "centerFlag = $centerFlag\n";
print LOG "maxGOtermGenes = $maxGOtermGenes\n";
print LOG "logxformFLag = $logxformFlag\n";
print LOG "pseudo = $pseudo\n";
print LOG "bgPeakFile = $bgPeakFile\n";
print LOG "lastSampleInput = $lastSampleInput\n";


if ($inputFile ne '' && ($gtrFile eq '' && $membershipFile eq '')) {
	`cp "$inputFile" "$normFile"`;
}

#cluster
if ($gtrFile eq '' && $membershipFile eq '' && $inputFile ne '') {

	if ($quantileFlag || $logxformFlag) {
		my $x = $annCols+1;
		`cut -f1,$x- $normFile > $tmpFile2`;
		`cut -f1-$annCols $normFile > $tmpFile3`;
		if ($logxformFlag) {
			print STDERR "\tLog2 tranforming the data...\n";
			print LOG "\tLog2 tranforming the data...\n";
			my $opt = '';
			if ($pseudo > 0) {
				$opt = " -add $pseudo ";
			}
			`logxform.pl $tmpFile2 $opt > $tmpFile4`;
			`mv $tmpFile4 $tmpFile2`;
		}
		if ($quantileFlag) {
			print STDERR "\tPerforming quantile normalization...\n";
			print LOG "\tPerforming quantile normalization...\n";
			#print STDERR "`quantileNorm.pl $tmpFile2 > $tmpFile4`;\n";
			`quantileNorm.pl $tmpFile2 > $tmpFile4`;
			`mv $tmpFile4 $tmpFile2`;
		}
		`addDataHeader.pl $tmpFile3 $tmpFile2 > $normFile`;
		`rm $tmpFile2 $tmpFile3`;
	}
	if ($centerFlag || $minExp > -1e20 || $minDiff > -1e20 || $minAbsExp ne '' || $lastSampleInput) {
		print STDERR "\tFiltering...\n";
		open IN, $normFile;
		open OUT, ">$tmpFile2";
		my $c = 0;
		my $total = 0;
		my $kept = 0;
		while (<IN>) {
			$c++;
			chomp;
			s/\r//g;
			my @line = split /\t/;
			if ($c == 1) {
				print OUT "$line[0]";
				for (my $i=1;$i<@line-$lastSampleInput;$i++) {
					print OUT "\t$line[$i]";
				}
				print OUT "\n";
				next;
			}
			my $rsum = 0;
			my $rmin = 1e20; 
			my $rmax = -1e20;
			my $rN = 0;
			for (my $i=$annCols;$i<@line-$lastSampleInput;$i++) {
#print STDERR "\t|$line[$i]|\n";
				$rmax = $line[$i] if ($line[$i] > $rmax);
				$rmin = $line[$i] if ($line[$i] < $rmin);
				$rsum += $line[$i];
				$rN++;
			}
			my $inputLevel = -1e10;
			if ($lastSampleInput) {
				$inputLevel = $line[@line-1];
			}
			$total++;
			if ($rmax < $minExp) {
					#print STDERR "\tFiltering rmax < minExp ($rmax,$rmin,$rsum,$rN, $minExp)\n";
				next;
			}
			if ($minAbsExp ne '') {
				if (abs($rmax) < $minAbsExp && abs($rmin) < $minAbsExp) {
					next;
				}
			}
			next if ($rmax-$rmin < $minDiff);
			next if ($rN < 1);
			my $avg = $rsum/$rN;
			if ($lastSampleInput) {
				if ($avg < $inputLevel) {
					next;
				}
			}
			if ($centerFlag) {
				for (my $i=$annCols;$i<@line;$i++) {
					$line[$i] -= $avg;
				}
			}
			$kept++;
			print OUT "$line[0]";
			for (my $i=1;$i<@line-$lastSampleInput;$i++) {
				print OUT "\t$line[$i]";
			}
			print OUT "\n";
		}
		close OUT;
		close IN;
		print STDERR "\tFiltering: kept $kept out of $total\n";
		print LOG "\tFiltering: kept $kept out of $total\n";
		`mv $tmpFile2 $normFile`;
	}

	if ($recastAsOligoFile) {
		recastAsOligos($recastAsOligoFile,$normFile, $oligoFile,$size, $genome, $norevopp);
	}

	my $lineCount = `wc -l < "$normFile"`;
	my $sub = '';
	if ($lineCount > $maxClusterCount) {
		$sub = "-sub $maxClusterCount";
	}
	`homerTools cluster -i "$normFile" -annCols $annCols $sub -o "$directory/out"`;
	$gtrFile = "$directory/out.gtr";

	open THRESH, ">$directory/thresholds.txt";
	print STDERR "\nCluster Threshold Analysis:\n\t\tThreshold\tValidClusters\tTotalClusters\tFraction of Genes\n";
	print THRESH "Threshold\tValidClusters\tTotalClusters\tFraction of Genes\n";
	for (my $i=$minThresh;$i<=$maxThresh;$i+=$incThresh) {
		my $checkThresh = $i;
		my ($nvc,$tc,$nvg,$tg) = checkThreshResults($checkThresh,$gtrFile);
		$tg = 1 if ($tg == 0);
		#print STDERR "$nvc\t$tc\t$nvg\t$tg\n";
		my $r = sprintf("%.2lf",$nvg/$tg*100);
		print STDERR "\t\t$checkThresh\t$nvc\t$tc\t$r%\n";
		print THRESH "$checkThresh\t$nvc\t$tc\t$r%\n";
	}
	close THRESH;

}

if ($membershipFile eq '') {
	`homerTools cluster -gtr "$gtrFile" -thresh $thresh > "$directory/clusterMembership.txt"`;
	$membershipFile = "$directory/clusterMembership.txt";
}


my %membership = ();
my %clusters = ();
my %counts = ();
open IN, $membershipFile;
my $numClusters = 0;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if ($line[1] eq '-1');
	$membership{$line[0]} = $line[1];
	if (!exists($clusters{$line[1]})) {
		$numClusters++;
		my @a = ();
		$clusters{$line[1]}=\@a;
	}
	push(@{$clusters{$line[1]}}, $line[0]);
	$counts{$line[1]}++;
}
close IN;

my $numValid = 0;
my @clusters = ();
my $total = 0;
my $validTotal = 0;
foreach(keys %counts) {
	$total += $counts{$_};
	if ($counts{$_} >= $minGenes) {
		$numValid++;
		$validTotal+=$counts{$_};
		push(@clusters, $_);
	}
}

$r = sprintf("%.2lf",$validTotal/$total*100);

print STDERR "\tFound $numValid (of $numClusters total) clusters (threh=$thresh, minPerCluster=$minGenes)\n";
print STDERR "\t\tComprising $validTotal of $total genes/regions ($r%)\n";
print STDERR "\tSleeping for 5 secs if you want to modify -thresh or -min\n";

print LOG "\tFound $numValid (of $numClusters total) clusters (threh=$thresh, minPerCluster=$minGenes)\n";
print LOG "\t\tComprising $validTotal of $total genes/regions ($r%)\n";
print LOG "\tSleeping for 5 secs if you want to modify -thresh or -min\n";

if ($stopFlag) {
	print STDERR "\t\"-stop\" flag used - quitting...\n";
	exit;
}

`sleep 5`;


my %data = ();
my @dataHeader = '';
if ($normFile ne '') {
	#print STDERR "Here? normFile = $normFile\n";
	open IN, $normFile;
	my $c = 0;
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($c==1) {
			@dataHeader = @line;;
			next;
		}
		$data{$line[0]} = \@line;
	}
	close IN;
}

#create input files
@clusters = sort {$a <=> $b} @clusters;
my @files = ();
foreach(@clusters) {
	my $c = $_;
	my $file = "$directory/$c.cluster.txt";
	push(@files, $file);
	open OUT, ">$file";
	foreach(@{$clusters{$c}}) {
		my $id = $_;
		if ($peakFlag) {
			print OUT "$data{$id}->[0]";
			for (my $i=1;$i<5;$i++) {
				print OUT "\t$data{$id}->[$i]";
			}
			print OUT "\n";
		} else {
			print OUT "$id\n";
		}
	}
	close OUT;
}

my $cpus = 0;


my @geneListFiles = ();

foreach(@files) {
	my $file = $_;
	my $pid = fork();
	$cpus++;
	if ($pid==0) {
		if ($peakFlag) {
			`annotatePeaks.pl $file $genome -size 600 -hist 1 -m ~/homer/motifs/sp1.motif ~/homer/motifs/nfy.motif ~/homer/motifs/ap1.motif > $file.mfreq.txt`; 
			my $gocmd = " -go \"$file.motifs/\" ";
			if ($tssGOfile ne '') {
				$gocmd = '';
				my $mopt = " -d $peakTSSGODist ";
				if ($peakTSSGOStrand eq 'same') {
					$mopt .= " -strand ";
				}
				`mergePeaks $mopt $tssGOfile $file -cobound 1 -prefix $file.tssGO`;
				open OUT, ">$file.tssGenes.txt";
				open IN, "$file.tssGO.coBoundBy1.txt";
				while (<IN>) {
					s/\-\-\d+\t/\t/;
					print OUT $_;
				}
				close IN;
				close OUT;
				`rm $file.tssGO.coBoundBy[01].txt`;
				push(@geneListFiles, "$file.tssGenes.txt");
				`findGO.pl $file.tssGenes.txt $org \"$file.motifs/\" -p $maxpCPUs`;

			}
			`annotatePeaks.pl $file $genome $gocmd -annStats $file.annStats.txt > $file.annotation.txt`; 
			`findMotifsGenome.pl $file $genome $file.motifs/ -p $maxpCPUs -size $size -len $lens -S $S $noDenovo $bgPeakFile`;
			if ($bgPeakFile ne '' && $includeNoBg) {
				`findMotifsGenome.pl $file $genome $file.motifsNoBg/ -p $maxpCPUs -size $size -len $lens -S $S $noDenovo`;
			}
			if ($peakCmpStr ne '') {
				`mergePeaks $file $peakCmpStr -d $peakCmpDist -cobound 1 -matrix $file.peakCmp`;
				`rm coBoundBy0.txt coBoundBy1.txt`;
			}
		} else {
			`findMotifs.pl $file $org $file.motifs/ -p $maxpCPUs -len $lens -S $S $noDenovo`;
		}
		exit;
	}	
	if ($cpus >= $maxCPUs) {
		wait();
		$cpus--;
	}
}
my $id = 0;
while ($id >= 0) {
	$id = wait();
}

if ($peakCmpStr ne '') {
	`addDataAll.pl $directory/*.peakCmp.logPvalue.txt > $tmpFile`;
	open IN, $tmpFile;
	open OUT, ">$directory/peakCmp.matrix.txt";
	my $z =0;
	while (<IN>) {
		$z++;
		chomp;
		my @line = split /\t/;
		if ($z == 1) {
			foreach(@line) {
				$_ =~ s/^.*\///;
				$_ =~ s/.txt.peakCmp.logPvalue.txt-vs Reference Peaks//;
			}
		}
		$line[0] =~ s/^.*\///;
		$line[0] =~ s/\.peaks\.txt$//;
		print OUT "$line[0]";
		for (my $i=1;$i<@line;$i++){ 
			print OUT "\t$line[$i]";
		}
		print OUT "\n";
	}
	close IN;
	close OUT;
	`rm $tmpFile`;

}

if (scalar(@geneListFiles) > 0) {

	my $longest = 0;
	my %geneLists = ();
	foreach(@geneListFiles) {
		my @a = ();
		my $file = $_;
		open IN, $file;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			next if ($line[0] eq '' || $line[0] =~ /^\#/);
			push(@a, $line[0]);
		}
		close IN;
		$longest = scalar(@a) if (scalar(@a) > $longest);
		$geneLists{$file}=\@a;
	}
			

	open OUT, ">$directory/metascape.multiInput.csv";
	my $z = 0;
	foreach(@geneListFiles) {
		$z++;
		my $name = $_;
		my $file = $_;
		$name =~ s/^.+\///;
		$name =~ s/.txt.tssGenes.txt//;
		print OUT "," if ($z > 1);
		print OUT "$name";
	}
	print OUT "\n";

	for (my $i=0;$i<$longest;$i++) {
		my $z = 0;
		foreach(@geneListFiles) {
			$z++;
			print OUT "," if ($z > 1);
			my $file = $_;
			if ($i < scalar(@{$geneLists{$file}})) {
				print OUT $geneLists{$file}->[$i];
			}
		}
		print OUT "\n";
	}	
	close OUT;
}


open OUT, ">$directory/summary.txt";
open OUT2, ">$directory/summary.clust.txt";
print OUT "Cluster\tCount\tBestGO\tBestMotif\tTopGOresults\tTopMotifResults\tTopPeakCmpResults";
print OUT2 "Cluster";
my $numExps=0;
if ($inputFile ne '') {
	for (my $i=$annCols;$i<@dataHeader;$i++) {
		print OUT "\t$dataHeader[$i]";
		print OUT2 "\t$dataHeader[$i]";
		$numExps++;
	}
}
print OUT "\n";
print OUT2 "\n";

for (my $i=0;$i<@files;$i++) {
	my $clust = $clusters[$i];
	my $file = $files[$i];

	my $topMotifStr = "";
	my $motifStr = "";

	open IN, "$file.motifs/knownResults.txt";
	$c=0;
	while (<IN>) {
		$c++;	
		chomp;
		s/\r//g;
		my @line = split /\t/;
		$line[0] =~ s/\/.*$//;
		if ($c == 2) {
			$topMotifStr .= $line[0] . "|" . $line[1] . "|" . $line[2];
		}
		if ($c < $topResults+2) {
			$motifStr .= "," if ($c > 2);
			$motifStr .= $line[0] . "|" . $line[1] . "|" . $line[2];
		}
	}
	close IN;
			
	my $topGoStr = "";	
	my $goStr = '';
	my $goCount = 0;
	open IN, "$file.motifs/$gofile";
	$c=0;
	while (<IN>) {
		$c++;	
		next if ($c < 2);
		chomp;
		s/\r//g;
		my @line = split /\t/;
		$line[2] = sprintf("%.2e",$line[2]);
		next if ($line[4] >= $maxGOtermGenes);
		if ($goCount < $topResults) {
			$goCount++;
			$goStr .= "," if ($goCount > 2);
			$goStr .= $line[1] . "|" . $line[0] . "|" . $line[2];
			if ($topGoStr eq '') {
				$topGoStr .= $line[1] . "|" . $line[0] . "|" . $line[2];
			}
		} else {
			last;
		}
	}
	close IN;

	#open IN, "$file.annStats.txt";
	#while (<IN>) {
	#}
	#close IN;

	my $peakCmpOutStr = "NA";
	if ($peakCmpStr ne '') {
		open IN, "$file.peakCmp.logPvalue.txt";
		my %res = ();
		my $z = 0;
		while (<IN>) {
			$z++;
			next if ($z < 2);
			chomp;
			s/\r//g;
			my @line = split /\t/;
			my $name = $line[0];
			$name =~ s/^.+\///;
			$name =~ s/\.txt//;
			$name =~ s/\.peaks//;
			$res{$name} = $line[1];
		}
		my @best = sort {$res{$a} <=> $res{$b}} keys %res;
		if (@best > 0) {
			$peakCmpOutStr = $best[0]  . "|" . $res{$best[0]};
		}
	}
		
	my $N = 0;	
	my @exp = ();
	for (my $j=0;$j<$numExps;$j++) {
		push(@exp, 0);
	}
	foreach(@{$clusters{$clust}}) {
		$N++;
		my $id = $_;
		if ($inputFile ne '') {
			for (my $j=$annCols;$j<$annCols+$numExps;$j++) {
				$exp[$j-$annCols] += $data{$id}->[$j];
			}
		}
	}
	if ($N > 0 && $inputFile ne '') {
		for (my $j=0;$j<$numExps;$j++) {
			$exp[$j] /= $N;
		}
	}
	print OUT "$clust\t$N\t$topGoStr\t$topMotifStr\t$goStr\t$motifStr\t$peakCmpOutStr";
	print OUT2 "C=$clust-N=$N-$topMotifStr-$topGoStr-$peakCmpOutStr";
	if ($inputFile ne '') {
		foreach(@exp) {
			print OUT "\t$_";
			print OUT2 "\t$_";
		}
	}
	print OUT "\n";
	print OUT2 "\n";
}
close OUT;
close OUT2;
close LOG;
`homerTools cluster -i "$directory/summary.clust.txt" -o "$directory/summary"`;
`zip -j $directory/s.zip $directory/summary.cdt $directory/summary.gtr $directory/out.cdt $directory/out.gtr`;


sub recastAsOligos {
	my ($oligoLen,$normFile, $oligoFile, $size, $genome, $norevopp) = @_;

	print STDERR "\tRecasting data as an oligo file... ($oligoFile)\n";


	`adjustPeakFile.pl $normFile -size $size > $tmpFile`;	
	`homerTools extract $tmpFile $genome > $tmpFile2`;
	

	my %regions = ();
	open OUT, ">$oligoFile";
	open IN, $normFile;
	my $c = 0;
	my $numSamples = 0;
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($c == 1) {
			print OUT "oligo";
			for (my $i=$annCols;$i<@line-$lastSampleInput;$i++) {
				print OUT "\t$line[$i]";
				$numSamples++;
			}
			print OUT "\n";
			next;
		}
		my @d = ();
		for (my $i=$annCols;$i<@line-$lastSampleInput;$i++) {
			push(@d, $line[$i]);
		}
		$regions{$line[0]} = \@d;
	}
	close IN;

	my %oligos = ();
	open IN, $tmpFile2;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $line[0];
		my $seq = $line[1];
		if (!exists($regions{$id})) {
			print STDERR "!!! Error - couldn't find region ID: $id\n";
			next;
		}
		my $len = length($seq);
		my %o = ();
		for (my $i=0;$i<$len-$oligoLen;$i++) {
			my $s = substr($seq,$i,$oligoLen);
			next if ($s =~ /N/i);
			if (!$norevopp) {
				if (!exists($o{$s})) {
					my $rv = HomerConfig::revopp($s);
					if (exists($o{$rv})) {
						$s = $rv;
					}
				}
			} 
			$o{$s}++;
		}
		
		foreach(keys %o) {
			my $s = $_;
			if (!$norevopp) {
				if (!exists($oligos{$s})) {
					my $rv = HomerConfig::revopp($s);
					if (exists($oligos{$rv})) {
						$s = $rv;
					}
				}
			}
			if (!exists($oligos{$s})) {
				my @a = ();
				for (my $i=0;$i<$numSamples;$i++) {
					push(@a, 0);
				}
				$oligos{$s} = {n=>0,v=>\@a};
			}

			for ($i=0;$i<$numSamples;$i++) {
				$oligos{$s}->{'v'}->[$i] += $regions{$id}->[$i];
			}
			$oligos{$s}->{'n'}++;
		}
	}
	close IN;

	foreach (keys %oligos) {
		my $s = $_;
		print OUT "$s";
		foreach(@{$oligos{$s}->{'v'}}) {
			$_ /= $oligos{$s}->{'n'};
			print OUT "\t$_";
		}
		print OUT "\n";
	}
	exit;

}




sub checkThreshResults {
	my ($t,$gtrFile) = @_;

	`homerTools cluster -gtr "$gtrFile" -thresh $t > "$tmpFile"`;

	my %membershipTmp = ();
	my %clusters = ();
	my %counts = ();

	open IN, $tmpFile;
	my $totalClusters = 0;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if ($line[1] eq '-1');
		$membership{$line[0]} = $line[1];
		if (!exists($clusters{$line[1]})) {
			$totalClusters++;
			my @a = ();
			$clusters{$line[1]}=\@a;
		}
		push(@{$clusters{$line[1]}}, $line[0]);
		$counts{$line[1]}++;
	}
	close IN;
	`rm "$tmpFile"`;

	my $numValidClusters = 0;
	my $totalGenes = 0;
	my $numValidGenes = 0;
	foreach(keys %counts) {
		$totalGenes += $counts{$_};
		if ($counts{$_} >= $minGenes) {
			$numValidClusters++;
			$numValidGenes+=$counts{$_};
		}
	}
	return ($numValidClusters,$totalClusters, $numValidGenes, $totalGenes);
}	

#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";

#
#
# Copyright 2009 - 2024 Christopher Benner <cbenner@ucsd.edu>
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
#


use POSIX;
use HomerConfig;
use Statistics;
use HomerSeqBias;

$scoreTestMode = "";
$tssStrandTestMode = "";

my $cmd = "";
foreach(@ARGV) {
	$cmd .= $_ . " ";
}

my $byIntervals=1;
my $rflag = 1;
my $randomizeVariantAssignments = 0;
my $preserveActiveStatusForRand = 1;

my $minDeltaChange = 1e-10;
my $outputDir = '';
my $maggieExe = "conda run -n maggie python /gpfs/data01/bennerlab/home/cbenner/software/maggie/maggie/bin/maggie_fasta_input.py";
my $maxCPUs = "56";
my $maxCPUs2 = "10";
my $maxCPUs3 = "50";
my $minSeqLen = 20;
$liftOverSize = 10;
my $strand = ".";
my $adjustMotifThreshold = '';
my $onlyIncludeSNV = 1;
my $onlyIncludeIndels = 0;
my $onlyIncludeBinaryVariants = 1;
my $useSignedPvalues = 0;
my $pvalueThreshold = 1;
my $mafThreshold = 0.05;
my $absSlopeMin = -1e100;
my $defmscore = -100;
my $epsPvalue = 1e-299;
my $sameNflag = 1;
my $zscore = 0;
my $motifRandOrder = 1;
my $mask = "";
my $markovOrder = 1;
my $randomizationsPerSequence = 5;
my $existingRandomizationFile = '';
my $maxIndelLength = 10;
my $subsetFile = '';
my $targetNumOfRandDiffs = 1e5;

my $defStart = -150;
my $defEnd = 150;
my $defStep = 10;
my $defWindow = 30;

my %nucIndex = ();
$nucIndex{'A'} = 0;
$nucIndex{'C'} = 1;
$nucIndex{'G'} = 2;
$nucIndex{'T'} = 3;

my $maggieMotifs = "  --motifPath motifs/jaspar/ ";

my $actSitesFile = "mutatedSites.activating.txt";
my $repSitesFile = "mutatedSites.repressing.txt";

sub printCMD {
	print STDERR "\n\tmaggiesPacifier.pl <genotype/phenotype/size inputs> <-m motifFile> <options>\n";

	print STDERR "\n\tIf using two genomes linked with liftOver files:\n";
	print STDERR "\t\t-g1 <genome fasta file for 1st strain>\n";
	print STDERR "\t\t-g2 <genome fasta file for 2nd strain>\n";
	print STDERR "\t\t-liftOver <1st to 2nd genome liftOver file, unzipped>\n";
	print STDERR "\t\t-p1 <regions with specific phenotype to 1st strain (coord. based on genome1)>\n";
	print STDERR "\t\t-p2 <regions with specific phenotype in 2nd strain (coord. also based on genome1)>\n";
	print STDERR "\t\t\t-allowN (by default, homolgous regions with different number of Ns are discarded)\n";

	print STDERR "\n\tIf using VCF file and peak files with regions of interest:\n";
	print STDERR "\t\t-g1 <genome fasta file for genome>\n";
	print STDERR "\t\t-vcf <vcf file with variants>\n";
	print STDERR "\t\t\t(filters, column specification?)\n";
	print STDERR "\t\t-p1 <regions with specific phenotype for reference alleles>\n";
	print STDERR "\t\t-p2 <regions with specific phenotype in alternate allele>\n";

	print STDERR "\n\tIf using QTL file and peak files with regions of interest:\n";
	print STDERR "\t\t-g1 <genome fasta file for genome>\n";
	print STDERR "\t\t-qtl <QTL mapping file>\n";
	print STDERR "\t\t\t-pvalueThreshold <#> (default: $pvalueThreshold)\n";
	print STDERR "\t\t\t-absSlopeMin <#> (default: $absSlopeMin)\n";
	print STDERR "\t\t\t-mafThreshold <#> (default: $mafThreshold)\n";
	print STDERR "\t\t\t-flipScore (treat negative scores as 'up' in the Alternate allele)\n";
	#print STDERR "\t\t\t-useSignedPvalues (default: use correlation)\n";
	print STDERR "\t\t-peaks <QTL feature position - peak/BED file>\n";

	print STDERR "\n\t--- General Options ---\n";
	print STDERR "\n\tSpecifying distance size intervals to analyze:\n";
	print STDERR "\t\t-size <#,#> [#,#] [#,#] (region(s) to analyze)\n";
	print STDERR "\t\t\t-or-\n";
	print STDERR "\t\t-start <#> (default: $defStart)\n";
	print STDERR "\t\t-end <#> (default: $defEnd)\n";
	print STDERR "\t\t-step <#> (default: $defStep)\n";
	print STDERR "\t\t-window <#> (default: $defWindow)\n";

	print STDERR "\n\tMotif options:\n";
	print STDERR "\t\t-m <motif(s)> (required, analyzed sequences must contain this motif)\n";
	print STDERR "\t\t\t-adjustMotifThreshold <#> (change all motif logodds thresholds by this)\n";

	print STDERR "\n\tVarients to consider (vcf/qtl analysis only, won't work with liftOver analysis):\n";
	print STDERR "\t\t-snvOnly | -allVariants (do not consider indels, default, or consider everything [may not work will all options])\n";
	print STDERR "\t\t-indelsOnly (only consider indels, not snvs)\n";
	print STDERR "\t\t-subset (only consider these variants)\n";

	print STDERR "\n\tOutput:\n";
	print STDERR "\t\t-o <output directory> (required)\n";
	print STDERR "\n\tOther options:\n";
	print STDERR "\t\t-strand <.|+|-> (default, \"$strand\" for both)\n";
	print STDERR "\t\t-cpu <#cpus> (number of CPUs/threads for interval analysis, def: $maxCPUs)\n";
	print STDERR "\t\t-p <#cpus> (number of CPUs per interval i.e. number of motifs to analyze in parallel, def: $maxCPUs2)\n";
	print STDERR "\t\t-pp <#cpus> (number of CPUs for sequence generation, liftover, def: $maxCPUs3)\n";
	print STDERR "\t\t-gwas <gwasCatalog.txt> (will annotate variant results associated with GWAS hits)\n";
	print STDERR "\t\t-mask (use repeatmasked sequences for analysis)\n";
	print STDERR "\t\t-randomize (Randomize variant assignments to features)\n";
	#print STDERR "\t\t-zscore <#> (output zscore, where # is number of randomized motifs to check e.g. 10)\n";
	print STDERR "\t\t\t-motifOrder <#> (when randomizing motif, 0(def)=swap columsn, 1=try to preserve dinuc, 2=trinuc)\n";
	print STDERR "\t\t-markov <#> (Markov order for randomization, add an i for position independent model e.g. \"-markov 2i\", def: $markovOrder)\n";
	print STDERR "\t\t-rps <#> (randomizations per sequence, default: $randomizationsPerSequence)\n";
	print STDERR "\t\t-loadExistingRandomizations <filename> (skip model generation and randomization)\n";

	#print STDERR "\n\tConfiguring maggie:\n";
	#print STDERR "\t\t-exe <maggie executable> (def: $maggieExe)\n";

	print STDERR "\n";
	exit;
}
my $peaks1 = "";
my $peaks2 = "";
my $genome1= "";
my $genome2= "";
my $peakPosFile = "";
my $vcfFile = "";
my $qtlFile = "";
my $snpFile = "";
my $size = '';
my $motifFile = '';
my $gwasFile = '';
my @sizes = ();
my $headerInfo = '';
my $liftOver = "";
my $flipScore = 0;


for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-p1') {
		$peaks1 = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-p2') {
		$peaks2 = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-peaks') {
		$peakPosFile = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-strand') {
		$strand = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-randomize') {
		$randomizeVariantAssignments = 1;
	} elsif ( $ARGV[$i] eq '-g1' || $ARGV[$i] eq '-g' || $ARGV[$i] eq '-genome') {
		$genome1 = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-g2') {
		$genome2 = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-allowN') {
		$sameNflag = 0;
	} elsif ( $ARGV[$i] eq '-markov') {
		$markovOrder = $ARGV[++$i];
		$markovOrder = '' if ($markovOrder < 0);
	} elsif ( $ARGV[$i] eq '-subset') {
		$subsetFile = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-rps') {
		$randomizationsPerSequence = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-loadExistingRandomizations') {
		$existingRandomizationFile = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-qtl') {
		$qtlFile = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-mask') {
		$mask = " -mask ";
	} elsif ( $ARGV[$i] eq '-gwas') {
		$gwasFile = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-absSlopeMin') {
		$absSlopeMin = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-vcf') {
		$vcfFile = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-zscore') {
		$zscore =  $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-motifOrder') {
		$motifRandOrder = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-snp') {
		$snpFile = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-flipScore') {
		$flipScore = 1;
	} elsif ( $ARGV[$i] eq '-liftOver') {
		$liftOver = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-adjustMotifThreshold') {
		$adjustMotifThreshold = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-p') {
		$maxCPUs2 = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-pp') {
		$maxCPUs3 = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-size') {
		print STDERR "\tSize Ranges to Analyze:\n";
		my $bail =0;
		while ($ARGV[++$i] !~ /^\-[a-zA-Z]/) {
			push(@sizes, $ARGV[$i]);
			print STDERR "\t\t$ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
	    }
		last if ($bail == 1);
		$i--;
	} elsif ( $ARGV[$i] eq '-exe') {
		$maggieExe = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-o') {
		$outputDir = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-m') {
		$motifFile = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-allVariants') {
		$onlyIncludeSNV = 0;
		$onlyIncludeIndels = 0;
	} elsif ( $ARGV[$i] eq '-snvOnly') {
		$onlyIncludeSNV = 1;
	} elsif ( $ARGV[$i] eq '-indelsOnly') {
		$onlyIncludeIndels = 1;
	} elsif ( $ARGV[$i] eq '-pvalueThreshold') {
		$pvalueThreshold = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-mafThreshold') {
		$mafThreshold = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-start') {
		$defStart = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-end') {
		$defEnd = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-step') {
		$defStep = $ARGV[++$i];
	} elsif ( $ARGV[$i] eq '-window') {
		$defWindow = $ARGV[++$i];
	} else {
		print STDERR "!!! What is \" $ARGV[$i] \"??\n";
		printCMD();
	}
}

if ($motifFile eq '') {
	print STDERR "!!! Error - must specify a motif file with -m\n";
	printCMD();
}
if ($genome1 eq '') {
	print STDERR "!!! At minimum genome1 (-g1) must be specified\n";
	printCMD();
}
if ($outputDir eq '') {
	print STDERR "!!! Please specify and output directory: -o <dir> required\n";
	printCMD();
}

if (length(@sizes < 1)) {
	my $halfWindowD = ceil($defWindow/2);
	my $halfWindowU = floor($defWindow/2);
	for (my $i=$defStart;$i<=$defEnd;$i+=$defStep) {
		my $s = $i-$halfWindowD;
		my $e = $i+$halfWindowU;
		my $x = "$s,$e";
		push(@sizes, $x);
	}	
}

print STDERR "\tWill analyze the following size ranges:\n";
my $minStart =0;
my $maxEnd = 0;
foreach(@sizes) {
	print STDERR "\t\t$_\n";
	my @x = split /\,/;
	for (my $i=0;$i<@x;$i++) {
		if ($x[$i] < $minStart) {
			$minStart = $x[$i];
		}
		if ($x[$i] > $maxEnd) {
			$maxEnd = $x[$i];
		}
	}
}
$minStart -= $maxIndelLength;
$maxEnd += $maxIndelLength;
print STDERR "\tExtremes: $minStart,$maxEnd\n\n";


my $mode = "";

if (($peaks1 ne '' || $peaks2 ne '') && $genome1 ne '' && $genome2 ne '' && $liftOver ne '') {
	print STDERR "\n\tWill analyze motif mutations between different genomes using liftOver\n";
	print STDERR "\tThis will for the program to assume variants should have no impact on phenotype\n";
	$mode = "liftOver";
	$markovOrder = '';
	$headerInfo  = "\tFeatureID\tchr\tstart\tend\tstrand\tscore";

} elsif ($genome1 ne '' && $vcfFile ne '' && ($peaks1 ne '' || $peaks2 ne '')) {
	print STDERR "\n\tWill analyze motif mutations based on a VCF file\n";
	$mode = "vcf";

} elsif ($genome1 ne '' && $qtlFile ne '' && $peakPosFile ne '') {
	print STDERR "\n\tWill analyze QTLs\n";
	$mode = "qtl";
	$headerInfo  = "\tQtlID\tVariantID\tFeatureID\tchr\tstart\tend\tstrand\tscore\tpvalue\tslope";

} else {
	print STDERR "!!! Did not detect the right combination of input files:\n";
	printCMD();
}

my $gwas = '';
if ($gwasFile ne '') {
	print STDERR "\tReading GWAS information\n";
	$gwas = readGWAS($gwasFile);
}


`mkdir -p "$outputDir"`;
open LOG, ">$outputDir/log.txt";
print LOG "cmd: $cmd\n";
if ($zscore > 0) {
	`mkdir -p "$outputDir/dataByTF/"`;
}

my $rand = rand();
my $tmpFile = "$outputDir/" . $rand . ".tmp";
my $tmpFile2 = "$outputDir/" . $rand . ".2.tmp";
my $tmpFile3 = "$outputDir/" . $rand . ".3.tmp";
my $tmpFile4 = "$outputDir/" . $rand . ".4.tmp";
my $tmpFile5 = "$outputDir/" . $rand . ".5.tmp";
my $tmpFile6 = "$outputDir/" . $rand . ".6.tmp";


my %dataZ = ();
my $seqs = \%dataZ;
my %qtlDataZ = ();
my $qtlData = \%qtlDataZ;

my $flip = 0;
if ($mode eq 'liftOver') {
	$flip = 0;
	$flip = 1 if ($flipScore);
	unless ($scoreTestMode eq "+") {
		if ($peaks1 ne '') {
			getSeqLiftOver($seqs,$peaks1,$flip,$minStart,$maxEnd);
		}
	}
	$flip = 1;
	$flip = 0 if ($flipScore);
	unless ($scoreTestMode eq '-') {
		if ($peaks2 ne '') {
			getSeqLiftOver($seqs,$peaks2,$flip,$minStart,$maxEnd);
		}
	}
	printSequences("$outputDir/seqs.tsv", $seqs);
} elsif ($mode eq 'qtl') {
	$qtlData = readQTLFile($qtlFile,$peakPosFile,$minStart,$maxEnd);
	getSeqQTL($seqs,$qtlData,$minStart,$maxEnd,$subsetFile);

	
} elsif ($mode eq 'vcf') {

	my $vcfData = readVCFFile($vcfFile);

	$flip = 0;
	unless ($scoreTestMode eq "+") {
		if ($peaks1 ne '') {
			assignVariantsToPeaks($qtlData,$vcfData,$peaks1,$flip,$minStart,$maxEnd);
		}
	}
	$flip = 1;
	unless ($scoreTestMode eq '-') {
		if ($peaks2 ne '') {
			assignVariantsToPeaks($qtlData,$vcfData,$peaks2,$flip,$minStart,$maxEnd);
		}
	}
	exportQTLData($qtlData);

	my $qfile = "$outputDir/qtls.tsv";
	my $pfile = "$outputDir/qtlFeatures.tsv";
	$qtlData = readQTLFile($qfile,$pfile,$minStart,$maxEnd);
	getSeqQTL($seqs,$qtlData,$minStart,$maxEnd,$subsetFile);
} else {
	print STDERR "$mode mode not supported just yet... (or was mistyped or made up)\n";
	exit;
}

my $motifs = '';
%results = ();
my $model = '';

$motifs = HomerConfig::readMotifFile($motifFile);

my $numSizes = scalar(@sizes);
my $cpus = 0;
my %cleanNames = ();

my @motifs = sort {$b->{'len'} <=> $a->{'len'}} @$motifs;
my $numMotifs = scalar(@motifs);

my $zz = 0;
foreach(@motifs) {
	$zz++;

	#first get a name that can be written as part of a filename
	my $motif = $_;
	my $name = $motif->{'name'};
	my $clean = cleanMotifName($name);
	if (exists($cleanNames{$clean})) {
		my $x = 2;
		my $n = $clean . "_" . $x;
		while (exists($cleanNames{$n})) {
			$x++;
			$n = $clean . "_" . $x;
		}
		$clean = $n;
	}
	$motif->{'cleanName'} = $clean;
	$cleanNames{$clean}=1;


	my $dir = "$outputDir/motif_$clean/";

	$results{$clean} = $dir;

	print STDERR "\tAnalyzing motif $clean ($zz of $numMotifs total motifs)\n";
	
	$cpus++;
	my $pid = fork();
	if ($pid == 0) {
		#print STDERR "\tAnalyzing motif ($clean)\n";
		analyzeMotif($seqs,$motif,\@sizes,$minStart,$dir);
		exit;
	}
	if ($cpus >= $maxCPUs) {
		my $id = wait();
		$cpus--;
	}
}
my $id = 0;
while ($id >= 0) {
	$id = wait();
}
print STDERR "\tFinished checking all motifs\n";


### ......



open LOGP, ">$outputDir/logp.table.txt";
open LOGP0, ">$outputDir/logp0.table.txt";
open FDR, ">$outputDir/fdr.table.txt";
open FDR0, ">$outputDir/fdr0.table.txt";
open SCORE, ">$outputDir/scoreDiff.table.txt";
print LOGP "Motif";
print LOGP0 "Motif";
print FDR "Motif";
print FDR0 "Motif";
print SCORE "Motif";
for (my $i=0;$i<@sizes;$i++) {
	print LOGP "\t$sizes[$i]";
	print LOGP0 "\t$sizes[$i]";
	print FDR "\t$sizes[$i]";
	print FDR0 "\t$sizes[$i]";
	print SCORE "\t$sizes[$i]";
}
print LOGP "\n";
print LOGP0 "\n";
print FDR "\n";
print FDR0 "\n";
print SCORE "\n";
print LOGP "midpoint";
print LOGP0 "midpoint";
print FDR "midpoint";
print FDR0 "midpoint";
print SCORE "midpoint";
for (my $i=0;$i<@sizes;$i++) {
	my @z = split /\,/,$sizes[$i];
	my $m = floor(($z[0]+$z[1])/2);
	print LOGP "\t$m";
	print LOGP0 "\t$m";
	print FDR "\t$m";
	print FDR0 "\t$m";
	print SCORE "\t$m";
}
print LOGP "\n";
print LOGP0 "\n";
print FDR "\n";
print FDR0 "\n";
print SCORE "\n";

open STATS, ">$outputDir/stats.tfs.txt";
print STATS "Motif\tRegion\tDistance to TSS\tTotal sequences\tTotal with differential motif"
			. "\tAvg motif score change\tStdev Motif score change"
			. "\tAvg random motif score change\tStdev random motif score change"
			. "\tZscore (vs. random)\tRankSum p-value (vs. random)\tRankSum signed logP (vs. random)"
			. "\tZscore (vs. zero)\tRankSum p-value (vs. zero)\tRankSum signed logP (vs. zero)\n";

open SITES, ">$outputDir/mutatedSites.txt";
print SITES "RegionName\tFeatureDistance\tRegionStart\tRegionEnd\tMotif\tRelativePositionToTSS\tRelativeStrandToTSS"
		. "\tControlAlleleLogOdds\tActiveAlleleLogOdds\tControlSite\tActiveSite\tDiffLogOdds\tGained/Lost\tPhenotype"
		. "\tchr(motif)\tstart(motif)\tend(motif)\tstrand(motif)"
		. "$headerInfo\tGWASinfo\n";

my %printedSites = ();

my %randResults = ();
my %TFfiles = ();

my %pvalues = ();
my %pvalues0 = ();


@motifs = sort {$a->{'cleanName'} cmp $b->{'cleanName'}} @motifs;
for (my $j=0;$j<@motifs;$j++) {
	my $mname = $motifs[$j]->{'name'};
	my $cleanName = $motifs[$j]->{'cleanName'};

	print LOGP "$cleanName";
	print LOGP0 "$cleanName";
	print SCORE "$cleanName";

	my $dir = $results{$cleanName};
	my $statFile = "$dir/$cleanName.stats.txt";
	my $siteFile = "$dir/$cleanName.mutatedSites.txt";


	open IN, $statFile or print STDERR "\tCould not open $statFile\n";
	my %statLines = ();
	while (<IN>) {
		print STATS $_;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $size = $line[1];
		$statLines{$size} = \@line;
	}
	close IN;

	for (my $i=0;$i<@sizes;$i++) {
		my $size = $sizes[$i];


		my $deltaAvg = '';
		my $rdeltaAvg = '';
		my $kwpvalue = '';
		my $logKWpvalue = "";
		my $kwpvalue0 = '';
		my $logKWpvalue0 = "";
		if (exists($statLines{$size})) {
			$deltaAvg = $statLines{$size}->[5];
			$rdeltaAvg = $statLines{$size}->[7];
			$kwpvalue = $statLines{$size}->[10];
			$logKWpvalue = $statLines{$size}->[11];
			$kwpvalue0 = $statLines{$size}->[13];
			$logKWpvalue0 = $statLines{$size}->[14];
		}

		my $pid = $cleanName . "_" . $size;
		$pvalues{$pid} = {p=>$kwpvalue};
		$pvalues0{$pid} = {p=>$kwpvalue0};

		my $deltaDiff = 0;
		if ($deltaAvg ne '') {
			$deltaDiff = $deltaAvg - $rdeltaAvg;
		}
		#my $logKWpvalue = log($kwpvalue)/log(10.0);
		#if ($deltaDiff < 0) {
		#	$logKWpvalue *= -1;
		#}


		open IN, $siteFile or print STDERR "Could not open $siteFile";
		while (<IN>) {
			chomp;
			my $og = $_;
			my @line = split /\t/;
			my $qtlid = $line[0];
			my $mname = $line[4];
			my $mstrand = $line[6];
			my $mstatus = $line[10];
			my $pheno = $line[11];
			my $uniqID = $qtlid . "|" . $mname . "|" . $mstrand . "|" . $mstatus . "|" . $pheno;
		
			if (exists($printedSites{$uniqID})) {
				next;
			} else {
				$printedSites{$uniqID}=1;
			}

			my $vid = $line[19];
			my $gwasinfo = '';
			if ($vid =~ /^(.*?)\_(\d+?)\_/) {
				my $chr = $1;
				my $p =  $2;
				my $gwasid = $chr . "-" . $p;
				if (exists($gwas->{$gwasid})) {
					$gwasinfo = $gwas->{$gwasid};
				}
			}
	
			my $tssDist=$line[5];
			if (exists($qtlData->{$qtlid})) { 
				$tssDist = $qtlData->{$qtlid}->{'tssDist'};
			} else {
				#print STDERR "!!! Couldn't find $qtlid\n";
			}

			print SITES "$qtlid";
			for (my $k=1;$k<@line;$k++) {
				if ($k==5) {
					print SITES "\t$tssDist";
				} else {
					print SITES "\t$line[$k]";
				}
			}
			print SITES "\t$gwasinfo\n";
		}
		close IN;


		my $numSeq = 0;
		my $givenName = '';
		my $medianPvalue = 0;
		my $medianScoreDiff = 0;
		my $avgScoreDiff = 0;
		my $posMut = 0;
		my $negMut = 0;

		#my $score = $avgScoreDiff;
		my $score = $medianPvalue;

		my $randAvg = 0;
		my $randStd = 0;
		my $zscoreValue = 0;
		if ($zscore > 0) {
			if (!exists($randResults{$mname}->{$size})) {
				print STDERR "!!! Error: Could not find randomized motif results for $mname (size=$size)\n";
				next;
			}
			($randAvg, $randStd) = Statistics::avevar($randResults{$mname}->{$size});
			$randStd = sqrt($randStd);
			if ($randStd < 1e-99) {
				$zscoreValue = 0;
			} else {
				$zscoreValue = ($score - $randAvg)/$randStd;
			}
		}

		$medianPvalue = 0 if ($medianPvalue eq 'na');


		if (1) {
			$logKWpvalue = 'NA' if ($logKWpvalue eq '');
			$logKWpvalue0 = 'NA' if ($logKWpvalue0 eq '');
			print LOGP "\t$logKWpvalue";
			print LOGP0 "\t$logKWpvalue0";
			print SCORE "\t$deltaDiff";
		} else {
			print LOGP "\t$medianPvalue";
			print LOGP0 "\t$medianPvalue";
			print SCORE "\t$avgScoreDiff";
		}

		#`rm -r "$dir"`;
	}
	print LOGP "\n";
	print LOGP0 "\n";
	print SCORE "\n";
}
close LOGP;
close LOGP0;
close SCORE;
close LOG;
close STATS;
close SITES;

#benjamini p-value correction
if (1) { 
	my @tests = keys %pvalues;
	my @p = ();
	my @p0 = ();
	foreach(@tests) {
		push(@p, $pvalues{$_}->{'p'});
		push(@p0, $pvalues0{$_}->{'p'});
	}		
	#my $correctedPvalues = Statistics::benjaminiFDR(\@p);
	#my $correctedPvalues0 = Statistics::benjaminiFDR(\@p0);
	my $correctedPvalues = UseRtoCorrectPvalues(\@p);
	my $correctedPvalues0 = UseRtoCorrectPvalues(\@p0);

	for (my $i=0;$i<@tests;$i++) {
		$pvalues{$tests[$i]}->{'q'} = $correctedPvalues->[$i];
		$pvalues0{$tests[$i]}->{'q'} = $correctedPvalues0->[$i];
	}
	open IN, "$outputDir/stats.tfs.txt";
	open STATS, ">$outputDir/stats.tfs.txt.new";

	$c=0;
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my $og = $_;
		my @line = split /\t/;

		if ($c ==1) {
			print STATS "$line[0]";
			for (my $i=1;$i<12;$i++) {
				print STATS "\t$line[$i]";
			}
			print STATS "\tq-value (vs. random)\tsigned log q-value (vs. random)";
			for (my $i=12;$i<@line;$i++) {
				print STATS "\t$line[$i]";
			}
			print STATS "\tq-value (vs. zero)\tsigned log q-value (vs. zero)";
			print STATS "\n";
			next;
		}
		my $pid = $line[0] . "_" . $line[1];
		my $str = "\tna\tna";
		my $str0 = "\tna\tna";
		if (!exists($pvalues{$pid})) {
			my @a = keys %pvalues;
			print STDERR "!!! Can't find entry in pvalue table! ($pid) (line=$og) (pvalues=@a)\n";
		} else {
			my $qvalue = $pvalues{$pid}->{'q'};
			my $logq= log($qvalue+$epsPvalue)/log(10.0);
			my $logp = $line[11];
			if ($logp > 0) {
				$logq = -1*$logq;
			}
			$pvalues{$pid}->{'logq'} = $logq;
			$str = "\t$qvalue\t$logq";
			my $qvalue0 = $pvalues0{$pid}->{'q'};
			my $logq0= log($qvalue0+$epsPvalue)/log(10.0);
			my $logp0 = $line[14];
			if ($logp0 > 0) {
				$logq0 = -1*$logq0;
			}
			$pvalues{$pid}->{'logq0'} = $logq0;
			$str2 = "\t$qvalue0\t$logq0";
		}
		print STATS "$line[0]";
		for (my $i=1;$i<12;$i++) {
			print STATS "\t$line[$i]";
		}
		print STATS $str;
		for (my $i=12;$i<@line;$i++) {
			print STATS "\t$line[$i]";
		}
		print STATS $str2;
		print STATS "\n";
	}
	close IN;
	close STATS;

	for (my $j=0;$j<@motifs;$j++) {
		my $mname = $motifs[$j]->{'cleanName'};
		next if ($mname =~ /^Rand_(\d+)_(.+)/);
		print FDR "$mname";
		print FDR0 "$mname";
		for (my $i=0;$i<@sizes;$i++) {
			my $size = $sizes[$i];
			my $pid = $mname . "_" . $size;
			my $logq = 0;
			my $logq0 = 0;
			if (!exists($pvalues{$pid})) {
				print STDERR "!!! Could not find $pid in fdr results\n";
			} else {
				if (exists($pvalues{$pid}->{'logq'})) {
					$logq = $pvalues{$pid}->{'logq'};
				}
				if (exists($pvalues{$pid}->{'logq0'})) {
					$logq0 = $pvalues{$pid}->{'logq0'};
				}
			}
			$logq = "NA" if ($logq eq ''); 
			$logq0 = "NA" if ($logq0 eq ''); 
			print FDR "\t$logq";
			print FDR0 "\t$logq0";
		}
		print FDR "\n";
		print FDR0 "\n";
	}
	`mv "$outputDir/stats.tfs.txt.new" "$outputDir/stats.tfs.txt"`;

}
close FDR;
close FDR0;

print STDERR "\n\tOutput will be found in the directory: $outputDir\n\n";

exit;

sub generateSequenceModel {
	my ($seqs, $markovOrder) = @_;

	print STDERR "\tCreating background sequence model...\n";
	#extract sequences to build model
	my %cleanSeqs = ();
	my $numSeqs = 0;
	my @seqIDs = keys %$seqs;
	foreach(@seqIDs) {
		$numSeqs++;
		my $id = $_;
		my $s1 = $seqs->{$id}->{'s1'};
		my $s2 = $seqs->{$id}->{'s2'};
		$cleanSeqs{$id . "-s1"} = $s1;
		$cleanSeqs{$id . "-s2"} = $s2;
	}
	my $model = '';
	my $posIndepFlag = 0;
	if ($markovOrder =~ s/i$//) {
		$posIndepFlag = 1;
		print STDERR "\t\tBuilding a position independent model (order=$markovOrder)\n";
		$model = HomerSeqBias::buildPosIndepMarkovModel(\%cleanSeqs, $markovOrder);
	} else {
		print STDERR "\t\tBuilding a position dependent model (order=$markovOrder)\n";
		$model = HomerSeqBias::buildPosDepMarkovModel(\%cleanSeqs, $markovOrder);
	}
	return $model;
}



sub generateRandomizedSequence {
	my ($seqs, $markovOrder, $randomizationsPerSequence) = @_;

	my $model = generateSequenceModel($seqs,$markovOrder);

	print STDERR "\tGenerating Randomized Control Sequences for Analysis... (n=$randomizationsPerSequence per input sequence)\n";

	my @files = ();
	my $filePrefix = $tmpFile . "_randSeq";

	my $seqsPerCPU = ceil($numSeqs/$maxCPUs3);
	my $index =0;
	my $cpus = 0;
	for (my $i=0;$i<$maxCPUs3;$i++) {
		my @ids=  ();
		for (;$index<(1+$i)*$seqsPerCPU;$index++) {
			last if ($index >= $numSeqs);
			#push(@ids, $seqIDs[$index]);
		}

		my $fname = $filePrefix . "$i";
		push(@files, $fname);
		$cpus++;
		my $pid = fork();

		if ($pid == 0) {

			#for each sequence, generate random sequences harboring variants from model
			open OUT, ">$fname";
			foreach(@ids) {
				my $id = $_;
		
				my $og1 = $seqs->{$id}->{'s1'};
				my $og2 = $seqs->{$id}->{'s2'};
		
				my $refSeq = $seqs->{$id}->{'refForRand'};
				my $altSeq = $seqs->{$id}->{'altForRand'};
				my $refPos = $seqs->{$id}->{'refAbsPos'};
				my $altPos = $seqs->{$id}->{'altAbsPos'};
				if ($seqs->{$id}->{'refFlip'}==1) {
					$refSeq = $seqs->{$id}->{'altForRand'};
					$altSeq = $seqs->{$id}->{'refForRand'};
					$refPos = $seqs->{$id}->{'altAbsPos'};
					$altPos = $seqs->{$id}->{'refAbsPos'};
				}

				my $L = length($og1);
		
				for (my $j=0;$j<$randomizationsPerSequence; $j++) {
					my $s1 = '';
					my $s2 = '';
					my $posIndepFlag=0;
					if ($posIndepFlag) {
						($s1,$s2) = HomerSeqBias::generateRandSequenceWithVariant($model,$refSeq,$refPos,$altSeq,$altPos,$L);
					} else {
						($s1,$s2) = HomerSeqBias::generateRandSequenceWithVariant($model,$refSeq,$refPos,$altSeq,$altPos);
					}
					#print STDERR "$refSeq\t$refPos\t$altSeq\t$altPos\n$s1\n$s2\n\n";
					next if ($s1 eq '' || $s2 eq '');
					print OUT "$id" . "_r$j" . "_1\t$s1\n";				
					print OUT "$id" . "_r$j" . "_2\t$s2\n";				
				}
			}
			close OUT;
			exit;
		}
		if ($cpus >= $maxCPUs3) {
			my $id = wait();
			$cpus--;
		}
	}
	my $pid = 0;
	while ($pid >=0) {
		$pid = wait();
	}

	foreach(@files) {
		my $f = $_;
		my $previousSeq = '';
		my $previousID = '';
		my $previousR = '';
		my $previousV = '';
		open IN, $f;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			my $seq = $line[1];
			my $id = '';
			my $r = '';
			my $v = '';
			if ($line[0] =~ /(.+?)_(r\d+)_([12])$/) {
				$id = $1;
				$r = $2;
				$v = $3;
				if ($previousSeq ne '' && $previousID eq $id && $previousR eq $r) {
					push(@{$seqs->{$id}->{'s1rand'}},$previousSeq);
					push(@{$seqs->{$id}->{'s2rand'}},$seq);
				}
				$previousID = $id;
				$previousR = $r;
				$previousV = $v;
			} else {
				print STDERR "\tWarning: couldn't parse id $line[0]\n";
			}

			#print STDERR "$id\t$refSeq $refPos $altSeq $altPos\n$og1\n$og2\n$s1\n$s2\n";

			$previousSeq = $seq;
		}
		close IN;
		`rm $f`;
	}
	print STDERR "\tFinished creating random control sequences\n";
}

sub cleanMotifName {
	my ($name) = @_;
	my $mname = $name;
	$mname =~ s/\/.*//;
	$mname =~ s/\,//g;
	$mname =~ s/\(/_/g;
	$mname =~ s/\)/_/g;
	$mname =~ s/\s/_/g;
	return $mname;
}



# this position could overlap a motif with different offsets.
# Find motif offset whose alt or ref produce the highest motif score
sub getBestMotifMatch {
	my ($scoresPos,$scoresNeg,$pdiffs,$ndiffs,$mlen,$pos,$strand) = @_;
	my $bestMotifPos = -1;
	my $bestMotifScore = -1e9;
	my $bestMotifStrand = '';
	my $bestMotifAlt = -1;
	my $bestMotifDiff = '';
	for (my $j=$pos-$mlen+1;$j<=$pos;$j++) {
		#skip if out of range
		next if ($j<0);
		next if ($j>=@$scoresPos);
		#skip if sequence has an 'N' - would produce a score of 'na' in both Pos and Neg strands
		next if ($scoresPos->[$j] eq 'na');
		my $varOffsetInMotif = $pos-$j;
#print STDERR "$varOffsetInMotif $pos $j $scoresPos->[$j] $scoresNeg->[$j]\n";
		my $altPosScore = $scoresPos->[$j]+ $pdiffs->[$varOffsetInMotif];
		my $altNegScore = $scoresNeg->[$j]+ $ndiffs->[$varOffsetInMotif];
		if (($scoresPos->[$j] > $bestMotifScore) && ($strand ne '-')) {
			$bestMotifScore = $scoresPos->[$j];
			$bestMotifPos = $j;
			$bestMotifStrand = "+";
			$bestMotifAlt = 0;
			$bestMotifDiff = $pdiffs->[$varOffsetInMotif];
		}
		if (($altPosScore > $bestMotifScore) && ($strand ne '-')) {
			$bestMotifScore = $altPosScore;
			$bestMotifPos = $j;
			$bestMotifStrand = "+";
			$bestMotifAlt = 1;
			$bestMotifDiff = $pdiffs->[$varOffsetInMotif];
		}
		if ($scoresNeg->[$j] > $bestMotifScore && ($strand ne '+')) {
			$bestMotifScore = $scoresNeg->[$j];
			$bestMotifPos = $j;
			$bestMotifStrand = "-";
			$bestMotifAlt = 0;
			$bestMotifDiff = $ndiffs->[$varOffsetInMotif];
		}
		if ($altNegScore > $bestMotifScore && ($strand ne '+')) {
			$bestMotifScore = $altNegScore;
			$bestMotifPos = $j;
			$bestMotifStrand = "-";
			$bestMotifAlt = 1;
			$bestMotifDiff = $ndiffs->[$varOffsetInMotif];
		}
	}
	return ($bestMotifPos, $bestMotifScore, $bestMotifStrand,$bestMotifAlt,$bestMotifDiff);
}

# seqs - has datastructure of all of the sequences and their annotation
# motif - motif (data structure) to analyze
# sizes - ref to array of size ranges to calculate enrichment on
# offset - reference position of the start of the sequences
# motifOutputDir - directory to dump results - will be read in by the main program (this will be exectuted on a fork)
sub analyzeMotif {
	my ($seqs,$motif,$sizes,$offset,$motifOutputDir) = @_;

	# create output directory
	`mkdir -p $motifOutputDir/`;

	# setup motif for analysis, get key information, open reporting files
	my $mname = $motif->{'name'};
	my $mlen = $motif->{'len'};
	my $cleanName = $motif->{'cleanName'};
	my $mthreshold = $motif->{'threshold'};
	my $detectThreshold = $mthreshold;
	if ($adjustMotifThreshold ne '') {
		$detectThreshold += $adjustMotifThreshold;
	}
	if (!defined($motif->{'lom'})) {
		HomerConfig::addLogOddsMatrix($motif);
	}
	my $sitesFile = "$motifOutputDir/$cleanName.mutatedSites.txt";	
	open SITES, ">$sitesFile" or die "Could not open $sitesFile for writing\n";
	my $psitesFile = "$motifOutputDir/$cleanName.potentialSites.txt";	
	open PSITES, ">$psitesFile" or die "Could not open $psitesFile for writing\n";
	my $statsFile = "$motifOutputDir/$cleanName.stats.txt";	
	open STATS, ">$statsFile" or die "Could not open $statsFile for writing\n";

	#setup reporting variables, datastructures
	my %motifVarTypes = ();
	my @motifVarPos = ();
	my @motifVarNeg = ();
	for (my $i=0;$i<$mlen;$i++) {
		my %a = ();
		my %b = ();
		push(@motifVarPos, \%a);
		push(@motifVarNeg, \%b);
	}
	my $totalSeqAnalyzed =0;
	my $totalSeqWithMotif = 0;
	my $totalSeqDiffMotif = 0;
	my $totalPotentialVarMotifs = 0;
	my %deltas = ();

	my %siteHist = ();
	my %varHist = ();
	my %possibleVars = ();
	my $seqdiff=0;
	my $seqmotif=0;	
	my $seqmotifdiff=0;

	# first lets collect the type of variants and their positional frequencies
	# the varTypes is indexed by the inactiveVar > activeVAR, not ref > alt
	my %varTypes = ();
	my $minPos=1e10;
	my $maxPos=-1e10;
	my $maxLen=-1;
	my $tSeqs = 0;
	foreach(keys %$seqs) {
		my $id = $_;
		my $varPos = $seqs->{$id}->{'refAbsPos'};
		$minPos = $varPos if ($varPos < $minPos);
		$maxPos = $varPos if ($varPos > $maxPos);
		$maxLen = length($seqs->{$id}->{'seq'});
		my $refVar = $seqs->{$id}->{'refForRand'};
		my $altVar = $seqs->{$id}->{'altForRand'};
		my $refActive = $seqs->{$id}->{'refActive'};
		my $var = $refVar . ">" . $altVar;
		my $flippedVar = $altVar . ">" . $refVar;
		if ($refActive) {
			$var = $altVar . ">" . $refVar;
			$flippedVar = $refVar . ">" . $altVar;
			my $x = $refVar;
			$refVar = $altVar;
			$altVar = $x;
		}
		if (!exists($varTypes{$var})) {
			my %a = ();
			my %b = ();
			my %c = ();
			my %d = ();
			my %e = ();
			my %f = ();
			my $x = {p=>\%a,d=>\%d,rap=>\%b,aap=>\%c,r=>$refVar,a=>$altVar,f=>$flippedVar,tp=>\%e,td=>\%f};
			$varTypes{$var} = $x;
		}
		$varTypes{$var}->{'p'}->{$varPos}++;
		if ($refActive) {
			$varTypes{$var}->{'rap'}->{$varPos}++;
		} else {
			$varTypes{$var}->{'aap'}->{$varPos}++;
		}
		$tSeqs++;
	}
	my @varTypes = sort {$a cmp $b} keys %varTypes;

	my %motifHist = ();
	for (my $i=0;$i<$maxLen;$i++) {
		my $a = {'+'=>0,'-'=>0};
		my $b = {'+'=>0,'-'=>0};
		my $c = {'+'=>0,'-'=>0};
		my $d = {'+'=>0,'-'=>0};
		my $e = {'+'=>0,'-'=>0};
		$motifHist{$i}={all=>$a,strainSpec=>$b,gained=>$c,lost=>$d,potential=>$e};
	}

my $totalPotential=0;
my $totalFound=0;


	# iterate through each sequence, scan for motif scores, consider variant impact, and tabulate changes in
	# motif scores, motif positions, and the variants that would cause it.
	foreach(keys %$seqs) {
		my $id = $_;
		my $sname = $id;
		my $seq = $seqs->{$id}->{'seq'};
		my $seqLen = length($seq);
		my $tssDist = $seqs->{$id}->{'tssDist'};
		my $varPos = $seqs->{$id}->{'refAbsPos'};
		my $refVar = $seqs->{$id}->{'refForRand'};
		my $altVar = $seqs->{$id}->{'altForRand'};
		my $refActive = $seqs->{$id}->{'refActive'};
		my $var = $refVar . ">" . $altVar;
		if ($refActive) {
			$var = $altVar . ">" . $refVar;
		}

		$totalSeqAnalyzed++;

		# update status every 1000 sequences processed
		if ($totalSeqAnalyzed % 1000 == 0) {
			my $r = 0;
			$r = sprintf("%.2f",100*$totalSeqAnalyzed/$tSeqs) if ($tSeqs>0);
			print STDERR "\t\t$cleanName $totalSeqAnalyzed/$tSeqs ($r%)\n";
		}

		# scan original sequence for its match to the motif at each position - scoresPos/Neg are arrays of these values
		my ($scoresPos, $scoresNeg) = HomerConfig::scoreMotifAtAllPositions($motif, $seq);

		# quantify motif distribution along sequence
		for (my $i=0;$i<@$scoresPos;$i++) {
			$motifHist{$i}->{'all'}->{'+'}++ if ($scoresPos->[$i] ne 'na' && $scoresPos->[$i] > $detectThreshold && $strand ne '-');
			$motifHist{$i}->{'all'}->{'-'}++ if ($scoresNeg->[$i] ne 'na' && $scoresNeg->[$i] > $detectThreshold && $strand ne '+');
		}

		# if a variant is in the sequence, analyze it's impact
		if ($var ne '' && $varPos >= 0) {
			my $pdiffs = HomerConfig::getPossibleScoreChanges($motif, $refVar,$altVar, "+");
			my $ndiffs = HomerConfig::getPossibleScoreChanges($motif, $refVar,$altVar, "-");

			my ($bestMotifPos,$bestMotifScore,$bestMotifStrand,$bestMotifAlt,$bestMotifDiff) 
							= getBestMotifMatch($scoresPos,$scoresNeg,$pdiffs,$ndiffs,$mlen,$varPos,$strand);

			# if variant is on a motif or creates a motif above the threshold
			if ($bestMotifScore > $detectThreshold) {
				$totalSeqDiffMotif++;

				# record motifs location for QC later
				if ($bestMotifAlt==1 && ($bestMotifScore-$bestMotifDiff < $detectThreshold)) {
					$motifHist{$bestMotifPos}->{'strainSpec'}->{$bestMotifStrand}++;
					$motifHist{$bestMotifPos}->{'all'}->{$bestMotifStrand}++; #would have missed in our scan of ref sequence
				} elsif ($bestMotifAlt==0 && ($bestMotifScore+$bestMotifDiff < $detectThreshold)) {
					$motifHist{$bestMotifPos}->{'strainSpec'}->{$bestMotifStrand}++;
				}
					   
				# assign motif scores for the inactive and active allele
				# default: bestMotifAlt==0, refActive==0
				my $mscore1 = $bestMotifScore;
				my $mscore2 = $bestMotifScore + $bestMotifDiff;
				if ($bestMotifAlt==1 && $refActive==0) {
					$mscore1 = $bestMotifScore - $bestMotifDiff;
					$mscore2 = $bestMotifScore;
				} elsif ($bestMotifAlt == 0 && $refActive==1) {
					$mscore1 = $bestMotifScore + $bestMotifDiff;
					$mscore2 = $bestMotifScore;
				} elsif ($bestMotifAlt==1 && $refActive==1) {
					$mscore1 = $bestMotifScore; 
					$mscore2 = $bestMotifScore - $bestMotifDiff;
				}	   

				# record differential motif score associated with the variant	
				if ($refActive) {
					$bestMotifDiff *= -1;
				}
				if (!exists($deltas{$varPos})) {
					my @a = ();
					$deltas{$varPos}=\@a;
				}
				push(@{$deltas{$varPos}}, $bestMotifDiff);
	
				# report the site
				my $status = "same";
				if ($bestMotifDiff < 0) {
					$status = "lost";
				} elsif ($bestMotifDiff > 0) {
					$status = "gained";
				}

				if ($status eq 'lost' || $status eq 'gained') {
					$motifHist{$bestMotifPos}->{$status}->{$bestMotifStrand}++;
				}

				# report site location
				my $info = $seqs->{$id}->{'info'};
				my $pheno = $seqs->{$id}->{'phenotype'};
				my $siteSeqRef = substr($seq,$bestMotifPos,$mlen);
				my $varPosInMotif = $varPos-$bestMotifPos;
				my $siteSeqAlt = substr($seq,$bestMotifPos,$varPosInMotif) . $altVar . substr($seq,$varPos+1,$mlen-$varPosInMotif-1);
				my $siteSeq1 = $siteSeqRef;
				my $siteSeq2 = $siteSeqAlt;
				if ($refActive) {
					$siteSeq1 = $siteSeqAlt;
					$siteSeq2 = $siteSeqRef;
				}
				my $siteChr = $seqs->{$id}->{'c'};
				my $siteStart = $seqs->{$id}->{'tp'}+$offset+$bestMotifPos;
				my $siteEnd = $siteStart + $mlen-1;
				my $siteStrand = $bestMotifStrand;
				if ($seqs->{$id}->{'td'} eq '-') {
					$siteEnd = $seqs->{$id}->{'tp'}-$offset-$bestMotifPos;
					$siteStart = $siteEnd-$mlen+1;
					if ($siteStrand eq '+') {
						$siteStrand = "-";
					} else {
						$siteStrand = "+";
					}
				}
				my $mpos = $bestMotifPos+$offset;
				print SITES "$id\t$tssDist\tna\tna\t$cleanName\t$mpos\t$bestMotifStrand\t$mscore1"
									. "\t$mscore2\t$siteSeq1\t$siteSeq2\t$bestMotifDiff"
									. "\t$status\t$pheno\t$siteChr\t$siteStart\t$siteEnd\t$siteStrand\t$info\n";

				#if ($diff > 0) {
				#	$motifVarPos[$varPosInMotif]->{$curCtrlVar."->".$curActiveVar}++;
				#	$motifVarTypes{$curCtrlVar."->".$curActiveVar}++;
				#} else {
				#	$motifVarNeg[$varPosInMotif]->{$curActiveVar."->".$curCtrlVar}++;
				#	$motifVarTypes{$curActiveVar."->".$curCtrlVar}++;
				#}
			}
		}

		# assess the change in motif score for potential variants found in each position to model background
		my %tmpMotifHist = ();
		my $motifFound = 0;
		foreach(@varTypes) {
			my $var = $_;
			my $flippedVar = $varTypes{$var}->{'f'};
			my $ref = $varTypes{$var}->{'r'};
			my $alt = $varTypes{$var}->{'a'};
			my $pdiffs = HomerConfig::getPossibleScoreChanges($motif, $ref,$alt, "+");
			my $ndiffs = HomerConfig::getPossibleScoreChanges($motif, $ref,$alt, "-");
			for (my $i=0;$i<$seqLen;$i++) {

				# next make sure variant is consistent with sequence - i.e. A>T should only happen where there is an A
				my $s = substr($seq,$i,1);
				next if ($s ne $ref);

				$tmpMotifHist{$i} = {'+'=>0,'-'=>0,'p'=>0} if (!exists($tmpMotifHist{$i}));
				$tmpMotifHist{$i}->{'p'}++;

				# if reference is the active allele, flip the variant order
				my $eVar = $var;
				if ($refActive) {
					$eVar = $flippedVar;
				}

				if (!exists($varTypes{$eVar}->{'tp'}->{$i})) {
					$varTypes{$eVar}->{'tp'}->{$i}=0;
					$varTypes{$eVar}->{'td'}->{$i}=0;
				}
				$varTypes{$eVar}->{'tp'}->{$i}++;

				# can skip if there are no variants of this type in this position
				# next if (!exists($varTypes{$var}->{'p'}->{$i}));

				# this position could overlap a motif with different offsets.
				# Find motif offset whose alt or ref produce the highest motif score
				my ($bestMotifPos,$bestMotifScore,$bestMotifStrand,$bestMotifAlt,$bestMotifDiff) 
							= getBestMotifMatch($scoresPos,$scoresNeg,$pdiffs,$ndiffs,$mlen,$i,$strand);


				# check to make sure a best motif was assigned and that it exceeds the detection threshold
				next if ($bestMotifScore < $detectThreshold);

				$tmpMotifHist{$bestMotifPos} = {'+'=>0,'-'=>0,'p'=>0} if (!exists($tmpMotifHist{$bestMotifPos}));
				$tmpMotifHist{$bestMotifPos}->{$bestMotifStrand}++;

				if ($refActive) {
					$bestMotifDiff *= -1;
				}

				if (!exists($varTypes{$eVar}->{'d'}->{$i})) {
					my @a = ();
					$varTypes{$eVar}->{'d'}->{$i} = \@a;
				}
				push(@{$varTypes{$eVar}->{'d'}->{$i}},$bestMotifDiff);

				$varTypes{$eVar}->{'td'}->{$i}++;
				$motifFound = 1;

				my $siteChr = $seqs->{$id}->{'c'};
				my $siteStart = $seqs->{$id}->{'tp'}+$offset+$bestMotifPos;
				my $siteEnd = $siteStart + $mlen-1;
				my $siteStrand = $bestMotifStrand;
				if ($seqs->{$id}->{'td'} eq '-') {
					$siteEnd = $seqs->{$id}->{'tp'}-$offset-$bestMotifPos;
					$siteStart = $siteEnd-$mlen+1;
					if ($siteStrand eq '+') {
						$siteStrand = "-";
					} else {
						$siteStrand = "+";
					}
				}

				print PSITES "$id\t$i\t$eVar\t$siteChr\t$siteStart\t$siteEnd\t$siteStrand\t$bestMotifScore\t$bestMotifDiff\t$bestMotifAlt\t$eVar\n";
				$totalPotentialVarMotifs++;
			}
		}
		# this helps us to avoid overcounting sites where variants may help identify a motif
		foreach(keys %tmpMotifHist) {
			my $p = $_;
			$totalPotential++ if ($tmpMotifHist{$p}->{'p'} > 0);
			$totalFound++ if ($tmpMotifHist{$p}->{'+'} > 0 || $tmpMotifHist{$p}->{'-'} > 0);

			$motifHist{$p}->{'potential'}->{'+'}++ if ($tmpMotifHist{$p}->{'+'} > 0);
			$motifHist{$p}->{'potential'}->{'-'}++ if ($tmpMotifHist{$p}->{'-'} > 0);
		}

		if ($motifFound) {
			$totalSeqWithMotif++;
		}
	}
	close SITES;
	close PSITES;

	#print STDERR "totalFound/vsPotential: $totalFound / $totalPotential\n";

	open MDELTAS, ">$motifOutputDir/motif.deltas.tsv";
	print MDELTAS "Position\tTotal differential sites\tAvg. motif difference\tAvg. random motif difference\n";
	my %reportedPos = ();

	#print STDERR "\t\tTotalSeqAnalyzed=$totalSeqAnalyzed\tTotalSeqDiffMotif=$totalSeqDiffMotif\ttotalPotentialVarMotifs=$totalPotentialVarMotifs\n";

	# now that the sequences and motifs have been identified, cycle through the different intervals and evaluate
	# the differences in motif scores found
	for (my $i=0;$i<@$sizes;$i++) {
		my $size = $sizes->[$i];
		my @z = split/\,/,$size;
		my $sizeStart = $z[0];
		my $sizeEnd = $z[1];
		$tssDist = floor(($sizeStart+$sizeEnd)/2);
		my $absSizeStart = $sizeStart - $offset;
		my $absSizeEnd = $sizeEnd - $offset;
		my $intervalDist = floor(($sizeStart+$sizeEnd)/2);
		my $sizeLen = $sizeEnd-$sizeStart;

		my @deltas = ();
		my @rdeltas = ();
		my %rdeltaWeight = ();
		my $numDeltas = 0;

		for (my $j=$absSizeStart;$j<$absSizeEnd;$j++) {
			my $posDeltaCount = 0;
			my $posAvg=0;
			my $curPos = $j+$offset;
			if (exists($deltas{$j})) {
				foreach(@{$deltas{$j}}) {
					push(@deltas, $_);
					$numDeltas++;
					$posAvg += $_;
					$posDeltaCount++;
				}
			}
			$posAvg /= $posDeltaCount if ($posDeltaCount > 0);

			my $posRAvg=0;
			my $numR = 0;
			foreach(@varTypes) {
				$var = $_;
				my $td = 0;
				my $tp = 0;
				$td = $varTypes{$var}->{'td'}->{$j} if (exists($varTypes{$var}->{'td'}->{$j}));
				$tp = $varTypes{$var}->{'tp'}->{$j} if (exists($varTypes{$var}->{'tp'}->{$j}));
				next if ($td < 1 || $tp < 1);
				my $rate = $td/$tp;
				my $numVar = 0;
				if (exists($varTypes{$var}->{'p'}->{$j})) {
					$numVar = $varTypes{$var}->{'p'}->{$j};
				}

				my $weight = $rate*$numVar;

				my $avg = 0;
				my $NNN = 0;
				if (exists($varTypes{$var}->{'d'}->{$j})) {
					foreach(@{$varTypes{$var}->{'d'}->{$j}}) {
						my $diff = $_;
						$avg += $diff;
						$NNN++;
					}
					if ($NNN > 0) {
						$avg /= $NNN;
						$weight /= $NNN;
					}
					foreach(@{$varTypes{$var}->{'d'}->{$j}}) {
						my $diff = $_;
						$rdeltaWeight{$diff} += $weight;
						$posRAvg += $diff*$weight;
						$numR += $weight;;
					}
				}
			}
			$posRAvg /= $numR if ($numR > 0);
			if (!exists($reportedPos{$curPos})) {
				print MDELTAS "$curPos\t$posDeltaCount\t$posAvg\t$posRAvg\n";
				$reportedPos{$curPos}=1;
			}
		}
		
		my $totalWeight = 0;
		foreach(values %rdeltaWeight) {
			$totalWeight+=$_;
		}
		my $scaleFactor = 1;
		my $targetNum = $targetNumOfRandDiffs;
		$targetNum = $numDeltas*10 if ($numDeltas*10 > $targetNum);
		if ($totalWeight > 0) {
			$scaleFactor = $targetNum/$totalWeight;
		}
		#print STDERR "\tscaleFactor=$scaleFactor (numDeltas=$numDeltas)\n";
		foreach(keys %rdeltaWeight) {
			my $diff = $_;
			my $v = $rdeltaWeight{$diff};
			for (my $i=0;$i<$v*$scaleFactor;$i++) {
				push(@rdeltas, $diff);
			}
		}

		my $nn = scalar(@deltas);
		my $rnn = scalar(@rdeltas);

		my $rpvalue = 1;
		my $rpvalue0 = 1;
	
		#i.e. make sure there is some data to analyze
		if ($nn > 3 && $rflag) {
			my $rscriptFile = "$motifOutputDir/$cleanName.rscript.R";
			open RSCRIPT, ">$rscriptFile"  or die "Could not open $rscriptFile for writing\n";
			my $rdataFile = "$motifOutputDir/$cleanName.rdata.txt";
			open RDATA, ">$rdataFile" or die "Could not open $rdataFile for writing\n";

			print RDATA "obs\trand\n";
			my $xnn = $nn;
			$xnn = $rnn if ($rnn > $xnn);
			for (my $k=0;$k<$xnn;$k++) {
				my $v1 = "";
				my $v2 = "";
				$v1 =$deltas[$k] if ($k < $nn);
				$v2 =$rdeltas[$k] if ($k < $rnn);
				print RDATA "$v1\t$v2\n";
			}
			close RDATA;
		
			print RSCRIPT "x <- read.delim(\"$rdataFile\")\n";
			print RSCRIPT "a <- x[,1]\n";
			print RSCRIPT "aa <- a[!is.na(a)]\n";
			print RSCRIPT "b <- x[,2]\n";
			print RSCRIPT "bb <- b[!is.na(b)]\n";
			print RSCRIPT "res <- wilcox.test(a,b)\n";
			print RSCRIPT "res2 <- wilcox.test(a)\n";
			print RSCRIPT "results <- data.frame(tf=c(\"$cleanName\"),pvalueVsRandom=c(res\$p.value),pvalueVsZero=c(res2\$p.value))\n";
			#	if ($firstFlag) {
			#	$firstFlag = 0;
			#} else {
			#	print RSCRIPT "results <- rbind(results, c(\"$mnameStr\",res\$p.value))\n";
			#}
			my $routput = "$motifOutputDir/$cleanName.routput.tsv";
			my $rErr = "$motifOutputDir/$cleanName.rErr.txt";
			print RSCRIPT "write.table(results, file=\"$routput\",sep=\"\t\",quote=FALSE)\n";
			close RSCRIPT;

			`R --no-save < "$rscriptFile" 2> "$rErr"`;
	
			open IN, "$routput" or die "Could not open $routput (routput) for reading\n";
			my $c = 0;
			while (<IN>) {
				$c++;
				next if ($c < 2);
				chomp;
				s/\r//g;
				my @line = split /\t/;
				my $tf = $line[1];
				my $p = $line[2];
				my $p0 = $line[3];
				$rpvalue = $p;
				$rpvalue0 = $p0;
				last;
			}
			close IN;
		}
	
		my $median=Statistics::median(\@deltas);
		my ($avg,$var)=Statistics::avevar(\@deltas);
		my $medianR=Statistics::median(\@rdeltas);
		my ($avgR,$varR) =Statistics::avevar(\@rdeltas);
		my $std = sqrt($var);
		my $stdR = sqrt($varR);
		my $zscore = 'na';
		my $zscore0 = 'na';
		if ($nn > 0) {
			my $sem = $std/sqrt($nn);
			$zscore = ($avg-$avgR)/$sem if ($sem > 0);
			$zscore0 = $avg/$sem if ($sem > 0);
		}
		my $deltaDiff = $avg - $avgR;
	
		my $sitesAboveZero = 0;
		my $sitesAboveAvg = 0;
		foreach(@deltas) {
			$sitesAboveZero++ if ($_ > 0);
			$sitesAboveAvg++ if ($_ > $avgR);
		}
	
		my $rlogpvalue = log($rpvalue+$epsPvalue)/log(10.0);
		$rlogpvalue *= -1 if ($deltaDiff > 0);
		my $rlogpvalue0 = log($rpvalue0+$epsPvalue)/log(10.0);
		$rlogpvalue0 *= -1 if ($avg > 0);
		#print STDERR "\t$cleanName\t$size\t$nn\t$rnn\t$rlogpvalue\t$rlogpvalue0\n";
			
		#print STDERR "\tavg=$avg, avgR=$avgR logp=$rlogpvalue logp0=$rlogpvalue0\n";			
		print STDERR "\t$cleanName\t$sizeStart,$sizeEnd\tavg=$avg\tavgR=$avgR\t(detla=$nn, rdelta=$rnn)\n";

		print STATS "$cleanName\t$size\t$tssDist\t$totalSeqAnalyzed\t$nn"
						. "\t$avg\t$std\t$avgR\t$stdR"
						. "\t$zscore\t$rpvalue\t$rlogpvalue"
						. "\t$zscore0\t$rpvalue0\t$rlogpvalue0\n";

	}
	close MDELTAS;
	close STATS; 

	#now print out a variant summary
	my $varSummaryFile = "$motifOutputDir/variant.summary.tsv";
	open VAR, ">$varSummaryFile" or die "Could not open $varSummaryFile for writing!\n";
	print VAR "Position/Variant";
	foreach(@varTypes) {
		my $var = $_;
		print VAR "\t$var";
	}
	print VAR "\n";
	for (my $i=$minPos;$i<=$maxPos;$i++) {
		my $p = $i + $offset;
		print VAR "$p";
		foreach(@varTypes) {
			my $var = $_;
			my $v = 0;
			if (exists($varTypes{$var}->{'p'}->{$i})) {
				$v = $varTypes{$var}->{'p'}->{$i};
			}
			print VAR "\t$v";
		}
		print VAR "\n";
	}
	close VAR;

	# print out a summary of motif positions
	my $motifDistSummaryFile = "$motifOutputDir/motif.distribution.summary.tsv";
	open MOTIF, ">$motifDistSummaryFile" or die "Could not open $motifDistSummaryFile for writing!\n";
	print MOTIF "Position\tall+\tall-\tstrainSpec+\tstrainSpec-\tgained+\tgained-\tlost+\tlost-\tpotential+\tpotential-\n";
	for (my $i=0;$i<$maxLen;$i++) {
		my $p = $i+$offset;
		print MOTIF "$p";
		print MOTIF "\t" . $motifHist{$i}->{'all'}->{'+'};
		print MOTIF "\t" . $motifHist{$i}->{'all'}->{'-'};
		print MOTIF "\t" . $motifHist{$i}->{'strainSpec'}->{'+'};
		print MOTIF "\t" . $motifHist{$i}->{'strainSpec'}->{'-'};
		print MOTIF "\t" . $motifHist{$i}->{'gained'}->{'+'};
		print MOTIF "\t" . $motifHist{$i}->{'gained'}->{'-'};
		print MOTIF "\t" . $motifHist{$i}->{'lost'}->{'+'};
		print MOTIF "\t" . $motifHist{$i}->{'lost'}->{'-'};
		print MOTIF "\t" . $motifHist{$i}->{'potential'}->{'+'};
		print MOTIF "\t" . $motifHist{$i}->{'potential'}->{'-'};
		print MOTIF "\n";
	}
	close MOTIF;

	if (0) {	
		my $varMatrixFile = "$motifOutputDir/$cleanName.varMatrix.txt";
		my @types = sort {$a cmp $b} keys %motifVarTypes;
		open MATRIX, ">$varMatrixFile";
		print MATRIX "GainedActive";
		foreach(@types) {
			print MATRIX "\t$_";
		}
		print MATRIX "\n";
		for (my $i=0;$i<$mlen;$i++) {
			print MATRIX "$i";
			foreach(@types) {
				my $v = 0;
				if (exists($motifVarPos[$i]->{$_})) {
					$v = $motifVarPos[$i]->{$_};
				}
				print MATRIX "\t$v";
			}
			print MATRIX "\n";
		}
		print MATRIX "LostInactive";
		foreach(@types) {
			print MATRIX "\t$_";
		}
		print MATRIX "\n";
		for (my $i=0;$i<$mlen;$i++) {
			print MATRIX "$i";
			foreach(@types) {
				my $v = 0;
				if (exists($motifVarNeg[$i]->{$_})) {
					$v = $motifVarNeg[$i]->{$_};
				}
				print MATRIX "\t$v";
			}
			print MATRIX "\n";
		}
		close MATRIX;
	}

}
	
sub readVCFFile {
	my ($vcfFile) = @_;
	print LOG "\tReading VCF file: $vcfFile\n";
	print STDERR "\tReading VCF file: $vcfFile\n";
	open IN, $vcfFile;
	my $totalVariantsInFile=0;
	my $totalVariantsConsidered=0;
	my $totalIndelTooLarge = 0;

	while (<IN>) {
		chomp;
		next if (/^#/);
		s/\r//g;


		$totalVariantsInFile++;
		my @line= split /\t/;

		my $chr = $line[0];
		my $pos = $line[1];
		my $ref= $line[3];
		my $alt = $line[4];

		#filtering
		next if ($line[6] ne 'PASS');

		my $vid = $chr . "_" . $pos . "_" . "$ref" . "_" . $alt;

		next if ($onlyIncludeSNV && (length($ref) != length($alt)));
		next if ($onlyIncludeIndels && (length($ref) == length($alt)));

		if (abs(length($ref)-length($alt)) > $maxIndelLength) {
			$totalIndelTooLarge++;
			next;
		}

		if (!exists($vcfData{$chr})) {
			my %a = ();
			$vcfData{$chr} = \%a;
		}
		$vcfData{$chr}->{$vid} = {c=>$chr,p=>$pos,ref=>$ref,alt=>$alt};
		$totalVariantsConsidered++;
	}
	print STDERR "\t\tTotal variants in file: $totalVariantsInFile\n";
	print STDERR "\t\tTotal variants considered: $totalVariantsConsidered\n";
	print LOG "\t\tTotal variants in file: $totalVariantsInFile\n";
	print LOG "\t\tTotal variants considered: $totalVariantsConsidered\n";
	if ($totalIndelTooLarge>0) {
		print STDERR "\t\tIndels that were too large (>$maxIndelLength): $totalIndelTooLarge\n";
		print LOG "\t\tIndels that were too large (>$maxIndelLength): $totalIndelTooLarge\n";
	}
	return \%vcfData;
}

sub assignVariantsToPeaks {
	my ($qtlData,$vcfData,$peaksFile,$flip,$minStart,$maxEnd) = @_;

	`bed2pos.pl "$peaksFile" -check > "$tmpFile"`;
	my %peaks = ();
	open IN, $tmpFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		my $og = $_;
		next if (/^#/);
		my @line = split /\t/;
		next if (@line < 5);
		my $m = floor(($line[2]+$line[3])/2);
		my $chr = $line[1];
		my $id = $line[0];
		my $d = $line[4];
		my $v = 'NA';
		$v = $line[5] if (@line > 5);

		if (!exists($peaks{$chr})) {
			my %a = ();
			$peaks{$chr} = \%a;
		}
		$peaks{$chr}->{$id} = {id=>$id,c=>$chr,p=>$m,s=>$line[2],e=>$line[3],d=>$d,v=>$v};
	}
	close IN;
	`rm "$tmpFile"`;

	my $totalAssociations = 0;
	my @chrs = sort {$a cmp $b} keys %peaks;
	for (my $z=0;$z<@chrs;$z++) {
		my $chr = $chrs[$z];
		print STDERR "\t\t$chr\n";

		my @vcfIDs = sort {$vcfData->{$chr}->{$a}->{'p'} <=> $vcfData->{$chr}->{$b}->{'p'}} keys %{$vcfData->{$chr}};
		my $vcfIndex = 0;
		my @peakIDs = sort {$peaks{$chr}->{$a}->{'p'} <=> $peaks{$chr}->{$b}->{'p'}} keys %{$peaks{$chr}};

		my $maxRange = abs($minStart);
		$maxRange = abs($maxEnd) if ($maxRange < abs($maxEnd));

		for (my $i=0;$i<@peakIDs;$i++) {
			my $tid = $peakIDs[$i];
			my $tp = $peaks{$chr}->{$tid}->{'p'};
			my $td = $peaks{$chr}->{$tid}->{'d'};
			my $v = $peaks{$chr}->{$tid}->{'v'};
			my $lower = $tp+$minStart;
			my $upper = $tp+$maxEnd;
			if ($td eq '-' || $td eq '1') {
				$lower = $tp-$maxEnd;
				$upper = $tp-$minStart;
			}

			for (my $j=$vcfIndex;$j<@vcfIDs;$j++) {
				my $vid = $vcfIDs[$j];
				my $vp = $vcfData->{$chr}->{$vid}->{'p'};
				if ($vp < $tp-$maxRange) {	
					$vcfIndex++;
					next;
				}
				last if ($vp > $tp+$maxRange);
				if (($vp >= $lower) && ($vp <= $upper)) {

					my $ref = $vcfData->{$chr}->{$vid}->{'ref'};
					my $alt = $vcfData->{$chr}->{$vid}->{'alt'};
					my $maf = 'NA';
					my $pvalue = 'NA';
					my $slope = $v;
					my $tssDist = $vp-$tp;
					if ($td eq '-' || $td eq '1') {
						$tssDist *= -1;
					}
					# score of -1 implies the reference allele is the 'active' phenotype, score of 1 means the alt is 'active'
					my $score = -1;
					if ($flip) {
						$score = 1;
					}


					my $qtlID = $vid . "-" . $tid;

					my $infoStr = "$qtlID\t$vid\t$tid\t$chr\t$tp\t$tp\t$td\t$score\t$pvalue\t$slope";

					$qtlData->{$qtlID} = {
						vid=>$vid,
						c=>$chr,
						vp=>$vp,
						tssID=>$tid,
						tp=>$tp,
						td=>$td,
						ref=>$ref,
						alt=>$alt,
						maf=>$maf,
						pvalue=>$pvalue,
						slope=>$slope,
						score=>$score,
						tssDist=>$tssDist,
						info=>$infoStr
					};
					$totalAssociations++;
				}
			}
		}
	}
	print STDERR "\t$totalAssociations total variant-feature associations found\n";
	print LOG "\t$totalAssociations total variant-feature associations found\n";
}

sub exportQTLData {
	my ($qtlData) = @_;

	open QTL, ">$outputDir/qtls.tsv";
	open PEAKS, ">$outputDir/qtlFeatures.tsv";
	open BED, ">$outputDir/variants.bed";
	print BED "track name=\"Genetic Variation: $outputDir\" visibility=3 useScore=0\n";

	my %printedVariants = ();

	my %peaks = ();
	print QTL "#variant_id\tphenotype_id\ttss_distance\tna_count\tna_samples\tmaf\tpvalue\tslope\n";
	foreach (values %$qtlData) {
		my $vid = $_->{'vid'};
		my $tid = $_->{'tssID'};
		my $tssDist = $_->{'tssDist'};
		my $count= 'NA';
		my $samples= 'NA';
		my $maf = $_->{'maf'};
		my $pvalue = $_->{'pvalue'};
		my $slope = $_->{'score'};

		# score/slope is used to assign which allele is 'active'
		# -1 (or negative slope) means the reference is the 'active' phenotype
		# 1 (or positive slope) means the alt is the 'active' phenotype
		print QTL "$vid\t$tid\t$tssDist\t$count\t$samples\t$maf\t$pvalue\t$slope\n";

		my $c = $_->{'c'};
		my $tp = $_->{'tp'};
		my $td = $_->{'td'};
		my $x = $c . "-" . $tp . "-" . $td;
		$peaks{$x} = "$tid\t$c\t$tp\t$tp\t$td";

		if (!exists($printedVariants{$vid})) {
			if ($vid =~ /^(.+?)\_(\d+?)\_(.+?)\_(.+?)$/) {
				my $c = $1;
				my $p = $2;
				my $r = $3;
				my $a = $4;
				my $s = $p-1;
				my $e = $s+length($r);
				print BED "$c\t$s\t$e\t$r" . "->" . "$a\n";
			} else {
				print STDERR " !! warning: couldn't parse vid: $vid \n";
			}
			$printedVariants{$vid}=1;
		}
	}
	
	print PEAKS "#phenotype_id\tchr\tstart\tend\tstrand\n";
	foreach(values %peaks) {
		print PEAKS "$_\n";
	}

	close QTL;
	close PEAKS;
	close BED;

}

sub readQTLFile {
	my ($qtlFile,$peakPosFile,$minStart,$maxEnd) = @_;

	my %qtl =();

	`bed2pos.pl "$peakPosFile" -check > "$tmpFile"`;
	my %info = ();
	open IN, $tmpFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		my $og = $_;
		next if (/^#/);
		my @line = split /\t/;
		next if (@line < 5);
		my $m = floor(($line[2]+$line[3])/2);
		$info{$line[0]} = {id=>$line[0],c=>$line[1],p=>$m,s=>$line[2],e=>$line[3],d=>$line[4],og=>$og};
	}
	close IN;
	`rm "$tmpFile"`;

	my %varPos = ();

	my $total = 0;
	my $filteredBadMAF=0;
	my $filteredBadPvalue = 0;
	my $filteredNonBinary = 0;
	my $filteredWildcard= 0;
	my $filteredBadSlope = 0;
	my $totalIndelTooLarge = 0;
	my $removedNotSNV = 0;
	my $removedNotIndel = 0;
	my $notInRange = 0;
	my $good = 0;
	my $totalPositiveScore =0;
	my $noTSSFound = 0;
	my $xx = 0;
	my $tempCount = 0;

my %counts = ();


	open IN, $qtlFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		$total++;
		my @line = split /\t/;

		my $vid = $line[0];
		my $tssID = $line[1];
		my $tssDist = $line[2];
		my $count = $line[3];
		my $sampleCount = $line[4];
		my $maf = $line[5];
		my $pvalue = $line[6];
		my $slope = $line[7];

		my $chr = $line[1];
		my $vp = $line[9];
		my $ref = $line[10];
		my $alt = $line[11];
		my $rate = $maf;

		if ($vid =~ /(chr.*?)\_(\d+)\_(.+?)\_(.+)$/) {
			$chr = $1;
			$vp = $2;
			$ref = $3;
			$alt = $4;
			if ($alt =~ /\_/ && $onlyIncludeBinaryVariants) {
				print STDERR "!!! Warning: multiple alt alleles: $vid\n";
				$filteredNonBinary++;
				next;
			}
			if ($ref eq '*' || $alt eq '*') {
				$filteredWildcard++;
				next;
			}
		} else {
			print STDERR "!!! Warning: Could not parse variant ID: $vid\n";
			next;
		}

		if (!exists($info{$tssID})) {
			print STDERR "!!! Warning: featureID ($tssID) not found in file: $peakPosFile\n" if ($noTSSFound < 30);
			print STDERR "Additional featureIDs not found...\n" if ($noTSSFound == 30);
			$noTSSFound++;
			next;
		}
		my $tc = $info{$tssID}->{'c'};
		if ($tc ne $chr) {
			print STDERR "!!! Warning feature chr different than variant: $vid ($info{$tssID}->{'c'}) vs. $tssID ($tc)\n";
			next;
		}
		my $tp = $info{$tssID}->{'p'};
		my $td = $info{$tssID}->{'d'};

		#calculate variant distance relative to feature/TSS
		$tssDist = $vp-$tp;

		$tssRef = $ref;
		$tssAlt = $alt;
		if ($td eq '-' || $td eq '1') {
			$tssDist *= -1;
			$tssRef = HomerConfig::revopp($ref);
			$tssAlt = HomerConfig::revopp($alt);
		}

		if ($tssDist < $minStart || $tssDist >= $maxEnd) {
			$notInRange++;
			next;
		}

		#filtering of variants:
		if ($onlyIncludeSNV && ((length($ref) != length($alt)) || $ref eq '-' || $alt eq '-')) {
			$removedNotSNV++;
			next;
		}
		if ($onlyIncludeIndels && (length($ref) == length($alt))) {
			$removedNotIndel++;
			next;
		}
		if (abs(length($ref)-length($alt)) > $maxIndelLength) {
			$totalIndelTooLarge++;
			next;
		}

		my $qtlID = $vid . "-" . $tssID;
		if (exists($qtl{$vid})) {
			print STDERR "!!! Warning - duplicate feature-variant entry ($vid, $tssID)\n";
			next;
		}

		if ($onlyIncludeBinaryVariants) {
			my $vp = $chr . ":" . $vp;
			if (exists($varPos{$vp})) {
				if ($varPos{$vp} eq 'bad') {
					$filteredNonBinary++;
					next;
				}
				if ($varPos{$vp} ne $vid) {
					$varPos{$vp} = 'bad';
					delete $qtl{$qtlID} if (exists($qtl{$qtlID}));
					$filteredNonBinary++;
					next;
				}			
			} else {
				$varPos{$vp} = $vid;
			}
		}
	
		my $score = $slope;

		if (($maf ne '' && $maf ne 'na' && $maf ne 'NA') && ($maf < $mafThreshold || $maf > 1-$mafThreshold)) {
			$filteredBadMAF++;
			next;
		}
		#if ($rate < $varMin || $rate > 67-$varMin) {
		#	$tooFewVariants++;
		#	next;
		#}
		if (abs($score) < $absSlopeMin) {
			$filteredBadSlope++;
			next;
		}
		my $logP = 0;
		if ($pvalue eq '' || $pvalue eq 'NA') {
			if ($pvalueThreshold < 1.0) {
				$filteredBadPvalue++;
				next;
			}
		} else {
			if ($pvalue > $pvalueThreshold) {
				$filteredBadPvalue++;
				next;
			}
			if ($pvalue < 1e-200) {
				$logP = log(1e-200)/log(10);
			} else {
				$logP = log($pvalue)/log(10);
			}
			if ($useSignedPvalues) {
				if ($score > 0) {
					$score = -1*$logP;
				} else {
					$score = $logP;
				}
			}
		}

		if ($tssStrandTestMode eq "+") {
			next if ($td eq '-');
		} elsif ($tssStrandTestMode eq "-") {
			next if ($td eq '+');
		}

$counts{$tssRef . ">" . $tssAlt}++;

		my $infoStr = "$qtlID\t$vid\t$tssID\t$chr\t$tp\t$tp\t$td\t$score\t$pvalue\t$slope";
		$qtl{$qtlID} = {
			vid=>$vid,
			c=>$chr,
			vp=>$vp,
			tssID=>$tssID,
			tp=>$tp,
			td=>$td,
			ref=>$ref,
			alt=>$alt,
			tssRef=>$tssRef,
			tssAlt=>$tssAlt,
			maf=>$maf,
			pvalue=>$pvalue,
			slope=>$slope,
			score=>$score,
			tssDist=>$tssDist,
			info=>$infoStr
		};
		$good++;
		$totalPositiveScore++ if ($score > 0);
	}
	close IN;


if (0) {
my @x = sort {$a cmp $b} keys %counts;
foreach(@x) {
	print STDERR "-- $_\t$counts{$_}\n";
}
}

	print LOG "\tFinished Processing QTL input files:\n";
	print LOG "\t\t$total total associations\n";
	print LOG "\t\t$noTSSFound filtered due to no matching feature/TSS ID (could be on purpose)\n";
	print LOG "\t\t$filteredBadMAF filtered for low MAF ($mafThreshold, 1-$mafThreshold)\n";
	print LOG "\t\t$filteredNonBinary filtered for not being binary variants\n";
	print LOG "\t\t$filteredWildcard filtered for having wildcard allele\n";
	print LOG "\t\t$removedNotSNV removed because they weren't a simple SNV\n";
	print LOG "\t\t$removedNotIndel removed because they weren't an indel\n";
	print LOG "\t\t$totalIndelTooLarge removed because indel was too large (> $maxIndelLength)\n";
	print LOG "\t\t$notInRange removed being outside of local distance range ($minStart, $maxEnd)\n";
	print LOG "\t\t$filteredBadPvalue removed for bad p-value: ($pvalueThreshold)\n";
	print LOG "\t\t$filteredBadSlope removed for bad slope/score (|slope| < $absSlopeMin)\n";
	print LOG "\tTotal tss-variant associations considered: $good\n";
	my $r = 0;
	if ($good > 0) {
		$r = $totalPositiveScore/$good;
	}
	print LOG "\tTotal tss-variant associations with + scores: $totalPositiveScore ($r)\n";

	print STDERR "\tTotal tss-variant associations considered: $good out of $total\n";
	print STDERR "\tTotal tss-variant associations with + scores: $totalPositiveScore ($r)\n";

	return \%qtl;
}

# getSeqQTL
#   $minStart is the offset for the beginnin of the sequence/end
sub getSeqQTL {
	my ($data,$qtlData,$minStart,$maxEnd,$subsetFile) = @_;

	my %valid = ();
	if ($subsetFile ne '') {
		open IN, $subsetFile;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			$valid{$line[0]} = 1;
		}
		close IN;
	}

	open OUT, ">$tmpFile";
	my %tssInfo = ();
	foreach(keys %{$qtlData}) {
		my $id = $_;
		my $vid = $qtlData->{$id}->{'vid'};
		next if ($subsetFile ne '' && !exists($valid{$vid}));
		my $c = $qtlData->{$id}->{'c'};
		my $tssID = $qtlData->{$id}->{'tssID'};
		my $tp = $qtlData->{$id}->{'tp'};
		my $td = $qtlData->{$id}->{'td'};
		my $score = $qtlData->{$id}->{'score'};
		$tssInfo{$id} = {tp=>$tp,td=>$td,s=>$score,tssID=>$tssID};
		my $s = $tp + $minStart;
		my $e = $tp + $maxEnd;
		if ($td eq '-' || $td eq '1') {
			$s = $tp - $maxEnd;
			$e = $tp - $minStart;
		}
		print OUT "$id\t$c\t$s\t$e\t$td\n";
	}
	close OUT;
	`homerTools extract $tmpFile $genome1 $mask > $tmpFile2`;

	my $numQTLs = scalar(keys %$qtlData);
	if ($randomizeVariantAssignments) {
		print STDERR "\tRandomizing variant-feature assignments (preserving position/nucleotide)\n";
		open IN, $tmpFile2;
		my @seqs = ();
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			my $x = {id=>$line[0],seq=>$line[1]};
			push(@seqs, $x);
		}
		my $numSeqs = scalar(@seqs);
		close IN;
print STDERR "numSeqs = $numSeqs\n";
		my %taken = ();
		my $maxTries = 10000;
		my $totalAssigned = 0;

		my @ogseqs = @seqs;

		open OUT, ">$tmpFile2";
		my $numAssigned = 0;
		my @qtlIDs = keys %$qtlData;
		my $numQtlIDs = scalar(@qtlIDs);
print STDERR "numSeqs = $numSeqs\tnumQtlIDs=$numQtlIDs\n";

		foreach(@qtlIDs) {
			my $qid = $_;
			my $tssDist = $qtlData->{$qid}->{'tssDist'};
			my $tssDir = $qtlData->{$qid}->{'td'};
			my $score = $qtlData->{$qid}->{'score'};
			my $ref = $qtlData->{$qid}->{'tssRef'};
			my $alt = $qtlData->{$qid}->{'tssAlt'};

			my $assigned = 0;
			my $tries = 1;

			while (!$assigned) {
				$tries++;
				if ($tries > $maxTries) {
print STDERR "!!!! Might be some duplicate sequences (at $totalAssigned of $numQtlIDs, $numSeqs left to be assigned)!!!!!\n";
					@seqs = @ogseqs;
					$numSeqs = scalar(@seqs);
					%taken = ();
					$tries = 0;
					next;
				}
				my $index = floor(rand()*$numSeqs);
				my $ptssID = $seqs[$index]->{'id'};

				#my $t = $ptssID . "_" . $tssDist;
				my $t = $ptssID;
				next if (exists($taken{$t}));

				# only assign if both scores are positive or both scores are negative
				if ($preserveActiveStatusForRand) {
					next if ($score * $qtlData->{$ptssID}->{'score'} < 0);
				}


				my $s = $seqs[$index]->{'seq'};

				#my $subx = $tssDist-$minStart;
				my $nuc = substr($s, $tssDist - $minStart,1);
				next if ($nuc ne $ref);

				splice(@seqs,$index,1);
				$numSeqs--;
				$taken{$t}=1;
				$assigned=1;
				$totalAssigned++;


#print STDERR "$qid\t$tssDist\t$score\t$ref\t$alt\t(nuc=$nuc, subx=$subx) $qtlData->{$qid}->{'tssID'} $qtlData->{$qid}->{'tp'} $qtlData->{$qid}->{'td'}\n";
#print STDERR "\t$og\n\t$s\n";

#print STDERR "$qtlData->{$qid}->{'tssID'} = $ptssID\n";
				$qtlData->{$qid}->{'tssID'} = $tssInfo{$ptssID}->{'tssID'};
				$qtlData->{$qid}->{'tp'} = $tssInfo{$ptssID}->{'tp'};
				$qtlData->{$qid}->{'td'} = $tssInfo{$ptssID}->{'td'};
				print OUT "$qid\t$s\n";


			}
			if ($assigned == 0) {
				print STDERR "!!! $totalAssigned of $numQTLs assigned (numSeqs=$numSeqs)\n";
				print STDERR "$qid\t$tssDist\t$score\t$ref\t$alt\n";
				print STDERR "!!! Error - could not assign a random feature for QTL $qid\n";
				exit;
			} else {
					#print STDERR "\t$qid assigned!\n";
			}
		}
		close OUT;
		print STDERR "\t\tFinished random assignment ($totalAssigned)\n";
	}

	my $indelTooLarge = 0;
	my $indelAtEdge = 0;
	my $tssDistNotRightRange = 0;

	my $N = '';
	for (my $i=0;$i<$maxIndelLength;$i++) {
		$N .= "N";
	}

	open OUT, ">$outputDir/seqs.tsv";
	my %varSummary = ();
	my %varKeys = ();

	open IN, $tmpFile2;
	while (<IN> ){
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $line[0];
		my $seq = $line[1];

		if (!exists($qtlData->{$id})) {
			print STDERR "!!! somethings wrong - didn't get correct qtlID...\n";
			exit;
		}
		#my $ref = $qtlData->{$id}->{'ref'};
		#my $alt = $qtlData->{$id}->{'alt'};
		my $ref = $qtlData->{$id}->{'tssRef'};
		my $alt = $qtlData->{$id}->{'tssAlt'};
		my $tssDist = $qtlData->{$id}->{'tssDist'};
		my $tssStrand = $qtlData->{$id}->{'td'};
		my $tp = $qtlData->{$id}->{'tp'};
		my $td = $qtlData->{$id}->{'td'};
		my $chr = $qtlData->{$id}->{'c'};

		my $overallLength = length($seq);

		if ($tssDist < $minStart || $tssDist >= $maxEnd) {
			#print STDERR "\t$tssDist ($minStart, $maxEnd)\n";
			$tssDistNotRightRange++;
			next;
		}
		my $lenRef = length($ref);
		my $lenAlt = length($alt);
		my $lenDiff = $lenAlt-$lenRef;
		my $maxSize = $lenRef;
		$maxSize = $lenAlt if ($lenAlt>$maxSize);
		my $sizeChange = abs($lenRef-$lenAlt);

		if ($sizeChange >= $maxIndelLength) {
			$indelTooLarge++;
			next;
		}

		my $aseq = '';
		my $pos = '';

		#for randomization leter:
		my $refAbsPos = -1;
		my $altAbsPos = -1;
		my $refForRand = $ref;
		my $altForRand = $alt;
		#print STDERR "$id\t$tssDist\t$tssStrand\t$ref\t$alt\n$seq\n";


		if ($tssStrand eq '-' || $tssStrand eq '1') {
			#my $revOppRef = HomerConfig::revopp($ref);
			#my $revOppAlt = HomerConfig::revopp($alt);
			my $revOppRef = $ref;
			my $revOppAlt = $alt;

			$refForRand = $revOppRef;
			$altForRand = $revOppAlt;

			$pos = $tssDist-$minStart+1;
			if ($pos-$sizeChange < 0 || $pos+$sizeChange >= $overallLength) {
				$indelAtEdge++;
				next;
			}
			$refAbsPos = $pos-length($revOppRef);

			if ($tssDist < 0) {
				if ($lenDiff < 0) {
					$aseq .= substr($N,0,$lenDiff*-1);
					$aseq .= substr($seq,0,$pos-$lenRef);
					$altAbsPos = length($aseq);
					$aseq .= $revOppAlt;
					$aseq .= substr($seq,$pos,$overallLength-$pos);
				} elsif ($lenDiff > 0) {
					$aseq .= substr($seq,$lenDiff,$pos-$lenDiff-$lenRef);
					$altAbsPos = length($aseq);
					$aseq .= $revOppAlt;
					$aseq .= substr($seq,$pos,$overallLength-$pos);
				} else {
					$aseq .= substr($seq,0,$pos-length($revOppRef));
					$altAbsPos = length($aseq);
					$aseq .= $revOppAlt;
					$aseq .= substr($seq,$pos,$overallLength-$pos);
				}
			} else {
				my $endLen = $overallLength-$pos;
				my $addedNs = '';
				if ($lenDiff > 0) {
					$endLen -= $lenDiff;
				} else {
					$addedNs = substr($N,0,abs($lenDiff));
				}
				$aseq .= substr($seq,0,$pos-$lenRef) . $revOppAlt . substr($seq,$pos,$endLen) . $addedNs;
				$altAbsPos = $pos-$lenRef;
			}
			#$aseq = substr($seq,0,$pos-length($revOppRef)) . $revOppAlt . substr($seq,$pos);
		} else {
			
			$pos = $tssDist-$minStart;
			if ($pos-$sizeChange < 0 || $pos+$sizeChange >= $overallLength) {
				$indelAtEdge++;
				next;
			}
			$refAbsPos = $pos;

			if ($tssDist < 0) {
				if ($lenDiff < 0) {
					$aseq .= substr($N,0,$lenDiff*-1);
					$aseq .= substr($seq,0,$pos);
					$altAbsPos = length($aseq);
					$aseq .= $alt;
					my $endLen = $overallLength-($pos+$lenRef);
					$aseq .= substr($seq,$pos+$lenRef,$endLen);
				} elsif ($lenDiff > 0) {
					$aseq .= substr($seq,$lenDiff,$pos-$lenDiff);
					$altAbsPos = length($aseq);
					$aseq .= $alt;
					my $endLen = $overallLength-($pos+$lenRef);
					$aseq .= substr($seq,$pos+$lenRef,$endLen);
				} else {
					$aseq .= substr($seq,0,$pos);
					$altAbsPos = length($aseq);
					$aseq .= $alt;
					my $endLen = $overallLength-($pos+$lenRef);
					$aseq .= substr($seq,$pos+$lenRef,$endLen);
				}
			} else {
				my $endLen = $overallLength-($pos+$lenRef);
				my $addedNs = '';
				if ($lenDiff > 0) {
					$endLen -= $lenDiff;
				} else {
					$addedNs = substr($N,0,abs($lenDiff));
				}
				$aseq .= substr($seq,0,$pos) . $alt . substr($seq,$pos+$lenRef,$endLen) . $addedNs;
				$altAbsPos = $pos;
			}
						
			#-length($ref);
		}
		#print STDERR "$aseq\n\n";

		#print STDERR "\t$id\t$tssStrand\t$pos ($minStart,$maxEnd)\t$tssDist\t$ref -> $alt\n";
		#print STDERR "\t\t$seq\n\t\t$aseq\n";
		print OUT "$id-REF-$ref\t$seq\n";
		print OUT "$id-ALT-$alt\t$aseq\n";

		my $score = $qtlData->{$id}->{'score'};
		$score *= -1 if ($flipScore);
		my $flip = 'na';

		if (0) {
			if ($scoreTestMode eq '') {
			} elsif ($scoreTestMode eq '+') {
				next if ($score < 0);
			} elsif ($scoreTestMode eq '-') {
				next if ($score > 0);
			}
		}

		my @s1randomizations = ();
		my @s2randomizations = ();

		# score: - scores mean the reference is 'active', + scores mean the alternative is 'active'
		# s1 = control (or inactive) sequence/allele
		# s2 = active sequence/allele
		# phenotype = which allele is active
		if ($score >= 0) {
			$data->{$id} = {name=>$id,s1=>$seq,s2=>$aseq,f=>$flip,info=>$qtlData->{$id}->{'info'},tssDist=>$tssDist,phenotype=>"alt",
				refAbsPos=>$refAbsPos, altAbsPos=>$altAbsPos, refForRand=>$refForRand,altForRand=>$altForRand,refFlip=>0,refActive=>0,
				s1rand=>\@s1randomizations,s2rand=>\@s2randomizations,score=>$score,c=>$chr,tp=>$tp,td=>$td,seq=>$seq};
		} else {
			$data->{$id} = {name=>$id,s1=>$aseq,s2=>$seq,f=>$flip,info=>$qtlData->{$id}->{'info'},tssDist=>$tssDist,phenotype=>"ref",
				refAbsPos=>$refAbsPos, altAbsPos=>$altAbsPos, refForRand=>$refForRand,altForRand=>$altForRand,refFlip=>1,refActive=>1,
				s1rand=>\@s1randomizations,s2rand=>\@s2randomizations,score=>$score,c=>$chr,tp=>$tp,td=>$td,seq=>$seq};
		}

		my $varKey = $altForRand . "->" . $refForRand;
		if ($score >= 0) {
			$varKey = $refForRand . "->" . $altForRand;
		}
		if (!exists($varSummary{$refAbsPos})) {
			my %a = ();
			$varSummary{$refAbsPos} = \%a;
		}
		$varSummary{$refAbsPos}->{$varKey}++;
		$varKeys{$varKey}=1;
	}
	close IN;
	close OUT;

	if (0) {
		print STDERR "\t\tremoved because:\n";
		print STDERR "\t\t\tIndel too large: $indelTooLarge\n";
		print STDERR "\t\t\tIndel at edge:   $indelAtEdge\n";
		print STDERR "\t\t\tVariant not correction range to TSS: $tssDistNotRightRange\n";
	}

	open VAR, ">$outputDir/variantSummary.tsv";
	my @keys = sort {$a cmp $b} keys %varKeys;
	print VAR "Position";
	foreach(@keys) {
		print VAR "\t$_";
	}
	print VAR "\n";
	my @pos = sort {$a <=> $b} keys %varSummary;
	foreach(@pos) {
		my $p = $_;
		my $pp = $p + $minStart;
		print VAR "$pp";
		foreach(@keys) {
			my $k = $_;
			my $v = 0;
			if (exists($varSummary{$p}->{$k})) {
				$v = $varSummary{$p}->{$k};
			}
			print VAR "\t$v";
		}
		print VAR "\n";
	}
	close VAR;


	`rm $tmpFile $tmpFile2`;
}


sub getSeqLiftOver {
	my ($data,$peaks,$flip,$minStart,$maxEnd) = @_;

	my %seq = ();
	if ($liftOver) {
		`adjustPeakFile.pl "$peaks" -size $liftOverSize > "$tmpFile"`;
		`convertCoordinates.pl "$liftOver" "$tmpFile" "$tmpFile3" -peaks -p $maxCPUs3`;
	} else {
		print STDERR "For now the program requires a liftOver file\n";
		exit;
	}
	my $len = $maxEnd-$minStart+1;
	`adjustPeakFile.pl "$peaks" -size $minStart,$maxEnd > "$tmpFile"`;
	`adjustPeakFile.pl "$tmpFile3" -size $minStart,$maxEnd > "$tmpFile2"`;

	my %peakInfo = ();
	open IN, $tmpFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		my $og = $_;
		my @line = split /\t/;
		my $strand = "+";
		if ($line[4] eq '-' || $line[4] eq '1') {
			$strand = "-";
		}
		if ($tssStrandTestMode eq "+") {
			next if ($strand eq '-');
		} elsif ($tssStrandTestMode eq "-") {
			next if ($strand eq '+');
		}
		$peakInfo{$line[0]} = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]";
	}
	close IN;


	`homerTools extract "$tmpFile" "$genome1" $mask > "$tmpFile3"`;
	`homerTools extract "$tmpFile2" "$genome2" $mask > "$tmpFile4"`;
	my %seq1 = ();
	my $totalG1 = 0;
	my $totalG2 = 0;
	my $totalCommon = 0;
	open IN, $tmpFile3;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if (@line < 2);
		my $id = $line[0];
		my $s = $line[1];
		next if ($s eq '');
		next if (!exists($peakInfo{$id}));
		$seq1{$id} = $s;
		$totalG1++;
	}
	close IN;

	my $differentNumNs=0;
	open OUT, ">$outputDir/seqs.tsv";
	open IN, $tmpFile4;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if (@line < 2);
		my $id = $line[0];
		my $s = $line[1];
		next if ($s eq '');
		next if (!exists($peakInfo{$id}));
		$totalG2++;
		if (!exists($seq1{$id})) {
			next;
		}
		my $s1len = length($seq1{$id});
	   	my $s2len =	 length($s);
		if ($len != $s1len || $len != $s2len) {
			print STDERR "\t!!! These sequences had different lengths: $id\t$len\t$s1len\t$s2len\n";
			next;
		}
		if ($sameNflag) {
			my $s1n = $seq1{$id};
			$s1n =~ s/N//i;
			my $s2n = $s;
			$s2n =~ s/N//i;
			if (length($s1n) != length($s2n)) {
				$differentNumNs++;
				next;
			}
		}


		$totalCommon++;
		print OUT "$id-g1\t$seq1{$id}\n";
		print OUT "$id-g2\t$s\n";

		my $pheno = "ref";
		my $d = {name=>$id,s1=>$seq1{$id},s2=>$s,f=>$flip,info=>$peakInfo{$id},tssDist=>".",phenotype=>$pheno};
		if ($flip) {
			$pheno = "alt";
			$d = {name=>$id,s1=>$s,s2=>$seq1{$id},f=>$flip,info=>$peakInfo{$id},tssDist=>".",phenotype=>$pheno};
		}
		$data->{$id} = $d;
	}
	close IN;
	print LOG "\tSequences extracted from g1: $totalG1\n";
	print LOG "\tSequences extracted from g2: $totalG2\n";
	print LOG "\tSequences with different Ns: $differentNumNs\n";
	print LOG "\tSequences extracted common : $totalCommon\n";
	close OUT;
	`rm $tmpFile3 $tmpFile4 $tmpFile2 $tmpFile`;
}
sub addToHist {
	my ($hist,$tfcat,$bin,$corr) = @_;
	if (!exists($hist->{$tfcat})) {
		my %a = ();
		my %b = ();
		$hist->{$tfcat} = {v=>\%a,n=>\%b};
	}
	$hist->{$tfcat}->{'v'}->{$bin}+=$corr;
	$hist->{$tfcat}->{'n'}->{$bin}+=1.0;
	$binMin = $bin if ($bin < $binMin);
	$binMax = $bin if ($bin >$binMax);
}

sub organizeMotifsByLength {
	my ($motifs) = @_;
	
	my %mset = ();
	foreach(@$motifs) {
		my $m = $_;
		my $mlen = $m->{'len'};
		if (!exists($mset{$mlen})) {
			my @a = ();
			$mset{$mlen} = \@a;
		}
		push(@{$mset{$mlen}}, $m);

	}
	return \%mset;
}
sub readGWAS {
	my ($file) = @_;
	my %gwas = ();
	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line =split /\t/;
		my $chr = $line[1];
		my $pos = $line[3];
		my $rsID = $line[4];
		my $pmid = $line[5];
		my $studyName = $line[9];
		my $trait = $line[10];
		my $variant = $line[15];
		my $pvalue = $line[17];
		my $str = $rsID . "|" . $pmid . "|" . $studyName . "|" . $trait . "|" . $variant . "|" . $pvalue;

		my $id = $chr . "-" . $pos;
		if (!exists($gwas{$id})) {
			$gwas{$id} = $str;
		} else {
			$gwas{$id} .= "," . $str;
		}
	}
	return \%gwas;
}
sub expandWithRandomizedMotifs {
	my ($motifs,$zscore,$motifRandOrder) = @_;

	my @newMotifs = ();
	foreach(@$motifs) {
		my $m = $_;
		push(@newMotifs, $m);
		HomerConfig::printMotifFile($m,$tmpFile);
		`randomizeMotif.pl "$tmpFile" -n $zscore -order $motifRandOrder > "$tmpFile2"`;
		my $rmotifs = HomerConfig::readMotifFile($tmpFile2);
		foreach(@$rmotifs) {
			push(@newMotifs, $_);
		}
	}
	my $n1 = scalar(@$motifs);
	my $n2 = scalar(@newMotifs);
	print STDERR "\tAdded randomized motifs (n=$zscore), expanded motifs to analyze from $n1 to $n2\n";
	return \@newMotifs;
}
sub printSequences {
	my ($filename, $seqs) = @_;
	open OUT, ">$filename";
	foreach(keys %$seqs) {
		my $id = $_;
		my $s1 = $seqs->{$id}->{'s1'};
		my $s2 = $seqs->{$id}->{'s2'};
		next if ($s1 eq $s2);
		print OUT "$id-g1\t$s1\n";
		print OUT "$id-g2\t$s2\n";
	}
	close OUT;
}
sub UseRtoCorrectPvalues {
	my ($pvalues) = @_;
	open OUT, ">$tmpFile";
	print OUT "pvalues\n";
	foreach(@$pvalues) {
		my $v = $_;
		if ($v eq '' || $v eq 'NA' || $v eq 'na') {
			$v = 1;
		}
		print OUT "$v\n";
	}
	close OUT;
	open OUT, ">$tmpFile2";
	print OUT "x <- read.delim(\"$tmpFile\")\n";
	print OUT "xx <- p.adjust(t(x), method=\"BH\")\n";
	print OUT "write.csv(xx,file=\"$tmpFile3\")\n";
	close OUT;
	`R --no-save < $tmpFile2`;

	open IN, $tmpFile3;
	my $index = 0;
	my @rv = ();
	my $z = 0;
	my $min = 1;
	while (<IN>) {
		$z++;
		next if ($z < 2);
		chomp;
		s/\r//g;
		my @line = split /\,/;
		my $p = $line[1];
		$min = $p if ($p > 0 && $p < $min);
		push(@rv, $p);
	}
	close IN;

	#`rm $tmpFile $tmpFile2 $tmpFile3`;


	for (my $i=0;$i<@rv;$i++) {
		$r = $rv[$i];
		if ($r < $min) {
			$r = $min;
		}
		$rv[$i] = $r;
	}
	return \@rv;
}

sub UseRtoCorrectLogPvalues {
	my ($pvalues) = @_;
	open OUT, ">$tmpFile";
	print OUT "pvalues\n";
	foreach(@$pvalues) {
		my $v = $_;
print STDERR "\t$v\n";
		if ($v > 0) {
			$v = 10**(-1*$v);
		} else {
			$v = 10**$v;
		}
print STDERR "\t$v\n";
		print OUT "$v\n";
	}
	close OUT;
	open OUT, ">$tmpFile2";
	print OUT "x <- read.delim(\"$tmpFile\")\n";
	print OUT "xx <- p.adjust(t(x), method=\"BH\")\n";
	print OUT "write.csv(xx,file=\"$tmpFile3\")\n";
	close OUT;
	`R --no-save < $tmpFile2`;

	open IN, $tmpFile3;
	my $index = 0;
	my @rv = ();
	my $z = 0;
	my $min = 1;
	while (<IN>) {
		$z++;
		next if ($z < 2);
		chomp;
		s/\r//g;
		my @line = split /\,/;
		my $p = $line[1];
		$min = $p if ($p > 0 && $p < $min);
		push(@rv, $p);
	}
	close IN;

	#`rm $tmpFile $tmpFile2 $tmpFile3`;


	for (my $i=0;$i<@rv;$i++) {
		$r = $rv[$i];
print STDERR "\t$i\t$r\n";
		if ($r < $min) {
			$r = $min;
		}
		$r = log($r)/log(10.0);
		if ($pvalues->[$i] > 0) {
			$r *= -1;
		}
print STDERR "\t$i\t$r\n";
		$rv[$i] = $r;
	}
	return \@rv;
}


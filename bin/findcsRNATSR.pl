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
# GNU General Public License for more details.


use POSIX;
use HomerConfig;
use Statistics;

#set key parameters, read command line arguments

my $config = HomerConfig::loadConfigFile();

my $reformatIDflag = 1;
my $includeSingleExonsFlag = 0;
my $include1stExonFlag = 0;
my $pseudoCount = 1.0;

my $maxCPUs = 3;
my $gtfTssSize = 400;
my $maxRNALog2Fold = 3e10;
my $maxInputLog2Fold = 3e10;
my $minDistDiff = 0.15;
my $defaultLog2Fold = 1.0;
my $localFold = 2.0;
my $tsrSize = 150;
my $minDist = '';
my $freqSize = 2000;
my $tagThreshold = 7;
my $mincsRNAforRatioCal = 7;
my $eps = 1e-7;
my $tbpLimitMetaGene = 3;
my $minMetaGeneSize = 300;
#my $exonFDR = 0.1;
my $senseSizeTSS = 150;
#my $senseUpstreamRNA = "-100,500";
my $senseUpstreamRNA = "-150,-10";
my $senseSizeRNA = "-100,500";
my $antisenseSize = "-500,100";
#$findPeaksTSSoptions = " -size 150 -L 2 -ntagThreshold 15 ";
$findPeaksTSSoptions = "";
my $tpRegions = "";
my $fpRegions = "";
my $metaGeneFlag = 0;
my $rnaseqStrand = '+';

my $reportDistDistance = 500;
my $reportRNAstableThresh = 2.0;
my $reportBidirectionalThresh = 2.0;

sub printCMD {
	print STDERR "\n\tfindcsRNATSR.pl <csRNA tag directory> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-o <prefix> \n";
	print STDERR "\t\t-i <csRNA input tag directory> (*see note below)\n";
	print STDERR "\t\t-rna <RNAseq tag directory> (*see note below)\n";
	#print STDERR "\t\t\t-rnaseqStrand <+|-|both> (strand relative to TSS associated with RNA-seq reads,def:$rnaseqStrand)\n";
	print STDERR "\t\t-genome <genome> (HOMER genome identifier or FASTA file)\n";
	print STDERR "\t\t-gtf <gtf file> (gene annotation, by default will use genome annotation if available)\n";
	print STDERR "\t\t\n";
	print STDERR "\t\t-ntagThreshold <#> (minimum number of norm. tags per TSR, default: $tagThreshold)\n";
	print STDERR "\t\t-size <#> (size of regions to check for significant TSR, default: $tsrSize)\n";
	print STDERR "\t\t-includeSingleExons (include single exon genes for TP/FP calculations, default: only use multiexon genes)\n";
	print STDERR "\t\t-metagene (Include a metagene plot of read distributions)\n";
	print STDERR "\t\t-cpu <#> (max CPUs to use)\n";

	print STDERR "\n\t\tEnrichment threshold options:\n";
	print STDERR "\t\t-maxInputLog2Fold <#> (maximum log2 fold enrichment vs. input or RNA, default: $maxInputLog2Fold)\n";
	print STDERR "\t\t-maxRNALog2Fold <#> (maximum log2 fold enrichment vs. input or RNA, default: $maxRNALog2Fold)\n";
	print STDERR "\t\t-minDistDiff <#> (minimum CDF difference to consider for TP/FP calculation, default: $minDistDiff)\n";
	print STDERR "\t\t-defaultLog2Fold <#> (if annotation-based TP/FP threshold calculation failes, default: $defaultLog2Fold)\n";
	print STDERR "\t\t-pseudoCount <#> (number of normalizaed reads added to input/target read counts, default: 1.0)\n";
	print STDERR "\t\t-ntagThreshold <#> (minimum number of norm. tags per TSS cluster, default: $tagThreshold)\n";
	print STDERR "\t\t-L <#> (local fold enrichment used for initial TSS identification, default: $localFold)\n";
	print STDERR "\t\t-noFilterRNA (only use RNAseq for stability calculations, not for filtering putative TSS)\n";

	print STDERR "\n\t\tRegion sizes:\n";
	print STDERR "\t\t-size <#> (size of regions to check for significant TSRs, default: $tsrSize)\n";
	print STDERR "\t\t-minDist <#> (minimum distance between adjacent TSRs, default: $minDist)\n";
	print STDERR "\t\t-senseRegion <#>,<#> (Size of region to calculate RNA/input read density,default: $senseSizeRNA)\n";
	print STDERR "\t\t-antisenseRegion <#>,<#> (Size of region to calculate RNA/input read density,default: $antisenseSize)\n";
	print STDERR "\t\t-rnaRegion <#>,<#> (Size of region to calculate RNA coverage for FP identification,default: $senseUpstreamRNA)\n";
	#print STDERR "\t\t-fisherExact (Use fisher exact to score enrichment in TSS vs. exons, default: max diff in CDF)\n";
	
	print STDERR "\n\t\tAnnotation options:\n";
	print STDERR "\t\t-stableThreshold <#> (minimum RNAseq read density for stable transcript, def: $reportRNAstableThresh)\n";
	print STDERR "\t\t-biThreshold <#> (minimum divergent csRNAseq signal for bi/divergent TSS, def: $reportBidirectionalThresh)\n";
	print STDERR "\t\t-dist <#> (minimum distance from annotated TSS to assign 'distal' annotation, def: $reportDistDistance)\n";
	print STDERR "\n\t* Tag Directories are assumed to represent the 5' ends of RNA in sense orientation\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 1) {
	printCMD();
}

my $prefix = "out";
my $tagDir = $ARGV[0];
my $inputDir = '';
my $rnaDir = '';
my $gtfFile = '';
my $gtfOptions = '';
my $genome = 'none';
my $fisherFlag = 0;
my $noFilterRNA = 0;
my $cmd = "findcsRNATSS.pl";
for (my $i=0;$i<@ARGV;$i++) {
	$cmd .= " " . $ARGV[$i];
}

for (my $i=1;$i<@ARGV;$i++) {

	if ($ARGV[$i] eq '-i') {
		$inputDir = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rna') {
		$rnaDir = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-o' || $ARGV[$i] eq '-prefix') {
		$prefix = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-gtf') {
		$gtfFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-metagene') {
		$metaGeneFlag = 1;
	} elsif ($ARGV[$i] eq '-minDistDiff') {
		$minDistDiff = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-defaultLog2Fold') {
		$defaultLog2Fold = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxRNALog2Fold') {
		$maxRNALog2Fold = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxInputLog2Fold') {
		$maxInputLog2Fold = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-ntagThreshold') {
		$tagThreshold = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-L') {
		$localFold = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-size') {
		$tsrSize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minDist') {
		$minDist = " -minDist " . $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rnaseqStrand') {
		$rnaseqStrand = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-senseRegion') {
		$senseSizeRNA = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-antisenseRegion') {
		$antisenseSize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rnaRegion') {
		$senseUpstreamRNA = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-fisherExact') {
		$fisherFlag = 1;
	} elsif ($ARGV[$i] eq '-includeSingleExons') {
		$includeSingleExonsFlag = 1;
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pseudoCount') {
		$pseudoCount = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-noFilterRNA') {
		$noFilterRNA = 0;
	} elsif ($ARGV[$i] eq '-stableThreshold') {
		$reportRNAstableThresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-biThreshold') {
		$reportBidirectionalThresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-dist') {
		$reportDistDistance = $ARGV[++$i];
	} else {
		print STDERR "\nERROR: Could not recognize option: \"$ARGV[$i]\"\n\n";
		printCMD();
	}
}

#stuff to work with HOMER's preconfigured genomes... probably best to just use GTF files... 
my $preloadedGenomeFlag = 0;
my $homerRNAfile = '';
my $homerTSSfile = '';
my $homerExonfile = '';
my $genomeDir = '';
if ($gtfFile eq '') {
	if ($genome eq 'none' || $genome eq '') {
		#not a valid genome
	} elsif (-f $genome) {
		#genome is likely a fasta file
	} elsif (exists($config->{'GENOMES'}->{$genome})) {
		$genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'};
		$homerRNAfile = $genomeDir . "/" . $genome . ".rna";
		unless (-f $homerRNAfile) {
			print STDERR "\tError - can't find $genome.rna file in $genomeDir\n";
			$homerRNAfile = '';
		}
		$homerExonfile = $genomeDir . "/annotations/basic/exons.ann.txt";
		unless (-f $homerExonfile) {
			print STDERR "\tError - can't find homer exon file: $homerExonfile\n";
			$homerExonfile = '';
		}
		$homerTSSfile = $genomeDir . "/" . $genome . ".tss";
		unless (-f $homerTSSfile) {
			print STDERR "\tError - can't find homer tss file: $homerTSSfile\n";
			$homerTSSfile = '';
		}
	} else {
		#not a valid genome
		print STDERR "\tWarrning - genome ($genome) is not a valid file\n";
	}
}


my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";
my $tmpFile4 = $rand . ".4.tmp";

my $initialTSSfile = $rand . ".tsr.tmp";
my $gtfTSSfile = $rand . ".gtfTSS.tmp";
my $gtfSingleExonFile = $rand . ".gtfSingleExon.tmp";
my $gtf1stExonFile = $rand . ".gtf1stExon.tmp";
my $gtfOtherExonsFile = $rand . ".gtfOtherExons.tmp";
my $divergentTSSfile = $rand . ".divergentTSS.txt";
my $antisenseTSSfile = $rand . ".antisenseTSS.txt";

if ($gtfFile ne '' || $homerRNAfile ne '') {
	$tpRegions = $rand . ".tpRegions.txt";
	$fpRegions = $rand . ".fpRegions.txt";
}


my $tagDirs = " \"$tagDir\"";
$tagDirs .= " \"$inputDir\"" if ($inputDir ne '');
$tagDirs .= " \"$rnaDir\"" if ($rnaDir ne '');


#Perform initial identification of TSS clusters (no input used)

$findPeaksTSSoptions .= " -L $localFold -size $tsrSize -ntagThreshold $tagThreshold $minDist ";

`findPeaks "$tagDir" -style tsr $findPeaksTSSoptions > $initialTSSfile`;


#From GTF file, identify key regions for QC and parameter optimization - (find annotated TSS, and where
#secondary exons are where we expect to find less TSS (and potentially many false positives).
if ($gtfFile ne '') {
	$gtfOptions = " -gtf \"$gtfFile\" ";
} elsif ($genome ne '' || $genome ne 'none') {

}

#Annotate file with the GTF regions (i.e. 1st exon, etd.), read in original TSS positions
`annotatePeaks.pl $initialTSSfile $genome -strand + -fragLength 1 -cpu $maxCPUs $gtfOptions > "$tmpFile2"`;
my %data = ();
my $c = 0;
open IN, $tmpFile2;
while (<IN>) {
	$c++;
	next if ($c < 2);
	chomp;
	s/\r//g;
	my $og = $_;
	my @line = split /\t/;
	my $id = $line[0];
	my $chr = $line[1];
	my $start = $line[2];
	my $end = $line[3];
	my $strand = $line[4];
	my $score = $line[5];
	my $focusRatio = $line[6];
	my $mid = floor(($start+$end)/2);

	my $printID = $id;
	if ($reformatIDflag) {
		$printID = "TSS-$chr:$start-$end,$strand";
	}
	my $str = "$printID\t$chr\t$start\t$end\t$strand\t$score\t$focusRatio";

	my $ann = 'other';
	my $csRNA = 0;
	my $csRNAinput = 0;
	my $rna =0;

	my $closestDist = 1000000;
	my $closestTrans = "";
	my $closestGene = "";
	my $closestSymbol = "";
	my $closestType = "";
	$closestDist = $line[9] if (@line > 9);
	$closestTrans = $line[10] if (@line > 10);
	$closestGene = $line[11] if (@line > 11);
	$closestSymbol = $line[15] if (@line > 15);
	$closestType = $line[18] if (@line > 18);

	$data{$id} = {str=>$str, csRNA=>$csRNA, csRNAinput=>$csRNAinput, rna=>$rna,ann=>$ann,anngene=>'na',assign=>'',
		closestDist=>$closestDist,closestTrans=>$closestTrans,closestGene=>$closestGene,closestSymbol=>$closestSymbol,
		closestType=>$closestType,p=>$mid,dist2donor=>'',donorSeq=>'',stablerna=>$rna,stablernafold=>'',
		foldRNA=>"na", foldInput=>"na",revrna=>'na',revcsRNA=>'na'};
}
close IN;
`rm $tmpFile2`;


#calculate reads in putative TSS clusters for primary experiment, input, RNA-seq, and annotate each region based
#on the GTF region (i.e. 1st exon, etd.)
`annotatePeaks.pl $initialTSSfile $genome -strand + -fragLength 1 -size $senseSizeTSS -d $tagDirs -cpu $maxCPUs > $tmpFile2`;

$c = 0;
open IN, $tmpFile2;
while (<IN>) {
	$c++;
	next if ($c < 2);
	chomp;
	s/\r//g;
	my $og = $_;
	my @line = split /\t/;
	my $id = $line[0];
	if (!exists($data{$id})) {
		print STDERR "!!! could not find $id\n";
		next;
	}

	if (@line < 20) {
		print STDERR "Warning... not enough columns in file\n";
		next;
	}
	my $csRNA = $line[19];
	my $csRNAinput = 'na';
	my $rna = 'na';

	my $ann = 'other';
	my $col = 20;
	if ($inputDir ne '' && @line > $col) {
		$csRNAinput = $line[$col];
		$col++;
	}
	if ($rnaDir ne '' && @line > $col) {
		$rna = $line[$col];
		$col++;
	}

	$data{$id}->{'csRNA'}=$csRNA;
	$data{$id}->{'csRNAinput'}=$csRNAinput;
	$data{$id}->{'rna'}=$rna;
	$data{$id}->{'stablerna'}=$rna;

	my $foldInput = 'na';
	$foldInput = log(($csRNA+$pseudoCount)/($csRNAinput+$pseudoCount))/log(2.0) if ($csRNAinput ne 'na');
	$data{$id}->{'foldInput'} = $foldInput;
	my $foldRNA = 'na';
	$foldRNA = log(($csRNA+$pseudoCount)/($rna+$pseudoCount))/log(2.0) if ($rna ne 'na');
	$data{$id}->{'foldRNA'} = $foldRNA;
}
close IN;
`rm $tmpFile2`;


if ($gtfFile ne '' || ($homerTSSfile ne '' && $homerRNAfile ne '' && $homerExonfile ne '')) {
	if ($gtfFile ne '') {
		`parseGTF.pl "$gtfFile" tss > $tmpFile`;
	} else {
		`cp "$homerTSSfile" "$tmpFile"`;
	}
	`adjustPeakFile.pl "$tmpFile" -size $gtfTssSize > "$gtfTSSfile"`;
	`adjustPeakFile.pl "$gtfTSSfile" -flipStrand > "$tmpFile"`;
	`adjustPeakFile.pl "$tmpFile" -size 600 > "$antisenseTSSfile"`;

	#my $Y1 = `wc -l $tmpFile $gtfTSSfile $antisenseTSSfile $homerTSSfile $gtfFile`;
	#print STDERR "Y1: $Y1\n";

	if ($gtfFile ne '') {
		`parseGTF.pl "$gtfFile" exons > "$tmpFile"`;
	} else {
		`cp "$homerExonfile" "$tmpFile"`;
	}
	#`cp $tmpFile $tmpFile.save`;
	
	my $grepInputFile1 = $rand . ".grepinput1.txt";
	my $grepInputFile2 = $rand . ".grepinput2.txt";
	my $grepInputFile3 = $rand . ".grepinput3.txt";
	`cp "$tmpFile" "$grepInputFile1"`;

	my $pattern = 'exon 1 of 1\)';
	`egrep "$pattern" "$grepInputFile1" > "$gtfSingleExonFile"`;
	`egrep -v "$pattern" "$grepInputFile1" > "$grepInputFile2"`;
	$pattern = 'exon 1 of ';
	`cp "$grepInputFile2" "$grepInputFile3"`;
	`egrep "$pattern" "$grepInputFile3" > "$gtf1stExonFile"`;
	`egrep -v "$pattern" "$grepInputFile3" > "$gtfOtherExonsFile"`;
	`rm $grepInputFile1 $grepInputFile2 $grepInputFile3`;
	
	#my $xxx = `wc -l $tmpFile $grepInputFile1 $grepInputFile2 $grepInputFile3 $gtfSingleExonFile $gtf1stExonFile $gtfOtherExonsFile`;
	#print STDERR "xxx=\n$xxx\n";

	`adjustPeakFile.pl "$gtf1stExonFile" -5p -size $gtfTssSize > "$tpRegions"`;
	`cp "$gtfOtherExonsFile" "$fpRegions"`;

	my $X1 = `wc -l $tpRegions`;
	my $X2 = `wc -l $fpRegions`;
	my $X3 = `wc -l $tpRegions`;
	my $X4 = `wc -l $fpRegions`;

	if ($includeSingleExonsFlag) {
		`adjustPeakFile.pl $gtfSingleExonFile -5p -size $gtfTssSize >> $tpRegions`;
		`cat $fpRegions $gtfSingleExonFile > $tmpFile`;
		`mv $tmpFile $fpRegions`;
		$X3 = `wc -l $tpRegions`;
		$X4 = `wc -l $fpRegions`;
	}
	#print STDERR "$X1\n$X2\n$X3\n$X4\n";


	`adjustPeakFile.pl $initialTSSfile -flipStrand > $tmpFile`;
	`adjustPeakFile.pl $tmpFile -move 200 -size 600 > $divergentTSSfile`;

	`adjustPeakFile.pl $gtf1stExonFile -3p -size -2,6 > $tmpFile`;
	if ($genomeDir ne '') {
		`homerTools extract $tmpFile $genomeDir > $tmpFile2`;
	} else {
		`homerTools extract $tmpFile $genome > $tmpFile2`;
	}
	my %junc = ();
	open IN, $tmpFile2;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		next if (@line < 2);
		my $id = $line[0];
		if ($id =~ /^exon \((.*?)\, exon 1 of/ || $id =~ /^non-coding \((.*?)\, exon 1 of/) {
			$id = $1;
			my $s = $line[1];
			$junc{$id} = {s=>$s};
		}
	}
	close IN;
	open IN, $gtf1stExonFile;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		next if (@line < 5);
		my $id = $line[0];
		if ($id =~ /^exon \((.*?)\, exon 1 of/ || $id =~ /^non-coding \((.*?)\, exon 1 of/) {
			$id = $1;
			next if (!exists($junc{$id}));
			my $chr = $line[1];
			my $p = $line[3];
			my $d = $line[4];
			if ($d eq '-') {
				$p = $line[2];
			}
			$junc{$id}->{'c'}=$chr;
			$junc{$id}->{'p'}=$p;
			$junc{$id}->{'d'}=$d;
		}
	}
	close IN;

	`rm $tmpFile $tmpFile2`;
	
	`mergePeaks -strand -cobound 1 -prefix $rand $initialTSSfile $gtfTSSfile $gtfSingleExonFile $gtf1stExonFile $gtfOtherExonsFile $divergentTSSfile $antisenseTSSfile`;
	open IN, "$rand.coBoundBy1.txt" or print STDERR "!!! Could not open $rand.coBoundBy1.txt";
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^\#/);
		my @line = split /\t/;
	
		my $id = $line[0];
		next if (!exists($data{$id}));
		my $tssPos = $data{$id}->{'p'};


		my @overlaps = ();
		for (my $i=7;$i<=12;$i++) {
			if ($i >= scalar(@line) || $line[$i] eq '') {
				push(@overlaps,'');
				next;
			}
			my @a = split /\,/, $line[$i];
			$a[0] =~ s/^.*?\(//;
			$a[0] =~ s/\)$//;
			push(@overlaps,$a[0]);
		}

		my $ann = '';
		my $gene= '';
		my $dist2donor = '';
		my $donorSeq = '';
		if ($overlaps[1] ne '') {
			$ann = 'singleExon';
			$gene = $overlaps[1];
		} elsif ($overlaps[0] ne '') {
			$ann = 'tss';
			$gene = $overlaps[0];
			if (exists($junc{$gene})) {
				$donorSeq = $junc{$gene}->{'s'};
				my $p = $junc{$gene}->{'p'};
				my $d = $junc{$gene}->{'d'};
				$dist2donor = $p-$tssPos;
				if ($d eq '-') {
					$dist2donor = $tssPos-$p;
				}
			}
		} else {
			if ($overlaps[3] ne '') {
				if ($overlaps[4] ne '') {
					$ann = 'otherExonBidirectional';
				} else {
					$ann = 'otherExon';
				}
				$gene = $overlaps[3];
			} elsif ($overlaps[2] ne '') {
				$ann = 'firstExon';
				$gene = $overlaps[2];
			} elsif ($overlaps[5] ne '') {
				$ann = 'promoterAntisense';
				$gene = $overlaps[5];
			} else {
				$ann = 'other';
				$gene = '';
			}
		}
		if (exists($data{$id})) {
			$data{$id}->{'ann'} = $ann;
			$data{$id}->{'anngene'} = $gene;
			$data{$id}->{'dist2donor'} = $dist2donor;
			$data{$id}->{'donorSeq'} = $donorSeq;
		}
			
	}
	close IN;

	if ($metaGeneFlag) {
		if ($gtfFile eq '') {
			`cp "$homerRNAfile" "$tmpFile"`;
		} else {
			`parseGTF.pl "$gtfFile" rna > "$tmpFile"`;
		}
		`makeMetaGeneProfile.pl "$tmpFile" none -size 5000 -bin 10 -min $minMetaGeneSize -gbin 100 -pc $tbpLimitMetaGene -d $tagDir $inputDir $rnaDir > "$prefix.metagene.txt"`;
		`rm "$tmpFile"`;
	}


	`rm "$rand.coBoundBy0.txt" "$rand.coBoundBy1.txt"`;
	`rm $gtfTSSfile $gtfSingleExonFile $gtf1stExonFile $gtfOtherExonsFile $divergentTSSfile $antisenseTSSfile`;
}

my $inputThresh = -1e10;
my $inputMaxDiff = -1e10;
my $rnaThresh = -1e10;
my $rnaMaxDiff = -1e10;
my $numTPinput = 0;
my $numTPrna = 0;
my $numFPinput = 0;
my $numFPrna = 0;
if ($fpRegions ne '' && $tpRegions ne '') {

	my %tpfpdata = ();
	my $halfGtfTssSize = floor($gtfTssSize/2);
	#recenter regions on the major TSS (just like the putative TSS
	`getPeakTags "$tpRegions" "$tagDir" -center -fragLength 1 -start -$halfGtfTssSize -end $halfGtfTssSize -strand + -tagAdjust 0 > "$tmpFile"`;
	`mv "$tmpFile" "$tpRegions"`;
	`getPeakTags "$fpRegions" "$tagDir" -center -fragLength 1 -start -$halfGtfTssSize -end $halfGtfTssSize -strand + -tagAdjust 0 > "$tmpFile"`;
	`mv "$tmpFile" "$fpRegions"`;

	#check for overlap between false positive regions and true positive regions, only considering FPs that don't overlap with TP
	#also discard regions that are potentially antisense to putative TSS clusters (i.e. might be bidirectional)
	`adjustPeakFile.pl "$initialTSSfile" -move -150 -flipStrand > "$tmpFile"`;
	`mergePeaks "$fpRegions" "$tpRegions" "$tmpFile" -cobound 1 -strand -d $gtfTssSize -prefix "$rand"`;

	my $numInitialFp = 0;
	my $numInitialTp = 0;
	open IN, "$rand.coBoundBy0.txt" or print STDERR "!! Warning: Couldn't open $rand.coBoundBy0.txt\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		my @line = split /\t/;
		$numInitialFp++;
		$tpfpdata{$line[0]} = {csRNA=>'na',csRNAinput=>'na',rna=>'na',assign=>'fp'};
	}
	close IN;
	open IN, $tpRegions or print STDERR "!! Warning: Couldn't open tpRegions ($tpRegions)\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		my @line = split /\t/;
		$numInitialTp++;
		$tpfpdata{$line[0]} = {csRNA=>'na',csRNAinput=>'na',rna=>'na',assign=>'tp'};
	}
	close IN;
	`rm "$rand.coBoundBy0.txt" "$rand.coBoundBy1.txt"`;
	
	`cat "$tpRegions" "$fpRegions" > "$tmpFile"`;
	`annotatePeaks.pl "$tmpFile" none -strand + -fragLength 1 -size $senseSizeTSS -d $tagDirs -cpu $maxCPUs > "$tmpFile2"`;
	$c = 0;

	open IN, $tmpFile2;
	while (<IN>) {
		$c++;
		next if ($c < 2);
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $line[0];
		next if (!exists($tpfpdata{$id}));
		if (@line < 20) {
			print STDERR "Warning... not enough columns in file\n";
			next;
		}
		my $csRNA = $line[19];
		my $csRNAinput = 'na';
		my $rna = 'na';

		my $ann = 'other';
		my $col = 20;
		if ($inputDir ne '' && @line > $col) {
			$csRNAinput = $line[$col];
			$col++;
		}
		if ($rnaDir ne '' && @line > $col) {
			$rna = $line[$col];
			$col++;
		}
		$tpfpdata{$id}->{'csRNA'}=$csRNA;
		$tpfpdata{$id}->{'csRNAinput'}=$csRNAinput;
		$tpfpdata{$id}->{'rna'}=$rna;
	}
	close IN;
	`rm "$tmpFile2"`;

	if ($rnaDir ne '') {
		#quantify RNA just upstream of TSS - signal here usually points to a false positive TSS from a 
		#highly expressed gene. 
		`adjustPeakFile.pl "$tmpFile" -size $senseUpstreamRNA > "$tmpFile3"`;
		`annotatePeaks.pl "$tmpFile3" none -strand + -fragLength 1 -d $tagDirs -cpu $maxCPUs > "$tmpFile2"`;
		$c = 0;
		open IN, $tmpFile2;
		while (<IN>) {
			$c++;
			next if ($c < 2);
			chomp;
			s/\r//g;
			my @line = split /\t/;
			my $id = $line[0];
			next if (!exists($tpfpdata{$id}));
			if (@line < 20) {
				print STDERR "Warning... not enough columns in file\n";
				next;
			}
			my $csRNA = $line[19];
			my $csRNAinput = 'na';
			my $rna = 'na';
			my $col = 20;
			if ($inputDir ne '' && @line > $col) {
				$csRNAinput = $line[$col];
				$col++;
			}
			if ($rnaDir ne '' && @line > $col) {
				$rna = $line[$col];
				$col++;
			}
			$tpfpdata{$id}->{'rna'} = $rna;
		}
		close IN;
		`rm "$tmpFile2" "$tmpFile3"`;
	}
	`rm "$tmpFile" "$fpRegions" "$tpRegions"`;


	# To establish fold-enrichment threshold for input and RNA vs. the true csRNA/Start/5'GRO/GROcap signal,
	# put putative TSS clusters into likely true-positive and false positive groups based on their
	# annotation near known TSS or in downstream exons.

	
	my @tssIDs = keys %tpfpdata;	
	my @tss = ();
	my @exons = ();
	my $totaltpfpdata = 0;
	foreach(@tssIDs) {
		my $id = $_;
		my $csRNA = $tpfpdata{$id}->{'csRNA'};
		my $csRNAinput = $tpfpdata{$id}->{'csRNAinput'};
		$totaltpfpdata++;
		next if ($csRNAinput eq 'na');
		my $foldInput = log(($csRNA+$pseudoCount)/($csRNAinput+$pseudoCount))/log(2.0);
		$tpfpdata{$id}->{'foldInput'} = $foldInput;
		next if ($csRNA+$csRNAinput < $mincsRNAforRatioCal);
		if ($tpfpdata{$id}->{'assign'} eq 'tp') {
			push(@tss, $foldInput);
		} elsif ($tpfpdata{$id}->{'assign'} eq 'fp') {
			push(@exons, $foldInput);
		}
	}


	#create differences in log2fold change between csRNA and input (or RNA) and calculate
	#threshold that results in the greatest difference in distributions (i.e. threshold
	#that distriminates between TSS and 'other exon' annotation, which are more likely 
	#to be noise.
	@tss = sort {$b <=> $a} @tss;
	@exons = sort {$b <=> $a} @exons;


	my $j=0;
	my $maxDiff = -1;
	my $minP = 1;
	my $thresh = 0;
	my $Nx = scalar(@tss);
	my $Ny = scalar(@exons);
	$numTPinput = $Nx;
	$numFPinput = $Ny;


	if ($Nx < 1 || $Ny < 1) {
		if ($inputDir ne '') {
			print STDERR "Warning - no TSS(true-positives) and/or Exons(false-positives) found (analysis of $tagDir(csRNA) and $inputDir(input))\n";
			print STDERR "\tTotal tpfpdata=$totaltpfpdata\n";
			print STDERR "\ttp=$numInitialTp fp=$numInitialFp\n";
			print STDERR "\ttpregions = $tpRegions, $fpRegions\n";
		}
	} else { 
		open INPUT, ">$prefix.inputDistribution.txt";
		print INPUT "csRNA/input log2 ratio\tTSS CDF\tExon CDF\tDifference\tHGpvalue\n";
		my $N = $Nx+$Ny;
		for (my $i=0;$i<@tss;$i++) {
			next if ($i+1<@tss && $tss[$i]-$tss[$i+1] < $eps); 
			while ($j < $Ny && $exons[$j]+$eps > $tss[$i]) {
				$j++;
			}
			my $fracX = ($i+1)/$Nx;
			my $fracY = ($j)/$Ny;
			my $diff = $fracX-$fracY;
			my $fdr = ($j+1)/($i+$j+1);
			my $p = 'NA';
			if ($fisherFlag) {
				$p = Statistics::loghypergeo($N,$Nx,$i+1+$j,$i+1);
				if ($p < $minP) {
					$minP = $p;
					$thresh = $tss[$i];
				}
			} else {
				if ($diff > $maxDiff) {
					$maxDiff = $diff;
					$thresh = $tss[$i];
				}
			}

			print INPUT "$tss[$i]\t$fracX\t$fracY\t$diff\t$p\n";
		}
		close INPUT;
		print STDERR "\n\tMaxDiff(input) = $maxDiff\n";
		if ($fisherFlag) {
			print STDERR "\n\tMinP(input) = $minP\n";
		}
		print STDERR "\tThreshold(input) = $thresh\n\n";
		if ($maxDiff < $minDistDiff) {
			print STDERR "\tDifference in distributions is very small - be careful! (vs. input)\n";
			#print STDERR "\tDifference in distributions is too small: setting fold enrichment to default\n";
			#$thresh = $defaultLog2Fold;
		}
		$inputThresh = $thresh;
		$inputMaxDiff = $maxDiff;
		$inputMaxDiff = $minP if ($fisherFlag);
	}
	if ($inputThresh > $maxInputLog2Fold) {
		print STDERR "\t!!!Warning: Automatic Log2 Fold enrichment vs. input ($inputThresh) greater than max ($maxInputLog2Fold)\n";
		print STDERR "\t\tSometimes this can happen if the data is very clean where most of the initial peaks are legit TSS\n";
		$inputThresh = $maxInputLog2Fold;
	}


	# Now do the same for RNAseq
	@tss = ();
	@exons = ();
	foreach(@tssIDs) {
		my $id = $_;
		my $csRNA = $tpfpdata{$id}->{'csRNA'};
		my $rna = $tpfpdata{$id}->{'rna'};
		next if ($rna eq 'na');
		my $foldRNA = log(($csRNA+$pseudoCount)/($rna+$pseudoCount))/log(2.0);
		$tpfpdata{$id}->{'foldRNA'} = $foldRNA;
		next if ($csRNA+$rna < $mincsRNAforRatioCal);
		if ($tpfpdata{$id}->{'assign'} eq 'tp') {
			push(@tss, $foldRNA);
		} elsif ($tpfpdata{$id}->{'assign'} eq 'fp') {
			push(@exons, $foldRNA);
		}
	}
	@tss = sort {$b <=> $a} @tss;
	@exons = sort {$b <=> $a} @exons;

	$j=0;
	$maxDiff = -1;
	$minP = 1;
	$thresh = 0;
	$Nx = scalar(@tss);
	$Ny = scalar(@exons);
	$numTPrna = $Nx;
	$numFPrna = $Ny;
	if ($Nx < 1 || $Ny < 1) {
		print STDERR "Warning - no TSS(true-positives) and/or Exons(false-positives) found (analysis of $tagDir(csRNA) and $rnaDir(RNAseq))\n" if ($rnaDir ne '');
	} else { 
		open RNA, ">$prefix.rnaDistribution.txt";
		print RNA "csRNA/RNA log2 ratio\tTSS CDF\tExon CDF\tDifference\tHGpvalue\n";
		my $N = $Nx+$Ny;
		for (my $i=0;$i<@tss;$i++) {
			next if ($i+1<@tss && $tss[$i]-$tss[$i+1] < $eps); 
			while ($j < $Ny && $exons[$j]+$eps >= $tss[$i]) {
				$j++;
			}
			my $fracX = ($i+1)/$Nx;
			my $fracY = ($j)/$Ny;
			my $diff = $fracX-$fracY;
			my $fdr = ($j+1)/($i+$j+1);
			my $p = 'NA';
			if ($fisherFlag) {
				$p = Statistics::loghypergeo($N,$Nx,$i+1+$j,$i+1);
				if ($p < $minP) {
					$minP = $p;
					$thresh = $tss[$i];
				}
			} else {
				if ($diff > $maxDiff) {
					$maxDiff = $diff;
					$thresh = $tss[$i];
				}
			}

			print RNA "$tss[$i]\t$fracX\t$fracY\t$diff\t$p\n";
		}
		close RNA;
		print STDERR "\n\tMaxDiff(rna) = $maxDiff\n";
		print STDERR "\n\tMinPvalue(rna) = $minP\n" if ($fisherFlag);
		print STDERR "\tThreshold(rna) = $thresh\n\n";
		if ($maxDiff < $minDistDiff) {
			print STDERR "\tDifference in distributions is very small - be careful! (vs. rna)\n";
			#print STDERR "\tDifference in distributions is too small: setting fold enrichment to default (rna)\n";
			#$thresh = $defaultLog2Fold;
		}
		$rnaThresh = $thresh;
		$rnaMaxDiff = $maxDiff;
		$rnaMaxDiff = $minP if ($fisherFlag);
	}
	if ($rnaThresh > $maxRNALog2Fold) {
		print STDERR "\t!!!Warning: Automatic Log2 Fold enrichment vs. rnaseq ($rnaThresh) greater than max ($maxRNALog2Fold)\n";
		print STDERR "\t\tSometimes this can happen if the data is very clean where most of the initial peaks are legit TSS\n";
		$rnaThresh = $maxRNALog2Fold;
	}

}
if ($inputThresh < -1e9 && $inputDir ne '') {
	print STDERR "\tUsing the default input threshold ($defaultLog2Fold)\n";
	$inputThresh = $defaultLog2Fold;
}
if ($rnaThresh < -1e9 && $rnaDir ne '') {
	print STDERR "\tUsing the default rnaseq threshold ($defaultLog2Fold)\n";
	$rnaThresh = $defaultLog2Fold;
}


if ($rnaDir ne '') {

	#quantify RNA just upstream of TSS - signal here usually points to a false positive TSS from a 
	#highly expressed gene. 
	`adjustPeakFile.pl $initialTSSfile -size $senseUpstreamRNA > $tmpFile`;
	`annotatePeaks.pl $tmpFile none -strand + -fragLength 1 -d $tagDirs -cpu $maxCPUs > $tmpFile2`;
	$c = 0;
	open IN, $tmpFile2;
	while (<IN>) {
		$c++;
		next if ($c < 2);
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $line[0];
		next if (!exists($data{$id}));
	
		if (@line < 20) {
			print STDERR "Warning... not enough columns in file\n";
			next;
		}
		my $csRNA = $data{$id}->{'csRNA'};
		my $csRNAinput = 'na';
		my $rna = 'na';
		my $ann = '';
		my $col = 20;
		if ($inputDir ne '' && @line > $col) {
			$csRNAinput = $line[$col];
			$col++;
		}
		if ($rnaDir ne '' && @line > $col) {
			$rna = $line[$col];
			$col++;
		}
		$data{$id}->{'rna'} = $rna;
		my $foldRNA = 'na';
		$foldRNA = log(($csRNA+$pseudoCount)/($rna+$pseudoCount))/log(2.0) if ($rna ne 'na');
		$data{$id}->{'foldRNA'} = $foldRNA;
		#print STDERR "$id\t$rna\t$foldRNA\n";
	}


	#quantify stable RNA if available
	`adjustPeakFile.pl $initialTSSfile -size $senseSizeRNA > $tmpFile`;
	`annotatePeaks.pl $tmpFile none -strand + -fragLength 1 -d $tagDirs -cpu $maxCPUs > $tmpFile2`;

	$c = 0;
	open IN, $tmpFile2;
	while (<IN>) {
		$c++;
		next if ($c < 2);
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $line[0];
	
		if (@line < 20) {
			print STDERR "Warning... not enough columns in file\n";
			next;
		}
		my $csRNA = $line[19];
		my $csRNAinput = 'na';
		my $rna = 'na';
		my $ann = '';
		my $col = 20;
		if ($inputDir ne '' && @line > $col) {
			$csRNAinput = $line[$col];
			$col++;
		}
		if ($rnaDir ne '' && @line > $col) {
			$rna = $line[$col];
			$col++;
		}
		my $logRatio = log(($rna+$pseudoCount)/($csRNA+$pseudoCount))/log(2.0);
		$data{$id}->{'stablerna'} = $rna;
		$data{$id}->{'stablernafold'} = $logRatio;
	}
	close IN;
	`rm $tmpFile $tmpFile2`;
}

#quantify negative strand signal in up-stream region (i.e. bidirectional transcription)
`adjustPeakFile.pl $initialTSSfile -flipStrand -size $antisenseSize > $tmpFile`;
`annotatePeaks.pl $tmpFile none -strand + -fragLength 1 -d $tagDirs -cpu $maxCPUs > $tmpFile2`;

$c = 0;
open IN, $tmpFile2;
while (<IN>) {
	$c++;
	next if ($c < 2);
	chomp;
	s/\r//g;
	my $og = $_;
	my @line = split /\t/;
	my $id = $line[0];

	if (@line < 20) {
		print STDERR "Warning... not enough columns in file\n";
		next;
	}
	my $csRNA = $line[19];
	my $csRNAinput = 'na';
	my $rna = 'na';
	my $ann = '';
	my $col = 20;
	if ($inputDir ne '' && @line > $col) {
		$csRNAinput = $line[$col];
		$col++;
	}
	if ($rnaDir ne '' && @line > $col) {
		$rna = $line[$col];
		$col++;
	}
	$data{$id}->{'revcsRNA'} = $csRNA;
	$data{$id}->{'revcsRNAinput'} = $csRNAinput;
	$data{$id}->{'revrna'} = $rna;
	my $logRatio = 'na';
	my $ogLevel = $data{$id}->{'csRNA'};
	if ($ogLevel ne 'na' && $csRNA ne 'na') {
		$logRatio = log(($ogLevel+$pseudoCount)/($csRNA+$pseudoCount))/log(2.0);
	}
   	$data{$id}->{'revfold'} = $logRatio;
}
close IN;
`rm $tmpFile $tmpFile2`;

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if ($noFilterRNA) {
	$noFilterRNA = $rnaThresh;
	$rnaThresh = -1e9;
}


open TSS, ">$prefix.tsr.txt";
open ALLTSS, ">$prefix.alltsr.txt";

my $total = 0;
my $Ntsr = 0;
my $Noverinput = 0;
my $Noverrna = 0;
my $Noverboth = 0;
my $Ninput = 0;
my $Nrna = 0;
my $Nann = 0;
my $Nexon = 0;
my $nDistal = 0;
my $nStable = 0;
my $nBidirectional = 0;
my %recodes=();

my $header = "#tsrID\tchr\tstart\tend\tstrand\tscore\tfocusRatio\tcsRNA\tcsRNAinput\tRNAseq\tannotation";
$header .= "\tgene associated with annotation\tLog2Ratio vs. Input\tLog2Ratio vs. RNAseq";
$header .= "\tstableRNA\tLog2Ratio vs. stableRNA";
$header .= "\trev-csRNA\trev-csRNAinput\trev-RNA\tBidirectionalRatio";
$header .= "\tClosest TSS Symbol\tClosest TSS distance\tClosest TSS transcript\tClosest TSS gene\tClosest TSS type";
$header .= "\tTSS status\tPromoter Proximal/Distal\tStable/Unstable\tBidirectional\tRegulatory Element Configuration";
$header .= "\tDistance to 1st Splice Donor\t1st Splice Donor Sequence";
$header .= "\n";
print TSS $header;
print ALLTSS $header;
foreach(keys %data) {
	$total++;
	my $id = $_;
	my $inputFold = $data{$id}->{'foldInput'};
	my $rnaFold = $data{$id}->{'foldRNA'};
	my $good = 1;
	my $status = 'FP';
	if ($inputFold eq 'na' || $inputFold > $inputThresh) {
		$Noverinput++;
	} else {
		$good = 0;
		$status .= "-overInput";
	}
	if ($rnaFold eq 'na' || $rnaFold > $rnaThresh) {
		$Noverrna++;
	} else {
		$good = 0;
		$status .= "-overRNA";
	}
	if ($good) {
		$Noverboth++;
		$status = 'valid';
	}

	my $str = $data{$id}->{'str'} . "\t" . $data{$id}->{'csRNA'} . "\t" . $data{$id}->{'csRNAinput'};
	$str .= "\t" . $data{$id}->{'rna'} . "\t" . $data{$id}->{'ann'} . "\t" . $data{$id}->{'anngene'};
	$str .= "\t" . $data{$id}->{'foldInput'};
	$str .= "\t" . $data{$id}->{'foldRNA'};
	$str .= "\t" . $data{$id}->{'stablerna'}; 
	$str .= "\t" . $data{$id}->{'stablernafold'}; 
	if (exists($data{$id}->{'revfold'})) {
		$str .= "\t" . $data{$id}->{'revcsRNA'} . "\t" . $data{$id}->{'revcsRNAinput'} . "\t" . $data{$id}->{'revrna'};
		$str .= "\t" . $data{$id}->{'revfold'};
	} else {
		$str .= "\tna\tna\tna\tna";
	}
	$str .= "\t" . $data{$id}->{'closestSymbol'};
	$str .= "\t" . $data{$id}->{'closestDist'};
	$str .= "\t" . $data{$id}->{'closestGene'};
	$str .= "\t" . $data{$id}->{'closestTrans'};
	$str .= "\t" . $data{$id}->{'closestType'};
	$str .= "\t" . $status;

	my $prox = 'na';
	my $stable = 'na';
	my $bi = 'na';
	my $re = 'na';
	
	if ($data{$id}->{'closestDist'} eq 'NA' || abs($data{$id}->{'closestDist'}) > $reportDistDistance) {
		$prox='distal';
	} else {
		$prox='proximal';
	}
	if ($data{$id}->{'stablerna'} ne 'na' && $data{$id}->{'stablerna'} > $reportRNAstableThresh) {
		$stable = 'stable';
		$re = 'S';
	} else {
		$stable = 'unstable';
		$re = 'U';
	}
	if ($data{$id}->{'revcsRNA'} ne 'na' && $data{$id}->{'revcsRNA'} > $reportBidirectionalThresh) {
		$bi = 'bidirectional';
		if ($data{$id}->{'revrna'} ne 'na' && $data{$id}->{'revrna'} > $reportRNAstableThresh) {
			$re .= 'S';
		} else {
			$re .= 'U';
		}
	} else {
		$bi = 'unidirectional';
	}
	
	$str .= "\t" . $prox;
	$str .= "\t" . $stable;
	$str .= "\t" . $bi;
	$str .= "\t" . $re;
	$str .= "\t" . $data{$id}->{'dist2donor'};
	$str .= "\t" . $data{$id}->{'donorSeq'};
	$str .= "\n";

	print ALLTSS $str;
	if ($good) {
		print TSS $str;
		if ($prox eq 'distal') {
			$nDistal++;
		} 
		if ($stable eq 'stable') {
			$nStable++;
		}
		if ($bi eq 'bidirectional') {
			$nBidirectional++;
		}
		$recodes{$re}++;
	}

}	
close TSS;
close ALLTSS;


if ($noFilterRNA ne '0') {
	$rnaThresh = $noFilterRNA;
	$noFilterRNA = 1;
}

#output stats file with parameters and QC numbers
open OUT, ">$prefix.stats.txt";
print OUT "cmd = $cmd\n\n";

my ($csTagTotal, $csPosTotal) = HomerConfig::readTagInfo($tagDir,0);
print OUT "Total csRNA reads: $csTagTotal ($tagDir)\n";
if ($inputDir ne '') {
	my ($inputTagTotal, $inputPosTotal) = HomerConfig::readTagInfo($inputDir,0);
	print OUT "Total input reads: $inputTagTotal ($inputDir)\n";
}
if ($rnaDir ne '') {
	my ($rnaTagTotal, $rnaPosTotal) = HomerConfig::readTagInfo($rnaDir,0);
	print OUT "Total rna reads: $rnaTagTotal ($rnaDir)\n";
}

$Ntss = $Noverboth;
print OUT "\ntotal putative TSS clusters\t$total\n";
print OUT "Valid TSRs\t$Ntss\n";
print OUT "\ttotal over input\t$Noverinput\n";
print OUT "\ttotal over rna  \t$Noverrna\n";
print OUT "\ttotal over both \t$Noverboth\n";

$r = sprintf("%.2lf", $nDistal/$Ntss*100);
print OUT "\nFraction Promoter-Distal TSRs: $r%\n";
$r = sprintf("%.2lf", $nStable/$Ntss*100);
print OUT "Fraction of stable transcript TSRs: $r%\n";
$r = sprintf("%.2lf", $nBidirectional/$Ntss*100);
print OUT "Fraction of bidirectional TSS clusters: $r%\n";

print OUT "Regulatory Element Configurations of each TSR:\n";
my @codes = ('SS','SU','S','US','UU','U');
foreach(@codes) {
	my $v = 0;
	if (exists($recodes{$_})) {
		$v = $recodes{$_};
	}
	$r = sprintf("%.2lf", $v/$Ntss*100);
	print OUT "\t$_: $v ($r%)\n";
}

print OUT "\nvs. Input (-i):\n";
print OUT "\tlog2 fold vs. input: $inputThresh\n";
print OUT "\tMaximum CDF difference: $inputMaxDiff\n";
print OUT "\tTotal TP (tss) regions: $numTPinput\n";
print OUT "\tTotal FP (exon) regions: $numFPinput\n";
print OUT "\nvs. RNA-seq (-rna):\n";
print OUT "\tlog2 fold vs. rna: $rnaThresh\n";
print OUT "\tMaximum CDF difference: $rnaMaxDiff\n";
print OUT "\tTotal TP (tss) regions: $numTPrna\n";
print OUT "\tTotal FP (exon) regions: $numFPrna\n";
if ($noFilterRNA) {
	print OUT "\tRNA-seq was not used to filter putative TSS clusters\n";
}
print OUT "Normalized tag threshold used for initial TSS cluster calls: $tagThreshold\n";
print OUT "\nRead threshold for input/RNA distribution: $mincsRNAforRatioCal\n";
print OUT "Thresholds for enrichment over input and rna determined by";
if ($fisherFlag) {
	print OUT " fisher exact test\n";
} else {
	print OUT " max CDF difference\n";
}
print OUT "PseudoCount for fold change stats: $pseudoCount\n";
print OUT "Size of TSS region for read counting: $senseSizeTSS\n";
print OUT "Size of region to assess RNAseq contamiation: $senseUpstreamRNA\n";
print OUT "Size of region to assess RNAseq stability: $senseSizeRNA\n";
print OUT "Size of region to assess antisense levels: $antisenseSize\n";
print OUT "Distance from annotated promoter considered distal: $reportDistDistance\n";
print OUT "Normalized RNAseq signal to be considered stable: $reportRNAstableThresh\n";
print OUT "Normalized antisense signal to be considered bidirectional: $reportBidirectionalThresh\n";
print OUT "\nRandom number used: $rand\n";

print OUT "\n";
close OUT;

open IN, "$prefix.stats.txt";
while (<IN>) {
	print STDERR "$_";
}
close IN;

my @files = (
	"$prefix.tsr.txt"
	#"$prefix.input.txt",
	#"$prefix.rna.txt",
	#"$prefix.anntss.txt",
	#"$prefix.annexons.txt"
);
if ($genome ne '' && $genome ne 'none') {
	foreach(@files) {
		my $f = $_;
		`annotatePeaks.pl $f $genome -size $freqSize -hist 1 -di > $f.freq.tsv 2> /dev/null`;
	}
}

`rm $initialTSSfile`;
#`rm $fpRegions $tpRegions`;
if ($gtfFile ne '') {
}

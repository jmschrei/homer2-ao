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

my $stitchFlag = 0;
my $stitchMinLen = 100000;
my $stitchMinPos = 0.0;
my $stitchMaxNeg = 0.0;
my $stitchBridgeLen = 50000;
my $normOnly = 0;

sub printCMD {
	print STDERR "\n\tThis program uses DESeq2/edgeR to find differential expression between sets of genes\n";
	print STDERR "\t(R must be installed in the executable path, and the DESeq2/edgeR package must be installed)\n";
	print STDERR "\t   Step 1:  Run analyzeRepeats.pl, but use -raw (or analyzeRNA.pl or annotatePeaks.pl)\n";
	print STDERR "\t   Step 2:  Run this program using that file (use -repeats/-rna/-peaks to match program)\n";
	print STDERR "\tThe output is sent to stdout - appends columns to original file containing data\n";
	print STDERR "\tGroup/Batch names below correspond to order of experiments given as tag directories\n";
	print STDERR "\n\tUsage: getDiffExpression.pl <data file> <group code...> [options]\n";
	print STDERR "\n\t\ti.e. getDiffExpression.pl inputFile.txt groupName1 groupName1 groupName2 > output.txt\n";
	print STDERR "\t\t\t(for 2 replicates in first group and 1 in second group)\n";
	print STDERR "\n\tTo normalize batch effects, add \"-batch <batch code 1> [batch code2] ...\" to the command\n";
	print STDERR "\t\tFor now batch normalization is only available with edgeR and DESeq2\n";
	print STDERR "\t\ti.e. getDiffExpression.pl inputFile.txt gName1 gName1 gName2 gName2 -batch 1 2 1 2 > output.txt\n";
	print STDERR "\t\t\t(for 2 replicates in each group, where the 1st in each and 2nd in are paired)\n";
	print STDERR "\n\t\tUpdate: Expression values are now variance stabalized with DESeq2/rlog by default\n";
	print STDERR "\n\tInput File format options (will try to autodetect):\n";
	print STDERR "\t\t-rna (for analyzeRNA.pl formatted input, default)\n";
	print STDERR "\t\t-repeats (for analyzeRepeats.pl formatted input file)\n";
	print STDERR "\t\t-peaks (for annotatePeaks.pl formatted input file)\n";
	print STDERR "\t\t-loop or -tad (for output files from findTADsAndLoops.pl)\n";
	print STDERR "\t\t-pc1 (for output files from annotatePeaks.pl looking at PC1 bedGraphs, try -pc1options)\n";
	print STDERR "\t\t-basic (for simple file with one column of gene identifiers and then the count data)\n";
	print STDERR "\n\tAdditional Options:\n";
	print STDERR "\t\t-dispersion <#> (edgeR common dispersion to use if no replicates, default: 0.05)\n";
	print STDERR "\t\t-norm2total (normalize using tag directory totals, default: normalize to gene totals(i.e.table)\n";
	print STDERR "\t\t-AvsA (compare each group vs. each group, default: compare 1st group vs. others)\n";
	print STDERR "\t\t-showR (Show R status updates, command output)\n";
	print STDERR "\n\tDifferential Expression program selection:\n";
	print STDERR "\t\t-DESeq2 (use DESeq2, now the default)\n";
	print STDERR "\t\t-DESeq (use DESeq instead of DESeq2 - doesn't support batch mode)\n";
	print STDERR "\t\t-edgeR (use edgeR, - doesn't support -norm2total for normalization to total mapped reads)\n";
	print STDERR "\t\t-limma (use limma, set to work on linear, non-integer values i.e. TAD IRs, PC1 values, etc.)\n";
	print STDERR "\n\tDE List Reporting:\n";
	print STDERR "\t\t-export <prefix> (output differential expression gene lists with filename prefix)\n";
	print STDERR "\t\t-fdr <#> (FDR threshold for diff. expression reporting, default: 0.05)\n";
	print STDERR "\t\t-pvalue <#> (p-value threshold for diff. expression reporting, default: uses -fdr)\n";
	print STDERR "\t\t-log2fold <#> (Log2 Fold threshold for diff. expression reporting, default: 1.0)\n";
	print STDERR "\n\tExpression value normalization/Variance Stabilization:\n";
	print STDERR "\t\t-rlog (Use DESeq2's rlog transform to normalize output matrix, default)\n";
	print STDERR "\t\t-vst (Use DESeq2's vst transform, quicker for large sample sets)\n";
	print STDERR "\t\t-simpleNorm (Normalize to experiment totals, i.e. basic normalization)\n"; 
	print STDERR "\t\t-raw (do not normalize read counts in output table, keep raw values)\n";
	print STDERR "\t\t-normOnly (only normalize data, don't find differentially expressed features)\n";
	print STDERR "\n";
	print STDERR "\tCurrent versions that work:\n";
	print STDERR "\t\tR: v3.3.2\n";
	#print STDERR "\t\tBioconductor: v3.2\n";
	print STDERR "\t\tDESeq v1.24.0\n";
	print STDERR "\t\tDESeq2 v1.14.1\n";
	print STDERR "\t\tedgeR v3.16.5\n";
	print STDERR "\t\tlimma v3.30.13\n";
	print STDERR "\n";
	exit;
}
sub printCMDpc1 {
	print STDERR "\n\tPC1 Hi-C analysis options for stitching differentially scoring regions:\n";
	print STDERR "\t\t-minRegionSize <#> (minimum size of stitched region to report, def: $stitchMinLen)\n";
	print STDERR "\t\t-bridgeDistance <#> (maximum distance to merge sig. regions, def: $stitchBridgeLen)\n";
	print STDERR "\t\t-minPosPC1 <#> (Minimum value for A compartment PC1 values, def: $stitchMinPos)\n";
	print STDERR "\t\t-maxNegPC1 <#> (Maximum value for B compartment PC1 values, def: $stitchMaxNeg)\n";
	print STDERR "\n";
	exit;
}
if (@ARGV < 3) {
	printCMD();
}


my $cmd = "getDiffExpression.pl";
for (my $i=0;$i<@ARGV;$i++) {
	$cmd .= " $ARGV[$i]";
}

#annotatePeaks
#my $col = 18;
my $normFlag = 1;
my $showRFlag = 0;
#annotateRNA.pl
my $col = -1;
my $commonDispersion = 0.05;
my $tagDirNormFlag = 0;
my $AvsAflag = 0;
my $tadFlag = 0;
my $cpFlag = 0;
my %stitchData = ();
my %chrAssignment = ();

my %og = ();
my @header = ();
my @names = ();
my %data = ();

my $mode = '';
#my $mode = 'DESeq2';
my %groups= ();
my %batchs= ();
my @groups = ();
my @batchs = ();
my @argvGroups = ();
my @argvBatchs = ();
my $batchID = 1;
my @tagDirTotals = ();
my $fdrCutoff=0.05;
my $pvalueCutoff='';
my $log2FoldCutoff=1.0;
my $export = "";
my $normMethod = '';



my $rand = rand();

my $batchFlag = 0;
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-batch') {
		$batchFlag = 1;
		next;
	}
	if ($ARGV[$i] eq '-export') {
		$export = $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-AvsA') {
		$AvsAflag = 1;
		next;
	}
	if ($ARGV[$i] eq '-rlog') {
		$normMethod = 'rlog';
		next;
	}
	if ($ARGV[$i] eq '-raw') {
		$normMethod = 'raw';
		next;
	}
	if ($ARGV[$i] eq '-vst') {
		$normMethod = 'vst';
		next;
	}
	if ($ARGV[$i] eq '-simpleNorm') {
		$normMethod = 'simple';
		next;
	}
	if ($ARGV[$i] eq '-normOnly') {
		$normOnly = 1;
		next;
	}
	if ($ARGV[$i] eq '-edgeR') {
		$mode = 'edgeR';
		next;
	}
	if ($ARGV[$i] eq '-limma') {
		$mode = 'limma';
		next;
	}
	if ($ARGV[$i] eq '-showR') {
		$showRFlag = 1;
		next;
	}
	if ($ARGV[$i] eq '-pc1options') {
		printCMDpc1();
		next;
	}
	if ($ARGV[$i] eq '-DESeq') {
		$mode = 'DESeq';
		next;
	}
	if ($ARGV[$i] eq '-minRegionSize') {
		$stitchMinLen = $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-bridgeDistance') {
		$stitchBridgeLen = $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-minPosPC1') {
		$stitchMinPos = $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-maxNegPC1') {
		$stitchMaxNeg = $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-DESeq2') {
		$mode = 'DESeq2';
		next;
	}
	if ($ARGV[$i] eq '-rna') {
		$col = 15;
		next;
	}
	if ($ARGV[$i] eq '-repeats') {
		$col = 8;
		next;
	}
	if ($ARGV[$i] eq '-fdr') {
		$fdrCutoff= $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-pvalue') {
		$pvalueCutoff= $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-log2fold') {
		$log2FoldCutoff= $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-peaks') {
		$col = 19;
		next;
	}
	if ($ARGV[$i] eq '-pc1') {
		$col = 19;
		$mode = "limma";
		$stitchFlag = 1;
		$normMethod = 'raw';
		next;
	}
	if ($ARGV[$i] eq '-tad') {
		$col = 10;
		$tadFlag = 1;
		$mode = "limma";
		$normMethod = 'raw';
		next;
	}
	if ($ARGV[$i] eq '-cp' || $ARGV[$i] eq '-loop') {
		$col = 10;
		next;
	}
	if ($ARGV[$i] eq '-basic') {
		$col = 1;
		next;
	}
	if ($ARGV[$i] eq '-norm2total') {
		$tagDirNormFlag = 1;
		next;
	}
	if ($ARGV[$i] eq '-dispersion') {
		$commonDispersion = $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] =~ /^\-/) {
		print STDERR "!!! $ARGV[$i] not recognized!\n";
		printCMD();
	}
	if ($batchFlag == 1) {
		if (!exists($batchs{$ARGV[$i]})) {
			$batchs{$ARGV[$i]} = $batchID++;
			push(@batchs, $ARGV[$i]);
		}
		push(@argvBatchs, $ARGV[$i]);
	} else {
		if (!exists($groups{$ARGV[$i]})) {
			push(@groups, $ARGV[$i]);
		}
		$groups{$ARGV[$i]} = 1;
		push(@argvGroups, $ARGV[$i]);
	}
}

# array of alpha names so that R doesn't through a fit or misparse things...
my @alpha = getAlphaNames();
my %inverseAlpha = ();
for (my $i=0;$i<@alpha;$i++) {
	$inverseAlpha{$alpha[$i]} = $i;
}

my $numGroups = scalar(@groups);
my $numBatchs = 0;
my $numComp = $numGroups*($numGroups-1);

if ($batchFlag) {
	if ($mode eq 'DESeq') {
		print STDERR "!!! Sorry - this script will only run -batch mode with edgeR or DESeq2\n";
		print STDERR "    You should probably be using DESeq2 instead anyway...\n";
		exit;
	}
}
if ($tagDirNormFlag) {
	if ($mode eq 'edgeR') {
		print STDERR "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		print STDERR "!!! Unfortunately -norm2total will only work with DESeq or DESeq2\n";
		print STDERR "    Analysis will continue, but be careful with the results!!!!\n";
		print STDERR "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		#exit;
	}
}



if ($col == 15) {
	print STDERR "\tTreating input as file generated by analyzeRNA.pl (-rna)\n";
} elsif ($col == 8) {
	print STDERR "\tTreating input as file generated by analyzeRepeats.pl (-repeats)\n";
} elsif ($col == 19) {
	print STDERR "\tTreating input as file generated by annotatePeaks.pl (-peaks)\n";
} elsif ($col == 10 && $tadFlag) {
	print STDERR "\tTreating input as file generated by findTADsAndLoops.pl (-tad)\n";
	$normMethod = 'raw' if ($normMethod eq '');
} elsif ($col == 10) {
	print STDERR "\tTreating input as file generated by findTADsAndLoops.pl (-loop)\n";
} elsif ($col == -1) {
	print STDERR "\tAutodetecting input file format...\n";
}


if ($batchFlag) {
	print STDERR "\tWill attempt to remove batch effects\n";
	if (scalar(@argvGroups) != scalar(@argvBatchs)) {
		my $ng = scalar(@argvGroups);
		my $nb = scalar(@argvBatchs);
		print STDERR "!!!! Entered different number of groups and batches on command line (groups=$ng,batch=$nb)!!!!\n";
		print STDERR "groups=@argvGroups\n";
		print STDERR "batch =@argvBatchs\n";
		exit;
	}
	$numBatchs = scalar(@groups);
	if ($numBatchs < 2) {
		print STDERR "!!!! Warning - you need at least 2 batch codes... otherwise omit \"-batch\" !!!\n";
		exit;
	}
}


my @colTotals = ();

my $defHeader = '';

open IN, $ARGV[0];
my $count = 0;
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $str = $line[0];

	if ($count == 1 && $col == -1) {
		if ($line[0] =~ /analyzeRepeats\.pl/) {
			$col = 8;
			print STDERR "\tAutodetected analyzeRepeats.pl file\n";
		} elsif ($line[0] =~ /analyzeRNA\.pl/) {
			$col = 15;
			print STDERR "\tAutodetected analyzeRNA.pl file\n";
		} elsif ($line[0] =~ /annotatePeaks\.pl/) {
			$col = 19;
			print STDERR "\tAutodetected annotatePeaks.pl file\n";
		} elsif ($line[0] =~ /^#TAD/) {
			$col = 10;
			$tadFlag = 1;
			print STDERR "\tAutodetected findTADsAndLoops.pl file (TADs)\n";
		} elsif ($line[0] =~ /^#Loop/ || $line[0] =~ /^#CP/) {
			$col = 10;
			$cpFlag = 1;
			print STDERR "\tAutodetected findTADsAndLoops.pl file (Loops)\n";
		} else {
			print STDERR "!!! Warning - could not autodetect the file format !!!\n";
			print STDERR "!!! Assuming it might be a simple/basic format with one column of geneIDs... !!\n";
			$col = 1;
		}
	}
	if ($count == 1) {
		$str .= " (cmd=$cmd)";
	}
			


	for (my $i=1;$i<@line;$i++) {
		$str .= "\t$line[$i]";
	}

	if ($count == 1) {
		print "$str";
		$defHeader = $str;
		for (my $i=$col;$i<@line;$i++) {
			$line[$i] =~ s/\"//g;
			$line[$i] =~ /^(.*?)\s/;
			push(@names, $1);
			push(@colTotals, 0);

			my $tagDirTotal=0;
			if ($line[$i] =~ /\(([\d\.]+) [tT]otal/) {
				$tagDirTotal = $1;
				#print STDERR "\t$line[$i] => $tagDirTotal\n";
			}
			push(@tagDirTotals, $tagDirTotal);

		}
		next;
	}	
	$og{$line[0]} = $str;

	my @d = ();
	for (my $i=$col;$i<@line;$i++) {
		push(@d,$line[$i]);
		$colTotals[$i-$col] += $line[$i];
	}
	$data{$line[0]} = \@d;



	my $name = $line[0];
	my $chr= $line[1];

	my $p = {c=>$line[1],s=>$line[2],e=>$line[3],d=>$line[4],n=>$line[0],diff=>0,fdr=>0,c1=>0,c2=>0};
	if (!exists($stitchData{$chr})) {
		my %a = ();
		$stitchData{$chr}=\%a;
	}
	$stitchData{$chr}->{$name} = $p;
	$chrAssignment{$name} = $chr;

}
close IN;

my $avgTotal = 0;
foreach(@colTotals) {
	$avgTotal+=$_;
}
if ($tagDirNormFlag == 1) {
	$avgTotal = 0;
	foreach(@tagDirTotals) {
		$avgTotal+=$_;
	}
}

if (@colTotals < 1) {
	print STDERR "!!! Error: Problem identifying columns totals - file format might not be correct!\n";
	exit;
}
$avgTotal /= scalar(@colTotals);



if ($tadFlag) {
	print STDERR "\tTAD format specified, using limma and skipping value normalization...\n";
	$mode = 'limma';
	$normMethod = 'raw' if ($normMethod eq '');
}
if ($cpFlag) {
	print STDERR "\tLoop format specified, using edgeR by performing simple normalization by default...\n";
	$mode = 'edgeR';
	$normMethod = 'simple' if ($normMethod eq '');
}

$normMethod = 'rlog' if ($normMethod eq '');

if ($normMethod eq 'simple') {
	foreach(keys %og) {
		my $id = $_;
		my $v = $og{$id};
		my @a = split /\t/,$v;
		for (my $i=$col;$i<@a;$i++) {
			if ($tagDirNormFlag == 1) {
				my $tdtotal = $tagDirTotals[$i-$col];
				my $name = $names[$i-$col];
				if ($tdtotal < 0.01) {
					print STDERR "!!! Error: Total reads for $name (column=$i) = $tdtotal\n";
					exit;
				}
				$a[$i] = sprintf("%.2lf",$a[$i]*$avgTotal/$tagDirTotals[$i-$col]);
			} else {
				my $colT = $colTotals[$i-$col];
				if ($colT < 0.01) {
					my $name = $names[$i-$col];
					print STDERR "!!! Error: column total for $name = $colT\n";
					exit;
				}
				$a[$i] = sprintf("%.2lf",$a[$i]*$avgTotal/$colT);
			}
		}
		$v = $a[0];
		for (my $i=1;$i<@a;$i++) {
			$v .= "\t$a[$i]";
		}
		$og{$id} = $v;
	}
} elsif ($normMethod eq 'rlog' || $normMethod eq 'vst') {

	my $inputFile = $rand . ".edgeR.in.data.txt";
	my $gFile = $rand . ".edgeR.groups.data.txt";
	my $normFile = $rand . ".edgeR.norm.data.txt";

	my $rOutputFile = $rand . ".edgeR.out.data.txt";
	my $rStderr = $rand . ".edgeR.out.stderr.txt";
	my $rScript = $rand . ".edgeR.R";

	open COUNTS, ">$inputFile";
	open GROUPS, ">$gFile";
	print GROUPS "Sample";
	if ($batchFlag) {
		print GROUPS "\tBatch";
	}
	print GROUPS "\tTreatment\n";
	open NORM, ">$normFile";
	print NORM "Sample\tNormFactor\tlibSize\n";

	my $z = 0;
	foreach(keys %og) {
		$z++;
		my $id = $_;
		my $v = $og{$id};
		my @a = split /\t/,$v;
		if ($z == 1) {
			#print STDERR "@a\n";
			print COUNTS "Gene";
			for (my $i=$col;$i<@a;$i++) {
				my $j = $i-$col;
				my $v = 'Sample' . $alpha[$j];
				print COUNTS "\t$v";

				print GROUPS "Sample" . $alpha[$j];
			   	if ($batchFlag) {
					my $b = $batchs{$argvBatchs[$j]};
					print GROUPS "\t$alpha[$b]";
				}
				#print STDERR " - $j - $argvGroups[$j] $a[$i]\n";
				if (defined($argvGroups[$j])) {
					print GROUPS "\t$argvGroups[$j]\n";
				} else {
					print GROUPS "\tSample" . "$alpha[$j]\n";
				}

				my $normF = 1;
				if ($tagDirTotals[$j] > 0 && $avgTotal > 0) {
					$normF = $tagDirTotals[$j]/$avgTotal;
				}
				print NORM "Sample" . $alpha[$j];
				print NORM "\t$normF\t$tagDirTotals[$j]\n";
			}
			print COUNTS "\n";
		}
		print COUNTS "$a[0]";
		for (my $i=$col;$i<@a;$i++) {
			my $v = floor($a[$i]+0.5);
			print COUNTS "\t$v";
		}
		print COUNTS "\n";
	}
	close COUNTS;
	close NORM;
	close GROUPS;

	open SCRIPT, ">$rScript";
	print SCRIPT "####### basic script for normalizing with DESeq2 (generated by HOMER) ########\n";
	print SCRIPT "library(DESeq2)\n";
	print SCRIPT "#Read Data in\n";
	print SCRIPT "countData <- read.delim(\"$inputFile\")\n";
	print SCRIPT "colData <- read.delim(\"$gFile\")\n\n";

	if ($batchFlag) {
		print SCRIPT "dds <- DESeqDataSetFromMatrix(countData, colData, design= ~Batch+Treatment, tidy=TRUE)\n";
	} else {
		print SCRIPT "dds <- DESeqDataSetFromMatrix(countData, colData,design=~Treatment,tidy=TRUE)\n";
	}
	if ($tagDirNormFlag == 1) {
		print SCRIPT "libData <- read.delim(\"$normFile\")\n\n";
		print SCRIPT "normFactors <- libData[,2]\n";
		print SCRIPT "libSizes <- libData[,3]\n";
		print SCRIPT "colData(dds)\$sizeFactor <- normFactors\n";
	}

	my $blind = 'FALSE';
	if ($normOnly) {
		$blind = 'TRUE';
	}
	if ($normMethod eq 'rlog') {
		print SCRIPT "norm <- rlog(dds,blind=$blind)\n";
	} elsif ($normMethod eq 'vst') {
		print SCRIPT "norm <- vst(dds,blind=$blind)\n";
	}

	print SCRIPT "norm_matrix <- assay(norm)\n";
	print SCRIPT "norm_df <- data.frame(Gene=rownames(norm_matrix), norm_matrix)\n";
	print SCRIPT "write.table(norm_df, \"$rOutputFile\", row.names = FALSE,sep=\"\\t\")\n";
	close SCRIPT;

	#print STDERR "`R --no-save < $rScript`;\n";
	#exit;
	print STDERR "\n\tPerforming variance stabalization ($normMethod)...\n\n";
	my $nsamples = scalar(@argvGroups); 
	if ($nsamples > 20 && $normMethod eq 'rlog') {
		print STDERR "\t!!! Warning, number of samples = $nsamples (>20). rlog transformation may take\n";
		print STDERR "\t    a long time.  Consider using -vst option instead if you are in a hurry.\n";     
	}
	if ($showRFlag) {
		print STDERR "R --no-save < $rScript\n";
		`R --no-save < "$rScript"`;
	} else {
		`R --no-save < "$rScript" 2> "$rStderr"`;
	}

	if (open IN, $rOutputFile) {
		$z = 0;
		while (<IN>) {	
			$z++;
			next if ($z < 2);
			chomp;
			s/\r//g;
			my @line = split /\t/;
			my $id  = $line[0];
			$id =~ s/\"//g;
	
			if (!exists($og{$id})) {
				print STDERR "! warning: what id = $id?\n";
				next;
			}
			my @a = split /\t/, $og{$id};
	
			$v = $a[0];
			for (my $i=1;$i<@a;$i++) {
				my $vv = $a[$i];
				if ($i >= $col) {
					my $j = $i-$col+1;
					$vv = $line[$j];
				}
				$v .= "\t$vv";
			}
			$og{$id} = $v;
		}
		close IN;
		`rm "$rStderr"`;
	} else {
		print STDERR "\tWarning: Likely a problem with variance stablizing transform\n";
		print STDERR "\tOften happens without replicates - consider specifying -simpleNorm\n";
		print STDERR "\tWill report raw counts in output table (no normalization)\n";
		print STDERR "\tSee $rStderr\n";
	}
	
	`rm "$inputFile" "$gFile" "$normFile" "$rOutputFile" "$rScript"`;

} elsif ($normMethod eq 'raw') {
	print STDERR "\tWill report raw counts in output table (no normalization)\n";
}
	
my $min = 1e10;
my $extra = '';
my @comps = ();

if ($normOnly == 0) {

if ($mode eq '') {
	$mode = 'DESeq2';
	print STDERR "\n\tWarning: getDiffExpression.pl now defaults to DESeq2 for calculations (instead of edgeR)\n\n";
}
print STDERR "\tUsing $mode to calculate differential expression/enrichment...\n";

for (my $i=0;$i<$numGroups;$i++) {
	last if ($i>0 && $AvsAflag == 0);
	for (my $j=$i+1;$j<$numGroups;$j++) {

		my $replicateFlag = 0;
		my $inputFile = $rand . ".edgeR.in.data.txt";
		my $gFile = $rand . ".edgeR.groups.data.txt";
		my $normFile = $rand . ".edgeR.norm.data.txt";

		my $rOutputFile = $rand . ".edgeR.out.data.txt";
		my $rStderr = $rand . ".edgeR.out.stderr.txt";
		my $rScript = $rand . ".edgeR.R";

		my %c1 = ();
		my %nc1 = ();
		my %c2 = ();
		my %nc2 = ();

		open COUNTS, ">$inputFile";
		open GROUPS, ">$gFile";
		open NORM, ">$normFile";
		my @dsamples = ();
		my %sampleCounts = ();
		my @totals = ();

		print COUNTS "Gene";
		print GROUPS "Sample";
		if ($batchFlag) {
			print GROUPS "\tBatch";
		}
		print GROUPS "\tTreatment\n";
		print NORM "Sample\tNormFactor\tLibSize\n";

		
		for (my $k=0;$k<@argvGroups;$k++) {
			# $k represents the sample index
			if ($argvGroups[$k] eq $groups[$i]) {
				print COUNTS "\tSample$alpha[$k]";
				push(@dsamples, "Sample" . $alpha[$k]);

				print GROUPS "Sample$alpha[$k]";
				if ($batchFlag) {
					my $b = $batchs{$argvBatchs[$k]};
					print GROUPS "\t$alpha[$b]";
				}
				$sampleCounts{$i}++;
				push(@totals, $tagDirTotals[$k]);
				print GROUPS "\t$alpha[$i]\n";
			} elsif ($argvGroups[$k] eq $groups[$j]) {
				print COUNTS "\tSample$alpha[$k]";
				push(@dsamples, "Sample" . $alpha[$k]);

				print GROUPS "Sample$alpha[$k]";
				if ($batchFlag) {
					my $b = $batchs{$argvBatchs[$k]};
					print GROUPS "\t$alpha[$b]";
				}
				$sampleCounts{$j}++;
				push(@totals, $tagDirTotals[$k]);
				print GROUPS "\t$alpha[$j]\n";
			}
		}
		print COUNTS "\n";
		close GROUPS;

		foreach(values %sampleCounts) {
			if ($_ > 1) {
				$replicateFlag = 1;
			}
		}

		if ($tagDirNormFlag == 1) {
			my $T = 0;
			$T+=$_ foreach(@totals);
			my $avg = $T/@totals if (@totals > 0);
			my $normLine="";
			my $sizeLine="";
			my $z = 0;
			foreach(@totals) {
				my $size = $_;
				my $normF = 1;
				if ($size > 0) {
					$normF = $size/$avg;
				}
				print NORM "$dsamples[$z]\t$normF\t$size\n";
			}
		}
		close NORM;

		foreach(keys %data) {
			my $gene = $_;
			print COUNTS "$gene";
			for (my $k=0;$k<@argvGroups;$k++) {
				if ($argvGroups[$k] eq $groups[$i] || $argvGroups[$k] eq $groups[$j]) {
					my $v = $data{$gene}->[$k];
					if ($mode eq 'limma') {
					} else {
						$v = floor($v+0.5);
					}
					print COUNTS "\t$v";
				}
				if ($argvGroups[$k] eq $groups[$i]) {
					$c1{$gene} +=  $data{$gene}->[$k];
					$nc1{$gene}++;
				}
				if ($argvGroups[$k] eq $groups[$j]) {
					$c2{$gene} +=  $data{$gene}->[$k];
					$nc2{$gene}++;
				}
			}
			print COUNTS "\n";
		}
		close COUNTS;

		open SCRIPT, ">$rScript";
		if ($mode eq 'limma') {
			print SCRIPT "####### basic script for running limma for PC1/TADs (generated by HOMER) ########\n";
			print SCRIPT "#### input is matrix data - first column is TAD id, others are TAD inclusion scores\n";
			print SCRIPT "#### <blank> Exp1 Exp2\n#### chr1:300000-400000  2.45  2.01\n";
			print SCRIPT "library(limma)\n";
			print SCRIPT "#Read Data in\n";
			print SCRIPT "countData <- read.delim(\"$inputFile\",row.names=1)\n";
			print SCRIPT "colData <- read.delim(\"$gFile\")\n\n";
			if ($batchFlag) {
				print SCRIPT "design <-model.matrix(~Batch+Treatment, data=colData)\n";
			} else {
				print SCRIPT "design <-model.matrix(~Treatment, data=colData)\n";
			}
			print SCRIPT "fit <- lmFit(countData, design)\n";
			print SCRIPT "fit2 <- eBayes(fit)\n";
			print SCRIPT "res <- topTable(fit2,n=Inf)\n";
			print SCRIPT "tidyResult <- data.frame(TAD=rownames(res), res)\n";
			print SCRIPT "write.table(tidyResult,file=\"$rOutputFile\",sep=\"\\t\",row.names=FALSE)\n\n";
		} elsif ($mode eq 'edgeR') {
			print SCRIPT "####### basic script for running edgeR (generated by HOMER) ########\n";
			print SCRIPT "#### input is matrix data - first column is gene id, others are tag un-normalized tag counts\n";
			print SCRIPT "#### <blank> Exp1 Exp2\n#### NM_123   1233  444\n##### NM_234    424  198\n";
			print SCRIPT "library(edgeR)\n";
			print SCRIPT "#Read Data in\n";
			print SCRIPT "countData <- read.delim(\"$inputFile\",row.names=1)\n";
			print SCRIPT "colData <- read.delim(\"$gFile\")\n\n";


			if ($tagDirNormFlag == 1) {
				print SCRIPT "libData <- read.delim(\"$normFile\")\n\n";
				print SCRIPT "normFactors <- libData[,2]\n";
				print SCRIPT "libSizes <- libData[,3]\n";
				print SCRIPT "y <- DGEList(countData, samples=colData, group=colData\$Treatment, lib.size=libSizes,norm.factors=normFactors)\n\n";
			} else {
				print SCRIPT "y <- DGEList(countData, samples=colData, group=colData\$Treatment)\n";
				print SCRIPT "y <- calcNormFactors(y)\n\n";
			}
	
			if ($batchFlag) {
				print SCRIPT "design <-model.matrix(~Batch+Treatment, data=colData)\n";
				print SCRIPT "y <- estimateDisp(y, design)\n";
				print SCRIPT "fit <- glmQLFit(y,design)\n";
				print SCRIPT "qlf <- glmQLFTest(fit,ncol(design))\n";
				print SCRIPT "res <- topTags(qlf,Inf)\n";
			} else {
				print SCRIPT "#Estimate Error Model\n";
				if ($replicateFlag) {
					print SCRIPT "design <-model.matrix(~Treatment, data=colData)\n";
					print SCRIPT "y <- estimateDisp(y, design)\n";
				} else {
					print STDERR "\tSetting commonDispersion to $commonDispersion\n";
					print SCRIPT "y\$common.dispersion <- $commonDispersion\n";
				}
				print SCRIPT "\n";
				print SCRIPT "#compute p-values, then output\n";
				print SCRIPT "et <- exactTest(y)\n";
				print SCRIPT "res <- topTags(et,Inf)\n";
			}
			print SCRIPT "tidyResult <- data.frame(Gene=rownames(res\$table), res\$table)\n";
			print SCRIPT "write.table(tidyResult,file=\"$rOutputFile\",sep=\"\\t\",row.names=FALSE)\n\n";
			print SCRIPT "write(\"Dispersion = \",stderr())\n";
			print SCRIPT "write(y\$common.dispersion,stderr())\n";

		} elsif ($mode eq 'DESeq') {
			print SCRIPT "####### basic script for running DESeq (generated by HOMER) ########\n";
			print SCRIPT "library(DESeq)\n";
			print SCRIPT "#Read Data in\n";
			print SCRIPT "countData <- read.delim(\"$inputFile\",row.names=1)\n";
			print SCRIPT "colData <- read.delim(\"$gFile\")\n\n";

			if ($batchFlag) {
					#print SCRIPT "bin <- as.matrix(read.table(\"$bFile\"))\n\n";
				#print SCRIPT "dds <- DESeqDataSetFromMatrix(countData = data, colData = gin,design=~batch+group)\n";
				print SCRIPT "colData <- read.delim(\"$gFile\")\n\n";
				print SCRIPT "cds = newCountDataSet(countData,colData)\n";
			} else {
				print SCRIPT "cds = newCountDataSet(countData,colData[,2])\n";
			}
			if ($tagDirNormFlag == 1) {
				print SCRIPT "libData <- read.delim(\"$normFile\")\n\n";
				print SCRIPT "normFactors <- libData[,2]\n";
				print SCRIPT "libSizes <- libData[,3]\n";
				print SCRIPT "pData(cds)\$sizeFactor <- normFactors\n";
			} else {
				print SCRIPT "cds = estimateSizeFactors(cds)\n";
			}

			if ($batchFlag) {
				print SCRIPT "cds = estimateDispersions(cds,\"pooled-CR\",modelFormula=count ~Batch+Treatment)\n";
				print SCRIPT "fit1 = fitNbinomGLMs(cds, count ~Batch + Treatment)\n";
				print SCRIPT "fit0 = fitNbinomGLMs(cds, count ~Batch)\n";
				print SCRIPT "pvalsGLM = nbinomGLMTest( fit1, fit0 )\n";
				print SCRIPT "padjGLM = p.adjust(pvalsGLM, method=\"BH\")\n";
				print SCRIPT "write.table(padjGLM,file=\"$rOutputFile\",sep=\"\\t\",row.names=FALSE)\n\n";
			} else {
				if ($replicateFlag) {
					print SCRIPT "cds = estimateDispersions(cds)\n";
				} else {
					print SCRIPT "cds = estimateDispersions(cds,method=\"blind\", sharingMode=\"fit-only\")\n";
				}
				print SCRIPT "res = nbinomTest(cds,\"A\",\"B\")\n";
				print SCRIPT "write.table(res,file=\"$rOutputFile\",sep=\"\\t\",row.names=FALSE)\n\n";
			}

		} elsif ($mode eq 'DESeq2') {
			print SCRIPT "####### basic script for running DESeq (generated by HOMER) ########\n";
			print SCRIPT "library(DESeq2)\n";
			print SCRIPT "#Read Data in\n";
			print SCRIPT "countData <- read.delim(\"$inputFile\")\n";
			print SCRIPT "colData <- read.delim(\"$gFile\")\n\n";

			if ($batchFlag) {
				print SCRIPT "dds <- DESeqDataSetFromMatrix(countData, colData, design= ~Batch+Treatment, tidy=TRUE)\n";
			} else {
				print SCRIPT "dds <- DESeqDataSetFromMatrix(countData, colData,design=~Treatment,tidy=TRUE)\n";
			}

			if ($tagDirNormFlag == 1) {
				print SCRIPT "libData <- read.delim(\"$normFile\")\n\n";
				print SCRIPT "normFactors <- libData[,2]\n";
				print SCRIPT "libSizes <- libData[,3]\n";

				print SCRIPT "colData(dds)\$sizeFactor <- normFactors\n";
				print SCRIPT "dds <- estimateDispersions(dds)\n";
				print SCRIPT "dds <- nbinomWaldTest(dds)\n";

			} else {
				print SCRIPT "dds <- DESeq(dds)\n";
			}
			print SCRIPT "res <- results(dds,tidy=TRUE)\n";
			print SCRIPT "write.table(res,file=\"$rOutputFile\",sep=\"\\t\",row.names=FALSE)\n\n";
		}
		close SCRIPT;

		#print STDERR "`R --no-save < $rScript`;\n";
		#exit;
		if ($showRFlag) {
			`R --no-save < "$rScript"`;
		} else {
			`R --no-save < "$rScript" 2> "$rStderr"`;
		}

		my $numUp = 0;
		my $numDown = 0;
		my $numTotal = 0;

		if ($export ne '') {
			open UP, ">$export.Up_$groups[$j]" . "_vs_" . $groups[$i] . ".txt";
			open DOWN, ">$export.Down_$groups[$j]" . "_vs_" . $groups[$i] . ".txt";
		}

		my $zzTotalLines=0;
		my $zzTotalWrongName=0;

		my %comp = ();
		open IN, $rOutputFile;
		$count = 0;
		while (<IN>) {
			$count++;
			chomp;
			s/\r//g;
			s/\"//g;
			my @line = split /\t/;
			foreach(@line) {
				s/^\s+//;
				s/\s+$//;
			}


			if ($count == 1) {
				if ($mode eq 'limma') {
					print "\t$groups[$i] vs. $groups[$j] Difference";
				} else {
					print "\t$groups[$i] vs. $groups[$j] Log2 Fold Change";
				}
				print "\t$groups[$i] vs. $groups[$j] p-value";
				print "\t$groups[$i] vs. $groups[$j] adj. p-value";
				if ($export ne '') {
					print UP "$defHeader";
					if ($mode eq 'limma') {
						print UP "\t$groups[$i] vs. $groups[$j] Difference";
					} else {
						print UP "\t$groups[$i] vs. $groups[$j] Log2 Fold Change";
					}
					print UP "\t$groups[$i] vs. $groups[$j] p-value";
					print UP "\t$groups[$i] vs. $groups[$j] adj. p-value";
					print UP "\n";
					print DOWN "$defHeader";
					if ($mode eq 'limma') {
						print DOWN "\t$groups[$i] vs. $groups[$j] Difference";
					} else {
						print DOWN "\t$groups[$i] vs. $groups[$j] Log2 Fold Change";
					}
					print DOWN "\t$groups[$i] vs. $groups[$j] p-value";
					print DOWN "\t$groups[$i] vs. $groups[$j] adj. p-value";
					print DOWN "\n";
				}
				next;
			}

			$zzTotalLines++;

			my $id = '';
			my $f = 0.0;
			my $p = 1.0;
			my $fdr = 1.0;
			my $str = "";

			if ($mode eq 'limma') {
				$id = $line[0];
				$f = $line[1];
				$p = $line[4];
				$fdr = $line[5];
				my $avg = $line[2];
				if ($batchFlag) {
					$f = $line[2];
					$avg = $line[3];
					$p = $line[5];
					$fdr = $line[6];
				}
				#$c1 = $avg - ($f/2);
				#$c2 = $avg + ($f/2);
			} elsif ($mode eq 'edgeR') {
				$id = $line[0];
				$f = $line[1];
				$p = $line[3];
				$fdr = $line[4];
				if ($batchFlag) {
					$p = $line[4];
					$fdr = $line[5];
				}
			} elsif ($mode eq 'DESeq') {
				$id = $line[0];
				$f = $line[5];
				$p = $line[6];
				$fdr = $line[7];
			} elsif ($mode eq 'DESeq2') {
				$id = $line[0];
				$f = $line[2];
				$p = $line[5];
				$fdr = $line[6];
			}

			if ($f eq '-Inf' || $f eq '-inf') {
				$f = -10;
			}
			if ($f eq 'Inf' || $f eq 'inf') {
				$f = 10;
			}
			if ($fdr eq 'NA' || $fdr > 1) {
				$fdr = 1;
			}
			if ($p eq 'NA' || $p > 1) {
				$p = 1;
			}

			$id =~ s/^\"//;
			$id =~ s/\"$//;
			if (!exists($og{$id})) {
				$zzTotalWrongName++;
				#print STDERR "id=$id\n";
				next;
			}

			$str = "\t$f\t$p\t$fdr";

			$numTotal++;
			if (($pvalueCutoff ne '' && $p < $pvalueCutoff) || ($fdr < $fdrCutoff)) {
				if ($f > $log2FoldCutoff) {
					$numUp++;
					print UP "$og{$id}" . "$str\n" if ($export ne '');
				} elsif ($f < -1*$log2FoldCutoff) {
					$numDown++;
					print DOWN "$og{$id}" . "$str\n" if ($export ne '');
				}
			}

			#$min = $line[$sIndex+1] if ($line[$sIndex+1] ne 'NA' && $line[$sIndex+1] < $min);
			$comp{$id}=$str;


			my $chr = '';
			next if (!exists($chrAssignment{$id}));
			$chr = $chrAssignment{$id};
			next if (!exists($stitchData{$chr}->{$id}));
			next if (!exists($c1{$id}));
			next if (!exists($c2{$id}));
				
			my $diff = $f;

			$stitchData{$chr}->{$id}->{'diff'} = $diff;
			$stitchData{$chr}->{$id}->{'fdr'} = $fdr;
			$stitchData{$chr}->{$id}->{'c1'} = $c1{$id}/$nc1{$id};
			$stitchData{$chr}->{$id}->{'c2'} = $c2{$id}/$nc2{$id};
		}
		close IN;

		if ($export ne '') {
			close UP;
			close DOWN;
		}

		push(@comps, \%comp);

		#print STDERR "totalLines=$zzTotalLines\ntotalWrongName=$zzTotalWrongName\n";

		if ($numTotal == 0) {
			print STDERR "!!! Warnining - something likely failed during R execution.\n";
			print STDERR "!!! R script that was used:\n";
			open IN1, "$rScript";	
			while (<IN1>) {
				print STDERR "\t$_";
			}
			close IN1;
			if ($showRFlag==0) {
				print STDERR "!!! R execution details:\n";
				open IN1, "$rStderr";	
				while (<IN1>) {
					print STDERR "\t$_";
				}
				close IN1;
			}
			exit;
		}
		if ($showRFlag==1) {
			open IN1, "$rScript";	
			while (<IN1>) {
				print STDERR "\t$_";
			}
			close IN1;

		}

		`rm "$inputFile" "$gFile" "$rOutputFile" "$rScript" "$normFile"`;
		if ($showRFlag==0) {
			`rm "$rStderr"`;
		}

		if ($numTotal == 0) {
			exit;
		}
		my $sigStr = "FDR<$fdrCutoff";
		if ($pvalueCutoff ne '') {
			$sigStr = "p<$pvalueCutoff";
		}

		print STDERR "\tOutput Stats $groups[$i] vs. $groups[$j]:\n";
		print STDERR "\t\tTotal Genes: $numTotal\n";
		my $r = sprintf("%.3f",$numUp/$numTotal*100);
		print STDERR "\t\tTotal Up-regulated in $groups[$j] vs. $groups[$i]: $numUp ($r%) [log2fold>$log2FoldCutoff, $sigStr]\n";
		$r = sprintf("%.3f",$numDown/$numTotal*100);
		print STDERR "\t\tTotal Dn-regulated in $groups[$j] vs. $groups[$i]: $numDown ($r%) [log2fold<-$log2FoldCutoff, $sigStr]\n";

		if ($stitchFlag==1 && $export ne '') {
			my $outputFile = "$export.Up_$groups[$j]" . "_vs_" . $groups[$i] . ".regions.txt";
			my $numUp = stitchRegions(\%stitchData, $outputFile, $fdrCutoff, $log2FoldCutoff, $stitchMinLen, 
															$stitchMinPos, $stitchMaxNeg, $stitchBridgeLen,0);
			$outputFile = "$export.Down_$groups[$j]" . "_vs_" . $groups[$i] . ".regions.txt";
			my $numDown = stitchRegions(\%stitchData, $outputFile, $fdrCutoff, $log2FoldCutoff, $stitchMinLen,
															$stitchMinPos, $stitchMaxNeg, $stitchBridgeLen,1);

			print STDERR "\n\t\tTotal continuous domains: $numUp (Up) and $numDown (Down)\n";
		}
	}
}

} # normOnly

print "\n";

foreach(keys %og) {
	print "$og{$_}";
	my $gene = $_;
	foreach(@comps) {
		if (exists($_->{$gene})) {
			print $_->{$gene};
		} else {
			print "\t0\t1\t1";
		}
	}
	print "\n";
}
print STDERR "\n";


sub getAlphaNames {
	my @a = ();
	for (my $i=65;$i<91;$i++) {
		push(@a, chr($i));
	}
	for (my $i=65;$i<91;$i++) {
		for (my $j=65;$j<91;$j++) {
			push(@a, chr($j) . chr($i));
		}
	}
	return @a;
}





sub stitchRegions {
	my ($stitchData, $outputFile, $fdrThresh, $diffThresh, $minLen, $minPos, $maxNeg, $bridgeLen, $reverseFlag) = @_;

	open OUT, ">$outputFile";
	print OUT "#PC1 region\tchr\tstart\tend\tstrand\tPC1-cond1\tPC1-cond2\tAvg. Difference\tAvg. FDR\tLength\n";

	my $numRegions = 0;
	foreach(keys %$stitchData) {
		my $chr = $_;
		my @p = sort {$a->{'s'} <=> $b->{'s'}} values %{$stitchData->{$chr}};
	
		
		my $good = 0;
		my $curStart = -1;
		my $curEnd = -1;
		my $avgDiff = 0;
		my $avgFDR = 0;
		my $avgC1 = 0;
		my $avgC2 = 0;
		my $N = 0;
		my $f = 1;
		if ($reverseFlag) {
			$f = -1;
		}
	
		for (my $i=0;$i<@p;$i++) {
			if ($p[$i]->{'fdr'} < $fdrThresh && $f*$p[$i]->{'diff'} > $diffThresh &&
					( ($f== 1 && $p[$i]->{'c1'} < $maxNeg && $p[$i]->{'c2'} > $minPos)
					 || ($f== -1 && $p[$i]->{'c2'} < $maxNeg && $p[$i]->{'c1'} > $minPos))) {
				if ($good) {
				} else {
					$good = 1;
					$curStart = $p[$i]->{'s'};
					$avgDiff = 0;
					$avgFDR = 0;
					$avgC1 = 0;
					$avgC2 = 0;
					$N = 0;
				}
				$curEnd = $p[$i]->{'e'};
				$avgDiff += $f*$p[$i]->{'diff'};
				$avgFDR += $p[$i]->{'fdr'};
				$avgC1 += $p[$i]->{'c1'};
				$avgC2 += $p[$i]->{'c2'};
				$N++;
			} else {
				if ($good && $p[$i]->{'s'}-$curEnd > $bridgeLen) {
	
					$avgDiff/=$N;
					$avgFDR/=$N;
					$avgC1/=$N;
					$avgC2/=$N;
					my $len = $curEnd-$curStart;
					$good = 0 if ($len < $minLen);
					if ($good) {
						print OUT "PC1-$chr:$curStart-$curEnd\t$chr\t$curStart\t$curEnd\t+\t$avgC1\t$avgC2\t$avgDiff\t$avgFDR\t$len\n";
						$numRegions++;
					}
					$good = 0;
				}
			}
	
		}
	}
	close OUT;
	return $numRegions;
	#print STDERR "\ttotal regions = $numRegions\n";
	#exit;
}

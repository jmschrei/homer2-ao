#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";


# Copyright 2009 - 2024 Christopher Benner <cbenner@ucsd.edu>
# #
# # This file is part of HOMER
# #
# # HOMER is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.
# #
# # HOMER is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.
#
use POSIX;
use HomerConfig;
use Statistics;
use HomerSeqBias;
use HomerSVGLogo;

my $maxCPUs = 56;
my $maxCPUs2 = 30;
my $expN = 10000;
my $expLen = 100;
my $kmer = 2;
my $logratioFile = '';
my $numRandSeq = -1;
my $pseudo = 0.00005;
my $l2ratioEdge = 20;
my $tgSeqFile = '';
my $bgSeqFile = '';
my $bgPosFile = '';
my $strand = 'both';
my $prefix = '';
my $outputDirectory = "";
my $mode = '';
my $denovoLength = 10;
my $homerOpts = "";
my $bitFlag =1;
my $denovoS = 5;
my $seqPropertyOrder = 1;
my $maxN = 0.50;
my $stratifiedModelFlag = "";
my @modelBinValues = ();
my $numBins = 10;
my $nsubBins = 10;
my $positionIndependentModelFlag=0;
my $offset = "default";
my $NN = 100000000;
my $allowTargetOverlap = 0;
my $allowBgOverlap = 0;
my $modelFlag = 0;
my @summaryRes = (3, 25, 100);
my $fourierFlag = 1;
my $fourierStart = -120;
my $fourierEnd = -40;
my $calcSeqPropertiesFlag = 0;
my $forceSummaryByWindow = 0;
my $enrichedFlag = 1;
my $depletedFlag = 1;

my $reduceThresh = 0.6;
my $matchThresh = "T10";	
my $rnaFlag = "";
my $bitsFlag = " -bits ";
my $knownMotifsFile = "/gpfs/data01/bennerlab/home/cbenner/software/homer/data/knownTFs/vertebrates/all.motifs";

my $seqBackgroundEXE = "conda run -n seqbackground python generate_background.py";

sub printCMD {
	print STDERR "\n\tusage:\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-o <output directory> (required)\n";
	print STDERR "\t\t-strand <+|-|separate|both> (default, both)\n";
	print STDERR "\t\t-cpu <#> (number of CPUs to use, def: $maxCPUs)\n";
	print STDERR "\t\t-cpu2 <#> (number of CPUs to per HOMER denovo analysis, def: $maxCPUs2)\n";
	print STDERR "\tProfiling enrichment of known/given motifs:\n";
	print STDERR "\t\t-m <motif file>\n";
	print STDERR "\tFinding motifs de novo\n";
	print STDERR "\t\t-denovo (look for motifs)\n";
	print STDERR "\t\t-len <#> (length of motif to find, def: $denovoLength)\n";
	print STDERR "\t\t-S <#> (number of de novo motifs to find, def: $denovoS)\n";
	print STDERR "\tEnrichment regions: (default: -windows 3)\n";
	print STDERR "\t\t-all (check all positions at all resolutions)\n";
	print STDERR "\t\t-regions <#,#> [#,#] ... (regions to search for enriched motifs)\n";
	print STDERR "\t\t-windows <#>[,#,#,...] | -wstart <#>, -wend <#>, -wstep <#> (size of enrichment region)\n";
	print STDERR "\t\t-positions <#>[,#,#,...] | -pstart <#>, -pend <#>, -pstep <#> (location of enrichment region)\n";
	print STDERR "\tTarget Sequences:\n";
	print STDERR "\t\t-p <peak/BED file> (i.e. tss positions)\n";
	print STDERR "\t\t-g <genome FASTA file> (for now you must provide a fasta file, not a 'homer genome'\n";
	print STDERR "\t\t-size <#> (resize regions to this size, -size #,# ok to use)\n";
	print STDERR "\t\t-mask (use repeat masked sequences)\n";
	print STDERR "\t\t-maxN <#> (Maximum fraction of Ns in sequence to still consider in analysis, default=$maxN)\n";
	print STDERR "\t\t\t-or-\n";
	print STDERR "\t\t-tg <fasta file> (target sequences in FASTA format)\n";
	print STDERR "\t\t-offset <#> (offset of fasta sequences for reporting, default: $offset)\n";
	print STDERR "\tControl Sequences:\n";
	print STDERR "\t\t-pkmer <#> (length of kmers to model sequence properties with by position, default: $kmer)\n";
	print STDERR "\t\t-ikmer <#> (length of kmers to model sequence properties with independent of postiion intsead)\n";
	print STDERR "\t\t-model (intead of using actual sequences, generate random sequences matching the desired kmer content)\n";
	print STDERR "\t\t-nbins <#> (number of GC bins to stratify sequences by, def: $numBins)\n";
	print STDERR "\t\t-nsubBins <#> (number of sub bins used to help stratify kmer content, def: $nsubBins)\n";
	print STDERR "\t\t-allowTargetOverlap (allow background regions to overlap target regions, by default not allowed)\n";
	print STDERR "\t\t-allowBgOverlap (allow background regions to overlap with other background regions, by default not allowed)\n";
	print STDERR "\t\t-N <#> (number of random sequences, def: same as target sequence set)\n";
	print STDERR "\t\t-NN <#> (number of initial sequences to extract from genome for selection, def: $NN)\n";
	print STDERR "\t\t\t-or- (these options will be used 'as is' and not normalized)\n";
	print STDERR "\t\t-bg <fasta file> (background sequences in FASTA format)\n";
	print STDERR "\t\t-bgp <peak/BED file> (background positions)\n";
	#print STDERR "\tregions:\n";
	print STDERR "\n\tOther options:\n";
	print STDERR "\t\t-fourierRange <#>,<#> (region to calculate frequency spectrum on, defaut: $fourierStart,$fourierEnd)\n";
	print STDERR "\t\t-calcSeqProperties (calculate sequence nucletodie pattern rates)\n";
	print STDERR "\n\tNotes:\n";
	print STDERR "\t\t-enrichment values are returned as natural log p values, where '+' values indicate depletion\n";
	print STDERR "\n";
	exit;
}

my $wfixed = '';
my @windows = (3);
my $wstart = '';
my $wend = '';
my $wstep = '';
my $pfixed = '';
my @positions = ();
my $pstart = '';
my $pend = '';
my $pstep = '';
my @regions = ();
my $allFlag = 0;

my $genomePath = '';
my $targetPosFile = '';
my $motifFile = '';
my $mask = '';
my $size = '';
my $empiricalFlag = 0;
my $cmd = '';
my $ogSize = '';
for (my $i=0;$i<@ARGV;$i++) {
	$cmd .= " " . $ARGV[$i];
}

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-p') {
		$targetPosFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-o') {
		$outputDirectory = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-m') {
		$motifFile = $ARGV[++$i];
		$mode = 'known';
	} elsif ($ARGV[$i] eq '-denovo') {
		$mode = 'denovo';
	} elsif ($ARGV[$i] eq '-len') {
		$denovoLength = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxN') {
		$maxN = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-stratifyByGC') {
		$stratifiedModelFlag = "gc";
	} elsif ($ARGV[$i] eq '-S') {
		$denovoS = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rna') {
		$rnaFlag = '-rna';
		$strand = "+";
	} elsif ($ARGV[$i] eq '-strand') {
		$strand = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-calcSeqProperties') {
		$calcSeqPropertiesFlag=1;
	} elsif ($ARGV[$i] eq '-fourierRange') {
		my @a = split /\,/,$ARGV[++$i];
		$fourierStart = $a[0];
		$fourierEnd = $a[1];
	} elsif ($ARGV[$i] eq '-pkmer') {
		$kmer = $ARGV[++$i];
		$positionIndependentModelFlag = 0;
	} elsif ($ARGV[$i] eq '-ikmer') {
		$kmer = $ARGV[++$i];
		$positionIndependentModelFlag = 1;
	} elsif ($ARGV[$i] eq '-model' || $ARGV[$i] eq '-modelBg') {
		$modelFlag = 1;
	} elsif ($ARGV[$i] eq '-bg') {
		$bgSeqFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bgp') {
		$bgPosFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-tg') {
		$tgSeqFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-offset') {
		$offset = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-empirical') {
		$empiricalFlag = 1;
	} elsif ($ARGV[$i] eq '-size') {
		$size = $ARGV[++$i];
		$ogSize = $size;
		# not sure if we need this
		if ($size eq 'given') {
			print STDERR "\tMust specify a size for this program (e.g. -size 200)\n";
			next;
		} elsif ($size =~ /\,/) {
			my @a = split /\,/, $size;
			my $sizeStart= $a[0];
			my $sizeEnd = $a[1];
			if ($sizeEnd < $sizeStart) {
				print STDERR "!!! Size end must be less than the size start range in -size $sizeStart,$sizeEnd\n";
				exit;
			}
			$size = ($sizeEnd - $sizeStart);
			$offset = $sizeStart;
		} else {
			$offset = -1*floor($size/2);
		}
		print STDERR "\tPeak Region set to $size (offset = $offset)\n";
	} elsif ($ARGV[$i] eq '-regions') {
		my @r = ();
		$i++;
		while ($i < @ARGV && $ARGV[$i] =~ /\,/) {
			my @a = split /\,/, $ARGV[$i];
			my $x = {s=>$a[0],e=>$a[1]};
			push(@regions, $x);
			$i++;
		}
		$i-- if ($i < @ARGV);
	} elsif ($ARGV[$i] eq '-nbins') {
		$numBins = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-nsubBins') {
		$nsubBins = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bins') {
		@modelBinValues = ();
		$i++;
		while ($i < @ARGV && $ARGV[$i] =~ /\-/) {
			push(@modelBinValues,$ARGV[$i]);
			$i++;
		}
		$i-- if ($i < @ARGV);
	} elsif ($ARGV[$i] eq '-N') {
		$numRandSeq = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-NN') {
		$NN = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-logratio') {
		$logratioFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-g') {
		$genomePath = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu2') {
		$maxCPUs2 = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pseudoRate') {
		$pseudo = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pstart') {
		$pstart = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-allowTargetOverlap') {
		$allowTargetOverlap=1;
	} elsif ($ARGV[$i] eq '-allowBgOverlap') {
		$allowBgOverlap=1;
	} elsif ($ARGV[$i] eq '-pend') {
		$pend = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pstep') {
		$pstep = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pfixed') {
		$pfixed = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-wstart') {
		$wstart = $ARGV[++$i];
		@windows=();
	} elsif ($ARGV[$i] eq '-wend') {
		$wend = $ARGV[++$i];
		@windows=();
	} elsif ($ARGV[$i] eq '-wstep') {
		$wstep = $ARGV[++$i];
		@windows=();
	} elsif ($ARGV[$i] eq '-wfixed') {
		$wfixed = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-windows') {
		@windows = split /\,/,$ARGV[++$i];;
	} elsif ($ARGV[$i] eq '-positions') {
		@positions = split /\,/,$ARGV[++$i];
	} elsif ($ARGV[$i] eq '-all') {
		$allFlag = 1;
	} elsif ($ARGV[$i] eq '-mask') {
		$mask = " -mask ";
	} else {
		print STDERR "What is \"$ARGV[$i]\"?\n";
		printCMD();
	}
}

if ($outputDirectory eq '') {
	print STDERR "!!! Error: Output directory (-o) required!!!\n";
	printCMD();
	exit;
}
if ($targetPosFile eq '' && $tgSeqFile eq '') {
	print STDERR "!!! -p <peak/BED file> or -tg <FASTA file> is required!\n";
	printCMD();
	exit;
}
if ($mode eq '') {
	print STDERR "!!! Either \"-denovo\" or \"-m <motifFile>\" are required!!!\n";
	printCMD();
}



`mkdir -p "$outputDirectory"`;
open OUT, ">$outputDirectory/cmd.txt" or die "!!! Error: Could not open file \"$$outputDirectory/cmd.txt\" for writing\n";
print OUT "cmd: $cmd\n";
close OUT;

my %mlen = ();


my $rand = rand();
my $tmpFile = $outputDirectory . "/" . $rand . ".tmp";
my $tmpFile2 = $outputDirectory . "/" . $rand . ".2.tmp";

my $targetSeqFile = $outputDirectory . "/targetSequences.tsv";


my %zzz = ();
my $seqs = \%zzz;
my $badLen = 0;
my $numSeqs =0;

my $targetPosResized = $outputDirectory . "/target.positions.resized.txt";
my $bgPosResized = $outputDirectory . "/target.positions.resized.txt";

if ($tgSeqFile eq '') {

	if ($ogSize ne '') {
		print STDERR "\tRegion size set to $ogSize\n";
		`adjustPeakFile.pl "$targetPosFile" -size $ogSize > "$targetPosResized"`;
		`homerTools extract "$targetPosResized" "$genomePath" $mask > "$tmpFile2"`;
	} else {
		print STDERR "\tRegion size set by file\n";
		`cp "$targetPosFile" "$targetPosResized"`;
		`homerTools extract "$targetPosResized" "$genomePath" $mask > "$tmpFile2"`;
	}
	
	my $foundSize = '';	
	open IN, $tmpFile2;
	while (<IN>){
		chomp;
		my @line = split /\t/;
		my $id = $line[0];
		my $seq = $line[1];
		my $slen = length($seq);
		if ($size ne '') {
			if ($slen != ($size+1)) {
				print STDERR "\t\t$id has the wrong sequence length (slen=$slen,size=$size+1)\n";
				next;
			}
		} else {
			$foundSize = $slen;
		}
		$seqs->{$id} = $seq;
		$numSeqs++;
	}
	close IN;
	`rm "$tmpFile2"`;
	if ($size eq '') {
		$size = $foundSize;
		print STDERR "\tSize of regions set to $size (based on target regions)\n";
	}
} else {
	unless (-e $tgSeqFile) {
		print STDERR "!!! Error: Could not open target FASTA file: $tgSeqFile\n";
		exit;
	}
	my $lenFlag = 1;
	$seqs = HomerSeqBias::readFASTA($tgSeqFile);
	my $len = -1;
	foreach(keys %$seqs) {
		my $L = length($seqs->{$_});
		$len = $L if ($len < 0);
		if ($L ne $len) {
			print STDERR "!!! Warning - target sequences are not the same length ($L != $len), seqID $_\n";
			next;
		}
		$numSeqs++;
	}
	print STDERR "\tInput target sequences read from file $tgSeqFile\n";
	print STDERR "\t\tlength/size = $len (offset = $offset)\n";
	print STDERR "\t\tN = $numSeqs\n";
	$size = $len;
	print STDERR "\t\tSetting size to length ($size)\n";
}	
if ($offset eq 'default') {
	$offset = -1*floor($size/2);
	print  STDERR "\tSetting default offset at half the sequence size: $offset\n";
}

if ($maxN < 1.0) {
	print STDERR "\tFiltering target sequences with a fraction of Ns greater than $maxN ";
	my $newNumSeqs =0;
	my $lenZero = 0;
	my $tooManyNs = 0;
	foreach(keys %$seqs) {
		my $x = $_;
		my $s = $seqs->{$x};
		my $L1 = length($s);
		if ($L1 < 1) {
			$lenZero++;
			delete $seqs->{$x};
			next;
		}
		$s =~ s/[Nn]//g;
		my $L2 = length($s);
		if (($L1-$L2)/$L1 > $maxN) {
			$tooManyNs++;
			delete $seqs->{$x};
			next;
		}
		$newNumSeqs++;
	}
	print STDERR " ($newNumSeqs of $numSeqs kept)\n";
	$numSeqs = $newNumSeqs;
}	
	
if ($numSeqs < 1) {
	print STDERR "\tProblem - number of sequences is less than 1!!!\n";
	exit;
}
print STDERR "\n\t$numSeqs target sequences extracted for analysis\n";


if (scalar(@regions) > 0) {
	print STDERR "\tWill analyze these regions for enrichment:\n";
	foreach(@regions) {
		print STDERR "\t\t$_->{'s'} to $_->{'e'}\n";
	}
} elsif (($pstart ne '' || $pfixed ne '') && ($wstart ne '' || $wfixed ne '')) {
	print STDERR "\tWill analyze a range of regions:\n";
	print STDERR "\t\tpstart=$pstart, pend=$pend, pfixed=$pfixed, pstep=$pstep\n";
	print STDERR "\t\twstart=$wstart, wend=$wend, wfixed=$wfixed, wstep=$wstep\n";
} elsif ($allFlag) {
	print STDERR "\tWill exhaustively check enrichment by position and window\n";
	@windows = ();
} elsif ($mode eq 'denovo') {
	my $s = $offset;
	my $e = $offset+$size;
	print STDERR "\tWill search for de novo motifs in given sequence from $s to $e\n";
	my $x = {s=>$s,e=>$e};
	push(@regions, $x);
}
if ($mode eq 'denovo' && $strand eq 'separate') {
	print STDERR "!!! Error - please use -strand with either +, -, or both when doing denovo enrichment\n";
	exit;
}

if ($calcSeqPropertiesFlag) {
	my $targetSeqPropertiesFile = $outputDirectory . "/targetSeqProperties.tsv";
	my $targetSeqProperties = calcSeqProperties($seqs, $seqPropertyOrder,$targetSeqPropertiesFile);
}

if ($numRandSeq < 0) {
	$numRandSeq = $numSeqs;
}
if ($bgSeqFile eq '' || $bgPosFile ne '') {
	print STDERR "\tWill use $numRandSeq background sequences for significance calculation\n";
}

# Generate Background Sequences

print STDERR "\n\tExtracting background sequences\n";


my $bgSeqFileToUse = $outputDirectory . "/background.fa";

if ($bgSeqFile ne '') {
	print STDERR "\tUsing background sequences provided by the user (file=$bgSeqFile)\n";
	unless (-e $bgSeqFile) {
		print STDERR "!!! Error: Could not open background FASTA file: $bgSeqFile\n";
		exit;
	}
	`cp "$bgSeqFile" "$bgSeqFileToUse"`;

} elsif ($bgPosFile ne '') {

	if ($ogSize ne '') {
		`adjustPeakFile.pl "$bgPosFile" -size $ogSize > "$bgPosResized"`;
		`homerTools extract "$bgPosResized" "$genomePath" $mask -fasta > "$bgSeqFileToUse"`;
	} else {
		`cp "$bgPosFile" "$bgPosResized"`;
		`homerTools extract "$bgPosResized" "$genomePath" $mask -fasta > "$bgSeqFileToUse"`;
	}

} else {
	# Ok - use the genome and/or a model to create background sequences.

	print STDERR "\tUsing 'homer2 bg' background selection (v5.0+)\n";
	
	my $opt = " -o \"$outputDirectory/homerBg\" ";
	$opt .= " $mask";
	$opt .= " -model" if ($modelFlag);
	$opt .= " -NN $NN -N $numRandSeq ";
	$opt .= " -nbins $numBins -nsubBins $nsubBins ";
	if ($positionIndependentModelFlag) {
		$opt .= " -ikmer $kmer ";
	} else {
		$opt .= " -pkmer $kmer ";
	}
	$opt .= " -allowTargetOverlaps" if ($allowTargetOverlap);
    $opt .= " -allowBgOverlaps" if ($allowBgOverlap);
	my $inputs = "";
	if ($tgSeqFile ne '') {
		$inputs = " -i \"$tgSeqFile\"";
	} else {
		$inputs = "-g \"$genomePath\" -p \"$targetPosResized\"";
	}
	print STDERR "\tExecuting homer2 bg:\n\t\thomer2 bg $inputs $opt\n\n";
	open OUTX, ">>$outputDirectory/cmd.txt" or die "!!! Error: Could not open file \"$$outputDirectory/cmd.txt\" for writing\n";
	print OUTX "homer2 cmd: homer2 bg $inputs $opt\n";
	close OUTX;
	`homer2 bg $inputs $opt`;
	`cp $outputDirectory/homerBg.bg.sequences.fasta $bgSeqFileToUse`;
}


$randSeqs = HomerSeqBias::readFASTA($bgSeqFileToUse);
$numRandSeq = 0;
my $wrongSizeCount = 0;
my $tooManyNsCount = 0;
foreach(keys %$randSeqs) {
	my $id = $_;
	my $s = $randSeqs->{$id};
	my $L = length($s);
	unless (($L == $size || $L == $size+1) && $L > 0) {
		$wrongSizeCount++;
		delete $randSeqs->{$id};
		print STDERR "\t\t!!! Warning - background sequences are not the same length as the target sequences ($L != $size), seqID $_\n";
		next;
	}
	$s =~ s/[Nn]//g;
	my $L2 = length($s);
	if (($L-$L2)/$L > $maxN) {
		$tooManyNsCount++;
		delete $randSeqs->{$id};
		next;
	}
	$numRandSeq++;
}
print STDERR "\tInput background sequences read from file $bgSeqFileToUse to use for analysis\n";
print STDERR "\t\tN = $numRandSeq ($tooManyNsCount had too many Ns, $wrongSizeCount were the wrong length)\n";

if ($calcSeqPropertiesFlag) {
	my $backgroundSeqPropertiesFile = $outputDirectory . "/backgroundSeqProperties.tsv";
	my $backgroundSeqProperties = calcSeqProperties($randSeqs, $seqPropertyOrder,$backgroundSeqPropertiesFile);
}

my $numTgSeq = $numSeqs;
my $numBgSeq = $numRandSeq;
my $pseudoFreq = 1/($numTgSeq+$numBgSeq);

if (scalar(@regions) < 1) {
	if (scalar(@windows) < 1) {
		if ($wfixed ne '') {
			$wstart = $wfixed;
			$wend = $wfixed;
			$wstep = 1;
		}
		if ($wstart eq '') {
			$wstart = 1;
		}
		if ($wend eq '') {
			$wend = $size;
		}
		if ($wstep eq '') {
			$wstep = 1;
		}
		print STDERR "\tWindow sizes of enrichment intervals to check:\n";
		print STDERR "\t\twstart:$wstart\n";
		print STDERR "\t\twend  :$wend\n";
		print STDERR "\t\twstep :$wstep\n";
		for (my $i=$wstart;$i<=$wend;$i+=$wstep) {
			push(@windows, $i);
		}
	} else {
		print STDERR "\tWill check the following window sizes for enrichment intervals:\n";
		@windows = sort {$a <=> $b} @windows;
		foreach(@windows) {
			print STDERR "\t\t$_\n";
		}
	}
	if (scalar(@positions) < 1) {
		if ($pfixed ne '') {
			$pstart = $pfixed;
			$pend = $pfixed;
			$pstep = 1;
		}
		if ($pstart eq '') {
			$pstart = $offset;
		}
		if ($pend eq '') {
			$pend = $offset+$size;
		}
		if ($pstep eq '') {
			$pstep = 1;
		}
		print STDERR "\tPositions (center of windows) where enrichment calculations will be made:\n";
		print STDERR "\t\tpstart:$pstart\n";
		print STDERR "\t\tpend  :$pend\n";
		print STDERR "\t\tpstep :$pstep\n";
		for (my $i=$pstart;$i<=$pend;$i+=$pstep) {
			push(@positions, $i);
		}
	} else {
		print STDERR "\tWill check the following positions when performing enrichment calculations:\n";
		@positions = sort {$a <=> $b} @positions;
		foreach(@positions) {
			print STDERR "\t\t$_\n";
		}
	}
} else {
	print STDERR "\tWill analyze the specified regions:\n";
	foreach(@regions) {
		my $s = $_->{'s'};
		my $e = $_->{'e'};
		print STDERR "\t\t$s,$e\n";
	}
}





if ($mode eq 'denovo') {
	`mkdir -p $outputDirectory/denovo`;

	my @denovoRegions = ();
	my @wpRegions = ();
	if (@regions > 0) {
		foreach(@regions) {
			if ($enrichedFlag) {
				my $r = {s=>$_->{'s'},e=>$_->{'e'},t=>'enriched'};
				push(@denovoRegions, $r);
				$_->{'enriched'} = $r;
			}
			if ($depletedFlag) {
				my $r = {s=>$_->{'s'},e=>$_->{'e'},t=>'depleted'};
				push(@denovoRegions, $r);
				$_->{'depleted'} = $r;
			}
		}
	} else {
		foreach(@windows) {
			my $w = $_;
			foreach(@positions) {
				my $p = $_;
				$wpRegions{$w . "_" . $p} = ();
				my $halfSizeLow = floor($w/2);
				my $halfSizeHigh = $halfSizeLow;
				$halfSizeHigh-- if ($w % 2== 0);
				my $s = $p-$halfSizeLow;
				my $e = $p+$halfSizeHigh;
				if ($enrichedFlag) {
					my $r = {s=>$s,e=>$e,t=>'enriched'};
					push(@denovoRegions, $r);
					$wpRegions{$w . "_" . $p}->{'enriched'} = $r;
				}
				if ($depletedFlag) {
					my $r = {s=>$s,e=>$e,t=>'depleted'};
					push(@denovoRegions, $r);
					$wpRegions{$w . "_" . $p}->{'depleted'} = $r;
				}
			}
		}
	}
	foreach(@denovoRegions) {
		my $s = $_->{'s'};
		my $e = $_->{'e'};
		my $t = $_->{'t'};
		print STDERR "\t\t$s\t$e\t$t\n";
	}

	foreach(@denovoRegions) {
		my $s = $_->{'s'};
		my $e = $_->{'e'};
		my $t = $_->{'t'};
		my $f = "motifs_$s" . "_$e" . ".$t.motifs";;
		my $outFile = "$outputDirectory/denovo/" . $f;
		$_->{'f'} = $outFile;

		my $pid = fork();
		$cpu++;
		if ($pid==0) {
			if ($t eq 'enriched') {
				findDenovoMotifs($outFile, $seqs, $randSeqs, $s, $e, $strand, $denovoLength);
			} else {
				findDenovoMotifs($outFile, $randSeqs, $seqs, $s, $e, $strand, $denovoLength);
			}
			exit;
		}
		if ($cpu >= $maxCPUs) {
			my $id = wait();
			$cpu--;
		}
	}
	my $id =0;
	while ($id>=0) {
		$id = wait();
	}

	open HTML, ">$outputDirectory/denovo/index.html";
	print HTML "<HTML><BODY>\n";
	if (@regions > 0) {
		print HTML "<H2>Regions Analyzed for de novo motif enrichment:</H2>\n";
		print HTML "<TABLE><TR><TH>Region</TH>\n";
		print HTML "\t<TH>Enrichment</TH><TH>Top Enriched</TH><TH>Top Enriched LogP</TH><TH>Top Match</TH>\n";
		print HTML "\t<TH>Depleted</TH><TH>Top Depleted</TH><TH>Top Depleted LogP</TH><TH>Top Match</TH>\n";
		print HTML "</TR>\n";
		foreach(@regions) {
			my $rr = $_;
			my $s = $rr->{'s'};
			my $e = $rr->{'e'};

			my $outFile1 = "";
			my $outFileLink1 = "";
			my $resultsDir1 = "";
			my $motifs1 = "";
			my $cons1 = "NA";
			my $match1 = "";
			my $logp1 = "NA";
			my $logo1 = '';
			if (exists($rr->{'enriched'})) {
				$outFile1 = $rr->{'enriched'}->{'f'};
				$outFileLink1 = $outFile1;
				$outFileLink1 =~ s/^.+\///;
				$resultsDir1 = "$outFileLink1" . ".results/homerResults.html";
				$motifs1 = HomerConfig::readMotifFile($outFile1);
				if (@$motifs1>0) {
					$logp1 = $motifs1->[0]->{'logp'};
					$cons1 = $motifs1->[0]->{'consensus'};
					$cons1 =~ s/^\>//;
					$match1 = parseBestGuess($motifs1->[0]->{'name'});
					$logo1 = HomerSVGLogo::motif2svg($motifs1->[0], $bitFlag, $rnaFlag);
				}
			}

			my $outFile2 = "";
			my $outFileLink2 = "";
			my $resultsDir2 = "";
			my $motifs2 = "";
			my $cons2 = "NA";
			my $match2 = "NA";
			my $logp2 = "NA";
			my $logo2 = '';
			if (exists($rr->{'depleted'})) {
				$outFile2 = $rr->{'depleted'}->{'f'};
				$outFileLink2 = $outFile2;
				$outFileLink2 =~ s/^.+\///;
				$resultsDir2 = "$outFileLink2" . ".results/homerResults.html";
				$motifs2 = HomerConfig::readMotifFile($outFile2);
				if (@$motifs2>0) {
					$logp2 = $motifs2->[0]->{'logp'};
					$cons2 = $motifs2->[0]->{'consensus'};
					$cons2 =~ s/^\>//;
					$match2 = parseBestGuess($motifs2->[0]->{'name'});
					$logo2 = HomerSVGLogo::motif2svg($motifs2->[0], $bitFlag, $rnaFlag);
				}
			}


			print HTML "<TR><TD>$s to $e</TD>";
			print HTML "<TD><A HREF=\"$resultsDir1\">results</A>|<A HREF=\"$outFileLink1\">motif file</A></TD>\n";
			print HTML "<TD>$logo1</TD>";
			print HTML "<TD>$logp1</TD>";
			print HTML "<TD>$match1</TD>";
			print HTML "<TD><A HREF=\"$resultsDir2\">results</A>|<A HREF=\"$outFileLink2\">motif file</A></TD>\n";
			print HTML "<TD>$logo2</TD>";
			print HTML "<TD>$logp2</TD>";
			print HTML "<TD>$match2</TD>";
			print HTML "</TR>\n";

		}
		print HTML "</TABLE>\n";
	} else {
		#for (my $w = $wstart;$w<=$wend;$w+=$wstep) {
		for (my $windex =0;$windex < scalar(@windows); $windex++) {
			
			my $w = $windows[$windex];
			open TABLE, ">$outputDirectory/denovo/window_$w.tsv";
			print HTML "<H2>Analyzing with window size of $w</H2>\n";

			print HTML "<TABLE><TR><TH>Position</TH><TH>Region</TH>\n";
			print HTML "\t<TH>Enrichment</TH><TH>Top Enriched Motif</TH><TH>Top Enriched LogP</TH><TH>Best Match</TH>\n";
			print HTML "\t<TH>Depleted</TH><TH>Top Depleted</TH><TH>Top Depleted LogP</TH><TH>Best Match</TH></TR>\n";
			#for (my $p = $pstart;$p<=$pend;$p+=$pstep) {
			for (my $pindex = 0;$pindex<scalar(@positions);$pindex++) {
				my $p = $positions[$pindex];
			
				my $halfSizeLow = floor($w/2);
				my $halfSizeHigh = $halfSizeLow;
				$halfSizeHigh-- if ($w % 2== 0);
				my $s = $p-$halfSizeLow;
				my $e = $p+$halfSizeHigh;

				my $code = $w . "_" . $p;

				my $outFile1 = "";
				my $outFileLink1 = "";
				my $resultsDir1 = "";
				my $motifs1 = "";
				my $cons1 = "NA";
				my $match1 = "NA";
				my $logp1 = "NA";
				my $logo1 = '';
				if (exists($wpRegions{$code}->{'enriched'})) {
					$outFile1 = $wpRegions{$code}->{'enriched'}->{'f'};
					$outFileLink1 = $outFile1;
					$outFileLink1 =~ s/^.+\///;
					$resultsDir1 = "$outFileLink1" . ".results/homerResults.html";
					$motifs1 = HomerConfig::readMotifFile($outFile1);
					if (@$motifs1>0) {
						$logp1 = $motifs1->[0]->{'logp'};
						$cons1 = $motifs1->[0]->{'consensus'};
						$cons1 =~ s/^\>//;
						$match1 = parseBestGuess($motifs1->[0]->{'name'});
						$logo1 = HomerSVGLogo::motif2svg($motifs1->[0], $bitFlag, $rnaFlag);
					}
				}

				my $outFile2 = "";
				my $outFileLink2 = "";
				my $resultsDir2 = "";
				my $motifs2 = "";
				my $cons2 = "NA";
				my $match2 = "NA";
				my $logp2 = "NA";
				my $logo2 = '';
				if (exists($wpRegions{$code}->{'depleted'})) {
					$outFile2 = $wpRegions{$code}->{'depleted'}->{'f'};
					$outFileLink2 = $outFile2;
					$outFileLink2 =~ s/^.+\///;
					$resultsDir2 = "$outFileLink1" . ".results/homerResults.html";
					$motifs2 = HomerConfig::readMotifFile($outFile2);
					if (@$motifs2>0) {
						$logp2 = $motifs2->[0]->{'logp'};
						$cons2 = $motifs2->[0]->{'consensus'};
						$cons2 =~ s/^\>//;
						$match2 = parseBestGuess($motifs2->[0]->{'name'});
						$logo2 = HomerSVGLogo::motif2svg($motifs2->[0], $bitFlag, $rnaFlag);
					}
				}	

				print HTML "<TR><TD>$p</TD><TD>$s to $e</TD>";
				print HTML "<TD><A HREF=\"$resultsDir1\">results</A>|<A HREF=\"$outFileLink1\">motif file</A></TD>\n";
				print HTML "<TD>$logo1</TD>";
				print HTML "<TD>$logp1</TD>";
				print HTML "<TD>$match1</TD>";
				print HTML "<TD><A HREF=\"$resultsDir2\">results</A>|<A HREF=\"$outFileLink2\">motif file</A></TD>\n";
				print HTML "<TD>$logo2</TD>";
				print HTML "<TD>$logp2</TD>";
				print HTML "<TD>$match2</TD>";
				print HTML "</TR>\n";

				print TABLE "$p\t$s,$e\t$logp1\t$cons1\t$logp2\t$cons2\n";
				
			}
			print HTML "</TABLE>\n";
			close TABLE;
		}
	}
	print HTML "<P>Command: $cmd</P>\n";
	print HTML "</BODY></HTML>\n";

	print STDERR "\n\tFinished with denovo motif finding\n";

} elsif ($mode eq 'known') {

	my %targetSites = ();
	print STDERR "\tScanning target sequences for motifs:\n";
	scanMotifsPosition(\%targetSites,$seqs,$motifFile,$offset,$strand);
	my %bgSites = ();
	print STDERR "\tScanning background sequences for motifs:\n";
	scanMotifsPosition(\%bgSites,$randSeqs,$motifFile,$offset,$strand);

	open LOGPALL, ">$outputDirectory/allmotifs.logp.txt";
	if (scalar(@regions) >0) {
		print LOGPALL "Motif";
		foreach(@regions) {
			my $s = $_->{'s'};
			my $e = $_->{'e'};
			print LOGPALL "\t$s,$e";
		}
		print LOGPALL "\n";
	} else {
		print LOGPALL "Motif\tWindowSize";
		#for (my $j=$pstart;$j<=$pend;$j+=$pstep) {
		foreach(@positions) {
			print LOGPALL "\t$_";
		}
		print LOGPALL "\n";
	}

	open INTERVALS, ">$outputDirectory/summary.bestIntervals.txt";
	print INTERVALS "Motif\tMost Enriched Start\tMost Enriched End\tEnrichment logp\tMost Depleted Start\tMost Depleted End\tDepletion logp\n";

	my %results = ();
	my %cleanNames = ();
	my %powerFiles = ();
	
	# this loop goes through each motif to output 
	my @motifs = sort {$a cmp $b} keys %targetSites;
	foreach(@motifs) {
		my $mname = $_;
		if (!exists($results{$mname})) {
			my %a = ();
			$results{$mname} = \%a;
		}
		my $clean = cleanMotifName($mname);
		if (exists($cleanNames{$clean})) {
			my $x = 2;
			my $n = $clean . "_" . $x;
			while (exists($cleanNames{$n})) {
				$x++;
				$n = $clean . "_" . $x;
			}
			$clean = $n;
		}
		$cleanNames{$clean}=1;
	
		my @files = ();
		my $cpu = 0;
		
		print STDERR "\tAnalyzing $clean\n";

		# Start with outputing the basic frequency counts
		open FREQ, ">$outputDirectory/$clean.freq.txt";
		print FREQ "$clean ($mname)\tTarget Frequency\tBackground Frequency\tLog2Ratio\n";
		my $tgFreq = getPositionCounts(\%targetSites, $mname,$offset,$size);
		my $bgFreq = getPositionCounts(\%bgSites, $mname,$offset,$size);
		for (my $i=$offset;$i<=$offset+$size;$i++) {
			my $v1 = 0;
			if (exists($tgFreq->{$i})) {
				$v1 = $tgFreq->{$i};
				$v1 /= $numTgSeq;
			}
			my $v2 = 0;
			if (exists($bgFreq->{$i})) {
				$v2 = $bgFreq->{$i};
				$v2 /= $numBgSeq;
			}
			my $ratio = log(($v1+$pseudoFreq)/($v2+$pseudoFreq))/log(2.0);
			print FREQ "$i\t$v1\t$v2\t$ratio\n";
		}
		close FREQ;

		# perform fourier analysis on region that normally exhibits helical density changes relative to TSS
		if ($fourierFlag) {
			if ($offset > $fourierStart) {
				$fourierStart = $offset;
			}
			if ($offset+$size < $fourierEnd) {
				$fourierEnd = $offset;
			}
			if ($fourierStart <= $fourierEnd+10) {
				my $fout = "$outputDirectory/$clean.fourier.$fourierStart" . 'to' . "$fourierEnd.txt";
				`powerSpec.pl "$outputDirectory/$clean.freq.txt" -col 4 -start $fourierStart -end $fourierEnd > "$fout"`;
				$powerFiles{$clean} = $fout;
			} else {
				print STDERR "\tSkiping fourier analysis - region outside of fourier range ($fourierStart,$fourierEnd)\n";
			}
		}

		if (scalar(@regions) > 0) {

			for (my $i=0;$i<@regions;$i++) {

				my $s = $regions[$i]->{'s'};
				my $e = $regions[$i]->{'e'};
				print STDERR "\t\tregion $s,$e ($clean)\n";
			#for (my $i=$size;$i>0;$i--) {
	
				my $resultFile = $tmpFile . "_region_" . $s . "_" . $e;
				push(@files, $resultFile);
	
				my $pid = fork();
				$cpu++;
				if ($pid==0) {
		
					open RES, ">$resultFile";
		
					my $ntg = 0;
					my $nbg = 0;
					my $tgRatio = 0;
					my $bgRatio = 0;
					my $v = 0;
					if ($s >= $offset && $e <= $offset+$size) {
						$ntg = getSeqCountWithMotif(\%targetSites,$mname,$s,$e);
						$nbg = getSeqCountWithMotif(\%bgSites,$mname,$s,$e);
						$tgRatio = $ntg/$numTgSeq;
						$bgRatio = $nbg/$numBgSeq;
						if ($tgRatio > $bgRatio) {
							$v = Statistics::loghypergeo($numTgSeq+$numBgSeq,$numTgSeq,$ntg+$nbg,$ntg);
						} else {
							$v = Statistics::iloghypergeo($numTgSeq+$numBgSeq,$numTgSeq,$ntg+$nbg,$ntg);
							$v *= -1;
						}
						$results{$mname}->{$s . "_" . $e} = $v;
						print RES "$mname\t$s\t$e\t$ntg\t$nbg\t$v\n";
						print STDERR "$mname\t$s\t$e\t$ntg\t$nbg\t$v\n";
						#print STDERR "\t$mname\tSize=$w,offset=$p ($s,$e)\t$ntg\t$nbg\t$v\n";
					}
					close RES;
					exit;
				}
				if ($cpu >= $maxCPUs) {
					my $id = wait();
					$cpu--;
				}
			}
			my $id =0;
			while ($id>=0) {
				$id = wait();
			}
	
			foreach(@files) {
				my $file =$_;
				open IN, $file;
				while (<IN>) {
					chomp;
					s/\r//g;
					my @line = split /\t/;
					my $mname = $line[0];
					my $i = $line[1];
					my $j = $line[2];
					my $v = $line[5];
					$results{$mname}->{$i . "_" . $j} = $v;
				}
				close IN;
				`rm "$file"`;
			}

			my $mname = $_;
			print LOGPALL "$mname";
	
			foreach(@regions) {
				my $s = $_->{'s'};	
				my $e = $_->{'e'};	
				my $v = 0;
				if (exists($results{$mname}->{$s . "_" . $e})) {
					$v = $results{$mname}->{$s . "_" . $e};
				}
				print LOGPALL "\t$v";
			}
			print LOGPALL "\n";

		} else {
			#for (my $w = $wstart;$w<=$wend;$w+=$wstep) {
			#for (my $i=$size;$i>0;$i--) {
			foreach(@windows) {
				my $w = $_;
				print STDERR "\t\tWindow = $w ($clean)\n";

				foreach(@positions) {
					my $p = $_;
					#print STDERR "\t\t\tposition = $p ($clean)\n";
					my $halfSizeLow = floor($w/2);
					my $halfSizeHigh = $halfSizeLow;
					$halfSizeHigh-- if ($w % 2== 0);
					$s = $p-$halfSizeLow;
					$e = $p+$halfSizeHigh;
					my $ntg = 0;
					my $nbg = 0;
					my $tgRatio = 0;
					my $bgRatio = 0;
					my $v = 0;

					if ($s >= $offset && $e <= $offset+$size) {
						$ntg = getSeqCountWithMotif(\%targetSites,$mname,$s,$e);
						$nbg = getSeqCountWithMotif(\%bgSites,$mname,$s,$e);
						$tgRatio = $ntg/$numTgSeq;
						$bgRatio = $nbg/$numBgSeq;
						if ($tgRatio > $bgRatio) {
							$v = Statistics::loghypergeo($numTgSeq+$numBgSeq,$numTgSeq,$ntg+$nbg,$ntg);
						} else {
							$v = Statistics::iloghypergeo($numTgSeq+$numBgSeq,$numTgSeq,$ntg+$nbg,$ntg);
							$v *= -1;
						}
						$results{$mname}->{$w . "_" . $p} = $v;
						#print RES "$mname\t$w\t$p\t$ntg\t$nbg\t$v\n";
						#print STDERR "\t$mname\tSize=$w,offset=$p ($s,$e)\t$ntg\t$nbg\t$v\n";
					}
				}
			}
	
			foreach(@files) {
				my $file =$_;
				open IN, $file;
				while (<IN>) {
					chomp;
					s/\r//g;
					my @line = split /\t/;
					my $mname = $line[0];
					my $i = $line[1];
					my $j = $line[2];
					my $v = $line[5];
					$results{$mname}->{$i . "_" . $j} = $v;
				}
				close IN;
				`rm "$file"`;
			}
		
		
			open LOGP, ">$outputDirectory/$clean.logp.txt";
		
			my $maxStr = '';
			my $maxStr2 = '\tNA\tNA\tNA';
			my $maxVal = -1e99;
			my $minStr = '';
			my $minStr2 = '\tNA\tNA\tNA';
			my $minVal = 1e99;

			print LOGP "$mname";
			#for (my $j=$pstart;$j<=$pend;$j+=$pstep) {
			foreach(@positions) {
				print LOGP "\t$_";
			}
			print LOGP "\n";
			#for (my $i = $wstart;$i<=$wend;$i+=$wstep) {
			foreach(@windows) {
				my $i=$_;
				print LOGP "$i";
				print LOGPALL "$mname\t$i";
				#for (my $j=$pstart;$j<=$pend;$j+=$pstep) {
				foreach(@positions) {
					my $j = $_;
					my $v = 0;
					if (exists($results{$mname}->{$i . "_" . $j})) {
						$v = $results{$mname}->{$i . "_" . $j};
					}
					if ($v ne '0') {
						my $ss = $j-floor($i/2);
						my $ee = $j+floor($i/2);
						if ($v > $maxVal) {
							$maxVal = $v;
							$maxStr = "tf=$clean\tposition=$j\twindow=$i\tlogp=$v";
							$maxStr2 = "\t$ss\t$ee\t$v";
						}
						if ($v < $minVal) {
							$minVal = $v;
							$minStr = "tf=$clean\tposition=$j\twindow=$i\tlogp=$v";
							$minStr2 = "\t$ss\t$ee\t$v";
						}
					}
					print LOGP "\t$v";
					print LOGPALL "\t$v";
				}
				print LOGP "\n";
				print LOGPALL "\n";
			}
			close LOGP;
		
			print STDERR "\n\tMin: $minStr\n";
			print STDERR "\tMax: $maxStr\n";
			print INTERVALS "$clean" . $minStr2 . $maxStr2 . "\n";
		
			open STATS, ">$outputDirectory/$clean.stats.txt";
			print STATS "$clean ($mname)\n";
			print STATS "cmd=$cmd\n";
			print STATS "Number of target sequences: $numTgSeq\n";
			print STATS "Number of background sequences: $numBgSeq\n";
			print STATS "Minimum Enrichment Interval:\t$minStr\n";
			print STATS "Maximum Enrichment Interval:\t$maxStr\n";
			close STATS;
	
		}
	}
	close LOGPALL;
	close INTERVALS;

	if (scalar(@windows) < 10 || $forceSummaryByWindow) {
		my %windowFiles = ();
		foreach(@windows) { 
			my $w = $_;
			my $f= "$outputDirectory/summary.window$w.logp.txt";
			open my $fh, ">$f";
			print $fh "Window=$w";
			$windowFiles{$w} = $fh;
		}
		open IN, "$outputDirectory/allmotifs.logp.txt";
		my $z = 0;
		while (<IN>) {
			$z++;
			chomp;
			s/\r//g;
			my @line = split /\t/;
			if ($z == 1) {
				my $header = "";
				for (my $i=2;$i<@line;$i++) {
					$header .= "\t$line[$i]";
				}
				foreach(values %windowFiles) {
					print {$_} $header . "\n";
				}
				next;
			}
			my $w = $line[1];
			if (exists($windowFiles{$w})) {
				my $fh = $windowFiles{$w};
				print $fh "$line[0]";
				for (my $i=2;$i<@line;$i++) {
					print $fh "\t$line[$i]";
				}
				print $fh "\n";
			}
		}
		close IN;
		foreach(@windows) { 
			my $w = $_;
			my $fh = $windowFiles{$w};
			close $fh;
		}

		foreach(@windows) { 
			my $w = $_;
			my $f1= "$outputDirectory/summary.window$w.logp.txt";
			my $f2= "$outputDirectory/summary.window$w.logq.txt";
			`adjustPvalueMatrix.pl -i "$f1" -logp -inbase e -outbase 10 -flip > "$f2"`;
		}
	}


	# summarize power spectrum results for all motifs, normalizing each one by their maximum signal
	my @pmotifs = keys %powerFiles;
	if (scalar(@pmotifs) > 0) {
		open POWER, ">$outputDirectory/summary.fourier.txt";
		print POWER "Period (Fourier analysis, max of power signals normalized to 1)";
		my %powerData = ();
		my %periods = ();
		for (my $z=0;$z<@pmotifs;$z++) {
			my $mname = $pmotifs[$z];
			print POWER "\t$mname";
			my $f = $powerFiles{$mname};
			open IN, $f;
			my $z = 0;
			my $max = -1e99;
			while (<IN> ){
				$z++;
				next if ($z < 2);
				chomp;
				s/\r//g;
				my @line = split /\t/;
				$periods{$line[0]}=1;
				if (!exists($powerData{$mname})) {
					my %a = ();
					$powerData{$mname} = \%a;
				}
				$powerData{$mname}->{$line[0]} = $line[1];
				if ($line[1] > $max) {
					$max = $line[1];
				}
			}
			close IN;
			if ($max == 0) {
				$max = 1;
			}
			foreach(values %{$powerData{$mname}}) {
				$_ /= $max;
			}
		}
		print POWER "\n";
		my @periods = sort {$a <=> $b} keys %periods;
		foreach(@periods) {
			my $p = $_;
			print POWER "$p";
			foreach(@pmotifs) {
				my $m = $_;
				my $v = 0;
				if (exists($powerData{$m}->{$p})) {
					$v = $powerData{$m}->{$p};
				}
				print POWER "\t$v";
			}
			print POWER "\n";
		}
		close POWER;
	}
	#`rm $tmpFile $tmpFile2`;
}
exit;


#Then scan at each position
sub getExpectedMotifFrequencies {
	my ($data,$seqs,$motifFile,$start,$end,$markov) = @_;
	
	my $set = "w=$window,m=$markov";
	#First generate sequence associated with that position
	#

	my $cpus = 0;

	for (my $p=$start;$p<=$end;$p++) {

		my $pid = fork();
		$cpus++;
		if ($pid==0) {

			my $forkprefix = $rand . ".pos" . $p;

			my $eWindow = $window;
			$eWindow = $markov+1 if ($markov+1 > $window);
	
			my $pp = $p-$offset;
			my $x = floor(($window-1)/2);
			my $s = $pp-$x;
		
			print STDERR "\tAnalyzing $p (w=$window,ew=$eWindow,m=$markov)\n";	
			open OUT, ">$forkprefix.in.fa";
			my $z = 1;
			foreach(values %$seqs) {
				my $s = substr($_,$s,$eWindow);
				print OUT ">S$z\n$s\n";
				$z++;
			}
			close OUT;
	
			`$seqBackgroundEXE --fa $forkprefix.in.fa --prefix $forkprefix --N $expN --slen $expLen --upper --order $markov`;
			`homer2 find -i $forkprefix.fa -m $motifFile > $forkprefix.sites.txt`;
			`rm $forkprefix.in.fa $forkprefix.fa $forkprefix.df.pkl`;

			exit(0);
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

	for (my $p=$start;$p<=$end;$p++) {

		my $forkprefix = $rand . ".pos" . $p;

		open IN, "$forkprefix.sites.txt";
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			my $mname = $line[3];
			my $d = $line[4];
			if ($d eq '+') {
				$d = 'p';
			} else {
				$d = 'n';
			}
			if (!exists($data->{$mname})) {
				my %a = ();
				my %b = ();
				$data->{$mname} = {p=>\%a,n=>\%b};
				$mlen{$mname} = length($line[2]);
			}
			if (!exists($data->{$mname}->{$d}->{$p})) {
				my %a = ();
				$data->{$mname}->{$d}->{$p}=\%a;
				$data->{$mname}->{$d}->{$p}->{$set}=0;
			}
			$data->{$mname}->{$d}->{$p}->{$set}++;
		}
		close IN;
		`rm $forkprefix.sites.txt`;
	}
}

#first scan in real sequences
sub scanMotifsPosition {
	my ($data,$seqs,$motifFile,$offset,$strand) = @_;

	open OUT, ">$tmpFile";
	foreach(keys %$seqs) {
		my $id = $_;
		my $seq = $seqs->{$id};
		print OUT "$id\t$seq\n";
	}
	close OUT;

	my $strandOpt = '';
	if ($strand eq '+' ){
		$strandOpt = ' -strand + ';
 	} elsif ($strand eq '-' ){
		$strandOpt = ' -strand - ';
	}
	`homer2 find -s $tmpFile -m $motifFile -p $maxCPUs -offset $offset $strandOpt > $tmpFile2`;
	
	open IN, $tmpFile2;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $seqID = $line[0];
		my $mname = $line[3];
		my $p = $line[1]+floor(length($line[2])/2);
		my $d = $line[4];
		if ($d eq '+') {
			$d = 'p';
		} else {
			$d = 'n';
			$p = $line[1]-floor(length($line[2])/2);
		}
		if ($strand eq 'separate') {
			$mname=$mname . "_$d";
		}
		if (!exists($data->{$mname})) {
			my %a = ();
			$data->{$mname} = \%a;
			$mlen{$mname} = length($line[2]);
		}
		if (!exists($data->{$mname}->{$seqID})) {
			my @a = ();
			$data->{$mname}->{$seqID}=\@a;
		}
		push(@{$data->{$mname}->{$seqID}},$p);
	}

	#finally, sort the motifs positions
	foreach(keys %$data) {
		my $mname = $_;
		foreach(keys %{$data->{$mname}}) {
			my $seqID = $_;
			my @a = sort {$a <=> $b} @{$data->{$mname}->{$seqID}};
			$data->{$mname}->{$seqID} = \@a;
		}
	}
	`rm "$tmpFile" "$tmpFile2"`;
}

sub getSeqCountWithMotif {
	my ($targetSites,$mname,$s,$e) = @_;
	my $n = 0;
	foreach(values %{$targetSites->{$mname}}) {
		my $x=0;
		foreach(@$_) {
			if ($_ >=$s && $_<=$e) {
				$x=1;
				last;
			}
			last if ($_ > $e);
		}	
		$n += $x;
	}
	return $n;
}
sub getPositionCounts {
	my ($targetSites, $mname,$s,$e) = @_;
	my %hist = ();
	foreach(values %{$targetSites->{$mname}}) {
		foreach(@$_) {
			if ($_ >=$s && $_<=$e) {
				$hist{$_}++;
			}
			last if ($_ > $e);
		}
	}
	return \%hist;
}

sub cleanMotifName {
	my ($name) = @_;
	my $mname = $name;
	my $s = "";
	if ($mname =~ /_+$/) {
		$s = "+";
	} elsif ($mname =~ /_-$/) {
		$s = "-";
	}
	$mname =~ s/\/.*//;
	$mname =~ s/\,//g;
	$mname =~ s/\(/_/g;
	$mname =~ s/\)/_/g;
	$mname =~ s/\s/_/g;
	$mname =~ s/\_+$//;
	if ($s ne '') {
		$mname .= "_$s";
	}
	#print STDERR "Cleaned \"$name\" -> \"$mname\"\n";
	return $mname;
}

sub findDenovoMotifs {
	my ($outFile, $seqs, $randSeqs, $start, $end, $strand, $len) = @_;

	my $s = $start-$offset;
	my $l = $end-$start+1;

	my $tmpSeqFile = "$outFile.seq.tsv";
	my $tmpGroupFile = "$outFile.group.tsv";
	open SEQ, ">$tmpSeqFile";
	open GROUP, ">$tmpGroupFile";
	foreach(keys %$seqs) {
		my $id =$_;
		my $seq = substr($seqs->{$id},$s,$l);
		print SEQ "$id\t$seq\n";
		print GROUP "$id\t1\n";
	}
	foreach(keys %$randSeqs) {
		my $id =$_;
		my $seq = substr($randSeqs->{$id},$s,$l);
		print SEQ "$id\t$seq\n";
		print GROUP "$id\t0\n";
	}
	close SEQ;
	close GROUP;
	
	my $opts = $homerOpts;
	$opts .= " -strand $strand ";
	$opts .= " -p $maxCPUs2 ";
	$opts .= " -S $denovoS ";

	my @lens = split /\,/, $len;
	my $motifFiles = "";
	foreach(@lens) {
		my $L = $_;
		my $OPTS = $opts . " -len $L";
		my $outFileLen = $outFile . ".length$L";
		$motifFiles .= " \"$outFileLen\"";
		`homer2 denovo -s "$tmpSeqFile" -g "$tmpGroupFile" $OPTS > "$outFileLen"`;
	}
	`cat $motifFiles > "$outFile.lengthAll"`;
	my $revopp = "";
	if ($strand eq '+' || $strand eq '-' || $strand eq 'separate') {
		$revopp = " -norevopp";
	}
	`compareMotifs.pl "$outFile.lengthAll" "$outFile.results/" -reduceThresh $reduceThresh -matchThresh $matchThresh -known $knownMotifsFile $bitsFlag -cpu $maxCPUs2 $rnaFlag $revopp`;
	`cp "$outFile.results/nonRedundant.motifs" "$outFile"`;
	`rm "$tmpSeqFile" "$tmpGroupFile"`;
}
sub calcSeqProperties {
	my ($seq, $order,$filename) = @_;

	my %properties = ();

	my @alpha = ('A','C','G','T');
	my ($kmerIndex,$kmerArray) = createKmerIndex($order,\@alpha);
	$properties{"_kmerIndex"} = $kmerIndex;
	$properties{"_kmerArray"} = $kmerArray;

	if (defined($filename)) {
		open OUT, ">$filename";
		print OUT "sequenceID";
		for (my $i=0;$i<@$kmerArray;$i++) {
			print OUT "\t$kmerArray->[$i]";
		}
		print OUT "\n";
	}

	foreach(keys %$seq) {
		my $name = $_;
		my $seq = $seq->{$name};
		my $res = getSeqProperties($seq,$order,$kmerIndex);
		$properties{$name} = $res;
		if (defined($filename)) {
			print OUT "$name";
			foreach(@$res) {
				print OUT "\t$_";
			}
			print OUT "\n";
		}
	}
	close OUT if (defined($filename));
	return \%properties;
}
sub getSeqProperties {
	my ($seq,$order,$kmerIndex) = @_;
	my @res = ();
	my @N = 0;
	for (my $j=0;$j<=$order;$j++) {
		my $f = 4**($j+1);
		for (my $i=0;$i<$f;$i++) {
			push(@res, 0);
		}
		push(@N,0);
	}
	for (my $i=0;$i<length($seq);$i++) {
		for (my $j=0;$j<=$order;$j++) {
			my $L = $j+1;
			next if ($i+$j > length($seq));
			my $s = substr($seq,$i,$L);
			next if ($L != length($s));
			if (exists($kmerIndex->{$s})) {
				my $index = $kmerIndex->{$s};
				$N[$j]++;
				$res[$index]++;
			}
		}
	}
	my $start = 0;
	for (my $j=0;$j<=$order;$j++) {
		my $f = 4**($j+1);
		next if ($N[$j] < 1);
		for (my $i=$start;$i<$start+$f;$i++) {
			$res[$i] /= $N[$j];
		}
		$start = $start+$f;
	}
	return \@res;
}


sub createKmerIndex {
	my ($order,$alpha) = @_;
	my $index = 0;
	my @kmerArray = ();
	my %kmerIndex = ();
	if ($order >= 0) {
		foreach(@$alpha) {
			my $s = $_;
			push(@kmerArray, $s);
			$kmerIndex{$s}=$index++;
		}
	}
	if ($order >= 1) {
		foreach(@$alpha) {
			my $s1 = $_;
			foreach(@$alpha) {
				my $s2 = $s1 . $_;
				push(@kmerArray, $s2);
				$kmerIndex{$s2}=$index++;
			}
		}
	}
	if ($order >= 2) {
		foreach(@$alpha) {
			my $s1 = $_;
			foreach(@$alpha) {
				my $s2 = $s1 . $_;
				foreach(@$alpha) {
					my $s3 = $s2 . $_;
					push(@kmerArray, $s3);
					$kmerIndex{$s3}=$index++;
				}
			}
		}
	}
	if ($order >= 3) {
		foreach(@$alpha) {
			my $s1 = $_;
			foreach(@$alpha) {
				my $s2 = $s1 . $_;
				foreach(@$alpha) {
					my $s3 = $s2 . $_;
					foreach(@$alpha) {
						my $s4 = $s3 . $_;
						push(@kmerArray, $s4);
						$kmerIndex{$s4}=$index++;
					}
				}
			}
		}
	}
	return (\%kmerIndex,\@kmerArray);
}
sub parseBestGuess {
	my ($name) = @_;
	$name =~ s/^.*BestGuess\://;
	$name =~ s/\/.*$//;
	return $name;
}
	

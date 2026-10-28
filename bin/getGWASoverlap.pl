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
use Statistics;
use HomerConfig;
use File::Spec;
use File::Basename;

my $maxBedNameLength = 200;
my $dist = 1000;
my $minSNPs =10;
my $maxCPUs = 1;
my $gregorPopulation = "EUR";
my $population = "EUR";
#my $population = "CEU";
my $LDthresh = 0.7;
my $seqOffset = -30;

sub printCMD {
	print STDERR "\n\tgetGWASoverlap.pl <gwas catolog file> -p <peak file1> [peak file2] ... [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-d <#> (Overlap distance, default: $dist)\n";
	print STDERR "\t\t-min <#> (minimum number of significant SNPs to consider study, default: $minSNPs)\n";
	print STDERR "\t\t-cpu <#> (number of threads to use, default: $maxCPUs)\n";
	#print STDERR "\t\t-empirical <#> (calculate FDR q-values empirically, <#> number of randomizations)\n";
	print STDERR "\t\t-GREGOR <path-to-GREGOR> <path-to-Reference LD info> (perform enrichment test with GREGOR)\n";
	print STDERR "\t\t\t-LD <#> (LD threshold for 'buddy SNPs', default 0.8, must be greater than 0.7)\n";
	print STDERR "\t\t-studyID <#> (only analyze this study)\n";
	print STDERR "\t\t-studyIDfile <#> (only analyze this study)\n";
	print STDERR "\t\t-snpOut <output file> (output overlapping snps & buddies as a BED file)\n";
	print STDERR "\t\t\t-snpOutAll (output all SNPs, even if they aren't overlapping)\n";
	#print STDERR "\t\tMotif scanning (all arguments are required):\n";
	#print STDERR "\t\t\t-m <motifs file> (Motifs to scan for)\n";
	#print STDERR "\t\t\t-genome <genome version> (genome version)\n";
	#print STDERR "\t\t\t-snps <parsed SNP definitions> (custom)\n";
	#print STDERR "\t-LDlink\n";
	#print STDERR "\t-localLDdb <snp.LD.table.tsv> (local proxy/buddy SNP table)\n";
	print STDERR "\n\tThe gwas catalog file can be downloaded from UCSC annotation database:\n";
	print STDERR "\t\ti.e.: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/gwasCatalog.txt.gz\n";
	print STDERR "\n";
	exit;
}



if (@ARGV < 1) {
	printCMD();
}


my $rand = rand();
my $gregorFlag = 0;
my $LDlinkFlag = 0;
my $localLDdbFile = '';

my $gwasFile = $ARGV[0];
my @peakFiles = ();
my $numEmpirical = 0;
my $pubmedIDflag = 0;
my $snpOutFile = '';
my $snpOutAllFlag = 0;
my $motifFile = '';
my $genome = '';
my $snpDefFile = '';
my %pubmedID = ();
%pop = ();
%buddy = ();
%buddyAllele = ();

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-d') {
		$dist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-m') {
		$motifFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-snps') {
		$snpDefFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$minSNPs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-LDlink') {
		$LDlinkFlag = 1;
	} elsif ($ARGV[$i] eq '-localLDdb') {
		$localLDdbFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-LD') {
		$LDthresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-snpOutAll') {
		$snpOutAllFlag = 1;
	} elsif ($ARGV[$i] eq '-snpOut') {
		$snpOutFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pubmedID' || $ARGV[$i] eq '-studyID') {
		$pubmedID{$ARGV[++$i]} = 1;
		$pubmedIDflag = 1;
	} elsif ($ARGV[$i] eq '-pubmedIDfile' || $ARGV[$i] eq '-studyIDfile') {
		open IN, $ARGV[++$i] or die "Couldn't open file specified after -pubmedIDfile/-studyIDfile\n";
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			#$line[0] =~ s/^.*\(//;
			#$line[0] =~ s/\).*$//;
			$pubmedID{$line[0]} = 1;
			print STDERR "Adding $line[0]\n";
		}
		close IN;
		$pubmedIDflag = 1;
	} elsif ($ARGV[$i] eq '-empirical') {
		$numEmpirical = $ARGV[++$i];
		print STDERR "!!! Error - '-empirical' is disabled in this version\n";
		exit;
	} elsif ($ARGV[$i] eq '-GREGOR' || $ARGV[$i] eq '-gregor') {
		$gregorFlag = 1;
		$gregorExe = $ARGV[++$i];
		$gregorRef = $ARGV[++$i];
		print STDERR "\tWill use GREGOR:\n";
		print STDERR "\t\tprogram location: $gregorExe\n";
		print STDERR "\t\tLD info location: $gregorRef\n";
	} elsif ($ARGV[$i] eq '-p') {
		$i++;
		while ($i<@ARGV && $ARGV[$i] !~ /^-/) {
			push(@peakFiles, $ARGV[$i++]);
		}
		$i--;
	} else {
		print STDERR "Couldn't recognize $ARGV[$i]\n";
		printCMD();
	}
}

if (0) {
my $ttmotifs = HomerConfig::readMotifFile($motifFile);
my $ttpos = -10;
my $ttseq = 'AGAACAGCCTGTTCTT';
my $ttstrand = "-";
my $ttref = 'G';
my $ttalt = 'A';
#my $ttseq = 'ARE(NR)/LNCAP-AR-ChIP-Seq(GSE27824)/Homer';
my $ttm = $ttmotifs->[0];
my $nnnn = $ttm->{'matrix'}->[0][0];
print STDERR "motif name = $nnnn   $ttref $ttalt\n";
my $ttdiff = scoreDifferenceSNP($ttpos,$ttseq,$ttstrand,$ttref,$ttalt,$ttm);
exit;
}




if ($gregorFlag) {
	$maxCPUs = ceil($maxCPUs/10.0);
	print STDERR "\tAdjusting simultaneous processes to $maxCPUs (GREGOR uses ~10 cpus per process)\n";
}
my %limits = ();

print STDERR "\n\tMinimum Risk SNPs per study: $minSNPs\n";
print STDERR "\tMaximum Distance between Peaks and risk SNPs: $dist\n\n";

print STDERR "\n\tAnalyzing peak files:\n";
my @fdrPvalues = ();

for (my $i=0;$i<@peakFiles;$i++) {
	print STDERR "\t\t$peakFiles[$i]\n";
	if ($numEmpirical > 0) {
		my $tmpFile = $rand . ".tmp";
		`bed2pos.pl "$peakFiles[$i]" -check > "$tmpFile"`;
		open IN, $tmpFile;
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			next if ($line[0] =~ /^#/);
			next if (@line < 5);
			if (!exists($limits{$line[1]})) {
				$limits{$line[1]} =0;
			}
			$limits{$line[1]} = $line[3] if ($limits{$line[1]} < $line[3]);
		}
		close IN;
		`rm $tmpFile`;
		my @a = ();
		push(@fdrPvalues, \@a);
	}
}
print STDERR "\n";

my @bedFiles = ();
my $gregorBedInputIndexFile = $rand . ".gregor.bed.index";
my $absBedIndexFile = "";
my $absGregorRef = "";
if ($gregorFlag) {
	open GREGINDEX, ">$gregorBedInputIndexFile";
	for (my $i=0;$i<@peakFiles;$i++) {
		my $tmpFile = $rand . ".tmp";
		`bed2pos.pl "$peakFiles[$i]" -check > "$tmpFile"`;
		my $bedFileName= $rand . "." . $peakFiles[$i] . ".bed";
		`pos2bed.pl "$tmpFile" | grep -v chrUn > "$bedFileName"`;
		push(@bedFiles, $bedFileName);
		my $absFile = File::Spec->rel2abs($bedFileName);
		print GREGINDEX "$absFile\n";
	}
	close GREGINDEX;
	$absBedIndexFile = File::Spec->rel2abs($gregorBedInputIndexFile);
	$absGregorRef = File::Spec->rel2abs($gregorRef);
}
my $nbed = scalar(@bedFiles);

#print out the header
print "Study";
print "\tTotal risk SNPs";
for (my $i=0;$i<@peakFiles;$i++) {
	print "\t$peakFiles[$i] Overlap";
	print "\t$peakFiles[$i] p-value";
	if ($numEmpirical > 0) {
		print "\t$peakFiles[$i] q-value FDR/Empirical(n=$numEmpirical)";
	} else {
		print "\t$peakFiles[$i] q-value FDR/Benjamini";
	}
	print "\t$peakFiles[$i] SNPs";
}
print "\n";

my $limitTotal = 0;
if ($numEmpirical > 0) {
	#print STDERR "\tChromosome Limits:\n";
	foreach(keys %limits) {
		#print STDERR "\t\t$_\t$limits{$_}\n";
		$limitTotal += $limits{$_};
	}
	#print STDERR "\tApprox Size: $limitTotal\n\n";
}

my $gwas = parseGWAScatalog($gwasFile, $pubmedIDflag, \%pubmedID);

if ($localLDdbFile ne '') {
	print STDERR "\tWill now look through local database for buddy SNPs\n";
	readLocalLDdbFile($localLDdbFile);
}


@results = ();
@studies = ();

my $numStudy = scalar(keys %$gwas);
my $count2 = 0;
my $cpus=0;
my %tmpIDs = ();
foreach(keys %$gwas) {
	$count2++;
	#last if ($count2 > 40);
	my $study = $_;
	my @snps = keys %{$gwas->{$study}};
	my $N = scalar(@snps);
	next if (@snps < $minSNPs);

	my $tmpID = $rand . "-$count2";
	my $pid = fork();
	if ($pid != 0) {
		#parent process
		#print STDERR "$pid\t$tmpID\n";
		$tmpIDs{$pid} = $tmpID;
		$cpus++;
		if ($cpus >= $maxCPUs) {
			my $id = wait();
			my $tid = $tmpIDs{$id};
			#print STDERR "\tProcessing tmpID: $tid\t process ID: $id\n";
			processResults($tid);
			$cpus--;
		}
		next;
	}


	print STDERR "\tAnalyzing Study $count2 of $numStudy ($study, $N total risk SNPs)\n";

	my @peakResults = ();

	my $resultsFile = $tmpID . ".results";
	open RESULTS, ">$resultsFile";

	my $numDsnps = $N;
	print RESULTS "$study\t$N\t$tmpID\n";
	#push(@studies, "$study\t$N");


	if ($gregorFlag) {
		my $snpIndexFile = "$tmpID" . ".snps.txt";
		open OUT, ">$snpIndexFile";
		foreach(@snps) {
			my $allele = $_;
			$allele =~ s/\-.+$//;
			#my $p = $gwas->{$study}->{$allele};
			print OUT "$allele\n";
		}
		close OUT;
		my $absSnpIndexFile = File::Spec->rel2abs($snpIndexFile);
		my $outputDir = $tmpID . "-outDir";
	   	$absOutputDir = File::Spec->rel2abs($outputDir);
		
		my $confFile = $tmpID . ".conf";
		open OUT, ">$confFile" or die "!!!Could not open $confFile for writing (GREGOR configuration file)\n";
		print OUT "#############################################################################\n";
		print OUT "# CHIPSEQ ENRICHMENT CONFIGURATION FILE\n";
		print OUT "# This configuration file contains run-time configuration of\n";
		print OUT "# CHIP_SEQ ENRICHMENT\n";
		print OUT "#############################################################################\n";
		print OUT "# KEY ELEMENTS TO CONFIGURE : NEED TO MODIFY\n";
		print OUT "#############################################################################\n";
		print OUT "INDEX_SNP_FILE = $absSnpIndexFile\n";
		print OUT "BED_FILE_INDEX = $absBedIndexFile\n";
		print OUT "REF_DIR = $absGregorRef\n";
		print OUT "R2THRESHOLD = $LDthresh ## must be greater than 0.7\n";
		print OUT "LDWINDOWSIZE = 1000000 ## must be less than 1MB; these two values define LD buddies\n";
		print OUT "OUT_DIR = $absOutputDir\n";
		print OUT "MIN_NEIGHBOR_NUM = 500 ## define the size of neighborhood\n";
		print OUT "BEDFILE_IS_SORTED = false  ## false, if the bed files are not sorted\n";
		print OUT "POPULATION = $gregorPopulation  ## define the population, you can specify EUR, AFR, AMR or ASN\n";
		print OUT "TOPNBEDFILES = 2000\n";
		print OUT "JOBNUMBER = 10\n";
		print OUT "#############################################################################\n";
		print OUT "#############################################################################\n";
		print OUT "BATCHTYPE = local ##  run jobs on local machine\n";
		close OUT;

		print STDERR "\tRunning GREGOR\n";
		print STDERR "perl $gregorExe --conf $confFile\n";
		`perl $gregorExe --conf $confFile`;
		print STDERR "\tGREGOR finished\n";

		if (-e "$outputDir/StatisticSummaryFile.txt") {
			open IN, "$outputDir/StatisticSummaryFile.txt";
			my $c = 0;
			my %peaks = ();
			while (<IN>) {
				$c++;
				next if ($c==1);
				chomp;
				s/\r//g;
				my @line = split /\t/;
				my $bedname = "$line[0]";
				my $oexp = "$line[1] ($line[2])";
				my $pvalue = $line[3];
				$peaks{$bedname} = {oexp=>$oexp,p=>$pvalue};
			}
			close IN;

			`ls -d $outputDir/index.SNP.and.LD.for.top.*/ > $outputDir/.dir`;
			open IN, "$outputDir/.dir";
			my $dd = $outputDir . "/index.SNP.and.LD.for.top." . $nbed . ".bed/";
			while (<IN>) {
				chomp;
				$dd = $_;
				last;
			}
			close IN;
	
			foreach(@bedFiles) {
				my $bedname = $_;
				my $N = 0;
				my $pvalue = 1;
				my $snpStr = '';
				my $fdr = 1;
				my $pvalueStr = '';
					if (exists($peaks{$bedname})) {
					$N = $peaks{$bedname}->{'oexp'};
					$pvalue = $peaks{$bedname}->{'p'};
				}
				my $overlapSNPsFile = $dd . "index.SNP.and.LD.in.bed." . $bedname . ".txt";
				if (-e $overlapSNPsFile) {
					print STDERR "\t\tCool - found overlap file.\n";
					open IN, $overlapSNPsFile;
					while (<IN>) {
						chomp;
						s/\r//g;
						my @line = split /\t/;
						next if (@line < 4);
						next if ($line[1] eq 'LD');
						if ($line[2] eq '-1') {
							# not overlapping a peak
							#$snpStr .= "chr" . $line[1] . "-n" . ",";
						} else {
							# overlapping a peak
							$snpStr .= "chr" . $line[1] . "-o" . ",";
						}
					}
					close IN;
				} else {
					print STDERR "Uh... $overlapSNPsFile is not right. Possible:\n";
					print STDERR "Listing of output Directory:\n";
					`ls $outputDir/`;
					`ls $outputDir/index.SNP.and.LD.for.top*`;
				}
	
				print RESULTS "$bedname\t$N\t$pvalue\t$snpStr\t$fdr\t$pvalueStr\n";
			}
			close RESULTS;
		} else {
			print STDERR "\t\tproblem, maybe no overlaping SNPs... ($study, $N total risk SNPs)\n";
		}
		#print STDERR "sleeping...\n";
		#`sleep 10000`;

		`rm "$snpIndexFile"`;
		`rm -r "$outputDir"` if ($outputDir ne '' || $outputDir ne '/');

		exit(0);
	}




	my @ldsnps = ();

	if ($LDlinkFlag) {
		my $zz = scalar(@snps);
		my $c=0;
		my $tmpFile = $tmpID . ".tmp";
		foreach(@snps) {
			$c++;
			$snp = $_;
			$snp =~ s/\-.*//;
			print STDERR "\t$c of $zz\t$snp\n";
			my $pop = $population;
			if (exists($pop{$study})) {
				$pop = $pop{$study};
			}
			`curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=$snp&pop=$pop&r2_d=r2' > $tmpFile`;
			open IN, $tmpFile;
			my $cc = 0;
			while (<IN>) {
				$cc++;
				next if ($cc < 2);
				chomp;
				s/\r//g;
				my @line = split /\t/;
				next if ($line[6] < $LDthresh);
				my $alt = $line[2];
				$alt =~ s/^.*\///;
				$alt =~ s/\)$//;
				my $allele = $line[0] . '-' . $alt;
				my @a = split /\:/,$line[1];
				my $c = $a[0];
				my $s = $a[1];
				my $e = $a[1];
				my $position = "$allele\t$c\t$s\t$e\t0";
				$gwas->{$study}->{$allele} = $position;
			}
			close IN;
		}
		`rm $tmpFile`;
		my $oldN = $N;
		@snps = keys %{$gwas->{$study}};
		$N = scalar(@snps);
		print STDERR "\t$study - $oldN -> $N with LD snps\n";
	} elsif ($localLDdbFile ne '') {
		#print STDERR "snps=@snps\n";
		foreach(@snps) {
			$allele = $_;
			my $snp = $allele;
			$snp =~ s/\-.*//;
			#print STDERR "Here??? $allele ($snp)\n";

			my $pop = $population;
			if (exists($pop{$study})) {
				$pop = $pop{$study};
			}

			if (!exists($buddy{$snp})) {
				print STDERR "!!! Warning: $snp not found in localLDdbFile\n";
				next;
			}
			foreach(keys %{$buddy{$snp}}) {
				my $allele2 = $_;
				my $position = $buddy{$snp}->{$allele2};
				#print STDERR "$allele\t->\t$allele2\t$position\n";
				$gwas->{$study}->{$allele2} = $position;
			}
		}
		my $oldN = $N;
		@snps = keys %{$gwas->{$study}};
		$N = scalar(@snps);
		print STDERR "\t$study - $oldN -> $N with LD snps\n";
		#@ldsnps = @snps;
	}

	# HOMER simplified enrichment scheme (instead of GREGOR)
	my $tmpFile = $tmpID . ".tmp";
	open OUT, ">$tmpFile";
	foreach(@snps) {
		my $allele = $_;
		my $p = $gwas->{$study}->{$allele};
		print OUT "$p\n";
	}
	close OUT;
	my $tmpFile2 = $tmpID . ".2.tmp";

	for (my $i=0;$i<@peakFiles;$i++) {
		`cp "$peakFiles[$i]" "$tmpFile2"`;

		`mergePeaks "$tmpFile" "$tmpFile2" -matrix $tmpID  -d $dist > /dev/null 2> /dev/null`;

		my $logp = 0;

		my $c = 0;
		open IN, "$tmpID.logPvalue.matrix.txt";
		while (<IN>) {
			$c++;
			next if ($c < 2);
			s/\r//g;
			my @line = split /\t/;
			$logp = $line[2];
			last;
		}
		close IN;
		my $pvalue = 1;
		if ($logp < 0.1) {
			$pvalue = exp($logp);
		} else {
			$pvalue = 1-exp(-1*$logp);
		}
		`rm $tmpID.*matrix.txt`;

		`mergePeaks "$tmpFile" "$tmpFile2" -cobound 1 -prefix $tmpID -d $dist > /dev/null 2> /dev/null`;

		my @bound = ();
		my $snpStr = "";
		open IN, "$tmpID.coBoundBy1.txt";	
		while (<IN>) {
			chomp;
			next if (/^#/);
			my @line = split /\t/;
			push(@bound, $line[0]);
			$snpStr .= "$line[1]:$line[2]-o" . ",";
			#$snpStr .= "$line[0]-o" . ",";
		}
		close IN;
		my $NN = scalar(@bound);	
		`rm $tmpID.coBoundBy*`;

		my $fdr = 1;
		my $pvalueStr = '';

		if ($numEmpirical > 0) {
			my $fileStr = "";
			my @chrs = keys %limits;
			for (my $j=0;$j<$numEmpirical;$j++) {
				my @rands = ();
				for (my $k=0;$k<$numDsnps;$k++) {
					push(@rands, floor($limitTotal*rand()));
					#print STDERR "$j $k $rands[$k]\n";
				}
				@rands = sort {$a <=> $b} @rands;
				my $f = "$tmpID.empirical.$j.tmp";
				$fileStr .= " $f";
				open OUT, ">$f";
				my $chrIndex = 0;
				my $totalIndex = 0;
				foreach(@rands) {
					my $r = $_;	
					while ($totalIndex+$limits{$chrs[$chrIndex]}< $r) {
						$totalIndex += $limits{$chrs[$chrIndex++]};
					}
					my $p = $r-$totalIndex;
					my $p2 = $p+1;
					print OUT "$r\t$chrs[$chrIndex]\t$p\t$p2\t+\n";
				}
				close OUT;
			}
			#print STDERR "`mergePeaks $tmpFile2 -cobound 1 -matrix $tmpID -d $dist $fileStr > /dev/null`;\n";
			`mergePeaks "$tmpFile2" -cobound 1 -matrix $tmpID -d $dist $fileStr -gsize $limitTotal > /dev/null 2> /dev/null`;

			my @pvalues = ();
			open IN, "$tmpID.logPvalue.txt";
			my $c = 0;
			while (<IN>) {
				$c++;
				next if ($c < 2);
				chomp;
				my @line = split /\t/;
				my $p = 1;
				if ($line[1] <= 0.1) {
					$p = exp($line[1]);
				} else {
					$p = 1-exp(-1*$line[1]);
				}
				push(@pvalues, $p);
				$pvalueStr .= "$p|";
				#print STDERR "\t$pvalue\t$p\n";
			}
			close IN;
			$pvalueStr .= ",";
			`rm -f $fileStr`;
			push(@{$fdrPvalues[$i]},\@pvalues);

		}

		#my $d = {n=>$N,p=>$pvalue,snps=>$snpStr,f=>$fdr};	
		#push(@peakResults, $d);
		print RESULTS "$peakFiles[$i]\t$NN\t$pvalue\t$snpStr\t$fdr\t$pvalueStr\n";
	}
	#push(@results, \@peakResults);
	`rm "$tmpFile" "$tmpFile2"`;

	close RESULTS;
	exit(0);
}

%snpDef = ();
if ($snpDefFile ne '') {
	print STDERR "\tReading dbSNP file\n";
	open IN, $snpDefFile;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $id = $line[1] . ":" . $line[2];
		my @a = ($line[5],$line[6]);
		$snpDef{$id} = \@a;
	}
	close IN;
}


$id = 0;
while ($id >= 0) {
	$id = wait();
	last if ($id < 0);
	my $tid = $tmpIDs{$id};
	#print STDERR "\tProcessing tmpID: $tid\t process ID: $id\n";
	processResults($tid);
}

for (my $i=0;$i<@peakFiles;$i++) {
	my @pvalues = ();
	foreach(@results) {
		push(@pvalues, $_->[$i]->{'p'});
	}
	my $fdrs = "";
	if ($numEmpirical < 1) {
		$fdrs = Statistics::benjaminiFDR(\@pvalues);
	} else {
		my $aa = scalar(@pvalues);
		my $bb = scalar(@{$fdrPvalues[$i]});
		$fdrs = Statistics::empiricalFDR(\@pvalues, $fdrPvalues[$i]);
	}
	for (my $j=0;$j<@results;$j++) {
		if (!exists($results[$j]->[$i]->{'p'})) {
			$results[$j]->[$i]->{'p'} = 'NA';
		}
		if ($results[$j]->[$i]->{'p'} eq 'NA') {
			$results[$j]->[$i]->{'f'} = 'NA';
		} else {
			$results[$j]->[$i]->{'f'} = $fdrs->[$j];
		}
	}
}
my %mutations = ();
if ($motifFile ne '' && $snpDefFile ne '' && $genome ne '') {
	my $tmpFile = $rand . ".in.tmp";
	my $tmpFile2 = $rand . ".in2.tmp";
	my $couldntfind = 0;
	open OUT, ">$tmpFile";
	for (my $i=0;$i<@results;$i++) {
		#print $studies[$i];
		my $z = 0;
		foreach(@{$results[$i]}) {
			#print "\t$_->{'n'}\t$_->{'p'}\t$_->{'f'}\t$_->{'snps'}";
			if ($snpOutFile ne '' && exists($_->{'snps'})) {
				my @snps = split /\,/,$_->{'snps'};
				#print STDERR "$_->{'snps'} = @snps\n";
				my $sname = $studies[$i];
				my $expName = $peakFiles[$z];
				$sname =~ s/\s/_/g;
				foreach(@snps) {
					my $snpCode =$_;
print STDERR "snpCode=$snpCode\n";
					$snpCode =~ /^(.*?)\:(\d+)\-(.)$/;
					my $chr = $1;
					my $pos = $2;
					my $overlap = $3;
					my $id = $chr . ":" . $pos;
					my $ref = '';
					my $alt = '';
					if (exists($snpDef{$id})) {
						$ref = $snpDef{$id}->[0];
						$alt = $snpDef{$id}->[1];
					} else {
						$couldntfind++;
						next;
					}
					my $name = $snpCode . "-" . $ref . "->" . $alt;
					$mutations{$snpCode} = "$ref>$alt";

					if ($overlap eq 'o') {
						# this is one to look at!
						my $s = $pos+$seqOffset;
						my $e = $pos-$seqOffset;
						print OUT "$name\t$chr\t$s\t$e\t+\n";
					} else {
						next;
					}
				}
			}
			$z++;
		}
	}
	close OUT;
	`annotatePeaks.pl $tmpFile $genome -nogene -noann -size 50 -m $motifFile > $tmpFile2`;
	my $motifs = HomerConfig::readMotifFile($motifFile);
	my %motifs = ();
	foreach(@$motifs) {
		$motifs{$_->{'name'}} = $_;
	}
	open IN, $tmpFile2;
	my $c=0;
	my @mcols = ();
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($c == 1) {
			for (my $i=0;$i<@line;$i++) {
				my $mname = $line[$i];
				$mname =~ s/ Distance.*$//;
				push(@mcols, $mname);
			}
			next;
		}
		my $pid =$line[0];
		$pid =~ /^(.*)\-(.*?)\-\>(.*)/;
		my $snpID = $1;
		my $ref=$2;
		my $alt=$3;
		for (my $i=9;$i<@line;$i++) {
			next if ($line[$i] eq '');
			if ($line[$i] =~ /^(.*?)\((.*?),([\+\-]),/) {
				my $pos = $1;
				my $seq = $2;
				my $strand = $3;
				if ($pos <= 0 && $pos+length($seq) >= 1) {
					#overlap!!
					my $diff = scoreDifferenceSNP($pos,$seq,$strand,$ref,$alt,$motifs{$mcols[$i]});
					my $msite = $mcols[$i] . "-$pos-$seq-$strand-delta:$diff";
					#print STDERR "msite=$msite\tsnpID=$snpID\tpid=$pid\n";
					if (!exists($snpMotifs{$snpID})) {
						my %m = ();
						$snpMotifs{$snpID} = \%m;
					}
					$snpMotifs{$snpID}->{$msite}=1;
				}
			}
		}
	}
	close IN;
	#`rm $tmpFile $tmpFile2`;	
}

if ($snpOutFile ne '') {
	open SNPBED, ">$snpOutFile";
	print SNPBED "track name=\"GWAS overlap SNPs\" itemRgb=\"On\"\n";
}
for (my $i=0;$i<@results;$i++) {
	print $studies[$i];
	my $z = 0;
	foreach(@{$results[$i]}) {
		print "\t$_->{'n'}\t$_->{'p'}\t$_->{'f'}\t$_->{'snps'}";
		if ($snpOutFile ne '' && exists($_->{'snps'})) {
			my @snps = split /\,/,$_->{'snps'};
			#print STDERR "$_->{'snps'} = @snps\n";
			my $sname = $studies[$i];
			my $expName = $peakFiles[$z];
			$sname =~ s/\s/_/g;
			foreach(@snps) {
				my $snpCode =$_;
				$snpCode =~ /^(.*?)\:(\d+)\-(.)$/;
				my $chr = $1;
				my $pos = $2;
				my $overlap = $3;
				$color = "180,180,180";
				if ($overlap eq 'o') {
					$color = "5,5,155";
					$overlap = "overlap";
				} elsif ($overlap eq 'n') {
					$overlap = "no-overlap";
				}
				if (exists($snpMotifs{$snpCode})) {
					$color = "0,100,0";
					my %dbd = ();
					foreach(keys %{$snpMotifs{$snpCode}}) {
						my $x = $_;
						$x =~ s/\).*?\/Homer\-/\)\-/;
						$x =~ s/\-.+\-delta/\-delta/;
						if ($x =~ /^.+\((.*?)\)/) {
							my $d = $1;
							if (exists($dbd{$d}) && $dbd{$d} > 3) {
								next;
							}
							$dbd{$d}++;
						}

						$overlap .= "|$x";
					}
				}
				my $mut = 'NA';
				if (exists($mutations{$snpCode})) {
					$mut = $mutations{$snpCode};
				}
				my $bname = $sname . "|" . $expName . "|$chr:$pos|" . "$mut|" . $overlap;
				my $start = $pos-1;
				if (length($bname) > $maxBedNameLength) {
					$bname = substr($bname,0,$maxBedNameLength);
				}
				print SNPBED "$chr\t$start\t$pos\t$bname\t0\t+\t$start\t$pos\t$color\n";
			}
		}
		$z++;
	}
	print "\n";
}	
if ($snpOutFile ne '') {
	close SNPBED;
}


if ($gregorFlag) {
	foreach(@bedFiles) {
		`rm "$_"`;
	}
}

`rm -rf $rand*`;
exit;

sub processResults {
	my ($tid) = @_;
	open IN, "$tid.results" or print STDERR "!!! Could not open results for $tid\n";
	my $c = 0;
	my $study = '';
	my $N = 0;
	my @r = ();
	my $peakIndex = 0;
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($c == 1) {
			$study = $line[0];
			$N = $line[1];
			next;
		}
		#print RESULTS "$peakFiles[$i]\t$N\t$pvalue\t$snpStr\t$fdr\t$pvalueStr\n";
		my $n = $line[1];
		my $pvalue = $line[2];
		my $snpStr = $line[3];
		my $fdr = $line[4];
		my $d = {n=>$n,p=>$pvalue,snps=>$snpStr,f=>$fdr};	
		push(@r, $d);

		my $epi = '';
		if (@line > 5) {
			$epi = $line[5];
			my @randomizations = split /\,/, $epi;
			if (@randomizations > 1) {
				pop @randomizations;
			}
			foreach(@randomizations) {
				my @ps = split /\|/, $_;
				pop @ps;
				push(@{$fdrPvalues[$peakIndex]},\@ps);
			}
		}
		$peakIndex++;

	}
	close IN;
	push(@studies, "$study\t$N");		
	push(@results, \@r);
	`rm "$tid.results"`;
}



sub parseGWAScatalog {
	my ($file, $pubmedIDflag, $pubmedID) = @_;
	my %gwas = ();
	my $format = 'gwas';


	my @givenIDs = keys %{$pubmedID};
	if (scalar(@givenIDs) > 0) {
		print STDERR "\n\tStudies: @givenIDs\n";	
	}


	open IN, $file;
	my $c = 0;
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($c == 1) {
			if (@line < 20) {
				print STDERR "\tAssuming SNP file foramt\n";
				$format = 'snps';
			} else {
				print STDERR "\tAssuming file is UCSC GWAS Catalog format...\n";
				$format = 'gwas';
			}
		}

		if ($format eq 'gwas') {
			my $chr = $line[1];
			my $start = $line[2];
			my $end = $line[3];
			my $snp = $line[4];
			my $pubmed = $line[5];
			my $study = $line[10] . "($pubmed)";
			if ($pubmedIDflag) {
				next if (!exists($pubmedID->{$study}));
				#print STDERR "\tFound $study in GWAS catalog file\n";
			} else {
			}
			my $allele = $line[15];
			my $pvalue = $line[17];
			my $pop = 'EUR';
			my $eur = 0;
			my $afr = 0;
			my $eas = 0;
			my $amr = 0;
			$eur = 1 if ($line[11] =~ /European/i || $line[11] =~ /Finnish/i);
			$afr = 1 if ($line[11] =~ /African/i);
			$eas = 1 if ($line[11] =~ /Asian/i || $line[11] =~ /Chinese/i || $line[11] =~ /Korea/i);
			$amr = 1 if ($line[11] =~ /Hispanic/i || $line[11] =~ /Latin/i);
			if ($eur + $afr + $eas + $amr > 1) {
				$pop = 'ALL';
			} elsif ($eas) {
				$pop = 'EAS';
			} elsif ($afr) {
				$pop = 'AFR';
			} elsif ($amr) {
				$pop = 'AMR';
			} else {
				$pop = 'EUR';
			}
			#print STDERR "\t$study\t$pop\t$line[11]\n";


			my $position = "$allele\t$chr\t$start\t$end\t0";
		
			if (!exists($gwas{$study})) {
				my %snps = ();
				$gwas{$study} = \%snps;
			}
			if (!exists($buddy{$snp})) {
				my %a = ();
				$buddy{$snp} = \%a;
			}
			if (!exists($buddyAllele{$snp})) {
				my %a = ();
				$buddyAllele{$snp} = $allele;
			}
			$gwas{$study}->{$allele} = $position;
			$pop{$study} = $pop;
			#print STDERR "$position\n";
		} elsif ($format eq 'snps') {
			my $study = $file;
			my $allele = $line[0];
			my $position = "$allele" . "\t" . $_;
			if (!exists($gwas{$study})) {
				my %snps = ();
				$gwas{$study} = \%snps;
			}
			$gwas{$study}->{$allele} = $position;
		}
	}
	close IN;
	return \%gwas;
}
sub readLocalLDdbFile {
	my ($file) = @_;
	#access %buddy as a global variable.
	my $numFound = 0;
	my $numBuddies = 0;
	my $notFound = 0;
	open IN, $file;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $snp = $line[0];
		my $allele = $line[0] . "-" . $line[3];
		if (!exists($buddy{$snp})) {
			$notFound++;
			next;
		}

		#example line: rs146153292     chr10:122963036 chr10:122958201:rs117707949:1,chr10:123225644:rs554623169:1
		my @a = split /\,/,$line[4];
		my $buddies = 0;
		foreach(@a) {
			my @b = split /\:/;
			my $s = $b[2] . "-" . $b[3];
			$buddyAllele{$b[2]}=$s;
			next if ($b[4] < $LDthresh);
			my $position = "$s\t$b[0]\t$b[1]\t$b[1]\t0";
			$buddy{$snp}->{$s} = $position;
			$buddies++;
		}
		$numBuddies += $buddies;
		$numFound++;
	}
	close IN;
	my $dem = $numFound;
	$dem = 1 if ($numFound < 1);
	my $avg = $numBuddies / $dem;
	my $total = $numFound + $notFound;
	print STDERR "\t$numBuddies buddies (in LD) found for $numFound (of $total) tag SNPs ($avg per)\n";
}

sub scoreDifferenceSNP {
	my ($pos,$seq,$strand,$ref,$alt,$motif) = @_;
	
	#print STDERR "MName= $motif->{'name'}  $seq  $pos  $strand $ref $alt\n";
	if ($ref eq '-' || $alt eq '-' || length($ref) != length($alt)) {
		#indel - assume the worst for now
		return -10;
	}

	my $diff = 0;

	my $p = 0;
	my $R = $ref;
	my $A = $alt;
	my $S = '';
	$p = abs($pos);
	$S = substr($seq,$p,length($ref));
	my @refIndex=();
	my @altIndex=();
	if ($strand eq '+') {
	} elsif ($strand eq '-') {
		$p = length($seq)+$pos-length($ref);
		$S = HomerConfig::revopp($S);
		$R = HomerConfig::revopp($R);
		$A = HomerConfig::revopp($A);
	}
	my $sI = seq2index($S);
	my $aI = seq2index($A);
	for (my $i=0;$i<scalar(@$sI);$i++) {
		my $offset = $p+$i;
		if ($sI->[$i] > 3 || $aI->[$i] > 3) {
			return 0;
		}
		my $ps = $motif->{'matrix'}->[$offset][$sI->[$i]];
		my $pa = $motif->{'matrix'}->[$offset][$aI->[$i]];
		my $lps = log($ps/0.25);
		my $lpa = log($pa/0.25);
		my $v = $lpa-$lps;
		#print STDERR "$i ($offset)\t$v\t$sI->[$i] ($ps)\t$aI->[$i] ($pa)\n";
		$diff += $v;
	}
	#print STDERR "diff=$diff\n";
	$diff = sprintf("%.2lf", $diff);
	return $diff;
}

sub seq2index {
	my ($seq) = @_;
	my @a = ();
	for (my $i=0;$i<length($seq);$i++) {
		my $x = substr($seq,$i,1);
		if ($x eq 'A' || $x eq 'a') {
			push(@a, 0);
		} elsif ($x eq 'C' || $x eq 'c') {
			push(@a, 1);
		} elsif ($x eq 'G' || $x eq 'g') {
			push(@a, 2);
		} elsif ($x eq 'T' || $x eq 't') {
			push(@a, 3);
		} else {
			push(@a, 4);
		}
	}
	return \@a;
}

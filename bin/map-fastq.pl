#!/usr/bin/env perl
use warnings;
use File::Basename;

#notes:
# hisat2 low quality: --mp 1,1
my $outSAMmultNmaxDef = 1;
my $outFilterMultimapNmaxDef = 10000;
#$outSAMmultNmaxDef = 100;
#$outFilterMultimapNmaxDef = 100;

my $sendErr2devnull = " 2> /dev/null ";
my $prefix = '';

sub printCMD {
	print STDERR "\n\tUsage: map-fastq.pl [options] <FASTQ file1> [FASTQ file2]...\n";
	print STDERR "\t\tpaired end: <FASTQ end1>,<FASTQ end2> [FASTQ file2 end1],[end2]...\n";
	print STDERR "\t\t\tor just specify first file and add \"-pe\" flag to auto detect second file\n";
	print STDERR "\n\tRequired Options:\n";
	print STDERR "\t\t-p <#> (Number of cpus per instance, default: 1)\n";
	print STDERR "\t\t-x <path-to-bt2 index> (path to index to use for mapping)\n";
	print STDERR "\t\t-pe (if only 1st read file is specified, this will look for 2nd read and do PE alignment)\n";
	print STDERR "\t\t-o <output folder, default: same location as alignment files>\n";
	print STDERR "\t\t-suffix <suffix to apply to mapped files, default genome/program>\n";
	print STDERR "\t\t-fastp (autotrim with fastp - requires fastp to be available in PATH)\n";
	print STDERR "\n\tAlignment Program:\n";
	print STDERR "\t\t-bowtie2 (map with bowtie2)\n";
	print STDERR "\t\t-bowtie (map with bowtie, add \"-C\" if using colorspace)\n";
	print STDERR "\t\t-bwa (map with bwa mem)\n";
	print STDERR "\t\t-hisat2 (map with hisat2)\n";
	print STDERR "\t\t-star (map with STAR)\n";
#	print STDERR "\t\t-tophat2 (map with tophat2)\n";
	print STDERR "\t\t\t-path <path-to-program-file> (executable file to run if not in path/diff name)\n";
#	print STDERR "\n\tCPU options:\n";
#	print STDERR "\t\t-cpu <#> (Number of instances to run at once, default:1)\n";
	print STDERR "\tGlobal Options:\n";
	print STDERR "\t\t-un (will output unaligned reads)\n";
	print STDERR "\t\t-bam (convert output files from sam to bam)\n";
	print STDERR "\n\tBowtie options:\n";
	print STDERR "\t\t--local (local alignment, default: global/end-to-end)\n";
	print STDERR "\n\tSTAR options:\n";
    print STDERR "\t\t-NoSharedMemory (Do not used shared memory for genome [slower])\n";
    print STDERR "\t\t-2pass <readlength> (Perform 2-pass alignment using --twopass1readsN option)\n";
	print STDERR "\t\t--outSAMmultNmax <#> (Sets STAR parameter for number of alignments in output, def: $outSAMmultNmaxDef)\n";
	print STDERR "\t\t--outFilterMultimapNmax <#> (Sets STAR parameter for multimapper filtering, def: $outFilterMultimapNmaxDef)\n";
	#print STDERR "\n\tTophat options:\n";
	#print STDERR "\t\t--library-type <type> (library type for tophat, default: fr-firststrand\n";
	#print STDERR "\t\t\t\tOther optoins: fr-unstranded, fr-secondstrand)\n";
	#print STDERR "\t\t-G <GTF file> (Use as guide for splice mapping)\n";
	#print STDERR "\t\t-mis <#> (Max number of mismatches, default: 2)\n";
	print STDERR "\n\tGeneral Options to pass along to alignment program:\n";
	print STDERR "\t\t-f (Input is FASTA files, default expects FASTQ)\n";
	print STDERR "\t\t-pass \"...\" (need to include quotes)\n";
	print STDERR "\t\t-showErrors (don't dump subprogram stderr to /dev/null)\n";
	#print STDERR "\t\t-remap (remap unaligned reads with bowtie afterwards, returning random position)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 2) {
	printCMD();
}
my $bamSortRam = 5000000000;
#my $bamSortRam =  50241897905;
my $bt2Index = "";
my $maxCPUs = 1;
my $pCPUs = 1;
my $remapFlag = 0;
my $fastpFlag = 0;
my @files = ();
my $libraryType = "RF";
my $maxMultihits = 20;
my $outBAMsortingThreads = 6;
my $guideGTF = "";
my $program = '';
my $exe = '';
my $local = '';
my $pass = "";
my $unFlag = 0;
my $peFlag = 0;
my $n = 3;
my $pass2Flag = 0;
my $readLength = 0;
my $removeSharedMemory = 1;
my $genomeLoad = "LoadAndKeep";
my $starOptions = "--outSAMstrandField intronMotif --outMultimapperOrder Random --limitOutSAMoneReadBytes 10000000 --outBAMsortingThreadN $outBAMsortingThreads";
my $hisat2Options = " --dta ";
my $maxMisMatches = " -n $n --genome-read-mismatches $n --read-mismatches $n ";
my $bamFlag = 0;
my $suffix = '';
my $colorSpaceFlag = 1;

for (my $i=0;$i<@ARGV;$i++) {
	#print STDERR "$ARGV[$i] $i\n";
	if ($ARGV[$i] eq '-index' || $ARGV[$i] eq '-x') {
		$bt2Index = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bam') {
		$bamFlag= 1;
	} elsif ($ARGV[$i] eq '-pe') {
		$peFlag = 1;
	} elsif ($ARGV[$i] eq '-2pass') {
		$pass2Flag = $ARGV[++$i];
		print STDERR "\tWill run STAR with 2-pass mode with read length of $pass2Flag\n";
		$pass2Flag--;
		print STDERR "\t\tRunning STAR with --sjdbOverhang $pass2Flag and NoSharedMemory\n";
		$starOptions .= " --sjdbOverhang $pass2Flag ";
        $genomeLoad = "NoSharedMemory";	
    } elsif ($ARGV[$i] eq '-NoSharedMemory') {
        $genomeLoad = "NoSharedMemory";	
	} elsif ($ARGV[$i] eq '-mis') {
		my $n = $ARGV[++$i];
		$maxMisMatches = " -n $n --genome-read-mismatches $n --read-mismatches $n ";
	} elsif ($ARGV[$i] eq '-prefix' || $ARGV[$i] eq '-o') {
		$prefix = $ARGV[++$i];
		print STDERR "\tOutput files will be placed in the directory: $prefix\n";
	} elsif ($ARGV[$i] eq '-bowtie2') {
		$program = 'bowtie2';
	} elsif ($ARGV[$i] eq '-bowtie') {
		$program = 'bowtie';
	} elsif ($ARGV[$i] eq '-fastp') {
		$fastpFlag = 1;
		print STDERR "\tWill trim files using fastp with default parameters\n";
	} elsif ($ARGV[$i] eq '-C') {
		$colorSpaceFlag = 1;
		$pass .= " -C ";
	} elsif ($ARGV[$i] eq '--outSAMmultNmax') {
		$outSAMmultNmaxDef = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '--outFilterMultimapNmax') {
		$outFilterMultimapNmaxDef = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-tophat2') {
		$program = 'tophat2';
	} elsif ($ARGV[$i] eq '-suffix') {
		$suffix = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bwa') {
		$program = 'bwa';
	} elsif ($ARGV[$i] eq '-star') {
		$program = 'STAR';
	} elsif ($ARGV[$i] eq '-hisat2') {
		$program = 'hisat2';
	} elsif ($ARGV[$i] eq '-f') {
		$pass .= " -f ";
	} elsif ($ARGV[$i] eq '--local') {
		$local = '--local';
	} elsif ($ARGV[$i] eq '-path') {
		$exe = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-un') {
		$unFlag = 1;
	} elsif ($ARGV[$i] eq '-p') {
		$pCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-remap') {
		$remapFlag = 1;
	} elsif ($ARGV[$i] eq '-showErrors') {
		$sendErr2devnull = "";
	} elsif ($ARGV[$i] eq '--library-type') {
		$libraryType = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pass') {
		my $stuff = $ARGV[++$i];
		$stuff=~ s/^\"//;
		$stuff=~ s/\"$//;
		$pass .= " " . $stuff
	} elsif ($ARGV[$i] eq '-G') {
		$guideGTF = " -G $ARGV[++$i]";
	} elsif ($ARGV[$i] =~ /^\-/) {
		print STDERR "\n!!! \"$ARGV[$i]\" not recognized...\n";
		printCMD();
	} else {
		push(@files, $ARGV[$i]);
	}
}


$hisat2Options .= " --rna-strandness $libraryType ";
$starOptions .= " --outSAMmultNmax $outSAMmultNmaxDef  --outFilterMultimapNmax $outFilterMultimapNmaxDef ";

if ($program eq '') {
	if ($bt2Index =~ /STAR/) {
		$program = 'STAR';
		print STDERR "\tAutodetecting program as STAR given the index\n";
	} elsif ($bt2Index =~ /bowtie2/) {
		$program = 'bowtie2';
		print STDERR "\tAutodetecting program as bowtie2 given the index\n";
	} elsif ($bt2Index =~ /bowtie/) {
		$program = 'bowtie';
		print STDERR "\tAutodetecting program as bowtie given the index\n";
	} elsif ($bt2Index =~ /bwa/) {
		$program = 'bwa';
		print STDERR "\tAutodetecting program as bwa given the index\n";
	} else {
		print STDERR "!!! Error: Must specify a program (-bowtie2, -star, -hisat2, or -bwa)\n";
		exit;
	}
}

#$program .= " -f ";
if ($exe eq '') {
	$exe = $program;
	print STDERR "\tWill run $exe ...\n";
}

if ($bt2Index eq '') {
	print STDERR "!!! Need to specify mapping index and/or CPUs !!!\n";
	printCMD();	
}

print STDERR "\tAlignment Program:\n\t";
my $a = `which $exe`;
print STDERR "$a\n";

my @file1 = ();
my @file2 = ();
foreach(@files) {
	my $f1 = $_;
	my $f2 = '';
	if ($_ =~ /\,/) {
		my @f = split /\,/;
		$f1 = $f[0];
		$f2 = $f[1];
	}
	if ($f2 eq '' && $peFlag) {
		my $testFile = $f1;
		if ($testFile =~ s/_R1_/_R2_/ || $testFile =~ s/\.R1\./\.R2\./ || $testFile =~ s/_R1\./_R2\./) {
			if (-e $testFile) {
				print STDERR "\tFound PE file pair:\n\t\t$f1\n\t\t$testFile\n";
				$f2 = $testFile;
			}
		} elsif ($testFile =~ /[ES]RR.*_1\.fastq/) {
			$testFile =~ s/_1\.fastq/_2\.fastq/;
			if (-e $testFile) {
				print STDERR "\tUsing PE file: $f1 & $testFile\n";
				$f2 = $testFile;
			}
		}
		if ($f2 eq '') {
			print STDERR "!!! Warning, Could not find 2nd read for file:\n\t\t$f1\n";
		}
	}
	push(@file1, $f1);
	push(@file2, $f2);
}
print STDERR "\t\n";
print STDERR "\n\tSleeping for 3 seconds...\n";
`sleep 3`;

if ($fastpFlag) {
	my $check = `which fastp`;
	if ($check eq "") {
		print STDERR "\n\t!!! Could not detect fastp program in the executable path !!!\n";
		print STDERR "\t     Either install fastp or remove the -fastp option from the command.\n\n";
		exit;
	}
	my @newF1 = ();
	my @newF2 = ();
	print STDERR "\tTrimming input files with fastp:\n";


	my $cpus = 0;
	for (my $i=0;$i<@file1;$i++) {
		my $pe = 0;
		my $f1 = $file1[$i];
		my $cf1 = $f1 . ".clean.fq";
		if ($prefix ne '') {
			$cf1 = $prefix . "/" . basename($f1) . ".clean.fq";
		}
		my $cmd = "fastp -i \"$f1\" -o \"$cf1\"";
		my $f2 = $file2[$i];
		my $cf2 = $f2;
		if ($f2 ne '') {
			$cf2 = $f2 . ".clean.fq";
			if ($prefix ne '') {
				$cf2 = $prefix . "/" . basename($f2) . ".clean.fq";
			}
			$cmd .= " -I \"$f2\" -O \"$cf2\" ";
		}
		my $html = $f1 . ".fastp.html";
		$cmd .= " -h \"$html\"";
		print STDERR "\t\ttrimming $f1\n";
		print STDERR "\t\t         $f2\n" if ($f2 ne '');
		push(@newF1, $cf1);
		push(@newF2, $cf2);
		
		my $pid = fork();
		$cpus++;
		if ($pid == 0) {
			`$cmd`;
			exit(0);
		}
		if ($cpus >= $pCPUs) {
			wait();
			$cpus--;
		}
	}
	my $id = 0;
	while ($id >= 0) {
		$id = wait();
	}
	print STDERR "\tDone trimmiing reads with fastp\n";

	@file1 = @newF1;
	@file2 = @newF2;
}

my $genomeName = $bt2Index;
$genomeName =~ s/\/+$//;
$genomeName =~ s/^.+\///;
$genomeName =~ s/\.fa//;
$genomeName =~ s/^Genome_//;

print STDERR "\n";
print STDERR "\tUsing $program to perform alignment\n";
print STDERR "\tIndex = $bt2Index ($genomeName)\n";
print STDERR "\tNumber of instances at once: $maxCPUs\n";
print STDERR "\tNumber of cpus per instance: $pCPUs\n";
for (my $i=0;$i<@file1;$i++){ 
	if ($file2[$i] eq '') {
		print STDERR "\tSE: $file1[$i]\n";
	} else {
		print STDERR "\tPE: $file1[$i]\n";
		print STDERR "\t    $file2[$i]\n";
	}
}

if (`which "$exe"`) {
} else {
    print STDERR "\n!!! Error: Could not execute $program installed ($exe)\n";
	print STDERR "!!! Make sure it is in the executable PATH or specify path with -path <...>\n";
    exit;
}

if ($program eq 'STAR') {
	if ($removeSharedMemory && $genomeLoad ne "NoSharedMemory") {
		`$exe --genomeLoad Remove --genomeDir $bt2Index $sendErr2devnull`;
		`rm -f Log.progress.out Aligned.out.sam Log.out`;
		`rm -fr _STARtmp/`;
	}
}




my @pids = ();
my $cpus = 0;
my @sams = ();
for (my $j=0;$j<@file1;$j++) {
	my $file1 = $file1[$j];
	my $file2 = $file2[$j];

	my %delete = ();
	my $outputPrefix = $file1;
	if ($prefix ne '') {
		$outputPrefix = $prefix . "/" . basename($file1);
	}


	$outputPrefix =~ s/\.gz$//;
	$outputPrefix =~ s/\.bz2$//;
	$outputPrefix =~ s/\.fq$//;
	$outputPrefix =~ s/\.fastq$//;
	if ($suffix eq '') {
		$outputPrefix .= "." . $genomeName . "." . $program;
	} else {
		$outputPrefix .= "." . $suffix;
	}

	my $outputSam = $outputPrefix . ".sam";
	my $outputBam = $outputPrefix . ".bam";
	my $outputLog = $outputPrefix . ".log";
	push(@sams, $outputSam);

	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		print STDERR "\tAligning: $file1\n";
		if ($file2 ne '') {
			print STDERR "\t          $file2\n";
		}
		#child process


		if ($program eq 'tophat2') {
			print STDERR "\tScript doesn't work with Tophat... try STAR\n";
			#my $outputDir = "$file1.$genomeName.tophat2";
			#my $outputFile = "$file1.$genomeName.tophat2.bam";
			#my $outputJunc = "$file1.$genomeName.tophat2.junc";
			#my $logFile = "$file1.$genomeName.tophat2.log";
			#`$exe --library-type $libraryType -p $pCPUs -g $maxMultihits $guideGTF -o $outputDir $pass $maxMisMatches "$bt2Index" $file 2> $logFile`;
			#`mv "$outputDir/accepted_hits.bam" "$outputFile"`;
        	#`parseTophatJunctions.pl "$outputDir/junctions.bed" > "$outputJunc"`;
			#`rm -r "$outputDir"`;
			#if ($unFlag) {
			#	my $samFile = "$file1.$genomeName.tophat2.sam";
			#	`samtools view -h "$outputFile" > "$samFile" 2>> $logFile`;
			#	my $unalignedFile = "$file1.$genomeName.tophat2.unaligned.fq";
			#	my $unalignedFile2 = "$file2.$genomeName.tophat2.unaligned.fq";
			#	`getUnalignedReadsSam.pl $samFile $file1 > $unalignedFile 2>> $logFile`;
			#	if ($peflag) {
			#		`getUnalignedReadsSam.pl $samFile $file2 > $unalignedFile2 2>> $logFile`;
			#	}
			#	`rm "$samFile"`;
			#}
		} elsif ($program eq 'bowtie') {
			my $unFile = "";
			if ($unFlag) {
				$unFile = "--un $outputPrefix.unaligned.fq";
			}
			my $input = "\"$file1\"";
			if ($file2 ne '') {
				$input = "-1 \"$file1\" -2 \"$file2\"";
			}
			`$exe $local $unFile --best -M 1 --sam -p $pCPUs $pass "$bt2Index" $input > "$outputSam" 2> "$outputLog"`;
			#print STDERR "`$exe $local $unFile --best -M 1 --sam -p $pCPUs $pass -x $bt2Index $input > $outputSam 2> $outputLog`;\n";
			push(@sams, $outputSam);
		} elsif ($program eq 'bowtie2') {
			my $unFile = " --no-unal";
			if ($unFlag) {
				$unFile = "--un $outputPrefix.unaligned.fq";
			}
			my $input = "-U \"$file1\"";
			if ($file2 ne '') {
				$input = "-1 \"$file1\" -2 \"$file2\"";
			}
			`$exe $local $unFile -p $pCPUs $pass -x "$bt2Index" $input > "$outputSam" 2> "$outputLog"`;
			push(@sams, $outputSam);
		} elsif ($program eq 'hisat2') {
			my $unFile = " --no-unal";
			if ($unFlag) {
				$unFile = "--un $outputPrefix.unaligned.fq";
			}
			my $input = "-U \"$file1\"";
			if ($file2 ne '') {
				$input = "-1 \"$file1\" -2 \"$file2\"";
			}
			#print STDERR "$exe $hisat2Options $local $unFile -p $pCPUs $pass -x $bt2Index $input > $outputSam 2> $outputLog\n";
			`$exe $hisat2Options $local $unFile -p $pCPUs $pass -x "$bt2Index" $input > "$outputSam" 2> "$outputLog"`;
			push(@sams, $outputSam);

		} elsif ($program eq 'bwa') {
		
			my $ff2 = "";
			if ($file2 ne '') {
				$ff2 = "\"$file2\"";
			}
			`$exe mem -t $pCPUs "$bt2Index" "$file1" $ff2 > "$outputSam" 2> "$outputLog"`;
			push(@sams, $outputSam);

		} elsif ($program eq 'STAR') {
			my $opt = $starOptions;
			if ($unFlag) {
				$opt .= ' --outReadsUnmapped Fastx ';
			}
			if ($file1 =~ /\.gz$/) {
				$opt .= ' --readFilesCommand zcat ';
			} elsif ($file1 =~ /\.bz2$/) {
				$opt .= ' --readFilesCommand bzcat ';
			}

			if ($bamFlag) {
				$opt .= " --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM $bamSortRam ";
			}
			if ($pass2Flag) {
				$opt .= " --twopass1readsN 10000000000 ";
			}
	
			`$exe --genomeLoad $genomeLoad --genomeDir "$bt2Index" --runThreadN $pCPUs --readFilesIn "$file1" "$file2" --outFileNamePrefix "$outputPrefix." $opt $pass`;

			`mv "$outputPrefix.Log.final.out" "$outputLog"`;
			if ($bamFlag) {
				`mv "$outputPrefix.Aligned.sortedByCoord.out.bam" "$outputBam"`;
				`samtools index "$outputBam"`;
			} else {
				`mv "$outputPrefix.Aligned.out.sam" "$outputSam"`;
			}
			if ($unFlag) {
				if ($file2 ne '') {
					`mv "$outputPrefix.Unmapped.out.mate1" "$outputPrefix.unaligned.R1.fq"`;
					`mv "$outputPrefix.Unmapped.out.mate2" "$outputPrefix.unaligned.R2.fq"`;
				} else {
					`mv "$outputPrefix.Unmapped.out.mate1" "$outputPrefix.unaligned.fq"`;
				}
			}
			`rm "$outputPrefix.Log.out" "$outputPrefix.Log.progress.out" "$outputPrefix.SJ.out.tab"`;
			`rm -rf "$outputPrefix._STARtmp/"`;
			if ($pass2Flag) {
				`rm -r "$outputPrefix._STARpass1/"`;
			}


		} else {
			print STDERR "!!! Error: Program not recognized.\n";
		}

		foreach(keys %delete) {
			`rm "$_"`;
		}
		exit(0);
	}
	push(@pids, $pid);
	if ($cpus >= $maxCPUs) {
		my $id = wait();
		$cpus--;
	}
}
my $id = 0;
while ($id >= 0) {
	$id = wait();
	if ($id == -1) {
	} else {
	}
}
print STDERR "\n";

if ($bamFlag && $program ne 'STAR') {
	print STDERR "\tCreating bam files from sam files...\n";
	foreach (@sams) {
		my $samFile = $_;
		my $bamFile = $samFile;
		$bamFile =~ s/\.sam$//;
		#print STDERR "\t\tConverting $samFile\n";
		print STDERR "\t\tSAM->BAM $samFile\n";
		my $tmpBamFile = $bamFile . ".tmp.bam";
		#print STDERR "samtools view -b -@ $pCPUs $samFile > $tmpBamFile 2> /dev/null\n";
		#print STDERR "samtools sort -m $bamSortRam -@ $pCPUs -O BAM $tmpBamFile > $bamFile.bam 2> /dev/null\n";
		#print STDERR "samtools index $bamFile.bam\n";
		`samtools view -b -@ $pCPUs "$samFile" > "$tmpBamFile" $sendErr2devnull`;
		`samtools sort -m $bamSortRam -@ $pCPUs -O BAM "$tmpBamFile" > "$bamFile.bam" $sendErr2devnull`;
		`samtools index "$bamFile.bam"`;
		`rm "$tmpBamFile"`;
		`rm "$samFile"`;
	}
}
print STDERR "\n\tFinished.\n\n";

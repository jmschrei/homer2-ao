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

my $mode='stack';
my $frac = 0.333;
my $maxCPUs = 1;

sub printCMD {
	print STDERR "\n\tbatchMakeHiCMatrix.pl -pos <chr:start-end> -res <#> -window <#> [etc.] -d <HiCtagDir1> [HiCtagDir2] ...\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-d <HiC TagDir1> [HiC TagDir2] ... (Tag Directories of Hi-C experiments to visulize)\n";
	print STDERR "\t\t-pos <chr:start-end> (genomic position to visualize)\n";
	print STDERR "\t\t-res <#> (resolution of step size to use for analysis)\n";
	print STDERR "\t\t-window <#> (resolution of window size for aggregating interactions)\n";
	print STDERR "\t\t-balance (balance resulting Hi-C matrix)\n";
	print STDERR "\t\t-stack (Stacks matricies on top of one another i.e. square and symetric, non-rotated, default)\n";
	print STDERR "\t\t-split (Creates split matricies i.e. square, non-symetric, non-rotated)\n";
	print STDERR "\t\t\t" . '(printed in order of directories: 1\2 3\4 5\6 ...)' . "\n";
	print STDERR "\t\t-rotate (Rotates matrices, default)\n";
	print STDERR "\t\t\t-frac <#> (fraction of square matrix to consider for rotating, default: $frac)\n";
	print STDERR "\t\t-cpu (number of different processes to use, def: $maxCPUs)\n";
	print STDERR "\n\t\tOther options are passed to analyzeHiC to control the creation of the matrices\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 2) {
	printCMD();
}

my $opt = '';
my @tagDirs = ();

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-frac') {
		$frac = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rotate') {
		$mode = 'rotate';
	} elsif ($ARGV[$i] eq '-stack') {
		$mode = 'stack';
	} elsif ($ARGV[$i] eq '-split') {
		$mode = 'split';
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
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
		$opt .= " " . $ARGV[$i];
	}
}

my $rand = rand();
my @outputFiles = ();

my $cpus = 0;
print STDERR "\n\tCreating contact matrices:\n";
for (my $i=0;$i<@tagDirs;$i++) {
	my $tagDir = $tagDirs[$i];
	my $outputFile = "$rand-$i.txt";
	push(@outputFiles, $outputFile);

	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		`analyzeHiC $tagDir $opt > "$outputFile"`;
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


if ($mode eq 'stack') {
	for (my $i=0;$i<@outputFiles;$i++) {
		my $outFile = $outputFiles[$i];
		open IN, $outFile;
		my $c = 0;
		while (<IN>) {
			$c++;
			next if ($i>0 && $c == 1);
			print $_;
		}
		close IN;
		`rm $outFile`;
	}
} elsif ($mode eq 'rotate') {
	my $tmpFile= $rand . ".tmp";

	for (my $i=0;$i<@outputFiles;$i++) {
		my $outFile = $outputFiles[$i];
		`rotateHiCmatrix.pl $outFile $frac > $tmpFile`;
		open IN, $tmpFile;
		my $c = 0;
		while (<IN>) {
			$c++;
			next if ($i>0 && $c == 1);
			print $_;
		}
		close IN;
		`rm $tmpFile $outFile`;
	}
} elsif ($mode eq 'split') {
	for (my $i=0;$i<@outputFiles;$i+=2) {
		my $outFile1 = $outputFiles[$i];
		my $outFile2 = $outputFiles[$i+1];
		my @data1=();
		open IN, $outFile1;
		my $c = 0;
		while (<IN>) {
			$c++;
			if ($i==0 && $c == 1) {
				print $_;
			}
			chomp;
			s/\r//g;
			my @line = split /\t/;
			push(@data1, \@line);
		}
		close IN;

		my $matrixSize = $c-1;

		open IN, $outFile2;
		$c = 0;
		while (<IN>) {
			$c++;
			next if ($c ==1);
			chomp;
			s/\r//g;
			my @line = split /\t/;
			next if (@line < 3);
			print "$line[0]\t$line[1]";
			for (my $j=2;$j<=$c+1;$j++) {
				next if ($j >= scalar(@{$data1[$c-1]}));
				print "\t" . $data1[$c-1][$j];
			}
			for (my $j=$c+2;$j<@line;$j++) {
				print "\t$line[$j]";
			}
			print "\n";
		}
		close IN;
		`rm $outFile1 $outFile2`;
	}
} else {
	print STDERR "!!! Error - mode ($mode) not recognized!\n";
	exit;
}
print STDERR "\n\tOrder of HiC tag directories:\n";
for (my $i=0;$i<@tagDirs;$i++) {
	my $v = $i+1;
	print STDERR "\t\t$v\t$tagDirs[$i]\n";
}
print STDERR "\n\n";


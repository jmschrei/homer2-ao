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
# GNU General Public License for more details.

my @suffix = (
	".PC1.bedGraph",
	".ICF.bedGraph",
	".DLR.bedGraph"
);
my $prefix = '';
my $dir1 =  "";
my $dir2 = "";

sub printCMD {
	print STDERR "\n\tsubtractBedGraphDirectory.pl <HiCTagDir1> <HiCTagDir2> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-prefix <prefixName> (will separate bedGraphs based on this name)\n";
	print STDERR "\t\t-center (center bedGraphs on average signal, removes global differences)\n";
	print STDERR "\n";
	exit;
}
if (@ARGV < 2) {
	printCMD();
}

my $centerFlag = '';
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-prefix') {
		$prefix = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-center') {
		$centerFlag = ' -center';
	} else {
		if ($ARGV[$i] =~ /^\-/) {
			printCMD();
		}
		if ($dir1 eq '') {
			$dir1 = $ARGV[$i];
		} elsif ($dir2 eq '') {
			$dir2 = $ARGV[$i];
		} else {
			print STDERR "!!! Warning, 2 tagdirs already specified (dir1=$dir1,dir2=$dir2)\n";
			exit;
		}
	}
}
if ($dir1 eq '' || $dir2 eq '') {
	print STDERR "!!! Error, need to specify two tag directories\n";
	exit;
}

my $rand = rand();
my $tmpfile = $rand . ".tmp";

$list1 = `ls $dir1/*bedGraph`;
$list2 = `ls $dir2/*bedGraph`;
my @files1 = split /\n/,$list1;
my @files2 = split /\n/,$list2;
print STDERR "list1=@files1\n";
print STDERR "list2=@files2\n";

my $f1cov = '';
foreach(@files1){
	if (/\.COV\.bedGraph/) {
		$f1cov = " -cov $_";
		last;
	}
}

my @used = ();

foreach(@suffix) {
	my $suffix = $_;
	my $f1 = '';
	my $f2 = '';
	foreach(@files1){
		if (/$suffix/) {
			$f1 = $_;
			last;
		}
	}
	foreach(@files2){
		if (/$suffix/) {
			$f2 = $_;
			last;
		}
	}
	if ($f1 eq '' || $f2 eq '') {
		print STDERR "\tNo match for $suffix (f1=$f1,f2=$f2)\n";
	}
	print STDERR "Using: $f1\n";
	print STDERR "       $f2\n";
	push(@used, $f1, $f2);
	if ($prefix eq '') {
		`subtractBedGraphs.pl $f1 $f2 $f1cov $centerFlag >> $tmpfile`;
	} else {
		my $outfile = $prefix . $suffix;
		`subtractBedGraphs.pl $f1 $f2 $f1cov $centerFlag > "$outfile"`;
	}

}

print STDERR "Used:\n";
foreach(@used) {
	print STDERR "\t$_\n";
}
if ($prefix eq '') {
	open IN, $tmpfile;
	while (<IN>) {
		print $_;
	}
	close IN;
	`rm $tmpfile`;
}

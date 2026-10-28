#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";
#
#
## Copyright 2009 - 2019 Christopher Benner <cbenner@ucsd.edu>
##
## This file is part of HOMER
##
## HOMER is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## HOMER is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.


use POSIX;
use HomerConfig;

my $N = 1000;
my $eps = 0.01;

sub printCMD {
	print STDERR "\n\tusage: homer2jaspar.pl <homerMotifFile> [homerMotifFile2] ...\n";
	print STDERR "\n\tWrites a JASPAR formatted motif file to stdout\n";
	exit;
}

if (@ARGV < 1) {
	printCMD();
}

unless (-e $ARGV[0]) {
	printCMD();
}

for (my $z=0;$z<@ARGV;$z++) {
	my $file = $ARGV[$z];
	my $motifs = HomerConfig::readMotifFile($file);
	foreach(@$motifs) {
		my $name = $_->{'name'};
		my $cons = $_->{'consensus'};
		print ">$name\t$name\n";
		my $m = $_->{'matrix'};
	
		for (my $i=0;$i<4;$i++) {
			if ($i==0) {
				print "A";
			} elsif ($i==1) {
				print "C";
			} elsif ($i==2) {
				print "G";
			} else {
				print "T";
			}
			print "  [";
			foreach(@$m) {
				my $v = $_->[$i];
				$v *= $N;
				my $t = sprintf("%d",$v);
				my $nn = ceil(log($t+$eps)/log(10));
				$nn = 1 if ($nn < 1);
				for (my $j=$nn;$j<6;$j++) {
					print " ";
				}
				print "$t ";
			}
			print "]\n";
		}
	}
}

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

my $frac = 0.333;
my $scaleFactor = 1;
if (@ARGV < 1) {
	print STDERR "<hic matrix> [fraction of height|$frac]\n";
	exit;
}
if (@ARGV > 1) {
	$frac = $ARGV[1];
}
if (@ARGV > 2) {
	$scaleFactor = $ARGV[2];
}
#print STDERR "\tFraction = $frac\n";
#print STDERR "\tScaleFactor = $scaleFactor\n";
my $size = -1;
my $chr = '';
my $start = -1;
my $end = -1;
my @data = ();


open IN, $ARGV[0];
my $c = 0;
while (<IN>) {
	$c++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($c==1) {
		$size = scalar(@line)-2;
		if ($line[2] =~ /(.*?)\-(\d+)/) {
			$chr = $1;
			$start = $2;
		}
		if ($line[@line-1] =~ /(.*?)\-(\d+)/) {
			$chr = $1;
			$end = $2;
		}
		if ($size < 0 || $chr eq '' || $start < 0 || $end < 0) {
			print STDERR "!!! Error - something was wrong reading the file:\n";
			print STDERR "\tsize=$size\n\tchr=$chr\n\tstart=$start\n\tend=$end\n\n";
			exit;
		}
		$bin = ($end-$start)/($size-1);
		$offsetStart = floor($size*$frac);
		$halfOffsetStart = ceil($offsetStart/2);
		$offsetEnd = $size-$halfOffsetStart;
		$newStart = $start+$halfOffsetStart*$bin;
		$newEnd = $start+$offsetEnd*$bin;
		$newSize = $offsetEnd-$halfOffsetStart;

		#print STDERR "\tsize=$size\n\tchr=$chr\n\tstart=$start\n\tend=$end\n\tbin=$bin\n";
		print STDERR "\tRotation range adjustment: New position is $chr:$newStart-$newEnd\n";
		next;
	}
	my $p = shift @line;
	$p = shift @line;
	push(@data, \@line);
}
close IN;

#print STDERR "\toffsetStart=$offsetStart\n";
#print STDERR "\thalfOffsetStart=$halfOffsetStart\n";
#print STDERR "\tnewSize = $newSize\n";

print "rotated\trotated";
my @newdata = ();
for (my $i=0;$i<=$offsetStart;$i++) {
	my $osx = $offsetStart-floor($i/2);
	my $osy = ceil($i/2);
	my $stagger = 0;
	if ($osy != $i/2) {
		$stagger=1;
	}
	#print STDERR "\ti=$i osx=$osx osy=$osy stagger=$stagger\n";

	my @d = ();
	for (my $j=0;$j<$newSize*2;$j++) {
		if ($i==0) {
			my $p = "$chr" . "-" . ($newStart+$j*floor($bin*0.5));
			print "\t$p";
		}
		push(@d, 0);
	}
	for (my $j=0;$j<$newSize;$j++) {
		#print STDERR "\t\tosx=$osx osy=$osy j=$j\n";
		$d[$j*2+$stagger] = $data[$osy+$j][$osx+$j];
	}
	my $startPos = $offsetStart-$i;
	my $nn = scalar(@d);
	push(@newdata, \@d);
}
print "\n";



if (0) {
for (my $i=0;$i<@newdata;$i++) {
		my $n = scalar(@{$newdata[$i]});
		for (my $j=1;$j<$n;$j++) {
		#for (my $j=0;$j<@{$newdata[$i]};$j++) {
			print STDERR "\t$newdata[$i][$j]";
		}
		print STDERR "\n";
	}
}

for (my $i=0;$i<@newdata;$i++) {
	my $n = scalar(@{$newdata[$i]});
	#print STDERR "$i\t$n\n";
	for (my $j=0;$j<@{$newdata[$i]};$j++) {
		if ((($i % 2 == 1) && ($j % 2 == 0)) || (($i % 2 == 0) && ($j % 2 == 1))) {
			my $sum = 0;
			my $N = 0;
			if ($i > 0) {
				#$sum += $newdata[$i-1][$j];
				#$N++;
			}
			if ($i < scalar(@newdata)-1) {
				#$sum += $newdata[$i+1][$j];
				#$N++;
			}
			if ($j > 0) {
				$sum += $newdata[$i][$j-1];
				$N++;
			}
			if ($j < scalar(@{$newdata[$i]})-1) {
				$sum += $newdata[$i][$j+1];
				$N++;
			}
			if ($N > 0) {
				$sum /= $N;
			}
			$newdata[$i][$j] = $sum;
		}
		
	}
}

$c = -1;
foreach(@newdata) {
	my $data = $_;
	$c++;
	my $d=  floor($halfOffsetStart)*$bin - $c*floor($bin/2);
	print "Dist-$d\tDist-$d";
	foreach(@$data) {
		my $v = $_ * $scaleFactor;
		print "\t$v";
	}
	print "\n";
}

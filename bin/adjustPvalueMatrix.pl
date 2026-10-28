#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";

#
##
## Copyright 2009 - 2024 Christopher Benner <cbenner@ucsd.edu>
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
##



use Statistics;

my $inputFile = '';
my $numHeader = 1;
my $numCol = 1;
my $logpFlag = 0;
my $inbase = 10;
my $outbase = 10;
my $homerFlag = 0;

sub printCMD {
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-i <inputFile>\n";
	print STDERR "\t\t-logp (input values are logp values)\n";
	print STDERR "\t\t\t-inbase <#> (log base for logp - def: 10)\n";
	print STDERR "\t\t\t-outbase <#> (log base for logp - def: 10)\n";
	print STDERR "\t\t-flip (flip direction on logp values)\n";
	print STDERR "\t\t-homer (Use homer's BH correction instead of R)\n";
	print STDERR "\n";
	exit;
}

my $flipFlag = 0;
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] !~ /^\-/) {
		$inputFile = $ARGV[$i];
	} elsif ($ARGV[$i] eq '-i') {
		$inputFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-logp') {
		$logpFlag = 1;
	} elsif ($ARGV[$i] eq '-flip') {
		$flipFlag = 1;
	} elsif ($ARGV[$i] eq '-homer') {
		$homerFlag = 1;
	} elsif ($ARGV[$i] eq '-inbase') {
		$inbase = $ARGV[++$i];
		if ($inbase eq 'e') {
			$inbase = exp(1);
		}
	} elsif ($ARGV[$i] eq '-outbase') {
		$outbase = $ARGV[++$i];
		if ($outbase eq 'e') {
			$outbase = exp(1);
		}
	} else {
		printCMD();
	}
}
				
if ($inputFile eq '') {
	print STDERR "!!! Error: Need an input file\n";	
	printCMD();
}

if ($logpFlag) {
	print STDERR "\tTreating p-values as logp (base=$inbase, will output base=$outbase)\n";
} else {
	print STDERR "\tTreating p-values as normal pvalues (i.e. not logp)\n";
}


my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";
open OUT, ">$tmpFile";
print OUT "pvalues\n";


my @pvalues = ();
open IN, $inputFile;
my $r = 0;
while (<IN>) {
	$r++;
	chomp;
	s/\r//g;
	if ($r <= $numHeader) {
		next;
	}
	my @line = split /\t/;
	for (my $i=$numCol;$i<@line;$i++) {
		my $v = $line[$i];
		if ($logpFlag) {
			if ($v > 0) {
				$v = $inbase**(-1*$v);
			} else {
				$v = $inbase**($v);
			}
		}
		push(@pvalues,$v);
		print OUT "$v\n";
	}
}
close IN;
close OUT;
open OUT, ">$tmpFile2";
print OUT "x <- read.delim(\"$tmpFile\")\n";
print OUT "xx <- p.adjust(t(x), method=\"BH\")\n";
print OUT "write.csv(xx,file=\"$tmpFile3\")\n";
close OUT;

my $min = 1;
my @adjusted = ();

if ($homerFlag) {
	my $xx = Statistics::benjaminiFDR(\@pvalues);
	@adjusted = @$xx;
	foreach(@adjusted) {
		if ($_ > 0 && $_ < $min) {
			$min = $_;
		}
	}
} else {
	`R --no-save < $tmpFile2`;
	my $z = 0;
	open IN, $tmpFile3;
	while (<IN>) {
		$z++;
		chomp;
		s/\r//g;
		next if ($z < 2);
		my @line = split /\,/;
		my $r = $line[1];
		if ($r > 0 && $r < $min) {
			$min = $r;
		}
		push(@adjusted, $r);
	}
	close IN;
}

my $warningFlag = 0;

my $x = 0;
open IN, $inputFile;
$r = 0;
while (<IN>) {
	my $og = $_;
	$r++;
	chomp;
	s/\r//g;
	if ($r <= $numHeader) {
		print $og;
		next;
	}
	my @line = split /\t/;
	for (my $i=0;$i<$numCol;$i++) {
		print "\t" if ($i>0);
		print "$line[$i]";
	}
	for (my $i=$numCol;$i<@line;$i++) {
		my $ov = $line[$i];
		my $v = $adjusted[$x++];
		if ($logpFlag) {
			if ($v < $min) {
				#$v = $min;
				if ($warningFlag == 0) {
					print STDERR "\tWarnning - some logp values are so significant that they won't be corrected\n";
					$warningFlag = 1;
				}
				$v = -1*abs($ov);
				$v *= log($inbase)/log($outbase);
			} else {
				$v = log($v) / log($outbase);
			}
			if ($ov > 0) {
				$v *= -1;
			}
		}
		print "\t" if ($i>0);

		if ($flipFlag) {
			$v *= -1;
		}
		print "$v";
	}
	print "\n";
}
close IN;
`rm $tmpFile $tmpFile2 $tmpFile3`;

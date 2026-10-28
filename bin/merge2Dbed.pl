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


sub printCMD {
	print STDERR "\n\tmerge2Dbed.pl [options] <2D BED file1> <2D BED file2> [2D BED file3]...\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-res <#> (maximum distance between endpoints to merge, default: 15000)\n";
	print STDERR "\t\t\tUsually for loops -res should be set to the window/superRes size, for TADs 2x window/superRes\n";
	print STDERR "\t\t-loop (treat 2D bed input files as loops, default)\n";
	print STDERR "\t\t-tad (treat 2D bed input files as TADs)\n";
	print STDERR "\t\t-prefix <filePrefix> (output venn diagram overlaps to separate files)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 2) {
	printCMD();
}
my @intFiles = ();
my $minRes = '';
$tadFlag = -1;
my $commonOnly=0;
my $prefix = '';
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-dist' || $ARGV[$i] eq '-res' || $ARGV[$i] eq '-window' || $ARGV[$i] eq '-superRes') {
		$minRes = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-tad' || $ARGV[$i] eq '-TAD') {
		$tadFlag = 1;
	} elsif ($ARGV[$i] eq '-cp' || $ARGV[$i] eq '-loop' || $ARGV[$i] eq '-loops') {
		$tadFlag = 0;
	} elsif ($ARGV[$i] eq '-commonOnly') {
		$commonOnly=1;
	} elsif ($ARGV[$i] eq '-prefix') {
		$prefix = $ARGV[++$i];
	} elsif ($ARGV[$i] =~ /^-/) {
		print STDERR "!!! Error: $ARGV[$i] not recognized!!!\n";
		printCMD();
	} else {
		push(@intFiles, $ARGV[$i]);
	}
}

my $numFiles = scalar(@intFiles);
if ($tadFlag == -1) {
	$tadFlag = 0;
	for (my $i=0;$i<@intFiles;$i++) {
		my $v = check2Dbed($intFiles[$i]);
		if ($v < 0) {
			print STDERR "\tWarning: $intFiles[$i] may not have any valid 2D entries...\n";
		} elsif ($v eq '1') {
			print STDERR "\tFile $intFiles[$i] looks like a TAD 2D BED file\n";
			$tadFlag++;
		} else {
			print STDERR "\tFile $intFiles[$i] looks like a loop 2D BED file\n";
		}
	}
	if ($tadFlag > 0) {
		if ($tadFlag < $numFiles) {
			print STDERR "\t!!! Warning: only $tadFlag of $numFiles look like TAD files\n";
		} else {
			print STDERR "\tTreating files as TAD 2D BED files ($tadFlag of $numFiles look like TAD files)\n";
		}
		$tadFlag = 1;
	} else {
		print STDERR "\tTreating files as loop 2D BED files\n";
		$tadFlag = 0;
	}
} elsif ($tadFlag == 1) {
	print STDERR "\tTreating files as TAD 2D BED files (-tad)\n";
} else {
	print STDERR "\tTreating files as loop 2D BED files (-loop)\n";
}

if ($minRes eq '') {
	if ($tadFlag) {
		$minRes = 30000;
		print STDERR "\n\tMerging resolution set to $minRes (default for TAD files)\n";
	} else {
		$minRes = 15000;
		print STDERR "\n\tMerging resolution set to $minRes (default for loop files)\n";
	}
} else {
	print STDERR "\n\tMerging resolution set to $minRes (user defined)\n";
}
print STDERR "\n";

my %peaks = ();
my $id1 = 1;
my $id2 = 1;
my $index = 0;


my %toCombine =();
$extraHeader = '';

my @bed = ();
for (my $i=0;$i<@intFiles;$i++) {
	my $bed = read2Dbed($intFiles[$i],$i,$numFiles);
	push(@bed, $bed);
}

my $newInt=0;
my @newInts = ();
foreach(keys %peaks) {
	# for each chr
	my $c = $_;
	my @peaks = sort {$a->{'s'} <=> $b->{'s'}} values %{$peaks{$c}};

	my $maxIter = 5;
	for (my $z = 0;$z<$maxIter;$z++) {
		my $numMerge = 0;
		for (my $i=0;$i<@peaks;$i++) {
			my $id1 = $peaks[$i]->{'id'};
			my $s1 = $peaks[$i]->{'s'};
			my $e1 = $peaks[$i]->{'e'};
			my $i1 = $peaks[$i]->{'i'};
			next if ($peaks[$i]->{'g'}==0);
			my $m1 = floor(($s1+$e1)/2);
			for (my $j=$i+1;$j<@peaks;$j++) {
				my $id2 = $peaks[$j]->{'id'};
				my $s2 = $peaks[$j]->{'s'};
				my $e2 = $peaks[$j]->{'e'};
				my $i2 = $peaks[$j]->{'i'};
				my $m2 = floor(($s2+$e2)/2);
				next if ($peaks[$j]->{'g'}==0);
				next if ($i1 == $i2);
				my $overlap =0;
				if ($e1 >= $s2) {
					#$overlap=1;
				}
				if (abs($m2-$m1) < $minRes) {
					$overlap=1;
				}
				if ($overlap==0 && $s2-$e1 > $minRes*4) {
					last;
				}
				if ($overlap) {
					my $int1=$bed[$i1]->{$c}->{$id1};
					my $int2=$bed[$i2]->{$c}->{$id2};
					my ($good,$int) = compareInts($int1,$int2,$minRes);
					if ($good) {
						$numMerge++;
						push(@newInts, $int);
						$newInt++;
	
						$bed[$i1]->{$c}->{$id1} = $int;
						$peaks{$c}->{$int1->{'p1'}}->{'g'}=1;
						$peaks{$c}->{$int1->{'p2'}}->{'g'}=1;
						$peaks{$c}->{$int1->{'p1'}}->{'s'} = $int->{'s1'};
						$peaks{$c}->{$int1->{'p1'}}->{'e'} = $int->{'e1'};
						$peaks{$c}->{$int1->{'p2'}}->{'s'} = $int->{'s2'};
						$peaks{$c}->{$int1->{'p2'}}->{'e'} = $int->{'e2'};
	
						$peaks{$c}->{$int2->{'p1'}}->{'g'}=0;
						$peaks{$c}->{$int2->{'p2'}}->{'g'}=0;
						if (0) {
							$peaks{$c}->{$int1->{'p1'}}->{'g'}=0;
							$peaks{$c}->{$int1->{'p2'}}->{'g'}=0;
							$peaks{$c}->{$int2->{'p1'}}->{'g'}=0;
							$peaks{$c}->{$int2->{'p2'}}->{'g'}=0;
						}
					}
				}
			}
		}
		#print STDERR "\t$c\t$z\t$numMerge\n";
		last if ($numMerge == 0);
	}
}
my %counts = ();
my $c = 0;
my %sets = ();
for (my $i=0;$i<@intFiles;$i++) {
	foreach(values %{$bed[$i]}) {
		foreach(values %$_) {
			if ($_->{'g'}) {
				my $memStr = getMemStr($_->{'membership'});
				$counts{$memStr}++;
				if (!exists($sets{$memStr})) {
					my @a = ();
					$sets{$memStr} = \@a;
				}
				push(@{$sets{$memStr}},$_);
				$c++;
			}
		}
	}
}
if ($prefix eq '') {
	print "#merged=$c chr1\tstart1\tend1\tchr2\tstart2\tend2" . $extraHeader . "\n";
}
print STDERR "\t$c total features after merging\n";
print STDERR "\nFeatures";
foreach(@intFiles) {
	print STDERR "\t$_";
}
print STDERR "\tName\n";
foreach(keys %counts) {

	my $name = '';	
	my $file = '';	
	print STDERR "$counts{$_}";
	for (my $i=0;$i<$numFiles;$i++) {
		my $substr = substr($_,$i,1);
		if ($substr eq '1') {
			print STDERR "\tX";
			$name .= "|" . $intFiles[$i];
			$file .= "_" . $intFiles[$i];
		} else {
			print STDERR "\t";
		}
	}
	$name =~ s/^\|//;
	print STDERR "\t$name\n";

	my $fh=*STDOUT;
	if ($prefix ne '') {
		my $outFile = $prefix . $file;
		$outFile =~ s/\//_/g;
		open ($fh,">",$outFile);
		print $fh "#merged=$counts{$_} chr1\tstart1\tend1\tchr2\tstart2\tend2" . $extraHeader . "\n";
	}
	foreach(@{$sets{$_}}) {
		printInt($_);
	}
	if ($prefix ne '') {
		close $fh;
	}
}
exit;









sub getMemStr {
	my ($mem)  = @_;
	my $str = '';
	foreach(@$mem) {
		$str .= $_;
	}
	return $str;
}
sub printInt {
	my ($int) = @_;
	my $s1 = $int->{'s1'};
	my $e1 = $int->{'e1'};
	my $s2 = $int->{'s2'};
	my $e2 = $int->{'e2'};
	if ($tadFlag) {
		$e1 = $e2;
		$s2 = $s1;
	}
	print "$int->{'c1'}\t$s1\t$e1\t$int->{'c2'}\t$s2\t$e2";
	foreach(@{$int->{'other'}}) {
		print "\t$_";
	}
	print "\n";
}

sub compareInts {
	my ($int1,$int2,$minRes) = @_;
	my $good = 1;
	my $flip=0;
	my $int = '';

	if ($int1->{'c1'} eq $int2->{'c1'}) {
		if ($int1->{'c2'} eq $int2->{'c2'}) {
			#TADs, Ints on same chromosome
			if ( (($int1->{'s1'} <= $int2->{'e1'} && $int1->{'e1'} >= $int2->{'s1'}) 
			  			|| abs($int1->{'m1'} - $int2->{'m1'}) <= $minRes)
					&& (($int1->{'s2'} <= $int2->{'e2'} && $int1->{'e2'} >= $int2->{'s2'}) 
			  			|| abs($int1->{'m2'} - $int2->{'m2'}) <= $minRes)) {
				$good = 1;
				my $s1 = $int1->{'s1'};
				$s1 = $int2->{'s1'} if ($int2->{'s1'} < $s1);
				my $e1 = $int1->{'e1'};
				$e1 = $int2->{'e1'} if ($int2->{'e1'} > $e1);
				my $s2 = $int1->{'s2'};
				$s2 = $int2->{'s2'} if ($int2->{'s2'} < $s2);
				my $e2 = $int1->{'e2'};
				$e2 = $int2->{'e2'} if ($int2->{'e2'} > $e2);
				$int = {c1=>$int1->{'c1'},s1=>$s1,e1=>$e1,m1=>floor(($s1+$e1)/2),c2=>$int1->{'c2'},s2=>$s2,e2=>$e2,
								m2=>floor(($s2+$e2)/2),other=>'',g=>1,
								p1=>$int1->{'p1'}, p2=>$int1->{'p2'},membership=>$int1->{'membership'}};
				$int1->{'g'} = 0;
				$int2->{'g'} = 0;
			} elsif ( (($int1->{'s1'} <= $int2->{'e2'} && $int1->{'e1'} >= $int2->{'s2'}) 
			  			|| abs($int1->{'m1'} - $int2->{'m2'}) <= $minRes)
					&& (($int1->{'s2'} <= $int2->{'e1'} && $int1->{'e2'} >= $int2->{'s1'}) 
			  			|| abs($int1->{'m2'} - $int2->{'m1'}) <= $minRes) && ($int1->{'c1'} eq $int1->{'c2'})) {
				$good = 1;
				$flip=1;
				my $s1 = $int1->{'s1'};
				$s1 = $int2->{'s2'} if ($int2->{'s2'} < $s1);
				my $e1 = $int1->{'e1'};
				$e1 = $int2->{'e2'} if ($int2->{'e2'} > $e1);
				my $s2 = $int1->{'s2'};
				$s2 = $int2->{'s1'} if ($int2->{'s1'} < $s2);
				my $e2 = $int1->{'e2'};
				$e2 = $int2->{'e1'} if ($int2->{'e1'} > $e2);
				$int = {c1=>$int1->{'c1'},s1=>$s1,e1=>$e1,m1=>floor(($s1+$e1)/2),c2=>$int1->{'c2'},s2=>$s2,e2=>$e2,
								m2=>floor(($s2+$e2)/2),other=>'',g=>1,
								p1=>$int1->{'p1'}, p2=>$int1->{'p2'},membership=>$int1->{'membership'}};
				$int1->{'g'} = 0;
				$int2->{'g'} = 0;
			} else {
				$good=0;
			}
		} else {
			$good=0;
		}
	} elsif ($int1->{'c1'} eq $int2->{'c2'}) {
		if ($int1->{'c2'} eq $int2->{'c1'}) {
			if ( (($int1->{'s1'} <= $int2->{'e2'} && $int1->{'e1'} >= $int2->{'s2'}) 
			  			|| abs($int1->{'m1'} - $int2->{'m2'}) <= $minRes)
					&& (($int1->{'s2'} <= $int2->{'e1'} && $int1->{'e2'} >= $int2->{'s1'}) 
			  			|| abs($int1->{'m2'} - $int2->{'m1'}) <= $minRes)) {
				$good = 1;
				$flip=1;
				my $s1 = $int1->{'s1'};
				$s1 = $int2->{'s1'} if ($int2->{'s1'} < $s1);
				my $e1 = $int1->{'e1'};
				$e1 = $int2->{'e1'} if ($int2->{'e1'} > $e1);
				my $s2 = $int1->{'s2'};
				$s2 = $int2->{'s2'} if ($int2->{'s2'} < $s2);
				my $e2 = $int1->{'e2'};
				$e2 = $int2->{'e2'} if ($int2->{'e2'} > $e2);
				$int = {c1=>$int1->{'c1'},s1=>$s1,e1=>$e1,m1=>floor(($s1+$e1)/2),c2=>$int1->{'c2'},s2=>$s2,e2=>$e2,
								m2=>floor(($s2+$e2)/2),other=>'',g=>1,
								p1=>$int1->{'p1'}, p2=>$int1->{'p2'},membership=>$int1->{'membership'}};
				$int1->{'g'} = 0;
				$int2->{'g'} = 0;
			} else {
				$good = 0;
			}
		} else {
			$good=0;
		}
	} else {
		$good=0;
	}
	if ($good) {
		for (my $i=0;$i<@{$int->{'membership'}};$i++) {
			if ($int2->{'membership'}->[$i] > 0) {
				$int->{'membership'}->[$i] = 1;
			}
		}
		$int->{'other'} = $int1->{'other'};
		for (my $i=1;$i<@{$int1->{'other'}};$i++) {
			my $max = $int1->{'other'}->[$i];
			#if (defined($int2->{'other'}->[$i]) && $max < $int2->{'other'}->[$i]) {
			#	$max = $int2->{'other'}->[$i];
			#}
			$int->{'other'}->[$i] = $max;
		}
	}
	return ($good,$int);
}





sub check2Dbed {
	my ($file) = @_;
	open IN, $file or die "Could not open file: $file\n";
	my $same=0;
	my $total=0;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if ($line[0] =~ /^#/);
		next if (@line < 6);
		if (($line[1] eq $line[4]) && ($line[0] eq $line[3]) && ($line[2] eq $line[5])) {
			$same++;
		}
		$total++;
	}
	close IN;
	return -1 if ($total < 1);
	return 1 if ($same/$total > 0.5);
	return 0;
}

sub read2Dbed {
	my ($file,$index,$numFiles) = @_;
	my %bed = ();
	$extraHeader = '';
	my $c = 0;
	open IN, $file or die "Could not open file: $file\n";
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($c == 1) {
			for (my $i=6;$i<@line;$i++){ 
				if ($line[0] =~ /^#/) {
					$extraHeader .= "\t" . $line[$i];
				} else {
					$extraHeader .= "\tinfo";
				}
			}
		}
		next if (/^#/);
		my $c1 = $line[0];
		my $s1 = $line[1];
		my $e1 = $line[2];
		my $m1 = floor(($s1+$e1)/2);
		my $c2 = $line[3];
		my $s2 = $line[4];
		my $e2 = $line[5];
		if ($c1 eq $c2) {
			#make sure first coordinate is lower one
			if ($s1 > $s2) {
				my $tmpS = $s2;
				$s2 = $s1;
				$s1 = $tmpS;
				my $tmpE = $e2;
				$e2 = $e1;
				$e1 = $tmpE;
			}
		}


		if ($tadFlag) {
			$e1 = $s1+$minRes;
			$s2 = $e2-$minRes;
		}
		my $m2 = floor(($s2+$e2)/2);
		my $id = "2d-" . $id1++;
		my $id1 = "p-" . $id2++;
		my $id2 = "p-" . $id2++;
		my @other = ();
		for (my $i=6;$i<@line;$i++) {
			push(@other, $line[$i]);
		}
		my $c = $c1;
		if (($c1 cmp $c2) > 0) {
			$c = $c2;
		}
		if (!exists($bed{$c})) {
			my %a = ();
			$bed{$c}=\%a;
		}
		my @membership = ();
		for (my $i=0;$i<$numFiles;$i++) {
			push(@membership,0);
		}
		$membership[$index] = 1;
		$bed{$c}->{$id} = {c1=>$c1,s1=>$s1,e1=>$e1,m1=>$m1,c2=>$c2,s2=>$s2,e2=>$e2,m2=>$m2,other=>\@other,g=>1,p1=>$id1,p2=>$id2,membership=>\@membership};
		if (!exists($peaks{$c1})) {
			my %a = ();
			$peaks{$c1}=\%a;
		}
		if (!exists($peaks{$c2})) {
			my %a = ();
			$peaks{$c2}=\%a;
		}
		$peaks{$c1}->{$id1} = {c=>$c1,s=>$s1,e=>$e1,id=>$id,i=>$index,g=>1};
		$peaks{$c2}->{$id2} = {c=>$c2,s=>$s2,e=>$e2,id=>$id,i=>$index,g=>1};
	}
	close IN;

	return \%bed;
}

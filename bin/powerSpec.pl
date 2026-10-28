#!/usr/bin/env perl
use warnings;

use POSIX;

if (@ARGV < 1) {
	print STDERR "\n\tUsage:\n";
	print STDERR "\t\tpowerSpec.pl <data file> [options]\n";
	print STDERR "\tThe input is a tab file (with header):\n";
	print STDERR "\t\tX1  Y1\n";
	print STDERR "\t\tX2  Y2\n";
	print STDERR "\tThe idea is to specify which column holds the y-value to analyze with the [col] paramter.\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-col <#> (Column with data to analyze, (1-based) 1,2,3,... default: 2)\n";
	print STDERR "\t\t-start <#> (minimum range, default: -100)\n";
	print STDERR "\t\t-end <#> (maximum range, default: 100)\n";
	print STDERR "\t\t-min <#> (minimum period to analyze, default: 1.5)\n";
	print STDERR "\t\t-max <#> (maximum period to analyze, default: 50)\n";
	print STDERR "\t\t-inc <#> (increment in period to analyze, default: 0.1)\n";
	print STDERR "\t\t-0 (Ignore data at the origin)\n";
	print STDERR "\t\t-vectorNorm (Normalize data vector to 1)\n";
	print STDERR "\t\t-totalNorm (Normalize total signal to 1)\n";
	print STDERR "\tDetails: Mean is subtracted from the data before analysis\n";
	print STDERR "\t\tUnits: input units^2\n";
	print STDERR "\n";
	exit;
}
my $useFullPeriods = 0;
my $window = 30;
my $autocorr = 100;

my $minPos= -100;
my $maxPos= 100;
my $col = 2;


my $minPeriod = 1.5;
my $maxPeriod = 50;
my $inc = 0.1;
my $PI = 3.14159265;
my $vectorNorm = 0;
my $totalNorm = 0;
my $zeroFlag = 0;

my $datafile = $ARGV[0];
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-col') {
		$col = $ARGV[++$i];
		print STDERR "\tAnalyzing frequency spectrum on column #$col\n";
	} elsif ($ARGV[$i] eq '-start') {
		$minPos = $ARGV[++$i];
		print STDERR "\tMinimum position analyzed for fourier analysis: $minPos\n";
	} elsif ($ARGV[$i] eq '-end') {
		$maxPos = $ARGV[++$i];
		print STDERR "\tMaximum position analyzed for fourier analysis: $maxPos\n";
	} elsif ($ARGV[$i] eq '-min') {
		$minPeriod = $ARGV[++$i];
		print STDERR "\tMinimum period analyzed: $minPeriod\n";
	} elsif ($ARGV[$i] eq '-max') {
		$maxPeriod = $ARGV[++$i];
		print STDERR "\tMaximum period analyzed: $maxPeriod\n";
	} elsif ($ARGV[$i] eq '-inc') {
		$inc = $ARGV[++$i];
		print STDERR "\tIncrement in period set to: $inc\n";
	} elsif ($ARGV[$i] eq '-vectorNorm') {
		$vectorNorm = 1;
		print STDERR "\tPerforming Vector Normalization\n";
	} elsif ($ARGV[$i] eq '-0') {
		$zeroFlag = 1;
		print STDERR "\tIgnoring data at the origin\n";
	} elsif ($ARGV[$i] eq '-totalNorm') {
		$totalNorm = 1;
		print STDERR "\tPerforming Total Signal Normalization\n";
	} else {
		print STDERR "!!!! Do not understand: $ARGV[$i]\n";
		exit;
	}
}

$col--;

my %data = ();
open IN, $ARGV[0];
my $sum = 0;
my $N = 0;
my $count = 0;
my $accMin = 1e100;
my $accMax = -1e100;
my $maxHeight = 1;
while (<IN>) {
	$count++;
	next if ($count ==1);
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($line[0] eq '0' && $zeroFlag == 1) {
		next;
	}
	next if ($line[0] < $minPos);
	next if ($line[0] > $maxPos);
	if ($line[0] < $accMin) {
		$accMin = $line[0];
	}
	if ($line[0] > $accMax) {
		$accMax = $line[0];
	}
	my $v=$line[$col];
	$data{$line[0]} = $v;
	$sum+=$v;
	$N++;
	$maxHeight = $v if ($v > $maxHeight);
}
close IN;


########## Normalization of signal

if ($totalNorm == 1) {
	my $sum = 0;
	foreach (values %data) {
		$sum += $_;
	}
	foreach (values %data) {
		$_ /= $sum;
	}
print STDERR "Hello!\n";
}

$sum = 1 if ($sum == 0);
$N = 1 if ($N == 0);


my $normSum = 0;
foreach(values %data) {
	#$_ /= $maxHeight;
	$normSum += $_;
}
my $avg = $normSum/$N;
foreach(values %data) {
	$_ -= $avg;
}
#print STDERR "\t$normSum";
if ($vectorNorm==1) {
	my $sum = 0;
	foreach(values %data) {
		$sum+=$_*$_;
	}
	$sum = sqrt($sum);
	foreach (values %data) {
		$_ /= $sum;
	}
}


# smoothing window
if (0) {
	my @x = sort {$a <=> $b} keys %data;
	my $win = ceil($window/2);
	my %newData = ();
	for (my $i=0;$i<@x;$i++) {
		my $avg = $data{$x[$i]};
		my $n = 1;
		for (my $j=$i-1;$j>$i-$win;$j--) {
			last if ($j<0);
			$avg += $data{$x[$j]};
			$n++;
		}
		for (my $j=$i+1;$j<$i+$win;$j++) {
			last if ($j >= @x);
			$avg += $data{$x[$j]};
			$n++;
		}
		$avg /= $n;
		$newData{$x[$i]} = $data{$x[$i]} - $avg;
	}
	%data = %newData;
}

#code for autocorrelation - not used
if (0) {
	@x = sort {$a <=> $b} keys %data;
	my %auto = ();
	my %autoN = ();
	for (my $i=0;$i<@x;$i++) {
		for (my $j=0;$j<$i;$j++) {
			my $diff = $i-$j;
			my $v = $data{$x[$i]}*$data{$x[$j]};
			$auto{$diff} += $v;
			$autoN{$diff} += $N;
		}
	}
	@x = sort {$a <=> $b} keys %auto;
	%data = ();
	foreach(@x) {
		$auto{$_} /= $autoN{$_};
		next if (abs($_) > $autocorr);
		$data{$_} = $auto{$_};
		print "$_\t$auto{$_}\t$autoN{$_}\n";
	}
}

@x = sort {$a <=> $b} keys %data;
my $min = $x[0];
my $max = $x[@x-1];
my $interval = $max-$min;


print "Period\tPower\tPhase\tFreq\n";
for (my $p=$minPeriod;$p<=$maxPeriod;$p+=$inc) {


	my $f = 2*$PI/$p;

	my $coverage = $interval/$p;
	my $curMax = 1e10;
	#my $curMax = $p*floor($coverage);
	if ($useFullPeriods == 1) {
		#$curMax = $max+1;
		#$curMax = floor($interval/$p)*$p;
	}

	#used for checking the effect of an offset
	my $theta = 0;

	my $F = 0;
	my $I = 0;
	my $maxX = -1e1000;

	my $P = 0;
	foreach(keys %data) {
		my $x = $_;
		next if ($x >= $curMax);
		if ($x > $maxX) {
			$maxX = $x;
		}
		my $v = $data{$_};
		my $vF = $v*cos($f*$x+$theta);
		my $vI = -1*$v*sin($f*$x+$theta);
		$F += $vF;
		$I += $vI;
	}
	my $steps = ($maxX-$min);

	$nP = $F*$F+$I*$I;

	$nP /= 2*$PI;
	my $Fm = $F;
	if ($F == 0) {
		$Fm = 1e-20;
	}
	my $phase = atan2($I,$Fm);
	$phase  = $phase/$PI*$p/2;
	$f = 1/$p;

	my $printP = $inc*floor($p/$inc+0.1);



	print "$printP\t$nP\t$phase\t$f\n";
}


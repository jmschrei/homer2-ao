#!/usr/bin/env perl
#

use POSIX;
my $I = 19;

my $upThresh = 0.58;
my $dnThresh = -1*$upThresh;
my $pseudo = 3.0;
my $logFlag = 0;
my $geneFlag= 0;

my @col1 = ();
my @col2 = ();

my $format = 'bed';
my $report = 'peak';
$format = 'peak';
my %genes = ();
$listFlag = 0;

sub printCMD {
	print STDERR "\n\tusage: calcScore.pl ... > output.tsv\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-i <input annotatePeaks.pl file> (required)\n";
	print STDERR "\t\t-1 #,#,... (column assignments for control condition, required)\n";
	print STDERR "\t\t-2 #,#,... (column assignments for treated condition, required)\n";
	print STDERR "\t\t\t-list (list column names/indexes)\n";
	print STDERR "\t\t-log (take log first, pseudo count: $pseudo)\n";
	print STDERR "\t\t\t-pseudo <#> (pseudo count to add before taking the log: $pseudo)\n";
	print STDERR "\t\t-gene (output gene symbols for nearest gene, taking most extreme regulation)\n";
	print STDERR "\t\t-thresh <#> (log2FC threshold to define up/down regulated genes)\n";
	print STDERR "\t\t-format <peak|gene> (def: $format)\n";
	print STDERR "\n";
	exit;
}


my $topCount = 25000;
my $inputFile = "";
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-i') {
		$inputFile = $ARGV[++$i];
		next;
	} elsif ($ARGV[$i] eq '-gene') {
		$report = 'gene';
	} elsif ($ARGV[$i] eq '-pseudo') {
		$pseudo = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-1') {
		@col1 = split /\,/, $ARGV[++$i];
		next;
	} elsif ($ARGV[$i] eq '-2') {
		@col2 = split /\,/, $ARGV[++$i];
		next;
	} elsif ($ARGV[$i] eq '-list') {
		$listFlag = 1;
	} elsif ($ARGV[$i] eq '-thresh') {
		$upThresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-format') {
		$format = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-log') {
		print STDERR "\tWill log data first (pseudo count=$pseudo)\n";
		$logFlag = 1;
	} elsif ($i == 0) {
		print STDERR "\tAssuming $ARGV[0] is the input file\n";
		$inputFile = $ARGV[0];
	} elsif ($i == 0 && $ARGV[$i] !~ /\-/) {
		$inputFile = $ARGV[$i];
		print STDERR "\tAssuming first argument is the input file ($inputFile)\n";
		next;
	} else {
		print STDERR "!!! \"$ARGV[$i]\" not recognized!\n";
		printCMD();
	}
}

if ($format eq 'peak') {
	$I = 19;
} elsif ($format eq 'gene') {
	$I = 8;
}

$dnThresh = -1*$upThresh;

if ($inputFile eq '') {
	print STDERR "!!! -i <inputfile> is required!!!\n";
	printCMD();
	exit;
}
if ($listFlag == 0 && (@col1 < 1 || @col2 < 1)) {
	print STDERR "!!! Must specify -1 <#,#,...> and -2 <#,#,...>\n";
	printCMD();
	exit;
}
print STDERR "\n\tfile: $inputFile\n";
print STDERR "\t\tCondition 1: @col1\n";
print STDERR "\t\tCondition 2: @col2\n";

open UP, ">up.txt";
open NC, ">nc.txt";
open DN, ">dn.txt";

my %data = ();
open IN, $inputFile or die "!!! Could not open \"$inputFile\"\n";
while (<IN>) {
	$c++;
	s/\r//g;
	my @line = split /\t/;
	if ($c < 2) {
		print STDERR "\tData columns:\n";
		for (my $i=0;$i+$I<scalar(@line);$i++) {	
			print STDERR "\t\t$i: $line[$i+$I]\n";
		}
		print STDERR "\n";
		if ($listFlag) {
			exit;
		}
		print STDERR "\tCondition1:\n";
		foreach(@col1) {
			my $x = $line[$I+$_];
			print STDERR "\t\t$x\n";
		}	
		print STDERR "\tCondition2:\n";
		foreach(@col2) {
			my $x = $line[$I+$_];
			print STDERR "\t\t$x\n";
		}	
		next;
	}
	my $p = $line[2];
	#my $p = floor(($line[2]+$line[3])/2);


	my $vv1 = 0;
	my $nn1 = 0;
	foreach(@col1) {
		my $x = $line[$I+$_];
		if ($logFlag) {
			$x = log($x+$pseudo)/log(2.0);
		}
		$vv1 += $x;
		$nn1++;
	}
	$vv1 /= $nn1;

	my $vv2 = 0;
	my $nn2 = 0;
	foreach(@col2) {
		my $x = $line[$I+$_];
		if ($logFlag) {
			$x = log($x+$pseudo)/log(2.0);
		}
		$vv2 += $x;
		$nn2++;
	}
	$vv2 /= $nn2;
	my $r = $vv2 - $vv1;

	#$r = ($line[$I+4]) - ($line[$I+0]);
	#$r = (($line[$I+2]+$line[$I+3])/2) - (($line[$I+0]+$line[$I+1])/2);
	#$r = (($line[$I+1]+$line[$I+2])/2) - ($line[$I+0]);
	#$r = (($line[$I+3]+$line[$I+4])/2) - (($line[$I+0]+$line[$I+1]+$line[$I+2])/3);
	#$r = (($line[$I+5]+$line[$I+4]+$line[$I+3])/3) - (($line[$I+0]+$line[$I+1]+$line[$I+2])/3);
	#$r = (($line[$I+2]+$line[$I+5]+$line[$I+8])/3) - (($line[$I+1]+$line[$I+4]+$line[$I+7])/3);
	#
	#
	#pseudo count (non rlog)
	#$r = log((($line[$I+1])+3)/($line[$I+0]+3))/log(2.0);
	#
	#$r = ($line[$I+10]) - ($line[$I+0]);
	#$r = log ((($line[$I+2])+3)/($line[$I+0]+3));
	#$r = (($line[$I+2]+$line[$I+3]+$line[$I+0]+$line[$I+1])/4);
	#$r = (($line[$I+1]+$line[$I+2])/2) - (($line[$I+0]+$line[$I+0])/2);
	#$r = (($line[$I+3]+$line[$I+3])/2) - (($line[$I+0]+$line[$I+0])/2);
	#$r = (($line[$I+5]+$line[$I+4]+$line[$I+3])/3) - (($line[$I+0]+$line[$I+1]+$line[$I+2])/3);
	#$r = (($line[$I+2]+$line[$I+2])/2) - (($line[$I+0]+$line[$I+0])/2);
	#$r = $line[$I+13] - $line[$I+0];
	#$r = (($line[$I+2]+$line[$I+5])/2) - (($line[$I+0]+$line[$I+3])/2);
	
	#NFY
	#$r = (($line[$I+3]+$line[$I+4])/2) - (($line[$I+0]+$line[$I+1]+$line[$I+2])/3);

	#focus ratio
	#my $den= ($line[$I+4]+$line[$I+5]/2);
	#next if ($den < 1e-10);
	#$r = ($line[$I+0]+$line[$I+1]/2) / ($line[$I+4]+$line[$I+5]/2);
	#
	my $s = $line[2];
	my $e = $line[3];
	my $id = $line[0];
	my $str = "$line[1]\t$s\t$e\t$line[0]\t$r\t$line[4]\n";
	if ($report eq 'peak') {
		$str = "$line[0]\t$line[1]\t$s\t$e\t$line[4]\t$r\n";
		$data{$id} = {s=>$str,v=>$r};
	} elsif ($report eq 'gene') {
		my $gene = $line[11];
		$str = "$gene\t$r\n";
		$id = $gene;
		if (exists($data{$id})) {
			if (abs($data{$id}->{'v'}) < abs($r)) {
				$data{$id} = {s=>$str,v=>$r};
			}
		} else {
			$data{$id} = {s=>$str,v=>$r};
		}
	} else {
		$data{$id} = {s=>$str,v=>$r};
	}

	if ($r >= $upThresh) {
		print UP $str;
		$numUP++;
	} elsif ($r < $dnThresh) {
		print DN $str;
		$numDN++;
	} else {
		print NC $str;
		$numNC++;
	}
}
close UP;
close DN;
close NC;
close IN;
print STDERR "\n\tNumber of regulated regions:\n";
print STDERR "\t\t$numUP\tup (>$upThresh)\n";
print STDERR "\t\t$numNC\tnc\n";
print STDERR "\t\t$numDN\tdn (<-$upThresh)\n";
print STDERR "\n";

my @ids = sort {$data{$a}->{'v'} <=> $data{$b}->{'v'}} keys %data;
my $totalTSS = scalar(@ids);

open TOP, ">top.txt";
open BOT, ">bottom.txt";

my $c = 0;
foreach(@ids) {
	$c++;
	my $id = $_;
	my $str = $data{$id}->{'s'};
	print $str;
	if ($c <= $topCount) {
		print BOT $str;
	}
	if ($c > $totalTSS-$topCount) {
		print TOP $str;
	}
}

close TOP;
close BOT;



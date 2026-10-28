#!/usr/bin/env perl
use warnings;


use POSIX;

sub printCMD() {

	print STDERR "\n\tUsage: removeRedundantPeaks.pl <peakfile> [options]\n";
	print STDERR "\n\tBy default, peaks with redundant positions [exact] will be removed\n";
	print STDERR "\tOutput is a new peak file sent to stdout. Additionally, peaks can be removed if they\n";
	print STDERR "\thave similar sequence composition (checked with BLAT).\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-palindromic (will remove peaks that overlap >50% and are on separate strands,\n";
	print STDERR "\t\t\t\tthat are likely representative of palindromic motif recognition)\n";
	print STDERR "\t\t-blat <#> (Where % is the percentage match to be considered redundant, default=0.33)\n";
	print STDERR "\t\t-size <#> (size of peaks to be used for sequence similarity, default=given)\n";
	print STDERR "\t\t-genome <directory of FASTA files> (genome for extracting sequence)\n";
	print STDERR "\t\t-mask (use repeat masked sequence)\n";
	print STDERR "\t\t-ignoreStrand\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 1) {
	printCMD();
}
my $peakFile = $ARGV[0];
my $minBLATscore = 0.33;
my $size = '';
my $genome = "";
my $maskFlag = "";
my $ignoreStrand = 0;
my $palindromicFlag = 0;
my $palindromicOverlap = 13;

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-blat') {
		$minBLATscore = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-mask') {
		$maskFlag = " -mask ";
	} elsif ($ARGV[$i] eq '-size') {
		$size = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-palindromic') {
		$palindromicFlag = 1;
	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-ignoreStrand') {
		$ignoreStrand = 1;
	} else {
		printCMD();
	}
}

my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";

`bed2pos.pl "$peakFile" -check > "$tmpFile"`;

my $initialPromoters = 0;
my %data = ();
my %data2 = ();
my $positionFiltered = 0;
my $idcounter = 0;
open IN, $tmpFile;
while (<IN>) {
	my $og = $_;
	chomp;
	s/\r//g;
	my @line  = split /\t/;

	next unless ($line[2] =~ /\d/);

	my $id = $line[0];
	my $chr = $line[1];
	my $mid = floor($line[2]+$line[3])/2;
	my $dir = '+';
	if ($line[4] eq '-' || $line[4] eq '1') {
		$dir = '-';
	}
	my $loc = $line[1] . '_' . $mid . '_' . $dir;
	if ($ignoreStrand) {
		$loc = $line[1] . '_' . $mid;
	}
	if (!exists($data{$loc})) {
		my %a = ();
		$data{$loc} = \%a;
	}
	my $v = 0;
	$v = $line[5] if (@line > 5);
	my $v2 = 0;
	$v2 = $line[6] if (@line > 6);

	$idcounter++;
	$id = "peakID-" . $idcounter;
	my $d = {og=>$og,s=>$line[2],e=>$line[3],d=>$dir,v=>$v,v2=>$v2,loc=>$loc};
	$data{$loc}->{$id}=$d;
	if (!exists($data2{$chr})) {
		my %a =();
		$data2{$chr}=\%a;
	}
	$data2{$chr}->{$id} = $d;
	$initialPromoters++;
}
close IN;
`rm "$tmpFile"`;

my %newData = ();
my %pos = ();
if ($palindromicFlag) {
	print STDERR "\tRemoving overlapping/palindromic matches. (palindromic overlap=$palindromicOverlap)\n";
	foreach(keys %data2) {	
		my $chr = $_;
		my @peaks = sort {$data2{$chr}->{$a}->{'s'} <=> $data2{$chr}->{$b}->{'s'}} keys %{$data2{$chr}};

		for (my $i=0;$i<@peaks;$i++) {
			my $good = 1;
			my $id1 = $peaks[$i];
			my $s1 = $data2{$chr}->{$id1}->{'s'};
			my $e1 = $data2{$chr}->{$id1}->{'e'};
			my $d1 = $data2{$chr}->{$id1}->{'d'};
			my $v1 = $data2{$chr}->{$id1}->{'v2'};
			my $mid1 = floor(($s1+$e1)/2);
			for (my $j=$i+1;$j<@peaks;$j++) {
				my $id2 = $peaks[$j];
				my $s2 = $data2{$chr}->{$id2}->{'s'};
				my $e2 = $data2{$chr}->{$id2}->{'e'};
				my $d2 = $data2{$chr}->{$id2}->{'d'};
				my $v2 = $data2{$chr}->{$id2}->{'v2'};
				my $mid2 = floor(($s2+$e2)/2);

				last if ($s2 > $e1+$palindromicOverlap);
				if (abs($mid1-$mid2) <= $palindromicOverlap) {
					if ($v2 >= $v1) {
						#print STDERR "bad+: $mid1\t$mid2\n";
						$good = 0;
						last;
					}
				}
			}
			for (my $j=$i-1;$j>=0;$j--) {
				my $id2 = $peaks[$j];
				my $s2 = $data2{$chr}->{$id2}->{'s'};
				my $e2 = $data2{$chr}->{$id2}->{'e'};
				my $d2 = $data2{$chr}->{$id2}->{'d'};
				my $v2 = $data2{$chr}->{$id2}->{'v2'};
				my $mid2 = floor(($s2+$e2)/2);
				last if ($e2 < $s1-$palindromicOverlap);
				if (abs($mid1-$mid2) <= $palindromicOverlap) {
					if ($v2 > $v1) {
						#print STDERR "bad-: $mid1\t$mid2\n";
						$good = 0;
						last;
					}
				}
			}
			if ($good) {
				my $best = $id1;
				my $og = $data2{$chr}->{$best}->{'og'};
				my $loc = $data2{$chr}->{$best}->{'loc'};
				$newData{$loc} = $og;
				$pos{$best} = $og;
				$positionFiltered++;
			}
		}
	}	
} else {
	print STDERR "\tOnly removing exact peak matches.\n";
	foreach(keys %data) {	
		my $loc = $_;
		my @a = sort {$data{$loc}->{$b}->{'v2'} <=> $data{$loc}->{$a}->{'v2'} || $a cmp $b} keys %{$data{$loc}};
		my $best = $a[0];
		my $og = $data{$loc}->{$best}->{'og'};
		$newData{$loc} = $og;
		$pos{$best} = $og;
		$positionFiltered++;
	}	
}
print STDERR "\n\tKept $positionFiltered of $initialPromoters based on genome position\n";

if ($genome eq '') {
	foreach(values %newData) {
		print $_;
	}
	exit;
}

open OUT, ">$tmpFile";
foreach(values %newData) {
	print OUT $_;
}
close OUT;
if ($size ne '') {
	print STDERR "\tResizing...\n";
	`adjustPeakFile.pl "$tmpFile" -size $size > "$tmpFile2"`;
	#`resizePosFile.pl "$tmpFile" $size > "$tmpFile2"`;
	`mv "$tmpFile2" "$tmpFile"`;
}
`makeBinaryFile.pl "$tmpFile" "$tmpFile" > "$tmpFile3"`;

`homerTools extract "$tmpFile" "$genome" $maskFlag > "$tmpFile2"`;
`findRedundantBLAT.pl "$tmpFile2" $minBLATscore > "$tmpFile"`;
`adjustRedunGroupFile.pl "$tmpFile3" "$tmpFile" > "$tmpFile2"`;
open IN, $tmpFile2;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if (exists($pos{$line[0]})) {
		print "$pos{$line[0]}";
	}
}
close IN;

`rm "$tmpFile" "$tmpFile2" "$tmpFile3"`;

#!/usr/bin/env perl
#use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";

use Statistics;

#
if (@ARGV < 2) {
	print STDERR "<input bedGraph> <match bedGraph> -> output to stdout\n";
	print STDERR "\n";
	exit;
}

my %factor = ();
my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my %data = ();
open OUT, ">$tmpFile";
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	if (/^track/) {
		print "$_\n";
		next;
	}
	my @line = split /\t/;
	my $id = "$line[0]:$line[1]-$line[2]";
	my $chr = $line[0];
	my $start = $line[1];
	my $end = $line[2];
	my $v = $line[3];
	my $p = {id=>$id,c=>$chr,s=>$start,e=>$end,v=>$v};
	push(@data, $p);
	print OUT "$id\t$chr\t$start\t$end\t+\t$v\n";
}
close IN;
close OUT;

`annotatePeaks.pl $tmpFile none -noblanks -bedGraph $ARGV[0] $ARGV[1] > $tmpFile2`;

my %data1 = ();
my %data2 = ();
my $c = 0;
open IN, $tmpFile2;
while (<IN>) {
	$c++;
	if ($c < 2) {
		next;
	}
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $chr= $line[1];
	my $v1 = $line[$numCols-2];
	my $v2 = $line[$numCols-1];
	if (!exists($data1{$chr})) {
		my @a = ();
		$data1{$chr}=\@a;
		my @b = ();
		$data2{$chr}=\@b;
	}
	push(@{$data1{$chr}},$v1);
	push(@{$data2{$chr}},$v2);
}
close IN;

foreach(keys %data1) {
	my $chr = $_;
	my ($r,$p) = Statistics::correlation($data1{$chr},$data2{$chr});
	print STDERR "\t$chr\t$p\t$r\n";
	if ($r > 0) {
		$factor{$chr}=1;
	} else {
		$factor{$chr}=-1;
	}
}
foreach(@data) {
	my $chr = $_->{'c'};
	next if (!exists($factor{$chr}));
	my $v = $_->{'v'} * $factor{$chr};
	print "$chr\t$_->{'s'}\t$_->{'e'}\t$v\n";
}
`rm $tmpFile $tmpFile2`;

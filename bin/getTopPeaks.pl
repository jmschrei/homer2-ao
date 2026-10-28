#!/usr/bin/env perl
use warnings;

if (@ARGV < 2) {
	print STDERR "<peak file> <# of peaks> [# header]\n";
	print STDERR "outputs to stdout - uses 6th column to find best peaks\n";
	exit;
}
my $header = 0;
if (@ARGV > 2) {
	$header = $ARGV[2];
}
my %peaks = ();
my $z= 0;
open IN, $ARGV[0];
while (<IN>) {
	$z++;
	chomp;
	s/\r//g;
	my $og = $_;
	if ($z <= $header) {
		print $og . "\n";
		next;
	}
	my @line = split /\t/;
	next if ($line[0] =~ /^#/);
	next if (@line < 6);
	next unless ($line[5] =~ /^\-*[\d\.e\+\-]+$/);
	$peaks{$line[0]} = {s=>$line[5], og=>$og};
	
}
close IN;

my @peaks = sort {$peaks{$b}->{'s'} <=> $peaks{$a}->{'s'}} keys %peaks;
for (my $i=0;$i<@peaks;$i++) {
	last if ($i >= $ARGV[1]);
	print "$peaks{$peaks[$i]}->{'og'}\n";
}

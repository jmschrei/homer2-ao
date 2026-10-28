#!/usr/bin/env perl
use warnings;


if (@ARGV < 2) {
	print STDERR "<data> <list of ids to remove> <# lines of header> <set to 1 if want to keep IDs> [1=ignore acc version]\n";
	exit;
}

my $header = 0;
if (@ARGV > 2) {
	$header = $ARGV[2];
}
my $oppFlag = 0;
if (@ARGV > 3) {
	$oppFlag = $ARGV[3];
}

my $accFlag = 0;
if (@ARGV > 4) {
	$accFlag = 1;
}
#print STDERR "\tFiltering:\n";
#print STDERR "\t\theaderLines set to $header\n";
#print STDERR "\t\tWill keep list instead of filtering: $oppFlag\n";
#print STDERR "\t\tWill ignore differences between accession number versions\n" if ($accFlag);

my %list = ();
open IN, $ARGV[1];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $id = shift @line;
	if ($accFlag) {
		$id =~ s/\.\d+$//;
	}
	$list{$id} = 1;
}
close IN;

my $c= 0;
open IN, $ARGV[0];
while (<IN>) {
	$c++;
	chomp;
	s/\r//g;
	if ($c <= $header) {
		print "$_\n";
		next;
	}
	my @line = split /\t/;
	my $id = shift @line;
	my $id2 = $id;
	if ($accFlag) {
		$id2 =~ s/\.\d+$//;
	}
	my $pflag = 0;
	if ($oppFlag) {
		if (exists($list{$id}) || exists($list{$id2})) {
			$pflag = 1;
		}
	} else {
		if (!exists($list{$id}) && !exists($list{$id2})) {
			$pflag = 1;
		}
	}
	if ($pflag) {
		print "$id";
		foreach(@line) {
			print "\t$_";
		}
		print "\n";
	}

}
close IN;
	




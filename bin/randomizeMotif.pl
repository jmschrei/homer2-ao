#!/usr/bin/perl -w

my $numMotifs = 1;
my $order=0;
my $originalFlag = 0;
my $prefix = "Rand";

sub printCMD {
	print STDERR "\n\tusage: randomizeMotif.pl <motif file> [options]\n";
	print STDERR "\t\tThis program will output motif files where the positions are randomly swapped\n";
	print STDERR "\t\tOutputs a new motif file (may contain multiple randomized motifs\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-n <#> (number of randomized motifs to produce, def: $numMotifs)\n";
	print STDERR "\t\t-order <#> (0=random swap, 1=tries to maintain dinuc, 2=trinuc, etc., def: $order)\n";
	print STDERR "\t\t-original (include original motif in output file)\n";
	print STDERR "\t\t-name <prefix> (New motifs will be named prefix_1_motifname, ... def: $prefix)\n";
	print STDERR "\n";
	exit;
}
if (@ARGV < 1) {
	printCMD();
}
my $motifFile = $ARGV[0];

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-n') {
		$numMotifs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-original') {
		$originalFlag = 1;
	} elsif ($ARGV[$i] eq '-order') {
		$order = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-prefix') {
		$prefix = $ARGV[++$i];
	} else {
		print STDERR "!!! Error: Did not recognize option: \"$ARGV[$i]\"\n";
		printCMD();
	}
}
unless (-e $motifFile) {
	print STDERR "!!! Error: Motif file \"$motifFile\" doesn't exist!\n";
	exit;
}
my @rows = ();
open IN, $motifFile;
my $z= 0;
while (<IN>) {
	$z++;
	my $og = $_;
	if ($originalFlag) {
		print $og;
	}
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($z==1) {
		@header = @line;
		@consensus = split "",$header[0];
		shift @consensus;
		next;
	} else {
		push(@rows, $og);
	}
}
close IN;

my $len = scalar(@rows);


for (my $z=0;$z<$numMotifs;$z++) {
	my $name = $prefix . "_" . $z . "_" . $header[1];

	my @set = ();
	my %used = ();
	for (my $i=0;$i<@rows;$i++) {
		my $x = {r=>rand(),i=>$i};
		$used{$i}=0;
		push(@set, $x);
	}
	@set = sort {$a->{'r'} <=> $b->{'r'}} @set;
	print ">random\t$name";
	for (my $i=2;$i<@header;$i++) {
		print "\t$header[$i]";
	}
	print "\n";
	
	foreach(@set) {
		my $index = $_->{'i'};
		next if ($used{$index});
		print $rows[$index];
		$used{$index}=1;
		for (my $i=0;$i<$order;$i++) {
			$index++;
			last if ($index >= $len);
			last if ($used{$index});
			print $rows[$index];
			$used{$index}=1;
		}
	}
}

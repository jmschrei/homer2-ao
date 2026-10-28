#!/usr/bin/perl -w
#

if (@ARGV < 2) {
	print STDERR "<RBP_id_encoding.csv> <dir with *.PWM files>\n";
	exit;
}

my %data = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\"//g;
	my @line = split /\, /;
	$line[1] =~ s/\ /\_/g;
	$line[1] =~ s/\(//g;
	$line[1] =~ s/\)//g;
	if ($line[0] < 10) {
		$line[0] = '00' . $line[0];
	} elsif ($line[0] < 100) {
		$line[0] = '0' . $line[0];
	}
	print STDERR "$line[0]\t$line[1]\n";
	open IN2, "$ARGV[1]/$line[0].PWM" or die "Could not open $ARGV[1]/$line[0].PWM\n";
	print ">$line[1]\t$line[1]/$line[1]-$line[0]/oRNAament\t0\n";
	my $z = 0;
	while (<IN2>) {
		$z++;
		chomp;
		my @line = split /\t/;
		next if ($z < 2);
		for (my $i=1;$i<@line;$i++){ 
			print "\t" if ($i > 1);
			print "$line[$i]";
		}
		print "\n";
	}
	close IN2;
}
close IN;



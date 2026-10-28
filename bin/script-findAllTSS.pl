#!/usr/bin/env perl
use warnings;

my $commonOpt = "";
#$commonOpt = " -ntagThreshold ";
#my $commonOpt = "-metagene";

if (@ARGV < 1) {
	print STDERR "<in.tss table> <gtf> <genome>\n";
	print STDERR "in.tss format:\n";
	print STDERR "<prefixname>TAB<csRNAdir>TAB<inputDir>TAB<rnaDir>TAB<gtf>TAB<genome>TAB<options>\n";
	exit;
}
my $gtf = "~/bennerlab/species/genomes/homo_sapiens/genes.gtf";
my $genome = "~/bennerlab/species/genomes/homo_sapiens/genome.fa";
if (@ARGV > 1) {
	$gtf = $ARGV[1];
	if ($gtf eq 'none') {
		$gtf = '';
	}
}
if (@ARGV > 2) {
	$genome = $ARGV[2];
	if ($genome eq 'none') {
		$genome = '';
	}
}
if (@ARGV > 3) {
	for (my $i=3;$i<@ARGV;$i++) {
		$commonOpt .= " " . $ARGV[$i];
	}
	print STDERR "\tWill use following options: $commonOpt\n";
}

my $numCPUs = 3;

print STDERR "\tGTF: $gtf\n";
print STDERR "\tGenome: $genome\n\n";

my @sets = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	next if ($_ eq '');
	my @line = split /\t/;
	next if ($line[0] eq '' || $line[0] =~ /^#/);
	my $name = $line[0];
	my $tagdir = '';
	my $input = '';
	my $rna = '';
	my $gtfX = $gtf;
	my $genomeX = $genome;
	my $opt = '';
	$tagdir = $line[1] if (@line > 1);
	$input = $line[2] if (@line > 2);
	$rna = $line[3] if (@line > 3);
	$gtfX = $line[4] if (@line > 4);
	$gtfX = $gtf if ($gtfX eq '');
	$genomeX = $line[5] if (@line > 5);
	$genomeX = $genome if ($genomeX eq '');
	$opt = $line[6] if (@line > 6);
	next if ($tagdir eq '' || $name eq '');

	my $t = {output=>$name, tagdir=>$tagdir, input=>$input, rna=>$rna,gtf=>$gtfX,genome=>$genomeX,opt=>$opt};
	push(@sets, $t);
	print STDERR "\t$name(name)\t$tagdir(tagdir)\t$input(input)\t$rna(RNA)\t$gtfX(gtf)\t$genomeX(genome)\t$opt\n";
}
close IN;

print STDERR "\tsleeping for 3 seconds...\n";
`sleep 3`;
my @cmds = ();
foreach(@sets) {
	my $output = $_->{'output'};
	my $tagdir = $_->{'tagdir'};
	my $input = $_->{'input'};
	my $rna = $_->{'rna'};
	my $gtf = $_->{'gtf'};
	my $genome = $_->{'genome'};
	my $opt = $_->{'opt'} . " $commonOpt ";
	my $cmd = "findcsRNATSS.pl $tagdir -o $output";
	if ($input ne '') {
		$cmd .= " -i $input ";
	}
	if ($rna ne '') {
		$cmd .= " -rna $rna ";
	}
	if ($gtf ne '') {
		$cmd .= " -gtf $gtf ";
	}
	$cmd .= " -genome $genome ";
	$cmd .= " -cpu $numCPUs ";
	$cmd .= " " . $opt;

	print "$cmd\n";
	push(@cmd, $cmd);
}

foreach(@cmd) {
	`sleep 1`;
	my $pid = fork();
	if ($pid==0) {
		`$_`;
		exit;
	}
}
my $id = 0;
while ($id>=0) {
	$id=wait();
}
print STDERR "\tAll FINISHED\n";

#!/usr/bin/env perl
use warnings;



if (@ARGV < 1) {
	print STDERR "\n\tbatchMakeTagDirectory.pl <key file> [makeTagDirectory options]...\n";
	print STDERR "\n\tKey File: Tag delimited: <directory>TAB<alignment file>\n";
	print STDERR "\t\tDuplicate directory name for multiple alignment files for same experiment\n";
	print STDERR "\t\tOption: -cpu <#> : to run more than one at a time...\n";
	print STDERR "\t\tOption: -r : TagDir names with \"-r1\", \"-r2\", etc. will also be merged\n";
	print STDERR "\t\t\t-ronly (only create merged directories, not individual replicate directories)\n";
	print STDERR "\n";
	exit;
}
my $maxCPUs = 1;
my $mergeReps = 0;
my $otherOptions = "";

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-r') {
		$mergeReps = 1;
		next;
	}
	if ($ARGV[$i] eq '-ronly') {
		$mergeReps = 2;
		next;
	}
	$otherOptions .= " $ARGV[$i]";
}


my %key = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line = split /\t/;
	if ($mergeReps != 2) {
		if (!exists($key{$line[0]})) {
			my @a = ();
			$key{$line[0]}= \@a;
		}
		push(@{$key{$line[0]}}, $line[1]);
	}
	if ($mergeReps) {
		my $mergeName = $line[0];
		$mergeName =~ s/\-r\d+\-.+/\-combo/;
		$mergeName =~ s/\-r\d+$/\-combo/;
		$mergeName =~ s/\_r\d+\_.+/\-combo/;
		$mergeName =~ s/\_r\d+$/\-combo/;
		if ($mergeName eq $line[0]) {
			print STDERR "\t! warning, $line[0] doesn't match -r \n";
			next;
		}
		print STDERR "\tmerging $line[0] -> $mergeName\n";
		if (!exists($key{$mergeName})) {
			my @a = ();
			$key{$mergeName}= \@a;
		}
		push(@{$key{$mergeName}}, $line[1]);
	}
}
close IN; 


print STDERR "\n\tWill create the following Tag Directories:\n";
my @keys = sort {$a cmp $b} keys %key;
foreach(@keys) {
	my $td = $_;
	print STDERR "\t$td\n";
	foreach(@{$key{$td}}) {
		print STDERR "\t\t$_\n";
	}
}
print STDERR "\n\tWaiting 3 seconds before starting...\n";
`sleep 3`;

my @pids = ();
my $cpus = 0;
foreach(keys %key) {
	my $dirName = $_;
	my $fileList = "";
	foreach(@{$key{$dirName}}) {
		$fileList .= " " . $_ . " ";
		#$fileList .= " \"" . $_ . "\"";
	}
	print STDERR "`makeTagDirectory $dirName $otherOptions $fileList`\n";

	my $pid = fork();
	$cpus++;

	if ($pid == 0) {
		`makeTagDirectory "$dirName" $otherOptions $fileList `;
		exit(0);
	}

	push(@pids, $pid);
	if ($cpus >= $maxCPUs) {
		my $id = wait();
		$cpus--;
	}
}
my $id = 0;
while ($id >= 0) {
	$id= wait();
}
	


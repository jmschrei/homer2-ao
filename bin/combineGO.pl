#!/usr/bin/env perl


sub printCMD {
	print STDERR "\n\tcombineGO.pl [options] -f <filename i.e. biological_process.txt> -d <directory1> [directory2] ...\n";
	print STDERR "\n\toptions:\n";
	print STDERR "\t\t-f <filename> (Filename to join enrichment values from)\n";
	print STDERR "\t\t\t-m (equivalent to -f knownResults.txt)\n";
	print STDERR "\t\t-d <directory1> [directory2] ... (findMotifs.pl output directories to join)\n";
	print STDERR "\t\t-top <#> (Only keep top # terms per directory, default: keep all)\n";
	print STDERR "\t\t-column <#> (Column in files to join, default: 4)\n";
	print STDERR "\t\t-ratio (return log2 ratio target/background if using motifs/knownResults.txt)\n";
	print STDERR "\t\t\t-pseudoFraction <#> (Default: 2%)\n";
	print STDERR "\n";
	exit;
}
if (@ARGV < 3) {
	printCMD();
}

#my $col = 10;
my $col = 4;
my $top = 0;
my $filename= $ARGV[0];
my @tagDirs = ();
my $ratioFlag = 0;
my $pseudoFraction = 0.02;

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-f') {
		$filename=$ARGV[++$i];
	} elsif ($ARGV[$i] eq '-m') {
		$filename= "knownResults.txt";
	} elsif ($ARGV[$i] eq '-d') {
        my $bail = 0;
        while ($ARGV[++$i] !~ /^\-/) {
            push(@tagDirs, $ARGV[$i]);
            print STDERR "\t\t$ARGV[$i]\n";
            if ($i>=@ARGV-1) {
                $bail = 1;
                last;
            }
        }
        last if ($bail == 1);
        $i--;
	} elsif ($ARGV[$i] eq '-top') {
		$top = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-column') {
		$col = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-ratio') {
		$ratioFlag = 1;
	} elsif ($ARGV[$i] eq '-pseudoFraction') {
		$pseudoFraction = $ARGV[++$i];
	} else {
		print STDERR "!! Can't recognize $ARGV[$i]\n";
		printCMD();
	}
}


my %top = ();

my $targetTotal = 1;
my $bgTotal = 1;

print "TermID\tTermDesc";
for (my $i=0;$i<@tagDirs;$i++) {
	my $dir = $tagDirs[$i];
	open IN, "$dir/$filename" or next;
	#open IN, "$dir/$filename" or die "Could not open $dir/$filename!!!\n";

	print "\t$dir";

	my $c=0;
	my %d = ();
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($c==1) {
			if ($filename eq 'knownResults.txt') {
				if (@line < 8) {
					print STDERR "!!! Something is wrong - knownResults.txt doesn't have enough columns!\n";
					exit;
				}
				$line[5] =~ /\(of (\d+)\)/;
				$targetTotal = $1;
				$line[7] =~ /\(of (\d+)\)/;
				$bgTotal = $1;
				if ($targetTotal < 2 || $bgTotal < 2) {
					print STDERR "!!!Problem figuring out target and bg totals\n";
					#exit;
				}
			}
			next;
		} 
		my $id = $line[0] . "\t" . $line[1];
		my $v = $line[$col-1];

		if ($ratioFlag && $filename eq 'knownResults.txt') {
			my $target = $line[5];
			my $bg  = $line[7];
			$bg = 1 if ($bg < 1);
			$target = 1 if ($target < 1);
			my $n = $target/$targetTotal;
			my $d = $bg/$bgTotal;
			$v = log(($n+$pseudoFraction)/($d+$pseudoFraction))/log(2.0);
				
		}
		$d{$id} = $v;
	}
	close IN;

	my @a = sort {$d{$a} <=> $d{$b}} keys %d;
	if ($top > 0) {
		for (my $j=0;$j<$top;$j++) {
			$top{$a[$j]}=1;
		}
	} else {
		foreach(keys %d) {
			$top{$_}=1;
		}
	}
	push(@data, \%d);
}
print "\n";
foreach(keys %top) {
	my $id = $_;
	my $str = $_;
	if ($filename eq 'knownResults.txt') {
		$str =~ s/\/.*\/Homer//;
	}
	print "$str";
	foreach(@data) {
		my $v = 0;
		if (exists($_->{$id})) {
			$v = $_->{$id};
		}
		print "\t$v";
	}
	print "\n";
}


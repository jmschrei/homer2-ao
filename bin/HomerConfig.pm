#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";

package HomerConfig;



sub loadConfigFile {

	my ($file) = @_;
	if (!defined($file)) {
		$file = $homeDir . "/config.txt";
	}
	my %a = ();
	my %b = ();
	my %c = ();
	my %d = ();
	my %e = ();
	my $config = {PROMOTERS=>\%a, GENOMES=>\%b, ORGANISMS=>\%c,SETTINGS=>\%d,SOFTWARE=>\%e};

	parseConfigFile($file,$config,$homeDir);

	my $localConfig = $ENV{"HOME"} . "/.homerConfig.txt";
	if (-e $localConfig) {
		#parseConfigFile($localConfig,$config,"");
	}
	return $config;
}
sub printConfigFile {
	my ($config,$outFile, $updateFlag) = @_;

	my $file = $homeDir . "/config.txt";
	if (defined($outFile)) {
		$file = $outFile;
		#print STDERR "\tPrinting Configuration to $outFile (may not be standard...)\n";
	}

	open CONFIG, ">$file";

	if (defined($updateFlag)) {
		print CONFIG "# Homer Update File\n";
	} else {
		print CONFIG "# Homer Configuration File (automatically generated)\n";
	}
	print CONFIG "#\n";
	print CONFIG "# This file is updated from the Homer website and contains information about data available for\n";
	print CONFIG "# use with the program.\n";
	print CONFIG "#\n";
	print CONFIG "# This file has been further modified by the local update script\n";
	print CONFIG "#\n";
	print CONFIG "# Each section has the same format, which is <tab> separated values specifying:\n";
	print CONFIG "# package name <tab> version <tab> description <tab> url <tab> optional parameters (, separated)\n";
	print CONFIG "#\n";

	my @groups = ("SOFTWARE","ORGANISMS","PROMOTERS","GENOMES");
	foreach(@groups) {
		my $group = $_;
		print CONFIG "\n$group\n";
		foreach(keys %{$config->{$group}}) {
			my $package = $_;
			my $version = $config->{$group}->{$package}->{'version'};
			my $desc = $config->{$group}->{$package}->{'desc'};
			my $url = $config->{$group}->{$package}->{'url'};
			my $location = $config->{$group}->{$package}->{'location'};
			my $params = "";
			my $c = 0;
			foreach(@{$config->{$group}->{$package}->{'parameters'}}) {
				$params .= "," if ($c > 0);
				$c++;
				$params .= $_;
			}
			print CONFIG "$package\t$version\t$desc\t$url\t$location\t$params\n";
		}
	}
	if (defined($updateFlag)) {
	} else {
		print CONFIG "SETTINGS\n";
		if (exists($config->{'SETTINGS'})) {
			foreach(keys %{$config->{'SETTINGS'}}) {
				my $var = $_;
				my $val = $config->{'SETTINGS'}->{$var};
				print CONFIG "$var=$val\n";
			}
		}
	}
	close CONFIG;
}
sub parseConfigFile {
	my ($file, $config,$homeDir) = @_;
	my $mode = '';	
	open IN, $file or die "Could not open configuration file ($file)\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		s/#.*$//;
		s/^\s+//;
		next if ($_ eq '');

		if (/^SOFTWARE/) {
			$mode = 'SOFTWARE';
			next;
		}
		if (/^ORGANISMS/) {
			$mode = 'ORGANISMS';
			next;
		}
		if (/^PROMOTERS/) {
			$mode = 'PROMOTERS';
			next;
		}
		if (/^GENOMES/) {
			$mode = 'GENOMES';
			next;
		}
		if (/^SETTINGS/) {
			$mode = 'SETTINGS';
			next;
		}
		#if (/^PREPARSED/) {
		#	$mode = 'PREPARSED';
		#	next;
		#}
		my @line = split /\t/;
		next if ($mode eq '');

		if (@line > 4) {
			my $package = $line[0];
			my $version = $line[1];
			my $description = $line[2];
			my $url = $line[3];
			my $location = $line[4];
			my @params = ();
			if (@line > 5) {
				@params = split /\,/, $line[5];
			}
			my $d = $homeDir . "/" . $location;
			if ($mode eq 'GENOMES') {
				$config->{'GENOMES'}->{$package} = {org=>$params[0],promoters=>$params[1], directory=>$d,
							location=>$location, version=>$version,url=>$url,desc=>$description, parameters=>\@params};
			} elsif ($mode eq 'PROMOTERS') {
				$config->{'PROMOTERS'}->{$package} = {org=>$params[0], directory=>$d, genome=>$params[1], 
							idtype=>$params[2], start=>$params[3], end=>$params[4],
							location=>$location, version=>$version,url=>$url,desc=>$description, parameters=>\@params};
			} elsif ($mode eq 'ORGANISMS') {
				my $taxid = 'NA';
				my $source = 'NA';
				$taxid = $params[0] if (@params > 0);
				$source = $params[1] if (@params > 1);
				$config->{'ORGANISMS'}->{$package} = {taxid=>$taxid,source=>$source, directory=>$d,
							location=>$location, version=>$version,url=>$url,desc=>$description, parameters=>\@params};
			} elsif ($mode eq 'SOFTWARE') {
				$config->{'SOFTWARE'}->{$package} = {
							location=>$location, version=>$version,url=>$url,desc=>$description, parameters=>\@params};
			#} elsif ($mode eq 'PREPARSED') {
			#	$config->{'PREPARSED'}->{$package} = {genome=>$params[0],size=>$params[1], directory=>$d,
			#				location=>$location, version=>$version,url=>$url,desc=>$description, parameters=>\@params};
			}
		} elsif ($mode eq 'SETTINGS') {
			my @var = split /\=/, $line[0];
			next if (@var < 2);
			$config->{'SETTINGS'}->{$var[0]} = $var[1];
		}
	}
	return $config;
}

sub readMotifFile {
	my ($file) = @_;
	my @motifs = ();
	my $index = -1;
	open MOTIFFILE, $file;
	while (<MOTIFFILE>) {
		my $og = $_;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($line[0] =~ /^>/) {
			my @a = ();
			my $logp = 0;
			if (@line > 3) {
				$logp = $line[3];
			}
			my $m = {name=>$line[1],consensus=>$line[0],len=>0,matrix=>\@a,text=>$og,threshold=>$line[2],logp=>$logp};
			push(@motifs, $m);
			$index++;
		} else {
			$motifs[$index]->{'len'}++;
			$motifs[$index]->{'text'} .= $og;
			push(@{$motifs[$index]->{'matrix'}}, \@line);
		}
	}
	close MOTIFFILE;
	return \@motifs;
}



# parameters:
# $motif - homer motif hash pointer
# $seq - sequence to scan
# Will not score motifs that overlap a 'N', will report 'na'
# return 2 arrays for + and - matching scores, of length sequence - motif length + 1
sub scoreMotifAtAllPositions {
	my ($motif, $seq) = @_;
	if (!defined($motif->{'lom'})) {
		addLogOddsMatrix($motif);
	}
	my $L = length($seq);
	my $mlen = $motif->{'len'};
	$seq =~ s/A/0/g;
	$seq =~ s/C/1/g;
	$seq =~ s/G/2/g;
	$seq =~ s/T/3/g;
	$seq =~ s/N/4/g;
	my @seq = split(//,$seq);
	my @scoresPos = ();
	my @scoresNeg = ();
	for (my $i=0;$i<$L-$mlen+1;$i++) {
		my $fscore = 0;
		my $rscore = 0;
		my $Nflag = 0;
		for (my $j=0;$j<$mlen;$j++) {
			my $x = $seq[$i+$j];
			if ($x > 3)  {
				$Nflag=1;
				last;
			}
			$fscore += $motif->{'lom'}->[$j][$x];
			$rscore += $motif->{'lom'}->[$mlen-1-$j][3-$x];
		}
		if ($Nflag) {
			$fscore = 'na';
			$rscore = 'na';
		}
		push(@scoresPos, $fscore);
		push(@scoresNeg, $rscore);
	}
	return (\@scoresPos, \@scoresNeg);
}

# parameters:
# $motif - homer motif hash pointer
# $seq - sequence to scan
# $strand - if "+" or "-" strand specific, otherwise both strands
# If no motifs are found, $bestPos will be negative
# Will not score motifs that overlap a 'N'
sub scoreMotifInSequence {
	my ($motif, $seq, $strand) = @_;
	if (!defined($motif->{'lom'})) {
		addLogOddsMatrix($motif);
	}
	my $L = length($seq);
	my $mlen = $motif->{'len'};
	$seq =~ s/A/0/g;
	$seq =~ s/C/1/g;
	$seq =~ s/G/2/g;
	$seq =~ s/T/3/g;
	$seq =~ s/N/4/g;
	my @seq = split(//,$seq);
	my $bestScore = -1e9;
	my $bestPos = -1;
	my $bestStrand = '.';
	my @scores = ();
	my @positions = ();
	my @strands = ();
	for (my $i=0;$i<$L-$mlen+1;$i++) {
		my $fscore = 0;
		my $rscore = 0;
		my $Nflag = 0;
		for (my $j=0;$j<$mlen;$j++) {
			my $x = $seq[$i+$j];
			if ($x > 3)  {
				$Nflag=1;
				last;
			}
			$fscore += $motif->{'lom'}->[$j][$x];
			$rscore += $motif->{'lom'}->[$mlen-1-$j][3-$x];
		}
		next if ($Nflag);
		if ($strand ne '-') {
			if ($fscore > $bestScore) {
				$bestScore = $fscore;
				$bestPos = $i;
				$bestStrand = "+";
			}
			if ($fscore >= $motif->{'threshold'}) {
				push(@scores, $fscore);
				push(@positions, $i);
				push(@strands, "+");
			}
		}
		if ($strand ne '+') {
			if ($rscore > $bestScore) {
				$bestScore = $rscore;
				$bestPos = $i;
				$bestStrand = '-';
			}
			if ($rscore >= $motif->{'threshold'}) {
				push(@scores, $rscore);
				push(@positions, $i);
				push(@strands, "-");
			}
		}
	}

	return ($bestScore, $bestPos, $bestStrand, \@scores, \@positions, \@strands);
}
sub addLogOddsMatrix {
	my ($motif) = @_;
	my $mlen = $motif->{'len'};
	my @m = ();
	for (my $i=0;$i<$mlen;$i++) {
		my @lodds = ();
		for (my $j=0;$j<4;$j++) {
			my $v = $motif->{'matrix'}->[$i][$j];
			$v = log($v/0.25);
			push(@lodds,$v);
		}
		push(@m, \@lodds);
	}
	$motif->{'lom'} = \@m;
}
sub getPossibleScoreChanges {
	my ($motif, $nuc1, $nuc2, $strand) = @_;
	my @scores = ();
	if (!defined($motif->{'lom'})) {
		addLogOddsMatrix($motif);
	}
	my $mlen = $motif->{'len'};
	my $i1 = 0;
	my $i2 = 0;
	if ($nuc1 eq 'A') {
		$i1 = 0;
	} elsif ($nuc1 eq 'C') {
		$i1 = 1;
	} elsif ($nuc1 eq 'G') {
		$i1 = 2;
	} elsif ($nuc1 eq 'T') {
		$i1 = 3;
	} else {
		return \@scores;
	}
	if ($nuc2 eq 'A') {
		$i2 = 0;
	} elsif ($nuc2 eq 'C') {
		$i2 = 1;
	} elsif ($nuc2 eq 'G') {
		$i2 = 2;
	} elsif ($nuc2 eq 'T') {
		$i2 = 3;
	} else {
		return \@scores;
	}
	for (my $i=0;$i<$mlen;$i++) {
		if ($strand ne '-') {
			my $fdiff = $motif->{'lom'}->[$i][$i2] - $motif->{'lom'}->[$i][$i1];
			push(@scores, $fdiff);
		}
		if ($strand ne '+') {
			my $rdiff = $motif->{'lom'}->[$mlen-1-$i][3-$i2] - $motif->{'lom'}->[$mlen-1-$i][3-$i1];
			push(@scores, $rdiff);
		}
	}
	return \@scores;
}

sub printMotifsFile {
	# i.e. 1 or more motifs (array)
	my ($motifs, $file) = @_;
	open MOTIFFILE, ">$file";
	foreach(@$motifs) {
		my $m = $_;		
		print MOTIFFILE $m->{'text'};
	}
	close MOTIFFILE;
}

sub printMotifFile {
	my ($m,$file) = @_;
	open MOTIFFILE, ">$file";
	print MOTIFFILE $m->{'text'};
	if (0) {
	for (my $i=0;$i<@{$m->{'matrix'}};$i++) {
		for (my $j=0;$j<4;$j++) {
			my $v = $m->{'matrix'}->[$i][$j];
			print MOTIFFILE "\t" if ($j > 0);
			print MOTIFFILE $v;
		}
		print MOTIFFILE "\n";
	}
	}
	close MOTIFFILE;
}

sub logOddsMotifMatrix {
	my ($m) = @_;
	my $minP = 0.001;
	my $mlen = scalar(@{$m->{'matrix'}});
	for (my $i=0;$i<$mlen;$i++) {
		my $sum = 0;
		for (my $j=0;$j<4;$j++) {
			$sum+=$m->{'matrix'}->[$i]->[$j];;
		}
	 	if ($sum>0) {
			my $sum2 = 0;
			for (my $j=0;$j<4;$j++) {
				$m->{'matrix'}->[$i]->[$j] /= $sum;
				$m->{'matrix'}->[$i]->[$j] = $minP if ($m->{'matrix'}->[$i]->[$j] < $minP);
				$sum2+=$m->{'matrix'}->[$i]->[$j];
			}
			for (my $j=0;$j<4;$j++) {
				$m->{'matrix'}->[$i]->[$j] /= $sum2;
				$m->{'matrix'}->[$i]->[$j] = log($m->{'matrix'}->[$i]->[$j]/0.25);
			}
		}
	}
	return $m;
}
sub scoreOligoWithMotif {
	my ($m, $s, $strand) = @_;
	my $score = 0;
	%a = ();
	$a{'A'} = 0;
	$a{'C'} = 1;
	$a{'G'} = 2;
	$a{'T'} = 3;
	$a{'a'} = 0;
	$a{'c'} = 1;
	$a{'g'} = 2;
	$a{'t'} = 3;
	my $mlen = scalar(@{$m->{'matrix'}});

	my $seq = $s;
	if ($strand eq '-' || $strand eq '1') {
		$seq = revopp($s);
	}
	for (my $i=0;$i<length($seq);$i++) {
		if ($i>=$mlen) {
			last;
		}
		my $nt = substr($seq,$i,1);
		if (exists($a{$nt})) {
			$index = $a{$nt};
			my $v = $m->{'matrix'}->[$i][$index];
			#print STDERR "\t$i\t$index\tmlen=$mlen\t$v\n";
			$score += $v;
		} else {
			my $ss = 0;
			for (my $j=0;$j<4;$j++) {
				$ss += $m->{'matrix'}->[$i][$j];
			}
			$ss /= 4.0;
			$score += $ss;
		}
	}
	return $score;
}

sub readTagInfo {
    my ($dir,$pc) = @_;

	my $file = "$dir/" . "/tagInfo.txt";
    my $p=0;
    my $t = 0;
    my $fragEst = "NA";
    my $peakEst = "NA";
	if (!defined($pc)) {
		$pc = 0;
	}
    open IN, $file or die "Could not open $file\n";
    my $count = 0;
    while (<IN>) {
        $count++;
        if ($count == 1) {
            if (/Estimated Fragment Length=(\d+?)[,\t]/) {
                $fragEst = $1;
            }
            if (/Estimated Peak Width=(\d+?)[,\t]/) {
                $peakEst = $1;
            }

        }
        next if ($count < 2);
        chomp;
        s/\r//g;
        my @line= split /\t/;
		if ($count == 2 || $line[0] eq 'genome') {
        	$p = $line[1];
        	$t = $line[2];
		} else {
			if ($line[0] =~ /fragmentLengthEstimate=(.+)$/) {
				$fragEst = $1;
			}
			if ($line[0] =~ /peakSizeEstimate=(.+)$/) {
				$peakEst = $1;
			}
		}
    }
    close IN;
	if ($pc == 1) {	
		$t = $p;
	} elsif ($pc > 1) {
		$file = "$dir/tagCountDistribution.txt";
		if (-e $file) {
			open IN, $file;	
			my $total = 0;
			my $count = 0;
			while (<IN>) {
				$count++;
				chomp;
				s/\r//g;
				my @line = split /\t/;
				next if ($count==1);
				$total += $line[1];
				if ($line[0] eq $pc) {
					last;
				}
			}
			close IN;
			#print STDERR "total=$total\n";
			$t *= $total;
		}
	}

    #print STDERR "$t, $p, $fragEst, $peakEst, $pc\n";
    return ($t, $p, $fragEst, $peakEst);
}

sub getHiCBgRes {
	my ($dir,$currentRes,$cpus) = @_;
	my $tmpFile = rand() . ".tmp";
	`ls -1 "$dir"/HiCbackground_* > "$tmpFile"`;
	open TMPIN, $tmpFile;
	my %res = ();
	while (<TMPIN>) {
		chomp;
		if (/HiCbackground_(\d+)_bp/) {
			$res{$1}=1;
		}
	}
	close TMPIN;
	`rm "$tmpFile"`;
	if (!exists($res{$currentRes})) {
		print STDERR "\t! Background for resolution $currentRes not found for Directory $dir\n";
		print STDERR "\t\tExisting Bg models:\n";
		foreach(keys %res) {
			print STDERR "\t\t\t$_\n";
		}
		print STDERR "\tWill create model for $currentRes in 5 seconds...\n";
		`sleep 5`;
		`analyzeHiC "$dir" -bgonly -res $currentRes -cpu $cpus`;
	}
	return \%res;

}

sub revopp {
	my ($seq) = @_;
	my $rv = reverse($seq);
	$rv =~ s/A/X/g;
	$rv =~ s/T/A/g;
	$rv =~ s/X/T/g;
	$rv =~ s/C/X/g;
	$rv =~ s/G/C/g;
	$rv =~ s/X/G/g;
	return $rv;
}


sub parseCustomGenome {
	my ($genome) = @_;

	my $genomeDir = "";
	my $genomeName = "";
	my $genomeParseDir = "";
	if (-e $genome) {
		$customGenome=1;
		$genomeDir = $genome;
		$genomeParseDir = $genome;
		if (-f $genomeParseDir) {
			if ($genomeParseDir =~ s/\/([^\/]*)$//) {
				$genomeName = $1;
			} else {
				$genomeName = $genome;
				$genomeParseDir = ".";
			}
		} else {
			my $a = $genomeParseDir;
			$a  =~ s/\/+$//;
			$a  =~ s/\/([^\/]*)$//;
			$genomeName = $1;
		}
		#print STDERR "$genomeName\n";
		$genomeParseDir .= "/preparsed/";
		print STDERR "\tUsing Custom Genome\n";
	} else {
		print STDERR "\n!!!!Genome $genome not found in $homeDir/config.txt\n\n";
		print STDERR "\tTo check if is available, run \"perl $homeDir/configureHomer.pl -list\"\n";
		print STDERR "\tIf so, add it by typing \"perl $homeDir/configureHomer.pl -install $genome\"\n";
		print STDERR "\n";
		exit;
	}
	return ($genomeName,$genomeDir,$genomeParseDir);
}

sub parseUCSCStr {
	my ($str) = @_;
	my $chr = "";
	my $start = 0;
	my $end = 1e10;
	$str =~ s/\,//g;
	$str =~ /(.*?)\:(\d+)\-(\d+)/;
	$chr = $1;
	$start = $2;
	$end = $3;
	return ($chr,$start,$end);

}

sub checkMSet {
	my ($mset, $org) = @_;
	if ($mset eq 'auto' && $org ne '' && $org ne 'null') {
		my %table = ();
		open ORGTABLE, "$homeDir/data/knownTFs/organism.table.txt" or die "Where is the file?\n";
		while (<ORGTABLE>) {
			next if (/^#/);
			chomp;
			s/\r//g;
			my @line = split /\t/;
			$table{$line[0]} = $line[1];
		}
		close ORGTABLE;
		if (exists($table{$org})) {
			$mset = $table{$org};
			print STDERR "\tFound mset for \"$org\", will check against $mset motifs\n";
		} else {
			print STDERR "\tCould not find mset for \"$org\", will check against all motifs\n";
			$mset = 'all';
		}
	} elsif ($mset eq 'auto') {
		$mset='all';
	}
	my $msetDir = $homeDir . "/data/knownTFs/$mset/";
	my $allMotifsFile = $msetDir . "all.motifs";
	my $knownMotifsFile = $msetDir . "known.motifs";
	unless (-e $allMotifsFile) {
		print STDERR "\tWarning, couldn't find $allMotifsFile for mset: \"$mset\"\n";
		print STDERR "\t\tUsing all motifs...\n";
		$mset = 'all';
		$msetDir = $homeDir . "/data/knownTFs/$mset/";
		$allMotifsFile = $msetDir . "all.motifs";
		$knownMotifsFile = $msetDir . "known.motifs";
	}
	#print STDERR "$allMotifsFile, $knownMotifsFile\n";
	return ($allMotifsFile, $knownMotifsFile);

}

sub readAnnTable {
	my ($file) = @_;

	my %table = ();
	my @header = ();

	open ANNTABLE, "$file" or die "!!! Error: Could not open file: $file\n";
	my $zz=0;
	while (<ANNTABLE>) {
		$zz++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($zz == 1) {
			@header = @line;
			next;
		}
		my $id= $line[0];
		my %d = ();
		for (my $i=0;$i<@line;$i++) {
			last if ($i >= @header);
			my $n = $header[$i];
			$d{$n} = $line[$i];
		}
		$table{$id} = \%d;
	}
	close ANNTABLE;
	return \%table;
}


1;

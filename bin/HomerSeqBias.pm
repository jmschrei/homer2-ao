#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";


package HomerSeqBias;

use POSIX;




sub generateSequences {
	my ($seqs,$markov,$numRandSeq,$maxCPUs) = @_;

	print STDERR "\tBuilding sequence model...\n";
	my $model = buildPosDepMarkovModel($seqs,$markov);
	
	print STDERR "\tGenerating sequences...\n";

	my $seqPerSet = ceil($numRandSeq/$maxCPUs);

	my $tmpFile = rand() . ".tmp";

	my @files = ();
	for (my $k=0;$k<$maxCPUs;$k++) {

		my $file = $k . "_" . $tmpFile;
		push(@files, $file);

		my $pid = fork();
		if ($pid == 0) {
			srand();
			open OUT, ">$file";
			for (my $i=0;$i<$seqPerSet;$i++) {
				my $seq = generateRandSequence($model);
				print OUT "$seq\n";
			}
			close OUT;
			exit;
		}
	}
	my $id = 0;
	while ($id >= 0) {
		$id = wait();
	}

	my %newSeqs = ();
	my $x = 0;
	foreach(@files) {
		my $file = $_;
		open IN, $_;
		while (<IN>) {
			chomp;
			my $id = "rand_" . ++$x;
			last if ($x >= $numRandSeq);
			$newSeqs{$id} = $_;
		}
		close IN;
		`rm "$file"`;
	}

	print STDERR "\tTotal Sequences created: $x\n";
	
	return \%newSeqs;
}



sub generateSequencesInParallel {
	my ($model,$posIndepFlag,$size,$numRandSeq,$maxCPUs) = @_;
	
	print STDERR "\tGenerating sequences using $maxCPUs CPUs...\n";

	my $seqPerSet = ceil($numRandSeq/$maxCPUs);

	my $tmpFile = rand() . ".tmp";

	my @files = ();
	for (my $k=0;$k<$maxCPUs;$k++) {
		my $file = $k . "_" . $tmpFile;
		push(@files, $file);

		my $pid = fork();
		if ($pid == 0) {
			srand();
			open OUT, ">$file";
			for (my $i=0;$i<$seqPerSet;$i++) {
				my $seq = ();
				if ($posIndepFlag) {
					$seq = HomerSeqBias::generatePosIndepRandSequence($model,$size);
				} else {
					$seq = HomerSeqBias::generateRandSequence($model);
				}
				print OUT "$seq\n";
			}
			close OUT;
			exit;
		}
	}
	my $id = 0;
	while ($id >= 0) {
		$id = wait();
	}

	my %newSeqs = ();
	my $x = 0;
	foreach(@files) {
		my $file = $_;
		open IN, $_;
		while (<IN>) {
			chomp;
			my $id = "rand_" . ++$x;
			last if ($x >= $numRandSeq);
			$newSeqs{$id} = $_;
		}
		close IN;
		`rm "$file"`;
	}
	print STDERR "\tTotal Sequences created: $x\n";
	
	return \%newSeqs;
}

sub readFASTA {
	my ($file) = @_;
	my %seqs = ();
	my $curName = '';
	my $curSeq = '';
	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		next if ($_ eq '');
		if (s/^>//) {
			if ($curName ne '' && $curSeq ne '') {
				$seqs{$curName} = $curSeq;
			}
			$curName = $_;
			$curSeq = '';
		} else {
			$curSeq .= $_;
		}
	}
	close IN;
	if ($curName ne '' && $curSeq ne '') {
		$seqs{$curName} = $curSeq;
	}
	return \%seqs;
}

# Takes a hashtable of sequences and the markov chain order (integer)
# and creates a model (custom data structure made up of hashtables)
# that can be used to generate random sequences with similar properties
# Sequences should all be the same length
# Arguments:
#   1. $seqs - reference to hashtable where keys are sequence names and values are sequences
#   2. $markov - order (i.e. 1 for dinucleotides, 2 for trinucleotides, etc.)
sub buildPosDepMarkovModel {
	my ($seqs,$markov) = @_;

	my @seqs = values %$seqs;

	my $seqLen = length($seqs[0]);
	my $len = $markov+1;

	my %model = ();
	$model{'len'} =$len;
	$model{'seqLen'} = $seqLen;
	$model{'type'} = 'posdep';
	print STDERR "\tsequence length=$seqLen\n";
	print STDERR "\tk-mer length for model=$len\n";
	for (my $i=0;$i<$seqLen;$i++) {
		my %a = ();
		my %b = ();
		$model{$i} = \%a;
		$model{$i}->{'total'} = 0;

		$model{$i}->{'reverse'} = \%b;
		$model{$i}->{'reverse'}->{'total'} = 0;
	}

	my $maxP = 0;
	for (my $i=0;$i<@seqs;$i++) {
		my $seq = $seqs[$i];
		if (length($seq) != $seqLen) {
			my $LL = length($seq);
			print STDERR "!!!Error: Sequences must be the same length for a position dependent model ($seqLen vs. $LL)\n";
			exit;
		}
		my $s0 = substr($seq,0,$len);
		for (my $j=1;$j<$seqLen-$len+2;$j++) {
			my $s1 = substr($seq,$j,$len);
			if (!defined($s1)) {
				print STDERR "$j $len $seq\n";
			}
			my $p = $j-1;
			$maxP = $p if ($p > $maxP);

			if (!exists($model{$p}->{$s0})) {
				my %a = ();
				$model{$p}->{$s0} = \%a;
				$model{$p}->{$s0}->{'total'} = 0;
			}
			$model{$p}->{$s0}->{$s1}++;
			$model{$p}->{$s0}->{'total'}++;
			$model{$p}->{'total'}++;

			my $pp = $p+1;
			if (!exists($model{$pp}->{'reverse'}->{$s1})) {
				my %a = ();
				$model{$pp}->{'reverse'}->{$s1} = \%a;
				$model{$pp}->{'reverse'}->{$s1}->{'total'} = 0;
			}
			$model{$pp}->{'reverse'}->{$s1}->{$s0}++;
			$model{$pp}->{'reverse'}->{$s1}->{'total'}++;
			$model{$pp}->{'reverse'}->{'total'}++;

			$s0 = $s1;
		}
	}
	print STDERR "\tmaxP=$maxP\n";
	$model{'maxP'} = $maxP;

	#now normalize model	
	for (my $i=0;$i<=$maxP;$i++) {
		my $p = $i;
		my @randArray = ();
		my $ptotal = $model{$p}->{'total'};
		my $sum =0;
		foreach (keys %{$model{$p}}) {
			my $key = $_;
			next if ($key eq 'total' || $key eq 'rarray' || $key eq 'reverse' || $key eq 'maxP');
			my $rate = $model{$p}->{$_}->{'total'} / $ptotal;
			$sum+=$rate;
			push(@randArray,{s=>$sum,k=>$_});
		}
		$model{$p}->{'rarray'} = \@randArray;

		foreach (keys %{$model{$p}}) {
			my $key = $_;
			my @randArray = ();
			next if ($key eq 'total' || $key eq 'rarray' || $key eq 'reverse' || $key eq 'maxP');
			my $ktotal = $model{$p}->{$key}->{'total'};
			my $sum =0;
			foreach (keys %{$model{$p}->{$key}}) {
				next if ($_ eq 'total' || $_ eq 'rarray' || $_ eq 'reverse' || $_ eq 'maxP');
				my $rate = $model{$p}->{$key}->{$_} / $ktotal;
				$sum+=$rate;
				push(@randArray,{s=>$sum,k=>$_});
			}
			$model{$p}->{$key}->{'rarray'} = \@randArray;
		}

		foreach (keys %{$model{$p}->{'reverse'}}) {
			my $key = $_;
			next if ($key eq 'total' || $key eq 'rarray' || $key eq 'reverse' || $key eq 'maxP');
			my @randArray = ();
			my $ktotal = $model{$p}->{'reverse'}->{$key}->{'total'};
			my $sum =0;
			foreach (keys %{$model{$p}->{'reverse'}->{$key}}) {
				next if ($_ eq 'total' || $_ eq 'rarray' || $_ eq 'reverse' || $_ eq 'maxP');
				my $rate = $model{$p}->{'reverse'}->{$key}->{$_} / $ktotal;
				$sum+=$rate;
				push(@randArray,{s=>$sum,k=>$_});
			}
			$model{$p}->{'reverse'}->{$key}->{'rarray'} = \@randArray;
		}
	}
	return \%model;
}

# Takes a hashtable of sequences and the markov chain order (integer)
# and creates a model (custom data structure made up of hashtables)
# that can be used to generate random sequences with similar properties
# Sequences do NOT need to be the same length for the position independent model
# Arguments:
#   1. $seqs - reference to hashtable where keys are sequence names and values are sequences
#   2. $markov - order (i.e. 1 for dinucleotides, 2 for trinucleotides, etc.)
sub buildPosIndepMarkovModel {
	my ($seqs,$markov) = @_;

	my @seqs = values %$seqs;

	my $len = $markov+1;

	my %model = ();
	$model{'len'} =$len;
	$model{'total'} = 0;
	$model{'type'} = 'posindep';
	my %zz = ();
	$model{'reverse'} = \%zz;
	$model{'reverse'}->{'total'} = 0;

	for (my $i=0;$i<@seqs;$i++) {
		my $seq = $seqs[$i];
		my $seqLen = length($seq);
		my $s0 = substr($seq,0,$len);
		for (my $j=1;$j<$seqLen-$len+1;$j++) {
			my $s1 = substr($seq,$j,$len);

			if (!exists($model{$s0})) {
				my %a = ();
				$model{$s0} = \%a;
				$model{$s0}->{'total'} = 0;
			}
			$model{$s0}->{$s1}++;
			$model{$s0}->{'total'}++;
			$model{'total'}++;

			if (!exists($model{'reverse'}->{$s1})) {
				my %a = ();
				$model{'reverse'}->{$s1} = \%a;
				$model{'reverse'}->{$s1}->{'total'} = 0;
			}
			$model{'reverse'}->{$s1}->{$s0}++;
			$model{'reverse'}->{$s1}->{'total'}++;
			$model{'reverse'}->{'total'}++;
			$s0 = $s1;
		}
	}

	#now normalize model	
	my @randArray = ();
	my $ptotal = $model{'total'};
	my $sum =0;
	foreach (keys %model) {
		my $key = $_;
		next if ($key eq 'total' || $key eq 'rarray' || $key eq 'reverse' || $key eq 'maxP' || $key eq 'len' || $key eq 'type');
		#print STDERR "$key $model{$key}->{'total'} / $ptotal\n";
		my $rate = $model{$key}->{'total'} / $ptotal;
		$sum+=$rate;
		push(@randArray,{s=>$sum,k=>$key});
	}
	$model{'rarray'} = \@randArray;

	foreach (keys %model) {
		my $key = $_;
		my @randArray = ();
		next if ($key eq 'total' || $key eq 'rarray' || $key eq 'reverse' || $key eq 'maxP' || $key eq 'len' || $key eq 'type');
		my $ktotal = $model{$key}->{'total'};
		my $sum =0;
		foreach (keys %{$model{$key}}) {
			next if ($_ eq 'total' || $_ eq 'rarray' || $_ eq 'reverse' || $_ eq 'maxP' || $_ eq 'len' || $_ eq 'type');
			my $rate = $model{$key}->{$_} / $ktotal;
			$sum+=$rate;
			push(@randArray,{s=>$sum,k=>$_});
		}
		$model{$key}->{'rarray'} = \@randArray;
	}

	foreach (keys %{$model{'reverse'}}) {
		my $key = $_;
		next if ($key eq 'total' || $key eq 'rarray' || $key eq 'reverse' || $key eq 'maxP' || $key eq 'len' || $key eq 'type');
		my @randArray = ();
		my $ktotal = $model{'reverse'}->{$key}->{'total'};
		my $sum =0;
		foreach (keys %{$model{'reverse'}->{$key}}) {
			next if ($_ eq 'total' || $_ eq 'rarray' || $_ eq 'reverse' || $_ eq 'maxP' || $_ eq 'len' || $_ eq 'type');
			my $rate = $model{'reverse'}->{$key}->{$_} / $ktotal;
			$sum+=$rate;
			push(@randArray,{s=>$sum,k=>$_});
		}
		$model{'reverse'}->{$key}->{'rarray'} = \@randArray;
	}
	return \%model;
}



# This function will generate a random sequence and return it
# argumeents: 
#  1. $model  - position depdennt sequence model.
sub generatePosIndepRandSequence {
	my ($model,$length) = @_;

	my $len = $model->{'len'};
	my $seqLen = $length;

	my $r = rand();
	my $s0 = '';
	my $seq = '';
	foreach(@{$model->{'rarray'}}) {
		if ($r < $_->{'s'}) {
			$s0 = $_->{'k'};
			last;
		}
	}
	$seq .= $s0;

	for (my $i=0;$i<$seqLen-$len+1;$i++) {
		$r = rand();
		$s1 = '';
		foreach(@{$model->{$s0}->{'rarray'}}) {
			if ($r < $_->{'s'}) {
				$s1 = $_->{'k'};
				last;
			}
		}	
		$seq .= substr($s1,$len-1,1);
		$s0 = $s1;
	}
	return $seq;
}




# This function will generate a random sequence and return it
# argumeents: 
#  1. $model  - position depdennt sequence model.
sub generateRandSequence {
	my ($model) = @_;

	my $len = $model->{'len'};
	my $seqLen = $model->{'seqLen'};
	my $maxP = $model->{'maxP'};
	my $r = rand();
	my $s0 = '';
	my $seq = '';
	foreach(@{$model->{0}->{'rarray'}}) {
		if ($r < $_->{'s'}) {
			$s0 = $_->{'k'};
			last;
		}
	}
	$seq .= $s0;

	for (my $i=0;$i<$seqLen-$len+1;$i++) {
		my $p = $i;
		$r = rand();
		$s1 = '';
		foreach(@{$model->{$p}->{$s0}->{'rarray'}}) {
			if ($r < $_->{'s'}) {
				$s1 = $_->{'k'};
				last;
			}
		}	
		$seq .= substr($s1,$len-1,1);
		$s0 = $s1;
	}
	return $seq;
}


# This function will generate a random sequence and return it
# argumeents: 
#  1. $model  - position depdennt sequence model.
sub generateRandSequenceReverse {
	my ($model) = @_;

	my $len = $model->{'len'};
	my $seqLen = $model->{'seqLen'};
	my $maxP = $model->{'maxP'};
	my $r = rand();
	my $s1 = '';
	my $seq = '';
	foreach(@{$model->{$maxP}->{'rarray'}}) {
		if ($r < $_->{'s'}) {
			$s1 = $_->{'k'};
			last;
		}
	}
	$seq = $s1 . $seq;

	for (my $i=$maxP;$i>=0;$i--) {
		my $p = $i;
		$r = rand();
		$s0 = '';
		foreach(@{$model->{$p}->{'reverse'}->{$s1}->{'rarray'}}) {
			if ($r < $_->{'s'}) {
				$s0 = $_->{'k'};
				last;
			}
		}	
		my $toadd = substr($s0,0,1);
		$seq = $toadd . $seq;
		$s1 = $s0;
	}
	return $seq;
}

# This function will generate a random sequence and return it
# argumeents: 
#  1. $model  - position depdennt sequence model.
sub generateRandSequencefromSeedSeq {
	my ($model, $seedSeq, $seedPos,$seqLen) = @_;


	my $len = $model->{'len'};
	if (!defined($seqLen)) {
		$seqLen = $model->{'seqLen'};
	}
	my $maxP = -1;
	if ($model->{'type'} eq 'posdep') {
		$maxP = $model->{'maxP'};
	} else {
		$maxP = $seqLen-$len;
	}

	my $seedLen = length($seedSeq);

	if ($seedPos + $seedLen-1 < 0 ) {
		return '';
	}
	if ($seedPos < 0) {
		my $s = substr($seedSeq,-1*$seedPos);
		$seedSeq = $s;
		$seedPos = 0;
	}

	if ($seedPos >= $seqLen) {
		return '';
	}
	if ($seedPos+$seedLen > $seqLen) {
		my $s = substr($seedSeq,0,$seqLen-$seedPos);
		$seedSeq = $s;
	}
	$seedLen = length($seedSeq);
	#prepart.
	# First determine the starting kmer that matches the seed seqeunce
	my $preseq = '';
	my $s0='';
	my $s1='';
	if ($seedLen < $len) {
		my $startPos = $seedPos;
		my $offset = 0;
		if ($seedPos >= $maxP) {
			$startPos = $maxP;
			$offset = $seedPos-$maxP;
			$seedPos = $maxP;
		}

		my $good = 0;
		while ($good == 0) {

			my $r = rand();
			my $rarray = ();
			if ($model->{'type'} eq 'posdep') {
				$rarray = $model->{$startPos}->{'rarray'};
			} else {
				$rarray = $model->{'rarray'};
			}
			foreach(@$rarray) {
				if ($r < $_->{'s'}) {
					$s1 = $_->{'k'};
					last;
				}
			}
			if (substr($s1,$offset,$seedLen) eq $seedSeq) {
				$good=1;
				last;
			}
		}
	} elsif ($seedLen == $len) {
		$s1 = $seedSeq;
	} else {
		$s1 = substr($seedSeq,0,$len);
	}
	my $seedSeqToAdd = $seedSeq;
	my $initialSeed = $s1;

	for (my $i=$seedPos;$i>0;$i--) {
		my $p = $i;
		$r = rand();
		$s0 = '';
		my $rarray = '';
		if ($model->{'type'} eq 'posdep') {
			$rarray = $model->{$p}->{'reverse'}->{$s1}->{'rarray'};
		} else {
			$rarray = $model->{'reverse'}->{$s1}->{'rarray'};
		}
		foreach(@$rarray) {
			if ($r < $_->{'s'}) {
				$s0 = $_->{'k'};
				last;
			}
		}	
		my $toadd = substr($s0,0,1);
		$preseq = $toadd . $preseq;
		$s1 = $s0;
	}


	#prepart.
	$postseq = '';
	$s0='';
	$s1='';
	my $postStart = $seedPos;
	if ($seedLen < $len) {
		$s0 = $initialSeed;
		$seedSeqToAdd = $initialSeed;
	} elsif ($seedLen == $len) {
		$s0 = $initialSeed; # i.e. seedSeq
		$seedSeqToAdd = $initialSeed;
	} else {
		$s0 = substr($seedSeq,$seedLen-$len,$len);
		$postStart = $seedPos + ($seedLen-$len);
	}

	for (my $i=$postStart;$i<$seqLen-$len;$i++) {
		my $p = $i;
		$r = rand();
		$s1 = '';
		my $rarray = '';
		if ($model->{'type'} eq 'posdep') {
			$rarray = $model->{$p}->{$s0}->{'rarray'};
		} else {
			$rarray = $model->{$s0}->{'rarray'};
		}
		foreach(@$rarray) {
			if ($r < $_->{'s'}) {
				$s1 = $_->{'k'};
				last;
			}
		}	
		if ($s1 eq '') {
			last;
		}
		$postseq .= substr($s1,$len-1,1);
		$s0 = $s1;
	}
	#my $seq = $preseq . "-" . $seedSeqToAdd . "-" . $postseq;
	my $seq = $preseq . $seedSeqToAdd . $postseq;
	#print STDERR "$seq\n";
	return $seq;
}



# This function will generate a random sequence and return it
# argumeents: 
#  1. $model  - position depdennt sequence model.
#  2. ref - sequence of the reference allele
#  3. alt - sequence of the alt allele
#  4. position relative to the beginning (0) where the reference starts.
#  This will return an array with 2 sequences (one for the ref, one for alt)
#  It is also assumed that these variants will overlap
sub generateRandSequenceWithVariant {
	my ($model, $ref, $refPos, $alt, $altPos,$seqLen) = @_;

	my $len = $model->{'len'};

	if (!defined($seqLen)) {
		if (!defined($model->{'seqLen'})) {
			print STDERR "!!! Error - need to specify $seqLen in generateRandSequenceWithVariant when using a position independent model\n";
			return ('','');
		}
		$seqLen = $model->{'seqLen'};
	}


	my $seedFlag = 0; #ref
	if (length($ref) < length($alt)) {
		$seedFlag = 0;
	} elsif (length($ref) > length($alt)) {
		$seedFlag = 1;
	} else {
		my $r = rand();
		if ($r > 0.5) {
			$seedFlag = 1;
		} else {
			$seedFlag = 0;
		}
	}
	my $rseq1 = '';
	my $rseq2 = '';

	#print STDERR "seedFlag=$seedFlag\n";	
	if ($seedFlag==0) {
		#use ref
		#print STDERR "ref=$ref, $refPos\n";
		$rseq1 = generateRandSequencefromSeedSeq($model, $ref, $refPos,$seqLen);
		my $preseq = substr($rseq1,$refPos-$altPos,$altPos);
		#print STDERR "rseq1=$rseq1\nrefPos=$refPos\nref=$ref\nseqLen=$seqLen\naltPos=$altPos\nalt=$alt\n";
		my $postseq = substr($rseq1,$refPos+length($ref),$seqLen-($altPos+length($alt)));
		#$rseq2 = $preseq . "|" .  $alt . "|"  . $postseq;
		$rseq2 = $preseq .  $alt . $postseq;
		
	} else {
		#use alt
		#print STDERR "alt=$alt, $altPos\n";
		$rseq2 = generateRandSequencefromSeedSeq($model, $alt, $altPos,$seqLen);
		my $preseq = substr($rseq2,$altPos-$refPos,$refPos);
		my $postseq = substr($rseq2,$altPos+length($alt),$seqLen-($refPos+length($ref)));
		#$rseq1 = $preseq . "|" . $ref . "|" . $postseq;
		$rseq1 = $preseq . $ref . $postseq;
	}

	return ($rseq1,$rseq2);
}

1;

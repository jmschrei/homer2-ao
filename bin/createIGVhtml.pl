#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/bennerlab/home/cbenner/software/homer/bin";
my $homeDir = "/gpfs/data01/bennerlab/home/cbenner/software/homer/";

use HomerConfig;
use HomerSVGLogo;

use POSIX;



sub printCMD {
	print STDERR "\n\tcreateIGVpage.pl [options] -i <keyFile> -o <outputFileDirectoryName>\n";
	print STDERR "\n\tGenerally speaking, it's probably best to run this with just a single parameter - the -i key file\n";
	print STDERR "\n\tRequired Options:\n";
	print STDERR "\t\t-i <keyFile> (format below)\n";
	print STDERR "\t\tAlthough required, options below can be encoded in the key file too (see below)\n";
	print STDERR "\t\t-o <Output Directory> (name of directory to place web files, i.e. /home/cbenner/www/igv)\n";
	print STDERR "\t\t-url <Output URL> (base URL to access web files, i.e. http://homer.ucsd.edu/cbenner/igv/)\n";
	print STDERR "\t\tShared resource files (for making multiple sets that can use the same genome/annotation files:\n";
	print STDERR "\t\t(by default these will simply be included in the output file and shared resources will not be used)\n";
	print STDERR "\t\t\t-resource <directory> (directory to put shared genome resource files i.e. genome fasta, gtf)\n";
	print STDERR "\t\t\t-resourceURL <URL> (base URL for shared genome files/annotation)\n";
	print STDERR "\t\t\t-org <organism name> (i.e. 'homo_sapiens' will use default genome/gtf from species data)\n";
	print STDERR "\t\t\tNote - organism name should have an underscore instead of a space\n";
	print STDERR "\t\tFor custom genomes:\n";
	print STDERR "\t\t\t-genome <genome FASTA file>\n";
	print STDERR "\t\t\t-gtf <gtf file>\n";
	print STDERR "\t\t-byChr (visualize by chromosome)\n";
	print STDERR "\t\tNormalization options:\n";
	print STDERR "\t\t\t-raw (raw read counts, default: normalized to 10^7)\n";

	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-force (force all files to be generate, otherwise checks if they exist)\n";
	print STDERR "\t\t-fragLength <#> (Force use this fragment length for chip/gro visualization)\n";
	print STDERR "\t\t-color <random|gradient|255,255,255> (if color is not specified use this,def:gradient)\n";

	print STDERR "\n\tKey File format columns:\n";
	print STDERR "\t\tColumn 1: Name\n";
	print STDERR "\t\tColumn 2: TagDirectory or File\n";
	print STDERR "\t\tColumn 3: Type (rna,tss,chip,gro,gtf,bed,bw,gff3)\n";
	print STDERR "\t\tColumn 4: Color (255,255,255)\n";
	print STDERR "\t\tColumn 5: other bigwig/makeUCSCfile parameters\n";
	print STDERR "\n\tExample key file:\n";
	print STDERR "\t\t#trackName<tab>TagDir<tab>Type<tab>Color<tab>other-parameters\n";
	print STDERR "\t\tEmbryo-csRNA	Embryo-Startseq/	tss	255,0,0\n";
	print STDERR "\t\tEmbryo-MNase	Embryo-MNase-seq/	chip	0,75,0	-fragLength 147\n";
	print STDERR "\n\t\t-- file can also contain key/value pairs for required parameters:\n";
	print STDERR "\t\to=/home/cbenner/www/igv/<outputDirName>\n";
	print STDERR "\t\turl=http://homer.ucsd.edu/cbenner/igv/<outputDirName>\n";
	print STDERR "\t\tresource=/home/cbenner/www/igv/genomes/\n";
	print STDERR "\t\tresourceURL=http://homer.ucsd.edu/cbenner/igv/genomes/\n";
	print STDERR "\t\torg=coccidioides_immitis\n";
	print STDERR "\t\tgenome=cocci.fasta\n";
	print STDERR "\t\tgtf=cocci.genes.gtf\n";
	print STDERR "\t\tbyChr=0\n";
	print STDERR "\n";
	exit;
}

my %igvcolors = ();
$igvcolors{"csRNA.1nt"} = "255,0,0";
$igvcolors{"csRNAinput.1nt"} = "55,5,5";
$igvcolors{"mSTART.1nt"} = "255,0,0";
$igvcolors{"mSTARTinput.1nt"} = "55,5,5";
$igvcolors{"5GRO.1nt"} = "100,100,0";
$igvcolors{"GRO.1nt"} = "30,30,5";
$igvcolors{"CAGE.1nt"} = "100,0,100";
$igvcolors{"CAGEinput.1nt"} = "30,5,30";
$igvcolors{"ATAC.cov"} = "0,255,0";
$igvcolors{"GRO"} = "150,150,0";
$igvcolors{"RNA"} = "0,0,255";
$igvcolors{"gtf"} = "0,200,200";
$igvcolors{"bed"} = "200,200,0";

my @colorStart = (150,0,0);
my @colorEnd = (0,0,150);

my $sharedResourceFlag = 0;
my $visualizeByChrFlag = 0;

my $keyFile = "";
my $outputDirectory = "";
my $org = '';
my $force = 0;
my $includeRNAbydefault = 0;
my $defaultFragLength = 75;
my $chipFragLength = -1e9;
my $colorScheme = "gradient";
my $customGenome = '';
my $customGTF = '';
my $normOpts = '';

my $baseDir = "/home/cbenner/www/igv/";
my $baseURL = "http://homer.ucsd.edu/cbenner/igv/";
my $resourceDir = "/home/cbenner/www/test/species-ref/";
my $resourceURL = "http://homer.ucsd.edu/cbenner/test/species-ref/";
my $orgSource = "/gpfs/data01/bennerlab/species/genomes/";
my $speciesDir = "/gpfs/data01/bennerlab/species/";


printCMD () if (@ARGV < 2);
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-i' || $ARGV[$i] eq '-k') {
		$keyFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-o' || $ARGV[$i] eq '-output') {
		$outputDirectory = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-url') {
		$baseURL = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-resource') {
		$resourceDir = $ARGV[++$i];
		$sharedResourceFlag = 1;
	} elsif ($ARGV[$i] eq '-resourceURL') {
		$resourceURL = $ARGV[++$i];
		$sharedResourceFlag = 1;
	} elsif ($ARGV[$i] eq '-org') {
		$org = $ARGV[++$i];
		#$sharedResourceFlag = 1;
	} elsif ($ARGV[$i] eq '-force') {
		$force = 1;
	} elsif ($ARGV[$i] eq '-genome') {
		$customGenome = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-gtf') {
		$customGTF = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-color') {
		$colorScheme = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-raw') {
		$normOpts .= " -raw ";
	} elsif ($ARGV[$i] eq '-byChr') {
		$visualizeByChrFlag = 1;
	} elsif ($ARGV[$i] eq '-fragLength') {
		$defaultFragLength = $ARGV[++$i];
		$chipFragLength = $defaultFragLength;
	} else {
		printCMD();
	}
}


if ($keyFile eq '') {
	print STDERR "!! Key file (-i) required!\n";
	exit;
}

my $rand = rand();
my $tmpFile = $rand . ".tmp";

my $numNeedColor = 0;
	
print STDERR "\n\tExperiments to visualize:\n";
open IN, $keyFile;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^\#/);
	my $og = $_;
	if (/org\s*=\s*(.*?)\s*$/) {
		$org = $1;
		#$sharedResourceFlag = 1;
		next;
	}
	if (/genome\s*=\s*(.*?)\s*$/) {
		$customGenome = $1;
		next;
	}
	if (/gtf\s*=\s*(.*?)\s*$/) {
		$customGTF = $1;
		next;
	}
	if (/output\s*=\s*(.*?)\s*$/) {
		$outputDirectory = $1;
		next;
	}
	if (/o\s*=\s*(.*?)\s*$/) {
		$outputDirectory = $1;
		next;
	}
	if (/url\s*=\s*(.*?)\s*$/) {
		$baseURL = $1;
		next;
	}
	if (/resource\s*=\s*(.*?)\s*$/) {
		$resourceDir = $1;
		$sharedResourceFlag = 1;
		next;
	}
	if (/resourceURL\s*=\s*(.*?)\s*$/) {
		$resourceURL = $1;
		$sharedResourceFlag = 1;
		next;
	}
	if (/byChr\s*=\s*(.*?)\s*$/) {
		$visualizeByChrFlag = 1;
		next;
	}

	if (/fragLength\s*=\s*(.*?)\s*$/) {
		$defaultFragLength = $1;
		$chipFragLength = $defaultFragLength;
		next;
	}
	my @line = split /\t/;
	if (@line < 2) {
		print STDERR "\t!!!Warning too few columns for $_\n";
		next;
	}
	my $name = $line[0];
	my $file = $line[1];
	my $type = "chip";
	$type = $line[2] if (@line > 2);
	my $color = "";
	$color = $line[3] if (@line > 3);
	my $params = '';
	$params = $line[4] if (@line > 4);
	my $sampleForce = 0;
	if ($params =~ s/-force//) {
		#print STDERR "Found force flag!\n";
		$sampleForce=1;
	}


	print STDERR "\t\t$name\t$file\t$type\t$color\t$params\n"; 
	my $t = {n=>$name,f=>$file,t=>$type,c=>$color,p=>$params,force=>$sampleForce};
	if ($type eq 'bed' || $type eq 'gtf' || $type eq 'gff3') {
		unless (-f $file) {
			print STDERR "!!Error - could find file for $og\n";
			next;
		}
		push(@anns, $t);
	} elsif ($type eq 'bw') {
		if ($file =~ /\,/) {
			print STDERR "\n";
		} else {
			unless (-f $file) {
				print STDERR "!!Error - could find file for $og\n";
				next;
			}
		}
		push(@data, $t);
	} else {
		unless (-d $file) {
			print STDERR "!!Error - could find directory for $og\n";
			next;
		}
		if ($color eq '') {
			$numNeedColor++;
		}
		push(@data, $t);
	}
}
close IN;


if ($outputDirectory eq '') {
	print STDERR "!! Output directory (-o) required!\n";
	exit;
}

if ($org eq '') {
	print STDERR "!! Organism required (-org)!\n";
	exit;
}

if ($customGenome eq '') {
	$sharedResourceFlag = 1;
	unless (-e $orgSource . $org) {
		print STDERR "!!! -org $org not recognized as an organism\n";
		$x = `ls -d $orgSource/*`;
		my @orgs = split /\s/,$x;
		foreach(@orgs) {
			$_ =~ s/\/+$//;
			$_ =~ s/^.*\///;
		}
		print STDERR "\tPossible -org: @orgs\n";
		exit;
	}
}



$baseDir = $outputDirectory;
#$baseDir .= $outputDirectory . "/";
#$baseURL .= $outputDirectory . "/";
$baseURL .= "/" unless ($baseURL =~ /\/$/);
$baseDir .= "/" unless ($baseDir =~ /\/$/);

if ($sharedResourceFlag==0) {
	$resourceDir = $baseDir;
	$resourceURL = $baseURL;
}

$resourceDir .= "/" unless ($resourceDir =~ /\/$/);
$resourceURL .= "/" unless ($resourceURL =~ /\/$/);
print STDERR "\n\tOutputDirectory: $baseDir\n";
print STDERR "\tURL: $baseURL\n";
print STDERR "\n\tResourceDirectory: $resourceDir\n";
print STDERR "\tResourceURL: $resourceURL\n";

my $orgDir = $resourceDir . $org . "/";
my $orgURL = $resourceURL . $org . "/";
print STDERR "\n\tSpecies Directory: $orgDir\n";
print STDERR "\tSpecies URL: $orgURL\n";

if ($customGenome eq '') {
	unless (-e $orgSource . "$org/genome.fa") {
		print STDERR "!!! Error - could not find genome file for $org\n";
		exit;
	}
}

if (-e $orgDir . "genome.fa") {
	print STDERR "\t\tFound genome already in webresource folder - skipping\n";
} else {

}
	


my $j = 0;
$numNeedColor--;
$numNeedColor=1 if ($numNeedColor < 1);
for (my $i=0;$i<@data;$i++) {
	if ($data[$i]->{'c'} eq '') {
		if ($colorScheme eq 'rand' || $colorScheme eq 'random') {
			$data[$i]->{'c'} = floor(255*rand()) . "," . floor(255*rand()) . "," . floor(255*rand());
		} elsif ($colorScheme eq 'gradient') {
			my $red = floor($colorStart[0] + ($colorEnd[0]-$colorStart[0])/$numNeedColor*$j);
			my $green = floor($colorStart[1] + ($colorEnd[1]-$colorStart[1])/$numNeedColor*$j);
			my $blue = floor($colorStart[2] + ($colorEnd[2]-$colorStart[2])/$numNeedColor*$j);
			$j++;
			$data[$i]->{'c'} = "$red,$green,$blue";
		}
	}
}
for (my $i=0;$i<@anns;$i++) {
	if ($anns[$i]->{'c'} eq '') {
		if ($colorScheme eq 'rand' || $colorScheme eq 'random') {
			$anns[$i]->{'c'} = floor(255*rand()) . "," . floor(255*rand()) . "," . floor(255*rand());
		} elsif ($colorScheme eq 'gradient') {
			my $red = floor($colorStart[0] + ($colorEnd[0]-$colorStart[0])/$numNeedColor*$j);
			my $green = floor($colorStart[1] + ($colorEnd[1]-$colorStart[1])/$numNeedColor*$j);
			my $blue = floor($colorStart[2] + ($colorEnd[2]-$colorStart[2])/$numNeedColor*$j);
			$j++;
			$anns[$i]->{'c'} = "$red,$green,$blue";
		}
	}
}






print STDERR "\n\tPausing for 3 seconds to make sure this looks ok...\n";
`sleep 3`;

`mkdir -p $baseDir`;
if (-e $orgDir) {
} else {
	print STDERR "\tCreating organism igv resource directory ($orgDir)\n";
	`mkdir -p $orgDir`;
}

if (-e $orgDir . "genome.fa" && $force == 0) {
} else {
	print STDERR "Can't find $orgDir" . "genome.fa\n\n";
	print STDERR "\n\tCopying organism files to shared folder (genome.fa, genes.gtf)\n";
	if ($customGenome eq '') {
		`cp "$orgSource/$org/genome.fa" "$orgDir/genome.fa"`;
		`cp "$orgSource/$org/genome.fa.fai" "$orgDir/genome.fa.fai"`;
		`cp "$orgSource/$org/genes.gtf" "$orgDir/genes.gtf"`;
		`cp "$orgSource/$org/genes.gtf.idx" "$orgDir/genes.gtf.idx"`;
	} else {
		`samtools faidx "$customGenome"`;
		`cp "$customGenome" "$orgDir/genome.fa"`;
		`cp "$customGenome.fai" "$orgDir/genome.fa.fai"`;


		if ($customGTF ne '') {
			`sortBed -i "$customGTF" > "$tmpFile.gtf"`;
			`igvtools index "$tmpFile.gtf"`;
			`mv "$tmpFile.gtf" "$orgDir/genes.gtf"`;
			`mv "$tmpFile.gtf.idx" "$orgDir/genes.gtf.idx"`;
		}
	}
}
$genomeURL = $resourceURL . "/$org/genome.fa";
$genomeIndexURL = $resourceURL . "/$org/genome.fa.fai";
$genesURL = $resourceURL . "/$org/genes.gtf";
$genesIndexURL = $resourceURL . "/$org/genes.gtf.idx";


my $chromSizes = "$orgSource/$org/chrom.sizes";
if ($customGenome ne '') {
	$chromSizes = "$orgDir/chrom.sizes";
	if (-e $chromSizes && $force == 0) {
	} else {
		`homerTools extract stats "$orgDir/genome.fa" > "$chromSizes.header"`;
		`sed '1d' "$chromSizes.header" > "$chromSizes"`;
		`rm "$chromSizes.header"`;
	}
}



open IGV, ">$baseDir/index.html";

print IGV "<!DOCTYPE html>\n";
print IGV "<HTML lang=\"en\">\n";

print IGV "<HEAD>\n";
#print IGV "<meta charset=\"utf-8\"><meta http-equiv=\"X-UA-Compatible\" content=\"IE=edge\">\n";
#print IGV "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no\">\n";
#print IGV "<meta name=\"description\" content=\"\"> <meta name=\"author\" content=\"cbenner\">\n";
#print IGV "<link rel=\"shortcut icon\" href=\"https://igv.org/web/img/favicon.ico\">\n";

print IGV "  <title>$org library IGV Session</title>\n";
print IGV "  <meta charset=\"utf-8\">\n";
print IGV "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no\">\n";
print IGV "  <!-- IGV JS-->\n";
print IGV '  <script src="https://cdn.jsdelivr.net/npm/igv@2.7.4/dist/igv.min.js"></script>' . "\n";

#
#print IGV "<script type=\"text/javascript\" src=\"https://igv.org/web/test/dist/igv.min.js\"></script>\n";
#print IGV "<script type=\"text/javascript\" src=\"https://igv.org/web/release/current/dist/igv.min.js\"></script>\n";
#print IGV "<script type=\"text/javascript\" src=\"https://igv.org/web/release/current/dist/igv.min.js\"></script>\n";

print IGV "</head>\n";


print IGV "\n<BODY>\n";
print IGV "<div id=\"myDiv\">\n";
print IGV "<h1>$org library IGV sessions</h1></div>\n";
print IGV "<script type=\"text/javascript\">\n";
print IGV "\tdocument.addEventListener(\"DOMContentLoaded\", function () {\n";
print IGV "\t\tvar div = document.getElementById(\"myDiv\");\n";
print IGV "\t\toptions = {\n";
print IGV "\t\t\treference: {\n";
print IGV "\t\t\t\tid: \"$org\",\n";
print IGV "\t\t\t\tfastaURL: \"$genomeURL\",\n";
print IGV "\t\t\t\tindexURL: \"$genomeIndexURL\",\n";
print IGV "\t\t\t},\n";
print IGV "\t\t\ttracks: [\n";

if (-e "$orgDir/genes.gtf") {
	print IGV "\t\t\t\t{\n";
	print IGV "\t\t\t\t\tname: \"Genes\",\n";
	print IGV "\t\t\t\t\ttype: \"annotation\",\n";
	print IGV "\t\t\t\t\tformat: \"gtf\",\n";
	print IGV "\t\t\t\t\tsourceType: \"file\",\n";
	print IGV "\t\t\t\t\turl: \"$genesURL\",\n";
	print IGV "\t\t\t\t\tindexURL: \"$genesIndexURL\",\n";
	print IGV "\t\t\t\t\tdisplayMode: \"EXPANDED\",\n";
	print IGV "\t\t\t\t\tsearchable: \"true\"\n";
	print IGV "\t\t\t\t}";
}

my $autoscaleGroup = 1;
my $srcRNA = $speciesDir . "data/$org/igv/RNA.pos.bw";
if (-e $srcRNA && $includeRNAbydefault) {
	my $color = "0,100,0";
	my $src = $speciesDir . "data/$org/igv/RNA.pos.bw";
	my $dest = $baseDir . "/RNA.pos.bigWig";
	`cp $src $dest`;
	my $bwURL = $baseURL . "/RNA.pos.bigWig";
	print IGV ",\n\t\t\t\t{\n";
	print IGV "\t\t\t\t\turl: \"$bwURL\",\n";
	print IGV "\t\t\t\t\tlabel: \"RNA.pos\",\n";
	print IGV "\t\t\t\t\tcolor: \"rgb($color)\",\n";
	#print IGV "\t\t\t\t\tmin: \"0\",\n";
	print IGV "\t\t\t\t\tautoscaleGroup: \"$autoscaleGroup\",\n";
	print IGV "\t\t\t\t}";
	$src = $speciesDir . "data/$org/igv/RNA.neg.bw";
	$dest = $baseDir . "/RNA.neg.bw";
	`cp $src $dest`;
	$bwURL = $baseURL . "/RNA.neg.bw";
	print IGV ",\n\t\t\t\t{\n";
	print IGV "\t\t\t\t\turl: \"$bwURL\",\n";
	print IGV "\t\t\t\t\tlabel: \"RNA.neg\",\n";
	print IGV "\t\t\t\t\tcolor: \"rgb($color)\",\n";
	#print IGV "\t\t\t\t\tmax: \"0\",\n";
	print IGV "\t\t\t\t\tautoscaleGroup: \"$autoscaleGroup\",\n";
	print IGV "\t\t\t\t}";
	$autoscaleGroup++;
}


my $totalN = scalar(@data);
my $z = 0;

foreach(@data) {
	$z++;
	my $name = $_->{'n'};
	my $dir = $_->{'f'};
	my $type = $_->{'t'};
	my $color = $_->{'c'};
	my $params = $_->{'p'};
	my $sampleForce = $_->{'force'};
	if ($force) {
		$sampleForce = $force;
	}

	$name =~ s/\///;

	my $premadebw = "";
	my $premadebw2 = "";
	if ($type eq 'bw') {
		if ($dir =~ /\,/) {
			#this means it was a dual/stranded bw submission
			my @files = split /\,/,$dir;
			$premadebw = $files[0];
			$premadebw2 = $files[1];
			$type = 'bw2';
		} else {
			$premadebw = $dir;
		}
	}

	my $bwFile = "";
	my $bwFile2 = "";
	if ($type eq 'chip' || $type eq 'atac' || $type eq 'chipseq' || $type eq 'atacseq' || $type eq 'bw') {
		$bwFile = "$name.cov.bw";
		my $fragOpt = '';
		if ($chipFragLength > -1e7) {
			$fragOpt .= " -fragLength $chipFragLength ";
		}
		if (-e "$baseDir/$bwFile" && $sampleForce == 0) {
		} elsif ($type eq 'bw1') {
			`cp "$premadebw" "$baseDir/$bwFile"`;
		} else {
			`makeUCSCfile $dir $fragOpt -noheader $params $normOpts > $baseDir/$name.cov.bedGraph`;
			`bedGraphToBigWig "$baseDir/$name.cov.bedGraph" "$chromSizes" "$baseDir/$bwFile"`;
			`rm "$baseDir/$name.cov.bedGraph"`;
		}

		$bwURL = $baseURL . $bwFile;
		print IGV ",\n\t\t\t\t{\n";
		print IGV "\t\t\t\t\turl: \"$bwURL\",\n";
		print IGV "\t\t\t\t\tlabel: \"$name\",\n";
		print IGV "\t\t\t\t\tcolor: \"rgb($color)\",\n";
		print IGV "\t\t\t\t\tautoscaleGroup: \"$autoscaleGroup\",\n";
		print IGV "\t\t\t\t}";
		$autoscaleGroup++;
	} elsif ($type eq 'gro' || $type eq 'tss' || $type eq 'rna' || $type eq 'bw2') {
		my $opt = " -fragLength 1 -normLength 0 ";
		if ($type eq 'gro') {
			$opt = " -fragLength $defaultFragLength ";
		} elsif ($type eq 'rna') {
			$opt = " -fragLength given ";
		}
		$bwFile = "$name.pos.bw";
		if (-e "$baseDir/$bwFile" && $sampleForce == 0) {
		} elsif ($type eq 'bw2') {
			`cp "$premadebw" "$baseDir/$bwFile"`;
		} else {
			`makeUCSCfile $dir $opt -strand + -noheader $params $normOpts > $baseDir/$name.pos.bedGraph`;
			`bedGraphToBigWig "$baseDir/$name.pos.bedGraph" "$chromSizes" "$baseDir/$bwFile"`;
			`rm "$baseDir/$name.pos.bedGraph"`;
		}
		$bwFile2 = "$name.neg.bw";
		if (-e "$baseDir/$bwFile2" && $sampleForce == 0) {
		} elsif ($type eq 'bw2') {
			`cp "$premadebw2" "$baseDir/$bwFile2"`;
		} else {
			`makeUCSCfile $dir $opt -strand - -neg -noheader $params $normOpts > $baseDir/$name.neg.bedGraph`;
			`bedGraphToBigWig "$baseDir/$name.neg.bedGraph" "$chromSizes" "$baseDir/$bwFile2"`;
			`rm "$baseDir/$name.neg.bedGraph"`;
		}

		my $bwURL = $baseURL . $bwFile;
		print IGV ",\n\t\t\t\t{\n";
		print IGV "\t\t\t\t\turl: \"$bwURL\",\n";
		print IGV "\t\t\t\t\tlabel: \"$name.pos\",\n";
		print IGV "\t\t\t\t\tcolor: \"rgb($color)\",\n";
		print IGV "\t\t\t\t\tautoscaleGroup: \"$autoscaleGroup\",\n";
		print IGV "\t\t\t\t}";

		$bwURL = $baseURL . $bwFile2;
		print IGV ",\n\t\t\t\t{\n";
		print IGV "\t\t\t\t\turl: \"$bwURL\",\n";
		print IGV "\t\t\t\t\tlabel: \"$name.neg\",\n";
		print IGV "\t\t\t\t\tcolor: \"rgb($color)\",\n";
		print IGV "\t\t\t\t\tautoscaleGroup: \"$autoscaleGroup\",\n";
		print IGV "\t\t\t\t}";
		$autoscaleGroup++;
	}
}
$z=0;
foreach(@anns) {
	$z++;
	my $name = $_->{'n'};
	my $file = $_->{'f'};
	my $type = $_->{'t'};
	my $color = $_->{'c'};
	my $params = $_->{'p'};

	$name =~ s/\///;

	#need to sort/index file
	my $dest = $baseDir . "$name.$type";
	my $bwURL = $baseURL . "$name.$type";

	`sortBed -i "$file" > "$dest"`;
	`igvtools index "$dest"`;
	#`mv "$tmpFile.gtf" "$orgDir/genes.gtf"`;
	#`mv "$tmpFile.gtf.idx" "$orgDir/genes.gtf.idx"`;
	#my $dest = $baseDir. "$name.bed";
	#`cp "$file" "$dest"`;
	#my $bwURL = $baseURL . "$name.bed";
	print IGV ",\n\t\t\t\t{\n";
	print IGV "\t\t\t\t\ttype: \"annotation\",\n";
	print IGV "\t\t\t\t\tsourceType: \"file\",\n";
	print IGV "\t\t\t\t\tformat: \"$type\",\n";
	print IGV "\t\t\t\t\tlabel: \"$name\",\n";
	print IGV "\t\t\t\t\turl: \"$bwURL\",\n";
	print IGV "\t\t\t\t\tcolor: \"rgb($color)\",\n";
	#print IGV "\t\t\t\t\tdisplayMode: \"EXPANDED\"\n";
	print IGV "\t\t\t\t}";
}

print IGV "\n";
print IGV "\t\t\t]\n";
print IGV "\t\t};\n";
print IGV "\t\tbrowser = igv.createBrowser(div, options);\n";
print IGV "\t\t})\n";
print IGV "</script>\n";
print IGV "</body>\n";
print IGV "</html>\n";

close IGV;

print STDERR "\n\tURL: $baseURL\n\n";


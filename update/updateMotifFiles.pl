#!/usr/bin/perl -w

print STDERR "\tNo input required...\n";
print STDERR "\n\t!!! If wget has trouble downloading JASPAR matrices, you can download them yourself\n";
print STDERR "\tPlace them in the current directory, filenames:\n";
print STDERR "\t\tpfm_vertebrates.txt\n";
print STDERR "\t\tpfm_plants.txt\n";
print STDERR "\t\tpfm_insects.txt\n";
print STDERR "\t\tpfm_nematodes.txt\n";
print STDERR "\t\tpfm_fungi.txt\n";
print STDERR "\tThis program will also copy over the context of the 'motifs/' folder\n";
print STDERR "\tPausing for 5 seconds...\n";
if (@ARGV < 1) {
#	print STDERR "\tNo input required...\n";
}
`sleep 5`;

my $opt = " --no-check-certificate ";

my $activeMotifs = "../motifs/";
my $homerMotifs = "../data/knownTFs/";

`mkdir -p $homerMotifs`;
`mkdir -p $homerMotifs/all`;
`mkdir -p $homerMotifs/motifs`;
`cp -r $activeMotifs/* $homerMotifs/motifs/`;

my @groups = (
	"vertebrates",
	"insects",
	"worms",
	"plants",
	"yeast"
);


#my $jasparServer = "http://jaspar.binf.ku.dk/";
my $jasparServer = "http://jaspar.genereg.net/";

# Update JASPAR motifs
unless (-e "pfm_vertebrates.txt") {
	`wget $opt -O pfm_vertebrates.txt https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt`;
	`./motifs/parseJasparMatrix.pl pfm_vertebrates.txt > motifs/vertebrates/jaspar.motifs`;
	`rm -f pfm_vertebrates.txt`;
} else {
	`./motifs/parseJasparMatrix.pl pfm_vertebrates.txt > motifs/vertebrates/jaspar.motifs`;
}

unless (-e "pfm_plants.txt") {
	`wget $opt -O pfm_plants.txt https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_plants_non-redundant_pfms_jaspar.txt`;
	`./motifs/parseJasparMatrix.pl pfm_plants.txt > motifs/plants/jaspar.motifs`;
	`rm -f pfm_plants.txt`;
} else {
	`./motifs/parseJasparMatrix.pl pfm_plants.txt > motifs/plants/jaspar.motifs`;
}
	

unless (-e "pfm_insects.txt") {
	`wget $opt -O pfm_insects.txt https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_insects_non-redundant_pfms_jaspar.txt`;
	`./motifs/parseJasparMatrix.pl pfm_insects.txt > motifs/insects/jaspar.motifs`;
	`rm -f pfm_insects.txt`;
} else {
	`./motifs/parseJasparMatrix.pl pfm_insects.txt > motifs/insects/jaspar.motifs`;
}

unless (-e "pfm_nematodes.txt") {
	`wget $opt -O pfm_nematodes.txt https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_nematodes_non-redundant_pfms_jaspar.txt`;
	`./motifs/parseJasparMatrix.pl pfm_nematodes.txt > motifs/worms/jaspar.motifs`;
	`rm -f pfm_nematodes.txt`;
} else {
	`./motifs/parseJasparMatrix.pl pfm_nematodes.txt > motifs/worms/jaspar.motifs`;
}

unless (-e "pfm_fungi.txt") {
	`wget $opt -O pfm_fungi.txt https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_fungi_non-redundant_pfms_jaspar.txt`;
	`./motifs/parseJasparMatrix.pl pfm_fungi.txt > motifs/yeast/jaspar.motifs`;
	`rm -f pfm_fungi.txt`;
} else {
	`./motifs/parseJasparMatrix.pl pfm_fungi.txt > motifs/yeast/jaspar.motifs`;
}

# known motifs first (ones that have degeneracy thresholds)
foreach(@groups) {
	my $g = $_;
	`mkdir -p $homerMotifs/$g`;
	if ($g eq 'vertebrates') {
		`cat $homerMotifs/motifs/*.motif > $homerMotifs/$g/known.motifs`;
	} else {
		`cat $homerMotifs/motifs/$g/*.motif > $homerMotifs/$g/known.motifs`;
	}
	`cat $homerMotifs/$g/known.motifs motifs/common.motifs motifs/$g/*.motifs > $homerMotifs/$g/all.motifs`;
}
`cat $homerMotifs/motifs/*.motif $homerMotifs/motifs/*/*.motif > $homerMotifs/all/known.motifs`;


`cat $homerMotifs/all/known.motifs motifs/common.motifs motifs/*/*.motifs > $homerMotifs/all/all.motifs`;
`cp $homerMotifs/all/known.motifs $homerMotifs/known.motifs`;
`cp $homerMotifs/all/all.motifs $homerMotifs/all.motifs`;

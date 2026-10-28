
// Copyright 2009 - 2014 Christopher Benner <cbenner@salk.edu>
// 
// This file is part of HOMER
//
// HOMER is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// HOMER is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include "Motif2.h"
#include <stdarg.h>
#include <algorithm>
#include <sys/stat.h>
#include <cmath>
#include <atomic>
#include <sys/mman.h>
#include <cerrno>
#include <condition_variable>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

int Motif::alphaSize = DNA_ALPHA;
char* Motif::alphaIndex = NULL;

#define BUFFER 100000
#define WHITE_SPACE 7
void split2(char* string, char** cols, int &numCols, char delim) {
	cols[0] = string;
	numCols=1;	
	char delim2 = 0;
	if (delim == WHITE_SPACE) {
		delim = '\t';
		delim2 = 32;
	}
	int len = strlen(string);
	for (int i=0;i<len;i++) {
		if (string[i] == delim || string[i] == delim2) {
			string[i] = '\0';
			cols[numCols] = &(string[i+1]);
			numCols++;
		} else if (string[i] == '\n') {
			string[i] = '\0';
		} else if (string[i] == '\r') {
			string[i] = '\0';
		}
	}
}




void printCMDhomer() {
	fprintf(stderr, "\n\thomer2 <command> [options]\n");
	fprintf(stderr, "\n\tCommands:\n");
	fprintf(stderr, "\t\tdenovo (discover motifs using de novo motif analysis)\n");
	fprintf(stderr, "\t\tknown (score known motifs for enrichment)\n");
	fprintf(stderr, "\t\tfind (find specific instances motifs in the data)\n");
	fprintf(stderr, "\t\tbackground (select background sequences, can also specify 'bg')\n");
	fprintf(stderr, "\t\tnorm (normalize background sequences)\n");
	fprintf(stderr, "\t\tmask (mask instances of motifs)\n");
	fprintf(stderr, "\n");
	exit(0);
}

void printCMDdenovo() {
	fprintf(stderr, "\n\thomer2 denovo -s <seq file> -g <group file> [options]\n");
	fprintf(stderr, "\n\tDiscover motif de novo in a set of sequences.\n");
	fprintf(stderr, "\n\tPrimary Inputs:\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t-g <group file> (sequence group and weight assignments)\n");
	fprintf(stderr, "\t\t\t--- or ---\n");
	fprintf(stderr, "\t\t-i <input FASTA file>\n");
	fprintf(stderr, "\t\t-b <background FASTA file>\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-o <output motif file> (default: sent to stdout)\n");
	//fprintf(stderr, "\t\t-sites <predicted sites file> (default: none)\n");
	fprintf(stderr, "\t\t-len <#> (length of motif to search for, default: 10)\n");
	fprintf(stderr, "\t\t-mis <#> (maximum number of mismatches in global search phase, default: 2)\n");
	fprintf(stderr, "\t\t-strand <+|-|both> (search for motifs on specific strand, default: both)\n");
	fprintf(stderr, "\t\t-stat <hypergeo|binomial> (enrichment statistic, default: binomial)\n");
	//fprintf(stderr, "\t\t-stat <hypergeo|binomial|zscore> (enrichment statistic, default: binomial)\n");
	fprintf(stderr, "\t\t-S <#> (Total motifs to find, default: %d)\n",MOTIF_DEFAULT_NUMBER_OF_MOTIFS);
	fprintf(stderr, "\t\t-olen <#> (length of lower-order oligos to normalize in oligo table)\n");
	fprintf(stderr, "\t\t\t-oout <filename> (output normalization weights to file)\n");
	fprintf(stderr, "\t\t\t-omax <#> (max oligo normalization iterations, default: %d)\n", 
																MAX_OLIGO_NORMALIZATION_ITERATIONS);
	fprintf(stderr, "\t\t-p <#> (Number of processers to use, default: 1)\n");
	fprintf(stderr, "\n\tAdvanced Options:\n");
	fprintf(stderr, "\t\t-nozoops (skip proper zoops scoring at end, default: -zoops)\n");
	fprintf(stderr, "\t\t-tmp <filename> (temporary results file, default: \".tmp.motifs\")\n");
	fprintf(stderr, "\t\t-oligos <filename> (print enrichment of individual oligos)\n");
	fprintf(stderr, "\t\t-opt <motif file> (expand/futher optimize these motifs/skip global phase)\n");

	fprintf(stderr, "\n\t\tSpeed vs. Sensitivity:\n");
	fprintf(stderr, "\t\t-fullMask | -quickMask (choose one:)\n");
	fprintf(stderr, "\t\t\t-fullMask (as motifs are found, mask them from original sequences.\n");
	fprintf(stderr, "\t\t\t\tRequires more memory but yields \"cleaner\" results, default)\n");
	fprintf(stderr, "\t\t\t-quickMask (as motifs are found, mask bound \"oligos\" only, old way)\n");
	fprintf(stderr, "\t\t-e <#> (maximum expected motif instances per bp, default: %.3lf)\n",
																	MOTIF_MAX_HIT_RATE_TO_CONSIDER);
	fprintf(stderr, "\t\t-T <#> (number of trial matrices when optimizing, default: %d)\n",
																	MOTIFS_NUM_POTENTIAL_MOTIFS);
	fprintf(stderr, "\t\t-Tcpus <#> (round -T to a multiple of this many CPUs instead of -p, default: -p)\n");
	fprintf(stderr, "\t\t-blen <#> (# of bp on either side to check for redundancy, default: 1)\n");
	fprintf(stderr, "\t\t-maxBack <#> (Max percentage of background that motifs may contain, default: 0.5)\n");
	fprintf(stderr, "\t\t-minlp <#> (minimum significance of seeds to optimize, default: %.3lf)\n",
																	MOTIF_SEED_MINIMUM_LOGP);

	fprintf(stderr, "\n\t\tMemory vs. Speed:\n");
	fprintf(stderr, "\t\t-cache <#> (size in MB of stat cache, helps for hypergeo, default: %d)\n",
																	MAX_STAT_CACHE_SIZE);
																	
	fprintf(stderr, "\n");
	exit(0);
}

void programDeNovo(int argc, char** argv) {

	DeNovoMotifs* denovo = new DeNovoMotifs();


	FILE* oligoOutputFile = NULL;
	FILE* outputFp = NULL;
	char* outputFileName = NULL;
	
	if (argc < 3) {
		printCMDdenovo();
	}	
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-s")==0) {
			denovo->seqFile = argv[++i];
		} else if (strcmp(argv[i],"-g")==0) {
			denovo->groupFile = argv[++i];
		} else if (strcmp(argv[i],"-i")==0) {
			denovo->inputFASTAfile = argv[++i];
		} else if (strcmp(argv[i],"-b")==0) {
			denovo->bgFASTAfile = argv[++i];
		} else if (strcmp(argv[i],"-sites")==0) {
			denovo->sitefp = fopen(argv[++i],"w");
			if (denovo->sitefp == NULL) {
				fprintf(stderr, "Could not open %s for writing sites!!!\n",argv[i]);
				exit(0);
			}
		} else if (strcmp(argv[i],"-len")==0) {
			sscanf(argv[++i],"%d", &(denovo->oligoLength));
		} else if (strcmp(argv[i],"-olen")==0) {
			sscanf(argv[++i],"%d", &(denovo->normLength));
		} else if (strcmp(argv[i],"-omax")==0) {
			sscanf(argv[++i],"%d", &(denovo->maxAutoNormIters));
		} else if (strcmp(argv[i],"-T")==0) {
			sscanf(argv[++i],"%d", &(denovo->numTrialMotifs));
		} else if (strcmp(argv[i],"-Tcpus")==0) {
			sscanf(argv[++i],"%d", &(denovo->trialCPUs));
		} else if (strcmp(argv[i],"-minlp")==0) {
			sscanf(argv[++i],"%lf", &(denovo->minimumSeedLogp));
		} else if (strcmp(argv[i],"-p")==0) {
			sscanf(argv[++i],"%d", &(denovo->numCPUs));
		} else if (strcmp(argv[i],"-blen")==0) {
			sscanf(argv[++i],"%d", &(denovo->bufferLength));
		} else if (strcmp(argv[i],"-cache")==0) {
			unsigned long long int sizeMB = 0;
			sscanf(argv[++i],"%lld", &sizeMB);
			denovo->maxCacheSize = sizeMB*1000000/sizeof(double);
		} else if (strcmp(argv[i],"-zoops")==0) {
			denovo->finalZoopsFlag=1;
		} else if (strcmp(argv[i],"-nozoops")==0) {
			denovo->finalZoopsFlag=0;
		} else if (strcmp(argv[i],"-rand")==0) {
			denovo->randFlag=1;
		} else if (strcmp(argv[i],"-o")==0) {
			outputFileName=argv[++i];
			outputFp=fopen(outputFileName,"w");
			if (outputFp==NULL) {
				fprintf(stderr, "!!! Could not open %s for writing output!!!\n",outputFileName);
				exit(1);
			}
		} else if (strcmp(argv[i],"-oout")==0) {
			denovo->normFile = argv[++i];
		} else if (strcmp(argv[i],"-tmp")==0) {
			denovo->tmpFile = argv[++i];
		} else if (strcmp(argv[i],"-opt")==0) {
			denovo->optimizeGivenFlag=1;
			fprintf(stderr, "\tOptimizing given motifs...\n");
			denovo->motifs = Motif::readMotifFile(argv[++i],denovo->numMotifs);
		} else if (strcmp(argv[i],"-mis")==0) {
			sscanf(argv[++i],"%d", &(denovo->maxMisMatches));
		} else if (strcmp(argv[i],"-maxBack")==0) {
			sscanf(argv[++i],"%lf", &(denovo->maxBackgroundPercent));
		} else if (strcmp(argv[i],"-fullMask")==0) {
			denovo->maskMethod = DENOVO_MASK_FULL;
		} else if (strcmp(argv[i],"-quickMask")==0) {
			denovo->maskMethod = DENOVO_MASK_QUICK;
		} else if (strcmp(argv[i],"-e")==0) {
			sscanf(argv[++i],"%lf", &(denovo->maximumExpectedPerBp));
		} else if (strcmp(argv[i],"-S")==0) {
			sscanf(argv[++i],"%d", &(denovo->totalMotifsToFind));
		} else if (strcmp(argv[i],"-oligos")==0) {
			oligoOutputFile = fopen(argv[++i],"w");
			if (oligoOutputFile == NULL) {
				fprintf(stderr, "!!! Could not open oligo output file (%s)!!!\n", argv[i]);
				exit(0);
			}
		} else if (strcmp(argv[i],"-speed")==0) {
			i++;
			if (strcmp(argv[i],"auto")==0) {
				denovo->speedFlag = DENOVO_SPEED_AUTOADJUST;
			} else if (strcmp(argv[i],"constant")==0) {
				denovo->speedFlag = DENOVO_SPEED_CONSTANT;
			} else {
				fprintf(stderr, "!!! option \"-speed %s\" not valid !!!\n", argv[i]);
				exit(0);
			}
		} else if (strcmp(argv[i],"-stat")==0) {
			i++;
			if (strcmp(argv[i],"hypergeo")==0) {
				denovo->scoringMethod = MOTIF_SCORING_HYPERGEOMETRIC;
			} else if (strcmp(argv[i],"binomial")==0) {
				denovo->scoringMethod = MOTIF_SCORING_BINOMIAL;
			} else if (strcmp(argv[i],"zscore")==0) {
				denovo->scoringMethod = MOTIF_SCORING_ZSCORE;
			} else {
				fprintf(stderr, "!!! option \"-stat %s\" not valid !!!\n", argv[i]);
				exit(0);
			}
		} else if (strcmp(argv[i],"-strand")==0) {
			i++;
			if (strcmp(argv[i],"+")==0) {
				denovo->strand = MOTIF_STRAND_POS;
			} else if (strcmp(argv[i],"-")==0) {
				denovo->strand = MOTIF_STRAND_NEG;
			} else if (strcmp(argv[i],"both")==0) {
				denovo->strand = MOTIF_STRAND_BOTH;
			} else {
				fprintf(stderr, "!!! Could not recognize -strand \"%s\" !!!\n", argv[i]);
				printCMDdenovo();
			}
		} else {
			printCMDdenovo();
		}
	}
	fprintf(stderr, "\n");
	
	if (outputFileName == NULL) outputFp = stdout;

	denovo->checkParameters();
	if (denovo->fastaFlag) {
		denovo->parseFasta();
	}

	denovo->readSequence2Tree();
	denovo->normalizeOligos();

	//setup for motif finding
	denovo->initializeEnrichmentScoring();

	if (denovo->optimizeGivenFlag == 0) {
		//normal operation
		//Global optimization
		denovo->oligoArray->initializeActiveOligoArray();
		denovo->oligoArray->initializeOligoPvalues(denovo->scorer);
		denovo->findMisMatchEnrichment(denovo->maxMisMatches, oligoOutputFile,HOMER_VERBOSE_LEVEL1);	
		//local optimization
		if (denovo->maskMethod == DENOVO_MASK_FULL) {
			denovo->loadFullSequences();
		}
		denovo->optimizeSeeds2Motifs();
	} else {
		for (int i=0;i<denovo->numMotifs;i++) denovo->motifs[i]->adjustLength(denovo->oligoLength);
		denovo->optimizeGivenMotifs();
	}

	denovo->freeOligos();

	if (denovo->finalZoopsFlag && denovo->maskMethod == DENOVO_MASK_QUICK) {
		fprintf(stderr, "\n\tFinalizing enrichment with proper zoops scoring...\n");
		denovo->scoreWithKnownMotifs(KNOWN_OPTIMIZE_THREASHOLD);
	} else {
		fprintf(stderr, "\n\tFinalizing Enrichment Statistics (new in v3.4)\n");
		denovo->scoreWithKnownMotifs(KNOWN_JUST_SCORE_ZOOPS);
	}
	denovo->printMotifs(outputFp);

	if (denovo->tmpFile != NULL) {
		(void)unlink(denovo->tmpFile);
	} else {	
		(void)unlink(".tmp.motifs");
	}
	if (outputFileName != NULL) {
		fclose(outputFp);	
		fprintf(stderr, "\tOutput in file: %s\n", outputFileName);
	}
	fprintf(stderr, "\n");


	if (denovo->fastaFlag) {
		fprintf(stderr, "\tCleaning up temporary sequence and group files:\n");
		fprintf(stderr, "\t\t%s\n",denovo->seqFile);
		fprintf(stderr, "\t\t%s\n",denovo->groupFile);
		(void)unlink(denovo->seqFile);
		(void)unlink(denovo->groupFile);
	}

	delete denovo;

}


void printCMDknown() {
	fprintf(stderr, "\n\thomer2 known -s <seq file> -g <group file> -k <known motif file> [options]\n");
	fprintf(stderr, "\n\tFind the enrichment of known motifs in a set of sequences.\n");
	fprintf(stderr, "\n\tPrimary Inputs:\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t-g <group file> (sequence group and weight assignments)\n");
	fprintf(stderr, "\t\t\t--- or ---\n");
	fprintf(stderr, "\t\t-i <input FASTA file>\n");
	fprintf(stderr, "\t\t-b <background FASTA file>\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-o <output enrichment file> (default: sent to stdout)\n");
	fprintf(stderr, "\t\t-m <known motif file> (Known motif file to check enrichment for)\n");
	fprintf(stderr, "\t\t\t-mout <output motif file> (Output updated motifs with statistics)\n");
	fprintf(stderr, "\t\t-strand <+|-|both> (search for motifs on specific strand, default: both)\n");
	fprintf(stderr, "\t\t-stat <hypergeo|binomial> (enrichment statistic, default: binomial)\n");
	fprintf(stderr, "\t\t-nlen <#> (length of lower-order oligos to normalize, default: 0)\n");
	fprintf(stderr, "\t\t\t-nout <filename> (output normalization weights to file)\n");
	fprintf(stderr, "\t\t\t-nmax <#> (max normalization iterations, default: %d\n", 
																MAX_OLIGO_NORMALIZATION_ITERATIONS);
	fprintf(stderr, "\t\t-cache <#> (size in MB of stat cache: %d)\n",MAX_STAT_CACHE_SIZE);
	fprintf(stderr, "\t\t-p <#> (Number of processers to use, default: 1)\n");
	fprintf(stderr, "\t\t-opt (Optimize degeneracy threshold to get best enrichment, use -mout to get motifs)\n");
	fprintf(stderr, "\t\t-siteReduce <#> (Eliminate redundant motifs sharing > %% of sites)\n");
	fprintf(stderr, "\t\t-maxBack <#> (Max percentage of background that motifs may contain, default: 0.5)\n");
	fprintf(stderr, "\n");
	exit(0);
}
void printCMDfind() {
	fprintf(stderr, "\n\thomer2 find -s <seq file> -m <motif file> [options]\n");
	fprintf(stderr, "\n\tFind the instances of motif(s) in a set of sequences.\n");
	fprintf(stderr, "\n\tPrimary Inputs:\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t\t--- or ---\n");
	fprintf(stderr, "\t\t-i <input FASTA file>\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-o <output file> (default: sent to stdout)\n");
	fprintf(stderr, "\t\t-m <motif file> (Motif(s) to find instances of)\n");
	fprintf(stderr, "\t\t-offset <#> (offset to report motif instances from, default: midpoint)\n");
	fprintf(stderr, "\t\t-strand <+|-|both> (search for motifs on specific strand, default: both)\n");
	fprintf(stderr, "\t\t-p <#> (Number of processers to use, default: 1)\n");
	fprintf(stderr, "\t\t-mscore (instead of reporting sites, report best motif score per sequence)\n");
	fprintf(stderr, "\n");
	exit(0);
}
void printCMDmask() {
	fprintf(stderr, "\n\thomer2 mask -s <seq file> -m <motif file> [options]\n");
	fprintf(stderr, "\n\tRemove instances of motif(s) in a set of sequences.\n");
	fprintf(stderr, "\n\tPrimary Inputs:\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t\t--- or ---\n");
	fprintf(stderr, "\t\t-i <input FASTA file>\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-o <output tsv sequence file> (default: sent to stdout)\n");
	fprintf(stderr, "\t\t-m <motif file> (Motif(s) to find instances of)\n");
	fprintf(stderr, "\t\t-strand <+|-|both> (search for motifs on specific strand, default: both)\n");
	fprintf(stderr, "\t\t-p <#> (Number of processers to use, default: 1)\n");
	fprintf(stderr, "\n");
	exit(0);
}
void printCMDnorm() {
	fprintf(stderr, "\n\thomer2 norm -s <seq file> -g <group file> [options]\n");
	fprintf(stderr, "\n\tNormalize background sequences to remove short oligo enrichment\n");
	fprintf(stderr, "\n\tCreates new group file to use with motif finding\n");
	fprintf(stderr, "\n\tPrimary Inputs:\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t-g <group file> (sequence group and weight assignments)\n");
	fprintf(stderr, "\t\t\t--- or ---\n");
	fprintf(stderr, "\t\t-i <input FASTA file>\n");
	fprintf(stderr, "\t\t-b <background FASTA file>\n");
	fprintf(stderr, "\t\t\tOutput will be sequence and group file, even with FASTA input files\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-o <output weighted group file> (default: sent to stdout)\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t-g <group file> (sequence group and weight assignments)\n");
	fprintf(stderr, "\t\t-strand <+|-|both> (search for motifs on specific strand, default: both)\n");
	fprintf(stderr, "\t\t-nlen <#> (length of lower-order oligos to normalize, default: 0)\n");
	fprintf(stderr, "\t\t\t-nout <filename> (output normalization weights to file)\n");
	fprintf(stderr, "\t\t\t-nmax <#> (max normalization iterations, default: %d)\n", 
																MAX_OLIGO_NORMALIZATION_ITERATIONS);
	fprintf(stderr, "\t\t-neutral (set target/background to neutral i.e. 25%%, 6.25%%, etc. frequencies)\n");
	fprintf(stderr, "\t\t-p <#> (Number of processers to use, default: 1)\n");
	fprintf(stderr, "\n");
	exit(0);
}




void programKnown(int argc, char** argv) {

	KnownMotifs* known = new KnownMotifs();

	FILE* outputFp = NULL;
	char* outputFileName = NULL;

	known->findFlag = 0;
	if (strcmp(argv[1],"find")==0) {
		known->findFlag =1;
		known->knownFlag =0;
	}
	known->maskFlag = 0;
	if (strcmp(argv[1],"mask")==0) {
		known->maskFlag =1;
		known->knownFlag =0;
	}
	known->normOnlyFlag = 0;
	if (strcmp(argv[1],"norm")==0) {
		known->normOnlyFlag =1;
		known->knownFlag =0;
	}
	
	if (argc < 3) {
		if (known->maskFlag) printCMDmask();
		if (known->findFlag) printCMDfind();
		if (known->normOnlyFlag) printCMDnorm();
		printCMDknown();
	}	
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-s")==0) {
			known->seqFile = argv[++i];
		} else if (strcmp(argv[i],"-g")==0) {
			known->groupFile = argv[++i];
		} else if (strcmp(argv[i],"-i")==0) {
			known->inputFASTAfile = argv[++i];
		} else if (strcmp(argv[i],"-b")==0) {
			known->bgFASTAfile = argv[++i];
		} else if (strcmp(argv[i],"-m")==0) {
			known->knownFile = argv[++i];
		} else if (strcmp(argv[i],"-p")==0) {
			sscanf(argv[++i],"%d", &(known->numCPUs));
		} else if (strcmp(argv[i],"-opt")==0) {
			known->threasholdOptimizationFlag = KNOWN_OPTIMIZE_THREASHOLD;
		} else if (strcmp(argv[i],"-mout")==0) {
			known->motifFilename = argv[++i];
		} else if (strcmp(argv[i],"-mscore")==0) {
			known->mscoreFlag = 1;
		} else if (strcmp(argv[i],"-neutral")==0) {
			known->neutralFlag = 1;
		} else if (strcmp(argv[i],"-siteReduce")==0) {
			sscanf(argv[++i],"%lf", &(known->siteReduceThreshold));
		} else if (strcmp(argv[i],"-maxBack")==0) {
			sscanf(argv[++i],"%lf", &(known->maxBackgroundPercent));
		} else if (strcmp(argv[i],"-offset")==0) {
			sscanf(argv[++i],"%d", &(known->offset));
		} else if (strcmp(argv[i],"-nlen")==0) {
			sscanf(argv[++i],"%d", &(known->normLength));
		} else if (strcmp(argv[i],"-nmax")==0) {
			sscanf(argv[++i],"%d", &(known->maxAutoNormIters));
		} else if (strcmp(argv[i],"-cache")==0) {
			unsigned long long int sizeMB = 0;
			sscanf(argv[++i],"%lld", &sizeMB);
			known->maxCacheSize = sizeMB*1000000/sizeof(double);
		} else if (strcmp(argv[i],"-o")==0) {
			outputFileName=argv[++i];
			outputFp=fopen(outputFileName,"w");
			if (outputFp==NULL) {
				fprintf(stderr, "!!! Could not open %s for writing output!!!\n",outputFileName);
				exit(1);
			}
		} else if (strcmp(argv[i],"-nout")==0) {
			known->normFile = argv[++i];
		} else if (strcmp(argv[i],"-stat")==0) {
			i++;
			if (strcmp(argv[i],"hypergeo")==0) {
				known->scoringMethod = MOTIF_SCORING_HYPERGEOMETRIC;
			} else if (strcmp(argv[i],"binomial")==0) {
				known->scoringMethod = MOTIF_SCORING_BINOMIAL;
			}
		} else if (strcmp(argv[i],"-strand")==0) {
			i++;
			if (strcmp(argv[i],"+")==0) {
				known->strand = MOTIF_STRAND_POS;
			} else if (strcmp(argv[i],"-")==0) {
				known->strand = MOTIF_STRAND_NEG;
			} else if (strcmp(argv[i],"both")==0) {
				known->strand = MOTIF_STRAND_BOTH;
			} else {
				fprintf(stderr, "!!! Could not recognize -strand \"%s\" !!!\n", argv[i]);
				if (known->maskFlag) printCMDmask();
				if (known->findFlag) printCMDfind();
				if (known->normOnlyFlag) printCMDnorm();
				printCMDknown();
			}
		} else {
			if (known->maskFlag) printCMDmask();
			if (known->findFlag) printCMDfind();
			if (known->normOnlyFlag) printCMDnorm();
			printCMDknown();
		}
	}
	fprintf(stderr, "\n");

	if (outputFileName == NULL) outputFp = stdout;
	if (known->findFlag) known->sitefp = outputFp;
	if (known->normOnlyFlag) known->normfp = outputFp;

	known->checkParameters();
	if (known->fastaFlag) {
		known->parseFasta();
	}
	known->loadSequence();

	if (known->normOnlyFlag == 0) {
		known->loadMotifs();
		known->getMotifEnrichment();
		if (known->knownFlag) {
			known->printMotifEnrichment(outputFp);
		}
		if (known->siteReduceThreshold > 0) known->findOverlappingMotifs();

		if (known->motifFilename != NULL) {
			FILE* mfp = fopen(known->motifFilename,"w");
			if (mfp != NULL) {
				known->printMotifs(mfp);
				fclose(mfp);
			}
		}
		if (known->maskFlag) {
			known->seqArray->printSequenceFile(outputFp);
		}
	}
	if (outputFileName != NULL) fclose(outputFp);

	if (known->fastaFlag) {
		fprintf(stderr, "\tCleaning up temporary sequence and group files:\n");
		fprintf(stderr, "\t\t%s\n",known->seqFile);
		fprintf(stderr, "\t\t%s\n",known->groupFile);
		(void)unlink(known->seqFile);
		(void)unlink(known->groupFile);
	}
}


void splitMotif2(char* string, char** cols, int &numCols, char delim) {
	cols[0] = string;
	numCols=1;
	int len = strlen(string);
	for (int i=0;i<len;i++) {
		if (string[i] == delim) {
			string[i] = '\0';
			cols[numCols] = &(string[i+1]);
			numCols++;
		} else if (string[i] == '\n') {
			string[i] = '\0';
		} else if (string[i] == '\r') {
			string[i] = '\0';
		}
	}
}


DeNovoMotifs::DeNovoMotifs() {
	init();
}
DeNovoMotifs::~DeNovoMotifs() {
	if (oligoArray != NULL) {
		delete oligoArray;
	}
	oligoArray = NULL;
	if (seedOligos != NULL) delete []seedOligos;
	if (scorer != NULL) delete scorer;
	if (fullSequences != NULL) delete fullSequences;
	// fill this in later
}
void DeNovoMotifs::freeOligos() {
	if (oligoArray != NULL) delete oligoArray;
	oligoArray = NULL;
}
void DeNovoMotifs::init() {
	seqFile = NULL;
	groupFile = NULL;
	inputFASTAfile = NULL;
	bgFASTAfile = NULL;
	numTargets = 0.0;
	numBackground = 0.0;
	maxMisMatches = 2;
	oligoArray = NULL;
	alphaSize = DNA_ALPHA;
	totalOligos = 0;
	strand = MOTIF_STRAND_BOTH;
	oligoLength = 10;
	normLength = 0;
	noDuplicateOligos=1;

	maskMethod = DENOVO_MASK_FULL;
	fullSequences = NULL;
	zoopsFlag = SCORER_APPROXIMATE_ZOOPS;
	seedOligos=NULL;
	normFile = NULL;
	tmpFile = NULL;
	sitefp = NULL;
	numSeedOligos = 0;
	bufferLength = -1;
	numOligoEstimate = 0;
	numCPUs=1;
	fastaFlag = 0;
	randFlag = 0;
	maximumExpectedPerBp = -1.0;
	maxCacheSize = MAX_STAT_CACHE_SIZE*1000000/sizeof(double);
	minimumSeedLogp = MOTIF_SEED_MINIMUM_LOGP;

	finalZoopsFlag = 1;
	scorer = NULL;
	siteReduceThreshold = MOTIF_DEFAULT_REDUCE_PERCENTAGE;

	optimizeGivenFlag = 0;
	maxBackgroundPercent = MOTIF_MAX_BACKGROUND_PERCENT;

	scoringMethod = MOTIF_SCORING_BINOMIAL;
	misMatchTargetOnlyFlag = 1;
	misMatchMinPvalue = MOTIF_MISMATCH_PVALUE_MINIMUM;
	minimumMisMatchToStartSkipping = MOTIF_MISMATCH_MINIMUM_SKIP;
	minimumCountsToStartSkipping = MOTIF_MISMATCH_MINIMUM_COUNTS;
	minimumFractionToStartSkipping = 0.000; // not used

	totalMotifsToFind = MOTIF_DEFAULT_NUMBER_OF_MOTIFS;
	numTrialMotifs = MOTIFS_NUM_POTENTIAL_MOTIFS;
	trialCPUs = 0;
	maxOptimizationIterations = MOTIFS_MAX_OPT_ITERATIONS;
	motifOptimizationMethod = MOTIF_OPTIMZATION_DIFFERENT_THRESHOLDS;
	
	maxAutoNormIters = MAX_OLIGO_NORMALIZATION_ITERATIONS;
	
}


void DeNovoMotifs::checkParameters() {
	if (numCPUs < 1) numCPUs = 1;
	// -Tcpus: the trial threads take matrices from a shared queue, so any -p gives
	// the same motifs for a given -T; this keeps the -T that -p <trialCPUs> would use.
	int tcpus = (trialCPUs > 0) ? trialCPUs : (int)numCPUs;
	if (tcpus > 1) {
		int og = numTrialMotifs;
		int trialsPerCPU = ((int)(((((double)numTrialMotifs)/((double)tcpus)))+0.5));
		if (trialsPerCPU < 1) trialsPerCPU = 1;
		numTrialMotifs = trialsPerCPU*tcpus;
		if (og != numTrialMotifs) {
			fprintf(stderr, "\tNumber of Trial motifs (-T) set to %d (from %d) to work well with %d CPUs\n",
							numTrialMotifs, og, tcpus);
		}
	}
	if (inputFASTAfile == NULL || bgFASTAfile == NULL) {
		if (seqFile == NULL ){
			fprintf(stderr, "!!! Need some type of sequence file for homer to do anything !!!\n");
			exit(0);
		}
		if (groupFile == NULL ){
			fprintf(stderr, "!!! Need a group file for homer to do anything !!!\n");
			exit(0);
		}
	} else {
		if (inputFASTAfile == NULL) {
			fprintf(stderr, "!!! Missing input FASTA file (-i) !!!\n");
			exit(0);
		}
		if (bgFASTAfile == NULL) {
			fprintf(stderr, "!!! Missing background FASTA file (-b) !!!\n");
			exit(0);
		}
		fastaFlag = 1;
		fprintf(stderr, "\tTreating input files as FASTA format\n");
	}
	if (bufferLength == -1) {
		if (oligoLength < 8) {
			bufferLength = 1;
		} else {
			bufferLength = 1 + (oligoLength-8)/4;
		}
		if (bufferLength > 1) {
			fprintf(stderr, "\t-blen automatically set to %d\n", bufferLength);
		}
		if (bufferLength > 1000) bufferLength = 0;
	}
	if (maximumExpectedPerBp < 0.0) {
		maximumExpectedPerBp = MOTIF_MAX_HIT_RATE_TO_CONSIDER;
		if (oligoLength == 10) {
			maximumExpectedPerBp = MOTIF_MAX_HIT_RATE_TO_CONSIDER*2;
		} else if (oligoLength == 9) {
			maximumExpectedPerBp = MOTIF_MAX_HIT_RATE_TO_CONSIDER*4;
		} else if (oligoLength < 9) {
			maximumExpectedPerBp = MOTIF_MAX_HIT_RATE_TO_CONSIDER*8;
		}
	}
}
void DeNovoMotifs::parseFasta() {
	SequenceArray* seqarray = new SequenceArray();
	seqFile = new char[10000];
	groupFile = new char[10000];

	struct timeval tv;
	gettimeofday(&tv,NULL);
	unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
	srand(time_in_micros);
	//fprintf(stderr, "seed=%ld\n", time_in_micros);

	int r = rand();
	sprintf(seqFile, "0%d.seq", r);
	sprintf(groupFile, "0%d.group", r);
	FILE* fptest = fopen(seqFile,"r");
	while (fptest != NULL) {
		fclose(fptest);
		r = rand();
		sprintf(seqFile, "0%d.seq", r);
		sprintf(groupFile, "0%d.group", r);
		fptest = fopen(seqFile,"r");
	}
	seqarray->parseFasta2SeqAndGroupFiles(inputFASTAfile,bgFASTAfile,seqFile,groupFile);
	delete seqarray;
}


void DeNovoMotifs::normalizeOligos() {
	oligoArray->maxAutoNormIters = maxAutoNormIters;
	if (normLength > 0) {
		oligoArray->normalizeOligos(normLength,normFile);
	}
}

void DeNovoMotifs::readSequence2Tree() {

	Hashtable* groupdata = new Hashtable();

	fprintf(stderr, "\tScanning input files...\n");
	
	FILE* fp = fopen(groupFile, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open group file: %s !!!\n", groupFile);
		exit(0);
	}
	char* buf = new char[MOTIF2_BUFFER];
	buf[MOTIF2_BUFFER-1]='\0';
	char** cols = new char*[1000];
	int numCols;
	char* seq = NULL;
	char* name = NULL;

	int totalTargetsGroup = 0;
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		if (buf[MOTIF2_BUFFER-1]!='\0') {
			fprintf(stderr, "!!! Input sequence is too long (>%d)!!!\n", MOTIF2_BUFFER);
		}
		buf[MOTIF2_BUFFER-1]='\0';
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 2) continue;
		GroupInfo* gi = new GroupInfo();
		gi->weight = 1.0;
		gi->group = 0;
		name = cols[0];
		if (numCols > 1) {
			sscanf(cols[1],"%d",&(gi->group));
			if (gi->group) totalTargetsGroup++;
		}
		if (numCols > 2) {
			sscanf(cols[2],"%lf",&(gi->weight));
		}
		groupdata->insert(gi, name);
	}
	fclose(fp);

	if (randFlag) {
		//randomize group labels
		double targetFrac = ((double)totalTargetsGroup)/((double)groupdata->total);
		srand(time(NULL));
		char** keys = groupdata->keys();
		for (int i=0;i<groupdata->total;i++) {
			GroupInfo* gi = (GroupInfo*)groupdata->search(keys[i]);
			if (((double)rand())/((double)RAND_MAX) < targetFrac) {
				gi->group = 1;
			} else {
				gi->group = 0;
			}
			delete [](keys[i]);
		}
		delete []keys;
	}


	fp = fopen(seqFile, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open seq file: %s !!!\n", seqFile);
		exit(0);
	}

	numOligoEstimate = 0;
	int totalSequences = 0;
	long long int totalLength = 0;
	buf[MOTIF2_BUFFER-1]='\0';
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		if (buf[MOTIF2_BUFFER-1]!='\0') {
			fprintf(stderr, "!!! Input sequence is too long (>%d)!!!\n", MOTIF2_BUFFER);
		}
		buf[MOTIF2_BUFFER-1]='\0';
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 2) continue;
		name = cols[0];
		seq = cols[1];
//fprintf(stderr, " %s",name);
		GroupInfo* gi = (GroupInfo*)groupdata->search(name);
		if (gi == NULL) {
			continue;
		}
		totalSequences++;
		int slen = strlen(seq);
		int contribution = slen-oligoLength+1;
		if (contribution < 0) contribution = 0;
		numOligoEstimate += contribution;
		totalLength += slen;
	}
	if (numOligoEstimate <= 0) {
		fprintf(stderr, "!!! Something is wrong... are you sure you chose the right length for motif finding?\n");
		fprintf(stderr, "!!! i.e. also check your sequence file!!!\n");
		exit(0);
	}
	double avgSeqLength = totalLength/totalSequences;

	oligoArray = new OligoArray();
	oligoArray->init(oligoLength,numOligoEstimate, alphaSize, strand, numCPUs);
	double totalTargetOligos = 0.0;
	double totalBackgroundOligos = 0.0;

	rewind(fp);	

	//now read in for real
	int curNumSequences = 0;
	int barLength = totalSequences/80;
	int curThreshold = barLength;
	fprintf(stderr, "\tParsing sequences...\n");
	fprintf(stderr, "\t|0%%                                   50%%                                  100%%|\n");
	fprintf(stderr, "\t");
	buf[MOTIF2_BUFFER-1]='\0';
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		if (buf[MOTIF2_BUFFER-1]!='\0') {
			fprintf(stderr, "!!! Input sequence is too long (>%d)!!!\n", MOTIF2_BUFFER);
		}
		buf[MOTIF2_BUFFER-1]='\0';
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 2) continue;
		GroupInfo* gi = (GroupInfo*)groupdata->search(cols[0]);
		if (gi == NULL) {
			continue;
		}
		name = cols[0];
		seq = cols[1];
		cleanSequence(seq);
		if (strand == MOTIF_STRAND_NEG) {
			Oligo::revopp(seq);
		}
			
		curNumSequences++;
		if (curNumSequences > curThreshold) {
			curThreshold += barLength;
			fprintf(stderr, "=");
		}
		double w = gi->weight;
		unsigned int g = gi->group;
		double pvalue = 0.0;
		double nvalue = 0.0;
		if (g==0) nvalue = w;
		if (g==1) pvalue = w;


		int length = strlen(seq);

		double numberOfGoodOligos = 0.0;
		for (int i=0;i<=length-(int)oligoLength;i++) {
			char tmp = seq[i+oligoLength];
			seq[i+oligoLength]='\0';
			char* oligo = &(seq[i]);
			int ok = checkSequence(oligo);
			if (ok) numberOfGoodOligos += 1.0;
			seq[i+oligoLength]=tmp;
		}
		if (numberOfGoodOligos < 1.0) {
			continue;
		}

		numTargets += pvalue;
		numBackground += nvalue;

		if (numberOfGoodOligos > avgSeqLength) {
			pvalue *= avgSeqLength/numberOfGoodOligos;
			nvalue *= avgSeqLength/numberOfGoodOligos;
		}
		for (int i=0;i<=length-(int)oligoLength;i++) {
			char tmp = seq[i+oligoLength];
			seq[i+oligoLength]='\0';
			char* oligo = &(seq[i]);
			int ok = checkSequence(oligo);
			if (ok) {
				totalOligos += oligoArray->addOligo(oligo,pvalue,nvalue,-1.0);
				//totalOligos += oligoArray->addOligo(oligo,pvalue,nvalue,(double)curNumSequences);
				totalTargetOligos += pvalue;
				totalBackgroundOligos += nvalue;
			}
			seq[i+oligoLength]=tmp;
		}
		
	}
	fprintf(stderr, "\n");
	fclose(fp);
	fprintf(stderr, "\tTotal number of Oligos: %lld\n", totalOligos);

	if (numTargets < 0.1 || numBackground < 0.1) {
		fprintf(stderr, "!!! Something is wrong - either no target sequences or no background sequences!!!\n");
		exit(0);
	}

	double oligoPerSeqTarget = totalTargetOligos/numTargets;
	//double oligoPerSeqBackround = totalBackgroundOligos/numBackground;
	double newNumBackground = totalBackgroundOligos/oligoPerSeqTarget;
	fprintf(stderr, "\tAutoadjustment for sequence coverage in background: %.2lfx\n",newNumBackground/numBackground);
	numBackground = newNumBackground;

	char** keys = groupdata->keys();
	for (int i=0;i<groupdata->total;i++) {
		GroupInfo* gi = (GroupInfo*)groupdata->search(keys[i]);
		delete gi;
		delete [](keys[i]);
	}
	delete []keys;
	delete groupdata;
	delete []buf;
	delete []cols;

	oligoArray->optimizeMemoryUsage();

}

void DeNovoMotifs::loadFullSequences() {

	if (fullSequences != NULL) delete fullSequences;
	fullSequences = new SequenceArray();
	fullSequences->readSequenceFiles(seqFile,groupFile);
	fullSequences->initializeScorer(scoringMethod,SCORER_STANDARD_ZOOPS,maxCacheSize);
	//if (normLength > 0 && !findFlag) {
	//	seqArray->normalizeSequence(normLength, strand, normFile);
	//}
}

int checkSequence(char* seq) {
	int i=0;
	while (seq[i] != '\0') {
		switch (seq[i]) {
			case 'a':
				seq[i] = 'A';
				break;
			case 'c':
				seq[i] = 'C';
				break;
			case 'g':
				seq[i] = 'G';
				break;
			case 't':
				seq[i] = 'T';
				break;
			case 'A':
			case 'C':
			case 'G':
			case 'T':
				break;
			default:
				return 0;
		}
		i++;
	}
	return 1;
}

ThreadArgs_findMisMatchEnrichmentThread::ThreadArgs_findMisMatchEnrichmentThread(
							DeNovoMotifs* obj,int cpuNum,int mis,int lasttotal,int vlevel) {
	denovo = obj;
	cpu = cpuNum;
	maxMisMatches = mis;
	totalChecked=0;
	numSkipped=0;
	dueToTargetTags=0;
    dueToLessThanLast=0;
   	dueToLowPvalue=0;
    lastChecked=0;
    lastCheckedTotal= lasttotal;
	verboseLevel = vlevel;
}


void DeNovoMotifs::findMisMatchEnrichment(int maxMisMatches,FILE* fp,int verboseLevel) {
	//oligoArray->resetOligoArray(); not needed with "activeOligos"
	int inc = (int)(((double)oligoArray->numActiveOligos)/80.0);
	int lastCheckedTotal = (int)oligoArray->numActiveOligos;
	if (verboseLevel != HOMER_VERBOSE_LEVEL0) {
		fprintf(stderr, "\n\tGlobal Optimization Phase: Looking for enriched oligos with up to %d mismatches...\n",maxMisMatches);
	}
	for (int j=0;j<=maxMisMatches;j++) {
		if (verboseLevel != HOMER_VERBOSE_LEVEL0) {
			fprintf(stderr, "\n\tScreening oligos %d (allowing %d mismatches):\n", oligoArray->numOligos, j);
			fprintf(stderr, "\t|0%%                                   50%%                                  100%%|\n");
			fprintf(stderr, "\t");
		}
		int totalChecked=0;
		int numSkipped=0;
		int dueToTargetTags=0;
		int dueToLessThanLast=0;
		int dueToLowPvalue=0;
		int lastChecked = 0;

		
		pthread_mutex_init(&mutex, NULL);
		mutexIndex = 0;
		mutexTotal = oligoArray->numActiveOligos;
		nextMeter = inc;
		incMeter = inc;
		ThreadArgs_findMisMatchEnrichmentThread** args = new ThreadArgs_findMisMatchEnrichmentThread*[numCPUs];

		pthread_t* thread = new pthread_t[numCPUs];
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		void* status;
		for (unsigned int i=0;i<numCPUs;i++) {
			args[i] = new ThreadArgs_findMisMatchEnrichmentThread(this,i,j,lastCheckedTotal,verboseLevel);
			(void)pthread_create(&(thread[i]), &attr, 
							&DeNovoMotifs_findMisMatchEnrichmentThread, args[i]);
		}
		pthread_attr_destroy(&attr);
		for (unsigned int i=0;i<numCPUs;i++) {
			(void)pthread_join(thread[i], &status);
			totalChecked += args[i]->totalChecked;
			numSkipped += args[i]->numSkipped;
			dueToTargetTags += args[i]->dueToTargetTags;
			dueToLessThanLast += args[i]->dueToLessThanLast;
			dueToLowPvalue += args[i]->dueToLowPvalue;
			if (args[i]->lastChecked > lastChecked) lastChecked = args[i]->lastChecked;
			delete args[i];
		}
		delete []thread;
		delete []args;
		pthread_mutex_destroy(&mutex);

		if (verboseLevel != HOMER_VERBOSE_LEVEL0) {
			fprintf(stderr, "\n\t\t%.2lf%% skipped, %.2lf%% checked (%d of %d), of those checked:\n", 
					(100.0*numSkipped)/oligoArray->numActiveOligos,
					(100.0*totalChecked)/oligoArray->numActiveOligos,
					totalChecked, oligoArray->numActiveOligos);
			fprintf(stderr, "\t\t%.2lf%% not in target, %.2lf%% increased p-value, %.2lf%% high p-value\n", 
					(100.0*dueToTargetTags)/oligoArray->numActiveOligos, 
					(100.0*dueToLessThanLast)/oligoArray->numActiveOligos,
					(dueToLowPvalue*100.0)/oligoArray->numActiveOligos);
		}
		inc = (int)(((double)lastChecked)/80.0);
		lastCheckedTotal = lastChecked;
	}
	oligoArray->sortActiveOligos();
	if (fp != NULL) {
		fprintf(fp, "OligoSequence\tTarget Occurences\tBackground Occurences\tEnrichment(ln p, up to %d mismatches)\n",
					maxMisMatches);
		for (int i=0;i<oligoArray->numActiveOligos;i++) {
			fprintf(fp, "%s\t%.1lf\t%.1lf\t%.3lf\n",oligoArray->activeOligos[i]->seq,oligoArray->activeOligos[i]->numTarget,
					oligoArray->activeOligos[i]->numBackground,oligoArray->activeOligos[i]->value);
		}
	}
}

void* DeNovoMotifs_findMisMatchEnrichmentThread(void* threadArgs) {

	ThreadArgs_findMisMatchEnrichmentThread* args = (ThreadArgs_findMisMatchEnrichmentThread*)threadArgs;
	args->denovo->findMisMatchEnrichmentThread(args);
	pthread_exit(threadArgs);
	return NULL;
}

void DeNovoMotifs::findMisMatchEnrichmentThread(ThreadArgs_findMisMatchEnrichmentThread* args) {

	int i=-1;

	while (i < mutexTotal) {

		pthread_mutex_lock(&mutex);
		i=mutexIndex++;
		if (i > nextMeter) {
			///////////////////Need to deal with "next"
			if (i <= args->lastCheckedTotal) {
				if (args->verboseLevel != HOMER_VERBOSE_LEVEL0) {
					fprintf(stderr, "=");
				}
			}
			//fprintf(stderr, "=%d %d",i,next);
			nextMeter += incMeter;
		}
		pthread_mutex_unlock(&mutex);
		if (i >= mutexTotal) break;

		double logp = 0.0;
		double vt = (double)oligoArray->activeOligos[i]->numTarget;
		double vbg = (double)oligoArray->activeOligos[i]->numBackground;
		int skipFlag = oligoArray->activeOligos[i]->flag2;
		if (!skipFlag && misMatchTargetOnlyFlag) {
			if (vt < 0.01) {
				skipFlag = 1;
				oligoArray->activeOligos[i]->flag2 = 1;
				args->dueToTargetTags++;
			}
		}
		if (!skipFlag) {
			if (args->maxMisMatches == 0) {
			} else {
				vt = 0.0;
				vbg = 0.0;
				oligoArray->searchOligo(oligoArray->activeOligos[i]->seq,
									args->maxMisMatches,vt,vbg,OLIGO_SEARCH_KEEP,args->cpu);
			}
	
			logp = scorer->scoreEnrichment(vt,vbg);
			if (oligoArray->activeOligos[i]->value > logp) {
				oligoArray->activeOligos[i]->value = logp;
			} else {
				if (minimumMisMatchToStartSkipping <= args->maxMisMatches 
							&& (vt+vbg > minimumCountsToStartSkipping)) {
					oligoArray->activeOligos[i]->flag2 = 1;
					args->dueToLessThanLast++;
				}
			}
			if ((minimumMisMatchToStartSkipping <= args->maxMisMatches)
						&& (oligoArray->activeOligos[i]->value > misMatchMinPvalue)
						&& (vt+vbg > minimumCountsToStartSkipping)) {
				oligoArray->activeOligos[i]->flag2 = 1;
				args->dueToLowPvalue++;
			}
			args->totalChecked++;
			args->lastChecked = i+1;
		} else {
			args->numSkipped++;
		}
	}
}


void DeNovoMotifs::getSeedOligos() {
	//assmes that oligoArray->oligos is sorted in terms of enrichment
	int maxSeedOligos = (int)MOTIF_MAX_SEED_OLIGOS_TOTAL;
	if (maxSeedOligos > oligoArray->numOligos) {
		maxSeedOligos = oligoArray->numOligos;
	}
	seedOligos = new Oligo*[maxSeedOligos];
	for (int i=0;i<maxSeedOligos;i++) {
		seedOligos[i] = oligoArray->activeOligos[i];
	}
	numSeedOligos = maxSeedOligos;
}

void DeNovoMotifs::initializeEnrichmentScoring() {
	scorer = new EnrichmentScorer();
	scorer->initializeEnrichmentScoring(scoringMethod,numTargets, numBackground,
															zoopsFlag,maxCacheSize);
	oligoArray->initializeOligoPvalues(scorer);
}

void DeNovoMotifs::optimizeMotifMatrix(Motif* &motif) {

	int cpu = 0;
	optimizeMotifThreshold(motif,cpu);
	if (motif->badFlag) return;

	if (motifOptimizationMethod == MOTIF_OPTIMZATION_DIFFERENT_THRESHOLDS) {

		Motif** newMotifs = new Motif*[numTrialMotifs];
		int newMotifIndex = 0;
		Motif* matrix = new Motif();
		double ogNumTarget = motif->numTarget;
		double ogNumBackground = motif->numBackground;
		double ogLogp = motif->logp;
		matrix->setBlank(oligoLength);

		Oligo** utilityOligos = oligoArray->utilityOligos[cpu];
		int numUtilityOligos = oligoArray->numUtilityOligos[cpu];

		int halfTrialMotifs = (int)((numTrialMotifs+1)/2.0);	
		if (halfTrialMotifs < 1) halfTrialMotifs=1;

		double inc = (utilityOligos[0]->value - motif->threshold)/(halfTrialMotifs);
		if (inc < (utilityOligos[0]->value)/((double)2*halfTrialMotifs)) {
			inc = (utilityOligos[0]->value)/((double)2*halfTrialMotifs);
		}

		double nextThreshold = utilityOligos[0]->value - inc;

		for (int i=0;i<numUtilityOligos;i++) {
			double currentThreshold = utilityOligos[i]->value;

			
			double nlogp = scorer->scoreEnrichment(ogNumTarget+utilityOligos[i]->numTarget,
													ogNumBackground+utilityOligos[i]->numBackground);
			utilityOligos[i]->value = 0;
			double diffPvalue = nlogp - ogLogp;
			//diffPvalue = -utilityOligos[i]->numTarget;
			if (diffPvalue > 0.0) {
				utilityOligos[i]->value = 0;
				//utilityOligos[i]->value = utilityOligos[i]->logp+diffPvalue;
				//if (utilityOligos[i]->value > 0) utilityOligos[i]->value = 0.0;
			} else {
				//utilityOligos[i]->value = diffPvalue;
				utilityOligos[i]->value = utilityOligos[i]->logp;
			}
				

			
			//set oligo contribution value
			matrix->addOligo(utilityOligos[i]);
			//fprintf(stderr, "======== %s %f\n", utilityOligos[i]->seq, utilityOligos[i]->logp);

			if (i+1<numUtilityOligos) {
				// need to keep going if a bunch of oligos have the same threshold
				double diff = currentThreshold - utilityOligos[i+1]->value;
				if (diff < FLOAT_EPS) {
					continue;
				}
				//fprintf(stderr, "diff=%le\n", diff);
			}

			if (currentThreshold < nextThreshold) {
				//matrix->print(stderr);
				newMotifs[newMotifIndex++] = new Motif(matrix);
				nextThreshold -= inc;
				if (newMotifIndex >= numTrialMotifs) break;
			}
		}


		//easily part to make parallel
		pthread_mutex_init(&mutex, NULL);
		mutexIndex = 0;
		mutexTotal = newMotifIndex;

		ThreadArgs_optimizeMotifThresholdThread** args = new ThreadArgs_optimizeMotifThresholdThread*[numCPUs];

		pthread_t* thread = new pthread_t[numCPUs];
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		void* status;
		for (unsigned int i=0;i<numCPUs;i++) {
			args[i] = new ThreadArgs_optimizeMotifThresholdThread(this,i,newMotifs);
			(void)pthread_create(&(thread[i]), &attr, 
							&DeNovoMotifs_optimizeMotifThresholdThread, args[i]);
		}
		pthread_attr_destroy(&attr);
		for (unsigned int i=0;i<numCPUs;i++) {
			(void)pthread_join(thread[i], &status);
			delete args[i];
		}
		delete []thread;
		delete []args;
		pthread_mutex_destroy(&mutex);


		int bestIndex = -1;
		double bestIndexScore  = motif->logp; 
		for (int i=0;i<newMotifIndex;i++) {
			//newMotifs[i]->print(stderr);
			if (newMotifs[i]->logp < bestIndexScore) {
				bestIndex = i;
				bestIndexScore = newMotifs[i]->logp;
			}
		}
		//fprintf(stderr, "Best Index = %d of %d (new: %lf\told: %lf)\n",bestIndex, newMotifIndex,bestIndexScore,motif->logp);
		if (bestIndex > -1) {
			Motif* tmp = motif;
			motif = newMotifs[bestIndex];
			newMotifs[bestIndex] = tmp;
		}
		for (int i=0;i<newMotifIndex;i++) {
			delete newMotifs[i];
		}
		delete matrix;	
		delete []newMotifs;
	}	
}

void* DeNovoMotifs_optimizeMotifThresholdThread(void* threadArgs) {
	ThreadArgs_optimizeMotifThresholdThread* args = (ThreadArgs_optimizeMotifThresholdThread*)threadArgs;
	args->denovo->optimizeMotifThresholdThread(args);
	pthread_exit(threadArgs);
	return NULL;
}
ThreadArgs_optimizeMotifThresholdThread::ThreadArgs_optimizeMotifThresholdThread(
								DeNovoMotifs* obj,int newcpu,Motif** newmotifs) {
	denovo = obj;
	cpu = newcpu;
	motifs=newmotifs;
}
void DeNovoMotifs::optimizeMotifThresholdThread(ThreadArgs_optimizeMotifThresholdThread* args) {
	int i=-1;
	while (i<mutexTotal) {
		pthread_mutex_lock(&mutex);
		i=mutexIndex++;
		pthread_mutex_unlock(&mutex);
		if (i >= mutexTotal) break;

		args->motifs[i]->normalize();
		optimizeMotifThreshold(args->motifs[i],args->cpu);
	}
}

void DeNovoMotifs::optimizeMotifThreshold(Motif* motif, int cpu) {

	motif->setMinimumScore(maximumExpectedPerBp);
	oligoArray->scoreOligosWithMotif(motif,motif->minimumScore,MOTIF_SCORE_OLIGOS,cpu);
	oligoArray->sortUtilityOligos(cpu);

	double bestThreshold = 0;
	double bestScore = 0;
	double bestTvalue = 0;
	double bestBGvalue = 0;
	double currentTvalue = 0;
	double currentBGvalue = 0;
	double currentThreshold = 0;
	double maxBackground = maxBackgroundPercent*numBackground;
	Oligo** utilityOligos = oligoArray->utilityOligos[cpu];
	int numUtilityOligos = oligoArray->numUtilityOligos[cpu];
	if (numUtilityOligos < 1) {
		// basically, this is no longer a valid motif... need to bail gracefully...
		motif->badFlag = 1;
		return;
	}
	//int bestIndex = 0;
	//int lastChecked = 0;
	double bestPossibleThreshold = utilityOligos[0]->cpuData[cpu];
	double lowestPossibleThreshold = utilityOligos[numUtilityOligos-1]->cpuData[cpu];
	//fprintf(stderr, "\t%lf\t%lf\t%lf\n", motif->minimumScore, bestPossibleThreshold, lowestPossibleThreshold);
	for (int i=0;i<numUtilityOligos;i++) {
		currentTvalue += utilityOligos[i]->numTarget;
		currentBGvalue += utilityOligos[i]->numBackground;
		currentThreshold = utilityOligos[i]->cpuData[cpu];
		if (i+1<numUtilityOligos) {
			// need to keep going if a bunch of oligos have the same threshold
			double diff = currentThreshold - utilityOligos[i+1]->cpuData[cpu];
			if (diff < FLOAT_EPS) {
				continue;
			}
			//fprintf(stderr, "diff=%le\n", diff);
		}
		double logp = scorer->scoreEnrichment(currentTvalue,currentBGvalue);
		//fprintf(stderr, "%lf\t%lf\t%lf\t%lf\n", currentThreshold, logp,currentTvalue,numTargets);
		if (logp < bestScore) {
			//bestIndex = i;
			bestScore = logp;
			bestThreshold = currentThreshold;
			bestTvalue = currentTvalue;
			bestBGvalue = currentBGvalue;
		}
		if (currentThreshold < motif->minimumScore) break;
		if (logp > bestScore/3.0 && currentTvalue > numTargets) break;
		if (currentBGvalue > maxBackground) break;
		//lastChecked = i;
		/*fprintf(stderr, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\n", utilityOligos[i]->seq, 
					currentTvalue, currentBGvalue, currentThreshold, logp, bestScore);*/
	}
	//fprintf(stderr, "%d of %d of %d\n", bestIndex, lastChecked,numUtilityOligos);
	double degfrac = (bestPossibleThreshold-bestThreshold)
							/(bestPossibleThreshold-lowestPossibleThreshold+0.01);
	if (degfrac > 0.9) {
		fprintf(stderr, "Optimal threshold is low - motif appears very degenerate!!!\n");
	}
	
	motif->threshold = bestThreshold-0.000001;
	motif->logp = bestScore;
	motif->numTarget = bestTvalue;
	motif->numBackground = bestBGvalue;
	motif->freqTarget = bestTvalue/numTargets;
	motif->freqBackground = bestBGvalue/numBackground;
	if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
		motif->freqTarget = scorer->getZoopsApproxTarget(bestTvalue)/numTargets;
		motif->freqBackground = scorer->getZoopsApproxBackground(bestBGvalue)/numBackground;
	}
//	fprintf(stderr, "\t%.3lf\t%.3lf\n", motif->minimumScore, motif->threshold);
}

void DeNovoMotifs::optimizeGivenMotifs() {
	int defaultCPU = 0;

	FILE* tmpfp = NULL;
	if (tmpFile == NULL) {
		tmpfp = fopen(".tmp.motifs","w");
		if (tmpfp == NULL) fprintf(stderr, "\t!! Could not open tmp results file (.tmp.motifs) !!\n");
	} else {
		tmpfp = fopen(tmpFile,"w");
		if (tmpfp == NULL) fprintf(stderr, "\t!! Could not open tmp results file (%s) !!\n", tmpFile);
	}

	fprintf(stderr, "\n\tLocal Optimization Phase:\n");
	for (int i=0;i<numMotifs;i++) {

		motifs[i]->getConsensus();
		fprintf(stderr, "\t\tInitial Consensus: %s\n", motifs[i]->consensus);
	
		double curLogp = 0.0;
		for (int j=0;j<maxOptimizationIterations;j++) {
			optimizeMotifMatrix(motifs[i]);
			motifs[i]->getConsensus();
			fprintf(stderr, "\t\tRound %d: %.2lf", j+1, motifs[i]->logp);
			fprintf(stderr, " %s ", motifs[i]->consensus);
			fprintf(stderr, "T:%.1lf(%.2f%%),B:%.1lf(%.2f%%),P:1e%d\n",motifs[i]->numTarget, 
									motifs[i]->freqTarget*100, motifs[i]->numBackground,
									motifs[i]->freqBackground*100,(int)(motifs[i]->logp/2.3026));
			if (curLogp-FLOAT_EPS < motifs[i]->logp) {
				break;
			}
			curLogp = motifs[i]->logp;
		}

		motifs[i]->getConsensus();
		fprintf(stderr, "\t\t=Final=: %.2lf", motifs[i]->logp);
		fprintf(stderr, " %s ", motifs[i]->consensus);
		fprintf(stderr, "T:%.1lf(%.2f%%),B:%.1lf(%.2f%%),P:1e%d\n",motifs[i]->numTarget, 
									motifs[i]->freqTarget*100, motifs[i]->numBackground,
									motifs[i]->freqBackground*100,(int)(motifs[i]->logp/2.3026));

		//do not remove oligos from the tree for this type of motif finding
		oligoArray->scoreOligosWithMotif(motifs[i],motifs[i]->minimumScore,MOTIF_SCORE_OLIGOS,defaultCPU);

		int totalBoundOligos=0;
		for (int j=0;j<oligoArray->numUtilityOligos[0];j++) {
			if (oligoArray->utilityOligos[0][j]->value > motifs[i]->threshold) {
				totalBoundOligos++;
			}
		}
		motifs[i]->getConsensus();
		fprintf(stderr, "\t\t%d of %d - %s Bound Oligos: %d\n",i+1,totalMotifsToFind,
											motifs[i]->consensus,totalBoundOligos);
		motifs[i]->print(tmpfp);
		if (tmpfp != NULL) fflush(tmpfp);
	}
	if (tmpfp != NULL) fclose(tmpfp);
}

void DeNovoMotifs::optimizeSeeds2Motifs() {

	//assmes that oligoArray->activeOligos is sorted in terms of enrichment
	if (oligoArray->numActiveOligos > oligoArray->numOligos*MOTIF_MAX_SEED_OLIGOS_FRACTION) {
		oligoArray->numActiveOligos = (int)(oligoArray->numOligos*MOTIF_MAX_SEED_OLIGOS_FRACTION);
	}
	if (oligoArray->numActiveOligos > MOTIF_MAX_SEED_OLIGOS_TOTAL) {
		oligoArray->numActiveOligos = (int)(MOTIF_MAX_SEED_OLIGOS_TOTAL);
	}
	if (minimumSeedLogp > MOTIF_SEED_MINIMUM_LOGP-0.01 
			&& minimumSeedLogp < MOTIF_SEED_MINIMUM_LOGP+0.01) {
		//default value - adjust for the target genes, the less of which lowers pvalues
		minimumSeedLogp = MOTIF_SEED_MINIMUM_LOGP+(4.0-log((double)numTargets)/log(10.0));
	}

	motifs = new Motif*[totalMotifsToFind];
	numMotifs = 0;
	int seedOligoIndex = 0;
	char* bufferSeq = new char[oligoLength+1];

	OligoArray* boundOligos = new OligoArray();
	boundOligos->init(oligoLength,numOligoEstimate/4,alphaSize,strand,1);
	int defaultCPU = 0;

	FILE* tmpfp = NULL;
	if (tmpFile == NULL) {
		tmpfp = fopen(".tmp.motifs","w");
		if (tmpfp == NULL) fprintf(stderr, "\t!! Could not open tmp results file (.tmp.motifs) !!\n");
	} else {
		tmpfp = fopen(tmpFile,"w");
		if (tmpfp == NULL) fprintf(stderr, "\t!! Could not open tmp results file (%s) !!\n", tmpFile);
	}


	fprintf(stderr, "\n\tLocal Optimization Phase:\n");
	for (int i=0;i<totalMotifsToFind;i++) {

		if (seedOligoIndex >= oligoArray->numActiveOligos) break;
		char* seedSeq = oligoArray->activeOligos[seedOligoIndex]->seq;
		double seedlogp = oligoArray->activeOligos[seedOligoIndex]->value;

		if (i > 4 && seedlogp > minimumSeedLogp && maskMethod == DENOVO_MASK_FULL) {
			fprintf(stderr, "\tRemaining seeds don't look promising (After initial 5 motifs, logp %.3lf > %.3lf)\n",
											seedlogp,minimumSeedLogp); 
			//fprintf(stderr, "\t\tTo analyze further, use -minlp <#> option\n"); 
			//fprintf(stderr, "\t\t%s\n", oligoArray->activeOligos[seedOligoIndex]->seq); 
			break;
		}
		oligoArray->activeOligos[seedOligoIndex++] = NULL;
		//fprintf(stderr, "\t%s",seedSeq);

		//first check if the oligo still exists and hasn't been masked out yet
		int there=0;
		for (int j=0;j<=bufferLength;j++) {
			double vt=0.0,vbg=0.0;
			for (int k=0;k<j;k++) bufferSeq[k]='N';
			for (int k=0;k<(int)oligoLength-j;k++) bufferSeq[k+j]=seedSeq[k];
			bufferSeq[oligoLength]='\0';
			there = boundOligos->searchOligo(bufferSeq,0,vt,vbg,OLIGO_SEARCH_KEEP,defaultCPU);
			if (there > 0) break;
			if (j==0) continue;
			vt=0.0,vbg=0.0;
			for (int k=0;k<j;k++) bufferSeq[oligoLength-k-1]='N';
			for (int k=0;k<(int)oligoLength-j;k++) bufferSeq[k]=seedSeq[k+j];
			bufferSeq[oligoLength]='\0';
			there = boundOligos->searchOligo(bufferSeq,0,vt,vbg,OLIGO_SEARCH_KEEP,defaultCPU);
			if (there > 0) break;
		}
		if (there > 0) { 
			i--;
			//fprintf(stderr, "\t%s... too similar, skipping\n",seedSeq);
			continue;
		}
		fprintf(stderr, "\t%d of %d Initial Sequence: %s... (%.3lf)\n",i+1,totalMotifsToFind,
															seedSeq,seedlogp);
		//fprintf(stderr, "...\n");


		motifs[i] = new Motif(seedSeq);
		motifs[i]->getConsensus();
		//fprintf(stderr, "\t\tInitial Consensus: %s\n", motifs[i]->consensus);
		numMotifs=i+1;
	
		double curLogp = 0.0;
		for (int j=0;j<maxOptimizationIterations;j++) {
			optimizeMotifMatrix(motifs[i]);
			if (motifs[i]->badFlag) {
				fprintf(stderr, "\t\t!! Motif no longer valid - masked out of existance...\n");
				break;
			}
			motifs[i]->getConsensus();
			fprintf(stderr, "\t\tRound %d: %.2lf", j+1, motifs[i]->logp);
			fprintf(stderr, " %s ", motifs[i]->consensus);
			fprintf(stderr, "T:%.1lf(%.2f%%),B:%.1lf(%.2f%%),P:1e%d\n",motifs[i]->numTarget, 
									motifs[i]->freqTarget*100, motifs[i]->numBackground,
									motifs[i]->freqBackground*100,(int)(motifs[i]->logp/2.3026));
			if (curLogp-FLOAT_EPS < motifs[i]->logp) {
				break;
			}
			curLogp = motifs[i]->logp;
		}
		if (motifs[i]->badFlag) {
			continue;
		}
		if (maskMethod == DENOVO_MASK_FULL) {
			//first optimize motif based on actual zoops scoring
			fullSequences->optimizeMotifThreshold(motifs[i],strand);
		}

		motifs[i]->getConsensus();
		fprintf(stderr, "\t\t=Final=: %.2lf", motifs[i]->logp);
		fprintf(stderr, " %s ", motifs[i]->consensus);
		fprintf(stderr, "T:%.1lf(%.2f%%),B:%.1lf(%.2f%%),P:1e%d\n",motifs[i]->numTarget, 
						motifs[i]->freqTarget*100, motifs[i]->numBackground,
						motifs[i]->freqBackground*100,(int)(motifs[i]->logp/2.3026));
		motifs[i]->print(tmpfp);
		if (tmpfp != NULL) fflush(tmpfp);


		if (maskMethod == DENOVO_MASK_FULL) {
			//remove sequences adjacent to motif from the original sequence that overlap
			fprintf(stderr, "\t\tPerforming exhaustive masking of motif...\n");
			fullSequences->maskMotifFromOligoArray(motifs[i],oligoArray,strand);

			//resort active Oligos to 
			fprintf(stderr, "\t\tReprioritizing potential motifs...\n");
			for (int i=0;i<oligoArray->numActiveOligos;i++) {
				if (oligoArray->activeOligos[i] == NULL) {
					oligoArray->activeOligos[i] = oligoArray->activeOligos[oligoArray->numActiveOligos-1];
					oligoArray->numActiveOligos--;
				} else {
					break;
				}
			}
			oligoArray->initializeOligoPvalues(scorer);
			findMisMatchEnrichment(maxMisMatches,NULL,HOMER_VERBOSE_LEVEL0);	
			seedOligoIndex=0;

		} else {

			//remove oligos from the tree to avoid finding the same motif again
			//remove bound oligos from data set
			oligoArray->scoreOligosWithMotif(motifs[i],motifs[i]->minimumScore,MOTIF_REMOVE_OLIGOS,defaultCPU);
	
			int totalBoundOligos=0;
			for (int j=0;j<oligoArray->numUtilityOligos[0];j++) {
				if (oligoArray->utilityOligos[0][j]->value > motifs[i]->threshold) {
					boundOligos->addOligo(oligoArray->utilityOligos[0][j]->seq,1.0,1.0,-1.0);
					totalBoundOligos++;
				}
			}
			fprintf(stderr, "\t\t%d of %d - %s Bound Oligos: %d\n",i+1,totalMotifsToFind,
											motifs[i]->consensus,totalBoundOligos);

		}


	}
	if (tmpfp != NULL) fclose(tmpfp);

	delete []bufferSeq;
	delete boundOligos;
}

void DeNovoMotifs::printMotifs(FILE* fp) {
	if (numMotifs < 1) return;

	char* name = new char[10000];	
	qsort(motifs, numMotifs, sizeof(Motif*), &(Motif::motifCmp));
	for (int i=0;i<numMotifs;i++) {
		motifs[i]->getConsensus();
		sprintf(name,"%d-%s", i+1, motifs[i]->consensus);
		motifs[i]->setName(name);
		motifs[i]->print(fp);
	}
	delete []name;
}
void DeNovoMotifs::scoreWithKnownMotifs(int optFlag) {

	KnownMotifs* known = new KnownMotifs();
	known->seqFile = seqFile;
	known->groupFile = groupFile;
	known->strand = strand;
	known->numCPUs = numCPUs;
	known->scoringMethod = scoringMethod;
	known->motifs = motifs;
	known->numMotifs = numMotifs;
	known->threasholdOptimizationFlag = optFlag;
	//known->siteReduceThreshold = siteReduceThreshold;

	known->loadSequence();
	known->getMotifEnrichment();
	if (known->siteReduceThreshold > 0) known->findOverlappingMotifs();

	motifs = known->motifs;
	numMotifs = known->numMotifs;
}

//============================ OligoArray ========================================

OligoArray::OligoArray() {
	tree = NULL;
	numTree = 0;
	maxTree = 0;
	alphaSize = DNA_ALPHA;
	oligoLength = 0;
	alpha = NULL;
	rvAlpha = NULL;
	rvMap = NULL;
	numCPUs = 1;

	pthread_mutex_init(&sortMutex, NULL);
	
	staticOligos = NULL;
	utilityOligos = NULL;
	numUtilityOligos = 0;
	numOligos = 0;
	activeOligos = NULL;
	numActiveOligos = 0;
	maxOligos = 0;

	rvFlag = 1;	

	indexStack = NULL;
	alphaStack = NULL;
	positionStack = NULL;
	fwdMisMatchStack = NULL;
	rvMisMatchStack = NULL;

	totalTarget = 0.0;	
	totalBackground = 0.0;	
	
	maxAutoNormIters = MAX_OLIGO_NORMALIZATION_ITERATIONS;
}

OligoArray::~OligoArray() {
	if (tree != NULL) delete []tree;
	if (staticOligos != NULL) {
		for (int i=0;i<numOligos;i++) {
			delete staticOligos[i];
		}
		delete []staticOligos;
	}
	if (alpha != NULL) delete []alpha;
	if (rvAlpha != NULL) delete []rvAlpha;
	if (rvMap != NULL) delete []rvMap;
	if (alphaIndex != NULL) delete []alphaIndex;
	if (rvAlphaIndex != NULL) delete []rvAlphaIndex;

	pthread_mutex_destroy(&sortMutex);

	if (numUtilityOligos != NULL) delete []numUtilityOligos;

	for (int i=0;i<numCPUs;i++) {	
		if (utilityOligos != NULL && utilityOligos[i] != NULL) delete [](utilityOligos[i]);
		if (indexStack != NULL && indexStack[i] != NULL) delete [](indexStack[i]);
		if (alphaStack != NULL && alphaStack[i] != NULL) delete [](alphaStack[i]);
		if (positionStack != NULL && positionStack[i] != NULL) delete [](positionStack[i]);
		if (fwdMisMatchStack != NULL && fwdMisMatchStack[i] != NULL) delete [](fwdMisMatchStack[i]);
		if (rvMisMatchStack != NULL && rvMisMatchStack[i] != NULL) delete [](rvMisMatchStack[i]);
		if (fwdScoreStack != NULL && fwdScoreStack[i] != NULL) delete [](fwdScoreStack[i]);
		if (rvScoreStack != NULL && rvScoreStack[i] != NULL) delete [](rvScoreStack[i]);
	}
	if (utilityOligos != NULL) delete []utilityOligos;
	if (indexStack != NULL) delete []indexStack;
	if (alphaStack != NULL) delete []alphaStack;
	if (positionStack != NULL) delete []positionStack;
	if (fwdMisMatchStack != NULL) delete []fwdMisMatchStack;
	if (rvMisMatchStack != NULL) delete []rvMisMatchStack;
	if (fwdScoreStack != NULL) delete []fwdScoreStack;
	if (rvScoreStack != NULL) delete []rvScoreStack;
}

void OligoArray::init(int length, int totalDataSetBp, int alphaSIZE, int strandFlag, int numcpu) {

	alphaSize = alphaSIZE;
	oligoLength = length;
	if (alphaSize == DNA_ALPHA) {
		alpha=new char[4];
		rvAlpha=new char[4];
		alpha[0] = 'A';
		alpha[1] = 'C';
		alpha[2] = 'G';
		alpha[3] = 'T';
		rvAlpha[0] = 'T';
		rvAlpha[1] = 'G';
		rvAlpha[2] = 'C';
		rvAlpha[3] = 'A';
		rvMap = new char[256];
		alphaIndex = new char[256];
		rvAlphaIndex = new char[256];
		rvMap[65]='T'; //rvMap['A']='T';
		rvMap[67]='G'; //rvMap['C']='G';
		rvMap[71]='C'; //rvMap['G']='C';
		rvMap[84]='A'; //rvMap['T']='A';
		alphaIndex[65]=0; //alphaIndex['A']=0;
		alphaIndex[67]=1; //alphaIndex['C']=1;
		alphaIndex[71]=2; //alphaIndex['G']=2;
		alphaIndex[84]=3; //alphaIndex['T']=3;
	}
	Motif::alphaIndex = alphaIndex;
	Motif::alphaSize = alphaSize;

	if (strandFlag == MOTIF_STRAND_BOTH) {
		rvFlag =1;
	} else {
		rvFlag = 0;
	}	
	
	maxOligos = INT_MAX;
	if (length < 20) {
		double dmax = pow(alphaSize,length)+10;
		if (rvFlag) {
			dmax = dmax/1.9;
		}
		if (dmax < (double)maxOligos) {
			maxOligos = (int) dmax;
		}
	}
	if (maxOligos > totalDataSetBp) {
		maxOligos = totalDataSetBp;
	}
	//int fullTreeEstimate = (int) (log(maxOligos)/log(4.0)+0.49);

	//maxTree = maxOligos*(3+oligoLength-fullTreeEstimate);
	maxTree = maxOligos*(5);
	numTree = 0;
	tree = new int[maxTree];
	for (int i=0;i<maxTree;i++) tree[i] = 0;

	staticOligos = new Oligo*[maxOligos+1];
	staticOligos[0] = NULL;
	for (int i=0;i<maxOligos;i++) {
		staticOligos[i] = NULL;
	}
	staticOligos[maxOligos] = NULL;
	numOligos = 0;

	numCPUs = numcpu;
	if (numCPUs < 1) numCPUs = 1;

	utilityOligos = new Oligo**[numCPUs];
	numUtilityOligos = new int[numCPUs];

	indexStack = new int*[numCPUs];
	alphaStack = new int*[numCPUs];
	positionStack = new int*[numCPUs];
	fwdMisMatchStack = new int*[numCPUs];
	rvMisMatchStack = new int*[numCPUs];
	fwdScoreStack = new double*[numCPUs];
	rvScoreStack = new double*[numCPUs];
	for (int i=0;i<numCPUs;i++) {
		if (i==0) {
			utilityOligos[i] = new Oligo*[maxOligos];
			for (int j=0;j<maxOligos;j++) utilityOligos[i][j] = NULL;
		} else {
			utilityOligos[i] = NULL;
		}

		numUtilityOligos[i]=0;
		indexStack[i] = new int[oligoLength];
		alphaStack[i] = new int[oligoLength];
		positionStack[i] = new int[oligoLength];
		fwdMisMatchStack[i] = new int[oligoLength];
		rvMisMatchStack[i] = new int[oligoLength];
		fwdScoreStack[i] = new double[oligoLength];
		rvScoreStack[i] = new double[oligoLength];
	}

	//fprintf(stderr, "\tEstimate of prefix tree depth = %d\n", fullTreeEstimate);
	//fprintf(stderr, "\tMax # of Oligos               = %d\n", maxOligos);
	//fprintf(stderr, "\tMax size of Tree              = %d\n", maxTree);

}
void OligoArray::optimizeMemoryUsage() {
	fprintf(stderr, "\n\tOligos: %d of %d max\n", numOligos, maxOligos);
	fprintf(stderr, "\tTree  : %d of %d max\n", numTree, maxTree);
	fprintf(stderr, "\tOptimizing memory usage...\n");

	maxTree = numTree+4;
	int* ntree = new int[maxTree];
	for (int i=0;i<maxTree;i++) {
		ntree[i] = tree[i];
	}
	delete []tree;
	tree = ntree;

	maxOligos = numOligos;
	Oligo** tmpOligos = new Oligo*[maxOligos];
	for (int i=0;i<maxOligos;i++) {
		tmpOligos[i] = staticOligos[i+1];
	}
	delete []staticOligos;
	staticOligos = new Oligo*[maxOligos+1];
	staticOligos[0]=0;
	for (int i=0;i<maxOligos;i++) {
		staticOligos[i+1] = tmpOligos[i];
		staticOligos[i+1]->cpuData = new double[numCPUs];;
	}

	for (int i=0;i<numCPUs;i++) {
		if (utilityOligos[i] != NULL) {
			delete [](utilityOligos[i]);
			utilityOligos[i] = NULL;
		}
		utilityOligos[i] = new Oligo*[maxOligos];
		for (int j=0;j<maxOligos;j++) {
			utilityOligos[i][j] = tmpOligos[j];
		}
	}
	delete []tmpOligos;
}


void OligoArray::readOligoCountFile(char* file,int strand, int ncpus) {
	fprintf(stderr, "\n\tReading oligos from file %s\n", file);
	FILE* fp = fopen(file,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Couldn't open oligo read count file: %s !!!\n", file);
		return;
	}
	char* buf = new char[MOTIF2_BUFFER];
	char** cols = new char*[1000];
	int numCols;

	int firstOligoSize = -1;
	int numLines = 0;
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 3) continue;
		int ok = checkSequence(cols[0]);
		if (!ok) continue;
		numLines++;
		int curLength = strlen(cols[0]);
		if (firstOligoSize == -1) {	
			fprintf(stderr, "\t\tDetected Oligos of length %d from file %s\n", curLength,file);
			firstOligoSize = curLength;
		} else if (curLength != firstOligoSize) {
			fprintf(stderr, "\t\t!!! Error - different sized oligos detected (first one %d != %d) (%s)!!!\n",
									firstOligoSize,curLength,cols[0]);
		}
	}
	init(firstOligoSize,firstOligoSize*numLines,DNA_ALPHA,strand,ncpus);
	rewind(fp);
	double TV=0.0;
	double BV=0.0;
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 3) continue;
		int ok = checkSequence(cols[0]);
		if (!ok) continue;
		sscanf(cols[1],"%lf",&TV);
		sscanf(cols[2],"%lf",&BV);
		addOligo(cols[0],TV,BV,0.0);
	}
	fclose(fp);

	optimizeMemoryUsage();
	initializeActiveOligoArray();
	calculateFoldEnrichment(0.5);
}



void OligoArray::printOligo(char* seq, FILE* fp) {
	double vt=0.0,vbg=0.0;
	int cpu = 0;
	int found = searchOligo(seq,0,vt,vbg,0,cpu);
	if (found > 0) {
		utilityOligos[cpu][0]->print(fp);
		fprintf(stderr, "\t0\t%lf\t%lf\n", vt,vbg);
		searchOligo(seq,1,vt,vbg,0,cpu);
		fprintf(stderr, "\t1\t%lf\t%lf\n", vt,vbg);
		searchOligo(seq,2,vt,vbg,0,cpu);
		fprintf(stderr, "\t2\t%lf\t%lf\n", vt,vbg);
	} else {
		fprintf(stderr, "Can't find %s\n", seq);
	}
}

void OligoArray::sortActiveOligos() {
	if (activeOligos == NULL || numActiveOligos < 1) return;
	qsort(activeOligos, numActiveOligos, sizeof(Oligo*), &(Oligo::oligoCmp));
}
void OligoArray::sortUtilityOligos(int cpu) {
	pthread_mutex_lock(&sortMutex);
	Oligo::cpu = cpu;
	qsort(utilityOligos[cpu], numUtilityOligos[cpu], sizeof(Oligo*), &(Oligo::oligoCmpReverseParallel));
	pthread_mutex_unlock(&sortMutex);
}
void OligoArray::initializeActiveOligoArray() {
	if (activeOligos != NULL) delete []activeOligos;
	activeOligos = new Oligo*[numOligos];
	numActiveOligos = numOligos;
	for (int i=0;i<numOligos;i++) {
		activeOligos[i] = staticOligos[i+1];
		activeOligos[i]->flag=0;
		activeOligos[i]->flag2=0;
	}
}

int OligoArray::addOligo(char* seq, double vt, double vbg, double addValue) {
	totalTarget += vt;
	totalBackground += vbg;
	double tmpT=0.0,tmpBg=0.0;
	int cpu = 0;
	int found = searchOligo(seq,0,tmpT,tmpBg,OLIGO_SEARCH_KEEP,cpu);
	int rvMode = 0;
	if (found == 0) {
		rvMode = 0;
	} else {
		//if (found > 1) fprintf(stderr, "!!!!Problem %s\t%s\n", seq, utilityOligos[cpu][0]->seq);
		if (strcmp(seq,utilityOligos[cpu][0]->seq)==0) {
			rvMode=0;
		} else {
			rvMode=1;
		}
		if (rvMode == 1 && rvFlag == 0) {
			//fprintf(stderr, "!!!!WTF %s\t%s\n", seq, utilityOligos[0]->seq);
		}
	}
	int added = addOligoTree(seq,vt,vbg,addValue,rvMode,0);
	rvMode = 0;
	return added;
}

int OligoArray::addOligoTree(char* seq, double vt, double vbg, double addValue, int rvMode, int knownIndex) {

	int index = 0;
	int lastIndex = 0;
	int same = 0;
	int knownIndexToAdd=0;

	for (int i=0;i<oligoLength;i++) {
		char bp = '\0';
		if (rvMode) {
			bp = seq[oligoLength-i-1];
		} else {
			bp = seq[i];
		}

		if (index == numTree) {
			if (numTree+alphaSize > maxTree) {
				fprintf(stderr, "Out of room in oligo TREE! %d of %d\n", numTree+alphaSize, maxTree);
			}
			for (int j=0;j<alphaSize;j++) {
				tree[numTree]=0;
				numTree++;
			}
		} else if (index > numTree) {
			fprintf(stderr, "ERROR IN TREE GENERATION!!\n");
		}

		int indexDelta = alphaIndex[(int)bp];
		if (rvMode) indexDelta = 3-indexDelta;

		lastIndex = index+indexDelta;
		index = tree[lastIndex];

		if (index < 0 && i< oligoLength-1) {
			//this is a pointer to an oligo
			//check if equal
			same=1;
			char* otherSeq = staticOligos[-index]->seq;
			//fprintf(stderr, "\tChecking if %s = %s rvMode=%d starting at %d\n", seq, otherSeq,rvMode,i);
			for (int j=i+1;j<oligoLength;j++) {
				if (rvMode) {
					if (rvMap[(int)seq[oligoLength-j-1]] == otherSeq[j]) same++;
					else break;
				} else {
					if (seq[j] == otherSeq[j]) same++;
					else break;
				}
			}
			if (i+same < oligoLength) {
				//	fprintf(stderr, " !!!!! NOT THE SAME\n");
				//this means they aren't equal 
				knownIndexToAdd=-index;
				index = 0; //set index to zero so that a new tree node is created
			} else {
				//	fprintf(stderr, " ARE THE SAME\n");
				//sequences are the same - break and add the counts below
				//fprintf(stderr, "Sequences are the same!\n");
				break;
			}
		}
		if (index == 0) {
			if (i<oligoLength-1 && same > 0) {
				tree[lastIndex] = numTree;
				tree[numTree] = 0;
				index = numTree;
				same--;
			} else {
				//sequences haven't been added yet-
				break;
			}
		}
	}

	int oligoIndex = -tree[lastIndex];
	int added = 0;
	if (knownIndex != 0) {
		oligoIndex = knownIndex;
		tree[lastIndex]= -knownIndex;
	} else {
		if (oligoIndex == 0) {
			oligoIndex = ++numOligos;
			if (oligoIndex > maxOligos-1) {
				fprintf(stderr, "Ran out of room for more oligos %d of %d!!!!\n",oligoIndex,maxOligos);
			}
			tree[lastIndex] = -oligoIndex;
			if (rvMode) {
				fprintf(stderr, "THIS SHOULDN'T Happen!!!! (%s)\n",seq);
			}
			staticOligos[oligoIndex] = new Oligo(seq);
			//oligos[oligoIndex-1] = staticOligos[oligoIndex]; not needed
			added=1;
		} else if (oligoIndex < 0) {
			fprintf(stderr, "ERROR WITH OLIGO CREATION\n");
		} else {
		}
		if (addValue <= FLOAT_EPS || staticOligos[oligoIndex]->value + FLOAT_EPS < addValue) {
			staticOligos[oligoIndex]->numTarget += vt;
			staticOligos[oligoIndex]->numBackground += vbg;
		}
		staticOligos[oligoIndex]->value = addValue;
		
	}

	if (knownIndexToAdd != 0) {
		addOligoTree(staticOligos[knownIndexToAdd]->seq,staticOligos[knownIndexToAdd]->numTarget,
						staticOligos[knownIndexToAdd]->numBackground,addValue,rvMode,knownIndexToAdd);
	}
		

	return added;
}

void OligoArray::adjustOligoInstances(char* seq, double vt, double vbg) {
	double tmp1, tmp2;
	int cpu = 0;
	int maxMisMatches = 0;
	int found = searchOligo(seq, maxMisMatches, tmp1, tmp2, OLIGO_SEARCH_KEEP,cpu);
	if (found < 1) return;
	if (found > 1) {
		fprintf(stderr, "!!! Somthing is probably wrong... duplicates in the oligoarray\n");
	}
	double normWeight = utilityOligos[cpu][0]->normWeight;
	utilityOligos[cpu][0]->numTarget += vt;
	utilityOligos[cpu][0]->numBackground += vbg*normWeight;
	int remove = 1;
	if (utilityOligos[cpu][0]->numTarget < FLOAT_EPS) {
		utilityOligos[cpu][0]->numTarget = 0.0;
	} else {
		remove = 0;
	}
	if (utilityOligos[cpu][0]->numBackground < FLOAT_EPS) {
		utilityOligos[cpu][0]->numBackground = 0.0;
	} else {
		remove = 0;
	}
	if (remove) {
		//all instances of the oligo have been removed, so get it out of the tree
		found = searchOligo(seq, maxMisMatches, tmp1, tmp2, OLIGO_SEARCH_REMOVE,cpu);
	}
}

int OligoArray::searchOligo(char* seq,int maxMisMatches, double& vt, double& vbg, 
															int removeFlag,int cpu) {
	
	numUtilityOligos[cpu] = 0;

	indexStack[cpu][0] = 0;
	alphaStack[cpu][0] = 0;
	fwdMisMatchStack[cpu][0] = 0;
	rvMisMatchStack[cpu][0] = 0;
	int stackPosition = 0;
	int valIndex = 0;
	int oligoIndex = 0;
	
	int alphaPos = 0;
	int fwdMisMatch = 0;
	int rvMisMatch = 0;

	while (stackPosition > -1) {


		alphaPos = alphaStack[cpu][stackPosition];

		if (alphaPos >= alphaSize) {
			if (removeFlag == OLIGO_SEARCH_REMOVE) {
				int numEmpty = 0;
				for (int i=0;i<alphaSize;i++) {
					if (tree[indexStack[cpu][stackPosition]+i] == 0) numEmpty++;
				}
				if (numEmpty == alphaSize && stackPosition > 0) {
					tree[indexStack[cpu][stackPosition-1]+alphaStack[cpu][stackPosition-1]]=0;
				}
			}
			stackPosition--;
			continue;
		}


		fwdMisMatch = fwdMisMatchStack[cpu][stackPosition];
		rvMisMatch = rvMisMatchStack[cpu][stackPosition];

		if (rvFlag) {
			if ('N' != seq[oligoLength-stackPosition-1]) {
				if (rvAlpha[alphaPos] != seq[oligoLength-stackPosition-1])
					rvMisMatch++;
			}
		}
		if ('N' != seq[stackPosition]) {
			if (alpha[alphaPos] != seq[stackPosition]) fwdMisMatch++;
		}

		if (fwdMisMatch > maxMisMatches) {
			if (!rvFlag || rvMisMatch > maxMisMatches) {
				alphaStack[cpu][stackPosition]++;
				continue;
			}
		}

		valIndex = tree[indexStack[cpu][stackPosition]+alphaPos];
		if (valIndex < 0) {
			oligoIndex = -valIndex;
			//this is a pointer to an oligo
			//check if equal
			char* otherSeq = staticOligos[oligoIndex]->seq;
	
			for (int i=stackPosition+1;i<oligoLength;i++) {
				if (seq[i] != 'N' && seq[i] != otherSeq[i]) {
					fwdMisMatch++;
					if (fwdMisMatch > maxMisMatches) break;
				}
			}
			if (rvFlag) {
				for (int i=stackPosition+1;i<oligoLength;i++) {
					if (seq[oligoLength-i-1] != 'N' && rvMap[(int)seq[oligoLength-i-1]] != otherSeq[i]) {
						rvMisMatch++;
						if (rvMisMatch > maxMisMatches) break;
					}
				}
			}

			if (fwdMisMatch <= maxMisMatches || (rvFlag && rvMisMatch <= maxMisMatches)) {
				vt += staticOligos[oligoIndex]->numTarget;
				vbg += staticOligos[oligoIndex]->numBackground;
				staticOligos[oligoIndex]->flag = 1;
				utilityOligos[cpu][numUtilityOligos[cpu]++] = staticOligos[oligoIndex];
//fprintf(stderr, "\t%s\t%lf\t%lf\n", otherSeq,vt,vbg);
				if (removeFlag == OLIGO_SEARCH_REMOVE) {
					tree[indexStack[cpu][stackPosition]+alphaPos] = 0;
				}
			}
			alphaStack[cpu][stackPosition]++;
			continue;
		} else if (valIndex == 0) {
			//dead end
			alphaStack[cpu][stackPosition]++;
			continue;
		} else {
			//next node
			alphaStack[cpu][stackPosition+1]=0;
			indexStack[cpu][stackPosition+1]=valIndex;
			fwdMisMatchStack[cpu][stackPosition+1]=fwdMisMatch;
			rvMisMatchStack[cpu][stackPosition+1]=rvMisMatch;

			alphaStack[cpu][stackPosition]++;
			stackPosition++;

			continue;
		}
	}
	return numUtilityOligos[cpu];
} 

void OligoArray::scoreOligosWithMotif(Motif* motif, double minimumScore, int removeFlag, int cpu) {

	motif->logxform();
	numUtilityOligos[cpu]=0;
	
	indexStack[cpu][0] = 0;
	alphaStack[cpu][0] = 0;
	positionStack[cpu][0] = oligoLength-1;
	fwdScoreStack[cpu][0] = 0.0;
	rvScoreStack[cpu][0] = 0.0;

	int stackPosition = 0;
	int valIndex = 0;
	int oligoIndex = 0;
	
	int alphaPos = 0;

	int fwdPosition=0;
	int rvPosition=0;
	double fwdScore=0.0;
	double rvScore=0.0;

	while (stackPosition > -1) {

		alphaPos = alphaStack[cpu][stackPosition];
		if (alphaPos >= alphaSize) {
			if (removeFlag == MOTIF_REMOVE_OLIGOS) {
				int numEmpty = 0;
				for (int i=0;i<alphaSize;i++) {
					if (tree[indexStack[cpu][stackPosition]+i] == 0) numEmpty++;
				}
				if (numEmpty == alphaSize && stackPosition > 0) {
					tree[indexStack[cpu][stackPosition-1]+alphaStack[cpu][stackPosition-1]]=0;
				}
			}
			stackPosition--;
			continue;
		}

		fwdPosition = stackPosition;
		rvPosition = oligoLength-fwdPosition-1;

		fwdScore = fwdScoreStack[cpu][stackPosition];
		rvScore = rvScoreStack[cpu][stackPosition];

		fwdScore += motif->m[fwdPosition][alphaPos];
		rvScore += motif->m[rvPosition][3-alphaPos];

		if (fwdPosition < oligoLength-1) {
			if (fwdScore+motif->fwdMaxPossible[fwdPosition+1] < minimumScore) {
				if (!rvFlag || rvScore+motif->rvMaxPossible[rvPosition-1] < minimumScore) {
					alphaStack[cpu][stackPosition]++;
					continue;
				}
			}
		}

		valIndex = tree[indexStack[cpu][stackPosition]+alphaPos];
		if (valIndex < 0) {
			oligoIndex = -valIndex;
			//this is a pointer to an oligo
			char* seq = staticOligos[oligoIndex]->seq;
			for (int i=fwdPosition+1;i<oligoLength;i++) {
				char bp = seq[i];
				fwdScore += motif->m[i][(int)alphaIndex[(int)bp]];
				rvScore += motif->m[oligoLength-i-1][3-alphaIndex[(int)bp]];
			}
			//fprintf(stderr, "%s\t%lf\t%lf\n", seq, fwdScore, rvScore);
			double bestScore = fwdScore;
			char direction = 0;
			if (fwdScore < rvScore && rvFlag) {
				bestScore = rvScore;
				direction = 1;
			}
			staticOligos[oligoIndex]->cpuData[cpu] = bestScore;
			staticOligos[oligoIndex]->value = bestScore;
			staticOligos[oligoIndex]->flag = direction;
	
			if (removeFlag == MOTIF_REMOVE_OLIGOS) {
				if (bestScore >= motif->threshold) {
					//fprintf(stderr, "\t\t\t%s\n", staticOligos[oligoIndex]->seq);
					tree[indexStack[cpu][stackPosition]+alphaPos] = 0;
				}
			}
			utilityOligos[cpu][numUtilityOligos[cpu]] = staticOligos[oligoIndex];
			numUtilityOligos[cpu]++;
			alphaStack[cpu][stackPosition]++;
			continue;
		} else if (valIndex == 0) {
			//dead end
			alphaStack[cpu][stackPosition]++;
			continue;
		} else {
			//next node
			alphaStack[cpu][stackPosition+1]=0;
			indexStack[cpu][stackPosition+1]=valIndex;
			fwdScoreStack[cpu][stackPosition+1]=fwdScore;
			rvScoreStack[cpu][stackPosition+1]=rvScore;

			alphaStack[cpu][stackPosition]++;

			stackPosition++;

			continue;
		}
	}
}

void OligoArray::calculateFoldEnrichment(double pseudoCount) {

	double ratio = 1;
	if (totalBackground > 0) {
		ratio = totalTarget / totalBackground;
	}
	double AAA = 0.0;
	fprintf(stderr, "\tCalculating Oligo Enrichment\n");
	fprintf(stderr, "\t\tTotal Target Oligos = %lf\n\t\tTotal Background Oligos = %lf\n", totalTarget, totalBackground);
	for (int i=0;i<numActiveOligos;i++) {
		double nt = activeOligos[i]->numTarget;
		double nb = activeOligos[i]->numBackground*ratio;
		AAA += activeOligos[i]->numBackground;
		if (nt < pseudoCount) nt = pseudoCount;
		if (nb < pseudoCount) nb = pseudoCount;
		activeOligos[i]->value = nt/nb;
	}
}

void OligoArray::normalizeOligos(int nlen,char* normFile) {
	if (nlen < 1) return;
	double lastError = 1e100;
	double cur = 0.0;
	double next = 1.0;
	for (int i=0;i<maxAutoNormIters;i++) {
		int printFlag = 0;
		if (i==0) printFlag = 1;
		double curError = normalizeOligosIteration(nlen,normFile,printFlag);
		double improvementRate = (lastError-curError)/lastError;
		//fprintf(stderr, "Round %d:\t%le\t%le\t%le\n", i+1,lastError,curError,improvementRate);
		if (curError >= lastError) break;
		if (improvementRate < MIN_OLIGO_NORMALIZATION_IMPROVEMENT) break;
		lastError = curError;
		if (printFlag) {
			fprintf(stderr, "\tConverging on autonormalization solution:\n\t");
		} else {
			cur += 80.0/(double)maxAutoNormIters;
			if (cur > next) {
				fprintf(stderr, ".");
				next += 1.0;
			}
		}
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "\n\tFinal normalization:\n");
	(void)normalizeOligosIteration(nlen,NULL,1);
}

double OligoArray::normalizeOligosIteration(int nlen,char* normFile, int printFlag) {

	FILE* fp = NULL;
	if (normFile != NULL && printFlag) {
		fp = fopen(normFile, "w");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for writing!!!\n",normFile);
		} else {
			fprintf(fp, "Oligo\tTargetCounts\tBackgroundCounts\tNormalizationFactor\n");
		}
	}

	double error = 0.0;

	char* mer = new char[1000];

	for (int curLen=1;curLen<=nlen;curLen++) {
		int total = (int) pow(4.0,(double)curLen);
		if (printFlag) fprintf(stderr, "\tAutonormalization: %d-mers (%d total)\n",curLen,total);
		double* targetLevels = new double[total];
		double* backgroundLevels = new double[total];
		double* normFactors = new double[total];
		int* pow4 = new int[curLen+1];
		pow4[0]=1;
		for (int i=1;i<=curLen;i++) pow4[i] = 4*pow4[i-1];
		int nbins = oligoLength-curLen+1;
		for (int i=0;i<total;i++) {
			targetLevels[i]=0.0;
			backgroundLevels[i]=0.0;
			normFactors[i] = 1.0;
		}

		double dnbins = (double)nbins;
		for (int i=1;i<numOligos+1;i++) {
			char* s = staticOligos[i]->seq;
			for (int j=0;j<nbins;j++) {
				int index = 0;
				for (int k=0;k<curLen;k++) {
					index += alphaIndex[(int)s[j+k]]*pow4[k];
				}
				targetLevels[index] += staticOligos[i]->numTarget/dnbins;
				backgroundLevels[index] += staticOligos[i]->numBackground/dnbins;
			}
		}

		double totalTarget = 0.0;
		double totalBackground = 0.0;
		for (int i=0;i<total;i++) {
			totalTarget += targetLevels[i];
			totalBackground += backgroundLevels[i];
		}
		double minimumBackgroundLevels = 0.5/totalBackground;
		double minimumTargetLevels = 0.5/totalTarget;
		for (int i=0;i<total;i++) {
			for (int j=0;j<curLen;j++) {
				int remainder = i % pow4[j+1];
				int asdf = (int)(remainder/pow4[j]);
				mer[j] = alpha[asdf];
			}
			mer[curLen]='\0';

			double tLevel = targetLevels[i];
			double bLevel = backgroundLevels[i];
			if (rvFlag) {
				int rindex = 0;
				for (int j=curLen-1;j>=0;j--) {
					rindex += (3-alphaIndex[(int)mer[j]])*pow4[curLen-j-1];
				}
				tLevel = (tLevel+targetLevels[rindex])/2.0;
				bLevel = (bLevel+backgroundLevels[rindex])/2.0;
			}

			if (bLevel < minimumBackgroundLevels) bLevel = minimumBackgroundLevels;
			if (tLevel < minimumTargetLevels) tLevel = minimumTargetLevels;

			double factor = totalBackground*tLevel/totalTarget/bLevel;
			if (printFlag) {
				if (curLen < 3) {
					fprintf(stderr, "\t\t%s\t%.2lf%%\t%.2lf%%\t%.3lf\n", mer,
									100.0*tLevel/totalTarget,
									100.0*bLevel/totalBackground,factor);
				}
				if (fp != NULL) {
					//fprintf(fp, "%s\t%.1lf\t%.1lf\t%.3lf\n", mer,tLevel,bLevel,factor);
					fprintf(fp, "\t\t%s\t%.2lf%%\t%.2lf%%\t%.3lf\n", mer,
									100.0*tLevel/totalTarget,
									100.0*bLevel/totalBackground,factor);
				}
			}
			normFactors[i] = factor;
			error += (factor-1.0)*(factor-1.0)/(double)total;
		}

		totalBackground = 0.0;
		for (int i=1;i<numOligos+1;i++) {
			char* s = staticOligos[i]->seq;
			double scoreUnit = (double)staticOligos[i]->numBackground;
			double newScore = 0;
			for (int j=0;j<nbins;j++) {
				int index = 0;
				for (int k=0;k<curLen;k++) {
					index += alphaIndex[(int)s[j+k]]*pow4[k];
				}
				newScore += normFactors[index];
			}
			newScore /= nbins;
			double score = (double)(newScore*scoreUnit);

			staticOligos[i]->numBackground = score;
			staticOligos[i]->normWeight *= newScore;
			totalBackground += staticOligos[i]->numBackground;
		}
		delete []targetLevels;
		delete []backgroundLevels;
		delete []normFactors;
	}
	if (fp != NULL) {
		fprintf(stderr, "\tNormalization weights can be found in file: %s\n", normFile);
		fclose(fp);
	}
	delete []mer;

	return error;
}
void OligoArray::initializeOligoPvalues(EnrichmentScorer* scorer) {
	//initialize individual oligo enrichments
	for (int i=1;i<numOligos+1;i++) {
		if (staticOligos[i] == NULL) {
			fprintf(stderr, "X");
			continue;
		}
		//fprintf(stderr, "- %x\n", activeOligos[i]);
		//activeOligos[i]->logp = -1*activeOligos[i]->numTarget;
		staticOligos[i]->value = 1.0;
		staticOligos[i]->flag2 = 0;
		staticOligos[i]->logp = (float)scorer->scoreEnrichment(staticOligos[i]->numTarget,
														staticOligos[i]->numBackground);
	}
}


// ----------------------- class PSSM -------------------------------

Motif::Motif() {
	init();
}
Motif::Motif(char* seq) {
	init();
	length = strlen(seq);
	consensus = new char[length+1];
	strcpy(consensus,seq);
	name = new char[length+1];
	strcpy(name,seq);
	m=new double*[length];
	for (int i=0;i<length;i++) {
		m[i]=new double[alphaSize];
		for (int j=0;j<alphaSize;j++) {
			m[i][j] = MOTIF_INITIAL_DEGENERACY;
			if (alphaIndex[(int)seq[i]] == j) {
				m[i][j] = 1.0 - MOTIF_INITIAL_DEGENERACY*3.0;
			}
		}
	}
	fwdMaxPossible = new double[length];
	rvMaxPossible = new double[length];
	for (int i=0;i<length;i++) {
		fwdMaxPossible[i]=0;
		rvMaxPossible[i]=0;
	}
	logFlag = 0;
}
Motif::Motif(Motif* motif) {
	init();
	motif->expxform();
	length = motif->length;
	m=new double*[length];
	for (int i=0;i<length;i++) {
		m[i]=new double[alphaSize];
		for (int j=0;j<alphaSize;j++) {
			m[i][j] = motif->m[i][j];
		}
	}
	fwdMaxPossible = new double[length];
	rvMaxPossible = new double[length];
	for (int i=0;i<length;i++) {
		fwdMaxPossible[i]=0;
		rvMaxPossible[i]=0;
	}
	logFlag = 0;
	numTarget = motif->numTarget;
	numBackground = motif->numBackground;
	freqTarget = motif->freqTarget;
	freqBackground = motif->freqBackground;
	normalize();
}
Motif::Motif(char* n, char* c, double t, double lp, double ** matrix, int len) {
	init();
	name = new char[strlen(n)+1];
	strcpy(name,n);
	consensus = new char[strlen(c)+1];
	strcpy(consensus,c);
	threshold = t; //-0.01;
	logp = lp;
	length = len;
	m = new double*[length];
	for (int i=0;i<length;i++) {
		m[i] = new double[alphaSize];
		for (int j=0;j<alphaSize;j++) {
			m[i][j] = matrix[i][j];
		}
	}
	logFlag = 0;
	//normalize();
}


Motif::~Motif() {
	if (m != NULL) {
		for (int i=0;i<length;i++) delete [](m[i]);
		delete []m;
	}
	if (seqStats != NULL) delete []seqStats;
	if (name != NULL) delete []name;
	if (consensus != NULL) delete []consensus;
	if (fwdMaxPossible != NULL) delete []fwdMaxPossible;
	if (rvMaxPossible != NULL) delete []rvMaxPossible;
}
void Motif::init() {
	m = NULL;
	name = NULL;
	consensus = NULL;
	length = 0;
	threshold = 0;
	logFlag = 0;
	numTarget = 0;
	numBackground = 0;
	freqTarget = 0;
	freqBackground = 0;
	seqStats = NULL;
	numSeqs = 0;
	logp=0;
	minimumScore = 0;
	badFlag =0;

	similarPercent = 0.0;
	similarPercentOfOther = 0.0;
	similarMotif = NULL;

	statAvgTarget = 0.0;
	statStdTarget = 0.0;
	statAvgBackground = 0.0;
	statStdBackground = 0.0;
	statStrandLogRatio = 0.0;
	siteMultiplicity = 0.0;
}
void Motif::setName(char* n) {
	if (name != NULL) delete []name;
	name = NULL;
	if (n == NULL) return;
	name = new char[strlen(n)+1];
	strcpy(name,n);
}
void Motif::addOligo(Oligo* oligo) {
	int rvFlag = oligo->flag;
	double value = -1*oligo->value;

	for (int i=0;i<length;i++) {
		int index = 0;
		if (rvFlag) index = 3-alphaIndex[(int)oligo->seq[length-i-1]];
		else index =alphaIndex[(int)oligo->seq[i]];
	
		m[i][index] += value;
	}
}
void Motif::initializeSeqStats(int num) {
	//if (seqStats != NULL) {
	if (seqStats != NULL && numSeqs != num) {
		delete []seqStats;
		seqStats = NULL;
	}
	numSeqs = num;
	if (seqStats == NULL) seqStats = new SequenceStats[numSeqs];
	for (int i=0;i<numSeqs;i++) {
		seqStats[i].initialize();
	}
}
void Motif::deleteSeqStats() {
	if (seqStats != NULL) delete []seqStats;
	seqStats = NULL;
}


void Motif::setBlank(int len) {
	length = len;
	m = new double*[length];
	for (int i=0;i<length;i++) {
		m[i] = new double[alphaSize];
		for (int j=0;j<alphaSize;j++) {
			m[i][j] = 0.0;
		}
	}
	fwdMaxPossible = new double[length];
	rvMaxPossible = new double[length];
	for (int i=0;i<length;i++) {
		fwdMaxPossible[i]=0;
		rvMaxPossible[i]=0;
	}
	logFlag = 0;
}

void Motif::logxform() {
	if (!logFlag) {
		for (int i=0;i<length;i++) {
			for (int j=0;j<alphaSize;j++) {
				if (m[i][j] < MOTIF_PROBABILITY_MINIMUM)
					m[i][j] = MOTIF_PROBABILITY_MINIMUM;
				m[i][j] = log(m[i][j]/0.25);
			}
		}
		logFlag =1;
	}
}
void Motif::expxform() {
	if (logFlag) {
		for (int i=0;i<length;i++) {
			for (int j=0;j<alphaSize;j++) {
				m[i][j] = exp(m[i][j])*0.25;
			}
		}
		logFlag = 0;
	}
	justNormalize();
}
void Motif::normalize() {
	justNormalize();
	logxform();
	expxform();
}
void Motif::justNormalize() {
	char tmp[1000];
	for (int i=0;i<length;i++) {
		double sum = 0;
		for (int j=0;j<alphaSize;j++) {
			sum += m[i][j];
		}
		for (int j=0;j<alphaSize;j++) {
			m[i][j] /= sum;
		}
		double totalerror = 0;
		double count = 0;
		for (int j=0;j<alphaSize;j++) {
			if (m[i][j] < MOTIF_PROBABILITY_MINIMUM) {
				totalerror += MOTIF_PROBABILITY_MINIMUM-m[i][j];
				count+=1.0;
			}
		}
		if (count > 0.0) {
			double adjust = totalerror/((double)alphaSize-count);
			for (int j=0;j<alphaSize;j++) {
				if (m[i][j] < MOTIF_PROBABILITY_MINIMUM) {
					m[i][j] = MOTIF_PROBABILITY_MINIMUM;
				} else {
					m[i][j] -= adjust;
				}
			}
		}
		for (int j=0;j<alphaSize;j++) {
			sprintf(tmp,"%.3lf",m[i][j]); // for rounding purposes
			sscanf(tmp,"%lf",&(m[i][j]));
		}
	}
}
void Motif::adjustLength(int newlen) {
	expxform();
	double** matrix = new double*[newlen];
	for (int i=0;i<newlen;i++) {
		matrix[i]=new double[alphaSize];
		for (int j=0;j<alphaSize;j++) matrix[i][j]=0.25;
	}
	if (length <= newlen) {
		int start = (newlen-length)/2;
		for (int i=start;i<length+start;i++) {
			for (int j=0;j<alphaSize;j++) matrix[i][j] = m[i-start][j];
		}
	} else {
		int start = (length-newlen)/2;
		for (int i=start;i<newlen+start;i++) {
			for (int j=0;j<alphaSize;j++) matrix[i-start][j] = m[i][j];
		}
	}
	for (int i=0;i<length;i++) delete [](m[i]);
	delete []m;
	m= matrix;
	length = newlen;

	if (fwdMaxPossible != NULL) delete []fwdMaxPossible;
	fwdMaxPossible = new double[length];
	for (int i=0;i<length;i++) fwdMaxPossible[i] = 0.0;

	if (rvMaxPossible != NULL) delete []rvMaxPossible;
	rvMaxPossible = new double[length];
	for (int i=0;i<length;i++) rvMaxPossible[i] = 0.0;
}



void Motif::print(FILE* fp) {
	if (fp == NULL) return;
	getConsensus();
	fprintf(fp, ">%s\t", consensus);
	if (name != NULL) {
		fprintf(fp, "%s", name);
	} else {
		fprintf(fp, "%s", consensus);
	}
	fprintf(fp, "\t%lf", threshold);
	fprintf(fp, "\t%lf", logp);
	fprintf(fp, "\t0");
	fprintf(fp, "\tT:%.1lf(%.2f%%),B:%.1lf(%.2f%%),P:1e%d",numTarget, freqTarget*100,  
									numBackground,freqBackground*100,(int)(logp/2.3026));
	if (siteMultiplicity > 0.5) {
		fprintf(fp,"\tTpos:%.1lf,Tstd:%.1lf,Bpos:%.1lf,Bstd:%.1lf,StrandBias:%.1lf,Multiplicity:%.2lf",
					statAvgTarget,statStdTarget,statAvgBackground,statStdBackground,
					statStrandLogRatio,siteMultiplicity);
	}
	if (similarMotif != NULL) {
		if (siteMultiplicity < 0.5) fprintf(fp, "\t");
		char* nn = similarMotif->name;
		if (nn == NULL) {
			if (similarMotif->consensus) similarMotif->getConsensus();
			nn = similarMotif->consensus;
		}
		fprintf(fp, "\tShares %.2lf%% of sites with %s(%.2lf%%)", similarPercent*100.0, nn,
														similarPercentOfOther*100.0);
	}
	fprintf(fp, "\n");
	for (int i=0;i<length;i++) {
		for (int j=0;j<alphaSize;j++) {
			if (j>0) fprintf(fp, "\t");
			fprintf(fp, "%.3lf",m[i][j]);
		}
		fprintf(fp, "\n");
	}
}

void Motif::getConsensus() {
	if (consensus != NULL) delete []consensus;
	consensus = new char[length+1];
	consensus[length] = '\0';
	expxform();
	char* alphaPresent = new char[alphaSize];
	for (int i=0;i<length;i++) {
		//int bestIndex = 0;
		double bestValue = 0;
		for (int j=0;j<alphaSize;j++) {
			alphaPresent[j]=0;
			if (m[i][j] > bestValue) {
				//bestIndex = j;
				bestValue = m[i][j];
			}
		}
		for (int j=0;j<alphaSize;j++) {
			double diff = bestValue - m[i][j];
			if (diff < MOTIF_IUPAC_THRESHOLD)
				alphaPresent[j]=1;
		}
		consensus[i] = getSymbol(alphaPresent);
	}
	delete []alphaPresent;
}
char Motif::getSymbol(char* p) {
	char bp = 'N';
	if (alphaSize == DNA_ALPHA) {
		if (p[0]) {
			if (p[1]) {
				if (p[2]) {
					if (p[3]) bp = 'N';
					else bp = 'V';
				} else {
					if (p[3]) bp = 'H';
					else bp = 'M';
				}
			} else {
				if (p[2]) {
					if (p[3]) bp = 'D';
					else bp = 'R';
				} else {
					if (p[3]) bp = 'W';
					else bp = 'A';
				}
			}
		} else {
			if (p[1]) {
				if (p[2]) {
					if (p[3]) bp = 'B';
					else bp = 'S';
				} else {
					if (p[3]) bp = 'Y';
					else bp = 'C';
				}
			} else {
				if (p[2]) {
					if (p[3]) bp = 'K';
					else bp = 'G';
				} else {
					if (p[3]) bp = 'T';
					else bp = 'N';
				}
			}
		}
	}
	return bp;
}

void Motif::setMinimumScore(double percent) {
	logxform();	
	int psize = (int) ceil(1.0/percent);
	if (psize < 2) {
		minimumScore = -1e60;
		return;
	}
	int psize4 = psize*alphaSize;
	double* lastResults=new double[psize];
	double* curResults=new double[psize4];
	for (int i=0;i<psize;i++) lastResults[i]=0.0;

	//set maximum scores, so you know when to quit
	for (int i=0;i<length;i++) {
		double max = 0;
		for (int j=0;j<alphaSize;j++) if (m[i][j] > max) max = m[i][j];
		rvMaxPossible[i]=max;
	}
	fwdMaxPossible[length-1]=rvMaxPossible[length-1];
	for (int i=length-2;i>=0;i--) fwdMaxPossible[i]=fwdMaxPossible[i+1]+rvMaxPossible[i];
	for (int i=1;i<length;i++) rvMaxPossible[i] += rvMaxPossible[i-1];

	//estimate percentile score
	for (int i=0;i<length;i++) {
		for (int j=0;j<psize;j++) {
			for (int k=0;k<alphaSize;k++) {
				curResults[j*alphaSize+k] = m[i][k]+lastResults[j];
			}
		}
		qsort(curResults, psize4, sizeof(double), &sortDouble);
		for (int j=0;j<psize;j++) {
			double sum = 0;
			for (int k=0;k<alphaSize;k++) {
				sum += curResults[j*alphaSize+k];
			}
			sum /= alphaSize;
			lastResults[j] = sum;
		}
	}
	minimumScore = lastResults[psize-2];
	//for (int i=0;i<psize;i++) {
	//	fprintf(stderr, "\t%lf\n", lastResults[i]);
	//}
	//fprintf(stderr, "\tMinimum Score set to %lf\n", minimumScore);
	delete []lastResults;
	delete []curResults;
}
 

int sortDouble(const void *a, const void* b) {
	if (*((double*)a) < *((double*)b)) return -1;
	if (*((double*)a) > *((double*)b)) return 1;
	return 0;
}
int Motif::nameID = 0;
int Motif::motifCmp(const void *a, const void *b) {
	double va=(*((Motif**)a))->logp;
	double vb=(*((Motif**)b))->logp;
	if (va < vb) return -1;
	if (va > vb) return 1;
	return 0;
}

// sscanf(s, "%lf", v): a field that is a plain decimal number (optional
// sign, digits with an optional point, an optional exponent with digits,
// nothing after) goes to strtod, which is what glibc's scanf converts such a
// field with; anything else goes to sscanf.
static void scanDouble(const char* s, double* v) {
	const char* p = s;
	if (*p == '+' || *p == '-') p++;
	int digits = 0;
	while (*p >= '0' && *p <= '9') {
		p++;
		digits++;
	}
	if (*p == '.') {
		p++;
		while (*p >= '0' && *p <= '9') {
			p++;
			digits++;
		}
	}
	if (digits > 0 && (*p == 'e' || *p == 'E')) {
		const char* e = p+1;
		if (*e == '+' || *e == '-') e++;
		int expDigits = 0;
		while (*e >= '0' && *e <= '9') {
			e++;
			expDigits++;
		}
		if (expDigits > 0) p = e;
		else digits = 0;
	}
	if (digits > 0 && *p == '\0') {
		*v = strtod(s, NULL);
		return;
	}
	sscanf(s, "%lf", v);
}

// The next line fgets(buf, MOTIF2_BUFFER, fp) would return from data[at, n):
// at most MOTIF2_BUFFER-1 bytes, ending after the first newline.
static int nextMotifLine(const std::vector<char>& data, size_t& at, const char*& line, size_t& len) {
	if (at >= data.size()) return 0;
	const size_t most = MOTIF2_BUFFER-1;
	size_t avail = data.size()-at;
	const char* p = data.data()+at;
	const char* nl = (const char*)memchr(p, '\n', avail < most ? avail : most);
	if (nl != NULL) len = nl-p+1;
	else len = avail < most ? avail : most;
	line = p;
	at += len;
	return 1;
}

Motif** Motif::readMotifFile(char* file, int &numMotifs) {

	numMotifs=0;
	if (file==NULL) return NULL;
	FILE* fp = fopen(file,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Couldn't open motif file: %s !!!\n", file);
		return NULL;
	}
	// the file once, in memory; both passes below read its fgets lines
	std::vector<char> data;
	{
		size_t got = 0;
		data.resize((size_t)1 << 20);
		while (1) {
			size_t r = fread(data.data()+got, 1, data.size()-got, fp);
			got += r;
			if (got < data.size()) break;
			data.resize(data.size()*2);
		}
		data.resize(got);
	}
	char* buf = new char[MOTIF2_BUFFER];
	char** cols = new char*[1000];
	int numCols;
	size_t at = 0;
	const char* line = NULL;
	size_t len = 0;

	int maxMotifs = 0;
	while (nextMotifLine(data, at, line, len)) {
		if (line[0] == '>') maxMotifs++;
	}

	int maxMotifLength = MAXIMUM_MOTIF_LENGTH;
	double** matrix = new double*[MAXIMUM_MOTIF_LENGTH];
	double* matrixRows = new double[(size_t)MAXIMUM_MOTIF_LENGTH*DNA_ALPHA];
	for (int i=0;i<maxMotifLength;i++) {
		matrix[i] = matrixRows + (size_t)i*DNA_ALPHA;
	}

	char* consensus = new char[MAXIMUM_MOTIF_LENGTH];
	char* name = new char[MAXIMUM_MOTIF_LENGTH];
	int curMotifLength = 0;
	consensus[0]='\0';
	name[0]='\0';
	double threshold = 0.0;
	double logp = 0.0;
	int started = 0;

	Motif** motifs = new Motif*[maxMotifs];

	// The second pass reads the file again from its start. A stream that
	// cannot seek (a pipe) gives it only what arrives after the first pass's
	// end of file, as rewind() + fgets() would.
	if (fseek(fp, 0, SEEK_SET) != 0) {
		clearerr(fp);
		size_t got = 0;
		data.resize((size_t)1 << 20);
		while (1) {
			size_t r = fread(data.data()+got, 1, data.size()-got, fp);
			got += r;
			if (got < data.size()) break;
			data.resize(data.size()*2);
		}
		data.resize(got);
	}
	at = 0;
	while (nextMotifLine(data, at, line, len)) {
		memcpy(buf, line, len);
		buf[len] = '\0';
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 3) continue;
		if (cols[0][0] == '>') {
			if (started) {
				motifs[numMotifs] = new Motif(name,consensus,threshold,
											logp,matrix,curMotifLength);
				numMotifs++;
			}
			curMotifLength = 0;
			if (numCols > 2) {
				started=1;
				strcpy(consensus,&(cols[0][1]));
				strcpy(name,cols[1]);
				scanDouble(cols[2],&threshold);
				logp = 0.0;
				if (numCols > 3) scanDouble(cols[3],&logp);
			} else {
				started=0;
			}
		} else if (numCols > 3) {
			for (int i=0;i<DNA_ALPHA;i++) {
				scanDouble(cols[i],&(matrix[curMotifLength][i]));
			}
			curMotifLength++;
		}
	}
	if (started) {
		motifs[numMotifs] = new Motif(name,consensus,threshold,
									logp,matrix,curMotifLength);
		numMotifs++;
	}
	
	fclose(fp);

	delete []buf;
	delete []cols;
	delete []consensus;
	delete []name;
 
	fprintf(stderr, "\t%d motifs loaded\n", numMotifs);

	//for (int i=0;i<numMotifs;i++) {
	//motifs[i]->print(stdout);
	//}
	//exit(0);


	return motifs;
}

void Motif::optimizeThreshold(EnrichmentScorer* scorer, double maxBackgroundPercent) {
	if (scorer == NULL) return;
	if (seqStats == NULL) return;
	qsort(seqStats, numSeqs, sizeof(SequenceStats),&(SequenceStats::scoreCmp)); 

	logp = 1.0;	
	double curNumTargets = 0;
	double curNumBackground = 0;
	double maxBackground =  scorer->numBackground*maxBackgroundPercent;
	for (int i=0;i<numSeqs-1;i++) {
		double curThresh = seqStats[i].score;
		double weight = seqStats[i].seq->w;
		if (seqStats[i].seq->g == 1) {
			curNumTargets+=weight;
		} else {
			curNumBackground+=weight;
		}
		if (curThresh < seqStats[i+1].score + 0.000001) {
			continue;
		}
		double curPvalue = scorer->scoreEnrichment(curNumTargets, curNumBackground);
		//fprintf(stderr, "%d\t%lf\t%lf\t%lf\t%lf\n", i,curThresh,curNumTargets,curNumBackground,curPvalue);
		if (curPvalue < logp) {
			logp = curPvalue;
			threshold = curThresh;
			numTarget = curNumTargets;
			numBackground = curNumBackground;
		}
		if (curNumBackground > maxBackground) break;
	}
	threshold -= 0.000001;
	freqTarget = numTarget/scorer->numTargets;
	freqBackground = numBackground/scorer->numBackground;
}

double Motif::scoreOverlap(Motif* m1, Motif* m2, double &percent1, double &percent2) {

	percent1 = 0;
	percent2 = 0;
	if (m1->numSeqs != m2->numSeqs) {
		fprintf(stderr, "!!! Different number of sequences for each motif !!!\n");
		exit(0);
	}
	qsort(m1->seqStats, m1->numSeqs, sizeof(SequenceStats),&(SequenceStats::defaultOrder)); 
	qsort(m2->seqStats, m2->numSeqs, sizeof(SequenceStats),&(SequenceStats::defaultOrder)); 

	int totalBp = 0;
	int totalSites1 = 0;
	int totalSites2 = 0;
	int overlap = 0;
	double expectedOverlap=0;

	int distance = 0;
	int maxMlen = m1->length;
	if (m2->length > maxMlen) maxMlen = m2->length;
	distance = (int) (((double)maxMlen)*MOTIF_OVERLAP_PERCENT);
	

	
	for (int i=0;i<m1->numSeqs;i++) {
		int seqlen = m1->seqStats[i].seq->length;
		int nsites1 = m1->seqStats[i].numSites;
		int nsites2 = m2->seqStats[i].numSites;
		totalBp += seqlen;
	
		int totalBp1 = nsites1*distance;
		int totalBp2 = nsites2*distance;
		int minSites = nsites1;
		if (nsites2 < minSites) minSites = nsites2;
		totalSites1 += nsites1;
		totalSites2 += nsites2;

		if (nsites1 == 0 || nsites2 == 0) continue;

		expectedOverlap += ((double)(totalBp1*totalBp2))/((double)seqlen*seqlen)*((double)(minSites));



		if (1) { //strand == MOTIF_STRAND_BOTH) {
			for (int j=0;j<nsites1;j++) {
				char d1 = m1->seqStats[i].strands[j];
				int p1 = m1->seqStats[i].pos[j];
				if (d1 == MOTIF_STRAND_POS) p1 += m1->length/2;
				else p1 -= m1->length/2;
				for (int k=0;k<nsites2;k++) {
					char d2 = m2->seqStats[i].strands[k];
					int p2 = m2->seqStats[i].pos[k];
					if (d2 == MOTIF_STRAND_NEG) p1 += m2->length/2;
					else p2 -= m2->length/2;

					if (p2-p1 > -distance && p2-p1 < distance) {
						overlap++;
					}
				}
			}
		} else { //if (strand != MOTIF_STRAND_BOTH) {
			//Strand specific code
			// The logic here is that sites are already ordered per sequence
			// first come + strand, in order i.e. 12+ 230+ 300+ then 309- 239- 34- 
			int startIndex2 = 0;
			for (int j=0;j<nsites1;j++) {
				char d1 = m1->seqStats[i].strands[j];
				int p1 = m1->seqStats[i].pos[j];
				for (int k=startIndex2;k<nsites2;k++) {
					char d2 = m2->seqStats[i].strands[k];
					int p2 = m2->seqStats[i].pos[k];
					if (d1 == MOTIF_STRAND_POS) {
						if (d2 == MOTIF_STRAND_NEG) {
							break;
						} else if (p2-p1 < -distance) {
							startIndex2 = k;
						} else if (p2-p1 < distance) {
							overlap++;
						} else {
							break;
						}
					} else { // d1 == MOTIF_STRAND_NEG
						if (d2 == MOTIF_STRAND_POS) {
							startIndex2 = k;
						} else if (p2-p1 > distance) {
							startIndex2 = k;
						} else if (p2-p1 > -distance) {
							overlap++;
						} else {
							break;
						}
					}
				}
			}
		}
	}
	int minNumSites = totalSites1;
	if (minNumSites > totalSites2) minNumSites = totalSites2;

	//fprintf(stderr, "%s = %d %s = %d Overlap = %d (%.2lf%%) (%.2lf%%) Expected = %lf\n",m1->name, totalSites1,m2->name, totalSites2, overlap, 100.0*overlap/totalSites1, 100.0*overlap/totalSites2,expectedOverlap);
	if (minNumSites == 0) return 0;

	percent1 = ((double)overlap)/(double)totalSites1;
	percent2 = ((double)overlap)/(double)totalSites2;

	double score = ((double)overlap)/(double)minNumSites;
	return score;
}

SequenceStats::SequenceStats() {
	initialize();
}
void SequenceStats::initialize() {
	seq = NULL;
	score = -1e100;
	numSites = 0;
	pos = NULL;
	strands = NULL;
}
SequenceStats::~SequenceStats() {
	if (pos != NULL) delete []pos;
	if (strands != NULL) delete []strands;
}
void SequenceStats::setSites(int* npos, char* nstrands, int nsites) {
	if (pos != NULL) {
		delete []pos;
		pos = NULL;
	}
	if (strands != NULL) {
		delete []strands;
		strands = NULL;
	}
	numSites = nsites;
	if (numSites < 1) return;
	pos = new int[numSites];
	strands = new char[numSites];
	for (int i=0;i<numSites;i++) {
		pos[i] = npos[i];
		strands[i] = nstrands[i];
	}
}

int SequenceStats::defaultOrder(const void *a, const void *b) {
	double va=((SequenceStats*)a)->seq->id;
	double vb=((SequenceStats*)b)->seq->id;
	if (va < vb) return 1;
	if (va > vb) return -1;
	return 0;
}
int SequenceStats::scoreCmp(const void *a, const void *b) {
	double va=((SequenceStats*)a)->score;
	double vb=((SequenceStats*)b)->score;
	if (va < vb) return 1;
	if (va > vb) return -1;
	return 0;
}



// ----------------------- class Oligo ------------------------------

int Oligo::cpu = 0;

Oligo::Oligo() {
	init();
}
void Oligo::init() {
	seq = NULL;
	numTarget = 0.0;
	numBackground = 0.0;
	value = 0.0;
	flag = 0;
	flag2 = 0;
	cpuData = NULL;
	normWeight = 1;
}
Oligo::~Oligo() {
	if (seq != NULL) {
		delete []seq;
	}
	if (cpuData != NULL) {
		delete []cpuData;
	}
}
Oligo::Oligo(char* s) {
	init();
	seq = new char[strlen(s)+1];
	strcpy(seq,s);
}
void Oligo::print(FILE* fp) {
	fprintf(fp, "%s\t%f\t%f\t%f\t%f\t%lf\t%d\t%d\n", seq, numTarget,numBackground,logp,normWeight,
				value,flag,flag2);
}
int Oligo::oligoCmp(const void *a, const void *b) {
	double va=(*((Oligo**)a))->value;
	double vb=(*((Oligo**)b))->value;
	if (va < vb) return -1;
	if (va > vb) return 1;
	return 0;
}
int Oligo::oligoCmpReverse(const void *a, const void *b) {
	double va=(*((Oligo**)a))->value;
	double vb=(*((Oligo**)b))->value;
	if (va < vb) return 1;
	if (va > vb) return -1;
	return 0;
}
int Oligo::oligoCmpReverseParallel(const void *a, const void *b) {
	double va=(*((Oligo**)a))->cpuData[cpu];
	double vb=(*((Oligo**)b))->cpuData[cpu];
	if (va < vb) return 1;
	if (va > vb) return -1;
	return 0;
}

void Oligo::revopp(char* seq) {
    if (seq == NULL) return;
    int len = strlen(seq);
    int mid = (int)((((double)len)/2.0)+0.5);
    for (int i=0;i<mid;i++) {
        int rIndex = len-1-i;
        char ogHere = seq[i];
        char ogThere = seq[rIndex];
        char there;
        char here;

        if (ogHere == 'A' || ogHere == 'a') {
            there = 'T';
        } else if (ogHere == 'C' || ogHere == 'c') {
            there = 'G';
        } else if (ogHere == 'G' || ogHere == 'g') {
            there = 'C';
        } else if (ogHere == 'T' || ogHere == 't') {
            there = 'A';
        } else {
            there = 'N';
        }
        if (ogThere == 'A' || ogThere == 'a') {
            here = 'T';
        } else if (ogThere == 'C' || ogThere == 'c') {
            here = 'G';
        } else if (ogThere == 'G' || ogThere == 'g') {
            here = 'C';
        } else if (ogThere == 'T' || ogThere == 't') {
            here = 'A';
        } else {
            here = 'N';
        }
		seq[i] = here;
        seq[rIndex] = there;
    }
}


//==================================================================================
Sequence::Sequence() {
	init();
}
Sequence::Sequence(char* n, char* seq, int group, double weight) {
	init();
	setSeq(n,seq,group,weight);
}
Sequence::~Sequence() {
	if (textShared) {
		s = NULL;
		name = NULL;
	}
	if (s != NULL) delete []s;
	if (name != NULL) delete []name;
	if (kmer5 != NULL) delete []kmer5;
	if (kmer6 != NULL) delete []kmer6;
	if (codes != NULL && !codesShared) delete []codes;
}
void Sequence::init() {
	s=NULL;
	name=NULL;
	kmer5=NULL;
	kmer6=NULL;
	codes=NULL;
	codesShared=0;
	textShared=0;
	length = 0;
	w=1.0;
	g=0;
	id=0;
}
void Sequence::setSeq(char* n, char* seq, int group, double weight) {
	if (textShared) {
		s = NULL;
		name = NULL;
		textShared = 0;
	}
	if (name != NULL) delete []name;
	if (n == NULL) {
		name = NULL;
	} else {
		name = new char[strlen(n)+1];
		strcpy(name,n);
	}
	w = weight;
	g = group;
	if (kmer5 != NULL) delete []kmer5;
	kmer5 = NULL;
	if (kmer6 != NULL) delete []kmer6;
	kmer6 = NULL;
	if (codes != NULL && !codesShared) delete []codes;
	codes = NULL;
	codesShared = 0;
	if (s != NULL) delete []s;
	if (seq == NULL) {
		s = NULL;
		length = 0;
	} else {
		length = strlen(seq);
		s = new char[length+1];
		strcpy(s,seq);
	}
}

// Base -> column of the flattened log-odds table: A,C,G,T = 0..3, and every
// other character (N, IUPAC letters, lowercase) = 4, whose entry is -inf.
static unsigned char motifScanCode[256];
static int motifScanCodeInit() {
	for (int i=0;i<256;i++) motifScanCode[i] = 4;
	motifScanCode[(unsigned char)'A'] = 0;
	motifScanCode[(unsigned char)'C'] = 1;
	motifScanCode[(unsigned char)'G'] = 2;
	motifScanCode[(unsigned char)'T'] = 3;
	return 1;
}
static int motifScanCodeReady = motifScanCodeInit();

// find and -mscore write each motif's lines in motif-file order, the order
// HOMER v5.1's homer2 gives at -p 1, whatever -p is. Lines of the first motif in
// file order not yet complete (the head) join the queue of lines ready to
// write; any other motif keeps its lines until every motif before it is
// complete. One thread at a time writes the queue, outside the lock. Threads
// take motifs in file order, which keeps little unwritten; with more than
// orderedCap bytes unwritten they take only the head, or wait for it.
static const size_t orderedCap = (size_t)64 << 20;
// Blocks smaller than this are copied; larger ones are taken over with the
// thread's buffer, which gets a spare back.
static const size_t orderedCopyMax = (size_t)1 << 16;
// A thread's output buffer: allocated whole (chunkArray, released
// with free()), flushed when it holds more than 1 MB.
static const size_t scanOutBlock = ((size_t)1 << 20) + ((size_t)64 << 10);
template <typename T>
static T* chunkArray(long n);
struct OrderedPiece {
	char* data;
	size_t len;
	size_t cap;
};
struct OrderedOutput {
	FILE* fp;
	int n;
	int fileCursor;
	int next;
	std::vector<char> started;
	std::vector<char> done;
	std::vector<std::vector<OrderedPiece> > held;
	std::vector<OrderedPiece> ready;
	std::vector<OrderedPiece> pool;
	size_t unwritten;
	bool writing;
	std::mutex lock;
	std::condition_variable advanced;

	// estimated cost per motif, and the start-early rule's settings
	std::vector<double> weight;
	double heavyMin;
	int lookAhead;

	OrderedOutput(FILE* f, int count) : fp(f), n(count), fileCursor(0), next(0),
			started(count, 0), done(count, 0), held(count), unwritten(0), writing(false),
			heavyMin(INFINITY), lookAhead(64) {}
	// a motif at 8x the mean estimate or more is heavy
	void setWeights(const std::vector<double>& w) {
		weight = w;
		double sum = 0.0;
		for (size_t k=0;k<w.size();k++) sum += w[k];
		const double mean = w.empty() ? 0.0 : sum/(double)w.size();
		heavyMin = (mean > 0.0) ? 8.0*mean : INFINITY;
	}
	~OrderedOutput() {
		for (size_t j=0;j<pool.size();j++) free(pool[j].data);
	}

	// The next motif to scan. Called once per motif, so one is always left.
	int take() {
		std::unique_lock<std::mutex> g(lock);
		while (1) {
			while (started[fileCursor]) fileCursor++;
			// When the head has not started, fileCursor is the head.
			if (unwritten <= orderedCap || !started[next]) {
				int k = fileCursor;
				// with the head under way, the heaviest motif within the
				// look-ahead goes first, so that it does not run alone at the
				// end; only its own output waits for its turn
				if (started[next] && !weight.empty()) {
					const int end = (fileCursor + lookAhead < n) ? fileCursor + lookAhead : n;
					double best = heavyMin;
					for (int j=fileCursor;j<end;j++) {
						if (!started[j] && weight[j] >= best) {
							best = weight[j];
							k = j;
						}
					}
				}
				started[k] = 1;
				return k;
			}
			// The head is running on another thread; its finish() wakes us.
			advanced.wait(g);
		}
	}
	// Writes the queue until it is empty, unless another thread is already
	// writing it. Called with g locked; returns with g locked.
	void drain(std::unique_lock<std::mutex>& g) {
		if (writing) return;
		writing = true;
		std::vector<OrderedPiece> batch;
		while (!ready.empty()) {
			batch.swap(ready);
			g.unlock();
			for (size_t j=0;j<batch.size();j++) fwrite(batch[j].data, 1, batch[j].len, fp);
			g.lock();
			for (size_t j=0;j<batch.size();j++) {
				unwritten -= batch[j].len;
				if (batch[j].cap > 0) {
					pool.push_back(batch[j]);
				} else {
					delete [](batch[j].data);
				}
			}
			batch.clear();
		}
		writing = false;
		advanced.notify_all();
	}
	void promote(int k) {
		std::vector<OrderedPiece>& h = held[k];
		ready.insert(ready.end(), h.begin(), h.end());
		std::vector<OrderedPiece>().swap(h);
	}
	// Lines of motif k in buf[0, len).
	void emit(int k, char*& buf, size_t len, size_t& cap) {
		OrderedPiece piece = {NULL, len, 0};
		if (len < orderedCopyMax) {
			piece.data = new char[len];
			memcpy(piece.data, buf, len);
		}
		std::unique_lock<std::mutex> g(lock);
		if (piece.data == NULL) {
			piece.data = buf;
			piece.cap = cap;
			if (pool.empty()) {
				buf = NULL;
				cap = 0;
			} else {
				buf = pool.back().data;
				cap = pool.back().cap;
				pool.pop_back();
			}
		}
		unwritten += len;
		if (k == next) {
			ready.push_back(piece);
			drain(g);
		} else {
			held[k].push_back(piece);
		}
	}
	void finish(int k) {
		std::unique_lock<std::mutex> g(lock);
		done[k] = 1;
		while (next < n && done[next]) {
			promote(next);
			next++;
		}
		// The new head may be running: what it has so far joins the queue,
		// the rest comes through emit().
		if (next < n) promote(next);
		advanced.notify_all();
		drain(g);
	}
};
static OrderedOutput* orderedOut = NULL;
// The index entries a motif's scan will visit, estimated as the pivot choice
// in MotifScanTables::buildSix does: the normal approximation of a 7-mer
// passing (base composition from the index's 7-mer counts), at the best
// pivot window, both strands, times the 7-mers indexed; 0 for a motif the
// 7-mer index does not serve. It orders the scan only.
// comp: the 7-mers' first-base composition, total: the 7-mers
// indexed (both from count7, computed once by the caller).
static double motifScanCost(Motif* motif, const double* comp, double total) {
	const int mlen = motif->length;
	if (mlen < 7 || total <= 0.0) return 0.0;
	// log-odds as logxform makes them, leaving the motif as it is
	std::vector<double> lo((size_t)mlen*4);
	for (int j=0;j<mlen;j++) {
		for (int b=0;b<4;b++) {
			double v = motif->m[j][b];
			if (!motif->logFlag) {
				if (v < MOTIF_PROBABILITY_MINIMUM) v = MOTIF_PROBABILITY_MINIMUM;
				v = log(v/0.25);
			}
			lo[j*4+b] = v;
		}
	}
	double allBest = 0.0;
	std::vector<double> best((size_t)mlen);
	for (int j=0;j<mlen;j++) {
		double v = lo[j*4];
		for (int b=1;b<4;b++) if (lo[j*4+b] > v) v = lo[j*4+b];
		best[j] = v;
		allBest += v;
	}
	const double cutEst = motif->threshold - 1e-6;
	double bestP = -1.0;
	for (int a=0;a+7<=mlen;a++) {
		double inside = 0.0;
		for (int t=0;t<7;t++) inside += best[a+t];
		const double need = cutEst - (allBest - inside);
		double pass = 0.0;
		for (int st=0;st<2;st++) {
			double mu = 0.0;
			double var = 0.0;
			for (int t=0;t<7;t++) {
				double m1 = 0.0;
				double m2 = 0.0;
				for (int b=0;b<4;b++) {
					const double v = (st == 0) ? lo[(a+t)*4+b] : lo[(a+t)*4+3-b];
					m1 += comp[b]*v;
					m2 += comp[b]*v*v;
				}
				mu += m1;
				var += m2 - m1*m1;
			}
			const double sd = sqrt(var > 0.0 ? var : 0.0) + 1e-12;
			pass += 0.5*erfc((need - mu)/sd/sqrt(2.0));
		}
		if (bestP < 0.0 || pass < bestP) bestP = pass;
	}
	return (bestP > 0.0) ? bestP*total : 0.0;
}
// The motif this thread is scanning for orderedOut, or -1.
static thread_local int orderedMotif = -1;

// Per-thread output for find and -mscore: lines are formatted here and
// written to the shared FILE in blocks, so threads do not take its lock once
// per line.
static thread_local char* scanOut = NULL;
static thread_local size_t scanOutLen = 0;
static thread_local size_t scanOutCap = 0;
static void scanOutFlush(FILE* fp) {
	if (scanOutLen > 0 && fp != NULL) {
		if (orderedMotif >= 0) {
			orderedOut->emit(orderedMotif, scanOut, scanOutLen, scanOutCap);
		} else {
			fwrite(scanOut, 1, scanOutLen, fp);
		}
	}
	scanOutLen = 0;
}
// printf("%.6lf", x) without printf, for the values a motif score takes.
// Writes into out (at least 512 bytes) and returns the length. Values that
// are not finite, at least 1e6 in magnitude, or within 1e-3 of a rounding tie
// in the sixth decimal go through snprintf, so the bytes are printf's.
static int formatFixed6(double x, char* out) {
	double ax = fabs(x);
	if (!(ax < 1e6)) return snprintf(out, 512, "%.6lf", x);
	// x*1e6 is within 1.2e-4 of the exact product for |x| < 1e6, so away
	// from a tie the nearest integer is the one printf rounds to.
	double y = ax * 1e6;
	double fl = floor(y);
	double frac = y - fl;
	if (frac > 0.499 && frac < 0.501) return snprintf(out, 512, "%.6lf", x);
	unsigned long long r = (unsigned long long)fl + (frac > 0.5 ? 1 : 0);
	unsigned long long ip = r / 1000000ULL;
	unsigned int fp6 = (unsigned int)(r % 1000000ULL);
	char tmp[24];
	int n = 0;
	do {
		tmp[n++] = (char)('0' + ip % 10);
		ip /= 10;
	} while (ip > 0);
	int len = 0;
	if (signbit(x)) out[len++] = '-';
	while (n > 0) out[len++] = tmp[--n];
	out[len++] = '.';
	for (int d=5;d>=0;d--) {
		out[len+d] = (char)('0' + fp6 % 10);
		fp6 /= 10;
	}
	len += 6;
	out[len] = '\0';
	return len;
}

// One site line, `name\tpos\tsite\tmotif\tstrand\tscore\n`, the bytes
// "%s\t%d\t%.*s\t%s\t%c\t%.6lf\n" gives; siteLen < 0 means site is a C string.
static void scanOutSite(FILE* fp, const char* name, int pos, const char* site, int siteLen,
		const char* motifName, char strand, double score) {
	size_t nameLen = strlen(name);
	size_t sLen = siteLen < 0 ? strlen(site) : strnlen(site, (size_t)siteLen);
	size_t mLen = strlen(motifName);
	size_t need = nameLen + sLen + mLen + 16 + 24 + 512;
	if (scanOutCap - scanOutLen < need) {
		size_t cap = scanOutBlock;
		if (cap < 2*scanOutCap) cap = 2*scanOutCap;
		if (cap < scanOutLen + need) cap = scanOutLen + need;
		char* grown = chunkArray<char>((long)cap);
		if (scanOutLen > 0) memcpy(grown, scanOut, scanOutLen);
		free(scanOut);
		scanOut = grown;
		scanOutCap = cap;
	}
	char* o = scanOut + scanOutLen;
	memcpy(o, name, nameLen);
	o += nameLen;
	*o++ = '\t';
	long long p = pos;
	if (p < 0) {
		*o++ = '-';
		p = -p;
	}
	char digits[24];
	int nd = 0;
	do {
		digits[nd++] = (char)('0' + p % 10);
		p /= 10;
	} while (p > 0);
	while (nd > 0) *o++ = digits[--nd];
	*o++ = '\t';
	memcpy(o, site, sLen);
	o += sLen;
	*o++ = '\t';
	memcpy(o, motifName, mLen);
	o += mLen;
	*o++ = '\t';
	*o++ = strand;
	*o++ = '\t';
	o += formatFixed6(score, o);
	*o++ = '\n';
	scanOutLen = o - scanOut;
	if (scanOutLen > (1<<20)) scanOutFlush(fp);
}

MotifScanTables::MotifScanTables() {
	length = 0;
	K = 0;
	restOfBest = 0.0;
	fwd = NULL;
	rev = NULL;
	sel = NULL;
	gap = NULL;
	margin = 0.0;
	size = 0;
	nblocks = 0;
	blockStart[0] = 0;
	blockStart[1] = 0;
	blockLen[0] = 0;
	blockLen[1] = 0;
	blockRest = 0.0;
	bfwd = new double[2*3125];
	brev = new double[2*3125];
	indexable = 0;
	pivot[0] = 0;
	pivot[1] = 0;
	npass[0] = 0;
	npass[1] = 0;
	pass[0] = new int[3125];
	pass[1] = new int[3125];
	six = 0;
	sixP = 6;
	s6n = 0;
	s6Rest = 0.0;
	for (int st=0;st<2;st++) {
		for (int q=0;q<3;q++) {
			t6[st][q] = NULL;
			o5[st][q] = NULL;
			o5x[st][q] = NULL;
			s6max[st][q] = 0.0;
		}
		s6pass[st] = NULL;
		s6passScore[st] = NULL;
		s6passGap[st] = NULL;
		s6npass[st] = 0;
		s6pivot[st] = 0;
		s6oth[st][0] = 1;
		s6oth[st][1] = 2;
	}
	for (int q=0;q<3;q++) {
		s6Start[q] = 0;
		s6Len[q] = 0;
		o5Width[q] = 0;
		o5Unused[q] = 1;
	}
}
MotifScanTables::~MotifScanTables() {
	delete []fwd;
	delete []rev;
	delete []sel;
	delete []gap;
	delete []bfwd;
	delete []brev;
	delete []pass[0];
	delete []pass[1];
	for (int st=0;st<2;st++) {
		for (int q=0;q<3;q++) {
			delete [](t6[st][q]);
			delete [](o5[st][q]);
			delete [](o5x[st][q]);
		}
		delete [](s6pass[st]);
		delete [](s6passScore[st]);
		delete [](s6passGap[st]);
	}
}
void MotifScanTables::build(Motif* motif, int blocks) {
	int mlen = motif->length;
	length = mlen;
	if (size < mlen+1) {
		delete []fwd;
		delete []rev;
		delete []sel;
		delete []gap;
		size = mlen+1;
		fwd = new double[size*5];
		rev = new double[size*5];
		sel = new int[size];
		gap = new double[size];
	}
	for (int j=0;j<mlen;j++) {
		for (int k=0;k<4;k++) {
			fwd[j*5+k] = motif->m[j][k];
			rev[j*5+k] = motif->m[j][3-k];
		}
		fwd[j*5+4] = -INFINITY;
		rev[j*5+4] = -INFINITY;
	}

	// The bound. For each position, the best column; positions are ranked by
	// how far a uniform base falls below that best (largest first), and the
	// first K are scored exactly for every window, the rest at their best.
	for (int j=0;j<mlen;j++) {
		double best = motif->m[j][0];
		double mean = 0.0;
		for (int k=0;k<4;k++) {
			if (motif->m[j][k] > best) best = motif->m[j][k];
			mean += motif->m[j][k];
		}
		gap[j] = best - mean*0.25;
		sel[j] = j;
	}
	for (int a=1;a<mlen;a++) {
		int x = sel[a];
		int b = a-1;
		while (b >= 0 && gap[sel[b]] < gap[x]) {
			sel[b+1] = sel[b];
			b--;
		}
		sel[b+1] = x;
	}
	// Four positions: in a sweep of K = 3 to 8 and an adaptive rule on 5,000
	// ATAC peaks x 472 known motifs, K = 4 and 5 were fastest.
	K = mlen < 4 ? mlen : 4;
	restOfBest = 0.0;
	for (int t=K;t<mlen;t++) {
		double best = motif->m[sel[t]][0];
		for (int k=1;k<4;k++) if (motif->m[sel[t]][k] > best) best = motif->m[sel[t]][k];
		restOfBest += best;
	}
	// a NaN in a first column or an infinite entry gives a margin that
	// filters nothing
	margin = 0.0;
	for (int j=0;j<mlen;j++) {
		double big = fabs(motif->m[j][0]);
		for (int k=1;k<4;k++) if (fabs(motif->m[j][k]) > big) big = fabs(motif->m[j][k]);
		margin += big;
	}
	margin = 1e-9*(1.0+margin);

	// find's blocks: the windows of the motif with the largest summed gap,
	// two that do not overlap (5+5 from width 10, 5+4 or 4+5 at 9, 4+4 at 8),
	// or one (5 for widths 5-7, 4 at width 4).
	nblocks = 0;
	if (!blocks) {
		indexable = 0;
		npass[0] = 0;
		npass[1] = 0;
		return;
	}
	{
		int options[3][2];
		int numOptions = 0;
		if (mlen >= 10) {
			options[numOptions][0] = 5; options[numOptions][1] = 5; numOptions++;
		} else if (mlen == 9) {
			options[numOptions][0] = 5; options[numOptions][1] = 4; numOptions++;
			options[numOptions][0] = 4; options[numOptions][1] = 5; numOptions++;
		} else if (mlen == 8) {
			options[numOptions][0] = 4; options[numOptions][1] = 4; numOptions++;
		}
		double bestPair = -1.0;
		for (int o=0;o<numOptions;o++) {
			int k0 = options[o][0];
			int k1 = options[o][1];
			for (int a=0;a+k0<=mlen;a++) {
				double ga = 0.0;
				for (int t=0;t<k0;t++) ga += gap[a+t];
				for (int b=a+k0;b+k1<=mlen;b++) {
					double gb = 0.0;
					for (int t=0;t<k1;t++) gb += gap[b+t];
					if (ga+gb > bestPair) {
						bestPair = ga+gb;
						blockStart[0] = a;
						blockLen[0] = k0;
						blockStart[1] = b;
						blockLen[1] = k1;
						nblocks = 2;
					}
				}
			}
		}
		if (numOptions == 0 && mlen >= 4) {
			int k = (mlen >= 5) ? 5 : 4;
			for (int a=0;a+k<=mlen;a++) {
				double ga = 0.0;
				for (int t=0;t<k;t++) ga += gap[a+t];
				if (ga > bestPair) {
					bestPair = ga;
					blockStart[0] = a;
					blockLen[0] = k;
					blockStart[1] = a;
					blockLen[1] = k;
					nblocks = 1;
				}
			}
		}
	}
	if (nblocks > 0) {
		blockRest = 0.0;
		for (int j=0;j<mlen;j++) {
			int inBlock = 0;
			for (int q=0;q<nblocks;q++) {
				if (j >= blockStart[q] && j < blockStart[q]+blockLen[q]) inBlock = 1;
			}
			if (inBlock) continue;
			double best = motif->m[j][0];
			for (int k=1;k<4;k++) if (motif->m[j][k] > best) best = motif->m[j][k];
			blockRest += best;
		}
		// Forward: motif positions a..a+k-1 read the 5-mer in order; reverse
		// complement: the 5-mer at i-a-k+1 holds positions a+k-1, ..., a of
		// the window that ends at i. Built by prefix sums over the first t
		// digits (the column for code 4 is -inf, so any N gives -inf); a
		// block of 4 repeats its 4-digit sums over the fifth digit.
		double pf[3125+625+125+25+5];
		double pr[3125+625+125+25+5];
		for (int q=0;q<2;q++) {
			if (q >= nblocks) {
				for (int idx=0;idx<3125;idx++) {
					bfwd[q*3125+idx] = 0.0;
					brev[q*3125+idx] = 0.0;
				}
				continue;
			}
			int a = blockStart[q];
			int k = blockLen[q];
			double* lf = pf;
			double* lr = pr;
			for (int d=0;d<5;d++) {
				lf[d] = 0.0 + fwd[a*5+d];
				lr[d] = 0.0 + rev[(a+k-1)*5+d];
			}
			int width = 5;
			for (int t=1;t<k;t++) {
				double* nf = lf + width;
				double* nr = lr + width;
				for (int x=0;x<width;x++) {
					for (int d=0;d<5;d++) {
						nf[x*5+d] = lf[x] + fwd[(a+t)*5+d];
						nr[x*5+d] = lr[x] + rev[(a+k-1-t)*5+d];
					}
				}
				lf = nf;
				lr = nr;
				width *= 5;
			}
			for (int idx=0;idx<3125;idx++) {
				int x = (k == 5) ? idx : idx/5;
				bfwd[q*3125+idx] = lf[x];
				brev[q*3125+idx] = lr[x];
			}
		}
	}

	// find's index. A window can survive the bound only if its pivot-block
	// 4-mer scores at least cut - blockRest - (the other block's best); the
	// 1e-9 slack covers rounding, so every window the bound keeps is visited.
	indexable = 0;
	npass[0] = 0;
	npass[1] = 0;
	if (nblocks > 0) {
		indexable = 1;
		for (int x=0;x<2*3125;x++) {
			if (bfwd[x] != bfwd[x] || brev[x] != brev[x]) indexable = 0;
		}
	}
	if (indexable) {
		const double cut = motif->threshold - 1e-6;
		for (int st=0;st<2;st++) {
			const double* T = (st == 0) ? bfwd : brev;
			double need[2] = {0.0, 0.0};
			int count[2] = {0, 0};
			for (int q=0;q<nblocks;q++) {
				const double* To = T + (1-q)*3125;
				double maxOther = -INFINITY;
				for (int x=0;x<3125;x++) if (To[x] > maxOther) maxOther = To[x];
				need[q] = cut - blockRest - maxOther - 1e-9;
				for (int x=0;x<3125;x++) if (!(T[q*3125+x] < need[q])) count[q]++;
			}
			int q = (nblocks == 2 && count[1] < count[0]) ? 1 : 0;
			pivot[st] = q;
			npass[st] = 0;
			for (int x=0;x<3125;x++) {
				if (!(T[q*3125+x] < need[q])) pass[st][npass[st]++] = x;
			}
		}
	}
}

int MotifScanTables::buildSix(Motif* motif, int P, const int* countP) {
	six = 0;
	sixP = P;
	const int mlen = length;
	if (mlen < P) return 0;
	s6n = 1;
	for (int q=0;q<3;q++) {
		s6Start[q] = 0;
		s6Len[q] = 0;
	}
	s6Len[0] = P;
	if (P == 7 && countP != NULL && mlen > P) {
		// The pivot: the window of 7 least likely to pass, from a normal
		// approximation of its score over random 7-mers with the sequences'
		// base composition (per-position mean and variance), both strands. A
		// window's 7-mer can pass only if it scores at least cut - (the best
		// of every position outside it), whatever the other blocks are. The
		// other blocks then tighten the bound: the largest-gap window of up
		// to 7 positions on each side of the pivot.
		static thread_local const int* compFor = NULL;
		static thread_local double comp[4];
		if (compFor != countP) {
			double c4[4] = {0.0, 0.0, 0.0, 0.0};
			for (int x=0;x<16384;x++) c4[x >> 12] += countP[x];
			const double sum = c4[0]+c4[1]+c4[2]+c4[3];
			for (int b=0;b<4;b++) comp[b] = (sum > 0.0) ? c4[b]/sum : 0.25;
			compFor = countP;
		}
		double allBest = 0.0;
		for (int j=0;j<mlen;j++) {
			double best = motif->m[j][0];
			for (int k=1;k<4;k++) if (motif->m[j][k] > best) best = motif->m[j][k];
			allBest += best;
		}
		const double cutEst = motif->threshold - 1e-6;
		double bestP = -1.0;
		int bestA = 0;
		for (int a=0;a+7<=mlen;a++) {
			double inside = 0.0;
			for (int t=0;t<7;t++) {
				double best = motif->m[a+t][0];
				for (int k=1;k<4;k++) if (motif->m[a+t][k] > best) best = motif->m[a+t][k];
				inside += best;
			}
			const double need = cutEst - (allBest - inside);
			double pass = 0.0;
			for (int st=0;st<2;st++) {
				double mu = 0.0;
				double var = 0.0;
				for (int t=0;t<7;t++) {
					const double* row = motif->m[a+t];
					double m1 = 0.0;
					double m2 = 0.0;
					for (int b=0;b<4;b++) {
						const double v = (st == 0) ? row[b] : row[3-b];
						m1 += comp[b]*v;
						m2 += comp[b]*v*v;
					}
					mu += m1;
					var += m2 - m1*m1;
				}
				const double sd = sqrt(var > 0.0 ? var : 0.0) + 1e-12;
				pass += 0.5*erfc((need - mu)/sd/sqrt(2.0));
			}
			if (bestP < 0.0 || pass < bestP) {
				bestP = pass;
				bestA = a;
			}
		}
		s6Start[0] = bestA;
		const int sides[2][2] = {{0, bestA}, {bestA+7, mlen}};
		for (int sd=0;sd<2;sd++) {
			const int from = sides[sd][0];
			const int to = sides[sd][1];
			const int k = (to-from < 7) ? to-from : 7;
			if (k <= 0) continue;
			double bestGap = -1.0;
			int bestB = from;
			for (int b=from;b+k<=to;b++) {
				double g = 0.0;
				for (int t=0;t<k;t++) g += gap[b+t];
				if (g > bestGap) {
					bestGap = g;
					bestB = b;
				}
			}
			s6Start[s6n] = bestB;
			s6Len[s6n] = k;
			s6n++;
		}
	} else if (mlen >= 2*P) {
		// the two windows of P with the largest summed gap
		double bestPair = -1.0;
		for (int a=0;a+P<=mlen;a++) {
			double ga = 0.0;
			for (int t=0;t<P;t++) ga += gap[a+t];
			for (int b=a+P;b+P<=mlen;b++) {
				double gb = 0.0;
				for (int t=0;t<P;t++) gb += gap[b+t];
				if (ga+gb > bestPair) {
					bestPair = ga+gb;
					s6Start[0] = a;
					s6Start[1] = b;
				}
			}
		}
		s6n = 2;
		s6Len[1] = P;
	} else if (mlen > P) {
		// a block of P at either end, whichever has the larger gap, and the rest
		double g0 = 0.0;
		double g1 = 0.0;
		for (int t=0;t<P;t++) {
			g0 += gap[t];
			g1 += gap[mlen-P+t];
		}
		s6n = 2;
		s6Len[1] = mlen-P;
		if (g0 >= g1) {
			s6Start[0] = 0;
			s6Start[1] = P;
		} else {
			s6Start[0] = mlen-P;
			s6Start[1] = 0;
		}
	}
	s6Rest = 0.0;
	for (int j=0;j<mlen;j++) {
		int inBlock = 0;
		for (int q=0;q<s6n;q++) {
			if (j >= s6Start[q] && j < s6Start[q]+s6Len[q]) inBlock = 1;
		}
		if (inBlock) continue;
		double best = motif->m[j][0];
		for (int k=1;k<4;k++) if (motif->m[j][k] > best) best = motif->m[j][k];
		s6Rest += best;
	}

	// Tables. Forward, digit t of a k-mer is motif position a+t; reverse
	// complement, a+k-1-t. A block of 6 or 7: its first three digits' sum
	// plus its last k-3 digits' sum (the bound only; the 1e-6 slack covers
	// the rounding).
	for (int q=0;q<s6n;q++) {
		const int a = s6Start[q];
		const int k = s6Len[q];
		if (k >= 6) {
			const int loDigits = k-3;
			const int loSize = 1 << (2*loDigits);
			double* af = halfHi[0][q];
			double* ar = halfHi[1][q];
			double* bf = halfLo[0][q];
			double* br = halfLo[1][q];
			for (int x=0;x<64;x++) {
				const int d0 = x >> 4, d1 = (x >> 2) & 3, d2 = x & 3;
				af[x] = (fwd[a*5+d0] + fwd[(a+1)*5+d1]) + fwd[(a+2)*5+d2];
				ar[x] = (rev[(a+k-1)*5+d0] + rev[(a+k-2)*5+d1]) + rev[(a+k-3)*5+d2];
			}
			for (int y=0;y<loSize;y++) {
				double sf = 0.0;
				double sr = 0.0;
				for (int t=0;t<loDigits;t++) {
					const int d = (y >> (2*(loDigits-1-t))) & 3;
					if (t == 0) {
						sf = fwd[(a+3)*5+d];
						sr = rev[(a+k-4)*5+d];
					} else {
						sf = sf + fwd[(a+3+t)*5+d];
						sr = sr + rev[(a+k-4-t)*5+d];
					}
				}
				bf[y] = sf;
				br[y] = sr;
			}
			for (int st=0;st<2;st++) {
				if (t6[st][q] == NULL) t6[st][q] = new double[16385];
			}
			double* tf = t6[0][q];
			double* tr = t6[1][q];
			for (int hi=0;hi<64;hi++) {
				const double hf = af[hi];
				const double hr = ar[hi];
				for (int lo=0;lo<loSize;lo++) {
					tf[hi*loSize+lo] = hf + bf[lo];
					tr[hi*loSize+lo] = hr + br[lo];
				}
			}
			tf[64*loSize] = -INFINITY;
			tr[64*loSize] = -INFINITY;
			if (P == 7 && k == 6) {
				// code >> 2: 4096 + v for a prefix of exactly 6
				// bases, 8192 for fewer
				for (int x=0;x<4096;x++) {
					tf[4096+x] = tf[x];
					tr[4096+x] = tr[x];
				}
				tf[8192] = -INFINITY;
				tr[8192] = -INFINITY;
			}
		} else {
			double lf[5+25+125+625+3125];
			double lr[5+25+125+625+3125];
			double* pf = lf;
			double* pr = lr;
			for (int d=0;d<5;d++) {
				pf[d] = 0.0 + fwd[a*5+d];
				pr[d] = 0.0 + rev[(a+k-1)*5+d];
			}
			int width = 5;
			for (int t=1;t<k;t++) {
				double* nf = pf + width;
				double* nr = pr + width;
				for (int x=0;x<width;x++) {
					for (int d=0;d<5;d++) {
						nf[x*5+d] = pf[x] + fwd[(a+t)*5+d];
						nr[x*5+d] = pr[x] + rev[(a+k-1-t)*5+d];
					}
				}
				pf = nf;
				pr = nr;
				width *= 5;
			}
			int unused = 1;
			for (int t=k;t<5;t++) unused *= 5;
			for (int st=0;st<2;st++) {
				if (o5[st][q] == NULL) o5[st][q] = new double[3125];
			}
			o5Width[q] = width;
			o5Unused[q] = unused;
			for (int x=0;x<width;x++) {
				for (int r=0;r<unused;r++) {
					o5[0][q][x*unused+r] = pf[x];
					o5[1][q][x*unused+r] = pr[x];
				}
			}
			if (P == 7) {
				// the slots of code >> 4: below 3072 a prefix of 5 or
				// more bases (its first five, x & 1023); from 3072, y = x - 3072
				// = 4w + j holds a prefix of 4 (j = 0, value w), 3 (j = 1, w =
				// 4v), 2 (j = 2, w = 16v) or 1 (j = 3, w = 64v); else none
				for (int st=0;st<2;st++) {
					if (o5x[st][q] == NULL) o5x[st][q] = new double[4096];
				}
				for (int x=0;x<4096;x++) {
					int r = 0;
					int w = 0;
					if (x < 3072) {
						r = 5;
						w = x & 1023;
					} else {
						const int y = x - 3072;
						const int ww = y >> 2;
						const int j = y & 3;
						if (j == 0) {
							r = 4;
							w = ww;
						} else if (j == 1 && (ww & 3) == 0) {
							r = 3;
							w = ww >> 2;
						} else if (j == 2 && (ww & 15) == 0) {
							r = 2;
							w = ww >> 4;
						} else if (j == 3 && (ww & 63) == 0) {
							r = 1;
							w = ww >> 6;
						}
					}
					double vf = -INFINITY;
					double vr = -INFINITY;
					if (r >= k) {
						// the prefix's first k bases, as the base-5 index of
						// the level-k sums
						const int v = w >> (2*(r-k));
						int idx = 0;
						for (int t=0;t<k;t++) idx = idx*5 + ((v >> (2*(k-1-t))) & 3);
						vf = pf[idx];
						vr = pr[idx];
					}
					o5x[0][q][x] = vf;
					o5x[1][q][x] = vr;
				}
			}
		}
	}
	for (int q=0;q<s6n;q++) {
		for (int st=0;st<2;st++) {
			const int k = s6Len[q];
			if (k >= 6) {
				// the sum of two finite halves is never NaN
				const int loSize = 1 << (2*(k-3));
				int finite = 1;
				for (int x=0;x<64;x++) if (!std::isfinite(halfHi[st][q][x])) finite = 0;
				for (int y=0;y<loSize;y++) if (!std::isfinite(halfLo[st][q][y])) finite = 0;
				if (finite) continue;
				const double* T = t6[st][q];
				for (int x=0;x<64*loSize;x++) if (T[x] != T[x]) return 0;
			} else {
				for (int x=0;x<o5Width[q];x++) {
					const double v = o5[st][q][x*o5Unused[q]];
					if (v != v) return 0;
				}
			}
		}
	}
	// the blocks' gaps, to order the other blocks (the larger first: it
	// rejects more entries alone)
	double blockGap[3] = {0.0, 0.0, 0.0};
	for (int q=0;q<s6n;q++) {
		for (int t=0;t<s6Len[q];t++) blockGap[q] += gap[s6Start[q]+t];
	}

	// The pivot: of the blocks of P, the one whose passing P-mers occur
	// least often in the sequences. A window survives the bound only if its
	// pivot P-mer scores at least cut - s6Rest - (the other blocks' best);
	// the 1e-9 slack covers rounding. A table's maximum is maxHi + maxLo
	// (rounding is monotone), and a row hi whose halfHi + maxLo is below need
	// holds no entry at or above it.
	const double cut = motif->threshold - 1e-6;
	const int loSizeP = 1 << (2*(P-3));
	for (int st=0;st<2;st++) {
		if (s6pass[st] == NULL) {
			s6pass[st] = new int[16384];
			s6passScore[st] = new double[16384];
			s6passGap[st] = new double[16384];
		}
		double maxLo[3] = {0.0, 0.0, 0.0};
		int pivotBlocks = 0;
		for (int q=0;q<s6n;q++) {
			s6max[st][q] = -INFINITY;
			const int k = s6Len[q];
			if (k >= 6) {
				if (k == P) pivotBlocks++;
				const int loSize = 1 << (2*(k-3));
				double maxHi = -INFINITY;
				maxLo[q] = -INFINITY;
				for (int x=0;x<64;x++) if (halfHi[st][q][x] > maxHi) maxHi = halfHi[st][q][x];
				for (int y=0;y<loSize;y++) if (halfLo[st][q][y] > maxLo[q]) maxLo[q] = halfLo[st][q][y];
				s6max[st][q] = maxHi + maxLo[q];
			} else {
				for (int x=0;x<o5Width[q];x++) {
					if (o5[st][q][x*o5Unused[q]] > s6max[st][q]) s6max[st][q] = o5[st][q][x*o5Unused[q]];
				}
			}
		}
		int bestQ = -1;
		long long bestCount = 0;
		double bestNeed = 0.0;
		for (int q=0;q<s6n;q++) {
			if (s6Len[q] != P) continue;
			double others = 0.0;
			for (int r=0;r<s6n;r++) if (r != q) others += s6max[st][r];
			const double need = cut - s6Rest - others - 1e-9;
			long long total = 0;
			if (pivotBlocks > 1) {
				const double* T = t6[st][q];
				for (int hi=0;hi<64;hi++) {
					if (halfHi[st][q][hi] + maxLo[q] < need) continue;
					for (int x=hi*loSizeP;x<hi*loSizeP+loSizeP;x++) {
						if (!(T[x] < need)) total += (countP != NULL) ? countP[x] : 1;
					}
				}
			}
			if (bestQ < 0 || total < bestCount) {
				bestQ = q;
				bestCount = total;
				bestNeed = need;
			}
		}
		s6pivot[st] = bestQ;
		int no = 0;
		for (int r=0;r<s6n;r++) if (r != bestQ) s6oth[st][no++] = r;
		if (no == 2 && blockGap[s6oth[st][1]] > blockGap[s6oth[st][0]]) {
			const int tmp = s6oth[st][0];
			s6oth[st][0] = s6oth[st][1];
			s6oth[st][1] = tmp;
		}
		s6npass[st] = 0;
		// The flanks: the bases before the pivot k-mer in sequence order
		// (nearest first) are motif positions ap-1, ap-2 (forward) or,
		// complemented, ap+P, ap+P+1 (reverse); the three after are ap+P,
		// ap+P+1, ap+P+2 (forward) or ap-1, ap-2, ap-3 (reverse). A base b there
		// can add at most m[pos][b] - (the position's best) to the bound,
		// wherever the position falls (rest or another block); outside the
		// motif it adds nothing. An entry reads an N, or a base past the
		// sequence, as A: inside the motif either makes the window fail
		// anyway (-inf, or out of its sequence), so the reading only lets
		// through windows the exact score drops.
		{
			const int ap = s6Start[bestQ];
			int lpos[2];
			int rpos[3];
			for (int j=0;j<2;j++) lpos[j] = (st == 0) ? ap-1-j : ap+P+j;
			for (int j=0;j<3;j++) rpos[j] = (st == 0) ? ap+P+j : ap-1-j;
			const double* tab = (st == 0) ? fwd : rev;
			double dl[2][4];
			double dr[3][4];
			for (int j=0;j<5;j++) {
				const int pos = (j < 2) ? lpos[j] : rpos[j-2];
				double* dd = (j < 2) ? dl[j] : dr[j-2];
				if (pos < 0 || pos >= mlen) {
					for (int bb=0;bb<4;bb++) dd[bb] = 0.0;
					continue;
				}
				double best = tab[pos*5];
				for (int bb=1;bb<4;bb++) if (tab[pos*5+bb] > best) best = tab[pos*5+bb];
				for (int bb=0;bb<4;bb++) dd[bb] = tab[pos*5+bb] - best;
			}
			for (int x=0;x<16;x++) s6flankL[st][x] = dl[0][x & 3] + dl[1][x >> 2];
			for (int x=0;x<64;x++) s6flankR[st][x] = (dr[0][x & 3] + dr[1][(x >> 2) & 3]) + dr[2][x >> 4];
		}
		const double* T = t6[st][bestQ];
		for (int hi=0;hi<64;hi++) {
			if (halfHi[st][bestQ][hi] + maxLo[bestQ] < bestNeed) continue;
			for (int x=hi*loSizeP;x<hi*loSizeP+loSizeP;x++) {
				if (!(T[x] < bestNeed)) {
					s6pass[st][s6npass[st]] = x;
					s6passScore[st][s6npass[st]] = T[x];
					s6passGap[st][s6npass[st]] = T[x] - bestNeed;
					s6npass[st]++;
				}
			}
		}
	}
	six = 1;
	return 1;
}

// Base-4 6-mers of codes[0, length) into k6: 4096 where the 6-mer holds an N
// or runs past the end. Base-5 5-mers into k5, first base most significant,
// 4 past the end. A sequence shorter than 6 (4) gets all 4096 (3124).
static void kmerCodes6(const unsigned char* codes, int length, unsigned short* k6) {
	if (length < 6) {
		for (int i=0;i<length;i++) k6[i] = 4096;
		return;
	}
	int v = 0;
	int lastN = -1;
	for (int j=0;j<length;j++) {
		v = ((v << 2) | (codes[j] & 3)) & 4095;
		if (codes[j] == 4) lastN = j;
		if (j >= 5) k6[j-5] = (unsigned short)((lastN >= j-5) ? 4096 : v);
	}
	for (int i=length-5;i<length;i++) k6[i] = 4096;
}
// A 7-mer chunk's one code array: each position's prefix code. With
// r = the bases from the position on before an N or the sequence's end (at
// most 7) and v the base-4 value of those r bases:
//   r = 7: v (0-16383)    r = 6: 16384 + 4v    r = 5: 32768 + 16v
//   r = 4: 49152 + 64v    r = 3: 49152 + 256v + 16
//   r = 2: 49152 + 1024v + 32    r = 1: 49152 + 4096v + 48    r = 0: CK7_NONE
// A block of 7 reads the code (16384 for any code past 16383), a block of 6
// the code >> 2 (8192 past 8191), a block of 5 or fewer the code >> 4
// (0-4095); every r >= the block's length has its own slots, which hold the
// block's score of the prefix's first bases, and the rest are -inf.
#define CK7_NONE 49232
static void kmerCodes7(const unsigned char* codes, int length, unsigned short* k7) {
	// by r: the offset, and the bits of the 7-base value kept
	static const unsigned int base[8] = {CK7_NONE, 49152+48, 49152+32, 49152+16, 49152, 32768, 16384, 0};
	static const unsigned int keep[8] = {0, 0x3000, 0x3C00, 0x3F00, 0x3FC0, 0x3FF0, 0x3FFC, 0x3FFF};
	int run = 0;
	unsigned int v = 0;
	for (int p=length-1;p>=0;p--) {
		const unsigned int b = codes[p];
		run = (b == 4) ? 0 : run+1;
		// the bases at p..p+6, p's the most significant (past the end: 0)
		v = (v >> 2) | ((b & 3) << 12);
		const int r = (run < 7) ? run : 7;
		k7[p] = (unsigned short)(base[r] + (v & keep[r]));
	}
}
static void kmerCodes5(const unsigned char* codes, int length, unsigned short* k5) {
	if (length < 4) {
		for (int i=0;i<length;i++) k5[i] = 3124;
		return;
	}
	int v = 0;
	for (int t=0;t<5;t++) v = v*5 + ((t < length) ? codes[t] : 4);
	k5[0] = (unsigned short)v;
	for (int i=1;i<length;i++) {
		v = (v % 625)*5 + ((i+4 < length) ? codes[i+4] : 4);
		k5[i] = (unsigned short)v;
	}
}
void Sequence::prepareKmers(int withKmer5, unsigned char* codesAt) {
	if (codes != NULL && !codesShared) delete []codes;
	if (codesAt != NULL) {
		codes = codesAt;
		codesShared = 1;
	} else {
		codes = new unsigned char[length+1];
		codesShared = 0;
	}
	for (int i=0;i<length;i++) codes[i] = motifScanCode[(unsigned char)s[i]];
	if (kmer6 != NULL) delete []kmer6;
	kmer6 = NULL;
	if (kmer5 != NULL) delete []kmer5;
	kmer5 = NULL;
	if (length < 4 || !withKmer5) return;
	kmer5 = new unsigned short[length];
	kmerCodes5(codes, length, kmer5);
}
KmerChunk::KmerChunk() {
	first = 0;
	count = 0;
	start = NULL;
	pos = NULL;
	sid = NULL;
	cstart = NULL;
	clen = NULL;
	ck7 = NULL;
	ck6 = NULL;
	ck5 = NULL;
	blockSeq = NULL;
	flanks = 0;
	codesBlock = NULL;
}
KmerChunk::~KmerChunk() {
	// a 6/7-mer chunk's pos and code arrays come from chunkArray
	const int fromChunkArray = (ck7 != NULL || ck6 != NULL);
	delete []start;
	if (fromChunkArray) {
		free(pos);
	} else {
		delete []pos;
	}
	delete []sid;
	delete []cstart;
	delete []clen;
	free(ck7);
	free(ck6);
	free(ck5);
	delete []blockSeq;
	free(codesBlock);
}
// The index holds a position when the first four bases of its 5-mer are
// not N (the fifth may be N, or past the end).
static int kmerHasN(int v) {
	return (v/5)%5 == 4 || (v/25)%5 == 4 || (v/125)%5 == 4 || v/625 == 4;
}
// Runs fn(i) for every i in [0, n) on up to `threads` threads, which take
// items in order from a shared counter.
template <typename F>
static void parallelFor(int n, int threads, F fn) {
	if (threads > n) threads = n;
	if (threads <= 1) {
		for (int i=0;i<n;i++) fn(i);
		return;
	}
	std::atomic<int> next(0);
	std::vector<std::thread> pool;
	for (int t=0;t<threads;t++) {
		pool.emplace_back([&]() {
			for (int i=next++;i<n;i=next++) fn(i);
		});
	}
	for (size_t t=0;t<pool.size();t++) pool[t].join();
}

static void buildKmerChunk(Sequence** seqs, KmerChunk* c) {
	const int first = c->first;
	const int last = c->first + c->count;
	c->start = new int[3126];
	for (int v=0;v<3126;v++) c->start[v] = 0;
	for (int q=first;q<last;q++) {
		const unsigned short* km = seqs[q]->kmer5;
		if (km == NULL) continue;
		for (int p=0;p+4<=seqs[q]->length;p++) {
			if (!kmerHasN(km[p])) c->start[km[p]+1]++;
		}
	}
	for (int v=0;v<3125;v++) c->start[v+1] += c->start[v];
	int total = c->start[3125];
	c->pos = new int[total > 0 ? total : 1];
	c->sid = new unsigned short[total > 0 ? total : 1];
	int* cursor = new int[3125];
	for (int v=0;v<3125;v++) cursor[v] = c->start[v];
	for (int q=first;q<last;q++) {
		const unsigned short* km = seqs[q]->kmer5;
		if (km == NULL) continue;
		for (int p=0;p+4<=seqs[q]->length;p++) {
			if (kmerHasN(km[p])) continue;
			int e = cursor[km[p]]++;
			c->pos[e] = p;
			c->sid[e] = (unsigned short)(q-first);
		}
	}
	delete []cursor;
}

// A chunk of the P-mer index (P = 6 or 7): its sequences' 7-mer (P = 7),
// 6-mer and 5-mer codes end to end, and the positions of each P-mer.
// n items for a 6/7-mer chunk's arrays, released with free(): 2 MB or
// more go in 2 MB-aligned memory advised MADV_HUGEPAGE, so the index visit's
// random reads into them miss the TLB less. Also the chunks' codes
// blocks and the sequences' text blocks.
template <typename T>
static T* chunkArray(long n) {
	size_t bytes = (size_t)(n > 0 ? n : 1) * sizeof(T);
	const size_t huge = (size_t)2 << 20;
	void* p = NULL;
	if (posix_memalign(&p, bytes >= huge ? huge : 64, bytes) != 0) {
		fprintf(stderr, "!!! out of memory for the k-mer index !!!\n");
		exit(1);
	}
#ifdef MADV_HUGEPAGE
	if (bytes >= huge) (void)madvise(p, bytes, MADV_HUGEPAGE);
#endif
	return (T*)p;
}
// An allocator whose blocks come from chunkArray (huge pages from
// 2 MB), for the FASTA load's piece texts.
template <typename T>
struct HugeAllocator {
	typedef T value_type;
	HugeAllocator() {}
	template <typename U> HugeAllocator(const HugeAllocator<U>&) {}
	T* allocate(size_t n) { return chunkArray<T>((long)n); }
	void deallocate(T* p, size_t) { free(p); }
	template <typename U> bool operator==(const HugeAllocator<U>&) const { return true; }
	template <typename U> bool operator!=(const HugeAllocator<U>&) const { return false; }
};
static void buildKmerChunkP(Sequence** seqs, KmerChunk* c, int pad, int P) {
	const int first = c->first;
	const int last = c->first + c->count;
	const int buckets = 1 << (2*P);
	c->cstart = new int[c->count > 0 ? c->count : 1];
	c->clen = new int[c->count > 0 ? c->count : 1];
	long span = pad;
	for (int q=first;q<last;q++) {
		c->cstart[q-first] = (int)span;
		c->clen[q-first] = seqs[q]->length;
		span += seqs[q]->length;
	}
	span += pad;
	// a 7-mer chunk has ck7 alone, in the prefix code
	if (P == 7) {
		c->ck7 = chunkArray<unsigned short>(span);
		for (long j=0;j<pad;j++) {
			c->ck7[j] = CK7_NONE;
			c->ck7[span-1-j] = CK7_NONE;
		}
	} else {
		c->ck6 = chunkArray<unsigned short>(span);
		c->ck5 = chunkArray<unsigned short>(span);
		for (long j=0;j<pad;j++) {
			c->ck6[j] = 4096;
			c->ck5[j] = 3124;
			c->ck6[span-1-j] = 4096;
			c->ck5[span-1-j] = 3124;
		}
	}
	const unsigned short* ckP = (P == 7) ? c->ck7 : c->ck6;
	c->start = new int[buckets+1];
	for (int v=0;v<=buckets;v++) c->start[v] = 0;
	for (int q=first;q<last;q++) {
		const int at = c->cstart[q-first];
		if (P == 7) {
			kmerCodes7(seqs[q]->codes, seqs[q]->length, c->ck7+at);
		} else {
			kmerCodes6(seqs[q]->codes, seqs[q]->length, c->ck6+at);
			kmerCodes5(seqs[q]->codes, seqs[q]->length, c->ck5+at);
		}
		const unsigned short* km = ckP+at;
		for (int p=0;p<seqs[q]->length;p++) {
			if (km[p] < buckets) c->start[km[p]+1]++;
		}
	}
	for (int v=0;v<buckets;v++) c->start[v+1] += c->start[v];
	int total = c->start[buckets];
	c->pos = chunkArray<int>(total);
	int* cursor = new int[buckets];
	for (int v=0;v<buckets;v++) cursor[v] = c->start[v];
	c->flanks = (span < (1L << 22));
	for (int q=first;q<last;q++) {
		const int at = c->cstart[q-first];
		const unsigned short* km = ckP+at;
		const unsigned char* code = seqs[q]->codes;
		const int len = seqs[q]->length;
		for (int p=0;p<len;p++) {
			if (km[p] >= buckets) continue;
			unsigned int v = (unsigned int)(at+p);
			if (c->flanks) {
				// 2 bits a base: an N (4) or past the sequence reads as A
				const unsigned int l1 = (p > 0) ? (code[p-1] & 3) : 0;
				const unsigned int l2 = (p > 1) ? (code[p-2] & 3) : 0;
				const unsigned int r1 = (p+P < len) ? (code[p+P] & 3) : 0;
				const unsigned int r2 = (p+P+1 < len) ? (code[p+P+1] & 3) : 0;
				const unsigned int r3 = (p+P+2 < len) ? (code[p+P+2] & 3) : 0;
				v |= ((l1 | (l2 << 2)) << 22) | ((r1 | (r2 << 2) | (r3 << 4)) << 26);
			}
			c->pos[cursor[km[p]]++] = (int)v;
		}
	}
	delete []cursor;
	const long blocks = (span + 63) / 64;
	c->blockSeq = new unsigned short[blocks > 0 ? blocks : 1];
	int qs = 0;
	for (long b=0;b<blocks;b++) {
		while (qs+1 < c->count && c->cstart[qs+1] <= b*64) qs++;
		c->blockSeq[b] = (unsigned short)qs;
	}
}

// mask: each thread takes a fixed run of sequences and masks it with every
// motif in file order, so each sequence sees the motifs in the order -p 1
// gives it, whatever -p is. (Threads over motifs would mask one sequence
// from several threads at once.) Sequences with a site are counted
// afterwards in sequence order, as the -p 1 loop counts them.
void SequenceArray::maskMotifs(Motif** motifs, int numMotifs, int strand, int offset, int threads) {
	for (int k=0;k<numMotifs;k++) motifs[k]->logxform();
	if (threads > numSeqs) threads = numSeqs;
	if (threads < 1) threads = 1;
	long total = 0;
	for (int i=0;i<numSeqs;i++) total += seqs[i]->length;
	std::vector<int> bounds(threads+1, numSeqs);
	bounds[0] = 0;
	long acc = 0;
	int t = 1;
	for (int i=0;i<numSeqs && t<threads;i++) {
		acc += seqs[i]->length;
		if (acc*threads >= total*t) bounds[t++] = i+1;
	}
	std::vector<std::vector<unsigned long long> > hit(threads);
	parallelFor(threads, threads, [&](int r) {
		const int lo = bounds[r];
		const int hi = bounds[r+1];
		const int words = (hi-lo+63)/64;
		hit[r].assign((size_t)numMotifs*words, 0);
		int* sitePositions = new int[longestSequence];
		char* siteStrands = new char[longestSequence];
		pthread_mutex_t notused;
		pthread_mutex_init(&notused,NULL);
		MotifScanTables* tables = new MotifScanTables();
		for (int k=0;k<numMotifs;k++) {
			tables->build(motifs[k]);
			for (int i=lo;i<hi;i++) {
				double score = 0;
				int nsites = seqs[i]->motifScan(motifs[k],strand,offset,KNOWN_MASK_MOTIFS,score,
							sitePositions,siteStrands,NULL,notused,tables);
				if (nsites > 0) hit[r][(size_t)k*words + (i-lo)/64] |= 1ULL << ((i-lo)%64);
			}
		}
		delete tables;
		pthread_mutex_destroy(&notused);
		delete []sitePositions;
		delete []siteStrands;
	});
	for (int k=0;k<numMotifs;k++) {
		Motif* m = motifs[k];
		m->numTarget = 0.0;
		m->numBackground = 0.0;
		for (int r=0;r<threads;r++) {
			const int words = (bounds[r+1]-bounds[r]+63)/64;
			for (int i=bounds[r];i<bounds[r+1];i++) {
				int j = i-bounds[r];
				if (!((hit[r][(size_t)k*words + j/64] >> (j%64)) & 1ULL)) continue;
				if (seqs[i]->g == 1) {
					m->numTarget += seqs[i]->w;
				} else {
					m->numBackground += seqs[i]->w;
				}
			}
		}
		if (numTarget > 0.0) m->freqTarget = m->numTarget / numTarget;
		if (numBackground > 0.0) m->freqBackground = m->numBackground / numBackground;
	}
}

void SequenceArray::prepareKmers(int threads, int need5, int need6, int maxWidth, int need7) {
	delete []chunks;
	chunks = new KmerChunk[numSeqs > 0 ? numSeqs : 1];
	numChunks = 0;
	// About one chunk per thread, at least 16,384 bases and at most 2 Mb: a
	// chunk's index has 4096 6-mer or 16,384 7-mer buckets, and each passing
	// k-mer costs a bucket lookup per chunk, so few large chunks visit their
	// entries with less overhead; past a few Mb its code arrays (2 bytes a
	// base each) leave the cache that the index visit reads them from at
	// random (1 Mb was best with the 6-mer index, 2 Mb with the 7-mer one).
	long totalBases = 0;
	for (int q=0;q<numSeqs;q++) totalBases += (seqs[q]->length > 3) ? seqs[q]->length-3 : 0;
	const int units = (threads > 0) ? threads : 1;
	long chunkBases = (totalBases + units - 1) / units;
	if (chunkBases < 16384) chunkBases = 16384;
	if (chunkBases > 2097152) chunkBases = 2097152;
	int i = 0;
	while (i < numSeqs) {
		int first = i;
		long bases = 0;
		while (i < numSeqs && i-first < 65535) {
			long add = seqs[i]->length - 3;
			if (add < 0) add = 0;
			if (i > first && bases + add > chunkBases) break;
			bases += add;
			i++;
		}
		KmerChunk* c = &(chunks[numChunks++]);
		c->first = first;
		c->count = i-first;
	}
	index5 = need5;
	delete []chunks6;
	chunks6 = NULL;
	delete []count6;
	count6 = NULL;
	delete []chunks7;
	chunks7 = NULL;
	delete []count7;
	count7 = NULL;
	if (need6) {
		chunks6 = new KmerChunk[numChunks > 0 ? numChunks : 1];
		for (int k=0;k<numChunks;k++) {
			chunks6[k].first = chunks[k].first;
			chunks6[k].count = chunks[k].count;
		}
	}
	if (need7) {
		chunks7 = new KmerChunk[numChunks > 0 ? numChunks : 1];
		for (int k=0;k<numChunks;k++) {
			chunks7[k].first = chunks[k].first;
			chunks7[k].count = chunks[k].count;
		}
	}
	// a chunk's sequences, then its indexes, on one thread
	parallelFor(numChunks, threads, [&](int k) {
		// the chunk's sequences' codes, in one block it owns
		long codeBytes = 0;
		for (int q=chunks[k].first;q<chunks[k].first+chunks[k].count;q++) codeBytes += seqs[q]->length + 1;
		chunks[k].codesBlock = chunkArray<unsigned char>(codeBytes);
		long codeAt = 0;
		for (int q=chunks[k].first;q<chunks[k].first+chunks[k].count;q++) {
			seqs[q]->prepareKmers(need5, chunks[k].codesBlock + codeAt);
			codeAt += seqs[q]->length + 1;
		}
		if (need5) buildKmerChunk(seqs, &(chunks[k]));
		if (need6) buildKmerChunkP(seqs, &(chunks6[k]), maxWidth+8, 6);
		if (need7) buildKmerChunkP(seqs, &(chunks7[k]), maxWidth+8, 7);
	});
	// maxChunk6Entries: the most entries in one chunk of either index
	maxChunk6Entries = 0;
	if (need6) {
		count6 = new int[4096];
		for (int x=0;x<4096;x++) count6[x] = 0;
		for (int k=0;k<numChunks;k++) {
			for (int x=0;x<4096;x++) count6[x] += chunks6[k].start[x+1] - chunks6[k].start[x];
			if (chunks6[k].start[4096] > maxChunk6Entries) maxChunk6Entries = chunks6[k].start[4096];
		}
	}
	if (need7) {
		count7 = new int[16384];
		for (int x=0;x<16384;x++) count7[x] = 0;
		for (int k=0;k<numChunks;k++) {
			for (int x=0;x<16384;x++) count7[x] += chunks7[k].start[x+1] - chunks7[k].start[x];
			if (chunks7[k].start[16384] > maxChunk6Entries) maxChunk6Entries = chunks7[k].start[16384];
		}
	}
}

// find with the index: sum the given candidate windows exactly, in the order
// Sequence::motifScan's find path visits them (forward ascending, reverse
// descending; both lists arrive ascending), with the same per-window logic.
int Sequence::scanCandidates(Motif* motif, int strand, int offset, double &bestScore,
					int* sitePositions, char* siteStrands, FILE* fp,
					const MotifScanTables* tables, const int* fwdCand, int nf,
					const int* revCand, int nr) {
	if (offset == MOTIF_OFFSET_NULL) {
		offset = length/-2;
	}
	const unsigned char* code = codes;
	const int mlen = tables->length;
	const double* fwd = tables->fwd;
	const double* rev = tables->rev;
	bestScore = -1e10;
	int numSites = 0;
	if (strand == MOTIF_STRAND_POS || strand==MOTIF_STRAND_BOTH) {
		for (int ci=0;ci<nf;ci++) {
			const int i = fwdCand[ci];
			const unsigned char* c = code+i;
			double score = 0.0;
			for (int j=0;j<mlen;j++) {
				score += fwd[j*5+c[j]];
			}
			if (score == -INFINITY) continue;
			if (score > bestScore) bestScore = score;
			if (score >= motif->threshold) {
				sitePositions[numSites] = i+offset;
				siteStrands[numSites] = MOTIF_STRAND_POS;
				numSites++;
				scanOutSite(fp,name,i+offset,&(s[i]),motif->length,motif->name,'+',score);
			}
		}
	}
	if (strand == MOTIF_STRAND_BOTH || strand==MOTIF_STRAND_NEG) {
		for (int ci=nr-1;ci>=0;ci--) {
			const int i = revCand[ci];
			const unsigned char* c = code+i;
			double score = 0.0;
			for (int j=0;j<mlen;j++) {
				score += rev[j*5+c[-j]];
			}
			if (score == -INFINITY) continue;
			if (score > bestScore) bestScore = score;
			if (score >= motif->threshold) {
				sitePositions[numSites] = i+offset;
				siteStrands[numSites] = MOTIF_STRAND_NEG;
				numSites++;
				scanOutSite(fp,name,i+offset,&(s[i-motif->length+1]),motif->length,motif->name,'-',score);
			}
		}
	}
	return numSites;
}

// find with the index: the windows of chunk c that survive find's bound,
// grouped by sequence. Visits only positions whose pivot-block 4-mer is in
// the motif's pass list, then applies the bound expression of
// Sequence::motifScan's find path unchanged. Fills, per sequence q of the
// chunk, candidates [off[q], off[q+1]) of cand, ascending.
static void orderCandidates(const long long* keys, int n, int count, int*& cand, int& candSize, int* off);
static void collectCandidates(Sequence** seqs, const KmerChunk* c, const MotifScanTables* t,
		int st, double cut, long long*& keys, int& keySize,
		int*& cand, int& candSize, int* off) {
	const int mlen = t->length;
	const int a0 = t->blockStart[0];
	const int a1 = t->blockStart[1];
	const int len0 = t->blockLen[0];
	const int len1 = t->blockLen[1];
	const double blockRest = t->blockRest;
	const double* T0 = (st == 0) ? t->bfwd : t->brev;
	const double* T1 = T0 + 3125;
	const int ap = t->blockStart[t->pivot[st]];
	const int lenp = t->blockLen[t->pivot[st]];
	int total = 0;
	for (int k=0;k<t->npass[st];k++) {
		int x = t->pass[st][k];
		total += c->start[x+1] - c->start[x];
	}
	if (keySize < total+1) {
		delete []keys;
		keySize = total+1;
		keys = new long long[keySize];
	}
	int n = 0;
	for (int k=0;k<t->npass[st];k++) {
		int x = t->pass[st][k];
		for (int e=c->start[x];e<c->start[x+1];e++) {
			int q = c->sid[e];
			const Sequence* sq = seqs[c->first+q];
			const unsigned short* km = sq->kmer5;
			int i;
			double bound;
			if (st == 0) {
				i = c->pos[e] - ap;
				if (i < 0 || i > sq->length-mlen) continue;
				bound = blockRest + T0[km[i+a0]] + T1[km[i+a1]];
			} else {
				i = c->pos[e] + ap + lenp - 1;
				if (i < mlen-1 || i > sq->length-1) continue;
				bound = blockRest + T0[km[i-a0-len0+1]] + T1[km[i-a1-len1+1]];
			}
			if (bound < cut) continue;
			keys[n++] = ((long long)q << 32) | (long long)i;
		}
	}
	orderCandidates(keys, n, c->count, cand, candSize, off);
}

// Candidates keys[0, n) as (sequence << 32 | position) into cand, grouped by
// sequence q of the chunk's count, [off[q], off[q+1]), positions ascending.
static void orderCandidates(const long long* keys, int n, int count, int*& cand, int& candSize, int* off) {
	if (candSize < n+1) {
		delete []cand;
		candSize = n+1;
		cand = new int[candSize];
	}
	// Order by (sequence, position): a counting sort by sequence, then an
	// insertion sort of each sequence's positions when they are few.
	for (int q=0;q<=count;q++) off[q] = 0;
	for (int k=0;k<n;k++) off[(int)(keys[k] >> 32)+1]++;
	for (int q=0;q<count;q++) off[q+1] += off[q];
	static thread_local std::vector<int> fillv;
	if ((int)fillv.size() < count+1) fillv.resize(count+1);
	int* fill = fillv.data();
	for (int q=0;q<count;q++) fill[q] = off[q];
	for (int k=0;k<n;k++) {
		cand[fill[(int)(keys[k] >> 32)]++] = (int)(keys[k] & 0xffffffffLL);
	}
	for (int q=0;q<count;q++) {
		if (off[q+1]-off[q] > 32) {
			// a long sequence: its candidates are many
			std::sort(cand+off[q], cand+off[q+1]);
			continue;
		}
		for (int a=off[q]+1;a<off[q+1];a++) {
			int x = cand[a];
			int b = a-1;
			while (b >= off[q] && cand[b] > x) {
				cand[b+1] = cand[b];
				b--;
			}
			cand[b+1] = x;
		}
	}
}

// find with the 6/7-mer index: as collectCandidates, for buildSix's blocks.
// Every entry of a passing pivot P-mer shares that P-mer's partial score, so
// an entry costs a table lookup per other block (the second only when the
// first, with the second's best, still reaches the cut). Forward, the window
// starting at i holds the pivot P-mer at i+ap and another block's k-mer at
// i+ao; reverse, the window ending at i holds them at i-ap-P+1 and
// i-ao-k+1. Positions are in the chunk's end-to-end code arrays, whose
// padding covers any motif width; an entry whose window leaves its sequence
// reads a neighbour's or the padding's codes, and is dropped by the range
// check once it passes the bound. O1, O2: the other blocks' kinds (0: none,
// 5: up to 5 positions through ck5, 6 and 7 through ck6 and ck7).
// The table slot of a block of kind O at code x: with the 7-mer
// pivot every block reads ck7's prefix code.
template <int P, int O>
static inline unsigned int blockSlot(unsigned int x) {
	if (P != 7) return x;
	if (O == 7) return (x < 16384) ? x : 16384;
	if (O == 6) return ((x >> 2) < 8192) ? (x >> 2) : 8192;
	return x >> 4;
}
template <int P, int O1, int O2, int FL>
static int gatherSix(Sequence** seqs, const KmerChunk* c, const MotifScanTables* t,
		int st, double cut, long long* keys, int atMode) {
	const int mlen = t->length;
	const int p = t->s6pivot[st];
	const int ap = t->s6Start[p];
	const int b1 = t->s6oth[st][0];
	const int b2 = t->s6oth[st][1];
	const int a1 = (O1 != 0) ? t->s6Start[b1] : 0;
	const int l1 = (O1 != 0) ? t->s6Len[b1] : 0;
	const int a2 = (O2 != 0) ? t->s6Start[b2] : 0;
	const int l2 = (O2 != 0) ? t->s6Len[b2] : 0;
	const double* T1 = (O1 >= 6) ? t->t6[st][b1] : ((O1 == 5) ? ((P == 7) ? t->o5x[st][b1] : t->o5[st][b1]) : NULL);
	const double* T2 = (O2 >= 6) ? t->t6[st][b2] : ((O2 == 5) ? ((P == 7) ? t->o5x[st][b2] : t->o5[st][b2]) : NULL);
	const double max2 = (O2 != 0) ? t->s6max[st][b2] : 0.0;
	const int pivotShift = (st == 0) ? -ap : ap+P-1;
	const int shift1 = (st == 0) ? a1 : -a1-l1+1;
	const int shift2 = (st == 0) ? a2 : -a2-l2+1;
	const int lowI = (st == 0) ? 0 : mlen-1;
	const int highLess = (st == 0) ? mlen : 1;
	const int np = t->s6npass[st];
	const int* px = t->s6pass[st];
	const double* ps = t->s6passScore[st];
	const unsigned short* ck1 = (P == 7 || O1 == 7) ? c->ck7 : ((O1 == 6) ? c->ck6 : c->ck5);
	const unsigned short* ck2 = (P == 7 || O2 == 7) ? c->ck7 : ((O2 == 6) ? c->ck6 : c->ck5);
	int n = 0;
	if (FL) {
		// Entries carry their flanking bases. In batches of up to FLAT
		// entries: phase 1, branch-free, writes the entries whose flanks pass
		// the flank test into a list on the stack (with their pass-list
		// index); phase 2 runs the lookups over the list, prefetched a fixed
		// distance ahead. A bucket may span batches.
		enum { FLAT = 4096 };
		int fp[FLAT];
		unsigned short fk[FLAT];
		const unsigned int* upos = (const unsigned int*)c->pos;
		const double* gap = t->s6passGap[st];
		const double* fL = t->s6flankL[st];
		const double* fR = t->s6flankR[st];
		int k = 0;
		int e = (np > 0) ? c->start[px[0]] : 0;
		while (k < np) {
			int m = 0;
			while (k < np && m < FLAT) {
				const int eEnd = c->start[px[k]+1];
				if (O1 == 0 && t->s6Rest + ps[k] < cut) {
					e = eEnd;
				} else {
					// a later bucket's head (bucket starts are scattered)
					if (k+4 < np) __builtin_prefetch(upos + c->start[px[k+4]]);
					const double g = gap[k];
					// m grows by at most one an entry, so FLAT - m entries fit
					int lim = e + (FLAT - m);
					if (lim > eEnd) lim = eEnd;
					for (;e<lim;e++) {
						const unsigned int u = upos[e];
						fp[m] = (int)(u & 0x3FFFFFu);
						fk[m] = (unsigned short)k;
						m += (g + fL[(u >> 22) & 15] + fR[u >> 26] >= 0.0);
					}
				}
				if (e == eEnd) {
					k++;
					if (k < np) e = c->start[px[k]];
				}
			}
			for (int j=0;j<m;j++) {
				const int ep = fp[j];
				const int at = ep + pivotShift;
				if (O1 != 0) {
					if (j + 32 < m) {
						__builtin_prefetch(ck1 + fp[j+32] + pivotShift + shift1);
						// the second block's code, and the blockSeq slot a
						// survivor reads
						if (O2 != 0) __builtin_prefetch(ck2 + fp[j+32] + pivotShift + shift2);
						if (!atMode) __builtin_prefetch(c->blockSeq + (fp[j+32] >> 6));
					}
					const double s1 = t->s6Rest + ps[fk[j]] + T1[blockSlot<P,O1>(ck1[at+shift1])];
					if (O2 != 0) {
						if (s1 + max2 < cut) continue;
						if (s1 + T2[blockSlot<P,O2>(ck2[at+shift2])] < cut) continue;
					} else if (s1 < cut) {
						continue;
					}
				}
				if (atMode) {
					// the window's position in the chunk; resolveAts
					// finds its sequence
					((unsigned int*)keys)[n++] = (unsigned int)at;
					continue;
				}
				int q = c->blockSeq[ep >> 6];
				while (q+1 < c->count && c->cstart[q+1] <= ep) q++;
				const int i = at - c->cstart[q];
				if (i < lowI || i > c->clen[q] - highLess) continue;
				keys[n++] = ((long long)q << 32) | (long long)i;
			}
		}
		return n;
	}
	for (int k=0;k<np;k++) {
		const double base = t->s6Rest + ps[k];
		if (O1 == 0 && base < cut) continue;
		const int x = px[k];
		if (O1 != 0 && k+1 < np) {
			// the next bucket's first entries: the in-bucket prefetch below
			// reaches only entries 32 past the current one
			const int nx = px[k+1];
			const int nEnd = (c->start[nx+1] < c->start[nx]+32) ? c->start[nx+1] : c->start[nx]+32;
			for (int j=c->start[nx];j<nEnd;j++) __builtin_prefetch(ck1 + c->pos[j] + pivotShift + shift1);
		}
		const int eEnd = c->start[x+1];
		for (int e=c->start[x];e<eEnd;e++) {
			const int at = c->pos[e] + pivotShift;
			if (O1 != 0) {
				// a bucket's entries step through the chunk's codes ~4 kB apart
				if (e + 32 < eEnd) __builtin_prefetch(ck1 + c->pos[e+32] + pivotShift + shift1);
				const double s1 = base + T1[blockSlot<P,O1>(ck1[at+shift1])];
				if (O2 != 0) {
					if (s1 + max2 < cut) continue;
					if (s1 + T2[blockSlot<P,O2>(ck2[at+shift2])] < cut) continue;
				} else if (s1 < cut) {
					continue;
				}
			}
			// the entry's sequence: the last one starting at or before its
			// position (zero-length sequences are stepped over)
			const int ep = c->pos[e];
			if (atMode) {
				((unsigned int*)keys)[n++] = (unsigned int)at;
				continue;
			}
			int q = c->blockSeq[ep >> 6];
			while (q+1 < c->count && c->cstart[q+1] <= ep) q++;
			const int i = at - c->cstart[q];
			if (i < lowI || i > c->clen[q] - highLess) continue;
			keys[n++] = ((long long)q << 32) | (long long)i;
		}
	}
	return n;
}
static int blockKind(int len) {
	if (len >= 7) return 7;
	if (len == 6) return 6;
	return 5;
}
// The keys gatherSix leaves for chunk c and strand st, unordered; keys holds
// at least the chunk's index entries. With atMode, keys holds instead
// 32-bit window positions in the chunk (start forward, end reverse), not
// yet checked against their sequence's ends (resolveAts does).
static int sixKeys(Sequence** seqs, const KmerChunk* c, const MotifScanTables* t,
		int st, double cut, long long* keys, int atMode) {
	const int k1 = (t->s6n >= 2) ? blockKind(t->s6Len[t->s6oth[st][0]]) : 0;
	const int k2 = (t->s6n >= 3) ? blockKind(t->s6Len[t->s6oth[st][1]]) : 0;
#define SIX_GATHER(PP, A, B) if (k1 == A && k2 == B) return c->flanks \
		? gatherSix<PP,A,B,1>(seqs, c, t, st, cut, keys, atMode) : gatherSix<PP,A,B,0>(seqs, c, t, st, cut, keys, atMode);
	if (t->sixP == 7) {
		SIX_GATHER(7,0,0)
		SIX_GATHER(7,5,0) SIX_GATHER(7,6,0) SIX_GATHER(7,7,0)
		SIX_GATHER(7,5,5) SIX_GATHER(7,5,6) SIX_GATHER(7,5,7)
		SIX_GATHER(7,6,5) SIX_GATHER(7,6,6) SIX_GATHER(7,6,7)
		SIX_GATHER(7,7,5) SIX_GATHER(7,7,6) SIX_GATHER(7,7,7)
	} else {
		SIX_GATHER(6,0,0)
		SIX_GATHER(6,5,0) SIX_GATHER(6,6,0) SIX_GATHER(6,7,0)
		SIX_GATHER(6,5,5) SIX_GATHER(6,5,6) SIX_GATHER(6,6,5) SIX_GATHER(6,6,6)
	}
#undef SIX_GATHER
	fprintf(stderr, "!!! no index path for blocks %d %d !!!\n", k1, k2);
	exit(1);
	return 0;
}
static int sixKeys(Sequence** seqs, const KmerChunk* c, const MotifScanTables* t,
		int st, double cut, long long* keys, int atMode);
// Window positions from sixKeys' atMode (n of them, in ats), sorted,
// then each turned into its chunk sequence (seq) and window position (pos),
// keeping only windows inside their sequence: the sequence is the last one
// starting at or before the position, found walking forward. Returns how
// many were kept.
// ats[0, n) ascending by an LSD radix sort, 8 bits a pass.
static void radixSortU32(unsigned int* ats, int n) {
	static thread_local std::vector<unsigned int> tmpv;
	if ((int)tmpv.size() < n) tmpv.resize(n);
	unsigned int* tmp = tmpv.data();
	unsigned int hi = 0;
	for (int k=0;k<n;k++) hi |= ats[k];
	const int passes = (hi < (1u << 24)) ? 3 : 4;
	unsigned int* src = ats;
	unsigned int* dst = tmp;
	for (int pass=0;pass<passes;pass++) {
		const int shift = 8*pass;
		int count[257];
		for (int b=0;b<257;b++) count[b] = 0;
		for (int k=0;k<n;k++) count[((src[k] >> shift) & 255)+1]++;
		for (int b=0;b<256;b++) count[b+1] += count[b];
		for (int k=0;k<n;k++) dst[count[(src[k] >> shift) & 255]++] = src[k];
		unsigned int* t = src;
		src = dst;
		dst = t;
	}
	if (src != ats) memcpy(ats, src, (size_t)n*sizeof(unsigned int));
}
static int resolveAts(const KmerChunk* c, int st, int mlen, unsigned int* ats, int n,
		int* seq, int* pos) {
	if (n > 64) {
		radixSortU32(ats, n);
	} else {
		for (int a=1;a<n;a++) {
			unsigned int x = ats[a];
			int b = a-1;
			while (b >= 0 && ats[b] > x) {
				ats[b+1] = ats[b];
				b--;
			}
			ats[b+1] = x;
		}
	}
	const int lowI = (st == 0) ? 0 : mlen-1;
	const int highLess = (st == 0) ? mlen : 1;
	int q = 0;
	int m = 0;
	for (int k=0;k<n;k++) {
		const int at = (int)ats[k];
		const int b = c->blockSeq[at >> 6];
		if (b > q) q = b;
		while (q+1 < c->count && c->cstart[q+1] <= at) q++;
		const int i = at - c->cstart[q];
		if (i < lowI || i > c->clen[q] - highLess) continue;
		seq[m] = q;
		pos[m] = i;
		m++;
	}
	return m;
}
static void collectCandidates6(Sequence** seqs, const KmerChunk* c, const MotifScanTables* t,
		int st, double cut, long long*& keys, int& keySize,
		int*& cand, int& candSize, int* off) {
	int total = 0;
	for (int k=0;k<t->s6npass[st];k++) {
		int x = t->s6pass[st][k];
		total += c->start[x+1] - c->start[x];
	}
	if (keySize < total+1) {
		delete []keys;
		keySize = total+1;
		keys = new long long[keySize];
	}
	int n = sixKeys(seqs, c, t, st, cut, keys, 0);
	orderCandidates(keys, n, c->count, cand, candSize, off);
}

int Sequence::motifScan(Motif* motif, int strand,int offset, int mode,double &bestScore,
										int* sitePositions, char* siteStrands,
										FILE* fp, pthread_mutex_t fileLock,
										const MotifScanTables* tables) {

	if (offset == MOTIF_OFFSET_NULL) {
		offset = length/-2;
	}
	// The sequence as table columns: prepared once for find, which never
	// edits them; otherwise a per-thread copy (mask writes 4s into it).
	static thread_local unsigned char* scratch = NULL;
	static thread_local int scratchSize = 0;
	unsigned char* code = codes;
	if (code == NULL || mode != KNOWN_FIND_ALL || fp == NULL) {
		if (scratchSize < length+5) {
			delete []scratch;
			scratchSize = length+5;
			scratch = new unsigned char[scratchSize];
		}
		for (int i=0;i<length;i++) scratch[i] = motifScanCode[(unsigned char)s[i]];
		// past the end, for the 5-mer codes below
		for (int i=length;i<length+4;i++) scratch[i] = 4;
		code = scratch;
	}

	if (tables == NULL) {
		static thread_local MotifScanTables* own = NULL;
		if (own == NULL) own = new MotifScanTables();
		own->build(motif);
		tables = own;
	}
	const int mlen = tables->length;
	const int K = tables->K;
	const double restOfBest = tables->restOfBest;
	const double* fwd = tables->fwd;
	const double* rev = tables->rev;
	const int* sel = tables->sel;

	// Bound and score differ only by rounding (~1e-13 here), far inside the
	// margin.
	const int prune = (mode == KNOWN_FIND_ALL && fp != NULL) || mode == KNOWN_MASK_MOTIFS;
	const double cut = motif->threshold - 1e-6;
	// find with prepared 4-mers bounds by two block lookups (the fast path
	// below); mask, which rewrites the sequence as it goes, by the 4 ranked
	// positions in the general loops.
	const int useBlocks = prune && mode == KNOWN_FIND_ALL && kmer5 != NULL && tables->nblocks > 0;
	const unsigned short* km = kmer5;
	const double blockRest = tables->blockRest;
	const int a0 = tables->blockStart[0];
	const int a1 = tables->blockStart[1];
	const int len0 = tables->blockLen[0];
	const int len1 = tables->blockLen[1];
	const double* bf0 = tables->bfwd;
	const double* bf1 = tables->bfwd+3125;
	const double* br0 = tables->brev;
	const double* br1 = tables->brev+3125;

	bestScore = -1e10;
	int numSites = 0;

	if (useBlocks) {
		// find's fast path: a branch-free pass bounds every window and keeps
		// those that could reach the threshold; only they are summed, in the
		// order the general loops below visit them.
		static thread_local int* cand = NULL;
		static thread_local int candSize = 0;
		if (candSize < length+1) {
			delete []cand;
			candSize = length+1;
			cand = new int[candSize];
		}
		if (strand == MOTIF_STRAND_POS || strand==MOTIF_STRAND_BOTH) {
			const int n = length-mlen+1;
			const unsigned short* k0 = km+a0;
			const unsigned short* k1 = km+a1;
			int nc = 0;
			for (int i=0;i<n;i++) {
				double bound = blockRest + bf0[k0[i]] + bf1[k1[i]];
				cand[nc] = i;
				nc += !(bound < cut);
			}
			for (int ci=0;ci<nc;ci++) {
				const int i = cand[ci];
				const unsigned char* c = code+i;
				double score = 0.0;
				for (int j=0;j<mlen;j++) {
					score += fwd[j*5+c[j]];
				}
				if (score == -INFINITY) continue;
				if (score > bestScore) bestScore = score;
				if (score >= motif->threshold) {
					sitePositions[numSites] = i+offset;
					siteStrands[numSites] = MOTIF_STRAND_POS;
					numSites++;
					scanOutSite(fp,name,i+offset,&(s[i]),motif->length,motif->name,'+',score);
				}
			}
		}
		if ((strand == MOTIF_STRAND_BOTH || strand==MOTIF_STRAND_NEG) && length >= mlen) {
			// Windows end at i = mlen-1 .. length-1, and block q's 5-mer
			// starts at i-a_q-len_q+1. Collected ascending, summed descending.
			const int first = mlen-1;
			const int n = length-first;
			const unsigned short* k0 = km+first-a0-len0+1;
			const unsigned short* k1 = km+first-a1-len1+1;
			int nc = 0;
			for (int t=0;t<n;t++) {
				double bound = blockRest + br0[k0[t]] + br1[k1[t]];
				cand[nc] = first+t;
				nc += !(bound < cut);
			}
			for (int ci=nc-1;ci>=0;ci--) {
				const int i = cand[ci];
				const unsigned char* c = code+i;
				double score = 0.0;
				for (int j=0;j<mlen;j++) {
					score += rev[j*5+c[-j]];
				}
				if (score == -INFINITY) continue;
				if (score > bestScore) bestScore = score;
				if (score >= motif->threshold) {
					sitePositions[numSites] = i+offset;
					siteStrands[numSites] = MOTIF_STRAND_NEG;
					numSites++;
					scanOutSite(fp,name,i+offset,&(s[i-motif->length+1]),motif->length,motif->name,'-',score);
				}
			}
		}
		return numSites;
	}

	// Without the threshold bound (the modes that keep the best score, and
	// ZOOPS), every window's block bound is computed first, and only windows
	// whose bound reaches need are summed, in the usual order. need is the
	// threshold for ZOOPS, which reads no best score; otherwise the lower of
	// the threshold and the score of the window with the highest bound, since
	// a window below that can be neither a site nor the best score, nor tie it.
	const int first = mlen-1;
	const int nwin = length-mlen+1;
	int nf = (nwin > 0) ? nwin : 0;
	int nr = nf;
	const int* fc = NULL;
	const int* rc = NULL;
	if (!prune && tables->nblocks > 0 && nwin > 0) {
		static thread_local unsigned short* k5 = NULL;
		static thread_local double* fb = NULL;
		static thread_local double* rb = NULL;
		static thread_local int* fcand = NULL;
		static thread_local int* rcand = NULL;
		static thread_local int filterSize = 0;
		if (filterSize < length+1) {
			delete []k5;
			delete []fb;
			delete []rb;
			delete []fcand;
			delete []rcand;
			filterSize = length+1;
			k5 = new unsigned short[filterSize];
			fb = new double[filterSize];
			rb = new double[filterSize];
			fcand = new int[filterSize];
			rcand = new int[filterSize];
		}
		// kmerCodes5's codes, each position on its own (code is the scratch
		// copy here, padded with 4s)
		for (int i=0;i<length;i++) {
			k5[i] = (unsigned short)(code[i]*625 + code[i+1]*125 + code[i+2]*25 + code[i+3]*5 + code[i+4]);
		}
		const int doF = strand == MOTIF_STRAND_POS || strand==MOTIF_STRAND_BOTH;
		const int doR = strand == MOTIF_STRAND_BOTH || strand==MOTIF_STRAND_NEG;
		if (doF) {
			const unsigned short* k0 = k5+a0;
			const unsigned short* k1 = k5+a1;
			for (int i=0;i<nwin;i++) fb[i] = blockRest + bf0[k0[i]] + bf1[k1[i]];
		}
		if (doR) {
			// as find's fast path: the window ending at first+t
			const unsigned short* k0 = k5+first-a0-len0+1;
			const unsigned short* k1 = k5+first-a1-len1+1;
			for (int t=0;t<nwin;t++) rb[t] = blockRest + br0[k0[t]] + br1[k1[t]];
		}
		double need = motif->threshold;
		if (mode != KNOWN_JUST_SCORE_ZOOPS) {
			// the highest bound (NaNs ignored), by four independent maxima, then
			// the first window that has it
			double tm[4] = {-INFINITY, -INFINITY, -INFINITY, -INFINITY};
			for (int st=0;st<2;st++) {
				if (st == 0 ? !doF : !doR) continue;
				const double* b = (st == 0) ? fb : rb;
				int i = 0;
				for (;i+4<=nwin;i+=4) {
					for (int u=0;u<4;u++) tm[u] = (b[i+u] > tm[u]) ? b[i+u] : tm[u];
				}
				for (;i<nwin;i++) tm[0] = (b[i] > tm[0]) ? b[i] : tm[0];
			}
			double top = tm[0];
			for (int u=1;u<4;u++) top = (tm[u] > top) ? tm[u] : top;
			int topAt = -1;
			int topRev = 0;
			if (top > -INFINITY) {
				if (doF) {
					for (int i=0;i<nwin;i++) {
						if (fb[i] == top) {
							topAt = i;
							break;
						}
					}
				}
				if (topAt < 0 && doR) {
					for (int t=0;t<nwin;t++) {
						if (rb[t] == top) {
							topAt = t;
							topRev = 1;
							break;
						}
					}
				}
			}
			double best = -INFINITY;
			if (topAt >= 0) {
				best = 0.0;
				if (topRev) {
					const unsigned char* c = code+first+topAt;
					for (int j=0;j<mlen;j++) best += rev[j*5+c[-j]];
				} else {
					const unsigned char* c = code+topAt;
					for (int j=0;j<mlen;j++) best += fwd[j*5+c[j]];
				}
			}
			// a best that is not a finite number filters nothing
			if (!(best > -INFINITY && best < INFINITY)) {
				need = -INFINITY;
			} else if (best < need) {
				need = best;
			}
		}
		need -= tables->margin;
		if (doF) {
			int n = 0;
			for (int i=0;i<nwin;i++) {
				fcand[n] = i;
				n += !(fb[i] < need);
			}
			nf = n;
			fc = fcand;
		}
		if (doR) {
			int n = 0;
			for (int t=0;t<nwin;t++) {
				rcand[n] = first+t;
				n += !(rb[t] < need);
			}
			nr = n;
			rc = rcand;
		}
	}

	if (strand == MOTIF_STRAND_POS || strand==MOTIF_STRAND_BOTH) {
		for (int ci=0;ci<nf;ci++) {
			const int i = (fc != NULL) ? fc[ci] : ci;
			const unsigned char* c = code+i;
			if (prune) {
				double bound = restOfBest;
				for (int t=0;t<K;t++) {
					bound += fwd[sel[t]*5+c[sel[t]]];
				}
				if (bound < cut) continue;
			}
			double score = 0.0;
			for (int j=0;j<mlen;j++) {
				score += fwd[j*5+c[j]];
			}
			if (score == -INFINITY) continue;
			if (score > bestScore) {
				bestScore = score;
				if (mode == KNOWN_MSCORE_MOTIFS) {
					sitePositions[0] = i+offset;
					siteStrands[0] = MOTIF_STRAND_POS;
					numSites=1;
				}
			}
			if (score >= motif->threshold) {
				sitePositions[numSites] = i+offset;
				siteStrands[numSites] = MOTIF_STRAND_POS;
				numSites++;
				if (mode==KNOWN_JUST_SCORE_ZOOPS && fp == NULL) {
					break;
				}
				if (mode == KNOWN_MASK_MOTIFS) {
					for (int j=0;j<motif->length;j++) {
						s[i+j] = 'N';
						code[i+j] = 4;
					}
				} else if (mode == KNOWN_MSCORE_MOTIFS) {
				} else if (fp != NULL) {
					scanOutSite(fp,name,i+offset,&(s[i]),motif->length,motif->name,'+',score);
				}
			}
		}
	}
	if (mode == KNOWN_JUST_SCORE_ZOOPS && numSites > 0) {
		return 1;
	}

	if (strand == MOTIF_STRAND_BOTH || strand==MOTIF_STRAND_NEG) {
		for (int ci=nr-1;ci>=0;ci--) {
			const int i = (rc != NULL) ? rc[ci] : first+ci;
			const unsigned char* c = code+i;
			if (prune) {
				double bound = restOfBest;
				for (int t=0;t<K;t++) {
					bound += rev[sel[t]*5+c[-sel[t]]];
				}
				if (bound < cut) continue;
			}
			double score = 0.0;
			for (int j=0;j<mlen;j++) {
				score += rev[j*5+c[-j]];
			}
			if (score == -INFINITY) continue;
			if (score > bestScore) {
				bestScore = score;
				if (mode == KNOWN_MSCORE_MOTIFS) {
					sitePositions[0] = i+offset;
					siteStrands[0] = MOTIF_STRAND_NEG;
					numSites=1;
				}
			}
			if (score >= motif->threshold) {
				sitePositions[numSites] = i+offset;
				siteStrands[numSites] = MOTIF_STRAND_NEG;
				numSites++;
				if (mode==KNOWN_JUST_SCORE_ZOOPS && fp == NULL) {
					break;
				}
				if (mode == KNOWN_MASK_MOTIFS) {
					for (int j=0;j<motif->length;j++) {
						s[i-j] = 'N';
						code[i-j] = 4;
					}
				} else if (mode == KNOWN_MSCORE_MOTIFS) {
				} else if (fp != NULL) {
					scanOutSite(fp,name,i+offset,&(s[i-motif->length+1]),motif->length,motif->name,'-',score);
				}
			}
		}
	}
	return numSites;
}





//=================--------------------------------
LinkedSequence::LinkedSequence(Sequence* seq, LinkedSequence* parent) {
	s = seq;
	link = parent;
}
LinkedSequence::~LinkedSequence() {
	//if (link != NULL) delete link;
}
SequenceArray::SequenceArray() {
	seqs = NULL;
	numSeqs= 0;
	numTarget=0.0;
	numBackground=0.0;
	bpTarget=0;
	bpBackground=0;
	link=NULL;
	longestSequence = 0;
	scorer = NULL;
	maxBackgroundPercent = MOTIF_MAX_BACKGROUND_PERCENT;
	maxAutoNormIters = MAX_OLIGO_NORMALIZATION_ITERATIONS;
	chunks = NULL;
	numChunks = 0;
	index5 = 0;
	chunks6 = NULL;
	count6 = NULL;
	chunks7 = NULL;
	count7 = NULL;
	maxChunk6Entries = 0;
	seqStatsUnread = 0;
	textBlocks = NULL;
	numTextBlocks = 0;
	objBlocks = NULL;
	objCounts = NULL;
	numObjBlocks = 0;
}
SequenceArray::~SequenceArray() {
	//if (link != NULL) delete link;
	delete []chunks;
	delete []chunks6;
	delete []count6;
	delete []chunks7;
	delete []count7;
	deleteLinks();
	for (int i=0;i<numSeqs;i++) {
		// a Sequence built in an object block is destroyed in place
		int inBlock = 0;
		for (int b=0;b<numObjBlocks;b++) {
			if (seqs[i] >= objBlocks[b] && seqs[i] < objBlocks[b]+objCounts[b]) inBlock = 1;
		}
		if (inBlock) {
			seqs[i]->~Sequence();
		} else {
			delete seqs[i];
		}
	}
	delete []seqs;
	for (int i=0;i<numTextBlocks;i++) free(textBlocks[i]);
	delete []textBlocks;
	for (int i=0;i<numObjBlocks;i++) free(objBlocks[i]);
	delete []objBlocks;
	delete []objCounts;
	if (scorer != NULL) delete scorer;
}
void SequenceArray::deleteLinks() {
	LinkedSequence* cur = link;
	while (cur != NULL) {
		LinkedSequence* next = cur->link;
		delete cur;
		cur  = next;
	}
}
void SequenceArray::initializeScorer(int scoringMethod,int zoopsFlag,unsigned long long int maxCacheSize) {
	scorer = new EnrichmentScorer();
	scorer->initializeEnrichmentScoring(scoringMethod,numTarget, numBackground,
															zoopsFlag,maxCacheSize);
}
void SequenceArray::addSequence(Sequence* s) {
	link = new LinkedSequence(s, link);
	numSeqs++;
	//fprintf(stderr, "%s\t%d\n", s->name, numSeqs);
}
void SequenceArray::optimizeArray() {
	seqs = new Sequence*[numSeqs];
	//fprintf(stderr, "numSeqs(optimize)=%d\n", numSeqs);
	numSeqs=0;
	numBackground = 0.0;
	numTarget = 0.0;
	bpTarget = 0;
	bpBackground = 0;
	longestSequence = 0;
	LinkedSequence* ls = link;
	while (ls != NULL) {
		seqs[numSeqs] = ls->s;
		seqs[numSeqs]->id = numSeqs;
		if (seqs[numSeqs]->length > longestSequence) {
			longestSequence = seqs[numSeqs]->length;
		}
		//fprintf(stderr, "New long: %d, %s %d\n", seqs[numSeqs]->length, seqs[numSeqs]->s,longestSequence);

		if (seqs[numSeqs]->g == 0) {
			numBackground += seqs[numSeqs]->w;
			bpBackground += seqs[numSeqs]->length;
		} else {
			numTarget += seqs[numSeqs]->w;
			bpTarget += seqs[numSeqs]->length;
		}
		ls = ls->link;
		numSeqs++;
	}
	longestSequence *= 2;
	//fprintf(stderr, "numSeqs(optimize)=%d\n", numSeqs);
}

void SequenceArray::parseFasta2SeqAndGroupFiles(char* inputFASTAfile,
							char* bgFASTAfile, char* seqFile, char* groupFile) {
	
	FILE* fpseq = fopen(seqFile,"w");
	if (fpseq == NULL) {
		fprintf(stderr, "!!! Could not open a temporary file (%s) for writing!!!\n", seqFile);
		exit(0);
	}
	FILE* fpgroup = fopen(groupFile,"w");
	if (fpgroup == NULL) {
		fprintf(stderr, "!!! Could not open a temporary file (%s) for writing!!!\n", groupFile);
		exit(0);
	}

	Inttable* names = new Inttable(1000000);	

	char* buf = new char[MOTIF2_BUFFER];

	char* curSeq = new char[1000000];
	char* curName = new char[MOTIF2_BUFFER];

	for (int i=0;i<2;i++) {

		if (i==0 && inputFASTAfile == NULL) continue;
		if (i==1 && bgFASTAfile == NULL) continue;
	
		FILE* fp = NULL;
		int g = 1;
		if (i==0) {		
			fp = fopen(inputFASTAfile,"r");
			if (fp == NULL) {
				fprintf(stderr, "Could not open FASTA file %s !!!\n", inputFASTAfile);
				exit(0);
			}
		} else {
			fp = fopen(bgFASTAfile,"r");
			if (fp == NULL) {
				fprintf(stderr, "Could not open FASTA file %s !!!\n", bgFASTAfile);
				exit(0);
			}
			g = 0;
		}
		curSeq[0] = '\0';
		curName[0] = '\0';
		buf[MOTIF2_BUFFER-1]='\0';
		while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
			if (buf[MOTIF2_BUFFER-1]!='\0') {
				fprintf(stderr, "!!!Error: Input sequence is too long (>%d)!!!\n", MOTIF2_BUFFER);
			}
			buf[MOTIF2_BUFFER-1]='\0';
			int lineLength = strlen(buf);
			//fprintf(stderr, "lineLength=%d\n", lineLength);
			if (lineLength > 0) buf[lineLength-1]='\0';
			//for carriage returns
			if (lineLength > 1 && buf[lineLength-2] == 13) buf[lineLength-2]='\0';
			

			if (buf[0] == '>') {
				if (curSeq[0] != '\0' && curName[0] != '\0') {
					fprintf(fpgroup,"%s\t%d\t1.0\n",curName,g);
					fprintf(fpseq,"%s\t%s\n",curName,curSeq);
				}
				curName[0] = '\0';
				curSeq[0] = '\0';
				if (lineLength > 1) {
					strcpy(curName,&(buf[1]));
					int c = names->search(curName);
					if (c != EMPTY_INT) {
						names->insert(++c,curName);
						sprintf(curName,"%s-%d",&(buf[1]),c);
					} else {
						names->insert(1,curName);
					}
				}
			} else {
				strcat(curSeq,buf);
			}
		}
		if (curSeq[0] != '\0' && curName[0] != '\0') {
			fprintf(fpgroup,"%s\t%d\t1.0\n",curName,g);
			fprintf(fpseq,"%s\t%s\n",curName,curSeq);
		}
		fclose(fp);

	}

	fclose(fpseq);
	fclose(fpgroup);

	delete []curSeq;
	delete []curName;
	delete names;
	delete []buf;
}


// fgets(buf, MOTIF2_BUFFER, fp) over a block buffer: the same lines (at most
// MOTIF2_BUFFER-1 bytes, ending after the first '\n' when there is one),
// without a library call per line. A line stays valid until the next call.
struct FastaLineReader {
	FILE* fp;
	std::vector<char> buf;
	size_t start;
	size_t end;
	bool eof;
	FastaLineReader(FILE* f) : fp(f), buf((size_t)1 << 20), start(0), end(0), eof(false) {}
	bool next(char*& line, size_t& len) {
		const size_t most = MOTIF2_BUFFER-1;
		while (1) {
			size_t avail = end-start;
			char* p = buf.data()+start;
			char* nl = (char*)memchr(p, '\n', avail < most ? avail : most);
			if (nl != NULL) {
				len = nl-p+1;
			} else if (avail >= most) {
				len = most;
			} else if (eof) {
				if (avail == 0) return false;
				len = avail;
			} else {
				if (start > 0) {
					memmove(buf.data(), p, avail);
					start = 0;
					end = avail;
				}
				if (buf.size()-end < ((size_t)1 << 20)) buf.resize(end + ((size_t)1 << 20));
				size_t got = fread(buf.data()+end, 1, buf.size()-end, fp);
				if (got == 0) eof = true;
				end += got;
				continue;
			}
			line = p;
			start += len;
			return true;
		}
	}
};
struct FastaRecord {
	std::string name;
	std::string seq;
	int g;
};

// -i/-b FASTA read in memory. The records, their names (duplicates renamed
// name-2, name-3, ...), groups and the stderr lines are those of
// parseFasta2SeqAndGroupFiles followed by readSequenceFiles. When no name or
// sequence holds a tab or CR and no name repeats, the .seq/.group round trip
// gives each record its own name, sequence, group and weight 1.0, so the
// Sequences are built directly; otherwise the lines that round trip would
// write are parsed with readSequenceFiles' code, as before.
static void copyCleaned(char* dst, const char* src, size_t n);
void Sequence::setSeqCleaned(const char* n, size_t nameLen, const char* seq, size_t seqLen,
		int group, double weight) {
	if (textShared) {
		s = NULL;
		name = NULL;
		textShared = 0;
	}
	if (name != NULL) delete []name;
	name = new char[nameLen+1];
	memcpy(name, n, nameLen);
	name[nameLen] = '\0';
	w = weight;
	g = group;
	if (kmer5 != NULL) delete []kmer5;
	kmer5 = NULL;
	if (kmer6 != NULL) delete []kmer6;
	kmer6 = NULL;
	if (codes != NULL && !codesShared) delete []codes;
	codes = NULL;
	codesShared = 0;
	if (s != NULL) delete []s;
	length = (int)seqLen;
	s = new char[seqLen+1];
	copyCleaned(s, seq, seqLen);
	s[seqLen] = '\0';
}
void Sequence::setSeqCleanedAt(const char* n, size_t nameLen, const char* seq, size_t seqLen,
		int group, double weight, char* at) {
	if (textShared) {
		s = NULL;
		name = NULL;
	}
	if (name != NULL) delete []name;
	if (s != NULL) delete []s;
	w = weight;
	g = group;
	if (kmer5 != NULL) delete []kmer5;
	kmer5 = NULL;
	if (kmer6 != NULL) delete []kmer6;
	kmer6 = NULL;
	if (codes != NULL && !codesShared) delete []codes;
	codes = NULL;
	codesShared = 0;
	name = at;
	memcpy(name, n, nameLen);
	name[nameLen] = '\0';
	length = (int)seqLen;
	s = at + nameLen + 1;
	copyCleaned(s, seq, seqLen);
	s[seqLen] = '\0';
	textShared = 1;
}

// One piece of a FASTA file in memory, starting at a header line or at the
// file's start: per header, its raw name and the kept bytes of the lines
// after it (as the serial reader keeps them: the last byte dropped, and a CR
// before it), back to back in `text`; lines before the first header form an
// entry with no header. ok is 0 when a line is longer than fgets(buf,
// MOTIF2_BUFFER) returns whole; tabOrCr when a kept byte is a tab or CR.
struct FastaSpan {
	size_t nameAt;
	size_t nameLen;
	size_t seqAt;
	size_t seqLen;
	int header;
	// FNV-1a of the raw name (never 0), for the renaming table
	unsigned long long nameHash;
};
static unsigned long long fastaNameHash(const char* p, size_t n) {
	unsigned long long h = 1469598103934665603ULL;
	for (size_t i=0;i<n;i++) {
		h ^= (unsigned char)p[i];
		h *= 1099511628211ULL;
	}
	return (h == 0) ? 1 : h;
}
// Raw name -> times seen, for renaming: open addressing on the parse's
// hashes; keys are views of the pieces' text, which outlive the table.
struct FastaNameTable {
	std::vector<unsigned long long> hashes;
	std::vector<std::string_view> keys;
	std::vector<int> counts;
	size_t used;
	FastaNameTable() : used(0) {
		hashes.assign(1024, 0);
		keys.resize(1024);
		counts.assign(1024, 0);
	}
	// room for more names without growing
	void reserve(size_t more) {
		size_t want = hashes.size();
		while (2*(used+more) > want) want *= 2;
		if (want != hashes.size()) growTo(want);
	}
	// The slot of name (hash h): found, or empty and claimed for it.
	size_t slot(unsigned long long h, std::string_view name, int& found) {
		if (2*(used+1) > hashes.size()) grow();
		const size_t mask = hashes.size()-1;
		size_t i = (size_t)h & mask;
		while (hashes[i] != 0) {
			if (hashes[i] == h && keys[i] == name) {
				found = 1;
				return i;
			}
			i = (i+1) & mask;
		}
		hashes[i] = h;
		keys[i] = name;
		counts[i] = 0;
		used++;
		found = 0;
		return i;
	}
	void grow() {
		growTo(hashes.size()*2);
	}
	// rehashed into n slots (a power of two)
	void growTo(size_t n) {
		std::vector<unsigned long long> oh;
		std::vector<std::string_view> ok;
		std::vector<int> oc;
		oh.swap(hashes);
		ok.swap(keys);
		oc.swap(counts);
		hashes.assign(n, 0);
		keys.resize(n);
		counts.assign(n, 0);
		const size_t mask = hashes.size()-1;
		for (size_t j=0;j<oh.size();j++) {
			if (oh[j] == 0) continue;
			size_t i = (size_t)oh[j] & mask;
			while (hashes[i] != 0) i = (i+1) & mask;
			hashes[i] = oh[j];
			keys[i] = ok[j];
			counts[i] = oc[j];
		}
	}
};
// A FASTA file's bytes.
struct FastaFileBuffer {
	char* data;
	size_t size;
	FastaFileBuffer() : data(NULL), size(0) {}
	~FastaFileBuffer() {
		free(data);
	}
	FastaFileBuffer(const FastaFileBuffer&) = delete;
	FastaFileBuffer& operator=(const FastaFileBuffer&) = delete;
	FastaFileBuffer(FastaFileBuffer&& o) : data(o.data), size(o.size) {
		o.data = NULL;
		o.size = 0;
	}
};
template <typename F>
static void parallelFor(int n, int threads, F fn);
// A regular file, read with one pread per thread into memory nothing
// initializes first (so the copies and the page faults run in parallel).
// Returns 0, reading nothing, when the file is not regular or not the size
// fstat gave.
static int readFileParallel(FILE* fp, FastaFileBuffer& out, int threads) {
	const int fd = fileno(fp);
	struct stat st;
	if (fstat(fd, &st) != 0 || !S_ISREG(st.st_mode) || st.st_size <= 0) return 0;
	const size_t size = (size_t)st.st_size;
	char* buf = chunkArray<char>((long)size);
	const int parts = (threads > 1) ? threads : 1;
	std::atomic<int> bad(0);
	parallelFor(parts, threads, [&](int t) {
		size_t lo = size*(size_t)t/parts;
		const size_t hi = size*(size_t)(t+1)/parts;
		while (lo < hi) {
			const ssize_t r = pread(fd, buf+lo, hi-lo, (off_t)lo);
			if (r < 0 && errno == EINTR) continue;
			if (r <= 0) {
				bad = 1;
				return;
			}
			lo += (size_t)r;
		}
	});
	char extra = 0;
	if (!bad && pread(fd, &extra, 1, (off_t)size) != 0) bad = 1;
	if (bad) {
		free(buf);
		return 0;
	}
	out.data = buf;
	out.size = size;
	return 1;
}
struct FastaPiece {
	std::vector<char, HugeAllocator<char> > text;
	std::vector<FastaSpan> spans;
	int ok;
	int tabOrCr;
};
static void parseFastaPiece(const char* p, const char* end, FastaPiece& out) {
	const size_t most = MOTIF2_BUFFER-1;
	out.text.clear();
	out.spans.clear();
	// room for a record every 32 bytes (reserved, not touched)
	out.spans.reserve((size_t)(end-p)/32 + 16);
	out.ok = 1;
	out.tabOrCr = 0;
	out.text.reserve(end-p);
	while (p < end) {
		const char* nl = (const char*)memchr(p, '\n', end-p);
		size_t len = (nl != NULL) ? (size_t)(nl-p+1) : (size_t)(end-p);
		if (len > most) {
			out.ok = 0;
			return;
		}
		size_t keep = len-1;
		if (len > 1 && p[len-2] == 13) keep = len-2;
		if (memchr(p, '\t', keep) != NULL || memchr(p, '\r', keep) != NULL) out.tabOrCr = 1;
		if (keep > 0 && p[0] == '>') {
			FastaSpan sp;
			sp.nameAt = out.text.size();
			sp.nameLen = keep-1;
			sp.nameHash = fastaNameHash(p+1, keep-1);
			out.text.insert(out.text.end(), p+1, p+keep);
			sp.seqAt = out.text.size();
			sp.seqLen = 0;
			sp.header = 1;
			out.spans.push_back(sp);
		} else {
			if (out.spans.empty()) {
				FastaSpan sp = {out.text.size(), 0, out.text.size(), 0, 0, 1};
				out.spans.push_back(sp);
			}
			out.text.insert(out.text.end(), p, p+keep);
			out.spans.back().seqLen += keep;
		}
		p += len;
	}
}

// A record to build: its final name (a view of a piece's text, or
// renamed[renamedAt] when it got a -N suffix), its kept sequence bytes, its
// group.
struct FastaKept {
	const char* name;
	size_t nameLen;
	const char* seq;
	size_t seqLen;
	int g;
	long renamedAt;
};

// Splits base[0, n) at header lines into pieces, parses them on `threads`
// threads and appends the kept records in file order. names has one table
// per thread, for the names whose hash falls in its partition. Returns 0
// when the serial reader must take every file instead.
static int readFastaPiecesParallel(const char* base, size_t n, int g, int threads,
		std::vector<FastaPiece>& pieces, std::vector<FastaNameTable>& names,
		std::vector<std::string>& renamed, std::vector<FastaKept>& kept, int& anyRenamed,
		int& tabOrCr) {
	if (n > 0 && memchr(base, '\0', n) != NULL) return 0;
	const int count = (threads > 1 ? threads : 1) * 4;
	std::vector<size_t> starts(count+1, n);
	starts[0] = 0;
	for (int t=1;t<count;t++) {
		size_t at = n*(size_t)t/count;
		if (at < starts[t-1]) at = starts[t-1];
		size_t s = n;
		if (at > 0) {
			const char* q = base+at-1;
			while (1) {
				q = (const char*)memchr(q, '\n', base+n-q);
				if (q == NULL || q+1 >= base+n) break;
				if (q[1] == '>') {
					s = q+1-base;
					break;
				}
				q++;
			}
		}
		starts[t] = s;
	}
	const size_t first = pieces.size();
	pieces.resize(first + count);
	parallelFor(count, threads, [&](int t) {
		FastaPiece& pc = pieces[first+t];
		parseFastaPiece(base+starts[t], base+starts[t+1], pc);
		// a later piece starts at a header line, unless the line is a bare
		// '>' at the end of the file, which keeps nothing
		if (pc.ok && t > 0 && !pc.spans.empty() && !pc.spans[0].header && pc.spans[0].seqLen > 0) pc.ok = 0;
	});
	for (int t=0;t<count;t++) {
		if (!pieces[first+t].ok) return 0;
		if (pieces[first+t].tabOrCr) tabOrCr = 1;
	}

	// Each header's count of its name so far, the count the serial
	// reader's table reaches there: a name's headers all fall in one hash
	// partition, and each thread takes its partition's headers in file
	// order (the tables carry over to the next file).
	const int parts = (int)names.size();
	std::vector<std::vector<int> > seen(count);
	std::vector<size_t> perPart(parts, 0);
	for (int t=0;t<count;t++) {
		const FastaPiece& pc = pieces[first+t];
		seen[t].assign(pc.spans.size(), 0);
		for (size_t k=0;k<pc.spans.size();k++) {
			if (pc.spans[k].header) perPart[(pc.spans[k].nameHash >> 40) % parts]++;
		}
	}
	parallelFor(parts, threads, [&](int part) {
		FastaNameTable& table = names[part];
		table.reserve(perPart[part]);
		for (int t=0;t<count;t++) {
			const FastaPiece& pc = pieces[first+t];
			for (size_t k=0;k<pc.spans.size();k++) {
				const FastaSpan& sp = pc.spans[k];
				if (!sp.header || (int)((sp.nameHash >> 40) % parts) != part) continue;
				std::string_view raw(pc.text.data()+sp.nameAt, sp.nameLen);
				int found = 0;
				const size_t at = table.slot(sp.nameHash, raw, found);
				seen[t][k] = found ? ++(table.counts[at]) : (table.counts[at] = 1);
			}
		}
	});

	{
		size_t headers = 0;
		for (int part=0;part<parts;part++) headers += perPart[part];
		kept.reserve(kept.size() + headers);
	}
	// The serial reader's order: a header renames its name whether or not its
	// record is kept; a record is kept when its name and sequence are not empty.
	FastaKept cur = {NULL, 0, NULL, 0, g, -1};
	int curNamed = 0;
	for (int t=0;t<count;t++) {
		FastaPiece& pc = pieces[first+t];
		for (size_t k=0;k<pc.spans.size();k++) {
			const FastaSpan& sp = pc.spans[k];
			if (!sp.header) continue;
			if (curNamed && cur.nameLen > 0 && cur.seqLen > 0) kept.push_back(cur);
			std::string_view raw(pc.text.data()+sp.nameAt, sp.nameLen);
			if (seen[t][k] > 1) {
				int c = seen[t][k];
				char suffix[32];
				sprintf(suffix, "-%d", c);
				renamed.push_back(std::string(raw) + suffix);
				anyRenamed = 1;
				cur.name = NULL;
				cur.nameLen = renamed.back().size();
				cur.renamedAt = (long)renamed.size()-1;
			} else {
				cur.name = raw.data();
				cur.nameLen = raw.size();
				cur.renamedAt = -1;
			}
			cur.seq = pc.text.data()+sp.seqAt;
			cur.seqLen = sp.seqLen;
			curNamed = 1;
		}
	}
	if (curNamed && cur.nameLen > 0 && cur.seqLen > 0) kept.push_back(cur);
	return 1;
}

void SequenceArray::readFastaDirect(char* inputFASTAfile, char* bgFASTAfile, int threads) {
	// Each file whole, in memory, parsed in parallel; a file the parallel
	// parse cannot take, or a tab or CR in a kept byte, sends every file
	// through the serial reader below.
	{
		std::vector<FastaFileBuffer> files;
		std::vector<FastaPiece> pieces;
		std::vector<FastaNameTable> names((threads > 1) ? threads : 1);
		std::vector<std::string> renamed;
		std::vector<FastaKept> kept;
		int anyRenamed = 0;
		int tabOrCr = 0;
		int ok = 1;
		files.reserve(2);
		for (int i=0;i<2 && ok;i++) {
			if (i==0 && inputFASTAfile == NULL) continue;
			if (i==1 && bgFASTAfile == NULL) continue;
			FILE* fp = fopen(i == 0 ? inputFASTAfile : bgFASTAfile, "r");
			if (fp == NULL) {
				fprintf(stderr, "Could not open FASTA file %s !!!\n", i == 0 ? inputFASTAfile : bgFASTAfile);
				exit(0);
			}
			files.push_back(FastaFileBuffer());
			FastaFileBuffer& file = files.back();
			if (!readFileParallel(fp, file, threads)) {
				std::vector<char> data;
				struct stat st;
				size_t size = (size_t)1 << 20;
				if (fstat(fileno(fp), &st) == 0 && S_ISREG(st.st_mode) && st.st_size > 0) size = (size_t)st.st_size + 1;
				data.resize(size);
				size_t got = 0;
				while (1) {
					size_t r = fread(data.data()+got, 1, data.size()-got, fp);
					got += r;
					if (got < data.size()) break;
					data.resize(data.size()*2);
				}
				file.data = chunkArray<char>((long)got);
				if (got > 0) memcpy(file.data, data.data(), got);
				file.size = got;
			}
			fclose(fp);
			if (!readFastaPiecesParallel(file.data, file.size, i == 0 ? 1 : 0, threads, pieces, names, renamed,
					kept, anyRenamed, tabOrCr)) ok = 0;
		}
		// renamed names live in `renamed`, which no longer grows
		for (size_t k=0;k<kept.size() && ok;k++) {
			if (kept[k].renamedAt >= 0) {
				kept[k].name = renamed[kept[k].renamedAt].data();
				kept[k].nameLen = renamed[kept[k].renamedAt].size();
			}
		}
		// final names are distinct unless a rename happened (a renamed name
		// can equal another record's)
		if (ok && !tabOrCr && anyRenamed) {
			std::unordered_set<std::string_view> seen;
			seen.reserve(kept.size()*2);
			for (size_t k=0;k<kept.size() && ok;k++) {
				if (!seen.insert(std::string_view(kept[k].name, kept[k].nameLen)).second) ok = 0;
			}
		}
		if (ok && !tabOrCr) {
			fprintf(stderr, "\tReading input files...\n");
			std::vector<Sequence*> built(kept.size(), (Sequence*)NULL);
			// every record's name and sequence in one block on huge
			// pages: the candidate walk reads them at random
			std::vector<size_t> textAt(kept.size()+1, 0);
			for (size_t k=0;k<kept.size();k++) textAt[k+1] = textAt[k] + kept[k].nameLen + 1 + kept[k].seqLen + 1;
			char* text = chunkArray<char>((long)textAt[kept.size()]);
			Sequence* objs = chunkArray<Sequence>((long)kept.size());
			{
				char** grown = new char*[numTextBlocks+1];
				for (int i=0;i<numTextBlocks;i++) grown[i] = textBlocks[i];
				grown[numTextBlocks++] = text;
				delete []textBlocks;
				textBlocks = grown;
				Sequence** og = new Sequence*[numObjBlocks+1];
				long* oc = new long[numObjBlocks+1];
				for (int i=0;i<numObjBlocks;i++) {
					og[i] = objBlocks[i];
					oc[i] = objCounts[i];
				}
				og[numObjBlocks] = objs;
				oc[numObjBlocks++] = (long)kept.size();
				delete []objBlocks;
				delete []objCounts;
				objBlocks = og;
				objCounts = oc;
			}
			const int blocks = (int)((kept.size() + 4095) / 4096);
			parallelFor(blocks, threads, [&](int b) {
				const size_t lo = (size_t)b*4096;
				const size_t hi = (lo+4096 < kept.size()) ? lo+4096 : kept.size();
				for (size_t k=lo;k<hi;k++) {
					Sequence* s = new (objs + k) Sequence();
					s->setSeqCleanedAt(kept[k].name, kept[k].nameLen, kept[k].seq, kept[k].seqLen,
							kept[k].g, 1.0, text + textAt[k]);
					built[k] = s;
				}
			});
			int totalSequences = 0;
			if (numSeqs == 0 && link == NULL) {
				// the array optimizeArray would build from the list
				// (newest first), with its sums in the same order, without
				// one list node per sequence
				const int n = (int)kept.size();
				totalSequences = n;
				seqs = new Sequence*[n > 0 ? n : 1];
				numBackground = 0.0;
				numTarget = 0.0;
				bpTarget = 0;
				bpBackground = 0;
				longestSequence = 0;
				for (int j=0;j<n;j++) {
					Sequence* sq = built[n-1-j];
					seqs[j] = sq;
					sq->id = j;
					if (sq->length > longestSequence) longestSequence = sq->length;
					if (sq->g == 0) {
						numBackground += sq->w;
						bpBackground += sq->length;
					} else {
						numTarget += sq->w;
						bpTarget += sq->length;
					}
				}
				numSeqs = n;
				longestSequence *= 2;
				fprintf(stderr, "\t%d total sequences read\n", totalSequences);
				return;
			}
			for (size_t k=0;k<kept.size();k++) {
				totalSequences++;
				addSequence(built[k]);
			}
			fprintf(stderr, "\t%d total sequences read\n", totalSequences);
			optimizeArray();
			return;
		}
	}

	std::vector<FastaRecord> records;
	std::unordered_map<std::string,int> names;
	std::string curName;
	std::string curSeq;

	for (int i=0;i<2;i++) {

		if (i==0 && inputFASTAfile == NULL) continue;
		if (i==1 && bgFASTAfile == NULL) continue;

		FILE* fp = NULL;
		int g = 1;
		if (i==0) {
			fp = fopen(inputFASTAfile,"r");
			if (fp == NULL) {
				fprintf(stderr, "Could not open FASTA file %s !!!\n", inputFASTAfile);
				exit(0);
			}
		} else {
			fp = fopen(bgFASTAfile,"r");
			if (fp == NULL) {
				fprintf(stderr, "Could not open FASTA file %s !!!\n", bgFASTAfile);
				exit(0);
			}
			g = 0;
		}
		curSeq.clear();
		curName.clear();
		FastaLineReader reader(fp);
		char* line = NULL;
		size_t len = 0;
		while (reader.next(line, len)) {
			// what strlen() and the two '\0' stores leave of the line
			size_t lineLength = strnlen(line, len);
			size_t keep = lineLength;
			if (lineLength > 0) keep = lineLength-1;
			if (lineLength > 1 && line[lineLength-2] == 13) keep = lineLength-2;

			if (keep > 0 && line[0] == '>') {
				if (!curSeq.empty() && !curName.empty()) {
					records.push_back(FastaRecord());
					records.back().name.swap(curName);
					records.back().seq.swap(curSeq);
					records.back().g = g;
				}
				curName.clear();
				curSeq.clear();
				if (lineLength > 1) {
					curName.assign(line+1, keep-1);
					std::unordered_map<std::string,int>::iterator it = names.find(curName);
					if (it != names.end()) {
						int c = ++(it->second);
						char suffix[32];
						sprintf(suffix, "-%d", c);
						curName += suffix;
					} else {
						names[curName] = 1;
					}
				}
			} else {
				curSeq.append(line, keep);
			}
		}
		if (!curSeq.empty() && !curName.empty()) {
			records.push_back(FastaRecord());
			records.back().name.swap(curName);
			records.back().seq.swap(curSeq);
			records.back().g = g;
		}
		fclose(fp);
	}

	int direct = 1;
	{
		std::unordered_set<std::string_view> seen;
		seen.reserve(records.size()*2);
		for (size_t k=0;k<records.size() && direct;k++) {
			const std::string& n = records[k].name;
			const std::string& q = records[k].seq;
			if (memchr(n.data(), '\t', n.size()) != NULL || memchr(n.data(), '\r', n.size()) != NULL) direct = 0;
			if (memchr(q.data(), '\t', q.size()) != NULL || memchr(q.data(), '\r', q.size()) != NULL) direct = 0;
			if (!seen.insert(std::string_view(n)).second) direct = 0;
		}
	}
	if (direct) {
		fprintf(stderr, "\tReading input files...\n");
		int totalSequences = 0;
		for (size_t k=0;k<records.size();k++) {
			char* seq = &(records[k].seq[0]);
			cleanSequence(seq);
			totalSequences++;
			Sequence* s = new Sequence(&(records[k].name[0]),seq,records[k].g,1.0);
			addSequence(s);
		}
		fprintf(stderr, "\t%d total sequences read\n", totalSequences);
		optimizeArray();
		return;
	}

	std::vector<std::string> groupLines;
	std::vector<std::string> seqLines;
	for (size_t k=0;k<records.size();k++) {
		char groupTail[64];
		sprintf(groupTail, "\t%d\t1.0\n", records[k].g);
		groupLines.push_back(records[k].name + groupTail);
		seqLines.push_back(records[k].name + "\t" + records[k].seq + "\n");
	}
	std::vector<FastaRecord>().swap(records);

	// readSequenceFiles, over the lines instead of the files
	Hashtable* groupdata = new Hashtable();

	fprintf(stderr, "\tReading input files...\n");

	char** cols = new char*[1000];
	int numCols;
	char* seq = NULL;
	char* name = NULL;
	std::vector<char> line;

	for (size_t k=0;k<groupLines.size();k++) {
		line.assign(groupLines[k].begin(), groupLines[k].end());
		line.push_back('\0');
		splitMotif2(&(line[0]),cols,numCols,'\t');
		if (numCols < 2) continue;
		GroupInfo* gi = new GroupInfo();
		gi->weight = 1.0;
		gi->group = 0;
		name = cols[0];
		if (numCols > 1) {
			sscanf(cols[1],"%d",&(gi->group));
		}
		if (numCols > 2) {
			sscanf(cols[2],"%lf",&(gi->weight));
		}
		groupdata->insert(gi, name);
	}

	int totalSequences = 0;
	for (size_t k=0;k<seqLines.size();k++) {
		line.assign(seqLines[k].begin(), seqLines[k].end());
		line.push_back('\0');
		splitMotif2(&(line[0]),cols,numCols,'\t');
		if (numCols < 2) continue;
		int group = 1;
		double weight = 1.0;
		GroupInfo* gi = (GroupInfo*)groupdata->search(cols[0]);
		if (gi == NULL) continue;
		group = gi->group;
		weight = gi->weight;
		name = cols[0];
		seq = cols[1];
		cleanSequence(seq);
		totalSequences++;
		Sequence* s = new Sequence(name,seq,group,weight);
		addSequence(s);
	}
	fprintf(stderr, "\t%d total sequences read\n", totalSequences);

	optimizeArray();

	char** keys = groupdata->keys();
	for (int i=0;i<groupdata->total;i++) {
		GroupInfo* gi = (GroupInfo*)groupdata->search(keys[i]);
		delete gi;
		delete [](keys[i]);
	}
	delete []keys;
	delete groupdata;
	delete []cols;
}

void SequenceArray::readSequenceFiles(char* seqFile, char* groupFile) {

	Hashtable* groupdata = new Hashtable();

	fprintf(stderr, "\tReading input files...\n");

	char* buf = new char[MOTIF2_BUFFER];
	char** cols = new char*[1000];
	int numCols;
	char* seq = NULL;
	char* name = NULL;


	FILE* fp = NULL;
	if (groupFile != NULL) {	
		fp = fopen(groupFile, "r");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open group file: %s !!!\n", groupFile);
			exit(0);
		}
		buf[MOTIF2_BUFFER-1]='\0';
		while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
			if (buf[MOTIF2_BUFFER-1]!='\0') {
				fprintf(stderr, "!!! Input sequence is too long (>%d)!!!\n", MOTIF2_BUFFER);
			}
			buf[MOTIF2_BUFFER-1]='\0';
			splitMotif2(buf,cols,numCols,'\t');
			if (numCols < 2) continue;
			GroupInfo* gi = new GroupInfo();
			gi->weight = 1.0;
			gi->group = 0;
			name = cols[0];
			if (numCols > 1) {
				sscanf(cols[1],"%d",&(gi->group));
			}
			if (numCols > 2) {
				sscanf(cols[2],"%lf",&(gi->weight));
			}
			groupdata->insert(gi, name);
		}
		fclose(fp);
	}

	fp = fopen(seqFile, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open seq file: %s !!!\n", seqFile);
		exit(0);
	}

	int totalSequences = 0;
	buf[MOTIF2_BUFFER-1]='\0';
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		if (buf[MOTIF2_BUFFER-1]!='\0') {
			fprintf(stderr, "!!! Input sequence is too long (>%d)!!!\n", MOTIF2_BUFFER);
		}
		buf[MOTIF2_BUFFER-1]='\0';
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 2) continue;
		int group = 1;
		double weight = 1.0;
		if (groupFile != NULL) {
			GroupInfo* gi = (GroupInfo*)groupdata->search(cols[0]);
			if (gi == NULL) continue;
			group = gi->group;
			weight = gi->weight;
		}
		name = cols[0];
		seq = cols[1];
		cleanSequence(seq);
		//if (strlen(seq) < 1) continue;
		totalSequences++;
		Sequence* s = new Sequence(name,seq,group,weight);
		//fprintf(stderr, "%s\t%s\t%d\n", name,seq,s->length);
		//fprintf(stderr, "add seqeuene.\n");
		addSequence(s);
	}
	fprintf(stderr, "\t%d total sequences read\n", totalSequences);

	optimizeArray();
	/*
	double oligoPerSeqTarget = totalTargetOligos/numTargets;
	double oligoPerSeqBackround = totalBackgroundOligos/numBackground;
	double newNumBackground = totalBackgroundOligos/oligoPerSeqTarget;

	fprintf(stderr, "\tAutoadjustment for sequence coverage in background: %.2lfx\n",newNumBackground/numBackground);
	numBackground = newNumBackground;
	*/

	char** keys = groupdata->keys();
	for (int i=0;i<groupdata->total;i++) {
		GroupInfo* gi = (GroupInfo*)groupdata->search(keys[i]);
		delete gi;
		delete [](keys[i]);
	}
	delete []keys;
	delete groupdata;
	delete []buf;
	delete []cols;

}
void SequenceArray::printGroupFile(FILE* fp) {
	if (fp == NULL) {
		fprintf(stderr, "!!! Cannot output normalized group file !!!\n");
		return;
	}
	for (int i=0;i<numSeqs;i++) {
		fprintf(fp, "%s\t%d\t%lf\n",seqs[i]->name,seqs[i]->g,seqs[i]->w);
	}
}
void SequenceArray::printSequenceFile(FILE* fp) {
	if (fp == NULL) {
		fprintf(stderr, "!!! Cannot output sequence file !!!\n");
		return;
	}
	for (int i=0;i<numSeqs;i++) {
		fprintf(fp, "%s\t%s\n",seqs[i]->name,seqs[i]->s);
	}
}

void SequenceArray::maskMotifFromOligoArray(Motif* motif, OligoArray* oligos, char strand) {
	motif->logxform();
	motif->initializeSeqStats(numSeqs);
	double score = 0;
	if (longestSequence < 1) longestSequence = 2;
	int *sitePositions = new int[longestSequence];
	char *siteStrands = new char[longestSequence];
	int numPositions2Mask = 0;
	int *positions2Mask = new int[longestSequence*motif->length*2];
	char *mask = new char[longestSequence];
	for (int i=0;i<longestSequence;i++) {
		mask[i] = 0;
	}
	int mlen = motif->length;
	int offset = 0;
	MotifScanTables tables;
	tables.build(motif);
	double totalTsites = 0;
	double totalBsites = 0;
	pthread_mutex_t notused;
	pthread_mutex_init(&notused,NULL);
	pthread_mutex_unlock(&notused);
	for (int i=0;i<numSeqs;i++) {

		int nsites = seqs[i]->motifScan(motif,strand,offset,KNOWN_FIND_ALL,score,
								sitePositions,siteStrands,NULL,notused,&tables);
		motif->seqStats[i].seq = seqs[i];
		motif->seqStats[i].score = score;
		motif->seqStats[i].setSites(sitePositions, siteStrands,nsites);

		numPositions2Mask = 0;
		int slen = seqs[i]->length;
		double targetValue = 0;
		double backgroundValue = 0;
		if (seqs[i]->g) {
			targetValue = -1*seqs[i]->w;
			totalTsites += seqs[i]->w*nsites;
		} else {
			backgroundValue = -1*seqs[i]->w;
			totalBsites += seqs[i]->w*nsites;
		}

		for (int j=0;j<nsites;j++) {
			int p = sitePositions[j];
			//fprintf(stderr, "%s site %d of %d. %d, strand=%d\n", seqs[i]->name, j, nsites,p,siteStrands[j]);
			if (siteStrands[j] == MOTIF_STRAND_POS) {
				for (int k=-mlen+1;k<mlen;k++) {
					int pos = p+k;
					if (pos < 0) continue;
					if (pos > slen-mlen) break;
					if (mask[pos] != 0) continue; // don't double delete
					char* s = &(seqs[i]->s[pos]);
					char tmp = seqs[i]->s[pos+mlen];
					seqs[i]->s[pos+mlen] = '\0';
					//fprintf(stderr, "%d\t%s\n", k,s);
					int ok = checkSequence(s);
					if (ok) {
						oligos->adjustOligoInstances(s, targetValue, backgroundValue);
					}
					seqs[i]->s[pos+mlen] = tmp;

					mask[pos] = 1;
					positions2Mask[numPositions2Mask++]=pos;
				}
			} else {
				for (int k=-mlen+1;k<mlen;k++) {
					int pos = p-mlen+1+k;
					if (pos < 0) continue;
					if (pos > slen-mlen) break;
					if (mask[pos] != 0) continue; // don't double delete

					char* s = &(seqs[i]->s[pos]);
					char tmp = seqs[i]->s[pos+mlen];
					seqs[i]->s[pos+mlen] = '\0';
					//fprintf(stderr, "%d\t%s\n", k,s);
					int ok = checkSequence(s);
					if (ok) {
						oligos->adjustOligoInstances(s, targetValue, backgroundValue);
					}
					seqs[i]->s[pos+mlen] = tmp;
				
					mask[pos] = 1;
					positions2Mask[numPositions2Mask++]=pos;
				}
			}
		}
		//now literally mask the sequence
		for (int j=0;j<nsites;j++) {
			int p = sitePositions[j];
			if (siteStrands[j] == MOTIF_STRAND_POS) {
				for (int k=0;k<mlen;k++) seqs[i]->s[p+k] = 'N';
			} else {
				for (int k=0;k<mlen;k++) seqs[i]->s[p-k] = 'N';
			}
		}
		for (int j=0;j<numPositions2Mask;j++) {
			int p  = positions2Mask[j];
			mask[p] = 0;
		}
	}
	motif->calculateDistributionStats();

	delete []mask;
	delete []sitePositions;
	delete []siteStrands;
	delete []positions2Mask;
}
void SequenceArray::optimizeMotifThreshold(Motif* m, char strand) {
	m->logxform();
	m->numTarget = 0.0;
	m->numBackground = 0.0;
	m->initializeSeqStats(numSeqs);
	int *sitePositions = new int[longestSequence];
	char *siteStrands = new char[longestSequence];
	int mode = KNOWN_OPTIMIZE_THREASHOLD;
	int offset = 0;
	MotifScanTables tables;
	tables.build(m);
	pthread_mutex_t notused;
	pthread_mutex_init(&notused,NULL);
	pthread_mutex_unlock(&notused);
	for (int i=0;i<numSeqs;i++) {
		double score = 0;
		//int nsites = 0;
		(void)seqs[i]->motifScan(m,strand,offset,mode,score,sitePositions,
											siteStrands,NULL,notused,&tables);
		if (mode == KNOWN_OPTIMIZE_THREASHOLD || mode==KNOWN_FIND_ALL) {
			m->seqStats[i].seq = seqs[i];
			m->seqStats[i].score = score;
		}
	}
	m->optimizeThreshold(scorer,maxBackgroundPercent);
	delete []sitePositions;
	delete []siteStrands;
}

void SequenceArray::motifScan(Motif* m, int strand,int offset, int mode,FILE* fp, pthread_mutex_t fileLock) {
	m->logxform();
	// find, width 6 and more, with the 6-mer index built: buildSix's blocks
	// replace the 5-mer ones
	int pivotP = 0;
	// known's ZOOPS count: only whether a sequence has a window at the
	// threshold matters, so the same index serves it
	const int zoops = mode == KNOWN_JUST_SCORE_ZOOPS && fp == NULL;
	if ((mode == KNOWN_FIND_ALL && fp != NULL) || zoops) {
		if (chunks7 != NULL && m->length >= 7) pivotP = 7;
		else if (chunks6 != NULL && m->length >= 6) pivotP = 6;
	}
	const int trySix = pivotP > 0;
	const KmerChunk* chunksP = (pivotP == 7) ? chunks7 : chunks6;
	MotifScanTables tables;
	tables.build(m, !trySix);
	const int six = trySix && tables.buildSix(m, pivotP, (pivotP == 7) ? count7 : count6);
	m->numTarget = 0.0;
	m->numBackground = 0.0;
	int *sitePositions = new int[longestSequence];
	char *siteStrands = new char[longestSequence];
	const int stats = mode == KNOWN_OPTIMIZE_THREASHOLD || (mode == KNOWN_FIND_ALL && !seqStatsUnread);
	if (stats) {
		m->initializeSeqStats(numSeqs);
	}
	// find with the prepared index visits only windows that can survive the
	// bound, chunk by chunk; sequences are still handled in order.
	const int indexed = six || (mode == KNOWN_FIND_ALL && fp != NULL && chunks != NULL && index5
			&& tables.indexable);
	static thread_local long long* keys = NULL;
	static thread_local int keySize = 0;
	static thread_local int* fcand = NULL;
	static thread_local int fcandSize = 0;
	static thread_local int* rcand = NULL;
	static thread_local int rcandSize = 0;
	static thread_local std::vector<int> foffv;
	static thread_local std::vector<int> roffv;
	if ((int)foffv.size() < numSeqs+1) {
		foffv.resize(numSeqs+1);
		roffv.resize(numSeqs+1);
	}
	int* foff = foffv.data();
	int* roff = roffv.data();
	int chunk = -1;
	int chunkEnd = 0;
	const double cut = m->threshold - 1e-6;
	//fprintf(stderr, "numSeqs = %d, offset = %d\n", numSeqs,offset);
	if (six && !stats) {
		// Plain find: nothing is printed or kept for a sequence without
		// candidates, so only those with some are visited, in sequence order,
		// and numTarget/numBackground add up in the same order.
		static thread_local long long* fkeys = NULL;
		static thread_local long long* rkeys = NULL;
		static thread_local int* fpos = NULL;
		static thread_local int* rpos = NULL;
		static thread_local int* fseq = NULL;
		static thread_local int* rseq = NULL;
		static thread_local int sixSize = 0;
		if (sixSize < maxChunk6Entries+1) {
			delete []fkeys;
			delete []rkeys;
			delete []fpos;
			delete []rpos;
			delete []fseq;
			delete []rseq;
			sixSize = maxChunk6Entries+1;
			fkeys = new long long[sixSize];
			rkeys = new long long[sixSize];
			fpos = new int[sixSize];
			rpos = new int[sixSize];
			fseq = new int[sixSize];
			rseq = new int[sixSize];
		}
		for (int ch=0;ch<numChunks;ch++) {
			const KmerChunk* c = &(chunksP[ch]);
			int nf = 0;
			int nr = 0;
			if (strand == MOTIF_STRAND_POS || strand==MOTIF_STRAND_BOTH) {
				nf = sixKeys(seqs,c,&tables,0,cut,fkeys,1);
			}
			if (strand == MOTIF_STRAND_BOTH || strand==MOTIF_STRAND_NEG) {
				nr = sixKeys(seqs,c,&tables,1,cut,rkeys,1);
			}
			if (nf == 0 && nr == 0) continue;
			// window positions -> (sequence, position), ascending
			nf = resolveAts(c, 0, m->length, (unsigned int*)fkeys, nf, fseq, fpos);
			nr = resolveAts(c, 1, m->length, (unsigned int*)rkeys, nr, rseq, rpos);
			if (nf == 0 && nr == 0) continue;
			if (zoops) {
				// a sequence counts when any of its candidates scores at the
				// threshold; summed exactly as Sequence::motifScan sums
				const int mlen = m->length;
				const double* fwd = tables.fwd;
				const double* rev = tables.rev;
				const double threshold = m->threshold;
				int a = 0;
				int b = 0;
				while (a < nf || b < nr) {
					const int qa = (a < nf) ? fseq[a] : INT_MAX;
					const int qb = (b < nr) ? rseq[b] : INT_MAX;
					const int q = (qa < qb) ? qa : qb;
					const Sequence* sq = seqs[c->first+q];
					int hit = 0;
					for (;a<nf && fseq[a]==q;a++) {
						if (hit) continue;
						const unsigned char* cc = sq->codes + fpos[a];
						double score = 0.0;
						for (int j=0;j<mlen;j++) score += fwd[j*5+cc[j]];
						if (score != -INFINITY && score >= threshold) hit = 1;
					}
					for (;b<nr && rseq[b]==q;b++) {
						if (hit) continue;
						const unsigned char* cc = sq->codes + rpos[b];
						double score = 0.0;
						for (int j=0;j<mlen;j++) score += rev[j*5+cc[-j]];
						if (score != -INFINITY && score >= threshold) hit = 1;
					}
					if (hit) {
						if (sq->g == 1) {
							m->numTarget += sq->w;
						} else {
							m->numBackground += sq->w;
						}
					}
				}
				continue;
			}
			// the sequences with candidates, in order: (q, first forward, first reverse)
			static thread_local std::vector<int> gq, ga, gb;
			gq.clear();
			ga.clear();
			gb.clear();
			{
				int a = 0;
				int b = 0;
				while (a < nf || b < nr) {
					const int qa = (a < nf) ? fseq[a] : INT_MAX;
					const int qb = (b < nr) ? rseq[b] : INT_MAX;
					const int q = (qa < qb) ? qa : qb;
					gq.push_back(q);
					ga.push_back(a);
					gb.push_back(b);
					while (a < nf && fseq[a] == q) a++;
					while (b < nr && rseq[b] == q) b++;
				}
				ga.push_back(nf);
				gb.push_back(nr);
			}
			const int ng = (int)gq.size();
			const int mlen = m->length;
			for (int g=0;g<ng;g++) {
				// prefetch: the Sequence eight groups ahead, its data four ahead:
				// the span of the group's windows, from the first one's
				// start to the last one's end, in codes and text, and the name
				// (which may cross a line)
				if (g+8 < ng) __builtin_prefetch(seqs[c->first+gq[g+8]]);
				if (g+4 < ng) {
					const int gg = g+4;
					const Sequence* p = seqs[c->first+gq[gg]];
					int lo = (ga[gg] < ga[gg+1]) ? fpos[ga[gg]] : rpos[gb[gg]] - mlen + 1;
					int hi = (gb[gg] < gb[gg+1]) ? rpos[gb[gg+1]-1] : fpos[ga[gg+1]-1] + mlen - 1;
					if (lo < 0) lo = 0;
					if (hi < lo) hi = lo;
					__builtin_prefetch(p->codes + lo);
					__builtin_prefetch(p->codes + hi);
					__builtin_prefetch(p->s + lo);
					__builtin_prefetch(p->s + hi);
					__builtin_prefetch(p->name);
					__builtin_prefetch(p->name + 32);
				}
				Sequence* sq = seqs[c->first+gq[g]];
				double score = 0;
				int nsites = sq->scanCandidates(m,strand,offset,score,sitePositions,siteStrands,fp,
							&tables, fpos+ga[g], ga[g+1]-ga[g], rpos+gb[g], gb[g+1]-gb[g]);
				if (nsites > 0) {
					if (sq->g == 1) {
						m->numTarget += sq->w;
					} else {
						m->numBackground += sq->w;
					}
				}
			}
		}
	} else for (int i=0;i<numSeqs;i++) {
		double score = 0;
		int nsites = 0;
		if (indexed) {
			if (i >= chunkEnd) {
				chunk++;
				const KmerChunk* c = six ? &(chunksP[chunk]) : &(chunks[chunk]);
				chunkEnd = c->first + c->count;
				for (int q=0;q<=c->count;q++) {
					foff[q] = 0;
					roff[q] = 0;
				}
				if (strand == MOTIF_STRAND_POS || strand==MOTIF_STRAND_BOTH) {
					if (six) {
						collectCandidates6(seqs,c,&tables,0,cut,keys,keySize,fcand,fcandSize,foff);
					} else {
						collectCandidates(seqs,c,&tables,0,cut,keys,keySize,fcand,fcandSize,foff);
					}
				}
				if (strand == MOTIF_STRAND_BOTH || strand==MOTIF_STRAND_NEG) {
					if (six) {
						collectCandidates6(seqs,c,&tables,1,cut,keys,keySize,rcand,rcandSize,roff);
					} else {
						collectCandidates(seqs,c,&tables,1,cut,keys,keySize,rcand,rcandSize,roff);
					}
				}
			}
			int q = i - chunks[chunk].first;
			nsites = seqs[i]->scanCandidates(m,strand,offset,score,sitePositions,siteStrands,fp,
						&tables, fcand+foff[q], foff[q+1]-foff[q], rcand+roff[q], roff[q+1]-roff[q]);
		} else {
			nsites = seqs[i]->motifScan(m,strand,offset,mode,score,sitePositions,
													siteStrands,fp,fileLock,&tables);
		}
		//fprintf(stderr, "%d\n", nsites);
		if (nsites > 0) {
			if (seqs[i]->g == 1) {
				m->numTarget += seqs[i]->w;
			} else {
				m->numBackground += seqs[i]->w;
			}
		}
		if (stats) {
			m->seqStats[i].seq = seqs[i];
			m->seqStats[i].score = score;
			m->seqStats[i].setSites(sitePositions, siteStrands,nsites);
		}
		if (mode == KNOWN_MSCORE_MOTIFS) {
			char curStrand = '+';
			if (siteStrands[0] == MOTIF_STRAND_NEG) {
				curStrand = '-';
			}
			scanOutSite(fp,seqs[i]->name,sitePositions[0],m->consensus,-1,m->name,curStrand,score);
		}
	}
	scanOutFlush(fp);
	if (numTarget > 0.0) m->freqTarget = m->numTarget / numTarget;
	if (numBackground > 0.0) m->freqBackground = m->numBackground / numBackground;

	if (stats) {
		m->calculateDistributionStats();
		if (mode == KNOWN_OPTIMIZE_THREASHOLD) {
			//m->deleteSites();
		}
	}

	delete []sitePositions;
	delete []siteStrands;
}
void Motif::calculateDistributionStats() {
	long double avgTarget= 0.0;
	double nTarget = 0;
	long double avgBack= 0.0;
	double nBack = 0;
	double nPos = 0.0;
	double nNeg = 0.0;
	double totalSites = 0.0;
	double totalSeqWithSites = 0.0;

	if (seqStats == NULL || numSeqs < 1) return;

	for (int i=0;i<numSeqs;i++) {
		double w = seqStats[i].seq->w;
		int g = seqStats[i].seq->g;
		if (seqStats[i].numSites < 1) continue;

		double uniqSites = 0;

		for (int j=0;j<seqStats[i].numSites;j++) {
			char d = seqStats[i].strands[j];
			int p = seqStats[i].pos[j];
			if (d == MOTIF_STRAND_POS) {
				p += length/2;
			} else {
				p -= length/2;
			}
			if (g == 1) {
				avgTarget += w*p;
				nTarget += w;
				if (d == MOTIF_STRAND_POS) nPos += w;
				else nNeg += w;
				int good = 1;
				for (int k=j+1;k<seqStats[i].numSites;k++) {
					char d2 = seqStats[i].strands[k];
					int p2 = seqStats[i].pos[k];
					if (d2 == MOTIF_STRAND_POS) {
						p2 += length/2;
					} else {
						p2 -= length/2;
					}
					if (p-p2>length/-2 && p-p2<length/2) {
						good = 0;
						break;
					}
				}
				if (good) uniqSites+=w;
			} else {
				avgBack += w*p;
				nBack += w;
			}
		}
		if (g == 1) {
			totalSeqWithSites += w;
			totalSites += uniqSites;
		}

	}

	if (nTarget > 0.0) avgTarget /= nTarget;
	if (nBack > 0.0) avgBack /= nBack;
	double stdTarget=0.0;
	double stdBack = 0.0;
	nTarget = 0.0;
	nBack = 0.0;

	for (int i=0;i<numSeqs;i++) {
		double w = seqStats[i].seq->w;
		int g = seqStats[i].seq->g;
		if (seqStats[i].numSites < 1) continue;

		for (int j=0;j<seqStats[i].numSites;j++) {
			char d = seqStats[i].strands[j];
			int p = seqStats[i].pos[j];
			if (d == MOTIF_STRAND_POS) {
				p += length/2;
			} else {
				p -= length/2;
			}
			if (g == 1) {
				double v = w*(p-avgTarget);
				stdTarget += v*v;
				nTarget += w;
			} else {
				double v = w*(p-avgBack);
				stdBack += v*v;
				nBack += w;
			}
		}
	}
	if (nTarget > 0.0) stdTarget /= nTarget;
	if (nBack > 0.0) stdBack /= nBack;
	stdTarget = sqrt(stdTarget);
	stdBack = sqrt(stdBack);

	statAvgTarget = avgTarget;
	statStdTarget = stdTarget;
	statAvgBackground = avgBack;
	statStdBackground = stdBack;

	statStrandLogRatio = 0.0;
	if (nNeg < 0.00001 && nPos > 0.00001) statStrandLogRatio = 10.0;
	else if (nNeg > 0.00001 && nPos < 0.00001) statStrandLogRatio = -10.0;
	else if (nNeg < 0.00001 && nPos < 0.00001) statStrandLogRatio = 0.0;
	else statStrandLogRatio = log(nPos/nNeg)/log(2);

	siteMultiplicity = 1.0;
	if (totalSeqWithSites > 0.01) siteMultiplicity = totalSites/totalSeqWithSites;

}


void SequenceArray::normalizeSequence(int nlen,int strand,char* normFile,int neutralFlag) {
	if (nlen < 1) return;

	double lastError = 1e100;
	double next = 1.0;
	double cur = 0.0;
	
	for (int i=0;i<maxAutoNormIters;i++) {
		int printFlag = 0;
		if (i==0) {
			printFlag = 1;
		}
		double curError = normalizeSequenceIteration(nlen,strand,normFile,printFlag,neutralFlag);
		//double improvementRate = (lastError-curError)/lastError;
		//fprintf(stderr, "Round %d:\tLastError: %le\tCurrentError:%le\tImprovment:%le\n", i+1,lastError,curError,improvementRate);
		if (curError >= lastError) break;
		//if (improvementRate < MIN_OLIGO_NORMALIZATION_IMPROVEMENT) break;
		lastError = curError;
		if (printFlag) {
			fprintf(stderr, "\tConverging on autonormalization solution:\n\t");
		} else {
			cur += 80.0/(double)maxAutoNormIters;
			if (cur > next) {
				fprintf(stderr, ".");
				next += 1.0;
			}
		}
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "\tFinal normalization:");
	(void)normalizeSequenceIteration(nlen,strand,NULL,1,neutralFlag);
}

double SequenceArray::normalizeSequenceIteration(int nlen,int strand,char* normFile,int printFlag,int neutralFlag) {

	FILE* fp = NULL;
	if (normFile != NULL && printFlag) {
		fp = fopen(normFile, "w");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for writing!!!\n",normFile);
		} else {
			fprintf(fp, "Oligo\tTargetCounts\tBackgroundCounts\tNormalizationFactor\n");
		}
	}
	char* alphaIndex = new char[256];
	alphaIndex[65]=0; //alphaIndex['A']=0;
	alphaIndex[67]=1; //alphaIndex['C']=1;
	alphaIndex[71]=2; //alphaIndex['G']=2;
	alphaIndex[84]=3; //alphaIndex['T']=3;
	char* alpha = new char[256];
	alpha[0]='A';
	alpha[1]='C';
	alpha[2]='G';
	alpha[3]='T';
	int* totalGood = new int[numSeqs];
	double* neutralFreq = new double[nlen+1];
	neutralFreq[0] = 1.0;
	for (int i=0;i<nlen+1;i++) {
		neutralFreq[i]=pow(4.0,-1.0*i);
	}

	double error = 0.0;

	char* mer = new char[1000];
	for (int curLen=1;curLen<=nlen;curLen++) {
		int total = (int) pow(4.0,(double)curLen);
		if (printFlag) fprintf(stderr, "\tAutonormalization: %d-mers (%d total)\n",curLen,total);
		double* targetLevels = new double[total];
		double* backgroundLevels = new double[total];
		double* normFactors = new double[total];
		double* normFactorsT = new double[total];
		int* pow4 = new int[curLen+1];
		pow4[0]=1;
		for (int i=1;i<=curLen;i++) pow4[i] = 4*pow4[i-1];
		for (int i=0;i<total;i++) {
			targetLevels[i]=0.0;
			backgroundLevels[i]=0.0;
			normFactors[i] = 1.0;
			normFactorsT[i] = 1.0;
		}

		for (int i=0;i<numSeqs;i++) {
			char* s = seqs[i]->s;
			int goodOligos = 0;
			int nbins = seqs[i]->length-curLen;
			for (int j=0;j<nbins;j++) {
				int bad = 0;
				for (int k=0;k<curLen;k++) {
					if (s[j+k] == 'N') {
						bad = 1;
						break;
					}
				}
				if (bad) continue;
				goodOligos++;
			}
			totalGood[i] = goodOligos;
			double gOligos = (double)totalGood[i];
			for (int j=0;j<nbins;j++) {
				int index = 0;
				int bad = 0;
				for (int k=0;k<curLen;k++) {
					if (s[j+k] == 'N') {
						bad = 1;
						break;
					}
					index += alphaIndex[(int)s[j+k]]*pow4[k];
				}
				if (bad) continue;
				if (seqs[i]->g == 1) {
					targetLevels[index] += seqs[i]->w/gOligos;
				} else {
					backgroundLevels[index] += seqs[i]->w/gOligos;
				}
			}
		}

		double totalTarget = 0.0;
		double totalBackground = 0.0;
		for (int i=0;i<total;i++) {
			totalTarget += targetLevels[i];
			totalBackground += backgroundLevels[i];
		}
		double minimumTargetLevels = 0.5/totalTarget;
		double minimumBackgroundLevels = 0.5/totalBackground;
		for (int i=0;i<total;i++) {
			for (int j=0;j<curLen;j++) {
				int remainder = i % pow4[j+1];
				int asdf = (int)(remainder/pow4[j]);
				mer[j] = alpha[asdf];
			}
			mer[curLen]='\0';

			double tLevel = targetLevels[i];
			double bLevel = backgroundLevels[i];
			if (strand == MOTIF_STRAND_BOTH) {
				int rindex = 0;
				for (int j=curLen-1;j>=0;j--) {
					rindex += (3-alphaIndex[(int)mer[j]])*pow4[curLen-j-1];
				}
				tLevel = (tLevel+targetLevels[rindex])/2.0;
				bLevel = (bLevel+backgroundLevels[rindex])/2.0;
			}

			if (bLevel < minimumBackgroundLevels) bLevel = minimumBackgroundLevels;
			if (tLevel < minimumTargetLevels) tLevel = minimumTargetLevels;
			
			if (neutralFlag == 0) {
				double factor = totalBackground*tLevel/totalTarget/bLevel;
				if (printFlag) {
					if (curLen < 3) {
						fprintf(stderr, "\t\t%s\t%.2lf%%\t%.2lf%%\t%.3lf\n", mer,
											100.0*tLevel/totalTarget,
											100.0*bLevel/totalBackground,factor);
					}
					if (fp != NULL) {
						//fprintf(fp, "%s\t%.1lf\t%.1lf\t%.3lf\n", mer,tLevel,bLevel,factor);
						fprintf(fp, "%s\t%.2lf%%\t%.2lf%%\t%.3lf\n", mer, 
											100.0*tLevel/totalTarget,
											100.0*bLevel/totalBackground,factor);
					}
				}
				normFactors[i] = factor;
				error += (factor-1.0)*(factor-1.0)/(double)total;
			} else {
				double factorNeutralT = totalTarget/tLevel*neutralFreq[curLen];
				double factorNeutralB = totalBackground/bLevel*neutralFreq[curLen];
				if (printFlag) {
					if (curLen < 3) {
						fprintf(stderr, "\t\t%s\t%.2lf%%\t%.2lf%%\t%.3lf,%.3lf\n", mer,
											100.0*tLevel/totalTarget,
											100.0*bLevel/totalBackground,factorNeutralT,factorNeutralB);
					}
					if (fp != NULL) {
						//fprintf(fp, "%s\t%.1lf\t%.1lf\t%.3lf\n", mer,tLevel,bLevel,factor);
						fprintf(fp, "%s\t%.2lf%%\t%.2lf%%\t%.3lf,%.3lf\n", mer, 
											100.0*tLevel/totalTarget,
											100.0*bLevel/totalBackground,factorNeutralT,factorNeutralB);
					}
				}
				normFactorsT[i] = factorNeutralT;
				normFactors[i] = factorNeutralB;
				error += (factorNeutralT-1.0)*(factorNeutralT-1.0)/(double)total;
				error += (factorNeutralB-1.0)*(factorNeutralB-1.0)/(double)total;
			}
		}

		totalBackground = 0.0;
		totalTarget = 0.0;
		for (int i=0;i<numSeqs;i++) {
			char* s = seqs[i]->s;
			double scoreUnit = (double)seqs[i]->w;
			if (neutralFlag == 0 && seqs[i]->g == 1) {
				totalTarget += scoreUnit;
				continue;
			}
			double newScore = 0;
			int nbins = seqs[i]->length-curLen;
			double gOligos = (double)totalGood[i];
			for (int j=0;j<nbins;j++) {
				int index = 0;
				int bad = 0;
				for (int k=0;k<curLen;k++) {
					if (s[j+k] == 'N') {
						bad=1;
						break;
					}
					index += alphaIndex[(int)s[j+k]]*pow4[k];
				}
				if (bad) continue;
				if (seqs[i]->g == 1) {
					newScore += normFactorsT[index];
				} else {
					newScore += normFactors[index];
				}
			}
			if (gOligos > 0.5) {
				newScore /= gOligos;
			}
			double newWeight = newScore*scoreUnit;
			if (newWeight < HOMER_MINIMUM_SEQ_WEIGHT) newWeight = HOMER_MINIMUM_SEQ_WEIGHT;
			if (newWeight > 1.0/HOMER_MINIMUM_SEQ_WEIGHT) newWeight = 1.0/HOMER_MINIMUM_SEQ_WEIGHT;

			double penalty = newWeight;
			if (newWeight < 1) penalty = 1/newWeight;
			penalty*=penalty;
		
			double delta = newWeight-seqs[i]->w;
			double newWeight1 = seqs[i]->w+delta;
			//double newWeight2 = seqs[i]->w+delta;
			if (penalty > 1 && ((delta > 0 && newWeight > 1) || (delta < 0 && newWeight < 1))) {
				newWeight1 = seqs[i]->w+delta/penalty;
			} else {
			}
			//fprintf(stdout, "1\t%lf\t%lf\n",newWeight2,seqs[i]->w);
			//fprintf(stdout, "%lf\n",penalty);
			//fprintf(stdout, "%lf\t%lf\n",newWeight2,seqs[i]->w);

			//if (newWeight1 > 4) newWeight1 = 4;
			//if (newWeight1 < 0.25) newWeight1 = 0.25;
			
			seqs[i]->w = newWeight1;
			if (seqs[i]->g == 1) {
				totalTarget += seqs[i]->w;
			} else {
				totalBackground += seqs[i]->w;
			}
		}
		numTarget = totalTarget;
		numBackground = totalBackground;
		delete []targetLevels;
		delete []backgroundLevels;
		delete []normFactors;
		delete []pow4;
	}
	if (fp != NULL) {
		fprintf(stderr, "\tNormalization weights can be found in file: %s\n", normFile);
		fclose(fp);
	}
	delete []mer;
	delete []alphaIndex;
	delete []alpha;
	delete []totalGood;
	delete []neutralFreq;

	return error;

}

// ACGTN stay, acgt become ACGT, every other byte becomes N
struct CleanSequenceTable {
	char t[256];
	CleanSequenceTable() {
		for (int i=0;i<256;i++) t[i] = 'N';
		t[(unsigned char)'A'] = 'A';
		t[(unsigned char)'C'] = 'C';
		t[(unsigned char)'G'] = 'G';
		t[(unsigned char)'T'] = 'T';
		t[(unsigned char)'a'] = 'A';
		t[(unsigned char)'c'] = 'C';
		t[(unsigned char)'g'] = 'G';
		t[(unsigned char)'t'] = 'T';
	}
};
static const CleanSequenceTable cleanSequenceTable;
static void copyCleaned(char* dst, const char* src, size_t n) {
	for (size_t i=0;i<n;i++) dst[i] = cleanSequenceTable.t[(unsigned char)src[i]];
}
void cleanSequence(char* seq) {
	if (seq == NULL) return;
	int len = strlen(seq);
	for (int i=0;i<len;i++) {
		seq[i] = cleanSequenceTable.t[(unsigned char)seq[i]];
	}
}

//=====================================================================
KnownMotifs::KnownMotifs() {
	directFasta = 0;
	seqArray = new SequenceArray();
	totalBp = 0;
	normLength = 0;
	alphaSize = DNA_ALPHA;
	strand = MOTIF_STRAND_BOTH;
	maxCacheSize = MAX_STAT_CACHE_SIZE*1000000/sizeof(double);
	offset = MOTIF_OFFSET_NULL;

	numTargets = 0.0;
	numBackground= 0.0;

	motifs = NULL;
	numMotifs = 0;
	maxAutoNormIters = MAX_OLIGO_NORMALIZATION_ITERATIONS;

	numCPUs = 1;
	nextMeter = 0.0;
	incMeter = 0.0;
	
	seqFile = NULL;
	groupFile = NULL;
	inputFASTAfile = NULL;
	bgFASTAfile = NULL;
	fastaFlag = 0;
	knownFile = NULL;
	normFile = NULL;
	sitefp = NULL;
	normfp = NULL;
	motifFilename = NULL;
	neutralFlag = 0;

	maxBackgroundPercent = MOTIF_MAX_BACKGROUND_PERCENT;
	siteReduceThreshold = -1.0;

	//scoringMethod = MOTIF_SCORING_HYPERGEOMETRIC;
	scoringMethod = MOTIF_SCORING_BINOMIAL;
	scorer = NULL;

	mscoreFlag = 0;
	threasholdOptimizationFlag = KNOWN_JUST_SCORE_ZOOPS;
	zoopsFlag = SCORER_STANDARD_ZOOPS;
	
	knownFlag = 1;	
	findFlag = 0;
	maskFlag = 0;
	normOnlyFlag = 0;
}

KnownMotifs::~KnownMotifs() {
	delete seqArray;
	delete scorer;
}
void KnownMotifs::checkParameters() {
	if (numCPUs < 1) numCPUs = 1;
	if (seqFile == NULL){
		if (inputFASTAfile != NULL) {
			fastaFlag=1;
		} else {
			fprintf(stderr, "!!! Need a sequence file for homer to do anything !!!\n");
			exit(0);
		}
	}
	if (fastaFlag) {
		if (findFlag == 0 && maskFlag == 0 && bgFASTAfile == NULL){
			fprintf(stderr, "!!! Need a background FASTA file for homer to do anything !!!\n");
			exit(0);
		}
	} else {
		if (findFlag == 0 && maskFlag == 0 && groupFile == NULL){
			fprintf(stderr, "!!! Need a group file for homer to do anything !!!\n");
			exit(0);
		}
	}
	if (normOnlyFlag == 0 && knownFile == NULL ){
		fprintf(stderr, "!!! Need a file of known motifs for homer to do anything !!!\n");
		exit(0);
	}
}
void KnownMotifs::parseFasta() {
	seqFile = new char[10000];
	groupFile = new char[10000];
	//srand(time(NULL));
	struct timeval tv;
	gettimeofday(&tv,NULL);
	unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
	srand(time_in_micros);
	//fprintf(stderr, "seed=%ld\n", time_in_micros);

	int r = rand();
	sprintf(seqFile, "0%d.seq", r);
	sprintf(groupFile, "0%d.group", r);
	FILE* fptest = fopen(seqFile,"r");
	while (fptest != NULL) {
		fclose(fptest);
		r = rand();
		sprintf(seqFile, "0%d.seq", r);
		sprintf(groupFile, "0%d.group", r);
		fptest = fopen(seqFile,"r");
	}

	// loadSequence reads the FASTA in memory; nothing is written here.
	directFasta = 1;
}


void KnownMotifs::loadMotifs() {
	motifs = Motif::readMotifFile(knownFile,numMotifs);
}

void KnownMotifs::loadSequence() {

	seqArray->maxAutoNormIters = maxAutoNormIters;
	if (directFasta) {
		seqArray->readFastaDirect(inputFASTAfile,bgFASTAfile,(int)numCPUs);
	} else {
		seqArray->readSequenceFiles(seqFile,groupFile);
	}
	if (normLength > 0 && !findFlag) {
		seqArray->normalizeSequence(normLength, strand, normFile,neutralFlag);
	}

	numTargets = seqArray->numTarget;
	numBackground = seqArray->numBackground;
	if (normOnlyFlag) {
		seqArray->printGroupFile(normfp);
	}
}


void KnownMotifs::getMotifEnrichment() {

	scorer = new EnrichmentScorer();
	if (!findFlag) {
		scorer->initializeEnrichmentScoring(scoringMethod,numTargets, numBackground,
												zoopsFlag,maxCacheSize);
	}

	if (maskFlag) {
		threasholdOptimizationFlag = KNOWN_MASK_MOTIFS;
	} else if (mscoreFlag) {
		threasholdOptimizationFlag = KNOWN_MSCORE_MOTIFS;
	} else if (siteReduceThreshold > 0) {
		threasholdOptimizationFlag = KNOWN_FIND_ALL;
	} else if (findFlag) {
		threasholdOptimizationFlag = KNOWN_FIND_ALL;
		fprintf(stderr, "\tFinding instances of %d motif(s)\n", numMotifs);
	} else {
		fprintf(stderr, "\tChecking enrichment of %d motif(s)\n", numMotifs);
	}

	fprintf(stderr, "\t|0%%                                    50%%                                  100%%|\n");
	fprintf(stderr, "\t");

	//parallel check
	pthread_mutex_init(&mutex, NULL);
	pthread_mutex_init(&fileLock, NULL);
	if (threasholdOptimizationFlag == KNOWN_FIND_ALL && sitefp != NULL) {
		// the 6-mer index serves widths 6 and more, the 5-mer one widths 4-5
		int need5 = 0;
		int need6 = 0;
		int need7 = 0;
		int maxWidth = 0;
		for (int k=0;k<numMotifs;k++) {
			if (motifs[k]->length >= 7) need7 = 1;
			else if (motifs[k]->length == 6) need6 = 1;
			else if (motifs[k]->length >= 4) need5 = 1;
			if (motifs[k]->length > maxWidth) maxWidth = motifs[k]->length;
		}
		seqArray->prepareKmers((int)numCPUs, need5, need6, maxWidth, need7);
	} else if (threasholdOptimizationFlag == KNOWN_JUST_SCORE_ZOOPS && sitefp == NULL) {
		// known: the 6/7-mer index finds, for motifs of width 6 and more,
		// the sequences with a window at the threshold
		int need6 = 0;
		int need7 = 0;
		int maxWidth = 0;
		for (int k=0;k<numMotifs;k++) {
			if (motifs[k]->length >= 7) need7 = 1;
			else if (motifs[k]->length == 6) need6 = 1;
			if (motifs[k]->length > maxWidth) maxWidth = motifs[k]->length;
		}
		if (need6 || need7) seqArray->prepareKmers((int)numCPUs, 0, need6, maxWidth, need7);
	}
	// The stats are printed by -mout (Motif::print) and compared by
	// -siteReduce (scoreOverlap); plain find prints neither.
	seqArray->seqStatsUnread = findFlag && threasholdOptimizationFlag == KNOWN_FIND_ALL
			&& siteReduceThreshold <= 0 && motifFilename == NULL;
	mutexIndex = 0;
	mutexTotal = numMotifs;
	incMeter = 1.0/((double)numMotifs);
	nextMeter = 0;

	if (threasholdOptimizationFlag == KNOWN_MASK_MOTIFS) {
		seqArray->maskMotifs(motifs, numMotifs, strand, offset, (int)numCPUs);
		for (int i=0;i<numMotifs;i++) {
			while (((double)i+1.0)/((double)numMotifs) >= nextMeter) {
				fprintf(stderr, "=");
				nextMeter += 1.0/80.0;
			}
			if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
				motifs[i]->freqTarget = scorer->getZoopsApproxTarget(motifs[i]->numTarget)/seqArray->numTarget;
				motifs[i]->freqBackground = scorer->getZoopsApproxBackground(motifs[i]->numBackground)
																			/seqArray->numBackground;
			}
		}
	} else {
	OrderedOutput* ordered = NULL;
	if (sitefp != NULL && (threasholdOptimizationFlag == KNOWN_FIND_ALL
			|| threasholdOptimizationFlag == KNOWN_MSCORE_MOTIFS)) {
		ordered = new OrderedOutput(sitefp, numMotifs);
		orderedOut = ordered;
		if (threasholdOptimizationFlag == KNOWN_FIND_ALL && seqArray->count7 != NULL) {
			double c4[4] = {0.0, 0.0, 0.0, 0.0};
			for (int x=0;x<16384;x++) c4[x >> 12] += seqArray->count7[x];
			const double total = c4[0]+c4[1]+c4[2]+c4[3];
			double comp[4];
			for (int b=0;b<4;b++) comp[b] = (total > 0.0) ? c4[b]/total : 0.25;
			std::vector<double> w(numMotifs, 0.0);
			for (int k=0;k<numMotifs;k++) w[k] = motifScanCost(motifs[k], comp, total);
			ordered->setWeights(w);
		}
	}

	ThreadArgs_scanMotifThread** args = new ThreadArgs_scanMotifThread*[numCPUs];

	pthread_t* thread = new pthread_t[numCPUs];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void* status;
	for (unsigned int i=0;i<numCPUs;i++) {
		args[i] = new ThreadArgs_scanMotifThread(this,i);
		(void)pthread_create(&(thread[i]), &attr, 
						&KnownMotifs_scanMotifThread, args[i]);
	}
	pthread_attr_destroy(&attr);
	for (unsigned int i=0;i<numCPUs;i++) {
		(void)pthread_join(thread[i], &status);
		delete args[i];
	}
	delete []thread;
	delete []args;
	orderedOut = NULL;
	delete ordered;
	}
	pthread_mutex_destroy(&mutex);
	pthread_mutex_destroy(&fileLock);

	fprintf(stderr, "\n");
}


void* KnownMotifs_scanMotifThread(void* threadArgs) {
	ThreadArgs_scanMotifThread* args = (ThreadArgs_scanMotifThread*)threadArgs;
	args->known->scanMotifThread(args);
	pthread_exit(threadArgs);
	return NULL;
}
ThreadArgs_scanMotifThread::ThreadArgs_scanMotifThread(KnownMotifs* km, int ncpu) {
	known = km;
	cpu = ncpu;
}

void KnownMotifs::scanMotifThread(ThreadArgs_scanMotifThread* args) {
	int i=-1;
	while (i<mutexTotal) {
		pthread_mutex_lock(&mutex);
		i=mutexIndex++;
		if (i < mutexTotal) {
			while (((double)i+1.0)/((double)numMotifs) >= nextMeter) {
				fprintf(stderr, "=");
				nextMeter += 1.0/80.0;
			}
		}
		pthread_mutex_unlock(&mutex);
		if (i >= mutexTotal) break;

		int k = i;
		if (orderedOut != NULL) k = orderedOut->take();
		Motif* motif = motifs[k];
		if (orderedOut != NULL) orderedMotif = k;
		seqArray->motifScan(motif,strand,offset,threasholdOptimizationFlag,sitefp,fileLock);
		if (orderedOut != NULL) {
			orderedOut->finish(k);
			orderedMotif = -1;
		}
		// Plain find reads a motif's per-sequence stats only inside its own
		// scan; free them so the next motif reuses the memory. -siteReduce
		// compares them across motifs afterwards.
		if (findFlag && threasholdOptimizationFlag == KNOWN_FIND_ALL && siteReduceThreshold <= 0) {
			motif->deleteSeqStats();
		}
		if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
			motif->freqTarget = scorer->getZoopsApproxTarget(motif->numTarget)/seqArray->numTarget;
			motif->freqBackground = scorer->getZoopsApproxBackground(motif->numBackground)
																			/seqArray->numBackground;
		}
		if (threasholdOptimizationFlag == KNOWN_OPTIMIZE_THREASHOLD) {
			motif->optimizeThreshold(scorer,maxBackgroundPercent);
		} else if (threasholdOptimizationFlag == KNOWN_JUST_SCORE_ZOOPS) {
			if (findFlag == 0 && numBackground > 0.5 && numTargets > 0.5) {
				motif->logp = scorer->scoreEnrichment(motif->numTarget, motif->numBackground);
			} else {
				motif->logp = 0.0;
			}
		}
	}

}



void KnownMotifs::printMotifEnrichment(FILE* fp) {
	if (fp == NULL) return;
	if (numMotifs < 1) {
		fprintf(stderr, "!!! No motifs to analyze !!!\n");
		return;
	}
	qsort(motifs, numMotifs, sizeof(Motif*),&(Motif::motifCmp));

	fprintf(fp, "Motif Name\tConsensus\tp-value\tLog p-value\tLog2 enrichment ratio");
	fprintf(fp, "\t# target sequences with Motif(total=%d)", (int) numTargets);
	fprintf(fp, "\t%% of target sequences with Motif");
	fprintf(fp, "\t# background sequences with Motif(total=%d)", (int) numBackground);
	fprintf(fp, "\t%% of background sequences with Motif");
	fprintf(fp, "\n");
	for (int i=0;i<numMotifs;i++) {
		fprintf(fp,"%s\t%s\t1e%d\t%.3lf",motifs[i]->name,motifs[i]->consensus,
					(int)(motifs[i]->logp/2.3026),motifs[i]->logp);
		double tfreq = motifs[i]->numTarget/numTargets;
		double bgfreq = motifs[i]->numBackground/numBackground;
		double ratio = tfreq;
		if (ratio < 1.0/numTargets) {
			ratio = 1.0/numTargets;
		}
		if (bgfreq < 1.0/numBackground) {
			ratio = ratio / (1.0/numBackground);
		} else {
			ratio = ratio / bgfreq;
		}
		ratio = log(ratio)/log(2.0);
		fprintf(fp, "\t%.3lf",ratio);
		fprintf(fp, "\t%d\t%.2lf%%",(int)motifs[i]->numTarget,tfreq*100);
		fprintf(fp, "\t%d\t%.2lf%%",(int)motifs[i]->numBackground,bgfreq*100);
		fprintf(fp, "\n");
	
	}
}
void KnownMotifs::printMotifs(FILE* fp) {
	if (numMotifs < 1) return;
	qsort(motifs, numMotifs, sizeof(Motif*), &(Motif::motifCmp));
	for (int i=0;i<numMotifs;i++) {
		motifs[i]->print(fp);
	}
}

void KnownMotifs::findOverlappingMotifs() {

	qsort(motifs, numMotifs, sizeof(Motif*), &(Motif::motifCmp));
	
	char* mask = new char[numMotifs];
	for (int i=0;i<numMotifs;i++) {
		mask[i] = 0;
	}
	for (int i=0;i<numMotifs;i++) {
		if (mask[i]) continue;
		for (int j=i+1;j<numMotifs;j++) {
			if (mask[j]) continue;
			double percentI = 0.0, percentJ=0.0;
			double score = Motif::scoreOverlap(motifs[i],motifs[j],percentI, percentJ);
			if (score > siteReduceThreshold) {
				mask[j] = 1;
				motifs[j]->similarPercent = percentJ;
				motifs[j]->similarPercentOfOther = percentI;
				motifs[j]->similarMotif = motifs[i];
			}
		}
	}
}

//=================================================================================

EnrichmentScorer::EnrichmentScorer() {
	numTargets = 0;
	numBackground = 0;
	zoopsFlag = 1;
	targetZoopsApprox = NULL;
	bgZoopsApprox = NULL;
	maxZoopsApprox = 0;
	enrichmentFunction = NULL;

	cache = NULL;
	maxCacheSize = MAX_STAT_CACHE_SIZE;
}
EnrichmentScorer::~EnrichmentScorer() {
	if (targetZoopsApprox != NULL) delete []targetZoopsApprox;
	if (bgZoopsApprox != NULL) delete []bgZoopsApprox;
	if (cache != NULL) {
		for (int i=0;i<cacheLength;i++) {
			double* row = cache[i].load(std::memory_order_relaxed);
			if (row != NULL) delete []row;
		}
		delete []cache;
	}
}

// A row is filled in private and published with a compare-and-swap; a thread
// that loses the race frees its copy and uses the winner's.
double* EnrichmentScorer::getCacheRow(unsigned int mtotal) {
	double* row = cache[mtotal].load(std::memory_order_acquire);
	if (row != NULL) return row;
	double* fresh = new double[mtotal+1];
	for (unsigned int i=0;i<mtotal+1;i++) fresh[i]=10.0;
	if (cache[mtotal].compare_exchange_strong(row, fresh, std::memory_order_acq_rel,
												std::memory_order_acquire)) {
		return fresh;
	}
	delete []fresh;
	return row;
}

double EnrichmentScorer::scoreEnrichmentHypergeometric(double pv, double nv) {
	unsigned int targetValue = (unsigned int)pv;
	unsigned int bgValue = (unsigned int)nv;
	if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
		if (targetValue >= maxZoopsApprox) targetValue = maxZoopsApprox-1;
		if (bgValue >= maxZoopsApprox) bgValue = maxZoopsApprox-1;
		targetValue = (unsigned int)targetZoopsApprox[targetValue];
		bgValue = (unsigned int)bgZoopsApprox[bgValue];
	}
	unsigned int mtotal = targetValue+bgValue;
	double logp = 0.0;

	if (mtotal < (unsigned int)cacheLength) {
		double* row = getCacheRow(mtotal);
		if (row[targetValue] > 1.0) {
			row[targetValue] = loghypergeoD((unsigned int)(numTargets+numBackground),
					(unsigned int)numTargets, mtotal, targetValue);
		}
		logp = row[targetValue];
	} else {
		logp = loghypergeoD((unsigned int)(numTargets+numBackground),
				(unsigned int)numTargets, mtotal, targetValue);
	}
	return logp;
}


double EnrichmentScorer::scoreEnrichmentZscore(double pv, double nv) {
	double targetValue = pv;
	double bgValue = nv;
	if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
		unsigned int pvInt = (unsigned int)pv;
		unsigned int nvInt = (unsigned int)nv;
		if (pvInt >= maxZoopsApprox) pvInt = maxZoopsApprox-1;
		if (nvInt >= maxZoopsApprox) nvInt = maxZoopsApprox-1;
		targetValue = targetZoopsApprox[pvInt];
		bgValue = bgZoopsApprox[nvInt];
	}
	if (bgValue < 1.0) bgValue = 1.0;

	double bgRatio = bgValue/numBackground;
	double std = sqrt(numBackground*bgRatio*(1-bgRatio));
	if (std < 1.0) std = 1.0;
	double expected = numTargets*bgRatio;
	return -1* (targetValue-expected)/std;

}
double EnrichmentScorer::scoreEnrichmentBinomial(double pv, double nv) {
	unsigned int targetValue = (unsigned int)pv;
	unsigned int bgValue = (unsigned int)nv;
	if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
		if (targetValue >= maxZoopsApprox) targetValue = maxZoopsApprox-1;
		if (bgValue >= maxZoopsApprox) bgValue = maxZoopsApprox-1;
		targetValue = (unsigned int)targetZoopsApprox[targetValue];
		bgValue = (unsigned int)bgZoopsApprox[bgValue];
	}
	unsigned mtotal = targetValue+bgValue;
	double ratio = ((double)bgValue)/numBackground;
	double logp = 0.0;
	if (mtotal < (unsigned int)cacheLength) {
		double* row = getCacheRow(mtotal);
		if (row[targetValue] > 1.0) {
			row[targetValue] = logbinomialD((unsigned int)(numTargets),
						(unsigned int)targetValue, ratio, (unsigned int)numBackground);
		} else {
			//fprintf(stderr, ".");
		}
		logp = row[targetValue];
	} else {
		logp = logbinomialD((unsigned int)(numTargets),
						(unsigned int)targetValue, ratio, (unsigned int)numBackground);
	}
	return logp;
}
void EnrichmentScorer::initializeEnrichmentScoring(int scoringMethod,double numT, 
									double numB,int zFlag, unsigned long long int newMaxCacheSize) {

	numTargets = (unsigned int)numT;
	numBackground = (unsigned int)numB;
	zoopsFlag = zFlag;
	maxCacheSize = newMaxCacheSize;
	cacheLength = (int)sqrt((double)maxCacheSize*2.0);
	fprintf(stderr, "\tCache length = %d\n", cacheLength);
	cache = new std::atomic<double*>[cacheLength];
	for (int i=0;i<cacheLength;i++) {
		cache[i].store(NULL, std::memory_order_relaxed);
	}

	if (scoringMethod == MOTIF_SCORING_HYPERGEOMETRIC) {
		enrichmentFunction = &EnrichmentScorer::scoreEnrichmentHypergeometric;
		fprintf(stderr, "\tUsing hypergeometric scoring\n");
	} else if (scoringMethod == MOTIF_SCORING_BINOMIAL) {
		enrichmentFunction = &EnrichmentScorer::scoreEnrichmentBinomial;
		fprintf(stderr, "\tUsing binomial scoring\n");
	} else if (scoringMethod == MOTIF_SCORING_ZSCORE) {
		enrichmentFunction = &EnrichmentScorer::scoreEnrichmentZscore;
		fprintf(stderr, "\tUsing z-score scoring\n");
	}

	if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
		if (targetZoopsApprox != NULL) delete []targetZoopsApprox;
		if (bgZoopsApprox != NULL) delete []bgZoopsApprox;
			
		maxZoopsApprox = (int) (5.0*(numT+numB));

		//fprintf(stderr, "MaxZoopsApprox=%d\n", maxZoopsApprox);
		targetZoopsApprox = new double[maxZoopsApprox];
		bgZoopsApprox = new double[maxZoopsApprox];
		targetZoopsApprox[0]=0.0;
		bgZoopsApprox[0]=0.0;
		targetZoopsApprox[1]=1.0;
		bgZoopsApprox[1]=1.0;
		for (unsigned int i=2;i<maxZoopsApprox;i++) {
			targetZoopsApprox[i]=targetZoopsApprox[i-1]+(numT-targetZoopsApprox[i-1])/numT;
			bgZoopsApprox[i]=bgZoopsApprox[i-1]+(numB-bgZoopsApprox[i-1])/numB;
		}
	}
}

double EnrichmentScorer::scoreEnrichment(double pv, double nv) {
	
	return (this->*enrichmentFunction)(pv, nv);
}
double EnrichmentScorer::getZoopsApproxTarget(double v) {
	if (v >= maxZoopsApprox) v = maxZoopsApprox-1;
	return targetZoopsApprox[(int)v];
}
double EnrichmentScorer::getZoopsApproxBackground(double v) {
	if (v >= maxZoopsApprox) v = maxZoopsApprox-1;
	return bgZoopsApprox[(int)v];
}




//=======================================================================================================


void programBackground(int argc, char** argv) {

	BgSequenceArray* targets = new BgSequenceArray();
	BgSequenceArray* bg = new BgSequenceArray();

	int repMaskFlag = 0;
	int strandFlag = MOTIF_STRAND_BOTH;
	int size = -1;
	int autoBinSizeFlag = 1;
	//int defaultSize = 200;
	long int numBgSeqs = 100000;
	long int numConsideredBgSeqs = 100000000;
	int numBins = 10;
	int numSubBins = 10;
	const char* outputFilePrefix = "out";
	char* inputTargetFastaFile = NULL;
	char* inputTargetBedFile = NULL;
	char* inputBackgroundFastaFile = NULL;
	char* inputBackgroundBedFile = NULL;
	char* inputBackgroundMaskBedFile = NULL;
	char* pscoreOutputFile = NULL;
	char* genomeFasta = NULL;
	double maxPercentNperSeq = 0.5;
	double decayRate = 0.75;
	long seed = -1;
	int bgOverlapFlag = 1;
	int targetOverlapFlag = 1;
	int kmer = 1;
	int modelFlag = 0;
	int maxKmerOptIterations = 10;
	int numkmerIters = 20;
	int overlapIteration = 5;
	int positionFlag = 0;
	int homerStyleOutputFlag = 0;
	int alphaLen = 0;
	int bgPosRegionFlag = 0;
	//double kmerChangeRate = 0.25;

	if (argc < 3) {
		printCMDbackground();
	}
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-h")==0) {
			printCMDbackground();
			//known->groupFile = argv[++i];
		} else if (strcmp(argv[i],"-i")==0) {
			inputTargetFastaFile = argv[++i];
		} else if (strcmp(argv[i],"-p")==0) {
			inputTargetBedFile = argv[++i];
		} else if (strcmp(argv[i],"-model")==0) {
			modelFlag = 1;
		} else if (strcmp(argv[i],"-excludeNs")==0) {
			alphaLen = 4;
		} else if (strcmp(argv[i],"-includeNs")==0) {
			alphaLen = 5;
		} else if (strcmp(argv[i],"-b")==0) {
			inputBackgroundFastaFile = argv[++i];
		} else if (strcmp(argv[i],"-bg")==0) {
			inputBackgroundBedFile = argv[++i];
			bgPosRegionFlag = 0;
		} else if (strcmp(argv[i],"-bgr")==0) {
			inputBackgroundBedFile = argv[++i];
			bgPosRegionFlag = 1;
		} else if (strcmp(argv[i],"-bgm")==0) {
			inputBackgroundMaskBedFile = argv[++i];
		} else if (strcmp(argv[i],"-g")==0) {
			genomeFasta = argv[++i];
		} else if (strcmp(argv[i],"-o")==0) {
			outputFilePrefix=argv[++i];
		} else if (strcmp(argv[i],"-pscore")==0) {
			pscoreOutputFile=argv[++i];
		} else if (strcmp(argv[i],"-size")==0) {
			sscanf(argv[++i],"%d", &(size));
		} else if (strcmp(argv[i],"-N")==0) {
			sscanf(argv[++i],"%ld", &(numBgSeqs));
		} else if (strcmp(argv[i],"-NN")==0) {
			sscanf(argv[++i],"%ld", &(numConsideredBgSeqs));
		} else if (strcmp(argv[i],"-pkmer")==0) {
			sscanf(argv[++i],"%d", &(kmer));
			positionFlag = 1;
		} else if (strcmp(argv[i],"-ikmer")==0) {
			sscanf(argv[++i],"%d", &(kmer));
			positionFlag = 0;
		} else if (strcmp(argv[i],"-maxFractionN")==0) {
			sscanf(argv[++i],"%lf", &(maxPercentNperSeq));
		} else if (strcmp(argv[i],"-mask")==0) {
			repMaskFlag = 1;
		} else if (strcmp(argv[i],"-seed")==0) {
			sscanf(argv[++i],"%ld", &(seed));
		} else if (strcmp(argv[i],"-nbins")==0) {
			sscanf(argv[++i],"%d", &(numBins));
			autoBinSizeFlag = 0;
		} else if (strcmp(argv[i],"-nsubBins")==0) {
			sscanf(argv[++i],"%d", &(numSubBins));
			autoBinSizeFlag = 0;
		} else if (strcmp(argv[i],"-maxIterations")==0) {
			sscanf(argv[++i],"%d", &maxKmerOptIterations);
		} else if (strcmp(argv[i],"-maxIterations")==0) {
			sscanf(argv[++i],"%d", &numkmerIters);
		} else if (strcmp(argv[i],"-overlapIteration")==0) {
			sscanf(argv[++i],"%d", &overlapIteration);
		} else if (strcmp(argv[i],"-gs")==0) {
			homerStyleOutputFlag = 1;
		} else if (strcmp(argv[i],"-decayRate")==0) {
			sscanf(argv[++i],"%lf", &decayRate);
		} else if (strcmp(argv[i],"-strand")==0) {
			strandFlag = MOTIF_STRAND_SEPARATE;
			/*
			i++;
			if (i < argc) {
				if (strcmp(argv[i],"+")==0) strandFlag = MOTIF_STRAND_POS;
				if (strcmp(argv[i],"-")==0) strandFlag = MOTIF_STRAND_NEG;
				if (strcmp(argv[i],".")==0) strandFlag = MOTIF_STRAND_BOTH;
				if (strcmp(argv[i],"both")==0) strandFlag = MOTIF_STRAND_BOTH;
			}*/
		} else if (strcmp(argv[i],"-allowTargetOverlaps")==0) {
			targetOverlapFlag = 0;
		} else if (strcmp(argv[i],"-allowBgOverlaps")==0) {
			bgOverlapFlag = 0;
		} else {
			fprintf(stderr, "!!! Error: Not sure what command line option \"%s\" is.\n\n", argv[i]);
			printCMDbackground();
		}
	}

	if (alphaLen < 1) {
		if (modelFlag) {
			alphaLen = 5;
		} else {
			alphaLen = 4;
		}
	}

	//Initialize random number generator
	if (seed < 0) {
		struct timeval tv;
		gettimeofday(&tv,NULL);
		unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
		srand(time_in_micros);
	} else {
		fprintf(stderr, "\tInitializing random number generator with seed: %ld\n", seed);
		srand((unsigned long) seed);
	}
	bg->overlapFlag = bgOverlapFlag;
	targets->alphaLen = alphaLen;
	bg->alphaLen = alphaLen;
	targets->setpKmer(kmer,positionFlag);
	bg->setpKmer(kmer,positionFlag);
	//targets->pkmer = kmer;
	//bg->pkmer = kmer;
	//targets->positionFlag = positionFlag;
	//bg->positionFlag = positionFlag;


	if (genomeFasta != NULL) {
		fprintf(stderr, "\n\tReading genome FASTA file (%s):\n", genomeFasta);
		bg->ReadFastaFile(genomeFasta, repMaskFlag);
		bg->shareGenome(targets);
	}

	// load target regions
	if (inputTargetFastaFile != NULL) {
		fprintf(stderr, "\n\tReading target FASTA file (%s):\n", inputTargetFastaFile);
		targets->ReadFastaFile(inputTargetFastaFile, repMaskFlag);
		targets->processSequencesAsIs(strandFlag,maxPercentNperSeq);
	} else if (inputTargetBedFile != NULL) {
		fprintf(stderr, "\n\tReading target genomic positions (%s):\n", inputTargetBedFile);
		targets->processGenomicPositions(inputTargetBedFile,strandFlag, maxPercentNperSeq);
	} else {
		fprintf(stderr, "!!! Error: Need to specify a target FASTA file (-i) or genomic coordinates (-p) file with a genome (-g)\n");
		printCMDbackground();
		exit(0);	
	}

	// set size
	if (size < 0) {
		size = targets->seqLen;
		bg->seqLen = size;
		fprintf(stderr, "\tSetting sequence length to the average target length: %d\n\n", size);
	}
	bg->seqLen = size;
	if (size < 2) {
		fprintf(stderr, "!!! Average Target length is %d - you must specify larger regions to calculate an appropriate background\n",size);
		exit(0);
	}
	if (autoBinSizeFlag) {
		double nseqs = (double)targets->numSeqs;		
		double binEstimate = nseqs/100.0;
		int suggestedBins = (int) floor(sqrt(binEstimate));
		if (suggestedBins < 3) {
			suggestedBins = 3;
		}
		if (suggestedBins > 10) suggestedBins = 10;
		fprintf(stderr, "\t%ld target sequences -> will set -nbins %d and -nsubBins %d\n",targets->numSeqs,suggestedBins,suggestedBins);
		numBins = suggestedBins;
		numSubBins = suggestedBins;
	}

	double* binThresholds = targets->generateBinThresholds(numBins);
	fprintf(stderr, "\n\tAssigning target sequences to bins:\n");
	targets->assignSeqsToBins(binThresholds,numBins);

	targets->targetFreqsByBin = targets->calculatePositionalKmerFrequenciesByBin();
	bg->targetFreqsByBin = targets->targetFreqsByBin;
	targets->scoreKmerProfiles();

	if (modelFlag) {
		// generate bg sequences using a model
		fprintf(stderr, "\tBuilding a background sequence models from target sequences\n");
		bg->numBins = numBins;
		long int numSeqPerBin = numBgSeqs/numBins;
		for (int bin=0;bin<numBins;bin++) {
			fprintf(stderr, "\t\tGenerating sequences for bin=%d (%ld)\n",bin+1,numSeqPerBin);
			BgModel* model = new BgModel();
			model->initializeModel(kmer, alphaLen, positionFlag, size);

			model->buildModel(targets->bins[bin], targets->binSizes[bin]);
			BgSequence** gseqs = model->generateSequences(numSeqPerBin);
			bg->addSequences(gseqs, numSeqPerBin, bin);
		}
		bg->selectAllSequences();

	} else {
		// extract real sequences for the background
		if (genomeFasta != NULL && inputBackgroundBedFile == NULL && inputBackgroundFastaFile == NULL) {
			//maybe load based on input BED file for background
			fprintf(stderr, "\tExtracting background sequences from the genome\n");
			bg->processSequencesAsGenome(size,strandFlag,maxPercentNperSeq,numConsideredBgSeqs,binThresholds,numBins);
		} else if (genomeFasta != NULL && inputBackgroundBedFile != NULL && inputBackgroundFastaFile == NULL) {
			if (bgPosRegionFlag) {
				fprintf(stderr, "\n\tReading background genomic ranges to select from (%s):\n", inputBackgroundBedFile);
				bg->processGenomicRanges(inputBackgroundBedFile,size,strandFlag,maxPercentNperSeq,numConsideredBgSeqs,binThresholds,numBins);
			} else {
				//regions are exact
				fprintf(stderr, "\n\tReading background genomic positions (%s):\n", inputBackgroundBedFile);
				bg->processGenomicPositions(inputBackgroundBedFile,strandFlag,maxPercentNperSeq);
			}
		} else if (inputBackgroundFastaFile != NULL) {
			fprintf(stderr, "\tExtracting background sequences from provided FASTA file (%s)\n", inputBackgroundFastaFile);
			bg->ReadFastaFile(inputBackgroundFastaFile, repMaskFlag);
			bg->processSequencesAsIs(strandFlag,maxPercentNperSeq);
		} else {
			fprintf(stderr, "Error: incorrect options selected for background sequence\n");
			exit(1);
		}

		if (genomeFasta != NULL && inputTargetBedFile != NULL && targetOverlapFlag != 0) {
			fprintf(stderr, "\tMasking target positions in the genome from consideration:\n");
			bg->maskGenomicPositionOverlaps(targets,2);
		}
		if (inputBackgroundMaskBedFile != NULL) {
			//bg->maskGenomicPositions(inputBackgroundMaskBedFile);
		}

		fprintf(stderr, "\tAssigning background sequences to bins:\n");
		bg->assignSeqsToBins(binThresholds,numBins);

		/*bg->scoreKmerProfiles();
		if (pscoreOutputFile != NULL) {
			fprintf(stderr, "\n\tOutputing profile scores to file %s\n", pscoreOutputFile);
			FILE* fp = fopen(pscoreOutputFile,"w");
			if (fp == NULL) {
				fprintf(stderr, "!!! Error, could not open %s to output profile scores\n", pscoreOutputFile);
				exit(0);
			}
			//for (long int i=0;i<targets->numSeqs;i++) targets->seqs[i]->pscore = targets->seqs[i]->stat;
			//for (long int i=0;i<bg->numSeqs;i++) bg->seqs[i]->pscore = bg->seqs[i]->stat;
			targets->printAllSequencesTable(fp);
			bg->printAllSequencesTable(fp);
			fclose(fp);
			exit(0);
		}*/

		// refine selected sequences to match position-specific kmer frequencies
		if (kmer > 0) {
			if (positionFlag) {
				fprintf(stderr, "\n\tSelecting sequences while matching positional nucleotide content (kmer=%d, N=%ld, subBins=%d)\n", kmer,numBgSeqs,numSubBins);
			} else {
				fprintf(stderr, "\n\tSelecting sequences while matching nucleotide content (kmer=%d, N=%ld, subBins=%d)\n", kmer,numBgSeqs,numSubBins);
			}
			matchKmerFrequencies(targets, bg, numBgSeqs, kmer, positionFlag,numSubBins, numkmerIters,overlapIteration,decayRate,pscoreOutputFile);
		} else {
			fprintf(stderr, "\n\tRandomly selecting sequences from bins (N=%ld)\n", numBgSeqs);
			bg->selectInitialSequences(numBgSeqs);
		}
	}

	//output
	long int tseqs = bg->getTotalSelectedSequences();
	fprintf(stderr, "\tTotal Sequences in output: %ld\n", tseqs);	
	if (tseqs < numBgSeqs) {
	}


	char* outputFilename = new char[100000];
	FILE* outfp = NULL;


	sprintf(outputFilename, "%s.bg.sequences.fasta", outputFilePrefix);
	outfp = fopen(outputFilename,"w");
	if (outfp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
		exit(1);
	}
	bg->printSelectedSequencesFASTA(outfp);
	fclose(outfp);

	sprintf(outputFilename, "%s.bg.stats.txt", outputFilePrefix);
	outfp = fopen(outputFilename,"w");
	if (outfp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
		exit(1);
	}
	bg->printSelectedSequencesTable(outfp);
	fclose(outfp);

	sprintf(outputFilename, "%s.bg.positions.bed", outputFilePrefix);
	outfp = fopen(outputFilename,"w");
	if (outfp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
		exit(1);
	}
	bg->printSelectedSequencesBED(outfp);
	fclose(outfp);



	sprintf(outputFilename, "%s.target.sequences.fasta", outputFilePrefix);
	outfp = fopen(outputFilename,"w");
	if (outfp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
		exit(1);
	}
	targets->printAllSequencesFASTA(outfp);
	fclose(outfp);

	sprintf(outputFilename, "%s.target.stats.txt", outputFilePrefix);
	outfp = fopen(outputFilename,"w");
	if (outfp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
		exit(1);
	}
	targets->printAllSequencesTable(outfp);
	fclose(outfp);

	sprintf(outputFilename, "%s.target.positions.bed", outputFilePrefix);
	outfp = fopen(outputFilename,"w");
	if (outfp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
		exit(1);
	}
	targets->printAllSequencesBED(outfp);
	fclose(outfp);

	targets->selectAllSequences();

	if (homerStyleOutputFlag) {
		sprintf(outputFilename, "%s.group.tsv", outputFilePrefix);
		outfp = fopen(outputFilename,"w");
		if (outfp == NULL) {
			fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
			exit(1);
		}
		targets->printGroupFile(outfp,1);
		bg->printGroupFile(outfp,0);
		fclose(outfp);

		sprintf(outputFilename, "%s.seq.tsv", outputFilePrefix);
		outfp = fopen(outputFilename,"w");
		if (outfp == NULL) {
			fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
			exit(1);
		}
		targets->printSeqFile(outfp);
		bg->printSeqFile(outfp);
		fclose(outfp);
	}

	if (kmer < 1) kmer =1;
	targets->setpKmer(kmer,positionFlag);
	bg->setpKmer(kmer,positionFlag);
	double **freq = bg->calculatePositionalKmerFrequencies(bg->sbins, bg->sbinSizes, bg->numBins);
	sprintf(outputFilename, "%s.bg.freq.txt", outputFilePrefix);
	outfp = fopen(outputFilename,"w");
	if (outfp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
		exit(1);
	}
	bg->printPositionalKmerFrequencies(outfp, freq);
	fclose(outfp);
	double **targetFreq = targets->calculatePositionalKmerFrequencies(targets->bins, targets->binSizes, targets->numBins);
	sprintf(outputFilename, "%s.target.freq.txt", outputFilePrefix);
	outfp = fopen(outputFilename,"w");
	if (outfp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
		exit(1);
	}
	targets->printPositionalKmerFrequencies(outfp, targetFreq);
	fclose(outfp);

	double totalDiff = 0.0;
	int flen = 1;
	if (positionFlag) flen = bg->seqLen;
	for (int j=0;j<flen;j++) {
		for (int k=0;k<bg->freqLength;k++) {
			//double d = targetFreqs[j][k]-bgFreqs[j][k];
			double d = freq[j][k]-targetFreq[j][k];
			totalDiff += fabs(d);
		}
	}
	fprintf(stderr, "\n\tFinal frequency difference = %lf\n", totalDiff);

	sprintf(outputFilename, "%s.info.txt", outputFilePrefix);
	outfp = fopen(outputFilename,"w");
	if (outfp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file %s for writing!\n", outputFilename);
		exit(1);
	}
	fprintf(outfp, "cmd =");
	for (int i=0;i<argc;i++) {
		fprintf(outfp, " %s", argv[i]);
	}
	fprintf(outfp, "\n");
	fprintf(outfp, "total bg sequences selected = %ld\n", tseqs);
	fprintf(outfp, "frequency difference = %lf\n", totalDiff);
	fclose(outfp);

	fprintf(stderr, "\n");

}


void printCMDbackground() {
	fprintf(stderr, "\n\thomer2 background -i <target sequences.fasta> [options]\n");
	fprintf(stderr, "\n\thomer2 bg -i <target sequences.fasta> [options]\n");
	fprintf(stderr, "\n\tGenerate/Select background sequences that match properties in a set of target sequences.\n");
	fprintf(stderr, "\n\tInputs:\n");
	fprintf(stderr, "\t  Target sequences you want to model:\n");
	fprintf(stderr, "\t\t-i <target sequences.fasta> (FASTA file)\n");
	fprintf(stderr, "\t\t-p <target positions.bed> (Alteratively, provide a BED or HOMER peak file with genomic coordinates)\n");
	fprintf(stderr, "\t  Background sequences to select from:\n");
	fprintf(stderr, "\t\t-model (generate sequences using a model, do not extract real background sequences)\n");
	fprintf(stderr, "\t\t-g <genome.fasta> (genome FASTA file or seqeunce resource to select sequences from)\n");
	fprintf(stderr, "\t\t-b <background sequences.fasta> (explicit set of background sequences to choose from FASTA file)\n");
	fprintf(stderr, "\t\t-bg <background positions.bed> (explicit set of background positions to choose from)\n");
	fprintf(stderr, "\t\t-bgr <background regions positions.bed> (regions of the genome to select bg sequences from)\n");
	//fprintf(stderr, "\t\t-bgm <masked positions.bed> (regions of the genome to mask - do not select regions from)\n");
	fprintf(stderr, "\tKey options:\n");
	fprintf(stderr, "\t\t-size <#> (size of regions to consider in background, default: avg of length of target sequences)\n");
	fprintf(stderr, "\t\t-N <#> (number of background sequences to select, default: 100000)\n");
	fprintf(stderr, "\t\t-NN <#> (number of background sequences initial screen from genome, default: 100000000)\n");
	fprintf(stderr, "\t\t-mask (mask lowercase sequence i.e. softmasked sequence, default: use all sequences)\n");
	fprintf(stderr, "\t\t-nbins <#> (number of bins to segregate sequences into for GC selection, def: 10)\n");
	fprintf(stderr, "\t\t-nsubBins <#> (number of bins to segregate sequences into for positional frequencies, def: 10)\n");
	fprintf(stderr, "\t\t\t(if not set, homer2 may attempt to adjust -nbins/-nsubBins automatically when using small datasets)\n");
	fprintf(stderr, "\t\t-maxFractionN <#> (Maximum fraction of sequence that can be N and still used, default: 0.5)\n");
	fprintf(stderr, "\t\t-allowTargetOverlaps (allow selected bg sequences from a genome to overlap targets, def: not allowed)\n");
	fprintf(stderr, "\t\t-allowBgOverlaps (allow selected bg sequences from a genome to overlap, def: not allowed)\n");
	fprintf(stderr, "\t\t-strand (allow sequences to overlap if on separate strands)\n");
	fprintf(stderr, "\t\t-pkmer <#> (match positional kmer content)\n");
	fprintf(stderr, "\t\t-ikmer <#> (match overall kmer content [position independent])\n");
	fprintf(stderr, "\t\t-excludeNs/-includeNs (by default, kmers with Ns are excluded when selecting bg sequences,\n");
	fprintf(stderr, "\t\t\t\tbut included when generating sequences with -model)\n");
	fprintf(stderr, "\t\t-pscore <outputBEDfile> (Report initial pscores)\n");
	fprintf(stderr, "\t\t-maxIterations <#> (maximum iterations, def: 20)\n");
	fprintf(stderr, "\t\t-overlapIteration <#> (iteration to start enforcing no overlaps, def: 5)\n");
	fprintf(stderr, "\t\t-decayRate <#> (selection rate per iteration, def: 0.75)\n");
	fprintf(stderr, "\t\t-seed <#> (seed for random number generator, def: uses time)\n");
	//fprintf(stderr, "\t\t-maxIterations <#> (maximum iterations, def: 10)\n");
	
	fprintf(stderr, "\tOutput:\n");
	fprintf(stderr, "\t\t-o <output prefix> (default: out)\n");
	fprintf(stderr, "\t\t-gs (include homer-style group and sequence output files)\n");
	fprintf(stderr, "\n");
	exit(0);
}


BgSequence::BgSequence() {
	seq = NULL;
	len = -1;
	name = NULL;
	pos = -1;
	strand=0;
	stat=0.0;
	mask=0;
	weight = 1.0;
	bin = -1;
	positionIndex=-1;
} 
BgSequence::~BgSequence() {
}

// returns a new string with dedicated memory
char* BgSequence::getSeq() {
	if (len < 1) return NULL;
	if (seq == NULL) return NULL;
	char* s = new char[len+1];
	memcpy(s,seq,len);
	s[len]='\0';
	return s;
}
void BgSequence::printAll(FILE* fp) {
	char* s = getSeq();
	fprintf(fp, "%s\n", s);
	delete []s;
	char* n = name;
	if (name == NULL) n = (char*) "Random";
	fprintf(fp, "len=%d, stat=%lf, mask=%d, positionIndex=%ld\n", len, stat, mask, positionIndex);
	fprintf(fp, "src=%s, %d, %d\n", n, pos, strand);
}
void BgSequence::printBED(FILE* fp) {
	printBED(fp,-1);
}
void BgSequence::printBED(FILE* fp,int index) {
	int start = pos;
	int end = pos+len;
	char s = '+';
	if (strand == MOTIF_STRAND_NEG) {
		s = '-';
	}
	char* n = name;
	if (name == NULL) n = (char*) "Random";
	if (index < 0) {
		fprintf(fp, "%s\t%d\t%d\t%s_%d_%d\t%lf\t%c\n",n,start,end,n,pos,strand,weight,s);
	} else {
		fprintf(fp, "%s\t%d\t%d\t%s_%d_%d_%d\t%lf\t%c\n",n,start,end,n,pos,strand,index,weight,s);
	}
}
void BgSequence::printTable(FILE* fp) {
	printTable(fp,-1);
}
void BgSequence::printTable(FILE* fp,int index) {
	char s = '+';
	if (strand == MOTIF_STRAND_NEG) s='-';
	double gc = getGC();	
	char* n = name;
	if (name == NULL) n = (char*) "Random";
	if (index < 0) {
		fprintf(fp, "%s_%d_%d\t%s\t%d\t%c\t%lf\t%lf\t%lf\t%d\n",n,pos,strand,n,pos,s,weight,pscore,gc,bin);
	} else {
		fprintf(fp, "%s_%d_%d_%d\t%s\t%d\t%c\t%lf\t%lf\t%lf\t%d\n",n,pos,strand,index,n,pos,s,weight,pscore,gc,bin);
	}
}
void BgSequence::printFASTA(FILE* fp) {
	printFASTA(fp,-1);
}
void BgSequence::printFASTA(FILE* fp, int index) {
	char* n = name;
	if (name == NULL) n = (char*) "Random";
	if (index >=0) {
		fprintf(fp, ">%s_%d_%d_%d %lf\n",n,pos,strand,index,weight);
	} else {
		fprintf(fp, ">%s_%d_%d %lf\n",n,pos,strand,weight);
	}
	char* s = getSeq();
	fprintf(fp, "%s\n", s);
	delete []s;
}
void BgSequence::printGroup(FILE* fp,int assignment) {
	printGroup(fp,assignment,-1);
}
void BgSequence::printGroup(FILE* fp, int assignment, int index) {
	char* n = name;
	if (name == NULL) n = (char*) "Random";
	if (index >=0) {
		fprintf(fp, "%s_%d_%d_%d\t%d\t%lf\n",n,pos,strand,index,assignment,weight);
	} else {
		fprintf(fp, "%s_%d_%d\t%d\t%lf\n",n,pos,strand,assignment,weight);
	}
}
void BgSequence::printSeq(FILE* fp) {
	printSeq(fp,-1);
}
void BgSequence::printSeq(FILE* fp,int index) {
	char* n = name;
	if (name == NULL) n = (char*) "Random";
	char* s = getSeq();
	if (index >=0) {
		fprintf(fp, "%s_%d_%d_%d\t%s\n",n,pos,strand,index,s);
	} else {
		fprintf(fp, "%s_%d_%d\t%s\n",n,pos,strand,s);
	}
	delete []s;
}

double BgSequence::getGC() {
	int gc = 0;
	int total = 0;
	for (int i=0;i<len;i++) {
		if (seq[i] == 'N') continue;
		total++;
		if (seq[i] == 'C' || seq[i] == 'G') {
			gc++;
		}
	}
	if (total < 1) return 0.0;
	return ((double)gc)/((double)total);
}

BgSequenceArray::BgSequenceArray() {
	seqs = NULL; //BgSequence**
	numSeqs = 0; //int
	seqLen = -1; //int
	ogSeqs = NULL; //char**
	ogSeqsRV = NULL; //char**
	ogLens = NULL; //char**
	ogNames = NULL; //char**
	numOGs = 0; //int
	posStrandArray = NULL; //BgSequence**
	negStrandArray = NULL; //BgSequence**
	shareFlag = 0;
	chr2Index = NULL;

	bins = NULL;
	binSizes = NULL;
	numBins = -1;
	sbins = NULL;
	sbinSizes = NULL;
	numBgSeqs = -1;
	numBgSeqPerBin = -1;

	targetFreqsByBin=NULL;

	positionFlag = 0;
	pkmer=0;
	alphaLen = 5;

	maxFailures = 100;
	overlapFlag = 1;
}
BgSequenceArray::~BgSequenceArray() {
	if (seqs != NULL) {
		for (int i=0;i<numSeqs;i++) {
			if (seqs[i] != NULL) delete seqs[i];
			seqs[i] = NULL;
		}
		delete []seqs;
		seqs = NULL;
	}
	if (shareFlag == 0) {
		for (int i=0;numOGs;i++) {
			if (ogSeqs!=NULL && ogSeqs[i] != NULL) delete [](ogSeqs[i]);
			if (ogSeqsRV!=NULL && ogSeqsRV[i] != NULL) delete [](ogSeqsRV[i]);
			if (ogNames!=NULL && ogNames[i] != NULL) delete [](ogNames[i]);
		}
		if (ogSeqs!=NULL) delete []ogSeqs;
		if (ogSeqsRV!=NULL) delete []ogSeqsRV;
		if (ogNames!=NULL) delete []ogNames;
		if (ogLens!=NULL) delete []ogLens;
		if (posStrandArray != NULL) delete []posStrandArray;
		if (negStrandArray != NULL) delete []negStrandArray;
		if (chr2Index != NULL) delete chr2Index;
	}
	numSeqs = 0; //int
	seqLen = -1; //int
}

void BgSequenceArray::ReadFastaFile(char* filename, char repMaskFlag) {

	numOGs = -1;
	totalOGSeqLength=0;

	FILE* fp = fopen(filename,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file \"%s\"\n", filename);
		exit(1);
	}

	char* curSeq = new char[MAX_CHR_LENGTH];
	long int curSeqLen = 0;
	char* buffer = new char[BUFFER];
	int maxOGs = 1000000;
	int maxOGsInc = 1000000;
	ogSeqs = new char*[maxOGs];
	ogLens = new long int[maxOGs];
	ogNames = new char*[maxOGs];
	char** cols = new char*[1000];
	int numCols = 0;

	while (fgets(buffer,BUFFER,fp)!=NULL) {
		int lineLen = strlen(buffer);
		if (lineLen > 0 && buffer[lineLen-1] == '\n') {
			buffer[lineLen-1] = '\0';
			lineLen--;
		}
		if (buffer[0] == '>') {
			if (numOGs >= 0) {
				ogSeqs[numOGs] = new char[curSeqLen+1];
				ogLens[numOGs] = curSeqLen;
				strcpy(ogSeqs[numOGs],curSeq);
				totalOGSeqLength+=curSeqLen;
			}
			numOGs++;

			//expand array if needed
			if (numOGs+2 >= maxOGs) {
				maxOGs += maxOGsInc;
				char** newOGs = new char*[maxOGs];
				long int* newOGLens = new long int[maxOGs];
				char** newOGNames = new char*[maxOGs];
				for (int i=0;i<numOGs-1;i++){ 
					newOGs[i] = ogSeqs[i];
					newOGLens[i] = ogLens[i];
					newOGNames[i] = ogNames[i];
				}
				delete []ogSeqs;
				delete []ogNames;
				delete []ogLens;
				ogSeqs = newOGs;
				ogLens = newOGLens;
				ogNames = newOGNames;
			}

			split2(buffer, cols, numCols,WHITE_SPACE);
			int nameLength = strlen(cols[0])+1;
			ogNames[numOGs] = new char[nameLength];
			strcpy(ogNames[numOGs],&(cols[0][1]));
			curSeqLen = 0;
			curSeq[0] = '\0';
			continue;
		} else {
			for (int i=0;i<lineLen;i++) {
				if (buffer[i] == 'a' || buffer[i] == 'A') {
					if (repMaskFlag && buffer[i] == 'a') buffer[i] = 'N';
					else buffer[i] = 'A';
				} else if (buffer[i] == 'c' || buffer[i] == 'C') {
					if (repMaskFlag && buffer[i] == 'c') buffer[i] = 'N';
					else buffer[i] = 'C';
				} else if (buffer[i] == 'g' || buffer[i] == 'G') {
					if (repMaskFlag && buffer[i] == 'g') buffer[i] = 'N';
					else buffer[i] = 'G';
				} else if (buffer[i] == 't' || buffer[i] == 'T') {
					if (repMaskFlag && buffer[i] == 't') buffer[i] = 'N';
					else buffer[i] = 'T';
				} else if (buffer[i] == 'n' || buffer[i] == 'N') {
					buffer[i] = 'N';
				} else {
					fprintf(stderr, "\t\t!!! Warning: Invalid sequence character: %d\n",buffer[i]);
					buffer[i] = 'N';
				}
			}
			if (lineLen+10 > MAX_CHR_LENGTH) {
				fprintf(stderr, "!!! Error: sequence length is greater than MAX_CHR_LENGTH (%d)\n", MAX_CHR_LENGTH);
			}
			memcpy(&(curSeq[curSeqLen]),buffer,lineLen);
			curSeqLen += lineLen;
			curSeq[curSeqLen]='\0';
		}
			
	}
	if (numOGs >= 0) {
		ogSeqs[numOGs] = new char[curSeqLen+1];
		ogLens[numOGs] = curSeqLen;
		strcpy(ogSeqs[numOGs],curSeq);
		totalOGSeqLength+=curSeqLen;
		numOGs++;
	}
	fclose(fp);

	//Setup Reverse opposite sequences
	ogSeqsRV = new char*[numOGs];
	for (int i=0;i<numOGs;i++) {
		ogSeqsRV[i] = revopp2(ogSeqs[i]);
	}

	//finish by creating a hashtable for lookup later
	fprintf(stderr, "\t\tSummary of sequences found (first 20):\n");
	fprintf(stderr, "\t\t\tName\tLength\n");
	chr2Index = new Inttable(numOGs*2);
	for (int i=0;i<numOGs;i++) {
		chr2Index->insert(i,ogNames[i]);
		if (i<20) fprintf(stderr, "\t\t\t%s\t%ld\n", ogNames[i], ogLens[i]);
	}

	delete []buffer;
	delete []curSeq;
	delete []cols;

	fprintf(stderr, "\t\tFound %d sequences (%ld total length) in file: %s\n",numOGs,totalOGSeqLength,filename);

}
void BgSequenceArray::shareGenome(BgSequenceArray* target) {
	target->ogSeqs = ogSeqs;
	target->ogSeqsRV = ogSeqsRV;
	target->ogNames = ogNames;
	target->ogLens = ogLens;
	target->numOGs = numOGs;
	target->totalOGSeqLength = totalOGSeqLength;
	target->posStrandArray = posStrandArray;
	target->negStrandArray = negStrandArray;
	target->numStrandArray = numStrandArray;
	target->chr2Index = chr2Index;
		
	target->shareFlag = 1;
}

void BgSequenceArray::addSequences(BgSequence** newSeqs, long int numNewSeqs, int bin) {


	BgSequence** nseqs = new BgSequence*[numSeqs+numNewSeqs];
	for (int i=0;i<numSeqs;i++) {
		nseqs[i] = seqs[i];
	}
	for (int i=0;i<numNewSeqs;i++) {
		nseqs[numSeqs+i] = newSeqs[i];
	}
	delete []seqs;
	numSeqs += numNewSeqs;
	seqs = nseqs;

	if (bins == NULL) {
		bins = new BgSequence**[numBins];
		binSizes = new long int[numBins];
		for (int i=0;i<numBins;i++) {
			bins[i] = NULL;
			binSizes[i] = 0;
		}
	}
	int newBinSize = binSizes[bin]+numNewSeqs;
	BgSequence** nbin = new BgSequence*[newBinSize];
	if (bins[bin] != NULL) {
		for (int i=0;i<binSizes[bin];i++) {
			nbin[i] = bins[bin][i];
		}
	}
	for (int i=0;i<numNewSeqs;i++) {
		nbin[binSizes[bin]+i] = newSeqs[i];
	}
	delete []bins[bin];
	binSizes[bin] += numNewSeqs;
	bins[bin] = nbin;
	
}

void BgSequenceArray::processGenomicPositions(char* bedFileName, char strandFlag, double maxNpercent) {

	strand = strandFlag;	
	FILE* fp = fopen(bedFileName,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file \"%s\"\n", bedFileName);
		exit(1);
	}

	char* buffer = new char[BUFFER];
	char** cols = new char*[1000];
	int numCols = 0;
	int totalPossible = 0; 
	int numBad = 0;
	int numRemovedDueToNs = 0;


	double avgLen = 0.0;
	double avgGC = 0.0;

	int numLines = 0;
	while (fgets(buffer,BUFFER,fp)!=NULL) {
		numLines++;
	}
	rewind(fp);

	seqs = new BgSequence*[numLines];
	numSeqs = 0;
	int format = -1;

	while (fgets(buffer,BUFFER,fp)!=NULL) {
		int lineLen = strlen(buffer);
		if (lineLen > 0 && buffer[lineLen-1] == '\n') {
			buffer[lineLen-1] = '\0';
			lineLen--;
		}
		if (buffer[0] == '#') {
			continue;
		}
		split2(buffer, cols, numCols,WHITE_SPACE);
		if (numCols < 3) {
			continue;
		}
		totalPossible++;
		char* chr = cols[0];
		char* startStr = cols[1];
		char* endStr = cols[2];
		char* strandStr = NULL;
		int startAdj = 0;
		int lenAdj = 0;
		if (numCols > 3 && (strcmp(cols[3],"-")==0 || strcmp(cols[3],"+")==0 || strcmp(cols[3],".")==0)) {
			// Assume BED file as chr, start, end, strand
			if (format != 1) {
				fprintf(stderr, "\t\tgussing file is BED format (chr <tab> start <tab> end <tab> strand)\n");
				format = 1;
			}
			strandStr = cols[3];
		} else if (numCols > 4 && (strcmp(cols[4],"-")==0 || strcmp(cols[4],"+")==0 || strcmp(cols[4],".")==0)) {
			// Assume HOMER Peak file as name, chr,  start,  end,  strand
			if (format != 2) {
				fprintf(stderr, "\t\tgussing file is HOMER peak format (name <tab> chr <tab> start <tab> end <tab> strand)\n");
				format = 2;
			}
			chr = cols[1];
			startStr =cols[2];
			endStr = cols[3];
			strandStr = cols[4];
			startAdj = -1;
			lenAdj = 0;
		} else if (numCols > 5 && (strcmp(cols[5],"-")==0 || strcmp(cols[5],"+")==0 || strcmp(cols[5],".")==0)) {
			// Assume BED file as chr,  start,  end, name, value, strand
			if (format != 2) {
				fprintf(stderr, "\t\tgussing file is BED format (chr <tab> start <tab> end <tab> name <tab> value <tab> strand)\n");
				format = 2;
			}
			strandStr = cols[5];
		} else {
			if (format < 0) {
				fprintf(stderr, "\t\tUnknown format\n");
				format = 0;
			}
		}
		int start = -1;
		int end = -1;
		sscanf(startStr,"%d",&start);
		sscanf(endStr,"%d",&end);
		char strand = MOTIF_STRAND_POS;
		if (strandStr != NULL) {
			if (strandStr[0] == '-' || strandStr[0] == '1') {
				strand = MOTIF_STRAND_NEG;
			}
		}
		//fprintf(stderr, "strand=%d\n", strand);
		start += startAdj; //for 1 vs. 0 based positions a la BED files
	
		//first check if chr is in the genome
		int index = chr2Index->search(chr);
		if (index == EMPTY_INT) {
			fprintf(stderr, "\t\t!!! Warning: chromosome %s not found in genome!\n", chr);
			numBad++;
			continue;
		}
		if (start < 0 || end < 0 || start > ogLens[index] || end > ogLens[index]) {
			fprintf(stderr, "\t\t!!! Warning: positions %d,%d are outside of range for chromosome %s (0-%ld)\n",start,end,chr,ogLens[index]);
			numBad++;
			continue;
		}

		//Ok - sequence should be extractable
		int len = end-start + lenAdj;

		if (len != seqLen && positionFlag) {
			//fprintf(stderr, "Error - background region not the correct size (needs to be exact for a positional model)\n");
			//exit(0);
		}

		char* seq = &(ogSeqs[index][start]);
		if (strand == MOTIF_STRAND_NEG) {
			seq = &(ogSeqs[index][start-len]);
			int rvPos = ogLens[index]-start-len;
			seq = &(ogSeqsRV[index][rvPos]);
		}

		int numNs = 0;
		int numGC = 0;
		int totalNonN = 0;
		for (int i=0;i<len;i++) {
			if (seq[i] == 'N') numNs++;
			else {
				if (seq[i] == 'G' || seq[i] == 'C') numGC++;
				totalNonN++;
			}
		}

		if (((double)numNs)/((double)len) > maxNpercent || totalNonN < 1) {
			numRemovedDueToNs++;
			continue;
		}

		double GC = ((double)numGC)/((double)totalNonN); //fraction GC
		avgGC+=GC;
		avgLen += end-start;

		BgSequence* s = new BgSequence();
		s->seq = seq;
		s->len = len;
		s->name = ogNames[index];
		s->pos = start;
		s->strand = strand;
		s->stat = -1*GC;
		s->mask = 0;
		s->weight = 1.0;
		s->bin = -1;
		s->positionIndex=-1;

		seqs[numSeqs++]=s;
	}
	fclose(fp);

	fprintf(stderr, "\t\t%d of %d removed (fraction of Ns > %.3lf)\n", numRemovedDueToNs, totalPossible, maxNpercent);
	if (numSeqs < 1) {
		fprintf(stderr, "!!! Error, no sequences were processed from BED/peak file!!!\n");
		return;
	}
	fprintf(stderr, "\t\t%d removed do to bad chromosome or bad coordinates\n", numBad);
	avgGC /= numSeqs;
	avgLen /= numSeqs;
	fprintf(stderr, "\t\tProcessed %ld target sequence positions from the genome (avg. length = %.1lf, GC = %lf)\n", numSeqs, avgLen, avgGC);
	seqLen = (int) avgLen;
	if (avgLen < 10) {
		fprintf(stderr, "!!!!! WARNING: Average target sequence length is < 10 bp (%.1lf). You may want to change the size of the\n", avgLen);
		fprintf(stderr, "               input regions (you can use adjustPeakFile.pl for this)\n");
	}


	delete []buffer;
	delete []cols;
	fprintf(stderr, "\n");
}
void BgSequenceArray::processSequencesAsGenome(int seqLength, char strandSetting, double maxNpercent,
												long int numConsideredBgSeqs,double *binThresholds,int nBins) {

	strand = strandSetting;
	seqLen = seqLength;

	double avgGC = 0.0;
	double avgGCselected = 0.0;
	int numRemovedDueToNs = 0;

	long int approxMaxSeqs = totalOGSeqLength*2; // for reverse opposites 

	seqs = new BgSequence*[approxMaxSeqs];
	posStrandArray = new BgSequence*[totalOGSeqLength];
	negStrandArray = new BgSequence*[totalOGSeqLength];
	numSeqs = 0;
	numStrandArray = 0;
	
	long int totalPossible = 0;
	long int totalPreselection = 0;

	int rateRebalanceCount = 100000;
	int curRebalanceCount = 0;
	long int* binTotals = NULL;
	long int* binTotalselected = NULL;
	double* binRates = NULL;
	if (numConsideredBgSeqs > 0) {
		binTotals = new long int[nBins];
		binTotalselected = new long int[nBins];
		binRates = new double[nBins];
		double initialRate = ((double)numConsideredBgSeqs)/((double)approxMaxSeqs);
		for (int i=0;i<nBins;i++) {
			binTotals[i] = 0;
			binTotalselected[i] = 0;
			binRates[i] = initialRate;
		}
		fprintf(stderr, "\t\tWill select %ld from the approximately %ld possible sequences (%lf) [showing first 20]:\n", numConsideredBgSeqs, approxMaxSeqs, initialRate);
	}


	for (int i=0;i<numOGs;i++) {
		if (i < 20) fprintf(stderr, "\t\tprocessing %s (length=%ld)...\n", ogNames[i], ogLens[i]);
		int numNs = 0;
		int numGC = 0;
		int totalNonN = 0;
		for (int j=0;j<ogLens[i]-seqLen;j++) {
			totalPossible++;
			if (j==0) {
				for (int k=j;k<j+seqLen;k++) {
					if (ogSeqs[i][k] == 'N') numNs++;
					else {
						if (ogSeqs[i][k] == 'G' || ogSeqs[i][k] == 'C') numGC++;
						totalNonN++;
					}
				}
			} else {
				int x = j-1;
				if (ogSeqs[i][x] == 'N') {
					numNs--;
				} else {
					if (ogSeqs[i][x] == 'G' || ogSeqs[i][x] == 'C') numGC--;
					totalNonN--;
				}
				x = j+seqLen-1;
				if (ogSeqs[i][x] == 'N') {
					numNs++;
				} else {
					if (ogSeqs[i][x] == 'G' || ogSeqs[i][x] == 'C') numGC++;
					totalNonN++;
				}
			}
			if (((double)numNs)/((double)seqLen) > maxNpercent) {
				numRemovedDueToNs++;
				continue;
			}
			double GC = ((double)numGC)/((double)totalNonN); //fraction GC
			avgGC+=GC;

			totalPreselection++;


			if (numConsideredBgSeqs > 0) {
				curRebalanceCount++;
				int bindex = 0;
				for (int k=0;k<nBins;k++) {
					if (-1*GC < binThresholds[k]) {
						bindex = k;
						break;
					}
				}
				binTotals[bindex]+=2;

				if (curRebalanceCount > rateRebalanceCount) {
					curRebalanceCount=0;
					//time to rebalance rates
					long int t = 0;
					for (int k=0;k<nBins;k++) {
						t += binTotals[k];
					}	
					for (int k=0;k<nBins;k++) {
						binRates[k] = ((double)numConsideredBgSeqs/(((double)nBins))) /
									(((double)approxMaxSeqs)*(((double)binTotals[k])/((double)t)));
					}
				}

				if (((double)rand())/((double)RAND_MAX) < binRates[bindex]) {
					//add
				} else {
					continue;
				}
				binTotalselected[bindex]+=2;
			}

			avgGCselected += GC;

			BgSequence* s = new BgSequence();
			s->seq = &(ogSeqs[i][j]);
			s->len = seqLen;
			s->name = ogNames[i];
			s->pos = j;
			s->strand = MOTIF_STRAND_POS;
			s->stat = -1*GC;
			s->mask = 0;
			s->weight = 1.0;
			s->bin = -1;
			s->positionIndex=numStrandArray;
			posStrandArray[numStrandArray]=s;

			seqs[numSeqs++]=s;
			
			BgSequence* ss = new BgSequence();
			int rvPos = ogLens[i]-j-seqLen;
			ss->seq = &(ogSeqsRV[i][rvPos]);
			ss->len = seqLen;
			ss->name = ogNames[i];
			ss->pos = j;
			ss->strand = MOTIF_STRAND_NEG;
			ss->stat = -1*GC;
			ss->mask = 0;
			ss->weight = 1.0;
			ss->bin = -1;
			ss->positionIndex=numStrandArray;
			negStrandArray[numStrandArray]=ss;
			seqs[numSeqs++]=ss;

			numStrandArray++;

		}
	}

	double r = ((double)numRemovedDueToNs)/((double)totalPossible);
	fprintf(stderr, "\t\t%d of %ld (%.2lf%%) removed (fraction of Ns > %.3lf)\n", numRemovedDueToNs, totalPossible, r, maxNpercent);
	if (numSeqs < 1) {
		fprintf(stderr, "!!! Error, no sequences were processed!!!\n");
		return;
	}
	avgGC /= totalPreselection;
	fprintf(stderr, "\t\tProcessed %ld sequences from input sequences (avg. GC = %lf)\n", totalPreselection*2, avgGC);
	avgGCselected /= numStrandArray;
	if (numConsideredBgSeqs > 0) {
		fprintf(stderr, "\t\tSelected %ld sequences (avg. GC = %lf)\n", numSeqs, avgGCselected);
		fprintf(stderr, "\t\t\tBin\tThreshold\tCount\tCount Selected\tFinal Bin Rate\n");
		for (int i=0;i<nBins;i++) {
			if (i==nBins-1) fprintf(stderr, "\t\t\t%d\tmaximum  \t%ld\t%ld\t%lf\n",i+1,binTotals[i],binTotalselected[i],binRates[i]);
			else fprintf(stderr, "\t\t\t%d\t%lf\t%ld\t%ld\t%lf\n",i+1,binThresholds[i],binTotals[i],binTotalselected[i],binRates[i]);
		}
		delete []binTotals;
		delete []binTotalselected;
		delete []binRates;
	}
	fprintf(stderr, "\n");

}



void BgSequenceArray::processGenomicRanges(char* bedFileName, int seqLength, char strandFlag, double maxNpercent,
												long int numConsideredBgSeqs,double *binThresholds,int nBins) {

	strand = strandFlag;	
	seqLen = seqLength;

	FILE* fp = fopen(bedFileName,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Error: Could not open file \"%s\"\n", bedFileName);
		exit(1);
	}

	char* buffer = new char[BUFFER];
	char** cols = new char*[1000];
	int numCols = 0;
	int numBad = 0;


	int numRegions = 0;
	while (fgets(buffer,BUFFER,fp)!=NULL) {
		numRegions++;
	}
	rewind(fp);


	int format = -1;
	BgPosition* regions = new BgPosition[numRegions];
	numRegions = 0;
	long int totalLength = 0;

	while (fgets(buffer,BUFFER,fp)!=NULL) {
		int lineLen = strlen(buffer);
		if (lineLen > 0 && buffer[lineLen-1] == '\n') {
			buffer[lineLen-1] = '\0';
			lineLen--;
		}
		if (buffer[0] == '#') {
			continue;
		}
		split2(buffer, cols, numCols,WHITE_SPACE);
		if (numCols < 3) {
			continue;
		}
		char* name = cols[0];
		char* chr = cols[0];
		char* startStr = cols[1];
		char* endStr = cols[2];
		char* strandStr = NULL;
		int startAdj = 0;
		int lenAdj = 0;
		if (numCols > 3 && (strcmp(cols[3],"-")==0 || strcmp(cols[3],"+")==0 || strcmp(cols[3],".")==0)) {
			// Assume BED file as chr, start, end, strand
			if (format != 1) {
				fprintf(stderr, "\t\tgussing file is BED format (chr <tab> start <tab> end <tab> strand)\n");
				format = 1;
			}
			strandStr = cols[3];
		} else if (numCols > 4 && (strcmp(cols[4],"-")==0 || strcmp(cols[4],"+")==0 || strcmp(cols[4],".")==0)) {
			// Assume HOMER Peak file as name, chr,  start,  end,  strand
			if (format != 2) {
				fprintf(stderr, "\t\tgussing file is HOMER peak format (name <tab> chr <tab> start <tab> end <tab> strand)\n");
				format = 2;
			}
			name = cols[0];
			chr = cols[1];
			startStr =cols[2];
			endStr = cols[3];
			strandStr = cols[4];
			startAdj = -1;
			lenAdj = 0;
		} else if (numCols > 5 && (strcmp(cols[5],"-")==0 || strcmp(cols[5],"+")==0 || strcmp(cols[5],".")==0)) {
			// Assume BED file as chr,  start,  end, name, value, strand
			if (format != 2) {
				fprintf(stderr, "\t\tgussing file is BED format (chr <tab> start <tab> end <tab> name <tab> value <tab> strand)\n");
				format = 2;
			}
			name = cols[3];
			strandStr = cols[5];
		}
		int start = -1;
		int end = -1;
		sscanf(startStr,"%d",&start);
		sscanf(endStr,"%d",&end);
		char regStrand = MOTIF_STRAND_POS;
		if (strandStr != NULL) {
			if (strandStr[0] == '-' || strandStr[0] == '1') {
				regStrand = MOTIF_STRAND_NEG;
			}
		}
		//fprintf(stderr, "strand=%d\n", strand);
		start += startAdj; //for 1 vs. 0 based positions a la BED files
	
		//first check if chr is in the genome
		int index = chr2Index->search(chr);
		if (index == EMPTY_INT) {
			fprintf(stderr, "\t\t!!! Warning: chromosome %s not found in genome!\n", chr);
			numBad++;
			continue;
		}
		if (start < 0 || end < 0 || start > ogLens[index] || end > ogLens[index]) {
			fprintf(stderr, "\t\t!!! Warning: positions %d,%d are outside of range for chromosome %s (0-%ld)\n",start,end,chr,ogLens[index]);
			numBad++;
			continue;
		}

		//Ok - sequence should be extractable
		int len = end-start + lenAdj;

		if (len < seqLen) {
			//region to short to consider
			continue;
		}
		regions[numRegions].set(name,chr,start,end,regStrand,len);
		regions[numRegions].chrIndex = index;
		numRegions++;
		totalLength += len;
	}
	fclose(fp);

	qsort(regions,numRegions,sizeof(BgPosition),&cmpBgPositions);

	//should probably look for regionst tht overlap... but out of time for that right now!

	double avgGC = 0.0;
	double avgGCselected = 0.0;
	int numRemovedDueToNs = 0;

	long int approxMaxSeqs = totalLength*2; // for reverse opposites 
	seqs = new BgSequence*[approxMaxSeqs];
	posStrandArray = new BgSequence*[totalLength];
	negStrandArray = new BgSequence*[totalLength];
	numSeqs = 0;
	numStrandArray = 0;
	
	long int totalPossible = 0;
	long int totalPreselection = 0;

	int rateRebalanceCount = 100000;
	int curRebalanceCount = 0;
	long int* binTotals = NULL;
	long int* binTotalselected = NULL;
	double* binRates = NULL;
	if (numConsideredBgSeqs > 0) {
		binTotals = new long int[nBins];
		binTotalselected = new long int[nBins];
		binRates = new double[nBins];
		double initialRate = ((double)numConsideredBgSeqs)/((double)approxMaxSeqs);
		for (int i=0;i<nBins;i++) {
			binTotals[i] = 0;
			binTotalselected[i] = 0;
			binRates[i] = initialRate;
		}
		fprintf(stderr, "\t\tWill select %ld from the approximately %ld possible sequences (%lf) [showing first 20]:\n", numConsideredBgSeqs, approxMaxSeqs, initialRate);
	}


	for (int i=0;i<numRegions;i++) {

		int chrIndex = regions[i].chrIndex;
		int start = regions[i].start;
		int end = regions[i].end;
		//char regStrand = regions[i].strand;

		int numNs = 0;
		int numGC = 0;
		int totalNonN = 0;
		for (int j=start;j<end-seqLen;j++) {
			if (j < 0) continue;
			if (j >= ogLens[chrIndex]-seqLen) break;

			totalPossible++;
			if (j==start) {
				for (int k=j;k<j+seqLen;k++) {
					if (ogSeqs[chrIndex][k] == 'N') numNs++;
					else {
						if (ogSeqs[chrIndex][k] == 'G' || ogSeqs[chrIndex][k] == 'C') numGC++;
						totalNonN++;
					}
				}
			} else {
				int x = j-1;
				if (ogSeqs[chrIndex][x] == 'N') {
					numNs--;
				} else {
					if (ogSeqs[chrIndex][x] == 'G' || ogSeqs[chrIndex][x] == 'C') numGC--;
					totalNonN--;
				}
				x = j+seqLen-1;
				if (ogSeqs[chrIndex][x] == 'N') {
					numNs++;
				} else {
					if (ogSeqs[chrIndex][x] == 'G' || ogSeqs[chrIndex][x] == 'C') numGC++;
					totalNonN++;
				}
			}
			if (((double)numNs)/((double)seqLen) > maxNpercent) {
				numRemovedDueToNs++;
				continue;
			}
			double GC = ((double)numGC)/((double)totalNonN); //fraction GC
			avgGC+=GC;

			totalPreselection++;


			if (numConsideredBgSeqs > 0) {
				curRebalanceCount++;
				int bindex = 0;
				for (int k=0;k<nBins;k++) {
					if (-1*GC < binThresholds[k]) {
						bindex = k;
						break;
					}
				}
				binTotals[bindex]+=2;

				if (curRebalanceCount > rateRebalanceCount) {
					curRebalanceCount=0;
					//time to rebalance rates
					long int t = 0;
					for (int k=0;k<nBins;k++) {
						t += binTotals[k];
					}	
					for (int k=0;k<nBins;k++) {
						binRates[k] = ((double)numConsideredBgSeqs/(((double)nBins))) /
									(((double)approxMaxSeqs)*(((double)binTotals[k])/((double)t)));
					}
				}

				if (((double)rand())/((double)RAND_MAX) < binRates[bindex]) {
					//add
				} else {
					continue;
				}
				binTotalselected[bindex]+=2;
			}

			avgGCselected += GC;

			BgSequence* s = new BgSequence();
			s->seq = &(ogSeqs[chrIndex][j]);
			s->len = seqLen;
			s->name = ogNames[chrIndex];
			s->pos = j;
			s->strand = MOTIF_STRAND_POS;
			s->stat = -1*GC;
			s->mask = 0;
			s->weight = 1.0;
			s->bin = -1;
			s->positionIndex=numStrandArray;
			posStrandArray[numStrandArray]=s;

			seqs[numSeqs++]=s;
			
			BgSequence* ss = new BgSequence();
			int rvPos = ogLens[chrIndex]-j-seqLen;
			ss->seq = &(ogSeqsRV[chrIndex][rvPos]);
			ss->len = seqLen;
			ss->name = ogNames[chrIndex];
			ss->pos = j;
			ss->strand = MOTIF_STRAND_NEG;
			ss->stat = -1*GC;
			ss->mask = 0;
			ss->weight = 1.0;
			ss->bin = -1;
			ss->positionIndex=numStrandArray;
			negStrandArray[numStrandArray]=ss;
			seqs[numSeqs++]=ss;

			numStrandArray++;

		}
	}

	double r = ((double)numRemovedDueToNs)/((double)totalPossible);
	fprintf(stderr, "\t\t%d of %ld (%.2lf%%) removed (fraction of Ns > %.3lf)\n", numRemovedDueToNs, totalPossible, r, maxNpercent);
	if (numSeqs < 1) {
		fprintf(stderr, "!!! Error, no sequences were processed!!!\n");
		return;
	}
	avgGC /= totalPreselection;
	fprintf(stderr, "\t\tProcessed %ld sequences from input sequences (avg. GC = %lf)\n", totalPreselection*2, avgGC);
	avgGCselected /= numStrandArray;
	if (numConsideredBgSeqs > 0) {
		fprintf(stderr, "\t\tSelected %ld sequences (avg. GC = %lf)\n", numSeqs, avgGCselected);
		fprintf(stderr, "\t\t\tBin\tThreshold\tCount\tCount Selected\tFinal Bin Rate\n");
		for (int i=0;i<nBins;i++) {
			if (i==nBins-1) fprintf(stderr, "\t\t\t%d\tmaximum  \t%ld\t%ld\t%lf\n",i+1,binTotals[i],binTotalselected[i],binRates[i]);
			else fprintf(stderr, "\t\t\t%d\t%lf\t%ld\t%ld\t%lf\n",i+1,binThresholds[i],binTotals[i],binTotalselected[i],binRates[i]);
		}
		delete []binTotals;
		delete []binTotalselected;
		delete []binRates;
	}
	fprintf(stderr, "\n");

}



void BgSequenceArray::maskGenomicPositionOverlaps(BgSequenceArray* targets, int level) {
	targets->sortBasedOnCoordinates();
	LongInttable* startCoords = new LongInttable(targets->numSeqs*5);
	char* curChr = NULL;
	int longest = 0;
	for (long int i=0;i<targets->numSeqs;i++) {
		char* c = targets->seqs[i]->name;
		if (targets->seqs[i]->len > longest) longest = targets->seqs[i]->len;
		if (curChr == NULL || strcmp(curChr,c) != 0) {
			curChr = c;
			startCoords->insert(i, c);
			//fprintf(stderr, "\t\t%s\t%d\n", c, i);
		}
	}

	long int curTargetIndex = -1;
	curChr = NULL;
	for (long int i=0;i<numStrandArray;i++) {
		char* chr = posStrandArray[i]->name;
		int p = posStrandArray[i]->pos;
		int len = posStrandArray[i]->len;

		if (curChr == NULL || strcmp(curChr,chr) != 0) {
			//new chromosome
			curChr = chr;
			curTargetIndex = startCoords->search(chr);
			if (curTargetIndex == EMPTY_INT) {
				curTargetIndex = -1;
			}
		}
		if (curTargetIndex < 0) continue;

		for (long int j=curTargetIndex;j<targets->numSeqs;j++) {
			if (strcmp(targets->seqs[j]->name,curChr) != 0) {
				curTargetIndex = -1;
				break;
			}
			int tp = targets->seqs[j]->pos;
			int tlen = targets->seqs[j]->len;
			if (tp+tlen < p) {
				curTargetIndex = j;
				continue;
			} else {
				break;
			}
		}
		if (curTargetIndex < 0) continue;
		for (long int j=curTargetIndex;j<targets->numSeqs;j++) {
			if (strcmp(targets->seqs[j]->name,curChr) != 0) {
				break;
			}
			int tp = targets->seqs[j]->pos;
			int tlen = targets->seqs[j]->len;
			if (tp+tlen < p) continue;
			if (tp <= p+len) {
				if (strand == MOTIF_STRAND_BOTH) {
					posStrandArray[i]->mask=level;
					negStrandArray[i]->mask=level;
				} else if (strand == MOTIF_STRAND_SEPARATE) {
					if (targets->seqs[j]->strand == MOTIF_STRAND_POS) {
						posStrandArray[i]->mask=level;
					} else if (targets->seqs[j]->strand == MOTIF_STRAND_NEG) {
						negStrandArray[i]->mask=level;
					}
				}
			}
			if (tp > p+len+longest) break;
		}
	}
	delete startCoords;
}
void BgSequenceArray::processSequencesAsIs(char strandSetting, double maxNpercent) {
	strand = strandSetting;
	seqLen = -1;

	double avgLen = 0.0;
	double avgGC = 0.0;
	int numRemovedDueToNs = 0;


	int approxMaxSeqs = numOGs;

	seqs = new BgSequence*[approxMaxSeqs];
	numSeqs = 0;

	for (int i=0;i<numOGs;i++) {
		int numNs = 0;
		int numGC = 0;
		int totalNonN = 0;
		for (int j=0;j<ogLens[i];j++) {
			if (ogSeqs[i][j] == 'N') numNs++;
			else {
				if (ogSeqs[i][j] == 'G' || ogSeqs[i][j] == 'C') numGC++;
				totalNonN++;
			}
		}
		if (((double)numNs)/((double)seqLen) > maxNpercent) {
			numRemovedDueToNs++;
			continue;
		}
		double GC = ((double)numGC)/((double)totalNonN); //fraction GC
		avgGC += GC;
		avgLen += ogLens[i];
		

		BgSequence* s = new BgSequence();
		s->seq = ogSeqs[i];
		s->len = ogLens[i];
		s->name = ogNames[i];
		s->pos = 0;
		s->strand = MOTIF_STRAND_POS;
		s->stat = -1*GC;
		s->weight = 1.0;
		s->bin = -1;
		s->mask = 0;
		s->positionIndex=-1;
		seqs[numSeqs++]=s;
//s->printAll(stderr);
	}

	fprintf(stderr, "\t%d of %d removed (fraction of Ns > %.3lf)\n", numRemovedDueToNs, numOGs, maxNpercent);
	if (numSeqs == 0) {
		fprintf(stderr, "\tError - no sequenced processed!!!\n");
		return;
	}
	avgLen /= (double)numSeqs;
	avgGC /= (double)numSeqs;
	seqLen = (int) avgLen;
	fprintf(stderr, "\t%ld sequences processed (avg length = %d, GC = %lf)\n", numSeqs, seqLen, avgGC);
	
}

void BgSequenceArray::sortBasedOnStat() {
	//fprintf(stderr, "sorting %ld sequences by the stat\n", numSeqs);
	qsort(seqs, numSeqs, sizeof(BgSequence*), &(BgSequenceArray::statCmp));
}
void BgSequenceArray::sortBasedOnCoordinates() {
	qsort(seqs, numSeqs, sizeof(BgSequence*), &(BgSequenceArray::coordCmp));
}
void BgSequenceArray::sortBasedOnStat(BgSequence** array, long int len) {
	qsort(array, len, sizeof(BgSequence*), &(BgSequenceArray::statCmp));
}
int BgSequenceArray::statCmp(const void *a, const void *b) {
	double va=(*((BgSequence**)a))->stat;
	double vb=(*((BgSequence**)b))->stat;
	if (va < vb) return -1;
	if (va > vb) return 1;
	return 0;
}
int BgSequenceArray::coordCmp(const void *a, const void *b) {
	char* va=(*((BgSequence**)a))->name;
	char* vb=(*((BgSequence**)b))->name;
	int c = strcmp(va, vb);
	if (c != 0) {
		return c;
	}
	int va2=(*((BgSequence**)a))->pos+(*((BgSequence**)a))->len;
	int vb2=(*((BgSequence**)b))->pos+(*((BgSequence**)b))->len;
	if (va2 < vb2) return -1;
	if (va2 > vb2) return 1;
	return 0;
}

int cmpBgPositions(const void *a, const void *b) {
	char* va=((BgPosition*)a)->chr;
	char* vb=((BgPosition*)b)->chr;
	int c = strcmp(va, vb);
	if (c != 0) {
		return c;
	}
	int va2 =((BgPosition*)a)->start;
	int vb2 =((BgPosition*)b)->start;
	if (va2 < vb2) return -1;
	if (va2 > vb2) return 1;
	return 0;
}

double* BgSequenceArray::generateBinThresholds(int numBins) {
	long int seqPerBin = (long int) ceil(((double)numSeqs)/((double)numBins));
	double* binThresholds = new double[numBins];
	fprintf(stderr, "\tEstablishing sequence bins (n=%d, %ld seqs per bin):\n", numBins, seqPerBin);
	fprintf(stderr, "\t\tSorting sequences by GC content:\n");
	sortBasedOnStat();

	fprintf(stderr, "\t\tBin\tStat Threshold\n");
	for (int i=0;i<numBins;i++) {
		long int index = (i+1)*seqPerBin;
		if (index >= numSeqs) index = numSeqs-1;
		binThresholds[i] = seqs[index]->stat;
		fprintf(stderr, "\t\t%d\t%lf\n", i+1, binThresholds[i]);
		if (i > 0 && binThresholds[i] == binThresholds[i-1]) {
			fprintf(stderr, "!!! Error - bin thresholds were the same - consider reducing the number of bins (-nbins)\n");
			fprintf(stderr, "!!! (Another common problem is that the input target regions are too short, e.g. 1bp long)\n");
			exit(0);
		}
	}
	binThresholds[numBins-1]=FLT_MAX;
	return binThresholds;
}
void BgSequenceArray::assignSeqsToBins(double* binThresholds,int numBINs) {
	numBins = numBINs;
	bins = new BgSequence**[numBins];
	binSizes = new long int[numBins];
	for (int i=0;i<numBins;i++) {
		bins[i] = NULL;
		binSizes[i] = 0;
	}
	for (long int i=0;i<numSeqs;i++) {
		double stat = seqs[i]->stat;
		for (int j=0;j<numBins;j++) {
			if (stat <= binThresholds[j]) {
				binSizes[j]++;
				break;
			}
		}
	}
	fprintf(stderr, "\t\tBin\tnumber of sequences assigned\n");
	for (int i=0;i<numBins;i++){ 
		fprintf(stderr, "\t\t%d\t%ld\n", i+1, binSizes[i]);
		bins[i] = new BgSequence*[binSizes[i]];
		binSizes[i] = 0;
	}
	fprintf(stderr, "\n");
	for (long int i=0;i<numSeqs;i++) {
		double stat = seqs[i]->stat;
		for (int j=0;j<numBins;j++) {
			if (stat <= binThresholds[j]) {
				bins[j][binSizes[j]] = seqs[i];
				seqs[i]->bin = j;
				binSizes[j]++;
				break;
			}
		}
	}
}
void BgSequenceArray::selectInitialSequences(long int numBGSEQs) {
	deleteSBins();
	clearMask();

	numBgSeqs = numBGSEQs;
	numBgSeqPerBin = (long int) ceil(((double)numBgSeqs)/((double)numBins));

	sbins = new BgSequence**[numBins];
	sbinSizes = new long int[numBins];
	sbinFixed = new char[numBins];
	fprintf(stderr, "\tSelecting sequences by bin:\n");
	for (int i=0;i<numBins;i++) {
		sbinSizes[i] = 0;
		sbinFixed[i] = 0;
		long int max = numBgSeqPerBin;
		if (numBgSeqPerBin > binSizes[i]) {
			//fprintf(stderr, "!!! Warning: Less sequences in bin %d (%d) than needed (%d).\n", i+1,binSizes[i],numBgSeqPerBin);
			max = binSizes[i];
		}
		sbins[i] = new BgSequence*[max];

		for (long int j=0;j<max;j++) {
			int fail = 1;
			while (fail) {
				long int index = rand() % binSizes[i];
				if (bins[i][index]->mask) {
					fail++;
					if (fail > maxFailures) break;
					continue;
				}
				sbins[i][sbinSizes[i]] = bins[i][index];
				mask(bins[i][index],1);
				sbinSizes[i]++;
				fail = 0;
			}
			if (fail > maxFailures) {
				fprintf(stderr, "\t\t\t!!! Reached maxFailures (%d)\n", maxFailures);
				break;
			}
		}
		fprintf(stderr, "\t\tBin=%d\tn=%ld\n", i+1, sbinSizes[i]);
	}
}
void BgSequenceArray::deleteSBins() {
	if (sbins != NULL) {
		for (int i=0;i<numBins;i++) {
			if (sbins[i] != NULL) delete []sbins[i];
		}
		delete []sbins;
		sbins = NULL;
	}
	if (sbinSizes != NULL) delete []sbinSizes;
	sbinSizes = NULL;
}
void BgSequenceArray::clearMask() {
	for (long int i=0;i<numSeqs;i++) {
		if (seqs[i]->mask < 2) seqs[i]->mask = 0;
	}
}
void BgSequenceArray::resetMask() {
	for (int i=0;i<numBins;i++) {
		for (long int j=0;j<sbinSizes[i];j++) {
			mask(sbins[i][j],1);
		}
	}
}
void BgSequenceArray::resetMask(int bin) {
	for (long int j=0;j<sbinSizes[bin];j++) {
		mask(sbins[bin][j],1);
	}
}
//Mask level 0, 1, 2(perm)
void BgSequenceArray::mask(BgSequence* bgSeq,int level) {
	if (bgSeq->mask < level) bgSeq->mask = level;
	if (overlapFlag==0) return;
		
	int L = bgSeq->len;
	long int pindex = bgSeq->positionIndex;
	if (pindex < 0) {
		// there is no 'positioning' to avoid overlap
		return;		
	}
	char curStrand = bgSeq->strand;
	for (long int i=pindex-L;i<=pindex+L;i++) {
		if (i<0) continue;
		if (i>=numStrandArray) break;
		if (strand == MOTIF_STRAND_BOTH || (curStrand == MOTIF_STRAND_POS)) {
			if (strcmp(posStrandArray[i]->name,bgSeq->name)==0) {
				int diff = posStrandArray[i]->pos - bgSeq->pos;
				if (diff < 0) diff *= -1;
				if (diff <= L) {
					if (posStrandArray[i]->mask < level) posStrandArray[i]->mask=level;
				}
			}
		}
		if (strand == MOTIF_STRAND_BOTH || (curStrand == MOTIF_STRAND_NEG)) {
			if (strcmp(negStrandArray[i]->name,bgSeq->name)==0) {
				int diff = negStrandArray[i]->pos - bgSeq->pos;
				if (diff < 0) diff *= -1;
				if (diff <= L) {
					if (negStrandArray[i]->mask < level) negStrandArray[i]->mask=level;
				}
			}
		}	
	}
}

void BgSequenceArray::optimizePositionalFrequncies(int maxIterations, double changeRate) {

	fprintf(stderr, "\tOptimizing positional kmer frequencies (kmer=%d)\n", pkmer);
	for (int i=0;i<maxIterations;i++) {
		double totalError = 0.0;
		for (int j=0;j<numBins;j++) {
			totalError += optimizePositionalFrequencyBin(j,changeRate);
		}
		clearMask();
		resetMask();
		fprintf(stderr, "\t\tIteration: %d\tError=%lf\n", i+1, totalError);
	}

}
double BgSequenceArray::optimizePositionalFrequencyBin(int bin,double changeRate) {
	//need to get target frequencies...
	double error = 0.0;

	double **targetFreq = targetFreqsByBin[bin];
	double **bgFreq = calculatePositionalKmerFrequencies(sbins[bin], sbinSizes[bin]);

	//fprintf(stderr, "Target:\n");
	//printPositionalKmerFrequencies(stderr, targetFreq);
	//fprintf(stderr, "Bg:\n");
	//printPositionalKmerFrequencies(stderr, bgFreq);

	double maxError = 0.0;
	for (int i=0;i<seqLen;i++) {
		double maxErrorPos = 0.0;
		for (int j=0;j<freqLength;j++) {
			double d = targetFreq[i][j]-bgFreq[i][j];
			bgFreq[i][j] = d;
			double fabsD = fabs(d);
			error += fabsD;
			if (fabsD > maxErrorPos) maxErrorPos = fabsD;
			//store the error back in d to check on later
		}
		maxError += maxErrorPos;
	}
	//fprintf(stderr, "Diff:\n");
	//printPositionalKmerFrequencies(stderr, bgFreq);
	fprintf(stderr, "Error: %lf\n\n",error);

	double *matchScores = new double[sbinSizes[bin]];
	double minMatch = 0;
	double maxMatch = 0;
	for (int i=0;i<sbinSizes[bin];i++) {
		double match = 0.0;
		BgSequence* bg = sbins[bin][i];
		for (int j=0;j<bg->len-pkmer+1;j++) {
			int index = seq2Index(&(bg->seq[j]),pkmer);
			if (index < 0) continue;
			match += bgFreq[j][index];
		}
		matchScores[i] = match;
		if (match>maxMatch)maxMatch=match;
		if (match<minMatch)minMatch=match;
	}
	//fprintf(stderr, "\tMinMatch: %lf\n", minMatch);
	//fprintf(stderr, "\tMaxMatch: %lf\n", maxMatch);
	double matchRange = maxMatch-minMatch;
	if (matchRange < 1e-10) matchRange = 1e-10;

	int totalReplaced = 0;
	for (int i=0;i<sbinSizes[bin];i++) {
		
		double r = (maxMatch-matchScores[i])/matchRange; // * changeRate;
		r *= r;

		//char* s = sbins[bin][i]->getSeq();
		//fprintf(stderr, "\t%s\t%lf\t%lf\n", s, matchScores[i],r);
		//delete []s;

		if (matchScores[i] > 0.0) continue;
		if (((double)rand()/(double)RAND_MAX) < r) {
			// replace 
	
			int fail = 1;
			while (fail) {
				int index = rand() % binSizes[bin];
				if (bins[bin][index]->mask) {
					fail++;
					if (fail > maxFailures) break;
					continue;
				}
				sbins[bin][i] = bins[bin][index];
				mask(bins[bin][index],1);
				fail = 0;
			}
			if (fail > maxFailures) {
				fprintf(stderr, "!!! Reached maxFailures (%d)\n", maxFailures);
				break;
			} else {
				totalReplaced++;
			}
		}
	}
	fprintf(stderr, "\ttotal replaced = %d\n", totalReplaced);

	return error;
}
void BgSequenceArray::setpKmer(int PKMER,int posFlag) {
	pkmer = PKMER;
	positionFlag = posFlag;
	freqLength = 1;
	for (int i=0;i<pkmer;i++) {
		freqLength *= alphaLen;
	}
}

double*** BgSequenceArray::calculatePositionalKmerFrequenciesByBin() {
	double*** freq = new double**[numBins];
	for (int i=0;i<numBins;i++) {
		freq[i] = calculatePositionalKmerFrequencies(bins[i], binSizes[i]);
	}
	return freq;
}
double** BgSequenceArray::calculatePositionalKmerFrequencies(BgSequence** seqArray, long int num) {

	double **freq = new double*[seqLen];
	double *N = new double[seqLen];
	for (int i=0;i<seqLen;i++) {
		freq[i] = new double[freqLength];
		N[i] = 0.0;
		for (int j=0;j<freqLength;j++) {
			freq[i][j]=0.0;
		}
	}
	for (long int i=0;i<num;i++) {
		BgSequence* bg = seqArray[i];
		for (int j=0;j<bg->len-pkmer+1;j++) {
			int index = seq2Index(&(bg->seq[j]),pkmer);
			if (index >= 0) {
				if (positionFlag) {
					freq[j][index]+=1.0;
					N[j]+=1.0;
				} else {
					freq[0][index]+=1.0;
					N[0]+=1.0;
				}
			}
		}
	}
	for (int i=0;i<seqLen;i++) {
		for (int j=0;j<freqLength;j++) {
			if (N[i] > 0.0) {
				freq[i][j] /= N[i];
			}
		}
	}
	delete []N;
	return freq;
}
double** BgSequenceArray::calculatePositionalKmerFrequencies(BgSequence*** seqArray, long int* sizes, int num) {

	double **freq = new double*[seqLen];
	double *N = new double[seqLen];
	for (int i=0;i<seqLen;i++) {
		freq[i] = new double[freqLength];
		N[i] = 0.0;
		for (int j=0;j<freqLength;j++) {
			freq[i][j]=0.0;
		}
	}
	for (int bin=0;bin<num;bin++) {
		for (long int i=0;i<sizes[bin];i++) {
			BgSequence* bg = seqArray[bin][i];
			for (int j=0;j<bg->len-pkmer+1;j++) {
				int index = seq2Index(&(bg->seq[j]),pkmer);
				if (index >= 0) {
					if (positionFlag) {
						freq[j][index]+=1.0;
						N[j]+=1.0;
					} else {
						freq[0][index]+=1.0;
						N[0]+=1.0;
					}
				}
			}
		}
	}
	for (int i=0;i<seqLen;i++) {
		for (int j=0;j<freqLength;j++) {
			if (N[i] > 0.0) {
				freq[i][j] /= N[i];
			}
		}
	}
	
	delete []N;
	return freq;
}

void BgSequenceArray::selectAllSequences() {
	deleteSBins();
	sbins = new BgSequence**[numBins];
	sbinSizes = new long int[numBins];
	sbinFixed = new char[numBins];
	for (int i=0;i<numBins;i++) {
		sbinSizes[i] = binSizes[i];
		sbinFixed[i] = 0;
		sbins[i] = new BgSequence*[sbinSizes[i]];
		for (long int j=0;j<sbinSizes[i];j++) {
			sbins[i][j] = bins[i][j];
		}
	}
}
void BgSequenceArray::scoreKmerProfiles() {
	double minKmerFreq = 1e-6;
	for (int i=0;i<numBins;i++) {
		double* freqAvg = new double[freqLength];
		for (int j=0;j<freqLength;j++) {
			freqAvg[j] = 0.0;
			double n = 1.0;
			if (positionFlag) {
				for (int k=0;k<seqLen;k++) {
					freqAvg[j] += targetFreqsByBin[i][k][j];	
					n += 1.0;
				}
			} else {
				freqAvg[j] += targetFreqsByBin[i][0][j];	
				n += 1.0;
			}
			if (n < 1.0) n=1.0;
			freqAvg[j] /= n;
		}
		for (long int k=0;k<binSizes[i];k++) {
			BgSequence* bg = bins[i][k];
			double s = 0.0;
			for (int j=0;j<bg->len-pkmer+1;j++) {
				int index = seq2Index(&(bg->seq[j]),pkmer);
				double v = 0.0;
				if (index >= 0) {
					if (positionFlag) {
						v = targetFreqsByBin[i][j][index];
					} else {
						v = targetFreqsByBin[i][0][index];
					}
				}
				if (v < minKmerFreq) v = minKmerFreq;
				if (positionFlag) {
					double vv = freqAvg[index];
					if (vv < minKmerFreq) vv = minKmerFreq;
					v /= vv;
				}
				s += log(v)/log(2.0);
			}
			bg->pscore = s;
		}
	}
}

void matchKmerFrequencies(BgSequenceArray* targets, BgSequenceArray* bg, long int numBgSeqs, int pkmer, int positionFlag, int numSubBins, int numIters,
							int overlapIteration, double decayRate, char* pscoreOutputFile) {
	bg->deleteSBins();
	bg->clearMask();
	bg->numBgSeqs = numBgSeqs;
	bg->numBgSeqPerBin = (long int) ceil(((double)numBgSeqs)/((double)bg->numBins));
	long int numBgSeqPerBin = bg->numBgSeqPerBin;
	long int numBgSeqPerBinTarget = bg->numBgSeqPerBin * 2;

	targets->numSubBins = numSubBins;
	bg->numSubBins = numSubBins;

	targets->targetFreqsByBin = targets->calculatePositionalKmerFrequenciesByBin();
	bg->targetFreqsByBin = targets->targetFreqsByBin;
	bg->selectAllSequences();

	// Idea is to start by allowing overlap to get the best sequences for positional matching
	// Then for the 2nd half of iterations enforce regions not to overlap
	int ogOverlapFlag = bg->overlapFlag;
	bg->overlapFlag=0;
	for (int z=0;z<numIters;z++) {
		fprintf(stderr, "\t\tStarting iteration %d of %d", z+1,numIters);
		if (ogOverlapFlag > 0 && z+1 >= overlapIteration) {
			fprintf(stderr, ", restricting overlap\n");
			bg->overlapFlag = 1;
		}
		fprintf(stderr, "\n");
		bg->clearMask();

		//set masked from the fixed bins first
		for (int i=0;i<targets->numBins;i++) {
			if (bg->sbinFixed[i]) {
				bg->resetMask(i);		
			}
		}
		int updateFlag = 0;
		for (int i=0;i<targets->numBins;i++) {

			if (bg->sbinFixed[i]) {
				continue;
			}
			updateFlag = 1;
			//long int currentCount = bg->binSizes[i];
			long int numSeqToSelect = numBgSeqPerBinTarget;
			if (z==0) {
				numSeqToSelect = ((long int)((bg->binSizes[i]-numBgSeqPerBinTarget)*decayRate))+numBgSeqPerBinTarget;
			} else {
				//currentCount = bg->sbinSizes[i];
				if (z==numIters-1) {
					numSeqToSelect = numBgSeqPerBinTarget;
				} else {
					numSeqToSelect = ((long int)((bg->sbinSizes[i]-numBgSeqPerBinTarget)*decayRate))+numBgSeqPerBinTarget;
				}
			}
			if (bg->sbinSizes[i] < numBgSeqPerBinTarget) {
				numSeqToSelect = bg->sbinSizes[i];
			}
			//fprintf(stderr, "\t\t\tSelecting from %ld -> %ld (%lf)\n", currentCount, numSeqToSelect, decayRate);
				

			double** targetFreqs = bg->targetFreqsByBin[i];		
			double** bgFreqs = bg->calculatePositionalKmerFrequencies(bg->sbins[i],bg->sbinSizes[i]);
			//double** bgFreqs = bg->calculatePositionalKmerFrequencies(sbins[i],sbinSizes[i]);
			double totalDiff = 0.0;
			int flen = 1;
			if (positionFlag) flen = bg->seqLen;
			for (int j=0;j<flen;j++) {
				for (int k=0;k<bg->freqLength;k++) {
					//double d = targetFreqs[j][k]-bgFreqs[j][k];
					double d = bgFreqs[j][k]-targetFreqs[j][k];
					bgFreqs[j][k] = d;
					totalDiff += fabs(d);
				}
			}
			double** freqDiff = bgFreqs;
			//bgFreq now has the differences between target and bg frequencies encoded for each position
			//Will now use that to score sequences in the target to get their distribution, then try to 
			//match that by selecting a similar distribution from the background.
			targets->scoreSequencesByKmerFreqDifference(targets->bins[i],targets->binSizes[i],freqDiff);
			double* thresholds = targets->generateSubBinThresholds(i);
			bg->scoreSequencesByKmerFreqDifference(bg->sbins[i],bg->sbinSizes[i],freqDiff);
			bg->selectSeqsBasedOnSubBins(i,thresholds,numSeqToSelect);

			if (bg->sbinSizes[i] < numBgSeqPerBinTarget) {
				if (bg->overlapFlag > 0 || ogOverlapFlag == 0) {
				//fprintf(stderr, "\t\t\t\t%d < %d (numBgSeqPerBinTarget)\n",bg->sbinSizes[i],numBgSeqPerBinTarget);
					bg->sbinFixed[i] = 1;
				}
			}
	
			//double** selectFreqs = bg->calculatePositionalKmerFrequencies(bg->sbins[i],bg->sbinSizes[i]);
	
			for (int j=0;j<bg->seqLen;j++) {
				delete []bgFreqs[j];
			}
			delete []bgFreqs;
			delete []thresholds;
			fprintf(stderr, "\t\t\tBin=%d numSeqs=%ld fixed=%d kmerDiff=%lf\n", i+1, bg->sbinSizes[i], bg->sbinFixed[i],totalDiff);
		}

		if (z==0) {
			//On the first iteration, store the profile score to estimate how well the sequence matches the positional freqs
			for (long int i=0;i<targets->numSeqs;i++) targets->seqs[i]->pscore = targets->seqs[i]->stat;
			for (long int i=0;i<bg->numSeqs;i++) bg->seqs[i]->pscore = bg->seqs[i]->stat;

 			if (pscoreOutputFile != NULL) {
				//output scores for all background - good if you want to see everything
				fprintf(stderr, "\n\tOutputing profile scores to file %s\n", pscoreOutputFile);
				FILE* fp = fopen(pscoreOutputFile,"w");
				if (fp == NULL) {
					fprintf(stderr, "!!! Error, could not open %s to output profile scores\n", pscoreOutputFile);
					exit(0);
				}
				targets->printAllSequencesTable(fp);
				bg->printAllSequencesTable(fp);
				fclose(fp);
			}
		}

		if (updateFlag == 0) {
			break;
		}
	}
	bg->equalizesBins(numBgSeqPerBin);
}
void BgSequenceArray::equalizesBins(long int numPerBin) {
	long int min = numSeqs;
	for (int i=0;i<numBins;i++) {
		if (sbinSizes[i] < min) min = sbinSizes[i];
	}
	fprintf(stderr, "\t\tEqualizing selected sequence counts by bin (requested=%ld, min=%ld)\n", numPerBin,min);
	if (min > numPerBin) min = numPerBin;
	fprintf(stderr, "\t\t\tBin\tOld Count\tNew Count\n");
	for (int i=0;i<numBins;i++) {
		long int curCount = sbinSizes[i];
		if (curCount > min) {
			for (long int j=0;j<sbinSizes[i];j++) {
				sbins[i][j]->stat = (double)(rand() % sbinSizes[i]);
			}
			qsort(sbins[i],sbinSizes[i],sizeof(BgSequence*),&(BgSequenceArray::statCmp));
			BgSequence** newset = new BgSequence*[min];
			for (long int j=0;j<min;j++) {
				newset[j] = sbins[i][j]; 
			}
			delete []sbins[i];
			sbins[i] = newset;
			sbinSizes[i] = min;
		}
		fprintf(stderr, "\t\t\t%d\t%ld\t%ld\n", i+1, curCount, min);
	}
}
void BgSequenceArray::scoreSequencesByKmerFreqDifference(BgSequence** seqSet,long int numSeqSet,double **freqDiff) {
	for (long int i=0;i<numSeqSet;i++) {
		double match = 0.0;
		BgSequence* s = seqSet[i];
		for (int j=0;j<s->len-pkmer+1;j++) {
			int index = seq2Index(&(s->seq[j]),pkmer);
			if (index < 0) continue;
			if (positionFlag) {
				match += freqDiff[j][index];
			} else {
				match += freqDiff[0][index];
			}
		}
		s->stat = match;
	}
}

double* BgSequenceArray::generateSubBinThresholds(int bin) {
	sortBasedOnStat(bins[bin],binSizes[bin]);
	int seqPerBin = (int) ceil(((double)binSizes[bin])/((double)numSubBins));
	double* binThresholds = new double[numSubBins];

	//fprintf(stderr, "\tEstablishing sequence subBins (%d seqs per subBin):\n", seqPerBin);
	for (int i=0;i<numSubBins;i++) {
		long int index = (i+1)*seqPerBin;
		if (index >= binSizes[bin]) index = binSizes[bin]-1;
		binThresholds[i] = bins[bin][index]->stat;
		//fprintf(stderr, "\t\t%d\t%lf\n", i+1, binThresholds[i]);
		if (i > 0 && binThresholds[i] == binThresholds[i-1]) {
			fprintf(stderr, "!!! Error - sub bin thresholds were the same - consider reducing the number of sub bins (-nsubBins)\n");
			for (int k=0;k<=i;k++) {
				fprintf(stderr, "\t\t%d\t%d\t%lf\n", i, k, binThresholds[k]);
			}
			if (i>1) {
				binThresholds[i-1] = (binThresholds[i-1]+binThresholds[i-2])/2.0;
				fprintf(stderr, "!!! Trying to recover...\n");
			} else {
				fprintf(stderr, "!!! Probably not recoverable, reduce -nsubBins\n");
				exit(0);
			}
			//exit(0);
		}
	}
	binThresholds[numSubBins-1]=FLT_MAX;
	return binThresholds;
}

void BgSequenceArray::selectSeqsBasedOnSubBins(int bin,double* thresholds, long int numSeqToSelect) {

	BgSequence** curBin = sbins[bin];
	long int curSize = sbinSizes[bin];

	sbins[bin] = new BgSequence*[numSeqToSelect];
	sbinSizes[bin] = 0;

	//numBgSeqs = numBGSEQs;
	//numBgSeqPerBin = (int) ceil(((double)numBgSeqs)/((double)numBins));
	long int numSeqPerSubBin = (long int) ceil(((double)numSeqToSelect)/((double)numSubBins));
	//fprintf(stderr, "numSeqPerSubBin=%ld\n", numSeqPerSubBin);

	BgSequence*** subBins = new BgSequence**[numSubBins];
	long int* subBinSizes = new long int[numSubBins];
	for (int i=0;i<numSubBins;i++) {
		subBins[i] = NULL;
		subBinSizes[i] = 0;
	}

	for (long int i=0;i<curSize;i++) {
		double stat = curBin[i]->stat;
		for (int j=0;j<numSubBins;j++) {
			if (stat <= thresholds[j]) {
				subBinSizes[j]++;
				break;
			}
		}
	}
	//fprintf(stderr, "\tNumber of sequences assigned to each subBin:\n");
	for (int i=0;i<numSubBins;i++){ 
		//fprintf(stderr, "\t\t%d\t%d\n", i, subBinSizes[i]);
		subBins[i] = new BgSequence*[subBinSizes[i]];
		subBinSizes[i] = 0;
	}
	for (long int i=0;i<curSize;i++) {
		double stat = curBin[i]->stat;
		for (int j=0;j<numSubBins;j++) {
			if (stat <= thresholds[j]) {
				subBins[j][subBinSizes[j]] = curBin[i];
				subBinSizes[j]++;
				break;
			}
		}
	}

	//fprintf(stderr, "\tSelecting sequences by sub bin:\n");
	for (int i=0;i<numSubBins;i++) {
		long int currSubBinTotal = 0;
		long int max = numSeqPerSubBin;
		if (numSeqPerSubBin > subBinSizes[i]) {
			//fprintf(stderr, "!!! Warning: Less sequences in subBin %d (%d) than needed (%d).\n", i+1,subBinSizes[i],numSeqPerSubBin);
			max = subBinSizes[i];
		}
		for (int j=0;j<max;j++) {
			int fail = 1;
			while (fail && sbinSizes[bin] < numSeqToSelect) {
				long int index = rand() % subBinSizes[i];
				if (subBins[i][index]->mask) {
					fail++;
					if (fail > maxFailures) break;
					continue;
				}
				sbins[bin][sbinSizes[bin]] = subBins[i][index];
				sbinSizes[bin]++;
				mask(subBins[i][index],1);
				currSubBinTotal++;
				fail = 0;
			}
			if (sbinSizes[bin] >= numSeqToSelect) {
				//fprintf(stderr, "=== sbinSizes[%d] %d >= %d\n", bin, sbinSizes[bin], numSeqToSelect);
				break;
			}
			if (fail > maxFailures) {
				//fprintf(stderr, "!!! Reached maxFailures (%d)\n", maxFailures);
				break;
			}
		}
		//fprintf(stderr, "\t\tsubBin=%d\tn=%ld\tmax=%ld\n", i+1, currSubBinTotal,max);

		delete []subBins[i];
	}

	delete []curBin;
	delete []subBins;
	delete []subBinSizes;

}




int BgSequenceArray::seq2Index(char* s, int len) {
	int index = 0;
	int scale = 1;
	for (int i=0;i<len;i++) {
		int x = 0;
		if (s[i] == 'N') {
			if (alphaLen == 5) x=4;
			else return -1;
		} else if (s[i] == 'A') {
			x=0;
		} else if (s[i] == 'C') {
			x=1;
		} else if (s[i] == 'G') {
			x=2;
		} else if (s[i] == 'T') {
			x=3;
		} else {
			return -1;
		}
		index += x * scale;
		scale *= alphaLen;
	}
	return index;
}

// creates new memory
char* BgSequenceArray::index2Seq(int index, int len) {

	char* rv = new char[len+1];
	rv[0] = '\0';
	rv[len] = '\0';
	int scale = alphaLen;
	for (int i=0;i<len;i++) {
		int x = index % scale;
		x /= (scale/alphaLen);
		if (x == 0) {
			rv[i] = 'A';
		} else if (x==1) {
			rv[i] = 'C';
		} else if (x==2) {
			rv[i] = 'G';
		} else if (x==3) {
			rv[i] = 'T';
		} else if (x==4) {
			rv[i] = 'N';
		} else {
			fprintf(stderr, "\tError something wrong with index2Seq - x=%d (input index=%d, len=%d)\n", x, index, len);
		}
		scale *= alphaLen;
	}
	return rv;
}


long int BgSequenceArray::getTotalSelectedSequences() {
	long int total = 0;
	if (sbinSizes != NULL) {
		for (int i=0;i<numBins;i++){ 
			total += sbinSizes[i];
		}
	}
	return total;
}


void BgSequenceArray::printPositionalKmerFrequencies(FILE* fp, double** freq) {
	fprintf(fp, "Position");
	for (int j=0;j<freqLength;j++) {
		char* kmer = index2Seq(j, pkmer);
		fprintf(fp, "\t%s", kmer);
		delete []kmer;
	}
	fprintf(fp, "\n");
	
	if (positionFlag) {	
		for (int i=0;i<seqLen;i++) {
			fprintf(fp, "%d",i);
			for (int j=0;j<freqLength;j++) {
				fprintf(fp, "\t%lf",freq[i][j]);
			}
			fprintf(fp,"\n");
		}
	} else {
		fprintf(fp, "PositionIndependent");
		for (int j=0;j<freqLength;j++) {
			fprintf(fp, "\t%lf",freq[0][j]);
		}
		fprintf(fp,"\n");
	}
}

void BgSequenceArray::printAllSequencesFASTA(FILE* fp) {
	int n = 0;
	for (long int i=0;i<numSeqs;i++) {
		seqs[i]->printFASTA(fp,n++);
	}
}
void BgSequenceArray::printAllSequencesBED(FILE* fp) {
	int n = 0;
	for (long int i=0;i<numSeqs;i++) {
		seqs[i]->printBED(fp,n++);
	}
}
void BgSequenceArray::printAllSequencesTable(FILE* fp) {
	fprintf(fp,"SeqID\tchr\tposition\tstrand\tweight\tProfileScore\tGC\tbin\n");
	int n = 0;
	for (long int i=0;i<numSeqs;i++) {
		seqs[i]->printTable(fp,n++);
	}
}
void BgSequenceArray::printSelectedSequencesFASTA(FILE* fp) {
	int n=0;
	for (int i=0;i<numBins;i++) {
		for (long int j=0;j<sbinSizes[i];j++) {
			sbins[i][j]->printFASTA(fp,n++);
		}
	}
}
void BgSequenceArray::printSelectedSequencesBED(FILE* fp) {
	int n = 0;
	for (int i=0;i<numBins;i++) {
		for (long int j=0;j<sbinSizes[i];j++) {
			sbins[i][j]->printBED(fp,n++);
		}
	}
}
void BgSequenceArray::printSelectedSequencesTable(FILE* fp) {
	fprintf(fp,"SeqID\tchr\tposition\tstrand\tweight\tProfileScore\tGC\tbin\n");
	int n= 0;
	for (int i=0;i<numBins;i++) {
		for (long int j=0;j<sbinSizes[i];j++) {
			sbins[i][j]->printTable(fp,n++);
		}
	}
}
void BgSequenceArray::printGroupFile(FILE* fp,int assignment) {
	int n=0;
	for (int i=0;i<numBins;i++) {
		for (long int j=0;j<sbinSizes[i];j++) {
			sbins[i][j]->printGroup(fp,assignment,n++);
		}
	}
}
void BgSequenceArray::printSeqFile(FILE* fp) {
	int n=0;
	for (int i=0;i<numBins;i++) {
		for (long int j=0;j<sbinSizes[i];j++) {
			sbins[i][j]->printSeq(fp,n++);
		}
	}
}

BgModel::BgModel() {
	freq = NULL;
	fwdFreq = NULL;
	revFreq = NULL;
	fwdIndex = NULL;
	revIndex = NULL;
	freqLength = 0;
	seqLen = 0;
	alphaLen = 5;
	positionFlag =0;
	alphaRef = NULL;
	freqRef = NULL;
	kmer = 0;
}
BgModel::~BgModel() {
}
void BgModel::initializeModel(int mer, int alphaLength, int pFlag, int slen) {
	kmer = mer;
	alphaLen = alphaLength;
	positionFlag = pFlag;
	seqLen = slen;
	modelLen = seqLen-kmer+1;
	if (positionFlag==0) {
		modelLen = 1;
	}

	initializeFreqIndex();

	freq = new double*[modelLen];
	counts = new double[modelLen];
	fwdFreq = new double**[modelLen];
	fwdCounts = new double*[modelLen];
	revFreq = new double**[modelLen];
	revCounts = new double*[modelLen];
	for (int i=0;i<modelLen;i++) {
		freq[i] = new double[freqLength];
		counts[i] = 0.0;
		fwdFreq[i] = new double*[freqLength];
		fwdCounts[i] = new double[freqLength];
		revFreq[i] = new double*[freqLength];
		revCounts[i] = new double[freqLength];
		for (int j=0;j<freqLength;j++) {
			freq[i][j] = 0.0;
			fwdFreq[i][j] = new double[alphaLen];
			fwdCounts[i][j] = 0.0;
			revFreq[i][j] = new double[alphaLen];
			revCounts[i][j] = 0.0;
			for (int k=0;k<alphaLen;k++) {
				fwdFreq[i][j][k] = 0.0;
				revFreq[i][j][k] = 0.0;
			}
		}
	}
}
void BgModel::initializeFreqIndex() {
	freqLength = 1;
	for (int i=0;i<kmer;i++) {
		freqLength *= alphaLen;
	}

	fwdIndex = new int*[freqLength];
	revIndex = new int*[freqLength];
	for (int i=0;i<freqLength;i++) {
		fwdIndex[i] = new int[alphaLen];
		revIndex[i] = new int[alphaLen];
		for (int j=0;j<alphaLen;j++) {
			fwdIndex[i][j] = -1;
			revIndex[i][j] = -1;
		}
	}
	char* s = new char[kmer+3];
	char* ss = &(s[1]);
	alphaRef = "ACGTN";
	freqRef = new char*[freqLength];

	s[kmer+2]='\0';
	for (int i=0;i<freqLength;i++) {
		index2Seq(i,kmer,ss);
		freqRef[i] = new char[kmer+1];
		ss[kmer]='\0';
		strcpy(freqRef[i], ss);
		for (int j=0;j<alphaLen;j++) {
			s[0] = alphaRef[j];
			s[kmer+1] = alphaRef[j];
			//fprintf(stderr, "\t%d\t%d\t%s\n",i,j,s);

			fwdIndex[i][j] = seq2Index(&(s[2]),kmer);
			revIndex[i][j] = seq2Index(s,kmer);
		}
	}
	delete []s;
}
void BgModel::buildModel(BgSequence** seqs, long int numSeqs) {

	for (long int i=0;i<numSeqs;i++) {
		char* seq = seqs[i]->seq;
		int len = seqs[i]->len;
		int lastIndex = -1;
		for (int j=0;j<len-kmer+1;j++) {
			int curIndex = seq2Index(&(seq[j]),kmer);
			if (curIndex < 0) {
				lastIndex = -1;
				continue;
			}

			int mIndex = j;
			if (positionFlag == 0) {
				mIndex = 0;
			} else {
				if (j >= modelLen) {
					//sequence is longer than sequenxce specific model
					fprintf(stderr, "!!! Warning: sequence is longer than sequenxce specific model (%d vs. %d)\n",len, modelLen);
					break;
				}
			}

			freq[mIndex][curIndex] += 1.0;
			counts[mIndex] += 1.0;
			if (lastIndex > -1) {
				int x = -1;
				char c = seq[j-1];
				if (c == 'A') {
					x=0;
				} else if (c == 'C') {
					x=1;
				} else if (c == 'G') {
					x=2;
				} else if (c == 'T') {
					x=3;
				} else if (c == 'N') {
					if (alphaLen == 5) x=4;
					else x = -1;
				}
				if (x > -1) {
					revFreq[mIndex][curIndex][x] += 1.0;
					revCounts[mIndex][curIndex] += 1.0;
				}
			}
			if (j<len-kmer) {
				int x = -1;
				char c = seq[j+kmer];
				if (c == 'A') {
					x=0;
				} else if (c == 'C') {
					x=1;
				} else if (c == 'G') {
					x=2;
				} else if (c == 'T') {
					x=3;
				} else if (c == 'N') {
					if (alphaLen == 5) x=4;
					else x = -1;
				}
				if (x > -1) {
					fwdFreq[mIndex][curIndex][x] += 1.0;
					fwdCounts[mIndex][curIndex] += 1.0;
				}

			}
		}
	}

	//normalize model, and turn frequencies into cumulative values to work better 
	// with random number generation later
	for (int i=0;i<modelLen;i++) {
		double N = counts[i];
		if (N < 1.0) N = 1.0;
		for (int j=0;j<freqLength;j++) {
			freq[i][j] /= N;
			if (j > 0) freq[i][j] += freq[i][j-1];

			double NF = fwdCounts[i][j];
			if (NF < 1.0) NF = 1.0;
			double NR = revCounts[i][j];
			if (NR < 1.0) NR = 1.0;

			for (int k=0;k<alphaLen;k++) {
				fwdFreq[i][j][k] /= NF;
				if (k > 0) fwdFreq[i][j][k] += fwdFreq[i][j][k-1];
				revFreq[i][j][k] /= NR;
				if (k > 0) revFreq[i][j][k] += revFreq[i][j][k-1];
			}
		}
	}
}
BgSequence** BgModel::generateSequences(long int numSeqs) {
	BgSequence** seqs = new BgSequence*[numSeqs];
	for (long int i=0;i<numSeqs;i++) {
		BgSequence* seq = new BgSequence();
		seq->seq = generateSequence(seqLen);
		seq->len = seqLen;
		seq->pos = -1;
		seq->strand = MOTIF_STRAND_POS;
		seq->stat = 0.0;
		seq->weight = 1.0;
		seq->bin = -1;
		seq->mask = 0;
		seq->positionIndex = -1;
		seqs[i] = seq;
	}
	return seqs;
}
char* BgModel::generateSequence(int len) {
	if (positionFlag) {
		if (len > seqLen) len = seqLen;
	}
	char* seq = new char[len+1];
	seq[len]='\0';

	int lastIndex = -1;
	for (int i=0;i<len-kmer+1;i++) {
		int mIndex = 0;
		if (positionFlag) {
			mIndex = i-1;
		}
		
		double r = ((double)rand())/((double)RAND_MAX);
//fprintf(stderr, "r=%lf\n", r);
		if (i==0) {
			//to start, pick rand kmer
			int index = -1;
			mIndex = 0;
			for (int j=0;j<freqLength;j++) {
//fprintf(stderr, "Start: %d %d %lf\n", i, j, freq[mIndex][j]);
				if (r < freq[mIndex][j]) {
					index = j;
					break;
				}
			}
			if (index == -1) index = freqLength-1;
			strcpy(seq,freqRef[index]);
			lastIndex = index;
		} else {
			int index = -1;
			for (int j=0;j<alphaLen;j++) {
//fprintf(stderr, "Start: %d %d %d %lf\n", i, j, lastIndex, fwdFreq[mIndex][lastIndex][j]);
				if (r < fwdFreq[mIndex][lastIndex][j]) {
					index = j;
					break;
				}
			}
			if (index == -1) index = alphaLen-1;
			seq[i+kmer-1] = alphaRef[index];
			lastIndex = fwdIndex[lastIndex][index];
		}
	}
	//int x = strlen(seq);
	//fprintf(stderr, "%d\t%d\t%s\n", x, len, seq);

	return seq;
}

int BgModel::seq2Index(char* s, int len) {
	int index = 0;
	int scale = 1;
	for (int i=0;i<len;i++) {
		int x = 0;
		if (s[i] == 'N') {
			if (alphaLen == 5) x=4;
			else return -1;
		} else if (s[i] == 'A') {
			x=0;
		} else if (s[i] == 'C') {
			x=1;
		} else if (s[i] == 'G') {
			x=2;
		} else if (s[i] == 'T') {
			x=3;
		} else {
			return -1;
		}
		index += x * scale;
		scale *= alphaLen;
	}
	return index;
}
// Does not set the NULL to end the string
void BgModel::index2Seq(int index, int len, char* rv) {

	//char* rv = new char[len+1];
	rv[0] = '\0';
	//rv[len] = '\0';
	int scale = alphaLen;
	for (int i=0;i<len;i++) {
		int x = index % scale;
		x /= (scale/alphaLen);
		if (x == 0) {
			rv[i] = 'A';
		} else if (x==1) {
			rv[i] = 'C';
		} else if (x==2) {
			rv[i] = 'G';
		} else if (x==3) {
			rv[i] = 'T';
		} else if (x==4) {
			rv[i] = 'N';
		} else {
			fprintf(stderr, "\tError something wrong with index2Seq - x=%d (input index=%d, len=%d)\n", x, index, len);
		}
		scale *= alphaLen;
	}
}

char* revopp2(char* seq) {
    if (seq == NULL) return NULL;
	int len = strlen(seq);
	if (len < 1) return NULL;
	char* newSeq = new char[len+1];
    for (int i=0;i<len;i++) {
        int j = len-1-i;

        if (seq[i] == 'A') {
            newSeq[j] = 'T';
        } else if (seq[i] == 'C') {
            newSeq[j] = 'G';
        } else if (seq[i] == 'G') {
            newSeq[j] = 'C';
        } else if (seq[i] == 'T') {
            newSeq[j] = 'A';
        } else if (seq[i] == 'N') {
            newSeq[j] = 'N';
		} else {
            newSeq[j] = 'N';
        }
    }
	return newSeq;
}

BgPosition::BgPosition() {
	name = NULL;
	chr = NULL;
	start = -1;
	end = -1;
	strand = -1;
	len = -1;
	chrIndex = -1;
}
BgPosition::~BgPosition() {
	if (name != NULL) {
		delete []name;
	}
	if (chr != NULL) {
		delete []chr;
	}
}
void BgPosition::set(char* n, char* c, int s, int e, char st, int l) {
	if (n != NULL) {
		name = new char[strlen(n)+1];
		strcpy(name,n);
	}
	if (c != NULL) {
		chr = new char[strlen(c)+1];
		strcpy(chr,c);
	}
	start = s;
	end = e;
	strand = st;
	len = l;
}

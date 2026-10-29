// Modified in October 2026 for homer2-ao: the bgbins, randbg and filterids commands.
// Copyright 2009 - 2016 Christopher Benner <cbenner@ucsd.edu>
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

#include "homerTools.h"
#define _FILE_OFFSET_BITS 64
#include <stdint.h>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>

void printCMD();
void printCMDbarcodes();
void printCMDtrim();
void printCMDfreq();
void printCMDextract();
void printCMDspecial();
void barcodesProgram(int argc, char** argv);
void truSeqProgram(int argc, char** argv);
void trimProgram(int argc, char** argv);
void freqProgram(int argc, char** argv);
void extractProgram(int argc, char** argv);
void clusterProgram(int argc, char** argv);
void matrix3DProgram(int argc, char** argv);
void specialProgram(int argc, char** argv);
void decontaminateProgram(int argc, char** argv);
void filterIdsProgram(int argc, char** argv);
void bgBinsProgram(int argc, char** argv);
void randBgProgram(int argc, char** argv);

int main(int argc, char** argv) {

	if (argc < 2) {
		printCMD();
	}

	char* program = argv[1];
	if (strcmp(program,"barcodes")==0) {
		barcodesProgram(argc,argv);
	} else if (strcmp(program,"trim")==0) {
		trimProgram(argc,argv);
	} else if (strcmp(program,"truseq")==0) {
		truSeqProgram(argc,argv);
	} else if (strcmp(program,"freq")==0) {
		freqProgram(argc,argv);
	} else if (strcmp(program,"extract")==0) {
		extractProgram(argc,argv);
	} else if (strcmp(program,"decontaminate")==0) {
		decontaminateProgram(argc,argv);
	} else if (strcmp(program,"cluster")==0) {
		clusterProgram(argc,argv);
	} else if (strcmp(program,"matrix3D")==0) {
		matrix3DProgram(argc,argv);
	} else if (strcmp(program,"special")==0) {
		specialProgram(argc,argv);
	} else if (strcmp(program,"filterids")==0) {
		filterIdsProgram(argc,argv);
	} else if (strcmp(program,"bgbins")==0) {
		bgBinsProgram(argc,argv);
	} else if (strcmp(program,"randbg")==0) {
		randBgProgram(argc,argv);
	} else if (strcmp(program,"--help")==0) {
		printCMD();
	} else {
		fprintf(stderr, "!!! Could not recognize \"%s\" as a command !!!\n", argv[1]);
		printCMD();
	}
	return 0;
}

void printCMD() {
	fprintf(stderr, "\n\tUsage: homerTools <command> [--help | options]\n");
	fprintf(stderr, "\n\tCollection of tools for sequence manipulation\n");
	fprintf(stderr, "\n\tCommands: [type \"homerTools <command>\" to see individual command options]\n");
	fprintf(stderr, "\t\tbarcodes - separate FASTQ file by barcodes\n");
	fprintf(stderr, "\t\ttruseq - process truseq barcodes from unidentified indexes (illumina)\n");
	fprintf(stderr, "\t\ttrim - trim adapter sequences or fixed sizes from FASTQ files(also splits)\n");
	fprintf(stderr, "\t\tfreq - calculate position-dependent nucleotide/dinucleotide frequencies\n");
	fprintf(stderr, "\t\textract - extract specific sequences from FASTA file(s)\n");
	fprintf(stderr, "\t\tdecontaminate - remove bad tags from a contaminated tag directory\n");
	fprintf(stderr, "\t\tcluster - hierarchical clustering of a NxN distance matrix\n");
//	fprintf(stderr, "\t\tmatrix3D - convert distance matrix to 3D coordinates\n");
	fprintf(stderr, "\t\tspecial - specialized routines (i.e. only really useful for chuck)\n");
	fprintf(stderr, "\t\tfilterids - keep the lines of a file whose first column is in an ID list\n");
	fprintf(stderr, "\t\tbgbins - combine target and background bins and mark targets (findMotifsGenome.pl)\n");
	fprintf(stderr, "\t\trandbg - randRemoveBackground.pl with Perl's rand() stream (findMotifsGenome.pl)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}


//------------------------------  trim -------------------------------------------------------------------


int matchAdapter5prime(char* seq, char* adapter,int maxMisMatches,int minMatchLength,int matchStart,int explicitN);
int matchAdapter3prime(char* seq, char* adapter,int maxMisMatches,int minMatchLength,int matchStart,int explicitN);
void printCMDtrim() {
	fprintf(stderr, "\n\tUsage: homerTools trim [options] <fastq file1> [file2] ...\n");
	//fprintf(stderr, "\n\tfor paired-end: homerTools trim [options] <fastq R1>,<fastq R2> ...\n");
	//fprintf(stderr, "\t\tno space between mate files and comma: <fastq R1>,<fastq R2>\n");
	fprintf(stderr, "\n\tProbably best to use only on option at a time...\n");
	fprintf(stderr, "\n\tOptions for command: trim\n");
	fprintf(stderr, "\t\t-3 <#|[ACGT]> (trim # bp or adapter sequence from 3' end of sequences)\n");
	fprintf(stderr, "\t\t-5 <#|[ACGT]> (trim # bp or adapter sequence from 5' end of sequences)\n");
	fprintf(stderr, "\t\t\t-mis <#> (Maximum allowed mismatches in adapter sequence, default: 0)\n");
	fprintf(stderr, "\t\t\t-minMatchLength <#> (minimum adapter sequence at edge to match, default: half adapter length)\n");
	fprintf(stderr, "\t\t\t-matchStart <#> (don't start searching for adapter until this position, default: 0)\n");
	fprintf(stderr, "\t\t-q <#> (Trim sequences once quality dips below threshold, default: none [range:0-40])\n");
	fprintf(stderr, "\t\t\t-qstart <#> (don't check quality until sequences are at least this long, default: 10)\n");
	fprintf(stderr, "\t\t\t-qwindow <#> (size of moving average to check for quality dropoff, default: 5)\n");
	fprintf(stderr, "\t\t-len <#> (Keep first # bp of sequence - i.e. make them the same length)\n");
	fprintf(stderr, "\t\t-stats <filename> (Output trimming statistics to filename, default: sent to stdout)\n");
	fprintf(stderr, "\t\t-min <#> (Minimum size of trimmed sequence to keep, default: 1)\n");
	fprintf(stderr, "\t\t-max <#> (Maximum read length, default: %d)\n", TRIM_MAXREADLENGTH);
	fprintf(stderr, "\t\t-suffix <filename suffix> (output is sent to InuptFileName.suffix, default: trimmed)\n");
	fprintf(stderr, "\t\t-lenSuffix <filename suffix> (length distribution is sent to InuptFileName.suffix, default: lengths)\n");
	fprintf(stderr, "\t\t-uniq (replace read names with unique IDs, for cases when there are problems with the names)\n");
	fprintf(stderr, "\t\t-pe (find look for 2nd mate pair file, matches _R1_ to _R2_ and _1.fastq to _2.fastq)\n");
	fprintf(stderr, "\t\t-explicitN (match Ns in adapter sequence directly to Ns in reads)\n");
	//fprintf(stderr, "\t\t-keepMatchOnly (when matching sequence to trim, only keep reads containing a match, default: keep all)\n");
	//fprintf(stderr, "\t\t-split <#> (Split reads into two reads at bp #, output to trimmed1 and trimmed2)\n");
	//fprintf(stderr, "\t\t-revopp <#> (Return reverse opposite of read [if used with -split, only the 2nd\n");
	//fprintf(stderr, "\t\t\t\thalf of the read will be retuned as reverse opposite])\n");
	fprintf(stderr, "\n\tCommon Examples:\n");
	fprintf(stderr, "\t\tIllumina TruSeq:\n");
	fprintf(stderr, "\t\t\thomerTools trim -3 AGATCGGAAGAGCACACGTCT -mis 2 -minMatchLength 4 -min 20\n");
	//fprintf(stderr, "\t\tIllumina TruSeq (PE - 2nd read):\n");
	//fprintf(stderr, "\t\t\thomerTools trim -3 AGATCGGAAGAGCGTCGTGT -mis 1 -minMatchLength 6\n");
	fprintf(stderr, "\t\tIllumina small RNA adapter:\n");
	fprintf(stderr, "\t\t\thomerTools trim -3 TCGTATGCCGTCTTCTGCTTG -mis 2 -minMatchLength 4 -min 20\n");
	fprintf(stderr, "\t\tOther common adapters::\n");
	fprintf(stderr, "\t\t\thomerTools trim -3 TGGAATTCTCGGGTGCCAAG -mis 2 -minMatchLength 4 -min 20\n");
	fprintf(stderr, "\t\tTypical ATAC-seq 3' adapter sequence:\n");
	fprintf(stderr, "\t\t\thomerTools trim -3 CTGTCTCTTATACACATCTC -mis 2 -minMatchLength 4 -min 20\n");
	fprintf(stderr, "\t\tHiC trimming:\n");
	fprintf(stderr, "\t\t\tMboI:    homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20\n");
	fprintf(stderr, "\t\t\tHindIII: homerTools trim -3 AAGCTAGCTT -mis 0 -matchStart 20 -min 20\n");
	fprintf(stderr, "\t\t\tNcoI:    homerTools trim -3 CCATGCATGG -mis 0 -matchStart 20 -min 20\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void trimProgram(int argc, char** argv) {

	char* suffix = NULL;
	char* lenSuffix = NULL;
	char* statsFile = NULL;
	char** files = new char*[100000];	
	int numfiles = 0;
	if (argc < 3) {
		printCMDtrim();
	}
	int p3len = -1;
	int p5len = -1;
	int p3len2 = -1;
	int p5len2 = -1;
	int matchStart = 0;
	int maxMisMatches = 0;
	int minMatchLength = -1;
	int splitPos = -1;
	int minSizeToKeep = 1;
	int qthresh = -1;
	int qwindow = 5;
	int qstart = 10;
	int uniqNamesFlag = 0;
	int findPEFlag = 0;
	int explicitNflag = 0;
	
	int revoppFlag = 0;
	int maxReadLength = TRIM_MAXREADLENGTH;
	char* p3adapter = NULL;
	char* p5adapter = NULL;
	char* p3adapter2 = NULL;
	char* p5adapter2 = NULL;
	int fixedLength = -1;
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-len")==0) {
			sscanf(argv[++i],"%d",&fixedLength);
			fprintf(stderr, "\tSequences will be output with a fixed length of %d\n",fixedLength);
		} else if (strcmp(argv[i],"-min")==0) {
			sscanf(argv[++i],"%d",&minSizeToKeep);
			fprintf(stderr, "\tMinimum length of sequence to keep %d\n",minSizeToKeep);
		} else if (strcmp(argv[i],"-keepMatchOnly")==0) {
			fprintf(stderr, "\tWill only keep sequences matching the sequence pattern\n");
		} else if (strcmp(argv[i],"-revopp")==0) {
			revoppFlag = 1;
			fprintf(stderr, "\tReturning the reverse opposite of the read\n");
		} else if (strcmp(argv[i],"-explicitN")==0) {
			explicitNflag = 1;
			fprintf(stderr, "\tTreating Ns as explicit bp\n");
		} else if (strcmp(argv[i],"-q")==0) {
			sscanf(argv[++i],"%d",&qthresh);
			fprintf(stderr, "\tWill trim reads once quality scores dip below %d\n", qthresh);
		} else if (strcmp(argv[i],"-qwindow")==0) {
			sscanf(argv[++i],"%d",&qwindow);
			fprintf(stderr, "\tWill trim reads based on quality with a %d bp moving window\n", qwindow);
		} else if (strcmp(argv[i],"-qstart")==0) {
			sscanf(argv[++i],"%d",&qstart);
			fprintf(stderr, "\tWill trim reads based on quality starting at position %d\n", qstart);
		} else if (strcmp(argv[i],"-split")==0) {
			sscanf(argv[++i],"%d",&splitPos);
			fprintf(stderr, "\tSplitting reads at bp %d\n", splitPos);
		} else if (strcmp(argv[i],"-pe")==0) {
			fprintf(stderr, "\tWill search for mate pair files that match for input files\n");
			findPEFlag =1;
		} else if (strcmp(argv[i],"-uniq")==0) {
			fprintf(stderr, "\tWill replace read names with unique names\n");
			uniqNamesFlag = 1;
		} else if (strcmp(argv[i],"-minMatchLength")==0) {
			sscanf(argv[++i],"%d",&minMatchLength);
			fprintf(stderr, "\tMinimum number of matches to adapter sequence at the edge: %d\n", minMatchLength);
		} else if (strcmp(argv[i],"-mis")==0) {
			sscanf(argv[++i],"%d",&maxMisMatches);
			fprintf(stderr, "\tMaximum number of mismatches in adapter sequence to consider %d\n", maxMisMatches);
		} else if (strcmp(argv[i],"-max")==0) {
			sscanf(argv[++i],"%d",&maxReadLength);
			fprintf(stderr, "\tMaximum length of sequence to keep %d\n", maxReadLength);
		} else if (strcmp(argv[i],"-matchStart")==0) {
			sscanf(argv[++i],"%d",&matchStart);
		} else if (strcmp(argv[i],"-q")==0) {
			sscanf(argv[++i],"%d",&qthresh);
			fprintf(stderr, "\tWill trim reads once quality scores dip below %d\n", qthresh);
		} else if (strcmp(argv[i],"-qwindow")==0) {
			sscanf(argv[++i],"%d",&qwindow);
			fprintf(stderr, "\tWill trim reads based on quality with a %d bp moving window\n", qwindow);
		} else if (strcmp(argv[i],"-qstart")==0) {
			sscanf(argv[++i],"%d",&qstart);
			fprintf(stderr, "\tWill trim reads based on quality starting at position %d\n", qstart);
		} else if (strcmp(argv[i],"-split")==0) {
			sscanf(argv[++i],"%d",&splitPos);
			fprintf(stderr, "\tSplitting reads at bp %d\n", splitPos);
		} else if (strcmp(argv[i],"-pe")==0) {
			fprintf(stderr, "\tWill search for mate pair files that match for input files\n");
			findPEFlag =1;
		} else if (strcmp(argv[i],"-uniq")==0) {
			fprintf(stderr, "\tWill replace read names with unique names\n");
			uniqNamesFlag = 1;
		} else if (strcmp(argv[i],"-minMatchLength")==0) {
			sscanf(argv[++i],"%d",&minMatchLength);
			fprintf(stderr, "\tMinimum number of matches to adapter sequence at the edge: %d\n", minMatchLength);
		} else if (strcmp(argv[i],"-mis")==0) {
			sscanf(argv[++i],"%d",&maxMisMatches);
			fprintf(stderr, "\tMaximum number of mismatches in adapter sequence to consider %d\n", maxMisMatches);
		} else if (strcmp(argv[i],"-max")==0) {
			sscanf(argv[++i],"%d",&maxReadLength);
			fprintf(stderr, "\tMaximum length of sequence to keep %d\n", maxReadLength);
		} else if (strcmp(argv[i],"-matchStart")==0) {
			sscanf(argv[++i],"%d",&matchStart);
			fprintf(stderr, "\tWill only start matching adapter sequence after %d [for 3', before for 5']\n", matchStart);
		} else if (strcmp(argv[i],"-suffix")==0) {
			suffix=argv[++i];
			fprintf(stderr, "\tEach output file will be given the original filename with .%s as a suffix\n", suffix);
		} else if (strcmp(argv[i],"-lenSuffix")==0) {
			lenSuffix=argv[++i];
			fprintf(stderr, "\tRead length distribution files will be given the original filename with .%s as a suffix\n", lenSuffix);
		} else if (strcmp(argv[i],"-stats")==0) {
			statsFile=argv[++i];
			fprintf(stderr, "\tTrimming stats will be sent to file: %s\n",statsFile);
		} else if (strcmp(argv[i],"-3")==0) {
			fprintf(stderr, "\tSequences will be trimmed at 3' end");
			i++;
			if (argv[i][0] > 47 && argv[i][0] < 58) {
				sscanf(argv[i],"%d",&p3len);	
				fprintf(stderr, " by %d bp\n", p3len);
			} else {
				fprintf(stderr, " with adapter: %s\n", argv[i]);
				p3adapter = argv[i];
				if (minMatchLength < 0) {
					minMatchLength = strlen(p3adapter)/2;
				}
			}
		} else if (strcmp(argv[i],"-3r2")==0) {
			fprintf(stderr, "\tSequences will be trimmed at 3' end of the 2nd read");
			i++;
			if (argv[i][0] > 47 && argv[i][0] < 58) {
				sscanf(argv[i],"%d",&p3len2);	
				fprintf(stderr, " by %d bp\n", p3len2);
			} else {
				fprintf(stderr, " with adapter: %s\n", argv[i]);
				p3adapter2 = argv[i];
				if (minMatchLength < 0) {
					minMatchLength = strlen(p3adapter2)/2;
				}
			}
		} else if (strcmp(argv[i],"-5")==0) {
			fprintf(stderr, "\tSequences will be trimmed at 5' end");
			i++;
			if (argv[i][0] > 47 && argv[i][0] < 58) {
				sscanf(argv[i],"%d",&p5len);	
				fprintf(stderr, " by %d bp\n", p5len);
			} else {
				fprintf(stderr, " with adapter: %s\n", argv[i]);
				p5adapter = argv[i];
				if (minMatchLength < 0) {
					minMatchLength = strlen(p5adapter)/2;
				}
			}
		} else if (strcmp(argv[i],"-5r2")==0) {
			fprintf(stderr, "\tSequences will be trimmed at 5' end of the 2nd read");
			i++;
			if (argv[i][0] > 47 && argv[i][0] < 58) {
				sscanf(argv[i],"%d",&p5len2);	
				fprintf(stderr, " by %d bp\n", p5len2);
			} else {
				fprintf(stderr, " with adapter: %s\n", argv[i]);
				p5adapter2 = argv[i];
				if (minMatchLength < 0) {
					minMatchLength = strlen(p5adapter2)/2;
				}
			}
		} else if (argv[i][0] == '-') {
			printCMDtrim();
		} else {
			fprintf(stderr, "\tWill parse file %s\n", argv[i]);
			files[numfiles++] = argv[i];
		}
	}
	if (numfiles < 1) {
		fprintf(stderr, "!!! No input files... Not much to do!!!\n");
		exit(0);
	}


	char* cmd = new char[BUFFER];
	char* file1 = new char[BUFFER];
	char* file2 = new char[BUFFER];
	char* buf = new char[BUFFER];
	char* buf2 = new char[BUFFER];
	char* current = new char[BUFFER];
	char* current2 = new char[BUFFER];
	char* curSeq = new char[BUFFER];
	char* curSeq2 = new char[BUFFER];
	char* curQual = new char[BUFFER];
	char* curQual2 = new char[BUFFER];
	char* curHeader = new char[BUFFER];
	char* curHeader2 = new char[BUFFER];
	char* newfilename = new char[10000];
	char* newfilename2 = new char[10000];
	int* readLengths = new int[maxReadLength+1];
	int* readLengths2 = new int[maxReadLength+1];
	current[0]='\0';

	for (int i=0;i<numfiles;i++) {

		//int peFlag = 0;

		file2[0]='\0';
		strcpy(file1,files[i]);
		for (int j=0;j<(int)strlen(file1);j++) {
			if (file1[j] == ',') {
				file1[j] = '\0';
				strcpy(file2,&(file1[j+1]));
				fprintf(stderr, "\tDetected PE files: %s\n", file1);
				fprintf(stderr, "\t                   %s\n", file2);
				//peFlag = 1;
				break;
			}
		}
		if (findPEFlag && file2[0] == '\0') {
			// Need to add auto detection of PE files
			strcpy(file2,file1);
			for (int j=0;j<(int)strlen(file2)-4;j++) {
				if (strncmp(&(file2[j]),"_R1.",4)==0 || 
								strncmp(&(file2[j]),"_R1_",4)==0) {
					file2[j+2]='2';
					break;
				}
				if (strncmp(&(file2[j]),"_R2.",4)==0 || 
								strncmp(&(file2[j]),"_R2_",4)==0) {
					file2[j+2]='1';
					break;
				}
				if (strncmp(&(file2[j]),"_1.f",4)==0) {
					file2[j+1]='2';
					break;
				}
				if (strncmp(&(file2[j]),"_2.f",4)==0) {
					file2[j+1]='1';
					break;
				}
			}
			FILE* check = fopen(file2, "r");
			if (check != NULL) {
				fprintf(stderr, "\tDetected PE files: %s\n", file1);
				fprintf(stderr, "\t                   %s\n", file2);
				fclose(check);
				//peFlag = 1;
			} else {
				fprintf(stderr, "\tNo PE file detected for %s\n", file1);
				file2[0]='\0';
			}
		}

		FILE* inputfp = NULL;	
		FILE* inputfp2 = NULL;	

		int gzipFlag1 = 0;
		int gzipFlag2 = 0;
		if (strlen(file1) > 3) {
			if (strcmp(&(file1[strlen(file1)-3]), ".gz")==0) {
				gzipFlag1 = 1;
			} else if (strlen(file1) > 4 && strcmp(&(file1[strlen(file1)-4]), ".bz2")==0) {
				gzipFlag1 = 2;
			}
		}
		if (strlen(file2) > 3) {
			if (strcmp(&(file2[strlen(file2)-3]), ".gz")==0) {
				gzipFlag2 = 1;
			} else if (strlen(file2) > 4 && strcmp(&(file2[strlen(file2)-4]), ".bz2")==0) {
				gzipFlag2 = 2;
			}
		}
	
		if (gzipFlag1) {
			if (gzipFlag1==1) sprintf(cmd, "zcat %s", file1);
			if (gzipFlag1==2) sprintf(cmd, "bunzip2 -c %s", file1);
			inputfp = popen(cmd,"r");
		} else {
			inputfp = fopen(file1,"r");
		}
		if (inputfp == NULL) {
			fprintf(stderr, "!!! Could not open input file: %s - skipping !!!\n", file1);
			continue;
		}
		if (file2[0] != '\0') {
			if (gzipFlag2) {
				if (gzipFlag2==1) sprintf(cmd, "zcat %s", file2);
				if (gzipFlag2==2) sprintf(cmd, "bunzip2 -c %s", file2);
				inputfp2 = popen(cmd,"r");
			} else {
				//fprintf(stderr, "cmd=%s\n", cmd);
				inputfp2 = fopen(file2,"r");
			}
			if (inputfp2 == NULL) {
				fprintf(stderr, "!!! Could not open 2nd read input file: %s - skipping !!!\n", file2);
				continue;
			}
		}

		FILE* outfp = stdout;
		FILE* outfp2 = stdout;

		strcpy(newfilename,file1);
		strcat(newfilename,".");
		strcpy(newfilename2,file2);
		strcat(newfilename2,".");

		if (suffix == NULL) {
			strcat(newfilename,"trimmed");
			strcat(newfilename2,"trimmed");
		} else {
			strcat(newfilename,suffix);
			if (splitPos > 0) {
				strcat(newfilename2,suffix);
				strcat(newfilename,"1");
				strcat(newfilename2,"2");
			} else {
				strcat(newfilename2,suffix);
			}
		}

		outfp= fopen(newfilename, "w");
		if (outfp == NULL) {
			fprintf(stderr, "!!! Could not open output file: %s - skipping !!!\n", newfilename);
			continue;
		}
		if (splitPos > 0 || inputfp2 != NULL) {
			//fprintf(stderr, "Outputing to %s\n", newfilename2);
			outfp2= fopen(newfilename2, "w");
			if (outfp2 == NULL) {
				fprintf(stderr, "!!! Could not open output file: %s - skipping !!!\n", newfilename2);
				continue;
			}
		}


		FILE* lenfp = stdout;
		FILE* lenfp2 = stdout;
		strcpy(newfilename,file1);
		strcat(newfilename,".");
		strcpy(newfilename2,file1);
		strcat(newfilename2,".");
		if (lenSuffix == NULL) {
			strcat(newfilename,"lengths");
			strcat(newfilename2,"lengths");
		} else {
			strcat(newfilename2,lenSuffix);
		}

		lenfp= fopen(newfilename, "w");
		if (lenfp == NULL) {
			fprintf(stderr, "!!! Could not open length distribution file: %s - skipping !!!\n", newfilename);
		}
		if (splitPos > 0 || inputfp2 != NULL) {
			lenfp2= fopen(newfilename2, "w");
			if (lenfp2 == NULL) {
				fprintf(stderr, "!!! Could not open length distribution file: %s - skipping !!!\n", newfilename2);
			}
		}

		long long int lastStart = 0;
		long long int currentline = 0;
		int lastLinePlus=0;
		current[0]='\0';
		current2[0]='\0';
		curSeq[0]='\0';
		curSeq2[0]='\0';
		curHeader[0]='\0';
		curHeader2[0]='\0';
		curQual[0]='\0';
		curQual2[0]='\0';
		int start = 0;
		int end = 1000000000;
		int start2 = 0;
		int end2 = 1000000000;
		int goodReads = 0;
		int totalReads = 0;
		int curLen = 0;
		int curLen2 = 0;
		int maxObservedLength = 0;
		int fastqFlag = 0;
		for (int i=0;i<maxReadLength+1;i++) readLengths[i]=0;
		for (int i=0;i<maxReadLength+1;i++) readLengths2[i]=0;

		while (fgets(buf, BUFFER, inputfp) != NULL) {

			currentline++;
			int lineLength = strlen(buf);
			int lineLength2 = 0;

			if (inputfp2 != NULL) {
				char* ok = fgets(buf2, BUFFER, inputfp2);
				if (ok == NULL) {
					fprintf(stderr, "!!! Error - Mate pair file %s has less lines than file %s\n",
								file2, file1);
					break;
				}
				lineLength2 = strlen(buf2);
			}

			if (lineLength > 1) {
				//trim off newline character
				if (lineLength > 0) buf[lineLength-1]='\0';
				if (lineLength2 > 0) buf2[lineLength2-1]='\0';

				if (((buf[0] == '@' && lastLinePlus==0) || (buf[0] == '>' && fastqFlag != 1))) {

							//&& currentline % 4 == 1) {
					if ((buf[0] == '@' && lastLinePlus==0)) fastqFlag = 1;

					int good1 = 0;
					int good2 = 1;

					if (curLen >= minSizeToKeep && curLen < maxReadLength) {
						if (maxObservedLength < curLen) maxObservedLength = curLen;
						readLengths[curLen]++;
						good1 = 1;
					} else {
						readLengths[0]++;
					}
					if (lineLength2 > 0) {
						good2 = 0;
						if (buf[0] != buf2[0]) {
							fprintf(stderr, "!!! Error, PE mate files are not in order:\n1: %s\n2: %s\n",
										buf, buf2);
							exit(0);
						}
						if (curLen2 >= minSizeToKeep && curLen2 < maxReadLength) {
							if (maxObservedLength < curLen2) maxObservedLength = curLen2;
							readLengths2[curLen2]++;
							good2 = 1;
						} else {
							readLengths2[0]++;
						}
					}
					if (good1 && good2) {
						fprintf(outfp, "%s\n%s\n+\n%s\n", curHeader,curSeq,curQual);
						if (lineLength2 > 0 || splitPos > 0) {
							fprintf(outfp2, "%s\n%s\n+\n%s\n", curHeader2,curSeq2,curQual2);
						}
						goodReads++;
					}
					totalReads++;

					current[0]='\0';
					current2[0]='\0';
					curSeq[0]='\0';
					curSeq2[0]='\0';
					curHeader[0]='\0';
					curHeader2[0]='\0';
					curQual[0]='\0';
					curQual2[0]='\0';
					strcpy(curHeader,buf);
					if (uniqNamesFlag) {
						sprintf(curHeader, "@Read%d", totalReads);
					}
					if (lineLength2 > 0) {
						strcpy(curHeader2,buf2);
						if (uniqNamesFlag) {
							sprintf(curHeader2, "@Read%d", totalReads);
						}
					} else if (splitPos > 0) {
						int x = strlen(curHeader);
						current[x-1] = '1';
						current[x] = '\0';
						strcpy(curHeader2,curHeader);
						current2[x-1] = '2';
					}
					lastStart = currentline;
					lastLinePlus=0;
					curLen = 0;
					curLen2 = 0;
				} else if (buf[0] == '+' && currentline-2==lastStart) {
					//quality header - ignore
					lastLinePlus=1;
				} else {
					if (currentline-1==lastStart) {
						//sequence
						char* seq = buf;
						curLen=lineLength-1;
						start = 0;
						end = curLen;
						if (p5len > 0) {
							if (p5len >= curLen) {
								seq[0]='\0';
								curLen=0;
								start=0;
								end=0;
							} else {
								seq = &(seq[p5len]);
								start=p5len;
								curLen-=p5len;
							}
						}
						if (p3len > 0) {
							if (p3len >= curLen) {
								seq[0]='\0';
								curLen=0;
								start=0;
								end=0;
							} else {
								seq[curLen-p3len] = '\0';
								end=start+curLen-p3len;
								curLen-=p3len;
							}
						}
						if (p5adapter != NULL && curLen > 0) {
							int offset =  matchAdapter5prime(seq, p5adapter,maxMisMatches,minMatchLength, matchStart,explicitNflag);
							if (offset >= 0) {
								seq = &(seq[offset+1]);
								curLen -= offset+1;
								start=offset+1;
								end = curLen;
							}
						}
						if (p3adapter != NULL && curLen > 0) {
							int offset =  matchAdapter3prime(seq, p3adapter,maxMisMatches,minMatchLength,matchStart,explicitNflag);
							if (offset >= 0) {
								seq[offset]='\0';
								curLen = offset;
								end = offset;
							}
						}
						if (fixedLength > 0) {
							if (curLen > fixedLength) {
								seq[fixedLength] = '\0';
								curLen = fixedLength;
								end = start + fixedLength;
							}
						}
						if (splitPos > 0) {
							if (splitPos < (int)strlen(seq)) {
								char tmp = seq[splitPos];
								seq[splitPos] = '\0';
								strcat(curSeq,seq);
								seq[splitPos] = tmp;
								if (revoppFlag) {
									int L = strlen(seq);
									seq[L-1] = '\0';
									revopp(&(seq[splitPos]));
									seq[L-1] = '\0';
								}
								strcat(curSeq2,&(seq[splitPos]));
							} else {
								strcat(curSeq,seq);
							}
						} else {
							if (revoppFlag) {
								seq[strlen(seq)-1] = '\0';
								revopp(seq);
								seq[strlen(seq)-1] = '\0';
							}
							strcat(curSeq,seq);
						}
						if (lineLength2>0) {
							//mate pair sequence
							char* seq = buf2;
							curLen2=lineLength2-1;
							start2 = 0;
							end2 = curLen2;
							if (p5len > 0 || p5len2 > 0) {
								int p5lenLocal = p5len;
								if (p5len2 > 0) p5lenLocal = p5len2;
								if (p5lenLocal >= curLen) {
									seq[0]='\0';
									curLen2=0;
									start2=0;
									end2=0;
								} else {
									seq = &(seq[p5lenLocal]);
									start2=p5lenLocal;
									curLen2-=p5lenLocal;
								}
							}
							if (p3len > 0 || p3len2 > 0) {
								int p3lenLocal = p3len;
								if (p3len2 > 0) p3lenLocal = p3len2;
								if (p3lenLocal >= curLen2) {
									seq[0]='\0';
									curLen2=0;
									start2=0;
									end2=0;
								} else {
									seq[curLen2-p3lenLocal] = '\0';
									end2=start2+curLen2-p3lenLocal;
									curLen2-=p3lenLocal;
								}
							}
							if ((p5adapter != NULL || p5adapter2 != NULL) && curLen > 0) {
								char* adapterLocal = p5adapter;
								if (p5adapter2 != NULL) {
									adapterLocal = p5adapter2;
								}
								int offset =  matchAdapter5prime(seq, adapterLocal,maxMisMatches,minMatchLength, matchStart,explicitNflag);
								if (offset >= 0) {
									seq = &(seq[offset+1]);
									curLen2 -= offset+1;
									start2=offset+1;
									end2 = curLen2;
								}
							}
							if ((p3adapter != NULL || p3adapter2 != NULL) && curLen2 > 0) {
								char* adapterLocal = p3adapter;
								if (p3adapter2 != NULL) {
									adapterLocal = p3adapter2;
								}
								int offset =  matchAdapter3prime(seq, adapterLocal,maxMisMatches,minMatchLength,matchStart,explicitNflag);
								if (offset >= 0) {
									seq[offset]='\0';
									curLen2 = offset;
									end2 = offset;
								}
							}
							if (fixedLength > 0) {
								if (curLen2 > fixedLength) {
									seq[fixedLength] = '\0';
									curLen2 = fixedLength;
									end2 = start2 + fixedLength;
								}
							}
							if (revoppFlag) {
								seq[strlen(seq)-1] = '\0';
								revopp(seq);
								seq[strlen(seq)-1] = '\0';
							}
							strcat(curSeq2,seq);
						}
					} else {
						//quality scores
						char* qual = buf;
						if (end > lineLength-1) {
							end=0;
						} else {
							qual = &(buf[start]);
						} 
						qual[end]='\0';

						//trim read based on quality score
						if (qthresh > 0) {
							int L = strlen(qual);
							int total = 0;
							int N = 0;
							for (int i=qstart;i<L;i++) {
								if (i==qstart) {
									for (int j=i;j<i+qwindow;j++) {
										if (j>=L) break;
										N++;
										total+=qual[j]-33;
									}
								} else {
									total -= qual[i-1]-33;
									if (i+qwindow < L) {
										total += qual[i+qwindow]-33;
									} else {
										N--;
									}
								}
								if (N > 0 && (total/N < qthresh)) {	
									qual[i] = '\0';
									curSeq[i] = '\0';
									curLen = i-1;
									break;
								} else if (N < 1) {
								}
							}
						}

						if (splitPos > 0) {
							if (splitPos < (int)strlen(qual)) {
								char tmp = qual[splitPos];
								qual[splitPos] = '\0';
								strcat(curQual,qual);
								qual[splitPos] = tmp;
								if (revoppFlag) {
									int L = strlen(qual);
									qual[L-1] = '\0';
									reverse(&(qual[splitPos]));
									qual[L-1] = '\0';
								}
								strcat(curQual2,&(qual[splitPos]));
							} else {
								strcat(curQual,qual);
							}
						} else {
							if (revoppFlag) {
								qual[strlen(qual)-1] = '\0';
								reverse(qual);
								qual[strlen(qual)-1] = '\0';
							}
							strcat(curQual,qual);
						}
						if (lineLength2 > 0) {
							//quality scores
							char* qual = buf2;
							if (end2 > lineLength2-1) {
								end2=0;
							} else {
								qual = &(buf2[start2]);
							} 
							qual[end2]='\0';

							//trim read based on quality score
							if (qthresh > 0) {
								int L = strlen(qual);
								int total = 0;
								int N = 0;
								for (int i=qstart;i<L;i++) {
									if (i==qstart) {
										for (int j=i;j<i+qwindow;j++) {
											if (j>=L) break;
											N++;
											total+=qual[j]-33;
										}
									} else {
										total -= qual[i-1]-33;
										if (i+qwindow < L) {
											total += qual[i+qwindow]-33;
										} else {
											N--;
										}
									}
									if (N > 0 && (total/N < qthresh)) {	
										qual[i] = '\0';
										curSeq2[i] = '\0';
										curLen2 = i-1;
										break;
									} else if (N < 1) {
									}
								}
							}
	
							if (revoppFlag) {
								qual[strlen(qual)-1] = '\0';
								reverse(qual);
								qual[strlen(qual)-1] = '\0';
							}
							strcat(curQual2,qual);
						}
					}
					lastLinePlus=0;
				}
			} else {
				lastLinePlus=0;
				if (currentline-1==lastStart || (currentline-3==lastStart && lastLinePlus==1)) {
					//blank sequence or quality information
				} else {
					//whitespace
				}
			}
		}
		//take care of last one
		if (lastLinePlus==0 || fastqFlag != 1) {
			int good1 = 0;
			int good2 = 1;
			if (curLen >= minSizeToKeep && curLen < maxReadLength) {
				if (maxObservedLength < curLen) maxObservedLength = curLen;
				readLengths[curLen]++;
				good1=1;
			} else {
				readLengths[0]++;
			}
			if (inputfp2 != NULL || splitPos > 0) {
				good2=0;
				if (curLen2 >= minSizeToKeep && curLen2 < maxReadLength) {
					if (maxObservedLength < curLen2) maxObservedLength = curLen2;
					readLengths2[curLen2]++;
					good2=1;
				} else {
					readLengths2[0]++;
				}
			}

			if (good1 && good2) {
				goodReads++;
				fprintf(outfp, "%s\n%s\n+\n%s\n", curHeader,curSeq,curQual);
				if (inputfp2 != NULL || splitPos > 0) {
					fprintf(outfp2, "%s\n%s\n+\n%s\n", curHeader2,curSeq2,curQual2);
				}
			}
		}

		fprintf(lenfp, "Length\t# reads\tFraction\n");
		for (int i=0;i<maxObservedLength+1;i++) {
			double ratio = ((double)(readLengths[i]))/((double)totalReads);
			fprintf(lenfp, "%d\t%d\t%lf%%\n", i, readLengths[i],ratio*100);
		}
		fclose(lenfp);

		if (inputfp2 != NULL || splitPos > 0) {
			fprintf(lenfp2, "Length\t# reads\tFraction\n");
			for (int i=0;i<maxObservedLength+1;i++) {
				double ratio = ((double)(readLengths2[i]))/((double)totalReads);
				fprintf(lenfp2, "%d\t%d\t%lf%%\n", i, readLengths2[i],ratio*100);
			}
			fclose(lenfp2);
		}
		fprintf(stderr, "\tTrimmed output: %d of %d reads\n", goodReads, totalReads);

		if (gzipFlag1) {
			pclose(inputfp);
		} else {
			fclose(inputfp);
		}
		if (inputfp2 != NULL) {
			if (gzipFlag2) {	
				pclose(inputfp2);
			} else {
				fclose(inputfp2);
			}
		}

		fclose(outfp);
		if (outfp2 != NULL) fclose(outfp2);

	}		

	delete []current;
	delete []newfilename;
	delete []newfilename2;
	delete []readLengths;
	delete []files;
	delete []buf;
	
	exit(0);
}

int matchAdapter5prime(char* seq, char* adapter, int maxMisMatches,int minMatchLength,int matchStart, int explicitNflag) {
	int seqLength = strlen(seq)-1; //seq is expected to have a \n character
	int adapterLength = strlen(adapter);
	for (int i=seqLength-1;i>=matchStart;i--) {
		int match=1;
		int numMis = 0;
		int numMatch = 0;
		for (int j=adapterLength-1;j>=0;j--) {
			int nindex = i-(adapterLength-1)+j;
			if (nindex<0) break;
			numMatch++;
			if (explicitNflag==0) {
				if (seq[nindex] == 'N') continue;
				if (adapter[j] == 'N') continue;
			}
			if (seq[nindex] == adapter[j]) {
				continue;
			} else {
				numMis++;
				if (numMis > maxMisMatches) {
					match=0;
					break;
				}
			}
		}
		if (numMatch < minMatchLength) match = 0;
		if (match==1) {
			return i;
		}
	}
	return -1;
}

int matchAdapter3prime(char* seq, char* adapter, int maxMisMatches, int minMatchLength, int matchStart, int explicitNflag) {
	int seqLength = strlen(seq)-1; //seq is expected to have a \n character
	int adapterLength = strlen(adapter);
	for (int i=matchStart;i<seqLength;i++) {
		int match=1;
		int numMis = 0;
		int numMatch = 0;
		for (int j=0;j<adapterLength;j++) {
			if (i+j>=seqLength) break;
			numMatch++;
			if (explicitNflag==0) {
				if (seq[i+j] == 'N') continue;
				if (adapter[j] == 'N') continue;
			}
			if (seq[i+j] == adapter[j]) {
				continue;
			} else {
				numMis++;
				if (numMis > maxMisMatches) {
					match=0;
					break;
				}
			}
		}
		if (numMatch < minMatchLength) match = 0;
		if (match==1) {
			return i;
		}
	}
	return -1;
}
//---------------------------- truseq -------------------------------------------------------------------

#define MAX_FREQ_TO_REPORT 30
void getTruSeqBarcode(char* buf,int lineLength,char* barcode);
int processBarcode(char* barcode, char* read, Hashtable* barcodeFPs, LongInttable* barcodeTotals,
							char* prefix,int printFlag);
void printCMDtruSeq() {
	fprintf(stderr, "\n\tUsage: homerTools truseq <fastq file1> [file2] ...\n");
	fprintf(stderr, "\n\tOptions for command: truseq\n");
	fprintf(stderr, "\t\t-min <#> (Minimum frequency of barcodes to keep: default=%.3f\n", BARCODE_MINFREQ);
	fprintf(stderr, "\t\t-freq <filename> (output file for barcode frequencies, default=file.freq.txt)\n");
	fprintf(stderr, "\t\t-qual <#> (Minimum quality score for barcode nucleotides, default=not used)\n");
	fprintf(stderr, "\t\t-qualBase <character> (Minimum quality character in FASTQ file, default=B)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void truSeqProgram(int argc, char** argv) {

	char* freqFileName = NULL;
	char** files = new char*[100000];	
	int numfiles = 0;
	int qualThreshold = -120;
	char qualBase = 'B';
	int printFlag = 0;

	if (argc < 3) printCMDtruSeq();

	//int barcodeSize = 0;
	double minPercentage = BARCODE_MINFREQ;
	fprintf(stderr, "\tDefault minimum barcode percentage is %.0lf%%\n", minPercentage*100);

	for (int i=2;i<argc;i++) { 
		if (strcmp(argv[i],"-min")==0) {
			sscanf(argv[++i],"%lf",&minPercentage);
			fprintf(stderr, "\tMinimum percentage for barcodes set at %.2lf%%\n", minPercentage*100);
		} else if (strcmp(argv[i],"-freq")==0) {
			freqFileName = argv[++i];
			fprintf(stderr, "\tBarcode frequencies will be recorded in *.%s files\n", freqFileName);
		} else if (strcmp(argv[i],"-qualBase")==0) {
			sscanf(argv[++i],"%c",&qualBase);
		} else {
			files[numfiles++] = argv[i];
		}
	}

	if (numfiles < 1) {
		fprintf(stderr, "!!! No input files... Not much to do!!!\n");
		exit(0);
	}


	char* prefix = NULL;
	char* buf = new char[BUFFER];
	char* current = new char[BUFFER];
	char* barcode = new char[BUFFER];
	char* name= new char[10000];
	current[0]='\0';
	barcode[0]='\0';

	for (int i=0;i<numfiles;i++) {
		FILE* fp = fopen(files[i], "r");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s\n", files[i]);
			continue;
		}
		prefix = files[i];
		Hashtable* barcodeFPs = new Hashtable();	
		LongInttable* barcodeTotals = new LongInttable();	
		long long int totalReads = 0;
		long long int lastStart = 0;
		long long int currentline = 0;
		int lastLinePlus=0;
		int format = SEQFILE_FORMAT_UNKNOWN ;
		char minQualityScore=0;
		current[0]='\0';
		barcode[0]='\0';
		while (fgets(buf, BUFFER, fp) != NULL) {
			currentline++;
			int lineLength = strlen(buf);
			if (lineLength > 1) {
				if (format == SEQFILE_FORMAT_UNKNOWN) {
					if (buf[0] == '@') {
						format = SEQFILE_FORMAT_FASTQ;
						fprintf(stderr, "\tFASTQ format detected\n");
					} else if (buf[0] == '>') {
						format = SEQFILE_FORMAT_FASTA;
						fprintf(stderr, "\tFASTA format detected\n");
					}
				}
				if ((format == SEQFILE_FORMAT_FASTQ && buf[0] == '@' && lastLinePlus==0)  
						|| (format == SEQFILE_FORMAT_FASTA && buf[0] == '>')) {
					if (minQualityScore >= qualThreshold) {
						totalReads += processBarcode(barcode,current,barcodeFPs,barcodeTotals,prefix,printFlag);
					}
					current[0]='\0';
					barcode[0]='\0';
					strcpy(current,buf);
					lastStart = currentline;
					lastLinePlus=0;
					minQualityScore=127;
					getTruSeqBarcode(buf,lineLength,barcode);
					//fprintf(stdout, "%s\n", barcode);

				} else if (buf[0] == '+' && currentline-2==lastStart) {
					//quality header
					lastLinePlus=1;
					strcat(current,buf);
				} else {
					lastLinePlus=0;
					strcpy(current,buf);
				}
			} else {
				lastLinePlus=0;
				if (currentline-1==lastStart || (currentline-3==lastStart && lastLinePlus==1)) {
					//blank sequence or quality information
					strcat(current,"\n");
				} else {
					//whitespace
				}
			}
		}
		if ((format == SEQFILE_FORMAT_FASTQ && lastLinePlus==0)  
								|| (format == SEQFILE_FORMAT_FASTA)) {
			if (minQualityScore >= qualThreshold) {
				totalReads += processBarcode(barcode,current,barcodeFPs,barcodeTotals,prefix,printFlag);
			}
		}
		fclose(fp);

		char** keys = barcodeTotals->keys();
		StrSort* data = new StrSort[barcodeTotals->total];
		for (int j=0;j<barcodeTotals->total;j++) {
			long long int barcodeTotal = barcodeTotals->search(keys[j]);
			if (printFlag) {
				FILE* bfp = (FILE*) barcodeFPs->search(keys[j]);
				fclose(bfp);
			}
			double ratio = ((double)barcodeTotal)/((double)totalReads);
			char* readInfo = new char[10000];
			sprintf(readInfo, "%s\t%lld\t%lf\n", keys[j], barcodeTotal,ratio);
			data[j].str = readInfo;
			data[j].v = ratio;
			if (ratio < minPercentage) {
				strcpy(name, prefix);
				strcat(name,".");
				strcat(name,keys[j]);
				remove(name);
			}
			delete [](keys[j]);
		}
		delete []keys;
	
		qsort(data,barcodeTotals->total,sizeof(StrSort),&decendStrSort);

		FILE* freqFP = NULL;
		strcpy(name,files[i]);
		strcat(name,".freq.txt");
		freqFP = fopen(name, "w");
		if (freqFP == NULL) {
			fprintf(stderr, "!!! Problem opening %s for writing (frequency file)!!!\n", name);
			exit(0);
		}
		fprintf(freqFP, "Barcode\tTotal Reads(of %lld)\tFrequency\n", totalReads);
		fprintf(stderr, "\tBarcode\tTotal Reads(of %lld)\tFrequency\n", totalReads);
		for (int i=0;i<barcodeTotals->total;i++) {
			fprintf(freqFP,"%s",data[i].str);
			if (i<MAX_FREQ_TO_REPORT) fprintf(stderr,"\t\t%s",data[i].str);
			delete [](data[i].str);
		}
		delete []data;
	
		if (freqFileName != NULL) {
			fclose(freqFP);
		}
		delete barcodeFPs;
		delete barcodeTotals;
	}

	delete []name;
	delete []files;
	delete []buf;
}
void getTruSeqBarcode(char* buf,int lineLength,char* barcode) {
	int lastColon = 0;
	int lastBp = 0;
	for (int i=0;i<lineLength;i++) {
		if (buf[i] == ':') {
			lastColon=i;
		} else if (buf[i] == 'A' || buf[i] == 'C' || buf[i]=='G' || buf[i]=='T' || buf[i] == 'N') {
			lastBp = i;
		}
	}
	if (lastBp-lastColon < 1) {
		barcode[0] = '\0';
	} else {
		strncpy(barcode,&(buf[lastColon+1]),lastBp-lastColon);
	}
}

//------------------------------  barcodes -------------------------------------------------------------------


void printCMDbarcodes() {
	fprintf(stderr, "\n\tUsage: homerTools barcodes <# bp in barcode> [options] <fastq file1> [file2] ...\n");
	fprintf(stderr, "\n\t3rd argument must be the number bp in the barcode\n");
	fprintf(stderr, "\n\tOptions for command: barcode\n");
	fprintf(stderr, "\t\t-min <#> (Minimum frequency of barcodes to keep: default=%.3f\n", BARCODE_MINFREQ);
	fprintf(stderr, "\t\t-freq <filename> (output file for barcode frequencies, default=file.freq.txt)\n");
	fprintf(stderr, "\t\t-qual <#> (Minimum quality score for barcode nucleotides, default=not used)\n");
	fprintf(stderr, "\t\t-qualBase <character> (Minimum quality character in FASTQ file, default=B)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void barcodesProgram(int argc, char** argv) {

	char* freqFileName = NULL;
	char** files = new char*[100000];	
	int numfiles = 0;
	int qualThreshold = -120;
	char qualBase = 'B';
	int printFlag = 1;

	if (argc < 4) printCMDbarcodes();

	int barcodeSize = 0;
	double minPercentage = BARCODE_MINFREQ;
	sscanf(argv[2],"%d",&barcodeSize);
	if (barcodeSize < 1 || barcodeSize > 100) {
		fprintf(stderr, "!!! Barcode of size %d seems odd - make sure you use the program correctly!!!\n", barcodeSize);
	} else {
		fprintf(stderr, "\tUsing barcode size of %d bp\n", barcodeSize);
	}
	fprintf(stderr, "\tDefault minimum barcode percentage is %.0lf%%\n", minPercentage*100);

	for (int i=3;i<argc;i++) { 
		if (strcmp(argv[i],"-min")==0) {
			sscanf(argv[++i],"%lf",&minPercentage);
			fprintf(stderr, "\tMinimum percentage for barcodes set at %.2lf%%\n", minPercentage*100);
		} else if (strcmp(argv[i],"-freq")==0) {
			freqFileName = argv[++i];
			fprintf(stderr, "\tBarcode frequencies will be recorded in *.%s files\n", freqFileName);
		} else if (strcmp(argv[i],"-qualBase")==0) {
			sscanf(argv[++i],"%c",&qualBase);
			fprintf(stderr, "\tQuality base character set to %c (%d)\n", qualBase, qualBase);
		} else if (strcmp(argv[i],"-qual")==0) {
			sscanf(argv[++i],"%d",&qualThreshold);
			fprintf(stderr, "\tBarcode quality threshold set at %d\n", qualThreshold);
		} else if (argv[i][0] == '-') {
			printCMDbarcodes();
		} else {
			fprintf(stderr, "\tWill parse file %s\n", argv[i]);
			files[numfiles++] = argv[i];
		}
	}
	if (numfiles < 1) {
		fprintf(stderr, "!!! No input files... Not much to do!!!\n");
		exit(0);
	}


	char* prefix = NULL;
	char* buf = new char[BUFFER];
	char* current = new char[BUFFER];
	char* barcode = new char[BUFFER];
	char* name= new char[10000];
	current[0]='\0';
	barcode[0]='\0';

	for (int i=0;i<numfiles;i++) {
		FILE* fp = fopen(files[i], "r");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s\n", files[i]);
			continue;
		}
		prefix = files[i];
		Hashtable* barcodeFPs = new Hashtable();	
		LongInttable* barcodeTotals = new LongInttable();	
		long long int totalReads = 0;
		long long int lastStart = 0;
		long long int currentline = 0;
		int lastLinePlus=0;
		int format = SEQFILE_FORMAT_UNKNOWN ;
		char minQualityScore=0;
		current[0]='\0';
		barcode[0]='\0';
		while (fgets(buf, BUFFER, fp) != NULL) {
			currentline++;
			int lineLength = strlen(buf);
			if (lineLength > 1) {
				if (format == SEQFILE_FORMAT_UNKNOWN) {
					if (buf[0] == '@') {
						format = SEQFILE_FORMAT_FASTQ;
						fprintf(stderr, "\tFASTQ format detected\n");
					} else if (buf[0] == '>') {
						format = SEQFILE_FORMAT_FASTA;
						fprintf(stderr, "\tFASTA format detected\n");
					}
				}
				if ((format == SEQFILE_FORMAT_FASTQ && buf[0] == '@' && lastLinePlus==0)  
						|| (format == SEQFILE_FORMAT_FASTA && buf[0] == '>')) {
					if (minQualityScore >= qualThreshold) {
						totalReads += processBarcode(barcode,current,barcodeFPs,barcodeTotals,prefix,printFlag);
					}
					current[0]='\0';
					barcode[0]='\0';
					strcpy(current,buf);
					lastStart = currentline;
					lastLinePlus=0;
					minQualityScore=127;
				} else if (buf[0] == '+' && currentline-2==lastStart) {
					//quality header
					lastLinePlus=1;
					strcat(current,buf);
				} else {
					if (lineLength-1 >= barcodeSize) {
						strcat(current,&(buf[barcodeSize]));
					} else {
						strcat(current,"\n");
					}
					if (currentline-1==lastStart) {
						//sequence - get barcode
						if (lineLength-1 >= barcodeSize) {
							strncpy(barcode,buf,barcodeSize);
						}
					} else {
						if (lineLength-1 < barcodeSize) {
							minQualityScore = -1;
						} else {
							minQualityScore = buf[0]-qualBase;
							for (int j=1;j<barcodeSize;j++) {
								if (buf[j]-qualBase < minQualityScore) minQualityScore = buf[j]-qualBase;
							}
						}
					}
					lastLinePlus=0;
				}
			} else {
				lastLinePlus=0;
				if (currentline-1==lastStart || (currentline-3==lastStart && lastLinePlus==1)) {
					//blank sequence or quality information
					strcat(current,"\n");
				} else {
					//whitespace
				}
			}
		}
		if ((format == SEQFILE_FORMAT_FASTQ && lastLinePlus==0)  
								|| (format == SEQFILE_FORMAT_FASTA)) {
			if (minQualityScore >= qualThreshold) {
				totalReads += processBarcode(barcode,current,barcodeFPs,barcodeTotals,prefix,printFlag);
			}
		}
		fclose(fp);

		char** keys = barcodeFPs->keys();
		StrSort* data = new StrSort[barcodeFPs->total];
		for (int j=0;j<barcodeFPs->total;j++) {
			long long int barcodeTotal = barcodeTotals->search(keys[j]);
			FILE* bfp = (FILE*) barcodeFPs->search(keys[j]);
			fclose(bfp);
			double ratio = ((double)barcodeTotal)/((double)totalReads);
			char* readInfo = new char[10000];
			sprintf(readInfo, "%s\t%lld\t%lf\n", keys[j], barcodeTotal,ratio);
			data[j].str = readInfo;
			data[j].v = ratio;
			if (ratio < minPercentage) {
				strcpy(name, prefix);
				strcat(name,".");
				strcat(name,keys[j]);
				remove(name);
			}
			delete [](keys[j]);
		}
		delete []keys;
	
		qsort(data,barcodeFPs->total,sizeof(StrSort),&decendStrSort);

		FILE* freqFP = NULL;
		strcpy(name,files[i]);
		strcat(name,".freq.txt");
		freqFP = fopen(name, "w");
		if (freqFP == NULL) {
			fprintf(stderr, "!!! Problem opening %s for writing (frequency file)!!!\n", name);
			exit(0);
		}
		fprintf(freqFP, "Barcode\tTotal Reads(of %lld)\tFrequency\n", totalReads);
		fprintf(stderr, "\tBarcode\tTotal Reads(of %lld)\tFrequency\n", totalReads);
		for (int i=0;i<barcodeFPs->total;i++) {
			fprintf(freqFP,"%s",data[i].str);
			fprintf(stderr,"\t\t%s",data[i].str);
			delete [](data[i].str);
		}
		delete []data;
	
		if (freqFileName != NULL) {
			fclose(freqFP);
		}
		delete barcodeFPs;
		delete barcodeTotals;
	}

	delete []name;
	delete []files;
	delete []buf;
}


int processBarcode(char* barcode, char* read, Hashtable* barcodeFPs, LongInttable* barcodeTotals,
												char* prefix,int printFlag) {
	if (barcode[0] == '\0') return 0;
	if (read[0] == '\0') return 0;
	long long int total = barcodeTotals->search(barcode);
	FILE* fp = NULL;
	if (total == EMPTY_INT) {
		barcodeTotals->insert(1,barcode);
	} else {
		barcodeTotals->insert(total+1,barcode);
	}	
	if (printFlag) {
		if (total == EMPTY_INT) {
			char* name= new char[10000];
			strcpy(name, prefix);
			strcat(name, ".");
			strcat(name, barcode);
			fp = fopen(name,"w");
			barcodeFPs->insert(fp,barcode);
			delete []name;
		} else {
			barcodeTotals->insert(total+1,barcode);
			fp = (FILE*) barcodeFPs->search(barcode);
		}
		fprintf(fp,"%s",read);
	}
	return 1;
}



int decendStrSort(const void* a, const void* b) {
	StrSort* aa = (StrSort*)a;
	StrSort* bb = (StrSort*)b;
	if (aa->v > bb->v) return -1;
	if (aa->v < bb->v) return 1;
	return 0;
}


//---------------------------- freq ---------------------------


void printCMDfreq() {
	fprintf(stderr, "\n\tUsage: homerTools freq [options] <sequence file>\n");
	fprintf(stderr, "\n\tOutputs nucleotide frequencies to stdout\n");
	fprintf(stderr, "\n\tOptions for command: freq\n");
	fprintf(stderr, "\t\t-format <tsv|fasta|fastq> (sequence file format, default: auto detect)\n");
	fprintf(stderr, "\t\t-offset <#> (offset of first base in output file, default: 0)\n");
	fprintf(stderr, "\t\t-maxlen <#> (Maximum length of sequences to consider, default: length of 1st seq)\n");
	fprintf(stderr, "\t\t-o <filename> (Output filename, default: output sent to stdout)\n");
	fprintf(stderr, "\t\t-gc <filename> (calculate CpG/GC content per sequence output to \"filename\")\n");
	fprintf(stderr, "\t\t\tOutputFormat: name<tab>CpG<tab>GC<tab>AG<tab>AC<tab>Length\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void freqProgram(int argc, char** argv) {

	char* outputFilename = NULL;
	char* gcFilename = NULL;
	char* seqFile = NULL;
	int fileFormat = SEQFILE_FORMAT_UNKNOWN;
	int maxLength = -1;
	int offset = 0;

	if (argc < 3) {
		printCMDfreq();
	}
	fprintf(stderr, "\n");
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-offset")==0) {
			sscanf(argv[++i],"%d",&offset);
			//fprintf(stderr, "\tOffset set to %d\n",offset);
		} else if (strcmp(argv[i],"-maxlen")==0) {
			sscanf(argv[++i],"%d",&maxLength);
			//fprintf(stderr, "\tOffset set to %d\n",offset);
		} else if (strcmp(argv[i],"-format")==0) {
			i++;
			if (strcmp(argv[i],"tsv")==0) {
				fileFormat = SEQFILE_FORMAT_TSV;
			} else if (strcmp(argv[i],"fasta")==0) {
				fileFormat = SEQFILE_FORMAT_FASTA;
			} else if (strcmp(argv[i],"fastq")==0) {
				fileFormat = SEQFILE_FORMAT_FASTQ;
			} else if (strcmp(argv[i],"fa")==0) {
				fileFormat = SEQFILE_FORMAT_FASTA;
			} else if (strcmp(argv[i],"fq")==0) {
				fileFormat = SEQFILE_FORMAT_FASTQ;
			} else {
				fprintf(stderr, "!!! Could not recognize format: %s !!!\n", argv[i]);
				printCMDfreq();
			}
		} else if (strcmp(argv[i],"-o")==0) {
			outputFilename = argv[++i];
		} else if (strcmp(argv[i],"-gc")==0) {
			gcFilename = argv[++i];
		} else {
			if (seqFile == NULL) {
				seqFile = argv[i];
			} else {
				fprintf(stderr, "!!! Error- multiple sequence files specified: %s & %s\n", seqFile, argv[i]);
				printCMDfreq();
			}
		}
	}

	int autoFormat = determineSeqFileFormat(seqFile);
	if (fileFormat != SEQFILE_FORMAT_UNKNOWN) {
		if (autoFormat != fileFormat) {
			fprintf(stderr, "!!! Auto-check for your file format disagrees with your selected parameters\n");
			fprintf(stderr, "!!! Check your file and your options! Continuing anyway...\n");
		}
	} else {
		fileFormat = autoFormat;
	}

	NucleotideFreq* nf = new NucleotideFreq(seqFile, fileFormat, offset, maxLength,gcFilename);

	FILE* fp = stdout;
	if (outputFilename != NULL) {
		fp = fopen(outputFilename, "w");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for writing!!!\n", outputFilename);
			exit(0);
		}
	}

	nf->print(fp);

	if (outputFilename != NULL) {
		fclose(fp);
	}
	fprintf(stderr, "\n");
}

//---------------------------- extract ---------------------------


void printCMDextract() {
	fprintf(stderr, "\n\tUsage: homerTools extract <peak file/BED file> <Directory of FASTA files> [options]\n");
	fprintf(stderr, "\n\tThe <Directory of FASTA files> can be a single FASTA file instead.\n");
	fprintf(stderr, "\tIf using a Directory, files should be named with chromosomes,\n");
	fprintf(stderr, "\t\ti.e. chr1.fa or chr1.fa.masked or genome.fa/genome.fa.masked\n");
	fprintf(stderr, "\t\tIf having trouble, place all FASTA entries in single file instead of a directory\n");
	fprintf(stderr, "\t\tFASTA format: >chrname ... (anything after whitespace will be ignored)\n");
	fprintf(stderr, "\tThis program will output sequences to stdout in tab-delimited format\n");
	fprintf(stderr, "\n\tOptions for command: extract\n");
	fprintf(stderr, "\t\t-fa (output sequences in FASTA format - default is tab-delimited format)\n"); 
	fprintf(stderr, "\t\t-mask (mask out lower case sequence from genome)\n"); 
	fprintf(stderr, "\n\tAlternate Usage: homerTools extract stats <Directory of FASTA files>\n");
	fprintf(stderr, "\t\tDisplays stats about the genome files (such as length)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void extractProgram(int argc, char** argv) {

	if (argc < 4) {
		printCMDextract();
	}

	char* posfile = argv[2];
	char* genomeDir = argv[3];
	int maskFlag = 0;
	struct stat *buf = NULL;

	if (stat(genomeDir,buf) != 0) {
		HomerConfig*  hc = new HomerConfig();
		char* genomeDirectory = hc->getGenomeDirectory(genomeDir);
		if (genomeDirectory == NULL) {
			fprintf(stderr, "!!!Error: %s not a valid file, directory or homer configured genome\n", genomeDir);
			exit(0);
		}
		genomeDir = genomeDirectory;
	}

	int fastaFlag = 0;
	int statsFlag = 0;
	if (strcmp(posfile,"stats")==0) {
		statsFlag = 1;
	}

	fprintf(stderr, "\n");
	for (int i=4;i<argc;i++) {
		if (strcmp(argv[i],"-fa")==0) {
			fastaFlag = 1;
		} else if (strcmp(argv[i],"-mask")==0) {
			maskFlag = 1;
		} else {
			printCMDextract();
		}
	}

	FILE* fpout = stdout;

	if (statsFlag) {
		PeakLibrary::extractSequenceStats(genomeDir,fpout);
	} else {
		PeakLibrary* peaks = new PeakLibrary(posfile, PEAK_READ_MODE_NORMAL);
		peaks->extractSequence(genomeDir,fpout,fastaFlag,maskFlag);
	}

	exit(0);
}



//---------------------------- decontaminate ---------------------------


void printCMDdecontaminate() {
	fprintf(stderr, "\n\tUsage: homerTools decontaminate <Contaminated Tag Directory>\n");
	fprintf(stderr, "\t\t\t\t<Contaminant(background) Tag Directory> [options]\n");
	fprintf(stderr, "\n\tRemoves reads from an experiment that contains reads originating from another.\n");
	fprintf(stderr, "\tEither specify the fraction of contaminated reads (-frac <#>), or the program\n");
	fprintf(stderr, "\twill attempt to estimate (only works for quanitatively different experiments)\n");
	fprintf(stderr, "\tCreates contaminationHistogram.txt and contaminationScatterPlot.txt in the\n");
	fprintf(stderr, "\toutput directory to help guage the extent of contamination when estimating.\n");
	fprintf(stderr, "\n\tTwo tag directories are required after \"homerTools decontaminate\":\n");
	fprintf(stderr, "\t\t<Contaminated Tag Directory> <Contaminant Tag Directory>\n");
	fprintf(stderr, "\t\tThe Contaminated Tag Directory will be modified (recommend to copy the original)\n");
	fprintf(stderr, "\n\tOptions for command: decontaminate\n");
	fprintf(stderr, "\t\t-frac <#> (Estimate fraction of sample that is contaminated, default: auto)\n"); 
	fprintf(stderr, "\t\t-estimateOnly (Only estimate the contamination, do not decontaminate)\n"); 
	fprintf(stderr, "\t\t-o <output tag directory> (default: overrites contaminated tag directory)\n"); 
	fprintf(stderr, "\t\t-size <#> (Peak size for estimating contamination/Max distance from contaminant\n");
	fprintf(stderr, "\t\t\treads to remove contaminated reads, default: 250)\n"); 
	fprintf(stderr, "\t\t-min <#> (Minimum tag count to consider when estimating contamination, default: 20)\n"); 
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void decontaminateProgram(int argc, char** argv) {

	if (argc < 4) {
		printCMDdecontaminate();
	}

	char* expDir = argv[2];
	char* inputDir = argv[3];
	char* outputDir = NULL;
	double fraction = -1.0;
	float minThreshold = 20.0;
	int estonly = 0;

	int size = 250;

	fprintf(stderr, "\n");
	for (int i=4;i<argc;i++) {
		if (strcmp(argv[i],"-frac")==0) {
			sscanf(argv[++i],"%lf",&fraction);
		} else if (strcmp(argv[i],"-o")==0) {
			outputDir= argv[++i];
		} else if (strcmp(argv[i],"-estimateOnly")==0) {
			estonly = 1;
		} else if (strcmp(argv[i],"-min")==0) {
			sscanf(argv[++i],"%f",&minThreshold);
		} else if (strcmp(argv[i],"-size")==0) {
			sscanf(argv[++i],"%d",&size);
		} else {
			printCMDdecontaminate();
		}
	}
	fprintf(stderr, "\n\tContaminated Directory: %s\n", expDir);
	fprintf(stderr, "\tContaminant Directory:  %s\n", inputDir);
	fprintf(stderr, "\n");
	if (outputDir == NULL) {
		if (estonly == 0) {
			fprintf(stderr, "\tWarning - this will overrite the data in %s\n", expDir);
			fprintf(stderr, "\tIt is recommended that you make a copy first\n");
			fprintf(stderr, "\tIf you wish to stop, hit CTRL+C NOW!!!\n");
			for (int i=5;i>0;i--) {
				fprintf(stderr, "\t\t%d\n",i);
				(void)system("sleep 1");
			}
		}
		outputDir = expDir;
	} else {
		fprintf(stderr, "\tOutput Directory:    %s\n", outputDir);
		char* command = new char[10000];
		sprintf(command, "mkdir -p \"%s\"",outputDir);
		(void)system(command);
		sprintf(command, "cp \"%s\"/* \"%s\"/",expDir, outputDir);
		(void)system(command);
	}

	TagLibrary* tags = new TagLibrary(outputDir);
	TagLibrary* input = new TagLibrary(inputDir);

	tags->readTagDirectory();
	input->readTagDirectory();
	if (fraction < 0.0) {
		fprintf(stderr, "\tEstimating fraction of reads in from contaminated experiement\n\n");
		fraction = tags->estimateContamination(input,size,minThreshold);
		if (estonly) exit(0);
	} else {
		fprintf(stderr, "\tFraction of reads in from contaminated experiement: %lf%% [user defined]\n",
											fraction*100.0);
	}
	tags->decontaminate(input,fraction,size);
	//should probably do more quality control...

	fprintf(stderr, "\n");
	exit(0);
}


//---------------------------- special ---------------------------

void specialUniqMapProgram(int argc, char** argv);

void printCMDspecial() {
	fprintf(stderr, "\n\tUsage: homerTools special <2nd command> ...\n");
	fprintf(stderr, "\n\tHighly specialized programs\n");
	fprintf(stderr, "\n\t2nd Commands:\n");
	fprintf(stderr, "\t\tuniqmap (parse Unique mappability files)\n"); 
	//fprintf(stderr, "\t\ttile (create peak file tiling the genome)\n"); 
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void specialProgram(int argc, char** argv) {

	if (argc < 3) {
		printCMDspecial();
	}
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"uniqmap")==0) {
			specialUniqMapProgram(argc,argv);
		} else {
			printCMDspecial();
		}
	}
	exit(0);
}

//---------------------------- special uniqmap ---------------------------


void printCMDspecialUniqMap() {
	fprintf(stderr, "\n\tUsage: homerTools special uniqmap <uniqmap directory(output)> <uniq map file(input)>\n");
	fprintf(stderr, "\n\tOptions for command: special uniqmap\n");
	fprintf(stderr, "\t\tnone\n"); 
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void specialUniqMapProgram(int argc, char** argv) {

	if (argc < 5) {
		printCMDspecialUniqMap();
	}

	char* directory = argv[3];
	char* inputFile = argv[4];

	char* command = new char[100000];
	strcpy(command, "mkdir -p \"");
	strcat(command, directory);
	strcat(command, "\"");
	(void)system(command);

	char* buf = new char[BUFFER];
	char** line = new char*[10000];
	int numCols = 0;

	int initializationSize = 400000000;

	FILE* fp = fopen(inputFile,"r");
	if (fp == NULL) {
		fprintf(stderr, "!! Could not open input file: %s\n", inputFile);
		exit(0);
	}

	Hashtable* uniqmapchr = new Hashtable();
	char* currentChr = new char[10000];
	currentChr[0] = '\0';


	UniqMapChrs* umc = NULL;
	int format = -1;

	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf, line, numCols, '\t');
		
		if (format == -1) {
			if (numCols < 3) continue;
			char* chr = line[0];
			int count = 0;
			while (chr[count] != '\0') {
				if (chr[count] == '.') {
					chr[count] = '\0';
					break;
				}
				count++;
			}
			if (strcmp(chr,currentChr) != 0) {
				strcpy(currentChr,chr);
				umc = (UniqMapChrs*) uniqmapchr->search(chr);
				if (umc == NULL) {
					umc = new UniqMapChrs(chr,directory,1); //1 to overwrite existing data
					umc->increaseMapSize(initializationSize);
					uniqmapchr->insert(umc,chr);
				}
			}
			int position = -1;
			int strand = -1;
			sscanf(line[1],"%d",&position);
			sscanf(line[2],"%d",&strand);
			umc->setMappable(position,strand);
		}
	}
	fclose(fp);

	sprintf(command, "%s/uniqMapStats.txt", directory);
	fp = fopen(command, "w");
	long long int gsize=0;
	long long int mappability = 0;
	fprintf(fp, "#Chr\tMappable bp (both strands)\tTotal bp (both strand)\n");

	char** keys = uniqmapchr->keys();
	for (int i=0;i<uniqmapchr->total;i++) {
		umc = (UniqMapChrs*) uniqmapchr->search(keys[i]);
		umc->saveFiles();
		gsize += umc->size;
		mappability += umc->numMappable;
		fprintf(fp, "%s\t%d\t%d\n",keys[i],umc->numMappable,umc->size);
		delete [](keys[i]);
	}
	delete []keys;
	fprintf(fp, "genome\t%lld\t%lld\n",mappability,gsize);
	fclose(fp);
	

	delete uniqmapchr;
	delete []buf;
	delete []currentChr;
	delete []line;
	delete []command;
	exit(0);
}


//---------------------------- special cluster ---------------------------


void printCMDcluster() {
	fprintf(stderr, "\n\tUsage: homerTools cluster [options]\n");
	fprintf(stderr, "\n\tOptions for command: cluster\n");
	fprintf(stderr, "\t\t-d <distance matrix file> (tab delimited)\n"); 
	fprintf(stderr, "\t\t\t-annCols <#> (number of annotation columns at beginning, def: 1)\n"); 
	fprintf(stderr, "\t\t-i <raw data matrix file> (tab delimited, i.e. gene expression instead of dist matrix)\n"); 
	fprintf(stderr, "\t\t\t-dist <correlation|euclidian> (distance measure to analyze raw data matrix with, def: corr)\n"); 
	fprintf(stderr, "\t\t\t-sub <#> (cluster using subset, then adds the rest of the data back,def=7000)\n"); 
	fprintf(stderr, "\t\t-o <outputfile prefix> (prefix for output files, def: out)\n"); 
	fprintf(stderr, "\t\t-r (reverse distance metric = high values are \"closer\")\n"); 
	fprintf(stderr, "\n\tExtract Clusters from output gtr or atr format file:\n"); 
	fprintf(stderr, "\t\t-gtr <gtr file> (output from clustering - cdt must be present too)\n"); 
	fprintf(stderr, "\t\t-thresh <#> (threshold to define clusters)\n"); 
	fprintf(stderr, "\n\tScore regions against clusters (use either -gtr/-d to supply distance matrix):\n");
	fprintf(stderr, "\t\t-c <cluster file> [cluster file 2] ... (First column is region name)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}


void clusterProgram(int argc, char** argv) {

	fprintf(stderr, "\n");
	if (argc < 3) {
		printCMDcluster();
	}

	char* distMatrixFile = NULL;
	char* dataMatrixFile = NULL;
	char* prefix = NULL;
	char** clusterFiles = new char*[10000];
	int numClusterFiles = 0;
	int reverseFlag = 0;
	int numAnnCols = 1;
	int distanceMetric = CLUSTER_DISTANCE_CORRELATION;
	int matrixOutputFlag = 1;
	int numSub = CLUSTER_SUB_INIT;

	char* gtrInputFile = NULL;
	double clusterThreshold = 0.0;

	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-d")==0) {
			distMatrixFile = argv[++i];
		} else if (strcmp(argv[i],"-i")==0) {
			dataMatrixFile = argv[++i];
			matrixOutputFlag = 0;
		} else if (strcmp(argv[i],"-o")==0) {
			prefix = argv[++i];
		} else if (strcmp(argv[i],"-gtr")==0) {
			gtrInputFile = argv[++i];
		} else if (strcmp(argv[i],"-thresh")==0) {
			sscanf(argv[++i],"%lf",&clusterThreshold);
		} else if (strcmp(argv[i],"-annCols")==0) {
			sscanf(argv[++i],"%d",&numAnnCols);
		} else if (strcmp(argv[i],"-sub")==0) {
			sscanf(argv[++i],"%d",&numSub);
			matrixOutputFlag = 0;
		} else if (strcmp(argv[i],"-c")==0) {
			i++;
			while (i<argc && argv[i][0] != '-') {
				clusterFiles[numClusterFiles++]=argv[i];
				i++;
			}
			if (i<argc && argv[i][0] == '-') i--;
		} else if (strcmp(argv[i],"-dist")==0 || strcmp(argv[i],"-distance")==0) {
			i++;
			if (strcmp(argv[i],"corr")==0 || strcmp(argv[i],"correlation")==0) {
				distanceMetric = CLUSTER_DISTANCE_CORRELATION;
			} else if (strcmp(argv[i],"euclidean")==0) {
				distanceMetric = CLUSTER_DISTANCE_EUCLIDEAN;
			} else {
				fprintf(stderr, "!!! Error - %s not a proper option for \"-distance\"\n",argv[i]);
				exit(0);
			}
		} else if (strcmp(argv[i],"-r")==0) {
			reverseFlag = 1;
		} else {
			fprintf(stderr, "!!! Invalid option: %s\n", argv[i]);
			printCMDcluster();
		}
	}

	FILE* clusterScoreFP = stdout;

	TreeCluster* cluster = NULL;

	if (gtrInputFile != NULL) {
		cluster = new TreeCluster();
		cluster->loadGTR(gtrInputFile);

		if (numClusterFiles < 1) {
			cluster->printClusters(stdout,NULL,clusterThreshold);
			exit(0);
		} else {
			cluster->scoreClusterFiles(clusterFiles, numClusterFiles,clusterScoreFP);
			exit(0);
		}
	}

	cluster = new TreeCluster();
	if (dataMatrixFile != NULL) {
		cluster->readDataMatrix(dataMatrixFile,numAnnCols);
		if (cluster->numMatrix > numSub) {
			cluster->sampleDataMatrix(numSub);
		}
	}
	if (distMatrixFile != NULL) {
		cluster->readDistMatrix(distMatrixFile,numAnnCols);
	}
	if (cluster->matrix == NULL) {
		cluster->calcDistMatrix(distanceMetric);
	}
	if (reverseFlag) {
		cluster->invertDistMatrix();
	}
	//cluster = new TreeCluster(distMatrix,numMatrix);
	if (numClusterFiles > 0) {
		cluster->scoreClusterFiles(clusterFiles, numClusterFiles,clusterScoreFP);
	} else {
		cluster->cluster();
		if (cluster->subFlag) {
			cluster->expandClusters();
		}
		cluster->printCDT(prefix,matrixOutputFlag);
	}
	exit(0);
}

//----------------------- matrix3D program

void printCMDmatrix3D() {
	fprintf(stderr, "\n\tUsage: homerTools matrix3D <distance matrix file> [options]\n");
	fprintf(stderr, "\n\tOptions for command: cluster\n");
	fprintf(stderr, "\t\t-r (invert the scores i.e. high scores = close together)\n");
	fprintf(stderr, "\t\t-2 (input matrix has two columns of identifiers [to work well with treeview])\n");
	fprintf(stderr, "\t\t-c <cluster assignment file> (color clusters)\n");
	fprintf(stderr, "\t\t-iter <#> (number of iterations)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}

void matrix3DProgram(int argc, char** argv) {

	fprintf(stderr, "\n");
	if (argc < 3) {
		printCMDmatrix3D();
	}

	srand(time(NULL));

	char* distMatrixFile = argv[2];
	//char* prefix = NULL;
	//int reverseFlag = 0;
	int col2Flag = 0;
	char* clusterFile = NULL;
	double iters = 1e6;

	for (int i=3;i<argc;i++) {
		if (strcmp(argv[i],"-r")==0) {
			//reverseFlag = 1;
		} else if (strcmp(argv[i],"-2")==0) {
			col2Flag = 1;
		} else if (strcmp(argv[i],"-c")==0) {
			clusterFile = argv[++i];
		} else if (strcmp(argv[i],"-iter")==0) {
			sscanf(argv[++i],"%lf",&iters);
		} else {
			fprintf(stderr, "!!! Invalid option: %s\n", argv[i]);
			printCMDcluster();
		}
	}

	int maxIterations = (int)iters; 
	char* buf = new char[BUFFER];
	char** line = new char*[10000];
	int numCols = 0;

	FILE* fp = fopen(distMatrixFile,"r");
	if (fp == NULL) {
		fprintf(stderr, "!! Could not open input file: %s\n", distMatrixFile);
		exit(0);
	}


	int numMatrix = 0;
	double** distMatrix = NULL;
	char** names = NULL;
	double fixedStep = -1;
	
	int lineCount = 0;
	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf, line, numCols, '\t');
		lineCount++;
		if (lineCount == 1) {
			numMatrix = numCols-1-col2Flag;
			distMatrix = new double*[numMatrix];
			names = new char*[numMatrix];
			for (int i=0;i<numMatrix;i++) {
				distMatrix[i] = new double[numMatrix];
				for (int j=0;j<numMatrix;j++) distMatrix[i][j]=FLT_MAX;
			}
			continue;
		}
		if (numCols < numMatrix+1) {
			fprintf(stderr, "Something is wrong with your distance matrix!!!\n");
			exit(0);
		}
		int index = lineCount -2;
		names[index] = new char[strlen(line[0])+1];
		strcpy(names[index], line[0]);
		for (int i=1+col2Flag;i<numCols;i++) {
			sscanf(line[i], "%lf", &(distMatrix[index][i-1-col2Flag]));
		}
	}




	Genome3D* genomeStruct = new Genome3D();
	genomeStruct->init(distMatrix,numMatrix,names,fixedStep);
	if (clusterFile != NULL)
		genomeStruct->addColorsFromClusters(clusterFile);

	genomeStruct->optimize(maxIterations);
	genomeStruct->print(stdout);


}



//---------------------------- filterids / bgbins ---------------------------
// Byte-exact C++ forms of filterListBy.pl <data> <list> 0 1 and of
// cat <target bins> <filtered background bins> | makeBinaryFile.pl - <target bins>,
// as findMotifsGenome.pl runs them. Perl's line handling is reproduced: a line
// ends at "\n" (a last line without one still counts), every "\r" in it is
// removed, and its ID is the text before the first tab (empty for an empty
// line). split /\t/ drops trailing empty fields, so filterListBy.pl prints a
// kept line with its trailing tabs removed. An input that cannot be opened
// reads as empty, as Perl's unchecked open does.

#if defined(__GNUC__) || defined(__clang__)
#define FF_PREFETCH(p) __builtin_prefetch(p)
#else
#define FF_PREFETCH(p) ((void)0)
#endif

class FFLineReader {
public:
	FILE* fp;
	std::vector<char> buf;
	size_t pos, end;
	bool eof;
	FFLineReader(const char* path) {
		fp = fopen(path, "rb");
		buf.resize(1 << 22);
		pos = end = 0;
		eof = (fp == NULL);
	}
	~FFLineReader() {
		if (fp != NULL) fclose(fp);
	}
	// the line stays valid until the next call
	bool next(const char* &line, size_t &len) {
		while (1) {
			if (pos < end) {
				const char* nl = (const char*)memchr(&(buf[pos]), '\n', end-pos);
				if (nl != NULL) {
					line = &(buf[pos]);
					len = nl - line;
					pos += len + 1;
					return true;
				}
			}
			if (eof) {
				if (pos < end) {
					line = &(buf[pos]);
					len = end-pos;
					pos = end;
					return true;
				}
				return false;
			}
			if (pos > 0) {
				if (end > pos) memmove(&(buf[0]), &(buf[pos]), end-pos);
				end -= pos;
				pos = 0;
			}
			if (end == buf.size()) buf.resize(buf.size()*2);
			size_t n = fread(&(buf[end]), 1, buf.size()-end, fp);
			end += n;
			if (n == 0) eof = true;
		}
	}
};

static size_t ffFileSize(const char* path) {
	struct stat st;
	if (stat(path, &st) != 0) return 0;
	return (size_t)st.st_size;
}

static void ffAppendFile(const char* path, std::string &out) {
	FILE* fp = fopen(path, "rb");
	if (fp == NULL) return;
	out.reserve(out.size() + ffFileSize(path) + 1);
	char rbuf[1 << 16];
	size_t n = 0;
	while ((n = fread(rbuf, 1, sizeof(rbuf), fp)) > 0) out.append(rbuf, n);
	fclose(fp);
}

static inline uint64_t ffHash(const char* p, size_t n) {
	uint64_t h = 0x9E3779B97F4A7C15ULL ^ ((uint64_t)n * 0xC2B2AE3D27D4EB4FULL);
	while (n >= 8) {
		uint64_t w;
		memcpy(&w, p, 8);
		h = (h ^ w) * 0xFF51AFD7ED558CCDULL;
		h ^= h >> 29;
		p += 8;
		n -= 8;
	}
	uint64_t w = 0;
	memcpy(&w, p, n);
	h = (h ^ w) * 0xC4CEB9FE1A85EC53ULL;
	h ^= h >> 32;
	h *= 0xFF51AFD7ED558CCDULL;
	h ^= h >> 29;
	return h;
}

struct FFSlot {
	uint32_t tag; // 0 = empty
	uint32_t len;
	const char* ptr;
};

// Set of byte strings (open addressing, linear probing). A key is copied into
// the set's own chunks, or referenced where it outlives the set. Slot numbers
// change when the table grows.
class FFKeySet {
public:
	std::vector<FFSlot> slots;
	std::vector<char*> chunks;
	size_t chunkUsed, chunkSize;
	size_t mask, count;
	FFKeySet() {
		slots.resize(1 << 10);
		memset(&(slots[0]), 0, slots.size()*sizeof(FFSlot));
		mask = slots.size()-1;
		count = chunkUsed = chunkSize = 0;
	}
	~FFKeySet() {
		for (size_t i=0;i<chunks.size();i++) delete [](chunks[i]);
	}
	static inline uint32_t tagOf(uint64_t h) {
		return (uint32_t)(h >> 32) | 1;
	}
	void reserve(size_t n) {
		size_t cap = slots.size();
		while (cap*7 < n*10) cap *= 2;
		if (cap != slots.size()) rehash(cap);
	}
	void rehash(size_t cap) {
		std::vector<FFSlot> old;
		old.swap(slots);
		slots.resize(cap);
		memset(&(slots[0]), 0, cap*sizeof(FFSlot));
		mask = cap-1;
		for (size_t i=0;i<old.size();i++) {
			if (old[i].tag == 0) continue;
			size_t s = ffHash(old[i].ptr, old[i].len) & mask;
			while (slots[s].tag != 0) s = (s+1) & mask;
			slots[s] = old[i];
		}
	}
	const char* store(const char* p, size_t n) {
		if (chunkUsed + n > chunkSize || chunks.size() == 0) {
			chunkSize = (n > (1 << 22)) ? n : (1 << 22);
			chunks.push_back(new char[chunkSize]);
			chunkUsed = 0;
		}
		char* d = chunks.back() + chunkUsed;
		if (n > 0) memcpy(d, p, n);
		chunkUsed += n;
		return d;
	}
	inline void prefetch(uint64_t h) {
		FF_PREFETCH(&(slots[h & mask]));
	}
	// the first slot at or after s that is empty or holds a key of this tag and length
	inline size_t probe(size_t s, uint32_t tag, size_t n) {
		while (slots[s].tag != 0 && (slots[s].tag != tag || slots[s].len != n)) s = (s+1) & mask;
		return s;
	}
	// slot of the key, or -1
	long long find(const char* p, size_t n, uint64_t h) {
		uint32_t tag = tagOf(h);
		size_t s = probe(h & mask, tag, n);
		while (slots[s].tag != 0) {
			if (memcmp(slots[s].ptr, p, n) == 0) return (long long)s;
			s = probe((s+1) & mask, tag, n);
		}
		return -1;
	}
	long long find(const char* p, size_t n) {
		return find(p, n, ffHash(p, n));
	}
	// slot of the key (valid until the table grows); isNew says whether it was added
	size_t insert(const char* p, size_t n, uint64_t h, bool copy, bool &isNew) {
		if (n > 0xFFFFFFFFULL) {
			fprintf(stderr, "!!! ID too long !!!\n");
			exit(1);
		}
		long long f = find(p, n, h);
		if (f >= 0) {
			isNew = false;
			return (size_t)f;
		}
		isNew = true;
		if ((count+1)*10 > slots.size()*7) rehash(slots.size()*2);
		size_t s = h & mask;
		while (slots[s].tag != 0) s = (s+1) & mask;
		slots[s].tag = tagOf(h);
		slots[s].len = (uint32_t)n;
		slots[s].ptr = copy ? store(p, n) : p;
		count++;
		return s;
	}
	size_t insert(const char* p, size_t n, bool copy, bool &isNew) {
		return insert(p, n, ffHash(p, n), copy, isNew);
	}
};

// The ID of a line: the text before the first tab, with every "\r" removed.
// When a "\r" is removed the ID is built in scratch.
static inline void ffLineKey(const char* line, size_t len, std::string &scratch,
						const char* &key, size_t &klen) {
	const char* tab = (const char*)memchr(line, '\t', len);
	klen = (tab != NULL) ? (size_t)(tab-line) : len;
	key = line;
	if (klen == 0 || memchr(line, '\r', klen) == NULL) return;
	scratch.clear();
	for (size_t i=0;i<klen;i++) if (line[i] != '\r') scratch.push_back(line[i]);
	key = scratch.data();
	klen = scratch.size();
}

// filterListBy.pl's printed form of a kept line, without the "\n": "\r"
// removed, trailing tabs dropped.
static inline void ffAppendFields(std::string &out, const char* line, size_t len) {
	size_t start = out.size();
	if (len > 0 && memchr(line, '\r', len) != NULL) {
		for (size_t i=0;i<len;i++) if (line[i] != '\r') out.push_back(line[i]);
	} else {
		out.append(line, len);
	}
	size_t e = out.size();
	while (e > start && out[e-1] == '\t') e--;
	out.resize(e);
}

static void ffFlush(std::string &out, FILE* fp, size_t keep) {
	if (out.size() < keep) return;
	if (out.size() > 0) fwrite(out.data(), 1, out.size(), fp);
	out.clear();
}

static void ffReadIds(const char* path, FFKeySet &ids) {
	ids.reserve(ffFileSize(path)/48);
	FFLineReader in(path);
	const int B = 16;
	const char* kp[B];
	size_t kn[B];
	uint64_t kh[B];
	std::string scratch;
	bool more = true, isNew = false;
	while (more) {
		int nb = 0;
		const char* line = NULL;
		size_t len = 0;
		while (nb < B && (more = in.next(line, len))) {
			const char* key = NULL;
			size_t klen = 0;
			ffLineKey(line, len, scratch, key, klen);
			kp[nb] = ids.store(key, klen);
			kn[nb] = klen;
			kh[nb] = ffHash(kp[nb], klen);
			ids.prefetch(kh[nb]);
			nb++;
		}
		for (int i=0;i<nb;i++) ids.insert(kp[i], kn[i], kh[i], false, isNew);
	}
}

void printCMDfilterids() {
	fprintf(stderr, "\n\tUsage: homerTools filterids <data file> <ID list file>\n");
	fprintf(stderr, "\n\tPrints to stdout the lines of <data file> whose first tab-separated column is the\n");
	fprintf(stderr, "\tfirst column of a line in <ID list file>, with carriage returns and trailing tabs\n");
	fprintf(stderr, "\tremoved: the output of filterListBy.pl <data file> <ID list file> 0 1\n");
	fprintf(stderr, "\n");
	exit(0);
}

void filterIdsProgram(int argc, char** argv) {
	if (argc < 4) printCMDfilterids();
	FFKeySet ids;
	ffReadIds(argv[3], ids);
	FFLineReader in(argv[2]);
	const char* line = NULL;
	size_t len = 0;
	std::string out, scratch;
	while (in.next(line, len)) {
		const char* key = NULL;
		size_t klen = 0;
		ffLineKey(line, len, scratch, key, klen);
		if (ids.find(key, klen) >= 0) {
			ffAppendFields(out, line, len);
			out.push_back('\n');
			ffFlush(out, stdout, 1 << 22);
		}
	}
	ffFlush(out, stdout, 0);
	fflush(stdout);
}

// The ID list of `mergePeaks -cobound 1 -prefix <prefix> [-strand] <pos> <peaks>`'s
// coBoundBy0.txt, as filterListBy.pl reads it: the IDs of the header line, and
// the name of every <pos> peak that overlaps no <peaks> peak (closed intervals,
// same chromosome, and same strand with -strand). Peak lines are classified
// as PeakLibrary::readPeakFile classifies them. Returns false, having read
// only, when mergePeaks could behave in a way this does not reproduce: an
// argument mergePeaks reads as an option or as the same file, an output it
// cannot open, a NUL or carriage return, a line longer than its fgets buffer,
// a BED-format line, a coordinate sscanf would not read exactly, a peak with
// start > end, a repeated name (mergePeaks renames it), target names longer in
// total than its overlap string, or a file with no peaks. ids holds every
// <pos> name and the header IDs; keep[slot] marks the coBoundBy0 IDs.

#define FF_PEAK_SKIP 0
#define FF_PEAK_OK 1
#define FF_PEAK_UNSUPPORTED 2
// the longest line SeqTag.cpp's fgets(buf, 100000) returns whole, with margin
#define FF_PEAK_MAX_LINE 99997

struct FFPeak {
	const char* name;
	size_t nameLen;
	const char* chr;
	size_t chrLen;
	int start, end, dir;
};

// SeqTag.cpp's checkInt on one column
static int ffCheckInt(const char* s, size_t n) {
	int bad = 0;
	for (size_t i=0;i<n;i++) {
		if (s[i] < 45) bad=1;
		if (s[i] == 46 || s[i] == 47) bad=1;
		if (s[i] > 57) bad=1;
	}
	return bad;
}

// SeqTag.cpp's checkStrand on one column
static int ffCheckStrand(const char* s, size_t n) {
	if (n != 1) return -1;
	if (s[0] == '0' || s[0] == '+' || s[0] == '.') return 0;
	if (s[0] == '1' || s[0] == '-') return 1;
	return -1;
}

// -?[0-9]{1,9}, which sscanf("%d") reads exactly
static bool ffPlainInt(const char* s, size_t n, int &v) {
	size_t i = 0;
	bool neg = false;
	if (n > 0 && s[0] == '-') {
		neg = true;
		i = 1;
	}
	if (n-i < 1 || n-i > 9) return false;
	int x = 0;
	for (;i<n;i++) {
		if (s[i] < '0' || s[i] > '9') return false;
		x = x*10 + (s[i]-'0');
	}
	v = neg ? -x : x;
	return true;
}

// One line (no "\n", "\r" or NUL) as PeakLibrary::readPeakFile in
// PEAK_READ_MODE_NORMAL reads it. Only columns 0-4 decide, and a line with
// five or more columns is treated alike.
static int ffPeakLine(const char* s, size_t len, FFPeak &pk) {
	const char* c[5];
	size_t cn[5];
	int numCols = 0;
	const char* e = s + len;
	const char* p = s;
	while (numCols < 5) {
		const char* t = (const char*)memchr(p, '\t', e-p);
		c[numCols] = p;
		cn[numCols] = ((t != NULL) ? t : e) - p;
		numCols++;
		if (t == NULL) break;
		p = t+1;
	}
	if (numCols < 3) return FF_PEAK_SKIP;
	if (cn[0] > 0 && c[0][0] == '#') return FF_PEAK_SKIP;
	char c2 = (cn[2] > 0) ? c[2][0] : '\0';
	if (c2 != '-' && (c2 < 48 || c2 > 57)) return FF_PEAK_SKIP;
	int col1notNum = ffCheckInt(c[1], cn[1]);
	if (ffCheckInt(c[2], cn[2])) return FF_PEAK_SKIP;
	int dir = -1;
	if (numCols == 3) {
		if (col1notNum) return FF_PEAK_SKIP;
		return FF_PEAK_UNSUPPORTED;
	}
	int col3notNum = ffCheckInt(c[3], cn[3]);
	if (col1notNum && col3notNum) return FF_PEAK_SKIP;
	if (numCols == 4) {
		if (col3notNum) return FF_PEAK_UNSUPPORTED;
	} else {
		int col4strand = ffCheckStrand(c[4], cn[4]);
		if (col4strand != -1 && col3notNum == 0) {
			dir = col4strand;
		} else if (col3notNum) {
			return FF_PEAK_UNSUPPORTED;
		} else {
			return FF_PEAK_SKIP;
		}
	}
	if (!ffPlainInt(c[2], cn[2], pk.start) || !ffPlainInt(c[3], cn[3], pk.end)) return FF_PEAK_UNSUPPORTED;
	if (pk.start > pk.end) return FF_PEAK_UNSUPPORTED;
	pk.name = c[0];
	pk.nameLen = cn[0];
	pk.chr = c[1];
	pk.chrLen = cn[1];
	pk.dir = (dir == -1) ? 0 : dir;
	return FF_PEAK_OK;
}

// the next line of buf from p ("\n" removed; a last line without one counts)
static inline bool ffNextLine(const std::string &buf, size_t &p, const char* &line, size_t &len) {
	if (p >= buf.size()) return false;
	line = buf.data() + p;
	const char* nl = (const char*)memchr(line, '\n', buf.size()-p);
	len = (nl != NULL) ? (size_t)(nl-line) : buf.size()-p;
	p += len + 1;
	return true;
}

struct FFTargetChr {
	std::vector<std::pair<int,int> > iv[2];
	std::vector<int> start[2];
	std::vector<int> maxEnd[2];
};

static bool ffCoboundIds(const char* posFile, const char* peaksFile, const char* prefix,
				bool strandSep, std::string &posBuf, FFKeySet &ids, std::vector<char> &keep) {
	if (posFile[0] == '-' || peaksFile[0] == '-' || strcmp(posFile, peaksFile) == 0) return false;

	std::string tbuf;
	ffAppendFile(peaksFile, tbuf);
	if (memchr(tbuf.data(), '\0', tbuf.size()) != NULL || memchr(tbuf.data(), '\r', tbuf.size()) != NULL) return false;
	std::unordered_map<std::string,int> chrIndex;
	std::vector<FFTargetChr> tchr;
	FFKeySet tnames;
	size_t nameBytes = 0, numTargets = 0, p = 0;
	const char* line = NULL;
	size_t len = 0;
	bool isNew = false;
	while (ffNextLine(tbuf, p, line, len)) {
		if (len > FF_PEAK_MAX_LINE) return false;
		FFPeak pk;
		int r = ffPeakLine(line, len, pk);
		if (r == FF_PEAK_SKIP) continue;
		if (r == FF_PEAK_UNSUPPORTED) return false;
		tnames.insert(pk.name, pk.nameLen, false, isNew);
		if (!isNew) return false;
		// getCoBoundPeaks joins a peak's overlapping target names into a 1,000,000-byte string
		nameBytes += pk.nameLen + 1;
		if (nameBytes > 900000) return false;
		std::string chr(pk.chr, pk.chrLen);
		std::unordered_map<std::string,int>::iterator it = chrIndex.find(chr);
		int ci = 0;
		if (it == chrIndex.end()) {
			ci = (int)tchr.size();
			chrIndex[chr] = ci;
			tchr.push_back(FFTargetChr());
		} else {
			ci = it->second;
		}
		tchr[ci].iv[strandSep ? pk.dir : 0].push_back(std::make_pair(pk.start, pk.end));
		numTargets++;
	}
	if (numTargets == 0) return false;
	for (size_t i=0;i<tchr.size();i++) {
		for (int k=0;k<2;k++) {
			std::vector<std::pair<int,int> > &iv = tchr[i].iv[k];
			std::sort(iv.begin(), iv.end());
			int m = INT_MIN;
			for (size_t j=0;j<iv.size();j++) {
				if (iv[j].second > m) m = iv[j].second;
				tchr[i].start[k].push_back(iv[j].first);
				tchr[i].maxEnd[k].push_back(m);
			}
		}
	}

	std::string cmd = std::string("mergePeaks -cobound 1 -prefix ") + prefix
					+ (strandSep ? " -strand " : " ") + posFile + " " + peaksFile;
	std::string header = "#Name (cmd = " + cmd + ")\tchr\tstart\tend\tstrand\tstat\tstat2\t"
					+ peaksFile + "\n";
	size_t headerLines = 0;
	for (size_t i=0;i<header.size();i++) if (header[i] == '\n') headerLines++;

	ffAppendFile(posFile, posBuf);
	if (memchr(posBuf.data(), '\0', posBuf.size()) != NULL || memchr(posBuf.data(), '\r', posBuf.size()) != NULL) return false;
	size_t numLines = 1;
	for (const char* q = posBuf.data(); (q = (const char*)memchr(q, '\n', posBuf.data()+posBuf.size()-q)) != NULL; q++) numLines++;
	ids.reserve(numLines + headerLines + 16);
	keep.assign(ids.slots.size(), 0);

	const int B = 16;
	const char* kp[B];
	size_t kn[B];
	uint64_t kh[B];
	char ovl[B];
	size_t numRef = 0;
	const char* lastChr = NULL;
	size_t lastChrLen = 0;
	int lastCi = -1;
	p = 0;
	bool more = true;
	while (more) {
		int nb = 0;
		while (nb < B && (more = ffNextLine(posBuf, p, line, len))) {
			if (len > FF_PEAK_MAX_LINE) return false;
			FFPeak pk;
			int r = ffPeakLine(line, len, pk);
			if (r == FF_PEAK_SKIP) continue;
			if (r == FF_PEAK_UNSUPPORTED) return false;
			if (lastChr == NULL || pk.chrLen != lastChrLen || memcmp(pk.chr, lastChr, lastChrLen) != 0) {
				std::unordered_map<std::string,int>::iterator it = chrIndex.find(std::string(pk.chr, pk.chrLen));
				lastCi = (it == chrIndex.end()) ? -1 : it->second;
				lastChr = pk.chr;
				lastChrLen = pk.chrLen;
			}
			char o = 0;
			if (lastCi >= 0) {
				int k = strandSep ? pk.dir : 0;
				const std::vector<int> &st = tchr[lastCi].start[k];
				size_t j = std::upper_bound(st.begin(), st.end(), pk.end) - st.begin();
				if (j > 0 && tchr[lastCi].maxEnd[k][j-1] >= pk.start) o = 1;
			}
			kp[nb] = pk.name;
			kn[nb] = pk.nameLen;
			kh[nb] = ffHash(pk.name, pk.nameLen);
			ovl[nb] = o;
			ids.prefetch(kh[nb]);
			nb++;
		}
		for (int i=0;i<nb;i++) {
			size_t s = ids.insert(kp[i], kn[i], kh[i], false, isNew);
			if (!isNew) return false;
			keep[s] = !ovl[i];
		}
		numRef += nb;
	}
	if (numRef == 0) return false;

	std::string scratch;
	size_t hp = 0;
	while (ffNextLine(header, hp, line, len)) {
		const char* key = NULL;
		size_t klen = 0;
		ffLineKey(line, len, scratch, key, klen);
		size_t s = ids.insert(key, klen, true, isNew);
		keep[s] = 1;
	}
	if (keep.size() != ids.slots.size()) {
		fprintf(stderr, "!!! bgbins: ID table grew past its reserve !!!\n");
		exit(1);
	}
	// mergePeaks creates its two outputs last, and fails when it cannot
	for (int i=0;i<2;i++) {
		std::string fname = std::string(prefix) + ((i == 0) ? ".coBoundBy0.txt" : ".coBoundBy1.txt");
		FILE* fp = fopen(fname.c_str(), "w");
		if (fp == NULL) return false;
		fclose(fp);
	}
	return true;
}

void printCMDbgbins() {
	fprintf(stderr, "\n\tUsage: homerTools bgbins -target <target bins> -bg <background bins> -ids <ID list file> -cat <output file>\n");
	fprintf(stderr, "\t   or: homerTools bgbins -target <target bins> -bg <background bins> -pos <pos file> -peaks <peak file>\n");
	fprintf(stderr, "\t                         -prefix <prefix> [-strand] -cat <output file>\n");
	fprintf(stderr, "\n\tWrites to -cat the target bins file followed by the lines of the background bins\n");
	fprintf(stderr, "\tfile whose ID is in -ids, as filterListBy.pl <bg> <ids> 0 1 prints them, and prints\n");
	fprintf(stderr, "\tto stdout what makeBinaryFile.pl <-cat file> <target bins> prints.\n");
	fprintf(stderr, "\tWith -pos/-peaks the ID list is the <prefix>.coBoundBy0.txt file that\n");
	fprintf(stderr, "\tmergePeaks -cobound 1 -prefix <prefix> [-strand] <pos file> <peak file> would write\n");
	fprintf(stderr, "\t(the two coBoundBy files are created empty); exits with status 2, having written\n");
	fprintf(stderr, "\tnothing, when the inputs need mergePeaks itself.\n");
	fprintf(stderr, "\n");
	exit(0);
}

void bgBinsProgram(int argc, char** argv) {
	const char* targetFile = NULL;
	const char* bgFile = NULL;
	const char* idFile = NULL;
	const char* catFile = NULL;
	const char* posFile = NULL;
	const char* peaksFile = NULL;
	const char* prefix = NULL;
	bool strandSep = false;
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-strand")==0) {
			strandSep = true;
			continue;
		}
		if (i+1 >= argc) printCMDbgbins();
		if (strcmp(argv[i],"-target")==0) {
			targetFile = argv[++i];
		} else if (strcmp(argv[i],"-bg")==0) {
			bgFile = argv[++i];
		} else if (strcmp(argv[i],"-ids")==0) {
			idFile = argv[++i];
		} else if (strcmp(argv[i],"-cat")==0) {
			catFile = argv[++i];
		} else if (strcmp(argv[i],"-pos")==0) {
			posFile = argv[++i];
		} else if (strcmp(argv[i],"-peaks")==0) {
			peaksFile = argv[++i];
		} else if (strcmp(argv[i],"-prefix")==0) {
			prefix = argv[++i];
		} else {
			printCMDbgbins();
		}
	}
	if (targetFile == NULL || bgFile == NULL || catFile == NULL) printCMDbgbins();
	bool coboundFlag = (posFile != NULL || peaksFile != NULL || prefix != NULL);
	if (coboundFlag && (posFile == NULL || peaksFile == NULL || prefix == NULL || idFile != NULL)) printCMDbgbins();
	if (!coboundFlag && idFile == NULL) printCMDbgbins();

	// built before anything is written, so that status 2 leaves no output
	FFKeySet ids;
	std::vector<char> keep;
	std::string posBuf;
	if (coboundFlag) {
		if (!ffCoboundIds(posFile, peaksFile, prefix, strandSep, posBuf, ids, keep)) exit(2);
	} else {
		ffReadIds(idFile, ids);
		keep.assign(ids.slots.size(), 1);
	}

	// makeBinaryFile.pl <target + kept background> <target> prints each first
	// occurrence of an ID in the combined lines, with 1 if it is a target ID
	// (the ID of any target line) and 0 if not, then the target IDs never met.
	// Target lines come first, so every target ID but the last line's is met
	// there. The last target line, when it has no "\n", merges with the first
	// kept background line.
	std::string combined, out, scratch;
	ffAppendFile(targetFile, combined);
	size_t targetBytes = combined.size();
	FFKeySet targetDone;
	bool isNew = false;
	size_t p = 0;
	while (p < targetBytes) {
		const char* nl = (const char*)memchr(combined.data()+p, '\n', targetBytes-p);
		if (nl == NULL) break;
		size_t len = nl - (combined.data()+p);
		const char* key = NULL;
		size_t klen = 0;
		ffLineKey(combined.data()+p, len, scratch, key, klen);
		targetDone.insert(key, klen, true, isNew);
		if (isNew) {
			out.append(key, klen);
			out.append("\t1\n");
		}
		p += len + 1;
	}
	ffFlush(out, stdout, 0);
	bool hasTail = (p < targetBytes);
	std::string tail, tailKey;
	bool tailSeen = false;
	if (hasTail) {
		tail.assign(combined.data()+p, targetBytes-p);
		const char* key = NULL;
		size_t klen = 0;
		ffLineKey(tail.data(), tail.size(), scratch, key, klen);
		tailKey.assign(key, klen);
		tailSeen = (targetDone.find(tailKey.data(), tailKey.size()) >= 0);
	}

	// done[slot]: the ID was already printed
	std::vector<char> done(ids.slots.size(), 0);
	for (size_t i=0;i<targetDone.slots.size();i++) {
		if (targetDone.slots[i].tag == 0) continue;
		long long s = ids.find(targetDone.slots[i].ptr, targetDone.slots[i].len);
		if (s >= 0) done[s] = 1;
	}

	std::string bg;
	ffAppendFile(bgFile, bg);
	combined.reserve(targetBytes + bg.size() + 1);
	bool first = true;
	const int B = 16;
	const char* lp[B];
	size_t ln[B];
	const char* kp[B];
	size_t kn[B];
	uint64_t kh[B];
	size_t cand[B];
	std::string kscratch[B];
	const char* bp = bg.data();
	const char* bend = bp + bg.size();
	while (bp < bend) {
		int nb = 0;
		while (nb < B && bp < bend) {
			const char* nl = (const char*)memchr(bp, '\n', bend-bp);
			size_t len = (nl != NULL) ? (size_t)(nl-bp) : (size_t)(bend-bp);
			lp[nb] = bp;
			ln[nb] = len;
			ffLineKey(bp, len, kscratch[nb], kp[nb], kn[nb]);
			kh[nb] = ffHash(kp[nb], kn[nb]);
			ids.prefetch(kh[nb]);
			bp += len + 1;
			nb++;
		}
		for (int i=0;i<nb;i++) {
			cand[i] = ids.probe(kh[i] & ids.mask, FFKeySet::tagOf(kh[i]), kn[i]);
			if (ids.slots[cand[i]].tag != 0) FF_PREFETCH(ids.slots[cand[i]].ptr);
		}
		for (int i=0;i<nb;i++) {
			size_t s = cand[i];
			while (ids.slots[s].tag != 0 && memcmp(ids.slots[s].ptr, kp[i], kn[i]) != 0) {
				s = ids.probe((s+1) & ids.mask, FFKeySet::tagOf(kh[i]), kn[i]);
			}
			if (ids.slots[s].tag == 0 || !keep[s]) continue;
			size_t start = combined.size();
			ffAppendFields(combined, lp[i], ln[i]);
			combined.push_back('\n');
			if (first && hasTail) {
				// the merged line: the last target line, then this one
				std::string merged = tail;
				merged.append(combined.data()+start, combined.size()-1-start);
				const char* key = NULL;
				size_t klen = 0;
				ffLineKey(merged.data(), merged.size(), scratch, key, klen);
				targetDone.insert(key, klen, true, isNew);
				if (isNew) {
					bool special = (klen == tailKey.size() && memcmp(key, tailKey.data(), klen) == 0);
					if (special) tailSeen = true;
					out.append(key, klen);
					out.append(special ? "\t1\n" : "\t0\n");
					long long m = ids.find(key, klen);
					if (m >= 0) done[m] = 1;
				}
			} else if (done[s] == 0) {
				done[s] = 1;
				bool special = (hasTail && kn[i] == tailKey.size() && memcmp(kp[i], tailKey.data(), kn[i]) == 0);
				if (special) tailSeen = true;
				out.append(kp[i], kn[i]);
				out.append(special ? "\t1\n" : "\t0\n");
			}
			first = false;
			ffFlush(out, stdout, 1 << 22);
		}
	}
	if (hasTail && first) {
		// no background line kept: the last target line ends the combined file
		if (targetDone.find(tailKey.data(), tailKey.size()) < 0) {
			out.append(tailKey);
			out.append("\t1\n");
		}
		tailSeen = true;
	}
	if (hasTail && !tailSeen) {
		out.append(tailKey);
		out.append("\t1\n");
	}
	ffFlush(out, stdout, 0);
	fflush(stdout);

	FILE* fp = fopen(catFile, "wb");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open %s for writing !!!\n", catFile);
		exit(1);
	}
	if (combined.size() > 0) fwrite(combined.data(), 1, combined.size(), fp);
	fclose(fp);
}



//---------------------------- randbg ---------------------------
// randRemoveBackground.pl <group file> <total> <bins file | null> in C++, drawing
// Perl's rand() stream. Perl 5.20 and later draw rand() from their own
// drand48 on every platform: state = (state*0x5DEECE66D + 0xB) mod 2^48, value =
// state/2^48. -perlrand gives the first two values a fresh perl process draws,
// times 2^48 (findMotifsGenome.pl asks perl for them), which is the whole state.
// randbg exits with status 2, having printed nothing, when it cannot reproduce
// the script: the two values are not consecutive drand48 states, <total> is
// not a plain integer, a group line whose ID is in the bins has a second column
// other than 0 or 1, such an ID repeats, or its bin is not a plain integer
// (Perl would order bins of equal value by hash order).

static const uint64_t FF_DRAND_A = 0x5DEECE66DULL;
static const uint64_t FF_DRAND_C = 0xBULL;
static const uint64_t FF_DRAND_MASK = (1ULL << 48) - 1;
static const long long FF_BIN_NONE = LLONG_MIN;

class FFPerlRand {
public:
	uint64_t x;
	bool first;
	FFPerlRand(uint64_t x1) {
		x = x1;
		first = true;
	}
	inline double next() {
		if (first) {
			first = false;
		} else {
			x = (x*FF_DRAND_A + FF_DRAND_C) & FF_DRAND_MASK;
		}
		return ldexp((double)x, -48);
	}
};

// [0-9]{1,15}
static bool ffUnsignedArg(const char* s, unsigned long long &v) {
	size_t n = strlen(s);
	if (n < 1 || n > 15) return false;
	v = 0;
	for (size_t i=0;i<n;i++) {
		if (s[i] < '0' || s[i] > '9') return false;
		v = v*10 + (unsigned long long)(s[i]-'0');
	}
	return true;
}

// 0 or -?[1-9][0-9]{0,14}: distinct strings have distinct numeric values
static long long ffCanonicalInt(const char* s, size_t n) {
	if (n == 1 && s[0] == '0') return 0;
	size_t i = 0;
	bool neg = false;
	if (n > 0 && s[0] == '-') {
		neg = true;
		i = 1;
	}
	if (n-i < 1 || n-i > 15 || s[i] == '0') return FF_BIN_NONE;
	long long v = 0;
	for (;i<n;i++) {
		if (s[i] < '0' || s[i] > '9') return FF_BIN_NONE;
		v = v*10 + (s[i]-'0');
	}
	return neg ? -v : v;
}

// Perl's chomp; s/\r//g; split /\t/ for the first two fields. Carriage returns
// are already removed from buf. has1 is false when split gives no second field.
static inline void ffTwoFields(const char* line, size_t len, const char* &f0, size_t &n0,
							const char* &f1, size_t &n1, bool &has1) {
	const char* tab = (const char*)memchr(line, '\t', len);
	f0 = line;
	n0 = (tab != NULL) ? (size_t)(tab-line) : len;
	has1 = false;
	f1 = NULL;
	n1 = 0;
	if (tab == NULL) return;
	f1 = tab+1;
	const char* e = line + len;
	const char* tab2 = (const char*)memchr(f1, '\t', e-f1);
	n1 = ((tab2 != NULL) ? tab2 : e) - f1;
	has1 = true;
}

static void ffReadNoCR(const char* path, std::string &buf) {
	ffAppendFile(path, buf);
	if (memchr(buf.data(), '\r', buf.size()) != NULL) {
		buf.erase(std::remove(buf.begin(), buf.end(), '\r'), buf.end());
	}
}

struct FFRandBin {
	long long value;
	long long numPos;
	std::vector<std::pair<const char*, size_t> > back;
};

void printCMDrandbg() {
	fprintf(stderr, "\n\tUsage: homerTools randbg <group file> <total # of sequences> <bins file | null> -perlrand <x1> <x2>\n");
	fprintf(stderr, "\n\tPrints what randRemoveBackground.pl <group file> <total> <bins file | null> prints when\n");
	fprintf(stderr, "\tPerl's rand() would first return x1/2^48 and x2/2^48 (perl -e 'printf(\"%%.0f %%.0f\",\n");
	fprintf(stderr, "\trand()*2**48, rand()*2**48)' gives them). Exits with status 2, printing nothing, when\n");
	fprintf(stderr, "\tit cannot reproduce the script on these inputs.\n");
	fprintf(stderr, "\n");
	exit(0);
}

static int ffRandBg(const char* groupFile, const char* maxArg, const char* binsFile,
					const char* x1Arg, const char* x2Arg, std::string &out) {
	unsigned long long x1 = 0, x2 = 0, maxv = 0;
	if (!ffUnsignedArg(x1Arg, x1) || !ffUnsignedArg(x2Arg, x2)) return 2;
	if (x1 > FF_DRAND_MASK || x2 > FF_DRAND_MASK) return 2;
	if (((x1*FF_DRAND_A + FF_DRAND_C) & FF_DRAND_MASK) != x2) return 2;
	if (!ffUnsignedArg(maxArg, maxv)) return 2;

	std::string groupBuf, binsBuf;
	ffReadNoCR(groupFile, groupBuf);
	bool nullBins = (strcmp(binsFile, "null") == 0);
	if (!nullBins) ffReadNoCR(binsFile, binsBuf);
	const std::string &setBuf = nullBins ? groupBuf : binsBuf;

	// $bins{ID}: the last line's second field wins
	FFKeySet bins;
	size_t numLines = 1;
	for (const char* q = setBuf.data(); (q = (const char*)memchr(q, '\n', setBuf.data()+setBuf.size()-q)) != NULL; q++) numLines++;
	bins.reserve(numLines + 16);
	std::vector<long long> binValue(bins.slots.size(), FF_BIN_NONE);
	std::vector<char> used(bins.slots.size(), 0);
	const int B = 16;
	const char* kp[B];
	size_t kn[B];
	uint64_t kh[B];
	long long kv[B];
	const char* fp1[B];
	size_t fn1[B];
	bool fh1[B];
	size_t cand[B];
	const char* line = NULL;
	size_t len = 0, p = 0;
	bool isNew = false, more = true;
	while (more) {
		int nb = 0;
		while (nb < B && (more = ffNextLine(setBuf, p, line, len))) {
			const char *f1;
			size_t n1;
			bool has1;
			ffTwoFields(line, len, kp[nb], kn[nb], f1, n1, has1);
			kv[nb] = nullBins ? 1 : (has1 ? ffCanonicalInt(f1, n1) : FF_BIN_NONE);
			kh[nb] = ffHash(kp[nb], kn[nb]);
			bins.prefetch(kh[nb]);
			nb++;
		}
		for (int i=0;i<nb;i++) {
			size_t s = bins.insert(kp[i], kn[i], kh[i], false, isNew);
			binValue[s] = kv[i];
		}
	}
	if (binValue.size() != bins.slots.size()) {
		fprintf(stderr, "!!! randbg: ID table grew past its reserve !!!\n");
		exit(1);
	}

	std::vector<FFRandBin> binList;
	std::unordered_map<long long,int> binIndex;
	std::vector<int> smallBin(4096, -1);
	long long n = 0;
	p = 0;
	more = true;
	while (more) {
		int nb = 0;
		while (nb < B && (more = ffNextLine(groupBuf, p, line, len))) {
			ffTwoFields(line, len, kp[nb], kn[nb], fp1[nb], fn1[nb], fh1[nb]);
			kh[nb] = ffHash(kp[nb], kn[nb]);
			bins.prefetch(kh[nb]);
			nb++;
		}
		for (int i=0;i<nb;i++) {
			cand[i] = bins.probe(kh[i] & bins.mask, FFKeySet::tagOf(kh[i]), kn[i]);
			if (bins.slots[cand[i]].tag != 0) FF_PREFETCH(bins.slots[cand[i]].ptr);
		}
		for (int i=0;i<nb;i++) {
			const char* f0 = kp[i];
			size_t n0 = kn[i];
			const char* f1 = fp1[i];
			size_t n1 = fn1[i];
			bool has1 = fh1[i];
			size_t s = cand[i];
			while (bins.slots[s].tag != 0 && memcmp(bins.slots[s].ptr, f0, n0) != 0) {
				s = bins.probe((s+1) & bins.mask, FFKeySet::tagOf(kh[i]), n0);
			}
			if (bins.slots[s].tag == 0) continue;
			if (!has1 || n1 != 1 || (f1[0] != '0' && f1[0] != '1')) return 2;
			if (used[s]) return 2;
			used[s] = 1;
			long long v = binValue[s];
			if (v == FF_BIN_NONE) return 2;
			int bi = (v >= 0 && v < 4096) ? smallBin[v] : -1;
			if (bi < 0) {
				std::unordered_map<long long,int>::iterator it = binIndex.find(v);
				if (it == binIndex.end()) {
					bi = (int)binList.size();
					binIndex[v] = bi;
					binList.push_back(FFRandBin());
					binList.back().value = v;
					binList.back().numPos = 0;
				} else {
					bi = it->second;
				}
				if (v >= 0 && v < 4096) smallBin[v] = bi;
			}
			if (f1[0] == '1') {
				n++;
				out.append(f0, n0);
				out.append("\t1\n");
				binList[bi].numPos++;
			} else {
				binList[bi].back.push_back(std::make_pair(f0, n0));
			}
		}
	}

	long long num = (long long)maxv - n;
	if (num < n) {
		fprintf(stderr, "Less background than target sequences!!!!!\n");
	}
	if (num < 1) {
		fprintf(stderr, "Something is wrong - not enough background sequences!!!!\n");
		return 0;
	}
	std::vector<std::pair<long long,int> > realBins;
	long long totalPos = 0;
	for (size_t i=0;i<binList.size();i++) {
		if (binList[i].numPos < 1 || binList[i].back.size() < 1) continue;
		totalPos += binList[i].numPos;
		realBins.push_back(std::make_pair(-binList[i].value, (int)i));
	}
	std::sort(realBins.begin(), realBins.end());

	FFPerlRand rnd(x1);
	double leftover = 0.0;
	for (size_t r=0;r<realBins.size();r++) {
		FFRandBin &b = binList[realBins[r].second];
		double num2add = floor((double)b.numPos/(double)totalPos*(double)num+0.49999999);
		num2add += leftover;
		leftover = 0.0;
		long long nids = (long long)b.back.size();
		long long z = 0;
		if (num2add > (double)nids/2.0) {
			// $ids{$_} = rand() for each ID, then a sort by value (all distinct)
			std::vector<std::pair<double,long long> > keyed(nids);
			for (long long i=0;i<nids;i++) keyed[i] = std::make_pair(rnd.next(), i);
			std::stable_sort(keyed.begin(), keyed.end());
			for (long long i=0;i<num2add;i++) {
				if (i >= nids) break;
				out.append(b.back[keyed[i].second].first, b.back[keyed[i].second].second);
				out.append("\t0\n");
				z++;
			}
		} else {
			std::vector<char> mask(nids, 0);
			while (z < num2add) {
				double index = floor(rnd.next()*(double)nids);
				if (index >= (double)nids) continue;
				long long k = (long long)index;
				if (mask[k]) continue;
				mask[k] = 1;
				out.append(b.back[k].first, b.back[k].second);
				out.append("\t0\n");
				z++;
			}
		}
		if (z < num2add) leftover += num2add - z;
	}
	return 0;
}

void randBgProgram(int argc, char** argv) {
	if (argc < 5) printCMDrandbg();
	const char* x1Arg = NULL;
	const char* x2Arg = NULL;
	for (int i=5;i<argc;i++) {
		if (strcmp(argv[i],"-perlrand")==0 && i+2 < argc) {
			x1Arg = argv[++i];
			x2Arg = argv[++i];
		} else {
			exit(2);
		}
	}
	if (x1Arg == NULL) exit(2);
	std::string out;
	if (ffRandBg(argv[2], argv[3], argv[4], x1Arg, x2Arg, out) != 0) exit(2);
	if (out.size() > 0) fwrite(out.data(), 1, out.size(), stdout);
	fflush(stdout);
}

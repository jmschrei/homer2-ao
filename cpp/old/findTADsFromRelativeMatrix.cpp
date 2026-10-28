#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include "SeqTag.h"

#define BUFFER 100000000
#define WHITE_SPACE 7

#define MAX_RELMATRIX_COLUMNS 1000000
#define MAX_RELMATRIX_ROWS 10000000

#define RELMATRIX_TRIAVG 0
#define RELMATRIX_TRIUPAVG 1
#define RELMATRIX_TRIDNAVG 2

#define DEFAULT_WINDOW_SIZE_PIXELS_TAD_REL_MAX 5
#define DEFAULT_PIXELS_FROM_DIAGONAL_TO_START_SCORING 3

#define RELMATRIX_TAD_RATIO_METHOD_AVG 0
#define RELMATRIX_TAD_RATIO_METHOD_MAX 1

#define RELMATRIX_DI_DEFAULT_DIST 1000000
#define RELMATRIX_INS_DEFAULT_DIST 200000
#define RELMATRIX_DI_PSEUDOCOUNT 5

#define TAD_MAXIMUM_OVERLAP 0.20

//void split(char* string, char** cols, int &numCols, char delim);
void findTADs(int argc, char** argv);
void scoreTADs(int argc, char** argv);

class RelMatrix {
public:
	double** m;
	int len; //goes with positions;
	int size; //goes with offsets
	int res;
	int *offsets;
	int *pos;
	int startCoordinate;
	char* chr;
	char* prefix;
	double* mstats;
	double* mThreshold;

	double** tri;
	int filterFlag;

	int minDist;
	int maxDist;
	int minDistIndex;
	int maxDistIndex;
	int diDist;
	int insDist;
	int domainMinSize;
	int domainMinSizeIndex;
	int windowSize;
	int windowSizeIndex;
	int minDistIntIndex;
	int intInnerSize;
	int intOuterSize;
	int innerIntSizeIndex;
	int outerIntSizeIndex;
	double minIntThreshFactor;
	double minIntScore;
	double minIntConstantFactor;
	int ratioMethod;
	double intPseudoFactor;
	double NAvalue;
	double pseudoCount;
	double poissonThreshold;

	RelMatrix(char* file);
	~RelMatrix();
	void init();
	void readRelMatrixFile(char* file);
	void zscoreMatrix();
	double* zscoreDI(double* di);
	void printRelMatrix(FILE* fp, double** matrix);
	void setMinMax(int min, int max, int minTadSize, int pixelWindow,int minIntDist,int innerSize,int outerSize);
	void getMatrixStats();
	double** calculateTriMatrix(int mode);
	double** normalizeMatrix(double** tri, double **upstream, double **downstream);
	//double** findLocalMaxs(double **matrix,int win,double thresh);
	PeakLibrary* findTADsFromNormMatrix(double **matrix,double thresh,int removeOverlap);
	double scoreTriangle(int indexI, int indexJ, int depth);
	void scoreTADs(PeakLibrary* peaks, double** scores, int sindex);
	void scoreInteractions(GenomeInteractionLibrary* peaks, double** scores, int sindex);
	double scoreInteraction(int indexI,int indexJ);

	GenomeInteractionLibrary* findInteractions(int windowSize,double localFoldThreshold);
	double* calculateDI(int length);
	double* calculateInsulationRatio(int insSize,double** triAvg,
							double** triUpAvg, double **triDownAvg);
	void printDIbedGraph(FILE* fp,double* di);
};

void printCMD() {
	fprintf(stderr, "\n\tusage: findTADsFromRelativeMatrix <cmd> [options]\n");
	fprintf(stderr, "\n\tAn internal program used by findTADsAndCPs.pl to calculate topological domains (TAD)\n");
	fprintf(stderr, "\tand contact points (CPs, i.e. interaction anchor points). Use findTADsAndCPs.pl instead.\n");
	fprintf(stderr, "\n\tfindTADsFromRelativeMatrix find <relative matrix> [options]\n");
	fprintf(stderr, "\t\t-prefix <filename> (Filename prefix for output files, default: out)\n");
	fprintf(stderr, "\t\t-minDist <#> (minimum contact distance to score (to avoid diagonal) default: %d x res)\n",
									DEFAULT_PIXELS_FROM_DIAGONAL_TO_START_SCORING);
	fprintf(stderr, "\t\t-diDist <#> (Distance used to calculate directionality index, default 1000000)\n");
	fprintf(stderr, "\t\t-insDist <#> (maximum length of interactions used to insulation ratio, default 200000)\n");
	fprintf(stderr, "\t\tTAD options:\n");
	fprintf(stderr, "\t\t  -minTADscore <#> (minimum inclusion ratio i.e. intra-TAD interactions vs. inter-TAD interaction, def: 1.75)\n");
	fprintf(stderr, "\t\t  -minTADsize <#> (minimum TAD size, default: 100000)\n");
	fprintf(stderr, "\t\t  -sizeTADoptima <#> (window size to look for local optima, default: 5xresolution)\n");
	fprintf(stderr, "\t\t  -noOverlap (remove overlapping TADs if they overlap more than %.2lf - keep best one)\n",TAD_MAXIMUM_OVERLAP);
	fprintf(stderr, "\t\t    -noOverlapStrict (remove overlapping TADs if they overlap at all - keep best one)\n");
	fprintf(stderr, "\t\t  -skipTADs (don't find TADs)\n");
	fprintf(stderr, "\t\tCP options:\n");
	fprintf(stderr, "\t\t  -minCPdist <#> (minimum interaction length for CPs)\n");
	fprintf(stderr, "\t\t  -sizeCP <#> (size of interaction region for CPs, default 5x res)\n");
	fprintf(stderr, "\t\t  -sizeCPlocal <#> (size of interaction bg region CPs, default 15x res)\n");
	fprintf(stderr, "\t\t  -foldCPlocal <#> (local fold threshold for interaction for CPs, def: 1.5)\n");
	fprintf(stderr, "\t\t  -foldCPglobal <#> (fold threshold vs. avg for interaction at given distance for CPs, def: 2.0)\n");
	fprintf(stderr, "\t\t  -poissonCP <#> (maximum poisson p-value for interaction score based on average interactions per dist)\n");
	fprintf(stderr, "\t\t  -minCPreads <#> (minimum number of interaction reads for CP, default: automatic)\n");
	fprintf(stderr, "\t\t  -pseudoCPreads <#> (pseudo count for interaction threshold for CPs, def: 5.0)\n");
	fprintf(stderr, "\t\t  -skipCPs (don't find contact points/anchors)\n");
	//fprintf(stderr, "\t\t-intFthreshold (fold threshold vs. avg for interaction for CPs, def: 2.0)\n");
	//fprintf(stderr, "\t\t-intCthreshold (pseudo count for interaction threshold for CPs, def: 5.0)\n");
	//fprintf(stderr, "\t\t-intLthreshold (local fold threshold for interaction for CPs, def: 1.5)\n");
	fprintf(stderr, "\n\tfindTADsFromRelativeMatrix score -tad <2D bed> -cp <2D bed> -matrix <rel matrix1> [rel matrix2] ...\n");
	fprintf(stderr, "\t\t-matrix <relative file1> [relative file2] ... (relative matricies to use for scoring)\n");
	fprintf(stderr, "\t\t-prefix <filename> (Filename prefix for output files, default: out)\n");
	fprintf(stderr, "\t\t-tad <TAD 2D BED file> (score TAD on inclusion ratio)\n");
	fprintf(stderr, "\t\t-cp <contact point 2D BED file> (score interactions found at contact points)\n");
	fprintf(stderr, "\t\t-minDist <#> (minimum contact distance to score (to avoid diagonal) default: %d x res)\n",
									DEFAULT_PIXELS_FROM_DIAGONAL_TO_START_SCORING);
	fprintf(stderr, "\t\t\n");
	exit(0);
}
int main (int argc, char** argv) {

	if (argc < 2) {
		//fprintf(stderr, "Problem asshole!\n");
		printCMD();
		exit(0);
	}
	if (strcmp("find",argv[1])==0) {
		findTADs(argc,argv);
	} else if (strcmp("score",argv[1])==0) {
		scoreTADs(argc,argv);
	} else {
		printCMD();
	}
}
void scoreTADs(int argc, char** argv) {
	int minDist = -1;
	int maxDist = INT_MAX;
	int windowSize = -1;
	int minIntDist = -1;
	//int minTADsize = -1;
	int minTADsize = 100000;
 	int intInnerSize = -1;
 	int intOuterSize = -1;
	const char* prefixConst = "out";
	char* prefix = (char*)prefixConst;
	char* int2dBEDfile = NULL;
	char* tad2dBEDfile = NULL;

	char** relMatrixFiles = new char*[argc];
	int numRelMatrixFiles = 0;

	for (int i=2;i<argc;i++) {
		if (strcmp("-matrix",argv[i])==0) {
			i++;
			for (;i<argc;i++) {
				if (argv[i][0] == '-') {
					i--;
					break;
				}
				relMatrixFiles[numRelMatrixFiles++] = argv[i];
			}
		} else if (strcmp(argv[i],"-prefix")==0) {
			prefix = argv[++i];
		} else if (strcmp(argv[i],"-cp")==0) {
			int2dBEDfile = argv[++i];
		} else if (strcmp(argv[i],"-tad")==0) {
			tad2dBEDfile = argv[++i];
		} else if (strcmp(argv[i],"-minDist")==0) {
			sscanf(argv[++i],"%d",&minDist);
		} else if (strcmp(argv[i],"-sizeCP")==0) {
			sscanf(argv[++i],"%d",&intInnerSize);
		} else if (strcmp(argv[i],"-sizeCPlocal")==0) {
			sscanf(argv[++i],"%d",&intOuterSize);
		} else {
			fprintf(stderr, "!!! CMD line Error: %s not recognized.\n", argv[i]);
			printCMD();
		}
	}
	char* filename = new char[strlen(prefix)+1000];

	PeakLibrary* tads = NULL;
	double **scores = NULL;
	if (tad2dBEDfile != NULL) {
		//tads = new PeakLibrary(tad2dBEDfile,PEAK_READ_MODE_NORMAL);
		tads = new PeakLibrary(tad2dBEDfile,PEAK_READ_MODE_2DBED);
		tads->setDefaultPeakOrder();
		fprintf(stderr, "\tFound %d TAD regions in file %s\n", tads->numPeaks, tad2dBEDfile);
		scores = new double*[tads->numPeaks];
		for (int i=0;i<tads->numPeaks;i++) {
			scores[i] = new double[numRelMatrixFiles];
			for (int j=0;j<numRelMatrixFiles;j++) {
				scores[i][j]=0.0;
			}
		}
	}
	
	double **intScores = NULL;
	GenomeInteractionLibrary* ints = NULL;
	if (int2dBEDfile != NULL) {
		ints = new GenomeInteractionLibrary();
		ints->read2dBED(int2dBEDfile);
		//fprintf(stderr, "\tScoring %d interactions\n", ints->numInteractions);
		intScores = new double*[ints->numInteractions];
		for (int i=0;i<ints->numInteractions;i++) {
			intScores[i] = new double[numRelMatrixFiles];
			for (int j=0;j<numRelMatrixFiles;j++) {
				intScores[i][j]=0.0;
			}
		}
	}


	for (int i=0;i<numRelMatrixFiles;i++) {
		RelMatrix* matrix = new RelMatrix(relMatrixFiles[i]);
		matrix->setMinMax(minDist,maxDist,minTADsize,windowSize,minIntDist,intInnerSize,intOuterSize);
		if (tads != NULL) matrix->scoreTADs(tads,scores,i);
		if (ints != NULL) matrix->scoreInteractions(ints,intScores,i);
	}

	if (tads != NULL) {
		sprintf(filename,"%s.tad.scores.txt",prefix);
		FILE* fp = fopen(filename,"w");

		fprintf(fp,"#TAD name\tchr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tinclusionRatio\tsize");
		for (int i=0;i<numRelMatrixFiles;i++) {
			fprintf(fp,"\t%s", relMatrixFiles[i]);
		}
		fprintf(fp,"\n");
		for (int i=0;i<tads->numPeaks;i++) {
			Peak* p = tads->peakOrder[i];
			//char dir = '+';
			//if (p->strand == 1 || p->strand == '-') dir = '-';
			fprintf(fp, "%s\t%s\t%d\t%d\t%s\t%d\t%d\t255,255,0",p->name,p->chr,p->start,p->end,p->chr,p->start,p->end);
			//fprintf(fp, "\t%f\t%d", p->v,p->end - p->start);
			fprintf(fp, "\t%f\t%d", p->v,p->end - p->start);
			for (int j=0;j<numRelMatrixFiles;j++) {
				fprintf(fp, "\t%lf",scores[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

		for (int i=0;i<tads->numPeaks;i++) {
			delete [](scores[i]);
		}
		delete []scores;
		delete tads;
	}

	if (ints != NULL) {
		sprintf(filename,"%s.cp.scores.txt",prefix);
		FILE* fp = fopen(filename,"w");

		fprintf(fp,"#interaction name\tchr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tscore\tlocalRatio");
		for (int i=0;i<numRelMatrixFiles;i++) {
			fprintf(fp,"\t%s", relMatrixFiles[i]);
		}
		fprintf(fp,"\n");
		for (int i=0;i<ints->numInteractions;i++) {
			GenomeInteraction* p = ints->interactions[i];
			fprintf(fp, "%s\t%s\t%d\t%d\t%s\t%d\t%d\t0,0,255",p->name,p->chr1,p->start1,p->end1,p->chr2,p->start2,p->end2);
			fprintf(fp, "\t%lf\t%lf",p->interactions,p->expected);
			for (int j=0;j<numRelMatrixFiles;j++) {
				fprintf(fp, "\t%lf",intScores[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

		for (int i=0;i<ints->numInteractions;i++) {
			delete [](intScores[i]);
		}
		delete []intScores;
		delete ints;
	}
	delete []filename;
}	
void findTADs(int argc, char** argv) {

	int minDist = -1;
	int maxDist = INT_MAX;
	int windowSize = -1;
	int minTADsize = 100000;
	int removeOverlap = 0;
	char* file = argv[2];
	double threshold = 1.75;
	double localFoldThreshold = 1.5;
	double minIntThreshFactor = 2.0;
	double minIntScore = 0.0;
	double minIntConstantFactor = 5.0;
	double poissonCP = 2.0;
	const char* defaultPrefix = "out";
	char* prefix = NULL;
	int filterFlag = 1;
	int skipTADs = 0;
	int minIntDist = -1;
	int skipInts = 0;
	int diSize = -1;
	int insSize = -1;
	int intInnerSize = -1;
	int intOuterSize = -1;
	for (int i=3;i<argc;i++) {
		if (strcmp(argv[i],"-prefix")==0) {
			prefix = argv[++i];
		} else if (strcmp(argv[i],"-nofilter")==0) {
			filterFlag = 0;
		} else if (strcmp(argv[i],"-minTADsize")==0) {
			sscanf(argv[++i],"%d",&minTADsize);
		} else if (strcmp(argv[i],"-minCPdist")==0) {
			sscanf(argv[++i],"%d",&minIntDist);
		} else if (strcmp(argv[i],"-minCPreads")==0) {
			sscanf(argv[++i],"%lf",&minIntScore);
		} else if (strcmp(argv[i],"-sizeTADoptima")==0) {
			sscanf(argv[++i],"%d",&windowSize);
		} else if (strcmp(argv[i],"-skipTADs")==0) {
			skipTADs = 1;
		} else if (strcmp(argv[i],"-skipCPs")==0) {
			skipInts = 1;
		} else if (strcmp(argv[i],"-noOverlap")==0) {
			removeOverlap = 1;
		} else if (strcmp(argv[i],"-noOverlapStrict")==0) {
			removeOverlap = 2;
		} else if (strcmp(argv[i],"-diDist")==0) {
			sscanf(argv[++i],"%d",&diSize);
		} else if (strcmp(argv[i],"-insDist")==0) {
			sscanf(argv[++i],"%d",&insSize);
		} else if (strcmp(argv[i],"-minTADscore")==0) {
			sscanf(argv[++i],"%lf",&threshold);
		} else if (strcmp(argv[i],"-minDist")==0) {
			sscanf(argv[++i],"%d",&minDist);
		} else if (strcmp(argv[i],"-foldCPlocal")==0) {
			sscanf(argv[++i],"%lf",&localFoldThreshold);
			//} else if (strcmp(argv[i],"-intCthreshold")==0) {
		} else if (strcmp(argv[i],"-pseudoCPreads")==0) {
			sscanf(argv[++i],"%lf",&minIntConstantFactor);
		} else if (strcmp(argv[i],"-poissonCP")==0) {
			sscanf(argv[++i],"%lf",&poissonCP);
		} else if (strcmp(argv[i],"-foldCPglobal")==0) {
			sscanf(argv[++i],"%lf",&minIntThreshFactor);
		} else if (strcmp(argv[i],"-sizeCP")==0) {
			sscanf(argv[++i],"%d",&intInnerSize);
		} else if (strcmp(argv[i],"-sizeCPlocal")==0) {
			sscanf(argv[++i],"%d",&intOuterSize);
		} else {
			fprintf(stderr, "!!! CMD line Error: %s not recognized.\n", argv[i]);
			printCMD();
		}
	}
	fprintf(stderr, "\tThreshold for TAD calling (-minTADscore)   : %.3lf\n", threshold);
	if (prefix == NULL) {
		prefix = (char*)defaultPrefix;
	}
	char* filename = new char[strlen(prefix)+1000];

	RelMatrix* matrix = new RelMatrix(file);
	matrix->setMinMax(minDist,maxDist,minTADsize,windowSize,minIntDist,intInnerSize,intOuterSize);

	PeakLibrary* tads = NULL;

	if (skipTADs == 0) {
		//matrix->zscoreMatrix();
		//matrix->printRelMatrix(stdout,matrix->m);
		fprintf(stderr, "\tFinished reading file\n");
		double** triAvg = matrix->calculateTriMatrix(RELMATRIX_TRIAVG);
		fprintf(stderr, "\tFinished calculating triangle density\n");
		double** triUpAvg = matrix->calculateTriMatrix(RELMATRIX_TRIUPAVG);
		fprintf(stderr, "\tFinished calculating upper control triangle density\n");
		double** triDnAvg = matrix->calculateTriMatrix(RELMATRIX_TRIDNAVG);
		fprintf(stderr, "\tFinished calculating lower control triangle density\n");
		double** norm = matrix->normalizeMatrix(triAvg, triUpAvg, triDnAvg);
	
		//matrix->printRelMatrix(stdout,norm);
		//exit(0);
	
		matrix->filterFlag = filterFlag;
		tads = matrix->findTADsFromNormMatrix(norm,threshold,removeOverlap);

		sprintf(filename, "%s.tad.2D.bed",prefix);
		FILE* fp = fopen(filename,"w");
		tads->printBED2D(fp);
		fclose(fp);

		double* di = matrix->calculateDI(diSize);
		di = matrix->zscoreDI(di);
		sprintf(filename, "%s.DI.bedGraph",prefix);
		fp = fopen(filename,"w");
		fprintf(fp, "track name=\"%s DI\" type=bedGraph\n",prefix);
		matrix->printDIbedGraph(fp,di);
		fclose(fp);

		double* ins = matrix->calculateInsulationRatio(insSize,triAvg,triUpAvg,triDnAvg);
		sprintf(filename, "%s.Insulation.bedGraph",prefix);
		fp = fopen(filename,"w");
		fprintf(fp, "track name=\"%s InsulationRatio\" type=bedGraph\n",prefix);
		matrix->printDIbedGraph(fp,ins);
		fclose(fp);

		//tads->print(stdout);
	} else {
		fprintf(stderr, "\tSkipping TAD calling...\n");
	}
	if (skipInts == 0) {
		//findInteractions
		matrix->minIntThreshFactor = minIntThreshFactor;
		matrix->minIntScore = minIntScore;
		matrix->minIntConstantFactor = minIntConstantFactor;
		matrix->poissonThreshold = poissonCP;
		matrix->getMatrixStats();
		if (matrix->mThreshold != NULL) {
			fprintf(stderr, "\tMin int reads at maxDist: %lf\n", matrix->mThreshold[matrix->size-1]);
		}
		fprintf(stderr, "\tFinding contact points/anchors:\n");
		fprintf(stderr, "\tFold over distance-average for contact point calling (-foldCPgloball): %.3lf\n", minIntThreshFactor);
		fprintf(stderr, "\tConstant Threshold for contact point calling (-pseudoCPreads): %.3lf\n", minIntConstantFactor);
		fprintf(stderr, "\tLocal Fold Threshold for contact point calling (-foldCPlocal): %.3lf\n", localFoldThreshold);
		fprintf(stderr, "\tMaximum poisson p-values threshold for minimum interaction score (-poissonCP): %.3le\n", poissonCP);
		fprintf(stderr, "\tMinimum Interaction score (-minCPreads): %.3lf\n", minIntScore);
		GenomeInteractionLibrary* interactions = matrix->findInteractions(windowSize,localFoldThreshold);
		sprintf(filename, "%s.cp.2D.bed",prefix);
		FILE* fp = fopen(filename,"w");
		interactions->print2dBED(fp);
		fclose(fp);
	} else {
		fprintf(stderr, "\tSkipping contact point calling...\n");
	}
}

RelMatrix::RelMatrix(char* filename) {
	init();
	readRelMatrixFile(filename);
}
RelMatrix::~RelMatrix() {
	if (m != NULL) {
		for (int i=0;i<len;i++) {
			if (m[i] != NULL) delete [](m[i]);
		}
		delete []m;
	}
	if (mstats != NULL) {
		delete []mstats;
	}
	if (mThreshold != NULL) {
		delete []mThreshold;
	}
}
void RelMatrix::init() {
	m = NULL;
	len = 0;
	size = 0;
	res = 0;
	chr = NULL;
	pos = NULL;
	offsets = NULL;
	NAvalue = 1.0;
	startCoordinate=0;
	intPseudoFactor = 1.0;
	minIntThreshFactor = 2.0;
	minIntConstantFactor = 20.0;
	innerIntSizeIndex = 5;
	outerIntSizeIndex = 15;
	pseudoCount = 0.0;
	mstats = NULL;
	poissonThreshold = 1.0;
	mThreshold = NULL;
	minDistIntIndex = 6;
	diDist = RELMATRIX_DI_DEFAULT_DIST;
	insDist = RELMATRIX_INS_DEFAULT_DIST;
	ratioMethod = RELMATRIX_TAD_RATIO_METHOD_AVG;

}
double RelMatrix::scoreInteraction(int i,int j) {

	double** matrix = m;
	
	//sig interactions must be at least greater than the average for that distance

	//double intScore = matrix[i][j];
	double innerSum=0.0;
	double outerSum=0.0;
	int innerN = 0;
	int outerN = 0;
	for (int k=-1*outerIntSizeIndex;k<outerIntSizeIndex;k++) {
		int x = i+k;
		if (x < 0) continue;
		if (x >= len) continue;
		for (int m= -1*outerIntSizeIndex;m<outerIntSizeIndex;m++) {
			int y = j+m;
			//if (y < minDistIntIndex) continue;
			if (y < 0) continue;
			if (y > maxDistIndex) continue;
			if ((k < innerIntSizeIndex && k >= -1*innerIntSizeIndex) 
						&& (m < innerIntSizeIndex && m >= -1*innerIntSizeIndex)) {
				innerSum+=matrix[x][y];
				innerN++;
			} else {
				outerSum+=matrix[x][y];
				outerN++;
			}
		}
	}

	if (innerN > 0) innerSum /= innerN;
	if (outerN > 0) outerSum /= outerN;
	//double pseudoInt = mstats[j]*intPseudoFactor;
	//double fscore = (innerSum+pseudoInt)/(outerSum+pseudoInt);
	//return fscore;
	return innerSum;
	//return intScore;

}
GenomeInteractionLibrary* RelMatrix::findInteractions(int windowSize,double localFoldThreshold) {



	double** matrix = m;
	double** mm = NULL;
	int smoothingFlag = 1;
	if (smoothingFlag) {
		mm = new double*[len];
		for (int i=0;i<len;i++) {
			mm[i] = new double[size];
			for (int j=0;j<size;j++) {
				mm[i][j] = 0;
			}
		}
		for (int i=0;i<len;i++) {
			for (int j=0;j<size;j++) {
				if (j < minDistIntIndex) continue;
				double intScore = m[i][j];
				//if (intScore < mstats[j]*minThreshFactor) continue;
				//if (intScore < mstats[j]*minIntThreshFactor+minIntConstantFactor) continue;
				if (intScore < mstats[j]*minIntThreshFactor) continue;
				if (intScore < mThreshold[j]) continue;
				double innerSum=0.0;
				int innerN = 0;
				for (int k=-1*innerIntSizeIndex;k<innerIntSizeIndex;k++) {
					int x = i+k;
					if (x < 0) continue;
					if (x >= len) continue;
					for (int l= -1*innerIntSizeIndex;l<innerIntSizeIndex;l++) {
						int y = j+l;
						//if (y < minDistIntIndex) continue;
						if (y < 0) continue;
						if (y > maxDistIndex) continue;
						if ((k < innerIntSizeIndex && k >= -1*innerIntSizeIndex) 
										&& (l < innerIntSizeIndex && l >= -1*innerIntSizeIndex)) {
							innerSum+=m[x][y];
							innerN++;
						}
					}
				}
				if (innerN > 0) innerSum /= innerN;
				mm[i][j] = innerSum;
				//fprintf(stderr, "%lf\t%d\t%d\t%d\n", mm[i][j],i,j,innerN);
			}
		}
		matrix = mm;
	}
	
	//sig interactions must be at least greater than the average for that distance
	//double constantFactor = 25.0;
	//int minDistIntIndex = 6;


	GenomeInteractionLibrary *ints = new GenomeInteractionLibrary();

	//int z= 0;
	for (int i=0;i<len;i++) {
		for (int j=0;j<size;j++) {
			if (j < minDistIntIndex) continue;
			double intScore = matrix[i][j];
			if (intScore <= 1e-20) continue;
			if (intScore < minIntScore) continue;
		//	if (smoothingFlag == 0 && intScore < mstats[j]*minIntThreshFactor+minIntConstantFactor) continue;
			if (intScore < mstats[j]*minIntThreshFactor) continue;
			if (intScore < mThreshold[j]) continue;

			int bad = 0;
			double innerSum=0.0;
			double outerSum=0.0;
			int innerN = 0;
			int outerN = 0;
			for (int k=-1*outerIntSizeIndex;k<outerIntSizeIndex;k++) {
				int x = i+k;
				if (x < 0) continue;
				if (x >= len) continue;
				for (int l= -1*outerIntSizeIndex;l<outerIntSizeIndex;l++) {
					int y = j+l;
					//if (y < minDistIntIndex) continue;
					if (y < 0) continue;
					if (y > maxDistIndex) continue;
					if (x!=i || y!=j) {
						if (matrix[x][y] > matrix[i][j]) {
							bad =1;
							break;
						}
						if (matrix[x][y] == matrix[i][j] && x+y >= i+j) {
							//fprintf(stderr, "%lf\t%lf\t%d\t%d\t%d\t%d\n", matrix[x][y],matrix[i][j],x,y,i,j);
							bad=1;
							break;
						}
					}
					if ((k < innerIntSizeIndex && k >= -1*innerIntSizeIndex) 
									&& (l < innerIntSizeIndex && l >= -1*innerIntSizeIndex)) {
						innerSum+=m[x][y];
						innerN++;
					} else {
						outerSum+=m[x][y];
						outerN++;
					}
				}
				if (bad) break;
			}
			if (bad) continue;

			if (innerN > 0) innerSum /= innerN;
			if (outerN > 0) outerSum /= outerN;

			double pseudoInt = mstats[j]*intPseudoFactor;
			//fprintf(stderr, "psuedo=%lf\t%lf\t%lf\t%lf\n", pseudoInt, innerSum, outerSum,intScore);
			double fscore = (innerSum+pseudoInt)/(outerSum+pseudoInt);

			if (bad==0 && fscore > localFoldThreshold) {
				int s = pos[i];
				int e = s+offsets[j];
				//float v = (float)matrix[i][j];
				char* str = new char[1000];
				//sprintf(str, "%s:%d-%d", chr, s,e);
				//fprintf(stderr, "\t%d\t%s\t%s\n", z++,chr,str);
				//
				GenomeInteraction* interaction = new GenomeInteraction(str,chr,s,s+res,0.0,i,NULL,
										chr,e,e+res,0.0,j,NULL,innerSum,fscore,0.0,0.0,1);
				ints->addInteraction(interaction);

				//candidatesX[numCandidates] = i;
				//candidatesY[numCandidates] = j;
				//numCandidates++;
			}
		}
	}

	//fprintf(stderr, "\tFound %d total interactions\n", ints->numInteractions);
	ints->optimizeInteractionArray();
	ints->sortInteractionArray();
	fprintf(stderr, "\tFound %d total interactions\n", ints->numInteractions);
	return ints;
}
void RelMatrix::scoreTADs(PeakLibrary* peaks, double** scores, int sindex) {

	for (int i=0;i<peaks->numPeaks;i++) {
		Peak* tad = peaks->peakOrder[i];
		if (strcmp(chr,tad->chr) != 0) continue;
		int indexI = (tad->start - startCoordinate)/res;
		int indexJ = (tad->end - tad->start)/res;
		double s = scoreTriangle(indexI,indexJ,indexJ);
		scores[i][sindex] = s;
	}
}
void RelMatrix::scoreInteractions(GenomeInteractionLibrary* ints,double** intScores,int sindex) {
	for (int i=0;i<ints->numInteractions;i++) {
		GenomeInteraction* gi = ints->interactions[i];
		if (strcmp(chr,gi->chr1) != 0) continue;
		if (strcmp(chr,gi->chr2) != 0) continue;
		int indexI = (gi->start1 - startCoordinate)/res;
		int indexJ = (gi->start2 - gi->start1)/res;
		if (gi->start2 < gi->start1) {
			indexI = (gi->start2 - startCoordinate)/res;
			indexJ = (gi->start1 - gi->start2)/res;
		}
		double s = 0.0;
		if (indexJ < size && indexJ >= 0 && indexI < len && indexI > 0) {
			s = scoreInteraction(indexI,indexJ);
		}
		intScores[i][sindex] = s;
	}
}
PeakLibrary* RelMatrix::findTADsFromNormMatrix(double **matrix,double thresh,int removeOverlap) {

	//int candidatesX = new int[len*2];
	//int candidatesY = new int[len*2];
	//int numCandidates = 0;
	PeakLibrary* putativeTADs = new PeakLibrary(100000);

	//int z= 0;
	for (int i=0;i<len;i++) {
		for (int j=0;j<size;j++) {
			if (matrix[i][j] < thresh) continue;
			int bad = 0;
			for (int k=-windowSizeIndex;k<windowSizeIndex;k++) {
				int x = i+k;
				if (x < 0) continue;
				if (x >= len) continue;
				for (int m=-windowSizeIndex;m<windowSizeIndex;m++) {
					int y = j+m;
					if (y < minDistIndex) continue;
					if (y > maxDistIndex) continue;
					if (matrix[x][y] > matrix[i][j]) {
						bad =1;
						break;
					}
				}
				if (bad) break;
			}
			if (bad==0 && matrix[i][j] >= thresh && j>=domainMinSizeIndex && (i+j < len)) {

				int s = pos[i];
				int e = s+offsets[j]+offsets[1]-1;
				float v = (float)matrix[i][j];
				char* str = new char[1000];
				sprintf(str, "%s:%d-%d", chr, s,e);
				//fprintf(stderr, "\t%d\t%s\t%s\n", z++,chr,str);
				Peak* p = putativeTADs->addPeak(str,chr,s,e,(s+e)/2,0,v,v,NULL,0,1);
				p->tagStart = i;
				p->tagEnd = j;
				//candidatesX[numCandidates] = i;
				//candidatesY[numCandidates] = j;
				//numCandidates++;
			}
		}
	}

	//fprintf(stderr, "Looking at start coordinats\n");
	putativeTADs->setDefaultPeakOrder(PEAK_SORT_BY_CHR_START_COORDINATE);
	//putativeTADs->printCurrentPeakOrder(stdout);
	//
	
	Peak** peakset = new Peak*[1000];
	int numPeakSet = 0;
	for (int i=0;i<putativeTADs->numPeaks;i++) {
		numPeakSet = 0;
		Peak* p =  putativeTADs->peakOrder[i];
		for (int j=i-1;j>=0;j--) {
			if (p->tagStart - putativeTADs->peakOrder[j]->tagStart <= windowSizeIndex) {
				peakset[numPeakSet++] = putativeTADs->peakOrder[j];
			} else {
				break;
			}
		}
		for (int j=i+1;j<putativeTADs->numPeaks;j++) {
			if (putativeTADs->peakOrder[j]->tagStart - p->tagStart <= windowSizeIndex) {
				peakset[numPeakSet++] = putativeTADs->peakOrder[j];
			} else {
				break;
			}
		}

		//double score = (double)p->v;
		//double rescore = scoreTriangle(p->tagStart,p->tagEnd,p->tagEnd);
		//fprintf(stderr, "\t%s\t%d\t%d\t%lf\t%lf\n", p->name, p->tagStart,p->tagEnd,score,rescore);
	
		//int psize = putativeTADs->peakOrder[i]->end - putativeTADs->peakOrder[i]->start;
		for (int j=0;j<numPeakSet;j++) {
			Peak* pp = peakset[j];

			int diff = p->tagEnd - pp->tagEnd;
			if (diff < 0) {
				//other pp/j domain is larger (superdomain)
				double subscore = scoreTriangle(pp->tagStart,pp->tagEnd,-1*diff);
				//fprintf(stderr, "\t\t+sup: %s\t%lf\tsupDiff=%lf (%d)\n", pp->name,(double)pp->v,subscore,diff);
				if (filterFlag) {
					//if (pp->v > p->v) p->priority = 0;
					//if (subscore > p->v) p->priority = 0;
					if (p->v > pp->v) {
						 if (subscore < thresh) pp->priority = 0;
					} else {
						if (subscore > p->v) p->priority = 0;
					}
				}
			} else {
				// pp/j domain is likely a subdomain
				double subscore = scoreTriangle(p->tagStart,p->tagEnd,diff);
				//fprintf(stderr, "\t\t-sub: %s\t%lf\tsubDiff=%lf (%d)\n", pp->name,(double)pp->v,subscore,diff);
				if (filterFlag) {
					//if (p->v > pp->v) pp->priority = 0;
					//if (subscore > pp->v) pp->priority = 0;
					if (pp->v > p->v) {
						if (subscore < thresh) p->priority = 0;
					} else {
						if (subscore > pp->v) pp->priority = 0;
					}
				}
			}
		}
	}

	//fprintf(stderr, "Looking at end coordinats\n");
	putativeTADs->setDefaultPeakOrder(PEAK_SORT_BY_CHR_END_POS_COORDINATE);
	//putativeTADs->setDefaultPeakOrder(PEAK_SORT_BY_CHR_START_COORDINATE);
	for (int i=0;i<putativeTADs->numPeaks;i++) {
		numPeakSet = 0;
		Peak* p =  putativeTADs->peakOrder[i];
		//p->print(stderr);
		//continue;
		for (int j=i-1;j>=0;j--) {
			Peak* pp = putativeTADs->peakOrder[j];
			if ((p->tagStart + p->tagEnd) - 
					(pp->tagStart+pp->tagEnd) <= windowSizeIndex) {
				peakset[numPeakSet++] = pp;
			} else {
				break;
			}
		}
		for (int j=i+1;j<putativeTADs->numPeaks;j++) {
			Peak* pp = putativeTADs->peakOrder[j];
			if ((pp->tagStart+pp->tagEnd) - (p->tagStart+p->tagEnd) <= windowSizeIndex) {
				peakset[numPeakSet++] = pp;
			} else {
				break;
			}
		}

		//double score = (double)p->v;
		//double rescore = scoreTriangle(p->tagStart,p->tagEnd,p->tagEnd);
		//fprintf(stderr, "\t%s\t%d\t%d\t%lf\t%lf\n", p->name, p->tagStart,p->tagEnd,score,rescore);

		//int psize = putativeTADs->peakOrder[i]->end - putativeTADs->peakOrder[i]->start;
		for (int j=0;j<numPeakSet;j++) {
			Peak* pp = peakset[j];

			//int diff = p->tagEnd - pp->tagEnd;
			int diff = p->tagStart - pp->tagStart;
			if (diff < 0) {
				//other j domain is larger (superdomain)
				double subscore = scoreTriangle(pp->tagStart,pp->tagEnd,-1*diff);
				if (filterFlag) {
					//if (pp->v > p->v) p->priority = 0;
					//if (subscore > p->v) p->priority = 0;
					if (p->v > pp->v) {
						 if (subscore < thresh) pp->priority = 0;
					} else {
						if (subscore > p->v) p->priority = 0;
					}
				}
				//fprintf(stderr, "\t\t+sup: %s\t%lf\tsupDiff=%lf (%d %d)\n", pp->name,(double)pp->v,subscore,p->priority,pp->priority);
			} else {
				// j domain is likely a subdomain
				double subscore = scoreTriangle(p->tagStart,p->tagEnd,diff);
				if (filterFlag) {
					//if (p->v > pp->v) pp->priority = 0;
					//if (subscore > pp->v) pp->priority = 0;
					//if (subscore < thresh) p->priority = 0;
					if (pp->v > p->v) {
						if (subscore < thresh) p->priority = 0;
					} else {
						if (subscore > pp->v) pp->priority = 0;
					}
				}
				//fprintf(stderr, "\t\t-sub: %s\t%lf\tsubDiff=%lf (%d %d)\n", pp->name,(double)pp->v,subscore,p->priority,pp->priority);
			}
		}
	}
	delete []peakset;


	PeakLibrary* TADs = new PeakLibrary(putativeTADs->numPeaks);
	for (int i=0;i<putativeTADs->numPeaks;i++) {
		Peak* p = putativeTADs->peakOrder[i];
		if (p->priority) {
			TADs->addPeak(p);
		}
	}
	fprintf(stderr, "\tFound %d total TADs\n", TADs->numPeaks);
	delete putativeTADs;

	if (removeOverlap) {
		fprintf(stderr, "\tRemoving Overlapping TADs\n");
		TADs->setDefaultPeakOrder(PEAK_SORT_BY_CHR_START_COORDINATE);
		for (int i=0;i<TADs->numPeaks;i++) {
			TADs->peakOrder[i]->priority=1;
		}
		for (int i=0;i<TADs->numPeaks;i++) {
			Peak* p1 = TADs->peakOrder[i];
			int p1size = TADs->peakOrder[i]->end - TADs->peakOrder[i]->start;
			for (int j=i+1;j<TADs->numPeaks;j++){ 
				Peak* p2 = TADs->peakOrder[j];
				int p2size = TADs->peakOrder[j]->end - TADs->peakOrder[j]->start;
				if (p1->end <= p2->start) break;
				int overlapEnd = TADs->peakOrder[i]->end;
				if (TADs->peakOrder[i]->end > TADs->peakOrder[j]->end) overlapEnd = TADs->peakOrder[j]->end;
				int overlapSize = overlapEnd - TADs->peakOrder[j]->start;
				double maxRatio = ((double)overlapSize)/((double)p1size);
				if (p2size < p1size) {
					maxRatio = ((double)overlapSize)/((double)p2size);
				}
				if (p1->v < p2->v) {
					if (removeOverlap == 2 || maxRatio > TAD_MAXIMUM_OVERLAP) {
						p1->priority=0;
						continue;
					}
				}
				if (p1->v >= p2->v) {
					if (removeOverlap == 2 || maxRatio > TAD_MAXIMUM_OVERLAP) {
						p2->priority=0;
						continue;
					}
				}
			}
		}
		PeakLibrary* nTADs = new PeakLibrary(TADs->numPeaks);
		for (int i=0;i<TADs->numPeaks;i++) {
			Peak* p = TADs->peakOrder[i];
			if (p->priority) {
				nTADs->addPeak(p);
			}
		}
		fprintf(stderr, "\tFound %d total TADs after removing overlapping TADs\n", nTADs->numPeaks);
		delete TADs;
		TADs = nTADs;
	}

	return TADs;
}
double** RelMatrix::normalizeMatrix(double** tri, double **upstream, double **downstream) {

	/*pseudoCount = 0.0;
	for (int i=0;i<len;i++) {
		pseudoCount += tri[i][minDistIndex];
	}
	pseudoCount /= 2*len;
	fprintf(stderr, "\tPseudoCount avg = %lf\n", pseudoCount);
	*/

	double** norm = new double*[len];
	for (int i=0;i<len;i++) {
		norm[i] = new double[size];
		for (int j=0;j<size;j++) {
			double avg = (upstream[i][j] + downstream[i][j])/2.0;
			if (ratioMethod == RELMATRIX_TAD_RATIO_METHOD_MAX) {
				avg = upstream[i][j];
				if (downstream[i][j] > upstream[i][j]) {
					avg = downstream[i][j];
				}
			}
			if (avg > 0) {
				//double pc = ((pseudoCount*minDistIndex)/j);
				double pc = mstats[minDistIndex];
			//	double pc = mstats[j];
				//fprintf(stderr, "pc=%lf\tpc1=%lf\n", pc, pc1);
				//double pc = mstats[minDistIndex];
				//fprintf(stderr, "%lf\t%lf\t%lf\t%d\n", tri[i][j],avg,pc,j);
				norm[i][j] = (tri[i][j]+pc)/(avg+pc);
			} else {
				norm[i][j]=0.0;
			}
		}
	}
	return norm;
}

double* RelMatrix::zscoreDI(double* di) {

	double avg = 0.0;
	double std = 0.0;
	double N = 0.0;
	for (int i=0;i<len;i++) {
		N+=1.0;
		avg+=di[i];
	}
	if (N > 0.0) avg /= N;
	for (int i=0;i<len;i++) {
		std += (avg-di[i])*(avg-di[i]);
	}
	if (N > 0.0) std /= N;
	std = sqrt(std);
	if (std == 0.0) std = 1.0;

	for (int i=0;i<len;i++) {
		di[i] = (di[i]-avg)/std;
	}
	return di;
}

void RelMatrix::zscoreMatrix() {

	double* avg = new double[size];
	double* std = new double[size];
	double* N = new double[size];
	for (int i=0;i<size;i++){ 
		avg[i] = 0.0;
		std[i] = 0.0;
		N[i] = 0.0;
	}
	for (int i=0;i<len;i++) {
		for (int j=0;j<size;j++) {
			if (m[i][j] < 1.0 + 1e-100 && m[i][j] > 1.0 - 1e-100) {
			} else {
				N[j]+=1.0;
				avg[j]+=m[i][j];
			}
		}
	}
	for (int j=0;j<size;j++) {
		if (N[j] > 0) {
			avg[j] /= N[j];
		}
	}
	for (int i=0;i<len;i++) {
		for (int j=0;j<size;j++) {
			if (m[i][j] < 1.0 + 1e-100 && m[i][j] > 1.0 - 1e-100) {
			} else {
				std[j] += (avg[j]-m[i][j])*(avg[j]-m[i][j]);
			}
		}
	}
	for (int j=0;j<size;j++) {
		if (N[j] > 0) {
			std[j] /= N[j];
		}
		std[j] = sqrt(std[j]);
		//fprintf(stderr, "\t\t%d\t%d\t%lf\t%lf\n", j, offsets[j], avg[j], std[j]);
	}
	for (int i=0;i<len;i++) {
		for (int j=0;j<size;j++) {
			if (N[j]>0 && std[j] > 0) {
				m[i][j] = (m[i][j]-avg[j])/std[j];
			}
		}
	}
	delete []avg;
	delete []std;
	delete []N;
}
			

double RelMatrix::scoreTriangle(int indexI, int indexJ, int depth) {

	if (indexI+indexJ >= len) return 0.0;
	if (indexJ >= size) return 0.0;

	int startIndex = minDistIndex;
	if (indexJ-depth > startIndex) {
		startIndex = indexJ-depth;
	}
	double triAvg = 0.0;
	double triN = 0.0;
	for (int j=startIndex;j<=indexJ;j++) {
		double colSum = 0.0;
		double colN = 0.0;
		int ilen = indexJ-j;
		for (int k=0;k<=ilen;k++) {
			int ii = indexI+k;
			if (ii >= len) continue;
			colSum += m[ii][j];
			colN+=1.0;
		}
		if (colN > 0.0) {
			triAvg += colSum/colN;
			triN++;
		}
	}
	if (triN > 0.0) {
		//triAvg /= triN;
	}

	double triAvgUp = 0.0;
	double triAvgDn = 0.0;
	double triNUp = 0.0;
	double triNDn = 0.0;
	int upOffset = -1;
	int dnOffset = indexJ;
	for (int j=startIndex;j<=indexJ;j++) {
		double colSumUp = 0.0;
		double colNUp = 0.0;
		double colSumDn = 0.0;
		double colNDn = 0.0;
		int ilen = j;
		for (int k=0;k<ilen;k++) {
			int iup = indexI-k+upOffset;
			if (iup >= 0) {
				colSumUp += m[iup][j];
				colNUp+=1.0;
			}
			int idn = indexI-k+dnOffset;
			if (idn < len) {
				colSumDn += m[idn][j];
				colNDn+=1.0;
			}
		}
		if (colNUp > 0.0) {
			triAvgUp += colSumUp/colNUp;
			triNUp++;
		}
		if (colNDn > 0.0) {
			triAvgDn += colSumDn/colNDn;
			triNDn++;
		}
	}

	//if (triNUp > 0) triAvgUp /= triNUp;
	//if (triNDn > 0) triAvgDn /= triNDn;
	double denominator = (triAvgUp+triAvgDn)/2.0;
	if (ratioMethod == RELMATRIX_TAD_RATIO_METHOD_MAX) {
		denominator = triAvgUp;
		if (triAvgDn > denominator) denominator = triAvgDn;
	}
	if (denominator < 1e-20) denominator = 1.0;
	//double pc = (pseudoCount*minDistIndex)/(double)(indexJ);
	double pc = mstats[startIndex];
	//double pc = mstats[indexJ];
	//double pc1 = mstats[indexJ];
	//fprintf(stderr, "pc=%lf\tpc1=%lf\n", pc, pc1);
	double ratio = (triAvg+pc)/(denominator+pc);

	return ratio;
}


double** RelMatrix::calculateTriMatrix(int mode) {
	double **tri = new double*[len];
	for (int i=0;i<len;i++) {
		tri[i] = new double[size];
		for (int j=0;j<size;j++) {
			tri[i][j] = 0.0;
		}
	}
	double *colSum = new double[size];
	double *colN = new double[size];

	if (mode == RELMATRIX_TRIAVG) {
		for (int i=0;i<len;i++) {
			//if (i%100==0) fprintf(stderr, "\t\t%d\n", i);
			for (int j=0;j<size;j++) {
				colSum[j] = 0.0;
				colN[j] = 0.0;
			}

			//fprintf(stderr, " %d %d \n", minDistIndex, maxDistIndex);
			for (int j=minDistIndex;j<maxDistIndex;j++) {
				double triAvg = 0.0;
				double triN = 0.0;
				for (int k=minDistIndex;k<=j;k++) {
					int ii = i+j-k;
					if (ii >= len) continue;
					colSum[k] += m[ii][k];
					colN[k]++;
					if (colN[k] > 0) {
						triAvg += colSum[k]/colN[k];
						triN++;
					}
				}
				if (triN > 0) {
					//triAvg /= triN;
				}
				if (i>=domainMinSizeIndex) tri[i][j] = triAvg;
			}
		}
	} else if (mode == RELMATRIX_TRIUPAVG) {

		int adj = -1;

		for (int i=0;i<len;i++) {
			for (int j=0;j<size;j++) {
				colSum[j] = 0.0;
				colN[j] = 0.0;
			}

			double triAvg = 0.0;
			double triN = 0.0;
			for (int j=minDistIndex;j<maxDistIndex;j++) {

				for (int k=0;k<j;k++) {
					int ii = i+adj-k;
					if (ii >= len) continue;
					if (ii <= 0) continue;
					colSum[j] += m[ii][j];
					colN[j]++;
				}
				if (colN[j] > 0) {
					triAvg += colSum[j]/colN[j];
					triN++;
				}
				double curavg = triAvg;
				if (triN > 0) {
					//curavg /= triN;
				}
				if (i>=domainMinSizeIndex) tri[i][j] = curavg;
			}
		}
	} else if (mode == RELMATRIX_TRIDNAVG) {

		for (int i=0;i<len+size;i++) {
			for (int j=0;j<size;j++) {
				colSum[j] = 0.0;
				colN[j] = 0.0;
			}
			double triAvg = 0.0;
			double triN = 0.0;
			for (int j=minDistIndex;j<maxDistIndex;j++) {
				int ii= i-j;
				if (ii >= len) continue;
				if (ii < 0) continue;

				for (int k=0;k<j;k++) {
					int iii = ii+1+k;
					if (iii >= len) continue;
					if (iii <= 0) continue;
					colSum[j] += m[iii][j];
					colN[j]++;
				}
				if (colN[j] > 0) {
					triAvg += colSum[j]/colN[j];
					triN++;
				}
				double curavg = triAvg;
				if (triN > 0) {
					//curavg /= triN;
				}
				if (ii>=domainMinSizeIndex) tri[ii][j] = curavg;
			}
		}
	}
	return tri;
}
void RelMatrix::setMinMax(int min, int max, int domainMin, int win, int minIntDist,
								int intInnerSize, int intOuterSize) {
	minDist = min;
	maxDist = max;
	domainMinSize = domainMin;
	windowSize = win;
	minDistIntIndex = 6;
	if (minIntDist > -1) {
		minDistIntIndex = (int)floor(minIntDist/res);
	}
	if (intInnerSize > -1) {
		innerIntSizeIndex = (int)ceil(intInnerSize/2.0/res);
		fprintf(stderr, "\tInteraction Size set to %d (index=%d)\n", intInnerSize,innerIntSizeIndex);
	}
	if (intOuterSize > -1) {
		outerIntSizeIndex = (int)ceil(intOuterSize/2.0/res);
		fprintf(stderr, "\tInteraction Size set to %d (index=%d)\n", intOuterSize,outerIntSizeIndex);
	}
	if (windowSize < 0) {
		windowSize = res*DEFAULT_WINDOW_SIZE_PIXELS_TAD_REL_MAX;
		//windowSize = 50000;
	}
	if (min < 0) {
		minDist = res * DEFAULT_PIXELS_FROM_DIAGONAL_TO_START_SCORING;
		//minDist = 75000;
	}
	if (domainMinSize < -1) {
		domainMinSize = minDist+win;
		//domainMinSize = 100000;
	} else if (min < 0) {
		//minDist = domainMinSize - windowSize;
		//if (minDist < 0) minDist = 0;
	}


	if (domainMinSize < minDist) domainMinSize = minDist;

	minDistIndex = -1;
	int actualMax = 0;
	for (int i=0;i<size;i++) {
		if (offsets[i] < minDist) {
			minDistIndex = i+1;
		}
		if (offsets[i] < domainMinSize) {
			domainMinSizeIndex = i+1;
		}
		if (offsets[i] < windowSize) {
			windowSizeIndex = i+1;
		}
		if (offsets[i] <= maxDist) {
			maxDistIndex = i;
			actualMax = offsets[i];
		}
	}
	maxDist=actualMax;
	pseudoCount = mstats[minDistIndex];
	fprintf(stderr, "\tMin Distance Index: %d (%d)\n",minDistIndex, minDist);
	fprintf(stderr, "\tMax Distance Index: %d (%d)\n",maxDistIndex, maxDist);
	fprintf(stderr, "\tMin TAD Size Index: %d (%d)\n",domainMinSizeIndex, domainMinSize);
	fprintf(stderr, "\tOptima size Index:  %d (%d)\n",windowSizeIndex,windowSize); 
	fprintf(stderr, "\tPseudoCount at minDist: %lf (%d)\n",pseudoCount,minDistIndex); 
	fprintf(stderr, "\tMin Int Dist Index: %d (%d)\n", minDistIntIndex,minIntDist);
	if (mThreshold != NULL) {
		fprintf(stderr, "\tMin int reads at maxDist: %lf\n", mThreshold[size-1]);
	}
}
void RelMatrix::getMatrixStats() {
	if (mstats == NULL) {
		mstats = new double[size];
	}
	for (int j=0;j<size;j++) {
		mstats[j] = 0.0;
		for (int i=0;i<len;i++) {
			mstats[j] += m[i][j];
		}
		mstats[j] /= len;
		//fprintf(stderr, "%d\t%lf\n", j, mstats[j]);
	}

	if (mThreshold == NULL) {
		mThreshold = new double[size];
	}
	for (int j=0;j<size;j++) {
		mThreshold[j] = 0.0;
		if (poissonThreshold < 1.0) {
			double logp = log(poissonThreshold);
			for (int i=0;i<1000000;i+=1) {
				double p = ilogCumulativePoisson(i,mstats[j]);
				//fprintf(stderr, "%d\t%lf\t%d\t%le\t%le\n", j, mstats[j], i , p , logp);
				if (p < logp) {
					mThreshold[j] = (double)i;
					break;
				}
			}
			//fprintf(stderr, "%d\t%lf\t%lf\n", j, mstats[j], mThreshold[j]);
			//fprintf(stderr, "\tMin interactions at maximum CP distance = %.1lf\n", mThreshold[maxDist-1]);
		}
	}
	for (int i=minDistIndex;i<minDistIndex+50;i++) {
		//fprintf(stderr, "\t%d\t%d\t%lf\t%lf\n", i, i*res, mstats[i],mThreshold[i]);
	}

}
void RelMatrix::readRelMatrixFile(char* file) {
	FILE* fp = fopen(file,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Error - could not open relative matrix file: %s\n", file);
		exit(0);
	}
	int lineCount = 0;
	char* buffer = new char[BUFFER];
	char** cols = new char*[MAX_RELMATRIX_COLUMNS];
	char** pcols = new char*[10];
	int numCols = 0;

	len = MAX_RELMATRIX_ROWS;
	m = new double*[len];
	pos = new int[len];
	for (int i=0;i<len;i++) {
		m[i] = NULL;
		pos[i] = -1;
	}
	int zeroCol = 0;
	int curLength = 0;
	while (fgets(buffer, BUFFER, fp) != NULL) {
		split(buffer, cols, numCols,'\t');
		lineCount++;
		if (lineCount == 1) {
			size = (numCols-2)/2+1;
			zeroCol = size+1;
			fprintf(stderr, "\tCheck: ZeroCol = %s\n", cols[zeroCol]);
			offsets = new int[size];
			for (int i=zeroCol;i<numCols;i++) {
				sscanf(cols[i],"%d",&(offsets[i-zeroCol]));
				//fprintf(stderr, "\t\t%d\t%s\n", offsets[i-zeroCol],cols[i]);
			}
			res = offsets[1]-offsets[0];
			fprintf(stderr, "\tres = %d\n", res);
			continue;
		}
		int mIndex = lineCount-2;

		double* row = new double[size];
		for (int i=zeroCol;i<numCols;i++) {
			sscanf(cols[i],"%lf",&(row[i-zeroCol]));
		}
		m[mIndex] = row;
		curLength++;

		split(cols[0],pcols,numCols,'-');
		sscanf(pcols[1],"%d",&(pos[mIndex]));
		if (chr == NULL) {
			chr = new char[strlen(pcols[0])+1];
			strcpy(chr,pcols[0]);
		}
		if (mIndex == 0) {
			startCoordinate = pos[mIndex];
		}
	}
	len = curLength;
	fprintf(stderr, "\tsize = %d\n", size);
	fprintf(stderr, "\tlen = %d\n", len);
	getMatrixStats();
	fclose(fp);
}
void RelMatrix::printRelMatrix(FILE* fp,double **matrix) {
	fprintf(fp, "RelativeMatrix\tPosition");
	for (int i=0;i<size;i++) {
		fprintf(fp, "\t%d", offsets[i]);
	}
	fprintf(fp, "\n");
	for (int i=0;i<len;i++) {
		fprintf(fp, "%s-%d\t%s-%d",chr,pos[i],chr,pos[i]);
		for (int j=0;j<size;j++) {
			fprintf(fp, "\t%.3lf",matrix[i][j]);
		}
		fprintf(fp, "\n");
	}

}
double* RelMatrix::calculateDI(int diSize) {
	int maxIndex = 0;
	if (diSize > 0) {
		diDist = diSize;
	}
	for (int i=0;i<size;i++) {
		if (offsets[i] > diDist) {
			maxIndex = i;
			break;
		}
	}
	if (maxIndex > size) maxIndex = size;
	fprintf(stderr, "\tCalculating DI with distance of %d\n", diDist);
	double* di = new double[len];
	for (int i=0;i<len;i++) {
		di[i] = 0.0;
		double up=RELMATRIX_DI_PSEUDOCOUNT,down=RELMATRIX_DI_PSEUDOCOUNT;
		for (int j=minDistIndex;j<maxIndex;j++) {
			down += m[i][j];
			int ii = i-j;
			if (ii >=0) {
				up += m[ii][j];
			}
		}
		di[i] = down-up;
		//di[i] = log(down/up)/log(2.0);
	}
	return di;
}
void RelMatrix::printDIbedGraph(FILE* fp, double* di) {
	//leave off last data point since it normally goes past the length of the chromosome...
	for (int i=0;i<len-1;i++) {
		fprintf(fp, "%s\t%d\t%d\t%lf\n", chr, pos[i], pos[i]+res,di[i]);
	}
}
double* RelMatrix::calculateInsulationRatio(int insSize,double** triAvg,
							double** triUpAvg, double **triDownAvg) {
	int maxIndex = 0;
	if (insSize > 0) {
		insDist = insSize;
	}
	for (int i=0;i<size;i++) {
		if (offsets[i] > insDist) {
			maxIndex = i;
			break;
		}
	}
	if (maxIndex > size) maxIndex = size;
	fprintf(stderr, "\tCalculating insulation ratio with a max distance of %d\n", insDist);
	double* ins = new double[len];
	for (int i=0;i<len;i++) {
		ins[i] = 0.0;
		if (i < maxIndex) continue;
		if (i > len-maxIndex) continue;
		double crossInteractions = triUpAvg[i][maxIndex];
		double upstream = triAvg[i-maxIndex][maxIndex];
		double downstream = triAvg[i][maxIndex];
		ins[i] = log((((upstream+downstream)/2.0)+1.0)/(crossInteractions+1.0))/log(2.0);
	}
	return ins;
}



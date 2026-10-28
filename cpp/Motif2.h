
// Copyright 2009, 2010, 2011, 2012 Christopher Benner <cbenner@gmail.com>
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


#include <stdio.h>
#include "Hashtable.h"
#include "statistics.h"
#include <limits.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <atomic>
#include <sys/time.h>

#ifndef MOTIF_H
#define MOTIF_H

#define MOTIF2_BUFFER 10000100

#define MOTIF_STRAND_POS 0
#define MOTIF_STRAND_NEG 1
#define MOTIF_STRAND_BOTH 2
#define MOTIF_STRAND_SEPARATE 3

#define MOTIF_SCORING_HYPERGEOMETRIC 0
#define MOTIF_SCORING_BINOMIAL 1
#define MOTIF_SCORING_ZSCORE 2
#define MOTIF_MISMATCH_PVALUE_MINIMUM -1.3863
//-0.6931 = 0.50
//-1.3863 = 0.25
//-2.3026 = 0.10
//-2.9957 = 0.05
//-4.6052 = 0.01
//-6.9078 = 0.001

#define MOTIF_MAX_BACKGROUND_PERCENT 0.5

#define MOTIF_SCORE_OLIGOS 0
#define MOTIF_REMOVE_OLIGOS 1

#define MOTIF_MISMATCH_MINIMUM_SKIP 1
#define MOTIF_MISMATCH_MINIMUM_COUNTS 1
#define MOTIF_MAX_SEED_OLIGOS_FRACTION 0.10
#define MOTIF_MAX_SEED_OLIGOS_TOTAL 1e5

#define MAX_OLIGO_NORMALIZATION_ITERATIONS 160
#define MIN_OLIGO_NORMALIZATION_IMPROVEMENT 0.01

#define MOTIF_SEED_MINIMUM_LOGP -10.0

#define OLIGO_SEARCH_KEEP 0
#define OLIGO_SEARCH_REMOVE 1

#define MOTIF_DEFAULT_NUMBER_OF_MOTIFS 25

#define MOTIF_INITIAL_DEGENERACY 0.1
#define MOTIF_PROBABILITY_MINIMUM 0.001
#define MOTIF_IUPAC_THRESHOLD 0.15

#define MOTIFS_NUM_POTENTIAL_MOTIFS 10
#define MOTIF_OPTIMZATION_DIFFERENT_THRESHOLDS 0
#define MOTIFS_MAX_OPT_ITERATIONS 10
#define MOTIF_MAX_HIT_RATE_TO_CONSIDER 0.005

#define MAXIMUM_MOTIF_LENGTH 10000
#define MAXIMUM_DENOVO_MOTIF_LENGTH 50

#define DENOVO_SPEED_AUTOADJUST 0
#define DENOVO_SPEED_CONSTANT 1
#define DENOVO_MASK_FULL 0
#define DENOVO_MASK_QUICK 1

#define SCORER_STANDARD_ZOOPS 0
#define SCORER_APPROXIMATE_ZOOPS 1

#define KNOWN_JUST_SCORE_ZOOPS 0
#define KNOWN_OPTIMIZE_THREASHOLD 1
#define KNOWN_MASK_MOTIFS 2
#define KNOWN_MSCORE_MOTIFS 3
#define KNOWN_FIND_ALL 4

#define MOTIF_OFFSET_NULL -12345678

#define MOTIF_OVERLAP_PERCENT 0.75

#define MOTIF_DEFAULT_REDUCE_PERCENTAGE 0.20

#define HOMER_MINIMUM_SEQ_WEIGHT 0.001

#define DNA_ALPHA 4

#define DOUBLE_EPS 1e-20
#define FLOAT_EPS 1e-10
#define TRIAL_CUT_SLACK 0.2
#define TRIAL_CUT_ESCALATE 0.5
// DeNovoMotifs::optimizeBaseThresholdAtCut cuts this fraction of the motif's score range below
// where optimizeMotifMatrix's loop can first stop.
#ifndef BASE_CUT_SLACK
#define BASE_CUT_SLACK 0.02
#endif

#define MAX_STAT_CACHE_SIZE 500

#define HOMER_VERBOSE_LEVEL0 0
#define HOMER_VERBOSE_LEVEL1 1

#define MAX_CHR_LENGTH 400000000


class DeNovoMotifs;
class KnownMotifs;
class OligoArray;
class Oligo;
class Motif;
class EnrichmentScorer;
struct ThreadArgs_findMisMatchEnrichmentThread;
struct ThreadArgs_optimizeMotifThresholdThread;
struct ThreadArgs_scanMotifThread;
class SequenceStats;

// One oligo's leaf-test data: the sequence as a 2-bit code (base i at bits
// 2i, A=0 C=1 G=2 T=3), its counts (kept equal to the Oligo's) and the Oligo.
struct OligoRecord {
	uint64_t code;
	double numTarget;
	double numBackground;
	Oligo* oligo;
};

// An oligo and its sort key, for sorting by key without reading the Oligo.
struct KeyedOligo {
	double key;
	Oligo* oligo;
};

// A trial thread's state for DeNovoMotifs::optimizeMotifThresholdAtCut: the fraction of the next
// matrix's score range above minimumScore at which to cut, and its hit buffer. Padded so CPUs
// do not share a cache line.
struct TrialCut {
	double frac;
	KeyedOligo* hits;
	int maxHits;
	char pad[64];
};

class OligoArray {
public:

	int* tree; //if tree instances are 0 or negative, they reference oligos, 1=tree root
	int numTree;
	int maxTree;
	int alphaSize;
	int oligoLength;

	char* alpha;
	char* rvAlpha;
	char* rvMap;
	char* alphaIndex;
	char* rvAlphaIndex;
	
	double totalTarget;
	double totalBackground;

	Oligo** staticOligos;
	int numOligos;

	Oligo** activeOligos;
	int numActiveOligos;

	int maxOligos;
	int rvFlag;
	int maxAutoNormIters;

	//stack variables for "recursive" looping - save's on new and delete time to keep them here

	int numCPUs;

	Oligo*** utilityOligos;
	int* numUtilityOligos;

	int** indexStack;
	int** seqIndex;
	int** alphaStack;
	int** positionStack;
	int** misMatchStack;
	int** fwdMisMatchStack;
	int** rvMisMatchStack;
	double** fwdScoreStack;
	double** rvScoreStack;
	// One allocation holding every CPU's stacks, each CPU's on its own cache lines.
	char* stackBlock;

	// With trackChanges set, adjustOligoInstances lists every oligo whose
	// counts it changes or that it removes from the tree, and clearTreeSlot
	// every oligo a clear unlinks.
	int trackChanges;
	Oligo** changedOligos;
	int numChangedOligos;
	int maxChangedOligos;
	// With snapContributions also set, changedSnap[2i] and [2i+1] hold the
	// listed oligo's contribution to the screen's sums (oligoContribution)
	// before its first change; NaN when it could not be computed.
	int snapContributions;
	double* changedSnap;
	// From scanCountSums on: countsUsable while every count is finite and
	// >= 0; countsExactT/Bg while that channel's counts are integers whose
	// total (countsTotalT/Bg, an upper bound) is below 2^52, so its sums are
	// exact in any order. adjustOligoInstances keeps them; they only clear.
	int countsUsable;
	int countsExactT;
	int countsExactBg;
	double countsTotalT;
	double countsTotalBg;

	// Built by optimizeMemoryUsage for DNA oligos of at most 32 bases, indexed
	// by oligo number; NULL otherwise, and searches compare the strings.
	OligoRecord* records;
	char* recordBlock;
	// remainMask[d]: the low bit of every 2-bit field after position d.
	uint64_t* remainMask;
	// firstHit[cpu][0]: oligo number of the last search's first hit.
	int** firstHit;

	// For the mismatch screen (records, -len 1-12): codeTable[k] is the number
	// of the oligo reachable in the tree whose sequence has lexicographic
	// 2-bit code k (base 0 in the top bits), else 0; codePresent has bit k set
	// for the same codes. clearTreeSlot keeps both equal to the tree.
	// codeTableState: 0 not built, 1 built, -1 not usable (the tree is searched).
	int* codeTable;
	uint64_t* codePresent;
	int codeTableState;
	// Per-CPU scratch for sumMisMatchCounts, 6 lists of codeScratchSize codes each.
	char* codeScratchBlock;
	uint32_t** codeScratch;
	int codeScratchSize;
	int codeScratchCPUs;
	// For the first screen (buildBallSums): per lexicographic code k,
	// ballLevels+1 (target, background) pairs: level 0 the counts of the oligo
	// with code k (0 if none), level r their sums over every code within r
	// mismatches of k. NULL when not built.
	double* ballCells;
	int ballLevels;
	// Per level, the relative bound on |HOMER v5.1's sum - ballSums' sum|; 0 for a
	// channel whose counts are integers with a total below 2^52 (exact sums).
	double ballErrT[3];
	double ballErrBg[3];


	OligoArray();
	~OligoArray();
	void init(int length, int totalDataSetBp, int alpha, int strandFlag,int numCPUs);
	void optimizeMotifThreshold(Motif* motif);
	void optimizeMemoryUsage();
	void printOligo(char* seq, FILE* fp);
	int addOligo(char* seq, double vt, double vbg, double addValue); //addValue = -1 for normal add
	int addSequenceWindows(const char* text, int numSeqs, const uint32_t* seqStart, const int* seqLen,
				const uint64_t* goodStart, const double* pv, const double* nv, int numThreads);
	int searchOligo(char* seq, int maxMisMatches, double& vt, double& vbg,int removeFlag,int cpu);
	void sortActiveOligos();
	void sortUtilityOligos(int cpu);
	void initializeActiveOligoArray();
	void initializeOligoPvalues(EnrichmentScorer* scorer);
	void normalizeOligos(int nlen, char* normFile);
	double normalizeOligosIteration(int nlen, char* normFile,int printFlag);
	void scoreOligosWithMotif(Motif* motif,double minimumScore, int removeFlag,int cpu);
	void scoreSortUtilityOligos(Motif* motif, double minimumScore, int cpu);
	int scoreOligosAtCut(Motif* motif, double pruneAt, double cut, int cpu, KeyedOligo* &hits, int &maxHits);
	void scoreOligosWithMotifRecords(Motif* motif,double minimumScore, int removeFlag,int cpu);
	int scoreOligosAtCutRecords(Motif* motif, double pruneAt, double cut, int cpu, KeyedOligo* &hits, int &maxHits);
	int scoreOligosAtCutTagged(Motif* motif, double pruneAt, double cut, KeyedOligo* &hits, int &maxHits);
	void calculateFoldEnrichment(double pseudoCount);
	void adjustOligoInstances(char* oligoSeq, double targetValue, double backgroundValue);
	void recordChangedOligo(Oligo* oligo);
	void clearChangedOligos();
	int oligoContribution(Oligo* oligo, double& vt, double& vbg);
	void scanCountSums();
	void readOligoCountFile(char* file, int strand, int cpus);
	void buildRecords();
	void dropRecords();
	int prepareCodeTable(int maxMisMatches);
	void dropCodeTable();
	void clearTreeSlot(int slot);
	void sumMisMatchCounts(char* seq, int maxMisMatches, double& vt, double& vbg, int cpu);
	int findExactOligo(char* seq, int cpu);
	int ballSumsWorthBuilding(int maxMisMatches, int fromMisMatches, long long numQueries);
	int buildBallSums(int maxMisMatches, int numThreads);
	void dropBallSums();
	int ballSums(const char* seq, int misMatches, double minCounts, double& vt, double& vbg, int& above);

private:
	int addOligoTree(char* seq, double vt, double vbg, double addvalue,int rvMode, int knownID);
	int searchOligoTree(char* seq,int maxMisMatches, double& vt, double& vbg);
	void buildCodeTable();
	void unsetCode(int oligoIndex);
	void dropTreeLeaf(int oligoIndex);
	void dropLeavesBelow(int node);

};


class Oligo {
public:
	char* seq;
	double numTarget;
	double numBackground;
	float normWeight; // used if -olen is used to normalize short oligos
	float logp;  // used to figure out contribution to motifs when optimizing matrix
	double value; //generic value - used for different things
					// pvalue during findMisMatchEnrichment
	char flag; //flag used to make sure we don't doulbe count oligos when traversing tree
	char flag2; //generic flag - used to mark oligos to skip during mismatch analysis
	char changed; // set while OligoArray::trackChanges records this oligo's counts changing
	int screenSlot; // index of this oligo's saved mismatch screen result (DeNovoMotifs::saveScreenResults)
	double* cpuData;
	Oligo();
	Oligo(char* seq);
	~Oligo();
	void init();
	void print(FILE*);
	static int oligoCmp(const void* a, const void* b);
	static int oligoCmpReverse(const void* a, const void* b);
	static void revopp(char* seq);
};

// What Sequence::motifScan needs from a motif, built once per motif: the
// log-odds as two flattened tables (forward and reverse complement, 5 entries
// per position, the fifth -inf for non-ACGT bases), and the pruning bound's
// positions and the best total of the others.
class MotifScanTables {
public:
	int length;
	int K;
	double restOfBest;
	double* fwd;
	double* rev;
	int* sel;
	double* gap;
	// how far rounding may move a window's score or bound: 1e-9 of the sum
	// of the largest magnitude at every position, plus 1e-9
	double margin;
	int size;
	// find's bound: up to two blocks of 4 or 5 positions, each with a
	// 3125-entry table of its partial log-odds per 5-mer (base 5) starting at
	// the block, forward and reverse complement, plus the best total of every
	// position outside the blocks. A block of 4 ignores the fifth base. With
	// one block the second table is all zeros at the same offset.
	int nblocks;
	int blockStart[2];
	int blockLen[2];
	double blockRest;
	double* bfwd;
	double* brev;
	// find's index: per strand (0 forward, 1 reverse), the block whose
	// passing 4-mers are fewest and those 4-mers; indexable is 0 when a
	// block table holds a NaN, and find then scans every window.
	int indexable;
	int pivot[2];
	int npass[2];
	int* pass[2];
	// find's 6-mer index, for widths 6 and more (buildSix): s6n blocks, block 0
	// of 6 positions, block 1 of 6 (width 12+) or of the width-6 positions
	// left (widths 7-11), plus the best total of every position outside them.
	// A block of 6 has a 4097-entry table per strand, indexed by kmer6 (entry
	// 4096, an N, is -inf); a shorter block a 3125-entry table per strand
	// indexed by kmer5, the unused digits ignored. Per strand, the pivot (a
	// block of 6) and its passing 6-mers with their partial scores.
	// The pivot has sixP = 6 or 7 positions; blocks of 6 or 7 have a
	// table of 4^k + 1 entries in t6 (the last, an N, -inf), read through
	// k-mer codes of their length; shorter blocks a 3125-entry o5 table.
	int six;
	int sixP;
	int s6n;
	int s6Start[3];
	int s6Len[3];
	double s6Rest;
	double* t6[2][3];
	double* o5[2][3];
	// with the 7-mer pivot, a block of 5 or fewer positions reads
	// ck7's prefix code >> 4: o5x has the block's partial score of each of
	// the 4096 slots (-inf for a prefix shorter than the block)
	double* o5x[2][3];
	int s6pivot[2];
	int s6npass[2];
	int* s6pass[2];
	double* s6passScore[2];
	// per passing pivot k-mer, its score less the need; a window
	// passes the flank test when that margin plus s6flankL[st][b1 + 4*b2] of
	// the two bases before the k-mer (b1 the nearer) plus
	// s6flankR[st][b1 + 4*b2 + 16*b3] of the three after is at least 0
	double* s6passGap[2];
	double s6flankL[2][16];
	double s6flankR[2][64];
	// a block of k = 6 or 7's table is halfHi[st][q][x >> 2(k-3)] +
	// halfLo[st][q][x & (4^(k-3) - 1)]; the shorter block's table repeats
	// o5Width values o5Unused times each
	double halfHi[2][3][64];
	double halfLo[2][3][256];
	int o5Width[3];
	int o5Unused[3];
	// per strand: each block's best partial score, and the blocks
	// other than the pivot, the larger-gap one first
	double s6max[2][3];
	int s6oth[2][2];

	MotifScanTables();
	~MotifScanTables();
	// blocks = 0 skips the 5-mer blocks and their index (find's 6-mer index
	// replaces them)
	void build(Motif* motif, int blocks = 1);
	int buildSix(Motif* motif, int P, const int* countP);
};

class Sequence {
public:
	char* s;
	int g; // group, i.e. target (=1) or background (=0)
	double w; // i.e. 1.0
	char* name;
	int length;
	int id;
	// Base-5 index of the 5-mer starting at each position (N past the end),
	// for find's block bound, and each base's table column; NULL until
	// SequenceArray::prepareKmers().
	unsigned short* kmer5;
	// base-4 6-mer starting at each position; 4096 when it holds an N or
	// runs past the end
	unsigned short* kmer6;
	unsigned char* codes;
	// codes points into a chunk's block (KmerChunk::codesBlock), which
	// the chunk frees
	int codesShared;
	// name and s point into a SequenceArray text block, which the
	// array frees
	int textShared;

	Sequence();
	Sequence(char* name, char* seq, int group, double weight);
	~Sequence();
	void init();
	void setSeq(char* name, char* seq, int group, double weight);
	// setSeq of a name and a sequence given by length, the sequence cleaned
	// (cleanSequence) as it is copied
	void setSeqCleaned(const char* name, size_t nameLen, const char* seq, size_t seqLen,
			int group, double weight);
	// setSeqCleaned with the name and then the sequence written at
	// `at` (nameLen + seqLen + 2 bytes), which the Sequence does not own
	void setSeqCleanedAt(const char* name, size_t nameLen, const char* seq, size_t seqLen,
			int group, double weight, char* at);
	// codes, and kmer5 when withKmer5 (the 5-mer index and the per-window
	// find path read it)
	// codes written at codesAt when it is not NULL (not owned)
	void prepareKmers(int withKmer5 = 1, unsigned char* codesAt = NULL);
	int motifScan(Motif* m,int strand, int offset,int mode, double &bestScore, 
					int* sitePositions, char* siteStrands, FILE* fp, pthread_mutex_t fileLock,
					const MotifScanTables* tables = NULL);
	int scanCandidates(Motif* motif, int strand, int offset, double &bestScore,
					int* sitePositions, char* siteStrands, FILE* fp,
					const MotifScanTables* tables, const int* fwdCand, int nf,
					const int* revCand, int nr);
};
// find's inverted index over consecutive sequences first..first+count-1:
// for each base-5 4-mer value v with no N, the entries start[v]..start[v+1]-1
// of pos (start of the 4-mer in its sequence) and sid (its sequence, minus
// first).
class KmerChunk {
public:
	int first;
	int count;
	int* start;
	int* pos;
	unsigned short* sid;
	// 6-mer chunks: the sequences' kmer6 and kmer5 codes end to end, from
	// cstart[q] (length clen[q]), with `pad` codes holding an N before the
	// first and after the last; pos then indexes these arrays
	int* cstart;
	int* clen;
	unsigned short* ck7;
	unsigned short* ck6;
	unsigned short* ck5;
	// 6/7-mer chunks keep no sid: blockSeq[p >> 6] is the last sequence that
	// starts at or before position p rounded down to 64 (code arrays' index)
	unsigned short* blockSeq;
	// pos packs the entry's flanking bases when flanks: bits 0-21 the
	// position, 22-25 the two bases before the k-mer (b1 + 4*b2, b1 the
	// nearer), 26-31 the three after (b1 + 4*b2 + 16*b3), 2 bits a base,
	// an N or a base past the sequence read as A (0); otherwise pos is the
	// position alone (a chunk of 2^22 positions or more)
	int flanks;
	// the codes of the chunk's sequences (base chunks only)
	unsigned char* codesBlock;

	KmerChunk();
	~KmerChunk();
};
class LinkedSequence {
public:
	Sequence* s;
	LinkedSequence* link;
	
	LinkedSequence(Sequence* seq, LinkedSequence* link);
	~LinkedSequence();
};

class SequenceArray {
public:
	Sequence** seqs;
	int numSeqs;
	LinkedSequence* link;
	double numTarget;
	double numBackground;
	int bpTarget;
	int bpBackground;
	int longestSequence;
	EnrichmentScorer* scorer;	
	double maxBackgroundPercent;
	int maxAutoNormIters;
	KmerChunk* chunks;
	int numChunks;
	// blocks of sequence names and texts, and of Sequence objects
	// built in place (objBlocks[i] holds objCounts[i]); freed with free()
	char** textBlocks;
	int numTextBlocks;
	Sequence** objBlocks;
	long* objCounts;
	int numObjBlocks;
	// chunks holds the 5-mer index when index5; chunks6 (same first/count)
	// the 6-mer index, with count6 its total per 6-mer
	int index5;
	KmerChunk* chunks6;
	int* count6;
	// the 7-mer index (widths 7 and more), with count7 its total per 7-mer
	KmerChunk* chunks7;
	int* count7;
	// the most index entries in one 6-mer chunk
	int maxChunk6Entries;
	// find: nothing reads the per-sequence site stats (no -mout, no
	// -siteReduce), so motifScan skips them
	int seqStatsUnread;
	// optimizeMotifThreshold's scan, for the mask scan that follows on
	// the same motif: each sequence's best window score, the forward log-odds
	// table and margin it scanned with; keptMotif NULL when nothing is kept
	Motif* keptMotif;
	int keptStrand;
	int keptNumSeqs;
	int keptLength;
	double keptMargin;
	double* keptBest;
	int keptBestSize;
	double* keptFwd;
	int keptFwdSize;

	SequenceArray();
	~SequenceArray();
	void deleteLinks();
	void addSequence(Sequence* s);
	void prepareKmers(int threads = 1, int need5 = 1, int need6 = 0, int maxWidth = 0, int need7 = 0);
	void maskMotifs(Motif** motifs, int numMotifs, int strand, int offset, int threads);
	void optimizeArray();
	void readSequenceFiles(char* seqFile, char* groupFile);
	void readFastaDirect(char* inputFASTAfile, char* bgFASTAfile, int threads = 1);
	void printGroupFile(FILE* fp);
	void printSequenceFile(FILE* fp);
	void motifScan(Motif* m,int strand, int offset,int mode,FILE* fp, pthread_mutex_t fileLock);
	void normalizeSequence(int nlen, int strand, char* normFile, int equalFlag, int threads = 1);
	double normalizeSequenceIteration(int nlen, int strand, char* normFile, int printFlag,int neutralFlag);
	void maskMotifFromOligoArray(Motif* motif,OligoArray* oligoArray,char strand);
	void initializeScorer(int scoringMethod,int zoopsFlag,unsigned long long int maxCacheSize);
	void optimizeMotifThreshold(Motif* motif, char strand);
	void parseFasta2SeqAndGroupFiles(char* inputFASTAfile,char* bgFASTAfile,char* sFile, char* gFile);

};


class Motif {
public:
	double** m;
	int length;
	double threshold;
	int logFlag;
	char* name;
	char* consensus;
	double numTarget;
	double numBackground;
	double freqTarget;
	double freqBackground;
	double numSeqTarget;
	double numSeqBackground;
	double logp;
	double minimumScore;
	int badFlag;

	double similarPercent;
	double similarPercentOfOther;
	Motif* similarMotif;

	SequenceStats* seqStats;
	int numSeqs;

	double statAvgTarget;
	double statAvgBackground;
	double statStdTarget;
	double statStdBackground;
	double statStrandLogRatio;
	double siteMultiplicity;


	double* fwdMaxPossible;
	double* rvMaxPossible;

	static int alphaSize;
	static char* alphaIndex;

	Motif();
	Motif(char* seq);
	Motif(Motif*);
	Motif(char* name, char* consensus, double threshold, double logp, double** matrix, int mlen);
	~Motif();
	void init();
	void setName(char*);
	void logxform();
	void expxform();
	void normalize();
	void justNormalize();
	void print(FILE*);
	void getConsensus();
	char getSymbol(char*);
	void setMinimumScore(double percent);
	void setBlank(int len);
	void addOligo(Oligo* o);
	void initializeSeqStats(int);
	void deleteSeqStats();
	void optimizeThreshold(EnrichmentScorer* scorer, double maxBackPercent);
	void adjustLength(int);
	void calculateDistributionStats();

	static double scoreOverlap(Motif* m1, Motif* m2, double &percent1, double &percent2);
	static int motifCmp(const void* a, const void* b);
	static int nameID;
	static Motif** readMotifFile(char* file, int &numMotifs);
};


class DeNovoMotifs {
public:
	OligoArray* oligoArray;
	unsigned long long int totalOligos;
	unsigned long long int numOligoEstimate;
	unsigned int oligoLength;
	unsigned int normLength;
	int bufferLength;
	unsigned char alphaSize;
	unsigned char strand;
	double maximumExpectedPerBp;
	unsigned long long int maxCacheSize;
	int noDuplicateOligos;

	unsigned int numCPUs;
	int mutexIndex;
	int mutexTotal;
	int nextMeter;
	int incMeter;
	pthread_mutex_t mutex;
	std::atomic<long long> screenNext; // next unclaimed screenOligos index
	Oligo** screenOligos; // the oligos findMisMatchEnrichmentThread screens
	int numScreenOligos;
	int screenBallFlag; // set for the first screen: screenMisMatches may use OligoArray::buildBallSums
	// Each active oligo's (value, flag2) after the last mismatch screen.
	Oligo** savedScreenOligos;
	double* savedScreenValues;
	char* savedScreenFlag2;
	int numSavedScreen;
	// Per saved slot and pass j = 1..sumPasses, the mismatch sums as of the
	// last screen: savedSums[(slot*sumPasses+j-1)*4 ..+3] = target, background,
	// and bounds on their distance from the exact sums of the current counts,
	// valid while bit j-1 of savedSumValid[slot] is set. savedCodes[slot] is
	// the oligo's 2-bit code when SAVED_CODE_BIT is set. sumPasses 0: not kept.
	double* savedSums;
	unsigned char* savedSumValid;
	unsigned long long* savedCodes;
	int sumPasses;
	int rescreenSumsFlag; // set during rescreenAfterMask's screen
	int markKeyBlocks; // blocks per key in markRescreenOligos, 1-4; 0 picks by cost
	// Per pass, the relative bound on HOMER v5.1's own rounding of a sum; 0 for an exact channel.
	double sumErrT[8];
	double sumErrBg[8];
	// One per CPU for the trial threads' optimizeMotifThreshold; NULL outside optimizeSeeds2Motifs' full masking.
	TrialCut* trialCuts;

	double numTargets;
	double numBackground;
	Oligo** seedOligos;
	int numSeedOligos;

	Motif** motifs;
	int numMotifs;

	SequenceArray* fullSequences;

	char* seqFile;
	char* groupFile;
	char* inputFASTAfile;
	char* bgFASTAfile;
	int fastaFlag;
	char* normFile;
	char* tmpFile;
	FILE* sitefp;

	int randFlag;

	int speedFlag;
	double siteReduceThreshold;

	int optimizeGivenFlag;
	int scoringMethod;
	int finalZoopsFlag;
	double maxBackgroundPercent;

	int maxMisMatches;
	int misMatchTargetOnlyFlag;
	double misMatchMinPvalue;
	int minimumMisMatchToStartSkipping;
	double minimumCountsToStartSkipping;
	double minimumFractionToStartSkipping;

	double minimumSeedLogp;
	int totalMotifsToFind;
	int numTrialMotifs;
	int trialCPUs;
	int numConcurrentMotifs;
	int maxOptimizationIterations;
	int motifOptimizationMethod;
	int maskMethod;
	int maxAutoNormIters;

	int zoopsFlag;
	EnrichmentScorer* scorer;	
	//double (DeNovoMotifs::*enrichmentFunction)(double pv, double nv);

	DeNovoMotifs();
	~DeNovoMotifs();
	
	void init();
	void readSequence2Tree();
	void loadFullSequences();
	void checkParameters();
	void getOligoEnrichment();
	void normalizeOligos();


	void findMisMatchEnrichment(int maxMisMatches,FILE* fp,int verboseLevel);
	void screenMisMatches(int maxMisMatches,int verboseLevel);
	void rescreenAfterMask(int maxMisMatches);
	int markRescreenOligos(char* dirty, int maxMisMatches);
	void saveScreenResults();
	void freeScreenResults();
	void setupStoredSums();
	void dropStoredSums();
	int storedSums(Oligo* o, int pass, double& vt, double& vbg, int& above);
	void storeSums(Oligo* o, int pass, double vt, double vbg);
	void findMisMatchEnrichmentThread(ThreadArgs_findMisMatchEnrichmentThread* threadArgs);
	void optimizeMotifMatrix(Motif* &motif);
	void optimizeMotifThreshold(Motif* m, int cpu, int trial=0, int loopFollows=0);
	int optimizeMotifThresholdAtCut(Motif* m, int cpu);
	int optimizeBaseThresholdAtCut(Motif* m);
	void setMotifThreshold(Motif* m, double bestThreshold, double bestScore, double bestTvalue, double bestBGvalue);
	void setTrialCut(Motif* m, int cpu, int breakIndex, int numHits, double bestKey, double nextKey,
				double bestThreshold);
	void optimizeMotifThresholdThread(ThreadArgs_optimizeMotifThresholdThread*);
	void initializeEnrichmentScoring();
	double scoreEnrichmentHypergeometric(double pv, double nv);
	double scoreEnrichmentBinomial(double pv, double nv);
	void getSeedOligos();
	void optimizeSeeds2Motifs();
	void optimizeGivenMotifs();
	void scoreWithKnownMotifs(int optFlag);
	void freeOligos();
	void parseFasta();
	void printMotifs(FILE* fp); // sorts motifs array
};

class EnrichmentScorer {
public:

	unsigned int numTargets;
	unsigned int numBackground;

	std::atomic<double*>* cache;
	int cacheLength;
	unsigned long long int maxCacheSize;

	int zoopsFlag;
	double* targetZoopsApprox;
	double* bgZoopsApprox;
	unsigned int maxZoopsApprox;

	double (EnrichmentScorer::*enrichmentFunction)(double pv, double nv);

	EnrichmentScorer();
	~EnrichmentScorer();
	double scoreEnrichment(double pv, double nv);
	void initializeEnrichmentScoring(int method, double numTarget, double numBack, int zoopsFlag, 
									unsigned long long int newMaxCacheSize);
	double scoreEnrichmentHypergeometric(double pv, double nv);
	double scoreEnrichmentBinomial(double pv, double nv);
	double scoreEnrichmentZscore(double pv, double nv);
	double getZoopsApproxTarget(double value);
	double getZoopsApproxBackground(double value);
	double* getCacheRow(unsigned int mtotal);
};


class KnownMotifs {
public:
	SequenceArray* seqArray;
	unsigned long long int totalBp;
	unsigned int normLength;
	unsigned char alphaSize;
	unsigned char strand;
	unsigned long long int maxCacheSize;
	int offset;

	double numTargets;
	double numBackground;

	Motif** motifs;
	int numMotifs;

	unsigned int numCPUs;
	int mutexIndex;
	int mutexTotal;
	double nextMeter;
	double incMeter;
	pthread_mutex_t mutex;
	pthread_mutex_t fileLock;

	char* seqFile;
	char* groupFile;
	char* inputFASTAfile;
	char* bgFASTAfile;
	int fastaFlag;
	// -i FASTA is read in memory by loadSequence, not through temp files
	int directFasta;
	char* knownFile;
	char* normFile;
	FILE* sitefp;
	FILE* normfp;
	char* motifFilename;

	double maxBackgroundPercent;
	double siteReduceThreshold;

	int scoringMethod;
	EnrichmentScorer* scorer;	

	int mscoreFlag;
	int threasholdOptimizationFlag;
	int zoopsFlag;

	int knownFlag;
	int findFlag;
	int maskFlag;
	int normOnlyFlag;
	int maxAutoNormIters;
	int neutralFlag;

	KnownMotifs();
	~KnownMotifs();
	void checkParameters();
	void loadMotifs();
	void loadSequence();
	void parseFasta();
	void getMotifEnrichment();
	void printMotifEnrichment(FILE* fp);
	void printMotifs(FILE* fp);
	void scanMotifThread(ThreadArgs_scanMotifThread* args);
	void findOverlappingMotifs();
	
};

class GroupInfo {
public:
	double weight;
	unsigned int group;
};

class SequenceStats {
public:
	Sequence* seq;
	double score;
	int numSites;
	int* pos;
	char* strands;
	SequenceStats();
	~SequenceStats();
	void initialize();
	void setSites(int* positions, char* strands, int num);
	static int scoreCmp(const void* a, const void* b);
	static int defaultOrder(const void* a, const void* b);
};

void splitMotif2(char* string, char** cols, int &numCols, char delim);
void cleanSequence(char* seq);
int checkSequence(char* seq);
int sortDouble(const void* a, const void* b);

void printCMDhomer();
void printCMDdenovo();
void printCMDbackground();
void programDeNovo(int argc, char** argv);
void programBackground(int argc, char** argv);
void programKnown(int argc, char** argv);


struct ThreadArgs_findMisMatchEnrichmentThread {
public:
	DeNovoMotifs* denovo;
	int cpu;
	int maxMisMatches;
	int totalChecked;
	int numSkipped;
	int dueToTargetTags;
	int dueToLessThanLast;
	int dueToLowPvalue;
	int lastChecked;
	int lastCheckedTotal;
	int verboseLevel;
	ThreadArgs_findMisMatchEnrichmentThread(DeNovoMotifs* obj,int cpu,int mis, int total,int vlevel);
};
void* DeNovoMotifs_findMisMatchEnrichmentThread(void* threadArgs);

struct ThreadArgs_optimizeMotifThresholdThread {
public:
	DeNovoMotifs* denovo;
	int cpu;
	Motif** motifs;
	ThreadArgs_optimizeMotifThresholdThread(DeNovoMotifs* obj,int cpu,Motif** motifs);
};
void* DeNovoMotifs_optimizeMotifThresholdThread(void* threadArgs);

struct ThreadArgs_scanMotifThread {
public:
	KnownMotifs* known;
	int cpu;
	ThreadArgs_scanMotifThread(KnownMotifs* obj,int cpu);
};
void* KnownMotifs_scanMotifThread(void* threadArgs);


class BgSequence {
// None of the char pointers in this class will be deleted.
public:
	char* seq;
	int len;
	char* name; //chr
	int pos; //position
	char strand;
	double stat; // i.e. gc content
	double pscore; // i.e. gc content
	double weight;
	int bin;
	char mask; //0=not masked, 1=masked, 2=perm mask
	long int positionIndex;

	BgSequence();
	~BgSequence();
	char* getSeq();
	void printAll(FILE*);
	void printFASTA(FILE*);
	void printFASTA(FILE*,int);
	void printBED(FILE*);
	void printBED(FILE*,int);
	void printTable(FILE*);
	void printTable(FILE*,int);
	void printGroup(FILE*,int assignment,int);
	void printGroup(FILE*,int assignment);
	void printSeq(FILE*,int);
	void printSeq(FILE*);
	double getGC();
};
class BgSequenceArray {
public:
	BgSequence** seqs;
	long int numSeqs;

	int seqLen;
	int strand;

	char** ogSeqs;
	char** ogSeqsRV;
	char** ogNames;
	long int* ogLens;
	int numOGs;
	long int totalOGSeqLength;
	BgSequence** posStrandArray;
	BgSequence** negStrandArray;
	long int numStrandArray;
	Inttable* chr2Index;
	int shareFlag;

	BgSequence*** bins;
	long int* binSizes;
	int numBins;

	BgSequence*** sbins;
	long int* sbinSizes;
	char* sbinFixed;
	long int numBgSeqs;
	long int numBgSeqPerBin;

	int pkmer;
	int positionFlag;
	int freqLength;
	int alphaLen;
	int numSubBins;

	double ***targetFreqsByBin;

	int maxFailures;
	int overlapFlag;

	BgSequenceArray();
	~BgSequenceArray();

	void ReadFastaFile(char* filename,char repMaskFlag);
	void processGenomicPositions(char* bedFileName, char strandSetting, double MaxNpercent);
	void processSequencesAsGenome(int seqLength, char strandSetting, double MaxNpercent,
								long int numToConsider, double* binThreshs, int nBins);
	void processGenomicRanges(char* inputBackgroundBedFile,int size,char strandFlag, 
								double maxPercentNperSeq,long int numConsideredBgSeqs,
								double* binThresholds,int numBins);
	void processSequencesAsIs(char strandSetting, double MaxNpercent);
	void maskGenomicPositionOverlaps(BgSequenceArray*, int level);
	void sortBasedOnStat();
	void sortBasedOnCoordinates();
	void sortBasedOnStat(BgSequence** array, long int len);
	static int statCmp(const void *a, const void *b);
	static int coordCmp(const void *a, const void *b);
	double* generateBinThresholds(int numBins);
	void assignSeqsToBins(double* binThresholds,int numBINs);
	void selectAllSequences();

	void selectInitialSequences(long int numBgSeqs);	
	void optimizePositionalFrequncies(int maxIters, double changeRate);
	double optimizePositionalFrequencyBin(int bin, double changeRate);
	void setpKmer(int kmer,int positionFlag);
	double*** calculatePositionalKmerFrequenciesByBin();
	double** calculatePositionalKmerFrequencies(BgSequence** seqArray, long int num);
	double** calculatePositionalKmerFrequencies(BgSequence*** seqArray, long int* sizes, int num);
	void scoreKmerProfiles();
	void scoreSequencesByKmerFreqDifference(BgSequence** bins,long int binSizes,double **freqDiff);
	double* generateSubBinThresholds(int bin);
	void selectSeqsBasedOnSubBins(int bin,double* thresholds, long int numSeqToSelect);
	void equalizesBins(long int numperbin);

	void printPositionalKmerFrequencies(FILE* fp, double** freq);
	long int getTotalSelectedSequences();

	void addSequences(BgSequence** s, long int numSeq, int bin);

	void shareGenome(BgSequenceArray*);

	void printAllSequencesFASTA(FILE*);
	void printAllSequencesTable(FILE*);
	void printAllSequencesBED(FILE*);
	void printSelectedSequencesFASTA(FILE*);
	void printSelectedSequencesTable(FILE*);
	void printSelectedSequencesBED(FILE*);
	void printGroupFile(FILE*,int assignment);
	void printSeqFile(FILE*);

	void addSeqs(BgSequence** seqs, int numSeqs);

	void deleteSBins();
	void mask(BgSequence* bgSeq,int level);
	void clearMask();
	void resetMask();
	void resetMask(int sBinIndex);
	int seq2Index(char* s, int len);
	char* index2Seq(int index, int len); // returns new memory allocated to a string

};
class BgPosition {
public:
	char *chr;
	char *name;
	int start;
	int end;
	char strand;
	int len;
	int chrIndex;
	BgPosition();
	~BgPosition();
	void set(char* n, char* c, int s, int e, char st, int l);
};

class BgModel {
public:
	double** freq; //seqLen x freqLength
	double* counts; //seqLen x freqLength
	double*** fwdFreq; //seqLen x freqLength x alphaLen
	double** fwdCounts; //seqLen x freqLength x alphaLen
	double*** revFreq; //seqLen x freqLength x alphaLen
	double** revCounts; //seqLen x freqLength x alphaLen
	int** fwdIndex; //freqLength x alphaLen
	int** revIndex; //freqLength x alphaLen
	int freqLength;
	int alphaLen;
	int seqLen;
	int modelLen;
	int positionFlag;
	int kmer;

	const char* alphaRef;
	char** freqRef;

	BgModel();
	~BgModel();

	void initializeModel(int kmer, int alphaLen, int positionFlag, int length);
	void initializeFreqIndex();
	void buildModel(BgSequence** seqs, long int numSeqs);
	BgSequence** generateSequences(long int numSeqs);
	char* generateSequence(int len);

	int seq2Index(char* s, int len);
	void index2Seq(int index, int len, char* rv);

};

void matchKmerFrequencies(BgSequenceArray* targets, BgSequenceArray* bg, long int numBgSeqs, int pkmer, int positionFlag, int numSubBins, 
									int numIters,int overlapIter, double decayRate, char* pscoreOutputFile);
void split2(char* string, char** cols, int &numCols, char delim);
char* revopp2(char* seq);
int cmpBgPositions(const void *a, const void *b);

#endif

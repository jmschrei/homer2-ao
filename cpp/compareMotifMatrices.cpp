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

// compareMotifMatrices: compareMotifs.pl's compareMotifsRehli with the
// correlation score, for many motif pairs at once (called by compareMotifs.pl).
//
//	compareMotifMatrices [-cpu <#>] < pairs.bin > scores.bin
//
// Input on stdin, native byte order: int32 mode (0: every pair of a query motif
// with a library motif, query-major; 1: every pair i<j of the query motifs),
// int32 norevopp, int32 query motif count, int32 library motif count (0 in
// mode 1), an int32 length per motif (query motifs first), then every motif's
// length*4 doubles, row by row, in the same order.
// Output on stdout, one record per pair: double score, int32 status (0 a score,
// 1 Perl's "NA", 2 not computed here: the caller compares the pair in Perl),
// int32 offset, int32 direction, int32 aligned length. Nothing is written until
// all of stdin has been read.
//
// Scores must equal Perl's to the bit: the sums run in Statistics::correlation's
// order with no fused multiply-add (the Makefile compiles this file with
// -ffp-contract=off). Perl adds and multiplies integer-valued numbers as
// integers, so a pair is computed here only when every value is finite and every
// sum and product stays below 2^53, where both give the same result.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <pthread.h>
#include <vector>

struct MotifMatrix {
	int len;
	const double* rows;  // len*4 values
	int exact;           // 1: every value finite, so the pair bound below applies
	double maxAbs;       // largest |value|, at least 1
};

struct PairResult {
	double score;
	int status;
	int offset;
	int dir;
	int len;
};

struct Job {
	const std::vector<MotifMatrix>* query;
	const std::vector<MotifMatrix>* library;
	int mode;
	int norevopp;
	long long numPairs;
	PairResult* results;
	int thread;
	int numThreads;
};

static const double PAD = 0.25;

// Statistics::correlation's r from its five sums; returns 1 for "NA".
static int correlationFromSums(double xysum, double xsum, double ysum, double x2sum,
		double y2sum, double n, double* r) {
	double numerator = xysum - xsum*ysum/n;
	double denominator = (x2sum - xsum*xsum/n)*(y2sum - ysum*ysum/n);
	if (denominator <= 0) return 1;
	denominator = std::sqrt(denominator);
	if (denominator > 0) {
		*r = numerator / denominator;
		return 0;
	}
	return 1;
}

// compareMotifsRehli(m1, m2): every offset with at least one overlapping column,
// both orientations of m2, 0.25 rows padding both sides to one length.
static void compareRehli(const MotifMatrix& m1, const MotifMatrix& m2, int norevopp,
		std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Yr,
		PairResult* res) {
	int len1 = m1.len;
	int len2 = m2.len;
	if (len1 < 1 || len2 < 1 || !m1.exact || !m2.exact
			|| 4.0*(len1+len2)*(m1.maxAbs > m2.maxAbs ? m1.maxAbs : m2.maxAbs) >= 67108864.0) {
		res->status = 2;
		res->score = 0;
		res->offset = 0;
		res->dir = 0;
		res->len = 0;
		return;
	}

	// X: len2 pad rows, m1, len2 pad rows; Y and Yr: len1 pad rows, m2 (or its
	// reverse opposite), len1 pad rows. Each offset reads a window of these.
	int nx = 4*(len2+len1+len2);
	int ny = 4*(len1+len2+len1);
	X.assign(nx, PAD);
	Y.assign(ny, PAD);
	Yr.assign(ny, PAD);
	memcpy(&(X[4*len2]), m1.rows, sizeof(double)*4*len1);
	memcpy(&(Y[4*len1]), m2.rows, sizeof(double)*4*len2);
	for (int i=0;i<len2;i++) {
		for (int j=0;j<4;j++) {
			Yr[4*(len1+i)+j] = m2.rows[4*(len2-1-i)+3-j];
		}
	}

	double bestScore = -1e10;
	int bestNA = 0;
	int bestOffset = 0;
	int bestDirection = 0;
	int bestLength = 0;
	for (int offset2 = len1-1; offset2 > -len2; offset2--) {
		int a = (offset2 < 0) ? -offset2 : 0;
		int b = offset2+len2-len1;
		if (b < 0) b = 0;
		int c = (offset2 > 0) ? offset2 : 0;
		int L = a+len1+b;
		const double* x = &(X[4*(len2-a)]);
		const double* y = &(Y[4*(len1-c)]);
		const double* yr = &(Yr[4*(len1-c)]);
		double xysum=0, xsum=0, ysum=0, x2sum=0, y2sum=0, xyrsum=0, yrsum=0, yr2sum=0;
		int n = 4*L;
		for (int k=0;k<n;k++) {
			double xv = x[k];
			double yv = y[k];
			double rv = yr[k];
			xysum += xv*yv;
			xsum += xv;
			ysum += yv;
			x2sum += xv*xv;
			y2sum += yv*yv;
			xyrsum += xv*rv;
			yrsum += rv;
			yr2sum += rv*rv;
		}
		double score = 0, rvScore = 0;
		int scoreNA = correlationFromSums(xysum, xsum, ysum, x2sum, y2sum, (double)n, &score);
		int rvNA = correlationFromSums(xyrsum, xsum, yrsum, x2sum, yr2sum, (double)n, &rvScore);
		// Perl compares "NA" as 0
		double s = scoreNA ? 0.0 : score;
		double rs = rvNA ? 0.0 : rvScore;
		double bs = bestNA ? 0.0 : bestScore;
		if (norevopp == 0 && rs > s) {
			if (rs > bs) {
				bestScore = rvScore;
				bestNA = rvNA;
				bestOffset = -offset2;
				bestDirection = 1;
				bestLength = L;
			}
		} else {
			if (s > bs) {
				bestScore = score;
				bestNA = scoreNA;
				bestOffset = -offset2;
				bestDirection = 0;
				bestLength = L;
			}
		}
	}
	res->score = bestNA ? 0.0 : bestScore;
	res->status = bestNA;
	res->offset = -bestOffset;
	res->dir = bestDirection;
	res->len = bestLength;
}

static void* compareThread(void* arg) {
	Job* job = (Job*)arg;
	std::vector<double> X, Y, Yr;
	const std::vector<MotifMatrix>& q = *(job->query);
	long long p = 0;
	if (job->mode == 0) {
		const std::vector<MotifMatrix>& lib = *(job->library);
		long long nl = (long long)lib.size();
		for (p=job->thread;p<job->numPairs;p+=job->numThreads) {
			compareRehli(q[p/nl], lib[p%nl], job->norevopp, X, Y, Yr, &(job->results[p]));
		}
	} else {
		int nq = (int)q.size();
		for (int i=0;i<nq;i++) {
			for (int j=i+1;j<nq;j++,p++) {
				if (p % job->numThreads != job->thread) continue;
				compareRehli(q[i], q[j], job->norevopp, X, Y, Yr, &(job->results[p]));
			}
		}
	}
	return NULL;
}

int main(int argc, char** argv) {
	int numThreads = 1;
	for (int i=1;i<argc;i++) {
		if (strcmp(argv[i],"-cpu")==0 && i+1 < argc) {
			numThreads = atoi(argv[++i]);
		} else {
			fprintf(stderr, "\n\tcompareMotifMatrices [-cpu <#>] < pairs.bin > scores.bin\n");
			fprintf(stderr, "\n\tScores motif pairs as compareMotifs.pl's compareMotifsRehli does"
								" (binary input and output, called by compareMotifs.pl)\n\n");
			return 1;
		}
	}
	if (numThreads < 1) numThreads = 1;

	std::vector<char> input;
	char chunk[1<<16];
	size_t got = 0;
	while ((got = fread(chunk, 1, sizeof(chunk), stdin)) > 0) {
		input.insert(input.end(), chunk, chunk+got);
	}
	size_t size = input.size();
	int header[4] = {0, 0, 0, 0};
	if (size < sizeof(header)) {
		fprintf(stderr, "!!! compareMotifMatrices: bad input header !!!\n");
		return 1;
	}
	memcpy(header, &(input[0]), sizeof(header));
	int mode = header[0], norevopp = header[1], nq = header[2], nl = header[3];
	if ((mode != 0 && mode != 1) || nq < 0 || nl < 0 || (mode == 1 && nl != 0)) {
		fprintf(stderr, "!!! compareMotifMatrices: bad input header !!!\n");
		return 1;
	}
	size_t numMotifs = (size_t)nq + (size_t)nl;
	size_t pos = sizeof(header);
	if ((size - pos)/sizeof(int) < numMotifs) {
		fprintf(stderr, "!!! compareMotifMatrices: truncated input !!!\n");
		return 1;
	}
	std::vector<int> lengths(numMotifs + 1);
	memcpy(&(lengths[0]), &(input[pos]), sizeof(int)*numMotifs);
	pos += sizeof(int)*numMotifs;
	size_t numValues = 0;
	for (size_t i=0;i<numMotifs;i++) {
		if (lengths[i] < 0) {
			fprintf(stderr, "!!! compareMotifMatrices: bad motif length !!!\n");
			return 1;
		}
		numValues += 4*(size_t)lengths[i];
	}
	if ((size - pos) != sizeof(double)*numValues) {
		fprintf(stderr, "!!! compareMotifMatrices: truncated or oversized input !!!\n");
		return 1;
	}
	std::vector<double> values(numValues + 1);
	if (numValues > 0) memcpy(&(values[0]), &(input[pos]), sizeof(double)*numValues);

	std::vector<MotifMatrix> query, library;
	size_t offset = 0;
	for (size_t i=0;i<numMotifs;i++) {
		MotifMatrix m;
		m.len = lengths[i];
		m.rows = &(values[offset]);
		offset += 4*(size_t)m.len;
		m.exact = 1;
		m.maxAbs = 1.0;
		for (int j=0;j<4*m.len;j++) {
			double v = m.rows[j];
			if (!std::isfinite(v)) {
				m.exact = 0;
			} else if (std::fabs(v) > m.maxAbs) {
				m.maxAbs = std::fabs(v);
			}
		}
		if (i < (size_t)nq) {
			query.push_back(m);
		} else {
			library.push_back(m);
		}
	}

	long long numPairs = (mode == 0) ? (long long)nq*nl : (long long)nq*(nq-1)/2;
	std::vector<PairResult> results(numPairs > 0 ? numPairs : 1);
	if (numThreads > numPairs) numThreads = (numPairs > 0) ? (int)numPairs : 1;
	std::vector<Job> jobs(numThreads);
	std::vector<pthread_t> threads(numThreads);
	for (int t=0;t<numThreads;t++) {
		jobs[t].query = &query;
		jobs[t].library = &library;
		jobs[t].mode = mode;
		jobs[t].norevopp = norevopp;
		jobs[t].numPairs = numPairs;
		jobs[t].results = &(results[0]);
		jobs[t].thread = t;
		jobs[t].numThreads = numThreads;
	}
	for (int t=1;t<numThreads;t++) {
		pthread_create(&(threads[t]), NULL, compareThread, &(jobs[t]));
	}
	compareThread(&(jobs[0]));
	for (int t=1;t<numThreads;t++) {
		pthread_join(threads[t], NULL);
	}

	const size_t recordSize = sizeof(double) + 4*sizeof(int);
	std::vector<char> out(recordSize*(numPairs > 0 ? numPairs : 1));
	for (long long p=0;p<numPairs;p++) {
		char* rec = &(out[recordSize*p]);
		int fields[4] = {results[p].status, results[p].offset, results[p].dir, results[p].len};
		memcpy(rec, &(results[p].score), sizeof(double));
		memcpy(rec + sizeof(double), fields, 4*sizeof(int));
	}
	if (numPairs > 0 && fwrite(&(out[0]), recordSize, numPairs, stdout) != (size_t)numPairs) {
		fprintf(stderr, "!!! compareMotifMatrices: could not write results !!!\n");
		return 1;
	}
	if (fflush(stdout) != 0) return 1;
	return 0;
}

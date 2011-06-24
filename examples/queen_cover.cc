
// $Id$

/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2009, Iowa State University Research Foundation, Inc.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
    Builds the set of solutions to the queen cover problem for 
    user-specified board size NxN.

    In other words, finds all possible ways to put queens onto
    an NxN chessboard so that every square either contains a queen,
    or is attached by a queen.
*/

#include <cstdio>
#include <cstdlib>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include "meddly.h"

using namespace MEDDLY;

int N;
FILE* outfile;

dd_edge** qic;
dd_edge** qidp;
dd_edge** qidm;


class timer {
  double start_time;
protected:
  inline double Time() const {
    rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_utime.tv_usec / 1000000.0;
  }
public:
  timer() { 
    reset();
  }
  void reset() { 
    start_time = Time(); 
  }
  double elapsed() const {
    return Time() - start_time;
  }
};


void printmem(long m)
{
  if (m<1024) {
    printf("%ld bytes", m);
    return;
  }
  double approx = m;
  approx /= 1024;
  if (approx < 1024) {
    printf("%3.2lf Kbytes", approx);
    return;
  }
  approx /= 1024;
  if (approx < 1024) {
    printf("%3.2lf Mbytes", approx);
    return;
  }
  approx /= 1024;
  if (approx < 1024) {
    printf("%3.2lf Gbytes", approx);
    return;
  }
  approx /= 1024;
  printf("%3.2lf Tbytes", approx);
}

forest* buildQueenForest()
{
  printf("Initializing domain and forest\n");
  int* vars = new int[N*N];
  for (int i=0; i<N*N; i++) {
    vars[i] = 2;
  }
  domain* d = createDomainBottomUp(vars, N*N);
  assert(d);
  forest* f = d->createForest(false, forest::INTEGER, forest::MULTI_TERMINAL);
  assert(f);

  // Set up MDD options
  f->setReductionRule(forest::FULLY_REDUCED);
  f->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE);
  f->setNodeDeletion(forest::PESSIMISTIC_DELETION);
  
  delete[] vars;
  return f;
}

inline int ijmap(int i, int j)
{
  return i*N+j+1;
}

void hasQueen(forest* f, int i, int j, dd_edge &e)
{
  f->createEdgeForVar(ijmap(i, j), false, e);
}

void queenInRow(forest* f, int r, dd_edge &e)
{
  f->createEdge(int(0), e);
  if (r<0 || r>=N) return;
  for (int j=0; j<N; j++) {
    dd_edge tmp(f);
    hasQueen(f, r, j, tmp);
    apply(MAXIMUM, e, tmp, e);
  }
}

void queenInCol(forest* f, int c, dd_edge &e)
{
  if (c>=0 && c<N && qic[c]) {
    e = *qic[c];
    return;
  }
  f->createEdge(int(0), e);
  if (c<0 || c>=N) return;
  for (int i=0; i<N; i++) {
    dd_edge tmp(f);
    hasQueen(f, i, c, tmp);
    apply(MAXIMUM, e, tmp, e);
  }
  qic[c] = new dd_edge(e);
}

void queenInDiagP(forest* f, int d, dd_edge &e)
{
  if (d>=0 && d<2*N-1 && qidp[d]) {
    e = *qidp[d];
    return;
  }
  f->createEdge(int(0), e);
  if (d<0 || d>=2*N-1) return;
  for (int i=0; i<N; i++) for (int j=0; j<N; j++) if (i+j==d) {
    dd_edge tmp(f);
    hasQueen(f, i, j, tmp);
    apply(MAXIMUM, e, tmp, e);
  }
  qidp[d] = new dd_edge(e);
}

void queenInDiagM(forest* f, int d, dd_edge &e)
{
  if (d>-N && d<N && qidm[d+N-1]) {
    e = *qidm[d+N-1];
    return;
  }
  f->createEdge(int(0), e);
  if (d<=-N || d>=N) return;
  for (int i=0; i<N; i++) for (int j=0; j<N; j++) if (i-j==d) {
    dd_edge tmp(f);
    hasQueen(f, i, j, tmp);
    apply(MAXIMUM, e, tmp, e);
  }
  qidm[d+N-1] = new dd_edge(e);
}

bool processArgs(int argc, const char** argv)
{
  if (argc<2) return false;
  if (argc>3) return false;
  N = atoi(argv[1]); 
  if (N<1) return false;
  if (argc>2) {
    outfile = fopen(argv[2], "w");
    if (0==outfile) {
      printf("Couldn't open %s for writing, no solutions will be written\n", argv[2]);
    }
  } else {
    outfile = 0;
  }
  return true;
}

int usage(const char* who)
{
  printf("Usage: %s N <outfile>\n\n\t        N:  board dimension\n", who);
  printf("\t<outfile>: if specified, we write all solutions to this file\n\n");
  return 1;
}

int main(int argc, const char** argv)
{
  if (!processArgs(argc, argv)) return usage(argv[0]);
  initialize();
  printf("Using %s\n", getLibraryInfo(0));

  timer stopwatch;
  printf("\nDetermining queen covers for %dx%d chessboard.\n", N, N);

  forest* f = buildQueenForest();

  dd_edge num_queens(f);
  f->createEdge(int(0), num_queens);
  for (int i=0; i<N; i++) for (int j=0; j<N; j++) {
    dd_edge tmp(f);
    hasQueen(f, i, j, tmp);
    apply(PLUS, num_queens, tmp, num_queens);
  } // for i,j

  // "By-hand" compute tables for queen in col, diag
  qic = new dd_edge*[N];
  for (int i=0; i<N; i++) qic[i] = 0;
  qidp = new dd_edge*[2*N-1];
  qidm = new dd_edge*[2*N-1];
  for (int i=0; i<2*N-1; i++) {
    qidp[i] = qidm[i] = 0;
  }

  dd_edge solutions(f);

  int q;
  for (q=1; q<=N; q++) {

    printf("\nTrying to cover with %d queens\n", q);

    // Build all possible placements of q queens:
    dd_edge qqueens(f);
    f->createEdge(q, qqueens);
    apply(EQUAL, qqueens, num_queens, qqueens);
    solutions = qqueens;

    printf("\tAdding constraints by row\n\t\t");
    for (int i=0; i<N; i++) {
      printf("%d", N-i);
      fflush(stdout);
      dd_edge rowcov = qqueens;
      dd_edge qir(f);
      queenInRow(f, i, qir);
      for (int j=0; j<N; j++) {
        dd_edge col(f);
        queenInCol(f, j, col);
        dd_edge dgp(f);
        queenInDiagP(f, i+j, dgp);
        dd_edge dgm(f);
        queenInDiagM(f, i-j, dgm);
        // "OR" together the possible attackers for this square
        apply(MAXIMUM, col, dgp, col);
        apply(MAXIMUM, col, dgm, col);
        apply(MAXIMUM, col, qir, col);
        // "AND" with this row
        apply(MULTIPLY, rowcov, col, rowcov);
      } // for j
      printf(", ");
      fflush(stdout);
      // "AND" directly with set of covers
      apply(MULTIPLY, solutions, rowcov, solutions);
    } // for i

    printf("\n");

    // Any solutions?
    if (0==solutions.getNode()) {
      printf("\tNo solutions\n");
      continue;
    } else {
      printf("\tSuccess\n");
      break;
    }
  } // for q

  // Cleanup
  for (int i=0; i<N; i++) delete qic[i];
  for (int i=0; i<2*N-1; i++) {
    delete qidp[i];
    delete qidm[i];
  }
  delete[] qic;
  delete[] qidp;
  delete[] qidm;

  printf("\n%lg seconds CPU time elapsed\n", stopwatch.elapsed());
  printf("Forest stats:\n");
  printf("\t%ld current nodes\n", f->getCurrentNumNodes());
  printf("\t%ld peak nodes\n", f->getPeakNumNodes());
  printf("\t");
  printmem(f->getCurrentMemoryUsed());
  printf(" current memory\n\t");
  printmem(f->getPeakMemoryUsed());
  printf(" peak memory\n");

  long c;
  apply(CARDINALITY, solutions, c);
  printf("\nFor a %dx%d chessboard, ", N, N);
  printf("there are %ld covers with %d queens\n\n", c, q);

  if (!outfile) return 0;


  printf("Writing solutions to file %s\n", argv[2]);
  
  fprintf(outfile, "%d # Board dimension\n\n", N);
  // show the solutions
  dd_edge::const_iterator iter = solutions.begin();
  long counter;
  for (counter = 1; iter; ++iter, ++counter) {
    fprintf(outfile, "solution %5ld:  ", counter);
    const int* minterm = iter.getAssignments();
    for (int i=0; i<N; i++) for (int j=0; j<N; j++) {
      if (minterm[ijmap(i,j)]) {
        fprintf(outfile, "(%2d, %2d) ", i+1, j+1);
      }
    }
    fprintf(outfile, "\n");
  } // for iter
  fprintf(outfile, "\n");
  cleanup();
  return 0;
}

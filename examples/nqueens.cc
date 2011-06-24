
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
    Builds the set of solutions to the N-queens problem for 
    user-specified N.

    In other words, finds all possible ways to put N queens onto
    an NxN chessboard so that no queen is attacking any other queen.
    Clearly, we cannot have two queens in any given row, therefore
    every such solution must have exactly one queen in each row.
    We use an MDD with N variables, with variable x_i indicating
    the column position (from 1 to N) for the queen in row i.
*/

#include <cstdio>

#include "meddly.h"
#include "timer.h"

using namespace MEDDLY;

int N;
int* scratch;
compute_manager* CM;

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


#define LINEAR_INTERSECTIONS
#ifdef  LINEAR_INTERSECTIONS

void intersect(dd_edge** A, int L)
{
  if (0==A[0]) {
    // find first non-zero and swap
    for (int i=1; i<L; i++) {
      if (A[i]) {
        A[0] = A[i];
        A[i] = 0;
        break;
      }
    }
    if (0==A[0]) return;
  }
  fprintf(stderr, "\t");
  for (int i=1; i<L; i++) {
    fprintf(stderr, "%d ", L-i);
    if (A[i]) {
      apply(MULTIPLY, *A[0], *A[i], *A[0]);
      delete A[i];
      A[i] = 0;
    }
  }
  fprintf(stderr, "\n");
}

#else 

void intersect(dd_edge** A, int L)
{
  while (L>1) {
    printf("\t%2d terms to combine ", L);
    // combine adjacent pairs
    for (int i=0; i<L; i+=2) {
      if (A[i] && A[i+1]) {
        apply(MULTIPLY, *A[i], *A[i+1], *A[i]);
        delete A[i+1];
        A[i+1] = 0;
        printf(".");
        fflush(stdout);
      }
    } // for i
    fprintf(stderr, "\n");
    // compact
    int p = 0;
    for (int i=0; i<L; i++) {
      if (A[i]) {
        if (i>p) {
          A[p] = A[i];
          A[i] = 0;
        }
        p++;
      }
    } // for i
    L = p;
  } // while
}

#endif


void createQueenNodes(forest* f, int q, dd_edge &col, dd_edge &cp, dd_edge &cm)
{
  assert(q>0);
  assert(q<=N);
  f->createEdgeForVar(q, false, col);
  for (int i=0; i<N; i++) {
    scratch[i] = i+q;
  }
  f->createEdgeForVar(q, false, scratch, cp);
  for (int i=0; i<N; i++) {
    scratch[i] = i-q;
  }
  f->createEdgeForVar(q, false, scratch, cm);
}

int main()
{
  timer watch;
  initialize();
  CM = getComputeManager();
  assert(CM);
  printf("Using %s\n", getLibraryInfo(0));
  printf("N-Queens solutions.  Enter the value for N:\n");
  scanf("%d", &N);
  if (N<1) return 0;
  scratch = new int[N+1];
  
  watch.note_time();
  printf("Initializing domain and forest\n");

  for (int i=0; i<N; i++) {
    scratch[i] = N;
  }
  domain* d = createDomainBottomUp(scratch, N);
  assert(d);
  forest* f = d->createForest(false, forest::INTEGER, forest::MULTI_TERMINAL);
  assert(f);

  // Set up MDD options
  f->setReductionRule(forest::FULLY_REDUCED);
  f->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE);
  f->setNodeDeletion(forest::PESSIMISTIC_DELETION);

  printf("Building nodes for queen column and diagonals\n");

  dd_edge** col = new dd_edge*[N];
  dd_edge** dgp = new dd_edge*[N];
  dd_edge** dgm = new dd_edge*[N];
  dd_edge** constr = new dd_edge*[N+1];

  for (int i=0; i<N; i++) {
    col[i] = new dd_edge(f);
    dgp[i] = new dd_edge(f);
    dgm[i] = new dd_edge(f);
    createQueenNodes(f, i+1, *col[i], *dgp[i], *dgm[i]);
    constr[i] = new dd_edge(f);
    f->createEdge(int(1), *constr[i]);
  }
  constr[N] = 0;

  for (int i=0; i<N-1; i++) {
    printf("Building queen %2d constraints\n", i+1);
    for (int j=N-1; j>i; j--) {
      dd_edge uniq_col(f);
      apply(NOT_EQUAL, *col[i], *col[j], uniq_col);
      dd_edge uniq_dgp(f);
      apply(NOT_EQUAL, *dgp[i], *dgp[j], uniq_dgp);
      dd_edge uniq_dgm(f);
      apply(NOT_EQUAL, *dgm[i], *dgm[j], uniq_dgm);
      // build overall "not attacking each other" set...
      apply(MULTIPLY, uniq_col, uniq_dgp, uniq_col);
      apply(MULTIPLY, uniq_col, uniq_dgm, uniq_col);
      int k = uniq_col.getLevel()-1;
      if (k<0) k=0;
      assert(k<N);
      apply(MULTIPLY, *constr[k], uniq_col, *constr[k]);
    } // for j
  } // for i
  printf("Building solutions\n");
  intersect(constr, N);
  assert(constr[0]);
  dd_edge solutions = *constr[0];
  watch.note_time();

  printf("Elapsed time: %lf seconds\n", watch.get_last_interval()/1000000.0);

  printf("Cleanup\n");
  for (int i=0; i<N; i++) {
    delete col[i];
    delete dgp[i];
    delete dgm[i];
    delete constr[i];
  }
  delete[] col;
  delete[] dgp;
  delete[] dgm;
  delete[] constr;

  printf("Done.\n\n");

  printf("Set of solutions requires %d nodes\n", solutions.getNodeCount());
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
  printf("\nThere are %ld solutions to the %d-queens problem\n\n", c, N);

  // show one of the solutions
  dd_edge::const_iterator first = solutions.begin();
  if (first) {
    const int* minterm = first.getAssignments();
    printf("One solution:\n");
    for (int i=1; i<=N; i++) {
      printf("\tQueen for row %2d in column %2d\n", i, minterm[i]+1);
    }
  }
  destroyDomain(d);
  cleanup();
  return 0;
}

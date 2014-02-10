
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
    Program to analyze the "Square 1" puzzle.

    Organized as follows.
    Top has 12 slots, but some pieces take up two consecutive slots.
    Bottom is the same way.
    We can do an "exchange" twist if there are no double-wide pieces
    in the critical slots.
*/

#include "meddly.h"
#include "meddly_expert.h"
#include "timer.h"

// #define DEBUG_EVENTS
// #define DUMP_NSF

using namespace MEDDLY;

/// "Global" array of minterms,
/// used by Equals and NotEquals

int** minterms;
int** minprime;
int mtlength;

void InitMinterms()
{
  minterms = 0;
  minprime = 0;
  mtlength = 0;
}

void ExpandMinterms(int n, int K)
{
  if (n<mtlength) return;
  minterms = (int**) realloc(minterms, n*sizeof(int*));
  minprime = (int**) realloc(minprime, n*sizeof(int*));
  if (0==minterms || 0==minprime) {
    fprintf(stderr, "Couldn't resize minterms %d\n", n);
    exit(1);
  }
  for (; mtlength < n; mtlength++) {
    minterms[mtlength] = new int[K+1];
    minprime[mtlength] = new int[K+1];
  }
}

/*
  Builds the "constraint" va == vb,
  for the given variables (levels).
  Use negative to denote primed.
*/
void Equals(int va, int vb, dd_edge &answer)
{
  expert_forest* EF = (expert_forest*) answer.getForest();
  int K = EF->getNumVariables();
  int asz = EF->getLevelSize(va);
  int bsz = EF->getLevelSize(vb);
  int size = MIN(asz, bsz);

  ExpandMinterms(size, K);

  for (int i=0; i<size; i++) {
    for (int k=1; k<=K; k++) {
      minterms[i][k] = DONT_CARE;
      minprime[i][k] = DONT_CARE;
    }
    if (va<0) minprime[i][-va] = i;
    else      minterms[i][va]  = i;
    if (vb<0) minprime[i][-vb] = i;
    else      minterms[i][vb]  = i;
  }

  EF->createEdge(minterms, minprime, size, answer);
}

/*
  Builds the "constraint" va != vb,
  for the given variables (levels).
  Use negative to denote primed.
*/
void NotEquals(int va, int vb, dd_edge &answer)
{
  expert_forest* EF = (expert_forest*) answer.getForest();
  int K = EF->getNumVariables();
  int asz = EF->getLevelSize(va);
  int bsz = EF->getLevelSize(vb);

  ExpandMinterms(asz * bsz, K);

  int mt = 0;
  for (int ai=0; ai<asz; ai++) {
    for (int bi=0; bi<bsz; bi++) if (ai != bi) {
      for (int k=1; k<=K; k++) {
        minterms[mt][k] = DONT_CARE;
        minprime[mt][k] = DONT_CARE;
      }
      if (va<0) minprime[mt][-va] = ai;
      else      minterms[mt][ va] = ai;
      if (vb<0) minprime[mt][-vb] = bi;
      else      minterms[mt][ vb] = bi;
      mt++;
    }
  }

  EF->createEdge(minterms, minprime, mt, answer);
}

/*
  "Rotates" a set of variables 
    x0' == x0+r mod n
    x1' == x1+r mod n
    x2' == x2+r mod n
    ...
    xn-1' == xn-1+r mod n
*/
void Rotate(const int* vars, int nv, int r, dd_edge &answer)
{
  expert_forest* EF = (expert_forest*) answer.getForest();

  int i=0;
  int j = (i+r) % nv;
  Equals(-vars[i], vars[j], answer);

  for (i=1; i<nv; i++) {
    j = (i+r) % nv;
    dd_edge temp(EF);
    Equals(-vars[i], vars[j], temp);
    answer *= temp;
  }
}

/*
  Builds next-state function for square1
*/
void Square1NSF(const int* top, const int* bottom, const int exch, dd_edge &answer)
{
  dd_edge topl(answer), topr(answer), botl(answer), botr(answer);


  printf("Building event: rotate top by %d\n", 1);
  fflush(stdout);
  Rotate(top, 12, 1, topl);
#ifdef DEBUG_EVENTS
  event.show(stdout, 2);
  fflush(stdout);
#endif

  printf("Building event: rotate top by %d\n", -1);
  fflush(stdout);
  Rotate(top, 12, 11, topr);
#ifdef DEBUG_EVENTS
  event.show(stdout, 2);
  fflush(stdout);
#endif

  printf("Building event: rotate bottom by %d\n", 1);
  fflush(stdout);
  Rotate(bottom, 12, 1, botl);
#ifdef DEBUG_EVENTS
  event.show(stdout, 2);
  fflush(stdout);
#endif

  printf("Building event: rotate bottom by %d\n", -1);
  fflush(stdout);
  Rotate(bottom, 12, 11, botr);
#ifdef DEBUG_EVENTS
  event.show(stdout, 2);
  fflush(stdout);
#endif

  printf("Accumulating events...");
  fflush(stdout);

  answer = topl + topr + botl + botr;

  printf("\n");
  fflush(stdout);

#ifdef DUMP_NSF
  printf("Overall NSF:\n");
  answer.show(stdout, 2);
  fflush(stdout);
#endif
}

void printStats(const char* who, const forest* f)
{
  printf("%s stats:\n", who);
  const expert_forest* ef = (expert_forest*) f;
  ef->reportStats(stdout, "\t",
    expert_forest::HUMAN_READABLE_MEMORY  |
    expert_forest::BASIC_STATS | expert_forest::EXTRA_STATS |
    expert_forest::STORAGE_STATS | expert_forest::HOLE_MANAGER_STATS
  );
}


int main()
{
  MEDDLY::initialize();

  // Variables for top piees
  const int top[] = { 25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5, 3 };
  // Variables for bottom pieces
  const int bottom[] = { 24, 22, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2 };
  // Variable for exchanged or not
  const int exchange = 1;

  //
  // Build the domain
  //
  int scratch[25];
  for (int i=0; i<25; i++) scratch[i] = 0;
  for (int i=0; i<12; i++) {
    // set top and bottom sizes
    if (scratch[top[i]-1]) {
      fprintf(stderr, "Error - variable %d duplicated\n", top[i]);
      return 1;
    }
    scratch[top[i]-1] = 16;   // There's 16 top/bottom pieces, total
    if (scratch[bottom[i]-1]) {
      fprintf(stderr, "Error - variable %d duplicated\n", bottom[i]);
      return 1;
    }
    scratch[bottom[i]-1] = 16;
  }
  if (scratch[exchange-1]) {
    fprintf(stderr, "Error - variable %d duplicated\n", exchange);
    return 1;
  }
  scratch[exchange-1] = 2;  // exchanged or not
  domain* D = createDomainBottomUp(scratch, 25);
 
  //
  // Build NSF for possible "1-step" moves
  //
  forest::policies p(true);
  // p.setPessimistic();
  // p.setFullyReduced();
  forest* F = D->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL, p);
  dd_edge nsf(F);
  Square1NSF(top, bottom, exchange, nsf);
  printStats("MxD", F);

  MEDDLY::cleanup();
  return 0;
}

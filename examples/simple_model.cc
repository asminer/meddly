
// $Id$

/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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

#include "meddly.h"
#include "simple_model.h"

#define VERBOSE

//#define DEBUG_GENERATE
//#define DEBUG_EVENTS
//#define DEBUG_EVENTS_LONG

inline int MAX(int a, int b) {
  return (a>b) ? a : b;
}

void buildNextStateFunction(const char* const* events, int nEvents,
  MEDDLY::forest* mxd, MEDDLY::dd_edge &nsf, int verb)
{
  using namespace MEDDLY;
  if (verb) fprintf(stderr, "Building next-state function\n");

  // set up auxiliary mtmxd forest and edges
  domain* d = mxd->useDomain();
  forest* mtmxd = d->createForest(
    true, forest::INTEGER, forest::MULTI_TERMINAL
  );
  int nVars = d->getNumVariables();
  int maxBound = d->getVariableBound(1, false);
  for (int i=2; i<=nVars; i++) {
    maxBound = MAX(maxBound, d->getVariableBound(i, false));
  }
  maxBound++;
  int* temp = new int[maxBound];
  int* minterm = new int[nVars+1];
  int* mtprime = new int[nVars+1];
  dd_edge** varP  = new dd_edge*[nVars+1];
  varP[0] = 0;
  dd_edge** inc   = new dd_edge*[nVars+1];
  inc[0] = 0;
  dd_edge** dec   = new dd_edge*[nVars+1];
  dec[0] = 0;

  //  Create edge for each variable xi'
  for (int i=1; i<=nVars; i++) {
    varP[i] = new dd_edge(mtmxd);
    mtmxd->createEdgeForVar(i, true, varP[i][0]);
  }

  // Create edge for each function xi+1
  for (int i=0; i<maxBound; i++) temp[i] = i+1;
  for (int i=1; i<=nVars; i++) {
    inc[i] = new dd_edge(mtmxd);
    mtmxd->createEdgeForVar(i, false, temp, inc[i][0]);
  }

  // Create edge for each function xi-1
  for (int i=0; i<maxBound; i++) temp[i] = i-1;
  for (int i=1; i<=nVars; i++) {
    dec[i] = new dd_edge(mtmxd);
    mtmxd->createEdgeForVar(i, false, temp, dec[i][0]);
  }

  mxd->createEdge(false, nsf);

  for (int e=0; e<nEvents; e++) {
    const char* ev = events[e];
    if (2==verb) fprintf(stderr, "%5d", e);
    if (verb>2) fprintf(stderr, "Event %5d", e);

    dd_edge nsf_ev(mxd);
    dd_edge term(mxd);

    //
    // build mask for this event
    //
    for (int i=1; i<=nVars; i++) {
      if ('.' == ev[i]) {
        minterm[i] = DONT_CARE;
        mtprime[i] = DONT_CHANGE;
      } else {
        minterm[i] = DONT_CARE;
        mtprime[i] = DONT_CARE;
      }
    }
    mxd->createEdge(&minterm, &mtprime, 1, nsf_ev);
#ifdef DEBUG_EVENTS
    printf("Initial nsf for event %d\n", e);
    nsf_ev.show(stdout, 2);
#endif
    if (verb>2) fprintf(stderr, " : ");
    
    //
    // 'and' with the "do care" levels
    //
    for (int i=1; i<=nVars; i++) {
      dd_edge docare(mtmxd);

      if ('.' == ev[i]) {
        if (verb>3) fputc('.', stderr);
        continue;
      } else {
        if (verb>2) fputc(ev[i], stderr);
      }
      switch (ev[i]) {
        case '+':   apply(EQUAL, varP[i][0], inc[i][0], docare);
                    break;

        case '-':   apply(EQUAL, varP[i][0], dec[i][0], docare);
                    break;

        default:    throw 1;
      } // switch
      apply(COPY, docare, term);
#ifdef DEBUG_EVENTS
      printf("Term for event %d, level %d\n", e, i);
      term.show(stdout, 2);
#endif
      nsf_ev *= term;
    } // for i

    //
    //  union with overall
    //
#ifdef DEBUG_EVENTS
    printf("Complete nsf for event %d:\n", e);
    nsf_ev.show(stdout, 2);
#endif
    if (verb>2) fputc(' ', stderr);
    nsf += nsf_ev;
#ifdef DEBUG_EVENTS_LONG
    printf("Complete after adding event %d:\n", e);
    nsf.show(stdout, 2);
#endif

    if (verb>2) fputc('\n', stderr);
  } // for e
  if (verb==2) fputc('\n', stderr);

  // cleanup
  delete[] mtprime;
  delete[] minterm;
  delete[] temp;
  for (int i=1; i<=nVars; i++) {
    delete varP[i];
    delete inc[i];
    delete dec[i];
  }
  delete[] varP;
  delete[] inc;
  delete[] dec;
  destroyForest(mtmxd);

#ifdef DEBUG_EVENTS
  printf("Complete NSF:\n");
  nsf.show(stdout, 2); 
#endif
}




//
//  Explicit RS construction
//

bool fireEvent(const char* event, const int* current, int* next, int nVars)
{
  for (int i=nVars; i; i--) {
    if ('.' == event[i]) {
      next[i] = current[i];
      continue;
    }
    if ('-' == event[i]) {
      next[i] = current[i] - 1;
      if (next[i] < 0) return false;
      continue;
    }
    if ('+' == event[i]) {
      next[i] = current[i] + 1;
      // TBD ... check for overflow
      continue;
    }
    throw 1;  // bad event string
  }
  return true;
}

void explicitReachset(const char* const* events, int nEvents, 
  MEDDLY::forest* f, MEDDLY::dd_edge &expl, MEDDLY::dd_edge &RS, int batchsize)
{
  int nVars = f->getDomain()->getNumVariables();
  if (batchsize < 1) batchsize = 256;

  // initialize batch memory
  int** minterms = new int*[batchsize];
  for (int b=0; b<batchsize; b++) {
    minterms[b] = new int[1+nVars];
  }
  
  // unexplored states
  MEDDLY::dd_edge unexplored(f);
  // batch of states
  MEDDLY::dd_edge batch(f);
  int b = 0; 
  // exploration loop.
  for (;;) {
    unexplored.clear();
    MEDDLY::enumerator I(expl);
    if (!I) break;    // nothing left to explore, bail out
    // explore everything in expl
    for (; I; ++I) {
      const int* curr = I.getAssignments();
#ifdef DEBUG_GENERATE
      printf("Exploring state: (%d", curr[1]);
      for (int n=2; n<=nVars; n++) printf(", %d", curr[n]);
      printf(")\n");
#endif
      // what's enabled?
      for (int e=0; e<nEvents; e++) {
        if (!fireEvent(events[e], curr, minterms[b], nVars)) continue;
#ifdef DEBUG_GENERATE
        printf("  -- (event %d) --> (%d", e, minterms[b][1]);
        for (int n=2; n<=nVars; n++) printf(", %d", minterms[b][n]);
        printf(")\n");
#endif
        bool seen;
        f->evaluate(RS, minterms[b], seen);
        if (seen) continue;     // already known in RS
        f->evaluate(unexplored, minterms[b], seen);
        if (seen) continue;     // already in unexplored list
        b++;
        if (b>=batchsize) {
          // Buffer is full; flush it
          f->createEdge(minterms, b, batch);
          unexplored += batch;
          RS += batch;
          b = 0;
        }
      }
    }
    // expl is empty.
    // Flush the buffer
    if (b) {
      f->createEdge(minterms, b, batch);
      unexplored += batch;
      b = 0;
    }
    RS += unexplored;
    expl = unexplored;
  }

  // cleanup batch memory
  for (int b=0; b<batchsize; b++) {
    delete[] minterms[b];
  }
  delete[] minterms;
}



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

#include "../src/meddly.h"
#include "simple_model.h"

#define VERBOSE
#define NOT_KNOWN     -2
#define OUT_OF_BOUNDS -1



// #define DEBUG_GENERATE
// #define DEBUG_EVENTS
// #define DEBUG_EVENTS_LONG

// #define SHOW_EVENT_HANDLES

// #define SAME_FOREST_OPERATIONS

inline int MAX(int a, int b) {
  return (a>b) ? a : b;
}

/*
    Common core for building the next state function.

      @param  events    Array of dimension \a nEvents.
                        Each entry is a string where
                        str[j] tells how this event affects variable j:
                        '+' means the event increases the variable by one
                            (unless this would violate the bound)
                        '-' means the event decreases the variable by one
                            (unless the value is already zero)
                        '.' means no change.
      @param  nEvents   Size of events array.
      @param  mxd       The forest to use; should be a boolean one for relations.

      @param  pnsf      If not null, store the relation for each event here.

      @param  mono      If not null, accumulate the monolithic relation here.

      @param  verb      Verbosity level.
*/

void buildNextStateFunction(const char* const* events, int nEvents,
  MEDDLY::forest* mxd,
  MEDDLY::satpregen_opname::pregen_relation* pnsf,
  MEDDLY::dd_edge* mono, int verb)
{
  using namespace MEDDLY;
  if (verb) fprintf(stderr, "Building next-state function\n");

  // set up auxiliary mtmxd forest and edges
  domain* d = mxd->useDomain();
  forest* mtmxd = forest::create(d,
    true, range_type::INTEGER, edge_labeling::MULTI_TERMINAL
  );
  int nVars = d->getNumVariables();
  int maxBound = d->getVariableBound(1, false);
  for (int i=2; i<=nVars; i++) {
    maxBound = MAX(maxBound, d->getVariableBound(i, false));
  }
  maxBound++;
  long* temp = new long[maxBound];
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

  //
  // Initialize accumulators
  //
  if (mono) mxd->createEdge(false, *mono);

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
#ifdef SAME_FOREST_OPERATIONS
      dd_edge docare(mtmxd);
#endif

      if ('.' == ev[i]) {
        if (verb>3) fputc('.', stderr);
        continue;
      } else {
        if (verb>2) fputc(ev[i], stderr);
      }
      switch (ev[i]) {
        case '+':
#ifdef SAME_FOREST_OPERATIONS
                    apply(EQUAL, varP[i][0], inc[i][0], docare);
#else
                    apply(EQUAL, varP[i][0], inc[i][0], term);
#endif
                    break;

        case '-':
#ifdef SAME_FOREST_OPERATIONS
                    apply(EQUAL, varP[i][0], dec[i][0], docare);
#else
                    apply(EQUAL, varP[i][0], dec[i][0], term);
#endif
                    break;

        default:    throw 1;
      } // switch
#ifdef SAME_FOREST_OPERATIONS
      apply(COPY, docare, term);
#endif
#ifdef DEBUG_EVENTS
      printf("Term for event %d, level %d\n", e, i);
      term.show(stdout, 2);
#endif
      nsf_ev *= term;
    } // for i

#ifdef DEBUG_EVENTS
    printf("Complete nsf for event %d:\n", e);
    nsf_ev.show(stdout, 2);
#endif
    if (verb>2) fputc(' ', stderr);
#ifdef SHOW_EVENT_HANDLES
    fprintf(stderr, "%d", nsf_ev.getNode());
#endif

    if (pnsf) {
      //
      //  Add event to relation
      //
      pnsf->addToRelation(nsf_ev);
    }

    if (mono) {
      //
      //  union with overall
      //
      *mono += nsf_ev;
#ifdef DEBUG_EVENTS_LONG
      printf("Complete after adding event %d:\n", e);
      mono->show(stdout, 2);
#endif
    }

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
  mono->show(stdout, 2);
#endif
}


void buildNextStateFunction(const char* const* events, int nEvents,
  MEDDLY::forest* mxd, MEDDLY::dd_edge &nsf, int verb)
{
  buildNextStateFunction(events, nEvents, mxd, 0, &nsf, verb);
}

void buildNextStateFunction(const char* const* events, int nEvents,
  MEDDLY::satpregen_opname::pregen_relation* pnsf, int verb)
{
  if (0==pnsf) return;
  buildNextStateFunction(events, nEvents, pnsf->getRelForest(), pnsf, 0, verb);
  pnsf->finalize();

#ifdef SHOW_EVENT_HANDLES
  using namespace MEDDLY;
  // check what we got
  for (int k=16; k; k--) {
    int len = pnsf->lengthForLevel(k);
    if (0==len) continue;
    printf("Events at level %d:\n\t", k);
    node_handle* List = pnsf->arrayForLevel(k);
    for (int i=0; i<len; i++)
      printf("%d ", List[i]);
    printf("\n");
  }
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
  int b;
  int nVars = f->getDomain()->getNumVariables();
  if (batchsize < 1) batchsize = 256;

  // initialize batch memory
  int** minterms = new int*[batchsize];
  for (b=0; b<batchsize; b++) {
    minterms[b] = new int[1+nVars];
  }

  // unexplored states
  MEDDLY::dd_edge unexplored(f);
  // batch of states
  MEDDLY::dd_edge batch(f);
  b = 0;
  // exploration loop.
  MEDDLY::enumerator I(expl);
  for (;;) {
    f->createEdge(false, unexplored);
    I.start(expl);
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
  for (b=0; b<batchsize; b++) {
    delete[] minterms[b];
  }
  delete[] minterms;
}


/*************************************************************/
int* nxtList;
class derRelNode : public MEDDLY::relation_node
{
  int en;
  int fr;
  int h;
public:
  derRelNode(MEDDLY::forest* mxdF, int level, int dwn, int e, int f, int inh): relation_node(141010, mxdF, level, dwn, e, f, inh)
  {
    en = e;
    fr = f;
    h = inh;
  }

  long nextOf(long i) override
  {

    if(i>=getPieceSize()) { //Array needs to be allocated
       expandTokenUpdate(i);
    }

    if(getTokenUpdate()[i]==NOT_KNOWN) //Array needs to be updated
    {
    	long result = ((i>=en) && (h==-1?true:i<h))? i+fr : OUT_OF_BOUNDS;
    	long val = result;
     	setTokenUpdateAtIndex(i,val);
    }
   return getTokenUpdate()[i];
  }
};


void buildImplicitRelation(const int* const* events, int nEvents,int nPlaces, int bounds, MEDDLY::forest* mddF, MEDDLY::forest* mxdF, MEDDLY::satimpl_opname::implicit_relation* T)
{

  unsigned node_count = 0;
  int* tops_of_events = (int*)malloc(size_t(nEvents)*sizeof(int));

  for(int e = 0;e < nEvents; e++)
    {
    bool done = false;
    for( int p = 1; p <= nPlaces; p++)
      {
      if(events[e][p]!='0') node_count +=1;
      if((events[e][nPlaces-p+1]!=0)&&(!done))
        {
        tops_of_events[e] = nPlaces-p+1;
        done = true;
        }
      }
    }


  derRelNode** rNode = (derRelNode**)malloc(node_count*sizeof(derRelNode*));

  // Add/Subtract Tokens
  nxtList = (int*)malloc((node_count+2)*sizeof(int));
  nxtList[0] = 0;
  nxtList[1] = 0;

  int rctr = 0;
  for(int e = 0;e < nEvents; e++)
    {
    unsigned long sign = 0;
    int previous_node_handle = 1;
    for( int p = 1; p <= nPlaces; p++)
      {
       sign = events[e][p]>=0?(sign*10)+events[e][p]:(sign*100)+events[e][p];
        if(events[e][p]!=0)
        {
          rNode[rctr] = new derRelNode(mxdF,p,previous_node_handle,events[e][p]<0?-events[e][p]:0,events[e][p],-1);
          previous_node_handle = T->registerNode((tops_of_events[e]==p),rNode[rctr]);

          rctr++;
          nxtList[previous_node_handle] = events[e][p];
        }
      }
    }
}


MEDDLY::sathyb_opname::event** buildHybridRelation(const int* const* events, int nEvents,int nPlaces, int bounds, MEDDLY::forest* mddF, MEDDLY::forest* mxdF) {

  unsigned node_count = 0;
  int* tops_of_events = (int*)malloc(size_t(nEvents)*sizeof(int));

  for(int e = 0;e < nEvents; e++)
    {
    bool done = false;
    for( int p = 1; p <= nPlaces; p++)
      {
      if(events[e][p]!='0') node_count +=1;
      if((events[e][nPlaces-p+1]!=0)&&(!done))
        {
        tops_of_events[e] = nPlaces-p+1;
        done = true;
        }
      }
    }

  derRelNode** rNode = (derRelNode**)malloc(node_count*sizeof(derRelNode*));
  MEDDLY::sathyb_opname::event** T = (MEDDLY::sathyb_opname::event**)malloc(nEvents*sizeof(MEDDLY::sathyb_opname::event*));

  // Add/Subtract Tokens
  nxtList = (int*)malloc((node_count+2)*sizeof(int));
  nxtList[0] = 0;
  nxtList[1] = 0;

  int rctr = 0;
  for(int e = 0;e < nEvents; e++)
    {
    unsigned long sign = 0;
    int previous_node_handle = 1;
    int e_rn = 0;
    for( int p = 1; p <= nPlaces; p++)
      {
       sign = events[e][p]>=0?(sign*10)+events[e][p]:(sign*100)+events[e][p];
        if(events[e][p]!=0)
        {
          e_rn++;
          rNode[rctr] = new derRelNode(mxdF,p,-1,events[e][p]<0?-events[e][p]:0,events[e][p],-1);
          rctr++;
          nxtList[previous_node_handle] = events[e][p];
        }
      }
       T[e] = new MEDDLY::sathyb_opname::event(NULL, 0, (MEDDLY::relation_node**)(&rNode[rctr-e_rn]), e_rn);
    }

    return T;

}



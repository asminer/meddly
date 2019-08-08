
#include "../src/meddly.h"
#include "simple_model.h"

#define VERBOSE

using namespace MEDDLY;

const char* mdl[] = { "X-+." };

int main()
{
  MEDDLY::initialize();

  // Initialize domain
  int* tmp = new int[4];
  int_extra* gtmp = new  int_extra[4];

  tmp[0] = 3;
  tmp[1] = 3;
  tmp[2] = 3;

  gtmp[0] = tmp[0];
  gtmp[1] = tmp[1];
  gtmp[2] = tmp[2];

  domain* d = createDomainBottomUp(variable::variableTypes::boundedClass,tmp, 3);

  // Initialize forests
  forest* mdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL);
  forest* mxd = d->createForest(1, forest::BOOLEAN, forest::MULTI_TERMINAL);

  // build initial state
  tmp[1] = 2;
  tmp[2] = 0;
  tmp[3] = 0;
  dd_edge init_state(mdd);
  gtmp[1] = tmp[1];
  gtmp[2] = tmp[2];
  gtmp[3] = tmp[3];
  mdd->createEdge(&gtmp, 1, init_state);

  // build next-state function
  dd_edge nsf(mxd);
  buildNextStateFunction(mdl, 1, mxd, nsf);
  tmp[3] = 0;
  tmp[2] = -1;
  tmp[1] = -1;

  gtmp[3] = tmp[3];
  gtmp[2] = tmp[2];
  gtmp[1] = tmp[1];
  dd_edge mask(mxd);
  mxd->createEdge(&gtmp, &gtmp, 1, mask);
  nsf *= mask;

  // build rs using traditional & saturation
  dd_edge reachable1(mdd);
  dd_edge reachable2(mdd);
  apply(REACHABLE_STATES_DFS, init_state, nsf, reachable1);
  apply(REACHABLE_STATES_BFS, init_state, nsf, reachable2);

  // Display everything
#ifdef VERBOSE
  FILE_output out(stdout);
  printf("Initial state:\n");
  init_state.show(out, 2);

  printf("NSF:\n");
  nsf.show(out, 2);

  printf("DFS states:\n");
  reachable1.show(out, 2);

  printf("BFS states:\n");
  reachable2.show(out, 2);
#endif

  int retval;
  if (reachable1 == reachable2) {
    retval = 0;
  } else {
    retval = 1;
  }

#ifdef VERBOSE
  if (retval) {
    printf("\nReachable states DO NOT match\n\n");
  } else {
    printf("\nReachable states match\n\n");
  }
#endif

  // cleanup
  MEDDLY::cleanup();
  delete[] tmp;
  return retval;
}


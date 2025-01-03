
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



/**
 * rebuild_mdd.cc
 *
 * Rebuild a mdd forest with a different variable ordr.
 */

#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <cassert>
#include "../src/meddly.h"
#include "reorder.h"

using namespace MEDDLY;

// Only include this file if referring to operations that are not
// available via the mddlib interface.
//
// Also, include this file if you would like to create a custom
// operation by deriving one of the existing operations.
// #include "operation_ext.h"

// Timer class
#include "../timing/timer.h"

#define USE_REALS 0

#if USE_REALS
  typedef float element_type;
#else
  typedef long element_type;
#endif


// verbose: 0: minimum, 2: maximum
const int verbose = 0;


void printStats(const char* who, const forest* f)
{
  printf("%s stats:\n", who);
  FILE_output mout(stdout);
  f->reportStats(mout, "\t",
    HUMAN_READABLE_MEMORY  |
    BASIC_STATS | EXTRA_STATS |
    STORAGE_STATS |
    HOLE_MANAGER_STATS | HOLE_MANAGER_DETAILED
  );
}



void printUsage(FILE *outputStream)
{
  fprintf(outputStream,
      "Usage: test_union_mtmdd <#Variables> <VariableBound> <#Elements>\
 [#terms]\n");
}


void reorderVariablesByRebuilding(dd_edge &e)
{
	printf("Initial: %ld\n", e.getForest()->getCurrentNumNodes());

	forest* f = e.getForest();
	int num_var = f->getNumVariables();
	int* level2var = new int[num_var + 1];
	level2var[0] = 0;
	for(int i=1; i<=num_var; i++) {
		level2var[i] = i;
	}
	shuffle(level2var, 1, num_var);

	forest* target = forest::create(f->getDomain(), f->isForRelations(),
			f->getRangeType(), f->getEdgeLabeling() , f->getPolicies());
	target->reorderVariables(level2var);
	delete[] level2var;

	dd_edge e2(target);

	{
	  global_rebuilder gr(f, target);

	  timer start;
	  start.note_time();
	  e2 = gr.rebuild(e);
	  start.note_time();
	  printf("Time interval: %.4e seconds\n", start.get_last_seconds());

	  printf("Source: %lu\n", e.getNodeCount());
	  printf("Final: %lu\n", e2.getNodeCount());
	}

	f->createConstant(0L, e);
	printf("Source Discarded: %lu\n", e.getNodeCount());

	{
	  global_rebuilder gr(target, f);

	  timer start;
	  start.note_time();
	  e = gr.rebuild(e2);
	  start.note_time();
	  printf("Time interval: %.4e seconds\n", start.get_last_seconds());

	  printf("Final: %lu\n", e2.getNodeCount());
	  printf("Source: %lu\n", e.getNodeCount());
	}
}

int main(int argc, char *argv[])
{
  if (argc != 4 && argc != 5) {
    printUsage(stdout);
    exit(1);
  }

  srand(time(NULL));

  // initialize number of variables, their bounds and the number of elements
  // to create

  int nVariables = 0;
  int variableBound = 0;
  int nElements = 0;
  int nTerms = 0;

  sscanf(argv[1], "%d", &nVariables);
  assert(nVariables > 0);

  sscanf(argv[2], "%d", &variableBound);
  assert(variableBound > 0);

  sscanf(argv[3], "%d", &nElements);
  assert(nElements > 0);

  if (argc == 5) {
    sscanf(argv[4], "%d", &nTerms);
    assert(nTerms > 0);
  } else {
    nTerms = variableBound;
  }

  printf("#variables: %d, variable bound: %d, #elements: %d\n",
      nVariables, variableBound, nElements);

  // initialize the variable bounds array to provide to the domain

  int* bounds = new int[nVariables];
  assert(bounds != 0);
  for (int i = 0; i < nVariables; ++i)
  {
    bounds[i] = variableBound;
  }

  initialize();

  // Create a domain
  domain *d = domain::createBottomUp(bounds, nVariables);
  assert(d != 0);

  // Create a MTMDD forest in this domain
  policies p(false);
  p.setPessimistic();
#if USE_REALS
  forest* mtmdd =
    forest::create(d, SET, range_type::REAL, edge_labeling::MULTI_TERMINAL, p);
#else
  forest* mtmdd =
    forest::create(d, SET, range_type::INTEGER, edge_labeling::MULTI_TERMINAL, p);
#endif
  assert(mtmdd != 0);

  // create a random minterm collection
  minterm_coll elements(nElements, mtmdd);
  while (elements.size() < nElements) {
      for (unsigned j=nVariables; j; --j) {
          const int v = int(variableBound * rand() / (RAND_MAX+1.0));
          elements.unused().setVar(j, v);
      }
#if USE_REALS
      elements.unused().setValue(nTerms * rand() / (RAND_MAX + 1.0));
#else
      elements.unused().setValue(int(nTerms * rand() / (RAND_MAX + 1.0)));
#endif
      elements.pushUnused();
  }

  // print elements
  if (verbose > 0) {
    FILE_output meddlyout(stdout);
    elements.show(meddlyout);
  }

  dd_edge result(mtmdd);
  elements.buildFunction(result);
  reorderVariablesByRebuilding(result);


  // Cleanup; in this case simply delete the domain
  domain::destroy(d);

  delete[] bounds;

  cleanup();

  return 0;
}



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



/**
 * test_evmdd.cc
 *
 * Testing EV+MDD operations for integers.
 * Operations: min, max, +, -, *, /.
 */

#include <iostream>

// meddly_expert.h include meddly.h
#include "../include/meddly_expert.h"

// Only include this file if referring to operations that are not
// available via the mddlib interface. 
// 
// Also, include this file if you would like to create a custom
// operation by deriving one of the existing operations.
#include "../src/operation_ext.h"

// Timer class
#include "../src/timer.h"

#define USE_REALS 1

#if USE_REALS
  typedef float element_type;
#else
  typedef int element_type;
#endif


// verbose: 0: minimum, 2: maximum
const int verbose = 1;

// Given a forest and an op_code returns the corresponding op_info.
// 
// This is only valid for operations of the form C = A op B,
// where A, B and C belong to the same forest.
op_info* getOp(forest* f, compute_manager::op_code op)
{
  static const int nForests = 3;
  static forest* forests[nForests];
  static expert_compute_manager* ecm = 
    static_cast<expert_compute_manager*>(MEDDLY_getComputeManager());
  assert(ecm != 0);
  assert(f != 0);

  forests[0] = f;
  forests[1] = f;
  forests[2] = f;
  return ecm->getOpInfo(op, forests, nForests);
}


// Given a forest and an instance of an operation returns
// the corresponding op_info.
// 
// This is only valid for operations of the form C = A op B,
// where A, B and C belong to the same forest.
op_info* getOp(forest* f, operation* op)
{
  static const int nForests = 3;
  static forest* forests[nForests];
  static expert_compute_manager* ecm = 
    static_cast<expert_compute_manager*>(MEDDLY_getComputeManager());
  assert(ecm != 0);
  assert(f != 0);
  assert(op != 0);

  forests[0] = f;
  forests[1] = f;
  forests[2] = f;
  return ecm->getOpInfo(op, forests, nForests);
}


// Tests a evmdd operation on the elements provided.
// This function assumes that each element[i] represents
// an element in the given MTMDD.
dd_edge test_evmdd(forest* evmdd, compute_manager::op_code opCode,
    int** element, element_type* terms, int nElements)
{
  // A = first nElements/2 elements combined using +.
  // B = second nElements/2 elements combined using +.
  // C = A op B

  static expert_compute_manager* ecm = 
    static_cast<expert_compute_manager*>(MEDDLY_getComputeManager());
  assert(ecm != 0);

  dd_edge A(evmdd);
  dd_edge B(evmdd);
  dd_edge C(evmdd);

  int half = nElements/2;

  assert(forest::SUCCESS ==
      evmdd->createEdge(element, terms, half, A));
  assert(forest::SUCCESS ==
      evmdd->createEdge(element + half,
        terms + half, nElements - half, B));

  op_info* op = getOp(evmdd, opCode);
  assert(op != NULL);
  assert(compute_manager::SUCCESS == ecm->apply(op, A, B, C));

  if (verbose > 0) {
    printf("A: ");
    A.show(stdout, 2);
    printf("\n\nB: ");
    B.show(stdout, 2);
    printf("\n\nC: ");
    C.show(stdout, 2);
    printf("\n\n");
  }

  return C;
}


void printUsage(FILE *outputStream)
{
  fprintf(outputStream,
      "Usage: test_evmdd <#Variables> <VariableBound> <#Elements>\n");
}


int main(int argc, char *argv[])
{
  if (argc != 4) {
    printUsage(stdout);
    exit(1);
  }

  srandom(1u);

  // initialize number of variables, their bounds and the number of elements
  // to create

  int nVariables = 0;
  int variableBound = 0;
  int nElements = 0;

  sscanf(argv[1], "%d", &nVariables);
  assert(nVariables > 0);

  sscanf(argv[2], "%d", &variableBound);
  assert(variableBound > 0);

  sscanf(argv[3], "%d", &nElements);
  assert(nElements > 0);

  printf("#variables: %d, variable bound: %d, #elements: %d\n",
      nVariables, variableBound, nElements);

  // create the elements randomly
  int** element = (int **) malloc(nElements * sizeof(int *));
  element_type* terms =
    (element_type *) malloc(nElements * sizeof(element_type));

  for (int i = 0; i < nElements; ++i)
  {
    element[i] = (int *) malloc((nVariables + 1) * sizeof(int));
    element[i][0] = 0;
    for (int j = nVariables; j > 0; --j)
    {
      element[i][j] = int(float(variableBound) * random() / (RAND_MAX + 1.0));
      assert(element[i][j] >= 0 && element[i][j] < variableBound);
    }
    terms[i] =
      element_type(float(variableBound) * random() / (RAND_MAX + 1.0));
  }

  // initialize the variable bounds array to provide to the domain

  int* bounds = (int *) malloc(nVariables * sizeof(int));
  assert(bounds != 0);
  for (int i = 0; i < nVariables; ++i)
  {
    bounds[i] = variableBound;
  }

  // Create a domain
  domain *d = MEDDLY_createDomain();
  assert(d != 0);
  assert(domain::SUCCESS == d->createVariablesBottomUp(bounds, nVariables));

  // Create a MTMDD forest in this domain
#if USE_REALS
  forest* evmdd = d->createForest(false, forest::REAL, forest::EVTIMES);
#else
  forest* evmdd = d->createForest(false, forest::INTEGER, forest::EVPLUS);
#endif
  assert(evmdd != 0);

  // print elements
  for (int i = 0; i < nElements; ++i)
  {
    printf("Element %d: [%d", i, element[i][0]);
    for (int j = 1; j <= nVariables; ++j)
    {
      printf(" %d", element[i][j]);
    }
#if USE_REALS
    printf(": %f]\n", terms[i]);
#else
    printf(": %d]\n", terms[i]);
#endif
  }

  assert(forest::SUCCESS ==
      evmdd->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE));
  assert(forest::SUCCESS ==
      evmdd->setNodeDeletion(forest::OPTIMISTIC_DELETION));

  timer start;
  start.note_time();
  dd_edge result = test_evmdd(evmdd, compute_manager::PLUS,
      element, terms, nElements);
  printf("Time interval: %.4e seconds\n",
      start.get_last_interval()/1000000.0);

  printf("Peak Nodes in MDD: %d\n", evmdd->getPeakNumNodes());
  printf("Nodes in compute table: %d\n",
      (MEDDLY_getComputeManager())->getNumCacheEntries());

  if (verbose > 1) {
    printf("\n\nForest Info:\n");
    evmdd->showInfo(stdout);
  }

  // Cleanup; in this case simply delete the domain
  delete d;

  free(bounds);
  for (int i = 0; i < nElements; ++i)
  {
    free(element[i]);
  }
  free(element);

  return 0;
}


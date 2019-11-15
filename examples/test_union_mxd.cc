
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
 * test_union_mxd.cc
 *
 * Testing MXD union.
 */

#include <iostream>
#include "../src/meddly.h"
#include "../src/timer.h"

#include <vector>

#define VERBOSE 1

#define CACHE_SIZE 262144u

using namespace MEDDLY;

void printUsage(FILE *outputStream)
{
  fprintf(outputStream,
      "Usage: test_union_mxd <#Variables> <VariableBound> <#Elements>\n");
}

int main(int argc, char *argv[])
{
  if (argc != 4) {
    printUsage(stdout);
    exit(1);
  }

  srand(1u);

  FILE_output meddlyout(stdout);
  
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

  timer mallocTimer;
  long mallocTime = 0;

  int** elements = (int **) malloc(nElements * sizeof(int *));
  int** pelements = (int **) malloc(nElements * sizeof(int *));
  for (int i = 0; i < nElements; ++i)
  {
    mallocTimer.note_time();
    elements[i] = (int *) malloc((nVariables + 1) * sizeof(int));
    pelements[i] = (int *) malloc((nVariables + 1) * sizeof(int));
    mallocTimer.note_time();
    mallocTime += mallocTimer.get_last_interval();

    elements[i][0] = 0;
    pelements[i][0] = 0;
    for (int j = nVariables; j >= 1; --j)
    {
      elements[i][j] = int(float(variableBound) * rand() / (RAND_MAX + 1.0));
      assert(elements[i][j] >= 0 && elements[i][j] < variableBound);
      pelements[i][j] = int(float(variableBound) * rand() / (RAND_MAX + 1.0));
      assert(pelements[i][j] >= 0 && pelements[i][j] < variableBound);
    }
   /*if(i==0)
    {
     pelements[i][2] =  elements[i][2];
    }
    if(i==1)
      {
        elements[i][2] = DONT_CARE;
        pelements[i][2] = DONT_CHANGE;
      }*/
    // print element[i]
#ifdef VERBOSE
    printf("Element %d: [%d", i, elements[i][0]);
    for (int j = 1; j <= nVariables; ++j) { printf(" %d", elements[i][j]); }
    printf("] --> [%d", pelements[i][0]);
    for (int j = 1; j <= nVariables; ++j) { printf(" %d", pelements[i][j]); }
    printf("]\n");
#endif
  }

  // initialize the variable bounds array to provide to the domain

  int* bounds = (int *) malloc(nVariables * sizeof(int));
  assert(bounds != 0);
  for (int i = 0; i < nVariables; ++i)
  {
#ifdef TESTING_AUTO_VAR_GROWTH
    bounds[i] = 2;
#else
    bounds[i] = variableBound;
#endif
  }

  initializer_list* L = defaultInitializerList(0);
#ifdef CACHE_SIZE
  ct_initializer::setMaxSize(CACHE_SIZE);
#endif
  initialize(L);

  // Create a domain
  domain *d = createDomainBottomUp(bounds, nVariables);
  assert(d != 0);

  // Create an MXD forest in this domain (to store a relation)
  forest* xd = d->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL);
  assert(xd != 0);

  dd_edge initial_state(xd);
  timer start;
  printf("Started... ");

  // Create a dd_edge per element and combine using the UNION operator.
  dd_edge** ddElements = new dd_edge*[nElements-1];
  for (int i = 0; i < nElements-1; ++i)
  {
    ddElements[i] = new dd_edge(xd);
    xd->createEdge(elements + i, pelements + i, 1, *(ddElements[i]));
  }

  // Combine dd_edges
  //start.note_time();
  int nDDElements = nElements-1;
  while (nDDElements > 1) {
    int nCombinations = nDDElements/2;
    for (int i = 0; i < nCombinations; i++)
    {
      // Combine i and (nDDElements-1-i)
      (*(ddElements[i])) += (*(ddElements[nDDElements-1-i]));
    }
    nDDElements = (nDDElements+1)/2;
  }
  initial_state = *(ddElements[0]);
  
#ifdef VERBOSE
  printf("\n\nInitial State:\n");
  initial_state.show(meddlyout, 2);
#endif
  
  double c0;
  apply(CARDINALITY, initial_state, c0);
  printf("Elements in result: %.4e\n", c0);
  printf("Peak Nodes in MXD: %ld\n", xd->getPeakNumNodes());
  
  
  printf("\n Will mark time for unioning a single minterm\n");
  start.note_time();
  dd_edge singleminterm(xd);
  xd->createEdge(elements+nElements-1, pelements+nElements-1,1,singleminterm);
  initial_state +=singleminterm;
  start.note_time();  

  for (int i = 0; i < nElements; ++i) delete ddElements[i];
  delete [] ddElements;

#ifdef VERBOSE
  printf("\n\nInitial State Union One Minterm:\n");
  initial_state.show(meddlyout, 2);
#endif
  printf("done. Time interval: %.4e seconds\n",
         start.get_last_seconds());

  double c;
  apply(CARDINALITY, initial_state, c);
  printf("Elements in result: %.4e\n", c);
  printf("Peak Nodes in MXD: %ld\n", xd->getPeakNumNodes());
  
  /* TBD: FIX
  printf("Nodes in compute table: %ld\n",
      (getComputeManager())->getNumCacheEntries());
  */

/*#ifdef VERBOSE
  printf("\n\nForest Info:\n");
  xd->showInfo(meddlyout);
#endif*/
  
  std::vector<int> f_term;
  std::vector<int> t_term;
  
  f_term.push_back(0);f_term.push_back(0);f_term.push_back(0);// f_term.push_back(0);f_term.push_back(0);f_term.push_back(1); f_term.push_back(0);f_term.push_back(0);f_term.push_back(1);f_term.push_back(1);
  t_term.push_back(1);t_term.push_back(1);t_term.push_back(1);//t_term.push_back(1);t_term.push_back(1);t_term.push_back(0);t_term.push_back(1);t_term.push_back(1);t_term.push_back(0);t_term.push_back(0);
  printf("\n Single element: [0 %d %d %d] --> [0 %d %d %d]",f_term.at(0),f_term.at(1),f_term.at(2),t_term.at(0),t_term.at(1),t_term.at(2));
  
  std::vector<std::vector<int>> terms;
  terms.push_back(f_term);
  terms.push_back(t_term);
  
  node_handle fhandle = initial_state.getNode();
  printf("Before  operation handle is: %d\n", fhandle);
  timer start1;
  printf("Started... ");
  start1.note_time();
  node_handle new_handle = xd->unionOneMinterm(fhandle, terms);
  printf("After New operation handle is: %d\n", new_handle);
  start1.note_time();
  printf("done1. Time interval: %.4e seconds\n",
         start1.get_last_seconds());
  dd_edge new_state(xd);
  new_state.set(new_handle);
  
#ifdef VERBOSE
  printf("\n\nNew State:\n");
  new_state.show(meddlyout, 2);
#endif
  
  double c1;
  apply(CARDINALITY, new_state, c1);
  printf("Elements in result: %.4e\n", c1);
  printf("Peak Nodes in MXD: %ld\n", xd->getPeakNumNodes());

  // Cleanup; in this case simply delete the domain
  destroyDomain(d);
  cleanup();

  free(bounds);
  for (int i = 0; i < nElements; ++i)
  {
    mallocTimer.note_time();
    free(elements[i]);
    free(pelements[i]);
    mallocTimer.note_time();
    mallocTime += mallocTimer.get_last_interval();
  }
  free(elements);
  free(pelements);

#ifdef VERBOSE
  printf("Malloc time: %.4e seconds\n", mallocTime/1000000.0);
#endif

  return 0;
}

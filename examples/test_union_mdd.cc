
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
 * test_union_mdd.cc
 *
 * Testing MDD union.
 */

//TODO: test union, intersection (optimize), diff (optimize)

#include <iostream>
#include <vector>
#include <algorithm>

#include "meddly.h"
#include "timer.h"

using namespace MEDDLY;

//#define TESTING_AUTO_VAR_GROWTH

//#define TESTING_TEMP_DD_EDGES
//#define TESTING_UNION_SPEED
//#define BUILD_INDEX_SET
//#define CHECK_ELEMENTS

//#define VERBOSE

void printUsage(FILE *outputStream)
{
  fprintf(outputStream,
      "Usage: test_union_mdd <#Variables> <VariableBound> <#Elements>\n");
}

void printState(const dd_edge& state)
{
	printf("State:\n");
	state.show(stdout, 2);
}

void printStat(const dd_edge& state)
{
	printf("Peak Nodes in MDD: %ld\n", state.getForest()->getPeakNumNodes());
	printf("Current Nodes in MDD: %ld\n", state.getForest()->getCurrentNumNodes());
}

typedef struct CompLt {
	bool operator() (const std::vector<int>& x, const std::vector<int>& y){
		MEDDLY_DCASSERT(x.size()==y.size());

		for(int i=x.size()-1; i>=0 ;i--) {
			if(x[i]<y[i]) {
				return true;
			}
			if(x[i]>y[i]) {
				return false;
			}
		}
		return false;
	}
}CompLt;

void printAssignments(const std::vector<std::vector<int> >& assign)
{
//	std::sort(assign.begin(), assign.end(), CompLt());
	for (std::vector<std::vector<int> >::const_iterator iter=assign.begin(); iter!=assign.end(); iter++) {
		for(std::vector<int>::const_iterator siter=iter->begin(); siter<iter->end(); siter++) {
			std::cout<<*siter<<" ";
		}
		std::cout<<std::endl;
	}
}

void saveAssignments(dd_edge& state, std::vector<std::vector<int> >& assign)
{
	assign.clear();

	int nVariables = state.getForest()->getDomain()->getNumVariables();
	for (enumerator iter(state); iter; ++iter) {
		assign.push_back(std::vector<int>(iter.getAssignments()+1, iter.getAssignments()+1+nVariables));
	}
	std::sort(assign.begin(), assign.end(), CompLt());
}

void verifyAssignments(const std::vector<std::vector<int> >& assign1, const std::vector<std::vector<int> >& assign2)
{
	if(assign1.size()!=assign2.size()) {
		std::cerr<<"Error in Size"<<std::endl;
		exit(1);
	}

	for(int i=0; i<assign1.size(); i++) {
		for(int j=0; j<assign1[i].size(); j++) {
			if(assign1[i][j] != assign2[i][j]) {
				std::cerr<<"Error in Value"<<std::endl;
				exit(1);
			}
		}
	}
}

void printInCount(forest* f)
{
	expert_forest* ef=static_cast<expert_forest*>(f);
	for(int i=1; i<ef->getPeakNumNodes()+1; i++) {
		if(ef->isValidNodeIndex(i) && (ef->isTerminalNode(i) || ef->getNode(i).offset > 0)) {
			std::cout<<"Node "<<i<<": "<<ef->getInCount(i)<<" @ "<<ef->getNodeLevel(i)<<std::endl;
		}
	}
}

//int main(int argc, char *argv[])
//{
//  if (argc != 4) {
//    printUsage(stdout);
//    exit(1);
//  }
//
//  unsigned int seed = time(NULL);
//  srandom(seed);
//
//  printf("Seed: %d\n", seed);
//
//  // initialize number of variables, their bounds and the number of elements
//  // to create
//
//  int nVariables = 0;
//  int variableBound = 0;
//  int nElements = 0;
//
//  sscanf(argv[1], "%d", &nVariables);
//  assert(nVariables > 0);
//
//  sscanf(argv[2], "%d", &variableBound);
//  assert(variableBound > 0);
//
//  sscanf(argv[3], "%d", &nElements);
//  assert(nElements > 0);
//
//  printf("#variables: %d, variable bound: %d, #elements: %d\n",
//      nVariables, variableBound, nElements);
//
//  // create the elements randomly
//
//  timer mallocTimer;
//  long mallocTime = 0;
//
//  int** elements = (int **) malloc(nElements * sizeof(int *));
//  for (int i = 0; i < nElements; ++i)
//  {
//    mallocTimer.note_time();
//    elements[i] = (int *) malloc((nVariables + 1) * sizeof(int));
//    mallocTimer.note_time();
//    mallocTime += mallocTimer.get_last_interval();
//
//    elements[i][0] = 0;
//    for (int j = nVariables; j >= 1; --j)
//    {
//      elements[i][j] = int(float(variableBound) * random() / (RAND_MAX + 1.0));
//      assert(elements[i][j] >= 0 && elements[i][j] < variableBound);
//    }
//    // print element[i]
//#ifdef VERBOSE
//    printf("Element %d: [%d", i, elements[i][0]);
//    for (int j = 1; j <= nVariables; ++j)
//    {
//      printf(" %d", elements[i][j]);
//    }
//    printf("]\n");
//#endif
//  }
//
//  // initialize the variable bounds array to provide to the domain
//
//  int* bounds = (int *) malloc(nVariables * sizeof(int));
//  assert(bounds != 0);
//  for (int i = 0; i < nVariables; ++i)
//  {
//#ifdef TESTING_AUTO_VAR_GROWTH
//    bounds[i] = 2;
//#else
//    bounds[i] = variableBound;
//#endif
//  }
//
//
//  settings s;
//#ifdef CACHE_SIZE
//  s.computeTable.maxSize = CACHE_SIZE;
//#endif
//  initialize(s);
//
//  // Create a domain
//  domain *d = createDomainBottomUp(bounds, nVariables);
//  assert(d != 0);
//
//  // Create an MDD forest in this domain (to store states)
//  forest* f = d->createForest(false, forest::BOOLEAN,
//      forest::MULTI_TERMINAL);
//  assert(f != 0);
//
//#if 0
//  assert(forest::SUCCESS ==
//      //states->setReductionRule(forest::FULLY_REDUCED));
//    f->setReductionRule(forest::QUASI_REDUCED));
//  assert(forest::SUCCESS ==
//      f->setNodeDeletion(forest::OPTIMISTIC_DELETION));
//  // states->setNodeDeletion(forest::PESSIMISTIC_DELETION));
//  if (variableBound < 4) {
//    assert(forest::SUCCESS ==
//        f->setNodeStorage(forest::FULL_STORAGE));
//  } else {
//    assert(forest::SUCCESS ==
//        f->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE));
//  }
//#endif
//
//  dd_edge state(f);
//
//  timer start;
//
//  printf("Started... ");
//
//  // Use Meddly's batch addition to combine all elements in one step.
//  f->createEdge(elements, nElements, state);
//
////  printState(state);
//
//  std::vector<std::vector<int> > assign1;
////  saveAssignments(state, assign1);
//
//  std::vector<std::vector<int> > assign2;
//
////  printf("\n==================== Bubble Reordering ====================\n");
//  int count=1;
////  for(int var=1; var<=nVariables; var++) {
//////	  printf("\n========================= Round %d =========================\n", var);
////	  for(int level=1; level<nVariables; level++) {
////		  printf("%d. Swap variables at level %d and %d\n", count, level, level+1);
////
////		  static_cast<expert_domain*>(d)->swapAdjacentVariables(level);
////
//////		  printState(state);
////
////		  count++;
////
//////		  printState(state);
//////		  printInCount(f);
////	  }
////  }
//
////  printState(state);
////  saveAssignments(state, assign2);
////  verifyAssignments(assign1, assign2);
//
//  printf("\n==================== Reverse Reordering ====================\n");
//  count=1;
//  for(int level=1; level<nVariables; level++) {
//	  printf("%d. Swap variables at level %d and %d\n", count, level, level+1);
//	  static_cast<expert_domain*>(d)->swapAdjacentVariables(level);
//	  count++;
//  }
//  for(int level=nVariables-1; level>=1; level--) {
//	  printf("%d. Swap variables at level %d and %d\n", count, level, level+1);
//	  static_cast<expert_domain*>(d)->swapAdjacentVariables(level);
//	  count++;
//  }
//
////  printState(state);
////  saveAssignments(state, assign2);
////  verifyAssignments(assign1, assign2);
//
//  // Cleanup; in this case simply delete the domain
//  destroyDomain(d);
//  MEDDLY::cleanup();
//
//  free(bounds);
//  for (int i = 0; i < nElements; ++i) {
//    mallocTimer.note_time();
//    free(elements[i]);
//    mallocTimer.note_time();
//    mallocTime += mallocTimer.get_last_interval();
//  }
//  free(elements);
//
//#ifdef VERBOSE
//  printf("Malloc time: %.4e seconds\n", mallocTime/1000000.0);
//#endif
//
//  return 0;
//}

//int main(int argc, char *argv[])
//{
//  if (argc != 4) {
//    printUsage(stdout);
//    exit(1);
//  }
//
//  unsigned int seed = time(NULL);
//  srandom(1408565459u);
//
//  printf("Seed: %d\n", seed);
//
//  // initialize number of variables, their bounds and the number of elements
//  // to create
//
//  int nVariables = 0;
//  int variableBound = 0;
//  int nElements = 0;
//
//  sscanf(argv[1], "%d", &nVariables);
//  assert(nVariables > 0);
//
//  sscanf(argv[2], "%d", &variableBound);
//  assert(variableBound > 0);
//
//  sscanf(argv[3], "%d", &nElements);
//  assert(nElements > 0);
//
//  printf("#variables: %d, variable bound: %d, #elements: %d\n",
//      nVariables, variableBound, nElements);
//
//  // create the elements randomly
//
//  int** elements = (int **) malloc(nElements * sizeof(int *));
//  for (int i = 0; i < nElements; ++i)
//  {
//    elements[i] = (int *) malloc((nVariables + 1) * sizeof(int));
//
//    elements[i][0] = 0;
//    for (int j = nVariables; j >= 1; --j)
//    {
//      elements[i][j] = int(float(variableBound+1) * random() / (RAND_MAX + 1.0))-1;
//      assert(elements[i][j] >= DONT_CARE && elements[i][j] < variableBound);
//    }
//    // print element[i]
//#ifdef VERBOSE
//    printf("Element %d: [%d", i, elements[i][0]);
//    for (int j = 1; j <= nVariables; ++j)
//    {
//      printf(" %d", elements[i][j]);
//    }
//    printf("]\n");
//#endif
//  }
//
//  // initialize the variable bounds array to provide to the domain
//
//  int* bounds = (int *) malloc(nVariables * sizeof(int));
//  assert(bounds != 0);
//  for (int i = 0; i < nVariables; ++i) {
//    bounds[i] = variableBound;
//  }
//
//
//  settings s;
//  initialize(s);
//
//  // Create a domain
//  domain *d = createDomainBottomUp(bounds, nVariables);
//  assert(d != 0);
//
//  // Create an MDD forest in this domain (to store states)
//  forest* f = d->createForest(false, forest::BOOLEAN,
//      forest::MULTI_TERMINAL);
//  assert(f != 0);
//
//#if 0
//  assert(forest::SUCCESS ==
//      //states->setReductionRule(forest::FULLY_REDUCED));
//    f->setReductionRule(forest::QUASI_REDUCED));
//  assert(forest::SUCCESS ==
//      f->setNodeDeletion(forest::OPTIMISTIC_DELETION));
//  // states->setNodeDeletion(forest::PESSIMISTIC_DELETION));
//  if (variableBound < 4) {
//    assert(forest::SUCCESS ==
//        f->setNodeStorage(forest::FULL_STORAGE));
//  } else {
//    assert(forest::SUCCESS ==
//        f->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE));
//  }
//#endif
//
//  dd_edge state(f);
//
//  timer start;
//
//  printf("Started... ");
//
//  // Use Meddly's batch addition to combine all elements in one step.
//  f->createEdge(elements, nElements, state);
//
//  static_cast<expert_forest*>(f)->removeAllComputeTableEntries();
//
//  printState(state);
//  printStat(state);
//
//  std::vector<std::vector<int> > assign1;
//  saveAssignments(state, assign1);
////  printAssignments(assign1);
//
//  std::cout<<std::endl;
//
//  std::vector<std::vector<int> > assign2;
//
////  std::cout<<"==================== Bubble Reordering ===================="<<std::endl;
////  for(int i=0; i<nVariables; i++) {
////	  int high = nVariables;
////	  int low = 1;
////
////	  std::cout<<i<<". Move the variable at level "<<high<<" down to level "<<low <<std::endl;
////
////	  static_cast<expert_domain*>(d)->moveDownVariable(high, low);
////
//////	  printState(state);
////	  printStat(state);
////
////	  saveAssignments(state, assign2);
//////	  printAssignments(assign2);
////	  verifyAssignments(assign1, assign2);
////
//////	  static_cast<expert_forest*>(f)->dumpUniqueTable(stdout);
////  }
//
////  printState(state);
////  printStat(state);
////  printInCount(f);
////  printAssignments(assign2);
//
//  printf("\n==================== Random Reordering ====================\n");
//  for(int i=0; i<1; i++) {
////	  int high = int(float(nVariables-1) * random() / (RAND_MAX + 1.0))+2;
////	  int low = int(float(high-1) * random() / (RAND_MAX + 1.0))+1;
//	  int high = 3;
//	  int low = 1;
//
//	  std::cout<<i<<". Move the variable at level "<<low<<" up to level "<<high <<std::endl;
//
//	  static_cast<expert_domain*>(d)->moveUpVariable(low, high);
//
//	  printState(state);
//	  printStat(state);
//
//	  saveAssignments(state, assign2);
////	  printAssignments(assign2);
//	  verifyAssignments(assign1, assign2);
//
////	  static_cast<expert_forest*>(f)->dumpUniqueTable(stdout);
//  }
//
////  printState(state);
////  printStat(state);
////  printInCount(f);
////  printAssignments(assign2);
//
//  // Cleanup; in this case simply delete the domain
//  destroyDomain(d);
//  MEDDLY::cleanup();
//
//  free(bounds);
//  for (int i = 0; i < nElements; ++i) {
//    free(elements[i]);
//  }
//  free(elements);
//
//  return 0;
//}

//int main(int argc, char *argv[])
//{
//  if (argc != 4) {
//    printUsage(stdout);
//    exit(1);
//  }
//
//  unsigned int seed = time(NULL);
//  srandom(seed);
//
//  printf("Seed: %d\n", seed);
//
//  // initialize number of variables, their bounds and the number of elements
//  // to create
//
//  int nVariables = 8;
//  int variableBound = 2;
//  int nElements = 15;
//
//  printf("#variables: %d, variable bound: %d, #elements: %d\n",
//      nVariables, variableBound, nElements);
//
//  // create the elements randomly
//
//  int** elements = (int **) malloc(nElements * sizeof(int *));
//
//  int es[][9] = {{1, 1, DONT_CARE, DONT_CARE, DONT_CARE, DONT_CARE, DONT_CARE, DONT_CARE},
//		  {1, 0, 1, 1, DONT_CARE, DONT_CARE, DONT_CARE, DONT_CARE},
//		  {0, DONT_CARE, 1, 1, DONT_CARE, DONT_CARE, DONT_CARE, DONT_CARE},
//		  {1, 0, 1, 0, 1, 1, DONT_CARE, DONT_CARE},
//		  {1, 0, 0, DONT_CARE, 1, 1, DONT_CARE, DONT_CARE},
//		  {0, DONT_CARE, 1, 0, 1, 1, DONT_CARE, DONT_CARE},
//		  {0, DONT_CARE, 0, DONT_CARE, 1, 1, DONT_CARE, DONT_CARE},
//		  {1, 0, 1, 0, 1, 0, 1, 1},
//		  {0, DONT_CARE, 1, 0, 1, 0, 1, 1},
//		  {1, 0, 0, DONT_CARE, 1, 0, 1, 1},
//		  {0, DONT_CARE, 0, DONT_CARE, 1, 0, 1, 1},
//		  {1, 0, 1, 0, 0, DONT_CARE, 1, 1},
//		  {0, DONT_CARE, 1, 0, 0, DONT_CARE, 1, 1},
//		  {1, 0, 0, DONT_CARE, 0, DONT_CARE, 1, 1},
//		  {0, DONT_CARE, 0, DONT_CARE, 0, DONT_CARE, 1, 1}
//    };
//
//  for(int i=0; i<nElements; i++) {
//	  elements[i] = (int *) malloc((nVariables + 1) * sizeof(int));
//	  for(int j=0; j<nVariables+1; j++) {
//		  elements[i][j] = es[i][j];
//	  }
//  }
//
//  // initialize the variable bounds array to provide to the domain
//
//  int* bounds = (int *) malloc(nVariables * sizeof(int));
//  assert(bounds != 0);
//  for (int i = 0; i < nVariables; ++i)
//  {
//#ifdef TESTING_AUTO_VAR_GROWTH
//    bounds[i] = 2;
//#else
//    bounds[i] = variableBound;
//#endif
//  }
//
//
//  settings s;
//#ifdef CACHE_SIZE
//  s.computeTable.maxSize = CACHE_SIZE;
//#endif
//  initialize(s);
//
//  // Create a domain
//  domain *d = createDomainBottomUp(bounds, nVariables);
//  assert(d != 0);
//
//  // Create an MDD forest in this domain (to store states)
//  forest* f = d->createForest(false, forest::BOOLEAN,
//      forest::MULTI_TERMINAL);
//  assert(f != 0);
//
//#if 0
//  assert(forest::SUCCESS ==
//      //states->setReductionRule(forest::FULLY_REDUCED));
//    f->setReductionRule(forest::QUASI_REDUCED));
//  assert(forest::SUCCESS ==
//      f->setNodeDeletion(forest::OPTIMISTIC_DELETION));
//  // states->setNodeDeletion(forest::PESSIMISTIC_DELETION));
//  if (variableBound < 4) {
//    assert(forest::SUCCESS ==
//        f->setNodeStorage(forest::FULL_STORAGE));
//  } else {
//    assert(forest::SUCCESS ==
//        f->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE));
//  }
//#endif
//
//  dd_edge state(f);
//
//  timer start;
//
//  printf("Started... ");
//
//  // Use Meddly's batch addition to combine all elements in one step.
//  f->createEdge(elements, nElements, state);
//
//  static_cast<expert_forest*>(f)->removeAllComputeTableEntries();
//
//  printState(state);
//  printStat(state);
//
//  std::vector<std::vector<int> > assign1;
//  saveAssignments(state, assign1);
////  printAssignments(assign1);
//
//  printInCount(f);
//
//  std::cout<<std::endl;
//
//  std::vector<std::vector<int> > assign2;
//
////  std::cout<<"==================== Bubble Reordering ===================="<<std::endl;
////  for(int i=0; i<nVariables; i++) {
////	  int high = nVariables;
////	  int low = 1;
////
////	  std::cout<<i<<". Move the variable at level "<<high<<" down to level "<<low <<std::endl;
////
////	  static_cast<expert_domain*>(d)->moveDownVariable(high, low);
////
//////	  printState(state);
////	  printStat(state);
////
////	  saveAssignments(state, assign2);
//////	  printAssignments(assign2);
////	  verifyAssignments(assign1, assign2);
////
//////	  static_cast<expert_forest*>(f)->dumpUniqueTable(stdout);
////  }
//
////  printState(state);
////  printStat(state);
////  printInCount(f);
////  printAssignments(assign2);
//
//  printf("\n==================== Random Reordering ====================\n");
//  for(int i=0; i<1; i++) {
//	  int high = int(float(nVariables-1) * random() / (RAND_MAX + 1.0))+2;
//	  int low = int(float(high-1) * random() / (RAND_MAX + 1.0))+1;
//
//	  std::cout<<i<<". Move the variable at level "<<high<<" down to level "<<low <<std::endl;
//
//	  static_cast<expert_domain*>(d)->moveUpVariable(low, high);
//
////	  printState(state);
//	  printStat(state);
//
//	  saveAssignments(state, assign2);
////	  printAssignments(assign2);
//	  verifyAssignments(assign1, assign2);
//
////	  static_cast<expert_forest*>(f)->dumpUniqueTable(stdout);
//  }
//
////  printState(state);
////  printStat(state);
////  printInCount(f);
////  printAssignments(assign2);
//
//  // Cleanup; in this case simply delete the domain
//  destroyDomain(d);
//  MEDDLY::cleanup();
//
//  free(bounds);
//  for (int i = 0; i < nElements; ++i) {
//    free(elements[i]);
//  }
//  free(elements);
//
//#ifdef VERBOSE
//  printf("Malloc time: %.4e seconds\n", mallocTime/1000000.0);
//#endif
//
//  return 0;
//}
//
//int main(int argc, char *argv[])
//{
//  if (argc != 4) {
//    printUsage(stdout);
//    exit(1);
//  }
//
//  unsigned int seed = time(NULL);
//  srandom(1408735565u);
//
//  printf("Seed: %d\n", seed);
//
//  timer runTimer;
//
//  // initialize number of variables, their bounds and the number of elements
//  // to create
//
//  int nVariables = 0;
//  int variableBound = 0;
//  int nElements = 0;
//
//  sscanf(argv[1], "%d", &nVariables);
//  assert(nVariables > 0);
//
//  sscanf(argv[2], "%d", &variableBound);
//  assert(variableBound > 0);
//
//  sscanf(argv[3], "%d", &nElements);
//  assert(nElements > 0);
//
//  printf("#variables: %d, variable bound: %d, #elements: %d\n",
//      nVariables, variableBound, nElements);
//
//  // create the elements randomly
//
//  int** elements = (int **) malloc(nElements * sizeof(int *));
//  for (int i = 0; i < nElements; ++i)
//  {
//    elements[i] = (int *) malloc((nVariables + 1) * sizeof(int));
//
//    elements[i][0] = 0;
//    for (int j = nVariables; j >= 1; --j)
//    {
//      elements[i][j] = int(float(variableBound+1) * random() / (RAND_MAX + 1.0))-1;
//      assert(elements[i][j] >= DONT_CARE && elements[i][j] < variableBound);
//    }
//    // print element[i]
//#ifdef VERBOSE
//    printf("Element %d: [%d", i, elements[i][0]);
//    for (int j = 1; j <= nVariables; ++j)
//    {
//      printf(" %d", elements[i][j]);
//    }
//    printf("]\n");
//#endif
//  }
//
//  // initialize the variable bounds array to provide to the domain
//
//  int* bounds = (int *) malloc(nVariables * sizeof(int));
//  assert(bounds != 0);
//  for (int i = 0; i < nVariables; ++i) {
//    bounds[i] = variableBound;
//  }
//
//
//  settings s;
//  initialize(s);
//
//  // Create a domain
//  domain *d = createDomainBottomUp(bounds, nVariables);
//  assert(d != 0);
//
//  // Create an MDD forest in this domain (to store states)
//  forest* f = d->createForest(false, forest::BOOLEAN,
//      forest::MULTI_TERMINAL);
//  assert(f != 0);
//
//#if 0
//  assert(forest::SUCCESS ==
//      //states->setReductionRule(forest::FULLY_REDUCED));
//    f->setReductionRule(forest::QUASI_REDUCED));
//  assert(forest::SUCCESS ==
//      f->setNodeDeletion(forest::OPTIMISTIC_DELETION));
//  // states->setNodeDeletion(forest::PESSIMISTIC_DELETION));
//  if (variableBound < 4) {
//    assert(forest::SUCCESS ==
//        f->setNodeStorage(forest::FULL_STORAGE));
//  } else {
//    assert(forest::SUCCESS ==
//        f->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE));
//  }
//#endif
//
//  dd_edge state(f);
//
//  timer start;
//
//  printf("Started... ");
//
//  // Use Meddly's batch addition to combine all elements in one step.
//  f->createEdge(elements, nElements, state);
//
//  static_cast<expert_forest*>(f)->removeAllComputeTableEntries();
//
////  printState(state);
//  printStat(state);
//
//  std::vector<std::vector<int> > assign1;
////  saveAssignments(state, assign1);
////  printAssignments(assign1);
//
//  std::cout<<std::endl;
//
//  std::vector<std::vector<int> > assign2;
//
////  std::cout<<"==================== Bubble Reordering ===================="<<std::endl;
////  for(int i=0; i<nVariables; i++) {
////	  int high = nVariables;
////	  int low = 1;
////
////	  std::cout<<i<<". Move the variable at level "<<high<<" down to level "<<low <<std::endl;
////
////	  static_cast<expert_domain*>(d)->moveDownVariable(high, low);
////
//////	  printState(state);
////	  printStat(state);
////
////	  saveAssignments(state, assign2);
//////	  printAssignments(assign2);
////	  verifyAssignments(assign1, assign2);
////
//////	  static_cast<expert_forest*>(f)->dumpUniqueTable(stdout);
////  }
//
////  printState(state);
////  printStat(state);
////  printInCount(f);
////  printAssignments(assign2);
//
////  printState(state);
//
//  printf("\n==================== Swapping ====================\n");
//
//  printf("Swap variables at level %d and %d\n", 1, nVariables);
//  runTimer.note_time();
//  for(int level=1; level<nVariables; level++) {
//  	  static_cast<expert_domain*>(d)->swapAdjacentVariables(level);
//  }
//  runTimer.note_time();
//  std::cout<<"Time: "<<runTimer.get_last_interval()/1000000.0<<" s"<<std::endl;
//  printStat(f);
//
//  printf("Swap variables at level %d and %d\n", nVariables, 1);
//  runTimer.note_time();
//  for(int level=nVariables-1; level>=1; level--) {
//  	  static_cast<expert_domain*>(d)->swapAdjacentVariables(level);
//  }
//  runTimer.note_time();
//  std::cout<<"Time: "<<runTimer.get_last_interval()/1000000.0<<" s"<<std::endl;
//  printStat(f);
//
//  printf("\n==================== Relocating ====================\n");
//
//  std::cout<<"Move the variable at level "<<1<<" up to level "<<nVariables<<std::endl;
//  runTimer.note_time();
//  static_cast<expert_domain*>(d)->moveUpVariable(1, nVariables);
////  for(int level=1; level<nVariables; level++) {
////  	  static_cast<expert_domain*>(d)->moveUpVariable(level, level+1);
////  }
//  runTimer.note_time();
//  std::cout<<"Time: "<<runTimer.get_last_interval()/1000000.0<<" s"<<std::endl;
//  printStat(f);
//
////  std::cout<<"Move the variable at level "<<nVariables<<" down to level "<<1<<std::endl;
////  runTimer.note_time();
//////  static_cast<expert_domain*>(d)->moveDownVariable(nVariables, 1);
////  for(int level=nVariables-1; level>=1; level--) {
////	  static_cast<expert_domain*>(d)->moveDownVariable(level+1, level);
////  }
////  runTimer.note_time();
////  std::cout<<"Time: "<<runTimer.get_last_interval()/1000000.0<<" s"<<std::endl;
////  printStat(f);
//
//
////    printf("\n==================== Reverse Reordering ====================\n");
////    for(int level=1; level<nVariables; level++) {
////  	  printf("Swap variables at level %d and %d\n", level, level+1);
////  	  static_cast<expert_domain*>(d)->swapAdjacentVariables(level);
////    }
////    for(int level=nVariables-1; level>=1; level--) {
////  	  printf("Swap variables at level %d and %d\n", level, level+1);
////  	  static_cast<expert_domain*>(d)->swapAdjacentVariables(level);
////    }
////
////    printf("\n==================== Reverse Reordering ====================\n");
////    for(int level=1; level<nVariables; level++) {
////  	  printf("Swap variables at level %d and %d\n", level, level+1);
////  	  static_cast<expert_domain*>(d)->moveUpVariable(level, level+1);
////    }
////    for(int level=nVariables-1; level>=1; level--) {
////  	  printf("Swap variables at level %d and %d\n", level, level+1);
////  	  static_cast<expert_domain*>(d)->moveDownVariable(level+1, level);
////    }
//
////  printState(state);
//
//
//
////  printState(state);
//
////  printf("\n==================== Random Reordering ====================\n");
////  for(int i=0; i<50; i++) {
////	  int high = int(float(nVariables-1) * random() / (RAND_MAX + 1.0))+2;
////	  int low = int(float(high-1) * random() / (RAND_MAX + 1.0))+1;
////
////	  int moveUp = int(float(2) * random() / (RAND_MAX + 1.0));
////	  if(moveUp==0) {
////		  std::cout<<i<<". Move the variable at level "<<high<<" down to level "<<low<<std::endl;
////		  static_cast<expert_domain*>(d)->moveDownVariable(high, low);
////	  }
////	  else {
////		  std::cout<<i<<". Move the variable at level "<<low<<" up to level "<<high<<std::endl;
////		  static_cast<expert_domain*>(d)->moveUpVariable(low, high);
////	  }
////
//////	  printState(state);
////	  printStat(state);
////
//////	  static_cast<expert_domain*>(d)->moveDownVariable(high, low);
////
//////	  printState(state);
//////	  printStat(state);
////
//////	  saveAssignments(state, assign2);
//////	  printAssignments(assign2);
//////	  verifyAssignments(assign1, assign2);
////
//////	  static_cast<expert_forest*>(f)->dumpUniqueTable(stdout);
////  }
//
////  printState(state);
////  printStat(state);
////  printInCount(f);
////  printAssignments(assign2);
//
//  // Cleanup; in this case simply delete the domain
//  destroyDomain(d);
//  MEDDLY::cleanup();
//
//  free(bounds);
//  for (int i = 0; i < nElements; ++i) {
//    free(elements[i]);
//  }
//  free(elements);
//
//  return 0;
//}

void shuffle(int* array, int start, int end)
{
	for(int i=start; i<end; i++) {
		int swap=i+rand()%(end-i+1);
		int value=array[i];
		array[i]=array[swap];
		array[swap]=value;
	}
}

int main(int argc, char *argv[])
{
  timer runTimer;

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

  const char* reorder = argv[4];

  unsigned int seed = 0;
  if(argc>5) {
	  sscanf(argv[5], "%u", &seed);
  }
  else {
	  seed = time(NULL);
  }
  srandom(seed);

  printf("Seed: %d\n", seed);

  // create the elements randomly

  int** elements = (int **) malloc(nElements * sizeof(int *));
  for (int i = 0; i < nElements; ++i)
  {
    elements[i] = (int *) malloc((nVariables + 1) * sizeof(int));

    elements[i][0] = 0;
    for (int j = nVariables; j >= 1; --j)
    {
    	int prob=rand()%8;
    	if(prob<=2){
    		elements[i][j]=0;
    	}
    					else if(prob<=5){
    						elements[i][j]=1;
    					}
    					else{
    						elements[i][j]=DONT_CARE;
    					}

      assert(elements[i][j] >= DONT_CARE && elements[i][j] < variableBound);
    }
    // print element[i]
#ifdef VERBOSE
    printf("Element %d: [%d", i, elements[i][0]);
    for (int j = 1; j <= nVariables; ++j)
    {
      printf(" %d", elements[i][j]);
    }
    printf("]\n");
#endif
  }

  // initialize the variable bounds array to provide to the domain

  int* bounds = (int *) malloc(nVariables * sizeof(int));
  assert(bounds != 0);
  for (int i = 0; i < nVariables; ++i) {
    bounds[i] = variableBound;
  }

  settings s;
  initialize(s);

  // Create a domain
  domain *d = createDomainBottomUp(bounds, nVariables);
  assert(d != 0);

  // Create an MDD forest in this domain (to store states)
  forest::policies p(false);
  if(strcmp(reorder, "LOWEST_INVERSION")==0) {
	  p.setLowestInversion();
  }
  else if(strcmp(reorder, "HIGHEST_INVERSION")==0) {
	  p.setHighestInversion();
  }
  else if(strcmp(reorder, "BUBBLE_DOWN")==0) {
	  p.setBubbleDown();
  }
  else if(strcmp(reorder, "BUBBLE_UP")==0) {
	  p.setBubbleUp();
  }
  else if(strcmp(reorder, "LOWEST_COST")==0) {
	  p.setLowestCost();
  }

  forest* f = d->createForest(false, forest::BOOLEAN,
      forest::MULTI_TERMINAL, p);
  assert(f != 0);

#if 0
  assert(forest::SUCCESS ==
      //states->setReductionRule(forest::FULLY_REDUCED));
    f->setReductionRule(forest::QUASI_REDUCED));
  assert(forest::SUCCESS ==
      f->setNodeDeletion(forest::OPTIMISTIC_DELETION));
  // states->setNodeDeletion(forest::PESSIMISTIC_DELETION));
  if (variableBound < 4) {
    assert(forest::SUCCESS ==
        f->setNodeStorage(forest::FULL_STORAGE));
  } else {
    assert(forest::SUCCESS ==
        f->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE));
  }
#endif

  dd_edge state(f);

  printf("Started...\n");

  // Use Meddly's batch addition to combine all elements in one step.
  f->createEdge(elements, nElements, state);

  static_cast<expert_forest*>(f)->removeAllComputeTableEntries();
  printStat(state);
//  printState(state);

  int* array=static_cast<int*>(calloc(sizeof(int), nVariables+1));
  for(int i=0; i<nVariables+1; i++) {
	  array[i]=i;
  }

  shuffle(array, 1, nVariables);

//  for(int i=1; i<nVariables+1; i++) {
//	printf("Var %d : Level %d\n", i, array[i]);
//  }

  runTimer.note_time();

  static_cast<expert_domain*>(d)->reorderVariables(array);

  runTimer.note_time();
  std::cout<<"Time: "<<runTimer.get_last_interval()/1000000.0<<" s"<<std::endl;

  printStat(state);

  // Cleanup; in this case simply delete the domain
  destroyDomain(d);
  MEDDLY::cleanup();

  free(array);

  free(bounds);
  for (int i = 0; i < nElements; ++i) {
    free(elements[i]);
  }
  free(elements);

  return 0;
}

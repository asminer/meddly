
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
 * test_temporary_nodes.cc
 *
 * Testing operations on Temporary nodes in Meddly's Expert interface.
 */

#include <iostream>
#include <string.h>
#include "meddly.h"
#include "meddly_expert.h"

#include "src/forests/mtmxd.h"

// Timer class
#include "timer.h"

#include <set>

//#define VERBOSE

const int trueNode = -1;

void findNodesInGraph(expert_forest* f, int node, std::set<int>& result);

void convertDDEdgeToTemporaryNode(const dd_edge& a, int& b);
int convertToTemporaryNode(std::map<int, int>& cache,
    expert_forest* f, int root);
void printUsage(FILE *outputStream);

// TestA:
// Builds temporary nodes by hand.
// Adds elements on at a time and prints the results.
void testA(expert_forest* f);

// TestB:
// Builds temporary nodes by hand (MXDs).
// Adds elements on at a time and prints the results.
void testB();

// PhaseI: use convertDDEdgeToTemporaryNode.
void doPhaseI(expert_forest* f, const int* const* elements,
    int start, int end, int& tempNode);
// PhaseII: use accumulate += minterms
void doPhaseII(expert_forest* f, const int* const* elements,
    int start, int end, int& tempNode);
// Phase III: use acccumulate += mdd
void doPhaseIII(expert_forest* f, const int* const* elements,
    int start, int end, int& tempNode);

void printIncounts(expert_forest* f, int node);

int main(int argc, char *argv[])
{
#if 1
  testB();
#else
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

  int** elements = (int **) malloc(nElements * sizeof(int *));
  for (int i = 0; i < nElements; ++i)
  {
    elements[i] = (int *) malloc((nVariables + 1) * sizeof(int));
    elements[i][0] = 0;
    for (int j = nVariables; j >= 1; --j)
    {
      elements[i][j] = int(float(variableBound) * random() / (RAND_MAX + 1.0));
      assert(elements[i][j] >= 0 && elements[i][j] < variableBound);
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

#ifdef VERBOSE
  printf("#########################################################\n\n");
#endif

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

  // Create an MDD forest in this domain (to store states)
  forest* states = d->createForest(false, forest::BOOLEAN,
      forest::MULTI_TERMINAL);
  assert(states != 0);

  assert(forest::SUCCESS ==
      states->setReductionRule(forest::FULLY_REDUCED));
    //states->setReductionRule(forest::QUASI_REDUCED));
#if 0
  assert(forest::SUCCESS ==
      states->setNodeDeletion(forest::OPTIMISTIC_DELETION));
  // states->setNodeDeletion(forest::PESSIMISTIC_DELETION));
  if (variableBound < 4) {
    assert(forest::SUCCESS ==
        states->setNodeStorage(forest::FULL_STORAGE));
  } else {
    assert(forest::SUCCESS ==
        states->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE));
  }
#endif

  expert_forest* expertStates = dynamic_cast<expert_forest*>(states);
  assert(0 != expertStates);

#if 0
  testA(expertStates);
#endif

  // Tests:
  // (1) temp_node += min_terms
  // (2) temp_node += mdd
  // (3) reduce(temp_node)
  //
  // (0) convertDDEdgeToTemporaryNode()
  // Use createEdge for nEments/3.
  // Convert the above to temp_node.
  //
  // (1) temp_node += min_terms:
  // Add nElements/3 to temp_node
  //
  // (2) temp_node += mdd
  // Add nElements/3 to mdd using createEdge.
  // Do temp_node += mdd
  //
  // (3) reduce(temp_node)
  // assert(reduce(temp_node) == createEdge(nElements))

  int tempNode = 0;
  int twoNBy3 = 2 * nElements / 3;

  // (0)
  doPhaseI(expertStates, elements, 0, nElements/3, tempNode);

  // (1)
  doPhaseII(expertStates, elements, nElements/3, twoNBy3, tempNode);

  // (2)
  doPhaseIII(expertStates, elements, twoNBy3, nElements, tempNode);

  // (3)
  int accumulatedNode = expertStates->reduceNode(tempNode);

#ifdef VERBOSE
  printf("Reduced tempNode\n");
  printf("----------------\n");
  expertStates->showNodeGraph(stdout, accumulatedNode);
  printf("\n");

  printf("#########################################################\n\n");
#endif

#ifdef DEBUG_ACCUMULATE_MDD
  printf("tempNode: %d, reducedNode: %d\n", tempNode, accumulatedNode);
  if (expertStates->isActiveNode(tempNode)) {
    printf("Incount(tempNode) after reduction: %d\n\n",
        expertStates->getInCount(tempNode));
  }
#endif

  dd_edge final(states);
  final.set(accumulatedNode, 0,
      expertStates->getNodeLevel(accumulatedNode));

  dd_edge nodeC(states);
  assert(forest::SUCCESS ==
      states->createEdge(elements, nElements, nodeC));

  if (final != nodeC) {
#ifdef VERBOSE
    printf("createEdge(nElements)\n");
    printf("---------------------\n");
    expertStates->showNodeGraph(stdout, nodeC.getNode());
    printf("\n");
    printf("#########################################################\n\n");
#endif
    printf("ERROR: Reduced tempNode != createEdge(nElements)\n");
    dd_edge diff = final - nodeC;
#if 0
    printf("\nAccumulate - createEdge:\n");
    expertStates->showNodeGraph(stdout, diff.getNode());
#else
    printf("\nAccumulate - createEdge: %d elements\n",
        (int)diff.getCardinality());
#endif
    diff = nodeC - final;
#if 0
    printf("\ncreateEdge - Accumulate:\n");
    expertStates->showNodeGraph(stdout, diff.getNode());
#else
    printf("\ncreateEdge - Accumulate: %d elements\n",
        (int)diff.getCardinality());
#endif
  }

  // Done with final and nodeC, clear them.
  final.clear();
  nodeC.clear();

  // Cleanup; in this case simply delete the domain
  delete d;

  free(bounds);
  for (int i = 0; i < nElements; ++i)
  {
    free(elements[i]);
  }
  free(elements);

#endif
  return 0;
}


void printUsage(FILE *outputStream)
{
  fprintf(outputStream,
      "Usage: test_temporary_nodes <#Variables> <VariableBound> <#Elements>\n");
}


int convertToTemporaryNode(std::map<int, int>& cache,
    expert_forest* f, int root)
{
  // Terminal nodes
  if (f->isTerminalNode(root)) return root;
  // Check cache for result
  std::map<int, int>::iterator iter = cache.find(root);
  if (iter != cache.end()) {
    f->linkNode(iter->second);
    return iter->second;
  }
  // Build new temporary node
  int level = f->getNodeLevel(root);
  int result = 0;
  if (f->isFullNode(root)) {
    int size = f->getFullNodeSize(root);
    result = f->createTempNode(level, size, false);
    for (int i = 0; i < size; i++)
    {
      int temp = convertToTemporaryNode(cache, f, f->getDownPtr(root, i));
      f->setDownPtrWoUnlink(result, i, temp);
      f->unlinkNode(temp);
    }
  }
  else {
    assert(f->isSparseNode(root));
    int sparseSize = f->getSparseNodeSize(root);
    int size = 1 + f->getSparseNodeIndex(root, sparseSize - 1);
    result = f->createTempNode(level, size, true);
    for (int i = 0; i < sparseSize; i++)
    {
      int temp = convertToTemporaryNode(cache, f,
          f->getSparseNodeDownPtr(root, i));
      f->setDownPtrWoUnlink(result, f->getSparseNodeIndex(root, i), temp);
      f->unlinkNode(temp);
    }
  }
  // Store result in cache
  cache[root] = result;
  return result;
}


void convertDDEdgeToTemporaryNode(const dd_edge& a, int& b)
{
  // Recursive procedure:
  // Start at root(a).
  // Build a temporary node for root(a).
  // -- Build a temporary node for each child(root(a)).
  // Keep track of duplicates via a compute cache which maps
  //   each reduced node to a temporary node.
  expert_forest* f = static_cast<expert_forest*>(a.getForest());
  int root = a.getNode();
  std::map<int, int> cache;
  b = convertToTemporaryNode(cache, f, root);
}


void findNodesInGraph(expert_forest* f, int node, std::set<int>& result)
{
  result.clear();
  std::set<int> toVisit;
  if (!f->isTerminalNode(node)) toVisit.insert(node);
  while (!toVisit.empty()) {
    std::set<int>::iterator iter = toVisit.begin();
    int next = *iter;
    toVisit.erase(iter);
    result.insert(next);
    if (f->isFullNode(next)) {
      int size = f->getFullNodeSize(next);
      for (int i = 0; i < size; ++i) {
        int dptr = f->getFullNodeDownPtr(next, i);
        if (!f->isTerminalNode(dptr))
          if (result.find(dptr) == result.end())
            toVisit.insert(dptr);
      }
    }
    else {
      assert(f->isSparseNode(next));
      int size = f->getSparseNodeSize(next);
      for (int i = 0; i < size; ++i) {
        int dptr = f->getSparseNodeDownPtr(next, i);
        if (!f->isTerminalNode(dptr))
          if (result.find(dptr) == result.end())
            toVisit.insert(dptr);
      }
    }
  }
}


void testA(expert_forest* f)
{
  // Domain must have 3 variables of size 2.

  int nVars = 3;
  assert(f->getDomain()->getNumVariables() == nVars);

  int varSize = 2;
  int level1 = f->getDomain()->getVariableAbove(domain::TERMINALS);
  int level2 = f->getDomain()->getVariableAbove(level1);
  int level3 = f->getDomain()->getVariableAbove(level2);

  assert(f->getDomain()->getVariableBound(level1) == varSize);
  assert(f->getDomain()->getVariableBound(level2) == varSize);
  assert(f->getDomain()->getVariableBound(level3) == varSize);

  {
    // --------------------------------------------------------------------
    // TEST #1
    //
    // Build a node for the set: {000, 010}
    int trueNode = f->getTerminalNode(true);
    int node1L1 = f->createTempNodeMaxSize(level1);
    int node2L2 = f->createTempNodeMaxSize(level2);
    int node3L3 = f->createTempNodeMaxSize(level3);

    // node1L1 = {0} --> T
    f->setDownPtr(node1L1, 0, trueNode);

    // node2L2 = {0,1} --> node1L1
    f->setDownPtr(node2L2, 0, node1L1);
    f->setDownPtr(node2L2, 1, node1L1);

    // node3L3 = {0} --> node2L2
    f->setDownPtr(node3L3, 0, node2L2);

    int tempNode = node3L3;
    f->linkNode(tempNode);

    f->unlinkNode(node3L3);
    f->unlinkNode(node2L2);
    f->unlinkNode(node1L1);
    f->unlinkNode(trueNode);

    printf("{000, 010}\n");
    f->showNodeGraph(stdout, tempNode);
    printIncounts(f, tempNode);

    // Add 001
    int* element = new int[nVars + 1];
    element[0] = 0;
    element[1] = 1;
    element[2] = 0;
    element[3] = 0;
    f->accumulate(tempNode, element);

    printf("Added {001}\n");
    f->showNodeGraph(stdout, tempNode);
    printIncounts(f, tempNode);

    // Add 011
    element[1] = 1;
    element[2] = 1;
    element[3] = 0;
    f->accumulate(tempNode, element);

    printf("Added {011}\n");
    f->showNodeGraph(stdout, tempNode);
    printIncounts(f, tempNode);

    f->unlinkNode(tempNode);

    delete [] element;
  }

  {
    // --------------------------------------------------------------------
    // TEST #2
    // Build a node for the set: {000, 100}
    int trueNode = f->getTerminalNode(true);
    int node1L1 = f->createTempNodeMaxSize(level1);
    int node2L2 = f->createTempNodeMaxSize(level2);
    int node3L3 = f->createTempNodeMaxSize(level3);

    // node1L1 = {0} --> T
    f->setDownPtr(node1L1, 0, trueNode);

    // node2L2 = {0} --> node1L1
    f->setDownPtr(node2L2, 0, node1L1);

    // node3L3 = {0,1} --> node2L2
    f->setDownPtr(node3L3, 0, node2L2);
    f->setDownPtr(node3L3, 1, node2L2);

    int tempNode = node3L3;
    f->linkNode(tempNode);

    f->unlinkNode(node3L3);
    f->unlinkNode(node2L2);
    f->unlinkNode(node1L1);
    f->unlinkNode(trueNode);

    printf("{000, 100}\n");
    f->showNodeGraph(stdout, tempNode);
    printIncounts(f, tempNode);

    // Add 001
    int* element = new int[nVars + 1];
    element[0] = 0;
    element[1] = 1;
    element[2] = 0;
    element[3] = 0;
    f->accumulate(tempNode, element);

    printf("Added {001}\n");
    f->showNodeGraph(stdout, tempNode);
    printIncounts(f, tempNode);

    // Add 101
    element[1] = 1;
    element[2] = 0;
    element[3] = 1;
    f->accumulate(tempNode, element);

    printf("Added {101}\n");
    f->showNodeGraph(stdout, tempNode);
    printIncounts(f, tempNode);

    f->unlinkNode(tempNode);

    delete [] element;
  }
}


void testB()
{
  int nVars = 3;
  int varSize = 2;

  // Initialize the variable bounds array to provide to the domain

  int* bounds = (int *) malloc(nVars * sizeof(int));
  assert(bounds != 0);
  for (int i = 0; i < nVars; ) { bounds[i++] = varSize; }

  // Create a domain
  domain *d = MEDDLY_createDomain();
  assert(d != 0);
  assert(domain::SUCCESS == d->createVariablesBottomUp(bounds, nVars));

  // Create an MDD forest in this domain (to store states)
  forest* states = d->createForest(true, forest::BOOLEAN,
      forest::MULTI_TERMINAL);
  assert(states != 0);

  mxd_node_manager* f = dynamic_cast<mxd_node_manager*>(states);
  assert(f != 0);

  int level1 = f->getDomain()->getVariableAbove(domain::TERMINALS);
  int level2 = f->getDomain()->getVariableAbove(level1);
  int level3 = f->getDomain()->getVariableAbove(level2);

  {
    // ----------------------------------------------------------------
    // TEST #1
    //
    // Build a node for the set: {000->000, 010->000}
    int trueNode = f->getTerminalNode(true);
    int node1L1P = f->createTempNodeMaxSize(-level1);
    int node2L1 = f->createTempNodeMaxSize(level1);
    int node3L2P = f->createTempNodeMaxSize(-level2);
    int node4L2 = f->createTempNodeMaxSize(level2);
    int node5L3P = f->createTempNodeMaxSize(-level3);
    int node6L3 = f->createTempNodeMaxSize(level3);

    // node1L1P = {0} --> T
    f->setDownPtr(node1L1P, 0, trueNode);

    // node2L1 = {0} --> node1L1P
    f->setDownPtr(node2L1, 0, node1L1P);

    // node3L2P = {0,1} --> node2L1
    f->setDownPtr(node3L2P, 0, node2L1);
    f->setDownPtr(node3L2P, 1, node2L1);

    // node4L2 = {0} --> node3L2P
    f->setDownPtr(node4L2, 1, node3L2P);

    // node5L3P = {0} --> node4L2
    f->setDownPtr(node5L3P, 0, node4L2);

    // node6L3 = {0} --> node5L3P
    f->setDownPtr(node6L3, 0, node5L3P);

    int tempNode = node6L3;
    f->linkNode(tempNode);

    f->unlinkNode(trueNode);
    f->unlinkNode(node1L1P);
    f->unlinkNode(node2L1);
    f->unlinkNode(node3L2P);
    f->unlinkNode(node4L2);
    f->unlinkNode(node5L3P);
    f->unlinkNode(node6L3);

    printf("{010->000, 010->010}\n");
    f->showNodeGraph(stdout, tempNode);
    printIncounts(f, tempNode);

    // Add 011->000
    int* element = new int[nVars + 1];
    int* pelement = new int[nVars + 1];
    element[0] = 0;
    element[1] = 1;
    element[2] = 1;
    element[3] = 0;
    pelement[0] = 0;
    pelement[1] = 0;
    pelement[2] = 0;
    pelement[3] = 0;
    f->accumulate(tempNode, element, pelement);

    printf("Added {011->000} i.e. {001010}\n");
    f->showNodeGraph(stdout, tempNode);
    printIncounts(f, tempNode);

    f->unlinkNode(tempNode);

    delete [] element;
    delete [] pelement;
  }


  // Cleanup; in this case simply delete the domain
  delete d;

  free(bounds);
}


void doPhaseI(expert_forest* f, const int* const* elements,
    int start, int end, int& tempNode)
{
  dd_edge nodeA(f);
  assert(forest::SUCCESS ==
      f->createEdge(elements + start, end - start, nodeA));

  convertDDEdgeToTemporaryNode(nodeA, tempNode);

#ifdef VERBOSE
  printf("Elements %d to %d\n\n", start, end - 1);

  printf("nodeA (createEdge)\n");
  printf("------------------\n");
  f->showNodeGraph(stdout, nodeA.getNode());
  printf("\n");

  printf("tempNode (convert)\n");
  printf("------------------\n");
  f->showNodeGraph(stdout, tempNode);
  printf("\n");

  printf("#########################################################\n\n");
#endif

}


// PhaseII: use accumulate += minterms
void doPhaseII(expert_forest* f, const int* const* elements,
    int start, int end, int& tempNode)
{
  for (int i = start; i < end; ++i) {
    //accumulate(f, tempNode, (int*)(elements[i]));
    f->accumulate(tempNode, (int*)(elements[i]));
  }

#ifdef VERBOSE
  printf("Elements %d to %d\n\n", start, end - 1);

  printf("tempNode (accumulate minterms)\n");
  printf("------------------------------\n");
  f->showNodeGraph(stdout, tempNode);
  printf("\n");

  printf("#########################################################\n\n");
#endif
}


// Phase III: use acccumulate += mdd
void doPhaseIII(expert_forest* f, const int* const* elements,
    int start, int end, int& tempNode)
{
  dd_edge nodeB(f);
  assert(forest::SUCCESS == f->createEdge(elements + start,
        end - start, nodeB));
#if 0
  accumulate(f, tempNode, nodeB.getNode());
#else
  f->accumulate(tempNode, nodeB.getNode());
#endif

#ifdef VERBOSE
  printf("Elements %d to %d\n\n", start, end - 1);

  printf("nodeB (createEdge)\n");
  printf("------------------\n");
  f->showNodeGraph(stdout, nodeB.getNode());
  printf("\n");

  printf("tempNode (accumulate mdd)\n");
  printf("-------------------------\n");
  f->showNodeGraph(stdout, tempNode);
  printf("\n");

  printf("#########################################################\n\n");
#endif

#ifdef DEBUG_ACCUMULATE_MDD
  printf("tempNode: %d\n", tempNode);
  if (f->isActiveNode(tempNode)) {
    printf("Incount(tempNode) before reduction: %d\n\n",
        f->getInCount(tempNode));
  }
  printIncounts(f, tempNode);
#endif

}


void printIncounts(expert_forest* f, int node)
{
  std::set<int> nodes;
  findNodesInGraph(f, node, nodes);
  for (std::set<int>::iterator iter = nodes.begin();
      iter != nodes.end(); ++iter) {
    printf("Incount(%d): %d\n", *iter, f->getInCount(*iter));
  }
}


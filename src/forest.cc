
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



// TODO: Testing

#include "defines.h"
#include <set>
#include <queue>
#include <vector>


// **************************** forest *********************************

forest::forest() {}
forest::~forest() {}


const char* forest::getErrorCodeName(forest::error e)
{
  switch (e) {
    case SUCCESS:
        return "Operation returned successfully";
    case NOT_IMPLEMENTED:
        return "Operation not implemented";
    case INSUFFICIENT_MEMORY:
        return "Operation failed -- lack of memory";
    case INVALID_OPERATION:
        return "Operation not supported for the given forest";
    case INVALID_VARIABLE:
        return "Operation failed -- invalid variable handle";
    case INVALID_ASSIGNMENT:
        return "Operation falied -- variable out of range";
    default:
        return "Unknown error code";
  }
}

// **************************** expert_forest *********************************

expert_forest::expert_forest(domain *d, bool rel, range_type t,
  edge_labeling ev, reduction_rule r, node_storage s, node_deletion_policy nd)
: d(d), isRelation(rel), rangeType(t), edgeLabel(ev),
  reductionRule(r), nodeStorage(s), nodeDeletionPolicy(nd),
  sz(256), firstFree(0), firstHole(-1)
{
  // firstHole < 0 indicates no holes.
  firstHole = -1;
  // Create an array to store pointers to dd_edges.
  edge = (edge_data *) malloc(sz * sizeof(edge_data));
  for (unsigned i = 0; i < sz; ++i) {
    edge[i].nextHole = -1;
    edge[i].edge = 0;
  }
}


void expert_forest::unregisterDDEdges() {
  // Go through the list of valid edges (value > 0), and set
  // the e.index to -1 (indicating unregistered edge).

  // ignore the NULLs; release the rest
  for (unsigned i = 0; i < firstFree; ++i) {
    if (edge[i].edge != 0) {
      DCASSERT(edge[i].nextHole == -1);
      int node = edge[i].edge->getNode();
      unlinkNode(node);
      edge[i].edge->setIndex(-1);
    }
  }

  // firstHole < 0 indicates no holes.
  for (unsigned i = 0; i < firstFree; ++i) {
    edge[i].nextHole = -1;
    edge[i].edge = 0;
  }
  firstHole = -1;
  firstFree = 0;
}


expert_forest::~expert_forest() {
  // Go through the list of valid edges (value > 0), and set
  // the e.index to -1 (indicating unregistered edge).
  // unregisterDDEdges();
  // No need to call this from here -- ~node_manager() calls it.

  // Delete the array.
  free(edge);

  // NOTE: since the user is provided with the dd_edges instances (as opposed
  // to a pointer), the user program will automatically call the
  // destructor for each dd_edge when the corresponding variable goes out of
  // scope. Therefore there is no need to destruct dd_edges from here.

  // Inform the domain that you are going away
  smart_cast<expert_domain*>(d)->unlinkForest(this);
}


void expert_forest::registerEdge(dd_edge& e) {
  // add to collection of edges for this forest.
  // change e.index to help find this edge at a later time.
  if (firstHole >= 0) {
    // hole available; fill it up
    int index = firstHole;
    firstHole = edge[firstHole].nextHole;
    edge[index].edge = &e;
    edge[index].nextHole = -1;
    e.setIndex(index);
  } else {
    // no holes available, add to end of array
    if (firstFree >= sz) {
      // expand edge[]
      int new_sz = sz * 2;
      edge_data* new_edge =
          (edge_data*) realloc(edge, new_sz * sizeof(edge_data));
      if (NULL == new_edge) outOfMemory();
      edge = new_edge;
      for (int i = sz; i < new_sz; ++i)
      {
        edge[i].nextHole = -1;
        edge[i].edge = 0;
      }
      sz = new_sz;
    }
    DCASSERT(firstFree < sz);
    edge[firstFree].nextHole = -1;
    edge[firstFree].edge = &e;
    e.setIndex(firstFree);
    ++firstFree;
  }
}


void expert_forest::unregisterEdge(dd_edge& e) {
  // remove this edge from the collection of edges for this forest.
  // change e.index to -1.
  DCASSERT(e.getIndex() >= 0);
  int index = e.getIndex();
  DCASSERT(edge[index].edge == &e);
  edge[index].edge = 0;
  edge[index].nextHole = firstHole;
  firstHole = index;
  e.setIndex(-1);
}


// TODO: make use of pointers to speed this up.
int expert_forest::getDownPtr(int p, int i) const {
  DCASSERT(isActiveNode(p));
  if (isTerminalNode(p)) return p;
  DCASSERT(i >= 0);
  if (isFullNode(p)) {
    // full or trunc-full node
    if (getFullNodeSize(p) > i) return getFullNodeDownPtr(p, i);

    // full or trunc-full node, but i lies after the last non-zero downpointer
    return 0;
  } else {
    // sparse node
    // binary search to find the index ptr corresponding to i
    int start = 0;
    int stop = getSparseNodeSize(p) - 1;
    
    // the index ptr corresponding to i has been compressed i.e. its a zero.
    if (getSparseNodeIndex(p, start) > i) return 0;
    if (getSparseNodeIndex(p, stop) < i) return 0;

    int mid = (start + stop)/2;
    while(start < stop) {
      if (getSparseNodeIndex(p, mid) == i)
        return getSparseNodeDownPtr(p, mid);
      if (getSparseNodeIndex(p, mid) > i)
        stop = mid - 1;
      else
        start = mid + 1;
      mid = (start + stop)/2;
    }
    if (getSparseNodeIndex(p, mid) == i)
      return getSparseNodeDownPtr(p, mid);
   
    // index ptr not found - it was compressed because it was a zero.
    return 0;
  }
}


bool expert_forest::getDownPtrs(int p, std::vector<int>& dptrs) const {
  if (!isActiveNode(p) || isTerminalNode(p)) return false;

  const int* ptrs = 0;
  assert(getDownPtrs(p, ptrs));

  if (isFullNode(p)) {
    int size = getFullNodeSize(p);
    if (dptrs.size() < unsigned(size)) dptrs.resize(size, 0);
    const int* end = ptrs + size;
    std::vector<int>::iterator iter = dptrs.begin();
    while (ptrs != end) { *iter++ = *ptrs++; }
  }
  else {
    int nnz = getSparseNodeSize(p);
    const int* index = 0;
    assert(getSparseNodeIndexes(p, index));
    const int* end = ptrs + nnz;

    int size = index[nnz-1] + 1;
    if (dptrs.size() < unsigned(size)) dptrs.resize(size, 0);

    while (ptrs != end) { dptrs[*index++] = *ptrs++; }
  }
  return true;
}


bool expert_forest::isStale(int h) const {
  return
    isTerminalNode(h)
    ? discardTemporaryNodesFromComputeCache()
    : isPessimistic()
    ? isZombieNode(h)
    : (getInCount(h) == 0);
}

#ifndef INLINED_REALS

float expert_forest::getReal(int term) const
{
  return (term == 0)? 0.0: *((float*)&(term <<= 1));
}

int expert_forest::getTerminalNode(float a) const
{
  return (a == 0.0)? 0: (*((int*)((void*)&a)) >> 1) | 0x80000000;
}

#endif


unsigned expert_forest::getNodeCount(int p) const
{
  std::set<int> discovered;
  std::queue<int> toExpand;

  if (p != 0 && p != -1) {
    toExpand.push(p);
    discovered.insert(p);
  }

  // expand the front of toExpand;
  // add newly discovered ones to discovered and toExpand

  while (!toExpand.empty()) {
    int p = toExpand.front();
    toExpand.pop();
    if (isTerminalNode(p)) continue;
    // expand
    if (isFullNode(p)) {
      const int sz = getFullNodeSize(p);
      for (int i = 0; i < sz; ++i)
      {
        int temp = getFullNodeDownPtr(p, i);
        if (temp == 0 || temp == -1) continue;
        // insert into discovered and toExpand if new
        if (discovered.find(temp) == discovered.end()) {
          toExpand.push(temp);
          discovered.insert(temp);
        }
      }
    }
    else {
      const int sz = getSparseNodeSize(p);
      for (int i = 0; i < sz; ++i)
      {
        int temp = getSparseNodeDownPtr(p, i);
        if (temp == 0 || temp == -1) continue;
        // insert into discovered and toExpand if new
        if (discovered.find(temp) == discovered.end()) {
          toExpand.push(temp);
          discovered.insert(temp);
        }
      }
    }
  }

  // Add 2 to discovered.size() for terminals 0 and -1.
  return discovered.size() + 2;
}


unsigned expert_forest::getEdgeCount(int p, bool countZeroes) const
{
  std::set<int> discovered;
  std::queue<int> toExpand;

  toExpand.push(p);
  discovered.insert(p);

  unsigned count = 0;

  // expand the front of toExpand;
  // add newly discovered ones to discovered and toExpand

  while (!toExpand.empty()) {
    int p = toExpand.front();
    toExpand.pop();
    if (isTerminalNode(p)) continue;
    // expand
    if (isFullNode(p)) {
      const int sz = getFullNodeSize(p);
      for (int i = 0; i < sz; ++i)
      {
        int temp = getFullNodeDownPtr(p, i);
        if (countZeroes) count++;
        else if (temp != 0) count++;
        if (temp == 0 || temp == -1)  continue;
        // insert into discovered and toExpand if new
        if (discovered.find(temp) == discovered.end()) {
          toExpand.push(temp);
          discovered.insert(temp);
        }
      }
    }
    else {
      const int sz = getSparseNodeSize(p);
      for (int i = 0; i < sz; ++i)
      {
        int temp = getSparseNodeDownPtr(p, i);
        if (countZeroes) count++;
        else if (temp != 0) count++;
        if (temp == 0 || temp == -1)  continue;
        // insert into discovered and toExpand if new
        if (discovered.find(temp) == discovered.end()) {
          toExpand.push(temp);
          discovered.insert(temp);
        }
      }
    }
  }

  return count;
}



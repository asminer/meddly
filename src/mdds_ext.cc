
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



// TODO: Completed implementing ALL functions in mdds_ext.h
// TODO: need to call node_manager:: functions?
//
// TODO: Test createNode variant to build a sparse node for:
//       mdds (done), mxd, evmdd, mtmdd

#include "mdds_ext.h"
#include <algorithm>
#include <limits.h>

//#define DEBUG_SORT_MATRIX
//#define DEBUG_SORT_BUILD

// N = domain->getNumVars() + 1
void sortVector(int** indexes, int* terms, int N, int nVars);
void sortMatrix(int** indexes, int** pindexes, int* terms, int N, int nVars);


// ********************************** MTMDDs **********************************

mtmdd_node_manager::mtmdd_node_manager(domain *d, forest::range_type t)
: node_manager(d, false, t,
      forest::MULTI_TERMINAL, forest::FULLY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{ }


mtmdd_node_manager::~mtmdd_node_manager()
{ }


int mtmdd_node_manager::reduceNode(int p)
{
  DCASSERT(isActiveNode(p));

  if (isReducedNode(p)) return p; 

  DCASSERT(!isTerminalNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(getInCount(p) == 1);

#if 0
  validateIncounts();
#endif

  int size = getFullNodeSize(p);
  int* ptr = getFullNodeDownPtrs(p);
  int node_level = getNodeLevel(p);

  decrTempNodeCount(node_level);

#ifdef DEVELOPMENT_CODE
  int node_height = getNodeHeight(p);
  for (int i=0; i<size; i++) {
    assert(isReducedNode(ptr[i]));
    assert(getNodeHeight(ptr[i]) < node_height);
  }
#endif

  // quick scan: is this node zero?
  int nnz = 0;
  int truncsize = 0;
  for (int i = 0; i < size; ++i) {
    if (0 != ptr[i]) {
      nnz++;
      truncsize = i;
    }
  }
  truncsize++;

  if (0 == nnz) {
    // duplicate of 0
    deleteTempNode(p);
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got 0\n", p);
#endif
    return 0;
  }

  // check for possible reductions
  if (reductionRule == forest::FULLY_REDUCED) {
    if (nnz == getLevelSize(node_level)) {
      int i = 1;
      for ( ; i < size && ptr[i] == ptr[0]; i++);
      if (i == size ) {
        // for all i, ptr[i] == ptr[0]
        int temp = sharedCopy(ptr[0]);
        deleteTempNode(p);
        return temp;
      }
    }
  } else {
    // Quasi-reduced -- no identity reduction for MDDs/MTMDDs
    DCASSERT(reductionRule == forest::QUASI_REDUCED);

    // ensure than all downpointers are pointing to nodes exactly one
    // level below or zero.
    for (int i = 0; i < size; ++i)
    {
      if (ptr[i] == 0) continue;
      if (getNodeLevel(ptr[i]) != (node_level - 1)) {
        int temp = ptr[i];
        ptr[i] = buildQuasiReducedNodeAtLevel(node_level - 1, ptr[i]);
        unlinkNode(temp);
      }
      DCASSERT(ptr[i] == 0 || (getNodeLevel(ptr[i]) == node_level - 1));
    }
  }

  // check unique table
  int q = find(p);
  if (getNull() != q) {
    // duplicate found
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got %d\n", p, q);
#endif
    deleteTempNode(p);
    return sharedCopy(q);
  }

  // insert into unique table
  insert(p);

#ifdef TRACE_REDUCE
  printf("\tReducing %d: unique, compressing\n", p);
#endif

  if (!areSparseNodesEnabled())
    nnz = size;

  // right now, tie goes to truncated full.
  if (2*nnz < truncsize) {
    // sparse is better; convert
    int newoffset = getHole(node_level, 4+2*nnz, true);
    // can't rely on previous ptr, re-point to p
    int* full_ptr = getNodeAddress(p);
    int* sparse_ptr = getAddress(node_level, newoffset);
    // copy first 2 integers: incount, next
    sparse_ptr[0] = full_ptr[0];
    sparse_ptr[1] = full_ptr[1];
    // size
    sparse_ptr[2] = -nnz;
    // copy index into address[]
    sparse_ptr[3 + 2*nnz] = p;
    // get pointers to the new sparse node
    int* indexptr = sparse_ptr + 3;
    int* downptr = indexptr + nnz;
    ptr = full_ptr + 3;
    // copy downpointers
    for (int i=0; i<size; i++, ++ptr) {
      if (*ptr) {
        *indexptr = i;
        *downptr = *ptr;
        ++indexptr;
        ++downptr;
      }
    }
    // trash old node
#ifdef MEMORY_TRACE
    int saved_offset = getNodeOffset(p);
    setNodeOffset(p, newoffset);
    makeHole(node_level, saved_offset, 4 + size);
#else
    makeHole(node_level, getNodeOffset(p), 4 + size);
    setNodeOffset(p, newoffset);
#endif
    // address[p].cache_count does not change
  } else {
    // full is better
    if (truncsize<size) {
      // truncate the trailing 0s
      int newoffset = getHole(node_level, 4+truncsize, true);
      // can't rely on previous ptr, re-point to p
      int* full_ptr = getNodeAddress(p);
      int* trunc_ptr = getAddress(node_level, newoffset);
      // copy first 2 integers: incount, next
      trunc_ptr[0] = full_ptr[0];
      trunc_ptr[1] = full_ptr[1];
      // size
      trunc_ptr[2] = truncsize;
      // copy index into address[]
      trunc_ptr[3 + truncsize] = p;
      // elements
      memcpy(trunc_ptr + 3, full_ptr + 3, truncsize * sizeof(int));
      // trash old node
#ifdef MEMORY_TRACE
      int saved_offset = getNodeOffset(p);
      setNodeOffset(p, newoffset);
      makeHole(node_level, saved_offset, 4 + size);
#else
      makeHole(node_level, getNodeOffset(p), 4 + size);
      setNodeOffset(p, newoffset);
#endif
      // address[p].cache_count does not change
    }
  }

  // address[p].cache_count does not change
  DCASSERT(getCacheCount(p) == 0);
  // Sanity check that the hash value is unchanged
  DCASSERT(find(p) == p);

#if 0
  validateIncounts();
#endif

  return p;
}


int mtmdd_node_manager::createNode(int k, int index, int dptr)
{
  DCASSERT(index >= -1);

#if 1
  if (index > -1 && getLevelSize(k) <= index) {
    //printf("level: %d, curr size: %d, index: %d, ",
    //k, getLevelSize(k), index);
    expertDomain->enlargeVariableBound(k, index + 1);
    //printf("new size: %d\n", getLevelSize(k));
  }
#endif


  if (dptr == 0) return 0;
  if (index == -1) {
    // all downpointers should point to dptr
    if (reductionRule == forest::FULLY_REDUCED) return sharedCopy(dptr);
    int curr = createTempNodeMaxSize(k, false);
    setAllDownPtrsWoUnlink(curr, dptr);
    return reduceNode(curr);
  }

#if 0

  // a single downpointer points to dptr
  int curr = createTempNode(k, index + 1);
  setDownPtrWoUnlink(curr, index, dptr);
  return reduceNode(curr);

#else

  // a single downpointer points to dptr
  if (nodeStorage == FULL_STORAGE ||
      (nodeStorage == FULL_OR_SPARSE_STORAGE && index < 2)) {
    // Build a full node
    int curr = createTempNode(k, index + 1);
    setDownPtrWoUnlink(curr, index, dptr);
    return reduceNode(curr);
  }
  else {
    DCASSERT (nodeStorage == SPARSE_STORAGE ||
        (nodeStorage == FULL_OR_SPARSE_STORAGE && index >= 2));
    // Build a sparse node
    int p = createTempNode(k, 2);
    int* nodeData = getNodeAddress(p);
    // For sparse nodes, size is -ve
    nodeData[2] = -1;
    // indexes followed by downpointers -- here we have one index and one dptr
    nodeData[3] = index;
    nodeData[4] = sharedCopy(dptr);
    // search in unique table
    int q = find(p);
    if (getNull() == q) {
      // no duplicate found; insert into unique table
      insert(p);
      DCASSERT(getCacheCount(p) == 0);
      DCASSERT(find(p) == p);
    }
    else {
      // duplicate found; discard this node and return the duplicate
      // revert to full temp node before discarding
      nodeData[2] = 2;
      nodeData[3] = 0;
      nodeData[4] = 0;
      unlinkNode(dptr);
      deleteTempNode(p);
      p = sharedCopy(q);
    }
    decrTempNodeCount(k);
    return p;
  }

#endif

}


void mtmdd_node_manager::createEdge(const int* v, int term, dd_edge& e)
{
  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;
  DCASSERT(isTerminalNode(term));
  int result = term;
  int curr = 0;
  for (int i=1; i<h_sz; i++) {
    curr = createNode(h2l_map[i], v[h2l_map[i]], result);
    unlinkNode(result);
    result = curr;
  }
  e.set(result, 0, getNodeLevel(result));
  // e.show(stderr, 2);
}


forest::error mtmdd_node_manager::createEdge(const int* const* vlist,
    const int* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || terms == 0 || N <= 0) return forest::INVALID_VARIABLE;

  createEdge(vlist[0], getTerminalNode(terms[0]), e);
  if (N > 1) {
    dd_edge curr(this);
    for (int i=1; i<N; i++) {
      createEdge(vlist[i], getTerminalNode(terms[i]), curr);
      e += curr;
    }
  }

  return forest::SUCCESS;
}


forest::error mtmdd_node_manager::createEdge(int term, dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  term = getTerminalNode(term);

  if (reductionRule == forest::FULLY_REDUCED) {
    e.set(term, 0, domain::TERMINALS);
    return forest::SUCCESS;
  }

  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;
  int result = term;
  int curr = 0;
  for (int i=1; i<h_sz; i++) {
    curr = createTempNodeMaxSize(h2l_map[i], false);
    setAllDownPtrsWoUnlink(curr, result);
    unlinkNode(result);
    result = reduceNode(curr);
  }
  e.set(result, 0, getNodeLevel(result));

  return forest::SUCCESS;
}


forest::error mtmdd_node_manager::evaluate(const dd_edge &f,
    const int* vlist, int &term) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  term = f.getNode();
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  while (!isTerminalNode(term)) {
    term = getDownPtr(term, vlist[h2l_map[getNodeHeight(term)]]);
  }
  term = getInteger(term);
  return forest::SUCCESS;
}


void mtmdd_node_manager::normalizeAndReduceNode(int& p, int& ev)
{
  assert(false);
}


void mtmdd_node_manager::normalizeAndReduceNode(int& p, float& ev)
{
  assert(false);
}


forest::error mtmdd_node_manager::createEdge(const int* const* vlist, int N,
    dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mtmdd_node_manager::createEdge(const int* const* vlist,
    const float* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || terms == 0 || N <= 0) return forest::INVALID_VARIABLE;

  createEdge(vlist[0], getTerminalNode(terms[0]), e);
  if (N > 1) {
    dd_edge curr(this);
    for (int i=1; i<N; i++) {
      createEdge(vlist[i], getTerminalNode(terms[i]), curr);
      e += curr;
    }
  }

  return forest::SUCCESS;
}


forest::error mtmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mtmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mtmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mtmdd_node_manager::createEdge(bool val, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mtmdd_node_manager::createEdge(float term, dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  int termNode = getTerminalNode(term);

  if (reductionRule == forest::FULLY_REDUCED) {
    e.set(termNode, 0, domain::TERMINALS);
    return forest::SUCCESS;
  }

  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;
  int result = termNode;
  int curr = 0;
  for (int i=1; i<h_sz; i++) {
    curr = createTempNodeMaxSize(h2l_map[i], false);
    setAllDownPtrsWoUnlink(curr, result);
    unlinkNode(result);
    result = reduceNode(curr);
  }
  e.set(result, 0, getNodeLevel(result));

  return forest::SUCCESS;
}


forest::error mtmdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mtmdd_node_manager::evaluate(const dd_edge &f,
    const int* vlist, float &term) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  int node = f.getNode();
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  while (!isTerminalNode(node)) {
    node = getDownPtr(node, vlist[h2l_map[getNodeHeight(node)]]);
  }
  term = getReal(node);
  return forest::SUCCESS;
}


forest::error mtmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mtmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mtmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  return forest::INVALID_OPERATION;
}


// *********************************** MXDs *********************************** 


mxd_node_manager::mxd_node_manager(domain *d)
: node_manager(d, true, forest::BOOLEAN,
      forest::MULTI_TERMINAL, forest::IDENTITY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{ }


mxd_node_manager::~mxd_node_manager()
{ }


int mxd_node_manager::reduceNode(int p)
{
  DCASSERT(isActiveNode(p));
  assert(reductionRule == forest::IDENTITY_REDUCED);

  if (isReducedNode(p)) return p; 

  DCASSERT(!isTerminalNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(getInCount(p) == 1);

#if 0
  validateIncounts();
#endif

  int size = getFullNodeSize(p);
  int* ptr = getFullNodeDownPtrs(p);
  int node_level = getNodeLevel(p);

  decrTempNodeCount(node_level);

#ifdef DEVELOPMENT_CODE
  int node_height = getNodeHeight(p);
  if (isUnprimedNode(p)) {
    // unprimed node
    for (int i=0; i<size; i++) {
      assert(isReducedNode(ptr[i]));
      assert(ptr[i] == 0 || (getNodeLevel(ptr[i]) == -node_level));
    }
  } else {
    // primed node
    for (int i=0; i<size; i++) {
      assert(isReducedNode(ptr[i]));
      assert(getNodeHeight(ptr[i]) < node_height);
      assert(isTerminalNode(ptr[i]) || isUnprimedNode(ptr[i]));
    }
  }
#endif

  // quick scan: is this node zero?
  int nnz = 0;
  int truncsize = 0;
  for (int i = 0; i < size; ++i) {
    if (0 != ptr[i]) {
      nnz++;
      truncsize = i;
    }
  }
  truncsize++;

  if (0 == nnz) {
    // duplicate of 0
    deleteTempNode(p);
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got 0\n", p);
#endif
    return 0;
  }

  // check for possible reductions
  if (isUnprimedNode(p) && nnz == getLevelSize(node_level)) {
    // check for identity matrix node
    int temp = getDownPtr(ptr[0], 0);
    if (temp != 0) {
      bool ident = true;
      for (int i = 0; i < size && ident; i++) {
        if (!singleNonZeroAt(ptr[i], temp, i)) ident = false;
      }
      if (ident) {
        // passed all tests for identity matrix node
        temp = sharedCopy(temp);
        deleteTempNode(p);
        return temp;
      }
    }
  }

  // check unique table
  int q = find(p);
  if (getNull() != q) {
    // duplicate found
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got %d\n", p, q);
#endif
    deleteTempNode(p);
    return sharedCopy(q);
  }

  // insert into unique table
  insert(p);

#ifdef TRACE_REDUCE
  printf("\tReducing %d: unique, compressing\n", p);
#endif

  if (!areSparseNodesEnabled())
    nnz = size;

  // right now, tie goes to truncated full.
  if (2*nnz < truncsize) {
    // sparse is better; convert
    int newoffset = getHole(node_level, 4+2*nnz, true);
    // can't rely on previous ptr, re-point to p
    int* full_ptr = getNodeAddress(p);
    int* sparse_ptr = getAddress(node_level, newoffset);
    // copy first 2 integers: incount, next
    sparse_ptr[0] = full_ptr[0];
    sparse_ptr[1] = full_ptr[1];
    // size
    sparse_ptr[2] = -nnz;
    // copy index into address[]
    sparse_ptr[3 + 2*nnz] = p;
    // get pointers to the new sparse node
    int* indexptr = sparse_ptr + 3;
    int* downptr = indexptr + nnz;
    ptr = full_ptr + 3;
    // copy downpointers
    for (int i=0; i<size; i++, ++ptr) {
      if (*ptr) {
        *indexptr = i;
        *downptr = *ptr;
        ++indexptr;
        ++downptr;
      }
    }
    // trash old node
#ifdef MEMORY_TRACE
    int saved_offset = getNodeOffset(p);
    setNodeOffset(p, newoffset);
    makeHole(node_level, saved_offset, 4 + size);
#else
    makeHole(node_level, getNodeOffset(p), 4 + size);
    setNodeOffset(p, newoffset);
#endif
    // address[p].cache_count does not change
  } else {
    // full is better
    if (truncsize<size) {
      // truncate the trailing 0s
      int newoffset = getHole(node_level, 4+truncsize, true);
      // can't rely on previous ptr, re-point to p
      int* full_ptr = getNodeAddress(p);
      int* trunc_ptr = getAddress(node_level, newoffset);
      // copy first 2 integers: incount, next
      trunc_ptr[0] = full_ptr[0];
      trunc_ptr[1] = full_ptr[1];
      // size
      trunc_ptr[2] = truncsize;
      // copy index into address[]
      trunc_ptr[3 + truncsize] = p;
      // elements
      memcpy(trunc_ptr + 3, full_ptr + 3, truncsize * sizeof(int));
      // trash old node
#ifdef MEMORY_TRACE
      int saved_offset = getNodeOffset(p);
      setNodeOffset(p, newoffset);
      makeHole(node_level, saved_offset, 4 + size);
#else
      makeHole(node_level, getNodeOffset(p), 4 + size);
      setNodeOffset(p, newoffset);
#endif
      // address[p].cache_count does not change
    }
  }

  // address[p].cache_count does not change
  DCASSERT(getCacheCount(p) == 0);
  // Sanity check that the hash value is unchanged
  DCASSERT(find(p) == p);

#if 0
  validateIncounts();
#endif

  return p;
}


int mxd_node_manager::createNode(int k, int index, int dptr)
{
  DCASSERT(index >= 0 && dptr >= -1);

  if (dptr == 0) return 0;

  // a single downpointer points to dptr
  if (nodeStorage == FULL_STORAGE ||
      (nodeStorage == FULL_OR_SPARSE_STORAGE && index < 2)) {
    // Build a full node
    int curr = createTempNode(k, index + 1);
    setDownPtrWoUnlink(curr, index, dptr);
    return reduceNode(curr);
  }
  else {
    DCASSERT (nodeStorage == SPARSE_STORAGE ||
        (nodeStorage == FULL_OR_SPARSE_STORAGE && index >= 2));
    // Build a sparse node
    int p = createTempNode(k, 2);
    int* nodeData = getNodeAddress(p);
    // For sparse nodes, size is -ve
    nodeData[2] = -1;
    // indexes followed by downpointers -- here we have one index and one dptr
    nodeData[3] = index;
    nodeData[4] = sharedCopy(dptr);
    // search in unique table
    int q = find(p);
    if (getNull() == q) {
      // no duplicate found; insert into unique table
      insert(p);
      DCASSERT(getCacheCount(p) == 0);
      DCASSERT(find(p) == p);
    }
    else {
      // duplicate found; discard this node and return the duplicate
      // revert to full temp node before discarding
      nodeData[2] = 2;
      nodeData[3] = 0;
      nodeData[4] = 0;
      unlinkNode(dptr);
      deleteTempNode(p);
      p = sharedCopy(q);
    }
    decrTempNodeCount(k);
    return p;
  }
}

int mxd_node_manager::createNode(int k, int index1, int index2, int dptr)
{
  DCASSERT(reductionRule == forest::IDENTITY_REDUCED);

  DCASSERT((index1 >= 0 && index2 >= 0) ||
      (index1 >= -1 && index2 >= -1) ||
      (index1 >= -2 && index2 >= -2 && index1 == index2));


  if (index1 == -2) {
    // "don't change"
    DCASSERT(index2 == -2);
    return sharedCopy(dptr);
  }

  int p = 0;
  if (index2 == -1) {
    // represents "don't care"
    p = createTempNodeMaxSize(-k, false);
    setAllDownPtrsWoUnlink(p, dptr);
    p = reduceNode(p);
  } else {
    p = createNode(-k, index2, dptr);
  }

  int curr = 0;
  if (index1 == -1) {
    // represents "don't care"
    curr = createTempNodeMaxSize(k, false);
    setAllDownPtrsWoUnlink(curr, p);
    curr = reduceNode(curr);
  } else {
    curr = createNode(k, index1, p);
  }

  unlinkNode(p);
  return curr;
}


void mxd_node_manager::createEdge(const int* v, const int* vp,  dd_edge& e)
{
  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;
  int curr = getTerminalNode(true);
  int prev = curr;
  for (int i=1; i<h_sz; i++) {
    prev = curr;
    curr = createNode(h2l_map[i], v[h2l_map[i]], vp[h2l_map[i]], prev);
    unlinkNode(prev);
  }
  e.set(curr, 0, getNodeLevel(curr));
}


forest::error mxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || vplist == 0 || N <= 0) return forest::INVALID_VARIABLE;

  createEdge(vlist[0], vplist[0], e);
  if (N > 1) {
    dd_edge curr(this);
    for (int i=1; i<N; i++) {
      createEdge(vlist[i], vplist[i], curr);
      e += curr;
    }
  }

  return forest::SUCCESS;
}


forest::error mxd_node_manager::createEdge(bool val, dd_edge &e)
{
  if (!val) {
    e.set(getTerminalNode(false), 0, domain::TERMINALS);
  }
  else {
    DCASSERT(val);
    // construct the edge bottom-up
    const int* h2l_map = expertDomain->getHeightsToLevelsMap();
    int h_sz = expertDomain->getNumVariables() + 1;
    int curr = getTerminalNode(true);
    int prev = curr;
    for (int i=1; i<h_sz; i++) {
      prev = curr;
      curr = createNode(h2l_map[i], -1, -1, prev);
      unlinkNode(prev);
    }
    e.set(curr, 0, getNodeLevel(curr));
  }
  return forest::SUCCESS;
}


forest::error mxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  int node = f.getNode();
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int level = 0;
  while (!isTerminalNode(node)) {
    level = h2l_map[getNodeHeight(node)];
    node = getDownPtr(node, vlist[level]);
    node = getDownPtr(node, vplist[level]);
  }
  term = getBoolean(node);
  return forest::SUCCESS;
}


forest::error
mxd_node_manager::
findFirstElement(const dd_edge& f, int* vlist, int* vplist) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  // vlist is based on level handles.
  int node = f.getNode();
  if (node == 0) return forest::INVALID_ASSIGNMENT;

  int currLevel = expertDomain->getTopVariable();
  DCASSERT(currLevel != domain::TERMINALS);
  while (currLevel != domain::TERMINALS)
  {
    DCASSERT(node != 0);
    DCASSERT(isUnprimedNode(node));
    if (currLevel != getNodeLevel(node)) {
      // currLevel is "higher" than node, and has been skipped.
      // Since this is a mxd, reduced nodes enable "don't change" paths
      // at the skipped level.
      vlist[currLevel] = 0;   // picking the first index
      vplist[currLevel] = 0;
    } else {
      // find a valid path at this unprime level
      if (isFullNode(node)) {
        int size = getFullNodeSize(node);
        for (int i = 0; i < size; i++)
        {
          int n = getFullNodeDownPtr(node, i);
          if (n != 0) {
            node = n;
            vlist[currLevel] = i;
            break;
          }
        }
      } else {
        vlist[currLevel] = getSparseNodeIndex(node, 0);
        node = getSparseNodeDownPtr(node, 0);
      }
      
      DCASSERT(!isTerminalNode(node));
      // can't be -1 because that violates MXD properties
      // can't be 0 because node cannot be set to 0 in the above construct.
      DCASSERT(isPrimedNode(node));
      // find a valid path at this prime level
      if (isFullNode(node)) {
        int size = getFullNodeSize(node);
        for (int i = 0; i < size; i++)
        {
          int n = getFullNodeDownPtr(node, i);
          if (n != 0) {
            node = n;
            vplist[currLevel] = i;
            break;
          }
        }
      } else {
        vplist[currLevel] = getSparseNodeIndex(node, 0);
        node = getSparseNodeDownPtr(node, 0);
      }
    }
    currLevel = expertDomain->getVariableBelow(currLevel);
  }

  return forest::SUCCESS;
}


void mxd_node_manager::normalizeAndReduceNode(int& p, int& ev)
{
  assert(false);
}


void mxd_node_manager::normalizeAndReduceNode(int& p, float& ev)
{
  assert(false);
}


forest::error mxd_node_manager::createEdge(const int* const* vlist, int N,
    dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::createEdge(const int* const* vlist,
    const int* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::createEdge(const int* const* vlist,
    const float* terms, int n, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::createEdge(int val, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::createEdge(float val, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    int &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    float &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  return forest::INVALID_OPERATION;
}


// ******************************** MTMXDs ******************************* 


mtmxd_node_manager::mtmxd_node_manager(domain *d, forest::range_type t)
: node_manager(d, true, t,
      forest::MULTI_TERMINAL, forest::IDENTITY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{ }


mtmxd_node_manager::~mtmxd_node_manager()
{ }


int mtmxd_node_manager::reduceNode(int p)
{
  DCASSERT(isActiveNode(p));
  assert(reductionRule == forest::IDENTITY_REDUCED);

  if (isReducedNode(p)) return p; 

  DCASSERT(!isTerminalNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(getInCount(p) == 1);

#if 0
  validateIncounts();
#endif

  int size = getFullNodeSize(p);
  int* ptr = getFullNodeDownPtrs(p);
  int node_level = getNodeLevel(p);

  decrTempNodeCount(node_level);

#ifdef DEVELOPMENT_CODE
  int node_height = getNodeHeight(p);
  if (isUnprimedNode(p)) {
    // unprimed node
    for (int i=0; i<size; i++) {
      assert(isReducedNode(ptr[i]));
      assert(ptr[i] == 0 || (getNodeLevel(ptr[i]) == -node_level));
    }
  } else {
    // primed node
    for (int i=0; i<size; i++) {
      assert(isReducedNode(ptr[i]));
      assert(getNodeHeight(ptr[i]) < node_height);
      assert(isTerminalNode(ptr[i]) || isUnprimedNode(ptr[i]));
    }
  }
#endif

  // quick scan: is this node zero?
  int nnz = 0;
  int truncsize = 0;
  for (int i = 0; i < size; ++i) {
    if (0 != ptr[i]) {
      nnz++;
      truncsize = i;
    }
  }
  truncsize++;

  if (0 == nnz) {
    // duplicate of 0
    deleteTempNode(p);
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got 0\n", p);
#endif
    return 0;
  }

  // check for possible reductions
  if (isUnprimedNode(p) && nnz == getLevelSize(node_level)) {
    // check for identity matrix node
    int temp = getDownPtr(ptr[0], 0);
    if (temp != 0) {
      bool ident = true;
      for (int i = 0; i < size && ident; i++) {
        if (!singleNonZeroAt(ptr[i], temp, i)) ident = false;
      }
      if (ident) {
        // passed all tests for identity matrix node
        temp = sharedCopy(temp);
        deleteTempNode(p);
        return temp;
      }
    }
  }

  // check unique table
  int q = find(p);
  if (getNull() != q) {
    // duplicate found
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got %d\n", p, q);
#endif
    deleteTempNode(p);
    return sharedCopy(q);
  }

  // insert into unique table
  insert(p);

#ifdef TRACE_REDUCE
  printf("\tReducing %d: unique, compressing\n", p);
#endif

  if (!areSparseNodesEnabled())
    nnz = size;

  // right now, tie goes to truncated full.
  if (2*nnz < truncsize) {
    // sparse is better; convert
    int newoffset = getHole(node_level, 4+2*nnz, true);
    // can't rely on previous ptr, re-point to p
    int* full_ptr = getNodeAddress(p);
    int* sparse_ptr = getAddress(node_level, newoffset);
    // copy first 2 integers: incount, next
    sparse_ptr[0] = full_ptr[0];
    sparse_ptr[1] = full_ptr[1];
    // size
    sparse_ptr[2] = -nnz;
    // copy index into address[]
    sparse_ptr[3 + 2*nnz] = p;
    // get pointers to the new sparse node
    int* indexptr = sparse_ptr + 3;
    int* downptr = indexptr + nnz;
    ptr = full_ptr + 3;
    // copy downpointers
    for (int i=0; i<size; i++, ++ptr) {
      if (*ptr) {
        *indexptr = i;
        *downptr = *ptr;
        ++indexptr;
        ++downptr;
      }
    }
    // trash old node
#ifdef MEMORY_TRACE
    int saved_offset = getNodeOffset(p);
    setNodeOffset(p, newoffset);
    makeHole(node_level, saved_offset, 4 + size);
#else
    makeHole(node_level, getNodeOffset(p), 4 + size);
    setNodeOffset(p, newoffset);
#endif
    // address[p].cache_count does not change
  } else {
    // full is better
    if (truncsize<size) {
      // truncate the trailing 0s
      int newoffset = getHole(node_level, 4+truncsize, true);
      // can't rely on previous ptr, re-point to p
      int* full_ptr = getNodeAddress(p);
      int* trunc_ptr = getAddress(node_level, newoffset);
      // copy first 2 integers: incount, next
      trunc_ptr[0] = full_ptr[0];
      trunc_ptr[1] = full_ptr[1];
      // size
      trunc_ptr[2] = truncsize;
      // copy index into address[]
      trunc_ptr[3 + truncsize] = p;
      // elements
      memcpy(trunc_ptr + 3, full_ptr + 3, truncsize * sizeof(int));
      // trash old node
#ifdef MEMORY_TRACE
      int saved_offset = getNodeOffset(p);
      setNodeOffset(p, newoffset);
      makeHole(node_level, saved_offset, 4 + size);
#else
      makeHole(node_level, getNodeOffset(p), 4 + size);
      setNodeOffset(p, newoffset);
#endif
      // address[p].cache_count does not change
    }
  }

  // address[p].cache_count does not change
  DCASSERT(getCacheCount(p) == 0);
  // Sanity check that the hash value is unchanged
  DCASSERT(find(p) == p);

#if 0
  validateIncounts();
#endif

  return p;
}


int mtmxd_node_manager::createNode(int k, int index, int dptr)
{
  DCASSERT(index >= 0 && index < getLevelSize(k) && isValidNodeIndex(dptr));

  if (dptr == 0) return 0;

  // a single downpointer points to dptr
  if (nodeStorage == FULL_STORAGE ||
      (nodeStorage == FULL_OR_SPARSE_STORAGE && index < 2)) {
    // Build a full node
    int curr = createTempNode(k, index + 1);
    setDownPtrWoUnlink(curr, index, dptr);
    return reduceNode(curr);
  }
  else {
    DCASSERT (nodeStorage == SPARSE_STORAGE ||
        (nodeStorage == FULL_OR_SPARSE_STORAGE && index >= 2));
    // Build a sparse node
    int p = createTempNode(k, 2);
    int* nodeData = getNodeAddress(p);
    // For sparse nodes, size is -ve
    nodeData[2] = -1;
    // indexes followed by downpointers -- here we have one index and one dptr
    nodeData[3] = index;
    nodeData[4] = sharedCopy(dptr);
    // search in unique table
    int q = find(p);
    if (getNull() == q) {
      // no duplicate found; insert into unique table
      insert(p);
      DCASSERT(getCacheCount(p) == 0);
      DCASSERT(find(p) == p);
    }
    else {
      // duplicate found; discard this node and return the duplicate
      // revert to full temp node before discarding
      nodeData[2] = 2;
      nodeData[3] = 0;
      nodeData[4] = 0;
      unlinkNode(dptr);
      deleteTempNode(p);
      p = sharedCopy(q);
    }
    decrTempNodeCount(k);
    return p;
  }
}

int mtmxd_node_manager::createNode(int k, int index1, int index2, int dptr)
{
  DCASSERT(reductionRule == forest::IDENTITY_REDUCED);

  DCASSERT((index1 >= 0 && index2 >= 0) ||
      (index1 >= -1 && index2 >= -1) ||
      (index1 >= -2 && index2 >= -2 && index1 == index2));



  if (index1 == -2) {
    // "don't change"
    return sharedCopy(dptr);
  }

  int p = 0;
  if (index2 == -1) {
    // represents "don't care"
    p = createTempNodeMaxSize(-k, false);
    setAllDownPtrsWoUnlink(p, dptr);
    p = reduceNode(p);
  } else {
    p = createNode(-k, index2, dptr);
  }

  int curr = 0;
  if (index1 == -1) {
    // represents "don't care"
    curr = createTempNodeMaxSize(k, false);
    setAllDownPtrsWoUnlink(curr, p);
    curr = reduceNode(curr);
  } else {
    curr = createNode(k, index1, p);
  }

  unlinkNode(p);
  return curr;
}


void mtmxd_node_manager::createEdge(const int* v, const int* vp, int term,
    dd_edge& e)
{
  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;
  DCASSERT(isTerminalNode(term));
  int curr = term;
  int prev = curr;
#if 1
  const int* end = h2l_map + h_sz;
  for (++h2l_map; h2l_map != end; ++h2l_map)
  {
    prev = curr;
    curr = createNode(*h2l_map, v[*h2l_map], vp[*h2l_map], prev);
    unlinkNode(prev);
  }
#else
  for (int i = 1; i < h_sz; ++i)
  {
    prev = curr;
    curr = createNode(h2l_map[i], v[h2l_map[i]], vp[h2l_map[i]], prev);
    unlinkNode(prev);
  }
#endif
  e.set(curr, 0, getNodeLevel(curr));
}


forest::error mtmxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge& e)
{
  DCASSERT(getRangeType() == forest::INTEGER);
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || vplist == 0 || terms == 0 || N <= 0)
    return forest::INVALID_VARIABLE;

  createEdge(vlist[0], vplist[0], getTerminalNode(terms[0]), e);
  if (N > 1) {
    dd_edge curr(this);
    for (int i=1; i<N; i++) {
      createEdge(vlist[i], vplist[i], getTerminalNode(terms[i]), curr);
      e += curr;
    }
  }

  return forest::SUCCESS;
}


forest::error mtmxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge& e)
{
  DCASSERT(getRangeType() == forest::REAL);
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || vplist == 0 || terms == 0 || N <= 0)
    return forest::INVALID_VARIABLE;

#ifndef SORT_BUILD
  createEdge(vlist[0], vplist[0], getTerminalNode(terms[0]), e);
  if (N > 1) {
    dd_edge curr(this);
    for (int i=1; i<N; i++) {
      createEdge(vlist[i], vplist[i], getTerminalNode(terms[i]), curr);
      e += curr;
    }
  }
#else
  if (N < 2) {
    createEdge(vlist[0], vplist[0], getTerminalNode(terms[0]), e);
  } else {
    // check for special cases
    bool specialCasesFound = false;
    for (int i = 0; i < N; i++)
    {
      int* curr = (int*)vlist[i];
      int* end = curr + expertDomain->getNumVariables() + 1;
      for ( ; curr != end; ) { if (*curr++ < 0) break; }
      if (curr != end) { specialCasesFound = true; break; }

      curr = (int*)vplist[i];
      end = curr + expertDomain->getNumVariables() + 1;
      for ( ; curr != end; ) { if (*curr++ < 0) break; }
      if (curr != end) { specialCasesFound = true; break; }
    }

    if (specialCasesFound) {
#if 0
      assert(N > 1);
      createEdge(vlist[0], vplist[0], getTerminalNode(terms[0]), e);
      dd_edge curr(this);
      for (int i=1; i<N; i++) {
        createEdge(vlist[i], vplist[i], getTerminalNode(terms[i]), curr);
        e += curr;
      }
#else
      // populate curr[]
      vector<dd_edge*> curr(N, (dd_edge*)0);
      for (vector<dd_edge*>::iterator iter = curr.begin();
          iter != curr.end(); ++iter, ++vlist, ++vplist, ++terms)
      {
        *iter = new dd_edge(this);
        createEdge(*vlist, *vplist, getTerminalNode(*terms), *(*iter));
      }

      // Iteratively halve curr[] by combining adjacent locations.
      // When curr is of size 1 curr[0] gives the result

      for (vector<dd_edge*> next; curr.size() > 1; curr = next)
      {
        // build next[] and then update curr[]
        next.resize((curr.size()+1)/2, (dd_edge*)0);

        vector<dd_edge*>::iterator currIter = curr.begin();
        vector<dd_edge*>::iterator nextIter = next.begin();
        for ( ; currIter != curr.end(); currIter += 2, ++nextIter)
        {
          *nextIter = *currIter;
        }

        currIter = curr.begin() + 1;
        nextIter = next.begin();
        for ( ; currIter != curr.end(); currIter += 2, ++nextIter)
        {
          *(*nextIter) += *(*currIter);
          delete *currIter;
        }
      }
      assert(curr.size() == 1);
      e = *(curr[0]);
      delete curr[0];
#endif
    }
    else {
      int* list[N];
      int* plist[N];
      int tlist[N];
      memcpy(list, vlist, N * sizeof(int*));
      memcpy(plist, vplist, N * sizeof(int*));
      int* curr = tlist;
      const float* currTerm = terms;
      for (int* last = tlist + N; curr != last; ++curr, ++currTerm)
        *curr = getTerminalNode(*currTerm);
      sortMatrix(list, plist, tlist, N, getDomain()->getNumVariables() + 1);

      int result = 
        sortBuild(list, plist, tlist, getDomain()->getNumVariables(), 0, N);
      e.set(result, 0, getNodeLevel(result));

      memset(list, 0, N * sizeof(int*));
      memset(plist, 0, N * sizeof(int*));
    }
  }
#endif
  return forest::SUCCESS;
}

#ifdef SORT_BUILD
int mtmxd_node_manager::sortBuild(int** list, int** plist, int* tlist,
    int height, int begin, int end)
{
  // [begin, end)

  // terminal condition
  if (height == 0)
  {
    assert(begin + 1 == end);
    return tlist[begin];
  }

  int** currList = 0;
  int nextHeight = 0;
  int level = 0;
  int absLevel = 0;

  if (height > 0) {
    currList = list;
    nextHeight = -height;
    level = expertDomain->getVariableWithHeight(height);
  } else {
    currList = plist;
    nextHeight = -height-1;
    level = -(expertDomain->getVariableWithHeight(-height));
  }
  absLevel = level < 0? -level: level;

  vector<int> nodes;
  int currBegin = begin;
  int currEnd = currBegin;
  for ( ; currEnd < end; currBegin = currEnd) 
  {
    int currIndex = currList[currBegin][absLevel];
    assert(currIndex >= 0);
    for (currEnd = currBegin + 1;
        currEnd < end && currIndex == currList[currEnd][absLevel];
        ++currEnd);
    // found new range
    // to be "unioned" and assigned to result[currIndex]
#ifdef DEBUG_SORT_BUILD
    printf("level: %d, currIndex: %d, currBegin: %d, currEnd: %d\n",
        level, currIndex, currBegin, currEnd);
    fflush(stdout);
#endif
    int node = sortBuild(list, plist, tlist, nextHeight, currBegin, currEnd);
#ifdef DEBUG_SORT_BUILD
    printf("level: %d, currIndex: %d, currBegin: %d, currEnd: %d\n",
        level, currIndex, currBegin, currEnd);
    fflush(stdout);
#endif
    nodes.push_back(currIndex);
    nodes.push_back(node);
  }

  assert(nodes.size() >= 2);
  if (nodes.size() == 2) {
    // single entry: store directly as sparse
    // TODO: this should be able to accept many more cases
    int result = createNode(level, nodes[0], nodes[1]);
    unlinkNode(nodes[1]);
    return result;
  }
  // full node
  int size = nodes[nodes.size() - 2] + 1;
  int result = createTempNode(level, size);
  for (vector<int>::iterator iter = nodes.begin();
      iter != nodes.end(); iter += 2)
  {
    if (getDownPtr(result, *iter) != 0) {
      for (unsigned i = 0u; i < nodes.size(); i+=2)
      {
        printf("%d:%d ", nodes[i], nodes[i+1]);
      }
      printf("\n");
      for (int i = begin; i < end; i++)
      {
        printf("%d:%d ", i, currList[i][absLevel]);
      }
      printf("\n");
      exit(1);
    }
    setDownPtrWoUnlink(result, *iter, *(iter+1));
    unlinkNode(*(iter+1));
  }
  return reduceNode(result);
}
#endif

int mtmxd_node_manager::createEdge(int dptr)
{
  DCASSERT(isTerminalNode(dptr));
  if (dptr == 0) return dptr;

  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;
  int curr = dptr;
  int prev = 0;
  for (int i=1; i<h_sz; i++) {
    prev = curr;
    curr = createNode(h2l_map[i], -1, -1, prev);
    unlinkNode(prev);
  }
  return curr;
}


forest::error mtmxd_node_manager::createEdge(int val, dd_edge &e)
{
  DCASSERT(getRangeType() == forest::INTEGER);
  if (e.getForest() != this) return forest::INVALID_OPERATION;

  int node = createEdge(getTerminalNode(val));
  e.set(node, 0, getNodeLevel(node));
  return forest::SUCCESS;
}


forest::error mtmxd_node_manager::createEdge(float val, dd_edge &e)
{
  DCASSERT(getRangeType() == forest::REAL);
  if (e.getForest() != this) return forest::INVALID_OPERATION;

  int node = createEdge(getTerminalNode(val));
  e.set(node, 0, getNodeLevel(node));
  return forest::SUCCESS;
}


forest::error mtmxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  DCASSERT(getRangeType() == forest::INTEGER);
  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  int node = f.getNode();
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int level = 0;
  while (!isTerminalNode(node)) {
    level = h2l_map[getNodeHeight(node)];
    node = getDownPtr(node, vlist[level]);
    node = getDownPtr(node, vplist[level]);
  }
  term = getInteger(node);
  return forest::SUCCESS;
}


forest::error mtmxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  DCASSERT(getRangeType() == forest::REAL);
  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  int node = f.getNode();
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int level = 0;
  while (!isTerminalNode(node)) {
    level = h2l_map[getNodeHeight(node)];
    node = getDownPtr(node, vlist[level]);
    node = getDownPtr(node, vplist[level]);
  }
  term = getReal(node);
  return forest::SUCCESS;
}


void mtmxd_node_manager::normalizeAndReduceNode(int& p, int& ev)
{
  assert(false);
}


void mtmxd_node_manager::normalizeAndReduceNode(int& p, float& ev)
{
  assert(false);
}


forest::error mtmxd_node_manager::createEdge(const int* const* vlist, int N,
    dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mtmxd_node_manager::createEdge(const int* const* vlist,
    const int* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mtmxd_node_manager::createEdge(const int* const* vlist,
    const float* terms, int n, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mtmxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mtmxd_node_manager::createEdge(bool val, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mtmxd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mtmxd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    int &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mtmxd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    float &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mtmxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  return forest::INVALID_OPERATION;
}


// ********************************** EVMDDs ********************************** 


evmdd_node_manager::evmdd_node_manager(domain *d, forest::edge_labeling el)
: node_manager(d, false, forest::BOOLEAN,
      el, forest::FULLY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{
}


evmdd_node_manager::~evmdd_node_manager()
{ }

int evmdd_node_manager::createTempNode(int k, int sz, bool clear)
{
  DCASSERT(k != 0);

  if (isTimeToGc()) { gc(); }

  CHECK_RANGE(1, mapLevel(k), l_size);
  DCASSERT(level[mapLevel(k)].data != NULL);
  CHECK_RANGE(1, sz, getLevelSize(k) + 1);

  // get a location in address[] to store the node
  int p = getFreeNode(k);

#ifdef DEBUG_MDD_H
  printf("%s: k: %d, sz: %d, new p: %d\n", __func__, k, sz, p);
  fflush(stdout);
#endif

  address[p].level = k;
  address[p].offset = getHole(k, 4 + 2 * sz, true);
  address[p].cache_count = 0;

#ifdef DEBUG_MDD_H
  printf("%s: offset: %d\n", __func__, address[p].offset);
  fflush(stdout);
#endif

  int* foo = level[mapLevel(k)].data + address[p].offset;
  foo[0] = 1;                   // #incoming
  foo[1] = temp_node;
  foo[2] = sz;                  // size
  foo[3 + sz + sz] = p;         // pointer to this node in the address array

  // initialize
  if (clear) {
    initDownPtrs(p);
    initEdgeValues(p);
  }

#ifdef TRACK_DELETIONS
  cout << "Creating node " << p << "\n";
  cout.flush();
#endif

  incrTempNodeCount(k);
  nodes_activated_since_gc++;

  return p;
}


// returns index with a[]; -1 if not found
int binarySearch(const int* a, int sz, int find)
{
  const int* begin = a;
  const int* end = a + sz;
  if (find < *begin || *(end - 1) < find) return -1;
  while (begin != end) {
    const int* mid = begin + (end - begin) / 2;
    if (*mid == find) return (mid - a);
    if (*mid < find) {
      begin = mid;
    } else {
      end = mid;
    }
  }
  return (*begin == find)? begin - a: -1;
}


int evmdd_node_manager::reduceNode(int p)
{
  assert(false);
  return 0;
}


forest::error
evmdd_node_manager::
createEdge(const int* const* vlist, const int* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error
evmdd_node_manager::
createEdge(int val, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error
evmdd_node_manager::
evaluate(const dd_edge &f, const int* vlist, int &term) const
{
  return forest::INVALID_OPERATION;
}



forest::error evmdd_node_manager::createEdge(const int* const* vlist, int N,
    dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::createEdge(const int* const* vlist,
    const float* terms, int n, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::createEdge(bool val, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::createEdge(float val, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    float &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error evmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  return forest::INVALID_OPERATION;
}


// ********************************* EV+MDDs ********************************** 

evplusmdd_node_manager::evplusmdd_node_manager(domain *d)
: evmdd_node_manager(d, forest::EVPLUS)
{ }

evplusmdd_node_manager::~evplusmdd_node_manager()
{ }

void evplusmdd_node_manager::initEdgeValues(int p) {
  DCASSERT(!isReducedNode(p));
  DCASSERT(isFullNode(p));
  int *edgeptr = getFullNodeEdgeValues(p);
  int *last = edgeptr + getFullNodeSize(p);
  for ( ; edgeptr != last; ++edgeptr) *edgeptr = INF;
}


// Similar to getDownPtrs() but for EV+MDDs
bool evplusmdd_node_manager::getDownPtrsAndEdgeValues(int p,
    std::vector<int>& dptrs, std::vector<int>& evs) const
{
  if (!isActiveNode(p) || isTerminalNode(p) || !isReducedNode(p))
    return false;

  if (isFullNode(p)) {
    int size = getFullNodeSize(p);
    if (dptrs.size() < unsigned(size)) {
      dptrs.resize(size, 0);
      evs.resize(size, INF);
    }
    const int* ptrs = getFullNodeDownPtrsReadOnly(p);
    const int* end = ptrs + size;
    const int* evptrs = getFullNodeEdgeValuesReadOnly(p);
    std::vector<int>::iterator iter = dptrs.begin();
    std::vector<int>::iterator eviter = evs.begin();
    while (ptrs != end)
    {
      *iter++ = *ptrs++;
      *eviter++ = *evptrs++;
    }
  }
  else {
    int nnz = getSparseNodeSize(p);
    int size = getLargestIndex(p) + 1;
    if (dptrs.size() < unsigned(size)) {
      dptrs.resize(size, 0);
      evs.resize(size, INF);
    }
    const int* ptrs = getSparseNodeDownPtrs(p);
    const int* index = getSparseNodeIndexes(p);
    const int* evptrs = getSparseNodeEdgeValues(p);
    const int* end = ptrs + nnz;
    while (ptrs != end)
    {
      evs[*index] = *evptrs++;
      dptrs[*index++] = *ptrs++;
    }
  }
  return true;
}


void evplusmdd_node_manager::normalizeAndReduceNode(int& p, int& ev)
{
  DCASSERT(isActiveNode(p));

  if (isReducedNode(p)) {
    if (p == 0) ev = INF;
    return;
  }

  DCASSERT(!isTerminalNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(getInCount(p) == 1);

#if 0
  validateIncounts();
#endif

  const int size = getFullNodeSize(p);
  int *dptr = getFullNodeDownPtrs(p);
  int *eptr = getFullNodeEdgeValues(p);
  const int node_level = getNodeLevel(p);

  decrTempNodeCount(node_level);

#ifdef DEVELOPMENT_CODE
  const int node_height = getNodeHeight(p);
  for (int i=0; i<size; i++) {
    assert(isReducedNode(dptr[i]));
    assert(getNodeHeight(dptr[i]) < node_height);
    assert((dptr[i] == 0 && eptr[i] == INF) ||
        (dptr[i] != 0 && eptr[i] != INF));
  }
#endif

  // quick scan: is this node zero?
  // find min for normalizing later
  int nnz = 0;
  int truncsize = 0;

  int min = INF;
  for (int i = 0; i < size; i++) {
    if (0 != dptr[i]) {
      nnz++;
      truncsize = i;
      DCASSERT(eptr[i] != INF);
      if (eptr[i] < min) min = eptr[i];
    }
  }
  truncsize++;

  if (0 == nnz) {
    // duplicate of 0
    deleteTempNode(p);
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got 0\n", p);
#endif
    p = 0;
    ev = INF;
    return;
  }

  // normalize -- there should be atleast one i s.t. eptr[i] == 0
  DCASSERT(min != INF);
  for (int i = 0; i < size; i++) {
    if (0 != dptr[i]) {
      eptr[i] -= min;
      DCASSERT(eptr[i] >= 0);
    } // else eptr[i] == INF
  }

  // after normalizing, residual edge-value (i.e. min) is pushed up
  // nothing needs to be added ev after this step
  ev = min;

  // check for possible reductions
  if (reductionRule == forest::FULLY_REDUCED &&
      nnz == getLevelSize(node_level) && eptr[0] == 0) {
    // if downpointers are the same and ev are same (i.e. 0 after
    // normalizing), eliminate node
    int i = 1;
    for ( ; i < size && dptr[i] == dptr[0] && eptr[i] == 0; i++);
    if (i == size ) {
      // for all i, dptr[i] == dptr[0] and eptr[i] == 0
      int temp = sharedCopy(dptr[0]);
      deleteTempNode(p);
      p = temp;
      return;
    }
  }

  // check unique table
  int q = find(p);
  if (getNull() != q) {
    // duplicate found
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got %d\n", p, q);
#endif
    linkNode(q);
    deleteTempNode(p);
    p = q;
    return;
  }

  // insert into unique table
  insert(p);

#ifdef TRACE_REDUCE
  printf("\tReducing %d: unique, compressing\n", p);
#endif

  if (!areSparseNodesEnabled())
    nnz = size;

  // right now, tie goes to truncated full.
  if (3*nnz < 2*truncsize) {
    // sparse is better; convert
    int newoffset = getHole(node_level, 4+3*nnz, true);
    // can't rely on previous dptr, re-point to p
    int* full_ptr = getNodeAddress(p);
    int* sparse_ptr = getAddress(node_level, newoffset);
    // copy first 2 integers: incount, next
    sparse_ptr[0] = full_ptr[0];
    sparse_ptr[1] = full_ptr[1];
    // size
    sparse_ptr[2] = -nnz;
    // copy index into address[]
    sparse_ptr[3 + 3*nnz] = p;
    // get pointers to the new sparse node
    int* indexptr = sparse_ptr + 3;
    int* downptr = indexptr + nnz;
    int* edgeptr = downptr + nnz;
    dptr = full_ptr + 3;
    eptr = dptr + size;
    // copy downpointers
    for (int i=0; i<size; i++, ++dptr, ++eptr) {
      if (0 != *dptr) {
        *indexptr = i;
        *downptr = *dptr;
        *edgeptr = *eptr;
        ++indexptr;
        ++downptr;
        ++edgeptr;
      }
    }
    // trash old node
#ifdef MEMORY_TRACE
    int saved_offset = getNodeOffset(p);
    setNodeOffset(p, newoffset);
    makeHole(node_level, saved_offset, 4 + 2 * size);
#else
    makeHole(node_level, getNodeOffset(p), 4 + 2 * size);
    setNodeOffset(p, newoffset);
#endif
    // address[p].cache_count does not change
  } else {
    // full is better
    if (truncsize < size) {
      // truncate the trailing 0s
      int newoffset = getHole(node_level, 4+2*truncsize, true);
      // can't rely on previous ptr, re-point to p
      int* full_ptr = getNodeAddress(p);
      int* trunc_ptr = getAddress(node_level, newoffset);
      // copy first 2 integers: incount, next
      trunc_ptr[0] = full_ptr[0];
      trunc_ptr[1] = full_ptr[1];
      // size
      trunc_ptr[2] = truncsize;
      // copy index into address[]
      trunc_ptr[3 + 2 * truncsize] = p;
      // elements
      memcpy(trunc_ptr + 3, full_ptr + 3, truncsize * sizeof(int));
      // edge values
      memcpy(trunc_ptr + 3 + truncsize, full_ptr + 3 + size,
          truncsize * sizeof(int));
      // trash old node
#ifdef MEMORY_TRACE
      int saved_offset = getNodeOffset(p);
      setNodeOffset(p, newoffset);
      makeHole(node_level, saved_offset, 4 + 2 * size);
#else
      makeHole(node_level, getNodeOffset(p), 4 + 2 * size);
      setNodeOffset(p, newoffset);
#endif
      // address[p].cache_count does not change
    }
  }

  // address[p].cache_count does not change
  DCASSERT(getCacheCount(p) == 0);
  // Sanity check that the hash value is unchanged
  DCASSERT(find(p) == p);

#if 0
  validateIncounts();
#endif

  return;
}


forest::error evplusmdd_node_manager::getElement(int a, int index, int* e)
{
  if (a == 0) return forest::INVALID_VARIABLE;
  if (a == -1) return forest::SUCCESS;

  int down = 0;
  int downIndex = -1;
  int ev = 0;

  if (isFullNode(a)) {
    int aSize = getFullNodeSize(a);
    for (int i = aSize - 1; i >= 0; i--)
    {
      getFullNodeEdgeValue(a, i, ev);
      if (index >= ev) {
        down = getFullNodeDownPtr(a, i);
        downIndex = i;
        break;
      }
    }
  }
  else {
    DCASSERT(isSparseNode(a));
    int aNnz = getSparseNodeSize(a);
    for (int i = aNnz - 1; i >= 0; i--)
    {
      getSparseNodeEdgeValue(a, i, ev);
      if (index >= ev) {
        down = getSparseNodeDownPtr(a, i);
        downIndex = getSparseNodeIndex(a, i);
        break;
      }
    }
  }

  DCASSERT(downIndex >= 0);
  int aLevel = getNodeLevel(a);
  DCASSERT(aLevel >= 0);
  e[aLevel] = downIndex;
  return getElement(down, index - ev, e);
}


forest::error evplusmdd_node_manager::getElement(const dd_edge& a,
    int index, int* e)
{
  assert(e != 0);
  e[0] = 0;
  if (index < 0) return forest::INVALID_VARIABLE;
  return getElement(a.getNode(), index, e);
}


forest::error
evplusmdd_node_manager::
createEdge(const int* const* vlist, const int* terms, int N, dd_edge &e)
{
  return createEdgeInternal(vlist, terms, N, e);
}


forest::error
evplusmdd_node_manager::
createEdge(int val, dd_edge &e)
{
  return createEdgeInternal(val, e);
}


forest::error
evplusmdd_node_manager::
evaluate(const dd_edge &f, const int* vlist, int &term) const
{
  return evaluateInternal(f, vlist, term);
}


int
evplusmdd_node_manager::
createTempNode(int lh, std::vector<int>& downPointers,
    std::vector<int>& edgeValues)
{
  int tempNode =
    evmdd_node_manager::createTempNode(lh, downPointers.size(), false);
  int* dptrs = getFullNodeDownPtrs(tempNode);
  int* evs = getFullNodeEdgeValues(tempNode);
  std::vector<int>::iterator dpiter = downPointers.begin();
  std::vector<int>::iterator eviter = edgeValues.begin();
  while (dpiter != downPointers.end())
  {
    *dptrs++ = *dpiter++;
    *evs++ = *eviter++;
  }
  return tempNode;
}


// ********************************* EV*MDDs ********************************** 

evtimesmdd_node_manager::evtimesmdd_node_manager(domain *d)
: evmdd_node_manager(d, forest::EVTIMES)
{ }

evtimesmdd_node_manager::~evtimesmdd_node_manager()
{ }

void evtimesmdd_node_manager::initEdgeValues(int p) {
  DCASSERT(!isReducedNode(p));
  DCASSERT(isFullNode(p));

  DCASSERT(sizeof(int) == sizeof(float));

#if 0
  int *edgeptr = getFullNodeEdgeValues(p);
  int *last = edgeptr + getFullNodeSize(p);
  for ( ; edgeptr != last; ++edgeptr) *edgeptr = toInt(NAN);
#else
  int *edgeptr = getFullNodeEdgeValues(p);
  float *fptr = (float *)edgeptr;
  float *end = fptr + getFullNodeSize(p);
  for ( ; fptr != end; ++fptr) *fptr = NAN;
#endif
}


// Similar to getDownPtrs() but for EV*MDDs
bool evtimesmdd_node_manager::getDownPtrsAndEdgeValues(int p,
    std::vector<int>& dptrs, std::vector<float>& evs) const
{
  if (!isActiveNode(p) || isTerminalNode(p) || !isReducedNode(p))
    return false;

  if (isFullNode(p)) {
    int size = getFullNodeSize(p);
    if (dptrs.size() < unsigned(size)) {
      dptrs.resize(size, 0);
      evs.resize(size, NAN);
    }
    const int* ptrs = getFullNodeDownPtrsReadOnly(p);
    const int* end = ptrs + size;
    const float* evptrs = (const float *)(getFullNodeEdgeValuesReadOnly(p));
    std::vector<int>::iterator iter = dptrs.begin();
    std::vector<float>::iterator eviter = evs.begin();
    while (ptrs != end)
    {
      *iter++ = *ptrs++;
      *eviter++ = *evptrs++;
    }
  }
  else {
    int nnz = getSparseNodeSize(p);
    int size = getLargestIndex(p) + 1;
    if (dptrs.size() < unsigned(size)) {
      dptrs.resize(size, 0);
      evs.resize(size, INF);
    }
    const int* ptrs = getSparseNodeDownPtrs(p);
    const int* index = getSparseNodeIndexes(p);
    const float* evptrs = (const float *)(getSparseNodeEdgeValues(p));
    const int* end = ptrs + nnz;
    while (ptrs != end)
    {
      evs[*index] = *evptrs++;
      dptrs[*index++] = *ptrs++;
    }
  }
  return true;
}


void evtimesmdd_node_manager::normalizeAndReduceNode(int& p, float& ev)
{
  DCASSERT(isActiveNode(p));

  if (isReducedNode(p)) {
    if (p == 0) ev = NAN;
    return;
  }

  DCASSERT(!isTerminalNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(getInCount(p) == 1);

#if 0
  validateIncounts();
#endif

  const int size = getFullNodeSize(p);
  int *dptr = getFullNodeDownPtrs(p);
  int *eptr = getFullNodeEdgeValues(p);
  float *fptr = (float *)eptr;
  const int node_level = getNodeLevel(p);

  decrTempNodeCount(node_level);

#ifdef DEVELOPMENT_CODE
  const int node_height = getNodeHeight(p);
  for (int i=0; i<size; i++) {
    assert(isReducedNode(dptr[i]));
    assert(getNodeHeight(dptr[i]) < node_height);
    assert((dptr[i] == 0 && isNan(fptr[i])) ||
        (dptr[i] != 0 && !isNan(fptr[i]) && fptr[i] >= 0.0));
  }
#endif

  // quick scan: is this node zero?
  // find max for normalizing later
  int nnz = 0;
  int truncsize = 0;

  float max = 0;
  for (int i = 0; i < size; i++) {
    if (0 != dptr[i]) {
      nnz++;
      truncsize = i;
      if (fptr[i] > max) { max = fptr[i]; }
    }
  }
  truncsize++;

  if (0 == nnz) {
    // duplicate of 0
    deleteTempNode(p);
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got 0\n", p);
#endif
    p = 0;
    ev = NAN;
    return;
  }

  // normalize -- all fptr[i] should be between 0 and 1 (or NAN)
  if (max != 0.0) {
    for (int i = 0; i < size; i++) {
      if (0 != dptr[i]) {
        fptr[i] /= max;
      }
    }
  }

  // after normalizing, residual edge-value (i.e. max) is pushed up
  // nothing needs to be added ev after this step
  ev = max;

  // check for possible reductions
  if (reductionRule == forest::FULLY_REDUCED &&
      nnz == getLevelSize(node_level)) {
    // if downpointers are the same and ev are same, eliminate node
    int i = 1;
    int src = dptr[0];
    for ( ; i < size && dptr[i] == src; i++);
    if (i == size) {
      src = eptr[0];
      for (i = 1; i < size && eptr[i] == src; i++);
      if (i == size ) {
        // for all i, dptr[i] == dptr[0] and eptr[i] == eptr[0]
        int temp = sharedCopy(dptr[0]);
        deleteTempNode(p);
        p = temp;
        return;
      }
    }
  }

  // check unique table
  int q = find(p);
  if (getNull() != q) {
    // duplicate found
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got %d\n", p, q);
#endif
    linkNode(q);
    deleteTempNode(p);
    p = q;
    return;
  }

  // insert into unique table
  insert(p);

#ifdef TRACE_REDUCE
  printf("\tReducing %d: unique, compressing\n", p);
#endif

  if (!areSparseNodesEnabled())
    nnz = size;

  // right now, tie goes to truncated full.
  if (3*nnz < 2*truncsize) {
    // sparse is better; convert
    int newoffset = getHole(node_level, 4+3*nnz, true);
    // can't rely on previous dptr, re-point to p
    int* full_ptr = getNodeAddress(p);
    int* sparse_ptr = getAddress(node_level, newoffset);
    // copy first 2 integers: incount, next
    sparse_ptr[0] = full_ptr[0];
    sparse_ptr[1] = full_ptr[1];
    // size
    sparse_ptr[2] = -nnz;
    // copy index into address[]
    sparse_ptr[3 + 3*nnz] = p;
    // get pointers to the new sparse node
    int* indexptr = sparse_ptr + 3;
    int* downptr = indexptr + nnz;
    int* edgeptr = downptr + nnz;
    dptr = full_ptr + 3;
    eptr = dptr + size;
    // copy downpointers
    for (int i=0; i<size; i++, ++dptr, ++eptr) {
      if (0 != *dptr) {
        *indexptr = i;
        *downptr = *dptr;
        *edgeptr = *eptr;
        ++indexptr;
        ++downptr;
        ++edgeptr;
      }
    }
    // trash old node
#ifdef MEMORY_TRACE
    int saved_offset = getNodeOffset(p);
    setNodeOffset(p, newoffset);
    makeHole(node_level, saved_offset, 4 + 2 * size);
#else
    makeHole(node_level, getNodeOffset(p), 4 + 2 * size);
    setNodeOffset(p, newoffset);
#endif
    // address[p].cache_count does not change
  } else {
    // full is better
    if (truncsize < size) {
      // truncate the trailing 0s
      int newoffset = getHole(node_level, 4+2*truncsize, true);
      // can't rely on previous ptr, re-point to p
      int* full_ptr = getNodeAddress(p);
      int* trunc_ptr = getAddress(node_level, newoffset);
      // copy first 2 integers: incount, next
      trunc_ptr[0] = full_ptr[0];
      trunc_ptr[1] = full_ptr[1];
      // size
      trunc_ptr[2] = truncsize;
      // copy index into address[]
      trunc_ptr[3 + 2 * truncsize] = p;
      // elements
      memcpy(trunc_ptr + 3, full_ptr + 3, truncsize * sizeof(int));
      // edge values
      memcpy(trunc_ptr + 3 + truncsize, full_ptr + 3 + size,
          truncsize * sizeof(int));
      // trash old node
#ifdef MEMORY_TRACE
      int saved_offset = getNodeOffset(p);
      setNodeOffset(p, newoffset);
      makeHole(node_level, saved_offset, 4 + 2 * size);
#else
      makeHole(node_level, getNodeOffset(p), 4 + 2 * size);
      setNodeOffset(p, newoffset);
#endif
      // address[p].cache_count does not change
    }
  }

  // address[p].cache_count does not change
  DCASSERT(getCacheCount(p) == 0);
  // Sanity check that the hash value is unchanged
  DCASSERT(find(p) == p);

#if 0
  validateIncounts();
#endif

  return;
}


forest::error
evtimesmdd_node_manager::
createEdge(const int* const* vlist, const float* terms, int N, dd_edge &e)
{
  return createEdgeInternal(vlist, terms, N, e);
}


forest::error
evtimesmdd_node_manager::
createEdge(float val, dd_edge &e)
{
  return createEdgeInternal(val, e);
}


forest::error
evtimesmdd_node_manager::
evaluate(const dd_edge &f, const int* vlist, float &term) const
{
  return evaluateInternal(f, vlist, term);
}


int
evtimesmdd_node_manager::
createTempNode(int lh, std::vector<int>& downPointers,
    std::vector<float>& edgeValues)
{
  int tempNode =
    evmdd_node_manager::createTempNode(lh, downPointers.size(), false);
  int* dptrs = getFullNodeDownPtrs(tempNode);
  int* ievs = getFullNodeEdgeValues(tempNode);
  float* evs = (float *)ievs;
  std::vector<int>::iterator dpiter = downPointers.begin();
  std::vector<float>::iterator eviter = edgeValues.begin();
  while (dpiter != downPointers.end())
  {
    *dptrs++ = *dpiter++;
    *evs++ = *eviter++;
  }
  return tempNode;
}


// *********************************** MDDs ***********************************

mdd_node_manager::mdd_node_manager(domain *d)
: node_manager(d, false, forest::BOOLEAN,
      forest::MULTI_TERMINAL, forest::FULLY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{ }


mdd_node_manager::~mdd_node_manager()
{ }


int mdd_node_manager::reduceNode(int p)
{
  DCASSERT(isActiveNode(p));

  if (isReducedNode(p)) return p; 

  DCASSERT(!isTerminalNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(getInCount(p) == 1);

#if 0
  validateIncounts();
#endif

  int size = getFullNodeSize(p);
  int* ptr = getFullNodeDownPtrs(p);
  int node_level = getNodeLevel(p);

  decrTempNodeCount(node_level);

#ifdef DEVELOPMENT_CODE
  int node_height = getNodeHeight(p);
  for (int i=0; i<size; i++) {
    assert(isReducedNode(ptr[i]));
    assert(getNodeHeight(ptr[i]) < node_height);
  }
#endif

  // quick scan: is this node zero?
  int nnz = 0;
  int truncsize = 0;
  for (int i = 0; i < size; ++i) {
    if (0 != ptr[i]) {
      nnz++;
      truncsize = i;
    }
  }
  truncsize++;

  if (0 == nnz) {
    // duplicate of 0
    deleteTempNode(p);
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got 0\n", p);
#endif
    return 0;
  }

  // check for possible reductions
  if (reductionRule == forest::FULLY_REDUCED) {
    if (nnz == getLevelSize(node_level)) {
      int i = 1;
      for ( ; i < size && ptr[i] == ptr[0]; i++);
      if (i == size ) {
        // for all i, ptr[i] == ptr[0]
        int temp = sharedCopy(ptr[0]);
        deleteTempNode(p);
        return temp;
      }
    }
  } else {
    // Quasi-reduced -- no identity reduction for MDDs/MTMDDs
    assert(reductionRule == forest::QUASI_REDUCED);

    // ensure than all downpointers are pointing to nodes exactly one
    // level below or zero.
    for (int i = 0; i < size; ++i)
    {
      if (ptr[i] == 0) continue;
      if (getNodeLevel(ptr[i]) != (node_level - 1)) {
        int temp = ptr[i];
        ptr[i] = buildQuasiReducedNodeAtLevel(node_level - 1, ptr[i]);
        unlinkNode(temp);
      }
      assert(ptr[i] == 0 || (getNodeLevel(ptr[i]) == node_level - 1));
    }
  }

  // check unique table
  int q = find(p);
  if (getNull() != q) {
    // duplicate found
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got %d\n", p, q);
#endif
    deleteTempNode(p);
    return sharedCopy(q);
  }

  // insert into unique table
  insert(p);

#ifdef TRACE_REDUCE
  printf("\tReducing %d: unique, compressing\n", p);
#endif

  if (!areSparseNodesEnabled())
    nnz = size;

  // right now, tie goes to truncated full.
  if (2*nnz < truncsize) {
    // sparse is better; convert
    int newoffset = getHole(node_level, 4+2*nnz, true);
    // can't rely on previous ptr, re-point to p
    int* full_ptr = getNodeAddress(p);
    int* sparse_ptr = getAddress(node_level, newoffset);
    // copy first 2 integers: incount, next
    sparse_ptr[0] = full_ptr[0];
    sparse_ptr[1] = full_ptr[1];
    // size
    sparse_ptr[2] = -nnz;
    // copy index into address[]
    sparse_ptr[3 + 2*nnz] = p;
    // get pointers to the new sparse node
    int* indexptr = sparse_ptr + 3;
    int* downptr = indexptr + nnz;
    ptr = full_ptr + 3;
    // copy downpointers
    for (int i=0; i<size; i++, ++ptr) {
      if (*ptr) {
        *indexptr = i;
        *downptr = *ptr;
        ++indexptr;
        ++downptr;
      }
    }
    // trash old node
#ifdef MEMORY_TRACE
    int saved_offset = getNodeOffset(p);
    setNodeOffset(p, newoffset);
    makeHole(node_level, saved_offset, 4 + size);
#else
    makeHole(node_level, getNodeOffset(p), 4 + size);
    setNodeOffset(p, newoffset);
#endif
    // address[p].cache_count does not change
  } else {
    // full is better
    if (truncsize<size) {
      // truncate the trailing 0s
      int newoffset = getHole(node_level, 4+truncsize, true);
      // can't rely on previous ptr, re-point to p
      int* full_ptr = getNodeAddress(p);
      int* trunc_ptr = getAddress(node_level, newoffset);
      // copy first 2 integers: incount, next
      trunc_ptr[0] = full_ptr[0];
      trunc_ptr[1] = full_ptr[1];
      // size
      trunc_ptr[2] = truncsize;
      // copy index into address[]
      trunc_ptr[3 + truncsize] = p;
      // elements
      memcpy(trunc_ptr + 3, full_ptr + 3, truncsize * sizeof(int));
      // trash old node
#ifdef MEMORY_TRACE
      int saved_offset = getNodeOffset(p);
      setNodeOffset(p, newoffset);
      makeHole(node_level, saved_offset, 4 + size);
#else
      makeHole(node_level, getNodeOffset(p), 4 + size);
      setNodeOffset(p, newoffset);
#endif
      // address[p].cache_count does not change
    }
  }

  // address[p].cache_count does not change
  DCASSERT(getCacheCount(p) == 0);
  // Sanity check that the hash value is unchanged
  DCASSERT(find(p) == p);

#if 0
  validateIncounts();
#endif

  return p;
}

int mdd_node_manager::createNode(int k, int index, int dptr)
{
  DCASSERT(index >= -1 && dptr >= -1);
  DCASSERT(false == getBoolean(0));

#if 1
  if (index > -1 && getLevelSize(k) <= index) {
    //printf("level: %d, curr size: %d, index: %d, ",
    //k, getLevelSize(k), index);
    expertDomain->enlargeVariableBound(k, index + 1);
    //printf("new size: %d\n", getLevelSize(k));
  }
#endif

  if (dptr == 0) return 0;
  if (index == -1) {
    // all downpointers should point to dptr
    if (reductionRule == forest::FULLY_REDUCED) return sharedCopy(dptr);
    int curr = createTempNodeMaxSize(k, false);
    setAllDownPtrsWoUnlink(curr, dptr);
    return reduceNode(curr);
  }

#if 0

  // a single downpointer points to dptr
  int curr = createTempNode(k, index + 1);
  setDownPtrWoUnlink(curr, index, dptr);
  return reduceNode(curr);

#else

  // a single downpointer points to dptr
  if (nodeStorage == FULL_STORAGE ||
      (nodeStorage == FULL_OR_SPARSE_STORAGE && index < 2)) {
    // Build a full node
    int curr = createTempNode(k, index + 1);
    setDownPtrWoUnlink(curr, index, dptr);
    return reduceNode(curr);
  }
  else {
    DCASSERT (nodeStorage == SPARSE_STORAGE ||
        (nodeStorage == FULL_OR_SPARSE_STORAGE && index >= 2));
    // Build a sparse node
    int p = createTempNode(k, 2);
    int* nodeData = getNodeAddress(p);
    // For sparse nodes, size is -ve
    nodeData[2] = -1;
    // indexes followed by downpointers -- here we have one index and one dptr
    nodeData[3] = index;
    nodeData[4] = sharedCopy(dptr);
    // search in unique table
    int q = find(p);
    if (getNull() == q) {
      // no duplicate found; insert into unique table
      insert(p);
      DCASSERT(getCacheCount(p) == 0);
      DCASSERT(find(p) == p);
    }
    else {
      // duplicate found; discard this node and return the duplicate
      // revert to full temp node before discarding
      nodeData[2] = 2;
      nodeData[3] = 0;
      nodeData[4] = 0;
      unlinkNode(dptr);
      deleteTempNode(p);
      p = sharedCopy(q);
    }
    decrTempNodeCount(k);
    return p;
  }

#endif

}


void mdd_node_manager::createEdge(const int* v, int term, dd_edge& e)
{
  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;

  int result = term;
  int curr = 0;
  for (int i=1; i<h_sz; i++) {
    curr = createNode(h2l_map[i], v[h2l_map[i]], result);
    unlinkNode(result);
    result = curr;
  }
  e.set(result, 0, getNodeLevel(result));
}


forest::error mdd_node_manager::createEdge(const int* const* vlist, int N,
    dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || N <= 0) return forest::INVALID_VARIABLE;

#ifndef SORT_BUILD
  int trueNode = getTerminalNode(true);
  createEdge(vlist[0], trueNode, e);
  if (N > 1) {
    dd_edge curr(this);
    for (int i=1; i<N; i++) {
      createEdge(vlist[i], trueNode, curr);
      e += curr;
    }
  }
#else

  if (N < 2) {
    createEdge(vlist[0], getTerminalNode(true), e);
  }
  else {
    // check for special cases
    bool specialCasesFound = false;
    for (int i = 0; i < N && !specialCasesFound; i++)
    {
      int* curr = (int*)vlist[i];
      int* end = curr + expertDomain->getNumVariables() + 1;
      for ( ; curr != end; )
      {
        if (*curr++ < 0) {
          specialCasesFound = true;
          break;
        }
      }
    }

    if (specialCasesFound) {
      // build using "standard" procedure
      int trueNode = getTerminalNode(true);
#if 0
      assert(N > 1);
      createEdge(vlist[0], trueNode, e);
      dd_edge curr(this);
      for (int i=1; i<N; i++) {
        createEdge(vlist[i], trueNode, curr);
        e += curr;
      }
#else
      // populate curr[]
      vector<dd_edge*> curr(N, (dd_edge*)0);
      for (vector<dd_edge*>::iterator iter = curr.begin();
          iter != curr.end(); ++iter, ++vlist)
      {
        *iter = new dd_edge(this);
        createEdge(*vlist, trueNode, *(*iter));
      }

      // Iteratively halve curr[] by combining adjacent locations.
      // When curr is of size 1 curr[0] gives the result

      for (vector<dd_edge*> next; curr.size() > 1; curr = next)
      {
        // build next[] and then update curr[]
        next.resize((curr.size()+1)/2, (dd_edge*)0);

        vector<dd_edge*>::iterator currIter = curr.begin();
        vector<dd_edge*>::iterator nextIter = next.begin();
        for ( ; currIter != curr.end(); currIter += 2, ++nextIter)
        {
          *nextIter = *currIter;
        }

        currIter = curr.begin() + 1;
        nextIter = next.begin();
        for ( ; currIter != curr.end(); currIter += 2, ++nextIter)
        {
          *(*nextIter) += *(*currIter);
          delete *currIter;
        }
      }
      assert(curr.size() == 1);
      e = *(curr[0]);
      delete curr[0];
#endif
    }
    else {
      // build using sort-based procedure
      int* list[N];
      memcpy(list, vlist, N * sizeof(int*));
      sortVector(list, 0, N, getDomain()->getNumVariables() + 1);

      int result = sortBuild(list, getDomain()->getNumVariables(), 0, N);
      e.set(result, 0, getNodeLevel(result));

      memset(list, 0, N * sizeof(int*));
    }
  }
#endif
  return forest::SUCCESS;
}


#ifdef SORT_BUILD
int mdd_node_manager::sortBuild(int** list, int height, int begin, int end)
{
  // [begin, end)

  // terminal condition
  if (height == 0)
  {
    assert(begin + 1 == end);
    return getTerminalNode(true);
  }

  int** currList = list;
  int nextHeight = height - 1;
  int level = expertDomain->getVariableWithHeight(height);

  vector<int> nodes;
  int currBegin = begin;
  int currEnd = currBegin;
  for ( ; currEnd < end; currBegin = currEnd) 
  {
    int currIndex = currList[currBegin][level];
    assert(currIndex >= 0);
    for (currEnd = currBegin + 1;
        currEnd < end && currIndex == currList[currEnd][level];
        ++currEnd);
    // found new range
    // to be "unioned" and assigned to result[currIndex]
#ifdef DEBUG_SORT_BUILD
    printf("level: %d, currIndex: %d, currBegin: %d, currEnd: %d\n",
        level, currIndex, currBegin, currEnd);
    fflush(stdout);
#endif
    int node = sortBuild(list, nextHeight, currBegin, currEnd);
#ifdef DEBUG_SORT_BUILD
    printf("level: %d, currIndex: %d, currBegin: %d, currEnd: %d\n",
        level, currIndex, currBegin, currEnd);
    fflush(stdout);
#endif
    nodes.push_back(currIndex);
    nodes.push_back(node);
  }

  assert(nodes.size() >= 2);
  if (nodes.size() == 2) {
    // single entry: store directly as sparse
    // TODO: this should be able to accept many more cases
    int result = createNode(level, nodes[0], nodes[1]);
    unlinkNode(nodes[1]);
    return result;
  }
  // full node
  int size = nodes[nodes.size() - 2] + 1;
  int result = createTempNode(level, size);
  for (vector<int>::iterator iter = nodes.begin();
      iter != nodes.end(); iter += 2)
  {
    if (getDownPtr(result, *iter) != 0) {
      for (unsigned i = 0u; i < nodes.size(); i+=2)
      {
        printf("%d:%d ", nodes[i], nodes[i+1]);
      }
      printf("\n");
      for (int i = begin; i < end; i++)
      {
        printf("%d:%d ", i, currList[i][level]);
      }
      printf("\n");
      exit(1);
    }
    setDownPtrWoUnlink(result, *iter, *(iter+1));
    unlinkNode(*(iter+1));
  }
  return reduceNode(result);
}
#endif


forest::error mdd_node_manager::createEdge(bool term, dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (term == false) return forest::INVALID_VARIABLE;

  if (reductionRule == forest::FULLY_REDUCED) {
    e.set(getTerminalNode(true), 0, domain::TERMINALS);
    return forest::SUCCESS;
  }

  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;
  int result = getTerminalNode(true);
  int curr = getTerminalNode(false);
  for (int i=1; i<h_sz; i++) {
    curr = createTempNodeMaxSize(h2l_map[i], false);
    setAllDownPtrsWoUnlink(curr, result);
    unlinkNode(result);
    result = reduceNode(curr);
  }
  e.set(result, 0, getNodeLevel(result));
  return forest::SUCCESS;
}


forest::error mdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  int node = f.getNode();
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  while (!isTerminalNode(node)) {
    node = getDownPtr(node, vlist[h2l_map[getNodeHeight(node)]]);
  }
  term = getBoolean(node);
  return forest::SUCCESS;
}


forest::error
mdd_node_manager::
findFirstElement(const dd_edge& f, int* vlist) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  // vlist is based on level handles.
  int node = f.getNode();
  if (node == 0) return forest::INVALID_ASSIGNMENT;

  int currLevel = expertDomain->getTopVariable();
  DCASSERT(currLevel != domain::TERMINALS);
  while (currLevel != domain::TERMINALS)
  {
    DCASSERT(node != 0);
    if (currLevel != getNodeLevel(node)) {
      // currLevel is "higher" than node, and has been skipped.
      // Since this is a mdd, reduced nodes enable all paths at the
      // skipped level.
      vlist[currLevel] = 0;   // picking the first index
    } else {
      // find a valid path at this level
      if (isFullNode(node)) {
        int size = getFullNodeSize(node);
        for (int i = 0; i < size; i++)
        {
          int n = getFullNodeDownPtr(node, i);
          if (n != 0) {
            node = n;
            vlist[currLevel] = i;
            break;
          }
        }
      } else {
        vlist[currLevel] = getSparseNodeIndex(node, 0);
        node = getSparseNodeDownPtr(node, 0);
      }
    }
    currLevel = expertDomain->getVariableBelow(currLevel);
  }

  return forest::SUCCESS;
}



void mdd_node_manager::normalizeAndReduceNode(int& p, int& ev)
{
  assert(false);
}


void mdd_node_manager::normalizeAndReduceNode(int& p, float& ev)
{
  assert(false);
}


forest::error mdd_node_manager::createEdge(const int* const* vlist,
    const int* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::createEdge(const int* const* vlist,
    const float* terms, int n, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::createEdge(int val, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::createEdge(float val, dd_edge &e)
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    int &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    float &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  return forest::INVALID_OPERATION;
}


forest::error mdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  return forest::INVALID_OPERATION;
}


// ********************************** Utils **********************************

class vector_sorter {
  public:
    vector_sorter(int* indexes, int terms, int *loc, int *size)
    : pri(indexes), sec(terms), ptr(loc), sz(size) {}
    int getPtr() const { return *ptr; }
    int getSize() const { return *sz; }
    int* pri;
    int sec;
    int* ptr;
    int* sz;
};

bool cmpVector(const vector_sorter *a, const vector_sorter *b)
{
  return (a->getPtr() == INT_MAX)? 
    // compare sec
    a->sec < b->sec:
    // compare pri[ptr]
    a->pri[a->getPtr()] < b->pri[a->getPtr()];
}

void print(FILE* out, const vector_sorter* a)
{
  assert(a != 0 && a->pri != 0);
  fprintf(out, "[%d", a->pri[0]);
  for (int i = 1; i < a->getSize(); ++i)
  {
    fprintf(out, " %d", a->pri[i]);
  }
  fprintf(out, ": %d]", a->sec);
}

#if 0
void sortArray(int** indexes, int* terms, int N, int nVars)
{
  // build objects
  int ptr = 0;
  int sz = nVars;

  // ptr and sz are passed to all objs via pointers.
  // any changes to ptr and sz are visible to the objs.

  vector<vector_sorter*> objs(N);
  if (terms == 0) {
    for (vector<vector_sorter*>::iterator iter = objs.begin();
        iter != objs.end(); ++iter, ++indexes)
    {
      *iter = new vector_sorter(*indexes, 1, &ptr, &sz);
    }
  }
  else {
    for (vector<vector_sorter*>::iterator iter = objs.begin();
        iter != objs.end(); ++iter, ++indexes, ++terms)
    {
      *iter = new vector_sorter(*indexes, *terms, &ptr, &sz);
    }
  }

  // if no terms then do not sort them
  ptr = (terms == 0? 0: INT_MAX);
  std::stable_sort(objs.begin(), objs.end(), cmpVector);

  for (vector<vector_sorter*>::iterator iter = objs.begin();
      iter != objs.end(); ++iter)
  {
    fprintf(stdout, "[%d]: ", iter - objs.begin());
    print(stdout, *iter);
    fprintf(stdout, "\n");
  }

  for (vector<vector_sorter*>::iterator iter = objs.begin();
      iter != objs.end(); ++iter)
  {
    delete *iter;
  }
}
#endif

void sortVector(int** indexes, int* terms, int N, int nVars)
{
  // build objects
  int ptr = 0;
  int sz = nVars;

  // ptr and sz are passed to all objs via pointers.
  // any changes to ptr and sz are visible to the objs.

  vector<vector_sorter*> objs(N);
  if (terms == 0) {
    int** levelA = indexes;
    for (vector<vector_sorter*>::iterator iter = objs.begin();
        iter != objs.end(); ++iter, ++levelA)
    {
      *iter = new vector_sorter(*levelA, 1, &ptr, &sz);
    }
  }
  else {
    int** levelA = indexes;
    int* levelT = terms;
    for (vector<vector_sorter*>::iterator iter = objs.begin();
        iter != objs.end(); ++iter, ++levelA, ++levelT)
    {
      *iter = new vector_sorter(*levelA, *levelT, &ptr, &sz);
    }
  }

#ifdef DEBUG_SORT_MATRIX
  for (vector<vector_sorter*>::iterator iter = objs.begin();
      iter != objs.end(); ++iter)
  {
    fprintf(stdout, "[%d]: ", iter - objs.begin());
    print(stdout, *iter);
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "\n");
#endif

  // if no terms then do not sort them
  if (terms != 0) {
    ptr = INT_MAX;
    std::stable_sort(objs.begin(), objs.end(), cmpVector);
  }

  // index 0 is ignored
  for (ptr = 1; ptr != nVars; ptr++)
  {
    std::stable_sort(objs.begin(), objs.end(), cmpVector);
  }

#ifdef DEBUG_SORT_MATRIX
  fprintf(stdout, "\n");
  for (vector<vector_sorter*>::iterator iter = objs.begin();
      iter != objs.end(); ++iter)
  {
    fprintf(stdout, "[%d]: ", iter - objs.begin());
    print(stdout, *iter);
    fprintf(stdout, "\n");
  }
#endif

  if (terms == 0) {
    for (vector<vector_sorter*>::iterator iter = objs.begin();
        iter != objs.end(); ++iter, ++indexes)
    {
      *indexes = (*iter)->pri;
    }
  } else {
    for (vector<vector_sorter*>::iterator iter = objs.begin();
        iter != objs.end(); ++iter, ++indexes, ++terms)
    {
      *indexes = (*iter)->pri;
      *terms = (*iter)->sec;
    }
  }

  for (vector<vector_sorter*>::iterator iter = objs.begin();
      iter != objs.end(); ++iter)
  {
    delete *iter;
  }
  // exit(1);
}



class matrix_sorter {
  public:
    matrix_sorter(int* indexes, int* pindexes, int terms, int *loc, int *size)
    : pri(indexes), ppri(pindexes), sec(terms), ptr(loc), sz(size) {}
    int getPtr() const { return *ptr; }
    int getSize() const { return *sz; }
    int* pri;
    int* ppri;
    int sec;
    int* ptr;
    int* sz;
};

bool cmpMatrix(const matrix_sorter *a, const matrix_sorter *b)
{
  return (a->getPtr() == INT_MAX)? 
    // compare sec
    a->sec < b->sec:
    // compare pri[ptr]
    (a->getPtr() < 0)?
      // prime
      a->ppri[-(a->getPtr())] < b->ppri[-(a->getPtr())]:
      // unprime
      a->pri[a->getPtr()] < b->pri[a->getPtr()];
}

void print(FILE* out, const matrix_sorter* a)
{
  assert(a != 0 && a->pri != 0 && a->ppri != 0);
  fprintf(out, "[%d", a->pri[0]);
  for (int i = 1; i < a->getSize(); ++i)
  {
    fprintf(out, " %d", a->pri[i]);
  }
  fprintf(out, " -> %d", a->ppri[0]);
  for (int i = 1; i < a->getSize(); ++i)
  {
    fprintf(out, " %d", a->ppri[i]);
  }
  fprintf(out, ": %d]", a->sec);
}

void sortMatrix(int** indexes, int** pindexes, int* terms, int N, int nVars)
{
  // build objects
  int ptr = 0;
  int sz = nVars;

  // ptr and sz are passed to all objs via pointers.
  // any changes to ptr and sz are visible to the objs.

  vector<matrix_sorter*> objs(N);
  if (terms == 0) {
    int** levelA = indexes;
    int** levelB = pindexes;
    for (vector<matrix_sorter*>::iterator iter = objs.begin();
        iter != objs.end(); ++iter, ++levelA, ++levelB)
    {
      *iter = new matrix_sorter(*levelA, *levelB, 1, &ptr, &sz);
    }
  }
  else {
    int** levelA = indexes;
    int** levelB = pindexes;
    int* levelT = terms;
    for (vector<matrix_sorter*>::iterator iter = objs.begin();
        iter != objs.end(); ++iter, ++levelA, ++levelB, ++levelT)
    {
      *iter = new matrix_sorter(*levelA, *levelB, *levelT, &ptr, &sz);
    }
  }

#ifdef DEBUG_SORT_MATRIX
  for (vector<matrix_sorter*>::iterator iter = objs.begin();
      iter != objs.end(); ++iter)
  {
    fprintf(stdout, "[%d]: ", iter - objs.begin());
    print(stdout, *iter);
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "\n");
#endif

  // if no terms then do not sort them
  if (terms != 0) {
    ptr = INT_MAX;
    std::stable_sort(objs.begin(), objs.end(), cmpMatrix);
  }

  // index 0 is ignored
  ptr = -1;
  do {
    std::stable_sort(objs.begin(), objs.end(), cmpMatrix);
    ptr = (ptr < 0)? -ptr: -ptr-1;
  } while (-ptr != nVars);

#ifdef DEBUG_SORT_MATRIX
  fprintf(stdout, "\n");
  for (vector<matrix_sorter*>::iterator iter = objs.begin();
      iter != objs.end(); ++iter)
  {
    fprintf(stdout, "[%d]: ", iter - objs.begin());
    print(stdout, *iter);
    fprintf(stdout, "\n");
  }
#endif

  if (terms == 0) {
    for (vector<matrix_sorter*>::iterator iter = objs.begin();
        iter != objs.end(); ++iter, ++indexes, ++pindexes)
    {
      *indexes = (*iter)->pri;
      *pindexes = (*iter)->ppri;
    }
  } else {
    for (vector<matrix_sorter*>::iterator iter = objs.begin();
        iter != objs.end(); ++iter, ++indexes, ++pindexes, ++terms)
    {
      *indexes = (*iter)->pri;
      *pindexes = (*iter)->ppri;
      *terms = (*iter)->sec;
    }
  }

  for (vector<matrix_sorter*>::iterator iter = objs.begin();
      iter != objs.end(); ++iter)
  {
    delete *iter;
  }
  // exit(1);
}



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


#include "mtmdd.h"

// TODO: Add option to create temporary nodes of max size when
//       accumulating minterms.

// ********************************** MTMDDs **********************************

mtmdd_node_manager::mtmdd_node_manager(int dsl, domain *d, forest::range_type t)
: node_manager(dsl, d, false, t,
      forest::MULTI_TERMINAL, forest::FULLY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION,
      mtmddDataHeaderSize)
{
  list = 0;
  termList = 0;
  listSize = 0;
  count = 0;
  slot = 0;
  countSize = 0;
}


mtmdd_node_manager::mtmdd_node_manager(int dsl, domain *d,
    bool relation, forest::range_type t,
    forest::edge_labeling e, forest::reduction_rule r,
    forest::node_storage s, forest::node_deletion_policy dp)
: node_manager(dsl, d, relation, t, e, r, s, dp, mtmddDataHeaderSize)
{
  list = 0;
  termList = 0;
  listSize = 0;
  count = 0;
  slot = 0;
  countSize = 0;
}



mtmdd_node_manager::~mtmdd_node_manager()
{
  if (list) free(list);
  if (termList) free(termList);
  if (count) free(count);
  if (slot) free(slot);
}


void mtmdd_node_manager::expandCountAndSlotArrays(int size)
{
  if (size <= countSize) return;

  int newCountSize = countSize == 0? 8: countSize;
  while (newCountSize < size) { newCountSize *= 2; }

  count = (int*) realloc(count, newCountSize * sizeof(int));
  slot = (int*) realloc(slot, newCountSize * sizeof(int));
  memset(count + countSize, 0, (newCountSize - countSize) * sizeof(int));
  memset(slot + countSize, 0, (newCountSize - countSize) * sizeof(int));
  countSize = newCountSize;
}


void mtmdd_node_manager::resizeNode(int p, int size)
{
  // This operation can only be performed on Temporary nodes.
  if (!isActiveNode(p) || isTerminalNode(p) || isReducedNode(p)) {
    throw error(error::INVALID_OPERATION);
  }

  // If node is already large enough, do nothing, and return SUCCESS.
  int nodeSize = getFullNodeSize(p);
  if (size <= nodeSize) return;

  MEDDLY_DCASSERT(size > nodeSize);

  // Expand node:
  // (a) Create array of desired size;
  // (b) Copy over data from the old array to the new array;
  // (c) Discard the old array
  // (d) Update pointers

  // Create array of desired size
  int oldDataArraySize = getDataHeaderSize() + nodeSize;
  int newDataArraySize = oldDataArraySize - nodeSize + size;
  int nodeLevel = getNodeLevel(p);
  int newOffset = getHole(nodeLevel, newDataArraySize, true);

  MEDDLY_DCASSERT(newDataArraySize > oldDataArraySize);

  // Pointers to old and new data arrays
  int* prev = getNodeAddress(p);
  int* curr = level[mapLevel(nodeLevel)].data + newOffset;

  // Copy old array to new array
  // Don't copy last location from the old array to the new array
  // -- last location stores the pointer to node p.
  memcpy(curr, prev, (oldDataArraySize - 1) * sizeof(int));

  // Clear out the trailing portion of the new array
  // -- except the last location.
  memset(curr + oldDataArraySize - 1, 0,
      (newDataArraySize - oldDataArraySize) * sizeof(int));

  // Discard the old array
  makeHole(nodeLevel, address[p].offset, oldDataArraySize);

  // Change the node size
  curr[2] = size;

  // Set the back-pointer in the new array
  curr[newDataArraySize - 1] = p;

  // Update the offset field to point to the new array
  address[p].offset = newOffset;

  MEDDLY_DCASSERT(size == getFullNodeSize(p));
  MEDDLY_DCASSERT(p == 
      getNodeAddress(p)[getDataHeaderSize() + getFullNodeSize(p) - 1]);
}


int mtmdd_node_manager::reduceNode(int p)
{
  MEDDLY_DCASSERT(isActiveNode(p));

  if (isReducedNode(p)) return p; 

  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(isFullNode(p));

  int size = getFullNodeSize(p);
  int* ptr = getFullNodeDownPtrs(p);
  int node_level = getNodeLevel(p);

#ifdef DEVELOPMENT_CODE
  validateDownPointers(p);
#endif

  // quick scan: is this node zero?
  int nnz = 0;
  int truncsize = 0;
  {
    int* curr = ptr;
    int* last = curr + size;
    while (curr != last) {
      if (!isReducedNode(*curr)) {
        MEDDLY_DCASSERT(getInCount(*curr) == 1);
        *curr = reduceNode(*curr);
      }
      MEDDLY_DCASSERT(isReducedNode(*curr));
      if (0 != *curr++) {
        ++nnz;
        truncsize = curr - ptr;
      }
    }
  }

  if (0 == nnz) {
    // duplicate of 0
    unlinkNode(p);
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got 0\n", p);
#endif
    return 0;
  }

  // check for possible reductions
  if (reductionRule == forest::QUASI_REDUCED) {
    // ensure than all downpointers are pointing to nodes exactly one
    // level below or zero.
    int nextLevel = node_level - 1;
    for (int i = 0; i < size; ++i)
    {
      if (ptr[i] == 0) continue;
      if (getNodeLevel(ptr[i]) != nextLevel) {
        int temp = ptr[i];
        ptr[i] = buildQuasiReducedNodeAtLevel(nextLevel, ptr[i]);
        unlinkNode(temp);
      }
      MEDDLY_DCASSERT(ptr[i] == 0 || (getNodeLevel(ptr[i]) == nextLevel));
    }
  } else {
    // check for possible reductions
    int temp = 0;
    if (checkForReductions(p, nnz, temp)) {
      linkNode(temp);
      unlinkNode(p);
      return temp;
    }
  }

  // check unique table
#if 0
  int q = find(p);
  if (getNull() != q) {
    // duplicate found
#ifdef TRACE_REDUCE
    printf("\tReducing %d, got %d\n", p, q);
#endif
    unlinkNode(p);
    return sharedCopy(q);
  }

  // insert into unique table
  insert(p);
#else
  int q = replace(p);
  if (q != p) {
    // duplicate found
    unlinkNode(p);
    return sharedCopy(q);
  }
#endif

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
  MEDDLY_DCASSERT(getCacheCount(p) == 0);
  // Sanity check that the hash value is unchanged
  MEDDLY_DCASSERT(find(p) == p);

  // Temporary node has been transformed to a reduced node; decrement
  // temporary node count.
  decrTempNodeCount(node_level);

  return p;
}


int mtmdd_node_manager::createNode(int lh,
    std::vector<int>& index, std::vector<int>& dptr)
{
  // last index in index[] should be the largest.
#ifdef DEVELOPMENT_CODE
  int max = 0;
  for (vector<int>::iterator iter = index.begin();
      iter != index.end(); ++iter)
  {
    if (max < *iter) max = *iter;
  }
  assert(max == index[index.size()-1]);
#endif

  int largestIndex = index[index.size()-1];
  int result = createTempNode(lh, largestIndex+1, true);
  int* ptr = getFullNodeDownPtrs(result);

  for (vector<int>::iterator iIter = index.begin(), dIter = dptr.begin();
      iIter != index.end(); )
  {
    // no need to for any linking because the links are "transferred"
    // from the vector
    ptr[*iIter++] = *dIter++;
  }

  return reduceNode(result);
}


int mtmdd_node_manager::createNode(int k, int index, int dptr)
{
  MEDDLY_DCASSERT(index >= -1);

  if (index > -1 && getLevelSize(k) <= index) {
    expertDomain->enlargeVariableBound(k, false, index + 1);
  }

  if (dptr == 0) return 0;
  if (index == -1) {
    // all downpointers should point to dptr
    if (reductionRule == forest::FULLY_REDUCED) return sharedCopy(dptr);
    int curr = createTempNodeMaxSize(k, false);
    setAllDownPtrsWoUnlink(curr, dptr);
    return reduceNode(curr);
  }

  // a single downpointer points to dptr
  if (nodeStorage == FULL_STORAGE ||
      (nodeStorage == FULL_OR_SPARSE_STORAGE && index < 2)) {
    // Build a full node
    int curr = createTempNode(k, index + 1);
    setDownPtrWoUnlink(curr, index, dptr);
    return reduceNode(curr);
  }
  else {
    MEDDLY_DCASSERT (nodeStorage == SPARSE_STORAGE ||
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
      MEDDLY_DCASSERT(getCacheCount(p) == 0);
      MEDDLY_DCASSERT(find(p) == p);
    }
    else {
      // duplicate found; discard this node and return the duplicate
      // revert to full temp node before discarding
      nodeData[2] = 2;
      nodeData[3] = 0;
      nodeData[4] = 0;
      unlinkNode(dptr);
      unlinkNode(p);
      p = sharedCopy(q);
    }
    return p;
  }
}


void mtmdd_node_manager::createEdge(const int* v, int term, dd_edge& e)
{
  // construct the edge bottom-up
  MEDDLY_DCASSERT(isTerminalNode(term));
  int result = term;
  int curr = 0;
  for (int i=1; i<=expertDomain->getNumVariables(); i++) {
    curr = createNode(i, v[i], result);
    unlinkNode(result);
    result = curr;
  }
  e.set(result, 0, getNodeLevel(result));
  // e.show(stderr, 2);
}


void mtmdd_node_manager::createEdge(const int* const* vlist,
    const int* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || terms == 0 || N <= 0) 
    throw error(error::INVALID_VARIABLE);

  createEdgeInternal(vlist, terms, N, e);
}


void mtmdd_node_manager::createEdgeHelper(int terminalNode, dd_edge& e)
{
  MEDDLY_DCASSERT(isTerminalNode(terminalNode));

  if (reductionRule == forest::FULLY_REDUCED || terminalNode == 0) {
    e.set(terminalNode, 0, 0);
    return;
  }

  // construct the edge bottom-up
  int result = terminalNode;
  int curr = 0;
  for (int i=1; i<=expertDomain->getNumVariables(); i++) {
    curr = createTempNodeMaxSize(i, false);
    setAllDownPtrsWoUnlink(curr, result);
    unlinkNode(result);
    result = reduceNode(curr);
  }
  e.set(result, 0, getNodeLevel(result));
}


void mtmdd_node_manager::createEdge(int term, dd_edge& e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  createEdgeHelper(getTerminalNode(term), e);
}


int mtmdd_node_manager::getTerminalNodeForEdge(int n, const int* vlist) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  while (!isTerminalNode(n)) {
    n = getDownPtr(n, vlist[getNodeHeight(n)]);
  }
  return n;
}


void mtmdd_node_manager::evaluate(const dd_edge &f,
    const int* vlist, int &term) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  term = getInteger(getTerminalNodeForEdge(f.getNode(), vlist));
}


void mtmdd_node_manager::normalizeAndReduceNode(int& p, int& ev)
{
  assert(false);
}


void mtmdd_node_manager::normalizeAndReduceNode(int& p, float& ev)
{
  assert(false);
}


void mtmdd_node_manager::createEdge(const int* const* vlist, int N,
    dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mtmdd_node_manager::createEdge(const int* const* vlist,
    const float* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || terms == 0 || N <= 0) 
    throw error(error::INVALID_VARIABLE);

  createEdgeInternal(vlist, terms, N, e);
}


void mtmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mtmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mtmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mtmdd_node_manager::createEdge(bool val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mtmdd_node_manager::createEdge(float term, dd_edge& e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  createEdgeHelper(getTerminalNode(term), e);
}


void mtmdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  throw error(error::INVALID_OPERATION);
}


void mtmdd_node_manager::evaluate(const dd_edge &f,
    const int* vlist, float &term) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  term = getReal(getTerminalNodeForEdge(f.getNode(), vlist));
}


void mtmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  throw error(error::INVALID_OPERATION);
}


void mtmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  throw error(error::INVALID_OPERATION);
}


void mtmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  throw error(error::INVALID_OPERATION);
}


void
mtmdd_node_manager::findFirstElement(const dd_edge& f, int* vlist) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  // vlist is based on level handles.
  int node = f.getNode();
  if (node == 0) 
    throw error(error::INVALID_ASSIGNMENT);

  for (int currLevel = expertDomain->getNumVariables(); currLevel; currLevel--)
  {
    MEDDLY_DCASSERT(node != 0);
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
  } // for currLevel
}



// *********************************** MDDs ***********************************

mdd_node_manager::mdd_node_manager(int dsl, domain *d)
: mtmdd_node_manager(dsl, d, false, forest::BOOLEAN,
      forest::MULTI_TERMINAL, forest::FULLY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{ }


mdd_node_manager::~mdd_node_manager()
{ }


void mdd_node_manager::createEdge(const int* const* vlist,
    int N, dd_edge &e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || N <= 0) 
    throw error(error::INVALID_VARIABLE);
  mtmdd_node_manager::createEdgeInternal(vlist, (bool*)0, N, e);
}


void mdd_node_manager::createEdge(bool term, dd_edge& e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  createEdgeHelper(getTerminalNode(term), e);
}


void mdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  term = getBoolean(getTerminalNodeForEdge(f.getNode(), vlist));
}


void mdd_node_manager::createEdge(const int* const* vlist,
    const int* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mdd_node_manager::createEdge(const int* const* vlist,
    const float* terms, int n, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mdd_node_manager::createEdge(int val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mdd_node_manager::createEdge(float val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    int &term) const
{
  throw error(error::INVALID_OPERATION);
}


void mdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    float &term) const
{
  throw error(error::INVALID_OPERATION);
}


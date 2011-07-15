
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


#include "mtmxd.h"


#define CREATE_TEMP_NODES_MAX_SIZE_ONLY

#define REDUCE_ACCUMULATE_MXD
  
// ******************************** MTMXDs ******************************* 


mtmxd_node_manager::mtmxd_node_manager(int dsl, domain *d, forest::range_type t)
: node_manager(dsl, d, true, t,
      forest::MULTI_TERMINAL, forest::IDENTITY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION,
      mtmxdDataHeaderSize)
{
  pList = 0;
  unpList = 0;
  tList = 0;
  listSize = 0;
  count = 0;
  slot = 0;
  countSize = 0;
}


mtmxd_node_manager::mtmxd_node_manager(int dsl, domain *d,
    bool relation, forest::range_type t,
    forest::edge_labeling e, forest::reduction_rule r,
    forest::node_storage s, forest::node_deletion_policy dp)
: node_manager(dsl, d, relation, t, e, r, s, dp, mtmxdDataHeaderSize)
{
  unpList = 0;
  pList = 0;
  tList = 0;
  listSize = 0;
  count = 0;
  slot = 0;
  countSize = 0;
}


mtmxd_node_manager::~mtmxd_node_manager()
{
  if (unpList) free(unpList);
  if (pList) free(pList);
  if (tList) free(tList);
  if (count) free(count);
  if (slot) free(slot);
}


void mtmxd_node_manager::expandCountAndSlotArrays(int size)
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


void mtmxd_node_manager::resizeNode(int p, int size)
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


int mtmxd_node_manager::reduceNode(int p)
{
  MEDDLY_DCASSERT(isActiveNode(p));

#ifdef DEVELOPMENT_CODE
  validateDownPointers(p);
#endif

  if (isReducedNode(p)) return p; 

  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(isFullNode(p));

  int size = getFullNodeSize(p);
  int* ptr = getFullNodeDownPtrs(p);
  int node_level = getNodeLevel(p);

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
  int temp = 0;
  if (checkForReductions(p, nnz, temp)) {
    linkNode(temp);
    unlinkNode(p);
    return temp;
  }

#if 0
  // check unique table
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
#ifndef TRUNCATED_REDUCE_IN_PLACE
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
#endif
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


int mtmxd_node_manager::createNode(int k, int index, int dptr)
{
  MEDDLY_DCASSERT(index >= 0 && index < getLevelSize(k) && isValidNodeIndex(dptr));

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

int mtmxd_node_manager::createNode(int k, int index1, int index2, int dptr)
{
  MEDDLY_DCASSERT((index1 >= 0 && index2 >= 0) ||
      (index1 >= -1 && index2 >= -1) ||
      (index1 >= -2 && index2 >= -2 && index1 == index2));

  int result = 0;

  if (reductionRule == forest::IDENTITY_REDUCED) {

    if (index1 == -2) {
      // "don't change"
      result = sharedCopy(dptr);
    }
    else {
      int p = 0;
      if (index2 == -1) {
        // represents "don't care"
        p = createTempNodeMaxSize(-k, false);
        setAllDownPtrsWoUnlink(p, dptr);
        p = reduceNode(p);
      } else {
        p = createNode(-k, index2, dptr);
      }
      if (index1 == -1) {
        // represents "don't care"
        result = createTempNodeMaxSize(k, false);
        setAllDownPtrsWoUnlink(result, p);
        result = reduceNode(result);
      } else {
        result = createNode(k, index1, p);
      }
      unlinkNode(p);
    }

  }
  else if (index1 == -2) {

    // "don't change"
    MEDDLY_DCASSERT(reductionRule == forest::QUASI_REDUCED ||
        reductionRule == forest::FULLY_REDUCED);
    int sz = getLevelSize(k);
    result = createTempNode(k, sz, false);
    int* unprimedDptrs = getFullNodeDownPtrs(result);
    for (int i = 0; i != sz; ++i)
    {
      unprimedDptrs[i] = createNode(-k, i, dptr);
    }
    result = reduceNode(result);

  }
  else if (reductionRule == forest::QUASI_REDUCED) {

    int p = 0;
    if (index2 == -1) {
      // represents "don't care"
      p = createTempNodeMaxSize(-k, false);
      setAllDownPtrsWoUnlink(p, dptr);
      p = reduceNode(p);
    } else {
      p = createNode(-k, index2, dptr);
    }
    if (index1 == -1) {
      // represents "don't care"
      result = createTempNodeMaxSize(k, false);
      setAllDownPtrsWoUnlink(result, p);
      result = reduceNode(result);
    } else {
      result = createNode(k, index1, p);
    }
    unlinkNode(p);

  }
  else {

    // deal with "don't care" for primed level
    int p = index2 == -1? sharedCopy(dptr): createNode(-k, index2, dptr);
    // deal with "don't care" for unprimed level
    result = index1 == -1? sharedCopy(p): createNode(k, index1, p);
    unlinkNode(p);

  }

  return result;
}


int mtmxd_node_manager::createNode(const int* v, const int* vp, int term,
    int startAtHeight, bool primedLevel)
{
  // construct the edge bottom-up
  for (int k=1; k<startAtHeight; k++) {
    int prev = term;
    term = createNode(k, v[k], vp[k], term);
    unlinkNode(prev);
  } // for k

  // deal with height == startAtHeight
  // handle primed level first
  if (primedLevel) {
    // only primed level to be handled at this height
    int prev = term;
    term = createNode(-startAtHeight, vp[startAtHeight], term);
    unlinkNode(prev);
  } else {
    // both primed and unprimed levels to be handled at this height
    int prev = term;
    term = createNode(startAtHeight, v[startAtHeight], vp[startAtHeight], term);
    unlinkNode(prev);
  }

  return term;
}


void mtmxd_node_manager::createEdge(const int* v, const int* vp, int term,
    int startAtHeight, bool primedLevel, dd_edge& e)
{
  term = createNode(v, vp, term, startAtHeight, primedLevel);
  e.set(term, 0, getNodeLevel(term));
}


void mtmxd_node_manager::createEdge(const int* v, const int* vp, int term,
    dd_edge& e)
{
  createEdge(v, vp, term, expertDomain->getNumVariables(), false, e);
}


void mtmxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0 || terms == 0 || N <= 0)
    throw error(error::INVALID_VARIABLE);

  createEdgeInternal(vlist, vplist, terms, N, e);
}


void mtmxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0 || terms == 0 || N <= 0)
    throw error(error::INVALID_VARIABLE);

  createEdgeInternal(vlist, vplist, terms, N, e);
}


int mtmxd_node_manager::createEdge(int dptr)
{
  MEDDLY_DCASSERT(isTerminalNode(dptr));
  if (dptr == 0) return 0;

  if (reductionRule == forest::FULLY_REDUCED) return sharedCopy(dptr);

  // construct the edge bottom-up
  int curr = dptr;
  int prev = 0;
  for (int i=1; i<=expertDomain->getNumVariables(); i++) {
    prev = curr;
    curr = createNode(i, -1, -1, prev);
    unlinkNode(prev);
  }
  return curr;
}


void mtmxd_node_manager::createEdge(int val, dd_edge &e)
{
  MEDDLY_DCASSERT(getRangeType() == forest::INTEGER);
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);

  int node = createEdge(getTerminalNode(val));
  e.set(node, 0, getNodeLevel(node));
}


void mtmxd_node_manager::createEdge(float val, dd_edge &e)
{
  MEDDLY_DCASSERT(getRangeType() == forest::REAL);
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);

  int node = createEdge(getTerminalNode(val));
  e.set(node, 0, getNodeLevel(node));
}


int mtmxd_node_manager::getTerminalNodeForEdge(int n, const int* vlist,
    const int* vplist) const
{
  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contain a single element.

  if (reductionRule == forest::IDENTITY_REDUCED) {

    while (!isTerminalNode(n)) {
      int nLevel = getNodeLevel(n);
      if (nLevel < 0) {
        // Primed Node
        int next = getDownPtr(n, vplist[-nLevel]);
        MEDDLY_DCASSERT(isTerminalNode(next) || isUnprimedNode(next));
        int currHeight = getNodeHeight(n) - 1;
        int nextHeight = getNodeHeight(next);
        if (nextHeight < currHeight) {
          // skipped levels
          for ( ; nextHeight != currHeight; --currHeight)
          {
            if (vlist[currHeight] != vplist[currHeight]) {
              next = 0;
              break;
            }
          }
        }
        n = next;
      }
      else {
        // Unprimed Node
        MEDDLY_DCASSERT(getDownPtr(n, vlist[nLevel]) == 0
            || -nLevel == getNodeLevel(getDownPtr(n, vlist[nLevel])));
        n = getDownPtr(n, vlist[nLevel]);
      }
    }

  }
  else {

    while (!isTerminalNode(n)) {
      n = isPrimedNode(n)
        ? getDownPtr(n, vplist[-getNodeLevel(n)])
        : getDownPtr(n, vlist[getNodeLevel(n)]);
    }

  }
  return n;
}


void mtmxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::INTEGER);
  if (f.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0) 
    throw error(error::INVALID_VARIABLE);

  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  term = getInteger(getTerminalNodeForEdge(f.getNode(), vlist, vplist));
}


void mtmxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::REAL);
  if (f.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0) 
    throw error(error::INVALID_VARIABLE);

  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  term = getReal(getTerminalNodeForEdge(f.getNode(), vlist, vplist));
}


void mtmxd_node_manager::normalizeAndReduceNode(int& p, int& ev)
{
  assert(false);
}


void mtmxd_node_manager::normalizeAndReduceNode(int& p, float& ev)
{
  assert(false);
}


void mtmxd_node_manager::createEdge(const int* const* vlist, int N,
    dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mtmxd_node_manager::createEdge(const int* const* vlist,
    const int* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mtmxd_node_manager::createEdge(const int* const* vlist,
    const float* terms, int n, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mtmxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mtmxd_node_manager::createEdge(bool val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mtmxd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  throw error(error::INVALID_OPERATION);
}


void mtmxd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    int &term) const
{
  throw error(error::INVALID_OPERATION);
}


void mtmxd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    float &term) const
{
  throw error(error::INVALID_OPERATION);
}


void mtmxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  throw error(error::INVALID_OPERATION);
}


void mtmxd_node_manager::findFirstElement(const dd_edge& f,
    int* vlist, int* vplist) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  // vlist is based on level handles.
  int node = f.getNode();
  if (node == 0) 
    throw error(error::INVALID_ASSIGNMENT);

  if (forest::IDENTITY_REDUCED == reductionRule) {

    for (int level = expertDomain->getNumVariables(); level; level--)
    {
      MEDDLY_DCASSERT(node != 0);
      MEDDLY_DCASSERT(isUnprimedNode(node));
      if (level != getNodeLevel(node)) {
        // level is "higher" than node, and has been skipped.
        // Since this is a mxd, reduced nodes enable "don't change" paths
        // at the skipped level.
        vlist[level] = 0;   // picking the first index
        vplist[level] = 0;
      } else {
        // find a valid path at this unprime level
        if (isFullNode(node)) {
          int size = getFullNodeSize(node);
          for (int i = 0; i < size; i++)
          {
            int n = getFullNodeDownPtr(node, i);
            if (n != 0) {
              node = n;
              vlist[level] = i;
              break;
            }
          }
        } else {
          vlist[level] = getSparseNodeIndex(node, 0);
          node = getSparseNodeDownPtr(node, 0);
        }

        MEDDLY_DCASSERT(!isTerminalNode(node));
        // can't be -1 because that violates MXD properties
        // can't be 0 because node cannot be set to 0 in the above construct.
        MEDDLY_DCASSERT(isPrimedNode(node));
        // find a valid path at this prime level
        if (isFullNode(node)) {
          int size = getFullNodeSize(node);
          for (int i = 0; i < size; i++)
          {
            int n = getFullNodeDownPtr(node, i);
            if (n != 0) {
              node = n;
              vplist[level] = i;
              break;
            }
          }
        } else {
          vplist[level] = getSparseNodeIndex(node, 0);
          node = getSparseNodeDownPtr(node, 0);
        }
      }
    } // for level

  }
  else {

    // !IDENTITY_REDUCED
    for (int currLevel = expertDomain->getNumVariables(); currLevel; )
    {
      MEDDLY_DCASSERT(node != 0);
      if (currLevel != getNodeLevel(node)) {
        // currLevel been skipped. !IDENTITY_REDUCED ==> reduced nodes
        // enable all paths at the skipped level. Pick the first index.
        if (currLevel < 0) {
          vplist[-currLevel] = 0;
        } else {
          vlist[currLevel] = 0;
        }
      } else {
        // find a valid path at this level
        if (isFullNode(node)) {
          int size = getFullNodeSize(node);
          for (int i = 0; i < size; i++)
          {
            int n = getFullNodeDownPtr(node, i);
            if (n != 0) {
              node = n;
              if (currLevel < 0) {
                vplist[-currLevel] = i;
              } else {
                vlist[currLevel] = i;
              }
              break;
            }
          }
          // Note: since n is not an empty node the for loop is
          // guaranteed to overwrite "node".
        } else {
          if (currLevel < 0) {
            vplist[-currLevel] = getSparseNodeIndex(node, 0);
          } else {
            vlist[currLevel] = getSparseNodeIndex(node, 0);
          }
          node = getSparseNodeDownPtr(node, 0);
        }
      }
      // Set next level
      currLevel = currLevel < 0
                    ? (-currLevel)-1
                    : -currLevel;
    } // currlevel


  }
}


// *********************************** MXDs *********************************** 


mxd_node_manager::mxd_node_manager(int dsl, domain *d)
: mtmxd_node_manager(dsl, d, true, forest::BOOLEAN,
      forest::MULTI_TERMINAL, forest::IDENTITY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{ }


mxd_node_manager::~mxd_node_manager()
{ }


void mxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0 || N <= 0) 
    throw error(error::INVALID_VARIABLE);
  createEdgeInternal(vlist, vplist, (bool*)0, N, e);
}


void mxd_node_manager::createEdge(bool val, dd_edge &e)
{
  MEDDLY_DCASSERT(getRangeType() == forest::BOOLEAN);
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);

  int node = createEdge(getTerminalNode(val));
  e.set(node, 0, getNodeLevel(node));
}


void mxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::BOOLEAN);
  if (f.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0) 
    throw error(error::INVALID_VARIABLE);

  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  term = getBoolean(getTerminalNodeForEdge(f.getNode(), vlist, vplist));
}


void mxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mxd_node_manager::createEdge(int val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mxd_node_manager::createEdge(float val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void mxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  throw error(error::INVALID_OPERATION);
}


void mxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  throw error(error::INVALID_OPERATION);
}


#if 1

int mxd_node_manager::buildQRIdentityNode(int node, int level)
{
  int size = getLevelSize(level);
  int result = createTempNode(level, size, false);
  for (int i = 0; i < size; )
  {
#ifndef CREATE_TEMP_NODES_MAX_SIZE_ONLY
    int temp = createTempNode(-level, i + 1, true);
#else
    int temp = createTempNode(-level, size, true);
#endif
    setDownPtrWoUnlink(temp, i, node);
    setDownPtrWoUnlink(result, i, temp);
    unlinkNode(temp);
  }
  return result;
}


#ifdef REDUCE_ACCUMULATE_MXD

int mxd_node_manager::accumulateExpandA(int a, int b, bool cBM)
{
  // a[i][i] += b
  // return reduceNode(a)

  MEDDLY_DCASSERT(getReductionRule() == forest::IDENTITY_REDUCED);
  MEDDLY_DCASSERT(getMappedNodeHeight(a) > getMappedNodeHeight(b));
  MEDDLY_DCASSERT(!isReducedNode(a));
  MEDDLY_DCASSERT(getInCount(a) == 1);
  MEDDLY_DCASSERT(cBM == false);
  MEDDLY_DCASSERT(b != 0);

  int aSize = getFullNodeSize(a);
  int aLevel = getNodeLevel(a);
  int levelSize = getLevelSize(aLevel);

  if (aSize < levelSize) {
    resizeNode(a, levelSize);
    aSize = getFullNodeSize(a);
  }

  MEDDLY_DCASSERT(aSize == levelSize);

  for (int i = aSize; --i >= 0; ) {
    int dptr = getFullNodeDownPtr(a, i);
    if (isReducedNode(dptr)) {
      MEDDLY_DCASSERT(-1 != dptr);
      int pDptr = dptr == 0? 0: getDownPtr(dptr, i);
      int acc = pDptr == 0? sharedCopy(b): addReducedNodes(pDptr, b);
      if (pDptr != acc) {
        int pNode = dptr == 0
          ? createTempNode(-aLevel, i + 1, true)
          : makeACopy(dptr, i + 1);
        setDownPtr(pNode, i, acc);
        pNode = reduceNode(pNode);
        setDownPtr(a, i, pNode);
        unlinkNode(pNode);
      }
      unlinkNode(acc);
    } else {
      MEDDLY_DCASSERT(getInCount(dptr) == 1);
      int pDptr = getFullNodeDownPtr(dptr, i);
      int acc =
        pDptr == 0
        ? sharedCopy(b)
        : accumulateMxd(pDptr, b,
            (pDptr == -1? false: (getInCount(pDptr) > 1)));
      MEDDLY_DCASSERT(isReducedNode(acc));
      if (pDptr != acc) {
        if (getFullNodeSize(dptr) <= i) { resizeNode(dptr, i + 1); }
        setDownPtr(dptr, i, acc);
        dptr = reduceNode(dptr);
        setDownPtr(a, i, dptr);
        unlinkNode(dptr);
      }
      unlinkNode(acc);
    }
  }

  return reduceNode(a);
}


void mxd_node_manager::accumulateMxdHelper(int& a, int b, bool cBM,
    bool needsToMakeACopy,
    int (mxd_node_manager::*function)(int, int, bool))
{
  // Expand both nodes. a is a full node, b can be either sparse or full.

  MEDDLY_DCASSERT(!isReducedNode(a));
  MEDDLY_DCASSERT(getInCount(a) == 1);
  MEDDLY_DCASSERT(isReducedNode(b));
  MEDDLY_DCASSERT(!needsToMakeACopy);
  MEDDLY_DCASSERT(!cBM);

  // Node b is Truncated-Full
  if (isFullNode(b)) {
    int size = getFullNodeSize(b);
    // Resize a.
    if (getFullNodeSize(a) < size) {
      resizeNode(a, size);
      MEDDLY_DCASSERT(getFullNodeSize(a) == size);
    }
    // Accumulate into a.
    int* aDptrs = getFullNodeDownPtrs(a);
    const int* bDptrs = getFullNodeDownPtrsReadOnly(b);
    for (int i = 0; i < size; ++i) {
      bool unlinkAOfI = isReducedNode(aDptrs[i]);
      int result = (this->*function)(aDptrs[i], bDptrs[i], cBM);
      MEDDLY_DCASSERT(isReducedNode(result));
      if (unlinkAOfI) unlinkNode(aDptrs[i]);
      aDptrs[i] = result;
    }
  }
  // Node b is Sparse
  else {
    MEDDLY_DCASSERT(isSparseNode(b));
    int nDptrs = getSparseNodeSize(b);
    int size = 1 + getSparseNodeIndex(b, nDptrs - 1);
    // Resize a.
    if (getFullNodeSize(a) < size) {
      resizeNode(a, size);
      MEDDLY_DCASSERT(getFullNodeSize(a) == size);
    }
    // Accumulate into a.
    int* aDptrs = getFullNodeDownPtrs(a);
    const int* bDptrs = getSparseNodeDownPtrs(b);
    const int* bIndexes = getSparseNodeIndexes(b);
    for (int i = 0; i < nDptrs; ++i) {
      int index = *bIndexes++;
      int* currA = aDptrs + index;
      bool unlinkCurrA = isReducedNode(*currA);
      int result = (this->*function)(*currA, *bDptrs++, cBM);
      MEDDLY_DCASSERT(isReducedNode(result));
      if (unlinkCurrA) unlinkNode(*currA);
      *currA = result;
    }
  }
}


int mxd_node_manager::accumulateMxdPrime(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(getReductionRule() == forest::IDENTITY_REDUCED);
  MEDDLY_DCASSERT(isReducedNode(b));
  MEDDLY_DCASSERT(a != -1);
  MEDDLY_DCASSERT(b != -1);

  // Terminal nodes
  if (a == 0 || b == 0) {
    int result = a + b;
    return isReducedNode(result)? sharedCopy(result): reduceNode(result);
  }

  MEDDLY_DCASSERT(getNodeLevel(a) == getNodeLevel(b));

  // a is a reduced node
  if (isReducedNode(a)) { return addPrimeReducedNodes(a, b); }

  if (getInCount(a) > 1) cBM = true;
  bool needsToMakeACopy = cBM;

  // Expand both nodes. a is a full node, b can be either sparse or full.
  accumulateMxdHelper(a, b, cBM, needsToMakeACopy, 
      &mxd_node_manager::accumulateMxd);

  return reduceNode(a);
}


int mxd_node_manager::accumulateMxd(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(getReductionRule() == forest::IDENTITY_REDUCED);
  MEDDLY_DCASSERT(isReducedNode(b));

  // Terminal nodes
  if (a == 0 || b == 0) { 
    int result = a + b;
    return isReducedNode(result)? sharedCopy(result): reduceNode(result);
  }

  // a is a reduced node
  if (isReducedNode(a)) { return addReducedNodes(a, b); }

  // a is a temporary node
  int aHeight = getMappedNodeHeight(a);
  int bHeight = getMappedNodeHeight(b);

  MEDDLY_DCASSERT(aHeight >= bHeight);

  if (getInCount(a) > 1) cBM = true;

  if (aHeight > bHeight) {
    // b's levels were skipped.
    // only identity-reduced Mxds.
    // a[i] += b
    return accumulateExpandA(a, b, cBM);
  }

  bool needsToMakeACopy = cBM;

  // Expand both nodes. a is a full node, b can be either sparse or full.
  accumulateMxdHelper(a, b, cBM, needsToMakeACopy,
      &mxd_node_manager::accumulateMxdPrime);

  return reduceNode(a);
}


#else

int mxd_node_manager::accumulateExpandA(int a, int b, bool cBM)
{
  // a[i][i] += b

  MEDDLY_DCASSERT(getReductionRule() == forest::IDENTITY_REDUCED);
  MEDDLY_DCASSERT(getMappedNodeHeight(a) > getMappedNodeHeight(b));
  MEDDLY_DCASSERT(!isReducedNode(a));

  bool needsToMakeACopy = cBM;
  int savedTempNode = a;

  int aSize = getFullNodeSize(a);
  int aLevel = getNodeLevel(a);
  int levelSize = getLevelSize(aLevel);

  if (aSize < levelSize) {
    if (needsToMakeACopy) {
      a = makeACopy(a, levelSize);
      needsToMakeACopy = false;
    } else {
      resizeNode(a, levelSize);
    }
    aSize = getFullNodeSize(a);
  }

  MEDDLY_DCASSERT(aSize == levelSize);

  for (int i = 0; i < aSize; i++) {
    int dptr = getFullNodeDownPtr(a, i);
    MEDDLY_DCASSERT(dptr != -1);
    int pdptr = getDownPtr(dptr, i);
    bool pcBM = dptr == 0? cBM: (cBM || (1 < getInCount(dptr)));
    int result = accumulateMxd(pdptr, b, pcBM);

    if (result != pdptr) {
      // Need to modify a[i]
      int pNode = 0;
#ifndef CREATE_TEMP_NODES_MAX_SIZE_ONLY
      int pNodeSize = i + 1;
#else
      int pNodeSize = getLevelSize(-aLevel);
#endif
      if (dptr == 0) {
        pNode = createTempNode(-aLevel, pNodeSize, true);
      }
      else if (isReducedNode(dptr) || pcBM) {
        pNode = makeACopy(dptr, pNodeSize);
      }
      else {
        if (getFullNodeSize(dptr) < pNodeSize)
          assert(forest::SUCCESS == resizeNode(dptr, pNodeSize));
        pNode = sharedCopy(dptr);
      }

      setDownPtr(pNode, i, result);

      if (pNode != dptr) {
        if (needsToMakeACopy) {
          a = makeACopy(a);
          needsToMakeACopy = false;
        }
        setDownPtr(a, i, pNode);
      }
      unlinkNode(pNode);
    }
    unlinkNode(result);
  }

  return savedTempNode == a? sharedCopy(a): a;
}


void mxd_node_manager::accumulateMxdHelper(int& a, int b, bool cBM,
    bool needsToMakeACopy,
    int (mxd_node_manager::*function)(int, int, bool))
{
  // Expand both nodes. a is a full node, b can be either sparse or full.

  // Node b is Truncated-Full
  if (isFullNode(b)) {
    int size = getFullNodeSize(b);
    // Resize a.
    if (getFullNodeSize(a) < size) {
      if (needsToMakeACopy) {
        a = makeACopy(a, size);
        needsToMakeACopy = false;
      } else {
        resizeNode(a, size);
      }
      MEDDLY_DCASSERT(getFullNodeSize(a) == size);
    }
    // Accumulate into a.
    int* aDptrs = getFullNodeDownPtrs(a);
    const int* bDptrs = getFullNodeDownPtrsReadOnly(b);
    for (int i = 0; i < size; ++i) {
      int result = (this->*function)(aDptrs[i], *bDptrs++, cBM);
      if (result == aDptrs[i]) { unlinkNode(result); continue; }
      if (needsToMakeACopy) {
        a = makeACopy(a);
        needsToMakeACopy = false;
        aDptrs = getFullNodeDownPtrs(a);
      }
      unlinkNode(aDptrs[i]);
      aDptrs[i] = result;
    }
  }
  // Node b is Sparse
  else {
    MEDDLY_DCASSERT(isSparseNode(b));
    int nDptrs = getSparseNodeSize(b);
    int size = 1 + getSparseNodeIndex(b, nDptrs - 1);
    // Resize a.
    if (getFullNodeSize(a) < size) {
      if (needsToMakeACopy) {
        a = makeACopy(a, size);
        needsToMakeACopy = false;
      } else {
        resizeNode(a, size);
      }
      MEDDLY_DCASSERT(getFullNodeSize(a) == size);
    }
    // Accumulate into a.
    int* aDptrs = getFullNodeDownPtrs(a);
    const int* bDptrs = getSparseNodeDownPtrs(b);
    const int* bIndexes = getSparseNodeIndexes(b);
    for (int i = 0; i < nDptrs; ++i) {
      int index = *bIndexes++;
      int* currA = aDptrs + index;
      int result = (this->*function)(*currA, *bDptrs++, cBM);
      if (result == *currA) { unlinkNode(result); continue; }
      if (needsToMakeACopy) {
        a = makeACopy(a);
        needsToMakeACopy = false;
        aDptrs = getFullNodeDownPtrs(a);
        currA = aDptrs + index;
      }
      unlinkNode(*currA);
      *currA = result;
    }
  }
}


// TODO: Make temporary nodes of max size.
int mxd_node_manager::accumulateMxdPrime(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(getReductionRule() == forest::IDENTITY_REDUCED);
  MEDDLY_DCASSERT(isReducedNode(b));
  MEDDLY_DCASSERT(a != -1);
  MEDDLY_DCASSERT(b != -1);

  // Terminal nodes
  if (a == 0 || b == 0) { return sharedCopy(a + b); }

  MEDDLY_DCASSERT(getNodeLevel(a) == getNodeLevel(b));

  // a is a reduced node
  if (isReducedNode(a)) { return addPrimeReducedNodes(a, b); }

  if (getInCount(a) > 1) cBM = true;
  bool needsToMakeACopy = cBM;
  int savedTempNode = a;

  // Expand both nodes. a is a full node, b can be either sparse or full.
  accumulateMxdHelper(a, b, cBM, needsToMakeACopy, 
      &mxd_node_manager::accumulateMxd);

  return savedTempNode == a? sharedCopy(a): a;
}


// TODO: Make temporary nodes of max size.
int mxd_node_manager::accumulateMxd(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(getReductionRule() == forest::IDENTITY_REDUCED);
  MEDDLY_DCASSERT(isReducedNode(b));

  // Terminal nodes
  if (a == 0 || b == 0) { return sharedCopy(a + b); }

  // a is a reduced node
  if (isReducedNode(a)) { return addReducedNodes(a, b); }

  // a is a temporary node
  int aHeight = getMappedNodeHeight(a);
  int bHeight = getMappedNodeHeight(b);

  if (getInCount(a) > 1) cBM = true;

  if (aHeight > bHeight) {
    // b's levels were skipped.
    // only identity-reduced Mxds.
    // a[i] += b
    return accumulateExpandA(a, b, cBM);
  }

  bool needsToMakeACopy = cBM;
  int savedTempNode = a;

  if (aHeight < bHeight) {
    // Build node c at the same level as b.
    // set all c[i] = a;
    a = buildQRIdentityNode(a, getNodeLevel(b));
    needsToMakeACopy = false;
  }

  // Expand both nodes. a is a full node, b can be either sparse or full.
  accumulateMxdHelper(a, b, cBM, needsToMakeACopy,
      &mxd_node_manager::accumulateMxdPrime);

  return savedTempNode == a? sharedCopy(a): a;
}


#endif


int mxd_node_manager::addPrimeReducedNodes(int a, int b)
{
  MEDDLY_DCASSERT(getNodeLevel(a) < 0);
  MEDDLY_DCASSERT(getNodeLevel(a) == getNodeLevel(b));
  MEDDLY_DCASSERT(isReducedNode(a));
  MEDDLY_DCASSERT(isReducedNode(b));

  int level = getNodeLevel(a);
  int size = getLevelSize(level);
  int result = makeACopy(a, size);

  if (isFullNode(b)) {
    int bSize = getFullNodeSize(b);
    for (int i = 0; i < bSize; ++i) {
      if (getFullNodeDownPtr(b, i) == 0) continue;
      int temp = accumulateMxd(getFullNodeDownPtr(result, i),
          getFullNodeDownPtr(b, i), true);
      setDownPtr(result, i, temp);
      unlinkNode(temp);
    }
  }
  else {
    int nDptrs = getSparseNodeSize(b);
    for (int i = 0; i < nDptrs; ++i) {
      int index = getSparseNodeIndex(b, i);
      int temp = accumulateMxd(getFullNodeDownPtr(result, index),
          getSparseNodeDownPtr(b, i), true);
      setDownPtr(result, index, temp);
      unlinkNode(temp);
    }
  }

  return reduceNode(result);
}


void mxd_node_manager::accumulate(int& a, int b)
{
  if (isActiveNode(a) && isActiveNode(b)) {
    // validateDownPointers(a, true);
    // validateDownPointers(b, true);
    bool unlinkA = isReducedNode(a);
    int result = accumulateMxd(a, b, false);
    if (unlinkA) unlinkNode(a);
    a = result;
    // validateDownPointers(a, true);
    return;
  }
  throw error(error::INVALID_OPERATION);
}



#endif


int mxd_node_manager::accumulateSkippedLevel(int tempNode,
    int* element, int* pelement, int level)
{
  MEDDLY_DCASSERT(getReductionRule() == forest::IDENTITY_REDUCED);
  MEDDLY_DCASSERT(level > 0);
  MEDDLY_DCASSERT(getNodeLevel(tempNode) >= 0);
  MEDDLY_DCASSERT(getNodeLevel(tempNode) != level);

  int index = element[level];
  int pindex = pelement[level];

  // Since the forest is Identity-Reduced, it follows:
  int dptr = index == pindex? tempNode: 0;

  int newDptr = accumulate(dptr, true, element, pelement, level-1);

  if (newDptr == dptr) {
    // Element got absorbed into dptr
    return tempNode;
  }

  // Since this is a skipped level node,
  // if newDptr != dptr,
  // a new node must be constructed regardless of cBM.
  int node = 0;
  if (tempNode == 0) {
#ifndef CREATE_TEMP_NODES_MAX_SIZE_ONLY
    int pNode = createTempNode(-level, pindex + 1, true);
#else
    int pNode = createTempNode(-level, getLevelSize(-level), true);
#endif
    setDownPtr(pNode, pindex, newDptr);
#ifndef CREATE_TEMP_NODES_MAX_SIZE_ONLY
    node = createTempNode(level, index + 1, true);
#else
    node = createTempNode(level, getLevelSize(level), true);
#endif
    setDownPtr(node, index, pNode);
    unlinkNode(pNode);
  } else {
    int size = getLevelSize(level);
#ifndef CREATE_TEMP_NODES_MAX_SIZE_ONLY
    int pSize = 0;
#else
    int pSize = getLevelSize(-level);
#endif
    node = createTempNode(level, size, false);
    for (int i = 0; i < size; ++i) {
#ifndef CREATE_TEMP_NODES_MAX_SIZE_ONLY
      pSize = 1 + (i == index? MAX( i , pindex ) : i);
#endif
      int pNode = createTempNode(-level, pSize, true);
      setDownPtrWoUnlink(pNode, i, tempNode);
      if (i == index) setDownPtr(pNode, pindex, newDptr);
      setDownPtr(node, i, pNode);
      unlinkNode(pNode);
    }
  }

  unlinkNode(newDptr);
  return node;
}


// Add an element to a temporary edge
// Start this recursion at the top level in the domain.
// Use expert_domain::getNumVariables() to obtain the topmost level in
// the domain.
// cBM: copy before modifying.
int mxd_node_manager::accumulate(int tempNode, bool cBM,
    int* element, int* pelement, int level)
{
  MEDDLY_DCASSERT(isMxd());
  MEDDLY_DCASSERT(this->getReductionRule() == forest::IDENTITY_REDUCED);

  if (level == 0) {
    if (tempNode == 0) accumulateMintermAddedElement = true;
    return -1;
  }
  MEDDLY_DCASSERT(level > 0);

  int nodeLevel = getNodeLevel(tempNode);
  MEDDLY_DCASSERT(nodeLevel >= 0);

  if (level != nodeLevel) {
    return accumulateSkippedLevel(tempNode, element, pelement, level);
  }

  int index = element[level];
  int pindex = pelement[level];

  // (1) Compute result for the next unprimed level.

  int dptr = getDownPtr(tempNode, index);
  MEDDLY_DCASSERT(!isTerminalNode(dptr) || dptr == 0);
  int pdptr = getDownPtr(dptr, pindex);

  int inCount = getInCount(tempNode);
  int pinCount = isTerminalNode(dptr)? 1: getInCount(dptr);

  // An incount > 1 indicates a need to duplicate the node before
  // modifying.
  if (inCount > 1) cBM = true;
  bool pcBM = cBM || (pinCount > 1);

  int newpDptr = accumulate(pdptr, pcBM, element, pelement, level-1);

  if (newpDptr == pdptr) {
    // Element got absorbed into pdptr
    return tempNode;
  }

  // (2) Create/update node at the primed level.

  // If dptr is 0, create a temporary node.
  // If dptr is a reduced node or if its incount > 1,
  //    create a copy (which is a temporary node).
  // Otherwise, use dptr (should be a temporary node with incount == 1).
  int pNode = 0;
#ifndef CREATE_TEMP_NODES_MAX_SIZE_ONLY
  int pNodeSize = pindex + 1;
#else
  int pNodeSize = getLevelSize(-level);
#endif
  if (dptr == 0) {
    pNode = createTempNode(-level, pNodeSize, true);
  } else if (isReducedNode(dptr) || pcBM) {
    pNode = makeACopy(dptr, pNodeSize);
  } else {
    pNode = dptr;
    MEDDLY_DCASSERT(!isReducedNode(pNode));
    if (pindex >= getFullNodeSize(pNode)) {
      resizeNode(pNode, pNodeSize);
    }
  }

  MEDDLY_DCASSERT(!isReducedNode(pNode));
  MEDDLY_DCASSERT(pindex < getFullNodeSize(pNode));
  setDownPtr(pNode, pindex, newpDptr);
  unlinkNode(newpDptr);
  
  // (3) Create/update node at the unprimed level.

  if (pNode == dptr) {
    // Element got absorbed into dptr
    return tempNode;
  }

  int node = 0;
#ifndef CREATE_TEMP_NODES_MAX_SIZE_ONLY
  int nodeSize = index + 1;
#else
  int nodeSize = getLevelSize(level);
#endif
  if (tempNode == 0) {
    node = createTempNode(level, nodeSize, true);
  } else if (isReducedNode(tempNode) || cBM) {
    node = makeACopy(tempNode, nodeSize);
  } else {
    node = tempNode;
    MEDDLY_DCASSERT(!isReducedNode(node));
    if (index >= getFullNodeSize(node)) {
      resizeNode(node, nodeSize);
    }
  }

  MEDDLY_DCASSERT(!isReducedNode(node));
  MEDDLY_DCASSERT(index < getFullNodeSize(node));
  setDownPtr(node, index, pNode);
  unlinkNode(pNode);

  return node;
}


#if 0

// Add an element to a temporary edge
bool mxd_node_manager::accumulate(int& tempNode,
    int* element, int* pelement)
{
  assert(isActiveNode(tempNode));
  assert(element != 0);
  assert(pelement != 0);

  // Enlarge variable bounds if necessary
  int level = domain::TERMINALS;
  while (-1 != (level = expertDomain->getVariableAbove(level))) {
    int sz = MAX( element[level] , pelement[level] ) + 1;
    if (sz > expertDomain->getVariableBound(level)) {
      expertDomain->enlargeVariableBound(level, false, sz);
    }
  }

  accumulateMintermAddedElement = false;
  int result = accumulate(tempNode, false,
      element, pelement, expertDomain->getNumVariables());
  if (tempNode != result) {
    // tempNode had to be copied into another node by accumulate().
    // This could be either because tempNode was a reduced node,
    // or because tempNode had incount > 1.
    unlinkNode(tempNode);
    tempNode = result;
  }
  // Note: tempNode == result indicates that the element was added
  // to the existing temporary node. Therefore, there is no need to
  // change incounts.
  return accumulateMintermAddedElement;
}


#else
// Add an element to a temporary edge
bool mxd_node_manager::accumulate(int& tempNode,
    int* element, int* pelement)
{
  assert(isActiveNode(tempNode));
  assert(element != 0);
  assert(pelement != 0);
  MEDDLY_DCASSERT(tempNode == 0 || !isReducedNode(tempNode));

  // Enlarge variable bounds if necessary
#if 0
  int level = domain::TERMINALS;
  while (-1 != (level = expertDomain->getVariableAbove(level))) {
    int sz = MAX( element[level] , pelement[level] ) + 1;
    if (sz > expertDomain->getVariableBound(level)) {
      expertDomain->enlargeVariableBound(level, false, sz);
    }
  }
#endif

  // Traverse Mdd till you find a 0.
  int parentNode = 0;
  int currNode = tempNode;
  int currLevel = expertDomain->getNumVariables();
  if (currNode != 0) {
    for ( ; currLevel ; currLevel-- )
    {
      // Skipped levels are not possible with temporary nodes.
      MEDDLY_DCASSERT(!isTerminalNode(currNode));
      MEDDLY_DCASSERT(getNodeLevel(currNode) == currLevel);

      // Unprimed level
      parentNode = currNode;
      int index = element[currLevel];
      if (index >= getFullNodeSize(currNode)) {
        // found a 0
#ifdef CREATE_TEMP_NODES_MAX_SIZE_ONLY
        resizeNode(currNode, getLevelSize(currLevel));
#else
        resizeNode(currNode, index + 1);
#endif
        currNode = 0;
        break;
      }
      currNode = getFullNodeDownPtr(currNode, index);
      if (currNode == 0) break;

      // Primed level
      parentNode = currNode;
      index = pelement[currLevel];
      if (index >= getFullNodeSize(currNode)) {
        // found a 0
#ifdef CREATE_TEMP_NODES_MAX_SIZE_ONLY
        resizeNode(currNode, getLevelSize(currLevel));
#else
        resizeNode(currNode, index + 1);
#endif
        currNode = 0;
        break;
      }
      currNode = getFullNodeDownPtr(currNode, index);
      if (currNode == 0) break;
    }

    // Element already exists!
    // Return false since do element was added.
    if (currLevel == 0) {
      MEDDLY_DCASSERT(currNode != 0);
      return false;
    }
  }

  // Build a node from the minterm starting at one level below currLevel.
  int unpNode = -1;
  int level;
  for (level=1; level != currLevel; level++) {
#ifdef CREATE_TEMP_NODES_MAX_SIZE_ONLY
    int prime = createTempNodeMaxSize(-level, true);
#else
    int prime = createTempNode(-level, pelement[level] + 1, true);
#endif
    getFullNodeDownPtrs(prime)[pelement[level]] = unpNode;
#ifdef CREATE_TEMP_NODES_MAX_SIZE_ONLY
    unpNode = createTempNodeMaxSize(level, true);
#else
    unpNode = createTempNode(level, element[level] + 1, true);
#endif
    getFullNodeDownPtrs(unpNode)[element[level]] = prime;
  }
  MEDDLY_DCASSERT(currLevel == level);
  MEDDLY_DCASSERT(getNodeLevel(unpNode) == level-1);

  // Deal with the currLevel
  int parentNodeLevel = getNodeLevel(parentNode);

  if (parentNodeLevel == 0) {
    MEDDLY_DCASSERT(tempNode == 0);
    // Build prime node with downpointer to unpNode.
    // Build unprime node with downpointer to prime node.
#ifdef CREATE_TEMP_NODES_MAX_SIZE_ONLY
    int prime = createTempNodeMaxSize(-level, true);
#else
    int prime = createTempNode(-level, pelement[level] + 1, true);
#endif
    getFullNodeDownPtrs(prime)[pelement[level]] = unpNode;
#ifdef CREATE_TEMP_NODES_MAX_SIZE_ONLY
    unpNode = createTempNodeMaxSize(level, true);
#else
    unpNode = createTempNode(level, element[level] + 1, true);
#endif
    getFullNodeDownPtrs(unpNode)[element[level]] = prime;
    // Update the root node.
    tempNode = unpNode;
  }
  else if (parentNodeLevel > 0) {
    // Parent is an unprimed node.
    // Build prime node with downpointer to unpNode.
    // Attach parent to prime node.
#ifdef CREATE_TEMP_NODES_MAX_SIZE_ONLY
    int prime = createTempNodeMaxSize(-level, true);
#else
    int prime = createTempNode(-level, pelement[level] + 1, true);
#endif
    getFullNodeDownPtrs(prime)[pelement[level]] = unpNode;
    MEDDLY_DCASSERT(element[level] < getFullNodeSize(parentNode));
    MEDDLY_DCASSERT(getFullNodeDownPtrs(parentNode)[element[level]] == 0);
    getFullNodeDownPtrs(parentNode)[element[level]] = prime;
  } else {
    // Parent is a primed node.
    // Attach parent to unpNode
    MEDDLY_DCASSERT(pelement[level] < getFullNodeSize(parentNode));
    MEDDLY_DCASSERT(getFullNodeDownPtrs(parentNode)[pelement[level]] == 0);
    getFullNodeDownPtrs(parentNode)[pelement[level]] = unpNode;
  }

  return true;
}
#endif


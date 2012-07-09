
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
#include "../unique_table.h"

// TODO: Add option to create temporary nodes of max size when
//       accumulating minterms.

// ********************************** MTMDDs **********************************

MEDDLY::mtmdd_forest::mtmdd_forest(int dsl, domain *d,
    bool relation, forest::range_type t,
    forest::edge_labeling e, const policies &p)
: mt_forest(dsl, d, relation, t, e, p)
{
  list = 0;
  termList = 0;
  listSize = 0;
  count = 0;
  slot = 0;
  countSize = 0;

  // TBD: this eventually belongs in mt_forest.
  // Initalize level data
  for (int k=getMinLevelIndex(); k<=getNumVariables(); k++) {
    levels[k].init(this, 0, 0, 0);
  }
}



MEDDLY::mtmdd_forest::~mtmdd_forest()
{
  if (list) free(list);
  if (termList) free(termList);
  if (count) free(count);
  if (slot) free(slot);
}


void MEDDLY::mtmdd_forest::expandCountAndSlotArrays(int size)
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


void MEDDLY::mtmdd_forest::resizeNode(int p, int size)
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
  int nodeLevel = getNodeLevel(p);
  int oldDataArraySize = levels[nodeLevel].slotsForNode(nodeSize);
  int newDataArraySize = levels[nodeLevel].slotsForNode(size);
  int newOffset = levels[nodeLevel].allocNode(size, p, false);

  MEDDLY_DCASSERT(newDataArraySize > oldDataArraySize);

  // Pointers to old and new data arrays
  int* prev = getNodeAddress(p);
  int* curr = levels[nodeLevel].data + newOffset;

  // Copy old array to new array
  // Don't copy last location from the old array to the new array
  // -- last location stores the pointer to node p.
  memcpy(curr, prev, (oldDataArraySize - 1) * sizeof(int));

  // Clear out the trailing portion of the new array
  // -- except the last location.
  memset(curr + oldDataArraySize - 1, 0,
      (newDataArraySize - oldDataArraySize) * sizeof(int));

  // Discard the old array
  levels[nodeLevel].recycleNode(address[p].offset);

  // Update the offset field to point to the new array
  address[p].offset = newOffset;

  MEDDLY_DCASSERT(size == getFullNodeSize(p));
}


int MEDDLY::mtmdd_forest::reduceNode(int p)
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
  if (isQuasiReduced()) {
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
  nodeFinder key(this, p);
  int q = unique->find(key);
  if (q) {
    unlinkNode(p);
    return linkNode(q);
  }
  unique->add(key.hash(), p);

#ifdef TRACE_REDUCE
  printf("\tReducing %d: unique, compressing\n", p);
#endif

  if (!areSparseNodesEnabled())
    nnz = size;

  int newoffset = 0;
  // right now, tie goes to truncated full.
  if (levels[node_level].slotsForNode(-nnz) 
    < 
    levels[node_level].slotsForNode(truncsize)) 
  {
    // sparse is better; convert
    newoffset = levels[node_level].allocNode(-nnz, p, false);
    // can't rely on previous ptr, re-point to p
    int* full_ptr = getNodeAddress(p);
    int* sparse_ptr = getAddress(node_level, newoffset);
    // copy first 2 integers: incount, next
    sparse_ptr[0] = full_ptr[0];
    sparse_ptr[1] = full_ptr[1];
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
  } else {
    // full is better
    if (truncsize<size) {
      // truncate the trailing 0s
      newoffset = levels[node_level].allocNode(truncsize, p, false);
      // can't rely on previous ptr, re-point to p
      int* full_ptr = getNodeAddress(p);
      int* trunc_ptr = getAddress(node_level, newoffset);
      // copy first 2 integers: incount, next
      trunc_ptr[0] = full_ptr[0];
      trunc_ptr[1] = full_ptr[1];
      // elements
      memcpy(trunc_ptr + 3, full_ptr + 3, truncsize * sizeof(int));
    }
  }
  // trash old node
  if (newoffset) {
#ifdef MEMORY_TRACE
    int saved_offset = getNodeOffset(p);
    setNodeOffset(p, newoffset);
    levels[node_level].recycleNode(saved_offset);
#else
    levels[node_level].recycleNode(getNodeOffset(p));
    setNodeOffset(p, newoffset);
#endif
  }

  // address[p].cache_count does not change
  MEDDLY_DCASSERT(getCacheCount(p) == 0);
  // Sanity check that the hash value is unchanged
  MEDDLY_DCASSERT(unique->find(key) == p);

  // Temporary node has been transformed to a reduced node; decrement
  // temporary node count.
  decrTempNodeCount(node_level);

  return p;
}


int MEDDLY::mtmdd_forest::createNode(int k, int index, int dptr)
{
  MEDDLY_DCASSERT(index >= -1);

  if (index > -1 && getLevelSize(k) <= index) {
    useExpertDomain()->enlargeVariableBound(k, false, index + 1);
  }

  if (dptr == 0) return 0;
  if (index == -1) {
    // all downpointers should point to dptr
    if (isFullyReduced()) return dptr;
    insertRedundantNode(k, dptr);
    return dptr;
  }

  nodeBuilder& nb = useSparseBuilder(k, 1);
  nb.d(0) = dptr;
  nb.i(0) = index;
  return createReducedNode(-1, nb);
}


void MEDDLY::mtmdd_forest::createEdge(const int* v, int term, dd_edge& e)
{
  // construct the edge bottom-up
  MEDDLY_DCASSERT(isTerminalNode(term));
  int result = term;
  int curr = 0;
  for (int i=1; i<=getExpertDomain()->getNumVariables(); i++) {
    result = createNode(i, v[i], result);
  }
  e.set(result, 0, getNodeLevel(result));
  // e.show(stderr, 2);
}


void
MEDDLY::mtmdd_forest::findFirstElement(const dd_edge& f, int* vlist) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  // vlist is based on level handles.
  int node = f.getNode();
  if (node == 0) 
    throw error(error::INVALID_ASSIGNMENT);

  for (int currLevel = getExpertDomain()->getNumVariables(); currLevel; currLevel--)
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



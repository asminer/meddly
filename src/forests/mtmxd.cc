
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
#include "../unique_table.h"

MEDDLY::mtmxd_forest::mtmxd_forest(int dsl, domain *d,
    bool relation, forest::range_type t,
    forest::edge_labeling e, const policies &p)
: mt_forest(dsl, d, relation, t, e, p)
{
  unpList = 0;
  pList = 0;
  tList = 0;
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


MEDDLY::mtmxd_forest::~mtmxd_forest()
{
  if (unpList) free(unpList);
  if (pList) free(pList);
  if (tList) free(tList);
  if (count) free(count);
  if (slot) free(slot);
}


void MEDDLY::mtmxd_forest::expandCountAndSlotArrays(int size)
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


void MEDDLY::mtmxd_forest::resizeNode(int p, int size)
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


int MEDDLY::mtmxd_forest::reduceNode(int p)
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
      //  int newoffset = levels[node_level].getHole(4+truncsize, true);
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


int MEDDLY::mtmxd_forest::createNode(int k, int index, int dptr)
{
  MEDDLY_DCASSERT(index >= 0 && index < getLevelSize(k) && isValidNodeIndex(dptr));

  if (dptr == 0) return 0;

  nodeBuilder& nb = useSparseBuilder(k, 1);
  nb.d(0) = dptr;
  nb.i(0) = index;
  return createReducedNode(-1, nb);
}

int MEDDLY::mtmxd_forest::createNode(int k, int index1, int index2, int dptr)
{
  MEDDLY_DCASSERT((index1 >= 0 && index2 >= 0) ||
      (index1 >= -1 && index2 >= -1) ||
      (index1 >= -2 && index2 >= -2 && index1 == index2));

  int result = 0;

  if (isIdentityReduced()) {

    if (index1 == -2) {
      // "don't change"
      result = dptr;
    }
    else {
      int p = 0;
      if (index2 == -1) {
        // represents "don't care"
        insertRedundantNode(-k, dptr);
        p = dptr;
      } else {
        p = createNode(-k, index2, dptr);
      }
      if (index1 == -1) {
        // represents "don't care"
        insertRedundantNode(k, p);
        result = p;
      } else {
        result = createNode(k, index1, p);
      }
    }

  }
  else if (index1 == -2) {

    // "don't change"
    MEDDLY_DCASSERT(isQuasiReduced() || isFullyReduced());
    int sz = getLevelSize(k);
    nodeBuilder &nb = useNodeBuilder(k, sz);
    for (int i=0; i<sz; i++) {
      nb.d(i) = createNode(-k, i, dptr);
    }
    result = createReducedNode(-1, nb);
  }
  else if (isQuasiReduced()) {

    int p = 0;
    if (index2 == -1) {
      // represents "don't care"
      insertRedundantNode(-k, dptr);
      p = dptr;
    } else {
      p = createNode(-k, index2, dptr);
    }
    if (index1 == -1) {
      // represents "don't care"
      insertRedundantNode(k, p);
      result = p;
    } else {
      result = createNode(k, index1, p);
    }
  }
  else {

    // deal with "don't care" for primed level
    int p = (index2 == -1) ? dptr : createNode(-k, index2, dptr);
    // deal with "don't care" for unprimed level
    result = (index1 == -1) ? p : createNode(k, index1, p);

  }

  return result;
}


int MEDDLY::mtmxd_forest::createNode(const int* v, const int* vp, int term,
    int startAtHeight, bool primedLevel)
{
  // construct the edge bottom-up
  for (int k=1; k<startAtHeight; k++) {
    term = createNode(k, v[k], vp[k], term);
  } // for k

  // deal with height == startAtHeight
  // handle primed level first
  if (primedLevel) {
    // only primed level to be handled at this height
    term = createNode(-startAtHeight, vp[startAtHeight], term);
  } else {
    // both primed and unprimed levels to be handled at this height
    term = createNode(startAtHeight, v[startAtHeight], vp[startAtHeight], term);
  }

  return term;
}


void MEDDLY::mtmxd_forest::createEdge(const int* v, const int* vp, int term,
    int startAtHeight, bool primedLevel, dd_edge& e)
{
  term = createNode(v, vp, term, startAtHeight, primedLevel);
  e.set(term, 0, getNodeLevel(term));
}


void MEDDLY::mtmxd_forest::createEdge(const int* v, const int* vp, int term,
    dd_edge& e)
{
  createEdge(v, vp, term, getExpertDomain()->getNumVariables(), false, e);
}


int MEDDLY::mtmxd_forest::createEdgeTo(int dptr)
{
  MEDDLY_DCASSERT(isTerminalNode(dptr));
  if (dptr == 0) return 0;

  if (isFullyReduced()) return linkNode(dptr);

  // construct the edge bottom-up
  int curr = dptr;
  for (int i=1; i<=getExpertDomain()->getNumVariables(); i++) {
    curr = createNode(i, -1, -1, curr);
  }
  return curr;
}


int MEDDLY::mtmxd_forest::getTerminalNodeForEdge(int n, const int* vlist,
    const int* vplist) const
{
  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contain a single element.

  if (isIdentityReduced()) {

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


void MEDDLY::mtmxd_forest::findFirstElement(const dd_edge& f,
    int* vlist, int* vplist) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  // vlist is based on level handles.
  int node = f.getNode();
  if (node == 0) 
    throw error(error::INVALID_ASSIGNMENT);

  if (isIdentityReduced()) {

    for (int level = getExpertDomain()->getNumVariables(); level; level--)
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
    for (int currLevel = getExpertDomain()->getNumVariables(); currLevel; )
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



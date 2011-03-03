
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

  
// ******************************** MTMXDs ******************************* 


mtmxd_node_manager::mtmxd_node_manager(domain *d, forest::range_type t)
: node_manager(d, true, t,
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


mtmxd_node_manager::mtmxd_node_manager(domain *d,
    bool relation, forest::range_type t,
    forest::edge_labeling e, forest::reduction_rule r,
    forest::node_storage s, forest::node_deletion_policy dp)
: node_manager(d, relation, t, e, r, s, dp, mtmxdDataHeaderSize)
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


forest::error mtmxd_node_manager::resizeNode(int p, int size)
{
  // This operation can only be performed on Temporary nodes.
  if (!isActiveNode(p) || isTerminalNode(p) || isReducedNode(p)) {
    return forest::INVALID_OPERATION;
  }

  // If node is already large enough, do nothing, and return SUCCESS.
  int nodeSize = getFullNodeSize(p);
  if (size <= nodeSize) return forest::SUCCESS;

  DCASSERT(size > nodeSize);

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

  DCASSERT(newDataArraySize > oldDataArraySize);

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

  DCASSERT(size == getFullNodeSize(p));
  DCASSERT(p == 
      getNodeAddress(p)[getDataHeaderSize() + getFullNodeSize(p) - 1]);

  return forest::SUCCESS;
}


int mtmxd_node_manager::reduceNode(int p)
{
  DCASSERT(isActiveNode(p));
  // assert(reductionRule == forest::IDENTITY_REDUCED);

  if (isReducedNode(p)) return p; 

  DCASSERT(!isTerminalNode(p));
  DCASSERT(isFullNode(p));

  int size = getFullNodeSize(p);
  int* ptr = getFullNodeDownPtrs(p);
  int node_level = getNodeLevel(p);

  decrTempNodeCount(node_level);

#ifdef DEVELOPMENT_CODE
  validateDownPointers(p);
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
      unlinkNode(p);
      p = sharedCopy(q);
    }
    decrTempNodeCount(k);
    return p;
  }
}

int mtmxd_node_manager::createNode(int k, int index1, int index2, int dptr)
{
  DCASSERT((index1 >= 0 && index2 >= 0) ||
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
    DCASSERT(reductionRule == forest::QUASI_REDUCED ||
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
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  DCASSERT(isTerminalNode(term));
  const int* end = h2l_map + startAtHeight;
  for (++h2l_map; h2l_map != end; ++h2l_map)
  {
    int prev = term;
    term = createNode(*h2l_map, v[*h2l_map], vp[*h2l_map], term);
    unlinkNode(prev);
  }

  // deal with height == startAtHeight
  // handle primed level first
  if (primedLevel) {
    // only primed level to be handled at this height
    int prev = term;
    term = createNode(-(*h2l_map), vp[*h2l_map], term);
    unlinkNode(prev);
  } else {
    // both primed and unprimed levels to be handled at this height
    int prev = term;
    term = createNode(*h2l_map, v[*h2l_map], vp[*h2l_map], term);
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


forest::error mtmxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || vplist == 0 || terms == 0 || N <= 0)
    return forest::INVALID_VARIABLE;

  return createEdgeInternal(vlist, vplist, terms, N, e);
}


forest::error mtmxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || vplist == 0 || terms == 0 || N <= 0)
    return forest::INVALID_VARIABLE;

  return createEdgeInternal(vlist, vplist, terms, N, e);
}


int mtmxd_node_manager::createEdge(int dptr)
{
  DCASSERT(isTerminalNode(dptr));
  if (dptr == 0) return 0;

  if (reductionRule == forest::FULLY_REDUCED) return sharedCopy(dptr);

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
        DCASSERT(isTerminalNode(next) || isUnprimedNode(next));
        int currHeight = getNodeHeight(n) - 1;
        int nextHeight = getNodeHeight(next);
        if (nextHeight < currHeight) {
          // skipped levels
          for ( ; nextHeight != currHeight; --currHeight)
          {
            int currLevel = expertDomain->getVariableWithHeight(currHeight);
            DCASSERT(currLevel != -1);
            if (vlist[currLevel] != vplist[currLevel]) {
              next = 0;
              break;
            }
          }
        }
        n = next;
      }
      else {
        // Unprimed Node
        DCASSERT(getDownPtr(n, vlist[nLevel]) == 0
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


forest::error mtmxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  DCASSERT(getRangeType() == forest::INTEGER);
  if (f.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || vplist == 0) return forest::INVALID_VARIABLE;

  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  term = getInteger(getTerminalNodeForEdge(f.getNode(), vlist, vplist));
  return forest::SUCCESS;
}


forest::error mtmxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  DCASSERT(getRangeType() == forest::REAL);
  if (f.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || vplist == 0) return forest::INVALID_VARIABLE;

  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  term = getReal(getTerminalNodeForEdge(f.getNode(), vlist, vplist));
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


forest::error mtmxd_node_manager::findFirstElement(const dd_edge& f,
    int* vlist, int* vplist) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  // vlist is based on level handles.
  int node = f.getNode();
  if (node == 0) return forest::INVALID_ASSIGNMENT;

  if (forest::IDENTITY_REDUCED == reductionRule) {

    for (int currLevel = expertDomain->getTopVariable();
        currLevel != domain::TERMINALS;
        currLevel = expertDomain->getVariableBelow(currLevel))
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
    }

  }
  else {

    // !IDENTITY_REDUCED
    for (int currLevel = expertDomain->getTopVariable();
        currLevel != domain::TERMINALS; )
    {
      DCASSERT(node != 0);
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
                    ? expertDomain->getVariableBelow(-currLevel)
                    : -currLevel;
    }


  }

  return forest::SUCCESS;
}


// *********************************** MXDs *********************************** 


mxd_node_manager::mxd_node_manager(domain *d)
: mtmxd_node_manager(d, true, forest::BOOLEAN,
      forest::MULTI_TERMINAL, forest::IDENTITY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{ }


mxd_node_manager::~mxd_node_manager()
{ }


forest::error mxd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || vplist == 0 || N <= 0) return forest::INVALID_VARIABLE;
  return createEdgeInternal(vlist, vplist, (bool*)0, N, e);
}


forest::error mxd_node_manager::createEdge(bool val, dd_edge &e)
{
  DCASSERT(getRangeType() == forest::BOOLEAN);
  if (e.getForest() != this) return forest::INVALID_OPERATION;

  int node = createEdge(getTerminalNode(val));
  e.set(node, 0, getNodeLevel(node));
  return forest::SUCCESS;
}


forest::error mxd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  DCASSERT(getRangeType() == forest::BOOLEAN);
  if (f.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || vplist == 0) return forest::INVALID_VARIABLE;

  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  term = getBoolean(getTerminalNodeForEdge(f.getNode(), vlist, vplist));
  return forest::SUCCESS;
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



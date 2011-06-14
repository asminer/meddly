
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


#include "evmdd.h"


// ********************************** EVMDDs ********************************** 


evmdd_node_manager
::evmdd_node_manager(domain *d, forest::range_type t, forest::edge_labeling el,
  int dataHeaderSize)
: node_manager(d, false, t,
      el, forest::FULLY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION,
      dataHeaderSize)
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
  address[p].offset = getHole(k, getDataHeaderSize() + 2 * sz, true);
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


void evmdd_node_manager::resizeNode(int p, int size)
{
  // This operation can only be performed on Temporary nodes.
  if (!isActiveNode(p) || isTerminalNode(p) || isReducedNode(p)) {
    throw error(error::INVALID_OPERATION);
  }

  // If node is already large enough, do nothing, and return SUCCESS.
  int oldSize = getFullNodeSize(p);
  if (size <= oldSize) return;

  DCASSERT(size > oldSize);

  // Expand node:
  // (0) Create array of desired size;
  // (1) Copy the header (change the new node size)
  // (2) Copy the downpointers (zero-out the trailing dptrs in the new array)
  // (3) Copy the edge-values (initialize the trailing evs in the new array)
  // (4) Copy the trailer -- part of the header stored at the end of the node
  // (5) Discard the old array
  // (6) Update the offset field

  // (0) Create array of desired size
  int oldDataArraySize = getDataHeaderSize() + 2 * oldSize ;
  int newDataArraySize = getDataHeaderSize() + 2 * size ;
  int nodeLevel = getNodeLevel(p);
  int newOffset = getHole(nodeLevel, newDataArraySize, true);

  DCASSERT(newDataArraySize > oldDataArraySize);

  // Pointers to old and new data arrays
  int* prev = getNodeAddress(p);
  int* curr = level[mapLevel(nodeLevel)].data + newOffset;

  // (1)+(2) Copy the first 3 ints (part of the header) and the downpointers
  memcpy(curr, prev, (3 + oldSize) * sizeof(int));
  // Set the new node size
  curr[2] = size;
  // Advance pointers
  prev += 3 + oldSize;
  curr += 3 + oldSize;
  // Zero-out the trailing dptrs in the new array
  memset(curr, 0, (size - oldSize) * sizeof(int));
  curr += (size - oldSize);

  // (3) Copy the edge-values (initialize the trailing evs in the new array)
  memcpy(curr, prev, oldSize * sizeof(int));
  prev += oldSize;
  curr += oldSize;
  // Initialize trailing edge-values in the new array
  DCASSERT(sizeof(int) == sizeof(float));
  int defaultEV = 0;
  getDefaultEdgeValue(defaultEV);
  for (int* last = curr + (size - oldSize); curr != last; )
  {
    *curr++ = defaultEV;
  }

  // (4) Copy the trailer -- part of the header stored at the end of the node.
  int trailerSize = getDataHeaderSize() - 3;
  memcpy(curr, prev, trailerSize * sizeof(int));
  DCASSERT(p == curr[trailerSize - 1]);

  // (5) Discard the old array
  makeHole(nodeLevel, address[p].offset, oldDataArraySize);

  // (6) Update the offset field
  address[p].offset = newOffset;

  DCASSERT(size == getFullNodeSize(p));
  DCASSERT(p == 
      getNodeAddress(p)[getDataHeaderSize() + 2 * getFullNodeSize(p) - 1]);
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
      if (mid == begin) {
        DCASSERT(begin + 1 == end);
        // find > *mid ==> find > *begin
        // therefore, compare with *end and quit
        // simply advance begin (loop will terminate)
        begin++;
      } else {
        begin = mid;
      }
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


void
evmdd_node_manager::
createEdge(const int* const* vlist, const int* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void
evmdd_node_manager::
createEdge(int val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void
evmdd_node_manager::
evaluate(const dd_edge &f, const int* vlist, int &term) const
{
  throw error(error::INVALID_OPERATION);
}



void evmdd_node_manager::createEdge(const int* const* vlist, int N,
    dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::createEdge(const int* const* vlist,
    const float* terms, int n, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::createEdge(bool val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::createEdge(float val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    float &term) const
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  throw error(error::INVALID_OPERATION);
}


void evmdd_node_manager::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  throw error(error::INVALID_OPERATION);
}


// ********************************* EV+MDDs ********************************** 

evplusmdd_node_manager::evplusmdd_node_manager(domain *d)
: evmdd_node_manager(d, forest::INTEGER, forest::EVPLUS,
  evplusmddDataHeaderSize)
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


int evplusmdd_node_manager::createTempNode(int k, int sz, bool clear)
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

  // Need 5 locations instead of 4 since we are also storing cardinality
  // of index sets
  // address[p].offset = getHole(k, 4 + 2 * sz, true);
  address[p].level = k;
  address[p].offset = getHole(k, getDataHeaderSize() + 2 * sz, true);
  address[p].cache_count = 0;

#ifdef DEBUG_MDD_H
  printf("%s: offset: %d\n", __func__, address[p].offset);
  fflush(stdout);
#endif

  int* foo  = level[mapLevel(k)].data + address[p].offset;
  *foo++    = 1;            // [0]: #incoming
  *foo++    = temp_node;    // [1]:
  *foo++    = sz;           // [2]: size
  foo       += 2 * sz;      // advance to [3 + sz + sz]
  *foo++    = -1;           // [3+sz+sz]: cardinality (-1: not been computed)
  *foo      = p;            // [4+sz+sz]: pointer back to the address array

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

  const int size = getFullNodeSize(p);
  int *dptr = getFullNodeDownPtrs(p);
  int *eptr = getFullNodeEdgeValues(p);
  const int node_level = getNodeLevel(p);

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
    unlinkNode(p);
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
      unlinkNode(p);
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
    unlinkNode(p);
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
    int newoffset = getHole(node_level, getDataHeaderSize()+3*nnz, true);
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
    makeHole(node_level, saved_offset, getDataHeaderSize() + 2 * size);
#else
    makeHole(node_level, getNodeOffset(p), getDataHeaderSize() + 2 * size);
    setNodeOffset(p, newoffset);
#endif
    // address[p].cache_count does not change
  } else {
    // full is better
    if (truncsize < size) {
      // truncate the trailing 0s
      int newoffset = getHole(node_level, getDataHeaderSize()+2*truncsize, true);
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
      makeHole(node_level, saved_offset, getDataHeaderSize() + 2 * size);
#else
      makeHole(node_level, getNodeOffset(p), getDataHeaderSize() + 2 * size);
      setNodeOffset(p, newoffset);
#endif
      // address[p].cache_count does not change
    }
  }

  // address[p].cache_count does not change
  DCASSERT(getCacheCount(p) == 0);
  // Sanity check that the hash value is unchanged
  DCASSERT(find(p) == p);

  // Temporary node has been transformed to a reduced node; decrement
  // temporary node count.
  decrTempNodeCount(node_level);

  return;
}


void evplusmdd_node_manager::getElement(int a, int index, int* e)
{
  if (a == 0) throw error(error::INVALID_VARIABLE);
  if (a == -1) return;

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


void evplusmdd_node_manager::getElement(const dd_edge& a,
    int index, int* e)
{
  if (e == 0 || index < 0) throw error(error::INVALID_VARIABLE);
  e[0] = 0;
  return getElement(a.getNode(), index, e);
}


void
evplusmdd_node_manager::
createEdge(const int* const* vlist, const int* terms, int N, dd_edge &e)
{
  return createEdgeInternal(vlist, terms, N, e);
}


void
evplusmdd_node_manager::
createEdge(int val, dd_edge &e)
{
  return createEdgeInternal(val, e);
}


void
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


void
evplusmdd_node_manager::
createNode(int lh, std::vector<int>& index, std::vector<int>& dptr,
    std::vector<int>& ev, int& result, int& resultEv)
{
#ifdef DEVELOPMENT_CODE
  // Sanity checks:
  // -  index[i] should be in the range [0, levelSize).
  // -  dptr[i] should be a reduced node at a level "below" lh.
  // -  edge-value of INF is reserved for edges pointing to terminal 0.
  for (unsigned i = 0; i < dptr.size(); i++)
  {
    CHECK_RANGE(0, index[i], getLevelSize(lh));
    assert(dptr[i] != 0 && ev[i] != INF);
    assert(isReducedNode(dptr[i]));
    assert(getNodeHeight(dptr[i]) < lh);
  }
#endif

  if (index.size() == 0) {
    result = 0; resultEv = INF;
    return;
  }

  if (isTimeToGc()) { garbageCollect(); }

  incrNodesActivatedSinceGc();

  // Compute minimum edge-value for normalization
  int minEv = INF;
  for (vector<int>::iterator iter = ev.begin(); iter != ev.end(); iter++)
  {
    if (*iter < minEv) minEv = *iter;
  }

  // Normalize edge-values
  for (vector<int>::iterator iter = ev.begin(); iter != ev.end(); iter++)
  {
    *iter -= minEv;
  }

  // ResultEv = minimum edge-value
  resultEv = minEv;

  // Check for possible reductions
  if (int(index.size()) == getLevelSize(lh) &&
      reductionRule == forest::FULLY_REDUCED) {
    // Check for fully-reduced: same dptr[i] and ev[i] == 0
    bool reducible = true;
    if (ev[0] == 0) {
      for (vector<int>::iterator iter = ev.begin(); iter != ev.end(); )
      {
        if (*iter++ != ev[0]) { reducible = false; break; }
      }
      if (reducible) {
        for (vector<int>::iterator iter = dptr.begin(); iter != dptr.end(); )
        {
          if (*iter++ != dptr[0]) { reducible = false; break; }
        }
        if (reducible) {
          // Reduce to dptr[0] and resultEv
          // Unlink all dptr[i], i = 1 to size - 1
          for (vector<int>::iterator iter = dptr.begin() + 1;
              iter != dptr.end(); )
          {
            unlinkNode(*iter++);
          }
          result = dptr[0];
          return;
        }
      }
    }
  }

  // Check if result node will be full or sparse.
  // For that you need to go through the indexes and find the largest index.

#ifdef DEVELOPMENT_CODE
  for (vector<int>::iterator iter = index.begin() + 1;
      iter != index.end(); ++iter)
  {
    assert(*iter > *(iter-1));
  }
#endif

  int largestIndex = index[index.size()-1];
  int fullNodeSize = (largestIndex + 1) * 2 + getDataHeaderSize();
  int sparseNodeSize = index.size() * 3 + getDataHeaderSize();
  int minNodeSize = MIN(fullNodeSize, sparseNodeSize);

  // Get a logical address for result (an index in address[]).
  result = getFreeNode(lh);

  // Fill in address[result].
  address[result].level = lh;
  address[result].offset = getHole(lh, minNodeSize, true);
  address[result].cache_count = 0;

  // Start filling in the actual node data
  int* nodeData = level[mapLevel(lh)].data + address[result].offset;
  nodeData[0] = 1;                      // in-count (# incoming pointers)
  nodeData[1] = getTempNodeId();
  nodeData[minNodeSize - 1] = result;   // pointer back to address[result]

  std::vector<int>::iterator inIter = index.begin();
  std::vector<int>::iterator dpIter = dptr.begin();
  std::vector<int>::iterator evIter = ev.begin();

  if (minNodeSize == fullNodeSize) {
    // Create full node
    // Size is +ve for full-nodes and -ve for sparse nodes.
    nodeData[2] = largestIndex + 1;
    int* resultDp = &nodeData[3];
    int* resultEvs = resultDp + nodeData[2];
    int* last = resultEvs + nodeData[2];

    int currIndex = 0;
    while (inIter != index.end())
    {
      if (currIndex == *inIter) {
        *resultDp++ = *dpIter++;
        *resultEvs++ = *evIter++;
        inIter++;
      } else {
        *resultDp++ = 0;
        *resultEvs++ = INF;
      }
      currIndex++;
    }
    while (resultEvs != last)
    {
      *resultDp++ = 0;
      *resultEvs++ = INF;
    }
  } else {
    // Create sparse node
    // Size is +ve for full-nodes and -ve for sparse nodes.
    nodeData[2] = -index.size();
    int* resultIn = &nodeData[3];
    int* resultDp = resultIn - nodeData[2];
    int* resultEvs = resultDp - nodeData[2];

    inIter = index.begin();
    dpIter = dptr.begin();
    evIter = ev.begin();
    while (dpIter != dptr.end())
    {
      *resultIn++ = *inIter++;
      *resultDp++ = *dpIter++;
      *resultEvs++ = *evIter++;
    }
  }

  // Search in unique table
  int found = find(result);
  if (getNull() == found) {
    // No duplicate found; insert into unique table
    insert(result);
    DCASSERT(getCacheCount(result) == 0);
    DCASSERT(find(result) == result);
  }
  else {
    // Duplicate found; unlink all dptr[] and return the duplicate
    for (dpIter = dptr.begin(); dpIter != dptr.end(); )
    {
      unlinkNode(*dpIter++);
    }
    // Code from deleteTempNode(result) adapted to work here
    {
      makeHole(lh, getNodeOffset(result), minNodeSize);
      freeNode(result);
      if (level[mapLevel(lh)].compactLevel) compactLevel(lh);
    }
    result = sharedCopy(found);
  }
}

// ********************************* EV*MDDs ********************************** 

evtimesmdd_node_manager::evtimesmdd_node_manager(domain *d)
: evmdd_node_manager(d, forest::REAL, forest::EVTIMES,
  evtimesmddDataHeaderSize)
{ }

evtimesmdd_node_manager::~evtimesmdd_node_manager()
{ }

void evtimesmdd_node_manager::initEdgeValues(int p) {
  DCASSERT(!isReducedNode(p));
  DCASSERT(isFullNode(p));

  DCASSERT(sizeof(int) == sizeof(float));

  float defaultEV;
  getDefaultEdgeValue(defaultEV);
  const int intDefaultEV = toInt(defaultEV);
  int *edgeptr = getFullNodeEdgeValues(p);
  for (int *last = edgeptr + getFullNodeSize(p); edgeptr != last; )
    *edgeptr++ = intDefaultEV;
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
      float dev;
      getDefaultEdgeValue(dev);
      evs.resize(size, dev);
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

  const int size = getFullNodeSize(p);
  int *dptr = getFullNodeDownPtrs(p);
  int *eptr = getFullNodeEdgeValues(p);
  float *fptr = (float *)eptr;
  const int node_level = getNodeLevel(p);

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
    unlinkNode(p);
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
        unlinkNode(p);
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
    unlinkNode(p);
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
    int newoffset = getHole(node_level, getDataHeaderSize()+3*nnz, true);
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
    makeHole(node_level, saved_offset, getDataHeaderSize() + 2 * size);
#else
    makeHole(node_level, getNodeOffset(p), getDataHeaderSize() + 2 * size);
    setNodeOffset(p, newoffset);
#endif
    // address[p].cache_count does not change
  } else {
    // full is better
    if (truncsize < size) {
      // truncate the trailing 0s
      int newoffset = getHole(node_level, getDataHeaderSize()+2*truncsize, true);
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
      makeHole(node_level, saved_offset, getDataHeaderSize() + 2 * size);
#else
      makeHole(node_level, getNodeOffset(p), getDataHeaderSize() + 2 * size);
      setNodeOffset(p, newoffset);
#endif
      // address[p].cache_count does not change
    }
  }

  // address[p].cache_count does not change
  DCASSERT(getCacheCount(p) == 0);
  // Sanity check that the hash value is unchanged
  DCASSERT(find(p) == p);

  // Temporary node has been transformed to a reduced node; decrement
  // temporary node count.
  decrTempNodeCount(node_level);

  return;
}


void
evtimesmdd_node_manager::
createEdge(const int* const* vlist, const float* terms, int N, dd_edge &e)
{
  return createEdgeInternal(vlist, terms, N, e);
}


void
evtimesmdd_node_manager::
createEdge(float val, dd_edge &e)
{
  return createEdgeInternal(val, e);
}


void
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


void
evtimesmdd_node_manager::
createNode(int lh, std::vector<int>& index, std::vector<int>& dptr,
    std::vector<float>& ev, int& result, float& resultEv)
{
#ifdef DEVELOPMENT_CODE
  // Sanity checks:
  // -  index[i] should be in the range [0, levelSize).
  // -  dptr[i] should be a reduced node at a level "below" lh.
  // -  edge-value of NAN is reserved for edges pointing to terminal 0.
  for (unsigned i = 0; i < dptr.size(); i++)
  {
    CHECK_RANGE(0, index[i], getLevelSize(lh));
    assert(dptr[i] != 0 && !isNan(ev[i]));
    assert(isReducedNode(dptr[i]));
    assert(getNodeHeight(dptr[i]) < lh);
  }
#endif

  if (index.size() == 0) {
    result = 0; resultEv = NAN;
    return;
  }

  if (isTimeToGc()) { garbageCollect(); }

  incrNodesActivatedSinceGc();

  // Compute minimum edge-value for normalization
  float maxEv = 0.0;
  for (vector<float>::iterator iter = ev.begin(); iter != ev.end(); iter++)
  {
    if (*iter > maxEv) maxEv = *iter;
  }

  // Normalize edge-values
  for (vector<float>::iterator iter = ev.begin(); iter != ev.end(); iter++)
  {
    *iter /= maxEv;
  }

  // ResultEv = minimum edge-value
  resultEv = maxEv;

  // Check for possible reductions
  if (int(index.size()) == getLevelSize(lh) &&
      reductionRule == forest::FULLY_REDUCED) {
    // Check for fully-reduced: same dptr[i] and ev[i] == 0.0
    bool reducible = true;
    if (ev[0] == 0.0) {
      for (vector<float>::iterator iter = ev.begin(); iter != ev.end(); )
      {
        if (*iter++ != ev[0]) { reducible = false; break; }
      }
      if (reducible) {
        for (vector<int>::iterator iter = dptr.begin(); iter != dptr.end(); )
        {
          if (*iter++ != dptr[0]) { reducible = false; break; }
        }
        if (reducible) {
          // Reduce to dptr[0] and resultEv
          // Unlink all dptr[i], i = 1 to size - 1
          for (vector<int>::iterator iter = dptr.begin() + 1;
              iter != dptr.end(); )
          {
            unlinkNode(*iter++);
          }
          result = dptr[0];
          return;
        }
      }
    }
  }

  // Check if result node will be full or sparse.
  // For that you need to go through the indexes and find the largest index.

#ifdef DEVELOPMENT_CODE
  for (vector<int>::iterator iter = index.begin() + 1;
      iter != index.end(); ++iter)
  {
    assert(*iter > *(iter-1));
  }
#endif

  int largestIndex = index[index.size()-1];
  int fullNodeSize = (largestIndex + 1) * 2 + getDataHeaderSize();
  int sparseNodeSize = index.size() * 3 + getDataHeaderSize();
  int minNodeSize = MIN(fullNodeSize, sparseNodeSize);

  // Get a logical address for result (an index in address[]).
  result = getFreeNode(lh);

  // Fill in address[result].
  address[result].level = lh;
  address[result].offset = getHole(lh, minNodeSize, true);
  address[result].cache_count = 0;

  // Start filling in the actual node data
  int* nodeData = level[mapLevel(lh)].data + address[result].offset;
  nodeData[0] = 1;                      // in-count (# incoming pointers)
  nodeData[1] = getTempNodeId();
  nodeData[minNodeSize - 1] = result;   // pointer back to address[result]

  std::vector<int>::iterator inIter = index.begin();
  std::vector<int>::iterator dpIter = dptr.begin();
  std::vector<float>::iterator evIter = ev.begin();

  if (minNodeSize == fullNodeSize) {
    // Create full node
    // Size is +ve for full-nodes and -ve for sparse nodes.
    nodeData[2] = largestIndex + 1;
    int* resultDp = &nodeData[3];
    float* resultEvs = (float*)(resultDp + nodeData[2]);
    float* last = resultEvs + nodeData[2];

    int currIndex = 0;
    while (inIter != index.end())
    {
      if (currIndex == *inIter) {
        *resultDp++ = *dpIter++;
        *resultEvs++ = *evIter++;
        inIter++;
      } else {
        *resultDp++ = 0;
        *resultEvs++ = NAN;;
      }
      currIndex++;
    }
    while (resultEvs != last)
    {
      *resultDp++ = 0;
      *resultEvs++ = NAN;
    }
  } else {
    // Create sparse node
    // Size is +ve for full-nodes and -ve for sparse nodes.
    nodeData[2] = -index.size();
    int* resultIn = &nodeData[3];
    int* resultDp = resultIn - nodeData[2];
    float* resultEvs = (float*)(resultDp - nodeData[2]);

    inIter = index.begin();
    dpIter = dptr.begin();
    evIter = ev.begin();
    while (dpIter != dptr.end())
    {
      *resultIn++ = *inIter++;
      *resultDp++ = *dpIter++;
      *resultEvs++ = *evIter++;
    }
  }

  // Search in unique table
  int found = find(result);
  if (getNull() == found) {
    // No duplicate found; insert into unique table
    insert(result);
    DCASSERT(getCacheCount(result) == 0);
    DCASSERT(find(result) == result);
  }
  else {
    // Duplicate found; unlink all dptr[] and return the duplicate
    for (dpIter = dptr.begin(); dpIter != dptr.end(); )
    {
      unlinkNode(*dpIter++);
    }
    // Code from deleteTempNode(result) adapted to work here
    {
      makeHole(lh, getNodeOffset(result), minNodeSize);
      freeNode(result);
      if (level[mapLevel(lh)].compactLevel) compactLevel(lh);
    }
    result = sharedCopy(found);
  }
}


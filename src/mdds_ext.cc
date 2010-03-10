
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



#include "mdds_ext.h"
#include <algorithm>
#include <limits.h>

//#define DEBUG_SORT_MATRIX
//#define DEBUG_SORT_BUILD

template <typename T>
class bucket { public: int dp; T ev; };


bool checkForReductions(node_manager *f, int p, int nnz, int& result)
{
  if (f->getReductionRule() == forest::QUASI_REDUCED) return false;
  if (nnz != f->getLevelSize(f->getNodeLevel(p))) return false;

  const int* ptr = f->getFullNodeDownPtrs(p);
  int size = f->getFullNodeSize(p);

  switch (f->getReductionRule()) {

    case forest::FULLY_REDUCED:
      result = ptr[0];
      for (int i = 1; i < size; ++i) {
        if (ptr[i] != result) return false;
      }
      break;

    case forest::IDENTITY_REDUCED:
      if (f->isForRelations()) {
        if (f->isPrimedNode(p)) return false;
        if (f->isFullNode(ptr[0])) {
          result = f->getFullNodeDownPtr(ptr[0], 0);
          if (result == 0) return false;
        } else {
          int index = f->getSparseNodeIndex(ptr[0], 0);
          if (index != 0) return false;
          result = f->getSparseNodeDownPtr(ptr[0], 0);
          DCASSERT(result != 0);
        }
        mtmxd_node_manager* xf = smart_cast<mtmxd_node_manager*>(f);
        if (xf == 0) {
          printf("Casting %p to mtmxd_node_manager failed. Terminating\n.", f);
          exit(1);
        }
        for (int i = 0; i < size; i++) {
          if (!xf->singleNonZeroAt(ptr[i], result, i)) return false;
        }
      }
      else {
        printf("Identity-Reduction is valid only for forests that ");
        printf("store relations.\n");
        printf("Either change reduction rule for forest %p or enable\n", f);
        printf("relations for it.\n");
        printf("Terminating.\n");
        exit(1);
      }
      break;

    default:
      return false;
  }

  return true;
}


void validateDownPointers(node_manager* f, int p)
{
  assert(f->isFullNode(p));
  int nodeHeight = f->getNodeHeight(p);
  int nodeLevel = f->getNodeLevel(p);
  int nodeSize = f->getFullNodeSize(p);
  const int* ptr = f->getFullNodeDownPtrs(p);
  switch (f->getReductionRule()) {
    case forest::FULLY_REDUCED:
      if (f->isUnprimedNode(p)) {
        // unprimed node
        for (int i = 0; i < nodeSize; ++i) {
          assert(f->isReducedNode(ptr[i]));
          assert(!f->isForRelations() ||
              f->isTerminalNode(ptr[i]) ||
              f->getNodeHeight(ptr[i]) < nodeHeight ||
              f->getNodeLevel(ptr[i]) == -nodeLevel);
        }
      } else {
        // primed node
        for (int i = 0; i < nodeSize; ++i) {
          assert(f->isReducedNode(ptr[i]));
          assert(f->isTerminalNode(ptr[i]) ||
              f->getNodeHeight(ptr[i]) < nodeHeight);
        }
      }
      break;

    case forest::QUASI_REDUCED:
      if (f->isUnprimedNode(p)) {
        // unprimed node
        for (int i = 0; i < nodeSize; ++i) {
          assert(f->isReducedNode(ptr[i]));
          assert(!f->isForRelations() ||
              f->isTerminalNode(ptr[i]) ||
              f->getNodeLevel(ptr[i]) == -nodeLevel);
        }
      } else {
        // primed node
        for (int i = 0; i < nodeSize; ++i) {
          assert(f->isReducedNode(ptr[i]));
          assert(f->isTerminalNode(ptr[i]) ||
              (f->getNodeHeight(ptr[i]) == (nodeHeight - 1) &&
               f->isUnprimedNode(ptr[i])));
        }
      }
      break;

    case forest::IDENTITY_REDUCED:
      assert(f->isForRelations());
      if (f->isUnprimedNode(p)) {
        // unprimed node
        for (int i = 0; i < nodeSize; ++i) {
          assert(f->isReducedNode(ptr[i]));
          assert(ptr[i] == 0 || (f->getNodeLevel(ptr[i]) == -nodeLevel));
        }
      } else {
        // primed node
        for (int i = 0; i < nodeSize; ++i) {
          assert(f->isReducedNode(ptr[i]));
          assert(f->getNodeHeight(ptr[i]) < nodeHeight);
          assert(f->isTerminalNode(ptr[i]) || f->isUnprimedNode(ptr[i]));
        }
      }
      break;

    default:
      break;
  }
}


// ********************************** MTMDDs **********************************

mtmdd_node_manager::mtmdd_node_manager(domain *d, forest::range_type t)
: node_manager(d, false, t,
      forest::MULTI_TERMINAL, forest::FULLY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{ }


mtmdd_node_manager::mtmdd_node_manager(domain *d,
    bool relation, forest::range_type t,
    forest::edge_labeling e, forest::reduction_rule r,
    forest::node_storage s, forest::node_deletion_policy dp)
: node_manager(d, relation, t, e, r, s, dp)
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

  int size = getFullNodeSize(p);
  int* ptr = getFullNodeDownPtrs(p);
  int node_level = getNodeLevel(p);

  decrTempNodeCount(node_level);

#ifdef DEVELOPMENT_CODE
  validateDownPointers(this, p);
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
      DCASSERT(ptr[i] == 0 || (getNodeLevel(ptr[i]) == nextLevel));
    }
  } else {
    // check for possible reductions
    int temp = 0;
    if (checkForReductions(this, p, nnz, temp)) {
      linkNode(temp);
      deleteTempNode(p);
      return temp;
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
  DCASSERT(index >= -1);

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

  return createEdgeInternal(vlist, terms, N, e);
}


forest::error mtmdd_node_manager::createEdgeHelper(int terminalNode, dd_edge& e)
{
  DCASSERT(isTerminalNode(terminalNode));

  if (reductionRule == forest::FULLY_REDUCED || terminalNode == 0) {
    e.set(terminalNode, 0, domain::TERMINALS);
    return forest::SUCCESS;
  }

  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;
  int result = terminalNode;
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


forest::error mtmdd_node_manager::createEdge(int term, dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  return createEdgeHelper(getTerminalNode(term), e);
}


int mtmdd_node_manager::getTerminalNodeForEdge(int n, const int* vlist) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  while (!isTerminalNode(n)) {
    n = getDownPtr(n, vlist[h2l_map[getNodeHeight(n)]]);
  }
  return n;
}


forest::error mtmdd_node_manager::evaluate(const dd_edge &f,
    const int* vlist, int &term) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  term = getInteger(getTerminalNodeForEdge(f.getNode(), vlist));
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

  return createEdgeInternal(vlist, terms, N, e);
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
  return createEdgeHelper(getTerminalNode(term), e);
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
  term = getReal(getTerminalNodeForEdge(f.getNode(), vlist));
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


forest::error
mtmdd_node_manager::findFirstElement(const dd_edge& f, int* vlist) const
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



// *********************************** MDDs ***********************************

mdd_node_manager::mdd_node_manager(domain *d)
: mtmdd_node_manager(d, false, forest::BOOLEAN,
      forest::MULTI_TERMINAL, forest::FULLY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{ }


mdd_node_manager::~mdd_node_manager()
{ }


forest::error mdd_node_manager::createEdge(const int* const* vlist,
    int N, dd_edge &e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || N <= 0) return forest::INVALID_VARIABLE;
  return mtmdd_node_manager::createEdgeInternal(vlist, (bool*)0, N, e);
}


forest::error mdd_node_manager::createEdge(bool term, dd_edge& e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  return createEdgeHelper(getTerminalNode(term), e);
}


forest::error mdd_node_manager::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  term = getBoolean(getTerminalNodeForEdge(f.getNode(), vlist));
  return forest::SUCCESS;
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



// ******************************** MTMXDs ******************************* 


mtmxd_node_manager::mtmxd_node_manager(domain *d, forest::range_type t)
: node_manager(d, true, t,
      forest::MULTI_TERMINAL, forest::IDENTITY_REDUCED,
      forest::FULL_OR_SPARSE_STORAGE, OPTIMISTIC_DELETION)
{ }


mtmxd_node_manager::mtmxd_node_manager(domain *d,
    bool relation, forest::range_type t,
    forest::edge_labeling e, forest::reduction_rule r,
    forest::node_storage s, forest::node_deletion_policy dp)
: node_manager(d, relation, t, e, r, s, dp)
{ }


mtmxd_node_manager::~mtmxd_node_manager()
{ }


bool mtmxd_node_manager::singleNonZeroAt(int p, int val, int index) const
{
  DCASSERT(isActiveNode(p));
  DCASSERT(!isTerminalNode(p));
  DCASSERT(!isZombieNode(p));
  DCASSERT(val != 0);
  if (isFullNode(p)) {
    const int* dptr = getFullNodeDownPtrsReadOnly(p);
    const int sz = getFullNodeSize(p);
    if (index >= sz || dptr[index] != val) return false;
    int i = 0;
    for ( ; i < index; ++i) { if (dptr[i] != 0) return false; }
    for (i = index + 1 ; i < sz; ++i) { if (dptr[i] != 0) return false; }
  } else {
    if (getSparseNodeSize(p) != 1) return false;
    if (getSparseNodeIndex(p, 0) != index) return false;
    if (getSparseNodeDownPtr(p, 0) != val) return false;
  }
  return true;
}


int mtmxd_node_manager::reduceNode(int p)
{
  DCASSERT(isActiveNode(p));
  // assert(reductionRule == forest::IDENTITY_REDUCED);

  if (isReducedNode(p)) return p; 

  DCASSERT(!isTerminalNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(getInCount(p) == 1);

  int size = getFullNodeSize(p);
  int* ptr = getFullNodeDownPtrs(p);
  int node_level = getNodeLevel(p);

  decrTempNodeCount(node_level);

#ifdef DEVELOPMENT_CODE
  validateDownPointers(this, p);
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
  int temp = 0;
  if (checkForReductions(this, p, nnz, temp)) {
    linkNode(temp);
    deleteTempNode(p);
    return temp;
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


void mtmxd_node_manager::createEdge(const int* v, const int* vp, int term,
    dd_edge& e)
{
  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;
  DCASSERT(isTerminalNode(term));
  const int* end = h2l_map + h_sz;
  for (++h2l_map; h2l_map != end; ++h2l_map)
  {
    int prev = term;
    term = createNode(*h2l_map, v[*h2l_map], vp[*h2l_map], term);
    unlinkNode(prev);
  }
  e.set(term, 0, getNodeLevel(term));
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
  const int nodeHeight = expertDomain->getVariableHeight(lh);
  for (unsigned i = 0; i < dptr.size(); i++)
  {
    CHECK_RANGE(0, index[i], getLevelSize(lh));
    assert(dptr[i] != 0 && ev[i] != INF);
    assert(isReducedNode(dptr[i]));
    assert(getNodeHeight(dptr[i]) < nodeHeight);
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
  int fullNodeSize = (largestIndex + 1) * 2 + 4;
  int sparseNodeSize = index.size() * 3 + 4;
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
: evmdd_node_manager(d, forest::EVTIMES)
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
  DCASSERT(getInCount(p) == 1);

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
  const int nodeHeight = expertDomain->getVariableHeight(lh);
  for (unsigned i = 0; i < dptr.size(); i++)
  {
    CHECK_RANGE(0, index[i], getLevelSize(lh));
    assert(dptr[i] != 0 && !isNan(ev[i]));
    assert(isReducedNode(dptr[i]));
    assert(getNodeHeight(dptr[i]) < nodeHeight);
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
  int fullNodeSize = (largestIndex + 1) * 2 + 4;
  int sparseNodeSize = index.size() * 3 + 4;
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
  DCASSERT(a->getPtr() != INT_MAX);
  return a->pri[a->getPtr()] < b->pri[a->getPtr()];
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

  // no need to sort terms

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
  const int index = a->getPtr();
  DCASSERT(index != INT_MAX);
  return index < 0?
    // prime
    a->ppri[-index] < b->ppri[-index]:
    // unprime
    a->pri[index] < b->pri[index];
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
}


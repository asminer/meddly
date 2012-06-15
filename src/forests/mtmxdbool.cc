
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


#include "mtmxdbool.h"


#define CREATE_TEMP_NODES_MAX_SIZE_ONLY

#define REDUCE_ACCUMULATE_MXD
  

MEDDLY::mt_mxd_bool::mt_mxd_bool(int dsl, domain *d, const policies &p)
: MEDDLY::mtmxd_forest(dsl, d, true, BOOLEAN, MULTI_TERMINAL, p)
{ }


MEDDLY::mt_mxd_bool::~mt_mxd_bool()
{ }


void MEDDLY::mt_mxd_bool::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0 || N <= 0) 
    throw error(error::INVALID_VARIABLE);
  createEdgeInternal(vlist, vplist, (bool*)0, N, e);
}


void MEDDLY::mt_mxd_bool::createEdge(bool val, dd_edge &e)
{
  MEDDLY_DCASSERT(getRangeType() == forest::BOOLEAN);
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);

  int node = createEdge(getTerminalNode(val));
  e.set(node, 0, getNodeLevel(node));
}


void MEDDLY::mt_mxd_bool::evaluate(const dd_edge& f, const int* vlist,
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


void MEDDLY::mt_mxd_bool::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void MEDDLY::mt_mxd_bool::createEdge(const int* const* vlist,
    const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void MEDDLY::mt_mxd_bool::createEdge(int val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void MEDDLY::mt_mxd_bool::createEdge(float val, dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}


void MEDDLY::mt_mxd_bool::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  throw error(error::INVALID_OPERATION);
}


void MEDDLY::mt_mxd_bool::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, float &term) const
{
  throw error(error::INVALID_OPERATION);
}


#if 1

int MEDDLY::mt_mxd_bool::buildQRIdentityNode(int node, int level)
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

int MEDDLY::mt_mxd_bool::accumulateExpandA(int a, int b, bool cBM)
{
  // a[i][i] += b
  // return reduceNode(a)

  MEDDLY_DCASSERT(isIdentityReduced());
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


void MEDDLY::mt_mxd_bool::accumulateMxdHelper(int& a, int b, bool cBM,
    bool needsToMakeACopy,
    int (MEDDLY::mt_mxd_bool::*function)(int, int, bool))
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


int MEDDLY::mt_mxd_bool::accumulateMxdPrime(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(isIdentityReduced());
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
      &MEDDLY::mt_mxd_bool::accumulateMxd);

  return reduceNode(a);
}


int MEDDLY::mt_mxd_bool::accumulateMxd(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(isIdentityReduced());
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
      &MEDDLY::mt_mxd_bool::accumulateMxdPrime);

  return reduceNode(a);
}


#else

int MEDDLY::mt_mxd_bool::accumulateExpandA(int a, int b, bool cBM)
{
  // a[i][i] += b

  MEDDLY_DCASSERT(isIdentityReduced());
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


void MEDDLY::mt_mxd_bool::accumulateMxdHelper(int& a, int b, bool cBM,
    bool needsToMakeACopy,
    int (MEDDLY::mt_mxd_bool::*function)(int, int, bool))
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
int MEDDLY::mt_mxd_bool::accumulateMxdPrime(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(isIdentityReduced());
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
      &MEDDLY::mt_mxd_bool::accumulateMxd);

  return savedTempNode == a? sharedCopy(a): a;
}


// TODO: Make temporary nodes of max size.
int MEDDLY::mt_mxd_bool::accumulateMxd(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(isIdentityReduced());
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
      &MEDDLY::mt_mxd_bool::accumulateMxdPrime);

  return savedTempNode == a? sharedCopy(a): a;
}


#endif


int MEDDLY::mt_mxd_bool::addPrimeReducedNodes(int a, int b)
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


void MEDDLY::mt_mxd_bool::accumulate(int& a, int b)
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


int MEDDLY::mt_mxd_bool::accumulateSkippedLevel(int tempNode,
    int* element, int* pelement, int level)
{
  MEDDLY_DCASSERT(isIdentityReduced());
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
int MEDDLY::mt_mxd_bool::accumulate(int tempNode, bool cBM,
    int* element, int* pelement, int level)
{
  MEDDLY_DCASSERT(isMxd());
  MEDDLY_DCASSERT(isIdentityReduced());

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
bool MEDDLY::mt_mxd_bool::accumulate(int& tempNode,
    int* element, int* pelement)
{
  assert(isActiveNode(tempNode));
  assert(element != 0);
  assert(pelement != 0);

  // Enlarge variable bounds if necessary
  int level = domain::TERMINALS;
  while (-1 != (level = getExpertDomain()->getVariableAbove(level))) {
    int sz = MAX( element[level] , pelement[level] ) + 1;
    if (sz > getExpertDomain()->getVariableBound(level)) {
      getExpertDomain()->enlargeVariableBound(level, false, sz);
    }
  }

  accumulateMintermAddedElement = false;
  int result = accumulate(tempNode, false,
      element, pelement, getExpertDomain()->getNumVariables());
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
bool MEDDLY::mt_mxd_bool::accumulate(int& tempNode,
    int* element, int* pelement)
{
  assert(isActiveNode(tempNode));
  assert(element != 0);
  assert(pelement != 0);
  MEDDLY_DCASSERT(tempNode == 0 || !isReducedNode(tempNode));

  // Enlarge variable bounds if necessary
#if 0
  int level = domain::TERMINALS;
  while (-1 != (level = getExpertDomain()->getVariableAbove(level))) {
    int sz = MAX( element[level] , pelement[level] ) + 1;
    if (sz > getExpertDomain()->getVariableBound(level)) {
      getExpertDomain()->enlargeVariableBound(level, false, sz);
    }
  }
#endif

  // Traverse Mdd till you find a 0.
  int parentNode = 0;
  int currNode = tempNode;
  int currLevel = getExpertDomain()->getNumVariables();
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


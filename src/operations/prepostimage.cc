
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


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "prepostimage.h"
#include "../compute_table.h"

namespace MEDDLY {
  class image_op;

  class preimage_mdd;
  class postimage_mdd;

  class preimage_opname;
  class postimage_opname;
};

// ******************************************************************
// *                                                                *
// *                         image_op class                         *
// *                                                                *
// ******************************************************************

/// Abstract base class for all pre/post image operations.
class MEDDLY::image_op : public binary_operation {
  public:
    image_op(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;

    inline bool findResult(int a, int b, int &c) {
      static int key[2];
      key[0] = a; key[1] = b;
      const int* cacheFind = CT->find(this, key);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[2]);
      return true;
    }
    inline int saveResult(int a, int b, int c) {
      static int cacheEntry[3];
      cacheEntry[0] = arg1F->cacheNode(a); 
      cacheEntry[1] = arg2F->cacheNode(b);
      cacheEntry[2] = resF->cacheNode(c);
      CT->add(this, cacheEntry);
      return c;
    }
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual int compute(int a, int b) = 0;
  protected:
    binary_operation* unionOp;
};

MEDDLY::image_op::image_op(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : binary_operation(oc, a1, a2, res)
{
  key_length = 2; 
  ans_length = 1;
  unionOp = 0;
}

bool MEDDLY::image_op::isEntryStale(const int* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}

void MEDDLY::image_op::discardEntry(const int* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void
MEDDLY::image_op::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d, %d): %d]", getName(), data[0], data[1], data[2]);
}

void MEDDLY::image_op
::compute(const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  if (resF->getRangeType() == forest::BOOLEAN) {
    unionOp = getOperation(UNION, a, b, c);
  } else {
    unionOp = getOperation(MAXIMUM, a, b, c);
  }
  DCASSERT(unionOp);
  int cnode = compute(a.getNode(), b.getNode());
  c.set(cnode, 0, resF->getNodeLevel(cnode));
  unionOp = 0;  // for sanity
}


// ******************************************************************
// *                                                                *
// *                       preimage_mdd class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::preimage_mdd : public image_op {
  public:
    preimage_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual int compute(int a, int b);
  protected:
    int expandMxdByOneLevel(int a, int b);
    int expandMxd(int a, int b);
    int expandMdd(int a, int b);
    int expandByOneLevel(int a, int b);
    int expand(int a, int b);
};

MEDDLY::preimage_mdd::preimage_mdd(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : image_op(oc, a1, a2, res)
{
}

int MEDDLY::preimage_mdd::compute(int mdd, int mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (arg2F->isTerminalNode(mxd)) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->getTerminalNode(true);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  int result = 0;
  if (findResult(mdd, mxd, result)) {
    return result;
  }

  // check if mxd and mdd are at the same level
  int mddLevel = arg1F->getNodeLevel(mdd);
  int mxdLevel = arg2F->getNodeLevel(mxd);

  if (mddLevel < mxdLevel) {
    // expand mxd
    result = expandMxd(mdd, mxd);
  } else if (mddLevel > mxdLevel) {
    // expand mdd
    result = expandMdd(mdd, mxd);
  } else {
    // same level
    DCASSERT(arg1F->getNodeLevel(mdd) == arg2F->getNodeLevel(mxd));
    result = expand(mdd, mxd);
  } // same level

  DCASSERT(resF->isReducedNode(result));

  // save result in compute cache and return it
  return saveResult(mdd, mxd, result);
}

// result += pre_image(mdd, iMxd[j])
int MEDDLY::preimage_mdd::expandMxdByOneLevel(int mdd, int iMxd)
{
  if (iMxd == 0) return 0;

  DCASSERT(!arg2F->isTerminalNode(iMxd));

  int result = 0;
  int tempResult = 0;
  if (arg2F->isFullNode(iMxd)) {
    int iMxdSize = arg2F->getFullNodeSize(iMxd);
    for (int j = 0; j < iMxdSize; j++) {
      int ijMxd = arg2F->getFullNodeDownPtr(iMxd, j);
      if (ijMxd == 0) continue;
      int jResult = compute(mdd, ijMxd);
      if (jResult == 0) continue;
      tempResult = result;
      result = unionOp->compute(result, jResult);
      resF->unlinkNode(jResult);
      resF->unlinkNode(tempResult);
    }
  } else {
    // iMxd is a sparse node
    int iMxdSize = arg2F->getSparseNodeSize(iMxd);
    for (int j = 0; j < iMxdSize; j++) {
      int jResult = compute(mdd, arg2F->getSparseNodeDownPtr(iMxd, j));
      if (jResult == 0) continue;
      tempResult = result;
      result = unionOp->compute(result, jResult);
      resF->unlinkNode(jResult);
      resF->unlinkNode(tempResult);
    }
  }
  return result;
}

// result[j] += pre_image(mdd[i], mxd[j][i])
// becomes
// result[j] += pre_image(mdd, mxd[j][i])
int MEDDLY::preimage_mdd::expandMxd(int mdd, int mxd)
{
  // create new node
  int result = 0;
  if (arg2F->isFullNode(mxd)) {
    // full mxd node
    int mxdSize = arg2F->getFullNodeSize(mxd);
    result = resF->createTempNode(arg2F->getNodeLevel(mxd), mxdSize, false);
    for (int i = 0; i < mxdSize; i++) {
      int temp = expandMxdByOneLevel(mdd, arg2F->getFullNodeDownPtr(mxd, i));
      resF->setDownPtrWoUnlink(result, i, temp);
      resF->unlinkNode(temp);
    } // for mxdSize
  } else {
    // sparse mxd node
    int mxdSize = arg2F->getSparseNodeSize(mxd);
    result = resF->createTempNode(arg2F->getNodeLevel(mxd),
        arg2F->getSparseNodeIndex(mxd, mxdSize-1) + 1
    );
    for (int i = 0; i < mxdSize; i++) {
      int tR = expandMxdByOneLevel(mdd, arg2F->getSparseNodeDownPtr(mxd, i));
      resF->setDownPtrWoUnlink(result, arg2F->getSparseNodeIndex(mxd, i), tR);
      resF->unlinkNode(tR);
    } // for mxdSize
  }
  return resF->reduceNode(result);
}

int MEDDLY::preimage_mdd::expandMdd(int mdd, int mxd)
{
  // create new node
  int result = resF->createTempNodeMaxSize(arg1F->getNodeLevel(mdd));
  if (arg1F->isFullNode(mdd)) {
    // mdd is a full node
    int mddSize = arg1F->getFullNodeSize(mdd);
    for (int i = 0; i < mddSize; i++) {
      int tempResult = compute(arg1F->getFullNodeDownPtr(mdd, i), mxd);
      resF->setDownPtr(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
  } else {
    // mdd is a sparse node
    int mddSize = arg1F->getSparseNodeSize(mdd);
    for (int i = 0; i < mddSize; i++) {
      int tempResult = compute(arg1F->getSparseNodeDownPtr(mdd, i), mxd);
      resF->setDownPtr(result, arg1F->getSparseNodeIndex(mdd, i), tempResult);
      resF->unlinkNode(tempResult);
    }
  }
  return resF->reduceNode(result);
}

// result += pre_image(mdd[j], iMxd[j])
int MEDDLY::preimage_mdd::expandByOneLevel(int mdd, int iMxd)
{
  if (iMxd == 0) return 0;
  DCASSERT(!arg2F->isTerminalNode(iMxd));
  int result = 0;
  int tempResult = 0;
  if (arg2F->isFullNode(iMxd)) {
    if (arg1F->isFullNode(mdd)) {
      int minSize =
        MIN(arg1F->getFullNodeSize(mdd), arg2F->getFullNodeSize(iMxd));
      for (int j = 0; j < minSize; j++) {
        int ijMxd = arg2F->getFullNodeDownPtr(iMxd, j);
        int jMdd = arg1F->getFullNodeDownPtr(mdd, j);
        if (ijMxd == 0 || jMdd == 0) continue;
        int jResult = compute(jMdd, ijMxd);
        if (jResult == 0) continue;
        tempResult = result;
        result = unionOp->compute(result, jResult);
        resF->unlinkNode(jResult);
        resF->unlinkNode(tempResult);
      }
    } // mdd is a full node
    else {
      DCASSERT(arg1F->isSparseNode(mdd));
      int iMxdSize = arg2F->getFullNodeSize(iMxd);
      int mddSize = arg1F->getSparseNodeSize(mdd);
      for (int j = 0; j < mddSize; j++) {
        int jIndex = arg1F->getSparseNodeIndex(mdd, j);
        if (jIndex >= iMxdSize) break;
        int ijMxd = arg2F->getFullNodeDownPtr(iMxd, jIndex);
        if (ijMxd == 0) continue;
        int jMdd = arg1F->getSparseNodeDownPtr(mdd, j);
        int jResult = compute(jMdd, ijMxd);
        if (jResult == 0) continue;
        tempResult = result;
        result = unionOp->compute(result, jResult);
        resF->unlinkNode(jResult);
        resF->unlinkNode(tempResult);
      }
    } // mdd is a sparse node
  } // iMxd is a full node
  else {
    DCASSERT(arg2F->isSparseNode(iMxd));
    if (arg1F->isFullNode(mdd)) {
      int iMxdSize = arg2F->getSparseNodeSize(iMxd);
      int mddSize = arg1F->getFullNodeSize(mdd);
      for (int j = 0; j < iMxdSize; j++) {
        int jIndex = arg2F->getSparseNodeIndex(iMxd, j);
        if (jIndex >= mddSize) break;
        int jMdd = arg1F->getFullNodeDownPtr(mdd, jIndex);
        if (jMdd == 0) continue;
        int ijMxd = arg2F->getSparseNodeDownPtr(iMxd, j);
        int jResult = compute(jMdd, ijMxd);
        if (jResult == 0) continue;
        tempResult = result;
        result = unionOp->compute(result, jResult);
        resF->unlinkNode(jResult);
        resF->unlinkNode(tempResult);
      }
    } // mdd is a full node
    else {
      DCASSERT(arg1F->isSparseNode(mdd));
      int iMxdLargestIndex = arg2F->getSparseNodeIndex(iMxd,
        arg2F->getSparseNodeSize(iMxd) - 1);
      int mddSize = arg1F->getSparseNodeSize(mdd);
      // k: traverses iMxd, j: traverses mdd
      int j = 0;
      int k = 0;
      int iMxdIndex = arg2F->getSparseNodeIndex(iMxd, k);
      for ( ; j < mddSize; ) {
        int mddIndex = arg1F->getSparseNodeIndex(mdd, j);
        if (mddIndex > iMxdLargestIndex) break;
        while (iMxdIndex < mddIndex) {
          DCASSERT((k+1) < arg2F->getSparseNodeSize(iMxd));
          iMxdIndex = arg2F->getSparseNodeIndex(iMxd, ++k);
        }
        if (mddIndex < iMxdIndex) { ++j; continue; }
        DCASSERT(mddIndex == iMxdIndex);
        int jResult = compute(
            arg1F->getSparseNodeDownPtr(mdd, j),
            arg2F->getSparseNodeDownPtr(iMxd, k)
        );
        ++j; ++k;
        if (jResult == 0) continue;
        tempResult = result;
        result = unionOp->compute(result, jResult);
        resF->unlinkNode(jResult);
        resF->unlinkNode(tempResult);
      }
    } // mdd is a sparse node
  } // iMxd is a sparse node

  return result;
}

// result[j] += pre_image(mdd[i], mxd[j][i])
int MEDDLY::preimage_mdd::expand(int mdd, int mxd)
{
  // create new node
  int result = 0;
  if (arg2F->isFullNode(mxd)) {
    // full mxd node
    int mxdSize = arg2F->getFullNodeSize(mxd);
    result = resF->createTempNode(arg2F->getNodeLevel(mxd), mxdSize, false);
    for (int i = 0; i < mxdSize; i++) {
      int tR = expandByOneLevel(mdd, arg2F->getFullNodeDownPtr(mxd, i));
      resF->setDownPtrWoUnlink(result, i, tR);
      resF->unlinkNode(tR);
    } // for mxdSize
  } else {
    // sparse mxd node
    int mxdSize = arg2F->getSparseNodeSize(mxd);
    result = resF->createTempNode(arg2F->getNodeLevel(mxd),
        arg2F->getSparseNodeIndex(mxd, mxdSize - 1) + 1
    );
    for (int i = 0; i < mxdSize; i++) {
      int tR = expandByOneLevel(mdd, arg2F->getSparseNodeDownPtr(mxd, i));
      resF->setDownPtrWoUnlink(result, arg2F->getSparseNodeIndex(mxd, i), tR);
      resF->unlinkNode(tR);
    } // for mxdSize
  }
  return resF->reduceNode(result);
}



// ******************************************************************
// *                                                                *
// *                      postimage_mdd  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::postimage_mdd : public image_op {
  public:
    postimage_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual int compute(int a, int b);
  protected:
    int fullFull(int mdd, int mxd);
    int fullSparse(int mdd, int mxd);
    int sparseFull(int mdd, int mxd);
    int sparseSparse(int mdd, int mxd);
    int expandMdd(int mdd, int mxd);
    int expandMxd(int mdd, int mxd);
    void expandMxdHelper(int mdd, int iMxd, int result);
};

MEDDLY::postimage_mdd::postimage_mdd(const binary_opname* oc, 
  expert_forest* a1, expert_forest* a2, expert_forest* res)
: image_op(oc, a1, a2, res)
{
}

int MEDDLY::postimage_mdd::compute(int mdd, int mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;

  // skipped level means identity matrix
  if (arg2F->isTerminalNode(mxd)) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->getTerminalNode(true);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  int result = 0;
  if (findResult(mdd, mxd, result)) {
    return result;
  }

  // check if mxd and mdd are at the same level
  int mddLevel = arg1F->getNodeLevel(mdd);
  int mxdLevel = arg2F->getNodeLevel(mxd);

  if (mddLevel < mxdLevel) {
    // expand mxd
    result = expandMxd(mdd, mxd);
  } else if (mddLevel > mxdLevel) {
    // expand mdd
    result = expandMdd(mdd, mxd);
  } else {
    // same level
    DCASSERT(arg1F->getNodeLevel(mdd) == arg2F->getNodeLevel(mxd));
    if (arg1F->isFullNode(mdd)) {
      if (arg2F->isFullNode(mxd))
        result = fullFull(mdd, mxd);
      else
        result = fullSparse(mdd, mxd);
    } else {
      if (arg2F->isFullNode(mxd))
        result = sparseFull(mdd, mxd);
      else
        result = sparseSparse(mdd, mxd);
    }
  } // same level

  DCASSERT(resF->isReducedNode(result));

  // save result in compute cache and return it
  return saveResult(mdd, mxd, result);
}

int MEDDLY::postimage_mdd::fullFull(int mdd, int mxd)
{
  DCASSERT(arg1F->isFullNode(mdd));
  DCASSERT(arg2F->isFullNode(mxd));

  // for each mxd[i]
  //   for each mxd[i][j]
  //     dptrs[j] += postimage(mdd[i],mxd[i][j])

  // create new node
  int minSize = MIN(arg1F->getFullNodeSize(mdd), arg2F->getFullNodeSize(mxd));
  // need to create result of max size since post image computation may
  // write to result[i] where i >= minSize
  int result = resF->createTempNodeMaxSize(arg1F->getNodeLevel(mdd));

  int iMdd = 0;
  int iMxd = 0;
  for (int i = 0; i < minSize; i++) {
    iMdd = arg1F->getFullNodeDownPtr(mdd, i);
    if (iMdd == 0) continue;
    iMxd = arg2F->getFullNodeDownPtr(mxd, i);
    expandMxdHelper(iMdd, iMxd, result);
  }
  return resF->reduceNode(result);
}


int MEDDLY::postimage_mdd::fullSparse(int mdd, int mxd)
{
  DCASSERT(arg1F->isFullNode(mdd));
  DCASSERT(arg2F->isSparseNode(mxd));

  // create new node
  int result = resF->createTempNodeMaxSize(arg1F->getNodeLevel(mdd));
  int mxdSize = arg2F->getSparseNodeSize(mxd);
  int mddSize = arg1F->getFullNodeSize(mdd);
  int iMdd = 0;
  int index = 0;
  for (int i = 0; i < mxdSize; i++) {
    index = arg2F->getSparseNodeIndex(mxd, i);
    if (index >= mddSize) break;
    iMdd = arg1F->getFullNodeDownPtr(mdd, index);
    expandMxdHelper(iMdd, arg2F->getSparseNodeDownPtr(mxd, i), result);
  }
  return resF->reduceNode(result);
}


int MEDDLY::postimage_mdd::sparseFull(int mdd, int mxd)
{
  DCASSERT(arg1F->isSparseNode(mdd));
  DCASSERT(arg2F->isFullNode(mxd));

  // create new node
  int result = resF->createTempNodeMaxSize(arg1F->getNodeLevel(mdd));
  int mddSize = arg1F->getSparseNodeSize(mdd);
  int mxdSize = arg2F->getFullNodeSize(mxd);
  int iMxd = 0;
  int index = 0;
  for (int i = 0; i < mddSize; i++) {
    index = arg1F->getSparseNodeIndex(mdd, i);
    if (index >= mxdSize) break;
    iMxd = arg2F->getFullNodeDownPtr(mxd, index);
    expandMxdHelper(arg1F->getSparseNodeDownPtr(mdd, i), iMxd, result);
  }
  return resF->reduceNode(result);
}


int MEDDLY::postimage_mdd::sparseSparse(int mdd, int mxd)
{
  DCASSERT(arg1F->isSparseNode(mdd));
  DCASSERT(arg2F->isSparseNode(mxd));

  // create new node
  int result = resF->createTempNodeMaxSize(arg1F->getNodeLevel(mdd));
  int mddSize = arg1F->getSparseNodeSize(mdd);
  int mxdSize = arg2F->getSparseNodeSize(mxd);
  int i = 0;
  int j = 0;
  int mddIndex = arg1F->getSparseNodeIndex(mdd, i);
  int mxdIndex = arg2F->getSparseNodeIndex(mxd, j);
  for ( ; i < mddSize && j < mxdSize; ) {
    if (mxdIndex < mddIndex) {
      ++j;
      if (j < mxdSize) mxdIndex = arg2F->getSparseNodeIndex(mxd, j);
    }
    else if (mxdIndex > mddIndex) {
      ++i;
      if (i < mddSize) mddIndex = arg1F->getSparseNodeIndex(mdd, i);
    }
    else {
      expandMxdHelper(arg1F->getSparseNodeDownPtr(mdd, i),
          arg2F->getSparseNodeDownPtr(mxd, j), result);
      ++i, ++j;
    }
  }
  return resF->reduceNode(result);
}


int MEDDLY::postimage_mdd::expandMdd(int mdd, int mxd)
{
  // create new node
  int result = resF->createTempNodeMaxSize(arg1F->getNodeLevel(mdd));
  if (arg1F->isFullNode(mdd)) {
    // mdd is a full node
    int mddSize = arg1F->getFullNodeSize(mdd);
    for (int i = 0; i < mddSize; i++) {
      int tempResult = compute(arg1F->getFullNodeDownPtr(mdd, i), mxd);
      resF->setDownPtr(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
  } else {
    // mdd is a sparse node
    int mddSize = arg1F->getSparseNodeSize(mdd);
    for (int i = 0; i < mddSize; i++) {
      int tempResult = compute(arg1F->getSparseNodeDownPtr(mdd, i), mxd);
      resF->setDownPtr(result, arg1F->getSparseNodeIndex(mdd, i), tempResult);
      resF->unlinkNode(tempResult);
    }
  }
  return resF->reduceNode(result);
}


int MEDDLY::postimage_mdd::expandMxd(int mdd, int mxd)
{
  // create new node
  int result = resF->createTempNodeMaxSize(arg2F->getNodeLevel(mxd));
  if (arg2F->isFullNode(mxd)) {
    // full mxd node
    int mxdSize = arg2F->getFullNodeSize(mxd);
    for (int i = 0; i < mxdSize; i++) {
      expandMxdHelper(mdd, arg2F->getFullNodeDownPtr(mxd, i), result);
    } // for mxdSize
  } else {
    // sparse mxd node
    int mxdSize = arg2F->getSparseNodeSize(mxd);
    for (int i = 0; i < mxdSize; i++) {
      expandMxdHelper(mdd, arg2F->getSparseNodeDownPtr(mxd, i), result);
    } // for mxdSize
  }
  return resF->reduceNode(result);
}


void MEDDLY::postimage_mdd::expandMxdHelper(int mdd, int iMxd, int result)
{
  if (iMxd == 0) return;

  DCASSERT(unionOp != 0);
  DCASSERT(!arg2F->isTerminalNode(iMxd));
  DCASSERT(!resF->isReducedNode(result));

  int tempResult = 0;
  if (arg2F->isFullNode(iMxd)) {
    int jResult = 0;
    int ijMxd = 0;
    int iMxdSize = arg2F->getFullNodeSize(iMxd);
    for (int j = 0; j < iMxdSize; j++) {
      ijMxd = arg2F->getFullNodeDownPtr(iMxd, j);
      if (ijMxd == 0) continue;
      tempResult = compute(mdd, ijMxd);
      if (tempResult == 0) continue;
      jResult = unionOp->compute( 
        arg1F->getFullNodeDownPtr(result, j), tempResult
      );
      resF->unlinkNode(tempResult);
      resF->setDownPtr(result, j, jResult);
      resF->unlinkNode(jResult);
    }
  } else {
    // iMxd is a sparse node
    int tempIndex = 0;
    int iMxdSize = arg2F->getSparseNodeSize(iMxd);
    int tempIndexResult = 0;
    for (int j = 0; j < iMxdSize; j++) {
      tempResult = compute(mdd, arg2F->getSparseNodeDownPtr(iMxd, j));
      if (tempResult == 0) continue;
      tempIndex = arg2F->getSparseNodeIndex(iMxd, j);
      tempIndexResult = unionOp->compute(
          arg1F->getFullNodeDownPtr(result, tempIndex), tempResult
      );
      resF->unlinkNode(tempResult);
      resF->setDownPtr(result, tempIndex, tempIndexResult);
      resF->unlinkNode(tempIndexResult);
    }
  }
}


// ******************************************************************
// *                                                                *
// *                     preimage_opname  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::preimage_opname : public binary_opname {
  public:
    preimage_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::preimage_opname::preimage_opname()
 : binary_opname("Pre-image")
{
}

MEDDLY::binary_operation* 
MEDDLY::preimage_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a2->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH);

  return new preimage_mdd(this, a1, a2, r);
}


// ******************************************************************
// *                                                                *
// *                     postimage_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::postimage_opname : public binary_opname {
  public:
    postimage_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::postimage_opname::postimage_opname()
 : binary_opname("Post-image")
{
}

MEDDLY::binary_operation* 
MEDDLY::postimage_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a2->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH);

  return new postimage_mdd(this, a1, a2, r);
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializePreImage(const settings &s)
{
  return new preimage_opname;
}

MEDDLY::binary_opname* MEDDLY::initializePostImage(const settings &s)
{
  return new postimage_opname;
}


/* Pre-image and Post-image operations */

// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************

#if 0

// -------------------------- MDD MXD Image operations ----------------------


class mdd_mxd_image_operation : public mdd_apply_operation {
  public:
    mdd_mxd_image_operation();
    virtual ~mdd_mxd_image_operation();

    virtual void typeCheck(const op_info* owner);
    virtual const char* getName() const { return "Mdd-Mxd Image Operation"; }
    virtual bool isCommutative() const { return false; }

    virtual int compute(op_info* owner, int a, int b) = 0;

  private:
    // Disabling this function
    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mdd_post_image : public mdd_mxd_image_operation {
  public:
    static mdd_post_image* getInstance();
    virtual int compute(op_info* owner, int a, int b);
    virtual const char* getName() const { return "Mdd-Mxd Post Image"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mdd_post_image();
    mdd_post_image(const mdd_post_image& copy);
    mdd_post_image& operator=(const mdd_post_image& copy);
    virtual ~mdd_post_image();

    virtual int compute(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int fullFull(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int fullSparse(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int sparseFull(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int sparseSparse(op_info* owner, op_info* unionOp,
        int mdd, int mxd);
    virtual int expandMdd(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int expandMxd(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual void expandMxdHelper (op_info* owner, op_info* unionOp,
        int mdd, int iMxd, int result);
};


class mdd_pre_image : public mdd_post_image {
  public:
    static mdd_pre_image* getInstance();
    virtual const char* getName() const { return "Mdd-Mxd Pre Image"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mdd_pre_image();
    mdd_pre_image(const mdd_pre_image& copy);
    mdd_pre_image& operator=(const mdd_pre_image& copy);
    virtual ~mdd_pre_image();

    virtual int compute(op_info* owner, op_info* unionOp, int mdd, int mxd);

    virtual int expandMxd (op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int expandMxdByOneLevel(op_info* owner, op_info* unionOp,
        int mdd, int iMxd); 

    virtual int expand(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int expandByOneLevel(op_info* owner, op_info* unionOp,
        int mdd, int iMxd);

    // No need for mdd_pre_image::expandMdd()
    // since we are reusing mdd_post_image::expandMdd()
    // This works correctly because when a mxd level is skipped
    // due to an identity reduced node:
    // result[j] += pre_image(mdd[i], mxd[j][i])
    // is only valid when i == j (those are the only non-zero entries for
    // an identity reduced node).
    // Therefore, the above becomes
    // result[i] += pre_image(mdd[i], mxd[i][i])
    // which is the same expansion used in mdd_post_image::expandMdd().
};


// And these?

class mtmdd_post_image : public mdd_post_image {
  public:
    static mtmdd_post_image* getInstance();
    virtual const char* getName() const { return "MtMdd Post-Image"; }
    virtual void typeCheck(const op_info* owner);
    virtual int compute(op_info* owner, int a, int b);

  protected:
    mtmdd_post_image();
    mtmdd_post_image(const mtmdd_post_image& copy);
    mtmdd_post_image& operator=(const mtmdd_post_image& copy);
    virtual ~mtmdd_post_image();

    virtual int compute(op_info* owner, op_info* plusOp, int mdd, int mxd);
};


class mtmdd_pre_image : public mdd_pre_image {
  public:
    static mtmdd_pre_image* getInstance();
    virtual const char* getName() const { return "MtMdd Pre-Image"; }
    virtual void typeCheck(const op_info* owner);
    virtual int compute(op_info* owner, int a, int b);

  protected:
    mtmdd_pre_image();
    mtmdd_pre_image(const mtmdd_pre_image& copy);
    mtmdd_pre_image& operator=(const mtmdd_pre_image& copy);
    virtual ~mtmdd_pre_image();

    virtual int compute(op_info* owner, op_info* plusOp, int mdd, int mxd);
};


// ---------------------- MDD MXD Image Operation -------------------


mdd_mxd_image_operation::mdd_mxd_image_operation()
{ }

mdd_mxd_image_operation::~mdd_mxd_image_operation() {}


bool mdd_mxd_image_operation::checkTerminals(op_info* op,
  int a, int b, int& c)
{
  assert(false);
  return false;
}


void
mdd_mxd_image_operation::typeCheck(const op_info* owner)
{
  // op1 == MDD, op2 == MXD, op3 = MDD
  if (owner == 0)
    throw error(error::UNKNOWN_OPERATION);
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    throw error(error::TYPE_MISMATCH);
  if (!owner->areAllForests())
    throw error(error::TYPE_MISMATCH);
  if (owner->nParams != 3)
    throw error(error::WRONG_NUMBER);
  if (owner->p[0] == owner->p[1] || owner->p[0] != owner->p[2])
    throw error(error::FOREST_MISMATCH);
  if (!getExpertForest(owner, 0)->isMdd() ||
      !getExpertForest(owner, 1)->isMxd())
    throw error(error::TYPE_MISMATCH);
}


// ------------------------------------------------------------------


// ---------------------- MDD Post-Image Operation -------------------


mdd_post_image* mdd_post_image::getInstance()
{
  static mdd_post_image instance;
  return &instance;
}


mdd_post_image::mdd_post_image()
{ }


mdd_post_image::~mdd_post_image() {}


int mdd_post_image::compute(op_info* owner, int mdd, int mxd)
{
  DCASSERT(owner->nParams == 3 && owner->p[0] == owner->p[2]);
  const int nOperands = 3;
  op_param plist[nOperands] = {owner->p[0], owner->p[0], owner->p[0]};
  op_info* unionOp = 
    smart_cast<expert_compute_manager*>(MEDDLY::getComputeManager())->
    getOpInfo(compute_manager::UNION, plist, nOperands);
  assert(unionOp != 0);

#if 0
  return smart_cast<expert_forest*>(owner->f[1])->getReductionRule() ==
    forest::IDENTITY_REDUCED
    ? compute(owner, unionOp, mdd, mxd)
    : computeNonIdent(owner, unionOp, mdd, mxd);
#else
  return compute(owner, unionOp, mdd, mxd);
#endif
}


// HERE: Test this
// TODO: Test
int mdd_post_image::compute(op_info* owner, op_info* unionOp, int mdd, int mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;

  // skipped level means identity matrix
  if (mxd == -1) {
    getExpertForest(owner, 0)->linkNode(mdd);
    return mdd;
  }

  // check the cache
  int result = 0;
  if (findResult(owner, mdd, mxd, result)) {
#ifdef POST_CACHE_HIT
    printf("%s: cache hit (%d, %d, %d)\n", __func__, mdd, mxd, result);
#endif
    return result;
  }

  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  // check if mxd and mdd are at the same level
  int mddLevel = mddNm->getNodeLevel(mdd);
  int mxdLevel = mxdNm->getNodeLevel(mxd);

  if (mddLevel < mxdLevel) {
    // expand mxd
    result = expandMxd(owner, unionOp, mdd, mxd);
  } else if (mddLevel > mxdLevel) {
    // expand mdd
    result = expandMdd(owner, unionOp, mdd, mxd);
  } else {
    // same level
    DCASSERT(mddNm->getNodeLevel(mdd) == mxdNm->getNodeLevel(mxd));
    if (mddNm->isFullNode(mdd)) {
      if (mxdNm->isFullNode(mxd))
        result = fullFull(owner, unionOp, mdd, mxd);
      else
        result = fullSparse(owner, unionOp, mdd, mxd);
    } else {
      if (mxdNm->isFullNode(mxd))
        result = sparseFull(owner, unionOp, mdd, mxd);
      else
        result = sparseSparse(owner, unionOp, mdd, mxd);
    }
  } // same level

  DCASSERT(mddNm->isReducedNode(result));

  // save result in compute cache and return it
  saveResult(owner, mdd, mxd, result);
  return result;
}


int mdd_post_image::expandMdd (op_info* owner, op_info* unionOp,
  int mdd, int mxd)
{
  expert_forest* mddNm = getExpertForest(owner, 0);

  // create new node
  int result = mddNm->createTempNodeMaxSize(mddNm->getNodeLevel(mdd));
  if (mddNm->isFullNode(mdd)) {
    // mdd is a full node
    int mddSize = mddNm->getFullNodeSize(mdd);
    for (int i = 0; i < mddSize; i++) {
      int tempResult = compute(owner, unionOp,
          mddNm->getFullNodeDownPtr(mdd, i), mxd);
      mddNm->setDownPtr(result, i, tempResult);
      mddNm->unlinkNode(tempResult);
    }
  } else {
    // mdd is a sparse node
    int mddSize = mddNm->getSparseNodeSize(mdd);
    for (int i = 0; i < mddSize; i++) {
      int tempResult = compute(owner, unionOp,
          mddNm->getSparseNodeDownPtr(mdd, i), mxd);
      mddNm->setDownPtr(result, mddNm->getSparseNodeIndex(mdd, i), tempResult);
      mddNm->unlinkNode(tempResult);
    }
  }
  return mddNm->reduceNode(result);
}

void mdd_post_image::expandMxdHelper (op_info* owner, op_info* unionOp,
  int mdd, int iMxd, int result)
{
  if (iMxd == 0) return;

  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);
#if 0
  expert_compute_manager* cm = 
    smart_cast<expert_compute_manager*>(MEDDLY::getComputeManager());
#endif

  mdd_union* unionOpPtr = smart_cast<mdd_union*>(unionOp->op);
  DCASSERT(unionOpPtr != 0);

  DCASSERT(!mxdNm->isTerminalNode(iMxd));
  DCASSERT(!mddNm->isReducedNode(result));

  int tempResult = 0;
  if (mxdNm->isFullNode(iMxd)) {
    int jResult = 0;
    int ijMxd = 0;
    int iMxdSize = mxdNm->getFullNodeSize(iMxd);
#if 0
    dd_edge tempResultEdge(mddNm);
    dd_edge jResultEdge(mddNm);
#endif
    for (int j = 0; j < iMxdSize; j++) {
      ijMxd = mxdNm->getFullNodeDownPtr(iMxd, j);
      if (ijMxd == 0) continue;
      tempResult = compute(owner, unionOp, mdd, ijMxd);
      if (tempResult == 0) continue;
#if 0
      tempResultEdge.set(tempResult, 0, mddNm->getNodeLevel(tempResult));
      jResult = mddNm->getFullNodeDownPtr(result, j);
      mddNm->linkNode(jResult);
      jResultEdge.set(jResult, 0, mddNm->getNodeLevel(jResult));
      cm->apply(unionOp, jResultEdge, tempResultEdge, jResultEdge);
      mddNm->setDownPtr(result, j, jResultEdge.getNode());
#else
      jResult = unionOpPtr->compute(unionOp,
          mddNm->getFullNodeDownPtr(result, j), tempResult);
      mddNm->unlinkNode(tempResult);
      mddNm->setDownPtr(result, j, jResult);
      mddNm->unlinkNode(jResult);
#endif
    }
  } else {
    // iMxd is a sparse node
    int tempIndex = 0;
    int iMxdSize = mxdNm->getSparseNodeSize(iMxd);
    int tempIndexResult = 0;
#if 0
    dd_edge tempResultEdge(mddNm);
    dd_edge tempIndexResultEdge(mddNm);
#endif
    for (int j = 0; j < iMxdSize; j++) {
      tempResult = compute(owner, unionOp, mdd,
          mxdNm->getSparseNodeDownPtr(iMxd, j));
      if (tempResult == 0) continue;
#if 0
      tempResultEdge.set(tempResult, 0, mddNm->getNodeLevel(tempResult));
      tempIndex = mxdNm->getSparseNodeIndex(iMxd, j);
      tempIndexResult = mddNm->getFullNodeDownPtr(result, tempIndex);
      mddNm->linkNode(tempIndexResult);
      tempIndexResultEdge.set(tempIndexResult, 0,
          mddNm->getNodeLevel(tempIndexResult));
      cm->apply(unionOp, tempIndexResultEdge, tempResultEdge,
          tempIndexResultEdge);
      mddNm->setDownPtr(result, tempIndex, tempIndexResultEdge.getNode());
#else
      tempIndex = mxdNm->getSparseNodeIndex(iMxd, j);
      tempIndexResult = unionOpPtr->compute(unionOp,
          mddNm->getFullNodeDownPtr(result, tempIndex), tempResult);
      mddNm->unlinkNode(tempResult);
      mddNm->setDownPtr(result, tempIndex, tempIndexResult);
      mddNm->unlinkNode(tempIndexResult);
#endif
    }
  }
}

int mdd_post_image::expandMxd (op_info* owner, op_info* unionOp,
  int mdd, int mxd)
{
  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  // create new node
  int result = mddNm->createTempNodeMaxSize(mxdNm->getNodeLevel(mxd));
  if (mxdNm->isFullNode(mxd)) {
    // full mxd node
    int mxdSize = mxdNm->getFullNodeSize(mxd);
    for (int i = 0; i < mxdSize; i++) {
      expandMxdHelper(owner, unionOp, mdd,
          mxdNm->getFullNodeDownPtr(mxd, i), result);
    } // for mxdSize
  } else {
    // sparse mxd node
    int mxdSize = mxdNm->getSparseNodeSize(mxd);
    for (int i = 0; i < mxdSize; i++) {
      expandMxdHelper(owner, unionOp, mdd,
          mxdNm->getSparseNodeDownPtr(mxd, i), result);
    } // for mxdSize
  }
  return mddNm->reduceNode(result);
}


int mdd_post_image::fullFull (op_info* owner, op_info* unionOp,
  int mdd, int mxd)
{
  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  DCASSERT(mddNm->isFullNode(mdd));
  DCASSERT(mxdNm->isFullNode(mxd));

  // for each mxd[i]
  //   for each mxd[i][j]
  //     dptrs[j] += postimage(mdd[i],mxd[i][j])

  // create new node
  int minSize = MIN(mddNm->getFullNodeSize(mdd), mxdNm->getFullNodeSize(mxd));
  // need to create result of max size since post image computation may
  // write to result[i] where i >= minSize
  int result = mddNm->createTempNodeMaxSize(mddNm->getNodeLevel(mdd));

  int iMdd = 0;
  int iMxd = 0;
  for (int i = 0; i < minSize; i++) {
    iMdd = mddNm->getFullNodeDownPtr(mdd, i);
    if (iMdd == 0) continue;
    iMxd = mxdNm->getFullNodeDownPtr(mxd, i);
    expandMxdHelper(owner, unionOp, iMdd, iMxd, result);
  }
  return mddNm->reduceNode(result);
}


int mdd_post_image::fullSparse (op_info* owner, op_info* unionOp,
  int mdd, int mxd)
{
  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  DCASSERT(mddNm->isFullNode(mdd));
  DCASSERT(mxdNm->isSparseNode(mxd));

  // create new node
  int result = mddNm->createTempNodeMaxSize(mddNm->getNodeLevel(mdd));
  int mxdSize = mxdNm->getSparseNodeSize(mxd);
  int mddSize = mddNm->getFullNodeSize(mdd);
  int iMdd = 0;
  int index = 0;
  for (int i = 0; i < mxdSize; i++) {
    index = mxdNm->getSparseNodeIndex(mxd, i);
    if (index >= mddSize) break;
    iMdd = mddNm->getFullNodeDownPtr(mdd, index);
    expandMxdHelper(owner, unionOp,
        iMdd, mxdNm->getSparseNodeDownPtr(mxd, i), result);
  }
  return mddNm->reduceNode(result);
}


int mdd_post_image::sparseFull (op_info* owner, op_info* unionOp,
  int mdd, int mxd)
{
  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  DCASSERT(mddNm->isSparseNode(mdd));
  DCASSERT(mxdNm->isFullNode(mxd));

  // create new node
  int result = mddNm->createTempNodeMaxSize(mddNm->getNodeLevel(mdd));
  int mddSize = mddNm->getSparseNodeSize(mdd);
  int mxdSize = mxdNm->getFullNodeSize(mxd);
  int iMxd = 0;
  int index = 0;
  for (int i = 0; i < mddSize; i++) {
    index = mddNm->getSparseNodeIndex(mdd, i);
    if (index >= mxdSize) break;
    iMxd = mxdNm->getFullNodeDownPtr(mxd, index);
    expandMxdHelper(owner, unionOp,
        mddNm->getSparseNodeDownPtr(mdd, i), iMxd, result);
  }
  return mddNm->reduceNode(result);
}


int mdd_post_image::sparseSparse (op_info* owner, op_info* unionOp,
  int mdd, int mxd)
{
  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  DCASSERT(mddNm->isSparseNode(mdd));
  DCASSERT(mxdNm->isSparseNode(mxd));

  // create new node
  int result = mddNm->createTempNodeMaxSize(mddNm->getNodeLevel(mdd));
  int mddSize = mddNm->getSparseNodeSize(mdd);
  int mxdSize = mxdNm->getSparseNodeSize(mxd);
  int i = 0;
  int j = 0;
  int mddIndex = mddNm->getSparseNodeIndex(mdd, i);
  int mxdIndex = mxdNm->getSparseNodeIndex(mxd, j);
  for ( ; i < mddSize && j < mxdSize; ) {
    if (mxdIndex < mddIndex) {
      ++j;
      if (j < mxdSize) mxdIndex = mxdNm->getSparseNodeIndex(mxd, j);
    }
    else if (mxdIndex > mddIndex) {
      ++i;
      if (i < mddSize) mddIndex = mddNm->getSparseNodeIndex(mdd, i);
    }
    else {
      expandMxdHelper(owner, unionOp, mddNm->getSparseNodeDownPtr(mdd, i),
          mxdNm->getSparseNodeDownPtr(mxd, j), result);
      ++i, ++j;
    }
  }
  return mddNm->reduceNode(result);
}


// ------------------------------------------------------------------


// ---------------------- MDD Pre-Image Operation -------------------


mdd_pre_image* mdd_pre_image::getInstance()
{
  static mdd_pre_image instance;
  return &instance;
}


mdd_pre_image::mdd_pre_image()
{ }


mdd_pre_image::~mdd_pre_image() {}


int mdd_pre_image::compute(op_info* owner, op_info* unionOp, int mdd, int mxd)
{
  // result[j] += pre_image(mdd[i], mxd[j][i])

  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;

  // skipped level means identity matrix
  if (mxd == -1) {
    getExpertForest(owner, 0)->linkNode(mdd);
    return mdd;
  }

  // check the cache
  int result = 0;
  if (findResult(owner, mdd, mxd, result)) {
#ifdef POST_CACHE_HIT
    printf("%s: cache hit (%d, %d, %d)\n", __func__, mdd, mxd, result);
#endif
    return result;
  }

  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  // check if mxd and mdd are at the same level
  int mddLevel = mddNm->getNodeLevel(mdd);
  int mxdLevel = mxdNm->getNodeLevel(mxd);

  if (mddLevel < mxdLevel) {
    // expand mxd
    result = expandMxd(owner, unionOp, mdd, mxd);
  } else if (mddLevel > mxdLevel) {
    // expand mdd
    result = expandMdd(owner, unionOp, mdd, mxd);
  } else {
    // same level
    DCASSERT(mddNm->getNodeLevel(mdd) == mxdNm->getNodeLevel(mxd));
    result = expand(owner, unionOp, mdd, mxd);
  } // same level

  DCASSERT(mddNm->isReducedNode(result));

  // save result in compute cache and return it
  saveResult(owner, mdd, mxd, result);
  return result;
}


int mdd_pre_image::expandMxdByOneLevel(op_info* owner, op_info* unionOp,
  int mdd, int iMxd)
{
  // result += pre_image(mdd, iMxd[j])

  if (iMxd == 0) return 0;

  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  mdd_union* unionOpPtr = smart_cast<mdd_union*>(unionOp->op);
  DCASSERT(unionOpPtr != 0);

  DCASSERT(!mxdNm->isTerminalNode(iMxd));

  int result = 0;
  int tempResult = 0;
  if (mxdNm->isFullNode(iMxd)) {
    int iMxdSize = mxdNm->getFullNodeSize(iMxd);
    for (int j = 0; j < iMxdSize; j++) {
      int ijMxd = mxdNm->getFullNodeDownPtr(iMxd, j);
      if (ijMxd == 0) continue;
      int jResult = compute(owner, unionOp, mdd, ijMxd);
      if (jResult == 0) continue;
      tempResult = result;
      result = unionOpPtr->compute(unionOp, result, jResult);
      mddNm->unlinkNode(jResult);
      mddNm->unlinkNode(tempResult);
    }
  } else {
    // iMxd is a sparse node
    int iMxdSize = mxdNm->getSparseNodeSize(iMxd);
    for (int j = 0; j < iMxdSize; j++) {
      int jResult = compute(owner, unionOp, mdd,
          mxdNm->getSparseNodeDownPtr(iMxd, j));
      if (jResult == 0) continue;
      tempResult = result;
      result = unionOpPtr->compute(unionOp, result, jResult);
      mddNm->unlinkNode(jResult);
      mddNm->unlinkNode(tempResult);
    }
  }

  return result;
}

int mdd_pre_image::expandMxd (op_info* owner, op_info* unionOp,
  int mdd, int mxd)
{
  // result[j] += pre_image(mdd[i], mxd[j][i])
  // becomes
  // result[j] += pre_image(mdd, mxd[j][i])
  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  // create new node
  int result = 0;
  if (mxdNm->isFullNode(mxd)) {
    // full mxd node
    int mxdSize = mxdNm->getFullNodeSize(mxd);
    result = mddNm->createTempNode(mxdNm->getNodeLevel(mxd), mxdSize, false);
    for (int i = 0; i < mxdSize; i++) {
      int tempResult = expandMxdByOneLevel(owner, unionOp, mdd,
          mxdNm->getFullNodeDownPtr(mxd, i));
      mddNm->setDownPtrWoUnlink(result, i, tempResult);
      mddNm->unlinkNode(tempResult);
    } // for mxdSize
  } else {
    // sparse mxd node
    int mxdSize = mxdNm->getSparseNodeSize(mxd);
    result = mddNm->createTempNode(mxdNm->getNodeLevel(mxd),
        mxdNm->getSparseNodeIndex(mxd, mxdSize-1) + 1);
    for (int i = 0; i < mxdSize; i++) {
      int tempResult = expandMxdByOneLevel(owner, unionOp, mdd,
          mxdNm->getSparseNodeDownPtr(mxd, i));
      mddNm->setDownPtrWoUnlink(result, mxdNm->getSparseNodeIndex(mxd, i),
          tempResult);
      mddNm->unlinkNode(tempResult);
    } // for mxdSize
  }
  return mddNm->reduceNode(result);
}


int mdd_pre_image::expand(op_info* owner, op_info* unionOp, int mdd, int mxd)
{
  // result[j] += pre_image(mdd[i], mxd[j][i])
  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  // create new node
  int result = 0;
  if (mxdNm->isFullNode(mxd)) {
    // full mxd node
    int mxdSize = mxdNm->getFullNodeSize(mxd);
    result = mddNm->createTempNode(mxdNm->getNodeLevel(mxd), mxdSize, false);
    for (int i = 0; i < mxdSize; i++) {
      int tempResult = expandByOneLevel(owner, unionOp, mdd,
          mxdNm->getFullNodeDownPtr(mxd, i));
      mddNm->setDownPtrWoUnlink(result, i, tempResult);
      mddNm->unlinkNode(tempResult);
    } // for mxdSize
  } else {
    // sparse mxd node
    int mxdSize = mxdNm->getSparseNodeSize(mxd);
    result = mddNm->createTempNode(mxdNm->getNodeLevel(mxd),
        mxdNm->getSparseNodeIndex(mxd, mxdSize - 1) + 1);
    for (int i = 0; i < mxdSize; i++) {
      int tempResult = expandByOneLevel(owner, unionOp, mdd,
          mxdNm->getSparseNodeDownPtr(mxd, i));
      mddNm->setDownPtrWoUnlink(result, mxdNm->getSparseNodeIndex(mxd, i),
          tempResult);
      mddNm->unlinkNode(tempResult);
    } // for mxdSize
  }
  return mddNm->reduceNode(result);
}


int mdd_pre_image::expandByOneLevel(op_info* owner, op_info* unionOp,
  int mdd, int iMxd)
{
  // result += pre_image(mdd[j], iMxd[j])

  if (iMxd == 0) return 0;

  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

  mdd_union* unionOpPtr = smart_cast<mdd_union*>(unionOp->op);
  DCASSERT(unionOpPtr != 0);

  DCASSERT(!mxdNm->isTerminalNode(iMxd));

  int result = 0;
  int tempResult = 0;
  if (mxdNm->isFullNode(iMxd)) {
    if (mddNm->isFullNode(mdd)) {
      int minSize =
        MIN(mddNm->getFullNodeSize(mdd), mxdNm->getFullNodeSize(iMxd));
      for (int j = 0; j < minSize; j++) {
        int ijMxd = mxdNm->getFullNodeDownPtr(iMxd, j);
        int jMdd = mddNm->getFullNodeDownPtr(mdd, j);
        if (ijMxd == 0 || jMdd == 0) continue;
        int jResult = compute(owner, unionOp, jMdd, ijMxd);
        if (jResult == 0) continue;
        tempResult = result;
        result = unionOpPtr->compute(unionOp, result, jResult);
        mddNm->unlinkNode(jResult);
        mddNm->unlinkNode(tempResult);
      }
    } // mdd is a full node
    else {
      DCASSERT(mddNm->isSparseNode(mdd));
      int iMxdSize = mxdNm->getFullNodeSize(iMxd);
      int mddSize = mddNm->getSparseNodeSize(mdd);
      for (int j = 0; j < mddSize; j++) {
        int jIndex = mddNm->getSparseNodeIndex(mdd, j);
        if (jIndex >= iMxdSize) break;
        int ijMxd = mxdNm->getFullNodeDownPtr(iMxd, jIndex);
        if (ijMxd == 0) continue;
        int jMdd = mddNm->getSparseNodeDownPtr(mdd, j);
        int jResult = compute(owner, unionOp, jMdd, ijMxd);
        if (jResult == 0) continue;
        tempResult = result;
        result = unionOpPtr->compute(unionOp, result, jResult);
        mddNm->unlinkNode(jResult);
        mddNm->unlinkNode(tempResult);
      }
    } // mdd is a sparse node
  } // iMxd is a full node
  else {
    DCASSERT(mxdNm->isSparseNode(iMxd));
    if (mddNm->isFullNode(mdd)) {
      int iMxdSize = mxdNm->getSparseNodeSize(iMxd);
      int mddSize = mddNm->getFullNodeSize(mdd);
      for (int j = 0; j < iMxdSize; j++) {
        int jIndex = mxdNm->getSparseNodeIndex(iMxd, j);
        if (jIndex >= mddSize) break;
        int jMdd = mddNm->getFullNodeDownPtr(mdd, jIndex);
        if (jMdd == 0) continue;
        int ijMxd = mxdNm->getSparseNodeDownPtr(iMxd, j);
        int jResult = compute(owner, unionOp, jMdd, ijMxd);
        if (jResult == 0) continue;
        tempResult = result;
        result = unionOpPtr->compute(unionOp, result, jResult);
        mddNm->unlinkNode(jResult);
        mddNm->unlinkNode(tempResult);
      }
    } // mdd is a full node
    else {
      DCASSERT(mddNm->isSparseNode(mdd));
      int iMxdLargestIndex = mxdNm->getSparseNodeIndex(iMxd,
        mxdNm->getSparseNodeSize(iMxd) - 1);
      int mddSize = mddNm->getSparseNodeSize(mdd);
      // k: traverses iMxd, j: traverses mdd
      int j = 0;
      int k = 0;
      int iMxdIndex = mxdNm->getSparseNodeIndex(iMxd, k);
      for ( ; j < mddSize; ) {
        int mddIndex = mddNm->getSparseNodeIndex(mdd, j);
        if (mddIndex > iMxdLargestIndex) break;
        while (iMxdIndex < mddIndex) {
          DCASSERT((k+1) < mxdNm->getSparseNodeSize(iMxd));
          iMxdIndex = mxdNm->getSparseNodeIndex(iMxd, ++k);
        }
        if (mddIndex < iMxdIndex) { ++j; continue; }
        DCASSERT(mddIndex == iMxdIndex);
        int jResult = compute(owner, unionOp,
            mddNm->getSparseNodeDownPtr(mdd, j),
            mxdNm->getSparseNodeDownPtr(iMxd, k));
        ++j; ++k;
        if (jResult == 0) continue;
        tempResult = result;
        result = unionOpPtr->compute(unionOp, result, jResult);
        mddNm->unlinkNode(jResult);
        mddNm->unlinkNode(tempResult);
      }
    } // mdd is a sparse node
  } // iMxd is a sparse node

  return result;
}


// ------------------------------------------------------------------


// ---------------------- MTMDD Post-Image Operation -------------------


mtmdd_post_image* mtmdd_post_image::getInstance()
{
  static mtmdd_post_image instance;
  return &instance;
}


mtmdd_post_image::mtmdd_post_image()
{ }


mtmdd_post_image::~mtmdd_post_image() {}


void
mtmdd_post_image::typeCheck(const op_info* owner)
{
  // op1 == MTMDD, op2 == MTMXD, op3 = MTMDD
  if (owner == 0)
    throw error(error::UNKNOWN_OPERATION);
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    throw error(error::TYPE_MISMATCH);
  if (!owner->areAllForests())
    throw error(error::TYPE_MISMATCH);
  if (owner->nParams != 3)
    throw error(error::WRONG_NUMBER);
  if (owner->p[0] == owner->p[1] || owner->p[0] != owner->p[2])
    throw error(error::FOREST_MISMATCH);
  if (!getExpertForest(owner, 0)->isMtMdd() ||
      !getExpertForest(owner, 1)->isMtMxd())
    throw error(error::TYPE_MISMATCH);
}


int mtmdd_post_image::compute(op_info* owner, int mdd, int mxd)
{
  DCASSERT(owner->nParams == 3 && owner->p[0] == owner->p[2]);
  const int nOperands = 3;
  op_param plist[nOperands] = {owner->p[0], owner->p[0], owner->p[0]};
  op_info* plusOp = 0;
/*
    smart_cast<expert_compute_manager*>(MEDDLY::getComputeManager())->
    getOpInfo(compute_manager::PLUS, plist, nOperands);
*/
  assert(plusOp != 0);
  return compute(owner, plusOp, mdd, mxd);
}


// HERE: Test this
// TODO: Test
int mtmdd_post_image::compute(op_info* owner, op_info* plusOp, int mdd, int mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;

  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

#if 0
  // skipped level means identity matrix
  if (mxd == -1) {
    getExpertForest(owner, 0)->linkNode(mdd);
    return mdd;
  }
#else
  if (mddNm->isTerminalNode(mdd) && mxdNm->isTerminalNode(mxd)) {
    return mxd;
  }
#endif

  // check the cache
  int result = 0;
  if (findResult(owner, mdd, mxd, result)) {
#ifdef POST_CACHE_HIT
    printf("%s: cache hit (%d, %d, %d)\n", __func__, mdd, mxd, result);
#endif
    return result;
  }

  // check if mxd and mdd are at the same level
  int mddLevel = mddNm->getNodeLevel(mdd);
  int mxdLevel = mxdNm->getNodeLevel(mxd);

  if (mddLevel < mxdLevel) {
    // expand mxd
    result = expandMxd(owner, plusOp, mdd, mxd);
  } else if (mddLevel > mxdLevel) {
    // expand mdd
    result = expandMdd(owner, plusOp, mdd, mxd);
  } else {
    // same level
    DCASSERT(mddNm->getNodeLevel(mdd) == mxdNm->getNodeLevel(mxd));
    if (mddNm->isFullNode(mdd)) {
      if (mxdNm->isFullNode(mxd))
        result = fullFull(owner, plusOp, mdd, mxd);
      else
        result = fullSparse(owner, plusOp, mdd, mxd);
    } else {
      if (mxdNm->isFullNode(mxd))
        result = sparseFull(owner, plusOp, mdd, mxd);
      else
        result = sparseSparse(owner, plusOp, mdd, mxd);
    }
  } // same level

  DCASSERT(mddNm->isReducedNode(result));

  // save result in compute cache and return it
  saveResult(owner, mdd, mxd, result);
  return result;
}


// ------------------------------------------------------------------


// ---------------------- MTMDD Pre-Image Operation -------------------


mtmdd_pre_image* mtmdd_pre_image::getInstance()
{
  static mtmdd_pre_image instance;
  return &instance;
}


mtmdd_pre_image::mtmdd_pre_image()
{ }


mtmdd_pre_image::~mtmdd_pre_image() {}


void
mtmdd_pre_image::typeCheck(const op_info* owner)
{
  // op1 == MTMDD, op2 == MTMXD, op3 = MTMDD
  if (owner == 0)
    throw error(error::UNKNOWN_OPERATION);
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    throw error(error::TYPE_MISMATCH);
  if (!owner->areAllForests())
    throw error(error::TYPE_MISMATCH);
  if (owner->nParams != 3)
    throw error(error::WRONG_NUMBER);
  if (owner->p[0] == owner->p[1] || owner->p[0] != owner->p[2])
    throw error(error::FOREST_MISMATCH);
  if (!getExpertForest(owner, 0)->isMtMdd() ||
      !getExpertForest(owner, 1)->isMtMxd())
    throw error(error::TYPE_MISMATCH);
}


int mtmdd_pre_image::compute(op_info* owner, int mdd, int mxd)
{
  DCASSERT(owner->nParams == 3 && owner->p[0] == owner->p[2]);
  const int nOperands = 3;
  op_param plist[nOperands] = {owner->p[0], owner->p[0], owner->p[0]};
  op_info* plusOp = 0;
/*
    smart_cast<expert_compute_manager*>(MEDDLY::getComputeManager())->
    getOpInfo(compute_manager::PLUS, plist, nOperands);
*/
  assert(plusOp != 0);
  return compute(owner, plusOp, mdd, mxd);
}


int mtmdd_pre_image::compute(op_info* owner, op_info* plusOp, int mdd, int mxd)
{
  // result[j] += pre_image(mdd[i], mxd[j][i])

  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;

  expert_forest* mddNm = getExpertForest(owner, 0);
  expert_forest* mxdNm = getExpertForest(owner, 1);

#if 0
  // skipped level means identity matrix
  if (mxd == -1) {
    getExpertForest(owner, 0)->linkNode(mdd);
    return mdd;
  }
#else
  if (mddNm->isTerminalNode(mdd) && mxdNm->isTerminalNode(mxd)) {
    return mxd;
  }
#endif

  // check the cache
  int result = 0;
  if (findResult(owner, mdd, mxd, result)) {
#ifdef POST_CACHE_HIT
    printf("%s: cache hit (%d, %d, %d)\n", __func__, mdd, mxd, result);
#endif
    return result;
  }

  // check if mxd and mdd are at the same level
  int mddLevel = mddNm->getNodeLevel(mdd);
  int mxdLevel = mxdNm->getNodeLevel(mxd);

  if (mddLevel < mxdLevel) {
    // expand mxd
    result = expandMxd(owner, plusOp, mdd, mxd);
  } else if (mddLevel > mxdLevel) {
    // expand mdd
    result = expandMdd(owner, plusOp, mdd, mxd);
  } else {
    // same level
    DCASSERT(mddNm->getNodeLevel(mdd) == mxdNm->getNodeLevel(mxd));
    result = expand(owner, plusOp, mdd, mxd);
  } // same level

  DCASSERT(mddNm->isReducedNode(result));

  // save result in compute cache and return it
  saveResult(owner, mdd, mxd, result);
  return result;
}



// ------------------------------------------------------------------

#endif



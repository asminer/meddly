
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
    virtual int compute(int a, int b);
  protected:
    binary_operation* unionOp;
    virtual int compute_rec(int a, int b) = 0;
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
  int cnode = compute(a.getNode(), b.getNode());
  c.set(cnode, 0, resF->getNodeLevel(cnode));
}

int MEDDLY::image_op::compute(int a, int b)
{
  if (resF->getRangeType() == forest::BOOLEAN) {
    unionOp = getOperation(UNION, resF, resF, resF);
  } else {
    unionOp = getOperation(MAXIMUM, resF, resF, resF);
  }
  DCASSERT(unionOp);
  return compute_rec(a, b);
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

  protected:
    virtual int compute_rec(int a, int b);
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

int MEDDLY::preimage_mdd::compute_rec(int mdd, int mxd)
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
      int jResult = compute_rec(mdd, ijMxd);
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
      int jResult = compute_rec(mdd, arg2F->getSparseNodeDownPtr(iMxd, j));
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
      int tempResult = compute_rec(arg1F->getFullNodeDownPtr(mdd, i), mxd);
      resF->setDownPtr(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
  } else {
    // mdd is a sparse node
    int mddSize = arg1F->getSparseNodeSize(mdd);
    for (int i = 0; i < mddSize; i++) {
      int tempResult = compute_rec(arg1F->getSparseNodeDownPtr(mdd, i), mxd);
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
        int jResult = compute_rec(jMdd, ijMxd);
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
        int jResult = compute_rec(jMdd, ijMxd);
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
        int jResult = compute_rec(jMdd, ijMxd);
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
        iMxdIndex = arg2F->getSparseNodeIndex(iMxd, k);
        while (iMxdIndex < mddIndex) {
          k++;
          DCASSERT(k < arg2F->getSparseNodeSize(iMxd));
          iMxdIndex = arg2F->getSparseNodeIndex(iMxd, k);
        }
        if (mddIndex < iMxdIndex) { ++j; continue; }
        DCASSERT(mddIndex == iMxdIndex);
        int jResult = compute_rec(
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

  protected:
    virtual int compute_rec(int a, int b);
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

int MEDDLY::postimage_mdd::compute_rec(int mdd, int mxd)
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
      int tempResult = compute_rec(arg1F->getFullNodeDownPtr(mdd, i), mxd);
      resF->setDownPtr(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
  } else {
    // mdd is a sparse node
    int mddSize = arg1F->getSparseNodeSize(mdd);
    for (int i = 0; i < mddSize; i++) {
      int tempResult = compute_rec(arg1F->getSparseNodeDownPtr(mdd, i), mxd);
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
      tempResult = compute_rec(mdd, ijMxd);
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
      tempResult = compute_rec(mdd, arg2F->getSparseNodeDownPtr(iMxd, j));
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


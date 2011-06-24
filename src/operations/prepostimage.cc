
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


#include "prepostimage.h"
#include "../compute_cache.h"


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



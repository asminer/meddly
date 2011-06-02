
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


#include "operation_ext.h"
#include "../compute_cache.h"

//#define IGNORE_TERMS 0
//#define IGNORE_INCOUNT 2

mdd_apply_operation::
mdd_apply_operation()
{}

mdd_apply_operation::
~mdd_apply_operation()
{}

compute_manager::error
mdd_apply_operation::
typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0] != owner->p[1] || owner->p[0] != owner->p[2])
    return compute_manager::FOREST_MISMATCH;
  if (!getExpertForest(owner, 0)->isMdd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mdd_apply_operation::
isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[i] <--> forest[i]
  // call isStale for each forest[i] and data[i]
  DCASSERT(owner->nParams == 3);
  return 
    getExpertForest(owner, 0)->isStale(data[0]) ||
    getExpertForest(owner, 1)->isStale(data[1]) ||
    getExpertForest(owner, 2)->isStale(data[2]);
}


void
mdd_apply_operation::
discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[i] <--> forest[i]
  // call uncacheNode for each forest[i] and data[i]
  DCASSERT(owner->nParams == 3);
  getExpertForest(owner, 0)->uncacheNode(data[0]);
  getExpertForest(owner, 1)->uncacheNode(data[1]);
  getExpertForest(owner, 2)->uncacheNode(data[2]);
}


void
mdd_apply_operation::
showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // data[] is of size owner.nParams
  // data[i] <--> forest[i]
  // call showNode for each forest[i] and data[i]
  DCASSERT(owner->nParams == 3);
  fprintf(strm, "[%s(%d, %d): %d]",
      owner->op->getName(), data[0], data[1], data[2]);
#if 0
  getExpertForest(owner, 0)->showNode(strm, data[0]);
  getExpertForest(owner, 1)->showNode(strm, data[1]);
  getExpertForest(owner, 2)->showNode(strm, data[2]);
#endif
}


bool
mdd_apply_operation::
findResult(op_info* owner, int a, int b, int& c)
{
  static int key[2];

  // create cache entry
  if (isCommutative() && a > b) {
    // sort the entry in ascending order
    key[0] = b; key[1] = a;
  } else {
    key[0] = a; key[1] = b;
  }

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  c = cacheEntry[2];
  getExpertForest(owner, 2)->linkNode(c);
  return true;
}


void
mdd_apply_operation::
saveResult(op_info* owner, int a, int b, int c)
{
  static int cacheEntry[3];

  // create cache entry
  if (isCommutative() && a > b) {
    // sort the entry in ascending order
    cacheEntry[0] = b; cacheEntry[1] = a;
  } else {
    cacheEntry[0] = a; cacheEntry[1] = b;
  }
  cacheEntry[2] = c;

  getExpertForest(owner, 0)->cacheNode(cacheEntry[0]);
  getExpertForest(owner, 1)->cacheNode(cacheEntry[1]);
  getExpertForest(owner, 2)->cacheNode(cacheEntry[2]);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}


compute_manager::error
mdd_apply_operation::
compute(op_info* owner, dd_edge** operands)
{
  if (operands == 0) return compute_manager::TYPE_MISMATCH;
  // compute(owner, dd_edge, dd_edge, dd_edge) checks for owner == 0
  return compute(owner, *operands[0], *operands[1], *operands[2]);
}

compute_manager::error
mdd_apply_operation::
compute(op_info* owner, const dd_edge& a, dd_edge& b)
{
  return compute_manager::TYPE_MISMATCH;
}

compute_manager::error
mdd_apply_operation::
compute(op_info* owner, const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  if (owner == 0) return compute_manager::TYPE_MISMATCH;
  int result = compute(owner, a.getNode(), b.getNode());
  c.set(result, 0, getExpertForest(owner, 2)->getNodeLevel(result));
  return compute_manager::SUCCESS;
}


void
mdd_apply_operation::
expandA(op_info* owner, expert_forest* expertForest,
  int a, int b, int result, int resultSize)
{
  // fill in result node and return
  // result[i] = compute(a[i], b)

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    int zeroB = compute(owner, 0, b);
    for (int i = 0; i < aSize; ++i)
    {
      int iA = expertForest->getFullNodeDownPtr(a, i);
      if (iA == 0) {
        expertForest->setDownPtrWoUnlink(result, i, zeroB);
      }
      else {
        int tempResult = compute(owner, iA, b);
        expertForest->setDownPtrWoUnlink(result, i, tempResult);
        expertForest->unlinkNode(tempResult);
      }
    }
    for (int i = aSize; i < resultSize; ++i)
    {
      expertForest->setDownPtrWoUnlink(result, i, zeroB);
    }
    expertForest->unlinkNode(zeroB);
  } else {
    DCASSERT(expertForest->isSparseNode(a));
    const int aSize = expertForest->getSparseNodeSize(a);
    int zeroB = compute(owner, 0, b);
    // j goes through every index (like a full-node index pointer)
    int j = 0;
    for (int i = 0; i < aSize; ++i, ++j)
    {
      // sparse-nodes skip indices which represent downpointer 0
      for (int index = expertForest->getSparseNodeIndex(a, i); j < index; ++j)
      {
        expertForest->setDownPtrWoUnlink(result, j, zeroB);
      }
      // done with skipped indices; deal with the next sparse node index
      DCASSERT(j == expertForest->getSparseNodeIndex(a, i));
      int tempResult = compute(owner,
          expertForest->getSparseNodeDownPtr(a, i), b);
      expertForest->setDownPtrWoUnlink(result, j, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    DCASSERT(j == expertForest->getSparseNodeIndex(a, aSize - 1) + 1);
    for ( ; j < resultSize; ++j)
    {
      expertForest->setDownPtrWoUnlink(result, j, zeroB);
    }
    expertForest->unlinkNode(zeroB);
  }
}


void
mdd_apply_operation::
expandB(op_info* owner, expert_forest* expertForest,
  int a, int b, int result, int resultSize)
{
  // fill in result node and return
  // result[i] = compute(a, b[i])

  if (expertForest->isFullNode(b)) {
    const int bSize = expertForest->getFullNodeSize(b);
    DCASSERT(bSize <= resultSize);
    int aZero = compute(owner, a, 0);
    for (int i = 0; i < bSize; ++i)
    {
      int iB = expertForest->getFullNodeDownPtr(b, i);
      if (iB == 0) {
        expertForest->setDownPtrWoUnlink(result, i, aZero);
      }
      else {
        int tempResult = compute(owner, a, iB);
        expertForest->setDownPtrWoUnlink(result, i, tempResult);
        expertForest->unlinkNode(tempResult);
      }
    }
    for (int i = bSize; i < resultSize; ++i)
    {
      expertForest->setDownPtrWoUnlink(result, i, aZero);
    }
    expertForest->unlinkNode(aZero);
  } else {
    DCASSERT(expertForest->isSparseNode(b));
    const int bSize = expertForest->getSparseNodeSize(b);
    // j goes through every index (like a full-node index pointer)
    int j = 0;
    int aZero = compute(owner, a, 0);
    for (int i = 0; i < bSize; ++i, ++j)
    {
      // sparse-nodes skip indices which represent downpointer 0
      for (int index = expertForest->getSparseNodeIndex(b, i); j < index; ++j)
      {
        expertForest->setDownPtrWoUnlink(result, j, aZero);
      }
      // done with skipped indices; deal with the next sparse node index
      DCASSERT(j == expertForest->getSparseNodeIndex(b, i));
      int tempResult = compute(owner, a,
          expertForest->getSparseNodeDownPtr(b, i));
      expertForest->setDownPtrWoUnlink(result, j, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    DCASSERT(j == expertForest->getSparseNodeIndex(b, bSize - 1) + 1);
    for ( ; j < resultSize; ++j)
    {
      expertForest->setDownPtrWoUnlink(result, j, aZero);
    }
    expertForest->unlinkNode(aZero);
  }
}


// Single-step compute
// result[i] = a[i] op b[i]
int
mdd_apply_operation::
compute(op_info* owner, int a, int b)
{
  int result = 0;
  expert_forest* expertForest = getExpertForest(owner, 0);
  if (checkTerminals(owner, a, b, result))
    return result;
  if (findResult(owner, a, b, result))
    return result;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = expertForest->getLevelSize(resultLevel);

#if 0

  result = expertForest->createTempNode(resultLevel, resultSize, false);

  if (aLevel < resultLevel) {
    // expand b
    // result[i] = a op b[i]
    expandB(owner, expertForest, a, b, result, resultSize);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    expandA(owner, expertForest, a, b, result, resultSize);
  }
  else {
    DCASSERT(aLevel == bLevel);
    // expand a and b
    // result[i] = a[i] op b[i]

    if (expertForest->isFullNode(a)) {
      if (expertForest->isFullNode(b)) {
        fullFull(owner, expertForest, a, b, result, resultSize);
      }
      else {
        DCASSERT(expertForest->isSparseNode(b));
        fullSparse(owner, expertForest, a, b, result, resultSize);
      }
    }
    else {
      DCASSERT(expertForest->isSparseNode(a));
      if (expertForest->isFullNode(b)) {
        sparseFull(owner, expertForest, a, b, result, resultSize);
      }
      else {
        DCASSERT(expertForest->isSparseNode(b));
        sparseSparse(owner, expertForest, a, b, result, resultSize);
      }
    }
  }

#else

  // Three vectors: operands a and b, and result c
  std::vector<int> A(resultSize, 0);
  std::vector<int> B(resultSize, 0);
  std::vector<int> C(resultSize, 0);

  if (aLevel < resultLevel) {
    // expand b
    // result[i] = a op b[i]
    int aZero = compute(owner, a, 0);
    expertForest->getDownPtrs(b, B);
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    for ( ; iterB != B.end(); iterB++, iterC++)
    {
      if (*iterB == 0) {
        *iterC = aZero;
        expertForest->linkNode(aZero);
      } else {
        *iterC = compute(owner, a, *iterB);
      }
    }
    expertForest->unlinkNode(aZero);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    int zeroB = compute(owner, 0, b);

    expertForest->getDownPtrs(a, A);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterC = C.begin();
    for ( ; iterA != A.end(); iterA++, iterC++)
    {
      if (*iterA == 0) {
        *iterC = zeroB;
        expertForest->linkNode(zeroB);
      } else {
        *iterC = compute(owner, *iterA, b);
      }
    }
    expertForest->unlinkNode(zeroB);
  }
  else {
    // expand both a and b
    // result[i] = a[i] op b[i]
    int zeroZero = compute(owner, 0, 0);

    expertForest->getDownPtrs(a, A);
    expertForest->getDownPtrs(b, B);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    for ( ; iterA != A.end(); iterA++, iterB++, iterC++)
    {
      if (*iterA == 0 && *iterB == 0) {
        *iterC = zeroZero;
        expertForest->linkNode(zeroZero);
      } else {
        *iterC = compute(owner, *iterA, *iterB); 
      }
    }

    expertForest->unlinkNode(zeroZero);
  }

  result = expertForest->createTempNode(resultLevel, C);

#endif

  // save result in compute cache and return it

#if 0
  printf("reduce(%d): ", result);
  result = expertForest->reduceNode(result);
  printf("%d  [", result);
  for (unsigned i = 0; i < C.size(); i++ )
  {
    printf("%d ", C[i]);
  }
  printf("]\n");
#else
  result = expertForest->reduceNode(result);
#endif

  saveResult(owner, a, b, result);
  return result;
}


void
mdd_apply_operation::
fullFull (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  DCASSERT(mddNm->isFullNode(a));
  DCASSERT(mddNm->isFullNode(b));
  DCASSERT(!mddNm->isReducedNode(result));

  int aSize = mddNm->getFullNodeSize(a);
  int bSize = mddNm->getFullNodeSize(b);
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;
  for ( ; i < minSize; ++i)
  {
    int tempResult = compute(owner, mddNm->getFullNodeDownPtr(a, i),
        mddNm->getFullNodeDownPtr(b, i));
    mddNm->setDownPtrWoUnlink(result, i, tempResult);
    mddNm->unlinkNode(tempResult);
  }
  for ( ; i < aSize; ++i)
  {
    int tempResult = compute(owner, mddNm->getFullNodeDownPtr(a, i), 0);
    mddNm->setDownPtrWoUnlink(result, i, tempResult);
    mddNm->unlinkNode(tempResult);
  }
  for ( ; i < bSize; ++i)
  {
    int tempResult = compute(owner, 0, mddNm->getFullNodeDownPtr(b, i));
    mddNm->setDownPtrWoUnlink(result, i, tempResult);
    mddNm->unlinkNode(tempResult);
  }
  if (i < resultSize) {
    int zeroZero = compute(owner, 0, 0);
    for ( ; i < resultSize; ++i)
    {
      mddNm->setDownPtrWoUnlink(result, i, zeroZero);
    }
    mddNm->unlinkNode(zeroZero);
  }
}


void
mdd_apply_operation::
sparseSparse (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  DCASSERT(mddNm->isSparseNode(a));
  DCASSERT(mddNm->isSparseNode(b));
  DCASSERT(!mddNm->isReducedNode(result));

  int aNnz = mddNm->getSparseNodeSize(a);
  int aSize = mddNm->getSparseNodeIndex(a, aNnz - 1) + 1;
  int bNnz = mddNm->getSparseNodeSize(b);
  int bSize = mddNm->getSparseNodeIndex(b, bNnz - 1) + 1;
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;  // index for a
  int j = 0;  // index for b
  int aIndex = mddNm->getSparseNodeIndex(a, i);
  int bIndex = mddNm->getSparseNodeIndex(b, j);
  int zeroZero = compute(owner, 0, 0);
  int kA = 0;
  int kB = 0;
  int k = 0;
  for ( ; k < minSize; ++k)
  {
    DCASSERT(k <= aIndex);
    DCASSERT(k <= bIndex);
    if (k < aIndex) {
      kA = 0;
    }
    else {
      kA = mddNm->getSparseNodeDownPtr(a, i);
      ++i;
      if (i < aNnz) { aIndex = mddNm->getSparseNodeIndex(a, i); }
    }
    if (k < bIndex) {
      kB = 0;
    }
    else {
      kB = mddNm->getSparseNodeDownPtr(b, j);
      ++j;
      if (j < bNnz) { bIndex = mddNm->getSparseNodeIndex(b, j); }
    }
    if (kA == 0 && kB == 0) {
      mddNm->setDownPtrWoUnlink(result, k, zeroZero);
    }
    else {
      int tempResult = compute(owner, kA, kB);
      mddNm->setDownPtrWoUnlink(result, k, tempResult);
      mddNm->unlinkNode(tempResult);
    }
  }
  for ( ; k < aSize; ++k, ++i)
  {
    aIndex = mddNm->getSparseNodeIndex(a, i);
    for ( ; k < aIndex; ++k)
    {
      // a has zeroes in these slots
      mddNm->setDownPtrWoUnlink(result, k, zeroZero);
    }
    DCASSERT(k == mddNm->getSparseNodeIndex(a, i));
    int tempResult = compute(owner, mddNm->getSparseNodeDownPtr(a, i), 0);
    mddNm->setDownPtrWoUnlink(result, k, tempResult);
    mddNm->unlinkNode(tempResult);
  }
  for ( ; k < bSize; ++k, ++j)
  {
    bIndex = mddNm->getSparseNodeIndex(b, j);
    for ( ; k < bIndex; ++k)
    {
      // b has zeroes in these slots
      mddNm->setDownPtrWoUnlink(result, k, zeroZero);
    }
    DCASSERT(k == mddNm->getSparseNodeIndex(b, j));
    int tempResult = compute(owner, 0, mddNm->getSparseNodeDownPtr(b, j));
    mddNm->setDownPtrWoUnlink(result, k, tempResult);
    mddNm->unlinkNode(tempResult);
  }
  for ( ; k < resultSize; ++k)
  {
    mddNm->setDownPtrWoUnlink(result, k, zeroZero);
  }
  mddNm->unlinkNode(zeroZero);
}


void
mdd_apply_operation::
fullSparse (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  DCASSERT(mddNm->isFullNode(a));
  DCASSERT(mddNm->isSparseNode(b));
  DCASSERT(!mddNm->isReducedNode(result));

  int aSize = mddNm->getFullNodeSize(a);
  int bNnz = mddNm->getSparseNodeSize(b);
  int bSize = mddNm->getSparseNodeIndex(b, bNnz - 1) + 1;
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;
  int j = 0; // j points to sparse node index
  int bIndex = mddNm->getSparseNodeIndex(b, j);
  for ( ; i < minSize; ++i)
  {
    if (i < bIndex) {
      // b has zeroes in these slots
      int tempResult = compute(owner, mddNm->getFullNodeDownPtr(a, i), 0);
      mddNm->setDownPtrWoUnlink(result, i, tempResult);
      mddNm->unlinkNode(tempResult);
    }
    else {
      DCASSERT(i == mddNm->getSparseNodeIndex(b, j));
      int tempResult = compute(owner, mddNm->getFullNodeDownPtr(a, i),
          mddNm->getSparseNodeDownPtr(b, j));
      mddNm->setDownPtrWoUnlink(result, i, tempResult);
      mddNm->unlinkNode(tempResult);
      ++j;
      if (j < bNnz) { bIndex = mddNm->getSparseNodeIndex(b, j); }
    }
  }

  int zeroZero = compute(owner, 0, 0);
  for ( ; i < aSize; ++i)
  {
    int iA = mddNm->getFullNodeDownPtr(a, i);
    if (iA == 0) {
      mddNm->setDownPtrWoUnlink(result, i, zeroZero);
    }
    else {
      int tempResult = compute(owner, iA, 0);
      mddNm->setDownPtrWoUnlink(result, i, tempResult);
      mddNm->unlinkNode(tempResult);
    }
  }
  for ( ; i < bSize; ++i, ++j)
  {
    for (int index = mddNm->getSparseNodeIndex(b, j); i < index; ++i)
    {
      // b has zeroes in these slots
      mddNm->setDownPtrWoUnlink(result, i, zeroZero);
    }
    DCASSERT(i == mddNm->getSparseNodeIndex(b, j));
    int tempResult = compute(owner, 0, mddNm->getSparseNodeDownPtr(b, j));
    mddNm->setDownPtrWoUnlink(result, i, tempResult);
    mddNm->unlinkNode(tempResult);
  }
  for ( ; i < resultSize; ++i)
  {
    mddNm->setDownPtrWoUnlink(result, i, zeroZero);
  }
  mddNm->unlinkNode(zeroZero);
}


void
mdd_apply_operation::
sparseFull (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  DCASSERT(mddNm->isSparseNode(a));
  DCASSERT(mddNm->isFullNode(b));
  DCASSERT(!mddNm->isReducedNode(result));

  int aNnz = mddNm->getSparseNodeSize(a);
  int aSize = mddNm->getSparseNodeIndex(a, aNnz - 1) + 1;
  int bSize = mddNm->getFullNodeSize(b);
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;
  int j = 0; // j points to sparse node index
  int aIndex = mddNm->getSparseNodeIndex(a, j);
  for ( ; i < minSize; ++i)
  {
    if (i < aIndex) {
      // a has zeroes in these slots
      int tempResult = compute(owner, 0, mddNm->getFullNodeDownPtr(b, i));
      mddNm->setDownPtrWoUnlink(result, i, tempResult);
      mddNm->unlinkNode(tempResult);
    }
    else {
      DCASSERT(i == mddNm->getSparseNodeIndex(a, j));
      int tempResult = compute(owner, mddNm->getSparseNodeDownPtr(a, j),
          mddNm->getFullNodeDownPtr(b, i));
      mddNm->setDownPtrWoUnlink(result, i, tempResult);
      mddNm->unlinkNode(tempResult);
      ++j;
      if (j < aNnz) { aIndex = mddNm->getSparseNodeIndex(a, j); }
    }
  }

  int zeroZero = compute(owner, 0, 0);
  for ( ; i < aSize; ++i, ++j)
  {
    for (int index = mddNm->getSparseNodeIndex(a, j); i < index; ++i)
    {
      // b has zeroes in these slots
      mddNm->setDownPtrWoUnlink(result, i, zeroZero);
    }
    DCASSERT(i < aSize && i == mddNm->getSparseNodeIndex(a, j));
    int tempResult = compute(owner, mddNm->getSparseNodeDownPtr(a, j), 0);
    mddNm->setDownPtrWoUnlink(result, i, tempResult);
    mddNm->unlinkNode(tempResult);
  }
  for ( ; i < bSize; ++i)
  {
    int iB = mddNm->getFullNodeDownPtr(b, i);
    if (iB == 0) {
      mddNm->setDownPtrWoUnlink(result, i, zeroZero);
    }
    else {
      int tempResult = compute(owner, 0, iB);
      mddNm->setDownPtrWoUnlink(result, i, tempResult);
      mddNm->unlinkNode(tempResult);
    }
  }
  for ( ; i < resultSize; ++i)
  {
    mddNm->setDownPtrWoUnlink(result, i, zeroZero);
  }
  mddNm->unlinkNode(zeroZero);
}


// ----------------------- MDD Union --------------------------------


mdd_union*
mdd_union::
getInstance()
{
  static mdd_union instance;
  return &instance;
}


mdd_union::
mdd_union()
{ }


mdd_union::
~mdd_union()
{ }


bool
mdd_union::
checkTerminals(op_info* op, int a, int b, int& c)
{
  if (a == 0 || a == b) {
    c = b;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }
  if (b == 0) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }
  if (a == -1 || b == -1) {
    c = -1;
    return true;
  }

  return false;
}


#if 0
#if 0

// BUG!!

// Single-step compute
// result[i] = a[i] op b[i]

int mdd_union::compute(op_info* owner, int a, int b)
{
  // std::cerr << "1-step mdd_union.\n";
  int result = 0;
  expert_forest* expertForest = smart_cast<expert_forest*>(owner->p[0]);
  if (checkTerminals(owner, a, b, result))
    return result;
  if (findResult(owner, a, b, result))
    return result;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  if (aLevel != bLevel) {
    if (aLevel < bLevel) {
      // expand b
      // result[i] = a op b[i]
      // commutative, use expandA(b, a) instead of expandB(a, b)
      int resultSize = expertForest->getLevelSize(bLevel);
      result = expertForest->createTempNode(bLevel, resultSize, false);
      expandA(owner, expertForest, b, a, result, resultSize);
    }
    else {
      DCASSERT(bLevel < aLevel);
      // expand a
      // result[i] = a[i] op b
      int resultSize = expertForest->getLevelSize(aLevel);
      result = expertForest->createTempNode(aLevel, resultSize, false);
      expandA(owner, expertForest, a, b, result, resultSize);
    }
  }
  else {
    DCASSERT(aLevel == bLevel);
    // expand a and b
    // result[i] = a[i] op b[i]

    if (expertForest->isFullNode(a)) {
      if (expertForest->isFullNode(b)) {
        int resultSize = MAX(
            expertForest->getFullNodeSize(a),
            expertForest->getFullNodeSize(b));
        result = expertForest->createTempNode(aLevel, resultSize, false);
        fullFull(owner, expertForest, a, b, result, resultSize);
      }
      else {
        DCASSERT(expertForest->isSparseNode(b));
        int resultSize = MAX(
            expertForest->getFullNodeSize(a),
            (expertForest->getSparseNodeIndex(b,
              expertForest->getSparseNodeSize(b)-1) + 1));
        result = expertForest->createTempNode(aLevel, resultSize, false);
        fullSparse(owner, expertForest, a, b, result, resultSize);
      }
    }
    else {
      DCASSERT(expertForest->isSparseNode(a));
      if (expertForest->isFullNode(b)) {
        // commutative, use fullSparse(b, a) instead of sparseFull(a, b)
        int resultSize = MAX(
            (expertForest->getSparseNodeIndex(a,
              expertForest->getSparseNodeSize(a)-1) + 1),
            expertForest->getFullNodeSize(b));
        result = expertForest->createTempNode(aLevel, resultSize, false);
        fullSparse(owner, expertForest, b, a, result, resultSize);
      }
      else {
        DCASSERT(expertForest->isSparseNode(b));
        int resultSize = MAX(
            (expertForest->getSparseNodeIndex(a,
              expertForest->getSparseNodeSize(a)-1) + 1),
            (expertForest->getSparseNodeIndex(b,
              expertForest->getSparseNodeSize(b)-1) + 1));
        result = expertForest->createTempNode(aLevel, resultSize, false);
        sparseSparse(owner, expertForest, a, b, result, resultSize);
      }
    }
  }

  // save result in compute cache and return it
  result = expertForest->reduceNode(result);
  saveResult(owner, a, b, result);
  return result;
}


void mdd_union::expandA(op_info* owner, expert_forest* expertForest,
  int a, int b, int result, int resultSize)
{
  // fill in result node and return
  // result[i] = compute(a[i], b)

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    for (int i = 0; i < aSize; ++i)
    {
      int aI = expertForest->getFullNodeDownPtr(a, i);
      if (aI == 0) {
        expertForest->setDownPtrWoUnlink(result, i, b);
        continue;
      }
      int tempResult = compute(owner, aI, b);
      expertForest->setDownPtrWoUnlink(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    for (int i = aSize; i < resultSize; ++i)
    {
      expertForest->setDownPtrWoUnlink(result, i, b);
    }
  } else {
    DCASSERT(expertForest->isSparseNode(a));
    const int aSize = expertForest->getSparseNodeSize(a);
    // j goes through every index (like a full-node index pointer)
    int j = 0;
    for (int i = 0; i < aSize; ++i, ++j)
    {
      // sparse-nodes skip indices which represent downpointer 0
      for (int index = expertForest->getSparseNodeIndex(a, i); j < index; ++j)
      {
        expertForest->setDownPtrWoUnlink(result, j, b);
      }
      // done with skipped indices; deal with the next sparse node index
      DCASSERT(j == expertForest->getSparseNodeIndex(a, i));
      int tempResult = compute(owner,
          expertForest->getSparseNodeDownPtr(a, i), b);
      expertForest->setDownPtrWoUnlink(result, j, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    DCASSERT(j == expertForest->getSparseNodeIndex(a, aSize - 1) + 1);
    for ( ; j < resultSize; ++j)
    {
      expertForest->setDownPtrWoUnlink(result, j, b);
    }
  }
}


void mdd_union::fullFull (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  DCASSERT(mddNm->isFullNode(a));
  DCASSERT(mddNm->isFullNode(b));
  DCASSERT(!mddNm->isReducedNode(result));

  int aSize = mddNm->getFullNodeSize(a);
  int bSize = mddNm->getFullNodeSize(b);
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;
  for ( ; i < minSize; ++i)
  {
    int tempResult = compute(owner, mddNm->getFullNodeDownPtr(a, i),
        mddNm->getFullNodeDownPtr(b, i));
    mddNm->setDownPtrWoUnlink(result, i, tempResult);
    mddNm->unlinkNode(tempResult);
  }
  for ( ; i < aSize; ++i)
  {
    // since union(a[i], 0) == a[i], set downptrs to a[i]
    mddNm->setDownPtrWoUnlink(result, i, mddNm->getFullNodeDownPtr(a, i));
  }
  for ( ; i < bSize; ++i)
  {
    // since union(0, b[i]) == b[i], set downptrs to b[i]
    mddNm->setDownPtrWoUnlink(result, i, mddNm->getFullNodeDownPtr(b, i));
  }
  for ( ; i < resultSize; ++i)
  {
    mddNm->setDownPtrWoUnlink(result, i, 0);
  }
}


void mdd_union::sparseSparse (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  DCASSERT(mddNm->isSparseNode(a));
  DCASSERT(mddNm->isSparseNode(b));
  DCASSERT(!mddNm->isReducedNode(result));

#if 0
  // Algo:
  // i: index to indexptrs of a
  // j: index to indexptrs of b
  // if (indexptr(i) < indexptr(j))
  //   result(indexptr(i)) = downptr(i); i++
  // else if (indexptr(i) > indexptr(j))
  //   result(indexptr(j)) = downptr(j); j++
  // else
  //   result(indexptr(i)) = union(downptr(i), downptr(j)); i++; j++

  int i = 0;
  int j = 0;

  int aNnz = mddNm->getSparseNodeSize(a);
  int bNnz = mddNm->getSparseNodeSize(b);

  int indexI, indexJ;

  while (i < aNnz && j < bNnz) {
    indexI = mddNm->getSparseNodeIndex(a, i);
    indexJ = mddNm->getSparseNodeIndex(b, j);
    if (indexI < indexJ) {
      mddNm->setDownPtr(result, indexI, mddNm->getSparseNodeDownPtr(a, i));
      ++i;
    }
    else if (indexI > indexJ) {
      mddNm->setDownPtr(result, indexJ, mddNm->getSparseNodeDownPtr(b, j));
      ++j;
    }
    else {
      // indexI == indexJ
      int tempResult = compute(owner, mddNm->getSparseNodeDownPtr(a, i),
          mddNm->getSparseNodeDownPtr(b, j));
      mddNm->setDownPtr(result, indexI, tempResult);
      mddNm->unlinkNode(tempResult);
      ++i; ++j;
    }
  }
  while (i < aNnz) {
    DCASSERT(mddNm->getFullNodeDownPtr(result, mddNm->getSparseNodeIndex(a, i))
        == 0);
    mddNm->setDownPtr(result, mddNm->getSparseNodeIndex(a, i),
        mddNm->getSparseNodeDownPtr(a, i));
    ++i;
  }
  while (j < bNnz) {
    DCASSERT(mddNm->getFullNodeDownPtr(result, mddNm->getSparseNodeIndex(b, j))
        == 0);
    mddNm->setDownPtr(result, mddNm->getSparseNodeIndex(b, j),
        mddNm->getSparseNodeDownPtr(b, j));
    ++j;
  }

#else

  // aIndex = 0; bIndex = 0;
  // aIndexIndex = index(a, aIndex);
  // bIndexIndex = index(b, bIndex);
  // For i = 0 to (resultSize-1)
  //   aIndexDown = aIndex == i? down(a, aIndex++): 0
  //   bIndexDown = bIndex == i? down(b, bIndex++): 0
  //   result[i] = union(aIndexDown, bIndexDown)

  int aIndex = 0;
  int bIndex = 0;
  int aIndexIndex = mddNm->getSparseNodeIndex(a, aIndex);
  int bIndexIndex = mddNm->getSparseNodeIndex(b, bIndex);
  int aSize = mddNm->getSparseNodeSize(a);
  int bSize = mddNm->getSparseNodeSize(b);
  int tempResult = 0;
  int i = 0;
  int max = 1 + MAX(mddNm->getSparseNodeIndex(a, aSize-1),
      mddNm->getSparseNodeIndex(b, bSize-1));
  for ( ; i < max; ++i)
  {
    if (i == aIndexIndex) {
      if (i == bIndexIndex) {
        // i == aIndexIndex == bIndexIndex
        // result[i] = union(a[aIndex], b[bIndex])
        tempResult = compute(owner, mddNm->getSparseNodeDownPtr(a, aIndex),
            mddNm->getSparseNodeDownPtr(b, bIndex));
        mddNm->setDownPtrWoUnlink(result, i, tempResult);
        mddNm->unlinkNode(tempResult);
        ++aIndex; ++bIndex;
        if (aIndex < aSize) aIndexIndex = mddNm->getSparseNodeIndex(a, aIndex);
        if (bIndex < bSize) bIndexIndex = mddNm->getSparseNodeIndex(b, bIndex);
      }
      else {
        // i == aIndexIndex && i != bIndexIndex
        // result[i] = union(a[aIndex], 0)
        mddNm->setDownPtrWoUnlink(result, i,
            mddNm->getSparseNodeDownPtr(a, aIndex));
        ++aIndex;
        if (aIndex < aSize) aIndexIndex = mddNm->getSparseNodeIndex(a, aIndex);
      }
    }
    else {
      if (i == bIndexIndex) {
        // i == bIndexIndex && i != aIndexIndex
        // result[i] = union(0, b[bIndex])
        mddNm->setDownPtrWoUnlink(result, i,
            mddNm->getSparseNodeDownPtr(b, bIndex));
        ++bIndex;
        if (bIndex < bSize) bIndexIndex = mddNm->getSparseNodeIndex(b, bIndex);
      }
      else {
        // i != bIndexIndex && i != aIndexIndex
        // result[i] = union(0, 0)
        mddNm->setDownPtrWoUnlink(result, i, 0);
      }
    }
  }
  for ( ; i < resultSize; ++i)
  {
    mddNm->setDownPtrWoUnlink(result, i, 0);
  }


#endif
}


void mdd_union::fullSparse (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  DCASSERT(mddNm->isFullNode(a));
  DCASSERT(mddNm->isSparseNode(b));
  DCASSERT(!mddNm->isReducedNode(result));

  int aSize = mddNm->getFullNodeSize(a);
  int bNnz = mddNm->getSparseNodeSize(b);
  int bSize = mddNm->getSparseNodeIndex(b, bNnz - 1) + 1;
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;
  int j = 0; // j points to sparse node index
  int bIndex = mddNm->getSparseNodeIndex(b, j);
  for ( ; i < minSize; ++i)
  {
    if (i < bIndex) {
      // b has zeroes in these slots
      mddNm->setDownPtrWoUnlink(result, i, mddNm->getFullNodeDownPtr(a, i));
    }
    else {
      DCASSERT(i == mddNm->getSparseNodeIndex(b, j));
      int tempResult = compute(owner, mddNm->getFullNodeDownPtr(a, i),
          mddNm->getSparseNodeDownPtr(b, j));
      mddNm->setDownPtrWoUnlink(result, i, tempResult);
      mddNm->unlinkNode(tempResult);
      ++j;
      if (j < bNnz) { bIndex = mddNm->getSparseNodeIndex(b, j); }
    }
  }
  for ( ; i < aSize; ++i)
  {
    mddNm->setDownPtrWoUnlink(result, i, mddNm->getFullNodeDownPtr(a, i));
  }
  for ( ; j < bNnz; ++j, ++i)
  {
    bIndex = mddNm->getSparseNodeIndex(b, j);
    for ( ; i < bIndex; ++i)
    {
      mddNm->setDownPtrWoUnlink(result, i, 0);
    }
    mddNm->setDownPtrWoUnlink(result, i, mddNm->getSparseNodeDownPtr(b, j));
  }
  for ( ; i < resultSize; ++i)
  {
    mddNm->setDownPtrWoUnlink(result, i, 0);
  }
}

#else

// Alternate single-step compute

int
mdd_union::
compute(op_info* owner, int a, int b)
{
  int result = 0;
  expert_forest* expertForest = getExpertForest(owner, 0);
  if (checkTerminals(owner, a, b, result))
    return result;
  if (findResult(owner, a, b, result))
    return result;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;

  // Three vectors: operands a and b, and result c
  if (aLevel < resultLevel) {
    // expand b
    // result[i] = a op b[i]
    std::vector<int> B;
    expertForest->getDownPtrs(b, B);
    std::vector<int> C(B.size(), 0);

    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();

    for ( ; iterB != B.end(); )
    {
      *iterC++ = compute(owner, a, *iterB++);
    }
    result = expertForest->createTempNode(resultLevel, C);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    std::vector<int> A();
    expertForest->getDownPtrs(a, A);
    std::vector<int> C(A.size(), 0);

    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterC = C.begin();
    for ( ; iterA != A.end(); )
    {
        *iterC++ = compute(owner, *iterA++, b);
    }
    result = expertForest->createTempNode(resultLevel, C);
  }
  else {
    // expand both a and b
    // result[i] = a[i] op b[i]
    std::vector<int> A, B;
    expertForest->getDownPtrs(a, A);
    expertForest->getDownPtrs(b, B);
    int min = MIN(A.size(), B.size());
    std::vector<int> C(MAX(A.size(), B.size()), 0);

    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    for (int i = 0; i < min; i++)
    {
      *iterC++ = compute(owner, *iterA++, *iterB++); 
    }
    while (iterA != A.end())
    {
      expertForest->linkNode(*iterA);
      *iterC++ = *iterA++;
    }
    while (iterB != B.end())
    {
      expertForest->linkNode(*iterB);
      *iterC++ = *iterB++;
    }
    result = expertForest->createTempNode(resultLevel, C);
  }


  // save result in compute cache and return it

#if 0
  printf("reduce(%d): ", result);
  result = expertForest->reduceNode(result);
  printf("%d  [", result);
  for (unsigned i = 0; i < C.size(); i++ )
  {
    printf("%d ", C[i]);
  }
  printf("]\n");
#else
  result = expertForest->reduceNode(result);
#endif

  saveResult(owner, a, b, result);
  return result;
}


#endif
#endif


// ------------------------------------------------------------------


// ----------------------- MDD Intersection --------------------------------


mdd_intersection* mdd_intersection::getInstance()
{
  static mdd_intersection instance;
  return &instance;
}


mdd_intersection::mdd_intersection()
{ }


mdd_intersection::~mdd_intersection() {}


bool
mdd_intersection::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (0 == a || 0 == b) {
    c = 0;
    return true;
  }

  if (-1 == a || a == b) {
    c = b;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (-1 == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  return false;
}


#if 0

#if 0

// Two-step compute
// - result[i] = a[i]
// - result[i] = result[i] union b[i]
// This makes for a simpler but less efficient implementation.

int mdd_intersection::compute(op_info* owner, int a, int b)
{
  // std::cerr << "2-step mdd_intersection.\n";
  int result = 0;
  expert_forest* expertForest = getExpertForest(owner, 0);
  
  // check terminals
  if (checkTerminals(owner, a, b, result))
    return result;

  // check cache for result
  if (findResult(owner, a, b, result))
    return result;
  
  // expand nodes
  // 0. initialize result
  // 1. copy node a to node result
  // 2. do operation between contents of result and node b

  // 0. initialize result
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = expertForest->getLevelSize(resultLevel);
  result = expertForest->createTempNode(resultLevel, resultSize);

  // 1. copy node a to node result
  if (aLevel < resultLevel) {
    // all down pointers of result point to node a
    for (int i = 0; i < resultSize; ++i)
      expertForest->setDownPtr(result, i, a);
  }
  else if (expertForest->isFullNode(a)) {
    // a is a full-node
    const int aSize = expertForest->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    for (int i = 0; i < aSize; ++i)
      expertForest->setDownPtr(result,
          i, expertForest->getFullNodeDownPtr(a, i));
  }
  else {
    // a is a sparse-node
    const int aSize = expertForest->getSparseNodeSize(a);
    for (int i = 0; i < aSize; ++i)
      expertForest->setDownPtr(result,
          expertForest->getSparseNodeIndex(a, i),
          expertForest->getSparseNodeDownPtr(a, i));
  }

  // 2. do operation between contents of result and node b
  if (bLevel < resultLevel) {
    for (int i = 0; i < resultSize; ++i)
    {
      int tempResult =
        compute(owner, expertForest->getFullNodeDownPtr(result, i), b);
      expertForest->setDownPtr(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
  }
  else if (expertForest->isFullNode(b)) {
    // b is a full-node
    const int bSize = expertForest->getFullNodeSize(b);
    DCASSERT(bSize <= resultSize);
    for (int i = 0; i < bSize; ++i)
    {
      int tempResult = compute(owner,
          expertForest->getFullNodeDownPtr(result, i),
          expertForest->getFullNodeDownPtr(b, i));
      expertForest->setDownPtr(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
  }
  else {
    // b is a sparse-node
    const int bSize = expertForest->getSparseNodeSize(b);
    DCASSERT(expertForest->getSparseNodeIndex(b, bSize - 1) <= resultSize);
    for (int i = 0; i < bSize; ++i)
    {
      int index = expertForest->getSparseNodeIndex(b, i);
      int tempResult = compute(owner,
          expertForest->getFullNodeDownPtr(result, index),
          expertForest->getSparseNodeDownPtr(b, i));
      expertForest->setDownPtr(result, index, tempResult);
      expertForest->unlinkNode(tempResult);
    }
  }

  // save result in compute cache and return it
  result = expertForest->reduceNode(result);
  saveResult(owner, a, b, result);
  return result;
}

#else

// Single-step compute
// result[i] = a[i] op b[i]

int mdd_intersection::compute(op_info* owner, int a, int b)
{
  // std::cerr << "1-step mdd_intersection.\n";
  int result = 0;
  expert_forest* expertForest = smart_cast<expert_forest*>(owner->p[0]);
  if (checkTerminals(owner, a, b, result))
    return result;
  if (findResult(owner, a, b, result))
    return result;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  if (aLevel != bLevel) {
    if (aLevel < bLevel) {
      // expand b
      // result[i] = a op b[i]
      // commutative, use expandA(b, a) instead of expandB(a, b)
      int resultSize = expertForest->isFullNode(b)?
        expertForest->getFullNodeSize(b):
        (expertForest->getSparseNodeIndex(b,
          expertForest->getSparseNodeSize(b)-1) + 1);
      result = expertForest->createTempNode(bLevel, resultSize);
      expandA(owner, expertForest, b, a, result, resultSize);
    }
    else {
      DCASSERT(bLevel < aLevel);
      // expand a
      // result[i] = a[i] op b
      int resultSize = expertForest->isFullNode(a)?
        expertForest->getFullNodeSize(a):
        (expertForest->getSparseNodeIndex(a,
          expertForest->getSparseNodeSize(a)-1) + 1);
      result = expertForest->createTempNode(aLevel, resultSize);
      expandA(owner, expertForest, a, b, result, resultSize);
    }
  }
  else {
    DCASSERT(aLevel == bLevel);
    // expand a and b
    // result[i] = a[i] op b[i]

    if (expertForest->isFullNode(a)) {
      if (expertForest->isFullNode(b)) {
        int resultSize = MIN(
            expertForest->getFullNodeSize(a),
            expertForest->getFullNodeSize(b));
        result = expertForest->createTempNode(aLevel, resultSize);
        fullFull(owner, expertForest, a, b, result, resultSize);
      }
      else {
        DCASSERT(expertForest->isSparseNode(b));
        int resultSize = MIN(
            expertForest->getFullNodeSize(a),
            (expertForest->getSparseNodeIndex(b,
              expertForest->getSparseNodeSize(b)-1) + 1));
        result = expertForest->createTempNode(aLevel, resultSize);
        fullSparse(owner, expertForest, a, b, result, resultSize);
      }
    }
    else {
      DCASSERT(expertForest->isSparseNode(a));
      if (expertForest->isFullNode(b)) {
        // commutative, use fullSparse(b, a) instead of sparseFull(a, b)
        int resultSize = MIN(
            (expertForest->getSparseNodeIndex(a,
              expertForest->getSparseNodeSize(a)-1) + 1),
            expertForest->getFullNodeSize(b));
        result = expertForest->createTempNode(aLevel, resultSize);
        fullSparse(owner, expertForest, b, a, result, resultSize);
      }
      else {
        DCASSERT(expertForest->isSparseNode(b));
        int resultSize = MIN(
            (expertForest->getSparseNodeIndex(a,
              expertForest->getSparseNodeSize(a)-1) + 1),
            (expertForest->getSparseNodeIndex(b,
              expertForest->getSparseNodeSize(b)-1) + 1));
        result = expertForest->createTempNode(aLevel, resultSize);
        sparseSparse(owner, expertForest, a, b, result, resultSize);
      }
    }
  }

  // save result in compute cache and return it
  result = expertForest->reduceNode(result);
  saveResult(owner, a, b, result);
  return result;
}

#endif


#if 1

void mdd_intersection::expandA(op_info* owner, expert_forest* expertForest,
  int a, int b, int result, int resultSize)
{
  // fill in result node and return
  // result[i] = compute(a[i], b)

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    for (int i = 0; i < aSize; ++i)
    {
      int iA = expertForest->getFullNodeDownPtr(a, i);
      if (iA != 0) {
        int tempResult = compute(owner, iA, b);
        expertForest->setDownPtr(result, i, tempResult);
        expertForest->unlinkNode(tempResult);
      }
    }
  } else {
    DCASSERT(expertForest->isSparseNode(a));
    const int aSize = expertForest->getSparseNodeSize(a);
    for (int i = 0; i < aSize; ++i)
    {
      int tempResult = compute(owner,
          expertForest->getSparseNodeDownPtr(a, i), b);
      expertForest->setDownPtr(result,
          expertForest->getSparseNodeIndex(a, i), tempResult);
      expertForest->unlinkNode(tempResult);
    }
  }
}


void mdd_intersection::expandB(op_info* owner, expert_forest* expertForest,
  int a, int b, int result, int resultSize)
{
  expandA(owner, expertForest, b, a, result, resultSize);
}


void mdd_intersection::fullFull (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  DCASSERT(mddNm->isFullNode(a));
  DCASSERT(mddNm->isFullNode(b));
  DCASSERT(!mddNm->isReducedNode(result));

  int aSize = mddNm->getFullNodeSize(a);
  int bSize = mddNm->getFullNodeSize(b);
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;
  for ( ; i < minSize; ++i)
  {
    int tempResult = compute(owner, mddNm->getFullNodeDownPtr(a, i),
        mddNm->getFullNodeDownPtr(b, i));
    mddNm->setDownPtr(result, i, tempResult);
    mddNm->unlinkNode(tempResult);
  }
}


void mdd_intersection::sparseSparse (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  DCASSERT(mddNm->isSparseNode(a));
  DCASSERT(mddNm->isSparseNode(b));
  DCASSERT(!mddNm->isReducedNode(result));

  // New Algo:
  // i: index to indexptrs of a
  // j: index to indexptrs of b
  // if (indexptr(i) < indexptr(j))
  //   result(indexptr(i)) = downptr(i); i++
  // else if (indexptr(i) > indexptr(j))
  //   result(indexptr(j)) = downptr(j); j++
  // else
  //   result(indexptr(i)) = union(downptr(i), downptr(j)); i++; j++

  int i = 0;
  int j = 0;
  
  int aNnz = mddNm->getSparseNodeSize(a);
  int bNnz = mddNm->getSparseNodeSize(b);

  int indexI, indexJ;

  while (i < aNnz && j < bNnz) {
    indexI = mddNm->getSparseNodeIndex(a, i);
    indexJ = mddNm->getSparseNodeIndex(b, j);
    if (indexI < indexJ) {
      ++i;
    }
    else if (indexI > indexJ) {
      ++j;
    }
    else {
      // indexI == indexJ
      int tempResult = compute(owner, mddNm->getSparseNodeDownPtr(a, i),
          mddNm->getSparseNodeDownPtr(b, j));
      mddNm->setDownPtr(result, indexI, tempResult);
      mddNm->unlinkNode(tempResult);
      ++i; ++j;
    }
  }
}


void mdd_intersection::fullSparse (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  DCASSERT(mddNm->isFullNode(a));
  DCASSERT(mddNm->isSparseNode(b));
  DCASSERT(!mddNm->isReducedNode(result));

  int aSize = mddNm->getFullNodeSize(a);
  int bNnz = mddNm->getSparseNodeSize(b);
  int bSize = mddNm->getSparseNodeIndex(b, bNnz - 1) + 1;
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;
  int bIndex = mddNm->getSparseNodeIndex(b, i);
  while (bIndex < minSize) {
    int tempResult = compute(owner, mddNm->getFullNodeDownPtr(a, bIndex),
        mddNm->getSparseNodeDownPtr(b, i));
    mddNm->setDownPtr(result, bIndex, tempResult);
    mddNm->unlinkNode(tempResult);

    // advance index
    ++i;
    bIndex = mddNm->getSparseNodeIndex(b, i);
  }
}


void mdd_intersection::sparseFull (op_info* owner, expert_forest* mddNm,
  int a, int b, int result, int resultSize)
{
  fullSparse(owner, mddNm, b, a, result, resultSize);
}

#endif

#endif


// ------------------------------------------------------------------


// ----------------------- MDD Difference --------------------------------


mdd_difference* mdd_difference::getInstance()
{
  static mdd_difference instance;
  return &instance;
}


mdd_difference::mdd_difference()
{ }


mdd_difference::~mdd_difference() {}


bool
mdd_difference::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (a == b || 0 == a || -1 == b) {
    c = 0;
    return true;
  }

  if (0 == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  return false;
}

// ------------------------------------------------------------------


// ------------------------- mxd_apply_operation --------------------


mxd_apply_operation::mxd_apply_operation()
{}


mxd_apply_operation::~mxd_apply_operation() {}


compute_manager::error
mxd_apply_operation::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0] != owner->p[1] || owner->p[0] != owner->p[2])
    return compute_manager::FOREST_MISMATCH;
  if (!owner->p[0].isMxd() ||
      !owner->p[0].isBoolForest() ||
      !owner->p[0].isMT())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


int mxd_apply_operation::compute(op_info* owner, int a, int b) {
  return owner->p[0].readForest()->getReductionRule() ==
      forest::IDENTITY_REDUCED
      ? computeIdent(owner, a, b)
      : computeNonIdent(owner, a, b);
}


int mxd_apply_operation::computeNonIdent(op_info* owner, int a, int b)
{
  expert_forest* expertForest = owner->p[0].getForest();

  DCASSERT(expertForest->getReductionRule() != forest::IDENTITY_REDUCED);

  int result = 0;
  if (checkTerminals(owner, a, b, result))
    return result;
  if (findResult(owner, a, b, result))
    return result;

  // expand nodes
  // 0. initialize result
  // 1. deal with special cases
  // 2. copy node a to node result
  // 3. do operation between contents of result and node b

  // 0. initialize result
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  int resultLevel = ABS(aLevel) == ABS(bLevel) ?
    MAX(aLevel, bLevel) :
    (ABS(aLevel) > ABS(bLevel)? aLevel: bLevel);
  int resultSize = expertForest->getLevelSize(resultLevel);
  // result = expertForest->createTempNode(resultLevel, resultSize);

  // get downpointers for a
  std::vector<int> A;
  if (aLevel != resultLevel) {
    if (int(A.size()) < resultSize) A.resize(resultSize);
    std::fill_n(A.begin(), resultSize, a);
  } else {
    expertForest->getDownPtrs(a, A);
  }

  // get downpointers for b
  std::vector<int> B;
  if (bLevel != resultLevel) {
    if (int(B.size()) < resultSize) B.resize(resultSize);
    std::fill_n(B.begin(), resultSize, b);
  } else {
    expertForest->getDownPtrs(b, B);
  }

  // compute C
  std::vector<int> C(resultSize, 0);
  int min = MIN(A.size(), B.size());
  std::vector<int>::iterator aIter = A.begin();
  std::vector<int>::iterator bIter = B.begin();
  std::vector<int>::iterator cIter = C.begin();

  for (std::vector<int>::iterator end = cIter + min; cIter != end; )
  {
    *cIter++ = computeNonIdent(owner, *aIter++, *bIter++);
  }

  for (std::vector<int>::iterator end = A.end(); aIter != end; )
  {
    DCASSERT(cIter != C.end());
    *cIter++ = computeNonIdent(owner, *aIter++, 0);
  }

  for (std::vector<int>::iterator end = B.end(); bIter != end; )
  {
    DCASSERT(cIter != C.end());
    *cIter++ = computeNonIdent(owner, 0, *bIter++);
  }

  DCASSERT(aIter == A.end() && bIter == B.end());

  if (cIter != C.end()) {
    int zeroZero = computeNonIdent(owner, 0, 0);
    if (zeroZero != 0) std::fill(cIter, C.end(), zeroZero);
    if (!expertForest->isTerminalNode(zeroZero)) {
      unsigned count = C.size() - (cIter - C.begin());
      expertForest->getInCount(zeroZero) += int(count);
    }
    expertForest->unlinkNode(zeroZero);
  }

  result = expertForest->createTempNode(resultLevel, C);
  result = expertForest->reduceNode(result);
  saveResult(owner, a, b, result);
  return result;
}


int mxd_apply_operation::computeIdent(op_info* owner, int a, int b)
{
  // Note:
  // The main difference between mxd::compute and mdd::compute is in the
  // manner in which skipped levels are dealt with.
  //
  // If aLevel < bLevel (or vice versa), MXD reduction rules specify that
  // a and b are both unprimed nodes. Therefore, when expanding b (since
  // it is at a higher level), expand the primed nodes also).

  int result = 0;
  expert_forest* expertForest = owner->p[0].getForest();
  if (checkTerminals(owner, a, b, result))
    return result;
  if (findResult(owner, a, b, result))
    return result;

  // expand nodes
  // 0. initialize result
  // 1. deal with special cases
  // 2. copy node a to node result
  // 3. do operation between contents of result and node b

  // 0. initialize result
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  int resultLevel = ABS(aLevel) == ABS(bLevel) ?
    MAX(aLevel, bLevel) :
    (ABS(aLevel) > ABS(bLevel)? aLevel: bLevel);
  int resultSize = expertForest->getLevelSize(resultLevel);
  result = expertForest->createTempNode(resultLevel, resultSize);

  // 1. deal with special cases
  if (aLevel != resultLevel) {
    // WARNING: this is not the same as bLevel == resultLevel as that will
    // include the case where aLevel == bLevel == resultLevel
    // a is lower; therefore expand B
    if (resultLevel < 0)
      singleExpandB(owner, a, b, expertForest, result, resultLevel, resultSize);
    else
      expandB(owner, a, b, expertForest, result, resultLevel, resultSize);
  }
  else if (bLevel != resultLevel) {
    // WARNING: this is not the same as aLevel == resultLevel as that will
    // include the case where aLevel == bLevel == resultLevel
    // b is lower; therefore expand A
    if (resultLevel < 0)
      singleExpandA(owner, a, b, expertForest, result, resultLevel, resultSize);
    else
      expandA(owner, a, b, expertForest, result, resultLevel, resultSize);
  }
  else {
    // 2. copy node a to node result
    if (expertForest->isFullNode(a)) {
      // a is a full-node
      const int aSize = expertForest->getFullNodeSize(a);
      DCASSERT(aSize <= resultSize);
      for (int i = 0; i < aSize; ++i)
        expertForest->setDownPtr(result,
            i, expertForest->getFullNodeDownPtr(a, i));
    }
    else {
      // a is a sparse-node
      const int aSize = expertForest->getSparseNodeSize(a);
      for (int i = 0; i < aSize; ++i)
        expertForest->setDownPtr(result,
            expertForest->getSparseNodeIndex(a, i),
            expertForest->getSparseNodeDownPtr(a, i));
    }

    // 3. do operation between contents of result and node b
    if (expertForest->isFullNode(b)) {
      // b is a full-node
      const int bSize = expertForest->getFullNodeSize(b);
      DCASSERT(bSize <= resultSize);
      for (int i = 0; i < bSize; ++i)
      {
        int tempResult = computeIdent(owner,
            expertForest->getFullNodeDownPtr(result, i),
            expertForest->getFullNodeDownPtr(b, i));
        expertForest->setDownPtr(result, i, tempResult);
        expertForest->unlinkNode(tempResult);
      }
      for (int i = bSize; i < resultSize; ++i)
      {
        int tempResult =
          computeIdent(owner, expertForest->getFullNodeDownPtr(result, i), 0);
        expertForest->setDownPtr(result, i, tempResult);
        expertForest->unlinkNode(tempResult);
      }
    }
    else {
      // b is a sparse-node
      const int bSize = expertForest->getSparseNodeSize(b);
      DCASSERT(expertForest->getSparseNodeIndex(b, bSize - 1) <= resultSize);
      // j goes through every index (like a full-node index pointer)
      int j = 0;
      for (int i = 0; i < bSize; ++i, ++j)
      {
        // sparse-nodes skip indices which represent downpointer 0
        // call compute of those skipped indices
        for (int index = expertForest->getSparseNodeIndex(b, i);
            j < index; ++j)
        {
          int tempResult = computeIdent(owner,
              expertForest->getFullNodeDownPtr(result, j), 0);
          expertForest->setDownPtr(result, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        // done with skipped indices; deal with the next sparse node index
        DCASSERT(j == expertForest->getSparseNodeIndex(b, i));
        int tempResult = computeIdent(owner,
            expertForest->getFullNodeDownPtr(result, j),
            expertForest->getSparseNodeDownPtr(b, i));
        expertForest->setDownPtr(result, j, tempResult);
        expertForest->unlinkNode(tempResult);
      }
      DCASSERT(j == expertForest->getSparseNodeIndex(b, bSize - 1) + 1);
      for ( ; j < resultSize; ++j)
      {
        int tempResult = computeIdent(owner,
            expertForest->getFullNodeDownPtr(result, j), 0);
        expertForest->setDownPtr(result, j, tempResult);
        expertForest->unlinkNode(tempResult);
      }
    }
  }

  result = expertForest->reduceNode(result);
  saveResult(owner, a, b, result);
  return result;
}


void mxd_apply_operation::singleExpandA(op_info* owner, int a, int b,
  expert_forest* expertForest, int result, int resultLevel, int resultSize)
{
  // result[i] = a[i] op b[i], but b[i] = b,
  // therefore, result[i] = a[i] op b

  // compute each result[i] using above before inserting into result node

  int zeroB = computeIdent(owner, 0, b);

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    int i = 0;
    for ( ; i < aSize; ++i)
    {
      int tempResult =
        computeIdent(owner, expertForest->getFullNodeDownPtr(a, i), b);
      expertForest->setDownPtr(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      expertForest->setDownPtr(result, i, zeroB);
    }
  }
  else {
    DCASSERT(expertForest->isSparseNode(a));
    const int aSize = expertForest->getSparseNodeSize(a);
    int i = 0;
    for (int aIndex = 0; aIndex < aSize; ++aIndex, ++i)
    {
      for (int index = expertForest->getSparseNodeIndex(a, aIndex);
          i < index; ++i)
      {
        DCASSERT(i < resultSize);
        expertForest->setDownPtr(result, i, zeroB);
      }
      DCASSERT(i == expertForest->getSparseNodeIndex(a, aIndex));
      DCASSERT(i < resultSize);
      int tempResult =
        computeIdent(owner, expertForest->getSparseNodeDownPtr(a, aIndex), b);
      expertForest->setDownPtr(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      expertForest->setDownPtr(result, i, zeroB);
    }
  }

  expertForest->unlinkNode(zeroB);
}


void mxd_apply_operation::singleExpandB(op_info* owner, int a, int b,
  expert_forest* expertForest, int result, int resultLevel, int resultSize)
{
  // result[i] = a[i] op b[i], but a[i] = b,
  // therefore, result[i] = a op b[i]

  // compute each result[i] using above before inserting into result node

  int aZero = computeIdent(owner, a, 0);

  if (expertForest->isFullNode(b)) {
    const int bSize = expertForest->getFullNodeSize(b);
    DCASSERT(bSize <= resultSize);
    int i = 0;
    for ( ; i < bSize; ++i)
    {
      int tempResult =
        computeIdent(owner, a, expertForest->getFullNodeDownPtr(b, i));
      expertForest->setDownPtr(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      expertForest->setDownPtr(result, i, aZero);
    }
  }
  else {
    DCASSERT(expertForest->isSparseNode(b));
    const int bSize = expertForest->getSparseNodeSize(b);
    int i = 0;
    for (int bIndex = 0; bIndex < bSize; ++bIndex, ++i)
    {
      for (int index = expertForest->getSparseNodeIndex(b, bIndex);
          i < index; ++i)
      {
        DCASSERT(i < resultSize);
        expertForest->setDownPtr(result, i, aZero);
      }
      DCASSERT(i == expertForest->getSparseNodeIndex(b, bIndex));
      DCASSERT(i < resultSize);
      int tempResult =
        computeIdent(owner, a, expertForest->getSparseNodeDownPtr(b, bIndex));
      expertForest->setDownPtr(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      expertForest->setDownPtr(result, i, aZero);
    }
  }

  expertForest->unlinkNode(aZero);
}


void mxd_apply_operation::expandA(op_info* owner, int a, int b,
  expert_forest* expertForest, int result, int resultLevel, int resultSize)
{
  // result[i][j] = a[i][j] op b[i][j]
  // but b[i][j] = b only when i == j, and 0 otherwise
  // result[i][i] = a[i][i] op b, and
  // result[i][j] = a[i][j] op 0, when i != j

  // compute each result[i] using above before inserting into result node

  // when a[i] == 0, a[i][j] = 0 for all j, and
  // result[i][i] = 0 op b
  // result[i][j] = 0 op 0 for i != j
  int zeroOpZero = computeIdent(owner, 0, 0);
  int zeroOpB = computeIdent(owner, 0, b);

  int pResultSize = expertForest->getLevelSize(-resultLevel);

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    for (int i = 0; i < aSize; ++i)
    {
      int iA = expertForest->getFullNodeDownPtr(a, i);
      int iResult = expertForest->createTempNode(-resultLevel, pResultSize);
      if (iA == 0) {
        for (int j = 0; j < pResultSize; ++j)
          expertForest->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
      }
      else if (expertForest->isFullNode(iA)) {
        const int iASize = expertForest->getFullNodeSize(iA);
        for (int j = 0; j < iASize; ++j)
        {
          int ijA = expertForest->getFullNodeDownPtr(iA, j);
          if (ijA == 0) {
            expertForest->setDownPtr(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(owner, ijA, b): computeIdent(owner, ijA, 0);
            expertForest->setDownPtr(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < pResultSize; ++j)
        {
          expertForest->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        }
      }
      else {
        DCASSERT(expertForest->isSparseNode(iA));
        const int iASize = expertForest->getSparseNodeSize(iA);
        int j = 0;
        for (int k = 0; k < iASize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = expertForest->getSparseNodeIndex(iA, k);
              j < index; ++j)
          {
            expertForest->setDownPtr(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iA, k));
          int ijA = expertForest->getSparseNodeDownPtr(iA, k);
          int tempResult = (i == j)?
            computeIdent(owner, ijA, b): computeIdent(owner, ijA, 0);
          expertForest->setDownPtr(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          expertForest->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        }
      }
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtr(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
    // TODO: can optimize when zeroOpZero == zeroOpB == 0
    for (int i = aSize; i < resultSize; ++i)
    {
      int iResult = expertForest->createTempNode(-resultLevel, pResultSize);
      for (int j = 0; j < pResultSize; ++j)
        expertForest->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtr(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
  }
  else {
    DCASSERT(expertForest->isSparseNode(a));
    const int aSize = expertForest->getSparseNodeSize(a);
    DCASSERT(expertForest->getSparseNodeIndex(a, aSize-1) < resultSize);
    int i = 0;
    for (int aIndex = 0; aIndex < aSize; ++aIndex, ++i)
    {
      for (int index = expertForest->getSparseNodeIndex(a, aIndex);
          i < index; ++i)
      {
        int iResult = expertForest->createTempNode(-resultLevel, pResultSize);
        for (int j = 0; j < pResultSize; ++j)
          expertForest->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        iResult = expertForest->reduceNode(iResult);
        expertForest->setDownPtr(result, i, iResult);
        expertForest->unlinkNode(iResult);
      }

      DCASSERT(i == expertForest->getSparseNodeIndex(a, aIndex));
      int iA = expertForest->getSparseNodeDownPtr(a, aIndex);
      int iResult = expertForest->createTempNode(-resultLevel, pResultSize);
      DCASSERT(iA != 0 && iA != -1);
      if (expertForest->isFullNode(iA)) {
        const int iASize = expertForest->getFullNodeSize(iA);
        for (int j = 0; j < iASize; ++j)
        {
          int ijA = expertForest->getFullNodeDownPtr(iA, j);
          if (ijA == 0) {
            expertForest->setDownPtr(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(owner, ijA, b): computeIdent(owner, ijA, 0);
            expertForest->setDownPtr(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < pResultSize; ++j)
        {
          expertForest->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        }
      }
      else {
        DCASSERT(expertForest->isSparseNode(iA));
        const int iASize = expertForest->getSparseNodeSize(iA);
        int j = 0;
        for (int k = 0; k < iASize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = expertForest->getSparseNodeIndex(iA, k);
              j < index; ++j)
          {
            expertForest->setDownPtr(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iA, k));
          int ijA = expertForest->getSparseNodeDownPtr(iA, k);
          int tempResult = (i == j)?
            computeIdent(owner, ijA, b): computeIdent(owner, ijA, 0);
          expertForest->setDownPtr(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          expertForest->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        }
      }
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtr(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }

    // TODO: can optimize when zeroOpZero == zeroOpB == 0
    DCASSERT(i == expertForest->getSparseNodeIndex(a, aSize - 1) + 1);
    for ( ; i < resultSize; ++i)
    {
      int iResult = expertForest->createTempNode(-resultLevel, pResultSize);
      for (int j = 0; j < pResultSize; ++j)
        expertForest->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtr(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
  }

  expertForest->unlinkNode(zeroOpZero);
  expertForest->unlinkNode(zeroOpB);
}


void mxd_apply_operation::expandB(op_info* owner, int a, int b,
  expert_forest* expertForest, int result, int resultLevel, int resultSize)
{
  // result[i][j] = a[i][j] op b[i][j]
  // but a[i][j] = a only when i == j, and 0 otherwise
  // result[i][i] = a op b[i][i], and
  // result[i][j] = 0 op b[i][j], when i != j

  // compute each result[i] using above before inserting into result node

  // when b[i] == 0, b[i][j] = 0 for all j, and
  // result[i][i] = a op 0
  // result[i][j] = 0 op 0 for i != j
  int zeroOpZero = computeIdent(owner, 0, 0);
  int aOpZero = computeIdent(owner, a, 0);

  int pResultSize = expertForest->getLevelSize(-resultLevel);

  if (expertForest->isFullNode(b)) {
    const int bSize = expertForest->getFullNodeSize(b);
    DCASSERT(bSize <= resultSize);
    for (int i = 0; i < bSize; ++i)
    {
      int iB = expertForest->getFullNodeDownPtr(b, i);
      int iResult = expertForest->createTempNode(-resultLevel, pResultSize);
      if (iB == 0) {
        for (int j = 0; j < pResultSize; ++j)
          expertForest->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
      }
      else if (expertForest->isFullNode(iB)) {
        const int iBSize = expertForest->getFullNodeSize(iB);
        for (int j = 0; j < iBSize; ++j)
        {
          int ijB = expertForest->getFullNodeDownPtr(iB, j);
          if (ijB == 0) {
            expertForest->setDownPtr(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(owner, a, ijB): computeIdent(owner, 0, ijB);
            expertForest->setDownPtr(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < pResultSize; ++j)
        {
          expertForest->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        }
      }
      else {
        DCASSERT(expertForest->isSparseNode(iB));
        const int iBSize = expertForest->getSparseNodeSize(iB);
        int j = 0;
        for (int k = 0; k < iBSize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = expertForest->getSparseNodeIndex(iB, k);
              j < index; ++j)
          {
            expertForest->setDownPtr(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iB, k));
          int ijB = expertForest->getSparseNodeDownPtr(iB, k);
          int tempResult = (i == j)?
            computeIdent(owner, a, ijB): computeIdent(owner, 0, ijB);
          expertForest->setDownPtr(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          expertForest->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        }
      }
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtr(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
    // TODO: can optimize when zeroOpZero == aOpZero == 0
    for (int i = bSize; i < resultSize; ++i)
    {
      int iResult = expertForest->createTempNode(-resultLevel, pResultSize);
      for (int j = 0; j < pResultSize; ++j)
        expertForest->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtr(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
  }
  else {
    DCASSERT(expertForest->isSparseNode(b));
    const int bSize = expertForest->getSparseNodeSize(b);
    int i = 0;
    for (int bIndex = 0; bIndex < bSize; ++bIndex, ++i)
    {
      for (int index = expertForest->getSparseNodeIndex(b, bIndex);
          i < index; ++i)
      {
        int iResult = expertForest->createTempNode(-resultLevel, pResultSize);
        for (int j = 0; j < pResultSize; ++j)
          expertForest->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        iResult = expertForest->reduceNode(iResult);
        expertForest->setDownPtr(result, i, iResult);
        expertForest->unlinkNode(iResult);
      }

      DCASSERT(i == expertForest->getSparseNodeIndex(b, bIndex));
      int iB = expertForest->getSparseNodeDownPtr(b, bIndex);
      int iResult = expertForest->createTempNode(-resultLevel, pResultSize);
      DCASSERT(iB != 0 && iB != -1);
      if (expertForest->isFullNode(iB)) {
        const int iBSize = expertForest->getFullNodeSize(iB);
        for (int j = 0; j < iBSize; ++j)
        {
          int ijB = expertForest->getFullNodeDownPtr(iB, j);
          if (ijB == 0) {
            expertForest->setDownPtr(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(owner, a, ijB): computeIdent(owner, 0, ijB);
            expertForest->setDownPtr(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < pResultSize; ++j)
        {
          expertForest->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        }
      }
      else {
        DCASSERT(expertForest->isSparseNode(iB));
        const int iBSize = expertForest->getSparseNodeSize(iB);
        int j = 0;
        for (int k = 0; k < iBSize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = expertForest->getSparseNodeIndex(iB, k);
              j < index; ++j)
          {
            expertForest->setDownPtr(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iB, k));
          int ijB = expertForest->getSparseNodeDownPtr(iB, k);
          int tempResult = (i == j)?
            computeIdent(owner, a, ijB): computeIdent(owner, 0, ijB);
          expertForest->setDownPtr(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          expertForest->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        }
      }
      iResult = expertForest->reduceNode(iResult);
      DCASSERT(i < resultSize);
      expertForest->setDownPtr(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }

    // TODO: can optimize when zeroOpZero == aOpZero == 0
    DCASSERT(i == expertForest->getSparseNodeIndex(b, bSize - 1) + 1);
    for ( ; i < resultSize; ++i)
    {
      int iResult = expertForest->createTempNode(-resultLevel, pResultSize);
      for (int j = 0; j < pResultSize; ++j)
        expertForest->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtr(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
  }

  expertForest->unlinkNode(zeroOpZero);
  expertForest->unlinkNode(aOpZero);
}


// ------------------------------------------------------------------

// ------------------------- mxd_alt_apply_operation --------------------


mxd_alt_apply_operation::mxd_alt_apply_operation()
{ }


mxd_alt_apply_operation::~mxd_alt_apply_operation() {}


compute_manager::error
mxd_alt_apply_operation::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0] != owner->p[1] || owner->p[0] != owner->p[2])
    return compute_manager::FOREST_MISMATCH;
  if (!(owner->p[0].isMxd()) ||
      !owner->p[0].isBoolForest() ||
      !owner->p[0].isMT())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mxd_alt_apply_operation::isEntryStale(const op_info* owner,
    const int* data)
{
  // data[] is of size owner.nParams + 1
  // data[0] is level
  // data[i] <--> forest[i-1]
  DCASSERT(owner->nParams == 3);
  return 
    getExpertForest(owner, 0)->isStale(data[1]) ||
    getExpertForest(owner, 1)->isStale(data[2]) ||
    getExpertForest(owner, 2)->isStale(data[3]);
}


void
mxd_alt_apply_operation::discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams + 1
  // data[0] is level
  // data[i] <--> forest[i-1]
  DCASSERT(owner->nParams == 3);
  getExpertForest(owner, 0)->uncacheNode(data[1]);
  getExpertForest(owner, 1)->uncacheNode(data[2]);
  getExpertForest(owner, 2)->uncacheNode(data[3]);
}


void
mxd_alt_apply_operation::showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // data[] is of size owner.nParams + 1
  // data[0] is level
  // data[i] <--> forest[i-1]
  DCASSERT(owner->nParams == 3);
  fprintf(strm, "[%s %d %d %d %d]",
      owner->op->getName(), data[0], data[1], data[2], data[3]);
}


bool
mxd_alt_apply_operation::findResult(op_info* owner,
    int k, int a, int b, int& c)
{
  static int key[3];

#ifdef IGNORE_TERMS
  if (getExpertForest(owner, 0)->isTerminalNode(a) ||
      getExpertForest(owner, 1)->isTerminalNode(b))
    return false;
#endif

  // create cache entry
  key[0] = k;
  if (isCommutative() && a > b) {
    // sort the entry in ascending order
    key[1] = b; key[2] = a;
  } else {
    key[1] = a; key[2] = b;
  }

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  c = cacheEntry[3];
  getExpertForest(owner, 2)->linkNode(c);
  return true;
}


void
mxd_alt_apply_operation::saveResult(op_info* owner, int k, int a, int b, int c)
{
  static int cacheEntry[4];

#ifdef IGNORE_TERMS
  if (getExpertForest(owner, 0)->isTerminalNode(a) ||
      getExpertForest(owner, 1)->isTerminalNode(b))
    return;
#endif

#ifdef IGNORE_INCOUNT
  if ((!getExpertForest(owner, 0)->isTerminalNode(a) &&
        getExpertForest(owner, 0)->getInCount(a) < IGNORE_INCOUNT) &&
      (!getExpertForest(owner, 0)->isTerminalNode(b) &&
       getExpertForest(owner, 0)->getInCount(b) < IGNORE_INCOUNT))
    return;
#endif

  // create cache entry
  cacheEntry[0] = k;
  if (isCommutative() && a > b) {
    // sort the entry in ascending order
    cacheEntry[1] = b; cacheEntry[2] = a;
  } else {
    cacheEntry[1] = a; cacheEntry[2] = b;
  }
  cacheEntry[3] = c;

#if 0
#ifdef DEVELOPMENT_CODE
  assert(!findResult(owner, k, a, b, c));
#endif
#endif

  getExpertForest(owner, 0)->cacheNode(cacheEntry[1]);
  getExpertForest(owner, 1)->cacheNode(cacheEntry[2]);
  getExpertForest(owner, 2)->cacheNode(cacheEntry[3]);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}


compute_manager::error
mxd_alt_apply_operation::compute(op_info* owner, dd_edge** operands)
{
  if (operands == 0) return compute_manager::TYPE_MISMATCH;
  // compute(owner, dd_edge, dd_edge, dd_edge) checks for owner == 0
  return compute(owner, *operands[0], *operands[1], *operands[2]);
}

compute_manager::error
mxd_alt_apply_operation::compute(op_info* owner, const dd_edge& a, dd_edge& b)
{
  return compute_manager::TYPE_MISMATCH;
}

compute_manager::error
mxd_alt_apply_operation::compute(op_info* owner, const dd_edge& a,
    const dd_edge& b, dd_edge& c)
{
  if (owner == 0) return compute_manager::TYPE_MISMATCH;
  int result = compute(owner,
      owner->p[0].getDomain()->getTopVariable(), a.getNode(), b.getNode());
  c.set(result, 0, getExpertForest(owner, 2)->getNodeLevel(result));
  return compute_manager::SUCCESS;
}


int mxd_alt_apply_operation::compute(op_info* owner, int resultLevel,
    int a, int b)
{
  return owner->p[0].readForest()->getReductionRule() ==
      forest::IDENTITY_REDUCED
      ? computeIdent(owner, resultLevel, a, b)
      : computeNonIdent(owner, resultLevel, a, b);
}


int mxd_alt_apply_operation::computeNonIdent(op_info* owner,
    int resultLevel, int a, int b)
{
  expert_forest* expertForest = owner->p[0].getForest();

  DCASSERT(expertForest->getReductionRule() != forest::IDENTITY_REDUCED);

  int result = 0;

  if (resultLevel == 0) {
    checkTerminals(owner, a, b, result);
    return result;
  }

  if (findResult(owner, resultLevel, a, b, result)) {
    return result;
  }


  // expand nodes
  // 0. initialize result
  // 1. deal with special cases
  // 2. copy node a to node result
  // 3. do operation between contents of result and node b

  // 0. initialize result
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  int resultSize = expertForest->getLevelSize(resultLevel);
  int nextLevel = (resultLevel > 0)? -resultLevel: -resultLevel-1;

  if (aLevel != resultLevel && bLevel != resultLevel) {

    // c[i] = a op b
    int c = computeNonIdent(owner, nextLevel, a, b);
    std::vector<int> C(resultSize, c);
    if (!expertForest->isTerminalNode(c))
      expertForest->getInCount(c) += resultSize;
    expertForest->unlinkNode(c);
    result = expertForest->createTempNode(resultLevel, C);
    result = expertForest->reduceNode(result);

  }
  else {

    // get downpointers for a
    std::vector<int> A;
    if (aLevel != resultLevel) {
      if (int(A.size()) < resultSize) A.resize(resultSize);
      std::fill_n(A.begin(), resultSize, a);
    } else {
      expertForest->getDownPtrs(a, A);
    }

    // get downpointers for b
    std::vector<int> B;
    if (bLevel != resultLevel) {
      if (int(B.size()) < resultSize) B.resize(resultSize);
      std::fill_n(B.begin(), resultSize, b);
    } else {
      expertForest->getDownPtrs(b, B);
    }

    // compute C
    std::vector<int> C(resultSize, 0);

    int min = MIN(A.size(), B.size());
    std::vector<int>::iterator aIter = A.begin();
    std::vector<int>::iterator bIter = B.begin();
    std::vector<int>::iterator cIter = C.begin();

    for (std::vector<int>::iterator end = cIter + min; cIter != end; )
    {
      *cIter++ = computeNonIdent(owner, nextLevel, *aIter++, *bIter++);
    }

    for (std::vector<int>::iterator end = A.end(); aIter != end; )
    {
      DCASSERT(cIter != C.end());
      *cIter++ = computeNonIdent(owner, nextLevel, *aIter++, 0);
    }

    for (std::vector<int>::iterator end = B.end(); bIter != end; )
    {
      DCASSERT(cIter != C.end());
      *cIter++ = computeNonIdent(owner, nextLevel, 0, *bIter++);
    }

    DCASSERT(aIter == A.end() && bIter == B.end());

    if (cIter != C.end()) {
      int zeroZero = computeNonIdent(owner, nextLevel, 0, 0);
      if (zeroZero != 0) std::fill(cIter, C.end(), zeroZero);
      if (!expertForest->isTerminalNode(zeroZero)) {
        unsigned count = C.size() - (cIter - C.begin());
        expertForest->getInCount(zeroZero) += int(count);
      }
      expertForest->unlinkNode(zeroZero);
    }

    result = expertForest->createTempNode(resultLevel, C);
    result = expertForest->reduceNode(result);

  }

  saveResult(owner, resultLevel, a, b, result);
  return result;
}


int mxd_alt_apply_operation::computeIdent(op_info* owner, int resultLevel,
    int a, int b)
{
  // Note 1:
  // The main difference between mxd::compute and mdd::compute is in the
  // manner in which skipped levels are dealt with.
  //
  // If aLevel < bLevel (or vice versa), MXD reduction rules specify that
  // a and b are both unprimed nodes. Therefore, when expanding b (since
  // it is at a higher level), expand the primed nodes also).
  //
  // Note 2:
  // Nodes in mxd_alt_apply_operation::compute() are built by expanding
  // one level at a time (i.e. no levels are skipped when building the result).
  // This is to ensure correctness for operations where op(0, 0) != 0.

  int result = 0;
  expert_forest* expertForest = owner->p[0].getForest();

  DCASSERT(expertForest->getReductionRule() == forest::IDENTITY_REDUCED);

  if (resultLevel == 0) {
    checkTerminals(owner, a, b, result);
    return result;
  }

  if (findResult(owner, resultLevel, a, b, result)) {
    return result;
  }

  // expand nodes
  // 0. initialize result
  // 1. deal with special cases
  // 2. copy node a to node result
  // 3. do operation between contents of result and node b

  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  // 0. initialize result
  int resultSize = expertForest->getLevelSize(resultLevel);
  result = expertForest->createTempNode(resultLevel, resultSize, false);

  if (aLevel != resultLevel && bLevel != resultLevel) {
    expandSkippedLevel(owner, a, b,
        expertForest, result, resultLevel, resultSize);
  }
  else if (aLevel != resultLevel) {
    // aLevel != resultLevel, bLevel == resultLevel
    if (resultLevel < 0)
      singleExpandB(owner, a, b, expertForest, result, resultLevel, resultSize);
    else
      expandB(owner, a, b, expertForest, result, resultLevel, resultSize);
  }
  else if (bLevel != resultLevel) {
    // aLevel == resultLevel, bLevel != resultLevel
    if (resultLevel < 0)
      singleExpandA(owner, a, b, expertForest, result, resultLevel, resultSize);
    else
      expandA(owner, a, b, expertForest, result, resultLevel, resultSize);
  }
  else {
    // aLevel == resultLevel, bLevel == resultLevel
    DCASSERT(aLevel == resultLevel && bLevel == resultLevel);
    // 2. copy node a to node result
    if (expertForest->isFullNode(a)) {
      // a is a full-node
      const int aSize = expertForest->getFullNodeSize(a);
      DCASSERT(aSize <= resultSize);
      for (int i = 0; i < aSize; ++i)
        expertForest->setDownPtrWoUnlink(result,
            i, expertForest->getFullNodeDownPtr(a, i));
      for (int i = aSize; i < resultSize; ++i)
        expertForest->setDownPtrWoUnlink(result, i, 0);
    }
    else {
      // a is a sparse-node
      const int aSize = expertForest->getSparseNodeSize(a);
      for (int i = 0; i < resultSize; ++i)
        expertForest->setDownPtrWoUnlink(result, i, 0);
      for (int i = 0; i < aSize; ++i)
        expertForest->setDownPtrWoUnlink(result,
            expertForest->getSparseNodeIndex(a, i),
            expertForest->getSparseNodeDownPtr(a, i));
    }

    // 3. do operation between contents of result and node b
    int nextLevel = (resultLevel > 0)? -resultLevel: -resultLevel-1;
    if (expertForest->isFullNode(b)) {
      // b is a full-node
      const int bSize = expertForest->getFullNodeSize(b);
      DCASSERT(bSize <= resultSize);
      for (int i = 0; i < bSize; ++i)
      {
        int tempResult = computeIdent(owner, nextLevel,
            expertForest->getFullNodeDownPtr(result, i),
            expertForest->getFullNodeDownPtr(b, i));
        expertForest->setDownPtr(result, i, tempResult);
        expertForest->unlinkNode(tempResult);
      }
      for (int i = bSize; i < resultSize; ++i)
      {
        int tempResult = computeIdent(owner, nextLevel,
            expertForest->getFullNodeDownPtr(result, i), 0);
        expertForest->setDownPtr(result, i, tempResult);
        expertForest->unlinkNode(tempResult);
      }
    }
    else {
      // b is a sparse-node
      const int bSize = expertForest->getSparseNodeSize(b);
      DCASSERT(expertForest->getSparseNodeIndex(b, bSize - 1) <= resultSize);
      // j goes through every index (like a full-node index pointer)
      int j = 0;
      for (int i = 0; i < bSize; ++i, ++j)
      {
        // sparse-nodes skip indices which represent downpointer 0
        // call compute of those skipped indices
        for (int index = expertForest->getSparseNodeIndex(b, i);
            j < index; ++j)
        {
          int tempResult = computeIdent(owner, nextLevel,
              expertForest->getFullNodeDownPtr(result, j), 0);
          expertForest->setDownPtr(result, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        // done with skipped indices; deal with the next sparse node index
        DCASSERT(j == expertForest->getSparseNodeIndex(b, i));
        int tempResult = computeIdent(owner, nextLevel,
            expertForest->getFullNodeDownPtr(result, j),
            expertForest->getSparseNodeDownPtr(b, i));
        expertForest->setDownPtr(result, j, tempResult);
        expertForest->unlinkNode(tempResult);
      }
      DCASSERT(j == expertForest->getSparseNodeIndex(b, bSize - 1) + 1);
      for ( ; j < resultSize; ++j)
      {
        int tempResult = computeIdent(owner, nextLevel,
            expertForest->getFullNodeDownPtr(result, j), 0);
        expertForest->setDownPtr(result, j, tempResult);
        expertForest->unlinkNode(tempResult);
      }
    }
  }

  result = expertForest->reduceNode(result);
  saveResult(owner, resultLevel, a, b, result);
  return result;
}


void mxd_alt_apply_operation::expandSkippedLevel(op_info* owner,
    int a, int b, expert_forest* expertForest,
    int result, int resultLevel, int resultSize)
{
  DCASSERT(expertForest->getNodeLevel(a) != resultLevel);
  DCASSERT(expertForest->getNodeLevel(b) != resultLevel);
  if (resultLevel > 0) {
    DCASSERT(expertForest->getNodeLevel(a) >= 0);
    DCASSERT(expertForest->getNodeLevel(b) >= 0);
    // both a and b are below result
    int zeroZeroAtOneLevelBelow = computeIdent(owner, resultLevel-1, 0, 0);
    int aBAtOneLevelBelow = computeIdent(owner, resultLevel-1, a, b);
    int pResultSize = expertForest->getLevelSize(-resultLevel);

    for (int i = 0; i < resultSize; ++i)
    {
      // create primed node at -resultLevel
      int p = expertForest->createTempNode(-resultLevel, pResultSize, false);
      for (int j = 0; j < pResultSize; ++j)
      {
        // create unprimed node at resultLevel-1
        expertForest->setDownPtrWoUnlink(p, j,
            (i == j)? aBAtOneLevelBelow: zeroZeroAtOneLevelBelow);
      }
      p = expertForest->reduceNode(p);
      expertForest->setDownPtrWoUnlink(result, i, p);
      expertForest->unlinkNode(p);
    }
    expertForest->unlinkNode(zeroZeroAtOneLevelBelow);
    expertForest->unlinkNode(aBAtOneLevelBelow);
  }
  else {
    DCASSERT(resultLevel < 0);
    DCASSERT(a == 0 && b == 0);
    int zeroZeroAtOneLevelBelow = computeIdent(owner, -resultLevel-1, 0, 0);
    for (int i = 0; i < resultSize; ++i)
    {
      expertForest->setDownPtrWoUnlink(result, i, zeroZeroAtOneLevelBelow);
    }
    expertForest->unlinkNode(zeroZeroAtOneLevelBelow);
  }
}


void mxd_alt_apply_operation::singleExpandA(op_info* owner, int a, int b,
  expert_forest* expertForest, int result, int resultLevel, int resultSize)
{
  // result[i] = a[i] op b[i], but b[i] = b,
  // therefore, result[i] = a[i] op b

  // compute each result[i] using above before inserting into result node

  DCASSERT(resultLevel < 0);
  int zeroB = computeIdent(owner, -resultLevel-1, 0, b);

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    int i = 0;
    for ( ; i < aSize; ++i)
    {
      int tempResult = computeIdent(owner, -resultLevel-1,
          expertForest->getFullNodeDownPtr(a, i), b);
      expertForest->setDownPtrWoUnlink(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      expertForest->setDownPtrWoUnlink(result, i, zeroB);
    }
  }
  else {
    DCASSERT(expertForest->isSparseNode(a));
    const int aSize = expertForest->getSparseNodeSize(a);
    int i = 0;
    for (int aIndex = 0; aIndex < aSize; ++aIndex, ++i)
    {
      for (int index = expertForest->getSparseNodeIndex(a, aIndex);
          i < index; ++i)
      {
        expertForest->setDownPtrWoUnlink(result, i, zeroB);
      }
      DCASSERT(i == expertForest->getSparseNodeIndex(a, aIndex));
      int tempResult = computeIdent(owner, -resultLevel-1,
          expertForest->getSparseNodeDownPtr(a, aIndex), b);
      expertForest->setDownPtrWoUnlink(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      expertForest->setDownPtrWoUnlink(result, i, zeroB);
    }
  }

  expertForest->unlinkNode(zeroB);
}


void mxd_alt_apply_operation::singleExpandB(op_info* owner, int a, int b,
  expert_forest* expertForest, int result, int resultLevel, int resultSize)
{
  // result[i] = a[i] op b[i], but a[i] = b,
  // therefore, result[i] = a op b[i]

  // compute each result[i] using above before inserting into result node

  DCASSERT(resultLevel < 0);
  int aZero = computeIdent(owner, -resultLevel-1, a, 0);

  if (expertForest->isFullNode(b)) {
    const int bSize = expertForest->getFullNodeSize(b);
    int i = 0;
    for ( ; i < bSize; ++i)
    {
      int tempResult = computeIdent(owner, -resultLevel-1,
          a, expertForest->getFullNodeDownPtr(b, i));
      expertForest->setDownPtrWoUnlink(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      expertForest->setDownPtrWoUnlink(result, i, aZero);
    }
  }
  else {
    DCASSERT(expertForest->isSparseNode(b));
    const int bSize = expertForest->getSparseNodeSize(b);
    int i = 0;
    for (int bIndex = 0; bIndex < bSize; ++bIndex, ++i)
    {
      for (int index = expertForest->getSparseNodeIndex(b, bIndex);
          i < index; ++i)
      {
        expertForest->setDownPtrWoUnlink(result, i, aZero);
      }
      DCASSERT(i == expertForest->getSparseNodeIndex(b, bIndex));
      int tempResult = computeIdent(owner, -resultLevel-1,
          a, expertForest->getSparseNodeDownPtr(b, bIndex));
      expertForest->setDownPtrWoUnlink(result, i, tempResult);
      expertForest->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      expertForest->setDownPtrWoUnlink(result, i, aZero);
    }
  }

  expertForest->unlinkNode(aZero);
}


void mxd_alt_apply_operation::expandA(op_info* owner, int a, int b,
  expert_forest* expertForest, int result, int resultLevel, int resultSize)
{
  // result[i][j] = a[i][j] op b[i][j]
  // but b[i][j] = b only when i == j, and 0 otherwise
  // result[i][i] = a[i][i] op b, and
  // result[i][j] = a[i][j] op 0, when i != j

  // compute each result[i] using above before inserting into result node

  // when a[i] == 0, a[i][j] = 0 for all j, and
  // result[i][i] = 0 op b
  // result[i][j] = 0 op 0 for i != j
  DCASSERT(resultLevel > 0);
  int zeroOpZero = computeIdent(owner, resultLevel-1, 0, 0);
  int zeroOpB = computeIdent(owner, resultLevel-1, 0, b);

  int pResultSize = expertForest->getLevelSize(-resultLevel);

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    for (int i = 0; i < aSize; ++i)
    {
      int iA = expertForest->getFullNodeDownPtr(a, i);
      int iResult =
        expertForest->createTempNode(-resultLevel, pResultSize, false);
      if (iA == 0) {
        for (int j = 0; j < pResultSize; ++j)
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
      }
      else if (expertForest->isFullNode(iA)) {
        const int iASize = expertForest->getFullNodeSize(iA);
        for (int j = 0; j < iASize; ++j)
        {
          int ijA = expertForest->getFullNodeDownPtr(iA, j);
          if (ijA == 0) {
            expertForest->setDownPtrWoUnlink(iResult, j, i == j
                ? zeroOpB
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(owner, resultLevel-1, ijA, b):
              computeIdent(owner, resultLevel-1, ijA, 0);
            expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < pResultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
        }
      }
      else {
        DCASSERT(expertForest->isSparseNode(iA));
        const int iASize = expertForest->getSparseNodeSize(iA);
        int j = 0;
        for (int k = 0; k < iASize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = expertForest->getSparseNodeIndex(iA, k);
              j < index; ++j)
          {
            expertForest->setDownPtrWoUnlink(iResult, j, i == j
                ? zeroOpB
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iA, k));
          int ijA = expertForest->getSparseNodeDownPtr(iA, k);
          int tempResult = (i == j)?
            computeIdent(owner, resultLevel-1, ijA, b):
            computeIdent(owner, resultLevel-1, ijA, 0);
          expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
        }
      }
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtrWoUnlink(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
    // TODO: can optimize when zeroOpZero == zeroOpB == 0
    for (int i = aSize; i < resultSize; ++i)
    {
      int iResult =
        expertForest->createTempNode(-resultLevel, pResultSize, false);
      for (int j = 0; j < pResultSize; ++j)
        expertForest->setDownPtrWoUnlink(iResult, j,
            i == j? zeroOpB: zeroOpZero);
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtrWoUnlink(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
  }
  else {
    DCASSERT(expertForest->isSparseNode(a));
    const int aSize = expertForest->getSparseNodeSize(a);
    int i = 0;
    for (int aIndex = 0; aIndex < aSize; ++aIndex, ++i)
    {
      for (int index = expertForest->getSparseNodeIndex(a, aIndex);
          i < index; ++i)
      {
        int iResult =
          expertForest->createTempNode(-resultLevel, pResultSize, false);
        for (int j = 0; j < pResultSize; ++j)
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
        iResult = expertForest->reduceNode(iResult);
        expertForest->setDownPtrWoUnlink(result, i, iResult);
        expertForest->unlinkNode(iResult);
      }

      DCASSERT(i == expertForest->getSparseNodeIndex(a, aIndex));
      int iA = expertForest->getSparseNodeDownPtr(a, aIndex);
      int iResult =
        expertForest->createTempNode(-resultLevel, pResultSize, false);
      DCASSERT(iA != 0 && iA != -1);
      if (expertForest->isFullNode(iA)) {
        const int iASize = expertForest->getFullNodeSize(iA);
        for (int j = 0; j < iASize; ++j)
        {
          int ijA = expertForest->getFullNodeDownPtr(iA, j);
          if (ijA == 0) {
            expertForest->setDownPtrWoUnlink(iResult, j,
                i == j? zeroOpB : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(owner, resultLevel-1, ijA, b):
              computeIdent(owner, resultLevel-1, ijA, 0);
            expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < pResultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
        }
      }
      else {
        DCASSERT(expertForest->isSparseNode(iA));
        const int iASize = expertForest->getSparseNodeSize(iA);
        int j = 0;
        for (int k = 0; k < iASize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = expertForest->getSparseNodeIndex(iA, k);
              j < index; ++j)
          {
            expertForest->setDownPtrWoUnlink(iResult, j, i == j
                ? zeroOpB
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iA, k));
          int ijA = expertForest->getSparseNodeDownPtr(iA, k);
          int tempResult = i == j?
            computeIdent(owner, resultLevel-1, ijA, b):
            computeIdent(owner, resultLevel-1, ijA, 0);
          expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
        }
      }
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtrWoUnlink(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }

    // TODO: can optimize when zeroOpZero == zeroOpB == 0
    DCASSERT(i == expertForest->getSparseNodeIndex(a, aSize - 1) + 1);
    for ( ; i < resultSize; ++i)
    {
      int iResult =
        expertForest->createTempNode(-resultLevel, pResultSize, false);
      for (int j = 0; j < pResultSize; ++j)
        expertForest->setDownPtrWoUnlink(iResult, j,
            i == j? zeroOpB: zeroOpZero);
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtrWoUnlink(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
  }

  expertForest->unlinkNode(zeroOpZero);
  expertForest->unlinkNode(zeroOpB);
}


void mxd_alt_apply_operation::expandB(op_info* owner, int a, int b,
  expert_forest* expertForest, int result, int resultLevel, int resultSize)
{
  // result[i][j] = a[i][j] op b[i][j]
  // but a[i][j] = a only when i == j, and 0 otherwise
  // result[i][i] = a op b[i][i], and
  // result[i][j] = 0 op b[i][j], when i != j

  // compute each result[i] using above before inserting into result node

  // when b[i] == 0, b[i][j] = 0 for all j, and
  // result[i][i] = a op 0
  // result[i][j] = 0 op 0 for i != j
  DCASSERT(resultLevel > 0);
  int zeroOpZero = computeIdent(owner, resultLevel-1, 0, 0);
  int aOpZero = computeIdent(owner, resultLevel-1, a, 0);

  int pResultSize = expertForest->getLevelSize(-resultLevel);

  if (expertForest->isFullNode(b)) {
    const int bSize = expertForest->getFullNodeSize(b);
    for (int i = 0; i < bSize; ++i)
    {
      int iB = expertForest->getFullNodeDownPtr(b, i);
      int iResult =
        expertForest->createTempNode(-resultLevel, pResultSize, false);
      if (iB == 0) {
        for (int j = 0; j < pResultSize; ++j)
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
      }
      else if (expertForest->isFullNode(iB)) {
        const int iBSize = expertForest->getFullNodeSize(iB);
        for (int j = 0; j < iBSize; ++j)
        {
          int ijB = expertForest->getFullNodeDownPtr(iB, j);
          if (ijB == 0) {
            expertForest->setDownPtrWoUnlink(iResult, j, i == j
                ? aOpZero
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(owner, resultLevel-1, a, ijB):
              computeIdent(owner, resultLevel-1, 0, ijB);
            expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < pResultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
        }
      }
      else {
        DCASSERT(expertForest->isSparseNode(iB));
        const int iBSize = expertForest->getSparseNodeSize(iB);
        int j = 0;
        for (int k = 0; k < iBSize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = expertForest->getSparseNodeIndex(iB, k);
              j < index; ++j)
          {
            expertForest->setDownPtrWoUnlink(iResult, j, i == j
                ? aOpZero
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iB, k));
          int ijB = expertForest->getSparseNodeDownPtr(iB, k);
          int tempResult = (i == j)?
            computeIdent(owner, resultLevel-1, a, ijB):
            computeIdent(owner, resultLevel-1, 0, ijB);
          expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
        }
      }
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtrWoUnlink(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
    // TODO: can optimize when zeroOpZero == aOpZero == 0
    for (int i = bSize; i < resultSize; ++i)
    {
      int iResult =
        expertForest->createTempNode(-resultLevel, pResultSize, false);
      for (int j = 0; j < pResultSize; ++j)
        expertForest->setDownPtrWoUnlink(iResult, j,
            i == j? aOpZero: zeroOpZero);
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtrWoUnlink(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
  }
  else {
    DCASSERT(expertForest->isSparseNode(b));
    const int bSize = expertForest->getSparseNodeSize(b);
    int i = 0;
    for (int bIndex = 0; bIndex < bSize; ++bIndex, ++i)
    {
      for (int index = expertForest->getSparseNodeIndex(b, bIndex);
          i < index; ++i)
      {
        int iResult =
          expertForest->createTempNode(-resultLevel, pResultSize, false);
        for (int j = 0; j < pResultSize; ++j)
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
        iResult = expertForest->reduceNode(iResult);
        expertForest->setDownPtrWoUnlink(result, i, iResult);
        expertForest->unlinkNode(iResult);
      }

      DCASSERT(i == expertForest->getSparseNodeIndex(b, bIndex));
      int iB = expertForest->getSparseNodeDownPtr(b, bIndex);
      int iResult =
        expertForest->createTempNode(-resultLevel, pResultSize, false);
      DCASSERT(iB != 0 && iB != -1);
      if (expertForest->isFullNode(iB)) {
        const int iBSize = expertForest->getFullNodeSize(iB);
        for (int j = 0; j < iBSize; ++j)
        {
          int ijB = expertForest->getFullNodeDownPtr(iB, j);
          if (ijB == 0) {
            expertForest->setDownPtrWoUnlink(iResult, j, i == j
                ? aOpZero
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(owner, resultLevel-1, a, ijB):
              computeIdent(owner, resultLevel-1, 0, ijB);
            expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < pResultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
        }
      }
      else {
        DCASSERT(expertForest->isSparseNode(iB));
        const int iBSize = expertForest->getSparseNodeSize(iB);
        int j = 0;
        for (int k = 0; k < iBSize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = expertForest->getSparseNodeIndex(iB, k);
              j < index; ++j)
          {
            expertForest->setDownPtrWoUnlink(iResult, j, i == j
                ? aOpZero
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iB, k));
          int ijB = expertForest->getSparseNodeDownPtr(iB, k);
          int tempResult = i == j?
            computeIdent(owner, resultLevel-1, a, ijB):
            computeIdent(owner, resultLevel-1, 0, ijB);
          expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
        }
      }
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtrWoUnlink(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }

    // TODO: can optimize when zeroOpZero == aOpZero == 0
    DCASSERT(i == expertForest->getSparseNodeIndex(b, bSize - 1) + 1);
    for ( ; i < resultSize; ++i)
    {
      int iResult =
        expertForest->createTempNode(-resultLevel, pResultSize, false);
      for (int j = 0; j < pResultSize; ++j)
        expertForest->setDownPtrWoUnlink(iResult, j,
            i == j? aOpZero: zeroOpZero);
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtrWoUnlink(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }
  }

  expertForest->unlinkNode(zeroOpZero);
  expertForest->unlinkNode(aOpZero);
}


// ------------------------------------------------------------------

// ----------------------- MXD Union --------------------------------


mxd_union* mxd_union::getInstance()
{
  static mxd_union instance;
  return &instance;
}


mxd_union::mxd_union()
{ }


mxd_union::~mxd_union() {}


bool
mxd_union::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (0 == a) {
    c = b;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (0 == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (a == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  return false;
}


int mxd_union::computeIdentExpandOneLevel(op_info* owner, int a, int b)
{
  expert_forest* expertForest = owner->p[0].getForest();
  DCASSERT(!expertForest->isTerminalNode(a));
  DCASSERT(!expertForest->isTerminalNode(b));
  DCASSERT(expertForest->getNodeLevel(b) == expertForest->getNodeLevel(a));

  int resultLevel = expertForest->getNodeLevel(a);
  int resultSize = expertForest->getLevelSize(resultLevel);

  // Set-up function pointer for call to the next level.
  // If a and b are at the unprime level, call computeIdentExpandOneLevel().
  // Otherwise, call computeIdent().
  int (mxd_union::*function)(op_info*, int, int) =
    resultLevel > 0
    ? &mxd_union::computeIdentExpandOneLevel
    : &mxd_union::computeIdent;

  // Copy a into result.
  int result = expertForest->makeACopy(a, resultSize);
  int* resultDptrs = 0;
  assert(expertForest->getDownPtrs(result, resultDptrs));

  bool noChange = true;

  // Add b to result.
  int* rDptrs = resultDptrs;
  const int* bDptrs = 0;
  assert(expertForest->getDownPtrs(b, bDptrs));
  if (expertForest->isFullNode(b)) {
    int bSize = expertForest->getFullNodeSize(b);
    for (const int* bEnd = bDptrs + bSize; bDptrs != bEnd; ++rDptrs, ++bDptrs)
    {
      // Terminal conditions.
      if (*rDptrs == *bDptrs || 0 == *bDptrs) continue;
      if (0 == *rDptrs) {
        *rDptrs = *bDptrs;
        expertForest->linkNode(*rDptrs);
        noChange = false;
        continue;
      }
      // Expand *rDptrs and *bDptrs.
      int pNode = (this->*function)(owner, *rDptrs, *bDptrs);
      if (*rDptrs != pNode) {
        expertForest->unlinkNode(*rDptrs);
        *rDptrs = pNode;
        noChange = false;
        continue;
      }
      expertForest->unlinkNode(pNode);
    }
  }
  else {
    DCASSERT(expertForest->isSparseNode(b));
    int nDptrs = expertForest->getSparseNodeSize(b);
    const int* bIndexes = 0;
    assert(expertForest->getSparseNodeIndexes(b, bIndexes));
    for (const int* bEnd = bDptrs + nDptrs; bDptrs != bEnd; ++bDptrs)
    {
      // Terminal conditions
      rDptrs = resultDptrs + *bIndexes++;
      DCASSERT(*bDptrs != 0);
      if (*rDptrs == *bDptrs) continue;
      if (0 == *rDptrs) {
        *rDptrs = *bDptrs;
        expertForest->linkNode(*rDptrs);
        noChange = false;
        continue;
      }
      // Expand *rDptrs and *bDptrs.
      int pNode = (this->*function)(owner, *rDptrs, *bDptrs);
      if (*rDptrs != pNode) {
        expertForest->unlinkNode(*rDptrs);
        *rDptrs = pNode;
        noChange = false;
        continue;
      }
      expertForest->unlinkNode(pNode);
    }

  }

  if (noChange) {
    // result is the same as node a.
    // Don't call reduce; discard result and return a.
    expertForest->linkNode(a);
    expertForest->unlinkNode(result);
    result = a;
  } else {
    result = expertForest->reduceNode(result);
  }

  return result;
}


#if 0
int mxd_union::computeIdentExpand(op_info* owner, int a, int b)
{
  expert_forest* expertForest = owner->p[0].getForest();
  int resultLevel = expertForest->getNodeLevel(a);
  int resultSize = expertForest->getLevelSize(resultLevel);

  // Copy a into result.
  int result = expertForest->makeACopy(a, resultSize);
  int* resultDptrs = 0;
  assert(expertForest->getDownPtrs(result, resultDptrs));

  bool noChange = true;

  // Add b to result.
  int* rDptrs = resultDptrs;
  const int* bDptrs = 0;
  assert(expertForest->getDownPtrs(b, bDptrs));
  if (expertForest->isFullNode(b)) {
    int bSize = expertForest->getFullNodeSize(b);
    for (const int* bEnd = bDptrs + bSize; bDptrs != bEnd; ++rDptrs, ++bDptrs)
    {
      // Terminal conditions.
      if (*rDptrs == *bDptrs || 0 == *bDptrs) continue;
      if (0 == *rDptrs) {
        *rDptrs = *bDptrs;
        expertForest->linkNode(*rDptrs);
        noChange = false;
        continue;
      }
      // Expand *rDptrs and *bDptrs.
      int pNode = computeIdentExpandPrimedLevel(owner, *rDptrs, *bDptrs);
      if (*rDptrs != pNode) {
        expertForest->unlinkNode(*rDptrs);
        *rDptrs = pNode;
        noChange = false;
      }
      expertForest->unlinkNode(pNode);
    }
  }
  else {
    DCASSERT(expertForest->isSparseNode(b));
    int nDptrs = expertForest->getSparseNodeSize(b);
    const int* bIndexes = 0;
    assert(expertForest->getSparseNodeIndexes(b, bIndexes));
    for (const int* bEnd = bDptrs + nDptrs; bDptrs != bEnd; ++bDptrs)
    {
      // Terminal conditions
      rDptrs = resultDptrs + *bIndexes++;
      DCASSERT(*bDptrs != 0);
      if (*rDptrs == *bDptrs) continue;
      if (0 == *rDptrs) {
        *rDptrs = *bDptrs;
        expertForest->linkNode(*rDptrs);
        noChange = false;
        continue;
      }
      // Expand *rDptrs and *bDptrs.
      int pNode = computeIdentExpandPrimedLevel(owner, *rDptrs, *bDptrs);
      if (*rDptrs != pNode) {
        expertForest->unlinkNode(*rDptrs);
        *rDptrs = pNode;
        noChange = false;
      }
      expertForest->unlinkNode(pNode);
    }

  }

  if (noChange) {
    // result is the same as node a.
    // Don't call reduce; discard result and return a.
    expertForest->linkNode(a);
    expertForest->unlinkNode(result);
    result = a;
  } else {
    result = expertForest->reduceNode(result);
  }

  return result;
}
#endif


int mxd_union::computeIdentExpandA(op_info* owner, int a, int b)
{
  expert_forest* expertForest = owner->p[0].getForest();
  int resultLevel = expertForest->getNodeLevel(a);
  int resultSize = expertForest->getLevelSize(resultLevel);

  // Copy a into result.
  int result = expertForest->makeACopy(a, resultSize);
  int* resultDptrs = 0;
  assert(expertForest->getDownPtrs(result, resultDptrs));

  bool noChange = true;

  // Add b to result[i][i]
  int* rDptrs = resultDptrs;
  for (int i = 0; i < resultSize; ++i, ++rDptrs) {
    if (*rDptrs == 0) {
      int pNode = expertForest->createTempNode(-resultLevel, i + 1, true);
      expertForest->setDownPtrWoUnlink(pNode, i, b);
      *rDptrs = expertForest->reduceNode(pNode);
      noChange = false;
    } else {
      int mxdII = expertForest->getDownPtr(*rDptrs, i);
      int temp = 0;
      if (mxdII == 0) {
        temp = b; expertForest->linkNode(b);
      } else {
        temp = computeIdent(owner, mxdII, b);
      }
      if (temp != mxdII) {
        int pNode = expertForest->makeACopy(*rDptrs, i + 1);
        expertForest->setDownPtr(pNode, i, temp);
        expertForest->unlinkNode(*rDptrs);
        *rDptrs = expertForest->reduceNode(pNode);
        noChange = false;
      }
      expertForest->unlinkNode(temp);
    }
  }

  if (noChange) {
    // result is the same as node a.
    // Don't call reduce; discard result and return a.
    expertForest->linkNode(a);
    expertForest->unlinkNode(result);
    result = a;
  } else {
    result = expertForest->reduceNode(result);
  }

  return result;
}


int mxd_union::computeIdent(op_info* owner, int a, int b)
{
  expert_forest* expertForest = owner->p[0].getForest();

  // Terminal conditions for recursion
  if (a == b) {
    expertForest->linkNode(a);
    return a;
  }
  if (0 == a || 0 == b) {
    int c = a + b;
    expertForest->linkNode(c);
    return c;
  }

  // Search compute table
  static int cacheEntry[3];
  if (a > b)  { cacheEntry[0] = b; cacheEntry[1] = a; }
  else        { cacheEntry[0] = a; cacheEntry[1] = b; }
  const int* ans = owner->cc->find(owner, (int*) cacheEntry);
  if (ans != 0) {
    expertForest->linkNode(ans[2]);
    return ans[2];
  }

  int aHeight = expertForest->getMappedNodeHeight(a);
  int bHeight = expertForest->getMappedNodeHeight(b);
  int result = 0;

  if (aHeight > bHeight) {
    // result[i][j] = a[i][j]
    // except result[i][i] = a[i][i] + b
    result = computeIdentExpandA(owner, a, b);
  } else if (aHeight < bHeight) {
    result = computeIdentExpandA(owner, b, a);
  } else {
    result = computeIdentExpandOneLevel(owner, b, a);
  }

  result = expertForest->reduceNode(result);

  // Save result to compute table
  if (a > b)  { cacheEntry[0] = b; cacheEntry[1] = a; }
  else        { cacheEntry[0] = a; cacheEntry[1] = b; }
  cacheEntry[2] = result;
  expertForest->cacheNode(a);
  expertForest->cacheNode(b);
  expertForest->cacheNode(result);
  owner->cc->add(owner, (int*) cacheEntry);

  return result;
}


// ------------------------------------------------------------------


// ----------------------- MXD Intersection --------------------------------


mxd_intersection* mxd_intersection::getInstance()
{
  static mxd_intersection instance;
  return &instance;
}


mxd_intersection::mxd_intersection()
{ }


mxd_intersection::~mxd_intersection() {}


bool
mxd_intersection::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (0 == a || 0 == b) {
    c = 0;
    return true;
  }

  // this comparison also covers a == b == -1
  if (a == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  return false;
}

// ------------------------------------------------------------------


// ----------------------- MXD Difference --------------------------------


mxd_difference* mxd_difference::getInstance()
{
  static mxd_difference instance;
  return &instance;
}


mxd_difference::mxd_difference()
{ }


mxd_difference::~mxd_difference() {}


bool
mxd_difference::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (a == b || 0 == a) {
    c = 0;
    return true;
  }

  if (0 == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  return false;
}

// ------------------------------------------------------------------


// ----------------------- MXD Complement --------------------------------


mxd_complement* mxd_complement::getInstance()
{
  static mxd_complement instance;
  return &instance;
}


mxd_complement::mxd_complement()
{ }


mxd_complement::~mxd_complement() {}


compute_manager::error
mxd_complement::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 0)->isMxd())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 1)->isMxd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mxd_complement::isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[0] --> level
  // data[1,2] <--> forest[0,1]
  // call isStale for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  return 
    getExpertForest(owner, 0)->isStale(data[1]) ||
    getExpertForest(owner, 1)->isStale(data[2]);
}


void
mxd_complement::discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[0] --> level
  // data[1,2] <--> forest[0,1]
  // call uncacheNode for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  getExpertForest(owner, 0)->uncacheNode(data[1]);
  getExpertForest(owner, 1)->uncacheNode(data[2]);
}


void
mxd_complement::showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // data[] is of size owner.nParams
  // data[0] --> level
  // data[1,2] <--> forest[0,1]
  // call showNode for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  fprintf(strm, "[%s %d %d %d]",
      owner->op->getName(), data[0], data[1], data[2]);
#if 0
  getExpertForest(owner, 0)->showNode(strm, data[1]);
  getExpertForest(owner, 1)->showNode(strm, data[2]);
#endif
}


bool
mxd_complement::findResult(op_info* owner, int level, int a, int& b)
{
  static int key[2];
  // create cache entry
  key[0] = level;
  key[1] = a;

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  b = cacheEntry[2];
  getExpertForest(owner, 1)->linkNode(b);
  return true;
}


void
mxd_complement::saveResult(op_info* owner, int level, int a, int b)
{
  static int cacheEntry[3];

#ifdef IGNORE_INCOUNT
  if (!getExpertForest(owner, 0)->isTerminalNode(a) &&
      getExpertForest(owner, 0)->getInCount(a) < IGNORE_INCOUNT)
    return;
#endif

  // create cache entry
  cacheEntry[0] = level; cacheEntry[1] = a; cacheEntry[2] = b;

#if 0
#ifdef DEVELOPMENT_CODE
  assert(!findResult(owner, cacheEntry[0], cacheEntry[1], cacheEntry[2]));
#endif
#endif

  getExpertForest(owner, 0)->cacheNode(cacheEntry[1]);
  getExpertForest(owner, 1)->cacheNode(cacheEntry[2]);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}


bool
mxd_complement::checkTerminals(op_info* op, int a, int& b)
{
  if (getExpertForest(op, 0)->isTerminalNode(a)) {
    b = (a == 0)? -1: 0;
    return true;
  }
  return false;
}


// result[i] = a[i]
int mxd_complement::compute(op_info* owner, int a)
{
  return compute(owner,
      getExpertForest(owner, 0)->getExpertDomain()->getTopVariable(), a);
}


int mxd_complement::compute(op_info* owner, int resultLevel, int a)
{
  int result = 0;
  if (checkTerminals(owner, a, result)) {
    if (resultLevel == 0)
      return result;
    if (getExpertForest(owner, 1)->getReductionRule() == forest::FULLY_REDUCED)
      return result;
  }
  if (findResult(owner, resultLevel, a, result))
    return result;

  // expand nodes
  // 0. initialize result
  // 1. copy node a to node result
  // 2. do operation between contents of result and node b

  // 0. initialize result
  // a and result belong to the same domain so the levels must have the
  // same handles and the same sizes
  expert_forest *f0 = getExpertForest(owner, 0);
  expert_forest *f1 = getExpertForest(owner, 1);
  const int nodeALevel = f0->getNodeLevel(a);
  int resultSize = f1->getLevelSize(resultLevel);
  result = f1->createTempNode(resultLevel, resultSize, false);

  int nextLevel = (resultLevel > 0)? -resultLevel:
      f0->getExpertDomain()->getVariableBelow(-resultLevel);

  if (nodeALevel != resultLevel) {
    // All result[i] = compute(nextLevel, a)
    // Special case: expand both levels at skipped identity-reduced level.
    if (resultLevel > 0 &&
        f1->getReductionRule() == forest::IDENTITY_REDUCED) {
      DCASSERT(resultLevel > 0);
      DCASSERT(nextLevel < 0);
      int nextNextLevel = f0->getExpertDomain()->getVariableBelow(-nextLevel);
      int zero = compute(owner, nextNextLevel, 0);
      int temp = compute(owner, nextNextLevel, a);
      for (int i = 0; i < resultSize; ++i) {
        // Build node at nextLevel
        // such that n[i==j] = compute(owner, nextNextLevel, a)
        //       and n[i!=j] = compute(owner, nextNextLevel, 0)
        int n = f1->createTempNode(nextLevel, resultSize, false);
        for (int j = 0; j < resultSize; ++j) {
          f1->setDownPtrWoUnlink(n, j, (i == j)? temp: zero);
        }
        n = f1->reduceNode(n);
        f1->setDownPtrWoUnlink(result, i, n);
        f1->unlinkNode(n);
      }
      f1->unlinkNode(temp);
      f1->unlinkNode(zero);
    }
    else {
      // For Fully and Quasi. Also for prime resultLevel for Identity reduced.
      int temp = compute(owner, nextLevel, a);
      for (int i = 0; i < resultSize; ++i) {
        f1->setDownPtrWoUnlink(result, i, temp);
      }
      f1->unlinkNode(temp);
    }
  }
  else if (f0->isFullNode(a)) {
    // a is a full-node
    const int aSize = f0->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    int i = 0;
    int temp = 0;
    for ( ; i < aSize; ++i) {
      temp = compute(owner, nextLevel, f0->getFullNodeDownPtr(a, i));
      f1->setDownPtrWoUnlink(result, i, temp);
      f1->unlinkNode(temp);
    }
    if (i < resultSize) {
      temp = compute(owner, nextLevel, 0);
      do {
        f1->setDownPtrWoUnlink(result, i++, temp);
      } while (i < resultSize);
      f1->unlinkNode(temp);
    }
  }
  else {
    // a is a sparse-node
    DCASSERT(f0->isSparseNode(a));
    const int aSize = f0->getSparseNodeSize(a);
    DCASSERT(f0->getSparseNodeIndex(a, aSize - 1) < resultSize);
    int i = 0;        // goes from 0..resultSize
    int aIndex = 0;   // j is the sparse index; aIndex is the equiv full index
    int temp = 0;
    int zero = compute(owner, nextLevel, 0);
    for (int j = 0; j < aSize; ++i, ++j) {
      aIndex = f0->getSparseNodeIndex(a, j);
      while (i < aIndex) {
        // a[i] is 0
        f1->setDownPtrWoUnlink(result, i++, zero);
      }
      DCASSERT(i == aIndex && i < resultSize);
      temp = compute(owner, nextLevel, f0->getSparseNodeDownPtr(a, j));
      f1->setDownPtrWoUnlink(result, i, temp);
      f1->unlinkNode(temp);
    }
    while (i < resultSize) {
      // a[i] is 0
      f1->setDownPtrWoUnlink(result, i++, zero);
    }
    f1->unlinkNode(zero);
  }

  // save result in compute cache and return it
  result = f1->reduceNode(result);
  saveResult(owner, resultLevel, a, result);
  return result;
}



// ------------------------------------------------------------------


// ----------------------- MTMDD Apply operation ----------------------


mtmdd_apply_operation::mtmdd_apply_operation()
{ }

mtmdd_apply_operation::~mtmdd_apply_operation() {}

compute_manager::error
mtmdd_apply_operation::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0] != owner->p[1] || owner->p[0] != owner->p[2])
    return compute_manager::FOREST_MISMATCH;
  if (owner->p[0].isMxd() ||
      owner->p[0].isBoolForest() ||
      // the above allows MTMDDs with INTEGERs and REALs
      // the line below allows only INTEGERs
      // owner->f[0]->getRangeType() != forest::INTEGER ||
      !owner->p[0].isMT())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


// ----------------------- MTMDD Divide ----------------------------


mtmdd_divide* mtmdd_divide::getInstance()
{
  static mtmdd_divide instance;
  return &instance;
}


mtmdd_divide::mtmdd_divide()
{ }


mtmdd_divide::~mtmdd_divide() {}


bool
mtmdd_divide::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getInteger(a) /
          getExpertForest(op, 0)->getInteger(b));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getReal(a) /
          getExpertForest(op, 0)->getReal(b));
    }
    return true;
  }
  return false;
}


// ----------------------- MTMDD Multiply ----------------------------


mtmdd_multiply* mtmdd_multiply::getInstance()
{
  static mtmdd_multiply instance;
  return &instance;
}


mtmdd_multiply::mtmdd_multiply()
{ }


mtmdd_multiply::~mtmdd_multiply() {}


bool
mtmdd_multiply::checkTerminals(op_info* op, int a, int b, int& c)
{

  if (a == 0 || b == 0) {
    c = 0;
    return true;
  }

  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getInteger(a) *
          getExpertForest(op, 0)->getInteger(b));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getReal(a) *
          getExpertForest(op, 0)->getReal(b));
    }
    return true;
  }
  return false;
}


// ----------------------- MTMDD Minus ----------------------------


mtmdd_minus* mtmdd_minus::getInstance()
{
  static mtmdd_minus instance;
  return &instance;
}


mtmdd_minus::mtmdd_minus()
{ }


mtmdd_minus::~mtmdd_minus() {}


bool
mtmdd_minus::checkTerminals(op_info* op, int a, int b, int& c)
{

  if (a == b) {
    c = 0;
    return true;
  }

  // TODO: if (a == 0) c = negate(b)

  if (b == 0) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getInteger(a) -
          getExpertForest(op, 0)->getInteger(b));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getReal(a) -
          getExpertForest(op, 0)->getReal(b));
    }
    return true;
  }
  return false;
}


// ----------------------- MTMDD Plus ----------------------------


mtmdd_plus* mtmdd_plus::getInstance()
{
  static mtmdd_plus instance;
  return &instance;
}


mtmdd_plus::~mtmdd_plus() {}


bool
mtmdd_plus::checkTerminals(op_info* op, int a, int b, int& c)
{

#if 1
  if (a == 0 || b == 0) {
    c = a | b;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }
#endif

  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getInteger(a) +
          getExpertForest(op, 0)->getInteger(b));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getReal(a) +
          getExpertForest(op, 0)->getReal(b));
    }
    return true;
  }

  return false;
}

#if 0

mtmdd_plus::mtmdd_plus()
{ }


#else

compute_manager::error
mtmdd_plus::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0] != owner->p[1] || owner->p[0] != owner->p[2])
    return compute_manager::FOREST_MISMATCH;
  if (owner->p[0].isMxd() ||
      owner->p[0].isBoolForest() ||
      // the above allows MTMDDs with INTEGERs and REALs
      // the line below allows only INTEGERs
      // owner->f[0]->getRangeType() != forest::INTEGER ||
      !owner->p[0].isMT())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


mtmdd_plus::mtmdd_plus()
{ }


#endif

// ----------------------- MTMDD Min ----------------------------


mtmdd_min* mtmdd_min::getInstance()
{
  static mtmdd_min instance;
  return &instance;
}


mtmdd_min::mtmdd_min()
{ }


mtmdd_min::~mtmdd_min() {}


bool
mtmdd_min::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = (getExpertForest(op, 0)->getInteger(a) <
          getExpertForest(op, 0)->getInteger(b))? a: b;
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = (getExpertForest(op, 0)->getReal(a) <
          getExpertForest(op, 0)->getReal(b))? a: b;
    }
    return true;
  }
  return false;
}


// ----------------------- MTMDD Max ----------------------------


mtmdd_max* mtmdd_max::getInstance()
{
  static mtmdd_max instance;
  return &instance;
}


mtmdd_max::mtmdd_max()
{ }


mtmdd_max::~mtmdd_max() {}


bool
mtmdd_max::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = (getExpertForest(op, 0)->getInteger(a) >
          getExpertForest(op, 0)->getInteger(b))? a: b;
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = (getExpertForest(op, 0)->getReal(a) >
          getExpertForest(op, 0)->getReal(b))? a: b;
#if 0
      fprintf(stderr, "a: %d, b: %d, c: %d\n", a, b, c);
      fprintf(stderr, "a: %f, b: %f, c: %f\n",
          getExpertForest(op, 0)->getReal(a),
          getExpertForest(op, 0)->getReal(b),
          getExpertForest(op, 0)->getReal(c));
#endif
    }
    return true;
  }

  return false;
}


// ----------------------- MTMDD Or-Min ----------------------------


mtmdd_or_min* mtmdd_or_min::getInstance()
{
  static mtmdd_or_min instance;
  return &instance;
}


mtmdd_or_min::mtmdd_or_min()
{ }


mtmdd_or_min::~mtmdd_or_min() {}


bool
mtmdd_or_min::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (0 == a) {
    c = b;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (0 == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (a == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = (getExpertForest(op, 0)->getInteger(a) <
          getExpertForest(op, 0)->getInteger(b))? a: b;
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = (getExpertForest(op, 0)->getReal(a) <
          getExpertForest(op, 0)->getReal(b))? a: b;
    }
    return true;
  }


  return false;
}


// ----------------------- MTMDD Or-Max ----------------------------


mtmdd_or_max* mtmdd_or_max::getInstance()
{
  static mtmdd_or_max instance;
  return &instance;
}


mtmdd_or_max::mtmdd_or_max()
{ }


mtmdd_or_max::~mtmdd_or_max() {}


bool
mtmdd_or_max::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (0 == a) {
    c = b;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (0 == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (a == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = (getExpertForest(op, 0)->getInteger(a) >
          getExpertForest(op, 0)->getInteger(b))? a: b;
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = (getExpertForest(op, 0)->getReal(a) >
          getExpertForest(op, 0)->getReal(b))? a: b;
    }
    return true;
  }

  return false;
}


// ----------------------- MTMDD And-Min ----------------------------


mtmdd_and_min* mtmdd_and_min::getInstance()
{
  static mtmdd_and_min instance;
  return &instance;
}


mtmdd_and_min::mtmdd_and_min()
{ }


mtmdd_and_min::~mtmdd_and_min() {}


bool
mtmdd_and_min::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (0 == a || 0 == b) {
    c = 0;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (a == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = (getExpertForest(op, 0)->getInteger(a) <
          getExpertForest(op, 0)->getInteger(b))? a: b;
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = (getExpertForest(op, 0)->getReal(a) <
          getExpertForest(op, 0)->getReal(b))? a: b;
    }
    return true;
  }

  return false;
}


// ----------------------- MTMDD And-Max ----------------------------


mtmdd_and_max* mtmdd_and_max::getInstance()
{
  static mtmdd_and_max instance;
  return &instance;
}


mtmdd_and_max::mtmdd_and_max()
{ }


mtmdd_and_max::~mtmdd_and_max() {}


bool
mtmdd_and_max::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (0 == a || 0 == b) {
    c = 0;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (a == b) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = (getExpertForest(op, 0)->getInteger(a) >
          getExpertForest(op, 0)->getInteger(b))? a: b;
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = (getExpertForest(op, 0)->getReal(a) >
          getExpertForest(op, 0)->getReal(b))? a: b;
    }
    return true;
  }

  return false;
}


// ----------------------- MTMDD Less-Than ----------------------------


mtmdd_less_than* mtmdd_less_than::getInstance()
{
  static mtmdd_less_than instance;
  return &instance;
}


mtmdd_less_than::mtmdd_less_than()
{ }


mtmdd_less_than::~mtmdd_less_than() {}


bool
mtmdd_less_than::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          (getExpertForest(op, 0)->getInteger(a) <
           getExpertForest(op, 0)->getInteger(b))? 1: 0);
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          (getExpertForest(op, 0)->getReal(a) <
           getExpertForest(op, 0)->getReal(b))? 1.0f: 0.0f);
    }
    return true;
  }

  return false;
}


// ---------------------- MTMDD Less-Than-Or-Equal-To ------------------------


mtmdd_less_than_equal* mtmdd_less_than_equal::getInstance()
{
  static mtmdd_less_than_equal instance;
  return &instance;
}


mtmdd_less_than_equal::mtmdd_less_than_equal()
{ }


mtmdd_less_than_equal::~mtmdd_less_than_equal() {}


bool
mtmdd_less_than_equal::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          (getExpertForest(op, 0)->getInteger(a) <=
           getExpertForest(op, 0)->getInteger(b))? 1: 0);
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          (getExpertForest(op, 0)->getReal(a) <=
           getExpertForest(op, 0)->getReal(b))? 1.0f: 0.0f);
    }
    return true;
  }

  return false;
}


// ----------------------- MTMDD Greater-Than ----------------------------


mtmdd_greater_than* mtmdd_greater_than::getInstance()
{
  static mtmdd_greater_than instance;
  return &instance;
}


mtmdd_greater_than::mtmdd_greater_than()
{ }


mtmdd_greater_than::~mtmdd_greater_than() {}


bool
mtmdd_greater_than::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          (getExpertForest(op, 0)->getInteger(a) >
           getExpertForest(op, 0)->getInteger(b))? 1: 0);
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          (getExpertForest(op, 0)->getReal(a) >
           getExpertForest(op, 0)->getReal(b))? 1.0f: 0.0f);
    }
    return true;
  }

  return false;
}


// -------------------- MTMDD Greater-Than-Or-Equal-To ----------------------


mtmdd_greater_than_equal* mtmdd_greater_than_equal::getInstance()
{
  static mtmdd_greater_than_equal instance;
  return &instance;
}


mtmdd_greater_than_equal::mtmdd_greater_than_equal()
{ }


mtmdd_greater_than_equal::~mtmdd_greater_than_equal() {}


bool
mtmdd_greater_than_equal::checkTerminals(op_info* op,
    int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() ==
        forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          (getExpertForest(op, 0)->getInteger(a) >=
           getExpertForest(op, 0)->getInteger(b))? 1: 0);
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() ==
          forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          (getExpertForest(op, 0)->getReal(a) >=
           getExpertForest(op, 0)->getReal(b))? 1.0f: 0.0f);
    }
    return true;
  }
  return false;
}


// -------------------- MTMDD Equal-To ----------------------


mtmdd_equal* mtmdd_equal::getInstance()
{
  static mtmdd_equal instance;
  return &instance;
}


mtmdd_equal::mtmdd_equal()
{ }


mtmdd_equal::~mtmdd_equal() {}


bool
mtmdd_equal::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(int((a == b)? 1: 0));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(float((a == b)? 1: 0));
    }
    return true;
  }
  return false;
}


// -------------------- MTMDD Not-Equal-To ----------------------


mtmdd_not_equal* mtmdd_not_equal::getInstance()
{
  static mtmdd_not_equal instance;
  return &instance;
}


mtmdd_not_equal::mtmdd_not_equal()
{ }


mtmdd_not_equal::~mtmdd_not_equal() {}


bool
mtmdd_not_equal::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(int((a != b)? 1: 0));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(float((a != b)? 1: 0));
    }
    return true;
  }
  return false;
}


// ----------------------- MTMDD-To-MDD Apply operation ----------------------


mtmdd_to_mdd_apply_operation::mtmdd_to_mdd_apply_operation()
: scratch(0)
{ }


void mtmdd_to_mdd_apply_operation::buildScratch(const expert_domain *d)
{
  // build new scratch:
  // each level in scratch corresponds to a level in the domain
  // each level[i] is of size levelSizes[i] == getVariableBound(i)

  int topVar = d->getTopVariable();

  // number of levels in the domain (incl. terminals)
  nLevels = topVar + 1;

  // allocate levelSizes int array
  levelSizes = (int *) malloc(nLevels * sizeof(int));
  assert(levelSizes != 0);
  memset(levelSizes, 0, nLevels * sizeof(int));

  // allocate scratch pointer array
  scratch = (int **) malloc(nLevels * sizeof(int*));
  assert(scratch != 0);
  memset(scratch, 0, nLevels * sizeof(int*));
  
  // scratch 2-D array: no need to alloc for each level,
  // just the ones that are valid (levelSizes[i] != 0) in this forest
  for (unsigned i = 0u; i < nLevels; ++i) {
    levelSizes[i] = d->getVariableBound(i);
    if (levelSizes[i] > 0)
      scratch[i] = (int *) malloc(levelSizes[i] * sizeof(int));
    // no need to clear the array: scratch == always assumed to be dirty!
  }
}


void mtmdd_to_mdd_apply_operation::initScratch(const expert_domain* d)
{
  // build scratch 2-D array based on domain

  if (scratch != 0) {
    // if its the same size, leave it alone otherwise rebuild it
    if (nLevels == (unsigned int)(1 + d->getTopVariable())) {
      // now check if the level sizes are the same
      unsigned i = 0;
      for ( ; i < nLevels && levelSizes[i] == d->getVariableBound(i); ++i);
      if (i == nLevels) return;
    }

    // get rid of previous scratch
    deleteScratch();
  }

  // build new scratch
  buildScratch(d);
}


void mtmdd_to_mdd_apply_operation::deleteScratch()
{
  if (scratch == 0) return;
  DCASSERT(levelSizes != 0);
  for (unsigned i = 0u; i < nLevels; ++i)
    if (scratch[i]) free(scratch[i]);
  free(scratch);
  free(levelSizes);
}


int* mtmdd_to_mdd_apply_operation::getScratch(int level)
{
  DCASSERT(level >= 0);
  DCASSERT((unsigned int)level < nLevels);
  DCASSERT(scratch[level] != 0);
  return scratch[level];
}
 

void mtmdd_to_mdd_apply_operation::setScratch(int level, int val)
{
  DCASSERT(level >= 0);
  DCASSERT((unsigned int)level < nLevels);
  int* data = getScratch(level);
  int* stop = data + levelSizes[level];
  while (data != stop) {
    *data = val;
    data++;
  }
}
 

mtmdd_to_mdd_apply_operation::~mtmdd_to_mdd_apply_operation()
{
  deleteScratch();
}


compute_manager::error
mtmdd_to_mdd_apply_operation::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0] != owner->p[1])
    return compute_manager::FOREST_MISMATCH;
  if (!getExpertForest(owner, 0)->isMtMdd() ||
      !getExpertForest(owner, 2)->isMdd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


#if 1

// Two-step compute
// - result[i] = a[i]
// - result[i] = result[i] op b[i]
// This makes for a simpler but less efficient implementation.
int mtmdd_to_mdd_apply_operation::compute(op_info* owner, int a, int b)
{
  initScratch(smart_cast<const expert_domain*>(owner->p[0].getDomain()));
  return computeHelper(owner, a, b);
}

int mtmdd_to_mdd_apply_operation::computeHelper(op_info* owner, int a, int b)
{
  DCASSERT(getExpertForest(owner, 0) == getExpertForest(owner, 1));

  int result = 0;
  if (checkTerminals(owner, a, b, result)) return result;
  if (findResult(owner, a, b, result)) return result;

  expert_forest* mtmdd = getExpertForest(owner, 0);
  expert_forest* mdd = getExpertForest(owner, 2);

  // expand nodes
  // 0. initialize result
  // 1. copy node a to node result
  // 2. do operation between contents of result and node b

  // 0. initialize result
  const int aLevel = mtmdd->getNodeLevel(a);
  const int bLevel = mtmdd->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = mtmdd->getLevelSize(resultLevel);
  int* scratch = getScratch(resultLevel);

  DCASSERT(mtmdd->getLevelSize(resultLevel) ==
      mdd->getLevelSize(resultLevel));

  // 1. copy node a to node result
  if (aLevel < resultLevel) {
    // all down pointers of result point to node a
    for (int i = 0; i < resultSize; ++i)
      scratch[i] = a;
  }
  else if (mtmdd->isFullNode(a)) {
    // a is a full-node
    const int aSize = mtmdd->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    for (int i = 0; i < aSize; ++i)
      scratch[i] = mtmdd->getFullNodeDownPtr(a, i);
    for (int i = aSize; i < resultSize; ++i)
      scratch[i] = 0;
  }
  else {
    // a is a sparse-node
    memset(scratch, 0, resultSize * sizeof(int));
    const int aSize = mtmdd->getSparseNodeSize(a);
    for (int i = 0; i < aSize; ++i)
      scratch[mtmdd->getSparseNodeIndex(a, i)] =
        mtmdd->getSparseNodeDownPtr(a, i);
  }

  // 2. do operation between contents of result and node b
  if (bLevel < resultLevel) {
    for (int i = 0; i < resultSize; ++i)
      scratch[i] = computeHelper(owner, scratch[i], b);
  }
  else if (mtmdd->isFullNode(b)) {
    // b is a full-node
    const int bSize = mtmdd->getFullNodeSize(b);
    DCASSERT(bSize <= resultSize);
    for (int i = 0; i < bSize; ++i)
      scratch[i] = computeHelper(owner, scratch[i],
          mtmdd->getFullNodeDownPtr(b, i));
    for (int i = bSize; i < resultSize; ++i)
      scratch[i] = computeHelper(owner, scratch[i], 0);
  }
  else {
    // b is a sparse-node
    const int bSize = mtmdd->getSparseNodeSize(b);
    DCASSERT(mtmdd->getSparseNodeIndex(b, bSize - 1) <= resultSize);
    // j goes through every index (like a full-node index pointer)
    int j = 0;
    for (int i = 0; i < bSize; ++i, ++j)
    {
      // sparse-nodes skip indices which represent downpointer 0
      // call compute of those skipped indices
      for (int index = mtmdd->getSparseNodeIndex(b, i); j < index; ++j)
        scratch[j] = computeHelper(owner, scratch[j], 0);
      // done with skipped indices; deal with the next sparse node index
      DCASSERT(j == mtmdd->getSparseNodeIndex(b, i));
      scratch[j] = computeHelper(owner, scratch[j],
          mtmdd->getSparseNodeDownPtr(b, i));
    }
    DCASSERT(j == mtmdd->getSparseNodeIndex(b, bSize - 1) + 1);
    for ( ; j < resultSize; ++j)
      scratch[j] = computeHelper(owner, scratch[j], 0);
  }

  // create node representing scratch[]
  result = mdd->createTempNode(resultLevel, resultSize);
  for (int i = 0; i < resultSize; ++i)
  {
    mdd->setDownPtrWoUnlink(result, i, scratch[i]);
    mdd->unlinkNode(scratch[i]);
  }
  result = mdd->reduceNode(result);

  // save result in compute cache and return it
  saveResult(owner, a, b, result);
  return result;
}

#endif

// ----------------------- MTMDD-MDD Less-Than ----------------------------


mtmdd_to_mdd_less_than* mtmdd_to_mdd_less_than::getInstance()
{
  static mtmdd_to_mdd_less_than instance;
  return &instance;
}


mtmdd_to_mdd_less_than::mtmdd_to_mdd_less_than()
{ }


mtmdd_to_mdd_less_than::~mtmdd_to_mdd_less_than() {}


bool
mtmdd_to_mdd_less_than::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    c = getExpertForest(op, 2)->getTerminalNode(
        (getExpertForest(op, 0)->getRangeType() == forest::INTEGER)
        ? (getExpertForest(op, 0)->getInteger(a) <
          getExpertForest(op, 0)->getInteger(b)? true: false)
        : (getExpertForest(op, 0)->getReal(a) <
          getExpertForest(op, 0)->getReal(b)? true: false));
    return true;
  }

  return false;
}


// ---------------------- MTMDD-MDD Less-Than-Or-Equal-To --------------------


mtmdd_to_mdd_less_than_equal* mtmdd_to_mdd_less_than_equal::getInstance()
{
  static mtmdd_to_mdd_less_than_equal
    instance;
  return &instance;
}


mtmdd_to_mdd_less_than_equal::mtmdd_to_mdd_less_than_equal()
{ }


mtmdd_to_mdd_less_than_equal::~mtmdd_to_mdd_less_than_equal() {}


bool
mtmdd_to_mdd_less_than_equal::checkTerminals(op_info* op,
    int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    c = getExpertForest(op, 2)->getTerminalNode(
        (getExpertForest(op, 0)->getRangeType() == forest::INTEGER)
        ? (getExpertForest(op, 0)->getInteger(a) <=
          getExpertForest(op, 0)->getInteger(b)? true: false)
        : (getExpertForest(op, 0)->getReal(a) <=
          getExpertForest(op, 0)->getReal(b)? true: false));
    return true;
  }
  return false;
}


// ----------------------- MTMDD-MDD Greater-Than --------------------------


mtmdd_to_mdd_greater_than* mtmdd_to_mdd_greater_than::getInstance()
{
  static mtmdd_to_mdd_greater_than instance;
  return &instance;
}


mtmdd_to_mdd_greater_than::mtmdd_to_mdd_greater_than()
{ }


mtmdd_to_mdd_greater_than::~mtmdd_to_mdd_greater_than() {}


bool
mtmdd_to_mdd_greater_than::checkTerminals(op_info* op,
    int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    c = getExpertForest(op, 2)->getTerminalNode(
        (getExpertForest(op, 0)->getRangeType() == forest::INTEGER)
        ? (getExpertForest(op, 0)->getInteger(a) >
          getExpertForest(op, 0)->getInteger(b)? true: false)
        : (getExpertForest(op, 0)->getReal(a) >
          getExpertForest(op, 0)->getReal(b)? true: false));
    return true;
  }

  return false;
}


// -------------------- MTMDD-MDD Greater-Than-Or-Equal-To ------------------


mtmdd_to_mdd_greater_than_equal*
mtmdd_to_mdd_greater_than_equal::getInstance()
{
  static mtmdd_to_mdd_greater_than_equal
    instance;
  return &instance;
}


mtmdd_to_mdd_greater_than_equal
::mtmdd_to_mdd_greater_than_equal()
{ }


mtmdd_to_mdd_greater_than_equal::~mtmdd_to_mdd_greater_than_equal() {}


bool
mtmdd_to_mdd_greater_than_equal::checkTerminals(op_info* op,
    int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    c = getExpertForest(op, 2)->getTerminalNode(
        (getExpertForest(op, 0)->getRangeType() == forest::INTEGER)
        ? (getExpertForest(op, 0)->getInteger(a) >=
          getExpertForest(op, 0)->getInteger(b)? true: false)
        : (getExpertForest(op, 0)->getReal(a) >=
          getExpertForest(op, 0)->getReal(b)? true: false));
    return true;
  }

  return false;
}


// -------------------- MTMDD-MDD Equal-To ------------------


mtmdd_to_mdd_equal* mtmdd_to_mdd_equal::getInstance()
{
  static mtmdd_to_mdd_equal instance;
  return &instance;
}


mtmdd_to_mdd_equal::mtmdd_to_mdd_equal()
{ }


mtmdd_to_mdd_equal::~mtmdd_to_mdd_equal() {}


bool mtmdd_to_mdd_equal::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    c = getExpertForest(op, 2)->getTerminalNode(bool(a == b? true: false));
    return true;
  }
return false;
}


// -------------------- MTMDD-MDD Not-Equal-To ------------------


mtmdd_to_mdd_not_equal* mtmdd_to_mdd_not_equal::getInstance()
{
  static mtmdd_to_mdd_not_equal instance;
  return &instance;
}


mtmdd_to_mdd_not_equal::mtmdd_to_mdd_not_equal()
{ }


mtmdd_to_mdd_not_equal::~mtmdd_to_mdd_not_equal() {}


bool mtmdd_to_mdd_not_equal::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    c = getExpertForest(op, 2)->getTerminalNode(bool(a != b? true: false));
    return true;
  }
  return false;
}


// -------------------- Conversion operations ------------------


conversion_operation::conversion_operation()
{}

conversion_operation::~conversion_operation()
{}

compute_manager::error
conversion_operation::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 2)
    return compute_manager::WRONG_NUMBER;
  // No assumption about what kind of forests these are, but they
  // have to be forests in the same domain
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool conversion_operation::isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[i] <--> forest[i]
  // call isStale for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  return 
    getExpertForest(owner, 0)->isStale(data[0]) ||
    getExpertForest(owner, 1)->isStale(data[1]);
}


void
conversion_operation::discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[i] <--> forest[i]
  // call uncacheNode for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  getExpertForest(owner, 0)->uncacheNode(data[0]);
  getExpertForest(owner, 1)->uncacheNode(data[1]);
}


void
conversion_operation::showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // TODO:
  // data[] is of size owner.nParams
  // data[i] <--> forest[i]
  // call showNode for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  fprintf(strm, "[%s %d %d]",
      owner->op->getName(), data[0], data[1]);
#if 0
  getExpertForest(owner, 0)->showNode(strm, data[0]);
  getExpertForest(owner, 1)->showNode(strm, data[1]);
#endif
}


bool
conversion_operation::findResult(op_info* owner, int a, int& b)
{
  static int key[1];
  // create cache entry
  key[0] = a;

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  b = cacheEntry[1];
  getExpertForest(owner, 1)->linkNode(b);
  return true;
}


void
conversion_operation::saveResult(op_info* owner, int a, int b)
{
  static int cacheEntry[2];

#ifdef IGNORE_INCOUNT
  if (!getExpertForest(owner, 0)->isTerminalNode(a) &&
      getExpertForest(owner, 0)->getInCount(a) < IGNORE_INCOUNT)
    return;
#endif

  // create cache entry
  cacheEntry[0] = a; cacheEntry[1] = b;

#if 0
#ifdef DEVELOPMENT_CODE
  assert(!findResult(owner, cacheEntry[0], cacheEntry[1]));
#endif
#endif

  getExpertForest(owner, 0)->cacheNode(cacheEntry[0]);
  getExpertForest(owner, 1)->cacheNode(cacheEntry[1]);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}


compute_manager::error
conversion_operation::compute(op_info* owner, dd_edge** operands)
{
  if (operands == 0) return compute_manager::TYPE_MISMATCH;
  // compute(owner, dd_edge, dd_edge, dd_edge) checks for owner == 0
  return compute(owner, *operands[0], *operands[1]);
}

compute_manager::error
conversion_operation::compute(op_info* owner, const dd_edge& a, dd_edge& b)
{
  if (owner == 0) return compute_manager::TYPE_MISMATCH;
  int result = compute(owner, a.getNode());
  b.set(result, 0, getExpertForest(owner, 1)->getNodeLevel(result));
  return compute_manager::SUCCESS;
}

compute_manager::error
conversion_operation::compute(op_info* owner, const dd_edge& a,
    const dd_edge& b, dd_edge& c)
{
  return compute_manager::TYPE_MISMATCH;
}


// result[i] = a[i]
int conversion_operation::compute(op_info* owner, int a)
{
  int result = 0;
  if (checkTerminals(owner, a, result))
    return result;
  if (findResult(owner, a, result))
    return result;

  // expand nodes
  // 0. initialize result
  // 1. copy node a to node result
  // 2. do operation between contents of result and node b

  // 0. initialize result
  // a and result belong to the same domain so the levels must have the
  // same handles and the same sizes
  const int resultLevel = getExpertForest(owner, 0)->getNodeLevel(a);
  int resultSize = getExpertForest(owner, 1)->getLevelSize(resultLevel);
  result =
    getExpertForest(owner, 1)->createTempNode(resultLevel, resultSize, false);

  // 1. copy node a to node result
  if (getExpertForest(owner, 0)->isFullNode(a)) {
    // a is a full-node
    const int aSize = getExpertForest(owner, 0)->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    int i = 0;
    int temp = 0;
    for ( ; i < aSize; ++i) {
      temp =
        compute(owner, getExpertForest(owner, 0)->getFullNodeDownPtr(a, i));
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, temp);
      getExpertForest(owner, 1)->unlinkNode(temp);
    }
    for ( ; i < resultSize; ++i) {
      temp = compute(owner, 0);
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, temp);
      getExpertForest(owner, 1)->unlinkNode(temp);
    }
  }
  else {
    // a is a sparse-node
    const int aSize =
      getExpertForest(owner, 0)->getSparseNodeSize(a);

    DCASSERT(getExpertForest(owner, 0)->getSparseNodeIndex(a, aSize - 1)
        < resultSize);

    int i = 0;        // goes from 0..resultSize
    int aIndex = 0;   // j is the sparse index; aIndex is the equiv full index
    int temp = 0;

    for (int j = 0; j < aSize; ++i, ++j) {
      aIndex = getExpertForest(owner, 0)->getSparseNodeIndex(a, j);
      for ( ; i < aIndex; ++i) {
        // a[i] is 0
        temp = compute(owner, 0);
        getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, temp);
        getExpertForest(owner, 1)->unlinkNode(temp);
      }
      DCASSERT(i == aIndex && i < resultSize);
      temp =
        compute(owner, getExpertForest(owner, 0)->getSparseNodeDownPtr(a, j));
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, temp);
      getExpertForest(owner, 1)->unlinkNode(temp);
    }

    for ( ; i < resultSize; ++i) {
      temp = compute(owner, 0);
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, temp);
      getExpertForest(owner, 1)->unlinkNode(temp);
    }
  }

  // save result in compute cache and return it
  result = getExpertForest(owner, 1)->reduceNode(result);
  saveResult(owner, a, result);
  return result;
}


// -------------------- MDD Complement ------------------


mdd_complement* mdd_complement::getInstance()
{
  static mdd_complement instance;
  return &instance;
}


mdd_complement::mdd_complement()
{ }


mdd_complement::~mdd_complement() {}


compute_manager::error
mdd_complement::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 0)->isMdd())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 1)->isMdd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mdd_complement::checkTerminals(op_info* op, int a, int& b)
{
  if (getExpertForest(op, 0)->isTerminalNode(a)) {
    b = (a == 0)? -1: 0;
    return true;
  }
  return false;
}


// -------------------- Convert MTMDD to MDD ------------------


mtmdd_to_mdd* mtmdd_to_mdd::getInstance()
{
  static mtmdd_to_mdd instance;
  return &instance;
}


mtmdd_to_mdd::mtmdd_to_mdd()
{ }


mtmdd_to_mdd::~mtmdd_to_mdd() {}


compute_manager::error
mtmdd_to_mdd::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 0)->isMtMdd())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 1)->isMdd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mtmdd_to_mdd::checkTerminals(op_info* op, int a, int& b)
{
  if (getExpertForest(op, 0)->isTerminalNode(a)) {
    b = getExpertForest(op, 1)->getTerminalNode(bool(
        (getExpertForest(op, 0)->getRangeType() == forest::INTEGER)
        ? (getExpertForest(op, 0)->getInteger(a) == 0 ? false: true)
        : (getExpertForest(op, 0)->getReal(a) == 0.0 ? false: true)));
    return true;
  }

  return false;
}


// -------------------- Convert MDD to MTMDD ------------------


mdd_to_mtmdd* mdd_to_mtmdd::getInstance()
{
  static mdd_to_mtmdd instance;
  return &instance;
}


mdd_to_mtmdd::mdd_to_mtmdd()
{ }


mdd_to_mtmdd::~mdd_to_mtmdd() {}


compute_manager::error
mdd_to_mtmdd::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 0)->isMdd())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 1)->isMtMdd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mdd_to_mtmdd::checkTerminals(op_info* op, int a, int& b)
{
  if (getExpertForest(op, 0)->isTerminalNode(a)) {
    b = getExpertForest(op, 1)->getRangeType() == forest::INTEGER
      ? getExpertForest(op, 1)->getTerminalNode(
          int(getExpertForest(op, 0)->getBoolean(a)? 1: 0))
      : getExpertForest(op, 1)->getTerminalNode(
          float(getExpertForest(op, 0)->getBoolean(a)? 1.0f: 0.0f));
    return true;
  }

  return false;
}


// -------------------- Convert MTMXD to MXD ------------------


mtmxd_to_mxd* mtmxd_to_mxd::getInstance()
{
  static mtmxd_to_mxd instance;
  return &instance;
}


mtmxd_to_mxd::mtmxd_to_mxd()
{ }


mtmxd_to_mxd::~mtmxd_to_mxd() {}


compute_manager::error
mtmxd_to_mxd::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 0)->isMtMxd())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 1)->isMxd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mtmxd_to_mxd::checkTerminals(op_info* op, int a, int& b)
{
  if (getExpertForest(op, 0)->isTerminalNode(a)) {
    b = getExpertForest(op, 1)->getTerminalNode(bool(
        (getExpertForest(op, 0)->getRangeType() == forest::INTEGER)
        ? (getExpertForest(op, 0)->getInteger(a) == 0 ? false: true)
        : (getExpertForest(op, 0)->getReal(a) == 0.0 ? false: true)));
    return true;
  }

  return false;
}


// -------------------- Convert MXD to MTMXD ------------------


mxd_to_mtmxd* mxd_to_mtmxd::getInstance()
{
  static mxd_to_mtmxd instance;
  return &instance;
}


mxd_to_mtmxd::mxd_to_mtmxd()
{ }


mxd_to_mtmxd::~mxd_to_mtmxd() {}


compute_manager::error
mxd_to_mtmxd::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 0)->isMxd())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 1)->isMtMxd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mxd_to_mtmxd::checkTerminals(op_info* op, int a, int& b)
{
  if (getExpertForest(op, 0)->isTerminalNode(a)) {
    b = getExpertForest(op, 1)->getRangeType() == forest::INTEGER
      ? getExpertForest(op, 1)->getTerminalNode(
          int(getExpertForest(op, 0)->getBoolean(a)? 1: 0))
      : getExpertForest(op, 1)->getTerminalNode(
          float(getExpertForest(op, 0)->getBoolean(a)? 1.0f: 0.0f));
    return true;
  }

  return false;
}


// -------------------- Convert MTMDD to EVMDD ------------------


mtmdd_to_evmdd* mtmdd_to_evmdd::getInstance()
{
  static mtmdd_to_evmdd instance;
  return &instance;
}


mtmdd_to_evmdd::mtmdd_to_evmdd()
{ }


mtmdd_to_evmdd::~mtmdd_to_evmdd() {}


compute_manager::error
mtmdd_to_evmdd::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 0)->isMtMdd())
    return compute_manager::TYPE_MISMATCH;
  if (owner->p[1].isMT())
    // !MULTI_TERMINAL is same as EV+ or EV*
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mtmdd_to_evmdd::isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nParams + 1 (last int is edge value)
  // data[i] <--> forest[i]
  // call isStale for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  return 
    getExpertForest(owner, 0)->isStale(data[0]) ||
    getExpertForest(owner, 1)->isStale(data[1]);
}


void
mtmdd_to_evmdd::discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams + 1 (last int is edge value)
  // data[i] <--> forest[i]
  // call uncacheNode for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  getExpertForest(owner, 0)->uncacheNode(data[0]);
  getExpertForest(owner, 1)->uncacheNode(data[1]);
}


void
mtmdd_to_evmdd::showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // data[] is of size owner.nParams + 1 (last int is edge value)
  // data[i] <--> forest[i]
  // call showNode for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  if (getExpertForest(owner, 1)->getRangeType() == forest::REAL) {
    fprintf(strm, "[%s %d %d %f]",
        owner->op->getName(), data[0], data[1], toFloat(data[2]));
  } else {
    fprintf(strm, "[%s %d %d %d]", owner->op->getName(),
        data[0], data[1], data[2]);
  }
}


bool
mtmdd_to_evmdd::findResult(op_info* owner, int a, int& b, int& ev)
{
  static int key[1];
  // create cache entry
  key[0] = a;

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  b = cacheEntry[1];
  ev = cacheEntry[2];
  getExpertForest(owner, 1)->linkNode(b);
  return true;
}


bool
mtmdd_to_evmdd::findResult(op_info* owner, int a, int& b, float& ev)
{
  static int key[1];
  // create cache entry
  key[0] = a;

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  b = cacheEntry[1];
  ev = toFloat(cacheEntry[2]);
  getExpertForest(owner, 1)->linkNode(b);
  return true;
}


void
mtmdd_to_evmdd::saveResult(op_info* owner, int a, int b, int c)
{
  static int cacheEntry[3];

#ifdef IGNORE_INCOUNT
  if (!getExpertForest(owner, 0)->isTerminalNode(a) &&
      getExpertForest(owner, 0)->getInCount(a) < IGNORE_INCOUNT)
    return;
#endif

  // create cache entry
  cacheEntry[0] = a; cacheEntry[1] = b; cacheEntry[2] = c;

  getExpertForest(owner, 0)->cacheNode(a);
  getExpertForest(owner, 1)->cacheNode(b);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}


void
mtmdd_to_evmdd::saveResult(op_info* owner, int a, int b, float c)
{
  static int cacheEntry[3];

#ifdef IGNORE_INCOUNT
  if (!getExpertForest(owner, 0)->isTerminalNode(a) &&
      getExpertForest(owner, 0)->getInCount(a) < IGNORE_INCOUNT)
    return;
#endif

  // create cache entry
  cacheEntry[0] = a; cacheEntry[1] = b; cacheEntry[2] = toInt(c);

  getExpertForest(owner, 0)->cacheNode(a);
  getExpertForest(owner, 1)->cacheNode(b);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}


bool mtmdd_to_evmdd::checkTerminals(op_info* op, int a, int& b, int& ev)
{
  if (getExpertForest(op, 0)->isTerminalNode(a)) {
    ev = getExpertForest(op, 0)->getInteger(a);
    if (ev == 0) {
      // replace with non-reachable node
      b = getExpertForest(op, 1)->getTerminalNode(false);
      ev = INF;
    } else {
      b = getExpertForest(op, 1)->getTerminalNode(true);
    }
    return true;
  }

  return false;
}


bool mtmdd_to_evmdd::checkTerminals(op_info* op, int a, int& b, float& ev)
{
  if (getExpertForest(op, 0)->isTerminalNode(a)) {
    b = getExpertForest(op, 1)->getTerminalNode(true);
    ev = getExpertForest(op, 0)->getReal(a);
    return true;
  }

  return false;
}


compute_manager::error
mtmdd_to_evmdd::compute(op_info* owner, dd_edge** operands)
{
  if (operands == 0) return compute_manager::TYPE_MISMATCH;
  // compute(owner, dd_edge, dd_edge, dd_edge) checks for owner == 0
  return compute(owner, *operands[0], *operands[1]);
}


compute_manager::error
mtmdd_to_evmdd::compute(op_info* owner, const dd_edge& a,
    const dd_edge& b, dd_edge& c)
{
  return compute_manager::TYPE_MISMATCH;
}


compute_manager::error
mtmdd_to_evmdd::compute(op_info* owner, const dd_edge& a, dd_edge& b)
{
  if (owner == 0) return compute_manager::TYPE_MISMATCH;
  int result = 0;
  int ev = 0;
  compute(owner, a.getNode(), result, ev);
  b.set(result, ev, getExpertForest(owner, 1)->getNodeLevel(result));
  return compute_manager::SUCCESS;
}


// result[i] = a[i]
void mtmdd_to_evmdd::compute(op_info* owner, int a, int& b, int &bev)
{
  b = 0; bev = 0;
  if (checkTerminals(owner, a, b, bev)) return;
  if (findResult(owner, a, b, bev)) return;

  // expand nodes
  // 0. initialize result
  // 1. copy node a to node result
  // 2. do operation between contents of result and node b

  // 0. initialize result
  // a and result belong to the same domain so the levels must have the
  // same handles and the same sizes
  const int resultLevel = getExpertForest(owner, 0)->getNodeLevel(a);
  int resultSize = getExpertForest(owner, 1)->getLevelSize(resultLevel);
  int result =
    getExpertForest(owner, 1)->createTempNode(resultLevel, resultSize, false);

  // 1. copy node a to node result
  if (getExpertForest(owner, 0)->isFullNode(a)) {
    // a is a full-node
    const int aSize = getExpertForest(owner, 0)->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    int i = 0;
    int node = 0;
    int ev = 0;
    for ( ; i < aSize; ++i) {
      node = ev = 0;
      compute(owner, getExpertForest(owner, 0)->getFullNodeDownPtr(a, i),
          node, ev);
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, node);
      getExpertForest(owner, 1)->setEdgeValue(result, i, ev);
      getExpertForest(owner, 1)->unlinkNode(node);
    }
    node = ev = 0;
    compute(owner, 0, node, ev);
    for ( ; i < resultSize; ++i) {
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, node);
      getExpertForest(owner, 1)->setEdgeValue(result, i, ev);
    }
    getExpertForest(owner, 1)->unlinkNode(node);
  }
  else {
    // a is a sparse-node
    const int aSize =
      getExpertForest(owner, 0)->getSparseNodeSize(a);

    DCASSERT(getExpertForest(owner, 0)->getSparseNodeIndex(a, aSize - 1)
        < resultSize);

    int i = 0;        // goes from 0..resultSize
    int aIndex = 0;   // j is the sparse index; aIndex is the equiv full index
    int node = 0;
    int ev = 0;

    for (int j = 0; j < aSize; ++i, ++j) {
      aIndex = getExpertForest(owner, 0)->getSparseNodeIndex(a, j);
      for ( ; i < aIndex; ++i) {
        // a[i] is 0
        node = ev = 0;
        compute(owner, 0, node, ev);
        getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, node);
        getExpertForest(owner, 1)->setEdgeValue(result, i, ev);
        getExpertForest(owner, 1)->unlinkNode(node);
      }
      DCASSERT(i == aIndex && i < resultSize);
      node = ev = 0;
      compute(owner, getExpertForest(owner, 0)->getSparseNodeDownPtr(a, j),
          node, ev);
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, node);
      getExpertForest(owner, 1)->setEdgeValue(result, i, ev);
      getExpertForest(owner, 1)->unlinkNode(node);
    }

    node = ev = 0;
    compute(owner, 0, node, ev);
    for ( ; i < resultSize; ++i) {
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, node);
      getExpertForest(owner, 1)->setEdgeValue(result, i, ev);
    }
    getExpertForest(owner, 1)->unlinkNode(node);
  }

  // save result in compute cache and return it
  b = result;
  getExpertForest(owner, 1)->normalizeAndReduceNode(b, bev);
  saveResult(owner, a, b, bev);
}


// -------------------- Convert MDD to EV+MDD ------------------


mdd_to_evplusmdd_index_set* mdd_to_evplusmdd_index_set::getInstance()
{
  static mdd_to_evplusmdd_index_set instance;
  return &instance;
}


mdd_to_evplusmdd_index_set::mdd_to_evplusmdd_index_set()
{ }


mdd_to_evplusmdd_index_set::~mdd_to_evplusmdd_index_set() {}


compute_manager::error
mdd_to_evplusmdd_index_set::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 0)->isMdd())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 1)->isEvplusMdd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


compute_manager::error
mdd_to_evplusmdd_index_set::compute(op_info* owner,
    const dd_edge& a, dd_edge& b)
{
  if (owner == 0) return compute_manager::TYPE_MISMATCH;
  int result = 0;
  int ev = 0;
  int nVars = owner->p[0].getDomain()->getNumVariables();
  compute(owner, a.getNode(), nVars, result, ev);
  // note that ev will be equal to the cardinality of the ev+mdd;
  // but we do not use that number here.
  b.set(result, 0, getExpertForest(owner, 1)->getNodeLevel(result));
  return compute_manager::SUCCESS;
}


void mdd_to_evplusmdd_index_set::compute(op_info* owner,
    int a, int height, int& b, int& bev)
{
  DCASSERT(getExpertForest(owner, 1)->getTerminalNode(false) == 0);
  DCASSERT(getExpertForest(owner, 0)->getTerminalNode(true) == -1);
  DCASSERT(getExpertForest(owner, 0)->getTerminalNode(false) == 0);

  b = 0; bev = INF;

  // Deal with terminals
  if (a == 0) return;
  if (height == 0) {
    DCASSERT(getExpertForest(owner, 0)->getTerminalNode(true) == a);
    b = getExpertForest(owner, 1)->getTerminalNode(true);
    bev = 1;
    return;
  }

  expert_forest* mddf = getExpertForest(owner, 0);
  expert_forest* evmddf = getExpertForest(owner, 1);
  DCASSERT(mddf->getDomain() == evmddf->getDomain());
  expert_domain* d = smart_cast<expert_domain*>(mddf->useDomain());

  int aHeight = mddf->getNodeHeight(a);
  DCASSERT(aHeight <= height);

  // Search compute table
  if (aHeight == height) {
    if (findResult(owner, a, b, bev)) return;
  }

  // Create node at appropriate height
  const int resultLevel = d->getVariableWithHeight(height);
  const int resultSize = d->getVariableBound(resultLevel);
  int result = evmddf->createTempNode(resultLevel, resultSize, true);

  // Total counts all the paths leading from result
  int total = 0;

  if (aHeight < height) {
    int node = 0;
    int ev = INF;
    compute(owner, a, height - 1, node, ev);
    // ev gives the cardinality of node

    // Set the down-pointers and edge-values
    // Down-pointers are the same for all indexes
    // Edge-value at index i is ev*i
    if (node != 0) {
      for (int i = 0; i < resultSize; i++)
      {
        evmddf->setDownPtrWoUnlink(result, i, node);
        // edge-value indicates the number of paths from all the previous
        // down-pointers
        evmddf->setEdgeValue(result, i, total);
        // update the number of paths
        total += ev;
      }
      evmddf->unlinkNode(node);
    }
  } // aHeight < height
  else {
    // create node at level corresponding to variable 'aHeight'
    // note that aHeight == height
    
    if (mddf->isFullNode(a)) {
      // a is full (or truncated-full)
      int aSize = mddf->getFullNodeSize(a);
      for (int i = 0; i < aSize; i++)
      {
        int node = 0;
        int ev = INF;
        compute(owner, mddf->getFullNodeDownPtr(a, i), height - 1, node, ev);
        if (node == 0) continue;

        evmddf->setDownPtrWoUnlink(result, i, node);
        // edge-value indicates the number of paths from all the previous
        // down-pointers
        evmddf->setEdgeValue(result, i, total);
        evmddf->unlinkNode(node);
        // update the number of paths
        total += ev;
      }
    } // a is full
    else {
      // a is sparse
      int aNnz = mddf->getSparseNodeSize(a);
      for (int i = 0; i < aNnz; i++)
      {
        int node = 0;
        int ev = INF;
        compute(owner, mddf->getSparseNodeDownPtr(a, i), height - 1, node, ev);
        if (node == 0) continue;

        int index = mddf->getSparseNodeIndex(a, i);
        evmddf->setDownPtrWoUnlink(result, index, node);
        // edge-value indicates the number of paths from all the previous
        // down-pointers
        evmddf->setEdgeValue(result, index, total);
        evmddf->unlinkNode(node);
        // update the number of paths
        total += ev;
      }
    } // a is sparse
  } // aHeight == height

  // Save result in compute cache and return it
  b = result;
  evmddf->normalizeAndReduceNode(b, bev);
  DCASSERT((b == 0 && bev == INF) || (b != 0 && bev == 0));
  bev = total;
  saveResult(owner, a, b, bev);
  // Store the cardinality along with node
  evmddf->setIndexSetCardinality(b, bev);
}


// ----------------------- MTMXD Apply operation ----------------------


mtmxd_apply_operation::mtmxd_apply_operation()
{ }

mtmxd_apply_operation::~mtmxd_apply_operation() {}

compute_manager::error
mtmxd_apply_operation::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0] != owner->p[1] || owner->p[0] != owner->p[2])
    return compute_manager::FOREST_MISMATCH;
  if (!(owner->p[0].isMxd()) ||
      owner->p[0].isBoolForest() ||
      // the above allows MTMXDs with INTEGERs and REALs
      // the line below allows only INTEGERs
      // owner->f[0]->getRangeType() != forest::INTEGER ||
      !owner->p[0].isMT())
    return compute_manager::TYPE_MISMATCH;
  // mtmxd_apply_operation::compute() works correctly only with
  // operations that return 0 when operating on 0s. For example,
  // 0+0, 0*0, 0-0 all result in 0.
  // 0<=0, 0>=0, 0==0 all result in 1 (assuming 1 == true).
  int result = 0;
  if (checkTerminals(const_cast<op_info*>(owner), 0, 0, result)) {
    if (result != 0) {
      const_cast<expert_forest*>(getExpertForest(owner, 0))->unlinkNode(result);
      return compute_manager::TYPE_MISMATCH;
    }
  }
  return compute_manager::SUCCESS;
}


// ----------------------- MTMXD Multiply ----------------------------


mtmxd_multiply* mtmxd_multiply::getInstance()
{
  static mtmxd_multiply instance;
  return &instance;
}


mtmxd_multiply::mtmxd_multiply()
{ }


mtmxd_multiply::~mtmxd_multiply() {}


bool
mtmxd_multiply::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (a == 0 || b == 0) {
    c = 0;
    return true;
  }

  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getInteger(a) *
          getExpertForest(op, 0)->getInteger(b));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getReal(a) *
          getExpertForest(op, 0)->getReal(b));
    }
    return true;
  }
  return false;
}


// ----------------------- MTMXD Minus ----------------------------


mtmxd_minus* mtmxd_minus::getInstance()
{
  static mtmxd_minus instance;
  return &instance;
}


mtmxd_minus::mtmxd_minus()
{ }


mtmxd_minus::~mtmxd_minus() {}


bool
mtmxd_minus::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (a == b) {
    c = 0;
    return true;
  }

  // TODO: if (a == 0) c = negate(b)

  if (b == 0) {
    c = a;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getInteger(a) -
          getExpertForest(op, 0)->getInteger(b));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getReal(a) -
          getExpertForest(op, 0)->getReal(b));
    }
    return true;
  }
  return false;
}


// ----------------------- MTMXD Plus ----------------------------


mtmxd_plus* mtmxd_plus::getInstance()
{
  static mtmxd_plus instance;
  return &instance;
}


mtmxd_plus::mtmxd_plus()
{ }


mtmxd_plus::~mtmxd_plus() {}


bool
mtmxd_plus::checkTerminals(op_info* op, int a, int b, int& c)
{
#if 1
  if (a == 0 || b == 0) {
    c = a + b;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }
#endif

  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getInteger(a) +
          getExpertForest(op, 0)->getInteger(b));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getReal(a) +
          getExpertForest(op, 0)->getReal(b));
    }
    return true;
  }

  return false;
}


// ----------------------- MTMXD Min ----------------------------


mtmxd_min* mtmxd_min::getInstance()
{
  static mtmxd_min instance;
  return &instance;
}


mtmxd_min::mtmxd_min()
{ }


mtmxd_min::~mtmxd_min() {}


bool
mtmxd_min::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = (getExpertForest(op, 0)->getInteger(a) <
          getExpertForest(op, 0)->getInteger(b))? a: b;
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = (getExpertForest(op, 0)->getReal(a) <
          getExpertForest(op, 0)->getReal(b))? a: b;
    }
    return true;
  }
  return false;
}


// ----------------------- MTMXD Max ----------------------------


mtmxd_max* mtmxd_max::getInstance()
{
  static mtmxd_max instance;
  return &instance;
}


mtmxd_max::mtmxd_max()
{ }


mtmxd_max::~mtmxd_max() {}


bool
mtmxd_max::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = (getExpertForest(op, 0)->getInteger(a) >
          getExpertForest(op, 0)->getInteger(b))? a: b;
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = (getExpertForest(op, 0)->getReal(a) >
          getExpertForest(op, 0)->getReal(b))? a: b;
#if 0
      fprintf(stderr, "a: %d, b: %d, c: %d\n", a, b, c);
      fprintf(stderr, "a: %f, b: %f, c: %f\n",
          getExpertForest(op, 0)->getReal(a),
          getExpertForest(op, 0)->getReal(b),
          getExpertForest(op, 0)->getReal(c));
#endif
    }
    return true;
  }

  return false;
}


// ----------------------- MTMXD Less-Than ----------------------------


mtmxd_less_than* mtmxd_less_than::getInstance()
{
  static mtmxd_less_than instance;
  return &instance;
}


mtmxd_less_than::mtmxd_less_than()
{ }


mtmxd_less_than::~mtmxd_less_than() {}


bool
mtmxd_less_than::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode((int)
          (getExpertForest(op, 0)->getInteger(a) <
           getExpertForest(op, 0)->getInteger(b))? 1: 0);
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode((float)
          (getExpertForest(op, 0)->getReal(a) <
           getExpertForest(op, 0)->getReal(b))? 1.0f: 0.0f);
    }
    return true;
  }

  return false;
}


// ----------------------- MTMXD Greater-Than ----------------------------


mtmxd_greater_than* mtmxd_greater_than::getInstance()
{
  static mtmxd_greater_than instance;
  return &instance;
}


mtmxd_greater_than::mtmxd_greater_than()
{ }


mtmxd_greater_than::~mtmxd_greater_than() {}


bool
mtmxd_greater_than::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode((int)
          (getExpertForest(op, 0)->getInteger(a) >
           getExpertForest(op, 0)->getInteger(b))? 1: 0);
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode((float)
          (getExpertForest(op, 0)->getReal(a) >
           getExpertForest(op, 0)->getReal(b))? 1.0f: 0.0f);
    }
    return true;
  }

  return false;
}


// -------------------- MTMXD Not-Equal-To ----------------------


mtmxd_not_equal* mtmxd_not_equal::getInstance()
{
  static mtmxd_not_equal instance;
  return &instance;
}


mtmxd_not_equal::mtmxd_not_equal()
{ }


mtmxd_not_equal::~mtmxd_not_equal() {}


bool
mtmxd_not_equal::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(int((a != b)? 1: 0));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(float((a != b)? 1: 0));
    }
    return true;
  }
  return false;
}


// --------------------- MTMXD Alternate Apply operation --------------------


mtmxd_alt_apply_operation::
mtmxd_alt_apply_operation()
{ }

mtmxd_alt_apply_operation::~mtmxd_alt_apply_operation() {}

compute_manager::error
mtmxd_alt_apply_operation::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0] != owner->p[1] || owner->p[0] != owner->p[2])
    return compute_manager::FOREST_MISMATCH;
  if (!(owner->p[0].isMxd()) ||
      owner->p[0].isBoolForest() ||
      // the above allows MTMXDs with INTEGERs and REALs
      // the line below allows only INTEGERs
      // owner->f[0]->getRangeType() != forest::INTEGER ||
      !owner->p[0].isMT())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


// ----------------------- MTMXD Divide ----------------------------


mtmxd_divide* mtmxd_divide::getInstance()
{
  static mtmxd_divide instance;
  return &instance;
}


mtmxd_divide::mtmxd_divide()
{ }


mtmxd_divide::~mtmxd_divide() {}


bool
mtmxd_divide::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(
          int(getExpertForest(op, 0)->getInteger(a) /
            float(getExpertForest(op, 0)->getInteger(b))));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(
          getExpertForest(op, 0)->getReal(a) /
          getExpertForest(op, 0)->getReal(b));
    }
    return true;
  }
  return false;
}


// ---------------------- MTMXD Less-Than-Or-Equal-To ------------------------


mtmxd_less_than_equal* mtmxd_less_than_equal::getInstance()
{
  static mtmxd_less_than_equal instance;
  return &instance;
}


mtmxd_less_than_equal::mtmxd_less_than_equal()
{ }


mtmxd_less_than_equal::~mtmxd_less_than_equal() {}


bool
mtmxd_less_than_equal::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode((int)
          (getExpertForest(op, 0)->getInteger(a) <=
           getExpertForest(op, 0)->getInteger(b))? 1: 0);
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode((float)
          (getExpertForest(op, 0)->getReal(a) <=
           getExpertForest(op, 0)->getReal(b))? 1.0f: 0.0f);
    }
    return true;
  }

  return false;
}


// -------------------- MTMXD Greater-Than-Or-Equal-To ----------------------


mtmxd_greater_than_equal* mtmxd_greater_than_equal::getInstance()
{
  static mtmxd_greater_than_equal instance;
  return &instance;
}


mtmxd_greater_than_equal::mtmxd_greater_than_equal()
{ }


mtmxd_greater_than_equal::~mtmxd_greater_than_equal() {}


bool
mtmxd_greater_than_equal::checkTerminals(op_info* op,
    int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode((int)
          (getExpertForest(op, 0)->getInteger(a) >=
           getExpertForest(op, 0)->getInteger(b))? 1: 0);
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode((float)
          (getExpertForest(op, 0)->getReal(a) >=
           getExpertForest(op, 0)->getReal(b))? 1.0f: 0.0f);
    }
    return true;
  }
  return false;
}


// -------------------- MTMXD Equal-To ----------------------


mtmxd_equal* mtmxd_equal::getInstance()
{
  static mtmxd_equal instance;
  return &instance;
}


mtmxd_equal::mtmxd_equal()
{ }


mtmxd_equal::~mtmxd_equal() {}


bool
mtmxd_equal::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 0)->isTerminalNode(b)) {
    if (getExpertForest(op, 0)->getRangeType() == forest::INTEGER) {
      c = getExpertForest(op, 0)->getTerminalNode(int((a == b)? 1: 0));
    } else {
      DCASSERT(getExpertForest(op, 0)->getRangeType() == forest::REAL);
      c = getExpertForest(op, 0)->getTerminalNode(float((a == b)? 1: 0));
    }
    return true;
  }
  return false;
}





// ---------------------- EVMDD Operations -------------------------



evmdd_apply_operation::
evmdd_apply_operation()
{ }


evmdd_apply_operation::
~evmdd_apply_operation()
{ }


bool
evmdd_apply_operation::
isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nParams * 2
  // data[2i] and data[2i+1] <--> forest[i]
  // call isStale for each forest[i] and data[2i]
  DCASSERT(owner->nParams == 3);
  return 
    getExpertForest(owner, 0)->isStale(data[0]) ||
    getExpertForest(owner, 1)->isStale(data[2]) ||
    getExpertForest(owner, 2)->isStale(data[4]);
}


void
evmdd_apply_operation::
discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams * 2
  // data[2i] and data[2i+1] <--> forest[i]
  // call uncacheNode for each forest[i] and data[2i]
  DCASSERT(owner->nParams == 3);
  getExpertForest(owner, 0)->uncacheNode(data[0]);
  getExpertForest(owner, 1)->uncacheNode(data[2]);
  getExpertForest(owner, 2)->uncacheNode(data[4]);
}


// Calls compute(op_info*, dd_edge, dd_edge, dd_edge)
compute_manager::error
evmdd_apply_operation::
compute(op_info* owner, dd_edge** operands)
{
  if (operands == 0) return compute_manager::TYPE_MISMATCH;
  // compute(owner, dd_edge, dd_edge, dd_edge) checks for owner == 0
  return compute(owner, *operands[0], *operands[1], *operands[2]);
}


// Returns an error
compute_manager::error
evmdd_apply_operation::
compute(op_info* owner, const dd_edge& a, dd_edge& b)
{
  return compute_manager::TYPE_MISMATCH;
}



// ---------------------- EV+MDD Operations -------------------------



evplusmdd_apply_operation::
evplusmdd_apply_operation()
{ }


evplusmdd_apply_operation::
~evplusmdd_apply_operation()
{ }


void
evplusmdd_apply_operation::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  // data[] is of size owner.nParams * 2
  // data[2i] and data[2i+1] <--> forest[i]
  DCASSERT(owner->nParams == 3);
  fprintf(strm, "[%s(%d:%d, %d:%d): %d:%d]",
      owner->op->getName(),
      data[0], data[1], data[2], data[3], data[4], data[5]);
}


compute_manager::error
evplusmdd_apply_operation::
typeCheck(const op_info* owner)
{
  // op1 == EV+MDD, op2 == EV+MDD, op3 = EV+MDD
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0] != owner->p[1] || owner->p[0] != owner->p[2])
    return compute_manager::FOREST_MISMATCH;
  if (!getExpertForest(owner, 0)->isEvplusMdd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool
evplusmdd_apply_operation::
findResult(op_info* owner, int a, int aev, int b, int bev, int& c, int& cev)
{
  static int key[4];

  // create cache entry
  if (isCommutative() && a > b) {
    // sort the entry in ascending order
    key[0] = b; key[2] = a;
    key[1] = bev; key[3] = aev;
  } else {
    key[0] = a; key[2] = b;
    key[1] = aev; key[3] = bev;
  }

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  c = cacheEntry[4];
  cev = cacheEntry[5];
  getExpertForest(owner, 2)->linkNode(c);
  return true;
}


void
evplusmdd_apply_operation::
saveResult(op_info* owner, int a, int aev, int b, int bev, int c, int cev)
{
  static int cacheEntry[6];

  getExpertForest(owner, 0)->cacheNode(a);
  getExpertForest(owner, 1)->cacheNode(b);
  getExpertForest(owner, 2)->cacheNode(c);

  // create cache entry
  if (isCommutative() && a > b) {
    // sort the entry in ascending order
    cacheEntry[0] = b; cacheEntry[2] = a;
    cacheEntry[1] = bev; cacheEntry[3] = aev;
  } else {
    cacheEntry[0] = a; cacheEntry[2] = b;
    cacheEntry[1] = aev; cacheEntry[3] = bev;
  }
  cacheEntry[4] = c;
  cacheEntry[5] = cev;

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}


// Implements APPLY operation -- calls checkTerminals to compute
// result for terminal nodes.
compute_manager::error
evplusmdd_apply_operation::
compute(op_info* owner, const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  if (owner == 0) return compute_manager::TYPE_MISMATCH;
  int result = 0;
  int ev = 0;
  int aev = 0;
  int bev = 0;
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(owner, a.getNode(), aev, b.getNode(), bev, result, ev);
  c.set(result, ev, getExpertForest(owner, 1)->getNodeLevel(result));
  return compute_manager::SUCCESS;
}


compute_manager::error
evplusmdd_apply_operation::
compute(op_info* owner, int a, int aev, int b, int bev, int& c, int& cev)
{
  DCASSERT(owner->p[0] == owner->p[1] && owner->p[1] == owner->p[2]);

  if (checkTerminals(owner, a, aev, b, bev, c, cev))
    return compute_manager::SUCCESS;
  if (findResult(owner, a, aev, b, bev, c, cev))
    return compute_manager::SUCCESS;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  expert_forest* expertForest = getExpertForest(owner, 0);
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = expertForest->getLevelSize(resultLevel);

  // Three vectors: operands a and b, and result c

  if (aLevel < resultLevel) {
    // expand b
    // result[i] = a op b[i]
    std::vector<int> B(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<int> Bev(resultSize, INF);
    std::vector<int> Cev(resultSize, INF);

    int aZero, aZeroev;
    compute(owner, a, aev, 0, INF, aZero, aZeroev);

    expertForest->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterBev = Bev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    for ( ; iterB != B.end(); iterB++, iterC++, iterBev++, iterCev++)
    {
      if (*iterB == 0) {
        *iterC = aZero;
        expertForest->linkNode(aZero);
        *iterCev = aZeroev;
      } else {
        compute(owner, a, aev, *iterB, *iterBev + bev, *iterC, *iterCev);
      }
    }
    expertForest->unlinkNode(aZero);
    c = expertForest->createTempNode(resultLevel, C, Cev);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    std::vector<int> A(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<int> Aev(resultSize, INF);
    std::vector<int> Cev(resultSize, INF);

    int zeroB, zeroBev;
    compute(owner, 0, INF, b, bev, zeroB, zeroBev);

    expertForest->getDownPtrsAndEdgeValues(a, A, Aev);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterAev = Aev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end(); iterA++, iterC++, iterAev++, iterCev++)
    {
      if (*iterA == 0) {
        *iterC = zeroB;
        expertForest->linkNode(zeroB);
        *iterCev = zeroBev;
      } else {
        compute(owner, *iterA, *iterAev + aev, b, bev, *iterC, *iterCev);
      }
    }
    expertForest->unlinkNode(zeroB);
    c = expertForest->createTempNode(resultLevel, C, Cev);
  }
  else {
    // expand both a and b
    // result[i] = a[i] op b[i]
    std::vector<int> A(resultSize, 0);
    std::vector<int> B(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<int> Aev(resultSize, INF);
    std::vector<int> Bev(resultSize, INF);
    std::vector<int> Cev(resultSize, INF);

    int zeroZero, zeroZeroEv;
    compute(owner, 0, INF, 0, INF, zeroZero, zeroZeroEv);

    expertForest->getDownPtrsAndEdgeValues(a, A, Aev);
    expertForest->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterAev = Aev.begin();
    std::vector<int>::iterator iterBev = Bev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end();
        iterA++, iterB++, iterC++, iterAev++, iterBev++, iterCev++)
    {
      if (*iterA == 0) {
        if (*iterB == 0) {
          *iterC = zeroZero;
          expertForest->linkNode(zeroZero);
          *iterCev = zeroZeroEv;
        } else {
          compute(owner, 0, INF, *iterB, *iterBev + bev, *iterC, *iterCev);
        }
      } else if (*iterB == 0) {
        compute(owner, *iterA, *iterAev + aev, 0, INF, *iterC, *iterCev);
      } else {
        compute(owner, *iterA, *iterAev + aev, *iterB, *iterBev + bev,
            *iterC, *iterCev);
      }
    }

    expertForest->unlinkNode(zeroZero);
    c = expertForest->createTempNode(resultLevel, C, Cev);
  }

  // save result in compute cache and return it

#if 0
  printf("reduce(%d): ", result);
  result = expertForest->reduceNode(result);
  printf("%d  [", result);
  for (unsigned i = 0; i < C.size(); i++ )
  {
    printf("%d ", C[i]);
  }
  printf("]\n");
#else
  expertForest->normalizeAndReduceNode(c, cev);
#endif

  saveResult(owner, a, aev, b, bev, c, cev);
  return compute_manager::SUCCESS;
}


// ---------------------- EV*MDD Operations -------------------------



evtimesmdd_apply_operation::
evtimesmdd_apply_operation()
{ }


evtimesmdd_apply_operation::
~evtimesmdd_apply_operation()
{ }


void
evtimesmdd_apply_operation::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  // data[] is of size owner.nParams * 2
  // data[2i] and data[2i+1] <--> forest[i]
  DCASSERT(owner->nParams == 3);
  fprintf(strm, "[%s(%d:%f, %d:%f): %d:%f]",
      owner->op->getName(),
      data[0], toFloat(data[1]), data[2], toFloat(data[3]),
      data[4], toFloat(data[5]));
}


compute_manager::error
evtimesmdd_apply_operation::
typeCheck(const op_info* owner)
{
  // op1 == EV+MDD, op2 == EV+MDD, op3 = EV+MDD
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (!owner->areAllForests())
    return compute_manager::TYPE_MISMATCH;
  if (owner->nParams != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->p[0] != owner->p[1] || owner->p[0] != owner->p[2])
    return compute_manager::FOREST_MISMATCH;
  if (!getExpertForest(owner, 0)->isEvtimesMdd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool
evtimesmdd_apply_operation::
findResult(op_info* owner, int a, float aev, int b, float bev,
    int& c, float& cev)
{
  static int key[4];

  // create cache entry
  if (isCommutative() && a > b) {
    // sort the entry in ascending order
    key[0] = b; key[2] = a;
    key[1] = toInt(bev); key[3] = toInt(aev);
  } else {
    key[0] = a; key[2] = b;
    key[1] = toInt(aev); key[3] = toInt(bev);
  }

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  c = cacheEntry[4];
  cev = toFloat(cacheEntry[5]);
  getExpertForest(owner, 2)->linkNode(c);
  return true;
}


void
evtimesmdd_apply_operation::
saveResult(op_info* owner, int a, float aev, int b, float bev, int c, float cev)
{
  static int cacheEntry[6];

  getExpertForest(owner, 0)->cacheNode(a);
  getExpertForest(owner, 1)->cacheNode(b);
  getExpertForest(owner, 2)->cacheNode(c);

  // create cache entry
  if (isCommutative() && a > b) {
    // sort the entry in ascending order
    cacheEntry[0] = b; cacheEntry[2] = a;
    cacheEntry[1] = toInt(bev); cacheEntry[3] = toInt(aev);
  } else {
    cacheEntry[0] = a; cacheEntry[2] = b;
    cacheEntry[1] = toInt(aev); cacheEntry[3] = toInt(bev);
  }
  cacheEntry[4] = c;
  cacheEntry[5] = toInt(cev);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}


// Implements APPLY operation -- calls checkTerminals to compute
// result for terminal nodes.
compute_manager::error
evtimesmdd_apply_operation::
compute(op_info* owner, const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  if (owner == 0) return compute_manager::TYPE_MISMATCH;
  int result = 0;
  float ev = 0;
  float aev = 0;
  float bev = 0;
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(owner, a.getNode(), aev, b.getNode(), bev, result, ev);
  c.set(result, ev, getExpertForest(owner, 1)->getNodeLevel(result));
  return compute_manager::SUCCESS;
}


compute_manager::error
evtimesmdd_apply_operation::
compute(op_info* owner, int a, float aev, int b, float bev, int& c, float& cev)
{
  DCASSERT(owner->p[0] == owner->p[1] && owner->p[1] == owner->p[2]);

  if (checkTerminals(owner, a, aev, b, bev, c, cev))
    return compute_manager::SUCCESS;
  if (findResult(owner, a, aev, b, bev, c, cev))
    return compute_manager::SUCCESS;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  expert_forest* expertForest = getExpertForest(owner, 0);
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = expertForest->getLevelSize(resultLevel);

  // Three vectors: operands a and b, and result c

  if (aLevel < resultLevel) {
    // expand b
    // result[i] = a op b[i]
    std::vector<int> B(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<float> Bev(resultSize, NAN);
    std::vector<float> Cev(resultSize, NAN);

    int aZero;
    float aZeroev;
    compute(owner, a, aev, 0, NAN, aZero, aZeroev);

    expertForest->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<float>::iterator iterBev = Bev.begin();
    std::vector<float>::iterator iterCev = Cev.begin();
    for ( ; iterB != B.end(); iterB++, iterC++, iterBev++, iterCev++)
    {
      if (*iterB == 0) {
        *iterC = aZero;
        expertForest->linkNode(aZero);
        *iterCev = aZeroev;
      } else {
        compute(owner, a, aev, *iterB, *iterBev * bev, *iterC, *iterCev);
      }
    }
    expertForest->unlinkNode(aZero);
    c = expertForest->createTempNode(resultLevel, C, Cev);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    std::vector<int> A(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<float> Aev(resultSize, NAN);
    std::vector<float> Cev(resultSize, NAN);

    int zeroB;
    float zeroBev;
    compute(owner, 0, NAN, b, bev, zeroB, zeroBev);

    expertForest->getDownPtrsAndEdgeValues(a, A, Aev);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<float>::iterator iterAev = Aev.begin();
    std::vector<float>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end(); iterA++, iterC++, iterAev++, iterCev++)
    {
      if (*iterA == 0) {
        *iterC = zeroB;
        expertForest->linkNode(zeroB);
        *iterCev = zeroBev;
      } else {
        compute(owner, *iterA, *iterAev * aev, b, bev, *iterC, *iterCev);
      }
    }
    expertForest->unlinkNode(zeroB);
    c = expertForest->createTempNode(resultLevel, C, Cev);
  }
  else {
    // expand both a and b
    // result[i] = a[i] op b[i]
    std::vector<int> A(resultSize, 0);
    std::vector<int> B(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<float> Aev(resultSize, NAN);
    std::vector<float> Bev(resultSize, NAN);
    std::vector<float> Cev(resultSize, NAN);

    int zeroZero;
    float zeroZeroEv;
    compute(owner, 0, NAN, 0, NAN, zeroZero, zeroZeroEv);

    expertForest->getDownPtrsAndEdgeValues(a, A, Aev);
    expertForest->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<float>::iterator iterAev = Aev.begin();
    std::vector<float>::iterator iterBev = Bev.begin();
    std::vector<float>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end();
        iterA++, iterB++, iterC++, iterAev++, iterBev++, iterCev++)
    {
      if (*iterA == 0) {
        if (*iterB == 0) {
          *iterC = zeroZero;
          expertForest->linkNode(zeroZero);
          *iterCev = zeroZeroEv;
        } else {
          compute(owner, 0, NAN, *iterB, *iterBev * bev, *iterC, *iterCev);
        }
      } else if (*iterB == 0) {
        compute(owner, *iterA, *iterAev * aev, 0, NAN, *iterC, *iterCev);
      } else {
        compute(owner, *iterA, *iterAev * aev, *iterB, *iterBev * bev,
            *iterC, *iterCev);
      }
    }

    expertForest->unlinkNode(zeroZero);
    c = expertForest->createTempNode(resultLevel, C, Cev);
  }

  // save result in compute cache and return it

#if 0
  int size = resultSize;
  float tempev = 0;
  printf("reduce(%d):\n", c);
  printf("  before %d:%f: [", c, cev);
  for (int i = 0; i < size; i++ ) {
    expertForest->getFullNodeEdgeValue(c, i, tempev);
    printf("%d:%0.6f ",
        expertForest->getFullNodeDownPtr(c, i), tempev);
  }
  printf("]\n");

  expertForest->normalizeAndReduceNode(c, cev);

  std::vector<int> C(resultSize, 0);
  std::vector<float> Cev(resultSize, NAN);
  expertForest->getDownPtrsAndEdgeValues(c, C, Cev);
  printf("  after %d:%f: [", c, cev);
  for (unsigned i = 0; i < C.size(); i++ ) {
    printf("%d:%0.6f ", C[i], Cev[i]);
  }
  printf("]\n\n");
#else
  expertForest->normalizeAndReduceNode(c, cev);
#endif

  saveResult(owner, a, aev, b, bev, c, cev);
  return compute_manager::SUCCESS;
}


evplusmdd_plus*
evplusmdd_plus::
getInstance()
{
  static evplusmdd_plus instance;
  return &instance;
}


bool
evplusmdd_plus::
checkTerminals(op_info* op, int a, int aev, int b, int bev, int& c, int& cev)
{
  if (a == 0) {
    c = b; cev = bev;
    getExpertForest(op, 1)->linkNode(b);
    return true;
  }
  if (b == 0) {
    c = a; cev = aev;
    getExpertForest(op, 0)->linkNode(a);
    return true;
  }
  if (a == -1 && b == -1) {
    c = -1; cev = aev + bev;
    return true;
  }
  return false;
}


#if 1

compute_manager::error
evplusmdd_plus::
compute(op_info* owner, int a, int aev, int b, int bev, int& c, int& cev)
{
#if 1
  DCASSERT(owner->p[0] == owner->p[1] && owner->p[1] == owner->p[2]);

  if (checkTerminals(owner, a, aev, b, bev, c, cev))
    return compute_manager::SUCCESS;
  if (findResult(owner, a, aev, b, bev, c, cev))
    return compute_manager::SUCCESS;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  expert_forest* expertForest = getExpertForest(owner, 0);
  const int aLevel = expertForest->getNodeLevel(a);
  const int bLevel = expertForest->getNodeLevel(b);
  int resultLevel = aLevel > bLevel? aLevel: bLevel;

  // Three vectors: operands a and b, and result c

  if (aLevel < resultLevel) {
    // expand b
    // result[i] = a op b[i]
    int aZero, aZeroev;
    compute(owner, a, aev, 0, INF, aZero, aZeroev);

    std::vector<int> B, Bev;
    expertForest->getDownPtrsAndEdgeValues(b, B, Bev);

    std::vector<int> C(B.size(), 0);
    std::vector<int> Cev(B.size(), INF);

    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterBev = Bev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    for ( ; iterB != B.end(); iterB++, iterC++, iterBev++, iterCev++)
    {
      if (*iterB == 0) {
        *iterC = aZero;
        expertForest->linkNode(aZero);
        *iterCev = aZeroev;
      } else {
        compute(owner, a, aev, *iterB, *iterBev + bev, *iterC, *iterCev);
      }
    }
    expertForest->unlinkNode(aZero);
    c = expertForest->createTempNode(resultLevel, C, Cev);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    int zeroB, zeroBev;
    compute(owner, 0, INF, b, bev, zeroB, zeroBev);

    std::vector<int> A, Aev;
    expertForest->getDownPtrsAndEdgeValues(a, A, Aev);

    std::vector<int> C(A.size(), 0);
    std::vector<int> Cev(A.size(), INF);

    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterAev = Aev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end(); iterA++, iterC++, iterAev++, iterCev++)
    {
      if (*iterA == 0) {
        *iterC = zeroB;
        expertForest->linkNode(zeroB);
        *iterCev = zeroBev;
      } else {
        compute(owner, *iterA, *iterAev + aev, b, bev, *iterC, *iterCev);
      }
    }
    expertForest->unlinkNode(zeroB);
    c = expertForest->createTempNode(resultLevel, C, Cev);
  }
  else {
    // expand both a and b
    // result[i] = a[i] op b[i]
    int zeroZero, zeroZeroEv;
    compute(owner, 0, INF, 0, INF, zeroZero, zeroZeroEv);

    std::vector<int> A, Aev, B, Bev;
    expertForest->getDownPtrsAndEdgeValues(a, A, Aev);
    expertForest->getDownPtrsAndEdgeValues(b, B, Bev);

    int resultSize = MAX(A.size(), B.size());
    std::vector<int> C(resultSize, 0);
    std::vector<int> Cev(resultSize, INF);

    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterAev = Aev.begin();
    std::vector<int>::iterator iterBev = Bev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();

    for ( ; iterA != A.end() && iterB != B.end();
        iterA++, iterB++, iterC++, iterAev++, iterBev++, iterCev++)
    {
      if (*iterA == 0) {
        if (*iterB == 0) {
          *iterC = zeroZero;
          expertForest->linkNode(zeroZero);
          *iterCev = zeroZeroEv;
        } else {
          compute(owner, 0, INF, *iterB, *iterBev + bev, *iterC, *iterCev);
        }
      } else if (*iterB == 0) {
        compute(owner, *iterA, *iterAev + aev, 0, INF, *iterC, *iterCev);
      } else {
        compute(owner, *iterA, *iterAev + aev, *iterB, *iterBev + bev,
            *iterC, *iterCev);
      }
    }
    for ( ; iterA != A.end(); iterA++, iterC++, iterAev++, iterCev++)
    {
      // C[i] = A[i], Cev[i] = Aev[i] + aev
      if (*iterA == 0) {
        *iterC = zeroZero;
        expertForest->linkNode(zeroZero);
        *iterCev = zeroZeroEv;
      } else {
        compute(owner, *iterA, *iterAev + aev, 0, INF, *iterC, *iterCev);
      }
    }
    for ( ; iterB != B.end(); iterB++, iterC++, iterBev++, iterCev++)
    {
      // C[i] = B[i], Cev[i] = Bev[i] + bev
      if (*iterB == 0) {
        *iterC = zeroZero;
        expertForest->linkNode(zeroZero);
        *iterCev = zeroZeroEv;
      } else {
        compute(owner, 0, INF, *iterB, *iterBev + bev, *iterC, *iterCev);
      }
    }


    expertForest->unlinkNode(zeroZero);
    c = expertForest->createTempNode(resultLevel, C, Cev);
  }

  // save result in compute cache and return it

#if 0
  printf("reduce(%d): ", result);
  result = expertForest->reduceNode(result);
  printf("%d  [", result);
  for (unsigned i = 0; i < C.size(); i++ )
  {
    printf("%d ", C[i]);
  }
  printf("]\n");
#else
  expertForest->normalizeAndReduceNode(c, cev);
#endif

  saveResult(owner, a, aev, b, bev, c, cev);
  return compute_manager::SUCCESS;
#else
  return evplusmdd_apply_operation::compute(owner, a, aev, b, bev, c, cev);
#endif
}

#else

compute_manager::error
evplusmdd_plus::
compute(op_info* owner, int a, int aev, int b, int bev, int& c, int& cev)
{
  DCASSERT(owner->p[0] == owner->p[1] && owner->p[1] == owner->p[2]);

  if (checkTerminals(owner, a, aev, b, bev, c, cev))
    return compute_manager::SUCCESS;
  if (findResult(owner, a, aev, b, bev, c, cev))
    return compute_manager::SUCCESS;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  expert_forest* expertForest = getExpertForest(owner, 0);
  int aLevel = expertForest->getNodeLevel(a);
  int bLevel = expertForest->getNodeLevel(b);
  int resultLevel = aLevel > bLevel? aLevel: bLevel;

  // Three vectors: operands a and b, and result c

  if (aLevel < resultLevel || bLevel < resultLevel) {
    if (bLevel < resultLevel) {
      // commutative operation; so swap and solve as (aLevel < resultLevel)
      SWAP(a, b);
      SWAP(aev, bev);
      SWAP(aLevel, bLevel);
    }
    DCASSERT(aLevel < resultLevel);

#if 1
    std::vector<int> B, Bev;
    expertForest->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int> C(expertForest->getLevelSize(resultLevel), a);
    // note: linknode is to be done later on
    std::vector<int> Cev(C.size(), aev);

    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterBev = Bev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    int noChangeCount = 0;
    for ( ; iterB != B.end(); iterB++, iterC++, iterBev++, iterCev++)
    {
      if (*iterB == 0) {
        noChangeCount++;
      } else {
        compute(owner, a, aev, *iterB, *iterBev + bev, *iterC, *iterCev);
      }
    }
    // update in-count for the rest of the nodes in C[]
    expertForest->getInCount(a) += (noChangeCount + C.size() - B.size());
    c = expertForest->createTempNode(resultLevel, C, Cev);
#else
    std::vector<int> C(expertForest->getLevelSize(resultLevel), a);
    // note: linknode is to be done later on
    std::vector<int> Cev(C.size(), aev);

    if (expertForest->isFullNode(b)) {
      int sz = expertForest->getFullNodeSize(b);
      int nnz = 0;
      for (int i = 0; i < sz; i++)
      {
        int dptr = expertForest->getFullNodeDownPtr(b, i);
        if (dptr != 0) {
          nnz++;
          int ev = INF;
          expertForest->getFullNodeEdgeValue(b, i, ev);
          compute(owner, a, aev, dptr, ev + bev, C[i], Cev[i]);
        }
      }
      expertForest->getInCount(a) += (C.size() - nnz);
    }
    else {
      DCASSERT(expertForest->isSparseNode(b));
      int nnz = expertForest->getSparseNodeSize(b);
      for (int i = 0; i < nnz; i++)
      {
        int dptr = expertForest->getSparseNodeDownPtr(b, i);
        int index = expertForest->getSparseNodeIndex(b, i);
        int ev = INF;
        expertForest->getSparseNodeEdgeValue(b, i, ev);
        compute(owner, a, aev, dptr, ev + bev, C[index], Cev[index]);
      }
      expertForest->getInCount(a) += (C.size() - nnz);
    }
    c = expertForest->createTempNode(resultLevel, C, Cev);
#endif
  }
  else {
    // expand both a and b
    // result[i] = a[i] op b[i]
#if 1
    std::vector<int> A, Aev, B, Bev;
    expertForest->getDownPtrsAndEdgeValues(a, A, Aev);
    expertForest->getDownPtrsAndEdgeValues(b, B, Bev);

    std::vector<int> *X, *Xev, *Y, *Yev;
    int x, y, xev, yev;
    if (B.size() > A.size()) {
      X = &B; Xev = &Bev; x = b; xev = bev;
      Y = &A; Yev = &Aev; y = a; yev = aev;
    } else {
      X = &A; Xev = &Aev; x = a; xev = aev;
      Y = &B; Yev = &Bev; y = b; yev = bev;
    }

    // result in X. Add Y to it.
    std::vector<int>::iterator iterX = (*X).begin();
    std::vector<int>::iterator iterY = (*Y).begin();
    std::vector<int>::iterator iterXev = (*Xev).begin();
    std::vector<int>::iterator iterYev = (*Yev).begin();

    for ( ; iterY != (*Y).end(); iterY++, iterX++, iterYev++, iterXev++)
    {
      if (*iterX == 0) {
        // if *iterY == 0, do nothing as *iterX is already 0
        if (*iterY != 0) {
          *iterX = *iterY;
          *iterXev = *iterYev + yev;
          expertForest->linkNode(*iterX);
        }
      } else {
        if (*iterY == 0) {
          *iterXev += xev;
          expertForest->linkNode(*iterX);
        } else {
          compute(owner, *iterX, *iterXev + xev, *iterY, *iterYev + yev,
              *iterX, *iterXev);
        }
      }
    }
    for ( ; iterX != (*X).end(); iterX++)
    {
      if (*iterX != 0) {
        *iterXev += xev;
        expertForest->linkNode(*iterX);
      }
    }

    c = expertForest->createTempNode(resultLevel, *X, *Xev);
#else
    std::vector<int> A, Aev;
    expertForest->getDownPtrsAndEdgeValues(a, A, Aev);

    if (expertForest->isFullNode(b)) {
      int sz = expertForest->getFullNodeSize(b);
      if (int(A.size()) < sz) { A.resize(sz, 0); Aev.resize(sz, INF); }
      for (int i = 0; i < sz; i++)
      {
        int dptr = expertForest->getFullNodeDownPtr(b, i);
        if (dptr == 0) {
          if (A[i] != 0) {
            expertForest->linkNode(A[i]);
            Aev[i] += aev;
          }
        }
        else {
          int ev = INF;
          expertForest->getFullNodeEdgeValue(b, i, ev);
          compute(owner, A[i], Aev[i] + aev, dptr, ev + bev, A[i], Aev[i]);
        }
      }
      for (int i = int(A.size()) - 1; i >= sz; i--)
      {
        if (A[i] != 0) {
          expertForest->linkNode(A[i]);
          Aev[i] += aev;
        }
      }
    }
    else {
      DCASSERT(expertForest->isSparseNode(b));
      for (int i = int(A.size()) - 1; i >= 0; i--)
      {
        if (A[i] != 0) {
          expertForest->linkNode(A[i]);
          Aev[i] += aev;
        }
      }
      int nnz = expertForest->getSparseNodeSize(b);
      int sz = expertForest->getSparseNodeIndex(b, nnz - 1) + 1;
      if (int(A.size()) < sz) { A.resize(sz, 0); Aev.resize(sz, INF); }
      for (int i = 0; i < nnz; i++)
      {
        int dptr = expertForest->getSparseNodeDownPtr(b, i);
        int ev = INF;
        expertForest->getSparseNodeEdgeValue(b, i, ev);
        int index = expertForest->getSparseNodeIndex(b, i);
        int temp = A[index];
        compute(owner, A[index], Aev[index], dptr, ev + bev,
          A[index], Aev[index]);
        if (temp != 0) expertForest->unlinkNode(temp);
      }
    }
    c = expertForest->createTempNode(resultLevel, A, Aev);
#endif
  }

  // save result in compute cache and return it

#if 0
  printf("reduce(%d): ", result);
  result = expertForest->reduceNode(result);
  printf("%d  [", result);
  for (unsigned i = 0; i < C.size(); i++ )
  {
    printf("%d ", C[i]);
  }
  printf("]\n");
#else
  expertForest->normalizeAndReduceNode(c, cev);
#endif

  saveResult(owner, a, aev, b, bev, c, cev);
  return compute_manager::SUCCESS;
}

#endif


evplusmdd_multiply*
evplusmdd_multiply::
getInstance()
{
  static evplusmdd_multiply instance;
  return &instance;
}


bool
evplusmdd_multiply::
checkTerminals(op_info* op, int a, int aev, int b, int bev, int& c, int& cev)
{
  if (a == 0 || b == 0) {
    c = 0; cev = INF;
    return true;
  }
  if (a == -1 && b == -1) {
    c = -1; cev = aev * bev;
    return true;
  }
  return false;
}


evplusmdd_minus*
evplusmdd_minus::
getInstance()
{
  static evplusmdd_minus instance;
  return &instance;
}


bool
evplusmdd_minus::
checkTerminals(op_info* op, int a, int aev, int b, int bev, int& c, int& cev)
{
  if (a == 0) {
    c = b; cev = -bev;
    getExpertForest(op, 1)->linkNode(b);
    return true;
  }
  if (b == 0) {
    c = a; cev = aev;
    getExpertForest(op, 0)->linkNode(a);
    return true;
  }
  if (a == -1 && b == -1) {
    c = -1; cev = aev - bev;
    return true;
  }
  return false;
}


evtimesmdd_plus*
evtimesmdd_plus::
getInstance()
{
  static evtimesmdd_plus instance;
  return &instance;
}


bool
evtimesmdd_plus::
checkTerminals(op_info* op, int a, float aev, int b, float bev,
    int& c, float& cev)
{
  if (a == 0) {
    c = b; cev = bev;
    getExpertForest(op, 1)->linkNode(b);
    return true;
  }
  if (b == 0) {
    c = a; cev = aev;
    getExpertForest(op, 0)->linkNode(a);
    return true;
  }
  if (a == -1 && b == -1) {
    c = -1; cev = aev + bev;
    return true;
  }
  return false;
}


evtimesmdd_multiply*
evtimesmdd_multiply::
getInstance()
{
  static evtimesmdd_multiply instance;
  return &instance;
}


bool
evtimesmdd_multiply::
checkTerminals(op_info* op, int a, float aev, int b, float bev,
    int& c, float& cev)
{
  if (a == 0 || b == 0) {
    c = 0; cev = NAN;
    return true;
  }
  if (a == -1 && b == -1) {
    c = -1; cev = aev * bev;
    return true;
  }
  return false;
}


evtimesmdd_minus*
evtimesmdd_minus::
getInstance()
{
  static evtimesmdd_minus instance;
  return &instance;
}


bool
evtimesmdd_minus::
checkTerminals(op_info* op, int a, float aev, int b, float bev,
    int& c, float& cev)
{
  if (a == 0) {
    c = b; cev = -bev;
    getExpertForest(op, 1)->linkNode(b);
    return true;
  }
  if (b == 0) {
    c = a; cev = aev;
    getExpertForest(op, 0)->linkNode(a);
    return true;
  }
  if (a == -1 && b == -1) {
    c = -1; cev = aev - bev;
    return true;
  }
  return false;
}


evtimesmdd_equal*
evtimesmdd_equal::
getInstance()
{
  static evtimesmdd_equal instance;
  return &instance;
}


bool
evtimesmdd_equal::
checkTerminals(op_info* op, int a, float aev, int b, float bev,
    int& c, float& cev)
{
  if (getExpertForest(op, 0)->isTerminalNode(a) &&
      getExpertForest(op, 1)->isTerminalNode(b)) {
    if (a == b && ((aev == bev) || (isNan(aev) && isNan(bev)))) {
      c = -1;
      cev = 1.0;
    } else {
      c = 0;
      cev = NAN;
    }
    return true;
  }

  return false;
}



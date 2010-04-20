
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
//#define DEBUG_DFS

inline expert_forest* getExpertForest(op_info* op, int index) {
  return smart_cast<expert_forest*>(op->f[index]);
}

inline const expert_forest* getExpertForest(const op_info* op, int index) {
  return smart_cast<const expert_forest*>(op->f[index]);
}

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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] != owner->f[1] || owner->f[0] != owner->f[2])
    return compute_manager::FOREST_MISMATCH;
  if (!getExpertForest(owner, 0)->isMdd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mdd_apply_operation::
isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nForests
  // data[i] <--> forest[i]
  // call isStale for each forest[i] and data[i]
  DCASSERT(owner->nForests == 3);
  return 
    getExpertForest(owner, 0)->isStale(data[0]) ||
    getExpertForest(owner, 1)->isStale(data[1]) ||
    getExpertForest(owner, 2)->isStale(data[2]);
}


void
mdd_apply_operation::
discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nForests
  // data[i] <--> forest[i]
  // call uncacheNode for each forest[i] and data[i]
  DCASSERT(owner->nForests == 3);
  getExpertForest(owner, 0)->uncacheNode(data[0]);
  getExpertForest(owner, 1)->uncacheNode(data[1]);
  getExpertForest(owner, 2)->uncacheNode(data[2]);
}


void
mdd_apply_operation::
showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // data[] is of size owner.nForests
  // data[i] <--> forest[i]
  // call showNode for each forest[i] and data[i]
  DCASSERT(owner->nForests == 3);
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
  expert_forest* expertForest = smart_cast<expert_forest*>(owner->f[0]);
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
  expert_forest* expertForest = smart_cast<expert_forest*>(owner->f[0]);
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] != owner->f[1] || owner->f[0] != owner->f[2])
    return compute_manager::FOREST_MISMATCH;
  if (!(owner->f[0]->isForRelations()) ||
      owner->f[0]->getRangeType() != forest::BOOLEAN ||
      owner->f[0]->getEdgeLabeling() != forest::MULTI_TERMINAL)
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


int mxd_apply_operation::compute(op_info* owner, int a, int b) {
  return smart_cast<expert_forest*>(owner->f[0])->getReductionRule() ==
      forest::IDENTITY_REDUCED
      ? computeIdent(owner, a, b)
      : computeNonIdent(owner, a, b);
}


int mxd_apply_operation::computeNonIdent(op_info* owner, int a, int b)
{
  expert_forest* expertForest = smart_cast<expert_forest*>(owner->f[0]);

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
  expert_forest* expertForest = smart_cast<expert_forest*>(owner->f[0]);
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] != owner->f[1] || owner->f[0] != owner->f[2])
    return compute_manager::FOREST_MISMATCH;
  if (!(owner->f[0]->isForRelations()) ||
      owner->f[0]->getRangeType() != forest::BOOLEAN ||
      owner->f[0]->getEdgeLabeling() != forest::MULTI_TERMINAL)
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mxd_alt_apply_operation::isEntryStale(const op_info* owner,
    const int* data)
{
  // data[] is of size owner.nForests + 1
  // data[0] is level
  // data[i] <--> forest[i-1]
  DCASSERT(owner->nForests == 3);
  return 
    getExpertForest(owner, 0)->isStale(data[1]) ||
    getExpertForest(owner, 1)->isStale(data[2]) ||
    getExpertForest(owner, 2)->isStale(data[3]);
}


void
mxd_alt_apply_operation::discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nForests + 1
  // data[0] is level
  // data[i] <--> forest[i-1]
  DCASSERT(owner->nForests == 3);
  getExpertForest(owner, 0)->uncacheNode(data[1]);
  getExpertForest(owner, 1)->uncacheNode(data[2]);
  getExpertForest(owner, 2)->uncacheNode(data[3]);
}


void
mxd_alt_apply_operation::showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // data[] is of size owner.nForests + 1
  // data[0] is level
  // data[i] <--> forest[i-1]
  DCASSERT(owner->nForests == 3);
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
      owner->f[0]->getDomain()->getTopVariable(), a.getNode(), b.getNode());
  c.set(result, 0, getExpertForest(owner, 2)->getNodeLevel(result));
  return compute_manager::SUCCESS;
}


int mxd_alt_apply_operation::compute(op_info* owner, int resultLevel,
    int a, int b)
{
  return smart_cast<expert_forest*>(owner->f[0])->getReductionRule() ==
      forest::IDENTITY_REDUCED
      ? computeIdent(owner, resultLevel, a, b)
      : computeNonIdent(owner, resultLevel, a, b);
}


int mxd_alt_apply_operation::computeNonIdent(op_info* owner,
    int resultLevel, int a, int b)
{
  expert_forest* expertForest = smart_cast<expert_forest*>(owner->f[0]);

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
  expert_forest* expertForest = smart_cast<expert_forest*>(owner->f[0]);

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


// ---------------------- MDD MXD Image Operation -------------------


mdd_mxd_image_operation::mdd_mxd_image_operation()
{ }

mdd_mxd_image_operation::~mdd_mxd_image_operation() {}


#if 0
bool
mdd_mxd_image_operation::findResult(op_info* owner, int a, int b, int& c)
{
  static int key[2];
  // create cache entry
  key[0] = a; key[1] = b;

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  c = cacheEntry[2];
  getExpertForest(owner, 2)->linkNode(c);
  return true;
}


void
mdd_mxd_image_operation::saveResult(op_info* owner, int a, int b, int c)
{
  static int cacheEntry[3];

  // create cache entry
  cacheEntry[0] = a; cacheEntry[1] = b; cacheEntry[2] = c;

#if 0
#ifdef DEVELOPMENT_CODE
  assert(!findResult(owner, cacheEntry[0], cacheEntry[1], cacheEntry[2]));
#endif
#endif

  getExpertForest(owner, 0)->cacheNode(cacheEntry[0]);
  getExpertForest(owner, 1)->cacheNode(cacheEntry[1]);
  getExpertForest(owner, 2)->cacheNode(cacheEntry[2]);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}
#endif


bool mdd_mxd_image_operation::checkTerminals(op_info* op,
  int a, int b, int& c)
{
  assert(false);
  return false;
}


compute_manager::error
mdd_mxd_image_operation::typeCheck(const op_info* owner)
{
  // op1 == MDD, op2 == MXD, op3 = MDD
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] == owner->f[1] || owner->f[0] != owner->f[2])
    return compute_manager::FOREST_MISMATCH;
  if (!getExpertForest(owner, 0)->isMdd() ||
      !getExpertForest(owner, 1)->isMxd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
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
  DCASSERT(owner->nForests == 3 && owner->f[0] == owner->f[2]);
  const int nOperands = 3;
  forest* forests[nOperands] = {owner->f[0], owner->f[0], owner->f[0]};
  op_info* unionOp = 
    smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
    getOpInfo(compute_manager::UNION, forests, nOperands);
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
    smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager());
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


// ---------------------- MDD Traditional Reachability -------------------


mdd_reachability_bfs* mdd_reachability_bfs::getInstance()
{
  static mdd_reachability_bfs instance;
  return &instance;
}


mdd_reachability_bfs::mdd_reachability_bfs()
{ }


mdd_reachability_bfs::~mdd_reachability_bfs() {}


int mdd_reachability_bfs::compute(op_info* owner, int mdd, int mxd)
{
  // set up aliases
  DCASSERT(owner->nForests == 3 && owner->f[0] == owner->f[2]);
  const int nOperands = 3;
  forest* forests[nOperands] = {owner->f[0], owner->f[0], owner->f[0]};
  expert_compute_manager* ecm = 
    smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager());
  assert(ecm != 0);
  op_info* unionOp =
    ecm->getOpInfo(compute_manager::UNION, forests, nOperands);
  assert(unionOp != 0);
  forests[1] = owner->f[1];
  op_info* postImageOp =
    ecm->getOpInfo(compute_manager::POST_IMAGE, forests, nOperands);
  assert(postImageOp != 0);
  expert_forest* mddNm = getExpertForest(owner, 0);

  // Traditional (breadth-first) reachability analysis

#if 1
  expert_forest* mxdNm = getExpertForest(owner, 1);
  dd_edge nsf(mxdNm);
  mxdNm->linkNode(mxd);
  nsf.set(mxd, 0, mxdNm->getNodeLevel(mxd));
  dd_edge reachableStates(mddNm);
  mddNm->linkNode(mdd);
  reachableStates.set(mdd, 0, mddNm->getNodeLevel(mdd));
  dd_edge prevReachableStates(mddNm);

  while(prevReachableStates != reachableStates)
  {
    prevReachableStates = reachableStates;
    // printf("\nPost-Image (mdd:%d, mxd:%d): ",
    //    reachableStates.getNode(), nsf.getNode());
    dd_edge postImage(mddNm);
    ecm->apply(postImageOp, reachableStates, nsf, postImage);
    // printf("%d\n", postImage.getNode());
    // postImage.show(stdout, 2);
    // printf("\nUnion (mdd:%d, mdd:%d): ",
    //    reachableStates.getNode(), postImage.getNode());
    ecm->apply(unionOp, reachableStates, postImage, reachableStates);
    // printf("%d\n", reachableStates.getNode());
  }

  int result = reachableStates.getNode();
  mddNm->linkNode(result);
#else

  mdd_union* unionOpPtr = smart_cast<mdd_union*>(unionOp->op);
  DCASSERT(unionOpPtr != 0);
  mdd_post_image* postImageOpPtr =
    smart_cast<mdd_post_image*>(postImageOp->op);
  DCASSERT(postImageOpPtr != 0);

  int nsf = mxd;
  int reachableStates = mdd;
  int prevReachableStates = 0;
  int postImage = mdd;

  mddNm->linkNode(reachableStates);
  mddNm->linkNode(postImage);

  do {
    int prevPostImage = postImage;
    postImage = postImageOpPtr->compute(postImageOp, postImage, nsf);
    mddNm->unlinkNode(prevPostImage);

    prevReachableStates = reachableStates;
    reachableStates = unionOpPtr->compute(unionOp, reachableStates, postImage);
    mddNm->unlinkNode(prevReachableStates);
  } while (reachableStates != prevReachableStates);

  mddNm->unlinkNode(postImage);

  int result = reachableStates;
  // no need for linkNode(reachableStates) because that has already been
  // called once.

#endif

  return result;
}


// ---------------------- MDD Saturation-based Reachability -------------------


mdd_reachability_dfs* mdd_reachability_dfs::getInstance()
{
  static mdd_reachability_dfs instance;
  return &instance;
}


mdd_reachability_dfs::mdd_reachability_dfs()
{ }


mdd_reachability_dfs::~mdd_reachability_dfs() {}


int mdd_reachability_dfs::compute(op_info* owner, int mdd, int mxd)
{
  DCASSERT(owner->nForests == 3 && owner->f[0] == owner->f[2]);

  // Initialize class members and helper operations
  initialize(owner);

  // Depth-first reachability analysis (Saturation)

#ifdef DEBUG_DFS
  printf("Consolidated Next-State Function:\n");
  xdf->showNodeGraph(stdout, mxd);
  printf("\n");

  printf("Initial State:\n");
  ddf->showNodeGraph(stdout, mdd);
  printf("\n");
#endif

  // Split the next-state function: each level has its own next-state function
  // The nsf is stored into the vector splits
  splitMxd(mxd);

#ifdef DEBUG_DFS
  printf("Split Next-State Function:\n");
  for (int i = splits.size() - 1; i >= 0; i--)
  {
    printf("Level %d, Node %d\n", i, splits[i]);
    xdf->showNodeGraph(stdout, splits[i]);
    printf("\n");
  }

  fflush(stdout);
#endif

#if 1
  // Saturate the node
  int result = saturate(mdd);

  // clear pointers to dd nodes, class members and operation pointers
  clear();

  return result;
#else
  clear();
  return 0;
#endif
}


void mdd_reachability_dfs::initialize(op_info* o)
{
  // set up aliases
  owner = o;
  ecm = smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager());
  assert(ecm != 0);
  ddf = getExpertForest(owner, 0);
  assert(ddf != 0);
  xdf = getExpertForest(owner, 1);
  assert(xdf != 0);
  DCASSERT(ddf->getDomain() == xdf->getDomain());
  ed = smart_cast<expert_domain*>(ddf->useDomain());
  assert(ed != 0);

  // set up mdd operation: union
  const int nOperands = 3;
  forest* forests[nOperands] = {owner->f[0], owner->f[0], owner->f[0]};

  mddUnionOp = ecm->getOpInfo(compute_manager::UNION, forests, nOperands);
  assert(mddUnionOp != 0);
  mddUnion = smart_cast<mdd_union*>(mddUnionOp->op);
  assert(mddUnion != 0);

  // set up mxd operations: intersection and difference
  forests[0] = owner->f[1]; forests[1] = owner->f[1]; forests[2] = owner->f[1];

  mxdIntersectionOp =
    ecm->getOpInfo(compute_manager::INTERSECTION, forests, nOperands);
  assert(mxdIntersectionOp != 0);
  mxdIntersection = smart_cast<mxd_intersection*>(mxdIntersectionOp->op);
  assert(mxdIntersection != 0);

  mxdDifferenceOp =
    ecm->getOpInfo(compute_manager::DIFFERENCE, forests, nOperands);
  assert(mxdDifferenceOp != 0);
  mxdDifference = smart_cast<mxd_difference*>(mxdDifferenceOp->op);
  assert(mxdDifference != 0);

  // Intialize the scratch 2-D vector (i.e. make it the correct size)
  int nLevels = ed->getTopVariable() + 1;
  scratch.clear();
  scratch.resize(nLevels);
  for (unsigned height = ed->getNumVariables(); height > 0; --height)
  {
    int lh = ed->getVariableWithHeight(height);
    int sz = ed->getVariableBound(lh);
    scratch[lh].resize(sz, 0);
  }

  // Initialize the splits vector
  splits.resize(nLevels, 0);

  // Initialize the boolean 2-D vectors
  curr.resize(nLevels);
  for (int i = 0; i < nLevels; ++i)
  {
    curr[i].resize(scratch[i].size(), false);
  }
  next = curr;
}


void mdd_reachability_dfs::clear()
{
  // clear pointer to dd nodes
  for (unsigned i = 0u; i < splits.size(); i++) xdf->unlinkNode(splits[i]);

  // clear class members and pointers to operations
  splits.clear();
  scratch.clear();
  curr.clear();
  next.clear();
  owner = 0;
  ecm = 0;
  ddf = 0;
  xdf = 0;
  ed = 0;
  mddUnionOp = 0;
  mxdIntersectionOp = 0;
  mxdDifferenceOp = 0;
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;
}


// split is used to split a mxd for the saturation algorithm
void mdd_reachability_dfs::splitMxd(int mxd)
{
#if 1
  DCASSERT(xdf != 0);

  int falseNode = xdf->getTerminalNode(false);
  int trueNode = xdf->getTerminalNode(true);

  // find intersection for all mxd[i][i]
  // -- if mxd is smaller than max level size, then some mxd[i] is zero,
  //    therefore intersection is 0.
  //    -- if node is sparse, then there is some mxd[i] = 0,
  //       therefore intersection is 0.
  //    -- if node is full, if size < max level size, intersection is 0.
  // 
  // if intersection == 0, add mxd to level_mxd[level], return.
  // otherwise, 
  // -- create mxdSize nodes at primed level with a copy of corresponding
  //    mxd[i].
  // -- for each new_mxd[i][i],
  //    -- mxd[i][i] = mxd[i][i] - intersection 
  // -- set new_mxd[i] after reducing the primed level nodes
  //    note that new_mxd will never be 0 since mxd is not an identity node
  //
  // add new_mxd to level_mxd[level]
  // 
  // repeat the above for intersection
  //
  int level = 0;
  int intersection = falseNode;
  int mxdSize = 0;
  int mxdI = falseNode;

  xdf->linkNode(mxd);

  while (!xdf->isTerminalNode(mxd)) {
    level = xdf->getNodeLevel(mxd);
    DCASSERT(level > 0); // we only deal with unprimed levels

    // Find intersection for all mxd[i][i]
    // Note: only do this if mxd is a full node; when it is sparse, some
    // mxd[i] is 0 therefore the intersection will always be 0 (falseNode).
    intersection = falseNode;
    if (xdf->isFullNode(mxd)) {
      mxdSize = xdf->getFullNodeSize(mxd);
      if (mxdSize == xdf->getLevelSize(level)) {
        // for all i, mxd[i] != 0
        intersection = trueNode;
        bool first = true;
        for (int i = 0; i < mxdSize && intersection != falseNode; ++i)
        {
          mxdI = xdf->getFullNodeDownPtr(mxd, i);

          // If any mxd[i] is a terminal (according to Identity Reduced rules)
          // it must be node 0, and mxd[i][i] is also 0. Therefore,
          // the intersection is 0. So check for this condition, and break
          // out of the loop it true.

          // if mxdI is a terminal node it must be a 0 (falseNode)
          DCASSERT((xdf->isTerminalNode(mxdI) && mxdI == falseNode) ||
              !xdf->isTerminalNode(mxdI));

          int mxdII = falseNode;

          if (!xdf->isTerminalNode(mxdI)) {
            if (xdf->isFullNode(mxdI)) {
              if (xdf->getFullNodeSize(mxdI) > i)
                mxdII = xdf->getFullNodeDownPtr(mxdI, i);
            } else {
              DCASSERT(xdf->isSparseNode(mxdI));
              // search for ith index
              int found = -1;
              int mxdINnz = xdf->getSparseNodeSize(mxdI);

              if (mxdINnz > 8) {
                // binary search
                int start = 0;
                int stop = mxdINnz - 1;

                while (start < stop) {
                  int mid = (start + stop) / 2;
                  int midIndex = xdf->getSparseNodeIndex(mxdI, mid);
                  if (midIndex < i) {
                    start = (mid == start)? mid + 1: mid;
                  } else {
                    stop = mid;
                  }
                }

                assert(start == stop);
                if (xdf->getSparseNodeIndex(mxdI, start) == i) {
                  found = start;
                }
              }
              else {
                // linear search
                for (int j = 0; j < mxdINnz; ++j)
                {
                  if (xdf->getSparseNodeIndex(mxdI, j) == i) {
                    found = j;
                    break;
                  }
                }
              }

              if (found != -1)
                mxdII = xdf->getSparseNodeDownPtr(mxdI, found);
            }
          }

          if (!first) {
            int temp = getMxdIntersection(intersection, mxdII);
            //mxdIntersection->compute(mxdIntersectionOp, intersection, mxdII);
            xdf->unlinkNode(intersection);
            intersection = temp;
          } else {
            first = false;
            xdf->linkNode(mxdII);
            xdf->unlinkNode(intersection);
            intersection = mxdII;
          }

#ifdef DEBUG_DFS
          printf("intersection: %d level: %d\n",
              intersection, xdf->getNodeLevel(intersection));
#endif
        }
      }
    }

    DCASSERT(splits[level] == falseNode);

    DCASSERT(intersection == falseNode ||
        xdf->getNodeLevel(mxd) > xdf->getNodeLevel(intersection));

    if (intersection != falseNode) {
      splits[level] = getMxdDifference(mxd, intersection);
        // mxdDifference->compute(mxdDifferenceOp, mxd, intersection);
    } else {
      splits[level] = mxd;
      xdf->linkNode(mxd);
    }

    // intersection becomes the mxd for the next iteration
    xdf->unlinkNode(mxd);
    mxd = intersection;
  }

  DCASSERT(xdf->isTerminalNode(mxd));
  xdf->unlinkNode(mxd);
#else
  xdf->linkNode(mxd);
  splits[xdf->getNodeLevel(mxd)] = mxd;
#endif
}

#ifdef ALT_SATURATE_HELPER

int mdd_reachability_dfs::saturate(int mdd)
{
#ifdef DEBUG_DFS
  printf("mdd: %d\n", mdd);
#endif

  // how does saturateHelper get called?
  // bottom-up i.e. call helper on children before calling helper for parent

  DCASSERT(ddf->isReducedNode(mdd));

  // TODO: need to define compute table for storing recFire results
  // mdd_mxd_image_operation already has a compute table, so use that.

  // terminal condition for recursion
  if (ddf->isTerminalNode(mdd)) return mdd;

  int k = ddf->getNodeLevel(mdd);      // level
  int sz = ddf->getLevelSize(k);       // size

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d\n", mdd, k, sz);
#endif

  DCASSERT(curr.size() > unsigned(k));
  DCASSERT(next.size() > unsigned(k));
  DCASSERT(curr[k].size() >= unsigned(sz));
  DCASSERT(next[k].size() >= unsigned(sz));

  int n = ddf->createTempNode(k, sz, true);  // new node

  // recursively call saturate for children of mdd
  if (ddf->isFullNode(mdd)) {
    int mddSize = ddf->getFullNodeSize(mdd);
    for (int i = 0; i < mddSize; ++i)
    {
      int mddI = ddf->getFullNodeDownPtr(mdd, i);
      if (mddI != 0) {
        int temp = saturate(mddI);
        ddf->setDownPtrWoUnlink(n, i, temp);
        ddf->unlinkNode(temp);
      }
    }
  } else {
    DCASSERT(ddf->isSparseNode(mdd));
    int mddNnz = ddf->getSparseNodeSize(mdd);
    for (int i = 0; i < mddNnz; ++i)
    {
      int mddIDownPtr = ddf->getSparseNodeDownPtr(mdd, i);
      int mddIIndex = ddf->getSparseNodeIndex(mdd, i);
      int temp = saturate(mddIDownPtr);
      ddf->setDownPtrWoUnlink(n, mddIIndex, temp);
      ddf->unlinkNode(temp);
    }
  }
  
  // call saturateHelper for n
#ifdef DEBUG_DFS
  printf("Calling saturate: level %d\n", k);
#endif
  saturateHelper(n);

  // reduce and return
  n = ddf->reduceNode(n);

#ifdef DEBUG_DFS
  ddf->showNodeGraph(stdout, n);
#endif

  return n;
}


void mdd_reachability_dfs::saturateHelperUnPrimeFull(int mdd, int mxd)
{
  DCASSERT(xdf->isFullNode(mxd));
  DCASSERT(ddf->isFullNode(mdd));

  int mddSize = ddf->getFullNodeSize(mdd);
  int mxdSize = xdf->getFullNodeSize(mxd);
  int minSize = MIN(mddSize, mxdSize);

  std::vector<bool> enabled(mddSize, true);
  bool repeat = true;

  while (repeat)
  {
    std::vector<bool> next(mddSize, false);
    for (int i = 0; i < minSize; i++)
    {
      if (!enabled[i]) continue;

      int mxdI = xdf->getFullNodeDownPtr(mxd, i);
      if (mxdI == 0) continue;

      if (xdf->isFullNode(mxdI)) {
        // mxdI is truncated-full
        saturateHelperPrimeFull(mdd, i, mxdI, next);
      }
      else {
        // mxdI is sparse
        saturateHelperPrimeSparse(mdd, i, mxdI, next);
      }
    }

    repeat = false;
    for (int i = 0; i < mddSize; i++)
    {
      if (next[i]) {
        enabled = next;
        repeat = true;
        break;
      }
    }
  } // while(repeat)
}


void mdd_reachability_dfs::saturateHelperUnPrimeSparse(int mdd, int mxd)
{
  DCASSERT(xdf->isSparseNode(mxd));
  DCASSERT(ddf->isFullNode(mdd));

  int mddSize = ddf->getFullNodeSize(mdd);
  int mxdNnz = xdf->getSparseNodeSize(mxd);

  // mdd must of the max size for this level
  DCASSERT(mddSize == ddf->getLevelSize(ddf->getNodeLevel(mdd)));

  std::vector<bool> enabled(mddSize, true);
  bool repeat = true;

  while (repeat)
  {
    std::vector<bool> next(mddSize, false);

    for (int index = 0; index < mxdNnz; index++)
    {
      int i = xdf->getSparseNodeIndex(mxd, index);
      
      if (!enabled[i]) continue;

      int mxdI = xdf->getSparseNodeDownPtr(mxd, index);

      if (xdf->isFullNode(mxdI)) {
        // mxdI is truncated-full
        saturateHelperPrimeFull(mdd, i, mxdI, next);
      }
      else {
        // mxdI is sparse
        saturateHelperPrimeSparse(mdd, i, mxdI, next);
      }
    }

    repeat = false;
    for (int i = 0; i < mddSize; i++)
    {
      if (next[i]) {
        enabled = next;
        repeat = true;
        break;
      }
    }
  } // while(repeat)
}


void mdd_reachability_dfs::saturateHelperPrimeFull(int mdd, int i,
    int mxdI, std::vector<bool>& next)
{
  DCASSERT(ddf->isFullNode(mdd));
  DCASSERT(!ddf->isReducedNode(mdd));
  DCASSERT(!xdf->isTerminalNode(mxdI));
  DCASSERT(xdf->isFullNode(mxdI));
  // Note: it does not matter if mddI is sparse or full since we do not expand
  // it here.

  int mddI = ddf->getFullNodeDownPtr(mdd, i);
  if (mddI == 0) return;

  int mxdISize = xdf->getFullNodeSize(mxdI);
  for (int j = 0; j < mxdISize; j++)
  {
    int mxdIJ = xdf->getFullNodeDownPtr(mxdI, j);
    if (mxdIJ == 0) continue;

    int f = recFire(mddI, mxdIJ);
    if (f == 0) continue;

    int mddJ = ddf->getFullNodeDownPtr(mdd, j);
    int u = getMddUnion(f, mddJ);
    ddf->unlinkNode(f);

    if (u != mddJ) {
      ddf->setDownPtr(mdd, j, u);
      if (i == j) mddI = u;
      next[j] = true;
    }
    ddf->unlinkNode(u);
  }
}


void mdd_reachability_dfs::saturateHelperPrimeSparse(int mdd, int i,
    int mxdI, std::vector<bool>& next)
{
  DCASSERT(ddf->isFullNode(mdd));
  DCASSERT(!ddf->isReducedNode(mdd));
  DCASSERT(!xdf->isTerminalNode(mxdI));
  DCASSERT(xdf->isSparseNode(mxdI));
  // Note: it doesi not matter if mddI is sparse or full since we do not expand
  // it here.

  int mddI = ddf->getFullNodeDownPtr(mdd, i);
  if (mddI == 0) return;

  int mxdINnz = xdf->getSparseNodeSize(mxdI);
  for (int j = 0; j < mxdINnz; j++)
  {
    int jDown = xdf->getSparseNodeDownPtr(mxdI, j);
    int jIndex = xdf->getSparseNodeIndex(mxdI, j);

    int f = recFire(mddI, jDown);
    if (f == 0) continue;

    int mddJIndex = ddf->getFullNodeDownPtr(mdd, jIndex);
    int u = getMddUnion(f, mddJIndex);
    ddf->unlinkNode(f);

    if (u != mddJIndex) {
      ddf->setDownPtr(mdd, jIndex, u);
      if (i == jIndex) mddI = u;
      next[jIndex] = true;
    }
    ddf->unlinkNode(u);
  }
}


#if 0

void mdd_reachability_dfs::saturateHelper(int mdd)
{
  // ============== BUG ===============
  DCASSERT(ddf->getTerminalNode(false) == 0);
  DCASSERT(ddf->getTerminalNode(true) == -1);

  DCASSERT(!ddf->isReducedNode(mdd));
  DCASSERT(ddf->isFullNode(mdd));
  DCASSERT(ddf->getLevelSize(ddf->getNodeLevel(mdd)) ==
      ddf->getFullNodeSize(mdd));

  int mddLevel = ddf->getNodeLevel(mdd);
  int mxd = splits[mddLevel];
  if (mxd == 0) return;

  // int *mdd_ptr = ddf->get_full_down_ptrs(mdd);
  int mddSize = ddf->getFullNodeSize(mdd);

  // curr[mddLevel];
  // next[mddLevel];
  fill(next[mddLevel].begin(), next[mddLevel].end(), true);
  bool repeat = true;

  if (xdf->isFullNode(mxd)) {
    int mxdSize = xdf->getFullNodeSize(mxd);
    DCASSERT(mxdSize <= mddSize);
    int min = MIN(mddSize, mxdSize);

    while (repeat)
    {
      curr[mddLevel] = next[mddLevel];
      fill(next[mddLevel].begin(), next[mddLevel].end(), false);
      repeat = false;
      // for each mxd[i] != 0
      for (int i = 0; i < min; i++)
      {
        int mddI = ddf->getFullNodeDownPtr(mdd, i);
        int mxdI = xdf->getFullNodeDownPtr(mxd, i);
        if (mddI == 0 || mxdI == 0 || !curr[mddLevel][i]) continue;

        DCASSERT(ddf->isReducedNode(mddI));
        DCASSERT(xdf->isReducedNode(mxdI));

        if (xdf->isFullNode(mxdI)) {
          int mxdISize = xdf->getFullNodeSize(mxdI);
          // for each mxd[i][j] != 0
          for (int j = 0; j < mxdISize; j++)
          {
            int mxdIJ = xdf->getFullNodeDownPtr(mxdI, j);
            if (mxdIJ == 0) continue;

            DCASSERT(xdf->isReducedNode(mxdIJ));
            int f = recFire(mddI, mxdIJ);
            if (f == 0) continue;

            int mddJ = ddf->getFullNodeDownPtr(mdd, j);
            DCASSERT(ddf->isReducedNode(mddJ));
            DCASSERT(ddf->isReducedNode(f));
            int u = mddUnion->compute(mddUnionOp, f, mddJ);
            DCASSERT(ddf->isReducedNode(u));
            ddf->unlinkNode(f);

            if (u != mddJ) {
              ddf->setDownPtr(mdd, j, u);
              if (i == j) mddI = u;
#if 1
              if (!ddf->isTerminalNode(u) && j < min) {
                if (j > i) {
                  // call recFire for it later in this iteration itself
                  curr[mddLevel][j] = true;
                }
                else {
                  // call recFire for it in the next iteration
                  next[mddLevel][j] = true;
                  repeat = true;
                }
              }
#else
              next[mddLevel][j] = true;
              repeat = true;
#endif
            }
            ddf->unlinkNode(u);
          } // for (j=0; j<mxdISize; j++)
        }
        else {
          DCASSERT(xdf->isSparseNode(mxdI));
          int mxdINnz = xdf->getSparseNodeSize(mxdI);
          // for each mxd[i][j] != 0
          for (int j = 0; j < mxdINnz; ++j)
          {
            int jIndex = xdf->getSparseNodeIndex(mxdI, j);
            int jDown = xdf->getSparseNodeDownPtr(mxdI, j);
            DCASSERT(xdf->isReducedNode(jDown));
            int f = recFire(mddI, jDown);
            if (f == 0) continue;

            int mddJIndex = ddf->getFullNodeDownPtr(mdd, jIndex);
            DCASSERT(ddf->isReducedNode(mddJIndex));
            DCASSERT(ddf->isReducedNode(f));
            int u = mddUnion->compute(mddUnionOp, f, mddJIndex);
            DCASSERT(ddf->isReducedNode(u));
            ddf->unlinkNode(f);

            if (u != mddJIndex) {
              // update mddJIndex
              ddf->setDownPtr(mdd, jIndex, u);
              if (i == jIndex) mddI = u;
#if 1
              if (!ddf->isTerminalNode(u) && jIndex < min) {
                if (jIndex > i) {
                  // call recFire for it later in this iteration itself
                  curr[mddLevel][jIndex] = true;
                }
                else {
                  // call recFire for it in the next iteration
                  next[mddLevel][jIndex] = true;
                  repeat = true;
                }
              }
#else
              next[mddLevel][j] = true;
              repeat = true;
#endif
            }
            ddf->unlinkNode(u);
          } // for (int j = 0; j < mxdINnz; ++j)
        } // mxdI is sparse
      } // for (i=0; i<mddSize && i<mxdSize; i++)
    } // while (repeat)
  }
  else {
    DCASSERT (xdf->isSparseNode(mxd));

    int mxdNnz = xdf->getSparseNodeSize(mxd);
    int mxdSize = xdf->getSparseNodeIndex(mxd, mxdNnz - 1) + 1;

    DCASSERT(mxdSize <= mddSize);
    int min = MIN(mddSize, mxdSize);

    while (repeat) {
      curr[mddLevel] = next[mddLevel];
      fill(next[mddLevel].begin(), next[mddLevel].end(), false);
      repeat = false;

      for (int index = 0; index < mxdNnz; index++) // same as i < min
      {
        int i = xdf->getSparseNodeIndex(mxd, index);
        if (i >= mddSize) break;
        int mxdI = xdf->getSparseNodeDownPtr(mxd, index);
        int mddI = ddf->getFullNodeDownPtr(mdd, i);
        if (mddI == 0 || !curr[mddLevel][i]) continue;

        // for each mxd[i] != 0
        if (xdf->isFullNode(mxdI)) {
          int mxdISize = xdf->getFullNodeSize(mxdI);
          // for each mxd[i][j] != 0
          for (int j = 0; j < mxdISize; j++)
          {
            int mxdIJ = xdf->getFullNodeDownPtr(mxdI, j);
            if (mxdIJ == 0) continue;
            DCASSERT(xdf->isReducedNode(mxdIJ));
            DCASSERT(ddf->isReducedNode(mddI));
            int f = recFire(mddI, mxdIJ);
            if (f == 0) continue;

            int mddJ = ddf->getFullNodeDownPtr(mdd, j);
            DCASSERT(ddf->isReducedNode(mddJ));
            DCASSERT(ddf->isReducedNode(f));
            int u = mddUnion->compute(mddUnionOp, f, mddJ);
            DCASSERT(ddf->isReducedNode(u));
            ddf->unlinkNode(f);

            if (u != mddJ) {
              ddf->setDownPtr(mdd, j, u);
              if (i == j) mddI = u;
#if 1
              if (!ddf->isTerminalNode(u) && j < min) {
                if (j > i) {
                  // call recFire for it later in this iteration itself
                  curr[mddLevel][j] = true;
                }
                else {
                  // call recFire for it in the next iteration
                  next[mddLevel][j] = true;
                  repeat = true;
                }
              }
#else
              next[mddLevel][j] = true;
              repeat = true;
#endif
            }
            ddf->unlinkNode(u);
          } // for (j=0; j<mxdISize; j++)
        }
        else {
          DCASSERT(xdf->isSparseNode(mxdI));
          int mxdINnz = xdf->getSparseNodeSize(mxdI);
          // for each mxd[i][j] != 0
          for (int j = 0; j < mxdINnz; ++j)
          {
            int jIndex = xdf->getSparseNodeIndex(mxdI, j);
            int jDown = xdf->getSparseNodeDownPtr(mxdI, j);
            DCASSERT(xdf->isReducedNode(jDown));
            int f = recFire(mddI, jDown);
            if (f == 0) continue;

            int mddJIndex = ddf->getFullNodeDownPtr(mdd, jIndex);
            DCASSERT(ddf->isReducedNode(mddJIndex));
            DCASSERT(ddf->isReducedNode(f));
            int u = mddUnion->compute(mddUnionOp, f, mddJIndex);
            DCASSERT(ddf->isReducedNode(u));
            ddf->unlinkNode(f);

            if (u != mddJIndex) {
              // update mddJIndex
              ddf->setDownPtr(mdd, jIndex, u);
              if (i == jIndex) mddI = u;
#if 1
              if (!ddf->isTerminalNode(u) && jIndex < min) {
                if (jIndex > i) {
                  // call recFire for it later in this iteration itself
                  curr[mddLevel][jIndex] = true;
                }
                else {
                  // call recFire for it in the next iteration
                  next[mddLevel][jIndex] = true;
                  repeat = true;
                }
              }
#else
              next[mddLevel][j] = true;
              repeat = true;
#endif
            }
            ddf->unlinkNode(u);
          } // for (int j = 0; j < mxdINnz; ++j)
        } // mxdI is sparse
      } // for (index=0; index<nnz; index++)
    } // while (repeat)
  } // mxd is sparse
}

#else

void mdd_reachability_dfs::saturateHelper(int mdd)
{
  // ============== BUG ===============
  DCASSERT(ddf->getTerminalNode(false) == 0);
  DCASSERT(ddf->getTerminalNode(true) == -1);

  DCASSERT(!ddf->isReducedNode(mdd));
  DCASSERT(ddf->isFullNode(mdd));
  DCASSERT(ddf->getLevelSize(ddf->getNodeLevel(mdd)) ==
      ddf->getFullNodeSize(mdd));

  int mddLevel = ddf->getNodeLevel(mdd);
  int mxd = splits[mddLevel];
  if (mxd == 0) return;

  if (xdf->isFullNode(mxd))
    saturateHelperUnPrimeFull(mdd, mxd);
  else
    saturateHelperUnPrimeSparse(mdd, mxd);
}


#endif


#if 0

int mdd_reachability_dfs::recFire(int mdd, int mxd)
{
#if 1
  DCASSERT(ddf->isReducedNode(mdd));
  DCASSERT(xdf->isReducedNode(mxd));

  DCASSERT(0 == xdf->getTerminalNode(false));
  DCASSERT(-1 == xdf->getTerminalNode(true));
  DCASSERT(0 == ddf->getTerminalNode(false));
  DCASSERT(-1 == ddf->getTerminalNode(true));

  if (mxd == -1) {
    ddf->linkNode(mdd);
    return mdd;
  }

  if (mxd == 0 || mdd == 0) {
    return 0;
  }

  int result = 0;
#if 1
  if (findResult(owner, mdd, mxd, result)) {
    return result;
  }
#endif

  int mxdHeight = xdf->getNodeHeight(mxd);
  int mddHeight = ddf->getNodeHeight(mdd);
  int newHeight = MAX(mxdHeight, mddHeight);
  int newLevel = ed->getVariableWithHeight(newHeight);

#if 0
  int mxdLevel = xdf->getNodeLevel(mxd);
  int mddLevel = ddf->getNodeLevel(mdd);
  int newLevel = MAX(mxdLevel, mddLevel);
#endif

  int newSize = ddf->getLevelSize(newLevel);
  int newNode = ddf->createTempNode(newLevel, newSize, true);

  if (mxdHeight < mddHeight) {
    if (ddf->isFullNode(mdd)) {
      int mddSize = ddf->getFullNodeSize(mdd);
      for (int i = 0; i < mddSize; i++)
      {
        int mddI = ddf->getFullNodeDownPtr(mdd, i);
        DCASSERT(ddf->isReducedNode(mddI));
        int temp = recFire(mddI, mxd); 
        DCASSERT(ddf->isReducedNode(temp));
        ddf->setDownPtrWoUnlink(newNode, i, temp);
        ddf->unlinkNode(temp);
      }
    } else {
      DCASSERT(ddf->isSparseNode(mdd));
      int mddNnz = ddf->getSparseNodeSize(mdd);
      for (int i = 0; i < mddNnz; i++)
      {
        int index = ddf->getSparseNodeIndex(mdd, i);
        int dptr = ddf->getSparseNodeDownPtr(mdd, i);
        int temp = recFire(dptr, mxd);
        DCASSERT(ddf->isReducedNode(temp));
        ddf->setDownPtrWoUnlink(newNode, index, temp);
        ddf->unlinkNode(temp);
      }
    }
  } else if (mxdHeight > mddHeight) {
    if (xdf->isFullNode(mxd)) {
      int mxdSize = xdf->getFullNodeSize(mxd);
      for (int i = 0; i < mxdSize; i++)
      {
        int mxdI = xdf->getFullNodeDownPtr(mxd, i);
        if (mxdI == 0) continue;
        if (xdf->isFullNode(mxdI)) {
          int mxdISize = xdf->getFullNodeSize(mxdI);
          for (int j = 0; j < mxdISize; j++)
          {
            int mxdIJ = xdf->getFullNodeDownPtr(mxdI, j);
            if (mxdIJ == 0) continue;

            int f = recFire(mdd, mxdIJ);
            if (f == 0) continue;

            int dptr = ddf->getFullNodeDownPtr(newNode, j);
            if (f != dptr) {
              int u = mddUnion->compute(mddUnionOp, f, dptr);
              ddf->setDownPtr(newNode, j, u);
              ddf->unlinkNode(u);
            }
            ddf->unlinkNode(f);
          }
        }
        else {
          DCASSERT(xdf->isSparseNode(mxdI));
          int mxdINnz = xdf->getSparseNodeSize(mxdI);
          for (int j = 0; j < mxdINnz; j++)
          {
            int mxdIJIndex = xdf->getSparseNodeIndex(mxdI, j);
            int mxdIJ = xdf->getSparseNodeDownPtr(mxdI, j);

            int f = recFire(mdd, mxdIJ);
            if (f == 0) continue;

            int dptr = ddf->getFullNodeDownPtr(newNode, mxdIJIndex);
            if (f != dptr) {
              int u = mddUnion->compute(mddUnionOp, f, dptr);
              ddf->setDownPtr(newNode, mxdIJIndex, u);
              ddf->unlinkNode(u);
            }
            ddf->unlinkNode(f);
          }
        }
      }
    }
    else {
      DCASSERT(xdf->isSparseNode(mxd));

      int mxdNnz = xdf->getSparseNodeSize(mxd);
      for (int index = 0; index < mxdNnz; index++)
      {
        // int i = xdf->getSparseNodeIndex(mxd, index);
        int mxdI = xdf->getSparseNodeDownPtr(mxd, index);

        if (xdf->isFullNode(mxdI)) {
          int mxdISize = xdf->getFullNodeSize(mxdI);
          for (int j = 0; j < mxdISize; j++)
          {
            int mxdIJ = xdf->getFullNodeDownPtr(mxdI, j);
            if (mxdIJ == 0) continue;

            int f = recFire(mdd, mxdIJ);
            if (f == 0) continue;

            int dptr = ddf->getFullNodeDownPtr(newNode, j);
            if (f != dptr) {
              int u = mddUnion->compute(mddUnionOp, f, dptr);
              ddf->setDownPtr(newNode, j, u);
              ddf->unlinkNode(u);
            }
            ddf->unlinkNode(f);
          }
        }
        else {
          DCASSERT(xdf->isSparseNode(mxdI));
          int mxdINnz = xdf->getSparseNodeSize(mxdI);
          for (int j = 0; j < mxdINnz; j++)
          {
            int mxdIJIndex = xdf->getSparseNodeIndex(mxdI, j);
            int mxdIJ = xdf->getSparseNodeDownPtr(mxdI, j);

            int f = recFire(mdd, mxdIJ);
            if (f == 0) continue;

            int dptr = ddf->getFullNodeDownPtr(newNode, mxdIJIndex);
            if (f != dptr) {
              int u = mddUnion->compute(mddUnionOp, f, dptr);
              ddf->setDownPtr(newNode, mxdIJIndex, u);
              ddf->unlinkNode(u);
            }
            ddf->unlinkNode(f);
          }
        }
      }
    }
  } else {
    DCASSERT(mxdHeight == mddHeight);

    if (ddf->isFullNode(mdd)) {
      int mddSize = ddf->getFullNodeSize(mdd);
      if (xdf->isFullNode(mxd)) {
        int mxdSize = xdf->getFullNodeSize(mxd);
        int min = MIN(mddSize, mxdSize);
        for (int i = 0; i < min; i++)
        {
          int mddI = ddf->getFullNodeDownPtr(mdd, i);
          int mxdI = xdf->getFullNodeDownPtr(mxd, i);
          if (mxdI == 0 || mddI == 0) continue;
          if (xdf->isFullNode(mxdI)) {
            int mxdISize = xdf->getFullNodeSize(mxdI);
            for (int j = 0; j < mxdISize; j++)
            {
              int mxdIJ = xdf->getFullNodeDownPtr(mxdI, j);
              if (mxdIJ == 0) continue;

              int f = recFire(mddI, mxdIJ);
              if (f == 0) continue;
              if (f == -1) {
                ddf->setDownPtr(newNode, j, f);
                continue;
              }

              int dptr = ddf->getFullNodeDownPtr(newNode, j);
              if (f != dptr) {
                int u = mddUnion->compute(mddUnionOp, f, dptr);
                ddf->setDownPtr(newNode, j, u);
                ddf->unlinkNode(u);
              }
              ddf->unlinkNode(f);
            }
          }
          else {
            DCASSERT(xdf->isSparseNode(mxdI));
            int mxdINnz = xdf->getSparseNodeSize(mxdI);
            for (int j = 0; j < mxdINnz; j++)
            {
              int mxdIJIndex = xdf->getSparseNodeIndex(mxdI, j);
              int mxdIJ = xdf->getSparseNodeDownPtr(mxdI, j);
              int f = recFire(mddI, mxdIJ);

              if (f == 0) continue;
              if (f == -1) {
                ddf->setDownPtr(newNode, mxdIJIndex, f);
                continue;
              }

              int dptr = ddf->getFullNodeDownPtr(newNode, mxdIJIndex);
              if (f != dptr) {
                int u = mddUnion->compute(mddUnionOp, f, dptr);
                ddf->setDownPtr(newNode, mxdIJIndex, u);
                ddf->unlinkNode(u);
              }
              ddf->unlinkNode(f);
            }
          }
        }
      }
      else {
        DCASSERT(xdf->isSparseNode(mxd));

        int mxdNnz = xdf->getSparseNodeSize(mxd);
        for (int index = 0; index < mxdNnz; index++)
        {
          int i = xdf->getSparseNodeIndex(mxd, index);
          if (i >= mddSize) break;
          int mddI = ddf->getFullNodeDownPtr(mdd, index);
          if (mddI == 0) continue;
          int mxdI = xdf->getSparseNodeDownPtr(mxd, index);

          if (xdf->isFullNode(mxdI)) {
            int mxdISize = xdf->getFullNodeSize(mxdI);
            for (int j = 0; j < mxdISize; j++)
            {
              int mxdIJ = xdf->getFullNodeDownPtr(mxdI, j);
              if (mxdIJ == 0) continue;

              int f = recFire(mddI, mxdIJ);
              if (f == 0) continue;
              if (f == -1) {
                ddf->setDownPtr(newNode, j, f);
                continue;
              }

              int dptr = ddf->getFullNodeDownPtr(newNode, j);
              if (f != dptr) {
                int u = mddUnion->compute(mddUnionOp, f, dptr);
                ddf->setDownPtr(newNode, j, u);
                ddf->unlinkNode(u);
              }
              ddf->unlinkNode(f);
            }
          }
          else {
            DCASSERT(xdf->isSparseNode(mxdI));
            int mxdINnz = xdf->getSparseNodeSize(mxdI);
            for (int j = 0; j < mxdINnz; j++)
            {
              int mxdIJIndex = xdf->getSparseNodeIndex(mxdI, j);
              int mxdIJ = xdf->getSparseNodeDownPtr(mxdI, j);
              int f = recFire(mddI, mxdIJ);

              if (f == 0) continue;
              if (f == -1) {
                ddf->setDownPtr(newNode, mxdIJIndex, f);
                continue;
              }

              int dptr = ddf->getFullNodeDownPtr(newNode, mxdIJIndex);
              if (f != dptr) {
                int u = mddUnion->compute(mddUnionOp, f, dptr);
                ddf->setDownPtr(newNode, mxdIJIndex, u);
                ddf->unlinkNode(u);
              }
              ddf->unlinkNode(f);
            }
          }
        }
      }
    }
    else {
      // ========= BUG?? ================

      DCASSERT(ddf->isSparseNode(mdd));
      int mddNnz = ddf->getSparseNodeSize(mdd);
      int mddSize = ddf->getSparseNodeIndex(mdd, mddNnz - 1) + 1;
      // pointer to traverse the sparse mdd node
      int mddIndex = 0;
      // the node index to which mddIndex refers to
      int mddIIndex = ddf->getSparseNodeIndex(mdd, mddIndex);

      if (xdf->isFullNode(mxd)) {
        int mxdSize = xdf->getFullNodeSize(mxd);
        int min = MIN(mddSize, mxdSize);
        for (int i = 0; i < min; i++)
        {
          if (i < mddIIndex) {
            i = mddIIndex - 1;
            continue;
          }
          if (i > mddIIndex) {
            mddIndex++;
            if (mddIndex >= mddNnz) break;
            mddIIndex = ddf->getSparseNodeIndex(mdd, mddIndex);
            i--;
            continue;
          }
          DCASSERT(i == mddIIndex);

          int mxdI = xdf->getFullNodeDownPtr(mxd, i);
          if (mxdI == 0) continue;
          int mddI = ddf->getSparseNodeDownPtr(mdd, mddIndex);

          if (xdf->isFullNode(mxdI)) {
            int mxdISize = xdf->getFullNodeSize(mxdI);
            for (int j = 0; j < mxdISize; j++)
            {
              int mxdIJ = xdf->getFullNodeDownPtr(mxdI, j);
              if (mxdIJ == 0) continue;

              int f = recFire(mddI, mxdIJ);
              if (f == 0) continue;

              int dptr = ddf->getFullNodeDownPtr(newNode, j);
              if (f != dptr) {
                int u = mddUnion->compute(mddUnionOp, f, dptr);
                ddf->setDownPtr(newNode, j, u);
                ddf->unlinkNode(u);
              }
              ddf->unlinkNode(f);
            }
          }
          else {
            DCASSERT(xdf->isSparseNode(mxdI));
            int mxdINnz = xdf->getSparseNodeSize(mxdI);
            for (int j = 0; j < mxdINnz; j++)
            {
              int mxdIJIndex = xdf->getSparseNodeIndex(mxdI, j);
              int mxdIJ = xdf->getSparseNodeDownPtr(mxdI, j);

              int f = recFire(mddI, mxdIJ);
              if (f == 0) continue;

              int dptr = ddf->getFullNodeDownPtr(newNode, mxdIJIndex);
              if (f != dptr) {
                int u = mddUnion->compute(mddUnionOp, f, dptr);
                ddf->setDownPtr(newNode, mxdIJIndex, u);
                ddf->unlinkNode(u);
              }
              ddf->unlinkNode(f);
            }
          }
        }
      }
      else {
        DCASSERT(xdf->isSparseNode(mxd));

        int mxdNnz = xdf->getSparseNodeSize(mxd);
        for (int index = 0; index < mxdNnz; index++)
        {
          int i = xdf->getSparseNodeIndex(mxd, index);
          if (i < mddIIndex) continue;
          if (i > mddIIndex) {
            mddIndex++;
            if (mddIndex >= mddNnz) break;
            mddIIndex = ddf->getSparseNodeIndex(mdd, mddIndex);
            index--;
            continue;
          }

          assert(i == mddIIndex);

          int mddI = ddf->getSparseNodeDownPtr(mdd, mddIndex);
          int mxdI = xdf->getSparseNodeDownPtr(mxd, index);

          if (xdf->isFullNode(mxdI)) {
            int mxdISize = xdf->getFullNodeSize(mxdI);
            for (int j = 0; j < mxdISize; j++)
            {
              int mxdIJ = xdf->getFullNodeDownPtr(mxdI, j);
              if (mxdIJ == 0) continue;

              int f = recFire(mddI, mxdIJ);
              if (f == 0) continue;

              int dptr = ddf->getFullNodeDownPtr(newNode, j);
              if (f != dptr) {
                int u = mddUnion->compute(mddUnionOp, f, dptr);
                ddf->setDownPtr(newNode, j, u);
                ddf->unlinkNode(u);
              }
              ddf->unlinkNode(f);
            }
          }
          else {
            DCASSERT(xdf->isSparseNode(mxdI));
            int mxdINnz = xdf->getSparseNodeSize(mxdI);
            for (int j = 0; j < mxdINnz; j++)
            {
              int mxdIJIndex = xdf->getSparseNodeIndex(mxdI, j);
              int mxdIJ = xdf->getSparseNodeDownPtr(mxdI, j);

              int f = recFire(mddI, mxdIJ);
              if (f == 0) continue;

              int dptr = ddf->getFullNodeDownPtr(newNode, mxdIJIndex);
              if (f != dptr) {
                int u = mddUnion->compute(mddUnionOp, f, dptr);
                ddf->setDownPtr(newNode, mxdIJIndex, u);
                ddf->unlinkNode(u);
              }
              ddf->unlinkNode(f);
            }
          }
        }
      }
    }
  }

  int newNode0 = ddf->getFullNodeDownPtr(newNode, 0);
  int i = 0;
  if (ddf->isTerminalNode(newNode0)) {
    for (i = 1;
        i < newSize && newNode0 == ddf->getFullNodeDownPtr(newNode, i); ++i);
  }
  if (i != newSize) {
    saturateHelper(newNode);
  }

  result = ddf->reduceNode(newNode);
  saveResult(owner, mdd, mxd, result);
  return result;
#else
  return 0;
#endif
}

#else


int mdd_reachability_dfs::recFire(int mdd, int mxd)
{
  if (xdf->isTerminalNode(mxd)) {
    if (mxd == 0)
      return 0;
    ddf->linkNode(mdd);
    return mdd;
  }

  // POSSIBLE BUG
  // when mdd == -1 && mxd is not a terminal, should this recursion stop?
  if (ddf->isTerminalNode(mdd)) {
    return mdd;
  }

  int result = 0;
  if (findResult(owner, mdd, mxd, result)) {
    return result;
  }

  int mddHeight = ddf->getNodeHeight(mdd);
  int mxdHeight = xdf->getNodeHeight(mxd);

  if (mddHeight < mxdHeight) {
    int resultLevel = xdf->getNodeLevel(mxd);
    result = ddf->createTempNodeMaxSize(resultLevel, true);
    recFireExpandMxd(mdd, mxd, result);
  }
  else if (mddHeight > mxdHeight) {
    int resultLevel = ddf->getNodeLevel(mdd);
    result = ddf->createTempNodeMaxSize(resultLevel, true);
    recFireExpandMdd(mdd, mxd, result);
  }
  else {
    // mddHeight == mxdHeight
    int resultLevel = ddf->getNodeLevel(mdd);
    result = ddf->createTempNodeMaxSize(resultLevel, true);

    // check for full/sparse and call appropriate function
    if (ddf->isFullNode(mdd))
      if (xdf->isFullNode(mxd))
        recFireFF(mdd, mxd, result);
      else
        recFireFS(mdd, mxd, result);
    else
      if (xdf->isFullNode(mxd))
        recFireSF(mdd, mxd, result);
      else
        recFireSS(mdd, mxd, result);
  }

  int resultSize = ddf->getFullNodeSize(result);
  int result0 = ddf->getFullNodeDownPtr(result, 0);
  int i = 0;
  if (ddf->isTerminalNode(result0)) {
    for (i = 1;
        i < resultSize && result0 == ddf->getFullNodeDownPtr(result, i); ++i);
  }
  if (i != resultSize) {
    saturateHelper(result);
  }

  result = ddf->reduceNode(result);
  saveResult(owner, mdd, mxd, result);
  return result;
}


void mdd_reachability_dfs::recFireExpandMdd(int mdd, int mxd, int result)
{
  DCASSERT(ddf->isFullNode(result));
  DCASSERT(!ddf->isReducedNode(result));

  if (ddf->isFullNode(mdd)) {
    int mddSize = ddf->getFullNodeSize(mdd);
    for (int i = 0; i < mddSize; i++)
    {
      int mddI = ddf->getFullNodeDownPtr(mdd, i);
      if (mddI != 0) {
        int f = recFire(mddI, mxd);
        if (f != 0) {
          // this corresponds to mxd[i][i] so add it to result[i]
          // but result[i] is 0 before now, so set result[i] to f.
          ddf->setDownPtr(result, i, f);
        }
        ddf->unlinkNode(f);
      }
    }
  }
  else {
    DCASSERT(ddf->isSparseNode(mdd));
    int mddNnz = ddf->getSparseNodeSize(mdd);
    for (int i = 0; i < mddNnz; i++)
    {
      int mddI = ddf->getSparseNodeDownPtr(mdd, i);
      int f = recFire(mddI, mxd);
      if (f != 0) {
        // this corresponds to mxd[i][i] so add it to result[i]
        // but result[i] is 0 before now, so set result[i] to f.
        ddf->setDownPtr(result, ddf->getSparseNodeIndex(mdd, i), f);
      }
      ddf->unlinkNode(f);
    }
  }
}


void mdd_reachability_dfs::recFireExpandMxd(int mdd, int mxd, int result)
{
  DCASSERT(ddf->isFullNode(result));
  DCASSERT(!ddf->isReducedNode(result));

  if (xdf->isFullNode(mxd)) {
    int mxdSize = xdf->getFullNodeSize(mxd);
    for (int i = 0; i < mxdSize; i++)
    {
      int mxdI = xdf->getFullNodeDownPtr(mxd, i);
      if (mxdI != 0) {
        DCASSERT(!xdf->isTerminalNode(mxdI));
        recFirePrime(mdd, mxdI, result);
      }
    }
  }
  else {
    DCASSERT(xdf->isSparseNode(mxd));

    int mxdNnz = xdf->getSparseNodeSize(mxd);
    for (int i = 0; i < mxdNnz; i++)
    {
      int mxdI = xdf->getSparseNodeDownPtr(mxd, i);
      DCASSERT(!xdf->isTerminalNode(mxdI));
      recFirePrime(mdd, mxdI, result);
    }
  }
}


void mdd_reachability_dfs::recFireFF(int mdd, int mxd, int result)
{
  DCASSERT(ddf->isFullNode(mdd));
  DCASSERT(xdf->isFullNode(mxd));
  DCASSERT(ddf->isFullNode(result));
  DCASSERT(!ddf->isReducedNode(result));

  int min = MIN(ddf->getFullNodeSize(mdd), xdf->getFullNodeSize(mxd));
  for (int i = 0; i < min; i++)
  {
    int mddI = ddf->getFullNodeDownPtr(mdd, i);
    int mxdI = xdf->getFullNodeDownPtr(mxd, i);
    if (mddI != 0 && mxdI != 0) {
      DCASSERT(!xdf->isTerminalNode(mxdI));
      recFirePrime(mddI, mxdI, result);
    }
  }
}


void mdd_reachability_dfs::recFireFS(int mdd, int mxd, int result)
{
  DCASSERT(ddf->isFullNode(mdd));
  DCASSERT(xdf->isSparseNode(mxd));
  DCASSERT(ddf->isFullNode(result));
  DCASSERT(!ddf->isReducedNode(result));

  // go through the sparse node since it will have fewer entries
  int mddSize = ddf->getFullNodeSize(mdd);
  int mxdNnz = xdf->getSparseNodeSize(mxd);
  for (int i = 0; i < mxdNnz; i++)
  {
    int index = xdf->getSparseNodeIndex(mxd, i);
    if (index >= mddSize) break;

    int mxdI = xdf->getSparseNodeDownPtr(mxd, i);
    int mddI = ddf->getFullNodeDownPtr(mdd, index);

    if (mddI != 0) {
      DCASSERT(!xdf->isTerminalNode(mxdI));
      recFirePrime(mddI, mxdI, result);
    }
  }
}


void mdd_reachability_dfs::recFireSF(int mdd, int mxd, int result)
{
  DCASSERT(ddf->isSparseNode(mdd));
  DCASSERT(xdf->isFullNode(mxd));
  DCASSERT(ddf->isFullNode(result));
  DCASSERT(!ddf->isReducedNode(result));

  // go through the sparse node since it will have fewer entries
  int mxdSize = xdf->getFullNodeSize(mxd);
  int mddNnz = ddf->getSparseNodeSize(mdd);
  for (int i = 0; i < mddNnz; i++)
  {
    int index = ddf->getSparseNodeIndex(mdd, i);
    if (index >= mxdSize) break;

    int mddI = ddf->getSparseNodeDownPtr(mdd, i);
    int mxdI = xdf->getFullNodeDownPtr(mxd, index);

    if (mxdI != 0) {
      DCASSERT(!xdf->isTerminalNode(mxdI));
      recFirePrime(mddI, mxdI, result);
    }
  }
}


void mdd_reachability_dfs::recFireSS(int mdd, int mxd, int result)
{
  DCASSERT(ddf->isSparseNode(mdd));
  DCASSERT(xdf->isSparseNode(mxd));
  DCASSERT(ddf->isFullNode(result));
  DCASSERT(!ddf->isReducedNode(result));

  int mddNnz = ddf->getSparseNodeSize(mdd);
  int mxdNnz = xdf->getSparseNodeSize(mxd);

  // i refers to mxd index and j refers to mdd index
  int i = 0; int j = 0;
  for ( ; i < mxdNnz && j < mddNnz; )
  {
    int mxdIndex = xdf->getSparseNodeIndex(mxd, i);
    int mddIndex = ddf->getSparseNodeIndex(mdd, j);
    if (mddIndex < mxdIndex) {
      j++;
    }
    else if (mddIndex > mxdIndex) {
      i++;
    }
    else {
      // mddIndex == mxdIndex
      // index(mdd, j) == index(mxd, i)
      int mxdI = xdf->getSparseNodeDownPtr(mxd, i);
      int mddI = ddf->getSparseNodeDownPtr(mdd, j);

      DCASSERT(!xdf->isTerminalNode(mxdI));
      recFirePrime(mddI, mxdI, result);

      i++; j++;
    }
  }
}


void mdd_reachability_dfs::recFirePrime(int mddI, int mxdI, int result)
{
  DCASSERT(ddf->isFullNode(result));
  DCASSERT(!ddf->isReducedNode(result));

  if (xdf->isFullNode(mxdI)) {
    int mxdISize = xdf->getFullNodeSize(mxdI);
    for (int j = 0; j < mxdISize; j++)
    {
      int mxdIJ = xdf->getFullNodeDownPtr(mxdI, j);
      if (mxdIJ != 0) {
        int f = recFire(mddI, mxdIJ);
        if (f != 0) {
          // update result[j]
          int u = getMddUnion(ddf->getFullNodeDownPtr(result, j), f);
          ddf->setDownPtr(result, j, u);
          ddf->unlinkNode(u);
        }
        ddf->unlinkNode(f);
      }
    }
  }
  else {
    DCASSERT(xdf->isSparseNode(mxdI));
    int mxdINnz = xdf->getSparseNodeSize(mxdI);
    for (int j = 0; j < mxdINnz; j++)
    {
      int mxdIJ = xdf->getSparseNodeDownPtr(mxdI, j);
      int f = recFire(mddI, mxdIJ);
      if (f != 0) {
        // update result[mxdIJIndex]
        int mxdIJIndex = xdf->getSparseNodeIndex(mxdI, j);
        int u = getMddUnion(ddf->getFullNodeDownPtr(result, mxdIJIndex), f);
        ddf->setDownPtr(result, mxdIJIndex, u);
        ddf->unlinkNode(u);
      }
      ddf->unlinkNode(f);
    }
  }
}

#endif

#else


int mdd_reachability_dfs::saturate(int mdd)
{
#ifdef DEBUG_DFS
  printf("mdd: %d\n", mdd);
#endif

  // how does saturateHelper get called?
  // bottom-up i.e. call helper on children before calling helper for parent

  DCASSERT(ddf->isReducedNode(mdd));

  // terminal condition for recursion
  if (ddf->isTerminalNode(mdd)) return mdd;

  int k = ddf->getNodeLevel(mdd);      // level
  int sz = ddf->getLevelSize(k);       // size

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d\n", mdd, k, sz);
#endif

  std::vector<int> node(sz, 0);
  std::vector<int> mddDptrs;
  ddf->getDownPtrs(mdd, mddDptrs);

  std::vector<int>::iterator nodeIter = node.begin();
  std::vector<int>::iterator mddIter = mddDptrs.begin();
  for ( ; mddIter != mddDptrs.end(); ++nodeIter, ++mddIter)
  {
    if (*mddIter != 0) *nodeIter = saturate(*mddIter);
  }
  
  // call saturateHelper for n
#ifdef DEBUG_DFS
  printf("Calling saturate: level %d\n", k);
#endif
  saturateHelper(k, node);

  // reduce and return
  int n = ddf->createTempNode(k, node);
  n = ddf->reduceNode(n);

#ifdef DEBUG_DFS
  ddf->showNodeGraph(stdout, n);
#endif

  return n;
}


void mdd_reachability_dfs::saturateHelper(int mddLevel, std::vector<int>& mdd)
{
  DCASSERT(unsigned(ddf->getLevelSize(mddLevel)) == mdd.size());

  int mxd = splits[mddLevel];
  if (xdf->isTerminalNode(mxd)) return;

  std::vector<int> mxdDptrs;
  if (!xdf->getDownPtrs(mxd, mxdDptrs)) return;

  std::vector<bool> curr(mdd.size(), false);
  std::vector<bool> next(mdd.size(), true);
  bool repeat = true;

  while (repeat)
  {
    curr = next;
    fill(next.begin(), next.end(), false);
    repeat = false;

    // for each mxd[i1 != 0
    for (unsigned i = 0u; i < mxdDptrs.size(); i++)
    {
      if (mxdDptrs[i] == 0 || mdd[i] == 0 || !curr[i]) continue;
      DCASSERT(!xdf->isTerminalNode(mxdDptrs[i]));

      std::vector<int> mxdIDptrs;
      xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs);

      // for each mxd[i][j] != 0
      for (unsigned j = 0u; j < mxdIDptrs.size(); j++)
      {
        if (mxdIDptrs[j] == 0) continue;
        int f = recFire(mdd[i], mxdIDptrs[j]);
        if (f == 0) continue;
        int u = getMddUnion(mdd[j], f);
        ddf->unlinkNode(f);
        if (u != mdd[j]) {
          // update mdd[j] and mark for next iteration
          ddf->unlinkNode(mdd[j]);
          mdd[j] = u;
          if (j > i) {
            curr[j] = true;
          } else {
            next[j] = true;
            repeat = true;
          }
        } else {
          ddf->unlinkNode(u);
        }
      }
    }
  }
}


int mdd_reachability_dfs::recFire(int mdd, int mxd)
{
  DCASSERT(ddf->isReducedNode(mdd));
  DCASSERT(xdf->isReducedNode(mxd));

  if (mxd == -1) {
    ddf->linkNode(mdd);
    return mdd;
  }
  if (mxd == 0 || mdd == 0) return 0;

  int result = 0;
  if (findResult(owner, mdd, mxd, result)) {
    return result;
  }

  int mxdHeight = xdf->getNodeHeight(mxd);
  int mddHeight = ddf->getNodeHeight(mdd);
  int nodeHeight = MAX(mxdHeight, mddHeight);
  int nodeLevel = ed->getVariableWithHeight(nodeHeight);
  int newSize = ddf->getLevelSize(nodeLevel);
  std::vector<int> node(newSize, 0);

  if (mxdHeight < mddHeight) {
    std::vector<int> mddDptrs;
    ddf->getDownPtrs(mdd, mddDptrs);
    for (unsigned i = 0; i < mddDptrs.size(); i++)
    {
      if (mddDptrs[i] != 0) node[i] = recFire(mddDptrs[i], mxd);
    }
  } else if (mxdHeight > mddHeight) {
    std::vector<int> mxdDptrs;
    xdf->getDownPtrs(mxd, mxdDptrs);
    for (unsigned i = 0; i < mxdDptrs.size(); i++)
    {
      if (mxdDptrs[i] == 0) continue;
      std::vector<int> mxdIDptrs;
      xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs);
      for (unsigned j = 0; j < mxdIDptrs.size(); j++)
      {
        if (mxdIDptrs[j] == 0) continue;
        int f = recFire(mdd, mxdIDptrs[j]);
        if (f == 0) continue;
        int u = getMddUnion(node[j], f);
        ddf->unlinkNode(f);
        ddf->unlinkNode(node[j]);
        node[j] = u;
      }
    }
  } else {
    DCASSERT(mxdHeight == mddHeight);
    std::vector<int> mddDptrs;
    std::vector<int> mxdDptrs;
    ddf->getDownPtrs(mdd, mddDptrs);
    xdf->getDownPtrs(mxd, mxdDptrs);
    unsigned min = MIN(mddDptrs.size(), mxdDptrs.size());
    for (unsigned i = 0; i < min; i++)
    {
      if (mxdDptrs[i] == 0 || mddDptrs[i] == 0) continue;
      std::vector<int> mxdIDptrs;
      xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs);
      for (unsigned j = 0; j < mxdIDptrs.size(); j++)
      {
        if (mxdIDptrs[j] == 0) continue;
        int f = recFire(mddDptrs[i], mxdIDptrs[j]);
        if (f == 0) continue;
        int u = getMddUnion(node[j], f);
        ddf->unlinkNode(f);
        ddf->unlinkNode(node[j]);
        node[j] = u;
      }
    }
  }

  unsigned i = 0u;
  for ( ; i < node.size() && node[i] == 0; i++);
  if (i != node.size()) saturateHelper(nodeLevel, node);
  int n = ddf->createTempNode(nodeLevel, node);

  result = ddf->reduceNode(n);
  saveResult(owner, mdd, mxd, result);
  return result;
}

#endif


int mdd_reachability_dfs::getMddUnion(int a, int b)
{
  return mddUnion->compute(mddUnionOp, a, b);
}

int mdd_reachability_dfs::getMxdIntersection(int a, int b)
{
  return mxdIntersection->compute(mxdIntersectionOp, a, b);
}

int mdd_reachability_dfs::getMxdDifference(int a, int b)
{
  return mxdDifference->compute(mxdDifferenceOp, a, b);
}


// ----------------------- MTMDD Apply operation ----------------------


mtmdd_apply_operation::mtmdd_apply_operation()
{ }

mtmdd_apply_operation::~mtmdd_apply_operation() {}

compute_manager::error
mtmdd_apply_operation::typeCheck(const op_info* owner)
{
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] != owner->f[1] || owner->f[0] != owner->f[2])
    return compute_manager::FOREST_MISMATCH;
  if (owner->f[0]->isForRelations() ||
      owner->f[0]->getRangeType() == forest::BOOLEAN ||
      // the above allows MTMDDs with INTEGERs and REALs
      // the line below allows only INTEGERs
      // owner->f[0]->getRangeType() != forest::INTEGER ||
      owner->f[0]->getEdgeLabeling() != forest::MULTI_TERMINAL)
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] != owner->f[1] || owner->f[0] != owner->f[2])
    return compute_manager::FOREST_MISMATCH;
  if (owner->f[0]->isForRelations() ||
      owner->f[0]->getRangeType() == forest::BOOLEAN ||
      // the above allows MTMDDs with INTEGERs and REALs
      // the line below allows only INTEGERs
      // owner->f[0]->getRangeType() != forest::INTEGER ||
      owner->f[0]->getEdgeLabeling() != forest::MULTI_TERMINAL)
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] != owner->f[1])
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
  initScratch(smart_cast<const expert_domain*>(owner->f[0]->getDomain()));
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 2)
    return compute_manager::WRONG_NUMBER;
  // No assumption about what kind of forests these are, but they
  // have to be forests in the same domain
  if (owner->f[0]->getDomain() != owner->f[1]->getDomain())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool conversion_operation::isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nForests
  // data[i] <--> forest[i]
  // call isStale for each forest[i] and data[i]
  DCASSERT(owner->nForests == 2);
  return 
    getExpertForest(owner, 0)->isStale(data[0]) ||
    getExpertForest(owner, 1)->isStale(data[1]);
}


void
conversion_operation::discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nForests
  // data[i] <--> forest[i]
  // call uncacheNode for each forest[i] and data[i]
  DCASSERT(owner->nForests == 2);
  getExpertForest(owner, 0)->uncacheNode(data[0]);
  getExpertForest(owner, 1)->uncacheNode(data[1]);
}


void
conversion_operation::showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // TODO:
  // data[] is of size owner.nForests
  // data[i] <--> forest[i]
  // call showNode for each forest[i] and data[i]
  DCASSERT(owner->nForests == 2);
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0]->getDomain() != owner->f[1]->getDomain())
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0]->getDomain() != owner->f[1]->getDomain())
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0]->getDomain() != owner->f[1]->getDomain())
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0]->getDomain() != owner->f[1]->getDomain())
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0]->getDomain() != owner->f[1]->getDomain())
    return compute_manager::TYPE_MISMATCH;
  if (!getExpertForest(owner, 0)->isMtMdd())
    return compute_manager::TYPE_MISMATCH;
  if (owner->f[1]->getEdgeLabeling() == forest::MULTI_TERMINAL)
    // !MULTI_TERMINAL is same as EV+ or EV*
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mtmdd_to_evmdd::isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nForests + 1 (last int is edge value)
  // data[i] <--> forest[i]
  // call isStale for each forest[i] and data[i]
  DCASSERT(owner->nForests == 2);
  return 
    getExpertForest(owner, 0)->isStale(data[0]) ||
    getExpertForest(owner, 1)->isStale(data[1]);
}


void
mtmdd_to_evmdd::discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nForests + 1 (last int is edge value)
  // data[i] <--> forest[i]
  // call uncacheNode for each forest[i] and data[i]
  DCASSERT(owner->nForests == 2);
  getExpertForest(owner, 0)->uncacheNode(data[0]);
  getExpertForest(owner, 1)->uncacheNode(data[1]);
}


void
mtmdd_to_evmdd::showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // data[] is of size owner.nForests + 1 (last int is edge value)
  // data[i] <--> forest[i]
  // call showNode for each forest[i] and data[i]
  DCASSERT(owner->nForests == 2);
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 2)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0]->getDomain() != owner->f[1]->getDomain())
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
  int nVars = owner->f[0]->getDomain()->getNumVariables();
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] != owner->f[1] || owner->f[0] != owner->f[2])
    return compute_manager::FOREST_MISMATCH;
  if (!(owner->f[0]->isForRelations()) ||
      owner->f[0]->getRangeType() == forest::BOOLEAN ||
      // the above allows MTMXDs with INTEGERs and REALs
      // the line below allows only INTEGERs
      // owner->f[0]->getRangeType() != forest::INTEGER ||
      owner->f[0]->getEdgeLabeling() != forest::MULTI_TERMINAL)
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] != owner->f[1] || owner->f[0] != owner->f[2])
    return compute_manager::FOREST_MISMATCH;
  if (!(owner->f[0]->isForRelations()) ||
      owner->f[0]->getRangeType() == forest::BOOLEAN ||
      // the above allows MTMXDs with INTEGERs and REALs
      // the line below allows only INTEGERs
      // owner->f[0]->getRangeType() != forest::INTEGER ||
      owner->f[0]->getEdgeLabeling() != forest::MULTI_TERMINAL)
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



// ---------------------- MTMDD Post-Image Operation -------------------


mtmdd_post_image* mtmdd_post_image::getInstance()
{
  static mtmdd_post_image instance;
  return &instance;
}


mtmdd_post_image::mtmdd_post_image()
{ }


mtmdd_post_image::~mtmdd_post_image() {}


compute_manager::error
mtmdd_post_image::typeCheck(const op_info* owner)
{
  // op1 == MTMDD, op2 == MTMXD, op3 = MTMDD
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] == owner->f[1] || owner->f[0] != owner->f[2])
    return compute_manager::FOREST_MISMATCH;
  if (!getExpertForest(owner, 0)->isMtMdd() ||
      !getExpertForest(owner, 1)->isMtMxd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


int mtmdd_post_image::compute(op_info* owner, int mdd, int mxd)
{
  DCASSERT(owner->nForests == 3 && owner->f[0] == owner->f[2]);
  const int nOperands = 3;
  forest* forests[nOperands] = {owner->f[0], owner->f[0], owner->f[0]};
  op_info* plusOp = 
    smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
    getOpInfo(compute_manager::PLUS, forests, nOperands);
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


compute_manager::error
mtmdd_pre_image::typeCheck(const op_info* owner)
{
  // op1 == MTMDD, op2 == MTMXD, op3 = MTMDD
  if (owner == 0)
    return compute_manager::UNKNOWN_OPERATION;
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] == owner->f[1] || owner->f[0] != owner->f[2])
    return compute_manager::FOREST_MISMATCH;
  if (!getExpertForest(owner, 0)->isMtMdd() ||
      !getExpertForest(owner, 1)->isMtMxd())
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


int mtmdd_pre_image::compute(op_info* owner, int mdd, int mxd)
{
  DCASSERT(owner->nForests == 3 && owner->f[0] == owner->f[2]);
  const int nOperands = 3;
  forest* forests[nOperands] = {owner->f[0], owner->f[0], owner->f[0]};
  op_info* plusOp = 
    smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
    getOpInfo(compute_manager::PLUS, forests, nOperands);
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
  // data[] is of size owner.nForests * 2
  // data[2i] and data[2i+1] <--> forest[i]
  // call isStale for each forest[i] and data[2i]
  DCASSERT(owner->nForests == 3);
  return 
    getExpertForest(owner, 0)->isStale(data[0]) ||
    getExpertForest(owner, 1)->isStale(data[2]) ||
    getExpertForest(owner, 2)->isStale(data[4]);
}


void
evmdd_apply_operation::
discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nForests * 2
  // data[2i] and data[2i+1] <--> forest[i]
  // call uncacheNode for each forest[i] and data[2i]
  DCASSERT(owner->nForests == 3);
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
  // data[] is of size owner.nForests * 2
  // data[2i] and data[2i+1] <--> forest[i]
  DCASSERT(owner->nForests == 3);
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] != owner->f[1] || owner->f[0] != owner->f[2])
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
  DCASSERT(owner->f[0] == owner->f[1] && owner->f[1] == owner->f[2]);

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
  // data[] is of size owner.nForests * 2
  // data[2i] and data[2i+1] <--> forest[i]
  DCASSERT(owner->nForests == 3);
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
  if (owner->op == 0 || owner->f == 0 || owner->cc == 0)
    return compute_manager::TYPE_MISMATCH;
  if (owner->nForests != 3)
    return compute_manager::WRONG_NUMBER;
  if (owner->f[0] != owner->f[1] || owner->f[0] != owner->f[2])
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
  DCASSERT(owner->f[0] == owner->f[1] && owner->f[1] == owner->f[2]);

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
  DCASSERT(owner->f[0] == owner->f[1] && owner->f[1] == owner->f[2]);

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
  DCASSERT(owner->f[0] == owner->f[1] && owner->f[1] == owner->f[2]);

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


#if 0

// ---------------------- MTMDD Saturation-based Reachability -----------------


mtmdd_reachability_dfs* mtmdd_reachability_dfs::getInstance()
{
  static mtmdd_reachability_dfs instance;
  return &instance;
}


mtmdd_reachability_dfs::mtmdd_reachability_dfs()
{ }


mtmdd_reachability_dfs::~mtmdd_reachability_dfs() {}


int mtmdd_reachability_dfs::compute(op_info* owner, int mtmdd, int mtmxd)
{
  DCASSERT(owner->nForests == 3 && owner->f[0] == owner->f[2]);

  // Initialize class members and helper operations
  initialize(owner);

  // Depth-first reachability analysis (Saturation)

#ifdef DEBUG_DFS
  printf("Consolidated Next-State Function:\n");
  xdf->showNodeGraph(stdout, mtmxd);
  printf("\n");

  printf("Initial State:\n");
  ddf->showNodeGraph(stdout, mtmdd);
  printf("\n");
#endif

  // Split the next-state function: each level has its own next-state function
  // The nsf is stored into the vector splits
  splitMxd(mtmxd);

#ifdef DEBUG_DFS
  printf("Split Next-State Function:\n");
  for (int i = splits.size() - 1; i >= 0; i--)
  {
    printf("Level %d, Node %d\n", i, splits[i]);
    xdf->showNodeGraph(stdout, splits[i]);
    printf("\n");
  }

  fflush(stdout);
#endif

  // Saturate the node
  int result = saturate(mtmdd);

  // clear pointers to dd nodes, class members and operation pointers
  clear();

  return result;
}


void mtmdd_reachability_dfs::initialize(op_info* o)
{
  // set up aliases
  owner = o;
  ecm = smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager());
  assert(ecm != 0);
  ddf = getExpertForest(owner, 0);
  assert(ddf != 0);
  xdf = getExpertForest(owner, 1);
  assert(xdf != 0);
  DCASSERT(ddf->getDomain() == xdf->getDomain());
  ed = smart_cast<expert_domain*>(ddf->useDomain());
  assert(ed != 0);

  // set up mtmdd operation: union
  const int nOperands = 3;
  forest* forests[nOperands] = {owner->f[0], owner->f[0], owner->f[0]};

  mtmddUnionOp = ecm->getOpInfo(compute_manager::PLUS, forests, nOperands);
  assert(mtmddUnionOp != 0);
  mtmddUnion = smart_cast<mtmdd_union*>(mtmddUnionOp->op);
  assert(mtmddUnion != 0);

  // set up mtmxd operations: intersection and difference
  forests[0] = owner->f[1]; forests[1] = owner->f[1]; forests[2] = owner->f[1];

  // Equivalent of MDD Intersection is MTMDD Min
  mtmxdIntersectionOp =
    ecm->getOpInfo(compute_manager::MIN, forests, nOperands);
  assert(mtmxdIntersectionOp != 0);
  mtmxdIntersection = smart_cast<mtmxd_intersection*>(mtmxdIntersectionOp->op);
  assert(mtmxdIntersection != 0);

  mtmxdDifferenceOp =
    ecm->getOpInfo(compute_manager::MINUS, forests, nOperands);
  assert(mtmxdDifferenceOp != 0);
  mtmxdDifference = smart_cast<mtmxd_difference*>(mtmxdDifferenceOp->op);
  assert(mtmxdDifference != 0);

  // Intialize the scratch 2-D vector (i.e. make it the correct size)
  int nLevels = ed->getTopVariable() + 1;
  scratch.clear();
  scratch.resize(nLevels);
  for (unsigned height = ed->getNumVariables(); height > 0; --height)
  {
    int lh = ed->getVariableWithHeight(height);
    int sz = ed->getVariableBound(lh);
    scratch[lh].resize(sz, 0);
  }

  // Initialize the splits vector
  splits.resize(nLevels, 0);

  // Initialize the boolean 2-D vectors
  curr.resize(nLevels);
  for (int i = 0; i < nLevels; ++i)
  {
    curr[i].resize(scratch[i].size(), false);
  }
  next = curr;
}


void mtmdd_reachability_dfs::clear()
{
  // clear pointer to dd nodes
  for (unsigned i = 0u; i < splits.size(); i++) xdf->unlinkNode(splits[i]);

  // clear class members and pointers to operations
  splits.clear();
  scratch.clear();
  curr.clear();
  next.clear();
  owner = 0;
  ecm = 0;
  ddf = 0;
  xdf = 0;
  ed = 0;
  mtmddUnionOp = 0;
  mtmxdIntersectionOp = 0;
  mtmxdDifferenceOp = 0;
  mtmddUnion = 0;
  mtmxdIntersection = 0;
  mtmxdDifference = 0;
}


// split is used to split a mtmxd for the saturation algorithm
void mtmdd_reachability_dfs::splitMxd(int mtmxd)
{
  DCASSERT(xdf != 0);

  int zeroNode;
  if (xdf->getRangeType() == forest::REAL) {
    zeroNode = xdf->getTerminalNode(0.0);
  } else {
    zeroNode = xdf->getTerminalNode(0);
  }

  // find intersection for all mtmxd[i][i]
  // -- if mtmxd is smaller than max level size, then some mtmxd[i] is zero,
  //    therefore intersection is 0.
  //    -- if node is sparse, then there is some mtmxd[i] = 0,
  //       therefore intersection is 0.
  //    -- if node is full, if size < max level size, intersection is 0.
  // 
  // if intersection == 0, add mtmxd to level_mtmxd[level], return.
  // otherwise, 
  // -- create mtmxdSize nodes at primed level with a copy of corresponding
  //    mtmxd[i].
  // -- for each new_mtmxd[i][i],
  //    -- mtmxd[i][i] = mtmxd[i][i] - intersection 
  // -- set new_mtmxd[i] after reducing the primed level nodes
  //    note that new_mtmxd will never be 0 since mtmxd is not an identity node
  //
  // add new_mtmxd to level_mtmxd[level]
  // 
  // repeat the above for intersection
  //
  int level = 0;
  int intersection = zeroNode;
  int mtmxdSize = 0;
  int mtmxdI = zeroNode;

  xdf->linkNode(mtmxd);

  while (!xdf->isTerminalNode(mtmxd)) {
    level = xdf->getNodeLevel(mtmxd);
    DCASSERT(level > 0); // we only deal with unprimed levels

    // Find intersection for all mtmxd[i][i]
    // Note: only do this if mtmxd is a full node; when it is sparse, some
    // mtmxd[i] is 0 therefore the intersection will always be 0 (zeroNode).
    intersection = zeroNode;
    if (xdf->isFullNode(mtmxd)) {
      mtmxdSize = xdf->getFullNodeSize(mtmxd);
      if (mtmxdSize == xdf->getLevelSize(level)) {
        // for all i, mtmxd[i] != 0
        intersection = zeroNode;
        bool first = true;
        for (int i = 0; i < mtmxdSize; ++i)
        {
          mtmxdI = xdf->getFullNodeDownPtr(mtmxd, i);

          // If any mtmxd[i] is a terminal (according to Identity Reduced rules)
          // it must be node 0, and mtmxd[i][i] is also 0. Therefore,
          // the intersection is 0. So check for this condition, and break
          // out of the loop it true.

          // if mtmxdI is a terminal node it must be a 0 (zeroNode)
          DCASSERT((xdf->isTerminalNode(mtmxdI) && mtmxdI == zeroNode) ||
              !xdf->isTerminalNode(mtmxdI));

          int mtmxdII = zeroNode;

          if (!xdf->isTerminalNode(mtmxdI)) {
            if (xdf->isFullNode(mtmxdI)) {
              if (xdf->getFullNodeSize(mtmxdI) > i)
                mtmxdII = xdf->getFullNodeDownPtr(mtmxdI, i);
            } else {
              DCASSERT(xdf->isSparseNode(mtmxdI));
              // search for ith index
              int found = -1;
              int mtmxdINnz = xdf->getSparseNodeSize(mtmxdI);

              if (mtmxdINnz > 8) {
                // binary search
                int start = 0;
                int stop = mtmxdINnz - 1;

                while (start < stop) {
                  int mid = (start + stop) / 2;
                  int midIndex = xdf->getSparseNodeIndex(mtmxdI, mid);
                  if (midIndex < i) {
                    start = (mid == start)? mid + 1: mid;
                  } else {
                    stop = mid;
                  }
                }

                assert(start == stop);
                if (xdf->getSparseNodeIndex(mtmxdI, start) == i) {
                  found = start;
                }
              }
              else {
                // linear search
                for (int j = 0; j < mtmxdINnz; ++j)
                {
                  if (xdf->getSparseNodeIndex(mtmxdI, j) == i) {
                    found = j;
                    break;
                  }
                }
              }

              if (found != -1)
                mtmxdII = xdf->getSparseNodeDownPtr(mtmxdI, found);
            }
          }

          if (!first) {
            int temp = getMxdIntersection(intersection, mtmxdII);
            xdf->unlinkNode(intersection);
            intersection = temp;
          } else {
            first = false;
            xdf->linkNode(mtmxdII);
            xdf->unlinkNode(intersection);
            intersection = mtmxdII;
          }
#ifdef DEBUG_DFS
          printf("intersection: %d level: %d\n",
              intersection, xdf->getNodeLevel(intersection));
#endif
          if (intersection == zeroNode) break;
        }
      }
    }

    DCASSERT(splits[level] == zeroNode);

    DCASSERT(intersection == zeroNode ||
        xdf->getNodeLevel(mtmxd) > xdf->getNodeLevel(intersection));

    if (intersection != zeroNode) {
      splits[level] = getMxdDifference(mtmxd, intersection);
        // mtmxdDifference->compute(mtmxdDifferenceOp, mtmxd, intersection);
    } else {
      splits[level] = mtmxd;
      xdf->linkNode(mtmxd);
    }

    // intersection becomes the mtmxd for the next iteration
    xdf->unlinkNode(mtmxd);
    mtmxd = intersection;
  }

  DCASSERT(xdf->isTerminalNode(mtmxd));
  xdf->unlinkNode(mtmxd);
}


int mtmdd_reachability_dfs::saturate(int mtmdd)
{
#ifdef DEBUG_DFS
  printf("mtmdd: %d\n", mtmdd);
#endif

  // how does saturateHelper get called?
  // bottom-up i.e. call helper on children before calling helper for parent

  DCASSERT(ddf->isReducedNode(mtmdd));

  // terminal condition for recursion
  if (ddf->isTerminalNode(mtmdd)) return mtmdd;

  int k = ddf->getNodeLevel(mtmdd);   // level
  int sz = ddf->getLevelSize(k);      // size

#ifdef DEBUG_DFS
  printf("mtmdd: %d, level: %d, size: %d\n", mtmdd, k, sz);
#endif

  std::vector<int> node(sz, 0);
  std::vector<int> mtmddDptrs;
  ddf->getDownPtrs(mtmdd, mtmddDptrs);

  int zeroNode = (ddf->getRangeType() == forest::REAL)?
      ddf->getTerminalNode(0.0): ddf->getTerminalNode(0);

  std::vector<int>::iterator nodeIter = node.begin();
  std::vector<int>::iterator mtmddIter = mtmddDptrs.begin();
  for ( ; mtmddIter != mtmddDptrs.end(); ++nodeIter, ++mtmddIter)
  {
    if (*mtmddIter != zeroNode) *nodeIter = saturate(*mtmddIter);
  }
  
  // call saturateHelper for n
#ifdef DEBUG_DFS
  printf("Calling saturate: level %d\n", k);
#endif
  saturateHelper(k, node);

  // reduce and return
  int n = ddf->createTempNode(k, node);
  n = ddf->reduceNode(n);

#ifdef DEBUG_DFS
  ddf->showNodeGraph(stdout, n);
#endif

  return n;
}


void mtmdd_reachability_dfs::replaceTerminalNode(int mtmdd, int term)
{
  DCASSERT(0 == (ddf->getRangeType() == forest::REAL?
    ddf->getTerminalNode(0.0): ddf->getTerminalNode(0)));
  DCASSERT(mtmdd != 0);

  if (ddf->isTerminalNode(mtmdd)) return term;

  std::vector<int> dptrs;
  ddf->getDownPtrs(mtmdd, dptrs);
  for (unsigned i = 0; i < dptrs.size(); i++)
  {
    if (dptrs[i] != 0) dptrs[i] = replaceTerminalNode(dptrs[i], term);
  }
  int result = ddf->createTempNode(ddf->getNodeLevel(mtmdd), dptrs);
  return ddf->reduceNode(result);
}


void mtmdd_reachability_dfs::saturateHelper(int mtmddLevel,
    std::vector<int>& mtmdd)
{
  DCASSERT(unsigned(ddf->getLevelSize(mtmddLevel)) == mtmdd.size());
  DCASSERT(0 == (ddf->getRangeType() == forest::REAL?
    ddf->getTerminalNode(0.0): ddf->getTerminalNode(0)));

  int mtmxd = splits[mtmddLevel];
  if (xdf->isTerminalNode(mtmxd)) {
    if (mtmxd == 0) return;

    // states don't change but the terminal values will change to mtmxd
    for (unsigned i = 0u; i < mtmdd.size(); i++)
    {
      if (mtmdd[i] != 0) {
        int temp = replaceTerminalNode(mtmdd[i], mtmxd);
        ddf->unlinkNode(mtmdd[i]);
        mtmdd[i] = temp;
      }
    }
    return;
  }
  
  std::vector<int> mtmxdDptrs;
  if (!xdf->getDownPtrs(mtmxd, mtmxdDptrs)) return;

  std::vector<bool> curr(mtmdd.size(), false);
  std::vector<bool> next(mtmdd.size(), true);
  bool repeat = true;

  while (repeat)
  {
    curr = next;
    fill(next.begin(), next.end(), false);
    repeat = false;

    // for each mtmxd[i1 != 0
    for (unsigned i = 0u; i < mtmxdDptrs.size(); i++)
    {
      if (mtmxdDptrs[i] == 0 || mtmdd[i] == 0 || !curr[i]) continue;
      DCASSERT(!xdf->isTerminalNode(mtmxdDptrs[i]));

      std::vector<int> mtmxdIDptrs;
      xdf->getDownPtrs(mtmxdDptrs[i], mtmxdIDptrs);

      // for each mtmxd[i][j] != 0
      for (unsigned j = 0u; j < mtmxdIDptrs.size(); j++)
      {
        if (mtmxdIDptrs[j] == 0) continue;
        int f = recFire(mtmdd[i], mtmxdIDptrs[j]);
        if (f == 0) continue;
        int u = getMddUnion(mtmdd[j], f);
        ddf->unlinkNode(f);
        if (u != mtmdd[j]) {
          // update mtmdd[j] and mark for next iteration
          ddf->unlinkNode(mtmdd[j]);
          mtmdd[j] = u;
          if (j > i) {
            curr[j] = true;
          } else {
            next[j] = true;
            repeat = true;
          }
        } else {
          ddf->unlinkNode(u);
        }
      }
    }
  }
}


int mtmdd_reachability_dfs::recFire(int mtmdd, int mtmxd)
{
  DCASSERT(ddf->isReducedNode(mtmdd));
  DCASSERT(xdf->isReducedNode(mtmxd));

  if (mtmxd == -1) {
    ddf->linkNode(mtmdd);
    return mtmdd;
  }
  if (mtmxd == 0 || mtmdd == 0) return 0;

  int result = 0;
  if (findResult(owner, mtmdd, mtmxd, result)) {
    return result;
  }

  int mtmxdHeight = xdf->getNodeHeight(mtmxd);
  int mtmddHeight = ddf->getNodeHeight(mtmdd);
  int nodeHeight = MAX(mtmxdHeight, mtmddHeight);
  int nodeLevel = ed->getVariableWithHeight(nodeHeight);
  int newSize = ddf->getLevelSize(nodeLevel);
  std::vector<int> node(newSize, 0);

  if (mtmxdHeight < mtmddHeight) {
    std::vector<int> mtmddDptrs;
    ddf->getDownPtrs(mtmdd, mtmddDptrs);
    for (unsigned i = 0; i < mtmddDptrs.size(); i++)
    {
      if (mtmddDptrs[i] != 0) node[i] = recFire(mtmddDptrs[i], mtmxd);
    }
  } else if (mtmxdHeight > mtmddHeight) {
    std::vector<int> mtmxdDptrs;
    xdf->getDownPtrs(mtmxd, mtmxdDptrs);
    for (unsigned i = 0; i < mtmxdDptrs.size(); i++)
    {
      if (mtmxdDptrs[i] == 0) continue;
      std::vector<int> mtmxdIDptrs;
      xdf->getDownPtrs(mtmxdDptrs[i], mtmxdIDptrs);
      for (unsigned j = 0; j < mtmxdIDptrs.size(); j++)
      {
        if (mtmxdIDptrs[j] == 0) continue;
        int f = recFire(mtmdd, mtmxdIDptrs[j]);
        if (f == 0) continue;
        int u = getMddUnion(node[j], f);
        ddf->unlinkNode(f);
        ddf->unlinkNode(node[j]);
        node[j] = u;
      }
    }
  } else {
    DCASSERT(mtmxdHeight == mtmddHeight);
    std::vector<int> mtmddDptrs;
    std::vector<int> mtmxdDptrs;
    ddf->getDownPtrs(mtmdd, mtmddDptrs);
    xdf->getDownPtrs(mtmxd, mtmxdDptrs);
    unsigned min = MIN(mtmddDptrs.size(), mtmxdDptrs.size());
    for (unsigned i = 0; i < min; i++)
    {
      if (mtmxdDptrs[i] == 0 || mtmddDptrs[i] == 0) continue;
      std::vector<int> mtmxdIDptrs;
      xdf->getDownPtrs(mtmxdDptrs[i], mtmxdIDptrs);
      for (unsigned j = 0; j < mtmxdIDptrs.size(); j++)
      {
        if (mtmxdIDptrs[j] == 0) continue;
        int f = recFire(mtmddDptrs[i], mtmxdIDptrs[j]);
        if (f == 0) continue;
        int u = getMddUnion(node[j], f);
        ddf->unlinkNode(f);
        ddf->unlinkNode(node[j]);
        node[j] = u;
      }
    }
  }

  unsigned i = 0u;
  for ( ; i < node.size() && node[i] == 0; i++);
  if (i != node.size()) saturateHelper(nodeLevel, node);
  int n = ddf->createTempNode(nodeLevel, node);

  result = ddf->reduceNode(n);
  saveResult(owner, mtmdd, mtmxd, result);
  return result;
}


int mtmdd_reachability_dfs::getMddUnion(int a, int b)
{
  return mtmddUnion->compute(mtmddUnionOp, a, b);
}

int mtmdd_reachability_dfs::getMxdIntersection(int a, int b)
{
  return mtmxdIntersection->compute(mtmxdIntersectionOp, a, b);
}

int mtmdd_reachability_dfs::getMxdDifference(int a, int b)
{
  return mtmxdDifference->compute(mtmxdDifferenceOp, a, b);
}

#endif

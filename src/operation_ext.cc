
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



#include "../src/operation_ext.h"
#include "../src/compute_cache.h"

//#define IGNORE_TERMS 0
//#define IGNORE_INCOUNT 2

inline expert_forest* getExpertForest(op_info* op, int index) {
  return smart_cast<expert_forest*>(op->f[index]);
}

inline const expert_forest* getExpertForest(const op_info* op, int index) {
  return smart_cast<const expert_forest*>(op->f[index]);
}

mdd_apply_operation::mdd_apply_operation(const char *name, bool commutative)
: operation(2, 1, name), commutative(commutative) {}

mdd_apply_operation::~mdd_apply_operation() {}

compute_manager::error
mdd_apply_operation::typeCheck(const op_info* owner)
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
      owner->f[0]->getRangeType() != forest::BOOLEAN ||
      owner->f[0]->getEdgeLabeling() != forest::MULTI_TERMINAL)
    return compute_manager::TYPE_MISMATCH;
  return compute_manager::SUCCESS;
}


bool mdd_apply_operation::isEntryStale(const op_info* owner, const int* data)
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
mdd_apply_operation::discardEntry(op_info* owner, const int* data)
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
mdd_apply_operation::showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // TODO:
  // data[] is of size owner.nForests
  // data[i] <--> forest[i]
  // call showNode for each forest[i] and data[i]
  DCASSERT(owner->nForests == 3);
  fprintf(strm, "[%s %d %d %d]",
      owner->op->getName(), data[0], data[1], data[2]);
#if 0
  getExpertForest(owner, 0)->showNode(strm, data[0]);
  getExpertForest(owner, 1)->showNode(strm, data[1]);
  getExpertForest(owner, 2)->showNode(strm, data[2]);
#endif
}


bool
mdd_apply_operation::findResult(op_info* owner, int a, int b, int& c)
{
#if 1

  static int key[2];

#ifdef IGNORE_TERMS
  if (getExpertForest(owner, 0)->isTerminalNode(a) ||
      getExpertForest(owner, 1)->isTerminalNode(b))
    return false;
#endif

  // create cache entry
  if (commutative && a > b) {
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

#else

#ifdef IGNORE_TERMS
  if (getExpertForest(owner, 0)->isTerminalNode(a) ||
      getExpertForest(owner, 1)->isTerminalNode(b))
    return false;
#endif

  // create cache entry

  if (commutative && a > b) {
    if (owner->cc->find(owner, b, a, c)) {
      getExpertForest(owner, 2)->linkNode(c);
      return true;
    }
    return false;
  }
  else {
    if (owner->cc->find(owner, a, b, c)) {
      getExpertForest(owner, 2)->linkNode(c);
      return true;
    }
    return false;
  }

#endif
}


void
mdd_apply_operation::saveResult(op_info* owner, int a, int b, int c)
{
#if 1

  static int cacheEntry[3];

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
  if (commutative && a > b) {
    // sort the entry in ascending order
    cacheEntry[0] = b; cacheEntry[1] = a;
  } else {
    cacheEntry[0] = a; cacheEntry[1] = b;
  }
  cacheEntry[2] = c;

#if 0
#ifdef DEVELOPMENT_CODE
  assert(!findResult(owner, cacheEntry[0], cacheEntry[1], cacheEntry[2]));
#endif
#endif

  getExpertForest(owner, 0)->cacheNode(cacheEntry[0]);
  getExpertForest(owner, 1)->cacheNode(cacheEntry[1]);
  getExpertForest(owner, 2)->cacheNode(cacheEntry[2]);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));

#else

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

  getExpertForest(owner, 0)->cacheNode(a);
  getExpertForest(owner, 1)->cacheNode(b);
  getExpertForest(owner, 2)->cacheNode(c);

  if (commutative && a > b)
    owner->cc->add(owner, b, a, c);
  else
    owner->cc->add(owner, a, b, c);

#endif
}


compute_manager::error
mdd_apply_operation::compute(op_info* owner, dd_edge** operands)
{
  if (operands == 0) return compute_manager::TYPE_MISMATCH;
  // compute(owner, dd_edge, dd_edge, dd_edge) checks for owner == 0
  return compute(owner, *operands[0], *operands[1], *operands[2]);
}

compute_manager::error
mdd_apply_operation::compute(op_info* owner, const dd_edge& a, dd_edge& b)
{
  return compute_manager::TYPE_MISMATCH;
}

compute_manager::error
mdd_apply_operation::compute(op_info* owner, const dd_edge& a,
    const dd_edge& b, dd_edge& c)
{
  if (owner == 0) return compute_manager::TYPE_MISMATCH;
  int result = compute(owner, a.getNode(), b.getNode());
  c.set(result, 0, getExpertForest(owner, 2)->getNodeLevel(result));
  return compute_manager::SUCCESS;
}


void mdd_apply_operation::expandA(op_info* owner, expert_forest* expertForest,
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


void mdd_apply_operation::expandB(op_info* owner, expert_forest* expertForest,
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


#if 0

// Two-step compute
// - result[i] = a[i]
// - result[i] = result[i] op b[i]
// This makes for a simpler but less efficient implementation.
int mdd_apply_operation::compute(op_info* owner, int a, int b)
{
  int result = 0;
  if (checkTerminals(owner, a, b, result))
    return result;
  if (findResult(owner, a, b, result))
    return result;

  // expand nodes
  // 0. initialize result
  // 1. copy node a to node result
  // 2. do operation between contents of result and node b

  // 0. initialize result
  const int aLevel = getExpertForest(owner, 0)->getNodeLevel(a);
  const int bLevel = getExpertForest(owner, 1)->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = getExpertForest(owner, 2)->getLevelSize(resultLevel);
  result = getExpertForest(owner, 2)->createTempNode(resultLevel, resultSize);

  // 1. copy node a to node result
  if (aLevel < resultLevel) {
    // all down pointers of result point to node a
    for (int i = 0; i < resultSize; ++i)
      getExpertForest(owner, 2)->setDownPtr(result, i, a);
  }
  else if (getExpertForest(owner, 0)->isFullNode(a)) {
    // a is a full-node
    const int aSize = getExpertForest(owner, 0)->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    for (int i = 0; i < aSize; ++i)
      getExpertForest(owner, 2)->setDownPtr(result, i,
          getExpertForest(owner, 0)->getFullNodeDownPtr(a, i));
  }
  else {
    // a is a sparse-node
    const int aSize =
      getExpertForest(owner, 0)->getSparseNodeSize(a);
    for (int i = 0; i < aSize; ++i)
      getExpertForest(owner, 2)->setDownPtr(result,
          getExpertForest(owner, 0)->getSparseNodeIndex(a, i),
          getExpertForest(owner, 0)->getSparseNodeDownPtr(a, i));
  }

  // 2. do operation between contents of result and node b
  if (bLevel < resultLevel) {
    for (int i = 0; i < resultSize; ++i)
    {
      int tempResult = compute(owner,
          getExpertForest(owner, 2)->getFullNodeDownPtr(result, i), b);
      getExpertForest(owner, 2)->setDownPtr(result, i, tempResult);
      getExpertForest(owner, 2)->unlinkNode(tempResult);
    }
  }
  else if (getExpertForest(owner, 1)->isFullNode(b)) {
    // b is a full-node
    const int bSize = getExpertForest(owner, 1)->getFullNodeSize(b);
    DCASSERT(bSize <= resultSize);
    for (int i = 0; i < bSize; ++i)
    {
      int tempResult = compute(owner,
          getExpertForest(owner, 2)->getFullNodeDownPtr(result, i),
          getExpertForest(owner, 1)->getFullNodeDownPtr(b, i));
      getExpertForest(owner, 2)->setDownPtr(result, i, tempResult);
      getExpertForest(owner, 2)->unlinkNode(tempResult);
    }
    for (int i = bSize; i < resultSize; ++i)
    {
      int tempResult = compute(owner,
          getExpertForest(owner, 2)->getFullNodeDownPtr(result, i), 0);
      getExpertForest(owner, 2)->setDownPtr(result, i, tempResult);
      getExpertForest(owner, 2)->unlinkNode(tempResult);
    }
  }
  else {
    // b is a sparse-node
    const int bSize = getExpertForest(owner, 1)->getSparseNodeSize(b);
    DCASSERT(getExpertForest(owner, 1)->getSparseNodeIndex(b, bSize - 1) <=
        resultSize);
    // j goes through every index (like a full-node index pointer)
    int j = 0;
    for (int i = 0; i < bSize; ++i, ++j)
    {
      // sparse-nodes skip indices which represent downpointer 0
      // call compute of those skipped indices
      for (int index = getExpertForest(owner, 1)->getSparseNodeIndex(b, i);
          j < index; ++j)
      {
        int tempResult = compute(owner,
            getExpertForest(owner, 2)->getFullNodeDownPtr(result, j), 0);
        getExpertForest(owner, 2)->setDownPtr(result, j, tempResult);
        getExpertForest(owner, 2)->unlinkNode(tempResult);
      }
      // done with skipped indices; deal with the next sparse node index
      DCASSERT(j == getExpertForest(owner, 1)->getSparseNodeIndex(b, i));
      int tempResult = compute(owner,
          getExpertForest(owner, 2)->getFullNodeDownPtr(result, j),
          getExpertForest(owner, 1)->getSparseNodeDownPtr(b, i));
      getExpertForest(owner, 2)->setDownPtr(result, j, tempResult);
      getExpertForest(owner, 2)->unlinkNode(tempResult);
    }
    DCASSERT(j ==
        getExpertForest(owner, 1)->getSparseNodeIndex(b, bSize - 1) + 1);
    for ( ; j < resultSize; ++j)
    {
      int tempResult = compute(owner,
          getExpertForest(owner, 2)->getFullNodeDownPtr(result, j), 0);
      getExpertForest(owner, 2)->setDownPtr(result, j, tempResult);
      getExpertForest(owner, 2)->unlinkNode(tempResult);
    }
  }

  // save result in compute cache and return it
  result = getExpertForest(owner, 2)->reduceNode(result);
  saveResult(owner, a, b, result);
  return result;
}

#else

// Single-step compute
// result[i] = a[i] op b[i]

int mdd_apply_operation::compute(op_info* owner, int a, int b)
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

  // save result in compute cache and return it
  result = expertForest->reduceNode(result);
  saveResult(owner, a, b, result);
  return result;
}

#endif


void mdd_apply_operation::fullFull (op_info* owner, expert_forest* mddNm,
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


void mdd_apply_operation::sparseSparse (op_info* owner, expert_forest* mddNm,
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


void mdd_apply_operation::fullSparse (op_info* owner, expert_forest* mddNm,
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


void mdd_apply_operation::sparseFull (op_info* owner, expert_forest* mddNm,
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


mdd_union* mdd_union::getInstance()
{
  static mdd_union instance("MDD Union");
  return &instance;
}


mdd_union::mdd_union(const char *name)
: mdd_apply_operation(name, true)
{ }


mdd_union::~mdd_union() {}


bool
mdd_union::checkTerminals(op_info* op, int a, int b, int& c)
{
#if 1
  if (a == 0 || b == 0 || a == b) {
    c = a | b;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }
  if (a == -1 || b == -1) {
    c = -1;
    return true;
  }

  return false;
#else
  if (-1 == a) {
    c = -1;
    return true;
  }

  if (-1 == b) {
    c = -1;
    return true;
  }

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
#endif
}

#if 1

// BUG!! in single-step compute

#if 0

// Two-step compute
// - result[i] = a[i]
// - result[i] = result[i] union b[i]
// This makes for a simpler but less efficient implementation.

int mdd_union::compute(op_info* owner, int a, int b)
{
  // std::cerr << "2-step mdd_union.\n";
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
  result = expertForest->createTempNode(resultLevel, resultSize, false);

  // 1. copy node a to node result
  if (aLevel < resultLevel) {
    // all down pointers of result point to node a
    for (int i = 0; i < resultSize; ++i)
      expertForest->setDownPtrWoUnlink(result, i, a);
  }
  else if (expertForest->isFullNode(a)) {
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
    const int lastIndex = expertForest->getSparseNodeIndex(a, aSize - 1);
    int i = 0;  // for traversing sparse node indices
    int index = expertForest->getSparseNodeIndex(a, i);
    int j = 0;
    for ( ; j < lastIndex; ++j) {
      if (j == index) {
        expertForest->setDownPtrWoUnlink(result, j,
            expertForest->getSparseNodeDownPtr(a, i));
        ++i;
        index = expertForest->getSparseNodeIndex(a, i);
      }
      else {
        expertForest->setDownPtrWoUnlink(result, j, 0);
      }
    }
    // deal with last index
    DCASSERT(j == lastIndex);
    expertForest->setDownPtrWoUnlink(result, j,
        expertForest->getSparseNodeDownPtr(a, i));
    for (++j; j < resultSize; ++j)
      expertForest->setDownPtrWoUnlink(result, j, 0);
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

#endif


#if 1

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


#endif

#endif


// ------------------------------------------------------------------


// ----------------------- MDD Intersection --------------------------------


mdd_intersection* mdd_intersection::getInstance()
{
  static mdd_intersection instance("MDD Intersection");
  return &instance;
}


mdd_intersection::mdd_intersection(const char *name)
: mdd_apply_operation(name, true)
{ }


mdd_intersection::~mdd_intersection() {}


bool
mdd_intersection::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (0 == a || 0 == b) {
    c = 0;
    return true;
  }

  if (-1 == a) {
    c = b;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (-1 == b) {
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


#if 1

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
  static mdd_difference instance("MDD Difference");
  return &instance;
}


mdd_difference::mdd_difference(const char *name)
: mdd_apply_operation(name, false)
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


mxd_apply_operation::mxd_apply_operation(const char *name, bool commutative)
: mdd_apply_operation(name, commutative) {}


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


int mxd_apply_operation::compute(op_info* owner, int a, int b)
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
        int tempResult = compute(owner,
            expertForest->getFullNodeDownPtr(result, i),
            expertForest->getFullNodeDownPtr(b, i));
        expertForest->setDownPtr(result, i, tempResult);
        expertForest->unlinkNode(tempResult);
      }
      for (int i = bSize; i < resultSize; ++i)
      {
        int tempResult =
          compute(owner, expertForest->getFullNodeDownPtr(result, i), 0);
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
          int tempResult = compute(owner,
              expertForest->getFullNodeDownPtr(result, j), 0);
          expertForest->setDownPtr(result, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        // done with skipped indices; deal with the next sparse node index
        DCASSERT(j == expertForest->getSparseNodeIndex(b, i));
        int tempResult = compute(owner,
            expertForest->getFullNodeDownPtr(result, j),
            expertForest->getSparseNodeDownPtr(b, i));
        expertForest->setDownPtr(result, j, tempResult);
        expertForest->unlinkNode(tempResult);
      }
      DCASSERT(j == expertForest->getSparseNodeIndex(b, bSize - 1) + 1);
      for ( ; j < resultSize; ++j)
      {
        int tempResult = compute(owner,
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

  int zeroB = compute(owner, 0, b);

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    int i = 0;
    for ( ; i < aSize; ++i)
    {
      int tempResult =
        compute(owner, expertForest->getFullNodeDownPtr(a, i), b);
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
        expertForest->setDownPtr(result, i, zeroB);
      }
      DCASSERT(i == expertForest->getSparseNodeIndex(a, aIndex));
      int tempResult =
        compute(owner, expertForest->getSparseNodeDownPtr(a, aIndex), b);
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

  int aZero = compute(owner, a, 0);

  if (expertForest->isFullNode(b)) {
    const int bSize = expertForest->getFullNodeSize(b);
    int i = 0;
    for ( ; i < bSize; ++i)
    {
      int tempResult =
        compute(owner, a, expertForest->getFullNodeDownPtr(b, i));
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
        expertForest->setDownPtr(result, i, aZero);
      }
      DCASSERT(i == expertForest->getSparseNodeIndex(b, bIndex));
      int tempResult =
        compute(owner, a, expertForest->getSparseNodeDownPtr(b, bIndex));
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
  int zeroOpZero = compute(owner, 0, 0);
  int zeroOpB = compute(owner, 0, b);

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    for (int i = 0; i < aSize; ++i)
    {
      int iA = expertForest->getFullNodeDownPtr(a, i);
      int iResult = expertForest->createTempNode(-resultLevel, resultSize);
      if (iA == 0) {
        for (int j = 0; j < resultSize; ++j)
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
              compute(owner, ijA, b): compute(owner, ijA, 0);
            expertForest->setDownPtr(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < resultSize; ++j)
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
            compute(owner, ijA, b): compute(owner, ijA, 0);
          expertForest->setDownPtr(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < resultSize; ++j)
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
      int iResult = expertForest->createTempNode(-resultLevel, resultSize);
      for (int j = 0; j < resultSize; ++j)
        expertForest->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtr(result, i, iResult);
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
        int iResult = expertForest->createTempNode(-resultLevel, resultSize);
        for (int j = 0; j < resultSize; ++j)
          expertForest->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        iResult = expertForest->reduceNode(iResult);
        expertForest->setDownPtr(result, i, iResult);
        expertForest->unlinkNode(iResult);
      }

      DCASSERT(i == expertForest->getSparseNodeIndex(a, aIndex));
      int iA = expertForest->getSparseNodeDownPtr(a, aIndex);
      int iResult = expertForest->createTempNode(-resultLevel, resultSize);
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
              compute(owner, ijA, b): compute(owner, ijA, 0);
            expertForest->setDownPtr(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < resultSize; ++j)
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
            compute(owner, ijA, b): compute(owner, ijA, 0);
          expertForest->setDownPtr(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < resultSize; ++j)
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
      int iResult = expertForest->createTempNode(-resultLevel, resultSize);
      for (int j = 0; j < resultSize; ++j)
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
  int zeroOpZero = compute(owner, 0, 0);
  int aOpZero = compute(owner, a, 0);

  if (expertForest->isFullNode(b)) {
    const int bSize = expertForest->getFullNodeSize(b);
    for (int i = 0; i < bSize; ++i)
    {
      int iB = expertForest->getFullNodeDownPtr(b, i);
      int iResult = expertForest->createTempNode(-resultLevel, resultSize);
      if (iB == 0) {
        for (int j = 0; j < resultSize; ++j)
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
              compute(owner, a, ijB): compute(owner, 0, ijB);
            expertForest->setDownPtr(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < resultSize; ++j)
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
            compute(owner, a, ijB): compute(owner, 0, ijB);
          expertForest->setDownPtr(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < resultSize; ++j)
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
      int iResult = expertForest->createTempNode(-resultLevel, resultSize);
      for (int j = 0; j < resultSize; ++j)
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
        int iResult = expertForest->createTempNode(-resultLevel, resultSize);
        for (int j = 0; j < resultSize; ++j)
          expertForest->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        iResult = expertForest->reduceNode(iResult);
        expertForest->setDownPtr(result, i, iResult);
        expertForest->unlinkNode(iResult);
      }

      DCASSERT(i == expertForest->getSparseNodeIndex(b, bIndex));
      int iB = expertForest->getSparseNodeDownPtr(b, bIndex);
      int iResult = expertForest->createTempNode(-resultLevel, resultSize);
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
              compute(owner, a, ijB): compute(owner, 0, ijB);
            expertForest->setDownPtr(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < resultSize; ++j)
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
            compute(owner, a, ijB): compute(owner, 0, ijB);
          expertForest->setDownPtr(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < resultSize; ++j)
        {
          expertForest->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        }
      }
      iResult = expertForest->reduceNode(iResult);
      expertForest->setDownPtr(result, i, iResult);
      expertForest->unlinkNode(iResult);
    }

    // TODO: can optimize when zeroOpZero == aOpZero == 0
    DCASSERT(i == expertForest->getSparseNodeIndex(b, bSize - 1) + 1);
    for ( ; i < resultSize; ++i)
    {
      int iResult = expertForest->createTempNode(-resultLevel, resultSize);
      for (int j = 0; j < resultSize; ++j)
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


mxd_alt_apply_operation::mxd_alt_apply_operation(const char *name,
    bool commutative) : operation(3, 1, name), commutative(commutative) {}


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
  if (commutative && a > b) {
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
  if (commutative && a > b) {
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

  DCASSERT(aLevel == resultLevel || bLevel == resultLevel);

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
        int tempResult = compute(owner, nextLevel,
            expertForest->getFullNodeDownPtr(result, i),
            expertForest->getFullNodeDownPtr(b, i));
        expertForest->setDownPtr(result, i, tempResult);
        expertForest->unlinkNode(tempResult);
      }
      for (int i = bSize; i < resultSize; ++i)
      {
        int tempResult = compute(owner, nextLevel,
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
          int tempResult = compute(owner, nextLevel,
              expertForest->getFullNodeDownPtr(result, j), 0);
          expertForest->setDownPtr(result, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        // done with skipped indices; deal with the next sparse node index
        DCASSERT(j == expertForest->getSparseNodeIndex(b, i));
        int tempResult = compute(owner, nextLevel,
            expertForest->getFullNodeDownPtr(result, j),
            expertForest->getSparseNodeDownPtr(b, i));
        expertForest->setDownPtr(result, j, tempResult);
        expertForest->unlinkNode(tempResult);
      }
      DCASSERT(j == expertForest->getSparseNodeIndex(b, bSize - 1) + 1);
      for ( ; j < resultSize; ++j)
      {
        int tempResult = compute(owner, nextLevel,
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
    int zeroZeroAtOneLevelBelow = compute(owner, resultLevel-1, 0, 0);
    int aBAtOneLevelBelow = compute(owner, resultLevel-1, a, b);
    for (int i = 0; i < resultSize; ++i)
    {
      // create primed node at -resultLevel
      int p = expertForest->createTempNode(-resultLevel, resultSize, false);
      for (int j = 0; j < resultSize; ++j)
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
    int zeroZeroAtOneLevelBelow = compute(owner, -resultLevel-1, 0, 0);
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
  int zeroB = compute(owner, -resultLevel-1, 0, b);

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    int i = 0;
    for ( ; i < aSize; ++i)
    {
      int tempResult = compute(owner, -resultLevel-1,
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
      int tempResult = compute(owner, -resultLevel-1,
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
  int aZero = compute(owner, -resultLevel-1, a, 0);

  if (expertForest->isFullNode(b)) {
    const int bSize = expertForest->getFullNodeSize(b);
    int i = 0;
    for ( ; i < bSize; ++i)
    {
      int tempResult = compute(owner, -resultLevel-1,
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
      int tempResult = compute(owner, -resultLevel-1,
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
  int zeroOpZero = compute(owner, resultLevel-1, 0, 0);
  int zeroOpB = compute(owner, resultLevel-1, 0, b);

  if (expertForest->isFullNode(a)) {
    const int aSize = expertForest->getFullNodeSize(a);
    for (int i = 0; i < aSize; ++i)
    {
      int iA = expertForest->getFullNodeDownPtr(a, i);
      int iResult =
        expertForest->createTempNode(-resultLevel, resultSize, false);
      if (iA == 0) {
        for (int j = 0; j < resultSize; ++j)
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? zeroOpB: zeroOpZero);
      }
      else if (expertForest->isFullNode(iA)) {
        const int iASize = expertForest->getFullNodeSize(iA);
        for (int j = 0; j < iASize; ++j)
        {
          int ijA = expertForest->getFullNodeDownPtr(iA, j);
          if (ijA == 0) {
            expertForest->setDownPtrWoUnlink(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              compute(owner, resultLevel-1, ijA, b):
              compute(owner, resultLevel-1, ijA, 0);
            expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < resultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? zeroOpB: zeroOpZero);
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
            expertForest->setDownPtrWoUnlink(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iA, k));
          int ijA = expertForest->getSparseNodeDownPtr(iA, k);
          int tempResult = (i == j)?
            compute(owner, resultLevel-1, ijA, b):
            compute(owner, resultLevel-1, ijA, 0);
          expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < resultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? zeroOpB: zeroOpZero);
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
        expertForest->createTempNode(-resultLevel, resultSize, false);
      for (int j = 0; j < resultSize; ++j)
        expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? zeroOpB: zeroOpZero);
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
          expertForest->createTempNode(-resultLevel, resultSize, false);
        for (int j = 0; j < resultSize; ++j)
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        iResult = expertForest->reduceNode(iResult);
        expertForest->setDownPtrWoUnlink(result, i, iResult);
        expertForest->unlinkNode(iResult);
      }

      DCASSERT(i == expertForest->getSparseNodeIndex(a, aIndex));
      int iA = expertForest->getSparseNodeDownPtr(a, aIndex);
      int iResult =
        expertForest->createTempNode(-resultLevel, resultSize, false);
      DCASSERT(iA != 0 && iA != -1);
      if (expertForest->isFullNode(iA)) {
        const int iASize = expertForest->getFullNodeSize(iA);
        for (int j = 0; j < iASize; ++j)
        {
          int ijA = expertForest->getFullNodeDownPtr(iA, j);
          if (ijA == 0) {
            expertForest->setDownPtrWoUnlink(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              compute(owner, resultLevel-1, ijA, b):
              compute(owner, resultLevel-1, ijA, 0);
            expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < resultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? zeroOpB: zeroOpZero);
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
            expertForest->setDownPtrWoUnlink(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iA, k));
          int ijA = expertForest->getSparseNodeDownPtr(iA, k);
          int tempResult = (i == j)?
            compute(owner, resultLevel-1, ijA, b):
            compute(owner, resultLevel-1, ijA, 0);
          expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < resultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? zeroOpB: zeroOpZero);
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
        expertForest->createTempNode(-resultLevel, resultSize, false);
      for (int j = 0; j < resultSize; ++j)
        expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? zeroOpB: zeroOpZero);
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
  int zeroOpZero = compute(owner, resultLevel-1, 0, 0);
  int aOpZero = compute(owner, resultLevel-1, a, 0);

  if (expertForest->isFullNode(b)) {
    const int bSize = expertForest->getFullNodeSize(b);
    for (int i = 0; i < bSize; ++i)
    {
      int iB = expertForest->getFullNodeDownPtr(b, i);
      int iResult =
        expertForest->createTempNode(-resultLevel, resultSize, false);
      if (iB == 0) {
        for (int j = 0; j < resultSize; ++j)
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? aOpZero: zeroOpZero);
      }
      else if (expertForest->isFullNode(iB)) {
        const int iBSize = expertForest->getFullNodeSize(iB);
        for (int j = 0; j < iBSize; ++j)
        {
          int ijB = expertForest->getFullNodeDownPtr(iB, j);
          if (ijB == 0) {
            expertForest->setDownPtrWoUnlink(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              compute(owner, resultLevel-1, a, ijB):
              compute(owner, resultLevel-1, 0, ijB);
            expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < resultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? aOpZero: zeroOpZero);
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
            expertForest->setDownPtrWoUnlink(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iB, k));
          int ijB = expertForest->getSparseNodeDownPtr(iB, k);
          int tempResult = (i == j)?
            compute(owner, resultLevel-1, a, ijB):
            compute(owner, resultLevel-1, 0, ijB);
          expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < resultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? aOpZero: zeroOpZero);
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
        expertForest->createTempNode(-resultLevel, resultSize, false);
      for (int j = 0; j < resultSize; ++j)
        expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? aOpZero: zeroOpZero);
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
          expertForest->createTempNode(-resultLevel, resultSize, false);
        for (int j = 0; j < resultSize; ++j)
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? aOpZero: zeroOpZero);
        iResult = expertForest->reduceNode(iResult);
        expertForest->setDownPtrWoUnlink(result, i, iResult);
        expertForest->unlinkNode(iResult);
      }

      DCASSERT(i == expertForest->getSparseNodeIndex(b, bIndex));
      int iB = expertForest->getSparseNodeDownPtr(b, bIndex);
      int iResult =
        expertForest->createTempNode(-resultLevel, resultSize, false);
      DCASSERT(iB != 0 && iB != -1);
      if (expertForest->isFullNode(iB)) {
        const int iBSize = expertForest->getFullNodeSize(iB);
        for (int j = 0; j < iBSize; ++j)
        {
          int ijB = expertForest->getFullNodeDownPtr(iB, j);
          if (ijB == 0) {
            expertForest->setDownPtrWoUnlink(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              compute(owner, resultLevel-1, a, ijB):
              compute(owner, resultLevel-1, 0, ijB);
            expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
            expertForest->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < resultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? aOpZero: zeroOpZero);
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
            expertForest->setDownPtrWoUnlink(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == expertForest->getSparseNodeIndex(iB, k));
          int ijB = expertForest->getSparseNodeDownPtr(iB, k);
          int tempResult = (i == j)?
            compute(owner, resultLevel-1, a, ijB):
            compute(owner, resultLevel-1, 0, ijB);
          expertForest->setDownPtrWoUnlink(iResult, j, tempResult);
          expertForest->unlinkNode(tempResult);
        }
        DCASSERT(j == expertForest->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < resultSize; ++j)
        {
          expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? aOpZero: zeroOpZero);
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
        expertForest->createTempNode(-resultLevel, resultSize, false);
      for (int j = 0; j < resultSize; ++j)
        expertForest->setDownPtrWoUnlink(iResult, j, (i == j)? aOpZero: zeroOpZero);
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
  static mxd_union instance("MXD Union");
  return &instance;
}


mxd_union::mxd_union(const char *name)
: mxd_apply_operation(name, true)
{ }


mxd_union::~mxd_union() {}


bool
mxd_union::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (-1 == a || -1 == b) {
    c = -1;
    return true;
  }

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
  static mxd_intersection instance("MXD Intersection");
  return &instance;
}


mxd_intersection::mxd_intersection(const char *name)
: mxd_apply_operation(name, true)
{ }


mxd_intersection::~mxd_intersection() {}


bool
mxd_intersection::checkTerminals(op_info* op, int a, int b, int& c)
{
  if (0 == a || 0 == b) {
    c = 0;
    return true;
  }

  if (-1 == a) {
    c = b;
    getExpertForest(op, 0)->linkNode(c);
    return true;
  }

  if (-1 == b) {
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


// ----------------------- MXD Difference --------------------------------


mxd_difference* mxd_difference::getInstance()
{
  static mxd_difference instance("MXD Difference");
  return &instance;
}


mxd_difference::mxd_difference(const char *name)
: mxd_apply_operation(name, false)
{ }


mxd_difference::~mxd_difference() {}


bool
mxd_difference::checkTerminals(op_info* op, int a, int b, int& c)
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


// ---------------------- MDD MXD Image Operation -------------------


mdd_mxd_image_operation::mdd_mxd_image_operation(const char* name)
: mdd_apply_operation(name, false) {}


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
  static mdd_post_image instance("MDD Post-Image");
  return &instance;
}


mdd_post_image::mdd_post_image(const char *name)
: mdd_mxd_image_operation(name)
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
  return compute(owner, unionOp, mdd, mxd);
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
  static mdd_pre_image instance("MDD Pre-Image");
  return &instance;
}


mdd_pre_image::mdd_pre_image(const char *name)
: mdd_post_image(name)
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
  static mdd_reachability_bfs instance("MDD Reachability BFS");
  return &instance;
}


mdd_reachability_bfs::mdd_reachability_bfs(const char *name)
: mdd_mxd_image_operation(name)
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
  dd_edge postImage(mddNm);

  while(prevReachableStates != reachableStates)
  {
    prevReachableStates = reachableStates;
    // printf("\nPost-Image (mdd:%d, mxd:%d): ",
    //    reachableStates.getNode(), nsf.getNode());
    ecm->apply(postImageOp, reachableStates, nsf, postImage);
    // printf("%d\n", postImage.getNode());
    // postImage.show(stdout, 2);
    // printf("\nUnion (mdd:%d, mdd:%d): ",
    //    reachableStates.getNode(), postImage.getNode());
    ecm->apply(unionOp, reachableStates, postImage, reachableStates);
    // printf("%d\n", reachableStates.getNode());
  }

  int result = postImage.getNode();
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
  int prevPostImage = 0;

  mddNm->linkNode(reachableStates);
  mddNm->linkNode(postImage);

  do {
    prevPostImage = postImage;
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


#if 0

// split is used to split a mxd for the saturation algorithm
vector<int> split(expert_forest* f, int mxd)
{
  HERE: use this with saturation

  DCASSERT(f != 0);

  expert_compute_manager* ecm =
    smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager());
  const int nForests = 3;
  forest* forests[nForests] = {f, f, f};
  op_info* intersectionOp = ecm->getOpInfo(compute_manager::INTERSECTION,
      forests, nForests);
  op_info* differenceOp = ecm->getOpInfo(compute_manager::DIFFERENCE,
      forests, nForests);

  int trueNode = f->getTerminalNode(true);
  int falseNode = f->getTerminalNode(false);
  vector<int> splits(f->getDomain()->getTopVariable() + 1, falseNode);

  mxd_union* unionOpPtr = smart_cast<mxd_union*>(unionOp->op);
  DCASSERT(unionOpPtr != 0);
  mxd_intersection* intersectionOpPtr =
    smart_cast<mxd_intersection*>(intersectionOp->op);
  DCASSERT(intersectionOpPtr != 0);
  mxd_difference* differenceOpPtr =
    smart_cast<mxd_difference*>(differenceOp->op);
  DCASSERT(differenceOpPtr != 0);

  // find intersection for all mxd[i][i]
  // -- if mxd is smaller than max level size, then some mxd[i] is zero,
  //    therefore intersection is 0.
  //    -- if node is sparse, then there is some mxd[i] = 0,
  //       therefore intersection is 0.
  //    -- if node is full, if size < max level size, intersection is 0.
  // 
  // if intersection == 0, add mxd to level_mxd[level], return.
  // otherwise, 
  // -- create mxd_size nodes at primed level with a copy of corresponding
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
  int temp = falseNode;

  f->linkNode(mxd);

  while (!f->isTerminalNode(mxd)) {
    level = f->getNodeLevel(mxd);
    DCASSERT(level != 0);
    DCASSERT(level > 0); // we only deal with unprimed levels

    // find intersection for all mxd[i][i]
    intersection = falseNode;
    if (f->isFullNode(mxd)) {
      mxdSize = f->getFullNodeSize(mxd);
      if (mxdSize == f->getLevelSize(level)) {
        // for all i, mxd[i] != 0
        intersection = trueNode;
        for (int i = 0; i < mxdSize; ++i)
        {
          mxdI = f->getFullNodeDownPtr(mxd, i);

          // If any mxd[i] is a terminal (according to Identity Reduced rules)
          // it must be node 0, and mxd[i][i] is also 0. Therefore,
          // the intersection is 0. So check for this condition, and break
          // out of the loop it true.
          if (f->isTerminalNode(mxdI)) {
            DCASSERT(mxdI == falseNode);
            f->unlinkNode(intersection);
            intersection = falseNode;
            break;
          }

          temp = intersectionOpPtr->compute(intersectionOp,
              intersection, f->getFullNodeDownPtr(mxdI, i));
          f->unlinkNode(intersection);
          intersection = temp;
        }
      }
    }

    DCASSERT(splits[level] == 0);

    DCASSERT(intersection == falseNode ||
        f->getNodeLevel(mxd) > f->getNodeLevel(intersection));

    splits[level] = differenceOpPtr->compute(differenceOp, mxd, intersection);

    // intersection becomes the mxd for the next iteration
    f->unlinkNode(mxd);
    mxd = intersection;
  }

  DCASSERT(f->isTerminalNode(mxd));
  f->unlinkNode(mxd);

  return splits;
}

#endif


// ----------------------- MTMDD Apply operation ----------------------


mtmdd_apply_operation::mtmdd_apply_operation(const char *name, bool comm)
: mdd_apply_operation(name, comm)
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
  static mtmdd_divide instance("MTMDD Divide");
  return &instance;
}


mtmdd_divide::mtmdd_divide(const char *name)
: mtmdd_apply_operation(name, false)
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
  static mtmdd_multiply instance("MTMDD Multiply");
  return &instance;
}


mtmdd_multiply::mtmdd_multiply(const char *name)
: mtmdd_apply_operation(name, true)
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
  static mtmdd_minus instance("MTMDD Minus");
  return &instance;
}


mtmdd_minus::mtmdd_minus(const char *name)
: mtmdd_apply_operation(name, false)
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
  static mtmdd_plus instance("MTMDD Plus");
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

mtmdd_plus::mtmdd_plus(const char *name)
: mtmdd_apply_operation(name, true)
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


mtmdd_plus::mtmdd_plus(const char *name)
: mdd_union(name)
{ }


#endif

// ----------------------- MTMDD Min ----------------------------


mtmdd_min* mtmdd_min::getInstance()
{
  static mtmdd_min instance("MTMDD Min");
  return &instance;
}


mtmdd_min::mtmdd_min(const char *name)
: mtmdd_apply_operation(name, true)
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
  static mtmdd_max instance("MTMDD Max");
  return &instance;
}


mtmdd_max::mtmdd_max(const char *name)
: mtmdd_apply_operation(name, true)
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
  static mtmdd_or_min instance("MTMDD Or-Min");
  return &instance;
}


mtmdd_or_min::mtmdd_or_min(const char *name)
: mtmdd_apply_operation(name, true)
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
  static mtmdd_or_max instance("MTMDD Or-Max");
  return &instance;
}


mtmdd_or_max::mtmdd_or_max(const char *name)
: mtmdd_apply_operation(name, true)
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
  static mtmdd_and_min instance("MTMDD And-Min");
  return &instance;
}


mtmdd_and_min::mtmdd_and_min(const char *name)
: mtmdd_apply_operation(name, true)
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
  static mtmdd_and_max instance("MTMDD And-Max");
  return &instance;
}


mtmdd_and_max::mtmdd_and_max(const char *name)
: mtmdd_apply_operation(name, true)
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
  static mtmdd_less_than instance("MTMDD Less-Than");
  return &instance;
}


mtmdd_less_than::mtmdd_less_than(const char *name)
: mtmdd_apply_operation(name, false)
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
  static mtmdd_less_than_equal instance("MTMDD Less-Than-Or-Equal-To");
  return &instance;
}


mtmdd_less_than_equal::mtmdd_less_than_equal(const char *name)
: mtmdd_apply_operation(name, false)
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
  static mtmdd_greater_than instance("MTMDD Greater-Than");
  return &instance;
}


mtmdd_greater_than::mtmdd_greater_than(const char *name)
: mtmdd_apply_operation(name, false)
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
  static mtmdd_greater_than_equal instance("MTMDD Greater-Than-Or-Equal-To");
  return &instance;
}


mtmdd_greater_than_equal::mtmdd_greater_than_equal(const char *name)
: mtmdd_apply_operation(name, false)
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
  static mtmdd_equal instance("MTMDD Equal-To");
  return &instance;
}


mtmdd_equal::mtmdd_equal(const char *name)
: mtmdd_apply_operation(name, true)
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
  static mtmdd_not_equal instance("MTMDD Not-Equal-To");
  return &instance;
}


mtmdd_not_equal::mtmdd_not_equal(const char *name)
: mtmdd_apply_operation(name, true)
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


mtmdd_to_mdd_apply_operation::mtmdd_to_mdd_apply_operation(const char *name,
    bool comm)
: mtmdd_apply_operation(name, comm), scratch(0)
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
  static mtmdd_to_mdd_less_than instance("MTMDD-MDD Less-Than");
  return &instance;
}


mtmdd_to_mdd_less_than::mtmdd_to_mdd_less_than(const char *name)
: mtmdd_to_mdd_apply_operation(name, false)
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
    instance("MTMDD-MDD Less-Than-Or-Equal-To");
  return &instance;
}


mtmdd_to_mdd_less_than_equal::mtmdd_to_mdd_less_than_equal(const char *name)
: mtmdd_to_mdd_apply_operation(name, false)
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
  static mtmdd_to_mdd_greater_than instance("MTMDD-MDD Greater-Than");
  return &instance;
}


mtmdd_to_mdd_greater_than::mtmdd_to_mdd_greater_than(const char *name)
: mtmdd_to_mdd_apply_operation(name, false)
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
    instance("MTMDD-MDD Greater-Than-Or-Equal-To");
  return &instance;
}


mtmdd_to_mdd_greater_than_equal
::mtmdd_to_mdd_greater_than_equal(const char *name)
: mtmdd_to_mdd_apply_operation(name, false)
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
  static mtmdd_to_mdd_equal instance("MTMDD-MDD Equal-To");
  return &instance;
}


mtmdd_to_mdd_equal::mtmdd_to_mdd_equal(const char *name)
: mtmdd_to_mdd_apply_operation(name, true)
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
  static mtmdd_to_mdd_not_equal instance("MTMDD-MDD Not-Equal-To");
  return &instance;
}


mtmdd_to_mdd_not_equal::mtmdd_to_mdd_not_equal(const char *name)
: mtmdd_to_mdd_apply_operation(name, true)
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


conversion_operation::conversion_operation(const char *name)
: operation(1, 1, name)
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
  static mtmdd_to_mdd instance("Convert MTMDD to MDD");
  return &instance;
}


mtmdd_to_mdd::mtmdd_to_mdd(const char *name)
: conversion_operation(name)
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
  static mdd_to_mtmdd instance("Convert MDD to MTMDD");
  return &instance;
}


mdd_to_mtmdd::mdd_to_mtmdd(const char *name)
: conversion_operation(name)
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
  static mtmxd_to_mxd instance("Convert MTMXD to MXD");
  return &instance;
}


mtmxd_to_mxd::mtmxd_to_mxd(const char *name)
: conversion_operation(name)
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
  static mxd_to_mtmxd instance("Convert MXD to MTMXD");
  return &instance;
}


mxd_to_mtmxd::mxd_to_mtmxd(const char *name)
: conversion_operation(name)
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
  static mtmdd_to_evmdd instance("Convert MTMDD to EVMDD");
  return &instance;
}


mtmdd_to_evmdd::mtmdd_to_evmdd(const char *name)
: operation(1, 2, name)
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
  static int cacheEntry[2];

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
  static int cacheEntry[2];

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


// ----------------------- MTMXD Apply operation ----------------------


mtmxd_apply_operation::mtmxd_apply_operation(const char *name, bool comm)
: mxd_apply_operation(name, comm)
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
  static mtmxd_multiply instance("MTMXD Multiply");
  return &instance;
}


mtmxd_multiply::mtmxd_multiply(const char *name)
: mtmxd_apply_operation(name, true)
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
  static mtmxd_minus instance("MTMXD Minus");
  return &instance;
}


mtmxd_minus::mtmxd_minus(const char *name)
: mtmxd_apply_operation(name, false)
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
  static mtmxd_plus instance("MTMXD Plus");
  return &instance;
}


mtmxd_plus::mtmxd_plus(const char *name)
: mtmxd_apply_operation(name, true)
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
  static mtmxd_min instance("MTMXD Min");
  return &instance;
}


mtmxd_min::mtmxd_min(const char *name)
: mtmxd_apply_operation(name, true)
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
  static mtmxd_max instance("MTMXD Max");
  return &instance;
}


mtmxd_max::mtmxd_max(const char *name)
: mtmxd_apply_operation(name, true)
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
  static mtmxd_less_than instance("MTMXD Less-Than");
  return &instance;
}


mtmxd_less_than::mtmxd_less_than(const char *name)
: mtmxd_apply_operation(name, false)
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
  static mtmxd_greater_than instance("MTMXD Greater-Than");
  return &instance;
}


mtmxd_greater_than::mtmxd_greater_than(const char *name)
: mtmxd_apply_operation(name, false)
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
  static mtmxd_not_equal instance("MTMXD Not-Equal-To");
  return &instance;
}


mtmxd_not_equal::mtmxd_not_equal(const char *name)
: mtmxd_apply_operation(name, true)
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
mtmxd_alt_apply_operation(const char *name, bool comm)
: mxd_alt_apply_operation(name, comm)
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
  static mtmxd_divide instance("MTMXD Divide");
  return &instance;
}


mtmxd_divide::mtmxd_divide(const char *name)
: mtmxd_alt_apply_operation(name, false)
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
  static mtmxd_less_than_equal instance("MTMXD Less-Than-Or-Equal-To");
  return &instance;
}


mtmxd_less_than_equal::mtmxd_less_than_equal(const char *name)
: mtmxd_alt_apply_operation(name, false)
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
  static mtmxd_greater_than_equal instance("MTMXD Greater-Than-Or-Equal-To");
  return &instance;
}


mtmxd_greater_than_equal::mtmxd_greater_than_equal(const char *name)
: mtmxd_alt_apply_operation(name, false)
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
  static mtmxd_equal instance("MTMXD Equal-To");
  return &instance;
}


mtmxd_equal::mtmxd_equal(const char *name)
: mtmxd_alt_apply_operation(name, true)
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
  static mtmdd_post_image instance("MTMDD Post-Image");
  return &instance;
}


mtmdd_post_image::mtmdd_post_image(const char *name)
: mdd_post_image(name)
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
  static mtmdd_pre_image instance("MTMDD Pre-Image");
  return &instance;
}


mtmdd_pre_image::mtmdd_pre_image(const char *name)
: mdd_pre_image(name)
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




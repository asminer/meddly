
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
#include "apply_base.h"

// ******************************************************************
// *                                                                *
// *                   generic_binary_mdd methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_mdd::generic_binary_mdd(const binary_opname* code,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : binary_operation(code, 2, 1, arg1, arg2, res)
{
  can_commute = false;
}

MEDDLY::generic_binary_mdd::~generic_binary_mdd()
{
}

bool MEDDLY::generic_binary_mdd::isStaleEntry(const int* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}

void MEDDLY::generic_binary_mdd::discardEntry(const int* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void
MEDDLY::generic_binary_mdd ::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d, %d): %d]", getName(), data[0], data[1], data[2]);
}

void MEDDLY::generic_binary_mdd::compute(const dd_edge &a, const dd_edge &b, 
  dd_edge &c)
{
  int cnode = compute(a.getNode(), b.getNode());
  c.set(cnode, 0, resF->getNodeLevel(cnode));
}


int MEDDLY::generic_binary_mdd::compute(int a, int b)
{
  int result = 0;
  if (checkTerminals(a, b, result))
    return result;
  if (findResult(a, b, result))
    return result;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = resF->getLevelSize(resultLevel);

#if 0

  result = resF->createTempNode(resultLevel, resultSize, false);

  if (aLevel < resultLevel) {
    // expand b
    // result[i] = a op b[i]
    expandB(a, b, result, resultSize);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    expandA(a, b, result, resultSize);
  }
  else {
    DCASSERT(aLevel == bLevel);
    // expand a and b
    // result[i] = a[i] op b[i]

    if (arg1F->isFullNode(a)) {
      if (arg2F->isFullNode(b)) {
        fullFull(a, b, result, resultSize);
      }
      else {
        DCASSERT(arg2F->isSparseNode(b));
        fullSparse(a, b, result, resultSize);
      }
    }
    else {
      DCASSERT(arg1F->isSparseNode(a));
      if (arg2F->isFullNode(b)) {
        sparseFull(a, b, result, resultSize);
      }
      else {
        DCASSERT(arg2F->isSparseNode(b));
        sparseSparse(a, b, result, resultSize);
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
    int aZero = compute(a, 0);
    arg2F->getDownPtrs(b, B);
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    for ( ; iterB != B.end(); iterB++, iterC++)
    {
      if (*iterB == 0) {
        *iterC = aZero;
        resF->linkNode(aZero);
      } else {
        *iterC = compute(a, *iterB);
      }
    }
    resF->unlinkNode(aZero);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    int zeroB = compute(0, b);

    arg1F->getDownPtrs(a, A);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterC = C.begin();
    for ( ; iterA != A.end(); iterA++, iterC++)
    {
      if (*iterA == 0) {
        *iterC = zeroB;
        resF->linkNode(zeroB);
      } else {
        *iterC = compute(*iterA, b);
      }
    }
    resF->unlinkNode(zeroB);
  }
  else {
    // expand both a and b
    // result[i] = a[i] op b[i]
    int zeroZero = compute(0, 0);

    arg1F->getDownPtrs(a, A);
    arg2F->getDownPtrs(b, B);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    for ( ; iterA != A.end(); iterA++, iterB++, iterC++)
    {
      if (*iterA == 0 && *iterB == 0) {
        *iterC = zeroZero;
        resF->linkNode(zeroZero);
      } else {
        *iterC = compute(*iterA, *iterB); 
      }
    }

    resF->unlinkNode(zeroZero);
  }

  result = resF->createTempNode(resultLevel, C);

#endif

  // save result in compute cache and return it

#if 0
  printf("reduce(%d): ", result);
  result = resF->reduceNode(result);
  printf("%d  [", result);
  for (unsigned i = 0; i < C.size(); i++ )
  {
    printf("%d ", C[i]);
  }
  printf("]\n");
#else
  result = resF->reduceNode(result);
#endif

  saveResult(a, b, result);
#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d) = %d\n", getName(), a, b, result);
#endif
  return result;
}

void
MEDDLY::generic_binary_mdd::
expandA(int a, int b, int result, int resultSize)
{
  // fill in result node and return
  // result[i] = compute(a[i], b)

  if (arg1F->isFullNode(a)) {
    const int aSize = arg1F->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    int zeroB = compute(0, b);
    for (int i = 0; i < aSize; ++i)
    {
      int iA = arg1F->getFullNodeDownPtr(a, i);
      if (iA == 0) {
        resF->setDownPtrWoUnlink(result, i, zeroB);
      }
      else {
        int tempResult = compute(iA, b);
        resF->setDownPtrWoUnlink(result, i, tempResult);
        resF->unlinkNode(tempResult);
      }
    }
    for (int i = aSize; i < resultSize; ++i)
    {
      resF->setDownPtrWoUnlink(result, i, zeroB);
    }
    resF->unlinkNode(zeroB);
  } else {
    DCASSERT(arg1F->isSparseNode(a));
    const int aSize = arg1F->getSparseNodeSize(a);
    int zeroB = compute(0, b);
    // j goes through every index (like a full-node index pointer)
    int j = 0;
    for (int i = 0; i < aSize; ++i, ++j)
    {
      // sparse-nodes skip indices which represent downpointer 0
      for (int index = arg1F->getSparseNodeIndex(a, i); j < index; ++j)
      {
        resF->setDownPtrWoUnlink(result, j, zeroB);
      }
      // done with skipped indices; deal with the next sparse node index
      DCASSERT(j == arg1F->getSparseNodeIndex(a, i));
      int tempResult = compute(arg1F->getSparseNodeDownPtr(a, i), b);
      resF->setDownPtrWoUnlink(result, j, tempResult);
      resF->unlinkNode(tempResult);
    }
    DCASSERT(j == arg1F->getSparseNodeIndex(a, aSize - 1) + 1);
    for ( ; j < resultSize; ++j)
    {
      resF->setDownPtrWoUnlink(result, j, zeroB);
    }
    resF->unlinkNode(zeroB);
  }
}


void
MEDDLY::generic_binary_mdd::
expandB(int a, int b, int result, int resultSize)
{
  // fill in result node and return
  // result[i] = compute(a, b[i])

  if (arg2F->isFullNode(b)) {
    const int bSize = arg2F->getFullNodeSize(b);
    DCASSERT(bSize <= resultSize);
    int aZero = compute(a, 0);
    for (int i = 0; i < bSize; ++i)
    {
      int iB = arg2F->getFullNodeDownPtr(b, i);
      if (iB == 0) {
        resF->setDownPtrWoUnlink(result, i, aZero);
      }
      else {
        int tempResult = compute(a, iB);
        resF->setDownPtrWoUnlink(result, i, tempResult);
        resF->unlinkNode(tempResult);
      }
    }
    for (int i = bSize; i < resultSize; ++i)
    {
      resF->setDownPtrWoUnlink(result, i, aZero);
    }
    resF->unlinkNode(aZero);
  } else {
    DCASSERT(arg2F->isSparseNode(b));
    const int bSize = arg2F->getSparseNodeSize(b);
    // j goes through every index (like a full-node index pointer)
    int j = 0;
    int aZero = compute(a, 0);
    for (int i = 0; i < bSize; ++i, ++j)
    {
      // sparse-nodes skip indices which represent downpointer 0
      for (int index = arg2F->getSparseNodeIndex(b, i); j < index; ++j)
      {
        resF->setDownPtrWoUnlink(result, j, aZero);
      }
      // done with skipped indices; deal with the next sparse node index
      DCASSERT(j == arg2F->getSparseNodeIndex(b, i));
      int tempResult = compute(a, arg2F->getSparseNodeDownPtr(b, i));
      resF->setDownPtrWoUnlink(result, j, tempResult);
      resF->unlinkNode(tempResult);
    }
    DCASSERT(j == arg2F->getSparseNodeIndex(b, bSize - 1) + 1);
    for ( ; j < resultSize; ++j)
    {
      resF->setDownPtrWoUnlink(result, j, aZero);
    }
    resF->unlinkNode(aZero);
  }
}

void
MEDDLY::generic_binary_mdd::
fullFull (int a, int b, int result, int resultSize)
{
  DCASSERT(arg1F->isFullNode(a));
  DCASSERT(arg2F->isFullNode(b));
  DCASSERT(!resF->isReducedNode(result));

  int aSize = arg1F->getFullNodeSize(a);
  int bSize = arg2F->getFullNodeSize(b);
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;
  for ( ; i < minSize; ++i)
  {
    int tempResult = compute(arg1F->getFullNodeDownPtr(a, i),
        arg2F->getFullNodeDownPtr(b, i));
    resF->setDownPtrWoUnlink(result, i, tempResult);
    resF->unlinkNode(tempResult);
  }
  for ( ; i < aSize; ++i)
  {
    int tempResult = compute(arg1F->getFullNodeDownPtr(a, i), 0);
    resF->setDownPtrWoUnlink(result, i, tempResult);
    resF->unlinkNode(tempResult);
  }
  for ( ; i < bSize; ++i)
  {
    int tempResult = compute(0, arg2F->getFullNodeDownPtr(b, i));
    resF->setDownPtrWoUnlink(result, i, tempResult);
    resF->unlinkNode(tempResult);
  }
  if (i < resultSize) {
    int zeroZero = compute(0, 0);
    for ( ; i < resultSize; ++i)
    {
      resF->setDownPtrWoUnlink(result, i, zeroZero);
    }
    resF->unlinkNode(zeroZero);
  }
}


void
MEDDLY::generic_binary_mdd::
sparseSparse (int a, int b, int result, int resultSize)
{
  DCASSERT(arg1F->isSparseNode(a));
  DCASSERT(arg2F->isSparseNode(b));
  DCASSERT(!resF->isReducedNode(result));

  int aNnz = arg1F->getSparseNodeSize(a);
  int aSize = arg1F->getSparseNodeIndex(a, aNnz - 1) + 1;
  int bNnz = arg2F->getSparseNodeSize(b);
  int bSize = arg2F->getSparseNodeIndex(b, bNnz - 1) + 1;
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;  // index for a
  int j = 0;  // index for b
  int aIndex = arg1F->getSparseNodeIndex(a, i);
  int bIndex = arg2F->getSparseNodeIndex(b, j);
  int zeroZero = compute(0, 0);
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
      kA = arg1F->getSparseNodeDownPtr(a, i);
      ++i;
      if (i < aNnz) { aIndex = arg1F->getSparseNodeIndex(a, i); }
    }
    if (k < bIndex) {
      kB = 0;
    }
    else {
      kB = arg2F->getSparseNodeDownPtr(b, j);
      ++j;
      if (j < bNnz) { bIndex = arg2F->getSparseNodeIndex(b, j); }
    }
    if (kA == 0 && kB == 0) {
      resF->setDownPtrWoUnlink(result, k, zeroZero);
    }
    else {
      int tempResult = compute(kA, kB);
      resF->setDownPtrWoUnlink(result, k, tempResult);
      resF->unlinkNode(tempResult);
    }
  }
  for ( ; k < aSize; ++k, ++i)
  {
    aIndex = arg1F->getSparseNodeIndex(a, i);
    for ( ; k < aIndex; ++k)
    {
      // a has zeroes in these slots
      resF->setDownPtrWoUnlink(result, k, zeroZero);
    }
    DCASSERT(k == arg1F->getSparseNodeIndex(a, i));
    int tempResult = compute(arg1F->getSparseNodeDownPtr(a, i), 0);
    resF->setDownPtrWoUnlink(result, k, tempResult);
    resF->unlinkNode(tempResult);
  }
  for ( ; k < bSize; ++k, ++j)
  {
    bIndex = arg2F->getSparseNodeIndex(b, j);
    for ( ; k < bIndex; ++k)
    {
      // b has zeroes in these slots
      resF->setDownPtrWoUnlink(result, k, zeroZero);
    }
    DCASSERT(k == arg2F->getSparseNodeIndex(b, j));
    int tempResult = compute(0, arg2F->getSparseNodeDownPtr(b, j));
    resF->setDownPtrWoUnlink(result, k, tempResult);
    resF->unlinkNode(tempResult);
  }
  for ( ; k < resultSize; ++k)
  {
    resF->setDownPtrWoUnlink(result, k, zeroZero);
  }
  resF->unlinkNode(zeroZero);
}


void
MEDDLY::generic_binary_mdd::
fullSparse (int a, int b, int result, int resultSize)
{
  DCASSERT(arg1F->isFullNode(a));
  DCASSERT(arg2F->isSparseNode(b));
  DCASSERT(!resF->isReducedNode(result));

  int aSize = arg1F->getFullNodeSize(a);
  int bNnz = arg2F->getSparseNodeSize(b);
  int bSize = arg2F->getSparseNodeIndex(b, bNnz - 1) + 1;
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;
  int j = 0; // j points to sparse node index
  int bIndex = arg2F->getSparseNodeIndex(b, j);
  for ( ; i < minSize; ++i)
  {
    if (i < bIndex) {
      // b has zeroes in these slots
      int tempResult = compute(arg1F->getFullNodeDownPtr(a, i), 0);
      resF->setDownPtrWoUnlink(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
    else {
      DCASSERT(i == arg2F->getSparseNodeIndex(b, j));
      int tempResult = compute(arg1F->getFullNodeDownPtr(a, i),
          arg2F->getSparseNodeDownPtr(b, j));
      resF->setDownPtrWoUnlink(result, i, tempResult);
      resF->unlinkNode(tempResult);
      ++j;
      if (j < bNnz) { bIndex = arg2F->getSparseNodeIndex(b, j); }
    }
  }

  int zeroZero = compute(0, 0);
  for ( ; i < aSize; ++i)
  {
    int iA = arg1F->getFullNodeDownPtr(a, i);
    if (iA == 0) {
      resF->setDownPtrWoUnlink(result, i, zeroZero);
    }
    else {
      int tempResult = compute(iA, 0);
      resF->setDownPtrWoUnlink(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
  }
  for ( ; i < bSize; ++i, ++j)
  {
    for (int index = arg2F->getSparseNodeIndex(b, j); i < index; ++i)
    {
      // b has zeroes in these slots
      resF->setDownPtrWoUnlink(result, i, zeroZero);
    }
    DCASSERT(i == arg2F->getSparseNodeIndex(b, j));
    int tempResult = compute(0, arg2F->getSparseNodeDownPtr(b, j));
    resF->setDownPtrWoUnlink(result, i, tempResult);
    resF->unlinkNode(tempResult);
  }
  for ( ; i < resultSize; ++i)
  {
    resF->setDownPtrWoUnlink(result, i, zeroZero);
  }
  resF->unlinkNode(zeroZero);
}


void
MEDDLY::generic_binary_mdd::
sparseFull (int a, int b, int result, int resultSize)
{
  DCASSERT(arg1F->isSparseNode(a));
  DCASSERT(arg2F->isFullNode(b));
  DCASSERT(!resF->isReducedNode(result));

  int aNnz = arg1F->getSparseNodeSize(a);
  int aSize = arg1F->getSparseNodeIndex(a, aNnz - 1) + 1;
  int bSize = arg2F->getFullNodeSize(b);
  int minSize = aSize < bSize? aSize: bSize;

  int i = 0;
  int j = 0; // j points to sparse node index
  int aIndex = arg1F->getSparseNodeIndex(a, j);
  for ( ; i < minSize; ++i)
  {
    if (i < aIndex) {
      // a has zeroes in these slots
      int tempResult = compute(0, arg2F->getFullNodeDownPtr(b, i));
      resF->setDownPtrWoUnlink(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
    else {
      DCASSERT(i == arg1F->getSparseNodeIndex(a, j));
      int tempResult = compute(arg1F->getSparseNodeDownPtr(a, j),
          arg2F->getFullNodeDownPtr(b, i));
      resF->setDownPtrWoUnlink(result, i, tempResult);
      resF->unlinkNode(tempResult);
      ++j;
      if (j < aNnz) { aIndex = arg1F->getSparseNodeIndex(a, j); }
    }
  }

  int zeroZero = compute(0, 0);
  for ( ; i < aSize; ++i, ++j)
  {
    for (int index = arg1F->getSparseNodeIndex(a, j); i < index; ++i)
    {
      // b has zeroes in these slots
      resF->setDownPtrWoUnlink(result, i, zeroZero);
    }
    DCASSERT(i < aSize && i == arg1F->getSparseNodeIndex(a, j));
    int tempResult = compute(arg1F->getSparseNodeDownPtr(a, j), 0);
    resF->setDownPtrWoUnlink(result, i, tempResult);
    resF->unlinkNode(tempResult);
  }
  for ( ; i < bSize; ++i)
  {
    int iB = arg2F->getFullNodeDownPtr(b, i);
    if (iB == 0) {
      resF->setDownPtrWoUnlink(result, i, zeroZero);
    }
    else {
      int tempResult = compute(0, iB);
      resF->setDownPtrWoUnlink(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
  }
  for ( ; i < resultSize; ++i)
  {
    resF->setDownPtrWoUnlink(result, i, zeroZero);
  }
  resF->unlinkNode(zeroZero);
}


// ******************************************************************
// *                                                                *
// *                   generic_binary_mxd methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_mxd::generic_binary_mxd(const binary_opname* code,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mdd(code, arg1, arg2, res)
{
}

MEDDLY::generic_binary_mxd::~generic_binary_mxd()
{
}

int MEDDLY::generic_binary_mxd::compute(int a, int b) {
  return arg1F->getReductionRule() == forest::IDENTITY_REDUCED
      ? computeIdent(a, b)
      : computeNonIdent(a, b);
}

int MEDDLY::generic_binary_mxd::computeIdent(int a, int b)
{
  // Note:
  // The main difference between mxd::compute and mdd::compute is in the
  // manner in which skipped levels are dealt with.
  //
  // If aLevel < bLevel (or vice versa), MXD reduction rules specify that
  // a and b are both unprimed nodes. Therefore, when expanding b (since
  // it is at a higher level), expand the primed nodes also).

  int result = 0;
  if (checkTerminals(a, b, result))
    return result;
  if (findResult(a, b, result))
    return result;

  // expand nodes
  // 0. initialize result
  // 1. deal with special cases
  // 2. copy node a to node result
  // 3. do operation between contents of result and node b

  // 0. initialize result
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultLevel = ABS(aLevel) == ABS(bLevel) ?
    MAX(aLevel, bLevel) :
    (ABS(aLevel) > ABS(bLevel)? aLevel: bLevel);
  int resultSize = resF->getLevelSize(resultLevel);
  result = resF->createTempNode(resultLevel, resultSize);

  // 1. deal with special cases
  if (aLevel != resultLevel) {
    // WARNING: this is not the same as bLevel == resultLevel as that will
    // include the case where aLevel == bLevel == resultLevel
    // a is lower; therefore expand B
    if (resultLevel < 0)
      singleExpandB(a, b, result, resultLevel, resultSize);
    else
      expandB(a, b, result, resultLevel, resultSize);
  }
  else if (bLevel != resultLevel) {
    // WARNING: this is not the same as aLevel == resultLevel as that will
    // include the case where aLevel == bLevel == resultLevel
    // b is lower; therefore expand A
    if (resultLevel < 0)
      singleExpandA(a, b, result, resultLevel, resultSize);
    else
      expandA(a, b, result, resultLevel, resultSize);
  }
  else {
    // 2. copy node a to node result
    if (arg1F->isFullNode(a)) {
      // a is a full-node
      const int aSize = arg1F->getFullNodeSize(a);
      DCASSERT(aSize <= resultSize);
      for (int i = 0; i < aSize; ++i)
        resF->setDownPtr(result, i, arg1F->getFullNodeDownPtr(a, i));
    }
    else {
      // a is a sparse-node
      const int aSize = arg1F->getSparseNodeSize(a);
      for (int i = 0; i < aSize; ++i)
        resF->setDownPtr(result,
            arg1F->getSparseNodeIndex(a, i),
            arg1F->getSparseNodeDownPtr(a, i));
    }

    // 3. do operation between contents of result and node b
    if (arg2F->isFullNode(b)) {
      // b is a full-node
      const int bSize = arg2F->getFullNodeSize(b);
      DCASSERT(bSize <= resultSize);
      for (int i = 0; i < bSize; ++i)
      {
        int tempResult = computeIdent(
          resF->getFullNodeDownPtr(result, i),
          arg2F->getFullNodeDownPtr(b, i)
        );
        resF->setDownPtr(result, i, tempResult);
        resF->unlinkNode(tempResult);
      }
      for (int i = bSize; i < resultSize; ++i)
      {
        int tempResult =
          computeIdent(resF->getFullNodeDownPtr(result, i), 0);
        resF->setDownPtr(result, i, tempResult);
        resF->unlinkNode(tempResult);
      }
    }
    else {
      // b is a sparse-node
      const int bSize = arg2F->getSparseNodeSize(b);
      DCASSERT(arg2F->getSparseNodeIndex(b, bSize - 1) <= resultSize);
      // j goes through every index (like a full-node index pointer)
      int j = 0;
      for (int i = 0; i < bSize; ++i, ++j)
      {
        // sparse-nodes skip indices which represent downpointer 0
        // call compute of those skipped indices
        for (int index = arg2F->getSparseNodeIndex(b, i);
            j < index; ++j)
        {
          int tempResult = computeIdent(
              resF->getFullNodeDownPtr(result, j), 0);
          resF->setDownPtr(result, j, tempResult);
          resF->unlinkNode(tempResult);
        }
        // done with skipped indices; deal with the next sparse node index
        DCASSERT(j == arg2F->getSparseNodeIndex(b, i));
        int tempResult = computeIdent(
            resF->getFullNodeDownPtr(result, j),
            arg2F->getSparseNodeDownPtr(b, i));
        resF->setDownPtr(result, j, tempResult);
        resF->unlinkNode(tempResult);
      }
      DCASSERT(j == arg2F->getSparseNodeIndex(b, bSize - 1) + 1);
      for ( ; j < resultSize; ++j)
      {
        int tempResult = computeIdent(
            resF->getFullNodeDownPtr(result, j), 0);
        resF->setDownPtr(result, j, tempResult);
        resF->unlinkNode(tempResult);
      }
    }
  }

  result = resF->reduceNode(result);
  saveResult(a, b, result);
#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d) = %d\n", getName(), a, b, result);
#endif
  return result;
}



int MEDDLY::generic_binary_mxd::computeNonIdent(int a, int b)
{
  DCASSERT(resF->getReductionRule() != forest::IDENTITY_REDUCED);

  int result = 0;
  if (checkTerminals(a, b, result))
    return result;
  if (findResult(a, b, result))
    return result;

  // expand nodes
  // 0. initialize result
  // 1. deal with special cases
  // 2. copy node a to node result
  // 3. do operation between contents of result and node b

  // 0. initialize result
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultLevel = ABS(aLevel) == ABS(bLevel) ?
    MAX(aLevel, bLevel) :
    (ABS(aLevel) > ABS(bLevel)? aLevel: bLevel);
  int resultSize = resF->getLevelSize(resultLevel);
  // result = resF->createTempNode(resultLevel, resultSize);

  // get downpointers for a
  std::vector<int> A;
  if (aLevel != resultLevel) {
    if (int(A.size()) < resultSize) A.resize(resultSize);
    std::fill_n(A.begin(), resultSize, a);
  } else {
    arg1F->getDownPtrs(a, A);
  }

  // get downpointers for b
  std::vector<int> B;
  if (bLevel != resultLevel) {
    if (int(B.size()) < resultSize) B.resize(resultSize);
    std::fill_n(B.begin(), resultSize, b);
  } else {
    arg2F->getDownPtrs(b, B);
  }

  // compute C
  std::vector<int> C(resultSize, 0);
  int min = MIN(A.size(), B.size());
  std::vector<int>::iterator aIter = A.begin();
  std::vector<int>::iterator bIter = B.begin();
  std::vector<int>::iterator cIter = C.begin();

  for (std::vector<int>::iterator end = cIter + min; cIter != end; )
  {
    *cIter++ = computeNonIdent(*aIter++, *bIter++);
  }

  for (std::vector<int>::iterator end = A.end(); aIter != end; )
  {
    DCASSERT(cIter != C.end());
    *cIter++ = computeNonIdent(*aIter++, 0);
  }

  for (std::vector<int>::iterator end = B.end(); bIter != end; )
  {
    DCASSERT(cIter != C.end());
    *cIter++ = computeNonIdent(0, *bIter++);
  }

  DCASSERT(aIter == A.end() && bIter == B.end());

  if (cIter != C.end()) {
    int zeroZero = computeNonIdent(0, 0);
    if (zeroZero != 0) std::fill(cIter, C.end(), zeroZero);
    if (!resF->isTerminalNode(zeroZero)) {
      unsigned count = C.size() - (cIter - C.begin());
      resF->getInCount(zeroZero) += int(count);
    }
    resF->unlinkNode(zeroZero);
  }

  result = resF->createTempNode(resultLevel, C);
  result = resF->reduceNode(result);
  saveResult(a, b, result);
  return result;
}


void MEDDLY::generic_binary_mxd::expandA(int a, int b, int result, 
  int resultLevel, int resultSize)
{
  // result[i][j] = a[i][j] op b[i][j]
  // but b[i][j] = b only when i == j, and 0 otherwise
  // result[i][i] = a[i][i] op b, and
  // result[i][j] = a[i][j] op 0, when i != j

  // compute each result[i] using above before inserting into result node

  // when a[i] == 0, a[i][j] = 0 for all j, and
  // result[i][i] = 0 op b
  // result[i][j] = 0 op 0 for i != j
  int zeroOpZero = computeIdent(0, 0);
  int zeroOpB = computeIdent(0, b);

  int pResultSize = resF->getLevelSize(-resultLevel);

  if (arg1F->isFullNode(a)) {
    const int aSize = arg1F->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    for (int i = 0; i < aSize; ++i)
    {
      int iA = arg1F->getFullNodeDownPtr(a, i);
      int iResult = resF->createTempNode(-resultLevel, pResultSize);
      if (iA == 0) {
        for (int j = 0; j < pResultSize; ++j)
          resF->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
      }
      else if (arg1F->isFullNode(iA)) {
        const int iASize = arg1F->getFullNodeSize(iA);
        for (int j = 0; j < iASize; ++j)
        {
          int ijA = arg1F->getFullNodeDownPtr(iA, j);
          if (ijA == 0) {
            resF->setDownPtr(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(ijA, b): computeIdent(ijA, 0);
            resF->setDownPtr(iResult, j, tempResult);
            resF->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < pResultSize; ++j)
        {
          resF->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        }
      }
      else {
        DCASSERT(arg1F->isSparseNode(iA));
        const int iASize = arg1F->getSparseNodeSize(iA);
        int j = 0;
        for (int k = 0; k < iASize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = arg1F->getSparseNodeIndex(iA, k);
              j < index; ++j)
          {
            resF->setDownPtr(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == arg1F->getSparseNodeIndex(iA, k));
          int ijA = arg1F->getSparseNodeDownPtr(iA, k);
          int tempResult = (i == j)?
            computeIdent(ijA, b): computeIdent(ijA, 0);
          resF->setDownPtr(iResult, j, tempResult);
          resF->unlinkNode(tempResult);
        }
        DCASSERT(j == arg1F->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          resF->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        }
      }
      iResult = resF->reduceNode(iResult);
      resF->setDownPtr(result, i, iResult);
      resF->unlinkNode(iResult);
    }
    // TODO: can optimize when zeroOpZero == zeroOpB == 0
    for (int i = aSize; i < resultSize; ++i)
    {
      int iResult = resF->createTempNode(-resultLevel, pResultSize);
      for (int j = 0; j < pResultSize; ++j)
        resF->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
      iResult = resF->reduceNode(iResult);
      resF->setDownPtr(result, i, iResult);
      resF->unlinkNode(iResult);
    }
  }
  else {
    DCASSERT(arg1F->isSparseNode(a));
    const int aSize = arg1F->getSparseNodeSize(a);
    DCASSERT(arg1F->getSparseNodeIndex(a, aSize-1) < resultSize);
    int i = 0;
    for (int aIndex = 0; aIndex < aSize; ++aIndex, ++i)
    {
      for (int index = arg1F->getSparseNodeIndex(a, aIndex);
          i < index; ++i)
      {
        int iResult = resF->createTempNode(-resultLevel, pResultSize);
        for (int j = 0; j < pResultSize; ++j)
          resF->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        iResult = resF->reduceNode(iResult);
        resF->setDownPtr(result, i, iResult);
        resF->unlinkNode(iResult);
      }

      DCASSERT(i == arg1F->getSparseNodeIndex(a, aIndex));
      int iA = arg1F->getSparseNodeDownPtr(a, aIndex);
      int iResult = resF->createTempNode(-resultLevel, pResultSize);
      DCASSERT(iA != 0 && iA != -1);
      if (arg1F->isFullNode(iA)) {
        const int iASize = arg1F->getFullNodeSize(iA);
        for (int j = 0; j < iASize; ++j)
        {
          int ijA = arg1F->getFullNodeDownPtr(iA, j);
          if (ijA == 0) {
            resF->setDownPtr(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(ijA, b): computeIdent(ijA, 0);
            resF->setDownPtr(iResult, j, tempResult);
            resF->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < pResultSize; ++j)
        {
          resF->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        }
      }
      else {
        DCASSERT(arg1F->isSparseNode(iA));
        const int iASize = arg1F->getSparseNodeSize(iA);
        int j = 0;
        for (int k = 0; k < iASize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = arg1F->getSparseNodeIndex(iA, k);
              j < index; ++j)
          {
            resF->setDownPtr(iResult, j, (i == j)
                ? zeroOpB
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == arg1F->getSparseNodeIndex(iA, k));
          int ijA = arg1F->getSparseNodeDownPtr(iA, k);
          int tempResult = (i == j)?
            computeIdent(ijA, b): computeIdent(ijA, 0);
          resF->setDownPtr(iResult, j, tempResult);
          resF->unlinkNode(tempResult);
        }
        DCASSERT(j == arg1F->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          resF->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
        }
      }
      iResult = resF->reduceNode(iResult);
      resF->setDownPtr(result, i, iResult);
      resF->unlinkNode(iResult);
    }

    // TODO: can optimize when zeroOpZero == zeroOpB == 0
    DCASSERT(i == arg1F->getSparseNodeIndex(a, aSize - 1) + 1);
    for ( ; i < resultSize; ++i)
    {
      int iResult = resF->createTempNode(-resultLevel, pResultSize);
      for (int j = 0; j < pResultSize; ++j)
        resF->setDownPtr(iResult, j, (i == j)? zeroOpB: zeroOpZero);
      iResult = resF->reduceNode(iResult);
      resF->setDownPtr(result, i, iResult);
      resF->unlinkNode(iResult);
    }
  }

  resF->unlinkNode(zeroOpZero);
  resF->unlinkNode(zeroOpB);
}


void MEDDLY::generic_binary_mxd::singleExpandA(int a, int b, int result, 
  int resultLevel, int resultSize)
{
  // result[i] = a[i] op b[i], but b[i] = b,
  // therefore, result[i] = a[i] op b

  // compute each result[i] using above before inserting into result node

  int zeroB = computeIdent(0, b);

  if (arg1F->isFullNode(a)) {
    const int aSize = arg1F->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    int i = 0;
    for ( ; i < aSize; ++i)
    {
      int tempResult =
        computeIdent(arg1F->getFullNodeDownPtr(a, i), b);
      resF->setDownPtr(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      resF->setDownPtr(result, i, zeroB);
    }
  }
  else {
    DCASSERT(arg1F->isSparseNode(a));
    const int aSize = arg1F->getSparseNodeSize(a);
    int i = 0;
    for (int aIndex = 0; aIndex < aSize; ++aIndex, ++i)
    {
      for (int index = arg1F->getSparseNodeIndex(a, aIndex);
          i < index; ++i)
      {
        DCASSERT(i < resultSize);
        resF->setDownPtr(result, i, zeroB);
      }
      DCASSERT(i == arg1F->getSparseNodeIndex(a, aIndex));
      DCASSERT(i < resultSize);
      int tempResult =
        computeIdent(arg1F->getSparseNodeDownPtr(a, aIndex), b);
      resF->setDownPtr(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      resF->setDownPtr(result, i, zeroB);
    }
  }

  resF->unlinkNode(zeroB);
}


void MEDDLY::generic_binary_mxd::expandB(int a, int b, int result, 
  int resultLevel, int resultSize)
{
  // result[i][j] = a[i][j] op b[i][j]
  // but a[i][j] = a only when i == j, and 0 otherwise
  // result[i][i] = a op b[i][i], and
  // result[i][j] = 0 op b[i][j], when i != j

  // compute each result[i] using above before inserting into result node

  // when b[i] == 0, b[i][j] = 0 for all j, and
  // result[i][i] = a op 0
  // result[i][j] = 0 op 0 for i != j
  int zeroOpZero = computeIdent(0, 0);
  int aOpZero = computeIdent(a, 0);

  int pResultSize = resF->getLevelSize(-resultLevel);

  if (arg2F->isFullNode(b)) {
    const int bSize = arg2F->getFullNodeSize(b);
    DCASSERT(bSize <= resultSize);
    for (int i = 0; i < bSize; ++i)
    {
      int iB = arg2F->getFullNodeDownPtr(b, i);
      int iResult = resF->createTempNode(-resultLevel, pResultSize);
      if (iB == 0) {
        for (int j = 0; j < pResultSize; ++j)
          resF->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
      }
      else if (arg2F->isFullNode(iB)) {
        const int iBSize = arg2F->getFullNodeSize(iB);
        for (int j = 0; j < iBSize; ++j)
        {
          int ijB = arg2F->getFullNodeDownPtr(iB, j);
          if (ijB == 0) {
            resF->setDownPtr(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(a, ijB): computeIdent(0, ijB);
            resF->setDownPtr(iResult, j, tempResult);
            resF->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < pResultSize; ++j)
        {
          resF->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        }
      }
      else {
        DCASSERT(arg2F->isSparseNode(iB));
        const int iBSize = arg2F->getSparseNodeSize(iB);
        int j = 0;
        for (int k = 0; k < iBSize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = arg2F->getSparseNodeIndex(iB, k);
              j < index; ++j)
          {
            resF->setDownPtr(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == arg2F->getSparseNodeIndex(iB, k));
          int ijB = arg2F->getSparseNodeDownPtr(iB, k);
          int tempResult = (i == j)?
            computeIdent(a, ijB): computeIdent(0, ijB);
          resF->setDownPtr(iResult, j, tempResult);
          resF->unlinkNode(tempResult);
        }
        DCASSERT(j == arg2F->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          resF->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        }
      }
      iResult = resF->reduceNode(iResult);
      resF->setDownPtr(result, i, iResult);
      resF->unlinkNode(iResult);
    }
    // TODO: can optimize when zeroOpZero == aOpZero == 0
    for (int i = bSize; i < resultSize; ++i)
    {
      int iResult = resF->createTempNode(-resultLevel, pResultSize);
      for (int j = 0; j < pResultSize; ++j)
        resF->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
      iResult = resF->reduceNode(iResult);
      resF->setDownPtr(result, i, iResult);
      resF->unlinkNode(iResult);
    }
  }
  else {
    DCASSERT(arg2F->isSparseNode(b));
    const int bSize = arg2F->getSparseNodeSize(b);
    int i = 0;
    for (int bIndex = 0; bIndex < bSize; ++bIndex, ++i)
    {
      for (int index = arg2F->getSparseNodeIndex(b, bIndex);
          i < index; ++i)
      {
        int iResult = resF->createTempNode(-resultLevel, pResultSize);
        for (int j = 0; j < pResultSize; ++j)
          resF->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        iResult = resF->reduceNode(iResult);
        resF->setDownPtr(result, i, iResult);
        resF->unlinkNode(iResult);
      }

      DCASSERT(i == arg2F->getSparseNodeIndex(b, bIndex));
      int iB = arg2F->getSparseNodeDownPtr(b, bIndex);
      int iResult = resF->createTempNode(-resultLevel, pResultSize);
      DCASSERT(iB != 0 && iB != -1);
      if (arg2F->isFullNode(iB)) {
        const int iBSize = arg2F->getFullNodeSize(iB);
        for (int j = 0; j < iBSize; ++j)
        {
          int ijB = arg2F->getFullNodeDownPtr(iB, j);
          if (ijB == 0) {
            resF->setDownPtr(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(a, ijB): computeIdent(0, ijB);
            resF->setDownPtr(iResult, j, tempResult);
            resF->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < pResultSize; ++j)
        {
          resF->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        }
      }
      else {
        DCASSERT(arg2F->isSparseNode(iB));
        const int iBSize = arg2F->getSparseNodeSize(iB);
        int j = 0;
        for (int k = 0; k < iBSize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = arg2F->getSparseNodeIndex(iB, k);
              j < index; ++j)
          {
            resF->setDownPtr(iResult, j, (i == j)
                ? aOpZero
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == arg2F->getSparseNodeIndex(iB, k));
          int ijB = arg2F->getSparseNodeDownPtr(iB, k);
          int tempResult = (i == j)?
            computeIdent(a, ijB): computeIdent(0, ijB);
          resF->setDownPtr(iResult, j, tempResult);
          resF->unlinkNode(tempResult);
        }
        DCASSERT(j == arg2F->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          resF->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
        }
      }
      iResult = resF->reduceNode(iResult);
      DCASSERT(i < resultSize);
      resF->setDownPtr(result, i, iResult);
      resF->unlinkNode(iResult);
    }

    // TODO: can optimize when zeroOpZero == aOpZero == 0
    DCASSERT(i == arg2F->getSparseNodeIndex(b, bSize - 1) + 1);
    for ( ; i < resultSize; ++i)
    {
      int iResult = resF->createTempNode(-resultLevel, pResultSize);
      for (int j = 0; j < pResultSize; ++j)
        resF->setDownPtr(iResult, j, (i == j)? aOpZero: zeroOpZero);
      iResult = resF->reduceNode(iResult);
      resF->setDownPtr(result, i, iResult);
      resF->unlinkNode(iResult);
    }
  }

  resF->unlinkNode(zeroOpZero);
  resF->unlinkNode(aOpZero);
}


void MEDDLY::generic_binary_mxd::singleExpandB(int a, int b, int result, 
  int resultLevel, int resultSize)
{
  // result[i] = a[i] op b[i], but a[i] = b,
  // therefore, result[i] = a op b[i]

  // compute each result[i] using above before inserting into result node

  int aZero = computeIdent(a, 0);

  if (arg2F->isFullNode(b)) {
    const int bSize = arg2F->getFullNodeSize(b);
    DCASSERT(bSize <= resultSize);
    int i = 0;
    for ( ; i < bSize; ++i)
    {
      int tempResult =
        computeIdent(a, arg2F->getFullNodeDownPtr(b, i));
      resF->setDownPtr(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      resF->setDownPtr(result, i, aZero);
    }
  }
  else {
    DCASSERT(arg2F->isSparseNode(b));
    const int bSize = arg2F->getSparseNodeSize(b);
    int i = 0;
    for (int bIndex = 0; bIndex < bSize; ++bIndex, ++i)
    {
      for (int index = arg2F->getSparseNodeIndex(b, bIndex);
          i < index; ++i)
      {
        DCASSERT(i < resultSize);
        resF->setDownPtr(result, i, aZero);
      }
      DCASSERT(i == arg2F->getSparseNodeIndex(b, bIndex));
      DCASSERT(i < resultSize);
      int tempResult =
        computeIdent(a, arg2F->getSparseNodeDownPtr(b, bIndex));
      resF->setDownPtr(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      resF->setDownPtr(result, i, aZero);
    }
  }

  resF->unlinkNode(aZero);
}


// ******************************************************************
// *                                                                *
// *                 generic_binbylevel_mxd methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binbylevel_mxd
::generic_binbylevel_mxd(const binary_opname* code, expert_forest* arg1, 
  expert_forest* arg2, expert_forest* res)
 : binary_operation(code, 3, 1, arg1, arg2, res)
{
  can_commute = false;
}

MEDDLY::generic_binbylevel_mxd::~generic_binbylevel_mxd()
{
}

bool MEDDLY::generic_binbylevel_mxd::isStaleEntry(const int* data)
{
  return arg1F->isStale(data[1]) ||
         arg2F->isStale(data[2]) ||
         resF->isStale(data[3]);
}

void MEDDLY::generic_binbylevel_mxd::discardEntry(const int* data)
{
  arg1F->uncacheNode(data[1]);
  arg2F->uncacheNode(data[2]);
  resF->uncacheNode(data[3]);
}

void
MEDDLY::generic_binbylevel_mxd::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d, %d, %d): %d]", getName(), 
    data[0], data[1], data[2], data[3]
  );
}

void MEDDLY::generic_binbylevel_mxd
::compute(const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  int result = compute(
    resF->getDomain()->getNumVariables(), a.getNode(), b.getNode()
  );
  c.set(result, 0, resF->getNodeLevel(result));
}

int MEDDLY::generic_binbylevel_mxd
::compute(int resultLevel, int a, int b)
{
  return resF->getReductionRule() == forest::IDENTITY_REDUCED
      ? computeIdent(resultLevel, a, b)
      : computeNonIdent(resultLevel, a, b);
}

int MEDDLY::generic_binbylevel_mxd
::computeIdent(int resultLevel, int a, int b)
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
  DCASSERT(resF->getReductionRule() == forest::IDENTITY_REDUCED);

  if (resultLevel == 0) {
    checkTerminals(a, b, result);
    return result;
  }

  if (findResult(resultLevel, a, b, result)) {
    return result;
  }

  // expand nodes
  // 0. initialize result
  // 1. deal with special cases
  // 2. copy node a to node result
  // 3. do operation between contents of result and node b

  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  // 0. initialize result
  int resultSize = resF->getLevelSize(resultLevel);
  result = resF->createTempNode(resultLevel, resultSize, false);

  if (aLevel != resultLevel && bLevel != resultLevel) {
    expandSkippedLevel(a, b, result, resultLevel, resultSize);
  }
  else if (aLevel != resultLevel) {
    // aLevel != resultLevel, bLevel == resultLevel
    if (resultLevel < 0)
      singleExpandB(a, b, result, resultLevel, resultSize);
    else
      expandB(a, b, result, resultLevel, resultSize);
  }
  else if (bLevel != resultLevel) {
    // aLevel == resultLevel, bLevel != resultLevel
    if (resultLevel < 0)
      singleExpandA(a, b, result, resultLevel, resultSize);
    else
      expandA(a, b, result, resultLevel, resultSize);
  }
  else {
    // aLevel == resultLevel, bLevel == resultLevel
    DCASSERT(aLevel == resultLevel && bLevel == resultLevel);
    // 2. copy node a to node result
    if (arg1F->isFullNode(a)) {
      // a is a full-node
      const int aSize = arg1F->getFullNodeSize(a);
      DCASSERT(aSize <= resultSize);
      for (int i = 0; i < aSize; ++i)
        resF->setDownPtrWoUnlink(result,
            i, arg1F->getFullNodeDownPtr(a, i));
      for (int i = aSize; i < resultSize; ++i)
        resF->setDownPtrWoUnlink(result, i, 0);
    }
    else {
      // a is a sparse-node
      const int aSize = arg1F->getSparseNodeSize(a);
      for (int i = 0; i < resultSize; ++i)
        resF->setDownPtrWoUnlink(result, i, 0);
      for (int i = 0; i < aSize; ++i)
        resF->setDownPtrWoUnlink(result,
            arg1F->getSparseNodeIndex(a, i),
            arg1F->getSparseNodeDownPtr(a, i));
    }

    // 3. do operation between contents of result and node b
    int nextLevel = (resultLevel > 0)? -resultLevel: -resultLevel-1;
    if (arg2F->isFullNode(b)) {
      // b is a full-node
      const int bSize = arg2F->getFullNodeSize(b);
      DCASSERT(bSize <= resultSize);
      for (int i = 0; i < bSize; ++i)
      {
        int tempResult = computeIdent(nextLevel,
            resF->getFullNodeDownPtr(result, i),
            arg2F->getFullNodeDownPtr(b, i));
        resF->setDownPtr(result, i, tempResult);
        resF->unlinkNode(tempResult);
      }
      for (int i = bSize; i < resultSize; ++i)
      {
        int tempResult = computeIdent(nextLevel,
            resF->getFullNodeDownPtr(result, i), 0);
        resF->setDownPtr(result, i, tempResult);
        resF->unlinkNode(tempResult);
      }
    }
    else {
      // b is a sparse-node
      const int bSize = arg2F->getSparseNodeSize(b);
      DCASSERT(arg2F->getSparseNodeIndex(b, bSize - 1) <= resultSize);
      // j goes through every index (like a full-node index pointer)
      int j = 0;
      for (int i = 0; i < bSize; ++i, ++j)
      {
        // sparse-nodes skip indices which represent downpointer 0
        // call compute of those skipped indices
        for (int index = arg2F->getSparseNodeIndex(b, i);
            j < index; ++j)
        {
          int tempResult = computeIdent(nextLevel,
              resF->getFullNodeDownPtr(result, j), 0);
          resF->setDownPtr(result, j, tempResult);
          resF->unlinkNode(tempResult);
        }
        // done with skipped indices; deal with the next sparse node index
        DCASSERT(j == arg2F->getSparseNodeIndex(b, i));
        int tempResult = computeIdent(nextLevel,
            resF->getFullNodeDownPtr(result, j),
            arg2F->getSparseNodeDownPtr(b, i));
        resF->setDownPtr(result, j, tempResult);
        resF->unlinkNode(tempResult);
      }
      DCASSERT(j == arg2F->getSparseNodeIndex(b, bSize - 1) + 1);
      for ( ; j < resultSize; ++j)
      {
        int tempResult = computeIdent(nextLevel,
            resF->getFullNodeDownPtr(result, j), 0);
        resF->setDownPtr(result, j, tempResult);
        resF->unlinkNode(tempResult);
      }
    }
  }

  result = resF->reduceNode(result);
  saveResult(resultLevel, a, b, result);
#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d, %d) = %d\n", getName(), 
    resultLevel, a, b, result
  );
#endif
  return result;
}

int MEDDLY::generic_binbylevel_mxd
::computeNonIdent(int resultLevel, int a, int b)
{
  DCASSERT(resF->getReductionRule() != forest::IDENTITY_REDUCED);

  int result = 0;

  if (resultLevel == 0) {
    checkTerminals(a, b, result);
    return result;
  }

  if (findResult(resultLevel, a, b, result)) {
    return result;
  }


  // expand nodes
  // 0. initialize result
  // 1. deal with special cases
  // 2. copy node a to node result
  // 3. do operation between contents of result and node b

  // 0. initialize result
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultSize = resF->getLevelSize(resultLevel);
  int nextLevel = (resultLevel > 0)? -resultLevel: -resultLevel-1;

  if (aLevel != resultLevel && bLevel != resultLevel) {

    // c[i] = a op b
    int c = computeNonIdent(nextLevel, a, b);
    std::vector<int> C(resultSize, c);
    if (!resF->isTerminalNode(c))
      resF->getInCount(c) += resultSize;
    resF->unlinkNode(c);
    result = resF->createTempNode(resultLevel, C);
    result = resF->reduceNode(result);

  }
  else {

    // get downpointers for a
    std::vector<int> A;
    if (aLevel != resultLevel) {
      if (int(A.size()) < resultSize) A.resize(resultSize);
      std::fill_n(A.begin(), resultSize, a);
    } else {
      arg1F->getDownPtrs(a, A);
    }

    // get downpointers for b
    std::vector<int> B;
    if (bLevel != resultLevel) {
      if (int(B.size()) < resultSize) B.resize(resultSize);
      std::fill_n(B.begin(), resultSize, b);
    } else {
      arg2F->getDownPtrs(b, B);
    }

    // compute C
    std::vector<int> C(resultSize, 0);

    int min = MIN(A.size(), B.size());
    std::vector<int>::iterator aIter = A.begin();
    std::vector<int>::iterator bIter = B.begin();
    std::vector<int>::iterator cIter = C.begin();

    for (std::vector<int>::iterator end = cIter + min; cIter != end; )
    {
      *cIter++ = computeNonIdent(nextLevel, *aIter++, *bIter++);
    }

    for (std::vector<int>::iterator end = A.end(); aIter != end; )
    {
      DCASSERT(cIter != C.end());
      *cIter++ = computeNonIdent(nextLevel, *aIter++, 0);
    }

    for (std::vector<int>::iterator end = B.end(); bIter != end; )
    {
      DCASSERT(cIter != C.end());
      *cIter++ = computeNonIdent(nextLevel, 0, *bIter++);
    }

    DCASSERT(aIter == A.end() && bIter == B.end());

    if (cIter != C.end()) {
      int zeroZero = computeNonIdent(nextLevel, 0, 0);
      if (zeroZero != 0) std::fill(cIter, C.end(), zeroZero);
      if (!resF->isTerminalNode(zeroZero)) {
        unsigned count = C.size() - (cIter - C.begin());
        resF->getInCount(zeroZero) += int(count);
      }
      resF->unlinkNode(zeroZero);
    }

    result = resF->createTempNode(resultLevel, C);
    result = resF->reduceNode(result);

  }

  saveResult(resultLevel, a, b, result);
  return result;
}


void MEDDLY::generic_binbylevel_mxd::expandSkippedLevel(int a, int b, 
  int result, int resultLevel, int resultSize)
{
  DCASSERT(arg1F->getNodeLevel(a) != resultLevel);
  DCASSERT(arg2F->getNodeLevel(b) != resultLevel);
  if (resultLevel > 0) {
    DCASSERT(arg1F->getNodeLevel(a) >= 0);
    DCASSERT(arg2F->getNodeLevel(b) >= 0);
    // both a and b are below result
    int zeroZeroAtOneLevelBelow = computeIdent(resultLevel-1, 0, 0);
    int aBAtOneLevelBelow = computeIdent(resultLevel-1, a, b);
    int pResultSize = resF->getLevelSize(-resultLevel);

    for (int i = 0; i < resultSize; ++i)
    {
      // create primed node at -resultLevel
      int p = resF->createTempNode(-resultLevel, pResultSize, false);
      for (int j = 0; j < pResultSize; ++j)
      {
        // create unprimed node at resultLevel-1
        resF->setDownPtrWoUnlink(p, j,
            (i == j)? aBAtOneLevelBelow: zeroZeroAtOneLevelBelow);
      }
      p = resF->reduceNode(p);
      resF->setDownPtrWoUnlink(result, i, p);
      resF->unlinkNode(p);
    }
    resF->unlinkNode(zeroZeroAtOneLevelBelow);
    resF->unlinkNode(aBAtOneLevelBelow);
  }
  else {
    DCASSERT(resultLevel < 0);
    DCASSERT(a == 0 && b == 0);
    int zeroZeroAtOneLevelBelow = computeIdent(-resultLevel-1, 0, 0);
    for (int i = 0; i < resultSize; ++i)
    {
      resF->setDownPtrWoUnlink(result, i, zeroZeroAtOneLevelBelow);
    }
    resF->unlinkNode(zeroZeroAtOneLevelBelow);
  }
}


void MEDDLY::generic_binbylevel_mxd::expandA(int a, int b,
  int result, int resultLevel, int resultSize)
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
  int zeroOpZero = computeIdent(resultLevel-1, 0, 0);
  int zeroOpB = computeIdent(resultLevel-1, 0, b);

  int pResultSize = resF->getLevelSize(-resultLevel);

  if (arg1F->isFullNode(a)) {
    const int aSize = arg1F->getFullNodeSize(a);
    for (int i = 0; i < aSize; ++i)
    {
      int iA = arg1F->getFullNodeDownPtr(a, i);
      int iResult =
        resF->createTempNode(-resultLevel, pResultSize, false);
      if (iA == 0) {
        for (int j = 0; j < pResultSize; ++j)
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
      }
      else if (arg1F->isFullNode(iA)) {
        const int iASize = arg1F->getFullNodeSize(iA);
        for (int j = 0; j < iASize; ++j)
        {
          int ijA = arg1F->getFullNodeDownPtr(iA, j);
          if (ijA == 0) {
            resF->setDownPtrWoUnlink(iResult, j, i == j
                ? zeroOpB
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(resultLevel-1, ijA, b):
              computeIdent(resultLevel-1, ijA, 0);
            resF->setDownPtrWoUnlink(iResult, j, tempResult);
            resF->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < pResultSize; ++j)
        {
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
        }
      }
      else {
        DCASSERT(arg1F->isSparseNode(iA));
        const int iASize = arg1F->getSparseNodeSize(iA);
        int j = 0;
        for (int k = 0; k < iASize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = arg1F->getSparseNodeIndex(iA, k);
              j < index; ++j)
          {
            resF->setDownPtrWoUnlink(iResult, j, i == j
                ? zeroOpB
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == arg1F->getSparseNodeIndex(iA, k));
          int ijA = arg1F->getSparseNodeDownPtr(iA, k);
          int tempResult = (i == j)?
            computeIdent(resultLevel-1, ijA, b):
            computeIdent(resultLevel-1, ijA, 0);
          resF->setDownPtrWoUnlink(iResult, j, tempResult);
          resF->unlinkNode(tempResult);
        }
        DCASSERT(j == arg1F->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
        }
      }
      iResult = resF->reduceNode(iResult);
      resF->setDownPtrWoUnlink(result, i, iResult);
      resF->unlinkNode(iResult);
    }
    // TODO: can optimize when zeroOpZero == zeroOpB == 0
    for (int i = aSize; i < resultSize; ++i)
    {
      int iResult =
        resF->createTempNode(-resultLevel, pResultSize, false);
      for (int j = 0; j < pResultSize; ++j)
        resF->setDownPtrWoUnlink(iResult, j,
            i == j? zeroOpB: zeroOpZero);
      iResult = resF->reduceNode(iResult);
      resF->setDownPtrWoUnlink(result, i, iResult);
      resF->unlinkNode(iResult);
    }
  }
  else {
    DCASSERT(arg1F->isSparseNode(a));
    const int aSize = arg1F->getSparseNodeSize(a);
    int i = 0;
    for (int aIndex = 0; aIndex < aSize; ++aIndex, ++i)
    {
      for (int index = arg1F->getSparseNodeIndex(a, aIndex);
          i < index; ++i)
      {
        int iResult =
          resF->createTempNode(-resultLevel, pResultSize, false);
        for (int j = 0; j < pResultSize; ++j)
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
        iResult = resF->reduceNode(iResult);
        resF->setDownPtrWoUnlink(result, i, iResult);
        resF->unlinkNode(iResult);
      }

      DCASSERT(i == arg1F->getSparseNodeIndex(a, aIndex));
      int iA = arg1F->getSparseNodeDownPtr(a, aIndex);
      int iResult =
        resF->createTempNode(-resultLevel, pResultSize, false);
      DCASSERT(iA != 0 && iA != -1);
      if (arg1F->isFullNode(iA)) {
        const int iASize = arg1F->getFullNodeSize(iA);
        for (int j = 0; j < iASize; ++j)
        {
          int ijA = arg1F->getFullNodeDownPtr(iA, j);
          if (ijA == 0) {
            resF->setDownPtrWoUnlink(iResult, j,
                i == j? zeroOpB : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(resultLevel-1, ijA, b):
              computeIdent(resultLevel-1, ijA, 0);
            resF->setDownPtrWoUnlink(iResult, j, tempResult);
            resF->unlinkNode(tempResult);
          }
        }
        for (int j = iASize; j < pResultSize; ++j)
        {
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
        }
      }
      else {
        DCASSERT(arg1F->isSparseNode(iA));
        const int iASize = arg1F->getSparseNodeSize(iA);
        int j = 0;
        for (int k = 0; k < iASize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = arg1F->getSparseNodeIndex(iA, k);
              j < index; ++j)
          {
            resF->setDownPtrWoUnlink(iResult, j, i == j
                ? zeroOpB
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == arg1F->getSparseNodeIndex(iA, k));
          int ijA = arg1F->getSparseNodeDownPtr(iA, k);
          int tempResult = i == j?
            computeIdent(resultLevel-1, ijA, b):
            computeIdent(resultLevel-1, ijA, 0);
          resF->setDownPtrWoUnlink(iResult, j, tempResult);
          resF->unlinkNode(tempResult);
        }
        DCASSERT(j == arg1F->getSparseNodeIndex(iA, iASize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? zeroOpB: zeroOpZero);
        }
      }
      iResult = resF->reduceNode(iResult);
      resF->setDownPtrWoUnlink(result, i, iResult);
      resF->unlinkNode(iResult);
    }

    // TODO: can optimize when zeroOpZero == zeroOpB == 0
    DCASSERT(i == arg1F->getSparseNodeIndex(a, aSize - 1) + 1);
    for ( ; i < resultSize; ++i)
    {
      int iResult =
        resF->createTempNode(-resultLevel, pResultSize, false);
      for (int j = 0; j < pResultSize; ++j)
        resF->setDownPtrWoUnlink(iResult, j,
            i == j? zeroOpB: zeroOpZero);
      iResult = resF->reduceNode(iResult);
      resF->setDownPtrWoUnlink(result, i, iResult);
      resF->unlinkNode(iResult);
    }
  }

  resF->unlinkNode(zeroOpZero);
  resF->unlinkNode(zeroOpB);
}


void MEDDLY::generic_binbylevel_mxd::singleExpandA(int a, int b,
  int result, int resultLevel, int resultSize)
{
  // result[i] = a[i] op b[i], but b[i] = b,
  // therefore, result[i] = a[i] op b

  // compute each result[i] using above before inserting into result node

  DCASSERT(resultLevel < 0);
  int zeroB = computeIdent(-resultLevel-1, 0, b);

  if (arg1F->isFullNode(a)) {
    const int aSize = arg1F->getFullNodeSize(a);
    int i = 0;
    for ( ; i < aSize; ++i)
    {
      int tempResult = computeIdent(-resultLevel-1,
          arg1F->getFullNodeDownPtr(a, i), b);
      resF->setDownPtrWoUnlink(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      resF->setDownPtrWoUnlink(result, i, zeroB);
    }
  }
  else {
    DCASSERT(arg1F->isSparseNode(a));
    const int aSize = arg1F->getSparseNodeSize(a);
    int i = 0;
    for (int aIndex = 0; aIndex < aSize; ++aIndex, ++i)
    {
      for (int index = arg1F->getSparseNodeIndex(a, aIndex);
          i < index; ++i)
      {
        resF->setDownPtrWoUnlink(result, i, zeroB);
      }
      DCASSERT(i == arg1F->getSparseNodeIndex(a, aIndex));
      int tempResult = computeIdent(-resultLevel-1,
          arg1F->getSparseNodeDownPtr(a, aIndex), b);
      resF->setDownPtrWoUnlink(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      resF->setDownPtrWoUnlink(result, i, zeroB);
    }
  }

  resF->unlinkNode(zeroB);
}


void MEDDLY::generic_binbylevel_mxd::expandB(int a, int b,
  int result, int resultLevel, int resultSize)
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
  int zeroOpZero = computeIdent(resultLevel-1, 0, 0);
  int aOpZero = computeIdent(resultLevel-1, a, 0);

  int pResultSize = resF->getLevelSize(-resultLevel);

  if (arg2F->isFullNode(b)) {
    const int bSize = arg2F->getFullNodeSize(b);
    for (int i = 0; i < bSize; ++i)
    {
      int iB = arg2F->getFullNodeDownPtr(b, i);
      int iResult =
        resF->createTempNode(-resultLevel, pResultSize, false);
      if (iB == 0) {
        for (int j = 0; j < pResultSize; ++j)
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
      }
      else if (arg2F->isFullNode(iB)) {
        const int iBSize = arg2F->getFullNodeSize(iB);
        for (int j = 0; j < iBSize; ++j)
        {
          int ijB = arg2F->getFullNodeDownPtr(iB, j);
          if (ijB == 0) {
            resF->setDownPtrWoUnlink(iResult, j, i == j
                ? aOpZero
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(resultLevel-1, a, ijB):
              computeIdent(resultLevel-1, 0, ijB);
            resF->setDownPtrWoUnlink(iResult, j, tempResult);
            resF->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < pResultSize; ++j)
        {
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
        }
      }
      else {
        DCASSERT(arg2F->isSparseNode(iB));
        const int iBSize = arg2F->getSparseNodeSize(iB);
        int j = 0;
        for (int k = 0; k < iBSize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = arg2F->getSparseNodeIndex(iB, k);
              j < index; ++j)
          {
            resF->setDownPtrWoUnlink(iResult, j, i == j
                ? aOpZero
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == arg2F->getSparseNodeIndex(iB, k));
          int ijB = arg2F->getSparseNodeDownPtr(iB, k);
          int tempResult = (i == j)?
            computeIdent(resultLevel-1, a, ijB):
            computeIdent(resultLevel-1, 0, ijB);
          resF->setDownPtrWoUnlink(iResult, j, tempResult);
          resF->unlinkNode(tempResult);
        }
        DCASSERT(j == arg2F->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
        }
      }
      iResult = resF->reduceNode(iResult);
      resF->setDownPtrWoUnlink(result, i, iResult);
      resF->unlinkNode(iResult);
    }
    // TODO: can optimize when zeroOpZero == aOpZero == 0
    for (int i = bSize; i < resultSize; ++i)
    {
      int iResult =
        resF->createTempNode(-resultLevel, pResultSize, false);
      for (int j = 0; j < pResultSize; ++j)
        resF->setDownPtrWoUnlink(iResult, j,
            i == j? aOpZero: zeroOpZero);
      iResult = resF->reduceNode(iResult);
      resF->setDownPtrWoUnlink(result, i, iResult);
      resF->unlinkNode(iResult);
    }
  }
  else {
    DCASSERT(arg2F->isSparseNode(b));
    const int bSize = arg2F->getSparseNodeSize(b);
    int i = 0;
    for (int bIndex = 0; bIndex < bSize; ++bIndex, ++i)
    {
      for (int index = arg2F->getSparseNodeIndex(b, bIndex);
          i < index; ++i)
      {
        int iResult =
          resF->createTempNode(-resultLevel, pResultSize, false);
        for (int j = 0; j < pResultSize; ++j)
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
        iResult = resF->reduceNode(iResult);
        resF->setDownPtrWoUnlink(result, i, iResult);
        resF->unlinkNode(iResult);
      }

      DCASSERT(i == arg2F->getSparseNodeIndex(b, bIndex));
      int iB = arg2F->getSparseNodeDownPtr(b, bIndex);
      int iResult =
        resF->createTempNode(-resultLevel, pResultSize, false);
      DCASSERT(iB != 0 && iB != -1);
      if (arg2F->isFullNode(iB)) {
        const int iBSize = arg2F->getFullNodeSize(iB);
        for (int j = 0; j < iBSize; ++j)
        {
          int ijB = arg2F->getFullNodeDownPtr(iB, j);
          if (ijB == 0) {
            resF->setDownPtrWoUnlink(iResult, j, i == j
                ? aOpZero
                : zeroOpZero);
          }
          else {
            int tempResult = (i == j)?
              computeIdent(resultLevel-1, a, ijB):
              computeIdent(resultLevel-1, 0, ijB);
            resF->setDownPtrWoUnlink(iResult, j, tempResult);
            resF->unlinkNode(tempResult);
          }
        }
        for (int j = iBSize; j < pResultSize; ++j)
        {
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
        }
      }
      else {
        DCASSERT(arg2F->isSparseNode(iB));
        const int iBSize = arg2F->getSparseNodeSize(iB);
        int j = 0;
        for (int k = 0; k < iBSize; ++k, ++j)
        {
          // sparse-nodes skip indices which represent downpointer 0
          // call compute of those skipped indices
          for (int index = arg2F->getSparseNodeIndex(iB, k);
              j < index; ++j)
          {
            resF->setDownPtrWoUnlink(iResult, j, i == j
                ? aOpZero
                : zeroOpZero);
          }
          // done with skipped indices; deal with the next sparse node index
          DCASSERT(j == arg2F->getSparseNodeIndex(iB, k));
          int ijB = arg2F->getSparseNodeDownPtr(iB, k);
          int tempResult = i == j?
            computeIdent(resultLevel-1, a, ijB):
            computeIdent(resultLevel-1, 0, ijB);
          resF->setDownPtrWoUnlink(iResult, j, tempResult);
          resF->unlinkNode(tempResult);
        }
        DCASSERT(j == arg2F->getSparseNodeIndex(iB, iBSize - 1) + 1);
        for ( ; j < pResultSize; ++j)
        {
          resF->setDownPtrWoUnlink(iResult, j,
              i == j? aOpZero: zeroOpZero);
        }
      }
      iResult = resF->reduceNode(iResult);
      resF->setDownPtrWoUnlink(result, i, iResult);
      resF->unlinkNode(iResult);
    }

    // TODO: can optimize when zeroOpZero == aOpZero == 0
    DCASSERT(i == arg2F->getSparseNodeIndex(b, bSize - 1) + 1);
    for ( ; i < resultSize; ++i)
    {
      int iResult =
        resF->createTempNode(-resultLevel, pResultSize, false);
      for (int j = 0; j < pResultSize; ++j)
        resF->setDownPtrWoUnlink(iResult, j,
            i == j? aOpZero: zeroOpZero);
      iResult = resF->reduceNode(iResult);
      resF->setDownPtrWoUnlink(result, i, iResult);
      resF->unlinkNode(iResult);
    }
  }

  resF->unlinkNode(zeroOpZero);
  resF->unlinkNode(aOpZero);
}


void MEDDLY::generic_binbylevel_mxd::singleExpandB(int a, int b,
  int result, int resultLevel, int resultSize)
{
  // result[i] = a[i] op b[i], but a[i] = b,
  // therefore, result[i] = a op b[i]

  // compute each result[i] using above before inserting into result node

  DCASSERT(resultLevel < 0);
  int aZero = computeIdent(-resultLevel-1, a, 0);

  if (arg2F->isFullNode(b)) {
    const int bSize = arg2F->getFullNodeSize(b);
    int i = 0;
    for ( ; i < bSize; ++i)
    {
      int tempResult = computeIdent(-resultLevel-1,
          a, arg2F->getFullNodeDownPtr(b, i));
      resF->setDownPtrWoUnlink(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      resF->setDownPtrWoUnlink(result, i, aZero);
    }
  }
  else {
    DCASSERT(arg2F->isSparseNode(b));
    const int bSize = arg2F->getSparseNodeSize(b);
    int i = 0;
    for (int bIndex = 0; bIndex < bSize; ++bIndex, ++i)
    {
      for (int index = arg2F->getSparseNodeIndex(b, bIndex);
          i < index; ++i)
      {
        resF->setDownPtrWoUnlink(result, i, aZero);
      }
      DCASSERT(i == arg2F->getSparseNodeIndex(b, bIndex));
      int tempResult = computeIdent(-resultLevel-1,
          a, arg2F->getSparseNodeDownPtr(b, bIndex));
      resF->setDownPtrWoUnlink(result, i, tempResult);
      resF->unlinkNode(tempResult);
    }
    for ( ; i < resultSize; ++i)
    {
      resF->setDownPtrWoUnlink(result, i, aZero);
    }
  }

  resF->unlinkNode(aZero);
}





// ******************************************************************
// *                                                                *
// *                   generic_binary_ev  methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_ev::generic_binary_ev(const binary_opname* code,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : binary_operation(code, 4, 2, arg1, arg2, res)
{
  can_commute = false;
}

MEDDLY::generic_binary_ev::~generic_binary_ev()
{
}

bool MEDDLY::generic_binary_ev::isStaleEntry(const int* data)
{
  return arg1F->isStale(data[1]) ||
         arg2F->isStale(data[3]) ||
         resF->isStale(data[5]);
}

void MEDDLY::generic_binary_ev::discardEntry(const int* data)
{
  arg1F->uncacheNode(data[1]);
  arg2F->uncacheNode(data[3]);
  resF->uncacheNode(data[5]);
}

// ******************************************************************
// *                                                                *
// *                 generic_binary_evplus  methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_evplus::generic_binary_evplus(const binary_opname* code,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_ev(code, arg1, arg2, res)
{
}

MEDDLY::generic_binary_evplus::~generic_binary_evplus()
{
}

void MEDDLY::generic_binary_evplus::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(<%d:%d>, <%d:%d>): <%d:%d>]",
      getName(),
      data[0], data[1], data[2], data[3], data[4], data[5]
  );
}

void MEDDLY::generic_binary_evplus
::compute(const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  int result, ev, aev, bev;
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(aev, a.getNode(), bev, b.getNode(), ev, result);
  c.set(result, ev, resF->getNodeLevel(result));
}


void MEDDLY::generic_binary_evplus
::compute(int aev, int a, int bev, int b, int& cev, int& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;
  if (findResult(aev, a, bev, b, cev, c))
    return;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = resF->getLevelSize(resultLevel);

  // Three vectors: operands a and b, and result c

  if (aLevel < resultLevel) {
    // expand b
    // result[i] = a op b[i]
    std::vector<int> B(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<int> Bev(resultSize, INF);
    std::vector<int> Cev(resultSize, INF);

    int aZero, aZeroev;
    compute(aev, a, INF, 0, aZeroev, aZero);

    arg2F->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterBev = Bev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    for ( ; iterB != B.end(); iterB++, iterC++, iterBev++, iterCev++)
    {
      if (*iterB == 0) {
        *iterC = aZero;
        resF->linkNode(aZero);
        *iterCev = aZeroev;
      } else {
        compute(aev, a, *iterBev + bev, *iterB, *iterCev, *iterC);
      }
    }
    resF->unlinkNode(aZero);
    c = resF->createTempNode(resultLevel, C, Cev);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    std::vector<int> A(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<int> Aev(resultSize, INF);
    std::vector<int> Cev(resultSize, INF);

    int zeroB, zeroBev;
    compute(INF, 0, bev, b, zeroBev, zeroB);

    arg1F->getDownPtrsAndEdgeValues(a, A, Aev);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterAev = Aev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end(); iterA++, iterC++, iterAev++, iterCev++)
    {
      if (*iterA == 0) {
        *iterC = zeroB;
        resF->linkNode(zeroB);
        *iterCev = zeroBev;
      } else {
        compute(*iterAev + aev, *iterA, bev, b, *iterCev, *iterC);
      }
    }
    resF->unlinkNode(zeroB);
    c = resF->createTempNode(resultLevel, C, Cev);
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
    compute(INF, 0, INF, 0, zeroZeroEv, zeroZero);

    arg1F->getDownPtrsAndEdgeValues(a, A, Aev);
    arg2F->getDownPtrsAndEdgeValues(b, B, Bev);
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
          resF->linkNode(zeroZero);
          *iterCev = zeroZeroEv;
        } else {
          compute(INF, 0, *iterBev + bev, *iterB, *iterCev, *iterC);
        }
      } else if (*iterB == 0) {
        compute(*iterAev + aev, *iterA, INF, 0, *iterCev, *iterC);
      } else {
        compute(*iterAev + aev, *iterA, *iterBev + bev, *iterB, 
          *iterCev, *iterC);
      }
    }

    resF->unlinkNode(zeroZero);
    c = resF->createTempNode(resultLevel, C, Cev);
  }

  // save result in compute cache and return it

#if 0
  printf("reduce(%d): ", result);
  result = resF->reduceNode(result);
  printf("%d  [", result);
  for (unsigned i = 0; i < C.size(); i++ )
  {
    printf("%d ", C[i]);
  }
  printf("]\n");
#else
  resF->normalizeAndReduceNode(c, cev);
#endif

  saveResult(aev, a, bev, b, cev, c);
}


// ******************************************************************
// *                                                                *
// *                 generic_binary_evtimes methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_evtimes
::generic_binary_evtimes(const binary_opname* code, expert_forest* arg1, 
  expert_forest* arg2, expert_forest* res)
: generic_binary_ev(code, arg1, arg2, res)
{
}

MEDDLY::generic_binary_evtimes::~generic_binary_evtimes()
{
}

void MEDDLY::generic_binary_evtimes
::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(<%f:%d>, <%f:%d>): <%f:%d>]",
      getName(),
      toFloat(data[0]), data[1], 
      toFloat(data[2]), data[3], 
      toFloat(data[4]), data[5]
  );
}

void MEDDLY::generic_binary_evtimes
::compute(const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  int result; 
  float ev, aev, bev;
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(aev, a.getNode(), bev, b.getNode(), ev, result);
  c.set(result, ev, resF->getNodeLevel(result));
}

void MEDDLY::generic_binary_evtimes
::compute(float aev, int a, float bev, int b, float& cev, int& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;
  if (findResult(aev, a, bev, b, cev, c))
    return;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = resF->getLevelSize(resultLevel);

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
    compute(aev, a, NAN, 0, aZeroev, aZero);

    arg2F->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<float>::iterator iterBev = Bev.begin();
    std::vector<float>::iterator iterCev = Cev.begin();
    for ( ; iterB != B.end(); iterB++, iterC++, iterBev++, iterCev++)
    {
      if (*iterB == 0) {
        *iterC = aZero;
        resF->linkNode(aZero);
        *iterCev = aZeroev;
      } else {
        compute(aev, a, *iterBev + bev, *iterB, *iterCev, *iterC);
      }
    }
    resF->unlinkNode(aZero);
    c = resF->createTempNode(resultLevel, C, Cev);
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
    compute(NAN, 0, bev, b, zeroBev, zeroB);

    arg1F->getDownPtrsAndEdgeValues(a, A, Aev);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<float>::iterator iterAev = Aev.begin();
    std::vector<float>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end(); iterA++, iterC++, iterAev++, iterCev++)
    {
      if (*iterA == 0) {
        *iterC = zeroB;
        resF->linkNode(zeroB);
        *iterCev = zeroBev;
      } else {
        compute(*iterAev + aev, *iterA, bev, b, *iterCev, *iterC);
      }
    }
    resF->unlinkNode(zeroB);
    c = resF->createTempNode(resultLevel, C, Cev);
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
    compute(NAN, 0, NAN, 0, zeroZeroEv, zeroZero);

    arg1F->getDownPtrsAndEdgeValues(a, A, Aev);
    arg2F->getDownPtrsAndEdgeValues(b, B, Bev);
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
          resF->linkNode(zeroZero);
          *iterCev = zeroZeroEv;
        } else {
          compute(NAN, 0, *iterBev + bev, *iterB, *iterCev, *iterC);
        }
      } else if (*iterB == 0) {
        compute(*iterAev + aev, *iterA, NAN, 0, *iterCev, *iterC);
      } else {
        compute(*iterAev + aev, *iterA, *iterBev + bev, *iterB, 
          *iterCev, *iterC);
      }
    }

    resF->unlinkNode(zeroZero);
    c = resF->createTempNode(resultLevel, C, Cev);
  }

  // save result in compute cache and return it

#if 0
  printf("reduce(%d): ", result);
  result = resF->reduceNode(result);
  printf("%d  [", result);
  for (unsigned i = 0; i < C.size(); i++ )
  {
    printf("%d ", C[i]);
  }
  printf("]\n");
#else
  resF->normalizeAndReduceNode(c, cev);
#endif

  saveResult(aev, a, bev, b, cev, c);
}




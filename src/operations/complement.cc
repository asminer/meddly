
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
#include "complement.h"
#include "../compute_table.h"

namespace MEDDLY {
  class compl_mdd;
  class compl_mxd;
  class compl_opname;
};

// #define DEBUG_MXD_COMPL

// ******************************************************************
// *                                                                *
// *                        compl_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_mdd : public unary_operation {
  public:
    compl_mdd(const unary_opname* oc, expert_forest* arg, expert_forest* res);

    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, dd_edge& b);

    int compute(int a);
};

MEDDLY::compl_mdd
::compl_mdd(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, arg, res)
{
  key_length = 1; 
  ans_length = 1;
  // ct entry 0: input node
  // ct entry 1: output node
}

bool MEDDLY::compl_mdd::isEntryStale(const int* data)
{
  return argF->isStale(data[0]) || resF->isStale(data[1]);
}

void MEDDLY::compl_mdd::discardEntry(const int* data)
{
  argF->uncacheNode(data[0]);
  resF->uncacheNode(data[1]);
}

void MEDDLY::compl_mdd::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d): %d]", getName(), data[0], data[1]);
}

void MEDDLY::compl_mdd::compute(const dd_edge& a, dd_edge& b) 
{
  int result = compute(a.getNode());
  b.set(result, 0, resF->getNodeLevel(result));
}

int MEDDLY::compl_mdd::compute(int a)
{
  if (argF->isTerminalNode(a)) {
    return resF->getTerminalNode(a==0);
  }
  // Check compute table
  const int* cacheFind = CT->find(this, &a);
  if (cacheFind) {
    return resF->linkNode(cacheFind[1]);
  }

  const int aLevel = argF->getNodeLevel(a);
  int bSize = resF->getLevelSize(aLevel);
  int b = resF->createTempNode(aLevel, bSize, false);

  if (argF->isFullNode(a)) {
    const int aSize = argF->getFullNodeSize(a);
    DCASSERT(aSize <= bSize);
    int i, temp;
    for (i=0; i<aSize; i++) {
      temp = compute(argF->getFullNodeDownPtr(a, i));
      resF->setDownPtrWoUnlink(b, i, temp);
      resF->unlinkNode(temp);
    }
    temp = resF->getTerminalNode(true);
    for (; i<bSize; i++) {
      resF->setDownPtrWoUnlink(b, i, temp);
    }
  } // a is full
  else {
    const int aNnz = argF->getSparseNodeSize(a);
    int i, zero;
    i=0;
    zero = resF->getTerminalNode(true);
    for (int z=0; z<aNnz; z++) {
      int aIndex = argF->getSparseNodeIndex(a, z);
      for (; i<aIndex; i++) {
        resF->setDownPtrWoUnlink(b, i, zero);
      }
      int temp = compute(argF->getSparseNodeDownPtr(a, z));
      resF->setDownPtrWoUnlink(b, i, temp);
      resF->unlinkNode(temp);
    } // for z
    for (; i<bSize; i++) {
      resF->setDownPtrWoUnlink(b, i, zero);
    }
  } // a is sparse

  // reduce, save in compute table
  b = resF->reduceNode(b);
  static int cacheEntry[2];
  cacheEntry[0] = argF->cacheNode(a);
  cacheEntry[1] = resF->cacheNode(b);
  CT->add(this, cacheEntry);
  return b;
}


// ******************************************************************
// *                                                                *
// *                        compl_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_mxd : public unary_operation {
  public:
    compl_mxd(const unary_opname* oc, expert_forest* arg, expert_forest* res);

    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, dd_edge& b);

    int compute(int ht, int a);
};

MEDDLY::compl_mxd
::compl_mxd(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, arg, res)
{
  key_length = 2; 
  ans_length = 1;
  // ct entry 0: level
  // ct entry 1: input node
  // ct entry 2: output node
}

bool MEDDLY::compl_mxd::isEntryStale(const int* data)
{
  return argF->isStale(data[1]) || resF->isStale(data[2]);
}

void MEDDLY::compl_mxd::discardEntry(const int* data)
{
  argF->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void MEDDLY::compl_mxd::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d, %d): %d]", getName(), data[0], data[1], data[2]);
}

void MEDDLY::compl_mxd::compute(const dd_edge& a, dd_edge& b) 
{
  int result = compute(argF->getDomain()->getNumVariables(), a.getNode());
  b.set(result, 0, resF->getNodeLevel(result));
}

int MEDDLY::compl_mxd::compute(int ht, int a)
{
  if (0==ht) {
    return resF->getTerminalNode(a==0);
  }
  if (argF->isTerminalNode(a) && 
      resF->getReductionRule() == forest::FULLY_REDUCED) 
  {
    return resF->getTerminalNode(a==0);
  }
  // Check compute table
  static int cacheEntry[3];
  cacheEntry[0] = ht;
  cacheEntry[1] = a;
  const int* cacheFind = CT->find(this, cacheEntry);
  if (cacheFind) {
#ifdef DEBUG_MXD_COMPL
    fprintf(stderr, "\tin CT:   compl_mxd(%d, %d) : %d\n", ht, a, cacheFind[2]);
#endif
    return resF->linkNode(cacheFind[2]);
  }

#ifdef DEBUG_MXD_COMPL
  fprintf(stderr, "\tstarting compl_mxd(%d, %d)\n", ht, a);
#endif

  const int aLevel = argF->getNodeLevel(a);
  int bSize = resF->getLevelSize(ht);
  int b = resF->createTempNode(ht, bSize, false);
  int htm1 = (ht > 0) ? -ht : (-ht)-1;

  if (aLevel != ht) {
    // All result[i] = compute(htm1, a)
    // Special case: expand both levels at skipped identity-reduced level.
    if (ht > 0 &&
        resF->getReductionRule() == forest::IDENTITY_REDUCED) {
      DCASSERT(ht > 0);
      DCASSERT(htm1 < 0);
      int htm2 = (-htm1)-1;
      int zero = compute(htm2, 0);
      int temp = compute(htm2, a);
      for (int i = 0; i < bSize; ++i) {
        // Build node at nextLevel
        // such that n[i==j] = compute(owner, nextNextLevel, a)
        //       and n[i!=j] = compute(owner, nextNextLevel, 0)
        int n = resF->createTempNode(htm1, bSize, false);
        for (int j = 0; j < bSize; ++j) {
          resF->setDownPtrWoUnlink(n, j, (i == j)? temp: zero);
        }
        n = resF->reduceNode(n);
        resF->setDownPtrWoUnlink(b, i, n);
        resF->unlinkNode(n);
      }
      resF->unlinkNode(temp);
      resF->unlinkNode(zero);
    }
    else {
      // For Fully and Quasi. Also for prime ht for Identity reduced.
      int temp = compute(htm1, a);
      for (int i = 0; i < bSize; ++i) {
        resF->setDownPtrWoUnlink(b, i, temp);
      }
      resF->unlinkNode(temp);
    }
  } // skipped level
  else if (argF->isFullNode(a)) {
    // a is a full-node
    const int aSize = argF->getFullNodeSize(a);
    DCASSERT(aSize <= bSize);
    int i = 0;
    for ( ; i < aSize; ++i) {
      int temp = compute(htm1, argF->getFullNodeDownPtr(a, i));
      resF->setDownPtrWoUnlink(b, i, temp);
      resF->unlinkNode(temp);
    }
    if (i < bSize) {
      int zero = compute(htm1, 0);
      do {
        resF->setDownPtrWoUnlink(b, i++, zero);
      } while (i < bSize);
      resF->unlinkNode(zero);
    }
  } // not skipped, a full
  else {
    // a is a sparse-node
    DCASSERT(argF->isSparseNode(a));
    const int aSize = argF->getSparseNodeSize(a);
    DCASSERT(argF->getSparseNodeIndex(a, aSize - 1) < bSize);
    int i = 0;        // goes from 0..bSize
    int aIndex = 0;   // j is the sparse index; aIndex is the equiv full index
    int zero = compute(htm1, 0);
    for (int j = 0; j < aSize; ++i, ++j) {
      aIndex = argF->getSparseNodeIndex(a, j);
      while (i < aIndex) {
        // a[i] is 0
        resF->setDownPtrWoUnlink(b, i++, zero);
      }
      DCASSERT(i == aIndex && i < bSize);
      int temp = compute(htm1, argF->getSparseNodeDownPtr(a, j));
      resF->setDownPtrWoUnlink(b, i, temp);
      resF->unlinkNode(temp);
    }
    while (i < bSize) {
      // a[i] is 0
      resF->setDownPtrWoUnlink(b, i++, zero);
    }
    resF->unlinkNode(zero);
  } // not skipped, a sparse

  // reduce, save in compute table
  b = resF->reduceNode(b);
  cacheEntry[0] = ht;
  cacheEntry[1] = argF->cacheNode(a);
  cacheEntry[2] = resF->cacheNode(b);
  CT->add(this, cacheEntry);

#ifdef DEBUG_MXD_COMPL
  fprintf(stderr, "\tfinished compl_mxd(%d, %d) : %d\n", ht, a, b);
#endif

  return b;
}



// ******************************************************************
// *                                                                *
// *                       compl_opname class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_opname : public unary_opname {
  public:
    compl_opname();
    virtual unary_operation* 
      buildOperation(const forest* ar, const forest* res) const;
};

MEDDLY::compl_opname::compl_opname()
 : unary_opname("Complement")
{
}

MEDDLY::unary_operation* 
MEDDLY::compl_opname
::buildOperation(const forest* arg, const forest* res) const
{
  if (0==arg || 0==res) return 0;

  if (arg->getRangeType() != forest::BOOLEAN ||
      arg->getEdgeLabeling() != forest::MULTI_TERMINAL ||
      res->getRangeType() != forest::BOOLEAN ||
      res->getEdgeLabeling() != forest::MULTI_TERMINAL ||
      arg->isForRelations() != res->isForRelations()
  ) throw error(error::TYPE_MISMATCH);

  if (arg->isForRelations())
    return new compl_mxd(this, (expert_forest*) arg, (expert_forest*) res);
  else
    return new compl_mdd(this, (expert_forest*) arg, (expert_forest*) res);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeComplement(const settings &s)
{
  return new compl_opname;
}


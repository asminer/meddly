
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
#include "copy.h"
#include "../compute_table.h"

namespace MEDDLY {
  class copy_MT;

  class copy_bool2MT;
  class copy_MT2bool;

  class copy_MT2Evplus;

  class copy_opname;
};

// ******************************************************************
// *                                                                *
// *                         copy_MT  class                         *
// *                                                                *
// ******************************************************************

/// Abstract base class for copying between multi-terminal DDs.
class MEDDLY::copy_MT : public unary_operation {
  public:
    copy_MT(const unary_opname* oc, expert_forest* arg, expert_forest* res);

    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge &arg, dd_edge &res);
  protected:
    virtual int compute(int a) = 0;
};

MEDDLY::copy_MT
:: copy_MT(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, arg, res)
{
  key_length = 1;
  ans_length = 1;
}

bool MEDDLY::copy_MT::isEntryStale(const int* entryData)
{
  return 
    argF->isStale(entryData[0]) ||
    resF->isStale(entryData[1]);
}

void MEDDLY::copy_MT::discardEntry(const int* entryData)
{
  argF->uncacheNode(entryData[0]);
  resF->uncacheNode(entryData[1]);
}

void MEDDLY::copy_MT::showEntry(FILE* strm, const int *entryData) const
{
  fprintf(strm, "[%s(%d) %d]", getName(), entryData[0], entryData[1]);
}

void MEDDLY::copy_MT::compute(const dd_edge &arg, dd_edge &res)
{
  int result = compute(arg.getNode());
  res.set(result, 0, resF->getNodeLevel(result));
}

// ******************************************************************
// *                                                                *
// *                       copy_bool2MT class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::copy_bool2MT : public copy_MT {
  public:
    copy_bool2MT(const unary_opname* N, expert_forest* A, expert_forest* R)
      : copy_MT(N, A, R) { }
  protected:
    virtual int compute(int a);
};

int MEDDLY::copy_bool2MT::compute(int a)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    int aTerm = argF->getBoolean(a) ? 1 : 0;
    if (resF->getRangeType() == forest::INTEGER) {
      return resF->getTerminalNode(aTerm);
    } else {
      return resF->getTerminalNode(float(aTerm));
    }
  }

  // Check compute table
  const int* cacheFind = CT->find(this, &a);
  if (cacheFind) {
    return resF->linkNode(cacheFind[1]);
  }

  // Make new node for answer
  const int aLevel = argF->getNodeLevel(a);
  int bSize = resF->getLevelSize(aLevel);
  int b = resF->createTempNode(aLevel, bSize, false);
  int zero = (resF->getRangeType()==forest::INTEGER) 
    ? resF->getTerminalNode(0)
    : resF->getTerminalNode(float(0));

  // Recurse
  int i;
  if (argF->isFullNode(a)) {
    const int aSize = argF->getFullNodeSize(a);
    DCASSERT(aSize <= bSize);
    for (i=0; i<aSize; i++) {
      int temp = compute(argF->getFullNodeDownPtr(a, i));
      resF->setDownPtrWoUnlink(b, i, temp);
      resF->unlinkNode(temp);
    }
  } // a is full
  else {
    const int aSize = argF->getSparseNodeSize(a);
    DCASSERT(argF->getSparseNodeIndex(a, aSize - 1) < bSize);
    i = 0;
    for (int z=0; z<aSize; z++) {
      int aIndex = argF->getSparseNodeIndex(a, z);
      for (; i<aIndex; i++) {
        resF->setDownPtrWoUnlink(b, i, zero);
      }
      DCASSERT(i == aIndex && i < bSize);
      int temp = compute(argF->getSparseNodeDownPtr(a, z));
      resF->setDownPtrWoUnlink(b, i, temp);
      resF->unlinkNode(temp);
    }
  } // a is sparse
  for (; i<bSize; i++) {
    resF->setDownPtrWoUnlink(b, i, zero);
  }

  // Reduce, add to compute table
  b = resF->reduceNode(b);
  static int cacheEntry[2];
  cacheEntry[0] = argF->cacheNode(a);
  cacheEntry[1] = resF->cacheNode(b);
  CT->add(this, cacheEntry);

  return b;
}

// ******************************************************************
// *                                                                *
// *                       copy_MT2bool class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::copy_MT2bool : public copy_MT {
  public:
    copy_MT2bool(const unary_opname* N, expert_forest* A, expert_forest* R)
      : copy_MT(N, A, R) { }
  protected:
    virtual int compute(int a);
};

int MEDDLY::copy_MT2bool::compute(int a)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    bool aTerm = (argF->getRangeType() == forest::INTEGER)
      ? (argF->getInteger(a) != 0)
      : (argF->getReal(a) !=0);
    return resF->getTerminalNode(aTerm);
  }

  // Check compute table
  const int* cacheFind = CT->find(this, &a);
  if (cacheFind) {
    return resF->linkNode(cacheFind[1]);
  }

  // Make new node for answer
  const int aLevel = argF->getNodeLevel(a);
  int bSize = resF->getLevelSize(aLevel);
  int b = resF->createTempNode(aLevel, bSize, false);
  int zero = resF->getTerminalNode(false);

  // Recurse
  int i;
  if (argF->isFullNode(a)) {
    const int aSize = argF->getFullNodeSize(a);
    DCASSERT(aSize <= bSize);
    for (i=0; i<aSize; i++) {
      int temp = compute(argF->getFullNodeDownPtr(a, i));
      resF->setDownPtrWoUnlink(b, i, temp);
      resF->unlinkNode(temp);
    }
  } // a is full
  else {
    const int aSize = argF->getSparseNodeSize(a);
    DCASSERT(argF->getSparseNodeIndex(a, aSize - 1) < bSize);
    i = 0;
    for (int z=0; z<aSize; z++) {
      int aIndex = argF->getSparseNodeIndex(a, z);
      for (; i<aIndex; i++) {
        resF->setDownPtrWoUnlink(b, i, zero);
      }
      DCASSERT(i == aIndex && i < bSize);
      int temp = compute(argF->getSparseNodeDownPtr(a, z));
      resF->setDownPtrWoUnlink(b, i, temp);
      resF->unlinkNode(temp);
      i++;
    }
  } // a is sparse
  for (; i<bSize; i++) {
    resF->setDownPtrWoUnlink(b, i, zero);
  }

  // Reduce, add to compute table
  b = resF->reduceNode(b);
  static int cacheEntry[2];
  cacheEntry[0] = argF->cacheNode(a);
  cacheEntry[1] = resF->cacheNode(b);
  CT->add(this, cacheEntry);

  return b;
}

// ******************************************************************
// *                                                                *
// *                      copy_MT2Evplus class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::copy_MT2Evplus : public unary_operation {
  public:
    copy_MT2Evplus(const unary_opname* oc, expert_forest* arg, 
      expert_forest* res);

    virtual bool isEntryStale(const int* entryData) {
      return 
        argF->isStale(entryData[0]) ||
        resF->isStale(entryData[1]);
    }
    virtual void discardEntry(const int* entryData) {
      argF->uncacheNode(entryData[0]);
      resF->uncacheNode(entryData[1]);
    }
    virtual void showEntry(FILE* strm, const int *entryData) const {
      fprintf(strm, "[%s(%d) <%d, %d>]", getName(), entryData[0], 
        entryData[2], entryData[1]);
    }
    virtual void compute(const dd_edge &arg, dd_edge &res) {
      int b, bev;
      compute(arg.getNode(), b, bev);
      res.set(b, bev, resF->getNodeLevel(b));
    }
    virtual void compute(int a, int &b, int &bev);
};

MEDDLY::copy_MT2Evplus::copy_MT2Evplus(const unary_opname* oc, 
  expert_forest* arg, expert_forest* res) : unary_operation(oc, arg, res)
{
  key_length = 1;
  ans_length = 2;
  // entry[0]: MT node
  // entry[1]: EV node
  // entry[2]: EV value
}

void MEDDLY::copy_MT2Evplus::compute(int a, int &b, int &bev)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    bev = argF->getInteger(a);
    b = resF->getTerminalNode(true);
    return;
  }

  // Check compute table
  const int* cacheFind = CT->find(this, &a);
  if (cacheFind) {
    b = resF->linkNode(cacheFind[1]);
    bev = cacheFind[2];
    return;
  }

  // Make new node for answer
  const int aLevel = argF->getNodeLevel(a);
  int bSize = resF->getLevelSize(aLevel);
  b = resF->createTempNode(aLevel, bSize, false);
  int zero, zev;
  compute(0, zero, zev);

  // Recurse
  int i;
  if (argF->isFullNode(a)) {
    const int aSize = argF->getFullNodeSize(a);
    DCASSERT(aSize <= bSize);
    for (i=0; i<aSize; i++) {
      int d, dev;
      compute(argF->getFullNodeDownPtr(a, i), d, dev);
      resF->setDownPtrWoUnlink(b, i, d);
      resF->unlinkNode(d);
      resF->setEdgeValue(b, i, dev);
    }
  } // a is full
  else {
    const int aSize = argF->getSparseNodeSize(a);
    DCASSERT(argF->getSparseNodeIndex(a, aSize - 1) < bSize);
    i = 0;
    for (int z=0; z<aSize; z++) {
      int aIndex = argF->getSparseNodeIndex(a, z);
      for (; i<aIndex; i++) {
        resF->setDownPtrWoUnlink(b, i, zero);
        resF->setEdgeValue(b, i, zev);
      }
      DCASSERT(i == aIndex && i < bSize);
      int d, dev;
      compute(argF->getSparseNodeDownPtr(a, z), d, dev);
      resF->setDownPtrWoUnlink(b, i, d);
      resF->unlinkNode(d);
      resF->setEdgeValue(b, i, dev);
    }
  } // a is sparse
  for (; i<bSize; i++) {
    resF->setDownPtrWoUnlink(b, i, zero);
    resF->setEdgeValue(b, i, zev);
  }

  // Reduce, add to compute table
  bev = 0;
  resF->normalizeAndReduceNode(b, bev);
  static int cacheEntry[3];
  cacheEntry[0] = argF->cacheNode(a);
  cacheEntry[1] = resF->cacheNode(b);
  cacheEntry[2] = bev;
  CT->add(this, cacheEntry);
}

// ******************************************************************
// *                                                                *
// *                       copy_opname  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::copy_opname : public unary_opname {
  public:
    copy_opname();
    virtual unary_operation* 
      buildOperation(expert_forest* ar, expert_forest* res) const;
};

MEDDLY::copy_opname::copy_opname()
 : unary_opname("Copy")
{
}

MEDDLY::unary_operation*
MEDDLY::copy_opname
::buildOperation(expert_forest* arg, expert_forest* res) const
{
  if (0==arg || 0==res) return 0;

  if (arg->getDomain() != res->getDomain())
    throw error(error::DOMAIN_MISMATCH);

  if (arg->isForRelations() != res->isForRelations())
    throw error(error::TYPE_MISMATCH);

  if (arg->getEdgeLabeling() != forest::MULTI_TERMINAL)
    throw error(error::NOT_IMPLEMENTED);

  if (res->getEdgeLabeling() != forest::MULTI_TERMINAL) {
    if (arg->getRangeType() != res->getRangeType())
      throw error(error::TYPE_MISMATCH);
   
    if (res->getEdgeLabeling() == forest::EVPLUS) {
      return new copy_MT2Evplus(this,  arg,  res);
    }

    throw error(error::NOT_IMPLEMENTED);
  }

  if (arg->getRangeType() == forest::BOOLEAN) {
    if (res->getRangeType() == forest::BOOLEAN)
      throw error(error::NOT_IMPLEMENTED);
    return new copy_bool2MT(this,  arg,  res);
  } // boolean
  else {
    if (res->getRangeType() != forest::BOOLEAN)
      throw error(error::NOT_IMPLEMENTED);
    return new copy_MT2bool(this,  arg,  res);
  }
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeCopy(const settings &s)
{
  return new copy_opname;
}


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

    virtual bool isStaleEntry(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge &arg, dd_edge &res);
  protected:
    virtual int compute(int a) = 0;

    inline bool findResult(int a, int &b) {
      CTsrch.key(0) = a;
      const int* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      b = resF->linkNode(cacheFind[1]);
      return true;
    }
    inline int saveResult(int a, int b) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      entry.key(0) = argF->cacheNode(a);
      entry.result(0) = resF->cacheNode(b);
      CT->addEntry();
      return b;
    }
};

MEDDLY::copy_MT
:: copy_MT(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, 1, 1, arg, res)
{
}

bool MEDDLY::copy_MT::isStaleEntry(const int* entryData)
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
    int compute(int in, int a);
};

int MEDDLY::copy_bool2MT::compute(int a)
{
  return compute(-1, a);
}

int MEDDLY::copy_bool2MT::compute(int in, int a)
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

  // See if we can ignore "in"
  if (!resF->isIdentityReduced()) in = -1;

  // Check compute table
  int b;
  if (findResult(a, b)) return b;

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  const int size = resF->getLevelSize(level);
  node_builder& nb = resF->useNodeBuilder(level, size);

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, true);

  // recurse
  for (int i=0; i<size; i++) {
    nb.d(i) = compute(i, A->d(i));
  }

  // Cleanup
  node_reader::recycle(A);

  // Reduce
  b = resF->createReducedNode(in, nb);

  // Add to compute table
  return saveResult(a, b);
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
    int compute(int in, int a);
    int compute(int in, int k, int a);
};

int MEDDLY::copy_MT2bool::compute(int a)
{
  if (resF->isQuasiReduced())
    return compute(-1, resF->getNumVariables(), a);
  else 
    return compute(-1, a);
}

int MEDDLY::copy_MT2bool::compute(int in, int a)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    bool aTerm = (argF->getRangeType() == forest::INTEGER)
      ? (argF->getInteger(a) != 0)
      : (argF->getReal(a) !=0);
    return resF->getTerminalNode(aTerm);
  }

  // See if we can ignore "in"
  if (!resF->isIdentityReduced()) in = -1;

  // Check compute table
  int b;
  if (findResult(a, b)) return b;

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  const int size = resF->getLevelSize(level);
  node_builder& nb = resF->useNodeBuilder(level, size);

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, true);

  // recurse
  for (int i=0; i<size; i++) {
    nb.d(i) = compute(i, A->d(i));
  }

  // Cleanup
  node_reader::recycle(A);

  // Reduce
  b = resF->createReducedNode(in, nb);

  // Add to compute table
  return saveResult(a, b);
}

int MEDDLY::copy_MT2bool::compute(int in, int k, int a)
{
  // Check terminals
  if (0==k) {
    bool aTerm = (argF->getRangeType() == forest::INTEGER)
      ? (argF->getInteger(a) != 0)
      : (argF->getReal(a) !=0);
    return resF->getTerminalNode(aTerm);
  }

  // See if we can ignore "in"
  if (resF->isIdentityReduced()) {
    if (k>0) in = -1;
  } else {
    in = -1;
  }

  // Get level number
  const int aLevel = argF->getNodeLevel(a);

  // Check compute table
  int b;
  if (k == aLevel) if (findResult(a, b)) return b;
  int nextk;
  if (resF->isForRelations()) {
    nextk = (k>0) ? -k : -k-1;
  } else {
    nextk = k-1;
  }

  // Initialize node builder
  const int size = resF->getLevelSize(k);
  node_builder& nb = resF->useNodeBuilder(k, size);

  // Initialize node reader
  node_reader* A;
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      A = argF->initIdentityReader(k, in, a, true);
    } else {
      A = argF->initRedundantReader(k, a, true);
    }
  } else {
    A = argF->initNodeReader(a, true);
  }

  // recurse
  for (int i=0; i<size; i++) {
    nb.d(i) = compute(i, nextk, A->d(i));
  }

  // Cleanup
  node_reader::recycle(A);

  // Reduce
  b = resF->createReducedNode(in, nb);

  // Add to compute table
  if (k == aLevel) saveResult(a, b);
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

    virtual bool isStaleEntry(const int* entryData) {
      return 
        argF->isStale(entryData[0]) ||
        resF->isStale(entryData[2]);
    }
    virtual void discardEntry(const int* entryData) {
      argF->uncacheNode(entryData[0]);
      resF->uncacheNode(entryData[2]);
    }
    virtual void showEntry(FILE* strm, const int *entryData) const {
      fprintf(strm, "[%s(%d) <%d, %d>]", getName(), entryData[0], 
        entryData[1], entryData[2]);
    }
    virtual void compute(const dd_edge &arg, dd_edge &res) {
      int b, bev;
      compute(arg.getNode(), b, bev);
      res.set(b, bev, resF->getNodeLevel(b));
    }
    virtual void compute(int a, int &b, int &bev);
};

MEDDLY::copy_MT2Evplus::copy_MT2Evplus(const unary_opname* oc, 
  expert_forest* arg, expert_forest* res)
: unary_operation(oc, 1, 2, arg, res)
{
  // entry[0]: MT node
  // entry[1]: EV value
  // entry[2]: EV node
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
  CTsrch.key(0) = a;
  const int* cacheFind = CT->find(CTsrch);
  if (cacheFind) {
    bev = cacheFind[1];
    b = resF->linkNode(cacheFind[2]);
    return;
  }

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  const int size = resF->getLevelSize(level);
  node_builder& nb = resF->useNodeBuilder(level, size);

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, true);

  // recurse
  for (int i=0; i<size; i++) {
    int d, dev;
    compute(A->d(i), d, dev);
    nb.d(i) = d;
    nb.ei(i) = dev;
  }

  // Cleanup
  node_reader::recycle(A);

  // Reduce
  resF->createReducedNode(-1, nb, bev, b);

  // Add to compute table
  compute_table::temp_entry &entry = CT->startNewEntry(this);
  entry.key(0) = argF->cacheNode(a);
  entry.result(0) = bev;
  entry.result(1) = resF->cacheNode(b);
  CT->addEntry();
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


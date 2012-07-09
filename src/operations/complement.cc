
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

    virtual bool isStaleEntry(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, dd_edge& b);

  protected:
    int compute(int a);

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

MEDDLY::compl_mdd
::compl_mdd(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, 1, 1, arg, res)
{
  // ct entry 0: input node
  // ct entry 1: output node
}

bool MEDDLY::compl_mdd::isStaleEntry(const int* data)
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
  // Check terminals
  if (argF->isTerminalNode(a)) {
    return resF->getTerminalNode(a==0);
  }

  // Check compute table
  int b;
  if (findResult(a, b)) return b;

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  const int size = resF->getLevelSize(level);
  expert_forest::nodeBuilder& nb = resF->useNodeBuilder(level, size);

  // Initialize node reader
  expert_forest::nodeReader* A = argF->initNodeReader(a);

  // recurse
  for (int i=0; i<size; i++) {
    nb.d(i) = compute((*A)[i]);
  }

  // Cleanup
  argF->recycle(A);

  // Reduce
  b = resF->createReducedNode(-1, nb);

  // Add to compute table
  return saveResult(a, b);
}


// ******************************************************************
// *                                                                *
// *                        compl_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_mxd : public unary_operation {
  public:
    compl_mxd(const unary_opname* oc, expert_forest* arg, expert_forest* res);

    virtual bool isStaleEntry(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, dd_edge& b);

    int compute(int in, int k, int a);
};

MEDDLY::compl_mxd
::compl_mxd(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, 2, 1, arg, res)
{
  // ct entry 0: level
  // ct entry 1: input node
  // ct entry 2: output node
}

bool MEDDLY::compl_mxd::isStaleEntry(const int* data)
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
  int result = compute(-1, argF->getDomain()->getNumVariables(), a.getNode());
  b.set(result, 0, resF->getNodeLevel(result));
}

int MEDDLY::compl_mxd::compute(int in, int k, int a)
{
  if (0==k) {
    return resF->getTerminalNode(a==0);
  }
  if (argF->isTerminalNode(a) && 
      resF->isFullyReduced()) 
  {
    return resF->getTerminalNode(a==0);
  }
  // Check compute table
  CTsrch.key(0) = k;
  CTsrch.key(1) = a;
  const int* cacheFind = CT->find(CTsrch);
  if (cacheFind) {
#ifdef DEBUG_MXD_COMPL
    fprintf(stderr, "\tin CT:   compl_mxd(%d, %d) : %d\n", ht, a, cacheFind[2]);
#endif
    return resF->linkNode(cacheFind[2]);
  }

#ifdef DEBUG_MXD_COMPL
  fprintf(stderr, "\tstarting compl_mxd(%d, %d)\n", ht, a);
#endif

  // Initialize node builder
  const int size = resF->getLevelSize(k);
  expert_forest::nodeBuilder& nb = resF->useNodeBuilder(k, size);

  // Initialize node reader
  const int aLevel = argF->getNodeLevel(a);
  MEDDLY_DCASSERT(!isLevelAbove(aLevel, k));
  expert_forest::nodeReader* A;
  bool canSave = true;
  if (aLevel == k) {
    A = argF->initNodeReader(a);
  } else if (k>0 || argF->isFullyReduced()) {
    A = argF->initRedundantReader(k, a);
  } else {
    MEDDLY_DCASSERT(in>=0);
    A = argF->initIdentityReader(k, in, a);
    canSave = false;
  }
  
  // recurse
  int nextLevel = (k>0) ? -k : -k-1;
  int nnz = 0;
  for (int i=0; i<size; i++) {
    int d = compute(i, nextLevel, (*A)[i]);
    nb.d(i) = d;
    if (d) nnz++;
  }

  // cleanup
  argF->recycle(A);

  // reduce, save in CT
  int result = resF->createReducedNode(in, nb);
  if (k<0 && 1==nnz) canSave = false;
  if (canSave) {
    compute_table::temp_entry &entry = CT->startNewEntry(this);
    entry.key(0) = k;
    entry.key(1) = argF->cacheNode(a);
    entry.result(0) = resF->cacheNode(result);
    CT->addEntry();
  }
  return result;
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
      buildOperation(expert_forest* ar, expert_forest* res) const;
};

MEDDLY::compl_opname::compl_opname()
 : unary_opname("Complement")
{
}

MEDDLY::unary_operation* 
MEDDLY::compl_opname
::buildOperation(expert_forest* arg, expert_forest* res) const
{
  if (0==arg || 0==res) return 0;

  if (arg->getDomain() != res->getDomain())
    throw error(error::DOMAIN_MISMATCH);

  if (arg->getRangeType() != forest::BOOLEAN ||
      arg->getEdgeLabeling() != forest::MULTI_TERMINAL ||
      res->getRangeType() != forest::BOOLEAN ||
      res->getEdgeLabeling() != forest::MULTI_TERMINAL ||
      arg->isForRelations() != res->isForRelations()
  ) throw error(error::TYPE_MISMATCH);

  if (arg->isForRelations())
    return new compl_mxd(this,  arg,  res);
  else
    return new compl_mdd(this,  arg,  res);
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


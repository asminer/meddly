
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
#include "mdd2evplus.h"

// #define TRACE_ALL_OPS

namespace MEDDLY {
  class mdd2evplus_operation;
  class mdd2evplus_opname;
};

// ******************************************************************
// *                                                                *
// *                   mdd2evplus_operation class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::mdd2evplus_operation : public unary_operation {
  public:
    mdd2evplus_operation(const unary_opname* oc, expert_forest* arg, 
      expert_forest* res);

    virtual bool isStaleEntry(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;

    virtual void compute(const dd_edge &arg, dd_edge &res);

    void compute(int k, int a, int &bdn, int &bcard);
};

MEDDLY::mdd2evplus_operation::mdd2evplus_operation(const unary_opname* oc, 
  expert_forest* arg, expert_forest* res)
: unary_operation(oc, 1, 2, arg, res)
{
  // answer[0] : pointer
  // answer[1] : cardinality
}

bool 
MEDDLY::mdd2evplus_operation
::isStaleEntry(const int* entryData)
{
  return 
    argF->isStale(entryData[0]) ||
    resF->isStale(entryData[1]);
}

void 
MEDDLY::mdd2evplus_operation
::discardEntry(const int* entryData)
{
  argF->uncacheNode(entryData[0]);
  resF->uncacheNode(entryData[1]);
}

void 
MEDDLY::mdd2evplus_operation
::showEntry(FILE* strm, const int *entryData) const
{
  fprintf(strm, "[%s %d %d (card %d)]", getName(), entryData[0], 
    entryData[1], entryData[2]);
}

void
MEDDLY::mdd2evplus_operation
::compute(const dd_edge &arg, dd_edge &res)
{
  MEDDLY_DCASSERT(arg.getForest() == argF);
  MEDDLY_DCASSERT(res.getForest() == resF);
  MEDDLY_DCASSERT(resF->getTerminalNode(false) == 0);
  MEDDLY_DCASSERT(argF->getTerminalNode(true) < 0);
  MEDDLY_DCASSERT(argF->getTerminalNode(false) == 0);
  int down, card;
  int nVars = argF->getDomain()->getNumVariables();
  compute(nVars, arg.getNode(), down, card);
  res.set(down, 0);
}

void
MEDDLY::mdd2evplus_operation
::compute(int k, int a, int &bdn, int &bcard)
{
  // Deal with terminals
  if (0 == a) {
    bdn = 0;
    bcard = 0;
    return;
  }
  if (0 == k) {
    bdn = resF->getTerminalNode(true);
    bcard = 1;
    return;
  }

  int aLevel = argF->getNodeLevel(a);
  MEDDLY_DCASSERT(aLevel <= k);

  // Check compute table
  if (aLevel == k) {
    CTsrch.key(0) = a;
    const int* cacheEntry = CT->find(CTsrch);
    if (cacheEntry) {
      bdn = resF->linkNode(cacheEntry[1]);
      bcard = cacheEntry[2];
      return;
    }
  }

#ifdef TRACE_ALL_OPS
  printf("calling mdd2evplus::compute(%d, %d)\n", height, a);
#endif

  // Initialize node builder
  const int size = resF->getLevelSize(k);
  node_builder& nb = resF->useNodeBuilder(k, size);
  
  // Initialize node reader
  node_reader* A = (aLevel < k)
    ? argF->initRedundantReader(k, a, true)
    : argF->initNodeReader(a, true);

  // recurse
  bcard = 0;
  for (int i=0; i<size; i++) {
    int ddn, dcard;
    compute(k-1, A->d(i), ddn, dcard);
    nb.d(i) = ddn;
    if (nb.d(i)) {
      nb.ei(i) = bcard;
      bcard += dcard;
    } else {
      MEDDLY_DCASSERT(0 == dcard);
      nb.ei(i) = 0;
    }
  }

  // Cleanup
  node_reader::recycle(A);

  // Reduce
  nb.uh(0) = bcard;
  int dummy;
  long bl;
  resF->createReducedNode(-1, nb, dummy, bl);
  bdn = bl;
  MEDDLY_DCASSERT(0==dummy);

  // Add to compute table
  compute_table::temp_entry &entry = CT->startNewEntry(this);
  entry.key(0) = argF->cacheNode(a);
  entry.result(0) = resF->cacheNode(bdn);
  entry.result(1) = bcard;
  CT->addEntry();
}

// ******************************************************************
// *                                                                *
// *                    mdd2evplus_opname  class                    *
// *                                                                *
// ******************************************************************

class MEDDLY::mdd2evplus_opname : public unary_opname {
  public:
    mdd2evplus_opname();
    virtual unary_operation* 
      buildOperation(expert_forest* ar, expert_forest* res) const;
};

MEDDLY::mdd2evplus_opname::mdd2evplus_opname()
 : unary_opname("ConvertToIndexSet")
{
}

MEDDLY::unary_operation*
MEDDLY::mdd2evplus_opname
::buildOperation(expert_forest* arg, expert_forest* res) const
{
  if (0==arg || 0==res) return 0;

  if (arg->getDomain() != res->getDomain())
    throw error(error::DOMAIN_MISMATCH);

  if (arg->isForRelations() || 
      arg->getRangeType() != forest::BOOLEAN ||
      arg->getEdgeLabeling() != forest::MULTI_TERMINAL ||
      res->isForRelations() ||
      res->getRangeType() != forest::INTEGER ||
      res->getEdgeLabeling() != forest::EVPLUS
  ) throw error(error::TYPE_MISMATCH);

  return new mdd2evplus_operation(
    this,  arg,  res
  );
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeMDD2EVPLUS(const settings &s)
{
  return new mdd2evplus_opname;
}


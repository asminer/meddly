
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
#include "../compute_table.h"

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

    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;

    virtual void compute(const dd_edge &arg, dd_edge &res);

    void compute(int ht, int a, int &bdn, int &bcard);
};

MEDDLY::mdd2evplus_operation::mdd2evplus_operation(const unary_opname* oc, 
  expert_forest* arg, expert_forest* res)
: unary_operation(oc, true, arg, res)
{
  key_length = 1;
  ans_length = 2; // pointer, cardinality
}

bool 
MEDDLY::mdd2evplus_operation
::isEntryStale(const int* entryData)
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
  DCASSERT(arg.getForest() == argF);
  DCASSERT(res.getForest() == resF);
  DCASSERT(resF->getTerminalNode(false) == 0);
  DCASSERT(argF->getTerminalNode(true) == -1);
  DCASSERT(argF->getTerminalNode(false) == 0);
  int down, card;
  int nVars = argF->getDomain()->getNumVariables();
  compute(nVars, arg.getNode(), down, card);
  res.set(down, 0, resF->getNodeLevel(down));
}

void
MEDDLY::mdd2evplus_operation
::compute(int height, int a, int &bdn, int &bcard)
{
  // Deal with terminals
  if (a == 0) {
    bdn = 0;
    bcard = 0;
    return;
  }
  if (height == 0) {
    DCASSERT(argF->getTerminalNode(true) == a);
    bdn = resF->getTerminalNode(true);
    bcard = 1;
    return;
  }

  int aHeight = argF->getNodeHeight(a);
  DCASSERT(aHeight <= height);

  // Check compute table
  if (aHeight == height) {
    const int* cacheEntry = CT->find(this, &a);
    if (cacheEntry) {
      bdn = cacheEntry[1];
      bcard = cacheEntry[2];
      resF->linkNode(bdn);
      return;
    }
  }

  // Create node at appropriate height 
  const int resultSize = resF->getLevelSize(height);
  bdn = resF->createTempNode(height, resultSize, true);
  bcard = 0;

  if (aHeight < height) {
    int ddn, dcard;
    compute(height-1, a, ddn, dcard);
    if (ddn) {
      for (int i=0; i<resultSize; i++) {
        resF->setDownPtrWoUnlink(bdn, i, ddn);
        resF->setEdgeValue(bdn, i, bcard);
        bcard += dcard;
      } // for i
      resF->unlinkNode(ddn);
    } // if ddn
  } // aHeight < height
  else {
    if (argF->isFullNode(a)) {
      int aSize = argF->getFullNodeSize(a);
      for (int i=0; i<aSize; i++) {
        int ddn, dcard;
        compute(height-1, argF->getFullNodeDownPtr(a, i), ddn, dcard);
        if (0==ddn) continue;
        resF->setDownPtrWoUnlink(bdn, i, ddn);
        resF->setEdgeValue(bdn, i, bcard);
        resF->unlinkNode(ddn);
        bcard += dcard;
      } // for i
    } // a is full
    else {
      int aNnz = argF->getSparseNodeSize(a);
      for (int i=0; i<aNnz; i++) {
        int ddn, dcard;
        compute(height-1, argF->getSparseNodeDownPtr(a, i), ddn, dcard);
        if (0==ddn) continue;
        int ii = argF->getSparseNodeIndex(a, i);
        resF->setDownPtrWoUnlink(bdn, ii, ddn);
        resF->setEdgeValue(bdn, ii, bcard);
        resF->unlinkNode(ddn);
        bcard += dcard;
      } // for i
    } // a is sparse
  } // aHeight == height

  // reduce, save in compute table
  int dummy;
  resF->normalizeAndReduceNode(bdn, dummy);
  DCASSERT((bdn==0 && dummy==INF) || (bdn!= 0 && dummy==0));
  static int cacheEntry[3];
  cacheEntry[0] = argF->cacheNode(a);
  cacheEntry[1] = resF->cacheNode(bdn);
  cacheEntry[2] = bcard;
  CT->add(this, cacheEntry);
  resF->setIndexSetCardinality(bdn, bcard);
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


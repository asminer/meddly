
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
#include "cross.h"
#include "../compute_table.h"

namespace MEDDLY {
  class cross_bool;
  class cross_opname;
};

// ******************************************************************
// *                                                                *
// *                        cross_bool class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::cross_bool : public binary_operation {
  public:
    cross_bool(const binary_opname* oc, expert_forest* a1,
      expert_forest* a2, expert_forest* res);

    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    int compute_pr(int ht, int a, int b);
    int compute_un(int ht, int a, int b);
  protected:
    compute_table::search_key CTsrch;
};

MEDDLY::cross_bool::cross_bool(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) 
: binary_operation(oc, true, a1, a2, res)
{
  key_length = 3; 
  ans_length = 1;
  // data[0] : height
  // data[1] : a
  // data[2] : b
  // data[3] : c
  CT->initializeSearchKey(CTsrch, this);
}

bool MEDDLY::cross_bool::isEntryStale(const int* data)
{
  // data[0] is the level number
  return arg1F->isStale(data[1]) ||
         arg2F->isStale(data[2]) ||
         resF->isStale(data[3]);
}

void MEDDLY::cross_bool::discardEntry(const int* data)
{
  // data[0] is the level number
  arg1F->uncacheNode(data[1]);
  arg2F->uncacheNode(data[2]);
  resF->uncacheNode(data[3]);
}

void
MEDDLY::cross_bool ::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d, %d, %d): %d]", 
    getName(), data[0], data[1], data[2], data[3]
  );
}

void
MEDDLY::cross_bool::compute(const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  int L = arg1F->getDomain()->getNumVariables();
  int cnode = compute_un(L, a.getNode(), b.getNode());
  c.set(cnode, 0, resF->getNodeLevel(cnode));
}

int MEDDLY::cross_bool::compute_un(int lh, int a, int b)
{
  if (0==lh) {
    return resF->getTerminalNode(
      arg1F->getBoolean(a)
    );
  }
  if (0==a || 0==b) return resF->getTerminalNode(0);

  // check compute table
  CTsrch.key(0) = lh;
  CTsrch.key(1) = a;
  CTsrch.key(2) = b;
  const int* cacheFind = CT->find(CTsrch);
  if (cacheFind) {
    return resF->linkNode(cacheFind[3]);
  }

  // build new result node
  int c = resF->createTempNodeMaxSize(lh, true);

  // recurse
  if (arg1F->getNodeHeight(a) < lh) {
    // skipped level
    int d = compute_pr(lh, a, b);
    for (int i=resF->getLevelSize(lh)-1; i>=0; i--) {
      resF->setDownPtrWoUnlink(c, i, d);
    } 
    resF->unlinkNode(d);
  } else {
    // not skipped level
    if (arg1F->isFullNode(a)) {
      // Full storage
      for (int i=arg1F->getFullNodeSize(a)-1; i>=0; i--) {
        int ai = arg1F->getFullNodeDownPtr(a, i);
        if (0==ai) continue;
        int d = compute_pr(lh, ai, b);
        resF->setDownPtrWoUnlink(c, i, d);
        resF->unlinkNode(d);
      }
    } else {
      // Sparse storage
      for (int z=arg1F->getSparseNodeSize(a)-1; z>=0; z--) {
        int i = arg1F->getSparseNodeIndex(a, z);
        int ai = arg1F->getSparseNodeDownPtr(a, z);
        int d = compute_pr(lh, ai, b);
        resF->setDownPtrWoUnlink(c, i, d);
        resF->unlinkNode(d);
      }
    }
  }

  // reduce, save in compute table
  c = resF->reduceNode(c);

  compute_table::temp_entry &entry = CT->startNewEntry(this);
  entry.key(0) = lh;
  entry.key(1) = arg1F->cacheNode(a);
  entry.key(2) = arg2F->cacheNode(b);
  entry.result(0) = resF->cacheNode(c);
  CT->addEntry();

  return c;
}

int MEDDLY::cross_bool::compute_pr(int ht, int a, int b)
{
  if (0==a || 0==b) return resF->getTerminalNode(0);

  // convert height to level handle
  int lh = -ht; 
  ht--;

  // check compute table
  CTsrch.key(0) = lh;
  CTsrch.key(1) = a;
  CTsrch.key(2) = b;
  const int* cacheFind = CT->find(CTsrch);
  if (cacheFind) {
    return resF->linkNode(cacheFind[3]);
  }

  // build new result node
  int c = resF->createTempNodeMaxSize(lh, true);

  // recurse
  if (arg2F->getNodeHeight(b) <= ht) {
    // skipped level
    int d = compute_un(ht, a, b);
    for (int i=resF->getLevelSize(lh)-1; i>=0; i--) {
      resF->setDownPtrWoUnlink(c, i, d);
    } 
    resF->unlinkNode(d);
  } else {
    // not skipped level
    if (arg2F->isFullNode(b)) {
      // Full storage
      for (int i=arg2F->getFullNodeSize(b)-1; i>=0; i--) {
        int bi = arg2F->getFullNodeDownPtr(b, i);
        if (0==bi) continue;
        int d = compute_un(ht, a, bi);
        resF->setDownPtrWoUnlink(c, i, d);
        resF->unlinkNode(d);
      }
    } else {
      // Sparse storage
      for (int z=arg2F->getSparseNodeSize(b)-1; z>=0; z--) {
        int i = arg2F->getSparseNodeIndex(b, z);
        int bi = arg2F->getSparseNodeDownPtr(b, z);
        int d = compute_un(ht, a, bi);
        resF->setDownPtrWoUnlink(c, i, d);
        resF->unlinkNode(d);
      }
    }
  }

  // reduce, save in compute table
  c = resF->reduceNode(c);

  compute_table::temp_entry &entry = CT->startNewEntry(this);
  entry.key(0) = lh;
  entry.key(1) = arg1F->cacheNode(a);
  entry.key(2) = arg2F->cacheNode(b);
  entry.result(0) = resF->cacheNode(c);
  CT->addEntry();

  return c;
}


// ******************************************************************
// *                                                                *
// *                       cross_opname class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::cross_opname : public binary_opname {
  public:
    cross_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::cross_opname::cross_opname()
 : binary_opname("Cross")
{
}

MEDDLY::binary_operation* 
MEDDLY::cross_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH);

  if (
    a1->isForRelations()  ||
    (a1->getRangeType() != forest::BOOLEAN) ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    a2->isForRelations()  ||
    (a2->getRangeType() != forest::BOOLEAN) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (!r->isForRelations())  ||
    (r->getRangeType() != forest::BOOLEAN) ||
    (r->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH);

  return new cross_bool(this, a1, a2, r);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeCross(const settings &s)
{
  return new cross_opname;
}


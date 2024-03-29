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

#include "../defines.h"
#include "prepostplus.h"
#include "apply_base.h"

namespace MEDDLY {
    class prepostplus_evplus;

    class preplus_evplus;
    class postplus_evplus;

    binary_list PREPLUS_cache;
    binary_list POSTPLUS_cache;
};

// ******************************************************************
// *                                                                *
// *                    prepostplus_evplus  class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::prepostplus_evplus : public generic_binary_evplus {
  public:
    prepostplus_evplus(binary_list &c, forest* arg1, forest* arg2, forest* res);

  protected:
    virtual ct_entry_key* findResult(long aev, node_handle a,
      long bev, node_handle b, long& cev, node_handle &c);
    virtual void saveResult(ct_entry_key* Key,
      long aev, node_handle a, long bev, node_handle b, long cev, node_handle c);

    virtual bool checkTerminals(long aev, node_handle a, long bev, node_handle b,
      long& cev, node_handle& c);
};

MEDDLY::prepostplus_evplus::prepostplus_evplus(binary_list &cache,
    forest* arg1, forest* arg2, forest* res)
  : generic_binary_evplus(cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);

    MEDDLY_DCASSERT(arg1->isForRelations());
    MEDDLY_DCASSERT(!arg2->isForRelations());
}

MEDDLY::ct_entry_key* MEDDLY::prepostplus_evplus::findResult(long aev, node_handle a,
  long bev, node_handle b, long& cev, node_handle &c)
{
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  MEDDLY_DCASSERT(!can_commute);
  CTsrch->writeL(0);
  CTsrch->writeN(a);
  CTsrch->writeL(0);
  CTsrch->writeN(b);

  CT0->find(CTsrch, CTresult[0]);
  if (!CTresult[0]) return CTsrch;
  cev = CTresult[0].readL();
  c = resF->linkNode(CTresult[0].readN());
  if (c != 0) {
    cev += aev + bev;
  }
  CT0->recycle(CTsrch);
  return 0;
}

void MEDDLY::prepostplus_evplus::saveResult(ct_entry_key* Key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  CTresult[0].reset();
  CTresult[0].writeL(c == 0 ? 0L : cev - aev - bev);
  CTresult[0].writeN(c);
  CT0->addEntry(Key, CTresult[0]);
}

bool MEDDLY::prepostplus_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == -1 && b == -1) {
    c = -1;
    cev = aev + bev;
    MEDDLY_DCASSERT(cev >= 0);
    return true;
  }
  if (a == 0 || b == 0) {
    c = 0;
    cev = 0;
    return true;
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                     preplus_evplus  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::preplus_evplus : public prepostplus_evplus {
  public:
    preplus_evplus(forest* arg1, forest* arg2, forest* res);

    virtual void compute(long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);
};

MEDDLY::preplus_evplus::preplus_evplus(forest* arg1, forest* arg2, forest* res)
  : prepostplus_evplus(PREPLUS_cache, arg1, arg2, res)
{
}

void MEDDLY::preplus_evplus::compute(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;

  ct_entry_key* Key = findResult(aev, a, bev, b, cev, c);
  if (0==Key) return;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  MEDDLY_DCASSERT(bLevel >= 0);

  const int resultLevel = ABS(aLevel) > bLevel ? ABS(aLevel): bLevel;
  const unsigned resultSize = unsigned(resF->getLevelSize(resultLevel));

  // Initialize result
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, 0L, a, FULL_ONLY)
    : arg1F->newUnpacked(a, FULL_ONLY)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, 0L, b, FULL_ONLY)
    : arg2F->newUnpacked(b, FULL_ONLY)
  ;

  // do computation
  for (unsigned i = 0; i < resultSize; i++) {
    if (A->down(i) == 0 || B->down(i) == 0) {
      nb->setFull(i, edge_value(0L), 0);
      // nb->d_ref(i) = 0;
      // nb->setEdge(i, 0L);
      continue;
    }

    int dLevel = arg1F->getNodeLevel(A->down(i));
    unpacked_node *D = (dLevel != -resultLevel)
      ? unpacked_node::newIdentity(arg1F, -resultLevel, i, 0L, A->down(i), FULL_ONLY)
      : arg1F->newUnpacked(A->down(i), FULL_ONLY);

    unpacked_node* nb2 = unpacked_node::newFull(resF, -resultLevel, resultSize);

    for (unsigned j = 0; j < resultSize; j++) {
      if (D->down(j) == 0) {
        nb2->setFull(j, edge_value(0L), 0);
        // nb2->d_ref(j) = 0;
        // nb2->setEdge(j, 0L);
        continue;
      }

      long ev = Inf<long>();
      node_handle ed = 0;
      compute(aev + A->edge_long(i) + D->edge_long(j), D->down(j),
        bev + B->edge_long(i), B->down(i),
        ev, ed);
      nb2->setFull(j, edge_value(ev), ed);
      // nb2->d_ref(j) = ed;
      // nb2->setEdge(j, ev);
    }

    unpacked_node::Recycle(D);

    // Reduce
    long dev = Inf<long>();
    node_handle d = 0;
    resF->createReducedNode(int(i), nb2, dev, d);

    nb->setFull(i, edge_value(dev), d);
    // nb->d_ref(i) = d;
    // nb->setEdge(i, dev);
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // Reduce
  node_handle cl;
  resF->createReducedNode(-1, nb, cev, cl);
  c = cl;

  // Add to CT
  saveResult(Key, aev, a, bev, b, cev, c);
}

// ******************************************************************
// *                                                                *
// *                     postplus_evplus  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::postplus_evplus : public prepostplus_evplus {
  public:
    postplus_evplus(forest* arg1, forest* arg2, forest* res);

    virtual void compute(long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);
};

MEDDLY::postplus_evplus
 ::postplus_evplus(forest* arg1, forest* arg2, forest* res)
  : prepostplus_evplus(POSTPLUS_cache, arg1, arg2, res)
{
}

void MEDDLY::postplus_evplus::compute(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;

  ct_entry_key* Key = findResult(aev, a, bev, b, cev, c);
  if (0==Key) return;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  MEDDLY_DCASSERT(bLevel >= 0);

  const int resultLevel = ABS(aLevel) > bLevel ? ABS(aLevel): bLevel;
  const unsigned resultSize = unsigned(resF->getLevelSize(resultLevel));

  // Initialize result
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, 0L, a, FULL_ONLY)
    : arg1F->newUnpacked(a, FULL_ONLY)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, 0L, b, FULL_ONLY)
    : arg2F->newUnpacked(b, FULL_ONLY)
  ;

  // do computation
  for (unsigned i = 0; i < resultSize; i++) {
    if (A->down(i) == 0) {
      nb->setFull(i, edge_value(0L), 0);
      // nb->d_ref(i) = 0;
      // nb->setEdge(i, 0L);
      continue;
    }

    int dLevel = arg1F->getNodeLevel(A->down(i));
    unpacked_node *D = (dLevel != -resultLevel)
      ? unpacked_node::newIdentity(arg1F, -resultLevel, i, 0L, A->down(i), FULL_ONLY)
      : arg1F->newUnpacked(A->down(i), FULL_ONLY);

    unpacked_node* nb2 = unpacked_node::newFull(resF, -resultLevel, resultSize);

    for (unsigned j = 0; j < resultSize; j++) {
      if (D->down(j) == 0 || B->down(j) == 0) {
        nb2->setFull(j, edge_value(0L), 0);
        // nb2->d_ref(j) = 0;
        // nb2->setEdge(j, 0L);
        continue;
      }

      long ev = Inf<long>();
      node_handle ed = 0;
      compute(aev + A->edge_long(i) + D->edge_long(j), D->down(j),
        bev + B->edge_long(j), B->down(j),
        ev, ed);
      nb2->setFull(j, edge_value(ev), ed);
      // nb2->d_ref(j) = ed;
      // nb2->setEdge(j, ev);
    }

    unpacked_node::Recycle(D);

    // Reduce
    long dev = Inf<long>();
    node_handle d = 0;
    resF->createReducedNode(int(i), nb2, dev, d);

    nb->setFull(i, edge_value(dev), d);
    // nb->d_ref(i) = d;
    // nb->setEdge(i, dev);
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // Reduce
  node_handle cl;
  resF->createReducedNode(-1, nb, cev, cl);
  c = cl;

  // Add to CT
  saveResult(Key, aev, a, bev, b, cev, c);
}

// ******************************************************************
// *                                                                *
// *                      preplus_opname  class                     *
// *                                                                *
// ******************************************************************

/*
class MEDDLY::preplus_opname : public binary_opname {
  public:
    preplus_opname();
    virtual binary_operation* buildOperation(binary_list &c, forest* a1,
      forest* a2, forest* r);
};

MEDDLY::preplus_opname::preplus_opname()
 : binary_opname("Pre-Plus")
{
}

MEDDLY::binary_operation*
MEDDLY::preplus_opname::buildOperation(binary_list &c, forest* a1, forest* a2,
  forest* r)
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH);


  // Exception: EV+MXD + EVMDD = EV+MXD
  if (
    a1->getEdgeLabeling() == edge_labeling::EVPLUS ||
    a2->getEdgeLabeling() == edge_labeling::EVPLUS ||
    r->getEdgeLabeling() == edge_labeling::EVPLUS ||
    a1->isForRelations() ||
    !a2->isForRelations() ||
    r->isForRelations()
  ) {
    return new preplus_evplus(c, a1, a2, r);
  }

  throw error(error::NOT_IMPLEMENTED);
}
*/

// ******************************************************************
// *                                                                *
// *                     postplus_opname  class                     *
// *                                                                *
// ******************************************************************

/*
class MEDDLY::postplus_opname : public binary_opname {
  public:
    postplus_opname();
    virtual binary_operation* buildOperation(binary_list &c, forest* a1,
      forest* a2, forest* r);
};

MEDDLY::postplus_opname::postplus_opname()
 : binary_opname("Post-Plus")
{
}

MEDDLY::binary_operation*
MEDDLY::postplus_opname::buildOperation(binary_list &c, forest* a1, forest* a2,
  forest* r)
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH);


  // Exception: EV+MXD + EVMDD = EV+MXD
  if (
    a1->getEdgeLabeling() == edge_labeling::EVPLUS ||
    a2->getEdgeLabeling() == edge_labeling::EVPLUS ||
    r->getEdgeLabeling() == edge_labeling::EVPLUS ||
    a1->isForRelations() ||
    !a2->isForRelations() ||
    r->isForRelations()
  ) {
    return new postplus_evplus(c, a1, a2, r);
  }

  throw error(error::NOT_IMPLEMENTED);
}
*/

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::PRE_PLUS(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;

    binary_operation* bop =  PREPLUS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if  (
            a->getEdgeLabeling() == edge_labeling::EVPLUS ||
            b->getEdgeLabeling() == edge_labeling::EVPLUS ||
            c->getEdgeLabeling() == edge_labeling::EVPLUS ||
            a->isForRelations() ||
            !b->isForRelations() ||
            c->isForRelations()
        )
    {
        return PREPLUS_cache.add(new preplus_evplus(a, b, c));
    }

    throw error(error::NOT_IMPLEMENTED);
}

void MEDDLY::PRE_PLUS_init()
{
    PREPLUS_cache.reset("Pre-plus");
}

void MEDDLY::PRE_PLUS_done()
{
    MEDDLY_DCASSERT(PREPLUS_cache.isEmpty());
}

MEDDLY::binary_operation* MEDDLY::POST_PLUS(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;

    binary_operation* bop =  POSTPLUS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if  (
            a->getEdgeLabeling() == edge_labeling::EVPLUS ||
            b->getEdgeLabeling() == edge_labeling::EVPLUS ||
            c->getEdgeLabeling() == edge_labeling::EVPLUS ||
            a->isForRelations() ||
            !b->isForRelations() ||
            c->isForRelations()
        )
    {
        return POSTPLUS_cache.add(new postplus_evplus(a, b, c));
    }

    throw error(error::NOT_IMPLEMENTED);
}

void MEDDLY::POST_PLUS_init()
{
    POSTPLUS_cache.reset("Post-plus");
}

void MEDDLY::POST_PLUS_done()
{
    MEDDLY_DCASSERT(POSTPLUS_cache.isEmpty());
}


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
#include "prepostplus.h"
#include "apply_base.h"

namespace MEDDLY {
  class prepostplus_evplus;

  class preplus_evplus;
  class postplus_evplus;

  class preplus_opname;
  class postplus_opname;
};

// ******************************************************************
// *                                                                *
// *                    prepostplus_evplus  class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::prepostplus_evplus : public generic_binary_evplus {
  public:
    prepostplus_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(long aev, node_handle a, long bev, node_handle b,
      long& cev, node_handle& c);
};

MEDDLY::prepostplus_evplus::prepostplus_evplus(const binary_opname* opcode,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_evplus(opcode, arg1, arg2, res)
{
  MEDDLY_DCASSERT(arg1->isForRelations());
  MEDDLY_DCASSERT(!arg2->isForRelations());

  operationCommutes();
}

bool MEDDLY::prepostplus_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == -1 && b == -1) {
    c = -1;
    MEDDLY_DCASSERT(aev != Inf<long>());
    MEDDLY_DCASSERT(bev != Inf<long>());
    cev = aev + bev;
    MEDDLY_DCASSERT(cev >= 0);
    return true;
  }
  if (a == 0 || b == 0) {
    c = 0;
    cev = Inf<long>();
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
    preplus_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual void compute(long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);
};

MEDDLY::preplus_evplus::preplus_evplus(const binary_opname* opcode,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : prepostplus_evplus(opcode, arg1, arg2, res)
{
}

void MEDDLY::preplus_evplus::compute(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;

  compute_table::search_key* Key = findResult(aev, a, bev, b, cev, c);
  if (0==Key) return;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  MEDDLY_DCASSERT(bLevel >= 0);

  const int resultLevel = ABS(aLevel) > bLevel ? ABS(aLevel): bLevel;
  const int resultSize = resF->getLevelSize(resultLevel);

  // Initialize result
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, 0L, a, true)
    : unpacked_node::newFromNode(arg1F, a, true)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, 0L, b, true)
    : unpacked_node::newFromNode(arg2F, b, true)
  ;

  // do computation
  for (int i = 0; i < resultSize; i++) {
    if (A->d(i) == 0 || B->d(i) == 0) {
      nb->d_ref(i) = 0;
      nb->setEdge(i, Inf<long>());
      continue;
    }

    int dLevel = arg1F->getNodeLevel(A->d(i));
    unpacked_node *D = (dLevel != -resultLevel)
      ? unpacked_node::newIdentity(arg1F, -resultLevel, i, 0L, A->d(i), true)
      : unpacked_node::newFromNode(arg1F, A->d(i), true);

    unpacked_node* nb2 = unpacked_node::newFull(resF, -resultLevel, resultSize);

    for (int j = 0; j < resultSize; j++) {
      if (D->d(j) == 0) {
        nb2->d_ref(j) = 0;
        nb2->setEdge(j, Inf<long>());
        continue;
      }

      long ev = Inf<long>();
      node_handle ed = 0;
      compute(aev + A->ei(i) + D->ei(j), D->d(j),
        bev + B->ei(i), B->d(i),
        ev, ed);
      nb2->d_ref(j) = ed;
      nb2->setEdge(j, ev);
    }

    unpacked_node::recycle(D);

    // Reduce
    long dev = Inf<long>();
    node_handle d = 0;
    resF->createReducedNode(i, nb2, dev, d);

    nb->d_ref(i) = d;
    nb->setEdge(i, dev);
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

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
    postplus_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual void compute(long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);
};

MEDDLY::postplus_evplus::postplus_evplus(const binary_opname* opcode,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : prepostplus_evplus(opcode, arg1, arg2, res)
{
}

void MEDDLY::postplus_evplus::compute(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;

  compute_table::search_key* Key = findResult(aev, a, bev, b, cev, c);
  if (0==Key) return;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  MEDDLY_DCASSERT(bLevel >= 0);

  const int resultLevel = ABS(aLevel) > bLevel ? ABS(aLevel): bLevel;
  const int resultSize = resF->getLevelSize(resultLevel);

  // Initialize result
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, 0L, a, true)
    : unpacked_node::newFromNode(arg1F, a, true)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, 0L, b, true)
    : unpacked_node::newFromNode(arg2F, b, true)
  ;

  // do computation
  for (int i = 0; i < resultSize; i++) {
    if (A->d(i) == 0) {
      nb->d_ref(i) = 0;
      nb->setEdge(i, Inf<long>());
      continue;
    }

    int dLevel = arg1F->getNodeLevel(A->d(i));
    unpacked_node *D = (dLevel != -resultLevel)
      ? unpacked_node::newIdentity(arg1F, -resultLevel, i, 0L, A->d(i), true)
      : unpacked_node::newFromNode(arg1F, A->d(i), true);

    unpacked_node* nb2 = unpacked_node::newFull(resF, -resultLevel, resultSize);

    for (int j = 0; j < resultSize; j++) {
      if (D->d(j) == 0 || B->d(j) == 0) {
        nb2->d_ref(j) = 0;
        nb2->setEdge(j, Inf<long>());
        continue;
      }

      long ev = Inf<long>();
      node_handle ed = 0;
      compute(aev + A->ei(i) + D->ei(j), D->d(j),
        bev + B->ei(j), B->d(j),
        ev, ed);
      nb2->d_ref(j) = ed;
      nb2->setEdge(j, ev);
    }

    unpacked_node::recycle(D);

    // Reduce
    long dev = Inf<long>();
    node_handle d = 0;
    resF->createReducedNode(i, nb2, dev, d);

    nb->d_ref(i) = d;
    nb->setEdge(i, dev);
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

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

class MEDDLY::preplus_opname : public binary_opname {
  public:
    preplus_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::preplus_opname::preplus_opname()
 : binary_opname("Pre-Plus")
{
}

MEDDLY::binary_operation*
MEDDLY::preplus_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH);


  // Exception: EV+MXD + EVMDD = EV+MXD
  if (
    a1->getEdgeLabeling() == forest::EVPLUS ||
    a2->getEdgeLabeling() == forest::EVPLUS ||
    r->getEdgeLabeling() == forest::EVPLUS ||
    a1->isForRelations() ||
    !a2->isForRelations() ||
    r->isForRelations()
  ) {
    return new preplus_evplus(this, a1, a2, r);
  }

  throw error(error::NOT_IMPLEMENTED);
}

// ******************************************************************
// *                                                                *
// *                     postplus_opname  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::postplus_opname : public binary_opname {
  public:
    postplus_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::postplus_opname::postplus_opname()
 : binary_opname("Post-Plus")
{
}

MEDDLY::binary_operation*
MEDDLY::postplus_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH);


  // Exception: EV+MXD + EVMDD = EV+MXD
  if (
    a1->getEdgeLabeling() == forest::EVPLUS ||
    a2->getEdgeLabeling() == forest::EVPLUS ||
    r->getEdgeLabeling() == forest::EVPLUS ||
    a1->isForRelations() ||
    !a2->isForRelations() ||
    r->isForRelations()
  ) {
    return new postplus_evplus(this, a1, a2, r);
  }

  throw error(error::NOT_IMPLEMENTED);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializePrePlus()
{
  return new preplus_opname();
}

MEDDLY::binary_opname* MEDDLY::initializePostPlus()
{
  return new postplus_opname();
}

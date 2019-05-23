
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
#include "cycle.h"

namespace MEDDLY {
  class cycle_opname;

  class cycle_EV2EV;
};

// ******************************************************************
// *                                                                *
// *                       cycle_EV2MT  class                       *
// *                                                                *
// ******************************************************************

// Extract cycles (EV+MDD) from transitive closure (EV+MxD).
class MEDDLY::cycle_EV2EV : public unary_operation {
  public:
    cycle_EV2EV(const unary_opname* oc, expert_forest* arg, expert_forest* res);

    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res);

  protected:
    virtual void compute_r(long aev, node_handle a, int k, long& bev, node_handle& b);

    inline compute_table::entry_key*
    findResult(long aev, node_handle a, long& bev, node_handle &b)
    {
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      bev = CTresult[0].readL();
      b = resF->linkNode(CTresult[0].readN());
      if (b != 0) {
        bev += aev;
      }
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(compute_table::entry_key* Key,
      long aev, node_handle a, long bev, node_handle b)
    {
      CTresult[0].reset();
      CTresult[0].writeL(b == 0 ? 0L : bev - aev);
      CTresult[0].writeN(b);
      CT0->addEntry(Key, CTresult[0]);
      return b;
    }
};

MEDDLY::cycle_EV2EV::cycle_EV2EV(const unary_opname* oc, expert_forest* arg, expert_forest* res)
  : unary_operation(oc, 1, arg, res)
{
  MEDDLY_DCASSERT(argF->isEVPlus() && argF->isForRelations());
  MEDDLY_DCASSERT(resF->isEVPlus() && !resF->isForRelations());

  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "LN:LN");
  et->setForestForSlot(1, arg);
  et->setForestForSlot(4, res);
  registerEntryType(0, et);
  buildCTs();
}

void MEDDLY::cycle_EV2EV::computeDDEdge(const dd_edge &arg, dd_edge &res)
{
  long aev = Inf<long>();
  arg.getEdgeValue(aev);
  long bev = Inf<long>();
  node_handle b = 0;
  compute_r(aev, arg.getNode(), argF->getNumVariables(), bev, b);
  res.set(b, bev);
}

void MEDDLY::cycle_EV2EV::compute_r(long aev, node_handle a, int k, long& bev, node_handle& b)
{
  if ((!resF->isQuasiReduced() || k == 0) && argF->isTerminalNode(a)) {
    if (a == 0) {
      bev = 0;
      b = 0;
    }
    else {
      bev = aev;
      b = a;
    }
    return;
  }

  if (resF->isQuasiReduced() && ABS(argF->getNodeLevel(a)) < k) {
    int size = resF->getLevelSize(k);
    unpacked_node* T = unpacked_node::newFull(resF, k, size);
    long tev = Inf<long>();
    node_handle t = 0;
    compute_r(aev, a, k - 1, tev, t);
    for (int i = 0; i < size; i++) {
      T->setEdge(i, tev);
      T->d_ref(i) = resF->linkNode(t);
    }
    resF->unlinkNode(t);
    resF->createReducedNode(-1, T, bev, b);
    return;
  }

  // check the cache
  compute_table::entry_key* key = findResult(aev, a, bev, b);
  if (key == 0) {
    return;
  }

  const int aLevel = argF->getNodeLevel(a);
  const int level = ABS(aLevel);
  const int size = resF->getLevelSize(level);

  unpacked_node* A = aLevel < 0
    ? unpacked_node::newRedundant(argF, level, 0L, a, true)
    : unpacked_node::newFromNode(argF, a, true);
  unpacked_node* T = unpacked_node::newFull(resF, level, size);
  for (int i = 0; i < size; i++) {
    unpacked_node* B = isLevelAbove(-level, argF->getNodeLevel(A->d(i)))
      ? (argF->isIdentityReduced()
        ? unpacked_node::newIdentity(argF, -level, i, 0L, A->d(i), true)
        : unpacked_node::newRedundant(argF, -level, 0L, A->d(i), true))
      : unpacked_node::newFromNode(argF, A->d(i), true);

    long tev = Inf<long>();
    node_handle t = 0;
    compute_r(aev + A->ei(i) + B->ei(i), B->d(i), level - 1, tev, t);
    T->setEdge(i, tev);
    T->d_ref(i) = t;

    unpacked_node::recycle(B);
  }

  unpacked_node::recycle(A);

  resF->createReducedNode(-1, T, bev, b);
  saveResult(key, aev, a, bev, b);
}

// ******************************************************************
// *                                                                *
// *                         cycle_opname                           *
// *                                                                *
// ******************************************************************

class MEDDLY::cycle_opname : public unary_opname {
public:
  cycle_opname();
  virtual unary_operation* buildOperation(expert_forest* arg, expert_forest* res) const;
};

MEDDLY::cycle_opname::cycle_opname()
  : unary_opname("Extract cycles")
{
}

MEDDLY::unary_operation* MEDDLY::cycle_opname::buildOperation(
  expert_forest* arg, expert_forest* res) const
{
  unary_operation* op = 0;
  if (arg->isEVPlus() && arg->isForRelations() && res->isEVPlus() && !res->isForRelations()) {
    op = new cycle_EV2EV(this, arg, res);
  }
  else {
    throw error(error::NOT_IMPLEMENTED);
  }
  return op;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeCycle()
{
  return new cycle_opname;
}


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

    virtual bool isStaleEntry(const node_handle* data);
    virtual void discardEntry(const node_handle* data);
    virtual void showEntry(output &strm, const node_handle* data) const;
    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res);

  protected:
    static const int NODE_INDICES_IN_KEY[2];

    virtual void compute_r(long aev, node_handle a, long& bev, node_handle& b);

    inline compute_table::search_key*
    findResult(long aev, node_handle a, long& bev, node_handle &b)
    {
      compute_table::search_key* CTsrch = useCTkey();
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->reset();
      CTsrch->writeNH(a);
      compute_table::search_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      cacheFind.read(bev);
      b = resF->linkNode(cacheFind.readNH());
      if (b != 0) {
        bev += aev;
      }
      doneCTkey(CTsrch);
      return 0;
    }
    inline node_handle saveResult(compute_table::search_key* Key,
      long aev, node_handle a, long bev, node_handle b)
    {
      argF->cacheNode(a);
      compute_table::entry_builder &entry = CT->startNewEntry(Key);
      entry.writeResult(b == 0 ? 0L : bev - aev);
      entry.writeResultNH(resF->cacheNode(b));
      CT->addEntry();
      return b;
    }
};

const int MEDDLY::cycle_EV2EV::NODE_INDICES_IN_KEY[2] = {
  0,
  (sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)
};

MEDDLY::cycle_EV2EV::cycle_EV2EV(const unary_opname* oc, expert_forest* arg, expert_forest* res)
  : unary_operation(oc,
      sizeof(node_handle) / sizeof(node_handle),
      (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle),
      arg, res)
{
  MEDDLY_DCASSERT(argF->isEVPlus() && argF->isForRelations());
  MEDDLY_DCASSERT(resF->isEVPlus() && !resF->isForRelations());
}

bool MEDDLY::cycle_EV2EV::isStaleEntry(const node_handle* data)
{
  return argF->isStale(data[NODE_INDICES_IN_KEY[0]])
    || resF->isStale(data[NODE_INDICES_IN_KEY[1]]);
}

void MEDDLY::cycle_EV2EV::discardEntry(const node_handle* data)
{
  argF->uncacheNode(data[NODE_INDICES_IN_KEY[0]]);
  resF->uncacheNode(data[NODE_INDICES_IN_KEY[1]]);
}

void MEDDLY::cycle_EV2EV::showEntry(output &strm, const node_handle* data) const
{
  strm << "[" << getName()
    << "(" << long(data[NODE_INDICES_IN_KEY[0]])
    << "): " << long(data[NODE_INDICES_IN_KEY[3]])
    << "]";
}

void MEDDLY::cycle_EV2EV::computeDDEdge(const dd_edge &arg, dd_edge &res)
{
  long aev = Inf<long>();
  arg.getEdgeValue(aev);
  long bev = Inf<long>();
  node_handle b = 0;
  compute_r(aev, arg.getNode(), bev, b);
  res.set(b, bev);
}

void MEDDLY::cycle_EV2EV::compute_r(long aev, node_handle a, long& bev, node_handle& b)
{
  if (argF->isTerminalNode(a)) {
    bev = aev;
    b = a;
    return;
  }

  // check the cache
  compute_table::search_key* key = findResult(aev, a, bev, b);
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
    compute_r(aev + A->ei(i) + B->ei(i), B->d(i), tev, t);
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
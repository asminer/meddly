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

#include "../oper_unary.h"
#include "../ct_vector.h"
#include "../forest_levels.h"

/*
#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
*/

namespace MEDDLY {
    class cycle_opname;
    class cycle_EV2EV;
    class CYCLE_factory;
};

// ******************************************************************
// *                                                                *
// *                       cycle_EV2MT  class                       *
// *                                                                *
// ******************************************************************

// Extract cycles (EV+MDD) from transitive closure (EV+MxD).
class MEDDLY::cycle_EV2EV : public unary_operation {
    public:
        cycle_EV2EV(forest* arg, forest* res);
        virtual ~cycle_EV2EV();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp);
    protected:
        void _compute(int L, long aev, node_handle a,
                long &bev, node_handle &b);

    private:
        ct_entry_type* ct;


        // OLD BELOW HERE

/*
        virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);

  protected:
        virtual void compute_r(long aev, node_handle a, int k, long& bev, node_handle& b);

        inline ct_entry_key*
        findResult(long aev, node_handle a, long& bev, node_handle &b)
        {
            ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
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
        inline node_handle saveResult(ct_entry_key* Key,
            long aev, node_handle a, long bev, node_handle b)
        {
            CTresult[0].reset();
            CTresult[0].writeL(b == 0 ? 0L : bev - aev);
            CTresult[0].writeN(b);
            CT0->addEntry(Key, CTresult[0]);
            return b;
        }
        */
};

MEDDLY::cycle_EV2EV::cycle_EV2EV(forest* arg, forest* res)
  : unary_operation(arg, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
    if (!arg->isForRelations() || res->isForRelations()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    ct = new ct_entry_type("cycle"); // CYCLE_cache.getName(), "LN:LN");
    ct->setFixed(arg->getEdgeType(), arg);
    ct->setResult(res->getEdgeType(), res);
    ct->doneBuilding();
/*
    et->setForestForSlot(1, arg);
    et->setForestForSlot(4, res);
    registerEntryType(0, et);
    buildCTs();
    */
}

MEDDLY::cycle_EV2EV::~cycle_EV2EV()
{
    ct->markForDestroy();
}

void MEDDLY::cycle_EV2EV::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp)
{
    long aev, cev;
    av.get(aev);
    _compute(L, aev, ap, cev, cp);
    cv.set(cev);
}

/*

void MEDDLY::cycle_EV2EV::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{
  long aev = Inf<long>();
  arg.getEdgeValue(aev);
  long bev = Inf<long>();
  node_handle b = 0;
  compute_r(aev, arg.getNode(), argF->getMaxLevelIndex(), bev, b);
  res.set(bev, b);
}
*/

void MEDDLY::cycle_EV2EV::_compute(int k, long aev, node_handle a,
                long &bev, node_handle &b)
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
        unpacked_node* T = unpacked_node::newWritable(resF, k, size, FULL_ONLY);
        long tev;
        node_handle t;
        _compute(k-1, aev, a, tev, t);
        for (int i = 0; i < size; i++) {
            T->setFull(i, edge_value(tev), resF->linkNode(t));
            // T->setEdge(i, tev);
            // T->d_ref(i) = resF->linkNode(t);
        }
        resF->unlinkNode(t);
        edge_value ev;
        resF->createReducedNode(T, ev, b);
        bev = ev.getLong();
        return;
    }

    // check the cache
    ct_vector key(2);
    ct_vector res(2);
    key[0].set(aev);
    key[1].setN(a);
    if (ct->findCT(key, res)) {
        bev = res[0].getL();
        b = resF->linkNode(res[1].getN());
        return;
    }

    const int aLevel = argF->getNodeLevel(a);
    const int level = ABS(aLevel);
    const int size = resF->getLevelSize(level);

    unpacked_node* A = aLevel < 0
        ? unpacked_node::newRedundant(argF, level, 0L, a, FULL_ONLY)
        : unpacked_node::newFromNode(argF, a, FULL_ONLY);
    unpacked_node* T = unpacked_node::newWritable(resF, level, size, FULL_ONLY);
    for (int i = 0; i < size; i++) {
        unpacked_node* B = isLevelAbove(-level, argF->getNodeLevel(A->down(i)))
            ? (argF->isIdentityReduced()
                ? unpacked_node::newIdentity(argF, -level, i, 0L, A->down(i), FULL_ONLY)
                : unpacked_node::newRedundant(argF, -level, 0L, A->down(i), FULL_ONLY))
            : unpacked_node::newFromNode(argF, A->down(i), FULL_ONLY);

        long tev;
        node_handle t;
        _compute(level - 1, aev + A->edgeval(i).getLong() + B->edgeval(i).getLong(), B->down(i), tev, t);
        T->setFull(i, edge_value(tev), t);
        // T->setEdge(i, tev);
        // T->d_ref(i) = t;

        unpacked_node::Recycle(B);
    }

    unpacked_node::Recycle(A);

    edge_value ev;
    resF->createReducedNode(T, ev, b);
    bev = ev.getLong();

    // Add to cache
    res[0].set(bev);
    res[1].setN(b);
    ct->addCT(key, res);
}

/*

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
    unpacked_node* T = unpacked_node::newWritable(resF, k, size, FULL_ONLY);
    long tev = Inf<long>();
    node_handle t = 0;
    compute_r(aev, a, k - 1, tev, t);
    for (int i = 0; i < size; i++) {
      T->setFull(i, edge_value(tev), resF->linkNode(t));
      // T->setEdge(i, tev);
      // T->d_ref(i) = resF->linkNode(t);
    }
    resF->unlinkNode(t);
    edge_value ev;
    resF->createReducedNode(T, ev, b);
    bev = ev.getLong();
    return;
  }

  // check the cache
  ct_entry_key* key = findResult(aev, a, bev, b);
  if (key == 0) {
    return;
  }

  const int aLevel = argF->getNodeLevel(a);
  const int level = ABS(aLevel);
  const int size = resF->getLevelSize(level);

  unpacked_node* A = aLevel < 0
    ? unpacked_node::newRedundant(argF, level, 0L, a, FULL_ONLY)
    : unpacked_node::newFromNode(argF, a, FULL_ONLY);
  unpacked_node* T = unpacked_node::newWritable(resF, level, size, FULL_ONLY);
  for (int i = 0; i < size; i++) {
    unpacked_node* B = isLevelAbove(-level, argF->getNodeLevel(A->down(i)))
      ? (argF->isIdentityReduced()
        ? unpacked_node::newIdentity(argF, -level, i, 0L, A->down(i), FULL_ONLY)
        : unpacked_node::newRedundant(argF, -level, 0L, A->down(i), FULL_ONLY))
      : unpacked_node::newFromNode(argF, A->down(i), FULL_ONLY);

    long tev = Inf<long>();
    node_handle t = 0;
    compute_r(aev + A->edgeval(i).getLong() + B->edgeval(i).getLong(), B->down(i), level - 1, tev, t);
    T->setFull(i, edge_value(tev), t);
    // T->setEdge(i, tev);
    // T->d_ref(i) = t;

    unpacked_node::Recycle(B);
  }

  unpacked_node::Recycle(A);

  edge_value ev;
  resF->createReducedNode(T, ev, b);
  bev = ev.getLong();
  saveResult(key, aev, a, bev, b);
}
*/

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

/*
MEDDLY::unary_operation* MEDDLY::CYCLE(forest* arg, forest* res)
{
    if (!arg) return nullptr;
    unary_operation* uop =  CYCLE_cache.find(arg, res);
    if (uop) {
        return uop;
    }
    if (arg->isEVPlus()) {
        return CYCLE_cache.add(new cycle_EV2EV(arg, res));
    }
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::CYCLE_init()
{
    CYCLE_cache.reset("ExtractCycles");
}

void MEDDLY::CYCLE_done()
{
    MEDDLY_DCASSERT(CYCLE_cache.isEmpty());
}
*/

// ******************************************************************
// *                                                                *
// *                      CYCLE_factory  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::CYCLE_factory : public unary_factory {
    public:
        virtual void setup();
        virtual unary_operation* build(forest* arg, forest* res);
};

// ******************************************************************

void MEDDLY::CYCLE_factory::setup()
{
    _setup(__FILE__, "CYCLE", "Cycle detection. The input forest should be an EV+MxD (relation) produced by a transitive closure operation, and the output should be an EV+MDD (set). TBD: needs testing");
}

MEDDLY::unary_operation* MEDDLY::CYCLE_factory::build(forest* arg, forest* res)
{
    if (!arg) return nullptr;
    unary_operation* uop =  cache_find(arg, res);
    if (uop) {
        return uop;
    }
    if (arg->isEVPlus()) {
        return cache_add(new cycle_EV2EV(arg, res));
    }
    return nullptr;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_factory& MEDDLY::CYCLE()
{
    static CYCLE_factory F;
    return F;
}


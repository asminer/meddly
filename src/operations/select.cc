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
#include "select.h"

#ifdef ALLOW_DEPRECATED_0_17_8

#include "../forest.h"
#include "../forest_levels.h"
#include "../oper_unary.h"

namespace MEDDLY {
    class select;
    class select_MT;
    class select_EVPlus;

    class SELECT_factory;
};

// ******************************************************************
// *                                                                *
// *                         select  class                          *
// *                                                                *
// ******************************************************************

/// Abstract base class for selecting one state randomly from a set of states.
class MEDDLY::select : public unary_operation {
    public:
        select(forest* arg, forest* res);
};

MEDDLY::select::select(forest* arg, forest* res)
    : unary_operation(arg, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
}

// ******************************************************************
// *                                                                *
// *                        select_MT class                         *
// *                                                                *
// ******************************************************************

// States are not selected with equal probability.
class MEDDLY::select_MT : public select {
  public:
    select_MT(forest* arg, forest* res);
    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);
    virtual node_handle _compute(node_handle node, int level);
};

MEDDLY::select_MT::select_MT(forest* arg, forest* res)
    : select(arg, res)
{
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

void MEDDLY::select_MT::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{
  node_handle result = _compute(arg.getNode(), resF->getMaxLevelIndex());
  res.set(result);
}

MEDDLY::node_handle MEDDLY::select_MT::_compute(node_handle a, int level)
{
  // Check terminals
  if (argF->isTerminalNode(a) && level == 0) {
    MEDDLY_DCASSERT(a != 0);
    return a;
  }

  // Initialize node reader
  unpacked_node* A = isLevelAbove(level, argF->getNodeLevel(a))
    ? unpacked_node::newRedundant(argF, level, a, SPARSE_ONLY)
    : unpacked_node::newFromNode(argF, a, SPARSE_ONLY);
  MEDDLY_DCASSERT(A->getSize() > 0);

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newWritable(resF, level, 1, SPARSE_ONLY);

  // recurse
  unsigned nz = rand() % A->getSize();
  nb->setSparse(0, A->index(nz), _compute(A->down(nz), level - 1));
  // nb->i_ref(0) = A->index(nz);
  // nb->d_ref(0) = _compute(A->down(nz), level - 1);

  // Cleanup
  unpacked_node::Recycle(A);

  edge_value ev;
  node_handle res;
  resF->createReducedNode(nb, ev, res, 0);
  return res;
}

// ******************************************************************
// *                                                                *
// *                      select_EVPlus class                       *
// *                                                                *
// ******************************************************************

// States are not selected with equal probability.
class MEDDLY::select_EVPlus : public select {
  public:
    select_EVPlus(forest* arg, forest* res);
    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);
    virtual void _compute(long aev, node_handle a, int level, long& bev, node_handle& b);
};

MEDDLY::select_EVPlus::select_EVPlus(forest* arg, forest* res)
    : select(arg, res)
{
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

void MEDDLY::select_EVPlus::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{
  long bev = 0;
  node_handle b = 0;
  long aev = 0;
  arg.getEdgeValue(aev);
  _compute(aev, arg.getNode(), resF->getMaxLevelIndex(), bev, b);
  res.set(bev, b);
}

void MEDDLY::select_EVPlus::_compute(long aev, node_handle a, int level, long& bev, node_handle& b)
{
  // Check terminals
  if (argF->isTerminalNode(a) && level == 0) {
    MEDDLY_DCASSERT(a != 0);
    b = a;
    bev = aev;
    return;
  }

  // Initialize node reader
  unpacked_node* A = isLevelAbove(level, argF->getNodeLevel(a))
    ? unpacked_node::newRedundant(argF, level, a, SPARSE_ONLY)
    : unpacked_node::newFromNode(argF, a, SPARSE_ONLY);
  MEDDLY_DCASSERT(A->getSize() > 0);

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newWritable(resF, level, 1, SPARSE_ONLY);

  // recurse
  // Zero-edge-valued indices
  int* izz = new int[A->getSize()];
  unsigned sz = 0;
  for (unsigned iz = 0; iz < A->getSize(); iz++) {
    if (0 == A->edgeval(iz).getLong()) {
      izz[sz] = iz;
      sz++;
    }
  }
  int nz = izz[rand() % sz];
  delete[] izz;

  long tev = 0;
  node_handle t = 0;
  _compute(aev + A->edgeval(nz).getLong(), A->down(nz), level - 1, tev, t);

  nb->setSparse(0, A->index(nz), tev, t);

  // nb->i_ref(0) = A->index(nz);
  // nb->d_ref(0) = t;
  // nb->setEdge(0, tev);

  // Cleanup
  unpacked_node::Recycle(A);

  edge_value ev;
  resF->createReducedNode(nb, ev, b);
  bev = ev.getLong();
}

// ******************************************************************
// *                                                                *
// *                      SELECT_factory class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::SELECT_factory : public unary_factory {
    public:
        virtual void setup();
        virtual unary_operation* build(forest* arg, forest* res);
};

// ******************************************************************

void MEDDLY::SELECT_factory::setup()
{
    _setup(__FILE__, "SELECT", "Select an element. Deprecated.");
}

MEDDLY::unary_operation* MEDDLY::SELECT_factory::build(forest* arg, forest* res)
{
    if (!arg || !res) return nullptr;
    unary_operation* uop =  cache_find(arg, res);
    if (uop) {
        return uop;
    }
    if (arg->isForRelations()) {
        return nullptr;
    }
    if (arg->isMultiTerminal()) {
        return cache_add(new select_MT(arg, res));
    }
    if (arg->isEVPlus()) {
        return cache_add(new select_EVPlus(arg, res));
    }
    return nullptr;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_factory& MEDDLY::SELECT()
{
    static SELECT_factory F;
    return F;
}


#endif


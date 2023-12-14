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

#include "../forest.h"
#include "../oper_unary.h"

namespace MEDDLY {
  class select;
  class select_MT;
  class select_EVPlus;

  class select_opname;
};

// ******************************************************************
// *                                                                *
// *                         select  class                          *
// *                                                                *
// ******************************************************************

/// Abstract base class for selecting one state randomly from a set of states.
class MEDDLY::select : public unary_operation {
  public:
    select(unary_opname* oc, expert_forest* arg, expert_forest* res);
};

MEDDLY::select
:: select(unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, 0, arg, res)
{
  MEDDLY_DCASSERT(!argF->isForRelations());
  MEDDLY_DCASSERT(!resF->isForRelations());
}

// ******************************************************************
// *                                                                *
// *                        select_MT class                         *
// *                                                                *
// ******************************************************************

// States are not selected with equal probability.
class MEDDLY::select_MT : public select {
  public:
    select_MT(unary_opname* oc, expert_forest* arg, expert_forest* res);
    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);
    virtual node_handle _compute(node_handle node, int level);
};

MEDDLY::select_MT
:: select_MT(unary_opname* oc, expert_forest* arg, expert_forest* res)
 : select(oc, arg, res)
{
  MEDDLY_DCASSERT(argF->isMultiTerminal());
  MEDDLY_DCASSERT(resF->isMultiTerminal());
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
    : argF->newUnpacked(a, SPARSE_ONLY);
  MEDDLY_DCASSERT(A->getSize() > 0);

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newSparse(resF, level, 1);

  // recurse
  unsigned nz = rand() % A->getSize();
  nb->i_ref(0) = A->index(nz);
  nb->d_ref(0) = _compute(A->down(nz), level - 1);

  // Cleanup
  unpacked_node::Recycle(A);

  return resF->createReducedNode(0, nb);
}

// ******************************************************************
// *                                                                *
// *                      select_EVPlus class                       *
// *                                                                *
// ******************************************************************

// States are not selected with equal probability.
class MEDDLY::select_EVPlus : public select {
  public:
    select_EVPlus(unary_opname* oc, expert_forest* arg, expert_forest* res);
    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);
    virtual void _compute(long aev, node_handle a, int level, long& bev, node_handle& b);
};

MEDDLY::select_EVPlus
:: select_EVPlus(unary_opname* oc, expert_forest* arg, expert_forest* res)
 : select(oc, arg, res)
{
  MEDDLY_DCASSERT(argF->isEVPlus());
  MEDDLY_DCASSERT(resF->isEVPlus());
}

void MEDDLY::select_EVPlus::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{
  long bev = 0;
  node_handle b = 0;
  long aev = 0;
  arg.getEdgeValue(aev);
  _compute(aev, arg.getNode(), resF->getMaxLevelIndex(), bev, b);
  res.set(b, bev);
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
    : argF->newUnpacked(a, SPARSE_ONLY);
  MEDDLY_DCASSERT(A->getSize() > 0);

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newSparse(resF, level, 1);

  // recurse
  // Zero-edge-valued indices
  int* izz = new int[A->getSize()];
  unsigned sz = 0;
  for (unsigned iz = 0; iz < A->getSize(); iz++) {
    if (0 == A->edge_long(iz)) {
      izz[sz] = iz;
      sz++;
    }
  }
  int nz = izz[rand() % sz];
  delete[] izz;
  nb->i_ref(0) = A->index(nz);

  long tev = 0;
  node_handle t = 0;
  _compute(aev + A->edge_long(nz), A->down(nz), level - 1, tev, t);
  nb->d_ref(0) = t;
  nb->setEdge(0, tev);

  // Cleanup
  unpacked_node::Recycle(A);

  resF->createReducedNode(-1, nb, bev, b);
}


// ******************************************************************
// *                                                                *
// *                       copy_opname  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::select_opname : public unary_opname {
  public:
    select_opname();
    virtual unary_operation*
      buildOperation(expert_forest* ar, expert_forest* res);
};

MEDDLY::select_opname::select_opname()
 : unary_opname("Select")
{
}

MEDDLY::unary_operation*
MEDDLY::select_opname
::buildOperation(expert_forest* arg, expert_forest* res)
{
  if (0==arg || 0==res) return 0;

  if (arg->getDomain() != res->getDomain())
    throw error(error::DOMAIN_MISMATCH);

  if (arg->isForRelations() || res->isForRelations())
    throw error(error::NOT_IMPLEMENTED);

  if (arg->isMultiTerminal() && res->isMultiTerminal()){
    return new select_MT(this, arg, res);
  }

  if (arg->isEVPlus() && res->isEVPlus()) {
    return new select_EVPlus(this, arg, res);
  }

  //
  // Catch all for any other cases
  //
  throw error(error::NOT_IMPLEMENTED);

}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeSelect()
{
  return new select_opname;
}


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
#include "../forests/mt.h"
#include "apply_base.h"

// #define TRACE_ALL_OPS
// #define DISABLE_CACHE


// ******************************************************************
// *                                                                *
// *                   generic_binary_ev  methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_ev::generic_binary_ev(binary_list& code,
  forest* arg1, forest* arg2, forest* res)
  : binary_operation(code, 1, arg1, arg2, res)
{
}

MEDDLY::generic_binary_ev::~generic_binary_ev()
{
}


// ******************************************************************
// *                                                                *
// *                 generic_binary_evplus  methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_evplus::generic_binary_evplus(binary_list& code,
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_ev(code, arg1, arg2, res)
{
  ct_entry_type* et = new ct_entry_type(code.getName(), "LNLN:LN");
  et->setForestForSlot(1, arg1);
  et->setForestForSlot(3, arg2);
  et->setForestForSlot(6, res);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::generic_binary_evplus::~generic_binary_evplus()
{
}

void MEDDLY::generic_binary_evplus
::computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge& c, bool userFlag)
{
  node_handle result;
  long ev = Inf<long>(), aev = Inf<long>(), bev = Inf<long>();
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(aev, a.getNode(), bev, b.getNode(), ev, result);
  c.set(ev, result);
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true, __FILE__, __LINE__);
#endif
}

void MEDDLY::generic_binary_evplus
::compute(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;

  ct_entry_key* Key = findResult(aev, a, bev, b, cev, c);
  if (0==Key) return;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  const int resultLevel = MAX(aLevel, bLevel);
  const unsigned resultSize = unsigned(resF->getLevelSize(resultLevel));

  // Initialize result
  unpacked_node* nb =
    unpacked_node::newWritable(resF, resultLevel, resultSize, FULL_ONLY);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, 0L, a, FULL_ONLY)
    : unpacked_node::newFromNode(arg1F, a, FULL_ONLY)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, 0L, b, FULL_ONLY)
    : unpacked_node::newFromNode(arg2F, b, FULL_ONLY)
  ;

  ASSERT(__FILE__, __LINE__, !A->isExtensible() && !B->isExtensible());


  // do computation
  for (unsigned i=0; i<resultSize; i++) {
    long ev = 0;
    node_handle ed = 0;
    compute(aev + A->edgeval(i).getLong(), A->down(i),
            bev + B->edgeval(i).getLong(), B->down(i),
            ev, ed);
    ASSERT(__FILE__, __LINE__, ed != 0 || ev == 0);
    nb->setFull(i, edge_value(ev), ed);
    // nb->d_ref(i) = ed;
    // nb->setEdge(i, ev);
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // Reduce
  edge_value ev;
  resF->createReducedNode(nb, ev, c);
  cev = ev.getLong();

  // Add to CT
  saveResult(Key, aev, a, bev, b, cev, c);
}


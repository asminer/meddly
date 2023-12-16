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
#include "mdd2index.h"

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"

// #define TRACE_ALL_OPS

namespace MEDDLY {
  class mdd2index_operation;
  class mdd2index_opname;
};

// ******************************************************************
// *                                                                *
// *                   mdd2index_operation  class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::mdd2index_operation : public unary_operation {
  public:
    mdd2index_operation(unary_opname* oc, forest* arg,
      forest* res);

    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);

    void compute_r(int k, node_handle a, node_handle &bdn, long &bcard);
};

MEDDLY::mdd2index_operation::mdd2index_operation(unary_opname* oc,
  forest* arg, forest* res)
  : unary_operation(oc, 1, arg, res)
{
  // answer[0] : node
  // answer[1] : cardinality
  ct_entry_type* et = new ct_entry_type(oc->getName(), "N:NL");
  et->setForestForSlot(0, arg);
  et->setForestForSlot(2, res);
  registerEntryType(0, et);
  buildCTs();
}

void
MEDDLY::mdd2index_operation
::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{
  MEDDLY_DCASSERT(arg.getForest() == argF);
  MEDDLY_DCASSERT(res.getForest() == resF);
  MEDDLY_DCASSERT(resF->handleForValue(false) == 0);
  MEDDLY_DCASSERT(argF->handleForValue(true) < 0);
  MEDDLY_DCASSERT(argF->handleForValue(false) == 0);
  node_handle down;
  long card = 0;
  int nVars = argF->getMaxLevelIndex();
  compute_r(nVars, arg.getNode(), down, card);
  res.set(down, long(0));
}

void
MEDDLY::mdd2index_operation
::compute_r(int k, node_handle a, node_handle &bdn, long &bcard)
{
  // Deal with terminals
  if (0 == a) {
    bdn = 0;
    bcard = 0;
    return;
  }
  if (0 == k) {
    terminal t(true);
    bdn = t.getHandle();
    bcard = 1;
    return;
  }

  int aLevel = argF->getNodeLevel(a);
  MEDDLY_DCASSERT(aLevel <= k);

  // Check compute table
  ct_entry_key* CTsrch = 0;
  if (aLevel == k) {
    CTsrch = CT0->useEntryKey(etype[0], 0);
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeN(a);
    CT0->find(CTsrch, CTresult[0]);
    if (CTresult[0]) {
      bdn = resF->linkNode(CTresult[0].readN());
      bcard = CTresult[0].readL();
      CT0->recycle(CTsrch);
      return;
    }
  }

#ifdef TRACE_ALL_OPS
  printf("calling mdd2index::compute(%d, %d)\n", k, a);
#endif

  // Initialize node builder
  const unsigned size = unsigned(resF->getLevelSize(k));
  unpacked_node* nb = unpacked_node::newFull(resF, k, size);

  // Initialize node reader
  unpacked_node *A = unpacked_node::New(argF);
  if (aLevel < k) {
    A->initRedundant(argF, k, a, FULL_ONLY);
  } else {
    argF->unpackNode(A, a, FULL_ONLY);
  }

  // recurse
  bcard = 0;
  for (unsigned i=0; i<size; i++) {
    node_handle ddn;
    long dcard = 0;
    compute_r(k-1, A->down(i), ddn, dcard);
    // nb->d_ref(i) = ddn;
    if (ddn) {
      nb->setFull(i, edge_value(bcard), ddn);
      // nb->setEdge(i, bcard);
      bcard += dcard;
    } else {
      MEDDLY_DCASSERT(0 == dcard);
      nb->setFull(i, edge_value(0L), ddn);
      // nb->setEdge(i, 0L);
    }
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Reduce
  nb->setUHdata(&bcard);
  // memcpy(nb->UHdata(), &bcard, sizeof(bcard));

  long dummy = 0;
  node_handle bl;
  resF->createReducedNode(-1, nb, dummy, bl);
  bdn = bl;
  MEDDLY_DCASSERT(0 == dummy);

  // Add to compute table
  if (CTsrch) {
    CTresult[0].reset();
    CTresult[0].writeN(bdn);
    CTresult[0].writeL(bcard);
    CT0->addEntry(CTsrch, CTresult[0]);
  }
}

// ******************************************************************
// *                                                                *
// *                     mdd2index_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::mdd2index_opname : public unary_opname {
  public:
    mdd2index_opname();
    virtual unary_operation*
      buildOperation(forest* ar, forest* res);
};

MEDDLY::mdd2index_opname::mdd2index_opname()
 : unary_opname("ConvertToIndexSet")
{
}

MEDDLY::unary_operation*
MEDDLY::mdd2index_opname
::buildOperation(forest* arg, forest* res)
{
  if (0==arg || 0==res) return 0;

  if (arg->getDomain() != res->getDomain())
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (arg->isForRelations() ||
      arg->getRangeType() != range_type::BOOLEAN ||
      arg->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL ||
      res->isForRelations() ||
      res->getRangeType() != range_type::INTEGER ||
      res->getEdgeLabeling() != edge_labeling::INDEX_SET
  ) throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  return new mdd2index_operation(this, arg, res);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeMDD2INDEX()
{
  return new mdd2index_opname;
}


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
#include "mm_mult.h"

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"

// #define TRACE_ALL_OPS

namespace MEDDLY {
  class mm_mult_op;

  class mm_mult_mxd;

  class mm_mult_opname;
};

// ************************************************************************
// *                                                                      *
// *                                                                      *
// *                                                                      *
// *                          actual  operations                          *
// *                                                                      *
// *                                                                      *
// *                                                                      *
// ************************************************************************

// ******************************************************************
// *                                                                *
// *                         mm_mult_op class                       *
// *                                                                *
// ******************************************************************

/// Abstract base class for all matrix-matrix multiplication operations.
class MEDDLY::mm_mult_op : public binary_operation {
  public:
    mm_mult_op(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res, binary_operation* acc);

    inline ct_entry_key*
    findResult(node_handle a, node_handle b, node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CTsrch->writeN(b);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(ct_entry_key* Key,
      node_handle a, node_handle b, node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
      return c;
    }
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);
    virtual node_handle compute(node_handle a, node_handle b);
  protected:
    binary_operation* accumulateOp;
    virtual node_handle compute_rec(node_handle a, node_handle b) = 0;
};

MEDDLY::mm_mult_op::mm_mult_op(binary_list& oc, forest* a1,
  forest* a2, forest* res, binary_operation* acc)
: binary_operation(oc, 1, a1, a2, res)
{
  accumulateOp = acc;

  if (!a1->isForRelations() || !a2->isForRelations())
    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);

  ct_entry_type* et = new ct_entry_type(oc.getName(), "NN:N");
  et->setForestForSlot(0, a1);
  et->setForestForSlot(1, a2);
  et->setForestForSlot(3, res);
  registerEntryType(0, et);
  buildCTs();
}

void MEDDLY::mm_mult_op
::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag)
{
  node_handle cnode;
  cnode = compute(a.getNode(), b.getNode());
  c.set(cnode);
}

MEDDLY::node_handle MEDDLY::mm_mult_op::compute(node_handle a, node_handle b)
{
  MEDDLY_DCASSERT(accumulateOp);
  return compute_rec(a, b);
}


// ******************************************************************
// *                                                                *
// *                         mm_mult_mxd  class                     *
// *                                                                *
// ******************************************************************

/** Generic base for matrix multiplication between two relations.
    Changing what happens at the terminals can give
    different meanings to this operation :^)
*/
class MEDDLY::mm_mult_mxd: public mm_mult_op {
  public:
    mm_mult_mxd(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res, binary_operation* acc);

  protected:
    virtual node_handle compute_rec(node_handle a, node_handle b);
    virtual node_handle processTerminals(node_handle a, node_handle b) = 0;
};

MEDDLY::mm_mult_mxd::mm_mult_mxd(binary_list& oc,
  forest* a1, forest* a2, forest* res,
  binary_operation* acc)
: mm_mult_op(oc, a1, a2, res, acc)
{
}

// TODO Test
MEDDLY::node_handle MEDDLY::mm_mult_mxd::compute_rec(node_handle a,
  node_handle b)
{
  // termination conditions
  if (a == 0 || b == 0) return 0;
  if (arg2F->isTerminalNode(b) && arg1F->isTerminalNode(a)) {
      return processTerminals(a, b);
  }

  // check the cache
  node_handle result = 0;
  ct_entry_key* Key = findResult(a, b, result);
  if (0==Key) return result;

  /**
   * Note: only one node builder can be used at a time for each level.
   *
   * Create a node builder for r. Call it nb.
   * Create a node reader for b.
   * Create a node reader for b[i] for all i in (0, rSize-1) s.t b[i] != 0.
   * Create a node reader for a.
   * Create a UseReader for a[i].
   *
   * For all i in rSize,
   *    If a[i] == 0, set nb[i] = 0 and continue.
   *    Use the UseReader for a[i].
   *    Create a node builder for r[i]. Call it nbi.
   *    For all j in (0, rSize-1),
   *      For all k in (0, rSize-1),
   *        nbi[k] += compute_rec( a[i][j], b[j][k] )
   *    Set nb[i] = reduce(i, nbi).
   * Recycle all node readers (a, b, a[i], b[i]).
   * Result = reduce(-1, nb)
   * Save and return Result.
   */

  // check if b and a are at the same level
  int aLevel = arg1F->getNodeLevel(a);
  int bLevel = arg2F->getNodeLevel(b);

  // No primed levels at this point
  MEDDLY_DCASSERT(aLevel >= 0);
  MEDDLY_DCASSERT(bLevel >= 0);

  // Create a node builder for the result.
  int rLevel = MAX(aLevel, bLevel);
  unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* nbr = unpacked_node::newFull(resF, rLevel, rSize);

  dd_edge resultik(resF), temp(resF);

  // Clear out result (important!)
  // nbr->clear(0, rSize); // newFull does this already

  /**
   * If a is identity reduced (i.e. lower level than b)
   * For all i and j, r[i][j]=compute_rec(a, b[i][j])
   *
   * If b is identity reduced (i.e. lower level than a)
   * For all i and j, r[i][j]=compute_rec(a[i][j], b)
   */

  if (aLevel < bLevel) {
    // For all i and j, r[i][j] = compute_rec(a, b[i][j])
    unpacked_node* nrb = arg2F->newUnpacked(b, SPARSE_ONLY);
    unpacked_node* nrbp = unpacked_node::New(arg2F);
    for (unsigned iz = 0; iz < nrb->getSize(); ++iz) {
      unpacked_node* nbri = unpacked_node::newFull(resF, -rLevel, rSize);
      // for (unsigned i = 0; i < rSize; ++i) nbri->d_ref(i) = 0;
      arg2F->unpackNode(nrbp, nrb->down(iz), SPARSE_ONLY);
      unsigned i = nrb->index(iz);
      for (unsigned jz = 0; jz < nrbp->getSize(); ++jz) {
        unsigned j = nrbp->index(jz);
        MEDDLY_DCASSERT(0 == nbri->down(j));
        // nbri->d_ref(j) = compute_rec(a, nrbp->down(jz));
        nbri->setFull(j, compute_rec(a, nrbp->down(jz)));
      }
      MEDDLY_DCASSERT(0 == nbr->down(i));
      // nbr->d_ref(i) = resF->createReducedNode(int(i), nbri);
      nbr->setFull(i, resF->createReducedNode(int(i), nbri));
    }
    unpacked_node::Recycle(nrbp);
    unpacked_node::Recycle(nrb);
  }
  else if (aLevel > bLevel) {
    // For all i and j, r[i][j]=compute_rec(a[i][j], b)
    unpacked_node* nra = arg1F->newUnpacked(a, SPARSE_ONLY);
    unpacked_node* nrap = unpacked_node::New(arg1F);
    for (unsigned iz = 0; iz < nra->getSize(); ++iz) {
      unpacked_node* nbri = unpacked_node::newFull(resF, -rLevel, rSize);
      // for (unsigned i = 0; i < rSize; ++i) nbri->d_ref(i) = 0;
      arg1F->unpackNode(nrap, nra->down(iz), SPARSE_ONLY);
      unsigned i = nra->index(iz);
      for (unsigned jz = 0; jz < nrap->getSize(); ++jz) {
        unsigned j = nrap->index(jz);
        MEDDLY_DCASSERT(0 == nbri->down(j));
        // nbri->d_ref(j) = compute_rec(nrap->down(jz), b);
        nbri->setFull(j, compute_rec(nrap->down(jz), b));
      }
      MEDDLY_DCASSERT(0 == nbr->down(i));
      // nbr->d_ref(i) = resF->createReducedNode(int(i), nbri);
      nbr->setFull(i, resF->createReducedNode(int(i), nbri));
    }
    unpacked_node::Recycle(nrap);
    unpacked_node::Recycle(nra);
  }
  else {
    // For all i, j and k, r[i][k] += compute_rec(a[i][j], b[j][k])
    MEDDLY_DCASSERT(aLevel == rLevel);
    MEDDLY_DCASSERT(bLevel == rLevel);

    // Node readers for a, b and all b[j].
    unpacked_node* nra = arg1F->newUnpacked(a, SPARSE_ONLY);
    unpacked_node* nrb = arg1F->newUnpacked(b, FULL_ONLY);

    unpacked_node* nrap = unpacked_node::New(arg1F);

    unpacked_node** nrbp = new unpacked_node*[nrb->getSize()];
    for (unsigned i = 0; i < nrb->getSize(); ++i) {
      if (0==nrb->down(i)) {
        nrbp[i] = nullptr;
        continue;
      }
      nrbp[i] = arg2F->newUnpacked(nrb->down(i), SPARSE_ONLY);
    }

    // For all i, j, and k:
    //    result[i][k] += compute_rec(a[i][j], b[j][k])
    for (unsigned iz = 0; iz < nra->getSize(); ++iz) {
      unpacked_node* nbri = unpacked_node::newFull(resF, -rLevel, rSize);
      // for (unsigned i = 0; i < rSize; ++i) nbri->d_ref(i) = 0;
      arg1F->unpackNode(nrap, nra->down(iz), SPARSE_ONLY);
      unsigned i = nra->index(iz);
      for (unsigned jz = 0; jz < nrap->getSize(); ++jz) {
        unsigned j = nrap->index(jz);
        if (0 == nrbp[j]) continue;
        for (unsigned kz = 0; kz < nrbp[j]->getSize(); ++kz) {
          node_handle res = compute_rec(nrap->down(jz), nrbp[j]->down(kz));
          if (0 == res) continue;
          unsigned k = nrbp[j]->index(kz);
          if (0 == nbri->down(k)) {
            nbri->setFull(k, res);
            // nbri->d_ref(k) = res;
            continue;
          }
          // Do the union
          resultik.set(nbri->down(k));
          temp.set(res);
          accumulateOp->computeTemp(resultik, temp, resultik);
          nbri->setFull(k, resultik);
        }
      }
      MEDDLY_DCASSERT(0 == nbr->down(i));
      // nbr->d_ref(i) = resF->createReducedNode(int(i), nbri);
      nbr->setFull(i, resF->createReducedNode(int(i), nbri));
    }
    unpacked_node::Recycle(nrap);
    unpacked_node::Recycle(nra);
    for (unsigned i = 0; i < nrb->getSize(); ++i) {
      if (nrbp[i]) unpacked_node::Recycle(nrbp[i]);
    }
    delete[] nrbp;
    unpacked_node::Recycle(nrb);
  }

  result = resF->createReducedNode(-1, nbr);
#ifdef TRACE_ALL_OPS
  printf("computed new mm_mult_mxd(%d, %d) = %d\n", a, b, result);
#endif
  return saveResult(Key, a, b, result);
}


// ******************************************************************
// *                                                                *
// *                        mm_mult_mt  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

  /** Matrix-Matrix multiplication.
      Matrices are stored as MTMXDs.
  */
  template <typename RTYPE>
  class mm_mult_mt : public mm_mult_mxd {
    public:
      mm_mult_mt(binary_list& opcode, forest* arg1,
        forest* arg2, forest* res, binary_operation* acc)
        : mm_mult_mxd(opcode, arg1, arg2, res, acc) { }

    protected:
      virtual node_handle processTerminals(node_handle a, node_handle b)
      {
        RTYPE aval;
        RTYPE bval;
        RTYPE rval;
        arg1F->getValueFromHandle(a, aval);
        arg2F->getValueFromHandle(b, bval);
        rval = aval * bval;
        return resF->handleForValue(rval);
      }

  };
};


// ************************************************************************
// *                                                                      *
// *                                                                      *
// *                                                                      *
// *                           operation  names                           *
// *                                                                      *
// *                                                                      *
// *                                                                      *
// ************************************************************************


// ******************************************************************
// *                                                                *
// *                      mm_mult_opname  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::mm_mult_opname : public binary_opname {
  public:
    mm_mult_opname();
    virtual binary_operation* buildOperation(binary_list &c, forest* a1,
      forest* a2, forest* r);
};

MEDDLY::mm_mult_opname::mm_mult_opname()
 : binary_opname("Matrix-Matrix multiplication")
{
}

MEDDLY::binary_operation*
MEDDLY::mm_mult_opname::buildOperation(binary_list &c, forest* a1, forest* a2,
  forest* r)
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    (a1->getRangeType() == range_type::BOOLEAN) ||
    (a2->getRangeType() == range_type::BOOLEAN) ||
    (r->getRangeType()  == range_type::BOOLEAN) ||
    !a1->isForRelations()   ||
    !a2->isForRelations()   ||
    !r->isForRelations()    ||
    (a1->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL) ||
    (r->getEdgeLabeling()  != edge_labeling::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

#ifdef USE_NEW_APPLY
  binary_operation* acc = PLUS(r, r, r);
#else
  binary_opname* accop = PLUS();
  MEDDLY_DCASSERT(accop);
  dd_edge er(r);
  binary_operation* acc = accop->getOperation(er, er, er);
#endif

  switch (r->getRangeType()) {
    case range_type::INTEGER:
      return new mm_mult_mt<int>(c, a1, a2, r, acc);

    case range_type::REAL:
      return new mm_mult_mt<float>(c, a1, a2, r, acc);

    default:
      throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
  }
}



// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeMMMultiply()
{
  return new mm_mult_opname;
}



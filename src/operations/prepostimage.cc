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
#include "prepostimage.h"

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"
#include <list>

// #define TRACE_ALL_OPS

namespace MEDDLY {
  class image_op;
  class relXset_mdd;
  class setXrel_mdd;

  class image_op_evplus;
  class relXset_evplus;
  class setXrel_evplus;
  class tcXrel_evplus;
  class covtcmxd;
  class covtc;

  class image_op_mxd;
  class mrrc;

  class preimage_opname;
  class postimage_opname;

  class transitive_closure_postimage_opname;
  class covtc_opname;
  class maintain_reachabilityrelation_cover_opname;

  class VMmult_opname;
  class MVmult_opname;
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
// *                         image_op class                         *
// *                                                                *
// ******************************************************************

/// Abstract base class for all MT-based pre/post image operations.
class MEDDLY::image_op : public binary_operation {
  public:
    image_op(binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

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
    virtual void computeDDEdgeSC(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag,std::list<int>* shouldConfirm);
    virtual node_handle compute(node_handle a, node_handle b);
    virtual node_handle computeSC(node_handle a, node_handle b,std::list<int>* shouldConfirm);
  protected:
    binary_operation* accumulateOp;
    virtual node_handle compute_rec(node_handle a, node_handle b) = 0;
    virtual node_handle compute_recSC(node_handle a, node_handle b,std::list<int>* shouldConfirm) = 0;

    expert_forest* argV;
    expert_forest* argM;
};

MEDDLY::image_op::image_op(binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res, binary_operation* acc)
: binary_operation(oc, 1, a1, a2, res)
{
  accumulateOp = acc;

  if (a1->isForRelations()) {
    argM = a1;
    argV = a2;
    if (a2->isForRelations()) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  } else {
    argM = a2;
    argV = a1;
    if (!a2->isForRelations()) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  }

  ct_entry_type* et = new ct_entry_type(oc->getName(), "NN:N");
  et->setForestForSlot(0, argV);
  et->setForestForSlot(1, argM);
  et->setForestForSlot(3, res);
  registerEntryType(0, et);
  buildCTs();
}

void MEDDLY::image_op
::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag)
{
  node_handle cnode;
  if (a.getForest() == argV) {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", a.getNode(), b.getNode());
#endif
    cnode = compute(a.getNode(), b.getNode());
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", a.getNode(), b.getNode(), cnode);
#endif
  } else {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", b.getNode(), a.getNode());
#endif
    cnode = compute(b.getNode(), a.getNode());
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", b.getNode(), a.getNode(), cnode);
#endif
  }
  c.set(cnode);
}
void MEDDLY::image_op
::computeDDEdgeSC(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag,std::list<int>* shouldConfirm)
{
  node_handle cnode;
  if (a.getForest() == argV) {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", a.getNode(), b.getNode());
#endif
    // printf("computeDDEdgeSC calling computeSC\n" );
    cnode = computeSC(a.getNode(), b.getNode(),shouldConfirm);
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", a.getNode(), b.getNode(), cnode);
#endif
  } else {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", b.getNode(), a.getNode());
#endif
    cnode = compute(b.getNode(), a.getNode());
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", b.getNode(), a.getNode(), cnode);
#endif
  }
  c.set(cnode);
}

MEDDLY::node_handle MEDDLY::image_op::compute(node_handle a, node_handle b)
{
  MEDDLY_DCASSERT(accumulateOp);
  return compute_rec(a, b);
}

MEDDLY::node_handle MEDDLY::image_op::computeSC(node_handle a, node_handle b,std::list<int>* shouldConfirm)
{
  MEDDLY_DCASSERT(accumulateOp);
  return compute_recSC(a, b,shouldConfirm);
}

// ******************************************************************
// *                                                                *
// *                       relXset_mdd  class                       *
// *                                                                *
// ******************************************************************

/** Generic base for relation multiplied by set.
    Changing what happens at the terminals can give
    different meanings to this operation :^)
*/
class MEDDLY::relXset_mdd : public image_op {
  public:
    relXset_mdd(binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

  protected:
    virtual node_handle compute_rec(node_handle a, node_handle b);
    virtual node_handle compute_recSC(node_handle a, node_handle b,std::list<int>* shouldConfirm);
    virtual node_handle processTerminals(node_handle mdd, node_handle mxd) = 0;
};

MEDDLY::relXset_mdd::relXset_mdd(binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res, binary_operation* acc)
: image_op(oc, a1, a2, res, acc)
{
}

MEDDLY::node_handle MEDDLY::relXset_mdd::compute_recSC(node_handle a, node_handle b,std::list<int>* shouldConfirm){printf("return 0\n" );return 0;}

MEDDLY::node_handle MEDDLY::relXset_mdd::compute_rec(node_handle mdd, node_handle mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (argM->isTerminalNode(mxd)) {
    if (argV->isTerminalNode(mdd)) {
      return processTerminals(mdd, mxd);
    }
    // mxd is identity
    if (argV == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  node_handle result = 0;
  ct_entry_key* Key = findResult(mdd, mxd, result);
  if (0==Key) return result;

  // check if mxd and mdd are at the same level
  const int mddLevel = argV->getNodeLevel(mdd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New();
  if (mddLevel < rLevel) {
    A->initRedundant(argV, rLevel, mdd, true);
  } else {
    argV->unpackNode(A, mdd, FULL_ONLY);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      C->d_ref(i) = compute_rec(A->d(i), mxd);
    }
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (unsigned i=0; i<rSize; i++) C->d_ref(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New();
    unpacked_node *Rp = unpacked_node::New();
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, false);
    } else {
      argM->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    dd_edge newstatesE(resF), cdi(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getNNZs(); iz++) {
      unsigned i = Ru->i(iz);
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->d(iz), false);
      } else {
        argM->unpackNode(Rp, Ru->d(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getNNZs(); jz++) {
        unsigned j = Rp->i(jz);
        MEDDLY_DCASSERT(0<=j && j < A->getSize());
        if (0==A->d(j))   continue;
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = compute_rec(A->d(j), Rp->d(jz));
        if (0==newstates) continue;
        if (0==C->d(i)) {
          C->d_ref(i) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        newstatesE.set(newstates);
        cdi.set(C->d(i));
        accumulateOp->computeTemp(newstatesE, cdi, cdi);
        C->set_d(i, cdi);
      } // for j

    } // for i

    unpacked_node::recycle(Rp);
    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);

  result = resF->createReducedNode(-1, C);
#ifdef TRACE_ALL_OPS
  printf("computed relXset(%d, %d) = %d\n", mdd, mxd, result);
#endif
  return saveResult(Key, mdd, mxd, result);
}


// ******************************************************************
// *                                                                *
// *                       setXrel_mdd  class                       *
// *                                                                *
// ******************************************************************

/** Generic base for set multiplied by relation.
    Changing what happens at the terminals can give
    different meanings to this operation :^)
*/
class MEDDLY::setXrel_mdd : public image_op {
  public:
    setXrel_mdd(binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

  protected:
    virtual node_handle compute_rec(node_handle a, node_handle b);
    virtual node_handle compute_recSC(node_handle a, node_handle b,std::list<int>* shouldConfirm);
    virtual node_handle processTerminals(node_handle mdd, node_handle mxd) = 0;
};

MEDDLY::setXrel_mdd::setXrel_mdd(binary_opname* oc,
  expert_forest* a1, expert_forest* a2, expert_forest* res, binary_operation* acc)
: image_op(oc, a1, a2, res, acc)
{
}

MEDDLY::node_handle MEDDLY::setXrel_mdd::compute_recSC(node_handle mdd, node_handle mxd,std::list<int>* shouldConfirm)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (argM->isTerminalNode(mxd)) {
    if (argV->isTerminalNode(mdd)) {
      return processTerminals(mdd, mxd);
    }
    // mxd is identity
    if (argV == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  node_handle result = 0;
  ct_entry_key* Key = findResult(mdd, mxd, result);
  if (0==Key) {
#ifdef TRACE_ALL_OPS
    printf("computing new setXrel(%d, %d), got %d from cache\n", mdd, mxd, result);
#endif

    return result;
  }

#ifdef TRACE_ALL_OPS
  printf("computing new setXrel(%d, %d)\n", mdd, mxd);
#endif


  // check if mxd and mdd are at the same level
  const int mddLevel = argV->getNodeLevel(mdd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New();
  if (mddLevel < rLevel) {
    A->initRedundant(argV, rLevel, mdd, true);
  } else {
    argV->unpackNode(A, mdd, FULL_ONLY);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      C->d_ref(i) = compute_recSC(A->d(i), mxd,shouldConfirm);
    }
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (unsigned i=0; i<rSize; i++) C->d_ref(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New();
    unpacked_node *Rp = unpacked_node::New();
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, false);
    } else {
      argM->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    dd_edge newstatesE(resF), cdj(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getNNZs(); iz++) {
      unsigned i = Ru->i(iz);
      if (0==A->d(i))   continue;
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->d(iz), false);
      } else {
        argM->unpackNode(Rp, Ru->d(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getNNZs(); jz++) {
        unsigned j = Rp->i(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = compute_recSC(A->d(i), Rp->d(jz),shouldConfirm);
        if (0==newstates) continue;
        shouldConfirm[mddLevel].push_back(j);
        if (0==C->d(j)) {
          C->d_ref(j) = newstates;
          continue;
        }

        // there's new states and existing states; union them.
        newstatesE.set(newstates);
        cdj.set(C->d(j));
        accumulateOp->computeTemp(newstatesE, cdj, cdj);
        C->set_d(j, cdj);
      } // for j

    } // for i

    unpacked_node::recycle(Rp);
    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);

  result = resF->createReducedNode(-1, C);
#ifdef TRACE_ALL_OPS
  printf("computed new setXrel(%d, %d) = %d\n", mdd, mxd, result);
#endif
  return saveResult(Key, mdd, mxd, result);
}

MEDDLY::node_handle MEDDLY::setXrel_mdd::compute_rec(node_handle mdd, node_handle mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (argM->isTerminalNode(mxd)) {
    if (argV->isTerminalNode(mdd)) {
      return processTerminals(mdd, mxd);
    }
    // mxd is identity
    if (argV == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  node_handle result = 0;
  ct_entry_key* Key = findResult(mdd, mxd, result);
  if (0==Key) {
#ifdef TRACE_ALL_OPS
    printf("computing new setXrel(%d, %d), got %d from cache\n", mdd, mxd, result);
#endif

    return result;
  }

#ifdef TRACE_ALL_OPS
  printf("computing new setXrel(%d, %d)\n", mdd, mxd);
#endif


  // check if mxd and mdd are at the same level
  const int mddLevel = argV->getNodeLevel(mdd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New();
  if (mddLevel < rLevel) {
    A->initRedundant(argV, rLevel, mdd, true);
  } else {
    argV->unpackNode(A, mdd, FULL_ONLY);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      C->d_ref(i) = compute_rec(A->d(i), mxd);
    }
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (unsigned i=0; i<rSize; i++) C->d_ref(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New();
    unpacked_node *Rp = unpacked_node::New();
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, false);
    } else {
      argM->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    dd_edge newstatesE(resF), cdj(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getNNZs(); iz++) {
      unsigned i = Ru->i(iz);
      if (0==A->d(i))   continue;
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->d(iz), false);
      } else {
        argM->unpackNode(Rp, Ru->d(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getNNZs(); jz++) {
        unsigned j = Rp->i(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = compute_rec(A->d(i), Rp->d(jz));
        if (0==newstates) continue;
        if (0==C->d(j)) {
          C->d_ref(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        newstatesE.set(newstates);
        cdj.set(C->d(j));
        accumulateOp->computeTemp(newstatesE, cdj, cdj);
        C->set_d(j, cdj);
      } // for j

    } // for i

    unpacked_node::recycle(Rp);
    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);

  result = resF->createReducedNode(-1, C);
#ifdef TRACE_ALL_OPS
  printf("computed new setXrel(%d, %d) = %d\n", mdd, mxd, result);
#endif
  return saveResult(Key, mdd, mxd, result);
}

// ******************************************************************
// *                                                                *
// *                      mtvect_mtmatr  class                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

  /** Matrix-vector multiplication.
      Matrices are stored using MTMXDs,
      and vectors are stored using MTMDDs.
      If the template type is boolean, then this
      is equivalent to pre-image computation.
  */
  template <typename RTYPE>
  class mtmatr_mtvect : public relXset_mdd {
    public:
      mtmatr_mtvect(binary_opname* opcode, expert_forest* arg1,
        expert_forest* arg2, expert_forest* res, binary_operation* acc)
        : relXset_mdd(opcode, arg1, arg2, res, acc) { }

    protected:
      virtual node_handle processTerminals(node_handle mdd, node_handle mxd)
      {
        RTYPE mddval;
        RTYPE mxdval;
        RTYPE rval;
        argV->getValueFromHandle(mdd, mddval);
        argM->getValueFromHandle(mxd, mxdval);
        rval = mddval * mxdval;
        return resF->handleForValue(rval);
      }

  };

  template <>
  class mtmatr_mtvect<bool> : public relXset_mdd {
    public:
      mtmatr_mtvect(binary_opname* opcode, expert_forest* arg1,
        expert_forest* arg2, expert_forest* res, binary_operation* acc)
        : relXset_mdd(opcode, arg1, arg2, res, acc) { }

    protected:
      virtual node_handle processTerminals(node_handle mdd, node_handle mxd)
      {
        bool mddval;
        bool mxdval;
        bool rval;
        argV->getValueFromHandle(mdd, mddval);
        argM->getValueFromHandle(mxd, mxdval);
        rval = mddval && mxdval;
        return resF->handleForValue(rval);
      }

  };
};


// ******************************************************************
// *                                                                *
// *                      mtvect_mtmatr  class                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

  /** Vector-matrix multiplication.
      Vectors are stored using MTMDDs, and matrices
      are stored using MTMXDs.
      If the template type is boolean, then this
      is equivalent to post-image computation.
  */
  template <typename RTYPE>
  class mtvect_mtmatr : public setXrel_mdd {
    public:
      mtvect_mtmatr(binary_opname* opcode, expert_forest* arg1,
        expert_forest* arg2, expert_forest* res, binary_operation* acc)
        : setXrel_mdd(opcode, arg1, arg2, res, acc) { }

    protected:
      virtual node_handle processTerminals(node_handle mdd, node_handle mxd)
      {
        RTYPE mddval;
        RTYPE mxdval;
        RTYPE rval;
        argV->getValueFromHandle(mdd, mddval);
        argM->getValueFromHandle(mxd, mxdval);
        rval = mddval * mxdval;
        return resF->handleForValue(rval);
      }

  };

  template <>
  class mtvect_mtmatr<bool> : public setXrel_mdd {
    public:
      mtvect_mtmatr(binary_opname* opcode, expert_forest* arg1,
        expert_forest* arg2, expert_forest* res, binary_operation* acc)
        : setXrel_mdd(opcode, arg1, arg2, res, acc) { }

    protected:
      virtual node_handle processTerminals(node_handle mdd, node_handle mxd)
      {
        bool mddval;
        bool mxdval;
        bool rval;
        argV->getValueFromHandle(mdd, mddval);
        argM->getValueFromHandle(mxd, mxdval);
        rval = mddval && mxdval;
        return resF->handleForValue(rval);
      }

  };
};

// ******************************************************************
// *                                                                *
// *                      image_op_evplus class                     *
// *                                                                *
// ******************************************************************

/// Abstract base class for all MT-based pre/post image operations.
class MEDDLY::image_op_evplus : public binary_operation {
  public:
    image_op_evplus(binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

    inline ct_entry_key*
    findResult(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle &resEvmdd)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(evmdd);
      CTsrch->writeN(mxd);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      resEv = CTresult[0].readL();
      resEvmdd = resF->linkNode(CTresult[0].readN());
      if (resEvmdd != 0) {
        resEv += ev;
      }
      CT0->recycle(CTsrch);
      return 0;
    }
    inline void saveResult(ct_entry_key* Key,
      long ev, node_handle evmdd, node_handle mxd, long resEv, node_handle resEvmdd)
    {
      CTresult[0].reset();
      CTresult[0].writeL(resEvmdd == 0 ? 0L : resEv - ev);
      CTresult[0].writeN(resEvmdd);
      CT0->addEntry(Key, CTresult[0]);
    }
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);
    virtual void compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
  protected:
    binary_operation* accumulateOp;
    virtual void compute_rec(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd) = 0;

    expert_forest* argV;
    expert_forest* argM;
};


MEDDLY::image_op_evplus::image_op_evplus(binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res, binary_operation* acc)
: binary_operation(oc, 1, a1, a2, res)
{
  accumulateOp = acc;

  argV = a1;
  argM = a2;

  ct_entry_type* et = new ct_entry_type(oc->getName(), "NN:LN");
  et->setForestForSlot(0, a1);
  et->setForestForSlot(1, a2);
  et->setForestForSlot(4, res);
  registerEntryType(0, et);
  buildCTs();
}

void MEDDLY::image_op_evplus::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag)
{
  long cev = Inf<long>();
  node_handle cnode = 0;
  if (a.getForest() == argV) {
    long aev = Inf<long>();
    a.getEdgeValue(aev);
    compute(aev, a.getNode(), b.getNode(), cev, cnode);
  } else {
    long bev = Inf<long>();
    b.getEdgeValue(bev);
    compute(bev, b.getNode(), a.getNode(), cev, cnode);
  }
  c.set(cnode, cev);
}

void MEDDLY::image_op_evplus::compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  MEDDLY_DCASSERT(accumulateOp);
  compute_rec(ev, evmdd, mxd, resEv, resEvmdd);
}


// ******************************************************************
// *                                                                *
// *                         image_op_mxd class                         *
// *                                                                *
// ******************************************************************

enum RV
{
 G=-1,
 E= -2,
 L= -3,
 N= -4

};
/// Abstract base class for all MT-based pre/post image operations.
class MEDDLY::image_op_mxd : public binary_operation {
  public:
    image_op_mxd(binary_opname* opcode, expert_forest* arg1,
        expert_forest* arg2, expert_forest* res,expert_forest* res2,expert_forest* res3,expert_forest* res4,binary_operation* acc);

    inline ct_entry_key*
    findResult(node_handle a, node_handle b, node_handle &c,node_handle &d,node_handle &e,node_handle &f,int above, int& below)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CTsrch->writeN(b);
       CTsrch->writeI(above);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      d = resF->linkNode(CTresult[0].readN());
      e = resF->linkNode(CTresult[0].readN());
      f = resF->linkNode(CTresult[0].readN());
      below = CTresult[0].readI();
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(ct_entry_key* Key,
        node_handle a, node_handle b, node_handle c, node_handle d,node_handle e,node_handle f,int above,int below)
    {
      CTresult[0].reset();
      CTresult[0].writeN(c);
      CTresult[0].writeN(d);
      CTresult[0].writeN(e);
      CTresult[0].writeN(f);
      CTresult[0].writeI(below);
      CT0->addEntry(Key, CTresult[0]);
      return c;
    }
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);
    virtual void computeDDEdgeSC(const dd_edge& a, const dd_edge& b, dd_edge &c,dd_edge& d, bool userFlag,std::list<int>* shouldConfirm);
    virtual void computeDDEdgeSC(const dd_edge& a, const dd_edge& b, dd_edge &c,dd_edge& d,dd_edge &e,dd_edge &f,int& g, bool userFlag,std::list<int>* shouldConfirm,markcmp* cmp);

    virtual void compute(node_handle a, node_handle b,node_handle& c);
    virtual void computeSC(node_handle a, node_handle b,node_handle& c,node_handle& d,node_handle& e,node_handle& g,int& h, std::list<int>* shouldConfirm,markcmp* cmp);
  protected:
    binary_operation* accumulateOp;
    virtual void compute_rec(node_handle a, node_handle b,node_handle& c) = 0;
    virtual void compute_recSC(node_handle a, node_handle b,node_handle& c,node_handle& d,node_handle& e,node_handle& f,node_handle&g, std::list<int>* shouldConfirm,markcmp* cmp, RV& above, RV& below) = 0;
    virtual void compute_recSC(node_handle a, node_handle b,node_handle& c,node_handle& d,std::list<int>* shouldConfirm) = 0;

    expert_forest* argV;
    expert_forest* argM;
};

MEDDLY::image_op_mxd::image_op_mxd(binary_opname* oc, expert_forest* a1,
    expert_forest* a2, expert_forest* res,expert_forest* res2,expert_forest* res3,expert_forest* res4, binary_operation* acc)
: binary_operation(oc, 1, a1, a2, res,res2,res3,res4)
{
    accumulateOp = acc;

    argV = a1;
    argM = a2;
    int above, below;
  ct_entry_type* et = new ct_entry_type(oc->getName(), "NNI:NNNNI");
  et->setForestForSlot(0, argV);
  et->setForestForSlot(1, argM);
  et->setForestForSlot(4, res);
  et->setForestForSlot(5, res2);
  et->setForestForSlot(6, res3);
  et->setForestForSlot(7, res4);
  registerEntryType(0, et);
  buildCTs();
}

void MEDDLY::image_op_mxd
::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag)
{
  node_handle cnode;
  if (a.getForest() == argV) {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", a.getNode(), b.getNode());
#endif
    compute(a.getNode(), b.getNode(),cnode);
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", a.getNode(), b.getNode(), cnode);
#endif
  } else {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", b.getNode(), a.getNode());
#endif
    compute(b.getNode(), a.getNode(),cnode);
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", b.getNode(), a.getNode(), cnode);
#endif
  }
  c.set(cnode);
}
void MEDDLY::image_op_mxd
::computeDDEdgeSC(const dd_edge &a, const dd_edge &b, dd_edge &c,dd_edge &d, dd_edge &e,dd_edge&f,int &g, bool userFlag,std::list<int>* shouldConfirm,markcmp* cmp)
{

   node_handle cnode=0;
   node_handle dnode=0;
   node_handle enode=0;
   node_handle fnode=0;
  if (a.getForest() == argV) {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", a.getNode(), b.getNode());
#endif
    computeSC(a.getNode(), b.getNode(),cnode,dnode,enode,fnode,g,shouldConfirm,cmp);
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", a.getNode(), b.getNode(), cnode);
#endif
  } else {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", b.getNode(), a.getNode());
#endif
    compute(b.getNode(), a.getNode(),cnode);
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", b.getNode(), a.getNode(), cnode);
#endif
  }
  c.set(cnode);
  d.set(dnode);
  e.set(enode);
  f.set(fnode);
}

void MEDDLY::image_op_mxd
::computeDDEdgeSC(const dd_edge &a, const dd_edge &b, dd_edge &c,dd_edge &d, bool userFlag,std::list<int>* shouldConfirm)
{
  node_handle cnode=0;
   node_handle dnode=0;
   node_handle enode=0;
   node_handle fnode=0;

   int g=0;
  if (a.getForest() == argV) {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", a.getNode(), b.getNode());
#endif
    computeSC(a.getNode(), b.getNode(),cnode,dnode,enode,fnode,g,shouldConfirm,0);
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", a.getNode(), b.getNode(), cnode);
#endif
  } else {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", b.getNode(), a.getNode());
#endif
    compute(b.getNode(), a.getNode(),cnode);
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", b.getNode(), a.getNode(), cnode);
#endif
  }
  c.set(cnode);
  d.set(dnode);
}

void MEDDLY::image_op_mxd::compute(node_handle a, node_handle b,node_handle& c)
{
  MEDDLY_DCASSERT(accumulateOp);
  compute_rec(a, b,c);
}

void MEDDLY::image_op_mxd::computeSC(node_handle a, node_handle b,node_handle& c,node_handle& d,node_handle& e,node_handle& f ,int& h,std::list<int>* shouldConfirm,markcmp* cmp)
{
  MEDDLY_DCASSERT(accumulateOp);
  if(cmp!=0){
      RV rva;
      RV rvb;
      rva=E;
      rvb=E;
      compute_recSC(a, b,c,d,e,f,h,shouldConfirm,cmp,rva,rvb);

    }
  else{

  compute_recSC(a, b,c,d,shouldConfirm);
    }
}

class MEDDLY::covtcmxd : public binary_operation {
  public:
    covtcmxd(binary_opname* opcode, expert_forest* arg1,
        expert_forest* arg2, expert_forest* res,binary_operation* acc);

    inline ct_entry_key*
    findResult(node_handle a, node_handle b, int above,node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CTsrch->writeN(b);
       CTsrch->writeI(above);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(ct_entry_key* Key,
        node_handle a, node_handle b, int above,node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
      return c;
    }
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);
    virtual void compute(node_handle a, node_handle b,node_handle& c);

  protected:
    binary_operation* accumulateOp;
    virtual void compute_rec(node_handle a, node_handle b,node_handle& c, RV & above) = 0;

    expert_forest* argV;
    expert_forest* argM;
};

MEDDLY::covtcmxd::covtcmxd(binary_opname* oc, expert_forest* a1,
    expert_forest* a2, expert_forest* res, binary_operation* acc)
: binary_operation(oc, 1, a1, a2, res)
{
    accumulateOp = acc;

    argV = a1;
    argM = a2;
    int above;
  ct_entry_type* et = new ct_entry_type(oc->getName(), "NNI:N");
  et->setForestForSlot(0, argV);
  et->setForestForSlot(1, argM);
  et->setForestForSlot(4, res);
  registerEntryType(0, et);
  buildCTs();
}

void MEDDLY::covtcmxd
::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag)
{
    printf("INSIDE COVTCMXD computeDDEdge\n");
  node_handle cnode;
  if (a.getForest() == argV) {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", a.getNode(), b.getNode());
#endif
    compute(a.getNode(), b.getNode(),cnode);
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", a.getNode(), b.getNode(), cnode);
#endif
  } else {
#ifdef TRACE_ALL_OPS
    printf("computing top-level product(%d, %d)\n", b.getNode(), a.getNode());
#endif
    compute(b.getNode(), a.getNode(),cnode);
#ifdef TRACE_ALL_OPS
    printf("computed top-level product(%d, %d) = %d\n", b.getNode(), a.getNode(), cnode);
#endif
  }
  c.set(cnode);
}

void MEDDLY::covtcmxd::compute(node_handle a, node_handle b,node_handle& c)
{
  MEDDLY_DCASSERT(accumulateOp);
  printf("INSIDE COVTCMXD compute\n" );
  RV above=E;
  compute_rec(a, b,c,above);
}

// ******************************************************************
// *                                                                *
// *                     relXset_evplus  class                      *
// *                                                                *
// ******************************************************************

/** Generic base for relation multiplied by set.
    Changing what happens at the terminals can give
    different meanings to this operation :^)
*/
class MEDDLY::relXset_evplus : public image_op_evplus {
  public:
    relXset_evplus(binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

  protected:
    virtual void compute_rec(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
    virtual void processTerminals(long ev, node_handle mdd, node_handle mxd, long& resEv, node_handle& resEvmdd) = 0;
};

MEDDLY::relXset_evplus::relXset_evplus(binary_opname* oc,
  expert_forest* a1, expert_forest* a2, expert_forest* res, binary_operation* acc)
: image_op_evplus(oc, a1, a2, res, acc)
{
}

void MEDDLY::relXset_evplus::compute_rec(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  // termination conditions
  if (mxd == 0 || evmdd == 0) {
    resEv = 0;
    resEvmdd = 0;
    return;
  }
  if (argM->isTerminalNode(mxd)) {
    if (argV->isTerminalNode(evmdd)) {
      processTerminals(ev, evmdd, mxd, resEv, resEvmdd);
      return;
    }
    // mxd is identity
    if (argV == resF) {
      resEv = ev;
      resEvmdd = resF->linkNode(evmdd);
      return;
    }
  }

  // check the cache
  ct_entry_key* Key = findResult(ev, evmdd, mxd, resEv, resEvmdd);
  if (0==Key) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int evmddLevel = argV->getNodeLevel(evmdd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), evmddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize evmdd reader
  unpacked_node *A = (evmddLevel < rLevel)
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmdd, true)
    : argV->newUnpacked(evmdd, FULL_ONLY);

  if (evmddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      long nev = Inf<long>();
      node_handle newstates = 0;
      compute_rec(A->ei(i), A->d(i), mxd, nev, newstates);

      C->setEdge(i, newstates == 0 ? 0L : ev + nev);
      C->d_ref(i) = newstates;
    }
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= evmddLevel);

    // clear out result (important!)
    for (unsigned i=0; i<rSize; i++) {
      C->setEdge(i, 0L);
      C->d_ref(i) = 0;
    }

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New();
    unpacked_node *Rp = unpacked_node::New();
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, false);
    } else {
      argM->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    dd_edge newstatesE(resF), cdi(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getNNZs(); iz++) {
      unsigned i = Ru->i(iz);
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->d(iz), false);
      } else {
        argM->unpackNode(Rp, Ru->d(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getNNZs(); jz++) {
        unsigned j = Rp->i(jz);
        if (0==A->d(j))   continue;
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        long nev = Inf<long>();
        node_handle newstates = 0;
        compute_rec(A->ei(j), A->d(j), Rp->d(jz), nev, newstates);
        if (0==newstates) continue;
        nev += ev;
        if (0==C->d(i)) {
          C->setEdge(i, nev);
          C->d_ref(i) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        newstatesE.set(newstates, nev);
        cdi.set(C->d(i), C->ei(i));
        accumulateOp->computeTemp(newstatesE, cdi, cdi);
        C->set_de(i, cdi);
      } // for j

    } // for i

    unpacked_node::recycle(Rp);
    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);

  resF->createReducedNode(-1, C, resEv, resEvmdd);
#ifdef TRACE_ALL_OPS
  printf("computed new relXset(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmdd, mxd, resEv, resEvmdd);
#endif
  saveResult(Key, ev, evmdd, mxd, resEv, resEvmdd);
}

// ******************************************************************
// *                                                                *
// *                     setXrel_evplus  class                      *
// *                                                                *
// ******************************************************************

/** Generic base for set multiplied by relation.
    Changing what happens at the terminals can give
    different meanings to this operation :^)
*/
class MEDDLY::setXrel_evplus : public image_op_evplus {
  public:
    setXrel_evplus(binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

  protected:
    virtual void compute_rec(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
    virtual void processTerminals(long ev, node_handle mdd, node_handle mxd, long& resEv, node_handle& resEvmdd) = 0;
};

MEDDLY::setXrel_evplus::setXrel_evplus(binary_opname* oc,
  expert_forest* a1, expert_forest* a2, expert_forest* res, binary_operation* acc)
: image_op_evplus(oc, a1, a2, res, acc)
{
}

void MEDDLY::setXrel_evplus::compute_rec(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  // termination conditions
  if (mxd == 0 || evmdd == 0) {
    resEv = 0;
    resEvmdd = 0;
    return;
  }
  if (argM->isTerminalNode(mxd)) {
    if (argV->isTerminalNode(evmdd)) {
      processTerminals(ev, evmdd, mxd, resEv, resEvmdd);
      return;
    }
    // mxd is identity
    if (argV == resF) {
      resEv = ev;
      resEvmdd = resF->linkNode(evmdd);
      return;
    }
  }

  // check the cache
  ct_entry_key* Key = findResult(ev, evmdd, mxd, resEv, resEvmdd);
  if (0==Key) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int evmddLevel = argV->getNodeLevel(evmdd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), evmddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize evmdd reader
  unpacked_node *A = (evmddLevel < rLevel)
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmdd, true)
    : argV->newUnpacked(evmdd, FULL_ONLY);

  if (evmddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      long nev = Inf<long>();
      node_handle newstates = 0;
      compute_rec(A->ei(i), A->d(i), mxd, nev, newstates);

      C->setEdge(i, newstates == 0 ? 0L : ev + nev);
      C->d_ref(i) = newstates;
    }
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= evmddLevel);

    // clear out result (important!)
    for (unsigned i=0; i<rSize; i++) {
      C->setEdge(i, 0L);
      C->d_ref(i) = 0;
    }

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New();
    unpacked_node *Rp = unpacked_node::New();
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, false);
    } else {
      argM->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    dd_edge newstatesE(resF), cdj(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getNNZs(); iz++) {
      unsigned i = Ru->i(iz);
      if (0==A->d(i))   continue;
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->d(iz), false);
      } else {
        argM->unpackNode(Rp, Ru->d(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getNNZs(); jz++) {
        unsigned j = Rp->i(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        long nev = Inf<long>();
        node_handle newstates = 0;
        compute_rec(A->ei(i), A->d(i), Rp->d(jz), nev, newstates);
        if (0==newstates) continue;
        nev += ev;
        if (0==C->d(j)) {
          C->setEdge(j, nev);
          C->d_ref(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        newstatesE.set(newstates, nev);
        cdj.set(C->d(j), C->ei(j));
        accumulateOp->computeTemp(newstatesE, cdj, cdj);
        C->set_de(j, cdj);
      } // for j

    } // for i

    unpacked_node::recycle(Rp);
    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);

  resF->createReducedNode(-1, C, resEv, resEvmdd);
#ifdef TRACE_ALL_OPS
  printf("computed new setXrel(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmdd, mxd, resEv, resEvmdd);
#endif
  saveResult(Key, ev, evmdd, mxd, resEv, resEvmdd);
}

// ******************************************************************
// *                                                                *
// *                    mtmatr_evplusvect  class                    *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

  /** Matrix-vector multiplication.
      Vectors are stored using EV+MDDs, and matrices
      are stored using MTMXDs.
      If the template type is boolean, then this
      is equivalent to post-image computation.
  */
  template <typename RTYPE>
  class mtmatr_evplusvect : public relXset_evplus {
    public:
    mtmatr_evplusvect(binary_opname* opcode, expert_forest* arg1,
        expert_forest* arg2, expert_forest* res, binary_operation* acc)
        : relXset_evplus(opcode, arg1, arg2, res, acc) { }

    protected:
      virtual void processTerminals(long ev, node_handle mdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
      {
        RTYPE evmddval;
        RTYPE mxdval;
        RTYPE rval;
        argV->getValueFromHandle(mdd, evmddval);
        argM->getValueFromHandle(mxd, mxdval);
        rval = evmddval * mxdval;
        resEv = ev;
        resEvmdd = resF->handleForValue(rval);
      }
  };
};

// ******************************************************************
// *                                                                *
// *                    evplusvect_mtmatr  class                    *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

  /** Vector-matrix multiplication.
      Vectors are stored using EV+MDDs, and matrices
      are stored using MTMXDs.
      If the template type is boolean, then this
      is equivalent to post-image computation.
  */
  template <typename RTYPE>
  class evplusvect_mtmatr : public setXrel_evplus {
    public:
    evplusvect_mtmatr(binary_opname* opcode, expert_forest* arg1,
        expert_forest* arg2, expert_forest* res, binary_operation* acc)
        : setXrel_evplus(opcode, arg1, arg2, res, acc) { }

    protected:
      virtual void processTerminals(long ev, node_handle mdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
      {
        RTYPE evmddval;
        RTYPE mxdval;
        RTYPE rval;
        argV->getValueFromHandle(mdd, evmddval);
        argM->getValueFromHandle(mxd, mxdval);
        rval = evmddval * mxdval;
        resEv = ev;
        resEvmdd = resF->handleForValue(rval);
      }
  };
};

// ******************************************************************
// *                                                                *
// *                     tcXrel_evplus  class                      *
// *                                                                *
// ******************************************************************

/** Generic base for transitive closure multiplied by relation.
    Changing what happens at the terminals can give
    different meanings to this operation :^)
*/
class MEDDLY::tcXrel_evplus : public image_op_evplus {
  public:
    tcXrel_evplus(binary_opname* opcode, expert_forest* tc,
      expert_forest* trans, expert_forest* res, binary_operation* acc);

  protected:
    virtual void compute_rec(long ev, node_handle evmxd, node_handle mxd, long& resEv, node_handle& resEvmxd);
    virtual void processTerminals(long ev, node_handle evmxd, node_handle mxd, long& resEv, node_handle& resEvmxd);
};

MEDDLY::tcXrel_evplus::tcXrel_evplus(binary_opname* oc,
  expert_forest* tc, expert_forest* trans, expert_forest* res, binary_operation* acc)
: image_op_evplus(oc, tc, trans, res, acc)
{
}

void MEDDLY::tcXrel_evplus::compute_rec(long ev, node_handle evmxd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  // termination conditions
  if (mxd == 0 || evmxd == 0) {
    resEv = 0;
    resEvmdd = 0;
    return;
  }
  if (argM->isTerminalNode(mxd)) {
    if (argV->isTerminalNode(evmxd)) {
      processTerminals(ev, evmxd, mxd, resEv, resEvmdd);
      return;
    }
    // mxd is identity
    if (argV == resF) {
      resEv = ev;
      resEvmdd = resF->linkNode(evmxd);
      return;
    }
  }

  // check the cache
  ct_entry_key* Key = findResult(ev, evmxd, mxd, resEv, resEvmdd);
  if (0==Key) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int evmxdLevel = argV->getNodeLevel(evmxd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), ABS(evmxdLevel));
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(rLevel, evmxdLevel)
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmxd, true)
    : argV->newUnpacked(evmxd, FULL_ONLY);

  for (unsigned i = 0; i < rSize; i++) {
    int pLevel = argV->getNodeLevel(A->d(i));
    unpacked_node* B = isLevelAbove(-rLevel, pLevel)
      ? unpacked_node::newIdentity(argV, -rLevel, i, 0L, A->d(i), true)
      : argV->newUnpacked(A->d(i), FULL_ONLY);

    unpacked_node* D = unpacked_node::newFull(resF, -rLevel, rSize);
    if (rLevel > ABS(mxdLevel)) {
      //
      // Skipped levels in the MXD,
      // that's an important special case that we can handle quickly.
      for (unsigned j = 0; j < rSize; j++) {
        long nev = Inf<long>();
        node_handle newstates = 0;
        compute_rec(A->ei(i) + B->ei(j), B->d(j), mxd, nev, newstates);

        D->setEdge(j, newstates == 0 ? 0L : ev + nev);
        D->d_ref(j) = newstates;
      }
    }
    else {
      //
      // Need to process this level in the MXD.
      MEDDLY_DCASSERT(ABS(mxdLevel) >= ABS(pLevel));

      // clear out result (important!)
      for (unsigned j = 0; j < rSize; j++) {
        D->setEdge(j, 0L);
        D->d_ref(j) = 0;
      }

      // Initialize mxd readers, note we might skip the unprimed level
      unpacked_node *Ru = unpacked_node::New();
      unpacked_node *Rp = unpacked_node::New();
      if (mxdLevel < 0) {
        Ru->initRedundant(argM, rLevel, mxd, false);
      } else {
        argM->unpackNode(Ru, mxd, SPARSE_ONLY);
      }

      dd_edge newstatesE(resF), djp(resF);

      // loop over mxd "rows"
      for (unsigned jz = 0; jz < Ru->getNNZs(); jz++) {
        unsigned j = Ru->i(jz);
        if (0 == B->d(j)) {
          continue;
        }

        if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(jz)))) {
          Rp->initIdentity(argM, rLevel, j, Ru->d(jz), false);
        } else {
          argM->unpackNode(Rp, Ru->d(jz), SPARSE_ONLY);
        }

        // loop over mxd "columns"
        for (unsigned jpz = 0; jpz < Rp->getNNZs(); jpz++) {
          unsigned jp = Rp->i(jpz);
          // ok, there is an i->j "edge".
          // determine new states to be added (recursively)
          // and add them
          long nev = Inf<long>();
          node_handle newstates = 0;
          compute_rec(A->ei(i) + B->ei(j), B->d(j), Rp->d(jpz), nev, newstates);
          if (0==newstates) {
            continue;
          }
          nev += ev;
          if (0 == D->d(jp)) {
            D->setEdge(jp, nev);
            D->d_ref(jp) = newstates;
            continue;
          }
          // there's new states and existing states; union them.
          newstatesE.set(newstates, nev);
          djp.set(D->d(jp), D->ei(jp));
          accumulateOp->computeTemp(newstatesE, djp, djp);
          D->set_de(jp, djp);
        } // for j

      } // for i

      unpacked_node::recycle(Rp);
      unpacked_node::recycle(Ru);
    } // else

    long cev = Inf<long>();
    node_handle cnode = 0;
    resF->createReducedNode(int(i), D, cev, cnode);
    C->setEdge(i, cev);
    C->d_ref(i) = cnode;

    unpacked_node::recycle(B);
  }

  // cleanup mdd reader
  unpacked_node::recycle(A);

  resF->createReducedNode(-1, C, resEv, resEvmdd);
#ifdef TRACE_ALL_OPS
  printf("computed new tcXrel(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmxd, mxd, resEv, resEvmdd);
#endif
  saveResult(Key, ev, evmxd, mxd, resEv, resEvmdd);
}

void MEDDLY::tcXrel_evplus::processTerminals(long ev, node_handle evmxd, node_handle mxd, long& resEv, node_handle& resEvmxd)
{
  long evmddval;
  long mxdval;
  long rval;
  argV->getValueFromHandle(evmxd, evmddval);
  argM->getValueFromHandle(mxd, mxdval);
  rval = evmddval * mxdval;
  resEv = ev;
  resEvmxd = resF->handleForValue(rval);
}

// ******************************************************************
// *                                                                *
// *                         covtc  class                           *
// *                                                                *
// ******************************************************************


class MEDDLY::covtc : public covtcmxd {
  public:
    covtc(binary_opname* opcode, expert_forest* tc,
      expert_forest* trans, expert_forest* res, binary_operation* acc);

  protected:
      virtual void compute_rec(node_handle a, node_handle b, node_handle&c, RV& above);
    virtual void processTerminals(node_handle evmxd, node_handle mxd, RV& above,node_handle& resEvmxd);
};

MEDDLY::covtc::covtc(binary_opname* oc,
  expert_forest* tc, expert_forest* trans, expert_forest* res, binary_operation* acc)
: covtcmxd(oc, tc, trans, res, acc)
{
}

void MEDDLY::covtc::compute_rec(node_handle cr, node_handle e, node_handle&res, RV& above)
{
    printf("COVTC START\n" );
//   // termination conditions
  if (cr == 0 || e == 0) {
    res=0;
    return;
  }
  if (argM->isTerminalNode(e)) {
    if (argV->isTerminalNode(cr)) {
      processTerminals(cr, e, above, res);
      return;
    }

  }
//
  // check the cache
  ct_entry_key* Key = findResult(cr,e, above,res);
  if (0==Key) {
    return;
  }
//
  // check if mxd and evmdd are at the same level
  const int crLevel = argV->getNodeLevel(cr);
  const int eLevel = argM->getNodeLevel(e);
  const int rLevel = MAX(ABS(eLevel), ABS(crLevel));
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize cr reader
  unpacked_node* A = isLevelAbove(rLevel, crLevel)
    ? unpacked_node::newRedundant(argV, rLevel, cr, true)
    : argV->newUnpacked(cr, FULL_ONLY);

  for (unsigned i = 0; i < rSize; i++) {
    int pLevel = argV->getNodeLevel(A->d(i));
    unpacked_node* B = isLevelAbove(-rLevel, pLevel)
      ? unpacked_node::newIdentity(argV, -rLevel, i, A->d(i), true)
      : argV->newUnpacked(A->d(i), FULL_ONLY);

    unpacked_node* D = unpacked_node::newFull(resF, -rLevel, rSize);
    if (rLevel > ABS(eLevel)) {
      //
      // Skipped levels in the MXD,
      // that's an important special case that we can handle quickly.
      for (unsigned j = 0; j < rSize; j++) {
        node_handle newstates = 0;
        RV newabove=above;
        if(i!=j){newabove=N;}
        compute_rec( B->d(j), e, newstates,newabove);

        D->d_ref(j) = newstates;
      }
    }
    else {
      //
      // Need to process this level in the MXD.
      MEDDLY_DCASSERT(ABS(eLevel) >= ABS(pLevel));

      // clear out result (important!)
      for (unsigned j = 0; j < rSize; j++) {
        D->d_ref(j) = 0;
      }

      // Initialize mxd readers, note we might skip the unprimed level
      unpacked_node *Ru = unpacked_node::New();
      unpacked_node *Rp = unpacked_node::New();
      if (eLevel < 0) {
        Ru->initRedundant(argM, rLevel, e, false);
      } else {
        argM->unpackNode(Ru, e, SPARSE_ONLY);
      }

      dd_edge newstatesE(resF), djp(resF);

      // loop over mxd "rows"
      for (unsigned jz = 0; jz < Ru->getNNZs(); jz++) {
        unsigned j = Ru->i(jz);
        if (0 == B->d(j)) {
          continue;
        }

        if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(jz)))) {
          Rp->initIdentity(argM, rLevel, j, Ru->d(jz), false);
        } else {
          argM->unpackNode(Rp, Ru->d(jz), SPARSE_ONLY);
        }

        // loop over mxd "columns"
        for (unsigned jpz = 0; jpz < Rp->getNNZs(); jpz++) {
          unsigned jp = Rp->i(jpz);
          // ok, there is an i->j "edge".
          // determine new states to be added (recursively)
          // and add them
          RV newabove=above;
          if(i!=j)newabove=N;
          node_handle newstates = 0;
          compute_rec(B->d(j), Rp->d(jpz), newstates,newabove);
          if (0==newstates) {
            continue;
          }
          if (0 == D->d(jp)) {
            D->d_ref(jp) = newstates;
            continue;
          }
          // there's new states and existing states; union them.
          newstatesE.set(newstates);
          djp.set(D->d(jp));
          accumulateOp->computeTemp(newstatesE, djp, djp);
          D->set_d(jp, djp);
            } // for jpz

        } // for jz

      unpacked_node::recycle(Rp);
      unpacked_node::recycle(Ru);
    } // else

    node_handle cnode = 0;
    cnode=resF->createReducedNode(int(i), D);
    C->d_ref(i) = cnode;

    unpacked_node::recycle(B);
  }

  // cleanup mdd reader
  unpacked_node::recycle(A);

  res=resF->createReducedNode(-1, C);
// #ifdef TRACE_ALL_OPS
//   printf("computed new tcXrel(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmxd, mxd, resEv, resEvmdd);
// #endif
  saveResult(Key, cr, e, above,res);
}

void MEDDLY::covtc::processTerminals( node_handle cr, node_handle e,RV& above, node_handle& res)
{
  long evmddval;
  long mxdval;
  long rval;
  argV->getValueFromHandle(cr, evmddval);
  argM->getValueFromHandle(e, mxdval);
  rval = evmddval * mxdval;
  if(above==E)
  {
      res=0;
  }else{
      res = resF->handleForValue(rval);
  }
}


// ******************************************************************
// *                                                                *
// *                   mrrc  class                                  *
// *                                                                *
// ******************************************************************

/** Generic base for maintain reachability relation.
    Changing what happens at the terminals can give
    different meanings to this operation :^)
*/

class MEDDLY::mrrc : public image_op_mxd {
  public:
    mrrc(binary_opname* opcode, expert_forest* tc,
        expert_forest* trans, expert_forest* res,expert_forest* res2,expert_forest* res3,expert_forest* res4, binary_operation* acc);

  protected:
      int l=-1;
      virtual void compute_rec(node_handle a, node_handle b, node_handle&c);
      virtual void compute_recSC(node_handle a, node_handle b, node_handle&c,node_handle&d,node_handle&e,node_handle& g, int&f,std::list<int>* shouldConfirm,markcmp* cmp, RV& above, RV& below);
      virtual void compute_recSC(node_handle a, node_handle b, node_handle&c,node_handle&d,std::list<int>* shouldConfirm);
      virtual void processTerminals(node_handle evmxd, node_handle mxd, node_handle& resEvmxd,node_handle& eq,node_handle& leq,node_handle& tleq,RV& above);
      virtual void processTerminals(node_handle evmxd, node_handle mxd, node_handle& resEvmxd,node_handle& eq,node_handle& leq,node_handle& tleq);
    void AddtoNH(unpacked_node* D, int i, node_handle s );
};

MEDDLY::mrrc::mrrc(binary_opname* oc,
expert_forest* tc, expert_forest* trans, expert_forest* res,expert_forest* res2,expert_forest* res3,expert_forest* res4, binary_operation* acc)
: image_op_mxd(oc, tc, trans, res,res2,res3,res4, acc)
{
}

void MEDDLY::mrrc::compute_rec(node_handle a, node_handle b,node_handle&c){printf("SHOULD Implemented!! return 0\n" );getchar();return ;}

void MEDDLY::mrrc::compute_recSC(node_handle evmxd, node_handle mxd, node_handle& resEvmdd,node_handle& eq,std::list<int>* shouldConfirm)
{
  // termination conditions
  if (mxd == 0 || evmxd == 0) {
    resEvmdd = 0;
    eq=0;
    return;
  }
  if (argM->isTerminalNode(mxd)) {
    if (argV->isTerminalNode(evmxd)) {
      processTerminals( evmxd, mxd, resEvmdd,eq,eq,eq);
      return;
    }
  }

  // check the cache
  int iabove,ibelow;
  ct_entry_key* Key = findResult( evmxd, mxd, resEvmdd,eq,eq,eq,iabove,ibelow);
  if (0==Key) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int evmxdLevel = argV->getNodeLevel(evmxd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), ABS(evmxdLevel));
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);
  unpacked_node* H = unpacked_node::newFull(resF, rLevel, rSize);


  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(rLevel, evmxdLevel)
    ? unpacked_node::newRedundant(argV, rLevel, evmxd, true)
    : argV->newUnpacked(evmxd, FULL_ONLY);
     unpacked_node** nst = new unpacked_node*[rSize];
     for (unsigned i = 0; i < rSize; i++) {
         nst[i]=unpacked_node::newFull(resF, -rLevel, rSize);
     }

  for (unsigned i = 0; i < rSize; i++) {
    int pLevel = argV->getNodeLevel(A->d(i));
    unpacked_node* B = isLevelAbove(-rLevel, pLevel)
      ? unpacked_node::newIdentity(argV, -rLevel, i,  A->d(i), true)
      : argV->newUnpacked(A->d(i), FULL_ONLY);

    unpacked_node* D = unpacked_node::newFull(resF, -rLevel, rSize);

    if (rLevel > ABS(mxdLevel)) {

      // Skipped levels in the MXD,
      // that's an important special case that we can handle quickly.
      for (unsigned j = 0; j < rSize; j++) {
        node_handle newstates = 0;
        node_handle eqnewstates = 0;

        compute_recSC( B->d(j), mxd, newstates,eqnewstates,shouldConfirm);

        if(newstates!=0 && eqnewstates==0){
            printf("ERROR! %d %d\n", newstates,eqnewstates);
        }
        D->d_ref(j) = newstates;
        if(eqnewstates!=0){
        nst[j]->d_ref(j) = eqnewstates;
        }
        // else{
        //     // printf("NST0 not set\n" );
        // }
      }
    }
    else {

      // Need to process this level in the MXD.
      MEDDLY_DCASSERT(ABS(mxdLevel) >= ABS(pLevel));

      // clear out result (important!)
      for (unsigned j = 0; j < rSize; j++) {
        D->d_ref(j) = 0;
      }

      // Initialize mxd readers, note we might skip the unprimed level
      unpacked_node *Ru = unpacked_node::New();
      unpacked_node *Rp = unpacked_node::New();
      if (mxdLevel < 0) {
        Ru->initRedundant(argM, rLevel, mxd, false);
      } else {
        argM->unpackNode(Ru, mxd, SPARSE_ONLY);
      }

      dd_edge newstatesE(resF), djp(resF);
      dd_edge newstatesEg(resF), djpg(resF);

      dd_edge newstatesEQ(resF), djpq(resF);

      // loop over mxd "rows"
      for (unsigned jz = 0; jz < Ru->getNNZs(); jz++) {
        unsigned j = Ru->i(jz);
        if (0 == B->d(j)) {
          continue;
        }

        if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(jz)))) {
          Rp->initIdentity(argM, rLevel, j, Ru->d(jz), false);
        } else {
          argM->unpackNode(Rp, Ru->d(jz), SPARSE_ONLY);
        }

        // loop over mxd "columns"
        for (unsigned jpz = 0; jpz < Rp->getNNZs(); jpz++) {
          unsigned jp = Rp->i(jpz);
          // ok, there is an i->j "edge".
          // determine new states to be added (recursively)
          // and add them
          node_handle newstates = 0;
          node_handle eqnewstates = 0;

          compute_recSC( B->d(j), Rp->d(jpz), newstates,eqnewstates,shouldConfirm);
          if(newstates!=0 && eqnewstates==0){
              printf("ERROR P2! %d %d\n", newstates,eqnewstates);
              getchar();
          }
          if (0==newstates) {
            continue;
          }
          shouldConfirm[mxdLevel].push_back(jp);
          if (0 == D->d(jp)) {
            D->d_ref(jp) = newstates;
            if(nst[jp]->d(jp)==0){
                if(eqnewstates!=0)
            nst[jp]->d_ref(jp)=eqnewstates;
            }
            continue;
          }
          // there's new states and existing states; union them.
          newstatesE.set(newstates);
          djp.set(D->d(jp));
          accumulateOp->computeTemp(newstatesE, djp, djp);
          D->set_d(jp, djp);

          if(eqnewstates!=0){
          dd_edge newstatesEnst(resF), djpnst(resF);
          newstatesEnst.set(eqnewstates);
          djpnst.set(nst[jp]->d(jp));
          accumulateOp->computeTemp(newstatesEnst, djpnst, djpnst);
          nst[jp]->set_d(jp, djpnst);
            }
        } // for j
    } // for jz

      unpacked_node::recycle(Rp);
      unpacked_node::recycle(Ru);
    } // else
    node_handle cnode = 0;
    cnode=resF->createReducedNode(int(i), D);
    C->d_ref(i) = cnode;
    ostream_output meddlyout(std::cout);
    unpacked_node::recycle(B);
} //for i
for (unsigned i = 0; i < rSize; i++) {
     dd_edge newstatesEG(resF);
     node_handle nstnode = 0;
     nstnode=resF->createReducedNode(int(i), nst[i]);
     ostream_output meddlyout(std::cout);
 // nst[i]->show(meddlyout,true);

     if(H->d(i)==0){
     H->d_ref(i) = nstnode;
    }else{
        printf("WRONGWRONG!SETING TWICE!\n" );
    }
}

for (unsigned i = 0; i < rSize; i++) {
    nst[i]=0;
}
delete[] nst;

  // cleanup mdd reader
  unpacked_node::recycle(A);

  resEvmdd=resF->createReducedNode(-1, C);
  eq=resF->createReducedNode(-1, H);
#ifdef TRACE_ALL_OPS
  printf("computed new tcXrel(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmxd, mxd, resEv, resEvmdd);
#endif
  saveResult(Key, evmxd, mxd, resEvmdd,eq,eq,0,0,0);

}


void MEDDLY::mrrc::compute_recSC(node_handle evmxd, node_handle mxd, node_handle& resEvmdd,node_handle& eq,node_handle& leq,node_handle& tleq, int& geq,std::list<int>* shouldConfirm,markcmp* cmp,RV& above, RV& below)
{
    ostream_output meddlyout(std::cout);
    bool d=true;
if(d)
    printf("call compute_recSC for evmxd %d, mxd %d, above %d, below %d\n",evmxd,mxd,above,below );

// termination conditions
if (mxd == 0 || evmxd == 0) {
        resEvmdd = 0;
        eq=0;
        leq=0;
        geq=0;
        tleq=0;
        below=N;
        if(above==G) { below=G;}

        return;
}
if (argM->isTerminalNode(mxd)) {
        if (argV->isTerminalNode(evmxd)) {
                if(above==G) { below=G; resEvmdd=0; eq=0; leq=0; tleq=0; return; }
                processTerminals( evmxd, mxd, resEvmdd,eq, leq,tleq, above);
                if (d)
                printf("terminal res%d eq %d  leq %d tleq %d above %d \n",resEvmdd,eq,leq,tleq, above );
                if(above==N) { below=N; }
                else if(above==L) { below=L; }
                else {below=E;}
                return;
        }
}

// check the cache
int iabove=(int)above, igeq=0;
ct_entry_key* Key = findResult( evmxd, mxd, resEvmdd,eq,leq,tleq,iabove,igeq);
if (0==Key) {
        if(d)
            printf("FOUND IN CACHE  above %d node %d, transitionRel %d resEvmdd%d,eq%d,leq%d,\n",above,evmxd, mxd, resEvmdd,eq,leq );
        if(igeq>0) geq=igeq;
        return;
}

// check if mxd and evmdd are at the same level
const int evmxdLevel = argV->getNodeLevel(evmxd);
const int mxdLevel = argM->getNodeLevel(mxd);
const int rLevel = MAX(ABS(mxdLevel), ABS(evmxdLevel));
const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);
unpacked_node* CLEQ = unpacked_node::newFull(resF, rLevel, rSize);
unpacked_node* CTLEQ = unpacked_node::newFull(resF, rLevel, rSize);
unpacked_node* CEQ = unpacked_node::newFull(resF, rLevel, rSize);

if(d)
printf("nodeLEVEL %d TransitionLEVEL %d\n", evmxdLevel,mxdLevel);
// Initialize evmdd reader
unpacked_node* A = isLevelAbove(rLevel, evmxdLevel)
    ? unpacked_node::newRedundant(argV, rLevel, evmxd, true)
    : argV->newUnpacked(evmxd, FULL_ONLY);
RV savedAbove=above;
for (unsigned i = 0; i < rSize; i++) {

        int pLevel = argV->getNodeLevel(A->d(i));
        unpacked_node* B = isLevelAbove(-rLevel, pLevel)
      ? unpacked_node::newIdentity(argV, -rLevel, i,  A->d(i), true)
      : argV->newUnpacked(A->d(i), FULL_ONLY);

        unpacked_node* D = unpacked_node::newFull(resF, -rLevel, rSize);
        unpacked_node* DLEQ = unpacked_node::newFull(resF, -rLevel, rSize);
        // unpacked_node* DTLEQ = unpacked_node::newFull(resF, -rLevel, rSize);
        above=savedAbove;
        if (rLevel > ABS(mxdLevel)) {
                if(d) {
                        printf("IF a %d b %d mxdl %d l %d \n",above, below,argM->getNodeLevel(evmxdLevel),l );
                }

                // Skipped levels in the MXD,
                // that's an important special case that we can handle quickly.
                for (unsigned j = 0; j < rSize; j++) {
                        above=savedAbove;
                        RV newabove=above;
                        RV newbelow=below;
                        if(d)
                        {
                            printf("i %d\n",i );
                            cmp->hasOmega(i,rLevel);
                            printf("j %d \n",j );
                            cmp->hasOmega(j,rLevel);
                            printf("XXXijXXX %d\n", B->d(j) );

                        }
                        int compareijresult=cmp->compare(i,j,evmxdLevel);
                        if(d){
                            printf("^^ IF CPMPARE IJ\n" );
                            ostream_output meddlyout(std::cout);
                            A->show(meddlyout,true);
                            B->show(meddlyout,true);
                            printf("^^ IF A,B IJ END!\n" );
                        }
                        // printf("i %d j %d cij %d\n",i,j,compareijresult );

                        if(i==j) {
                                ;
                        }else if(compareijresult>0) { //i<=j
                                if(above==G || above==N) {newabove=N; }//printf("above set to N1\n" );}
                                else if(above==L || above==E) {newabove=L; if(d) printf("above set to L1\n" ); }
                        }else if(compareijresult==-1) { //i>j
                                if(above==G || above==E) newabove=G;
                                else if(above==L || above==N) {newabove=N;}//printf("above set to N2\n" );}
                        }else if(compareijresult==-4) { // i is not comparable to j!
                                newabove=N;//printf("above set to N3\n" );
                        }
                        node_handle newstates = 0;
                        node_handle eqnewstates = 0;
                        node_handle leqnewstates = 0;
                        node_handle tleqnewstates = 0;
                        int geq=0;
                        compute_recSC( B->d(j), mxd, newstates,eqnewstates,leqnewstates,tleqnewstates,geq,shouldConfirm,cmp,newabove,newbelow);

                        if(newstates!=0 || eqnewstates!=0 || leqnewstates!=0) {
                                shouldConfirm[evmxdLevel].push_back(j);
                                if(d)
                                printf("shouldConfirm [%d].pushback(%d) If\n",evmxdLevel,j );
                                if(d)
                                printf("IF LVL%d  cij %d i:%d , j: %d, above: %d, below:%d leqnewstates: %d \n",evmxdLevel,compareijresult,i,j, above,below,leqnewstates );
                                if((leqnewstates!=0)&&((above==E&&compareijresult>0)||(above==L&&compareijresult>0))) {
                                    geq=1;

                                        if(compareijresult>=rSize) {
                                                shouldConfirm[evmxdLevel].push_back(compareijresult);
                                                if(d)
                                                printf("shouldConfirm [%d].pushback(%d) IfOmega\n",evmxdLevel,compareijresult );
                                                printf("IF OMEGA2  set HERE!!\n" );
                                                getchar();
                                        }else{
                                                // if(d)
                                                printf("IF OMEGA2 inside HERE!!\n" );
                                                getchar();
                                                if(newstates!=0){
                                                    AddtoNH(D,j,newstates);
                                                    if(d)
                                                    printf("D[%d]=%d\n",j,newstates );
                                                    unpacked_node* DEQ = unpacked_node::newFull(resF, -rLevel, rSize);
                                                    AddtoNH(DEQ,j,eqnewstates);
                                                    if(d)
                                                    printf("DEQ[%d]=%d\n",j,eqnewstates );
                                                    node_handle denode = 0;
                                                    denode=resF->createReducedNode(int(j), DEQ);//, cev, cnode);
                                                    AddtoNH(CEQ,j,denode);
                                                }
                                                if(leqnewstates!=0){

                                                AddtoNH(DLEQ,compareijresult,leqnewstates);
                                                unpacked_node* DTLEQ = unpacked_node::newFull(resF, -rLevel, rSize);
                                                AddtoNH(DTLEQ,compareijresult,tleqnewstates);
                                                node_handle tlnode=0;
                                                tlnode=resF->createReducedNode(int(compareijresult), DTLEQ);
                                                AddtoNH(CTLEQ,j,tlnode);

                                                if(d)
                                                {
                                                    printf("DLEQ[%d]=%d\n",compareijresult,leqnewstates );
                                                    printf("DTLEQ[%d]=%d\n",compareijresult,leqnewstates );
                                                    printf("CTLEQ[%d]=%d\n",j,tlnode );
                                                }
                                                unpacked_node* DEQ = unpacked_node::newFull(resF, -rLevel, rSize);
                                                AddtoNH(DEQ,compareijresult,eqnewstates);
                                                if(d)
                                                printf("DEQ[%d]=%d\n",compareijresult,eqnewstates );
                                                node_handle denode = 0;
                                                denode=resF->createReducedNode(int(compareijresult), DEQ);//, cev, cnode);
                                                AddtoNH(CEQ,compareijresult,denode);
                                                }
                                        }
                                }else{
                                    if(d)
                                    printf("CAme to IF AddtoNH !\n" );
                                        AddtoNH(D,j,newstates);
                                        if(d)
                                        printf("D[%d]=%d\n",j,newstates );
                                        unpacked_node* DEQ = unpacked_node::newFull(resF, -rLevel, rSize);
                                        AddtoNH(DEQ,j,eqnewstates);
                                        if(d)
                                        printf("DEQ[%d]=%d\n",j,eqnewstates );
                                        node_handle denode = 0;
                                        denode=resF->createReducedNode(int(j), DEQ);
                                        AddtoNH(CEQ,j,denode);
                                        if((above==L &&  compareijresult==-2)||(above==E&& compareijresult==-2/* &&leqnewstates>0*/) ) {
                                            if(d)
                                            printf("CAme to IF AddtoNH ADDED!\n" );

                                                AddtoNH(DLEQ,j,leqnewstates);
                                                unpacked_node* DTLEQ = unpacked_node::newFull(resF, -rLevel, rSize);
                                                AddtoNH(DTLEQ,j,tleqnewstates);
                                                node_handle tlnode=0;
                                                tlnode=resF->createReducedNode(int(j), DTLEQ);
                                                AddtoNH(CTLEQ,j,tlnode);
                                                if(d)
                                                {
                                                    printf("DLEQ[%d]=%d\n",j,leqnewstates );
                                                    printf("DTLEQ[%d]=%d\n",j,leqnewstates );
                                                    printf("CTLEQ[%d]=%d\n",j,tlnode );
                                                }
                                        }
                                }
                        }
                }
        }else {
                if(d) {
                        printf("ELSE a %d b %d mxdl %d \n",above, below,evmxdLevel);//argM->getNodeLevel(evmxdLevel),l );
                }

                // Need to process this level in the MXD.
                MEDDLY_DCASSERT(ABS(mxdLevel) >= ABS(pLevel));

                // clear out result (important!)
                for (unsigned j = 0; j < rSize; j++) {
                        D->d_ref(j) = 0;
                        DLEQ->d_ref(j) = 0;
                }

                // Initialize mxd readers, note we might skip the unprimed level
                unpacked_node *Ru = unpacked_node::New();
                unpacked_node *Rp = unpacked_node::New();
                if (mxdLevel < 0) {
                        Ru->initRedundant(argM, rLevel, mxd, false);
                } else {
                        argM->unpackNode(Ru, mxd, SPARSE_ONLY);
                }

                dd_edge newstatesE(resF), djp(resF);
                dd_edge newstatesEg(resF), djpg(resF);

                dd_edge newstatesEQ(resF), djpq(resF);

                // loop over mxd "rows"
                for (unsigned jz = 0; jz < Ru->getNNZs(); jz++) {
                        unsigned j = Ru->i(jz);
                        if (0 == B->d(j)) {
                                continue;
                        }

                        if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(jz)))) {
                                Rp->initIdentity(argM, rLevel, j, Ru->d(jz), false);
                        } else {
                                argM->unpackNode(Rp, Ru->d(jz), SPARSE_ONLY);
                        }

                        // loop over mxd "columns"
                        for (unsigned jpz = 0; jpz < Rp->getNNZs(); jpz++) {
                                unsigned jp = Rp->i(jpz);
                                // ok, there is an i->j "edge".
                                // determine new states to be added (recursively)
                                // and add them
                                // printf("i %d jp %d\n",i,jp );
                                if(d)
                                {
                                    printf("ELSE i\n" );
                                    cmp->hasOmega(i,mxdLevel);
                                    printf("ELSE j\n" );
                                    cmp->hasOmega(jp,mxdLevel);
                                    printf("XXXELSE ijXXX\n" );

                                }
                                int compareijresult=cmp->compare(i,jp,mxdLevel);
                                if(d){
                                    printf("^^ ELSE CPMPARE IJ\n" );
                                }
                                if(Rp->d(jpz)==0) {
                                        printf("NEED TO CHECK JP!\n" );
                                }
                                above=savedAbove;
                                RV newabove=above;
                                RV newbelow=below;
                                if(d)
                                printf("above %d compareijresult %d\n",above, compareijresult );
                                if(i==jp) {
                                        newabove=above;
                                }else if(compareijresult>0) {
                                        if(above==G || above==N) {newabove=N; if(d)printf("above set to N4\n" );}
                                        else if(above==L || above==E) {newabove=L; if(d) printf("HERE! above set to L4\n" ); }
                                        if(d)
                                        printf("Inside if else newAbove %d\n",newabove );
                                }else if(compareijresult==-1) {
                                    if(d)
                                    printf("CAME IN compareijresult==-1 \n");
                                        if(above==G || above==E) newabove=G;
                                        else if(above==L || above==N) {newabove=N; if(d)printf("above set to N5\n" );}
                                }else if(compareijresult==-4) {
                                        newabove=N; if(d)printf("above set to N6\n" );
                                }
                                if(d)
                                printf("ELSE newAbove %d\n",newabove );
                                node_handle newstates = 0;
                                node_handle eqnewstates = 0;
                                node_handle leqnewstates = 0;
                                node_handle tleqnewstates = 0;
                                int geq=0;


                                compute_recSC( B->d(j), Rp->d(jpz), newstates,eqnewstates,leqnewstates,tleqnewstates,geq, shouldConfirm,cmp,newabove,newbelow);

                                if (newstates==0 && eqnewstates==0 && leqnewstates==0) {
                                        continue;
                                }
                                shouldConfirm[mxdLevel].push_back(jp);
                                if(d)
                                printf("shouldConfirm [%d].pushback(%d)Else\n",mxdLevel,jp );
                                if(d)
                                printf("ELSE LVL%d cij %d i:%d , j:%d,  jp: %d, above: %d, below:%d leqnewstates: %d \n",mxdLevel,compareijresult, i,j, jp, above,below,leqnewstates );
                                        // printf("LVL%d  cij %d i:%d , j: %d, above: %d, below:%d \n",evmxdLevel,compareijresult,i,j, above,below );
                                        if((leqnewstates!=0)&&((above==E&&compareijresult>0)||(above==L&&compareijresult>0))) {
                                            geq=1;
                                                if(compareijresult>=rSize) {
                                                        shouldConfirm[evmxdLevel].push_back(compareijresult);
                                                        if(d)
                                                        printf("shouldConfirm [%d].pushback(%d) ElseOmega!\n",mxdLevel,compareijresult );
                                                        // if(d)
                                                        printf(" OMEGA2  set HERE!!\n" );
                                                        getchar();
                                                }else{
                                                        // if(d)
                                                        printf("ELSE OMEGA2 inside HERE %d !!\n",eqnewstates );
                                                        if(d)
                                                        printf("Omega part newstates%d,leqnewstates%d,eqnewstates%d\n",newstates,leqnewstates,eqnewstates );
                                                        getchar();
                                                        if(newstates!=0){
                                                            AddtoNH(D,jp,newstates);
                                                            if(d)
                                                            printf("D[%d]=%d\n",jp,newstates );
                                                            unpacked_node* DEQ = unpacked_node::newFull(resF, -rLevel, rSize);
                                                            AddtoNH(DEQ,jp,eqnewstates);
                                                            if(d)
                                                            printf("DEQ[%d]=%d\n",jp,eqnewstates );

                                                            node_handle denode = 0;
                                                            denode=resF->createReducedNode(int(jp), DEQ);//, cev, cnode);
                                                            AddtoNH(CEQ,jp,denode);
                                                        }
                                                        if(leqnewstates!=0){
                                                        AddtoNH(DLEQ,compareijresult,leqnewstates);
                                                        unpacked_node* DTLEQ = unpacked_node::newFull(resF, -rLevel, rSize);
                                                        AddtoNH(DTLEQ,compareijresult,tleqnewstates);
                                                        node_handle tlnode=0;
                                                        tlnode=resF->createReducedNode(int(jp), DTLEQ);
                                                        AddtoNH(CTLEQ,jp,tlnode);
                                                        if(d)
                                                        {
                                                            printf("DLEQ[%d]=%d\n",compareijresult,leqnewstates );
                                                            printf("DTLEQ[%d]=%d\n",compareijresult,tleqnewstates );
                                                            printf("CTLEQ[%d]=%d\n",jp,tlnode );
                                                        }

                                                        unpacked_node* DEQ = unpacked_node::newFull(resF, -rLevel, rSize);
                                                        AddtoNH(DEQ,compareijresult,eqnewstates);
                                                        if(d)
                                                        printf("DEQ[%d]=%d\n",compareijresult,eqnewstates );

                                                        node_handle denode = 0;
                                                        denode=resF->createReducedNode(int(compareijresult), DEQ);//, cev, cnode);
                                                        AddtoNH(CEQ,compareijresult,denode);
                                                        }
                                                }
                                        }else{
                                            if(d)
                                            printf("CAme to else AddtoNH\n" );
                                            if(d)
                                            printf("newstates%d,leqnewstates%d,eqnewstates%d\n",newstates,leqnewstates,eqnewstates );

                                                AddtoNH(D,jp,newstates);
                                                if(d)
                                                printf("D[%d]=%d\n",jp,newstates );

                                                unpacked_node* DEQ = unpacked_node::newFull(resF, -rLevel, rSize);
                                                AddtoNH(DEQ,jp,eqnewstates);
                                                if(d)
                                                printf("DEQ[%d]=%d\n",jp,eqnewstates );

                                                node_handle denode = 0;
                                                denode=resF->createReducedNode(int(jp), DEQ);
                                                AddtoNH(CEQ,jp,denode);
                                                if((above==L &&  compareijresult==-2)||(above==E&& compareijresult==-2 ) ) {
                                                        if(d)
                                                        printf("CAme to else AddtoNH ADDED!\n" );
                                                        AddtoNH(DLEQ,jp,leqnewstates);
                                                        unpacked_node* DTLEQ = unpacked_node::newFull(resF, -rLevel, rSize);
                                                        AddtoNH(DTLEQ,jp,tleqnewstates);
                                                        node_handle tlnode=0;
                                                        tlnode=resF->createReducedNode(int(jp), DTLEQ);
                                                        AddtoNH(CTLEQ,jp,tlnode);
                                                        if(d)
                                                        {
                                                            printf("DLEQ[%d]=%d\n",jp,leqnewstates );
                                                            printf("DTLEQ[%d]=%d\n",jp,tleqnewstates );
                                                            printf("CTLEQ[%d]=%d\n",jp,tlnode );
                                                        }

                                                }
                                        }
                        }
                } // for jz

                unpacked_node::recycle(Rp);
                unpacked_node::recycle(Ru);
        } // else
        node_handle cnode = 0;
        cnode=resF->createReducedNode(int(i), D);
        C->d_ref(i) = cnode;

        node_handle cleqnode = 0;
        cleqnode=resF->createReducedNode(int(i), DLEQ);
        if(d)
        printf("cleqnode %d\n",cleqnode );
        CLEQ->d_ref(i) = cleqnode;

        // ostream_output meddlyout(std::cout);
        unpacked_node::recycle(B);

} //for i
// cleanup mdd reader
unpacked_node::recycle(A);

resEvmdd=resF->createReducedNode(-1, C);
eq=resF->createReducedNode(-1, CEQ);

leq=resF->createReducedNode(-1, CLEQ);
tleq=resF->createReducedNode(-1, CTLEQ);
if(d)
printf("leq %d\n",leq );
#ifdef TRACE_ALL_OPS
printf("computed new tcXrel(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmxd, mxd, resEv, resEvmdd);
#endif
if(evmxd!=-1&& mxd!=-1)
saveResult(Key, evmxd, mxd, resEvmdd,eq,leq,tleq,above,geq);
if(d)
printf("SAVED!evmxd %d, mxd %d, resEvmdd%d,eq%d,leq%d,above%d,geq%d \n", evmxd, mxd, resEvmdd,eq,leq,above,geq);
}
void MEDDLY::mrrc::AddtoNH(unpacked_node* UN, int i, node_handle s ){
    if(0!=s)
    if (0 == UN->d(i)) {
        UN->d_ref(i) = s;

    }else{ // there's new states and existing states; union them.
        dd_edge newstatesE(resF), djp(resF);
        newstatesE.set(s);
        djp.set(UN->d(i));
        accumulateOp->computeTemp(newstatesE, djp, djp);
        UN->set_d(i, djp);
    }
}

void MEDDLY::mrrc::processTerminals( node_handle evmxd, node_handle mxd, node_handle& resEvmxd,node_handle& eq,node_handle& leq,node_handle& tleq)
{
  long evmddval;
  long mxdval;
  long rval;
  argV->getValueFromHandle(evmxd, evmddval);
  argM->getValueFromHandle(mxd, mxdval);
  rval = evmddval * mxdval;
  resEvmxd = resF->handleForValue(rval);
  eq = resF->handleForValue(rval);
  leq = resF->handleForValue(rval);
  tleq=resF->handleForValue(rval);
}

void MEDDLY::mrrc::processTerminals( node_handle evmxd, node_handle mxd, node_handle& resEvmxd,node_handle& eq,node_handle& leq,node_handle& tleq, RV& above)
{
  long evmddval;
  long mxdval;
  long rval;
  argV->getValueFromHandle(evmxd, evmddval);
  argM->getValueFromHandle(mxd, mxdval);
  rval = evmddval * mxdval;
  if(above==N){
      resEvmxd = resF->handleForValue(rval);
      eq = resF->handleForValue(rval);
      leq = 0;
      tleq = 0;
  }else if(above==L){
      resEvmxd =0;
      eq = resF->handleForValue(rval);
      leq = resF->handleForValue(rval);
      tleq = resF->handleForValue(rval);
  }else if(above==E){
      resEvmxd =0;
      eq = resF->handleForValue(rval);
      leq = 0;
      tleq=0;
  }else if(above==G){
      resEvmxd =0;
      eq =0;
      leq = 0;
      tleq=0;
  }
}


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
// *                     preimage_opname  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::preimage_opname : public binary_opname {
  public:
    preimage_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r);
};

MEDDLY::preimage_opname::preimage_opname()
 : binary_opname("Pre-image")
{
}

MEDDLY::binary_operation*
MEDDLY::preimage_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r)
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  binary_opname* accop = nullptr;
  if (a1->getEdgeLabeling() == edge_labeling::EVPLUS || r->getRangeType() == range_type::BOOLEAN) {
    accop = UNION();
  } else {
    accop = MAXIMUM();
  }
  MEDDLY_DCASSERT(accop);

  dd_edge er(r);
  binary_operation* acc = accop->getOperation(er, er, er);

  if (a1->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
    return new mtmatr_mtvect<bool>(this, a1, a2, r, acc);
  }
  else if (a1->getEdgeLabeling() == edge_labeling::EVPLUS) {
    return new mtmatr_evplusvect<int>(this, a1, a2, r, acc);
  }
  else {
    throw error(error::TYPE_MISMATCH);
  }
}


// ******************************************************************
// *                                                                *
// *                     postimage_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::postimage_opname : public binary_opname {
  public:
    postimage_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r);
};

MEDDLY::postimage_opname::postimage_opname()
 : binary_opname("Post-image")
{
}

MEDDLY::binary_operation*
MEDDLY::postimage_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r)
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  binary_opname* accop = nullptr;
  if (a1->getEdgeLabeling() == edge_labeling::EVPLUS || r->getRangeType() == range_type::BOOLEAN) {
    accop = UNION();
  } else {
    accop = MAXIMUM();
  }
  MEDDLY_DCASSERT(accop);

  dd_edge er(r);
  binary_operation* acc = accop->getOperation(er, er, er);

  if (a1->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
    return new mtvect_mtmatr<bool>(this, a1, a2, r, acc);
  }
  else if(a1->getEdgeLabeling() == edge_labeling::EVPLUS) {
      return new evplusvect_mtmatr<int>(this, a1, a2, r, acc);
  }
  else {
    throw error(error::TYPE_MISMATCH);
  }
}

// ******************************************************************
// *                                                                *
// *           transitive_closure_postimage_opname class            *
// *                                                                *
// ******************************************************************

class MEDDLY::transitive_closure_postimage_opname : public binary_opname {
  public:
  transitive_closure_postimage_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r);
};

MEDDLY::transitive_closure_postimage_opname::transitive_closure_postimage_opname()
 : binary_opname("Transitive Closure Post-image")
{
}

MEDDLY::binary_operation*
MEDDLY::transitive_closure_postimage_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r)
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH);

  if (
    !a1->isForRelations()    ||
    !a2->isForRelations()   ||
    !r->isForRelations()     ||
    (a1->getEdgeLabeling() != edge_labeling::EVPLUS) ||
    (a2->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != edge_labeling::EVPLUS)
  )
    throw error(error::TYPE_MISMATCH);

  binary_opname* accop = UNION();
  MEDDLY_DCASSERT(accop);
  dd_edge er(r);
  binary_operation* acc = accop->getOperation(er, er, er);
  return new tcXrel_evplus(this, a1, a2, r, acc);
}

// ******************************************************************
// *                                                                *
// *                   covtc_opname class                           *
// *                                                                *
// ******************************************************************

class MEDDLY::covtc_opname : public binary_opname {
  public:
  covtc_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r);
};

MEDDLY::covtc_opname::covtc_opname()
 : binary_opname("Coverability Transitive Closure")
{
}

MEDDLY::binary_operation*
MEDDLY::covtc_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r)
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH);

  if (
    !a1->isForRelations()    ||
    !a2->isForRelations()   ||
    !r->isForRelations()
  )
    throw error(error::TYPE_MISMATCH);

  binary_opname* accop = UNION();
  MEDDLY_DCASSERT(accop);
  dd_edge er(r);
  binary_operation* acc = accop->getOperation(er, er, er);
  return new covtc(this, a1, a2, r, acc);
}

// ******************************************************************
// *                                                                *
// *    maintain_reachabilityrelation_cover_opname class            *
// *                                                                *
// ******************************************************************

class MEDDLY::maintain_reachabilityrelation_cover_opname : public binary_opname {
  public:
  maintain_reachabilityrelation_cover_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r);
};

MEDDLY::maintain_reachabilityrelation_cover_opname::maintain_reachabilityrelation_cover_opname()
 : binary_opname("Maintain RR and cover")
{
}

MEDDLY::binary_operation*
MEDDLY::maintain_reachabilityrelation_cover_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r)
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH);
  if (
    !a1->isForRelations()    ||
    !a2->isForRelations()   ||
    !r->isForRelations()
    )
    {

        throw error(error::TYPE_MISMATCH);
    }

  binary_opname* accop = UNION();
  MEDDLY_DCASSERT(accop);
  dd_edge er(r);
  binary_operation* acc = accop->getOperation(er,er,er);
  return new mrrc(this, a1, a2, r,r,r,r, acc);
}

// ******************************************************************
// *                                                                *
// *                      VMmult_opname  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::VMmult_opname : public binary_opname {
  public:
    VMmult_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r);
};

MEDDLY::VMmult_opname::VMmult_opname()
 : binary_opname("Vector-matrix multiply")
{
}

MEDDLY::binary_operation*
MEDDLY::VMmult_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r)
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
    (r->getRangeType() == range_type::BOOLEAN) ||
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  binary_opname* accop = PLUS();
  MEDDLY_DCASSERT(accop);
  dd_edge er(r);
  binary_operation* acc = accop->getOperation(er, er, er);

  switch (r->getRangeType()) {
    case range_type::INTEGER:
      return new mtvect_mtmatr<int>(this, a1, a2, r, acc);

    case range_type::REAL:
      return new mtvect_mtmatr<float>(this, a1, a2, r, acc);

    default:
      throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
  }
}

// ******************************************************************
// *                                                                *
// *                      MVmult_opname  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::MVmult_opname : public binary_opname {
  public:
    MVmult_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r);
};

MEDDLY::MVmult_opname::MVmult_opname()
 : binary_opname("Matrix-vector multiply")
{
}

MEDDLY::binary_operation*
MEDDLY::MVmult_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r)
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
    (r->getRangeType() == range_type::BOOLEAN) ||
    !a1->isForRelations()    ||
    a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  binary_opname* accop = PLUS();
  MEDDLY_DCASSERT(accop);
  dd_edge er(r);
  binary_operation* acc = accop->getOperation(er, er, er);

  //
  // We're switching the order of the arguments
  //

  switch (r->getRangeType()) {
    case range_type::INTEGER:
      return new mtmatr_mtvect<int>(this, a2, a1, r, acc);

    case range_type::REAL:
      return new mtmatr_mtvect<float>(this, a2, a1, r, acc);

    default:
      throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
  }
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializePreImage()
{
  return new preimage_opname;
}

MEDDLY::binary_opname* MEDDLY::initializePostImage()
{
  return new postimage_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeTCPostImage()
{
  return new transitive_closure_postimage_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeCOVTC()
{
  return new covtc_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeMRCPostImage()
{
  return new maintain_reachabilityrelation_cover_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeVMmult()
{
  return new VMmult_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeMVmult()
{
  return new MVmult_opname;
}

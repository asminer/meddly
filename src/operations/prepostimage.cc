
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
#include "old_meddly.h"
#include "old_meddly.hh"
#include "old_meddly_expert.h"
#include "old_meddly_expert.hh"
#include "prepostimage.h"

// #define TRACE_ALL_OPS

namespace MEDDLY {
  class image_op;
  class relXset_mdd;
  class setXrel_mdd;

  class image_op_evplus;
  class relXset_evplus;
  class setXrel_evplus;
  class tcXrel_evplus;

  class preimage_opname;
  class postimage_opname;

  class transitive_closure_postimage_opname;

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
    image_op(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

    inline compute_table::entry_key*
    findResult(node_handle a, node_handle b, node_handle &c)
    {
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CTsrch->writeN(b);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(compute_table::entry_key* Key,
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

    expert_forest* argV;
    expert_forest* argM;
};

MEDDLY::image_op::image_op(const binary_opname* oc, expert_forest* a1,
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

  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "NN:N");
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

MEDDLY::node_handle MEDDLY::image_op::compute(node_handle a, node_handle b)
{
  MEDDLY_DCASSERT(accumulateOp);
  return compute_rec(a, b);
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
    relXset_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

  protected:
    virtual node_handle compute_rec(node_handle a, node_handle b);
    virtual node_handle processTerminals(node_handle mdd, node_handle mxd) = 0;
};

MEDDLY::relXset_mdd::relXset_mdd(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res, binary_operation* acc)
: image_op(oc, a1, a2, res, acc)
{
}

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
  compute_table::entry_key* Key = findResult(mdd, mxd, result);
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
    setXrel_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

  protected:
    virtual node_handle compute_rec(node_handle a, node_handle b);
    virtual node_handle processTerminals(node_handle mdd, node_handle mxd) = 0;
};

MEDDLY::setXrel_mdd::setXrel_mdd(const binary_opname* oc,
  expert_forest* a1, expert_forest* a2, expert_forest* res, binary_operation* acc)
: image_op(oc, a1, a2, res, acc)
{
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
  compute_table::entry_key* Key = findResult(mdd, mxd, result);
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
      mtmatr_mtvect(const binary_opname* opcode, expert_forest* arg1,
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
      mtmatr_mtvect(const binary_opname* opcode, expert_forest* arg1,
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
      mtvect_mtmatr(const binary_opname* opcode, expert_forest* arg1,
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
      mtvect_mtmatr(const binary_opname* opcode, expert_forest* arg1,
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
    image_op_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

    inline compute_table::entry_key*
    findResult(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle &resEvmdd)
    {
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
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
    inline void saveResult(compute_table::entry_key* Key,
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


MEDDLY::image_op_evplus::image_op_evplus(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res, binary_operation* acc)
: binary_operation(oc, 1, a1, a2, res)
{
  accumulateOp = acc;

  argV = a1;
  argM = a2;

  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "NN:LN");
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
// *                     relXset_evplus  class                      *
// *                                                                *
// ******************************************************************

/** Generic base for relation multiplied by set.
    Changing what happens at the terminals can give
    different meanings to this operation :^)
*/
class MEDDLY::relXset_evplus : public image_op_evplus {
  public:
    relXset_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

  protected:
    virtual void compute_rec(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
    virtual void processTerminals(long ev, node_handle mdd, node_handle mxd, long& resEv, node_handle& resEvmdd) = 0;
};

MEDDLY::relXset_evplus::relXset_evplus(const binary_opname* oc,
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
  compute_table::entry_key* Key = findResult(ev, evmdd, mxd, resEv, resEvmdd);
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
    setXrel_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res, binary_operation* acc);

  protected:
    virtual void compute_rec(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
    virtual void processTerminals(long ev, node_handle mdd, node_handle mxd, long& resEv, node_handle& resEvmdd) = 0;
};

MEDDLY::setXrel_evplus::setXrel_evplus(const binary_opname* oc,
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
  compute_table::entry_key* Key = findResult(ev, evmdd, mxd, resEv, resEvmdd);
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
    mtmatr_evplusvect(const binary_opname* opcode, expert_forest* arg1,
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
    evplusvect_mtmatr(const binary_opname* opcode, expert_forest* arg1,
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
    tcXrel_evplus(const binary_opname* opcode, expert_forest* tc,
      expert_forest* trans, expert_forest* res, binary_operation* acc);

  protected:
    virtual void compute_rec(long ev, node_handle evmxd, node_handle mxd, long& resEv, node_handle& resEvmxd);
    virtual void processTerminals(long ev, node_handle evmxd, node_handle mxd, long& resEv, node_handle& resEvmxd);
};

MEDDLY::tcXrel_evplus::tcXrel_evplus(const binary_opname* oc,
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
  compute_table::entry_key* Key = findResult(ev, evmxd, mxd, resEv, resEvmdd);
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
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::preimage_opname::preimage_opname()
 : binary_opname("Pre-image")
{
}

MEDDLY::binary_operation*
MEDDLY::preimage_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r) const
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
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  binary_operation* acc = 0;
  if (a1->getEdgeLabeling() == forest::EVPLUS || r->getRangeType() == forest::BOOLEAN) {
    acc = getOperation(UNION, r, r, r);
  } else {
    acc = getOperation(MAXIMUM, r, r, r);
  }

  if (a1->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    return new mtmatr_mtvect<bool>(this, a1, a2, r, acc);
  }
  else if (a1->getEdgeLabeling() == forest::EVPLUS) {
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
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::postimage_opname::postimage_opname()
 : binary_opname("Post-image")
{
}

MEDDLY::binary_operation*
MEDDLY::postimage_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r) const
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
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  binary_operation* acc = 0;
  if (a1->getEdgeLabeling() == forest::EVPLUS || r->getRangeType() == forest::BOOLEAN) {
    acc = getOperation(UNION, r, r, r);
  } else {
    acc = getOperation(MAXIMUM, r, r, r);
  }

  if (a1->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    return new mtvect_mtmatr<bool>(this, a1, a2, r, acc);
  }
  else if(a1->getEdgeLabeling() == forest::EVPLUS) {
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
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::transitive_closure_postimage_opname::transitive_closure_postimage_opname()
 : binary_opname("Transitive Closure Post-image")
{
}

MEDDLY::binary_operation*
MEDDLY::transitive_closure_postimage_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r) const
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
    (a1->getEdgeLabeling() != forest::EVPLUS) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != forest::EVPLUS)
  )
    throw error(error::TYPE_MISMATCH);

  binary_operation* acc = getOperation(UNION, r, r, r);
  return new tcXrel_evplus(this, a1, a2, r, acc);
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
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::VMmult_opname::VMmult_opname()
 : binary_opname("Vector-matrix multiply")
{
}

MEDDLY::binary_operation*
MEDDLY::VMmult_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    (a1->getRangeType() == forest::BOOLEAN) ||
    (a2->getRangeType() == forest::BOOLEAN) ||
    (r->getRangeType() == forest::BOOLEAN) ||
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  binary_operation* acc = getOperation(PLUS, r, r, r);

  switch (r->getRangeType()) {
    case forest::INTEGER:
      return new mtvect_mtmatr<int>(this, a1, a2, r, acc);

    case forest::REAL:
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
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::MVmult_opname::MVmult_opname()
 : binary_opname("Matrix-vector multiply")
{
}

MEDDLY::binary_operation*
MEDDLY::MVmult_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    (a1->getRangeType() == forest::BOOLEAN) ||
    (a2->getRangeType() == forest::BOOLEAN) ||
    (r->getRangeType() == forest::BOOLEAN) ||
    !a1->isForRelations()    ||
    a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  binary_operation* acc = getOperation(PLUS, r, r, r);

  //
  // We're switching the order of the arguments
  //

  switch (r->getRangeType()) {
    case forest::INTEGER:
      return new mtmatr_mtvect<int>(this, a2, a1, r, acc);

    case forest::REAL:
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

MEDDLY::binary_opname* MEDDLY::initializeVMmult()
{
  return new VMmult_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeMVmult()
{
  return new MVmult_opname;
}


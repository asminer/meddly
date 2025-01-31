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

#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"

// #define TRACE_ALL_OPS

namespace MEDDLY {
    class image_op;
    class relXset_mdd;
    class setXrel_mdd;

    class image_op_evplus;
    class relXset_evplus;
    class setXrel_evplus;
    class tcXrel_evplus;

    binary_list PRE_IMAGE_cache;
    binary_list POST_IMAGE_cache;

    binary_list TC_POST_IMAGE_cache;
    binary_list VM_MULTIPLY_cache;
    binary_list MV_MULTIPLY_cache;
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
    image_op(binary_list& opcode, forest* arg1,
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

    forest* argV;
    forest* argM;
};

MEDDLY::image_op::image_op(binary_list& oc, forest* a1,
  forest* a2, forest* res, binary_operation* acc)
: binary_operation(oc, 1, a1, a2, res)
{
  accumulateOp = acc;

  checkDomains(__FILE__, __LINE__);
  checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

  if (a1->isForRelations()) {
    argM = a1;
    argV = a2;
    if (a2->isForRelations()) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  } else {
    argM = a2;
    argV = a1;
    if (!a2->isForRelations()) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  }

  ct_entry_type* et = new ct_entry_type(oc.getName(), "NN:N");
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
    relXset_mdd(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res, binary_operation* acc);

  protected:
    virtual node_handle compute_rec(node_handle a, node_handle b);
    virtual node_handle processTerminals(node_handle mdd, node_handle mxd) = 0;
};

MEDDLY::relXset_mdd::relXset_mdd(binary_list& oc, forest* a1,
  forest* a2, forest* res, binary_operation* acc)
: image_op(oc, a1, a2, res, acc)
{
    checkRelations(__FILE__, __LINE__, SET, RELATION, SET);
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
  ct_entry_key* Key = findResult(mdd, mxd, result);
  if (0==Key) return result;

  // check if mxd and mdd are at the same level
  const int mddLevel = argV->getNodeLevel(mdd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(argV);
  if (mddLevel < rLevel) {
    A->initRedundant(argV, rLevel, mdd, FULL_ONLY);
  } else {
    argV->unpackNode(A, mdd, FULL_ONLY);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      C->setFull(i, compute_rec(A->down(i), mxd));
      // C->d_ref(i) = compute_rec(A->down(i), mxd);
    }
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    C->clear(0, rSize);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New(argM);
    unpacked_node *Rp = unpacked_node::New(argM);
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, SPARSE_ONLY);
    } else {
      argM->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    dd_edge newstatesE(resF), cdi(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      unsigned i = Ru->index(iz);
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->down(iz), SPARSE_ONLY);
      } else {
        argM->unpackNode(Rp, Ru->down(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        unsigned j = Rp->index(jz);
        MEDDLY_DCASSERT(0<=j && j < A->getSize());
        if (0==A->down(j))   continue;
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = compute_rec(A->down(j), Rp->down(jz));
        if (0==newstates) continue;
        if (0==C->down(i)) {
          C->setFull(i, newstates);
          continue;
        }
        // there's new states and existing states; union them.
        newstatesE.set(newstates);
        cdi.set(C->down(i));
        accumulateOp->computeTemp(newstatesE, cdi, cdi);
        C->setFull(i, cdi);
      } // for j

    } // for i

    unpacked_node::Recycle(Rp);
    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);

  edge_value ev;
  resF->createReducedNode(C, ev, result);
  MEDDLY_DCASSERT(ev.isVoid());
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
    setXrel_mdd(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res, binary_operation* acc);

  protected:
    virtual node_handle compute_rec(node_handle a, node_handle b);
    virtual node_handle processTerminals(node_handle mdd, node_handle mxd) = 0;
};

MEDDLY::setXrel_mdd::setXrel_mdd(binary_list& oc,
  forest* a1, forest* a2, forest* res, binary_operation* acc)
: image_op(oc, a1, a2, res, acc)
{
    checkRelations(__FILE__, __LINE__, SET, RELATION, SET);
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
  unpacked_node *A = unpacked_node::New(argV);
  if (mddLevel < rLevel) {
    A->initRedundant(argV, rLevel, mdd, FULL_ONLY);
  } else {
    argV->unpackNode(A, mdd, FULL_ONLY);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      C->setFull(i, compute_rec(A->down(i), mxd));
      // C->d_ref(i) = compute_rec(A->down(i), mxd);
    }
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    C->clear(0, rSize);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New(argM);
    unpacked_node *Rp = unpacked_node::New(argM);
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, SPARSE_ONLY);
    } else {
      argM->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    dd_edge newstatesE(resF), cdj(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      unsigned i = Ru->index(iz);
      if (0==A->down(i))   continue;
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->down(iz), SPARSE_ONLY);
      } else {
        argM->unpackNode(Rp, Ru->down(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        unsigned j = Rp->index(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = compute_rec(A->down(i), Rp->down(jz));
        if (0==newstates) continue;
        if (0==C->down(j)) {
          C->setFull(j, newstates);
          continue;
        }
        // there's new states and existing states; union them.
        newstatesE.set(newstates);
        cdj.set(C->down(j));
        accumulateOp->computeTemp(newstatesE, cdj, cdj);
        C->setFull(j, cdj);
      } // for j

    } // for i

    unpacked_node::Recycle(Rp);
    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);

  edge_value ev;
  resF->createReducedNode(C, ev, result);
  MEDDLY_DCASSERT(ev.isVoid());
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
      mtmatr_mtvect(binary_list& opcode, forest* arg1,
        forest* arg2, forest* res, binary_operation* acc)
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
      mtmatr_mtvect(binary_list& opcode, forest* arg1,
        forest* arg2, forest* res, binary_operation* acc)
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
      mtvect_mtmatr(binary_list& opcode, forest* arg1,
        forest* arg2, forest* res, binary_operation* acc)
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
      mtvect_mtmatr(binary_list& opcode, forest* arg1,
        forest* arg2, forest* res, binary_operation* acc)
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
    image_op_evplus(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res, binary_operation* acc);

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

    forest* argV;
    forest* argM;
};


MEDDLY::image_op_evplus::image_op_evplus(binary_list& oc, forest* a1,
  forest* a2, forest* res, binary_operation* acc)
: binary_operation(oc, 1, a1, a2, res)
{
  accumulateOp = acc;

    checkDomains(__FILE__, __LINE__);
    checkLabelings(__FILE__, __LINE__,
        edge_labeling::EVPLUS,
        edge_labeling::MULTI_TERMINAL,
        edge_labeling::EVPLUS
    );

  argV = a1;
  argM = a2;

  ct_entry_type* et = new ct_entry_type(oc.getName(), "NN:LN");
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
  c.set(cev, cnode);
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
    relXset_evplus(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res, binary_operation* acc);

  protected:
    virtual void compute_rec(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
    virtual void processTerminals(long ev, node_handle mdd, node_handle mxd, long& resEv, node_handle& resEvmdd) = 0;
};

MEDDLY::relXset_evplus::relXset_evplus(binary_list& oc,
  forest* a1, forest* a2, forest* res, binary_operation* acc)
: image_op_evplus(oc, a1, a2, res, acc)
{
    checkRelations(__FILE__, __LINE__, RELATION, SET, SET);
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
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmdd, FULL_ONLY)
    : argV->newUnpacked(evmdd, FULL_ONLY);

  if (evmddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      long nev = Inf<long>();
      node_handle newstates = 0;
      compute_rec(A->edge_long(i), A->down(i), mxd, nev, newstates);

      C->setFull(i, edge_value( newstates ? ev + nev : 0L ), newstates);
      // C->setEdge(i, newstates == 0 ? 0L : ev + nev);
      // C->d_ref(i) = newstates;
    }
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= evmddLevel);

    // clear out result (important!)
    C->clear(0, rSize);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New(argM);
    unpacked_node *Rp = unpacked_node::New(argM);
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, SPARSE_ONLY);
    } else {
      argM->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    dd_edge newstatesE(resF), cdi(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      unsigned i = Ru->index(iz);
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->down(iz), SPARSE_ONLY);
      } else {
        argM->unpackNode(Rp, Ru->down(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        unsigned j = Rp->index(jz);
        if (0==A->down(j))   continue;
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        long nev = Inf<long>();
        node_handle newstates = 0;
        compute_rec(A->edge_long(j), A->down(j), Rp->down(jz), nev, newstates);
        if (0==newstates) continue;
        nev += ev;
        if (0==C->down(i)) {
          C->setFull(i, edge_value(nev), newstates);
          // C->setEdge(i, nev);
          // C->d_ref(i) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        newstatesE.set(nev, newstates);
        cdi.set(C->edge_long(i), C->down(i));
        accumulateOp->computeTemp(newstatesE, cdi, cdi);
        C->setFull(i, cdi);
      } // for j

    } // for i

    unpacked_node::Recycle(Rp);
    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);

  edge_value Cev;
  resF->createReducedNode(C, Cev, resEvmdd);
  resEv = Cev.getLong();
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
    setXrel_evplus(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res, binary_operation* acc);

  protected:
    virtual void compute_rec(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
    virtual void processTerminals(long ev, node_handle mdd, node_handle mxd, long& resEv, node_handle& resEvmdd) = 0;
};

MEDDLY::setXrel_evplus::setXrel_evplus(binary_list& oc,
  forest* a1, forest* a2, forest* res, binary_operation* acc)
: image_op_evplus(oc, a1, a2, res, acc)
{
    checkRelations(__FILE__, __LINE__, SET, RELATION, SET);
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
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmdd, FULL_ONLY)
    : argV->newUnpacked(evmdd, FULL_ONLY);

  if (evmddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      long nev = Inf<long>();
      node_handle newstates = 0;
      compute_rec(A->edge_long(i), A->down(i), mxd, nev, newstates);

      C->setFull(i, edge_value( newstates ? ev + nev : 0L), newstates);
      // C->setEdge(i, newstates == 0 ? 0L : ev + nev);
      // C->d_ref(i) = newstates;
    }
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= evmddLevel);

    // clear out result (important!)
    C->clear(0, rSize);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New(argM);
    unpacked_node *Rp = unpacked_node::New(argM);
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, SPARSE_ONLY);
    } else {
      argM->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    dd_edge newstatesE(resF), cdj(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      unsigned i = Ru->index(iz);
      if (0==A->down(i))   continue;
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->down(iz), SPARSE_ONLY);
      } else {
        argM->unpackNode(Rp, Ru->down(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        unsigned j = Rp->index(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        long nev = Inf<long>();
        node_handle newstates = 0;
        compute_rec(A->edge_long(i), A->down(i), Rp->down(jz), nev, newstates);
        if (0==newstates) continue;
        nev += ev;
        if (0==C->down(j)) {
          C->setFull(j, nev, newstates);
          // C->setEdge(j, nev);
          // C->d_ref(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        newstatesE.set(nev, newstates);
        cdj.set(C->edge_long(j), C->down(j));
        accumulateOp->computeTemp(newstatesE, cdj, cdj);
        C->setFull(j, cdj);
      } // for j

    } // for i

    unpacked_node::Recycle(Rp);
    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);

  edge_value Cev;
  resF->createReducedNode(C, Cev, resEvmdd);
  resEv = Cev.getLong();
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
    mtmatr_evplusvect(binary_list& opcode, forest* arg1,
        forest* arg2, forest* res, binary_operation* acc)
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
    evplusvect_mtmatr(binary_list& opcode, forest* arg1,
        forest* arg2, forest* res, binary_operation* acc)
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
    tcXrel_evplus(binary_list& opcode, forest* tc,
      forest* trans, forest* res, binary_operation* acc);

  protected:
    virtual void compute_rec(long ev, node_handle evmxd, node_handle mxd, long& resEv, node_handle& resEvmxd);
    virtual void processTerminals(long ev, node_handle evmxd, node_handle mxd, long& resEv, node_handle& resEvmxd);
};

MEDDLY::tcXrel_evplus::tcXrel_evplus(binary_list& oc,
  forest* tc, forest* trans, forest* res, binary_operation* acc)
: image_op_evplus(oc, tc, trans, res, acc)
{
    checkRelations(__FILE__, __LINE__, RELATION, RELATION, RELATION);
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
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmxd, FULL_ONLY)
    : argV->newUnpacked(evmxd, FULL_ONLY);

  for (unsigned i = 0; i < rSize; i++) {
    int pLevel = argV->getNodeLevel(A->down(i));
    unpacked_node* B = isLevelAbove(-rLevel, pLevel)
      ? unpacked_node::newIdentity(argV, -rLevel, i, 0L, A->down(i), FULL_ONLY)
      : argV->newUnpacked(A->down(i), FULL_ONLY);

    unpacked_node* D = unpacked_node::newFull(resF, -rLevel, rSize);
    if (rLevel > ABS(mxdLevel)) {
      //
      // Skipped levels in the MXD,
      // that's an important special case that we can handle quickly.
      for (unsigned j = 0; j < rSize; j++) {
        long nev = Inf<long>();
        node_handle newstates = 0;
        compute_rec(A->edge_long(i) + B->edge_long(j), B->down(j), mxd, nev, newstates);

        D->setFull(j, edge_value(newstates ? ev + nev : 0L), newstates);
        // D->setEdge(j, newstates == 0 ? 0L : ev + nev);
        // D->d_ref(j) = newstates;
      }
    }
    else {
      //
      // Need to process this level in the MXD.
      MEDDLY_DCASSERT(ABS(mxdLevel) >= ABS(pLevel));

      // clear out result (important!)
      D->clear(0, rSize);

      // Initialize mxd readers, note we might skip the unprimed level
      unpacked_node *Ru = unpacked_node::New(argM);
      unpacked_node *Rp = unpacked_node::New(argM);
      if (mxdLevel < 0) {
        Ru->initRedundant(argM, rLevel, mxd, SPARSE_ONLY);
      } else {
        argM->unpackNode(Ru, mxd, SPARSE_ONLY);
      }

      dd_edge newstatesE(resF), djp(resF);

      // loop over mxd "rows"
      for (unsigned jz = 0; jz < Ru->getSize(); jz++) {
        unsigned j = Ru->index(jz);
        if (0 == B->down(j)) {
          continue;
        }

        if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->down(jz)))) {
          Rp->initIdentity(argM, rLevel, j, Ru->down(jz), SPARSE_ONLY);
        } else {
          argM->unpackNode(Rp, Ru->down(jz), SPARSE_ONLY);
        }

        // loop over mxd "columns"
        for (unsigned jpz = 0; jpz < Rp->getSize(); jpz++) {
          unsigned jp = Rp->index(jpz);
          // ok, there is an i->j "edge".
          // determine new states to be added (recursively)
          // and add them
          long nev = Inf<long>();
          node_handle newstates = 0;
          compute_rec(A->edge_long(i) + B->edge_long(j), B->down(j), Rp->down(jpz), nev, newstates);
          if (0==newstates) {
            continue;
          }
          nev += ev;
          if (0 == D->down(jp)) {
            D->setFull(jp, edge_value(nev), newstates);
            // D->setEdge(jp, nev);
            // D->d_ref(jp) = newstates;
            continue;
          }
          // there's new states and existing states; union them.
          newstatesE.set(nev, newstates);
          djp.set(D->edge_long(jp), D->down(jp));
          accumulateOp->computeTemp(newstatesE, djp, djp);
          D->setFull(jp, djp);
        } // for j

      } // for i

      unpacked_node::Recycle(Rp);
      unpacked_node::Recycle(Ru);
    } // else

    edge_value cev;
    node_handle cnode;
    resF->createReducedNode(D, cev, cnode, int(i));
    C->setFull(i, cev, cnode);
    // C->setEdge(i, cev);
    // C->d_ref(i) = cnode;

    unpacked_node::Recycle(B);
  }

  // cleanup mdd reader
  unpacked_node::Recycle(A);

  edge_value Cev;
  resF->createReducedNode(C, Cev, resEvmdd);
  resEv = Cev.getLong();
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
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::PRE_IMAGE(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  PRE_IMAGE_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    binary_operation *acc = nullptr;
    if  (
            a->getEdgeLabeling() == edge_labeling::EVPLUS ||
            c->getRangeType() == range_type::BOOLEAN
        )
    {
        acc = UNION(c, c, c);
    } else {
        acc = MAXIMUM(c, c, c);
    }
    MEDDLY_DCASSERT(acc);

    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        return PRE_IMAGE_cache.add(
            new mtmatr_mtvect<bool>(PRE_IMAGE_cache, a, b, c, acc)
        );
    }
    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {
        return PRE_IMAGE_cache.add(
            new mtmatr_evplusvect<int>(PRE_IMAGE_cache, a, b, c, acc)
        );
    }
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::PRE_IMAGE_init()
{
    PRE_IMAGE_cache.reset("Pre-image");
}

void MEDDLY::PRE_IMAGE_done()
{
    MEDDLY_DCASSERT(PRE_IMAGE_cache.isEmpty());
}

// ******************************************************************

MEDDLY::binary_operation* MEDDLY::POST_IMAGE(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  POST_IMAGE_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    binary_operation *acc = nullptr;
    if  (
            a->getEdgeLabeling() == edge_labeling::EVPLUS ||
            c->getRangeType() == range_type::BOOLEAN
        )
    {
        acc = UNION(c, c, c);
    } else {
        acc = MAXIMUM(c, c, c);
    }
    MEDDLY_DCASSERT(acc);

    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        return POST_IMAGE_cache.add(
            new mtvect_mtmatr<bool>(POST_IMAGE_cache, a, b, c, acc)
        );
    }
    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {
        return POST_IMAGE_cache.add(
            new evplusvect_mtmatr<int>(POST_IMAGE_cache, a, b, c, acc)
        );
    }
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::POST_IMAGE_init()
{
    POST_IMAGE_cache.reset("Post-image");
}

void MEDDLY::POST_IMAGE_done()
{
    MEDDLY_DCASSERT(POST_IMAGE_cache.isEmpty());
}

// ******************************************************************

MEDDLY::binary_operation* MEDDLY::TC_POST_IMAGE(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  TC_POST_IMAGE_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    return TC_POST_IMAGE_cache.add(
        new tcXrel_evplus(TC_POST_IMAGE_cache, a, b, c, UNION(c, c, c))
    );
}

void MEDDLY::TC_POST_IMAGE_init()
{
    TC_POST_IMAGE_cache.reset("TC-Post-image");
}

void MEDDLY::TC_POST_IMAGE_done()
{
    MEDDLY_DCASSERT(TC_POST_IMAGE_cache.isEmpty());
}

// ******************************************************************

MEDDLY::binary_operation* MEDDLY::VM_MULTIPLY(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  VM_MULTIPLY_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if  (
            (a->getRangeType() == range_type::BOOLEAN) ||
            (b->getRangeType() == range_type::BOOLEAN) ||
            (c->getRangeType() == range_type::BOOLEAN)
        )
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    binary_operation* acc = PLUS(c, c, c);

    switch (c->getRangeType()) {
        case range_type::INTEGER:
            return VM_MULTIPLY_cache.add(
                new mtvect_mtmatr<int>(VM_MULTIPLY_cache, a, b, c, acc)
            );

        case range_type::REAL:
            return VM_MULTIPLY_cache.add(
                new mtvect_mtmatr<float>(VM_MULTIPLY_cache, a, b, c, acc)
            );

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
}

void MEDDLY::VM_MULTIPLY_init()
{
    VM_MULTIPLY_cache.reset("VM-multiply");
}

void MEDDLY::VM_MULTIPLY_done()
{
    MEDDLY_DCASSERT(VM_MULTIPLY_cache.isEmpty());
}

// ******************************************************************

MEDDLY::binary_operation* MEDDLY::MV_MULTIPLY(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  MV_MULTIPLY_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if  (
            (a->getRangeType() == range_type::BOOLEAN) ||
            (b->getRangeType() == range_type::BOOLEAN) ||
            (c->getRangeType() == range_type::BOOLEAN)
        )
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    binary_operation* acc = PLUS(c, c, c);

    //
    // We're switching the order of the arguments
    //

    switch (c->getRangeType()) {
        case range_type::INTEGER:
            return MV_MULTIPLY_cache.add(
                new mtmatr_mtvect<int>(MV_MULTIPLY_cache, b, a, c, acc)
            );

        case range_type::REAL:
            return MV_MULTIPLY_cache.add(
                new mtmatr_mtvect<float>(MV_MULTIPLY_cache, b, a, c, acc)
            );

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
}

void MEDDLY::MV_MULTIPLY_init()
{
    MV_MULTIPLY_cache.reset("MV-multiply");
}

void MEDDLY::MV_MULTIPLY_done()
{
    MEDDLY_DCASSERT(MV_MULTIPLY_cache.isEmpty());
}


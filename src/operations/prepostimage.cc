
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

#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
#endif
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData) const;

    inline compute_table::search_key* 
    findResult(node_handle a, node_handle b, node_handle &c) 
    {
      compute_table::search_key* CTsrch = useCTkey();
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->reset();
      CTsrch->writeNH(a);
      CTsrch->writeNH(b);
      compute_table::search_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      c = resF->linkNode(cacheFind.readNH());
      doneCTkey(CTsrch);
      return 0;
    }
    inline node_handle saveResult(compute_table::search_key* Key, 
      node_handle a, node_handle b, node_handle c) 
    {
      argV->cacheNode(a);
      argM->cacheNode(b);
      compute_table::entry_builder &entry = CT->startNewEntry(Key);
      entry.writeResultNH(resF->cacheNode(c));
      CT->addEntry();
      return c;
    }
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual node_handle compute(node_handle a, node_handle b);
  protected:
    binary_operation* accumulateOp;
    virtual node_handle compute_rec(node_handle a, node_handle b) = 0;

    expert_forest* argV;
    expert_forest* argM;
};

MEDDLY::image_op::image_op(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res, binary_operation* acc)
: binary_operation(oc, 2, 1, a1, a2, res)
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
}

#ifndef USE_NODE_STATUS
bool MEDDLY::image_op::isStaleEntry(const node_handle* data)
{
  return argV->isStale(data[0]) ||
         argM->isStale(data[1]) ||
         resF->isStale(data[2]);
}
#else
MEDDLY::forest::node_status
MEDDLY::image_op::getStatusOfEntry(const node_handle* data)
{
  MEDDLY::forest::node_status a = argV->getNodeStatus(data[0]);
  MEDDLY::forest::node_status b = argM->getNodeStatus(data[1]);
  MEDDLY::forest::node_status c = resF->getNodeStatus(data[2]);

  if (a == MEDDLY::forest::DEAD ||
      b == MEDDLY::forest::DEAD ||
      c == MEDDLY::forest::DEAD)
    return MEDDLY::forest::DEAD;
  else if (a == MEDDLY::forest::RECOVERABLE ||
      b == MEDDLY::forest::RECOVERABLE ||
      c == MEDDLY::forest::RECOVERABLE)
    return MEDDLY::forest::RECOVERABLE;
  else
    return MEDDLY::forest::ACTIVE;
}
#endif

void MEDDLY::image_op::discardEntry(const node_handle* data)
{
  argV->uncacheNode(data[0]);
  argM->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void
MEDDLY::image_op::showEntry(output &strm, const node_handle* data) const
{
  strm  << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1])
        << "): " << long(data[2]) << "]";
}

void MEDDLY::image_op
::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  node_handle cnode;
  if (a.getForest() == argV) {
    cnode = compute(a.getNode(), b.getNode());
  } else {
    cnode = compute(b.getNode(), a.getNode());
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
  compute_table::search_key* Key = findResult(mdd, mxd, result);
  if (0==Key) return result;

  // check if mxd and mdd are at the same level
  const int mddLevel = argV->getNodeLevel(mdd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const int rSize = resF->getLevelSize(rLevel);
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::useUnpackedNode();
  if (mddLevel < rLevel) {
    A->initRedundant(argV, rLevel, mdd, true);
  } else {
    A->initFromNode(argV, mdd, true);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (int i=0; i<rSize; i++) {
      C->d_ref(i) = compute_rec(A->d(i), mxd);
    }
  } else {
    // 
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (int i=0; i<rSize; i++) C->d_ref(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::useUnpackedNode();
    unpacked_node *Rp = unpacked_node::useUnpackedNode();
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, false);
    } else {
      Ru->initFromNode(argM, mxd, false);
    }

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->d(iz), false);
      } else {
        Rp->initFromNode(argM, Ru->d(iz), false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
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
        node_handle oldi = C->d(i);
        C->d_ref(i) = accumulateOp->compute(newstates, oldi);
        resF->unlinkNode(oldi);
        resF->unlinkNode(newstates);
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
  compute_table::search_key* Key = findResult(mdd, mxd, result);
  if (0==Key) return result;

  // check if mxd and mdd are at the same level
  const int mddLevel = argV->getNodeLevel(mdd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const int rSize = resF->getLevelSize(rLevel);
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::useUnpackedNode();
  if (mddLevel < rLevel) {
    A->initRedundant(argV, rLevel, mdd, true);
  } else {
    A->initFromNode(argV, mdd, true);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (int i=0; i<rSize; i++) {
      C->d_ref(i) = compute_rec(A->d(i), mxd);
    }
  } else {
    // 
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (int i=0; i<rSize; i++) C->d_ref(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::useUnpackedNode();
    unpacked_node *Rp = unpacked_node::useUnpackedNode();
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, false);
    } else {
      Ru->initFromNode(argM, mxd, false);
    }

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      if (0==A->d(i))   continue; 
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->d(iz), false);
      } else {
        Rp->initFromNode(argM, Ru->d(iz), false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
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
        node_handle oldj = C->d(j);
        C->d_ref(j) = accumulateOp->compute(newstates, oldj);
        resF->unlinkNode(oldj);
        resF->unlinkNode(newstates);
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

    virtual bool isStaleEntry(const node_handle* entryData);
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData) const;

    inline compute_table::search_key*
    findResult(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle &resEvmdd)
    {
      compute_table::search_key* CTsrch = useCTkey();
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->reset();
      CTsrch->writeNH(evmdd);
      CTsrch->writeNH(mxd);
      compute_table::search_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      cacheFind.read(resEv);
      resEvmdd = resF->linkNode(cacheFind.readNH());
      if (resEvmdd != 0) {
        resEv += ev;
      }
      doneCTkey(CTsrch);
      return 0;
    }
    inline void saveResult(compute_table::search_key* Key,
      long ev, node_handle evmdd, node_handle mxd, long resEv, node_handle resEvmdd)
    {
      argV->cacheNode(evmdd);
      argM->cacheNode(mxd);
      compute_table::entry_builder &entry = CT->startNewEntry(Key);
      entry.writeResult(resEvmdd == 0 ? 0L : resEv - ev);
      entry.writeResultNH(resF->cacheNode(resEvmdd));
      CT->addEntry();
    }
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual void compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
  protected:
    binary_operation* accumulateOp;
    virtual void compute_rec(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd) = 0;

    expert_forest* argV;
    expert_forest* argM;
};

MEDDLY::image_op_evplus::image_op_evplus(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res, binary_operation* acc)
: binary_operation(oc,
    (sizeof(node_handle) + sizeof(node_handle)) / sizeof(node_handle),
    (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle),
    a1, a2, res)
{
  accumulateOp = acc;

  argV = a1;
  argM = a2;

//  if (a1->isForRelations()) {
//    argM = a1;
//    argV = a2;
//    if (a2->isForRelations()) throw error(error::MISCELLANEOUS);
//  } else {
//    argM = a2;
//    argV = a1;
//    if (!a2->isForRelations()) throw error(error::MISCELLANEOUS);
//  }
}

bool MEDDLY::image_op_evplus::isStaleEntry(const node_handle* data)
{
  return argV->isStale(data[0]) ||
         argM->isStale(data[sizeof(node_handle) / sizeof(node_handle)]) ||
         resF->isStale(data[(2 * sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)]);
}

void MEDDLY::image_op_evplus::discardEntry(const node_handle* data)
{
  argV->uncacheNode(data[0]);
  argM->uncacheNode(data[sizeof(node_handle) / sizeof(node_handle)]);
  resF->uncacheNode(data[(2 * sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)]);
}

void MEDDLY::image_op_evplus::showEntry(output &strm, const node_handle* data) const
{
  strm  << "[" << getName()
        << "(" << long(data[0])
        << ", " << long(data[sizeof(node_handle) / sizeof(node_handle)])
        << "): " << long(data[(2 * sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)])
        << "]";
}

void MEDDLY::image_op_evplus::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c)
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
  compute_table::search_key* Key = findResult(ev, evmdd, mxd, resEv, resEvmdd);
  if (0==Key) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int evmddLevel = argV->getNodeLevel(evmdd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), evmddLevel);
  const int rSize = resF->getLevelSize(rLevel);
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize evmdd reader
  unpacked_node *A = (evmddLevel < rLevel)
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmdd, true)
    : unpacked_node::newFromNode(argV, evmdd, true);

  if (evmddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (int i=0; i<rSize; i++) {
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
    for (int i=0; i<rSize; i++) {
      C->setEdge(i, 0L);
      C->d_ref(i) = 0;
    }

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::useUnpackedNode();
    unpacked_node *Rp = unpacked_node::useUnpackedNode();
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, false);
    } else {
      Ru->initFromNode(argM, mxd, false);
    }

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->d(iz), false);
      } else {
        Rp->initFromNode(argM, Ru->d(iz), false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
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
        node_handle oldi = C->d(i);
        long cev = Inf<long>();
        node_handle cnode = 0;
        accumulateOp->compute(nev, newstates, C->ei(i), oldi, cev, cnode);
        C->setEdge(i, cev);
        C->d_ref(i) = cnode;

        resF->unlinkNode(oldi);
        resF->unlinkNode(newstates);
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
  compute_table::search_key* Key = findResult(ev, evmdd, mxd, resEv, resEvmdd);
  if (0==Key) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int evmddLevel = argV->getNodeLevel(evmdd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), evmddLevel);
  const int rSize = resF->getLevelSize(rLevel);
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize evmdd reader
  unpacked_node *A = (evmddLevel < rLevel)
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmdd, true)
    : unpacked_node::newFromNode(argV, evmdd, true);

  if (evmddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (int i=0; i<rSize; i++) {
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
    for (int i=0; i<rSize; i++) {
      C->setEdge(i, 0L);
      C->d_ref(i) = 0;
    }

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::useUnpackedNode();
    unpacked_node *Rp = unpacked_node::useUnpackedNode();
    if (mxdLevel < 0) {
      Ru->initRedundant(argM, rLevel, mxd, false);
    } else {
      Ru->initFromNode(argM, mxd, false);
    }

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      if (0==A->d(i))   continue;
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(argM, rLevel, i, Ru->d(iz), false);
      } else {
        Rp->initFromNode(argM, Ru->d(iz), false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
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
        node_handle oldj = C->d(j);
        long cev = Inf<long>();
        node_handle cnode = 0;
        accumulateOp->compute(nev, newstates, C->ei(j), oldj, cev, cnode);
        C->setEdge(j, cev);
        C->d_ref(j) = cnode;

        resF->unlinkNode(oldj);
        resF->unlinkNode(newstates);
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
  compute_table::search_key* Key = findResult(ev, evmxd, mxd, resEv, resEvmdd);
  if (0==Key) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int evmxdLevel = argV->getNodeLevel(evmxd);
  const int mxdLevel = argM->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), ABS(evmxdLevel));
  const int rSize = resF->getLevelSize(rLevel);
  unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(rLevel, evmxdLevel)
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmxd, true)
    : unpacked_node::newFromNode(argV, evmxd, true);

  for (int i = 0; i < rSize; i++) {
    int pLevel = argV->getNodeLevel(A->d(i));
    unpacked_node* B = isLevelAbove(-rLevel, pLevel)
      ? unpacked_node::newIdentity(argV, -rLevel, i, 0L, A->d(i), true)
      : unpacked_node::newFromNode(argV, A->d(i), true);

    unpacked_node* D = unpacked_node::newFull(resF, -rLevel, rSize);
    if (rLevel > ABS(mxdLevel)) {
      //
      // Skipped levels in the MXD,
      // that's an important special case that we can handle quickly.
      for (int j = 0; j < rSize; j++) {
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
      for (int j = 0; j < rSize; j++) {
        D->setEdge(j, 0L);
        D->d_ref(j) = 0;
      }

      // Initialize mxd readers, note we might skip the unprimed level
      unpacked_node *Ru = unpacked_node::useUnpackedNode();
      unpacked_node *Rp = unpacked_node::useUnpackedNode();
      if (mxdLevel < 0) {
        Ru->initRedundant(argM, rLevel, mxd, false);
      } else {
        Ru->initFromNode(argM, mxd, false);
      }

      // loop over mxd "rows"
      for (int jz = 0; jz < Ru->getNNZs(); jz++) {
        int j = Ru->i(jz);
        if (0 == B->d(j)) {
          continue;
        }

        if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->d(jz)))) {
          Rp->initIdentity(argM, rLevel, j, Ru->d(jz), false);
        } else {
          Rp->initFromNode(argM, Ru->d(jz), false);
        }

        // loop over mxd "columns"
        for (int jpz = 0; jpz < Rp->getNNZs(); jpz++) {
          int jp = Rp->i(jpz);
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
          node_handle oldjp = D->d(jp);
          long dev = Inf<long>();
          node_handle dnode = 0;
          accumulateOp->compute(nev, newstates, D->ei(jp), oldjp, dev, dnode);
          D->setEdge(jp, dev);
          D->d_ref(jp) = dnode;

          resF->unlinkNode(oldjp);
          resF->unlinkNode(newstates);
        } // for j

      } // for i

      unpacked_node::recycle(Rp);
      unpacked_node::recycle(Ru);
    } // else

    long cev = Inf<long>();
    node_handle cnode = 0;
    resF->createReducedNode(i, D, cev, cnode);
    C->setEdge(i, cev);
    C->d_ref(i) = cnode;

    unpacked_node::recycle(B);
  }

  // cleanup mdd reader
  unpacked_node::recycle(A);

  resF->createReducedNode(-1, C, resEv, resEvmdd);
#ifdef TRACE_ALL_OPS
  printf("computed new tcXrel(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmdd, mxd, resEv, resEvmdd);
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


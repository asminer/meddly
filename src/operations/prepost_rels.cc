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
#include "prepost_rels.h"

#include "../ops_builtin.h"
#include "../ct_vector.h"
#include "../oper_unary.h"
#include "../oper_binary.h"

#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../forest_levels.h"
#include "../forest_edgerules.h"

// #define TRACE_ALL_OPS
// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

namespace MEDDLY {
    class image_op_evplus2;
    class tcXrel_evplus;

    binary_list TC_POST_IMAGE_cache;
};

// ******************************************************************
// *                                                                *
// *                     image_op_evplus2  class                    *
// *                                                                *
// ******************************************************************

/// Abstract base class for all MT-based pre/post image operations.
class MEDDLY::image_op_evplus2 : public binary_operation {
  public:
    image_op_evplus2(binary_list& opcode, forest* arg1,
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


MEDDLY::image_op_evplus2::image_op_evplus2(binary_list& oc, forest* a1,
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

void MEDDLY::image_op_evplus2::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag)
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

void MEDDLY::image_op_evplus2::compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  MEDDLY_DCASSERT(accumulateOp);
  compute_rec(ev, evmdd, mxd, resEv, resEvmdd);
}

// ******************************************************************
// *                                                                *
// *                      tcXrel_evplus class                       *
// *                                                                *
// ******************************************************************

/** Generic base for transitive closure multiplied by relation.
    Changing what happens at the terminals can give
    different meanings to this operation :^)
*/
class MEDDLY::tcXrel_evplus : public image_op_evplus2 {
  public:
    tcXrel_evplus(binary_list& opcode, forest* tc,
      forest* trans, forest* res, binary_operation* acc);

  protected:
    virtual void compute_rec(long ev, node_handle evmxd, node_handle mxd, long& resEv, node_handle& resEvmxd);
    virtual void processTerminals(long ev, node_handle evmxd, node_handle mxd, long& resEv, node_handle& resEvmxd);
};

MEDDLY::tcXrel_evplus::tcXrel_evplus(binary_list& oc,
  forest* tc, forest* trans, forest* res, binary_operation* acc)
: image_op_evplus2(oc, tc, trans, res, acc)
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
  unpacked_node* C = unpacked_node::newWritable(resF, rLevel, rSize, FULL_ONLY);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(rLevel, evmxdLevel)
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmxd, FULL_ONLY)
    : unpacked_node::newFromNode(argV, evmxd, FULL_ONLY);

  for (unsigned i = 0; i < rSize; i++) {
    int pLevel = argV->getNodeLevel(A->down(i));
    unpacked_node* B = isLevelAbove(-rLevel, pLevel)
      ? unpacked_node::newIdentity(argV, -rLevel, i, 0L, A->down(i), FULL_ONLY)
      : unpacked_node::newFromNode(argV, A->down(i), FULL_ONLY);

    unpacked_node* D = unpacked_node::newWritable(resF, -rLevel, rSize, FULL_ONLY);
    if (rLevel > ABS(mxdLevel)) {
      //
      // Skipped levels in the MXD,
      // that's an important special case that we can handle quickly.
      for (unsigned j = 0; j < rSize; j++) {
        long nev = Inf<long>();
        node_handle newstates = 0;
        compute_rec(A->edgeval(i).getLong() + B->edgeval(j).getLong(), B->down(j), mxd, nev, newstates);

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
      unpacked_node *Ru = unpacked_node::New(argM, SPARSE_ONLY);
      unpacked_node *Rp = unpacked_node::New(argM, SPARSE_ONLY);
      if (mxdLevel < 0) {
        Ru->initRedundant(rLevel, mxd);
      } else {
        Ru->initFromNode(mxd);
      }

      dd_edge newstatesE(resF), djp(resF);

      // loop over mxd "rows"
      for (unsigned jz = 0; jz < Ru->getSize(); jz++) {
        unsigned j = Ru->index(jz);
        if (0 == B->down(j)) {
          continue;
        }

        if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->down(jz)))) {
          Rp->initIdentity(rLevel, j, Ru->down(jz));
        } else {
          Rp->initFromNode(Ru->down(jz));
        }

        // loop over mxd "columns"
        for (unsigned jpz = 0; jpz < Rp->getSize(); jpz++) {
          unsigned jp = Rp->index(jpz);
          // ok, there is an i->j "edge".
          // determine new states to be added (recursively)
          // and add them
          long nev = Inf<long>();
          node_handle newstates = 0;
          compute_rec(A->edgeval(i).getLong() + B->edgeval(j).getLong(), B->down(j), Rp->down(jpz), nev, newstates);
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
          djp.set(D->edgeval(jp).getLong(), D->down(jp));
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

MEDDLY::binary_operation* MEDDLY::TC_POST_IMAGE(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  TC_POST_IMAGE_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    return TC_POST_IMAGE_cache.add(
        new tcXrel_evplus(TC_POST_IMAGE_cache, a, b, c, build(UNION, c, c, c))
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


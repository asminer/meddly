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
// ************************************************************************

/*
    Template class for MT-based pre-image and post-image operations.
    This includes vector-matrix and matrix-vector multiplications.

    Template parameters:

        EOP: one of the EdgeOp classes in forests_edgerules.h.

        FORWD: if true, do post-image (one step forward), otherwise
                        do pre-image  (one step backward).

        ATYPE: arithmetic type class, must provide the following methods.

            /// Get the operation name, for display purposes
            static const char* name(bool forwd);

            /// Get the accumulate operation for result nodes.
            static binary_operation* accumulateOp(const forest* resF);

            /// Apply the operation when b is a terminal node.
            static void apply(
                const forest* fa, const edge_value &av, node_handle a,
                const forest* fb, node_handle b,
                const forest* fc, edge_value &cv, node_handle &c
            );

*/

namespace MEDDLY {

    template <class EOP, bool FORWD, class ATYPE>
    class prepost_set_op : public binary_operation {
        public:
            prepost_set_op(forest* arg1, forest* arg2,
                    forest* res, bool swap=false);

            virtual ~prepost_set_op();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    const edge_value &bv, node_handle bp,
                    edge_value &cv, node_handle &cp);

        protected:
            void _compute(int L, const edge_value &av, node_handle A,
                    node_handle B, edge_value &cv, node_handle &C);

        private:
            inline const edge_value &edgeval(unpacked_node *U, unsigned i) const
            {
                if (EOP::hasEdgeValues()) {
                    return U->edgeval(i);
                } else {
                    return nothing;
                }
            }

        private:
            ct_entry_type* ct;
            binary_operation* accumulateOp;
            const bool swap_opnds;
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
            edge_value nothing;

    }; // class prepost_op
}; // namespace MEDDLY

// ************************************************************************

template <class EOP, bool FORWD, class ATYPE>
MEDDLY::prepost_set_op<EOP, FORWD, ATYPE>
    ::prepost_set_op(forest* arg1, forest* arg2, forest* res, bool swap)
    : binary_operation(
            swap ? arg2 : arg1,
            swap ? arg1 : arg2, res), swap_opnds(swap)
#ifdef TRACE
      , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkRelations(__FILE__, __LINE__, SET, RELATION, SET);
    checkLabelings(__FILE__, __LINE__,
        res->getEdgeLabeling(),
        edge_labeling::MULTI_TERMINAL,
        res->getEdgeLabeling()
    );

    if (range_type::BOOLEAN == arg2F->getRangeType()) {
        if (arg1F->getRangeType() != resF->getRangeType()) {
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        }
    } else {
        checkAllRanges(__FILE__, __LINE__, res->getRangeType());
    }

    //
    // Addition operation for the vector-matrix multiply.
    //  union for pre/post image, addition otherwise.
    //
    accumulateOp = ATYPE::accumulateOp(res);
    MEDDLY_DCASSERT(accumulateOp);

    //
    // Build compute table key and result types.
    //
    ct = new ct_entry_type(ATYPE::name(FORWD));
    if (swap_opnds) {
        ct->setFixed(arg2, arg1);
    } else {
        ct->setFixed(arg1, arg2);
    }
    if (EOP::hasEdgeValues()) {
        ct->setResult(EOP::edgeValueTypeLetter(), res);
    } else {
        ct->setResult(res);
    }
    ct->doneBuilding();
}

template <class EOP, bool FORWD, class ATYPE>
MEDDLY::prepost_set_op<EOP, FORWD, ATYPE>::~prepost_set_op()
{
    ct->markForDestroy();
}

template <class EOP, bool FORWD, class ATYPE>
void MEDDLY::prepost_set_op<EOP, FORWD, ATYPE>
    ::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << ATYPE::name(FORWD) << " #" << top_count << " begin\n";
#endif

    MEDDLY_DCASSERT(bv.isVoid());

    if (swap_opnds) {
        _compute(L, bv, bp, ap, cv, cp);
    } else {
        _compute(L, av, ap, bp, cv, cp);
    }

#ifdef TRACE
    out << ATYPE::name(FORWD) << " #" << top_count << " end\n";
#endif
}

template <class EOP, bool FORWD, class ATYPE>
void MEDDLY::prepost_set_op<EOP, FORWD, ATYPE>::_compute(int L,
        const edge_value &av, node_handle A, node_handle B,
        edge_value &cv, node_handle &C)
{
    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************
    if (0==A || 0==B) {
        cv = av;
        C = 0;
        return;
    }

    if (arg2F->isTerminalNode(B)) {
        //
        // We're either at the bottom,
        // or the matrix is an identity (or a scalar times identity).
        // Treat that case quickly.
        //
        ATYPE::apply(arg1F, av, A, arg2F, B, resF, cv, C);
        return;
    }

    // **************************************************************
    //
    // Determine level information
    //
    // **************************************************************
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = ABS(arg2F->getNodeLevel(B));
    const int Clevel = MAX(Alevel, Blevel);
    const int Cnextlevel = MDD_levels::downLevel(Clevel);

#ifdef TRACE
    out << ATYPE::name(FORWD) << " prepost_set_op::compute("
        << L << ", " << A << ", " << B << ")\n";
    out << "A: #" << A << " ";
    arg1F->showNode(out, A, SHOW_DETAILS);
    out << "\n";
    out << "B: #" << B << " ";
    arg2F->showNode(out, B, SHOW_DETAILS);
    out << "\n";
    // out << A << " level " << Alevel << "\n";
    // out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << "\n";
#endif

    // **************************************************************
    //
    // Check the compute table
    //
    // **************************************************************
    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    key[0].setN(A);
    key[1].setN(B);

    if (ct->findCT(key, res)) {
        //
        // compute table hit
        //
        if (EOP::hasEdgeValues()) {
            res[0].get(cv);
            EOP::accumulateOp(cv, av);
            C = resF->linkNode(res[1].getN());
        } else {
            EOP::clear(cv);
            C = resF->linkNode(res[0].getN());
        }
#ifdef TRACE
        out << "CT hit ";
        key.show(out);
        out << " -> ";
        res.show(out);
        out << "\n";
#endif
        C = resF->makeRedundantsTo(C, Clevel, L);
        return;
        //
        // done compute table hit
        //
    }

    // **************************************************************
    //
    // Compute table 'miss'; do computation
    //
    // **************************************************************

    //
    // Set up unpacked nodes
    //

    unpacked_node* Au = unpacked_node::New(arg1F, SPARSE_ONLY);
    if (Alevel != Clevel) {
        Au->initRedundant(Clevel, A);
    } else {
        Au->initFromNode(A);
    }

    rel_node* Brn;
    if (Blevel != Clevel) {
        Brn = nullptr;
    } else {
        Brn = arg2F->buildRelNode(B);
    }

    unpacked_node* Cu = nullptr;

#ifdef TRACE
    edge_value nothing;
    out << "A: ";
    Au->show(out, true);
    out << "\nB: ";
    if (Brn) {
        Brn->show(out);
    } else {
        out << "identity to " << B;
    }
    out.indent_more();
    out.put('\n');
#endif

    //
    // Recurse
    //
    if (!Brn) {
        //
        // Identity level
        //
        Cu = unpacked_node::newWritable(resF, Clevel, SPARSE_ONLY);
        unsigned zc = 0;
        for (unsigned z=0; z<Au->getSize(); z++) {
            edge_value  cdv;
            node_handle cdp;
            _compute(Cnextlevel, edgeval(Au, z), Au->down(z), B, cdv, cdp);
            if (!resF->isTransparentEdge(cdv, cdp)) {
                Cu->setSparse(zc, Au->index(z), cdv, cdp);
                ++zc;
            }
        }
        Cu->resize(zc);
    } else {
        //
        // Non-identity level
        //
        Cu = unpacked_node::newWritable(resF, Clevel, FULL_ONLY);
        // clear out result (important!)
        Cu->clear(0, Cu->getSize());

        if (FORWD) {
            unpacked_node* Bu = unpacked_node::New(arg2F, SPARSE_ONLY);

            for (unsigned zi=0; zi<Au->getSize(); zi++) {
                const unsigned i = Au->index(zi);
                if (Brn->outgoing(i, *Bu)) {
                    for (unsigned zj=0; zj<Bu->getSize(); zj++) {
                        const unsigned j = Bu->index(zj);
                        // C[j] = C[j] + A[i] * B[i,j]
                        node_handle cdp;
                        edge_value  cdv;
                        _compute(Cnextlevel, edgeval(Au, zi), Au->down(zi),
                                Bu->down(zj), cdv, cdp);
                        if (EOP::isZeroFunction(cdv, cdp)) {
                            continue;
                        }
                        if (EOP::isZeroFunction(edgeval(Cu, j), Cu->down(j)))
                        {
                            Cu->setFull(j, cdv, cdp);
                            continue;
                        }

                        edge_value  newdv;
                        node_handle newdp;
                        accumulateOp->compute(L, ~0,
                            edgeval(Cu, j), Cu->down(j),
                            cdv, cdp,
                            newdv, newdp
                        );
                        resF->unlinkNode(cdp);
                        resF->unlinkNode(Cu->down(j));
                        Cu->setFull(j, newdv, newdp);
                    } // for zj
                } // if brn[i]
            } // for zi
        } else {
            const unsigned kSize = unsigned(resF->getLevelSize(Clevel));
            unpacked_node* Bu = unpacked_node::New(arg2F, FULL_ONLY);
            for (unsigned i=0; i<kSize; i++) {
                if (Brn->outgoing(i, *Bu)) {
                    for (unsigned zj=0; zj<Au->getSize(); zj++) {
                        const unsigned j = Au->index(zj);
                        if (Bu->down(j)) {
                            // C[i] = C[i] + B[i,j] * A[j]
                            node_handle cdp;
                            edge_value  cdv;
                            _compute(Cnextlevel, edgeval(Au, zj), Au->down(zj),
                                    Bu->down(j), cdv, cdp);
                            if (EOP::isZeroFunction(cdv, cdp)) {
                                continue;
                            }
                            if (EOP::isZeroFunction(edgeval(Cu, i), Cu->down(i)))
                            {
                                Cu->setFull(i, cdv, cdp);
                                continue;
                            }
                            edge_value  newdv;
                            node_handle newdp;
                            accumulateOp->compute(L, ~0,
                                edgeval(Cu, i), Cu->down(i),
                                cdv, cdp,
                                newdv, newdp
                            );
                            resF->unlinkNode(cdp);
                            resF->unlinkNode(Cu->down(i));
                            Cu->setFull(i, newdv, newdp);
                        } // if bu[j]
                    } // for zj
                } // if brn[i]
            } // for i
        }
    }

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << ATYPE::name(FORWD) << " prepost_set_op::compute("
        << L << ", " << A << ", " << B << ") done\n";
    out << "  A: #" << A << ": ";
    arg1F->showNode(out, A, SHOW_DETAILS);
    out << "\n  B: #" << B << ": ";
    if (Brn) {
        Brn->show(out);
    } else {
        out << "identity to " << B;
    }
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce
    //
    resF->createReducedNode(Cu, cv, C);
#ifdef TRACE
    out << "reduced to ";
    resF->showEdge(out, cv, C);
    out << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT
    //
    if (EOP::hasEdgeValues()) {
        res[0].set(cv);
        res[1].setN(C);
    } else {
        res[0].setN(C);
    }
    ct->addCT(key, res);

    //
    // Cleanup
    //
    arg2F->doneRelNode(Brn);
    unpacked_node::Recycle(Au);

    //
    // Adjust result
    //
    C = resF->makeRedundantsTo(C, Clevel, L);
    EOP::accumulateOp(cv, av);
}

// ******************************************************************
// *                                                                *
// *                       mt_prepost  struct                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    struct mt_prepost {
        inline static const char* name(bool forwd)
        {
            return forwd ? "post-image" : "pre-image";
        }

        /// Get the accumulate operation for result nodes.
        inline static binary_operation* accumulateOp(forest* resF)
        {
            return UNION(resF, resF, resF);
        }

        /// Apply the operation when b is a terminal node.
        static void apply(forest* fa, const edge_value &av, node_handle a,
                          const forest* fb, node_handle b,
                          forest* fc, edge_value &cv, node_handle &c)
        {
            MEDDLY_DCASSERT(b<0);
            unary_operation* copy = COPY(fa, fc);
            MEDDLY_DCASSERT(copy);
            edge_value dummy;
            dummy.set();
            copy->compute(fa->getNodeLevel(a), ~0,
                av, a, cv, c);
        }

    };
};

// ******************************************************************
// *                                                                *
// *                      mt_vectXmatr  struct                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class TTYPE>
    struct mt_vectXmatr {
        inline static const char* name(bool forwd)
        {
            return forwd ? "VM-multiply" : "MV-multiply";
        }

        /// Get the accumulate operation for result nodes.
        inline static binary_operation* accumulateOp(forest* resF)
        {
            return PLUS(resF, resF, resF);
        }

        /// Apply the operation when b is a terminal node.
        static void apply(forest* fa, const edge_value &av, node_handle a,
                          const forest* fb, node_handle b,
                          forest* fc, edge_value &cv, node_handle &c)
        {
            MEDDLY_DCASSERT(b<0);
            TTYPE mxdval;
            fb->getValueFromHandle(b, mxdval);
            binary_operation* mult = MULTIPLY(fa, fa, fc);
            node_handle fa_mxdval = fa->handleForValue(mxdval);

            edge_value dummy;
            dummy.set();
            mult->compute(fa->getNodeLevel(a), ~0,
                av, a, dummy, fa_mxdval,
                cv, c);
        }

    };
};



// ************************************************************************
// ************************************************************************
// ************************************************************************
// *                                                                      *
// *                                                                      *
// *                                                                      *
// *                      OLD actual  operations OLD                      *
// *                                                                      *
// *                                                                      *
// *                                                                      *
// ************************************************************************


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
  unpacked_node* C = unpacked_node::newWritable(resF, rLevel, rSize, FULL_ONLY);

  // Initialize evmdd reader
  unpacked_node *A = (evmddLevel < rLevel)
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmdd, FULL_ONLY)
    : unpacked_node::newFromNode(argV, evmdd, FULL_ONLY);

  if (evmddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      long nev = Inf<long>();
      node_handle newstates = 0;
      compute_rec(A->edgeval(i).getLong(), A->down(i), mxd, nev, newstates);

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
    unpacked_node *Ru = unpacked_node::New(argM, SPARSE_ONLY);
    unpacked_node *Rp = unpacked_node::New(argM, SPARSE_ONLY);
    if (mxdLevel < 0) {
      Ru->initRedundant(rLevel, mxd);
    } else {
      Ru->initFromNode(mxd);
    }

    dd_edge newstatesE(resF), cdi(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      unsigned i = Ru->index(iz);
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(rLevel, i, Ru->down(iz));
      } else {
        Rp->initFromNode(Ru->down(iz));
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
        compute_rec(A->edgeval(j).getLong(), A->down(j), Rp->down(jz), nev, newstates);
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
        cdi.set(C->edgeval(i).getLong(), C->down(i));
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
  unpacked_node* C = unpacked_node::newWritable(resF, rLevel, rSize, FULL_ONLY);

  // Initialize evmdd reader
  unpacked_node *A = (evmddLevel < rLevel)
    ? unpacked_node::newRedundant(argV, rLevel, 0L, evmdd, FULL_ONLY)
    : unpacked_node::newFromNode(argV, evmdd, FULL_ONLY);

  if (evmddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      long nev = Inf<long>();
      node_handle newstates = 0;
      compute_rec(A->edgeval(i).getLong(), A->down(i), mxd, nev, newstates);

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
    unpacked_node *Ru = unpacked_node::New(argM, SPARSE_ONLY);
    unpacked_node *Rp = unpacked_node::New(argM, SPARSE_ONLY);
    if (mxdLevel < 0) {
      Ru->initRedundant(rLevel, mxd);
    } else {
      Ru->initFromNode(mxd);
    }

    dd_edge newstatesE(resF), cdj(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      unsigned i = Ru->index(iz);
      if (0==A->down(i))   continue;
      if (isLevelAbove(-rLevel, argM->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(rLevel, i, Ru->down(iz));
      } else {
        Rp->initFromNode(Ru->down(iz));
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        unsigned j = Rp->index(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        long nev = Inf<long>();
        node_handle newstates = 0;
        compute_rec(A->edgeval(i).getLong(), A->down(i), Rp->down(jz), nev, newstates);
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
        cdj.set(C->edgeval(j).getLong(), C->down(j));
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
// ******************************************************************
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
            new prepost_set_op<EdgeOp_none, false, mt_prepost>(a, b, c)
        );
    }
    /*
    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {
        return PRE_IMAGE_cache.add(
            new mtmatr_evplusvect<int>(PRE_IMAGE_cache, a, b, c, acc)
        );
    }
    */
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
            new prepost_set_op<EdgeOp_none, true, mt_prepost>(a, b, c)
        );
    }
    /*
    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {
        return POST_IMAGE_cache.add(
            new evplusvect_mtmatr<int>(POST_IMAGE_cache, a, b, c, acc)
        );
    }
    */
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

    switch (c->getRangeType()) {
        case range_type::INTEGER:
            return VM_MULTIPLY_cache.add(
                new prepost_set_op<EdgeOp_none, true, mt_vectXmatr<int> >(a, b, c)
            );

        case range_type::REAL:
            return VM_MULTIPLY_cache.add(
                new prepost_set_op<EdgeOp_none, true, mt_vectXmatr<float> >(a, b, c)
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

    switch (c->getRangeType()) {
        case range_type::INTEGER:
            return MV_MULTIPLY_cache.add(
                new prepost_set_op<EdgeOp_none, false, mt_vectXmatr<int> >(a, b, c, true)
            );

        case range_type::REAL:
            return MV_MULTIPLY_cache.add(
                new prepost_set_op<EdgeOp_none, false, mt_vectXmatr<float> >(a, b, c, true)
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


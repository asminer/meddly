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
#include "prepost_sets.h"

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
    class VM_MULTIPLY_factory;
    class MV_MULTIPLY_factory;
    // class image_op_evplus;
    // class relXset_evplus;
    // class setXrel_evplus;
    // class tcXrel_evplus;

    // binary_list TC_POST_IMAGE_cache;
    // binary_list VM_MULTIPLY_cache;
    // binary_list MV_MULTIPLY_cache;
};

// ************************************************************************
// ************************************************************************

/*
    Template class for pre-image and post-image operations,
    when the relation is MT-based.
    This includes vector-matrix and matrix-vector multiplications.

    Template parameters:

        EOP: one of the EdgeOp classes in forests_edgerules.h.

        FORWD: if true, do post-image (one step forward), otherwise
                        do pre-image  (one step backward).

        ATYPE: arithmetic type class, must provide the following methods.

            /// Get the operation name, for display purposes
            static const char* name(bool forwd);

            /// Return true if the given edge is unreachable
            static bool isUnreachable(const edge_value &cv, node_handle c);

            /// Set an edge to be unreachable
            static void setUnreachable(edge_value &cv, node_handle &c);

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
    class prepost_set_mtrel : public binary_operation {
        public:
            prepost_set_mtrel(forest* arg1, forest* arg2,
                    forest* res, bool swap=false);

            virtual ~prepost_set_mtrel();

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

            //
            // Correctly do C[i] = C[i] + <v, p>
            //
            inline void addToCi(int nextL, unpacked_node *C, unsigned i,
                    const edge_value &v, node_handle p)
            {
                // This case should be caught already
                MEDDLY_DCASSERT( !ATYPE::isUnreachable(v, p) );
                if (ATYPE::isUnreachable(edgeval(C, i), C->down(i)))
                {
                    C->setFull(i, v, p);
                    return;
                }
                edge_value  newdv;
                node_handle newdp;
                accumulateOp->compute(nextL, ~0,
                    edgeval(C, i), C->down(i), v, p, newdv, newdp
                );
                resF->unlinkNode(p);
                resF->unlinkNode(C->down(i));
                C->setFull(i, newdv, newdp);
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
            bool forced_by_levels;

    }; // class prepost_op
}; // namespace MEDDLY

// ************************************************************************

template <class EOP, bool FORWD, class ATYPE>
MEDDLY::prepost_set_mtrel<EOP, FORWD, ATYPE>
    ::prepost_set_mtrel(forest* arg1, forest* arg2, forest* res, bool swap)
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
    // Do we need to recurse by levels and store level info in the CT?
    // YES, if the set and relation are both fully-reduced.
    // (If the set is quasi reduced, we will recurse by levels anyway.)
    // (If the relation is identity-reduced, we can skip levels.)
    //
    forced_by_levels = arg1->isFullyReduced() && arg2->isFullyReduced();

    //
    // Build compute table key and result types.
    //
    ct = new ct_entry_type(ATYPE::name(FORWD));
    if (forced_by_levels) {
        if (swap_opnds) {
            ct->setFixed('I', arg2, arg1);
        } else {
            ct->setFixed('I', arg1, arg2);
        }
    } else {
        if (swap_opnds) {
            ct->setFixed(arg2, arg1);
        } else {
            ct->setFixed(arg1, arg2);
        }
    }
    if (EOP::hasEdgeValues()) {
        ct->setResult(EOP::edgeValueTypeLetter(), res);
    } else {
        ct->setResult(res);
    }
    ct->doneBuilding();
}

template <class EOP, bool FORWD, class ATYPE>
MEDDLY::prepost_set_mtrel<EOP, FORWD, ATYPE>::~prepost_set_mtrel()
{
    ct->markForDestroy();
}

template <class EOP, bool FORWD, class ATYPE>
void MEDDLY::prepost_set_mtrel<EOP, FORWD, ATYPE>
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
void MEDDLY::prepost_set_mtrel<EOP, FORWD, ATYPE>::_compute(int L,
        const edge_value &av, node_handle A, node_handle B,
        edge_value &cv, node_handle &C)
{
    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************
    if (0==B || ATYPE::isUnreachable(av, A)) {
        ATYPE::setUnreachable(cv, C);
        EOP::accumulateOp(cv, av);
        return;
    }

    if (arg2F->isTerminalNode(B) && (0==L || arg2F->isIdentityReduced())) {
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
    const int Clevel = forced_by_levels ? L : MAX(Alevel, Blevel);
    const int nextL = MDD_levels::downLevel(Clevel);

#ifdef TRACE
    out << ATYPE::name(FORWD) << " prepost_set_mtrel::compute(" << L << ", ";
    arg1F->showEdge(out, av, A);
    out << ", " << B << ")\n";
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
    if (forced_by_levels) {
        key[0].setI(L);
        key[1].setN(A);
        key[2].setN(B);
    } else {
        key[0].setN(A);
        key[1].setN(B);
    }

    if (ct->findCT(key, res)) {
        //
        // compute table hit
        //
        if (EOP::hasEdgeValues()) {
            res[0].get(cv);
            EOP::accumulateOp(cv, av);
            C = resF->linkNode(res[1].getN());
            EOP::normalize(cv, C);
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
        edge_value zero;
        EOP::clear(zero);
        Au->initRedundant(Clevel, zero, A);
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
        if (arg2F->isFullyReduced()) {
            //
            // Skipped Fully level(s)
            //
            // Forward:  C[j] = C[j] + A[i] * B, for all i,j
            // Backward: C[i] = C[i] + B * A[j], for all i,j
            //
            // That's really the same thing, so we'll do:
            // C = 0
            // for all i s.t. A[i] != 0
            //     tmp = A[i] * B
            //     for all j
            //         C[j] = C[j] + tmp
            //
            Cu = unpacked_node::newWritable(resF, Clevel, FULL_ONLY);
            // clear out result (important!)
            Cu->clear(0, Cu->getSize());
            for (unsigned zi=0; zi<Au->getSize(); zi++) {
                node_handle ab_p;
                edge_value  ab_v;
                _compute(nextL, edgeval(Au, zi), Au->down(zi),
                        B, ab_v, ab_p);
                if (ATYPE::isUnreachable(ab_v, ab_p)) {
                    continue;
                }
                for (unsigned j=0; j<Cu->getSize(); j++) {
                    addToCi(nextL, Cu, j, ab_v, resF->linkNode(ab_p));
                }
                resF->unlinkNode(ab_p);
            } // for zi
        } else {
            //
            // Skipped Identity level(s)
            //
            MEDDLY_DCASSERT(arg2F->isIdentityReduced());
            Cu = unpacked_node::newWritable(resF, Clevel, SPARSE_ONLY);
            unsigned zc = 0;
            for (unsigned z=0; z<Au->getSize(); z++) {
                edge_value  ab_v;
                node_handle ab_p;
                _compute(nextL, edgeval(Au, z), Au->down(z), B, ab_v, ab_p);
                if (!resF->isTransparentEdge(ab_v, ab_p)) {
                    Cu->setSparse(zc, Au->index(z), ab_v, ab_p);
                    ++zc;
                }
            }
            Cu->resize(zc);
        }
    } else {
        //
        // Non-identity level
        //
        Cu = unpacked_node::newWritable(resF, Clevel, FULL_ONLY);
        // clear out result (important!)
        Cu->clear(0, Cu->getSize());

        if (FORWD) {
            /*
             * Go forward one step.
             */
            unpacked_node* Bu = unpacked_node::New(arg2F, SPARSE_ONLY);

            for (unsigned zi=0; zi<Au->getSize(); zi++) {
                const unsigned i = Au->index(zi);
                if (Brn->outgoing(i, *Bu)) {
                    for (unsigned zj=0; zj<Bu->getSize(); zj++) {
                        const unsigned j = Bu->index(zj);
                        // C[j] = C[j] + A[i] * B[i,j]
                        node_handle cdp;
                        edge_value  cdv;
                        _compute(nextL, edgeval(Au, zi), Au->down(zi),
                                Bu->down(zj), cdv, cdp);
                        if (ATYPE::isUnreachable(cdv, cdp)) {
                            continue;
                        }
                        addToCi(nextL, Cu, j, cdv, cdp);
#ifdef TRACE
                        out << A << "x" << B << " completed "
                            << i << "->" << j << "; C is now ";

                        Cu->show(out, false);
                        out << "\n";
#endif
                    } // for zj
                } // if brn[i]
            } // for zi
        } else {
            /*
             * Go backward one step.
             */
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
                            _compute(nextL, edgeval(Au, zj), Au->down(zj),
                                    Bu->down(j), cdv, cdp);
                            if (ATYPE::isUnreachable(cdv, cdp)) {
                                continue;
                            }
                            addToCi(nextL, Cu, i, cdv, cdp);
                        } // if bu[j]
#ifdef TRACE
                        out << A << "x" << B << " completed "
                            << i << "<-" << j << "; C is now ";

                        Cu->show(out, false);
                        out << "\n";
#endif
                    } // for zj
                } // if brn[i]
            } // for i
        }
    }

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << ATYPE::name(FORWD) << " prepost_set_mtrel::compute(" << L << ", ";
    arg1F->showEdge(out, av, A);
    out << ", " << B << ") done\n";
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
    EOP::normalize(cv, C);
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
            return build(UNION, resF, resF, resF);
        }

        /// Does the edge correspond to an unreachable state?
        inline static bool isUnreachable(const edge_value &ev, node_handle p)
        {
            MEDDLY_DCASSERT(ev.isVoid());
            return 0==p;
        }
        /// Set edge to be unreachable states
        inline static void setUnreachable(edge_value &v, node_handle &p)
        {
            v.set();
            p = 0;
        }

        /// Apply the operation when b is a terminal node.
        static void apply(forest* fa, const edge_value &av, node_handle a,
                          const forest* fb, node_handle b,
                          forest* fc, edge_value &cv, node_handle &c)
        {
            MEDDLY_DCASSERT(b<0);
            unary_operation* copy = build(COPY, fa, fc);
            MEDDLY_DCASSERT(copy);
            copy->compute(fa->getNodeLevel(a), ~0, av, a, cv, c);
        }

    };
};

// ******************************************************************
// *                                                                *
// *                       ev_prepost  struct                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <typename INT>
    struct ev_prepost {
        inline static const char* name(bool forwd)
        {
            return forwd ? "post-image" : "pre-image";
        }

        /// Get the accumulate operation for result nodes.
        inline static binary_operation* accumulateOp(forest* resF)
        {
            return build(MINIMUM, resF, resF, resF);
        }

        /// Does the edge correspond to an unreachable state?
        inline static bool isUnreachable(const edge_value &ev, node_handle p)
        {
            return OMEGA_INFINITY == p;
        }
        /// Set edge to be unreachable states
        inline static void setUnreachable(edge_value &v, node_handle &p)
        {
            p = OMEGA_INFINITY;
            v = INT(0);
        }


        /// Apply the operation when b is a terminal node.
        static void apply(forest* fa, const edge_value &av, node_handle a,
                          const forest* fb, node_handle b,
                          forest* fc, edge_value &cv, node_handle &c)
        {
            MEDDLY_DCASSERT(b<0);
            unary_operation* copy = build(COPY, fa, fc);
            MEDDLY_DCASSERT(copy);
            copy->compute(fa->getNodeLevel(a), ~0, av, a, cv, c);

            if (OMEGA_INFINITY != c) {
                cv.add(INT(1));
            }
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
            return build(PLUS, resF, resF, resF);
        }

        /// Does the edge correspond to an unreachable state?
        inline static bool isUnreachable(const edge_value &ev, node_handle p)
        {
            MEDDLY_DCASSERT(ev.isVoid());
            return 0==p;
        }
        inline static void setUnreachable(edge_value &v, node_handle &p)
        {
            v.set();
            p = 0;
        }

        /// Apply the operation when b is a terminal node.
        static void apply(forest* fa, const edge_value &av, node_handle a,
                          const forest* fb, node_handle b,
                          forest* fc, edge_value &cv, node_handle &c)
        {
            MEDDLY_DCASSERT(b<0);
            TTYPE mxdval;
            fb->getValueFromHandle(b, mxdval);
            binary_operation* mult = build(MULTIPLY, fa, fa, fc);
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

#if 0

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

#endif

// ******************************************************************
// *                                                                *
// *                     relXset_evplus  class                      *
// *                                                                *
// ******************************************************************

#if 0

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
      compute_rec(A->edgeval(i), A->down(i), mxd, nev, newstates);

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
        compute_rec(A->edgeval(j), A->down(j), Rp->down(jz), nev, newstates);
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
        cdi.set(C->edgeval(i), C->down(i));
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
  resEv = long(Cev);
#ifdef TRACE_ALL_OPS
  printf("computed new relXset(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmdd, mxd, resEv, resEvmdd);
#endif
  saveResult(Key, ev, evmdd, mxd, resEv, resEvmdd);
}

#endif

// ******************************************************************
// *                                                                *
// *                     setXrel_evplus  class                      *
// *                                                                *
// ******************************************************************

#if 0

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
      compute_rec(A->edgeval(i), A->down(i), mxd, nev, newstates);

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
        compute_rec(A->edgeval(i), A->down(i), Rp->down(jz), nev, newstates);
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
        cdj.set(C->edgeval(j), C->down(j));
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
  resEv = long(Cev);
#ifdef TRACE_ALL_OPS
  printf("computed new setXrel(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmdd, mxd, resEv, resEvmdd);
#endif
  saveResult(Key, ev, evmdd, mxd, resEv, resEvmdd);
}

#endif

// ******************************************************************
// *                                                                *
// *                    mtmatr_evplusvect  class                    *
// *                                                                *
// ******************************************************************

#if 0

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

#endif

// ******************************************************************
// *                                                                *
// *                    evplusvect_mtmatr  class                    *
// *                                                                *
// ******************************************************************

#if 0

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

#endif

// ******************************************************************
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                 _IMAGE_factory  template class                 *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <bool FWD>
    class _IMAGE_factory : public binary_factory {
        public:
            virtual void setup();
            virtual binary_operation*
                build_new(forest* a, forest* b, forest* c);
        private:
            char docs[1024];
    };
};

// ******************************************************************

template <bool FWD>
void MEDDLY::_IMAGE_factory <FWD>::setup()
{
    snprintf(docs, 1024, "Follow edges %s in a transition relation. The first operand should be a set/vector. The second operand should be a relation/matrix. The result should be a set/vector, and should have the same range as the first operand. All forests should be over the same domain.\nIf the output function has type boolean, then the output is the set of states %s one of the input states.\nTBD...",
            FWD ? "forward" : "backward",
            FWD ? "reachable from" : "that can reach");

    _setup(__FILE__, FWD ? "POST_IMAGE" : "PRE_IMAGE", docs);
}

template <bool FWD>
MEDDLY::binary_operation*
MEDDLY::_IMAGE_factory <FWD>::build_new(forest* a, forest* b, forest* c)
{
    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        return new prepost_set_mtrel<EdgeOp_none, FWD, mt_prepost>(a, b, c);
    }

    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {

        switch (a->getEdgeType()) {
            case edge_type::INT:
                return new prepost_set_mtrel<EdgeOp_plus<int>, FWD,
                       ev_prepost<int> > (a, b, c);

            case edge_type::LONG:
                return new prepost_set_mtrel<EdgeOp_plus<long>, FWD,
                       ev_prepost<long> > (a, b, c);

            default:
                return nullptr;
        };
    }
    return nullptr;
}

// ******************************************************************
// *                                                                *
// *                   VM_MULTIPLY_factory  class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::VM_MULTIPLY_factory : public binary_factory {
    public:
        virtual void setup();
        virtual binary_operation* build_new(forest* a, forest* b, forest* c);
};

// ******************************************************************

void MEDDLY::VM_MULTIPLY_factory::setup()
{
    _setup(__FILE__, "VM_MULTIPLY", "Vector matrix multiply. The first operand is a vector (MDD), the second operand is a matrix (MxD), and the result is a vector (MDD). All forests must multi-terminal and over the same domain.");
}

MEDDLY::binary_operation*
MEDDLY::VM_MULTIPLY_factory::build_new(forest* a, forest* b, forest* c)
{
    switch (c->getRangeType()) {
        case range_type::INTEGER:
            return new prepost_set_mtrel<EdgeOp_none, true,
                            mt_vectXmatr<int> >(a, b, c);

        case range_type::REAL:
            return new prepost_set_mtrel<EdgeOp_none, true,
                            mt_vectXmatr<float> >(a, b, c);

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
}

// ******************************************************************
// *                                                                *
// *                   MV_MULTIPLY_factory  class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::MV_MULTIPLY_factory : public binary_factory {
    public:
        virtual void setup();
        virtual binary_operation* build_new(forest* a, forest* b, forest* c);
};

// ******************************************************************

void MEDDLY::MV_MULTIPLY_factory::setup()
{
    _setup(__FILE__, "MV_MULTIPLY", "Matrix vector multiply. The first operand is a matrix (MxD), the second operand is a vector (MDD), and the result is a vector (MDD). All forests must multi-terminal and over the same domain.");
}

MEDDLY::binary_operation*
MEDDLY::MV_MULTIPLY_factory::build_new(forest* a, forest* b, forest* c)
{
    switch (c->getRangeType()) {
        case range_type::INTEGER:
            return new prepost_set_mtrel<EdgeOp_none, false,
                            mt_vectXmatr<int> >(a, b, c, true);

        case range_type::REAL:
            return new prepost_set_mtrel<EdgeOp_none, false,
                            mt_vectXmatr<float> >(a, b, c, true);

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_factory& MEDDLY::PRE_IMAGE()
{
    static _IMAGE_factory<false> F;
    return F;
}

MEDDLY::binary_factory& MEDDLY::POST_IMAGE()
{
    static _IMAGE_factory<true> F;
    return F;
}

MEDDLY::binary_factory& MEDDLY::VM_MULTIPLY()
{
    static VM_MULTIPLY_factory F;
    return F;
}

MEDDLY::binary_factory& MEDDLY::MV_MULTIPLY()
{
    static MV_MULTIPLY_factory F;
    return F;
}


// ******************************************************************

#if 0

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
                new prepost_set_mtrel<EdgeOp_none, true, mt_vectXmatr<int> >(a, b, c)
            );

        case range_type::REAL:
            return VM_MULTIPLY_cache.add(
                new prepost_set_mtrel<EdgeOp_none, true, mt_vectXmatr<float> >(a, b, c)
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
                new prepost_set_mtrel<EdgeOp_none, false, mt_vectXmatr<int> >(a, b, c, true)
            );

        case range_type::REAL:
            return MV_MULTIPLY_cache.add(
                new prepost_set_mtrel<EdgeOp_none, false, mt_vectXmatr<float> >(a, b, c, true)
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

#endif

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
#include "reach_satur.h"

#include "../forest.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"
#include "../ct_vector.h"
#include "../forest_levels.h"
#include "../forest_edgerules.h"

#include "prepost_common.h"

// #define TRACE
// #define DEBUG_SPLIT

#ifdef TRACE
#include "../operators.h"
#endif


// ************************************************************************
// ************************************************************************

/*
    Template class for saturation (version 1) operations,
    when the relation is MT-based.

    Template parameters:

        EOP: one of the EdgeOp classes in forests_edgerules.h.

        FORWD: if true, do forward reachability, otherwise
                        do backward reachability.

        ATYPE: arithmetic type class, must provide the following methods.

            /// Get the operation name, for display purposes
            static const char* name(bool forwd);

            /// Return true if the given edge is unreachable
            static bool isUnreachable(const edge_value &cv, node_handle c);

            /// Set an edge to be unreachable
            static void setUnreachable(edge_value &cv, node_handle &c);

            /// Set all edges of an unpacked node to unreachable
            static void setAllUnreachable(unpacked_node *u);

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
    class saturation_set_mtrel : public binary_operation {
        public:
            saturation_set_mtrel(forest* arg1, forest* arg2,
                    forest* res);

            virtual ~saturation_set_mtrel();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    const edge_value &bv, node_handle bp,
                    edge_value &cv, node_handle &cp);

            virtual void showNameAndOptions(output &s) const;
            virtual option getOption(unsigned i) const;
        protected:
            virtual bool _setOption(unsigned i, option x);

            /* For an input set (av, A) from arg1F,
             * determine saturated node (cv, C) in resF,
             * with respect to level L relation(s).
             */
            void saturate(int L, const edge_value &av, node_handle A,
                    edge_value &cv, node_handle &C);

            void recFire(int L, const edge_value &av, node_handle A,
                    node_handle B, edge_value &cv, node_handle &c);

            static inline const char* opName() {
                return FORWD ? "fwd-sat" : "bck-sat";
            }

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
            /// Split relation: events whose top is exactly k
            std::vector <dd_edge> top_exactly;
            /// Split relation: events whose top is at k or below
            std::vector <dd_edge> top_at_or_below;

            ct_entry_type* firect;
            ct_entry_type* satct;
            unary_operation*  copyOp;
            binary_operation* accumulateOp;
            binary_operation* mxdIntersection;
            binary_operation* mxdDifference;
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
            // options
            int version;
            char index_order;

            edge_value nothing;
            bool forced_by_levels;
            bool sat_cache_level;

    }; // class
}; // namespace MEDDLY

// ************************************************************************

template <class EOP, bool FORWD, class ATYPE>
MEDDLY::saturation_set_mtrel<EOP, FORWD, ATYPE>
    ::saturation_set_mtrel(forest* arg1, forest* arg2, forest* res)
    : binary_operation(arg1, arg2, res)
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

    if (arg1F->getRangeType() != resF->getRangeType()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    //
    // Set the options for saturation
    //   option1: version number
    //   option2: index explore order
    //
    setOptionTypes("ic");
    setOption(0, 0L);
    setOption(1, ' ');


    // TBD: move split stuff; only use it for saturation v1
    //
    // Set up space for relation split by top levels
    //
    top_exactly.resize(arg2F->getNumVariables()+1);
    top_at_or_below.resize(arg2F->getNumVariables()+1);

    for (unsigned i=1; i<=arg2F->getNumVariables(); i++) {
        top_exactly[i].attach(arg2F);
        top_at_or_below[i].attach(arg2F);
    }

    //
    // Helper operations
    //
    copyOp = build(COPY, arg1F, resF);
    MEDDLY_DCASSERT(copyOp);

    accumulateOp = ATYPE::accumulateOp(res);
    MEDDLY_DCASSERT(accumulateOp);

    mxdIntersection = build(INTERSECTION, arg2F, arg2F, arg2F);
    MEDDLY_DCASSERT(mxdIntersection);

    mxdDifference = build(DIFFERENCE, arg2F, arg2F, arg2F);
    MEDDLY_DCASSERT(mxdDifference);

    //
    // Do we need to recurse by levels and store level info in the CT?
    // YES, if the set and relation are both fully-reduced.
    // (If the set is quasi reduced, we will recurse by levels anyway.)
    // (If the relation is identity-reduced, we can skip levels.)
    //
    forced_by_levels = arg1->isFullyReduced() && arg2->isFullyReduced();


    // TBD

    //
    // Build compute table key and result types.
    //
    firect = new ct_entry_type("satfire");
    if (forced_by_levels) {
        firect->setFixed('I', arg1, arg2);
    } else {
        firect->setFixed(arg1, arg2);
    }
    if (EOP::hasEdgeValues()) {
        firect->setResult(EOP::edgeValueTypeLetter(), res);
    } else {
        firect->setResult(res);
    }
    firect->doneBuilding();


    //
    // If the set is fully-reduced, we need to store the level in
    // the saturation cache.
    //
    satct = new ct_entry_type("saturate");
    sat_cache_level = arg1->isFullyReduced();
    if (sat_cache_level) {
        satct->setFixed('I', arg1, arg2);
    } else {
        satct->setFixed(arg1, arg2);
    }
    if (EOP::hasEdgeValues()) {
        satct->setResult(EOP::edgeValueTypeLetter(), res);
    } else {
        satct->setResult(res);
    }
    satct->doneBuilding();

}

template <class EOP, bool FORWD, class ATYPE>
MEDDLY::saturation_set_mtrel<EOP, FORWD, ATYPE>::~saturation_set_mtrel()
{
    firect->markForDestroy();
    satct->markForDestroy();
}


template <class EOP, bool FORWD, class ATYPE>
void
MEDDLY::saturation_set_mtrel<EOP, FORWD, ATYPE>
    ::showNameAndOptions(output &s) const
{
    if (FORWD) {
        s.put("Forward saturation version ");
    } else {
        s.put("Backward saturation version ");
    }
    s.put(version);
    s.put(", ");
    switch (index_order) {
        case 'q':
            s.put("queue index order");
            break;

        // TBD: other index orders

        default:
            s.put("unknown index order");
    }
}

template <class EOP, bool FORWD, class ATYPE>
MEDDLY::operation::option
MEDDLY::saturation_set_mtrel<EOP, FORWD, ATYPE>
    ::getOption(unsigned i) const
{
    if (0==i) {
        return option(long(version));
    }
    if (1==i) {
        return option(index_order);
    }
    throw error(error::INVALID_OPTION, __FILE__, __LINE__);
}

template <class EOP, bool FORWD, class ATYPE>
bool
MEDDLY::saturation_set_mtrel<EOP, FORWD, ATYPE>
    ::_setOption(unsigned i, option x)
{
    if (0==i) {
        // Version number
        if (x.opt_i < 0 || x.opt_i > 3) {
            return false;
        }
        if (0==x.opt_i) {
            // Use default
            version = 1;
            return true;
        }
        version = x.opt_i;
        return true;
    }
    if (1==i) {
        // index order
        if (' ' == x.opt_c) {
            // Use default
        }
        switch (x.opt_c) {

            case ' ':
                // Use default
                index_order = 'q';
                return true;

            case 'q':
                index_order = x.opt_c;
                return true;

            // tbd: other orders here


            default:
                // unknown or unsupported character
                return false;
        }
    }
    return false;
}

template <class EOP, bool FORWD, class ATYPE>
void MEDDLY::saturation_set_mtrel<EOP, FORWD, ATYPE>
    ::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
    //
    // Split the relation
    //
    MEDDLY_DCASSERT(bv.isVoid());

    dd_edge diag(arg2F);
    dd_edge mxd(arg2F);
    mxd.set(arg2F->linkNode(bp));
    for (;;)
    {
        node_handle resp;
        //
        // Get relation node for current mxd
        //
        const node_handle mxdn = mxd.getNode();
        const int k = ABS(arg2F->getNodeLevel(mxdn));
        if (0 == k) break;

        top_at_or_below[k] = mxd;
        rel_node* Brn = arg2F->buildRelNode(mxdn);

        // Determine common diagonal
        diag.set(Brn->getDiagonal(0));
        const unsigned maxi = arg2F->getLevelSize(k);
        for (unsigned i=1; i<maxi; i++) {
            mxdIntersection->compute(k, ~0,
                    nothing, diag.getNode(),
                    nothing, Brn->getDiagonal(i),
                    diag.setEdgeValue(), resp
            );
            diag.set(resp);
        }

        // Set relation with top=k to relation minus common diagonal
        // and continue the iteration with the common diagonal
        mxdDifference->compute(k, ~0,
            nothing, mxd.getNode(),
            nothing, diag.getNode(),
            top_exactly[k].setEdgeValue(), resp
        );
        top_exactly[k].set(resp);
        mxd = diag;

        // cleanup
        arg2F->doneRelNode(Brn);
    }
    top_at_or_below[0].set(0);
    top_exactly[0].set(0);

#ifdef DEBUG_SPLIT
    ostream_output splout(std::cout);
    std::cout << "After splitting monolithic event in " << opName() << "\n";
    for (unsigned k=0; k <= arg2F->getNumVariables(); k++) {
        std::cout << "Relation with top level<=" << k << ": ";
        top_at_or_below[k].show(splout);
        std::cout << "\n";
        top_at_or_below[k].showGraph(splout);
        std::cout << "======================================================================\n";
    }
    for (unsigned k=1; k <= arg2F->getNumVariables(); k++) {
        std::cout << "Relation with top level=" << k << ": ";
        top_exactly[k].show(splout);
        std::cout << "\n";
        top_exactly[k].showGraph(splout);
        std::cout << "======================================================================\n";
    }
#endif

#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << opName() << " #" << top_count << " begin\n";
#endif

    saturate(L, av, ap, cv, cp);

#ifdef TRACE
    out << opName() << " #" << top_count << " end\n";
#endif

    //
    // clear out split relation
    //
    for (int k=arg2F->getMaxLevelIndex(); k; --k) {
        top_exactly[k].set(0);
        top_at_or_below[k].set(0);
    }
}

template <class EOP, bool FORWD, class ATYPE>
void MEDDLY::saturation_set_mtrel<EOP, FORWD, ATYPE>::saturate(int L,
        const edge_value &av, node_handle A,
        edge_value &cv, node_handle &C)
{
    const node_handle B = top_at_or_below[L].getNode();

    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************
    if (ATYPE::isUnreachable(av, A)) {
        ATYPE::setUnreachable(cv, C);
        C = resF->makeRedundantsTo(C, 0, L);
        return;
    }
    if (0==B) {
        copyOp->compute(L, ~0, av, A, cv, C);
        return;
    }


#ifdef TRACE
    out << ATYPE::name(FORWD) << " saturate(" << L << ", ";
    arg1F->showEdge(out, av, A);
    out << ", " << B << ")\n";
#endif

    // **************************************************************
    //
    // Check the compute table
    //
    // **************************************************************
    ct_vector key(satct->getKeySize());
    ct_vector res(satct->getResultSize());
    if (sat_cache_level) {
        key[0].setI(L);
        key[1].setN(A);
        key[2].setN(B);
    } else {
        key[0].setN(A);
        key[1].setN(B);
    }

    if (satct->findCT(key, res)) {
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
    // Copy A to C, saturating children as we go
    //
    unpacked_node* Au = unpacked_node::New(arg1F, SPARSE_ONLY);
    const int Alevel = arg1F->getNodeLevel(A);
    if (Alevel < L) {
        edge_value zero;
        EOP::clear(zero);
        Au->initRedundant(L, zero, A);
    } else {
        Au->initFromNode(A);
    }

    unpacked_node* Cu = unpacked_node::newWritable(resF, L, FULL_ONLY);
    ATYPE::setAllUnreachable(Cu);

#ifdef TRACE
    out << "saturating children, node A: ";
    Au->show(out, true);
    out << "\n";
#endif

    for (unsigned z = 0; z<Au->getSize(); z++) {
        node_handle cdp;
        edge_value cdv;
        saturate(L-1, edgeval(Au, z), Au->down(z), cdv, cdp);
        const unsigned i = Au->index(z);
        Cu->setFull(i, cdv, cdp);
    }

    unpacked_node::Recycle(Au);

#ifdef TRACE
    out << "done saturating children, node C: ";
    Cu->show(out, true);
    out << "\n";
#endif


    // TBD
    MEDDLY_DCASSERT(false);

    //
    // OLD
    //

#ifdef OLD_SAT

#ifdef DEBUG_DFS
  printf("mdd: %d, k: %d\n", mdd, k);
#endif

  // terminal condition for recursion
  if (argF->isTerminalNode(mdd)) return mdd;

  // search compute table
  node_handle n = 0;
  ct_entry_key* Key = findSaturateResult(mdd, k, n);
  if (0==Key) return n;

  const unsigned sz = unsigned(argF->getLevelSize(k));    // size
  const int mdd_level = argF->getNodeLevel(mdd);          // mdd level

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d, mdd_level: %d\n",
      mdd, k, sz, mdd_level);
#endif

  unpacked_node* C = unpacked_node::newWritable(resF, k, sz, FULL_ONLY);
  // Initialize mdd reader
  unpacked_node *mddDptrs = unpacked_node::New(argF, FULL_ONLY);
  if (mdd_level < k) {
    mddDptrs->initRedundant(k, mdd);
  } else {
    mddDptrs->initFromNode(mdd);
  }

  // Do computation
  for (unsigned i=0; i<sz; i++) {
    C->setFull(i,
      mddDptrs->down(i) ? saturate(mddDptrs->down(i), k-1) : 0
    );
  }

  // Cleanup
  unpacked_node::Recycle(mddDptrs);

  parent->saturateHelper(*C);
  edge_value ev;
  resF->createReducedNode(C, ev, n);
  MEDDLY_DCASSERT(ev.isVoid());

  // save in compute table
  saveSaturateResult(Key, mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  return n;
#endif // OLD_SAT

}

template <class EOP, bool FORWD, class ATYPE>
void MEDDLY::saturation_set_mtrel<EOP, FORWD, ATYPE>::recFire(int L,
        const edge_value &av, node_handle A,
        node_handle B, edge_value &cv, node_handle &C)
{
    // TBD
    MEDDLY_DCASSERT(false);
}

// ******************************************************************
// *                                                                *
// *                  reachset_satur_factory class                  *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <bool FWD>
    class reachset_satur_factory : public binary_factory {
        public:
            virtual void setup();

            virtual binary_operation*
                build_new(forest* a, forest* b, forest* c);
    };
};

// ******************************************************************

template <bool FWD>
void MEDDLY::reachset_satur_factory <FWD>::setup()
{
    if (FWD) {
        _setup(__FILE__, "REACHABLE_SATUR(true)", "Build forward reachability set using saturation (citation? TBD). The first argument is the set of initial states, and the second argument is the transition relation.");
    } else {
        _setup(__FILE__, "REACHABLE_SATUR(false)", "Build backward reachability set using saturation (citation? TBD). The first argument is the set of initial states, and the second argument is the transition relation.");
    }
}

template <bool FWD>
MEDDLY::binary_operation*
MEDDLY::reachset_satur_factory <FWD>::build_new(forest* a, forest* b, forest* c)
{
    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {

        switch (c->getRangeType()) {
            case range_type::BOOLEAN:
                return new saturation_set_mtrel<EdgeOp_none, FWD,
                            mt_prepost>(a, b, c);

            case range_type::INTEGER:
                if (c->isFullyReduced())  {
                    return new saturation_set_mtrel<EdgeOp_none, FWD,
                            mt_distance>(a, b, c);
                }

            default:
                return nullptr;
        }
    }

    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {

        switch (a->getEdgeType()) {
            case edge_type::INT:
                return new saturation_set_mtrel<EdgeOp_plus<int>, FWD,
                            ev_prepost<int> > (a, b, c);

            case edge_type::LONG:
                return new saturation_set_mtrel<EdgeOp_plus<long>, FWD,
                            ev_prepost<long> > (a, b, c);

            default:
                return nullptr;
        };
    }
    return nullptr;
}




// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_factory& MEDDLY::REACHABLE_SATUR(bool fwd)
{
    static reachset_satur_factory<true>  forwd;
    static reachset_satur_factory<false> bckwd;

    if (fwd) return forwd;
    else     return bckwd;
}


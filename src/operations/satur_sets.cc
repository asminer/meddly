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
#include "satur_sets.h"
#include "satur_index.h"

#include "../forest.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"
#include "../ct_vector.h"
#include "../forest_levels.h"
#include "../forest_edgerules.h"

#include "prepost_common.h"
#include "../operators.h"

// #define RECFIRE_THEN_SAT

// #define TRACE
// #define DEBUG_SPLIT
// #define DEBUG_SPLIT_FULL
// #define TRACE_RECFIRE

// #define COUNT_CALLS

// ************************************************************************
// ************************************************************************

/*
    Template class for saturation (version 1) operations,
    when the relation is MT-based.

    Template parameters:

        EOP: one of the EdgeOp classes in forests_edgerules.h.

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

    template <class EOP, class ATYPE>
    class saturation_set_mtrel : public binary_operation {
        public:
            saturation_set_mtrel(const int version, bool fwd,
                    forest* arg1, forest* arg2, forest* res);

            virtual ~saturation_set_mtrel();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    const edge_value &bv, node_handle bp,
                    edge_value &cv, node_handle &cp);

        protected:

            /* For an input set (av, A) in resF,
             * determine saturated node (cv, C) in resF,
             * with respect to level L relation(s).
             * will saturate children.
             */
            void saturate_1(int L, const edge_value &av, node_handle A,
                    edge_value &cv, node_handle &C);

            /*
             * core of saturate_1. will NOT saturate children.
             */
            inline void saturate_1(unpacked_node *C)
            {
                MEDDLY_DCASSERT(C);
                if (top_exactly[C->getLevel()].getNode()) {
                    _saturate_1(C);
                }
            }

            void _saturate_1(unpacked_node *C);

            /* For an input set (av, A) in resF,
             * fire relation B, and in some cases saturate the result
             * with respect to level L relation(s),
             * and store the result in (cv, C) in resF.
             */
            void recFire(int L, const edge_value &av, node_handle A,
                    node_handle B, edge_value &cv, node_handle &c);

            inline const char* opName() {
                return FORWD ? "fwd-sat" : "bck-sat";
            }

        private:
            void initSplit();
            void fillSplit(int L, node_handle top);

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
            // Return true iff C[i] was updated
            //
            inline bool addToCi(int nextL, unpacked_node *C, unsigned i,
                    const edge_value &v, node_handle p)
            {
                if (ATYPE::isUnreachable(v, p)) return false;
                if (ATYPE::isUnreachable(edgeval(C, i), C->down(i)))
                {
                    C->setFull(i, v, p);
                    return true;
                }
                edge_value  newdv;
                node_handle newdp;
                accumulateOp->compute(nextL, ~0,
                    edgeval(C, i), C->down(i), v, p, newdv, newdp
                );
                resF->unlinkNode(p);
                bool changed;
                if (EOP::hasEdgeValues()) {
                    changed = (newdp != C->down(i)) || (newdv != edgeval(C,i));
                } else {
                    changed = (newdp != C->down(i));
                }
                resF->unlinkNode(C->down(i));
                C->setFull(i, newdv, newdp);
                return changed;
            }

        private:
            /// Split relation: events whose top is exactly k
            std::vector <dd_edge> top_exactly;
            /// Split relation: events whose top is at k or below
            std::vector <dd_edge> top_at_or_below;

            /// Helper for exploring indexes, by level.
            /// TBD: make this a template
            std::vector <satur_index_basic> explorers;

            ct_entry_type* fire_ct;
            ct_entry_type* sat_ct;
            unary_operation*  copyOp;
            binary_operation* accumulateOp;
            binary_operation* mxdIntersection;
            binary_operation* mxdDifference;
            const bool FORWD;
            const int VERSN;
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif

            edge_value nothing;
            // bool fire_cache_level;
            // bool sat_cache_level;

            // should we store levels in the cache
            bool store_levels;

#ifdef COUNT_CALLS
            size_t sat_calls;
            size_t recfire_calls;
#endif

    }; // class
}; // namespace MEDDLY

// ************************************************************************
// constructor
// ************************************************************************

template <class EOP, class ATYPE>
MEDDLY::saturation_set_mtrel<EOP, ATYPE>
    ::saturation_set_mtrel(const int version, bool fwd,
            forest* arg1, forest* arg2, forest* res)
    : binary_operation(arg1, arg2, res), FORWD(fwd), VERSN(version)
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

    initSplit();

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
    // Build keys and results for the compute tables.
    // If the set is fully-reduced, then we need to store the level in the CT.
    //
    store_levels = resF->isFullyReduced();

    fire_ct = new ct_entry_type("satfire");
    sat_ct  = new ct_entry_type("saturate");

    if (store_levels) {
        fire_ct->setFixed('I', resF, arg2F);
        sat_ct->setFixed('I', resF, arg2F);
    } else {
        fire_ct->setFixed(resF, arg2F);
        sat_ct->setFixed(resF, arg2F);
    }

    if (EOP::hasEdgeValues()) {
        fire_ct->setResult(EOP::edgeValueTypeLetter(), resF);
        sat_ct->setResult(EOP::edgeValueTypeLetter(), resF);
    } else {
        fire_ct->setResult(resF);
        sat_ct->setResult(resF);
    }

    fire_ct->doneBuilding();
    sat_ct->doneBuilding();

    //
    // Initialize explorers
    //
    explorers.resize(arg2F->getNumVariables()+1);
    for (unsigned i=1; i<explorers.size(); i++) {
        explorers[i].attach(arg2F, i, FORWD);
    }
}

// ************************************************************************
// destructor
// ************************************************************************

template <class EOP, class ATYPE>
MEDDLY::saturation_set_mtrel<EOP, ATYPE>::~saturation_set_mtrel()
{
    fire_ct->markForDestroy();
    sat_ct->markForDestroy();
}

// ************************************************************************
//
// compute()
//
// ************************************************************************

template <class EOP, class ATYPE>
void MEDDLY::saturation_set_mtrel<EOP, ATYPE>
    ::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(bv.isVoid());

    //
    // Split the relation
    //
    if (1==VERSN) fillSplit(L, bp);

#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << opName() << " #" << top_count << " begin\n";
#endif
#ifdef COUNT_CALLS
    sat_calls = 0;
    recfire_calls = 0;
#endif

    //
    // Copy (av, ap) from arg1F to resF, then operate entirely in resF
    //
    edge_value acv;
    node_handle acp;
    copyOp->compute(L, ~0, av, ap, acv, acp);

    //
    // Saturate
    //
    if (1==VERSN) {
        saturate_1(L, acv, acp, cv, cp);
    } else {
        MEDDLY_DCASSERT(false);
    }
    resF->unlinkNode(acp);

#ifdef TRACE
    out << opName() << " #" << top_count << " end\n";
#endif
#ifdef COUNT_CALLS
    std::cout << "#saturate calls: " << sat_calls << "\n";
    std::cout << "#recfire  calls: " << recfire_calls << "\n";
#endif

    // if (1==VERSN) fillSplit(L, 0);
}

// ************************************************************************
//
// saturate_1
//
// ************************************************************************

template <class EOP, class ATYPE>
void MEDDLY::saturation_set_mtrel<EOP, ATYPE>::saturate_1(int L,
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
        cv = av;
        C = resF->linkNode(A);
        return;
    }

#ifdef COUNT_CALLS
    ++sat_calls;
#endif

#ifdef TRACE
    out << ATYPE::name(FORWD) << " saturate_1(" << L << ", ";
    resF->showEdge(out, av, A);
    out << ", " << B << ")\n";
#endif

    // **************************************************************
    //
    // Check the compute table
    //
    // **************************************************************
    ct_vector key(sat_ct->getKeySize());
    ct_vector res(sat_ct->getResultSize());
    if (store_levels) {
        key[0].setI(L);
        key[1].setN(A);
        key[2].setN(B);
    } else {
        key[0].setN(A);
        key[1].setN(B);
    }

    if (sat_ct->findCT(key, res)) {
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
    unpacked_node* Au = unpacked_node::New(resF, SPARSE_ONLY);
    const int Alevel = resF->getNodeLevel(A);
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
        saturate_1(L-1, edgeval(Au, z), Au->down(z), cdv, cdp);
        const unsigned i = Au->index(z);
        Cu->setFull(i, cdv, cdp);
    }

    unpacked_node::Recycle(Au);
#ifdef DEVELOPMENT_CODE
    Au = nullptr;
#endif

#ifdef TRACE
    out << "done saturating children, node C: ";
    Cu->show(out, true);
    out << "\nsaturating this node";
    out.indent_more();
    out.put('\n');
#endif

    saturate_1(Cu);

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "done saturating node C: ";
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
    sat_ct->addCT(key, res);

    //
    // Adjust result
    //
    EOP::accumulateOp(cv, av);
    EOP::normalize(cv, C);
}

// ************************************************************************
//
// _saturate_1()
//
// ************************************************************************
template <class EOP, class ATYPE>
void MEDDLY::saturation_set_mtrel<EOP, ATYPE>::
    _saturate_1(unpacked_node *Cu)
{
    //
    // Initialize explorer
    //
    const int L = Cu->getLevel();
    unsigned i, j;
    node_handle d;
    for (i=0; i<Cu->getSize(); i++) {
        if (Cu->down(i)) {
            explorers[L].wasUpdated(i);
        }
    }
#ifdef TRACE
    explorers[L].show(out);
#endif
    //
    // Saturation loop :)
    //
    while (explorers[L].nextEdge(i, j, d)) {
#ifdef TRACE
        out << "firing " << i << "->" << j << " down " << d << "\n";
#endif
        if (ATYPE::areAllReachable(edgeval(Cu, j), Cu->down(j))) {
#ifdef TRACE
            out << "    target index all reachable; skipping\n";
#endif
            continue;
        }
        node_handle rfp;
        edge_value  rfv;
        recFire(L-1, edgeval(Cu, i), Cu->down(i), d, rfv, rfp);
        if (addToCi(L-1, Cu, j, rfv, rfp)) {
#ifdef TRACE
            out << "element " << j << " was updated\n";
#endif
            explorers[L].wasUpdated(j);
        }
#ifdef TRACE
        out << "after firing " << i << "->" << j << " down " << d
            << " node C is ";
        Cu->show(out, false);
        out.put('\n');
#endif
    }

}

// ************************************************************************
//
// recFire()
//
// ************************************************************************


template <class EOP, class ATYPE>
void MEDDLY::saturation_set_mtrel<EOP, ATYPE>::recFire(int L,
        const edge_value &av, node_handle A,
        node_handle B, edge_value &cv, node_handle &C)
{
    // **************************************************************
    //
    // Determine level information
    //
    // **************************************************************
    const int Alevel = resF->getNodeLevel(A);
    const int Blevel = ABS(arg2F->getNodeLevel(B));
    const int Clevel = L; // fire_cache_level ? L : MAX(Alevel, Blevel);
    const int nextL = MDD_levels::downLevel(Clevel);

    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************
    if (0==B || ATYPE::isUnreachable(av, A)) {
        ATYPE::setUnreachable(cv, C);
        C = resF->makeRedundantsTo(C, Clevel, L);
        return;
    }

    if (arg2F->isTerminalNode(B) && (0==L || arg2F->isIdentityReduced())) {
        //
        // We're either at the bottom,
        // or the matrix is an identity (or a scalar times identity).
        // Treat that case quickly.
        //
        ATYPE::apply(resF, av, A, arg2F, B, resF, cv, C);
        return;
    }

#ifdef COUNT_CALLS
    ++recfire_calls;
#endif

#ifdef TRACE
    out << ATYPE::name(FORWD) << " recFire(" << L << ", ";
    resF->showEdge(out, av, A);
    out << ", " << B << ")\n";
    out << "A: #" << A << " ";
    resF->showNode(out, A, SHOW_DETAILS);
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
    ct_vector key(fire_ct->getKeySize());
    ct_vector res(fire_ct->getResultSize());
    if (store_levels) {
        key[0].setI(L);
        key[1].setN(A);
        key[2].setN(B);
    } else {
        key[0].setN(A);
        key[1].setN(B);
    }

    if (fire_ct->findCT(key, res)) {
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

#ifdef RECFIRE_THEN_SAT
        node_handle oldC = C;
        saturate_1(L, cv, C, cv, C);
        resF->unlinkNode(oldC);
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

#ifdef TRACE_RECFIRE
    std::cout << "starting recfire(" << A << ", " << B << ")\n";
#endif

    //
    // Set up unpacked nodes
    //

    unpacked_node* Au = unpacked_node::New(resF, FULL_ONLY);
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
    // Initialize result
    //
    unpacked_node* Cu = unpacked_node::newWritable(resF, Clevel, FULL_ONLY);
    ATYPE::setAllUnreachable(Cu);

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
            // for all i s.t. A[i] != 0
            //     tmp = A[i] * B
            //     for all j
            //         C[j] = C[j] + tmp
            //
            for (unsigned i=0; i<Au->getSize(); i++) {
                if (ATYPE::isUnreachable(edgeval(Au, i), Au->down(i))) {
                    continue;
                }
                node_handle ab_p;
                edge_value  ab_v;
                recFire(nextL, edgeval(Au, i), Au->down(i), B, ab_v, ab_p);
                if (ATYPE::isUnreachable(ab_v, ab_p)) {
                    continue;
                }
                for (unsigned j=0; j<Cu->getSize(); j++) {
                    if (ATYPE::areAllReachable(edgeval(Cu, j), Cu->down(j))) {
                        continue;
                    }
                    addToCi(nextL, Cu, j, ab_v, resF->linkNode(ab_p));
                }
                resF->unlinkNode(ab_p);
            } // for zi
        } else {
            //
            //  Skipped Identity level(s)
            //
            //  For both forward and backward, compute
            //      C[i] = A[i] * B
            //  for all i.
            //
            MEDDLY_DCASSERT(arg2F->isIdentityReduced());
            for (unsigned i=0; i<Au->getSize(); i++) {
                if (ATYPE::isUnreachable(edgeval(Au, i), Au->down(i))) {
                    continue;
                }
                edge_value  ab_v;
                node_handle ab_p;
                recFire(nextL, edgeval(Au, i), Au->down(i), B, ab_v, ab_p);
                Cu->setFull(i, ab_v, ab_p);
            }
        }
    } else {
        //
        // Non-identity level
        //
        if (FORWD) {
            /*
             * Go forward one step.
             */
            unpacked_node* Bu = unpacked_node::New(arg2F, SPARSE_ONLY);

            for (unsigned i=0; i<Au->getSize(); i++) {
                if (ATYPE::isUnreachable(edgeval(Au, i), Au->down(i))) {
                    continue;
                }
                if (Brn->outgoing(i, *Bu)) {
                    for (unsigned zj=0; zj<Bu->getSize(); zj++) {
                        const unsigned j = Bu->index(zj);
                        if (ATYPE::areAllReachable(edgeval(Cu, j), Cu->down(j))) {
                            continue;
                        }
#ifdef TRACE
                        out << A << "x" << B << " computes  "
                            << i << ">->" << j << ": "
                            << Au->down(i) << "x" << Bu->down(zj) << "\n";
#endif
                        // C[j] = C[j] + A[i] * B[i,j]
                        node_handle cdp;
                        edge_value  cdv;
                        recFire(nextL, edgeval(Au, i), Au->down(i),
                                Bu->down(zj), cdv, cdp);
                        if (!ATYPE::isUnreachable(cdv, cdp)) {
                            addToCi(nextL, Cu, j, cdv, cdp);
                        }
#ifdef TRACE
                        out << A << "x" << B << " completed "
                            << i << ">->" << j << "; C is now ";

                        Cu->show(out, false);
                        out << "\n";
#endif
                    } // for zj
                } // if brn[i]
            } // for zi
            unpacked_node::Recycle(Bu);
        } else {
            /*
             * Go backward one step.
             */
            const unsigned kSize = unsigned(resF->getLevelSize(Clevel));
            unpacked_node* Bu = unpacked_node::New(arg2F, FULL_ONLY);
            for (unsigned i=0; i<kSize; i++) {
                if (Brn->outgoing(i, *Bu)) {
                    for (unsigned j=0; j<Au->getSize(); j++) {
                        if (ATYPE::isUnreachable(edgeval(Au, j), Au->down(j))) {
                            continue;
                        }
                        if (ATYPE::areAllReachable(edgeval(Cu, i), Cu->down(i))) {
                            break;
                        }
#ifdef TRACE
                        out << A << "x" << B << " computes  "
                            << i << "<-<" << j << ": "
                            << Au->down(j) << "x" << Bu->down(j) << "\n";
#endif
                        if (Bu->down(j)) {
                            // C[i] = C[i] + B[i,j] * A[j]
                            node_handle cdp;
                            edge_value  cdv;
                            recFire(nextL, edgeval(Au, j), Au->down(j),
                                    Bu->down(j), cdv, cdp);
                            if (!ATYPE::isUnreachable(cdv, cdp)) {
                                addToCi(nextL, Cu, i, cdv, cdp);
                            }
                        } // if bu[j]
#ifdef TRACE
                        out << A << "x" << B << " completed "
                            << i << "<-<" << j << "; C is now ";

                        Cu->show(out, false);
                        out << "\n";
#endif
                    } // for zj
                } // if brn[i]
            } // for i
            unpacked_node::Recycle(Bu);
        }
    }

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << ATYPE::name(FORWD) << " recFire(" << L << ", ";
    resF->showEdge(out, av, A);
    out << ", " << B << ") done\n";
    out << "  A: #" << A << ": ";
    resF->showNode(out, A, SHOW_DETAILS);
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
#ifdef TRACE_RECFIRE
    std::cout << "saturate recfire(" << A << ", " << B << ")\n";
#endif

#ifndef RECFIRE_THEN_SAT
    //
    // Saturate the unpacked node
    //
    saturate_1(Cu);
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
    // Cleanup
    //
    arg2F->doneRelNode(Brn);
    unpacked_node::Recycle(Au);

    //
    // Save result in CT
    //
    if (EOP::hasEdgeValues()) {
        res[0].set(cv);
        res[1].setN(C);
    } else {
        res[0].setN(C);
    }
    fire_ct->addCT(key, res);

    //
    // Adjust result
    //
    C = resF->makeRedundantsTo(C, Clevel, L);
    EOP::accumulateOp(cv, av);
    EOP::normalize(cv, C);

#ifdef RECFIRE_THEN_SAT
    //
    // Saturate result
    //
    node_handle oldC = C;
    saturate_1(L, cv, C, cv, C);
    resF->unlinkNode(oldC);
#endif

#ifdef TRACE_RECFIRE
    std::cout << "computed recfire(" << A << ", " << B << ") = " << C << "\n";
#endif
}

// ************************************************************************
// initSplit()
// ************************************************************************

template <class EOP, class ATYPE>
void MEDDLY::saturation_set_mtrel<EOP, ATYPE>::initSplit()
{
    if (VERSN != 1) {
        //
        // Clear relation split by top levels
        //
        top_exactly.clear();
        top_at_or_below.clear();
        return;
    }

    //
    // Set up space for relation split by top levels
    //
    top_exactly.resize(arg2F->getNumVariables()+1);
    top_at_or_below.resize(arg2F->getNumVariables()+1);

    for (unsigned i=1; i<=arg2F->getNumVariables(); i++) {
        top_exactly[i].attach(arg2F);
        top_at_or_below[i].attach(arg2F);
    }
}

// ************************************************************************
// fillSplit()
// ************************************************************************

template <class EOP, class ATYPE>
void MEDDLY::saturation_set_mtrel<EOP, ATYPE>::fillSplit(int L, node_handle bp)
{
    MEDDLY_DCASSERT(1==VERSN);
#ifdef DEBUG_SPLIT
    ostream_output splout(std::cout);
#endif

    //
    // fill split relation by levels
    //
    dd_edge diag(arg2F);
    dd_edge mxd(arg2F);
    mxd.set(arg2F->linkNode(bp));
    for (int k=L; k; --k)
    {
#ifdef DEBUG_SPLIT_FULL
        splout << "Splitting relation, at level " << k << "\n";
#endif
        node_handle resp;
        //
        // Get relation node for current mxd
        //
        const node_handle mxdn = mxd.getNode();
        top_at_or_below[k] = mxd;
#ifdef DEBUG_SPLIT_FULL
        splout << "    at or below: ";
        top_at_or_below[k].show(splout);
        splout << "\n";
#endif

        if (ABS(arg2F->getNodeLevel(mxdn)) < k) {
#ifdef DEBUG_SPLIT_FULL
            splout << "    no dependency on this level\n";
            splout << "    exactly: 0\n";
#endif
            top_exactly[k].set(0);
            continue;
        }

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

#ifdef DEBUG_SPLIT_FULL
        splout << "    diagonals: ";
        for (unsigned i=0; i<maxi; i++) {
            if (i) splout << ", ";
            splout << Brn->getDiagonal(i);
        }
        splout << "\n";
        splout << "    common  : ";
        diag.show(splout);
        splout << "\n";
#endif


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


#ifdef DEBUG_SPLIT_FULL
    splout << "After splitting monolithic event in " << opName() << "\n";
    for (unsigned k=0; k <= arg2F->getNumVariables(); k++) {
        splout << "Relation with top level<=" << k << ": ";
        top_at_or_below[k].show(splout);
        splout << "\n";
        if (ABS(arg2F->getNodeLevel(top_at_or_below[k].getNode())) == k) {
            top_at_or_below[k].showGraph(splout);
        }
        splout << "======================================================================\n";
    }
#endif
#ifdef DEBUG_SPLIT
    for (unsigned k=1; k <= arg2F->getNumVariables(); k++) {
        splout << "Relation with top level=" << k << ": ";
        top_exactly[k].show(splout);
        splout << "\n";
#ifdef DEBUG_SPLIT_FULL
        top_exactly[k].showGraph(splout);
        splout << "======================================================================\n";
#endif
    }
#endif


    /*
     * Set the explorers. We only need to do this once per split.
     */

    for (unsigned k=1; k<=arg2F->getNumVariables(); k++) {
        const node_handle n = top_exactly[k].getNode();
        explorers[k].restart(n);
    }
}

// ******************************************************************
// *                                                                *
// *                  reachset_satur_factory class                  *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <bool FWD, int VER>
    class reachset_satur_factory : public binary_factory {
        public:
            virtual void setup();

            virtual binary_operation*
                build_new(forest* a, forest* b, forest* c);
    };
};

// ******************************************************************

template <bool FWD, int VER>
void MEDDLY::reachset_satur_factory <FWD, VER>::setup()
{
    //
    // TBD: rewrite these
    //
    if (FWD) {
        _setup(__FILE__, "REACHABLE_SATUR(true)", "Build forward reachability set using saturation (citation? TBD). The first argument is the set of initial states, and the second argument is the transition relation. The operation has two options. Option 0 is an integer, to indicate the version number of saturation to use (0 for default, 1 for original). Option 1 is a character, to indicate how to select the next index to explore (' ' for default, 'q' whatever is next in the queue).");
    } else {
        _setup(__FILE__, "REACHABLE_SATUR(false)", "Build backward reachability set using saturation (citation? TBD). The first argument is the set of initial states, and the second argument is the transition relation. The operation has two options. Option 0 is an integer, to indicate the version number of saturation to use (0 for default, 1 for original). Option 1 is a character, to indicate how to select the next index to explore (' ' for default, 'q' whatever is next in the queue).");
    }
}

template <bool FWD, int VER>
MEDDLY::binary_operation*
MEDDLY::reachset_satur_factory <FWD, VER>::build_new(forest* a, forest* b, forest* c)
{
    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {

        switch (c->getRangeType()) {
            case range_type::BOOLEAN:
                return new saturation_set_mtrel<EdgeOp_none,
                            mt_prepost>(VER, FWD, a, b, c);

            case range_type::INTEGER:
                if (c->isFullyReduced())  {
                    return new saturation_set_mtrel<EdgeOp_none,
                            mt_distance>(VER, FWD, a, b, c);
                }

            default:
                return nullptr;
        }
    }

    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {

        switch (a->getEdgeType()) {
            case edge_type::INT:
                return new saturation_set_mtrel<EdgeOp_plus<int>,
                            ev_prepost<int> > (VER, FWD, a, b, c);

            case edge_type::LONG:
                return new saturation_set_mtrel<EdgeOp_plus<long>,
                            ev_prepost<long> > (VER, FWD, a, b, c);

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

MEDDLY::binary_factory& MEDDLY::REACHABLE_SATUR(bool fwd, int version)
{
    //
    // Version 1
    //
    static reachset_satur_factory<true, 1>  forwd1;
    static reachset_satur_factory<false, 1> bckwd1;

    //
    // Version 2
    //

    // TBD

    switch (version) {
        case 1:
            if (fwd) return forwd1;
            else     return bckwd1;

        default:
            return BOGUS_BINARY();
    }
}


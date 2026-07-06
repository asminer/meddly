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

#ifndef SATUR_SETS_V2_H
#define SATUR_SETS_V2_H

#include "../defines.h"
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

// ************************************************************************
// ************************************************************************

/*
    Template class for saturation (version 2) operations,
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
    class saturation2_set_mtrel : public binary_operation {
        public:
            saturation2_set_mtrel(bool fwd,
                    forest* arg1, forest* arg2, forest* res);

            virtual ~saturation2_set_mtrel();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    const edge_value &bv, node_handle bp,
                    edge_value &cv, node_handle &cp);

        protected:

            /* For an input set (av, A) in resF,
             * determine saturated node (cv, C) in resF,
             * with respect to relation R.
             * will saturate children.
             */
            void satbelow(int L, const edge_value &av, node_handle A, node_handle R,
                    edge_value &cv, node_handle &C);

            /*
             * core of saturate. will NOT saturate children.
             */
            inline void saturate(unpacked_node *C, node_handle R)
            {
                MEDDLY_DCASSERT(C);
                if (R) {
                    _saturate(C, R);
                }
            }

            void _saturate(unpacked_node *C, node_handle R);

            /* For an input set (av, A) in resF,
             * fire relation B, saturate the result using relation R,
             * and store the result in (cv, C) in resF.
             */
            void recFire(int L, const edge_value &av, node_handle A,
                    node_handle B, node_handle R, edge_value &cv, node_handle &c);

            inline const char* opName() {
                return FORWD ? "fwd-sat2" : "bck-sat2";
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
MEDDLY::saturation2_set_mtrel<EOP, ATYPE>
    ::saturation2_set_mtrel(bool fwd, forest* arg1, forest* arg2, forest* res)
    : binary_operation(arg1, arg2, res), FORWD(fwd)
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
        fire_ct->setFixed('I', resF, arg2F, arg2F);
        sat_ct->setFixed('I', resF, arg2F);
    } else {
        fire_ct->setFixed(resF, arg2F, arg2F);
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
MEDDLY::saturation2_set_mtrel<EOP, ATYPE>::~saturation2_set_mtrel()
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
void MEDDLY::saturation2_set_mtrel<EOP, ATYPE>
    ::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(bv.isVoid());

#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << opName() << " #" << top_count << " begin\n";

    out << "Relation:\n";
    dd_edge B(arg2F);
    B.set(bv, bp);
    B.showGraph(out);
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
    satbelow(L, acv, acp, bp, cv, cp);
    resF->unlinkNode(acp);

#ifdef TRACE
    out << opName() << " #" << top_count << " end\n";
#endif
#ifdef COUNT_CALLS
    std::cout << "#saturate calls: " << sat_calls << "\n";
    std::cout << "#recfire  calls: " << recfire_calls << "\n";
#endif
}

// ************************************************************************
//
// satbelow
//
// ************************************************************************

template <class EOP, class ATYPE>
void MEDDLY::saturation2_set_mtrel<EOP, ATYPE>::satbelow(int L,
        const edge_value &av, node_handle A, node_handle R,
        edge_value &cv, node_handle &C)
{
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
    if (0==R) {
        cv = av;
        C = resF->linkNode(A);
        return;
    }

#ifdef COUNT_CALLS
    ++sat_calls;
#endif

#ifdef TRACE
    out << ATYPE::name(FORWD) << " saturate_2(" << L << ", ";
    resF->showEdge(out, av, A);
    out << ", " << R << ")\n";
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
        key[2].setN(R);
    } else {
        key[0].setN(A);
        key[1].setN(R);
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
    // We'll need this eventually anyway, so let's use it
    // here for getting diagonals of R.
    //
    explorers[L].restart(R);
#ifdef TRACE
    out.indent_more();
    out << "Level " << L << " explorer:\n";
    explorers[L].show(out);
    out.indent_less();
    out.put('\n');
#endif

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
    bool Cempty = true;

#ifdef TRACE
    out << "saturating children, node A: ";
    Au->show(out, true);
    out << "\n";
#endif

    for (unsigned z = 0; z<Au->getSize(); z++) {
        node_handle cdp;
        edge_value cdv;
        const unsigned i = Au->index(z);
        explorers[L].confirm(i);
        satbelow(L-1, edgeval(Au, z), Au->down(z),
                explorers[L].getDiagonal(i), cdv, cdp);

        if (!ATYPE::isUnreachable(cdv, cdp)) {
            Cu->setFull(i, cdv, cdp);
            Cempty = false;
        }
    }

    unpacked_node::Recycle(Au);
#ifdef DEVELOPMENT_CODE
    Au = nullptr;
#endif

#ifdef TRACE
    out << "done saturating children, node C: ";
    Cu->show(out, true);
    out.put('\n');
#endif

    if (!Cempty && R) {

#ifdef TRACE
        out << "saturating this node";
        out.indent_more();
        out.put('\n');
#endif

        saturate(Cu, R);

#ifdef TRACE
        out.indent_less();
        out.put('\n');
        out << "done saturating node C: ";
        Cu->show(out, true);
        out << "\n";
#endif
    }

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
// _saturate()
//
// ************************************************************************
template <class EOP, class ATYPE>
void MEDDLY::saturation2_set_mtrel<EOP, ATYPE>::
    _saturate(unpacked_node *Cu, node_handle R)
{
    MEDDLY_DCASSERT(R);

    //
    // Initialize explorer
    //
    const int L = Cu->getLevel();
    unsigned i, j;
    node_handle d;
    explorers[L].restart(R);
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
        if (i == j) continue;

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
        explorers[L].confirm(j);
        recFire(L-1, edgeval(Cu, i), Cu->down(i), d, explorers[L].getDiagonal(j), rfv, rfp);
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
void MEDDLY::saturation2_set_mtrel<EOP, ATYPE>::recFire(int L,
        const edge_value &av, node_handle A, node_handle B,
        node_handle R, edge_value &cv, node_handle &C)
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
    out << ", " << B << ", " << R << ")\n";
    out << "A: #" << A << " ";
    resF->showNode(out, A, SHOW_DETAILS);
    out << "\n";
    out << "B: #" << B << " ";
    arg2F->showNode(out, B, SHOW_DETAILS);
    out << "\n";
    out << "R: #" << R << " ";
    arg2F->showNode(out, R, SHOW_DETAILS);
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
        key[3].setN(R);
    } else {
        key[0].setN(A);
        key[1].setN(B);
        key[1].setN(R);
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
    std::cout << "starting recfire(" << A << ", " << B << ", " << R << ")\n";
#endif

    //
    // We'll need this eventually anyway, so let's use it
    // here for getting diagonals of R.
    //
    explorers[L].restart(R);

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
    bool Cempty = true;

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
                explorers[L].confirm(i);
                recFire(nextL, edgeval(Au, i), Au->down(i), B,
                            explorers[L].getDiagonal(i), ab_v, ab_p);

                if (ATYPE::isUnreachable(ab_v, ab_p)) {
                    continue;
                }
                for (unsigned j=0; j<Cu->getSize(); j++) {
                    if (ATYPE::areAllReachable(edgeval(Cu, j), Cu->down(j))) {
                        continue;
                    }
                    addToCi(nextL, Cu, j, ab_v, resF->linkNode(ab_p));
                    Cempty = false;
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
                explorers[L].confirm(i);
                recFire(nextL, edgeval(Au, i), Au->down(i), B,
                            explorers[L].getDiagonal(i), ab_v, ab_p);

                if (!ATYPE::isUnreachable(ab_v, ab_p)) {
                    Cu->setFull(i, ab_v, ab_p);
                    Cempty = false;
                }
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
                        explorers[L].confirm(j);
                        recFire(nextL, edgeval(Au, i), Au->down(i), Bu->down(zj),
                                    explorers[L].getDiagonal(j), cdv, cdp);
                        if (!ATYPE::isUnreachable(cdv, cdp)) {
                            addToCi(nextL, Cu, j, cdv, cdp);
                            Cempty = false;
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
                            explorers[L].confirm(i);
                            recFire(nextL, edgeval(Au, j), Au->down(j),
                                    Bu->down(j), explorers[L].getDiagonal(i),
                                    cdv, cdp);
                            if (!ATYPE::isUnreachable(cdv, cdp)) {
                                addToCi(nextL, Cu, i, cdv, cdp);
                                Cempty = false;
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

    if (!Cempty && R) {

#ifdef TRACE_RECFIRE
        std::cout << "saturate2 recfire(" << A << ", " << B << ", " << R << ")\n";
#endif

        //
        // Saturate the unpacked node
        //
        saturate(Cu, R);
    }

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

#ifdef TRACE_RECFIRE
    std::cout << "computed recfire(" << A << ", " << B << ") = " << C << "\n";
#endif
}



#endif // include guard

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
#include "compare.h"

#include "../oper_item.h"
#include "../oper_binary.h"
#include "../ct_vector.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

//
// Operation instance caches
//
namespace MEDDLY {
    binary_list EQUAL_cache;
    binary_list NEQ_cache;

    binary_list GT_cache;
    binary_list GE_cache;

    binary_list LT_cache;
    binary_list LE_cache;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *            Template operations  for all comparisons            *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                        compare_mt class                        *
// *                                                                *
// ******************************************************************

/*
    Required methods for CTYPE classes (should be inlined)
    for multi-terminal comparisons:

        /// Get the operation name, for display purposes
        static const char* name();

        /// Is the comparison relation symmetric?
        /// Should be true for == and !=, false for the others.
        static bool isSymmetric();

        /// Is the comparison relation reflexive?
        /// I.e., does x ~ x always evaluate to true, or to false?
        static bool isReflexive();

        /// Compare terminals a and b.
        static bool compare(const forest* fa, node_handle a,
                            const forest* fb, node_handle b);

 */

namespace MEDDLY {
    template <class CTYPE>
    class compare_mt : public binary_operation {
        public:
            compare_mt(forest* arg1, forest* arg2, forest* res);
            virtual ~compare_mt();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    const edge_value &bv, node_handle bp,
                    edge_value &cv, node_handle &cp);

        protected:
            void _compute(int L, unsigned in, node_handle A, node_handle B,
                    node_handle &C);

        private:

            inline int topLevelOf(int L, int alevel, int blevel) const
            {
                if (forced_by_levels) return L;
                if (resF->isForRelations()) {
                    return MXD_levels::topLevel(alevel, blevel);
                } else {
                    return MDD_levels::topLevel(alevel, blevel);
                }
            }

        private:
            ct_entry_type* ct;
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
            bool forced_by_levels;
    };
};

// ******************************************************************

template <class CTYPE>
MEDDLY::compare_mt<CTYPE>::compare_mt(forest* arg1, forest* arg2,
        forest* res) : binary_operation(arg1, arg2, res)
#ifdef TRACE
      , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    if (res->isForRelations()) {
        checkAllRelations(__FILE__, __LINE__, RELATION);
    } else {
        checkAllRelations(__FILE__, __LINE__, SET);
    }
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    if (arg1->getRangeType() != arg2->getRangeType()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    //
    // Do we need to recurse by levels and store level info in the CT?
    // YES, if either forest is identity-reduced
    //      (even if the other is quasi-reduced, because we can
    //       still jump to terminal 0)
    //
    forced_by_levels = arg1->isIdentityReduced() || arg2->isIdentityReduced();

    // Build compute table key and result types.
    // If we recurse by levels, then we need the level as part of the key.
    ct = new ct_entry_type(CTYPE::name());
    if (forced_by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();
}

template <class CTYPE>
MEDDLY::compare_mt<CTYPE>::~compare_mt()
{
    ct->markForDestroy();
}

template <class CTYPE>
void MEDDLY::compare_mt<CTYPE>::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << CTYPE::name() << " #" << top_count << " begin\n";
#endif
    ASSERT(__FILE__, __LINE__, av.isVoid());
    ASSERT(__FILE__, __LINE__, bv.isVoid());

    _compute(L, in, ap, bp, cp);

#ifdef TRACE
    out << CTYPE::name() << " #" << top_count << " end\n";
#endif
    cv.set();
}


template <class CTYPE>
void MEDDLY::compare_mt<CTYPE>::_compute(int L, unsigned in,
        node_handle A, node_handle B, node_handle &C)
{
    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************

    //
    // Both are the zero function.
    //
    if (0==A && 0==B) {
        terminal tt(CTYPE::isReflexive(), resF->getTerminalType());
        C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
        return;
    }

    //
    // Arguments are equal.
    //
    if (A == B && arg1F == arg2F) {
        terminal tt(CTYPE::isReflexive(), resF->getTerminalType());
        C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
        return;
    }

    //
    // Both are constant functions.
    //
    if (arg1F->isTerminalNode(A) && arg2F->isTerminalNode(B)) {
        if (0==L || !forced_by_levels) {
            terminal tt( CTYPE::compare(arg1F, A, arg2F, B), resF->getTerminalType() );
            C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
            return;
        }
    }

    //
    // Reorder A and B if the comparison is symmetric
    // and A,B are from the same forest
    //
    if (CTYPE::isSymmetric() && arg1F == arg2F) {
        if (A > B) {
            SWAP(A, B);
        }
    }

    // **************************************************************
    //
    // Determine level information
    //
    // **************************************************************
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = topLevelOf(L, Alevel, Blevel);
    const int Cnextlevel = resF->isForRelations()
        ? MXD_levels::downLevel(Clevel)
        : MDD_levels::downLevel(Clevel)
    ;

#ifdef TRACE
    out << CTYPE::name() << " compare::compute(" << L << ", " << in << ", "
        << A << ", " << B << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
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
        C = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "CT hit ";
        key.show(out);
        out << " -> ";
        res.show(out);
        out << "\n";
#endif
        if (L == resF->getNodeLevel(C)) {
            // Make sure we don't point to a singleton
            // from the same index.
            C = resF->redirectSingleton(in, C);
        } else {
            C = resF->makeRedundantsTo(C, Clevel, L);
        }
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

    unpacked_node* Au = unpacked_node::New(arg1F, FULL_ONLY);
    if (Alevel != Clevel) {
        if (arg1F->isIdentityReduced() && Clevel<0) {
            ASSERT(__FILE__, __LINE__, Clevel == L);
            Au->initIdentity(Clevel, in, A);
            ASSERT(__FILE__, __LINE__, Au->wasIdentity());
        } else {
            Au->initRedundant(Clevel, A);
            ASSERT(__FILE__, __LINE__, !Au->wasIdentity());
        }
    } else {
        Au->initFromNode(A);
        ASSERT(__FILE__, __LINE__, !Au->wasIdentity());
    }

    unpacked_node* Bu = unpacked_node::New(arg2F, FULL_ONLY);
    if (Blevel != Clevel) {
        if (arg2F->isIdentityReduced() && Clevel<0) {
            ASSERT(__FILE__, __LINE__, Clevel == L);
            Bu->initIdentity(Clevel, in, B);
            ASSERT(__FILE__, __LINE__, Bu->wasIdentity());
        } else {
            Bu->initRedundant(Clevel, B);
            ASSERT(__FILE__, __LINE__, !Bu->wasIdentity());
        }
    } else {
        Bu->initFromNode(B);
        ASSERT(__FILE__, __LINE__, !Bu->wasIdentity());
    }

    unpacked_node* Cu = unpacked_node::newWritable(resF, Clevel, FULL_ONLY);

    ASSERT(__FILE__, __LINE__, Cu->getSize() == Au->getSize());
    ASSERT(__FILE__, __LINE__, Cu->getSize() == Bu->getSize());


#ifdef TRACE
    out << "A: ";
    Au->show(out, true);
    out << "\nB: ";
    Bu->show(out, true);
    out.indent_more();
    out.put('\n');
#endif

    //
    // Recurse
    //
    for (unsigned i=0; i<Cu->getSize(); i++) {
        node_handle cd;
        _compute(Cnextlevel, i, Au->down(i), Bu->down(i), cd);
        Cu->setFull(i, cd);
    }

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << CTYPE::name() << " compare::compute(" << L << ", " << in << ", "
        << A << ", " << B << ") done\n";
    out << "  A: ";
    Au->show(out, true);
    out << "\n  B: ";
    Bu->show(out, true);
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce
    //
    edge_value dummy;
    resF->createReducedNode(Cu, dummy, C);
    ASSERT(__FILE__, __LINE__, dummy.isVoid());
#ifdef TRACE
    out << "reduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT, if we can
    //
    if (Au->wasIdentity() || Bu->wasIdentity()) {
        ct->noaddCT(key);
    } else {
        res[0].setN(C);
        ct->addCT(key, res);
    }

    //
    // Cleanup
    //
    unpacked_node::Recycle(Bu);
    unpacked_node::Recycle(Au);

    //
    // Adjust result
    //
    if (L == resF->getNodeLevel(C)) {
        //
        // We don't need to add identities or redundants;
        // just check if we need to avoid a singleton.
        //
        C = resF->redirectSingleton(in, C);
    } else {
        //
        // Add nodes but only from Clevel;
        // if the actual level of C is below Clevel it means
        // nodes were eliminated in the result forest.
        //
        C = resF->makeRedundantsTo(C, Clevel, L);
    }
}

// ******************************************************************
// *                                                                *
// *                        compare_ev class                        *
// *                                                                *
// ******************************************************************

/*
    Required methods for CTYPE classes (should be inlined)
    for edge-valued comparisons:

        /// Get the operation name, for display purposes
        static const char* name();

        /// Is the comparison relation symmetric?
        /// Should be true for == and !=, false for the others.
        static bool isSymmetric();

        /// Is the comparison relation reflexive?
        /// I.e., does x ~ x always evaluate to true, or to false?
        static bool isReflexive();

        /// Is the comparison always true or false?
        /// E.g., x < infinity
        static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer);

        /// Compare constant functions <av, ap> and <bv, bp>
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp);


    Required methods for FACTOR classes (should be inlined)
    for edge-valued comparisons:

        ///
        /// Factor for <c,d> OP <e,f>.
        ///     On output, c will be -a + c, and e will be -a + e.
        ///
        static void factor(edge_value &c, node_handle d, edge_value &e,
            node_handle f, edge_value &a);

        ///
        /// Will factoring always produce the identity edge value
        /// for the first element?
        ///
        static bool alwaysFactorsToIdentity();



 */

namespace MEDDLY {
    template <class EOP, class FACTOR, class CTYPE>
    class compare_ev : public binary_operation {
        public:
            compare_ev(forest* arg1, forest* arg2, forest* res);
            virtual ~compare_ev();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    const edge_value &bv, node_handle bp,
                    edge_value &cv, node_handle &cp);

        protected:
            void _compute(int L, unsigned in,
                    edge_value av, node_handle A,
                    edge_value bv, node_handle B,
                    node_handle &C);

        private:

            inline int topLevelOf(int L, int alevel, int blevel) const
            {
                if (forced_by_levels) return L;
                if (resF->isForRelations()) {
                    return MXD_levels::topLevel(alevel, blevel);
                } else {
                    return MDD_levels::topLevel(alevel, blevel);
                }
            }

        private:
            ct_entry_type* ct;
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
            bool forced_by_levels;
    };
};

// ******************************************************************

template <class EOP, class FACTOR, class CTYPE>
MEDDLY::compare_ev<EOP, FACTOR, CTYPE>::compare_ev(forest* arg1,
        forest* arg2, forest* res) : binary_operation(arg1, arg2, res)
#ifdef TRACE
      , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    if (res->isForRelations()) {
        checkAllRelations(__FILE__, __LINE__, RELATION);
    } else {
        checkAllRelations(__FILE__, __LINE__, SET);
    }
    if (arg1->getRangeType() != arg2->getRangeType()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
    if (arg1->getEdgeLabeling() != arg2->getEdgeLabeling()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
    if (arg1->getEdgeType() != arg2->getEdgeType()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    //
    // Do we need to recurse by levels and store level info in the CT?
    // YES, if either forest is identity-reduced
    //      (even if the other is quasi-reduced, because we can
    //       still jump to terminal 0)
    //
    forced_by_levels = arg1->isIdentityReduced() || arg2->isIdentityReduced();

    // Build compute table key and result types.
    // If we recurse by levels, then we need the level as part of the key.
    ct = new ct_entry_type(CTYPE::name());
    if (forced_by_levels) {
        ct->appendFixed('I');
    }
    // If we can't always factor to identity, then we need to store
    // the edge value for the first operand.
    if (!FACTOR::alwaysFactorsToIdentity()) {
        ct->appendFixed(EOP::edgeValueTypeLetter());
    }
    ct->appendFixed(arg1);
    ct->appendFixed(EOP::edgeValueTypeLetter());
    ct->appendFixed(arg2);
    ct->setResult(res);
    ct->doneBuilding();
}

template <class EOP, class FACTOR, class CTYPE>
MEDDLY::compare_ev<EOP, FACTOR, CTYPE>::~compare_ev()
{
    ct->markForDestroy();
}

template <class EOP, class FACTOR, class CTYPE>
void MEDDLY::compare_ev<EOP, FACTOR, CTYPE>::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << CTYPE::name() << " #" << top_count << " begin\n";
#endif
    ASSERT(__FILE__, __LINE__, !av.isVoid());
    ASSERT(__FILE__, __LINE__, !bv.isVoid());

    _compute(L, in, av, ap, bv, bp, cp);

#ifdef TRACE
    out << CTYPE::name() << " #" << top_count << " end\n";
#endif
    cv.set();
}

template <class EOP, class FACTOR, class CTYPE>
void MEDDLY::compare_ev<EOP, FACTOR, CTYPE>::_compute(int L, unsigned in,
        edge_value av, node_handle A,
        edge_value bv, node_handle B, node_handle &C)
{
    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************

    //
    // Both are the zero function.
    //
    if (EOP::isZeroFunction(av, A) && EOP::isZeroFunction(bv, B))
    {
        terminal tt(CTYPE::isReflexive(), resF->getTerminalType());
        C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
        return;
    }

    //
    // Arguments are equal.
    //
    if ((arg1F == arg2F) && (A == B) && (av == bv)) {
        terminal tt(CTYPE::isReflexive(), resF->getTerminalType());
        C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
        return;
    }

    //
    // Both are constant functions.
    //
    if (arg1F->isTerminalNode(A) && arg2F->isTerminalNode(B)) {
        if (0==L || !forced_by_levels) {
            terminal tt( CTYPE::compare(av, A, bv, B), resF->getTerminalType() );
            C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
            return;
        }
    }

    //
    // Other special cases, like everything is less than infinity
    //
    bool answer;
    if (CTYPE::isSpecialCase(av, A, bv, B, answer)) {
        terminal tt(answer, resF->getTerminalType());
        C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
        return;
    }

    //
    // Reorder A and B if the comparison is symmetric
    // and A,B are from the same forest
    //
    if (CTYPE::isSymmetric() && arg1F == arg2F) {
        if (A > B) {
            SWAP(A, B);
            SWAP(av, bv);
        }
    }

    //
    // "normalize" the edge values
    //
    edge_value fac;
    FACTOR::factor(av, A, bv, B, fac);

    // **************************************************************
    //
    // Determine level information
    //
    // **************************************************************
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = topLevelOf(L, Alevel, Blevel);
    const int Cnextlevel = resF->isForRelations()
        ? MXD_levels::downLevel(Clevel)
        : MDD_levels::downLevel(Clevel)
    ;

#ifdef TRACE
    out << CTYPE::name() << " compare::compute(" << L << ", " << in << ", ";
    arg1F->showEdge(out, av, A);
    out << ", ";
    arg2F->showEdge(out, bv, B);
    out << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << "\n";
#endif

    // **************************************************************
    //
    // Check the compute table
    //
    // **************************************************************
    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    unsigned keyind=0;
    if (forced_by_levels) {
        key[keyind++].setI(L);
    }
    if (!FACTOR::alwaysFactorsToIdentity()) {
        key[keyind++].set(av);
    }
    key[keyind++].setN(A);
    key[keyind++].set(bv);
    key[keyind]  .setN(B);

    if (ct->findCT(key, res)) {
        //
        // compute table hit
        //
        C = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "CT hit ";
        key.show(out);
        out << " -> ";
        res.show(out);
        out << "\n";
#endif
        if (L == resF->getNodeLevel(C)) {
            // Make sure we don't point to a singleton
            // from the same index.
            C = resF->redirectSingleton(in, C);
        } else {
            C = resF->makeRedundantsTo(C, Clevel, L);
        }
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

    edge_value zero;
    EOP::clear(zero);
    unpacked_node* Au = unpacked_node::New(arg1F, FULL_ONLY);
    if (Alevel != Clevel) {
        if (arg1F->isIdentityReduced() && Clevel<0) {
            ASSERT(__FILE__, __LINE__, Clevel == L);
            Au->initIdentity(Clevel, in, zero, A);
            ASSERT(__FILE__, __LINE__, Au->wasIdentity());
        } else {
            Au->initRedundant(Clevel, zero, A);
            ASSERT(__FILE__, __LINE__, !Au->wasIdentity());
        }
    } else {
        Au->initFromNode(A);
        ASSERT(__FILE__, __LINE__, !Au->wasIdentity());
    }

    unpacked_node* Bu = unpacked_node::New(arg2F, FULL_ONLY);
    if (Blevel != Clevel) {
        if (arg2F->isIdentityReduced() && Clevel<0) {
            ASSERT(__FILE__, __LINE__, Clevel == L);
            Bu->initIdentity(Clevel, in, zero, B);
            ASSERT(__FILE__, __LINE__, Bu->wasIdentity());
        } else {
            Bu->initRedundant(Clevel, zero, B);
            ASSERT(__FILE__, __LINE__, !Bu->wasIdentity());
        }
    } else {
        Bu->initFromNode(B);
        ASSERT(__FILE__, __LINE__, !Bu->wasIdentity());
    }

    unpacked_node* Cu = unpacked_node::newWritable(resF, Clevel, FULL_ONLY);

    ASSERT(__FILE__, __LINE__, Cu->getSize() == Au->getSize());
    ASSERT(__FILE__, __LINE__, Cu->getSize() == Bu->getSize());

#ifdef TRACE
    out << "A: ";
    Au->show(out, true);
    out << "\nB: ";
    Bu->show(out, true);
    out.indent_more();
    out.put('\n');
#endif

    //
    // Recurse
    //
    for (unsigned i=0; i<Cu->getSize(); i++) {
        node_handle cd;
        _compute(Cnextlevel, i,
                EOP::applyOp(av, Au->edgeval(i)), Au->down(i),
                EOP::applyOp(bv, Bu->edgeval(i)), Bu->down(i),
                cd);
        Cu->setFull(i, cd);
    }

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << CTYPE::name() << " compare::compute(" << L << ", " << in << ", ";
    arg1F->showEdge(out, av, A);
    out << ", ";
    arg2F->showEdge(out, bv, B);
    out << ") done\n";
    out << "  A: ";
    Au->show(out, true);
    out << "\n  B: ";
    Bu->show(out, true);
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce
    //
    edge_value dummy;
    resF->createReducedNode(Cu, dummy, C);
    ASSERT(__FILE__, __LINE__, dummy.isVoid());
#ifdef TRACE
    out << "reduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT, if we can
    //
    if (Au->wasIdentity() || Bu->wasIdentity()) {
        ct->noaddCT(key);
    } else {
        res[0].setN(C);
        ct->addCT(key, res);
    }

    //
    // Cleanup
    //
    unpacked_node::Recycle(Bu);
    unpacked_node::Recycle(Au);

    //
    // Adjust result
    //
    if (L == resF->getNodeLevel(C)) {
        //
        // We don't need to add identities or redundants;
        // just check if we need to avoid a singleton.
        //
        C = resF->redirectSingleton(in, C);
    } else {
        //
        // Add nodes but only from Clevel;
        // if the actual level of C is below Clevel it means
        // nodes were eliminated in the result forest.
        //
        C = resF->makeRedundantsTo(C, Clevel, L);
    }
}

// ******************************************************************
// *                                                                *
// *                      evplus_factor  class                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EDGETYPE>
        struct evplus_factor {
            static void factor(edge_value &c, node_handle d, edge_value &e,
                    node_handle f, edge_value &a)
            {
                if (OMEGA_INFINITY == d) {
                    //
                    // left operand is infinity;
                    // factor using the right operand.
                    //
                    if (OMEGA_INFINITY == f) {
                        a = edge_value(EDGETYPE(0));
                    } else {
                        a = e;
                        e = edge_value(EDGETYPE(0));
                    }
                } else {
                    //
                    // Factor on left operand
                    //
                    a = c;
                    c = edge_value(EDGETYPE(0));
                    EDGETYPE av;
                    a.get(av);
                    e.subtract(av);
                }
            }
            static bool alwaysFactorsToIdentity()
            {
                return true;
            }
        };
};

// ******************************************************************
// *                                                                *
// *                      evstar_factor  class                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EDGETYPE>
        struct evstar_factor {
            inline static void factor(
                    edge_value &c, node_handle d,
                    edge_value &e, node_handle f,
                    edge_value &a)
            {
                if (OMEGA_ZERO == d) {
                    //
                    // left operand is zero;
                    // factor using the right operand.
                    //
                    if (OMEGA_ZERO == f) {
                        a = edge_value(EDGETYPE(1));
                    } else {
                        EDGETYPE av;
                        e.get(av);
                        if (av < 0) {
                            av = -av;
                            e = EDGETYPE(-1);
                        } else {
                            e = EDGETYPE(1);
                        }
                        a = av;
                    }
                } else {
                    //
                    // Factor on left operand
                    //
                    EDGETYPE av;
                    c.get(av);
                    // normalize left operand
                    if (av < 0) {
                        av = -av;
                        c = EDGETYPE(-1);
                    } else {
                        c = EDGETYPE(1);
                    }
                    a = av;
                    // normalize right operand
                    e.divide(av);
                }
            }
            static bool alwaysFactorsToIdentity()
            {
                return false;
            }
        };
};

// ******************************************************************
// *                                                                *
// *                          eq_  classes                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    struct eq_base {
        inline static const char* name() {
            return "==";
        }
        inline static bool isSymmetric() {
            return true;    // x == y is the same as y == x
        }
        inline static bool isReflexive() {
            return true;    // x == x is true
        }
    };
    template <class RANGE>
    struct eq_mt : public eq_base {
        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av == bv );
        }
    };
    template <class EDGETYPE>
    struct eq_evplus : public eq_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            if ((OMEGA_INFINITY == ap) && (OMEGA_INFINITY == bp))
            {
                return true;
            }
            if ((OMEGA_NORMAL == ap) && (OMEGA_NORMAL == bp))
            {
                EDGETYPE avv, bvv;
                av.get(avv);
                bv.get(bvv);
                return avv == bvv;
            }
            return false;
        }
    };
    template <class EDGETYPE>
    struct eq_evstar : public eq_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            if ((OMEGA_ZERO == ap) && (OMEGA_ZERO == bp))
            {
                return true;
            }
            if ((OMEGA_NORMAL == ap) && (OMEGA_NORMAL == bp))
            {
                EDGETYPE avv, bvv;
                av.get(avv);
                bv.get(bvv);
                return avv == bvv;
            }
            return false;
        }
    };
};

// ******************************************************************
// *                                                                *
// *                          ne_  classes                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    struct ne_base {
        inline static const char* name() {
            return "!=";
        }
        inline static bool isSymmetric() {
            return true;    // x != y is the same as y != x
        }
        inline static bool isReflexive() {
            return false;   // x != x is false
        }
    };
    template <class RANGE>
    struct ne_mt : public ne_base {
        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av != bv );
        }
    };
    template <class EDGETYPE>
    struct ne_evplus : public ne_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            if ((OMEGA_INFINITY == ap) && (OMEGA_INFINITY == bp))
            {
                return false;
            }
            if ((OMEGA_NORMAL == ap) && (OMEGA_NORMAL == bp))
            {
                EDGETYPE avv, bvv;
                av.get(avv);
                bv.get(bvv);
                return avv != bvv;
            }
            return true;
        }
    };
    template <class EDGETYPE>
    struct ne_evstar : public ne_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            if ((OMEGA_ZERO == ap) && (OMEGA_ZERO == bp))
            {
                return false;
            }
            if ((OMEGA_NORMAL == ap) && (OMEGA_NORMAL == bp))
            {
                EDGETYPE avv, bvv;
                av.get(avv);
                bv.get(bvv);
                return avv != bvv;
            }
            return true;
        }
    };
};

// ******************************************************************
// *                                                                *
// *                          gt_  classes                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    struct gt_base {
        inline static const char* name() {
            return ">";
        }
        inline static bool isSymmetric() {
            return false;   // x > y is not the same as y > x
        }
        inline static bool isReflexive() {
            return false;   // x > x is false
        }
    };
    template <class RANGE>
    struct gt_mt : public gt_base {
        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av > bv );
        }
    };
    template <class EDGETYPE>
    struct gt_evplus : public gt_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            // We can assume <av, ap> and <bv, bp>
            // are not both infinity, as that case
            // would have been handled already
            if (OMEGA_INFINITY == bp)
            {
                answer = false;
                return true;
            }
            if (OMEGA_INFINITY == ap)
            {
                answer = true;
                return true;
            }
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            if (OMEGA_INFINITY == bp)
            {
                return false;
            }
            if (OMEGA_INFINITY == ap)
            {
                return true;
            }
            EDGETYPE avv, bvv;
            av.get(avv);
            bv.get(bvv);
            return avv > bvv;
        }
    };
    template <class EDGETYPE>
    struct gt_evstar : public gt_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            EDGETYPE avv, bvv;
            if (OMEGA_ZERO == ap) {
                avv = 0;
            } else {
                av.get(avv);
            }
            if (OMEGA_ZERO == bp) {
                bvv = 0;
            } else {
                bv.get(bvv);
            }
            return avv > bvv;
        }
    };
};

// ******************************************************************
// *                                                                *
// *                          ge_  classes                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    struct ge_base {
        inline static const char* name() {
            return ">=";
        }
        inline static bool isSymmetric() {
            return false;   // x >= y is not the same as y >= x
        }
        inline static bool isReflexive() {
            return true;    // x >= x is true
        }
    };
    template <class RANGE>
    struct ge_mt : public ge_base {
        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av >= bv );
        }
    };
    template <class EDGETYPE>
    struct ge_evplus : public ge_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            // We can assume <av, ap> and <bv, bp>
            // are not both infinity, as that case
            // would have been handled already
            if (OMEGA_INFINITY == ap)
            {
                answer = true;
                return true;
            }
            if (OMEGA_INFINITY == bp)
            {
                answer = false;
                return true;
            }
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            if (OMEGA_INFINITY == ap)
            {
                return true;
            }
            if (OMEGA_INFINITY == bp)
            {
                return false;
            }
            EDGETYPE avv, bvv;
            av.get(avv);
            bv.get(bvv);
            return avv >= bvv;
        }
    };
    template <class EDGETYPE>
    struct ge_evstar : public ge_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            EDGETYPE avv, bvv;
            if (OMEGA_ZERO == ap) {
                avv = 0;
            } else {
                av.get(avv);
            }
            if (OMEGA_ZERO == bp) {
                bvv = 0;
            } else {
                bv.get(bvv);
            }
            return avv >= bvv;
        }
    };
};

// ******************************************************************
// *                                                                *
// *                          lt_  classes                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    struct lt_base {
        inline static const char* name() {
            return "<";
        }
        inline static bool isSymmetric() {
            return false;   // x < y is not the same as y < x
        }
        inline static bool isReflexive() {
            return false;   // x < x is false
        }
    };
    template <class RANGE>
    struct lt_mt : public lt_base {
        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av < bv );
        }
    };
    template <class EDGETYPE>
    struct lt_evplus : public lt_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            // We can assume <av, ap> and <bv, bp>
            // are not both infinity, as that case
            // would have been handled already
            if (OMEGA_INFINITY == ap)
            {
                answer = false;
                return true;
            }
            if (OMEGA_INFINITY == bp)
            {
                answer = true;
                return true;
            }
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            if (OMEGA_INFINITY == ap)
            {
                return false;
            }
            if (OMEGA_INFINITY == bp)
            {
                return true;
            }
            EDGETYPE avv, bvv;
            av.get(avv);
            bv.get(bvv);
            return avv < bvv;
        }
    };
    template <class EDGETYPE>
    struct lt_evstar : public lt_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            EDGETYPE avv, bvv;
            if (OMEGA_ZERO == ap) {
                avv = 0;
            } else {
                av.get(avv);
            }
            if (OMEGA_ZERO == bp) {
                bvv = 0;
            } else {
                bv.get(bvv);
            }
            return avv < bvv;
        }
    };
};

// ******************************************************************
// *                                                                *
// *                          le_  classes                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    struct le_base {
        inline static const char* name() {
            return "<=";
        }
        inline static bool isSymmetric() {
            return false;   // x <= y is not the same as y <= x
        }
        inline static bool isReflexive() {
            return true;    // x <= x is true
        }
    };
    template <class RANGE>
    struct le_mt : public le_base {
        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av <= bv );
        }
    };
    template <class EDGETYPE>
    struct le_evplus : public le_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            // We can assume <av, ap> and <bv, bp>
            // are not both infinity, as that case
            // would have been handled already
            if (OMEGA_INFINITY == ap)
            {
                answer = false;
                return true;
            }
            if (OMEGA_INFINITY == bp)
            {
                answer = true;
                return true;
            }
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            if (OMEGA_INFINITY == bp)
            {
                return true;
            }
            if (OMEGA_INFINITY == ap)
            {
                return false;
            }
            EDGETYPE avv, bvv;
            av.get(avv);
            bv.get(bvv);
            return avv <= bvv;
        }
    };
    template <class EDGETYPE>
    struct le_evstar : public le_base {
        inline static bool isSpecialCase(
                const edge_value& av, node_handle ap,
                const edge_value& bv, node_handle bp,
                bool &answer)
        {
            return false;
        }
        static bool compare(const edge_value& av, node_handle ap,
                            const edge_value& bv, node_handle bp)
        {
            EDGETYPE avv, bvv;
            if (OMEGA_ZERO == ap) {
                avv = 0;
            } else {
                av.get(avv);
            }
            if (OMEGA_ZERO == bp) {
                bvv = 0;
            } else {
                bv.get(bvv);
            }
            return avv <= bvv;
        }
    };
};


// ******************************************************************
// *                                                                *
// *                                                                *
// *                           Front ends                           *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                             EQUAL                              *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::EQUAL(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  EQUAL_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (a->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            bop = new compare_mt<eq_mt<float> > (a,b,c);
        } else {
            bop = new compare_mt<eq_mt<long> > (a,b,c);
        }
        return EQUAL_cache.add(bop);
    }

    if (a->isEVTimes()) {
        switch (a->getEdgeType()) {
            case edge_type::FLOAT:
                bop = new compare_ev<EdgeOp_times<float>,
                    evstar_factor<float>, eq_evstar<float> > (a,b,c);
                break;

            case edge_type::DOUBLE:
                bop = new compare_ev<EdgeOp_times<double>,
                    evstar_factor<double>, eq_evstar<double> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    } else {
        switch (a->getEdgeType()) {
            case edge_type::INT:
                bop = new compare_ev<EdgeOp_plus<int>, evplus_factor<int>,
                    eq_evplus<int> > (a,b,c);
                break;

            case edge_type::LONG:
                bop = new compare_ev<EdgeOp_plus<long>, evplus_factor<long>,
                    eq_evplus<long> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    }
    return EQUAL_cache.add( bop );
}

void MEDDLY::EQUAL_init()
{
    EQUAL_cache.reset("Equal");
}

void MEDDLY::EQUAL_done()
{
    ASSERT(__FILE__, __LINE__, EQUAL_cache.isEmpty());
}

// ******************************************************************
// *                            NOT_EQUAL                           *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::NOT_EQUAL(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  NEQ_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (a->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            bop = new compare_mt<ne_mt<float> > (a,b,c);
        } else {
            bop = new compare_mt<ne_mt<long> > (a,b,c);
        }
        return NEQ_cache.add(bop);
    }

    if (a->isEVTimes()) {
        switch (a->getEdgeType()) {
            case edge_type::FLOAT:
                bop = new compare_ev<EdgeOp_times<float>,
                    evstar_factor<float>, ne_evstar<float> > (a,b,c);
                break;

            case edge_type::DOUBLE:
                bop = new compare_ev<EdgeOp_times<double>,
                    evstar_factor<double>, ne_evstar<double> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    } else {
        switch (a->getEdgeType()) {
            case edge_type::INT:
                bop = new compare_ev<EdgeOp_plus<int>, evplus_factor<int>,
                    ne_evplus<int> > (a,b,c);
                break;

            case edge_type::LONG:
                bop = new compare_ev<EdgeOp_plus<long>, evplus_factor<long>,
                    ne_evplus<long> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    }
    return NEQ_cache.add( bop );
}

void MEDDLY::NOT_EQUAL_init()
{
    NEQ_cache.reset("Unequal");
}

void MEDDLY::NOT_EQUAL_done()
{
    ASSERT(__FILE__, __LINE__, NEQ_cache.isEmpty());
}

// ******************************************************************
// *                         GREATER  THAN                          *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::GREATER_THAN(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  GT_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (a->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            bop = new compare_mt<gt_mt<float> > (a,b,c);
        } else {
            bop = new compare_mt<gt_mt<long> > (a,b,c);
        }
        return GT_cache.add(bop);
    }

    if (a->isEVTimes()) {
        switch (a->getEdgeType()) {
            case edge_type::FLOAT:
                bop = new compare_ev<EdgeOp_times<float>,
                    evstar_factor<float>, gt_evstar<float> > (a,b,c);
                break;

            case edge_type::DOUBLE:
                bop = new compare_ev<EdgeOp_times<double>,
                    evstar_factor<double>, gt_evstar<double> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    } else {
        switch (a->getEdgeType()) {
            case edge_type::INT:
                bop = new compare_ev<EdgeOp_plus<int>, evplus_factor<int>,
                    gt_evplus<int> > (a,b,c);
                break;

            case edge_type::LONG:
                bop = new compare_ev<EdgeOp_plus<long>, evplus_factor<long>,
                    gt_evplus<long> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    }
    return GT_cache.add( bop );
}

void MEDDLY::GREATER_THAN_init()
{
    GT_cache.reset("MoreThan");
}

void MEDDLY::GREATER_THAN_done()
{
    ASSERT(__FILE__, __LINE__, GT_cache.isEmpty());
}

// ******************************************************************
// *                       GREATER  OR EQUAL                        *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::GREATER_THAN_EQUAL(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  GE_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (a->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            bop = new compare_mt<ge_mt<float> > (a,b,c);
        } else {
            bop = new compare_mt<ge_mt<long> > (a,b,c);
        }
        return GE_cache.add(bop);
    }

    if (a->isEVTimes()) {
        switch (a->getEdgeType()) {
            case edge_type::FLOAT:
                bop = new compare_ev<EdgeOp_times<float>,
                    evstar_factor<float>, ge_evstar<float> > (a,b,c);
                break;

            case edge_type::DOUBLE:
                bop = new compare_ev<EdgeOp_times<double>,
                    evstar_factor<double>, ge_evstar<double> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    } else {
        switch (a->getEdgeType()) {
            case edge_type::INT:
                bop = new compare_ev<EdgeOp_plus<int>, evplus_factor<int>,
                    ge_evplus<int> > (a,b,c);
                break;

            case edge_type::LONG:
                bop = new compare_ev<EdgeOp_plus<long>, evplus_factor<long>,
                    ge_evplus<long> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    }
    return GE_cache.add( bop );
}

void MEDDLY::GREATER_THAN_EQUAL_init()
{
    GE_cache.reset("MoreEqual");
}

void MEDDLY::GREATER_THAN_EQUAL_done()
{
    ASSERT(__FILE__, __LINE__, GE_cache.isEmpty());
}

// ******************************************************************
// *                           LESS THAN                            *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::LESS_THAN(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  LT_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (a->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            bop = new compare_mt<lt_mt<float> > (a,b,c);
        } else {
            bop = new compare_mt<lt_mt<long> > (a,b,c);
        }
        return LT_cache.add(bop);
    }

    if (a->isEVTimes()) {
        switch (a->getEdgeType()) {
            case edge_type::FLOAT:
                bop = new compare_ev<EdgeOp_times<float>,
                    evstar_factor<float>, lt_evstar<float> > (a,b,c);
                break;

            case edge_type::DOUBLE:
                bop = new compare_ev<EdgeOp_times<double>,
                    evstar_factor<double>, lt_evstar<double> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    } else {
        switch (a->getEdgeType()) {
            case edge_type::INT:
                bop = new compare_ev<EdgeOp_plus<int>, evplus_factor<int>,
                    lt_evplus<int> > (a,b,c);
                break;

            case edge_type::LONG:
                bop = new compare_ev<EdgeOp_plus<long>, evplus_factor<long>,
                    lt_evplus<long> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    }
    return LT_cache.add( bop );
}

void MEDDLY::LESS_THAN_init()
{
    LT_cache.reset("LessThan");
}

void MEDDLY::LESS_THAN_done()
{
    ASSERT(__FILE__, __LINE__, LT_cache.isEmpty());
}

// ******************************************************************
// *                         LESS OR EQUAL                          *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::LESS_THAN_EQUAL(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  LE_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (a->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            bop = new compare_mt<le_mt<float> > (a,b,c);
        } else {
            bop = new compare_mt<le_mt<long> > (a,b,c);
        }
        return LE_cache.add(bop);
    }

    if (a->isEVTimes()) {
        switch (a->getEdgeType()) {
            case edge_type::FLOAT:
                bop = new compare_ev<EdgeOp_times<float>,
                    evstar_factor<float>, le_evstar<float> > (a,b,c);
                break;

            case edge_type::DOUBLE:
                bop = new compare_ev<EdgeOp_times<double>,
                    evstar_factor<double>, le_evstar<double> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    } else {
        switch (a->getEdgeType()) {
            case edge_type::INT:
                bop = new compare_ev<EdgeOp_plus<int>, evplus_factor<int>,
                    le_evplus<int> > (a,b,c);
                break;

            case edge_type::LONG:
                bop = new compare_ev<EdgeOp_plus<long>, evplus_factor<long>,
                    le_evplus<long> > (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
    }
    return LE_cache.add( bop );
}

void MEDDLY::LESS_THAN_EQUAL_init()
{
    LE_cache.reset("LessEqual");
}

void MEDDLY::LESS_THAN_EQUAL_done()
{
    ASSERT(__FILE__, __LINE__, LE_cache.isEmpty());
}


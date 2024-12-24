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
#include "arith.h"

#include "../ops_builtin.h" // for COPY
#include "../oper_item.h"
#include "../oper_binary.h"
#include "../oper_unary.h"
#include "../ct_vector.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

//
// Operation instance caches
//
namespace MEDDLY {
    binary_list MAXIMUM_cache;
    binary_list MINIMUM_cache;
};

// #define USE_MT_OPS

// ******************************************************************
// ******************************************************************
// ******************************************************************
// ***                                                            ***
// ***                                                            ***
// ***                 Multi-terminal  operations                 ***
// ***                                                            ***
// ***                                                            ***
// ******************************************************************
// ******************************************************************
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                                                                *
// *                         arith_mt class                         *
// *                                                                *
// *                                                                *
// ******************************************************************

/*
    Single template class for all multi-terminal,
    element-wise arithmetic operations.

    Required methods for ATYPE classes (should be inlined):

        /// Get the operation name, for display purposes
        static const char* name();

        /// Does the operation commute?
        static bool commutes();

        /// Is the second operand the identity element?
        /// I.e., does x OP b = x, for all x?
        static bool isIdentityValue(node_handle b);

        ///
        /// Does x OP x = x for all x?
        ///
        static bool isIdempotent();

        ///
        /// Does there exist a constant c such that
        ///    x OP x = c, for all x?
        /// If yes, set the value of terminal c.
        ///
        static bool isConstOnSameOpnds(const forest*, node_handle &c)

        /// Apply the operation on terminals a and b,
        /// to obtain the result terminal c.
        static void apply(const forest* fa, node_handle a,
                          const forest* fb, node_handle b,
                          const forest* fc, node_handle &c);

 */

namespace MEDDLY {
    template <class ATYPE>
    class arith_mt : public binary_operation {
        public:
            arith_mt(forest* arg1, forest* arg2, forest* res);
            virtual ~arith_mt();

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
            unary_operation* copy_arg1res;
            unary_operation* copy_arg2res;
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
            bool forced_by_levels;
    };
};

// ******************************************************************

template <class ATYPE>
MEDDLY::arith_mt<ATYPE>::arith_mt(forest* arg1, forest* arg2,
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
    // Quick copying, even across forests, for terminal cases :)
    //
    copy_arg1res = COPY(arg1, res);
    copy_arg2res = COPY(arg2, res);

    //
    // Do we need to recurse by levels and store level info in the CT?
    // YES, if either forest is identity-reduced
    //      (even if the other is quasi-reduced, because we can
    //       still jump to terminal 0)
    //
    forced_by_levels = arg1->isIdentityReduced() || arg2->isIdentityReduced();

    // Build compute table key and result types.
    // If we recurse by levels, then we need the level as part of the key.
    ct = new ct_entry_type(ATYPE::name());
    if (forced_by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();
}

template <class ATYPE>
MEDDLY::arith_mt<ATYPE>::~arith_mt()
{
    ct->markForDestroy();
}

template <class ATYPE>
void MEDDLY::arith_mt<ATYPE>::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << ATYPE::name() << " #" << top_count << " begin\n";
#endif
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());

    _compute(L, in, ap, bp, cp);

#ifdef TRACE
    out << ATYPE::name() << " #" << top_count << " end\n";
#endif
    cv.set();
}


template <class ATYPE>
void MEDDLY::arith_mt<ATYPE>::_compute(int L, unsigned in,
        node_handle A, node_handle B, node_handle &C)
{
    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************

    if ( ATYPE::isIdentityValue(B) )
    {
        //
        // Result is A
        //
        edge_value dummy;
        dummy.set();
        MEDDLY_DCASSERT(copy_arg1res);
        copy_arg1res->compute(L, in, dummy, A, dummy, C);
        MEDDLY_DCASSERT(dummy.isVoid());
        return;
    }

    if ( ATYPE::commutes() && ATYPE::isIdentityValue(A) )
    {
        //
        // Result is B
        //
        edge_value dummy;
        dummy.set();
        MEDDLY_DCASSERT(copy_arg2res);
        copy_arg2res->compute(L, in, dummy, B, dummy, C);
        MEDDLY_DCASSERT(dummy.isVoid());
        return;
    }

    if ( ATYPE::isIdempotent() && A == B && arg1F == arg2F )
    {
        //
        // Result is A
        //
        edge_value dummy;
        dummy.set();
        MEDDLY_DCASSERT(copy_arg1res);
        copy_arg1res->compute(L, in, dummy, A, dummy, C);
        MEDDLY_DCASSERT(dummy.isVoid());
        return;
    }

    if ( ATYPE::isConstOnSameOpnds(resF, C) && A == B && arg1F == arg2F )
    {
        //
        // Result is a constant C
        //
        C = resF->makeRedundantsTo(C, 0, L);
        return;
    }

    if ( arg1F->isTerminalNode(A) && arg2F->isTerminalNode(B) )
    {
        // A and B are both terminal nodes.
        // We can stop the recursion if L==0 or
        // we're not forced to proceed by levels.
        if ((0==L) || (!forced_by_levels))
        {
            ATYPE::apply(arg1F, A, arg2F, B, resF, C);
            C = resF->makeRedundantsTo(C, 0, L);
            return;
        }
    }

    //
    // Reorder A and B if the operation commutes
    // and A,B are from the same forest
    //
    if ( ATYPE::commutes() && arg1F == arg2F )
    {
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
    out << ATYPE::name() << " arith_mt::compute(" << L << ", " << in << ", "
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
    ct_vector key(forced_by_levels ? 3 : 2);
    ct_vector res(1);
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

    unpacked_node* Au;
    if (Alevel != Clevel) {
        if (arg1F->isIdentityReduced() && Clevel<0) {
            MEDDLY_DCASSERT(Clevel == L);
            Au = unpacked_node::newIdentity(arg1F, Clevel, in, A, FULL_ONLY);
            MEDDLY_DCASSERT(Au->wasIdentity());
        } else {
            Au = unpacked_node::newRedundant(arg1F, Clevel, A, FULL_ONLY);
            MEDDLY_DCASSERT(!Au->wasIdentity());
        }
    } else {
        Au = arg1F->newUnpacked(A, FULL_ONLY);
        MEDDLY_DCASSERT(!Au->wasIdentity());
    }

    unpacked_node* Bu;
    if (Blevel != Clevel) {
        if (arg2F->isIdentityReduced() && Clevel<0) {
            MEDDLY_DCASSERT(Clevel == L);
            Bu = unpacked_node::newIdentity(arg2F, Clevel, in, B, FULL_ONLY);
            MEDDLY_DCASSERT(Bu->wasIdentity());
        } else {
            Bu = unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY);
            MEDDLY_DCASSERT(!Bu->wasIdentity());
        }
    } else {
        Bu = arg2F->newUnpacked(B, FULL_ONLY);
        MEDDLY_DCASSERT(!Bu->wasIdentity());
    }

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());
    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

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
    out << ATYPE::name() << " arith_mt::compute(" << L << ", " << in << ", "
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
    MEDDLY_DCASSERT(dummy.isVoid());
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
// *                        mt_maximum class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct mt_maximum {
        inline static const char* name() {
            return "max";
        }

        inline static bool commutes() {
            return true;
        }

        inline static bool isIdentityValue(node_handle b) {
            // only true for negative infinity
            return false;
        }

        inline static bool isIdempotent() {
            return true;
            // MAX (x, x) = x
        }

        inline static bool isConstOnSameOpnds(const forest*, node_handle &)
        {
            return false;
        }

        inline static void apply(const forest* fa, node_handle a,
                const forest* fb, node_handle b,
                const forest* fc, node_handle &c)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            c = fc->handleForValue( MAX(av, bv) );
        }
    };
};

// ******************************************************************
// *                                                                *
// *                        mt_minimum class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct mt_minimum {
        inline static const char* name() {
            return "min";
        }

        inline static bool commutes() {
            return true;
        }

        inline static bool isIdentityValue(node_handle b) {
            // only true for positive infinity
            return false;
        }

        inline static bool isIdempotent() {
            return true;
            // MIN (x, x) = x
        }

        inline static bool isConstOnSameOpnds(const forest*, node_handle &)
        {
            return false;
        }

        inline static void apply(const forest* fa, node_handle a,
                const forest* fb, node_handle b,
                const forest* fc, node_handle &c)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            c = fc->handleForValue( MIN(av, bv) );
        }
    };
};



// ******************************************************************
// ******************************************************************
// ******************************************************************
// ***                                                            ***
// ***                                                            ***
// ***                   Edge-valued operations                   ***
// ***                                                            ***
// ***                                                            ***
// ******************************************************************
// ******************************************************************
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                                                                *
// *                       arith_compat class                       *
// *                                                                *
// *                                                                *
// ******************************************************************

/*
    Template class for "compatible" edge-valued,
    element-wise arithmetic operations.
    For an edge operator +, "compatible" means that the
    arithmetic operator OP satisfies
        (a + f()) OP (b + g()) = (a+b) + (f() OP g())

    Examples include edge operator + with OP + or -,
    edge operator * with OP * or /,
    and any multi-termianl operator.

    Required methods for ATYPE classes (should be inlined):

        /// Get the operation name, for display purposes
        static const char* name();

        /// Does the operation commute?
        static bool commutes();

        /// Is the second operand the identity element?
        /// I.e., does x OP b = x, for all x?
        static bool isIdentityValue(node_handle b);

        ///
        /// Does x OP x = x for all x?
        ///
        static bool isIdempotent();

        ///
        /// Does there exist a terminal node c such that
        ///    x OP x = <e, c>, for all x?
        /// Here, e is the identity edge value for the edge operator.
        /// If yes, set the value of terminal c.
        ///
        static bool isConstOnSameOpnds(const forest*, node_handle &c)

        /// Apply the operation on terminals a and b,
        /// to obtain the result terminal c.
        static void apply(const forest* fa, node_handle a,
                          const forest* fb, node_handle b,
                          const forest* fc, node_handle &c);

 */

namespace MEDDLY {
    template <class EOP, class ATYPE>
    class arith_compat : public binary_operation {
        public:
            arith_compat(forest* arg1, forest* arg2, forest* res);
            virtual ~arith_compat();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    const edge_value &bv, node_handle bp,
                    edge_value &cv, node_handle &cp);

        protected:
            void _compute(int L, unsigned in, node_handle A, node_handle B,
                    edge_value &cv, node_handle &C);

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
            unary_operation* copy_arg1res;
            unary_operation* copy_arg2res;
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
            bool forced_by_levels;
    };
};

// ******************************************************************

template <class EOP, class ATYPE>
MEDDLY::arith_compat<EOP, ATYPE>::arith_compat(forest* arg1, forest* arg2,
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
    // Quick copying, even across forests, for terminal cases :)
    //
    copy_arg1res = COPY(arg1, res);
    copy_arg2res = COPY(arg2, res);

    //
    // Do we need to recurse by levels and store level info in the CT?
    // YES, if either forest is identity-reduced
    //      (even if the other is quasi-reduced, because we can
    //       still jump to terminal 0)
    //
    forced_by_levels = arg1->isIdentityReduced() || arg2->isIdentityReduced();

    // Build compute table key and result types.
    // If we recurse by levels, then we need the level as part of the key.
    ct = new ct_entry_type(ATYPE::name());
    if (forced_by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();
}

template <class EOP, class ATYPE>
MEDDLY::arith_compat<EOP, ATYPE>::~arith_compat()
{
    ct->markForDestroy();
}

template <class EOP, class ATYPE>
void MEDDLY::arith_compat<EOP, ATYPE>::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << ATYPE::name() << " #" << top_count << " begin\n";
#endif

    _compute(L, in, ap, bp, cv, cp);
    EOP::accumulateOp(cv, av);
    EOP::accumulateOp(cv, bv);

#ifdef TRACE
    out << ATYPE::name() << " #" << top_count << " end\n";
#endif

}


template <class EOP, class ATYPE>
void MEDDLY::arith_compat<EOP, ATYPE>::_compute(int L, unsigned in,
        node_handle A, node_handle B, edge_value &cv, node_handle &C)
{
    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************

    if ( ATYPE::isIdentityValue(B) )
    {
        //
        // Result is A
        //
        EOP::clear(cv);
        MEDDLY_DCASSERT(copy_arg1res);
        copy_arg1res->compute(L, in, cv, A, cv, C);
        return;
    }

    if ( ATYPE::commutes() && ATYPE::isIdentityValue(A) )
    {
        //
        // Result is B
        //
        EOP::clear(cv);
        MEDDLY_DCASSERT(copy_arg2res);
        copy_arg2res->compute(L, in, cv, B, cv, C);
        return;
    }

    if ( ATYPE::isIdempotent() && A == B && arg1F == arg2F )
    {
        //
        // Result is A
        //
        EOP::clear(cv);
        MEDDLY_DCASSERT(copy_arg1res);
        copy_arg1res->compute(L, in, cv, A, cv, C);
        return;
    }

    if ( ATYPE::isConstOnSameOpnds(resF, C) && A == B && arg1F == arg2F )
    {
        //
        // Result is a constant C
        //
        EOP::clear(cv);
        C = resF->makeRedundantsTo(C, 0, L);
        return;
    }

    if ( arg1F->isTerminalNode(A) && arg2F->isTerminalNode(B) )
    {
        // A and B are both terminal nodes.
        // We can stop the recursion if L==0 or
        // we're not forced to proceed by levels.
        if ((0==L) || (!forced_by_levels))
        {
            ATYPE::apply(arg1F, A, arg2F, B, resF, C);
            C = resF->makeRedundantsTo(C, 0, L);
            EOP::clear(cv);
            return;
        }
    }

    //
    // Reorder A and B if the operation commutes
    // and A,B are from the same forest
    //
    if ( ATYPE::commutes() && arg1F == arg2F )
    {
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
    out << ATYPE::name() << " arith_mt::compute(" << L << ", " << in << ", "
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
    ct_vector key(forced_by_levels ? 3 : 2);
    ct_vector res(1);
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
        EOP::clear(cv);
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

    unpacked_node* Au;
    if (Alevel != Clevel) {
        if (arg1F->isIdentityReduced() && Clevel<0) {
            MEDDLY_DCASSERT(Clevel == L);
            Au = unpacked_node::newIdentity(arg1F, Clevel, in, A, FULL_ONLY);
            MEDDLY_DCASSERT(Au->wasIdentity());
        } else {
            Au = unpacked_node::newRedundant(arg1F, Clevel, A, FULL_ONLY);
            MEDDLY_DCASSERT(!Au->wasIdentity());
        }
    } else {
        Au = arg1F->newUnpacked(A, FULL_ONLY);
        MEDDLY_DCASSERT(!Au->wasIdentity());
    }

    unpacked_node* Bu;
    if (Blevel != Clevel) {
        if (arg2F->isIdentityReduced() && Clevel<0) {
            MEDDLY_DCASSERT(Clevel == L);
            Bu = unpacked_node::newIdentity(arg2F, Clevel, in, B, FULL_ONLY);
            MEDDLY_DCASSERT(Bu->wasIdentity());
        } else {
            Bu = unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY);
            MEDDLY_DCASSERT(!Bu->wasIdentity());
        }
    } else {
        Bu = arg2F->newUnpacked(B, FULL_ONLY);
        MEDDLY_DCASSERT(!Bu->wasIdentity());
    }

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());
    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

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
        edge_value x;
        _compute(Cnextlevel, i, Au->down(i), Bu->down(i), x, cd);
        if (EOP::hasEdgeValues())
        {
            EOP::accumulateOp(x, Au->edgeval(i));
            EOP::accumulateOp(x, Bu->edgeval(i));
            Cu->setFull(i, x, cd);
        }
        else
        {
            Cu->setFull(i, cd);
        }
    }

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << ATYPE::name() << " arith_mt::compute(" << L << ", " << in << ", "
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
    resF->createReducedNode(Cu, cv, C);
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
// *                                                                *
// *                           Front ends                           *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                            MAXIMUM                             *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::MAXIMUM(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  MAXIMUM_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
#ifdef USE_MT_OPS
        if (use_reals) {
            bop = new arith_mt<mt_maximum<float> > (a,b,c);
        } else {
            bop = new arith_mt<mt_maximum<long> > (a,b,c);
        }
#else
        if (use_reals) {
            bop = new arith_compat<EdgeOp_none, mt_maximum<float> > (a,b,c);
        } else {
            bop = new arith_compat<EdgeOp_none, mt_maximum<long> > (a,b,c);
        }
#endif
        return MAXIMUM_cache.add( bop );
    }


    // throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    return nullptr;
}

void MEDDLY::MAXIMUM_init()
{
    MAXIMUM_cache.reset("Maximum");
}

void MEDDLY::MAXIMUM_done()
{
    MEDDLY_DCASSERT(MAXIMUM_cache.isEmpty());
}

// ******************************************************************
// *                            MINIMUM                             *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::MINIMUM(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  MINIMUM_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
#ifdef USE_MT_OPS
        if (use_reals) {
            bop = new arith_mt<mt_minimum<float> > (a,b,c);
        } else {
            bop = new arith_mt<mt_minimum<long> > (a,b,c);
        }
#else
        if (use_reals) {
            bop = new arith_compat<EdgeOp_none, mt_minimum<float> > (a,b,c);
        } else {
            bop = new arith_compat<EdgeOp_none, mt_minimum<long> > (a,b,c);
        }
#endif
        return MINIMUM_cache.add( bop );
    }


    // throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    return nullptr;
}

void MEDDLY::MINIMUM_init()
{
    MINIMUM_cache.reset("Maximum");
}

void MEDDLY::MINIMUM_done()
{
    MEDDLY_DCASSERT(MINIMUM_cache.isEmpty());
}


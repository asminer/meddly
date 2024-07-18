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

// #define OLD_OPS

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
// *         Massive template operation for all comparisons         *
// *                                                                *
// *                                                                *
// ******************************************************************

/*
    Required methods for CTYPE classes (should be inlined):

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

// ******************************************************************
// *                                                                *
// *                        compare_op class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class CTYPE>
    class compare_op : public binary_operation {
        public:
            compare_op(forest* arg1, forest* arg2, forest* res);
            virtual ~compare_op();

            virtual void compute(const edge_value &av, node_handle ap,
                    const edge_value &bv, node_handle bp,
                    int L,
                    edge_value &cv, node_handle &cp);

        protected:
            node_handle _compute(node_handle A, node_handle B, int L);

            node_handle _compute_un(node_handle A, node_handle B, int L);
            node_handle _compute_pr(int in, node_handle A, node_handle B, int L);

        private:
            ct_entry_type* ct;
#ifdef TRACE
            ostream_output out;
#endif
            bool go_by_levels;
    };
};

// ******************************************************************

template <class CTYPE>
MEDDLY::compare_op<CTYPE>::compare_op(forest* arg1, forest* arg2,
        forest* res) : binary_operation(arg1, arg2, res)
#ifdef TRACE
      , out(std::cout)
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
    // Do we need to recurse by levels?
    // If either input is an identity-reduced relation, then YES.
    //
    go_by_levels = arg1->isIdentityReduced() || arg2->isIdentityReduced();

    // Build compute table key and result types.
    // If we recurse by levels, then we need the level as part of the key.
    ct = new ct_entry_type(CTYPE::name());
    if (go_by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();
}

template <class CTYPE>
MEDDLY::compare_op<CTYPE>::~compare_op()
{
    ct->markForDestroy();
}

template <class CTYPE>
void MEDDLY::compare_op<CTYPE>::compute(
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        int L,
        edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
    out << "********************************************************\n";
    out << "Starting top-level call to " << CTYPE::name() << " compute\n";
#endif
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
    cv.set();
    if (resF->isForRelations()) {
        cp = _compute_un(ap, bp, L);
    } else {
        cp = _compute(ap, bp, L);
    }
#ifdef TRACE
    out.indentation(0);
    out << "Finished top-level call to " << CTYPE::name() << " compute\n";
    out << "********************************************************\n";
#endif
}

//
// SET VERSION
//

template <class CTYPE>
MEDDLY::node_handle
MEDDLY::compare_op<CTYPE>::_compute(node_handle A, node_handle B, int L)
{
    MEDDLY_DCASSERT(!go_by_levels);
    MEDDLY_DCASSERT(L >= 0);

    //
    // Terminal cases
    //
    if (A == B && arg1F == arg2F) {
        terminal tt(CTYPE::isReflexive(), resF->getTerminalType());
        return resF->makeRedundantsTo(tt.getHandle(), 0, L);
    }

    if (arg1F->isTerminalNode(A) && arg2F->isTerminalNode(B)) {
        terminal tt( CTYPE::compare(arg1F, A, arg2F, B), resF->getTerminalType() );
        return resF->makeRedundantsTo(tt.getHandle(), 0, L);
    }

    //
    // Reorder A and B if the relation is symmetric,
    // and they're in the same forest.
    //

    if (CTYPE::isSymmetric() && arg1F == arg2F) {
        if (A > B) {
            SWAP(A, B);
        }
    }

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = MAX(Alevel, Blevel);

    MEDDLY_DCASSERT(Alevel <= L);
    MEDDLY_DCASSERT(Blevel <= L);

#ifdef TRACE
    out << CTYPE::name() << " compute("
        << A << ", " << B << ", " << L << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << " before chain\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(2);
    ct_vector res(1);
    key[0].setN(A);
    key[1].setN(B);
    if (ct->findCT(key, res)) {
        node_handle C = resF->makeRedundantsTo(
                resF->linkNode(res[0].getN()), Clevel, L)
            ;
#ifdef TRACE
        out << "\tCT hit " << res[0].getN() << "\n";
        out << "\tafter chain " << C << "\n";
#endif
        return C;
    }

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel == Clevel)
        ?   arg1F->newUnpacked(A, FULL_ONLY)
        :   unpacked_node::newRedundant(arg1F, Clevel, A, FULL_ONLY);

    unpacked_node* Bu = (Blevel == Clevel)
        ?   arg2F->newUnpacked(B, FULL_ONLY)
        :   unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY);

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    for (unsigned i=0; i<Au->getSize(); i++) {
        Cu->setFull(i, _compute(Au->down(i), Bu->down(i), Clevel-1) );
    }
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << CTYPE::name() << " compute(" << A << ", " << B << ", " << L
              << ") done\n";
    out << "  A: ";
    Au->show(out, true);
    out << "\n  B: ";
    Bu->show(out, true);
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Bu);
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C);
    MEDDLY_DCASSERT(dummy.isVoid());
#ifdef TRACE
    out << "reduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT
    //
    if (Clevel > 0) {
        res[0].setN(C);
        ct->addCT(key, res);
    }
    C = resF->makeRedundantsTo(C, Clevel, L);

#ifdef TRACE
    out << "chain to level " << L << " = " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif
    return C;
}

//
// RELATION VERSION: unprimed levels
//

template <class CTYPE>
MEDDLY::node_handle
MEDDLY::compare_op<CTYPE>::_compute_un(node_handle A, node_handle B, int L)
{
    MEDDLY_DCASSERT(L >= 0);

    //
    // Terminal cases
    //
    if (0==A && 0==B) {
        terminal tt(CTYPE::isReflexive(), resF->getTerminalType());
        return resF->makeRedundantsTo(tt.getHandle(), 0, L);
    }

    if (A == B && arg1F == arg2F) {
        terminal tt(CTYPE::isReflexive(), resF->getTerminalType());
        return resF->makeRedundantsTo(tt.getHandle(), 0, L);
    }

    if (arg1F->isTerminalNode(A) && arg2F->isTerminalNode(B)) {
        if (0==L || !go_by_levels) {
            terminal tt( CTYPE::compare(arg1F, A, arg2F, B), resF->getTerminalType() );
            return resF->makeRedundantsTo(tt.getHandle(), 0, L);
        }
    }

    //
    // Reorder A and B if the relation is symmetric,
    // and they're in the same forest.
    //

    if (CTYPE::isSymmetric() && arg1F == arg2F) {
        if (A > B) {
            SWAP(A, B);
        }
    }

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = go_by_levels
        ? L
        : MAX(ABS(Alevel), ABS(Blevel));

    MEDDLY_DCASSERT(Clevel>0);
    MEDDLY_DCASSERT(ABS(Alevel) <= L);
    MEDDLY_DCASSERT(ABS(Blevel) <= L);

#ifdef TRACE
    out << CTYPE::name() << " compute_un("
        << A << ", " << B << ", " << L << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << " before chain\n";
#endif

    //
    // Check compute table
    //
    ct_vector key( go_by_levels ? 3 : 2);
    ct_vector res(1);
    if (go_by_levels) {
        key[0].setI(L);
        key[1].setN(A);
        key[2].setN(B);
    } else {
        key[0].setN(A);
        key[1].setN(B);
    }
    if (ct->findCT(key, res)) {
        node_handle C = resF->makeRedundantsTo(
                resF->linkNode(res[0].getN()), Clevel, L)
        ;
#ifdef TRACE
        out << "\tCT hit " << res[0].getN() << "\n";
        out << "\tafter chain " << C << "\n";
#endif
        return C;
    }

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel == Clevel)
        ?   arg1F->newUnpacked(A, FULL_ONLY)
        :   unpacked_node::newRedundant(arg1F, Clevel, A, FULL_ONLY);

    unpacked_node* Bu = (Blevel == Clevel)
        ?   arg2F->newUnpacked(B, FULL_ONLY)
        :   unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY);

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    for (unsigned i=0; i<Au->getSize(); i++) {
        Cu->setFull(i,
            _compute_pr(int(i), Au->down(i), Bu->down(i), -Clevel)
        );
    }
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << CTYPE::name() << " compute_un(" << A << ", " << B << ", " << L
              << ") done\n";
    out << "  A: ";
    Au->show(out, true);
    out << "\n  B: ";
    Bu->show(out, true);
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Bu);
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C);
    MEDDLY_DCASSERT(dummy.isVoid());
#ifdef TRACE
    out << "reduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT
    //
    res[0].setN(C);
    ct->addCT(key, res);
    C = resF->makeRedundantsTo(C, Clevel, L);

#ifdef TRACE
    out << "chain to level " << L << " = " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif
    return C;
}


template <class CTYPE>
MEDDLY::node_handle
MEDDLY::compare_op<CTYPE>::_compute_pr(int in, node_handle A, node_handle B,
        int L)
{
    MEDDLY_DCASSERT(L < 0);

    //
    // Terminal cases
    //
    if (0==A && 0==B) {
        terminal tt(CTYPE::isReflexive(), resF->getTerminalType());
        return resF->makeRedundantsTo(tt.getHandle(), 0, L);
    }

    if (A == B && arg1F == arg2F) {
        terminal tt(CTYPE::isReflexive(), resF->getTerminalType());
        return resF->makeRedundantsTo(tt.getHandle(), 0, L);
    }

    //
    // Reorder A and B if the relation is symmetric,
    // and they're in the same forest.
    //

    if (CTYPE::isSymmetric() && arg1F == arg2F) {
        if (A > B) {
            SWAP(A, B);
        }
    }

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = L;

#ifdef TRACE
    out << CTYPE::name() << " compute_pr("
        << A << ", " << B << ", " << L << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << " before chain\n";
#endif

    //
    // Check compute table
    //
    /*
    ct_vector key( go_by_levels ? 3 : 2);
    ct_vector res(1);
    if (go_by_levels) {
        key[0].setI(L);
        key[1].setN(A);
        key[2].setN(B);
    } else {
        key[0].setN(A);
        key[1].setN(B);
    }
    if (ct->findCT(key, res)) {
        node_handle C = resF->makeRedundantsTo(
                resF->linkNode(res[0].getN()), Clevel, L)
        ;
#ifdef TRACE
        out << "\tCT hit " << res[0].getN() << "\n";
        out << "\tafter chain " << C << "\n";
#endif
        return C;
    }
    */

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel == Clevel)
        ?   arg1F->newUnpacked(A, FULL_ONLY)
        :   MXD_levels::isFullyLevel(arg1F, Clevel)
            ?   unpacked_node::newRedundant(arg1F, Clevel, A, FULL_ONLY)
            :   unpacked_node::newIdentity(arg1F, Clevel, in, A, FULL_ONLY);

    unpacked_node* Bu = (Blevel == Clevel)
        ?   arg2F->newUnpacked(B, FULL_ONLY)
        :   MXD_levels::isFullyLevel(arg2F, Clevel)
            ?   unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY)
            :   unpacked_node::newIdentity(arg2F, Clevel, in, B, FULL_ONLY);

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    unsigned nnz = 0;
    for (unsigned i=0; i<Au->getSize(); i++) {
        node_handle d = _compute_un(Au->down(i), Bu->down(i), -Clevel-1);
        Cu->setFull(i, d);
        if (d) ++nnz;
    }
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << CTYPE::name() << " compute_pr(" << A << ", " << B << ", " << L
              << ") done\n";
    out << "  A: ";
    Au->show(out, true);
    out << "\n  B: ";
    Bu->show(out, true);
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Bu);
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C, in);
    MEDDLY_DCASSERT(dummy.isVoid());
#ifdef TRACE
    out << "reduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT when it is safe to do so:
    //  (1) A is not an identity node
    //  (2) B is not an identity node
    //  (3) Result cannot possibly be an identity node
    //
    /*
    bool canSaveResult = !Au->wasIdentity() && !Bu->wasIdentity();
    if (resF->isIdentityReduced() && (1==nnz)) canSaveResult = false;

    if (canSaveResult) {
        res[0].setN(C);
        ct->addCT(key, res);
    } else {
        ct->noaddCT(key);
    }
    */

    return C;
}




// ******************************************************************
// *                                                                *
// *                         eq_templ class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct eq_templ {
        inline static const char* name() {
            return "==";
        }

        inline static bool isSymmetric() {
            return true;    // x == y is the same as y == x
        }

        inline static bool isReflexive() {
            return true;    // x == x is true
        }

        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av == bv );
        }
    };

};

// ******************************************************************
// *                                                                *
// *                         ne_templ class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct ne_templ {
        inline static const char* name() {
            return "!=";
        }

        inline static bool isSymmetric() {
            return true;    // x != y is the same as y != x
        }

        inline static bool isReflexive() {
            return false;   // x != x is false
        }

        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av != bv );
        }
    };
};

// ******************************************************************
// *                                                                *
// *                         gt_templ class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct gt_templ {
        inline static const char* name() {
            return ">";
        }

        inline static bool isSymmetric() {
            return false;   // x > y is not the same as y > x
        }

        inline static bool isReflexive() {
            return false;   // x > x is false
        }

        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av > bv );
        }
    };
};

// ******************************************************************
// *                                                                *
// *                         ge_templ class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct ge_templ {
        inline static const char* name() {
            return ">=";
        }

        inline static bool isSymmetric() {
            return false;   // x >= y is not the same as y >= x
        }

        inline static bool isReflexive() {
            return true;    // x >= x is true
        }

        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av >= bv );
        }
    };
};

// ******************************************************************
// *                                                                *
// *                         lt_templ class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct lt_templ {
        inline static const char* name() {
            return "<";
        }

        inline static bool isSymmetric() {
            return false;   // x < y is not the same as y < x
        }

        inline static bool isReflexive() {
            return false;   // x < x is false
        }

        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av < bv );
        }
    };
};

// ******************************************************************
// *                                                                *
// *                         le_templ class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct le_templ {
        inline static const char* name() {
            return "<=";
        }

        inline static bool isSymmetric() {
            return false;   // x <= y is not the same as y <= x
        }

        inline static bool isReflexive() {
            return true;    // x <= x is true
        }

        inline static bool compare(const forest* fa, node_handle a,
                const forest* fb, node_handle b)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            return ( av <= bv );
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
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            if (c->isForRelations())
                return EQUAL_cache.add(
                    new compare_op<eq_templ<float> > (a,b,c)
                );
            else
                return EQUAL_cache.add(
                    new compare_op<eq_templ<float> > (a,b,c)
                );
        } else {
            if (c->isForRelations())
                return EQUAL_cache.add(
                    new compare_op<eq_templ<long> > (a,b,c)
                );
            else
                return EQUAL_cache.add(
                    new compare_op<eq_templ<long> > (a,b,c)
                );
        }
    }

    /*
    if (c->getEdgeLabeling() == edge_labeling::EVTIMES) {
        return EQUAL_cache.add(new equal_evtimes(a, b, c));
    }
    */

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::EQUAL_init()
{
    EQUAL_cache.reset("Equal");
}

void MEDDLY::EQUAL_done()
{
    MEDDLY_DCASSERT(EQUAL_cache.isEmpty());
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
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            if (c->isForRelations())
                return NEQ_cache.add(
                    new compare_op<ne_templ<float> > (a,b,c)
                );
            else
                return NEQ_cache.add(
                    new compare_op<ne_templ<float> > (a,b,c)
                );
        } else {
            if (c->isForRelations())
                return NEQ_cache.add(
                    new compare_op<ne_templ<long> > (a,b,c)
                );
            else
                return NEQ_cache.add(
                    new compare_op<ne_templ<long> > (a,b,c)
                );
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::NOT_EQUAL_init()
{
    NEQ_cache.reset("Unequal");
}

void MEDDLY::NOT_EQUAL_done()
{
    MEDDLY_DCASSERT(NEQ_cache.isEmpty());
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
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            if (c->isForRelations())
                return GT_cache.add(
                    new compare_op<gt_templ<float> > (a,b,c)
                );
            else
                return GT_cache.add(
                    new compare_op<gt_templ<float> > (a,b,c)
                );
        } else {
            if (c->isForRelations())
                return GT_cache.add(
                    new compare_op<gt_templ<long> > (a,b,c)
                );
            else
                return GT_cache.add(
                    new compare_op<gt_templ<long> > (a,b,c)
                );
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::GREATER_THAN_init()
{
    GT_cache.reset("MoreThan");
}

void MEDDLY::GREATER_THAN_done()
{
    MEDDLY_DCASSERT(GT_cache.isEmpty());
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
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            if (c->isForRelations())
                return GE_cache.add(
                    new compare_op<ge_templ<float> > (a,b,c)
                );
            else
                return GE_cache.add(
                    new compare_op<ge_templ<float> > (a,b,c)
                );
        } else {
            if (c->isForRelations())
                return GE_cache.add(
                    new compare_op<ge_templ<long> > (a,b,c)
                );
            else
                return GE_cache.add(
                    new compare_op<ge_templ<long> > (a,b,c)
                );
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::GREATER_THAN_EQUAL_init()
{
    GE_cache.reset("MoreEqual");
}

void MEDDLY::GREATER_THAN_EQUAL_done()
{
    MEDDLY_DCASSERT(GE_cache.isEmpty());
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
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            if (c->isForRelations())
                return LT_cache.add(
                    new compare_op<lt_templ<float> > (a,b,c)
                );
            else
                return LT_cache.add(
                    new compare_op<lt_templ<float> > (a,b,c)
                );
        } else {
            if (c->isForRelations())
                return LT_cache.add(
                    new compare_op<lt_templ<long> > (a,b,c)
                );
            else
                return LT_cache.add(
                    new compare_op<lt_templ<long> > (a,b,c)
                );
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::LESS_THAN_init()
{
    LT_cache.reset("LessThan");
}

void MEDDLY::LESS_THAN_done()
{
    MEDDLY_DCASSERT(LT_cache.isEmpty());
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
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            if (c->isForRelations())
                return LE_cache.add(
                    new compare_op<le_templ<float> > (a,b,c)
                );
            else
                return LE_cache.add(
                    new compare_op<le_templ<float> > (a,b,c)
                );
        } else {
            if (c->isForRelations())
                return LE_cache.add(
                    new compare_op<le_templ<long> > (a,b,c)
                );
            else
                return LE_cache.add(
                    new compare_op<le_templ<long> > (a,b,c)
                );
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::LESS_THAN_EQUAL_init()
{
    LE_cache.reset("LessEqual");
}

void MEDDLY::LESS_THAN_EQUAL_done()
{
    MEDDLY_DCASSERT(LE_cache.isEmpty());
}


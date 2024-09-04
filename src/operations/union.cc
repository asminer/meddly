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
#include "union.h"
#include "apply_base.h" // remove this when we can

#include "../ops_builtin.h" // for COPY
#include "../oper_binary.h"
#include "../oper_unary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class union_mdd;
    class union_mxd;

    class union_mt;

    class union_min_evplus;
    class union_min_evplus_mxd;

    binary_list UNION_cache;
};

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

// #define NEW_UNION

// ******************************************************************
// *                                                                *
// *                         union_mt class                         *
// *                                                                *
// ******************************************************************

/**
    Union operation for multi-terminal forests.
*/

class MEDDLY::union_mt : public binary_operation {
    public:
        union_mt(forest* arg1, forest* arg2, forest* res);
        virtual ~union_mt();

        virtual void compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp);

    private:
        /*
            Recursively compute.

                @param  in  Incoming edge label. Used ONLY when we skip
                            a primed level in an identity reduced forest.
                            And if that happens, we CANNOT cache the result.

                @param  L   Level to process.
                            If by_levels is false, then this parameter
                            is ignored.

                @param  A   First argument node
                @param  B   Second argument node

                @param  creqlev     Level that topOpndLevel says result
                                    will be at. Actual level could be lower,
                                    for example if the answer is redundant.
                                    If by_levels is true, this will equal L.

                @param  C   On output: answer node handle.
                            C = A u B.

         */
        void _compute(unsigned in, int L, node_handle A, node_handle B,
                int &creqlev, node_handle &C);

    private:
        // inlined helpers

        inline int topOpndLevel(int alevel, int blevel, int L) const
        {
            if (by_levels) return L;
            if (resF->isForRelations()) {
                return MXD_levels::topLevel(alevel, blevel);
            } else {
                return MDD_levels::topLevel(alevel, blevel);
            }
        }

        inline unpacked_node* unpackEdge(forest* F, int atL, int in,
                node_handle A, int alevel, node_storage_flags fs)
        {
            MEDDLY_DCASSERT(A);
            MEDDLY_DCASSERT(F->getNodeLevel(A) == alevel);

            if (alevel == atL) {
                return F->newUnpacked(A, fs);
            }
            if (atL >= 0 || 0==A || F->isFullyReduced()) {
                return unpacked_node::newRedundant(F, atL, A, fs);
            }
            if (F->isIdentityReduced()) {
                return unpacked_node::newIdentity(F, atL, in, A, fs);
            }
            std::cout << "unpackEdge fall through\n";
            std::cout << "    atL: " << atL << "\n";
            std::cout << "    A: " << A << "\n";
            std::cout << "    a level: " << alevel << "\n";

            MEDDLY_DCASSERT(false);

            return nullptr;
        }

    private:
        ct_entry_type* ct;
        unary_operation* copy_arg1res;
        unary_operation* copy_arg2res;

#ifdef TRACE
        ostream_output out;
#endif

        bool by_levels;
        bool both_fully;
        bool both_identity;
};

// ******************************************************************

MEDDLY::union_mt::union_mt(forest* arg1, forest* arg2, forest* res)
    : binary_operation(arg1, arg2, res)
#ifdef TRACE
      , out(std::cout)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    //
    // Quick copying, even across forests, for terminal cases :)
    //
    copy_arg1res = COPY(arg1, res);
    copy_arg2res = COPY(arg2, res);

    //
    // Decide if we need to recurse by levels or not.
    // If yes, then we store the level in the compute table and
    // force the recursion to go by levels.
    // If no, we do not store the level in the compute table
    // and we recurse by the top level of the operands.
    //
    //  Cases that naturally have NO level skipping
    //
    //      quasi, quasi    :   by_levels = false
    //      quasi, fully    :   by_levels = false
    //      quasi, ident    :   by_levels = false
    //      fully, quasi    :   by_levels = false
    //      ident, quasi    :   by_levels = false
    //
    //  Cases with arbitrary level skipping allowed
    //
    //      fully, fully    :   by_levels = false
    //      ident, ident    :   by_levels = false
    //
    //  Cases that must proceed by levels, and store levels
    //  in the compute table entry
    //
    //      fully, ident    :   by_levels = true;
    //      ident, fully    :   by_levels = true;
    //

    by_levels =
        (arg1F->isFullyReduced() && arg2F->isIdentityReduced())
        ||
        (arg2F->isFullyReduced() && arg1F->isIdentityReduced())
    ;

    both_fully = arg1F->isFullyReduced() && arg2F->isFullyReduced();
    both_identity = arg1F->isIdentityReduced() && arg2F->isIdentityReduced();

    ct = new ct_entry_type("union");
    if (by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();
}


MEDDLY::union_mt::~union_mt()
{
    ct->markForDestroy();
}

void MEDDLY::union_mt::compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
#endif
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
    cv.set();
    int ablevel;
    _compute(0, L, ap, bp, ablevel, cp);
    if (both_fully) {
        cp = resF->makeRedundantsTo(cp, ablevel, L);
    }
    if (both_identity) {
        cp = resF->makeIdentitiesTo(cp, ablevel, L, -1);
    }
}

void MEDDLY::union_mt::_compute(unsigned in, int L,
        node_handle A, node_handle B, int &Clevel, node_handle &C)
{
    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    Clevel = topOpndLevel(Alevel, Blevel, L);

    //
    // Terminal case: both zero
    //
    if ((A == 0) && (B == 0)) {
        // This is correct, even if we're forced to go by levels.
        C = 0;
        return;
    }

    if (A == 0) {
        //
        // Result is B
        //
        edge_value bv, cv;
        bv.set();
        MEDDLY_DCASSERT(copy_arg2res);
        copy_arg2res->compute(bv, B, Blevel, cv, C);
        if (arg2F->isIdentityReduced()) {
            C = resF->makeIdentitiesTo(C, Blevel, Clevel, in);
        } else {
            C = resF->makeRedundantsTo(C, Blevel, Clevel);
        }
        return;
    }

    if ( (B == 0) || ((A==B)&&(arg1F==arg2F)) )
    {
        //
        // Result is A
        //
        edge_value av, cv;
        av.set();
        MEDDLY_DCASSERT(copy_arg1res);
        copy_arg1res->compute(av, A, Alevel, cv, C);
        if (arg1F->isIdentityReduced()) {
            C = resF->makeIdentitiesTo(C, Alevel, Clevel, in);
        } else {
            C = resF->makeRedundantsTo(C, Alevel, Clevel);
        }
        return;
    }

    //
    // Both terminal one; result is fully-fully one unless
    // both argument forests are identity reduced; then its I.
    //
    if (A<0 && B<0) {
        terminal tt(true, resF->getTerminalType());
        if (both_identity) {
            C = resF->makeIdentitiesTo(tt.getHandle(), 0, Clevel, in);
        } else {
            C = resF->makeRedundantsTo(tt.getHandle(), 0, Clevel);
        }
        return;
    }

    //
    // One argument is a fully-fully ONE.
    //
    if ( (A<0 && arg1F->isFullyReduced()) || (B<0 && arg2F->isFullyReduced()) )
    {
        terminal tt(true, resF->getTerminalType());
        C = resF->makeRedundantsTo(tt.getHandle(), 0, Clevel);
        return;
    }


    //
    // Reorder A and B if they commute (same forest)
    //
    if (arg1F == arg2F) {
        if (A > B) {
            SWAP(A, B);
        }
    }

#ifdef TRACE
    out << "union_mt compute(" << L << ", " << A << ", " << B << ")";
    if (by_levels) {
        out << " by levels\n";
    } else {
        out << " by topmost\n";
    }
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << "\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(by_levels ? 3 : 2);
    ct_vector res(1);
    if (by_levels) {
        key[0].setI(L);
        key[1].setN(A);
        key[2].setN(B);
    } else {
        key[0].setN(A);
        key[1].setN(B);
    }
    if (ct->findCT(key, res)) {
        C = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "\tCT hit " << res[0].getN() << "\n";
        out << "\tlevel is " << resF->getNodeLevel(C) << "\n";
#endif
        return;
    }


    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = unpackEdge(arg1F, Clevel, in, A, Alevel, FULL_ONLY);
    unpacked_node* Bu = unpackEdge(arg2F, Clevel, in, B, Blevel, FULL_ONLY);
    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());
    const bool depends_on_in = Au->wasIdentity() || Bu->wasIdentity();

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    const int nextlevel = resF->isForRelations()
            ? MXD_levels::downLevel(Clevel)
            : MDD_levels::downLevel(Clevel);

// HERE TBD

    /*
    for (unsigned i=0; i<Cu->getSize(); i++) {

        node_handle d = _compute(int(i), Au->down(i), Bu->down(i), nextlevel);
        Cu->setFull(i, d);
        if (d) ++nnzs;
    }



    // OLD BELOW HERE TBD


    unsigned nnzs = 0;
    for (unsigned i=0; i<Cu->getSize(); i++) {
        node_handle d = _compute(int(i), Au->down(i), Bu->down(i), nextlevel);
        Cu->setFull(i, d);
        if (d) ++nnzs;
    }
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "union_mt compute(" << in << ", "
        << A << ", " << B << ", " << L << ") done\n";
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
    node_handle C;
    resF->createReducedNode(Cu, dummy, C, in);
    MEDDLY_DCASSERT(dummy.isVoid());
#ifdef TRACE
    out << "reduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    if (depends_on_in) {
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

    return resF->makeRedundantsTo(C, Clevel, L);
    */
}

// ******************************************************************
// *                                                                *
// *                        union_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::union_mdd : public binary_operation {
    public:
        union_mdd(forest* arg1, forest* arg2, forest* res);
        virtual ~union_mdd();

        virtual void compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        node_handle _compute(node_handle A, node_handle B, int L);

    private:
        ct_entry_type* ct;
};

// ******************************************************************

MEDDLY::union_mdd::union_mdd(forest* arg1, forest* arg2, forest* res)
  : binary_operation(arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("union");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    ct->setFixed(arg1, arg2);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::union_mdd::~union_mdd()
{
    ct->markForDestroy();
}

void MEDDLY::union_mdd::compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
    cv.set();
    cp = _compute(ap, bp, L);
}

MEDDLY::node_handle
MEDDLY::union_mdd::_compute(node_handle A, node_handle B, int L)
{
    //
    // Check terminals
    //
    if (A < 0 || B < 0) {
        terminal tt(true);
        return resF->makeRedundantsTo(tt.getHandle(), 0, L);
    }

    if (A == 0) {
        if (B == 0) {
            return 0;
        }

        //
        // Return B if we can
        //
        if (arg2F == resF) {
            return resF->linkNode(B);
            // Don't need to make redundant chain b/c same forest
        }
    } // zero A

    if (B == 0) {
        //
        // Return A if we can
        //
        MEDDLY_DCASSERT(A);
        if (arg1F == resF) {
            return resF->linkNode(A);
        }
    } // zero B

    if (A == B) {
        if ((arg1F == arg2F) && (arg1F == resF)) {
            return resF->linkNode(A);
        }
    }

    //
    // Reorder A and B if they commute (same forest)
    //

    if (arg1F == arg2F) {
        if (A > B) {
            SWAP(A, B);
        }
    }

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = resF->isForRelations()
                        ? MXD_levels::topLevel(Alevel, Blevel)
                        : MDD_levels::topLevel(Alevel, Blevel);

    //
    // Check compute table
    //
    ct_vector key(2);
    ct_vector res(1);
    key[0].setN(A);
    key[1].setN(B);
    if (ct->findCT(key, res)) {
        return resF->makeRedundantsTo(resF->linkNode(res[0].getN()), Clevel, L);
    }

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   unpacked_node::newRedundant(arg1F, Clevel, A, FULL_ONLY)
        :   arg1F->newUnpacked(A, FULL_ONLY);

    unpacked_node* Bu = (Blevel != Clevel)
        ?   unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY)
        :   arg2F->newUnpacked(B, FULL_ONLY);

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i, _compute(Au->down(i), Bu->down(i), Clevel-1));
    }

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Bu);
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C);
    MEDDLY_DCASSERT(dummy.isVoid());

    //
    // Save result in CT
    //
    res[0].setN(C);
    ct->addCT(key, res);

    return resF->makeRedundantsTo(C, Clevel, L);
}


// ******************************************************************
// *                                                                *
// *                        union_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::union_mxd : public binary_operation {
    public:
        union_mxd(forest* arg1, forest* arg2, forest* res);
        virtual ~union_mxd();

        virtual void compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        inline node_handle makeChainTo(node_handle p, int K, int L)
        {
            if (go_by_levels) return p;
            if (identity_chains) {
                MEDDLY_DCASSERT(!redundant_chains);
                return resF->makeIdentitiesTo(p, K, L, -1);
            }
            if (redundant_chains) {
                return resF->makeRedundantsTo(p, K, L);
            }
            return p;
        }

        inline unpacked_node* patternNode(forest* F, int L, int in,
                node_handle A, node_storage_flags fs)
        {
            MEDDLY_DCASSERT(L<0);
            if (0==A || F->isFullyReduced()) {
                return unpacked_node::newRedundant(F, L, A, fs);
            }
            if (F->isIdentityReduced()) {
                return unpacked_node::newIdentity(F, L, in, A, fs);
            }
            std::cout << "L: " << L << "\n";
            std::cout << "A: " << A << "\n";
            std::cout << "A.level: " << F->getNodeLevel(A) << "\n";
            MEDDLY_DCASSERT(false);
        }

    protected:
        node_handle _compute(node_handle A, node_handle B, int L);
        node_handle _compute_primed(int in, node_handle A, node_handle B,
                int L);


    private:
        ct_entry_type* ct;

        bool identity_chains;
        bool redundant_chains;

        bool go_by_levels;
};

// ******************************************************************

MEDDLY::union_mxd::union_mxd(forest* arg1, forest* arg2, forest* res)
  : binary_operation(arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    //
    // Determine if we can 'skip levels', and how to add chains
    // of nodes when needed.
    //
    identity_chains = arg1->isIdentityReduced() || arg2->isIdentityReduced();
    redundant_chains = arg1->isFullyReduced() || arg2->isFullyReduced();
    go_by_levels = identity_chains && redundant_chains;

    //
    // Possible input forests:
    //
    //      quasi, quasi    :   no level skipping will occur
    //      quasi, fully    :   (same)
    //      quasi, ident    :   (same)
    //      fully, quasi    :   (same)
    //      ident, quasi    :   (same)
    //
    //      fully, fully    :   can skip levels, build redundant chains
    //      ident, ident    :   can skip levels, build identity chains
    //
    //      fully, ident    :
    //      ident, fully    :   must go by levels to handle unioning
    //                          of skipped redundant/identity nodes.
    //                          compute table entries add level info.
    //                          go_by_levels will be true only for
    //                          these two cases.
    //

    ct = new ct_entry_type("union_mxd");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    if (go_by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::union_mxd::~union_mxd()
{
    ct->markForDestroy();
}

void MEDDLY::union_mxd::compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
    cv.set();
    cp = _compute(ap, bp, L);
}

MEDDLY::node_handle
MEDDLY::union_mxd::_compute(node_handle A, node_handle B, int L)
{
    MEDDLY_DCASSERT(L>=0);

    //
    // Check terminal cases
    //

    if (A == 0) {
        // A is empty set
        if (B == 0) {
            return 0;
        }
        if (B < 0) {
            terminal tt(true);
            if (arg2F->isIdentityReduced()) {
                return resF->makeIdentitiesTo(tt.getHandle(), 0, L, -1);
            } else {
                return resF->makeRedundantsTo(tt.getHandle(), 0, L);
            }
        }

        //
        // Return B if we can
        //
        if (arg2F == resF) {
            return resF->linkNode(B);
            // Don't need to make redundant chain b/c same forest
        }
    } // zero A

    if (B == 0) {
        // B is empty set
        if (A < 0) {
            terminal tt(true);
            if (arg1F->isIdentityReduced()) {
                return resF->makeIdentitiesTo(tt.getHandle(), 0, L, -1);
            } else {
                return resF->makeRedundantsTo(tt.getHandle(), 0, L);
            }
        }
        //
        // Return A if we can
        //
        MEDDLY_DCASSERT(A);
        if (arg1F == resF) {
            return resF->linkNode(A);
        }
    } // zero B

    if (A == B) {
        if ((arg1F == arg2F) && (arg1F == resF)) {
            return resF->linkNode(A);
        }
    }

    // Both terminal one
    if (A < 0 && B < 0) {
        // if both are identity, return an identity chain to 1
        // otherwise, return a redundant chain to 1
        // (if either one is quasi, the chain length will be 0)
        terminal tt(true);
        if (arg1F->isIdentityReduced() && arg2F->isIdentityReduced()) {
            return resF->makeIdentitiesTo(tt.getHandle(), 0, L, -1);
        } else {
            return resF->makeRedundantsTo(tt.getHandle(), 0, L);
        }
    }
    // Just A is terminal one
    if (A < 0) {
        if (arg1F->isFullyReduced()) {
            terminal tt(true);
            return resF->makeRedundantsTo(tt.getHandle(), 0, L);
        }
    }
    // Just B is terminal one
    if (B < 0) {
        if (arg2F->isFullyReduced()) {
            terminal tt(true);
            return resF->makeRedundantsTo(tt.getHandle(), 0, L);
        }
    }


    //
    // Reorder A and B if they commute (same forest)
    //

    if (arg1F == arg2F) {
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

#ifdef TRACE
    std::cout << "union_mxd::_compute(" << A << ", " << B << ", " << L << ")\n";
    std::cout << "\t" << A << " level " << Alevel << "\n";
    std::cout << "\t" << B << " level " << Blevel << "\n";
    std::cout << "\tresult level " << Clevel << " before chain\n";
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
        node_handle C = makeChainTo(resF->linkNode(res[0].getN()), Clevel, L);
#ifdef TRACE
        std::cout << "\tCT hit " << res[0].getN() << "\n";
        std::cout << "\tafter chain " << C << "\n";
#endif
        return C;
    }

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   unpacked_node::newRedundant(arg1F, Clevel, A, FULL_ONLY)
        :   arg1F->newUnpacked(A, FULL_ONLY);

    unpacked_node* Bu = (Blevel != Clevel)
        ?   unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY)
        :   arg2F->newUnpacked(B, FULL_ONLY);

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i, _compute_primed(i, Au->down(i), Bu->down(i),
                    MXD_levels::downLevel(Clevel)));
    }
#ifdef TRACE
    std::cout << "union_mxd::_compute(" << A << ", " << B << ", " << L
              << ") = " << "\n\t";
    ostream_output out(std::cout);
    Cu->show(out, true);
    std::cout << "\n";
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
    std::cout << "\treduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    std::cout << "\n";
#endif

    //
    // Save result in CT
    //
    res[0].setN(C);
    ct->addCT(key, res);
    C = makeChainTo(C, Clevel, L);

#ifdef TRACE
    std::cout << "\tchain to level " << L << " = " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    std::cout << "\n";
#endif
    return C;
}

MEDDLY::node_handle
MEDDLY::union_mxd::_compute_primed(int in, node_handle A, node_handle B,
        const int Clevel)
{
    MEDDLY_DCASSERT(Clevel<0);

    //
    // Terminal cases
    //

    if (A == 0) {
        if (B == 0) {
            return 0;
        }

        //
        // Return B if we can
        //
        if (arg2F == resF) {
            return resF->linkNode(B);
            // Don't need to make redundant chain b/c same forest
        }
    } // zero A

    if (B == 0) {
        //
        // Return A if we can
        //
        MEDDLY_DCASSERT(A);
        if (arg1F == resF) {
            return resF->linkNode(A);
        }
    } // zero B

    if (A == B) {
        if ((arg1F == arg2F) && (arg1F == resF)) {
            return resF->linkNode(A);
        }
    }

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
#ifdef TRACE
    std::cout << "union_mxd::_compute_primed(" << in << ", " << A << ", " << B << ", " << Clevel << ")\n";
    std::cout << "\t" << A << " level " << Alevel << "\n";
    std::cout << "\t" << B << " level " << Blevel << "\n";
#endif

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   patternNode(arg1F, Clevel, in, A, FULL_ONLY)
        :   arg1F->newUnpacked(A, FULL_ONLY);

    unpacked_node* Bu = (Blevel != Clevel)
        ?   patternNode(arg2F, Clevel, in, B, FULL_ONLY)
        :   arg2F->newUnpacked(B, FULL_ONLY);

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i,
            _compute(Au->down(i), Bu->down(i),
                MXD_levels::downLevel(Cu->getLevel()))
        );
    }

#ifdef TRACE
    ostream_output out(std::cout);
    std::cout << "\t";
    Cu->show(out, true);
    std::cout << "\n";
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
    std::cout << "union_mxd::_compute_primed(" << in << ", " << A << ", "
              << B << ", " << Clevel << ") = " << C << "\n\t";
    resF->showNode(out, C, SHOW_DETAILS);
    out.put('\n');
#endif
    return C;
}


// ******************************************************************
// *                                                                *
// *                    union_min_evplus  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::union_min_evplus : public generic_binary_evplus {
    public:
        union_min_evplus(forest* arg1, forest* arg2, forest* res);

        virtual ct_entry_key* findResult(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle &c);
        virtual void saveResult(ct_entry_key* key, long aev, node_handle a,
            long bev, node_handle b, long cev, node_handle c);

        virtual bool checkTerminals(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle& c);
};

MEDDLY::union_min_evplus::union_min_evplus(forest* a, forest* b, forest* c)
  : generic_binary_evplus(UNION_cache, a, b, c)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

MEDDLY::ct_entry_key* MEDDLY::union_min_evplus::findResult(long aev, node_handle a,
  long bev, node_handle b, long& cev, node_handle &c)
{
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  if (canCommute() && a > b) {
    CTsrch->writeL(0);
    CTsrch->writeN(b);
    CTsrch->writeL(aev - bev);
    CTsrch->writeN(a);
  } else {
    CTsrch->writeL(0);
    CTsrch->writeN(a);
    CTsrch->writeL(bev - aev);
    CTsrch->writeN(b);
  }
  CT0->find(CTsrch, CTresult[0]);
  if (!CTresult[0]) return CTsrch;
  cev = CTresult[0].readL();
  MEDDLY_DCASSERT(cev == 0);
  c = resF->linkNode(CTresult[0].readN());
  if (c != 0) {
    cev = MIN(aev, bev);
  }
  CT0->recycle(CTsrch);
  return 0;
}

void MEDDLY::union_min_evplus::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  MEDDLY_DCASSERT(c == 0 || cev == MIN(aev, bev));
  CTresult[0].reset();
  CTresult[0].writeL(0);   //   Why always 0?
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

bool MEDDLY::union_min_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == 0) {
    if (b == 0) {
      cev = 0;
      c = 0;
      return true;
    }
    else if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg1F->isTerminalNode(a) && aev <= bev) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg2F->isTerminalNode(b) && bev <= aev) {
    if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (a == b) {
    if (arg1F == arg2F && arg2F == resF) {
      cev = MIN(aev, bev);
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                  union_min_evplus_mxd  class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::union_min_evplus_mxd : public generic_binary_evplus_mxd {
    public:
        union_min_evplus_mxd(forest* arg1, forest* arg2, forest* res);

        virtual ct_entry_key* findResult(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle &c);
        virtual void saveResult(ct_entry_key* key, long aev, node_handle a,
            long bev, node_handle b, long cev, node_handle c);

        virtual bool checkTerminals(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle& c);
};

MEDDLY::union_min_evplus_mxd::union_min_evplus_mxd(
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_evplus_mxd(UNION_cache, arg1, arg2, res)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

MEDDLY::ct_entry_key* MEDDLY::union_min_evplus_mxd::findResult(long aev, node_handle a,
  long bev, node_handle b, long& cev, node_handle &c)
{
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  if (canCommute() && a > b) {
    CTsrch->writeL(0);
    CTsrch->writeN(b);
    CTsrch->writeL(aev - bev);
    CTsrch->writeN(a);
  } else {
    CTsrch->writeL(0);
    CTsrch->writeN(a);
    CTsrch->writeL(bev - aev);
    CTsrch->writeN(b);
  }
  CT0->find(CTsrch, CTresult[0]);
  if (!CTresult[0]) return CTsrch;
  cev = CTresult[0].readL();
  MEDDLY_DCASSERT(cev == 0);
  c = resF->linkNode(CTresult[0].readN());
  if (c != 0) {
    cev = MIN(aev, bev);
  }
  CT0->recycle(CTsrch);
  return 0;
}

void MEDDLY::union_min_evplus_mxd::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  MEDDLY_DCASSERT(c == 0 || cev == MIN(aev, bev));
  CTresult[0].reset();
  CTresult[0].writeL(0);   // why always 0?
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

bool MEDDLY::union_min_evplus_mxd::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == 0) {
    if (b == 0) {
      cev = 0;
      b = 0;
      return true;
    }
    else if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg1F->isTerminalNode(a) && aev <= bev) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg2F->isTerminalNode(b) && bev <= aev) {
    if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (a == b) {
    if (arg1F == arg2F && arg2F == resF) {
      cev = MIN(aev, bev);
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation*
MEDDLY::UNION(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }
    binary_operation* bop =  UNION_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
#ifdef NEW_UNION
        return UNION_cache.add(new union_mt(a, b, c));
#else
        if (c->isForRelations()) {
            return UNION_cache.add(new union_mxd(a, b, c));
        } else {
            return UNION_cache.add(new union_mdd(a, b, c));
        }
#endif
    }

    if (c->getEdgeLabeling() == edge_labeling::EVPLUS) {
        if (c->isForRelations()) {
            return UNION_cache.add(new union_min_evplus_mxd(a, b, c));
        } else {
            return UNION_cache.add(new union_min_evplus(a, b, c));
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::UNION_init()
{
    UNION_cache.reset("Union");
}

void MEDDLY::UNION_done()
{
    MEDDLY_DCASSERT(UNION_cache.isEmpty());
}


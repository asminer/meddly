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

#include "../ops_builtin.h" // for COPY
#include "../oper_binary.h"
#include "../oper_unary.h"
#include "../ct_vector.h"
#include "../forest_levels.h"

namespace MEDDLY {
    class union_mt;
    class UNION_factory;
};

// #define TRACE
#define USE_PRIMED_CACHE

#ifdef TRACE
#include "../operators.h"
#endif

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

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

    private:
        void _compute(int L, unsigned in,
                node_handle A, node_handle B, node_handle &C);

    private:
        // inlined helpers

        inline int topLevelOf(int L, int alevel, int blevel) const
        {
            if (by_levels) return L;
            if (resF->isForRelations()) {
                if (both_identity) {
                    const int t = MXD_levels::topUnprimed(alevel, blevel);
                    return (t == MXD_levels::upLevel(L)) ? L : t;
                } else {
                    return MXD_levels::topLevel(alevel, blevel);
                }
            } else {
                return MDD_levels::topLevel(alevel, blevel);
            }
        }

        inline void chainToLevel(node_handle &C, int Clevel, int L, unsigned in)
        {
            //
            // Add nodes from Clevel to L
            //
            if (both_identity) {
#ifdef TRACE
                out << "I chain to " << C << ", levels " << Clevel
                    << " to " << L << "\n";
#endif
                C = resF->makeIdentitiesTo(C, Clevel, L, in);
            } else {
#ifdef TRACE
                out << "X chain to " << C << ", levels " << Clevel
                    << " to " << L << "\n";
#endif
                C = resF->makeRedundantsTo(C, Clevel, L);
            }
        }

    private:
        ct_entry_type* ct;
#ifdef USE_PRIMED_CACHE
        ct_entry_type* ct_primed;
#endif
        unary_operation* copy_arg1res;
        unary_operation* copy_arg2res;

#ifdef TRACE
        ostream_output out;
        unsigned top_count;
#endif

        bool by_levels;
        bool forced_by_levels;
        bool both_fully;
        bool both_identity;
};

// ******************************************************************

MEDDLY::union_mt::union_mt(forest* arg1, forest* arg2, forest* res)
    : binary_operation(arg1, arg2, res)
#ifdef TRACE
      , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    //
    // Quick copying, even across forests, for terminal cases :)
    //
    copy_arg1res = build(COPY, arg1, res);
    copy_arg2res = build(COPY, arg2, res);

    //
    // How to handle different reductions
    //
    //      quasi v quasi   :   no level skipping will occur
    //      quasi v fully   :
    //      quasi v ident   :
    //      fully v quasi   :
    //      ident v quasi   :
    //
    //      ident v ident   :   identity pattern: skip by unprimed
    //
    //                          [p 0 0]   [q 0 0]   [ pvq  0   0  ]
    //                          [0 p 0] v [0 q 0] = [  0  pvq  0  ]
    //                          [0 0 p]   [0 0 q]   [  0   0  pvq ]
    //
    //      ident v fully   :   force by levels, and store levels in
    //      fully v ident   :   the compute table entry
    //
    //                          [p p p]   [q 0 0]   [ pvq  p   p  ]
    //                          [p p p] v [0 q 0] = [  p  pvq  p  ]
    //                          [p p p]   [0 0 q]   [  p   p  pvq ]
    //
    //      fully v fully   :   fully pattern: skip by top level
    //
    //                          [p p p]   [q q q]   [ pvq pvq pvq ]
    //                          [p p p] v [q q q] = [ pvq pvq pvq ]
    //                          [p p p]   [q q q]   [ pvq pvq pvq ]
    //
    //

    both_fully = arg1F->isFullyReduced() && arg2F->isFullyReduced();
    both_identity = arg1F->isIdentityReduced() && arg2F->isIdentityReduced();

    by_levels = !(both_fully || both_identity);
    forced_by_levels =
        (arg1F->isFullyReduced() && arg2F->isIdentityReduced())
        ||
        (arg2F->isFullyReduced() && arg1F->isIdentityReduced())
    ;

    ct = new ct_entry_type("union");
    if (forced_by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();

#ifdef USE_PRIMED_CACHE
    //
    // If we're skipping by unprimed levels,
    // keep a second CT for the primed level computations
    // we are able to save. This is needed to differentiate
    // the unprimed level result and the primed level result.
    ct_primed = nullptr;
    if (both_identity) {
        ct_primed = new ct_entry_type("union_pr");
        ct_primed->setFixed(arg1, arg2);
        ct_primed->setResult(res);
        ct_primed->doneBuilding();
    }
#endif // USE_PRIMED_CACHE
}

MEDDLY::union_mt::~union_mt()
{
    ct->markForDestroy();
#ifdef USE_PRIMED_CACHE
    if (ct_primed) ct_primed->markForDestroy();
#endif
}

void MEDDLY::union_mt::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << "Union #" << top_count << " begin\n";
#endif
    _compute(L, in, ap, bp, cp);
#ifdef TRACE
    out << "Union #" << top_count << " end\n";
#endif
    cv.set();
}

void MEDDLY::union_mt::_compute(int L, unsigned in,
        node_handle A, node_handle B, node_handle &C)
{
    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************
    if (0==A && 0==B) {
        // This is correct, even if we're forced to go by levels.
        C = 0;
        return;
    }

    if (0==A) {
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

    if ( (0 == B) || ((A==B)&&(arg1F==arg2F)) )
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

    //
    // Both terminal one; result is fully-fully one unless
    // both argument forests are identity reduced; then it is I.
    //
    if (arg1F->isTerminalNode(A) && arg2F->isTerminalNode(B)) {
        terminal tt(true, resF->getTerminalType());
        if (both_identity) {
            C = resF->makeIdentitiesTo(tt.getHandle(), 0, L, in);
        } else {
            C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
        }
        return;
    }

    //
    // One argument is a fully-fully ONE.
    // Result is fully-fully one.
    //
    if ( (arg1F->isTerminalNode(A) && arg1F->isFullyReduced())
        || (arg2F->isTerminalNode(B) && arg2F->isFullyReduced()) )
    {
        terminal tt(true, resF->getTerminalType());
        C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
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
    const bool useCT = !both_identity || Clevel>0;
#ifdef USE_PRIMED_CACHE
    const bool useCTpr = !useCT && (L<0) &&
        (Alevel == L) && (Blevel == L) && (Clevel == L);
    MEDDLY_DCASSERT(!useCTpr || ct_primed);
#endif

#ifdef TRACE
    out << "union_mt::_compute(" << L << ", " << in << ", " << A << ", "
        << B << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << "\n";
#endif

    // **************************************************************
    //
    // Check the compute table or primed compute table
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

    if (useCT && ct->findCT(key, res)) {
        //
        // 'main' compute table hit
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
            chainToLevel(C, Clevel, L, in);
        }
        return;
        //
        // done 'main' compute table hit
        //
    }
#ifdef USE_PRIMED_CACHE
    if (useCTpr && ct_primed->findCT(key, res)) {
        //
        // 'primed' compute table hit
        //
        C = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "CT' hit ";
        key.show(out);
        out << " -> ";
        res.show(out);
        out << "\n";
#endif
        if (L == resF->getNodeLevel(C)) {
            // Make sure we don't point to a singleton
            // from the same index.
            C = resF->redirectSingleton(in, C);
        }
        return;
        //
        // done 'primed' compute table hit
        //
    }
#endif // USE_PRIMED_CACHE


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
            MEDDLY_DCASSERT(Clevel == L);
            // ^ if we skip too much, in is irrelevant
            Au->initIdentity(Clevel, in, A);
            MEDDLY_DCASSERT(Au->wasIdentity());
        } else {
            Au->initRedundant(Clevel, A);
            MEDDLY_DCASSERT(!Au->wasIdentity());
        }
    } else {
        Au->initFromNode(A);
        MEDDLY_DCASSERT(!Au->wasIdentity());
    }

    unpacked_node* Bu = unpacked_node::New(arg2F, FULL_ONLY);
    if (Blevel != Clevel) {
        if (arg2F->isIdentityReduced() && Clevel<0) {
            MEDDLY_DCASSERT(Clevel == L);
            Bu->initIdentity(Clevel, in, B);
            MEDDLY_DCASSERT(Bu->wasIdentity());
        } else {
            Bu->initRedundant(Clevel, B);
            MEDDLY_DCASSERT(!Bu->wasIdentity());
        }
    } else {
        Bu->initFromNode(B);
        MEDDLY_DCASSERT(!Bu->wasIdentity());
    }

    unpacked_node* Cu = unpacked_node::newWritable(resF, Clevel, FULL_ONLY);
    MEDDLY_DCASSERT(Cu->getSize() == Au->getSize());
    MEDDLY_DCASSERT(Cu->getSize() == Bu->getSize());

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
    out << "union_mt::_compute(" << L << ", " << in << ", " << A << ", "
        << B << ") done\n";
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
    if (useCT) {
        if (Au->wasIdentity() || Bu->wasIdentity()) {
            ct->noaddCT(key);
        } else {
            res[0].setN(C);
            ct->addCT(key, res);
        }
#ifdef USE_PRIMED_CACHE
    } else if (useCTpr) {
        MEDDLY_DCASSERT(!Au->wasIdentity());
        MEDDLY_DCASSERT(!Bu->wasIdentity());
        res[0].setN(C);
        ct_primed->addCT(key, res);
#endif // USE_PRIMED_CACHE
    } else {
        // recycle the key
        ct->noaddCT(key);
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
        chainToLevel(C, Clevel, L, in);
    }

}


// ******************************************************************
// *                                                                *
// *                      UNION_factory  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::UNION_factory : public binary_factory {
    public:
        virtual void setup();
        virtual binary_operation* build(forest* a, forest* b, forest* c);
};

// ******************************************************************

void MEDDLY::UNION_factory::setup()
{
    _setup(__FILE__, "UNION", "Set union, meaning logical OR of the inputs. Operands and the result may be within the same forest or across forests, but all forests must be over the same domain. Forests must be multi-terminal.");
}

MEDDLY::binary_operation*
MEDDLY::UNION_factory::build(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }
    binary_operation* bop =  cache_find(a, b, c);
    if (bop) {
        return bop;
    }

    return cache_add(new union_mt(a, b, c));
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_factory& MEDDLY::UNION()
{
    static UNION_factory F;
    return F;
}


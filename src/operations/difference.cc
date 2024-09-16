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
#include "difference.h"

#include "../ops_builtin.h" // for COPY
#include "../oper_binary.h"
#include "../oper_unary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class diffr_mt;

    class diffr_mdd;
    class diffr_mxd;

    binary_list DIFFR_cache;
};

// #define TRACE

#define NEW_DIFF

#ifdef TRACE
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                         diffr_mt class                         *
// *                                                                *
// ******************************************************************

//
// TBD: adjust to do the following.
//
//      quasi - quasi   :   no level skipping will occur
//      quasi - fully   :   (same)
//      quasi - ident   :   (same)
//      fully - quasi   :   (same)
//      ident - quasi   :   (same)
//
//      ident - ident   :   identity pattern
//                          skip to top unprimed level
//                          [p 0 0]   [q 0 0]   [ p-q  0   0  ]
//                          [0 p 0] - [0 q 0] = [  0  p-q  0  ]
//                          [0 0 p]   [0 0 q]   [  0   0  p-q ]
//
//      ident - fully   :   identity pattern
//                          skip to top unprimed level
//                          [p 0 0]   [q q q]   [ p-q  0   0  ]
//                          [0 p 0] - [q q q] = [  0  p-q  0  ]
//                          [0 0 p]   [q q q]   [  0   0  p-q ]
//
//      fully - fully   :   fully pattern
//                          skip to top level
//                          [p p p]   [q q q]   [ p-q p-q p-q ]
//                          [p p p] - [q q q] = [ p-q p-q p-q ]
//                          [p p p]   [q q q]   [ p-q p-q p-q ]
//
//      fully - ident   :   Go by levels to build result pattern.
//                          Level information included in CT entries.
//                          [p p p]   [q 0 0]   [ p-q  p   p  ]
//                          [p p p] - [0 q 0] = [  p  p-q  p  ]
//                          [p p p]   [0 0 q]   [  p   p  p-q ]
//
//  skip to top unprimed:
//      use CT for unprimed levels, and have a separate entry
//      for the primed levels?
//      (usable when both nodes are at the same primed level).
//
//  skip to top:
//      works just like 2L level fully/quasi reduced.
//

class MEDDLY::diffr_mt : public binary_operation {
    public:
        diffr_mt(forest* arg1, forest* arg2, forest* res);
        virtual ~diffr_mt();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

    protected:
        void _compute(bool useCT, int L, unsigned in,
                node_handle A, node_handle B, node_handle &C);

    private:
        //
        //  Level we should build, based on reduction types of arguments.
        //
        //  @param  L       topmost level of computation, from the recursion
        //
        //  @param  Alevel  level of the A node (argument 1)
        //  @param  Blevel  level of the B node (argument 2)
        //
        //  @param  force_unprimed
        //                  on output, true if we're forcing computation
        //                  at the unprimed level above the top of A and B.
        //
        //  @return         What level should we build
        //
        inline int topLevelOf(int L, int Alevel, int Blevel,
                bool &force_unprimed) const
        {
            force_unprimed = false;
            if (force_by_levels) {
                //
                // No level skipping allowed
                //
                return L;
            }
            if (! resF->isForRelations()) {
                //
                // MDD; skip as much as we can
                //
                return MDD_levels::topLevel(Alevel, Blevel);
            }

            const int topAB = MXD_levels::topLevel(Alevel, Blevel);

            if ( (topAB == L) || (topAB>0) || !arg1F->isIdentityReduced())
            {
                return topAB;
            }

            //
            //  Still here? Then
            //      topAB is primed
            //      topAB is below level L
            //      arg1F is identity reduced
            //
            //  Now, determine if this operation must proceed at the
            //  level just above topAB, which is the case if the operation
            //  involves identity reduced nodes:
            //
            //      (1) Alevel is below topAB
            //  or
            //      (2) Blevel is below topAB
            //
            //  If that's the case, set force_unprimed to true
            //  and use the topUnprimed.
            //
            //  Otherwise we can proceed as usual?
            //

            MEDDLY_DCASSERT(topAB < 0);
            MEDDLY_DCASSERT(arg1F->isIdentityReduced());

            force_unprimed = (Alevel != topAB) || (Blevel != topAB);

            if (force_unprimed) {
                return -topAB;      // the level above topAB
            }

            return topAB;

        }

    private:
        ct_entry_type* ct;
        unary_operation* copy_arg1res;
        bool force_by_levels;

#ifdef TRACE
        ostream_output out;
#endif
};


// ******************************************************************

MEDDLY::diffr_mt::diffr_mt(forest* arg1, forest* arg2, forest* res)
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

    //
    // How to handle different reductions
    //
    //      quasi - quasi   :   no level skipping will occur
    //      quasi - fully   :   (same)
    //      quasi - ident   :   (same)
    //      fully - quasi   :   (same)
    //      ident - quasi   :   (same)
    //
    //      ident - ident   :   identity pattern (can skip)
    //                          [p 0 0]   [q 0 0]   [ p-q  0   0  ]
    //                          [0 p 0] - [0 q 0] = [  0  p-q  0  ]
    //                          [0 0 p]   [0 0 q]   [  0   0  p-q ]
    //
    //      ident - fully   :   identity pattern (can skip)
    //                          [p 0 0]   [q q q]   [ p-q  0   0  ]
    //                          [0 p 0] - [q q q] = [  0  p-q  0  ]
    //                          [0 0 p]   [q q q]   [  0   0  p-q ]
    //
    //      fully - fully   :   fully pattern (can skip)
    //                          [p p p]   [q q q]   [ p-q p-q p-q ]
    //                          [p p p] - [q q q] = [ p-q p-q p-q ]
    //                          [p p p]   [q q q]   [ p-q p-q p-q ]
    //
    //      fully - ident   :   Go by levels to build result pattern.
    //                          Level information included in CT entries.
    //                          [p p p]   [q 0 0]   [ p-q  p   p  ]
    //                          [p p p] - [0 q 0] = [  p  p-q  p  ]
    //                          [p p p]   [0 0 q]   [  p   p  p-q ]
    //
    force_by_levels = arg1->isFullyReduced() && arg2->isIdentityReduced();

    ct = new ct_entry_type("difference");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    if (force_by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::diffr_mt::~diffr_mt()
{
    ct->markForDestroy();
}

void MEDDLY::diffr_mt::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
#ifdef TRACE
    out.indentation(0);
#endif
    _compute(true, L, in, ap, bp, cp);
    cv.set();
}

void MEDDLY::diffr_mt::_compute(bool useCT, int L, unsigned in,
            node_handle A, node_handle B, node_handle &C)
{
    //
    // Check terminal cases
    //
    if (A==0) {
        // 0 - B = 0
        C = 0;
        return;
    }

    if (B < 0) {
        if (arg2F->isFullyReduced() || 0==L) {
            // A -1 = 0
            C = 0;
            return;
        }
        //
        // Must be identity reduced, so B=I
        //
        MEDDLY_DCASSERT(arg2F->isIdentityReduced());
        if (A < 0) {
            if (arg1F->isIdentityReduced()) {
                // I - I
                C = 0;
                return;
            }

            // 1 - I, need to compute it
        }
    }

    if (B==0) {
        // A - 0 = A
        edge_value dummy;
        dummy.set();
        MEDDLY_DCASSERT(copy_arg1res);
        copy_arg1res->compute(L, in,  dummy, A, dummy, C);
        MEDDLY_DCASSERT(dummy.isVoid());
        return;
    }

    if (A == B) {
        if (arg1F == arg2F && !force_by_levels) {
            C = 0;
            return;
        }
    }

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    bool force_unprimed;
    const int Clevel = topLevelOf(L, Alevel, Blevel, force_unprimed);
    const int Cnextlevel = resF->isForRelations()
        ? MXD_levels::downLevel(Clevel)
        : MDD_levels::downLevel(Clevel)
    ;

#ifdef TRACE
    out << "diffr_mt::_compute(" << L << ", " << in << ", " << A << ", "
        << B << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << "; force_unprimed: "
        << force_unprimed << "\n";
#endif

    //
    // Check compute table, unless we can't
    // (we're recursing on the primed part of
    // a recursion that forces us to unprimed levels).
    //
    ct_vector key( force_by_levels ? 3 : 2);
    ct_vector res(1);
    if (force_by_levels) {
        key[0].setI(L);
        key[1].setN(A);
        key[2].setN(B);
    } else {
        key[0].setN(A);
        key[1].setN(B);
    }
    if (useCT && ct->findCT(key, res)) {
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

        if (resF->isQuasiReduced()) {
            if (C && (Clevel != resF->getNodeLevel(C))) {
                std::cout << "\nLevel fail\n";
                std::cout << "         L: " << L << "\n";
                std::cout << "    Alevel: " << Alevel << "\n";
                std::cout << "    Blevel: " << Blevel << "\n";
                std::cout << "    Clevel: " << Clevel << "\n";
                std::cout << "    force_unprimed: " << force_unprimed << "\n";

                std::cout << "    actual level: " << resF->getNodeLevel(C) << "\n";
            }
        }

        //
        // done compute table hit
        //
    } else {
        //
        // compute table 'miss'; do computation
        //

        //
        // Set up unpacked nodes
        //

        unpacked_node* Au;
        if (Alevel != Clevel) {
            if (arg1F->isIdentityReduced() && Clevel<0) {
                MEDDLY_DCASSERT(Clevel == L);
                // ^ if we skip too much, in is irrelevant
                Au = unpacked_node::newIdentity(arg1F, Clevel, in,
                        A, SPARSE_ONLY);
            } else {
                Au = unpacked_node::newRedundant(arg1F, Clevel,
                        A, SPARSE_ONLY);
            }
        } else {
            Au = arg1F->newUnpacked(A, SPARSE_ONLY);
        }

        unpacked_node* Bu;
        if (Blevel != Clevel) {
            if (arg2F->isIdentityReduced() && Clevel<0) {
                MEDDLY_DCASSERT(Clevel == L);
                Bu = unpacked_node::newIdentity(arg2F, Clevel, in,
                        B, FULL_ONLY);
            } else {
                Bu = unpacked_node::newRedundant(arg2F, Clevel,
                        B, FULL_ONLY);
            }
        } else {
            Bu = arg2F->newUnpacked(B, FULL_ONLY);
        }

        unpacked_node* Cu =
            unpacked_node::newSparse(resF, Clevel, Au->getSize());

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

        unsigned zc = 0;
        for (unsigned z=0; z<Au->getSize(); z++) {
            const unsigned i = Au->index(z);
            node_handle cd;
            _compute(!force_unprimed, Cnextlevel, int(i),
                        Au->down(z), Bu->down(i), cd);
            if (cd) {
                Cu->setSparse(zc++, i, cd);
            }
        }
        Cu->resize(zc);

#ifdef TRACE
        out.indent_less();
        out.put('\n');
        out << "diffr_mt::_compute(" << L << ", " << in << ", " << A << ", "
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
        if (!useCT || Au->wasIdentity() || Bu->wasIdentity()) {
            //
            // DON'T save, because result depends on in,
            // or we were told not to use the CT.
            //
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
        // done compute table 'miss'
        //
    }

    //
    // Adjust result for singletons, added identities/redundants.
    // Do this for both CT hits and misses.
    //
    if (L == resF->getNodeLevel(C)) {
        //
        // We don't need to add identities or redundants;
        // just check if we need to avoid a singleton.
        //
        C = resF->redirectSingleton(in, C);
        return;
    }

    //
    // Add nodes but only from Clevel;
    // if the actual level of C is below Clevel it means
    // nodes were eliminated in the result forest.
    //
    if (arg1F->isIdentityReduced()) {
#ifdef TRACE
        out << "I chain to " << C << ", levels " << Clevel << " to " << L << "\n";
#endif
        C = resF->makeIdentitiesTo(C, Clevel, L, in);
    } else {
#ifdef TRACE
        out << "X chain to " << C << ", levels " << Clevel << " to " << L << "\n";
#endif
        C = resF->makeRedundantsTo(C, Clevel, L);
    }
}

//
// TBD: OLD STUFF BELOW
//

// ******************************************************************
// *                                                                *
// *                        diffr_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::diffr_mdd : public binary_operation {
    public:
        diffr_mdd(forest* arg1, forest* arg2, forest* res);
        virtual ~diffr_mdd();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

    protected:
        node_handle _compute(node_handle A, node_handle B, int L);

    private:
        ct_entry_type* ct;
};

// ******************************************************************

MEDDLY::diffr_mdd::diffr_mdd(forest* arg1, forest* arg2, forest* res)
  : binary_operation(arg1, arg2, res)
{
    ct = new ct_entry_type("difference");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    ct->setFixed(arg1, arg2);
    ct->setResult(res);
    ct->doneBuilding();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

MEDDLY::diffr_mdd::~diffr_mdd()
{
    ct->markForDestroy();
}

void MEDDLY::diffr_mdd::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
    cv.set();
    cp = _compute(ap, bp, L);
}

MEDDLY::node_handle
MEDDLY::diffr_mdd::_compute(node_handle A, node_handle B, int L)
{
    //
    // Check terminal cases
    //
    if (A==0) {
        // 0 - B = 0
        return 0;
    }
    if (B < 0) {
        // A - 1 = 0
        return 0;
    }

    if (B==0) {
        if (A < 0) {
            // 1 - 0 = 1
            terminal tt(true);
            return resF->makeRedundantsTo(tt.getHandle(), 0, L);
        }
        // A - 0 = A
        // return A if we can (same forest as result)
        if (arg1F == resF) {
            return resF->linkNode(A);
            // Shouldn't need to make redundant chain b/c same forest
        }
    }

    if (A == B) {
        if (arg1F == arg2F) {
            return 0;
        }
    }

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = MAX(Alevel, Blevel);

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
    // Initialize unpacked nodes.
    // Sparse for A, full for B
    //
    unpacked_node* Au = (Alevel < Clevel)
        ?   unpacked_node::newRedundant(arg1F, Clevel, A, SPARSE_ONLY)
        :   arg1F->newUnpacked(A, SPARSE_ONLY);

    unpacked_node* Bu = (Blevel <Clevel)
        ?   unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY)
        :   arg2F->newUnpacked(B, FULL_ONLY);

    unpacked_node* Cu = unpacked_node::newSparse(resF, Clevel, Au->getSize());

    //
    // Build result node
    // Scan through (sparse) entries of A,
    // and subtract the corresponding entries in B.
    //
    unsigned zc = 0;
    for (unsigned z=0; z<Au->getSize(); z++) {
        const unsigned i = Au->index(z);
        const node_handle cd = _compute(Au->down(z), Bu->down(i), Clevel-1);
        if (cd) {
            Cu->setSparse(zc++, i, cd);
        }
    }
    Cu->resize(zc);

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
// *                        diffr_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::diffr_mxd : public binary_operation {
    public:
        diffr_mxd(forest* arg1, forest* arg2, forest* res);
        virtual ~diffr_mxd();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

    protected:
        node_handle _compute(node_handle A, node_handle B, int L);
        node_handle _compute_primed(int in, node_handle A, node_handle B,
                int L);

    protected:
        inline node_handle makeChainTo(node_handle p, int K, int L)
        {
            if (arg1F->isIdentityReduced()) {
                return resF->makeIdentitiesTo(p, K, L, -1);
            }
            if (arg1F->isFullyReduced() && !force_by_levels) {
                return resF->makeRedundantsTo(p, K, L);
            }
            MEDDLY_DCASSERT(K == L);
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


    private:
        ct_entry_type* ct;

        bool force_by_levels;

#ifdef TRACE
        ostream_output out;
#endif
};

// ******************************************************************

MEDDLY::diffr_mxd::diffr_mxd(forest* arg1, forest* arg2, forest* res)
  : binary_operation(arg1, arg2, res)
#ifdef TRACE
    , out(std::cout)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    //
    // How to handle different reductions
    //
    //      quasi - quasi   :   no level skipping will occur
    //      quasi - fully   :   (same)
    //      quasi - ident   :   (same)
    //      fully - quasi   :   (same)
    //      ident - quasi   :   (same)
    //
    //      ident - ident   :   identity pattern (can skip)
    //                          [p 0 0]   [q 0 0]   [ p-q  0   0  ]
    //                          [0 p 0] - [0 q 0] = [  0  p-q  0  ]
    //                          [0 0 p]   [0 0 q]   [  0   0  p-q ]
    //
    //      ident - fully   :   identity pattern (can skip)
    //                          [p 0 0]   [q q q]   [ p-q  0   0  ]
    //                          [0 p 0] - [q q q] = [  0  p-q  0  ]
    //                          [0 0 p]   [q q q]   [  0   0  p-q ]
    //
    //      fully - fully   :   fully pattern (can skip)
    //                          [p p p]   [q q q]   [ p-q p-q p-q ]
    //                          [p p p] - [q q q] = [ p-q p-q p-q ]
    //                          [p p p]   [q q q]   [ p-q p-q p-q ]
    //
    //      fully - ident   :   Go by levels to build result pattern.
    //                          Level information included in CT entries.
    //                          [p p p]   [q 0 0]   [ p-q  p   p  ]
    //                          [p p p] - [0 q 0] = [  p  p-q  p  ]
    //                          [p p p]   [0 0 q]   [  p   p  p-q ]
    //
    force_by_levels = arg1->isFullyReduced() && arg2->isIdentityReduced();

    ct = new ct_entry_type("difference");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    if (force_by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::diffr_mxd::~diffr_mxd()
{
    ct->markForDestroy();
}

void MEDDLY::diffr_mxd::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
    cv.set();
#ifdef TRACE
    out.indentation(0);
#endif
    cp = _compute(ap, bp, L);
}

MEDDLY::node_handle
MEDDLY::diffr_mxd::_compute(node_handle A, node_handle B, int L)
{
    MEDDLY_DCASSERT(L>=0);

    //
    // Check terminal cases
    //
    if (A==0) {
        // 0 - B = 0
        return 0;
    }

    if (B < 0) {
        if (arg2F->isFullyReduced() || 0==L) {
            // A -1 = 0
            return 0;
        }
        //
        // Must be identity reduced, so B=I
        //
        MEDDLY_DCASSERT(arg2F->isIdentityReduced());
        if (A < 0) {
            if (arg1F->isIdentityReduced()) {
                // I - I
                return 0;
            }

            // 1 - I, need to compute it
        }
    }

    if (B==0) {
        // A - 0 = A
        // return A if we can (same forest as result)
        if (arg1F == resF) {
            return resF->linkNode(A);
            // Shouldn't need to make any chains b/c same forest
        }
        if (A < 0) {
            // Treat the arg1F != resF case when A=1
            // 1 - 0 = 1
            terminal tt(true);
            if (0==L) return tt.getHandle();
            if (arg1F->isIdentityReduced()) {
                return resF->makeIdentitiesTo(tt.getHandle(), 0, L, -1);
            }
            if (arg1F->isFullyReduced()) {
                return resF->makeRedundantsTo(tt.getHandle(), 0, L);
            }
            MEDDLY_DCASSERT(false);
        }
    }

    if (A == B) {
        if (arg1F == arg2F && !force_by_levels) {
            return 0;
        }
    }

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = force_by_levels
        ? L
        : MAX(ABS(Alevel), ABS(Blevel));

#ifdef TRACE
    out << "diffr_mxd::_compute(" << A << ", " << B << ", " << L << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << " before chain\n";
#endif

    //
    // Check compute table
    //
    ct_vector key( force_by_levels ? 3 : 2);
    ct_vector res(1);
    if (force_by_levels) {
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
        out << "CT hit " << res[0].getN() << "\n";
        out << "after chain " << C << "\n";
#endif
        return C;
    }

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    // Sparse for A, full for B
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   unpacked_node::newRedundant(arg1F, Clevel, A, SPARSE_ONLY)
        :   arg1F->newUnpacked(A, SPARSE_ONLY);

    unpacked_node* Bu = (Blevel != Clevel)
        ?   unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY)
        :   arg2F->newUnpacked(B, FULL_ONLY);

    unpacked_node* Cu = unpacked_node::newSparse(resF, Clevel, Au->getSize());


    //
    // Build result node
    // Scan through (sparse) entries of A,
    // and subtract the corresponding entries in B.
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    unsigned zc = 0;
    for (unsigned z=0; z<Au->getSize(); z++) {
        const unsigned i = Au->index(z);
        const node_handle cd = _compute_primed(int(i), Au->down(z),
                Bu->down(i), MXD_levels::downLevel(Clevel));
        if (cd) {
            Cu->setSparse(zc++, i, cd);
        }
    }
    Cu->resize(zc);
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "diffr_mxd::_compute(" << A << ", " << B << ", " << L
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
    C = makeChainTo(C, Clevel, L);

#ifdef TRACE
    out << "chain to level " << L << " = " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif
    return C;
}

MEDDLY::node_handle
MEDDLY::diffr_mxd::_compute_primed(int in, node_handle A, node_handle B,
        const int Clevel)
{
    MEDDLY_DCASSERT(Clevel<0);
    MEDDLY_DCASSERT(A != 0);

    //
    // Terminal cases
    //
    if (A == B) {
        if (arg1F == arg2F) {
            return 0;
        }
    }

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
#ifdef TRACE
    out << "diffr_mxd::_compute_primed(" << in << ", " << A << ", " << B << ", " << Clevel << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
#endif

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   patternNode(arg1F, Clevel, in, A, SPARSE_ONLY)
        :   arg1F->newUnpacked(A, SPARSE_ONLY);

    unpacked_node* Bu = (Blevel != Clevel)
        ?   patternNode(arg2F, Clevel, in, B, FULL_ONLY)
        :   arg2F->newUnpacked(B, FULL_ONLY);

    unpacked_node* Cu = unpacked_node::newSparse(resF, Clevel, Au->getSize());

    //
    // Build result node
    // Scan through (sparse) entries of A,
    // and subtract the corresponding entries in B.
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    unsigned zc = 0;
    for (unsigned z=0; z<Au->getSize(); z++) {
        const unsigned i = Au->index(z);
        const node_handle cd = _compute(Au->down(z), Bu->down(i),
                MXD_levels::downLevel(Clevel));
        if (cd) {
            Cu->setSparse(zc++, i, cd);
        }
    }
    Cu->resize(zc);
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "diffr_mxd::_compute_primed(" << in << ", " << A << ", " << B << ", " << Clevel << ") done\n";
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
    out << "diffr_mxd::_compute_primed(" << in << ", " << A << ", "
              << B << ", " << Clevel << ") = " << C << "\n\t";
    resF->showNode(out, C, SHOW_DETAILS);
    out.put('\n');
#endif
    return C;

}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation*
MEDDLY::DIFFERENCE(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }
    binary_operation* bop =  DIFFR_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

#ifdef NEW_DIFF
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        return DIFFR_cache.add(new diffr_mt(a, b, c));
    }
#else
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations()) {
            return DIFFR_cache.add(new diffr_mxd(a, b, c));
        } else {
            return DIFFR_cache.add(new diffr_mdd(a, b, c));
        }
    }
#endif

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::DIFFERENCE_init()
{
    DIFFR_cache.reset("Difference");
}

void MEDDLY::DIFFERENCE_done()
{
    MEDDLY_DCASSERT(DIFFR_cache.isEmpty());
}

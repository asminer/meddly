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

class MEDDLY::diffr_mt : public binary_operation {
    public:
        diffr_mt(forest* arg1, forest* arg2, forest* res);
        virtual ~diffr_mt();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

    protected:
        void _compute(int L, unsigned in,
                node_handle A, node_handle B, node_handle &C);

    private:
        inline int topLevelOf(int L, int Alevel, int Blevel) const
        {
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

            if (force_by_unprimed) {
                const int t = MXD_levels::topUnprimed(Alevel, Blevel);
                return (t == MXD_levels::upLevel(L)) ? L : t;
            } else {
                return MXD_levels::topLevel(Alevel, Blevel);
            }
        }

        inline void chainToLevel(node_handle &C, int Clevel, int L, unsigned in)
        {
            //
            // Add nodes from Clevel to L
            //
            if (arg1F->isIdentityReduced()) {
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
        ct_entry_type* ct_primed;
        unary_operation* copy_arg1res;
        bool force_by_levels;
        bool force_by_unprimed;

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
    //      quasi - fully   :
    //      quasi - ident   :
    //      fully - quasi   :
    //      ident - quasi   :
    //
    //      ident - ident   :   identity pattern: skip by unprimed
    //
    //                          [p 0 0]   [q 0 0]   [ p-q  0   0  ]
    //                          [0 p 0] - [0 q 0] = [  0  p-q  0  ]
    //                          [0 0 p]   [0 0 q]   [  0   0  p-q ]
    //
    //      ident - fully   :   identity pattern: skip by unprimed
    //
    //                          [p 0 0]   [q q q]   [ p-q  0   0  ]
    //                          [0 p 0] - [q q q] = [  0  p-q  0  ]
    //                          [0 0 p]   [q q q]   [  0   0  p-q ]
    //
    //      fully - fully   :   fully pattern: skip by top level
    //
    //                          [p p p]   [q q q]   [ p-q p-q p-q ]
    //                          [p p p] - [q q q] = [ p-q p-q p-q ]
    //                          [p p p]   [q q q]   [ p-q p-q p-q ]
    //
    //      fully - ident   :   Go by levels to build result pattern.
    //                          Level information included in CT entries.
    //
    //                          [p p p]   [q 0 0]   [ p-q  p   p  ]
    //                          [p p p] - [0 q 0] = [  p  p-q  p  ]
    //                          [p p p]   [0 0 q]   [  p   p  p-q ]
    //
    force_by_levels = arg1->isFullyReduced() && arg2->isIdentityReduced();
    force_by_unprimed = arg1->isIdentityReduced() && !arg2->isQuasiReduced();

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

    //
    // If we're force by unprimed levels,
    // keep a second CT for the primed level computations
    // we are able to save. This is needed to differentiate
    // the unprimed level result and the primed level result.
    ct_primed = nullptr;
    if (force_by_unprimed) {
        ct_primed = new ct_entry_type("difference_pr");
        ct_primed->setFixed(arg1, arg2);
        ct_primed->setResult(res);
        ct_primed->doneBuilding();
    }
}

MEDDLY::diffr_mt::~diffr_mt()
{
    ct->markForDestroy();
    if (ct_primed) ct_primed->markForDestroy();
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
    _compute(L, in, ap, bp, cp);
    cv.set();
}

void MEDDLY::diffr_mt::_compute(int L, unsigned in,
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
    const int Clevel = topLevelOf(L, Alevel, Blevel);
    const int Cnextlevel = resF->isForRelations()
        ? MXD_levels::downLevel(Clevel)
        : MDD_levels::downLevel(Clevel)
    ;
    const bool useCT = !force_by_unprimed || Clevel>0;
    const bool useCTpr = !useCT && (L<0) &&
        (Alevel == L) && (Blevel == L) && (Clevel == L);
    MEDDLY_DCASSERT(!useCTpr || ct_primed);


#ifdef TRACE
    out << "diffr_mt::_compute(" << L << ", " << in << ", " << A << ", "
        << B << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << "\n";
#endif

    //
    // Check the compute table, unless we can't
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
            Au = unpacked_node::newIdentity(arg1F, Clevel, in, A, SPARSE_ONLY);
            MEDDLY_DCASSERT(Au->wasIdentity());
        } else {
            Au = unpacked_node::newRedundant(arg1F, Clevel, A, SPARSE_ONLY);
            MEDDLY_DCASSERT(!Au->wasIdentity());
        }
    } else {
        Au = arg1F->newUnpacked(A, SPARSE_ONLY);
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
        _compute(Cnextlevel, int(i), Au->down(z), Bu->down(i), cd);
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
    if (useCT) {
        if (Au->wasIdentity() || Bu->wasIdentity()) {
            ct->noaddCT(key);
        } else {
            res[0].setN(C);
            ct->addCT(key, res);
        }
    } else if (useCTpr) {
        MEDDLY_DCASSERT(!Au->wasIdentity());
        MEDDLY_DCASSERT(!Bu->wasIdentity());
        res[0].setN(C);
        ct_primed->addCT(key, res);
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

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
#include "apply_base.h" // remove this when we can

#include "../oper_binary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class diffr_mdd;
    class diffr_mxd;

    binary_list DIFFR_cache;
};

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                        diffr_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::diffr_mdd : public binary_operation {
    public:
        diffr_mdd(forest* arg1, forest* arg2, forest* res);
        virtual ~diffr_mdd();

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

void MEDDLY::diffr_mdd::compute(const edge_value &av, node_handle ap,
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

        virtual void compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        node_handle _compute(node_handle A, node_handle B, int L);
        node_handle _compute_primed(int in, node_handle A, node_handle B,
                int L);

    protected:
        inline node_handle makeChainTo(node_handle p, int K, int L)
        {
            if (arg1F->isIdentityReduced()) {
                return resF->makeIdentitiesTo(p, K, L);
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

        ostream_output out;
};

// ******************************************************************

MEDDLY::diffr_mxd::diffr_mxd(forest* arg1, forest* arg2, forest* res)
  : binary_operation(arg1, arg2, res), out(std::cout)
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

void MEDDLY::diffr_mxd::compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
    cv.set();
    out.indentation(0);
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
                return resF->makeIdentitiesTo(tt.getHandle(), 0, L);
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
                Bu->down(i), forest::downLevel(Clevel));
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
                forest::downLevel(Clevel));
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

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations()) {
            return DIFFR_cache.add(new diffr_mxd(a, b, c));
        } else {
            return DIFFR_cache.add(new diffr_mdd(a, b, c));
        }
    }

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

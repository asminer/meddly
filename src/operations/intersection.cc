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
#include "intersection.h"

#include "apply_base.h" // remove this when we can

#include "../ops_builtin.h" // for COPY
#include "../oper_binary.h"
#include "../oper_unary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class inter_mt;

    binary_list INTER_cache;
};

// #define TRACE
#define USE_PRIMED_CACHE

#ifdef TRACE
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                         inter_mt class                         *
// *                                                                *
// ******************************************************************

class MEDDLY::inter_mt : public binary_operation {
    public:
        inter_mt(forest* arg1, forest* arg2, forest* res);
        virtual ~inter_mt();

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
            if (! resF->isForRelations()) {
                //
                // MDD; skip as much as we can
                //
                return MDD_levels::topLevel(Alevel, Blevel);
            }

            if (identity_pattern) {
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
            if (identity_pattern) {
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
        bool identity_pattern;

#ifdef TRACE
        ostream_output out;
        unsigned top_count;
#endif
};


// ******************************************************************

MEDDLY::inter_mt::inter_mt(forest* arg1, forest* arg2, forest* res)
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
    copy_arg1res = COPY(arg1, res);
    copy_arg2res = COPY(arg2, res);

    //
    // How to handle different reductions
    //
    //      quasi ^ quasi   :   no level skipping will occur
    //      quasi ^ fully   :
    //      quasi ^ ident   :
    //      fully ^ quasi   :
    //      ident ^ quasi   :
    //
    //      ident ^ ident   :   identity pattern: skip by unprimed
    //
    //                          [p 0 0]   [q 0 0]   [ p^q  0   0  ]
    //                          [0 p 0] ^ [0 q 0] = [  0  p^q  0  ]
    //                          [0 0 p]   [0 0 q]   [  0   0  p^q ]
    //
    //      ident ^ fully   :   identity pattern: skip by unprimed
    //
    //                          [p 0 0]   [q q q]   [ p^q  0   0  ]
    //                          [0 p 0] ^ [q q q] = [  0  p^q  0  ]
    //                          [0 0 p]   [q q q]   [  0   0  p^q ]
    //
    //      fully ^ ident   :   identity pattern: skip by unprimed
    //
    //                          [p p p]   [q 0 0]   [ p^q  0   0  ]
    //                          [p p p] ^ [0 q 0] = [  0  p^q  0  ]
    //                          [p p p]   [0 0 q]   [  0   0  p^q ]
    //
    //      fully ^ fully   :   fully pattern: skip by top level
    //
    //                          [p p p]   [q q q]   [ p^q p^q p^q ]
    //                          [p p p] ^ [q q q] = [ p^q p^q p^q ]
    //                          [p p p]   [q q q]   [ p^q p^q p^q ]
    //
    //
    identity_pattern = (arg1->isIdentityReduced() && !arg2->isQuasiReduced())
                        ||
                       (arg2->isIdentityReduced() && !arg1->isQuasiReduced());

    ct = new ct_entry_type("intersection");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    ct->setFixed(arg1, arg2);
    ct->setResult(res);
    ct->doneBuilding();

#ifdef USE_PRIMED_CACHE
    //
    // If we're skipping by unprimed levels,
    // keep a second CT for the primed level computations
    // we are able to save. This is needed to differentiate
    // the unprimed level result and the primed level result.
    ct_primed = nullptr;
    if (identity_pattern) {
        ct_primed = new ct_entry_type("intersection_pr");
        ct_primed->setFixed(arg1, arg2);
        ct_primed->setResult(res);
        ct_primed->doneBuilding();
    }
#endif // USE_PRIMED_CACHE
}

MEDDLY::inter_mt::~inter_mt()
{
    ct->markForDestroy();
#ifdef USE_PRIMED_CACHE
    if (ct_primed) ct_primed->markForDestroy();
#endif
}

void MEDDLY::inter_mt::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << "Intersection #" << top_count << " begin\n";
#endif
    _compute(L, in, ap, bp, cp);
#ifdef TRACE
    out << "Intersection #" << top_count << " end\n";
#endif
    cv.set();
}

void MEDDLY::inter_mt::_compute(int L, unsigned in,
            node_handle A, node_handle B, node_handle &C)
{
    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************
    if (A==0 || B==0) {
        C = 0;
        return;
    }

    if (arg1F->isTerminalNode(A)) {
        if (L==0 || arg1F->isFullyReduced()) {
            // TRUE and B = B
            edge_value dummy;
            dummy.set();
            MEDDLY_DCASSERT(copy_arg2res);
            copy_arg2res->compute(L, in, dummy, B, dummy, C);
            MEDDLY_DCASSERT(dummy.isVoid());
            return;
        }
    }

    if (arg2F->isTerminalNode(B)) {
        if (L==0 || arg2F->isFullyReduced()) {
            // A and TRUE = A
            edge_value dummy;
            dummy.set();
            MEDDLY_DCASSERT(copy_arg1res);
            copy_arg1res->compute(L, in, dummy, A, dummy, C);
            MEDDLY_DCASSERT(dummy.isVoid());
            return;
        }
    }

    if ((A == B) && (arg1F==arg2F)) {
        // A and A = A
        edge_value dummy;
        dummy.set();
        MEDDLY_DCASSERT(copy_arg1res);
        copy_arg1res->compute(L, in, dummy, A, dummy, C);
        MEDDLY_DCASSERT(dummy.isVoid());
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
    const bool useCT = !identity_pattern || Clevel>0;
#ifdef USE_PRIMED_CACHE
    const bool useCTpr = !useCT && (L<0) &&
        (Alevel == L) && (Blevel == L) && (Clevel == L);
    MEDDLY_DCASSERT(!useCTpr || ct_primed);
#endif

#ifdef TRACE
    out << "inter_mt::_compute(" << L << ", " << in << ", " << A << ", "
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
    key[0].setN(A);
    key[1].setN(B);


#ifdef TRACE
    if (useCT) {
        out << "Check CT ";
        key.show(out);
        out << "\n";
    }
#endif
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
#ifdef TRACE
    if (useCTpr) {
        out << "Check CT' ";
        key.show(out);
        out << "\n";
    }
#endif
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
            Bu = unpacked_node::newIdentity(arg2F, Clevel, in, B, SPARSE_ONLY);
            MEDDLY_DCASSERT(Bu->wasIdentity());
        } else {
            Bu = unpacked_node::newRedundant(arg2F, Clevel, B, SPARSE_ONLY);
            MEDDLY_DCASSERT(!Bu->wasIdentity());
        }
    } else {
        Bu = arg2F->newUnpacked(B, SPARSE_ONLY);
        MEDDLY_DCASSERT(!Bu->wasIdentity());
    }

    unpacked_node* Cu = unpacked_node::newSparse(resF, Clevel,
            MIN(Au->getSize(), Bu->getSize()));

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

    unsigned za=0, zb=0, zc=0;
    while ( (za < Au->getSize()) && (zb < Bu->getSize()) ) {
        if (Au->index(za) < Bu->index(zb)) {
            ++za;
            continue;
        }
        if (Au->index(za) > Bu->index(zb)) {
            ++zb;
            continue;
        }
        const unsigned i = Au->index(za);
        MEDDLY_DCASSERT(i == Bu->index(zb));

        node_handle cd;
        _compute(Cnextlevel, i, Au->down(za), Bu->down(zb), cd);
        if (cd) {
            Cu->setSparse(zc++, i, cd);
        }
        za++;
        zb++;
    }
    Cu->resize(zc);
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "inter_mt::_compute(" << L << ", " << in << ", " << A << ", "
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
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation*
MEDDLY::INTERSECTION(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }

    binary_operation* bop =  INTER_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    return INTER_cache.add(new inter_mt(a, b, c));
}

void MEDDLY::INTERSECTION_init()
{
    INTER_cache.reset("Intersection");
}

void MEDDLY::INTERSECTION_done()
{
    MEDDLY_DCASSERT(INTER_cache.isEmpty());
}


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
#include "complement.h"

#include "../ct_vector.h"
#include "../compute_table.h"
#include "../oper_unary.h"

namespace MEDDLY {
    class compl_mt;

    class compl_mdd;
    class compl_mxd;

    unary_list COMPL_cache;
};

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                         compl_mt class                         *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_mt : public unary_operation {
    public:
        compl_mt(forest* arg, forest* res);
        virtual ~compl_mt();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp);

    protected:
        /// Implement compute(), recursively
        void _compute(int L, unsigned in, node_handle A, node_handle &C);

        /*
            Build a complemented identity pattern:
                [ p 1 1 ... 1 ]
                [ 1 p 1 ... 1 ]
                [ 1 1 p ... 1 ]
                [ :         : ]
                [ 1 1 1 ... p ]


                @param  p   Diagonal pointer
                @param  K   Start the pattern above this level
                @param  L   Stop the pattern at this level
                @param  in  Incoming edge index; needed only if L is primed
        */
        node_handle _identity_complement(node_handle p, int K, int L,
                unsigned in);

        inline node_handle identity_complement(node_handle p, int K, int L,
                unsigned in)
        {
            if (0==L) return p;
            if (K==L) return p;
            MEDDLY_DCASSERT(MXD_levels::topLevel(L, K) == L);
            return _identity_complement(p, K, L, in);
        }

    private:
        ct_entry_type* ct;

#ifdef TRACE
        ostream_output out;
        unsigned top_count;
#endif
};

// ******************************************************************

MEDDLY::compl_mt::compl_mt(forest* arg, forest* res)
    : unary_operation(arg, res)
#ifdef TRACE
      , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("mdd_complement");
    ct->setFixed(arg);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::compl_mt::~compl_mt()
{
    ct->markForDestroy();
}

void MEDDLY::compl_mt::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    ++top_count;
    out.indentation(0);
    out << "Starting top-level compl_mt::compute #" << top_count << "\n";
#endif
    MEDDLY_DCASSERT(av.isVoid());
    cv.set();
    _compute(L, in, ap, cp);
#ifdef TRACE
    out.indentation(0);
    out << "Finishing top-level compl_mt::compute #" << top_count << "\n";
#endif
}

void MEDDLY::compl_mt::_compute(int L, unsigned in,
        node_handle A, node_handle &cp)
{
    //
    // Terminal cases.
    //
    if (argF->isTerminalNode(A)) {
        bool ta;
        argF->getValueFromHandle(A, ta);
        cp = resF->handleForValue(!ta);
        //
        // Add nodes up to level L.
        //
        if (argF->isIdentityReduced()) {
            cp = identity_complement(cp, 0, L, in);
        } else {
            cp = resF->makeRedundantsTo(cp, 0, L);
        }
        return;
    }

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);
#ifdef TRACE
    out << "compl_mt::_compute(" << A << ")\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(1);
    ct_vector res(1);
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        //
        // compute table 'hit'
        //
        cp = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "  CT hit " << cp << "\n";
        out << "  at level " << resF->getNodeLevel(cp) << "\n";
#endif
        //
        // done: compute table 'hit'
        //
    } else {
        //
        // compute table 'miss'; do computation
        //

        //
        // Initialize unpacked nodes
        //
        unpacked_node* Au = argF->newUnpacked(A, FULL_ONLY);
        unpacked_node* Cu = unpacked_node::newFull(resF, Alevel, Au->getSize());
#ifdef TRACE
        out << "A: ";
        Au->show(out, true);
        out.indent_more();
        out.put('\n');
#endif
        const int Cnextlevel = resF->isForRelations()
            ? MXD_levels::downLevel(Alevel)
            : MDD_levels::downLevel(Alevel);

        //
        // Recurse over child edges
        //
        for (unsigned i=0; i<Cu->getSize(); i++) {
            node_handle d;
            _compute(Cnextlevel, i, Au->down(i), d);
            Cu->setFull(i, d);
        }

#ifdef TRACE
        out.indent_less();
        out.put('\n');
        out << "compl_mt::_compute(" << A << ") done\n";
        out << "A: ";
        Au->show(out, true);
        out << "\nC: ";
        Cu->show(out, true);
        out << "\n";
#endif

        //
        // Reduce
        //
        edge_value cv;
        resF->createReducedNode(Cu, cv, cp);
        MEDDLY_DCASSERT(cv.isVoid());
#ifdef TRACE
        out << "reduced to " << cp << ": ";
        resF->showNode(out, cp, SHOW_DETAILS);
        out << "\n";
#endif

        //
        // Add to CT
        //
        res[0].setN(cp);
        ct->addCT(key, res);

        //
        // Cleanup
        //
        unpacked_node::Recycle(Au);
        // Cu is recycled when we reduce it :)

        //
        // done: compute table 'miss'
        //
    }

    //
    // Adjust result for singletons, added identities/redundants.
    // Do this for both CT hits and misses.
    //
    const int Clevel = resF->getNodeLevel(cp);
    if (Clevel == L) {
        //
        // We don't need to add identities or redundants;
        // just check if we need to avoid a singleton.
        //
        cp = resF->redirectSingleton(in, cp);
        return;
    }

    //
    // Add nodes but only from Alevel;
    // if Clevel is below Alevel it means nodes were eliminated
    // in the result forest.
    //
    if (argF->isIdentityReduced()) {
        cp = identity_complement(cp, Alevel, L, in);
    } else {
        cp = resF->makeRedundantsTo(cp, Alevel, L);
    }
}

MEDDLY::node_handle MEDDLY::compl_mt::_identity_complement(node_handle p,
    int K, int L, unsigned in)
{
    MEDDLY_DCASSERT(L!=0);
    unpacked_node* Uun;
    unpacked_node* Upr;
    edge_value voidedge;
    voidedge.set();

    if (K<0) {
        //
        // Add a redundant layer as needed;
        // this also handles the only case when we need to check
        // if p is a singleton node :)
        //
        p = resF->makeRedundantsTo(p, K, MXD_levels::upLevel(K));
        K = MXD_levels::upLevel(K);
    }

    MEDDLY_DCASSERT(K>=0);
    const int Lstop = (L<0) ? MXD_levels::downLevel(L) : L;

    terminal ONE(true);
    node_handle chain_to_one = resF->makeRedundantsTo(ONE.getHandle(), 0, K);

    //
    // Proceed in unprimed, primed pairs
    //
    for (K++; K<=Lstop; K++) {
        Uun = unpacked_node::newFull(resF, K, resF->getLevelSize(K));

        for (unsigned i=0; i<Uun->getSize(); i++) {
            Upr = unpacked_node::newFull(resF, -K, Uun->getSize());
            for (unsigned j=0; j<Uun->getSize(); j++) {
                Upr->setFull(j, resF->linkNode(chain_to_one));
            }
            Upr->setFull(i, (i ? resF->linkNode(p) : p));
            resF->unlinkNode(chain_to_one);
            node_handle h;
            resF->createReducedNode(Upr, voidedge, h);
            MEDDLY_DCASSERT(voidedge.isVoid());
            Uun->setFull(i, h);
        }

        resF->createReducedNode(Uun, voidedge, p);
        MEDDLY_DCASSERT(voidedge.isVoid());

        resF->unlinkNode(chain_to_one);
        chain_to_one = resF->makeRedundantsTo(chain_to_one, K-1, K);
    } // for k

    //
    // Add top primed node, if L is negative
    //
    if (L<0) {
        MEDDLY_DCASSERT(-K == L);
        Upr = unpacked_node::newFull(resF, -K, resF->getLevelSize(-K));
        for (unsigned j=0; j<Upr->getSize(); j++) {
            Upr->setFull(j, resF->linkNode(chain_to_one));
        }
        Upr->setFull(in, p);
        resF->unlinkNode(chain_to_one);
        resF->createReducedNode(Upr, voidedge, p);
        MEDDLY_DCASSERT(voidedge.isVoid());
    }

    resF->unlinkNode(chain_to_one);
    return p;
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_operation* MEDDLY::COMPLEMENT(forest* arg, forest* res)
{
    if (!arg || !res) return nullptr;

    unary_operation* uop = COMPL_cache.find(arg, res);
    if (uop) {
        return uop;
    }
    return COMPL_cache.add(new compl_mt(arg, res));
}

void MEDDLY::COMPLEMENT_init()
{
    COMPL_cache.reset("Complement");
}

void MEDDLY::COMPLEMENT_done()
{
    MEDDLY_DCASSERT( COMPL_cache.isEmpty() );
}


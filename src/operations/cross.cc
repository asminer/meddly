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

// #define NEW_CT

#include "../defines.h"
#include "cross.h"

#ifdef NEW_CT
#include "../ct_vector.h"
#else
#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#endif
#include "../compute_table.h"
#include "../oper_binary.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

namespace MEDDLY {
    class cross_bool;

    binary_list CROSS_cache;
};

// ******************************************************************
// *                                                                *
// *                        cross_bool class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::cross_bool : public binary_operation {
    public:
        cross_bool(forest* a1, forest* a2, forest* res);

#ifdef NEW_CT
        virtual ~cross_bool();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);
#else
        virtual void computeDDEdge(const dd_edge& a, const dd_edge& b,
                dd_edge &c, bool userFlag);
#endif

        node_handle compute_pr(unsigned in, int L, node_handle a,
                node_handle b);

        node_handle compute_un(int L, node_handle a, node_handle b);

    private:

#ifdef NEW_CT
        ct_entry_type* ct;
#endif

#ifdef TRACE
        ostream_output out;
        unsigned top_count;
#endif

};

MEDDLY::cross_bool::cross_bool(forest* arg1, forest* arg2, forest* res)
#ifdef NEW_CT
: binary_operation(arg1, arg2, res)
#else
: binary_operation(CROSS_cache, 1, arg1, arg2, res)
#endif
#ifdef TRACE
    , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    checkRelations(__FILE__, __LINE__, SET, SET, RELATION);

#ifdef NEW_CT
    ct = new ct_entry_type(CROSS_cache.getName());
    ct.setFixed('I', arg1, arg2);
    st->setResult(res);
    ct->doneBuilding();
#else
    ct_entry_type* et = new ct_entry_type(CROSS_cache.getName(), "INN:N");
    et->setForestForSlot(1, arg1);
    et->setForestForSlot(2, arg2);
    et->setForestForSlot(4, res);
    registerEntryType(0, et);
    buildCTs();
#endif
}

#ifdef NEW_CT

MEDDLY::cross_bool::~cross_bool()
{
    ct->markForDestroy();
}

void MEDDLY::cross_bool::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << "Cross #" << top_count << " begin\n";
#endif
    if (L<0) {
        cp = compute_pr(in, L, ap, bp);
    } else {
        cp = compute_un(L, ap, bp, cp);
    }
#ifdef TRACE
    out << "Cross #" << top_count << " end\n";
#endif
    cv.set();
}

#else

void
MEDDLY::cross_bool::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag)
{
    int L = arg1F->getMaxLevelIndex();
    node_handle cnode = compute_un(L, a.getNode(), b.getNode());
    c.set(cnode);
}

#endif

MEDDLY::node_handle MEDDLY::cross_bool::compute_un(int k, node_handle a,
        node_handle b)
{
    MEDDLY_DCASSERT(k>=0);

    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************
    if (0==a || 0==b) {
        return 0;
    }

    if (0==k) {
        return a;
    }

#ifdef TRACE
    out << "cross::compute_un(" << k << ", " << a << ", " << b << ")\n";
#endif

    //
    // TBD: if both a and b skip level k, can we just build
    // a redundant node at level k in the mxd?
    //

    // **************************************************************
    //
    // Check the compute table
    //
    // **************************************************************

#ifdef NEW_CT
    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    key[0].setI(k);
    key[1].setN(a);
    key[2].setN(b);
    if (ct->find(key, res)) {
        //
        // Hit
        //
#ifdef TRACE
        out << "CT hit ";
        key.show(out);
        out << " -> ";
        res.show(out);
        out << "\n";
#endif
        return resF->linkNode(res[0].getN());
    }
#else
    // check compute table
    ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeI(k);
    CTsrch->writeN(a);
    CTsrch->writeN(b);
    CT0->find(CTsrch, CTresult[0]);
    if (CTresult[0]) {
        CT0->recycle(CTsrch);
        return resF->linkNode(CTresult[0].readN());
    }
#endif

    // **************************************************************
    //
    // Compute table 'miss'; do computation
    //
    // **************************************************************

    // Initialize unpacked node
    unpacked_node *A = unpacked_node::New(arg1F, SPARSE_ONLY);
    if (arg1F->getNodeLevel(a) < k) {
        A->initRedundant(k, a);
    } else {
        A->initFromNode(a);
    }

    unpacked_node *C =
        unpacked_node::newWritable(resF, k, A->getSize(), SPARSE_ONLY);

#ifdef TRACE
    out << "A: ";
    Au->show(out, true);
    out.indent_more();
    out.put('\n');
#endif

    //
    // Recurse
    //
    for (unsigned z=0; z<A->getSize(); z++) {
        const unsigned i = A->index(z);
        C->setSparse(z, i, compute_pr(i, -k, A->down(z), b));
    }

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "cross::compute_un(" << k << ", " << a << ", " << b << ")\n";
    out << "  A: ";
    A->show(out, true);
    out << "\n  C: ";
    C->show(out, true);
    out << "\n";
#endif

    //
    // Reduce
    //
    edge_value ev;
    node_handle c;
    resF->createReducedNode(C, ev, c);
    MEDDLY_DCASSERT(ev.isVoid());
#ifdef TRACE
    out << "reduced to " << c << ": ";
    resF->showNode(out, c, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT
    //

#ifdef NEW_CT
    res[0].setN(c);
    ct->addCT(key, res);
#else
    CTresult[0].reset();
    CTresult[0].writeN(c);
    CT0->addEntry(CTsrch, CTresult[0]);
#endif

    //
    // Cleanup
    //
    unpacked_node::Recycle(A);
    return c;
}


MEDDLY::node_handle MEDDLY::cross_bool::compute_pr(unsigned in, int k, node_handle a, node_handle b)
{
    MEDDLY_DCASSERT(k<=0);

    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************
    if (0==a || 0==b) {
        return 0;
    }

#ifdef TRACE
    out << "cross::compute_pr(" << in << ", " << k << ", " << a << ", " << b << ")\n";
#endif

    // **************************************************************
    //
    // Check the compute table
    //
    // **************************************************************

#ifdef NEW_CT
    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    key[0].setI(k);
    key[1].setN(a);
    key[2].setN(b);
    if (ct->find(key, res)) {
        //
        // Hit
        //
#ifdef TRACE
        out << "CT hit ";
        key.show(out);
        out << " -> ";
        res.show(out);
        out << "\n";
#endif
        node_handle c = resF->linkNode(res[0].getN());

        if (k == resF->getNodeLevel(c)) {
            // Make sure we don't point to a singleton
            // from the same index.
            c = resF->redirectSingleton(in, c);
        }

        return c;
    }
#else

    // DON'T check compute table

#endif

    // **************************************************************
    //
    // Compute table 'miss'; do computation
    //
    // **************************************************************

    // Initialize unpacked node
    unpacked_node *B = unpacked_node::New(arg2F, SPARSE_ONLY);
    if (arg2F->getNodeLevel(b) < -k) {
        B->initRedundant(-k, b);
    } else {
        B->initFromNode(b);
    }

    unpacked_node *C =
        unpacked_node::newWritable(resF, k, B->getSize(), SPARSE_ONLY);

#ifdef TRACE
    out << "B: ";
    Bu->show(out, true);
    out.indent_more();
    out.put('\n');
#endif

    //
    // Recurse
    //
    for (unsigned z=0; z<B->getSize(); z++) {
        C->setSparse(z, B->index(z), compute_un(-(k+1), a, B->down(z)));
    }

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "cross::compute_pr(" << k << ", " << a << ", " << b << ")\n";
    out << "  B: ";
    B->show(out, true);
    out << "\n  C: ";
    C->show(out, true);
    out << "\n";
#endif

    //
    // Reduce
    //
    edge_value ev;
    node_handle c;
    resF->createReducedNode(C, ev, c);
    MEDDLY_DCASSERT(ev.isVoid());
#ifdef TRACE
    out << "reduced to " << c << ": ";
    resF->showNode(out, c, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT
    //

#ifdef NEW_CT
    res[0].setN(c);
    ct->addCT(key, res);
#else

    // DON'T save in compute table

#endif

    //
    // Make sure we don't point to a singleton from the same index
    //
    if (k == resF->getNodeLevel(c)) {
        c = resF->redirectSingleton(in, c);
    }

    //
    // Cleanup
    //
    unpacked_node::Recycle(B);
    return c;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation*
MEDDLY::CROSS(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }
    binary_operation* bop =  CROSS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    return CROSS_cache.add(new cross_bool(a, b, c));
}

void MEDDLY::CROSS_init()
{
    CROSS_cache.reset("Cross");
}

void MEDDLY::CROSS_done()
{
    MEDDLY_DCASSERT(CROSS_cache.isEmpty());
}


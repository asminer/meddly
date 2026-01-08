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
#include "../ct_vector.h"
#include "../compute_table.h"
#include "../oper_binary.h"
#include "../forest_levels.h"
#include "cross.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

namespace MEDDLY {
    class cross_bool;
    class CROSS_factory;
};

// ******************************************************************
// *                                                                *
// *                        cross_bool class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::cross_bool : public binary_operation {
    public:
        cross_bool(forest* a1, forest* a2, forest* res);
        virtual ~cross_bool();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

        node_handle compute_pr(unsigned in, int L, node_handle a,
                node_handle b);

        node_handle compute_un(int L, node_handle a, node_handle b);

    private:
        ct_entry_type* ct;

#ifdef TRACE
        ostream_output out;
        unsigned top_count;
#endif

};

MEDDLY::cross_bool::cross_bool(forest* arg1, forest* arg2, forest* res)
: binary_operation(arg1, arg2, res)
#ifdef TRACE
    , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    checkRelations(__FILE__, __LINE__, SET, SET, RELATION);

    ct = new ct_entry_type("cross");
    ct->setFixed('I', arg1, arg2);
    ct->setResult(res);
    ct->doneBuilding();
}


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
        cp = compute_un(L, ap, bp);
    }
#ifdef TRACE
    out << "Cross #" << top_count << " end\n";
#endif
    cv.set();
}


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

    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    key[0].setI(k);
    key[1].setN(a);
    key[2].setN(b);
    if (ct->findCT(key, res)) {
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
    res[0].setN(c);
    ct->addCT(key, res);

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

    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    key[0].setI(k);
    key[1].setN(a);
    key[2].setN(b);
    if (ct->findCT(key, res)) {
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

    res[0].setN(c);
    ct->addCT(key, res);

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
// *                      CROSS_factory  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::CROSS_factory : public binary_factory {
    public:
        virtual void setup();
        virtual binary_operation* build_new(forest* a, forest* b, forest* c);
};

// ******************************************************************

void MEDDLY::CROSS_factory::setup()
{
    _setup(__FILE__, "CROSS", "Cross product of sets (boolean functions). The result is a relation. The input sets may be in different forests, but must have the same domain. The output forest must have the same domain, but be a relation.");
}

MEDDLY::binary_operation*
MEDDLY::CROSS_factory::build_new(forest* a, forest* b, forest* c)
{
    return new cross_bool(a, b, c);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_factory& MEDDLY::CROSS()
{
    static CROSS_factory F;
    return F;
}



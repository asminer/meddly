/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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

#include "../src/meddly.h"
#include <iostream>

using namespace MEDDLY;


void show_my_node(const char* where, unsigned who, const unpacked_node* un)
{
    ostream_output out(std::cout);

    out << "          " << where << " node " << char('A'+who) << ": ";
    un->show(out, false);
    out << "\n";
}


void compare(const unpacked_node* fnode, const unpacked_node* bnode)
{
    if (fnode->isFull() != bnode->isFull()) {
        throw "Compare: node storage mismatch";
    }

    if (fnode->isFull()) {
        unsigned i;
        for (i=0; i<fnode->getSize(); i++) {
            if (fnode->down(i) != bnode->down(i)) {
                throw "Full nodes DIFFER";
            }
        }
        for (; i<bnode->getSize(); i++) {
            if (bnode->down(i)) {
                throw "Full nodes DIFFER (extra in built)";
            }
        }
        std::cout << "      Full nodes MATCH\n";
        return;
    }

    // Both sparse
    if (fnode->getSize() != bnode->getSize()) {
        throw "Sparse nodes DIFFER (size mismatch)";
    }
    unsigned z;
    for (z=0; z<fnode->getSize(); z++) {
        if (fnode->down(z) != bnode->down(z)) {
            throw "Sparse nodes DIFFER";
        }
        if (fnode->index(z) != bnode->index(z)) {
            throw "Sparse nodes DIFFER";
        }
    }
    std::cout << "      Sparse nodes MATCH\n";
}


void compare(long ev, const unpacked_node* fnode, const unpacked_node* bnode)
{
    if (fnode->isFull() != bnode->isFull()) {
        throw "Compare: node storage mismatch";
    }

    long fev, bev;

    if (fnode->isFull()) {
        unsigned i;
        for (i=0; i<fnode->getSize(); i++) {
            if (fnode->down(i) != bnode->down(i)) {
                throw "Full nodes DIFFER";
            }
            fnode->edgeval(i).get(fev);
            bnode->edgeval(i).get(bev);
            //
            // Infinity special case
            //
            if (0==fnode->down(i)) {
                if ((0==fev) && (0==bev)) continue;
            }
            //
            // Normal edge
            //
            if (fev+ev != bev) {
                throw "Full nodes DIFFER (edge values)";
            }
        }
        for (; i<bnode->getSize(); i++) {
            if (bnode->down(i)) {
                throw "Full nodes DIFFER (extra in built)";
            }
        }
        std::cout << "      Full nodes MATCH\n";
        return;
    }

    // Both sparse
    if (fnode->getSize() != bnode->getSize()) {
        throw "Sparse nodes DIFFER (size mismatch)";
    }
    unsigned z;
    for (z=0; z<fnode->getSize(); z++) {
        if (fnode->down(z) != bnode->down(z)) {
            throw "Sparse nodes DIFFER";
        }
        if (fnode->index(z) != bnode->index(z)) {
            throw "Sparse nodes DIFFER";
        }
    }
    std::cout << "      Sparse nodes MATCH\n";
}

inline bool within_epsilon(double v1, double v2, bool show_delta)
{
    double delta;
    if (0.0 == v1) {
        delta = v2;
    } else if (0.0 == v2) {
        delta = v1;
    } else {
        delta = (v1-v2) / v1;
    }

    if (show_delta) {
        std::cerr << "Delta is " << delta << "\n";
    }

    if (delta > 1e-6) return false;
    return (delta > -1e-6);
}

void compare(double ev, const unpacked_node* fnode, const unpacked_node* bnode)
{
    if (fnode->isFull() != bnode->isFull()) {
        throw "Compare: node storage mismatch";
    }

    float fev, bev;

    if (fnode->isFull()) {
        unsigned i;
        for (i=0; i<fnode->getSize(); i++) {
            if (fnode->down(i) != bnode->down(i)) {
                throw "Full nodes DIFFER";
            }
            fnode->edgeval(i).get(fev);
            bnode->edgeval(i).get(bev);
            //
            // Infinity special case
            //
            if (0==fnode->down(i)) {
                if ((0==fev) && (0==bev)) continue;
            }
            //
            // Normal edge
            //
            if (!within_epsilon(fev * ev, bev, false)) {
                std::cerr << "fev*ev: " << fev*ev << "\n";
                std::cerr << "bev   : " << bev << "\n";
                within_epsilon(fev * ev, bev, true);
                throw "Full nodes DIFFER (edge values)";
            }
        }
        for (; i<bnode->getSize(); i++) {
            if (bnode->down(i)) {
                throw "Full nodes DIFFER (extra in built)";
            }
        }
        std::cout << "      Full nodes MATCH\n";
        return;
    }

    // Both sparse
    if (fnode->getSize() != bnode->getSize()) {
        throw "Sparse nodes DIFFER (size mismatch)";
    }
    unsigned z;
    for (z=0; z<fnode->getSize(); z++) {
        if (fnode->down(z) != bnode->down(z)) {
            throw "Sparse nodes DIFFER";
        }
        if (fnode->index(z) != bnode->index(z)) {
            throw "Sparse nodes DIFFER";
        }
    }
    std::cout << "      Sparse nodes MATCH\n";
}


unpacked_node* build_node(forest* f, unsigned who, bool full)
{
    dd_edge e[7];

    for (unsigned i=0; i<7; i++) {
        e[i].attach(f);
    }

    switch (f->getRangeType()) {
        case range_type::BOOLEAN:
                f->createEdge(true, e[0]);
                f->createEdge(true, e[1]);
                f->createEdge(true, e[2]);
                f->createEdge(true, e[3]);
                f->createEdge(true, e[4]);
                f->createEdge(true, e[5]);
                f->createEdge(true, e[6]);
                break;

        case range_type::INTEGER:
                f->createEdge(1L, e[0]);
                f->createEdge(2L, e[1]);
                f->createEdge(3L, e[2]);
                f->createEdge(44L, e[3]);
                f->createEdge(5L, e[4]);
                f->createEdge(666L, e[5]);
                f->createEdge(7777777L, e[6]);
                break;

        case range_type::REAL:
                f->createEdge(1.0f, e[0]);
                f->createEdge(2.2f, e[1]);
                f->createEdge(0.33333333f, e[2]);
                f->createEdge(4.4f, e[3]);
                f->createEdge(5.05f, e[4]);
                f->createEdge(6.0606f, e[5]);
                f->createEdge(0.007f, e[6]);
                break;
    }

    /*
        Strings of position, edge
     */
    const char* initstrs[] = {
        "10",
        "114283",
        "142536607182",
        "83",
        "04",
        "657680",
        "112233445566",
        "70"
    };
    if (who > 7) throw "Bad node index";

    const char* init = initstrs[who];

    unpacked_node* un = nullptr;

    if (full) {
        un = unpacked_node::newFull(f, who%4+1, 9);
        for (unsigned i=0; init[i]; i+=2) {
            unsigned pos = init[i]-'0';
            unsigned dwn = init[i+1]-'0';

            un->setFull(pos, e[dwn].getEdgeValue(), f->linkNode(e[dwn].getNode()));

        }
    } else {
        un = unpacked_node::newSparse(f, who%4+1, 9);
        unsigned nnzs = 0;
        for (unsigned i=0; init[i]; i+=2) {
            unsigned pos = init[i]-'0';
            unsigned dwn = init[i+1]-'0';

            un->setSparse(nnzs, pos, e[dwn].getEdgeValue(), f->linkNode(e[dwn].getNode()));

            ++nnzs;
        }
        un->shrink(nnzs);
    }

    return un;
}



void test_nodes(domain* d, range_type r, edge_labeling e, policies &p)
{
    forest* f = forest::create(d,  (e==edge_labeling::EVTIMES), r, e, p);

    node_handle nh;
    long int_ev;
    float real_ev;

    if (f->isMultiTerminal()) {

        for (unsigned w=0; w<8; w++) {
            unpacked_node* un = build_node(f, w, w<4);

            show_my_node("Built ", w, un);

            nh = f->createReducedNode(-1, un);

            unpacked_node* fn = f->newUnpacked(nh, FULL_OR_SPARSE);

            show_my_node("Pulled", w, fn);
            unpacked_node* bn = build_node(f, w, fn->isFull());
            show_my_node("Reblt ", w, bn);

            compare(fn, bn);
            unpacked_node::Recycle(fn);
            unpacked_node::Recycle(bn);
        }
    }

    if (f->isEVPlus()) {
        for (unsigned w=0; w<8; w++) {
            unpacked_node* un = build_node(f, w, w<4);

            show_my_node("Built ", w, un);

            f->createReducedNode(-1, un, int_ev, nh);

            unpacked_node* fn = f->newUnpacked(nh, FULL_OR_SPARSE);

            show_my_node("Pulled", w, fn);
            unpacked_node* bn = build_node(f, w, fn->isFull());
            show_my_node("Reblt ", w, bn);

            compare(int_ev, fn, bn);
            unpacked_node::Recycle(fn);
            unpacked_node::Recycle(bn);
        }
    }

    if (f->isEVTimes()) {
        for (unsigned w=0; w<8; w++) {
            unpacked_node* un = build_node(f, w, w<4);

            show_my_node("Built ", w, un);

            f->createReducedNode(-1, un, real_ev, nh);

            unpacked_node* fn = f->newUnpacked(nh, FULL_OR_SPARSE);

            show_my_node("Pulled", w, fn);
            unpacked_node* bn = build_node(f, w, fn->isFull());
            show_my_node("Reblt ", w, bn);

            compare(real_ev, fn, bn);
            unpacked_node::Recycle(fn);
            unpacked_node::Recycle(bn);
        }
    }

    forest::destroy(f);

}

int main()
{
    using namespace std;
    try {
        MEDDLY::initialize();

        const int bounds[] = { 9, 9, 9, 9, 9 };
        domain* d = domain::createBottomUp(bounds, 5);
        policies p;
        p.useDefaults(false);

        cout << "Testing full/sparse storage\n";

        p.setFullOrSparse();

        cout << "  Testing MT boolean\n";
        test_nodes(d, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL, p);
        cout << "  Testing MT integer\n";
        test_nodes(d, range_type::INTEGER, edge_labeling::MULTI_TERMINAL, p);
        cout << "  Testing MT real\n";
        test_nodes(d, range_type::REAL,    edge_labeling::MULTI_TERMINAL, p);
        cout << "  Testing EV+ integer\n";
        test_nodes(d, range_type::INTEGER, edge_labeling::EVPLUS,         p);
        cout << "  Testing EV* real\n";
        test_nodes(d, range_type::REAL,    edge_labeling::EVTIMES,        p);

        cout << "Passed\n";
        cout << "Testing full storage only\n";

        p.setFullStorage();

        cout << "  Testing MT boolean\n";
        test_nodes(d, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL, p);
        cout << "  Testing MT integer\n";
        test_nodes(d, range_type::INTEGER, edge_labeling::MULTI_TERMINAL, p);
        cout << "  Testing MT real\n";
        test_nodes(d, range_type::REAL,    edge_labeling::MULTI_TERMINAL, p);
        cout << "  Testing EV+ integer\n";
        test_nodes(d, range_type::INTEGER, edge_labeling::EVPLUS,         p);
        cout << "  Testing EV* real\n";
        test_nodes(d, range_type::REAL,    edge_labeling::EVTIMES,        p);

        cout << "Passed\n";
        cout << "Testing sparse storage only\n";

        p.setSparseStorage();

        cout << "  Testing MT boolean\n";
        test_nodes(d, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL, p);
        cout << "  Testing MT integer\n";
        test_nodes(d, range_type::INTEGER, edge_labeling::MULTI_TERMINAL, p);
        cout << "  Testing MT real\n";
        test_nodes(d, range_type::REAL,    edge_labeling::MULTI_TERMINAL, p);
        cout << "  Testing EV+ integer\n";
        test_nodes(d, range_type::INTEGER, edge_labeling::EVPLUS,         p);
        cout << "  Testing EV* real\n";
        test_nodes(d, range_type::REAL,    edge_labeling::EVTIMES,        p);

        cout << "Passed\n";

        MEDDLY::cleanup();
        return 0;
    }
    catch (MEDDLY::error e) {
        cerr << "Caught meddly error '" << e.getName() << "'\n";
        cerr << "    thrown in " << e.getFile() << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        cerr << "Caught our own error: " << e << "\n";
        return 2;
    }
    cerr << "Some other error?\n";
    return 3;
}

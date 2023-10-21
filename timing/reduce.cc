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
#include "timer.h"

#include <iostream>

// #define REPORTING

using namespace MEDDLY;

const unsigned BATCHSIZE = 16384;
const unsigned NUMBATCHES = 1024;
const unsigned TOTAL = BATCHSIZE * NUMBATCHES;

const unsigned VARSIZE = 10;

long seed = -1;

inline double Random()
{
    static const long MODULUS = 2147483647L;
    static const long MULTIPLIER = 48271L;
    static const long Q = MODULUS / MULTIPLIER;
    static const long R = MODULUS % MULTIPLIER;

    long t = MULTIPLIER * (seed % Q) - R * (seed / Q);
    if (t > 0) {
        seed = t;
    } else {
        seed = t + MODULUS;
    }
    return ((double) seed / MODULUS);
}

inline int Equilikely(int a, int b)
{
    return (a + (int) ((b - a + 1) * Random()));
}


void test_MT_full_Reductions(expert_forest* f, const char* what)
{
    dd_edge zero(f), one(f), two(f), x1(f), x2(f);

    f->createEdge(0L, zero);
    f->createEdge(1L, one);
    f->createEdge(2L, two);

    f->createEdgeForVar(1, false, x1);
    f->createEdgeForVar(2, false, x2);

    node_handle fixed[5];
    fixed[0] = zero.getNode();
    fixed[1] = one.getNode();
    fixed[2] = two.getNode();
    fixed[3] = x1.getNode();
    fixed[4] = x2.getNode();

    std::cerr << "    Building " << what << " MT nodes from full ";

    timer T;

    node_handle roots[BATCHSIZE];

    for (unsigned n=0; n<NUMBATCHES; n++) {
        if (0==n%64) std::cerr << '.';
        for (unsigned b=0; b<BATCHSIZE; b++) {
            unpacked_node* un = unpacked_node::New();
            un->bind_to_forest(f, 3, VARSIZE, true);
            for (unsigned i=0; i<VARSIZE; i++) {
                int d = Equilikely(0, 4);
                un->d_ref(i) = f->linkNode(fixed[d]);
            }
            roots[b] = f->createReducedNode(-1, un);
        } // for b
        for (unsigned b=0; b<BATCHSIZE; b++) {
            f->unlinkNode(roots[b]);
        }
    } // for n

    T.note_time();

    std::cerr << ' ' << T.get_last_seconds() << " seconds\n";

#ifdef REPORTING
    ostream_output out(std::cout);
    f->reportStats(out, "    ", ~0u);
#endif
}


void test_EV_full_Reductions(expert_forest* f, const char* what)
{
    dd_edge zero(f), x1(f), x2(f);

    f->createEdge(0L, zero);

    f->createEdgeForVar(1, false, x1);
    f->createEdgeForVar(2, false, x2);

    node_handle fixed[3];
    fixed[0] = zero.getNode();
    fixed[1] = x1.getNode();
    fixed[2] = x2.getNode();

    std::cerr << "    Building " << what << " EV nodes from full ";

    timer T;

    node_handle roots[BATCHSIZE];

    for (unsigned n=0; n<NUMBATCHES; n++) {
        if (0==n%64) std::cerr << '.';
        for (unsigned b=0; b<BATCHSIZE; b++) {
            unpacked_node* un = unpacked_node::New();
            un->bind_to_forest(f, 3, VARSIZE, true);
            for (unsigned i=0; i<VARSIZE; i++) {
                int d = Equilikely(0, 2);
                un->d_ref(i) = f->linkNode(fixed[d]);
                un->setEdge(i, long(Equilikely(0, 10)));
            }
            long x;
            f->createReducedNode(-1, un, x, roots[b]);
        } // for b
        for (unsigned b=0; b<BATCHSIZE; b++) {
            f->unlinkNode(roots[b]);
        }
    } // for n

    T.note_time();

    std::cerr << ' ' << T.get_last_seconds() << " seconds\n";

#ifdef REPORTING
    ostream_output out(std::cout);
    f->reportStats(out, "    ", ~0u);
#endif
}


int main()
{
    try {
        std::cerr   << "Timing experiments for building "
                    << TOTAL << " nodes, in batches of " << BATCHSIZE << "\n\n";

        MEDDLY::initialize();

        seed = 123456789;

        int bounds[3];
        bounds[0] = bounds[1] = bounds[2] = VARSIZE;

        domain *D = createDomainBottomUp(bounds, 3);
        forest* F;
        policies p;
        p.useDefaults(false);

        p.setFullStorage();
        F = D->createForest(false, range_type::INTEGER,
                        edge_labeling::MULTI_TERMINAL, p);
        test_MT_full_Reductions((expert_forest*) F, "  full");
        destroyForest(F);

        p.setSparseStorage();
        F = D->createForest(false, range_type::INTEGER,
                        edge_labeling::MULTI_TERMINAL, p);
        test_MT_full_Reductions((expert_forest*) F, "sparse");
        destroyForest(F);

        p.setFullStorage();
        F = D->createForest(false, range_type::INTEGER,
                        edge_labeling::EVPLUS, p);
        test_EV_full_Reductions((expert_forest*) F, "  full");
        destroyForest(F);

        p.setSparseStorage();
        F = D->createForest(false, range_type::INTEGER,
                        edge_labeling::EVPLUS, p);
        test_EV_full_Reductions((expert_forest*) F, "sparse");
        destroyForest(F);

        std::cerr << "\n";

        MEDDLY::cleanup();
        return 0;
    }

    catch (MEDDLY::error e) {
        std::cerr   << "\nCaught meddly error '" << e.getName()
                    << "'\n    thrown in " << e.getFile()
                    << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        std::cerr << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    std::cerr << "\nSome other error?\n";
    return 3;
}

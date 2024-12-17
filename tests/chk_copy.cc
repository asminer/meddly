
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

/*
    Exhaustively test the copy operator.
*/


#include <cstdlib>
#include <time.h>
#include <cassert>

#include "../src/meddly.h"

#define CHECK_ALL_MTMDDS
#define CHECK_ALL_MTMXDS
#define CHECK_ALL_INTMDDS
#define CHECK_ALL_REALMDDS
#define CHECK_ALL_REALMXDS

// #define DEBUG_BUILD

const int TERMS = 20;

using namespace MEDDLY;

long seed = -1;

double Random()
{
    const long MODULUS = 2147483647L;
    const long MULTIPLIER = 48271L;
    const long Q = MODULUS / MULTIPLIER;
    const long R = MODULUS % MULTIPLIER;

    long t = MULTIPLIER * (seed % Q) - R * (seed / Q);
    if (t > 0) {
        seed = t;
    } else {
        seed = t + MODULUS;
    }
    return ((double) seed / MODULUS);
}

int Equilikely(int a, int b)
{
    return (a + (int) ((b - a + 1) * Random()));
}

void randomizeMinterm(minterm &m, int max)
{
    if (m.isForRelations()) {
        for (unsigned i=1; i<=m.getNumVars(); i++) {
            int from = Equilikely(-1, max);
            int to = Equilikely(-2, max);
            if (DONT_CHANGE == to) {
                from = DONT_CARE;
            }
            m.setVars(i, from, to);
        }
    } else {
        for (unsigned i=1; i<=m.getNumVars(); i++) {
            m.setVar(i, Equilikely(-1, max));
        }
    }
}

void printMinterm(FILE* fout, const minterm& m)
{
    FILE_output out(fout);
    out << "    ";
    m.show(out);
    out << "\n";
}

void buildRandomFunc(long s, int terms, dd_edge &out, FILE* fout)
{
    seed = s;
    forest* f = out.getForest();
    minterm mint(f);

    switch (f->getRangeType()) {
        case range_type::BOOLEAN:
            f->createEdge(false, out);
            break;

        case range_type::INTEGER:
            f->createEdge(long(0), out);
            break;

        case range_type::REAL:
            f->createEdge(float(0), out);
            break;
    }
    if (fout) {
        fprintf(fout, "Function(s) based on %d terms:\n", terms);
    }
    for (int i=0; i<terms; i++) {
        randomizeMinterm(mint, 4);

        long i_value = Equilikely(1, 5);
        switch (f->getRangeType()) {
            case range_type::BOOLEAN:
                mint.setTerm(true);
                break;

            case range_type::INTEGER:
                mint.setTerm(i_value);
                break;

            case range_type::REAL:
                mint.setTerm(double(i_value));
                break;
        }

        dd_edge temp(out);
        mint.buildFunction(temp);

        if (fout) {
            if (f->getRangeType() != range_type::BOOLEAN) {
                fprintf(fout, "    value %ld\n", i_value);
            }
            printMinterm(fout, mint);
        }

        out += temp;

#ifdef DEBUG_BUILD
        FILE_output meddlyout(stdout);
        printf("new term: ");
        temp.show(meddlyout);
        printf("\nnew func: ");
        out.show(meddlyout);
        printf("\n");
#endif

    } // for i
}

inline void writeLongReduction(const forest* f)
{
    printf("%s", nameOf(f->getReductionRule()));
    /*
    switch(f->getReductionRule()) {

        case reduction_rule::FULLY_REDUCED:
            printf("fully-reduced ");
            break;

        case reduction_rule::QUASI_REDUCED:
            printf("quasi-reduced ");
            break;

        case reduction_rule::IDENTITY_REDUCED:
            printf("identity-reduced ");
            break;

        default:
            printf("unknown-reduced ");
            break;
    }
    */
}

void writeType(const forest* f)
{
    switch (f->getRangeType()) {
        case range_type::BOOLEAN:
            printf("bool ");
            break;

        case range_type::INTEGER:
            printf("int. ");
            break;

        case range_type::REAL:
            printf("real ");
            break;

        default:
            printf("unk. ");
            break;
    }
    switch(f->getReductionRule()) {

        case reduction_rule::FULLY_REDUCED:
            printf("FR ");
            break;

        case reduction_rule::QUASI_REDUCED:
            printf("QR ");
            break;

        case reduction_rule::IDENTITY_REDUCED:
            printf("IR ");
            break;

        default:
            printf("?R ");
            break;
    }

    switch (f->getEdgeLabeling()) {
        case edge_labeling::MULTI_TERMINAL:
            printf(" mt");
            break;

        case edge_labeling::EVPLUS:
            printf("ev+");
            break;

        case edge_labeling::INDEX_SET:
            printf("ind");
            break;

        case edge_labeling::EVTIMES:
            printf("ev*");
            break;

        default:
            printf(" ??");
            break;
    }

    if (f->isForRelations())  printf("mxd");
    else                      printf("mdd");
}

void testCopy(forest* srcF, forest* destF)
{
    printf("\t");
    writeType(srcF);
    printf(" -> ");
    writeType(destF);
    printf("      ");
    fflush(stdout);

    dd_edge srcE(srcF), dummy(srcF);
    dd_edge destE(destF);

    try {
        for (int t=1; t<=TERMS; t++) {

            long save_seed = seed;
#ifdef DEBUG_BUILD
            printf("\n\nBuilding source\n");
#endif
            buildRandomFunc(save_seed, t, srcE, nullptr);
#ifdef DEBUG_BUILD
            printf("\n\nBuilding destination\n");
#endif
            buildRandomFunc(save_seed, t, destE, nullptr);

            if (srcF->getRangeType() == range_type::BOOLEAN) {
                if (destF->getRangeType() == range_type::INTEGER) {
                    // convert destE to boolean
                    dd_edge zero(destF);
                    destF->createEdge(long(0), zero);
                    apply(NOT_EQUAL, destE, zero, destE);
                }
                if (destF->getRangeType() == range_type::REAL) {
                    // convert destE to boolean
                    dd_edge zero(destF);
                    destF->createEdge(float(0), zero);
                    apply(NOT_EQUAL, destE, zero, destE);
                }
            }

#ifdef DEBUG_BUILD
            printf("\n\nBuilding copy\n");
#endif

            dd_edge copyE(destF);
            apply(COPY, srcE, copyE);
            printf(".");
            fflush(stdout);

            if (copyE != destE) {
                FILE_output meddlyout(stdout);

                printf("failed!\n\n");

                printf("Source (first forest):\n");
                writeLongReduction(srcE.getForest());
                srcE.showGraph(meddlyout);

                printf("Destination (should get):\n");
                writeLongReduction(destE.getForest());
                destE.showGraph(meddlyout);

                printf("Copy (built from source):\n");
                writeLongReduction(copyE.getForest());
                copyE.showGraph(meddlyout);

                buildRandomFunc(save_seed, t, dummy, stdout);
                exit(1);
            }

        } // for t

        printf("OK\n");
    }
    catch (MEDDLY::error e) {
        printf("%s\n", e.getName());
    }

}



int processArgs(int argc, const char** argv)
{
    if (argc>2) {
        /* Strip leading directory, if any: */
        const char* name = argv[0];
        for (const char* ptr=name; *ptr; ptr++) {
            if ('/' == *ptr) name = ptr+1;
        }
        printf("Usage: %s <seed>\n", name);
        return 0;
    }
    if (argc>1) {
        seed = atol(argv[1]);
    }
    if (seed < 1) {
        seed = time(0);
    }

    printf("Using rng seed %ld\n", seed);
    return 1;
}

void addMDDforests(domain* D, forest** list, int &i,
        edge_labeling ev, range_type type)
{
    policies fr(false);
    fr.setFullyReduced();
    policies qr(false);
    qr.setQuasiReduced();

    list[i++] = forest::create(D, 0, type, ev, fr);
    list[i++] = forest::create(D, 0, type, ev, qr);
}

void addMXDforests(domain* D, forest** list, int &i,
        edge_labeling ev, range_type type)
{
    policies ir(true);
    ir.setIdentityReduced();
    policies fr(true);
    fr.setFullyReduced();
    policies qr(true);
    qr.setQuasiReduced();

    list[i++] = forest::create(D, 1, type, ev, ir);
    list[i++] = forest::create(D, 1, type, ev, fr);
    list[i++] = forest::create(D, 1, type, ev, qr);
}

int makeMTMDDforests(domain* D, forest** list)
{
    int i = 0;
    addMDDforests(D, list, i, edge_labeling::MULTI_TERMINAL, range_type::BOOLEAN);
    addMDDforests(D, list, i, edge_labeling::MULTI_TERMINAL, range_type::INTEGER);
    addMDDforests(D, list, i, edge_labeling::MULTI_TERMINAL, range_type::REAL);
    return i;
}

int makeMTMXDforests(domain* D, forest** list)
{
    int i = 0;
    addMXDforests(D, list, i, edge_labeling::MULTI_TERMINAL, range_type::BOOLEAN);
    addMXDforests(D, list, i, edge_labeling::MULTI_TERMINAL, range_type::INTEGER);
    addMXDforests(D, list, i, edge_labeling::MULTI_TERMINAL, range_type::REAL);
    return i;
}

int makeIntegerMDDforests(domain* D, forest** list)
{
    int i = 0;
    addMDDforests(D, list, i, edge_labeling::MULTI_TERMINAL, range_type::INTEGER);
    addMDDforests(D, list, i, edge_labeling::EVPLUS, range_type::INTEGER);
    return i;
}

void addRealMDDforests(int& i, domain* D, forest** list)
{
    addMDDforests(D, list, i, edge_labeling::MULTI_TERMINAL, range_type::REAL);
    // TBD - these are not supported yet
    //  addMDDforests(D, list, i, edge_labeling::EVPLUS, range_type::REAL);
    //  addMDDforests(D, list, i, edge_labeling::EVTIMES, range_type::REAL);
    // TBD
}

void addRealMXDforests(int& i, domain* D, forest** list)
{
    addMXDforests(D, list, i, edge_labeling::EVTIMES, range_type::REAL);
    addMXDforests(D, list, i, edge_labeling::MULTI_TERMINAL, range_type::REAL);
    // TBD - these are not supported yet
    //  addMXDforests(D, list, i, edge_labeling::EVPLUS, range_type::REAL);
    // TBD
}

int makeRealMDDforests(domain* D, forest** list)
{
    int i = 0;
    addRealMDDforests(i, D, list);
    return i;
}

int makeRealMXDforests(domain* D, forest** list)
{
    int i = 0;
    addRealMXDforests(i, D, list);
    return i;
}

void clearForests(forest** list, int N)
{
    for (int i=0; i<N; i++) {
        forest::destroy(list[i]);
        list[i] = 0;
    }
}

int main(int argc, const char** argv)
{
    if (!processArgs(argc, argv)) return 1;

    initialize();

    int vars[] = {5, 5, 5, 5, 5, 5};
    domain* myd = domain::createBottomUp(vars, 6);
    assert(myd);

    //
    // For later - arrays of forests (!)
    //

    forest* srcs[9];
    forest* dests[9];
    int slen = 9;
    int dlen = 9;
    for (int i=0; i<9; i++) {
        srcs[i] = 0;
        dests[i] = 0;
    }

    //
    // MTMDD tests, build all possible types of forests
    //

#ifdef CHECK_ALL_MTMDDS
    printf("Checking all possible copies between MTMDD forests\n");
    slen = makeMTMDDforests(myd, srcs);
    dlen = makeMTMDDforests(myd, dests);

    for (int i=0; i<slen; i++)
        for (int j=0; j<dlen; j++)
            testCopy(srcs[i], dests[j]);

    clearForests(srcs, slen);
    clearForests(dests, dlen);
    printf("\n");
#endif

    //
    // MTMXD tests, build all possible types
    //

#ifdef CHECK_ALL_MTMXDS
    printf("Checking all possible copies between MTMXD forests\n");
    slen = makeMTMXDforests(myd, srcs);
    dlen = makeMTMXDforests(myd, dests);

    for (int i=0; i<slen; i++)
        for (int j=0; j<dlen; j++)
            testCopy(srcs[i], dests[j]);

    clearForests(srcs, slen);
    clearForests(dests, dlen);
    printf("\n");
#endif

    //
    // MT and EV part 1
    // (integer to integer or real, mt or ev)
    //

#ifdef CHECK_ALL_INTMDDS
    printf("Checking all possible copies from integer to integer/real MDD forests\n");
    slen = makeIntegerMDDforests(myd, srcs);
    dlen = makeIntegerMDDforests(myd, dests);
    addRealMDDforests(dlen, myd, dests);

    for (int i=0; i<slen; i++) {
        for (int j=0; j<dlen; j++) {
            if (srcs[i]->isEVPlus() != dests[j]->isEVPlus()) {
                // Do not compare EV+MDD with others
                // Check the semantics of PLUS operation with EV+MDD
                printf("\tDo not compare ");
                writeType(srcs[i]);
                printf(" and ");
                writeType(dests[j]);
                printf("\n");
                fflush(stdout);
                continue;
            }
            testCopy(srcs[i], dests[j]);
        }
    }

    clearForests(srcs, slen);
    clearForests(dests, dlen);
    printf("\n");
#endif

    //
    // MT and EV part 2
    // (real only)
    //

#ifdef CHECK_ALL_REALMDDS
    printf("Checking all possible copies for real MDD forests\n");
    slen = makeRealMDDforests(myd, srcs);
    dlen = makeRealMDDforests(myd, dests);

    for (int i=0; i<slen; i++)
        for (int j=0; j<dlen; j++)
            testCopy(srcs[i], dests[j]);

    clearForests(srcs, slen);
    clearForests(dests, dlen);
    printf("\n");
#endif

    //
    // MT and EV part 3
    // (real only)
    //

#ifdef CHECK_ALL_REALMXDS
    printf("Checking all possible copies for real MXD forests\n");
    slen = makeRealMXDforests(myd, srcs);
    dlen = makeRealMXDforests(myd, dests);

    for (int i=0; i<slen; i++)
        for (int j=0; j<dlen; j++)
            testCopy(srcs[i], dests[j]);

    clearForests(srcs, slen);
    clearForests(dests, dlen);
    printf("\n");
#endif


    cleanup();
    return 0;
}



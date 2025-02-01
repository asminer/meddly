
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
    Tests the cross-product operator.
*/


#include <cstdlib>
#include <time.h>
#include <cassert>

#include "../src/meddly.h"
#include "randomize.h"

// #define DEBUG_RANDSET

int vars[] = {4, 4, 4, 4, 4, 4};

using namespace MEDDLY;

void randomizeMinterm(minterm &m, bool rows=true)
{
    if (m.isForSets()) {
        for (int i=1; i<=m.getNumVars(); i++) {
            m.setVar(i, vectorgen::Equilikely_I(-1, 3));
        }
    } else {
        if (rows) {
            for (int i=1; i<=m.getNumVars(); i++) {
                m.setVars(i, vectorgen::Equilikely_I(-1, 3), DONT_CARE);
            }
        } else {
            for (int i=1; i<=m.getNumVars(); i++) {
                m.setVars(i, DONT_CARE, vectorgen::Equilikely_I(-1, 3));
            }
        }
    }
#ifdef DEBUG_RANDSET
    FILE_output out(stdout);
    out << "Random minterm: ";
    m.show(out);
    out << "\n";
#endif
}


void makeRandomSet(forest* f, int nmt, dd_edge &x)
{
    dd_edge tmp(f);
    minterm mint(f);
    for (; nmt; nmt--) {
        randomizeMinterm(mint);
        mint.buildFunction(false, tmp);
        x += tmp;
    }
}

void makeRandomRows(forest* f, int nmt, dd_edge &x)
{
    dd_edge tmp(f);
    minterm mint(f);
    for (; nmt; nmt--) {
        randomizeMinterm(mint, true);
        mint.buildFunction(false, tmp);
        x += tmp;
    }
}

void makeRandomCols(forest* f, int nmt, dd_edge &x)
{
    dd_edge tmp(f);
    minterm mint(f);
    for (; nmt; nmt--) {
        randomizeMinterm(mint, false);
        mint.buildFunction(false, tmp);
        x += tmp;
    }
}

void test(forest* mdd, forest* mxd, int nmt)
{
    dd_edge rs(mdd), cs(mdd);
    dd_edge one(mdd);
    mdd->createConstant(true, one);

    dd_edge rr(mxd), cr(mxd), rcr(mxd), tmp(mxd);

    long saveseed = vectorgen::getSeed();
    makeRandomSet(mdd, nmt, rs);
    vectorgen::setSeed(saveseed, false);
    makeRandomRows(mxd, nmt, rr);

#ifdef DEBUG_RANDSET
    FILE_output out(stdout);
    printf("Generated random set:\n");
    rs.showGraph(out);
    printf("Generated random rows:\n");
    rr.showGraph(out);
#endif

    // check: generate rr from rs, make sure they match
    apply(CROSS, rs, one, tmp);
#ifdef DEBUG_RANDSET
    printf("rs x 1:\n");
    tmp.showGraph(out);
#endif
    if (tmp != rr) {
        FILE_output out(stdout);
        printf("mismatch on rs x 1\n");
        printf("rs:\n");
        rs.showGraph(out);
        printf("rr:\n");
        rr.showGraph(out);
        printf("rs x 1, should equal rr:\n");
        tmp.showGraph(out);
        throw "mismatch";
    }

    saveseed = vectorgen::getSeed();
    makeRandomSet(mdd, nmt, cs);
    vectorgen::setSeed(saveseed, false);
    makeRandomCols(mxd, nmt, cr);

#ifdef DEBUG_RANDSET
    printf("Generated random set:\n");
    cs.showGraph(out);
    printf("Generated random cols:\n");
    cr.showGraph(out);
#endif

    // check: generate cr from cs, make sure they match
    apply(CROSS, one, cs, tmp);
    if (tmp != cr) {
        FILE_output out(stdout);
        printf("mismatch on 1 x cs\n");
        printf("cs:\n");
        cs.showGraph(out);
        printf("cr:\n");
        cr.showGraph(out);
        printf("1 x cs, should equal cr:\n");
        tmp.showGraph(out);
        throw "mismatch";
    }

    // intersection of rr and cr should equal rs x cs.
    apply(CROSS, rs, cs, rcr);
    tmp = rr * cr;
    if (tmp != rcr) {
        FILE_output out(stdout);
        printf("mismatch on rs x cs\n");
        printf("rs:\n");
        rs.showGraph(out);
        printf("rr:\n");
        rr.showGraph(out);
        printf("cs:\n");
        cs.showGraph(out);
        printf("cr:\n");
        cr.showGraph(out);
        printf("rr * cr:\n");
        tmp.showGraph(out);
        printf("Cross product, should equal rr*cr:\n");
        rcr.showGraph(out);
        throw "mismatch";
    }
}


int processArgs(int argc, const char** argv)
{
    long seed = 0;
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

    vectorgen::setSeed(seed, true);
    return 1;
}

int main(int argc, const char** argv)
{
    if (!processArgs(argc, argv)) return 1;

    try {
        initialize();

        domain* myd = domain::createBottomUp(vars, 6);
        assert(myd);

        forest* mdd = forest::create(myd, SET, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);
        assert(mdd);
        forest* mxd = forest::create(myd, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);
        assert(mxd);

        for (int m=1; m<=20; m++) {
            printf("\tChecking cross-product for %2d random minterms\n", m);
            test(mdd, mxd, m);
        }
        domain::destroy(myd);
        cleanup();
        return 0;
    }
    catch (MEDDLY::error e) {
        std::cerr   << "\nCaught meddly error " << e.getName()
                    << "\n    thrown in " << e.getFile()
                    << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        std::cerr << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    std::cerr << "\nSome other error?\n";
    return 4;
}


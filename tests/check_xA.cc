
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
    Tests vector matrix multiplication.
*/


#include <string.h>
#include <cassert>

#include "../src/meddly.h"

#define SHOW_INDEXES
#define SHOW_MATRIX
// #define TEST_MATRIX_ELEMS

using namespace MEDDLY;

const char* JOB = nullptr;

inline void trying(const char* j)
{
    JOB = j;
}
inline void success()
{
    JOB = nullptr;
}

const int vars[] = {10, 10, 10};

const int R[] = {0, 3, 1, 4};
const int N[] = {0, 1, 4, 9};
const int S[] = {0, 2, 6, 5};

// State indexes:
//    R : 0
//    N : 1
//    S : 2

void build_oz(forest* indf, forest* mxd, dd_edge &ss, dd_edge &P)
{
    //
    // Build state indexes
    //
    trying("build indexes");
    minterm_coll sslist(3, indf);
    sslist.unused().setAll(R, 0L);
    sslist.pushUnused();
    sslist.unused().setAll(N, 1L);
    sslist.pushUnused();
    sslist.unused().setAll(S, 2L);
    sslist.pushUnused();
    rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);
    sslist.buildFunctionMin(infty, ss);
    success();

    FILE_output out(stdout);
#ifdef SHOW_INDEXES
    out << "Indexes (minterms):\n";
    sslist.show(out);
    out << "Indexes:\n";
    ss.showGraph(out);
#endif

    //
    // Build matrix elements
    //
    trying("build matrix");
    minterm_coll matrix(8, mxd);
    matrix.unused().setAll(R, R, 0.5f);
    matrix.pushUnused();
    matrix.unused().setAll(R, N, 0.25f);
    matrix.pushUnused();
    matrix.unused().setAll(R, S, 0.25f);
    matrix.pushUnused();

    matrix.unused().setAll(N, R, 0.5f);
    matrix.pushUnused();
    matrix.unused().setAll(N, S, 0.5f);
    matrix.pushUnused();

    matrix.unused().setAll(S, R, 0.25f);
    matrix.pushUnused();
    matrix.unused().setAll(S, N, 0.25f);
    matrix.pushUnused();
    matrix.unused().setAll(S, S, 0.5f);
    matrix.pushUnused();
    matrix.buildFunctionMax(0.0, P);
    success();

#ifdef TEST_MATRIX_ELEMS
    for (unsigned i=0; i<matrix.size(); i++) {
        out << "Building Mxd for minterm:\n";
        matrix.at(i).show(out);
        out << "\n";
        dd_edge ij(mxd);
        matrix.at(i).buildFunction(ij);
        ij.showGraph(out);
    }
#endif

#ifdef SHOW_MATRIX
    out << "Matrix (minterms):\n";
    matrix.show(out);
    out << "Matrix:\n";
    P.showGraph(out);
#endif
}

void by_hand_oz_xA(double* y, const double* const x)
{
    y[0] += x[0]/2 + x[1]/2 + x[2]/4;
    y[1] += x[0]/4 +          x[2]/4;
    y[2] += x[0]/4 + x[1]/2 + x[2]/2;
}

void by_hand_oz_Ax(double* y, const double* const x)
{
    y[0] += x[0]/2 + x[1]/4 + x[2]/4;
    y[1] += x[0]/2 +          x[2]/2;
    y[2] += x[0]/4 + x[1]/4 + x[2]/2;
}

bool equals(double *xgot, double *xans, int size)
{
    for (int i=0; i<size; i++) {
        double delta = xgot[i] - xans[i];
        if (xgot[i]) delta /= xgot[i];
        if (delta > -1e-5 && delta < 1e-5) continue;
        // not equal, print an error
        printf("Got     : [%lf", xgot[0]);
        for (int j=1; j<size; j++) printf(", %lf", xgot[j]);
        printf("]\nExpected: [%lf", xans[0]);
        for (int j=1; j<size; j++) printf(", %lf", xans[j]);
        printf("]\n");
        return false;
    }
    return true;
}

void readVector(const dd_edge &x, double* ans)
{
    trying("read vector");
    forest* f = x.getForest();
    rangeval temp;
    minterm Rm(f), Nm(f), Sm(f);
    Rm.setAll(R, 1.0);
    Nm.setAll(N, 1.0);
    Sm.setAll(S, 1.0);
    x.evaluate(Rm, temp);
    ans[0] = temp;
    x.evaluate(Nm, temp);
    ans[1] = temp;
    x.evaluate(Sm, temp);
    ans[2] = temp;
    success();
}

void impl_xA_check(dd_edge &x, const dd_edge &P)
{
    printf("xA multiplications (implicit):\n");

    forest* f = x.getForest();
    minterm init(f);
    init.setAll(N, 1.0f);
    init.buildFunction(0.0, x);

    for (int i=0; ; i++) {
        double ex[3];
        double exP[3];
        double test[3];
        readVector(x, ex);
        printf("p%d: [%lf, %lf, %lf]\n", i, ex[0], ex[1], ex[2]);
        if (i>=9) break;
        trying("multiply");
        apply(VM_MULTIPLY, x, P, x);
        success();

        // determine product by hand
        exP[0] = exP[1] = exP[2] = 0;
        by_hand_oz_xA(exP, ex);

        // Compare!
        readVector(x, test);
        if (!equals(test, exP, 3)) throw "vectors not equal";

    }
}

void expl_xA_check(const dd_edge &ss, const dd_edge &P)
{
    int i;
    double p[3];
    double q[3];
    double q_alt[3];
    p[0] = 0; p[1] = 1; p[2] = 0;
    printf("xA multiplications (explicit):\n");
    numerical_operation* VM = EXPLVECT_MATR_MULT(ss, P, ss);
    for (i=0; i<9; i++) {
        printf("p%d: [%lf, %lf, %lf]\n", i, p[0], p[1], p[2]);
        q[0] = q[1] = q[2] = 0;
        trying("multiply");
        VM->compute(q, p);
        success();

        // determine product by hand
        q_alt[0] = q_alt[1] = q_alt[2] = 0;
        by_hand_oz_xA(q_alt, p);

        // Compare!
        if (!equals(q, q_alt, 3)) throw "vectors not equal";
        p[0] = q[0];
        p[1] = q[1];
        p[2] = q[2];
    }
    printf("p%d: [%lf, %lf, %lf]\n", i, p[0], p[1], p[2]);
    operation::destroy(VM);
}

void impl_Ax_check(dd_edge &x, const dd_edge &P)
{
    printf("Ax multiplications (implicit):\n");

    forest* f = x.getForest();
    minterm init(f);
    init.setAll(N, 1.0f);
    init.buildFunction(0.0, x);

    for (int i=0; ; i++) {
        double ex[3];
        double exP[3];
        double test[3];
        readVector(x, ex);
        printf("p%d: [%lf, %lf, %lf]\n", i, ex[0], ex[1], ex[2]);
        if (i>=9) break;
        trying("multiply");
        apply(MV_MULTIPLY, P, x, x);
        success();

        // determine product by hand
        exP[0] = exP[1] = exP[2] = 0;
        by_hand_oz_Ax(exP, ex);

        // Compare!
        readVector(x, test);
        if (!equals(test, exP, 3)) throw "vectors not equal";

    }
}

void expl_Ax_check(const dd_edge &ss, const dd_edge &P)
{
    int i;
    double p[3];
    double q[3];
    double q_alt[3];
    p[0] = 0; p[1] = 1; p[2] = 0;
    printf("Ax multiplications (explicit):\n");
    numerical_operation* MV = MATR_EXPLVECT_MULT(ss, P, ss);
    for (i=0; i<9; i++) {
        printf("p%d: [%lf, %lf, %lf]\n", i, p[0], p[1], p[2]);
        q[0] = q[1] = q[2] = 0;
        trying("multiply");
        MV->compute(q, p);
        success();
        q_alt[0] = q_alt[1] = q_alt[2] = 0;
        by_hand_oz_Ax(q_alt, p);
        if (!equals(q, q_alt, 3)) throw "vectors not equal";
        p[0] = q[0];
        p[1] = q[1];
        p[2] = q[2];
    }
    printf("p%d: [%lf, %lf, %lf]\n", i, p[0], p[1], p[2]);
    operation::destroy(MV);
}

int main(int argc, const char** argv)
{
    try {
        initialize();

        domain* ozd = domain::createBottomUp(vars, 3);
        assert(ozd);
        forest* evpmdds = forest::create(ozd, SET, range_type::INTEGER,
                edge_labeling::EVPLUS);
        assert(evpmdds);
        forest* mtmxds = forest::create(ozd, RELATION, range_type::REAL,
                edge_labeling::MULTI_TERMINAL);
        assert(mtmxds);
        forest* mtmdds = forest::create(ozd, SET, range_type::REAL,
                edge_labeling::MULTI_TERMINAL);
        assert(mtmdds);

        dd_edge ss(evpmdds);
        dd_edge P(mtmxds);
        dd_edge x(mtmdds);

        build_oz(evpmdds, mtmxds, ss, P);

        expl_xA_check(ss, P);
        impl_xA_check(x, P);

        expl_Ax_check(ss, P);
        impl_Ax_check(x, P);

        // Avoid active node warning
        evpmdds->createConstant(0L, ss);
        mtmxds->createConstant(0.0, P);

        cleanup();
        return 0;
    }
    catch (MEDDLY::error e) {
        printf("\nCouldn't %s: %s\n\tThrown at %s line %u\n",
                JOB, e.getName(), e.getFile(), e.getLine()
        );
        return 1;
    }
    catch (const char* e) {
        printf("\nError: %s\n", e);
        return 2;
    }

    printf("\nUnknown error\n");
    return 3;
}

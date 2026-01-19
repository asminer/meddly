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
    Tests pre-image and post-image operations.
*/


#include <cstdlib>
#include <time.h>
#include <cassert>
#include <iostream>
#include <iomanip>

#include "../src/meddly.h"

#include "randomize.h"

// #define SHOW_INITIAL
// #define SHOW_GRAPH

vectorgen SG(MEDDLY::SET,      5, 3);
vectorgen RG(MEDDLY::RELATION, 5, 3);

const unsigned MIN_SET_CARD = 1;
const unsigned MAX_SET_CARD = 16;
const unsigned MULT_SET_CARD = 4;

const unsigned MIN_REL_CARD = 8;
const unsigned MAX_REL_CARD = 4096;
const unsigned MULT_REL_CARD = 8;

using namespace MEDDLY;

// Forest type, set in main()
//
//  'b' : multi terminal, boolean
//  'i' : multi terminal, integer
//  'p' : EV+, integer
//
char ftype;
MEDDLY::rangeval unreachable;

void usage(const char* arg0)
{
    /* Strip leading directory, if any: */
    const char* name = arg0;
    for (const char* ptr=arg0; *ptr; ptr++) {
        if ('/' == *ptr) name = ptr+1;
    }
    std::cerr << "\nUsage: " << name << "options seed\n\n";
    std::cerr << "Options:\n";
    std::cerr << "    --set:    Test on sets (default).\n";
    std::cerr << "    --rel:    Test on relations\n";
    std::cerr << "\n";
    std::cerr << "    --inF:    Input/output set/rel is fully-reduced\n";
    std::cerr << "    --inI:    Input/output set/rel is identity-reduced\n";
    std::cerr << "    --inQ:    Input/output set/rel is quasi-reduced\n";
    std::cerr << "\n";

    std::cerr << "    --relF:   Transition relation is fully-reduced\n";
    std::cerr << "    --relI:   Transition relation is identity-reduced\n";
    std::cerr << "    --relQ:   Transition relation is quasi-reduced\n";
    std::cerr << "\n";

    std::cerr << "    --MTb:    Input/output set/rel is MT with boolean range\n";
    std::cerr << "    --MTi:    Input/output set/rel is MT with integer range\n";
    std::cerr << "    --EVp:    Input/output set/rel is EV+ with integer range\n";
    std::cerr << "\n";

    exit(1);
}


struct intlist {
    unsigned index;
    intlist* next;

    intlist(unsigned i, intlist* n) {
        index = i;
        next = n;
    }
};

void showList(intlist* L)
{
    while (L) {
        std::cout << L->index;
        std::cout << " -> ";
        L = L->next;
    }
    std::cout << " null\n";
}

void makeExplicitGraph(const dd_edge &E, std::vector <intlist*> &elist)
{
    elist.resize(SG.potential());
    for (unsigned i=0; i<elist.size(); i++) {
        elist[i] = nullptr;
    }
#ifdef SHOW_GRAPH
    ostream_output out(std::cout);
    out << "Explicit graph:\n";
#endif
    for (dd_edge::iterator I = E.begin(); I != E.end(); I++) {
        unsigned f, t;
#ifdef SHOW_GRAPH
        out << "    ";
        (*I).show(out);
        out << "\n";
#endif
        RG.minterm2indexes(*I, f, t);
        // std::cout << "    " << f << " -> " << t << "\n";

        elist[f] = new intlist(t, elist[f]);
    }
}

void showExplicitGraph(const std::vector <intlist*> &elist)
{
    for (unsigned i=0; i<elist.size(); i++) {
        if (elist[i]) {
            std::cout << i << ": ";
            showList(elist[i]);
        }
    }
}

inline void dist_update(MEDDLY::rangeval &x, long d)
{
    if (x.isPlusInfinity()) {
        x = d;
        return;
    }
    if (long(x) < 0) {
        x = d;
        return;
    }
    if (d < long(x)) {
        x = d;
    }
}

void preImage_b(const std::vector <MEDDLY::rangeval> &s0,
        const std::vector <intlist*> &E, std::vector <MEDDLY::rangeval> &s1)
{
    if (s0.size() != E.size()) {
        throw "preImage size mismatch s0, E";
    }
    if (s1.size() != E.size()) {
        throw "preImage size mismatch s1, E";
    }
    for (unsigned i=0; i<s1.size(); i++) {
        s1[i] = unreachable;
    }
    for (unsigned i=0; i<E.size(); i++) {
        const intlist* curr = E[i];
        while (curr) {
            if (s0[curr->index]) {
                s1[i] = true;
                break;
            }
            curr = curr->next;
        }
    }
}

void preImage_i(const std::vector <MEDDLY::rangeval> &s0,
        const std::vector <intlist*> &E, std::vector <MEDDLY::rangeval> &s1)
{
    if (s0.size() != E.size()) {
        throw "preImage size mismatch s0, E";
    }
    if (s1.size() != E.size()) {
        throw "preImage size mismatch s1, E";
    }
    for (unsigned i=0; i<s1.size(); i++) {
        s1[i] = unreachable;
    }

    for (unsigned i=0; i<E.size(); i++) {
        const intlist* curr = E[i];
        while (curr) {
            if (s0[curr->index] != unreachable) {
                dist_update(s1[i], long(s0[curr->index]) + 1);
            }
            curr = curr->next;
        }
    }
}

void postImage_b(const std::vector <MEDDLY::rangeval> &s0,
        const std::vector <intlist*> &E, std::vector <MEDDLY::rangeval> &s1)
{
    if (s0.size() != E.size()) {
        throw "postImage size mismatch s0, E";
    }
    if (s1.size() != E.size()) {
        throw "postImage size mismatch s1, E";
    }
    for (unsigned i=0; i<s1.size(); i++) {
        s1[i] = unreachable;
    }
    for (unsigned i=0; i<E.size(); i++) {
        if (!s0[i]) continue;
        const intlist* curr = E[i];
        while (curr) {
            s1[curr->index] = true;
            curr = curr->next;
        }
    }
/*
    unsigned card = 0;
    for (unsigned i=0; i<s1.size(); i++) {
        if (s1[i]) card++;
    }
    std::cout << " card " << card << "\n";
    */
}

void postImage_i(const std::vector <MEDDLY::rangeval> &s0,
        const std::vector <intlist*> &E, std::vector <MEDDLY::rangeval> &s1)
{
    if (s0.size() != E.size()) {
        throw "postImage size mismatch s0, E";
    }
    if (s1.size() != E.size()) {
        throw "postImage size mismatch s1, E";
    }
    for (unsigned i=0; i<s1.size(); i++) {
        s1[i] = unreachable;
    }
    for (unsigned i=0; i<E.size(); i++) {
        if (unreachable == s0[i]) continue;
        long newdist = long(s0[i]) + 1;
        const intlist* curr = E[i];
        while (curr) {
            dist_update(s1[curr->index], newdist);
            curr = curr->next;
        }
    }
/*
    unsigned card = 0;
    for (unsigned i=0; i<s1.size(); i++) {
        if (s1[i]) card++;
    }
    std::cout << " card " << card << "\n";
    */
}


void checkEqual(const char* what, const dd_edge &e1, const dd_edge &e2)
{
    if (e1 == e2) return;

    ostream_output out(std::cout);

    out << "\nMismatch on " << what << "\n";
    out << "Expected DD:\n";
    e2.showGraph(out);
    out << "Obtained DD:\n";
    e1.showGraph(out);

    out << "\nExpected minterms:\n";
    for (dd_edge::iterator s = e2.begin(); s; ++s)
    {
        out << "    ";
        (*s).show(out);
        out << "\n";
    }

    out << "Obtained minterms:\n";
    for (dd_edge::iterator s = e1.begin(); s; ++s)
    {
        out << "    ";
        (*s).show(out);
        out << "\n";
    }


    throw "mismatch";
}

void setCompare(const std::vector <MEDDLY::rangeval> &S0,
        const dd_edge &ddR, forest* Fset)
{
    std::vector <MEDDLY::rangeval> post(S0.size());
    std::vector <MEDDLY::rangeval> pre(S0.size());

    std::vector <intlist*> R;
    makeExplicitGraph(ddR, R);

    if ('b' == ftype) {
        postImage_b(S0, R, post);
        preImage_b(S0, R, pre);
    } else {
        postImage_i(S0, R, post);
        preImage_i(S0, R, pre);
    }

    dd_edge dd0(Fset);
    if ('p' == ftype) {
        SG.explicit2edgeMin(S0, dd0, unreachable);
    } else {
        SG.explicit2edgeMax(S0, dd0, unreachable);
    }
#ifdef SHOW_INITIAL
    std::cout << "Initial minterm set:\n";
    SG.showMinterms(std::cout, Fset->getDomain(), S0, unreachable);
    std::cout << "\n";
#endif

    //
    // Check post-image
    //
    dd_edge ddpost(Fset), sympost(Fset);
    if ('p' == ftype) {
        SG.explicit2edgeMin(post, ddpost, unreachable);
    } else {
        SG.explicit2edgeMax(post, ddpost, unreachable);
    }
    apply(POST_IMAGE, dd0, ddR, sympost);
    checkEqual("post image", sympost, ddpost);

    //
    // Check pre-image
    //
    dd_edge ddpre(Fset), sympre(Fset);
    if ('p' == ftype) {
        SG.explicit2edgeMin(pre, ddpre, unreachable);
    } else {
        SG.explicit2edgeMax(pre, ddpre, unreachable);
    }
    apply(PRE_IMAGE, dd0, ddR, sympre);
    checkEqual("pre image", sympre, ddpre);
}

void setTestsOnForests(unsigned scard, forest* Fset,
        unsigned rcard, forest* Frel)
{
    if (!Fset) throw "null Fset";
    if (!Frel) throw "null frel";

    std::vector <MEDDLY::rangeval> expset(SG.potential());
    std::vector <MEDDLY::rangeval> exprel(RG.potential());
    dd_edge Rel(Frel);

    std::vector <MEDDLY::rangeval> rgvals(2);
    rgvals[0] = false;
    rgvals[1] = true;

    std::vector <MEDDLY::rangeval> values;

    if (Fset->isRangeType(MEDDLY::range_type::BOOLEAN)) {
        values.resize(2);
        values[0] = unreachable;
        values[1] = true;
    } else {
        values.resize(4);
        values[0] = unreachable;
        values[1] = 0;
        values[2] = 1;
        values[3] = 2;
    }

    const unsigned N = 8;

    std::cerr << "    " << std::setw(5) << scard << " x "
              << std::setw(5) << rcard << " ";

    for (unsigned i=0; i<N; i++) {
        std::cerr << '.';
        RG.randomizeVector(exprel, rcard, rgvals);
        RG.explicit2edgeMax(exprel, Rel, false);

        SG.randomizeVector(expset, scard, values);
        setCompare(expset, Rel, Fset);

        SG.randomizeFully(expset, scard, values);
        setCompare(expset, Rel, Fset);
    }
    for (unsigned i=0; i<N; i++) {
        std::cerr << "x";
        RG.randomizeFully(exprel, rcard, rgvals);
        RG.explicit2edgeMax(exprel, Rel, false);

        SG.randomizeVector(expset, scard, values);
        setCompare(expset, Rel, Fset);

        SG.randomizeFully(expset, scard, values);
        setCompare(expset, Rel, Fset);
    }
    for (unsigned i=0; i<N; i++) {
        std::cerr << "i";
        RG.randomizeIdentity(exprel, rcard, rgvals);
        RG.explicit2edgeMax(exprel, Rel, false);

        SG.randomizeVector(expset, scard, values);
        setCompare(expset, Rel, Fset);

        SG.randomizeFully(expset, scard, values);
        setCompare(expset, Rel, Fset);
    }
    std::cerr << std::endl;
}


void testSets(domain* D, range_type RT, edge_labeling EL, const policies &sp,
        const policies &rp)
{
    forest* Fset = forest::create(D, SET, RT, EL, sp);

    forest* Frel = forest::create(D, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, rp);


    std::cerr << "Input/output forest:\n"
              << "    " << nameOf(EL) << " " << nameOf(RT)
              << " " << nameOf(Fset->getReductionRule()) << "\n";

    std::cerr << "Transition relation forest:\n"
              << "    " << nameOf(Frel->getReductionRule()) << "\n";


    std::cerr << "Running tests for input cardinality x transition relation cardinality\n";

    for (unsigned rc=MIN_REL_CARD; rc<=MAX_REL_CARD; rc*=MULT_REL_CARD)
    {
        for (unsigned sc=MIN_SET_CARD; sc<=MAX_SET_CARD; sc*=MULT_SET_CARD)
        {
            setTestsOnForests(sc, Fset, rc, Frel);
        }
    }

}

void testRels(domain* D, range_type RT, edge_labeling EL, policies sp, policies rp)
{
}

int main(int argc, const char** argv)
{
    /*
     * Command-line options
     */
    bool sets = true;
    long seed = 0;
    char setred = 'F';
    char relred = 'I';
    ftype = 'b';
    edge_labeling EL = edge_labeling::MULTI_TERMINAL;
    range_type    RT = range_type::BOOLEAN;

    for (int i=1; i<argc; i++) {

        if (0==strcmp("--set", argv[i])) {
            sets = true;
            continue;
        }
        if (0==strcmp("--rel", argv[i])) {
            sets = false;
            continue;
        }

        if (0==strcmp("--inF", argv[i])) {
            setred = 'F';
            continue;
        }
        if (0==strcmp("--inI", argv[i])) {
            setred = 'I';
            continue;
        }
        if (0==strcmp("--inQ", argv[i])) {
            setred = 'Q';
            continue;
        }

        if (0==strcmp("--relF", argv[i])) {
            relred = 'F';
            continue;
        }
        if (0==strcmp("--relI", argv[i])) {
            relred = 'I';
            continue;
        }
        if (0==strcmp("--relQ", argv[i])) {
            relred = 'Q';
            continue;
        }

        if (0==strcmp("--MTb", argv[i])) {
            EL = edge_labeling::MULTI_TERMINAL;
            RT = range_type::BOOLEAN;
            ftype = 'b';
            continue;
        }
        if (0==strcmp("--MTi", argv[i])) {
            EL = edge_labeling::MULTI_TERMINAL;
            RT = range_type::INTEGER;
            ftype = 'i';
            continue;
        }
        if (0==strcmp("--EVp", argv[i])) {
            EL = edge_labeling::EVPLUS;
            RT = range_type::INTEGER;
            ftype = 'p';
            continue;
        }

        if ((argv[i][0] < '0') || (argv[i][0] > '9')) {
            usage(argv[0]);
        }

        seed = atol(argv[i]);
    }
    vectorgen::setSeed(seed);

    /*
     * Full unreachable value
     */
    switch (ftype) {
        case 'p':
            unreachable = MEDDLY::rangeval(MEDDLY::range_special::PLUS_INFINITY,
                            MEDDLY::range_type::INTEGER);
            break;

        case 'i':
            unreachable = -1;
            break;

        default:
            unreachable = false;
    }

    /*
     * Run requested test
     */

    try {
        MEDDLY::initialize();
        domain* D = RG.makeDomain();

        policies sp;
        policies rp;

        sp.useDefaults(sets ? SET : RELATION);
        rp.useDefaults(RELATION);

        switch (setred) {
            case 'Q':
                sp.setQuasiReduced();
                break;

            case 'F':
                sp.setFullyReduced();
                break;

            case 'I':
                sp.setIdentityReduced();
                break;

            default:
                throw "internal error 1";
        }
        switch (relred) {
            case 'Q':
                rp.setQuasiReduced();
                break;

            case 'F':
                rp.setFullyReduced();
                break;

            case 'I':
                rp.setIdentityReduced();
                break;

            default:
                throw "internal error 2";
        }


        if (sets) {
            testSets(D, RT, EL, sp, rp);
        } else {
            testRels(D, RT, EL, sp, rp);
        }

        domain::destroy(D);
        MEDDLY::cleanup();
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



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

#include <iostream>
#include <iomanip>

#include "../src/meddly.h"

using namespace MEDDLY;

const int VARSIZE = 2;

// #define DEBUG

class myop : public MEDDLY::operation {
    public:
        myop(forest* f1, forest* f2);

        virtual bool checkForestCompatibility() const {
            return true;
        }

        inline compute_table* ct0() {
            return CT[0];
        }
        inline compute_table* ct1() {
            return CT[1];
        }
        inline compute_table* ct2() {
            return CT[2];
        }

    public:
        ct_entry_type* et0;
        ct_entry_type* et1;
        ct_entry_type* et2;
};

myop::myop(forest* F1, forest* F2) : operation("myop", 3)
{
    et0 = new ct_entry_type("test0", "IN:I");
    et0->setForestForSlot(1, F1);
    registerEntryType(0, et0);

    et1 = new ct_entry_type("test1", "NN:I");
    et1->setForestForSlot(0, F1);
    et1->setForestForSlot(1, F1);
    registerEntryType(1, et1);

    et2 = new ct_entry_type("test2", "NN:I");
    et2->setForestForSlot(0, F1);
    et2->setForestForSlot(1, F2);
    registerEntryType(2, et2);

    buildCTs();
}

void initEdges(forest* f, std::vector <dd_edge> &E)
{
    ostream_output out(std::cout);

    long terms[VARSIZE];
    terms[0] = 0;
    for (unsigned i=0; i<E.size(); i++) {
        terms[1] = i+1;
        E[i].attach(f);
        f->createEdgeForVar(1, false, terms, E[i]);
#ifdef DEBUG
        out << "Built edge #" << i << ": " << E[i] << "\n";
#endif
    }
}

//
// Add entries of type  IN:I
// returns the number of entries added.
//
unsigned addEntries(const ct_entry_type* CTE, compute_table* CT,
        const std::vector <dd_edge> &E)
{
    unsigned cnt=0;
    ct_entry_result res;
    res.initialize(CTE);
    for (int i=0; i<E.size(); i++) {
        for (unsigned j=0; j<E.size(); j++) {
            ct_entry_key* key = compute_table::useEntryKey(CTE, 0);

            key->writeI(i);
            key->writeN(E[j].getNode());

            CT->find(key, res);
            if (res) throw "Found result in CT; shouldn't have.";

            res.reset();
            res.writeI(i+E[j].getNode());
            CT->addEntry(key, res);
            ++cnt;
        } // for j
    } // for i

    return cnt;
}

//
// Check entries of type  IN:I
// Returns the number of CT hits
//
unsigned checkEntries(const ct_entry_type* CTE, compute_table* CT,
        const std::vector <dd_edge> &E)
{
    unsigned hits = 0;
    ct_entry_result res;
    res.initialize(CTE);
    for (int i=0; i<E.size(); i++) {
        for (unsigned j=0; j<E.size(); j++) {
            ct_entry_key* key = compute_table::useEntryKey(CTE, 0);

            key->writeI(i);
            key->writeN(E[j].getNode());

            CT->find(key, res);
            CT->recycle(key);

            if (!res) continue;
            int answer = res.readI();
            if (answer != i + E[j].getNode()) {
                std::cerr << "Wrong CT entry.\n";
                std::cerr << "    got     : " << i << ", " << E[j].getNode()
                          << " : " << answer << "\n";
                std::cerr << "    expected: " << i << ", " << E[j].getNode()
                          << " : " << i+E[j].getNode() << "\n";
                throw "cache error";
            }
            hits++;

        } // for j
    } // for i

    return hits;
}

//
// Add entries of type  NN:I
//
unsigned addEntries(const ct_entry_type* CTE, compute_table* CT,
        const std::vector <dd_edge> &E1, const std::vector <dd_edge> &E2)
{
    unsigned cnt = 0;
    ct_entry_result res;
    res.initialize(CTE);
    for (unsigned i=0; i<E1.size(); i++) {
        for (unsigned j=0; j<E2.size(); j++) {
            ct_entry_key* key = compute_table::useEntryKey(CTE, 0);

            key->writeN(E1[i].getNode());
            key->writeN(E2[j].getNode());

            CT->find(key, res);
            if (res) throw "Found result in CT; shouldn't have.";

            res.reset();
            res.writeI(E1[i].getNode() + E2[j].getNode());
            CT->addEntry(key, res);

            ++cnt;
        } // for j
    } // for i

    return cnt;
}

//
// Check entries of type  NN:I
// Returns the number of CT hits
//
unsigned checkEntries(const ct_entry_type* CTE, compute_table* CT,
        const std::vector <dd_edge> &E1, const std::vector <dd_edge> &E2)
{
    unsigned hits = 0;
    ct_entry_result res;
    res.initialize(CTE);
    for (unsigned i=0; i<E1.size(); i++) {
        for (unsigned j=0; j<E2.size(); j++) {
            ct_entry_key* key = compute_table::useEntryKey(CTE, 0);

            key->writeN(E1[i].getNode());
            key->writeN(E2[j].getNode());

            CT->find(key, res);
            CT->recycle(key);

            if (!res) continue;
            int answer = res.readI();
            if (answer != E1[i].getNode() + E2[j].getNode()) {
                std::cerr << "Wrong CT entry.\n";
                std::cerr << "    got     : " << E1[i].getNode()
                          << ", " << E2[j].getNode()
                          << " : " << answer << "\n";
                std::cerr << "    expected: " << E1[i].getNode()
                          << ", " << E2[j].getNode()
                          << " : " << E1[i].getNode()+E2[j].getNode() << "\n";
                throw "cache error";
            }
            hits++;

        } // for j
    } // for i

    return hits;
}

//
// Run CT tests.
//
void check_CT(bool monolithic)
{
    using namespace std;

    //
    // Set up domain
    //
    int bounds[3];
    bounds[0] = bounds[1] = bounds[2] = VARSIZE;
    domain *D = domain::createBottomUp(bounds, 3);

    //
    // Set up two (identical) forests
    //
    forest *F1, *F2;
    policies p;
    p.useDefaults(false);
    F1 = forest::create(D, false, range_type::INTEGER,
                        edge_labeling::MULTI_TERMINAL, p);
    F2 = forest::create(D, false, range_type::INTEGER,
                        edge_labeling::MULTI_TERMINAL, p);

    //
    // Make some dd_edges
    //
    vector <dd_edge> E1(100);
    vector <dd_edge> E2(100);

    initEdges(F1, E1);
    initEdges(F2, E2);

    //
    // Add some CT entries
    //

    unsigned add, entries1, entries2, entries3, hits1, hits2, hits3;
    myop* foo = new myop(F1, F2);

    cout << "    Creating CT entries\n";
    add = addEntries(foo->et0, foo->ct0(), E1);
    cout << "        " << setw(5) << add << " op1 entries added\n";

    add = addEntries(foo->et1, foo->ct1(), E1, E1);
    cout << "        " << setw(5) << add << " op2 entries added\n";

    add = addEntries(foo->et2, foo->ct2(), E1, E2);
    cout << "        " << setw(5) << add << " op3 entries added\n";

    entries1 = foo->ct0()->getStats().numEntries;
    entries2 = foo->ct1()->getStats().numEntries;
    entries3 = foo->ct2()->getStats().numEntries;

    if (monolithic) {
        cout << "        " << setw(5) << entries1 << " entries in CT\n";
    } else {
        cout << "        " << setw(5) << entries1 << " entries in op1 CT\n";
        cout << "        " << setw(5) << entries2 << " entries in op2 CT\n";
        cout << "        " << setw(5) << entries3 << " entries in op3 CT\n";
    }

    cout << "    Checking CT entries\n";

    hits1 = checkEntries(foo->et0, foo->ct0(), E1);
    cout << "        " << setw(5) << hits1 << " op1 entries matched\n";
    hits2 = checkEntries(foo->et1, foo->ct1(), E1, E1);
    cout << "        " << setw(5) << hits2 << " op2 entries matched\n";
    hits3 = checkEntries(foo->et2, foo->ct2(), E1, E2);
    cout << "        " << setw(5) << hits3 << " op3 entries matched\n";

    if (monolithic) {
        cout << "        " << setw(5) << hits1+hits2+hits3
             << " total entries matched\n";
    }

    if (monolithic) {
        if (hits1+hits2+hits3 != entries1) {
            throw "monolithic CT hit count mismatch";
        }
    } else {
        if (hits1 != entries1) {
            throw "CT1 hit count mismatch";
        }
        if (hits2 != entries2) {
            throw "CT2 hit count mismatch";
        }
        if (hits3 != entries3) {
            throw "CT3 hit count mismatch";
        }
    }


    //
    // Destroy forest 2 and re-count!
    // TBD
    //

}

int main(int argc, const char** argv)
{
    using namespace std;
    try {

        // go through all CT styles
        //
        MEDDLY::initialize();

        check_CT(true);

        MEDDLY::cleanup();

        cout << "Done\n";
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
    return 4;
}

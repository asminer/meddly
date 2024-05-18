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

// const unsigned EDGES = 3;
const unsigned EDGES = 100;

// #define DEBUG

// #define DEBUG_CT

inline node_handle MAX(node_handle a, node_handle b)
{
    return (a>b) ? a : b;
}

class myop : public MEDDLY::operation {
    public:
        myop(forest* f1, forest* f2);

        virtual bool checkForestCompatibility() const {
            return true;
        }

    public:
        ct_entry_type* et0;
        ct_entry_type* et1;
        ct_entry_type* et2;
};

myop::myop(forest* F1, forest* F2) : operation("myop", 3)
{
    et0 = new ct_entry_type("op1", "IN:I");
    et0->setForestForSlot(1, F1);
    registerEntryType(0, et0);

    et1 = new ct_entry_type("op2", "NN:N");
    et1->setForestForSlot(0, F1);
    et1->setForestForSlot(1, F1);
    et1->setForestForSlot(3, F1);
    registerEntryType(1, et1);

    et2 = new ct_entry_type("op3", "NN:N");
    et2->setForestForSlot(0, F1);
    et2->setForestForSlot(1, F2);
    et2->setForestForSlot(3, F2);
    registerEntryType(2, et2);

    buildCTs();

    /*
    ostream_output out(std::cout);
    out << "Built entry types:\n";
    et0->show(out);
    et1->show(out);
    et2->show(out);
    */
}

void initEdges(forest* f, std::vector <dd_edge> &E)
{
    ostream_output out(std::cout);

    long terms[VARSIZE];
    terms[0] = 0;
    for (unsigned i=0; i<E.size(); i++) {
        E[i].attach(f);
        if (i) {
            terms[1] = i+1;
            f->createEdgeForVar(1, false, terms, E[i]);
        } else {
            f->createEdge(42L, E[i]);
        }
#ifdef DEBUG
        out << "Built edge #" << i << ": " << E[i] << "\n";
#endif
    }
}

//
// Add entries of type  IN:I
// returns the number of entries added.
//
unsigned addEntries(ct_entry_type* CTE, const std::vector <node_handle> &N)
{
    unsigned cnt=0;
    ct_entry_result res;
    res.initialize(CTE);
    for (int i=0; i<N.size(); i++) {
        for (unsigned j=0; j<N.size(); j++) {
            ct_entry_key* key = compute_table::useEntryKey(CTE, 0);

            key->writeI(i);
            key->writeN(N[j]);

            CTE->findCT(key, res);
            if (res) throw "Found result in CT; shouldn't have.";

            res.reset();
            res.writeI(i+N[j]);
            CTE->addCT(key, res);
            ++cnt;
        } // for j
    } // for i

    return cnt;
}

//
// Check entries of type  IN:I
//
unsigned checkEntries(const char* name, ct_entry_type* CTE, const std::vector <node_handle> &N)
{
    using namespace std;
    unsigned hits = 0;
    ct_entry_result res;
    res.initialize(CTE);
    for (int i=0; i<N.size(); i++) {
        for (unsigned j=0; j<N.size(); j++) {
            ct_entry_key* key = compute_table::useEntryKey(CTE, 0);

            key->writeI(i);
            key->writeN(N[j]);

            CTE->findCT(key, res);
            CTE->noaddCT(key);

            if (!res) continue;
            int answer = res.readI();
            if (answer != i + N[j]) {
                std::cerr << "Wrong CT entry.\n";
                std::cerr << "    got     : " << i << ", " << N[j]
                          << " : " << answer << "\n";
                std::cerr << "    expected: " << i << ", " << N[j]
                          << " : " << i+N[j] << "\n";
                throw "cache error";
            }
            hits++;

        } // for j
    } // for i

    cout << "        " << setw(5) << hits << " " << name << " entries matched\n";
    const unsigned actual = CTE->getNumEntries();

    cout << "        " << setw(5) << actual << " " << name << " entries present\n";

    if (hits != actual) throw "hit count mismatch";

    return hits;
}

//
// Add entries of type  NN:N
//
unsigned addEntries(ct_entry_type* CTE,
        const std::vector <node_handle> &N1,
        const std::vector <node_handle> &N2)
{
    unsigned cnt = 0;
    ct_entry_result res;
    res.initialize(CTE);
    for (unsigned i=0; i<N1.size(); i++) {
        for (unsigned j=0; j<N2.size(); j++) {
            ct_entry_key* key = compute_table::useEntryKey(CTE, 0);

            key->writeN(N1[i]);
            key->writeN(N2[j]);

            CTE->findCT(key, res);
            if (res) throw "Found result in CT; shouldn't have.";

            res.reset();
            res.writeN(MAX(N1[i],N2[j]));
            CTE->addCT(key, res);

            ++cnt;
        } // for j
    } // for i

    return cnt;
}

//
// Check entries of type  NN:N
// Returns the number of CT hits
//
unsigned checkEntries(const char* name, ct_entry_type* CTE,
        const std::vector <node_handle> &N1,
        const std::vector <node_handle> &N2)
{
    using namespace std;

    unsigned hits = 0;
    ct_entry_result res;
    res.initialize(CTE);
    for (unsigned i=0; i<N1.size(); i++) {
        for (unsigned j=0; j<N2.size(); j++) {
            ct_entry_key* key = compute_table::useEntryKey(CTE, 0);

            key->writeN(N1[i]);
            key->writeN(N2[j]);

            CTE->findCT(key, res);
            CTE->noaddCT(key);

            if (!res) continue;
            node_handle answer = res.readN();
            node_handle max = MAX(N1[i], N2[j]);
            if (answer != max) {
                std::cerr << "Wrong CT entry.\n";
                std::cerr << "    got     : " << N1[i]
                          << ", " << N2[j]
                          << " : " << answer << "\n";
                std::cerr << "    expected: " << N1[i]
                          << ", " << N2[j]
                          << " : " << max << "\n";
                throw "cache error";
            }
            hits++;

        } // for j
    } // for i

    cout << "        " << setw(5) << hits << " " << name << " entries matched\n";
    const unsigned actual = CTE->getNumEntries();

    cout << "        " << setw(5) << actual << " " << name << " entries present\n";

    if (hits != actual) throw "hit count mismatch";
    return hits;
}

inline void show_entries(unsigned entries, unsigned adds, const char* ctname)
{
    unsigned retain = entries*1000 / adds;
    std::cout   << "        " << std::setw(5) << entries << " entries in "
                << ctname << " (" << retain/10 << "." << retain%10
                << "% retention)\n";

    if (retain < 800) throw "not enough hits?";
}

inline void check_new_hits(unsigned oldhits, unsigned newhits)
{
    if (oldhits == newhits) return;
    if (newhits > oldhits) {
        std::cout << "              too many hits?\n";
        throw "old vs new mismatch";
    }
    // we can drop, but not by too many
    unsigned retain = newhits * 100 / oldhits;
    std::cout << "              (" << retain << "% retention)\n";

    if (retain < 90) {
        throw "not enough new hits?";
    }
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
    vector <dd_edge> E1(EDGES);
    vector <dd_edge> E2(EDGES);

    initEdges(F1, E1);
    initEdges(F2, E2);

    vector <node_handle> N1(EDGES);
    vector <node_handle> N2(EDGES);

    for (unsigned i=0; i<EDGES; i++) {
        N1[i] = E1[i].getNode();
        N2[i] = E2[i].getNode();
    }

    //
    // Add some CT entries
    //

    myop* foo = new myop(F1, F2);

    cout << "    Creating CT entries\n";
    const unsigned add1 = addEntries(foo->et0, N1);
    cout << "        " << setw(5) << add1 << " op1 entries added\n";

    const unsigned add2 = addEntries(foo->et1, N1, N1);
    cout << "        " << setw(5) << add2 << " op2 entries added\n";

    const unsigned add3 = addEntries(foo->et2, N1, N2);
    cout << "        " << setw(5) << add3 << " op3 entries added\n";

    const unsigned entries1 = foo->et0->getNumEntries();
    const unsigned entries2 = foo->et1->getNumEntries();
    const unsigned entries3 = foo->et2->getNumEntries();

    show_entries(entries1, add1, "op1");
    show_entries(entries2, add2, "op2");
    show_entries(entries3, add3, "op3");

#ifdef DEBUG_CT
    ostream_output out(cout);
    foo->ct2()->show(out, 5);
    out << "Forest 1:\n";
    F1->dump(out, SHOW_DETAILS);
    out << "Forest 2:\n";
    F2->dump(out, SHOW_DETAILS);
    out << foo->ct2()->getStats().numEntries << " entries\n";
#endif


    cout << "    Checking CT entries\n";

    const unsigned hits1 = checkEntries("op1", foo->et0, N1);
    const unsigned hits2 = checkEntries("op2", foo->et1, N1, N1);
    checkEntries("op3", foo->et2, N1, N2);

    //
    // Destroy forest 2 and re-count!
    //
    cout << "    Killing op3 forest (FID " << F2->FID() << ")\n";
    forest::destroy(F2);
#ifdef DEBUG_CT
    cout << "Updated CT\n";
    foo->ct2()->show(out, 5);
    out << foo->ct2()->getStats().numEntries << " entries\n";
#endif
    cout << "    Re-checking CT entries\n";

    const unsigned nhits1 = checkEntries("op1", foo->et0, N1);
    const unsigned nhits2 = checkEntries("op2", foo->et1, N1, N1);
    checkEntries("op3", foo->et2, N1, N2);

    check_new_hits(hits1, nhits1);
    check_new_hits(hits2, nhits2);

#ifdef DEBUG_CT
    cout << "All CT entry types\n";
    ct_entry_type::showAll(out);
#endif

}

bool setStyleNumber(unsigned i)
{
    switch (i) {
        case 0:
            std::cout << "Monolithic chained ";
            ct_initializer::setBuiltinStyle(
                    ct_initializer::MonolithicChainedHash);
            return true;

        case 1:
            std::cout << "Monolithic unchained ";
            ct_initializer::setBuiltinStyle(
                    ct_initializer::MonolithicUnchainedHash);
            return true;

        case 2:
            std::cout << "Operation chained ";
            ct_initializer::setBuiltinStyle(
                    ct_initializer::OperationChainedHash);
            return false;

        case 3:
            std::cout << "Operation unchained ";
            ct_initializer::setBuiltinStyle(
                    ct_initializer::OperationUnchainedHash);
            return false;

        default:
            throw "Bad style number";
    }
}

void setCompressionNumber(unsigned i)
{
    switch (i) {
        case 0:
            std::cout << "(uncompressed)";
            ct_initializer::setCompression(compressionOption::None);
            return;

        case 1:
            std::cout << "(type compressed)";
            ct_initializer::setCompression(compressionOption::TypeBased);
            return;

        default:
            throw "Bad compression number";

    }
}

int main(int argc, const char** argv)
{
    using namespace std;
    try {
        bool mono;
        // unsigned st=2;
        // unsigned cmp=0;
        for (unsigned st=0; st<4; st++) {
            for (unsigned cmp=0; cmp<2; cmp++) {
                MEDDLY::initializer_list* IL = defaultInitializerList(nullptr);
                cout << "=================================================================\n";
                cout << "    Testing ";
                mono = setStyleNumber(st);
                setCompressionNumber(cmp);
                cout << " style CTs\n";
                cout << "=================================================================\n";
                MEDDLY::initialize(IL);
                check_CT(mono);
                MEDDLY::cleanup();

            } // compression
        } // style

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

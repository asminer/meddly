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
//
void addEntries(const ct_entry_type* CTE, compute_table* CT,
        const std::vector <dd_edge> &E)
{
    ct_entry_result res;
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
        } // for j
    } // for i

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
void addEntries(const ct_entry_type* CTE, compute_table* CT,
        const std::vector <dd_edge> &E1, const std::vector <dd_edge> &E2)
{
    ct_entry_result res;
    for (unsigned i=0; i<E1.size(); i++) {
        for (unsigned j=0; j<E2.size(); j++) {
            ct_entry_key* key = compute_table::useEntryKey(CTE, 0);

            key->writeN(E1[i].getNode());
            key->writeN(E2[j].getNode());

            CT->find(key, res);
            if (res) throw "Found result in CT; shouldn't have.";

            res.reset();
            res.writeI(E1[i].getNode() + E2[i].getNode());
            CT->addEntry(key, res);
        } // for j
    } // for i

}


int main(int argc, const char** argv)
{
    try {
        std::vector <dd_edge> E(100);

        int bounds[3];
        bounds[0] = bounds[1] = bounds[2] = VARSIZE;

        MEDDLY::initialize();

        domain *D = domain::createBottomUp(bounds, 3);
        forest* F;
        policies p;
        p.useDefaults(false);
        F = forest::create(D, false, range_type::INTEGER,
                        edge_labeling::MULTI_TERMINAL, p);

        initEdges(F, E);

        ct_entry_type* et = new ct_entry_type("test1", "IN:I");
        et->setForestForSlot(1, F);
        compute_table::registerEntryType(first_etid, et);

        compute_table CT = ...



        std::cout << "Done\n";

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
    return 4;
}

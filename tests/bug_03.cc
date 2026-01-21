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
    Tests union, intersection, set difference operations.
*/


#include <cstdlib>
#include <cassert>

#include "../src/meddly.h"

const unsigned VARS = 5;
const unsigned RELDOM = 2;
const unsigned SETBITS = 2;
const unsigned POTENTIAL = 1024;    // SETDOM ^ VARS

using namespace MEDDLY;

void showSet(std::ostream &out, const std::vector <bool> &elems)
{
    out << "{ ";
    bool printed = false;
    for (unsigned i=0; i<elems.size(); i++) {
        if (!elems[i]) continue;
        if (printed) out << ", ";
        out << i;
        printed = true;
    }
    out << " }";
}

void showRelMinterms(std::ostream &out, const std::vector <bool> &elems)
{
    out << "{ ";
    bool printed = false;
    for (unsigned i=0; i<elems.size(); i++) {
        if (!elems[i]) continue;
        if (printed) out << ",\n            ";
        out << "[";
        printed = true;

        unsigned x = i;
        for (unsigned j=1; j<=VARS; j++) {
            unsigned un = x % RELDOM;
            x /= RELDOM;
            unsigned pr = x % RELDOM;
            x /= RELDOM;
            if (j>1) out << ", ";
            out << un << "->" << pr;
        }
        out << "]";
    }
    out << " }";
}


inline void fillMinterm(unsigned x, int* un, int* pr)
{
#ifdef DEBUG_MINTERMS
    std::cout << std::hex << x << ": [";
#endif
    for (unsigned j=1; j<=VARS; j++) {
        un[j] = x % RELDOM;
        x /= RELDOM;
        pr[j] = x % RELDOM;
        x /= RELDOM;
    }
#ifdef DEBUG_MINTERMS
    for (unsigned j=VARS; j; j--) {
        if (j<VARS) std::cout << ", ";
        std::cout << un[j] << "->" << pr[j];
    }
    std::cout << "]\n";
#endif
}

inline unsigned whichMinterm(const int* un, const int* pr)
{
    unsigned x=0;
    for (unsigned j=VARS; j; j--) {
        x <<= SETBITS;
        x |= (un[j] << (SETBITS/2)) | pr[j];
    }
    return x;
}


void set2mxd(const std::vector<bool> &S, forest *F, dd_edge &s)
{
    if (!F) throw "null forest";

    // Determine S cardinality
    unsigned card = 0;
    for (unsigned i=0; i<S.size(); i++) {
        if (S[i]) ++card;
    }

    // Special case - empty set
    if (0==card) {
        F->createEdge(false, s);
        return;
    }

    // Convert set S to list of minterms
    int** unlist = new int* [card];
    int** prlist = new int* [card];
    card = 0;
    for (unsigned i=0; i<S.size(); i++) {
        if (!S[i]) continue;
        unlist[card] = new int[VARS+1];
        prlist[card] = new int[VARS+1];
        fillMinterm(i, unlist[card], prlist[card]);
        ++card;
    }

    F->createEdge(unlist, prlist, card, s);

    // Cleanup
    for (unsigned i=0; i<card; i++) {
        delete[] unlist[i];
        delete[] prlist[i];
    }
    delete[] unlist;
    delete[] prlist;
}


int main()
{
    using namespace MEDDLY;

    try {
        MEDDLY::initialize();
        policies p;

        int br[VARS];
        for (unsigned i=0; i<VARS; i++) {
            br[i] = RELDOM;
        }
        domain* Dr = domain::createBottomUp(br, VARS);

        p.useDefaults(RELATION);

        p.setFullyReduced();

        forest* F = forest::create(Dr, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);

        std::vector <bool> Aset(POTENTIAL, false);

        Aset[512] = 1;
        Aset[513] = 1;
        Aset[514] = 1;
        Aset[515] = 1;

        std::cout << "Relation: ";
        showSet(std::cout, Aset);
        std::cout << "\n";

        std::cout << "Minterms: ";
        showRelMinterms(std::cout, Aset);
        std::cout << "\n";

        dd_edge Add(F);

        std::cout << "Mxd: ";
        set2mxd(Aset, F, Add);
        ostream_output out(std::cout);
        Add.showGraph(out);

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

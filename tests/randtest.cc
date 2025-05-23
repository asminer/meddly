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

//
// Test the randomization stuff in randomize.h/cc
//

#include "randomize.h"

const unsigned VARS = 5;
const unsigned RELDOM = 2;
const unsigned SETDOM = RELDOM * RELDOM;

//
// OLD IMPLEMENTATION HERE, for comparison
//

const unsigned SETBITS = 2;

inline void fillMinterm(unsigned x, int* mt)
{
    for (unsigned j=1; j<=VARS; j++) {
        mt[j] = x % SETDOM;
        x /= SETDOM;
    }
/*
    std::cout << "got digits " << mt[VARS];
    for (unsigned j=VARS-1; j; --j) {
        std::cout << ", " << mt[j];
    }
    std::cout << "\n";
    */
}

inline void fillMinterm(unsigned x, int* un, int* pr)
{
    for (unsigned j=1; j<=VARS; j++) {
        pr[j] = x % RELDOM;
        x /= RELDOM;
        un[j] = x % RELDOM;
        x /= RELDOM;
    }
/*
    std::cout << "got digits ";
    for (unsigned j=VARS; j; --j) {
        if (j < VARS) std::cout << ", ";
        std::cout << un[j] << " -> " << pr[j];
    }
    std::cout << "\n";
*/
}

inline unsigned whichMinterm(const int* mt)
{
    unsigned x=0;
    for (unsigned j=VARS; j; j--) {
        x <<= SETBITS;
        x |= mt[j];
    }
/*
    std::cout << "digits " << mt[VARS];
    for (unsigned j=VARS-1; j; --j) {
        std::cout << ", " << mt[j];
    }
    std::cout << " --> " << x << "\n";
*/
    return x;
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

void flipFullyElements(std::vector <bool> &elems, unsigned seed, unsigned k,
        bool on)
{
    // Convert 'seed' to minterm
    int tmp[VARS+1];
    fillMinterm(seed, tmp);
    // Replace variable k with all possible values in minterm,
    // and convert back to unsigned
    for (tmp[k]=0; tmp[k]<SETDOM; tmp[k]++) {
        unsigned x = whichMinterm(tmp);
        elems[x] = on;
    }
}

void flipFullyElements(std::vector <bool> &elems, unsigned seed, unsigned k1,
        unsigned k2, bool on)
{
    // Convert 'seed' to minterm
    int tmp[VARS+1];
    fillMinterm(seed, tmp);
    // Replace variable k with all possible values in minterm,
    // and convert back to unsigned
    for (tmp[k1]=0; tmp[k1]<SETDOM; tmp[k1]++) {
        for (tmp[k2]=0; tmp[k2]<SETDOM; tmp[k2]++) {
            unsigned x = whichMinterm(tmp);
            elems[x] = on;
        }
    }
}


void flipIdentityElements(std::vector <bool> &elems, unsigned seed,
        unsigned k, bool on)
{
    // Convert 'seed' to minterm
    int untmp[VARS+1];
    int prtmp[VARS+1];
    fillMinterm(seed, untmp, prtmp);
    for (unsigned v=0; v<RELDOM; v++) {
        untmp[k] = v;
        prtmp[k] = v;
        unsigned x = whichMinterm(untmp, prtmp);
        elems[x] = on;
    }
}


void flipIdentityElements(std::vector <bool> &elems, unsigned seed,
        unsigned k1, unsigned k2, bool on)
{
    // Convert 'seed' to minterm
    int untmp[VARS+1];
    int prtmp[VARS+1];
    fillMinterm(seed, untmp, prtmp);
    for (unsigned v1=0; v1<RELDOM; v1++) {
        untmp[k1] = v1;
        prtmp[k1] = v1;
        for (unsigned v2=0; v2<RELDOM; v2++) {
            untmp[k2] = v2;
            prtmp[k2] = v2;
            unsigned x = whichMinterm(untmp, prtmp);
            elems[x] = on;
        }
    }
}


bool equals(const std::vector <bool> &set1, const std::vector <bool> &set2)
{
    if (set1.size() != set2.size()) return false;
    for (unsigned i=0; i<set1.size(); i++) {
        if (set1[i] != set2[i]) return false;
    }
    return true;
}

void showset(std::ostream &out, const std::vector <bool> &elems)
{
    out << "{";
    bool printed = false;
    for (unsigned i=0; i<elems.size(); i++) {
        if (!elems[i]) continue;
        if (printed) {
            out << ", ";
        } else {
            printed = true;
        }
        out << i;
    }
    out << "}";
}

inline void clrset(std::vector <bool> &set)
{
    for (unsigned i=0; i<set.size(); i++) {
        set[i] = false;
    }
}

void list2set(const std::vector <unsigned> &list, std::vector <bool> &set)
{
    clrset(set);
    for (unsigned i=0; i<list.size(); i++) {
        set[ list[i] ] = true;
    }
}

void show(std::ostream &out, const std::vector <unsigned> &list)
{
    out << "[";
    for (unsigned i=0; i<list.size(); i++) {
        if (i) out << ", ";
        out << list[i];
    }
    out << "]";
}

int main(int argc, const char** argv)
{
    using namespace std;

    //
    // First argument: seed
    //
    long seed = 0;
    if (argv[1]) {
        seed = atol(argv[1]);
    }

    try {
        MEDDLY::initialize();
        vectorgen::setSeed(seed);

        //
        // Test minterm indexing for sets
        //
        vectorgen Gs(MEDDLY::SET, VARS, SETDOM);
        MEDDLY::domain* Ds = Gs.makeDomain();
        MEDDLY::minterm mts(Ds, MEDDLY::SET);
        cout << "Testing set indexing\n    ";
        for (unsigned i=0; i<Gs.potential(); i++) {
            Gs.index2minterm(i, mts);
            unsigned _i = Gs.minterm2index(mts);
            if (i != _i) {
                cout << "\n";
                cout << "Index " << i << "\nproduced minterm ";
                MEDDLY::ostream_output out(cout);
                mts.show(out);
                cout << "\nand got index " << _i << "\n";
                throw "index mismatch";
            }
            cout << ".";
        }
        cout << "\n";


        //
        // Test minterm indexing for relations
        //
        vectorgen Gr(MEDDLY::RELATION, VARS, RELDOM);
        MEDDLY::domain* Dr = Gr.makeDomain();
        MEDDLY::minterm mtr(Dr, MEDDLY::RELATION);
        cout << "Testing relation indexing\n    ";
        for (unsigned i=0; i<Gr.potential(); i++) {
            Gr.index2minterm(i, mtr);
            unsigned _i = Gr.minterm2index(mtr);
            if (i != _i) {
                cout << "\n";
                cout << "Index " << i << "\nproduced minterm ";
                MEDDLY::ostream_output out(cout);
                mts.show(out);
                cout << "\nand got index " << _i << "\n";
                throw "index mismatch";
            }
            cout << ".";
        }
        cout << "\n";


        //
        // Test buildFullyFromIndex
        //

        vector <unsigned> L;

        cout << "Testing buildFullyFromIndex\n    ";

        vector <bool> gset;
        vector <bool> oldset;

        gset.resize(Gr.potential());
        oldset.resize(Gr.potential());

        // unsigned x = 37;

        for (unsigned x=0; x<Gr.potential(); x++) {
            cout << ".";

            for (unsigned k1=1; k1 <= VARS; k1++) {
                for (unsigned k2=k1; k2 <= VARS; k2++) {
                    Gr.buildIdentityFromIndex(x, k1, k2, L);
                    list2set(L, gset);

                    clrset(oldset);
                    if (k1 == k2) {
                        flipIdentityElements(oldset, x, k1, true);
                    } else {
                        flipIdentityElements(oldset, x, k1, k2, true);
                    }

                    if (!equals(gset, oldset)) {
                        cerr << "Mismatch: x=" << x << ", k1 = " << k1
                            << ", k2 = " << k2 << "\n";

                        cerr << "gset: ";
                        showset(cerr, gset);
                        cerr << "\noldset: ";
                        showset(cerr, oldset);
                        cerr << "\nL: ";
                        show(cerr, L);
                        cerr << "\n";
                        throw "mismatch";
                    }
                }
            }

        } // for x
        cout << "\n";

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
        cerr << "\nCaught our own error: " << e << "\n";
        return 2;
    }
}



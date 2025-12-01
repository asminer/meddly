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

#include <iostream>
#include <vector>

#include "../src/meddly.h"

// Number of rings.
int N;

// **********************************************************************
// build an event to move a ring
//    @param ring2level     Mapping from rings to level
//    @param ring           Which ring to move
//    @param src            Source peg (0, 1, 2)
//    @param dest           Destination peg (0, 1, 2)
//    @param aux            Unused peg (0, 1, 2)
//
//    @param E              in:  blank dd_edge attached to relation forest
//                          out: dd_edge holding the event
// **********************************************************************
void moveRing(std::vector <int> &ring2level, int ring,
        int src, int dest, int aux, MEDDLY::dd_edge &E)
{
    assert(src != dest);
    assert(src != aux);
    assert(dest != aux);

    MEDDLY::minterm delta(E.getForest());

    for (unsigned i=N; i>0; i--) {
        const int k = ring2level[i];
        if (i > ring) {
            delta.setVars(k, MEDDLY::DONT_CARE, MEDDLY::DONT_CHANGE);
            continue;
        }
        if (i < ring) {
            delta.setVars(k, aux, aux);
            continue;
        }
        delta.setVars(k, src, dest);
    }

    delta.buildFunction(0, E);
}

// **********************************************************************
// print usage instructions
// **********************************************************************

int usage(const char* exe)
{
    using namespace std;
    const char* base = exe;
    for (; *exe; ++exe)
    {
        if ('/' == *exe)
        {
            base = exe+1;
        }
    }

    cerr << "\n";
    cerr << "Usage: " << base << " [switches] #rings\n";
    cerr << "\n";
    cerr << "Generate reachable states for the Towers of Hanoi puzzle.\n";
    cerr << "The number of rings given should be a natural number.\n";
    cerr << "\n";
    cerr << "Switches:\n";
    cerr << "    -h:    This help\n\n";
    cerr << "    -o:    Order the variables so that the smallest ring is the bottom-most\n";
    cerr << "           variable in the MDD. (Default)\n";
    cerr << "    -O:    Order the variables so that the largest ring is the bottom-most\n";
    cerr << "           variable in the MDD.\n";

    cerr << "\n";
    return 1;
}

inline int naturalArg(const char* arg)
{
    if (0==arg) return -1;
    return atol(arg);
}

// **********************************************************************
// main
// **********************************************************************

int main(int argc, const char** argv)
{
    using namespace std;

    N = -1;
    bool reverse_order = false;
    //
    // Process command line
    //
    for (unsigned i=1; i<argc; i++)
    {
        const char *arg = argv[i];
        if ('-' == arg[0])
        {
            // single character switches
            if (arg[1] && 0==arg[2]) {
                switch (arg[1])
                {
                    case 'o':
                                reverse_order = false;
                                continue;

                    case 'O':
                                reverse_order = true;
                                continue;

                    default:
                                if (arg[1] != 'h') {
                                    cerr << "Unknown switch \"-" << arg[1] << "\"\n";
                                }
                                return usage(argv[0]);
                }
            }
            cerr << "Unknown switch \"" << arg << "\"\n";
            return usage(argv[0]);
        }
        if (N>=0) {
            cerr << "More than one #rings specified?\n";
            return usage(argv[0]);
        }
        N = naturalArg(arg);
        if (N<0) {
            return usage(argv[0]);
        }
    }
    if (N<0) {
        cerr << "#rings not specified\n";
        return usage(argv[0]);
    }

    //
    // Build ring to level mappings
    //
    cout << "Towers of Hanoi with " << N << " rings.\n";
    if (0==N) return 0;

    vector <int> ring2level(N+1);
    vector <int> level2ring(N+1);
    ring2level[0] = 0;
    level2ring[0] = 0;
    if (reverse_order) {
        for (unsigned i=N; i>0; i--) {
            ring2level[i] = N-i+1;
            level2ring[N-i+1] = i;
        }
    } else {
        for (unsigned i=N; i>0; i--) {
            ring2level[i] =i;
            level2ring[i] = i;
        }
    }

    cout << "Variable order is TOP, ";
    for (unsigned i=N; i>0; i--) {
        cout << "ring " << level2ring[i] << ", ";
        if (i>3) {
            cout << "..., ";
            i = 3;
        }
    }
    cout << "BOTTOM\n";

    return 0;
}

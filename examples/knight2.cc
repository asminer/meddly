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
#include <iomanip>
#include <cstdlib>

#include "../src/meddly.h"

// using namespace MEDDLY;

// Board dimensions

#define MAXDIM 10

bool verbose, debug;
int N, M;
int start_n, start_m;


// Forest for multi-terminal MDDs
MEDDLY::forest* mtF;
// Forest for boolean constraints
MEDDLY::forest* boolF;

void buildConstraints(const coord &sq, constraint* C)
{
    if (0==C) return;

    using namespace std;
    using namespace MEDDLY;

    if (verbose) {
        cout << "--------------------------------------\n";
        cout << "Constraints for ";
        showVarOf(cout, sq) << "\n";
    }

    C->vardeps[B.get_var(sq)] = true;

    unsigned i, j;
    coord_list moves(8);
    MEDDLY::ostream_output sout(cout);

    knight_moves(sq, moves);
    if (verbose) {
        cout << "    Moves:\n";
        for (i=0; i<moves.length(); i++) {
            cout << "\t-> ";
            showVarOf(cout, moves.get(i)) << " = " << moves.get(i) << "\n";
        }
    }
    if (verbose) {
        cout << "    Disjunction of:\n";
        cout << "\t(";
        showVarOf(cout, sq);
        cout << "==" << B.getN()*B.getM() << ")\n";
    }
    dd_edge xsq(mtF);
    mtF->createEdgeForVar(B.get_var(sq), false, xsq);
    dd_edge one(mtF);
    mtF->createEdge(long(1), one);
    xsq += one;

    // xsq holds 1+ visit number for the given square

    dd_edge NxM(mtF);
    mtF->createEdge(long(B.getN()*B.getM()), NxM);
    apply(EQUAL, xsq, NxM, C->the_dd);


    for (i=0; i<moves.length(); i++) {
        dd_edge ands(boolF);
        boolF->createEdge(true, ands);
        if (verbose) cout << "\t";
        for (j=0; j<moves.length(); j++) {
            if (verbose) {
                if (j) cout << " ^ ";
            }
            if (verbose || debug) {
                cout << "(1+";
                showVarOf(cout, sq);
                if (i == j) {
                    cout << "==";
                } else {
                    cout << "!=";
                }
                showVarOf(cout, moves.get(j)) << ")";
            }
            // build constraint here...
            //

            dd_edge xmv(mtF);
            int vm = B.get_var(moves.get(j));
            mtF->createEdgeForVar(vm, false, xmv);
            C->vardeps[vm] = true;

            dd_edge term(boolF);
            if (i == j) {
                apply(EQUAL, xsq, xmv, term);
            } else {
                apply(NOT_EQUAL, xsq, xmv, term);
            }

            if (debug) {
                cout << "\n";
                term.show(sout, 2);
            }
            ands *= term;
        }
        if (verbose) cout << "\n";
        if (debug) {
            cout << "\nANDS:\n\n";
            ands.show(sout, 2);
        }
        // OR with the others here
        C->the_dd += ands;
    }
    if (debug) {
        cout << "\nEntire constraint:\n\n";
        C->the_dd.show(sout, 2);
    }
}

void show_solution(std::ostream &s, const int* minterm)
{
    using namespace std;
    coord c;
    for (c.n=1; c.n<=B.getN(); c.n++) {
        for (c.m=1; c.m<=B.getM(); c.m++) {
            cout << setfill(' ') << setw(3) << minterm[B.get_var(c)] << ' ';
        }
        cout << '\n';
    }
}

void usage(const char* arg)
{
    using namespace std;
    cout << "\nUsage: " << arg << " [switches]\n\n";
    cout << "Build and count solutions to the Knight's tour on an NxM chessboard.\n";
    cout << "\nLegal switches:\n";
    cout << "\t-h:\tThis help\n\n";
    cout << "\t-b n m\tBeginning position (default is 1 1)\n";
    cout << "\t-d:\tToggle debugging (default: off)\n";
    cout << "\t-n N\tSet number of rows\n";
    cout << "\t-m M\tSet number of columns\n";
    cout << "\t-r:\tReverse variable order\n";
    cout << "\t-v:\tToggle verbosity (default: off)\n";
    cout << "\n";
    exit(1);
}

inline long intFrom(const char* arg)
{
    if (0==arg) return 0;
    return atol(arg);
}

void process_args(int argc, const char** argv)
{
    // Set defaults for everything
    N = 4;
    M = 4;
    start_n = 1;
    start_m = 1;
    reverse = false;
    verbose = false;
    debug = false;
    long L, L2;

    // Go through args
    unsigned i;
    for (i=1; i<argc; i++) {
        char sw = ' ';
        const char* arg = argv[i];
        if (('-' == arg[0]) || (0 != arg[1]) || (0 == arg[2])) {
            sw = arg[1];
        }

        switch (sw) {
            case 'h':   usage(argv[0]);
                        continue;

            case 'b':   if (i+2<argc) {
                            L = intFrom(argv[++i]);
                            L2 = intFrom(argv[++i]);
                            if (L<1 || L>MAXDIM || L2<1 || L2>MAXDIM) {
                                std::cerr << "Value for -b out of range; ignoring\n";
                            } else {
                                start_n = L;
                                start_m = L2;
                            }
                        } else {
                            std::cerr << "Switch -b requires two arguments\n";
                            exit(2);
                        }
                        continue;

            case 'd':   debug = !debug;
                        continue;

            case 'n':   if (i+1<argc) {
                            L = intFrom(argv[++i]);
                            if (L<1 || L>MAXDIM) {
                                std::cerr << "Value " << L << " for -n out of range; ignoring\n";
                            } else {
                                N = L;
                            }
                        } else {
                            std::cerr << "Switch -n requires an argument\n";
                            exit(2);
                        }
                        continue;

            case 'm':   if (i+1<argc) {
                            L = intFrom(argv[++i]);
                            if (L<1 || L>MAXDIM) {
                                std::cerr << "Value " << L << " for -m out of range; ignoring\n";
                            } else {
                                M = L;
                            }
                        } else {
                            std::cerr << "Switch -m requires an argument\n";
                            exit(2);
                        }
                        continue;

            case 'r':   reverse = !reverse;
                        continue;

            case 'v':   verbose = !verbose;
                        continue;

            default:    std::cerr << "\nUnknown option: " << arg << "\n";
                        std::cerr << "Run with -h for help\n\n";
                        exit(2);
        }
    }

    if (start_n<1 || start_n>N || start_m<1 || start_m>M) {
        std::cerr << "Starting position " << start << " is not on the board\n";
        exit(2);
    }

}

int main(int argc, const char** argv)
{
    using namespace MEDDLY;

    process_args(argc, argv);

    const int numvars = N * M * 2;
    int* bounds = new int[numvars];
    for (int i=0; i<numvars; i++) {
        bounds[i] = (i%2) ? M : N;
    }

    //
    // Set up forests.
    //
    // Variables are, from the bottom up,
    //      col of first (or last if reversed) position
    //      row of first position
    //      col of second position
    //      row of second position
    //      ...
    //      col of last position
    //      row of last position
    //

    initialize();
    domain* D = createDomainBottomUp(bounds, numvars);

    policies p;
    p.useDefaults(false);
    p.useReferenceCounts = false;
    // p.setPessimistic();

    mtF = D->createForest(false, range_type::INTEGER,
                    edge_labeling::MULTI_TERMINAL, p);

    boolF = D->createForest(false, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    //
    // TBD HERE
    //

    constraint** clist = new constraint* [B.getN()*B.getM()];

    i = NxM;
    while (i) {
        --i;
        clist[i] = new constraint(boolF);
    }

    std::cout << "Building constraints per square" << std::endl;

    /*
     * Build constraints for each square, determining
     * what must come after
     */
    coord c;
    for (c.n=1; c.n<=B.getN(); c.n++)
        for (c.m=1; c.m<=B.getM(); c.m++) {
            buildConstraints(c, clist[(c.n-1)*B.getM()+c.m-1]);
        }

    /*
     * Build initial position constraint
     *
     */
    dd_edge all(boolF);
    dd_edge xst(mtF);
    dd_edge zero(mtF);

    mtF->createEdgeForVar(B.get_var(start), false, xst);
    mtF->createEdge(long(0), zero);
    apply(EQUAL, xst, zero, all);

    std::cout << "Sorting constraints" << std::endl;
    i = NxM;
    while (i) {
        i--;
        // Put 'smallest' at end
        unsigned min = i;
        for (unsigned j=0; j<i; j++) {
            if (less_than(clist[j], clist[min])) {
                min = j;
            }
        }
        if (min != i) {
            SWAP(clist[i], clist[min]);
        }
    }

    std::cout << "Anding constraints" << std::endl;
    i=NxM;
    try {
        while (i) {
            std::cout << "    " << std::setw(2) << i << " : ";
            i--;
            show_deps(clist[i], NxM);
            std::cout.flush();

            all *= clist[i]->the_dd;
            delete clist[i];
            clist[i] = nullptr;

            std::cout << "  (" << boolF->getCurrentNumNodes() << " nodes)"
                << std::endl;
        }

        long allcard;
        apply(CARDINALITY, all, allcard);
        std::cout << allcard << " tours total\n";

        // Show solutions
        enumerator sol(all);
        for (; sol; ++sol) {
            std::cout << "Solution:\n";
            show_solution(std::cout, sol.getAssignments());
        }

        // Memory stats
        std::cout << "Forest stats:\n";
        ostream_output myout(std::cout);
        expert_forest* ef = (expert_forest*)boolF;
        ef->reportStats(myout, "\t",
            expert_forest::HUMAN_READABLE_MEMORY  |
            expert_forest::BASIC_STATS | expert_forest::EXTRA_STATS |
            expert_forest::STORAGE_STATS | expert_forest::STORAGE_DETAILED |
            expert_forest::HOLE_MANAGER_STATS | expert_forest::HOLE_MANAGER_DETAILED
        );
    }
    catch (MEDDLY::error e) {
        std::cerr << "\n\nCaught MEDDLY error: " << e.getName() << "\n";
        std::cerr << "    thrown at " << e.getFile() << " line " << e.getLine() << "\n";
    }



    MEDDLY::cleanup();
    return 0;
}

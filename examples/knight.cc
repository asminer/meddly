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
#define SHOWSOLS 10

bool topdown, reverse, verbose, debug;
int N, M;
long maxsols;
long start_n, start_m;


// Forest for multi-terminal MDDs
MEDDLY::forest* mtF;
// Forest for boolean constraints
MEDDLY::forest* boolF;

const char* with_sign(long delta)
{
    static char buffer[3];
    buffer[2] = 0;
    if (delta<0) {
        buffer[0] = '-';
        buffer[1] = '0' - delta;
    } else {
        buffer[0] = '+';
        buffer[1] = '0' + delta;
    }
    return buffer;
}

void setStart(int moveno, MEDDLY::dd_edge& constr)
{
    using namespace std;
    using namespace MEDDLY;

    if (verbose) {
        cout << "--------------------------------------\n";
        cout << "Constraint for starting position\n";
    }
    const int thisrow = moveno*2;
    const int thiscol = thisrow-1;
    dd_edge xtrow(mtF), xtcol(mtF), strow(mtF), stcol(mtF);

    if (verbose) {
        cout << "    (x_" << thiscol << " == " << start_m << ") ^ ";
        cout << "(x_" << thisrow << " == " << start_n << ")\n";
    }

    mtF->createEdgeForVar(thisrow, false, xtrow);
    mtF->createEdgeForVar(thiscol, false, xtcol);

    mtF->createEdge(start_n, strow);
    mtF->createEdge(start_m, stcol);

    dd_edge term2(boolF);

    apply(EQUAL, xtrow, strow, constr);
    apply(EQUAL, xtcol, stcol, term2);
    constr *= term2;

    if (debug) {
        cout << "Forest for start\n";
        ostream_output myout(cout);
        constr.showGraph(myout);
    }

}

void buildConstraint(int movea, int moveb, MEDDLY::dd_edge& constr)
{
    using namespace std;
    using namespace MEDDLY;

    if (verbose) {
        cout << "--------------------------------------\n";
        cout << "Constraints for move from position " << movea << "->" << moveb << "\n";
        cout << "Disjunction of:\n";
    }

    const int thisrow = movea*2;
    const int thiscol = thisrow-1;
    const int nextrow = moveb*2;
    const int nextcol = nextrow-1;

    const long drow[] = {  1,  1, -1, -1,  2,  2, -2, -2, 0 };
    const long dcol[] = {  2, -2,  2, -2,  1, -1,  1, -1, 0 };

    boolF->createEdge(false, constr);

    for (unsigned i=0; drow[i]; i++) {
        if (verbose) {
            cout << "    (x_" << thiscol << with_sign(dcol[i]) << " == x_" << nextcol << ") ^ ";
            cout << "(x_" << thisrow << with_sign(drow[i]) << " == x_" << nextrow << ")\n";
        }

        dd_edge move(boolF);
        dd_edge rowch(boolF), colch(boolF);

        dd_edge xtrow(mtF), xtcol(mtF), xnrow(mtF), xncol(mtF);
        dd_edge dr(mtF), dc(mtF);

        mtF->createEdgeForVar(thisrow, false, xtrow);
        mtF->createEdgeForVar(thiscol, false, xtcol);
        mtF->createEdgeForVar(nextrow, false, xnrow);
        mtF->createEdgeForVar(nextcol, false, xncol);

        mtF->createEdge(drow[i], dr);
        mtF->createEdge(dcol[i], dc);

        xtrow += dr;
        xtcol += dc;

        apply(EQUAL, xtrow, xnrow, rowch);
        apply(EQUAL, xtcol, xncol, colch);
        apply(INTERSECTION, rowch, colch, move);

        constr += move;
    }

    if (debug) {
        cout << "Forest for move\n";
        ostream_output myout(cout);
        constr.showGraph(myout);
    }
}

void uniqueDown(int pnum, MEDDLY::dd_edge &constr)
{
    using namespace std;
    using namespace MEDDLY;

    if (verbose) {
        cout << "\n--------------------------------------\n";
        cout << "Position " << pnum << " different from all previous\n";
        cout << "Conjunction of:\n";
    }

    const int rowa = pnum*2;
    const int cola = rowa-1;

    dd_edge xrowa(mtF), xcola(mtF), xrowb(mtF), xcolb(mtF);
    dd_edge disj1(boolF), disj2(boolF);

    mtF->createEdgeForVar(rowa, false, xrowa);
    mtF->createEdgeForVar(cola, false, xcola);

    boolF->createEdge(true, constr);

    for (int i=pnum-1; i>0; i--) {
        const int rowb = i*2;
        const int colb = rowb-1;

        if (verbose) {
            cout << "    (x_" << rowa << " != x_" << rowb << ") v ";
            cout << "(x_" << cola << " != x_" << colb << ")\n";
        }

        mtF->createEdgeForVar(rowb, false, xrowb);
        mtF->createEdgeForVar(colb, false, xcolb);

        apply(NOT_EQUAL, xrowa, xrowb, disj1);
        apply(NOT_EQUAL, xcola, xcolb, disj2);
        disj1 += disj2;


        constr *= disj1;
    }

    if (debug) {
        cout << "Forest for constraint:\n";
        ostream_output myout(cout);
        constr.showGraph(myout);
    }
}

void usage(const char* arg)
{
    using namespace std;
    cout << "\nUsage: " << arg << " [switches]\n\n";
    cout << "Build and count solutions to the Knight's tour on an NxM chessboard.\n";
    cout << "\nLegal switches:\n";
    cout << "\t-h:\tThis help\n\n";
    cout << "\t-a t\tAccumulate/remove from top down\n";
    cout << "\t-a b\tAccumulate/remove from bottom up (default)\n";
    cout << "\t-b n m\tBeginning position (default is 0 0)\n";
    cout << "\t-d:\tToggle debugging (default: off)\n";
    cout << "\t-n N\tSet number of rows\n";
    cout << "\t-m M\tSet number of columns\n";
    cout << "\t-r:\tReverse variable order\n";
    cout << "\t-s #\tMax number of solutions to display\n";
    cout << "\t-v:\tToggle verbosity (default: off)\n";
    cout << "\n";
    exit(1);
}

inline long intFrom(const char* arg)
{
    if (0==arg) return 0;
    return atol(arg);
}

inline char charFrom(const char* arg)
{
    if (0==arg) return 0;
    return arg[0];
}


void process_args(int argc, const char** argv)
{
    // Set defaults for everything
    N = 4;
    M = 4;
    start_n = 0;
    start_m = 0;
    reverse = false;
    verbose = false;
    debug = false;
    topdown = false;
    maxsols = 1;
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

            case 'a':   if (i+1<argc) {
                            sw = charFrom(argv[++i]);
                            if ('t' == sw) {
                                topdown = true;
                                continue;
                            }
                            if ('b' == sw) {
                                topdown = false;
                                continue;
                            }
                            std::cerr << "Invalid argument " << sw << " for -a\n";
                            exit(2);
                        } else {
                            std::cerr << "Switch -a requires an argument\n";
                            exit(2);
                        }
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
                            if (L<0 || L>=MAXDIM) {
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
                            if (L<0 || L>=MAXDIM) {
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

            case 's':   if (i+1<argc) {
                            L = intFrom(argv[++i]);
                            if (L<0) {
                                L=0;
                            }
                            maxsols = L;
                        } else {
                            std::cerr << "Switch -s requires an argument\n";
                            exit(2);
                        }
                        continue;

            case 'v':   verbose = !verbose;
                        continue;

            default:    std::cerr << "\nUnknown option: " << arg << "\n";
                        std::cerr << "Run with -h for help\n\n";
                        exit(2);
        }
    }

    if (start_n<0 || start_n>=N || start_m<0 || start_m>=M) {
        std::cerr << "Starting position " << start_n << "," << start_m
                  << " is not on the board\n";
        exit(2);
    }

}

void showSolution(std::ostream& s, const int* minterm, const char* sep=nullptr)
{
    using namespace std;
    bool first = true;
    for (int i=N*M; i; i--) {
        if (!first) {
            if (sep) s << sep;
            else s << endl;
        } else {
            first = false;
        }
        s << '(' << setfill(' ') << setw(2) << minterm[i*2];
        s << ',' << setfill(' ') << setw(2) << minterm[i*2-1] << ')';
    }
}

void showAtMost(std::ostream& s, const MEDDLY::dd_edge &sols, long count)
{
    MEDDLY::enumerator iter(sols);
    long c = 0;
    for (; iter; ++iter) {
        ++c;
        if (c>count) return;
        s << "#" << c << ":  ";
        showSolution(s, iter.getAssignments(), " -> ");
        s << std::endl;
    }
}

void generate()
{
    using namespace MEDDLY;

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

    domain* D = domain::createBottomUp(bounds, numvars);

    policies p;
    p.useDefaults(false);
    // p.useReferenceCounts = false;
    p.setPessimistic();

    mtF = forest::create(D, false, range_type::INTEGER,
                    edge_labeling::MULTI_TERMINAL, p);

    boolF = forest::create(D, false, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);


    dd_edge* moves = new dd_edge[N*M+1];
    for (int i=1; i<N*M; i++) {
        moves[i].attach(boolF);
        buildConstraint(i, i+1, moves[i]);
    }

    moves[0].attach(boolF);
    moves[N*M].attach(boolF);

    if (reverse) {
        boolF->createEdge(true, moves[0]);
        setStart(N*M, moves[N*M]);
    } else {
        boolF->createEdge(true, moves[N*M]);
        setStart(1, moves[0]);
    }

    dd_edge tours(boolF);
    boolF->createEdge(true, tours);

    std::cout << "\nDetermining all " << N*M-1 << " knight moves\n";
    if (topdown) {
        for (int i=N*M; i>=0; i--) {
            std::cout << "  " << i;
            std::cout.flush();
            tours *= moves[i];
        }
    } else {
        for (int i=0; i<=N*M; i++) {
            std::cout << "  " << i;
            std::cout.flush();
            tours *= moves[i];
        }
    }
    delete[] moves;
    moves = nullptr;
    std::cout << "\n";

    double card;
    apply(CARDINALITY, tours, card);
    std::cout << "Approx. " << card << " total moves\n";

    std::cout << "\nRestricting to tours\n";
    std::cout << "\tBuilding unique constraint\n";

    dd_edge unique(boolF), nodup(boolF);
    boolF->createEdge(true, unique);
    if (topdown) {
        for (int i=N*M; i>0; i--) {
            std::cout << "    " << i;
            std::cout.flush();
            uniqueDown(i, nodup);
            std::cout << "..1";
            std::cout.flush();
            unique *= nodup;

            std::cout << " \t(" << boolF->getCurrentNumNodes()
                      << " forest nodes)"
                      << std::endl;
        }
    } else {
        for (int i=1; i<=N*M; i++) {
            std::cout << "    " << i;
            std::cout.flush();
            uniqueDown(i, nodup);
            std::cout << "..1";
            std::cout.flush();
            unique *= nodup;

            std::cout << " \t(" << boolF->getCurrentNumNodes()
                      << " forest nodes)"
                      << std::endl;
        }
    }

    std::cout << "\tApplying unique constraint\n";
    tours *= unique;

    apply(CARDINALITY, tours, card);
    std::cout << "Approx. " << card << " tours (no repeated squares)\n";

    dd_edge hamil(boolF);
    buildConstraint(1, N*M, hamil);
    hamil *= tours;

    apply(CARDINALITY, hamil, card);
    std::cout << "Approx. " << card << " Hamiltonian tours\n";

    dd_edge unsat(boolF);
    boolF->createEdge(false, unsat);

    if (maxsols) {
        //
        // Display some solutions
        //

        if (hamil != unsat) {
            std::cout << "First " << maxsols << " Hamiltonian tours:\n";
            showAtMost(std::cout, hamil, maxsols);
        } else if (tours != unsat) {
            std::cout << "First " << maxsols << " tours:\n";
            showAtMost(std::cout, tours, maxsols);
        }
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

//    boolF->showInfo(myout, 1);
}

int main(int argc, const char** argv)
{
    process_args(argc, argv);

    try {
        MEDDLY::initialize();
        generate();
        MEDDLY::cleanup();
        return 0;
    }
    catch (MEDDLY::error e) {
        std::cerr << "\n\nCaught MEDDLY error: " << e.getName() << "\n";
        std::cerr << "    thrown at " << e.getFile() << " line " << e.getLine() << "\n";
        return 1;
    }

    /*
    std::cout << "Anding constraints" << std::endl;
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

    */
}

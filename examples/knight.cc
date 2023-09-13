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


struct coord {
    signed char n, m;
    inline void set(int _n, int _m) {
        n = _n;
        m = _m;
    }
};

std::ostream& operator<< (std::ostream &s, coord c)
{
    return s << '(' << int(c.n) << ", " << int(c.m) << ')';
}

class board {
        int N, M;

        unsigned* var_at;

    public:
        board();
        ~board();

        void init(int n, int m);

        inline int getN() const { return N; }
        inline int getM() const { return M; }

        inline bool in_bounds(int n, int m) const {
            if (n < 1) return false;
            if (n > N) return false;
            if (m < 1) return false;
            if (m > M) return false;
            return true;
        }

        inline unsigned get_var(const coord &c) const {
            return var_at[(c.n-1)*M+c.m-1];
        }
        inline void set_var(const coord &c, unsigned v) {
            var_at[(c.n-1)*M+c.m-1] = v;
        }
        void reverse_vars();
};

board::board()
{
    var_at = nullptr;
    N = 0;
    M = 0;
}

board::~board()
{
    delete[] var_at;
}

void board::init(int n, int m)
{
    delete[] var_at;
    N = n;
    M = m;
    var_at = new unsigned[N*M];
    for (int i=0; i<N*M; i++) {
        var_at[i] = 0;
    }
}

void board::reverse_vars()
{
    unsigned adj = unsigned(N*M+1);
    for (int i=0; i<N*M; i++) {
        var_at[i] = adj - var_at[i];
    }
}

board B;
bool verbose, debug;

class coord_list {
        coord* coords;
        unsigned coords_size;
        unsigned coords_len;
    public:
        coord_list(unsigned sz);
        ~coord_list();

        inline void clear() {
            coords_len = 0;
        }

        inline void append(int n, int m) {
            if (!B.in_bounds(n, m)) return;
            if (coords_len >= coords_size) throw "ugh";
            coords[coords_len++].set(n, m);
        }
        inline unsigned length() const { return coords_len; }

        inline const coord& get(unsigned i) const { return coords[i]; }

        inline void reach(const coord &c) {
            if (coords_len >= coords_size) throw "ugh";
            coords[coords_len++].set(c.n, c.m);
            B.set_var(c, coords_len);
        }
};


coord_list::coord_list(unsigned sz)
{
    coords = new coord[sz];
    coords_size = sz;
    coords_len = 0;
}

coord_list::~coord_list()
{
    delete[] coords;
}

void knight_moves(const coord &p, coord_list &L)
{
    L.clear();

    L.append(p.n+1, p.m+2);
    L.append(p.n-1, p.m+2);
    L.append(p.n+1, p.m-2);
    L.append(p.n-1, p.m-2);

    L.append(p.n+2, p.m+1);
    L.append(p.n-2, p.m+1);
    L.append(p.n+2, p.m-1);
    L.append(p.n-2, p.m-1);
}

//
// Put vars in reachable order starting from square n, m
//
void initVars(char order, int n, int m)
{
    coord_list queue(unsigned(B.getN() * B.getM()));

    coord c;
    if ('r' == order) {
        for (c.n=1; c.n<=B.getN(); c.n++)
            for (c.m=1; c.m<=B.getM(); c.m++)
                queue.reach(c);
        return;
    }
    if ('c' == order) {
        for (c.m=1; c.m<=B.getM(); c.m++)
            for (c.n=1; c.n<=B.getN(); c.n++)
                queue.reach(c);
        return;
    }

    // Move based.

    c.set(n, m);
    queue.reach(c);

    coord_list moves(8);


    unsigned i = 0;
    while (i < queue.length()) {
        moves.clear();
        knight_moves(queue.get(i), moves);

        for (unsigned j=0; j<moves.length(); j++) {
            if (B.get_var(moves.get(j))) continue;
            queue.reach(moves.get(j));
        }
        i++;
    }

}

void showVars()
{
    using namespace std;

    cout << "Using variable order:\n";

    coord c;
    for (c.n=1; c.n<=B.getN(); c.n++) {
        for (c.m=1; c.m<=B.getM(); c.m++) {
            cout << setw(3) << B.get_var(c) << ' ';
        }
        cout << '\n';
    }
}

inline std::ostream& showVarOf(std::ostream &s, const coord &sq)
{
    return s << "x_" << std::setfill('0') << std::setw(2) << B.get_var(sq);
}


// Forest for multi-terminal MDDs
MEDDLY::forest* mtF;
// Forest for boolean constraints
MEDDLY::forest* boolF;

void buildConstraints(const coord &sq, MEDDLY::dd_edge &cons)
{
    using namespace std;
    using namespace MEDDLY;

    if (verbose) {
        cout << "--------------------------------------\n";
        cout << "Constraints for ";
        showVarOf(cout, sq) << "\n";
    }

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
    if (verbose) cout << "    Disjunction of:\n";
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
            dd_edge xsq(mtF);
            dd_edge one(mtF);
            mtF->createEdgeForVar(B.get_var(sq), false, xsq);
            mtF->createEdge(long(1), one);
            xsq += one;

            dd_edge xmv(mtF);
            mtF->createEdgeForVar(B.get_var(moves.get(j)), false, xmv);

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
        cons += ands;
    }
    if (debug) {
        cout << "\nEntire constraint:\n\n";
        cons.show(sout, 2);
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
    cout << "\t-o x\tSet ordering heuristic:\n";
    cout << "\t\t    -o r: by rows\n";
    cout << "\t\t    -o c: by columns\n";
    cout << "\t\t    -o m: by moves\n";
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

inline char charFrom(const char* arg)
{
    if (0==arg) return 0;
    return arg[0];
}

void process_args(int argc, const char** argv, coord &start)
{
    // Set defaults for everything
    int N = 4;
    int M = 4;
    char order='m';
    bool reverse = false;
    verbose = false;
    debug = false;
    start.n = 1;
    start.m = 1;
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
                            if (L<1 || L>10 || L2<1 || L2>10) {
                                std::cerr << "Value for -b out of range; ignoring\n";
                            } else {
                                start.n = L;
                                start.m = L2;
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
                            if (L<1 || L>10) {
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
                            if (L<1 || L>10) {
                                std::cerr << "Value " << L << " for -m out of range; ignoring\n";
                            } else {
                                M = L;
                            }
                        } else {
                            std::cerr << "Switch -m requires an argument\n";
                            exit(2);
                        }
                        continue;

            case 'o':   order = charFrom(argv[++i]);
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

    std::cout << "Knights tour on " << N << "x" << M << " chessboard\n";
    B.init(N, M);

    if (!B.in_bounds(start.n, start.m)) {
        std::cerr << "Starting position " << start << " is not on the board\n";
        exit(2);
    }

    initVars(order, start.n, start.m);
    if (reverse) B.reverse_vars();
    showVars();
}

int main(int argc, const char** argv)
{
    coord start;
    process_args(argc, argv, start);

    unsigned i = unsigned(B.getN() * B.getM());
    int* bounds = new int[i];
    while (i) {
        --i;
        bounds[i] = B.getN() * B.getM();
    }

    using namespace MEDDLY;

    initialize();
    domain* D = createDomainBottomUp(bounds, B.getN()*B.getM());

    mtF = D->createForest(false, range_type::INTEGER,
                    edge_labeling::MULTI_TERMINAL);

    boolF = D->createForest(false, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL);

    dd_edge* constr = new dd_edge[B.getN()*B.getM()];
    i = unsigned(B.getN() * B.getM());
    while (i) {
        --i;
        constr[i].setForest(boolF);

    }

    std::cout << "Building constraints per square" << std::endl;

    /*
     * Build constraints for each square, determining
     * what must come after
     */
    coord c;
    for (c.n=1; c.n<=B.getN(); c.n++)
        for (c.m=1; c.m<=B.getM(); c.m++) {
            buildConstraints(c, constr[(c.n-1)*B.getM()+c.m-1]);
        }



    MEDDLY::cleanup();
    return 0;
}

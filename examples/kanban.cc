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

#include "meddly/meddly.h"
#include "utils/rsgen_model.h"

class kanban : public rsgen_model {
    public:
        kanban();

        virtual bool process_arg(unsigned &i, const char** argv);
        virtual void show_usage(std::ostream &out, const char* arg0);

        virtual bool build_model(MEDDLY::output &out);
        virtual void build_initial(MEDDLY::minterm &init);

    protected:
        virtual bool known_cardinality(double &card) const;
#ifdef HAVE_LIBGMP
        virtual bool known_cardinality(mpz_t card) const;

        inline void axn(mpz_t item, unsigned a, unsigned pow) const {
            mpz_set_ui(item, a);
            for (unsigned i=0; i<pow; i++) {
                mpz_mul_ui(item, item, N);
            }
        }
#endif

    private:
        int N;
};

// ******************************************************************
// *                         kanban methods                         *
// ******************************************************************

kanban::kanban()
{
    N = -1;
}

// ******************************************************************

bool kanban::process_arg(unsigned &i, const char** argv)
{
    if (rsgen_model::process_arg(i, argv)) return true;

    if ('-' == argv[i][0])  return false;

    // assume this is the size
    N = atoi(argv[i]);
    i++;
    return true;
}

// ******************************************************************

void kanban::show_usage(std::ostream &out, const char* arg0)
{
    if (arg0) {
        out << "\nUsage: " << basename(arg0) << " [options] size\n\n";
    }
    out << "\tsize   : number of parts\n";

    rsgen_model::show_usage(out);
}

// ******************************************************************

bool kanban::build_model(MEDDLY::output &out)
{
    if (N<0) {
        std::cerr << "Missing or bad size.\n";
        return false;
    }
    out << "+-------------------------------------------+\n";
    out << "|   Initializing Kanban model for N = ";
    out.put(N, 4);
    out << "  |\n";
    out << "+-------------------------------------------+\n";
    out.flush();

    //
    // Initialize domain
    //
    int* sizes = new int[16];
    for (int i=15; i>=0; i--) sizes[i] = N+1;
    set_domain(MEDDLY::domain::createBottomUp(sizes, 16));
    delete[] sizes;

    //
    // Make the events
    //
    add_event(new simple_event("Tin1",    "X-+.............."));
    add_event(new simple_event("Tr1",     "X.-+............."));
    add_event(new simple_event("Tb1",     "X.+-............."));
    add_event(new simple_event("Tg1",     "X.-.+............"));
    add_event(new simple_event("Tr2",     "X.....-+........."));
    add_event(new simple_event("Tb2",     "X.....+-........."));
    add_event(new simple_event("Tg2",     "X.....-.+........"));
    add_event(new simple_event("Ts1_23",  "X+..--+..-+......"));
    add_event(new simple_event("Tr3",     "X.........-+....."));
    add_event(new simple_event("Tb3",     "X.........+-....."));
    add_event(new simple_event("Tg3",     "X.........-.+...."));
    add_event(new simple_event("Ts23_4",  "X....+..-+..--+.."));
    add_event(new simple_event("Tr4",     "X.............-+."));
    add_event(new simple_event("Tb4",     "X.............+-."));
    add_event(new simple_event("Tout4",   "X............+..-"));
    add_event(new simple_event("Tg4",     "X.............-.+"));

    return true;
}

// ******************************************************************

void kanban::build_initial(MEDDLY::minterm &init)
{
    init.setAllVars(0);
    init.setVar(1, N);
    init.setVar(5, N);
    init.setVar(9, N);
    init.setVar(13, N);
}

// ******************************************************************

bool kanban::known_cardinality(double &card) const
{
    double n = N;
    double n2 = n*n;
    double n3 = n2*n;
    double n4 = n3*n;
    double n5 = n4*n;

    double p1 = (n3 + 6*n2 + 11*n + 6) / 6;
    double p2 = (3*n5 + 30*n4 + 115*n3 + 210*n2 + 182*n + 60) / 60;

    card = p1*p1*p2;
    return true;
}

// ******************************************************************

#ifdef HAVE_LIBGMP

bool kanban::known_cardinality(mpz_t card) const
{
    mpz_t term;
    mpz_t p1;
    mpz_t p2;

    mpz_init(term);
    mpz_init(p1);
    mpz_init(p2);

    // p1 = (N^3 + 6 N^2 + 11 N + 6) / 6
    axn(p1, 1, 3);
    axn(term, 6, 2);
    mpz_add(p1, p1, term);
    axn(term, 11, 1);
    mpz_add(p1, p1, term);
    axn(term, 6, 0);
    mpz_add(p1, p1, term);
    mpz_div_ui(p1, p1, 6);

    // p2 = (3 N^5 + 30 N^4 + 115 N^3 + 210 N^2 + 182 N + 60) / 60;
    axn(p2, 3, 5);
    axn(term, 30, 4);
    mpz_add(p2, p2, term);
    axn(term, 115, 3);
    mpz_add(p2, p2, term);
    axn(term, 210, 2);
    mpz_add(p2, p2, term);
    axn(term, 182, 1);
    mpz_add(p2, p2, term);
    axn(term, 60, 0);
    mpz_add(p2, p2, term);
    mpz_div_ui(p2, p2, 60);

    mpz_mul(p1, p1, p1);
    mpz_mul(card, p1, p2);
    return true;
}

#endif

// ******************************************************************
// *                          main program                          *
// ******************************************************************

int main(int argc, const char** argv)
{
    kanban K;

    unsigned nArgs = (unsigned)argc;
    unsigned i=1;

    if (1 == argc) {
        K.show_usage(std::cerr, argv[0]);
        return 1;
    }

    while (i < nArgs) {
        if (!K.process_arg(i, argv)) {
            K.show_usage(std::cerr, argv[0]);
            return 1;
        }
    }

    try {
        MEDDLY::initialize();
        MEDDLY::ostream_output out(std::cout);

        if (!K.build_model(out)) {
            K.show_usage(std::cerr, argv[0]);
            return 1;
        }
        K.generate(out, true, true);

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
        std::cerr   << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    std::cerr << "\nSome other error?\n";
    return 4;

}


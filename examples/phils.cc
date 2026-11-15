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

class phils : public rsgen_model {
    public:
        phils();

        virtual bool process_arg(unsigned &i, const char** argv);
        virtual void show_usage(std::ostream &out, const char* arg0);

        virtual bool build_model(MEDDLY::output &out);
        virtual void build_initial(MEDDLY::minterm &init);

    private:
        /*
            Philosopher states:
                0   Idle
                1   Waiting on both forks
                2   Has left fork, wants right
                3   Has right fork, wants left
                4   Eating

            Fork states:
                0   Taken
                1   Available
         */

        rsgen_model::kron_event* Idle2Wait(int ph);
        rsgen_model::kron_event* Wait2HaveL(int ph);
        rsgen_model::kron_event* Wait2HaveR(int ph);
        rsgen_model::kron_event* HaveL2Eat(int ph);
        rsgen_model::kron_event* HaveR2Eat(int ph);
        rsgen_model::kron_event* Eat2Idle(int ph);

        /*
            TBD: allow different variable orders?

            Top, fork, phil, fork, phil, ...
        */

        // Get level number for philosopher i,  0 <= i < N.
        inline int Philosopher(int i) const
        {
            return i*2+1;
        }
        // Get level number for fork i,  0 <= i < N.
        inline int Fork(int i) const
        {
            return i*2+2;
        }

    private:
        std::vector <rsgen_model::level_delta> delta;
        int N;
};

// ******************************************************************
// *                         phils  methods                         *
// ******************************************************************

phils::phils()
{
    N = -1;
}

// ******************************************************************

bool phils::process_arg(unsigned &i, const char** argv)
{
    if (rsgen_model::process_arg(i, argv)) return true;

    if ('-' == argv[i][0])  return false;

    // assume this is the size
    N = atoi(argv[i]);
    i++;
    return true;
}

// ******************************************************************

void phils::show_usage(std::ostream &out, const char* arg0)
{
    show_synopsis(out, arg0, "[options] N");

    show_model_info(out, "Dining philosophers model.");

    show_section(out, "SIZE PARAMETERS");
    show_switch(out, "N", "Number of philosophers.");

    rsgen_model::show_usage(out);
}

// ******************************************************************

bool phils::build_model(MEDDLY::output &out)
{
    if (N<=0) {
        std::cerr << "Missing or bad size.\n";
        return false;
    }
    out << "+----------------------------------------------------+\n";
    out << "|   Initializing ";
    out.put(N, 6);
    out << " dining philosophers model   |\n";
    out << "+----------------------------------------------------+\n";
    out << "Variable order: Top, fork, phil, fork, phil, ...\n";
    out.flush();

    //
    // Initialize domain
    //
    int* sizes = new int[2*N];
    for (int i=2*N-1; i>=0; i--) sizes[i] = (i%2) ? 2 : 5;
    set_domain(MEDDLY::domain::createBottomUp(sizes, 2*N));
    delete[] sizes;

    //
    // Make the events
    //

    for (int i=0; i<N; i++) {
        add_event(Idle2Wait(i));
        add_event(Wait2HaveL(i));
        add_event(Wait2HaveR(i));
        add_event(HaveL2Eat(i));
        add_event(HaveR2Eat(i));
        add_event(Eat2Idle(i));
    }

    return true;
}

// ******************************************************************

void phils::build_initial(MEDDLY::minterm &init)
{
    for (int i=0; i<N; i++) {
        init.setVar( Philosopher(i), 0 );
        init.setVar( Fork(i), 1 );
    }
}

// ******************************************************************

rsgen_model::kron_event* phils::Idle2Wait(int ph)
{
    delta.resize(1);
    delta[0].set( Philosopher(ph), 0, 1);
    return new rsgen_model::kron_event("idle2wait", delta);
}

// ******************************************************************

rsgen_model::kron_event* phils::Wait2HaveL(int ph)
{
    delta.resize(2);
    delta[0].set( Philosopher(ph), 1, 2);
    delta[1].set( Fork(ph), 1, 0);
    return new rsgen_model::kron_event("wait2haveL", delta);
}

// ******************************************************************

rsgen_model::kron_event* phils::Wait2HaveR(int ph)
{
    delta.resize(2);
    delta[0].set( Philosopher(ph), 1, 3);
    delta[1].set( Fork((ph+N-1)%N), 1, 0);
    return new rsgen_model::kron_event("wait2haveR", delta);
}

// ******************************************************************

rsgen_model::kron_event* phils::HaveL2Eat(int ph)
{
    delta.resize(2);
    delta[0].set( Philosopher(ph), 2, 4);
    delta[1].set( Fork((ph+N-1)%N), 1, 0);
    return new rsgen_model::kron_event("haveL2eat", delta);
}

// ******************************************************************

rsgen_model::kron_event* phils::HaveR2Eat(int ph)
{
    delta.resize(2);
    delta[0].set( Philosopher(ph), 3, 4);
    delta[1].set( Fork(ph), 1, 0);
    return new rsgen_model::kron_event("haveR2eat", delta);
}

// ******************************************************************

rsgen_model::kron_event* phils::Eat2Idle(int ph)
{
    delta.resize(3);
    delta[0].set( Philosopher(ph), 4, 0);
    delta[1].set( Fork(ph),       0, 1);
    delta[2].set( Fork((ph+N-1)%N), 0, 1);
    return new rsgen_model::kron_event("eat2idle", delta);
}

// ******************************************************************
// *                          main program                          *
// ******************************************************************

int main(int argc, const char** argv)
{
    phils P;

    unsigned nArgs = (unsigned)argc;
    unsigned i=1;

    if (1 == argc) {
        P.show_usage(std::cerr, argv[0]);
        return 1;
    }

    while (i < nArgs) {
        if (!P.process_arg(i, argv)) {
            P.show_usage(std::cerr, argv[0]);
            return 1;
        }
    }

    try {
        MEDDLY::initialize();
        MEDDLY::ostream_output out(std::cout);

        if (!P.build_model(out)) {
            P.show_usage(std::cerr, argv[0]);
            return 1;
        }
        P.generate(out, true, true);

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

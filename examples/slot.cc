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

class slotted : public rsgen_model {
    public:
        slotted();

        virtual bool process_arg(unsigned &i, const char** argv);
        virtual void show_usage(std::ostream &out, const char* arg0);

        virtual bool build_model(MEDDLY::output &out);
        virtual void build_initial(MEDDLY::minterm_coll &init);

    private:
        rsgen_model::simple_event* Get(int i);
        rsgen_model::simple_event* Put(int i);
        rsgen_model::simple_event* Free(int i);
        rsgen_model::simple_event* Used(int i);
        rsgen_model::simple_event* Other(int i);
        rsgen_model::simple_event* Owner(int i);
        rsgen_model::simple_event* Write(int i);
        rsgen_model::simple_event* Go(int i);

    private:
        std::vector <rsgen_model::level_axb> delta;
        int N;
};

// ******************************************************************
// *                        slotted  methods                        *
// ******************************************************************

slotted::slotted()
{
    N = -1;
}

// ******************************************************************

bool slotted::process_arg(unsigned &i, const char** argv)
{
    if (rsgen_model::process_arg(i, argv)) return true;

    if ('-' == argv[i][0])  return false;

    // assume this is the size
    N = atoi(argv[i]);
    i++;
    return true;
}

// ******************************************************************

void slotted::show_usage(std::ostream &out, const char* arg0)
{
    show_synopsis(out, arg0, "[options] N");

    show_description(out, "MODEL DESCRIPTION", "(Safe) Petri net model of a slotted ring network protocol. See E. Pastor, O. Roig, J. Cortadella, and R. Badia, \"Petri net analysis using boolean manipulation\", 1994.");


    show_section(out, "SIZE PARAMETERS");
    show_switch(out, "N", "Number of network nodes. The overall Petri net has 8N places and 8N transitions.");

    rsgen_model::show_usage(out);
}

// ******************************************************************

bool slotted::build_model(MEDDLY::output &out)
{
    if (N<=0) {
        std::cerr << "Missing or bad size.\n";
        return false;
    }
    out << "+--------------------------------------------------+\n";
    out << "|   Initializing Slotted ring model for N = ";
    out.put(N, 4);
    out << "   |\n";
    out << "+--------------------------------------------------+\n";
    out.flush();

    //
    // Initialize domain
    //
    int* sizes = new int[8*N];
    for (int i=8*N-1; i>=0; i--) sizes[i] = 2;
    set_domain(MEDDLY::domain::createBottomUp(sizes, 8*N));
    delete[] sizes;

    //
    // Make the events
    //
    for (int i=0; i<N; i++) {
        add_event(Other(i));
        add_event(Owner(i));
        add_event(Write(i));
        add_event(Go(i));
        add_event(Get(i));
        add_event(Put(i));
        add_event(Used(i));
        add_event(Free(i));
    }

    return true;
}

// ******************************************************************

void slotted::build_initial(MEDDLY::minterm_coll &init)
{
    init.unused().setAllVars(0);
    for (int i=0; i<N; i++) {
        init.unused().setVar(i*8+3, 1);
        init.unused().setVar(i*8+5, 1);
    }
    init.pushUnused();
}

// ******************************************************************

rsgen_model::simple_event* slotted::Get(int i)
{
    delta.resize(4);
    delta[0].set(8*i+3, 1, +1);
    delta[1].set(8*i+5, 1, +1);
    delta[2].set(8*i+6, 1, -1);
    delta[3].set(8*i+8, 1, -1);
    return new rsgen_model::simple_event("Get", delta);
}

// ******************************************************************

rsgen_model::simple_event* slotted::Put(int i)
{
    delta.resize(4);
    delta[0].set(8*i+3, 1, +1);
    delta[1].set(8*i+4, 1, -1);
    delta[2].set(8*i+6, 1, -1);
    delta[3].set(8*i+7, 1, +1);
    return new rsgen_model::simple_event("Put", delta);
}

// ******************************************************************

rsgen_model::simple_event* slotted::Free(int i)
{
    delta.resize(4);
    delta[0].set(8*i+5, 1, -1);
    delta[1].set(8*i+6, 1, +1);
    int im1 = i ? (i-1) : (N-1);
    delta[2].set(8*im1+2, 1, +1);
    delta[3].set(8*im1+3, 1, -1);
    return new rsgen_model::simple_event("Free", delta);
}

// ******************************************************************

rsgen_model::simple_event* slotted::Used(int i)
{
    delta.resize(4);
    delta[0].set(8*i+7, 1, -1);
    delta[1].set(8*i+6, 1, +1);
    int im1 = i ? (i-1) : (N-1);
    delta[2].set(8*im1+3, 1, -1);
    delta[3].set(8*im1+1, 1, +1);
    return new rsgen_model::simple_event("Used", delta);
}

// ******************************************************************

rsgen_model::simple_event* slotted::Other(int i)
{
    delta.resize(2);
    delta[0].set(8*i+1, 1, -1);
    delta[1].set(8*i+4, 1, +1);
    return new rsgen_model::simple_event("Other", delta);
}

// ******************************************************************

rsgen_model::simple_event* slotted::Owner(int i)
{
    delta.resize(2);
    delta[0].set(8*i+1, 1, -1);
    delta[1].set(8*i+2, 1, +1);
    return new rsgen_model::simple_event("Owner", delta);
}

// ******************************************************************

rsgen_model::simple_event* slotted::Write(int i)
{
    delta.resize(2);
    delta[0].set(8*i+2, 1, -1);
    delta[1].set(8*i+4, 1, +1);
    return new rsgen_model::simple_event("Write", delta);
}


// ******************************************************************

rsgen_model::simple_event* slotted::Go(int i)
{
    delta.resize(2);
    delta[0].set(8*i+2, 1, -1);
    delta[1].set(8*i+8, 1, +1);
    return new rsgen_model::simple_event("Go", delta);
}

// ******************************************************************
// *                          main program                          *
// ******************************************************************

int main(int argc, const char** argv)
{
    slotted K;

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


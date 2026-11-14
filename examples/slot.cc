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
        virtual void build_initial(MEDDLY::minterm &init);

    private:
        int N;
};

// ******************************************************************
// *                         kanban methods                         *
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
    if (arg0) {
        out << "\nUsage: " << basename(arg0) << " [options] size\n\n";
    }
    out << "\tsize   : number of network nodes\n";

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

    // TBD

    return true;
}

// ******************************************************************

void slotted::build_initial(MEDDLY::minterm &init)
{
    init.setAllVars(0);
    for (int i=0; i<N; i++) {
        init.setVar(i*8+3, 1);
        init.setVar(i*8+5, 1);
    }
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


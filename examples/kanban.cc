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

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string.h>

#include "meddly/meddly.h"
#include "utils/rsgen_model.h"

class kanban : public rsgen_model {
    public:
        kanban();

        virtual bool process_arg(unsigned &i, const char** argv);
        virtual void show_usage(std::ostream &out, const char* arg0);

        virtual bool build_model();
        virtual void build_initial(MEDDLY::minterm &init);

    private:
        int N;

};

kanban::kanban()
{
    N = -1;
}

bool kanban::process_arg(unsigned &i, const char** argv)
{
    if (rsgen_model::process_arg(i, argv)) return true;

    if ('-' == argv[i][0])  return false;

    // assume this is the size
    N = atoi(argv[i]);
    i++;
    return true;
}

void kanban::show_usage(std::ostream &out, const char* arg0)
{
    if (arg0) {
        out << "\nUsage: " << basename(arg0) << " [options] size\n\n";
    }
    out << "\tsize   : number of parts\n";

    rsgen_model::show_usage(out);
}

bool kanban::build_model()
{
    if (N<0) {
        std::cerr << "Missing or bad size.\n";
        return false;
    }
    std::cout << "+-------------------------------------------+\n";
    std::cout << "|   Initializing Kanban model for N = "
              << std::left << std::setw(4) << N << "  |\n";
    std::cout << "+-------------------------------------------+" << std::endl;

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

    std::vector <level_axb> uplist;
    // "X-+..............",    // Tin1
    uplist.resize(2);
    uplist[0].set(1, 1, -1);
    uplist[1].set(2, 1, +1);
    add_event(new simple_event("Tin1", uplist));
    // "X.-+.............",    // Tr1
    uplist[0].set(2, 1, -1);
    uplist[1].set(3, 1, +1);
    add_event(new simple_event("Tr1", uplist));
    // "X.+-.............",    // Tb1
    uplist[0].set(2, 1, +1);
    uplist[1].set(3, 1, -1);
    add_event(new simple_event("Tb1", uplist));
    // "X.-.+............",    // Tg1
    uplist[0].set(2, 1, -1);
    uplist[1].set(4, 1, +1);
    add_event(new simple_event("Tg1", uplist));
    // "X.....-+.........",    // Tr2
    uplist[0].set(6, 1, -1);
    uplist[1].set(7, 1, +1);
    add_event(new simple_event("Tr2", uplist));
    // "X.....+-.........",    // Tb2
    uplist[0].set(6, 1, +1);
    uplist[1].set(7, 1, -1);
    add_event(new simple_event("Tb2", uplist));
    // "X.....-.+........",    // Tg2
    uplist[0].set(6, 1, -1);
    uplist[1].set(8, 1, +1);
    add_event(new simple_event("Tg2", uplist));
    // "X.........-+.....",    // Tr3
    uplist[0].set(10, 1, -1);
    uplist[1].set(11, 1, +1);
    add_event(new simple_event("Tr3", uplist));
    // "X.........+-.....",    // Tb3
    uplist[0].set(10, 1, +1);
    uplist[1].set(11, 1, -1);
    add_event(new simple_event("Tb3", uplist));
    // "X.........-.+....",    // Tg3
    uplist[0].set(10, 1, -1);
    uplist[1].set(12, 1, +1);
    add_event(new simple_event("Tg3", uplist));
    // "X.............-+.",    // Tr4
    uplist[0].set(14, 1, -1);
    uplist[1].set(15, 1, +1);
    add_event(new simple_event("Tr4", uplist));
    // "X.............+-.",    // Tb4
    uplist[0].set(14, 1, +1);
    uplist[1].set(15, 1, -1);
    add_event(new simple_event("Tb4", uplist));
    // "X............+..-",    // Tout4
    uplist[0].set(13, 1, +1);
    uplist[1].set(16, 1, -1);
    add_event(new simple_event("Tout4", uplist));
    // "X.............-.+"     // Tg4 TP
    uplist[0].set(14, 1, -1);
    uplist[1].set(16, 1, +1);
    add_event(new simple_event("Tg4", uplist));
    // "X+..--+..-+......",    // Ts1_23
    uplist.resize(6);
    uplist[0].set(1,  1, +1);
    uplist[1].set(4,  1, -1);
    uplist[2].set(5,  1, -1);
    uplist[3].set(6,  1, +1);
    uplist[4].set(9,  1, -1);
    uplist[5].set(10, 1, +1);
    add_event(new simple_event("Ts1_23", uplist));
    // "X....+..-+..--+..",    // Ts23_4
    uplist[0].set(5,  1, +1);
    uplist[1].set(8,  1, -1);
    uplist[2].set(9,  1, +1);
    uplist[3].set(12, 1, -1);
    uplist[4].set(13, 1, -1);
    uplist[5].set(14, 1, +1);
    add_event(new simple_event("Ts23_4", uplist));

    return true;
}

void kanban::build_initial(MEDDLY::minterm &init)
{
    init.setAllVars(0);
    init.setVar(1, N);
    init.setVar(5, N);
    init.setVar(9, N);
    init.setVar(13, N);
}


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

        if (!K.build_model()) {
            K.show_usage(std::cerr, argv[0]);
            return 1;
        }
        K.generate(true, true);

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


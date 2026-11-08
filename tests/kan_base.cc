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

#include "kan_base.h"

// ============================================

kanban_base::kanban_base()
{
    N = -1;
}

// ============================================

bool kanban_base::build_model(MEDDLY::output &out)
{
    if (N<0) {
        std::cerr << "Missing or bad size.\n";
        return false;
    }
    out << "Checking Kanban " << get_set_name() << ", N=" << N << "\n";

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

// ============================================

void kanban_base::build_initial(MEDDLY::minterm &init)
{
    init.setAllVars(0);
    init.setVar(1, N);
    init.setVar(5, N);
    init.setVar(9, N);
    init.setVar(13, N);
}


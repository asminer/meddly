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

//
// Test the randomization stuff in randomize.h/cc
//

#include "randomize.h"

void show(std::ostream &out, const std::vector <unsigned> &list)
{
    out << "[";
    for (unsigned i=0; i<list.size(); i++) {
        if (i) out << ", ";
        out << list[i];
    }
    out << "]";
}

int main(int argc, const char** argv)
{
    using namespace std;

    //
    // First argument: seed
    //
    long seed = 0;
    if (argv[1]) {
        seed = atol(argv[1]);
    }

    try {
        vectorgen_base::setSeed(seed);

        //
        // Test buildFullyFromIndex
        //

        vector <unsigned> L;

        vectorgen_base G(5, 2);
        cout << "Generator has potential " << G.potential() << "\n";

        G.buildFullyFromIndex(37, 0, 5, L);

        cout << "Got list ";
        show(cout, L);
        cout << "\n";

        return 0;
    }
    catch (const char* e) {
        cerr << "\nCaught our own error: " << e << "\n";
        return 1;
    }
}



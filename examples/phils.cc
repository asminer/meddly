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

// ******************************************************************
// *                          main program                          *
// ******************************************************************

int main(int argc, const char** argv)
{
    // test new kronecker events :)
    //

    std::vector <rsgen_model::level_delta> update;

    update.resize(10);

    update[ 0].set(3, 1, 4);
    update[ 1].set(1, 5, 9);
    update[ 2].set(2, 6, 5);
    update[ 3].set(3, 5, 8);
    update[ 4].set(9, 7, 9);
    update[ 5].set(3, 2, 3);
    update[ 6].set(8, 4, 6);
    update[ 7].set(6, 4, 2);
    update[ 8].set(3, 3, 8);
    update[ 9].set(2, 7, 3);

    rsgen_model::kron_event E("test", update);

    return 0;
}

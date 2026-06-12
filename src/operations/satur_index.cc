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

#include "satur_index.h"
#include "../variable.h"
#include "../domain.h"
#include "../forest.h"
#include "../rel_node.h"
#include "../io.h"

#include "satur_queue.h"

// #define TRACE

// **********************************************************************
// *                                                                    *
// *                     satur_index_basic  methods                     *
// *                                                                    *
// **********************************************************************

MEDDLY::satur_index_basic::satur_index_basic()
{
    curr = -1;
}

bool MEDDLY::satur_index_basic::
nextEdge(unsigned &i, unsigned &j, node_handle &down)
{
    for (;;) {
        if (queue.isEmpty()) return false;

        int ii = queue.front();
        MEDDLY_DCASSERT(ii>=0);
        i = unsigned(ii);

        if (curr < 0) {
            //
            // Starting fresh
            //
            curr = graph.getRowStart(i);
            down = graph.getDiagonal(i);
            if (down) {
                j = i;
                break;
            }
        }

        if (graph.getElement(curr, j, down)) {
            //
            // We got the next element; advance pointer
            //
            ++curr;
            break;
        }
        //
        // We're done with the row. Reset and loop
        //
        queue.remove();
        curr = -1;
    }

#ifdef TRACE
    std::cout << "level " << level << ": " << i << " -> " << j
              << " (down " << down << ")\n";
#endif
    return true;
}

void MEDDLY::satur_index_basic::show(output &s) const
{
    graph.show(s);
    s.put("updated: ");
    queue.show(s);
    s.put('\n');
}


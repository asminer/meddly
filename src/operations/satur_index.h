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

#ifndef MEDDLY_SATUR_INDEX_H
#define MEDDLY_SATUR_INDEX_H

#include "../defines.h"
#include <vector>

#include "satur_queue.h"
#include "satur_graph.h"

/*
 *  Interface for index explorer objects used in saturation.
 *
 *      constructor with no parameters
 *      destructor if needed
 *
 *      void attach(forest* F, int level, bool forwd):
 *
 *          called once, to attach the explorer to the MXD forest and level
 *
 *
 *      void restart(node handle n):
 *
 *          initialize the explorer from MXD node n
 *
 *
 *      void wasUpdated(unsigned i):
 *
 *          Indicate that index i has been updated, and its outgoing edges
 *          need to be re-visited.
 *
 *
 *      bool nextEdge(unsigned &i, unsigned &j, node_handle &down):
 *
 *          Get the next edge to explore, if there is one. Returns
 *          false if there are no more edges, otherwise returns true.
 *          On true return, we should explore edge [i, j] which has
 *          the specified downward pointer.
 *
 *
 *      node_handle getDiagonal(unsigned i):
 *
 *          Obtain edge [i,i]
 *
 *
 *      void show(output &s):
 *
 *          Display information, for debugging.
 */

namespace MEDDLY {
    class satur_index_nothing;
    class satur_index_basic;
}

// TBD: satur_index_old class that does NOT build the graph
// ******************************************************************
// *                                                                *
// *                   satur_index_nothing  class                   *
// *                                                                *
// ******************************************************************

// Just a queue; we explore by examining MxD nodes
class MEDDLY::satur_index_nothing {
        index_queue queue;

        rel_node* RN;
        unpacked_node* U;
        forest* For;
        node_handle node;
        int curr;

        int level;
        bool forwd;

    public:
        satur_index_nothing();
        ~satur_index_nothing();

        void attach(forest* F, int level, bool forwd);
        inline void restart(node_handle n) {
            if (n != node) _restart(n);
        }

        inline void wasUpdated(unsigned i) {
            queue.add(int(i));
        }

        bool nextEdge(unsigned &i, unsigned &j, node_handle &down);
        inline node_handle getDiagonal(unsigned i) {
            MEDDLY_DCASSERT(false);
            return 0;
        }
        void show(output &s) const;

    private:
        void _restart(node_handle n);
};

// ******************************************************************
// *                                                                *
// *                    satur_index_basic  class                    *
// *                                                                *
// ******************************************************************

// Basic graph and queue
class MEDDLY::satur_index_basic {
        satur_graph graph;
        index_queue queue;
        int curr;
    public:
        satur_index_basic();

        inline void attach(forest* F, int level, bool forwd) {
            graph.attach(F, level, forwd);
        }
        inline void restart(node_handle n) {
            curr = -1;
            graph.restart(n);
            queue.clear();
        }
        inline void wasUpdated(unsigned i) {
            graph.ensureRowExplored(i);
            queue.add(int(i));
        }
        bool nextEdge(unsigned &i, unsigned &j, node_handle &down);
        inline node_handle getDiagonal(unsigned i) {
            return graph.getDiagonal(i);
        }
        void show(output &s) const;
};

#endif

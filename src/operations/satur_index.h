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

// #define OLD_INDEX_INTERFACE

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
    class satur_index_basic;
}

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

// ******************************************************************
// *                                                                *
// *                         OLD  INTERFACE                         *
// *                                                                *
// ******************************************************************
#ifdef OLD_INDEX_INTERFACE

// Interface for index selection in saturation
//
//

namespace MEDDLY {
    class sat_index_explorer;
    class forest;
    class rel_node;
    class variable;
    class unpacked_node;
    class output;
};

class MEDDLY::sat_index_explorer {
    public:
        sat_index_explorer(forest* _F, int _level, bool _forwd);

        virtual ~sat_index_explorer();

        void restart(node_handle n);

        /**
            Get the next edge to fire, if there is one.
                @param  i       On output, source index.
                @param  j       On output, destination index. Might equal i.
                @param  down    On output, downward pointer.

                @return false, if there is nothing to fire;
                        true otherwise.
         */
        virtual bool nextEdge(unsigned &i, unsigned &j, node_handle &down) = 0;

        /**
            Get diagonal element [i,i].
            For forward exploration, if row i has not been built,
            we go ahead and build it.
         */
        inline node_handle getDiagonal(unsigned i) {
#ifdef DEVELOPMENT_CODE
            if (! rows.at(i).explored) {
                exploreRow(i);
            }
            return diagonals.at(i);
#else
            if (! rows[i].explored) {
                exploreRow(i);
            }
            return diagonals[i];
#endif
        }

        /**
            Indicate that index i has been updated, and its outgoing edges
            need to be re-visited.
         */
        inline void wasUpdated(unsigned i) {
#ifdef DEVELOPMENT_CODE
            if (! rows.at(i).explored) {
                exploreRow(i);
            }
#else
            if (! rows[i].explored) {
                exploreRow(i);
            }
#endif
            finishUpdate(i);
        }

        // For debugging
        virtual void show(output &s) const;

    protected:
        /**
            Clear out internal structures as necessary to start over
            from another relation node at the same level.
            Override this in derived classes as necessary. The default
            behavior does nothing.
            This method is called by method restart() and otherwise
            should not be called directly.
         */
        virtual void clear();

        void exploreRow(unsigned i);

        /**
            Update internal structures as necessary after row i is built.
            Override this in derived classes as necessary. The default
            behavior does nothing.
            This method is called by method exploreRow() and otherwise
            should not be called directly.
         */
        virtual void finishRow(unsigned i);

        /**
            Update internal structures as necessary after all rows are built.
            Override this in derived classes as necessary. The default
            behavior does nothing.
            This method is called by method restart() for backward
            exploration and otherwise should not be called directly.
         */
        virtual void finishAllRows();


        void expandRows(unsigned newsz);

        /**
            Update internal structures as necessary when the number
            of rows increases.
            Override this in derived classes as necessary. The default
            behavior does nothing.
            This method is called by method expandRows() and otherwise
            should not be called directly.
         */
        virtual void finishExpandRows(unsigned oldsz, unsigned newsz);

        /**
            Update internal structures as necessary when index i is updated.
            Override this in derived classes as necessary. The default
            behavior does nothing.
            This method is called by method wasUpdated() and otherwise
            should not be called directly.
         */
        virtual void finishUpdate(unsigned i);

    protected:
        struct row_element {
            unsigned index;
            node_handle down;
        };
        struct row_info {
            std::vector <row_element> elements;
            bool explored;

            row_info() {
                explored = false;
            }
        };

    protected:
        forest* For;
        int level;
        bool forwd;
        const variable* var;
        node_handle node;
        rel_node* RN;
        unpacked_node* U;

        std::vector <row_info> rows;
        std::vector <node_handle> diagonals;
};

// **********************************************************************
// *                                                                    *
// *                             Front  end                             *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {
    sat_index_explorer* makeSatIndexExplorer(char which, forest* F,
            int level, bool forwd);
};
#endif // ifdef OLD_INDEX_INTERFACE

#endif

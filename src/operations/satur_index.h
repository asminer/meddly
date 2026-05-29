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


#endif

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

#ifndef MEDDLY_REL_NODE
#define MEDDLY_REL_NODE

#include "defines.h"

namespace MEDDLY {
    class rel_node;
    class unpacked_node;
    class forest;
};

/**
    Abstract base class for a relation node.

    Allows different relation implementations to use the same
    implementation for relational product and saturation.

*/
class MEDDLY::rel_node {
    public:
        /**
            Initialize this node from handle n.
            Can be called more than once.
            Use node handle 0 to 'clear' a node.
         */
        virtual void initFromHandle(node_handle n) = 0;


        /**
            Get the outgoing edges from index i.
            When viewed as a matrix, this obtains row i.
                @param  i   (unprimed) variable value
                @param  u   will be filled with outgoing edges

                @return true, iff there is at least one outgoing edge
                        false, otherwise (in which case, u is left unchanged).
        */
        virtual bool outgoing(unsigned i, unpacked_node &u) = 0;

    protected:
        rel_node(forest* _parent);
        virtual ~rel_node();

        inline forest* getParent() const {
            return parent;
        };

        /// Check if the parent matches the given one
        bool hasParent(const forest* f) const {
            return (parent == f);
        }

        /// Delete a list
        static void deleteList(rel_node* p);

    private:
        forest* parent;
        rel_node* next;

        friend class forest;
};

#endif // include guard

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
};

/**
    Abstract base class for a relation node.

    Allows different relation implementations to use the same
    implementation for relational product and saturation.

*/
class MEDDLY::rel_node {
    public:
        /**
            Get the outgoing edges from index i.
            When viewed as a matrix, this obtains row i.
                @param  i   (unprimed) variable value
                @param  u   will be filled with outgoing edges

                @return true, iff there is at least one outgoing edge
        */
        virtual bool outgoing(unsigned i, unpacked_node &u) = 0;


        /**
            Get the incoming edges to index i.
            When viewed as a matrix, this obtains column i.
                @param  i   (primed) variable value
                @param  u   will be filled with outgoing edges

                @return true, iff there is at least one incoming edge
        */
        virtual bool incoming(unsigned i, unpacked_node &u) = 0;

    protected:
        rel_node();
        virtual ~rel_node();
};

#endif // include guard

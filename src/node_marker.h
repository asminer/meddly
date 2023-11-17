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

#ifndef MEDDLY_NODE_MARKER_H
#define MEDDLY_NODE_MARKER_H

#include "forest.h"

namespace MEDDLY {
    class node_marker;
};


/**
    Helper object used for anything requiring node marking.
*/
class MEDDLY::node_marker {
    public:
        node_marker(node_headers *H, const node_storage *nm);
        ~node_marker();

        inline void expand(size_t ns) {
            marked.expand(ns);
        }

        inline void shrink(size_t ns) {
            marked.shrink(ns);
        }

        inline void unmarkAll() {
            marked.clearAll();
        }

        inline bool isMarked(node_handle p) const {
            if (p<0) return false;
            return marked.get(p);
        }

        inline void mark(node_handle p) {
            if (p>0 && !marked.get(p)) _mark(p);
        }

        /// Return smallest non-terminal node handle >= i that is marked.
        inline node_handle nextMarked(node_handle i) const {
            if (i<1) return marked.firstOne(1);
            return (node_handle) marked.firstOne( (size_t) i );
        }

    private:
        void _mark(node_handle p);

    private:
        bitvector marked;
        node_headers* nodeHead;
        const node_storage* nodeMan;
};


#endif

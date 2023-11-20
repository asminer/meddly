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
        node_marker(bool permanent, node_headers &H, const node_storage *nm,
                expert_forest* f);
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

        /// Add a node to the queue to be explored.
        inline void addToQueue(node_handle p) {
            // TBD: try out an actual explore queue
            if (p>0 && !marked.get(p)) _mark(p);
        }

        /// Return smallest non-terminal node handle >= i that is marked.
        inline node_handle nextMarked(node_handle i) const {
            if (i<1) return marked.firstOne(1);
            return (node_handle) marked.firstOne( (size_t) i );
        }

        inline size_t getSize() const {
            return marked.getSize();
        }

        inline size_t countMarked() const {
            return marked.count();
        }

        /// Count the number of outgoing edges for marked nodes.
        size_t countEdges() const;

        /// Count the number of non-transparent outgoing edges for marked nodes.
        size_t countNonzeroEdges() const;

        /// Display all marked nodes, by levels.
        void showByLevels(output &s) const;

        /// TBD: writeNodeGraphPicture()

    private:
        void _mark(node_handle p);

    private:
        bitvector marked;
        node_headers& nodeHead;
        const node_storage* nodeMan;
        expert_forest* For;
};


#endif

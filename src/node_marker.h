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
#include "terminal.h"
#include <set>

namespace MEDDLY {
    class node_marker;
    class array_watcher;
};


/**
    Helper object used for anything requiring node marking.
*/
class MEDDLY::node_marker {
    public:
        node_marker(const forest &F, array_watcher* w = nullptr);
        ~node_marker();

        inline bool hasParent(const forest* f) const {
            return (&For) == f;
        }

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

        inline void setMarked(node_handle p) {
            if (p>0) {
                marked.set(p, true);
            }
        }

        inline void mark(node_handle p) {
            if (p>0 && !marked.get(p)) _mark(p);
        }

        /// Add a node to the queue to be explored.
        inline void addToQueue(node_handle p) {
            if (p>0 && !marked.get(p)) {
                MEDDLY_DCASSERT(S_top);
                marked.set(p, true);
                S_top->queue.push_back(p);
            }
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

        /**
            Add marked nodes at the specified level, to the list.
                @param  k       Level number to match.
                @param  v       Vector of nodes to add to.
        */
        void getNodesAtLevel(int k, std::vector <node_handle> &v) const;

        /**
            Determine all reached terminals.
            Visit all marked nodes, and for any non-transparent
            edge to a terminal node, add the terminal node to the set.
                @param  v       Vector of terminals to add to.
        */
        void getTerminals(std::set <node_handle> &v) const;

    private:
        void _mark(node_handle p);

    private:
        struct mystack {
            mystack* next;
            std::vector <node_handle> queue;
        };

    private:
        inline void push() {
            mystack* n;
            if (S_free) {
                n = S_free;
                S_free = S_free->next;
            } else {
                n = new mystack;
            }
            n->next = S_top;
            S_top = n;
        }
        inline void pop() {
            MEDDLY_DCASSERT(S_top);
            mystack* n = S_top->next;
            S_top->next = S_free;
            S_free = S_top;
            S_top = n;
        }

        void debug(const mystack *s);

    private:
        bitvector marked;
        const forest &For;

        mystack* S_top;
        mystack* S_free;

    private:
        inline bitvector* linkBits() {
            return &marked;
        }

        friend class MEDDLY::forest;

};


#endif

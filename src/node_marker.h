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

#include "defines.h"
#include "arrays.h"
#include <set>

namespace MEDDLY {
    class forest;
    class unpacked_node;

    class node_voyager;
    class node_counter;
    class edge_counter;
    class forest_writer;

    void exploreNodes(node_voyager &x,  node_handle p);
    void exploreNodes(node_counter &x,  node_handle p);
    void exploreNodes(edge_counter &x,  node_handle p);
    void exploreNodes(forest_writer &x, node_handle p);

#ifdef ALLOW_DEPRECATED_0_18_2
    class node_marker;
#endif
};

/*
    Objects used in node exploration should provide the following methods:

        const forest* getParent();

        static bool needsPreVisit();
        static bool needsPostVisit();

        bool isMarked(node_handle p) const;
        void setMarked(node_handle p);

        void preVisit(node_handle p, unpacked_node &u);
        void postVisit(node_handle p, unpacked_node &u);

    The above methods must support only the case p>0.

*/

/**
    Class to just visit (and mark) nodes.
*/
class MEDDLY::node_voyager {
    public:
        node_voyager(const forest* F, array_watcher* w = nullptr);
        ~node_voyager();

        inline const forest* getParent() const {
            return For;
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
            MEDDLY_DCASSERT(p>0);
            return marked.get(size_t(p));
        }

        inline void setMarked(node_handle p) {
            MEDDLY_DCASSERT(p>0);
            marked.set(size_t(p), true);
        }

        static inline bool needsPreVisit() {
            return false;
        }

        static inline bool needsPostVisit() {
            return false;
        }

        inline void preVisit(node_handle p, unpacked_node &u) {
            MEDDLY_DCASSERT(false);
        }

        inline void postVisit(node_handle p, unpacked_node &u) {
            MEDDLY_DCASSERT(false);
        }

        inline bitvector* linkBits() {
            return &marked;
        }

    private:
        bitvector marked;
        const forest* For;
};

// ==========================================================================

/*
 * OLD node_marker class below here
 */

#ifdef ALLOW_DEPRECATED_0_18_2

/**
    Helper object used for anything requiring node marking.
*/
class MEDDLY::node_marker {
    public:
        node_marker(const forest* F, array_watcher* w = nullptr);
        ~node_marker();

        inline bool hasParent(const forest* f) const {
            return f == For;
        }

        inline const forest* getParent() const {
            return For;
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
            return marked.get(size_t(p));
        }

        inline void setMarked(node_handle p) {
            if (p>0) {
                marked.set(size_t(p), true);
            }
        }

        inline void mark(node_handle p) {
            if (p>0 && !marked.get(size_t(p))) _mark(p);
        }

        /// Add a node to the queue to be explored.
        inline void addToQueue(node_handle p) {
            if (p>0 && !marked.get(size_t(p))) {
                MEDDLY_DCASSERT(S_top);
                marked.set(size_t(p), true);
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
        const forest* For;

        mystack* S_top;
        mystack* S_free;

    private:
        inline bitvector* linkBits() {
            return &marked;
        }

        friend class MEDDLY::forest;

};

#endif // allow_deprecated_0_18_2

#endif

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

#ifndef MEDDLY_DD_EDGE_H
#define MEDDLY_DD_EDGE_H

#include "defines.h"

namespace MEDDLY {
    class dd_edge;
    class forest;

    class input;
    class output;
};


/** New minimalist structure for handles to edges in forests.

    A dd_edge is a handle for user manipulation of functions stored in
    forests.  Based on the forest type, edges may or may not include
    edge values.

    There are a few useful operations that can be applied directly
    to a dd_edge; all the rest are done either through the "parent" forest,
    or through operations.  These include:

    - Deletion of a dd_edge.  This will cause the parent forest to recycle
      nodes as appropriate.

    - Checking for equality of two dd_edges, using the method equals().
*/
class MEDDLY::dd_edge {
    public:
        /// Construct and attach to a forest.
        dd_edge(forest* p=nullptr);

        /// Copy Constructor.
        dd_edge(const dd_edge &e);

        /// Assignment operator.
        dd_edge& operator=(const dd_edge &e);

        /// Destructor.  Will notify parent as appropriate.
        ~dd_edge();

        /// Attach to a forest.
        void attach(forest* p);

        /// Detach from the forest.
        inline void detach() { attach(nullptr); }

        /// Clear the edge but keep the forest
        inline void clear() {
            set(0);
            raw_value = 0;
        }

        /// Check if edges have the same parent forest
        inline bool sameForest(const dd_edge &e) const {
            return parentFID == e.parentFID;
        }

        /// Get the forest that contains us, or null
        forest* getForest() const;

        /// Check our parent
        inline bool isAttachedTo(const forest* p) const {
            return getForest() == p;
        }

        /// Get this dd_edge's label
        inline const char* getLabel() const { return label; }

        /** Set the edge's label.
            @param  L   Label to use; will be copied.
        */
        void setLabel(const char* L);

        inline long getNode() const { return node; }
        inline unsigned long getEdgeRaw() const { return raw_value; }
        inline const unsigned long* getEdgePtr() const { return &raw_value; }
        inline long getEdgeInt() const { return edge_int; }
        inline float getEdgeFloat() const { return edge_float; }

        inline void getEdgeValue(long &v) const { v = edge_int; }
        inline void getEdgeValue(float &v) const { v = edge_float; }

        /** Counts the number of unique nodes below this edge.
            @return     The number of unique nodes starting at the root node
                        of this dd_edge.
        */
        unsigned long getNodeCount() const;

        /** Counts the number of unique edges in this decision diagram.
            @param  countZeroes
                        if true, the stored zero edges are also counted
                        (sparse nodes do not store zero edges, so this
                        does not effect them; truncated nodes do store
                        some zero edges, so those edges will be counted).

            @return     the number of unique edges starting at the root node
                        of this dd_edge.
        */
        unsigned long getEdgeCount(bool countZeroes = false) const;


        /** Get this dd_edge's level.
            @return         the level.
        */
        int getLevel() const;

        /** Check for equality.
            @return true    iff this edge has the same parent and refers to
                            the same edge as \a e.
        */
        inline bool operator==(const dd_edge& e) const {
            return equals(e);
        }

        /** Check for inequality.
            @return true    iff this edge does not refer to the
                            same edge as \a e.
        */
        inline bool operator!=(const dd_edge& e) const {
            return !equals(e);
        }

        /**
            Display the edge information, compactly.
        */
        void show(output &s) const;

        /**
            Display the graph rooted at this node.
        */
        void showGraph(output &s) const;

        /** Draws a pictographical representation of the graph with
            this node as the root.
                @param  filename  Name of output file (without extension)
                @param  extension File extension (without ".").
                        E.g. "dot", "pdf", "ps"
        */
        void writePicture(const char* filename, const char* extension) const;

        //
        // Methods that will soon be replaced?
        // Or at least, made private?
        // Added here to ease the transition
        //
        void set(node_handle n);
        void set_and_link(node_handle n);

        inline void setEdgeValue(long value) {
            edge_int = value;
        }

        inline void setEdgeValue(float value) {
            edge_float = value;
        }

        inline void set(node_handle n, int value) {
            set(n);
            setEdgeValue(long(value));
        }
        inline void set(node_handle n, long value) {
            set(n);
            setEdgeValue(value);
        }
        inline void set(node_handle n, float value) {
            set(n);
            setEdgeValue(value);
        }

    private:
        void init(const dd_edge &e);

        inline bool equals(const dd_edge e) const {
            return
                (parentFID == e.parentFID) &&
                (node == e.node) &&
                (raw_value == e.raw_value);
        }

    private:
        char* label;            // for displaying
        unsigned parentFID;     // ID of parent forest
        unsigned index;         // our index in the parent forest
        node_handle node;
        union {
            long edge_int;
            float edge_float;

            unsigned long raw_value;    // must be at least as large
                                        // as the largest union element.
        };

        friend class forest;
        // Add any edge-valued forests here
        friend class evmdd_pluslong;
        friend class evmxd_pluslong;
        friend class evmdd_timesreal;
        friend class evmxd_timesreal;

        friend class unpacked_node; // maybe?
};

#endif // include guard

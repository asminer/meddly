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
#include "edge_value.h"

#include <string>

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

        /*
        /// Clear the edge but keep the forest
        inline void clear() {
            set(0);
            edgeval.set();
        }
        */

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
        inline const char* getLabel() const { return label.c_str(); }

        /** Set the edge's label.
            @param  L   Label to use; will be copied.
        */
        inline void setLabel(const char* L) {
            label = L;
        }

        inline long getNode() const { return node; }
        inline const edge_value& getEdgeValue() const { return edgeval; }
        inline long getEdgeInt() const { return edgeval.getLong(); }
        inline float getEdgeFloat() const { return edgeval.getFloat(); }

        inline void getEdgeValue(long &v) const { edgeval.get(v); }
        inline void getEdgeValue(float &v) const { edgeval.get(v); }

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

#ifdef ALLOW_DEPRECATED_0_17_3

        /** DEPRECATED; use dot_maker object instead (see io_dot.h).
            Draws a pictographical representation of the graph with
            this node as the root.
                @param  filename  Name of output file (without extension)
                @param  extension File extension (without ".").
                        E.g. "dot", "pdf", "ps"
        */
        void writePicture(const char* filename, const char* extension) const;
#endif

        /**
            Write the edge information to a file (stream).
                @param  s       Stream to write to
                @param  map     Mapping for non-terminal handles,
                                or null for identity mapping.
        */
        void write(output &s, node_handle* map) const;

        /**
            Read the edge information from a file (stream).
                @param  s       Stream to read from
                @param  map     Mapping for non-terminal handles,
                                or null for identity mapping.
        */
        void read(input &s, node_handle* map);

        //
        // Methods that will soon be replaced?
        // Or at least, made private?
        // Added here to ease the transition
        //
        void set(node_handle n);
        void set_and_link(node_handle n);

        inline void setEdgeValue(long value) {
            edgeval.set(value);
        }

        inline void setEdgeValue(float value) {
            edgeval.set(value);
        }

        inline edge_value& setEdgeValue() {
            return edgeval;
        }

        inline void set(node_handle n, long value) {
            set(n);
            setEdgeValue(value);
        }
        inline void set(node_handle n, float value) {
            set(n);
            setEdgeValue(value);
        }
        inline void set(node_handle n, edge_value v) {
            set(n);
            edgeval = v;
        }

    private:
        void init(const dd_edge &e);

        inline bool equals(const dd_edge e) const {
            return
                (parentFID == e.parentFID) &&
                (node == e.node) &&
                (edgeval == e.edgeval);
        }

    //
    // Actual edge information
    //
    private:
        std::string label;            // for displaying
        unsigned parentFID;     // ID of parent forest
        node_handle node;
        edge_value edgeval;

        friend class unpacked_node;

    //
    // for the dd_edge registry in the parent forest
    //
    private:
        dd_edge* prev;          // previous edge in forest registry
        dd_edge* next;          // next edge in forest registry
        friend class forest;
        friend class expert_forest; // until we move 'markAllRoots'

};

#endif // include guard

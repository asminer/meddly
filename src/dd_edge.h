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
#include "rangeval.h"

#include <string>
#include <vector>

namespace MEDDLY {
    class dd_edge;
    class iterator_helper;

    class forest;

    class input;
    class output;

    class minterm;
    class unpacked_node;
    class rangeval;
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

        class iterator {
            private:
                void init_with_forest(const forest* _F);

            public:
                iterator();

                iterator(const dd_edge &E, const minterm* mask=nullptr);

                // iterate from a random location.
                // rngfunc(a) returns a random integer in the range
                // [0, ..., a-1]
                iterator(const dd_edge &E, unsigned (*rngfunc)(unsigned));

                // Move constructor
                iterator(iterator &&I);

                ~iterator();

                //
                // Restart on a new edge, in the same forest though.
                //
                void restart(const dd_edge &E, const minterm* mask=nullptr);

                inline operator bool() const {
                    return !atEnd;
                }

                inline void operator++() {
                    if (atEnd)      return;
                    next();
                }
                inline void operator++(int) {
                    if (atEnd)      return;
                    next();
                }
                inline bool operator!=(const iterator& I) const {
                    if (I.atEnd && atEnd) return false;
                    if (I.atEnd || atEnd) return true;
                    return !equals(I);
                }
                inline bool operator==(const iterator& I) const {
                    if (I.atEnd && atEnd) return true;
                    if (I.atEnd || atEnd) return false;
                    return equals(I);
                }

                inline const minterm& operator*() const {
                    if (atEnd) {
                        throw error(error::INVALID_ITERATOR, __FILE__, __LINE__);
                    }
                    return *M;
                }


            private:
                iterator(const iterator &i) = delete;
                void operator=(const iterator &i) = delete;

                //
                // Increment
                //
                void next();

                bool equals(const iterator &I) const;

                friend class MEDDLY::dd_edge;
                friend class MEDDLY::iterator_helper;

            private:
                /*
                 * Keep a set of unpacked nodes at each level.
                 * If the variable is free (DONT_CARE), then
                 * the nodes are in sparse format for quick traversal.
                 * Otherwise we have a null pointer.
                 */
                unpacked_node** U_from;
                unpacked_node** U_to;

                /*
                 * Current index at each level.
                 * For free variables, this is the nonzero index
                 * into the sparse unpacked node.
                 * Otherwise, it is the index into the full unpacked
                 * node, and the index should equal the variable value.
                 */
                unsigned* Z_from;
                unsigned* Z_to;

                /*
                 * Accumulated edge values through each level.
                 */
                edge_value* ev_from;
                edge_value* ev_to;

                /*
                 * Root edge
                 */
                edge_value  root_ev;
                node_handle root_node;

                const forest* F;
                minterm* M;
                const minterm* mask;

                bool atEnd;
        };

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

        inline node_handle getNode() const { return node; }
        inline void xferNode(node_handle &h) {
            h = node; node = 0;
        }
        inline const edge_value& getEdgeValue() const { return edgeval; }
        inline int getEdgeInt() const { return edgeval.getInt(); }
        inline long getEdgeLong() const { return edgeval.getLong(); }
        inline float getEdgeFloat() const { return edgeval.getFloat(); }
        inline double getEdgeDouble() const { return edgeval.getDouble(); }

        inline void getEdgeValue(int &v) const { edgeval.get(v); }
        inline void getEdgeValue(long &v) const { edgeval.get(v); }
        inline void getEdgeValue(float &v) const { edgeval.get(v); }
        inline void getEdgeValue(double &v) const { edgeval.get(v); }

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

        /**
            Write the edge information to a file (stream).
                @param  s       Stream to write to
                @param  map     Translation from node handle to file node#
        */
        void write(output &s, const std::vector <unsigned> &map) const;

        /**
            Read the edge information from a file (stream).
                @param  s       Stream to read from
                @param  map     Translation from file node# to node handle
        */
        void read(input &s, const std::vector <node_handle> &map);

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

        /*
        inline void set(node_handle n, long value) {
            set(n);
            setEdgeValue(value);
        }
        inline void set(node_handle n, float value) {
            set(n);
            setEdgeValue(value);
        }
        */

        inline void set(edge_value v, node_handle n) {
            set(n);
            edgeval = v;
        }

        ///
        /// Build an end iterator
        ///
        inline iterator end() const
        {
            return iterator();
        }

        ///
        /// Build an iterator through this function
        /// The mask tells which values we can iterate over
        /// (set to DONT_CARE), while the others remain fixed.
        /// A NULL mask corresponds to all DONT_CAREs.
        ///
        inline iterator begin(const minterm *mask = nullptr) const
        {
            return iterator(*this, mask);
        }

        ///
        /// Default random number generator.
        ///
        /// Returns a "random" integer in the range [0, ..., a-1].
        ///
        static unsigned deflt_RNG(unsigned a);


        ///
        /// Build an iterator that starts at a random path
        /// but only among random paths that do not reach terminal 0.
        ///
        inline iterator random(unsigned (*RNG)(unsigned) = deflt_RNG) const
        {
            return iterator(*this, RNG);
        }

        /**
            For an index set only,
            find the set of assignments corresponding to the given index.

                @param  index   Desired index.
                @param  m       Output: the minterm such that f(m) = index.

                @throws     INVALID_OPERATION, if this is not an
                            Index Set EV+MDD.

                @return true,   if we found a minterm; false otherwise.
        */
        inline bool getElement(long index, minterm &m) const
        {
            switch (edgeval.getType()) {
                case edge_type::INT:
                    return  getElemInt(index, m);

                case edge_type::LONG:
                    return  getElemLong(index, m);

                default:
                    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
            }
        }


        ///
        /// Evaluate the function for variable assignments given in m,
        /// and store the result in val.
        /// This is the safest version, as it will correctly handle
        /// special values like infinity.
        ///
        void evaluate(const minterm &m, rangeval &val) const;

        ///
        /// Evaluate and convert to a 'standard' range.
        ///
        template <typename T>
        inline void evaluate(const minterm &m, T &val) const
        {
            rangeval v;
            evaluate(m, v);
            val = v;
        }
    private:
        bool getElemInt(long index, minterm &m) const;
        bool getElemLong(long index, minterm &m) const;

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

    //
    // for the dd_edge registry in the parent forest
    //
    private:
        dd_edge* prev;          // previous edge in forest registry
        dd_edge* next;          // next edge in forest registry
        friend class forest;

};

#endif // include guard

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
    class binary_operation;
};

// #define NEW_DD_EDGES

#ifdef NEW_DD_EDGES

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

        /// Check our parent
        inline bool isAttachedTo(const forest* p) const {
            return p == parent;
        }

        /// For now - get our parent
        /// (try to avoid the need for this)
        inline forest* getParent() const {
            return parent;
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

        //
        // Methods that will soon be replaced?
        // Or at least, made private?
        // Added here to ease the transition
        //
        void set(node_handle n);

        inline void setEdgeValue(long value) {
            edge_int = value;
        }

        inline void setEdgeValue(float value) {
            edge_float = value;
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
                (parent == e.parent) &&
                (node == e.node) &&
                (raw_value == e.raw_value);
        }


    private:
        char* label;    // for displaying
        forest *parent;
        node_handle node;
        union {
            long edge_int;
            float edge_float;

            unsigned long raw_value;    // must be at least as large
                                        // as the largest union element.
        };

        friend class forest;
        friend class unpacked_node; // maybe?
};

#else // not NEW_DD_EDGES

/** Structure for handles to edges in forests.

    A dd_edge is a handle for user manipulation of functions stored in
    forests.

    There are a few useful operations that can be applied directly
    to a dd_edge; all the rest are done either through the "parent" forest,
    or through operations.  These include:

    - Deletion of a dd_edge.  This will cause the parent forest to recycle
      nodes as appropriate.

    - Checking for equality of two dd_edges, using the method equals().
*/
class MEDDLY::dd_edge {
  public:
    /// Empty constructor.
    dd_edge();

    /** Constructor.
        Creates an empty edge in forest \a p.
        @param  p     forest to which this dd_edge will belong to.
    */
    dd_edge(forest* p);

    /** Copy Constructor.
        @param  e       dd_edge to copy.
    */
    dd_edge(const dd_edge &e);

    /** Assignment operator.
        @param  e       dd_edge to copy.
        @return         the new dd_edge.
    */
    dd_edge& operator=(const dd_edge &e);

    /// Destructor.  Will notify parent as appropriate.
    ~dd_edge();

  private:
    void init(const dd_edge &e);
    void destroy();


  public:

    /** Clears the contents of this edge. It will belong to the same
        forest as before.
    */
    void clear();

    /** Obtain a modifiable copy of the forest owning this edge.
        @return         the forest to which this edge belongs to.
    */
    forest* getForest() const;

        /// Check our parent
        inline bool isAttachedTo(const forest* p) const {
            return p == parent;
        }


    /// Set the forest owning this edge.
    void setForest(forest* f);

    /** Get this dd_edge's node handle.
        @return         the node handle.
    */
    node_handle getNode() const;

    /// Get this dd_edge's edge value (only valid for edge-valued MDDs).
//    void getEdgeValue(int& ev) const;
    void getEdgeValue(long& ev) const;

    /// Get this dd_edge's edge value (only valid for edge-valued MDDs).
    void getEdgeValue(float& ev) const;

    /// Get this dd_edge's label
    const char* getLabel() const;

    /** Get this dd_edge's level.
        @return         the level.
    */
    int getLevel() const;

    /** Get node cardinality.
        Provided for backward compatibility.
        Use apply(CARDINALITY, ...) instead.
        @return         the cardinality of the node.
    */
    double getCardinality() const;

    /** Counts the number of unique nodes in this decision diagram.
        @return       the number of unique nodes starting at the root node
                      of this dd_edge.
    */
    unsigned getNodeCount() const;

    /** Counts the number of unique edges in this decision diagram.
        @param  countZeroes
                      if true, the stored zero edges are also counted
                      (sparse nodes do not store zero edges, so this
                      does not effect them; truncated nodes do store
                      some zero edges, so those edges will be counted).
        @return       the number of unique edges starting at the root node
                      of this dd_edge.
    */
    unsigned getEdgeCount(bool countZeroes = false) const;

    /** Modifies the dd_edge fields.
        The dd_edge is cleared (it will still belong to the same forest),
        and the dd_edge data is filled with the data provided.
        @param  node    node handle.
    */
    void set(node_handle node);

    /** Modifies the dd_edge fields.
        The dd_edge is cleared (it will still belong to the same forest),
        and the dd_edge data is filled with the data provided.
        @param  node    node handle.
        @param  value   value of edge coming into the node (only useful
                        for edge-valued MDDs)
    */
    void set(node_handle node, int value);
    void set(node_handle node, long value);

    /** Modifies the dd_edge fields.
        The dd_edge is cleared (it will still belong to the same forest),
        and the dd_edge data is filled with the data provided.
        @param  node    node handle.
        @param  value   value of edge coming into the node (only useful
                        for edge-valued MDDs)
    */
    void set(node_handle node, float value);

    /** Modifies the edge value only.
        @param  value  value of edge coming into the node (only useful
                       for edge-valued MDDs)
     */
    void setEdgeValue(int value);
    void setEdgeValue(long value);
    void setEdgeValue(float value);

    /** Set the edge's label.
        @param  L   Label to use; will be copied.
    */
    void setLabel(const char* L);

    /** Check for equality.
        @return true    iff this edge has the same parent and refers to
                        the same edge as \a e.
    */
    bool operator==(const dd_edge& e) const;

    /** Check for inequality.
        @return true    iff this edge does not refer to the same edge as \a e.
    */
    bool operator!=(const dd_edge& e) const;

    /** Plus operator.
        BOOLEAN forests: Union; INTEGER/REAL forests: Addition.
        @param  e       dd_edge to Union/Add with this dd_edge.
        @return         \a this + \a e.
    */
    const dd_edge operator+(const dd_edge& e) const;

    /** Compound Plus operator.
        BOOLEAN forests: Union; INTEGER/REAL forests: Addition.
        This edge is overwritten with the result of the operation.
        @param  e       dd_edge to Union/Add with this dd_edge.
        @return         \a this + \a e.
    */
    dd_edge& operator+=(const dd_edge &e);

    /** Star operator.
        BOOLEAN forests: Intersection; INTEGER/REAL forests: Multiplication.
        @param  e       dd_edge to Intersection/Multiply with this dd_edge.
        @return         \a this * \a e.
    */
    const dd_edge operator*(const dd_edge& e) const;

    /** Compound Star operator.
        BOOLEAN forests: Intersection; INTEGER/REAL forests: Multiplication.
        This edge is overwritten with the result of the operation.
        @param  e       dd_edge to Intersection/Multiply with this dd_edge.
        @return         \a this * \a e.
    */
    dd_edge& operator*=(const dd_edge &e);

    /** Minus operator.
        BOOLEAN forests: Difference; INTEGER/REAL forests: Subtraction.
        @param  e       dd_edge for difference/subtract.
        @return         \a this - \a e.
    */
    const dd_edge operator-(const dd_edge& e) const;

    /** Compound Minus operator.
        BOOLEAN forests: Difference; INTEGER/REAL forests: Subtraction.
        This edge is overwritten with the result of the operation.
        @param  e       dd_edge for difference/subtract.
        @return         \a this - \a e.
    */
    dd_edge& operator-=(const dd_edge &e);

    /** Divide operator.
        BOOLEAN forests: INVALID; INTEGER/REAL forests: Division.
        @param  e       dd_edge for division.
        @return         \a this / \a e.
    */
    const dd_edge operator/(const dd_edge& e) const;

    /** Compound Divide operator.
        BOOLEAN forests: INVALID; INTEGER/REAL forests: Division.
        This edge is overwritten with the result of the operation.
        @param  e       dd_edge for division.
        @return         \a this / \a e.
    */
    dd_edge& operator/=(const dd_edge &e);

    /** Display the edge information.
        This is primarily for aid in debugging.
        Note that cardinality only works for MDDs, MTMDDs, MXDs and MTMXDs.
        @param  strm      File stream to write to.
        @param  verbosity 0: default
                          1: default + cardinality
                          2: default + displays graph rooted at this node.
                          3: default + cardinality + graph.
    */
    void show(output &s, int verbosity = 0) const;

    /** Draws a pictographical representation of the graph with this node as the root.
        @param  filename  Name of output file (without extension)
        @param  extension File extension (without "."). E.g. "pdf", "ps"
    */
    void writePicture(const char* filename, const char* extension) const;

    /// Write to a file
    void write(output &s, const node_handle* map) const;

    /// Read from a file
    void read(forest* p, input &s, const node_handle* map);

  private:
    friend class forest;
    friend class unpacked_node;

    void setIndex(unsigned ind);
    unsigned getIndex() const;

    forest *parent;
    unsigned index;  //  our slot number in the parent forest's list, or 0 for no slot.

    node_handle node;
    long raw_value;

    binary_operation* opPlus;
    binary_operation* opStar;
    binary_operation* opMinus;
    binary_operation* opDivide;

    // called when the parent is destroyed
    void orphan();

    // Label; used only for display purposes
    char* label;
};


// ******************************************************************
// *                                                                *
// *                    inlined dd_edge  methods                    *
// *                                                                *
// ******************************************************************

inline void MEDDLY::dd_edge::clear() {
  MEDDLY_DCASSERT(index);
  set(0);
  raw_value = 0;
}

inline MEDDLY::forest* MEDDLY::dd_edge::getForest() const {
  return parent;
}

inline MEDDLY::node_handle MEDDLY::dd_edge::getNode() const {
  return node;
}

inline const char* MEDDLY::dd_edge::getLabel() const {
  return label;
}

inline bool MEDDLY::dd_edge::operator==(const MEDDLY::dd_edge& e) const {
  if (parent != e.parent) return false;
  return (node == e.node) && (raw_value == e.raw_value);
}

inline bool MEDDLY::dd_edge::operator!=(const MEDDLY::dd_edge& e) const {
  return !(*this == e);
}

inline const MEDDLY::dd_edge
MEDDLY::dd_edge::operator+(const MEDDLY::dd_edge& e) const {
  return dd_edge(*this) += e;
}

inline const MEDDLY::dd_edge MEDDLY::dd_edge::operator*(const MEDDLY::dd_edge& e) const {
  return dd_edge(*this) *= e;
}

inline const MEDDLY::dd_edge MEDDLY::dd_edge::operator-(const MEDDLY::dd_edge& e) const {
  return dd_edge(*this) -= e;
}

inline void MEDDLY::dd_edge::setIndex(unsigned ind) {
  index = ind;
}

inline unsigned MEDDLY::dd_edge::getIndex() const {
  return index;
}

inline void MEDDLY::dd_edge::orphan() {
  parent = 0;
}

#endif

#endif // include guard

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

#ifndef MEDDLY_UNPACKED_NODE_H
#define MEDDLY_UNPACKED_NODE_H

#include "defines.h"
#include "edge_value.h"
#include "dd_edge.h"
#include "terminal.h"
#include "policies.h"

#include <cstring> // for memcpy

namespace MEDDLY {
    class forest;
    class initializer_list;
    class node_marker;

    class unpacked_node;
    struct unpacked_lists;

    class unreduced_node;
    struct unreduced_lists;
}

#define ALLOW_SET_FROM_DDEDGE
#define ALLOW_SET_FROM_POINTER
// #define DEBUG_UNPACKED_HASH

// ******************************************************************
// *                                                                *
// *                      unpacked_node  class                      *
// *                                                                *
// ******************************************************************

/**
    Class for unpacked nodes, i.e., copies of nodes outside a forest.
    Ideally - used anywhere we want to read node data, or create nodes.
    Unpacked nodes may be "full" or "sparse", independent of how
    the actual node is stored in the forest.

    Currently, an unpacked node is:
      array of node_handles, for the downward pointers
      array of integers, for indexes (sparse only)
      array of edge_values, for edge values
            (unless there are no edge values)
    All arrays are the same size, which could be larger than
    the requested size.

    If a node is extensible, the last element of the arrays
    is the extensible edge. In sparse storage, the index gives
    the starting index for the extensible edges.

    Unpacked nodes may be created and destroyed using the
    constructor and destructor as usual, or with static
    methods New() and Recycle(), which will pull from
    per-forest pools of Recycled nodes.

*/
class MEDDLY::unpacked_node {
        friend class initializer_list;
        friend struct unpacked_lists;
    protected:
        /** Constructor.
            The class must be "filled" by a forest before
            it can be used, however.
                @param  f   Parent forest.
                @param  fs  Should sparse storage be allowed

            Should not be called directly; use static
            method New() instead.
        */
        unpacked_node(const forest* f, node_storage_flags fs);

    public:
#ifdef ALLOW_DEPRECATED_0_17_8
        /** Constructor.
            The class must be "filled" by a forest before
            it can be used, however.
                @param  f   Parent forest.
                            If null, we must call attach()
                            before use.
        */
        unpacked_node(const forest* f = nullptr);
#endif

        /// Destructor.
        ~unpacked_node();

#ifdef ALLOW_DEPRECATED_0_17_8
        /// Attach the parent forest.
        void attach(const forest* f);
#endif

        /// Are we attached to f
        inline bool isAttachedTo(const forest* f) {
            return f == parent;
        }

    public:
        /* Initialization methods, primarily for reading */

        /**
            Initialize from a reduced node in our forest.
                @param  node    Node to unpack from
         */
        void initFromNode(node_handle node);

        /**
            Initialize as a redundant node.
                @param  k       Level of this node
                @param  node    All the down pointers
        */
        void initRedundant(int k, node_handle node);

        /**
            Initialize as a redundant node.
                @param  k       Level of this node
                @param  ev      All edge values
                @param  node    All the down pointers
        */
        inline void initRedundant(int k, const edge_value &ev, node_handle node)
        {
            if (ev.isVoid()) {
                initRedundant(k, node);
            } else {
                _initRedundant(k, ev, node);
            }
        }

    protected:
        void _initRedundant(int k, const edge_value &ev, node_handle node);

    public:

#ifdef ALLOW_DEPRECATED_0_17_8
        inline void initRedundant(const forest *f, int k,
                node_handle node, node_storage_flags fs)
        {
            ASSERT(__FILE__, __LINE__, isAttachedTo(f));
            ASSERT(__FILE__, __LINE__, (fs == nodestor) || (FULL_OR_SPARSE == nodestor));
            if (FULL_OR_SPARSE == nodestor) {
                is_full = (FULL_ONLY == fs);
            }
            initRedundant(k, node);
        }

        inline void initRedundant(const forest *f, int k,
                const edge_value &ev, node_handle node, node_storage_flags fs)
        {
            ASSERT(__FILE__, __LINE__, isAttachedTo(f));
            ASSERT(__FILE__, __LINE__, (fs == nodestor) || (FULL_OR_SPARSE == nodestor));
            if (FULL_OR_SPARSE == nodestor) {
                is_full = (FULL_ONLY == fs);
            }
            initRedundant(k, ev, node);
        }
#endif

        //
        // Build identity nodes
        //

        /**
            Initialize as an identity (singleton) node.
                @param  k       Level of this node; should be negative.
                @param  i       Index of non-zero edge
                @param  node    Down pointer for index i
        */
        void initIdentity(int k, unsigned i, node_handle node);

        /**
            Initialize as an identity (singleton) node.
                @param  k       Level of this node; should be negative.
                @param  i       Index of non-zero edge
                @param  ev      Edge value for index i
                @param  node    Down pointer for index i
        */
        inline void initIdentity(int k, unsigned i, const edge_value &ev,
                node_handle node)
        {
            if (ev.isVoid()) {
                initIdentity(k, i, node);
            } else {
                _initIdentity(k, i, ev, node);
            }
        }

    protected:
        void _initIdentity(int k, unsigned i, const edge_value &ev,
                node_handle node);

    public:

#ifdef ALLOW_DEPRECATED_0_17_8
        inline void initIdentity(const forest *f, int k, unsigned i,
                node_handle node, node_storage_flags fs)
        {
            ASSERT(__FILE__, __LINE__, isAttachedTo(f));
            ASSERT(__FILE__, __LINE__, (fs == nodestor) || (FULL_OR_SPARSE == nodestor));
            if (FULL_OR_SPARSE == nodestor) {
                is_full = (FULL_ONLY == fs);
            }
            initIdentity(k, i, node);
        }

        inline void initIdentity(const forest *f, int k, unsigned i,
                const edge_value &ev, node_handle node, node_storage_flags fs)
        {
            ASSERT(__FILE__, __LINE__, isAttachedTo(f));
            ASSERT(__FILE__, __LINE__, (fs == nodestor) || (FULL_OR_SPARSE == nodestor));
            if (FULL_OR_SPARSE == nodestor) {
                is_full = (FULL_ONLY == fs);
            }
            initIdentity(k, i, ev, node);
        }
#endif

    public:
        //
        // For convenience: get recycled instance and initialize
        //

        static inline unpacked_node* newFromNode(const forest *f,
                node_handle node, node_storage_flags fs)
        {
            unpacked_node* U = New(f, fs);
            ASSERT(__FILE__, __LINE__, U);
            U->initFromNode(node);
            return U;
        }

        static inline unpacked_node* newRedundant(const forest *f, int k,
                node_handle node, node_storage_flags fs)
        {
            unpacked_node* U = New(f, fs);
            ASSERT(__FILE__, __LINE__, U);
            U->initRedundant(k, node);
            return U;
        }

        static inline unpacked_node* newRedundant(const forest *f, int k,
                const edge_value &ev, node_handle node, node_storage_flags fs)
        {
            unpacked_node* U = New(f, fs);
            ASSERT(__FILE__, __LINE__, U);
            if (ev.isVoid()) {
                U->initRedundant(k, node);
            } else {
                U->_initRedundant(k, ev, node);
            }
            return U;
        }

        static inline unpacked_node* newIdentity(const forest *f, int k,
                unsigned i, node_handle node, node_storage_flags fs)
        {
            unpacked_node* U = New(f, fs);
            ASSERT(__FILE__, __LINE__, U);
            U->initIdentity(k, i, node);
            return U;
        }

        static inline unpacked_node* newIdentity(const forest *f, int k,
                unsigned i, const edge_value &ev, node_handle node,
                node_storage_flags fs)
        {
            unpacked_node* U = New(f, fs);
            ASSERT(__FILE__, __LINE__, U);
            if (ev.isVoid()) {
                U->initIdentity(k, i, node);
            } else {
                U->_initIdentity(k, i, ev, node);
            }
            return U;
        }

        /**
            Build a new, zeroed-out node for writing.
                @param  f       Forest the node will (maybe) get reduced in
                @param  lvl     Level of the node
                @param  tsz     (Initial) size of the node.
                @param  fs      Storage type (just here; the forest
                                may have other plans.)

                @return A pointer to a fresh node
         */
        static unpacked_node* newWritable(forest* f, int lvl,
                unsigned tsz, node_storage_flags fs);

        /**
            Build a new, zeroed-out node for writing.
                @param  f       Forest the node will (maybe) get reduced in
                @param  lvl     Level of the node. The node's size will
                                be the current size of the level.
                @param  fs      Storage type (just here; the forest
                                may have other plans.)

                @return A pointer to a fresh node
         */
        static unpacked_node* newWritable(forest* f, int lvl,
                node_storage_flags fs);

#ifdef ALLOW_DEPRECATED_0_17_8
        /** Create a zeroed-out full node of a given size */
        template <class T>
        static inline unpacked_node* newFull(forest *f,
                int levl, T tsz)
        {
            unpacked_node* U = New(f, FULL_ONLY);
            ASSERT(__FILE__, __LINE__, U);
            ASSERT(__FILE__, __LINE__, tsz >= 0);
            U->level = levl;
            U->resize(unsigned(tsz));
            U->clear(0, unsigned(tsz));
            U->setFull();
            U->allowWrites(f);
            return U;
        }

        /** Create a sparse node */
        template <class T>
        static inline unpacked_node* newSparse(forest *f,
                int levl, T nnzs)
        {
            unpacked_node* U = New(f, SPARSE_ONLY);
            ASSERT(__FILE__, __LINE__, U);
            U->level = levl;
            U->resize(unsigned(nnzs));
            U->setSparse();
            U->allowWrites(f);
            return U;
        }
#endif

    public:
        //
        // Display methods
        //

        /** Write a node in human-readable format.

            @param  s       Output stream.
            @param  details Should we show "details" or not.
        */
        void show(output &s, bool details) const;


        /** Write a node in machine-readable format.

            @param  s       Output stream.
            @param  map     Translation to use from node handle to file node#.
                            Allows us to renumber nodes as we write them.
        */
        void write(output &s, const std::vector <unsigned> &map) const;


        /** Read a node in machine-readable format.

            @param  s       Input stream.
            @param  map     Translation from file node# to node handles.
                            Allows us to renumber nodes as we read them.
        */
        void read(input &s, const std::vector <node_handle> &map);

    public:
        //
        // Node Access methods (inlined)
        //

        /// Get a pointer to the unhashed header data.
        inline const void* UHptr() const
        {
            ASSERT(__FILE__, __LINE__, extra_unhashed);
            return extra_unhashed;
        }

        /// Set the unhashed header data
        inline void setUHdata(const void* p)
        {
            ASSERT(__FILE__, __LINE__, p);
            memcpy(extra_unhashed, p, extra_unhashed_size);
        }

        /// Get the unhashed header data
        inline void getUHdata(void* p) const
        {
            ASSERT(__FILE__, __LINE__, p);
            memcpy(p, extra_unhashed, extra_unhashed_size);
        }

        /// Get the number of bytes of unhashed header data.
        inline unsigned UHbytes() const
        {
            return extra_unhashed_size;
        }

        /// Get a pointer to the hashed header data.
        inline const void* HHptr() const
        {
            ASSERT(__FILE__, __LINE__, extra_hashed);
            return extra_hashed;
        }

        /// Set the hashed header data
        inline void setHHdata(const void* p)
        {
            ASSERT(__FILE__, __LINE__, p);
            memcpy(extra_hashed, p, extra_hashed_size);
        }

        /// Get the hashed header data
        inline void getHHdata(void* p) const
        {
            ASSERT(__FILE__, __LINE__, p);
            memcpy(p, extra_hashed, extra_hashed_size);
        }

        /// Get the number of bytes of hashed header data.
        inline unsigned HHbytes() const
        {
            return extra_hashed_size;
        }

        /** Get a downward pointer.
            @param  n   Which pointer.
            @return     If this is a full node,
                        return pointer with index n.
                        If this is a sparse node,
                        return the nth non-zero pointer.
        */
        inline node_handle down(unsigned n) const
        {
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _down);
            return _down[n];
        }

        /** Get a downward pointer reference.
            @param  n   Which pointer.
            @return     If this is a full node,
                        return pointer with index n.
                        If this is a sparse node,
                        return the nth non-zero pointer.
        */
        inline node_handle& down(unsigned n)
        {
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _down);
            return _down[n];
        }


#ifdef ALLOW_DEPRECATED_0_17_8
        /** Get a downward pointer.
            @param  n   Which pointer.
            @return     If this is a full node,
                        return pointer with index n.
                        If this is a sparse node,
                        return the nth non-zero pointer.
        */
        inline node_handle down(int n) const
        {
            CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
            ASSERT(__FILE__, __LINE__, _down);
            return _down[n];
        }

        /** Get a downward pointer reference.
            @param  n   Which pointer.
            @return     If this is a full node,
                        return pointer with index n.
                        If this is a sparse node,
                        return the nth non-zero pointer.
        */
        inline node_handle& down(int n)
        {
            CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
            ASSERT(__FILE__, __LINE__, _down);
            return _down[n];
        }
#endif

        /** Get the index of the nth non-zero pointer.
            Use only for sparse nodes.
            @param  n   Which pointer
            @return     The index of the pointer
        */
        inline unsigned index(unsigned n) const
        {
            ASSERT(__FILE__, __LINE__, !is_full);
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _index);
            return _index[n];
        }

        /** Get a reference to the index of the nth non-zero pointer.
            Use only for sparse nodes.
            @param  n   Which pointer
            @return     The index of the pointer
        */
        inline unsigned& index(unsigned n)
        {
            ASSERT(__FILE__, __LINE__, !is_full);
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _index);
            return _index[n];
        }

#ifdef ALLOW_DEPRECATED_0_17_8
        /** Get the index of the nth non-zero pointer.
            Use only for sparse nodes.
            @param  n   Which pointer
            @return     The index of the pointer
        */
        inline unsigned index(int n) const
        {
            ASSERT(__FILE__, __LINE__, !is_full);
            CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
            ASSERT(__FILE__, __LINE__, _index);
            return _index[n];
        }

        /** Get a reference to the index of the nth non-zero pointer.
            Use only for sparse nodes.
            @param  n   Which pointer
            @return     The index of the pointer
        */
        inline unsigned& index(int n)
        {
            ASSERT(__FILE__, __LINE__, !is_full);
            CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
            ASSERT(__FILE__, __LINE__, _index);
            return _index[n];
        }
#endif

        /** Get the nth edge value.
            @param  n   Which pointer
            @return     The edge value
        */
        inline const edge_value& edgeval(unsigned n) const {
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _edge);
            return _edge[n];
        }

        /** Get a modifiable reference to the nth edge value.
            @param  n   Which pointer
            @return     The edge value
        */
        inline edge_value& edgeval(unsigned n) {
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _edge);
            return _edge[n];
        }

#ifdef ALLOW_DEPRECATED_0_17_8
        /** Get the nth edge value.
            @param  n   Which pointer
            @return     The edge value
        */
        inline const edge_value& edgeval(int n) const {
            CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
            ASSERT(__FILE__, __LINE__, _edge);
            return _edge[n];
        }

        /** Get a modifiable reference to the nth edge value.
            @param  n   Which pointer
            @return     The edge value
        */
        inline edge_value& edgeval(int n) {
            CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
            ASSERT(__FILE__, __LINE__, _edge);
            return _edge[n];
        }

        /** Subtract from an edge value.
            @param  n       Which pointer
            @param  v       Value to subtract
        */
        template <class T>
        inline void subtractFromEdge(unsigned n, T v) {
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _edge);
            _edge[n].subtract(v);
        }


        /** Divide an edge value.
            @param  n       Which pointer
            @param  v       Value to divide by
        */
        template <class T>
        inline void divideEdge(unsigned n, T v) {
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _edge);
            _edge[n].divide(v);
        }

        /// Get the nth edge value, as a long.
        inline long edge_long(unsigned n) const {
            return edgeval(n).getLong();
        }

        /// Get the edge value, as a float.
        inline float edge_float(unsigned n) const {
            return edgeval(n).getFloat();
        }

        /**
            Set just the edge value.
            If you're setting the down pointer also,
            probably you want to use setFull() or setSparse().
            @param  n   Which pointer
            @param  ev  New edge value
        */
        inline void setEdgeval(unsigned n, const edge_value &ev) {
            ASSERT(__FILE__, __LINE__, ev.hasType(the_edge_type));
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _edge);
            _edge[n] = ev;
        }
#endif

        /**
            Set a full edge.
                @param  n       Which pointer
                @param  h       Node handle
        */
        inline void setFull(unsigned n, node_handle h)
        {
            ASSERT(__FILE__, __LINE__, !hasEdges());
            ASSERT(__FILE__, __LINE__, isFull());
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _down);
            _down[n] = h;
        }

        /**
            Set a full edge.
                @param  n       Which pointer
                @param  v       Edge value
                @param  h       Node handle
        */
        inline void setFull(unsigned n, const edge_value &v, node_handle h)
        {
            ASSERT(__FILE__, __LINE__, v.hasType(the_edge_type));
            ASSERT(__FILE__, __LINE__, isFull());
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _down);
            _down[n] = h;
            if (_edge) {
                _edge[n] = v;
            } else {
                ASSERT(__FILE__, __LINE__, v.isVoid());
            }
        }

#ifdef ALLOW_SET_FROM_DDEDGE
        /**
            Set a full edge from E, and destroy E.
                @param  n       Which pointer
                @param  E       Edge value (if needed) and node
        */
        inline void setFull(unsigned n, dd_edge &E)
        {
            ASSERT(__FILE__, __LINE__, E.isAttachedTo(parent));
            ASSERT(__FILE__, __LINE__, E.getEdgeValue().hasType(the_edge_type));
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            E.xferNode(_down[n]);
            if (_edge) {
                _edge[n] = E.getEdgeValue();
            } else {
                ASSERT(__FILE__, __LINE__, E.getEdgeValue().isVoid());
            }
        }
#endif


#ifdef ALLOW_SET_FROM_POINTER
        /**
            Set a full edge.
                @param  n       Which pointer
                @param  p       Raw pointer for edge value
                @param  h       Node handle
        */
        inline void setFull(unsigned n, const void* p, node_handle h)
        {
            ASSERT(__FILE__, __LINE__, isFull());
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _down);
            ASSERT(__FILE__, __LINE__, _edge);
            _down[n] = h;
            _edge[n].set(the_edge_type, p);
        }
#endif


        /**
            Set a sparse edge.
                @param  n       Which nonzero edge
                @param  i       Index of the edge
                @param  h       Node handle
        */
        inline void setSparse(unsigned n, unsigned i, node_handle h)
        {
            ASSERT(__FILE__, __LINE__, !hasEdges());
            ASSERT(__FILE__, __LINE__, isSparse());
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _down);
            ASSERT(__FILE__, __LINE__, _index);
            _index[n] = i;
            _down[n] = h;
        }


        /**
            Set a sparse edge.
                @param  n       Which nonzero edge
                @param  i       Index of the edge
                @param  v       Edge value
                @param  h       Node handle
        */
        inline void setSparse(unsigned n, unsigned i, const edge_value &v,
                node_handle h)
        {
            ASSERT(__FILE__, __LINE__, v.hasType(the_edge_type));
            ASSERT(__FILE__, __LINE__, isSparse());
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _down);
            ASSERT(__FILE__, __LINE__, _index);
            _index[n] = i;
            _down[n] = h;
            if (_edge) {
                _edge[n] = v;
            } else {
                ASSERT(__FILE__, __LINE__, v.isVoid());
            }
        }


#ifdef ALLOW_SET_FROM_DDEDGE
        /**
            Set a sparse edge.
                @param  n       Which nonzero edge
                @param  i       Index of the edge
                @param  E       Edge value (if needed) and node
        */
        inline void setSparse(unsigned n, unsigned i, dd_edge &E)
        {
            ASSERT(__FILE__, __LINE__, E.getEdgeValue().hasType(the_edge_type));
            ASSERT(__FILE__, __LINE__, isSparse());
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _down);
            ASSERT(__FILE__, __LINE__, _index);
            E.xferNode(_down[n]);
            _index[n] = i;
            if (_edge) {
                _edge[n] = E.getEdgeValue();
            } else {
                ASSERT(__FILE__, __LINE__, E.getEdgeValue().isVoid());
            }
        }
#endif


#ifdef ALLOW_SET_FROM_POINTER
        /**
            Set a sparse edge.
                @param  n       Which nonzero edge
                @param  i       Index of the edge
                @param  p       Raw pointer for edge value
                @param  h       Node handle
        */
        inline void setSparse(unsigned n, unsigned i, const void* p,
                node_handle h)
        {
            ASSERT(__FILE__, __LINE__, isSparse());
            CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            ASSERT(__FILE__, __LINE__, _down);
            ASSERT(__FILE__, __LINE__, _index);
            ASSERT(__FILE__, __LINE__, _edge);
            _index[n] = i;
            _down[n] = h;
            _edge[n].set(the_edge_type, p);
        }
#endif


        /**
            Hack for mark and sweep.
            Add an extra, temporary node to the list of roots
            for mark and sweep.
            Do this for example in an "APPLY" operation that
            needs to create a temporary result.
            We can hold one such result and make sure it's marked,
            by calling this method.
                @param  t   Node handle for temporary result.
         */
        inline void setTempRoot(node_handle t)
        {
            mark_extra = t;
        }

    public:
        //
        // Methods to access the extensible portion of the node
        //
        // Note: Only nodes that belong to an extensible level can
        // be extensible.
        //
        // Extensible node  : (Regular part, Extensible part).
        // Regular Part     : same as non-extensible nodes.
        // Extensible Part  : a single edge <index, node, edge-value>,
        //                    for all indices in [extensible_index, +infinity].
        //
        // Example:
        // The node
        //    [0:n0:ev0, 1:n1:ev1, 2:n2:ev2, 3:n2:ev2, ..., +infinity:n2:ev2]
        // is represented as the following extensible node:
        //    ([0:n0:ev0, 1:n1:ev1], Extensible: [2:n2:ev2])
        // with,
        //    size of node           : 2
        //    extensible index       : 2
        //    extensible node_handle : n2
        //    extensible edge-value  : ev2
        //



        /// Does this node have an extensible edge?
        inline bool isExtensible() const
        {
#ifdef ALLOW_EXTENSIBLE
            return is_extensible;
#else
            return false;
#endif
        }

#ifdef ALLOW_EXTENSIBLE
        /// Set this node as extensible
        inline void markAsExtensible()
        {
            ASSERT(__FILE__, __LINE__, can_be_extensible);
            is_extensible = true;
        }

        /// Set this node as not extensible
        inline void markAsNotExtensible()
        {
            is_extensible = false;
        }

        /// Get extensible down pointer
        inline node_handle ext_d() const
        {
            ASSERT(__FILE__, __LINE__, isExtensible());
            return down(size-1);
        }

        /// Get the extensible index
        inline unsigned ext_i() const
        {
            ASSERT(__FILE__, __LINE__, isExtensible());
            return index(size - 1);
        }

        /// Get the extensible edge value
        inline const edge_value& ext_ev() const
        {
            ASSERT(__FILE__, __LINE__, isExtensible());
            return edgeval(size - 1);
        }
#endif

    public:
        //
        // Access of other information
        //

        /// Get the level number of this node.
        inline int getLevel() const
        {
            return level;
        }

        /// Set the level number of this node.
        inline void setLevel(int k)
        {
            level = k;
        }

        /// Get the size of this node.
        inline unsigned getSize() const
        {
            return size;
        }

        /// Is this a sparse node?
        inline bool isSparse() const
        {
            return !is_full;
        }

        /// Is this a full node?
        inline bool isFull() const
        {
            return is_full;
        }

        /// Does this node have edge values?
        inline bool hasEdges() const
        {
            return edge_type::VOID != the_edge_type;
        }

        /// Edge type
        inline edge_type getEdgeType() const {
            return the_edge_type;
        }

    public:
        /// Get the node's hash
        inline unsigned hash() const
        {
#ifdef DEVELOPMENT_CODE
            ASSERT(__FILE__, __LINE__, has_hash);
#endif
            return the_hash;
        }

        /// Compute the node's hash
        void computeHash();

#ifdef DEBUG_UNPACKED_HASH
        void debugHash(output &s) const;
#endif

        /// Checks if the node has no trailing redundant edges.
        /// I.e., the node does NOT have trailing edges that can
        /// be collapsed into the extensible edge.
        inline bool isTrim() const
        {
#ifdef ALLOW_EXTENSIBLE
            if (!isExtensible()) return true;
            if (size < 2) {
                // We have only the extensible edge.
                return true;
            }
            if (_down[size-1] != _down[size-2]) {
                // extensible pointer != trailing pointer
                return true;
            }
            // TBD: check edge values
            if (isSparse() && (_index[size-1] != _index[size-2]+1) ) {
                // There's a gap between the last normal pointer
                // and the extensible pointer
                return true;
            }
            // Still here? we can truncate edges.
            return false;
#else
            return true;
#endif
        }


        /// Removes redundant trailing edges.
        /// If the unpacked node is sparse, it assumes its indices
        /// to be in ascending order.
        void trim();

        /// If the unpacked node is sparse, it is sorted so that the
        /// indices are in ascending order.
        void sort();

        /// checks if the node indices are in ascending order
        bool isSorted() const;
    public:

        /// Change the size of a node
        inline void resize(unsigned ns) {
            size = ns;
            if (ns > alloc) expand(ns);
        }

        inline void resize(int ns) {
            ASSERT(__FILE__, __LINE__, ns>=0);
            resize(unsigned(ns));
        }

        /// Shrink the size of a node
        inline void shrink(unsigned ns)
        {
            ASSERT(__FILE__, __LINE__, ns <= getSize());
            size = ns;
        }

        /// Set the node as sparse.
        inline void setSparse() {
            is_full = false;
        }

        /// Set the node as truncated full.
        inline void setFull() {
            is_full = true;
        }

        /// Allow writing to this node.
        inline void allowWrites(forest* mp) {
            if (modparent == mp) return;
            ASSERT(__FILE__, __LINE__, nullptr == modparent);
            ASSERT(__FILE__, __LINE__, parent == mp);
            modparent = mp;
            AddToBuildList(this);
            mark_extra = 0;
        }


        // Set edges to transparent
        void clear(unsigned low, unsigned high);

        /// Was this node initialized from an identity node
        inline bool wasIdentity() const {
            return orig_was_identity;
        }

    protected:
        //
        // This is the ugliest hack that ever hacked
        // but it allows us to keep forest "hidden"
        // until later
        //
        template <class FORST>
        inline void _clear(const FORST &F, unsigned low, unsigned high)
        {
            CHECK_RANGE(__FILE__, __LINE__, 0u, low, alloc);
            CHECK_RANGE(__FILE__, __LINE__, 0u, high, alloc+1);
            ASSERT(__FILE__, __LINE__, _down);

            if (hasEdges()) {
                ASSERT(__FILE__, __LINE__, _edge);
                for (unsigned i=low; i<high; i++) {
                    F.getTransparentEdge(_edge[i], _down[i]);
                }
            } else {
                for (unsigned i=low; i<high; i++) {
                    _down[i] = F.getTransparentNode();
                }
            }
#ifdef DEVELOPMENT_CODE
            has_hash = false;
#endif
        }
        //
        // Special case of _clear for all edges
        //
        template <class FORST>
        inline void _clear(const FORST &F)
        {
            ASSERT(__FILE__, __LINE__, _down);

            if (hasEdges()) {
                ASSERT(__FILE__, __LINE__, _edge);
                for (unsigned i=getSize(); i; ) {
                    --i;
                    F.getTransparentEdge(_edge[i], _down[i]);
                }
            } else {
                for (unsigned i=getSize(); i; ) {
                    --i;
                    _down[i] = F.getTransparentNode();
                }
            }
#ifdef DEVELOPMENT_CODE
            has_hash = false;
#endif
        }

    public:
        //
        // Centralized per-forest (using FIDs) recycling
        //

        /// Pull a recycled node off of f's free list,
        /// or create a new one if needed.
        static unpacked_node* New(const forest* f, node_storage_flags ns);

        /// Add a node to the build list for the appropriate forest.
        static void AddToBuildList(unpacked_node* n);

        /// Mark children in writable nodes.
        ///     @param  M   node marker we should use to mark nodes;
        ///                 this also gives us the forest to check.
        ///
        static void MarkWritable(node_marker &M);

        /// Remove r from its forest's build list if needed,
        /// and add r back to its forest's recycle list.
        static void Recycle(unpacked_node* r);

        /// Update counts of children in writable nodes.
        ///     @param  F           Forest we care about.
        ///     @param  incounts    Vector of counts for each
        ///                         nonterminal node.
        ///
        static void AddToIncomingCounts(const forest* F,
                        std::vector <unsigned> &incounts);

        /// forest should call this in its constructor, after it has its FID.
        /// Initializes empty lists for forest f.
        static void initForest(const forest* f);

        /// forest should call this in its destructor.
        /// deletes all recycled nodes for forest f.
        static void doneForest(const forest* f);

    private:
        static void initStatics();
        static void doneStatics();

    private:
        /// Per-forest free and build lists.
        /// Stored as an array, indexed by the forest's FID.
        static unpacked_lists* ForLists;

        /// Allocated size of ForLists
        static unsigned ForListsAlloc;


    private:
        void expand(unsigned ns);

        static inline void deleteList(unpacked_node* &p) {
            while (p) {
                unpacked_node* n = p->next;
                delete p;
                p = n;
            }
        }

        static void showSingly(const unpacked_node* list);
        static void showDoubly(const unpacked_node* list);

    private:
        /// Next in list
        unpacked_node* next;
        /// Previous in list (for build lists only)
        unpacked_node* prev;

        /// Forest where the node belongs
        const forest* parent;
        /// Modifiable parent forest; required for writable nodes
        forest* modparent;

#ifdef ALLOW_DEPRECATED_0_17_8
        /// FID of the parent
        unsigned pFID;
#else
        /// FID of the parent
        const unsigned pFID;
#endif


        /// Down pointers
        node_handle* _down;

        /// Indexes, for sparse; otherwise unused
        unsigned* _index;

        /// Edge values; or null if void edges
        edge_value* _edge;

        /// Extra, temporary, node handle for marking.
        node_handle mark_extra;

        /// Allocated sizes of arrays
        unsigned alloc;

        /// Used sizes of arrays
        unsigned size;

        /// Extra header information that is not hashed
        void* extra_unhashed;

#ifdef ALLOW_DEPRECATED_0_17_8
        /// Number of bytes in extra unhashed header
        unsigned extra_unhashed_size;
#else
        /// Number of bytes in extra unhashed header
        const unsigned extra_unhashed_size;
#endif

        /// Extra header information that is hashed
        void* extra_hashed;

#ifdef ALLOW_DEPRECATED_0_17_8
        /// Number of bytes in extra hashed header
        unsigned extra_hashed_size;
#else
        /// Number of bytes in extra hashed header
        const unsigned extra_hashed_size;
#endif

        /// Level of the node
        int level;

        /// Hash of the node
        unsigned the_hash;

        /// Allowed storage
        /// For FULL_OR_SPARSE, we allocate _index
        /// but don't require the node to be sparse.
        const node_storage_flags nodestor;

        /// Are we using full storage?
        bool is_full;

#ifdef ALLOW_EXTENSIBLE
        /// Can this node be extensible?
        /// Determined by the parent forest.
        bool can_be_extensible;

        /// Is this node extensible?
        bool is_extensible;
#endif

#ifdef DEVELOPMENT_CODE
        /// only node pointers built by New() should be Recycle()d.
        bool can_be_recycled;

        /// Has the hash been computed
        bool has_hash;
#endif

#ifdef ALLOW_DEPRECATED_0_17_8
        edge_type the_edge_type;
#else
        const edge_type the_edge_type;
#endif

        /// True iff this was expanded from an identity-reduced edge
        bool orig_was_identity;

};

#endif

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
    class unpacked_node;
    class forest;
    class initializer_list;
    class node_marker;

    struct unpacked_lists;
}

#define ALLOW_SET_FROM_DDEDGE
#define ALLOW_SET_FROM_POINTER
// #define DEBUG_UNPACKED_HASH

#define USE_STRUCT

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
    public:

        /** Constructor.
            The class must be "filled" by a forest before
            it can be used, however.
                @param  f   Parent forest.
                            If null, we must call attach()
                            before use.
        */
        unpacked_node(const forest* f = nullptr);

        /// Destructor.
        ~unpacked_node();

        /// Attach the parent forest.
        void attach(const forest* f);

        /// Are we attached to f
        inline bool isAttachedTo(const forest* f) {
            return f == parent;
        }

    public:
        /* Initialization methods, primarily for reading */

        //
        // Build redundant nodes
        //

        void initRedundant(const forest *f, int k,
                node_handle node, node_storage_flags fs);

        void initRedundant(const forest *f, int k, const edge_value &ev,
                node_handle node, node_storage_flags fs);

        template <class T>
        inline void initRedundant(const forest *f, int k, T _ev, node_handle node,
                node_storage_flags fs)
        {
            initRedundant(f, k, edge_value(_ev), node, fs);
        }


        //
        // Build identity nodes
        //

        void initIdentity(const forest *f, int k, unsigned i,
                node_handle node, node_storage_flags fs);

        void initIdentity(const forest *f, int k, unsigned i,
                const edge_value &ev, node_handle node, node_storage_flags fs);

        template <class T>
        inline void initIdentity(const forest *f, int k, unsigned i, T _ev,
                node_handle node, node_storage_flags fs)
        {
            initIdentity(f, k, i, edge_value(_ev), node, fs);
        }


    public:
        //
        // For convenience: get recycled instance and initialize
        //

        static inline unpacked_node* newRedundant(const forest *f, int k,
                node_handle node, node_storage_flags fs)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->initRedundant(f, k, node, fs);
            return U;
        }

        static inline unpacked_node* newRedundant(const forest *f, int k,
                const edge_value &ev, node_handle node, node_storage_flags fs)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            if (ev.isVoid()) {
                U->initRedundant(f, k, node, fs);
            } else {
                U->initRedundant(f, k, ev, node, fs);
            }
            return U;
        }

        template <class T>
        static inline unpacked_node* newRedundant(const forest *f, int k,
                T _ev, node_handle node, node_storage_flags fs)
        {
            return newRedundant(f, k, edge_value(_ev), node, fs);
        }



        static inline unpacked_node* newIdentity(const forest *f, int k,
                unsigned i, node_handle node, node_storage_flags fs)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->initIdentity(f, k, i, node, fs);
            return U;
        }

        static inline unpacked_node* newIdentity(const forest *f, int k,
                unsigned i, const edge_value &ev, node_handle node,
                node_storage_flags fs)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            if (ev.isVoid()) {
                U->initIdentity(f, k, i, node, fs);
            } else {
                U->initIdentity(f, k, i, ev, node, fs);
            }
            return U;
        }

        template <class T>
        static inline unpacked_node* newIdentity(const forest *f, int k,
                unsigned i, T _ev, node_handle node, node_storage_flags fs)
        {
            return newIdentity(f, k, i, edge_value(_ev), node, fs);
        }

        /** Create a zeroed-out full node of a given size */
        template <class T>
        static inline unpacked_node* newFull(forest *f,
                int levl, T tsz)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            MEDDLY_DCASSERT(tsz >= 0);
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
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->level = levl;
            U->resize(unsigned(nnzs));
            U->setSparse();
            U->allowWrites(f);
            return U;
        }

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
            MEDDLY_DCASSERT(extra_unhashed);
            return extra_unhashed;
        }

        /// Set the unhashed header data
        inline void setUHdata(const void* p)
        {
            MEDDLY_DCASSERT(p);
            memcpy(extra_unhashed, p, extra_unhashed_size);
        }

        /// Get the unhashed header data
        inline void getUHdata(void* p) const
        {
            MEDDLY_DCASSERT(p);
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
            MEDDLY_DCASSERT(extra_hashed);
            return extra_hashed;
        }

        /// Set the hashed header data
        inline void setHHdata(const void* p)
        {
            MEDDLY_DCASSERT(p);
            memcpy(extra_hashed, p, extra_hashed_size);
        }

        /// Get the hashed header data
        inline void getHHdata(void* p) const
        {
            MEDDLY_DCASSERT(p);
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
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            return _idev[n].down;
#else
            MEDDLY_DCASSERT(_down);
            return _down[n];
#endif
        }

        /** Get a downward pointer.
            @param  n   Which pointer.
            @return     If this is a full node,
                        return pointer with index n.
                        If this is a sparse node,
                        return the nth non-zero pointer.
        */
        inline node_handle down(int n) const
        {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            return _idev[n].down;
#else
            MEDDLY_DCASSERT(_down);
            return _down[n];
#endif
        }

        /** Get the index of the nth non-zero pointer.
            Use only for sparse nodes.
            @param  n   Which pointer
            @return     The index of the pointer
        */
        inline unsigned index(unsigned n) const
        {
            MEDDLY_DCASSERT(!is_full);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            return _idev[n].index;
#else
            MEDDLY_DCASSERT(_index);
            return _index[n];
#endif
        }

        /** Get the index of the nth non-zero pointer.
            Use only for sparse nodes.
            @param  n   Which pointer
            @return     The index of the pointer
        */
        inline unsigned index(int n) const
        {
            MEDDLY_DCASSERT(!is_full);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            return _idev[n].index;
#else
            MEDDLY_DCASSERT(_index);
            return _index[n];
#endif
        }

        /** Get the nth edge value.
            @param  n   Which pointer
            @return     The edge value
        */
        inline const edge_value& edgeval(unsigned n) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            return _idev[n].edgeval;
#else
            MEDDLY_DCASSERT(_edge);
            return _edge[n];
#endif
        }

        /** Get the nth edge value.
            @param  n   Which pointer
            @return     The edge value
        */
        inline const edge_value& edgeval(int n) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            return _idev[n].edgeval;
#else
            MEDDLY_DCASSERT(_edge);
            return _edge[n];
#endif
        }


        /** Subtract from an edge value.
            @param  n       Which pointer
            @param  v       Value to subtract
        */
        template <class T>
        inline void subtractFromEdge(unsigned n, T v) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            _idev[n].edgeval.subtract(v);
#else
            MEDDLY_DCASSERT(_edge);
            _edge[n].subtract(v);
#endif
        }


        /** Divide an edge value.
            @param  n       Which pointer
            @param  v       Value to divide by
        */
        template <class T>
        inline void divideEdge(unsigned n, T v) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            _idev[n].edgeval.divide(v);
#else
            MEDDLY_DCASSERT(_edge);
            _edge[n].divide(v);
#endif
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
            MEDDLY_DCASSERT(ev.hasType(the_edge_type));
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            _idev[n].edgeval = ev;
#else
            MEDDLY_DCASSERT(_edge);
            _edge[n] = ev;
#endif
        }

        /**
            Set a full edge.
                @param  n       Which pointer
                @param  h       Node handle
        */
        inline void setFull(unsigned n, node_handle h)
        {
            MEDDLY_DCASSERT(!hasEdges());
            MEDDLY_DCASSERT(isFull());
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            _idev[n].down = h;
#else
            MEDDLY_DCASSERT(_down);
            _down[n] = h;
#endif
        }

        /**
            Set a full edge.
                @param  n       Which pointer
                @param  v       Edge value
                @param  h       Node handle
        */
        inline void setFull(unsigned n, const edge_value &v, node_handle h)
        {
            MEDDLY_DCASSERT(v.hasType(the_edge_type));
            MEDDLY_DCASSERT(isFull());
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            _idev[n].down = h;
            _idev[n].edgeval = v;
#else
            MEDDLY_DCASSERT(_down);
            _down[n] = h;
            if (_edge) {
                _edge[n] = v;
            }
#endif
        }

#ifdef ALLOW_SET_FROM_DDEDGE
        /**
            Set a full edge from E, and destroy E.
                @param  n       Which pointer
                @param  E       Edge value (if needed) and node
        */
        inline void setFull(unsigned n, dd_edge &E)
        {
            MEDDLY_DCASSERT(E.isAttachedTo(parent));
            MEDDLY_DCASSERT(E.getEdgeValue().hasType(the_edge_type));
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            E.xferNode(_idev[n].down);
            _idev[n].edgeval = E.getEdgeValue();
#else
            E.xferNode(_down[n]);
            if (_edge) {
                _edge[n] = E.getEdgeValue();
            }
#endif
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
            MEDDLY_DCASSERT(isFull());
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            _idev[n].down = h;
            _idev[n].edgeval.set(the_edge_type, p);
#else
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_edge);
            _down[n] = h;
            _edge[n].set(the_edge_type, p);
#endif
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
            MEDDLY_DCASSERT(!hasEdges());
            MEDDLY_DCASSERT(isSparse());
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            _idev[n].index = i;
            _idev[n].down = h;
#else
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_index);
            _index[n] = i;
            _down[n] = h;
#endif
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
            MEDDLY_DCASSERT(v.hasType(the_edge_type));
            MEDDLY_DCASSERT(isSparse());
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            _idev[n].index = i;
            _idev[n].down = h;
            _idev[n].edgeval = v;
#else
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_index);
            _index[n] = i;
            _down[n] = h;
            if (_edge) {
                _edge[n] = v;
            }
#endif
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
            MEDDLY_DCASSERT(E.getEdgeValue().hasType(the_edge_type));
            MEDDLY_DCASSERT(isSparse());
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            _idev[n].index = i;
            E.xferNode(_idev[n].down);
            _idev[n].edgeval = E.getEdgeValue();
#else
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_index);
            E.xferNode(_down[n]);
            _index[n] = i;
            if (_edge) {
                _edge[n] = E.getEdgeValue();
            }
#endif
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
            MEDDLY_DCASSERT(isSparse());
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
#ifdef USE_STRUCT
            MEDDLY_DCASSERT(_idev);
            _idev[n].index = i;
            _idev[n].down = h;
            _idev[n].edgeval.set(the_edge_type, p);
#else
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_index);
            MEDDLY_DCASSERT(_edge);
            _index[n] = i;
            _down[n] = h;
            _edge[n].set(the_edge_type, p);
#endif
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
            MEDDLY_DCASSERT(can_be_extensible);
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
            MEDDLY_DCASSERT(isExtensible());
            return down(size-1);
        }

        /// Get the extensible index
        inline unsigned ext_i() const
        {
            MEDDLY_DCASSERT(isExtensible());
            return index(size - 1);
        }

        /// Get the extensible edge value
        inline const edge_value& ext_ev() const
        {
            MEDDLY_DCASSERT(isExtensible());
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
            MEDDLY_DCASSERT(has_hash);
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
            MEDDLY_DCASSERT(ns>=0);
            resize(unsigned(ns));
        }

        /// Shrink the size of a node
        inline void shrink(unsigned ns)
        {
            MEDDLY_DCASSERT(ns <= getSize());
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
            MEDDLY_DCASSERT(nullptr == modparent);
            MEDDLY_DCASSERT(parent == mp);
            modparent = mp;
            AddToBuildList(this);
            mark_extra = 0;
        }

        // Set edges to transparent
        void clear(unsigned low, unsigned high);

    public:
        /// Was this node initialized from a redundant node
        inline bool wasRedundant() const {
            return orig_was_fully;
        }
        /// Was this node initialized from an identity node
        inline bool wasIdentity() const {
            return orig_was_identity;
        }

    public:
        //
        // Centralized per-forest (using FIDs) recycling
        //

        /// Pull a recycled node off of f's free list,
        /// or create a new one if needed.
        static unpacked_node* New(const forest* f);

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
        inline void setRegular() {
            orig_was_fully = orig_was_identity = false;
        }
        inline void setRedundant() {
            orig_was_fully = true;
            orig_was_identity = false;
        }
        inline void setIdentity() {
            orig_was_fully = false;
            orig_was_identity = true;
        }

#ifdef USE_STRUCT
    private:
        struct edgeinfo {
            unsigned index;
            node_handle down;
            edge_value edgeval;
        };
#endif

    private:
        /// Next in list
        unpacked_node* next;
        /// Previous in list (for build lists only)
        unpacked_node* prev;

        /// Forest where the node belongs
        const forest* parent;
        /// Modifiable parent forest; required for writable nodes
        forest* modparent;

#ifdef USE_STRUCT
        /// Array of indexes, down pointers, and edge values.
        edgeinfo* _idev;
#else
        /// Down pointers
        node_handle* _down;

        /// Indexes, for sparse; otherwise unused
        unsigned* _index;

        /// Edge values; or null if void edges
        edge_value* _edge;
#endif

        /// Extra, temporary, node handle for marking.
        node_handle mark_extra;

        /// Allocated sizes of arrays
        unsigned alloc;

        /// Used sizes of arrays
        unsigned size;

        /// Extra header information that is not hashed
        void* extra_unhashed;
        /// Number of bytes in extra unhashed header
        unsigned extra_unhashed_size;

        /// Extra header information that is hashed
        void* extra_hashed;
        /// Number of bytes in extra hashed header
        unsigned extra_hashed_size;

        /// Level of the node
        int level;

        /// Hash of the node
        unsigned the_hash;

        /// FID of the parent
        unsigned pFID;

        /// Are we assuming full storage?
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

        edge_type the_edge_type;

        /// True iff this was expanded from a fully-reduced edge
        bool orig_was_fully;

        /// True iff this was expanded from an identity-reduced edge
        bool orig_was_identity;

};


#endif

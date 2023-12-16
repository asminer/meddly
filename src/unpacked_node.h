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

namespace MEDDLY {
    class unpacked_node;
    class forest;
    class initializer_list;
    class node_marker;

    struct unpacked_lists;
}

// #define REMOVE_OLD

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
            U->initRedundant(f, k, ev, node, fs);
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
            U->initIdentity(f, k, i, ev, node, fs);
            return U;
        }

        template <class T>
        static inline unpacked_node* newIdentity(const forest *f, int k,
                unsigned i, T _ev, node_handle node, node_storage_flags fs)
        {
            return newIdentity(f, k, i, edge_value(_ev), node, fs);
        }

        /** Create a zeroed-out full node */
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
            MEDDLY_DCASSERT(_down);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            return _down[n];
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
            MEDDLY_DCASSERT(_down);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
            return _down[n];
        }

        /** Get the index of the nth non-zero pointer.
            Use only for sparse nodes.
            @param  n   Which pointer
            @return     The index of the pointer
        */
        inline unsigned index(unsigned n) const
        {
            MEDDLY_DCASSERT(_index);
            MEDDLY_DCASSERT(!is_full);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            return _index[n];
        }

        /** Get the index of the nth non-zero pointer.
            Use only for sparse nodes.
            @param  n   Which pointer
            @return     The index of the pointer
        */
        inline unsigned index(int n) const
        {
            MEDDLY_DCASSERT(_index);
            MEDDLY_DCASSERT(!is_full);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
            return _index[n];
        }

        /** Get the nth edge value.
            @param  n   Which pointer
            @return     The edge value
        */
        inline const edge_value& edgeval(unsigned n) const {
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            return _edge[n];
        }

        /** Get the nth edge value.
            @param  n   Which pointer
            @return     The edge value
        */
        inline const edge_value& edgeval(int n) const {
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, int(size));
            return _edge[n];
        }

        /** Subtract from an edge value.
            @param  n       Which pointer
            @param  v       Value to subtract
        */
        template <class T>
        inline void subtractFromEdge(unsigned n, T v) {
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            _edge[n].subtract(v);
        }


        /** Divide an edge value.
            @param  n       Which pointer
            @param  v       Value to divide by
        */
        template <class T>
        inline void divideEdge(unsigned n, T v) {
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
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
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            MEDDLY_DCASSERT(ev.hasType(the_edge_type));
            _edge[n] = ev;
        }

        /**
            Set a full edge.
                @param  n       Which pointer
                @param  h       Node handle
        */
        inline void setFull(unsigned n, node_handle h)
        {
            MEDDLY_DCASSERT(_down);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            MEDDLY_DCASSERT(!hasEdges());
            MEDDLY_DCASSERT(isFull());

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
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            MEDDLY_DCASSERT(v.hasType(the_edge_type));
            MEDDLY_DCASSERT(isFull());

            _down[n] = h;
            _edge[n] = v;
        }

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
            E.xferNode(_down[n]);
            if (_edge) {
                _edge[n] = E.getEdgeValue();
            }
        }


        /**
            Set a full edge.
                @param  n       Which pointer
                @param  p       Raw pointer for edge value
                @param  h       Node handle
        */
        inline void setFull(unsigned n, const void* p, node_handle h)
        {
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            MEDDLY_DCASSERT(isFull());

            _down[n] = h;
            _edge[n].set(the_edge_type, p);
        }


        /**
            Set a sparse edge.
                @param  n       Which nonzero edge
                @param  i       Index of the edge
                @param  h       Node handle
        */
        inline void setSparse(unsigned n, unsigned i, node_handle h)
        {
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_index);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            MEDDLY_DCASSERT(!hasEdges());
            MEDDLY_DCASSERT(isSparse());

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
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_index);
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            MEDDLY_DCASSERT(v.hasType(the_edge_type));
            MEDDLY_DCASSERT(isSparse());

            _index[n] = i;
            _down[n] = h;
            _edge[n] = v;
        }


        /**
            Set a sparse edge.
                @param  n       Which nonzero edge
                @param  i       Index of the edge
                @param  E       Edge value (if needed) and node
        */
        inline void setSparse(unsigned n, unsigned i, dd_edge &E)
        {
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_index);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            MEDDLY_DCASSERT(E.getEdgeValue().hasType(the_edge_type));
            MEDDLY_DCASSERT(isSparse());

            E.xferNode(_down[n]);
            _index[n] = i;
            if (_edge) {
                _edge[n] = E.getEdgeValue();
            }
        }


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
            MEDDLY_DCASSERT(_down);
            MEDDLY_DCASSERT(_index);
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            MEDDLY_DCASSERT(isSparse());

            _index[n] = i;
            _down[n] = h;
            _edge[n].set(the_edge_type, p);
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
            return is_extensible;
        }

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

        /// Checks if the node has no trailing redundant edges.
        /// I.e., the node does NOT have trailing edges that can
        /// be collapsed into the extensible edge.
        inline bool isTrim() const
        {
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
            MEDDLY_DCASSERT(ns <= size);
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
        }

        // Set edges to transparent
        void clear(unsigned low, unsigned high);

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

        /// Down pointers
        node_handle* _down;

        /// Indexes, for sparse; otherwise unused
        unsigned* _index;

        /// Edge values; or null if void edges
        edge_value* _edge;

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

        /// Can this node be extensible?
        /// Determined by the parent forest.
        bool can_be_extensible;

        /// Is this node extensible?
        bool is_extensible;

#ifdef DEVELOPMENT_CODE
        /// only node pointers built by New() should be Recycle()d.
        bool can_be_recycled;

        /// Has the hash been computed
        bool has_hash;
#endif

        edge_type the_edge_type;


/*
    OLD interface below here
*/

#ifdef ALLOW_DEPRECATED_0_17_4

public:
        inline void initRedundant(const forest *f, int k,
                node_handle node, bool full)
        {
            initRedundant(f, k, node, full ? FULL_ONLY : SPARSE_ONLY);
        }

        inline void initIdentity(const forest *f, int k, unsigned i,
                node_handle node, bool full)
        {
            initIdentity(f, k, i, node, full ? FULL_ONLY : SPARSE_ONLY);
        }

#ifdef REMOVE_OLD
        template <class T>
        inline void initRedundant(const forest *f, int k, T _ev, node_handle node,
                bool full)
        {
            edge_value ev(_ev);
            initRedundant(f, k, ev, node, full ? FULL_ONLY : SPARSE_ONLY);
        }

        template <class T>
        inline void initIdentity(const forest *f, int k, unsigned i, T _ev,
                node_handle node, bool full)
        {
            edge_value ev(_ev);
            initIdentity(f, k, i, ev, node, full ? FULL_ONLY : SPARSE_ONLY);
        }
#else
        void initRedundant(const forest *f, int k, int ev, node_handle node,
                bool full);

        void initRedundant(const forest *f, int k, long ev, node_handle node,
                bool full);

        void initRedundant(const forest *f, int k, float ev, node_handle node,
                bool full);

        void initIdentity(const forest *f, int k, unsigned i,
                int ev, node_handle node, bool full);

        void initIdentity(const forest *f, int k, unsigned i,
                long ev, node_handle node, bool full);

        void initIdentity(const forest *f, int k, unsigned i,
                float ev, node_handle node, bool full);
#endif

        static inline unpacked_node* newRedundant(const forest *f,
                int k, node_handle node, bool full)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->initRedundant(f, k, node, full);
            return U;
        }

        static inline unpacked_node* newRedundant(const forest *f,
                int k, long ev, node_handle node, bool full)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->initRedundant(f, k, ev, node, full);
            return U;
        }

        static inline unpacked_node* newRedundant(const forest *f,
                int k, float ev, node_handle node, bool full)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->initRedundant(f, k, ev, node, full);
            return U;
        }

        static inline unpacked_node* newIdentity(const forest *f,
                int k, unsigned i, node_handle node, bool full)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->initIdentity(f, k, i, node, full);
            return U;
        }

        static inline unpacked_node* newIdentity(const forest *f,
                int k, unsigned i, long ev, node_handle node, bool full)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->initIdentity(f, k, i, ev, node, full);
            return U;
        }

        static inline unpacked_node* newIdentity(const forest *f,
                int k, unsigned i, float ev, node_handle node, bool full)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->initIdentity(f, k, i, ev, node, full);
            return U;
        }

        /// Modify a pointer to the unhashed header data
        inline void* UHdata()
        {
            MEDDLY_DCASSERT(extra_unhashed);
            return extra_unhashed;
        }

        /// Modify a pointer to the hashed header data.
        inline void* HHdata()
        {
            MEDDLY_DCASSERT(extra_hashed);
            return extra_hashed;
        }


        inline node_handle d(unsigned n) const
        {
            return down(n);
        }

        /** Reference to a downward pointer.
            @param  n   Which pointer.
            @return     If this is a full node,
                        modify pointer with index n.
                        If this is a sparse node,
                        modify the nth non-zero pointer.
        */
        inline node_handle& d_ref(unsigned n)
        {
            MEDDLY_DCASSERT(_down);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            return _down[n];
        }

        /// Set the nth pointer from E, and destroy E.
        inline void set_d(unsigned n, dd_edge &E)
        {
            MEDDLY_DCASSERT(E.isAttachedTo(parent));
            MEDDLY_DCASSERT(!hasEdges());
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            E.xferNode(_down[n]);
        }

        inline unsigned i(unsigned n) const {
            return index(n);
        }

        /** Modify the index of the nth non-zero pointer.
            Use only for sparse nodes.
        */
        inline unsigned& i_ref(unsigned n)
        {
            MEDDLY_DCASSERT(_index);
            MEDDLY_DCASSERT(!is_full);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            return _index[n];
        }

        /// Set the nth pointer and edge value from E, and destroy E.
        inline void set_de(unsigned n, dd_edge &E)
        {
            MEDDLY_DCASSERT(E.isAttachedTo((forest*) parent));
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            E.xferNode(_down[n]);
            _edge[n] = E.getEdgeValue();
        }

        /// Get the edge value
        template <class T>
        inline void getEdge(unsigned i, T& val) const
        {
            edgeval(i).get(val);
        }

        /// Get the edge value, generic
        inline const edge_value& getEdge(unsigned i) const
        {
            return edgeval(i);
        }

        /** Set the nth edge value.
            @param  n   Which pointer
            @param  ev  Edge value
        */
        inline void set_edgeval(unsigned n, const edge_value &v) {
            MEDDLY_DCASSERT(parent);
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            MEDDLY_DCASSERT(v.hasType(the_edge_type));
            _edge[n] = v;
        }

        /** Set the nth edge value from a pointer.
            @param  n   Which edge
            @param  p   Pointer to edge value
        */
        inline void setEdgeRaw(unsigned n, const void* p) {
            MEDDLY_DCASSERT(parent);
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, size);
            _edge[n].set(the_edge_type, p);
        }

        /// Set the edge value
        template <class T>
        inline void setEdge(unsigned i, T val)
        {
            edge_value ev(val);
            set_edgeval(i, ev);
        }

        // Set the edge value, generic
        inline void setEdge(unsigned i, const edge_value &ev) {
            set_edgeval(i, ev);
        }

        /// Get the edge value, as an integer.
        inline long ei(unsigned i) const
        {
            return edge_long(i);
        }

        /// Get the edge value, as a float.
        inline float ef(unsigned i) const
        {
            return edge_float(i);
        }


        /// Get the extensible edge value, as a long
        inline long ext_ei() const
        {
            MEDDLY_DCASSERT(isExtensible());
            return ei(size - 1);
        }

        /// Get the extensible edge value, as a float
        inline float ext_ef() const
        {
            MEDDLY_DCASSERT(isExtensible());
            return ef(size - 1);
        }

        /// Get the number of nonzeroes of this node (sparse nodes only).
        inline unsigned getNNZs() const
        {
            return size;
        }


        /// Shrink the size of a (truncated) full node
        inline void shrinkFull(unsigned ns)
        {
            shrink(ns);
        }


        /// Shrink the size of a sparse node
        inline void shrinkSparse(unsigned ns)
        {
            shrink(ns);
        }


        /// Called by node_storage when building an unpacked
        /// node based on how it's stored.
        inline void bind_as_full(bool full)
        {
            is_full = full;
        }

#endif








};


#endif

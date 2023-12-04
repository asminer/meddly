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


    public:
        /* Initialization methods, primarily for reading */

        //
        // Build redundant nodes
        //

        void initRedundant(const forest *f, int k, const edge_value &ev,
                node_handle node, node_storage_flags fs);

#ifdef ALLOW_DEPRECATED_0_17_4
#ifdef REMOVE_OLD
        inline void initRedundant(const forest *f, int k, node_handle node,
                bool full)
        {
            edge_value ev;
            initRedundant(f, k, ev, node, full ? FULL_ONLY : SPARSE_ONLY);
        }

        template <class T>
        inline void initRedundant(const forest *f, int k, T _ev, node_handle node,
                bool full)
        {
            edge_value ev(_ev);
            initRedundant(f, k, ev, node, full ? FULL_ONLY : SPARSE_ONLY);
        }
#else
        void initRedundant(const forest *f, int k, node_handle node,
                bool full);

        void initRedundant(const forest *f, int k, int ev, node_handle node,
                bool full);

        void initRedundant(const forest *f, int k, long ev, node_handle node,
                bool full);

        void initRedundant(const forest *f, int k, float ev, node_handle node,
                bool full);
#endif
#endif

        //
        // Build identity nodes
        //

        void initIdentity(const forest *f, int k, unsigned i,
                const edge_value &ev, node_handle node, node_storage_flags fs);

#ifdef ALLOW_DEPRECATED_0_17_4
#ifdef REMOVE_OLD
        inline void initIdentity(const forest *f, int k, unsigned i,
                node_handle node, bool full)
        {
            edge_value ev;
            initIdentity(f, k, i, ev, node, full ? FULL_ONLY : SPARSE_ONLY);
        }

        template <class T>
        inline void initIdentity(const forest *f, int k, unsigned i, T _ev,
                node_handle node, bool full)
        {
            edge_value ev(_ev);
            initIdentity(f, k, i, ev, node, full ? FULL_ONLY : SPARSE_ONLY);
        }
#else
        void initIdentity(const forest *f, int k, unsigned i,
                node_handle node, bool full);

        void initIdentity(const forest *f, int k, unsigned i,
                int ev, node_handle node, bool full);

        void initIdentity(const forest *f, int k, unsigned i,
                long ev, node_handle node, bool full);

        void initIdentity(const forest *f, int k, unsigned i,
                float ev, node_handle node, bool full);
#endif
#endif

        //
        // Create blank nodes, primarily for writing
        //


/*
        inline void initFull(const forest *f, int levl, unsigned tsz)
        {
            MEDDLY_DCASSERT(f);
            bind_to_forest(f, levl, tsz, true);
        }

        inline void initSparse(const forest *f, int levl, unsigned nnz)
        {
            MEDDLY_DCASSERT(f);
            bind_to_forest(f, levl, nnz, false);
        }
*/

    public:
        //
        // For convenience: get recycled instance and initialize
        //

        static inline unpacked_node* newRedundant(const forest *f, int k,
                const edge_value &ev, node_handle node, node_storage_flags fs)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->initRedundant(f, k, ev, node, fs);
            return U;
        }

#ifdef ALLOW_DEPRECATED_0_17_4
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
#endif

        static inline unpacked_node* newIdentity(const forest *f, int k,
                unsigned i, const edge_value &ev, node_handle node,
                node_storage_flags fs)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->initIdentity(f, k, i, ev, node, fs);
            return U;
        }

#ifdef ALLOW_DEPRECATED_0_17_4
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
#endif

        /** Create a zeroed-out full node */
        static inline unpacked_node* newFull(forest *f,
                int levl, unsigned tsz)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->level = levl;
            U->resize(tsz);
            U->setFull();
            U->allowWrites(f);
            return U;
        }

        /** Create a zeroed-out sparse node */
        static inline unpacked_node* newSparse(forest *f,
                int levl, unsigned nnzs)
        {
            unpacked_node* U = New(f);
            MEDDLY_DCASSERT(U);
            U->level = levl;
            U->resize(nnzs);
            U->setSparse();
            U->allowWrites(f);
            return U;
        }

    public:
        //
        // Display methods
        //

// HERE <<<------------

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

        /// Modify a pointer to the unhashed header data
        inline void* UHdata()
        {
            MEDDLY_DCASSERT(extra_unhashed);
            return extra_unhashed;
        }

        /// Get a pointer to the unhashed header data.
        inline const void* UHptr() const
        {
            MEDDLY_DCASSERT(extra_unhashed);
            return extra_unhashed;
        }

        /// Get the number of bytes of unhashed header data.
        inline unsigned UHbytes() const
        {
            return ext_uh_size;
        }

        /// Get a pointer to the hashed header data.
        inline const void* HHptr() const
        {
            MEDDLY_DCASSERT(extra_hashed);
            return extra_hashed;
        }

        /// Modify a pointer to the hashed header data.
        inline void* HHdata()
        {
            MEDDLY_DCASSERT(extra_hashed);
            return extra_hashed;
        }

        /// Get the number of bytes of hashed header data.
        inline unsigned HHbytes() const
        {
            return ext_h_size;
        }

        /** Get a downward pointer.
            @param  n   Which pointer.
            @return     If this is a full reader,
                        return pointer with index n.
                        If this is a sparse reader,
                        return the nth non-zero pointer.
        */
        inline node_handle d(unsigned n) const
        {
            MEDDLY_DCASSERT(_down);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n,
                    is_full ? size : nnzs);
            return _down[n];
        }

        /** Reference to a downward pointer.
            @param  n   Which pointer.
            @return     If this is a full reader,
                        modify pointer with index n.
                        If this is a sparse reader,
                        modify the nth non-zero pointer.
        */
        inline node_handle& d_ref(unsigned n)
        {
            MEDDLY_DCASSERT(_down);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n,
                    is_full ? size : nnzs);
            return _down[n];
        }

        /// Set the nth pointer from E, and destroy E.
        inline void set_d(unsigned n, dd_edge &E)
        {
            MEDDLY_DCASSERT(E.isAttachedTo((forest*) parent));
            MEDDLY_DCASSERT(edge_type::VOID == the_edge_type);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n,
                    is_full ? size : nnzs);
            _down[n] = E.node;
            E.node = 0; // avoid having to adjust the link count
        }

        /** Get the index of the nth non-zero pointer.
            Use only for sparse readers.
        */
        inline unsigned i(unsigned n) const
        {
            MEDDLY_DCASSERT(_index);
            MEDDLY_DCASSERT(!is_full);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, nnzs);
            return _index[n];
        }

        /** Modify the index of the nth non-zero pointer.
            Use only for sparse readers.
        */
        inline unsigned& i_ref(unsigned n)
        {
            MEDDLY_DCASSERT(_index);
            MEDDLY_DCASSERT(!is_full);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, nnzs);
            return _index[n];
        }

        /*
        /// Get a pointer to an edge
        inline const void* eptr(unsigned i) const
        {
            MEDDLY_DCASSERT(edge);
            MEDDLY_DCASSERT(edge_bytes);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i,
                    is_full ? size : nnzs);
            return ((char*) edge) + i * edge_bytes;
        }

        /// Modify pointer to an edge
        inline void* eptr_write(unsigned i)
        {
            MEDDLY_DCASSERT(edge);
            MEDDLY_DCASSERT(edge_bytes);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i,
                    is_full ? size : nnzs);
            return ((char*) edge) + i * edge_bytes;
        }
        */

        /// Set the nth pointer and edge value from E, and destroy E.
        inline void set_de(unsigned n, dd_edge &E)
        {
            MEDDLY_DCASSERT(E.isAttachedTo((forest*) parent));
            MEDDLY_DCASSERT(_edge);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n,
                    is_full ? size : nnzs);
            _down[n] = E.node;
            _edge[n] = E.getEdgeValue();
            E.node = 0; // avoid having to adjust the link count
        }

        /// Get the edge value, as an integer.
        inline void getEdge(unsigned i, long& val) const
        {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i,
                    is_full ? size : nnzs);
            MEDDLY_DCASSERT(edge_type::LONG == the_edge_type);
            _edge[i].get(val);
        }

        /// Get the edge value, as a float.
        inline void getEdge(unsigned i, float& val) const
        {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i,
                    is_full ? size : nnzs);
            MEDDLY_DCASSERT(edge_type::FLOAT == the_edge_type);
            _edge[i].get(val);
        }

        /// Get the edge value, generic
        inline const edge_value& getEdge(unsigned i) const
        {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i,
                    is_full ? size : nnzs);
            return _edge[i];
        }

        /// Set the edge value, as an integer.
        inline void setEdge(unsigned i, long ev)
        {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i,
                    is_full ? size : nnzs);
            MEDDLY_DCASSERT(edge_type::LONG == the_edge_type);
            _edge[i].set(ev);
        }

        /// Set the edge value, as a float.
        inline void setEdge(unsigned i, float ev)
        {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i,
                    is_full ? size : nnzs);
            MEDDLY_DCASSERT(edge_type::FLOAT == the_edge_type);
            _edge[i].set(ev);
        }

        // Set the edge value, generic
        inline void setEdge(unsigned i, const edge_value &ev) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i,
                    is_full ? size : nnzs);
            MEDDLY_DCASSERT(ev.hasType(the_edge_type));
            _edge[i] = ev;
        }

        // Set the edge value, generic
        inline edge_value& setEdge(unsigned i) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i,
                    is_full ? size : nnzs);
            return _edge[i];
        }


        /// Get the edge value, as an integer.
        inline long ei(unsigned i) const
        {
            long ev;
            getEdge(i, ev);
            return ev;
        }

        /// Get the edge value, as a float.
        inline float ef(unsigned i) const
        {
            float ev;
            getEdge(i, ev);
            return ev;
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
            return d( (isSparse()? getNNZs(): getSize()) - 1 );
        }

        /// Get the extensible index
        inline unsigned ext_i() const
        {
            MEDDLY_DCASSERT(isExtensible());
            return isSparse()? i(getNNZs() - 1): getSize() - 1;
        }

        /// Get the extensible edge value, as a long
        inline long ext_ei() const
        {
            MEDDLY_DCASSERT(isExtensible());
            return ei( (isSparse()? getNNZs(): getSize()) - 1 );
        }

        /// Get the extensible edge value, as a float
        inline float ext_ef() const
        {
            MEDDLY_DCASSERT(isExtensible());
            return ef( (isSparse()? getNNZs(): getSize()) - 1 );
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

#ifdef ALLOW_DEPRECATED_0_17_4
        /// Get the number of nonzeroes of this node (sparse readers only).
        inline unsigned getNNZs() const
        {
            return nnzs;
        }
#endif

        /// Is this a sparse reader?
        inline bool isSparse() const
        {
            return !is_full;
        }

        /// Is this a full reader?
        inline bool isFull() const
        {
            return is_full;
        }

        /// Does this node have edge values?
        inline bool hasEdges() const
        {
            return edge_type::VOID != parent->getEdgeType();
        }

        /// Edge type
        inline edge_type getEdgeType() const {
            return parent->getEdgeType();
        }

    public:
        /// Get the node's hash
        inline unsigned hash() const
        {
#ifdef DEVELOPMENT_CODE
            MEDDLY_DCASSERT(has_hash);
#endif
            return h;
        }


        /// Set the node's hash
        inline void setHash(unsigned H)
        {
#ifdef DEVELOPMENT_CODE
            MEDDLY_DCASSERT(!has_hash);
            has_hash = true;
#endif
            h = H;
        }

        /// Compute the node's hash
        void computeHash();

        /// Removes redundant trailing edges.
        /// If the unpacked node is sparse, it assumes its indices
        /// to be in ascending order.
        void trim();

        /// If the unpacked node is sparse, it is sorted so that the
        /// indices are in ascending order.
        void sort();

        /// Checks if the node is has no trailing redundant edges
        bool isTrim() const;

        // checks if the node indices are in ascending order
        bool isSorted() const;

        // Is this a "build" node?  Important for mark and sweep.
        inline bool isBuildNode() const
        {
            return is_in_build_list;
        }


    public:

        /// Change the size of a node
        inline void resize(unsigned ns) {
            size = ns;
            if (ns > alloc) enlarge(ns);
        }

        /// Shrink the size of a node
        inline void shrink(unsigned ns)
        {
            MEDDLY_DCASSERT(ns <= size);
            size = ns;
        }

#ifdef ALLOW_DEPRECATED_0_17_4
        /// Shrink the size of a (truncated) full node
        inline void shrinkFull(unsigned ns)
        {
            MEDDLY_DCASSERT(isFull());
            MEDDLY_DCASSERT(ns >= 0);
            MEDDLY_DCASSERT(ns <= size);
            size = ns;
        }


        /// Shrink the size of a sparse node
        inline void shrinkSparse(unsigned ns)
        {
            MEDDLY_DCASSERT(isSparse());
            MEDDLY_DCASSERT(ns >= 0);
            MEDDLY_DCASSERT(ns <= nnzs);
            nnzs = ns;
        }
#endif

        /// Called within forest to allocate space.
        ///   @param  p     Parent, for building nodes.
        ///   @param  k     Level number.
        ///   @param  ns    Size of node.
        ///   @param  full  If true, we'll be filling a full reader.
        ///                 Otherwise it is a sparse one.
        // void bind_to_forest(const forest* p, int k, unsigned ns, bool full);

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

    public:
        /// Set down edges from 0..stop-1 to the transparent edge
        void clearEdges(unsigned stop);

        /// Set all down edges (and values) to 0, for full storage
        inline void clearFullEdges() {
            MEDDLY_DCASSERT(isFull());
            clearEdges(size);
        }

        /// Set all down edges (and values) to 0, for sparse storage
        inline void clearSparseEdges() {
            MEDDLY_DCASSERT(isSparse());
            clearEdges(nnzs);
        }


    public:
        //
        // Centralized per-forest (using FIDs) recycling
        //

        /// Pull a recycled node off of f's free list,
        /// or create a new one if needed.
        static unpacked_node* New(const forest* f);

        /// Add r to its forest's build list
        static void AddToBuildList(unpacked_node* r);

        /// Remove r from its forest's build list if needed,
        /// and add r back to its forest's recycle list.
        static void Recycle(unpacked_node* r);


/*
        // Pull a recycled node off the free list
        // static inline unpacked_node* useUnpackedNode()
        static inline unpacked_node* New(const forest* f)
        {
            unpacked_node* nr;
            if (freeList) {
                nr = freeList;
                freeList = nr->next;
            } else {
                nr = new unpacked_node;
            }
#ifdef DEVELOPMENT_CODE
            nr->has_hash = false;
#endif
            nr->is_in_build_list = false;
            return nr;
        }


        // Add a recycled node to the free list
        static inline void recycle(unpacked_node* r)
        {
            if (r) {
                if (r->is_in_build_list) {
                    removeFromBuildList(r);
                }
                r->next = freeList;
                freeList = r;
            }
        }

        // delete all nodes in the free list
        static void freeRecycled();

        // Lists so we can treat nodes under construction
        // as root nodes for mark and sweep

        // Add a node to the build list
        // Signifies that we are building this node, so we link/unlink
        // the down pointers.
        static inline void addToBuildList(unpacked_node* b, forest *mp)
        {
#ifdef DEBUG_BUILDLIST
            printf("Adding unpacked node at level %d to build list\n", b->getLevel());
#endif
            MEDDLY_DCASSERT(b);
            MEDDLY_DCASSERT(b->parent == mp);
            b->modparent = mp;

            b->is_in_build_list = true;
            b->next = buildList;
            buildList = b;
        }

        static void removeFromBuildList(unpacked_node* b);
        static void markBuildListChildren(node_marker* M);
*/


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

    private:
        /// Next in list
        unpacked_node* next;
        /// Previous in list
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
        unsigned h;

        /// only node pointers built by new() should be recycle()d.
        bool can_be_recycled;

        /// Are we assuming full storage?
        bool is_full;

        /// Can this node be extensible?
        /// Determined by the parent forest.
        bool can_be_extensible;

        /// Is this node extensible?
        bool is_extensible;

#ifdef DEVELOPMENT_CODE
        /// Has the hash been computed
        bool has_hash;
#endif

        // TBD: do we need these?

        edge_type the_edge_type;
        terminal_type the_terminal_type;


};


#endif

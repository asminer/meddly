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
#include "encoders.h"
#include "dd_edge.h"

namespace MEDDLY {
    class unpacked_node;
    class forest;
    class expert_forest;
    class initializer_list;
}

// ******************************************************************
// *                                                                *
// *                      unpacked_node  class                      *
// *                                                                *
// ******************************************************************

/** Class for reading nodes.
    Ideally - used anywhere we want to read node data.
    Backend implementation may change :^)
    Implemented in node_wrappers.cc.
    Readers may be "full" or "sparse",
    regardless of how the actual node is stored.

    Currently, a node is:
      array of node_handles, for the downward pointers
      array of integers, for indexes (sparse only)
      "compact" chunk of memory for edge values
*/
class MEDDLY::unpacked_node {
        friend class initializer_list;
    public:
        /** Constructor.
            The class must be "filled" by a forest before
            it can be used, however.

            TBD: make this private
        */
        unpacked_node();

        /// Destructor.  TBD: make this private
        ~unpacked_node();

        /// Free memory, but don't delete.
        void clear();

    public:
        /* Initialization methods, primarily for reading */

        //
        // Build redundant nodes
        //

        void initRedundant(const expert_forest *f, int k, node_handle node,
                bool full);

        void initRedundant(const expert_forest *f, int k, int ev,
                node_handle node, bool full);

        void initRedundant(const expert_forest *f, int k, long ev,
                node_handle node, bool full);

        void initRedundant(const expert_forest *f, int k, float ev,
                node_handle node, bool full);

        //
        // Build identity nodes
        //

        void initIdentity(const expert_forest *f, int k, unsigned i,
                node_handle node, bool full);

        void initIdentity(const expert_forest *f, int k, unsigned i,
                int ev, node_handle node, bool full);

        void initIdentity(const expert_forest *f, int k, unsigned i,
                long ev, node_handle node, bool full);

        void initIdentity(const expert_forest *f, int k, unsigned i,
                float ev, node_handle node, bool full);

        //
        // Create blank nodes, primarily for writing
        //

        inline void initFull(const expert_forest *f, int levl, unsigned tsz)
        {
            MEDDLY_DCASSERT(f);
            bind_to_forest(f, levl, tsz, true);
        }

        inline void initSparse(const expert_forest *f, int levl, unsigned nnz)
        {
            MEDDLY_DCASSERT(f);
            bind_to_forest(f, levl, nnz, false);
        }

    public:
        //
        // For convenience: get recycled instance and initialize
        //

        static inline unpacked_node* newRedundant(const expert_forest *f,
                int k, node_handle node, bool full)
        {
            unpacked_node* U = New();
            MEDDLY_DCASSERT(U);
            U->initRedundant(f, k, node, full);
            return U;
        }

        static inline unpacked_node* newRedundant(const expert_forest *f,
                int k, long ev, node_handle node, bool full)
        {
            unpacked_node* U = New();
            MEDDLY_DCASSERT(U);
            U->initRedundant(f, k, ev, node, full);
            return U;
        }

        static inline unpacked_node* newRedundant(const expert_forest *f,
                int k, float ev, node_handle node, bool full)
        {
            unpacked_node* U = New();
            MEDDLY_DCASSERT(U);
            U->initRedundant(f, k, ev, node, full);
            return U;
        }

        static inline unpacked_node* newIdentity(const expert_forest *f,
                int k, unsigned i, node_handle node, bool full)
        {
            unpacked_node* U = New();
            MEDDLY_DCASSERT(U);
            U->initIdentity(f, k, i, node, full);
            return U;
        }

        static inline unpacked_node* newIdentity(const expert_forest *f,
                int k, unsigned i, long ev, node_handle node, bool full)
        {
            unpacked_node* U = New();
            MEDDLY_DCASSERT(U);
            U->initIdentity(f, k, i, ev, node, full);
            return U;
        }

        static inline unpacked_node* newIdentity(const expert_forest *f,
                int k, unsigned i, float ev, node_handle node, bool full)
        {
            unpacked_node* U = New();
            MEDDLY_DCASSERT(U);
            U->initIdentity(f, k, i, ev, node, full);
            return U;
        }

        /** Create a zeroed-out full node */
        static inline unpacked_node* newFull(expert_forest *f,
                int levl, unsigned tsz)
        {
            unpacked_node* U = New();
            MEDDLY_DCASSERT(U);
            U->initFull(f, levl, tsz);
            U->clearFullEdges();
            addToBuildList(U, f);
            return U;
        }

        /** Create a zeroed-out sparse node */
        static inline unpacked_node* newSparse(expert_forest *f,
                int levl, unsigned nnzs)
        {
            unpacked_node* U = New();
            MEDDLY_DCASSERT(U);
            U->initSparse(f, levl, nnzs);
            U->clearSparseEdges();
            addToBuildList(U, f);
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
            @param  map     Translation to use on node handles.
                            Allows us to renumber nodes as we write them.
        */
        void write(output &s, const node_handle* map) const;

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
            MEDDLY_DCASSERT(down);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n,
                    is_full ? size : nnzs);
            return down[n];
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
            MEDDLY_DCASSERT(down);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n,
                    is_full ? size : nnzs);
            return down[n];
        }

        /// Set the nth pointer from E, and destroy E.
        inline void set_d(unsigned n, dd_edge &E)
        {
            MEDDLY_DCASSERT(parent == E.parent);
            MEDDLY_DCASSERT(0==edge_bytes);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n,
                    is_full ? size : nnzs);
            down[n] = E.node;
            E.node = 0; // avoid having to adjust the link count
        }

        /** Get the index of the nth non-zero pointer.
            Use only for sparse readers.
        */
        inline unsigned i(unsigned n) const
        {
            MEDDLY_DCASSERT(index);
            MEDDLY_DCASSERT(!is_full);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, nnzs);
            return index[n];
        }

        /** Modify the index of the nth non-zero pointer.
            Use only for sparse readers.
        */
        inline unsigned& i_ref(unsigned n)
        {
            MEDDLY_DCASSERT(index);
            MEDDLY_DCASSERT(!is_full);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n, nnzs);
            return index[n];
        }

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

        /// Set the nth pointer and edge value from E, and destroy E.
        inline void set_de(unsigned n, dd_edge &E)
        {
            MEDDLY_DCASSERT(parent == E.parent);
            MEDDLY_DCASSERT(edge);
            MEDDLY_DCASSERT(edge_bytes);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, n,
                    is_full ? size : nnzs);
            down[n] = E.node;
            memcpy(
                ((char*) edge) + n * edge_bytes, & (E.raw_value), edge_bytes
            );
            E.node = 0; // avoid having to adjust the link count
        }

        /// Get the edge value, as an integer.
        inline void getEdge(unsigned i, long& val) const
        {
            MEDDLY_DCASSERT(sizeof(long) == edge_bytes);
            MEDDLY::EVencoder<long>::readValue(eptr(i), val);
        }

        /// Get the edge value, as a float.
        inline void getEdge(unsigned i, float& val) const
        {
            MEDDLY_DCASSERT(sizeof(float) == edge_bytes);
            MEDDLY::EVencoder<float>::readValue(eptr(i), val);
        }

        /// Set the edge value, as an integer.
        inline void setEdge(unsigned i, long ev)
        {
            MEDDLY_DCASSERT(sizeof(long) == edge_bytes);
            MEDDLY::EVencoder<long>::writeValue(eptr_write(i), ev);
#ifdef DEVELOPMENT_CODE
            long test_ev = 256;
            getEdge(i, test_ev);
            MEDDLY_DCASSERT(test_ev == ev);
#endif
        }

        /// Set the edge value, as a float.
        inline void setEdge(unsigned i, float ev)
        {
            MEDDLY_DCASSERT(sizeof(float) == edge_bytes);
            MEDDLY::EVencoder<float>::writeValue(eptr_write(i), ev);
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

        /// Get the size of this node (full readers only).
        inline unsigned getSize() const
        {
            MEDDLY_DCASSERT(is_full);
            return size;
        }

        /// Get the number of nonzeroes of this node (sparse readers only).
        inline unsigned getNNZs() const
        {
            MEDDLY_DCASSERT(!is_full);
            return nnzs;
        }

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
            return edge_bytes;
        }

        /// Number of bytes per edge
        inline unsigned edgeBytes() const
        {
            return edge_bytes;
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
        void resize(unsigned ns);

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


        /// Called within expert_forest to allocate space.
        ///   @param  p     Parent, for building nodes.
        ///   @param  k     Level number.
        ///   @param  ns    Size of node.
        ///   @param  full  If true, we'll be filling a full reader.
        ///                 Otherwise it is a sparse one.
        void bind_to_forest(const expert_forest* p, int k, unsigned ns, bool full);

        /// Called by node_storage when building an unpacked
        /// node based on how it's stored.
        inline void bind_as_full(bool full)
        {
            is_full = full;
        }

    protected:
        /// Set all down edges (and values) to 0, for full storage
        inline void clearFullEdges()
        {
            MEDDLY_DCASSERT(isFull());
            memset(down, 0, unsigned(size) * sizeof(node_handle));
            if (edge_bytes) {
                memset(edge, 0, unsigned(size) * edge_bytes);
            }
        }

        /// Set all down edges (and values) to 0, for sparse storage
        inline void clearSparseEdges()
        {
            MEDDLY_DCASSERT(isSparse());
            memset(down, 0, unsigned(nnzs) * sizeof(node_handle));
            if (edge_bytes) {
                memset(edge, 0, unsigned(nnzs) * edge_bytes);
            }
        }


    public:
        // Centralized recycling

        // Pull a recycled node off the free list
        // static inline unpacked_node* useUnpackedNode()
        static inline unpacked_node* New()
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
        static inline void addToBuildList(unpacked_node* b, expert_forest *mp)
        {
#ifdef DEBUG_BUILDLIST
            printf("Adding unpacked node at level %d to build list\n", b->getLevel());
#endif
            MEDDLY_DCASSERT(b);
            MEDDLY_DCASSERT(b->eparent == mp);
            b->modparent = mp;

            b->is_in_build_list = true;
            b->next = buildList;
            buildList = b;
        }

        static void removeFromBuildList(unpacked_node* b);
        static void markBuildListChildren(expert_forest* F);
    private:
        static void initStatics();

    private:
        const forest* parent;   // TBD: eventually remove
        const expert_forest* eparent;
        expert_forest* modparent;

        static unpacked_node* freeList;
        static unpacked_node* buildList;

        unpacked_node* next; // for recycled list, and list of nodes being built
        bool is_in_build_list;

        /*
            TBD - extra info that is not hashed
        */
        void* extra_unhashed;
        unsigned ext_uh_alloc;
        unsigned ext_uh_size;
        /*
            Extra info that is hashed
        */
        void* extra_hashed;
        unsigned ext_h_alloc;
        unsigned ext_h_size;
        /*
            Down pointers, indexes, edge values.
        */
        node_handle* down;
        unsigned* index;
        void* edge;
        bool can_be_extensible;
        bool is_extensible;
        unsigned alloc;
        unsigned ealloc;
        unsigned size;
        unsigned nnzs;
        int level;
        unsigned h;
        unsigned char edge_bytes; // number of bytes for an edge value.
        bool is_full;
#ifdef DEVELOPMENT_CODE
        bool has_hash;
#endif
};


#endif

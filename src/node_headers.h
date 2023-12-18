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

#ifndef MEDDLY_NODE_HEADERS_H
#define MEDDLY_NODE_HEADERS_H


#include "arrays.h"
#include "policies.h"

namespace MEDDLY {
    class node_headers;
    class forest;
    class output;

    class memstats;
    class statset;
};

// ******************************************************************
// *                                                                *
// *                       node_headers class                       *
// *                                                                *
// ******************************************************************

/**
    Node header information, for a collection of nodes.

    This is used within a forest to store meta-information
    about each node.  The actual node itself (i.e., the downward
    pointers and any edge values) is stored in a separate data
    structure (the node_storage class) using various encoding schemes.

    Conceptually, this class is simply an array of structs,
    plus a mechanism for recycling unused elements.
    Node handle 0 is reserved for special use, so
    "real" nodes must have a non-zero handle.

    For each node, its header contains:
      node_address  address     The "address" of the node in the
                                node_storage class.

      integer       level       Level number of the node.  0 indicates that
                                the node has been deleted, but we cannot
                                recycle the node handle yet (probably because
                                it might be contained in a compute table
                                somewhere).

      natural       cache_count Optional (i.e., can be turned on/off for all
                                nodes).  Number of compute table references
                                to this handle.

      natural       incoming    Optional (i.e., can be turned on/off for all
                                nodes).  Number of incoming edges to the node
                                referred to by this handle.

    In practice, we may use a different data structure for speed
    or (more likely) to save space.

*/
class MEDDLY::node_headers : public array_watcher {
    public:
        node_headers(forest &P, memstats &ms, statset &ss);
        virtual ~node_headers();

        /// Delay initialization, so we can set up the forest first.
        void initialize();

    public:
        //
        // array_watcher overloads
        //
        virtual void expandElementSize(unsigned oldbits, unsigned newbits);
        virtual void shrinkElementSize(unsigned oldbits, unsigned newbits);


    public:
        //
        // Misc. Inlined getters/setters
        //

        /**
            Get largest handle of active nodes.
        */
        inline node_handle lastUsedHandle() const {
            return a_last;
        }

    public:
        //
        // Node status
        //

        /// Is this a deleted node
        inline bool isDeleted(node_handle p) const {
            MEDDLY_DCASSERT(p>0);
            MEDDLY_DCASSERT(levels);
            return (0==levels->get(size_t(p)));
        }

        /// Is this node active?
        inline bool isActive(node_handle p) const {
            return !isDeleted(p);
        }

        /// Set the given node as inactive/deleted
        inline void deactivate(node_handle p) {
            MEDDLY_DCASSERT(isActive(p));
            MEDDLY_DCASSERT(levels);
            levels->set(size_t(p), 0);
        }


    public:
        //
        // Inlined node attribute setters/getters:
        //      address,
        //      level,
        //      cache count,    (if using reference counts)
        //      in some cache,  (if using mark and sweep)
        //      incoming count, (if using reference counts)
        //      reachable,      (if using mark and sweep)
        //      implicit        (for implicit nodes)
        //

        /// Get the address for node p.
        inline node_address getNodeAddress(node_handle p) const {
            MEDDLY_DCASSERT(p>0);
            MEDDLY_DCASSERT(addresses);
            return addresses->get(size_t(p));
        }

        /// Set the address for node p to a.
        inline void setNodeAddress(node_handle p, node_address a) {
            MEDDLY_DCASSERT(p>0);
            MEDDLY_DCASSERT(addresses);
            addresses->set(size_t(p), a);
        }

        /** Change the address for node p.
            Same as setNodeAddress, except we can verify the old address
            as a sanity check.
        */
        inline void moveNodeAddress(node_handle p, node_address old_addr,
                node_address new_addr)
        {
            MEDDLY_DCASSERT(p>0);
            MEDDLY_DCASSERT(addresses);
            MEDDLY_DCASSERT(old_addr == addresses->get(size_t(p)));
            addresses->set(size_t(p), new_addr);
        }

        // ----------------------------------------------------------

        /// Get the level for node p.
        inline int getNodeLevel(node_handle p) const {
            MEDDLY_DCASSERT(p>0);
            MEDDLY_DCASSERT(levels);
            return levels->get(size_t(p));
        }

        /// Set the level for node p to k.
        inline void setNodeLevel(node_handle p, int k) {
            MEDDLY_DCASSERT(p>0);
            MEDDLY_DCASSERT(levels);
            levels->set(size_t(p), k);
        }

        // ----------------------------------------------------------

        /// Get the cache count for node p.
        inline unsigned long getNodeCacheCount(node_handle p) const {
            MEDDLY_DCASSERT(p>0);
            MEDDLY_DCASSERT(cache_counts);
            return cache_counts->get(size_t(p));
        }


        /// Increment the cache count for node p and return p.
        inline void cacheNode(node_handle p) {
            if (p<1) return;    // terminal node

            MEDDLY_DCASSERT(isActive(p));
            MEDDLY_DCASSERT(cache_counts);
            cache_counts->increment(size_t(p));
#ifdef TRACK_CACHECOUNT
            std::cerr << "\t+Node " << p << " is in " <<
                cache_counts->get(size_t(p)) << " caches" << std::endl;
#endif
        }

        /// Decrement the cache count for node p.
        inline void uncacheNode(node_handle p) {
            if (p<1) return;

            MEDDLY_DCASSERT(cache_counts);

            if (cache_counts->isPositiveAfterDecrement(size_t(p))) {
#ifdef TRACK_CACHECOUNT
                std::cerr << "\t-Node " << p << " is in " <<
                    cache_counts->get(size_t(p)) << " caches" << std::endl;
#endif
                return;
            }

#ifdef TRACK_CACHECOUNT
            std::cerr << "\t-Node " << p << " is not in any caches" << std::endl;
#endif
            lastUncache(p);
        }

        // ----------------------------------------------------------

        /// Indicate that node p is (or might be) in some cache entry.
        inline void setInCacheBit(node_handle p) {
            if (p<1) return;    // terminal node

            MEDDLY_DCASSERT(is_in_cache);
            is_in_cache->set(size_t(p), 1);
        }

        /// Clear cache entry bit for all nodes
        inline void clearAllInCacheBits() {
            MEDDLY_DCASSERT(is_in_cache);
            is_in_cache->clearAll();
        }

        // ----------------------------------------------------------

        /// Get the incoming count for node p.
        inline unsigned long getIncomingCount(node_handle p) const {
            MEDDLY_DCASSERT(p>0);
            MEDDLY_DCASSERT(incoming_counts);
            return incoming_counts->get(size_t(p));
        }

        /// Increment the incoming count for node p and return p.
        inline node_handle linkNode(node_handle p) {
            if (p<1) return p;    // terminal node

            MEDDLY_DCASSERT(isActive(p));
            MEDDLY_DCASSERT(incoming_counts);

            if (incoming_counts->isZeroBeforeIncrement(size_t(p))) {
                reviveNode(p);
            }

#ifdef TRACK_DELETIONS
            std::cerr << "\t+Node " << p << " count now "
                  << incoming_counts->get(size_t(p)) << std::endl;
#endif

            return p;
        }

        /// Decrement the incoming count for node p.
        inline void unlinkNode(node_handle p) {
            if (p<1) return;    // terminal node

            MEDDLY_DCASSERT(isActive(p));
            MEDDLY_DCASSERT(incoming_counts);

            if (incoming_counts->isPositiveAfterDecrement(size_t(p))) {
#ifdef TRACK_DELETIONS
                std::cerr << "\t-Node " << p << " count now "
                      << incoming_counts->get(size_t(p)) << std::endl;
#endif
                return;
            }

#ifdef TRACK_DELETIONS
            std::cerr << "\t-Node " << p << " count now zero" << std::endl;
#endif

            lastUnlink(p);
        }

        // ----------------------------------------------------------

        /// Set the node marker for reachable nodes.
        inline void linkReachable(bitvector* R) {
            is_reachable = R;
        }

        /*
        /// Indicate that node p is (or might be) reachable in the forest.
        inline void setReachableBit(node_handle p) {
            if (p<1) return;    // terminal node
            MEDDLY_DCASSERT(isActive(p));
            MEDDLY_DCASSERT(reachable);
            reachable->setMarked(p);
        }

        /// Mark a node
        inline void setReachable(node_handle p) {
            MEDDLY_DCASSERT(reachable);
            reachable->mark(p);
        }

        /// Is a node marked?
        inline bool isReachable(node_handle p) const {
            return reachable->isMarked(p);
        }

        /// Clear reachable bit for all nodes
        inline void clearAllReachableBits() {
            MEDDLY_DCASSERT(reachable);
            reachable->unmarkAll();
        }

        inline node_marker* getReachableMarker() {
            return reachable;
        }
        */

        // ----------------------------------------------------------

        /// Get whether node p is implicit.
        inline int getNodeImplicitFlag(node_handle p) const {
            MEDDLY_DCASSERT(p>0);
            return (implicit_bits) ? implicit_bits->get(size_t(p)) : false;
        }

        /// Set as true if node p is implicit.
        inline void setNodeImplicitFlag(node_handle p,bool flag) {
            MEDDLY_DCASSERT(p>0);
            if (implicit_bits) {
                implicit_bits->set(size_t(p), flag);
            } else {
                if (!flag) return;
                implicit_bits = new bitvector(this);
                implicit_bits->expand(a_size);
                implicit_bits->set(size_t(p), flag);
            }
        }


    public:
        //
        // Substantial methods
        //

        /// Called when done setting InCacheBits
        void sweepAllInCacheBits();


        //
        // Node handle management
        //

        /**
            Get an unused node handle.
            This is either a recycled one or
            the next one in the available pool
            (which will be expanded if necessary).
        */
        node_handle getFreeNodeHandle();

        /**
            Recycle a used node handle.
            The recycled handle can eventually be
            reused when returned by a call to
            getFreeNodeHandle().
        */
        void recycleNodeHandle(node_handle p);


        /**
            Swap two nodes.
            Used for in-place forest reordering.
                @param  p               First node
                @param  q               Second node
                @param  swap_incounts   If true, swap everything;
                                        if false, don't swap incoming counts.
        */
        void swapNodes(node_handle p, node_handle q, bool swap_incounts);


        //
        // I/O related methods
        //

        /** Show various memory stats.
                @param  s       Output stream to write to
                @param  pad     Padding string, written at the start of
                                each output line.
                @param  flags   Which stats to display, as "flags" (see
                                policies.h); use bitwise or to combine values.
        */
        void reportStats(output &s, const char* pad, display_flags flags) const;

        /** Display header information, primarily for debugging.
                @param  s     Output stream to write to
                @param  p     Node to display
        */
        void showHeader(output &s, node_handle p) const;

        /// Show internal data structures for debugging.
        void dumpInternal(output &s) const;


    private:
        //
        // helper methods
        //

        // for debugging
        void validateFreeLists() const;

        /// Increase the number of node handles.
        void expandHandleList();

        /// Decrease the number of node handles.
        void shrinkHandleList();

        /// Called when the link count reaches zero.
        void lastUnlink(node_handle p);

        /// Called when linking to an unreachable node.
        void reviveNode(node_handle p);

        /// Called when the cache count reaches zero.
        void lastUncache(node_handle p);

        /// Get next deleted node header.
        inline size_t getNextOf(size_t p) const {
            MEDDLY_DCASSERT(isDeleted(p));
            MEDDLY_DCASSERT(addresses);
            return addresses->get(size_t(p));
        }

        /// Set next deleted node header.
        inline void setNextOf(size_t p, size_t n) {
            MEDDLY_DCASSERT(isDeleted(p));
            MEDDLY_DCASSERT(n>=0);
            MEDDLY_DCASSERT(addresses);
            addresses->set(size_t(p), node_address(n));
        }


    private:
        address_array* addresses;
        level_array* levels;
        counter_array* cache_counts;
        bitvector* is_in_cache;
        counter_array* incoming_counts;
        bitvector* implicit_bits;
        bitvector* is_reachable;

        /// Last used address.
        size_t a_last;

        /// Allocated sizes of arrays.
        size_t a_size;

        /// Next time we shink the address list.
        size_t a_next_shrink;

        //
        // Data structure for free lists,
        // used with cache reference counts
        //

        /// Number of addresses in free lists
        size_t a_freed;
        /// Pointer to unused address lists, based on size
        size_t a_unused[8];  // number of bytes per handle
        /// Lowest non-empty address list
        unsigned a_lowest_index;

        //
        // Place to search next for free items,
        // used with mark & sweep cache bits.
        //
        size_t a_sweep;


        /// Current size of node "header"
        unsigned h_bits;

        /// Are we using the pessimistic strategy?
        bool pessimistic;

        /// Parent forest, needed for recycling
        forest &parent;

        /// Memory stats to update
        memstats &mstats;

        /// Other stats to update
        statset &stats;
};



#endif

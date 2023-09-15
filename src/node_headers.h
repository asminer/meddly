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

// #define OLD_NODE_HEADERS

#include "defines.h"
#include <cstring>  // for memset

namespace MEDDLY {
    class node_headers;
    class expert_forest;
    class output;
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
class MEDDLY::node_headers {
    public:
        node_headers(expert_forest &P);
        ~node_headers();


        /** Show various memory stats.
                @param  s       Output stream to write to
                @param  pad     Padding string, written at the start of
                                each output line.
                @param  flags   Which stats to display, as "flags";
                                use bitwise or to combine values.
        */
        void reportStats(output &s, const char* pad, unsigned flags) const;

        /** Display header information, primarily for debugging.
                @param  s     Output stream to write to
                @param  p     Node to display
        */
        void showHeader(output &s, node_handle p) const;

        /**
            Indicate that we don't need to track cache counts.
            The default is that we will track cache counts.
            For efficiency, this should be called soon after
            construction (otherwise we may allocate space for nothing).
        */
        // void turnOffCacheCounts();

        /**
            Indicate that we don't need to track incoming counts.
            The default is that we will track incoming counts.
            For efficiency, this should be called soon after
            construction (otherwise we may allocate space for nothing).
        */
        // void turnOffIncomingCounts();

        /**
            Set node recycling to pessimistic or optimistic.
            Pessimistic:    disconnected nodes are recycled immediately.
            Optimistic:     disconnected nodes are recycled only when
                            they no longer appear in any caches.
        */
        void setPessimistic(bool pess);

    public: // node handle management

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
            Get largest handle of active nodes.
        */
        node_handle lastUsedHandle() const;

        /**
            Swap two nodes.
            Used for in-place forest reordering.
                @param  p               First node
                @param  q               Second node
                @param  swap_incounts   If true, swap everything;
                                        if false, don't swap incoming counts.
        */
        void swapNodes(node_handle p, node_handle q, bool swap_incounts);

    public: // node status
        /// Is this node active?
        bool isActive(node_handle p) const;

        /// Is this a zombie node (dead but not able to be deleted yet)
        // bool isZombie(node_handle p) const;

        /// Is this a deleted node
        bool isDeleted(node_handle p) const;

        /// Deactivated: 0 level
        // bool isDeactivated(node_handle p) const;

        void deactivate(node_handle p);


    public: // address stuff

        /// Get the address for node p.
        node_address getNodeAddress(node_handle p) const;

        /// Set the address for node p to a.
        void setNodeAddress(node_handle p, node_address a);

        /** Change the address for node p.
            Same as setNodeAddress, except we can verify the old address
            as a sanity check.
        */
        void moveNodeAddress(node_handle p, node_address old_addr,
                node_address new_addr);

    public: // level stuff

        /// Get the level for node p.
        int getNodeLevel(node_handle p) const;

        /// Set the level for node p to k.
        void setNodeLevel(node_handle p, int k);

    public: // cache count stuff

        /// Are we tracking cache counts
        // bool trackingCacheCounts() const;

        /// Get the cache count for node p.
        unsigned long getNodeCacheCount(node_handle p) const;

        /// Increment the cache count for node p and return p.
        void cacheNode(node_handle p);

        /// Decrement the cache count for node p.
        void uncacheNode(node_handle p);

    public: // cache bit stuff (if we're not using cache counts)

        /// Indicate that node p is (or might be) in some cache entry.
        void setInCacheBit(node_handle p);

        /// Clear cache entry bit for all nodes
        void clearAllInCacheBits();

        /// Called when done setting InCacheBits
        void sweepAllInCacheBits();

    public: // incoming count stuff

        /// Are we tracking incoming counts
        // bool trackingIncomingCounts() const;

        /// Get the incoming count for node p.
        unsigned long getIncomingCount(node_handle p) const;

        /// Increment the incoming count for node p and return p.
        node_handle linkNode(node_handle p);

        /// Decrement the incoming count for node p.
        void unlinkNode(node_handle p);

    public: // reachable bit stuff (if we're not using incoming counts)

        /// Indicate that node p is (or might be) reachable in the forest.
        void setReachableBit(node_handle p);

        /// Does node p have its reachable bit set?
        bool hasReachableBit(node_handle p) const;

        /// Clear reachable bit for all nodes
        void clearAllReachableBits();

    public: // implicit stuff

        /// Get whether node p is implicit.
        int getNodeImplicitFlag(node_handle p) const;

        /// Set as true if node p is implicit.
        void setNodeImplicitFlag(node_handle p,bool flag);


    public: // for debugging

        void dumpInternal(output &s) const;


    public: // interface for node header size changes
        void changeHeaderSize(unsigned oldbits, unsigned newbits);

    private:  // helper methods

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

        size_t getNextOf(size_t p) const;
        void setNextOf(size_t p, size_t n);

    private:

        /*
            Stores the level of each node,
            as an array. This allows us to
            shrink the array elements based on
            the number of levels, to one of
            char, short, or int.
         */
        class level_array {
                node_headers &parent;
                char* data8;
                short* data16;
                int* data32;
                size_t size;
                unsigned char bytes;
            public:
                level_array(node_headers &p, int max_level);
                ~level_array();

                void expand(size_t ns);
                void shrink(size_t ns);

                void show(output &s, size_t first, size_t last, int width)
                    const;

                inline size_t entry_bits() const {
                    return size_t(bytes) * 8;
                }

                inline int get(size_t i) const
                {
                    MEDDLY_DCASSERT(i<size);
                    if (data8) {
                        MEDDLY_DCASSERT(!data16);
                        MEDDLY_DCASSERT(!data32);
                        return data8[i];
                    }
                    if (data16) {
                        MEDDLY_DCASSERT(!data32);
                        return data16[i];
                    }
                    MEDDLY_DCASSERT(data32);
                    return data32[i];
                }

                inline void set(size_t i, int v)
                {
                    MEDDLY_DCASSERT(i<size);
                    if (data8) {
                        MEDDLY_DCASSERT(!data16);
                        MEDDLY_DCASSERT(!data32);
                        MEDDLY_DCASSERT(v>-128);
                        MEDDLY_DCASSERT(v<128);
                        data8[i] = v;
                        return;
                    }
                    if (data16) {
                        MEDDLY_DCASSERT(!data32);
                        MEDDLY_DCASSERT(v>-32768);
                        MEDDLY_DCASSERT(v<32768);
                        data16[i] = v;
                        return;
                    }
                    MEDDLY_DCASSERT(data32);
                    data32[i] = v;
                }

                inline void swap(size_t i, size_t j)
                {
                    MEDDLY_DCASSERT(i<size);
                    MEDDLY_DCASSERT(j<size);
                    if (data8) {
                        MEDDLY_DCASSERT(!data16);
                        MEDDLY_DCASSERT(!data32);
                        SWAP(data8[i], data8[j]);
                        return;
                    }
                    if (data16) {
                        MEDDLY_DCASSERT(!data32);
                        SWAP(data16[i], data16[j]);
                        return;
                    }
                    MEDDLY_DCASSERT(data32);
                    SWAP(data32[i], data32[j]);
                }

        };

        /**
            Stores incoming or cache counts for each node,
            as an array.  This allows us to disable counts,
            or to shrink the array elements based on the
            maximum count.
        */
        class counter_array {
                node_headers &parent;
                unsigned char* data8;
                unsigned short* data16;
                unsigned int* data32;
                size_t size;
                size_t counts_09bit;  // #counts requiring at least 9 bits
                size_t counts_17bit;  // #counts requiring at least 17 bits
                unsigned char bytes;
            public:
                counter_array(node_headers &p);
                ~counter_array();

                void expand(size_t ns);
                void shrink(size_t ns);

                unsigned int get(size_t i) const;
                void swap(size_t i, size_t j);
                void increment(size_t i);
                void decrement(size_t i);

                bool isZeroBeforeIncrement(size_t i);
                bool isPositiveAfterDecrement(size_t i);

                void show(output &s, size_t first, size_t last, int width)
                    const;
                inline size_t entry_bits() const { return size_t(bytes) * 8; }


                // Expand from 8-bit to 16-bit entries because of element i
                void expand8to16(size_t i);
                // Expand from 16-bit to 32-bit entries because of element i
                void expand16to32(size_t i);

                void shrink16to8(size_t ns);
                void shrink32to16(size_t ns);
                void shrink32to8(size_t ns);
            };

        /**
            Stores the address of each node, as either an
            array of ints or an array of longs.
        */
        class address_array {
                node_headers &parent;
                unsigned int* data32;
                unsigned long* data64;
                size_t size;
                size_t num_large_elements;
                unsigned char bytes;
            public:
                address_array(node_headers &p);
                ~address_array();

                void expand(size_t ns);
                void shrink(size_t ns);

                unsigned long get(size_t i) const;
                void set(size_t i, unsigned long v);
                void swap(size_t i, size_t j);

                void show(output &s, size_t first, size_t last, int width)
                    const;
                inline size_t entry_bits() const { return size_t(bytes) * 8; }


                void expand32to64();
                void shrink64to32(size_t ns);
        };

        /**
            Array of booleans.
            Used for mark and sweep, or any other
            cases where we need a bit per node.
        */
        class bitvector {
                node_headers &parent;
                bool* data;
                size_t size;
            public:
                bitvector(node_headers &p);
                ~bitvector();

                void expand(size_t ns);
                void shrink(size_t ns);

                inline bool get(size_t i) const
                {
                    MEDDLY_DCASSERT(data);
                    MEDDLY_DCASSERT(i<size);
                    return data[i];
                }

                inline void set(size_t i, bool v)
                {
                    MEDDLY_DCASSERT(data);
                    MEDDLY_DCASSERT(i<size);
                    data[i] = v;
                }

                inline void clearAll()
                {
                    if (size) {
                        MEDDLY_DCASSERT(data);
                        memset(data, 0, size * sizeof(bool));
                    }
                }

                inline void swap(size_t i, size_t j)
                {
                    MEDDLY_DCASSERT(data);
                    MEDDLY_DCASSERT(i<size);
                    MEDDLY_DCASSERT(j<size);
                    SWAP(data[i], data[j]);
                }

                inline size_t entry_bits() const { return sizeof(bool) * 8; }

                /// Return smallest index i >= start with bit i cleared.
                inline size_t firstZero(size_t start) const
                {
                    for (; start < size; start++) {
                        if (0==data[start]) return start;
                    }
                    return size;
                }

        };

    private:
        address_array* addresses;
        level_array* levels;
        counter_array* cache_counts;
        bitvector* is_in_cache;
        counter_array* incoming_counts;
        bitvector* is_reachable;

        bitvector* implicit_bits;

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
        char a_lowest_index;

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
        expert_forest &parent;
};

// ******************************************************************
// *                                                                *
// *          inlined node_headers::counter_array  methods          *
// *                                                                *
// ******************************************************************

inline unsigned int MEDDLY::node_headers::counter_array::get(size_t i) const
{
    MEDDLY_DCASSERT(i<size);
    if (data8) {
        MEDDLY_DCASSERT(!data16);
        MEDDLY_DCASSERT(!data32);
        return data8[i];
    }
    if (data16) {
        MEDDLY_DCASSERT(!data32);
        return data16[i];
    }
    MEDDLY_DCASSERT(data32);
    return data32[i];
}

inline void MEDDLY::node_headers::counter_array::swap(size_t i, size_t j)
{
    MEDDLY_DCASSERT(i<size);
    MEDDLY_DCASSERT(j<size);
    if (data8) {
        MEDDLY_DCASSERT(!data16);
        MEDDLY_DCASSERT(!data32);
        SWAP(data8[i], data8[j]);
        return;
    }
    if (data16) {
        MEDDLY_DCASSERT(!data32);
        SWAP(data16[i], data16[j]);
        return;
    }
    MEDDLY_DCASSERT(data32);
    SWAP(data32[i], data32[j]);
}

inline void MEDDLY::node_headers::counter_array::increment(size_t i)
{
    MEDDLY_DCASSERT(i<size);
    if (data8) {
        MEDDLY_DCASSERT(!data16);
        MEDDLY_DCASSERT(!data32);
        if (0 == ++data8[i]) expand8to16(i);
        return;
    }
    if (data16) {
        MEDDLY_DCASSERT(!data32);
        ++data16[i];
        if (256 == data16[i]) ++counts_09bit;
        if (0 == data16[i]) expand16to32(i);
        return;
    }
    MEDDLY_DCASSERT(data32);
    data32[i]++;
    if (256 == data32[i]) ++counts_09bit;
    if (65536 == data32[i]) ++counts_17bit;
}

inline void MEDDLY::node_headers::counter_array::decrement(size_t i)
{
    MEDDLY_DCASSERT(i<size);
    if (data8) {
        MEDDLY_DCASSERT(!data16);
        MEDDLY_DCASSERT(!data32);
        MEDDLY_DCASSERT(data8[i]);
        --data8[i];
        return;
    }
    if (data16) {
        MEDDLY_DCASSERT(!data32);
        MEDDLY_DCASSERT(data16[i]);
        if (256 == data16[i]) {
            MEDDLY_DCASSERT(counts_09bit);
            --counts_09bit;
        }
        --data16[i];
        return;
    }
    MEDDLY_DCASSERT(data32);
    MEDDLY_DCASSERT(data32[i]);
    if (256 == data32[i]) {
        MEDDLY_DCASSERT(counts_09bit);
        --counts_09bit;
    }
    if (65536 == data32[i]) {
        MEDDLY_DCASSERT(counts_17bit);
        --counts_17bit;
    }
    --data32[i];
}

inline bool MEDDLY::node_headers::counter_array::isZeroBeforeIncrement(size_t i)
{
    MEDDLY_DCASSERT(i<size);
    if (data8) {
        MEDDLY_DCASSERT(!data16);
        MEDDLY_DCASSERT(!data32);
        if (0==data8[i]) {
            data8[i] = 1;
            return true;
        }
        if (0 == ++data8[i]) expand8to16(i);
        return false;
    }
    if (data16) {
        MEDDLY_DCASSERT(!data32);
        if (0==data16[i]) {
            data16[i] = 1;
            return true;
        }
        ++data16[i];
        if (256 == data16[i]) ++counts_09bit;
        if (0 == data16[i]) expand16to32(i);
        return false;
    }
    MEDDLY_DCASSERT(data32);
    if (0==data32[i]) {
        data32[i] = 1;
        return true;
    }
    data32[i]++;
    if (256 == data32[i]) ++counts_09bit;
    if (65536 == data32[i]) ++counts_17bit;
    return false;
}

inline bool MEDDLY::node_headers::counter_array
    ::isPositiveAfterDecrement(size_t i)
{
    MEDDLY_DCASSERT(i<size);
    if (data8) {
        MEDDLY_DCASSERT(!data16);
        MEDDLY_DCASSERT(!data32);
        MEDDLY_DCASSERT(data8[i]);
        return --data8[i] > 0;
    }
    if (data16) {
        MEDDLY_DCASSERT(!data32);
        MEDDLY_DCASSERT(data16[i]);
        if (256 == data16[i]) {
            MEDDLY_DCASSERT(counts_09bit);
            --counts_09bit;
        }
        return --data16[i] > 0;
    }
    MEDDLY_DCASSERT(data32);
    MEDDLY_DCASSERT(data32[i]);
    if (256 == data32[i]) {
        MEDDLY_DCASSERT(counts_09bit);
        --counts_09bit;
    }
    if (65536 == data32[i]) {
        MEDDLY_DCASSERT(counts_17bit);
        --counts_17bit;
    }
    return --data32[i] > 0;
}


// ******************************************************************
// *                                                                *
// *          inlined node_headers::address_array  methods          *
// *                                                                *
// ******************************************************************

inline unsigned long MEDDLY::node_headers::address_array::get(size_t i) const
{
    MEDDLY_DCASSERT(i<size);
    if (4==bytes) {
        MEDDLY_DCASSERT(data32);
        MEDDLY_DCASSERT(!data64);
        MEDDLY_DCASSERT(0==num_large_elements);
        return data32[i];
    }
    MEDDLY_DCASSERT(8==bytes);
    MEDDLY_DCASSERT(!data32);
    MEDDLY_DCASSERT(data64);
    return data64[i];
}

inline void MEDDLY::node_headers::address_array::set(size_t i, unsigned long v)
{
    MEDDLY_DCASSERT(i<size);
    if (4==bytes) {
        MEDDLY_DCASSERT(data32);
        MEDDLY_DCASSERT(!data64);
        MEDDLY_DCASSERT(0==num_large_elements);

        if (v & 0xffffffff00000000) {
            // v won't fit in 32 bits
            expand32to64();
            MEDDLY_DCASSERT(data64);
            data64[i] = v;
        } else {
            // v will fit in 32 bits
            data32[i] = v;
        }
        return;
    }
    MEDDLY_DCASSERT(8==bytes);
    MEDDLY_DCASSERT(!data32);
    MEDDLY_DCASSERT(data64);
    if (v & 0xffffffff00000000) {
        // v is large
        if (0 == (data64[i] & 0xffffffff00000000)) {
            // replacing small
            num_large_elements++;
        }
    } else {
        // v is small
        if (data64[i] & 0xffffffff00000000) {
            // replacing large
            MEDDLY_DCASSERT(num_large_elements);
            num_large_elements--;
        }
    }
    data64[i] = v;
}

inline void MEDDLY::node_headers::address_array::swap(size_t i, size_t j)
{
    MEDDLY_DCASSERT(i<size);
    MEDDLY_DCASSERT(j<size);
    if (4==bytes) {
        MEDDLY_DCASSERT(data32);
        MEDDLY_DCASSERT(!data64);
        SWAP(data32[i], data32[j]);
        return;
    }
    MEDDLY_DCASSERT(8==bytes);
    MEDDLY_DCASSERT(!data32);
    MEDDLY_DCASSERT(data64);
    SWAP(data64[i], data64[j]);
}

// ******************************************************************
// *                                                                *
// *                  inlined node_headers methods                  *
// *                                                                *
// ******************************************************************

inline void MEDDLY::node_headers::setPessimistic(bool pess)
{
  pessimistic = pess;
}

// ******************************************************************

inline MEDDLY::node_handle
MEDDLY::node_headers::lastUsedHandle() const
{
  return a_last;
}

// ******************************************************************

inline bool
MEDDLY::node_headers::isActive(node_handle p) const
{
  return !isDeleted(p);
}

// ******************************************************************

/*
inline bool
MEDDLY::node_headers::isZombie(node_handle p) const
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  MEDDLY_DCASSERT(addresses);
  MEDDLY_DCASSERT(levels);
  return (0==addresses->get(size_t(p))) && (0!=levels->get(size_t(p)));
}
*/

// ******************************************************************

inline bool
MEDDLY::node_headers::isDeleted(node_handle p) const
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(levels);
    return (0==levels->get(size_t(p)));
}

// ******************************************************************

/*
inline bool
MEDDLY::node_headers::isDeactivated(node_handle p) const
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  MEDDLY_DCASSERT(levels);
  return (0==levels->get(size_t(p)));
}
*/

// ******************************************************************

inline void
MEDDLY::node_headers::deactivate(node_handle p)
{
    MEDDLY_DCASSERT(isActive(p));
    MEDDLY_DCASSERT(levels);
    levels->set(size_t(p), 0);
}


// ******************************************************************

inline MEDDLY::node_address
MEDDLY::node_headers::getNodeAddress(node_handle p) const
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(addresses);
    return addresses->get(size_t(p));
}

// ******************************************************************

inline void MEDDLY::node_headers::setNodeAddress(node_handle p, node_address a)
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(addresses);
    addresses->set(size_t(p), a);
}

// ******************************************************************

inline void MEDDLY::node_headers::moveNodeAddress(node_handle p,
  node_address old_addr, node_address new_addr)
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(addresses);
    MEDDLY_DCASSERT(old_addr == addresses->get(size_t(p)));
    addresses->set(size_t(p), new_addr);
}

// ******************************************************************

inline int
MEDDLY::node_headers::getNodeLevel(node_handle p) const
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(levels);
    return levels->get(size_t(p));
}

// ******************************************************************

inline void
MEDDLY::node_headers::setNodeLevel(node_handle p, int k)
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(levels);
    levels->set(size_t(p), k);
}

// ******************************************************************

/*
inline bool
MEDDLY::node_headers::trackingCacheCounts() const
{
  return cache_counts;
}
*/

// ******************************************************************

inline unsigned long
MEDDLY::node_headers::getNodeCacheCount(node_handle p) const
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(cache_counts);
    return cache_counts->get(size_t(p));
}

// ******************************************************************

inline void
MEDDLY::node_headers::cacheNode(node_handle p)
{
    if (p<1) return;    // terminal node

    MEDDLY_DCASSERT(isActive(p));
    MEDDLY_DCASSERT(cache_counts);
    cache_counts->increment(size_t(p));
#ifdef TRACK_CACHECOUNT
    fprintf(stdout, "\t+Node %d is in %u caches\n",
            p, cache_counts->get(size_t(p)));
    fflush(stdout);
#endif
}

// ******************************************************************

inline void
MEDDLY::node_headers::uncacheNode(MEDDLY::node_handle p)
{
    if (p<1) return;

    MEDDLY_DCASSERT(addresses);
    MEDDLY_DCASSERT(cache_counts);

    if (cache_counts->isPositiveAfterDecrement(size_t(p))) {
#ifdef TRACK_CACHECOUNT
        fprintf(stdout, "\t-Node %d is in %lu caches\n",
                p, cache_counts->get(size_t(p)));
        fflush(stdout);
#endif
        return;
    }

#ifdef TRACK_CACHECOUNT
    fprintf(stdout, "\t-Node %d is not in any caches\n", p);
    fflush(stdout);
#endif

    lastUncache(p);
}

// ******************************************************************

inline void
MEDDLY::node_headers::setInCacheBit(node_handle p)
{
    if (p<1) return;    // terminal node

    MEDDLY_DCASSERT(is_in_cache);
    is_in_cache->set(size_t(p), 1);
}

// ******************************************************************

inline void
MEDDLY::node_headers::clearAllInCacheBits()
{
    MEDDLY_DCASSERT(is_in_cache);
    is_in_cache->clearAll();
}

// ******************************************************************

/*
inline bool
MEDDLY::node_headers::trackingIncomingCounts() const
{
  return incoming_counts;
}
*/

// ******************************************************************

inline unsigned long
MEDDLY::node_headers::getIncomingCount(node_handle p) const
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(incoming_counts);
    return incoming_counts->get(size_t(p));
}

// ******************************************************************

inline MEDDLY::node_handle
MEDDLY::node_headers::linkNode(node_handle p)
{
    if (p<1) return p;    // terminal node

    MEDDLY_DCASSERT(isActive(p));
    MEDDLY_DCASSERT(incoming_counts);
    if (incoming_counts->isZeroBeforeIncrement(size_t(p))) {
        reviveNode(p);
    }

#ifdef TRACK_DELETIONS
    std::cout << "\t+Node " << p << " count now " << incoming_counts->get(size_t(p)) << std::endl;
#endif

    return p;
}

// ******************************************************************

inline void
MEDDLY::node_headers::unlinkNode(node_handle p)
{
    if (p<1) return;    // terminal node

    MEDDLY_DCASSERT(isActive(p));

    MEDDLY_DCASSERT(incoming_counts);
    MEDDLY_DCASSERT(addresses);

    if (incoming_counts->isPositiveAfterDecrement(size_t(p))) {
#ifdef TRACK_DELETIONS
        std::cout << "\t-Node " << p << " count now " << incoming_counts->get(size_t(p)) << std::endl;
#endif
        return;
    }

#ifdef TRACK_DELETIONS
    std::cout << "\t-Node " << p << " count now zero" << std::endl;
#endif

    lastUnlink(p);
}

// ******************************************************************

inline void
MEDDLY::node_headers::setReachableBit(node_handle p)
{
    if (p<1) return;    // terminal node
    MEDDLY_DCASSERT(isActive(p));
    MEDDLY_DCASSERT(is_reachable);
    is_reachable->set(size_t(p), 1);
}

// ******************************************************************

inline bool
MEDDLY::node_headers::hasReachableBit(node_handle p) const
{
    if (p<1) return 1;    // terminal node

    MEDDLY_DCASSERT(is_reachable);
    return is_reachable->get(size_t(p));
}

// ******************************************************************

inline void
MEDDLY::node_headers::clearAllReachableBits()
{
    MEDDLY_DCASSERT(is_reachable);
    is_reachable->clearAll();
}

// ******************************************************************

inline int
MEDDLY::node_headers::getNodeImplicitFlag(node_handle p) const
{
    MEDDLY_DCASSERT(p>0);
    return (implicit_bits) ? implicit_bits->get(size_t(p)) : false;
}

// ******************************************************************

inline void
MEDDLY::node_headers::setNodeImplicitFlag(node_handle p, bool flag)
{
    MEDDLY_DCASSERT(p>0);
    if (implicit_bits) {
        implicit_bits->set(size_t(p), flag);
    } else {
        if (!flag) return;
        implicit_bits = new bitvector(*this);
        implicit_bits->expand(a_size);
        implicit_bits->set(size_t(p), flag);
    }
}

// ******************************************************************

inline size_t
MEDDLY::node_headers::getNextOf(size_t p) const
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(addresses);
    MEDDLY_DCASSERT(levels);
    MEDDLY_DCASSERT(0==levels->get(size_t(p)));
    return addresses->get(size_t(p));
}

// ******************************************************************

inline void
MEDDLY::node_headers::setNextOf(size_t p, size_t n)
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(n>=0);
    MEDDLY_DCASSERT(addresses);
    MEDDLY_DCASSERT(levels);
    MEDDLY_DCASSERT(0==levels->get(size_t(p)));
    addresses->set(size_t(p), node_address(n));
}

#endif

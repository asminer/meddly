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


#include "defines.h"
#include "node_headers.h"
#include "io.h"
#include "operators.h"
#include "forest.h"
#include "storage/bytepack.h"
#include <climits>

// #define DEBUG_HANDLE_MGT

// #define DEBUG_HANDLE_FREELIST

// #define ALWAYS_CHECK_FREELISTS

// #define DEBUG_SWEEP
// #define DEBUG_SWEEP_DETAIL
// #define DEBUG_SWEEP_DETAIL_CT

// ******************************************************************
// *                                                                *
// *                        Helper functions                        *
// *                                                                *
// ******************************************************************


// ******************************************************************

#ifdef DEBUG_HANDLE_FREELIST
void print_sequence(long a)
{
    static bool printed;
    static long first = 0;
    static long last = -1;
    if (a<=0) {
        if (first>last) return;
        if (printed) printf(", ");
        printed = false;
        if (first < last) {
            if (first+1<last) {
                printf("%ld ... %ld", first, last);
            } else {
                printf("%ld, %ld", first, last);
            }
        } else {
            printf("%ld", first);
        }
        first = 0;
        last = -1;
        return;
    }
    // a > 0
    if (0==first) {
        first = last = a;
        return;
    }
    if (last+1 == a) {
        last++;
        return;
    }
    // break in the sequence, we need to print
    if (printed) printf(", "); else printed = true;
    if (first < last) {
        if (first+1<last) {
            printf("%ld ... %ld", first, last);
        } else {
            printf("%ld, %ld", first, last);
        }
    } else {
        printf("%ld", first);
    }
    first = last = a;
}

// ******************************************************************

inline void dump_handle_info(const MEDDLY::node_headers &NH, long size)
{
    printf("Used handles:  ");
    print_sequence(0);
    for (long i=1; i<size; i++) if (!NH.isDeleted(i)) {
        print_sequence(i);
    }
    print_sequence(0);
    printf("\nFree handles:  ");
    for (long i=1; i<size; i++) if (NH.isDeleted(i)) {
        print_sequence(i);
    }
    print_sequence(0);
    printf("\n");
}
#endif

// ******************************************************************
//
// Helpers for dealing with node header array sizes
//
// ******************************************************************

const size_t START_SIZE = 512;
// const size_t MAX_ADD = 65536;
const size_t MAX_ADD = 16777216;

inline size_t next_size(size_t s)
{
    if (0==s) return START_SIZE;
    if (s < MAX_ADD) return s*2;
    return s+MAX_ADD;
}

inline size_t prev_size(size_t s)
{
    if (s <= MAX_ADD) return s/2;
    return s-MAX_ADD;
}

inline size_t next_check(size_t s)
{
    if (0==s) return START_SIZE/4;
    if (s < MAX_ADD) return s*2;
    return s+MAX_ADD;
}

inline size_t prev_check(size_t s)
{
    if (s <= MAX_ADD) return s/2;
    return s-MAX_ADD;
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_headers methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_headers::node_headers(forest &P, memstats &ms, statset &ss)
  : parent(P), mstats(ms), stats(ss)
{
    addresses = nullptr;
    levels = nullptr;
#ifdef REFCOUNTS_ON
    cache_counts = nullptr;
#endif
    is_in_cache = nullptr;
#ifdef REFCOUNTS_ON
    incoming_counts = nullptr;
#endif
    implicit_bits = nullptr;
    is_reachable = nullptr;
    pessimistic = true;

    //
    // Set up handle information and recycling
    //
    a_last = 0;
    a_size = 0;
    a_next_shrink = 0;
    a_freed = 0;
    for (unsigned i=0; i<8; i++) a_unused[i] = 0;
    a_lowest_index = 8;
    a_sweep = SIZE_MAX;
    h_bits = 0;

}

// ******************************************************************

MEDDLY::node_headers::~node_headers()
{
    delete addresses;
    delete levels;
#ifdef REFCOUNTS_ON
    delete cache_counts;
#endif
    delete is_in_cache;
#ifdef REFCOUNTS_ON
    delete incoming_counts;
#endif
    delete implicit_bits;
    // DON'T delete is_reachable, we don't own it
    // delete is_reachable;
}

// ******************************************************************

void MEDDLY::node_headers::initialize()
{
    addresses = new address_array(this);
    levels = new level_array(parent.getNumVariables(), this);
    if (parent.getPolicies().useReferenceCounts) {
#ifdef REFCOUNTS_ON
        cache_counts = new counter_array(this);
        incoming_counts = new counter_array(this);
#else
        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
#endif
    } else {
        is_in_cache = new bitvector(this);
    }

    pessimistic = parent.getPolicies().isPessimistic();
}

// ******************************************************************

void MEDDLY::node_headers::expandElementSize(unsigned oldbits, unsigned newbits)
{
    MEDDLY_DCASSERT(oldbits <= newbits);

    mstats.incMemUsed((a_last - a_freed)*(newbits-oldbits)/8);
    mstats.incMemAlloc(a_size*(newbits-oldbits)/8);

    h_bits += newbits;
    h_bits -= oldbits;
}

// ******************************************************************

void MEDDLY::node_headers::shrinkElementSize(unsigned oldbits, unsigned newbits)
{
    MEDDLY_DCASSERT(oldbits >= newbits);

    mstats.decMemUsed((a_last - a_freed)*(oldbits-newbits)/8);
    mstats.decMemAlloc(a_size*(oldbits-newbits)/8);

    h_bits += newbits;
    h_bits -= oldbits;
}

// ******************************************************************

void MEDDLY::node_headers::sweepAllInCacheBits()
{
    a_sweep = 0;

#ifdef DEBUG_SWEEP
    std::cerr << "Forest " << parent.FID() << " finished cache scan\n";
#endif
#ifdef DEBUG_SWEEP_DETAIL
    std::cerr << "Done cache scan; nodes that will be reclaimed:\n\t";

    MEDDLY_DCASSERT(is_in_cache);
    size_t p = 0;
    for (;;) {
        p = is_in_cache->firstZero(++p);
        if (is_reachable->get(p)) continue;
        std::cerr << p << ". ";
    } // infinite loop
    std::cerr << ".\n";

#ifdef DEBUG_SWEEP_DETAIL_CT
    std::cerr << "Compute tables:\n";
    stream_output s(stderr);
    operation::showAllComputeTables(s, 5);
#endif
#endif

    // TBD - scan backwards and collapse a_last, and then shrink array if needed
}


// ******************************************************************

MEDDLY::node_handle MEDDLY::node_headers::getFreeNodeHandle()
{
    mstats.incMemUsed(h_bits/8);
    node_handle found = 0;
    for (unsigned i=a_lowest_index; i<8; i++) {
        // try the various lists
        while (a_unused[i] > a_last) {
            a_unused[i] = getNextOf(a_unused[i]);
        }
        if (a_unused[i]) {  // get a recycled one from list i
            found = a_unused[i];
            a_unused[i] = getNextOf(a_unused[i]);
            break;
        } else {
            if (i == a_lowest_index) a_lowest_index++;
        }
    }

    if (found) {
        //
        // We got something from a free list.
        //
#ifdef DEBUG_HANDLE_FREELIST
        // address[found].setNotDeleted();
        dump_handle_info(*this, a_last+1);
#endif
#ifdef DEBUG_HANDLE_MGT
        std::cerr << "Forest " << parent.FID()
            << " using recycled handle " << found << std::endl;
#endif
        a_freed--;
        return found;
    }

    //
    // Free lists are empty.
    // Use sweep index if it's small enough (this is set
    // if we're using mark and sweep).
    //
    if (a_sweep < a_last) {
        MEDDLY_DCASSERT(addresses);
        MEDDLY_DCASSERT(is_in_cache);
        MEDDLY_DCASSERT(is_reachable);
        for (;;) {
            a_sweep = is_in_cache->firstZero(++a_sweep);
            if (a_sweep > a_last) break;
            if (is_reachable->get(a_sweep)) continue;
            found = a_sweep;
            break;
        }
        if (found) {
            MEDDLY_DCASSERT(0==addresses->get((size_t)found));
        } else {
            // we're done sweeping
            a_sweep = SIZE_MAX;
        }
    }

    if (found) {
        //
        // We're using a swept node handle.
        //
#ifdef DEBUG_HANDLE_FREELIST
        // address[found].setNotDeleted();
        dump_handle_info(*this, a_last+1);
#endif
#ifdef DEBUG_HANDLE_MGT
        std::cerr << "Forest " << parent.FID()
            << " using swept handle " << found << std::endl;
#endif
        return found;
    }

    //
    // Nothing free; pull off the end of the list,
    // expanding the arrays if necessary.
    //

    if (a_last+1 >= a_size) {
        expandHandleList();
    }
    a_last++;
    MEDDLY_DCASSERT(a_last < a_size);

#ifdef DEBUG_HANDLE_FREELIST
    dump_handle_info(*this, a_last+1);
#endif
#ifdef DEBUG_HANDLE_MGT
    std::cerr << "Forest " << parent.FID()
        << " using end handle " << a_last << std::endl;
#endif
    return a_last;
}

// ******************************************************************

void MEDDLY::node_headers::recycleNodeHandle(node_handle p)
{
    MEDDLY_DCASSERT(p>0);
    const size_t pp = size_t(p);

    MEDDLY_DCASSERT(0==getNodeCacheCount(p));

#ifdef DEBUG_HANDLE_MGT
    std::cerr << "Forest " << parent.FID()
        << " recycling handle " << p << std::endl;
#endif

    mstats.decMemUsed(h_bits/8);
    a_freed++;

    // Determine which list to add this into,
    // and add it there
    unsigned i = bytesRequiredForDown(p) -1;
    setNextOf(pp, a_unused[i]);
    a_unused[i] = pp;
    a_lowest_index = MIN(a_lowest_index, i);

    // if this was the last node, collapse nodes into the
    // "not yet allocated" pile.  But, we don't remove them
    // from the free list(s); we simply discard any too-large
    // ones when we pull from the free list(s).
    if (pp == a_last) {
        while (a_last && isDeleted(a_last) && (0==getNodeCacheCount(a_last))) {
            a_last--;
            a_freed--;
        }

#ifdef DEBUG_HANDLE_MGT
        std::cerr << "Forest " << parent.FID()
            << " collapsing handle end from " << p
            << " to " << a_last << std::endl;
#endif

        if (a_last < a_next_shrink) {
#ifdef DEBUG_HANDLE_MGT
            std::cerr << "Forest " << parent.FID()
                << " triggered shrink check " << a_next_shrink << std::endl;
#endif
            shrinkHandleList();
        }

    }

#ifdef DEBUG_HANDLE_FREELIST
    dump_handle_info(*this, a_last+1);
#endif
#ifdef ALWAYS_CHECK_FREELISTS
    validateFreeLists();
#endif
}

// ******************************************************************

void MEDDLY::node_headers::swapNodes(node_handle p, node_handle q,
        bool swap_incounts)
{
    MEDDLY_DCASSERT(p!=q);
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(q>0);

    const size_t sp = size_t(p);
    const size_t sq = size_t(q);
    if (addresses)                          addresses->swap(sp, sq);
    if (levels)                             levels->swap(sp, sq);
#ifdef REFCOUNTS_ON
    if (cache_counts)                       cache_counts->swap(sp, sq);
    if (swap_incounts && incoming_counts)   incoming_counts->swap(sp, sq);
#endif
    if (implicit_bits)                      implicit_bits->swap(sp, sq);
}

// ******************************************************************

template <class A>
inline void show_num_bits(MEDDLY::output &s, const char* pad, const char* what, const A* array)
{
    if (array) {
        s << pad << "    " << what;
        s.put(array->entry_bits(), 3);
        s << " bits\n";
    }
}

void MEDDLY::node_headers
::reportStats(output &s, const char* pad, display_flags flags) const
{
    if (flags & (STORAGE_STATS | STORAGE_DETAILED) ) {
        s << pad << "Node headers: array based\n";
    }
    if (flags & STORAGE_STATS) {
        s << pad << "  Node header : ";
        s.put(h_bits, 3);
        s << " bits\n";
    }

    if (flags & STORAGE_DETAILED) {
        show_num_bits(s, pad, "address   : ", addresses);
        show_num_bits(s, pad, "level     : ", levels);
#ifdef REFCOUNTS_ON
        show_num_bits(s, pad, "#caches   : ", cache_counts);
#endif
        show_num_bits(s, pad, "in_cache? : ", is_in_cache);
#ifdef REFCOUNTS_ON
        show_num_bits(s, pad, "#incoming : ", incoming_counts);
#endif
        show_num_bits(s, pad, "reachable?: ", is_reachable);
        show_num_bits(s, pad, "implicit? : ", implicit_bits);
    }

    if (flags & STORAGE_STATS) {
        s << pad << "  " << a_size << " headers allocated\n";
        s << pad << "  " << a_last << " last header\n";
        s << pad << "  " << a_freed << " recycled headers\n";
        s << pad << "  " << a_last - a_freed << " headers in use\n";
    }

}

// ******************************************************************

template <class A>
inline void show_element(MEDDLY::output &s, const char* what, const A* array,
            MEDDLY::node_handle p)
{
    if (array) {
        s << what << (unsigned long)array->get((size_t)p);
    }
}

void MEDDLY::node_headers::showHeader(output &s, node_handle p) const
{
    MEDDLY_DCASSERT(p>0);

    int k = ABS(getNodeLevel(p));
    const variable* v = parent.getDomain()->getVar(
            unsigned(parent.getVarByLevel(k))
    );

    if (v->hasName()) {
        s << " level: " << v->getName();
    } else {
        s << " level: " <<  k;
    }
    s.put( (getNodeLevel(p) < 0) ? '\'' : ' ' );

#ifdef REFCOUNTS_ON
    show_element(s, " in: ", incoming_counts, p);
    show_element(s, " cc: ", cache_counts, p);
#endif
    show_element(s, " reach: ", is_reachable, p);
    show_element(s, " cache: ", is_in_cache, p);
}

// ******************************************************************

template <class A>
inline void show_array(MEDDLY::output &s, const char* what, const A* array,
        size_t a_last, int awidth)
{
    s << "    " << what;
    if (array) {
        array->show(s, 1, a_last, awidth);
        s.put('\n');
    } else {
        s << "(null)\n";
    }
}

void MEDDLY::node_headers::dumpInternal(output &s) const
{
    s << "Node headers and management:\n";
    for (unsigned i=0; i<8; i++) {
        s << "    First " << i << "-byte unused node index: "
          << a_unused[i] << "\n";
    }
    int awidth = (int) digits(a_last);
    s << "    Node# :  ";
    for (unsigned long p=1; p<=a_last; p++) {
        if (p>1) s.put(' ');
        s.put(p, awidth);
    }
    s.put('\n');

    show_array(s, "Level  : ", levels, a_last, awidth);
    show_array(s, "Offset : ", addresses, a_last, awidth);
#ifdef REFCOUNTS_ON
    show_array(s, "Incount: ", incoming_counts, a_last, awidth);
    show_array(s, "Cache  : ", cache_counts, a_last, awidth);
#endif
}

// ******************************************************************

void MEDDLY::node_headers::validateFreeLists() const
{
    const unsigned MAX_ERRORS = 25;
    unsigned errors = 0;
    std::cerr << "  Validating free lists\n";
    bool* inlist = new bool[a_last+1];
    for (size_t i=0; i<=a_last; i++) {
        inlist[i] = false;
    }
    for (unsigned i=0; i<8; i++) {
        for (size_t curr = a_unused[i]; curr; curr = getNextOf(curr)) {
            if (curr > a_last) continue;
            inlist[curr] = true;
            if (isDeleted(curr)) continue;
            std::cerr << "\tBAD FREE LIST: item " << curr
                << " in free list#" << i << " is not deleted?\n",
            ++errors;
            if (errors >= MAX_ERRORS) {
                std::cerr << "\tToo many errors, not printing the rest\n";
                exit(1);
            }
        }
    }
    for (size_t i=1; i<=a_last; i++) {
        if (!isDeleted(i)) continue;
        if (inlist[i]) continue;
        std::cerr << "\tMISSING: item " << i
            << " is deleted, not in any free list\n";
        ++errors;
        if (errors >= MAX_ERRORS) {
            std::cerr << "\tToo many errors, not printing the rest\n";
            exit(1);
        }
    }
    std::cerr << "  Done validating free lists\n";
    delete[] inlist;
    if (errors>0) exit(1);
}

// ******************************************************************

void MEDDLY::node_headers::expandHandleList()
{
    //
    // If we're not using reference counts,
    // determine which nodes are reachable.
    //
    if (!parent.getPolicies().useReferenceCounts) {

#ifdef DEBUG_SWEEP
        std::cerr << "Forest " << parent.FID() << " starting unreachable scan\n";
#endif

        parent.markAllRoots();

#ifdef DEBUG_SWEEP
        std::cerr << "Forest " << parent.FID() << " finished unreachable scan\n";
#endif
#ifdef DEBUG_SWEEP_DETAIL
        std::cerr << "Unreachable nodes:\n\t";

        for (node_handle p=1; p<=a_last; p++) {
            MEDDLY_DCASSERT(is_reachable);
            if (is_reachable->get(p)) continue;
            std::cerr << p << ", ";
        }
        std::cerr << ".\n";
#endif  // DEBUG_SWEEP_DETAIL

        //
        // Because there's no way for us to recover
        // an unreachable node, go ahead and delete them now.
        //
        MEDDLY_DCASSERT(is_reachable);
        size_t unrch = 0;
        for (;;) {
            unrch = is_reachable->firstZero(++unrch);
            if (unrch >= a_size) break;
            if (!isDeleted(unrch)) parent.deleteNode(unrch);
        }

    } // no incoming counts

    //
    // Done with GC, now work on expanding
    //

    size_t old_size = a_size;
    do {
        a_next_shrink = next_check(a_next_shrink);
        a_size = next_size(a_size);
    } while (a_last+1 >= a_size);

    mstats.incMemAlloc( ((a_size-old_size)*h_bits)/8 );

    if (addresses)        addresses->expand(a_size);
    if (levels)           levels->expand(a_size);
#ifdef REFCOUNTS_ON
    if (cache_counts)     cache_counts->expand(a_size);
#endif
    if (is_in_cache)      is_in_cache->expand(a_size);
#ifdef REFCOUNTS_ON
    if (incoming_counts)  incoming_counts->expand(a_size);
#endif
    if (is_reachable)     is_reachable->expand(a_size);
    if (implicit_bits)    implicit_bits->expand(a_size);

#ifdef DEBUG_ADDRESS_RESIZE
    std::cerr << "Expanded node header arrays, new size "
        << a_size << " (shrink check " << a_next_shrink << ")\n";
#endif
}


// ******************************************************************

void MEDDLY::node_headers::shrinkHandleList()
{
#ifdef DEBUG_ADDRESS_RESIZE
    std::cerr << "About to shrink node headers arrays (current size "
        << a_size << ")...\n";
    validateFreeLists();
#endif

    // clean out free lists, because we're about
    // to shrink the address list which holds the pointers
    for (unsigned i=0; i<8; i++) {
        //
        // clean list i
        //
        size_t prev = 0;
        size_t curr;
        for (curr = a_unused[i]; curr; curr=getNextOf(curr))
        {
            if (curr > a_last) continue;  // don't add to the list
            if (prev) {
                setNextOf(prev, curr);
            } else {
                a_unused[i] = curr;
            }
            prev = curr;
        }
        if (prev) {
            setNextOf(prev, 0);
        } else {
            a_unused[i] = 0;
        }
    } // for i

    size_t old_size = a_size;
    do {
        a_size = prev_size(a_size);
        a_next_shrink = (a_size > START_SIZE) ? prev_check(a_next_shrink) : 0;
    } while (a_last < a_next_shrink);

    mstats.decMemAlloc( ((old_size-a_size)*h_bits)/8 );

    if (addresses)        addresses->shrink(a_size);
    if (levels)           levels->shrink(a_size);
#ifdef REFCOUNTS_ON
    if (cache_counts)     cache_counts->shrink(a_size);
    if (incoming_counts)  incoming_counts->shrink(a_size);
#endif
    if (is_in_cache)      is_in_cache->shrink(a_size);
    if (is_reachable)     is_reachable->shrink(a_size);
    if (implicit_bits)    implicit_bits->shrink(a_size);

#ifdef DEBUG_ADDRESS_RESIZE
    std::cerr << "Shrank node header arrays, new size "
        << a_size << " (shrink check " << a_next_shrink << ")\n";
    std::cerr << "Done shrinking node headers arrays (size is now "
        << a_size << ")...\n";
    validateFreeLists();
#endif
}


// ******************************************************************

void MEDDLY::node_headers::lastUnlink(node_handle p)
{
    //
    // Node p is unreachable.  See if we need to do some cleanup
    //

    if (
#ifdef REFCOUNTS_ON
        (cache_counts && (0==cache_counts->get(size_t(p))))
        ||
#endif
        (is_in_cache && (0==is_in_cache->get(size_t(p))))
    ) {
        //
        // Unreachable and not in any caches.  Delete and recycle handle.
        //
        parent.deleteNode(p);
        recycleNodeHandle(p);
        return;
    }

    if (pessimistic) {
        //
        // Delete; keep handle until caches are cleared.
        //
        parent.deleteNode(p);
    } else {
        //
        // Optimistic.  Keep unreachables around.
        //
#ifdef TRACK_UNREACHABLE_NODES
        stats.unreachable_nodes++;
#endif
    }
}

// ******************************************************************

void MEDDLY::node_headers::reviveNode(node_handle p)
{
    // Reclaim an unreachable
    stats.reclaimed_nodes++;
#ifdef TRACK_UNREACHABLE_NODES
    stats.unreachable_nodes--;
#endif
}

// ******************************************************************

void MEDDLY::node_headers::lastUncache(node_handle p)
{
    //
    // Cache count for p is now zero; might need to clean up
    //

    if (isDeleted(p)) {
        //
        // Must be using pessimistic
        //
#ifdef TRACK_UNREACHABLE_NODES
        stats.unreachable_nodes--;
#endif
        recycleNodeHandle(p);
    } else {
        // we were active.
        // See if we're now completely disconnected
        // and if so, tell parent to recycle node storage
        if (
#ifdef REFCOUNTS_ON
          (incoming_counts && (0==incoming_counts->get(size_t(p))))
          ||
#endif
          (is_reachable && (0==is_reachable->get(size_t(p))))
        ) {
#ifdef TRACK_UNREACHABLE_NODES
            stats.unreachable_nodes--;
#endif
            parent.deleteNode(p);
            recycleNodeHandle(p);
        }
    }
}

// ******************************************************************

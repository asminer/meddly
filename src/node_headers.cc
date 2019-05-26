
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
#include "storage/bytepack.h"


// #define DEBUG_HANDLE_FREELIST

// #define DEBUG_ADDRESS_RESIZE

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
// *                                                                *
// *               node_headers::level_array  methods               *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

MEDDLY::node_headers::level_array::level_array(memstats &ms) : MS(ms)
{
  data32 = 0;
  size = 0;
}

MEDDLY::node_headers::level_array::~level_array()
{
  MS.decMemAlloc(size * sizeof(int));
  free(data32);
}

void MEDDLY::node_headers::level_array::expand(size_t ns)
{
  if (ns <= size) return;

  int* d = (int*) realloc(data32, ns * sizeof(int));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(d + size, 0, (ns-size) * sizeof(int) );
  MS.incMemAlloc( (ns - size) * sizeof(int) );
  size = ns;
  data32 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged level array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::level_array::shrink(size_t ns)
{
  if (ns >= size) return;
  int* d = (int*) realloc(data32, ns * sizeof(int));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  MS.decMemAlloc( (size - ns) * sizeof(int) );
  size = ns;
  data32 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Reduced level array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::level_array
  ::show(output &s, size_t first, size_t last, int width) const
{
  s << "[";
  s.put(long(data32[first]), width);
  for (size_t i=first+1; i<=last; i++) {
    s.put('|');
    s.put(long(data32[i]), width);
  }
  s << "]";
}

#endif // OLD_NODE_HEADERS

// ******************************************************************
// *                                                                *
// *              node_headers::counter_array  methods              *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

MEDDLY::node_headers::counter_array::counter_array(memstats &ms) : MS(ms)
{
  data32 = 0;
  size = 0;
}

MEDDLY::node_headers::counter_array::~counter_array()
{
  MS.decMemAlloc(size * sizeof(unsigned int));
  free(data32);
}

void MEDDLY::node_headers::counter_array::expand(size_t ns)
{
  if (ns <= size) return;

  unsigned int* d = (unsigned int*) realloc(data32, ns * sizeof(unsigned int));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(d + size, 0, (ns-size) * sizeof(unsigned int) );
  MS.incMemAlloc( (ns - size) * sizeof(unsigned int) );
  size = ns;
  data32 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged counter array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::counter_array::shrink(size_t ns)
{
  if (ns >= size) return;
  unsigned int* d = (unsigned int*) realloc(data32, ns * sizeof(unsigned int));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  MS.decMemAlloc( (size - ns) * sizeof(unsigned int) );
  size = ns;
  data32 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Reduced counter array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::counter_array
  ::show(output &s, size_t first, size_t last, int width) const
{
  s << "[";
  s.put((unsigned long) data32[first], width);
  for (size_t i=first+1; i<=last; i++) {
    s.put('|');
    s.put((unsigned long) data32[i], width);
  }
  s << "]";
}

#endif

// ******************************************************************
// *                                                                *
// *              node_headers::address_array  methods              *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

MEDDLY::node_headers::address_array::address_array(memstats &ms) : MS(ms)
{
  data64 = 0;
  size = 0;
}

MEDDLY::node_headers::address_array::~address_array()
{
  MS.decMemAlloc(size * sizeof(unsigned long));
  free(data64);
}

void MEDDLY::node_headers::address_array::expand(size_t ns)
{
  if (ns <= size) return;

  unsigned long* d = (unsigned long*) realloc(data64, ns * sizeof(unsigned long));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(d + size, 0, (ns-size) * sizeof(unsigned long) );
  MS.incMemAlloc( (ns - size) * sizeof(unsigned long) );
  size = ns;
  data64 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged address array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::address_array::shrink(size_t ns)
{
  if (ns >= size) return;
  unsigned long* d = (unsigned long*) realloc(data64, ns * sizeof(unsigned long));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  MS.decMemAlloc( (size - ns) * sizeof(unsigned long) );
  size = ns;
  data64 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Reduced address array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::address_array
  ::show(output &s, size_t first, size_t last, int width) const
{
  s << "[";
  s.put(data64[first], width);
  for (size_t i=first+1; i<=last; i++) {
    s.put('|');
    s.put(data64[i], width);
  }
  s << "]";
}

#endif

// ******************************************************************
// *                                                                *
// *                node_headers::bitvector  methods                *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

MEDDLY::node_headers::bitvector::bitvector(memstats &ms) : MS(ms)
{
  data = 0;
  size = 0;
}

MEDDLY::node_headers::bitvector::~bitvector()
{
  MS.decMemAlloc(size * sizeof(bool));
  free(data);
}

void MEDDLY::node_headers::bitvector::expand(size_t ns)
{
  if (ns <= size) return;

  bool* d = (bool*) realloc(data, ns * sizeof(bool));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(d + size, 0, (ns-size) * sizeof(bool) );
  MS.incMemAlloc( (ns - size) * sizeof(bool) );
  size = ns;
  data = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged bitvector, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::bitvector::shrink(size_t ns)
{
  if (ns >= size) return;
  bool* d = (bool*) realloc(data, ns * sizeof(bool));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  MS.decMemAlloc( (size - ns) * sizeof(bool) );
  size = ns;
  data = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Reduced bitvector, new size %lu\n", size);
#endif
}

#endif

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_headers methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

const size_t START_SIZE = 512;
const size_t MAX_ADD = 65536;

inline size_t next_size(size_t s)
{
  if (0==s) return START_SIZE; 
  if (s < MAX_ADD) return s*2;
  return s+MAX_ADD;
}

inline static size_t prev_size(size_t s) 
{
  if (s < MAX_ADD) return s/2;
  return s-MAX_ADD;
}

MEDDLY::node_headers::node_headers(expert_forest &P)
  : parent(P)
{
#ifdef OLD_NODE_HEADERS
  //
  // Inltialize address array
  //
  a_size = a_min_size;
  address = (node_header *) malloc(a_size * sizeof(node_header));
  if (0 == address) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  parent.mstats.incMemAlloc(a_size * sizeof(node_header));
  memset(address, 0, a_size * sizeof(node_header));
  a_last = a_next_shrink = 0;
  for (int i=0; i<8; i++) a_unused[i] = 0;
  a_lowest_index = 8;

  //
  // hooks for later
  //
  usesCacheCounts = true;
  usesIncomingCounts = true;
#else
  addresses = new address_array(parent.mstats);
  levels = new level_array(parent.mstats);
  cache_counts = new counter_array(parent.mstats);
  incoming_counts = new counter_array(parent.mstats);
  implicit_bits = 0;
  a_last = 0;
  a_size = 0;
  a_next_shrink = START_SIZE / 4;
  for (int i=0; i<8; i++) a_unused[i] = 0;
  a_lowest_index = 8;
#endif
}

// ******************************************************************

MEDDLY::node_headers::~node_headers()
{
#ifdef OLD_NODE_HEADERS
  parent.mstats.decMemAlloc(a_size * sizeof(node_header));

  // Address array
  free(address);
#else
  delete addresses;
  delete levels;
  delete cache_counts;
  delete incoming_counts;
  delete implicit_bits;
#endif
}

// ******************************************************************

void MEDDLY::node_headers::turnOffCacheCounts()
{
#ifdef OLD_NODE_HEADERS
  usesCacheCounts = false;
  // For now - nothing else to do
#else
  delete cache_counts;
  cache_counts = 0;
#endif
}

// ******************************************************************

void MEDDLY::node_headers::turnOffIncomingCounts()
{
#ifdef OLD_NODE_HEADERS
  usesIncomingCounts = false;
  // For now - nothing else to do
#else
  delete incoming_counts;
  incoming_counts = 0;
#endif
}

// ******************************************************************

MEDDLY::node_handle MEDDLY::node_headers::getFreeNodeHandle() 
{
#ifdef OLD_NODE_HEADERS
  parent.mstats.incMemUsed(sizeof(node_header));
#endif
  node_handle found = 0;
  for (int i=a_lowest_index; i<8; i++) {
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
#ifdef DEBUG_HANDLE_FREELIST
    // address[found].setNotDeleted();
    dump_handle_info(*this, a_last+1);
#endif
    return found;
  }
  a_last++;
  if (a_last >= a_size) {
    expandHandleList();
  }
  MEDDLY_DCASSERT(a_last < a_size);

#ifdef DEBUG_HANDLE_FREELIST
  dump_handle_info(*this, a_last+1);
#endif
  return a_last;
}

// ******************************************************************

void MEDDLY::node_headers::recycleNodeHandle(node_handle p) 
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  MEDDLY_DCASSERT(0==getNodeCacheCount(p));
#ifdef OLD_NODE_HEADERS
  parent.mstats.decMemUsed(sizeof(node_header));
#endif
  deactivate(p);

  // Determine which list to add this into,
  // and add it there
  int i = bytesRequiredForDown(p) -1;
#ifdef OLD_NODE_HEADERS
  setNextOf(p, a_unused[i]);
  a_unused[i] = p;
#else
  setNextOf(size_t(p), a_unused[i]);
  a_unused[i] = size_t(p);
#endif
  a_lowest_index = MIN(a_lowest_index, (char)i);

  // if this was the last node, collapse nodes into the
  // "not yet allocated" pile.  But, we don't remove them
  // from the free list(s); we simply discard any too-large
  // ones when we pull from the free list(s).
  if (p == a_last) {
    while (a_last && isDeactivated(a_last)) {
      a_last--;
    }

    if (a_last < a_next_shrink) shrinkHandleList();

    /*
#else
    if (addresses) {
      if (addresses->resize_will_shrink(a_last+1)) {
        cleanFreeLists();
        addresses->resize(a_last+1);
      }
    }
    if (levels)           levels->resize(a_last+1);
    if (cache_counts)     cache_counts->resize(a_last+1);
    if (incoming_counts)  incoming_counts->resize(a_last+1);
    if (implicit_bits)    implicit_bits->resize(a_last+1);
#endif
    */
  }
#ifdef DEBUG_HANDLE_FREELIST
  dump_handle_info(*this, a_last+1);
#endif
}

// ******************************************************************

void MEDDLY::node_headers::swapNodes(node_handle p, node_handle q, bool swap_incounts)
{
  MEDDLY_DCASSERT(p!=q);
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  MEDDLY_DCASSERT(q>0);
  MEDDLY_DCASSERT(q<=a_last);

#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  SWAP(address[p], address[q]);

  if (!swap_incounts) {
    SWAP(address[p].incoming_count, address[q].incoming_count);
  } 
#else
  size_t sp = size_t(p);
  size_t sq = size_t(q);
  if (addresses)                          addresses->swap(sp, sq);
  if (levels)                             levels->swap(sp, sq);
  if (cache_counts)                       cache_counts->swap(sp, sq);
  if (swap_incounts && incoming_counts)   incoming_counts->swap(sp, sq);
  if (implicit_bits)                      implicit_bits->swap(sp, sq);
#endif
}

// ******************************************************************

void MEDDLY::node_headers::dumpInternal(output &s) const
{
  s << "Node headers and management:\n";
  for (int i=0; i<8; i++) {
    s << "    First " << i << "-byte unused node index: " << a_unused[i] << "\n";
  }
  int awidth = digits(a_last);
  s << "    Node# :  ";
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put(' ');
    s.put(long(p), awidth);
  }
#if OLD_NODE_HEADERS
  s << "\n    Level  : [";
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put('|');
    s.put(long(address[p].level), awidth);
  }
  s << "]\n";
  s << "\n    Offset : [";
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put('|');
    s.put(long(address[p].offset), awidth);
  }
  s << "]\n";
  s << "\n   Incount : [";
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put('|');
    s.put(long(address[p].incoming_count), awidth);
  }
  s << "]\n";
  s << "\n    Cache  : [";
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put('|');
    s.put(long(address[p].cache_count), awidth);
  }
  s << "]\n\n";
#else
  s << "\n    Level  : ";
  if (levels) {
    levels->show(s, 1, a_last, awidth);
    s << "\n";
  } else {
    s << "(null)\n";
  }
  s << "\n    Offset : ";
  if (addresses) {
    addresses->show(s, 1, a_last, awidth);
    s << "\n";
  } else {
    s << "(null)\n";
  }
  s << "\n   Incount : ";
  if (incoming_counts) {
    incoming_counts->show(s, 1, a_last, awidth);
    s << "\n";
  } else {
    s << "(null)\n";
  }
  s << "\n    Cache  : ";
  if (cache_counts) {
    cache_counts->show(s, 1, a_last, awidth);
    s << "\n";
  } else {
    s << "(null)\n";
  }
  s << "\n";
#endif
}

// ******************************************************************

void MEDDLY::node_headers::expandHandleList()
{
#ifdef OLD_NODE_HEADERS
  // increase size by 50%
  int delta = a_size / 2;
  MEDDLY_DCASSERT(delta>=0);
  node_header* new_address = (node_header*) 
    realloc(address, (a_size+delta) * sizeof(node_header));
  if (0==new_address) {
    /*
    fprintf(stderr, "Error in allocating array of size %lu at %s, line %d\n",
        (a_size+delta) * sizeof(node_header), __FILE__, __LINE__);
    */
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  address = new_address;
  parent.mstats.incMemAlloc(delta * sizeof(node_header));
  memset(address + a_size, 0, delta * sizeof(node_header));
  a_size += delta;
  a_next_shrink = a_size / 2;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged address array, new size %d\n", a_size);
#endif

#else // OLD_NODE_HEADERS

  a_next_shrink = next_size(a_next_shrink);
  a_size = next_size(a_size);

  if (addresses)        addresses->expand(a_size);
  if (levels)           levels->expand(a_size);
  if (cache_counts)     cache_counts->expand(a_size);
  if (incoming_counts)  incoming_counts->expand(a_size);
  if (implicit_bits)    implicit_bits->expand(a_size);

#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged address array, new size %lu\n", a_size);
#endif
#endif
}


// ******************************************************************

void MEDDLY::node_headers::shrinkHandleList()
{
#ifdef OLD_NODE_HEADERS

  // Determine new size
  int new_size = a_min_size;
  while (a_last >= new_size) new_size += new_size/2;
  int delta = a_size - new_size;
  if (0==delta) {
    a_next_shrink = 0;
    return;
  }

  // clean out free lists, because we're about
  // to trash memory beyond new_size.
  for (int i=0; i<8; i++) {
    //
    // clean list i
    //
    node_handle prev = 0;
    node_handle curr;
    for (curr = a_unused[i]; curr; curr=getNextOf(curr)) 
    {
      if (curr > a_last) continue;  // don't add to the list
      if (prev) {
        setNextOf(prev, curr);
      } else {
        a_unused[i] = curr;
      }
    } 
    if (prev) {
      setNextOf(prev, 0);
    } else {
      a_unused[i] = 0;
    }
  } // for i

  // shrink the array
  MEDDLY_DCASSERT(delta>=0);
  MEDDLY_DCASSERT(a_size-delta>=a_min_size);
  node_header* new_address = (node_header*) 
    realloc(address, new_size * sizeof(node_header));
  if (0==new_address) {
    /*
    fprintf(stderr, "Error in allocating array of size %lu at %s, line %d\n",
        new_size*sizeof(node_header), __FILE__, __LINE__);
    */
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  address = new_address;
  parent.mstats.decMemAlloc(delta * sizeof(node_header));
  a_size -= delta;
  a_next_shrink = a_size / 2;
  MEDDLY_DCASSERT(a_last < a_size);
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Shrank address array, new size %d\n", a_size);
#endif

#else // OLD_NODE_HEADERS

  // clean out free lists, because we're about
  // to shrink the address list which holds the pointers
  for (int i=0; i<8; i++) {
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
    } 
    if (prev) {
      setNextOf(prev, 0);
    } else {
      a_unused[i] = 0;
    }
  } // for i

  a_size = next_size(a_size);
  a_next_shrink = (a_size > START_SIZE) ? next_size(a_next_shrink) : 0;

  if (addresses)        addresses->shrink(a_size);
  if (levels)           levels->shrink(a_size);
  if (cache_counts)     cache_counts->shrink(a_size);
  if (incoming_counts)  incoming_counts->shrink(a_size);
  if (implicit_bits)    implicit_bits->shrink(a_size);

#ifdef DEBUG_ADDRESS_RESIZE
  printf("Shrank address array, new size %lu\n", a_size);
#endif

#endif
}


// ******************************************************************



// $Id$

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



// TODO: Testing

#include "compact.h"
#include "bytepack.h"

#define MERGE_AND_SPLIT_HOLES
// #define DEBUG_COMPACTION
// #define DEBUG_SLOW
// #define MEMORY_TRACE
// #define DEEP_MEMORY_TRACE




#if 0

// ******************************************************************
// *                                                                *
// *                                                                *
// *                     node_compacted methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************

/*
    Data is stored in a char* array, so we have to copy stuff off
    (since things might not be aligned properly).
    Nodes are stored as follows (numbers are offsets into the "address").
    Hole storage details follow.

      long: >=0:  The incoming count; this is an active node.
            < 0:  This is a hole; number of bytes in the hole.

      long: >=0:  Next pointer in the unique table.
            < 0:  This node is not in the unique table yet;
                  some special value is here to indicate status.
       
      uint(32b):  Size.  Number of downward pointers stored.

      uchar(8b):  Compression info.
                  High 4 bits:  dpbytes, the number of bytes used 
                                to store downward pointers.
                                Range of possible values is 0-15, 
                                but currently, only the values 1-8 are used.

                  Low 4 bits:   ixbytes, number of bytes used to store indexes.
                                Range of possible values is 0-15,
                                but currently, only the values 0-4 are used.
                                (Max node size is 2^31-1).
                                A count of 0 means the node is stored in full;
                                otherwise, the node is stored sparsely.

      uhbytes:    Unhashed header info; might be 0 bytes.
      hhbytes:    Hashed header info; might be 0 bytes.

      dpbytes:    Down pointer [0]
        :
        :
      dpbytes:    Down pointer [size-1]


      ixbytes:    Index value [0]
        :
        :
      ixbytes:    Index value [size-1]


      evbytes:    Edge value [0]
        :
        :
      evbytes:    Edge value [size-1]

      uchar(8b):  Padding length (can be 0).
                  Number of unused extra bytes, at most 255.

                  (extra bytes here, if any)

      long: > 0:  The very last slot is a non-negative long, giving
                  the forest node number.
                  (A negative value here indicates a hole.)


      The absolute smallest number of bytes for a node,
      assuming a zero size is allowed, is therefore:
        2 longs + 1 int + 1 char + 1 char + 1 long
      

      Holes are stored using a grid.
      Currently, holes require more storage than the smallest possible
      node.  If we switch to lists, we can save one long per hole,
      and then they will be very close.

      Holes are stored as follows.

        long: < 0:  -number of bytes in the hole

        long: >=0:  This is an index hole,
                    and this is an upward pointer in the grid.
              < 0:  This is a non-index hole.  Value is ignored.

        long:       down or previous pointer.
        long:       next pointer.
          :
          :         (unused bytes)
          :
        long: < 0:  -number of bytes in the hole



*/

MEDDLY::node_compacted::node_compacted()
{
  parent = 0;
}

MEDDLY::node_compacted::~node_compacted()
{
  if (parent) parent->changeStats().decMemAlloc(size);
  free(data);
}

void MEDDLY::node_compacted::init(expert_forest* p)
{
  MEDDLY_DCASSERT(0==parent);
  parent = p;
  data = 0;
  data_down = 0;
  size = 0;
  last = 0;
  max_request = 0;
  large_holes = 0;
  holes_top = 0;
  holes_bottom = 0;
  hole_slots = 0;
  evbytes = parent->edgeBytes();
  uhbytes = parent->unhashedHeaderBytes();
  hhbytes = parent->hashedHeaderBytes();
}

void MEDDLY::node_compacted::compact(bool shrink)
{
  //
  // Should we even bother?
  //
  if (0==data || 0==hole_slots) return;

  if (hole_slots <= parent->getPolicies().compact_min)  return;

  if (hole_slots <  parent->getPolicies().compact_max) {

    // If percentage of holes is below trigger, then don't compact

    if (100 * hole_slots < last * parent->getPolicies().compact_frac) return;

  }

#ifdef DEBUG_SLOW
  fprintf(stderr, "Compacting forest level\n");
#endif
#ifdef MEMORY_TRACE
  printf("Compacting\n");
#endif
#ifdef DEBUG_COMPACTION
  printf("Before compaction:\n");
  dumpInternal(stdout);
  printf("\n");
#endif

  //
  // Scan the whole array of data, copying over itself and skipping holes.
  // 
  long node = 1;  // since we leave [0] empty
  long curr = 1;

  while (node < last) {
    long L;
    // read first slot; determine if we are a hole or not.
    dataToLong<sizeof(long)>(data+node, L);
    if (L<0) {
      //
      // This is a hole, skip it
      // 
      MEDDLY_DCASSERT(
        0==memcmp(data+node, data+node-L-sizeof(long), sizeof(long))
      );
      node -= L;
      continue;
    }

    //
    // A real node, move it
    //
    long old_off = node;
    long new_off = curr;
    MEDDLY_DCASSERT(!parent->isPessimistic() || L>0);

    int size;
    dataToInt<sizeof(int)>(data_size + node, size);

    long nodebytes = size;
    nodebytes *= dpbytes(node) + ixbytes(node) + evbytes;
    nodebytes += uhbytes + hhbytes + commonHeaderBytes;
    
    // move everything up to the padding value
    if (node != curr) {
      memmove(data + curr, data + node, nodebytes);
    }
    node += nodebytes;
    curr += nodebytes;

    // skip padding
    data[curr] = 0;
    node += data[node];
    node++;

    // copy trailer
    memmove(data + curr, data + node, commonTailBytes);
    node += commonTailBytes;
    curr += commonTailBytes;
    
    // update node header
    dataToLong<sizeof(long)>(data + curr - sizeof(long), L);
    moveNodeOffset(L, old_off, new_off);
  } // while
  MEDDLY_DCASSERT(node == 1+last);
  last = curr-1;

  // set up hole pointers and such
  holes_top = holes_bottom = 0;
  hole_slots = 0;
  large_holes = 0;

  parent->changeStats().num_compactions++;

  if (shrink && size > min_size && last < size/2) {
    long new_size = size/2;
    while (new_size > min_size && new_size > last * 3) { new_size /= 2; }
    resize(new_size);
  }

#ifdef DEBUG_COMPACTION
  printf("After compaction:\n");
  dumpInternal(stdout);
  printf("\n");
#endif
}


void MEDDLY::node_compacted::dumpInternal(FILE *s) const
{
  if (0==data) return; // nothing to display
  
  fprintf(s, "Last slot used: %ld\n", last);
  fprintf(s, "large_holes: %ld\n", large_holes);
  fprintf(s, "Grid: top = %ld bottom = %ld\n", holes_top, holes_bottom);

  fprintf(s, "Data array by record: \n");
  int awidth = digits(parent->getLastNode());
  long a;
  for (a=1; a<=last; ) {
    fflush(s);
    a = dumpNode(s, a);
  } // for a
  fprintf(s, "%*ld : free slots\n", awidth, a);
  fflush(s);
  MEDDLY_DCASSERT(a == (last)+1);
}

void MEDDLY::node_compacted::dumpInternal(FILE *s, long a) const
{
  MEDDLY_DCASSERT(data);
  
  fprintf(s, "Last slot used: %ld\n", last);
  fprintf(s, "large_holes: %ld\n", large_holes);
  fprintf(s, "Grid: top = %ld bottom = %ld\n", holes_top, holes_bottom);

  dumpNode(s, a);
}


//
//
// Private
//
//

void MEDDLY::node_compacted::resize(long new_size)
{
  unsigned char* new_data = (unsigned char *) realloc(data, new_size);
  if (0 == new_data) throw error(error::INSUFFICIENT_MEMORY);
  if (0== data) new_data[0] = 0;
  if (new_size > size) {
    parent->changeStats().incMemAlloc(new_size - size);
  } else {
    parent->changeStats().decMemAlloc(size - new_size);
  }
#ifdef MEMORY_TRACE
    printf("Resized data[]. Old size: %ld, New size: %ld, Last: %ld.\n", 
      size, new_size, last
    );
#endif
  size = new_size;
  if (data != new_data) {
    // update pointers
    data = new_data;
    data_size = data + 2*sizeof(long);
    data_info = data_size + sizeof(int);
    data_down = data_info + 1 + uhbytes + hhbytes;
  }
}

long MEDDLY::node_compacted::dumpNode(FILE* s, long a) const
{
  MEDDLY_DCASSERT(data);
  if (0==a) return 0;
  int awidth = digits(parent->getLastNode());
  fprintf(s, "%*ld : [", awidth, a);

  const unsigned char* d = data + a;
  const unsigned char* stop;

  // first slot: definitely a long
  long L;
  dataToLong<sizeof(long)>(d, L);
  fprintf(s, "%ld", L);
  bool is_hole = (L<0);

  if (is_hole) {
    //
    // Hole.
    //
    stop = d - L;
    // get next 3 longs
    for (int i=3; i; i--) {
      d += sizeof(long);
      dataToLong<sizeof(long)>(d, L);
      fprintf(s, "|%ld", L);
    }
    fprintf(s, "| ... ");
  } else {
    //
    // Node.
    //
    
    // get next: long
    d += sizeof(long);
    dataToLong<sizeof(long)>(d, L);
    fprintf(s, "|n=%ld", L);
    d += sizeof(long);

    // get size: int
    int size;
    dataToInt<sizeof(int)>(d, size);
    fprintf(s, "|s=%d", size);
    d += sizeof(int);

    // get compression info
    int dpb = dpbytesOf(d[0]);
    MEDDLY_DCASSERT(dpbytes(a) == dpb);
    int ixb = ixbytesOf(d[0]);
    MEDDLY_DCASSERT(ixbytes(a) == ixb);
    d++;
    fprintf(s, "|dpb=%d, ixb=%d", dpb, ixb);

    // show unhashed header
    if (uhbytes) {
      fprintf(s, "|uh");
      for (int i=0; i<uhbytes; i++) {
        fprintf(s, " %02x", d[0]);
        d++;
      }
    }

    // show hashed header
    if (hhbytes) {
      fprintf(s, "|hh");
      for (int i=0; i<hhbytes; i++) {
        fprintf(s, " %02x", d[0]);
        d++;
      }
    }

    // show next chunk: downward pointers
    for (int i=0; i<size; i++) {
      switch (dpb) {
        case 1: dataToDown<1>(d, L);  break;
        case 2: dataToDown<2>(d, L);  break;
        case 3: dataToDown<3>(d, L);  break;
        case 4: dataToDown<4>(d, L);  break;
        case 5: dataToDown<5>(d, L);  break;
        case 6: dataToDown<6>(d, L);  break;
        case 7: dataToDown<7>(d, L);  break;
        case 8: dataToDown<8>(d, L);  break;
        default:
          throw error(error::MISCELLANEOUS);
      }
      fprintf(s, "|d %ld", L);
      d += dpb;
    }

    // show next chunk: indexes
    if (ixb) for (int i=0; i<size; i++) {
      int I;
      switch (ixb) {
        case 1: dataToInt<1>(d, I); break;
        case 2: dataToInt<2>(d, I); break;
        case 3: dataToInt<3>(d, I); break;
        case 4: dataToInt<4>(d, I); break;
        default:
          throw error(error::MISCELLANEOUS);
      }
      fprintf(s, "|i %d", I);
      d += ixb;
    }

    // show next chunk: edge values (raw)
    if (evbytes) for (int i=0; i<size; i++) {
      fprintf(s, "|e");
      for (int b=0; b<evbytes; b++) {
        fprintf(s, " %02x", d[0]);
        d++;
      }
    }

    // get the padding
    unsigned int padding = d[0];
    d++;
    fprintf(s, "| ... %d unused ... ", padding);
    stop = d + padding + sizeof(long);
  }
  // show the last, long slot
  dataToLong<sizeof(long)>(stop-sizeof(long), L);
  fprintf(s, "|%ld]\n", L);
  return stop - (data+a);
}


#endif  // if 0



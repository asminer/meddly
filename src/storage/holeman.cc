
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

#include "holeman.h"

// ******************************************************************
// *                                                                *
// *                        holeman  methods                        *
// *                                                                *
// ******************************************************************

MEDDLY::holeman::holeman(int smallestHole, node_storage* p)
{
  parent = p;
  num_holes = 0;
  max_holes = 0;
  hole_slots = 0;
  max_hole_slots = 0;
  fragment_slots = 0;
  max_fragment_slots = 0;
  last_slot = 0;
  MEDDLY_DCASSERT(parent);
  smallest = MAX(smallestHole, parent->smallestNode());
  data = 0;
  size = 0;
  parent->updateData(data);
}

// ******************************************************************

MEDDLY::holeman::~holeman()
{
  // don't delete our parent
  decMemAlloc(size*sizeof(node_handle));
  free(data);
}

// ******************************************************************

MEDDLY::node_address MEDDLY::holeman::chunkAfterHole(node_address addr) const
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(data[addr]<0);
  MEDDLY_DCASSERT(data[addr-data[addr]-1] == data[addr]);
  return addr - data[addr];
}

// ******************************************************************

void MEDDLY::holeman
::reportStats(FILE* s, const char* pad, unsigned flags) const
{
  if (flags & expert_forest::HOLE_MANAGER_STATS) {
    fprintf(s, "%s    %d holes currently\n", pad, num_holes);
    fprintf(s, "%s    %d max holes seen\n", pad, max_holes);

    unsigned long holemem = hole_slots * sizeof(node_handle);
    fprintf(s, "%s    ", pad);
    fprintmem(s, holemem, flags & expert_forest::HUMAN_READABLE_MEMORY);
    fprintf(s, " wasted in holes\n");

    holemem = max_hole_slots * sizeof(node_handle);
    fprintf(s, "%s    ", pad);
    fprintmem(s, holemem, flags & expert_forest::HUMAN_READABLE_MEMORY);
    fprintf(s, " max in holes\n");

    unsigned long fragmem = fragment_slots * sizeof(node_handle);
    fprintf(s, "%s    ", pad);
    fprintmem(s, fragmem, flags & expert_forest::HUMAN_READABLE_MEMORY);
    fprintf(s, " wasted in fragments\n");

    fragmem = max_fragment_slots * sizeof(node_handle);
    fprintf(s, "%s    ", pad);
    fprintmem(s, fragmem, flags & expert_forest::HUMAN_READABLE_MEMORY);
    fprintf(s, " max in fragments\n");
  }
}

// ******************************************************************

void MEDDLY::holeman::clearHolesAndShrink(node_address new_last, bool shrink)
{
  last_slot = new_last;
  num_holes = 0;
  hole_slots = 0;
  fragment_slots = 0;

  if (shrink && size > min_size && last_slot < size/2) {
    node_address new_size = size/2;
    while (new_size > min_size && new_size > last_slot * 3) { new_size /= 2; }
    resize(new_size);
  }
}

// ******************************************************************

void MEDDLY::holeman::resize(node_address new_slots)
{
  node_handle* new_data 
    = (node_handle*) realloc(data, new_slots * sizeof(node_handle));

  if (0 == new_data) throw error(error::INSUFFICIENT_MEMORY);
  if (0 == data) new_data[0] = 0;
  if (new_slots > size) {
    incMemAlloc((new_slots - size) * sizeof(node_handle));
  } else {
    decMemAlloc((size - new_slots) * sizeof(node_handle));
  }
#ifdef MEMORY_TRACE
    printf("Resized data[]. Old size: %ld, New size: %ld, Last: %ld.\n", 
      size, new_slots, last
    );
#endif
  size = new_slots;
  if (data != new_data) {
    // update pointers
    data = new_data;
    parent->updateData(data);
  }
}


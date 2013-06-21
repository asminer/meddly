
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

#include "unique_table.h"

#include <limits.h>

MEDDLY::unique_table::unique_table(expert_forest* ef)
{
  parent = ef;
  size = minSize;
  num_entries = 0;
  table = (node_handle*) malloc(sizeof(node_handle) * size);
  if (0==table) throw error(error::INSUFFICIENT_MEMORY);
  for (unsigned i=0; i<size; i++) table[i] = 0;
  next_expand = 2*size;
  next_shrink = 0;
}

MEDDLY::unique_table::~unique_table()
{
  free(table);
}

void MEDDLY::unique_table::show(FILE* s) const
{
  fprintf(s, "Unique table:\n");
  for (unsigned i=0; i<size; i++) {
    fprintf(s, "[%d] : ", i);
    for (int index = table[i]; index; index = parent->getNext(index)) {
      fprintf(s, "%d ", index);
    }
    fprintf(s, "\n");
  } 
  fflush(s);
}

//
// Helpers (private)
//

MEDDLY::node_handle MEDDLY::unique_table::convertToList()
{
  /*
  printf("Converting to list\n");
  show(stdout);
*/
  node_handle front = 0;
  node_handle next = 0;
  node_handle curr = 0;
  for (unsigned i = 0; i < size; i++) {
    for (curr = table[i]; curr; curr = next) {
      // add to head of list
      next = parent->getNext(curr);
      parent->setNext(curr, front);
      front = curr;
    }
    table[i] = 0;
  } // for i
  num_entries = 0;
  return front;
}

void MEDDLY::unique_table::buildFromList(node_handle front)
{
  node_handle next = 0;
  for ( ; front; front = next) {
    next = parent->getNext(front);
    unsigned h = parent->hash(front) % size;
    MEDDLY_CHECK_RANGE(0, h, size);
    parent->setNext(front, table[h]);
    table[h] = front;
    num_entries++;
  }
}

void MEDDLY::unique_table::expand()
{
  MEDDLY_DCASSERT(size < maxSize);
#ifdef DEBUG_SLOW
  fprintf(stderr, "Enlarging unique table (current size: %d)\n", size);
#endif
  node_handle ptr = convertToList();
  // length will the same as num_entries previously
  unsigned newSize = size * 2;
  node_handle *temp = (node_handle*) realloc(table, sizeof(node_handle) * newSize);
  if (0==temp) throw error(error::INSUFFICIENT_MEMORY);
  table = temp;
  for (unsigned i = newSize-1; i>=size; i--) table[i] = 0;
  next_shrink = size;
  size = newSize;
  if (size >= maxSize) next_expand = UINT_MAX;
  else                 next_expand = size * 2;
  buildFromList(ptr);
}

void MEDDLY::unique_table::shrink()
{
  MEDDLY_DCASSERT(size > minSize);
#ifdef DEBUG_SLOW
  fprintf(stderr, "Shrinking unique table (current size: %d)\n", size);
#endif
  node_handle ptr = convertToList();
  // length will the same as num_entries previously
  unsigned newSize = size / 2;
  node_handle *temp = (node_handle*) realloc(table, sizeof(node_handle) * newSize);
  if (0==temp) throw error(error::INSUFFICIENT_MEMORY);
  table = temp;
  next_expand = size;
  size = newSize;
  if (size <= minSize) next_shrink = 0;
  else                 next_shrink = size / 2;
  buildFromList(ptr);
}


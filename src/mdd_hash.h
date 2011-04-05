
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



/*
	MDD Hash table template class.

	Operations are performed on a "node manager" object, where nodes are indexed
  by handles.  Several hash tables can share a single node manager if desired.

  This class is used meant to be used by node managers only since this does
  not discard any nodes unless explicitly asked to do so by calling remove().

	The node manager must provide the following methods 
	(preferably inlined for speed):

	int getNull();                   // which integer handle to use for NULL.

  int getNext(int h);              // next node in this chain

  void setNext(int h);             // set the "next node" field for h

	unsigned hash(int h);             // compute hash value

	bool equals(int h1, int h2);      // is h1 == h2?

	for debugging:
	show(FILE *s, int h);
*/

#ifndef MDD_HASH_H
#define MDD_HASH_H

#include <sys/time.h>
#include <sys/resource.h>
#include "defines.h"

template <typename MANAGER>
class mdd_hash_table {

  protected:
    // size grows by a factor of 2
    unsigned size;
    unsigned num_entries;
    int* table;
    MANAGER *nodes;

    bool out_of_mem;      // no memory system memory available, can't expand

  public:
    mdd_hash_table(MANAGER *n) {
      num_entries = 0;
      size = getMinSize();
      nodes = n;
      table = (int*) malloc(sizeof(int) * size);
      if (NULL == table) outOfMemory();
      out_of_mem = false;
      int *end = table + size;
      for (int* curr = table; curr < end; curr++) {
        *curr = nodes->getNull();
      }
    }

    ~mdd_hash_table() {
      free(table);
    }

    inline unsigned getSize() const { return size; }
    inline unsigned getMaxSize() const { return 1073741824; }
    inline unsigned getMinSize() const { return 8; }
    inline unsigned getEntriesCount() const { return num_entries; }

    void show(FILE *s) const {
      fprintf(s, "%s :\n", __func__);
      for (unsigned i = 0; i < size; i++) {
        fprintf(s, "[%d] : ", i);
        for (int index = table[i];
            index != nodes->getNull();
            index = nodes->getNext(index)) {
          nodes->show(s, index);
          fprintf(s, " ");
        }
        fprintf(s, "\n");
      }
      fflush(s);
    }

    /// Empty the hash table into a list; returns the list.
    inline int convertToList(unsigned& listlength) {
      // fprintf(stdout, "mdd_hash: %s\n", __func__);
      listlength = 0;
      int front = nodes->getNull();
      int next = nodes->getNull();
      int curr = nodes->getNull();
      for (unsigned i = 0; i < size; i++) {
        for (curr = table[i]; curr != nodes->getNull(); curr = next) {
          // add to head of list
          next = nodes->getNext(curr);
          nodes->setNext(curr, front);
          front = curr;
          listlength++;
        }
        table[i] = nodes->getNull();
      } // for i
      DCASSERT(listlength == num_entries);
      num_entries = 0;
      return front;
    }

    /// A series of inserts; doesn't check for duplicates or expand.
    inline void buildFromList(int front) {
      // fprintf(stdout, "mdd_hash: %s\n", __func__);
      int next = nodes->getNull();
      unsigned h = 0;
      for ( ; front != nodes->getNull(); front = next) {
        next = nodes->getNext(front);
        h = hash(front);
        CHECK_RANGE(0, h, size);
        nodes->setNext(front, table[h]);
        table[h] = front;
        num_entries++;
      }
    }

    inline unsigned hash(int h) {
      return nodes->hash(h) % size;
    }

    /// Expand the hash table (if possible)
    inline void expand() {
      if (size >= getMaxSize() || out_of_mem) return;
#ifdef DEBUG_MDD_HASH_EXPAND_H
      fprintf(stderr, "%s: Trying to enlarge unique table...\n", __func__);
      // fprintf(stderr, "%s: Old table:\n", __func__);
      // show(stderr);
      fflush(stderr);
#endif
      unsigned length = 0;
      int ptr = convertToList(length);
      // length will the same as num_entries previously
      int newSize = size << 1;
#ifdef DEBUG_MDD_HASH_EXPAND_H
      fprintf(stderr, "%s: Enlarging table...\n", __func__);
      fprintf(stderr, "\tnum_entries = %d, new size = %d\n", length, newSize);
      fflush(stderr);
#endif
      int *temp = (int*) realloc(table, sizeof(int) * newSize);
      if (NULL == temp) {
        fprintf(stderr, "Memory allocation error in unique table expand.\n");
        out_of_mem = true;
        buildFromList(ptr);
        return;
      }
      table = temp;
      size = newSize;
      int *last = table + size;
      for (int *curr = table ; curr < last; ++curr) {
        *curr = nodes->getNull();
      }
      buildFromList(ptr);
#ifdef DEBUG_MDD_HASH_EXPAND_H
      // fprintf(stderr, "%s: New table:\n", __func__);
      // show(stderr);
      fprintf(stderr, "%s:    finished enlarging unique table\n", __func__);
      fflush(stderr);
#endif
    }


    inline void shrink() {
      if (size <= getMinSize()) return;
#ifdef DEBUG_MDD_HASH_H
      // fprintf(stderr, "Old table:\n");
      // show(stderr);
      fprintf(stderr, "Shrinking table from size %d", size);
      fflush(stderr);
#endif
      unsigned length = 0;
      int front = convertToList(length);
      int newSize = size >> 1;
#ifdef DEBUG_MDD_HASH_H
      fprintf(stderr, " to size %d\n", newSize);
      fflush(stderr);
#endif
      int *temp = (int*) realloc(table, sizeof(int) * newSize);
      if (NULL == temp) {
        fprintf(stderr, "Memory allocation error in unique table shrink.\n");
        exit(1);
      } else {
        table = temp;
        size = newSize;
      }
      int *last = table + size;
      for (int *curr = table; curr < last; ++curr) {
        *curr = nodes->getNull();
      }
      out_of_mem = false;
      buildFromList(front);
#ifdef DEBUG_MDD_HASH_H
      fprintf(stderr, "   finished shrinking\n");
      // fprintf(stderr, "New table:\n");
      // show(stderr);
      fflush(stderr);
#endif
    }

    /** If table contains key, move it to the front of the list.
      Otherwise, do nothing.
      Returns index of the item if found, nodes->getNull() otherwise.
     */
    int find(int key) {
#ifdef DEBUG_MDD_HASH_H
      fprintf(stderr, "%s: key = %d, size = %d ", __func__, key, size);
      nodes->show(stderr, key);
      fprintf(stderr, "\n");
#endif

      unsigned h = hash(key);

#ifdef DEBUG_MDD_HASH_H
      fprintf(stderr, "%s: h = %d\n", __func__, h);
#endif
      CHECK_RANGE(0, h, size);
      int null = nodes->getNull();
      if (table[h] == null) return null;
      int parent = null;
      for (int ptr = table[h];
          ptr != null;
          ptr = nodes->getNext(ptr)) {
        if (nodes->equals(key, ptr)) {
          if (ptr != table[h]) {
            // not at front; move it to front
            DCASSERT(parent != null);
            nodes->setNext(parent, nodes->getNext(ptr));
            nodes->setNext(ptr, table[h]);
            table[h] = ptr;
          }
#ifdef DEBUG_MDD_HASH_H
          fprintf(stderr, "%s: found match = %d\n", __func__, ptr);
#endif
          return ptr;
        }
        parent = ptr;
      } // for ptr
#ifdef DEBUG_MDD_HASH_H
      fprintf(stderr, "%s: did not find match\n", __func__);
#endif
      return null;
    }


    /** If table contains key, remove it and return it.
      I.e., the exact key.
      Otherwise, return nodes->getNull().
     */
    inline int remove(int key) {
      unsigned h = hash(key);
      CHECK_RANGE(0, h, size);
      int parent = nodes->getNull();
      int ptr = table[h];
      for ( ; ptr != nodes->getNull(); ptr = nodes->getNext(ptr)) {
        if (key == ptr) {
          if (ptr != table[h]) {
            // remove from middle
            nodes->setNext(parent, nodes->getNext(ptr));
          } else {
            // remove from head
            table[h] = nodes->getNext(ptr);
          }
          num_entries--;
          DCASSERT(table[h] != 0);
          return ptr;
        }
        parent = ptr;
      }
      DCASSERT(table[h] != 0);
      return nodes->getNull();
    }


    /** If table contains key, move it to the front of the list.
      Otherwise, add key to the front of the list.
      Returns the (new) front of the list.
     */
    inline int insert(int key) {
      DCASSERT(key != 0);
      DCASSERT(key >= 1);
      DCASSERT(find(key) == nodes->getNull());
      if (num_entries >= (size << 1)) expand();
      unsigned h = hash(key);
      CHECK_RANGE(0, h, size);
      // insert at head of list
      nodes->setNext(key, table[h]);
      table[h] = key;
      num_entries++;
      return table[h];
    }


    int replace(int key) {
      unsigned h = hash(key);
      CHECK_RANGE(0, h, size);

#ifdef DEBUG_MDD_HASH_H
      fprintf(stderr, "%s: h = %d\n", __func__, h);
#endif

      // search this bucket
      int null = nodes->getNull();
      if (table[h] != null) {
        int parent = null;
        for (int ptr = table[h];
            ptr != null;
            ptr = nodes->getNext(ptr)) {
          if (nodes->equals(key, ptr)) {
#ifdef DEBUG_MDD_HASH_H
            fprintf(stderr, "%s: found match = %d\n", __func__, ptr);
#endif
            return ptr;
          }
          parent = ptr;
        } // for ptr
      }

      // not found; now insert
      if (num_entries >= (size << 1)) {
        expand();
        h = hash(key);
      }
      nodes->setNext(key, table[h]);
      table[h] = key;
      ++num_entries;
#ifdef DEBUG_MDD_HASH_H
      fprintf(stderr, "%s: did not find match, inserted.\n", __func__);
#endif
      return key;
    }
};

#endif

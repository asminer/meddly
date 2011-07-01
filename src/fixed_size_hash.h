
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
	Fixed-size Hash table template class.

  This is similar to hash_table described in hash.h except that it cannot
  grow or shink. It stays at a size indicated by the initial size.

	Operations are performed on a "node manager" object, where
	nodes are indexed by handles.  Several hash tables can
	share a single node manager if desired.
	The node manager must provide the following methods 
	(preferably inlined for speed):

	int getNull();                   // which integer handle to use for NULL.

  int getNext(int h);              // next node in this chain

  void setNext(int h);             // set the "next node" field for h

	unsigned hash(int h, unsigned M); // compute hash value

	bool equals(int h1, int h2);      // is h1 == h2?

  bool isStale(int h);             // is h a stale (invalid/unwanted) node

  void uncacheNode(int h);         // node h has been removed from the cache,
                                    // perform clean up of node (decrement
                                    // cache count, release node, etc.)

  void cacheNode(int h);            // node h has been added to the cache,
                                    // perform updates to external data-
                                    // structures (increment cache count etc.)

	for debugging:
	show(FILE *s, int h);
*/

#ifndef FS_HASH_H
#define FS_HASH_H

#include <sys/time.h>
#include <sys/resource.h>
#include "defines.h"

#define CHECK_STALES_BEFORE_EQUAL 1
#define ENABLE_STALE_CHECKING 1
#define AGGRESIVE_GC 0
#define START_AT_MAX_SIZE 0

#define EXP_FACTOR 1.50       // new_size = old_size * 1.50

template <typename MANAGER>
class fixed_size_hash_table {

  protected:
    unsigned  table_size;
    unsigned  max_table_size;
    unsigned  min_table_size;
    float     max_load_factor;
    float     min_load_factor;
    unsigned  max_load_size;
    unsigned  min_load_size;
    unsigned  num_entries;
    int*      table;
    MANAGER*  nodes;

  public:
    // table_size = # of hash_table buckets
    fixed_size_hash_table(MANAGER *n, unsigned max_table_size = 16777216u,
        float max_load_factor = 0.7) {
      nodes = n;
      num_entries = 0;
      min_table_size = 1024u;
      this->max_table_size = max_table_size;
#if START_AT_MAX_SIZE
      table_size = max_table_size;
#else
      table_size = min_table_size;
#endif
      this->max_load_factor = max_load_factor;
      min_load_factor = 0.25;
      max_load_size = unsigned(table_size * max_load_factor);
      min_load_size = unsigned(table_size * min_load_factor);
      table = (int*) malloc(sizeof(int) * getSize());
      if (NULL == table) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);
      int *end = table + getSize();
      for (int* curr = table; curr < end; curr++) {
        *curr = nodes->getNull();
      }
    }

    ~fixed_size_hash_table() {
      // Do not delete nodes -- why? Because it performs needless
      // uncacheNode()'s incase you just want to discard this table
#if 0
      // Deleting nodes!
      for (unsigned i = 0; i < getSize(); i++) {
        int next = nodes->getNull();
        for (int curr = table[i]; curr != nodes->getNull(); curr = next) {
          next = nodes->getNext(curr);
          // remove from cache
          nodes->uncacheNode(curr);
        }
      } // for i
#endif
      free(table);
    }

    inline unsigned getMemoryUsage() const {
      return table_size * sizeof(int);
    }

    inline unsigned getSize() const {
      return table_size;
    }

    inline unsigned getEntriesCount() const { return num_entries; }
    /*
       inline unsigned get_max_entries_count() const { return max_entries; }
       inline unsigned get_max_chain_count() const { return maxchain; }
     */

    void show(FILE *s, bool verbose = false) const {
      char filler[] = "\t";
      fprintf(s, "%sNumber of slots:        %d\n", filler, getSize());
      fprintf(s, "%sNumber of entries:      %d\n", filler, getEntriesCount());
      if (verbose) {
        fprintf(s, "Hash table entries:\n");
        for (unsigned i = 0; i < getSize(); i++) {
          fprintf(s, "[%d] : ", i);
          if (table[i] != nodes->getNull()) {
            nodes->show(s, table[i]);
          }
          fprintf(s, "\n");
        }
      }
      fflush(s);
    }


    /// Empty the hash table into a list; returns the list.
    inline int convertToList(unsigned& listlength) {
#if AGGRESIVE_GC
      unsigned stale_count = 0;
      do {
        stale_count = 0;
        for (unsigned i = 0; i < getSize(); i++) {
          if (table[i] != nodes->getNull() && nodes->isStale(table[i])) {
            // remove from cache
            nodes->uncacheNode(table[i]);
            table[i] = nodes->getNull();
            stale_count++;
          }
        }
      } while (stale_count > 0);
#endif
      listlength = 0;
      int front = nodes->getNull();
      for (unsigned i = 0; i < getSize(); i++) {
        if (table[i] != nodes->getNull()) {
          DCASSERT(nodes->getNext(table[i]) == nodes->getNull());
#if ENABLE_STALE_CHECKING
          if (nodes->isStale(table[i])) {
            // remove from cache
            nodes->uncacheNode(table[i]);
          } else {
            // add to head of list
            nodes->setNext(table[i], front);
            front = table[i];
            listlength++;
          }
#else
          // add to head of list
          nodes->setNext(table[i], front);
          front = table[i];
          listlength++;
#endif
          table[i] = nodes->getNull();
        }
      } // for i
      num_entries = 0;
      // note that there may still be stale entries in this list
      return front;
    }

    /// A series of inserts; doesn't check for duplicates or expand.
    inline void buildFromList(int front) {
      unsigned h;
      DCASSERT(num_entries == 0);
      while (front != nodes->getNull()) {
        h = nodes->hash(front, getSize());
        CHECK_RANGE(0, h, getSize());
        if (table[h] != nodes->getNull()) {
          // discard the old entry
          nodes->uncacheNode(table[h]);
          num_entries--;
        }
        // add
        table[h] = front;
        front = nodes->getNext(front);
        nodes->setNext(table[h], nodes->getNull());
        num_entries++;
      }
    }

    /// Expand the hash table size and rehash the entries
    inline void expand() {
      // convert old entries to a list
      unsigned length = 0;
      int front = convertToList(length);

      // expand hash table
#ifdef EXP_FACTOR
      table_size = int(table_size * float(EXP_FACTOR));
#else
      table_size *= 2;
#endif
      if (table_size > max_table_size) table_size = max_table_size;
      max_load_size = int(table_size * max_load_factor);
      min_load_size = int(table_size * min_load_factor);
      table = (int *) realloc(table, table_size * sizeof(int));
      if (NULL == table) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);

      for (unsigned i = 0; i < table_size; i++) {
        table[i] = nodes->getNull();
      }

      // build from list
      buildFromList(front);
    }

    /// Shrink the hash table and rehash entries
    inline void shrink() {
      // convert old entries to a list
      unsigned length = 0;
      int front = convertToList(length);

      // shrink hash table
      table_size /= 2;
      max_load_size = int(table_size * max_load_factor);
      min_load_size = int(table_size * min_load_factor);
      table = (int *) realloc(table, table_size * sizeof(int));
      if (NULL == table) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);

      for (unsigned i = 0; i < getSize(); i++) {
        table[i] = nodes->getNull();
      }

      // build from list
      buildFromList(front);
    }


    /// Go through the hash table and remove all entries
    inline void clear() {
      for (unsigned i = 0; i < getSize(); i++) {
        if (table[i] != nodes->getNull()) {
          nodes->uncacheNode(table[i]);
          table[i] = nodes->getNull();
        }
      }

      num_entries = 0;
#if START_AT_MAX_SIZE
      table_size = max_table_size;
#else
      table_size = min_table_size;
#endif
      max_load_size = unsigned(table_size * max_load_factor);
      min_load_size = unsigned(table_size * min_load_factor);

      table = (int*) realloc(table, sizeof(int) * getSize());
      if (NULL == table) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);
      int *end = table + getSize();
      for (int* curr = table; curr < end; curr++) {
        *curr = nodes->getNull();
      }
    }


    /// Go through the hash table and remove all stale entries
    inline unsigned removeStaleEntries(bool shrinkable = true) {
      // return 0;
      unsigned stale_count = 0;
      for (unsigned i = 0; i < getSize(); i++) {
        if (table[i] != nodes->getNull()) {
          if (nodes->isStale(table[i])) {
            stale_count++;
            nodes->uncacheNode(table[i]);
            table[i] = nodes->getNull();
          }
        }
      } // for i
      num_entries -= stale_count;
      return stale_count;
    }

    /** If table contains key, return the index of the item found
     *  Otherwise, return getNull()
     */
    int find(int key) {
#ifdef DEBUG_HASH_H
      printf("%s: key = %d, size = %d ", __func__, key, getSize());
      nodes->show(stdout, key);
      printf("\n");
#endif

      unsigned h = nodes->hash(key, getSize());

#ifdef DEBUG_HASH_H
      printf("%s: key = %d, h = %d, size = %d ", __func__, key, h, getSize());
      nodes->show(stdout, key);
      printf("\n");
#endif
      CHECK_RANGE(0, h, getSize());

#if CHECK_STALES_BEFORE_EQUAL
      if (table[h] != nodes->getNull()) {
        if (nodes->isStale(table[h])) {
          nodes->uncacheNode(table[h]);
          table[h] = nodes->getNull();
          num_entries--;
        } else if (nodes->equals(key, table[h])) {
          return table[h];
        }
      }
#else
      if (table[h] != nodes->getNull()) {
        if (nodes->equals(key, table[h])) {
          return table[h];
#if ENABLE_STALE_CHECKING
        } else if (nodes->isStale(table[h])) {
          nodes->uncacheNode(table[h]);
          table[h] = nodes->getNull();
          num_entries--;
#endif
        }
      }
#endif
      return nodes->getNull();
    }


    /** If table contains key, remove it and return it; otherwise return
        getNull(). Note that the key in this case represents a logical
        node address -- the contents of the node are not investigated
        to determine equality.
    */
    inline int remove(int key) {
#ifdef DEVELOPMENT_CODE
      unsigned h = nodes->hash(key, getSize());
      CHECK_RANGE(0, h, getSize());
      int& head = table[h];
#else
      int& head = table[nodes->hash(key, getSize())];
#endif

      if (head != nodes->getNull() && head == key) {
        nodes->uncacheNode(head);
        head = nodes->getNull();
        num_entries--;
        return key;
      }

      return nodes->getNull();
    }


    /** If table contains key, do nothing.
     *  Otherwise, add it to the front of the list.
     *  Returns the front of the list.
     */
    inline int insert(int key) {
      DCASSERT(find(key) == nodes->getNull());
#if START_AT_MAX_SIZE
#else
      if (num_entries >= max_load_size) {
        if (table_size < max_table_size) expand();
      }
      else if (num_entries <= min_load_size) {
        if (table_size > min_table_size) shrink();
      }
#endif

      unsigned h = nodes->hash(key, getSize());
      CHECK_RANGE(0, h, getSize());
      // insert at head of list
      if (table[h] != nodes->getNull()) {
        // regardless of this nodes state (could be stale),
        // remove it from the cache
        DCASSERT(nodes->getNext(table[h]) == nodes->getNull());
        nodes->uncacheNode(table[h]);
        table[h] = nodes->getNull();
        num_entries--;
      }
      nodes->setNext(key, table[h]);
      table[h] = key;
      num_entries++;
      // nodes->cacheNode(table[h]);
      return table[h];
    }
};

#endif

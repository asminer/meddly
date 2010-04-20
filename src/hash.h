
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
	Hash table template class.

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

#ifndef HASH_H
#define HASH_H

#include <sys/time.h>
#include <sys/resource.h>
#include "defines.h"

#define CHAIN_LENGTH_INFO
#ifdef CHAIN_LENGTH_INFO
#include <vector>
#endif

#define DISCARD_STALES 1            // 1: check for stale entries
#define STALE_BEFORE_EQUAL 1        // 1: check for stales before equality
                                    //    this MUST be 1 for pessimistic caches
#define AGGRESIVE_GC 0

#define MAX_DEPTH 2
#define EXPANSION_FACTOR 2

//#define DEBUG_SIZE

template <typename MANAGER>
class hash_table {

  protected:
    unsigned  size;         // size growth: (* EXPANSION_FACTOR) or (* 2)
    unsigned  sizeBy4;
    unsigned  nEntries;
    unsigned  maxTableSize;
    unsigned  minTableSize;
    int*      table;
#ifdef MAX_DEPTH
    int*      depth;
    int       discarded;
#endif
    MANAGER*  nodes;
    bool      noMem;        // no memory system memory available, can't expand

  public:
    hash_table(MANAGER *n, unsigned maxTableSize = 16777216u) {
      nEntries = 0;
      minTableSize = 1024u;

#ifdef MAX_DEPTH
      this->maxTableSize = maxTableSize/int(MAX_DEPTH);
      discarded = 0;
#else
      this->maxTableSize = maxTableSize;
#endif

      setSize(getMinSize());
      nodes = n;
      table = (int*) malloc(sizeof(int) * getSize());
      assert(NULL != table);

#ifdef MAX_DEPTH
      depth = (int*) malloc(sizeof(int) * getSize());
      assert(NULL != depth);
      memset(depth, 0, sizeof(int) * getSize());
#endif

      noMem = false;
      int *end = table + getSize();
      for (int* curr = table; curr < end; curr++) {
        *curr = nodes->getNull();
      }

#ifdef ACCESS_STATS
      pings = 0;
      hits = 0;
#endif
    }

    ~hash_table() {
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
#ifdef MAX_DEPTH
      free(depth);
#endif
    }

    inline void setSize(unsigned sz) {
      DCASSERT(sz > 0);
      size = sz;
#ifdef MAX_DEPTH
      sizeBy4 = int(sz * int(MAX_DEPTH) * 0.7f);
#else
      sizeBy4 = sz * 2;
#endif
#ifdef DEBUG_SIZE
      printf("max size: %d, min size: %d, current size: %d, sizeby4: %d\n",
          getMaxSize(), getMinSize(), getSize(), getSizeBy4());
#endif
    }
    inline unsigned getSizeBy4() const { return sizeBy4; }
    inline unsigned getSize() const { return size; }
    inline unsigned getMemoryUsage() const { return size * sizeof(int); }
    inline unsigned getMaxSize() const { return maxTableSize;}
    inline unsigned getMinSize() const { return minTableSize; }
    inline unsigned getEntriesCount() const { return nEntries; }
    /*
       inline unsigned get_max_entries_count() const { return max_entries; }
       inline unsigned get_max_chain_count() const { return maxchain; }
       */

    void show(FILE *s, bool verbose = false) const {
      char filler[] = "\t";
      fprintf(s, "%sNumber of slots:        %d\n", filler, getSize());
      fprintf(s, "%sNumber of entries:      %d\n", filler, getEntriesCount());
#ifdef CHAIN_LENGTH_INFO
      fprintf(s, "%sChains:\n", filler); 
      fprintf(s, "%s   Depth    Count\n", filler); 
      // compute chain distribution
      std::vector<unsigned> chains;
      for (unsigned i = 0; i < getSize(); i++) {
        unsigned depth = 0;
        for (int index = table[i];
            index != nodes->getNull();
            index = nodes->getNext(index)) {
          depth++;
        }
        if (chains.size() <= depth) chains.resize(depth + 1);
        chains[depth]++;
      }
      for (unsigned i = 0u; i < chains.size(); i++)
      {
        fprintf(s, "%s   %d       %d\n", filler, i, chains[i]); 
      }
#endif
      if (verbose) {
        fprintf(s, "Hash table entries:\n");
        for (unsigned i = 0; i < getSize(); i++) {
          fprintf(s, "[%d] : ", i);
          for (int index = table[i];
              index != nodes->getNull();
              index = nodes->getNext(index)) {
            nodes->show(s, index);
            fprintf(s, " ");
          }
          fprintf(s, "\n");
        }
      }
      fflush(s);
    }

    /// Empty the hash table into a list; returns the list.
    // there could still be stale nodes in this list...
    // call removeStaleEntries() before this to remove as many
    // stale nodes as possible
    inline int convertToList(unsigned& listlength) {
      // printf("hash: %s\n", __func__);
      listlength = 0;
      int front = nodes->getNull();
      for (unsigned i = 0; i < getSize(); i++) {
        int next = nodes->getNull();
        for (int curr = table[i]; curr != nodes->getNull(); curr = next) {
          next = nodes->getNext(curr);
#if DISCARD_STALES
          if (nodes->isStale(curr)) {
            // remove from cache
            nodes->uncacheNode(curr);
          } else {
            // add to head of list
            nodes->setNext(curr, front);
            front = curr;
            listlength++;
          }
#else
          // add to head of list
          nodes->setNext(curr, front);
          front = curr;
          listlength++;
#endif
        }
        table[i] = nodes->getNull();
      } // for i
      nEntries = 0;
      return front;
    }

    /// A series of inserts; doesn't check for duplicates or expand.
    inline void buildFromList(int front) {
      // printf("hash: %s\n", __func__);
      // note that there could still be stale entries after this
      // note insertWithoutExpand() will increment depth[i].

#ifdef MAX_DEPTH
      int* temp = (int*) realloc(depth, sizeof(int) * getSize());
      assert(temp != NULL);
      depth = temp;
      memset(depth, 0, sizeof(int) * getSize());
#endif

      int next;
      for ( ; front != nodes->getNull(); front = next) {
        next = nodes->getNext(front);
#if DISCARD_STALES
        if (nodes->isStale(front)) {
          nodes->uncacheNode(front);
        } else {
          insertWithoutExpand(front);
        }
#else
        insertWithoutExpand(front);
#endif
      }
    }

    inline int insertWithoutExpand(int key) {
      unsigned h = nodes->hash(key, getSize());
      CHECK_RANGE(0, h, getSize());

      // insert at head of list

#ifndef MAX_DEPTH
      nodes->setNext(key, table[h]);
      table[h] = key;
      nEntries++;
#else
      nodes->setNext(key, nodes->getNull());
      // Two situations:
      // MAX_DEPTH reached: remove tail, insert at tail
      // MAX_DEPTH not reached: insert at tail
      DCASSERT(depth[h] <= MAX_DEPTH);
      if (depth[h] == int(MAX_DEPTH)) {
        DCASSERT(int(MAX_DEPTH) > 0);
        DCASSERT(depth[h] > 0);
        // find last entry
        // there is at least 1 entry in this list
        // remove the last entry from cache
        // insert at tail
#if MAX_DEPTH == 1
        nodes->uncacheNode(table[h]);
        table[h] = key;
#elif MAX_DEPTH == 2
        nodes->uncacheNode(nodes->getNext(table[h]));
        nodes->setNext(table[h], key);
#elif MAX_DEPTH == 3
        int prev = nodes->getNext(table[h]);
        nodes->uncacheNode(nodes->getNext(prev));
        nodes->setNext(prev, key);
#elif MAX_DEPTH == 4
        int prev = nodes->getNext(nodes->getNext(table[h]));
        nodes->uncacheNode(nodes->getNext(prev));
        nodes->setNext(prev, key);
#else
        assert(false);
        int prev = nodes->getNext(nodes->getNext(nodes->getNext(table[h])));
        for (int i = MAX_DEPTH - 5; i > 0; i--)
          prev = nodes->getNext(prev);
        nodes->uncacheNode(nodes->getNext(prev));
        nodes->setNext(prev, key);
#endif
        // nEntries and depth[h] remain the same
        discarded++;
      } else {
        // insert at tail
        DCASSERT(depth[h] < int(MAX_DEPTH));
        switch (depth[h]) {
          case 0:
            table[h] = key;
            break;
          case 1:
            nodes->setNext(table[h], key);
            break;
          default:
            int last = nodes->getNext(table[h]);
            for (int i = depth[h] - 2; i > 0; i--)
              last = nodes->getNext(last);
            nodes->setNext(last, key);
        }
        nEntries++;
        depth[h]++;
      }
#endif

      return key;
    }

    /// Go through the hash table and remove all entries
    inline void clear() {
      int temp;
      for (unsigned i = 0; i < getSize(); i++) {
        while (table[i] != nodes->getNull()) {
          temp = table[i];
          table[i] = nodes->getNext(table[i]);
          nodes->uncacheNode(temp);
        }
      }

      nEntries = 0;
      setSize(getMinSize());

      table = (int*) realloc(table, sizeof(int) * getSize());
      assert(NULL != table);
      int *end = table + getSize();
      for (int* curr = table; curr < end; curr++) {
        *curr = nodes->getNull();
      }

#ifdef MAX_DEPTH
      depth = (int*) realloc(depth, sizeof(int) * getSize());
      assert(NULL != depth);
      memset(depth, 0, sizeof(int) * getSize());
#endif

      noMem = false;
    }


    /// Go through the hash table and remove all stale entries
    /// If shrinkable is true, it will reduce the size of the
    /// hash table when there is a certain amount of wasted space in it.
    inline unsigned removeStaleEntries(bool shrinkable = true) {
      unsigned total_stale_count = 0;
      unsigned stale_count = 0;
      int prev, curr, next;
      do {
        stale_count = 0;
        for (unsigned i = 0; i < getSize(); i++) {
#if 1
          if (table[i] == nodes->getNull()) continue;
          while (table[i] != nodes->getNull() && nodes->isStale(table[i])) {
            next = nodes->getNext(table[i]);
            nodes->uncacheNode(table[i]);
#ifdef MAX_DEPTH
            depth[i]--;
#endif
            table[i] = next;
            stale_count++;
          }
          if (table[i] == nodes->getNull()) continue;

          for (prev = table[i], curr = nodes->getNext(prev);
              curr != nodes->getNull();
              curr = next) {
            next = nodes->getNext(curr);
            if (nodes->isStale(curr)) {
              nodes->setNext(prev, next);
              nodes->uncacheNode(curr);
#ifdef MAX_DEPTH
              depth[i]--;
#endif
              stale_count++;
              // prev remains the same
            } else {
              // advance prev
              prev = curr;
            }
          }
#else
          int front = nodes->getNull();  // empty "list"
          while (table[i] != nodes->getNull()) {
            next = nodes->getNext(table[i]);
            if (nodes->isStale(table[i])) {
              // don't add to list
              // remove from cache
              nodes->uncacheNode(table[i]);
#ifdef MAX_DEPTH
              depth[i]--;
#endif
              stale_count++;
            } else {
              nodes->setNext(table[i], front);
              front = table[i];
            }
            table[i] = next;
          }
          table[i] = front;  // list is reversed but that's ok
#endif
        } // for i
        total_stale_count += stale_count;
#if AGGRESIVE_GC
      } while (stale_count > 0);
#else
      } while (false);
#endif

      nEntries -= total_stale_count;
      if (shrinkable && getSize() > getMinSize() &&
          nEntries < (getSize() >> 1))
        shrink();
      return total_stale_count;
    }


    // Is is time to expand the hash table?
    inline bool isTimeToExpand() {
      return (getSize() < getMaxSize() &&
          nEntries > getSizeBy4() && !noMem);
    }

    // Expands the hash table (the number of buckets) if necessary
    // (based on isTimeToExpand()). After expanding it rehashes the entries.
    inline void expand() {
#ifdef DEBUG_HASH_EXPAND_H
      fprintf(stdout, "%s: Trying to enlarge table...\n", __func__);
      // fprintf(stdout, "%s: Old table:\n", __func__);
      // show(stdout);
#endif
      // removeStaleEntries(false);
      if (!isTimeToExpand()) return;
      unsigned length = 0;
      int ptr = convertToList(length);
#ifdef EXPANSION_FACTOR
      unsigned newSize = int(size * float(EXPANSION_FACTOR));
#else
      unsigned newSize = size << 1;
#endif
      if (newSize > maxTableSize) newSize = maxTableSize;
#ifdef DEBUG_HASH_EXPAND_H
      printf("%s: Enlarging table...\n", __func__);
      printf(" nEntries = %d", length);
      printf(", new size = %d\n", newSize);
#endif
      int *temp = (int*) realloc(table, sizeof(int) * newSize);
      if (NULL == temp) {
        fprintf(stderr, "Memory allocation error in unique table expand.\n");
        noMem = true;
        buildFromList(ptr);
        return;
      }
      
      table = temp;
      setSize(newSize);
      int *last = table + size;
      for (int *curr = table; curr < last; ++curr) {
        *curr = nodes->getNull();
      }

      buildFromList(ptr);
#ifdef DEBUG_HASH_EXPAND_H
      // fprintf(stdout, "%s: New table:\n", __func__);
      // show(stdout);
#endif
    }

    inline void shrink() {
#ifdef DEBUG_HASH_H
      fprintf(stdout, "Shrinking table.  Old table:\n");
      show(stdout);
#endif
      if (getSize() <= getMinSize()) return;
#if 0
      fprintf(stdout, "Shrinking table from size %d", size);
      fflush(stdout);
#endif
      unsigned length = 0;
      int front = convertToList(length);
      unsigned newSize = size >> 1;
      while (newSize > 0u && length < newSize) {
        newSize = newSize >> 1;
      }
      newSize = newSize << 1;  // otherwise list won't fit
      if (newSize < getMinSize()) newSize = getMinSize();
#if 0
      fprintf(stdout, "to size %d.", newSize);
      fflush(stdout);
#endif
      int *temp = (int*) realloc(table, sizeof(int) * newSize);
      if (NULL == temp) {
        fprintf(stderr, "Memory allocation error in hash table resize.\n");
        exit(1);
      } else {
        table = temp;
        setSize(newSize);
      }
      int *last = table + size;
      for (int *curr = table; curr < last; ++curr) {
        *curr = nodes->getNull();
      }
      noMem = false;
      buildFromList(front);
#if 0
      fprintf(stdout, "Done resizing\n");
      fflush(stdout);
#endif
#ifdef DEBUG_HASH_H
      fprintf(stdout, "New table:\n");
      show(stdout);
#endif
    }

    /** If table contains key, move it to the front of the list.
        Otherwise, do nothing.
        Returns index of the item if found, getNull() otherwise.
    */
    int find(int key) {
      unsigned h = nodes->hash(key, getSize());
      CHECK_RANGE(0, h, getSize());
      if (table[h] == nodes->getNull()) { return nodes->getNull(); }
      int ptr = nodes->getNull();
#if STALE_BEFORE_EQUAL
      while (table[h] != nodes->getNull()) {
        if (nodes->isStale(table[h])) {
          ptr = table[h];
          table[h] = nodes->getNext(table[h]);
          nodes->uncacheNode(ptr);
          nEntries--;
#ifdef MAX_DEPTH
          depth[h]--;
#endif
        } else if (nodes->equals(table[h], key)) {
          return table[h];
        } else {
          break;
        }
      }
      if (table[h] == nodes->getNull()) { return nodes->getNull(); }
      int prev = table[h];
      ptr = nodes->getNext(table[h]);
      int next = nodes->getNull();
      for ( ; ptr != nodes->getNull(); ptr = next) {
        next = nodes->getNext(ptr);
        if (nodes->isStale(ptr)) {
          nodes->setNext(prev, next);
          nodes->uncacheNode(ptr);
          // prev stays the same
          nEntries--;
#ifdef MAX_DEPTH
          depth[h]--;
#endif
        } else if (nodes->equals(ptr, key)) {
          // move it to front and return
          nodes->setNext(prev, next);
          nodes->setNext(ptr, table[h]);
          table[h] = ptr;
          return ptr;
        } else {
          prev = ptr;
        }
      }
      return nodes->getNull();
#else
      while (table[h] != nodes->getNull()) {
        if (nodes->equals(table[h], key)) {
          return table[h];
#if DISCARD_STALES
        } else if (nodes->isStale(table[h])) {
          ptr = table[h];
          table[h] = nodes->getNext(table[h]);
          nodes->uncacheNode(ptr);
          nEntries--;
#ifdef MAX_DEPTH
          depth[h]--;
#endif
#endif
        } else {
          break;
        }
      }
      if (table[h] == nodes->getNull()) { return nodes->getNull(); }
      int prev = table[h];
      ptr = nodes->getNext(table[h]);
      int next = nodes->getNull();
      for ( ; ptr != nodes->getNull(); ptr = next) {
        next = nodes->getNext(ptr);
        if (nodes->equals(ptr, key)) {
          // move it to front and return
          nodes->setNext(prev, next);
          nodes->setNext(ptr, table[h]);
          table[h] = ptr;
          return ptr;
#if DISCARD_STALES
        } else if (nodes->isStale(ptr)) {
          nodes->setNext(prev, next);
          nodes->uncacheNode(ptr);
          // prev stays the same
          nEntries--;
#ifdef MAX_DEPTH
          depth[h]--;
#endif
#endif
        } else {
          prev = ptr;
        }
      }
      return nodes->getNull();
#endif
    }


    /** If table contains key, remove it and return it; otherwise return
        getNull(). Note that the key in this case represents a logical
        node address -- the contents of the node are not investigated
        to determine equality.
    */
    inline int remove(int key) {
      unsigned h = nodes->hash(key, getSize());
      CHECK_RANGE(0, h, getSize());
      int& head = table[h];

      int result = nodes->getNull();

      if (head != nodes->getNull()) {
        if (head == key) {
          int next = nodes->getNext(head);
          nodes->uncacheNode(head);
          head = next;
          nEntries--;
          result = key;
        }
        else {
          for (int prev = head, curr = nodes->getNext(head),
              end = nodes->getNull(); curr != end; )
          {
            if (curr == key) {
              nodes->setNext(prev, nodes->getNext(curr));
              nodes->uncacheNode(curr);
              nEntries--;
              result = key;
              break;
            }
            prev = curr;
            curr = nodes->getNext(curr);
          }
        }
      }

#ifdef MAX_DEPTH
      if (result != nodes->getNull()) depth[h]--;
#endif

      return result;
    }


    /** Add key to the front of the list.
        Returns the (new) front of the list.

        Note: It is up to the user to search for duplicates
        (using find(..)) before inserting.
    */
    inline int insert(int key) {
      // DCASSERT(find(key) == nodes->getNull());
      if (isTimeToExpand()) expand();
      return insertWithoutExpand(key);
    }
};

#endif

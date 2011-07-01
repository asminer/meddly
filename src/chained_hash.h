
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

#ifndef CHAINED_HASH_H
#define CHAINED_HASH_H

//#define CHECK_STALE_ONLY_ON_EQUAL
//#define DEBUG_CT_ENTRIES
// #define DEBUG_HASH_EXPAND_H

#include <sys/time.h>
#include <sys/resource.h>
#include <vector>
#include "defines.h"

const unsigned  minTableSize = 1024u;
const unsigned maxLoadFactor = 4u;  // i.e. nEntries / tableSize < 10

template <typename MANAGER>
class chained_hash_table {

  protected:
    unsigned    tableSize;
    unsigned    nEntries;
    unsigned    maxEntries;
    unsigned    maxTableSize;
    std::vector<int> table;
    MANAGER*    nodes;

  public:
    chained_hash_table(MANAGER *n, unsigned max = 16777216u)
    : nEntries(0u), maxEntries(max/2), maxTableSize(max/2), nodes(n)
    {
      // nEntries can increase to 2 x maxEntries, therefore, maxEntries is
      // reduced by a factor of 2.

      table.resize(minTableSize, nodes->getNull());

      std::vector<int>::iterator iter = table.begin();
      for ( ; iter != table.end(); )
      {
        *iter++ = nodes->getNull();
      }
    }

    ~chained_hash_table() { }

    inline unsigned getTableSize() const { return table.size(); }
    inline unsigned getMemoryUsage() const {
      return table.size() * sizeof(int);
    }
    inline unsigned getEntriesCount() const { return nEntries; }

    void show(FILE *s, bool verbose = false) const {
      char filler[] = "\t";
      fprintf(s, "%sNumber of slots:        %lu\n",
          filler, (unsigned long)(table.size()));
      fprintf(s, "%sNumber of entries:      %d\n", filler, nEntries);
#ifdef CHAIN_LENGTH_INFO
      fprintf(s, "%sChains:\n", filler); 
      fprintf(s, "%s   Depth    Count\n", filler); 
      // compute chain distribution
      std::vector<unsigned> chains;
      for (unsigned i = 0; i < table.size(); i++) {
        unsigned depth = 0;
        for (int index = table[i];
            index != nodes->getNull();
            index = nodes->getNext(index), depth++);
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
        for (unsigned i = 0u; i < table.size(); i++) {
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
      for (unsigned i = 0; i < table.size(); i++) {
        int next = nodes->getNull();
        for (int curr = table[i]; curr != nodes->getNull(); curr = next) {
          next = nodes->getNext(curr);
          if (nodes->isStale(curr)) {
            // remove from cache
            nodes->uncacheNode(curr);
#ifdef DEBUG_CT_ENTRIES
            printf("Removed stale CT entry %d\n", curr);
#endif
          } else {
            // add to head of list
            nodes->setNext(curr, front);
            front = curr;
            listlength++;
          }
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

      for ( ; front != nodes->getNull(); ) {
        int next = nodes->getNext(front);
        if (nodes->isStale(front)) {
          nodes->uncacheNode(front);
#ifdef DEBUG_CT_ENTRIES
          printf("Removed stale CT entry %d\n", front);
#endif
        } else {
          insertWithoutExpand(front);
        }
        front = next;
      }
    }

    /// Removes an entry from the list at table[h].
    /// Either a stale entry or the last entry in the list is removed.
    /// If table[h] is empty, the nothing is removed.
    inline void removeAnEntry(int h) {
      // remove one element:
      // either a stale entry or an entry from the end of the list

      if (table[h] == nodes->getNull()) return;

      if (nodes->isStale(table[h])) {
        int stale = table[h];
        table[h] = nodes->getNext(table[h]);
        nodes->uncacheNode(stale);
        nEntries--;
#ifdef DEBUG_CT_ENTRIES
        printf("Removed stale CT entry %d (#entries %d)\n", stale, nEntries);
#endif
        return;
      }

      int prev = table[h];
      int curr = nodes->getNext(table[h]);

      if (curr == nodes->getNull()) {
        // remove table[h]
        nodes->uncacheNode(table[h]);
        nEntries--;
#ifdef DEBUG_CT_ENTRIES
        printf("Removed stale CT entry %d (#entries %d)\n", table[h], nEntries);
#endif
        table[h] = nodes->getNull();
        return;
      }

      if (nodes->isStale(curr)) {
        nodes->setNext(table[h], nodes->getNext(curr));
        nodes->uncacheNode(curr);
        nEntries--;
#ifdef DEBUG_CT_ENTRIES
        printf("Removed stale CT entry %d (#entries %d)\n", curr, nEntries);
#endif
        return;
      }

      int next = nodes->getNext(curr);

      while (next != nodes->getNull()) {
        if (nodes->isStale(next)) {
          nodes->setNext(curr, nodes->getNext(next));
          nodes->uncacheNode(next);
          nEntries--;
#ifdef DEBUG_CT_ENTRIES
          printf("Removed stale CT entry %d (#entries %d)\n", next, nEntries);
#endif
          return;
        }
        prev = curr;
        curr = next;
        next = nodes->getNext(next);
      }

      // If you get here then there are atleast two entries remaining
      // in this list (prev and curr). Remove curr.
      nodes->setNext(prev, next);
      nodes->uncacheNode(curr);
      nEntries--;
#ifdef DEBUG_CT_ENTRIES
      printf("Removed stale CT entry %d (#entries %d)\n", curr, nEntries);
#endif
    }

    /// If table contains key, remove it and return it; otherwise return
    /// getNull(). Note that the key in this case represents a logical
    /// node address -- the contents of the node are not investigated
    /// to determine equality.
    inline int remove(int key) {
      CHECK_RANGE(0, nodes->hash(key, table.size()), table.size());
      int& head = table[nodes->hash(key, table.size())];
      int result = nodes->getNull();

      if (head != nodes->getNull()) {
        if (head == key) {
          int next = nodes->getNext(head);
          nodes->uncacheNode(head);
          nEntries--;
#ifdef DEBUG_CT_ENTRIES
          printf("Removed stale CT entry %d (#entries %d)\n", head, nEntries);
#endif
          head = next;
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
#ifdef DEBUG_CT_ENTRIES
              printf("Removed stale CT entry %d (#entries %d)\n", curr, nEntries);
#endif
              result = key;
              break;
            }
            prev = curr;
            curr = nodes->getNext(curr);
          }
        }
      }

      return result;
    }


    /// Insert new entry without expanding the hash table
    inline int insertWithoutExpand(int key) {
      // find hash slot
      unsigned h = nodes->hash(key, table.size());
      CHECK_RANGE(0, h, table.size());

      if (nEntries >= maxEntries) removeAnEntry(h);

      // insert at head of list
      nodes->setNext(key, table[h]);
      table[h] = key;
      nEntries++;

#ifdef DEBUG_CT_ENTRIES
      printf("Added CT entry %d (#entries = %d)\n", key, nEntries);
#endif

      return key;
    }

    /// Go through the hash table and remove all entries
    inline void clear() {
      for (std::vector<int>::iterator iter = table.begin();
          iter != table.end(); ++iter) {
        while (*iter != nodes->getNull()) {
          int curr = *iter;
          *iter = nodes->getNext(*iter);
          nodes->uncacheNode(curr);
        }
      }

      nEntries = 0;
      table.resize(minTableSize, nodes->getNull());
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
        for (unsigned i = 0; i < table.size(); i++) {
          if (table[i] == nodes->getNull()) continue;
          while (table[i] != nodes->getNull() && nodes->isStale(table[i])) {
            next = nodes->getNext(table[i]);
            nodes->uncacheNode(table[i]);
#ifdef DEBUG_CT_ENTRIES
            printf("Removed stale CT entry %d\n", table[i]);
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
#ifdef DEBUG_CT_ENTRIES
              printf("Removed stale CT entry %d\n", curr);
#endif
              stale_count++;
              // prev remains the same
            } else {
              // advance prev
              prev = curr;
            }
          }
        } // for i
        total_stale_count += stale_count;
#if AGGRESIVE_GC
      } while (stale_count > 0);
#else
      } while (false);
#endif

      nEntries -= total_stale_count;
      if (shrinkable && table.size() > minTableSize &&
        nEntries / table.size() < 0.5)
        shrink();
      return total_stale_count;
    }


    // Is is time to expand the hash table?
    inline bool isTimeToExpand() {
      return (nEntries > maxLoadFactor * table.size()
          && table.size() < maxTableSize);
    }

    // Expands the hash table (the number of buckets) if necessary
    // (based on isTimeToExpand()). After expanding it rehashes the entries.
    inline void expand() {

#ifdef DEBUG_SLOW
      fprintf(stdout, "Enlarging compute table (current size %ld)\n", table.size());
#endif
#ifdef DEBUG_HASH_EXPAND_H
      fprintf(stdout, "%s: Trying to enlarge table...\n", __func__);
#endif

      unsigned length = 0;
      int ptr = convertToList(length);
      unsigned newSize = table.size() * 2;
      if (newSize > maxTableSize) newSize = maxTableSize;

#ifdef DEBUG_HASH_EXPAND_H
      printf("%s: Enlarging table...\n", __func__);
      printf(" nEntries = %d", length);
      printf(", new size = %d\n", newSize);
#endif

      table.resize(newSize, nodes->getNull());
      buildFromList(ptr);
    }

    inline void shrink() {
      if (table.size() <= minTableSize) return;
      if (nEntries / table.size() > 0.5) return;
#ifdef DEBUG_SLOW
      fprintf(stderr, "Shrinking compute table (current size %ld)\n", table.size());
#endif
#ifdef DEBUG_HASH_H
      fprintf(stdout, "Shrinking table.  Old table:\n");
      show(stdout);
#endif

      unsigned length = 0u;
      int front = convertToList(length);
      unsigned newSize = table.size() / 2;
      while (newSize > minTableSize && length/newSize < 0.5) { newSize /= 2; }
      newSize *= 2;
      if (newSize < minTableSize) newSize = minTableSize;

      table.resize(newSize, nodes->getNull());
      buildFromList(front);
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
      unsigned h = nodes->hash(key, table.size());
      CHECK_RANGE(0, h, table.size());

      if (table[h] == nodes->getNull()) { return nodes->getNull(); }

#ifdef CHECK_STALE_ONLY_ON_EQUAL
      // Check the head of the list for match.
      // If a match is found, check if it is a stale entry.
      // If the entry is stale, discard it and return null,
      // otherwise, return the matching entry.
      if (nodes->equals(table[h], key)) {
        if (!nodes->isStale(table[h])) return table[h];
        // Remove stale entry and return null.
        int ptr = table[h];
        table[h] = nodes->getNext(table[h]);
        nodes->uncacheNode(ptr);
        nEntries--;
#ifdef DEBUG_CT_ENTRIES
        printf("Removed stale CT entry %d (#entries %d)\n", ptr, nEntries);
#endif
        return nodes->getNull();
      }

      // Check the rest of the list for match.
      int prev = table[h];
      int ptr = nodes->getNext(table[h]);
      for ( ; ptr != nodes->getNull(); ) {
        if (nodes->equals(ptr, key)) {
          if (!nodes->isStale(ptr)) {
            // Move it to front and return.
            nodes->setNext(prev, nodes->getNext(ptr));
            nodes->setNext(ptr, table[h]);
            table[h] = ptr;
            return ptr;
          }
          // Remove stale entry and return null.
          nodes->setNext(prev, nodes->getNext(ptr));
          nodes->uncacheNode(ptr);
          nEntries--;
#ifdef DEBUG_CT_ENTRIES
          printf("Removed stale CT entry %d (#entries %d)\n", ptr, nEntries);
#endif
          return nodes->getNull();
        }
        // Advance pointers
        prev = ptr;
        ptr = nodes->getNext(ptr);
      }
#else
      // Check the head of the list for match.
      // Remove stale entries before checking for match.
      int ptr = nodes->getNull();
      while (table[h] != nodes->getNull()) {
        if (nodes->isStale(table[h])) {
          ptr = table[h];
          table[h] = nodes->getNext(table[h]);
          nodes->uncacheNode(ptr);
          nEntries--;
#ifdef DEBUG_CT_ENTRIES
          printf("Removed stale CT entry %d (#entries %d)\n", ptr, nEntries);
#endif
        } else if (nodes->equals(table[h], key)) {
          return table[h];
        } else {
          break;
        }
      }
      if (table[h] == nodes->getNull()) { return nodes->getNull(); }

      // Check the rest of the list for match.
      // Remove stale entries before checking for match.
      int prev = table[h];
      ptr = nodes->getNext(table[h]);
      int next = nodes->getNull();
      for ( ; ptr != nodes->getNull(); ptr = next) {
        next = nodes->getNext(ptr);
        if (nodes->isStale(ptr)) {
          nodes->setNext(prev, next);
          nodes->uncacheNode(ptr);
          nEntries--;
#ifdef DEBUG_CT_ENTRIES
          printf("Removed stale CT entry %d (#entries %d)\n", ptr, nEntries);
#endif
        } else if (nodes->equals(ptr, key)) {
          // Move it to front and return.
          nodes->setNext(prev, next);
          nodes->setNext(ptr, table[h]);
          table[h] = ptr;
          return ptr;
        } else {
          prev = ptr;
        }
      }
#endif

      // Not found, return null
      return nodes->getNull();
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

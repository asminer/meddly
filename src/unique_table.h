
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


#ifndef UNIQUE_TABLE_H
#define UNIQUE_TABLE_H

#include "defines.h"

namespace MEDDLY {
  class unique_table;
};

/** Unique table for discovering duplicate nodes.

    This is now a stand-alone, non-template class
    designed specifically for expert_forests.
*/
class MEDDLY::unique_table {
private:
    class subtable {
    public:
    	subtable();
    	~subtable();

        inline unsigned getSize() const         { return size; }
        inline unsigned getNumEntries() const   { return num_entries; }
        inline unsigned getMemUsed() const      { return size * sizeof(node_handle); }

        void reportStats(FILE* s, const char* pad, unsigned flags) const;

        /// For debugging
        void show(FILE *s) const;

        /**
         * Initialize the sub table. Must be called before use.
         */
    	void init(expert_forest *ef);

        /** If table contains key, move it to the front of the list.
            Otherwise, do nothing.
            Returns the item if found, 0 otherwise.

            Class T must have the following methods:
              unsigned hash():    return the hash value for this item.
              bool equals(int p): return true iff this item equals node p.
         */
        template <typename T>
        int find(const T &key);

        /** Add the item to the front of the list.
            Used when we KNOW that the item is not in the unique table already.
        */
        void add(unsigned hash, node_handle item);

        /** If table contains key, remove it and return it.
          I.e., the exact key.
          Otherwise, return 0.
        */
        int remove(unsigned hash, node_handle item);

        /**
         * Retrieve the items in the table.
         * The argument sz specifies the number of items to be retrieved.
         * If sz is larger than the number of items in the table, only the number
         * of items in the table are to be retrieved.
         * Returns the actual number of items retrieved.
         */
        int getItems(node_handle* items, int sz) const;

    private:  // helper methods
      /// Empty the hash table into a list; returns the list.
      node_handle convertToList();

      /// A series of inserts; doesn't check for duplicates or expand.
      void buildFromList(node_handle front);

      /// Expand the hash table (if possible)
      void expand();

      /// Shrink the hash table
      void shrink();

    private:
        expert_forest* parent;
        unsigned size;
        unsigned num_entries;
        unsigned next_expand;
        unsigned next_shrink;
        node_handle* table;

        static const unsigned MAX_SIZE = 1073741824;
        static const unsigned MIN_SIZE = 8;
    };

public:
    unique_table(expert_forest *ef);
    ~unique_table();

    unsigned getSize() const;
    unsigned getNumEntries() const;
    unsigned getMemUsed() const;

    void reportStats(FILE* s, const char* pad, unsigned flags) const;

    /// For debugging
    void show(FILE *s) const;

    /** If the table at the given level contains key, move it to the front of the list.
        Otherwise, do nothing.
        Returns the item if found, 0 otherwise.

        Class T must have the following methods:
          unsigned hash():    return the hash value for this item.
          bool equals(int p): return true iff this item equals node p.
     */
    template <typename T>
    node_handle find(const T &key, int level);

    /** Add the item to the front of the list at the corresponding level.
        Used when we KNOW that the item is not in the unique table already.
    */
    void add(unsigned h, node_handle item);

    /** If the table at the corresponding level contains key, remove it and return it.
      I.e., the exact key.
      Otherwise, return 0.
    */
    node_handle remove(unsigned h, node_handle item);

    /**
     * Retrieve the items in the table at the given level.
     * The argument sz specifies the number of items to be retrieved.
     * If sz is larger than the number of items in the table, only the number
     * of items in the table are to be retrieved.
     * Returns the actual number of items retrieved.
     */
    int getItems(int level, node_handle* items, int sz) const;

private:
    expert_forest* parent;
    subtable* tables;
};

template <typename T>
MEDDLY::node_handle MEDDLY::unique_table::subtable::find(const T &key)
{
	unsigned h = key.hash() % size;
	MEDDLY_CHECK_RANGE(0, h, size);
	node_handle prev = 0;
	for (node_handle ptr = table[h]; ptr!=0; ptr = parent->getNext(ptr)) {
		if (parent->areDuplicates(ptr, key)) { // key.equals(ptr)) {
			// MATCH
			if (ptr != table[h]) {
				// Move to front
				MEDDLY_DCASSERT(prev);
				parent->setNext(prev, parent->getNext(ptr));
				parent->setNext(ptr, table[h]);
				table[h] = ptr;
			}
			MEDDLY_DCASSERT(table[h] == ptr);
			return ptr;
		}
		prev = ptr;
	}

	return 0;
}

template <typename T>
inline MEDDLY::node_handle MEDDLY::unique_table::find(const T &key, int level)
{
	return tables[level].find(key);
}

inline void MEDDLY::unique_table::add(unsigned hash, node_handle item)
{
	int level=parent->getNodeLevel(item);
	tables[level].add(hash, item);
}

inline MEDDLY::node_handle MEDDLY::unique_table::remove(unsigned hash, node_handle item)
{
	int level=parent->getNodeLevel(item);
	return tables[level].remove(hash, item);
}

inline int MEDDLY::unique_table::getItems(int level, node_handle* items, int sz) const
{
	return tables[level].getItems(items, sz);
}

#endif

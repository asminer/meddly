
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



/*! \file compute_table.h

    Compute table interface.

    This interface is for "expert" interface users who wish to implement
    user-defined operations.

    A compute table is used to cache the results of operations on MDDs.
    An expert-user wishing to cache results of user-defined operations may
    need to use this interface.

*/

#ifndef COMPUTE_TABLE_H
#define COMPUTE_TABLE_H

#include "defines.h"

namespace MEDDLY {
  class compute_table;
  struct settings;

  /** Build a new, monolithic table.
      Monolithic means that the table stores entries for several
      (ideally, all) operations.
  */
  compute_table* createMonolithicTable(const settings &s); 

  /** Build a new table for a single operation.
  */
  compute_table* createOperationTable(const settings &s, operation* op);
}


class MEDDLY::compute_table {
    public:
      /// If true, use a table with chaining; otherwise, don't chain.
      bool chaining;
      /// The maximum size of the hash table.
      unsigned maxSize;
      /// Aggressively try to eliminate atale entries.
      bool checkStales;

      struct stats {
        long numEntries;
        unsigned hits;
        unsigned pings;
        static const int chainHistogramSize = 256;
        long chainHistogram[chainHistogramSize];
        long numLargeChains;
        int maxChainLength;
      };

      class search_key {
          friend class base_table;
          friend class monolithic_table;
          friend class operation_table;
          int hashLength;
          int* data;
          int* key_data;
#ifdef DEVELOPMENT_CODE
          int keyLength;  // used for range checking
#endif
        public:
          search_key();
          ~search_key();
          inline int& key(int i) { 
#ifdef DEVELOPMENT_CODE
            assert(i>=0);
            assert(i<keyLength);
#endif
            return key_data[i]; 
          }
      };

      class temp_entry {
          friend class base_table;
          friend class monolithic_table;
          friend class operation_table;
          int handle;
          int hashLength;
          int* entry;
          int* key_entry;
          int* res_entry;
#ifdef DEVELOPMENT_CODE
          int keyLength;
          int resLength;
#endif
        public:
          inline int& key(int i) { 
#ifdef DEVELOPMENT_CODE
            assert(i>=0);
            assert(i<keyLength);
#endif
            return key_entry[i]; 
          }
          inline int& result(int i) { 
#ifdef DEVELOPMENT_CODE
            assert(i>=0);
            assert(i<resLength);
#endif
            return res_entry[i]; 
          }
          // inline void copyResult(int i, void* data, size_t bytes) {
          void copyResult(int i, void* data, size_t bytes) {
#ifdef DEVELOPMENT_CODE
            assert(i>=0);
            assert(i+bytes<=resLength*sizeof(int));
#endif
            memcpy(res_entry+i, data, bytes);
          }
      };

    public:
      /// Constructor
      compute_table(const settings &s);

      /** Destructor. 
          Does NOT properly discard all table entries;
          use \a clear() for this.
      */
      virtual ~compute_table();

      /// Is this a per-operation compute table?
      virtual bool isOperationTable() const = 0;

      /// Initialize a search key for a given operation.
      virtual void initializeSearchKey(search_key &key, operation* op) = 0;

      /** Find an entry in the compute table based on the key provided.
          @param  key   Key to search for.
          @return       0, if not found;
                        otherwise, an integer array of size 
                        op->getCacheEntryLength()
      */
      virtual const int* find(const search_key &key) = 0;

      /** Start a new compute table entry.
          The operation should "fill in" the values for the entry,
          then call \a addEntry().
      */
      virtual temp_entry& startNewEntry(operation* op) = 0;

      /** Add the "current" new entry to the compute table.
          The entry may be specified by filling in the values 
          for the struct returned by \a startNewEntry().
      */
      virtual void addEntry() = 0;

      /** Remove all stale entries.
          Scans the table for entries that are no longer valid (i.e. they are
          stale, according to operation::isEntryStale) and removes them. This
          can be a time-consuming process (proportional to the number of cached
          entries).
      */
      virtual void removeStales() = 0;

      /** Removes all entries.
      */
      virtual void removeAll() = 0;

      /// Get performance stats for the table.
      inline const stats& getStats() {
        return perf;
      }

      /// For debugging.
      virtual void show(FILE *s, int verbLevel = 0) const = 0;

    protected:
      stats perf;
      temp_entry currEntry;
};

#endif


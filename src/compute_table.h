
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

  class compute_table {
    public:
      struct settings {
        /// If true, use a table with chaining; otherwise, don't chain.
        bool chaining;
        /// The maximum size of the hash table.
        unsigned maxSize;

        settings(); // to set default values, from MEDDLY::settings
      };

      struct stats {
        long numEntries;
        unsigned hits;
        unsigned pings;
      };

    public:
      /// Constructor
      compute_table(settings s);

      /** Destructor. 
          Does NOT properly discard all table entries;
          use \a clear() for this.
      */
      virtual ~compute_table();

      /// Is this a per-operation compute table?
      virtual bool isOperationTable() const = 0;

      /** Add an entry to the compute table. Note that this table allows for
          duplicate entries for the same key. The user may use find() before
          using add() to prevent such duplication. A copy of the data
          in entry[] is stored in the table.
          @param  op    Operation associated with this table entry
          @param  entry integer array of size op->getCacheEntryLength(),
                        containing the operands and the result to be stored
      */    
      virtual void add(operation* op, const int* entry) = 0;

      /** Find an entry in the compute table based on the key provided.
          If more than an entry with the same key exists, this will return
          the first matching entry.
          @param  op    Operation associated with this table entry
          @param  entry integer array of size op->getKeyLength(),
                        containing the key to the table entry to look for
          @return       integer array of size op->getCacheEntryLength()
      */
      virtual const int* find(operation* op, const int* entryKey) = 0;

      /** Remove stale entries.
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
        updateStats();
        return perf;
      }

    protected:
      settings opts;
      stats perf;
      static unsigned raw_hash(const int* data, int length);

      virtual void updateStats() = 0;
  };


  /** Build a new, monolithic table.
      Monolithic means that the table stores entries for several
      (ideally, all) operations.
  */
  compute_table* createMonolithicTable(compute_table::settings s); 

  /** Build a new table for a single operation.
  */
  compute_table* createOperationTable(compute_table::settings, operation* op);

} // namespace MEDDLY

#endif

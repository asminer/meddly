
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


#ifndef CT_NONE_H
#define CT_NONE_H

#include "../hash_stream.h"

// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                           ct_none  class                           *
// *                                                                    *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {
  template <bool MONOLITHIC, bool CHAINED>
  class ct_none : public compute_table {
    public:
      ct_none(const ct_initializer::settings &s, operation* op, unsigned slot);
      virtual ~ct_none();

      /**
          Find an entry.
          Used by find() and updateEntry().
            @param  key   Key to search for.
            @return Pointer to the result portion of the entry, or null if not found.
      */
      entry_item* findEntry(entry_key* key);

      // required functions

      virtual bool isOperationTable() const   { return !MONOLITHIC; }
      virtual void find(entry_key* key, entry_result &res);
      virtual void addEntry(entry_key* key, const entry_result& res);
      virtual void updateEntry(entry_key* key, const entry_result& res);
      virtual void removeStales();
      virtual void removeAll();
      virtual void show(output &s, int verbLevel = 0);

    private:  // helper methods

      unsigned long convertToList(bool removeStales);
      void listToTable(unsigned long h);

      void scanForStales();
      void rehashTable(const unsigned long* oldT, unsigned long oldS);

      /// Grab space for a new entry
      node_address newEntry(unsigned size);

      /*
          Use our own, built-in, specialized hash
          instead of hash_stream because it's faster.
      */

      // TBD - migrate to 64-bit hash

      static inline unsigned raw_hash(const void* x, unsigned bytes) {
        hash_stream H;
        H.start(12345);
        H.push(x, bytes);
        return H.finish();
      }

      inline unsigned hash(const void* x, unsigned bytes) const {
        return raw_hash(x, bytes) % tableSize;
      }
    
      inline void incMod(unsigned long &h) {
        h++;
        if (h>=tableSize) h=0;
      }

      /// Update stats: we just searched through c items
      inline void sawSearch(unsigned c) {
        if (c>=stats::searchHistogramSize) {
          perf.numLargeSearches++;
        } else {
          perf.searchHistogram[c]++;
        }
        if (c>perf.maxSearchLength) perf.maxSearchLength = c;
      }

      /**
          Try to set table[h] to curr.
          If the slot is full, check ahead the next couple for empty.
          If none are empty, then recycle current table[h]
          and overwrite it.
      */
      inline void setTable(unsigned long h, unsigned long curr) {
        MEDDLY_DCASSERT(!CHAINED);
        unsigned long hfree = h;
        for (int i=maxCollisionSearch; i>=0; i--, incMod(hfree)) {
          // find a free slot
          if (0==table[hfree]) {
            table[hfree] = curr;
            return;
          }
        }
        MEDDLY_DCASSERT(table[h]);
        // full; remove entry at our slot.
#ifdef DEBUG_CT
        printf("Collision; removing CT entry ");
        FILE_output out(stdout);
        showEntry(out, table[h]);
        printf(" handle %lu in slot %lu\n", table[h], h);
        fflush(stdout);
#endif
        collisions++;    

        discardAndRecycle(table[h]);
        table[h] = curr;
      }


      /// Check equality
      static inline bool equal_sw(const entry_item* a, const entry_item* b, int N) {
        /*
        switch (N) {  // note: cases 8 - 2 fall through
          case  8:    if (a[7] != b[7]) return false;
          case  7:    if (a[6] != b[6]) return false;
          case  6:    if (a[5] != b[5]) return false;
          case  5:    if (a[4] != b[4]) return false;
          case  4:    if (a[3] != b[3]) return false;
          case  3:    if (a[2] != b[2]) return false;
          case  2:    if (a[1] != b[1]) return false;
          case  1:    return a[0] == b[0];
          case  0:    return true;
          default:    return (0==memcmp(a, b, N*sizeof(entry_item)));
        };
        */
        return (0==memcmp(a, b, N*sizeof(entry_item)));
      }


      /**
          Check if the key portion of an entry equals key and should be discarded.
            @param  entry   Complete entry in CT to check.
            @param  key     Key to compare against. 
            @param  discard On output: should this entry be discarded

            @return Result portion of the entry, if they key portion matches;
                    0 if the entry does not match the key.
      */
      inline entry_item* checkEqualityAndStatus(entry_item* entry, const entry_key* key, bool &discard)
      {
        entry_item* entry_without_next = CHAINED ? (entry+1) : entry;
        const unsigned keyslots = key->numTempBytes() / sizeof(entry_item);
        // TBD HERE
        if (equal_sw(entry_without_next, (const entry_item*) key->readTempData(), keyslots))
        {
          //
          // Equal.
          //
          entry_item* result = entry_without_next + keyslots;
          discard = isDead(result, key->getET());
          return result;
        } else {
          //
          // Not equal.
          //
          if (checkStalesOnFind) {
            discard = isStale(entry);
          } else {
            discard = false;
          } // if checkStalesOnFind
          return 0;
        }
      }

      /**
          Check if an entry is stale.
            @param  entry   Pointer to complete entry to check.
            @return         true, if the entry should be discarded.
      */
      bool isStale(const entry_item* entry) const;

      /**
          Check if a result is dead (unrecoverable).
            @param  res   Pointer to result portion of entry.
            @param  et    Entry type information.
            @return       true, if the entry should be discarded.
      */
      bool isDead(const entry_item* res, const entry_type* et) const;

      /**
          Copy a result into an entry.
      */
      inline void setResult(entry_item* respart, const entry_result &res, const entry_type* et) {
        const entry_item* resdata = res.rawData();
        memcpy(respart, resdata, res.dataLength()*sizeof(entry_item));
      }

      /// Display a chain 
      inline void showChain(output &s, unsigned long L) const {
        s << L;
        if (CHAINED) {
          while (L) {
#ifdef INTEGRATED_MEMMAN
            const entry_item* entry = entries+L;
#else
            const entry_item* entry = (const int*) MMAN->getChunkAddress(L);
#endif
            L = entry[0].U;
            s << "->" << L;
          } // while L
        }
        s << "\n";
      }

      /**
          Discard an entry and recycle the memory it used.
            @param  h   Handle of the entry.
      */
      void discardAndRecycle(unsigned long h);

      /**
          Display an entry.
          Used for debugging, and by method(s) to display the entire CT.
            @param  s         Stream to write to.
            @param  h         Handle of the entry.
      */
      void showEntry(output &s, unsigned long h) const;

      /// Display a key.
      void showKey(output &s, const entry_key* k) const;


    private:
      /// Global entry type.  Ignored when MONOLITHIC is true.
      const entry_type* global_et;

      /// Hash table
      unsigned long* table;

      /// Hash table size
      unsigned long tableSize;

      /// When to next expand the table
      unsigned long tableExpand;

      /// When to next shrink the table
      unsigned long tableShrink;

#ifdef INTEGRATED_MEMMAN
      /// Memory space for entries
      entry_item*  entries;
      /// Used entries
      unsigned long entriesSize;
      /// Memory allocated for entries
      unsigned long entriesAlloc;

      static const unsigned maxEntrySize = 15;
      static const unsigned maxEntryBytes = sizeof(entry_item) * maxEntrySize;

      /// freeList[i] is list of all unused i-sized entries.
      unsigned long* freeList;
#else
  
      memory_manager* MMAN;
#endif
  
      /// Memory statistics
      memstats mstats;

      /// How many slots to search in unchained tables
      static const unsigned maxCollisionSearch = 2;

      /// Stats: how many collisions
      unsigned long collisions;
  }; // class ct_none
} // namespace


// **********************************************************************
// *                                                                    *
// *                          ct_none  methods                          *
// *                                                                    *
// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
MEDDLY::ct_none<MONOLITHIC, CHAINED>::ct_none(
  const ct_initializer::settings &s, operation* op, unsigned slot)
: compute_table(s)
{
  if (MONOLITHIC) {
    MEDDLY_DCASSERT(0==op);
    MEDDLY_DCASSERT(0==slot);
  } else {
    MEDDLY_DCASSERT(op);
  }
  global_et = op ? getEntryType(op, slot) : 0;

  /*
      Initialize memory management for entries.
  */
#ifdef INTEGRATED_MEMMAN
  freeList = new unsigned long[1+maxEntrySize];
  for (unsigned i=0; i<=maxEntrySize; i++) {
    freeList[i] = 0;
  }
  mstats.incMemUsed( (1+maxEntrySize) * sizeof(unsigned long) );
  mstats.incMemAlloc( (1+maxEntrySize) * sizeof(unsigned long) );

  entriesAlloc = 1024;
  entriesSize = 1;    
  entries = (entry_item*) malloc(entriesAlloc * sizeof(entry_item) );
  if (0==entries) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  entries[0].L = 0;     // NEVER USED; set here for sanity.

  mstats.incMemUsed( entriesSize * sizeof(entry_item) );
  mstats.incMemAlloc( entriesAlloc * sizeof(entry_item) );
#else
  MEDDLY_DCASSERT(s.MMS);
  MMAN = s.MMS->initManager(sizeof(entry_item), 2, mstats);
#endif

  /*
      Initialize hash table
  */
  tableSize = 1024;
  tableExpand = CHAINED ? 4*1024 : 512;
  tableShrink = 0;
  table = (unsigned long*) malloc(tableSize * sizeof(unsigned long));
  if (0==table) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  for (unsigned i=0; i<tableSize; i++) table[i] = 0;

  mstats.incMemUsed(tableSize * sizeof(unsigned long));
  mstats.incMemAlloc(tableSize * sizeof(unsigned long));

  collisions = 0;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
MEDDLY::ct_none<MONOLITHIC, CHAINED>::~ct_none()
{
  /*
    Clean up hash table
  */
  free(table);

  /*
    Clean up memory manager for entries
  */
#ifdef INTEGRATED_MEMMAN
  free(entries);
  delete[] freeList;
#else
  delete MMAN;
#endif

  /*
    Update stats: important for global usage
  */
  mstats.zeroMemUsed();
  mstats.zeroMemAlloc();
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
inline MEDDLY::compute_table::entry_item* MEDDLY::ct_none<MONOLITHIC, CHAINED>
::findEntry(entry_key* key)
{
  MEDDLY_DCASSERT(key);

  // |-----------||--------|-------------|---------|
  // | status    || active | recoverable | dead    |
  // |-----------||--------|-------------|---------|
  // | equal     || use    | use         | discard |
  // |-----------||--------|-------------|---------|
  // | not equal || keep   | discard     | discard |
  // |-----------||--------|-------------|---------|
  //
  // if (equal)
  //    if (dead) discard
  //    else use
  // else
  //    if (!active) discard
  //    else do nothing

  unsigned chain;

  entry_item* answer = 0;
  entry_item* preventry = 0;
  unsigned long hcurr = key->getHash() % tableSize;
  unsigned long curr = table[hcurr];

#ifdef DEBUG_CT_SEARCHES
  printf("Searching for CT entry ");
  FILE_output out(stdout);
  showKey(out, key);
  printf(" in hash slot %lu\n", hcurr);
  fflush(stdout);
#endif

  for (chain=0; ; ) {

    if (0==curr) {
      if (CHAINED) break;
      if (++chain > maxCollisionSearch) break;
      incMod(hcurr);
      curr = table[hcurr];
      continue;
    } else {
      if (CHAINED) chain++;
    }

#ifdef INTEGRATED_MEMMAN
    entry_item* entry = entries + curr;
#else
    entry_item* entry = (entry_item*) MMAN->getChunkAddress(curr);
#endif

    bool discard;
    answer = checkEqualityAndStatus(entry, key, discard);

    if (discard) {
        //
        // Delete the entry.
        //
#ifdef DEBUG_CT
        if (answer) printf("Removing stale CT hit   ");
        else        printf("Removing stale CT entry ");
        FILE_output out(stdout);
        showEntry(out, curr);
        printf(" handle %lu in slot %lu\n", curr, hcurr);
        fflush(stdout);
#endif
        if (CHAINED) {
          if (preventry) {
            preventry[0].UL = entry[0].UL;
          } else {
            table[hcurr] = entry[0].UL;
          }
        } else {
          table[hcurr] = 0;
        }

        discardAndRecycle(curr);

        if (answer) {
          answer = 0;
          // Since there can NEVER be more than one match
          // in the table, we're done!
          break;
        }
    } // if discard

    if (answer) {
        //
        // "Hit"
        //
        if (CHAINED) {
          //
          // If the matching entry is not already at the front of the list,
          // move it there
          //
          if (preventry) {
            preventry[0].UL = entry[0].UL;
            entry[0].UL = table[hcurr];
            table[hcurr] = curr;
          }
        }
#ifdef DEBUG_CT
        printf("Found CT entry ");
        FILE_output out(stdout);
        showEntry(out, curr);
        printf(" handle %lu in slot %lu\n", curr, hcurr);
        fflush(stdout);
#endif
        break;
    } // if equal

    //
    // Advance
    //
    if (CHAINED) {
      curr = entry[0].UL;
      preventry = entry;
    } else {
      if (++chain > maxCollisionSearch) break;
      incMod(hcurr);
      curr = table[hcurr];
    }
    
  } // for chain

  sawSearch(chain);

  return answer;
}

// **********************************************************************

// HERE

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>
::find(entry_key *key, entry_result& res)
{
  //
  // Allocate temporary space for key preprocessing.
  //
  const entry_type* et = key->getET();
  MEDDLY_DCASSERT(et);
  unsigned temp_bytes = 
    et->getKeySize( key->numRepeats() )
    +
    (MONOLITHIC ? 1 : 0)
    +
    (et->isRepeating() ? 1 : 0);
  temp_bytes *= sizeof(entry_item);

  entry_item* temp_entry = (entry_item*) key->allocTempData(temp_bytes);

  //
  // Copy the operation index if we're monolithic
  //
  int tptr = 0;
  if (MONOLITHIC) {
    temp_entry[0].UL = 0;
    temp_entry[0].U = et->getID();
    tptr++;
  }

  //
  // Copy the number of repeats if we're a repeating entry
  //
  if (et->isRepeating()) {
    temp_entry[tptr].UL = 0;
    temp_entry[tptr].U = key->numRepeats();
    tptr++;
  }

  //
  // Copy the key into temp_entry
  //
  memcpy(temp_entry+tptr, key->rawData(), key->dataLength() * sizeof(entry_item) );

  // 
  // TBD - probably shouldn't hash the floats, doubles.
  //
  
  //
  // Hash the key
  //
  setHash(key, raw_hash(temp_entry, temp_bytes));

  entry_item* entry_result = findEntry(key);
  perf.pings++;

  if (entry_result) {
    perf.hits++;
    //
    // Fill res
    //
    res.reset();
    res.setValid(entry_result); // TBD THIS
  } else {
    res.setInvalid();
  }
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::addEntry(entry_key* key, const entry_result &res)
{
  MEDDLY_DCASSERT(key);
  if (!MONOLITHIC) {
    if (key->getET() != global_et)
      throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
  }

  //
  // Increment cache counters for nodes in the key and result
  //
  key->cacheNodes();
  res.cacheNodes();

  unsigned long h = key->getHash() % tableSize;

  const entry_type* et = key->getET();
  MEDDLY_DCASSERT(et);

  //
  // Allocate an entry
  //

  const unsigned num_slots = 
    ( key->numTempBytes() / sizeof(entry_item) )
    + key->getET()->getResultSize()
    + (CHAINED ? 1 : 0)
  ;
  unsigned long curr = newEntry(num_slots);

#ifdef INTEGRATED_MEMMAN
  entry_item* entry = entries + curr;
#else
  entry_item* entry = (entry_item*) MMAN->getChunkAddress(curr);
#endif

  //
  // Copy into the entry
  //
  entry_item* key_portion = entry + (CHAINED ? 1 : 0);
  memcpy(key_portion, key->readTempData(), key->numTempBytes());
  entry_item* res_portion = key_portion + key->numTempBytes() / sizeof(entry_item);
  setResult(res_portion, res, et);

  //
  // Recycle key
  //
  recycle(key);


  //
  // Add entry to CT
  //
  if (CHAINED) {
    // Add this to front of chain
    entry[0].UL = table[h];
    table[h] = curr;
  } else {
    setTable(h, curr);
  }

#ifdef DEBUG_CT
  printf("Added CT entry ");
  FILE_output out(stdout);
  showEntry(out, curr);
  printf(" handle %lu in slot %lu (%u slots long)\n", curr, h, num_slots);
  fflush(stdout);
#endif

  if (perf.numEntries < tableExpand) return;

  //
  // Time to GC and maybe resize the table
  //

#ifdef DEBUG_SLOW
  fprintf(stdout, "Running GC in compute table (size %d, entries %u)\n", 
    tableSize, perf.numEntries
  );
#endif

  unsigned long list = 0;
  if (CHAINED) {
    list = convertToList(checkStalesOnResize);
    if (perf.numEntries < tableSize) {
      listToTable(list);
#ifdef DEBUG_SLOW
      fprintf(stdout, "Done CT GC, no resizing (now entries %u)\n", perf.numEntries);
#endif
      return;
    }
  } else {
    scanForStales();
    if (perf.numEntries < tableExpand / 4) {
#ifdef DEBUG_SLOW
      fprintf(stdout, "Done CT GC, no resizing (now entries %u)\n", perf.numEntries);
#endif
      return;
    }
  }

  unsigned long newsize = tableSize * 2;
  if (newsize > maxSize) newsize = maxSize;

  if (CHAINED) {
    if (newsize != tableSize) {
      //
      // Enlarge table
      //
    
      unsigned long* newt = (unsigned long*) realloc(table, newsize * sizeof(unsigned long));
      if (0==newt) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      }

      for (unsigned long i=tableSize; i<newsize; i++) newt[i] = 0;

      MEDDLY_DCASSERT(newsize > tableSize);
      mstats.incMemUsed( (newsize - tableSize) * sizeof(unsigned long) );
      mstats.incMemAlloc( (newsize - tableSize) * sizeof(unsigned long) );

      table = newt;
      tableSize = newsize;
    }

    if (tableSize == maxSize) {
      tableExpand = INT_MAX;      // TBD THIS
    } else {
      tableExpand = 4*tableSize;
    }
    tableShrink = tableSize / 2;

    listToTable(list);
  } else {  // not CHAINED
    if (newsize != tableSize) {
      //
      // Enlarge table
      //
      unsigned long* oldT = table;
      unsigned long oldSize = tableSize;
      tableSize = newsize;
      table = (unsigned long*) malloc(newsize * sizeof(unsigned long));
      if (0==table) {
        table = oldT;
        tableSize = oldSize;
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      }
      for (unsigned long i=0; i<newsize; i++) table[i] = 0;

      mstats.incMemUsed(newsize * sizeof(unsigned long));
      mstats.incMemAlloc(newsize * sizeof(unsigned long));

      rehashTable(oldT, oldSize);
      free(oldT);

      mstats.decMemUsed(oldSize * sizeof(unsigned long));
      mstats.decMemAlloc(oldSize * sizeof(unsigned long));

      if (tableSize == maxSize) {
        tableExpand = INT_MAX;    // TBD THIS
      } else {
        tableExpand = tableSize / 2;
      }
      tableShrink = tableSize / 8;
    }
  } // if CHAINED

#ifdef DEBUG_SLOW
  fprintf(stdout, "CT enlarged to size %lu\n", tableSize);
#endif

}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::updateEntry(entry_key* key, const entry_result &res)
{
  MEDDLY_DCASSERT(key->getET()->isResultUpdatable());
  entry_item* entry_result = findEntry(key);
  if (!entry_result) {
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }

  //
  // decrement cache counters for old result,
  //

  const entry_type* et = key->getET();
  for (unsigned i=0; i<et->getResultSize(); i++) {
#ifdef DEVELOPMENT_CODE
    typeID t;
    expert_forest* f;
    et->getResultType(i, t, f);
#else
    expert_forest* f = et->getResultForest(i);
#endif
    if (f) {
      MEDDLY_DCASSERT(NODE==t);
      f->uncacheNode( entry_result[i].N );
    }
  } // for i

  //
  // increment cache counters for new result.
  //
  res.cacheNodes();

  //
  // Overwrite result
  //
  setResult(entry_result, res, key->getET());
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::removeStales()
{
#ifdef DEBUG_REMOVESTALES
  fprintf(stdout, "Removing stales in CT (size %d, entries %u)\n", 
        tableSize, perf.numEntries
  );
#endif

  if (CHAINED) {

    //
    // Chained 
    //
    unsigned long list = convertToList(true);

    if (perf.numEntries < tableShrink) {
      //
      // Time to shrink table
      //
      unsigned long newsize = tableSize / 2;
      if (newsize < 1024) newsize = 1024;
      unsigned long* newt = (unsigned long*) realloc(table, newsize * sizeof(unsigned long));
      if (0==newt) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__); 
      }

      MEDDLY_DCASSERT(tableSize > newsize);
      mstats.decMemUsed( (tableSize - newsize) * sizeof(unsigned long) );
      mstats.decMemAlloc( (tableSize - newsize) * sizeof(unsigned long) );

      table = newt;
      tableSize = newsize;
      tableExpand = 4*tableSize;
      if (1024 == tableSize) {
        tableShrink = 0;
      } else {
        tableShrink = tableSize / 2;
      }
    }

    listToTable(list);

  } else {

    //
    // Unchained
    //

    scanForStales();

    if (perf.numEntries < tableShrink) {
      //
      // Time to shrink table
      //
      unsigned long newsize = tableSize / 2;
      if (newsize < 1024) newsize = 1024;
      if (newsize < tableSize) {
          unsigned long* oldT = table;
          unsigned long oldSize = tableSize;
          tableSize = newsize;
          table = (unsigned long*) malloc(newsize * sizeof(unsigned long));
          if (0==table) {
            table = oldT;
            tableSize = oldSize;
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
          }
          for (unsigned long i=0; i<newsize; i++) table[i] = 0;
          mstats.incMemUsed(newsize * sizeof(unsigned long));
          mstats.incMemAlloc(newsize * sizeof(unsigned long));
    
          rehashTable(oldT, oldSize);
          free(oldT);
      
          mstats.decMemUsed(oldSize * sizeof(unsigned long));
          mstats.decMemAlloc(oldSize * sizeof(unsigned long));
    
          tableExpand = tableSize / 2;
          if (1024 == tableSize) {
            tableShrink = 0;
          } else {
            tableShrink = tableSize / 8;
          }
      } // if different size
    }
  } // if CHAINED
    
#ifdef DEBUG_REMOVESTALES
  fprintf(stdout, "Done removing CT stales (size %lu, entries %lu)\n", 
    tableSize, perf.numEntries
  );
  FILE_output out(stderr);
  out << "CT after removing stales:\n";
  show(out, 9);
  fflush(stdout);
#endif
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::removeAll()
{
#ifdef DEBUG_REMOVESTALES
  fprintf(stdout, "Removing all CT entries\n");
#endif
  for (unsigned i=0; i<tableSize; i++) {
    while (table[i]) {
      int curr = table[i];
#ifdef INTEGRATED_MEMMAN
      const entry_item* entry = entries + curr;
#else
      const entry_item* entry = (const int*) MMAN->getChunkAddress(curr);
#endif
      if (CHAINED) {
        table[i] = entry[0].UL;
      } else {
        table[i] = 0;
      }
      discardAndRecycle(curr);
    } // while
  } // for i
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>
::show(output &s, int verbLevel)
{
  if (verbLevel < 1) return;

  if (MONOLITHIC) {
    s << "Monolithic compute table\n";
  } else {
    s << "Compute table for " << global_et->getName() << " (index " 
      << long(global_et->getID()) << ")\n";
  }

  s.put("", 6);
  s << "Current CT memory   :\t" << mstats.getMemUsed() << " bytes\n";
  s.put("", 6);
  s << "Peak    CT memory   :\t" << mstats.getPeakMemUsed() << " bytes\n";
  s.put("", 6);
  s << "Current CT alloc'd  :\t" << mstats.getMemAlloc() << " bytes\n";
  s.put("", 6);
  s << "Peak    CT alloc'd  :\t" << mstats.getPeakMemAlloc() << " bytes\n";
  if (!CHAINED) {
    s.put("", 6);
    s << "Collisions          :\t" << long(collisions) << "\n";
  }
  s.put("", 6);
  s << "Hash table size     :\t" << long(tableSize) << "\n";
  s.put("", 6);
  s << "Number of entries   :\t" << long(perf.numEntries) << "\n";

  if (--verbLevel < 1) return;

  s.put("", 6);
  s << "Pings               :\t" << long(perf.pings) << "\n";
  s.put("", 6);
  s << "Hits                :\t" << long(perf.hits) << "\n";

  if (--verbLevel < 1) return;

  s.put("", 6);
  s << "Search length histogram:\n";
  for (int i=0; i<stats::searchHistogramSize; i++) {
    if (perf.searchHistogram[i]) {
      s.put("", 10);
      s.put(long(i), 3);
      s << ": " << long(perf.searchHistogram[i]) << "\n";
    }
  }
  if (perf.numLargeSearches) {
    s.put("", 6);
    s << "Searches longer than " << long(stats::searchHistogramSize-1)
      << ": " << long(perf.numLargeSearches) << "\n";
  }
  s.put("", 6);
  s << "Max search length: " << long(perf.maxSearchLength) << "\n";


  if (--verbLevel < 1) return;

  s << "Hash table:\n";

  for (unsigned i=0; i<tableSize; i++) {
    if (0==table[i]) continue;
    s << "table[";
    s.put(long(i), 9);
    s << "]: ";
    showChain(s, table[i]);
  }

  if (--verbLevel < 1) return;

  s << "\nHash table nodes:\n";
  
  for (unsigned i=0; i<tableSize; i++) {
    int curr = table[i];
    while (curr) {
      s << "\tNode ";
      s.put(long(curr), 9);
      s << ":  ";
      showEntry(s, curr);
      s.put('\n');
      if (CHAINED) {
#ifdef INTEGRATED_MEMMAN
        curr = entries[curr].UL;
#else
        const int* entry = (const int*) MMAN->getChunkAddress(curr);
        curr = entry[0].UL;
#endif
      } else {
        curr = 0;
      }
    } // while curr
  } // for i
  s.put('\n');

  if (--verbLevel < 1) return;

#ifdef INTEGRATED_MEMMAN
  if (0==freeList) {
    s << "Free: null\n";
  } else {
    for (int i=1; i<=maxEntrySize; i++) {
      if (freeList[i]) {
        s << "freeList[" << i << "]: ";

        int L = freeList[i];
        s << L;
        while (L) {
          L = entries[L].UL;
          s << "->" << L;
        }
        s << "\n";
      }
    }
  }

  if (0==entries) {
    s << "Entries: null\n";
  } else {
    s << "Entries: [" << long(entries[0].L);
    for (int i=1; i<entriesSize; i++) {
      s << ", " << long(entries[i].L);
    }
    s << "]\n";
  }
#else
  if (MMAN) MMAN->dumpInternal(s);
#endif

}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
unsigned long MEDDLY::ct_none<MONOLITHIC, CHAINED>::convertToList(bool removeStales)
{
  MEDDLY_DCASSERT(CHAINED);
  // const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

  unsigned long list = 0;
  for (unsigned long i=0; i<tableSize; i++) {
    while (table[i]) {
      unsigned long curr = table[i];
#ifdef INTEGRATED_MEMMAN
      entry_item* entry = entries + curr;
#else
      entry_item* entry = (entry_item*) MMAN->getChunkAddress(curr);
#endif
      table[i] = entry[0].UL;
      if (removeStales && isStale(entry)) {

#ifdef DEBUG_TABLE2LIST
          printf("\tstale ");
          FILE_output out(stdout);
          showEntry(out, curr);
          printf(" (handle %lu table slot %lu)\n", curr, i);
#endif

          discardAndRecycle(curr);
          continue;
        // }
      } // if removeStales
      //
      // Not stale, move to list
      //
#ifdef DEBUG_TABLE2LIST
      printf("\tkeep  ");
      FILE_output out(stdout);
      showEntry(out, curr);
      printf(" (handle %lu table slot %lu)\n", curr, i);
#endif
      entry[0].UL = list;
      list = curr;
    } // while
  } // for i
#ifdef DEBUG_TABLE2LIST
  printf("Built list: ");
  showChain(stdout, list);
#endif
  return list;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::listToTable(unsigned long L)
{
  const int M = MONOLITHIC ? 1 : 0;
  MEDDLY_DCASSERT(CHAINED);
#ifdef DEBUG_LIST2TABLE
  printf("Recovering  list");
  // showChain(stdout, L);
#endif
  while (L) {
#ifdef INTEGRATED_MEMMAN
    entry_item* entry = entries+L;
#else
    entry_item* entry = (entry_item*) MMAN->getChunkAddress(L);
#endif
    const unsigned long curr = L;
    L = entry[0].UL;

    const entry_type* et = MONOLITHIC
      ?   getEntryType(entry[1].U)
      :   global_et;
    MEDDLY_DCASSERT(et);
    const unsigned reps = (et->isRepeating()) ? entry[ M+1 ].U : 0;
    // HERE
    const unsigned hashlength = M + (et->isRepeating() ? 1 : 0) + ( et->getKeyBytes(reps) / sizeof(entry_item) );  

    const unsigned h = hash(entry + 1, hashlength * sizeof(entry_item));
    entry[0].UL = table[h];
    table[h] = curr;
#ifdef DEBUG_LIST2TABLE
    printf("\tsave  ");
    FILE_output out(stdout);
    showEntry(out, curr);
    printf(" (handle %lu table slot %lu hashlength %u)\n", curr, h, hashlength);
#endif
  }
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::scanForStales()
{
  MEDDLY_DCASSERT(!CHAINED);
  for (unsigned long i=0; i<tableSize; i++) {
    if (0==table[i]) continue;

#ifdef INTEGRATED_MEMMAN
    const entry_item* entry = entries + table[i];
#else
    const entry_item* entry = (const entry_item*) MMAN->getChunkAddress(table[i]);
#endif

    if (isStale(entry)) {

#ifdef DEBUG_CT
      printf("Removing CT stale entry ");
      FILE_output out(stdout);
      showEntry(out, table[i]);
      printf(" in table slot %lu\n", i);
#endif  

      discardAndRecycle(table[i]);
      table[i] = 0;

    } // if isStale
  } // for i
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::rehashTable(const unsigned long* oldT, unsigned long oldS)
{
#ifdef DEBUG_REHASH
  printf("Rebuilding hash table\n");
#endif
  MEDDLY_DCASSERT(!CHAINED);
  const int M = MONOLITHIC ? 1 : 0;
  for (unsigned long i=0; i<oldS; i++) {
    unsigned long curr = oldT[i];
    if (0==curr) continue;
#ifdef INTEGRATED_MEMMAN
    const entry_item* entry = entries + curr;
#else
    const entry_item* entry = (const entry_item*) MMAN->getChunkAddress(curr);
#endif

    const entry_type* et = MONOLITHIC
      ?   getEntryType(entry[0].UL)
      :   global_et;
    MEDDLY_DCASSERT(et);
    const unsigned reps = (et->isRepeating()) ? entry[ M ].U : 0;
    const unsigned hashlength = M + (et->isRepeating() ? 1 : 0) + ( et->getKeyBytes(reps) / sizeof(entry_item) );

    unsigned long h = hash(entry, hashlength*sizeof(entry_item));
    setTable(h, curr);

#ifdef DEBUG_REHASH
    printf("\trehash  ");
    FILE_output out(stdout);
    showEntry(out, curr);
    printf(" (handle %lu table slot %lu hashlength %u)\n", curr, h, hashlength);
#endif
  }
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
MEDDLY::node_address MEDDLY::ct_none<MONOLITHIC, CHAINED>
::newEntry(unsigned size)
{
#ifdef INTEGRATED_MEMMAN
  // check free list
  if (size > maxEntrySize) {
    fprintf(stderr, "MEDDLY error: request for compute table entry larger than max size\n");
    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);  // best we can do
  }
  if (size<1) return 0;
  perf.numEntries++;
  if (freeList[size]) {
    node_address h = freeList[size];
    freeList[size] = entries[h].UL;
#ifdef DEBUG_CTALLOC
    fprintf(stderr, "Re-used entry %ld size %d\n", h, size);
#endif
    mstats.incMemUsed( size * sizeof(entry_item) );
    return h;
  }
  if (entriesSize + size > entriesAlloc) {
    // Expand by a factor of 1.5
    unsigned neA = entriesAlloc + (entriesAlloc/2);
    entry_item* ne = (entry_item*) realloc(entries, neA * sizeof(entry_item));
    if (0==ne) {
      fprintf(stderr,
          "Error in allocating array of size %lu at %s, line %d\n",
          neA * sizeof(int), __FILE__, __LINE__);
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    mstats.incMemAlloc( (entriesAlloc / 2) * sizeof(entry_item) );
    entries = ne;
    entriesAlloc = neA;
  }
  MEDDLY_DCASSERT(entriesSize + size <= entriesAlloc);
  node_address h = entriesSize;
  entriesSize += size;
#ifdef DEBUG_CTALLOC
  fprintf(stderr, "New entry %lu size %u\n", h, size);
#endif
  mstats.incMemUsed( size * sizeof(entry_item) );
  return h;

#else
  perf.numEntries++;
  size_t the_size = size;
  return MMAN->requestChunk(the_size);
#endif
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
bool MEDDLY::ct_none<MONOLITHIC, CHAINED> 
::isStale(const entry_item* entry) const
{
  const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

  //
  // Check entire entry for staleness
  //
  const entry_type* et = MONOLITHIC
      ?   getEntryType(entry[CHAINED?1:0].U)
      :   global_et;
  MEDDLY_DCASSERT(et);

  if (et->isMarkedForDeletion()) return true;

  entry += SHIFT;
  const unsigned reps = (et->isRepeating()) ? (*entry++).U : 0;
  const unsigned klen = et->getKeySize(reps);

  //
  // Key portion
  //
  for (unsigned i=0; i<klen; i++) {
#ifdef DEVELOPMENT_CODE
    typeID t;
    expert_forest* f;
    et->getKeyType(i, t, f);
#else
    expert_forest* f = et->getKeyForest(i);
#endif
    if (f) {
      MEDDLY_DCASSERT(NODE==t);
      if (MEDDLY::forest::ACTIVE != f->getNodeStatus(entry[i].N)) {
        return true;
      }
    }
  } // for i
  entry += klen;

  // 
  // Result portion
  //
  for (unsigned i=0; i<et->getResultSize(); i++) {
#ifdef DEVELOPMENT_CODE
    typeID t;
    expert_forest* f;
    et->getResultType(i, t, f);
#else
    expert_forest* f = et->getResultForest(i);
#endif
    if (f) {
      MEDDLY_DCASSERT(NODE==t);
      if (MEDDLY::forest::ACTIVE != f->getNodeStatus(entry[i].N)) {
        return true;
      }
    }
  } // for i

  return false;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
bool MEDDLY::ct_none<MONOLITHIC, CHAINED> 
::isDead(const entry_item* result, const entry_type* et) const
{
  MEDDLY_DCASSERT(et);
  MEDDLY_DCASSERT(result);
  //
  // Check result portion for dead nodes - cannot use result in that case
  //
  for (unsigned i=0; i<et->getResultSize(); i++) {
#ifdef DEVELOPMENT_CODE
    typeID t;
    expert_forest* f;
    et->getResultType(i, t, f);
#else
    expert_forest* f = et->getResultForest(i);
#endif
    if (f) {
      MEDDLY_DCASSERT(NODE == t);
      if (MEDDLY::forest::DEAD == f->getNodeStatus(result[i].N)) {
        return true;
      }
    }
  } // for i
  return false;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>
::discardAndRecycle(unsigned long h)
{
  const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

#ifdef INTEGRATED_MEMMAN
  entry_item* entry = entries + h;
#else
  entry_item* entry = (entry_item*) MMAN->getChunkAddress(table[h]);
#endif

  unsigned slots;

  //
  // Discard.
  // For the portion of the entry that are nodes,
  // we need to notify the forest to decrement
  // the CT counter.
  //

  const entry_type* et = MONOLITHIC
    ?   getEntryType(entry[CHAINED ? 1 : 0].U)
    :   global_et;
  MEDDLY_DCASSERT(et);

  const entry_item* ptr = entry + SHIFT;
  unsigned reps;
  if (et->isRepeating()) {
    reps = (*ptr).U;
    ptr++;
  } else {
    reps = 0;
  }

  //
  // Key portion
  //
  const unsigned stop = et->getKeySize(reps);
  for (unsigned i=0; i<stop; i++) {
    typeID t;
    expert_forest* f;
    et->getKeyType(i, t, f);
    if (f) {
      MEDDLY_DCASSERT(NODE == t);
      f->uncacheNode( ptr[i].N );
      continue;
    }
    if (POINTER == t) {
      ct_object* P = (ct_object*) ptr[i].P;
      delete P;
      continue;
    }
    MEDDLY_DCASSERT( t != NODE );
  } // for i
  ptr += stop;

  //
  // Result portion
  //
  for (unsigned i=0; i<et->getResultSize(); i++) {
    typeID t;
    expert_forest* f;
    et->getResultType(i, t, f);
    if (f) {
      MEDDLY_DCASSERT(NODE == t);
      f->uncacheNode( ptr[i].N );
      continue;
    }
    if (POINTER == t) {
      ct_object* P = (ct_object*) ptr[i].P;
      delete P;
      continue;
    }
    MEDDLY_DCASSERT( t != NODE );
  } // for i
  
  slots = SHIFT + (et->isRepeating() ? 1 : 0);
  slots += stop;
  slots += et->getResultSize();

  //
  // Recycle
  //
#ifdef DEBUG_CTALLOC
  fprintf(stderr, "Recycling entry %lu size %u\n", h, slots);
#endif
#ifdef INTEGRATED_MEMMAN
  entries[h].UL = freeList[slots];
  freeList[slots] = h;
  mstats.decMemUsed( slots * sizeof(entry_item) );
#else
  MMAN->recycleChunk(h, slots);
#endif
  perf.numEntries--;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>
::showEntry(output &s, unsigned long h) const
{
#ifdef INTEGRATED_MEMMAN
  const entry_item* entry = entries+h;
#else
  const entry_item* entry = (const entry_item*) MMAN->getChunkAddress(h);
#endif

  const entry_type* et = MONOLITHIC
    ?   getEntryType(entry[CHAINED?1:0].U)
    :   global_et;
  MEDDLY_DCASSERT(et);

  const entry_item* ptr = entry + (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);
  unsigned reps;
  if (et->isRepeating()) {
    reps = (*ptr).U;
    ptr++;
  } else {
    reps = 0;
  }
  s << "[" << et->getName() << "(";
  unsigned stop = et->getKeySize(reps);
  for (unsigned i=0; i<stop; i++) {
      if (i) s << ", ";
      switch (et->getKeyType(i)) {
        case NODE:
                        s.put(long(ptr[i].N));
                        break;
        case INTEGER:
                        s.put(long(ptr[i].I));
                        break;
        case LONG:
                        s.put(ptr[i].L);
                        break;
        case FLOAT:
                        s.put(ptr[i].F, 0, 0, 'e');
                        break;
        case DOUBLE:
                        s.put(ptr[i].D, 0, 0, 'e');
                        break;
        case POINTER:
                        s.put_hex((unsigned long)ptr[i].P);
                        break;
        default:
                        MEDDLY_DCASSERT(0);
      } // switch et->getKeyType()
  } // for i
  ptr += stop;
  s << "): ";
  for (unsigned i=0; i<et->getResultSize(); i++) {
      if (i) s << ", ";
      switch (et->getResultType(i)) {
        case NODE:
                        s.put(long(ptr[i].N));
                        break;
        case INTEGER:
                        s.put(long(ptr[i].I));
                        break;
        case LONG:
                        s.put(ptr[i].L);
                        s << "(L)";
                        break;
        case FLOAT:
                        s.put(ptr[i].F, 0, 0, 'e');
                        break;
        case DOUBLE:
                        s.put(ptr[i].D, 0, 0, 'e');
                        break;
                            
        case POINTER:
                        s.put_hex((unsigned long)ptr[i].P);
                        break;
        default:
                        MEDDLY_DCASSERT(0);
      } // switch et->getResultType()
  } // for i
  s << "]";
}


// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>
::showKey(output &s, const entry_key* key) const
{
  MEDDLY_DCASSERT(key);

  const entry_type* et = key->getET();
  MEDDLY_DCASSERT(et);

  unsigned reps = key->numRepeats();
  s << "[" << et->getName() << "(";
  unsigned stop = et->getKeySize(reps);

  for (unsigned i=0; i<stop; i++) {
      entry_item item = key->rawData()[i];
      if (i) s << ", ";
      switch (et->getKeyType(i)) {
        case NODE:
                        s.put(long(item.N));
                        break;
        case INTEGER:
                        s.put(long(item.I));
                        break;
        case LONG:
                        s.put(item.L);
                        break;
        case FLOAT:
                        s.put(item.F);
                        break;
        case DOUBLE:
                        s.put(item.D);
                        break;
        case POINTER:
                        s.put_hex((unsigned long)item.P);
                        break;
        default:
                        MEDDLY_DCASSERT(0);
      } // switch et->getKeyType()
  } // for i
  s << "): ?]";
}


#endif  // include guard


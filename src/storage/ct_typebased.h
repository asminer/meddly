
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


#ifndef CT_TYPEBASED_H
#define CT_TYPEBASED_H

// #define EQUAL_USES_MEMCPY

#define USE_EQUAL_SW

// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                         ct_typebased class                         *
// *                                                                    *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {
  template <bool MONOLITHIC, bool CHAINED>
  class ct_typebased : public compute_table {
    public:
      ct_typebased(const ct_initializer::settings &s, operation* op, unsigned slot);
      virtual ~ct_typebased();

      /**
          Find an entry.
          Used by find() and updateEntry().
            @param  key   Key to search for.
            @return Pointer to the result portion of the entry, or null if not found.
      */
      int* findEntry(entry_key* key);

      // required functions

      virtual bool isOperationTable() const   { return !MONOLITHIC; }
      virtual void find(entry_key* key, entry_result &res);
      virtual void addEntry(entry_key* key, const entry_result& res);
      virtual void updateEntry(entry_key* key, const entry_result& res);
      virtual void removeStales();
      virtual void removeAll();
      virtual void show(output &s, int verbLevel = 0);

    private:  // helper methods

      int convertToList(bool removeStales);
      void listToTable(int h);

      void scanForStales();
      void rehashTable(const int* oldT, unsigned oldS);

      /// Grab space for a new entry
      node_address newEntry(unsigned size);

      static inline unsigned raw_hash(const int* k, int length) {
        return hash_stream::raw_hash( (const unsigned*) k, length );
      }

      inline unsigned hash(const int* k, int length) const {
        return raw_hash(k, length) % tableSize;
      }
    
      inline void incMod(unsigned &h) {
        h++;
        if (h>=tableSize) h=0;
      }

      /// Update stats: we just searched through c items
      inline void sawSearch(int c) {
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
      inline void setTable(unsigned h, int curr) {
        MEDDLY_DCASSERT(!CHAINED);
        unsigned hfree = h;
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
#ifdef DEBUG_CT_SLOTS
        printf(" handle %d in slot %u\n", table[h], h);
#else
        printf("\n");
#endif
        fflush(stdout);
#endif
        collisions++;    

        discardAndRecycle(table[h]);
        table[h] = curr;
      }


      /// Check equality
      static inline bool equal_sw(const int* a, const int* b, int N) {
#ifdef EQUAL_USES_MEMCPY
        return (0==memcmp(a, b, N*sizeof(int)));
#else
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
          default:    return (0==memcmp(a, b, N*sizeof(int)));
        };
#endif
      }



      /** 
        Check equality.
        We advance pointer a and return the result portion if equal.
        If unequal we return 0.
      */
      static inline int* equal(int* a, const entry_key* key) {
        const entry_type* et = key->getET();
        MEDDLY_DCASSERT(et);
        //
        // Compare operator, if needed
        //
        if (MONOLITHIC) {
          if (et->getID() != *a) return 0;
          a++;
        }
        //
        // Compare #repetitions, if needed
        //
        if (et->isRepeating()) {
          if (key->numRepeats() != *a) return 0;
          a++;
        }
        //
        // Compare key portion
        //
        const entry_item* b = key->rawData();
        const unsigned klen = et->getKeySize(key->numRepeats());
        for (unsigned i=0; i<klen; i++) {
          const typeID t = et->getKeyType(i);
          switch (t) {
              case FLOAT:   // Hack: treat like integer
              case INTEGER:
                        if (b[i].I != *a) return 0;
                        a++;
                        continue;
              case NODE:
                        if (b[i].N != *a) return 0;
                        a++;
                        continue;
              case LONG:
                        if (b[i].L != *((long*)a) ) return 0;
                        a += 2;
                        continue;
              default:
                        // Other cases shouldn't be in the key
                        MEDDLY_DCASSERT(0);
          } // switch
        } // for i
        return a;
      }

      /**
          Check if the key portion of an entry equals key and should be discarded.
            @param  entry   Complete entry in CT to check.
            @param  key     Key to compare against. 
            @param  discard On output: should this entry be discarded

            @return Result portion of the entry, if they key portion matches;
                    0 if the entry does not match the key.
      */
      inline int* checkEqualityAndStatus(int* entry, const entry_key* key, bool &discard)
      {
        int* entry_without_next = CHAINED ? (entry+1) : entry;
#ifdef USE_EQUAL_SW
        const unsigned keyslots = key->numTempBytes() / sizeof(int);
        if (equal_sw(entry_without_next, (const int*) key->readTempData(), keyslots))
        {
          //
          // Equal.
          //
          int* result = entry_without_next + keyslots;
#else // USE_EQUAL_SW
        int* result = equal(entry_without_next, key);
        if (result) 
        {
#endif  // USE_EQUAL_SW

#ifdef DEBUG_ISDEAD
          printf("Checking entry result for deadness: ");
          FILE_output out(stdout);
          showEntry(out, entry);
          printf("\n");
          fflush(stdout);
#endif
          discard = isDead(result, key->getET());
#ifdef DEBUG_ISDEAD
          if (discard) {
            printf("\tdead\n");
          } else {
            printf("\tnot dead\n");
          }
          fflush(stdout);
#endif
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
      bool isStale(const int* entry) const;

      /**
          Check if a result is dead (unrecoverable).
            @param  res   Pointer to result portion of entry.
            @param  et    Entry type information.
            @return       true, if the entry should be discarded.
      */
      bool isDead(const int* res, const entry_type* et) const;


      /**
          Copy a result into an entry.
      */
      inline void setResult(int* respart, const entry_result &res, const entry_type* et) {
        const entry_item* resdata = res.rawData();
        for (unsigned i=0; i<res.dataLength(); i++) {
            typeID t = et->getResultType(i);
            switch (t) {
              case NODE:    *respart = resdata[i].N;
                            respart++;
                            continue;

              case INTEGER: *respart = resdata[i].I;
                            respart++;
                            continue;

              case FLOAT:   *((float*)respart) = resdata[i].F;
                            respart++;
                            continue;

              case DOUBLE:  *((double*)respart) = resdata[i].D;
                            respart += sizeof(double) / sizeof(int);
                            continue;

              case LONG:    *((long*)respart) = resdata[i].L;
                            respart += sizeof(long) / sizeof(int);
                            continue;

              case POINTER: *((void**)respart) = resdata[i].P;
                            respart += sizeof(void*) / sizeof(int);
                            continue;

              default:      MEDDLY_DCASSERT(0);
            } // switch t
        } // for i
      }

      /// Display a chain 
      inline void showChain(output &s, int L) const {
        s << L;
        if (CHAINED) {
          while (L) {
#ifdef INTEGRATED_MEMMAN
            const int* entry = entries+L;
#else
            const int* entry = (const int*) MMAN->getChunkAddress(L);
#endif
            L = entry[0];
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
            @param  entry     Pointer to the entry.
      */
      void showEntry(output &s, const int* entry) const;

      /**
          Display an entry.
          Used for debugging, and by method(s) to display the entire CT.
            @param  s         Stream to write to.
            @param  h         Handle of the entry.
      */
      inline void showEntry(output &s, unsigned long h) const 
      {
#ifdef INTEGRATED_MEMMAN
        showEntry(s, entries+h);
#else
        showEntry(s, (const int*) MMAN->getChunkAddress(h));
#endif
      }

      /// Display a key.
      void showKey(output &s, const entry_key* k) const;


    private:
      /// Global entry type.  Ignored when MONOLITHIC is true.
      const entry_type* global_et;

      /// Hash table
      int*  table;

      /// Hash table size
      unsigned int tableSize;

      /// When to next expand the table
      unsigned int tableExpand;

      /// When to next shrink the table
      unsigned int tableShrink;

      /// Space to build an entry
      // our_temp_entry currEntry;

#ifdef INTEGRATED_MEMMAN
      /// Memory space for entries
      int*  entries;
      /// Used entries
      int entriesSize;
      /// Memory allocated for entries
      int entriesAlloc;

      static const int maxEntrySize = 15;
      static const int maxEntryBytes = sizeof(int) * maxEntrySize;

      /// freeList[i] is list of all unused i-sized entries.
      int* freeList;
#else
  
      memory_manager* MMAN;
#endif
  
      /// Memory statistics
      memstats mstats;

      /// How many slots to search in unchained tables
      static const int maxCollisionSearch = 2;

      /// Stats: how many collisions
      long collisions;
  }; // class ct_typebased
} // namespace


// **********************************************************************
// *                                                                    *
// *                        ct_typebased methods                        *
// *                                                                    *
// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
MEDDLY::ct_typebased<MONOLITHIC, CHAINED>::ct_typebased(
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
  freeList = new int[1+maxEntrySize];
  for (int i=0; i<=maxEntrySize; i++) {
    freeList[i] = 0;
  }
  mstats.incMemUsed( (1+maxEntrySize) * sizeof(int) );
  mstats.incMemAlloc( (1+maxEntrySize) * sizeof(int) );

  entriesAlloc = 1024;
  entriesSize = 1;    
  entries = (int*) malloc(entriesAlloc * sizeof(int) );
  if (0==entries) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  entries[0] = 0;     // NEVER USED; set here for sanity.

  mstats.incMemUsed( entriesSize * sizeof(int) );
  mstats.incMemAlloc( entriesAlloc * sizeof(int) );
#else
  MEDDLY_DCASSERT(s.MMS);
  MMAN = s.MMS->initManager(sizeof(int), 1, mstats);
#endif

  /*
      Initialize hash table
  */
  tableSize = 1024;
  tableExpand = CHAINED ? 4*1024 : 512;
  tableShrink = 0;
  table = (int*) malloc(tableSize * sizeof(int));
  if (0==table) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  for (unsigned i=0; i<tableSize; i++) table[i] = 0;

  mstats.incMemUsed(tableSize * sizeof(int));
  mstats.incMemAlloc(tableSize * sizeof(int));

  collisions = 0;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
MEDDLY::ct_typebased<MONOLITHIC, CHAINED>::~ct_typebased()
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
inline int* MEDDLY::ct_typebased<MONOLITHIC, CHAINED>
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

  int chain;

  int* answer = 0;
  int* preventry = 0;
  unsigned hcurr = key->getHash() % tableSize;
  int curr = table[hcurr];

#ifdef DEBUG_CT_SEARCHES
  printf("Searching for CT entry ");
  FILE_output out(stdout);
  showKey(out, key);
#ifdef DEBUG_CT_SLOTS
  printf(" in hash slot %u\n", hcurr);
#else
  printf("\n");
#endif
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
    int* entry = entries + curr;
#else
    int* entry = (int*) MMAN->getChunkAddress(curr);
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
#ifdef DEBUG_CT_SLOTS
        printf(" handle %d in slot %u\n", curr, hcurr);
#else
        printf("\n");
#endif
        fflush(stdout);
#endif
        if (CHAINED) {
          if (preventry) {
            preventry[0] = entry[0];
          } else {
            table[hcurr] = entry[0];
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
            preventry[0] = entry[0];
            entry[0] = table[hcurr];
            table[hcurr] = curr;
          }
        }
#ifdef DEBUG_CT
        printf("Found CT entry ");
        FILE_output out(stdout);
        showEntry(out, curr);
#ifdef DEBUG_CT_SLOTS
        printf(" handle %d in slot %u\n", curr, hcurr);
#else
        printf("\n");
#endif
        fflush(stdout);
#endif
        break;
    } // if equal

    //
    // Advance
    //
    if (CHAINED) {
      curr = entry[0];
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

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>
::find(entry_key *key, entry_result& res)
{
  //
  // Allocate temporary space for key preprocessing.
  //
  const entry_type* et = key->getET();
  MEDDLY_DCASSERT(et);
  unsigned temp_bytes =
    et->getKeyBytes( key->numRepeats() ) 
    +
    (MONOLITHIC ? sizeof(int) : 0)
    +
    (et->isRepeating() ? sizeof(int) : 0);

  int* temp_entry = (int*) key->allocTempData(temp_bytes);
  const entry_item* data = key->rawData();

  //
  // Copy the operation index if we're monolithic
  //
  int tptr = 0;
  if (MONOLITHIC) {
    temp_entry[0] = et->getID();
    tptr++;
  }

  //
  // Copy the number of repeats if we're a repeating entry
  //
  if (et->isRepeating()) {
    temp_entry[tptr] = key->numRepeats();
    tptr++;
  }

  //
  // Copy the key into temp_entry
  //
  unsigned datalen = key->dataLength();
  for (unsigned i=0; i<datalen; i++) {
    typeID t = et->getKeyType(i);
    switch (t) {
      case NODE:
                  temp_entry[tptr] = data[i].N;
                  tptr++;
                  continue;
      case INTEGER:
                  temp_entry[tptr] = data[i].I;
                  tptr++;
                  continue;
      case LONG: 
                  memcpy(temp_entry+tptr, &(data[i].L), sizeof(long));
                  tptr += sizeof(long) / sizeof(int);
                  continue;
      case FLOAT:
                  temp_entry[tptr] = data[i].I;   // Hack!
                  tptr++;
                  continue;
      default:
                  MEDDLY_DCASSERT(0);
    } // switch t
  } // for i

  // 
  // TBD - probably shouldn't hash the floats, doubles.
  //
  
  //
  // Hash the key
  //
  setHash(key, raw_hash(temp_entry, temp_bytes / sizeof(int)));


  int* entry_result = findEntry(key);
  perf.pings++;

  if (entry_result) {
    perf.hits++;
    //
    // Fill res
    //
    res.reset();
    for (unsigned i=0; i<key->getET()->getResultSize(); i++) {
      typeID t = et->getResultType(i);
      switch (t) {
        case NODE:    res.writeN( *entry_result++ );
                      continue;

        case INTEGER: res.writeI( *entry_result++ );
                      continue;

        case FLOAT:   res.writeF( *((float*)entry_result) );
                      entry_result++;
                      continue;

        case DOUBLE:  res.writeD( *((double*)entry_result) );
                      entry_result += sizeof(double) / sizeof(int);
                      continue;

        case LONG:    res.writeL( *((long*)entry_result) );
                      entry_result += sizeof(long) / sizeof(int);
                      continue;

        case POINTER: res.writeP( *((void**)entry_result) );
                      entry_result += sizeof(void*) / sizeof(int);
                      continue;
                    
        default:      MEDDLY_DCASSERT(0);
      } // switch t
    } // for i
    res.reset();
    res.setValid();
  } else {
    res.setInvalid();
  }
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>::addEntry(entry_key* key, const entry_result &res)
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

  unsigned h = key->getHash() % tableSize;

  const entry_type* et = key->getET();
  MEDDLY_DCASSERT(et);

  //
  // Allocate an entry
  //

  const unsigned num_slots = 
    ( key->numTempBytes() + key->getET()->getResultBytes() ) / sizeof(int)
    + (CHAINED ? 1 : 0)
  ;
  node_address curr = newEntry(num_slots);

#ifdef INTEGRATED_MEMMAN
  int* entry = entries + curr;
#else
  int* entry = (int*) MMAN->getChunkAddress(curr);
#endif

  //
  // Copy into the entry
  //
  int* key_portion = entry + (CHAINED ? 1 : 0);
  memcpy(key_portion, key->readTempData(), key->numTempBytes());
  int* res_portion = key_portion + key->numTempBytes() / sizeof(int);
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
    entry[0] = table[h];
    table[h] = curr;
  } else {
    setTable(h, curr);
  }

#ifdef DEBUG_CT
  printf("Added CT entry ");
  FILE_output out(stdout);
  showEntry(out, curr);
#ifdef DEBUG_CT_SLOTS
  printf(" handle %lu in slot %u (%u slots long)\n", curr, h, num_slots);
#else
  printf("\n");
#endif
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

  int list = 0;
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
    
      int* newt = (int*) realloc(table, newsize * sizeof(int));
      if (0==newt) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      }

      for (unsigned i=tableSize; i<newsize; i++) newt[i] = 0;

      MEDDLY_DCASSERT(newsize > tableSize);
      mstats.incMemUsed( (newsize - tableSize) * sizeof(int) );
      mstats.incMemAlloc( (newsize - tableSize) * sizeof(int) );

      table = newt;
      tableSize = newsize;
    }

    if (tableSize == maxSize) {
      tableExpand = INT_MAX;
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
      int* oldT = table;
      unsigned oldSize = tableSize;
      tableSize = newsize;
      table = (int*) malloc(newsize * sizeof(int));
      if (0==table) {
        table = oldT;
        tableSize = oldSize;
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      }
      for (unsigned i=0; i<newsize; i++) table[i] = 0;

      mstats.incMemUsed(newsize * sizeof(int));
      mstats.incMemAlloc(newsize * sizeof(int));

      rehashTable(oldT, oldSize);
      free(oldT);

      mstats.decMemUsed(oldSize * sizeof(int));
      mstats.decMemAlloc(oldSize * sizeof(int));

      if (tableSize == maxSize) {
        tableExpand = INT_MAX;
      } else {
        tableExpand = tableSize / 2;
      }
      tableShrink = tableSize / 8;
    }
  } // if CHAINED

#ifdef DEBUG_SLOW
  fprintf(stdout, "CT enlarged to size %d\n", tableSize);
#endif

}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>::updateEntry(entry_key* key, const entry_result &res)
{
  MEDDLY_DCASSERT(key->getET()->isResultUpdatable());
  int* entry_result = findEntry(key);
  if (!entry_result) {
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }

  //
  // decrement cache counters for old result,
  //

  const unsigned slots_for_type[] = { 1, 1, 1, 2, 1, 2, 2 };

  MEDDLY_DCASSERT(1 == slots_for_type[ERROR]);
  MEDDLY_DCASSERT(1 == slots_for_type[INTEGER]);
  MEDDLY_DCASSERT(sizeof(node_handle) / sizeof(int) == slots_for_type[NODE]);
  MEDDLY_DCASSERT(sizeof(long) / sizeof(int) == slots_for_type[LONG]);
  MEDDLY_DCASSERT(sizeof(float) / sizeof(int) == slots_for_type[FLOAT]);
  MEDDLY_DCASSERT(sizeof(double) / sizeof(int) == slots_for_type[DOUBLE]);
  MEDDLY_DCASSERT(sizeof(void*) / sizeof(int) == slots_for_type[POINTER]);

  const entry_type* et = key->getET();
  int* ptr = entry_result;
  for (unsigned i=0; i<et->getResultSize(); i++) {
    typeID t;
    expert_forest* f;
    et->getResultType(i, t, f);
    MEDDLY_CHECK_RANGE(0, t, 7);
    if (f) {
      MEDDLY_DCASSERT(NODE == t);
      f->uncacheNode( *ptr );
      ptr++;
    } else {
      MEDDLY_DCASSERT(NODE != t);
      ptr += slots_for_type[t];
    }
    /*
    switch (t) {
        case NODE:
                        MEDDLY_DCASSERT(f);
                        f->uncacheNode( *ptr );
                        ptr++;
                        continue;
        case INTEGER:
                        ptr++;
                        continue;
        case LONG:
                        ptr += sizeof(long) / sizeof(int);
                        continue;
        case FLOAT:
                        ptr += sizeof(float) / sizeof(int);
                        continue;
        case DOUBLE:
                        ptr += sizeof(double) / sizeof(int);
                        continue;
        case POINTER: {
                        ct_object* P = *((ct_object**)(ptr));
                        delete P;
                        ptr += sizeof(void*) / sizeof(int);
                        continue;
                      }
        default:
                        MEDDLY_DCASSERT(0);
    }
    */
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
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>::removeStales()
{
#ifdef DEBUG_REMOVESTALES
  fprintf(stdout, "Removing stales in CT (size %d, entries %lu)\n", 
        tableSize, perf.numEntries
  );
#endif

  if (CHAINED) {

    //
    // Chained 
    //
    int list = convertToList(true);

    if (perf.numEntries < tableShrink) {
      //
      // Time to shrink table
      //
      int newsize = tableSize / 2;
      if (newsize < 1024) newsize = 1024;
      int* newt = (int*) realloc(table, newsize * sizeof(int));
      if (0==newt) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__); 
      }

      MEDDLY_DCASSERT(tableSize > newsize);
      mstats.decMemUsed( (tableSize - newsize) * sizeof(int) );
      mstats.decMemAlloc( (tableSize - newsize) * sizeof(int) );

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
      unsigned newsize = tableSize / 2;
      if (newsize < 1024) newsize = 1024;
      if (newsize < tableSize) {
          int* oldT = table;
          unsigned oldSize = tableSize;
          tableSize = newsize;
          table = (int*) malloc(newsize * sizeof(int));
          if (0==table) {
            table = oldT;
            tableSize = oldSize;
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
          }
          for (unsigned i=0; i<newsize; i++) table[i] = 0;
          mstats.incMemUsed(newsize * sizeof(int));
          mstats.incMemAlloc(newsize * sizeof(int));
    
          rehashTable(oldT, oldSize);
          free(oldT);
      
          mstats.decMemUsed(oldSize * sizeof(int));
          mstats.decMemAlloc(oldSize * sizeof(int));
    
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
  fprintf(stdout, "Done removing CT stales (size %d, entries %lu)\n", 
    tableSize, perf.numEntries
  );
#ifdef DEBUG_REMOVESTALES_DETAILS
  FILE_output out(stderr);
  out << "CT after removing stales:\n";
  show(out, 9);
#endif
  fflush(stdout);
#endif
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>::removeAll()
{
#ifdef DEBUG_REMOVESTALES
  fprintf(stdout, "Removing all CT entries\n");
#endif
  for (unsigned i=0; i<tableSize; i++) {
    while (table[i]) {
      int curr = table[i];
#ifdef INTEGRATED_MEMMAN
      const int* entry = entries + curr;
#else
      const int* entry = (const int*) MMAN->getChunkAddress(curr);
#endif
      if (CHAINED) {
        table[i] = entry[0];
      } else {
        table[i] = 0;
      }
      discardAndRecycle(curr);
    } // while
  } // for i
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>
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
        curr = entries[curr];
#else
        const int* entry = (const int*) MMAN->getChunkAddress(curr);
        curr = entry[0];
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
          L = entries[L];
          s << "->" << L;
        }
        s << "\n";
      }
    }
  }

  if (0==entries) {
    s << "Entries: null\n";
  } else {
    s << "Entries: [" << long(entries[0]);
    for (int i=1; i<entriesSize; i++) {
      s << ", " << long(entries[i]);
    }
    s << "]\n";
  }
#else
  if (MMAN) MMAN->dumpInternal(s);
#endif

}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
int MEDDLY::ct_typebased<MONOLITHIC, CHAINED>::convertToList(bool removeStales)
{
  MEDDLY_DCASSERT(CHAINED);

  int list = 0;
  for (unsigned i=0; i<tableSize; i++) {
    while (table[i]) {
      int curr = table[i];
#ifdef INTEGRATED_MEMMAN
      int* entry = entries + curr;
#else
      int* entry = (int*) MMAN->getChunkAddress(curr);
#endif
      table[i] = entry[0];
      if (removeStales && isStale(entry)) {

#ifdef DEBUG_TABLE2LIST
          printf("\tstale ");
          FILE_output out(stdout);
          showEntry(out, curr);
          printf(" (handle %d table slot %d)\n", curr, i);
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
      printf(" (handle %d table slot %d)\n", curr, i);
#endif
      entry[0] = list;
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
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>::listToTable(int L)
{
  const int M = MONOLITHIC ? 1 : 0;
  MEDDLY_DCASSERT(CHAINED);
#ifdef DEBUG_LIST2TABLE
  printf("Recovering  list");
  // showChain(stdout, L);
#endif
  while (L) {
#ifdef INTEGRATED_MEMMAN
    int* entry = entries+L;
#else
    int* entry = (int*) MMAN->getChunkAddress(L);
#endif
    const int curr = L;
    L = entry[0];

    const entry_type* et = MONOLITHIC
      ?   getEntryType(entry[1])
      :   global_et;
    MEDDLY_DCASSERT(et);
    const unsigned reps = (et->isRepeating()) ? entry[ M+1 ] : 0;
    const unsigned hashlength = M + (et->isRepeating() ? 1 : 0) + ( et->getKeyBytes(reps) / sizeof(int) );

    const unsigned h = hash(entry + 1, hashlength);
    entry[0] = table[h];
    table[h] = curr;
#ifdef DEBUG_LIST2TABLE
    printf("\tsave  ");
    FILE_output out(stdout);
    showEntry(out, curr);
    printf(" (handle %d table slot %d hashlength %u)\n", curr, h, hashlength);
#endif
  }
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>::scanForStales()
{
  MEDDLY_DCASSERT(!CHAINED);
  for (unsigned i=0; i<tableSize; i++) {
    if (0==table[i]) continue;

#ifdef INTEGRATED_MEMMAN
    const int* entry = entries + table[i];
#else
    const int* entry = (const int*) MMAN->getChunkAddress(table[i]);
#endif

    if (isStale(entry)) {

#ifdef DEBUG_CT
      printf("Removing CT stale entry ");
      FILE_output out(stdout);
      showEntry(out, table[i]);
#ifdef DEBUG_CT_SLOTS
      printf(" in table slot %u\n", i);
#else
      printf("\n");
#endif
#endif  

      discardAndRecycle(table[i]);
      table[i] = 0;

    } // if isStale
  } // for i
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>::rehashTable(const int* oldT, unsigned oldS)
{
#ifdef DEBUG_REHASH
  printf("Rebuilding hash table\n");
#endif
  MEDDLY_DCASSERT(!CHAINED);
  const int M = MONOLITHIC ? 1 : 0;
  for (unsigned i=0; i<oldS; i++) {
    int curr = oldT[i];
    if (0==curr) continue;
#ifdef INTEGRATED_MEMMAN
    const int* entry = entries + curr;
#else
    const int* entry = (const int*) MMAN->getChunkAddress(curr);
#endif

    const entry_type* et = MONOLITHIC
      ?   getEntryType(entry[0])
      :   global_et;
    MEDDLY_DCASSERT(et);
    const unsigned reps = (et->isRepeating()) ? entry[ M ] : 0;
    const unsigned hashlength = M + (et->isRepeating() ? 1 : 0) + ( et->getKeyBytes(reps) / sizeof(int) );

    unsigned h = hash(entry, hashlength);
    setTable(h, curr);

#ifdef DEBUG_REHASH
    printf("\trehash  ");
    FILE_output out(stdout);
    showEntry(out, curr);
    printf(" (handle %d table slot %d hashlength %u)\n", curr, h, hashlength);
#endif
  }
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
MEDDLY::node_address MEDDLY::ct_typebased<MONOLITHIC, CHAINED>
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
    freeList[size] = entries[h];
#ifdef DEBUG_CTALLOC
    fprintf(stderr, "Re-used entry %ld size %d\n", h, size);
#endif
    mstats.incMemUsed( size * sizeof(int) );
    return h;
  }
  if (entriesSize + size > entriesAlloc) {
    // Expand by a factor of 1.5
    int neA = entriesAlloc + (entriesAlloc/2);
    int* ne = (int*) realloc(entries, neA * sizeof(int));
    if (0==ne) {
      fprintf(stderr,
          "Error in allocating array of size %lu at %s, line %d\n",
          neA * sizeof(int), __FILE__, __LINE__);
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    mstats.incMemAlloc( (entriesAlloc / 2) * sizeof(int) );
    entries = ne;
    entriesAlloc = neA;
  }
  MEDDLY_DCASSERT(entriesSize + size <= entriesAlloc);
  node_address h = entriesSize;
  entriesSize += size;
#ifdef DEBUG_CTALLOC
  fprintf(stderr, "New entry %ld size %d\n", h, size);
#endif
  mstats.incMemUsed( size * sizeof(int) );
  return h;

#else
  perf.numEntries++;
  size_t the_size = size;
  return MMAN->requestChunk(the_size);
#endif
}

// **********************************************************************

inline bool YES_stale() 
{
#ifdef DEBUG_ISSTALE
  printf("stale\n");
  fflush(stdout);
#endif
  return true;
}

inline bool NO_stale()
{
#ifdef DEBUG_ISSTALE
  printf("not stale\n");
  fflush(stdout);
#endif
  return false;
}

template <bool MONOLITHIC, bool CHAINED>
bool MEDDLY::ct_typebased<MONOLITHIC, CHAINED> 
::isStale(const int* entry) const
{
#ifdef DEBUG_ISSTALE
  printf("Checking entry for staleness: ");
  FILE_output out(stdout);
  showEntry(out, entry);
  printf("\n");
  fflush(stdout);
#endif
  const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

  //
  // Check entire entry for staleness
  //
  const entry_type* et = MONOLITHIC
      ?   getEntryType(entry[CHAINED?1:0])
      :   global_et;
  MEDDLY_DCASSERT(et);

#ifdef DEBUG_ISSTALE
  printf("\tChecking et marked?\n");
#endif
  if (et->isMarkedForDeletion()) return YES_stale();

  entry += SHIFT;
  const unsigned reps = (et->isRepeating()) ? *entry++ : 0;
  const unsigned klen = et->getKeySize(reps);


  const unsigned slots_for_type[] = { 1, 1, 1, 2, 1, 2, 2 };

  MEDDLY_DCASSERT(1 == slots_for_type[ERROR]);
  MEDDLY_DCASSERT(1 == slots_for_type[INTEGER]);
  MEDDLY_DCASSERT(sizeof(node_handle) / sizeof(int) == slots_for_type[NODE]);
  MEDDLY_DCASSERT(sizeof(long) / sizeof(int) == slots_for_type[LONG]);
  MEDDLY_DCASSERT(sizeof(float) / sizeof(int) == slots_for_type[FLOAT]);
  MEDDLY_DCASSERT(sizeof(double) / sizeof(int) == slots_for_type[DOUBLE]);
  MEDDLY_DCASSERT(sizeof(void*) / sizeof(int) == slots_for_type[POINTER]);

  //
  // Key portion
  //
  for (unsigned i=0; i<klen; i++) {
    typeID t;
    expert_forest* f;
    et->getKeyType(i, t, f);
    MEDDLY_CHECK_RANGE(0, t, 7);
    if (f) {
#ifdef DEBUG_ISSTALE
      printf("\tchecking key item %u\n", i);
#endif
      MEDDLY_DCASSERT(NODE == t);
      if (MEDDLY::forest::ACTIVE != f->getNodeStatus(*entry)) {
        return YES_stale();
      }
      entry++;
    } else {
#ifdef DEBUG_ISSTALE
      printf("\tskipping key item %u, %u slots\n", i, slots_for_type[t]);
#endif
      MEDDLY_DCASSERT(NODE != t);
      entry += slots_for_type[t];
    }

    /*
    switch (t) {
        case NODE:
                        MEDDLY_DCASSERT(f);
                        if (MEDDLY::forest::ACTIVE != f->getNodeStatus(*entry)) {
                          return true;
                        }
                        entry++;
                        continue;
        case INTEGER:
                        entry++;
                        continue;
        case LONG:
                        entry += sizeof(long) / sizeof(int);
                        continue;
        case FLOAT:
                        entry += sizeof(float) / sizeof(int);
                        continue;
        case DOUBLE:
                        entry += sizeof(double) / sizeof(int);
                        continue;
        case POINTER: 
                        entry += sizeof(void*) / sizeof(int);
                        continue;
        default:
                        MEDDLY_DCASSERT(0);
    }
    */
  } // for i

  // 
  // Result portion
  //
  for (unsigned i=0; i<et->getResultSize(); i++) {
    typeID t;
    expert_forest* f;
    et->getResultType(i, t, f);
    MEDDLY_CHECK_RANGE(0, t, 7);
    if (f) {
#ifdef DEBUG_ISSTALE
      printf("\tchecking result item %u\n", i);
#endif
      MEDDLY_DCASSERT(NODE == t);
      if (MEDDLY::forest::ACTIVE != f->getNodeStatus(*entry)) {
        return YES_stale();
      }
      entry++;
    } else {
#ifdef DEBUG_ISSTALE
      printf("\tskipping result item %u, %u slots\n", i, slots_for_type[t]);
#endif
      MEDDLY_DCASSERT(NODE != t);
      entry += slots_for_type[t];
    }
    /*
    switch (t) {
        case NODE:
                        MEDDLY_DCASSERT(f);
                        if (MEDDLY::forest::ACTIVE != f->getNodeStatus(*entry)) {
                          return true;
                        }
                        entry++;
                        continue;
        case INTEGER:
                        entry++;
                        continue;
        case LONG:
                        entry += sizeof(long) / sizeof(int);
                        continue;
        case FLOAT:
                        entry += sizeof(float) / sizeof(int);
                        continue;
        case DOUBLE:
                        entry += sizeof(double) / sizeof(int);
                        continue;
        case POINTER: 
                        entry += sizeof(void*) / sizeof(int);
                        continue;
        default:
                        MEDDLY_DCASSERT(0);
    }
    */
  } // for i

  return NO_stale();
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
bool MEDDLY::ct_typebased<MONOLITHIC, CHAINED> 
::isDead(const int* result, const entry_type* et) const
{
  MEDDLY_DCASSERT(et);
  MEDDLY_DCASSERT(result);

  const unsigned slots_for_type[] = { 1, 1, 1, 2, 1, 2, 2 };

  MEDDLY_DCASSERT(1 == slots_for_type[ERROR]);
  MEDDLY_DCASSERT(1 == slots_for_type[INTEGER]);
  MEDDLY_DCASSERT(sizeof(node_handle) / sizeof(int) == slots_for_type[NODE]);
  MEDDLY_DCASSERT(sizeof(long) / sizeof(int) == slots_for_type[LONG]);
  MEDDLY_DCASSERT(sizeof(float) / sizeof(int) == slots_for_type[FLOAT]);
  MEDDLY_DCASSERT(sizeof(double) / sizeof(int) == slots_for_type[DOUBLE]);
  MEDDLY_DCASSERT(sizeof(void*) / sizeof(int) == slots_for_type[POINTER]);

  //
  // Check result portion for dead nodes - cannot use result in that case
  //
  for (unsigned i=0; i<et->getResultSize(); i++) {
    typeID t;
    expert_forest* f;
    et->getResultType(i, t, f);
    MEDDLY_CHECK_RANGE(0, t, 7);
    if (f) {
#ifdef DEBUG_ISDEAD
      printf("\tchecking result item %u\n", i);
#endif
      MEDDLY_DCASSERT(NODE == t);
      if (MEDDLY::forest::DEAD == f->getNodeStatus(*result)) {
        return true;
      }
      result++;
    } else {
#ifdef DEBUG_ISDEAD
      printf("\tskipping result item %u, %u slots\n", i, slots_for_type[t]);
#endif
      MEDDLY_DCASSERT(NODE != t);
      result += slots_for_type[t];
    }
    /*
    switch (t) {
      case NODE:  
                  MEDDLY_DCASSERT(f);
                  if (MEDDLY::forest::DEAD == f->getNodeStatus(*result)) {
                    return true;
                  }
                  result++;
                  continue;
      case INTEGER:
                  result++;
                  continue;
      case FLOAT:
                  result += sizeof(float) / sizeof(int);
                  continue;
      case LONG:
                  result += sizeof(long) / sizeof(int);
                  continue;
      case DOUBLE:
                  result += sizeof(double) / sizeof(int);
                  continue;
      case POINTER:
                  result += sizeof(void*) / sizeof(int);
                  continue;
      default:
                  MEDDLY_DCASSERT(0);

    }
    */
  } // for i
  return false;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>
::discardAndRecycle(unsigned long h)
{
  const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

#ifdef INTEGRATED_MEMMAN
  int* entry = entries + h;
#else
  int* entry = (int*) MMAN->getChunkAddress(table[h]);
#endif

  unsigned slots;

  //
  // Discard.
  // For the portion of the entry that are nodes,
  // we need to notify the forest to decrement
  // the CT counter.
  //

  const entry_type* et = MONOLITHIC
    ?   getEntryType(entry[CHAINED ? 1 : 0])
    :   global_et;
  MEDDLY_DCASSERT(et);

  const int* ptr = entry + SHIFT;
  unsigned reps;
  if (et->isRepeating()) {
    reps = *((const unsigned*)(ptr));
    ptr++;
  } else {
    reps = 0;
  }

  const unsigned slots_for_type[] = { 1, 1, 1, 2, 1, 2, 2 };

  MEDDLY_DCASSERT(1 == slots_for_type[ERROR]);
  MEDDLY_DCASSERT(1 == slots_for_type[INTEGER]);
  MEDDLY_DCASSERT(sizeof(node_handle) / sizeof(int) == slots_for_type[NODE]);
  MEDDLY_DCASSERT(sizeof(long) / sizeof(int) == slots_for_type[LONG]);
  MEDDLY_DCASSERT(sizeof(float) / sizeof(int) == slots_for_type[FLOAT]);
  MEDDLY_DCASSERT(sizeof(double) / sizeof(int) == slots_for_type[DOUBLE]);
  MEDDLY_DCASSERT(sizeof(void*) / sizeof(int) == slots_for_type[POINTER]);

  //
  // Key portion
  //
  const unsigned stop = et->getKeySize(reps);
  for (unsigned i=0; i<stop; i++) {
    typeID t;
    expert_forest* f;
    et->getKeyType(i, t, f);
    MEDDLY_CHECK_RANGE(0, t, 7);
    switch (t) {
        case NODE:
                        MEDDLY_DCASSERT(f);
                        f->uncacheNode( *ptr );
                        ptr++;
                        continue;
        /*
        case INTEGER:
                        ptr++;
                        continue;
        case LONG:
                        ptr += sizeof(long) / sizeof(int);
                        continue;
        case FLOAT:
                        ptr += sizeof(float) / sizeof(int);
                        continue;
        case DOUBLE:
                        ptr += sizeof(double) / sizeof(int);
                        continue;
        */
        case POINTER: {
                        ct_object* P = *((ct_object**)(ptr));
                        delete P;
                        ptr += sizeof(void*) / sizeof(int);
                        continue;
                      }
        default:
                        MEDDLY_DCASSERT(ERROR != t);
                        ptr += slots_for_type[t];
    }
  } // for i

  //
  // Result portion
  //
  for (unsigned i=0; i<et->getResultSize(); i++) {
    typeID t;
    expert_forest* f;
    et->getResultType(i, t, f);
    MEDDLY_CHECK_RANGE(0, t, 7);
    switch (t) {
        case NODE:
                        MEDDLY_DCASSERT(f);
                        f->uncacheNode( *ptr );
                        ptr++;
                        continue;
        /*
        case INTEGER:
                        ptr++;
                        continue;
        case LONG:
                        ptr += sizeof(long) / sizeof(int);
                        continue;
        case FLOAT:
                        ptr += sizeof(float) / sizeof(int);
                        continue;
        case DOUBLE:
                        ptr += sizeof(double) / sizeof(int);
                        continue;
        */
        case POINTER: {
                        ct_object* P = *((ct_object**)(ptr));
                        delete P;
                        ptr += sizeof(void*) / sizeof(int);
                        continue;
                      }
        default:
                        MEDDLY_DCASSERT(ERROR != t);
                        ptr += slots_for_type[t];
    }
  } // for i
  slots = ptr - entry;

  //
  // Recycle
  //
#ifdef DEBUG_CTALLOC
  fprintf(stderr, "Recycling entry %ld size %u\n", long(h), slots);
#endif
#ifdef INTEGRATED_MEMMAN
  entries[h] = freeList[slots];
  freeList[slots] = h;
  mstats.decMemUsed( slots * sizeof(int) );
#else
  MMAN->recycleChunk(h, slots);
#endif
  perf.numEntries--;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>
::showEntry(output &s, const int* entry) const
{
  MEDDLY_DCASSERT(entry);

  const entry_type* et = MONOLITHIC
    ?   getEntryType(entry[CHAINED?1:0])
    :   global_et;
  MEDDLY_DCASSERT(et);

  const int* ptr = entry + (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);
  unsigned reps;
  if (et->isRepeating()) {
    reps = *((const unsigned*)(ptr));
    ptr++;
  } else {
    reps = 0;
  }
  s << "[" << et->getName() << "(";
  unsigned stop = et->getKeySize(reps);
  for (unsigned i=0; i<stop; i++) {
      entry_item item;
      if (i) s << ", ";
      switch (et->getKeyType(i)) {
        case NODE:
                        item.N = *ptr;
                        s.put(long(item.N));
                        ptr++;
                        break;
        case INTEGER:
                        item.I = *ptr;
                        s.put(long(item.I));
                        ptr++;
                        break;
        case LONG:
                        item.L = *((const long*)(ptr));
                        s.put(item.L);
                        ptr += sizeof(long) / sizeof(int);
                        break;
        case FLOAT:
                        item.F = *((const float*)(ptr));
                        s.put(item.F, 0, 0, 'e');
                        ptr += sizeof(float) / sizeof(int);
                        break;
        case DOUBLE:
                        item.D = *((const double*)(ptr));
                        s.put(item.D, 0, 0, 'e');
                        ptr += sizeof(double) / sizeof(int);
                        break;
        case POINTER:
                        item.P = *((void**)(ptr));
                        s.put_hex((unsigned long)item.P);
                        ptr += sizeof(void*) / sizeof(int);
                        break;
        default:
                        MEDDLY_DCASSERT(0);
      } // switch et->getKeyType()
  } // for i
  s << "): ";
  for (unsigned i=0; i<et->getResultSize(); i++) {
      entry_item item;
      if (i) s << ", ";
      switch (et->getResultType(i)) {
        case NODE:
                        item.N = *ptr;
                        s.put(long(item.N));
                        ptr++;
                        break;
        case INTEGER:
                        item.I = *ptr;
                        s.put(long(item.I));
                        ptr++;
                        break;
        case LONG:
                        item.L = *((const long*)(ptr));
                        s.put(item.L);
                        s << "(L)";
                        ptr += sizeof(long) / sizeof(int);
                        break;
        case FLOAT:
                        item.F = *((const float*)(ptr));
                        s.put(item.F, 0, 0, 'e');
                        ptr += sizeof(float) / sizeof(int);
                        break;
        case DOUBLE:
                        item.D = *((const double*)(ptr));
                        s.put(item.D, 0, 0, 'e');
                        ptr += sizeof(double) / sizeof(int);
                        break;
                            
        case POINTER:
                        item.P = *((void**)(ptr));
                        s.put_hex((unsigned long)item.P);
                        ptr += sizeof(void*) / sizeof(int);
                        break;
        default:
                        MEDDLY_DCASSERT(0);
      } // switch et->getResultType()
  } // for i
  s << "]";
}


// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_typebased<MONOLITHIC, CHAINED>
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


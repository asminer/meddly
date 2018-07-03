
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


#include "ct_classic.h"

#include <map>  // for operation_map
#include <limits.h>

// #define DEBUG_SLOW
// #define DEBUG_CT
// #define DEBUG_TABLE2LIST
// #define DEBUG_LIST2TABLE
// #define DEBUG_CTALLOC

// #define DEBUG_REMOVESTALES
// #define SUMMARY_STALES

#define INTEGRATED_MEMMAN



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                         ct_template  class                         *
// *                                                                    *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {
  template <bool MONOLITHIC, bool CHAINED>
  class ct_template : public compute_table {
    public:
      ct_template(const ct_initializer::settings &s, operation* op, unsigned slot);
      virtual ~ct_template();

      /**
          Find an entry.
          Used by find() and updateEntry().
            @param  key   Key to search for.
            @return Pointer to the entry, or null if not found.
      */
      int* findEntry(entry_key* key);

      // required functions

      virtual bool isOperationTable() const   { return !MONOLITHIC; }
      virtual entry_result& find(entry_key* key);
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

      inline unsigned hash(const int* k, int length) const {
        return raw_hash(k, length) % tableSize;
      }

      /*
          Use our own, built-in, specialized hash
          instead of hash_stream because it's faster.
      */

      static inline unsigned rot(unsigned x, int k) {
          return (((x)<<(k)) | ((x)>>(32-(k))));
      }
      static inline void mix(unsigned &a, unsigned &b, unsigned &c) {
          a -= c;  a ^= rot(c, 4);  c += b; 
          b -= a;  b ^= rot(a, 6);  a += c;
          c -= b;  c ^= rot(b, 8);  b += a;
          a -= c;  a ^= rot(c,16);  c += b;
          b -= a;  b ^= rot(a,19);  a += c;
          c -= b;  c ^= rot(b, 4);  b += a;
      }
      static inline void final(unsigned &a, unsigned &b, unsigned &c) {
          c ^= b; c -= rot(b,14);
          a ^= c; a -= rot(c,11);
          b ^= a; b -= rot(a,25);
          c ^= b; c -= rot(b,16);
          a ^= c; a -= rot(c,4); 
          b ^= a; b -= rot(a,14);
          c ^= b; c -= rot(b,24);
      }

      static inline unsigned raw_hash(const int* k, int length) {
        unsigned a, b, c;
        a = b = c = 0xdeadbeef;

        // handle most of the key
        while (length > 3)
        {
          a += *k++;
          b += *k++;
          c += *k++;
          mix(a,b,c);
          length -= 3;
        }

        // handle the last 3 uint32_t's
        switch(length)
        { 
          // all the case statements fall through
          case 3: c += k[2];
          case 2: b += k[1];
          case 1: a += k[0];
                  final(a,b,c);
          case 0: // nothing left to add
                  break;
        }

        return c;
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
        printf(" in slot %u\n", h);
#endif
        collisions++;    
        //
        // Remove table[h]
        //
#ifdef INTEGRATED_MEMMAN
        int* entry = entries + table[h];
#else
        int* entry = (int*) MMAN->getChunkAddress(table[h]);
#endif
        operation* op = MONOLITHIC
          ?   operation::getOpWithIndex(entry[0])
          :   global_op;
        MEDDLY_DCASSERT(op);
        op->discardEntry(entry + (MONOLITHIC?1:0) );  
        recycleEntry(table[h], op->getCacheEntryLength() + (MONOLITHIC?1:0));
        //
        table[h] = curr;
      }

      /// Recycle an entry
      inline void recycleEntry(node_address h, unsigned size) {
#ifdef DEBUG_CTALLOC
        fprintf(stderr, "Recycling entry %ld size %u\n", long(h), size);
#endif
#ifdef INTEGRATED_MEMMAN
        entries[h] = freeList[size];
        freeList[size] = h;
        mstats.decMemUsed( size * sizeof(int) );
#else
        MMAN->recycleChunk(h, size);
#endif
        perf.numEntries--;
      }

      /// Check equality
      static inline bool equal_sw(const int* a, const int* b, int N) {
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


      inline void showEntry(output &s, int h) const
      {
#ifdef INTEGRATED_MEMMAN
        const int* entry = entries+h;
#else
        const int* entry = (const int*) MMAN->getChunkAddress(h);
#endif
        operation* op = MONOLITHIC
            ?   operation::getOpWithIndex(entry[CHAINED ? 1 : 0])
            :   global_op;
        MEDDLY_DCASSERT(op);
        op->showEntry(s, entry + (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0));
      }




    private:
      /// Global operation.  Ignored when MONOLITHIC is true.
      operation* global_op; 

      /// Global slot number (with operation).  Ignored when MONOLITHIC is true.
      unsigned global_slot;


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
  }; // class ct_template
} // namespace


// **********************************************************************
// *                                                                    *
// *                        ct_template  methods                        *
// *                                                                    *
// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
MEDDLY::ct_template<MONOLITHIC, CHAINED>::ct_template(
  const ct_initializer::settings &s, operation* op, unsigned slot)
: compute_table(s)
{
  if (MONOLITHIC) {
    MEDDLY_DCASSERT(0==op);
    MEDDLY_DCASSERT(0==slot);
  } else {
    MEDDLY_DCASSERT(op);
  }
  global_op = op;
  global_slot = slot;

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
MEDDLY::ct_template<MONOLITHIC, CHAINED>::~ct_template()
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
inline int* MEDDLY::ct_template<MONOLITHIC, CHAINED>
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

  const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

  /*
#ifdef DEBUG_CT
  printf("Searching for CT entry ");
  FILE_output out(stdout);
  key->show(out, key->dataLength());
  printf("\n");
#endif
  */

  int chain;

  int* answer = 0;
  int* preventry = 0;
  unsigned hcurr = key->getHash() % tableSize;
  int curr = table[hcurr];
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

    //
    // Check for match
    //
    if (equal_sw(entry + (CHAINED?1:0), key->rawData(MONOLITHIC), key->dataLength(MONOLITHIC))) {
      if (key->getOp()->shouldStaleCacheHitsBeDiscarded()) {
#ifndef USE_NODE_STATUS
        const bool stale = key->getOp()->isEntryStale(entry+SHIFT);
#else
        const bool stale = (MEDDLY::forest::node_status::DEAD ==
             key->getOp()->getEntryStatus(entry+SHIFT) );
#endif
        if (stale) {
          //
          // The match is stale.
          // Delete the entry.
          //
#ifdef DEBUG_CT
          printf("Removing stale CT hit ");
          FILE_output out(stdout);
          key->getOp()->showEntry(out, entry + SHIFT);
          printf(" in slot %u\n", hcurr);
#endif
          key->getOp()->discardEntry(entry+SHIFT);
          if (CHAINED) {
            if (preventry) {
              preventry[0] = entry[0];
            } else {
              table[hcurr] = entry[0];
            }
          } else {
            table[hcurr] = 0;
          }
          recycleEntry(curr, key->getOp()->getCacheEntryLength() + SHIFT);

          // Since there can NEVER be more than one match
          // in the table, we're done!
          break;
        }
      } // shouldStaleCacheHitsBeDiscarded

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
      key->getOp()->showEntry(out, entry + SHIFT);
      printf(" in slot %u\n", hcurr);
#endif
      answer = entry;
      break;
    } // equal_sw


    // 
    // Not a cache hit.
    //
    if (checkStalesOnFind) {
      operation* op = MONOLITHIC
          ?   operation::getOpWithIndex(entry[CHAINED ? 1 : 0])
          :   global_op;
      MEDDLY_DCASSERT(op);

#ifndef USE_NODE_STATUS
      const bool stale = key->getOp()->isEntryStale(entry+SHIFT);
#else
      const bool stale = (MEDDLY::forest::node_status::ACTIVE !=
             key->getOp()->getEntryStatus(entry+SHIFT) );
#endif
      if (stale) {
        //
        // Delete the entry
        //
#ifdef DEBUG_CT
        printf("Removing stale CT entry ");
        FILE_output out(stdout);
        op->showEntry(out, entry + SHIFT);
        printf(" in slot %u\n", hcurr);
#endif
        op->discardEntry(entry+SHIFT);
        if (CHAINED) {
          if (preventry) {
            preventry[0] = entry[0];
          } else {
            table[hcurr] = entry[0];
          }
        } else {
          table[hcurr] = 0;
        }
        recycleEntry(curr, op->getCacheEntryLength() + SHIFT);
      }

    } // if checkStalesOnFind


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
MEDDLY::compute_table::entry_result& MEDDLY::ct_template<MONOLITHIC, CHAINED>
::find(entry_key *key)
{
  static entry_result ANS;

  setHash(key, raw_hash(key->rawData(MONOLITHIC), key->dataLength(MONOLITHIC)));
  int* entry = findEntry(key);
  perf.pings++;

  if (entry) {
    perf.hits++;
    ANS.setResult(
      entry+(CHAINED?1:0)+key->dataLength(MONOLITHIC), 
      key->getOp()->getAnsLength()
    );
  } else {
    ANS.setInvalid();
  }
  return ANS;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_template<MONOLITHIC, CHAINED>::addEntry(entry_key* key, const entry_result &res)
{
  MEDDLY_DCASSERT(key);
  if (!MONOLITHIC) {
    if (key->getOp() != global_op)
      throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
  }

  unsigned h = key->getHash() % tableSize;

#ifdef DEBUG_CT
  printf("Adding CT entry ");
  FILE_output out(stdout);
  showEntry(out, currEntry.readHandle());
  printf(" in slot %u\n", h);
#endif

  operation* op = key->getOp();
  MEDDLY_DCASSERT(op);

  //
  // Allocate an entry
  //

  node_address curr = newEntry(
    op->getCacheEntryLength() + 
    (CHAINED ? 1 : 0) + 
    (MONOLITHIC ? 1 : 0)
  );

#ifdef INTEGRATED_MEMMAN
  int* entry = entries + curr;
#else
  int* entry = (int*) MMAN->getChunkAddress(curr);
#endif

  //
  // Copy into the entry
  //
  int* key_portion = entry + (CHAINED ? 1 : 0);
  memcpy(key_portion, key->rawData(MONOLITHIC), key->dataLength(MONOLITHIC)*sizeof(node_handle));
  int* res_portion = key_portion + key->dataLength(MONOLITHIC);
  memcpy(res_portion, res.rawData(), res.dataLength()*sizeof(node_handle));

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
void MEDDLY::ct_template<MONOLITHIC, CHAINED>::updateEntry(entry_key* key, const entry_result &res)
{
  int* entry = findEntry(key);
  if (entry) {
    int* key_portion = entry + (CHAINED ? 1 : 0);
    int* res_portion = key_portion + key->dataLength(MONOLITHIC);
    memcpy(res_portion, res.rawData(), res.dataLength()*sizeof(node_handle));
  } else {
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_template<MONOLITHIC, CHAINED>::removeStales()
{
#ifdef DEBUG_SLOW
  fprintf(stdout, "Removing stales in CT (size %d, entries %u)\n", 
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
    
#ifdef DEBUG_SLOW
  fprintf(stdout, "Done removing CT stales (size %d, entries %u)\n", 
    tableSize, perf.numEntries
  );
#endif
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_template<MONOLITHIC, CHAINED>::removeAll()
{
  const int SHIFT = (MONOLITHIC?1:0) + (CHAINED?1:0);

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

      operation* currop = MONOLITHIC
          ?   operation::getOpWithIndex(entry[CHAINED ? 1 : 0])
          :   global_op;
      MEDDLY_DCASSERT(currop);

      currop->discardEntry(entry + SHIFT);
      recycleEntry(curr, SHIFT+currop->getCacheEntryLength());
    } // while
  } // for i
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_template<MONOLITHIC, CHAINED>
::show(output &s, int verbLevel)
{
  if (verbLevel < 1) return;

  if (MONOLITHIC) {
    s << "Monolithic compute table\n";
  } else {
    s << "Compute table for " << global_op->getName() << " (index " 
      << long(global_op->getIndex()) << ")\n";
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
  if (0==entries) {
    s << "Entries: null\n";
  } else {
    s << "Entries: [" << long(entries[0]);
    for (int i=1; i<entriesSize; i++) {
      s << ", " << long(entries[i]);
    }
    s << "]\n";
  }
  if (0==freeList) {
    s << "Free: null\n";
  } else {
    s << "Free: [" << long(freeList[0]);
    for (int i=1; i<=maxEntrySize; i++) {
      s << ", " << long(freeList[i]);
    }
    s << "]\n";
  }
#else
  if (MMAN) MMAN->dumpInternal(s);
#endif

}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
int MEDDLY::ct_template<MONOLITHIC, CHAINED>::convertToList(bool removeStales)
{
  MEDDLY_DCASSERT(CHAINED);
  const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

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
      if (removeStales) {
        operation* currop = MONOLITHIC
          ?   operation::getOpWithIndex(entry[1])
          :   global_op;
        MEDDLY_DCASSERT(currop);
        //
        // Check for stale
        //
#ifndef USE_NODE_STATUS
        if (currop->isEntryStale(entry+SHIFT)) {
#else
        if (currop->getEntryStatus(entry+SHIFT)) {
#endif
#ifdef DEBUG_TABLE2LIST
          printf("\tstale ");
          FILE_output out(stdout);
          currop->showEntry(out, entry+SHIFT);
          printf(" (handle %d slot %d)\n", curr, i);
#endif
          currop->discardEntry(entry+SHIFT);
          recycleEntry(curr, SHIFT+currop->getCacheEntryLength());
          continue;
        }
      } // if removeStales
      //
      // Not stale, move to list
      //
#ifdef DEBUG_TABLE2LIST
      printf("\tkeep  ");
      FILE_output out(stdout);
      currop->showEntry(out, entry+SHIFT);
      printf(" (handle %d slot %d)\n", curr, i);
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
void MEDDLY::ct_template<MONOLITHIC, CHAINED>::listToTable(int L)
{
  MEDDLY_DCASSERT(CHAINED);
  const int M = MONOLITHIC ? 1 : 0;
#ifdef DEBUG_LIST2TABLE
  printf("Recovering  list: ");
  showChain(stdout, L);
#endif
  while (L) {
#ifdef INTEGRATED_MEMMAN
    int* entry = entries+L;
#else
    int* entry = (int*) MMAN->getChunkAddress(L);
#endif
    const int curr = L;
    L = entry[0];
    operation* currop = MONOLITHIC
      ?   operation::getOpWithIndex(entry[1])
      :   global_op;
    MEDDLY_DCASSERT(currop);
    const int hashlength = M+currop->getKeyLength();
    const unsigned h = hash(entry + 1, hashlength);
    entry[0] = table[h];
    table[h] = curr;
#ifdef DEBUG_LIST2TABLE
    printf("\tsave  ");
    FILE_output out(stdout);
    currop->showEntry(out, entry+M+1);
    printf(" (handle %d slot %d)\n", curr, h);
#endif
  }
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_template<MONOLITHIC, CHAINED>::scanForStales()
{
  MEDDLY_DCASSERT(!CHAINED);
  const int M = MONOLITHIC ? 1 : 0;
  for (unsigned i=0; i<tableSize; i++) {
    if (0==table[i]) continue;
#ifdef INTEGRATED_MEMMAN
    const int* entry = entries + table[i];
#else
    const int* entry = (const int*) MMAN->getChunkAddress(table[i]);
#endif
    operation* currop = MONOLITHIC
      ?   operation::getOpWithIndex(entry[0])
      :   global_op;
    MEDDLY_DCASSERT(currop);

#ifndef USE_NODE_STATUS
    if (!currop->isEntryStale(entry+M)) {
      continue;
    }
#else
    if (MEDDLY::forest::node_status::ACTIVE == currop->getEntryStatus(entry+M) {
      continue;
    }
#endif

#ifdef DEBUG_CT
    printf("Removing CT stale entry ");
    FILE_output out(stdout);
    currop->showEntry(out, entry + M );
    printf(" in slot %u\n", i);
#endif  
    currop->discardEntry(entry + M );
    int length = currop->getCacheEntryLength();
    recycleEntry(table[i], length + M );
    table[i] = 0;
  } // for i
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_template<MONOLITHIC, CHAINED>::rehashTable(const int* oldT, unsigned oldS)
{
  MEDDLY_DCASSERT(!CHAINED);
  for (unsigned i=0; i<oldS; i++) {
    int curr = oldT[i];
    if (0==curr) continue;
#ifdef INTEGRATED_MEMMAN
    const int* entry = entries + curr;
#else
    const int* entry = (const int*) MMAN->getChunkAddress(curr);
#endif
    operation* currop = MONOLITHIC
      ?   operation::getOpWithIndex(entry[0])
      :   global_op;
    MEDDLY_DCASSERT(currop);
    int hashlength = currop->getKeyLength() + (MONOLITHIC ? 1 : 0);
    unsigned h = hash(entry, hashlength);
    setTable(h, curr);
  }
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
MEDDLY::node_address MEDDLY::ct_template<MONOLITHIC, CHAINED>
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
    fprintf(stderr, "Re-used entry %d size %d\n", h, size);
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
  fprintf(stderr, "New entry %d size %d\n", h, size);
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
// *                                                                    *
// *                  monolithic_chained_style methods                  *
// *                                                                    *
// **********************************************************************

MEDDLY::monolithic_chained_style::monolithic_chained_style() 
{ 
}

MEDDLY::compute_table* 
MEDDLY::monolithic_chained_style::create(const ct_initializer::settings &s) const 
{
  return new ct_template<true, true>(s, 0, 0);
}

bool MEDDLY::monolithic_chained_style::usesMonolithic() const 
{
  return true;
}

// **********************************************************************
// *                                                                    *
// *                 monolithic_unchained_style methods                 *
// *                                                                    *
// **********************************************************************


MEDDLY::monolithic_unchained_style::monolithic_unchained_style() 
{ 
}

MEDDLY::compute_table* 
MEDDLY::monolithic_unchained_style::create(const ct_initializer::settings &s) const 
{
  return new ct_template<true, false>(s, 0, 0);
}

bool MEDDLY::monolithic_unchained_style::usesMonolithic() const 
{
  return true;
}

// **********************************************************************
// *                                                                    *
// *                  operation_chained_style  methods                  *
// *                                                                    *
// **********************************************************************

MEDDLY::operation_chained_style::operation_chained_style() 
{ 
}

MEDDLY::compute_table* 
MEDDLY::operation_chained_style::create(const ct_initializer::settings &s, operation* op, unsigned slot) const 
{
  return new ct_template<false, true>(s, op, slot);
}

bool MEDDLY::operation_chained_style::usesMonolithic() const 
{
  return false;
}


// **********************************************************************
// *                                                                    *
// *                 operation_unchained_style  methods                 *
// *                                                                    *
// **********************************************************************


MEDDLY::operation_unchained_style::operation_unchained_style() 
{ 
}

MEDDLY::compute_table* 
MEDDLY::operation_unchained_style::create(const ct_initializer::settings &s, operation* op, unsigned slot) const 
{
  return new ct_template<false, false>(s, op, slot);
}

bool MEDDLY::operation_unchained_style::usesMonolithic() const 
{
  return false;
}



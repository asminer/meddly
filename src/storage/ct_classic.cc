
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

// #include "../hash_stream.h"

// #include <map>  // for operation_map
// #include <limits.h>

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
            @return Pointer to the result portion of the entry, or null if not found.
      */
      int* findEntry(entry_key* key);

      // required functions

      virtual bool isOperationTable() const   { return !MONOLITHIC; }
#ifdef OLD_OP_CT
      virtual entry_result& find(entry_key* key);
#else
      virtual void find(entry_key* key, entry_result &res);
#endif
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
        printf(" in slot %u\n", h);
#endif
        collisions++;    

        discardAndRecycle(table[h]);
        table[h] = curr;
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

#ifndef OLD_OP_CT
      /**
          Check if the key portion of an entry equals key and should be discarded.
            @param  entry   Complete entry in CT to check.
            @param  key     Key to compare against. 
            @param  discard On output: should this entry be discarded

            @return Result portion of the entry, if they key portion matches;
                    0 if the entry does not match the key.
      */
      int* checkEqualityAndStatus(const int* entry, const entry_key* key, bool &discard);
#endif

#ifdef OLD_OP_CT
      inline bool isStale(unsigned long h) const 
      {
          const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

#ifdef INTEGRATED_MEMMAN
          const int* entry = entries + h;
#else
          const int* entry = (const int*) MMAN->getChunkAddress(h);
#endif
          operation* currop = MONOLITHIC
            ?   operation::getOpWithIndex(entry[CHAINED?1:0])
            :   global_op;
          MEDDLY_DCASSERT(currop);

#ifndef USE_NODE_STATUS
          return currop->isEntryStale(entry+SHIFT);
#else
          return MEDDLY::forest::node_status::ACTIVE != currop->getEntryStatus(entry+SHIFT);
#endif
      }

#else // OLD_OP_CT

      /**
          Check if an entry is stale.
            @param  h   Handle of entry to check.
            @return true, if the entry should be discarded.
      */
      bool isStale(unsigned long h) const;

#endif // OLD_OP_CT

      /**
          Copy a result into an entry.
      */
#ifdef OLD_OP_CT
      inline void setResult(int* respart, const entry_result &res) {
        memcpy(respart, res.rawData(), res.dataLength()*sizeof(node_handle));
      }
#else
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
#endif

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

#ifndef OLD_OP_CT
      /**
          Determine if an entry is unrecoverable.
          This means the entry is unusable because the
          result node(s) are DEAD.
            @param  h   Handle of the entry.
      */
      // bool entryUnrecoverable(unsigned long h) const;
#endif

      /**
          Discard an entry and recycle the memory it used.
            @param  h   Handle of the entry.
      */
      void discardAndRecycle(unsigned long h);

      /**
          Display an entry.
          Used for debugging, and by method(s) to display the entire CT.
            @param  s   Stream to write to.
            @param  h   Handle of the entry.
      */
      void showEntry(output &s, unsigned long h) const;


    private:
#ifdef OLD_OP_CT
      /// Global operation.  Ignored when MONOLITHIC is true.
      operation* global_op; 
#else
      /// Global entry type.  Ignored when MONOLITHIC is true.
      const entry_type* global_et;
#endif

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
#ifdef OLD_OP_CT
  global_op = op;
#else
  global_et = getEntryType(op, slot);
#endif

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

#ifdef OLD_OP_CT
  const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);
#endif

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


#ifdef OLD_OP_CT  ////////////////////////////////////////////////////////////// NOT OLD_OP_CT

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
#endif // USE_NODE_STATUS
        if (stale) {
          //
          // The match cannot be used.
          // Delete the entry.
          //
#ifdef DEBUG_CT
          printf("Removing stale CT hit ");
          FILE_output out(stdout);
          showEntry(out, hcurr);
          printf(" in slot %u\n", hcurr);
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
      showEntry(out, hcurr);
      printf(" in slot %u\n", hcurr);
#endif
      answer = entry + (CHAINED ? 1 : 0) + key->dataLength(MONOLITHIC);
      break;
    } // equal_sw


    // 
    // Not a cache hit.
    //
    if (checkStalesOnFind) {

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
        showEntry(out, hcurr);
        printf(" in slot %u\n", hcurr);
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
      }
    }

#else  ////////////////////////////////////////////////////////////// NOT OLD_OP_CT

    bool discard;
    answer = checkEqualityAndStatus(entry, key, discard);

    if (discard) {
        //
        // Delete the entry.
        //
#ifdef DEBUG_CT
        if (equal)  printf("Removing stale CT hit   ");
        else        printf("Removing stale CT entry ");
        FILE_output out(stdout);
        showEntry(out, hcurr);
        printf(" in slot %u\n", hcurr);
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
        showEntry(out, hcurr);
        printf(" in slot %u\n", hcurr);
#endif
        break;
    } // if equal

#endif ////////////////////////////////////////////////////////////// OLD_OP_CT

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

#ifdef OLD_OP_CT
template <bool MONOLITHIC, bool CHAINED>
MEDDLY::compute_table::entry_result& 
MEDDLY::ct_template<MONOLITHIC, CHAINED>
::find(entry_key *key)
#else
template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_template<MONOLITHIC, CHAINED>
::find(entry_key *key, entry_result& res)
#endif
{
#ifdef OLD_OP_CT
  static entry_result res;
  setHash(key, raw_hash(key->rawData(MONOLITHIC), key->dataLength(MONOLITHIC)));
#else
  //
  // Allocate temporary space for key preprocessing.
  //
  int* temp_entry = (int*) key->allocTempData( key->dataBytes(MONOLITHIC) );

  //
  // Copy the operation index if we're monolithic
  //
  const entry_item* data = key->rawData(MONOLITHIC);
  int tptr = 0;
  if (MONOLITHIC) {
    temp_entry[0] = data[0].U;
    data++;
    tptr++;
  }

  //
  // Copy the key into temp_entry
  //
  const entry_type* et = key->getET();
  MEDDLY_DCASSERT(et);
  unsigned datalen = key->dataLength(false);
  for (unsigned i=0; i<datalen; i++) {
    typeID t = et->getKeyType(i);
    switch (t) {
      case NODE:
                  temp_entry[tptr++] = data[i].N;
                  continue;
      case INTEGER:
                  temp_entry[tptr++] = data[i].I;
                  continue;
      case LONG: 
                  memcpy(temp_entry+tptr, &(data[i].L), sizeof(long));
                  tptr += sizeof(long) / sizeof(int);
                  continue;
      case FLOAT:
                  temp_entry[tptr++] = data[i].I;   // Hack!
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
  setHash(key, raw_hash(temp_entry, key->dataBytes(MONOLITHIC) / sizeof(int)));

#endif


  int* entry_result = findEntry(key);
  perf.pings++;

  if (entry_result) {
    perf.hits++;
#ifdef OLD_OP_CT
    res.setResult(entry_result, key->getOp()->getAnsLength());
#else
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
#endif
  } else {
    res.setInvalid();
  }
#ifdef OLD_OP_CT
  return res;
#endif
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_template<MONOLITHIC, CHAINED>::addEntry(entry_key* key, const entry_result &res)
{
  MEDDLY_DCASSERT(key);
  if (!MONOLITHIC) {
#ifdef OLD_OP_CT
    if (key->getOp() != global_op)
      throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
#else
    if (key->getET() != global_et)
      throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
#endif
  }

  unsigned h = key->getHash() % tableSize;

#ifdef DEBUG_CT
  printf("Adding CT entry ");
  FILE_output out(stdout);
  showEntry(out, currEntry.readHandle());
  printf(" in slot %u\n", h);
#endif

#ifdef OLD_OP_CT
  operation* op = key->getOp();
  MEDDLY_DCASSERT(op);
#else
  const entry_type* et = key->getET();
  MEDDLY_DCASSERT(et);
#endif

  //
  // Allocate an entry
  //

#ifdef OLD_OP_CT
  node_address curr = newEntry(
    op->getCacheEntryLength() + 
    (CHAINED ? 1 : 0) + 
    (MONOLITHIC ? 1 : 0)
  );
#else
  node_address curr = newEntry(
    key->dataBytes(MONOLITHIC) / sizeof(int) +
    res.dataBytes() / sizeof(int) +
    (CHAINED ? 1 : 0)
  );
#endif

#ifdef INTEGRATED_MEMMAN
  int* entry = entries + curr;
#else
  int* entry = (int*) MMAN->getChunkAddress(curr);
#endif

  //
  // Copy into the entry
  //
  int* key_portion = entry + (CHAINED ? 1 : 0);
#ifdef OLD_OP_CT
  memcpy(key_portion, key->rawData(MONOLITHIC), key->dataLength(MONOLITHIC)*sizeof(node_handle));
  int* res_portion = key_portion + key->dataLength(MONOLITHIC);
  setResult(res_portion, res);
#else
  memcpy(key_portion, key->readTempData(), key->dataBytes(MONOLITHIC));
  int* res_portion = key_portion + key->dataBytes(MONOLITHIC) / sizeof(int);
  setResult(res_portion, res, et);
#endif

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
#ifndef OLD_OP_CT
  MEDDLY_DCASSERT(key->getET()->isResultUpdatable());
#endif
  int* entry_result = findEntry(key);
  if (!entry_result) {
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }

#ifdef OLD_OP_CT
  setResult(entry_result, res);
#else
  setResult(entry_result, res, key->getET());
#endif
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
void MEDDLY::ct_template<MONOLITHIC, CHAINED>
::show(output &s, int verbLevel)
{
  if (verbLevel < 1) return;

  if (MONOLITHIC) {
    s << "Monolithic compute table\n";
  } else {
#ifdef OLD_OP_CT
    s << "Compute table for " << global_op->getName() << " (index " 
      << long(global_op->getIndex()) << ")\n";
#else
    s << "Compute table for " << global_et->getName() << " (index " 
      << long(global_et->getID()) << ")\n";
#endif
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
  // const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

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
      if (removeStales && isStale(curr)) {

        /*
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
          */

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
    showEntry(out, curr);
    printf(" (handle %d table slot %d)\n", curr, h);
#endif
  }
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_template<MONOLITHIC, CHAINED>::scanForStales()
{
  MEDDLY_DCASSERT(!CHAINED);
  for (unsigned i=0; i<tableSize; i++) {
    if (0==table[i]) continue;

    if (isStale(table[i])) {

#ifdef DEBUG_CT
      printf("Removing CT stale entry ");
      FILE_output out(stdout);
      showEntry(out, table[i]);
      printf(" in table slot %u\n", i);
#endif  

      discardAndRecycle(table[i]);
      table[i] = 0;

    } // if isStale
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

#ifndef OLD_OP_CT

template <bool MONOLITHIC, bool CHAINED>
int* MEDDLY::ct_template<MONOLITHIC, CHAINED>
::checkEqualityAndStatus(const int* entry, const entry_key* k, bool &equal, bool &discard)
{

  equal = equal_sw(entry, (const int*) k->readTempData(), k->dataBytes(MONOLITHIC) / sizeof(int) );

  //
  // Now, scan entry and check for staleness
  //

  // TBD!
  MEDDLY_DCASSERT(0);
  return 0;
}

#endif

// **********************************************************************

#ifndef OLD_OP_CT

template <bool MONOLITHIC, bool CHAINED>
bool MEDDLY::ct_template<MONOLITHIC, CHAINED> 
::isStale(unsigned long h) const
{
}

#endif

// **********************************************************************

/*
#ifndef OLD_OP_CT

template <bool MONOLITHIC, bool CHAINED>
bool MEDDLY::ct_template<MONOLITHIC, CHAINED>
::entryUnrecoverable(unsigned long h) const
{
  const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

#ifdef INTEGRATED_MEMMAN
  int* entry = entries + h;
#else
  int* entry = (int*) MMAN->getChunkAddress(table[h]);
#endif

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


  MEDDLY_DCASSERT(0);
  return true;
}

#endif
*/

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_template<MONOLITHIC, CHAINED>
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

#ifdef OLD_OP_CT

  operation* op = MONOLITHIC
    ?   operation::getOpWithIndex(entry[CHAINED ? 1 : 0])
    :   global_op;
  MEDDLY_DCASSERT(op);
  op->discardEntry(entry + SHIFT);    // Old way - OP handles decrementing counters
  slots = op->getCacheEntryLength() + SHIFT;

#else

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

  //
  // Key portion
  //
  slots = et->getKeySize(reps);
  for (unsigned i=0; i<slots; i++) {
    typeID t;
    forest* f;
    et->getKeyType(i, t, f);
    if (f) {
      MEDDLY_DCASSERT(NODE == t);
      f->uncacheNode( *ptr );
      ptr++;
      continue;
    }
    switch (t) {
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
        default:
                        MEDDLY_DCASSERT(0);
    }
  } // for i

  //
  // Result portion
  //
  slots += et->getResultSize();
  for (unsigned i=0; i<et->getResultSize(); i++) {
    typeID t;
    forest* f;
    et->getResultType(i, t, f);
    if (f) {
      MEDDLY_DCASSERT(NODE == t);
      f->uncacheNode( *ptr );
      ptr++;
      continue;
    }
    switch (t) {
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
        default:
                        MEDDLY_DCASSERT(0);
    }
  } // for i

#endif

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
void MEDDLY::ct_template<MONOLITHIC, CHAINED>
::showEntry(output &s, unsigned long h) const
{
#ifdef INTEGRATED_MEMMAN
  const int* entry = entries+h;
#else
  const int* entry = (const int*) MMAN->getChunkAddress(h);
#endif

#ifdef OLD_OP_CT
  operation* op = MONOLITHIC
      ?   operation::getOpWithIndex(entry[CHAINED ? 1 : 0])
      :   global_op;
  MEDDLY_DCASSERT(op);
  op->showEntry(s, entry + (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0));

#else

  const entry_type* et = MONOLITHIC
    ?   getEntryType(entry[0])
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
                        s.put(item.F);
                        ptr += sizeof(float) / sizeof(int);
                        break;
        case DOUBLE:
                        item.D = *((const double*)(ptr));
                        s.put(item.D);
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
                        ptr += sizeof(long) / sizeof(int);
                        break;
        case FLOAT:
                        item.F = *((const float*)(ptr));
                        s.put(item.F);
                        ptr += sizeof(float) / sizeof(int);
                        break;
        case DOUBLE:
                        item.D = *((const double*)(ptr));
                        s.put(item.D);
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
#endif  // ifdef OLD_OP_CT
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



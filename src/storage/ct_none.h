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


#ifndef MEDDLY_CT_NONE_H
#define MEDDLY_CT_NONE_H

#include "../node_storage.h"
#include "../compute_table.h"
#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../ct_vector.h"
#include "../ct_initializer.h"

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
      ct_none(const ct_settings &s, operation* op, unsigned slot);
      virtual ~ct_none();

    protected:
      /**
          Find an entry.
          Used by find() and updateEntry().
            @param  key   Key to search for.
            @return Pointer to the result portion of the entry, or null if not found.
      */
      ct_entry_item* findEntry(ct_entry_key* key);

      /**
          Find an entry.
          Used by find() and updateEntry().
            @param  key   Key to search for.
            @return Pointer to the result portion of the entry, or null if not found.
      */
      ct_entry_item* findEntry(const ct_entry_type &ET, const ct_vector &key);

    public:

      // required functions

      virtual void find(ct_entry_key* key, ct_entry_result &res);
      virtual void addEntry(ct_entry_key* key, const ct_entry_result& res);
      virtual void updateEntry(ct_entry_key* key, const ct_entry_result& res);

      // new functions

      virtual bool find(const ct_entry_type &ET, ct_vector &k, ct_vector &r);
      virtual void addEntry(const ct_entry_type &ET, ct_vector &key,
              const ct_vector &res);


      virtual void removeStales();
      virtual void removeAll();
      virtual void show(output &s, int verbLevel = 0);
      virtual void countNodeEntries(const forest* f, size_t* counts) const;

    private:  // helper methods

      inline void scanForStales(unsigned long i) {
          //
          // Check entry i (only) for staleness
          //
          MEDDLY_DCASSERT(!CHAINED);

          if (0==table[i]) return;

#ifdef INTEGRATED_MEMMAN
          const ct_entry_item* entry = entries + table[i];
#else
          const ct_entry_item* entry = (const ct_entry_item*) MMAN->getChunkAddress(table[i]);
#endif
          if (!isStale(entry)) return;

#ifdef DEBUG_CT_SCAN
          printf("SCAN removing CT stale entry ");
          FILE_output out(stdout);
          showEntry(out, table[i]);
          printf(" in table slot %u\n", i);
#endif
          discardAndRecycle(table[i]);
          table[i] = 0;
      }

      void scanListForStales(unsigned long i) {
          //
          // Check list at slot i for staleness
          //
          MEDDLY_DCASSERT(CHAINED);

          unsigned long curr = table[i];
          table[i] = 0;
          ct_entry_item* preventry = 0;

          while (curr) {
#ifdef INTEGRATED_MEMMAN
            ct_entry_item* entry = entries + curr;
#else
            ct_entry_item* entry = (ct_entry_item*) MMAN->getChunkAddress(curr);
#endif
            unsigned long next = entry[0].UL;

            if (isStale(entry)) {

#ifdef DEBUG_CT_SCAN
              printf("SCAN removing CT stale entry ");
              FILE_output out(stdout);
              showEntry(out, curr);
              printf(" in table list %u\n", i);
#endif
              discardAndRecycle(curr);

            } else {
              // Add entry back to list
              if (preventry) {
                preventry[0].UL = curr;
              } else {
                table[i] = curr;
              }
              preventry = entry;
            }

            curr = next;
          } // while curr
          if (preventry) preventry[0].UL = 0;
      }

      unsigned long convertToList(bool removeStales);
      void listToTable(unsigned long h);

      void scanForStales();
      void rehashTable(const unsigned long* oldT, unsigned long oldS);

      /// Grab space for a new entry
      node_address newEntry(unsigned size);

      // TBD - migrate to 64-bit hash

      static unsigned hash(const ct_entry_key* key);
      static void hash(const ct_entry_type &ET, ct_vector &key);
      static unsigned hash(const ct_entry_type* et, const ct_entry_item* entry);

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
#ifdef DEBUG_CT_SLOTS
        printf(" handle %lu in slot %lu\n", table[h], h);
#else
        printf("\n");
#endif
        fflush(stdout);
#endif
        collisions++;

        discardAndRecycle(table[h]);
        table[h] = curr;
      }


      /**
        Check equality.
        We advance pointer a and return the result portion if equal.
        If unequal we return 0.
      */
      static inline ct_entry_item* equal(ct_entry_item* a, const ct_entry_key* key) {
        const ct_entry_type* et = key->getET();
        MEDDLY_DCASSERT(et);
        //
        // Compare operator, if needed
        //
        if (MONOLITHIC) {
          if (et->getID() != (*a).U) return 0;
          a++;
        }
        //
        // compare #repetitions, if needed
        //
        if (et->isRepeating()) {
          if (key->numRepeats() != (*a).U) return 0;
          a++;
        }
        //
        // Compare key portion only
        //
        const ct_entry_item* b = key->rawData();
        const unsigned klen = et->getKeySize(key->numRepeats());
        for (unsigned i=0; i<klen; i++) {
          const ct_typeID t = et->getKeyType(i).getType();
          switch (t) {
              case ct_typeID::FLOAT:
                        if (a[i].F != b[i].F) return 0;
                        continue;
              case ct_typeID::NODE:
                        if (a[i].N != b[i].N) return 0;
                        continue;
              case ct_typeID::INTEGER:
                        if (a[i].I != b[i].I) return 0;
                        continue;
              case ct_typeID::DOUBLE:
                        if (a[i].D != b[i].D) return 0;
                        continue;
              case ct_typeID::GENERIC:
                        if (a[i].G != b[i].G) return 0;
                        continue;
              case ct_typeID::LONG:
                        if (a[i].L  != b[i].L) return 0;
                        continue;
              default:
                        MEDDLY_DCASSERT(0);
          } // switch t
        } // for i
        return a + klen;
      }



      /**
          Check if the key portion of an entry equals key and should be discarded.
            @param  entry   Complete entry in CT to check.
            @param  key     Key to compare against.
            @param  discard On output: should this entry be discarded

            @return Result portion of the entry, if they key portion matches;
                    0 if the entry does not match the key.
      */
      inline ct_entry_item* checkEqualityAndStatus(ct_entry_item* entry, const ct_entry_key* key, bool &discard)
      {
        ct_entry_item* answer = equal( CHAINED ? (entry+1) : entry, key);
        if (answer)
        {
          //
          // Equal.
          //
          discard = isDead(answer, key->getET());
        } else {
          //
          // Not equal.
          //
          if (checkStalesOnFind) {
            discard = isStale(entry);
          } else {
            discard = false;
          } // if checkStalesOnFind
        }
        return answer;
      }


      /**
            Check equality.
            We advance pointer a and return the result portion if equal.
            If unequal we return null.
      */
      static inline ct_entry_item* equal(ct_entry_item* a,
              const ct_entry_type &ET, const ct_vector &key)
        {
            //
            // Compare operator, if needed
            //
            if (MONOLITHIC) {
                if (ET.getID() != (*a).U) return nullptr;
                a++;
            }
            //
            // compare size, ff needed
            //
            if (ET.isRepeating()) {
                if (key.getSize() != (*a).U) return nullptr;
                a++;
            }
            //
            // Compare key portion only
            //
            for (unsigned i=0; i<key.getSize(); i++) {
                MEDDLY_DCASSERT( ET.getKeyType(i).hasType(key[i].getType()) );
                if (! key[i].equals(a[i])) return nullptr;
            } // for i
            return a + key.getSize();
      }

      /**
          Check if the key portion of an entry equals key and should be discarded.
            @param  entry   Complete entry in CT to check.
            @param  key     Key to compare against.
            @param  ET      Entry type
            @param  discard On output: should this entry be discarded

            @return Result portion of the entry, if they key portion matches;
                    0 if the entry does not match the key.
      */
      inline ct_entry_item* checkEqualityAndStatus(ct_entry_item* entry,
              const ct_entry_type &ET, const ct_vector &key, bool &discard)
        {
            ct_entry_item* answer
                = equal( CHAINED ? (entry+1) : entry, ET, key);
            if (answer)
            {
                //
                // Equal.
                //
                discard = isDead(answer, &ET);
            } else {
                //
                // Not equal.
                //
                if (checkStalesOnFind) {
                    discard = isStale(entry);
                } else {
                    discard = false;
                } // if checkStalesOnFind
            }
            return answer;
        }


      /**
          Check if an entry is stale.
            @param  entry   Pointer to complete entry to check.
            @return         true, if the entry should be discarded.
      */
      bool isStale(const ct_entry_item* entry) const;

      /**
          Check if a result is dead (unrecoverable).
            @param  res   Pointer to result portion of entry.
            @param  et    Entry type information.
            @return       true, if the entry should be discarded.
      */
      bool isDead(const ct_entry_item* res, const ct_entry_type* et) const;

      /**
          Copy a result into an entry.
      */
      inline void setResult(ct_entry_item* respart, const ct_entry_result &res, const ct_entry_type* et) {
        const ct_entry_item* resdata = res.rawData();
        memcpy(respart, resdata, res.dataLength()*sizeof(ct_entry_item));
      }

      /// Display a chain
      inline void showChain(output &s, unsigned long L) const {
        s << L;
        if (CHAINED) {
          while (L) {
#ifdef INTEGRATED_MEMMAN
            const ct_entry_item* entry = entries+L;
#else
            const ct_entry_item* entry = (const int*) MMAN->getChunkAddress(L);
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
      void showKey(output &s, const ct_entry_key* k) const;


    private:
      /// TBD: use vector
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
      ct_entry_item*  entries;
      /// Used entries
      unsigned long entriesSize;
      /// Memory allocated for entries
      unsigned long entriesAlloc;

      static const unsigned maxEntrySize = 15;
      static const unsigned maxEntryBytes = sizeof(ct_entry_item) * maxEntrySize;

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
  const ct_settings &s, operation* op, unsigned slot)
: compute_table(s, op, slot)
{
  if (MONOLITHIC) {
    MEDDLY_DCASSERT(0==op);
    MEDDLY_DCASSERT(0==slot);
  } else {
    MEDDLY_DCASSERT(op);
  }

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
  entries = (ct_entry_item*) malloc(entriesAlloc * sizeof(ct_entry_item) );
  if (0==entries) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  entries[0].L = 0;     // NEVER USED; set here for sanity.

  mstats.incMemUsed( entriesSize * sizeof(ct_entry_item) );
  mstats.incMemAlloc( entriesAlloc * sizeof(ct_entry_item) );
#else
  MEDDLY_DCASSERT(s.MMS);
  MMAN = s.MMS->initManager(sizeof(ct_entry_item), 2, mstats);
#endif

  /*
      Initialize hash table
  */
  tableSize = 1024;
  tableExpand = CHAINED ? 4*1024 : 512;
  tableShrink = 0;
  table = (unsigned long*) malloc(tableSize * sizeof(unsigned long));
  if (0==table) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  for (unsigned long i=0; i<tableSize; i++) table[i] = 0;

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
inline MEDDLY::ct_entry_item* MEDDLY::ct_none<MONOLITHIC, CHAINED>
::findEntry(ct_entry_key* key)
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

  ct_entry_item* answer = 0;
  ct_entry_item* preventry = 0;
  unsigned long hcurr = key->getHash() % tableSize;
  unsigned long curr = table[hcurr];

#ifdef DEBUG_CT_SEARCHES
  printf("Searching for CT entry ");
  FILE_output out(stdout);
  showKey(out, key);
#ifdef DEBUG_CT_SLOTS
  printf(" in hash slot %lu\n", hcurr);
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
    ct_entry_item* entry = entries + curr;
#else
    ct_entry_item* entry = (ct_entry_item*) MMAN->getChunkAddress(curr);
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
        printf(" handle %lu in slot %lu\n", curr, hcurr);
#else
        printf("\n");
#endif
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
#ifdef DEBUG_CT_SLOTS
        printf(" handle %lu in slot %lu\n", curr, hcurr);
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

template <bool MONOLITHIC, bool CHAINED>
inline MEDDLY::ct_entry_item* MEDDLY::ct_none<MONOLITHIC, CHAINED>
::findEntry(const ct_entry_type &ET, const ct_vector &key)
{
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

    ct_entry_item* answer = nullptr;
    ct_entry_item* preventry = nullptr;
    unsigned long hcurr = key.getHash() % tableSize;
    unsigned long curr = table[hcurr];

#ifdef DEBUG_CT_SEARCHES
    ostream_output out(std::cout);
    std::cout << "Searching for CT entry ";
    showKey(out, key);
#ifdef DEBUG_CT_SLOTS
    std::cout << " in hash slot " << hcurr;
#endif
    std::cout << std::endl;
#endif

    for (chain=0; ; ) {

        if (!curr) {
            if (CHAINED) break;
            if (++chain > maxCollisionSearch) break;
            incMod(hcurr);
            curr = table[hcurr];
            continue;
        } else {
            if (CHAINED) chain++;
        }

#ifdef INTEGRATED_MEMMAN
        ct_entry_item* entry = entries + curr;
#else
        ct_entry_item* entry = (ct_entry_item*) MMAN->getChunkAddress(curr);
#endif

        bool discard;
        answer = checkEqualityAndStatus(entry, ET, key, discard);

        if (discard) {
            //
            // Delete the entry.
            //
#ifdef DEBUG_CT
            if (answer) std::cout << "Removing stale CT hit   ";
            else        std::cout << "Removing stale CT entry ";
            ostream_output out(std::cout);
            showEntry(out, curr);
#ifdef DEBUG_CT_SLOTS
            std::cout << " handle " << curr << " in slot " << hcurr;
#endif
            std::cout << std::endl;
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
                answer = nullptr;
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
                // If the matching entry is not already at the front
                // of the list, move it there
                //
                if (preventry) {
                    preventry[0].UL = entry[0].UL;
                    entry[0].UL = table[hcurr];
                    table[hcurr] = curr;
                }
            }
#ifdef DEBUG_CT
            std::cout << "Found CT entry ";
            ostream_output out(std::cout);
            showEntry(out, curr);
#ifdef DEBUG_CT_SLOTS
            std::cout << " handle " << curr << " in slot " << hcurr;
#endif
            std::cout << std::endl;
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

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>
::find(ct_entry_key *key, ct_entry_result& res)
{
  //
  // TBD - probably shouldn't hash the floats, doubles.
  //

  //
  // Hash the key
  //
  setHash(key, hash(key));

  ct_entry_item* entry_result = findEntry(key);
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
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::addEntry(ct_entry_key* key, const ct_entry_result &res)
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

  const ct_entry_type* et = key->getET();
  MEDDLY_DCASSERT(et);

  //
  // Allocate an entry
  //

  const unsigned num_slots =
    et->getKeySize(key->numRepeats()) +
    + et->getResultSize()
    + (CHAINED ? 1 : 0)
    + (MONOLITHIC ? 1 : 0)
    + (et->isRepeating() ? 1 : 0)
  ;
  unsigned long curr = newEntry(num_slots);

#ifdef INTEGRATED_MEMMAN
  ct_entry_item* entry = entries + curr;
#else
  ct_entry_item* entry = (ct_entry_item*) MMAN->getChunkAddress(curr);
#endif

  //
  // Copy into the entry
  //
  ct_entry_item* key_portion = entry + (CHAINED ? 1 : 0);
  if (MONOLITHIC) {
    (*key_portion).U = et->getID();
    key_portion++;
  }
  if (et->isRepeating()) {
    (*key_portion).U = key->numRepeats();
    key_portion++;
  }
  const unsigned key_slots = et->getKeySize(key->numRepeats());
  memcpy(key_portion, key->rawData(), key_slots * sizeof(ct_entry_item));
  ct_entry_item* res_portion = key_portion + key_slots;
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
#ifdef DEBUG_CT_SLOTS
  printf(" handle %lu in slot %lu (%u slots long)\n", curr, h, num_slots);
#else
  printf("\n");
#endif
  fflush(stdout);
#endif

  if (perf.numEntries < tableExpand) return;

  //
  // Time to GC and maybe resize the table
  //
  perf.resizeScans++;

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
      tableExpand = std::numeric_limits<int>::max();
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
        tableExpand = std::numeric_limits<int>::max();
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
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::updateEntry(ct_entry_key* key, const ct_entry_result &res)
{
  MEDDLY_DCASSERT(key->getET()->isResultUpdatable());
  ct_entry_item* entry_result = findEntry(key);
  if (!entry_result) {
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }

  //
  // decrement cache counters for old result,
  //

  const ct_entry_type* et = key->getET();
  for (unsigned i=0; i<et->getResultSize(); i++) {
    const ct_itemtype &item = et->getResultType(i);
    if (item.hasForest()) {
        item.getForest()->uncacheNode( entry_result[i].N );
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
bool MEDDLY::ct_none<MONOLITHIC, CHAINED>
::find(const ct_entry_type &ET, ct_vector &key, ct_vector &res)
{
    //
    // TBD - probably shouldn't hash the floats, doubles.
    //

    //
    // Hash the key
    //
    hash(ET, key);

    ct_entry_item* entry_result = findEntry(ET, key);
    perf.pings++;
    if (!entry_result) return false;

    perf.hits++;
    //
    // Fill res
    //
    MEDDLY_DCASSERT(res.getSize() == ET.getResultSize());
    for (unsigned i=0; i<ET.getResultSize(); i++) {
        res[i].set(ET.getResultType(i).getType(), entry_result[i]);
    }
    return true;
}


// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::addEntry(const ct_entry_type &ET,
        ct_vector &key, const ct_vector &res)
{
    if (!MONOLITHIC) {
        if (&ET != global_et)
            throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }

    //
    // Increment cache counters for nodes in the entry
    //
    for (unsigned i=0; i<key.getSize(); i++) {
        const ct_itemtype &item = ET.getKeyType(i);
        if (item.hasForest()) {
            item.getForest()->cacheNode(key[i].getN());
        }
    }
    for (unsigned i=0; i<res.getSize(); i++) {
        const ct_itemtype &item = ET.getResultType(i);
        if (item.hasForest()) {
            item.getForest()->cacheNode(res[i].getN());
        }
    }

    unsigned long h = key.getHash() % tableSize;

    //
    // Allocate an entry
    //

    const unsigned num_slots =
        key.getSize() +
        + res.getSize()
        + (CHAINED ? 1 : 0)
        + (MONOLITHIC ? 1 : 0)
        + (ET.isRepeating() ? 1 : 0)
    ;
    unsigned long curr = newEntry(num_slots);
#ifdef INTEGRATED_MEMMAN
    ct_entry_item* entry = entries + curr;
#else
    ct_entry_item* entry = (ct_entry_item*) MMAN->getChunkAddress(curr);
#endif

    //
    // Copy into the entry
    //
    ct_entry_item* key_portion = entry + (CHAINED ? 1 : 0);
    if (MONOLITHIC) {
        (*key_portion).U = ET.getID();
        key_portion++;
    }
    if (ET.isRepeating()) {
        (*key_portion).U = key.getSize();
        key_portion++;
    }
    ct_entry_item* res_portion = key_portion + key.getSize();
    //
    // copy key
    //
    for (unsigned i=0; i<key.getSize(); i++) {
        key[i].get( key_portion[i] );
    }
    //
    // copy result
    //
    for (unsigned i=0; i<res.getSize(); i++) {
        res[i].get( res_portion[i] );
    }

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
    std::cout << "Added CT entry ";
    ostream_output out(std::cout);
    showEntry(out, curr);
#ifdef DEBUG_CT_SLOTS
    std::cout   << " handle " << curr << " ins slot " << h
                << " (" << num_slots << " slots long";
#endif
    std::cout << std::endl;
#endif

    if (perf.numEntries < tableExpand) return;

    //
    // Time to GC and maybe resize the table
    //
    perf.resizeScans++;

#ifdef DEBUG_SLOW
    std::cout   << "Running GC in compute table (size "
                << tableSize << ", entries " << perf.numEntries << ")\n";
#endif

    unsigned long list = 0;
    if (CHAINED) {
        list = convertToList(checkStalesOnResize);
        if (perf.numEntries < tableSize) {
            listToTable(list);
#ifdef DEBUG_SLOW
            std::cout   << "Done CT GC, no resizing (now entries "
                        << perf.numEntries << ")\n";
#endif
            return;
        }
    } else {
        scanForStales();
        if (perf.numEntries < tableExpand / 4) {
#ifdef DEBUG_SLOW
            std::cout   << "Done CT GC, no resizing (now entries "
                        << perf.numEntries << ")\n";
#endif
            return;
        }
    }

    //
    // Resize needed
    //

    unsigned long newsize = tableSize * 2;
    if (newsize > maxSize) newsize = maxSize;

    if (CHAINED) {
        if (newsize != tableSize) {
            //
            // Enlarge table
            //

            unsigned long* newt = (unsigned long*)
                realloc(table, newsize * sizeof(unsigned long));
            if (!newt) {
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
            tableExpand = std::numeric_limits<int>::max();
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
            if (!table) {
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
                tableExpand = std::numeric_limits<int>::max();
            } else {
                tableExpand = tableSize / 2;
            }
            tableShrink = tableSize / 8;
        }
    } // if CHAINED

#ifdef DEBUG_SLOW
    std::cout << "CT enlarged to size " << tableSize << std::endl;
#endif

}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::removeStales()
{
#ifdef DEBUG_REMOVESTALES
  fprintf(stdout, "Removing stales in CT (size %lu, entries %lu)\n",
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
void MEDDLY::ct_none<MONOLITHIC, CHAINED>::removeAll()
{
#ifdef DEBUG_REMOVESTALES
  fprintf(stdout, "Removing all CT entries\n");
#endif
  for (unsigned long i=0; i<tableSize; i++) {
    while (table[i]) {
      unsigned long curr = table[i];
#ifdef INTEGRATED_MEMMAN
      const ct_entry_item* entry = entries + curr;
#else
      const ct_entry_item* entry = (const int*) MMAN->getChunkAddress(curr);
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
  s.put("", 6);
  s << "Resize (GC) scans   :\t" << long(perf.resizeScans) << "\n";

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

  for (unsigned long i=0; i<tableSize; i++) {
    if (0==table[i]) continue;
    s << "table[";
    s.put(long(i), 9);
    s << "]: ";
    showChain(s, table[i]);
  }

  if (--verbLevel < 1) return;

  s << "\nHash table nodes:\n";

  for (unsigned long i=0; i<tableSize; i++) {
    unsigned long curr = table[i];
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
void MEDDLY::ct_none<MONOLITHIC, CHAINED>
::countNodeEntries(const forest* f, size_t* counts) const
{
#ifdef DEBUG_VALIDATE_COUNTS
  printf("    Counting in ct_none\n");
#endif
  const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

  for (unsigned long i=0; i<tableSize; i++) {
    //
    // Process table[i]
    //
    unsigned long curr = table[i];
    while (curr) {
      //
      // Loop through the chain
      // (unchained: done properly)
      //

#ifdef INTEGRATED_MEMMAN
      ct_entry_item* entry = entries + curr;
#else
      ct_entry_item* entry = (ct_entry_item*) MMAN->getChunkAddress(curr);
#endif
      const ct_entry_type* et = MONOLITHIC
        ?   getEntryType(entry[CHAINED ? 1 : 0].U)
        :   global_et;
      MEDDLY_DCASSERT(et);

      const ct_entry_item* ptr = entry + SHIFT;
      unsigned reps;
      if (et->isRepeating()) {
        reps = (*ptr).U;
        ptr++;
      } else {
        reps = 0;
      }

      //
      // Count the key portion
      //
      const unsigned stop = et->getKeySize(reps);
      for (unsigned i=0; i<stop; i++) {
        if (! et->getKeyType(i).hasForest(f)) continue;
        if (ptr[i].N > 0) {
          ++counts[ ptr[i].N ];
        }
      } // for i
      ptr += stop;

      //
      // Count the result portion
      //
      for (unsigned i=0; i<et->getResultSize(); i++) {
        if (! et->getResultType(i).hasForest(f)) continue;
        if (ptr[i].N > 0) {
          ++counts[ ptr[i].N ];
        }
      } // for i;

      //
      // Next in the chain
      //
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
      ct_entry_item* entry = entries + curr;
#else
      ct_entry_item* entry = (ct_entry_item*) MMAN->getChunkAddress(curr);
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
  MEDDLY_DCASSERT(CHAINED);
#ifdef DEBUG_LIST2TABLE
  printf("Recovering  list");
  // showChain(stdout, L);
#endif
  while (L) {
#ifdef INTEGRATED_MEMMAN
    ct_entry_item* entry = entries+L;
#else
    ct_entry_item* entry = (ct_entry_item*) MMAN->getChunkAddress(L);
#endif
    const unsigned long curr = L;
    L = entry[0].UL;

    const ct_entry_type* et = MONOLITHIC
      ?   getEntryType(entry[1].U)
      :   global_et;
    MEDDLY_DCASSERT(et);

    const unsigned h = hash(et, entry + 1) % tableSize;
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
    scanForStales(i);
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
  for (unsigned long i=0; i<oldS; i++) {
    unsigned long curr = oldT[i];
    if (0==curr) continue;
#ifdef INTEGRATED_MEMMAN
    const ct_entry_item* entry = entries + curr;
#else
    const ct_entry_item* entry = (const ct_entry_item*) MMAN->getChunkAddress(curr);
#endif

    const ct_entry_type* et = MONOLITHIC
      ?   getEntryType(entry[0].UL)
      :   global_et;
    MEDDLY_DCASSERT(et);

    unsigned long h = hash(et, entry) % tableSize;
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
    mstats.incMemUsed( size * sizeof(ct_entry_item) );
    return h;
  }
  if (entriesSize + size > entriesAlloc) {
    // Expand by a factor of 1.5
    unsigned neA = entriesAlloc + (entriesAlloc/2);
    ct_entry_item* ne = (ct_entry_item*) realloc(entries, neA * sizeof(ct_entry_item));
    if (0==ne) {
      fprintf(stderr,
          "Error in allocating array of size %lu at %s, line %d\n",
          neA * sizeof(int), __FILE__, __LINE__);
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    mstats.incMemAlloc( (entriesAlloc / 2) * sizeof(ct_entry_item) );
    entries = ne;
    entriesAlloc = neA;
  }
  MEDDLY_DCASSERT(entriesSize + size <= entriesAlloc);
  node_address h = entriesSize;
  entriesSize += size;
#ifdef DEBUG_CTALLOC
  fprintf(stderr, "New entry %lu size %u\n", h, size);
#endif
  mstats.incMemUsed( size * sizeof(ct_entry_item) );
  return h;

#else
  perf.numEntries++;
  size_t the_size = size;
  return MMAN->requestChunk(the_size);
#endif
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
unsigned MEDDLY::ct_none<MONOLITHIC, CHAINED>
::hash(const ct_entry_key* key)
{
  const ct_entry_type* et = key->getET();
  MEDDLY_DCASSERT(et);

  hash_stream H;
  H.start();

  //
  // Operator, if needed
  //
  if (MONOLITHIC) {
    H.push(et->getID());
  }
  //
  // #repetitions, if needed
  //
  if (et->isRepeating()) {
    H.push(key->numRepeats());
  }

  //
  // Hash key portion only
  //

  const ct_entry_item* entry = key->rawData();
  const unsigned klen = et->getKeySize(key->numRepeats());
  for (unsigned i=0; i<klen; i++) {    // i initialized earlier
    const ct_typeID t = et->getKeyType(i).getType();
    switch (t) {
        case ct_typeID::FLOAT:
                        MEDDLY_DCASSERT(sizeof(entry[i].F) == sizeof(entry[i].U));

        case ct_typeID::NODE:
                        MEDDLY_DCASSERT(sizeof(entry[i].N) == sizeof(entry[i].U));
        case ct_typeID::INTEGER:
                        MEDDLY_DCASSERT(sizeof(entry[i].I) == sizeof(entry[i].U));
                        H.push(entry[i].U);
                        continue;

        case ct_typeID::DOUBLE:
                        MEDDLY_DCASSERT(sizeof(entry[i].D) == sizeof(entry[i].L));
        case ct_typeID::GENERIC:
                        MEDDLY_DCASSERT(sizeof(entry[i].G) == sizeof(entry[i].L));
        case ct_typeID::LONG:
                        {
                          unsigned* hack = (unsigned*) (& (entry[i].L));
                          H.push(hack[0], hack[1]);
                        }
                        continue;
        default:
                        MEDDLY_DCASSERT(0);
    } // switch t
  } // for i

  return H.finish();
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
void MEDDLY::ct_none<MONOLITHIC, CHAINED>
::hash(const ct_entry_type &ET, ct_vector &key)
{
    hash_stream H;
    H.start();
    //
    // Operator, if needed
    //
    if (MONOLITHIC) {
        H.push(ET.getID());
    }
    //
    // key size, if needed
    //
    if (ET.isRepeating()) {
        H.push(key.getSize());
    }

    //
    // Hash key
    //
    for (unsigned i=0; i<key.getSize(); i++) {
        MEDDLY_DCASSERT(ET.getKeyType(i).hasType(key[i].getType()));
        key[i].hash(H);
    }

    //
    // Set hash
    //
    key.setHash( H.finish() );
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
unsigned MEDDLY::ct_none<MONOLITHIC, CHAINED>
::hash(const ct_entry_type* et, const ct_entry_item* entry)
{
  hash_stream H;
  H.start();
  //
  // Operator, if needed
  //
  if (MONOLITHIC) {
    MEDDLY_DCASSERT(et->getID() == entry[0].U);
    H.push(et->getID());
    entry++;
  }
  //
  // #repetitions, if needed
  //
  unsigned reps = 0;
  if (et->isRepeating()) {
    H.push(entry[0].U);
    reps = entry[0].U;
    entry++;
  }

  //
  // Hash key portion only
  //

  const unsigned klen = et->getKeySize(reps);
  for (unsigned i=0; i<klen; i++) {
    const ct_typeID t = et->getKeyType(i).getType();
    switch (t) {
        case ct_typeID::FLOAT:
                        MEDDLY_DCASSERT(sizeof(entry[i].F) == sizeof(entry[i].U));
        case ct_typeID::NODE:
                        MEDDLY_DCASSERT(sizeof(entry[i].N) == sizeof(entry[i].U));
        case ct_typeID::INTEGER:
                        MEDDLY_DCASSERT(sizeof(entry[i].I) == sizeof(entry[i].U));
                        H.push(entry[i].U);
                        continue;

        case ct_typeID::DOUBLE:
                        MEDDLY_DCASSERT(sizeof(entry[i].D) == sizeof(entry[i].L));
        case ct_typeID::GENERIC:
                        MEDDLY_DCASSERT(sizeof(entry[i].G) == sizeof(entry[i].L));
        case ct_typeID::LONG:      {
                          unsigned* hack = (unsigned*) (& (entry[i].L));
                          H.push(hack[0], hack[1]);
                        }
                        continue;
        default:
                        MEDDLY_DCASSERT(0);
    } // switch t
  } // for i

  return H.finish();
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
bool MEDDLY::ct_none<MONOLITHIC, CHAINED>
::isStale(const ct_entry_item* entry) const
{
  const int SHIFT = (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);

  //
  // Check entire entry for staleness
  //
  const ct_entry_type* et = MONOLITHIC
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
    const ct_itemtype &item = et->getKeyType(i);
    if (item.hasForest()) {
      if (item.getForest()->isStaleEntry(entry[i].N)) {
        return true;
      } else {
        item.getForest()->setCacheBit(entry[i].N);
      }
    }
  } // for i

  entry += klen;

  //
  // Result portion
  //
  for (unsigned i=0; i<et->getResultSize(); i++) {
    const ct_itemtype &item = et->getResultType(i);
    if (item.hasForest()) {
      if (item.getForest()->isStaleEntry(entry[i].N)) {
        return true;
      } else {
        item.getForest()->setCacheBit(entry[i].N);
      }
    }
  } // for i

  return false;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED>
bool MEDDLY::ct_none<MONOLITHIC, CHAINED>
::isDead(const ct_entry_item* result, const ct_entry_type* et) const
{
  MEDDLY_DCASSERT(et);
  MEDDLY_DCASSERT(result);
  //
  // Check result portion for dead nodes - cannot use result in that case
  //
  for (unsigned i=0; i<et->getResultSize(); i++) {
    const ct_itemtype &item = et->getResultType(i);
    if (item.hasForest()) {
      if (item.getForest()->isDeadEntry(result[i].N)) {
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
  ct_entry_item* entry = entries + h;
#else
  ct_entry_item* entry = (ct_entry_item*) MMAN->getChunkAddress(table[h]);
#endif

  unsigned slots;

  //
  // Discard.
  // For the portion of the entry that are nodes,
  // we need to notify the forest to decrement
  // the CT counter.
  //

  const ct_entry_type* et = MONOLITHIC
    ?   getEntryType(entry[CHAINED ? 1 : 0].U)
    :   global_et;
  MEDDLY_DCASSERT(et);

  const ct_entry_item* ptr = entry + SHIFT;
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
    const ct_itemtype &item = et->getKeyType(i);
    if (item.hasForest()) {
      item.getForest()->uncacheNode( ptr[i].N );
      continue;
    }
    if (item.hasType(ct_typeID::GENERIC)) {
      delete ptr[i].G;
      continue;
    }
    MEDDLY_DCASSERT(! item.hasType(ct_typeID::NODE) );
  } // for i

  ptr += stop;

  //
  // Result portion
  //
  for (unsigned i=0; i<et->getResultSize(); i++) {
    const ct_itemtype &item = et->getResultType(i);
    if (item.hasForest()) {
      item.getForest()->uncacheNode( ptr[i].N );
      continue;
    }
    if (item.hasType(ct_typeID::GENERIC)) {
      delete ptr[i].G;
      continue;
    }
    MEDDLY_DCASSERT(! item.hasType(ct_typeID::NODE) );
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
  mstats.decMemUsed( slots * sizeof(ct_entry_item) );
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
  const ct_entry_item* entry = entries+h;
#else
  const ct_entry_item* entry = (const ct_entry_item*) MMAN->getChunkAddress(h);
#endif

  const ct_entry_type* et = MONOLITHIC
    ?   getEntryType(entry[CHAINED?1:0].U)
    :   global_et;
  MEDDLY_DCASSERT(et);

  const ct_entry_item* ptr = entry + (MONOLITHIC ? 1 : 0) + (CHAINED ? 1 : 0);
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
      switch (et->getKeyType(i).getType()) {
          case ct_typeID::NODE:
                        s.put(long(ptr[i].N));
                        break;
          case ct_typeID::INTEGER:
                        s.put(long(ptr[i].I));
                        break;
          case ct_typeID::LONG:
                        s.put(ptr[i].L);
                        break;
          case ct_typeID::FLOAT:
                        s.put(ptr[i].F, 0, 0, 'e');
                        break;
          case ct_typeID::DOUBLE:
                        s.put(ptr[i].D, 0, 0, 'e');
                        break;
          case ct_typeID::GENERIC:
                        s.put_hex((unsigned long)ptr[i].G);
                        break;
        default:
                        MEDDLY_DCASSERT(0);
      } // switch et->getKeyType()
  } // for i
  ptr += stop;
  s << "): ";
  for (unsigned i=0; i<et->getResultSize(); i++) {
      if (i) s << ", ";
      switch (et->getResultType(i).getType()) {
          case ct_typeID::NODE:
                        s.put(long(ptr[i].N));
                        break;
          case ct_typeID::INTEGER:
                        s.put(long(ptr[i].I));
                        break;
          case ct_typeID::LONG:
                        s.put(ptr[i].L);
                        s << "(L)";
                        break;
          case ct_typeID::FLOAT:
                        s.put(ptr[i].F, 0, 0, 'e');
                        break;
          case ct_typeID::DOUBLE:
                        s.put(ptr[i].D, 0, 0, 'e');
                        break;

          case ct_typeID::GENERIC:
                        s.put_hex((unsigned long)ptr[i].G);
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
::showKey(output &s, const ct_entry_key* key) const
{
  MEDDLY_DCASSERT(key);

  const ct_entry_type* et = key->getET();
  MEDDLY_DCASSERT(et);

  unsigned reps = key->numRepeats();
  s << "[" << et->getName() << "(";
  unsigned stop = et->getKeySize(reps);

  for (unsigned i=0; i<stop; i++) {
      ct_entry_item item = key->rawData()[i];
      if (i) s << ", ";
      switch (et->getKeyType(i).getType()) {
          case ct_typeID::NODE:
                        s.put(long(item.N));
                        break;
          case ct_typeID::INTEGER:
                        s.put(long(item.I));
                        break;
          case ct_typeID::LONG:
                        s.put(item.L);
                        break;
          case ct_typeID::FLOAT:
                        s.put(item.F);
                        break;
          case ct_typeID::DOUBLE:
                        s.put(item.D);
                        break;
          case ct_typeID::GENERIC:
                        s.put_hex((unsigned long)item.G);
                        break;
        default:
                        MEDDLY_DCASSERT(0);
      } // switch et->getKeyType()
  } // for i
  s << "): ?]";
}


#endif  // include guard


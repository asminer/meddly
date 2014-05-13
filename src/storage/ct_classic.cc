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


#include "ct_classic.h"

// #define DEBUG_SLOW
// #define DEBUG_CT
// #define DEBUG_TABLE2LIST
// #define DEBUG_LIST2TABLE
// #define DEBUG_CTALLOC

// #define DEBUG_REMOVESTALES
// #define SUMMARY_STALES


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                 base_table::old_search_key methods                 *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::base_table::old_search_key::old_search_key(operation* op)
 : search_key(op)
{
  // hashLength = 0;
  data = 0;
  key_data = 0;
  killData = false;
}

MEDDLY::base_table::old_search_key::~old_search_key()
{
  if (killData) delete[] data;
}

// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                         base_table methods                         *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::base_table::base_table(const settings::computeTableSettings &s)
 : compute_table(s)
{
  entriesAlloc = 1024;
  entries = (int*) malloc(entriesAlloc * sizeof(int));
  entriesSize = 1;
  // entries[0] is never, ever, used.
  if (0==entries) throw error(error::INSUFFICIENT_MEMORY);
  // for recycling entries
  freeList = new int[1+maxEntrySize];
  for (int i=0; i<=maxEntrySize; i++) freeList[i] = 0;

  currMemory = entriesAlloc * sizeof(int) + (1+maxEntrySize) * sizeof(int);
  peakMemory = currMemory;
}

MEDDLY::base_table::~base_table()
{
  free(entries);
  delete[] freeList;
}

int MEDDLY::base_table::newEntry(int size)
{
  // check free list
  if (size > maxEntrySize) {
    fprintf(stderr, "MEDDLY error: request for compute table entry larger than max size\n");
    throw error(error::MISCELLANEOUS);  // best we can do
  }
  if (size<1) return 0;
  perf.numEntries++;
  if (freeList[size]) {
    int h = freeList[size];
    freeList[size] = entries[h];
#ifdef DEBUG_CTALLOC
    fprintf(stderr, "Re-used entry %d size %d\n", h, size);
#endif
    return h;
  }
  if (entriesSize + size > entriesAlloc) {
    // Expand by a factor of 1.5
    int neA = entriesAlloc + (entriesAlloc/2);
    int* ne = (int*) realloc(entries, neA * sizeof(int));
    if (0==ne) throw error(error::INSUFFICIENT_MEMORY);
    currMemory += (neA - entriesAlloc) * sizeof(int);
    if (currMemory > peakMemory) peakMemory = currMemory;
    entries = ne;
    entriesAlloc = neA;
  }
  MEDDLY_DCASSERT(entriesSize + size <= entriesAlloc);
  int h = entriesSize;
  entriesSize += size;
#ifdef DEBUG_CTALLOC
  fprintf(stderr, "New entry %d size %d\n", h, size);
#endif
  return h;
}

void MEDDLY::base_table::dumpInternal(FILE* s, int verbLevel) const
{
  if (verbLevel < 1) return;
  if (0==entries) fprintf(s, "Entries: null\n");
  else {
    fprintf(s, "Entries: [%d", entries[0]);
    for (int i=1; i<entriesSize; i++) 
      fprintf(s, ", %d", entries[i]);
    fprintf(s, "]\n");
  }
  if (0==freeList) fprintf(s, "Free: null\n");
  else {
    fprintf(s, "Free: [%d", freeList[0]);
    for (int i=1; i<=maxEntrySize; i++) 
      fprintf(s, ", %d", freeList[i]);
    fprintf(s, "]\n");
  }
}

void MEDDLY::base_table
::report(FILE* s, int indent, int &level) const
{
  if (level < 1) return;
  fprintf(s, "%*sNumber of entries :\t%u\n", indent, "", perf.numEntries);
  fprintf(s, "%*sEntry array size  :\t%d\n", indent, "", entriesSize);
  fprintf(s, "%*sEntry array alloc :\t%d\n", indent, "", entriesAlloc);

  if (--level < 1) return;

  fprintf(s, "%*sPings             :\t%d\n", indent, "", perf.pings);
  fprintf(s, "%*sHits              :\t%d\n", indent, "", perf.hits);

  if (--level < 1) return;

  fprintf(s, "%*sSearch length histogram:\n", indent, "");
  for (int i=0; i<stats::searchHistogramSize; i++) {
    if (perf.searchHistogram[i]) {
      fprintf(s, "%*s%3d: %ld\n", indent+4, "", i, perf.searchHistogram[i]);
    }
  }
  if (perf.numLargeSearches)
    fprintf(s, "%*sSearches longer than %d: %ld\n", indent, "",
            stats::searchHistogramSize-1, perf.numLargeSearches
    );
  fprintf(s, "%*sMax search length: %d\n", indent, "", perf.maxSearchLength);
}

// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                         base_hash  methods                         *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::base_hash::base_hash(const settings::computeTableSettings &s, 
  int initTsz, int initTex) : base_table(s)
{
  tableSize = initTsz;
  tableExpand = initTex;
  tableShrink = 0;
  table = (int*) malloc(tableSize * sizeof(int));
  if (0==table) throw error(error::INSUFFICIENT_MEMORY);
  for (unsigned i=0; i<tableSize; i++) table[i] = 0;

  currMemory += tableSize * sizeof(int);
  peakMemory = currMemory;
}

MEDDLY::base_hash::~base_hash()
{
  free(table);
}

void MEDDLY::base_hash::dumpInternal(FILE* s, int verbLevel) const
{
  if (verbLevel < 1) return;
  if (0==table) fprintf(s, "Table: null\n");
  else {
    fprintf(s, "Table: [%d", table[0]);
    for (unsigned i=1; i<tableSize; i++) 
      fprintf(s, ", %d", table[i]);
    fprintf(s, "]\n");
  }
  base_table::dumpInternal(s, verbLevel-1);
}

void MEDDLY::base_hash
::report(FILE* s, int indent, int &level) const
{
  if (level < 1) return;
  fprintf(s, "%*sHash table size   :\t%d\n", indent, "", tableSize);
  base_table::report(s, indent, level);
}

unsigned MEDDLY::base_hash::raw_hash(const int* k, int length)
{
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



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                        base_chained methods                        *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::base_chained::base_chained(const settings::computeTableSettings &s)
 : base_hash(s, 1024, 4*1024)
{
}

MEDDLY::base_chained::~base_chained()
{
}

void MEDDLY::base_chained::addEntry()
{
  unsigned h = currEntry.getHash() % tableSize;
  currEntry.data(0) = table[h];
  table[h] = currEntry.readHandle();

#ifdef DEBUG_CT
  printf("Adding CT entry ");
  showEntry(stdout, currEntry.readHandle());
  // fprintf(stderr, " to slot %u", h);
  printf("\n");
#endif

  if (perf.numEntries < tableExpand) return;

#ifdef DEBUG_SLOW
  fprintf(stdout, "Running GC in compute table (size %d, entries %ld)\n", 
    tableSize, perf.numEntries
  );
#endif

  int list = convertToList(checkStalesOnResize);
  if (perf.numEntries < tableSize) {
    // Don't need to expand
    listToTable(list);
#ifdef DEBUG_SLOW
    fprintf(stdout, "Done CT GC, no resizing (now entries %ld)\n", 
      perf.numEntries
    );
#endif
    return;
  } 

  unsigned newsize = tableSize*2;
  if (newsize > maxSize) newsize = maxSize;

  int* newt = (int*) realloc(table, newsize * sizeof(int));
  if (0==newt) throw error(error::INSUFFICIENT_MEMORY);

  for (unsigned i=tableSize; i<newsize; i++) newt[i] = 0;

  currMemory += (newsize - tableSize) * sizeof(int);
  if (currMemory > peakMemory) peakMemory = currMemory;

  table = newt;
  tableSize = newsize;
  if (tableSize == maxSize) {
    tableExpand = INT_MAX;
  } else {
    tableExpand = 4*tableSize;
  }
  tableShrink = tableSize / 2;

  listToTable(list);
#ifdef DEBUG_SLOW
  fprintf(stdout, "CT enlarged to size %d\n", tableSize);
#endif
}

void MEDDLY::base_chained::removeStales()
{
#ifdef DEBUG_SLOW
  fprintf(stdout, "Removing stales in CT (size %d, entries %ld)\n", 
    tableSize, perf.numEntries
  );
#endif
  int list = convertToList(true);
  if (perf.numEntries < tableShrink) {
    // shrink table
    int newsize = tableSize / 2;
    if (newsize < 1024) newsize = 1024;
    int* newt = (int*) realloc(table, newsize * sizeof(int));
    if (0==newt) throw error(error::INSUFFICIENT_MEMORY); 

    currMemory -= (tableSize - newsize) * sizeof(int);  

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
#ifdef DEBUG_SLOW
  fprintf(stdout, "Done removing CT stales (size %d, entries %ld)\n", 
    tableSize, perf.numEntries
  );
#endif
}

void MEDDLY::base_chained::dumpInternal(FILE* s, int verbLevel) const
{
  if (verbLevel < 1) return;

  fprintf(s, "Hash table chains:\n");

  for (unsigned i=0; i<tableSize; i++) {
    if (0==table[i]) continue;
    fprintf(s, "table[%9d]: %d", i, table[i]);
    int curr = entries[table[i]];
    while (curr) {
      fprintf(s, " -> %d", curr);
      curr = entries[curr];
    }
    fprintf(s, "\n");
  }

  fprintf(s, "\nHash table nodes:\n");
  
  for (unsigned i=0; i<tableSize; i++) {
    int curr = table[i];
    while (curr) {
      fprintf(s, "\tNode %9d:  ", curr);
      showEntry(s, curr);
      fprintf(s, "\n");
      curr = entries[curr];
    }
  }
  fprintf(s, "\n");

  base_hash::dumpInternal(s, verbLevel-1);
}


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                     monolithic_chained methods                     *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::monolithic_chained
::monolithic_chained(const settings::computeTableSettings &s)
 : base_chained(s)
{
}

MEDDLY::monolithic_chained::~monolithic_chained()
{
}

MEDDLY::compute_table::search_key* 
MEDDLY::monolithic_chained::initializeSearchKey(operation* op)
{
  return init(op, 1);
}


MEDDLY::compute_table::search_result& 
MEDDLY::monolithic_chained::find(search_key *k)
{
  static old_search_result ANS;
  old_search_key* key = smart_cast <old_search_key*>(k);
  MEDDLY_DCASSERT(key);
  perf.pings++;
  unsigned h = hash(key);
  int prev = 0;
  int curr = table[h];
  int chain = 0;
  ANS.setInvalid();
  while (curr) {
    chain++;
    //
    // Check for match
    //
    if (equal_sw(entries+curr+1, key->rawData(), key->dataLength())) {
      if (key->getOp()->shouldStaleCacheHitsBeDiscarded()) {
        if (checkStale(h, prev, curr)) {
          // The match is stale.
          // Since there can NEVER be more than one match
          // in the table, we're done!
          break;
        }
      }
      // "Hit"
      perf.hits++;
      if (prev) {
        // not at the front; move it there
        entries[prev] = entries[curr];
        entries[curr] = table[h];
        table[h] = curr;
      }
#ifdef DEBUG_CT
      printf("Found CT entry ");
      key->getOp()->showEntry(stdout, entries + curr + 2);
      // fprintf(stderr, " in slot %u", h);
      printf("\n");
#endif
      ANS.setResult(entries+curr+1+key->dataLength(), key->getOp()->getAnsLength());
      break;
    };
    //
    // No match; maybe check stale
    //
    if (checkStalesOnFind) {
      if (checkStale(h, prev, curr)) continue;
    }
    // advance pointers
    prev = curr;
    curr = entries[curr];
  }
  sawSearch(chain);
  return ANS;
}


MEDDLY::compute_table::entry_builder& 
MEDDLY::monolithic_chained::startNewEntry(search_key *key)
{
  MEDDLY_DCASSERT(key);
  startIndexedEntry(smart_cast<old_search_key*>(key), 1, 1);
  return currEntry;
}

void MEDDLY::monolithic_chained::removeAll()
{
  for (unsigned i=0; i<tableSize; i++) {
    while (table[i]) {
      int curr = table[i];
      table[i] = entries[curr];
      operation* currop = operation::getOpWithIndex(entries[curr+1]);
      MEDDLY_DCASSERT(currop);
      currop->discardEntry(entries + curr + 2);
      recycleEntry(curr, 2+currop->getCacheEntryLength());
    } // while
  } // for i
}

void MEDDLY::monolithic_chained::show(FILE *s, int verbLevel) 
{
  if (verbLevel < 1) return;
  fprintf(s, "Monolithic compute table\n");
  fprintf(s, "%*sCurrent CT memory :\t%lu bytes\n", 6, "", currMemory);
  fprintf(s, "%*sPeak    CT memory :\t%lu bytes\n", 6, "", peakMemory);
  // verbLevel--;
  report(s, 6, verbLevel);
  verbLevel--;
  dumpInternal(s, verbLevel);
}


int MEDDLY::monolithic_chained::convertToList(bool removeStales)
{
  int list = 0;
  for (unsigned i=0; i<tableSize; i++) {
    while (table[i]) {
      int curr = table[i];
      table[i] = entries[curr];
      if (removeStales) {
        operation* currop = operation::getOpWithIndex(entries[curr+1]);
        MEDDLY_DCASSERT(currop);
        const int* entry = entries + curr + 2;
        //
        // Check for stale
        //
        if (currop->isEntryStale(entry)) {
#ifdef DEBUG_TABLE2LIST
          printf("\tstale ");
          currop->showEntry(stdout, entry);
          printf(" (handle %d slot %d)\n", curr, i);
#endif
          currop->discardEntry(entry);
          recycleEntry(curr, 2+currop->getCacheEntryLength());
          continue;
        }
      } // if removeStales
      //
      // Not stale, move to list
      //
#ifdef DEBUG_TABLE2LIST
      printf("\tkeep  ");
      currop->showEntry(stdout, entry);
      printf(" (handle %d slot %d)\n", curr, i);
#endif
      entries[curr] = list;
      list = curr;
    } // while
  } // for i
#ifdef DEBUG_TABLE2LIST
  printf("Built list: %d", list);
  if (list ) for (int L=entries[list]; L; L=entries[L])
    printf("->%d", L);
  printf("\n");
#endif
  return list;
}

void MEDDLY::monolithic_chained::listToTable(int L)
{
#ifdef DEBUG_LIST2TABLE
  printf("Recovering  list: %d", L);
  if (L) for (int i=entries[L]; i; i=entries[i])
    printf("->%d", i);
  printf("\n");
#endif
  while (L) {
    int curr = L;
    L = entries[L];
    operation* currop = operation::getOpWithIndex(entries[curr+1]);
    MEDDLY_DCASSERT(currop);
    int hashlength = 1+currop->getKeyLength();
    unsigned h = hash(entries + curr + 1, hashlength);
    entries[curr] = table[h];
    table[h] = curr;
#ifdef DEBUG_LIST2TABLE
    printf("\tsave  ");
    currop->showEntry(stdout, entries + curr + 2);
    printf(" (handle %d slot %d)\n", curr, h);
#endif
  }
}

void MEDDLY::monolithic_chained::showEntry(FILE *s, int curr) const
{ 
  operation* op = operation::getOpWithIndex(entries[curr+1]);
  op->showEntry(s, entries + curr + 2);
}


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                     operation_chained  methods                     *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::operation_chained
::operation_chained(const settings::computeTableSettings &s, operation* op)
 : base_chained(s)
{
  global_op = op;
}

MEDDLY::operation_chained::~operation_chained()
{
}

MEDDLY::compute_table::search_key* 
MEDDLY::operation_chained::initializeSearchKey(operation* op)
{
  if (op != global_op)
    throw error(error::UNKNOWN_OPERATION);
  return init(op, 0);
}

MEDDLY::compute_table::entry_builder& 
MEDDLY::operation_chained::startNewEntry(search_key *key)
{
  MEDDLY_DCASSERT(key);
  if (key->getOp() != global_op)
    throw error(error::UNKNOWN_OPERATION);
  startIndexedEntry(smart_cast<old_search_key*>(key), 1, 0);
  return currEntry;
}

void MEDDLY::operation_chained::removeAll()
{
  for (unsigned i=0; i<tableSize; i++) {
    while (table[i]) {
      int curr = table[i];
      table[i] = entries[curr];
      global_op->discardEntry(entries + curr + 1);
      recycleEntry(curr, 1+global_op->getCacheEntryLength());
    } // while
  } // for i
}

void MEDDLY::operation_chained::show(FILE *s, int verbLevel)
{
  if (verbLevel < 1) return;
  fprintf(s, "Compute table for %s (index %d)\n", 
    global_op->getName(), global_op->getIndex()
  );
  fprintf(s, "%*sCurrent CT memory :\t%lu bytes\n", 6, "", currMemory);
  fprintf(s, "%*sPeak    CT memory :\t%lu bytes\n", 6, "", peakMemory);
  // verbLevel--;
  report(s, 6, verbLevel);
  verbLevel--;
  dumpInternal(s, verbLevel);
}


int MEDDLY::operation_chained::convertToList(bool removeStales)
{
  int list = 0;
  for (unsigned i=0; i<tableSize; i++) {
    while (table[i]) {
      int curr = table[i];
      table[i] = entries[curr];
      if (removeStales) {
        const int* entry = entries + curr + 1;
        //
        // Check for stale
        //
        if (global_op->isEntryStale(entry)) {
          global_op->discardEntry(entry);
          recycleEntry(curr, 1+global_op->getCacheEntryLength());
          continue;
        }
      }
      //
      // Not stale, move to list
      //
      entries[curr] = list;
      list = curr;
    } // while
  } // for i
  return list;
}

void MEDDLY::operation_chained::listToTable(int L)
{
  while (L) {
    int curr = L;
    L = entries[L];
    int hashlength = global_op->getKeyLength();
    int h = hash(entries + curr + 1, hashlength);
    entries[curr] = table[h];
    table[h] = curr;
  }
}

// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                   operation_chained_fast methods                   *
// *                                                                    *
// *                                                                    *
// **********************************************************************


template <int N>
MEDDLY::compute_table::search_result& 
MEDDLY::operation_chained_fast<N>::find(search_key *k)
{
  static old_search_result ANS;
  old_search_key* key = smart_cast <old_search_key*>(k);
  MEDDLY_DCASSERT(key);
  perf.pings++;
  unsigned h = hash(key);
  int prev = 0;
  int curr = table[h];
  int chain = 0;
  ANS.setInvalid();
  while (curr) {
    chain++;
    //
    // Check for match
    //
    if (equal_sw(entries+curr+1, key->rawData(), N)) {
      if (global_op->shouldStaleCacheHitsBeDiscarded()) {
        if (checkStale(h, prev, curr)) {
          // The match is stale.
          // Since there can NEVER be more than one match
          // in the table, we're done!
          break;
        }
      } 
      // "Hit"
      perf.hits++;
      if (prev) {
        // not at the front; move it there
        entries[prev] = entries[curr];
        entries[curr] = table[h];
        table[h] = curr;
      }
#ifdef DEBUG_CT
      printf("Found CT entry ");
      global_op->showEntry(stdout, entries + curr + 1);
      printf("\n");
#endif
      ANS.setResult(entries+curr+1+key->dataLength(), global_op->getAnsLength());
      break;
    };
    //
    // No match; maybe check stale
    //
    if (checkStalesOnFind) {
      if (checkStale(h, prev, curr)) continue;
    }
    // advance pointers
    prev = curr;
    curr = entries[curr];
  }
  sawSearch(chain);
  return ANS;
}

// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                       base_unchained methods                       *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::base_unchained
::base_unchained(const settings::computeTableSettings &s)
 : base_hash(s, 1024, 512)
{
  collisions = 0;
  global_op = 0;
}

MEDDLY::base_unchained::~base_unchained()
{
}

void MEDDLY::base_unchained::show(FILE *s, int verbLevel) 
{
  if (verbLevel < 1) return;
  showTitle(s);
  fprintf(s, "%*sCurrent CT memory :\t%lu bytes\n", 6, "", currMemory);
  fprintf(s, "%*sPeak    CT memory :\t%lu bytes\n", 6, "", peakMemory);
  verbLevel--;
  if (verbLevel < 1) return;
  fprintf(s, "%*sCollisions        :\t%ld\n", 6, "", collisions);
  report(s, 6, verbLevel);
  verbLevel--;
  if (verbLevel < 1) return;

  fprintf(s, "\nHash table:\n");
  
  for (unsigned i=0; i<tableSize; i++) {
    int curr = table[i];
    if (0==curr) continue;
    fprintf(s, "\t%9u:  node %9d: ", i, curr);
    showEntry(s, curr);
    fprintf(s, "\n");
  }

  base_hash::dumpInternal(s, verbLevel-1);
}


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                    monolithic_unchained methods                    *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::monolithic_unchained
::monolithic_unchained(const settings::computeTableSettings &s)
 : base_unchained(s)
{
}

MEDDLY::monolithic_unchained::~monolithic_unchained()
{
}

MEDDLY::compute_table::search_key* 
MEDDLY::monolithic_unchained::initializeSearchKey(operation* op)
{
  return init(op, 1);
}


MEDDLY::compute_table::search_result&
MEDDLY::monolithic_unchained::find(search_key *k)
{
  static old_search_result ANS;
  old_search_key* key = smart_cast <old_search_key*>(k);
  MEDDLY_DCASSERT(key);
  perf.pings++;
  unsigned h = hash(key);
  unsigned hcurr = h;
  int chain;
  ANS.setInvalid();
  for (chain=0; chain<=maxCollisionSearch; chain++, incMod(hcurr) ) {
    int curr = table[hcurr];
    if (0==curr) continue;
    //
    // Check for match
    //
    if (equal_sw(entries+curr, key->rawData(), key->dataLength())) {
      if (key->getOp()->shouldStaleCacheHitsBeDiscarded()) {
        if (checkStale<1>(hcurr, curr)) {
          // The match is stale.
          // Since there can NEVER be more than one match
          // in the table, we're done!
          break;
        }
      }
      // "Hit"
      perf.hits++;
#ifdef DEBUG_CT
      printf("Found CT entry ");
      key->getOp()->showEntry(stdout, entries + curr + 1);
      // fprintf(stderr, " in slot %u", h);
      printf("\n");
#endif
      ANS.setResult(entries+curr+key->dataLength(), key->getOp()->getAnsLength());
      break;
    };
    //
    // No match; maybe check stale
    //
    if (checkStalesOnFind) {
      checkStale<1>(hcurr, curr);
    }
  } // for chain
  sawSearch(chain);
  return ANS;
}

MEDDLY::compute_table::entry_builder& 
MEDDLY::monolithic_unchained::startNewEntry(search_key *key)
{
  MEDDLY_DCASSERT(key);
  startIndexedEntry(smart_cast<old_search_key*>(key), 0, 1);
  return currEntry;
}

void MEDDLY::monolithic_unchained::addEntry()
{
  addEntryT<1>();
}

void MEDDLY::monolithic_unchained::removeStales()
{
  removeStalesT<1>();
}
   
void MEDDLY::monolithic_unchained::removeAll()
{
  for (unsigned i=0; i<tableSize; i++) {
    if (0==table[i]) continue;
    remove<1>(table[i]);
    table[i] = 0;
  }
}

void MEDDLY::monolithic_unchained::showTitle(FILE* s) const
{
  fprintf(s, "Monolithic compute table\n");
}

void MEDDLY::monolithic_unchained::showEntry(FILE *s, int curr) const 
{ 
  operation* op = operation::getOpWithIndex(entries[curr]);
  op->showEntry(s, entries + curr + 1);
}

// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                    operation_unchained  methods                    *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::operation_unchained
::operation_unchained(const settings::computeTableSettings &s, operation* op)
 : base_unchained(s)
{
  global_op = op;
}

MEDDLY::operation_unchained::~operation_unchained()
{
}

MEDDLY::compute_table::search_key* 
MEDDLY::operation_unchained::initializeSearchKey(operation* op)
{
  if (op != global_op)
    throw error(error::UNKNOWN_OPERATION);
  return init(op, 0);
}

MEDDLY::compute_table::entry_builder& 
MEDDLY::operation_unchained::startNewEntry(search_key *key)
{
  MEDDLY_DCASSERT(key);
  startIndexedEntry(smart_cast<old_search_key*>(key), 0, 0);
  return currEntry;
}

void MEDDLY::operation_unchained::addEntry()
{
  addEntryT<0>();
}

void MEDDLY::operation_unchained::removeStales()
{
  removeStalesT<0>();
}
   
void MEDDLY::operation_unchained::removeAll()
{
  for (unsigned i=0; i<tableSize; i++) {
    if (0==table[i]) continue;
    remove<0>(table[i]);
    table[i] = 0;
  }
}

void MEDDLY::operation_unchained::showTitle(FILE* s) const
{
  fprintf(s, "Compute table for %s (index %d)\n", 
    global_op->getName(), global_op->getIndex()
  );
}

void MEDDLY::operation_unchained::showEntry(FILE *s, int curr) const 
{ 
  global_op->showEntry(s, entries + curr);
}


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                  operation_unchained_fast methods                  *
// *                                                                    *
// *                                                                    *
// **********************************************************************


template <int N>
MEDDLY::compute_table::search_result&
MEDDLY::operation_unchained_fast<N>::find(search_key *k)
{
  static old_search_result ANS;
  old_search_key* key = smart_cast <old_search_key*>(k);
  MEDDLY_DCASSERT(key);
  perf.pings++;
  unsigned h = hash(key);
  unsigned hcurr = h;
  int chain;
  ANS.setInvalid();
  for (chain=0; chain<=maxCollisionSearch; chain++, incMod(hcurr) ) {
    int curr = table[hcurr];
    if (0==curr) continue;
    //
    // Check for match
    //
    if (equal_sw(entries+curr, key->rawData(), N)) {
      if (key->getOp()->shouldStaleCacheHitsBeDiscarded()) {
        if (checkStale<0>(hcurr, curr)) {
          // The match is stale.
          // Since there can NEVER be more than one match
          // in the table, we're done!
          break;
        }
      }
      // "Hit"
      perf.hits++;
#ifdef DEBUG_CT
      printf("Found CT entry ");
      global_op->showEntry(stdout, entries + curr);
      // fprintf(stderr, " in slot %u", h);
      printf("\n");
#endif
      ANS.setResult(entries+curr+key->dataLength(), global_op->getAnsLength());
      break;
    };
    //
    // No match; maybe check stale
    //
    if (checkStalesOnFind) {
      checkStale<0>(hcurr, curr);
    }
  } // for chain
  sawSearch(chain);
  return ANS;
}

// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                          base_map methods                          *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::base_map
::base_map(const settings::computeTableSettings &s, operation* op)
 : base_table(s)
{
  global_op = op;
  current = 0;
}

MEDDLY::base_map::~base_map()
{
}

MEDDLY::compute_table::search_key* 
MEDDLY::base_map::initializeSearchKey(operation* op)
{
  if (op != global_op)
    throw error(error::UNKNOWN_OPERATION);
  return init(op, 0);
}

MEDDLY::compute_table::entry_builder& 
MEDDLY::base_map::startNewEntry(search_key *key)
{
  MEDDLY_DCASSERT(key);
  if (key->getOp() != global_op)
    throw error(error::UNKNOWN_OPERATION);
  current = startPtrEntry(smart_cast<old_search_key*>(key), 0, 0);
  return currEntry;
}


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                       operation_map  methods                       *
// *                                                                    *
// *                                                                    *
// **********************************************************************

template <int K>
MEDDLY::operation_map<K>
::operation_map(const settings::computeTableSettings &s, operation* op)
 : base_map(s, op)
{
}

template <int K>
MEDDLY::operation_map<K>
::~operation_map()
{
}


template <int K>
MEDDLY::compute_table::search_result&
MEDDLY::operation_map<K>::find(search_key *k) 
{
  static old_search_result ANS;
  old_search_key* key = smart_cast <old_search_key*>(k);
  MEDDLY_DCASSERT(key);
  perf.pings++;
  typename std::map <int*, int*, less<K> >::iterator 
    ans = ct.find(key->rawData());
      
  ANS.setInvalid();

  if (ans == ct.end()) {
    return ANS;
  }
  int* h = ans->second;
  if (global_op->shouldStaleCacheHitsBeDiscarded()) {
    if (isStale(h)) {
      ct.erase(ans);
      removeEntry(h);
      return ANS;
    }
  }
#ifdef DEBUG_CT
  printf("Found CT entry ");
  global_op->showEntry(stdout, h);
  printf("\n");
#endif
  perf.hits++;
  ANS.setResult(h+K, global_op->getAnsLength());
  return ANS;
}


template <int K>
void MEDDLY::operation_map<K>::addEntry()
{
#ifdef DEBUG_CT
  printf("Adding CT entry ");
  showEntry(stdout, current);
  printf("\n");
#endif
#ifdef DEVELOPMENT_CODE
  assert(
#endif
  (ct.insert(std::make_pair(current, current)))
#ifdef DEVELOPMENT_CODE
    .second
  )
#endif
  ;
  current = 0;
}

template <int K>
void MEDDLY::operation_map<K>::removeStales()
{
#ifdef DEBUG_SLOW
  fprintf(stdout, "Removing stales in CT (entries %ld)\n", perf.numEntries);
#endif
#ifdef SUMMARY_STALES
  int stales = 0;
#endif

  typename std::map<int*, int*, less<K> >::iterator curr = ct.begin();
  typename std::map<int*, int*, less<K> >::iterator end = ct.end();
  while (curr != end) {
    int* h = curr->second;
    if (isStale(h)) {
      ct.erase(curr++);
      removeEntry(h);
#ifdef SUMMARY_STALES
      stales++;
#endif
    } else {
      ++curr;
    }
  }

#ifdef DEBUG_SLOW
  fprintf(stdout, "Done removing CT stales (entries %ld)\n", perf.numEntries);
#endif
#ifdef SUMMARY_STALES
  printf("CT %s (index %d) removed %d stales\n", 
    global_op->getName(), global_op->getIndex(), stales);
#endif
}

template <int K>
void MEDDLY::operation_map<K>::removeAll()
{
  typename std::map<int*, int*, less<K> >::iterator curr = ct.begin();
  typename std::map<int*, int*, less<K> >::iterator end = ct.end();
  while (curr != end) {
    int* h = curr->second;
    ct.erase(curr++);
    removeEntry(h);
  }
  MEDDLY_DCASSERT(ct.empty());
}

template <int K>
void MEDDLY::operation_map<K>::show(FILE *s, int verbLevel) 
{
  if (verbLevel < 1) return;
  fprintf(s, "Compute table for %s (index %d)\n", 
    global_op->getName(), global_op->getIndex()
  );

  fprintf(s, "\tMap size: %ld\n", long(ct.size()));
  verbLevel -= 4;
  if (verbLevel<1) return;

  fprintf(s, "Map entries, in order:\n\t");
  typename std::map<int*, int*, less<K> >::iterator curr = ct.begin();
  typename std::map<int*, int*, less<K> >::iterator end = ct.end();
  const char* comma = "";
  while (curr != end) {
    showEntry(s, curr->second);
    fputs(comma, s);
    comma = ", ";
    ++curr;
  }
  fprintf(s, "\n\n");

  verbLevel--;
  dumpInternal(s, verbLevel);
}



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                           Hidden  styles                           *
// *                                                                    *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {
  class monolithic_chained_style;
  class monolithic_unchained_style;
  class operation_chained_style;
  class operation_unchained_style;
  class operation_map_style;
};

// **********************************************************************

class MEDDLY::monolithic_chained_style : public compute_table_style {
  public:
    monolithic_chained_style() { }
    virtual ~monolithic_chained_style() { }

    virtual compute_table* create(const settings::computeTableSettings &s)
      const 
    {
      return new monolithic_chained(s);
    }

    virtual bool usesMonolithic() const {
      return true;
    }
};

// **********************************************************************

class MEDDLY::monolithic_unchained_style : public compute_table_style {
  public:
    monolithic_unchained_style() { }
    virtual ~monolithic_unchained_style() { }

    virtual compute_table* create(const settings::computeTableSettings &s)
      const 
    {
      return new monolithic_unchained(s);
    }

    virtual bool usesMonolithic() const {
      return true;
    }
};

// **********************************************************************

class MEDDLY::operation_chained_style : public compute_table_style {
  public:
    operation_chained_style() { }
    virtual ~operation_chained_style() { }

    virtual compute_table* create(const settings::computeTableSettings &s,
      operation* op) const 
    {
      switch (op->getKeyLength()) {
        case 8:   return new operation_chained_fast<8>(s, op);
        case 7:   return new operation_chained_fast<7>(s, op);
        case 6:   return new operation_chained_fast<6>(s, op);
        case 5:   return new operation_chained_fast<5>(s, op);
        case 4:   return new operation_chained_fast<4>(s, op);
        case 3:   return new operation_chained_fast<3>(s, op);
        case 2:   return new operation_chained_fast<2>(s, op);
        case 1:   return new operation_chained_fast<1>(s, op);
        default:  assert(0);
                  return 0;
      }
      assert(0);
      return 0;
    }

    virtual bool usesMonolithic() const {
      return false;
    }
};


// **********************************************************************

class MEDDLY::operation_unchained_style : public compute_table_style {
  public:
    operation_unchained_style() { }
    virtual ~operation_unchained_style() { }

    virtual compute_table* create(const settings::computeTableSettings &s,
      operation* op) const 
    {
      switch (op->getKeyLength()) {
        case 8:   return new operation_unchained_fast<8>(s, op);
        case 7:   return new operation_unchained_fast<7>(s, op);
        case 6:   return new operation_unchained_fast<6>(s, op);
        case 5:   return new operation_unchained_fast<5>(s, op);
        case 4:   return new operation_unchained_fast<4>(s, op);
        case 3:   return new operation_unchained_fast<3>(s, op);
        case 2:   return new operation_unchained_fast<2>(s, op);
        case 1:   return new operation_unchained_fast<1>(s, op);
        default:  assert(0);
                  return 0;
      }
      assert(0);
      return 0;
    }

    virtual bool usesMonolithic() const {
      return false;
    }
};


// **********************************************************************

class MEDDLY::operation_map_style : public compute_table_style {
  public:
    operation_map_style() { }
    virtual ~operation_map_style() { }

    virtual compute_table* create(const settings::computeTableSettings &s,
      operation* op) const 
    {
      switch (op->getKeyLength()) {
        case 8:   return new operation_map<8>(s, op);
        case 7:   return new operation_map<7>(s, op);
        case 6:   return new operation_map<6>(s, op);
        case 5:   return new operation_map<5>(s, op);
        case 4:   return new operation_map<4>(s, op);
        case 3:   return new operation_map<3>(s, op);
        case 2:   return new operation_map<2>(s, op);
        case 1:   return new operation_map<1>(s, op);
        default:  assert(0);
                  return 0;
      }
      assert(0);
      return 0;
    }

    virtual bool usesMonolithic() const {
      return false;
    }
};


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                             Front  End                             *
// *                                                                    *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {
  monolithic_chained_style THE_MCS;
  monolithic_unchained_style THE_MUS;
  operation_chained_style THE_OCS;
  operation_unchained_style THE_OUS;
  operation_map_style THE_OMS;

  const compute_table_style* MonolithicChainedHash = &THE_MCS;
  const compute_table_style* MonolithicUnchainedHash = &THE_MUS;
  const compute_table_style* OperationChainedHash = &THE_OCS;
  const compute_table_style* OperationUnchainedHash = &THE_OUS;
  const compute_table_style* OperationMap = &THE_OMS;
};



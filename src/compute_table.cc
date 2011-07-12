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


#include "compute_table.h"

#include <map>  // for operation_map

// #define DEBUG_CT
// #define DEBUG_TABLE2LIST
// #define DEBUG_LIST2TABLE

namespace MEDDLY {
  /// base class for all compute tables; 
  /// handles allocation of entries :^)
  class base_table;

  /// base class for hash tables (gives hash function)
  class base_hash;

  class base_chained;

  class monolithic_chained;
  class operation_chained;

  // const float expansionFactor = 1.5;


  // extern settings meddlySettings;
}

// **********************************************************************
// *                                                                    *
// *                       compute_table  methods                       *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table::compute_table(const settings::computeTableSettings &s)
{
  maxSize = s.maxSize;
  if (0==maxSize)
    throw error(error::INVALID_ASSIGNMENT);

  switch (s.staleRemoval) {
    case settings::computeTableSettings::Aggressive:
            checkStalesOnFind = true;
            checkStalesOnResize = true;
            break;
    case settings::computeTableSettings::Moderate:
            checkStalesOnFind = false;
            checkStalesOnResize = true;
            break;
    case settings::computeTableSettings::Lazy:
            checkStalesOnFind = false;
            checkStalesOnResize = false;
            break;
  }

  perf.numEntries = 0;
  perf.hits = 0;
  perf.pings = 0;
  perf.numLargeSearches = 0;
  perf.maxSearchLength = 0;
  for (int i=0; i<perf.searchHistogramSize; i++)
    perf.searchHistogram[i] = 0;
}

MEDDLY::compute_table::~compute_table()
{
}

MEDDLY::compute_table::search_key::search_key()
{
  op = 0;
  hashLength = 0;
  data = 0;
  key_data = 0;
}

MEDDLY::compute_table::search_key::~search_key()
{
  delete[] data;
}

// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                          base_table class                          *
// *                                                                    *
// *                                                                    *
// **********************************************************************

class MEDDLY::base_table : public compute_table {
  public:
    base_table(const settings::computeTableSettings &s);
    virtual ~base_table();

  protected:
    inline void sawSearch(int c) {
      if (c>=stats::searchHistogramSize) {
        perf.numLargeSearches++;
      } else {
        perf.searchHistogram[c]++;
      }
      if (c>perf.maxSearchLength) perf.maxSearchLength = c;
    }
    int newEntry(int size);
    inline void recycleEntry(int h, int size) {
      entries[h] = freeList[size];
      freeList[size] = h;
      perf.numEntries--;
    }
    void dumpInternal(FILE* s, int verbLevel) const;
    void report(FILE* s, int indent, int &level, unsigned long mem) const;

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

    // M is 1 if we need a slot for the operation index, 0 otherwise.
    static inline void init(search_key &key, operation* op, int M) {
      DCASSERT(0==key.data);
      key.hashLength = op->getKeyLength() + M;
      key.data = new int[key.hashLength];
      key.key_data = key.data+M;
      if (M) key.data[0] = op->getIndex();
      key.op = op;
#ifdef DEVELOPMENT_CODE
      key.keyLength = op->getKeyLength();
#endif
    }
    
    // C is the number of slots we need for "chaining".
    // M is 1 if we need a slot for the operation index, 0 otherwise.
    inline void startNewEntry(operation* op, int C, int M) {
      currEntry.handle = newEntry(op->getCacheEntryLength()+C+M);
      currEntry.hashLength = op->getKeyLength()+M;
      currEntry.entry = entries + currEntry.handle;
      if (M) currEntry.entry[C] = op->getIndex();
      currEntry.key_entry = currEntry.entry + C + M;
      currEntry.res_entry = currEntry.key_entry + op->getKeyLength();
#ifdef DEVELOPMENT_CODE
      currEntry.keyLength = op->getKeyLength();
      currEntry.resLength = op->getAnsLength();
#endif
    }
  protected:
    int*  entries;
    int entriesSize;
    int entriesAlloc;
  private:
    static const int maxEntrySize = 15;
    static const int maxEntryBytes = sizeof(int) * maxEntrySize;
    int* freeList;
};

// **********************************************************************
// *                         base_table methods                         *
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
}

MEDDLY::base_table::~base_table()
{
  free(entries);
  delete[] freeList;
}

//
// TBD: old method was to expand by a factor of 1.5
//
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
    return h;
  }
  if (entriesSize + size > entriesAlloc) {
    // need to expand
    int neA = MIN(entriesAlloc * 2, entriesAlloc + 4096);
    int* ne = (int*) realloc(entries, neA * sizeof(int));
    if (0==ne) throw error(error::INSUFFICIENT_MEMORY);
    entries = ne;
    entriesAlloc = neA;
  }
  DCASSERT(entriesSize + size <= entriesAlloc);
  int h = entriesSize;
  entriesSize += size;
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
::report(FILE* s, int indent, int &level, unsigned long mem) const
{
  if (level < 1) return;
  fprintf(s, "%*sNumber of entries :\t%ld\n", indent, "", perf.numEntries);
  fprintf(s, "%*sEntry array size  :\t%d\n", indent, "", entriesSize);
  fprintf(s, "%*sEntry array alloc :\t%d\n", indent, "", entriesAlloc);
  mem += sizeof(int) * entriesAlloc;
  fprintf(s, "%*sTotal Memory usage:\t%lu\n", indent, "", mem);

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
// *                          base_hash  class                          *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/** Abstract base class for hashing compute tables.
  * All use the same hash function:
  *
  *   Bob Jenkin's Hash
  *   Free to use for educational or commerical purposes
  *   http://burtleburtle.net/bob/hash/doobs.html
  */
class MEDDLY::base_hash : public base_table {
  public:
    base_hash(const settings::computeTableSettings &s, int initTS, int initTE);
    virtual ~base_hash();

  protected:
    inline unsigned hash(const int* k, int length) const {
      return raw_hash(k, length) % tableSize;
    }
    inline unsigned hash(const search_key &key) const {
      return hash(key.rawData(), key.dataLength());
    }
    
    void dumpInternal(FILE* s, int verbLevel) const;
    void report(FILE* s, int indent, int &level, unsigned long mem) const;
  protected:
    int*  table;
    unsigned int tableSize;
    unsigned int tableExpand;
  private:
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
    static unsigned raw_hash(const int* k, int length);
};

// **********************************************************************
// *                         base_hash  methods                         *
// **********************************************************************

MEDDLY::base_hash::base_hash(const settings::computeTableSettings &s, 
  int initTS, int initTE) : base_table(s)
{
  tableSize = initTS;
  tableExpand = initTE;
  table = (int*) malloc(tableSize * sizeof(int));
  if (0==table) throw error(error::INSUFFICIENT_MEMORY);
  for (unsigned i=0; i<tableSize; i++) table[i] = 0;
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
::report(FILE* s, int indent, int &level, unsigned long mem) const
{
  if (level < 1) return;
  fprintf(s, "%*sHash table size   :\t%d\n", indent, "", tableSize);
  mem += sizeof(int) * tableSize;
  base_table::report(s, indent, level, mem);
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
// *                         base_chained class                         *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/// Abstract base class for hash tables with chaining.
class MEDDLY::base_chained : public base_hash {
  public:
    base_chained(const settings::computeTableSettings &s);
    virtual ~base_chained();

    virtual void addEntry();
    virtual void removeStales();
  protected:
    void dumpInternal(FILE* s, int verbLevel) const;
    virtual int convertToList(bool removeStales) = 0;
    virtual void listToTable(int h) = 0;
    virtual void showEntry(FILE* s, int h) const = 0;
};

// **********************************************************************
// *                        base_chained methods                        *
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
  unsigned h = hash(currEntry.readEntry(1), currEntry.readLength());
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

  if (tableSize >= maxSize) {
    // Cannot expand
    listToTable(list);
    tableExpand = perf.numEntries + tableSize;
    if (tableExpand < 0) tableExpand = INT_MAX;
#ifdef DEBUG_SLOW
    fprintf(stdout, "CT already at maximum size\n");
#endif
    return;
  }

  unsigned newsize = tableSize*2;
  if (newsize > maxSize) newsize = maxSize;

  int* newt = (int*) realloc(table, newsize * sizeof(int));
  if (0==newt) throw error(error::INSUFFICIENT_MEMORY);

  for (unsigned i=tableSize; i<newsize; i++) newt[i] = 0;

  table = newt;
  tableSize = newsize;
  tableExpand = 4*tableSize;

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
  if ((tableSize > 1024) && (perf.numEntries * 2 < tableSize)) {
    // shrink table
    int newsize = tableSize / 2;
    int* newt = (int*) realloc(table, newsize * sizeof(int));
    if (0==newt) throw error(error::INSUFFICIENT_MEMORY); 
    table = newt;
    tableSize = newsize;
    tableExpand = 4*tableSize;
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
// *                      monolithic_chained class                      *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/*
    Anatomy of an entry:

      entries[h]      : next pointer in chain
      entries[h+1]    : operation index
      entries[h+2]    : first "payload" item
      ...
      entries[h+L+1]  : last "payload" item, L is cache entry length.
*/
class MEDDLY::monolithic_chained : public base_chained {
  public:
    monolithic_chained(const settings::computeTableSettings &s);
    virtual ~monolithic_chained();

    // required functions

    virtual bool isOperationTable() const   { return false; }
    virtual void initializeSearchKey(search_key &key, operation* op);
    virtual const int* find(const search_key &key);
    virtual temp_entry& startNewEntry(operation* op);
    virtual void removeAll();

    void show(FILE *s, int verbLevel) const;

  protected:
    virtual int convertToList(bool removeStales);
    virtual void listToTable(int h);
    virtual void showEntry(FILE* s, int h) const;

    inline bool checkStale(unsigned h, int prev, int &curr) {
        operation* currop = operation::getOpWithIndex(entries[curr+1]);
        DCASSERT(currop);
        if (currop->isEntryStale(entries+curr+2)) {
          currop->discardEntry(entries+curr+2);
          int next = entries[curr];
          if (prev) entries[prev] = next;
          else      table[h] = next;
          int length = currop->getCacheEntryLength();
          recycleEntry(curr, length+2);
          curr = next;
          return true;
        }
        return false;
    }
};

// **********************************************************************
// *                     monolithic_chained methods                     *
// **********************************************************************

MEDDLY::monolithic_chained
::monolithic_chained(const settings::computeTableSettings &s)
 : base_chained(s)
{
}

MEDDLY::monolithic_chained::~monolithic_chained()
{
}

void MEDDLY::monolithic_chained
::initializeSearchKey(search_key &key, operation* op)
{
  init(key, op, 1);
}

const int* MEDDLY::monolithic_chained::find(const search_key &key)
{
  perf.pings++;
  unsigned h = hash(key);
  int prev = 0;
  int curr = table[h];
  int chain = 0;
  while (curr) {
    chain++;
    //
    // Check for match
    //
    if (equal_sw(entries+curr+1, key.rawData(), key.dataLength())) {
      sawSearch(chain);
      if (key.getOp()->shouldStaleCacheHitsBeDiscarded()) {
        if (checkStale(h, prev, curr)) {
          // The match is stale.
          // Since there can NEVER be more than one match
          // in the table, we're done!
          return 0;
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
      key.getOp()->showEntry(stdout, entries + curr + 2);
      // fprintf(stderr, " in slot %u", h);
      printf("\n");
#endif
      return entries + curr + 2;
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
  return 0;
}

MEDDLY::compute_table::temp_entry& 
MEDDLY::monolithic_chained::startNewEntry(operation* op)
{
  base_table::startNewEntry(op, 1, 1);
  return currEntry;
}

void MEDDLY::monolithic_chained::removeAll()
{
  for (unsigned i=0; i<tableSize; i++) {
    while (table[i]) {
      int curr = table[i];
      table[i] = entries[curr];
      operation* currop = operation::getOpWithIndex(entries[curr+1]);
      DCASSERT(currop);
      currop->discardEntry(entries + curr + 2);
      recycleEntry(curr, 2+currop->getCacheEntryLength());
    } // while
  } // for i
}

void MEDDLY::monolithic_chained::show(FILE *s, int verbLevel) const
{
  if (verbLevel < 1) return;
  fprintf(s, "Monolithic compute table\n");
  report(s, 6, verbLevel, 0);
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
        DCASSERT(currop);
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
    DCASSERT(currop);
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
// *                      operation_chained  class                      *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/*
    Anatomy of an entry:

      entries[h]      : next pointer in chain
      entries[h+1]    : first "payload" item
      ...
      entries[h+L]    : last "payload" item, L is cache entry length.
*/
class MEDDLY::operation_chained : public base_chained {
  public:
    operation_chained(const settings::computeTableSettings &s, operation* op);
    virtual ~operation_chained();

    // required functions

    virtual bool isOperationTable() const   { return true; }
    virtual void initializeSearchKey(search_key &key, operation* op);
    virtual temp_entry& startNewEntry(operation* op);
    virtual void removeAll();

    void show(FILE *s, int verbLevel) const;

  protected:
    virtual int convertToList(bool removeStales);
    virtual void listToTable(int h);
    virtual void showEntry(FILE* s, int h) const { 
      global_op->showEntry(s, entries + h + 1);
    }

    inline bool checkStale(unsigned h, int prev, int &curr) {
        if (global_op->isEntryStale(entries+curr+1)) {
          global_op->discardEntry(entries+curr+1);
          int next = entries[curr];
          if (prev) entries[prev] = next;
          else      table[h] = next;
          int length = global_op->getCacheEntryLength();
          recycleEntry(curr, length+1);
          curr = next;
          return true;
        }
        return false;
    }
  protected:
    operation* global_op;
};

// **********************************************************************
// *                     operation_chained  methods                     *
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

void MEDDLY::operation_chained
::initializeSearchKey(search_key &key, operation* op)
{
  if (op != global_op)
    throw error(error::UNKNOWN_OPERATION);
  init(key, op, 0);
}

MEDDLY::compute_table::temp_entry& 
MEDDLY::operation_chained::startNewEntry(operation* op)
{
  if (op != global_op)
    throw error(error::UNKNOWN_OPERATION);
  base_table::startNewEntry(op, 1, 0);
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

void MEDDLY::operation_chained::show(FILE *s, int verbLevel) const
{
  if (verbLevel < 1) return;
  fprintf(s, "Compute table for %s (index %d)\n", 
    global_op->getName(), global_op->getIndex()
  );
  report(s, 6, verbLevel, 0);
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
// *                    operation_chained_fast class                    *
// *                                                                    *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {
  template <int N>
  class operation_chained_fast : public operation_chained {
    public:
      operation_chained_fast(const settings::computeTableSettings &s, 
        operation* op) : operation_chained(s, op) { }
      virtual ~operation_chained_fast() { }
      virtual const int* find(const search_key &key);
  };
};

// **********************************************************************
// *                     operation_chained  methods                     *
// **********************************************************************

template <int N>
const int* MEDDLY::operation_chained_fast<N>::find(const search_key &key)
{
  perf.pings++;
  unsigned h = hash(key);
  int prev = 0;
  int curr = table[h];
  int chain = 0;
  while (curr) {
    chain++;
    //
    // Check for match
    //
    if (equal_sw(entries+curr+1, key.rawData(), N)) {
      sawSearch(chain);
      if (global_op->shouldStaleCacheHitsBeDiscarded()) {
        if (checkStale(h, prev, curr)) {
          // The match is stale.
          // Since there can NEVER be more than one match
          // in the table, we're done!
          return 0;
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
      return entries + curr + 1;
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
  return 0;
}


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                        operation_map  class                        *
// *                                                                    *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {

/*
    Anatomy of an entry:

      entries[h  ]    : first "payload" item
      ...
      entries[h+L-1]  : last "payload" item, L is cache entry length.

    Pure payload :^)
    Don't worry, we pay for it in the "map" I'm sure.
*/
template <int N>
class operation_map : public base_table {
  protected:
    // for now ... key for our map.
    friend class key_type;

    class key_type {
        int handle;
        union {
          operation_map* parent;
          const int* search;
        };
      public:
        key_type() {
          handle = 0;
          search = 0;
        }
        key_type(operation_map* p, int h) {
          parent = p;
          handle = h;
        }
        void set(const int* s) {
          handle = 0;
          search = s;
        }
        inline const int* getData() const {
          if (handle) return parent->entries + handle;
          return search;
        }
    };

    // For comparing entries in our map.
    class key_less {
      public:
        inline bool operator() (const key_type &p, const key_type &q) const {
          const int* a = p.getData();
          const int* b = q.getData();
          switch (N) {
              case 8:   if (a[7] < b[7]) return true;
                        if (a[7] > b[7]) return false;
              case 7:   if (a[6] < b[6]) return true;
                        if (a[6] > b[6]) return false;
              case 6:   if (a[5] < b[5]) return true;
                        if (a[5] > b[5]) return false;
              case 5:   if (a[4] < b[4]) return true;
                        if (a[4] > b[4]) return false;
              case 4:   if (a[3] < b[3]) return true;
                        if (a[3] > b[3]) return false;
              case 3:   if (a[2] < b[2]) return true;
                        if (a[2] > b[2]) return false;
              case 2:   if (a[1] < b[1]) return true;
                        if (a[1] > b[1]) return false;
              case 1:   return a[0] < b[0];
              case 0:   return false;
              default:  return memcmp(a, b, N*sizeof(int)) < 0;
          }
        }
    };  

  public:
    operation_map(const settings::computeTableSettings &s, operation* op);
    virtual ~operation_map();

    // required functions

    virtual bool isOperationTable() const   { return true; }
    virtual void initializeSearchKey(search_key &key, operation* op);
    virtual const int* find(const search_key &key);

    virtual temp_entry& startNewEntry(operation* op);
    virtual void addEntry();
    virtual void removeStales();
    virtual void removeAll();

    virtual void show(FILE *s, int verbLevel = 0) const;
  protected:
    virtual void showEntry(FILE* s, int h) const {
      global_op->showEntry(s, entries + h);
    }
    inline bool isStale(const int* entry) {
      return global_op->isEntryStale(entry);
    }
    inline void removeEntry(int h) {
      global_op->discardEntry(entries+h);
      recycleEntry(h, global_op->getCacheEntryLength());
    }
  private:
    std::map<key_type, int, key_less> *ct;
    operation* global_op;
    key_type search;
};

};

// **********************************************************************
// *                       operation_map  methods                       *
// **********************************************************************

template <int N>
MEDDLY::operation_map<N>
::operation_map(const settings::computeTableSettings &s, operation* op)
 : base_table(s) 
{
  global_op = op;
  ct = new std::map<key_type, int, key_less>;
}

template <int N>
MEDDLY::operation_map<N>::~operation_map()
{
  delete ct;
}

template <int N>
void MEDDLY::operation_map<N>
::initializeSearchKey(search_key &key, operation* op)
{
  if (op != global_op)
    throw error(error::UNKNOWN_OPERATION);
  init(key, op, 0);
  search.set(key.rawData()); // really neat trick if it works...
}

template <int N>
const int* MEDDLY::operation_map<N>
::find(const search_key &key) 
{
  DCASSERT(ct);
  perf.pings++;
  typename std::map <key_type, int, key_less>::iterator 
    ans = ct->find(search);
      
  if (ans == ct->end()) return 0;
  int h = ans->second;
  if (global_op->shouldStaleCacheHitsBeDiscarded()) {
    if (isStale(entries+h)) {
      ct->erase(ans);
      removeEntry(h);
      return 0;
    }
  }
#ifdef DEBUG_CT
  printf("Found CT entry ");
  global_op->showEntry(stdout, entries+h);
  printf("\n");
#endif
  perf.hits++;
  return entries+h;
}

template <int N>
MEDDLY::compute_table::temp_entry& 
MEDDLY::operation_map<N>::startNewEntry(operation* op)
{
  if (op != global_op)
    throw error(error::UNKNOWN_OPERATION);
  base_table::startNewEntry(op, 0, 0);
  return currEntry;
}

template <int N>
void MEDDLY::operation_map<N>::addEntry()
{
  DCASSERT(ct);
#ifdef DEBUG_CT
  printf("Adding CT entry ");
  showEntry(stdout, currEntry.readHandle());
  printf("\n");
#endif
#ifdef DEVELOPMENT_CODE
  assert(
#endif
  (ct->insert(
    std::make_pair(
      key_type(this, currEntry.readHandle()),
      currEntry.readHandle()
    )
  ))
#ifdef DEVELOPMENT_CODE
  .second)
#endif
  ;
}

template <int N>
void MEDDLY::operation_map<N>::removeStales()
{
  if (0==ct) return;
#ifdef DEBUG_SLOW
  fprintf(stdout, "Removing stales in CT (entries %ld)\n", perf.numEntries);
#endif

  typename std::map<key_type, int, key_less>::iterator curr = ct->begin();
  typename std::map<key_type, int, key_less>::iterator end = ct->end();
  while (curr != end) {
    int h = curr->second;
    if (isStale(entries+h)) {
      ct->erase(curr++);
      removeEntry(h);
    } else {
      ++curr;
    }
  }

#ifdef DEBUG_SLOW
  fprintf(stdout, "Done removing CT stales (entries %ld)\n", perf.numEntries);
#endif
}

template <int N>
void MEDDLY::operation_map<N>::removeAll()
{
  if (0==ct) return;
  typename std::map<key_type, int, key_less>::iterator curr = ct->begin();
  typename std::map<key_type, int, key_less>::iterator end = ct->end();
  while (curr != end) {
    int h = curr->second;
    ct->erase(curr++);
    removeEntry(h);
  }
  DCASSERT(ct->empty());
}

template <int N>
void MEDDLY::operation_map<N>::show(FILE *s, int verbLevel) const
{
  if (verbLevel < 1) return;
  fprintf(s, "Compute table for %s (index %d)\n", 
    global_op->getName(), global_op->getIndex()
  );
  report(s, 6, verbLevel, 0);
  verbLevel--;
  if (verbLevel < 1) return;

  if (0==ct) {
    fprintf(s, "Empty map\n");
  } else {
    fprintf(s, "Map entries, in order:\n\t");
    typename std::map<key_type, int, key_less>::iterator curr = ct->begin();
    typename std::map<key_type, int, key_less>::iterator end = ct->end();
    const char* comma = "";
    while (curr != end) {
      fprintf(s, "%s%d", comma, curr->second);
      comma = ", ";
      ++curr;
    }
    fprintf(s, "\n\n");
  }

  verbLevel--;
  dumpInternal(s, verbLevel);
}



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                             Front  End                             *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table*
MEDDLY::createMonolithicTable(const settings::computeTableSettings  &s)
{
  switch (s.style) {

    case settings::computeTableSettings::MonolithicChainedHash:
        return new monolithic_chained(s);

    default:
        throw error(error::INVALID_ASSIGNMENT);
  }
}

MEDDLY::compute_table*
MEDDLY::createOperationTable(const settings::computeTableSettings  &s, 
  operation* op)
{
  switch (s.style) {

    case settings::computeTableSettings::OperationChainedHash:
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

    case settings::computeTableSettings::OperationMap:
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

    default:
        throw error(error::INVALID_ASSIGNMENT);
  }
}


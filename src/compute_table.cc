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

// #define DEBUG_CT
// #define DEBUG_TABLE2LIST
// #define DEBUG_LIST2TABLE

namespace MEDDLY {

  class base_table;

  class chained_table;

  class monolithic_table;
  class operation_table;
  const float expansionFactor = 1.5;

  extern settings meddlySettings;
}

// **********************************************************************
// *                                                                    *
// *                       compute_table  methods                       *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table::compute_table(const settings &s)
{
  chaining = s.doComputeTablesUseChaining;
  maxSize = s.maxComputeTableSize;

  if (0==maxSize)
    throw error(error::INVALID_ASSIGNMENT);

  perf.numEntries = 0;
  perf.hits = 0;
  perf.pings = 0;
  perf.numLargeChains = 0;
  perf.maxChainLength = 0;
  for (int i=0; i<perf.chainHistogramSize; i++)
    perf.chainHistogram[i] = 0;
}

MEDDLY::compute_table::~compute_table()
{
}

MEDDLY::compute_table::search_key::search_key()
{
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
    base_table(const settings &s);
    virtual ~base_table();

    virtual void addEntry();
    virtual void removeStales();

    virtual void show(FILE *s, int verbLevel = 0) const;

  protected:
    unsigned raw_hash(const int* k, int length); 

    int newEntry(int size);
    inline void recycleEntry(int h, int size) {
      entries[h] = freeList[size];
      freeList[size] = h;
      perf.numEntries--;
    }

    inline void sawChain(int c) {
      if (c>=stats::chainHistogramSize) {
        perf.numLargeChains++;
      } else {
        perf.chainHistogram[c]++;
      }
      if (c>perf.maxChainLength) perf.maxChainLength = c;
    }

    int*  table;
    unsigned int tableSize;
    unsigned int tableExpand;

    int*  entries;
    int entriesSize;
    int entriesAlloc;

    // Anatomy of an active compute table entry h:
    //  Note: C is 1 if we use chaining, 0 otherwise.
    //  Note: M is 1 if this is a monolithic table, 0 otherwise.
    //
    //  if (C),   entries[h]        : next pointer.
    //  if (M),   entries[h+C]      : operation index.
    //            entries[h+C+M]    : first item in entry.
    //            entries[h+C+M+1]  : second item in entry.
    //            ...
    //
    //  Anatomy of a recycled compute table entry h:
    //
    //            entries[h]        : next pointer in free list
    //

    void dumpInternal(FILE* s) const;

    virtual int convertToList() = 0;
    virtual void listToTable(int h) = 0;
    virtual void showTitle(FILE* s) const = 0;
    virtual void showEntry(FILE* s, int h) const = 0;

  private:
    static const int maxEntrySize = 15;
    int* freeList;
};

// **********************************************************************
// *                         base_table methods                         *
// **********************************************************************

MEDDLY::base_table::base_table(const settings &s)
 : compute_table(s)
{
  tableSize = 1024;
  tableExpand = 4*tableSize;
  table = (int*) malloc(tableSize * sizeof(int));
  if (0==table) throw error(error::INSUFFICIENT_MEMORY);
  for (unsigned i=0; i<tableSize; i++) table[i] = 0;
  
  entriesAlloc = 1024;
  entries = (int*) malloc(entriesAlloc * sizeof(int));
  entriesSize = 1;
  // entries[0] is never, ever, used.
  if (0==entries) throw error(error::INSUFFICIENT_MEMORY);

  freeList = new int[1+maxEntrySize];
  for (int i=0; i<=maxEntrySize; i++) freeList[i] = 0;
}

MEDDLY::base_table::~base_table()
{
  free(entries);
  free(table);
}

void MEDDLY::base_table::addEntry()
{
  unsigned h = raw_hash(currEntry.entry+1, currEntry.hashLength) % tableSize;
  currEntry.entry[0] = table[h];
  table[h] = currEntry.handle;

#ifdef DEBUG_CT
  printf("Adding CT entry ");
  showEntry(stdout, currEntry.handle);
  // fprintf(stderr, " to slot %u", h);
  printf("\n");
#endif

  if (perf.numEntries < tableExpand) return;

#ifdef DEBUG_SLOW
  fprintf(stdout, "Running GC in compute table (size %d, entries %ld)\n", tableSize, perf.numEntries);
#endif

  int list = convertToList();
  if (perf.numEntries < tableSize) {
    // Don't need to expand
    listToTable(list);
#ifdef DEBUG_SLOW
    fprintf(stdout, "Done CT GC, no resizing (now entries %ld)\n", perf.numEntries);
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

void MEDDLY::base_table::removeStales()
{
#ifdef DEBUG_SLOW
  fprintf(stdout, "Removing stales in CT (size %d, entries %ld)\n", tableSize, perf.numEntries);
#endif
  int list = convertToList();
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
  fprintf(stdout, "Done removing CT stales (size %d, entries %ld)\n", tableSize, perf.numEntries);
#endif
}

void MEDDLY::base_table::show(FILE *s, int verbLevel) const
{ 
  showTitle(s);
  char filler[] = "\t";
  fprintf(s, "%sNumber of entries :\t%ld\n", filler, perf.numEntries);
  fprintf(s, "%sHash table size   :\t%d\n", filler, tableSize);
  fprintf(s, "%sEntry array size  :\t%d\n", filler, entriesSize);
  fprintf(s, "%sEntry array alloc :\t%d\n", filler, entriesAlloc);
  unsigned long usage = 0;
  usage += sizeof(int) * tableSize;
  usage += sizeof(int) * entriesAlloc;
  fprintf(s, "%sTotal Memory usage:\t%lu\n", filler, usage);
  
  if (verbLevel < 1) return;

  fprintf(s, "%sPings             :\t%d\n", filler, perf.pings);
  fprintf(s, "%sHits              :\t%d\n", filler, perf.hits);

  if (verbLevel<2) return;

  fprintf(s, "%sChain length histogram:\n", filler);
  for (int i=0; i<stats::chainHistogramSize; i++) {
    if (perf.chainHistogram[i]) {
      fprintf(s, "%s%s%3d: %ld\n", filler, filler, i, perf.chainHistogram[i]);
    }
  }
  if (perf.numLargeChains)
    fprintf(s, "%sChains longer than %d: %ld\n", filler, 
            stats::chainHistogramSize-1, perf.numLargeChains
    );
  fprintf(s, "%sMax chain length seen: %d\n", filler, perf.maxChainLength);

  if (verbLevel<5) return;

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
  fprintf(s, "\n");

  if (verbLevel<6) return;

  fprintf(s, "Hash table nodes:\n");
  
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

  if (verbLevel<8) return;

  fprintf(s, "Raw hash table info:\n");
  dumpInternal(s);
}
/*
 * Bob Jenkin's Hash
 * Free to use for educational or commerical purposes
 * http://burtleburtle.net/bob/hash/doobs.html
 */
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

inline void mix(unsigned &a, unsigned &b, unsigned &c)
{
  a -= c;  a ^= rot(c, 4);  c += b; 
  b -= a;  b ^= rot(a, 6);  a += c;
  c -= b;  c ^= rot(b, 8);  b += a;
  a -= c;  a ^= rot(c,16);  c += b;
  b -= a;  b ^= rot(a,19);  a += c;
  c -= b;  c ^= rot(b, 4);  b += a;
}

inline void final(unsigned &a, unsigned &b, unsigned &c)
{
  c ^= b; c -= rot(b,14);
  a ^= c; a -= rot(c,11);
  b ^= a; b -= rot(a,25);
  c ^= b; c -= rot(b,16);
  a ^= c; a -= rot(c,4); 
  b ^= a; b -= rot(a,14);
  c ^= b; c -= rot(b,24);
}

unsigned MEDDLY::base_table::raw_hash(const int* k, int length)
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

  return c % tableSize;
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

void MEDDLY::base_table::dumpInternal(FILE* s) const
{
  if (0==table) fprintf(s, "Table: null\n");
  else {
    fprintf(s, "Table: [%d", table[0]);
    for (unsigned i=1; i<tableSize; i++) 
      fprintf(s, ", %d", table[i]);
    fprintf(s, "]\n");
  }
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


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                       monolithic_table class                       *
// *                                                                    *
// *                                                                    *
// **********************************************************************

class MEDDLY::monolithic_table : public base_table {
  public:
    monolithic_table(const settings &s);
    virtual ~monolithic_table();

    // required functions

    virtual bool isOperationTable() const   { return false; }
    virtual void initializeSearchKey(search_key &key, operation* op);
    virtual const int* find(const search_key &key);
    virtual temp_entry& startNewEntry(operation* op);
    virtual void removeAll();

  protected:
    virtual int convertToList();
    virtual void listToTable(int h);
    virtual void showTitle(FILE* s) const;
    virtual void showEntry(FILE* s, int h) const;
};

// **********************************************************************
// *                      monolithic_table methods                      *
// **********************************************************************

MEDDLY::monolithic_table::monolithic_table(const settings &s)
 : base_table(s)
{
}

MEDDLY::monolithic_table::~monolithic_table()
{
}

void MEDDLY::monolithic_table
::initializeSearchKey(search_key &key, operation* op)
{
  DCASSERT(0==key.data);
  key.hashLength = op->getKeyLength()+1;
  key.data = new int[key.hashLength];
  key.key_data = key.data+1;
  key.data[0] = op->getIndex();
#ifdef DEVELOPMENT_CODE
  key.keyLength = op->getKeyLength();
#endif
}

const int* MEDDLY::monolithic_table::find(const search_key &key)
{
  perf.pings++;
  unsigned h = raw_hash(key.data, key.hashLength);
  int prev = 0;
  int curr = table[h];
  int chain = 0;
  while (curr) {
    chain++;
    operation* currop = operation::getOpWithIndex(entries[curr+1]);
    DCASSERT(currop);
    //
    // Check for stale
    //
    if (currop->isMarkedForDeletion() || currop->isEntryStale(entries+curr+2)) {
      currop->discardEntry(entries+curr+2);
      int next = entries[curr];
      if (prev) entries[prev] = next;
      else      table[h] = next;
      int length = 1+currop->getCacheEntryLength();
      recycleEntry(curr, length+1);
      curr = next;
      continue;
    }
    //
    // Check for match
    //
    if (memcmp(entries+curr+1, key.data, key.hashLength*sizeof(int))==0) {
      // "Hit"
      perf.hits++;
      sawChain(chain);
      if (prev) {
        // not at the front; move it there
        entries[prev] = entries[curr];
        entries[curr] = table[h];
        table[h] = curr;
      }
#ifdef DEBUG_CT
      printf("Found CT entry ");
      currop->showEntry(stdout, entries + curr + 2);
      // fprintf(stderr, " in slot %u", h);
      printf("\n");
#endif
      return entries + curr + 2;
    };
    // advance pointers
    prev = curr;
    curr = entries[curr];
  }
  sawChain(chain);
  return 0;
}

MEDDLY::compute_table::temp_entry& 
MEDDLY::monolithic_table::startNewEntry(operation* op)
{
  currEntry.handle = newEntry(2+op->getCacheEntryLength());
  currEntry.hashLength = op->getKeyLength()+1;
  currEntry.entry = entries + currEntry.handle;
  currEntry.entry[1] = op->getIndex();
  currEntry.key_entry = currEntry.entry + 2;
  currEntry.res_entry = currEntry.key_entry + op->getKeyLength();
#ifdef DEVELOPMENT_CODE
  currEntry.keyLength = op->getKeyLength();
  currEntry.resLength = op->getAnsLength();
#endif
  return currEntry;
}

void MEDDLY::monolithic_table::removeAll()
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


int MEDDLY::monolithic_table::convertToList()
{
  int list = 0;
  for (unsigned i=0; i<tableSize; i++) {
    while (table[i]) {
      int curr = table[i];
      table[i] = entries[curr];
      operation* currop = operation::getOpWithIndex(entries[curr+1]);
      DCASSERT(currop);
      const int* entry = entries + curr + 2;
      //
      // Check for stale
      //
      if (currop->isMarkedForDeletion() || currop->isEntryStale(entry)) {
#ifdef DEBUG_TABLE2LIST
        printf("\tstale ");
        currop->showEntry(stdout, entry);
        printf(" (handle %d slot %d)\n", curr, i);
#endif
        currop->discardEntry(entry);
        recycleEntry(curr, 2+currop->getCacheEntryLength());
        continue;
      }
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

void MEDDLY::monolithic_table::listToTable(int L)
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
    unsigned h = raw_hash(entries + curr + 1, hashlength);
    entries[curr] = table[h];
    table[h] = curr;
#ifdef DEBUG_LIST2TABLE
    printf("\tsave  ");
    currop->showEntry(stdout, entries + curr + 2);
    printf(" (handle %d slot %d)\n", curr, h);
#endif
  }
}

void MEDDLY::monolithic_table::showTitle(FILE *s) const
{
  fprintf(s, "Monolithic compute table\n");
}

void MEDDLY::monolithic_table::showEntry(FILE *s, int curr) const
{ 
  operation* op = operation::getOpWithIndex(entries[curr+1]);
  op->showEntry(s, entries + curr + 2);
}


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                       operation_table  class                       *
// *                                                                    *
// *                                                                    *
// **********************************************************************

class MEDDLY::operation_table : public base_table {
  public:
    operation_table(const settings &s, operation* op);
    virtual ~operation_table();

    // required functions

    virtual bool isOperationTable() const   { return true; }
    virtual void initializeSearchKey(search_key &key, operation* op);
    virtual const int* find(const search_key &key);
    virtual temp_entry& startNewEntry(operation* op);
    virtual void removeAll();

  protected:
    virtual int convertToList();
    virtual void listToTable(int h);
    virtual void showTitle(FILE* s) const;
    virtual void showEntry(FILE* s, int h) const;
  private:
    operation* global_op;
};

// **********************************************************************
// *                      operation_table  methods                      *
// **********************************************************************

MEDDLY::operation_table::operation_table(const settings &s, operation* op)
 : base_table(s)
{
  global_op = op;
}

MEDDLY::operation_table::~operation_table()
{
}

void MEDDLY::operation_table
::initializeSearchKey(search_key &key, operation* op)
{
  if (op != global_op)
    throw error(error::UNKNOWN_OPERATION);
  DCASSERT(0==key.data);
  key.hashLength = op->getKeyLength();
  key.data = new int[key.hashLength];
  key.key_data = key.data;
#ifdef DEVELOPMENT_CODE
  key.keyLength = op->getKeyLength();
#endif
}

const int* MEDDLY::operation_table::find(const search_key &key)
{
  perf.pings++;
  unsigned h = raw_hash(key.data, key.hashLength);
  int prev = 0;
  int curr = table[h];
  int chain = 0;
  while (curr) {
    chain++;
    //
    // Check for stale
    //
    if (global_op->isMarkedForDeletion() 
          || global_op->isEntryStale(entries+curr+1)) 
    {
      global_op->discardEntry(entries+curr+1);
      int next = entries[curr];
      if (prev) entries[prev] = next;
      else      table[h] = next;
      int length = global_op->getCacheEntryLength();
      recycleEntry(curr, length+1);
      curr = next;
      continue;
    }
    //
    // Check for match
    //
    if (memcmp(entries+curr+1, key.data, key.hashLength*sizeof(int))==0) {
      // "Hit"
      perf.hits++;
      sawChain(chain);
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
    // advance pointers
    prev = curr;
    curr = entries[curr];
  }
  sawChain(chain);
  return 0;
}

MEDDLY::compute_table::temp_entry& 
MEDDLY::operation_table::startNewEntry(operation* op)
{
  if (op != global_op)
    throw error(error::UNKNOWN_OPERATION);
  currEntry.handle = newEntry(1+op->getCacheEntryLength());
  currEntry.hashLength = op->getKeyLength();
  currEntry.entry = entries + currEntry.handle;
  currEntry.key_entry = currEntry.entry + 1;
  currEntry.res_entry = currEntry.key_entry + op->getKeyLength();
#ifdef DEVELOPMENT_CODE
  currEntry.keyLength = op->getKeyLength();
  currEntry.resLength = op->getAnsLength();
#endif
  return currEntry;
}

void MEDDLY::operation_table::removeAll()
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

int MEDDLY::operation_table::convertToList()
{
  int list = 0;
  for (unsigned i=0; i<tableSize; i++) {
    while (table[i]) {
      int curr = table[i];
      table[i] = entries[curr];
      const int* entry = entries + curr + 1;
      //
      // Check for stale
      //
      if (global_op->isMarkedForDeletion() || global_op->isEntryStale(entry)) {
        global_op->discardEntry(entry);
        recycleEntry(curr, 1+global_op->getCacheEntryLength());
        continue;
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

void MEDDLY::operation_table::listToTable(int L)
{
  while (L) {
    int curr = L;
    L = entries[L];
    int hashlength = global_op->getKeyLength();
    int h = raw_hash(entries + curr + 1, hashlength);
    entries[curr] = table[h];
    table[h] = curr;
  }
}

void MEDDLY::operation_table::showTitle(FILE *s) const
{
  fprintf(s, "Compute table for %s (index %d)\n", 
    global_op->getName(), global_op->getIndex()
  );
}

void MEDDLY::operation_table::showEntry(FILE *s, int curr) const
{ 
  global_op->showEntry(s, entries + curr + 1);
}



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                             Front  End                             *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table*
MEDDLY::createMonolithicTable(const settings &s)
{
  return new monolithic_table(s);
}

MEDDLY::compute_table*
MEDDLY::createOperationTable(const settings &s, operation* op)
{
  return new operation_table(s, op);
}


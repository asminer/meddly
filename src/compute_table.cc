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
  fprintf(stderr, "Adding CT entry ");
  showEntry(stderr, currEntry.handle);
  // fprintf(stderr, " to slot %u", h);
  fprintf(stderr, "\n");
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
  char filler[] = "\t";
  fprintf(s, "Monolithic compute table\n");
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
    int hashLength = 1+currop->getKeyLength();
    if (memcmp(entries+curr+1, key.data, hashLength*sizeof(int))==0) {
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
      fprintf(stderr, "Found CT entry ");
      currop->showEntry(stderr, entries + curr + 2);
      // fprintf(stderr, " in slot %u", h);
      fprintf(stderr, "\n");
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
        currop->discardEntry(entry);
        recycleEntry(curr, 2+currop->getCacheEntryLength());
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

void MEDDLY::monolithic_table::listToTable(int L)
{
  while (L) {
    int curr = L;
    L = entries[L];
    operation* currop = operation::getOpWithIndex(entries[curr+1]);
    DCASSERT(currop);
    int hashlength = 1+currop->getCacheEntryLength();
    int h = raw_hash(entries + curr + 1, hashlength);
    entries[L] = table[h];
    table[h] = L;
  }
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
    int length = global_op->getCacheEntryLength();
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
      recycleEntry(curr, length+1);
      curr = next;
      continue;
    }
    //
    // Check for match
    //
    if (memcmp(entries+curr+1, key.data, length)==0) {
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
      fprintf(stderr, "Found CT entry ");
      global_op->showEntry(stderr, entries + curr + 1);
      fprintf(stderr, "\n");
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
    int hashlength = global_op->getCacheEntryLength();
    int h = raw_hash(entries + curr + 1, hashlength);
    entries[L] = table[h];
    table[h] = L;
  }
}

void MEDDLY::operation_table::showEntry(FILE *s, int curr) const
{ 
  global_op->showEntry(s, entries + curr + 1);
}



// **********************************************************************
// **********************************************************************
// **********************************************************************
#ifdef OLD_TABLE
// **********************************************************************
// **********************************************************************
// **********************************************************************



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                       monolithic_table class                       *
// *                                                                    *
// *                                                                    *
// **********************************************************************

class MEDDLY::monolithic_table : public compute_table {
  public:
    monolithic_table(settings s);
    virtual ~monolithic_table();

    // required functions

    virtual bool isOperationTable() const   { return false; }
    virtual void add(operation* op, const int* entry);
    virtual const int* find(operation* op, const int* entryKey);
    virtual void removeStales();
    virtual void removeAll();
    virtual void updateStats();

  private:

    // ******************************************************************
    // *                    Internal data-structures                    *
    // ******************************************************************

    /*  The compute table maintains an array of table entries which can
        be identified by {owner, data}.
    */
    struct table_entry {
      operation* op;                // operation to which entry belongs to
      int dataOffset;               // entry data (usually: {key, result})
      int next;                     // for hash table
    };

    /*  Cache entries that have been discard (or of no use anymore) are
        recycled and placed on a recycled_list corresponding to the size of the
        table entry (size of table entry is determined by the size of its data
        array).
    */
    struct recycled_list {
      int size;                     // size of nodes in this list
      int front;                    // first node in list, -1 terminates
      recycled_list* next;          // points to next recycled node list
    };

    // ******************************************************************
    // *                      Internal data                             *
    // ******************************************************************

    table_entry* nodes;             // Array of nodes
    int nodeCount;                  // size of array nodes
    int lastNode;                   // last used index in array nodes
    int* data;                      // Array of ints for storing data
    int dataCount;                  // size of array data
    int lastData;                   // last used index in array data
#ifdef RECYCLED_LIST
    int recycledNodes;              // -1 indicates no recycled nodes
    recycled_list* recycledFront;
#else
    std::vector<int> holes;
#endif
    fixed_size_hash_table<monolithic_table>* fsht;
#ifdef USE_CHAINED_HASH_TABLE
    chained_hash_table<monolithic_table>* ht;
#else
    hash_table<monolithic_table>* ht;
#endif

    // ******************************************************************
    // *                      Internal functions                        *
    // ******************************************************************

    // Expand the nodes array (grows by expansionFactor).
    void expandNodes();
    // Expand the data array (grows by expansionFactor).
    void expandData();

    // Return the address of the data associated with this entry
    inline int* getDataAddress(const table_entry& n) const {
      return data + n.dataOffset;
    }

    // Set the address of the data associated with this entry
    inline void setDataAddress(table_entry& n, const int* address) {
      n.dataOffset = address - data;
    }

    // Set the address of the data associated with this entry using offset
    inline void setDataOffset(table_entry& n, int offset) {
      n.dataOffset = offset;
    }

    // Is this a free node (i.e. stores no valid data)
    inline bool isFreeNode(const table_entry& n) const {
      return 0 == n.op;
    }

    // Is the nodes array full?
    inline bool isNodesStorageFull() const {
      return (lastNode + 1) >= nodeCount;
    }

    // Can the data array accomodate (size * sizeof(int)) more bytes?
    inline bool isDataStorageFull(int size) const {
      return (lastData + size) >= dataCount;
    }

    // Helper for getFreeNode; tries to find a previously recycled node.
    inline int getRecycledNode(int size) {
#ifdef RECYCLED_LIST
      if (recycledFront) {
        recycled_list* curr = recycledFront;
        while (curr && curr->size != size) { curr = curr->next; }
        if (curr && curr->front != -1) {
          int node = curr->front;
          curr->front = nodes[node].next;
          nodes[node].next = getNull();
          // no need to clean up data array here; look at recycleNode(..)
          return node;
        }
      }
      return -1;                        // did not find recycled node
#else
      if (int(holes.size()) <= size) return -1;
      if (holes[size] == -1) return -1;
      int node = holes[size];
      holes[size] = nodes[node].next;
      nodes[node].next = getNull();
      return node;
#endif
    }

    // Returns a free node with a data array of (size * sizeof(int)) bytes.
    inline int getFreeNode(int size) {
      // try to find a recycled node that fits
      int newNode = getRecycledNode(size);

      // if not found, get a new node
      if (newNode == -1) {
        // expand nodes and data arrays if necessary
        if (isNodesStorageFull()) expandNodes();
        if (isDataStorageFull(size)) expandData();
        newNode = ++lastNode;
        setDataOffset(nodes[newNode], lastData + 1);
        // nodes[newNode].data = data + lastData + 1;
        lastData += size;
      }

      return newNode;
    }

    // Places the node in the recyled nodes list.
    // Note: Usually called after owner->discardEntry(..).
    inline void recycleNode(int h) {
      int nodeSize = nodes[h].op->getCacheEntryLength();

#ifdef RECYCLED_LIST

      // Holes stored in a list of lists, where each list contains nodes of
      // a certain size

      // find correct list
      recycled_list* curr = recycledFront;
      while(curr && curr->size != nodeSize) {
        curr = curr->next;
      }

      // if corresponding list does not exist, create one
      if (0 == curr) {
        // create a new recycled list for nodes of this size
        curr = (recycled_list *) malloc(sizeof(recycled_list));
        curr->size = nodeSize;
        curr->next = recycledFront;
        recycledFront = curr;
        curr->front = -1;
      }

      // add node to front of corresponding list
      nodes[h].op = 0;
      nodes[h].next = curr->front;
      curr->front = h;

#else

      // Holes stored in a vector of lists, where each list contains nodes of a
      // certain size. The index of the vector gives the size of the nodes in
      // the corresponding list.
  
      // enlarge vector if necessary
      if (int(holes.size()) <= nodeSize) {
        holes.resize(nodeSize + 1, -1);
      }

      // add node to front of corresponding list
      nodes[h].owner = 0;
      nodes[h].next = holes[nodeSize];
      holes[nodeSize] = h;
#endif
    }

  public:
    // ******************************************************************
    // *                    functions for  hash table                   *
    // ******************************************************************

    // which integer handle to use for NULL.
    inline int getNull() const { 
      return -1; 
    }

    // next node in this chain 
    inline int getNext(int h) const { 
      return nodes[h].next;
    }

    // set the "next node" field for h 
    inline void setNext(int h, int n) { 
      nodes[h].next = n; 
    }

    // compute hash value
    inline unsigned hash(int h, unsigned n) const {
      return raw_hash(getDataAddress(nodes[h]), nodes[h].op->getKeyLength()) % n;
    }

    // is key(h1) == key(h2)?
    inline bool equals(int h1, int h2) const {
      if (nodes[h1].op != nodes[h2].op) return false;
      return (0 == memcmp(
        getDataAddress(nodes[h1]), 
        getDataAddress(nodes[h2]), 
        nodes[h1].op->getKeyLengthInBytes()
      ));
    }

    // is h a stale (invalid/unwanted) node
    inline bool isStale(int h) const {
      return
        nodes[h].op->isMarkedForDeletion() ||
        nodes[h].op->isEntryStale(getDataAddress(nodes[h]));
    }

    // node h has been removed from the table, perform clean up of node
    // (decrement table count, release node, etc.)
    inline void uncacheNode(int h) {
      nodes[h].op->discardEntry(getDataAddress(nodes[h]));
      recycleNode(h);
    }

    // for debugging
    virtual void show(FILE* s, int h) const;
    virtual void show(FILE *s, bool verbose = false) const;
};

// **********************************************************************
// *                                                                    *
// *                      monolithic_table methods                      *
// *                                                                    *
// **********************************************************************

MEDDLY::monolithic_table::monolithic_table(settings s)
 : compute_table(s)
{
  // Initialize node array
  nodeCount = 1024;
  lastNode = -1;
  nodes = (table_entry *) malloc(nodeCount * sizeof(table_entry));
  if (0==nodes)
    throw error(error::INSUFFICIENT_MEMORY);
  for (int i=0; i<nodeCount; i++) {
    nodes[i].op = 0;
    nodes[i].next = getNull();
    setDataOffset(nodes[i], -1);
  }

  // Initialize data array
  dataCount = 1024;
  lastData = -1;
  data = (int *) malloc(dataCount * sizeof(int));
  if (0==data) 
    throw error(error::INSUFFICIENT_MEMORY);
  memset(data, 0, dataCount * sizeof(int));

  // Initialize recycled lists
  recycledNodes = -1;
  recycledFront = 0;

  // Initialize actual tables
  if (chaining) {
    // create hash table with chaining
#ifdef USE_CHAINED_HASH_TABLE
    ht = new chained_hash_table<monolithic_table>(this, maxSize);
#else
    ht = new hash_table<monolithic_table>(this, maxSize);
#endif
    fsht = 0;
  } else {
    // create hash table with no chaining
    fsht = new fixed_size_hash_table<monolithic_table>(this, maxSize);
    ht = 0;
  }
}

MEDDLY::monolithic_table:: ~monolithic_table()
{
  // delete hash table
  if (ht) { delete ht; ht = 0; }
  if (fsht) { delete fsht; fsht = 0; }

  // free data and nodes arrays
  free(data);
  free(nodes);
}

void MEDDLY::monolithic_table::add(operation* op, const int* entry)
{
#ifdef DEBUG_CT
  fprintf(stderr, "MT adding entry ");
  op->showEntry(stderr, entry);
  fprintf(stderr, "\n");
#endif
  // copy entry data
  int node = getFreeNode(op->getCacheEntryLength());
  nodes[node].op = op;
  memcpy(getDataAddress(nodes[node]), entry, op->getCacheEntryLengthInBytes());
  nodes[node].next = getNull();
  // insert node into hash table
  if (ht) ht->insert(node); else fsht->insert(node);
}

const int* MEDDLY::monolithic_table::find(operation* op, const int* entryKey)
{
  perf.pings++;
  // build key node
  int node = getFreeNode(op->getCacheEntryLength());
  nodes[node].op = op;
  // copying only keyLength data because hash() ignores the rest anyway.
  // moreover the key is all the user is reqd to provide
#ifdef DEVELOPMENT_CODE
  memset(getDataAddress(nodes[node]), 0, op->getCacheEntryLengthInBytes());
#endif
  memcpy(getDataAddress(nodes[node]), entryKey, op->getKeyLengthInBytes());
  nodes[node].next = getNull();

  // search for key node
  int h = (ht)? ht->find(node): fsht->find(node);

  // recycle key node
  recycleNode(node);

  if (h != getNull()) {
    // found entry
    perf.hits++;
#ifdef DEBUG_CT
    fprintf(stderr, "MT found entry ");
    op->showEntry(stderr, getDataAddress(nodes[h]));
    fprintf(stderr, "\n");
#endif
    return getDataAddress(nodes[h]);
  }
  // did not find entry
  return 0;
}

void MEDDLY::monolithic_table::removeStales()
{
  static bool removingStales = false;
  if (!removingStales) {
    removingStales = true;
    if (ht) ht->removeStaleEntries();
    if (fsht) fsht->removeStaleEntries();
    removingStales = false;
  }
}

void MEDDLY::monolithic_table::removeAll()
{
  // go through all nodes and remove each valid node
  table_entry* end = nodes + lastNode + 1;
  for (table_entry* curr = nodes; curr != end; ++curr) {
    if (!isFreeNode(*curr)) {
      DCASSERT(curr->dataOffset != -1);
      if (ht) ht->remove(curr - nodes);
      else    fsht->remove(curr - nodes);
    }
  }
}

void MEDDLY::monolithic_table::updateStats()
{
  DCASSERT(ht != 0 || fsht != 0);
  perf.numEntries = (ht)? ht->getEntriesCount(): fsht->getEntriesCount();
}

void MEDDLY::monolithic_table::expandNodes()
{
  DCASSERT(nodeCount != 0);
  int newNodeCount = int(nodeCount * expansionFactor);
  table_entry* tempNodes =
    (table_entry *) realloc(nodes, newNodeCount * sizeof(table_entry));
  if (0 == tempNodes)
    throw error(error::INSUFFICIENT_MEMORY);
  nodes = tempNodes;
  for (int i = nodeCount; i < newNodeCount; ++i) {
    nodes[i].op = 0;
    // nodes[i].data = 0;
    nodes[i].next = getNull();
    setDataOffset(nodes[i], -1);
  }
  nodeCount = newNodeCount;
}


void MEDDLY::monolithic_table::expandData()
{
  int newDataCount = int(dataCount * expansionFactor);
  data = (int *) realloc(data, newDataCount * sizeof(int));
  if (0 == data) 
    throw error(error::INSUFFICIENT_MEMORY);
  memset(data + dataCount, 0, (newDataCount - dataCount) * sizeof(int));
  dataCount = newDataCount;
}


void MEDDLY::monolithic_table::show(FILE* s, int h) const
{
  nodes[h].op->showEntry(s, getDataAddress(nodes[h]));
}


void MEDDLY::monolithic_table::show(FILE *s, bool verbose) const
{ 
  char filler[] = "\t";
  fprintf(s, "Monolithic compute table\n");
  fprintf(s, "%sNumber of slots:\t%d\n", filler, nodeCount);
  fprintf(s, "%sMemory usage:   \t%lu\n",
      filler, (unsigned long)(
        dataCount * sizeof(int) +
        nodeCount * sizeof(table_entry) +
        (ht == 0? 0: ht->getMemoryUsage()) +
        (fsht == 0? 0: fsht->getMemoryUsage())));
  if (verbose) {
    fprintf(s, "%s  Nodes[]:      \t%lu\n",
        filler, (unsigned long)(nodeCount * sizeof(table_entry)));
    fprintf(s, "%s  Data[]:       \t%lu\n",
        filler, (unsigned long)(dataCount * sizeof(int)));
  }
  fprintf(s, "%sPings:          \t%d\n", filler, perf.pings);
  fprintf(s, "%sHits:           \t%d\n", filler, perf.hits);
  fprintf(s, "Internal hash table info:\n");
  DCASSERT(ht == 0 || fsht == 0);
  if (ht != 0) ht->show(s, verbose);
  if (fsht != 0) fsht->show(s, verbose);
}


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                       operation_table  class                       *
// *                                                                    *
// *                                                                    *
// **********************************************************************

class MEDDLY::operation_table : public compute_table {
  public:
    operation_table(settings s, operation* op);
    virtual ~operation_table();

    // required functions

    virtual bool isOperationTable() const   { return false; }
    virtual void add(operation* op, const int* entry);
    virtual const int* find(operation* op, const int* entryKey);
    virtual void removeStales();
    virtual void removeAll();
    virtual void updateStats();

  private:

    // ******************************************************************
    // *                    Internal data-structures                    *
    // ******************************************************************

    /*  The compute table maintains an array of table entries which can
        be identified by {owner, data}.
    */
    struct table_entry {
      int dataOffset;               // entry data (usually: {key, result})
      int next;                     // for hash table
    };

    // ******************************************************************
    // *                      Internal data                             *
    // ******************************************************************

    operation* global_op;           // Operation everywhere
    table_entry* nodes;             // Array of nodes
    int nodeCount;                  // size of array nodes
    int lastNode;                   // last used index in array nodes
    int* data;                      // Array of ints for storing data
    int dataCount;                  // size of array data
    int lastData;                   // last used index in array data
    int recycledFront;              // first recycled node
    fixed_size_hash_table<operation_table>* fsht;
#ifdef USE_CHAINED_HASH_TABLE
    chained_hash_table<operation_table>* ht;
#else
    hash_table<operation_table>* ht;
#endif

    // ******************************************************************
    // *                      Internal functions                        *
    // ******************************************************************

    // Expand the nodes array (grows by expansionFactor).
    void expandNodes();
    // Expand the data array (grows by expansionFactor).
    void expandData();

    // Return the address of the data associated with this entry
    inline int* getDataAddress(const table_entry& n) const {
      return data + n.dataOffset;
    }

    // Set the address of the data associated with this entry
    inline void setDataAddress(table_entry& n, const int* address) {
      n.dataOffset = address - data;
    }

    // Set the address of the data associated with this entry using offset
    inline void setDataOffset(table_entry& n, int offset) {
      n.dataOffset = offset;
    }

    // Is this a free node (i.e. stores no valid data)
    inline bool isFreeNode(const table_entry& n) const {
      return n.dataOffset >= 0;
    }

    // Is the nodes array full?
    inline bool isNodesStorageFull() const {
      return (lastNode + 1) >= nodeCount;
    }

    // Can the data array accomodate (size * sizeof(int)) more bytes?
    inline bool isDataStorageFull(int size) const {
      return (lastData + size) >= dataCount;
    }

    // Helper for getFreeNode; tries to find a previously recycled node.
    inline int getRecycledNode() {
      if (recycledFront >= 0) {
        int ans = recycledFront;
        recycledFront = nodes[recycledFront].next;
        return ans;
      }
      return -1;                        // did not find recycled node
    }

    // Returns a free node 
    inline int getFreeNode() {
      // try to find a recycled node that fits
      int newNode = getRecycledNode();

      // if not found, get a new node
      if (newNode == -1) {
        int size = global_op->getCacheEntryLength();
        // expand nodes and data arrays if necessary
        if (isNodesStorageFull()) expandNodes();
        if (isDataStorageFull(size)) expandData();
        newNode = ++lastNode;
        setDataOffset(nodes[newNode], lastData + 1);
        lastData += size;
      }

      return newNode;
    }

    // Places the node in the recyled nodes list.
    // Note: Usually called after owner->discardEntry(..).
    inline void recycleNode(int h) {
      nodes[h].next = recycledFront;
      recycledFront = h;
    }

  public:
    // ******************************************************************
    // *                    functions for  hash table                   *
    // ******************************************************************

    // which integer handle to use for NULL.
    inline int getNull() const { 
      return -1; 
    }

    // next node in this chain 
    inline int getNext(int h) const { 
      return nodes[h].next;
    }

    // set the "next node" field for h 
    inline void setNext(int h, int n) { 
      nodes[h].next = n; 
    }

    // compute hash value
    inline unsigned hash(int h, unsigned n) const {
      return raw_hash(getDataAddress(nodes[h]), global_op->getKeyLength()) % n;
    }

    // is key(h1) == key(h2)?
    inline bool equals(int h1, int h2) const {
      return (0 == memcmp(
        getDataAddress(nodes[h1]), 
        getDataAddress(nodes[h2]), 
        global_op->getKeyLengthInBytes()
      ));
    }

    // is h a stale (invalid/unwanted) node
    inline bool isStale(int h) const {
      return
        global_op->isMarkedForDeletion() ||
        global_op->isEntryStale(getDataAddress(nodes[h]));
    }

    // node h has been removed from the table, perform clean up of node
    // (decrement table count, release node, etc.)
    inline void uncacheNode(int h) {
      global_op->discardEntry(getDataAddress(nodes[h]));
      recycleNode(h);
    }

    // for debugging
    virtual void show(FILE* s, int h) const;
    virtual void show(FILE *s, bool verbose = false) const;
};

// **********************************************************************
// *                                                                    *
// *                      operation_table  methods                      *
// *                                                                    *
// **********************************************************************

MEDDLY::operation_table::operation_table(settings s, operation* op)
 : compute_table(s)
{
  if (0==op) 
    throw error(error::INVALID_OPERATION);
  global_op = op;

  // Initialize node array
  nodeCount = 1024;
  lastNode = -1;
  nodes = (table_entry *) malloc(nodeCount * sizeof(table_entry));
  if (0==nodes)
    throw error(error::INSUFFICIENT_MEMORY);
  for (int i=0; i<nodeCount; i++) {
    nodes[i].next = getNull();
    setDataOffset(nodes[i], -1);
  }

  // Initialize data array
  dataCount = 1024;
  lastData = -1;
  data = (int *) malloc(dataCount * sizeof(int));
  if (0==data) 
    throw error(error::INSUFFICIENT_MEMORY);
  memset(data, 0, dataCount * sizeof(int));

  // Initialize recycled lists
  recycledFront = -1;

  // Initialize actual tables
  if (chaining) {
    // create hash table with chaining
#ifdef USE_CHAINED_HASH_TABLE
    ht = new chained_hash_table<operation_table>(this, maxSize);
#else
    ht = new hash_table<operation_table>(this, maxSize);
#endif
    fsht = 0;
  } else {
    // create hash table with no chaining
    fsht = new fixed_size_hash_table<operation_table>(this, maxSize);
    ht = 0;
  }
}

MEDDLY::operation_table:: ~operation_table()
{
  // delete hash table
  if (ht) { delete ht; ht = 0; }
  if (fsht) { delete fsht; fsht = 0; }

  // free data and nodes arrays
  free(data);
  free(nodes);
}

void MEDDLY::operation_table::add(operation* op, const int* entry)
{
  if (op != global_op)
    throw error(error::INVALID_OPERATION);
  // copy entry data
  int node = getFreeNode();
  memcpy(getDataAddress(nodes[node]), entry, op->getCacheEntryLengthInBytes());
  nodes[node].next = getNull();
#ifdef DEBUG_CT
  fprintf(stderr, "OT adding entry ");
  op->showEntry(stderr, entry);
  // fprintf(stderr, " to slot %d", node);
  fprintf(stderr, "\n");
#endif
  // insert node into hash table
  if (ht) ht->insert(node); else fsht->insert(node);
}

const int* MEDDLY::operation_table::find(operation* op, const int* entryKey)
{
  if (op != global_op)
    throw error(error::INVALID_OPERATION);
  perf.pings++;
  // build key node
  int node = getFreeNode();
  // copying only keyLength data because hash() ignores the rest anyway.
  // moreover the key is all the user is reqd to provide
#ifdef DEVELOPMENT_CODE
  memset(getDataAddress(nodes[node]), 0, op->getCacheEntryLengthInBytes());
#endif
  memcpy(getDataAddress(nodes[node]), entryKey, op->getKeyLengthInBytes());
  nodes[node].next = getNull();

  // search for key node
  int h = (ht)? ht->find(node): fsht->find(node);

  // recycle key node
  recycleNode(node);

  if (h != getNull()) {
    // found entry
    perf.hits++;
#ifdef DEBUG_CT
    fprintf(stderr, "OT found entry ");
    op->showEntry(stderr, getDataAddress(nodes[h]));
    fprintf(stderr, " at slot %d\n", h);
#endif
    return getDataAddress(nodes[h]);
  }
  // did not find entry
  return 0;
}

void MEDDLY::operation_table::removeStales()
{
  static bool removingStales = false;
  if (!removingStales) {
    removingStales = true;
    if (ht) ht->removeStaleEntries();
    if (fsht) fsht->removeStaleEntries();
    removingStales = false;
  }
}

void MEDDLY::operation_table::removeAll()
{
  // go through all nodes and remove each valid node
  table_entry* end = nodes + lastNode + 1;
  for (table_entry* curr = nodes; curr != end; ++curr) {
    if (!isFreeNode(*curr)) {
      DCASSERT(curr->dataOffset != -1);
      if (ht) ht->remove(curr - nodes);
      else    fsht->remove(curr - nodes);
    }
  }
}

void MEDDLY::operation_table::updateStats()
{
  DCASSERT(ht != 0 || fsht != 0);
  perf.numEntries = (ht)? ht->getEntriesCount(): fsht->getEntriesCount();
}

void MEDDLY::operation_table::expandNodes()
{
  DCASSERT(nodeCount != 0);
  int newNodeCount = int(nodeCount * expansionFactor);
  table_entry* tempNodes =
    (table_entry *) realloc(nodes, newNodeCount * sizeof(table_entry));
  if (0 == tempNodes)
    throw error(error::INSUFFICIENT_MEMORY);
  nodes = tempNodes;
  for (int i = nodeCount; i < newNodeCount; ++i) {
    // nodes[i].data = 0;
    nodes[i].next = getNull();
    setDataOffset(nodes[i], -1);
  }
  nodeCount = newNodeCount;
}


void MEDDLY::operation_table::expandData()
{
  int newDataCount = int(dataCount * expansionFactor);
  data = (int *) realloc(data, newDataCount * sizeof(int));
  if (0 == data) 
    throw error(error::INSUFFICIENT_MEMORY);
  memset(data + dataCount, 0, (newDataCount - dataCount) * sizeof(int));
  dataCount = newDataCount;
}


void MEDDLY::operation_table::show(FILE* s, int h) const
{
  global_op->showEntry(s, getDataAddress(nodes[h]));
  fprintf(s, "[slot %d] ", h);
}


void MEDDLY::operation_table::show(FILE *s, bool verbose) const
{ 
  char filler[] = "\t";
  fprintf(s, "Compute table for %s operation\n", global_op->getName());
  fprintf(s, "%sNumber of slots:\t%d\n", filler, nodeCount);
  fprintf(s, "%sMemory usage:   \t%lu\n",
      filler, (unsigned long)(
        dataCount * sizeof(int) +
        nodeCount * sizeof(table_entry) +
        (ht == 0? 0: ht->getMemoryUsage()) +
        (fsht == 0? 0: fsht->getMemoryUsage())));
  if (verbose) {
    fprintf(s, "%s  Nodes[]:      \t%lu\n",
        filler, (unsigned long)(nodeCount * sizeof(table_entry)));
    fprintf(s, "%s  Data[]:       \t%lu\n",
        filler, (unsigned long)(dataCount * sizeof(int)));
  }
  fprintf(s, "%sPings:          \t%d\n", filler, perf.pings);
  fprintf(s, "%sHits:           \t%d\n", filler, perf.hits);
  fprintf(s, "Internal hash table info:\n");
  DCASSERT(ht == 0 || fsht == 0);
  if (ht != 0) ht->show(s, verbose);
  if (fsht != 0) fsht->show(s, verbose);
}

#endif

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


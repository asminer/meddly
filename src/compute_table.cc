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

#include "hash.h"
#include "fixed_size_hash.h"
#include "chained_hash.h"
#include <map>

#define RECYCLED_LIST
#define USE_CHAINED_HASH_TABLE

namespace MEDDLY {
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

MEDDLY::compute_table::settings::settings()
{
  chaining = meddlySettings.doComputeTablesUseChaining;
  maxSize = meddlySettings.maxComputeTableSize;
}

MEDDLY::compute_table::compute_table(settings s)
{
  opts = s;

  if (0==opts.maxSize)
    throw error(error::INVALID_ASSIGNMENT);

  perf.numEntries = 0;
  perf.hits = 0;
  perf.pings = 0;
}

MEDDLY::compute_table::~compute_table()
{
}

/*
 * Bob Jenkin's Hash
 * Free to use for educational or commerical purposes
 * http://burtleburtle.net/bob/hash/doobs.html
 */
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

template <class T>
inline void mix(T &a, T &b, T &c)
{
  a -= c;  a ^= rot(c, 4);  c += b; 
  b -= a;  b ^= rot(a, 6);  a += c;
  c -= b;  c ^= rot(b, 8);  b += a;
  a -= c;  a ^= rot(c,16);  c += b;
  b -= a;  b ^= rot(a,19);  a += c;
  c -= b;  c ^= rot(b, 4);  b += a;
}

template <class T>
inline void final(T &a, T &b, T &c)
{
  c ^= b; c -= rot(b,14);
  a ^= c; a -= rot(c,11);
  b ^= a; b -= rot(a,25);
  c ^= b; c -= rot(b,16);
  a ^= c; a -= rot(c,4); 
  b ^= a; b -= rot(a,14);
  c ^= b; c -= rot(b,24);
}

unsigned MEDDLY::compute_table::raw_hash(const int* k, int length)
{
  unsigned long a, b, c;
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
    // Note: Usually called after owner->op->discardEntry(..).
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
  if (opts.chaining) {
    // create hash table with chaining
#ifdef USE_CHAINED_HASH_TABLE
    ht = new chained_hash_table<monolithic_table>(this, opts.maxSize);
#else
    ht = new hash_table<monolithic_table>(this, opts.maxSize);
#endif
    fsht = 0;
  } else {
    // create hash table with no chaining
    fsht = new fixed_size_hash_table<monolithic_table>(this, opts.maxSize);
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
    int recycledNodes;              // -1 indicates no recycled nodes
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
    // Note: Usually called after owner->op->discardEntry(..).
    inline void recycleNode(int h) {
      nodes[h].dataOffset = -1;
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
  recycledNodes = -1;
  recycledFront = 0;

  // Initialize actual tables
  if (opts.chaining) {
    // create hash table with chaining
#ifdef USE_CHAINED_HASH_TABLE
    ht = new chained_hash_table<operation_table>(this, opts.maxSize);
#else
    ht = new hash_table<operation_table>(this, opts.maxSize);
#endif
    fsht = 0;
  } else {
    // create hash table with no chaining
    fsht = new fixed_size_hash_table<operation_table>(this, opts.maxSize);
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
}


void MEDDLY::operation_table::show(FILE *s, bool verbose) const
{ 
  char filler[] = "\t";
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
// *                             Front  End                             *
// *                                                                    *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table*
MEDDLY::createMonolithicTable(compute_table::settings s)
{
  return new monolithic_table(s);
}

MEDDLY::compute_table*
MEDDLY::createOperationTable(compute_table::settings s, operation* op)
{
  return new operation_table(s, op);
}

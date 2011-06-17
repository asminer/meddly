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

MEDDLY::compute_table::settings::settings()
{
  chaining = true;
  maxSize = INT_MAX;
}

MEDDLY::compute_table::compute_table(settings s)
{
}

MEDDLY::compute_table::~compute_table()
{
}


// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************
// **********************************************************************

#include "hash.h"
#include "fixed_size_hash.h"
#include "chained_hash.h"
#include <map>

// const int maxEntries = 262144 * 2;
#define RECYCLED_LIST
#define USE_CHAINED_HASH_TABLE


const float expansionFactor = 1.5;

/* compute table methods */

MEDDLY::compute_table::compute_table()
: nodes(0), nodeCount(1024), lastNode(-1),
  data(0), dataCount(1024), lastData(-1),
  recycledNodes(-1), recycledFront(0),
  ht(0), hits(0), pings(0)
{
  // initialize node and data arrays
  nodes = (table_entry *) malloc(nodeCount * sizeof(table_entry));
  if (nodes == NULL) outOfMemory();
  for (int i = 0; i < nodeCount; ++i)
  {
    nodes[i].owner = 0;
    nodes[i].next = getNull();
    setDataOffset(nodes[i], -1);
  }
  data = (int *) malloc(dataCount * sizeof(int));
  if (data == NULL) outOfMemory();
  memset(data, 0, dataCount * sizeof(int));

  // create new hash table
#ifdef USE_CHAINED_HASH_TABLE
  ht = new chained_hash_table<compute_table>(this, 262144*4);
#else
  ht = new hash_table<compute_table>(this, 262144*4);
#endif
  fsht = 0;
}


MEDDLY::compute_table::~compute_table()
{
  // delete hash table
  if (ht) { delete ht; ht = 0; }
  if (fsht) { delete fsht; fsht = 0; }

  // go through all nodes and call discardEntry for each valid node
  table_entry* end = nodes + lastNode + 1;
  for (table_entry* curr = nodes; curr != end; ++curr)
  {
    if (!isFreeNode(*curr)) {
      DCASSERT(curr->dataOffset != -1);
      curr->owner->op->discardEntry(curr->owner, getDataAddress(*curr));
    }
    // note: we are not recycling the nodes here since we are just going
    // to delete this structure
  }

  // free data and nodes arrays
  free(data);
  free(nodes);
}


void MEDDLY::compute_table::setPolicy(bool chaining, unsigned maxSize)
{
  if (0==maxSize)
    throw error(error::INVALID_ASSIGNMENT);
  // some data is already in table; abort
  if (lastData > -1) 
    throw error(error::MISCELLANEOUS);

  // delete existing hash tables
  clear();
  if (fsht != 0) { delete fsht; fsht = 0; }
  if (ht != 0) { delete ht; ht = 0; }

  if (chaining) {
    // create hash table with chaining
#ifdef USE_CHAINED_HASH_TABLE
    ht = new chained_hash_table<compute_table>(this, maxSize);
#else
    ht = new hash_table<compute_table>(this, maxSize);
#endif
  }
  else {
    // create hash table with no chaining
    fsht = new fixed_size_hash_table<compute_table>(this, maxSize);
  }
}


void MEDDLY::compute_table::show(FILE* s, int h) const
{
  nodes[h].owner->op->showEntry(nodes[h].owner, s, getDataAddress(nodes[h]));
}


void MEDDLY::compute_table::show(FILE *s, bool verbose) const
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
  fprintf(s, "%sPings:          \t%d\n", filler, pings);
  fprintf(s, "%sHits:           \t%d\n", filler, hits);
  fprintf(s, "Internal hash table info:\n");
  DCASSERT(ht == 0 || fsht == 0);
  if (ht != 0) ht->show(s, verbose);
  if (fsht != 0) fsht->show(s, verbose);
}


void MEDDLY::compute_table::expandNodes()
{
  DCASSERT(nodeCount != 0);
  int newNodeCount = int(nodeCount * expansionFactor);
  table_entry* tempNodes =
    (table_entry *) realloc(nodes, newNodeCount * sizeof(table_entry));
  if (tempNodes == NULL) outOfMemory();
  nodes = tempNodes;
  for (int i = nodeCount; i < newNodeCount; ++i)
  {
    nodes[i].owner = 0;
    // nodes[i].data = 0;
    nodes[i].next = getNull();
    setDataOffset(nodes[i], -1);
  }
  nodeCount = newNodeCount;
}


void MEDDLY::compute_table::expandData()
{
  int newDataCount = int(dataCount * expansionFactor);
  data = (int *) realloc(data, newDataCount * sizeof(int));
  if (data == NULL) outOfMemory();
  memset(data + dataCount, 0, (newDataCount - dataCount) * sizeof(int));
  dataCount = newDataCount;
}


void MEDDLY::compute_table::removeStales(op_info* owner)
{
  static bool removingStales = false;
  if (!removingStales) {
    DCASSERT(ht != 0 || fsht != 0);
    if (owner) {
      // for each entry belonging to owner, call isStale() and if necessary
      // hash-table's remove() (which will call uncacheNode() and which
      // in-turn will call discardEntry()).
      if (ht) {
        table_entry* end = nodes + lastNode + 1;
        for (table_entry* current = nodes; current != end; ++current)
        {
          DCASSERT(!removingStales);
          if (current->owner == owner &&
              owner->op->isEntryStale(owner, getDataAddress(*current)))
            ht->remove(current - nodes);
        }
      }
      if (fsht) {
        table_entry* end = nodes + lastNode + 1;
        for (table_entry* current = nodes; current != end; ++current)
        {
          DCASSERT(!removingStales);
          if (current->owner == owner &&
              owner->op->isEntryStale(owner, getDataAddress(*current)))
            fsht->remove(current - nodes);
        }
      }
    }
    else {
      removingStales = true;
      if (ht) ht->removeStaleEntries();
      if (fsht) fsht->removeStaleEntries();
      removingStales = false;
    }
  }
}


void MEDDLY::compute_table::removeEntries(op_info* owner)
{
  static bool removingEntries = false;
  if (!removingEntries) {
    DCASSERT(ht != 0 || fsht != 0);
    if (owner) {
      // for each entry belonging to owner, call hash-table's remove()
      // which will call uncacheNode(), which in-turn will call discardEntry()).
      if (ht) {
        table_entry* end = nodes + lastNode + 1;
        for (table_entry* current = nodes; current != end; ++current)
        {
          DCASSERT(!removingEntries);
          if (current->owner == owner) ht->remove(current - nodes);
        }
      }
      if (fsht) {
        table_entry* end = nodes + lastNode + 1;
        for (table_entry* current = nodes; current != end; ++current)
        {
          DCASSERT(!removingEntries);
          if (current->owner == owner) fsht->remove(current - nodes);
        }
      }
    }
    else {
      removingEntries = true;
      clear();
      removingEntries = false;
    }
  }
}


void MEDDLY::compute_table::clear()
{
  DCASSERT(ht != 0 || fsht != 0);
  if (ht) ht->clear();
  if (fsht) fsht->clear();
}

int MEDDLY::compute_table::getNumEntries() const
{
  DCASSERT(ht != 0 || fsht != 0);
  return (ht)? ht->getEntriesCount(): fsht->getEntriesCount();
}



// ****************************************************************************
//
//                          Binary Compute Cache
//
// ****************************************************************************


MEDDLY::binary_compute_cache::binary_compute_cache()
: hits(0), pings(0), adds(0), inserts(0), op(0), f0(0), f1(0), f2(0),
  checkForStales(true)
{
  compute_cache::clear();
}


MEDDLY::binary_compute_cache::binary_compute_cache(const old_operation* op,
  const op_param* plist, int n)
: hits(0), pings(0), adds(0), inserts(0)
{
  assert(n == 3);
  for (int i = 0; i < 3; ++i) assert(plist[i].isForest());
  this->op = op;
  f0 = const_cast<expert_forest*>(plist[0].readForest());
  f1 = const_cast<expert_forest*>(plist[1].readForest());
  f2 = const_cast<expert_forest*>(plist[2].readForest());

  // Set the checkForStales flag based on f2's node deletion scheme.
  // Ignoring f0 and f1 since the result belongs to f2, and the result
  // is what will be accessed by the compute table user.
  checkForStales =
    f2->getNodeDeletion() == forest::PESSIMISTIC_DELETION;
  
  compute_cache::clear();
}


void MEDDLY::binary_compute_cache::set(const old_operation* op,
  expert_forest* f0, expert_forest* f1, expert_forest* f2)
{
  // op is allowed to be null.
  assert(f0 != 0 && f1 != 0 && f2 != 0);
  clear();
  hits = pings = adds = inserts = 0;
  this->op = op;
  this->f0 = f0;
  this->f1 = f1;
  this->f2 = f2;

  // See constructor for the logic used for setting checkForStales.
  checkForStales =
    f2->getNodeDeletion() == forest::PESSIMISTIC_DELETION;
}


MEDDLY::binary_compute_cache::~binary_compute_cache()
{
  clear();
}


void MEDDLY::binary_compute_cache::setPolicy(bool chaining, unsigned maxSize)
{ 
  clear();
}


void MEDDLY::binary_compute_cache::add(op_info* owner, const int* entry)
{
  assert(false);
  add(owner, entry[0], entry[1], entry[2]);
}


const int* MEDDLY::binary_compute_cache::find(op_info* owner, const int* entryKey) {
  assert(false);
  static int result[3];
  result[0] = entryKey[0];
  result[1] = entryKey[1];
  return find(owner, result[0], result[1], result[2])? result: 0;
}


const char* MEDDLY::binary_compute_cache::getOpName() const
{
  return op? op->getName(): "Unnamed Operation";
}


void MEDDLY::binary_compute_cache::removeStales(op_info* owner)
{
  adds = 0;
  if (ct.empty()) return;
  if (owner) { assert(owner->cc == this); }
  int staleCount = 0;
  int count = 0;
  std::map< key_type, ans_type >::iterator curr = ct.begin();
  std::map< key_type, ans_type >::iterator end = ct.end();
  while (curr != end) {
    if (f0->isStale((curr->first).a) ||
        f1->isStale((curr->first).b) ||
        f2->isStale(curr->second)) {
      f0->uncacheNode((curr->first).a);
      f1->uncacheNode((curr->first).b);
      f2->uncacheNode(curr->second);
      ct.erase(curr++);
      ++staleCount;
    } else {
      ++curr;
    }
    ++count;
  }
  // fprintf(stderr, "Removed %d stale entries out of %d entries in %s\n",
  //    staleCount, count, getOpName());
}


void MEDDLY::binary_compute_cache::clear()
{
  adds = 0;
  if (ct.empty()) return;
  // fprintf(stderr, "Removing all entries in %s\n", getOpName());
  std::map< key_type, ans_type >::iterator curr = ct.begin();
  std::map< key_type, ans_type >::iterator end = ct.end();
  while (curr != end) {
    f0->uncacheNode((curr->first).a);
    f1->uncacheNode((curr->first).b);
    f2->uncacheNode(curr->second);
    ct.erase(curr++);
  }
  DCASSERT(ct.empty());
}


void MEDDLY::binary_compute_cache::removeEntries(op_info* owner)
{
  if (owner) { assert(owner->cc == this); }
  clear();
}


// for debugging
void MEDDLY::binary_compute_cache::show(FILE *s, bool verbose) const
{
  fprintf(s, "Compute table for %s\n", getOpName());
  fprintf(s, "Inserts: %d, Pings: %d, Hits: %d, Adds (max %d): %d\n",
      inserts, pings, hits, MEDDLY::binary_compute_cache::maxAdds, adds);
  if (ct.empty()) return;
  if (verbose) {
    std::map< key_type, ans_type >::const_iterator curr = ct.begin();
    std::map< key_type, ans_type >::const_iterator end = ct.end();
    while (curr != end) {
      fprintf(s, "[%d, %d]: %d\n",
          (curr->first).a, (curr->first).b, curr->second);
      ++curr;
    }
  }
}



//////////////////////////////////////////////////////////////////////


    private:

      // ******************************************************************
      // *                    Internal data-structures                    *
      // ******************************************************************

      /* The compute table maintains an array of table entries which can
      * be identified by {owner, data}.
      */
      struct table_entry {
        operation* op;                // operation to which entry belongs to
        int dataOffset;               // entry data (usually: {key, result})
        int next;                     // for hash table
      };

      /* Cache entries that have been discard (or of no use anymore) are
      * recycled and placed on a recycled_list corresponding to the size
      * of the table entry (size of table entry is determined by the size
      * of its data array).
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
      int recycledNodes;              // -1 indicates no recycled nodes
      recycled_list* recycledFront;
      std::vector<int> holes;
      fixed_size_hash_table<compute_table>* fsht;
#ifdef USE_CHAINED_HASH_TABLE
      chained_hash_table<compute_table>* ht;
#else
      hash_table<compute_table>* ht;
#endif
      unsigned hits;
      unsigned pings;

      // ******************************************************************
      // *                      Internal functions                        *
      // ******************************************************************

      // Return the address of the data associated with this entry
      int* getDataAddress(const table_entry& n) const;
      // Set the address of the data associated with this entry
      void setDataAddress(table_entry& n, int* address);
      // Set the address of the data associated with this entry using offset
      void setDataOffset(table_entry& n, int offset);
      // Is this a free node (i.e. stores no valid data)
      bool isFreeNode(const table_entry& n) const;
      // Returns a free node with a data array of (size * sizeof(int)) bytes.
      int getFreeNode(int size);
      // Helper for getFreeNode; tries of find a previously recycled node.
      int getRecycledNode(int size);
      // Places the node in the recyled nodes list.
      // Note: Usually called after owner->op->discardEntry(..).
      void recycleNode(int h);

      // Is the nodes array full?
      bool isNodesStorageFull() const;
      // Expand the nodes array (grows by a factor of 2).
      void expandNodes();
      // Can the data array accomodate (size * sizeof(int)) more bytes?
      bool isDataStorageFull(int size) const;
      // Expand the data array (grows by a factor of 2).
      void expandData();

    public:
      // ******************************************************************
      // *                    functions for  hash table                   *
      // ******************************************************************

      // which integer handle to use for NULL.
      inline int getNull() const          { return -1; }

      // next node in this chain 
      inline int getNext(int h) const     { return nodes[h].next; }

      // set the "next node" field for h 
      inline void setNext(int h, int n)   { nodes[h].next = n; }

      // compute hash value
      unsigned hash(int h, unsigned n) const;

      // is key(h1) == key(h2)?
      inline bool equals(int h1, int h2) const {
          if (nodes[h1].op != nodes[h2].op) return false;
          return 0 == memcmp(
            getDataAddress(nodes[h1]), 
            getDataAddress(nodes[h2]), 
            nodes[h1].op->getKeyLengthInBytes()
          );
      }

      // is h a stale (invalid/unwanted) node
      inline bool isStale(int h) const {
          return nodes[h].op->isEntryStale(nodes[h].owner, getDataAddress(nodes[h]));

inline
bool MEDDLY::compute_table::isStale(int h) const
{
}


      // node h has been removed from the table, perform clean up of node
      // (decrement table count, release node, etc.)
      void uncacheNode(int h);

      // for debugging
      virtual void show(FILE* s, int h) const;
      virtual void show(FILE *s, bool verbose = false) const;

  };


// **********************************************************************
//
//                    Inlined methods for compute_table
//
// **********************************************************************

inline
void MEDDLY::compute_table::add(op_info* owner, const int* entry)
{
  // copy entry data
  int node = getFreeNode(owner->op->getCacheEntryLength());
  nodes[node].owner = owner;
  memcpy(getDataAddress(nodes[node]), entry,
      owner->op->getCacheEntryLengthInBytes());
  nodes[node].next = getNull();
  // insert node into hash table
  if (ht) ht->insert(node); else fsht->insert(node);
}

inline
void MEDDLY::compute_table::add(op_info* owner, int a, int b, int c)
{
  // copy entry data
  DCASSERT(owner->op->getKeyLength() == 2 &&
      owner->op->getCacheEntryLength() == 3);
  int node = getFreeNode(3);
  nodes[node].owner = owner;
  int* data = getDataAddress(nodes[node]);
  // { *data++ = a; } is the same as { *data = a; data++; }
  *data++ = a;
  *data++ = b;
  *data++ = c;
  nodes[node].next = getNull();
  // insert node into hash table
  if (ht) ht->insert(node); else fsht->insert(node);
}

inline
const int* MEDDLY::compute_table::find(op_info* owner, const int* entryKey)
{
  pings++;
  int keyLength = owner->op->getKeyLength();
  int entryLength = owner->op->getCacheEntryLength();
  // build key node
  int node = getFreeNode(entryLength);
  nodes[node].owner = owner;
  // copying only keyLength data because hash() ignores the rest anyway.
  // moreover the key is all the user is reqd to provide
  memcpy(getDataAddress(nodes[node]), entryKey, sizeof(int) * keyLength);
#ifdef DEVELOPMENT_CODE
  memset(getDataAddress(nodes[node]) + keyLength, 0,
      sizeof(int) * (entryLength - keyLength));
#endif
  nodes[node].next = getNull();

  // search for key node
  int h = (ht)? ht->find(node): fsht->find(node);

  // recycle key node
  recycleNode(node);

  if (h != getNull()) {
    // found entry
    hits++;
    return getDataAddress(nodes[h]);
  }
  // did not find entry
  return 0;
}

inline
bool MEDDLY::compute_table::find(op_info* owner, int a, int b, int& c)
{
  DCASSERT(owner->op->getKeyLength() == 2 &&
      owner->op->getCacheEntryLength() == 3);
  pings++;

  // build key node
  int node = getFreeNode(3);
  nodes[node].owner = owner;
  // copying only keyLength data because hash() ignores the rest anyway.
  // moreover the key is all the user is reqd to provide
  int* data = getDataAddress(nodes[node]);
  // { *data++ = a; } is the same as { *data = a; data++; }
  *data++ = a;
  *data++ = b;
#ifdef DEVELOPMENT_CODE
  *data++ = 0;
#endif
  nodes[node].next = getNull();

  // search for key node
  int h = (ht)? ht->find(node): fsht->find(node);

  // recycle key node
  recycleNode(node);

  if (h != getNull()) {
    // found entry
    hits++;
    c = getDataAddress(nodes[h])[2];
    return true;
  }
  // did not find entry
  return false;
}


/*
 * Hash table functions
 */


/*
 * Bob Jenkin's Hash
 * Free to use for educational or commerical purposes
 * http://burtleburtle.net/bob/hash/doobs.html
 */
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))
#define mix(a,b,c) \
  { \
      a -= c;  a ^= rot(c, 4);  c += b; \
      b -= a;  b ^= rot(a, 6);  a += c; \
      c -= b;  c ^= rot(b, 8);  b += a; \
      a -= c;  a ^= rot(c,16);  c += b; \
      b -= a;  b ^= rot(a,19);  a += c; \
      c -= b;  c ^= rot(b, 4);  b += a; \
  }
#define final(a,b,c) \
  { \
      c ^= b; c -= rot(b,14); \
      a ^= c; a -= rot(c,11); \
      b ^= a; b -= rot(a,25); \
      c ^= b; c -= rot(b,16); \
      a ^= c; a -= rot(c,4);  \
      b ^= a; b -= rot(a,14); \
      c ^= b; c -= rot(b,24); \
  }
inline
unsigned smallFinal(unsigned a, unsigned b, unsigned c)
{
  c ^= b; c -= rot(b,14);
  a ^= c; a -= rot(c,11);
  b ^= a; b -= rot(a,25);
  c ^= b; c -= rot(b,16);
#if 1
  a ^= c; a -= rot(c,4);
  b ^= a; b -= rot(a,14);
  c ^= b; c -= rot(b,24);
#endif
  return c;
}

inline
unsigned MEDDLY::compute_table::hash(int h, unsigned n) const
{
  int length = nodes[h].owner->op->getKeyLength();
  int* k = getDataAddress(nodes[h]);

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

  return c % n;
}

inline
void MEDDLY::compute_table::uncacheNode(int h)
{
  nodes[h].owner->op->discardEntry(nodes[h].owner, getDataAddress(nodes[h]));
  recycleNode(h);
}

inline
int* MEDDLY::compute_table::getDataAddress(const table_entry& n) const
{
  return data + n.dataOffset;
}

inline
void MEDDLY::compute_table::setDataOffset(table_entry& n, int offset)
{
  n.dataOffset = offset;
}

inline
void MEDDLY::compute_table::setDataAddress(table_entry& n, int* address)
{
  n.dataOffset = address - data;
}

inline
bool MEDDLY::compute_table::isFreeNode(const table_entry& n) const
{
  return n.owner == NULL;
}

inline
void MEDDLY::compute_table::recycleNode(int h)
{
  int nodeSize = nodes[h].owner->op->getCacheEntryLength();

#ifdef RECYCLED_LIST

  // Holes stored in a list of lists, where each list contains nodes of
  // a certain size

  // find correct list
  recycled_list* curr = recycledFront;
  while(curr != NULL && curr->size != nodeSize)
  {
    curr = curr->next;
  }

  // if corresponding list does not exist, create one
  if (curr == NULL) {
    // create a new recycled list for nodes of this size
    curr = (recycled_list *) malloc(sizeof(recycled_list));
    curr->size = nodeSize;
    curr->next = recycledFront;
    recycledFront = curr;
    curr->front = -1;
  }

  // add node to front of corresponding list
  nodes[h].owner = 0;
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

inline
int MEDDLY::compute_table::getRecycledNode(int size)
{
#ifdef RECYCLED_LIST
  if (recycledFront != NULL) {
    recycled_list* curr = recycledFront;
    while (curr != NULL && curr->size != size) { curr = curr->next; }
    if (curr != NULL && curr->front != -1) {
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

inline
bool MEDDLY::compute_table::isNodesStorageFull() const
{
  return (lastNode + 1) >= nodeCount;
}

inline
bool MEDDLY::compute_table::isDataStorageFull(int size) const
{
  return (lastData + size) >= dataCount;
}

inline
int MEDDLY::compute_table::getFreeNode(int size)
{
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




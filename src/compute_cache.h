
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



/*! \file compute_cache.h

    Compute Cache interface.

    This interface is for "expert" interface users who wish to implement
    user-defined operations.

    A compute cache is used to cache the results of operations on MDDs.
    An expert-user wishing to cache results of user-defined operations may
    need to use this interface.

*/

#ifndef COMPUTE_CACHE_H
#define COMPUTE_CACHE_H

#include "defines.h"
#include "hash.h"
#include "fixed_size_hash.h"
#include "chained_hash.h"
#include <map>

// const int maxEntries = 262144 * 2;
#define RECYCLED_LIST
#define USE_CHAINED_HASH_TABLE

namespace MEDDLY {

class compute_cache {
  public:
    // ******************************************************************
    // *                     expert-user interface                      *
    // ******************************************************************

    /// Constructor
    compute_cache();

    /** Destructor. Will discard all cache entries rendering all pointers
        to data within the cache invalid.
    */
    virtual ~compute_cache();

    /** Sets the policy for the hash table used by the cache.
        @param  chaining
                        if true, use cache with chaining,
                        otherwise, use cache without chaining.
        @param  maxSize the maximum size of the hash table.
        @return         true, if the operation completed successfully.
    */
    virtual bool setPolicy(bool chaining, unsigned maxSize);

    /** Add an entry to the compute cache. Note that this cache allows for
        duplicate entries for the same key. The user may use find() before
        using add() to prevent such duplication. A copy of the data
        in entry[] is stored in the cache.
        @param  owner op_info associated with this cache entry
        @param  entry integer array of size owner->op->getCacheEntryLength(),
                      containing the operands and the result to be stored
    */    
    virtual void add(op_info* owner, const int* entry);
    // void add(op_info* owner, int a, int b);
    virtual void add(op_info* owner, int a, int b, int c);

    /** Find an entry in the compute cache based on the key provided.
        If more than an entry with the same key exists, this will return
        the first matching entry.
        @param  owner op_info associated with this cache entry
        @param  entry integer array of size owner->op->getKeyLength(),
                      containing the key to the cache entry to look for
        @return       integer array of size owner->op->getCacheEntryLength()
    */
    virtual const int* find(op_info* owner, const int* entryKey);
    // int find(op_info* owner, int a);
    virtual bool find(op_info* owner, int a, int b, int& c);

    /** Remove stale entries.
        Scans the cache for entries that are no longer valid (i.e. they are
        stale) and removes them. This can be a time-consuming process
        (proportional to the number of cached entries).
        If owner is non-null, all stale entries associated with a particular
        operation are removed.
        @param  owner op_info of the operation whose stale entries are
                      to be removed from the cache. If owner is null,
                      all the stale entries in the cache are removed.
    */
    virtual void removeStales(op_info* op = 0);

    /** Remove operation entries.

        Scans the cache for entries that belong to the specified
        operation and removes them.
        
        This can be a time-consuming process (proportional to
        the number of cached entries).

        If owner is non-null, all entries in the compute table are removed.

        @param  owner op_info of the operation whose entries are
                      to be removed from the cache. If owner is null,
                      all the entries in the cache are removed.
    */
    virtual void removeEntries(op_info* op = 0);

    /** Removes all cached entries.
    */
    virtual void clear();

    /** Get the number of entries in the cache.
        @return     The number of entries in the cache.
    */
    virtual int getNumEntries() const;

    // ******************************************************************
    // *                 end of expert-user interface                   *
    // ******************************************************************


    // ******************************************************************
    // *         functions for hash table (inlined when possible)       *
    // ******************************************************************

    // which integer handle to use for NULL.
    virtual int getNull() const;

    // next node in this chain 
    virtual int getNext(int h) const;

    // set the "next node" field for h 
    virtual void setNext(int h, int n);

    // compute hash value
    virtual unsigned hash(int h, unsigned n) const;

    // is key(h1) == key(h2)?
    virtual bool equals(int h1, int h2) const;

    // is h a stale (invalid/unwanted) node
    virtual bool isStale(int h) const;

    // node h has been removed from the cache, perform clean up of node
    // (decrement cache count, release node, etc.)
    virtual void uncacheNode(int h);

    // for debugging
    virtual void show(FILE* s, int h) const;
    virtual void show(FILE *s, bool verbose = false) const;

  private:

    // ******************************************************************
    // *                    Internal data-structures                    *
    // ******************************************************************

    /* The compute cache maintains an array of cache entries which can
     * be identified by {owner, data}.
     */
    struct cache_entry {
      op_info* owner;                // op_info to which entry belongs to
      int dataOffset;                // entry data (usually: {key, result})
      int next;                      // for hash table
    };

    /* Cache entries that have been discard (or of no use anymore) are
     * recycled and placed on a recycled_list corresponding to the size
     * of the cache entry (size of cache entry is determined by the size
     * of its data array).
     */
    struct recycled_list {
      int size;                         // size of nodes in this list
      int front;                        // first node in list, -1 terminates
      recycled_list* next;              // points to next recycled node list
    };

    // ******************************************************************
    // *                      Internal functions                        *
    // ******************************************************************

    // Return the address of the data associated with this entry
    int* getDataAddress(const cache_entry& n) const;
    // Set the address of the data associated with this entry
    void setDataAddress(cache_entry& n, int* address);
    // Set the address of the data associated with this entry using offset
    void setDataOffset(cache_entry& n, int offset);
    // Is this a free node (i.e. stores no valid data)
    bool isFreeNode(const cache_entry& n) const;
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

    // ******************************************************************
    // *                      Internal data                             *
    // ******************************************************************

    cache_entry* nodes;                 // Array of nodes
    int nodeCount;                      // size of array nodes
    int lastNode;                       // last used index in array nodes
    int* data;                          // Array of ints for storing data
    int dataCount;                      // size of array data
    int lastData;                       // last used index in array data
    int recycledNodes;                  // -1 indicates no recycled nodes
    recycled_list* recycledFront;
    std::vector<int> holes;
    fixed_size_hash_table<compute_cache>* fsht;
#ifdef USE_CHAINED_HASH_TABLE
    chained_hash_table<compute_cache>* ht;
#else
    hash_table<compute_cache>* ht;
#endif
    unsigned hits;
    unsigned pings;
};



// **********************************************************************
//
//  Binary Compute Cache:
//  
//  This class is meant to store the results of operations that follow
//  the form:
//  operator(operandA, operandB) = operandC
//
//  where,
//  operandA belongs to forest0,
//  operandB belongs to forest1, and
//  operandC belongs to forest2.
//
//  There are no restrictions on the types of forests.
//  Also, forest0, forest1 and forest2 may or may not represent the
//  same forest.
//
//  Deprecated methods:
//  (1) bool setPolicy(bool, unsigned)
//  (2) void add(op_info*, const int*)
//  (3) const int* find(op_info*, const int*)
//
// **********************************************************************
class binary_compute_cache : public compute_cache {

  public:

    // Default Constructor: to be used with set().
    // Note that set() is not an inherited method.
    binary_compute_cache();
    virtual void set(const operation* op, expert_forest* f0,
        expert_forest* f1, expert_forest* f2);

    // Constructor
    binary_compute_cache(const operation* op, const op_param* plist, int n);

    // Destructor
    virtual ~binary_compute_cache();

    virtual const char* getOpName() const;
    virtual int getNumEntries() const;

    virtual void add(op_info* owner, int a, int b, int c);
    virtual bool find(op_info* owner, int a, int b, int& c);

    virtual void clear();
    virtual void removeStales(op_info* owner = 0);
    virtual void removeEntries(op_info* owner = 0);

    virtual void show(FILE *s, bool verbose = false) const;

    // Deprecated Methods
    virtual bool setPolicy(bool chaining, unsigned maxSize);
    virtual void add(op_info* owner, const int* entry);
    virtual const int* find(op_info* owner, const int* entryKey);

  private:

    class key {
      public:
        key() : a(0), b(0) {}
        key(int A, int B) : a(A), b(B) {}
        key(const key& k) : a(k.a), b(k.b) {}
        key& operator=(const key& k) {
          if (this != &k) { a = k.a; b = k.b; }
          return *this;
        }
        void set(int A, int B) { a = A; b = B; }
        bool operator<(const key& k) const {
          if (a < k.a) return true;
          if (a > k.a) return false;
          return (b < k.b);
        }

        int a;
        int b;
    };

    typedef key key_type;
    typedef int ans_type;

    std::map< key_type, ans_type > ct;

    unsigned hits;
    unsigned pings;
    unsigned adds;
    unsigned inserts;
    static const unsigned maxAdds = 1000000;
    static const unsigned maxStaleCount = maxAdds/10;

    const operation* op;
    expert_forest* f0;
    expert_forest* f1;
    expert_forest* f2;
    bool checkForStales;
};

} // namespace MEDDLY























// **********************************************************************
//
//                    Inlined methods for compute_cache
//
// **********************************************************************

inline
void MEDDLY::compute_cache::add(op_info* owner, const int* entry)
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
void MEDDLY::compute_cache::add(op_info* owner, int a, int b, int c)
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
const int* MEDDLY::compute_cache::find(op_info* owner, const int* entryKey)
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
bool MEDDLY::compute_cache::find(op_info* owner, int a, int b, int& c)
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

inline
int MEDDLY::compute_cache::getNull() const
{
  return -1;
}

inline
int MEDDLY::compute_cache::getNext(int h) const
{
  return nodes[h].next;
}

inline
void MEDDLY::compute_cache::setNext(int h, int n)
{
  nodes[h].next = n;
}


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
unsigned MEDDLY::compute_cache::hash(int h, unsigned n) const
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
bool MEDDLY::compute_cache::equals(int h1, int h2) const
{
  if (nodes[h1].owner != nodes[h2].owner) return false;
  return 
    (0 == memcmp(getDataAddress(nodes[h1]),
                 getDataAddress(nodes[h2]),
                 nodes[h1].owner->op->getKeyLengthInBytes()));
}

inline
bool MEDDLY::compute_cache::isStale(int h) const
{
  return
    nodes[h].owner->op->isEntryStale(nodes[h].owner, getDataAddress(nodes[h]));
}

inline
void MEDDLY::compute_cache::uncacheNode(int h)
{
  nodes[h].owner->op->discardEntry(nodes[h].owner, getDataAddress(nodes[h]));
  recycleNode(h);
}

inline
int* MEDDLY::compute_cache::getDataAddress(const cache_entry& n) const
{
  return data + n.dataOffset;
}

inline
void MEDDLY::compute_cache::setDataOffset(cache_entry& n, int offset)
{
  n.dataOffset = offset;
}

inline
void MEDDLY::compute_cache::setDataAddress(cache_entry& n, int* address)
{
  n.dataOffset = address - data;
}

inline
bool MEDDLY::compute_cache::isFreeNode(const cache_entry& n) const
{
  return n.owner == NULL;
}

inline
void MEDDLY::compute_cache::recycleNode(int h)
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
int MEDDLY::compute_cache::getRecycledNode(int size)
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
bool MEDDLY::compute_cache::isNodesStorageFull() const
{
  return (lastNode + 1) >= nodeCount;
}

inline
bool MEDDLY::compute_cache::isDataStorageFull(int size) const
{
  return (lastData + size) >= dataCount;
}

inline
int MEDDLY::compute_cache::getFreeNode(int size)
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

// **********************************************************************
//
//                Inlined methods for binary_compute_cache
//
// **********************************************************************

inline
void MEDDLY::binary_compute_cache::add(op_info* owner, int a, int b, int c)
{
  DCASSERT(owner == 0 || owner->cc == this);
  static key k;
  k.a = a; k.b = b;
  f0->cacheNode(a);
  f1->cacheNode(b);
  f2->cacheNode(c);

  // "ct[k] = c" will not overwrite an existing entry.
  // Using map::insert() gives us a means to check whether the new
  // entry was added to the ct.
  // Note that insert() returns pair<iterator, bool>, and the bool is
  // true if the pair was inserted.
  assert((ct.insert(std::make_pair(k, c))).second);

  ++inserts;
#if 0
  if (++adds > binary_compute_cache::maxAdds) {
    removeStales();
  }
#endif
}

inline
bool MEDDLY::binary_compute_cache::find(op_info* owner, int a, int b, int& c)
{
  DCASSERT(owner->cc == this);
  static key k;
  k.a = a; k.b = b;
  ++pings;
  std::map< key_type, ans_type >::iterator ans = ct.find(k);
  if (ans == ct.end()) return false;
  c = ans->second;
  if (checkForStales) {
    if (f2->isStale(c)) {
      ct.erase(ans);
      f0->uncacheNode(a);
      f1->uncacheNode(b);
      f2->uncacheNode(c);
      return false;
    }
  }
  ++hits;
  return true;
}

inline
int MEDDLY::binary_compute_cache::getNumEntries() const
{
  return ct.size();
}


#endif

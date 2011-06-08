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


#include "compute_cache.h"


const float expansionFactor = 1.5;

/* compute cache methods */

MEDDLY::compute_cache::compute_cache()
: nodes(0), nodeCount(1024), lastNode(-1),
  data(0), dataCount(1024), lastData(-1),
  recycledNodes(-1), recycledFront(0),
  ht(0), hits(0), pings(0)
{
  // initialize node and data arrays
  nodes = (cache_entry *) malloc(nodeCount * sizeof(cache_entry));
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
  ht = new chained_hash_table<compute_cache>(this, 262144*4);
#else
  ht = new hash_table<compute_cache>(this, 262144*4);
#endif
  fsht = 0;
}


MEDDLY::compute_cache::~compute_cache()
{
  // delete hash table
  if (ht) { delete ht; ht = 0; }
  if (fsht) { delete fsht; fsht = 0; }

  // go through all nodes and call discardEntry for each valid node
  cache_entry* end = nodes + lastNode + 1;
  for (cache_entry* curr = nodes; curr != end; ++curr)
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


bool MEDDLY::compute_cache::setPolicy(bool chaining, unsigned maxSize)
{
  // some data is already in cache; abort
  if (lastData > -1) return false;

  // delete existing hash tables
  clear();
  if (fsht != 0) { delete fsht; fsht = 0; }
  if (ht != 0) { delete ht; ht = 0; }

  if (chaining) {
    // create hash table with chaining
#ifdef USE_CHAINED_HASH_TABLE
    ht = new chained_hash_table<compute_cache>(this, maxSize);
#else
    ht = new hash_table<compute_cache>(this, maxSize);
#endif
  }
  else {
    // create hash table with no chaining
    fsht = new fixed_size_hash_table<compute_cache>(this, maxSize);
  }
  return true;
}


void MEDDLY::compute_cache::show(FILE* s, int h) const
{
  nodes[h].owner->op->showEntry(nodes[h].owner, s, getDataAddress(nodes[h]));
}


void MEDDLY::compute_cache::show(FILE *s, bool verbose) const
{ 
  char filler[] = "\t";
  fprintf(s, "%sNumber of slots:\t%d\n", filler, nodeCount);
  fprintf(s, "%sMemory usage:   \t%lu\n",
      filler, (unsigned long)(
        dataCount * sizeof(int) +
        nodeCount * sizeof(cache_entry) +
        (ht == 0? 0: ht->getMemoryUsage()) +
        (fsht == 0? 0: fsht->getMemoryUsage())));
  if (verbose) {
    fprintf(s, "%s  Nodes[]:      \t%lu\n",
        filler, (unsigned long)(nodeCount * sizeof(cache_entry)));
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


void MEDDLY::compute_cache::expandNodes()
{
  DCASSERT(nodeCount != 0);
  int newNodeCount = int(nodeCount * expansionFactor);
  cache_entry* tempNodes =
    (cache_entry *) realloc(nodes, newNodeCount * sizeof(cache_entry));
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


void MEDDLY::compute_cache::expandData()
{
  int newDataCount = int(dataCount * expansionFactor);
  data = (int *) realloc(data, newDataCount * sizeof(int));
  if (data == NULL) outOfMemory();
  memset(data + dataCount, 0, (newDataCount - dataCount) * sizeof(int));
  dataCount = newDataCount;
}


void MEDDLY::compute_cache::removeStales(op_info* owner)
{
  static bool removingStales = false;
  if (!removingStales) {
    DCASSERT(ht != 0 || fsht != 0);
    if (owner) {
      // for each entry belonging to owner, call isStale() and if necessary
      // hash-table's remove() (which will call uncacheNode() and which
      // in-turn will call discardEntry()).
      if (ht) {
        cache_entry* end = nodes + lastNode + 1;
        for (cache_entry* current = nodes; current != end; ++current)
        {
          DCASSERT(!removingStales);
          if (current->owner == owner &&
              owner->op->isEntryStale(owner, getDataAddress(*current)))
            ht->remove(current - nodes);
        }
      }
      if (fsht) {
        cache_entry* end = nodes + lastNode + 1;
        for (cache_entry* current = nodes; current != end; ++current)
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


void MEDDLY::compute_cache::removeEntries(op_info* owner)
{
  static bool removingEntries = false;
  if (!removingEntries) {
    DCASSERT(ht != 0 || fsht != 0);
    if (owner) {
      // for each entry belonging to owner, call hash-table's remove()
      // which will call uncacheNode(), which in-turn will call discardEntry()).
      if (ht) {
        cache_entry* end = nodes + lastNode + 1;
        for (cache_entry* current = nodes; current != end; ++current)
        {
          DCASSERT(!removingEntries);
          if (current->owner == owner) ht->remove(current - nodes);
        }
      }
      if (fsht) {
        cache_entry* end = nodes + lastNode + 1;
        for (cache_entry* current = nodes; current != end; ++current)
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


void MEDDLY::compute_cache::clear()
{
  DCASSERT(ht != 0 || fsht != 0);
  if (ht) ht->clear();
  if (fsht) fsht->clear();
}

int MEDDLY::compute_cache::getNumEntries() const
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


MEDDLY::binary_compute_cache::binary_compute_cache(const operation* op,
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


void MEDDLY::binary_compute_cache::set(const operation* op,
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


bool MEDDLY::binary_compute_cache::setPolicy(bool chaining, unsigned maxSize)
{ 
  clear();
  return true;
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


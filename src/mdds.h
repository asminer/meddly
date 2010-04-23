
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



//TODO: need to keep track of memory allocated as well.
//TODO: revamp the node counting mechanism.
//TODO: add a mechanism to node_manager so that reduction rule can
//      be set after instantiation of the node_manager.
//TODO: rename node_manager?
//TODO: IMPORTANT: find a faster way of doing this:
//      if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
//TODO:

#ifndef MDDS_H
#define MDDS_H

#include <fstream>
#include <iostream>
#include <vector>
//#include <map>  // for getCardinality()

#include "defines.h"
#include "mdd_hash.h"

using namespace std;

/// Entities that use node_manager should inherit and implement the following
/// functions which will enable the node_manager to remove nodes that have
/// gone stale (zombie nodes).
class node_manager_user {
  public:

    node_manager_user() {}
    virtual ~node_manager_user() {}
    virtual bool removeStaleEntries() {return false;}
    virtual long getEntriesCount() const {return 0;}
    virtual void getCacheCounts(unsigned count[], unsigned sz) const {}
    virtual unsigned getCacheCount(int p) const {return 0;}
};

/**
 * \var mdd_node
 * An mdd node. For now its handle is a simple integer.
 */
typedef int mdd_node;

typedef struct _evmdd_node
{
	mdd_node node;
	int value;
  _evmdd_node(): node(0), value(INF) {}
  _evmdd_node(mdd_node a, int v): node(a), value(v) {}
  void operator=(const _evmdd_node& a)
  {
    node = a.node;
    value = a.value;
  }
  bool operator!=(const _evmdd_node& a)
  {
    return ((node != a.node) || (value != a.value));
  }
  bool operator==(const _evmdd_node& a)
  {
    return ((node == a.node) && (value == a.value));
  }
  void print(FILE *s)
  {
    fprintf(s, "[n: %d, v: %d]", node, value);
  }
} evmdd_node;


/**
  If active, the node data is as follows:

  [0] incoming pointer count (>=0)
  [1] If this node is in the unique table,
        the next pointer (>=-1)
        Otherwise, a value < -1. 
  [2] Size.
        If > 0, full storage is used.
        If < 0, sparse storage is used.
        If = 0, the node is deleted but still in caches.

  Full storage:
  [3..2+Size]      Downward pointers (>=0)
  [3+Size] index ptr to this node in address array (>=0)

  Sparse storage:
  [3..2+nnz]       Indexes (>=0)
  [3+nnz..2+2*nnz] Downward pointers (>=0)
  [3+2*nnz] index ptr to this node in address array (>=0)
*/
typedef struct {
  int incount;
  int next_ptr;
  int size;
  int *down_ptrs;
  int index;
} mdd_active_node;

/**
  If deleted, the node becomes a hole. There are two kinds depending on their
  location in the "hole grid". The hole grid is a list of holes:

  (holes_bottom)
  holes_of_size_0 (index) -- holes_of_size_0 (non_index) -- (non_index) -- NULL
  |
  holes_of_size_1 -- ""
  |
  :
  :
  (holes_top)
*/

/**
  Index hole data is as follows:
  [0] -size (number of slots in hole)     
  [1] up
  [2] down 
  [3] next pointer (nodes of same size)
  [4..size-2] Unused
  :
  :
  [size-1] -size
*/
typedef struct {
  int size;
  int up_ptr_offset;
  int down_ptr_offset;
  int next_ptr_offset;
  int *data;
} mdd_index_hole;

/**
  Non-index hole data is as follows:
  [0] -size (number of slots in hole)     
  [1] flag (<0, indicates non-index node)
  [3] prev pointer (nodes of same size)
  [3] next pointer (nodes of same size)
  [4..size-2] Unused
  :
  :
  [size-1] -size
*/
typedef struct {
  int size;
  int flag;
  int prev_ptr_offset;
  int next_ptr_offset;
  int *data;
} mdd_non_index_hole;


class node_manager : public expert_forest {
  public:
    node_manager(domain *d, bool rel, range_type t, edge_labeling ev,
        reduction_rule r, node_storage s, node_deletion_policy nd);
    virtual ~node_manager();

    // ------------- inherited from expert_forest ---------------------------

  public:
    /// Refer to meddly.h
    forest::error createEdgeForVar(int vh, bool primedLevel, dd_edge& result);
    forest::error createEdgeForVar(int vh, bool primedLevel,
       bool* terms, dd_edge& a);
    forest::error createEdgeForVar(int vh, bool primedLevel,
       int* terms, dd_edge& a);
    forest::error createEdgeForVar(int vh, bool primedLevel,
       float* terms, dd_edge& a);

    forest::error createSubMatrix(const bool* const* vlist,
        const bool* const* vplist, const dd_edge a, dd_edge& b);
    forest::error createSubMatrix(const dd_edge& rows, const dd_edge& cols,
        const dd_edge& a, dd_edge& b);

    virtual forest::error getElement(const dd_edge& a, int index, int* e);

    virtual forest::error findFirstElement(const dd_edge& f, int* vlist) const;
    virtual forest::error findFirstElement(const dd_edge& f, int* vlist,
        int* vplist) const;

    /// Create a temporary node -- a node that can be modified by the user
    virtual int createTempNode(int lh, int size, bool clear = true);

    /// Create a temporary node with the given downpointers. Note that
    /// downPointers[i] corresponds to the downpointer at index i.
    /// IMPORTANT: The incounts for the downpointers are not incremented.
    /// The returned value is the handle for the temporary node.
    virtual int createTempNode(int lh, std::vector<int>& downPointers);

    /// Same as createTempNode(int, vector<int>) except this is for EV+MDDs.
    virtual int createTempNode(int lh, std::vector<int>& downPointers,
        std::vector<int>& edgeValues) { return 0; }

    /// Same as createTempNode(int, vector<int>) except this is for EV*MDDs.
    virtual int createTempNode(int lh, std::vector<int>& downPointers,
        std::vector<float>& edgeValues) { return 0; }


    /// Apply reduction rule to the temporary node and finalize it. Once
    /// a node is reduced, its contents cannot be modified.
    virtual int reduceNode(int node) = 0;
    // Dummy version available here.
    // In mdds_ext.h

    /// Reduce and finalize an node with an incoming edge value
    virtual void normalizeAndReduceNode(int& node, int& ev) = 0;
    // Dummy version available here.
    // In mdds_ext.h

    // Copy the downpointers into dptrs where dptrs[i] corresponds
    // to the downpointer at index i. For sparse nodes,
    // dptrs[index(i)] corresponds to the downpointer at index i.
    virtual bool getDownPtrs(int node, std::vector<int>& dptrs) const;

    // Similar to getDownPtrs() but for EV+MDDs
    virtual bool getDownPtrsAndEdgeValues(int node,
        std::vector<int>& dptrs, std::vector<int>& evs) const
    { return false; }

    // Similar to getDownPtrs() but for EV*MDDs
    virtual bool getDownPtrsAndEdgeValues(int node,
        std::vector<int>& dptrs, std::vector<float>& evs) const
    { return false; }

    /// Has the node been reduced
    bool isReducedNode(int node) const;

    bool isValidNodeIndex(int node) const;
    void reclaimOrphanNode(int node);     // for linkNode()
    void handleNewOrphanNode(int node);   // for unlinkNode()
    void deleteOrphanNode(int node);      // for uncacheNode()
    void freeZombieNode(int node);        // for uncacheNode()

    bool isStale(int h) const;

    void showNode(FILE *s, int p, int verbose = 0) const;
    void showNodeGraph(FILE *s, int p) const;
    void showInfo(FILE* strm, int verbosity);
    double getCardinality(int p) const;
    unsigned getNodeCount(int p) const;
    unsigned getEdgeCount(int p, bool countZeroes) const;

    long getCurrentNumNodes() const;
    long getCurrentMemoryUsed() const;
    long getCurrentMemoryAllocated() const;
    long getPeakNumNodes() const;
    long getPeakMemoryUsed() const;
    long getPeakMemoryAllocated() const;

    // Compaction threshold is a percentage value.
    // To set compaction threshold to 45%, call setCompactionThreshold(45).
    // Compaction will occur if
    // level[i].hole_slots > (level[i].size * compactionThreshold / 100) 
    unsigned getCompactionThreshold() const;
    forest::error setCompactionThreshold(unsigned p);
    forest::error compactMemory();
    forest::error garbageCollect();

    // *************** override expert_forest class -- done ***************

    bool areHolesRecycled() const;
    void setHoleRecycling(bool policy);
    int sharedCopy(int p);

    // use this to find out which level the node maps to
    int getNodeLevelMapping(int p) const;

    // Dealing with slot 2 (node size)
    int getLargestIndex(int p) const;

    // Dealing with entries

    // full node: entries start in the 4th slot (location 3, counting from 0)
    int* getFullNodeDownPtrs(int p);
    const int* getFullNodeDownPtrsReadOnly(int p) const;
    const int* getSparseNodeIndexes(int p) const;
    const int* getSparseNodeDownPtrs(int p) const;
    int getSparseNodeLargestIndex(int p) const;

    void setAllDownPtrs(int p, int value);
    void setAllDownPtrsWoUnlink(int p, int value);
    void initDownPtrs(int p);

    // for EVMDDs
    int* getFullNodeEdgeValues(int p);
    const int* getFullNodeEdgeValuesReadOnly(int p) const;
    const int* getSparseNodeEdgeValues(int p) const;
    void setAllEdgeValues(int p, int value);
    void setAllEdgeValues(int p, float fvalue);

    // p: node
    // i: the ith downpointer.
    // note: for sparse nodes this may not be the same as the ith index pointer.
    int getDownPtr(int p, int i) const;
    int getDownPtrAfterIndex(int p, int i, int &index) const;

    int getMddLevelMaxBound(int k) const;
    int getMxdLevelMaxBound(int k) const;
    int getLevelMaxBound(int k) const;

    bool isPrimedNode(int p) const;
    bool isUnprimedNode(int p) const;
    int buildQuasiReducedNodeAtLevel(int k, int p);

    void showLevel(FILE *s, int k) const;
    void showAll(FILE *s, int verb) const;

    void showNode(int p) const;
    void showAll() const;

    void reportMemoryUsage(FILE * s, const char filler=' ');

    void compareCacheCounts(int p = -1);
    void validateIncounts();


    // Remove zombies if more than max
    void removeZombies(int max = 100);

    // For uniqueness table
    int getNull() const;
    int getNext(int h) const;
    void setNext(int h, int n);
    // void cacheNode(int p);
    // void uncacheNode(int p);
    void show(FILE *s, int h) const;
#if ALT_HASH_CALL
    unsigned hash(int h) const;
#else
    unsigned hash(int h, unsigned M) const;
#endif
    bool equals(int h1, int h2) const;

    bool isCounting();

  protected:

    // Building level nodes
    int getLevelNode(int lh) const;
    int buildLevelNodeHelper(int lh, int* terminalNodes, int sz);
    void buildLevelNode(int lh, int* terminalNodes, int sz);
    void clearLevelNode(int lh);
    void clearLevelNodes();
    void clearAllNodes();

    // Building custom level nodes
    int* getTerminalNodes(int n);
    int* getTerminalNodes(int n, bool* terms);
    int* getTerminalNodes(int n, int* terms);
    int* getTerminalNodes(int n, float* terms);

    bool isValidVariable(int vh) const;
    bool doesLevelNeedCompaction(int k);

    // Dealing with node addressing
    void setNodeOffset(int p, int offset);

    // Dealing with node status
    bool isDeletedNode(int p) const;

    // Debug output
    void dump(FILE *s) const; 
    void dumpInternal(FILE *s) const; 
    void dumpInternalLevel(FILE *s, int k) const; 

    int setLevelBoundsAndHeights();
    int setLevelBoundAndHeight(int k, int sz, int h);

    long getUniqueTableMemoryUsed() const;
    long getTempNodeCount() const;
    long getZombieNodeCount() const;
    long getOrphanNodeCount() const;
    void updateMemoryAllocated(long bytes);
    long getHoleMemoryUsage() const;

    int getMaxHoleChain() const;
    int getCompactionsCount() const;

    // level based operations
    /// number of levels in the current mdd
    int getLevelCount() const;
    
    /// Move nodes so that all holes are at the end.
    void compactAllLevels();
    void compactLevel(int k);

    /// garbage collect
    bool gc();
    bool isTimeToGc();

    // zombify node p
    void zombifyNode(int p);

    // delete node p
    void deleteNode(int p);
    void deleteTempNode(int p);

    // free node p
    void freeNode(int p);

    // find the next free node in address[]
    int getFreeNode(int k);

    // returns offset to the hole found in level
    int getHole(int k, int slots, bool search_holes);

    // makes a hole of size == slots, at the specified level and offset
    void makeHole(int k, int p_offset, int slots);

    // int p in these functions is p's offset in data
    // change it to mean the mdd node handle

    // add a hole to the hole grid
    void gridInsert(int k, int p_offset);

    bool isHoleNonIndex(int k, int p_offset);

    // remove a non-index hole from the hole grid
    void midRemove(int k, int p_offset);

    // remove an index hole from the hole grid
    void indexRemove(int k, int p_offset);

  protected:
  
    // modify temp nodes count for level k as well as the global count
    void incrTempNodeCount(int k);
    void decrTempNodeCount(int k);
    
    // increment the count for "nodes activated since last garbage collection"
    void incrNodesActivatedSinceGc();

    // find, insert nodes into the unique table
    int find(int node);
    int insert(int node);

    // are sparse nodes enabled
    bool areSparseNodesEnabled() const;

    // is k a level that has been initialized
    bool isValidLevel(int k) const;

    // get the id used to indicate temporary nodes
    int getTempNodeId() const;

    // Cardinality for Mdds and MtMdds
    // k is the height (>= 0)
    double cardinality(int p, int k, std::map<int, double>& visited) const;

    // Cardinality for Mxds and MtMxds
    // k is height (always >=0). k combined with primeLevel indicates
    // the current primed or unprimed level.
    double cardinalityForRelations(int p, int k, bool primeLevel,
        std::map<int, double>& visited) const;

    // Pointer to expert_domain
    expert_domain* expertDomain;

    // Special next values
    static const int temp_node = -5;
    static const int non_index_hole = -2;

    /// Should we try to recycle holes.
    bool holeRecycling;
    
    /**
      Number of hole slots that trigger a compaction.
      This is a number between 0 and 1 indicating percentage.
      Compaction will occur if
      level[i].hole_slots > (level[i].size * compactionThreshold) 
      */
    float compactionThreshold;

    /// Size of address/next array.
    int a_size;
    /// Last used address.
    int a_last;
    /// Pointer to unused address list.
    int a_unused;
    /// Peak nodes;
    int peak_nodes;

    /// Number of levels. This is not the size of all the levels put together.
    int l_size;
    
    // performance stats

    /// For peak memory.
    long max_slots;
    /// Also for peak memory.
    long curr_slots;
    /// Number of alive nodes.
    long active_nodes;
    /// Largest traversed height of holes grid
    int max_hole_chain;
    /// Number of zombie nodes
    long zombie_nodes;
    /// These are just like zombies but they have not been zombified --
    /// exist only in non-pessimistic caches
    long orphan_nodes;
    /// Number of temporary nodes -- nodes that have not been reduced
    long temp_nodes;
    /// Number reclaimed nodes
    long reclaimed_nodes;
    /// Total number of compactions
    int num_compactions;
    /// Switch to turn on garbage collection; should be turned off when done
    bool enable_garbageCollection;
    /// Count of nodes created since last gc
    unsigned nodes_activated_since_gc;

    // Garbage collection in progress
    bool performing_gc;

    // Deleting terminal nodes (used in isStale() -- this enables
    // the removal of compute cache entries which refer to terminal nodes)
    bool delete_terminal_nodes;

    /// Uniqueness table
    mdd_hash_table <node_manager> *unique;

    long curr_mem_alloc;
    long max_mem_alloc;

    bool counting;

    // scratch pad for buildLevelNode and getTerminalNodes
    int* dptrs;
    int dptrsSize;
};


/// Inline functions implemented here

inline bool node_manager::isValidNodeIndex(int node) const {
  return node <= a_last;
}

inline void node_manager::reclaimOrphanNode(int p) {
  DCASSERT(!isPessimistic() || !isZombieNode(p));
  DCASSERT(isActiveNode(p));
  DCASSERT(!isTerminalNode(p));
  DCASSERT(isReducedNode(p));
  DCASSERT(getInCount(p) == 0);

#if 0
  validateIncounts();
#endif

  // increase incount
  ++(*(getNodeAddress(p)));

  reclaimed_nodes++;
  orphan_nodes--;
  
#ifdef TRACK_DELETIONS
  cout << "\t+Node " << p << " count now " << *(getNodeAddress(p)) << "\n";
  cout.flush();
#endif
}  

inline void node_manager::handleNewOrphanNode(int p) {
  DCASSERT(!isPessimistic() || !isZombieNode(p));
  DCASSERT(isActiveNode(p));
  DCASSERT(!isTerminalNode(p));
  DCASSERT(isReducedNode(p));
  DCASSERT(getInCount(p) == 0);

  // insted of orphan_nodes++ here; do it only when the orphan is not going
  // to get deleted or converted into a zombie

  if (getCacheCount(p) == 0) {
    // delete node
#ifdef TRACK_DELETIONS
    cout << "Deleting node " << p << " from unlinkNode\t";
    showNode(stdout, p);
    cout << "\n";
    cout.flush();
#endif
    deleteNode(p);
  } else if (isPessimistic()) {
    // zombify node
    zombifyNode(p);
  } else {
    orphan_nodes++;
  }

#if 0
  if (getOrphanNodeCount() > 100000)
    smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())
      ->removeStales();
#endif
}

inline void node_manager::deleteOrphanNode(int p) {
  DCASSERT(!isPessimistic());
  DCASSERT(getCacheCount(p) == 0 && getInCount(p) == 0);
#ifdef TRACK_DELETIONS
  cout << "Deleting node " << p << " from uncacheNode\t";
  showNode(stdout, p);
  cout << "\n";
  cout.flush();
#endif
  orphan_nodes--;
  deleteNode(p);
}

inline int node_manager::getDownPtrAfterIndex(int p, int i, int &index)
  const {
  DCASSERT(isActiveNode(p));
  DCASSERT(i >= 0);
  if (isTerminalNode(p)) return p;
  DCASSERT(i < getLevelSize(getNodeLevel(p)));
  if (isFullNode(p)) {
    // full or trunc-full node
    // full or trunc-full node, but i lies after the last non-zero downpointer
    // index = i;
    return (i < getFullNodeSize(p))? getFullNodeDownPtr(p, i): 0;
  } else {
    // sparse node
    // binary search to find the index ptr corresponding to i
    int stop = getSparseNodeSize(p);
    while (index < stop && i > getSparseNodeIndex(p, index)) index++;
    return (index < stop && i == getSparseNodeIndex(p, index))?
        getSparseNodeDownPtr(p, index): 0;
  }
}


inline int node_manager::getDownPtr(int p, int i) const {
  DCASSERT(isActiveNode(p));
  if (isTerminalNode(p)) return p;
  DCASSERT(i >= 0);
  if (isFullNode(p)) {
    // full or trunc-full node
    if (getFullNodeSize(p) > i) return getFullNodeDownPtr(p, i);

    // full or trunc-full node, but i lies after the last non-zero downpointer
    return 0;
  } else {
    // sparse node
    // binary search to find the index ptr corresponding to i
    int start = 0;
    int stop = getSparseNodeSize(p) - 1;
    
    // the index ptr corresponding to i has been compressed i.e. its a zero.
    if (getSparseNodeIndex(p, start) > i) return 0;
    if (getSparseNodeIndex(p, stop) < i) return 0;

    int mid = (start + stop)/2;
    while(start < stop) {
      if (getSparseNodeIndex(p, mid) == i)
        return getSparseNodeDownPtr(p, mid);
      if (getSparseNodeIndex(p, mid) > i)
        stop = mid - 1;
      else
        start = mid + 1;
      mid = (start + stop)/2;
    }
    if (getSparseNodeIndex(p, mid) == i)
      return getSparseNodeDownPtr(p, mid);
   
    // index ptr not found - it was compressed because it was a zero.
    return 0;
  }
}


#if 0
inline int node_manager::createVarNode(int lh)
{
  int node = level[mapLevel(k)].varNode;
  linkNode(node);
  return node;
}
#endif

inline int node_manager::createTempNode(int k, int sz, bool clear)
{
  DCASSERT(k != 0);

  if (isTimeToGc()) { gc(); }

  CHECK_RANGE(1, mapLevel(k), l_size);
  DCASSERT(level[mapLevel(k)].data != NULL);
  CHECK_RANGE(1, sz, getLevelSize(k) + 1);

  // get a location in address[] to store the node
  int p = getFreeNode(k);

#ifdef DEBUG_MDD_H
  printf("%s: k: %d, sz: %d, new p: %d\n", __func__, k, sz, p);
  fflush(stdout);
#endif

  // fill in the location with p's address info
  DCASSERT(getEdgeLabeling() == forest::MULTI_TERMINAL);
  address[p].level = k;
  address[p].offset = getHole(k, 4 + sz, true);
  address[p].cache_count = 0;

#ifdef DEBUG_MDD_H
  printf("%s: offset: %d\n", __func__, address[p].offset);
  fflush(stdout);
#endif

  int* foo = level[mapLevel(k)].data + address[p].offset;
  foo[0] = 1;                   // #incoming
  foo[1] = temp_node;
  foo[2] = sz;                  // size
  foo[3 + sz] = p;              // pointer to this node in the address array

  // initialize
  if (clear) initDownPtrs(p);

#ifdef TRACK_DELETIONS
  cout << "Creating node " << p << "\n";
  cout.flush();
#endif

  incrTempNodeCount(k);
  nodes_activated_since_gc++;

  return p;
}


inline
int
node_manager::createTempNode(int lh, std::vector<int>& downPointers)
{
  int tempNode = createTempNode(lh, downPointers.size(), false);
  int* dptrs = getFullNodeDownPtrs(tempNode);
  std::vector<int>::iterator iter = downPointers.begin();
  while (iter != downPointers.end())
  {
    *dptrs++ = *iter++;
  }
  return tempNode;
}


// TODO: HERE: should we move these next few inlines up to mddinternal?
inline bool node_manager::areHolesRecycled() const {
  return holeRecycling;
}
inline unsigned node_manager::getCompactionThreshold() const {
  return unsigned(compactionThreshold * 100.0);
}

inline forest::error node_manager::setCompactionThreshold(unsigned t) {
  if (t > 100) return forest::INVALID_ASSIGNMENT;
  compactionThreshold = t/100.0;
  return forest::SUCCESS;
}

// Dealing with cache count
// Dealing with node level

// use this to find out which level the node maps to
inline int node_manager::getNodeLevelMapping(int p) const {
  return mapLevel(getNodeLevel(p));
}

// Dealing with slot 0 (incount)

// linkNodeing and unlinkNodeing nodes

inline int node_manager::sharedCopy(int p) {
  linkNode(p);
  return p;
}

// Dealing with slot 1 (next pointer)

inline bool node_manager::isReducedNode(int p) const {
#ifdef DEBUG_MDD_H
  printf("%s: p: %d\n", __func__, p);
#endif
  DCASSERT(isActiveNode(p));
  return (isTerminalNode(p) || (getNext(p) >= getNull()));
}

// Dealing with slot 2 (node size)
inline int node_manager::getLargestIndex(int p) const {
  DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return isFullNode(p)? getFullNodeSize(p) - 1: getSparseNodeLargestIndex(p);
}

// Dealing with entries

inline
bool node_manager::getDownPtrs(int p, std::vector<int>& dptrs) const {
  if (!isActiveNode(p) || isTerminalNode(p) || !isReducedNode(p))
    return false;

  if (isFullNode(p)) {
    int size = getFullNodeSize(p);
    if (dptrs.size() < unsigned(size))
      dptrs.resize(size, 0);
    const int* ptrs = getFullNodeDownPtrsReadOnly(p);
    const int* end = ptrs + size;
    std::vector<int>::iterator iter = dptrs.begin();
    while (ptrs != end)
    {
      *iter++ = *ptrs++;
    }
  }
  else {
    int nnz = getSparseNodeSize(p);
    int size = getLargestIndex(p) + 1;
    if (dptrs.size() < unsigned(size))
      dptrs.resize(size, 0);
    const int* ptrs = getSparseNodeDownPtrs(p);
    const int* index = getSparseNodeIndexes(p);
    const int* end = ptrs + nnz;
    while (ptrs != end)
    {
      dptrs[*index++] = *ptrs++;
    }
  }
  return true;
}

// full node: entries start in the 4th slot (location 3, counting from 0)
inline int* node_manager::getFullNodeDownPtrs(int p) {
#ifdef DEBUG_MDD_H
  printf("%s: p: %d\n", __func__, p);
#endif
  DCASSERT(isFullNode(p));
  DCASSERT(!isReducedNode(p));
  return (getNodeAddress(p) + 3);
}

inline const int* node_manager::getFullNodeDownPtrsReadOnly(int p) const {
  DCASSERT(isFullNode(p));
  return (getNodeAddress(p) + 3);
}

inline const int* node_manager::getFullNodeEdgeValuesReadOnly(int p) const {
  DCASSERT(isFullNode(p));
  return (getNodeAddress(p) + 3 + getFullNodeSize(p));
}

inline int* node_manager::getFullNodeEdgeValues(int p) {
  DCASSERT(isFullNode(p));
  DCASSERT(!isReducedNode(p));
  return (getNodeAddress(p) + 3 + getFullNodeSize(p));
}

inline const int* node_manager::getSparseNodeIndexes(int p) const {
  DCASSERT(isSparseNode(p));
  return (getNodeAddress(p) + 3);
}

inline const int* node_manager::getSparseNodeDownPtrs(int p) const {
  DCASSERT(isSparseNode(p));
  return (getNodeAddress(p) + 3 + getSparseNodeSize(p));
}

inline const int* node_manager::getSparseNodeEdgeValues(int p) const {
  DCASSERT(isSparseNode(p));
  return (getNodeAddress(p) + 3 + 2 * getSparseNodeSize(p));
}

inline int node_manager::getSparseNodeLargestIndex(int p) const {
  DCASSERT(isSparseNode(p));
  return getSparseNodeIndex(p, getSparseNodeSize(p) - 1);
}

inline void node_manager::setAllDownPtrs(int p, int value) {
  DCASSERT(!isReducedNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(isReducedNode(value));
  int* curr = getFullNodeDownPtrs(p);
  int size = getFullNodeSize(p);
  for (int* end = curr + size; curr != end; )
  {
    unlinkNode(*curr);
    *curr++ = value;
  }
  if (!isTerminalNode(value)) getInCount(value) += size;
}

inline void node_manager::setAllDownPtrsWoUnlink(int p, int value) {
  DCASSERT(!isReducedNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(isReducedNode(value));
  int* curr = getFullNodeDownPtrs(p);
  int size = getFullNodeSize(p);
  for (int* end = curr + size; curr != end; )
  {
    *curr++ = value;
  }
  if (!isTerminalNode(value)) getInCount(value) += size;
}

inline void node_manager::initDownPtrs(int p) {
  DCASSERT(!isReducedNode(p));
  DCASSERT(isFullNode(p));
  memset(getFullNodeDownPtrs(p), 0, sizeof(int) * getFullNodeSize(p));
}

inline void node_manager::setAllEdgeValues(int p, int value) {
  DCASSERT(edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES);
  DCASSERT(!isReducedNode(p));
  DCASSERT(isFullNode(p));
  int *edgeptr = getFullNodeEdgeValues(p);
  int *last = edgeptr + getFullNodeSize(p);
  for ( ; edgeptr != last; ++edgeptr) *edgeptr = value;
}


inline void node_manager::setAllEdgeValues(int p, float fvalue) {
  DCASSERT(edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES);
  DCASSERT(!isReducedNode(p));
  DCASSERT(isFullNode(p));
  int *edgeptr = getFullNodeEdgeValues(p);
  int *last = edgeptr + getFullNodeSize(p);
  int value = toInt(fvalue);
  for ( ; edgeptr != last; ++edgeptr) *edgeptr = value;
}


inline bool node_manager::isPrimedNode(int p) const {
  return (getNodeLevel(p) < 0);
}
inline bool node_manager::isUnprimedNode(int p) const {
  return (getNodeLevel(p) > 0);
}

// For uniqueness table
inline int node_manager::getNull() const { return -1; }
inline int node_manager::getNext(int h) const { 
  DCASSERT(isActiveNode(h));
  DCASSERT(!isTerminalNode(h));
  // next pointer is at slot 1 (counting from 0)
  return *(getNodeAddress(h) + 1);
}
inline void node_manager::setNext(int h, int n) { 
  DCASSERT(isActiveNode(h));
  DCASSERT(!isTerminalNode(h));
  *(getNodeAddress(h) + 1) = n; 
}
inline bool node_manager::isStale(int h) const {
#if 1
  return
    isTerminalNode(h)?
      delete_terminal_nodes:
      isPessimistic()?
        isZombieNode(h):
        (enable_garbageCollection && getInCount(h) == 0);
#else
  return isTerminalNode(p)? false: isZombieNode(h);
#endif
}
inline bool node_manager::isCounting() { return counting; }

// Dealing with node addressing

inline void node_manager::setNodeOffset(int p, int offset) {
  CHECK_RANGE(1, p, a_last+1);
  address[p].offset = offset;
}

// Dealing with node status

inline bool node_manager::isDeletedNode(int p) const {
  CHECK_RANGE(1, p, a_last+1);
  return !(isActiveNode(p) || isZombieNode(p));
}

inline long node_manager::getUniqueTableMemoryUsed() const {
  return (unique->getSize() * sizeof(int));
}
inline long node_manager::getPeakNumNodes() const {
  return peak_nodes;    // terminal nodes not included
}
inline long node_manager::getCurrentNumNodes() const {
#if 1
  return active_nodes; // excludes terminal nodes
#else
  return active_nodes - orphan_nodes - temp_nodes;
#endif
}
inline long node_manager::getTempNodeCount() const {
  return temp_nodes;
}
inline long node_manager::getZombieNodeCount() const {
  return zombie_nodes;
}
inline long node_manager::getOrphanNodeCount() const {
  return orphan_nodes;
}
inline long node_manager::getPeakMemoryUsed() const {
#if 0
  int sum = 0;
  for (int i = 0; i < l_size; i++) sum += level[i].max_slots;
  return sum * sizeof(int);
#else
  return max_slots * sizeof(int);
#endif
}

inline void node_manager::updateMemoryAllocated(long bytes) {
  curr_mem_alloc += bytes;
  if (curr_mem_alloc > max_mem_alloc) max_mem_alloc = curr_mem_alloc;
#if 0
  printf("%p: Curr: %d, Peak: %d\n", this, curr_mem_alloc, max_mem_alloc);
#endif
}

/// number of levels in the current mdd
inline int node_manager::getLevelCount() const { return l_size; }

/// garbage collect
inline bool node_manager::isTimeToGc()
{
#if 1
  return (isPessimistic())? (getZombieNodeCount() > 1000000):
      (getOrphanNodeCount() > 500000);
  // return (nodes_activated_since_gc > 5000000); //10,000,000
#else
  return false;
#endif
}

inline bool node_manager::isHoleNonIndex(int k, int p_offset) {
  return (level[mapLevel(k)].data[p_offset + 1] == non_index_hole);
}

inline bool node_manager::doesLevelNeedCompaction(int k)
{
#if 0
  return (level[mapLevel(k)].hole_slots >
      (level[mapLevel(k)].size * compactionThreshold));
#else
  return ((level[mapLevel(k)].hole_slots > 10000) ||
      ((level[mapLevel(k)].hole_slots > 100) && 
       (level[mapLevel(k)].hole_slots >
        (level[mapLevel(k)].last * compactionThreshold))));
#endif
}

inline void node_manager::midRemove(int k, int p_offset) {
  DCASSERT(isHoleNonIndex(k, p_offset));
  int p_level = mapLevel(k);
  int left = level[p_level].data[p_offset+2];
  DCASSERT(left);
  int right = level[p_level].data[p_offset+3];

  level[p_level].data[left + 3] = right;
  if (right) level[p_level].data[right + 2] = left;
}

inline void node_manager::incrTempNodeCount(int k) {
  level[mapLevel(k)].temp_nodes++;
  temp_nodes++;
}


inline void node_manager::decrTempNodeCount(int k) {
  level[mapLevel(k)].temp_nodes--;
  temp_nodes--;
}


inline void node_manager::incrNodesActivatedSinceGc() {
  nodes_activated_since_gc++;
}


inline int node_manager::find(int node) {
  return unique->find(node);
}


inline int node_manager::insert(int node) {
  return unique->insert(node);
}


inline bool node_manager::areSparseNodesEnabled() const {
  return nodeStorage == forest::FULL_OR_SPARSE_STORAGE ||
      nodeStorage == forest::SPARSE_STORAGE;
}


inline bool node_manager::isValidLevel(int k) const {
  int mapped_level = mapLevel(k);
  return (1 <= mapped_level && mapped_level < l_size &&
    level[mapped_level].data != NULL);
}


inline int node_manager::getTempNodeId() const {
  return temp_node;
}

// ****************** override expert_forest class ****************** 


inline int node_manager::reduceNode(int node) {
  fprintf(stderr, "Error:\n");
  fprintf(stderr, "  Called reduceNode() for a forest with edge-values.\n");
  fprintf(stderr, "  Call normalizeAndReduceNode() instead.\n");
  fflush(stderr);
  return 0;
}

inline void node_manager::normalizeAndReduceNode(int& node, int& ev) {
  fprintf(stderr, "Error:\n");
  fprintf(stderr, "  Called normalizeAndReduceNode() for a forest\n");
  fprintf(stderr, "  without edge-values. Call reduceNode() instead.\n");
  fflush(stderr);
  node = 0;
  ev = INF;
}

inline forest::error node_manager::compactMemory() {
  compactAllLevels();
  return forest::SUCCESS;
}

inline long node_manager::getCurrentMemoryAllocated() const {
  return curr_mem_alloc;
}

inline long node_manager::getPeakMemoryAllocated() const {
  return max_mem_alloc;
}

inline void node_manager::showInfo(FILE* strm, int verbosity) {
  showAll(strm, verbosity);
  fprintf(strm, "DD stats:\n");
  reportMemoryUsage(strm, '\t');
}

// *************** override expert_forest class -- done ***************

inline int node_manager::getLevelNode(int k) const {
  return level[mapLevel(k)].levelNode;
}

inline bool node_manager::isValidVariable(int vh) const {
  return expertDomain->getVariableHeight(vh) != -1;
}

inline forest::error
node_manager::findFirstElement(const dd_edge& f, int* vlist) const
{
  return forest::INVALID_OPERATION;
}

inline forest::error
node_manager::findFirstElement(const dd_edge& f, int* vlist, int* vplist) const
{
  return forest::INVALID_OPERATION;
}

#endif

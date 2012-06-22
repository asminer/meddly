
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


// TODO: go through this file & cleanup

//TODO: add a mechanism to mt_forest so that reduction rule can
//      be set after instantiation of the mt_forest.

#ifndef MT_FOREST
#define MT_FOREST

#include <fstream>
#include <iostream>
#include <vector>

#include "../defines.h"
#include "../mdd_hash.h"

namespace MEDDLY {
  class mt_forest;
};


/*
 * The MDD node handle is an integer.
 * The MDD node handle represents the offset in address[] where the node data
 * is stored.
 */


/**

  Node states: Active, Hole.

  Active Node:
  ------------

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
  [3..2+Size]       Downward pointers (>=0)
  [3+Size]          index ptr to this node in address array (>=0)

  Sparse storage:
  [3..2+nnz]        Indexes (>=0)
  [3+nnz..2+2*nnz]  Downward pointers (>=0)
  [3+2*nnz]         index ptr to this node in address array (>=0)


  Hole:
  -----

  When deleted, the node becomes a hole.
  A hole must satisfy hole[0] == -(size_of_hole) == hole[-(size_of_hole)-1]
  i.e. the first and last locations in hole[] store -(size_of_hole).

  There are two kinds oh holes depending on their location in the "hole grid":
  Index Holes and Non-index Holes.

  The hole grid structure:
  ------------------------
  (holes_bottom)
  holes_of_size_0 (index) -- holes_of_size_0 (non_index) -- (non_index) -- NULL
  |
  holes_of_size_1 -- ""
  |
  :
  :
  (holes_top)

  Index holes are represented as follows:
  ---------------------------------------
  [0] -size (number of slots in hole)     
  [1] up
  [2] down 
  [3] next pointer (nodes of same size)
  [4..size-2] Unused
  :
  :
  [size-1] -size

  Non-index holes are represented as follows:
  [0] -size (number of slots in hole)     
  [1] flag (<0, indicates non-index node)
  [2] prev pointer (nodes of same size)
  [3] next pointer (nodes of same size)
  [4..size-2] Unused
  :
  :
  [size-1] -size
*/


class MEDDLY::mt_forest : public expert_forest {
  public:
    mt_forest(int dsl, domain *d, bool rel, range_type t, edge_labeling ev, 
      const policies &p);
    virtual ~mt_forest();

  // ------------------------------------------------------------
  // |
  // |  Required for uniqueness hash table
  // |
  public:
    inline int getNull() const { return -1; }
    inline int getNext(int h) const {
      MEDDLY_DCASSERT(isActiveNode(h));
      MEDDLY_DCASSERT(!isTerminalNode(h));
      // next pointer is at slot 1 (counting from 0)
      MEDDLY_DCASSERT(getNodeAddress(h));
      return *(getNodeAddress(h) + level_data::next_index);
    }

    inline void setNext(int h, int n) {
      MEDDLY_DCASSERT(isActiveNode(h));
      MEDDLY_DCASSERT(!isTerminalNode(h));
      *(getNodeAddress(h) + level_data::next_index) = n; 
    }

    // For uniqueness table
    // void cacheNode(int p);
    // void uncacheNode(int p);
    void show(FILE *s, int h) const;
    unsigned hash(int h) const;
    bool equals(int h1, int h2) const;
    
    bool equalsFF(int h1, int h2) const;
    bool equalsSS(int h1, int h2) const;
    bool equalsFS(int h1, int h2) const;


  public:
    using expert_forest::getDownPtrsAndEdgeValues;
    using expert_forest::getSparseNodeIndexes;

    /// Refer to meddly.h
    void createEdgeForVar(int vh, bool primedLevel,
        bool* terms, dd_edge& a);
    void createEdgeForVar(int vh, bool primedLevel,
        int* terms, dd_edge& a);
    void createEdgeForVar(int vh, bool primedLevel,
        float* terms, dd_edge& a);

    void createSubMatrix(const bool* const* vlist,
        const bool* const* vplist, const dd_edge a, dd_edge& b);
    void createSubMatrix(const dd_edge& rows, const dd_edge& cols,
        const dd_edge& a, dd_edge& b);

    virtual void getElement(const dd_edge& a, int index, int* e);

    virtual void findFirstElement(const dd_edge& f, int* vlist) const;
    virtual void findFirstElement(const dd_edge& f, int* vlist,
        int* vplist) const;

    virtual void accumulate(int& a, int b);
    virtual bool accumulate(int& tempNode, int* element);

    // cBM: Copy before modifying.
    virtual int accumulateMdd(int a, int b, bool cBM);
    virtual int addReducedNodes(int a, int b);
    virtual int accumulateExpandA(int a, int b, bool cBM);
    int accumulate(int tempNode, bool cBM, int* element, int level);
    virtual int makeACopy(int node, int size = 0);

    /// Create a temporary node -- a node that can be modified by the user
    virtual int createTempNode(int lh, int size, bool clear);

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

    // Helpers for reduceNode().
    // These method assume that the top leve node is a temporary node
    // whose children may be either reduced or temporary nodes
    // (with an incount >= 1).
    // Since there is a possibility of a temporary node being referred
    // to multiple times, these methods use a cache to ensure that each
    // temporary node is reduced only once.
    // Note: The same cache can be used across consecutive reduce
    // operations by specifying clearCache to false.
    int recursiveReduceNode(int tempNode, bool clearCache = true);
    int recursiveReduceNode(std::map<int, int>& cache, int root);

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

    virtual void dumpUniqueTable(FILE* s) const;
    void reclaimOrphanNode(int node);     // for linkNode()
    void handleNewOrphanNode(int node);   // for unlinkNode()
    void deleteOrphanNode(int node);      // for uncacheNode()
    void freeZombieNode(int node);        // for uncacheNode()

    bool discardTemporaryNodesFromComputeCache() const;   // for isStale()

    void showNode(FILE *s, int p, int verbose = 0) const;
    void showNodeGraph(FILE *s, int p) const;
    void showInfo(FILE* strm, int verbosity);

    // Compaction threshold is a percentage value.
    // To set compaction threshold to 45%, call setCompactionThreshold(45).
    // Compaction will occur if
    // level[i].hole_slots > (level[i].size * compactionThreshold / 100) 
    // unsigned getCompactionThreshold() const;
    // void setCompactionThreshold(unsigned p);
    void compactMemory();
    void garbageCollect();

    // *************** override expert_forest class -- done ***************

    // bool areHolesRecycled() const;
    // void setHoleRecycling(bool policy);
    int sharedCopy(int p);

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

    // for EVMDDs
    int* getFullNodeEdgeValues(int p);
    const int* getFullNodeEdgeValuesReadOnly(int p) const;
    const int* getSparseNodeEdgeValues(int p) const;
    void setAllEdgeValues(int p, int value);
    void setAllEdgeValues(int p, float fvalue);

    // p: node
    // i: the ith downpointer.
    // note: for sparse nodes this may not be the same as the ith index pointer.
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
    // int* getTerminalNodes(int n);
    int* getTerminalNodes(int n, bool* terms);
    int* getTerminalNodes(int n, int* terms);
    int* getTerminalNodes(int n, float* terms);

    bool isValidVariable(int vh) const;
    bool doesLevelNeedCompaction(int k);

    // Dealing with node addressing
    void setNodeOffset(int p, int offset);

    // Dealing with node status
    bool isDeletedNode(int p) const;

    long getUniqueTableMemoryUsed() const;
    long getHoleMemoryUsage() const;

    /// Move nodes so that all holes are at the end.
    void compactAllLevels();

    /// garbage collect
    bool gc(bool zombifyOrphanNodes = false);
    bool isTimeToGc();

    // zombify node p
    void zombifyNode(int p);

    // delete node p
    void deleteNode(int p);

    // free node p
    void freeNode(int p);

  protected:

    // modify temp nodes count for level k as well as the global count
    // the temporary node count should be incremented only within
    // createTempNode() or variants.
    // decrTempNodeCount() should be called by any method that changes a
    // temporary node to a reduced node.
    // Note: deleting a temp node automatically calls decrTempNodeCount().
    void incrTempNodeCount(int k);
    void decrTempNodeCount(int k);

    // increment the count for "nodes activated since last garbage collection"
    void incrNodesActivatedSinceGc();

    // find, insert nodes into the unique table
    int find(int node);
    int insert(int node);
    int replace(int node);

    // get the id used to indicate temporary nodes
    int getTempNodeId() const;

    // Checks if one of the following reductions is satisfied:
    // Fully, Quasi, Identity Reduced.
    // If the node can be reduced to an existing node, the existing node
    // is returned.
    bool checkForReductions(int p, int nnz, int& result);

    // Checks if the node has a single downpointer enabled and at
    // the given index.
    bool singleNonZeroAt(int p, int val, int index) const;

    // Checks if the node satisfies the forests reduction rules.
    // If it does not, an assert violation occurs.
    void validateDownPointers(int p, bool recursive = false);

    // Pointer to expert_domain
    // expert_domain* expertDomain;

    // Special next values
    static const int temp_node = -5;

    // performance stats

    /// Count of nodes created since last gc
    unsigned nodes_activated_since_gc;

    // Garbage collection in progress
    bool performing_gc;

    // Deleting terminal nodes (used in isStale() -- this enables
    // the removal of compute cache entries which refer to terminal nodes)
    bool delete_terminal_nodes;

    /// Uniqueness table
    mdd_hash_table <mt_forest> *unique;

    // long curr_mem_alloc;
    // long max_mem_alloc;

    bool counting;

    // scratch pad for buildLevelNode and getTerminalNodes
    int* dptrs;
    int dptrsSize;

    // Place holder for accumulate-minterm result.
    bool accumulateMintermAddedElement;

  private:
    // Cache for recursiveReduceNode()
    std::map<int, int> recursiveReduceCache;

    // Persistant variables used in addReducedNodes()
    dd_edge* nodeA;
    dd_edge* nodeB;
};


/// Inline functions implemented here

inline void MEDDLY::mt_forest::reclaimOrphanNode(int p) {
  MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));
  MEDDLY_DCASSERT(isActiveNode(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(isReducedNode(p));
  stats.reclaimed_nodes++;
  stats.orphan_nodes--;
}  

inline void MEDDLY::mt_forest::deleteOrphanNode(int p) {
  MEDDLY_DCASSERT(!isPessimistic());
  MEDDLY_DCASSERT(getCacheCount(p) == 0 && getInCount(p) == 0);
#ifdef TRACK_DELETIONS
  cout << "Deleting node " << p << " from uncacheNode\t";
  showNode(stdout, p);
  cout << "\n";
  cout.flush();
#endif
  stats.orphan_nodes--;
  deleteNode(p);
}

inline int MEDDLY::mt_forest::getDownPtrAfterIndex(int p, int i, int &index)
  const {
  MEDDLY_DCASSERT(isActiveNode(p));
  MEDDLY_DCASSERT(i >= 0);
  if (isTerminalNode(p)) return p;
  MEDDLY_DCASSERT(i < getLevelSize(getNodeLevel(p)));
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


inline
int MEDDLY::mt_forest::createTempNode(int lh, std::vector<int>& downPointers)
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


// Dealing with slot 0 (incount)

// linkNodeing and unlinkNodeing nodes

inline int MEDDLY::mt_forest::sharedCopy(int p) {
  linkNode(p);
  return p;
}

// Dealing with slot 1 (next pointer)

inline bool MEDDLY::mt_forest::isReducedNode(int p) const {
#ifdef DEBUG_MDD_H
  printf("%s: p: %d\n", __func__, p);
#endif
  MEDDLY_DCASSERT(isActiveNode(p));
  return (isTerminalNode(p) || (getNext(p) >= getNull()));
}

// Dealing with slot 2 (node size)
inline int MEDDLY::mt_forest::getLargestIndex(int p) const {
  MEDDLY_DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return isFullNode(p)? getFullNodeSize(p) - 1: getSparseNodeLargestIndex(p);
}

// Dealing with entries

// full node: entries start in the 4th slot (location 3, counting from 0)
inline int* MEDDLY::mt_forest::getFullNodeDownPtrs(int p) {
#ifdef DEBUG_MDD_H
  printf("%s: p: %d\n", __func__, p);
#endif
  MEDDLY_DCASSERT(isFullNode(p));
  MEDDLY_DCASSERT(!isReducedNode(p));
  return (getNodeAddress(p) + 3);
}

inline const int* MEDDLY::mt_forest::getFullNodeDownPtrsReadOnly(int p) const {
  MEDDLY_DCASSERT(isFullNode(p));
  return (getNodeAddress(p) + 3);
}

inline const int* MEDDLY::mt_forest::getFullNodeEdgeValuesReadOnly(int p) const {
  MEDDLY_DCASSERT(isFullNode(p));
  return (getNodeAddress(p) + 3 + getFullNodeSize(p));
}

inline int* MEDDLY::mt_forest::getFullNodeEdgeValues(int p) {
  MEDDLY_DCASSERT(isFullNode(p));
  MEDDLY_DCASSERT(!isReducedNode(p));
  return (getNodeAddress(p) + 3 + getFullNodeSize(p));
}

inline const int* MEDDLY::mt_forest::getSparseNodeIndexes(int p) const {
  MEDDLY_DCASSERT(isSparseNode(p));
  return (getNodeAddress(p) + 3);
}

inline const int* MEDDLY::mt_forest::getSparseNodeDownPtrs(int p) const {
  MEDDLY_DCASSERT(isSparseNode(p));
  return (getNodeAddress(p) + 3 + getSparseNodeSize(p));
}

inline const int* MEDDLY::mt_forest::getSparseNodeEdgeValues(int p) const {
  MEDDLY_DCASSERT(isSparseNode(p));
  return (getNodeAddress(p) + 3 + 2 * getSparseNodeSize(p));
}

inline int MEDDLY::mt_forest::getSparseNodeLargestIndex(int p) const {
  MEDDLY_DCASSERT(isSparseNode(p));
  return getSparseNodeIndex(p, getSparseNodeSize(p) - 1);
}

inline void MEDDLY::mt_forest::setAllDownPtrs(int p, int value) {
  MEDDLY_DCASSERT(!isReducedNode(p));
  MEDDLY_DCASSERT(isFullNode(p));
  MEDDLY_DCASSERT(isActiveNode(value));
  int* curr = getFullNodeDownPtrs(p);
  int size = getFullNodeSize(p);
  for (int* end = curr + size; curr != end; )
  {
    unlinkNode(*curr);
    *curr++ = value;
  }
  if (!isTerminalNode(value)) getInCount(value) += size;
}

inline void MEDDLY::mt_forest::setAllDownPtrsWoUnlink(int p, int value) {
  MEDDLY_DCASSERT(!isReducedNode(p));
  MEDDLY_DCASSERT(isFullNode(p));
  MEDDLY_DCASSERT(isActiveNode(value));
  int* curr = getFullNodeDownPtrs(p);
  int size = getFullNodeSize(p);
  for (int* end = curr + size; curr != end; )
  {
    *curr++ = value;
  }
  if (!isTerminalNode(value)) getInCount(value) += size;
}

inline void MEDDLY::mt_forest::setAllEdgeValues(int p, int value) {
  MEDDLY_DCASSERT(isEVPlus() || isEVTimes());
  MEDDLY_DCASSERT(!isReducedNode(p));
  MEDDLY_DCASSERT(isFullNode(p));
  int *edgeptr = getFullNodeEdgeValues(p);
  int *last = edgeptr + getFullNodeSize(p);
  for ( ; edgeptr != last; ++edgeptr) *edgeptr = value;
}


inline void MEDDLY::mt_forest::setAllEdgeValues(int p, float fvalue) {
  MEDDLY_DCASSERT(isEVPlus() || isEVTimes());
  MEDDLY_DCASSERT(!isReducedNode(p));
  MEDDLY_DCASSERT(isFullNode(p));
  int *edgeptr = getFullNodeEdgeValues(p);
  int *last = edgeptr + getFullNodeSize(p);
  int value = toInt(fvalue);
  for ( ; edgeptr != last; ++edgeptr) *edgeptr = value;
}


inline bool MEDDLY::mt_forest::isPrimedNode(int p) const {
  return (getNodeLevel(p) < 0);
}
inline bool MEDDLY::mt_forest::isUnprimedNode(int p) const {
  return (getNodeLevel(p) > 0);
}

inline bool MEDDLY::mt_forest::discardTemporaryNodesFromComputeCache() const {
  return delete_terminal_nodes;
}

inline bool MEDDLY::mt_forest::isCounting() { return counting; }

// Dealing with node addressing

inline void MEDDLY::mt_forest::setNodeOffset(int p, int offset) 
{
  MEDDLY_DCASSERT(isValidNonterminalIndex(p));
  address[p].offset = offset;
}

// Dealing with node status

inline bool MEDDLY::mt_forest::isDeletedNode(int p) const 
{
  MEDDLY_DCASSERT(isValidNonterminalIndex(p));
  return !(isActiveNode(p) || isZombieNode(p));
}

inline long MEDDLY::mt_forest::getUniqueTableMemoryUsed() const {
  return (unique->getSize() * sizeof(int));
}

/// garbage collect
inline bool MEDDLY::mt_forest::isTimeToGc()
{
#if 1
  // const int zombieTrigger = 1000;  // use for debugging
  // const int orphanTrigger = 500;  // use for debugging
  // const int zombieTrigger = 1000000;
  // const int orphanTrigger = 500000;
  return isPessimistic() 
    ? (stats.zombie_nodes > deflt.zombieTrigger)
    : (stats.orphan_nodes > deflt.orphanTrigger);
#elif 0
  return false;
#else
  return isPessimistic()? false: (getOrphanNodeCount() > 500000);
#endif
}

inline void MEDDLY::mt_forest::incrTempNodeCount(int k) {
  levels[k].incrTempNodeCount();
  stats.temp_nodes++;
}


inline void MEDDLY::mt_forest::decrTempNodeCount(int k) {
  levels[k].decrTempNodeCount();
  stats.temp_nodes--;
}


inline void MEDDLY::mt_forest::incrNodesActivatedSinceGc() {
  nodes_activated_since_gc++;
}


inline int MEDDLY::mt_forest::find(int node) {
  return unique->find(node);
}


inline int MEDDLY::mt_forest::insert(int node) {
  return unique->insert(node);
}


inline int MEDDLY::mt_forest::replace(int node) {
  return unique->replace(node);
}

inline int MEDDLY::mt_forest::getTempNodeId() const {
  return temp_node;
}

// ****************** override expert_forest class ****************** 


inline void MEDDLY::mt_forest::compactMemory() {
  compactAllLevels();
}

inline void MEDDLY::mt_forest::showInfo(FILE* strm, int verbosity) {
  showAll(strm, verbosity);
  fprintf(strm, "DD stats:\n");
  reportMemoryUsage(strm, '\t');
  fprintf(strm, "Unique table stats:\n");
  unique->showInfo(strm);
}

// *************** override expert_forest class -- done ***************

inline int MEDDLY::mt_forest::getLevelNode(int k) const {
  return levels[k].levelNode;
}

inline bool MEDDLY::mt_forest::isValidVariable(int vh) const {
  return (vh > 0) && (vh <= getExpertDomain()->getNumVariables());
  //return expertDomain->getVariableHeight(vh) != -1;
}

inline void
MEDDLY::mt_forest::findFirstElement(const dd_edge& f, int* vlist) const
{
  throw error(error::INVALID_OPERATION);
}

inline void
MEDDLY::mt_forest::findFirstElement(const dd_edge& f, int* vlist, int* vplist) const
{
  throw error(error::INVALID_OPERATION);
}


#endif

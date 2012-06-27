
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

// #define CIARDO_IDENTITY


namespace MEDDLY {
  class mt_forest;
};


/*
 * The MDD node handle is an integer.
 * The MDD node handle represents the offset in address[] where the node data
 * is stored.
 */

class MEDDLY::mt_forest : public expert_forest {
  public:
    mt_forest(int dsl, domain *d, bool rel, range_type t, edge_labeling ev, 
      const policies &p);
    virtual ~mt_forest();


  // ------------------------------------------------------------
  // virtual and overriding default behavior
  public: 

    virtual void createEdgeForVar(int vh, bool pr, bool* terms, dd_edge& a);
    virtual void createEdgeForVar(int vh, bool pr, int* terms, dd_edge& a);
    virtual void createEdgeForVar(int vh, bool pr, float* terms, dd_edge& a);

  // ------------------------------------------------------------
  // virtual and overriding default behavior
  protected:
    virtual int createReducedHelper(int i, const nodeBuilder& nb);

  



  // ------------------------------------------------------------
  // still to be organized:
  public:
    using expert_forest::getDownPtrsAndEdgeValues;
    using expert_forest::getSparseNodeIndexes;

    /// Refer to meddly.h
    void createSubMatrix(const bool* const* vlist,
        const bool* const* vplist, const dd_edge a, dd_edge& b);
    void createSubMatrix(const dd_edge& rows, const dd_edge& cols,
        const dd_edge& a, dd_edge& b);

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

    // virtual void dumpUniqueTable(FILE* s) const;
    void reclaimOrphanNode(int node);     // for linkNode()
    void handleNewOrphanNode(int node);   // for unlinkNode()
    void deleteOrphanNode(int node);      // for uncacheNode()
    // void freeZombieNode(int node);        // for uncacheNode()

    bool discardTemporaryNodesFromComputeCache() const;   // for isStale()

    void showNode(FILE *s, int p, int verbose = 0) const;
    void showNodeGraph(FILE *s, int p) const;

    // *************** override expert_forest class -- done ***************

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

    void showNode(int p) const;
    void showAll() const;

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

    // long getUniqueTableMemoryUsed() const;
    long getHoleMemoryUsage() const;

    // zombify node p
    void zombifyNode(int p);

    // delete node p
    void deleteNode(int p);

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

    // get the id used to indicate temporary nodes
    int getTempNodeId() const;

    // Checks if one of the following reductions is satisfied:
    // Fully, Quasi, Identity Reduced.
    // If the node can be reduced to an existing node, the existing node
    // is returned.
    // OBSOLETE
    bool checkForReductions(int p, int nnz, int& result);

#ifndef CIARDO_IDENTITY
    bool checkForReductions(const nodeBuilder &nb, int nnz, int& result);
#endif

    // Checks if the node has a single downpointer enabled and at
    // the given index.
    bool singleNonZeroAt(int p, int val, int index) const;

    // Returns zero, or the downward pointer if this is a
    // node with one downward pointer at the given index.
    int singleNonZeroAt(int p, int index) const;

    // Checks if the node satisfies the forests reduction rules.
    // If it does not, an assert violation occurs.
    void validateDownPointers(int p, bool recursive = false);

    // Special next values
    static const int temp_node = -5;

    // performance stats

    /// Count of nodes created since last gc
    unsigned nodes_activated_since_gc;

    // Deleting terminal nodes (used in isStale() -- this enables
    // the removal of compute cache entries which refer to terminal nodes)
    bool delete_terminal_nodes;

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
  return (isTerminalNode(p) || (getNext(p)>=0));
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


inline int MEDDLY::mt_forest::getTempNodeId() const {
  return temp_node;
}

inline int MEDDLY::mt_forest::getLevelNode(int k) const {
  return levels[k].levelNode;
}

inline bool MEDDLY::mt_forest::isValidVariable(int vh) const {
  return (vh > 0) && (vh <= getExpertDomain()->getNumVariables());
  //return expertDomain->getVariableHeight(vh) != -1;
}


#endif

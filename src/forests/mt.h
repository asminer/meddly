
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

    
    virtual char edgeSize(int k) const {
        return 0;
    }
    virtual char unhashedHeaderSize(int k) const {
        return 0;
    }
    virtual char hashedHeaderSize(int k) const {
        return 0;
    }
    virtual bool areEdgeValuesHashed(int k) const {
        return false;
    }
    virtual bool areDuplicates(int node, const node_builder &nb) const;
    virtual bool areDuplicates(int node, const node_reader &nr) const;

  private:
    template <class T>
    inline bool areDupsInternal(int p, const T &nb) const {
        const nodeData &node = getNode(p);
        if (node.level != nb.getLevel()) return false;
        const level_data &ld = levels[node.level];
        if (ld.isFull(node.offset)) {
          //
          // p is full
          //
          int fs = ld.fullSizeOf(node.offset);
          const int* pd = ld.fullDownOf(node.offset);
          if (nb.isFull()) {
            int i;
            for (i=0; i<fs; i++) if (pd[i] != nb.d(i)) return false;
            for (; i<nb.getSize(); i++) if (nb.d(i)) return false;
            return true;
          }
          MEDDLY_DCASSERT(nb.isSparse()); 
          int i = 0;
          for (int z=0; z<nb.getNNZs(); z++) {
            if (nb.i(z) >= fs) return false;
            for (; i<nb.i(z); i++) if (pd[i]) return false;
            if (pd[i] != nb.d(z)) return false;
            i++;
          } // for z
          for (; i<fs; i++) if (pd[i]) return false;
          return true;
        }
        //
        // p is sparse
        //
        int nnz = ld.sparseSizeOf(node.offset);
        const int* pd = ld.sparseDownOf(node.offset);
        const int* pi = ld.sparseIndexesOf(node.offset);

        if (nb.isSparse()) {
          if (nnz != nb.getNNZs()) return false;
          for (int z=0; z<nnz; z++) {
            if (nb.d(z) != pd[z]) return false;
            if (nb.i(z) != pi[z]) return false;
          }
          return true;
        }
        MEDDLY_DCASSERT(nb.isFull()); 
        int i = 0;
        for (int z=0; z<nnz; z++) {
          for (; i<pi[z]; i++) if (nb.d(i)) return false;
          if (nb.d(i) != pd[z]) return false;
          i++;
        }
        for (; i<nb.getSize(); i++) if (nb.d(i)) return false;
        return true;
    }

  // ------------------------------------------------------------
  // virtual and overriding default behavior
  protected:
    // virtual int createReducedHelper(int in, const node_builder& nb, bool &u);

  
  // ------------------------------------------------------------
  // Helpers for this and derived classes
  protected:
    /// Add a redundant node at level k.
    /// On input: d is the node "below" us, to point to.
    /// On output: d is the redundant node.
    inline void insertRedundantNode(int k, int& d) {
      MEDDLY_DCASSERT(!isFullyReduced());
      bool useIdentity = false;
      if (isIdentityReduced()) {
        useIdentity = getNodeLevel(d) < 0; 
      }
      int sz = getLevelSize(k);
      node_builder& nb = useNodeBuilder(k, sz);
      if (useIdentity) {
        // check for identity reductions with d
        int sd;
        int si = getSingletonIndex(d, sd);
        for (int i=0; i<sz; i++) {
          nb.d(i) = linkNode(
            (i==si) ? sd : d
          );
        }
      } else {
        // easy - just set all the pointers
        for (int i=0; i<sz; i++) nb.d(i) = linkNode(d);
      }
      unlinkNode(d);
      d = createReducedNode(-1, nb);
    }

    // Sanity check for development code.
    // void validateDownPointers(const node_builder &nb) const;


  // ------------------------------------------------------------
  // still to be organized:
  public:
    // using expert_forest::getDownPtrsAndEdgeValues;
    // using expert_forest::getSparseNodeIndexes;

    /// Refer to meddly.h
    /*
    void createSubMatrix(const bool* const* vlist,
        const bool* const* vplist, const dd_edge a, dd_edge& b);
    void createSubMatrix(const dd_edge& rows, const dd_edge& cols,
        const dd_edge& a, dd_edge& b);
*/

#ifdef ACCUMULATE_ON
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

#endif

    /// Has the node been reduced
    bool isReducedNode(int node) const;

    // virtual void dumpUniqueTable(FILE* s) const;
    // void reclaimOrphanNode(int node);     // for linkNode()
    // void handleNewOrphanNode(int node);   // for unlinkNode()
    // void deleteOrphanNode(int node);      // for uncacheNode()
    // void freeZombieNode(int node);        // for uncacheNode()

    // virtual void showNode(FILE *s, int p, int verbose = 0) const;

    // *************** override expert_forest class -- done ***************
#if 0
    // Dealing with slot 2 (node size)
    int getLargestIndex(int p) const;

    // Dealing with entries

    // full node: entries start in the 4th slot (location 3, counting from 0)
    int* getFullNodeDownPtrs(int p);
    const int* getFullNodeDownPtrsReadOnly(int p) const;
    const int* getSparseNodeIndexes(int p) const;
    const int* getSparseNodeDownPtrs(int p) const;
    int getSparseNodeLargestIndex(int p) const;

    // for EVMDDs
    int* getFullNodeEdgeValues(int p);
    const int* getFullNodeEdgeValuesReadOnly(int p) const;
    const int* getSparseNodeEdgeValues(int p) const;


    // p: node
    // i: the ith downpointer.
    // note: for sparse nodes this may not be the same as the ith index pointer.
    // int getDownPtrAfterIndex(int p, int i, int &index) const;

    int getMddLevelMaxBound(int k) const;
    int getMxdLevelMaxBound(int k) const;
    int getLevelMaxBound(int k) const;
#endif
    int buildQuasiReducedNodeAtLevel(int k, int p);

    void compareCacheCounts(int p = -1);
    void validateIncounts();


    // Remove zombies if more than max
    // void removeZombies(int max = 100);

    // bool isCounting();

  protected:
    // Building level nodes
    int getLevelNode(int lh) const;
    int buildLevelNodeHelper(int lh, int* terminalNodes, int sz);
    void buildLevelNode(int lh, int* terminalNodes, int sz);
    void clearLevelNode(int lh);
    void clearLevelNodes();
    // void clearAllNodes();

    // Building custom level nodes
    // int* getTerminalNodes(int n);
    int* getTerminalNodes(int n, bool* terms);
    int* getTerminalNodes(int n, int* terms);
    int* getTerminalNodes(int n, float* terms);

    bool isValidVariable(int vh) const;
    bool doesLevelNeedCompaction(int k);

    // Dealing with node addressing
    // void setNodeOffset(int p, int offset);

    // Dealing with node status
    // bool isDeletedNode(int p) const;

    // long getUniqueTableMemoryUsed() const;
    long getHoleMemoryUsage() const;



  protected:

    // modify temp nodes count for level k as well as the global count
    // the temporary node count should be incremented only within
    // createTempNode() or variants.
    // decrTempNodeCount() should be called by any method that changes a
    // temporary node to a reduced node.
    // Note: deleting a temp node automatically calls decrTempNodeCount().

    // increment the count for "nodes activated since last garbage collection"
    void incrNodesActivatedSinceGc();

    // get the id used to indicate temporary nodes
    int getTempNodeId() const;

    // Checks if one of the following reductions is satisfied:
    // Fully, Quasi, Identity Reduced.
    // If the node can be reduced to an existing node, the existing node
    // is returned.
    // OBSOLETE
    // bool checkForReductions(int p, int nnz, int& result);

    // Checks if the node has a single downpointer enabled and at
    // the given index.
    // bool singleNonZeroAt(int p, int val, int index) const;

    // Checks if the node satisfies the forests reduction rules.
    // If it does not, an assert violation occurs.
    // void validateDownPointers(int p, bool recursive = false);

    // Special next values
    static const int temp_node = -5;

    // performance stats

    /// Count of nodes created since last gc
    unsigned nodes_activated_since_gc;

    bool counting;

    // scratch pad for buildLevelNode and getTerminalNodes
    int* dptrs;
    int dptrsSize;

    // Place holder for accumulate-minterm result.
    bool accumulateMintermAddedElement;

  private:
#ifdef ACCUMULATE_ON
    // Persistant variables used in addReducedNodes()
    dd_edge* nodeA;
    dd_edge* nodeB;
#endif
};


/// Inline functions implemented here

/*
inline void MEDDLY::mt_forest::reclaimOrphanNode(int p) {
  MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));
  MEDDLY_DCASSERT(isActiveNode(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(isReducedNode(p));
  stats.reclaimed_nodes++;
  stats.orphan_nodes--;
}  
*/

/*
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
*/


/*
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
*/

// Dealing with slot 0 (incount)

// Dealing with slot 1 (next pointer)

inline bool MEDDLY::mt_forest::isReducedNode(int p) const {
#ifdef DEBUG_MDD_H
  printf("%s: p: %d\n", __func__, p);
#endif
  MEDDLY_DCASSERT(isActiveNode(p));
  return (isTerminalNode(p) || (getNext(p)>=0));
}

/*

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
*/

// inline bool MEDDLY::mt_forest::isCounting() { return counting; }

// Dealing with node addressing

/*
inline void MEDDLY::mt_forest::setNodeOffset(int p, int offset) 
{
  MEDDLY_DCASSERT(isValidNonterminalIndex(p));
  address[p].offset = offset;
}
*/

// Dealing with node status

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

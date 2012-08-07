
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

// #include <fstream>
// #include <iostream>
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

    
    virtual char edgeSize() const {
        return 0;
    }
    virtual char unhashedHeaderSize() const {
        return 0;
    }
    virtual char hashedHeaderSize() const {
        return 0;
    }
    virtual bool areEdgeValuesHashed() const {
        return false;
    }
    virtual bool areDuplicates(long node, const node_builder &nb) const;
    virtual bool areDuplicates(long node, const node_reader &nr) const;

    virtual bool isRedundant(const node_builder &nb) const;
    virtual bool isIdentityEdge(const node_builder &nb, int i) const;

  private:
    template <class T>
    inline void edgeForVarInternal(int vh, bool pr, T* terms, dd_edge& result) {
        if (vh < 0 || vh > getNumVariables())
            throw error(error::INVALID_VARIABLE);
        if (result.getForest() != this) 
            throw error(error::INVALID_OPERATION);
        int k = pr ? -vh: vh;
        MEDDLY_DCASSERT(isValidLevel(k));

        if (!isForRelations() && pr) 
            throw error(error::INVALID_ASSIGNMENT);
        if (getEdgeLabeling() != MULTI_TERMINAL)
            throw error(error::INVALID_OPERATION);
        long *terminalNodes = getTerminalNodes(getLevelSize(vh), terms);
        long node = buildLevelNodeHelper(k, terminalNodes, getLevelSize(vh));

        result.set(node, 0);
    }

    template <class T>
    inline bool areDupsInternal(long p, const T &nb) const {
        const node_header &node = getNode(p);
        if (node.level != nb.getLevel()) return false;
#ifdef NODE_STORAGE_PER_LEVEL
        const node_storage &nodeMan = levels[node.level];
#endif
        node_storage::scanner np;
        nodeMan.initScanner(np, node.offset);
        if (np.isFull()) {
          //
          // p is full
          //
          if (nb.isFull()) {
            for (; np; ++np) {
              if (np.d() != nb.d(np.count())) return false;
            }
            for (int i=np.count(); i<nb.getSize(); i++) 
              if (nb.d(i)) return false;
            return true;
          }
          MEDDLY_DCASSERT(nb.isSparse()); 
          for (int z=0; z<nb.getNNZs(); z++) {
            if (nb.i(z) >= np.fullSize()) return false;
            for (; np.count()<nb.i(z); ++np) if (np.d()) return false;
            if (np.d() != nb.d(z)) return false;
            ++np;
          } // for z
          for (; np; ++np) if (np.d()) return false;
          return true;
        } 
        //
        // p is sparse
        //
        if (nb.isSparse()) {
          if (np.sparseSize() != nb.getNNZs()) return false;
          for (int z=0; z<np.sparseSize(); ++z, ++np) {
            if (nb.d(z) != np.d()) return false;
            if (nb.i(z) != np.i()) return false;
          }
          return true;
        }
        MEDDLY_DCASSERT(nb.isFull()); 
        int i = 0;
        for (; np; ++np) {
          for (; i<np.i(); i++) if (nb.d(i)) return false;
          if (nb.d(i) != np.d()) return false;
          i++;
        }
        for (; i<nb.getSize(); i++) if (nb.d(i)) return false;
        return true;
    }

  // ------------------------------------------------------------
  // Helpers for this and derived classes
  protected:
    /// Add a redundant node at level k.
    /// On input: d is the node "below" us, to point to.
    /// On output: d is the redundant node.
    inline void insertRedundantNode(int k, long& d) {
      MEDDLY_DCASSERT(!isFullyReduced());
      bool useIdentity = false;
      if (isIdentityReduced()) {
        useIdentity = getNodeLevel(d) < 0; 
      }
      int sz = getLevelSize(k);
      node_builder& nb = useNodeBuilder(k, sz);
      if (useIdentity) {
        // check for identity reductions with d
        long sd;
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

  // ------------------------------------------------------------
  // still to be organized:
  public:
    // void compareCacheCounts(int p = -1);


  protected:
    // Building level nodes
    long buildLevelNodeHelper(int lh, long* terminalNodes, int sz);

    // Building custom level nodes
    long* getTerminalNodes(int n, bool* terms);
    long* getTerminalNodes(int n, int* terms);
    long* getTerminalNodes(int n, float* terms);

  protected:

    bool counting;

    // scratch pad for buildLevelNode and getTerminalNodes
    long* dptrs;
    int dptrsSize;

};


#endif

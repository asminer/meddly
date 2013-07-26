
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
  protected:
    mt_forest(int dsl, domain *d, bool rel, range_type t, edge_labeling ev, 
      const policies &p);
    virtual ~mt_forest();


  // ------------------------------------------------------------
  // virtual and overriding default behavior
  public: 

    virtual void createEdgeForVar(int vh, bool pr, bool* terms, dd_edge& a);
    virtual void createEdgeForVar(int vh, bool pr, int* terms, dd_edge& a);
    virtual void createEdgeForVar(int vh, bool pr, float* terms, dd_edge& a);

    
    virtual void writeNode(FILE* s, const node_reader& nr, 
      const node_handle* map) const;
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
        node_handle *terminalNodes = getTerminalNodes(getLevelSize(vh), terms);
        long node = buildLevelNodeHelper(k, terminalNodes, getLevelSize(vh));

        result.set(node, 0);
    }

  // ------------------------------------------------------------
  // Helpers for this and derived classes
  protected:
    /// Add a redundant node at level k.
    /// On input: d is the node "below" us, to point to.
    /// On output: d is the redundant node.
    inline void insertRedundantNode(int k, node_handle& d) {
      MEDDLY_DCASSERT(!isFullyReduced());
      bool useIdentity = false;
      if (isIdentityReduced()) {
        useIdentity = getNodeLevel(d) < 0; 
      }
      int sz = getLevelSize(k);
      node_builder& nb = useNodeBuilder(k, sz);
      if (useIdentity) {
        // check for identity reductions with d
        node_handle sd;
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

  protected:
    // Building level nodes
    node_handle buildLevelNodeHelper(int lh, node_handle* terminalNodes, int sz);

    // Building custom level nodes
    node_handle* getTerminalNodes(int n, bool* terms);
    node_handle* getTerminalNodes(int n, int* terms);
    node_handle* getTerminalNodes(int n, float* terms);

  private:  

    bool counting;

    // scratch pad for buildLevelNode and getTerminalNodes
    node_handle* dptrs;
    int dptrsSize;

};


#endif

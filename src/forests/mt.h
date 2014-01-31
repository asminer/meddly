
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

#include "../defines.h"

namespace MEDDLY {
  class mt_base_forest;
};


/**
    Base class for all multi-terminal DDs.
    Common things, that do not depend on the terminal type,
    are implemented here.
*/
class MEDDLY::mt_base_forest : public expert_forest {
  protected:
    mt_base_forest(int dsl, domain *d, bool rel, range_type t, const policies &p);

  public:
    virtual bool isRedundant(const node_builder &nb) const;
    virtual bool isIdentityEdge(const node_builder &nb, int i) const;
  
  // ------------------------------------------------------------
  // Helpers for this and derived classes
  protected:
    /**
        Add redundant nodes from level k to the given node.
    */
    inline node_handle makeNodeAtLevel(int k, node_handle d) {
      MEDDLY_DCASSERT(abs(k) >= abs(getNodeLevel(d)));
      if (0==d) return d;
      if (isFullyReduced()) return d;
      int dk = getNodeLevel(d); 
      while (dk != k) {
        int up;
        if (dk<0) up = -dk;
        else up = isForRelations() ? -(dk+1) : dk+1;

        // make node at level "up"
        int sz = getLevelSize(up);
        node_builder& nb = useNodeBuilder(up, sz);

        if (isIdentityReduced() && (dk<0)) {
          // make identity reductions below as necessary
          node_handle sd;
          int si = isTerminalNode(d) ? -1 : getSingletonIndex(d, sd);
          for (int i=0; i<sz; i++) {
            nb.d(i) = linkNode( (i==si) ? sd : d );
          }
        } else {
          // don't worry about identity reductions
          for (int i=0; i<sz; i++) {
            nb.d(i) = linkNode(d);
          }
        }
        unlinkNode(d);
        d = createReducedNode(-1, nb);
        dk = up;
      } // while
      return d;
    }

    /// make a node at the top level
    inline node_handle makeNodeAtTop(node_handle d) {
      return makeNodeAtLevel(getDomain()->getNumVariables(), d);
    }

    /**
        Special case for createEdge(), with only one minterm.
        For sets; used in mtmdd.h 
    */
    inline node_handle createEdgePath(int k, int* vlist, node_handle bottom) {
        MEDDLY_DCASSERT(!isForRelations());
        if (0==bottom) return bottom;
        for (int i=1; i<=k; i++) {
          if (DONT_CARE == vlist[i]) {
            // make a redundant node
            if (isFullyReduced()) continue; 
            int sz = getLevelSize(i);
            node_builder& nb = useNodeBuilder(i, sz);
            nb.d(0) = bottom;
            for (int v=1; v<sz; v++) {
              nb.d(v) = linkNode(bottom);
            }
            bottom = createReducedNode(-1, nb);
          } else {
            // make a singleton node
            node_builder& nb = useSparseBuilder(i, 1);
            nb.i(0) = vlist[i];
            nb.d(0) = bottom;
            bottom = createReducedNode(-1, nb);
          }
        } // for i
        return bottom;
    }

    /**
        Special case for createEdge(), with only one minterm.
        For relations; used in mtmxd.h
    */
    inline 
    node_handle createEdgePath(int k, int* vlist, int* vplist, node_handle next) 
    {
        MEDDLY_DCASSERT(isForRelations());
        if (0==next) return next;
        for (int i=1; i<=k; i++) {
          if (DONT_CHANGE == vplist[i]) {
            //
            // Identity node
            //
            MEDDLY_DCASSERT(DONT_CARE == vlist[i]);
            if (isIdentityReduced()) continue;
            // Build an identity node by hand
            int sz = getLevelSize(i);
            node_builder& nb = useNodeBuilder(i, sz);
            for (int v=0; v<sz; v++) {
              node_builder& nbp = useSparseBuilder(-i, 1);
              nbp.i(0) = v;
              nbp.d(0) = linkNode(next);
              nb.d(v) = createReducedNode(v, nbp);
            }
            unlinkNode(next);
            next = createReducedNode(-1, nb);
            continue;
          }
          //
          // process primed level
          //
          node_handle nextpr;
          if (DONT_CARE == vplist[i]) {
            if (isFullyReduced()) {
              // DO NOTHING
              nextpr = next;
            } else {
              // build redundant node
              int sz = getLevelSize(-i);
              node_builder& nb = useNodeBuilder(-i, sz);
              for (int v=0; v<sz; v++) {
                nb.d(v) = linkNode(next);
              }
              unlinkNode(next);
              nextpr = createReducedNode(-1, nb);
            }
          } else {
            // sane value
            node_builder& nb = useSparseBuilder(-i, 1);
            nb.i(0) = vplist[i];
            nb.d(0) = next;
            nextpr = createReducedNode(vlist[i], nb);
          }
          //
          // process unprimed level
          //
          if (DONT_CARE == vlist[i]) {
            if (isFullyReduced()) continue;
            // build redundant node
            int sz = getLevelSize(i);
            node_builder& nb = useNodeBuilder(i, sz);
            if (isIdentityReduced()) {
              // Below is likely a singleton, so check for identity reduction
              // on the appropriate v value
              for (int v=0; v<sz; v++) {
                node_handle dpr = (v == vplist[i]) ? next : nextpr;
                nb.d(v) = linkNode(dpr);
              }
            } else {
              // Doesn't matter what happened below
              for (int v=0; v<sz; v++) {
                nb.d(v) = linkNode(nextpr);
              }
            }
            unlinkNode(nextpr);
            next = createReducedNode(-1, nb);
          } else {
            // sane value
            node_builder& nb = useSparseBuilder(i, 1);
            nb.i(0) = vlist[i];
            nb.d(0) = nextpr;
            next = createReducedNode(-1, nb);
          }
        } // for i
        return next;
    }
    
    /**
        Enlarge variables to include all given minterms.
    */
    inline node_handle enlargeVariables(int** const vlist, int N, bool primed) {
      for (int k=1; k<=getDomain()->getNumVariables(); k++) {
        int maxv = vlist[0][k];
        for (int i=1; i<N; i++) {
          maxv = MAX(maxv, vlist[i][k]);
        }
        if (maxv < 1) continue;
        if (maxv >= getDomain()->getVariableBound(k, primed)) {
          useExpertDomain()->enlargeVariableBound(k, primed, maxv+1);
        }
      }
    }


  protected:

    // set, and used, in derived classes for createEdge()
    binary_operation* unionOp;
};


namespace MEDDLY {

/**
    Base class for all multi-terminal DDs.
    It's a template class, based on the terminal type,
    so we can implement certain things here, quickly.
*/
template <class ENCODER>
class mt_forest : public mt_base_forest {
  protected:
    mt_forest(int dsl, domain *d, bool rel, range_type t, const policies &p);

    template <class T>
    void createEdgeForVarTempl(int vh, bool pr, const T* terms, dd_edge& result);

    template <class T>
    inline void createEdgeTempl(T term, dd_edge& e) {
      if (e.getForest() != this) throw error(error::INVALID_OPERATION);
      e.set(makeNodeAtTop(ENCODER::value2handle(term)), 0);
    }

  public:
    virtual void showTerminal(FILE* s, node_handle tnode) const;
    virtual void writeTerminal(FILE* s, node_handle tnode) const;
    virtual node_handle readTerminal(FILE* s);
};

//
// Template implementation
//
template <class ENCODER>
mt_forest<ENCODER>
::mt_forest(int dsl, domain *d, bool rel, range_type t, const policies &p)
 : mt_base_forest(dsl, d, rel, t, p)
{
  // nothing to construct
}


template <class ENCODER>
template <class T>
void mt_forest<ENCODER>::
createEdgeForVarTempl(int vh, bool pr, const T* termvals, dd_edge& result)
{
    /*
        Sanity checks
    */
    if (vh < 0 || vh > getNumVariables())
        throw error(error::INVALID_VARIABLE);
    if (result.getForest() != this) 
        throw error(error::INVALID_OPERATION);
    if (!isForRelations() && pr) 
        throw error(error::INVALID_ASSIGNMENT);

    /*
        Get info for node we're building
    */
    int k = pr ? -vh: vh;
    int km1;
    if (isForRelations()) {
      km1 = (k<0) ? (-k)-1 : -k;
    } else {
      km1 = k-1;
    }
    int sz = getLevelSize(vh);

    /*
        Make this node
    */
    node_builder &nb = useNodeBuilder(k, sz);
    for (int i=0; i<sz; i++) {
      nb.d(i) = makeNodeAtLevel(km1, 
        ENCODER::value2handle(termvals ? termvals[i] : i)
      );
    }

    /*
        Reduce, add redundant as necessary, and set answer
    */
    node_handle node = createReducedNode(-1, nb);
    node = makeNodeAtLevel(getDomain()->getNumVariables() , node); 
    result.set(node, 0);
}



template <class ENCODER>
void mt_forest<ENCODER>::showTerminal(FILE* s, node_handle tnode) const
{
  ENCODER::show(s, tnode);
}

template <class ENCODER>
void mt_forest<ENCODER>::writeTerminal(FILE* s, node_handle tnode) const
{
  ENCODER::write(s, tnode);
}

template <class ENCODER>
node_handle mt_forest<ENCODER>::readTerminal(FILE* s)
{
  return ENCODER::read(s);
}


};  // end MEDDLY namespace

#endif

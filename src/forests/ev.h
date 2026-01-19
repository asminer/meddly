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

#ifndef MEDDLY_EV_FOREST
#define MEDDLY_EV_FOREST

#include "../defines.h"
#include "../minterms.h"
#include "../dd_edge.h"
#include "../forest.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"

namespace MEDDLY {
  class ev_forest;
};

/**
    Base class for all edge-valued DDs.

    Derived classes probably want to build an OPERATION class
    that contains inlined static methods:
        void setEdge(void*, type v);
*/
class MEDDLY::ev_forest : public forest
{
  protected:
    ev_forest(domain *d, bool rel, range_type t, edge_labeling ev,
      const policies &p);

  public:

  // ------------------------------------------------------------
  // Helpers for this and derived classes

#ifdef ALLOW_DEPRECATED_0_18_0
    /** Add redundant nodes from level k to the given node.
          @param  k   Top level we want
          @param  ev  In/Out: Edge value
          @param  ed  In/Out: Edge node pointer
    */
    template <class OPERATION, typename TYPE>
    void makeNodeAtLevel(int k, TYPE &ev, node_handle &ed);
#endif

  protected:
#ifdef ALLOW_DEPRECATED_0_18_0
    /// make a node at the top level
    template <class OPERATION, typename TYPE>
    inline void makeNodeAtTop(TYPE &ev, node_handle &ed) {
      makeNodeAtLevel<OPERATION>((int) getDomain()->getNumVariables(), ev, ed);
    }
#endif

#ifdef ALLOW_DEPRECATED_0_17_7
    /**
        Enlarge variables to include all given minterms.
    */
    inline void enlargeVariables(const int* const* vlist, int N, bool primed) {
      for (unsigned k=1; k<=getDomain()->getNumVariables(); k++) {
        int maxv = vlist[0][k];
        for (int i=1; i<N; i++) {
          maxv = MAX(maxv, vlist[i][k]);
        }
        if (maxv < 1) continue;
        if (maxv >= getDomain()->getVariableBound(k, primed)) {
          variable* vh = getDomain()->getVar(k);
          vh->enlargeBound(primed, (maxv+1));
        }
      }
    }
#endif

#ifdef ALLOW_DEPRECATED_0_18_0

    template <class OPERATION, typename T>
    inline
    void createEdgeForVarTempl(int vh, bool pr, const T* vals, dd_edge& result)
    {
      /*
          Sanity checks
      */
      if (vh < 0 || vh > getNumVariables())
        throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);
      if (!result.isAttachedTo(this))
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
      if (!isForRelations() && pr)
        throw error(error::INVALID_ASSIGNMENT, __FILE__, __LINE__);

      int level = getLevelByVar(vh);

      /*
          Get info for node we're building
      */
      int k = pr ? -level : level;
      int km1;
      if (isForRelations()) {
        km1 = (k<0) ? (-k)-1 : -k;
      } else {
        km1 = k-1;
      }

      /*
          Make this node
      */
      unpacked_node *nb = unpacked_node::newWritable(this, k, FULL_ONLY);
      for (unsigned i=0; i<nb->getSize(); i++) {
        T ev = vals ? vals[i] : i;
        terminal t;
        t.setOmega();
        node_handle ed = t.getHandle();
        makeNodeAtLevel<OPERATION, T>(km1, ev, ed);
        nb->setFull(i, ev, ed);
      }

      /*
          Reduce, add redundant as necessary, and set answer
      */
      T ev;
      node_handle node;
      createReducedNode(-1, nb, ev, node);
      makeNodeAtTop<OPERATION, T>(ev, node);
      result.set(ev, node);
    }


    template <class OPERATION, class T>
    inline void createEdgeTempl(T term, dd_edge& e) {
        if (!e.isAttachedTo(this)) {
          throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
        }
        terminal t;
        t.setOmega();
        node_handle ed = t.getHandle();
        makeNodeAtTop<OPERATION, T>(term, ed);
        e.set(term, ed);
    }
#endif


  // statics

#ifdef ALLOW_DEPRECATED_0_17_7
  public:
    static void initStatics();
    static void enlargeStatics(int n);
    static void clearStatics();

  protected:
    static int* order;
    static int order_size;

#endif
};






#ifdef ALLOW_DEPRECATED_0_18_0

namespace MEDDLY {

  template <class OPERATION, typename TYPE>
  void ev_forest::makeNodeAtLevel(int k, TYPE &ev, node_handle &ed)
  {
    MEDDLY_DCASSERT(abs(k) >= abs(getNodeLevel(ed)));
    if (0==ed) return;
    if (isFullyReduced()) return;

    // prevSize == 1 enables getSingletonIndex for the first run.
    unsigned prevSize = 1;
    int dk = getNodeLevel(ed);
    while (dk != k) {
      // make node at one level up
      int up = (dk < 0) ? -dk : isForRelations() ? -(dk + 1) : (dk + 1);


      // unsigned sz = unsigned(getLevelSize(up));
      // unpacked_node* nb = unpacked_node::newFull(this, up, sz);
      unpacked_node* nb = unpacked_node::newWritable(this, up, FULL_ONLY);
      unsigned sz = nb->getSize();

      if (isIdentityReduced() && (dk < 0) && (1 == prevSize)) {
        // Build unprimed node with check for identity reduction.
        // Note 0: ed cannot be a terminal node (dk < 0).
        // Note 1: if prevSize > 1, and ed is not the original argument,
        //         and if it is a primed node, it cannot be identity reduced.
        MEDDLY_DCASSERT(!isTerminalNode(ed));
        node_handle sd;
        unsigned si;
        if (!isSingletonNode(ed, si, sd))  {
            si = sz;
        }
        for (unsigned i = 0; i < sz; i++) {
          nb->setFull(i, ev, linkNode( (si == i) ? sd : ed ) );
        }
      } else {
        // No identity reduction possible.
        for (unsigned i = 0; i < sz; i++) {
          nb->setFull(i, ev, linkNode(ed));
        }
      }

      unlinkNode(ed);
      createReducedNode(-1, nb, ev, ed);
      dk = up;
      prevSize = sz;
    } // while
  }

};  // namespace MEDDLY

#endif

#endif


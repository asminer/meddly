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

#ifndef MEDDLY_MT_FOREST
#define MEDDLY_MT_FOREST

#include "../defines.h"
#include "../minterms.h"
#include "../dd_edge.h"
#include "../forest.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"

namespace MEDDLY {
  class mt_forest;
};

/**
    Base class for all multi-terminal DDs.
    Common things, that do not depend on the terminal type,
    are implemented here.
*/
class MEDDLY::mt_forest : public expert_forest {
  protected:
    mt_forest(domain *d, bool rel, range_type t, const policies &p, int* level_reduction_rule=NULL);

  public:
    virtual bool isRedundant(const unpacked_node &nb) const;
    virtual bool isIdentityEdge(const unpacked_node &nb, int i) const;

  // ------------------------------------------------------------
  // Helpers for this and derived classes

    /// Add redundant nodes from level k to the given node.
    node_handle makeNodeAtLevel(int k, node_handle d);

        // No edge values right now, so these are no-ops
        // That could change with complement edges.
        /*
        virtual void readEdgeValue(input &s, dd_edge &E) const;
        virtual void writeEdgeValue(output &s, const dd_edge &E) const;
        virtual void showEdgeValue(output &s, const dd_edge &E) const;
        */

  protected:
    /// make a node at the top level
    inline node_handle makeNodeAtTop(node_handle p) {
      return makeNodeAtLevel(getDomain()->getNumVariables(), p);
    }

    /**
        Enlarge variables to include all given minterms.
    */
    inline void enlargeVariables(const int* const* vlist, int N, bool primed) {
      for (int k=1; k<=getDomain()->getNumVariables(); k++) {
        int maxv = vlist[0][k];
        for (int i=1; i<N; i++) {
          maxv = MAX(maxv, vlist[i][k]);
        }
        if (maxv < 1) continue;
        if (maxv >= getDomain()->getVariableBound(k, primed)) {
          variable* vh = getDomain()->getVar(k);
          if (vh->isExtensible())
            vh->enlargeBound(primed, -(maxv+1));
          else
            vh->enlargeBound(primed, (maxv+1));
        }
      }
    }

    // template <class ENCODER, typename T>
    template <typename T>
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
      unsigned sz = unsigned(getLevelSize(level));

      /*
          Make this node
      */
      unpacked_node* nb = unpacked_node::newFull(this, k, sz);
      for (unsigned i=0; i<sz; i++) {
        terminal t(vals ? vals[i] : T(i));
        nb->d_ref(i) = makeNodeAtLevel(km1, t.getHandle());
      }

      /*
          Reduce, add redundant as necessary, and set answer
      */
      node_handle node = createReducedNode(-1, nb);
      node = makeNodeAtTop(node);
      result.set(node);
    }


    // template <class ENCODER, class T>
    template <typename T>
    inline void createEdgeTempl(T term, dd_edge& e) {
        if (!e.isAttachedTo(this)) {
            throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
        }
        // e.set(makeNodeAtTop(ENCODER::value2handle(term)));
        terminal t(term);
        e.set(makeNodeAtTop(t.getHandle()));
    }

  protected:
    // iterator base class for MT
    class mt_iterator : public enumerator::iterator {
      public:
        mt_iterator(const expert_forest* F);
        virtual ~mt_iterator();

        virtual void getValue(int &termVal) const;
        virtual void getValue(float &termVal) const;
    };


  // statics

  public:
    static void initStatics();
    static void enlargeStatics(int n);
    static void clearStatics();

  protected:
    static int* order;
    static int order_size;
};


#endif

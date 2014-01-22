
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

#define NEW_MT

namespace MEDDLY {
  class int_terminal;
  class real_terminal;
};


/** Wrapper around integers as terminal values.
*/
class MEDDLY::int_terminal {
    int value;
  public:
    inline void setFromValue(int v) {
      // Sanity check
      MEDDLY_DCASSERT(4 == sizeof(node_handle));
      value = v; 
      if (v < -1073741824 || v > 1073741823) {
        // Can't fit in 31 bits (signed)
        throw error(error::OVERFLOW);
      }
    }
    inline void setFromHandle(node_handle h) {
      // << 1 kills the sign bit
      // >> 1 puts us back, and extends the (new) sign bit
      value = (h << 1) >> 1;
    }
    inline void unionValue(bool v) {
      if (v) value = v;
    }
    inline void unionValue(int v) {
      setFromValue(value + v);
    }
    inline node_handle toHandle() const {
      // set the sign bit, unless we're 0
      return (0==value) ? 0 : value | 0x80000000;
    }
    inline operator bool() const {
      return value;
    }
    inline operator int() const {
      return value;
    }
    inline void show(FILE* s, forest::range_type r) const {
      switch (r) {
        case forest::BOOLEAN:
          th_fputc(value ? 'T' : 'F', s);
          return;

        case forest::INTEGER:
          th_fprintf(s, "t%d", value);
          return;

        default:
          throw error(error::TYPE_MISMATCH);
      }
    }
    inline void write(FILE* s, forest::range_type r) const {
      show(s, r);
    }
    inline void read(FILE* s, forest::range_type r) {
      stripWS(s);
      char c = fgetc(s);
      switch (r) {
        case forest::BOOLEAN: 
            if ('T' == c || 'F' == c) {
              value = ('T' == c);
              return;
            }
            throw error(error::INVALID_FILE);

        case forest::INTEGER: 
            if ('t' == c) {
              th_fscanf(1, s, "%d", &value);
              return;
            }
            throw error(error::INVALID_FILE);
          
        default:
          throw error(error::TYPE_MISMATCH);
      }
    }
};


/** Wrapper around reals as terminal values.
*/
class MEDDLY::real_terminal {
    union {
      float value;
      int ivalue;
    };
  public:
    inline void setFromValue(float v) {
      // Sanity checks
      MEDDLY_DCASSERT(4 == sizeof(node_handle));
      MEDDLY_DCASSERT(sizeof(float) <= sizeof(node_handle));
      value = v;
    }
    inline void setFromHandle(node_handle h) {
      // remove sign bit
      ivalue = (h << 1);
    }
    inline void unionValue(float v) {
      value += v;
    }
    inline node_handle toHandle() const {
      if (0.0 == value) return 0;
      // Strip lsb in fraction, and add sign bit
      return (ivalue >> 1) | 0x80000000;
    }
    inline operator float() const {
      return value;
    }
    inline void show(FILE* s, forest::range_type r) const {
      if (r != forest::REAL)
          throw error(error::TYPE_MISMATCH);
      th_fprintf(s, "t%f", value);
    }
    inline void write(FILE* s, forest::range_type r) const {
      if (r != forest::REAL)
          throw error(error::TYPE_MISMATCH);
      th_fprintf(s, "t%8e", value);
    }
    inline void read(FILE* s, forest::range_type r) {
      if (r != forest::REAL)
          throw error(error::TYPE_MISMATCH);
      stripWS(s);
      char c = fgetc(s);
      if ('t' == c) {
        th_fscanf(1, s, "%8e", &value);
      }
      throw error(error::INVALID_FILE);
    }
};


#ifdef NEW_MT

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
    /// Add redundant nodes from level k to the given node.
    inline node_handle makeNodeAtLevel(int k, node_handle d) {
      MEDDLY_DCASSERT(abs(k) >= abs(getNodeLevel(d)));
      if (isFullyReduced()) return d;
      int dk;
      while ((dk = getNodeLevel(d)) != k) {
        int up;
        if (dk<0) up = -dk;
        else up = isForRelations() ? -(dk+1) : dk+1;

        // make node at level "up"
        int sz = getLevelSize(up);
        node_builder& nb = useNodeBuilder(up, sz);

        if (isIdentityReduced() && (dk<0)) {
          // make identity reductions below as necessary
          node_handle sd;
          int si = getSingletonIndex(d, sd);
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
      } // while
      return d;
    }

    /// make a node at the top level
    inline node_handle makeNodeAtTop(node_handle d) {
      return makeNodeAtLevel(getDomain()->getNumVariables(), d);
    }


};


namespace MEDDLY {

/**
    Base class for all multi-terminal DDs.
    It's a template class, based on the terminal type,
    so we can implement certain things here, quickly.
*/
template <class TTYPE>
class mt_forest : public mt_base_forest {
  protected:
    mt_forest(int dsl, domain *d, bool rel, range_type t, const policies &p);

    template <class T>
    void createEdgeForVarTempl(int vh, bool pr, const T* terms, dd_edge& result);

    template <class T>
    inline void createEdgeTempl(T term, dd_edge& e) {
      if (e.getForest() != this) throw error(error::INVALID_OPERATION);
      TTYPE tnode;
      tnode.setFromValue(term);
      e.set(makeNodeAtTop(tnode.toHandle()), 0);
    }

  public:
    virtual void showTerminal(FILE* s, node_handle tnode) const;
    virtual void writeTerminal(FILE* s, node_handle tnode) const;
    virtual node_handle readTerminal(FILE* s);

};

//
// Template implementation
//
template <class TTYPE>
mt_forest<TTYPE>
::mt_forest(int dsl, domain *d, bool rel, range_type t, const policies &p)
 : mt_base_forest(dsl, d, rel, t, p)
{
  // nothing to construct
}


template <class TTYPE>
template <class T>
void mt_forest<TTYPE>::
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
    int sz = getLevelSize(vh);

    /*
        Make this node
    */
    node_builder &nb = useNodeBuilder(k, sz);
    for (int i=0; i<sz; i++) {
      TTYPE termnode;
      if (termvals) {
        termnode.setFromValue(termvals[i]);
      } else {
        termnode.setFromValue(i);
      }
      nb.d(i) = makeNodeAtLevel(k, termnode.toHandle());
    }

    /*
        Reduce, add redundant as necessary, and set answer
    */
    node_handle node = createReducedNode(-1, nb);
    node = makeNodeAtLevel(getDomain()->getNumVariables() , node); 
    result.set(node, 0);
}



template <class TTYPE>
void mt_forest<TTYPE>::showTerminal(FILE* s, node_handle tnode) const
{
  TTYPE term;
  term.setFromHandle(tnode);
  term.show(s, getRangeType());
}

template <class TTYPE>
void mt_forest<TTYPE>::writeTerminal(FILE* s, node_handle tnode) const
{
  TTYPE term;
  term.setFromHandle(tnode);
  term.write(s, getRangeType());
}

template <class TTYPE>
node_handle mt_forest<TTYPE>::readTerminal(FILE* s)
{
  TTYPE term;
  term.read(s, getRangeType());
  return term.toHandle();
}


};  // end MEDDLY namespace





#else

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

    virtual void createEdgeForVar(int vh, bool pr, const bool* terms, dd_edge& a);
    virtual void createEdgeForVar(int vh, bool pr, const int* terms, dd_edge& a);
    virtual void createEdgeForVar(int vh, bool pr, const float* terms, dd_edge& a);

    
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
    node_handle* getTerminalNodes(int n, const bool* terms);
    node_handle* getTerminalNodes(int n, const int* terms);
    node_handle* getTerminalNodes(int n, const float* terms);

  private:  

    bool counting;

    // scratch pad for buildLevelNode and getTerminalNodes
    node_handle* dptrs;
    int dptrsSize;

};

#endif

#endif

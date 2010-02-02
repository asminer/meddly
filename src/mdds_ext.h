
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



// TODO: Testing
// TODO: mdd_node_manager
// TODO: mxd_node_manager
// TODO: evmdd_node_manager
// TODO: mtmdd_node_manager
// TODO: mtmxd_node_manager (??)

/* 
  TODO: ensure this rule
  All extensions must over-ride either reduceNode() or normalizeAndReduceNode().
  Normalizing is only for edge-valued decision diagrams.
*/
  
#ifndef MDDS_EXTENSIONS_H
#define MDDS_EXTENSIONS_H

#include "mdds.h"

#define SORT_BUILD

// N = domain->getNumVars() + 1
void sortVector(int** indexes, int* terms, int N, int nVars);
void sortMatrix(int** indexes, int** pindexes, int* terms, int N, int nVars);


class mdd_node_manager : public node_manager {
  public:

    mdd_node_manager(domain *d);
    ~mdd_node_manager();

    // Refer to meddly.h
    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(bool val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual error findFirstElement(const dd_edge& f, int* vlist) const;
    
    // Refer to meddly_expert.h
    virtual int reduceNode(int p);

  protected:

    // Create edge representing f(vlist[]) = term and store it in curr.
    void createEdge(const int* vlist, int term, dd_edge& e);

    // Create a node, at level k, whose ith index points to dptr.
    // If i is -1, all indices of the node will point to dptr.
    int createNode(int k, int i, int dptr);

#ifdef SORT_BUILD
    int sortBuild(int** list, int height, int begin, int end);
#endif

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev);

    // Refer to meddly.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist, 
        const float* terms, int N, dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, int &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, float &term) const;
};


class mtmdd_node_manager : public node_manager {
  public:

    mtmdd_node_manager(domain *d, forest::range_type t);
    ~mtmdd_node_manager();

    // Refer to meddly.h
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    
    // Refer to meddly_expert.h
    virtual int reduceNode(int p);

  protected:

    // Create edge representing f(vlist[]) = term and store it in curr.
    void createEdge(const int* vlist, int term, dd_edge& e);

    // Create a node, at level k, whose ith index points to dptr.
    // If i is -1, all indices of the node will point to dptr.
    int createNode(int k, int i, int dptr);

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev);

    // Refer to meddly.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.

    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist, 
        const float* terms, int N, dd_edge &e);
    virtual error createEdge(bool val, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, int &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, float &term) const;
};


class mxd_node_manager : public node_manager {
  // TODO: mxds can only be forest::IDENTITY_REDUCED
  public:

    mxd_node_manager(domain *d);
    ~mxd_node_manager();

    // Refer to meddly.h
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge& e);
    virtual error createEdge(bool val, dd_edge &e);
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
    virtual error findFirstElement(const dd_edge& f, int* vlist,
        int* vplist) const;

    // Refer to meddly_expert.h
    virtual int reduceNode(int p);

  protected:
  
    // Create edge representing vlist[] and store it in curr.
    void createEdge(const int* vlist, const int* vplist, dd_edge& e);

    // Create a node n, at level -k, whose jth index points to dptr.
    // Create a node m, at level +k, whose ith index points to n.
    // If i or j is -1, all indices of n will point to dptr and all of m will
    //    point to n.
    // If i or j is -2, simply returns dptr.
    int createNode(int k, int i, int j, int dptr);

    // Create a node, at level k, whose ith index points to dptr.
    // 0 <= i < level bound
    // Used by createNode(k, i, j, dptr)
    int createNode(int k, int i, int dptr);

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev);

    // Refer to meddly.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.

    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist, 
        const float* terms, int N, dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, int &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, float &term) const;
};


class mtmxd_node_manager : public node_manager {
  // TODO: mtmxds can only be forest::IDENTITY_REDUCED
  public:

    mtmxd_node_manager(domain *d, forest::range_type t);
    ~mtmxd_node_manager();

    // Refer to meddly.h
    error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge& e);
    error createEdge(const int* const* vlist, const int* const* vplist,
        const float* terms, int N, dd_edge& e);

    error createEdge(int val, dd_edge &e);
    error createEdge(float val, dd_edge &e);

    error evaluate(const dd_edge& f, const int* vlist, const int* vplist,
        int &term) const;
    error evaluate(const dd_edge& f, const int* vlist, const int* vplist,
        float &term) const;

    // Refer to meddly_expert.h
    virtual int reduceNode(int p);

  protected:

    // Creates an edge representing v[] vp[] = terminal node (not value),
    // and stores it in e.
    void createEdge(const int* v, const int* vp, int termNode, dd_edge& e);

    // Creates a top-level node representing {-1, -1, ..., -1} = terminal node
    // (not value), and returns it (returned node is already linked to.
    int createEdge(int termNode);

    // Create a node n, at level -k, whose jth index points to dptr.
    // Create a node m, at level +k, whose ith index points to n.
    // If i or j is -1, all indices of n will point to dptr and all of m will
    //    point to n.
    // If i or j is -2, simply returns dptr.
    int createNode(int k, int i, int j, int dptr);

    // Create a node, at level k, whose ith index points to dptr.
    // 0 <= i < level bound
    // Used by createNode(k, i, j, dptr)
    int createNode(int k, int i, int dptr);

#ifdef SORT_BUILD
    int sortBuild(int** list, int** plist, int* tlist,
        int height, int begin, int end);
#endif

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev);

    // Refer to meddly.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.

    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);

    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e);
    virtual error createEdge(bool val, dd_edge &e);

    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
};


int binarySearch(const int* a, int sz, int find);

class evmdd_node_manager : public node_manager {
  public:
    evmdd_node_manager(domain *d, forest::edge_labeling el);
    ~evmdd_node_manager();

    virtual void getDefaultEdgeValue(int& n) const { n = INF; }
    virtual void getIdentityEdgeValue(int& n) const { n = 0; }
    virtual void getDefaultEdgeValue(float& n) const { n = NAN; }
    virtual void getIdentityEdgeValue(float& n) const { n = 1; }

    virtual void initEdgeValues(int p) = 0;

  protected:
    using node_manager::createTempNode;
    virtual int createTempNode(int k, int sz, bool clear = true);

    virtual int doOp(int a, int b) const = 0;
    virtual float doOp(float a, float b) const = 0;

    // Create edge representing vlist[], ending in a terminal value term
    // (or you could think of it as a vector with an edge-value equal to term),
    // and store it in curr.
    template <typename T>
    void createEdgeInternal(const int* vlist, T term, dd_edge &e);

    // Similar to above except that all paths lead to term.
    template <typename T>
    forest::error createEdgeInternal(T term, dd_edge &e);

    // This create a EVMDD from a collection of edges (represented 
    // as vectors).
    template <typename T>
    forest::error createEdgeInternal(const int* const* vlist,
        const T* terms, int N, dd_edge &e);

    // Create a node, at level k, whose ith index points to dptr with an
    // edge value ev. If i is -1, all indices of the node will point to dptr
    // with an edge value ev.
    template <typename T>
    void createNode(int k, int index, int dptr, T ev, int& res, T& resEv);

    // Create a sparse node, at level k, whose ith index points to dptr with
    // edge value ev.
    // 0 <= i < level bound.
    // Used by createNode(k, i, dptr, ev, res, resEv).
    template<typename T>
    void createSparseNode(int k, int index, int dptr, T ev,
        int& res, T& resEv);

    // Create a reduced node, at level k such that ith non-zero edge is
    // represented by {index[i], dptr[i] and ev[i]}.

    virtual void createNode(int lh, std::vector<int>& index,
        std::vector<int>& dptr, std::vector<int>& ev,
        int& result, int& resultEv) { }

    virtual void createNode(int lh, std::vector<int>& index,
        std::vector<int>& dptr, std::vector<float>& ev,
        int& result, float& resultEv) { }

    // If the element represented by vlist exists, return the value
    // corresponding to it.
    template <typename T>
    forest::error evaluateInternal(const dd_edge &f, const int* vlist,
        T &term) const;

    // list and tList are sorted lists containing elements to be added.
    // height is the height of the top variable in each element
    // (usually equal to the number of variables in the domain).
    // begin and end indicate the range of elements within the lists
    // that are to be added.
    // The sum of these edges is returned via node and edgeValue.
#if 1
    template <typename T>
    void sortBuild(int** list, T* tList, int height, int begin, int end,
        int& node, T& edgeValue);
#else
    template <typename T>
    void sortBuild(vector<int*>& list, vector<T>& tList, int height,
        int& node, T& edgeValue);
#endif

    virtual void normalizeAndReduceNode(int& p, int& ev) = 0;
    virtual void normalizeAndReduceNode(int& p, float& ev) = 0;

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual int reduceNode(int p);

    // Refer to meddly.h:
    // The following will either abort or return an error since they are not
    // applicable to this forest.

    virtual error createEdge(const int* const* vlist, const int* terms,
        int N, dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;

    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;

    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const float* terms, int N, dd_edge &e);
    virtual error createEdge(bool val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, int &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, float &term) const;
};


class evplusmdd_node_manager : public evmdd_node_manager {
  public:
    evplusmdd_node_manager(domain *d);
    ~evplusmdd_node_manager();

    virtual int createTempNode(int lh, std::vector<int>& downPointers,
        std::vector<int>& edgeValues);

    virtual void createNode(int lh, std::vector<int>& index,
        std::vector<int>& dptr, std::vector<int>& ev,
        int& result, int& resultEv);

    virtual error createEdge(const int* const* vlist, const int* terms,
        int N, dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;

    virtual forest::error getElement(const dd_edge& a, int index, int* e);
    virtual forest::error getElement(int a, int index, int* e);

    virtual int doOp(int a, int b) const { return a + b; }
    virtual float doOp(float a, float b) const {
      assert(false); return a + b;
    }

    virtual void initEdgeValues(int p);
    virtual bool getDownPtrsAndEdgeValues(int node,
        std::vector<int>& dptrs, std::vector<int>& evs) const;
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev) { }
};


class evtimesmdd_node_manager : public evmdd_node_manager {
  public:
    evtimesmdd_node_manager(domain *d);
    ~evtimesmdd_node_manager();
    
    virtual int createTempNode(int lh, std::vector<int>& downPointers,
        std::vector<float>& edgeValues);

    virtual void createNode(int lh, std::vector<int>& index,
        std::vector<int>& dptr, std::vector<float>& ev,
        int& result, float& resultEv);

    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;

    virtual int doOp(int a, int b) const { assert(false); return a * b; }
    virtual float doOp(float a, float b) const { return a * b; }

    virtual void initEdgeValues(int p);
    virtual bool getDownPtrsAndEdgeValues(int node,
        std::vector<int>& dptrs, std::vector<float>& evs) const;
    virtual void normalizeAndReduceNode(int& p, int& ev) { }
    virtual void normalizeAndReduceNode(int& p, float& ev);
};

























// ------------------------ Inline methods -----------------------------------

template <typename T>
forest::error
evmdd_node_manager::createEdgeInternal(T term, dd_edge &e)
{
  if (e.getForest() != this) return forest::INVALID_OPERATION;

  if (reductionRule == forest::FULLY_REDUCED) {
    e.set(getTerminalNode(true), term, domain::TERMINALS);
  } else {
    // construct the edge bottom-up
    const int* h2l_map = expertDomain->getHeightsToLevelsMap();
    int h_sz = expertDomain->getNumVariables() + 1;
    int prev = getTerminalNode(false);
    T prevEv = 0;
    int curr = getTerminalNode(true);
    T currEv = term;
    for (int i=1; i<h_sz; i++) {
      prev = curr;
      prevEv = currEv;
      createNode(h2l_map[i], -1, prev, prevEv, curr, currEv);
      unlinkNode(prev);
    }
    e.set(curr, currEv, getNodeLevel(curr));
  }
  return forest::SUCCESS;
}


template <typename T>
forest::error
evmdd_node_manager::createEdgeInternal(const int* const* vlist,
    const T* terms, int N, dd_edge &e)
{

  // TODO: modifying this to make buffered updates to the EVMDD
  // TODO: test modifications

  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || terms == 0 || N <= 0) return forest::INVALID_VARIABLE;

  // check if the vlist contains valid indexes
  for (int i = 0; i < N; i++)
  {
    const int* h2l_map = expertDomain->getHeightsToLevelsMap();
    int currHeight = expertDomain->getNumVariables();
    int currLevel = h2l_map[currHeight];
    while (currHeight > 0)
    {
      if (vlist[i][currLevel] >= getLevelSize(currLevel))
        return forest::INVALID_VARIABLE;
      currHeight--;
      currLevel = h2l_map[currHeight];
    }
  }

#ifndef SORT_BUILD
  // build using "standard" procedure
  createEdgeInternal(vlist[0], terms[0], e);
  dd_edge curr(this);
  for (int i=1; i<N; i++) {
    createEdgeInternal(vlist[i], terms[i], curr);
    e += curr;
  }
#else

  if (N < 2) {
    createEdgeInternal(vlist[0], terms[0], e);
  }
  else {
    // check for special cases
    bool specialCasesFound = false;
    for (int i = 0; i < N && !specialCasesFound; i++)
    {
      int* curr = (int*)vlist[i];
      int* end = curr + expertDomain->getNumVariables() + 1;
      for ( ; curr != end; )
      {
        if (*curr++ < 0) {
          specialCasesFound = true;
          break;
        }
      }
    }

    if (specialCasesFound) {
      // build using "standard" procedure
      createEdgeInternal(vlist[0], terms[0], e);
      dd_edge curr(this);
      for (int i=1; i<N; i++) {
        createEdgeInternal(vlist[i], terms[i], curr);
        e += curr;
      }
    }
    else {
      // build using sort-based procedure
#if 0
      int* list[N];
      memcpy(list, vlist, N * sizeof(int*));
      int intList[N];
      memcpy(intList, terms, N * sizeof(int));
      // int* intList = (int*)terms;
      sortVector(list, intList, N, getDomain()->getNumVariables() + 1);
      int node;
      T ev;
      int* temp = (int*)intList;
      T* tList = (T*)temp;
      sortBuild(list, tList, getDomain()->getNumVariables(), 0, N, node, ev);
      e.set(node, ev, getNodeLevel(node));
#else
      // put terms into vlist[i][0]
      int* list[N];
      memcpy(list, vlist, N * sizeof(int*));
      int intList[N];
      memcpy(intList, terms, N * sizeof(int));
      int* temp = (int*)intList;
      T* tList = (T*)temp;
      int node;
      T ev;
      sortBuild(list, tList, getDomain()->getNumVariables(), 0, N, node, ev);
      e.set(node, ev, getNodeLevel(node));
#endif
    }
  }
#endif
  return forest::SUCCESS;
}


template <typename T>
void
evmdd_node_manager::createEdgeInternal(const int* v, T term, dd_edge &e)
{
  // construct the edge bottom-up
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();
  int h_sz = expertDomain->getNumVariables() + 1;
  int prev = getTerminalNode(false);
  T prevEv;
  getIdentityEdgeValue(prevEv);
  int curr = getTerminalNode(true);
  T currEv = term;
  for (int i=1; i<h_sz; i++) {
    prev = curr;
    prevEv = currEv;
    createNode(h2l_map[i], v[h2l_map[i]], prev, prevEv, curr, currEv);
    unlinkNode(prev);
  }
  e.set(curr, currEv, getNodeLevel(curr));
}

template <typename T>
void evmdd_node_manager::createNode(int k, int index, int dptr, T ev,
    int& res, T& resEv)
{
  DCASSERT(dptr >= -1 && ev >= 0);
  DCASSERT(index >= -1);

  if (index == -1) {
    // if all edge values are the same for a evmdd node, it is
    // equivalent to say that they are all 0 (since we subtract the minimum
    // anyway).
    if (reductionRule == forest::FULLY_REDUCED) {
      res = sharedCopy(dptr);
      resEv = ev;
      return;
    }
    res = createTempNodeMaxSize(k, false);
    setAllDownPtrsWoUnlink(res, dptr);
    setAllEdgeValues(res, ev);
    normalizeAndReduceNode(res, resEv);
    return;
  }

  // a single downpointer points to dptr
  if (nodeStorage == FULL_STORAGE ||
      (nodeStorage == FULL_OR_SPARSE_STORAGE && index < 2)) {
    // Build a full node
    res = createTempNode(k, index + 1);
    setDownPtrWoUnlink(res, index, dptr);
    setEdgeValue(res, index, ev);
    normalizeAndReduceNode(res, resEv);
  }
  else {
    DCASSERT (nodeStorage == SPARSE_STORAGE ||
        (nodeStorage == FULL_OR_SPARSE_STORAGE && index >= 2));
    // Build a sparse node
    createSparseNode(k, index, dptr, ev, res, resEv);
    // res, resEv set by createSparseNode(..); do not call normalizeAndReduce
  }
}


template <typename T>
void evmdd_node_manager::createSparseNode(int k, int index,
    int dptr, T ev, int& res, T& resEv)
{
  DCASSERT(k != 0);

  if (isTimeToGc()) { garbageCollect(); }

  DCASSERT(isValidLevel(k));
  CHECK_RANGE(0, index, getLevelSize(k));

  // get a location in address[] to store the node
  int p = getFreeNode(k);

#ifdef DEBUG_MDD_H
  printf("%s: k: %d, index: %d, new p: %d\n", __func__, k, index, p);
  fflush(stdout);
#endif

  // fill in the location with p's address info
  address[p].level = k;
  address[p].offset = getHole(k, 7 /*4 + 3*/, true);
  address[p].cache_count = 0;

#ifdef DEBUG_MDD_H
  printf("%s: offset: %d\n", __func__, address[p].offset);
  fflush(stdout);
#endif

  int* foo = level[mapLevel(k)].data + address[p].offset;
  foo[0] = 1;                     // #incoming
  foo[1] = getTempNodeId();
  foo[2] = -1;                    // size
  foo[3] = index;                 // index
  foo[4] = sharedCopy(dptr);      // downpointer
  T identEv;
  getIdentityEdgeValue(identEv);
  foo[5] = toInt(identEv);        // this is the only ev, set resEv = ev
  foo[6] = p;                     // pointer to this node in the address array

  resEv = ev;

#ifdef TRACK_DELETIONS
  cout << "Creating node " << p << "\n";
  cout.flush();
#endif

  incrNodesActivatedSinceGc();

  // search in unique table
  int q = find(p);
  if (getNull() == q) {
    // no duplicate found; insert into unique table
    insert(p);
    DCASSERT(getCacheCount(p) == 0);
    DCASSERT(find(p) == p);
    res = p;
  }
  else {
    // duplicate found; discard this node and return the duplicate
    unlinkNode(dptr);
    // code from deleteTempNode(p) adapted to work here
    {
      makeHole(k, getNodeOffset(p), 7);
      freeNode(p);
      if (level[mapLevel(k)].compactLevel) compactLevel(k);
    }
    res = sharedCopy(q);
  }
}


template <typename T>
forest::error evmdd_node_manager::evaluateInternal(const dd_edge &f,
    const int* vlist, T &term) const
{
  if (f.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0) return forest::INVALID_VARIABLE;

  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  int node = f.getNode();
  const int* h2l_map = expertDomain->getHeightsToLevelsMap();

  f.getEdgeValue(term);
  while (!isTerminalNode(node)) {
    T ev;
    getDefaultEdgeValue(ev);
    int n = 0;
    if (isFullNode(node)) {
      int index = vlist[h2l_map[getNodeHeight(node)]];
      if (index < getFullNodeSize(node)) {
        getFullNodeEdgeValue(node, index, ev);
        n = getFullNodeDownPtr(node, index);
      }
    } else {
      // find the index
      // binary search
      int index = binarySearch(getSparseNodeIndexes(node),
          getSparseNodeSize(node), vlist[h2l_map[getNodeHeight(node)]]);
      if (index != -1 ) {
        getSparseNodeEdgeValue(node, index, ev);
        n = getSparseNodeDownPtr(node, index);
      }
    }
    term = (n == 0)? ev: doOp(term, ev);
    node = n;
  }
  return forest::SUCCESS;
}


#ifdef SORT_BUILD

#if 0

template <typename T>
void evmdd_node_manager::sortBuild(int** list, T* tList,
    int height, int begin, int end, int& node, T& edgeValue)
{
  // [begin, end)

  // terminal condition
  if (height == 0)
  {
    node = -1;
    edgeValue = tList[begin++];
    while (begin != end) { edgeValue += tList[begin++]; }
    return;
  }

  int** currList = list;
  int nextHeight = height - 1;
  int level = expertDomain->getVariableWithHeight(height);

  vector<int> indices;
  vector<int> dptrs;
  vector<T> edgeValues;
  int currBegin = begin;
  int currEnd = currBegin;
  for ( ; currEnd < end; currBegin = currEnd) 
  {
    int currIndex = currList[currBegin][level];
    assert(currIndex >= 0);
    for (currEnd = currBegin + 1;
        currEnd < end && currIndex == currList[currEnd][level];
        ++currEnd);
    // found new range
    // to be "unioned" and assigned to node[currIndex]
#ifdef DEBUG_SORT_BUILD
    printf("level: %d, currIndex: %d, currBegin: %d, currEnd: %d\n",
        level, currIndex, currBegin, currEnd);
    fflush(stdout);
#endif
    int n;
    T ev;
    sortBuild(list, tList, nextHeight, currBegin, currEnd, n, ev);
#ifdef DEBUG_SORT_BUILD
    printf("level: %d, currIndex: %d, currBegin: %d, currEnd: %d\n",
        level, currIndex, currBegin, currEnd);
    fflush(stdout);
#endif
    indices.push_back(currIndex);
    dptrs.push_back(n);
    edgeValues.push_back(ev);
  }

  assert(dptrs.size() >= 1);
  if (dptrs.size() == 1) {
    // single entry: store directly as sparse
    // TODO: this should be able to accept many more cases
    createNode(level, indices[0], dptrs[0], edgeValues[0], node, edgeValue);
    unlinkNode(dptrs[0]);
    return;
  }
  // full node
  int size = indices[indices.size() - 1] + 1;
  node = createTempNode(level, size);
  vector<int>::iterator indexIter = indices.begin();
  vector<int>::iterator dptrIter = dptrs.begin();
  typename vector<T>::iterator evIter = edgeValues.begin();
  for ( ; indexIter != indices.end(); indexIter++, dptrIter++, evIter++)
  {
    if (getDownPtr(node, *indexIter) != 0) {
      for (unsigned i = 0u; i < dptrs.size(); i++)
      {
        printf("%d:%d ", indices[i], dptrs[i+1]);
      }
      printf("\n");
      for (int i = begin; i < end; i++)
      {
        printf("%d:%d ", i, currList[i][level]);
      }
      printf("\n");
      exit(1);
    }
    setDownPtrWoUnlink(node, *indexIter, *dptrIter);
    setEdgeValue(node, *indexIter, *evIter);
    unlinkNode(*dptrIter);
  }
  normalizeAndReduceNode(node, edgeValue);
}

#else

template <typename T>
void evmdd_node_manager::sortBuild(int** list, T* tList,
    int height, int begin, int end, int& node, T& edgeValue)
{
  // [begin, end)

  // terminal condition
  if (height == 0)
  {
    node = -1;
    edgeValue = tList[begin++];
    while (begin != end) { edgeValue += tList[begin++]; }
    return;
  }

  int N = end - begin;
  int level = expertDomain->getVariableWithHeight(height);
  int nextHeight = height - 1;

  if (N == 1) {
    // nothing to sort; just build a node starting at this level
    int index = list[begin][level];
    int n;
    T ev;
    sortBuild(list, tList, nextHeight, begin, end, n, ev);
    createNode(level, index, n, ev, node, edgeValue);
    unlinkNode(n);
    return;
  }

  // do radix sort for this level

  int levelSize = getLevelSize(level);
#if 0
  vector<int> count(levelSize+1, 0);
#else
  vector<int> count(1, 0);
#endif

  // Curly braces here to limit the scope of the vectors defined within.
  // Without limiting the scope, memory usage will increase significantly
  // -- especially if there are a lot of variables in the domain.
  {
    vector<int*> sortedList(N, (int*)0);
    vector<T> sortedtList(N, 0);

    // determine size for count[]
    levelSize = 0;
    for (int i = begin; i < end; i++) {
      int index = list[i][level];
      if (index > levelSize) { levelSize = index; }
    }
    // levelSize refers to the maximum index found so far,
    // add 1 to convert to maximum size.
    levelSize++;
    // an extra space is needed at the end for the radix sort algorithm
    count.resize(levelSize+1, 0);

    // go through list and count the number of entries in each "bucket"
    for (int i = begin; i < end; i++) { count[list[i][level]]++; }

    // find starting index for each "bucket" in sorted lists
    // levelSize == number of buckets
#if 0
    for (int i = 0; i < levelSize; i++)
    {
      count[i+1] += count[i];
      count[i]--;
    }
    count[levelSize]--;
#else
    vector<int>::iterator last = count.begin() + levelSize;
    for (vector<int>::iterator iter = count.begin(); iter != last; iter++)
    {
      *(iter+1) += *iter;
      (*iter)--;
    }
    (*last)--;
#endif

    // insert into correct positions in sorted lists
    // go from last to first to preserve order
#if 0
    for (int i = end - 1; i >= begin; i--)
    {
      // getting index and get count[] ready for next insert
      int index = count[list[i][level]];
      count[list[i][level]]--;
      // insert at index
      sortedList[index] = list[i];
      sortedtList[index] = tList[i];
    }
#else
    int** listPtr = list + end - 1;
    int** firstListPtr = list + begin - 1;
    T* tListPtr = tList + end - 1;
    for ( ; listPtr != firstListPtr; )
    {
      // getting index and get count[] ready for next insert
      int index = count[(*listPtr)[level]]--;
      // insert at index
      sortedList[index] = *listPtr--;
      sortedtList[index] = *tListPtr--;
    }
#endif

    // write sorted lists to the original lists
#if 0
    for (int i = 0; i < N; i++)
    {
      list[begin + i] = sortedList[i];
      tList[begin + i] = sortedtList[i];
    }
#else
    listPtr = list + begin;
    tListPtr = tList + begin;
    vector<int*>::iterator sortedListIter = sortedList.begin();
    typename vector<T>::iterator sortedtListIter = sortedtList.begin();
    for ( ; sortedListIter != sortedList.end(); )
    {
      *listPtr++ = *sortedListIter++;
      *tListPtr++ = *sortedtListIter++;
    }
#endif
  }

  // after insertion, range for bucket[i] is [count[i]+1, count[i+1]+1)

  // call sortBuild for each index and store result
  // as (index, node, edge-value).
  vector<int> indices;
  vector<int> dptrs;
  vector<T> edgeValues;
  for (int i = 0; i < levelSize; i++)
  {
    if (count[i+1] > count[i]) {
      int n;
      T ev;
      sortBuild(list, tList, nextHeight, begin + count[i] + 1,
          begin + count[i+1] + 1, n, ev);
      indices.push_back(i);
      dptrs.push_back(n);
      edgeValues.push_back(ev);
    }
  }

  // build node from indices, dptrs and edgeValues

  DCASSERT(dptrs.size() > 0);

#if 0
  if (dptrs.size() == 1) {
    createNode(level, indices[0], dptrs[0], edgeValues[0], node, edgeValue);
    unlinkNode(dptrs[0]);
    return;
  }

  // build full node
  int size = indices[indices.size() - 1] + 1;
  node = createTempNode(level, size);
  vector<int>::iterator indexIter = indices.begin();
  vector<int>::iterator dptrIter = dptrs.begin();
  typename vector<T>::iterator evIter = edgeValues.begin();
  for ( ; indexIter != indices.end(); )
  {
    setDownPtrWoUnlink(node, *indexIter, *dptrIter);
    setEdgeValue(node, *indexIter++, *evIter++);
    unlinkNode(*dptrIter++);
  }
  normalizeAndReduceNode(node, edgeValue);
#else
  createNode(level, indices, dptrs, edgeValues, node, edgeValue);
#endif
}

#endif

#endif

#endif


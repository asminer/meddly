
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

int binarySearch(const int* a, int sz, int find);

template <typename T>
T handleMultipleTerminalValues(T* tList, int begin, int end);
template <>
bool handleMultipleTerminalValues(bool* tList, int begin, int end);

/*
TODO: modifying mtmxd_node_manager to work like mtmdd_node_manager; and
mxd_node_manager to work like mdd_node_manager.
TODO: adding real-valued terminals to mtmdd_node_manager.
*/

class mtmdd_node_manager : public node_manager {
  public:

    mtmdd_node_manager(domain *d, forest::range_type t);
    ~mtmdd_node_manager();

    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);

    virtual error createEdge(int val, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);

    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;

    virtual error findFirstElement(const dd_edge& f, int* vlist) const;

    virtual int reduceNode(int p);

    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev);
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

  protected:

    // Used by derived classes for initialization
    mtmdd_node_manager(domain *d, bool relation, forest::range_type t,
        forest::edge_labeling e, forest::reduction_rule r,
        forest::node_storage s, forest::node_deletion_policy dp);

    // This create a MTMDD from a collection of edges (represented 
    // as vectors).
    template <typename T>
    forest::error createEdgeInternal(const int* const* vlist,
        const T* terms, int N, dd_edge &e);

    // Create edge representing f(vlist[]) = term and store it in curr.
    void createEdge(const int* vlist, int term, dd_edge& e);

    // Create a node, at level k, whose ith index points to dptr.
    // If i is -1, all indices of the node will point to dptr.
    int createNode(int k, int i, int dptr);

    // Create a node at level k, such that dptr[i] is the downpointer
    // corresponding to index[i].
    // Important note: no node linking performed because the links
    // are "transferred" from the dptr vector.
    int createNode(int lh, std::vector<int>& index, std::vector<int>& dptr);

    // Creates an edge representing the terminal node given by
    // terminalNode.
    // Note: terminalNode is usually a terminal value converted into
    // the equivalent node representation. This makes this method useful
    // for createEdge(int, e) or (float, e) or (bool, e) as long as the
    // terminal value can be converted into an equivalent node.
    error createEdgeHelper(int terminalNode, dd_edge& e);

    // Get the terminal node at the bottom of the edge with root n
    // and vlist representing the indexes for the levels.
    // Used by evaluate()
    int getTerminalNodeForEdge(int n, const int* vlist) const;

    template <typename T>
    int sortBuild(int** list, T* tList, int height, int begin, int end);

};


class mdd_node_manager : public mtmdd_node_manager {
  public:

    mdd_node_manager(domain *d);
    ~mdd_node_manager();

    using mtmdd_node_manager::createEdge;
    using mtmdd_node_manager::evaluate;

    // Refer to meddly.h
    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(bool val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;

    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
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

    virtual error findFirstElement(const dd_edge& f, int* vlist, int* vplist)
      const;

    virtual int reduceNode(int p);

    bool singleNonZeroAt(int p, int val, int index) const;

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

  protected:

    // Used by derived classes for initialization
    mtmxd_node_manager(domain *d, bool relation, forest::range_type t,
        forest::edge_labeling e, forest::reduction_rule r,
        forest::node_storage s, forest::node_deletion_policy dp);

    // This create a MTMXD from a collection of edges (represented 
    // as vectors vlist and vplist).
    template <typename T>
      forest::error createEdgeInternal(const int* const* vlist,
          const int* const* vplist, const T* terms, int N, dd_edge &e);

    // Creates an edge representing v[] vp[] = terminal node (not value),
    // and stores it in e.
    void createEdge(const int* v, const int* vp, int termNode, dd_edge& e);

    // Starting at height given by {startAtHeight, primedLevel},
    // creates an edge representing v[] vp[] = terminal node (not value),
    // and stores it in e.
    void createEdge(const int* v, const int* vp, int termNode,
        int startAtHeight, bool primedLevel, dd_edge& e);

    // Creates a top-level node representing {-1, -1, ..., -1} = terminal node
    // (not value), and returns it (returned node is already linked to.
    int createEdge(int termNode);

    // Starting at height given by {startAtHeight, primedLevel},
    // creates a node representing v[] vp[] = terminal node (not value)
    // and returns it. Used by createEdge().
    int createNode(const int* v, const int* vp, int termNode,
        int startAtHeight, bool primedLevel);

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

    // Get the terminal node at the bottom of the edge with root n
    // and vlist and vplist representing the indexes for the levels.
    // Used by evaluate()
    int getTerminalNodeForEdge(int n, const int* vlist, const int* vplist)
      const;

    template <typename T>
      int sort(int** list, int** otherList, T* tList,
          int absLevel, int begin, int end);

    template <typename T>
      int sortBuild(int** unpList, int** pList, T* tList,
          int height, int begin, int end);

    void addToTree(int* unp, int* p, int terminalNode);
    int convertTreeToMtMxd();
    int convertToMtMxd(int addr, int height);

    int root;
};


class mxd_node_manager : public mtmxd_node_manager {
  // TODO: mxds can only be forest::IDENTITY_REDUCED
  public:

    mxd_node_manager(domain *d);
    ~mxd_node_manager();

    using mtmxd_node_manager::createEdge;
    using mtmxd_node_manager::evaluate;

    // Refer to meddly.h
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge& e);
    virtual error createEdge(bool val, dd_edge &e);
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;

    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const float* terms, int N, dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist,
        const int* vplist, int &term) const;
    virtual error evaluate(const dd_edge &f, const int* vlist,
        const int* vplist, float &term) const;
};


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
    template <typename T>
    void sortBuild(int** list, T* tList, int height, int begin, int end,
        int& node, T& edgeValue);

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
  if (e.getForest() != this) return forest::INVALID_OPERATION;
  if (vlist == 0 || terms == 0 || N <= 0) return forest::INVALID_VARIABLE;

  // check if the vlist contains valid indexes
  bool specialCasesFound = false;
  for (int i = 0; i < N; i++)
  {
    const int* h2l_map = expertDomain->getHeightsToLevelsMap();
    int currHeight = expertDomain->getNumVariables();
    int currLevel = h2l_map[currHeight];
    while (currHeight > 0)
    {
      int bound = vlist[i][currLevel] + 1;
      if (bound >= getLevelSize(currLevel))
        expertDomain->enlargeVariableBound(currLevel, false, bound);
      else if (bound < 0)
        specialCasesFound = true;
      currHeight--;
      currLevel = h2l_map[currHeight];
    }
  }

  if (N == 1 || specialCasesFound) {
    createEdgeInternal(vlist[0], terms[0], e);
    if (N > 1) {
      dd_edge curr(this);
      for (int i=1; i<N; i++) {
        createEdgeInternal(vlist[i], terms[i], curr);
        e += curr;
      }
    }
  }
  else {
    // build using sort-based procedure
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
  }
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


template <typename T>
void evmdd_node_manager::sortBuild(int** list, T* tList,
    int height, int begin, int end, int& node, T& edgeValue)
{
  // [begin, end)

  // terminal condition
  if (height == 0)
  {
    node = -1;
#if 0
    edgeValue = tList[begin++];
    while (begin != end) { edgeValue += tList[begin++]; }
#else
    edgeValue = handleMultipleTerminalValues(tList, begin, end);
#endif
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
  createNode(level, indices, dptrs, edgeValues, node, edgeValue);
}


template <typename T>
forest::error
mtmdd_node_manager::createEdgeInternal(const int* const* vlist,
    const T* terms, int N, dd_edge &e)
{
  // check if the vlist contains valid indexes
  bool specialCasesFound = false;
  for (int i = 0; i < N; i++)
  {
    const int* h2l_map = expertDomain->getHeightsToLevelsMap();
    int currHeight = expertDomain->getNumVariables();
    int currLevel = h2l_map[currHeight];
    while (currHeight > 0)
    {
      int bound = vlist[i][currLevel] + 1;
      if (bound >= getLevelSize(currLevel))
        expertDomain->enlargeVariableBound(currLevel, false, bound);
      else if (bound < 0)
        specialCasesFound = true;
      currHeight--;
      currLevel = h2l_map[currHeight];
    }
  }

  if (N == 1 || specialCasesFound) {
    // build using "standard" procedure
    if (terms == 0) {
      int trueNode = getTerminalNode(true);
      createEdge(vlist[0], trueNode, e);
      if (N > 1) {
        dd_edge curr(this);
        for (int i=1; i<N; i++) {
          createEdge(vlist[i], trueNode, curr);
          e += curr;
        }
      }
    }
    else {
      createEdge(vlist[0], getTerminalNode(terms[0]), e);
      if (N > 1) {
        dd_edge curr(this);
        for (int i=1; i<N; i++) {
          createEdge(vlist[i], getTerminalNode(terms[i]), curr);
          e += curr;
        }
      }
    }
  }
  else {
    // build using sort-based procedure
    // put terms into vlist[i][0]
    int* list[N];
    memcpy(list, vlist, N * sizeof(int*));

    int result = 0;
    if (terms == 0) {
      result =
        sortBuild(list, (T*)0, getDomain()->getNumVariables(), 0, N);
    } else {
      T tList[N];
      memcpy((void*)tList, (void*)terms, N * sizeof(T));
      result =
        sortBuild(list, tList, getDomain()->getNumVariables(), 0, N);
    }

    e.set(result, 0, getNodeLevel(result));
  }

  return forest::SUCCESS;
}


template <typename T>
int mtmdd_node_manager::sortBuild(int** list, T* tList,
    int height, int begin, int end)
{
  // [begin, end)

  // terminal condition
  if (height == 0)
  {
    return getTerminalNode(handleMultipleTerminalValues(tList, begin, end));
  }

  int N = end - begin;
  int level = expertDomain->getVariableWithHeight(height);
  int nextHeight = height - 1;

  if (N == 1) {
    // nothing to sort; just build a node starting at this level
    int n = sortBuild(list, tList, nextHeight, begin, end);
    int index = list[begin][level];
    int result = createNode(level, index, n);
    unlinkNode(n);
    return result;
  }

  // do radix sort for this level
  int levelSize = 0;
  vector<int> count(1, 0);

  // Curly braces here to limit the scope of the vectors defined within.
  // Without limiting the scope, memory usage will increase significantly
  // -- especially if there are a lot of variables in the domain.
  if (tList != 0) {
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
    vector<int>::iterator last = count.begin() + levelSize;
    for (vector<int>::iterator iter = count.begin(); iter != last; iter++)
    {
      *(iter+1) += *iter;
      (*iter)--;
    }
    (*last)--;

    // insert into correct positions in sorted lists
    // go from last to first to preserve order
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

    // write sorted lists to the original lists
    listPtr = list + begin;
    tListPtr = tList + begin;
    vector<int*>::iterator sortedListIter = sortedList.begin();
    typename vector<T>::iterator sortedtListIter = sortedtList.begin();
    for ( ; sortedListIter != sortedList.end(); )
    {
      *listPtr++ = *sortedListIter++;
      *tListPtr++ = *sortedtListIter++;
    }
  }
  else {
    // same as tList != 0, except that there is no tList to deal with

    vector<int*> sortedList(N, (int*)0);

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
    vector<int>::iterator last = count.begin() + levelSize;
    for (vector<int>::iterator iter = count.begin(); iter != last; iter++)
    {
      *(iter+1) += *iter;
      (*iter)--;
    }
    (*last)--;

    // insert into correct positions in sorted lists
    // go from last to first to preserve order
    int** listPtr = list + end - 1;
    int** firstListPtr = list + begin - 1;
    for ( ; listPtr != firstListPtr; )
    {
      // getting index and get count[] ready for next insert
      int index = count[(*listPtr)[level]]--;
      // insert at index
      sortedList[index] = *listPtr--;
    }

    // write sorted lists to the original lists
    listPtr = list + begin;
    for (vector<int*>::iterator sortedListIter = sortedList.begin();
        sortedListIter != sortedList.end(); )
    {
      *listPtr++ = *sortedListIter++;
    }
  }

  // after insertion, range for bucket[i] is [count[i]+1, count[i+1]+1)

  // call sortBuild for each index and store result as (index, node).
  vector<int> indices;
  vector<int> dptrs;
  for (int i = 0; i < levelSize; i++)
  {
    if (count[i+1] > count[i]) {
      int n = sortBuild(list, tList, nextHeight, begin + count[i] + 1,
          begin + count[i+1] + 1);
      indices.push_back(i);
      dptrs.push_back(n);
    }
  }

  // build node from indices, dptrs and edgeValues

  DCASSERT(dptrs.size() > 0);

  return createNode(level, indices, dptrs);
}


template <typename T>
forest::error
mtmxd_node_manager::createEdgeInternal(const int* const* vlist,
    const int* const* vplist, const T* terms, int N, dd_edge &e)
{
  // check if the vlist contains valid indexes
  bool specialCasesFound = false;
  for (int i = 0; i < N; i++)
  {
    for (int currLevel = expertDomain->getTopVariable();
        currLevel != domain::TERMINALS;
        currLevel = expertDomain->getVariableBelow(currLevel))
    {
      // unprimed level
      int bound = vlist[i][currLevel] + 1;
      if (bound >= expertDomain->getVariableBound(currLevel, false))
        expertDomain->enlargeVariableBound(currLevel, false, bound);
      else if (bound < 0)
        specialCasesFound = true;
      // primed level
      bound = vplist[i][currLevel] + 1;
      if (bound >= expertDomain->getVariableBound(currLevel, true))
        expertDomain->enlargeVariableBound(currLevel, true, bound);
      else if (bound < 0)
        specialCasesFound = true;
    }
  }

  if (N == 1 || specialCasesFound) {
    // build using "standard" procedure
    if (terms == 0) {
      int trueNode = getTerminalNode(true);
      createEdge(vlist[0], vplist[0], trueNode, e);
      if (N > 1) {
        dd_edge curr(this);
        for (int i=1; i<N; i++) {
          createEdge(vlist[i], vplist[i], trueNode, curr);
          e += curr;
        }
      }
    }
    else {
      createEdge(vlist[0], vplist[0], getTerminalNode(terms[0]), e);
      if (N > 1) {
        dd_edge curr(this);
        for (int i=1; i<N; i++) {
          createEdge(vlist[i], vplist[i], getTerminalNode(terms[i]), curr);
          e += curr;
        }
      }
    }
  }
  else {
#if 1
    int result = 0;
    // build using sort-based procedure
    // put terms into vlist[i][0]
    int* list[N];
    memcpy(list, vlist, N * sizeof(int*));
    int* pList[N];
    memcpy(pList, vplist, N * sizeof(int*));
    if (terms == 0) {
      result =
        sortBuild(list, pList, (T*)0, getDomain()->getNumVariables(), 0, N);
    } else {
      T tList[N];
      memcpy((void*)tList, (void*)terms, N * sizeof(T));
      result =
        sortBuild(list, pList, tList, getDomain()->getNumVariables(), 0, N);
    }
#else
    if (terms == 0) {
      int terminalNode = getTerminalNode(
          handleMultipleTerminalValues((T*)0, 0, N));
      for (int i = 0; i < N; ++i)
      {
        addToTree((int*)vlist[i], (int*)vplist[i], terminalNode);
      }
    } else {
      for (int i = 0; i < N; ++i)
      {
        addToTree((int*)vlist[i], (int*)vplist[i], getTerminalNode(terms[i]));
      }
    }
    int result = convertTreeToMtMxd();
#endif
    e.set(result, 0, getNodeLevel(result));
  }

  return forest::SUCCESS;
}


template <typename T>
int mtmxd_node_manager::sort(int** list, int** otherList, T* tList,
    int absLevel, int begin, int end)
{
  DCASSERT(tList != 0);

  static vector<int*> sortedList;
  static vector<int*> sortedOtherList;
  static vector<T> sortedtList;
  static vector<int> count;

  int N = end - begin;
  int levelSize = 0;

  if (int(sortedList.size()) < N) {
    sortedList.resize(N);
    sortedOtherList.resize(N);
    sortedtList.resize(N);
  }

  // determine size for count[]
  for (int i = begin; i < end; i++) {
    int index = list[i][absLevel];
    if (index > levelSize) { levelSize = index; }
  }
  // levelSize refers to the maximum index found so far,
  // add 1 to convert to maximum size.
  levelSize++;
  // an extra space is needed at the end for the radix sort algorithm
  if (int(count.size()) < 1 + levelSize) {
    count.resize(levelSize + 1);
  }
  fill_n(count.begin(), 1 + levelSize, 0);

  // go through list and count the number of entries in each "bucket"
  for (int i = begin; i < end; i++) { count[list[i][absLevel]]++; }

  // find starting index for each "bucket" in sorted lists
  // levelSize == number of buckets
  vector<int>::iterator last = count.begin() + levelSize;
  for (vector<int>::iterator iter = count.begin(); iter != last; ++iter)
  {
    *(iter+1) += *iter;
    --*iter;
  }
  --*last;

  // insert into correct positions in sorted lists
  // go from last to first to preserve order
  int** listPtr = list + end - 1;
  int** firstListPtr = list + begin - 1;
  int** otherListPtr = otherList + end - 1;
  T* tListPtr = tList + end - 1;
  for ( ; listPtr != firstListPtr; )
  {
    // getting index and get count[] ready for next insert
    int index = count[(*listPtr)[absLevel]]--;
    // insert at index
    sortedList[index] = *listPtr--;
    sortedOtherList[index] = *otherListPtr--;
    sortedtList[index] = *tListPtr--;
  }

  // write sorted lists to the original lists
  memcpy(list + begin, &sortedList[0], N * sizeof(int*));
  memcpy(otherList + begin, &sortedOtherList[0], N * sizeof(int*));
  memcpy(tList + begin, &sortedtList[0], N * sizeof(T));
  return levelSize;
}


template <>
inline
int mtmxd_node_manager::sort(int** list, int** otherList, bool* tList,
    int absLevel, int begin, int end)
{
  int N = end - begin;
  int levelSize = 0;

  // same as tList != 0, except that there is no tList to deal with

  static vector<int*> sortedList;
  static vector<int*> sortedOtherList;
  static vector<int> count;

  if (int(sortedList.size()) < N) {
    sortedList.resize(N);
    sortedOtherList.resize(N);
  }

  // determine size for count[]
  for (int i = begin; i < end; i++) {
    int index = list[i][absLevel];
    if (index > levelSize) { levelSize = index; }
  }
  // levelSize refers to the maximum index found so far,
  // add 1 to convert to maximum size.
  levelSize++;
  // an extra space is needed at the end for the radix sort algorithm
  if (int(count.size()) < 1 + levelSize) {
    count.resize(levelSize + 1);
  }
  fill_n(count.begin(), 1 + levelSize, 0);

  // go through list and count the number of entries in each "bucket"
  for (int i = begin; i < end; i++) { count[list[i][absLevel]]++; }

  // find starting index for each "bucket" in sorted lists
  // levelSize == number of buckets
  vector<int>::iterator last = count.begin() + levelSize;
  for (vector<int>::iterator iter = count.begin(); iter != last; ++iter)
  {
    *(iter+1) += *iter;
    --*iter;
  }
  --*last;

  // insert into correct positions in sorted lists
  // go from last to first to preserve order
  int** listPtr = list + end - 1;
  int** firstListPtr = list + begin - 1;
  int** otherListPtr = otherList + end - 1;
  for ( ; listPtr != firstListPtr; )
  {
    // getting index and get count[] ready for next insert
    int index = count[(*listPtr)[absLevel]]--;
    // insert at index
    sortedList[index] = *listPtr--;
    sortedOtherList[index] = *otherListPtr--;
  }

  // write sorted lists to the original lists
  memcpy(list + begin, &sortedList[0], N * sizeof(int*));
  memcpy(otherList + begin, &sortedOtherList[0], N * sizeof(int*));
  return levelSize;
}


template <typename T>
int mtmxd_node_manager::sortBuild(int** unpList, int** pList, T* tList,
    int height, int begin, int end)
{
  // [begin, end)

  // terminal condition
  if (height == 0) {
    return getTerminalNode(handleMultipleTerminalValues(tList, begin, end));
  }

  if (begin + 1 == end) {
    return createNode(unpList[begin], pList[begin],
        getTerminalNode(handleMultipleTerminalValues(tList, begin, end)),
        ABS(height), height < 0);
  }

  int** list = 0;
  int** otherList = 0;
  int nextHeight = 0;
  int level = 0;
  if (height > 0) {
    list = unpList;
    otherList = pList;
    nextHeight = -height;
    level = expertDomain->getVariableWithHeight(height);
  } else {
    list = pList;
    otherList = unpList;
    nextHeight = -height-1;
    level = -(expertDomain->getVariableWithHeight(-height));
  }
  int absLevel = level < 0? -level: level;

  int levelSize = sort(list, otherList, tList, absLevel, begin, end);

#ifdef DEVELOPMENT_CODE
  // find largest index
  int largestIndex = list[begin][absLevel];
  for (int i = begin + 1; i < end; )
  {
    DCASSERT(largestIndex <= list[i][absLevel]);
    largestIndex = list[i][absLevel];
    // skip the elements with the same index at this level
    for (++i; i < end && list[i][absLevel] == largestIndex; ++i);
  }
  if (largestIndex + 1 != levelSize) {
    printf("largest index: %d, levelSize: %d\n", largestIndex, levelSize);
    assert(false);
  }
#endif

  // build node
  int result = createTempNode(level, levelSize, true);
  int* ptr = getFullNodeDownPtrs(result);
  for (int i = begin; i < end; )
  {
    int index = list[i][absLevel];
    int start = i++;
    // skip the elements with the same index at this level
    for ( ; i < end && list[i][absLevel] == index; ++i);
    ptr[index] = sortBuild(unpList, pList, tList, nextHeight, start, i);
  }

  return reduceNode(result);
}


inline
int mtmxd_node_manager::convertToMtMxd(int addr, int height)
{
  DCASSERT(height != 0);

  int level = height > 0
    ? expertDomain->getVariableWithHeight(height)
    : -expertDomain->getVariableWithHeight(-height);
  int* ptr = getAddress(level, addr);
  DCASSERT(ptr[0] >= 5);

  // find largest index
  int largestIndex = -1;
  for (int* i = ptr + 1 + *ptr, *begin = ptr + 1; i != begin; )
  {
    --i;
    if (*i != 0) { largestIndex = i - begin; break; }
  }

  DCASSERT(largestIndex != -1);

  // build node at this level
  int result = createTempNode(level, largestIndex + 1, false);
  int* resultPtr = getFullNodeDownPtrs(result);

  // there is a chance that creating the temp node at the same level
  // will invalidate ptr[], so reset this pointer.
  ptr = getAddress(level, addr);

  // non-terminal case: build children
  if (height != -1) {
    int nextHeight = height > 0? -height: -height-1;
    for (int* i = ptr + 1, *end = i + largestIndex + 1; i != end; ++i)
    {
      *resultPtr++ = (0 == *i? 0: convertToMtMxd(*i, nextHeight));
    }
  } else {
    for (int* i = ptr + 1, *end = i + largestIndex + 1; i != end; )
    {
      *resultPtr++ = *i++;
    }
  }

  makeHole(level, addr, *ptr + 1);

  return reduceNode(result);
}


inline
int mtmxd_node_manager::convertTreeToMtMxd()
{
  int result =
    root == 0
    ? 0
    : convertToMtMxd(root, expertDomain->getNumVariables());
  root = 0;
  return result;
}


inline
void mtmxd_node_manager::addToTree(int* unpList, int* pList, int terminalNode)
{
  if (terminalNode == 0) return;

  // Tree starts at root. root stores an offset which together with
  // the root level (usually the top variable in the domain),
  // represents a unique address.
  // All "nodes" store offsets to their "downpointers", except for
  // the nodes at level -1 (i.e. the first primed level above terminals).
  // The nodes at level -1 store the actual terminal nodes.

  // getHole() will not accept a hole size smaller than 5.
  const int minNodeSize = 5;

  const int* h2lMap = expertDomain->getHeightsToLevelsMap();
  int currHeight = expertDomain->getNumVariables();
  int currLevel = h2lMap[currHeight];

  if (root == 0) {
    root = getHole(currLevel, minNodeSize + 1, true);
    int* nodePtr = getAddress(currLevel, root);
    *nodePtr = minNodeSize;
    memset(nodePtr + 1, 0, *nodePtr * sizeof(int));
  }

  int* currNode = &root;
  while (true) {
    DCASSERT(*currNode != 0);

    // unprimed level
    int currIndex = unpList[currLevel];
    int* ptr = getAddress(currLevel, *currNode);
    // expand currNode if necessary
    if (*ptr <= currIndex) {
      DCASSERT(currIndex + 2 >= minNodeSize);
      int newSize =
        MIN( MAX( currIndex + 1, *ptr * 2 ) , getLevelSize(currLevel) );
      int temp = getHole(currLevel, newSize + 1, true);
      int* tempPtr = getAddress(currLevel, temp);
      *tempPtr = newSize;
      ptr = getAddress(currLevel, *currNode);
      memcpy(tempPtr + 1, ptr + 1, *ptr * sizeof(int));
      memset(tempPtr + 1 + *ptr, 0, (*tempPtr - *ptr) * sizeof(int));
      makeHole(currLevel, *currNode, *ptr + 1);
      *currNode = temp;
      ptr = tempPtr;
    }
    // create a new branch if necessary
    if (ptr[currIndex + 1] == 0) {
      ptr[currIndex + 1] = getHole(-currLevel, minNodeSize + 1, true);
      int* nodePtr = getAddress(-currLevel, ptr[currIndex + 1]);
      *nodePtr = minNodeSize;
      memset(nodePtr + 1, 0, *nodePtr * sizeof(int));
    }
    currNode = ptr + 1 + currIndex;

    DCASSERT(*currNode != 0);
    // primed level
    currIndex = pList[currLevel];
    ptr = getAddress(-currLevel, *currNode);
    // expand currNode if necessary
    if (*ptr <= currIndex) {
      DCASSERT(currIndex +2 >= minNodeSize);
      int newSize =
        MIN( MAX( currIndex + 1, *ptr * 2 ) , getLevelSize(-currLevel) );
      int temp = getHole(-currLevel, newSize + 1, true);
      int* tempPtr = getAddress(-currLevel, temp);
      *tempPtr = newSize;
      ptr = getAddress(-currLevel, *currNode);
      memcpy(tempPtr + 1, ptr + 1, *ptr * sizeof(int));
      memset(tempPtr + 1 + *ptr, 0, (*tempPtr - *ptr) * sizeof(int));
      makeHole(-currLevel, *currNode, *ptr + 1);
      *currNode = temp;
      ptr = tempPtr;
    }
    // deal with terminal case (height == -1, i.e. first level
    // above the terminal nodes).
    if (currHeight == 1) {
      // store the terminals and exit
      ptr[1 + currIndex] = terminalNode;
      break;
    }
    // create a new branch if necessary
    currLevel = h2lMap[--currHeight];
    if (ptr[currIndex + 1] == 0) {
      ptr[currIndex + 1] = getHole(currLevel, minNodeSize + 1, true);
      int* nodePtr = getAddress(currLevel, ptr[currIndex + 1]);
      *nodePtr = minNodeSize;
      memset(nodePtr + 1, 0, *nodePtr * sizeof(int));
    }
    currNode = ptr + 1 + currIndex;
  }
}


template <typename T>
inline
T handleMultipleTerminalValues(T* tList, int begin, int end)
{
  DCASSERT(begin < end);
  T result = tList[begin++];
  while (begin != end) result += tList[begin++];
  return result;
}


template <>
inline
bool handleMultipleTerminalValues(bool* tList, int begin, int end)
{
  DCASSERT(begin < end);
  return true;
}


#endif


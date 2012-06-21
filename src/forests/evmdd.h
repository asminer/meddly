
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
// TODO: evmdd_node_manager

// TODO: inPlaceSortBuild() must be modified to deal with don't care and
//       don't change while building the node instead of deal with them
//       separately (before the call to inPlaceSortBuild()).
//       For this purpose, verify that compute_manager::UNION and PLUS
//       work with nodes that are not at the top-level correctly.


/* 
  TODO: ensure this rule
  All extensions must over-ride either reduceNode() or normalizeAndReduceNode().
  Normalizing is only for edge-valued decision diagrams.
*/
  
#ifndef EVMDD_H
#define EVMDD_H

#include "mt.h"

#define SORT_BUILD

//#define TREE_SORT
#define IN_PLACE_SORT

namespace MEDDLY {
  class evmdd_forest;
  class evp_mdd_int;    // TBD - split into its own file
  class evt_mdd_real;   // TBD - split into its own file
};

// ******************************************************************

const int evplusmddDataHeaderSize = 5;
const int evtimesmddDataHeaderSize = 4;

int binarySearch(const int* a, int sz, int find);

class MEDDLY::evmdd_forest : public mt_forest {
  public:
    evmdd_forest(int dsl, domain *d, range_type t,
        edge_labeling el, const policies &p, int dataHeaderSize);
    ~evmdd_forest();

    virtual void getDefaultEdgeValue(int& n) const { n = INF; }
    virtual void getIdentityEdgeValue(int& n) const { n = 0; }
    virtual void getDefaultEdgeValue(float& n) const { n = NAN; }
    virtual void getIdentityEdgeValue(float& n) const { n = 1; }

    virtual void initEdgeValues(int p) = 0;

  protected:
    using mt_forest::createTempNode;
    virtual int createTempNode(int k, int sz, bool clear = true);

    // Enlarges a temporary node, if new size is greater than old size.
    virtual void resizeNode(int node, int size);

    virtual int doOp(int a, int b) const = 0;
    virtual float doOp(float a, float b) const = 0;

    // Create edge representing vlist[], ending in a terminal value term
    // (or you could think of it as a vector with an edge-value equal to term),
    // and store it in curr.
    template <typename T>
    void createEdgeInternal(const int* vlist, T term, dd_edge &e);

    // Similar to above except that all paths lead to term.
    template <typename T>
    void createEdgeInternal(T term, dd_edge &e);

    // This create a EVMDD from a collection of edges (represented 
    // as vectors).
    template <typename T>
    void createEdgeInternal(const int* const* vlist,
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
    void evaluateInternal(const dd_edge &f, const int* vlist,
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

    template <typename T>
    T handleMultipleTerminalValues(const T* tList, int begin, int end);

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual int reduceNode(int p);

    // Refer to meddly.h:
    // The following will either abort or return an error since they are not
    // applicable to this forest.

    virtual void createEdge(const int* const* vlist, const int* terms,
        int N, dd_edge &e);
    virtual void createEdge(int val, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, int &term)
      const;

    virtual void createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual void createEdge(float val, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, float &term)
      const;

    virtual void createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        const float* terms, int N, dd_edge &e);
    virtual void createEdge(bool val, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual void evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
    virtual void evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, int &term) const;
    virtual void evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, float &term) const;
};


class MEDDLY::evp_mdd_int : public evmdd_forest {
  // EV+MDD data header:
  // { incount, next (unique table), size, ..., cardinality, logical address }
  // Data Header Size: 5

  public:
    evp_mdd_int(int dsl, domain *d, const policies &p);
    ~evp_mdd_int();

    virtual int createTempNode(int k, int sz, bool clear = true);
    virtual int createTempNode(int lh, std::vector<int>& downPointers,
        std::vector<int>& edgeValues);

    virtual void createNode(int lh, std::vector<int>& index,
        std::vector<int>& dptr, std::vector<int>& ev,
        int& result, int& resultEv);

    virtual void createEdge(const int* const* vlist, const int* terms,
        int N, dd_edge &e);
    virtual void createEdge(int val, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, int &term)
      const;

    virtual void getElement(const dd_edge& a, int index, int* e);
    virtual void getElement(int a, int index, int* e);

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


class MEDDLY::evt_mdd_real : public evmdd_forest {
  // EV*MDD data header:
  // { incount, next (unique table), size, ..., logical address }
  // Data Header Size: 4

  public:
    evt_mdd_real(int dsl, domain *d, const policies &p);
    ~evt_mdd_real();

    using evmdd_forest::getDefaultEdgeValue;
    using evmdd_forest::getIdentityEdgeValue;
    virtual void getDefaultEdgeValue(int& n) const {
      static int ev = toInt(NAN);
      n = ev;
    }
    virtual void getIdentityEdgeValue(int& n) const {
      static int ev = toInt(1.0f);
      n = ev;
    }

    virtual int createTempNode(int lh, std::vector<int>& downPointers,
        std::vector<float>& edgeValues);

    virtual void createNode(int lh, std::vector<int>& index,
        std::vector<int>& dptr, std::vector<float>& ev,
        int& result, float& resultEv);

    virtual void createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual void createEdge(float val, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, float &term)
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
void
MEDDLY::evmdd_forest::createEdgeInternal(T term, dd_edge &e)
{
  if (e.getForest() != this) throw error(error::INVALID_OPERATION);

  if (isFullyReduced()) {
    e.set(getTerminalNode(true), term, 0);
  } else {
    // construct the edge bottom-up
    int prev = getTerminalNode(false);
    T prevEv = 0;
    int curr = getTerminalNode(true);
    T currEv = term;
    for (int i=1; i<=getExpertDomain()->getNumVariables(); i++) {
      prev = curr;
      prevEv = currEv;
      createNode(i, -1, prev, prevEv, curr, currEv);
      unlinkNode(prev);
    }
    e.set(curr, currEv, getNodeLevel(curr));
  }
}


template <typename T>
void
MEDDLY::evmdd_forest::createEdgeInternal(const int* const* vlist,
    const T* terms, int N, dd_edge &e)
{
  if (e.getForest() != this) throw error(error::INVALID_OPERATION);
  if (vlist == 0 || terms == 0 || N <= 0) throw error(error::INVALID_VARIABLE);

  // check if the vlist contains valid indexes
  bool specialCasesFound = false;
  for (int i = 0; i < N; i++)
  {
    for (int level = getExpertDomain()->getNumVariables(); level; level--) {
      int bound = vlist[i][level] + 1;
      if (bound >= getLevelSize(level))
        useExpertDomain()->enlargeVariableBound(level, false, bound);
      else if (bound < 0)
        specialCasesFound = true;
    } // for level
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
}


template <typename T>
void
MEDDLY::evmdd_forest::createEdgeInternal(const int* v, T term, dd_edge &e)
{
  // construct the edge bottom-up
  int prev = getTerminalNode(false);
  T prevEv;
  getIdentityEdgeValue(prevEv);
  int curr = getTerminalNode(true);
  T currEv = term;
  for (int i=1; i<=getExpertDomain()->getNumVariables(); i++) {
    prev = curr;
    prevEv = currEv;
    createNode(i, v[i], prev, prevEv, curr, currEv);
    unlinkNode(prev);
  }
  e.set(curr, currEv, getNodeLevel(curr));
}

template <typename T>
void MEDDLY::evmdd_forest::createNode(int k, int index, int dptr, T ev,
    int& res, T& resEv)
{
  MEDDLY_DCASSERT(dptr >= -1 && ev >= 0);
  MEDDLY_DCASSERT(index >= -1);

  if (index == -1) {
    // if all edge values are the same for a evmdd node, it is
    // equivalent to say that they are all 0 (since we subtract the minimum
    // anyway).
    if (isFullyReduced()) {
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
  if (!areSparseNodesEnabled() || (areFullNodesEnabled() && index < 2)) {
    // Build a full node
    res = createTempNode(k, index + 1);
    setDownPtrWoUnlink(res, index, dptr);
    setEdgeValue(res, index, ev);
    normalizeAndReduceNode(res, resEv);
  }
  else {
    MEDDLY_DCASSERT (!areFullNodesEnabled() ||
        (areSparseNodesEnabled() && index >= 2));
    // Build a sparse node
    createSparseNode(k, index, dptr, ev, res, resEv);
    // res, resEv set by createSparseNode(..); do not call normalizeAndReduce
  }
}


template <typename T>
void MEDDLY::evmdd_forest::createSparseNode(int k, int index,
    int dptr, T ev, int& res, T& resEv)
{
  MEDDLY_DCASSERT(k != 0);

  if (isTimeToGc()) { garbageCollect(); }

  MEDDLY_DCASSERT(isValidLevel(k));
  MEDDLY_CHECK_RANGE(0, index, getLevelSize(k));

  // get a location in address[] to store the node
  int p = getFreeNode(k);

#ifdef DEBUG_MDD_H
  printf("%s: k: %d, index: %d, new p: %d\n", __func__, k, index, p);
  fflush(stdout);
#endif

  // fill in the location with p's address info

  const int nodeSize = getDataHeaderSize() + 3;
  address[p].level = k;
  // address[p].offset = levels[k].getHole(nodeSize, true);
  address[p].offset = levels[k].allocNode(-1, p, false);
  address[p].cache_count = 0;

#ifdef DEBUG_MDD_H
  printf("%s: offset: %d\n", __func__, address[p].offset);
  fflush(stdout);
#endif

  int* foo = levels[k].data + address[p].offset;
  // foo[0] = 1;                     // #incoming
  // foo[1] = getTempNodeId();
  // foo[2] = -1;                    // size
  foo[3] = index;                 // index
  foo[4] = sharedCopy(dptr);      // downpointer
  T identEv;
  getIdentityEdgeValue(identEv);
  foo[5] = toInt(identEv);        // this is the only ev, set resEv = ev
  foo[6] = -1;                    // cardinality (-1: not been computed)
  // foo[7] = p;                     // pointer to this node in the address array

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
    MEDDLY_DCASSERT(getCacheCount(p) == 0);
    MEDDLY_DCASSERT(find(p) == p);
    res = p;
  }
  else {
    // duplicate found; discard this node and return the duplicate
    unlinkNode(dptr);
    // code from deleteTempNode(p) adapted to work here
    {
      levels[k].makeHole(getNodeOffset(p), nodeSize);
      freeNode(p);
      if (levels[k].compactLevel) levels[k].compact(address);
    }
    res = sharedCopy(q);
  }
}


template <typename T>
void MEDDLY::evmdd_forest::evaluateInternal(const dd_edge &f,
    const int* vlist, T &term) const
{
  if (f.getForest() != this) throw error(error::INVALID_OPERATION);
  if (vlist == 0) throw error(error::INVALID_VARIABLE);

  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  int node = f.getNode();

  f.getEdgeValue(term);
  while (!isTerminalNode(node)) {
    T ev;
    getDefaultEdgeValue(ev);
    int n = 0;
    if (isFullNode(node)) {
      int index = vlist[getNodeHeight(node)];
      if (index < getFullNodeSize(node)) {
        getFullNodeEdgeValue(node, index, ev);
        n = getFullNodeDownPtr(node, index);
      }
    } else {
      // find the index
      // binary search
      int index = binarySearch(getSparseNodeIndexes(node),
          getSparseNodeSize(node), vlist[getNodeHeight(node)]);
      if (index != -1 ) {
        getSparseNodeEdgeValue(node, index, ev);
        n = getSparseNodeDownPtr(node, index);
      }
    }
    term = (n == 0)? ev: doOp(term, ev);
    node = n;
  }
}


template <typename T>
void MEDDLY::evmdd_forest::sortBuild(int** list, T* tList,
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
  int nextHeight = height - 1;

  if (N == 1) {
    // nothing to sort; just build a node starting at this level
    int index = list[begin][height];
    int n;
    T ev;
    sortBuild(list, tList, nextHeight, begin, end, n, ev);
    createNode(height, index, n, ev, node, edgeValue);
    unlinkNode(n);
    return;
  }

  // do radix sort for this level

  int levelSize = getLevelSize(height);
#if 0
  vector<int> count(levelSize+1, 0);
#else
  std::vector<int> count(1, 0);
#endif

  // Curly braces here to limit the scope of the vectors defined within.
  // Without limiting the scope, memory usage will increase significantly
  // -- especially if there are a lot of variables in the domain.
  {
    std::vector<int*> sortedList(N, (int*)0);
    std::vector<T> sortedtList(N, 0);

    // determine size for count[]
    levelSize = 0;
    for (int i = begin; i < end; i++) {
      int index = list[i][height];
      if (index > levelSize) { levelSize = index; }
    }
    // levelSize refers to the maximum index found so far,
    // add 1 to convert to maximum size.
    levelSize++;
    // an extra space is needed at the end for the radix sort algorithm
    count.resize(levelSize+1, 0);

    // go through list and count the number of entries in each "bucket"
    for (int i = begin; i < end; i++) { count[list[i][height]]++; }

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
    std::vector<int>::iterator last = count.begin() + levelSize;
    for (std::vector<int>::iterator iter = count.begin(); iter != last; iter++)
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
      int index = count[(*listPtr)[height]]--;
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
    std::vector<int*>::iterator sortedListIter = sortedList.begin();
    typename std::vector<T>::iterator sortedtListIter = sortedtList.begin();
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
  std::vector<int> indices;
  std::vector<int> dptrs;
  std::vector<T> edgeValues;
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

  MEDDLY_DCASSERT(dptrs.size() > 0);
  createNode(height, indices, dptrs, edgeValues, node, edgeValue);
}

namespace MEDDLY {

template <typename T>
inline
T evmdd_forest::handleMultipleTerminalValues(const T* tList,
    int begin, int end)
{
  MEDDLY_DCASSERT(begin < end);
  T result = tList[begin++];
  while (begin != end) result += tList[begin++];
  return result;
}


template <>
inline
bool evmdd_forest::handleMultipleTerminalValues(const bool* tList,
    int begin, int end)
{
  MEDDLY_DCASSERT(begin < end);
  return true;
}

} // namespace

#endif


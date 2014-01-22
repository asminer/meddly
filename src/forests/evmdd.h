
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


#ifndef EVMDD_H
#define EVMDD_H

#include <vector>

#define SORT_BUILD

//#define TREE_SORT

namespace MEDDLY {
  class evmdd_forest;
  class evp_mdd_int;    // TBD - split into its own file
  class evt_mdd_real;   // TBD - split into its own file
};

// ******************************************************************

// int binarySearch(const int* a, int sz, int find);

class MEDDLY::evmdd_forest : public expert_forest {
  public:
    evmdd_forest(int dsl, domain *d, range_type t,
        edge_labeling el, const policies &p);
    virtual ~evmdd_forest();

  // ------------------------------------------------------------
  // still to be organized:
    virtual void getDefaultEdgeValue(int& n) const { n = INF; }
    virtual void getIdentityEdgeValue(int& n) const { n = 0; }
    virtual void getDefaultEdgeValue(float& n) const { n = NAN; }
    virtual void getIdentityEdgeValue(float& n) const { n = 1; }


  protected:
    virtual void showTerminal(FILE* s, node_handle tnode) const;
    virtual void writeTerminal(FILE* s, node_handle tnode) const;
    virtual node_handle readTerminal(FILE* s);

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
    void createNode(int k, int index, node_handle dptr, T ev, node_handle& res, T& resEv);

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
        node_handle& node, T& edgeValue);

    // virtual void normalizeAndReduceNode(int& p, int& ev) = 0;
    // virtual void normalizeAndReduceNode(int& p, float& ev) = 0;

    template <typename T>
    T handleMultipleTerminalValues(const T* tList, int begin, int end);

};


class MEDDLY::evp_mdd_int : public evmdd_forest {
  public:
    evp_mdd_int(int dsl, domain *d, const policies &p);
    ~evp_mdd_int();

    virtual void createEdge(int** vlist, int* terms, int N, dd_edge &e);
    virtual void createEdge(int val, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, int &term)
      const;

    virtual void getElement(const dd_edge& a, int index, int* e);

    virtual int doOp(int a, int b) const { return a + b; }
    virtual float doOp(float a, float b) const {
      assert(false); return a + b;
    }

    virtual bool areEdgeValuesEqual(const void* eva, const void* evb) const;

    virtual bool isRedundant(const node_builder &nb) const;
    virtual bool isIdentityEdge(const node_builder &nb, int i) const;

  protected:
    virtual const char* codeChars() const;
    virtual void normalize(node_builder &nb, int& ev) const;
    virtual void showEdgeValue(FILE* s, const void* edge) const;
    virtual void writeEdgeValue(FILE* s, const void* edge) const;
    virtual void readEdgeValue(FILE* s, void* edge);
    virtual void showUnhashedHeader(FILE* s, const void* uh) const;
    virtual void writeUnhashedHeader(FILE* s, const void* uh) const;
    virtual void readUnhashedHeader(FILE* s, node_builder &nb) const;
};


class MEDDLY::evt_mdd_real : public evmdd_forest {
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

    virtual void createEdge(int** vlist, float* terms, int N, dd_edge &e);
    virtual void createEdge(float val, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, float &term)
      const;

    virtual int doOp(int a, int b) const { assert(false); return a * b; }
    virtual float doOp(float a, float b) const { return a * b; }

    virtual bool areEdgeValuesEqual(const void* eva, const void* evb) const;

    virtual bool isRedundant(const node_builder &nb) const;
    virtual bool isIdentityEdge(const node_builder &nb, int i) const;

  protected:
    virtual const char* codeChars() const;
    virtual void normalize(node_builder &nb, float& ev) const;
    virtual void showEdgeValue(FILE* s, const void* edge) const;
    virtual void writeEdgeValue(FILE* s, const void* edge) const;
    virtual void readEdgeValue(FILE* s, void* edge);

  private:
    static inline bool notClose(float a, float b) {
      if (a) {
        double diff = a-b;
        return ABS(diff/a) > 1e-6;
      } else {
        return ABS(b) > 1e-10;
      }
    }
    
};

























// ------------------------ Inline methods -----------------------------------

template <typename T>
void
MEDDLY::evmdd_forest::createEdgeInternal(T term, dd_edge &e)
{
  if (e.getForest() != this) throw error(error::INVALID_OPERATION);

  if (isFullyReduced()) {
    e.set(getTerminalNode(true), term);
  } else {
    // construct the edge bottom-up
    node_handle prev = getTerminalNode(false);
    T prevEv = 0;
    node_handle curr = getTerminalNode(true);
    T currEv = term;
    for (int i=1; i<=getExpertDomain()->getNumVariables(); i++) {
      prev = curr;
      prevEv = currEv;
      createNode(i, -1, prev, prevEv, curr, currEv);
      // unlinkNode(prev);
    }
    e.set(curr, currEv);
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
    node_handle node;
    T ev;
    sortBuild(list, tList, getDomain()->getNumVariables(), 0, N, node, ev);
    e.set(node, ev);
  }
}


template <typename T>
void
MEDDLY::evmdd_forest::createEdgeInternal(const int* v, T term, dd_edge &e)
{
  // construct the edge bottom-up
  node_handle prev = getTerminalNode(false);
  T prevEv;
  getIdentityEdgeValue(prevEv);
  node_handle curr = getTerminalNode(true);
  T currEv = term;
  for (int i=1; i<=getExpertDomain()->getNumVariables(); i++) {
    prev = curr;
    prevEv = currEv;
    createNode(i, v[i], prev, prevEv, curr, currEv);
  }
  e.set(curr, currEv);
}

template <typename T>
void MEDDLY::evmdd_forest::createNode(int k, int index, node_handle dptr, T ev,
    node_handle& res, T& resEv)
{
  MEDDLY_DCASSERT(dptr >= -1 && ev >= 0);
  MEDDLY_DCASSERT(index >= -1);

  if (index == -1) {
    // if all edge values are the same for a evmdd node, it is
    // equivalent to say that they are all 0 (since we subtract the minimum
    // anyway).
    if (isFullyReduced()) {
      res = dptr;
      resEv = ev;
      return;
    }
    int sz = getLevelSize(k);
    node_builder& nb = useNodeBuilder(k, sz);
    for (int i=0; i<sz; i++) {
      nb.d(i) = linkNode(dptr);
      nb.setEdge(i, ev);
    }
    unlinkNode(dptr);
    createReducedNode(-1, nb, ev, res);
    return;
  }

  // a single downpointer points to dptr
  node_builder& nb = useSparseBuilder(k, 1);
  nb.d(0) = dptr;
  nb.i(0) = index;
  nb.setEdge(0, ev);
  createReducedNode(-1, nb, resEv, res);
}


template <typename T>
void MEDDLY::evmdd_forest::evaluateInternal(const dd_edge &f,
    const int* vlist, T &term) const
{
  if (f.getForest() != this) throw error(error::INVALID_OPERATION);
  if (vlist == 0) throw error(error::INVALID_VARIABLE);

  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  node_handle node = f.getNode();

  f.getEdgeValue(term);
  while (!isTerminalNode(node)) {
    T ev;
    getDownPtr(node, vlist[getNodeHeight(node)], ev, node);
    term = (node) ? doOp(term, ev) : ev;
  }
}


template <typename T>
void MEDDLY::evmdd_forest::sortBuild(int** list, T* tList,
    int height, int begin, int end, node_handle& node, T& edgeValue)
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
    node_handle n;
    T ev;
    sortBuild(list, tList, nextHeight, begin, end, n, ev);
    createNode(height, index, n, ev, node, edgeValue);
    // unlinkNode(n);
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

  node_builder& nb = useSparseBuilder(height, levelSize);
  int z = 0;
  for (int i=0; i<levelSize; i++) {
    if (count[i+1] > count[i]) {
      node_handle n;
      T ev;
      sortBuild(list, tList, nextHeight, begin + count[i] + 1,
          begin + count[i+1] + 1, n, ev);
      nb.i(z) = i;
      nb.d(z) = n;
      nb.setEdge(z, ev);
      z++;
    }
  }
  nb.shrinkSparse(z);
  createReducedNode(-1, nb, edgeValue, node);
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


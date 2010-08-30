
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
// TODO: mtmdd_node_manager

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
  
#ifndef MTMDD_H
#define MTMDD_H

#include "mdds.h"

#define IN_PLACE_SORT

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
    T handleMultipleTerminalValues(const T* tList, int begin, int end);

#ifndef IN_PLACE_SORT
    template <typename T>
    int sortBuild(int** list, T* tList, int height, int begin, int end);
#else
    template <typename T>
    int inPlaceSort(int level, int begin, int end);
    template <typename T>
    int inPlaceSortBuild(int height, int begin, int end);
#endif

    // Methods and data for batch addition via sorting
    template <typename T> void copyLists(const int* const* vlist,
        const T* terms, int nElements);

    void expandCountAndSlotArrays(int size);

    int** list;
    int*  termList;
    int   listSize;

    int* count;
    int* slot;
    int countSize;
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



// ------------------------ Inline methods -----------------------------------


template <typename T>
inline
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
    DCASSERT(N > 0);

    // copy elements into internal volatile storage
    copyLists(vlist, terms, N);

    // call sort-based procedure for building the DD
#ifdef IN_PLACE_SORT
    int result = inPlaceSortBuild<T>(expertDomain->getNumVariables(), 0, N);
#else
    int result = sortBuild(list, (T*)(terms == 0? 0: termList),
        expertDomain->getNumVariables(), 0, N);
#endif

    e.set(result, 0, getNodeLevel(result));
  }

  return forest::SUCCESS;
}


template <typename T>
inline
void mtmdd_node_manager::copyLists(const int* const* vlist,
    const T* terms, int nElements)
{
  if (listSize < nElements) {
    list = (int**) realloc(list, sizeof(int*) * nElements);
    assert(list);
    if (terms) {
      termList = (int*) realloc(termList, sizeof(T) * nElements);
      assert(termList);
    }
    listSize = nElements;
  }

  memcpy(list, vlist, nElements * sizeof(int*));
  if (terms) {
    T* tempTList = (T*)termList;
    // Not doing memcpy in case T is a class
    for (int i = 0; i < nElements; i++) { tempTList[i] = terms[i]; }
  }
}



#ifndef IN_PLACE_SORT

template <typename T>
inline
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


#else


template<typename T>
inline
int mtmdd_node_manager::inPlaceSort(int level, int begin, int end)
{
  // Determine range of values
  int min = list[begin][level];
  int max = min;
  for (int i = begin + 1; i < end; ++i) {
    max = MAX(max, list[i][level]);
    min = MIN(min, list[i][level]);
  }

  // Prepare arrays (expand them as necessary and clear them as necessary).
  // TODO: move the array definitions and initialization out of here.
  expandCountAndSlotArrays(max + 1 - min);

#ifdef DEVELOPMENT_CODE
  for (int i = 0; i < countSize; i++) { assert(0 == count[i]); }
#endif

  // c and s reduce the number of subtractions in indexes
  int* c = count - min;
  int* s = slot - min;

  // Count the number of entries for each value
  for (int i = begin; i < end; i++) {
    c[list[i][level]]++;
  }

  // Determine the initial slot positions
  s[min] = begin;
  for (int i = min + 1; i <= max; ++i) {
    s[i] = s[i-1] + c[i-1];
  }

  // We have the correct bucket sizes, now move items into
  // appropriate buckets.

  T* terms = (T*)termList;
  
  for (int i = min; i < max; ++i) {
    // Move elements in bucket i to the correct slots.
    // Repeat this until all the elements in bucket i belong in bucket i.
    while (c[i] > 0) {
      // Find appropriate slot for list[s[i]]
      int* elem = list[s[i]];
      int elemIndex = elem[level];
      if (i == elemIndex) {
        // Already in the correct slot
        --c[i];
        ++s[i];
      }
      else {
        // Move elem to correct slot
        DCASSERT(elemIndex > i);
        while (c[elemIndex] > 0 && elemIndex == list[s[elemIndex]][level]) {
          // These elements are already in the correct slots; advance pointers.
          --c[elemIndex];
          ++s[elemIndex];
        }
        // At correct slot for elem
        DCASSERT(c[elemIndex] > 0);
        CHECK_RANGE(begin, s[elemIndex], end);
        SWAP(list[s[i]], list[s[elemIndex]]);
        if (terms) { SWAP(terms[s[i]], terms[s[elemIndex]]); }
        --c[elemIndex];
        ++s[elemIndex];
        // list[s[elemIndex]] now contains the correct element.
        // Also, list[s[i]] now contains an unknown and this 
        // will be handled in the next iteration.
        // Note that we do not advance c[i] and s[i].
      }
    }
    // Bucket i now contains only elements that belong in it.
  }

  c[max] = 0;

#ifdef DEVELOPMENT_CODE
  // Check if all buckets have been dealt with
  for (int i = min; i <= max; i++) { assert(0 == c[i]); }
#endif

#ifdef DEVELOPMENT_CODE
  // Check if sorted
  for (int i = begin + 1; i < end; i++) {
    assert(list[i-1][level] <= list[i][level]);
  }
#endif

  // max represents the largest index
  return max;
}


template<>
inline
int mtmdd_node_manager::inPlaceSort<bool>(int level, int begin, int end)
{
  // Determine range of values
  int min = list[begin][level];
  int max = min;
  for (int i = begin + 1; i < end; ++i) {
    max = MAX(max, list[i][level]);
    min = MIN(min, list[i][level]);
  }

  // Prepare arrays (expand them as necessary and clear them as necessary).
  expandCountAndSlotArrays(max + 1 - min);

#ifdef DEVELOPMENT_CODE
  for (int i = 0; i < countSize; i++) { assert(0 == count[i]); }
#endif

  // c and s reduce the number of subtractions in indexes
  int* c = count - min;
  int* s = slot - min;

  // Count the number of entries for each value
  for (int i = begin; i < end; i++) {
    ++c[list[i][level]];
  }

  // Determine the initial slot positions
  s[min] = begin;
  for (int i = min + 1; i <= max; ++i) {
    s[i] = s[i-1] + c[i-1];
  }

  // We have the correct bucket sizes, now move items into
  // appropriate buckets.

  for (int i = min; i < max; ++i) {
    // Move elements in bucket i to the correct slots.
    // Repeat this until all the elements in bucket i belong in bucket i.
    while (c[i] > 0) {
      // Find appropriate slot for list[s[i]]
      int* elem = list[s[i]];
      int elemIndex = elem[level];
      if (i == elemIndex) {
        // Already in the correct slot
        --c[i];
        ++s[i];
      }
      else {
        // Move elem to correct slot
        DCASSERT(elemIndex > i);
        while (c[elemIndex] > 0 && elemIndex == list[s[elemIndex]][level]) {
          // These elements are already in the correct slots; advance pointers.
          --c[elemIndex];
          ++s[elemIndex];
        }
        // At correct slot for elem
        DCASSERT(c[elemIndex] > 0);
        CHECK_RANGE(begin, s[elemIndex], end);
        SWAP(list[s[i]], list[s[elemIndex]]);
        --c[elemIndex];
        ++s[elemIndex];
        // list[s[elemIndex]] now contains the correct element.
        // Also, list[s[i]] now contains an unknown and this 
        // will be handled in the next iteration.
        // Note that we do not advance c[i] and s[i].
      }
    }
    // Bucket i now contains only elements that belong in it.
  }

  c[max] = 0;

#ifdef DEVELOPMENT_CODE
  // Check if all buckets have been dealt with
  for (int i = min; i <= max; i++) { assert(0 == c[i]); }
#endif

#ifdef DEVELOPMENT_CODE
  // Check if sorted
  for (int i = begin + 1; i < end; i++) {
    assert(list[i-1][level] <= list[i][level]);
  }
#endif

  // max represents the largest index
  return max;
}


template <typename T>
inline
int mtmdd_node_manager::inPlaceSortBuild(int height, int begin, int end)
{
  // [begin, end)

  // terminal condition
  if (height == 0) {
    return getTerminalNode(
        handleMultipleTerminalValues((T*)termList, begin, end));
  }

  int nextHeight = height - 1;
  int level = expertDomain->getVariableWithHeight(height);
  DCASSERT(level > 0);

  if (begin + 1 == end) {
    // nothing to sort; just build a node starting at this level
    int n = inPlaceSortBuild<T>(nextHeight, begin, end);
    int index = list[begin][level];
    int result = createNode(level, index, n);
    unlinkNode(n);
    return result;
  }

  // Sort elements at this level
  int nodeSize = 1 + inPlaceSort<T>(level, begin, end);

  // build node
  int result = createTempNode(level, nodeSize, true);
  int* ptr = getFullNodeDownPtrs(result);
  for (int i = begin; i < end; )
  {
    int index = list[i][level];
    int start = i++;
    // skip the elements with the same index at this level
    for ( ; i < end && list[i][level] == index; ++i);
    ptr[index] = inPlaceSortBuild<T>(nextHeight, start, i);
  }

  return reduceNode(result);
}

#endif


template <typename T>
inline
T mtmdd_node_manager::handleMultipleTerminalValues(const T* tList,
    int begin, int end)
{
  DCASSERT(begin < end);
  T result = tList[begin++];
  while (begin != end) result += tList[begin++];
  return result;
}


template <>
inline
bool mtmdd_node_manager::handleMultipleTerminalValues(const bool* tList,
    int begin, int end)
{
  DCASSERT(begin < end);
  return true;
}


#endif


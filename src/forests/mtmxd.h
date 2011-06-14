
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
// TODO: mxd_node_manager
// TODO: mtmxd_node_manager (??)

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
  
#ifndef MTMXD_H
#define MTMXD_H

#include "mdds.h"

#define IN_PLACE_SORT

const int mtmxdDataHeaderSize = 4;

class mtmxd_node_manager : public node_manager {
  // TODO: mtmxds can only be forest::IDENTITY_REDUCED
  // MTMDD data header:
  // { incount, next (unique table), size, ..., logical address}
  // Data Header Size: 4

  public:

    mtmxd_node_manager(domain *d, forest::range_type t);
    ~mtmxd_node_manager();

    // Refer to meddly.h
    void createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge& e);
    void createEdge(const int* const* vlist, const int* const* vplist,
        const float* terms, int N, dd_edge& e);

    void createEdge(int val, dd_edge &e);
    void createEdge(float val, dd_edge &e);

    void evaluate(const dd_edge& f, const int* vlist, const int* vplist,
        int &term) const;
    void evaluate(const dd_edge& f, const int* vlist, const int* vplist,
        float &term) const;

    virtual void findFirstElement(const dd_edge& f, int* vlist, int* vplist)
      const;

    // Enlarges a temporary node, if new size is greater than old size.
    virtual void resizeNode(int node, int size);

    virtual int reduceNode(int p);

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev);

    // Refer to meddly.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.

    virtual void createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual void createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);

    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e);
    virtual void createEdge(bool val, dd_edge &e);

    virtual void evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual void evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    virtual void evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
    virtual void evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;

  protected:

    // Used by derived classes for initialization
    mtmxd_node_manager(domain *d, bool relation, forest::range_type t,
        forest::edge_labeling e, forest::reduction_rule r,
        forest::node_storage s, forest::node_deletion_policy dp);

    // This create a MTMXD from a collection of edges (represented 
    // as vectors vlist and vplist).
    template <typename T>
      void createEdgeInternal(const int* const* vlist,
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
      T handleMultipleTerminalValues(const T* tList, int begin, int end);

#ifndef IN_PLACE_SORT
    template <typename T>
      int sortBuild(int** unpList, int** pList, T* tList,
          int height, int begin, int end);

    template <typename T>
      int sort(int** list, int** otherList, T* tList,
          int absLevel, int begin, int end);

    // Persistant storage for sort()
    vector<int*> sortedList;
    vector<int*> sortedOtherList;
    vector<T> sortedtList;
    vector<int> countVector;

#else
    template <typename T>
      int inPlaceSort(int level, bool isPrime, int begin, int end);
    template <typename T>
      int inPlaceSortBuild(int height, bool isPrime, int begin, int end);
#endif

    // Methods and data for batch addition via sorting
    template <typename T> void copyLists(const int* const* vlist,
        const int* const* vplist, const T* terms, int nElements);

    void expandCountAndSlotArrays(int size);

    int** unpList;
    int** pList;
    void* tList;
    int   listSize;

    int* count;
    int* slot;
    int countSize;
};



class mxd_node_manager : public mtmxd_node_manager {
  // TODO: mxds can only be forest::IDENTITY_REDUCED
  public:

    mxd_node_manager(domain *d);
    ~mxd_node_manager();

    using mtmxd_node_manager::createEdge;
    using mtmxd_node_manager::evaluate;

    virtual bool accumulate(int& tempNode, int* element, int* pelement);
    virtual void accumulate(int& a, int b);

    virtual int accumulate(int tempNode, bool cBM,
        int* element, int* pelement, int level);
    virtual int accumulateSkippedLevel(int tempNode,
        int* element, int* pelement, int level);

    virtual int accumulateMxd(int a, int b, bool cBM);
    virtual int accumulateMxdPrime(int a, int b, bool cBM);
    virtual int addPrimeReducedNodes(int a, int b);
    virtual int accumulateExpandA(int a, int b, bool cBM);
    virtual int buildQRIdentityNode(int node, int level);

    virtual void accumulateMxdHelper(int& a, int b, bool cBM,
        bool needsToMakeACopy,
        int (mxd_node_manager::*function)(int, int, bool));

    // Refer to meddly.h
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge& e);
    virtual void createEdge(bool val, dd_edge &e);
    virtual void evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;

    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        const float* terms, int N, dd_edge &e);
    virtual void createEdge(int val, dd_edge &e);
    virtual void createEdge(float val, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist,
        const int* vplist, int &term) const;
    virtual void evaluate(const dd_edge &f, const int* vlist,
        const int* vplist, float &term) const;
};



// ------------------------ Inline methods -----------------------------------

template <typename T>
inline
void mtmxd_node_manager::copyLists(const int* const* vlist,
    const int* const* vplist, const T* terms, int nElements)
{
  if (listSize < nElements) {
    unpList = (int**) realloc(unpList, sizeof(int*) * nElements);
    if (NULL == unpList) outOfMemory();
    pList = (int**) realloc(pList, sizeof(int*) * nElements);
    if (NULL == pList) outOfMemory();
    tList = (void*) realloc(tList, sizeof(T) * nElements);
    if (NULL == tList) outOfMemory();
    listSize = nElements;
  }

  memcpy(unpList, vlist, nElements * sizeof(int*));
  memcpy(pList, vplist, nElements * sizeof(int*));
  if (terms != 0) {
    T* tempTList = (T*)tList;
    for (int i = 0; i < nElements; i++) { tempTList[i] = terms[i]; }
  }
}


template <typename T>
void
mtmxd_node_manager::createEdgeInternal(const int* const* vlist,
    const int* const* vplist, const T* terms, int N, dd_edge &e)
{
  // check if the vlist contains valid indexes
  bool specialCasesFound = false;
  for (int i = 0; i < N; i++)
  {
    for (int level = expertDomain->getNumVariables(); level; level--) {
      // unprimed level
      int bound = vlist[i][level] + 1;
      if (bound >= expertDomain->getVariableBound(level, false))
        expertDomain->enlargeVariableBound(level, false, bound);
      else if (bound < 0)
        specialCasesFound = true;
      // primed level
      bound = vplist[i][level] + 1;
      if (bound >= expertDomain->getVariableBound(level, true))
        expertDomain->enlargeVariableBound(level, true, bound);
      else if (bound < 0)
        specialCasesFound = true;
    } // for level
  } // for i

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
    // build using sort-based procedure
    DCASSERT(N > 0);

    // copy elements into internal volatile storage
    copyLists(vlist, vplist, terms, N);

    // call sort-based procedure for building the DD
#ifdef IN_PLACE_SORT
    int result = inPlaceSortBuild<T>(expertDomain->getNumVariables(),
        false, 0, N);
#else
    int result = sortBuild(unpList, pList, (T*)(terms == 0? 0: tList),
        expertDomain->getNumVariables(), 0, N);
#endif

    e.set(result, 0, getNodeLevel(result));
  }
}


#ifndef IN_PLACE_SORT

template <typename T>
int mtmxd_node_manager::sort(int** list, int** otherList, T* tList,
    int absLevel, int begin, int end)
{
  DCASSERT(tList != 0);

  vector<int>& count = countVector;

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

  vector<int>& count = countVector;

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

#endif


template<typename T>
inline
int mtmxd_node_manager::inPlaceSort(int level, bool isPrime,
    int begin, int end)
{
  // plist, unpList, tList
  // list, otherList, tList

  int **list = isPrime? pList: unpList;
  int** otherList = isPrime? unpList: pList;

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
    c[list[i][level]]++;
  }

  // Determine the initial slot positions
  s[min] = begin;
  for (int i = min + 1; i <= max; ++i) {
    s[i] = s[i-1] + c[i-1];
  }

  // We have the correct bucket sizes, now move items into
  // appropriate buckets.

  DCASSERT(tList);
  T* terms = (T*)tList;
  DCASSERT(terms);
  
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
        SWAP(otherList[s[i]], otherList[s[elemIndex]]);
        SWAP(terms[s[i]], terms[s[elemIndex]]);
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
int mtmxd_node_manager::inPlaceSort<bool>(int level, bool isPrime,
    int begin, int end)
{
  // plist, unpList, tList
  // list, otherList, tList

  int **list = isPrime? pList: unpList;
  int** otherList = isPrime? unpList: pList;

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
    c[list[i][level]]++;
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
        SWAP(otherList[s[i]], otherList[s[elemIndex]]);
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
int mtmxd_node_manager::inPlaceSortBuild(int height, bool isPrime,
    int begin, int end)
{
  // [begin, end)

  // terminal condition
  if (height == 0) {
    return getTerminalNode(
        handleMultipleTerminalValues((T*)tList, begin, end));
  }

  if (begin + 1 == end) {
    return createNode(unpList[begin], pList[begin],
        getTerminalNode(handleMultipleTerminalValues((T*)tList, begin, end)),
        height, isPrime);
  }

  // Sort elements at this level
  int nodeSize = 1 + inPlaceSort<T>(height, isPrime, begin, end);

  int nextHeight = isPrime? height - 1: height;
  bool nextIsPrime = !isPrime;
  int** list = isPrime? pList: unpList;

  // build node
  int result = createTempNode(isPrime? -height: height, nodeSize, true);
  int* ptr = getFullNodeDownPtrs(result);
  for (int i = begin; i < end; )
  {
    int index = list[i][height];
    int start = i++;
    // skip the elements with the same index at this level
    for ( ; i < end && list[i][height] == index; ++i);
    ptr[index] = inPlaceSortBuild<T>(nextHeight, nextIsPrime, start, i);
  }

  return reduceNode(result);
}


template <typename T>
inline
T mtmxd_node_manager::handleMultipleTerminalValues(const T* tList,
    int begin, int end)
{
  DCASSERT(begin < end);
  T result = tList[begin++];
  while (begin != end) result += tList[begin++];
  return result;
}


template <>
inline
bool mtmxd_node_manager::handleMultipleTerminalValues(const bool* tList,
    int begin, int end)
{
  DCASSERT(begin < end);
  return true;
}


#endif


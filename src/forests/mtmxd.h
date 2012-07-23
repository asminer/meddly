
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

#include "mt.h"

namespace MEDDLY {
  class mtmxd_forest;
};

// ******************************************************************

class MEDDLY::mtmxd_forest : public mt_forest {
  // TODO: mtmxds can only be forest::IDENTITY_REDUCED
  // MTMDD data header:
  // { incount, next (unique table), size, ..., logical address}
  // Data Header Size: 4

  public:
    ~mtmxd_forest();

  protected:

    // Used by derived classes for initialization
    mtmxd_forest(int dsl, domain *d, bool relation, range_type t,
        edge_labeling e, const policies &p);

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
    int createEdgeTo(int termNode);

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

    template <typename T>
      int inPlaceSort(int k, int begin, int end);
    template <typename T>
      int inPlaceSortBuild(int in, int k, int begin, int end);

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

// ------------------------ Inline methods -----------------------------------

template <typename T>
inline
void MEDDLY::mtmxd_forest::copyLists(const int* const* vlist,
    const int* const* vplist, const T* terms, int nElements)
{
  if (listSize < nElements) {
    unpList = (int**) realloc(unpList, sizeof(int*) * nElements);
    if (NULL == unpList) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);
    pList = (int**) realloc(pList, sizeof(int*) * nElements);
    if (NULL == pList) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);
    tList = (void*) realloc(tList, sizeof(T) * nElements);
    if (NULL == tList) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);
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
MEDDLY::mtmxd_forest::createEdgeInternal(const int* const* vlist,
    const int* const* vplist, const T* terms, int N, dd_edge &e)
{
  // check if the vlist contains valid indexes
  bool specialCasesFound = false;
  for (int i = 0; i < N; i++)
  {
    for (int level = getExpertDomain()->getNumVariables(); level; level--) {
      // unprimed level
      int bound = vlist[i][level] + 1;
      if (bound >= getExpertDomain()->getVariableBound(level, false))
        useExpertDomain()->enlargeVariableBound(level, false, bound);
      else if (bound < 0)
        specialCasesFound = true;
      // primed level
      bound = vplist[i][level] + 1;
      if (bound >= getExpertDomain()->getVariableBound(level, true))
        useExpertDomain()->enlargeVariableBound(level, true, bound);
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
    MEDDLY_DCASSERT(N > 0);

    // copy elements into internal volatile storage
    copyLists(vlist, vplist, terms, N);

    // call sort-based procedure for building the DD
    int result = inPlaceSortBuild<T>(-1, getExpertDomain()->getNumVariables(), 0, N);

    e.set(result, 0, getNodeLevel(result));
  }
}



namespace MEDDLY {

template<typename T>
inline
int mtmxd_forest::inPlaceSort(int k, int begin, int end)
{
  // plist, unpList, tList
  // list, otherList, tList

  int **list = (k<0)? pList: unpList;
  int** otherList = (k<0)? unpList: pList;
  int level = ABS(k);

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

  MEDDLY_DCASSERT(tList);
  T* terms = (T*)tList;
  MEDDLY_DCASSERT(terms);
  
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
        MEDDLY_DCASSERT(elemIndex > i);
        while (c[elemIndex] > 0 && elemIndex == list[s[elemIndex]][level]) {
          // These elements are already in the correct slots; advance pointers.
          --c[elemIndex];
          ++s[elemIndex];
        }
        // At correct slot for elem
        MEDDLY_DCASSERT(c[elemIndex] > 0);
        MEDDLY_CHECK_RANGE(begin, s[elemIndex], end);
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
int mtmxd_forest::inPlaceSort<bool>(int k, int begin, int end)
{
  // plist, unpList, tList
  // list, otherList, tList

  int **list = (k<0)? pList: unpList;
  int** otherList = (k<0)? unpList: pList;
  int level = ABS(k);

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
        MEDDLY_DCASSERT(elemIndex > i);
        while (c[elemIndex] > 0 && elemIndex == list[s[elemIndex]][level]) {
          // These elements are already in the correct slots; advance pointers.
          --c[elemIndex];
          ++s[elemIndex];
        }
        // At correct slot for elem
        MEDDLY_DCASSERT(c[elemIndex] > 0);
        MEDDLY_CHECK_RANGE(begin, s[elemIndex], end);
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
int mtmxd_forest::inPlaceSortBuild(int in, int k, int begin, int end)
{
  // [begin, end)

  // terminal condition
  if (0 == k) {
    return getTerminalNode(
        handleMultipleTerminalValues((T*)tList, begin, end));
  }

  if (begin + 1 == end) {
    return createNode(unpList[begin], pList[begin],
        getTerminalNode(handleMultipleTerminalValues((T*)tList, begin, end)),
        ABS(k), (k<0));
  }

  // Sort elements at this level
  int nodeSize = 1 + inPlaceSort<T>(k, begin, end);

  int nextK = (k>0) ? (-k) : (-k-1);
  int** list = (k<0)? pList: unpList;

  // build node
  node_builder& nb = useSparseBuilder(k, nodeSize);
  int height = ABS(k);
  int z = 0;
  for (int i = begin; i < end; )
  {
    int index = list[i][height];
    int start = i++;
    // skip the elements with the same index at this level
    for ( ; i < end && list[i][height] == index; ++i);
    nb.i(z) = index;
    nb.d(z) = inPlaceSortBuild<T>(index, nextK, start, i);
    z++;
  }
  nb.shrinkSparse(z);
  return createReducedNode(in, nb);
}


template <typename T>
inline
T mtmxd_forest::handleMultipleTerminalValues(const T* tList,
    int begin, int end)
{
  MEDDLY_DCASSERT(begin < end);
  T result = tList[begin++];
  while (begin != end) result += tList[begin++];
  return result;
}


template <>
inline
bool mtmxd_forest::handleMultipleTerminalValues(const bool* tList,
    int begin, int end)
{
  MEDDLY_DCASSERT(begin < end);
  return true;
}

} // namespace

#endif


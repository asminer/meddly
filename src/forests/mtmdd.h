
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

#include "mt.h"

#define IN_PLACE_SORT

namespace MEDDLY {
  class mtmdd_forest;
};

// ******************************************************************

class MEDDLY::mtmdd_forest : public mt_forest {

  protected:

    // Used by derived classes for initialization
    mtmdd_forest(int dsl, domain *d, bool relation, range_type t,
        edge_labeling e, const policies &p);

    ~mtmdd_forest();

  protected:
    // Creates an edge representing the terminal node given by
    // terminalNode.
    // Note: terminalNode is usually a terminal value converted into
    // the equivalent node representation. This makes this method useful
    // for createEdge(int, e) or (float, e) or (bool, e) as long as the
    // terminal value can be converted into an equivalent node.
    inline void createEdgeHelper(node_handle terminalNode, dd_edge& e) {
        MEDDLY_DCASSERT(isTerminalNode(terminalNode));

        if (isFullyReduced() || terminalNode == 0) {
          e.set(terminalNode, 0);
          return;
        }

        // construct the edge bottom-up
        node_handle result = terminalNode;
        for (int i=1; i<=getExpertDomain()->getNumVariables(); i++) {
          insertRedundantNode(i, result);
        }
        e.set(result, 0);
    }

    // Get the terminal node at the bottom of the edge with root n
    // and vlist representing the indexes for the levels.
    // Used by evaluate()
    inline node_handle getTerminalNodeForEdge(int n, const int* vlist) const {
        // assumption: vlist does not contain any special values (-1, -2, etc).
        // vlist contains a single element.
        while (!isTerminalNode(n)) {
          n = getDownPtr(n, vlist[getNodeHeight(n)]);
        }
        return n;
    }

    // This create a MTMDD from a collection of edges (represented 
    // as vectors).
    template <typename T>
      void createEdgeInternal(const int* const* vlist,
          const T* terms, int N, dd_edge &e);

  protected: // still to be organized

    // Create edge representing f(vlist[]) = term and store it in e
    void createEdgeTo(const int* vlist, node_handle term, dd_edge& e);

    // Create a node, at level k, whose ith index points to dptr.
    // If i is -1, all indices of the node will point to dptr.
    node_handle createNode(int k, int i, node_handle dptr);

    template <typename T>
    T handleMultipleTerminalValues(const T* tList, int begin, int end);

    template <typename T>
    node_handle inPlaceSort(int level, int begin, int end);
    template <typename T>
    node_handle inPlaceSortBuild(int height, int begin, int end);

    // Methods and data for batch addition via sorting
    template <typename T> void copyLists(const int* const* vlist,
        const T* terms, int nElements);

    void expandCountAndSlotArrays(int size);

  private:
    int** list;
    node_handle*  termList;
    int   listSize;

    int* count;
    int* slot;
    int countSize;
};


// ------------------------ Inline methods -----------------------------------


template <typename T>
inline
void
MEDDLY::mtmdd_forest::createEdgeInternal(const int* const* vlist,
    const T* terms, int N, dd_edge &e)
{
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
  } // for i

  if (N == 1 || specialCasesFound) {
    // build using "standard" procedure
    if (terms == 0) {
      node_handle trueNode = getTerminalNode(true);
      createEdgeTo(vlist[0], trueNode, e);
      if (N > 1) {
        dd_edge curr(this);
        for (int i=1; i<N; i++) {
          createEdgeTo(vlist[i], trueNode, curr);
          e += curr;
        }
      }
    }
    else {
      createEdgeTo(vlist[0], getTerminalNode(terms[0]), e);
      if (N > 1) {
        dd_edge curr(this);
        for (int i=1; i<N; i++) {
          createEdgeTo(vlist[i], getTerminalNode(terms[i]), curr);
          e += curr;
        }
      }
    }
  }
  else {
    // build using sort-based procedure
    MEDDLY_DCASSERT(N > 0);

    // copy elements into internal volatile storage
    copyLists(vlist, terms, N);

    // call sort-based procedure for building the DD
    node_handle result = inPlaceSortBuild<T>(getExpertDomain()->getNumVariables(), 0, N);

    e.set(result, 0);
  }
}


template <typename T>
inline
void MEDDLY::mtmdd_forest::copyLists(const int* const* vlist,
    const T* terms, int nElements)
{
  if (listSize < nElements) {
    list = (int**) realloc(list, sizeof(void*) * nElements);
    if (NULL == list) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);
    if (terms) {
      termList = (node_handle*) realloc(termList, sizeof(node_handle) * nElements);
      if (NULL == termList) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);
    }
    listSize = nElements;
  }

  memcpy(list, vlist, nElements * sizeof(void*));
  if (terms) {
    T* tempTList = (T*)termList;
    // Not doing memcpy in case T is a class
    for (int i = 0; i < nElements; i++) { tempTList[i] = terms[i]; }
  }
}




namespace MEDDLY {

template<typename T>
inline
node_handle mtmdd_forest::inPlaceSort(int level, int begin, int end)
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
    c[list[i][level]]++;
  }

  // Determine the initial slot positions
  s[min] = begin;
  for (int i = min + 1; i <= max; ++i) {
    s[i] = s[i-1] + c[i-1];
  }

  // We have the correct bucket sizes, now move items into
  // appropriate buckets.

  MEDDLY_DCASSERT(termList);
  T* terms = (T*)termList;
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
node_handle mtmdd_forest::inPlaceSort<bool>(int level, int begin, int end)
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
MEDDLY::node_handle mtmdd_forest::inPlaceSortBuild(int height, int begin, int end)
{
  // [begin, end)

  // terminal condition
  if (height == 0) {
    return getTerminalNode(
        handleMultipleTerminalValues((T*)termList, begin, end));
  }

  int nextHeight = height - 1;

  if (begin + 1 == end) {
    // nothing to sort; just build a node starting at this level
    node_handle n = inPlaceSortBuild<T>(nextHeight, begin, end);
    node_handle index = list[begin][height];
    node_handle result = createNode(height, index, n);
    return result;
  }

  // Sort elements at this level
  int nodeSize = 1 + inPlaceSort<T>(height, begin, end);

  // build node
  node_builder& nb = useSparseBuilder(height, nodeSize);
  int z = 0;
  for (int i = begin; i < end; ) {
    int index = list[i][height];
    int start = i++;
    // skip the elements with the same index at this level
    for ( ; i < end && list[i][height] == index; ++i);
    // set next downward pointer
    nb.i(z) = index;
    nb.d(z) = inPlaceSortBuild<T>(nextHeight, start, i);
    z++;
  }
  nb.shrinkSparse(z);
  return createReducedNode(-1, nb);
}


template <typename T>
inline
T mtmdd_forest::handleMultipleTerminalValues(const T* tList,
    int begin, int end)
{
  MEDDLY_DCASSERT(begin < end);
  T result = tList[begin++];
  while (begin != end) result += tList[begin++];
  return result;
}


template <>
inline
bool mtmdd_forest::handleMultipleTerminalValues(const bool* tList,
    int begin, int end)
{
  MEDDLY_DCASSERT(begin < end);
  return true;
}

} // namespace

#endif


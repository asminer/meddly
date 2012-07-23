
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


#include "mtmdd.h"
#include "../unique_table.h"

// TODO: Add option to create temporary nodes of max size when
//       accumulating minterms.

// ********************************** MTMDDs **********************************

MEDDLY::mtmdd_forest::mtmdd_forest(int dsl, domain *d,
    bool relation, forest::range_type t,
    forest::edge_labeling e, const policies &p)
: mt_forest(dsl, d, relation, t, e, p)
{
  list = 0;
  termList = 0;
  listSize = 0;
  count = 0;
  slot = 0;
  countSize = 0;
}



MEDDLY::mtmdd_forest::~mtmdd_forest()
{
  if (list) free(list);
  if (termList) free(termList);
  if (count) free(count);
  if (slot) free(slot);
}


void MEDDLY::mtmdd_forest::expandCountAndSlotArrays(int size)
{
  if (size <= countSize) return;

  int newCountSize = countSize == 0? 8: countSize;
  while (newCountSize < size) { newCountSize *= 2; }

  count = (int*) realloc(count, newCountSize * sizeof(int));
  slot = (int*) realloc(slot, newCountSize * sizeof(int));
  memset(count + countSize, 0, (newCountSize - countSize) * sizeof(int));
  memset(slot + countSize, 0, (newCountSize - countSize) * sizeof(int));
  countSize = newCountSize;
}

int MEDDLY::mtmdd_forest::createNode(int k, int index, int dptr)
{
  MEDDLY_DCASSERT(index >= -1);

  if (index > -1 && getLevelSize(k) <= index) {
    useExpertDomain()->enlargeVariableBound(k, false, index + 1);
  }

  if (dptr == 0) return 0;
  if (index == -1) {
    // all downpointers should point to dptr
    if (isFullyReduced()) return dptr;
    insertRedundantNode(k, dptr);
    return dptr;
  }

  node_builder& nb = useSparseBuilder(k, 1);
  nb.d(0) = dptr;
  nb.i(0) = index;
  return createReducedNode(-1, nb);
}


void MEDDLY::mtmdd_forest::createEdge(const int* v, int term, dd_edge& e)
{
  // construct the edge bottom-up
  MEDDLY_DCASSERT(isTerminalNode(term));
  int result = term;
  int curr = 0;
  for (int i=1; i<=getExpertDomain()->getNumVariables(); i++) {
    result = createNode(i, v[i], result);
  }
  e.set(result, 0, getNodeLevel(result));
  // e.show(stderr, 2);
}


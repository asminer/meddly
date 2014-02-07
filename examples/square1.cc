
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

/*
    Program to analyze the "Square 1" puzzle.

    Organized as follows.
    Top has 12 slots, but some pieces take up two consecutive slots.
    Bottom is the same way.
    We can do an "exchange" twist if there are no double-wide pieces
    in the critical slots.
*/

#include "meddly.h"
#include "meddly_expert.h"
#include "timer.h"

using namespace MEDDLY;

/*
  Builds the "constraint" va == vb,
  for the given variables (levels).
  Use negative to denote primed.
  Assumes FULLY reduced.
*/
void Equals(int va, int vb, dd_edge &answer)
{
  if (ABS(va) < ABS(vb)) SWAP(va, vb);

  expert_forest* EF = (expert_forest*) answer.getForest();

  int asz = EF->getLevelSize(va);
  int bsz = EF->getLevelSize(vb);
  int size = MIN(asz, bsz);

  node_builder& na = EF->useNodeBuilder(va, size);

  for (int i=0; i<size; i++) {
    node_builder& nb = EF->useSparseBuilder(vb, 1);
    nb.i(0) = i;
    nb.d(0) = EF->handleForValue(1);
    na.d(i) = EF->createReducedNode(-1, nb);
  } // for i

  answer.set(EF->createReducedNode(-1, na), 0);
}

/*
  Builds the "constraint" va == vb,
  for the given variables (levels).
  Use negative to denote primed.
  Assumes FULLY reduced.
*/
void NotEquals(int va, int vb, dd_edge &answer)
{
  if (ABS(va) < ABS(vb)) SWAP(va, vb);

  expert_forest* EF = (expert_forest*) answer.getForest();

  int asz = EF->getLevelSize(va);
  int bsz = EF->getLevelSize(vb);

  node_builder& na = EF->useNodeBuilder(va, asz);

  for (int i=0; i<asz; i++) {
    node_builder& nb = EF->useNodeBuilder(vb, bsz);
    for (int j=0; j<asz; j++) {
      nb.d(j) = EF->handleForValue(i != j);
    }
    na.d(i) = EF->createReducedNode(-1, nb);
  } // for i

  answer.set(EF->createReducedNode(-1, na), 0);
}

/*
  "Rotates" a set of variables 
    x0' == x0+r mod n
    x1' == x1+r mod n
    x2' == x2+r mod n
    ...
    xn-1' == xn-1+r mod n

  Assumes FULLY reduced.
*/
void Rotate(const int* vars, int nv, int r, dd_edge &answer)
{
  expert_forest* EF = (expert_forest*) answer.getForest();

  answer.set(EF->handleForValue(1), 0);

  for (int i=0; i<nv; i++) {
    int j = (i+r) % nv;
    dd_edge temp(EF);
    Equals(-vars[i], vars[j], temp);
    answer *= temp;
  }
}


int main()
{
  MEDDLY::initialize();

  // Build the forest

  // Testing from here...

  int* scratch = new int[8+1];
  for (int i=0; i<=8; i++) scratch[i] = 5;

  domain* D = createDomainBottomUp(scratch, 8);
  forest::policies p(true);
  p.setFullyReduced();
  forest* F = D->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL, p);

  dd_edge test(F);

  int evens[] = {2, 4, 6, 8};
  Rotate(evens, 4, 1, test);

  test.show(stdout, 2);

  MEDDLY::cleanup();
  return 0;
}

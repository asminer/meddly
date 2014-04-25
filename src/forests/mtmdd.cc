
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

MEDDLY::mtmdd_forest
::mtmdd_forest(int dsl, domain* d, range_type t, const policies &p)
 : mt_forest(dsl, d, false, t, p)
{
  // anything to construct?
}

// ******************************************************************
// *                                                                *
// *              mtmdd_forest::mtmdd_iterator methods              *
// *                                                                *
// ******************************************************************

MEDDLY::mtmdd_forest::mtmdd_iterator::mtmdd_iterator(const expert_forest *F)
 : mt_iterator(F)
{
}

MEDDLY::mtmdd_forest::mtmdd_iterator::~mtmdd_iterator()
{
}

bool MEDDLY::mtmdd_forest::mtmdd_iterator::start(const dd_edge &e)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH);
  }
  return first(maxLevel, e.getNode());
}

bool MEDDLY::mtmdd_forest::mtmdd_iterator::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(!F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  int k;
  node_handle down = 0;
  for (k=1; k<=maxLevel; k++) {
    nzp[k]++;
    if (nzp[k] < path[k].getNNZs()) {
      index[k] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      break;
    }
  }
  level_change = k;
  if (k>maxLevel) {
    return false;
  }

  return first(k-1, down);
}

bool MEDDLY::mtmdd_forest::mtmdd_iterator::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(!F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  if (0==down) return false;

  for ( ; k; k--) {
    MEDDLY_DCASSERT(down);
    int kdn = F->getNodeLevel(down);
    MEDDLY_DCASSERT(kdn <= k);
    if (kdn < k)  F->initRedundantReader(path[k], k, down, false);
    else          F->initNodeReader(path[k], down, false);
    nzp[k] = 0;
    index[k] = path[k].i(0);
    down = path[k].d(0);
  }
  // save the terminal value
  index[0] = down;
  return true;
}

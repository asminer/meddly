
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

#include "mtmxd.h"

MEDDLY::mtmxd_forest
::mtmxd_forest(int dsl, domain* d, range_type t, const policies &p)
 : mt_forest(dsl, d, true, t, p)
{
  // anything to construct?
}

// ******************************************************************
// *                                                                *
// *              mtmxd_forest::mtmxd_iterator methods              *
// *                                                                *
// ******************************************************************

MEDDLY::mtmxd_forest::mtmxd_iterator::mtmxd_iterator(const expert_forest *F)
 : mt_iterator(F)
{
}

void MEDDLY::mtmxd_forest::swapAdjacentVariables(int level)
{
	throw error(error::NOT_IMPLEMENTED);
}

void MEDDLY::mtmxd_forest::moveDownVariable(int high, int low)
{
	throw error(error::NOT_IMPLEMENTED);
}

MEDDLY::mtmxd_forest::mtmxd_iterator::~mtmxd_iterator()
{
}

bool MEDDLY::mtmxd_forest::mtmxd_iterator::start(const dd_edge &e)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH);
  }
  return first(maxLevel, e.getNode());
}

bool MEDDLY::mtmxd_forest::mtmxd_iterator::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  int k = -1;
  node_handle down = 0;
  for (;;) { 
    nzp[k]++;
    if (nzp[k] < path[k].getNNZs()) {
      index[k] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      break;
    }
    if (k<0) {
      k = -k;
    } else {
      if (maxLevel == k) {
        level_change = k;
        return false;
      }
      k = -k-1;
    }
  }
  level_change = k;

  return first( (k>0) ? -k : -k-1, down);
}

bool MEDDLY::mtmxd_forest::mtmxd_iterator::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  if (0==down) return false;

  bool isFully = F->isFullyReduced();

  for ( ; k; k = downLevel(k) ) {
    MEDDLY_DCASSERT(down);
    int kdn = F->getNodeLevel(down);
    MEDDLY_DCASSERT(!isLevelAbove(kdn, k));

    if (isLevelAbove(k, kdn)) {
      if (k>0 || isFully) {
        F->initRedundantReader(path[k], k, down, false);
      } else {
        F->initIdentityReader(path[k], k, index[-k], down, false);
      }
    } else {
      F->initNodeReader(path[k], down, false);
    }
    nzp[k] = 0;
    index[k] = path[k].i(0);
    down = path[k].d(0);
  }
  // save the terminal value
  index[0] = down;
  return true;
}

// ******************************************************************
// *                                                                *
// *           mtmxd_forest::mtmxd_fixedrow_iter  methods           *
// *                                                                *
// ******************************************************************

MEDDLY::mtmxd_forest::
mtmxd_fixedrow_iter::mtmxd_fixedrow_iter(const expert_forest *F)
 : mt_iterator(F)
{
}

MEDDLY::mtmxd_forest::mtmxd_fixedrow_iter::~mtmxd_fixedrow_iter()
{
}

bool MEDDLY::mtmxd_forest::mtmxd_fixedrow_iter
::start(const dd_edge &e, const int* minterm)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH);
  }
  for (int k=1; k<=maxLevel; k++) {
    index[k] = minterm[k];
  }
  return first(maxLevel, e.getNode());
}

bool MEDDLY::mtmxd_forest::mtmxd_fixedrow_iter::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  node_handle down = 0;
  // Only try to advance the column, because the row is fixed.
  for (int k=-1; k>=-maxLevel; k--) { 
    for (nzp[k]++; nzp[k] < path[k].getNNZs(); nzp[k]++) {
      index[k] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      level_change = k;
      if (first(downLevel(k), down)) return true;
    }
  } // for

  return false;
}


bool MEDDLY::mtmxd_forest::mtmxd_fixedrow_iter::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  if (0==k) {
    index[0] = down;
    return true;
  }

  // Check that this "row" node has a non-zero pointer
  // for the fixed index.
  MEDDLY_DCASSERT(k>0);
  int cdown;
  if (isLevelAbove(k, F->getNodeLevel(down))) {
    // skipped unprimed level, must be "fully" reduced
    cdown = down;
  } else {
    cdown = F->getDownPtr(down, index[k]);
  }
  if (0==cdown) return false;

  //
  // Ok, set up the "column" node below
  k = downLevel(k);
  MEDDLY_DCASSERT(k<0);

  if (isLevelAbove(k, F->getNodeLevel(cdown))) {
    // Skipped level, we can be fast about this.
    // first, recurse.
    if (!first(downLevel(k), cdown)) return false;
    // Ok, there is a valid path.
    // Set up this level.
    nzp[k] = 0;
    if (F->isFullyReduced()) {
      F->initRedundantReader(path[k], k, cdown, false);
      index[k] = 0;
    } else {
      index[k] = index[upLevel(k)];
      F->initIdentityReader(path[k], k, index[k], cdown, false);
    }
    return true;
  } 

  // Proper node here.
  // cycle through it and recurse... 

  F->initNodeReader(path[k], cdown, false);

  for (int z=0; z<path[k].getNNZs(); z++) {
    if (first(downLevel(k), path[k].d(z))) {
      nzp[k] = z;
      index[k] = path[k].i(z);
      return true;
    }
  }

  return false;
}

// ******************************************************************
// *                                                                *
// *           mtmxd_forest::mtmxd_fixedcol_iter  methods           *
// *                                                                *
// ******************************************************************

MEDDLY::mtmxd_forest::
mtmxd_fixedcol_iter::mtmxd_fixedcol_iter(const expert_forest *F)
 : mt_iterator(F)
{
}

MEDDLY::mtmxd_forest::mtmxd_fixedcol_iter::~mtmxd_fixedcol_iter()
{
}

bool MEDDLY::mtmxd_forest::mtmxd_fixedcol_iter
::start(const dd_edge &e, const int* minterm)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH);
  }
  
  for (int k=1; k<=maxLevel; k++) {
    index[-k] = minterm[k];
  }

  return first(maxLevel, e.getNode());
}

bool MEDDLY::mtmxd_forest::mtmxd_fixedcol_iter::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  node_handle down = 0;
  // Only try to advance the row, because the column is fixed.
  for (int k=1; k<=maxLevel; k++) { 
    for (nzp[k]++; nzp[k] < path[k].getNNZs(); nzp[k]++) {
      index[k] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      level_change = k;
      if (first(downLevel(k), down)) return true;
    }
  } // for

  return false;
}


bool MEDDLY::mtmxd_forest::mtmxd_fixedcol_iter::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  if (0==k) {
    index[0] = down;
    return true;
  }

  if (k<0) {
    // See if this "column node" has a path
    // at the specified index.
    if (isLevelAbove(k, F->getNodeLevel(down))) {
      if (!F->isFullyReduced()) {
        // Identity node here - check index
        if (index[k] != index[upLevel(k)]) return false;
      }
      return first(downLevel(k), down);
    }
    int cdown = F->getDownPtr(down, index[k]);
    if (0==cdown) return false;
    return first(downLevel(k), cdown);
  }

  // Row node.  Find an index, if any,
  // such that there is a valid path below.
  MEDDLY_DCASSERT(k>0);
  int kdn = F->getNodeLevel(down);
  if (isLevelAbove(k, kdn)) {
    // Skipped level, handle quickly
    int kpr = downLevel(k);
    if (isLevelAbove(kpr, F->getNodeLevel(kdn))) {
      // next level is also skipped.
      // See if there is a valid path below.
      if (!first(downLevel(kpr), down)) return false;
      // There's one below, set up the one at these levels.
      F->initRedundantReader(path[k], k, down, false);
      if (F->isFullyReduced()) {
        nzp[k] = 0;
        index[k] = 0;
      } else {
        nzp[k] = index[kpr];
        index[k] = index[kpr];
      }
      return true;
    }
    // next level is not skipped.
    // See if there is a valid path below.
    int cdown = F->getDownPtr(down, index[kpr]);
    if (0==cdown) return false;
    if (!first(kpr, cdown)) return false;
    F->initRedundantReader(path[k], k, down, false);
    nzp[k] = 0;
    index[k] = 0;
    return true;
  }

  // Level is not skipped.
  F->initNodeReader(path[k], down, false);
  
  for (int z=0; z<path[k].getNNZs(); z++) {
    index[k] = path[k].i(z);
    if (first(downLevel(k), path[k].d(z))) {
      nzp[k] = z;
      return true;
    }
  }
  return false;
}



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



#include "defines.h"

// #define DEBUG_CLEANUP

// #define DEBUG_ITER_BEGIN

// ******************************************************************
// *                                                                *
// *                                                                *
// *                  enumerator::iterator methods                  *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::enumerator::iterator::iterator(const expert_forest* f)
{
  prindex = 0;
  F = f;
  if (0==f) {
    rawpath = path = 0;
    rawnzp = nzp = 0;
    rawindex = index = 0;
    minLevel = 0;
    maxLevel = 0;
    level_change = 0;
    return;
  }
  int N = f->getNumVariables();
  maxLevel = N;
  if (f->isForRelations()) {
    rawpath = new node_reader[2*N+1];
    rawnzp = new int[2*N+1];
    rawindex = new int[2*N+1];
    path = rawpath + N;
    nzp = rawnzp + N;
    index = rawindex + N;
    minLevel = -N;
  } else {
    rawpath = new node_reader[N+1];
    rawnzp = new int[N+1];
    rawindex = new int[N+1];
    path = rawpath;
    nzp = rawnzp;
    index = rawindex;
    minLevel = 1;
  }
  level_change = N+1;
}

MEDDLY::enumerator::iterator::~iterator()
{
  delete[] rawpath;
  delete[] rawnzp;
  delete[] rawindex;
  delete[] prindex;
}

bool MEDDLY::enumerator::iterator::start(const dd_edge &e)
{
  throw error(error::INVALID_OPERATION);
}

bool MEDDLY::enumerator::iterator::start(const dd_edge &e, const int* m)
{
  throw error(error::INVALID_OPERATION);
}

const int* MEDDLY::enumerator::iterator::getPrimedAssignments()
{
  if (0==F) return 0;
  if (!F->isForRelations()) return 0;
  if (0==prindex) {
    prindex = new int[1+maxLevel];
    prindex[0] = 0;
  }
  MEDDLY_DCASSERT(index);
  for (int k=maxLevel; k; k--) {
    prindex[k] = index[-k];
  }
  return prindex; 
}

void MEDDLY::enumerator::iterator::getValue(int &) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::enumerator::iterator::getValue(float &) const
{
  throw error(error::TYPE_MISMATCH);
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                       enumerator methods                       *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::enumerator::enumerator()
{
  I = 0;
  is_valid = false;
  T = EMPTY;
}

MEDDLY::enumerator::enumerator(type t, const forest* F)
{
  I = 0;
  init(t, F);
  is_valid = false;
}

MEDDLY::enumerator::enumerator(const dd_edge &e)
{
  I = 0;
  init(FULL, e.getForest());
  start(e);
}

MEDDLY::enumerator::~enumerator()
{
  delete I;
}

void MEDDLY::enumerator::init(type t, const forest* f)
{
  delete I;
  is_valid = false;

  T = t;
  const expert_forest* F = smart_cast<const expert_forest*>(f);
  if (0==F) {
    I = 0;
    return;
  }

  switch (t) {
    case FULL:
      I = F->makeFullIter();
      break;

    case ROW_FIXED:
      I = F->makeFixedRowIter();
      break;

    case COL_FIXED:
      I = F->makeFixedColumnIter();
      break;

    default:
      I = 0;
      return;
  }
  if (I->build_error()) {
    delete I;
    I = 0;
  } 
}

void MEDDLY::enumerator::start(const dd_edge &e)
{
  if (0==I) return;
  if (FULL != T) throw error(error::MISCELLANEOUS);
  MEDDLY_DCASSERT(I);
  is_valid = I->start(e);
}

void MEDDLY::enumerator::startFixedRow(const dd_edge &e, const int* minterm)
{
  if (0==I) return;
  if (ROW_FIXED != T) throw error(error::MISCELLANEOUS);
  MEDDLY_DCASSERT(I);
  is_valid = I->start(e, minterm);
}

void MEDDLY::enumerator::startFixedColumn(const dd_edge &e, const int* minterm)
{
  if (0==I) return;
  if (COL_FIXED != T) throw error(error::MISCELLANEOUS);
  MEDDLY_DCASSERT(I);
  is_valid = I->start(e, minterm);
}

#if 0

// OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
// OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
// OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
// OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
// OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
// OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
// OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
// OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
// OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD

// ******************************************************************
// *                                                                *
// *                                                                *
// *                       enumerator methods                       *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::enumerator::enumerator() 
{
  prindex = 0;
  initEmpty();
}

MEDDLY::enumerator::enumerator(const dd_edge &e)
{
  prindex = 0;
  initEmpty();
  start(e);
}

MEDDLY::enumerator::enumerator(const dd_edge &e, const int* allvars)
{
  prindex = 0;
  initEmpty();
  startFixed(e, allvars);
}

MEDDLY::enumerator::~enumerator()
{
  destroy();
}

void MEDDLY::enumerator::destroy()
{
  delete[] rawpath;
  delete[] rawnzp;
  delete[] rawindex;
  delete[] prindex;
  rawpath = 0;
  rawnzp = 0;
  rawindex = 0;
  prindex = 0;
}

void MEDDLY::enumerator::initEmpty()
{
  F = 0;
  type = EMPTY;
  rawpath = path = 0;
  rawnzp = nzp = 0;
  rawindex = index = 0;
  minLevel = maxLevel = 0;
  isValid = false;
  incr = 0;
  level_change = 0;
}

void MEDDLY::enumerator::newForest(expert_forest* f)
{
  MEDDLY_DCASSERT(f);
  if (f==F) return;
  destroy();
  F = f;
  int N = f->getNumVariables();
  maxLevel = N;
  if (f->isForRelations()) {
    rawpath = new node_reader[2*N+1];
    rawnzp = new int[2*N+1];
    rawindex = new int[2*N+1];
    path = rawpath + N;
    nzp = rawnzp + N;
    index = rawindex + N;
    minLevel = -N;
  } else {
    rawpath = new node_reader[N+1];
    rawnzp = new int[N+1];
    rawindex = new int[N+1];
    path = rawpath;
    nzp = rawnzp;
    index = rawindex;
    minLevel = 1;
  }
  level_change = N+1;
}

void MEDDLY::enumerator::start(const dd_edge &_e)
{
  e = _e;
  if (0==e.getForest() || 0==e.getNode()) {
    isValid = false;
    return;
  }
  newForest(smart_cast <expert_forest*>(e.getForest()));
  
  if (F->isForRelations()) {
    type = RELATION;
    isValid = firstRelElement(maxLevel, e.getNode());
    incr = &enumerator::incrRelation;
  } else {
    type = SET;
    isValid = firstSetElement(maxLevel, e.getNode());
    incr = &enumerator::incrNonRelation;
  }
}

void MEDDLY::enumerator::startFixedRow(const dd_edge &_e, const int* minterm)
{
  e = _e;
  if (0==e.getForest() || 0==e.getNode()) {
    isValid = false;
    return;
  }
  newForest(smart_cast <expert_forest*>(e.getForest()));
  
  if (!F->isForRelations()) throw error(error::TYPE_MISMATCH);

  type = ROW;
  for (int k=1; k<=maxLevel; k++) {
    index[k] = minterm[k];
  }

  isValid = firstColumn(maxLevel, e.getNode());
  incr = &enumerator::incrColumn;
}

void MEDDLY::enumerator::startFixedColumn(const dd_edge &_e, const int* minterm)
{
  e = _e;
  if (0==e.getForest() || 0==e.getNode()) {
    isValid = false;
    return;
  }
  newForest(smart_cast <expert_forest*>(e.getForest()));
  
  if (!F->isForRelations()) throw error(error::TYPE_MISMATCH);

  type = COLUMN;
  for (int k=1; k<=maxLevel; k++) {
    index[-k] = minterm[k];
  }

  isValid = firstRow(maxLevel, e.getNode());
  incr = &enumerator::incrRow;
}

void MEDDLY::enumerator::startFixed(const dd_edge &_e, const int* allvars)
{
  e = _e;
  throw error(error::NOT_IMPLEMENTED);
}

const int* MEDDLY::enumerator::getPrimedAssignments()
{
  if (!F->isForRelations()) return 0;
  if (0==prindex) {
    prindex = new int[1+maxLevel];
    prindex[0] = 0;
  }
  for (int k=maxLevel; k; k--) {
    prindex[k] = index[-k];
  }
  return prindex; 
}

void MEDDLY::enumerator::getValue(int& val) const
{
  MEDDLY_DCASSERT(F);
  if (F->isMultiTerminal()) {
    val = expert_forest::int_Tencoder::handle2value(index[0]);
  } else {
    // return the value using evaluate
    F->evaluate(e, index, val);
  }
}


void MEDDLY::enumerator::getValue(float& val) const
{
  MEDDLY_DCASSERT(F);
  if (F->isMultiTerminal()) {
    val = expert_forest::float_Tencoder::handle2value(index[0]);
  } else {
    // return the value using evaluate
    F->evaluate(e, index, val);
  }
}



bool MEDDLY::enumerator::incrNonRelation()
{
  MEDDLY_DCASSERT(SET == type);
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(!F->isForRelations());

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

  return firstSetElement(k-1, down);
}


bool MEDDLY::enumerator::incrRelation()
{
  MEDDLY_DCASSERT(RELATION == type);
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());

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

  return firstRelElement( (k>0) ? -k : -k-1, down);
}

bool MEDDLY::enumerator::incrRow()
{
  MEDDLY_DCASSERT(COLUMN == type);  // Column is fixed
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());

  node_handle down = 0;
  // Only try to advance the row, because the column is fixed.
  for (int k=1; k<=maxLevel; k++) { 
    for (nzp[k]++; nzp[k] < path[k].getNNZs(); nzp[k]++) {
      index[k] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      level_change = k;
      if (firstRow(downLevel(k), down)) return true;
    }
  } // for

  return false;
}

bool MEDDLY::enumerator::incrColumn()
{
  MEDDLY_DCASSERT(ROW == type);   // Row is fixed
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());

  node_handle down = 0;
  // Only try to advance the column, because the row is fixed.
  for (int k=-1; k>=-maxLevel; k--) { 
    for (nzp[k]++; nzp[k] < path[k].getNNZs(); nzp[k]++) {
      index[k] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      level_change = k;
      if (firstColumn(downLevel(k), down)) return true;
    }
  } // for

  return false;
}

bool MEDDLY::enumerator::firstSetElement(int k, node_handle down)
{
  if (0==down) return false;
  MEDDLY_DCASSERT(SET == type);
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(!F->isForRelations());

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

bool MEDDLY::enumerator::firstRelElement(int k, node_handle down)
{
  if (0==down) return false;
  MEDDLY_DCASSERT(RELATION == type);
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());

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

bool MEDDLY::enumerator::firstRow(int k, node_handle down)
{
  if (0==k) {
    index[0] = down;
    return true;
  }
  MEDDLY_DCASSERT(COLUMN == type);  // Column is fixed
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());

  if (k<0) {
    // See if this "column node" has a path
    // at the specified index.
    if (isLevelAbove(k, F->getNodeLevel(down))) {
      if (!F->isFullyReduced()) {
        // Identity node here - check index
        if (index[k] != index[upLevel(k)]) return false;
      }
      return firstRow(downLevel(k), down);
    }
    int cdown = F->getDownPtr(down, index[k]);
    if (0==cdown) return false;
    return firstRow(downLevel(k), cdown);
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
      if (!firstRow(downLevel(kpr), down)) return false;
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
    if (!firstRow(kpr, cdown)) return false;
    F->initRedundantReader(path[k], k, down, false);
    nzp[k] = 0;
    index[k] = 0;
    return true;
  }

  // Level is not skipped.
  F->initNodeReader(path[k], down, false);
  
  for (int z=0; z<path[k].getNNZs(); z++) {
    index[k] = path[k].i(z);
    if (firstRow(downLevel(k), path[k].d(z))) {
      nzp[k] = z;
      return true;
    }
  }
  return false;
}

bool MEDDLY::enumerator::firstColumn(int k, node_handle down)
{
  if (0==k) {
    index[0] = down;
    return true;
  }
  MEDDLY_DCASSERT(ROW == type);   // Row is fixed
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());

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
    if (!firstColumn(downLevel(k), cdown)) return false;
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
    if (firstColumn(downLevel(k), path[k].d(z))) {
      nzp[k] = z;
      index[k] = path[k].i(z);
      return true;
    }
  }

  return false;
}

#endif


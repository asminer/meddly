
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

// Helper functions
inline void linkNode(MEDDLY::forest* p, int node)
{
  MEDDLY_DCASSERT(p);
  MEDDLY_DCASSERT(smart_cast<MEDDLY::expert_forest*>(p));
  smart_cast<MEDDLY::expert_forest*>(p)->linkNode(node);
}

inline void unlinkNode(MEDDLY::forest* p, int node)
{
  if (p) {
    MEDDLY_DCASSERT(smart_cast<MEDDLY::expert_forest*>(p));
    smart_cast<MEDDLY::expert_forest*>(p)->unlinkNode(node);
  }
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        dd_edge  methods                        *
// *                                                                *
// *                                                                *
// ******************************************************************


// Constructor.
MEDDLY::dd_edge::dd_edge(forest* p)
: parent(p),
  node(0), value(0), level(0), index(-1),
  opPlus(0), opStar(0), opMinus(0), opDivide(0),
  updateNeeded(true), beginIterator(0)
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating dd_edge %x\n", this);
#endif
  MEDDLY_DCASSERT(p != NULL);
  parent->registerEdge(*this);
  MEDDLY_DCASSERT(index != -1);
}


// Copy Constructor.
MEDDLY::dd_edge::dd_edge(const dd_edge& e)
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating dd_edge %x\n", this);
#endif
  init(e);
  MEDDLY_DCASSERT(index != -1);
}


// Assignment operator.
MEDDLY::dd_edge& MEDDLY::dd_edge::operator=(const dd_edge& e)
{
  if (&e != this) {
    destroy();
    init(e);
  }
  MEDDLY_DCASSERT(index != -1);
  return *this;
}

// Destructor.  Will notify parent as appropriate.
MEDDLY::dd_edge::~dd_edge()
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Deleting dd_edge %x\n", this);
#endif
  destroy();
}

//
// Initialization helper
//
void MEDDLY::dd_edge::init(const dd_edge &e)
{
  parent = e.parent;
  node = e.node;
  value = e.value;
  level = e.level;

  linkNode(parent, node);

  updateNeeded = e.updateNeeded;
  
  opPlus = e.opPlus;
  opStar = e.opStar;
  opMinus = e.opMinus;
  opDivide = e.opDivide;

  if (updateNeeded) {
    beginIterator = 0;
  } else {
    beginIterator = new const_iterator(*(e.beginIterator));
  }

  if (parent) parent->registerEdge(*this);
  MEDDLY_DCASSERT(index != -1);
}

//
// destruction helper
//
void MEDDLY::dd_edge::destroy()
{
  if (beginIterator != 0) delete beginIterator;
  if (index != -1) {
    // still registered; unregister before discarding
    int old = node;
    node = 0;
    unlinkNode(parent, old);
    if (parent) parent->unregisterEdge(*this);
  }
}

void MEDDLY::dd_edge::getEdgeValue(float& ev) const
{
  ev = toFloat(value);
}


void MEDDLY::dd_edge::set(int n, int v, int l)
{
  if (node != n) { updateNeeded = true; }
  int old = node;
  node = n;
  unlinkNode(parent, old);
  value = v;
  level = l;
}


void MEDDLY::dd_edge::set(int n, float v, int l)
{
  if (node != n) { updateNeeded = true; }
  int old = node;
  node = n;
  unlinkNode(parent, old);
  value = toInt(v);
  level = l;
}


unsigned MEDDLY::dd_edge::getNodeCount() const
{
  return smart_cast<expert_forest*>(parent)->getNodeCount(node);
}

unsigned MEDDLY::dd_edge::getEdgeCount(bool countZeroes) const
{
  return smart_cast<expert_forest*>(parent)->getEdgeCount(node, countZeroes);
}

//
// Operator +=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator+=(const dd_edge& e)
{
  if (opPlus == 0) {
    if (parent->getRangeType() == forest::BOOLEAN)
      opPlus = getOperation(UNION, *this, e, *this);
    else
      opPlus = getOperation(PLUS, *this, e, *this);
    MEDDLY_DCASSERT(opPlus != 0);
  }
  opPlus->compute(*this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  MEDDLY_DCASSERT(updateNeeded == true);
  return *this;
}


// Operator *=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator*=(const dd_edge& e)
{
  if (opStar == 0) {
    if (parent->getRangeType() == forest::BOOLEAN)
      opStar = getOperation(INTERSECTION, *this, e, *this);
    else
      opStar = getOperation(MULTIPLY, *this, e, *this);
    MEDDLY_DCASSERT(opStar != 0);
  }
  opStar->compute(*this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  MEDDLY_DCASSERT(updateNeeded == true);
  return *this;
}


// Operator -=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator-=(const dd_edge& e)
{
  if (opMinus == 0) {
    if (parent->getRangeType() == forest::BOOLEAN)
      opMinus = getOperation(DIFFERENCE, *this, e, *this);
    else
      opMinus = getOperation(MINUS, *this, e, *this);
    MEDDLY_DCASSERT(opMinus != 0);
  }
  opMinus->compute(*this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  MEDDLY_DCASSERT(updateNeeded == true);
  return *this;
}


// Operator /=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator/=(const dd_edge& e)
{
  if (opDivide == 0) {
    opDivide = getOperation(DIVIDE, *this, e, *this);
  }
  opDivide->compute(*this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  MEDDLY_DCASSERT(updateNeeded == true);
  return *this;
}


// Display the edge information.
void MEDDLY::dd_edge::show(FILE* strm, int verbosity) const
{
  expert_forest* eParent = smart_cast<expert_forest*>(parent);
  fprintf(strm, "(Forest Addr: %p, ", parent);
  if (eParent->isTerminalNode(node)) {
    if (eParent->getRangeType() == forest::REAL) {
      fprintf(strm, "node: %f*, ", eParent->getReal(node));
    }
    else if (eParent->getRangeType() == forest::INTEGER) {
      fprintf(strm, "node: %d*, ", eParent->getInteger(node));
    }
    else {
      MEDDLY_DCASSERT(eParent->getRangeType() == forest::BOOLEAN);
      fprintf(strm, "node: %s*, ",(eParent->getBoolean(node)? "true": "false"));
    }
  }
  else {
    fprintf(strm, "node: %d, ", node);
  }
  if (eParent->getEdgeLabeling() == forest::EVTIMES) {
    // fprintf(strm, "value: %f, level: %d)\n", toFloat(value), level);
    float ev = 0;
    getEdgeValue(ev);
    fprintf(strm, "value: %f, level: %d)\n", ev, level);
  } else {
    fprintf(strm, "value: %d, level: %d)\n", value, level);
  }
  if (verbosity == 2 || verbosity == 3) {
    fprintf(strm, "MDD rooted at this node:\n");
    smart_cast<expert_forest*>(parent)->showNodeGraph(strm, node);
  }
  if (verbosity == 1 || verbosity == 3) {
    fprintf(strm, "Cardinality of node %d: %0.8e\n", node, getCardinality());
  }
}


MEDDLY::dd_edge::iterator MEDDLY::dd_edge::begin()
{
  if (updateNeeded) {
    updateIterators();
    updateNeeded = false;
  }
  MEDDLY_DCASSERT(beginIterator != 0);
  return *beginIterator;
}


MEDDLY::dd_edge::iterator MEDDLY::dd_edge::beginRow(const int* minterm)
{
  if (updateNeeded) {
    updateIterators();
    updateNeeded = false;
  }
  if (this->parent->isForRelations())
    return iterator(this, iterator::ROW, minterm);
  else
    return iterator();
}


MEDDLY::dd_edge::iterator MEDDLY::dd_edge::beginColumn(const int* minterm)
{
  if (updateNeeded) {
    updateIterators();
    updateNeeded = false;
  }
  // return iterator(this, false, minterm);
  return iterator(this, iterator::COLUMN, minterm);
}


void MEDDLY::dd_edge::updateIterators()
{
  // update beginIterator
  if (beginIterator != 0) delete beginIterator;
  
  if (parent->isForRelations()) {
    beginIterator = new iterator(this, iterator::RELATION, 0);
  } else {
    beginIterator = new iterator(this, iterator::SET, 0);
  }
}

// ******************************************************************
// *                                                                *
// *                   dd_edge::iterator  methods                   *
// *                                                                *
// ******************************************************************


MEDDLY::dd_edge::iterator::iterator()
{
  prindex = 0;
  initEmpty();
}

MEDDLY::dd_edge::iterator::iterator(dd_edge* _e, iter_type t, const int* minterm)
{
  prindex = 0;
  if (_e == 0) {
    initEmpty();
    return;
  }
  if (t == ROW || t == COLUMN) {
    if (!_e->parent->isForRelations()) {
      throw error(error::TYPE_MISMATCH);
    } 
  }
  e = _e;
  type = t;

  // Set up arrays
  int N = e->parent->getDomain()->getNumVariables();
  maxLevel = N;
  if (e->parent->isForRelations()) {
    rawpath = new void*[2*N+1];
    rawnzp = new int[2*N+1];
    rawindex = new int[2*N+1];
    for (int i=2*N; i>=0; i--) {
      rawpath[i] = 0;
      rawnzp[i] = -1;
      rawindex[i] = -1;
    }
    path = rawpath + N;
    nzp = rawnzp + N;
    index = rawindex + N;
    minLevel = -N;
  } else {
    rawpath = new void*[N+1];
    rawnzp = new int[N+1];
    rawindex = new int[N+1];
    for (int i=N; i>=0; i--) {
      rawpath[i] = 0;
      rawnzp[i] = -1;
      rawindex[i] = -1;
    }
    path = rawpath;
    nzp = rawnzp;
    index = rawindex;
    minLevel = 1;
  }

  // Go to the first element.
  switch (type) {
    case SET:
      isValid = firstSetElement(maxLevel, e->node);
      break;

    case RELATION:
      isValid = firstRelElement(maxLevel, e->node);
      break;

    case ROW:
      MEDDLY_DCASSERT(e->parent->isForRelations());
      // Row is fixed
      for (int k=1; k<=N; k++) {
        index[k] = minterm[k];
      }
      isValid = firstColumn(maxLevel, e->node);
      break;

    case COLUMN:
      MEDDLY_DCASSERT(e->parent->isForRelations());
      // Column is fixed
      for (int k=1; k<=N; k++) {
        index[-k] = minterm[k];
      }
      isValid = firstRow(maxLevel, e->node);
      break;

    default:
      throw error(error::MISCELLANEOUS);
  }
}

MEDDLY::dd_edge::iterator::iterator(const iterator& iter)
{
  prindex = 0;
  init(iter);
}

MEDDLY::dd_edge::iterator::~iterator()
{
  destroy();
}

MEDDLY::dd_edge::iterator& 
MEDDLY::dd_edge::iterator::operator=(const iterator& iter)
{
  if (this == &iter) return *this;
  destroy();
  init(iter);
  return *this;
}

void MEDDLY::dd_edge::iterator::destroy()
{
  if (e) {
    expert_forest* f = (expert_forest*) e->parent;
    MEDDLY_DCASSERT(f);
    for (int k=minLevel; k<=maxLevel; k++) {
      expert_forest::nodeReader* nlk = (expert_forest::nodeReader*) path[k];
      f->recycle(nlk);
      path[k] = 0;
    }
  }
  delete[] rawpath;
  delete[] rawnzp;
  delete[] rawindex;
  delete[] prindex;
  rawpath = 0;
  rawnzp = 0;
  rawindex = 0;
  prindex = 0;
}

void MEDDLY::dd_edge::iterator::initEmpty()
{
  e = 0;
  type = EMPTY;
  rawpath = path = 0;
  rawnzp = nzp = 0;
  rawindex = index = 0;
  minLevel = maxLevel = 0;
  isValid = false;
}

void MEDDLY::dd_edge::iterator::init(const iterator& iter)
{
  if (0==iter.e) {
    initEmpty();
    return;
  }
  MEDDLY_DCASSERT(iter.type != EMPTY);
  e = iter.e;
  type = iter.type;
  isValid = iter.isValid;

  // Set up arrays
  expert_forest* f = smart_cast <expert_forest*>(e->parent);
  MEDDLY_DCASSERT(f);
  int N = f->getNumVariables();
  maxLevel = N;
  if (f->isForRelations()) {
    rawpath = new void*[2*N+1];
    rawnzp = new int[2*N+1];
    rawindex = new int[2*N+1];
    for (int i=2*N; i>=0; i--) {
      rawpath[i] = f->copyNodeReader((expert_forest::nodeReader*)iter.rawpath[i]);
    }
    memcpy(rawnzp, iter.rawnzp, (2*N+1)*sizeof(int));
    memcpy(rawindex, iter.rawindex, (2*N+1)*sizeof(int));
    path = rawpath + N;
    nzp = rawnzp + N;
    index = rawindex + N;
    minLevel = -N;
  } else {
    rawpath = new void*[N+1];
    rawnzp = new int[N+1];
    rawindex = new int[N+1];
    for (int i=N; i>=0; i--) {
      rawpath[i] = f->copyNodeReader((expert_forest::nodeReader*)iter.rawpath[i]);
    }
    memcpy(rawnzp, iter.rawnzp, (N+1)*sizeof(int));
    memcpy(rawindex, iter.rawindex, (N+1)*sizeof(int));
    path = rawpath;
    nzp = rawnzp;
    index = rawindex;
    minLevel = 1;
  }
  MEDDLY_DCASSERT(iter.maxLevel == maxLevel);
  MEDDLY_DCASSERT(iter.minLevel == minLevel);
}

void MEDDLY::dd_edge::iterator::operator++()
{
  switch (type) {
    case EMPTY:
        return;

    case SET:
        MEDDLY_DCASSERT(e);
        MEDDLY_DCASSERT(isValid);
        isValid &= incrNonRelation();
        return;

    case RELATION:
        MEDDLY_DCASSERT(e);
        MEDDLY_DCASSERT(isValid);
        isValid &= incrRelation();
        return;

    case ROW:
        MEDDLY_DCASSERT(isValid);
        isValid &= incrColumn();
        return;

    case COLUMN:
        MEDDLY_DCASSERT(isValid);
        isValid &= incrRow();
        return;

    default:
        throw error(error::MISCELLANEOUS);
  };
}

bool MEDDLY::dd_edge::iterator::operator==(const iterator& iter) const
{
  if (e != iter.e) return false;
  if (type != iter.type) return false;
  if (e == 0) return true;
  MEDDLY_DCASSERT(iter.minLevel == minLevel);
  MEDDLY_DCASSERT(iter.maxLevel == maxLevel);
  for (int k=minLevel; k<=maxLevel; k++) {
    if (iter.index[k] != index[k]) return false;
  }
  return true;
}

const int* MEDDLY::dd_edge::iterator::getPrimedAssignments()
{
  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  if (!f->isForRelations()) return 0;
  if (0==prindex) {
    prindex = new int[1+maxLevel];
    prindex[0] = 0;
  }
  for (int k=maxLevel; k; k--) {
    prindex[k] = index[-k];
  }
  return prindex; 
}

void MEDDLY::dd_edge::iterator::getValue(int& val) const
{
  if (e == 0) return;
  if (e->parent->isMultiTerminal()) {
    val = static_cast<expert_forest*>(e->parent)->getInteger(index[0]);
  } else {
    // return the value using evaluate
    e->parent->evaluate(*e, index, val);
  }
}


void MEDDLY::dd_edge::iterator::getValue(float& val) const
{
  if (e == 0) return;
  if (e->parent->isMultiTerminal()) {
    val = static_cast<expert_forest*>(e->parent)->getReal(index[0]);
  } else {
    // return the value using evaluate
    e->parent->evaluate(*e, index, val);
  }
}



bool MEDDLY::dd_edge::iterator::incrNonRelation()
{
  MEDDLY_DCASSERT(e != 0);
  MEDDLY_DCASSERT(e->node != 0);
  MEDDLY_DCASSERT(type == SET);

  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  MEDDLY_DCASSERT(f);
  MEDDLY_DCASSERT(!f->isForRelations());

  int k;
  int down = 0;
  for (k=1; k<=maxLevel; k++) {
    nzp[k]++;
    expert_forest::nodeReader* nlk = (expert_forest::nodeReader*) path[k];
    if (nzp[k] < nlk->getNNZs()) {
      index[k] = nlk->i(nzp[k]);
      down = nlk->d(nzp[k]);
      MEDDLY_DCASSERT(down);
      break;
    }
    f->recycle(nlk);
    path[k] = 0;
  }
  if (k>maxLevel) {
    return false;
  }

  return firstSetElement(k-1, down);
}


bool MEDDLY::dd_edge::iterator::incrRelation()
{
  MEDDLY_DCASSERT(e != 0);
  MEDDLY_DCASSERT(e->node != 0);
  MEDDLY_DCASSERT(type == RELATION);

  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  MEDDLY_DCASSERT(f);
  MEDDLY_DCASSERT(f->isForRelations());

  int k = -1;
  int down = 0;
  for (;;) { 
    nzp[k]++;
    expert_forest::nodeReader* nlk = (expert_forest::nodeReader*) path[k];
    if (nzp[k] < nlk->getNNZs()) {
      index[k] = nlk->i(nzp[k]);
      down = nlk->d(nzp[k]);
      MEDDLY_DCASSERT(down);
      break;
    }
    f->recycle(nlk);
    path[k] = 0;
    if (k<0) {
      k = -k;
    } else {
      if (maxLevel == k) return false;
      k = -k-1;
    }
  }

  return firstRelElement( (k>0) ? -k : -k-1, down);
}

bool MEDDLY::dd_edge::iterator::incrRow()
{
  MEDDLY_DCASSERT(e != 0);
  MEDDLY_DCASSERT(e->node != 0);
  MEDDLY_DCASSERT(type == COLUMN);  // Column is fixed

  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  MEDDLY_DCASSERT(f);
  MEDDLY_DCASSERT(f->isForRelations());

  int down = 0;
  // Only try to advance the row, because the column is fixed.
  for (int k=1; k<=maxLevel; k++) { 
    expert_forest::nodeReader* nlk = (expert_forest::nodeReader*) path[k];
    for (nzp[k]++; nzp[k] < nlk->getNNZs(); nzp[k]++) {
      index[k] = nlk->i(nzp[k]);
      down = nlk->d(nzp[k]);
      MEDDLY_DCASSERT(down);
      if (firstRow(downLevel(k), down)) return true;
    }
    f->recycle(nlk);
    path[k] = 0;
  } // for

  return false;
}

bool MEDDLY::dd_edge::iterator::incrColumn()
{
  MEDDLY_DCASSERT(e != 0);
  MEDDLY_DCASSERT(e->node != 0);
  MEDDLY_DCASSERT(type == ROW);   // Row is fixed

  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  MEDDLY_DCASSERT(f);
  MEDDLY_DCASSERT(f->isForRelations());

  int down = 0;
  // Only try to advance the column, because the row is fixed.
  for (int k=-1; k>=-maxLevel; k--) { 
    expert_forest::nodeReader* nlk = (expert_forest::nodeReader*) path[k];
    for (nzp[k]++; nzp[k] < nlk->getNNZs(); nzp[k]++) {
      index[k] = nlk->i(nzp[k]);
      down = nlk->d(nzp[k]);
      MEDDLY_DCASSERT(down);
      if (firstColumn(downLevel(k), down)) return true;
    }
    f->recycle(nlk);
    path[k] = 0;
  } // for

  return false;
}

bool MEDDLY::dd_edge::iterator::firstSetElement(int k, int down)
{
  if (0==down) return false;
  MEDDLY_DCASSERT(type == SET);

  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  MEDDLY_DCASSERT(f);
  MEDDLY_DCASSERT(!f->isForRelations());

  for ( ; k; k--) {
    MEDDLY_DCASSERT(down);
    int kdn = f->getNodeLevel(down);
    MEDDLY_DCASSERT(kdn <= k);
    expert_forest::nodeReader* nlk = (kdn < k)
      ? f->initRedundantReader(k, down, false)
      : f->initNodeReader(down, false);
    path[k] = nlk;
    nzp[k] = 0;
    index[k] = nlk->i(0);
    down = nlk->d(0);
  }
  // save the terminal value
  index[0] = down;
  return true;
}

bool MEDDLY::dd_edge::iterator::firstRelElement(int k, int down)
{
  if (0==down) return false;
  MEDDLY_DCASSERT(type == RELATION);

  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  MEDDLY_DCASSERT(f);
  MEDDLY_DCASSERT(f->isForRelations());
  bool isFully = f->isFullyReduced();

  for ( ; k; k = downLevel(k) ) {
    MEDDLY_DCASSERT(down);
    int kdn = f->getNodeLevel(down);
    MEDDLY_DCASSERT(!isLevelAbove(kdn, k));

    expert_forest::nodeReader* nlk;
    if (isLevelAbove(k, kdn)) {
      if (k>0 || isFully) {
        nlk = f->initRedundantReader(k, down, false);
      } else {
        nlk = f->initIdentityReader(k, index[-k], down, false);
      }
    } else {
      nlk = f->initNodeReader(down, false);
    }
    path[k] = nlk;
    nzp[k] = 0;
    index[k] = nlk->i(0);
    down = nlk->d(0);
  }
  // save the terminal value
  index[0] = down;
  return true;
}

bool MEDDLY::dd_edge::iterator::firstRow(int k, int down)
{
  if (0==k) {
    index[0] = down;
    return true;
  }
  MEDDLY_DCASSERT(e);
  MEDDLY_DCASSERT(e->node);
  MEDDLY_DCASSERT(type == COLUMN);  // Column is fixed

  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  MEDDLY_DCASSERT(f);
  MEDDLY_DCASSERT(f->isForRelations());

  if (k<0) {
    // See if this "column node" has a path
    // at the specified index.
    if (isLevelAbove(k, f->getNodeLevel(down))) {
      if (!f->isFullyReduced()) {
        // Identity node here - check index
        if (index[k] != index[upLevel(k)]) return false;
      }
      return firstRow(downLevel(k), down);
    }
    int cdown = f->getDownPtr(down, index[k]);
    if (0==cdown) return false;
    return firstRow(downLevel(k), cdown);
  }

  // Row node.  Find an index, if any,
  // such that there is a valid path below.
  MEDDLY_DCASSERT(k>0);
  int kdn = f->getNodeLevel(down);
  if (isLevelAbove(k, kdn)) {
    // Skipped level, handle quickly
    int kpr = downLevel(k);
    if (isLevelAbove(kpr, f->getNodeLevel(kdn))) {
      // next level is also skipped.
      // See if there is a valid path below.
      if (!firstRow(downLevel(kpr), down)) return false;
      // There's one below, set up the one at these levels.
      path[k] = f->initRedundantReader(k, down, false);
      if (f->isFullyReduced()) {
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
    int cdown = f->getDownPtr(down, index[kpr]);
    if (0==cdown) return false;
    if (!firstRow(kpr, cdown)) return false;
    path[k] = f->initRedundantReader(k, down, false);
    nzp[k] = 0;
    index[k] = 0;
    return true;
  }

  // Level is not skipped.
  expert_forest::nodeReader* nlk = f->initNodeReader(down, false);
  
  for (int z=0; z<nlk->getNNZs(); z++) {
    index[k] = nlk->i(z);
    if (firstRow(downLevel(k), nlk->d(z))) {
      path[k] = nlk;
      nzp[k] = z;
      return true;
    }
  }

  f->recycle(nlk);
  return false;
}

bool MEDDLY::dd_edge::iterator::firstColumn(int k, int down)
{
  if (0==k) {
    index[0] = down;
    return true;
  }
  MEDDLY_DCASSERT(e);
  MEDDLY_DCASSERT(e->node);
  MEDDLY_DCASSERT(type == ROW);   // Row is fixed

  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  MEDDLY_DCASSERT(f);
  MEDDLY_DCASSERT(f->isForRelations());

  // Check that this "row" node has a non-zero pointer
  // for the fixed index.
  MEDDLY_DCASSERT(k>0);
  int cdown;
  if (isLevelAbove(k, f->getNodeLevel(down))) {
    // skipped unprimed level, must be "fully" reduced
    cdown = down;
  } else {
    cdown = f->getDownPtr(down, index[k]);
  }
  if (0==cdown) return false;

  //
  // Ok, set up the "column" node below
  k = downLevel(k);
  MEDDLY_DCASSERT(k<0);

  if (isLevelAbove(k, f->getNodeLevel(cdown))) {
    // Skipped level, we can be fast about this.
    // first, recurse.
    if (!firstColumn(downLevel(k), cdown)) return false;
    // Ok, there is a valid path.
    // Set up this level.
    nzp[k] = 0;
    if (f->isFullyReduced()) {
      path[k] = f->initRedundantReader(k, cdown, false);
      index[k] = 0;
    } else {
      index[k] = index[upLevel(k)];
      path[k] = f->initIdentityReader(k, index[k], cdown, false);
    }
    return true;
  } 

  // Proper node here.
  // cycle through it and recurse... 

  expert_forest::nodeReader* nlk = f->initNodeReader(cdown, false);

  for (int z=0; z<nlk->getNNZs(); z++) {
    if (firstColumn(downLevel(k), nlk->d(z))) {
      path[k] = nlk;
      nzp[k] = z;
      index[k] = nlk->i(z);
      return true;
    }
  }

  f->recycle(nlk);
  return false;
}


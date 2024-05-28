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

#include "evmxd_pluslong.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                    evmxd_pluslong  methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::evmxd_pluslong::evmxd_pluslong(domain *d, const policies &p, int* level_reduction_rule)
 : evmxd_forest(d, range_type::INTEGER, edge_labeling::EVPLUS,
         p, level_reduction_rule)
{
    setLongEdges();
    setTransparentEdge(0, long(0));
    initializeStorage();
}

MEDDLY::evmxd_pluslong::~evmxd_pluslong()
{ }

void MEDDLY::evmxd_pluslong::createEdge(long val, dd_edge &e)
{
  createEdgeTempl<OP, long>(val, e);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::evmxd_pluslong
::createEdge(const int* const* vlist, const int* const* vplist,
  const long* terms, int N, dd_edge &e)
{
  // XXX: Requires UnionPlus
  binary_operation* unionOp = PLUS(this, this, this);
  MEDDLY_DCASSERT(unionOp);
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);
  enlargeVariables(vplist, N, true);

  evmxd_edgemaker<OP, long>
  EM(this, vlist, vplist, terms, order, N,
    getDomain()->getNumVariables(), unionOp);

  long ev;
  node_handle ep;
  EM.createEdge(ev, ep);
  e.set(ep, ev);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::evmxd_pluslong
::createEdgeForVar(int vh, bool vp, const long* terms, dd_edge& a)
{
  createEdgeForVarTempl<OP, long>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::evmxd_pluslong
::evaluate(const dd_edge &f, const int* vlist, const int* vplist,
  long &term) const
{
  evaluateT<OP, long>(f, vlist, vplist, term);
}

void MEDDLY::evmxd_pluslong::showEdge(output &s, const edge_value &ev,
        node_handle d) const
{
    if (d == 0) {
        s.put("<oo, w>");
    } else {
        s.put('<');
        s.put(ev.getLong());
        s.put(", ");
        if (d < 0) {
            s.put('w');
        } else {
            s.put('#');
            s.put(d);
        }
        s.put('>');
    }
}

#ifdef ALLOW_DEPRECATED_0_17_3

const char* MEDDLY::evmxd_pluslong::codeChars() const
{
  return "dd_etxr";
}

#endif

// ******************************************************************
// *                                                                *
// *                                                                *
// *           evmxd_pluslong::evtrmxd_baseiter  methods           *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::evmxd_pluslong::evtrmxd_baseiter::evtrmxd_baseiter(const forest *F)
: iterator(F)
{
  int N = F->getNumVariables();
  raw_acc_evs = new long[2*N+1];
  acc_evs = raw_acc_evs + N;
}

MEDDLY::evmxd_pluslong::evtrmxd_baseiter::~evtrmxd_baseiter()
{
  delete[] raw_acc_evs;
}

void MEDDLY::evmxd_pluslong::evtrmxd_baseiter::getValue(long &tv) const
{
  MEDDLY_DCASSERT(acc_evs);
  tv = acc_evs[0];
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *           evmxd_pluslong::evtrmxd_iterator  methods           *
// *                                                                *
// *                                                                *
// ******************************************************************

bool MEDDLY::evmxd_pluslong::evtrmxd_iterator::start(const dd_edge &e)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
  }

  MEDDLY_DCASSERT(acc_evs);
  long ev;
  e.getEdgeValue(ev);
  acc_evs[maxLevel] = ev;

  return first(maxLevel, e.getNode());
}

bool MEDDLY::evmxd_pluslong::evtrmxd_iterator::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);
  MEDDLY_DCASSERT(acc_evs);

  int k = -1;
  node_handle down = 0;
  for (;;) {
    nzp[k]++;
    if (nzp[k] < path[k].getSize()) {
      index[k] = path[k].index(nzp[k]);
      down = path[k].down(nzp[k]);
      MEDDLY_DCASSERT(down);
      const long ev = path[k].edgeval(nzp[k]).getLong();
      acc_evs[downLevel(k)] = acc_evs[k] + ev;
      break;
    }
    if (maxLevel == k) {
      level_change = k+1;
      return false;
    }
    k = upLevel(k);
  } // infinite loop
  level_change = k;

  return first( (k>0) ? -k : -k-1, down);
}

bool MEDDLY::evmxd_pluslong::evtrmxd_iterator::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);
  MEDDLY_DCASSERT(acc_evs);

  if (0==down) return false;

  bool isFully = F->isFullyReduced();

  for ( ; k; k = downLevel(k) ) {
    MEDDLY_DCASSERT(down);
    int kdn = F->getNodeLevel(down);
    MEDDLY_DCASSERT(!isLevelAbove(kdn, k));

    if (isLevelAbove(k, kdn)) {
      if (k>0 || isFully) {
        path[k].initRedundant(F, k, 0L, down, SPARSE_ONLY);
      } else {
        path[k].initIdentity(F, k, index[-k], 0L, down, SPARSE_ONLY);
      }
    } else {
      F->unpackNode(path+k, down, SPARSE_ONLY);
    }
    nzp[k] = 0;
    index[k] = path[k].index(0);
    down = path[k].down(0);
    const long ev = path[k].edgeval(0).getLong();
    acc_evs[downLevel(k)] = acc_evs[k] + ev;
  }
  // save the terminal value
  index[0] = down;
  return true;
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *         evmxd_pluslong::evtrmxd_fixedrow_iter methods         *
// *                                                                *
// *                                                                *
// ******************************************************************

bool MEDDLY::evmxd_pluslong::evtrmxd_fixedrow_iter
::start(const dd_edge &e, const int* minterm)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
  }

  MEDDLY_DCASSERT(acc_evs);
  long ev;
  e.getEdgeValue(ev);
  acc_evs[maxLevel] = ev;

  for (int k=1; k<=maxLevel; k++) {
    index[k] = minterm[k];
  }
  return first(maxLevel, e.getNode());
}

bool MEDDLY::evmxd_pluslong::evtrmxd_fixedrow_iter::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);
  MEDDLY_DCASSERT(acc_evs);

  node_handle down = 0;
  // Only try to advance the column, because the row is fixed.
  for (int k=-1; k>=-maxLevel; k--) {
    for (nzp[k]++; nzp[k] < path[k].getSize(); nzp[k]++) {
      index[k] = path[k].index(nzp[k]);
      down = path[k].down(nzp[k]);
      MEDDLY_DCASSERT(down);
      level_change = k;
      const long ev = path[k].edgeval(nzp[k]).getLong();
      acc_evs[downLevel(k)] = acc_evs[k] + ev;
      if (first(downLevel(k), down)) return true;
    }
  } // for

  return false;
}


bool MEDDLY::evmxd_pluslong::evtrmxd_fixedrow_iter::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);
  MEDDLY_DCASSERT(acc_evs);

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
    if (0==cdown) return false;
    acc_evs[-k] = acc_evs[k];
  } else {
    long ev;
    F->getDownPtr(down, index[k], ev, cdown);
    if (0==cdown) return false;
    acc_evs[-k] = acc_evs[k] + ev;
  }

  //
  // Ok, set up the "column" node below
  k = downLevel(k);
  MEDDLY_DCASSERT(k<0);

  if (isLevelAbove(k, F->getNodeLevel(cdown))) {
    // Skipped level, we can be fast about this.
    acc_evs[downLevel(k)] = acc_evs[k];
    // first, recurse.
    if (!first(downLevel(k), cdown)) return false;
    // Ok, there is a valid path.
    // Set up this level.
    nzp[k] = 0;
    if (F->isFullyReduced()) {
      path[k].initRedundant(F, k, 0L, cdown, SPARSE_ONLY);
      index[k] = 0;
    } else {
      index[k] = index[upLevel(k)];
      path[k].initIdentity(F, k, index[k], 0L, cdown, SPARSE_ONLY);
    }
    return true;
  }

  // Proper node here.
  // cycle through it and recurse...

  F->unpackNode(path+k, cdown, SPARSE_ONLY);

  for (unsigned z=0; z<path[k].getSize(); z++) {
    const long ev = path[k].edgeval(z).getLong();
    acc_evs[downLevel(k)] = acc_evs[k] + ev;
    if (first(downLevel(k), path[k].down(z))) {
      nzp[k] = z;
      index[k] = path[k].index(z);
      return true;
    }
  }

  return false;
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *         evmxd_pluslong::evtrmxd_fixedcol_iter methods         *
// *                                                                *
// *                                                                *
// ******************************************************************

bool MEDDLY::evmxd_pluslong::evtrmxd_fixedcol_iter
::start(const dd_edge &e, const int* minterm)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
  }

  MEDDLY_DCASSERT(acc_evs);
  long ev;
  e.getEdgeValue(ev);
  acc_evs[maxLevel] = ev;

  for (int k=1; k<=maxLevel; k++) {
    index[-k] = minterm[k];
  }
  return first(maxLevel, e.getNode());
}

bool MEDDLY::evmxd_pluslong::evtrmxd_fixedcol_iter::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);
  MEDDLY_DCASSERT(acc_evs);

  node_handle down = 0;
  // Only try to advance the row, because the column is fixed.
  for (int k=1; k<=maxLevel; k++) {
    for (nzp[k]++; nzp[k] < path[k].getSize(); nzp[k]++) {
      index[k] = path[k].index(nzp[k]);
      down = path[k].down(nzp[k]);
      MEDDLY_DCASSERT(down);
      const long ev = path[k].edgeval(nzp[k]).getLong();
      acc_evs[downLevel(k)] = acc_evs[k] + ev;
      level_change = k;
      if (first(downLevel(k), down)) return true;
    }
  } // for

  return false;
}

bool MEDDLY::evmxd_pluslong::evtrmxd_fixedcol_iter::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);
  MEDDLY_DCASSERT(acc_evs);

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
      acc_evs[downLevel(k)] = acc_evs[k];
      return first(downLevel(k), down);
    }
    long ev;
    int cdown;
    F->getDownPtr(down, index[k], ev, cdown);
    if (0==cdown) return false;
    acc_evs[downLevel(k)] = acc_evs[k] + ev;
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
      acc_evs[kpr] = acc_evs[k];
      // See if there is a valid path below.
      if (!first(downLevel(kpr), down)) return false;
      // There's one below, set up the one at these levels.
      path[k].initRedundant(F, k, 0L, down, SPARSE_ONLY);
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
    long ev;
    int cdown;
    F->getDownPtr(down, index[kpr], ev, cdown);
    if (0==cdown) return false;
    acc_evs[downLevel(kpr)] = acc_evs[kpr] * ev;
    if (!first(kpr, cdown)) return false;
    path[k].initRedundant(F, k, 0L, down, SPARSE_ONLY);
    nzp[k] = 0;
    index[k] = 0;
    return true;
  }

  // Level is not skipped.
  F->unpackNode(path+k, down, SPARSE_ONLY);

  for (unsigned z=0; z<path[k].getSize(); z++) {
    index[k] = path[k].index(z);
    const long ev = path[k].edgeval(z).getLong();
    acc_evs[downLevel(k)] = acc_evs[k] + ev;
    if (first(downLevel(k), path[k].down(z))) {
      nzp[k] = z;
      return true;
    }
  }
  return false;
}

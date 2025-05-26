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

#include "evmdd_timesreal.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                    evmdd_timesreal  methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::evmdd_timesreal::evmdd_timesreal(domain *d, const policies &p)
 : evmdd_forest(d, range_type::REAL, edge_labeling::EVTIMES,
         p)
{
    setFloatEdges();
    setTransparentEdge(0, float(0));
    initializeStorage();
}

MEDDLY::evmdd_timesreal::~evmdd_timesreal()
{ }

#ifdef ALLOW_DEPRECATED_0_17_7

void MEDDLY::evmdd_timesreal::createEdge(float val, dd_edge &e)
{
  createEdgeTempl<OP, float>(val, e);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true, __FILE__, __LINE__);
#endif
}

void MEDDLY::evmdd_timesreal::createEdge(double val, dd_edge &e)
{
  createEdgeTempl<OP, double>(val, e);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true, __FILE__, __LINE__);
#endif
}


void MEDDLY::evmdd_timesreal
::createEdge(const int* const* vlist, const float* terms, int N, dd_edge &e)
{
  binary_operation* unionOp = 0;  // for now
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);

  evmdd_edgemaker<OP, float>
  EM(this, vlist, terms, order, N, getDomain()->getNumVariables(), unionOp);

  float ev;
  node_handle ep;
  EM.createEdge(ev, ep);
  e.set(ev, ep);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true, __FILE__, __LINE__);
#endif
}

void MEDDLY::evmdd_timesreal
::evaluate(const dd_edge &f, const int* vlist, float &term) const
{
  evaluateT<OP, float>(f, vlist, term);
}

#endif

#ifdef ALLOW_DEPRECATED_0_17_9

void MEDDLY::evmdd_timesreal
::createEdgeForVar(int vh, bool vp, const float* terms, dd_edge& a)
{
  createEdgeForVarTempl<OP, float>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true, __FILE__, __LINE__);
#endif
}

#endif

#ifdef VIRTUAL_IO_METHODS

void MEDDLY::evmdd_timesreal::showEdge(output &s, const edge_value &ev,
        node_handle d) const
{
    if (d == 0) {
        s.put("<0, w>");
    } else {
        s.put('<');
        s.put(ev.getFloat());
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

#endif

// ******************************************************************
// *                                                                *
// *                                                                *
// *           evmdd_timesreal::evtrmdd_iterator  methods           *
// *                                                                *
// *                                                                *
// ******************************************************************

#ifdef ALLOW_DEPRECATED_0_17_7

MEDDLY::evmdd_timesreal::evtrmdd_iterator::evtrmdd_iterator(const forest *F)
: iterator(F)
{
  int N = F->getNumVariables();
  acc_evs = new double[N+1];
}

MEDDLY::evmdd_timesreal::evtrmdd_iterator::~evtrmdd_iterator()
{
  delete[] acc_evs;
}

void MEDDLY::evmdd_timesreal::evtrmdd_iterator::getValue(float &tv) const
{
  MEDDLY_DCASSERT(acc_evs);
  tv = acc_evs[0];
}

bool MEDDLY::evmdd_timesreal::evtrmdd_iterator::start(const dd_edge &e)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
  }

  MEDDLY_DCASSERT(acc_evs);
  float ev;
  e.getEdgeValue(ev);
  acc_evs[maxLevel] = ev;

  return first(maxLevel, e.getNode());
}

bool MEDDLY::evmdd_timesreal::evtrmdd_iterator::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(!F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);
  MEDDLY_DCASSERT(acc_evs);

  int k;
  node_handle down = 0;
  for (k=1; k<=maxLevel; k++) {
    nzp[k]++;
    if (nzp[k] < path[k].getSize()) {
      index[k] = path[k].index(nzp[k]);
      down = path[k].down(nzp[k]);
      MEDDLY_DCASSERT(down);
      const float ev = path[k].edgeval(nzp[k]).getFloat();
      acc_evs[k-1] = acc_evs[k] * ev;
      break;
    }
  }
  level_change = k;
  if (k>maxLevel) {
    return false;
  }

  return first(k-1, down);
}

bool MEDDLY::evmdd_timesreal::evtrmdd_iterator::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(!F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);
  MEDDLY_DCASSERT(acc_evs);

  if (0==down) return false;

  for ( ; k; k--) {
    MEDDLY_DCASSERT(down);
    int kdn = F->getNodeLevel(down);
    MEDDLY_DCASSERT(kdn <= k);
    if (kdn < k)  {
      path[k].initRedundant(F, k, 1.0f, down, SPARSE_ONLY);
    } else {
      F->unpackNode(path+k, down, SPARSE_ONLY);
    }
    nzp[k] = 0;
    index[k] = path[k].index(0);
    down = path[k].down(0);
    const float ev = path[k].edgeval(0).getFloat();
    acc_evs[k-1] = acc_evs[k] * ev;
  }
  // save the terminal value
  index[0] = down;
  return true;
}

#endif

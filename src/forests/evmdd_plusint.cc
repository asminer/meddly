
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

#include "evmdd_plusint.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                     evmdd_plusint  methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::evmdd_plusint
 ::evmdd_plusint(int dsl, domain *d, const policies &p, bool index_set)
 : evmdd_forest(dsl, d, INTEGER, index_set ? INDEX_SET : EVPLUS, p)
{
  setEdgeSize(sizeof(node_handle), true);
  if (index_set) setUnhashedSize(sizeof(node_handle));
  initializeForest();
}

MEDDLY::evmdd_plusint::~evmdd_plusint()
{ }

void MEDDLY::evmdd_plusint::createEdge(int val, dd_edge &e)
{
  createEdgeTempl<OP, int>(val, e);
}

void MEDDLY::evmdd_plusint
::createEdge(const int* const* vlist, const int* terms, int N, dd_edge &e)
{
  // binary_operation* unionOp = getOperation(PLUS, this, this, this);
  binary_operation* unionOp = 0;  // for now
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);

  evmdd_edgemaker<OP, int>
  EM(this, vlist, terms, order, N, getDomain()->getNumVariables(), unionOp);

  int ev;
  node_handle ep;
  EM.createEdge(ev, ep);
  e.set(ep, ev);
}

void MEDDLY::evmdd_plusint
::createEdgeForVar(int vh, bool vp, const int* terms, dd_edge& a)
{
  createEdgeForVarTempl<OP, int>(vh, vp, terms, a);
}

void MEDDLY::evmdd_plusint
::evaluate(const dd_edge &f, const int* vlist, int &term) const
{
  evaluateT<OP, int>(f, vlist, term);
}

bool MEDDLY::evmdd_plusint
::areEdgeValuesEqual(const void* eva, const void* evb) const
{
  int val1, val2;
  OP::readValue(eva, val1);
  OP::readValue(evb, val2);
  return (val1 == val2);
}

bool MEDDLY::evmdd_plusint::isRedundant(const node_builder &nb) const
{
  return isRedundantTempl<OP>(nb);
}

bool MEDDLY::evmdd_plusint::isIdentityEdge(const node_builder &nb, int i) const
{
  return isIdentityEdgeTempl<OP>(nb, i); 
}

bool MEDDLY::evmdd_plusint::isRedundant(const unpacked_node &nb) const
{
  return isRedundantTempl<OP>(nb);
}

bool MEDDLY::evmdd_plusint::isIdentityEdge(const unpacked_node &nb, int i) const
{
  return isIdentityEdgeTempl<OP>(nb, i); 
}




void MEDDLY::evmdd_plusint::normalize(node_builder &nb, int& ev) const
{
  int minindex = -1;
  int stop = nb.isSparse() ? nb.getNNZs() : nb.getSize();
  for (int i=0; i<stop; i++) {
    if (0==nb.d(i)) continue;
    if ((minindex < 0) || (nb.ei(i) < nb.ei(minindex))) {
      minindex = i;
    }
  }
  if (minindex < 0) return; // this node will eventually be reduced to "0".
  ev = nb.ei(minindex);
  for (int i=0; i<stop; i++) {
    if (0==nb.d(i)) continue;
    int temp;
    nb.getEdge(i, temp);
    temp -= ev;
    nb.setEdge(i, temp);
  }
}

void MEDDLY::evmdd_plusint::showEdgeValue(output &s, const void* edge) const
{
  OP::show(s, edge);
}

void MEDDLY::evmdd_plusint::writeEdgeValue(output &s, const void* edge) const
{
  OP::write(s, edge);
}

void MEDDLY::evmdd_plusint::readEdgeValue(input &s, void* edge)
{
  OP::read(s, edge);
}

void MEDDLY::evmdd_plusint::showUnhashedHeader(output &s, const void* uh) const
{
  s.put(" card: ");
  s.put(long( ((const node_handle*)uh)[0]  ));
  // fprintf(s, " card: %d", ((const node_handle*)uh)[0]);
}

void MEDDLY::evmdd_plusint::writeUnhashedHeader(output &s, const void* uh) const
{
  s.put("\t ");
  s.put(long( ((const node_handle*)uh)[0]  ));
  s.put('\n');
  // th_fprintf(s, "\t %d\n", ((const node_handle*)uh)[0]);
}

void MEDDLY::evmdd_plusint::readUnhashedHeader(input &s, node_builder &nb) const
{
  ((node_handle*)nb.UHptr())[0] = s.get_integer();
  // th_fscanf(1, s, "%d", (node_handle*)nb.UHptr());
}

const char* MEDDLY::evmdd_plusint::codeChars() const
{
  return "dd_epvi";
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *            evmdd_plusint::evpimdd_iterator  methods            *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::evmdd_plusint::evpimdd_iterator::evpimdd_iterator(const expert_forest *F)
: iterator(F)
{
  int N = F->getNumVariables();
  acc_evs = new long[N+1];
}

MEDDLY::evmdd_plusint::evpimdd_iterator::~evpimdd_iterator()
{
  delete[] acc_evs;
}

void MEDDLY::evmdd_plusint::evpimdd_iterator::getValue(int &tv) const
{
  MEDDLY_DCASSERT(acc_evs);
  tv = acc_evs[0];
}

bool MEDDLY::evmdd_plusint::evpimdd_iterator::start(const dd_edge &e)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH);
  }

  int ev;
  e.getEdgeValue(ev);
  acc_evs[maxLevel] = ev;

  return first(maxLevel, e.getNode());
}

bool MEDDLY::evmdd_plusint::evpimdd_iterator::next()
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
    if (nzp[k] < path[k].getNNZs()) {
      index[k] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      int ev;
      path[k].getEdge(nzp[k], ev);
      acc_evs[k-1] = acc_evs[k] + ev;
      break;
    }
  }
  level_change = k;
  if (k>maxLevel) {
    return false;
  }

  return first(k-1, down);
}

bool MEDDLY::evmdd_plusint::evpimdd_iterator::first(int k, node_handle down)
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
    if (kdn < k)  path[k].initRedundant(F, k, 0, down, false);
    else          path[k].initFromNode(F, down, false);
    nzp[k] = 0;
    index[k] = path[k].i(0);
    down = path[k].d(0);
    int ev;
    path[k].getEdge(0, ev);
    acc_evs[k-1] = acc_evs[k] + ev;
  }
  // save the terminal value
  index[0] = down;
  return true;
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                    evmdd_index_set  methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::evmdd_index_set::evmdd_index_set(int dsl, domain *d, const policies &p)
 : evmdd_plusint(dsl, d, p, true)
{ }

MEDDLY::evmdd_index_set::~evmdd_index_set()
{ }

void MEDDLY::evmdd_index_set::getElement(const dd_edge &a, int index, int* e)
{
  if (e == 0) throw error(error::INVALID_VARIABLE);
  if (index < 0) {
    e[0] = 0;
    return;
  }
  int p = a.getNode();
  unpacked_node* R = unpacked_node::useUnpackedNode();
  for (int k=getNumVariables(); k; k--) {
    MEDDLY_DCASSERT(index >= 0);
    if (p <= 0) {
      e[k] = 0;
      continue;
    }
    R->initFromNode(this, p, false);
    MEDDLY_DCASSERT(R->getLevel() <= k);
    if (R->getLevel() < k) {
      e[k] = 0;
      continue; 
    }
    // Find largest i such that edge value i is not greater than index
    e[k] = 0;
    p = 0;
    for (int z=R->getNNZs()-1; z>=0; z--) {
      if (index < R->ei(z)) continue;
      e[k] = R->i(z);
      p = R->d(z);
      index -= R->ei(z);
      break;
    } // for z
  } // for k
  if (index)    e[0] = 0;
  else          e[0] = -p;
  unpacked_node::recycle(R);
}


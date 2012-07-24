
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


#include "../unique_table.h"
#include "evmdd.h"


// ********************************** EVMDDs ********************************** 


MEDDLY::evmdd_forest
::evmdd_forest(int dsl, domain *d, range_type t, edge_labeling el, 
  const policies &p)
: expert_forest(dsl, d, false, t, el, p)
{
}


MEDDLY::evmdd_forest::~evmdd_forest()
{ }


void MEDDLY::evmdd_forest::showTerminal(FILE* s, int tnode) const
{
  fprintf(s, "t%d", -tnode);
}



// ********************************* EV+MDDs ********************************** 

MEDDLY::evp_mdd_int::evp_mdd_int(int dsl, domain *d, const policies &p)
: evmdd_forest(dsl, d, forest::INTEGER, forest::EVPLUS, p)
{ 
  initializeForest();
}

MEDDLY::evp_mdd_int::~evp_mdd_int()
{ }


void MEDDLY::evp_mdd_int::getElement(const dd_edge& a, int index, int* e)
{
  if (e == 0) throw error(error::INVALID_VARIABLE);
  if (index < 0) {
    e[0] = 0;
    return;
  }
  int p = a.getNode();
  node_reader* R = node_reader::useReader();
  for (int k=getNumVariables(); k; k--) {
    MEDDLY_DCASSERT(index >= 0);
    if (p <= 0) {
      e[k] = 0;
      continue;
    }
    initNodeReader(*R, p, false);
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
  node_reader::recycle(R);
}

void
MEDDLY::evp_mdd_int::
createEdge(const int* const* vlist, const int* terms, int N, dd_edge &e)
{
  return createEdgeInternal(vlist, terms, N, e);
}


void
MEDDLY::evp_mdd_int::
createEdge(int val, dd_edge &e)
{
  return createEdgeInternal(val, e);
}


void
MEDDLY::evp_mdd_int::
evaluate(const dd_edge &f, const int* vlist, int &term) const
{
  return evaluateInternal(f, vlist, term);
}

bool MEDDLY::evp_mdd_int::areDuplicates(int node, const node_builder &nb) const
{
  return areDupsInternal(node, nb);
}

bool MEDDLY::evp_mdd_int::areDuplicates(int node, const node_reader &nr) const
{
  return areDupsInternal(node, nr);
}

void MEDDLY::evp_mdd_int::normalize(node_builder &nb, int &ev) const
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
    nb.ei(i) -= ev;
  }
}

void MEDDLY::evp_mdd_int::showEdgeValue(FILE* s, const void* edge, int i) const
{
  fprintf(s, "%d", ((const int*)edge)[i]);
}

void MEDDLY::evp_mdd_int::showUnhashedHeader(FILE* s, const int* uh) const
{
  fprintf(s, " card: %d", uh[0]);
}

bool MEDDLY::evp_mdd_int::isRedundant(const node_builder &nb) const
{
  if (isQuasiReduced()) return false;
  int common = nb.d(0);
  for (int i=1; i<nb.rawSize(); i++) {
    if (nb.d(i) != common)  return false;
    if (nb.ei(i) != 0)      return false;
  }
  return true;
}

bool MEDDLY::evp_mdd_int::isIdentityEdge(const node_builder &nb, int i) const
{
  return false;
}


// ********************************* EV*MDDs ********************************** 

MEDDLY::evt_mdd_real::evt_mdd_real(int dsl, domain *d, const policies &p)
: evmdd_forest(dsl, d, forest::REAL, forest::EVTIMES, p)
{ 
  initializeForest();
}

MEDDLY::evt_mdd_real::~evt_mdd_real()
{ }

void
MEDDLY::evt_mdd_real::
createEdge(const int* const* vlist, const float* terms, int N, dd_edge &e)
{
  return createEdgeInternal(vlist, terms, N, e);
}


void
MEDDLY::evt_mdd_real::
createEdge(float val, dd_edge &e)
{
  return createEdgeInternal(val, e);
}


void
MEDDLY::evt_mdd_real::
evaluate(const dd_edge &f, const int* vlist, float &term) const
{
  return evaluateInternal(f, vlist, term);
}


void MEDDLY::evt_mdd_real::normalize(node_builder &nb, float &ev) const
{
  int maxindex = -1;
  int stop = nb.isSparse() ? nb.getNNZs() : nb.getSize();
  for (int i=0; i<stop; i++) {
    if (0==nb.d(i)) continue;
    if ((maxindex < 0) || (nb.ef(i) < nb.ef(maxindex))) {
      maxindex = i;
    }
  }
  if (maxindex < 0) return; // this node will eventually be reduced to "0".
  ev = nb.ef(maxindex);
  MEDDLY_DCASSERT(ev);
  for (int i=0; i<stop; i++) {
    if (0==nb.d(i)) continue;
    nb.ef(i) /= ev;
  }
}


void MEDDLY::evt_mdd_real::showEdgeValue(FILE* s, const void* edge, int i) const
{
  fprintf(s, "%f", ((const float*)edge)[i]);
}

bool MEDDLY::evt_mdd_real::areDuplicates(int node, const node_builder &nb) const
{
  return areDupsInternal(node, nb);
}

bool MEDDLY::evt_mdd_real::areDuplicates(int node, const node_reader &nr) const
{
  return areDupsInternal(node, nr);
}

bool MEDDLY::evt_mdd_real::isRedundant(const node_builder &nb) const
{
  if (isQuasiReduced()) return false;
  int common = nb.d(0);
  for (int i=1; i<nb.rawSize(); i++) {
    if (nb.d(i) != common)  return false;
    if (nb.ef(i) != 1)      return false;
  }
  return true;
}

bool MEDDLY::evt_mdd_real::isIdentityEdge(const node_builder &nb, int i) const
{
  return false;
}



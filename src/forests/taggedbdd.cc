
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

#include "taggedbdd.h"

#include "../unique_table.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      taggedbdd  methods                        *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::taggedbdd
 ::taggedbdd(int dsl, domain *d, const policies &p, int* level_reduction_rule)
 : evmdd_forest(dsl, d, BOOLEAN, TAGGED, p, level_reduction_rule)
{
  // The reduction can only be edge-specific
  deflt.reduction = policies::EDGE_SPECIFIC;

  setEdgeSize(sizeof(long), true);
  initializeForest();
}

MEDDLY::taggedbdd::~taggedbdd()
{ }

void MEDDLY::taggedbdd::createEdge(long reduction, dd_edge &e)
{
  createEdgeTempl<OP, long>(reduction, e);
}

//void MEDDLY::evmdd_pluslong
//::createEdge(const int* const* vlist, const long* terms, int N, dd_edge &e)
//{
//  // binary_operation* unionOp = getOperation(PLUS, this, this, this);
//  binary_operation* unionOp = 0;  // for now
//  enlargeStatics(N);
//  enlargeVariables(vlist, N, false);
//
//  int num_vars = getNumVariables();
//
//  // Create vlist following the mapping between variable and level
//  int** ordered_vlist = static_cast<int**>(malloc(N * sizeof(int*) + (num_vars + 1) * N * sizeof(int)));
//  if (ordered_vlist == nullptr) {
//	  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
//  }
//
//  ordered_vlist[0] = reinterpret_cast<int*>(&ordered_vlist[N]);
//  for (int i = 1; i < N; i++) {
//	  ordered_vlist[i] = (ordered_vlist[i - 1] + num_vars + 1);
//  }
//  for (int i = 0; i <= num_vars; i++) {
//	  int level = getLevelByVar(i);
//	  for (int j = 0; j < N; j++) {
//		  ordered_vlist[j][level] = vlist[j][i];
//	  }
//  }
//
//  long* terms_long = static_cast<long*>(malloc(N * sizeof(long)));
//  for (int i = 0; i < N; i++) {
//    terms_long[i] = terms[i];
//  }
//
//  evmdd_edgemaker<OP, long>
//  EM(this, ordered_vlist, terms_long, order, N, num_vars, unionOp);
//
//  long ev = Inf<long>();
//  node_handle ep = 0;
//  EM.createEdge(ev, ep);
//  e.set(ep, ev);
//
//  free(ordered_vlist);
//  free(terms_long);
//}

//void MEDDLY::evmdd_pluslong
//::createEdgeForVar(int vh, bool vp, const long* terms, dd_edge& a)
//{
//  int sz = this->getVariableSize(vh);
//
//  long* terms_long = new long[sz];
//  for (int i = 0; i < sz; i++) {
//    terms_long[i] = terms[i];
//  }
//
//  createEdgeForVarTempl<OP, long>(vh, vp, terms_long, a);
//
//  delete[] terms_long;
//}


void MEDDLY::taggedbdd
::evaluate(const dd_edge &f, const int* vlist, long &term) const
{
  evaluateT<OP, long>(f, vlist, term);
}

bool MEDDLY::taggedbdd
::isTransparentEdge(node_handle ep, const void* ev) const
{
  MEDDLY_DCASSERT(ep != 0 || OP::isTransparentEdge(ev));
  if (ep != getTransparentNode()) return false;
  return OP::isTransparentEdge(ev);
}

void MEDDLY::taggedbdd
::getTransparentEdge(node_handle &ep, void* ev) const
{
  ep = 0;
  OP::setEdge(ev, OP::getRedundantEdge());
}

bool MEDDLY::taggedbdd
::areEdgeValuesEqual(const void* eva, const void* evb) const
{
  long val1, val2;
  OP::readValue(eva, val1);
  OP::readValue(evb, val2);
  return (val1 == val2);
}

bool MEDDLY::taggedbdd::isRedundant(const unpacked_node &nb) const
{
  if (nb.isSparse()) {
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
  }
  else {
    int common_child = nb.d(0);
    long common_child_tag;
    nb.getEdge(0, common_child_tag);
    for (int i = 1; i < nb.getSize(); i++) {
      if (nb.d(i) != common_child || nb.ei(i) != common_child_tag) {
        return false;
      }
    }
    return true;
  }
}

bool MEDDLY::taggedbdd::isIdentityEdge(const unpacked_node &nb, int i) const
{
  return isIdentityEdgeTempl<OP>(nb, i);
}

bool MEDDLY::taggedbdd::isZeroSuppressed(const unpacked_node &nb) const
{
  if (nb.isSparse()) {
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
  }
  else {
    for (int i = 1; i < nb.getSize(); i++) {
      if (nb.d(i) != getTransparentNode() || nb.ei(i) != 0L) {
        return false;
      }
    }
    return true;
  }
}

//bool MEDDLY::cbdd::isOneSuppressed(const unpacked_node &nb) const
//{
//  if (getLevelSize(nb.getLevel()) < 2) {
//    return false;
//  }
//
//  if (nb.isSparse()) {
//    if (nb.getNNZs() > 1) {
//      return false;
//    }
//    if (nb.getNNZs() == 1) {
//      if (nb.i(0) != 1) {
//        return false;
//      }
//    }
//    return true;
//  }
//  else {
//    for (int i = 0; i < nb.getSize(); i++) {
//      if (i == 1) {
//        continue;
//      }
//
//      if (nb.d(i) != getTransparentNode()) {
//        return false;
//      }
//    }
//    return true;
//  }
//}

void MEDDLY::taggedbdd::normalize(unpacked_node &nb, long& ev) const
{
//  assert(nb.getSize() == 2);
//
//  if (nb.isSparse()) nb.sort();
//
//  // get sparse, truncated full sizes and check
//  // for redundant / identity reductions.
//  int nnz;
//  if (nb.isSparse()) {
//    // Reductions for sparse nodes
//    nnz = nb.getNNZs();
//#ifdef DEVELOPMENT_CODE
//    for (int z=0; z<nnz; z++) {
//      MEDDLY_DCASSERT(nb.d(z)!=getTransparentNode());
//    } // for z
//#endif
//
//    // Check for identity nodes
//    if (0==nnz) {
//        return nb.d(0);
//      }
//    }
//
//    // Check for redundant nodes
//    if (nnz == getLevelSize(nb.getLevel()) && !isExtensibleLevel(nb.getLevel())) {
//      if (isRedundant(nb)) {
//        // unlink downward pointers, except the one we're returning.
//        for (int i = 1; i<nnz; i++)  unlinkNode(nb.d(i));
//#ifdef DEBUG_CREATE_REDUCED
//        printf("Redundant node ");
//        showNode(stdout, nb.d(0), SHOW_DETAILS | SHOW_INDEX);
//        printf("\n");
//#endif
//        return nb.d(0);
//      }
//    }
//
//  }
}

void MEDDLY::taggedbdd::showEdgeValue(output &s, const void* edge) const
{
  OP::show(s, edge);
}

void MEDDLY::taggedbdd::writeEdgeValue(output &s, const void* edge) const
{
  OP::write(s, edge);
}

void MEDDLY::taggedbdd::readEdgeValue(input &s, void* edge)
{
  OP::read(s, edge);
}

void MEDDLY::taggedbdd::showUnhashedHeader(output &s, const void* uh) const
{
  s.put(" card: ");
  s.put(static_cast<const long*>(uh)[0]);
  // fprintf(s, " card: %d", ((const node_handle*)uh)[0]);
}

void MEDDLY::taggedbdd::writeUnhashedHeader(output &s, const void* uh) const
{
  s.put("\t ");
  s.put(static_cast<const long*>(uh)[0]);
  s.put('\n');
  // th_fprintf(s, "\t %d\n", ((const node_handle*)uh)[0]);
}

void MEDDLY::taggedbdd::readUnhashedHeader(input &s, unpacked_node &nb) const
{
  static_cast<long*>(nb.UHdata())[0] = s.get_integer();
  // th_fscanf(1, s, "%d", (node_handle*)nb.UHptr());
}

const char* MEDDLY::taggedbdd::codeChars() const
{
  return "dd_epvi";
}

void MEDDLY::taggedbdd::createReducedHelper(int in, unpacked_node& nb, long& tag, node_handle& node)
{
#ifdef DEVELOPMENT_CODE
  validateDownPointers(nb);
#endif

  MEDDLY_DCASSERT(isEdgeSpecificReduced());
  MEDDLY_DCASSERT(getLevelSize(nb.getLevel()) == 2);

  // check is the node is written in order,
  // if not rearrange it in ascending order of indices.
  if (nb.isSparse()) nb.sort();

  tag = nb.getLevel();

  //if (nb.isExtensible()) return createReducedExtensibleNodeHelper(in, nb);

  if (nb.isSparse()) {
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
  }
  else {
    // Reductions for full nodes
    MEDDLY_DCASSERT(nb->getSize() == getLevelSize(nb.getLevel()));
    // Check for redundant nodes
    if (isRedundant(nb)) {
      // Rule 1
      nb.getEdge(0, tag);
      for (int i = 1; i < nb.getSize(); i++) {
        unlinkNode(nb.d(i));
      }
      node = nb.d(0);
      return;
    }
    // Check for zero-suppressed nodes
    else if (isZeroSuppressed(nb)) {
      node_handle child = nb.d(0);
      long child_tag;
      nb.getEdge(0, child_tag);
      if (child_tag == nb.getLevel() - 1) {
        // Rule 2a
        node = child;
        return;
      }
      else {
        // Rule 2b
        MEDDLY_DCASSERT(child_tag < nb->getLevel() - 1);

        int nextk = nb.getLevel() - 1;
        nb.initFull(this, nextk, getLevelSize(nextk));
        // Create an extra redundant node
        for (int i = 0; i < nb.getSize(); i++) {
          nb.d_ref(i) = linkNode(child);
          nb.setEdge(i, child_tag);
        }
        unlinkNode(child);
      }
    }
  }

  nb.computeHash();

  // check for duplicates in unique table
  node_handle q = unique->find(nb, getVarByLevel(nb.getLevel()));
  if (q) {
    // unlink all downward pointers
    int rawsize = nb.isSparse() ? nb.getNNZs() : nb.getSize();
    for (int i = 0; i < rawsize; i++)  unlinkNode(nb.d(i));
    node = linkNode(q);
    return;
  }

  //
  // Not eliminated by reduction rule.
  // Not a duplicate.
  //
  // We need to create a new node for this.

  // NOW is the best time to run the garbage collector, if necessary.
//#ifndef GC_OFF
//  if (isTimeToGc()) garbageCollect();
//#endif

  // Grab a new node
  node_handle p = nodeHeaders.getFreeNodeHandle();
  nodeHeaders.setNodeLevel(p, nb.getLevel());
  MEDDLY_DCASSERT(0 == nodeHeaders.getNodeCacheCount(p));
  MEDDLY_DCASSERT(0 == nodeHeaders.getIncomingCount(p));

  stats.incActive(1);
  if (theLogger && theLogger->recordingNodeCounts()) {
    theLogger->addToActiveNodeCount(this, nb.getLevel(), 1);
  }

  // All of the work is in nodeMan now :^)
  nodeHeaders.setNodeAddress(p, nodeMan->makeNode(p, nb, getNodeStorage()));
  node = linkNode(p);

  // add to UT
  unique->add(nb.hash(), p);

#ifdef DEVELOPMENT_CODE
  unpacked_node* key = unpacked_node::newFromNode(this, p, false);
  key->computeHash();
  MEDDLY_DCASSERT(key->hash() == nb.hash());
  node_handle f = unique->find(*key, getVarByLevel(key->getLevel()));
  MEDDLY_DCASSERT(f == p);
  unpacked_node::recycle(key);
#endif
#ifdef DEBUG_CREATE_REDUCED
  printf("Created node ");
  showNode(stdout, p, SHOW_DETAILS | SHOW_INDEX);
  printf("\n");
#endif
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *               esrbdd::evpimdd_iterator  methods                *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::taggedbdd::evpimdd_iterator::evpimdd_iterator(const expert_forest *F)
: iterator(F)
{
  int N = F->getNumVariables();
  acc_evs = new long[N+1];
}

MEDDLY::taggedbdd::evpimdd_iterator::~evpimdd_iterator()
{
  delete[] acc_evs;
}

void MEDDLY::taggedbdd::evpimdd_iterator::getValue(long &tv) const
{
  MEDDLY_DCASSERT(acc_evs);
  tv = acc_evs[0];
}

bool MEDDLY::taggedbdd::evpimdd_iterator::start(const dd_edge &e)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
  }

  long ev = Inf<long>();
  e.getEdgeValue(ev);
  acc_evs[maxLevel] = ev;

  return first(maxLevel, e.getNode());
}

bool MEDDLY::taggedbdd::evpimdd_iterator::next()
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
      long ev = Inf<long>();
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

bool MEDDLY::taggedbdd::evpimdd_iterator::first(int k, node_handle down)
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
    long ev = Inf<long>();
    path[k].getEdge(0, ev);
    acc_evs[k-1] = acc_evs[k] + ev;
  }
  // save the terminal value
  index[0] = down;
  return true;
}

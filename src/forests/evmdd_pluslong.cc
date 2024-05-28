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

#include "evmdd_pluslong.h"

#include "../unique_table.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                     evmdd_pluslong  methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::evmdd_pluslong
 ::evmdd_pluslong(domain *d, const policies &p, int* level_reduction_rule, bool index_set)
 : evmdd_forest(d, range_type::INTEGER,
         index_set ? edge_labeling::INDEX_SET : edge_labeling::EVPLUS,
         p, level_reduction_rule)
{
    setLongEdges();
    setTransparentEdge(0, long(0));
    if (index_set) setUnhashedSize(sizeof(long));
    initializeStorage();
}

MEDDLY::evmdd_pluslong::~evmdd_pluslong()
{ }

void MEDDLY::evmdd_pluslong::createEdge(long val, dd_edge &e)
{
  createEdgeTempl<OP, long>(val, e);
}

void MEDDLY::evmdd_pluslong
::createEdge(const int* const* vlist, const long* terms, int N, dd_edge &e)
{
  binary_operation* unionOp = 0;  // for now
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);

  int num_vars = getNumVariables();

  // Create vlist following the mapping between variable and level
  int** ordered_vlist = static_cast<int**>(malloc(N * sizeof(int*) + (num_vars + 1) * N * sizeof(int)));
  if (ordered_vlist == nullptr) {
	  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }

  ordered_vlist[0] = reinterpret_cast<int*>(&ordered_vlist[N]);
  for (int i = 1; i < N; i++) {
	  ordered_vlist[i] = (ordered_vlist[i - 1] + num_vars + 1);
  }
  for (int i = 0; i <= num_vars; i++) {
	  int level = getLevelByVar(i);
	  for (int j = 0; j < N; j++) {
		  ordered_vlist[j][level] = vlist[j][i];
	  }
  }

  long* terms_long = static_cast<long*>(malloc(N * sizeof(long)));
  for (int i = 0; i < N; i++) {
    terms_long[i] = terms[i];
  }

  evmdd_edgemaker<OP, long>
  EM(this, ordered_vlist, terms_long, order, N, num_vars, unionOp);

  long ev = Inf<long>();
  node_handle ep = 0;
  EM.createEdge(ev, ep);
  e.set(ep, ev);

  free(ordered_vlist);
  free(terms_long);
}

void MEDDLY::evmdd_pluslong
::createEdgeForVar(int vh, bool vp, const long* terms, dd_edge& a)
{
  int sz = this->getVariableSize(vh);

  long* terms_long = new long[sz];
  if (terms) {
    for (int i = 0; i < sz; i++) {
        terms_long[i] = terms[i];
    }
  } else {
    for (int i = 0; i < sz; i++) {
        terms_long[i] = i;
    }
  }

  createEdgeForVarTempl<OP, long>(vh, vp, terms_long, a);

  delete[] terms_long;
}


void MEDDLY::evmdd_pluslong
::evaluate(const dd_edge &f, const int* vlist, long &term) const
{
  evaluateT<OP, long>(f, vlist, term);
}

void MEDDLY::evmdd_pluslong::swapAdjacentVariables(int level)
{
  MEDDLY_DCASSERT(level >= 1);
  MEDDLY_DCASSERT(level < getNumVariables());

  int hvar = getVarByLevel(level + 1);  // The variable at the higher level
  int lvar = getVarByLevel(level);    // The variable at the lower level
  int hsize = getVariableSize(hvar);  // The size of the variable at the higher level
  int lsize = getVariableSize(lvar);  // The size of the variable at the lower level

  int hnum = unique->getNumEntries(hvar); // The number of nodes associated with the variable at the higher level
  node_handle* hnodes = new node_handle[hnum];
  unique->getItems(hvar, hnodes, hnum);

  int lnum = unique->getNumEntries(lvar); // The nubmer of nodes associated with the variable at the lower level
  node_handle* lnodes = new node_handle[lnum];
  unique->getItems(lvar, lnodes, lnum);

  //    printf("Before: Level %d : %d, Level %d : %d\n",
  //            level+1, hnum,
  //            level, lnum);

  int num = 0;
  // Renumber the level of nodes for the variable to be moved down
  for (int i = 0; i < hnum; i++) {
    unpacked_node* nr = newUnpacked(hnodes[i], FULL_ONLY);

    MEDDLY_DCASSERT(nr->getLevel() == level + 1);
    MEDDLY_DCASSERT(nr->getSize() == hsize);

    for (int j = 0; j < hsize; j++) {
      if (isLevelAbove(getNodeLevel(nr->down(j)), level - 1)) {
        // Remove the nodes corresponding to functions that
        // are independent of the variable to be moved up
        hnodes[num++] = hnodes[i];
        break;
      }
    }
    unpacked_node::Recycle(nr);

    setNodeLevel(hnodes[i], level);
  }
  hnum = num;

  // Renumber the level of nodes for the variable to be moved up
  for (int i = 0; i < lnum; i++) {
    setNodeLevel(lnodes[i], level + 1);
  }

  // Update the variable order
  std::const_pointer_cast<variable_order>(var_order)->exchange(hvar, lvar);

  node_handle** children = new node_handle*[hsize];
  long** sum_evs = new long*[hsize];
  for (int i = 0; i < hsize; i++) {
    children[i] = new node_handle[lsize];
    sum_evs[i] = new long[lsize];
  }

  // Process the rest of nodes for the variable to be moved down
  for (int i = 0; i < hnum; i++) {
    unpacked_node* high_nr = newUnpacked(hnodes[i], FULL_ONLY);

    unpacked_node* high_nb = unpacked_node::newFull(this, level + 1, lsize);
    for (int j = 0; j < hsize; j++) {
      long ev1 = high_nr->edgeval(j).getLong();
      MEDDLY_DCASSERT(ev1 >= 0);

      if (isLevelAbove(level, getNodeLevel(high_nr->down(j)))) {
        for (int k = 0; k < lsize; k++) {
          children[j][k] = high_nr->down(j);
          sum_evs[j][k] = ev1;
        }
      }
      else {
        unpacked_node* nr = newUnpacked(high_nr->down(j), FULL_ONLY);

        MEDDLY_DCASSERT(nr->getSize() == lsize);
        for (int k = 0; k < lsize; k++) {
          children[j][k] = nr->down(k);

          long ev2 = nr->edgeval(k).getLong();
          MEDDLY_DCASSERT(ev2 >= 0);

          sum_evs[j][k] = ev1 + ev2;
        }
        unpacked_node::Recycle(nr);
      }
    }

    for (int j = 0; j < lsize; j++) {
      unpacked_node* low_nb = unpacked_node::newFull(this, level, hsize);
      for (int k = 0; k < hsize; k++) {
        low_nb->setFull(k, sum_evs[k][j], linkNode(children[k][j]));
        // low_nb->d_ref(k) = linkNode(children[k][j]);
        // low_nb->setEdge(k, sum_evs[k][j]);
      }
      node_handle node = 0;
      long ev = 0;
      createReducedNode(-1, low_nb, ev, node);
      high_nb->setFull(j, ev, node);
      // high_nb->d_ref(j) = node;
      // high_nb->setEdge(j, ev);
    }

    unpacked_node::Recycle(high_nr);

    // The reduced node of high_nb must be at level+1
    // Assume the reduced node is at level
    // Then high_nodes[i] corresponds to a function that
    // is independent of the variable to be moved up
    // This is a contradiction
    modifyReducedNodeInPlace(high_nb, hnodes[i]);
  }

  for (int i = 0; i < hsize; i++) {
    delete[] children[i];
    delete[] sum_evs[i];
  }
  delete[] children;
  delete[] sum_evs;

  delete[] hnodes;
  delete[] lnodes;

  //    printf("After: Level %d : %d, Level %d : %d\n",
  //            level+1, unique->getNumEntries(lvar),
  //            level, unique->getNumEntries(hvar));
  //    printf("#Node: %d\n", getCurrentNumNodes());
}

void MEDDLY::evmdd_pluslong::showEdge(output &s, const edge_value &ev,
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

const char* MEDDLY::evmdd_pluslong::codeChars() const
{
  return "dd_epvi";
}

#endif

// ******************************************************************
// *                                                                *
// *                                                                *
// *            evmdd_pluslong::evpimdd_iterator  methods            *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::evmdd_pluslong::evpimdd_iterator::evpimdd_iterator(const forest *F)
: iterator(F)
{
  int N = F->getNumVariables();
  acc_evs = new long[N+1];
}

MEDDLY::evmdd_pluslong::evpimdd_iterator::~evpimdd_iterator()
{
  delete[] acc_evs;
}

void MEDDLY::evmdd_pluslong::evpimdd_iterator::getValue(long &tv) const
{
  MEDDLY_DCASSERT(acc_evs);
  tv = acc_evs[0];
}

bool MEDDLY::evmdd_pluslong::evpimdd_iterator::start(const dd_edge &e)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
  }

  long ev = Inf<long>();
  e.getEdgeValue(ev);
  acc_evs[maxLevel] = ev;

  return first(maxLevel, e.getNode());
}

bool MEDDLY::evmdd_pluslong::evpimdd_iterator::next()
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
      const long ev = path[k].edgeval(nzp[k]).getLong();
      // long ev = Inf<long>();
      // path[k].getEdge(nzp[k], ev);
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

bool MEDDLY::evmdd_pluslong::evpimdd_iterator::first(int k, node_handle down)
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
        edge_value ev(0L);
        path[k].initRedundant(F, k, ev, down, SPARSE_ONLY);
    }
    else {
        F->unpackNode(path+k, down, SPARSE_ONLY);
    }
    nzp[k] = 0;
    index[k] = path[k].index(0);
    down = path[k].down(0);
    // long ev = Inf<long>();
    // path[k].getEdge(0, ev);
    const long ev = path[k].edgeval(0).getLong();
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

MEDDLY::evmdd_index_set_long::evmdd_index_set_long(domain *_d, const policies &p, int* level_reduction_rule)
 : evmdd_pluslong(_d, p, level_reduction_rule, true)
{ }

MEDDLY::evmdd_index_set_long::~evmdd_index_set_long()
{ }

void MEDDLY::evmdd_index_set_long::getElement(const dd_edge &a, int index, int* e)
{
  // TODO: Redesign interface
  getElement(a, (long)index, e);
}

void MEDDLY::evmdd_index_set_long::getElement(const dd_edge &a, long index, int* e)
{
  if (e == nullptr) {
    throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);
  }
  if (index < 0) {
    e[0] = 0;
    return;
  }
  int p = a.getNode();
  unpacked_node* R = unpacked_node::New(this);
  for (int k = getNumVariables(); k > 0; k--) {
	int var = getVarByLevel(k);

    MEDDLY_DCASSERT(index >= 0);
    if (p <= 0) {
      e[var] = 0;
      continue;
    }
    unpackNode(R, p, SPARSE_ONLY);
    MEDDLY_DCASSERT(R->getLevel() <= k);
    if (R->getLevel() < k) {
      e[var] = 0;
      continue;
    }
    // Find largest i such that edge value i is not greater than index
    e[var] = 0;
    p = 0;
    for (int z = R->getSize() - 1; z >= 0; z--) {
      if (index < R->edgeval(z).getLong()) {
        continue;
      }
      e[var] = R->index(z);
      p = R->down(z);
      index -= R->edgeval(z).getLong();
      break;
    } // for z
  } // for k
  e[0] = (index > 0 ? 0 : -p);
  unpacked_node::Recycle(R);
}

void MEDDLY::evmdd_index_set_long::showHeaderInfo(output &s,
        const unpacked_node &nb) const
{

  s.put(" card: ");
  s.put(static_cast<const long*>(nb.UHptr())[0]);
  // fprintf(s, " card: %d", ((const node_handle*)uh)[0]);
}

void MEDDLY::evmdd_index_set_long::writeHeaderInfo(output &s,
        const unpacked_node &nb) const
{
  s.put('\t');
  s.put(static_cast<const long*>(nb.UHptr())[0]);
  s.put('\n');
  // th_fprintf(s, "\t %d\n", ((const node_handle*)uh)[0]);
}

void MEDDLY::evmdd_index_set_long::readHeaderInfo(input &s,
        unpacked_node &nb) const
{
    long card = s.get_integer();
    nb.setUHdata(&card);
#ifdef DEBUG_READ_DD
    std::cerr << "    got cardinality " << card << "\n";
#endif
}


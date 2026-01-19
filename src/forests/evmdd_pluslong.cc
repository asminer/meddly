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
 ::evmdd_pluslong(domain *d, const policies &p, bool index_set)
 : evmdd_forest(d, range_type::INTEGER,
         index_set ? edge_labeling::INDEX_SET : edge_labeling::EVPLUS, p)
{
}

MEDDLY::evmdd_pluslong::~evmdd_pluslong()
{ }

#ifdef ALLOW_DEPRECATED_0_18_0
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
#endif


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
    // unpacked_node* nr = newUnpacked(hnodes[i], FULL_ONLY);
    unpacked_node* nr = unpacked_node::newFromNode(this, hnodes[i], FULL_ONLY);

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
    unpacked_node* high_nr =
        unpacked_node::newFromNode(this, hnodes[i], FULL_ONLY);
    // unpacked_node* high_nr = newUnpacked(hnodes[i], FULL_ONLY);

    unpacked_node* high_nb =
        unpacked_node::newWritable(this, level + 1, lsize, FULL_ONLY);
    for (int j = 0; j < hsize; j++) {
      long ev1 = long(high_nr->edgeval(j));
      MEDDLY_DCASSERT(ev1 >= 0);

      if (isLevelAbove(level, getNodeLevel(high_nr->down(j)))) {
        for (int k = 0; k < lsize; k++) {
          children[j][k] = high_nr->down(j);
          sum_evs[j][k] = ev1;
        }
      }
      else {
        unpacked_node* nr =
            unpacked_node::newFromNode(this, high_nr->down(j), FULL_ONLY);
        // unpacked_node* nr = newUnpacked(high_nr->down(j), FULL_ONLY);

        MEDDLY_DCASSERT(nr->getSize() == lsize);
        for (int k = 0; k < lsize; k++) {
          children[j][k] = nr->down(k);

          long ev2 = long(nr->edgeval(k));
          MEDDLY_DCASSERT(ev2 >= 0);

          sum_evs[j][k] = ev1 + ev2;
        }
        unpacked_node::Recycle(nr);
      }
    }

    for (int j = 0; j < lsize; j++) {
      unpacked_node* low_nb =
          unpacked_node::newWritable(this, level, hsize, FULL_ONLY);
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


// ******************************************************************
// *                                                                *
// *                                                                *
// *                    evmdd_index_set  methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::evmdd_index_set_long::evmdd_index_set_long(domain *_d, const policies &p)
 : evmdd_pluslong(_d, p, true)
{ }

MEDDLY::evmdd_index_set_long::~evmdd_index_set_long()
{ }


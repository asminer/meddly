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


#include "mt.h"
#include "../unique_table.h"

#define ENABLE_CACHE_COUNTING 0
#define ENABLE_IN_COUNTING 0


// ******************************************************************
// *                                                                *
// *                       mt_forest  methods                       *
// *                                                                *
// ******************************************************************

#ifdef ALLOW_DEPRECATED_0_17_7
int* MEDDLY::mt_forest::order;
int  MEDDLY::mt_forest::order_size;
#endif

MEDDLY::mt_forest::mt_forest(domain *d, bool rel,
  range_type t, const policies &p)
: forest(d, rel, t, edge_labeling::MULTI_TERMINAL, p)
{
}

#ifdef ALLOW_DEPRECATED_0_18_0

MEDDLY::node_handle MEDDLY::mt_forest::_makeNodeAtLevel(int k, node_handle d)
{
  int dk = getNodeLevel(d);
  while (dk != k) {
    int up;
    if (dk<0) up = -dk;
    else up = isForRelations() ? -(dk+1) : dk+1;

    // make node at level "up"
    unpacked_node* nb = 0;

    if (isIdentityReduced() && (dk<0)) {
      //
      // make identity reductions below as necessary

      node_handle sd;
      unsigned indx;
      int si;
      if (isSingletonNode(d, indx, sd)) {
        si = int(indx);
      } else {
        si = -1;
        sd = 0;
      }

      int sz = getLevelSize(up);
      MEDDLY_DCASSERT(si < sz);
      const bool add_edge = (si+1 == sz && isExtensibleLevel(up));
      nb = unpacked_node::newWritable(this, up, (add_edge? (1+sz): sz), FULL_ONLY);

      for (int i=0; i<sz; i++) {
        nb->setFull(i, linkNode( (i==si) ? sd : d ));
      }
#ifdef ALLOW_EXTENSIBLE
      if (isExtensibleLevel(up)) {
        nb->markAsExtensible();
        if (add_edge) {
            nb->setFull(sz, linkNode(d));
        }
      }
#endif
    } else {
      //
      // don't worry about identity reductions

#ifdef ALLOW_EXTENSIBLE
      if (isExtensibleLevel(up)) {
        nb = unpacked_node::newFull(this, up, 1);
        nb->setFull(0, linkNode(d));
        nb->markAsExtensible();
      } else {
#endif
        nb = unpacked_node::newWritable(this, up, FULL_ONLY);
        for (int i=0; i<nb->getSize(); i++) {
          nb->setFull(i, linkNode(d));
        }
#ifdef ALLOW_EXTENSIBLE
      }
#endif
    }

    unlinkNode(d);
    d = createReducedNode(-1, nb);
    dk = up;
  } // while
  return d;
}

#endif

#ifdef ALLOW_DEPRECATED_0_17_7

void MEDDLY::mt_forest::initStatics()
{
  order = 0;
  order_size = 0;
}

void MEDDLY::mt_forest::enlargeStatics(int n)
{
  MEDDLY_DCASSERT(n>=0);
  if (n>order_size) {
    order = (int*) realloc(order, n*sizeof(int));
    //terminals = (node_handle*) realloc(terminals, n*sizeof(node_handle));
    //if (0==order || 0==terminals) {
    if (0==order) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    order_size = n;
  }
  for (int i=0; i<n; i++) order[i] = i;
}

void MEDDLY::mt_forest::clearStatics()
{
  free(order);
  order = 0;
  order_size = 0;
}

#endif

// ******************************************************************
// *                                                                *
// *                 mt_forest::mt_iterator methods                 *
// *                                                                *
// ******************************************************************

#ifdef ALLOW_DEPRECATED_0_17_7

MEDDLY::mt_forest::mt_iterator::mt_iterator(const forest *F)
 : iterator(F)
{
}

MEDDLY::mt_forest::mt_iterator::~mt_iterator()
{
}

void MEDDLY::mt_forest::mt_iterator::getValue(int &termVal) const
{
  MEDDLY_DCASSERT(index);
  terminal t(terminal_type::INTEGER, index[0]);
  termVal = t.getInteger();
  // termVal = int_Tencoder::handle2value(index[0]);
}

void MEDDLY::mt_forest::mt_iterator::getValue(float &termVal) const
{
  MEDDLY_DCASSERT(index);
  terminal t(terminal_type::REAL, index[0]);
  termVal = t.getReal();
}

#endif

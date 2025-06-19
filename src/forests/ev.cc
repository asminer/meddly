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


#include "ev.h"
#include "../unique_table.h"

// ******************************************************************
// *                                                                *
// *                       ev_forest  methods                       *
// *                                                                *
// ******************************************************************

#ifdef ALLOW_DEPRECATED_0_17_7

int* MEDDLY::ev_forest::order;
int  MEDDLY::ev_forest::order_size;

#endif

MEDDLY::ev_forest::ev_forest(domain *d, bool rel,
  range_type t, edge_labeling ev, const policies &p)
: forest(d, rel, t, ev, p)
{
  MEDDLY_DCASSERT(ev != edge_labeling::MULTI_TERMINAL);
}

// statics

#ifdef ALLOW_DEPRECATED_0_17_7

void MEDDLY::ev_forest::initStatics()
{
  order = 0;
  order_size = 0;
}

void MEDDLY::ev_forest::enlargeStatics(int n)
{
  MEDDLY_DCASSERT(n>0);
  if (n>order_size) {
    order = (int*) realloc(order, unsigned(n)*sizeof(int));
    //terminals = (node_handle*) realloc(terminals, n*sizeof(node_handle));
    //if (0==order || 0==terminals) {
    if (0==order) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    order_size = n;
  }
  for (int i=0; i<n; i++) order[i] = i;
}

void MEDDLY::ev_forest::clearStatics()
{
  free(order);
  order = 0;
  order_size = 0;
}

#endif

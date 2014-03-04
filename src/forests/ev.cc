
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


#include "ev.h"
#include "../unique_table.h"

// ******************************************************************
// *                                                                *
// *                        ev_cleanup class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::ev_cleanup : public cleanup_procedure {
  public:
    ev_cleanup() : cleanup_procedure() { 
      // nothing
    }
    virtual void execute() {
      ev_forest::clearStatics();
    }
};

// ******************************************************************
// *                                                                *
// *                       ev_forest  methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::ev_cleanup* MEDDLY::ev_forest::the_ev_cleaner = 0;
int* MEDDLY::ev_forest::order = 0;
int  MEDDLY::ev_forest::order_size = 0;

MEDDLY::ev_forest::ev_forest(int dsl, domain *d, bool rel,
  range_type t, edge_labeling ev, const policies &p)
: expert_forest(dsl, d, rel, t, ev, p)
{
  MEDDLY_DCASSERT(ev != MULTI_TERMINAL);
}

void MEDDLY::ev_forest::showTerminal(FILE* s, node_handle tnode) const
{
  fprintf(s, "t%d", -tnode);
}

void MEDDLY::ev_forest::writeTerminal(FILE* s, node_handle tnode) const
{
  th_fprintf(s, "t%d", -tnode);
}

MEDDLY::node_handle MEDDLY::ev_forest::readTerminal(FILE* s)
{
  stripWS(s);
  char c = fgetc(s);
  if ('t' == c) {
    int N;
    if (1==fscanf(s, "%d", &N)) {
      if (N>=0) {
        return -N;
      }
    }
  }
  throw error(error::INVALID_FILE);
}

// statics

void MEDDLY::ev_forest::enlargeStatics(int n)
{
  MEDDLY_DCASSERT(n>0);
  if (0==the_ev_cleaner) {
    the_ev_cleaner = new ev_cleanup();
    // DO NOT EVER DELETE the_ev_cleaner, it is done automatically :^)
  }
  if (n>order_size) {
    order = (int*) realloc(order, n*sizeof(int));
    //terminals = (node_handle*) realloc(terminals, n*sizeof(node_handle));
    //if (0==order || 0==terminals) {
    if (0==order) {
      throw error(error::INSUFFICIENT_MEMORY);
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
  // DO NOT delete the_ev_cleaner
  the_ev_cleaner = 0; 
}


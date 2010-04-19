
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



// TODO: Testing

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_LIBGMP
#include <gmp.h>
#endif

#include "../src/defines.h"

// **************************** forest *********************************

forest::forest() {}
forest::~forest() {}


const char* forest::getErrorCodeName(forest::error e)
{
  switch (e) {
    case SUCCESS:
        return "Operation returned successfully";
    case NOT_IMPLEMENTED:
        return "Operation not implemented";
    case INSUFFICIENT_MEMORY:
        return "Operation failed -- lack of memory";
    case INVALID_OPERATION:
        return "Operation not supported for the given forest";
    case INVALID_VARIABLE:
        return "Operation failed -- invalid variable handle";
    case INVALID_ASSIGNMENT:
        return "Operation falied -- variable out of range";
    default:
        return "Unknown error code";
  }
}

// **************************** expert_forest *********************************

expert_forest::expert_forest(domain *d, bool rel, range_type t,
  edge_labeling ev, reduction_rule r, node_storage s, node_deletion_policy nd)
: d(d), isRelation(rel), rangeType(t), edgeLabel(ev),
  reductionRule(r), nodeStorage(s), nodeDeletionPolicy(nd),
  sz(256), firstFree(0), firstHole(-1)
{
  // firstHole < 0 indicates no holes.
  firstHole = -1;
  // Create an array to store pointers to dd_edges.
  edge = (edge_data *) malloc(sz * sizeof(edge_data));
  for (unsigned i = 0; i < sz; ++i) {
    edge[i].nextHole = -1;
    edge[i].edge = 0;
  }
}


void expert_forest::unregisterDDEdges() {
  // Go through the list of valid edges (value > 0), and set
  // the e.index to -1 (indicating unregistered edge).

  // ignore the NULLs; release the rest
  for (unsigned i = 0; i < firstFree; ++i) {
    if (edge[i].edge != 0) {
      DCASSERT(edge[i].nextHole == -1);
      int node = edge[i].edge->getNode();
      unlinkNode(node);
      edge[i].edge->setIndex(-1);
    }
  }

  // firstHole < 0 indicates no holes.
  for (unsigned i = 0; i < firstFree; ++i) {
    edge[i].nextHole = -1;
    edge[i].edge = 0;
  }
  firstHole = -1;
  firstFree = 0;
}


expert_forest::~expert_forest() {
  // Go through the list of valid edges (value > 0), and set
  // the e.index to -1 (indicating unregistered edge).
  // unregisterDDEdges();
  // No need to call this from here -- ~node_manager() calls it.

  // Delete the array.
  free(edge);

  // NOTE: since the user is provided with the dd_edges instances (as opposed
  // to a pointer), the user program will automatically call the
  // destructor for each dd_edge when the corresponding variable goes out of
  // scope. Therefore there is no need to destruct dd_edges from here.

  // Inform the domain that you are going away
  smart_cast<expert_domain*>(d)->unlinkForest(this);
}


void expert_forest::registerEdge(dd_edge& e) {
  // add to collection of edges for this forest.
  // change e.index to help find this edge at a later time.
  if (firstHole >= 0) {
    // hole available; fill it up
    int index = firstHole;
    firstHole = edge[firstHole].nextHole;
    edge[index].edge = &e;
    edge[index].nextHole = -1;
    e.setIndex(index);
  } else {
    // no holes available, add to end of array
    if (firstFree >= sz) {
      // expand edge[]
      int new_sz = sz * 2;
      edge_data* new_edge =
          (edge_data*) realloc(edge, new_sz * sizeof(edge_data));
      assert(new_edge != NULL);
      edge = new_edge;
      for (int i = sz; i < new_sz; ++i)
      {
        edge[i].nextHole = -1;
        edge[i].edge = 0;
      }
      sz = new_sz;
    }
    DCASSERT(firstFree < sz);
    edge[firstFree].nextHole = -1;
    edge[firstFree].edge = &e;
    e.setIndex(firstFree);
    ++firstFree;
  }
}


void expert_forest::unregisterEdge(dd_edge& e) {
  // remove this edge from the collection of edges for this forest.
  // change e.index to -1.
  DCASSERT(e.getIndex() >= 0);
  int index = e.getIndex();
  DCASSERT(edge[index].edge == &e);
  edge[index].edge = 0;
  edge[index].nextHole = firstHole;
  firstHole = index;
  e.setIndex(-1);
}


#ifdef INLINED_REALS

#else

float expert_forest::getReal(int term) const
{
  return (term == 0)? 0.0: *((float*)&(term <<= 1));
}

int expert_forest::getTerminalNode(float a) const
{
#if 1
  // if ((a & 0xffffffff) == 0x00000000) return 0;
  // if (a == 0.0) return 0;
  return (a == 0.0)? 0: (*((int*)((void*)&a)) >> 1) | 0x80000000;

#else
  unsigned int node = *((unsigned int*)((void*)(&a)));
  // printf("%x\n", node);
  return a == 0.0
    ? 0
    // : ((*((int*)((void*)&a))) >> 1) | 0x80000000;
    : (node >> 1) | 0x80000000;
#endif
}

#endif

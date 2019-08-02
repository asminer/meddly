
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"

#include "sccgraph.h"

#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif

// ******************************************************************
// *                                                                *
// *                  sccgraph::edge_label methods                  *
// *                                                                *
// ******************************************************************

MEDDLY::sccgraph::edge_label::edge_label()
{
}

MEDDLY::sccgraph::edge_label::~edge_label()
{
}

// ******************************************************************
// *                                                                *
// *                sccgraph::edge_iterator  methods                *
// *                                                                *
// ******************************************************************

MEDDLY::sccgraph::edge_iterator::edge_iterator(const edge_iterator& ei)
 : parent(ei.parent)
{
  edgeptr = ei.edgeptr;
}

MEDDLY::sccgraph::edge_iterator::edge_iterator(const sccgraph &p, unsigned ep)
 : parent(p)
{
  edgeptr = ep;
}

// ******************************************************************
// *                                                                *
// *                        sccgraph methods                        *
// *                                                                *
// ******************************************************************

MEDDLY::sccgraph::sccgraph()
{
  graph_edges = 0;
  graph_edges_used = 1;   // 0 reserved for null
  graph_edges_alloc = 0;

  graph_from = 0;
  graph_vertices_used = 0;
  graph_vertices_alloc = 0;

  scc_edges = 0;
  scc_edges_used = 1;     // 0 reserved for null
  scc_edges_alloc = 0;

  scc_from = 0;
  scc_vertices_used = 0;
  scc_vertices_alloc = 0;

  vertex_to_scc = 0;
  vertices_by_scc = 0;
  scc_vertex_offset = 0;
}

MEDDLY::sccgraph::~sccgraph()
{
#ifdef DEBUG
  fprintf(stderr, "Entering sccgraph destructor\n");
#endif
  clear();
  free(graph_edges);
  free(graph_from);
  free(scc_edges);
  free(scc_from);
#ifdef DEBUG
  fprintf(stderr, "Exiting  sccgraph destructor\n");
#endif
}

void MEDDLY::sccgraph::clear()
{
  // Clear out edge labels
  if (graph_edges) {
    for (unsigned i=1; i<graph_edges_used; i++) {
      delete graph_edges[i].edge;
      graph_edges[i].edge = 0;
    }
  }
  // Reset sizes
  graph_edges_used = 1;
  graph_vertices_used = 0;
  scc_edges_used = 1;
  scc_vertices_used = 0;
}

void MEDDLY::sccgraph::dumpGraph(MEDDLY::output &out) const
{
  out << "sccgraph structure\n";
  out << "  Graph:\n";
  out << "    " << graph_vertices_used << " vertices\n";
  out << "    Edges:\n";
  for (unsigned i=0; i<graph_vertices_used; i++) {
    out << "      From vertex " << i << ":\n";
    unsigned ptr = graph_from[i];
    if (0==ptr) continue;
    do {
      ptr = graph_edges[ptr].next;
      out << "\tTo vertex " << graph_edges[ptr].to;
      out << " edge label ";
      graph_edges[ptr].edge->show(out);
      out << "\n";
    } while (ptr != graph_from[i]);
  } // for i

  out << "  TBD...\n";
}

void MEDDLY::sccgraph::add_edge(unsigned I, unsigned J, edge_label* L)
{
  MEDDLY_CHECK_RANGE(0, I, graph_vertices_used);
  MEDDLY_CHECK_RANGE(0, J, graph_vertices_used);
  MEDDLY_DCASSERT(L);

  addGraphEdge(I, J, L);

  // TBD - scc edge
}

void MEDDLY::sccgraph::update_SCCs()
{
  // TBD
}

void MEDDLY::sccgraph::expand_vertices(unsigned I)
{ 
  if (I > graph_vertices_alloc) {
    // enlarge
    if (0==graph_vertices_alloc) {
      graph_vertices_alloc = 16;
    } else if (graph_vertices_alloc > 1024) {
      graph_vertices_alloc += 1024;
    } else {
      graph_vertices_alloc *= 2;
    }
    graph_from = (unsigned*) 
      realloc(graph_from, graph_vertices_alloc * sizeof(unsigned));

    if (0==graph_from) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
  }

  // TBD - new sccs also, expand here if needed

  for (unsigned i=graph_vertices_used; i<I; i++) {
    graph_from[i] = 0;
  }

  graph_vertices_used = I;
}

void MEDDLY::sccgraph::addGraphEdge(unsigned I, unsigned J, edge_label* L)
{
  //
  // Put this edge into the list for I.
  // We use a circular list and keep a pointer to the last element;
  // that way "add to front" and "add to end" are both fast.
  //

  MEDDLY_DCASSERT(graph_from);
  MEDDLY_CHECK_RANGE(0, I, graph_vertices_used);

  //
  // Find the place in the list before we need to insert this edge.
  // If we find an existing edge, append L
  //
  // Insertion point will be after prev and before curr.
  //
  unsigned prev = 0;
  unsigned curr = 0;

  if (graph_from[I]) {

    unsigned last = graph_from[I];
    if (J == graph_edges[last].to) {
      // merge edges
      MEDDLY_DCASSERT(graph_edges[last].edge);
      graph_edges[last].edge->append_and_recycle(L);
      return;
    }
    if (J > graph_edges[last].to) {
      prev = last;
    } else {
      // Search the list until curr > J
      for (curr = graph_edges[last].next; curr != last; curr = graph_edges[curr].next) {
        if (graph_edges[curr].to == J ) {
          MEDDLY_DCASSERT(graph_edges[curr].edge);
          graph_edges[curr].edge->append_and_recycle(L);
          return;
        }
        if (graph_edges[curr].to > J ) {
          // found the spot
          break;
        }
        // Still here? curr < J so insertion point must be here or after.
        prev = curr;
      } // for curr
    }
  } // list not empty

  //
  // Need to create a new edge.
  //
  unsigned newedge = graph_edges_used;
  graph_edges_used++;
  if (graph_edges_used > graph_edges_alloc) {
    // enlarge
    if (0==graph_edges_alloc) {
      graph_edges_alloc = 16;
    } else if (graph_edges_alloc > 1024) {
      graph_edges_alloc += 1024;
    } else {
      graph_edges_alloc *= 2;
    }
    graph_edges = (graph_list_node*) 
      realloc(graph_edges, graph_edges_alloc * sizeof(graph_list_node));

    if (0==graph_edges) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    graph_edges[0].edge = 0;
  }

#ifdef DEBUG
  printf("Adding edge# %u  (prev %u curr %u)\n", newedge, prev, curr);
#endif

  //
  // Add the new edge to the list in the appropriate spot
  //

  graph_edges[newedge].to = J;
  graph_edges[newedge].edge = L;
  if (0==prev) {
    MEDDLY_DCASSERT(0==curr);
    MEDDLY_DCASSERT(0==graph_from[I]);
    graph_from[I] = newedge;
    graph_edges[newedge].next = newedge;
    return;
  }

  MEDDLY_DCASSERT(graph_edges[prev].to < J);
  if (curr) {
    MEDDLY_DCASSERT(graph_edges[prev].next == curr);
    MEDDLY_DCASSERT(J < graph_edges[curr].to);

    // We're adding anywhere except the very end of the list.

    graph_edges[newedge].next = curr;
  } else {
    MEDDLY_DCASSERT(last == prev);

    // We're adding at the very end of the list, so we also
    // need to update the list pointer.

    graph_from[I] = newedge;
    graph_edges[newedge].next = graph_edges[prev].next;
  }
  graph_edges[prev].next = newedge;
}




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

  sccs_to_update = 0;
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
  free(vertex_to_scc);
  free(vertices_by_scc);
  free(scc_vertex_offset);
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

inline void show_array(MEDDLY::output &out, const char* what, const unsigned* A, unsigned n)
{
  out << what;
  if (0==A) {
    out << "null\n";
    return;
  }
  out << "[";
  if (n) {
    out << A[0];
    for (unsigned i=1; i<n; i++) out << ", " << A[i];
  }
  out << "]\n";
}

void MEDDLY::sccgraph::dumpGraph(MEDDLY::output &out) const
{
  out << "sccgraph structure\n";
  out << "  Graph:\n";
  out << "    " << graph_vertices_used << " vertices\n";
  out << "    Edges:\n";
  for (unsigned i=0; i<graph_vertices_used; i++) {
    out << "      From vertex " << i << ":\t";
    show_graph_list(out, graph_from[i]);
    out << "\n";
  } // for i

  out << "  Vertex to SCC mapping:\n";
  for (unsigned i=0; i<graph_vertices_used; i++) {
    out << "      Vertex " << i << " is in SCC " << vertex_to_scc[i] << "\n";
  }
  for (unsigned i=0; i<scc_vertices_used; i++) {
    out << "  SCC #" << i << " contains vertices:  ";
    bool commas = false;
    for (unsigned j=scc_vertex_offset[i]; j<scc_vertex_offset[i+1]; j++) {
      if (commas) out << ", ";
      out << vertices_by_scc[j];
      commas = true;
    }
    out << "\n";
  }
  out << "  Edges between SCCs:\n";
  for (unsigned i=0; i<scc_vertices_used; i++) {
    out << "      From SCC#" << i << ":\t";
    show_scc_list(out, scc_from[i]);
    out << "\n";
  } // for i

  out << "  SCCs to explore on update:\n\t";
  show_scc_list(out, sccs_to_update);

  out << "\ninternal structure\n";
  show_array(out, "   graph_from: ", graph_from, graph_vertices_used); 
  show_array(out, "   scc_from: ", scc_from, scc_vertices_used); 
  show_array(out, "   vertex_to_scc: ", vertex_to_scc, graph_vertices_used);
  show_array(out, "   scc_vertex_offset: ", scc_vertex_offset, 1+scc_vertices_used);
  show_array(out, "   vertices_by_scc: ", vertices_by_scc, graph_vertices_used);
  out << "end\n";
}

void MEDDLY::sccgraph::add_edge(unsigned I, unsigned J, edge_label* L)
{
  MEDDLY_CHECK_RANGE(0, I, graph_vertices_used);
  MEDDLY_CHECK_RANGE(0, J, graph_vertices_used);
  MEDDLY_DCASSERT(L);
  MEDDLY_DCASSERT(vertex_to_scc);

  add_graph_edge(I, J, L);
  add_scc_edge(vertex_to_scc[I], vertex_to_scc[J]);
}

void MEDDLY::sccgraph::update_SCCs()
{
  // TBD
}

void MEDDLY::sccgraph::expand_vertices(unsigned I)
{ 
  //
  // Add graph vertices as needed
  //

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

    vertex_to_scc = (unsigned*)
      realloc(vertex_to_scc, graph_vertices_alloc * sizeof(unsigned));

    if (0==vertex_to_scc) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }

    vertices_by_scc = (unsigned*)
      realloc(vertices_by_scc, graph_vertices_alloc * sizeof(unsigned));

    if (0==vertices_by_scc) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
  }

  for (unsigned i=graph_vertices_used; i<I; i++) {
    graph_from[i] = 0;
  }
  const unsigned added_vertices = I - graph_vertices_used;

  //
  // Add SCC vertices as needed
  //

  if (scc_vertices_used + added_vertices > scc_vertices_alloc) {
    // enlarge
    if (0==scc_vertices_alloc) {
      scc_vertices_alloc = 16;
    } else if (scc_vertices_alloc > 1024) {
      scc_vertices_alloc += 1024;
    } else {
      scc_vertices_alloc *= 2;
    }
    scc_from = (unsigned*) 
      realloc(scc_from, scc_vertices_alloc * sizeof(unsigned));

    if (0==scc_from) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }

    scc_vertex_offset = (unsigned*) 
      realloc(scc_vertex_offset, (1+scc_vertices_alloc) * sizeof(unsigned));

    scc_vertex_offset[0] = 0;

    if (0==scc_vertex_offset) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
  }

  //
  // Add a new SCC for each graph vertex added
  //

  for (unsigned i=0; i<added_vertices; i++) {
    scc_from[scc_vertices_used+i] = 0;
    vertex_to_scc[graph_vertices_used+i] = scc_vertices_used+i;

    vertices_by_scc[ scc_vertex_offset[ scc_vertices_used+i ] ] = scc_vertices_used+i;

    scc_vertex_offset[scc_vertices_used+i+1] = scc_vertex_offset[scc_vertices_used+i]+1;
  }

  graph_vertices_used = I;
  scc_vertices_used += added_vertices;
}

void MEDDLY::sccgraph::add_graph_edge(unsigned I, unsigned J, edge_label* L)
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
  // Create a new edge, add it to the list in the appropriate spot
  //
  unsigned newedge = new_graph_edge(J, L);

#ifdef DEBUG
  printf("Adding graph edge# %u  (prev %u curr %u)\n", newedge, prev, curr);
#endif

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

unsigned MEDDLY::sccgraph::new_graph_edge(unsigned J, edge_label* L)
{
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

  graph_edges[newedge].to = J;
  graph_edges[newedge].edge = L;

  return newedge;
}

void MEDDLY::sccgraph::add_scc_edge(unsigned I, unsigned J)
{
  if (I==J) return;

  MEDDLY_DCASSERT(scc_from);
  MEDDLY_CHECK_RANGE(0, I, scc_vertices_used);

  scc_from[I] = add_to_scc_list(scc_from[I], J);

  //
  // If we need to update SCC I, add it to the (ordered) list
  //

  if (I>J) {
    sccs_to_update = add_to_scc_list(sccs_to_update, I);
  }
}

unsigned MEDDLY::sccgraph::new_scc_edge(unsigned J)
{
  unsigned newedge = scc_edges_used;
  scc_edges_used++;
  if (scc_edges_used > scc_edges_alloc) {
    // enlarge
    if (0==scc_edges_alloc) {
      scc_edges_alloc = 16;
    } else if (scc_edges_alloc > 1024) {
      scc_edges_alloc += 1024;
    } else {
      scc_edges_alloc *= 2;
    }
    scc_edges = (scc_list_node*) 
      realloc(scc_edges, scc_edges_alloc * sizeof(scc_list_node));

    if (0==scc_edges) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
  }

  scc_edges[newedge].to = J;

  return newedge;
}

unsigned MEDDLY::sccgraph::add_to_scc_list(unsigned ptr, unsigned J)
{
  //
  // Find the place in the circular list before we need to insert this edge.
  // If we find an existing edge, do nothing and return.
  //
  // Insertion point will be after prev and before curr.
  //
  unsigned prev = 0;
  unsigned curr = 0;

  if (ptr) {

    unsigned last = ptr;
    if (J == scc_edges[last].to) {
      // edge already exists
      MEDDLY_DCASSERT(scc_edges[last].edge);
      return ptr;
    }
    if (J > scc_edges[last].to) {
      prev = last;
    } else {
      // Search the list until curr > J
      for (curr = scc_edges[last].next; curr != last; curr = scc_edges[curr].next) {
        if (scc_edges[curr].to == J ) {
          MEDDLY_DCASSERT(scc_edges[curr].edge);
          // Already in the list, do nothing
          return ptr;
        }
        if (scc_edges[curr].to > J ) {
          // found the spot
          break;
        }
        // Still here? curr < J so insertion point must be here or after.
        prev = curr;
      } // for curr
    }
  } // list not empty

  //
  // Create a new edge, add it to the list in the appropriate spot
  //
  unsigned newedge = new_scc_edge(J);

#ifdef DEBUG
  printf("Adding to SCC list item# %u  (prev %u curr %u)\n", newedge, prev, curr);
#endif

  if (0==prev) {
    // We had an empty list to start with
    MEDDLY_DCASSERT(0==curr);
    MEDDLY_DCASSERT(0==ptr);
    scc_edges[newedge].next = newedge;
    return newedge;
  }

  MEDDLY_DCASSERT(scc_edges[prev].to < J);
  if (curr) {
    MEDDLY_DCASSERT(scc_edges[prev].next == curr);
    MEDDLY_DCASSERT(J < scc_edges[curr].to);

    // We're adding anywhere except the very end of the list.

    scc_edges[newedge].next = curr;
  } else {
    MEDDLY_DCASSERT(last == prev);

    // We're adding at the very end of the list, so we also
    // need to update the list pointer.

    ptr = newedge;
    scc_edges[newedge].next = scc_edges[prev].next;
  }
  scc_edges[prev].next = newedge;
  return ptr;
}



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
#define DEBUG_SCC
// #define DEBUG_SCC_MIN

#ifdef DEBUG
#include <cstdio>
#endif

#ifdef DEBUG_SCC
#include <cstdio>
#endif

#ifdef DEBUG_SCC_MIN
#include <cstdio>
#endif


/**
        Add node to a circularly-linked list, in order of "to" field.
          @param  List    Pointer to LAST element in the list.
                          It's circular, so this allows us to add quickly
                          at both the beginning and end of the list.
                          Will be updated if needed.

          @param  node    Pointer to the new node.

          @param  edges   Edge array that holds the actual edges
                          (our "pointers" are indexes into this array).

          @return         Pointer to a node.  If there was already an element
                          in the list with the same "to" value, return a
                          pointer to that element so the caller can merge or
                          whatever.  Otherwise, returns node.
*/
template <class LNT>
unsigned add_to_circular_list(unsigned& List, unsigned node, LNT* edges)
{
      MEDDLY_DCASSERT(edges);
      MEDDLY_DCASSERT(node);

      if (0==List) {
        // Trivial case - empty list
        edges[node].next = node;
        List = node;
        return node;
      }
      if (edges[node].to == edges[List].to) {
        //
        // Trivial case - node matches end of the list
        //
        return List;
      }
      if (edges[node].to > edges[List].to) {
        //
        // Trivial case - node goes at the end of the list
        // Note - this is the ONLY case where List is updated.
        //
        edges[node].next = edges[List].next;
        edges[List].next = node;
        List = node;
        return node;
      }

      //
      // There must be some element in the list with to value not less than than node's.
      // Use pointer curr to find the smallest one, with prev before that.
      //
      unsigned prev = List;
      unsigned curr = edges[prev].next;

      for (;;) {
        if (edges[node].to == edges[curr].to) {
          // 
          // We found a match, don't update anything in the list
          //
          return curr;
        }
        if (edges[node].to < edges[curr].to) {
          //
          // Insert here.  We want the list to become:
          //
          // prev -> node -> curr
          //
          edges[node].next = curr;
          edges[prev].next = node;
          return node;
        }

        //
        // Advance the pointers
        //
        prev = curr;
        curr = edges[curr].next;
        MEDDLY_DCASSERT(prev != List);    // we should never loop back around
      } // infinite loop
}

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
  graph_edges_freelist = 0;

  graph_from = 0;
  graph_vertices_used = 0;
  graph_vertices_alloc = 0;

  scc_edges = 0;
  scc_edges_used = 1;     // 0 reserved for null
  scc_edges_alloc = 0;
  scc_edges_freelist = 0;

  scc_from = 0;
  scc_vertices_used = 0;
  scc_vertices_alloc = 0;

  vertex_to_scc = 0;
  vertices_by_scc = 0;
  scc_vertex_offset = 0;

  scc_updatelist = 0;

  visit_stack = 0;
  visit_index = 0;
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
  free(visit_stack);
  free(visit_index);
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
  show_scc_list(out, scc_updatelist);

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

  //
  // Add edge to the main graph
  //
  unsigned newge = new_graph_edge(J, L);
  unsigned ge = add_to_circular_list(graph_from[I], newge, graph_edges);
  if (ge != newge) {
    // Existing edge; merge them!
    graph_edges[ge].edge->append_and_recycle(L);
    graph_edges[newge].edge = 0; 
    recycle_graph_edge(newge);
  }

  //
  // Add edge to the scc graph
  //
  unsigned newse = new_scc_edge(vertex_to_scc[J]);
  if (newse != add_to_circular_list(scc_from[vertex_to_scc[I]], newse, scc_edges)) {
    // Existing edge
    recycle_scc_edge(newse);
  }

  //
  // Do we need to update the scc graph?
  //
  if (vertex_to_scc[I]>vertex_to_scc[J]) {
    // Add SCC I to the update list
    unsigned nu = new_scc_edge(vertex_to_scc[I]);
    if (nu != add_to_circular_list(scc_updatelist, nu, scc_edges)) {
      recycle_scc_edge(nu);
    }
  }
}

void MEDDLY::sccgraph::update_SCCs()
{
  if (0==scc_updatelist) return;

  //
  // Initialize 'globals' before scc recursions
  //

  stack_top = 0;
  for (unsigned i=0; i<scc_vertices_used; i++) {
    visit_index[i] = 0;
  }
  curr_index = 0;

  //
  // Reverse the update list because we need
  // to visit them from highest to lowest.
  //
  unsigned ptr = scc_updatelist;
  unsigned next = 0;
  unsigned revlist = 0;
  next = scc_edges[ptr].next;
  do {
    ptr = next;
    next = scc_edges[ptr].next;
    scc_edges[ptr].next = revlist;
    revlist = ptr;
  } while (ptr != scc_updatelist);
  scc_updatelist = 0;

  //
  // Traverse list of sccs we know need updating
  // and do the scc traversal algorithm on them.
  // Also recycle the list.
  //
  scc_number = 2*scc_vertices_used;
  while (revlist) {
    ptr = revlist;
    revlist = scc_edges[revlist].next;

#ifdef DEBUG_SCC
    printf("Updating scc %u\n", ptr);
#endif
    unsigned v = scc_edges[ptr].to;
    if (visit_index[v]) continue;
    scc_visit(v);
    recycle_scc_edge(ptr);
  }
  
#ifdef DEBUG_SCC
  printf("Done exploring sccs\n");
  FILE_output cout(stdout);
  show_array(cout, "  visit_index: ", visit_index, scc_vertices_used);
#endif

  //
  // Determine renumbering
  //

  //
  // Pass 1: anything with visit_index=0 is numbered first, in order.
  // These vertices were never touched during exploration.
  // Note - if visit_index != 0 then we know visit_index >= scc_vertices_used.
  //
  unsigned up_counter = 0;
  for (unsigned i=0; i<scc_vertices_used; i++) {
    if (0==visit_index[i]) {
      visit_index[i] = up_counter++;
    }
  }

  //
  // There will be a gap of exactly scc_number - up_counter
  // between the ones numbered upwards (that were 0) and the
  // SCCs that were discovered counting downwards.
  // Pass 2: collapse the gap.
  //
  const unsigned gap = 1+scc_number - up_counter;
  for (unsigned i=0; i<scc_vertices_used; i++) {
    if (visit_index[i] > up_counter) {
      visit_index[i] -= gap;
    }
  }

#ifdef DEBUG_SCC
  show_array(cout, "  renumbering: ", visit_index, scc_vertices_used);
#endif

  //
  // Renumbering has been determined!
  // visit_index is now used for SCC renumbering.
  // Next - renumber everything :)
  //

  //
  // Move the scc_graph adjacency lists to visit_stack, temporarily.
  //
  for (unsigned i=0; i<scc_vertices_used; i++) {
    visit_stack[i] = scc_from[i];
    scc_from[i] = 0;
  }
  
  //
  // Copy edges back where they belong.
  // During the copy, renumber the "to" value.
  // Note - we can't simply do that in place, because the lists
  // will be out of order and some elements may need to be merged.
  // Also, eliminate any self loops.
  //
  for (unsigned i=0; i<scc_vertices_used; i++) {
    if (0==visit_stack[i]) continue;
    const unsigned newi = visit_index[i];
    // convert circular list to linear list by breaking the cycle
    unsigned list = scc_edges[visit_stack[i]].next;
    scc_edges[visit_stack[i]].next = 0;
    // transfer from list to the correct list, renumbering to values as we go
    while (list) {
      next = scc_edges[list].next;
      const unsigned newj = visit_index[ scc_edges[list].to ];
      scc_edges[list].to = newj;

      // Check our invariant
      MEDDLY_DCASSERT(newj >= newi);
      
      if ((newi == newj) || (list != add_to_circular_list(scc_from[newi], list, scc_edges))) {
        // this edge collapses now
        recycle_scc_edge(list);
      }
      list = next;
    } // while
  } // for i

  //
  // Update vertex_to_scc mapping
  //
  for (unsigned i=0; i<graph_vertices_used; i++) {
    MEDDLY_CHECK_RANGE(0, vertex_to_scc[i], scc_vertices_used);
    vertex_to_scc[i] = visit_index[ vertex_to_scc[i] ];
  }
  
  // TBD:
  // rebuild the vertices_by_scc and scc_vertex_offset arrays

  // TBD:
  // shrink scc_vertices_used as needed (first scc# with no vertices)
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

    visit_stack = (unsigned*)
      realloc(visit_stack, scc_vertices_alloc * sizeof(unsigned));

    if (0==visit_stack) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }

    visit_index = (unsigned*)
      realloc(visit_index, scc_vertices_alloc * sizeof(unsigned));
      
    if (0==visit_index) {
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


unsigned MEDDLY::sccgraph::new_graph_edge(unsigned J, edge_label* L)
{
  unsigned newedge = 0;

  if (graph_edges_freelist) {
    //
    // Use a recycled edge
    //
    newedge = graph_edges_freelist;
    graph_edges_freelist = graph_edges[graph_edges_freelist].next;
  }
  else {
    //
    // Use the next edge in the array
    //
    newedge = graph_edges_used;
    graph_edges_used++;
    if (graph_edges_used > graph_edges_alloc) {
      //
      // enlarge the edge array
      //
      if (0==graph_edges_alloc) {
        graph_edges_alloc = 16;
      } else if (graph_edges_alloc > 1024) {
        graph_edges_alloc += 1024;
      } else {
        graph_edges_alloc *= 2;
      }
      graph_edges = (labeled_list_node*) 
        realloc(graph_edges, graph_edges_alloc * sizeof(labeled_list_node));

      if (0==graph_edges) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      }
      graph_edges[0].edge = 0;
    }
  }

  graph_edges[newedge].to = J;
  graph_edges[newedge].edge = L;
  return newedge;
}


unsigned MEDDLY::sccgraph::new_scc_edge(unsigned J)
{
  unsigned newedge = 0;

  if (scc_edges_freelist) {
    //
    // Use a recycled edge
    //
    newedge = scc_edges_freelist;
    scc_edges_freelist = scc_edges[scc_edges_freelist].next;
  }
  else {
    newedge = scc_edges_used;
    scc_edges_used++;
    if (scc_edges_used > scc_edges_alloc) {
      //
      // Enlarge edge array
      //
      if (0==scc_edges_alloc) {
        scc_edges_alloc = 16;
      } else if (scc_edges_alloc > 1024) {
        scc_edges_alloc += 1024;
      } else {
        scc_edges_alloc *= 2;
      }
      scc_edges = (unlabeled_list_node*) 
        realloc(scc_edges, scc_edges_alloc * sizeof(unlabeled_list_node));

      if (0==scc_edges) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      }
    }
  }

  scc_edges[newedge].to = J;

  return newedge;
}

unsigned MEDDLY::sccgraph::scc_visit(unsigned v)
{
#ifdef DEBUG_SCC
  printf("entering scc_visit(%u)\n", v);
#endif

  visit_index[v] = ++curr_index;
  unsigned min = visit_index[v];
  visit_push(v);

#ifdef DEBUG_SCC_MIN
  printf("    min for %u is initially %u\n", v, min);
#endif

  // go through outgoing edges from v
  if (scc_from[v]) {
    unsigned ptr = scc_from[v];
    do {
      ptr = scc_edges[ptr].next;
      unsigned w = scc_edges[ptr].to;

      if (visit_index[w]) {
        // already visited
        min = MIN(min, visit_index[w]);
      } else {
        min = MIN(min, scc_visit(w));
      }
#ifdef DEBUG_SCC_MIN
      printf("    min for %u is now %u\n", v, min);
#endif
    } while (ptr != scc_from[v]);
  }

#ifdef DEBUG_SCC_MIN
  printf("    min for %u is finally %u\n", v, min);
#endif

  // check if v is the root of an SCC
  if (min == visit_index[v]) {
    // Pop everything off the stack until v;
    // these vertices must be merged
    // But first, determine the minimum of everything that will be popped;
    // use this as the new scc number.
  //  const unsigned sccnum = min_visit_stack_until(v);
#ifdef DEBUG_SCC
    printf("  %u is an scc root\n", v);
    printf("  Create new SCC# %u from vertices: ", scc_number);
   // printf("  Create new SCC# %u from vertices:  ", sccnum);
#endif
    unsigned w;
    do {
      w = visit_pop();
#ifdef DEBUG_SCC
      printf("%u ", w);
#endif
      visit_index[w] = scc_number;
    } while (w != v);
#ifdef DEBUG_SCC
    printf("\n");
#endif
    scc_number--;   // we get a reverse topological sort for free; number backwards
  }

#ifdef DEBUG_SCC
  printf("exiting  scc_visit(%u) returning %u\n", v, min);
#endif

  return min;
}


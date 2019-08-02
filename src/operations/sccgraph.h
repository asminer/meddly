
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

#ifndef SCCGRAPH_H
#define SCCGRAPH_H

namespace MEDDLY {
  class sccgraph;
};

/**
  Explicit graph representation with SCCs.
  Useful for fixed point computations such as saturation.

  To use this class, derive a class from subclass edge_label
  and implement append_and_recycle() which is used to merge
  edge labels.
*/
class MEDDLY::sccgraph {

  public:
    class edge_label {
      public:
        edge_label();
        virtual ~edge_label();
        virtual void append_and_recycle(edge_label*) = 0;
        virtual void show(MEDDLY::output &out) = 0;
    };

  public:
    class edge_iterator {
        const sccgraph& parent;
        unsigned edgeptr;
      public:
        edge_iterator(const edge_iterator& ei);
        
        unsigned to() const;              // Inlined
        const edge_label* edge() const;   // Inlined

        edge_iterator& operator++();      // Inlined
        bool operator==(const edge_iterator& ei) const;   // Inlined

      private:
        edge_iterator(const sccgraph &parent, unsigned ep);
        friend class sccgraph;
    };

  public:
    sccgraph();
    ~sccgraph();

    /// Reset the graph but keep the memory.
    void clear();

    /// For debugging and testing.
    void dumpGraph(MEDDLY::output &out) const; 

  public:
    //
    // Graph info
    // 

    /**
        Add vertices 0,...,I to the graph, unless already present.
        Inlined, implemented below.
    */
    void add_vertex(unsigned I);

    /**
        Add edge I->J with label L to the graph.
        If there's already an edge I->J, merge the edge labels
        using append_and_recycle().
    */
    void add_edge(unsigned I, unsigned J, edge_label* L);
    
    /// Update SCCs if necessary
    void update_SCCs();

  private:
    void expand_vertices(unsigned I);

  public:
    //
    // Iterator for graph edges
    //

    /**
        Start iterating graph edges from vertex I
    */
    edge_iterator begin_from(unsigned I) const;   // Inlined
    edge_iterator end() const;                    // Inlined

  public:
    //
    // SCC info
    //

    /// Current number of SCCs.
    unsigned get_num_SCCs() const;    // Inlined

    /// Return the SCC that contains graph vertex I.
    unsigned get_SCC_containing(unsigned I) const;    // Inlined

    /**
        Return the size of an SCC.
          @param  s   SCC, should be between 0 and get_num_SCCs() - 1.
          @return     Number of vertexes contained in SCC s.
    */
    unsigned get_SCC_size(unsigned s) const;    // Inlined
    
    /**
        Return the list of vertexes in an SCC.
          @param  s   SCC, should be between 0 and get_num_SCCs() - 1.
          @return     "Array" A such that A[i] gives a vertex
                      in SCC s, for all i between 0 and get_SCC_size(s) - 1.
    */
    const unsigned* get_SCC_vertices(unsigned s) const;   // Inlined

  private:
    // Complete graph

    struct graph_list_node {
      edge_label* edge;
      unsigned to;
      unsigned next;
    };

    graph_list_node* graph_edges;
    unsigned graph_edges_used;
    unsigned graph_edges_alloc;

    unsigned* graph_from;
    unsigned graph_vertices_used;
    unsigned graph_vertices_alloc;

  private:
    void addGraphEdge(unsigned I, unsigned J, edge_label* L);

  private:
    // SCC graph

    struct scc_list_node {
      unsigned to;
      unsigned next;
    };

    scc_list_node* scc_edges;
    unsigned scc_edges_used;
    unsigned scc_edges_alloc;

    unsigned* scc_from;
    unsigned scc_vertices_used;
    unsigned scc_vertices_alloc;

  private:
    // mappings from vertices to SCCs

    unsigned* vertex_to_scc;  // size is graph_vertices_alloc

    unsigned* vertices_by_scc;  // size is graph_vertices_alloc
    unsigned* scc_vertex_offset;  // size is 1+scc_vertices_alloc

  private:
    friend class sccgraph::edge_iterator;
};

//
// Inlined stuff
//

unsigned MEDDLY::sccgraph::edge_iterator::to() const
{
  MEDDLY_DCASSERT(edgeptr);
  return parent.graph_edges[edgeptr].to;
}

const MEDDLY::sccgraph::edge_label* MEDDLY::sccgraph::edge_iterator::edge() const
{
  MEDDLY_DCASSERT(edgeptr);
  return parent.graph_edges[edgeptr].edge;
}

MEDDLY::sccgraph::edge_iterator& MEDDLY::sccgraph::edge_iterator::operator++()
{
  if (edgeptr) {
    // TBD - sanity checks
    edgeptr = parent.graph_edges[edgeptr].next;
  }
  return *this;
}

bool MEDDLY::sccgraph::edge_iterator::operator==(const edge_iterator& ei) const
{
  MEDDLY_DCASSERT(parent == ei.parent);
  return edgeptr == ei.edgeptr;
}

MEDDLY::sccgraph::edge_iterator MEDDLY::sccgraph::begin_from(unsigned I) const
{
  MEDDLY_CHECK_RANGE(0, I, graph_vertices_used);
  MEDDLY_DCASSERT(graph_from);
  return edge_iterator(*this, graph_from[I]);
}

MEDDLY::sccgraph::edge_iterator MEDDLY::sccgraph::end() const
{
  return edge_iterator(*this, 0);
}

inline void MEDDLY::sccgraph::add_vertex(unsigned I)
{
  if (I >= graph_vertices_used) {
    expand_vertices(I+1);
  }
}

inline unsigned MEDDLY::sccgraph::get_num_SCCs() const
{
  return scc_vertices_used;
}

inline unsigned MEDDLY::sccgraph::get_SCC_containing(unsigned I) const
{
  MEDDLY_DCASSERT(vertex_to_scc);
  MEDDLY_CHECK_RANGE(0, I, graph_vertices_used);
  return vertex_to_scc[I];
}

inline unsigned MEDDLY::sccgraph::get_SCC_size(unsigned s) const
{
  MEDDLY_DCASSERT(scc_vertex_offset);
  MEDDLY_CHECK_RANGE(0, s, scc_vertices_used);
  return scc_vertex_offset[s+1] - scc_vertex_offset[s];
}

inline const unsigned* MEDDLY::sccgraph::get_SCC_vertices(unsigned s) const
{
  MEDDLY_DCASSERT(vertices_by_scc);
  MEDDLY_DCASSERT(scc_vertex_offset);
  MEDDLY_CHECK_RANGE(0, s, scc_vertices_used);
  return vertices_by_scc + scc_vertex_offset[s];
}

#endif  // ifndef SCCGRAPH_H


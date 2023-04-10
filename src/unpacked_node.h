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

#ifndef MEDDLY_UNPACKED_NODE_H
#define MEDDLY_UNPACKED_NODE_H

#include "defines.h"
#include "encoders.h"

namespace MEDDLY {
    class unpacked_node;
    class expert_forest;
}

// ******************************************************************
// *                                                                *
// *                      unpacked_node  class                      *
// *                                                                *
// ******************************************************************

/** Class for reading nodes.
    Ideally - used anywhere we want to read node data.
    Backend implementation may change :^)
    Implemented in node_wrappers.cc.
    Readers may be "full" or "sparse",
    regardless of how the actual node is stored.

    Currently, a node is:
      array of node_handles, for the downward pointers
      array of integers, for indexes (sparse only)
      "compact" chunk of memory for edge values
*/
class MEDDLY::unpacked_node {
  public:
    /**
        Options for filling an unpacked node from an existing one.
    */
      /*
    enum storage_style {
      /// Unpacked node should be stored as truncated full
      FULL_NODE,
      /// Unpacked node should be stored sparsely
      SPARSE_NODE,
      /// Unpacked node should be stored same as packed node
      AS_STORED
    };
    */

  public:
    /** Constructor.
     The class must be "filled" by a forest before
     it can be used, however.
     */
    unpacked_node();

    /// Destructor.
    ~unpacked_node();

    /// Free memory, but don't delete.
    void clear();

  public:
  /* Initialization methods, primarily for reading */

    void initFromNode(const expert_forest *f, node_handle node, bool full);
    void initFromNode(const expert_forest *f, node_handle node, storage_style st2);

    void initRedundant(const expert_forest *f, int k, node_handle node, bool full);
    void initRedundant(const expert_forest *f, int k, int ev, node_handle node, bool full);
    void initRedundant(const expert_forest *f, int k, long ev, node_handle node, bool full);
    void initRedundant(const expert_forest *f, int k, float ev, node_handle node, bool full);

    void initIdentity(const expert_forest *f, int k, unsigned i, node_handle node, bool full);
    void initIdentity(const expert_forest *f, int k, unsigned i, int ev, node_handle node, bool full);
    void initIdentity(const expert_forest *f, int k, unsigned i, long ev, node_handle node, bool full);
    void initIdentity(const expert_forest *f, int k, unsigned i, float ev, node_handle node, bool full);

  /* Create blank node, primarily for writing */

    void initFull(const expert_forest *f, int level, unsigned tsz);
    void initSparse(const expert_forest *f, int level, unsigned nnz);

  public:
  /* For convenience: get recycled instance and initialize */

    static unpacked_node* newFromNode(const expert_forest *f, node_handle node, bool full);
    static unpacked_node* newFromNode(const expert_forest *f, node_handle node, storage_style st2);

    static unpacked_node* newRedundant(const expert_forest *f, int k, node_handle node, bool full);
    static unpacked_node* newRedundant(const expert_forest *f, int k, long ev, node_handle node, bool full);
    static unpacked_node* newRedundant(const expert_forest *f, int k, float ev, node_handle node, bool full);

    static unpacked_node* newIdentity(const expert_forest *f, int k, unsigned i, node_handle node, bool full);
    static unpacked_node* newIdentity(const expert_forest *f, int k, unsigned i, long ev, node_handle node, bool full);
    static unpacked_node* newIdentity(const expert_forest *f, int k, unsigned i, float ev, node_handle node, bool full);

    /** Create a zeroed-out full node */
    static unpacked_node* newFull(const expert_forest *f, int level, unsigned tsz);

    /** Create a zeroed-out sparse node */
    static unpacked_node* newSparse(const expert_forest *f, int level, unsigned nnz);

  public:
  /* Display / access methods */

    /** Write a node in human-readable format.

        @param  s       Output stream.
        @param  details Should we show "details" or not.
    */
    void show(output &s, bool details) const;

    /** Write a node in machine-readable format.

        @param  s       Output stream.
        @param  map     Translation to use on node handles.
                        Allows us to renumber nodes as we write them.
    */
    void write(output &s, const node_handle* map) const;

    /// Get a pointer to the unhashed header data.
    const void* UHptr() const;

    /// Modify a pointer to the unhashed header data
    void* UHdata();

    /// Get the number of bytes of unhashed header data.
    unsigned UHbytes() const;

    /// Get a pointer to the hashed header data.
    const void* HHptr() const;

    /// Modify a pointer to the hashed header data.
    void* HHdata();

    /// Get the number of bytes of hashed header data.
    unsigned HHbytes() const;

    /** Get a downward pointer.
          @param  n   Which pointer.
          @return     If this is a full reader,
                      return pointer with index n.
                      If this is a sparse reader,
                      return the nth non-zero pointer.
    */
    node_handle d(unsigned n) const;

    /** Reference to a downward pointer.
          @param  n   Which pointer.
          @return     If this is a full reader,
                      modify pointer with index n.
                      If this is a sparse reader,
                      modify the nth non-zero pointer.
    */
    node_handle& d_ref(unsigned n);

    /// Set the nth pointer from E, and destroy E.
    void set_d(unsigned n, dd_edge &E);

    /** Get the index of the nth non-zero pointer.
        Use only for sparse readers.
     */
    unsigned i(unsigned n) const;

    /** Modify the index of the nth non-zero pointer.
        Use only for sparse readers.
    */
    unsigned& i_ref(unsigned n);

    /// Get a pointer to an edge
    const void* eptr(unsigned i) const;

    /// Modify pointer to an edge
    void* eptr_write(unsigned i);

    /// Set the nth pointer and edge value from E, and destroy E.
    void set_de(unsigned n, dd_edge &E);

    /// Get the edge value, as an integer.
    void getEdge(unsigned i, long& ev) const;

    /// Get the edge value, as a float.
    void getEdge(unsigned i, float& ev) const;

    /// Set the edge value, as an integer.
    void setEdge(unsigned i, long ev);

    /// Set the edge value, as a float.
    void setEdge(unsigned i, float ev);

    /// Get the edge value, as an integer.
    long ei(unsigned i) const;

    /// Get the edge value, as a float.
    float ef(unsigned i) const;

    // -------------------------------------------------------------------------
    // Methods to access the extensible portion of the node
    //
    // Note: Only nodes that belong to an extensible level can be extensible.
    //
    // Extensible node  : (Regular part, Extensible part).
    // Regular Part     : same as non-extensible nodes.
    // Extensible Part  : a single edge, i.e. a tuple <index, node, edge-value>,
    //                    for all indices in [extensible_index, +infinity].
    //
    // Example:
    // The node
    //    [0:n0:ev0, 1:n1:ev1, 2:n2:ev2, 3:n2:ev2, ..., +infinity:n2:ev2]
    // is represented as the following extensible node:
    //    ([0:n0:ev0, 1:n1:ev1], Extensible: [2:n2:ev2])
    // with,
    //    size of node           : 2
    //    extensible index       : 2
    //    extensible node_handle : n2
    //    extensible edge-value  : ev2
    // -------------------------------------------------------------------------

    /// Does this reader store an extensible edge?
    /// Note: in an extensible node, all edges starting from
    ///       index to +infinity refer to the same edge, i.e. (node, edge-value).
    bool isExtensible() const;

    void markAsExtensible();
    void markAsNotExtensible();

    node_handle ext_d() const;
    unsigned ext_i() const;
    long ext_ei() const;
    float ext_ef() const;

    // -------------------- End of extensible portion --------------------------

    /// Get the level number of this node.
    int getLevel() const;

    /// Set the level number of this node.
    void setLevel(int k);

    /// Get the size of this node (full readers only).
    unsigned getSize() const;

    /// Get the number of nonzeroes of this node (sparse readers only).
    unsigned getNNZs() const;

    /// Is this a sparse reader?
    bool isSparse() const;

    /// Is this a full reader?
    bool isFull() const;

    /// Does this node have edge values?
    bool hasEdges() const;

    /// Number of bytes per edge
    unsigned edgeBytes() const;

    // For debugging unique table.
    unsigned hash() const;

    void setHash(unsigned H);

    void computeHash();

    /// Removes redundant trailing edges.
    /// If the unpacked node is sparse, it assumes its indices to be in ascending order.
    void trim();

    /// If the unpacked node is sparse, it is sorted so that the
    /// indices are in ascending order.
    void sort();

    /// Checks if the node is has no trailing redundant edges
    bool isTrim() const;

    // checks if the node indices are in ascending order
    bool isSorted() const;

    // Is this a "build" node?  Important for mark and sweep.
    bool isBuildNode() const;

  public:

    /// Change the size of a node
    void resize(unsigned ns);

    /// Shrink the size of a (truncated) full node
    void shrinkFull(unsigned ns);

    /// Shrink the size of a sparse node
    void shrinkSparse(unsigned ns);

    /// Called within expert_forest to allocate space.
    ///   @param  p     Parent.
    ///   @param  k     Level number.
    ///   @param  ns    Size of node.
    ///   @param  full  If true, we'll be filling a full reader.
    ///                 Otherwise it is a sparse one.
    void bind_to_forest(const expert_forest* p, int k, unsigned ns, bool full);

    /// Called by node_storage when building an unpacked
    /// node based on how it's stored.
    void bind_as_full(bool full);

  protected:
    void clearFullEdges();
    void clearSparseEdges();

  public:
    // Centralized recycling
    static unpacked_node* useUnpackedNode();
    static void recycle(unpacked_node* r);
    static void freeRecycled();

    static void addToBuildList(unpacked_node* b);
    static void removeFromBuildList(unpacked_node* b);
    static void markBuildListChildren(expert_forest* F);

  private:
    const expert_forest* parent;
    static unpacked_node* freeList;
    static unpacked_node* buildList;

    unpacked_node* next; // for recycled list, and list of nodes being built
    bool is_in_build_list;

    /*
      TBD - extra info that is not hashed
    */
    void* extra_unhashed;
    unsigned ext_uh_alloc;
    unsigned ext_uh_size;
    /*
     Extra info that is hashed
     */
    void* extra_hashed;
    unsigned ext_h_alloc;
    unsigned ext_h_size;
    /*
     Down pointers, indexes, edge values.
     */
    node_handle* down;
    unsigned* index;
    void* edge;
    bool is_extensible;
    unsigned alloc;
    unsigned ealloc;
    unsigned size;
    unsigned nnzs;
    int level;
    unsigned h;
    unsigned char edge_bytes; // number of bytes for an edge value.
    bool is_full;
#ifdef DEVELOPMENT_CODE
    bool has_hash;
#endif
};

// ******************************************************************
// *                                                                *
// *                 inlined  unpacked_node methods                 *
// *                                                                *
// ******************************************************************


inline void
MEDDLY::unpacked_node::initFromNode(const expert_forest *f,
  node_handle node, bool full)
{
  MEDDLY_DCASSERT(f);
  f->fillUnpacked(*this, node, full ? FULL_NODE : SPARSE_NODE);
}

inline void
MEDDLY::unpacked_node::initFromNode(const expert_forest *f,
  node_handle node, storage_style st2)
{
  MEDDLY_DCASSERT(f);
  f->fillUnpacked(*this, node, st2);
}

inline void MEDDLY::unpacked_node::initFull(const expert_forest *f, int levl, unsigned tsz)
{
  MEDDLY_DCASSERT(f);
  bind_to_forest(f, levl, tsz, true);
}

inline void MEDDLY::unpacked_node::initSparse(const expert_forest *f, int levl, unsigned nnz)
{
  MEDDLY_DCASSERT(f);
  bind_to_forest(f, levl, nnz, false);
}

// ****************************************************************************

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::newFromNode(const expert_forest *f, node_handle node, bool full)
{
  unpacked_node* U = useUnpackedNode();
  MEDDLY_DCASSERT(U);
  U->initFromNode(f, node, full);
  return U;
}

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::newFromNode(const expert_forest *f, node_handle node, storage_style st2)
{
  unpacked_node* U = useUnpackedNode();
  MEDDLY_DCASSERT(U);
  U->initFromNode(f, node, st2);
  return U;
}

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::newRedundant(const expert_forest *f, int k, node_handle node, bool full)
{
  unpacked_node* U = useUnpackedNode();
  MEDDLY_DCASSERT(U);
  U->initRedundant(f, k, node, full);
  return U;
}

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::newRedundant(const expert_forest *f, int k, long ev, node_handle node, bool full)
{
  unpacked_node* U = useUnpackedNode();
  MEDDLY_DCASSERT(U);
  U->initRedundant(f, k, ev, node, full);
  return U;
}

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::newRedundant(const expert_forest *f, int k, float ev, node_handle node, bool full)
{
  unpacked_node* U = useUnpackedNode();
  MEDDLY_DCASSERT(U);
  U->initRedundant(f, k, ev, node, full);
  return U;
}

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::newIdentity(const expert_forest *f, int k, unsigned i, node_handle node, bool full)
{
  unpacked_node* U = useUnpackedNode();
  MEDDLY_DCASSERT(U);
  U->initIdentity(f, k, i, node, full);
  return U;
}

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::newIdentity(const expert_forest *f, int k, unsigned i, long ev, node_handle node, bool full)
{
  unpacked_node* U = useUnpackedNode();
  MEDDLY_DCASSERT(U);
  U->initIdentity(f, k, i, ev, node, full);
  return U;
}

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::newIdentity(const expert_forest *f, int k, unsigned i, float ev, node_handle node, bool full)
{
  unpacked_node* U = useUnpackedNode();
  MEDDLY_DCASSERT(U);
  U->initIdentity(f, k, i, ev, node, full);
  return U;
}

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::newFull(const expert_forest *f, int level, unsigned tsz)
{
  unpacked_node* U = useUnpackedNode();
  MEDDLY_DCASSERT(U);
  U->initFull(f, level, tsz);
  U->clearFullEdges();
  addToBuildList(U);
  return U;
}

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::newSparse(const expert_forest *f, int level, unsigned nnzs)
{
  unpacked_node* U = useUnpackedNode();
  MEDDLY_DCASSERT(U);
  U->initSparse(f, level, nnzs);
  U->clearSparseEdges();
  addToBuildList(U);
  return U;
}

// ****************************************************************************

inline const void*
MEDDLY::unpacked_node::UHptr() const
{
  MEDDLY_DCASSERT(extra_unhashed);
  return extra_unhashed;
}

inline void*
MEDDLY::unpacked_node::UHdata()
{
  MEDDLY_DCASSERT(extra_unhashed);
  return extra_unhashed;
}

inline unsigned
MEDDLY::unpacked_node::UHbytes() const
{
  return ext_uh_size;
}

inline const void*
MEDDLY::unpacked_node::HHptr() const
{
  MEDDLY_DCASSERT(extra_hashed);
  return extra_hashed;
}

inline void*
MEDDLY::unpacked_node::HHdata()
{
  MEDDLY_DCASSERT(extra_hashed);
  return extra_hashed;
}

inline unsigned
MEDDLY::unpacked_node::HHbytes() const
{
  return ext_h_size;
}

inline MEDDLY::node_handle
MEDDLY::unpacked_node::d(unsigned n) const
{
  MEDDLY_DCASSERT(down);
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, (is_full ? size : nnzs));
  return down[n];
}

inline MEDDLY::node_handle&
MEDDLY::unpacked_node::d_ref(unsigned n)
{
  MEDDLY_DCASSERT(down);
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, (is_full ? size : nnzs));
  return down[n];
}

inline void
MEDDLY::unpacked_node::set_d(unsigned n, dd_edge &E)
{
  MEDDLY_DCASSERT(parent == E.parent);
  MEDDLY_DCASSERT(0==edge_bytes);
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, (is_full ? size : nnzs));
  down[n] = E.node;
  E.node = 0; // avoid having to adjust the link count
}

inline unsigned
MEDDLY::unpacked_node::i(unsigned n) const
{
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(!is_full);
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, nnzs);
  return index[n];
}

inline unsigned&
MEDDLY::unpacked_node::i_ref(unsigned n)
{
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(!is_full);
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, nnzs);
  return index[n];
}

inline const void*
MEDDLY::unpacked_node::eptr(unsigned i) const
{
  MEDDLY_DCASSERT(edge);
  MEDDLY_DCASSERT(edge_bytes);
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, (is_full ? size : nnzs));
  return ((char*) edge) + i * edge_bytes;
}

inline void*
MEDDLY::unpacked_node::eptr_write(unsigned i)
{
  MEDDLY_DCASSERT(edge);
  MEDDLY_DCASSERT(edge_bytes);
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, (is_full ? size : nnzs));
  return ((char*) edge) + i * edge_bytes;
}

inline void
MEDDLY::unpacked_node::set_de(unsigned n, dd_edge &E)
{
  MEDDLY_DCASSERT(parent == E.parent);
  MEDDLY_DCASSERT(edge);
  MEDDLY_DCASSERT(edge_bytes);
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, n, (is_full ? size : nnzs));
  down[n] = E.node;
  memcpy( ((char*) edge) + n * edge_bytes, & (E.raw_value), edge_bytes );
  E.node = 0; // avoid having to adjust the link count
}

inline void
MEDDLY::unpacked_node::getEdge(unsigned n, long &val) const
{
  MEDDLY_DCASSERT(sizeof(long) == edge_bytes);
  MEDDLY::EVencoder<long>::readValue(eptr(n), val);
}

inline void
MEDDLY::unpacked_node::getEdge(unsigned n, float &val) const
{
  MEDDLY_DCASSERT(sizeof(float) == edge_bytes);
  MEDDLY::EVencoder<float>::readValue(eptr(n), val);
}

inline void
MEDDLY::unpacked_node::setEdge(unsigned n, long ev)
{
  MEDDLY_DCASSERT(sizeof(long) == edge_bytes);
  MEDDLY::EVencoder<long>::writeValue(eptr_write(n), ev);

  long test_ev = 256;
  getEdge(n, test_ev);
  MEDDLY_DCASSERT(test_ev == ev);
}


inline void
MEDDLY::unpacked_node::setEdge(unsigned n, float ev)
{
  MEDDLY_DCASSERT(sizeof(float) == edge_bytes);
  MEDDLY::EVencoder<float>::writeValue(eptr_write(n), ev);
}

inline long
MEDDLY::unpacked_node::ei(unsigned i) const
{
  long ev;
  getEdge(i, ev);
  return ev;
}

inline float
MEDDLY::unpacked_node::ef(unsigned i) const
{
  float ev;
  getEdge(i, ev);
  return ev;
}

inline int
MEDDLY::unpacked_node::getLevel() const
{
  return level;
}

inline void
MEDDLY::unpacked_node::setLevel(int k)
{
  level = k;
}

// ---------------------------------------------
// Extensible portion of the node
// ---------------------------------------------

inline bool
MEDDLY::unpacked_node::isExtensible() const
{
  return is_extensible;
}

inline void
MEDDLY::unpacked_node::markAsExtensible()
{
  MEDDLY_DCASSERT(parent->isExtensibleLevel(getLevel()));
  is_extensible = true;
}

inline void
MEDDLY::unpacked_node::markAsNotExtensible()
{
  is_extensible = false;
}

inline unsigned
MEDDLY::unpacked_node::ext_i() const
{
  MEDDLY_DCASSERT(isExtensible());
  return isSparse()? i(getNNZs() - 1): getSize() - 1;
}

inline MEDDLY::node_handle
MEDDLY::unpacked_node::ext_d() const
{
  MEDDLY_DCASSERT(isExtensible());
  return d( (isSparse()? getNNZs(): getSize()) - 1 );
}

inline long
MEDDLY::unpacked_node::ext_ei() const
{
  MEDDLY_DCASSERT(isExtensible());
  return ei( (isSparse()? getNNZs(): getSize()) - 1 );
}

inline float
MEDDLY::unpacked_node::ext_ef() const
{
  MEDDLY_DCASSERT(isExtensible());
  return ef( (isSparse()? getNNZs(): getSize()) - 1 );
}

// --- End of Extensible portion of the node ---

inline unsigned
MEDDLY::unpacked_node::getSize() const
{
  MEDDLY_DCASSERT(is_full);
  return size;
}

inline unsigned
MEDDLY::unpacked_node::getNNZs() const
{
  MEDDLY_DCASSERT(!is_full);
  return nnzs;
}

inline bool
MEDDLY::unpacked_node::isSparse() const
{
  return !is_full;
}

inline bool
MEDDLY::unpacked_node::isFull() const
{
  return is_full;
}

inline bool
MEDDLY::unpacked_node::hasEdges() const
{
  return edge_bytes;
}
inline unsigned
MEDDLY::unpacked_node::edgeBytes() const
{
  return edge_bytes;
}

inline unsigned
MEDDLY::unpacked_node::hash() const
{
#ifdef DEVELOPMENT_CODE
  MEDDLY_DCASSERT(has_hash);
#endif
  return h;
}

inline void
MEDDLY::unpacked_node::setHash(unsigned H)
{
#ifdef DEVELOPMENT_CODE
  MEDDLY_DCASSERT(!has_hash);
  has_hash = true;
#endif
  h = H;
}

inline bool
MEDDLY::unpacked_node::isBuildNode() const
{
  return is_in_build_list;
}

inline void
MEDDLY::unpacked_node::shrinkFull(unsigned ns)
{
  MEDDLY_DCASSERT(isFull());
  MEDDLY_DCASSERT(ns >= 0);
  MEDDLY_DCASSERT(ns <= size);
  size = ns;
}

inline void
MEDDLY::unpacked_node::shrinkSparse(unsigned ns)
{
  MEDDLY_DCASSERT(isSparse());
  MEDDLY_DCASSERT(ns >= 0);
  MEDDLY_DCASSERT(ns <= nnzs);
  nnzs = ns;
}

inline void
MEDDLY::unpacked_node::bind_as_full(bool full)
{
  is_full = full;
}

inline void
MEDDLY::unpacked_node::clearFullEdges()
{
  MEDDLY_DCASSERT(isFull());
  memset(down, 0, unsigned(size) * sizeof(node_handle));
  if (edge_bytes) {
    memset(edge, 0, unsigned(size) * edge_bytes);
  }
}

inline void
MEDDLY::unpacked_node::clearSparseEdges()
{
  MEDDLY_DCASSERT(isSparse());
  memset(down, 0, unsigned(nnzs) * sizeof(node_handle));
  if (edge_bytes) {
    memset(edge, 0, unsigned(nnzs) * edge_bytes);
  }
}

inline MEDDLY::unpacked_node*
MEDDLY::unpacked_node::useUnpackedNode()
{
  unpacked_node* nr;
  if (freeList) {
    nr = freeList;
    freeList = nr->next;
  }
  else {
    nr = new unpacked_node;
  }
#ifdef DEVELOPMENT_CODE
  nr->has_hash = false;
#endif
  nr->is_in_build_list = false;
  return nr;
}

inline void
MEDDLY::unpacked_node::recycle(MEDDLY::unpacked_node* r)
{
  if (r) {
    if (r->is_in_build_list) {
      removeFromBuildList(r);
    }
    r->next = freeList;
    freeList = r;
  }
}

inline void
MEDDLY::unpacked_node::freeRecycled()
{
  while (freeList) {
    MEDDLY::unpacked_node* n = freeList->next;
    delete freeList;
    freeList = n;
  }
}

inline void
MEDDLY::unpacked_node::addToBuildList(unpacked_node* b)
{
#ifdef DEBUG_BUILDLIST
  printf("Adding unpacked node at level %d to build list\n", b->getLevel());
#endif
  MEDDLY_DCASSERT(b);
  b->is_in_build_list = true;
  b->next = buildList;
  buildList = b;
}





#endif

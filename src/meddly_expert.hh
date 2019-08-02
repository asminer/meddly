
/*
 Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
 Copyright (C) 2009, Iowa State University Research Foundation, Inc.

 This library is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published
 by the Free Software Foundation, either version 3 of the License, orf
 (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library.  If not, see <http://www.gnu.org/licenses/>.
 */

/*! \file meddly_expert.hh

 Implementation details for interface in meddly_expert.h.
 */

#ifndef MEDDLY_EXPERT_HH
#define MEDDLY_EXPERT_HH

// ******************************************************************
// *                                                                *
// *                 inlined  expert_domain methods                 *
// *                                                                *
// ******************************************************************

inline MEDDLY::expert_variable*
MEDDLY::expert_domain::getExpertVar(int lev) const
{
  return (expert_variable*) vars[lev];
}

inline const MEDDLY::expert_variable*
MEDDLY::expert_domain::readExpertVar(int lev) const
{
  return (expert_variable*) vars[lev];
}

inline void
MEDDLY::expert_domain::enlargeVariableBound(int vh, bool prime, int b)
{
  getExpertVar(vh)->enlargeBound(prime, b);
}

inline void
MEDDLY::expert_domain::shrinkVariableBound(int vh, int b, bool force)
{
  getExpertVar(vh)->shrinkBound(b, force);
}


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
  MEDDLY_CHECK_RANGE(0, n, (is_full ? size : nnzs));
  return down[n];
}

inline MEDDLY::node_handle&
MEDDLY::unpacked_node::d_ref(unsigned n) 
{
  MEDDLY_DCASSERT(down);
  MEDDLY_CHECK_RANGE(0, n, (is_full ? size : nnzs));
  return down[n];
}

inline void
MEDDLY::unpacked_node::set_d(unsigned n, dd_edge &E)
{
  MEDDLY_DCASSERT(parent == E.parent);
  MEDDLY_DCASSERT(0==edge_bytes);
  MEDDLY_CHECK_RANGE(0, n, (is_full ? size : nnzs));
  down[n] = E.node;
  E.node = 0; // avoid having to adjust the link count
}

inline unsigned
MEDDLY::unpacked_node::i(unsigned n) const
{
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(!is_full);
  MEDDLY_CHECK_RANGE(0, n, nnzs);
  return index[n];
}

inline unsigned&
MEDDLY::unpacked_node::i_ref(unsigned n)
{
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(!is_full);
  MEDDLY_CHECK_RANGE(0, n, nnzs);
  return index[n];
}

inline const void*
MEDDLY::unpacked_node::eptr(unsigned i) const
{
  MEDDLY_DCASSERT(edge);
  MEDDLY_DCASSERT(edge_bytes);
  MEDDLY_CHECK_RANGE(0, i, (is_full ? size : nnzs));
  return ((char*) edge) + i * edge_bytes;
}

inline void*
MEDDLY::unpacked_node::eptr_write(unsigned i)
{
  MEDDLY_DCASSERT(edge);
  MEDDLY_DCASSERT(edge_bytes);
  MEDDLY_CHECK_RANGE(0, i, (is_full ? size : nnzs));
  return ((char*) edge) + i * edge_bytes;
}

inline void
MEDDLY::unpacked_node::set_de(unsigned n, dd_edge &E)
{
  MEDDLY_DCASSERT(parent == E.parent);
  MEDDLY_DCASSERT(edge);
  MEDDLY_DCASSERT(edge_bytes);
  MEDDLY_CHECK_RANGE(0, n, (is_full ? size : nnzs));
  down[n] = E.node;
  memcpy( ((char*) edge) + n * edge_bytes, & (E.raw_value), edge_bytes );
  E.node = 0; // avoid having to adjust the link count
}

inline void
MEDDLY::unpacked_node::getEdge(unsigned n, long &val) const
{
  MEDDLY_DCASSERT(sizeof(long) == edge_bytes);
  MEDDLY::expert_forest::EVencoder<long>::readValue(eptr(n), val);
}

inline void
MEDDLY::unpacked_node::getEdge(unsigned n, float &val) const
{
  MEDDLY_DCASSERT(sizeof(float) == edge_bytes);
  MEDDLY::expert_forest::EVencoder<float>::readValue(eptr(n), val);
}

inline void
MEDDLY::unpacked_node::setEdge(unsigned n, long ev)
{
  MEDDLY_DCASSERT(sizeof(long) == edge_bytes);
  MEDDLY::expert_forest::EVencoder<long>::writeValue(eptr_write(n), ev);

  long test_ev = 256;
  getEdge(n, test_ev);
  MEDDLY_DCASSERT(test_ev == ev);
}


inline void
MEDDLY::unpacked_node::setEdge(unsigned n, float ev)
{
  MEDDLY_DCASSERT(sizeof(float) == edge_bytes);
  MEDDLY::expert_forest::EVencoder<float>::writeValue(eptr_write(n), ev);
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
  printf("\n Setting hash as %ul",H);
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

// ******************************************************************
// *                                                                *
// *           inlined node_headers::level_array  methods           *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

inline int MEDDLY::node_headers::level_array::get(size_t i) const
{
  MEDDLY_DCASSERT(i<size);
  if (data8) {
    MEDDLY_DCASSERT(0==data16);
    MEDDLY_DCASSERT(0==data32);
    return data8[i];
  }
  if (data16) {
    MEDDLY_DCASSERT(0==data8);
    MEDDLY_DCASSERT(0==data32);
    return data16[i];
  }
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data16);
  MEDDLY_DCASSERT(data32);
  return data32[i];
}

inline void MEDDLY::node_headers::level_array::set(size_t i, int v)
{
  MEDDLY_DCASSERT(i<size);
  if (data8) {
    MEDDLY_DCASSERT(0==data16);
    MEDDLY_DCASSERT(0==data32);
    MEDDLY_DCASSERT(v>-128);
    MEDDLY_DCASSERT(v<128);
    data8[i] = v;
    return;
  }
  if (data16) {
    MEDDLY_DCASSERT(0==data8);
    MEDDLY_DCASSERT(0==data32);
    MEDDLY_DCASSERT(v>-32768);
    MEDDLY_DCASSERT(v<32768);
    data16[i] = v;
    return;
  }
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data16);
  MEDDLY_DCASSERT(data32);
  data32[i] = v;
}

inline void MEDDLY::node_headers::level_array::swap(size_t i, size_t j)
{
  MEDDLY_DCASSERT(i<size);
  MEDDLY_DCASSERT(j<size);
  if (data8) {
    char tmp = data8[i];
    data8[i] = data8[j];
    data8[j] = tmp;
    return;
  }
  if (data16) {
    short tmp = data16[i];
    data16[i] = data16[j];
    data16[j] = tmp;
    return;
  }
  MEDDLY_DCASSERT(data32);
  int tmp = data32[i];
  data32[i] = data32[j];
  data32[j] = tmp;
}

inline size_t MEDDLY::node_headers::level_array::entry_bits() const
{
  return size_t(bytes) * 8;
}

#endif

// ******************************************************************
// *                                                                *
// *          inlined node_headers::counter_array  methods          *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

inline unsigned int MEDDLY::node_headers::counter_array::get(size_t i) const
{
  MEDDLY_DCASSERT(i<size);
  if (data8) {
    MEDDLY_DCASSERT(0==data16);
    MEDDLY_DCASSERT(0==data32);
    return data8[i];
  }
  if (data16) {
    MEDDLY_DCASSERT(0==data8);
    MEDDLY_DCASSERT(0==data32);
    return data16[i];
  }
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data16);
  MEDDLY_DCASSERT(data32);
  return data32[i];
}

inline void MEDDLY::node_headers::counter_array::swap(size_t i, size_t j)
{
  MEDDLY_DCASSERT(i<size);
  MEDDLY_DCASSERT(j<size);
  if (data8) {
    unsigned char tmp = data8[i];
    data8[i] = data8[j];
    data8[j] = tmp;
    return;
  }
  if (data16) {
    unsigned short tmp = data16[i];
    data16[i] = data16[j];
    data16[j] = tmp;
    return;
  }
  MEDDLY_DCASSERT(data32);
  unsigned int tmp = data32[i];
  data32[i] = data32[j];
  data32[j] = tmp;
}

inline void MEDDLY::node_headers::counter_array::increment(size_t i)
{
  MEDDLY_DCASSERT(i<size);
  if (data8) {
    MEDDLY_DCASSERT(0==data16);
    MEDDLY_DCASSERT(0==data32);
    if (0 == ++data8[i]) expand8to16(i);
    return;
  }
  if (data16) {
    MEDDLY_DCASSERT(0==data8);
    MEDDLY_DCASSERT(0==data32);
    ++data16[i];
    if (256 == data16[i]) ++counts_09bit;
    if (0 == data16[i]) expand16to32(i);
    return;
  }
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data16);
  MEDDLY_DCASSERT(data32);
  data32[i]++;
  if (256 == data32[i]) ++counts_09bit;
  if (65536 == data32[i]) ++counts_17bit;
}

inline void MEDDLY::node_headers::counter_array::decrement(size_t i)
{
  MEDDLY_DCASSERT(i<size);
  if (data8) {
    MEDDLY_DCASSERT(0==data16);
    MEDDLY_DCASSERT(0==data32);
    MEDDLY_DCASSERT(data8[i]);
    --data8[i];
    return;
  }
  if (data16) {
    MEDDLY_DCASSERT(0==data8);
    MEDDLY_DCASSERT(0==data32);
    MEDDLY_DCASSERT(data16[i]);
    if (256 == data16[i]) {
      MEDDLY_DCASSERT(counts_09bit);
      --counts_09bit;
    }
    --data16[i];
    return;
  }
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data16);
  MEDDLY_DCASSERT(data32);
  MEDDLY_DCASSERT(data32[i]);
  if (256 == data32[i]) {
    MEDDLY_DCASSERT(counts_09bit);
    --counts_09bit;
  }
  if (65536 == data32[i]) {
    MEDDLY_DCASSERT(counts_17bit);
    --counts_17bit;
  }
  --data32[i];
}

inline bool MEDDLY::node_headers::counter_array::isZeroBeforeIncrement(size_t i)
{
  MEDDLY_DCASSERT(i<size);
  if (data8) {
    MEDDLY_DCASSERT(0==data16);
    MEDDLY_DCASSERT(0==data32);
    if (0==data8[i]) {
      data8[i] = 1;
      return true;
    }
    if (0 == ++data8[i]) expand8to16(i);
    return false;
  }
  if (data16) {
    MEDDLY_DCASSERT(0==data8);
    MEDDLY_DCASSERT(0==data32);
    if (0==data16[i]) {
      data16[i] = 1;
      return true;
    }
    ++data16[i];
    if (256 == data16[i]) ++counts_09bit;
    if (0 == data16[i]) expand16to32(i);
    return false;
  }
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data16);
  MEDDLY_DCASSERT(data32);
  if (0==data32[i]) {
    data32[i] = 1;
    return true;
  }
  data32[i]++;
  if (256 == data32[i]) ++counts_09bit;
  if (65536 == data32[i]) ++counts_17bit;
  return false;
}

inline bool MEDDLY::node_headers::counter_array::isPositiveAfterDecrement(size_t i)
{
  MEDDLY_DCASSERT(i<size);
  if (data8) {
    MEDDLY_DCASSERT(0==data16);
    MEDDLY_DCASSERT(0==data32);
    MEDDLY_DCASSERT(data8[i]);
    return 0 < --data8[i];
  }
  if (data16) {
    MEDDLY_DCASSERT(0==data8);
    MEDDLY_DCASSERT(0==data32);
    MEDDLY_DCASSERT(data16[i]);
    if (256 == data16[i]) {
      MEDDLY_DCASSERT(counts_09bit);
      --counts_09bit;
    }
    return 0 < --data16[i];
  }
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data16);
  MEDDLY_DCASSERT(data32);
  MEDDLY_DCASSERT(data32[i]);
  if (256 == data32[i]) {
    MEDDLY_DCASSERT(counts_09bit);
    --counts_09bit;
  }
  if (65536 == data32[i]) {
    MEDDLY_DCASSERT(counts_17bit);
    --counts_17bit;
  }
  return 0<--data32[i];
}

inline size_t MEDDLY::node_headers::counter_array::entry_bits() const
{
  return bytes * 8;
}


#endif

// ******************************************************************
// *                                                                *
// *          inlined node_headers::address_array  methods          *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

inline unsigned long MEDDLY::node_headers::address_array::get(size_t i) const
{
  MEDDLY_DCASSERT(i<size);
  if (4==bytes) {
    MEDDLY_DCASSERT(data32);
    MEDDLY_DCASSERT(0==data64);
    MEDDLY_DCASSERT(0==num_large_elements);
    return data32[i];
  }
  MEDDLY_DCASSERT(0==data32);
  MEDDLY_DCASSERT(8==bytes);
  MEDDLY_DCASSERT(data64);
  return data64[i];
}

inline void MEDDLY::node_headers::address_array::set(size_t i, unsigned long v)
{
  MEDDLY_DCASSERT(i<size);
  if (4==bytes) {
    MEDDLY_DCASSERT(data32);
    MEDDLY_DCASSERT(0==data64);
    MEDDLY_DCASSERT(0==num_large_elements);

    if (v & 0xffffffff00000000) {
      // v won't fit in 32 bits
      expand32to64();
      MEDDLY_DCASSERT(data64);
      data64[i] = v;
    } else {
      // v will fit in 32 bits
      data32[i] = v;
    }
    return;
  }
  MEDDLY_DCASSERT(0==data32);
  MEDDLY_DCASSERT(8==bytes);
  MEDDLY_DCASSERT(data64);
  if (v & 0xffffffff00000000) {
    // v is large
    if (0 == (data64[i] & 0xffffffff00000000)) {
      // replacing small
      num_large_elements++;
    }
  } else {
    // v is small
    if (data64[i] & 0xffffffff00000000) {
      // replacing large
      MEDDLY_DCASSERT(num_large_elements);
      num_large_elements--;
    }
  }
  data64[i] = v;
}

inline void MEDDLY::node_headers::address_array::swap(size_t i, size_t j)
{
  MEDDLY_DCASSERT(i<size);
  MEDDLY_DCASSERT(j<size);
  if (4==bytes) {
    MEDDLY_DCASSERT(data32);
    unsigned int tmp = data32[i];
    data32[i] = data32[j];
    data32[j] = tmp;
    return;
  }
  MEDDLY_DCASSERT(8==bytes);
  MEDDLY_DCASSERT(data64);
  unsigned long tmp = data64[i];
  data64[i] = data64[j];
  data64[j] = tmp;
}

inline size_t MEDDLY::node_headers::address_array::entry_bits() const
{
  return bytes * 8;
}

#endif

// ******************************************************************
// *                                                                *
// *            inlined node_headers::bitvector  methods            *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

inline bool MEDDLY::node_headers::bitvector::get(size_t i) const
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(i<size);
  return data[i];
}

inline void MEDDLY::node_headers::bitvector::set(size_t i, bool v)
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(i<size);
  data[i] = v;
}

inline void MEDDLY::node_headers::bitvector::clearAll()
{
  if (size) {
    MEDDLY_DCASSERT(data);
    memset(data, 0, size * sizeof(bool));
  }
}

inline void MEDDLY::node_headers::bitvector::swap(size_t i, size_t j)
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(i<size);
  MEDDLY_DCASSERT(j<size);
  bool tmp = data[i];
  data[i] = data[j];
  data[j] = tmp;
}

inline size_t MEDDLY::node_headers::bitvector::entry_bits() const
{
  return sizeof(bool) * 8;
}

inline size_t MEDDLY::node_headers::bitvector::firstZero(size_t start) const
{
  for (; start < size; start++) {
    if (0==data[start]) return start;
  }
  return size;
}

#endif

// ******************************************************************
// *                                                                *
// *                  inlined node_headers methods                  *
// *                                                                *
// ******************************************************************

inline void MEDDLY::node_headers::setPessimistic(bool pess)
{
  pessimistic = pess;
}

// ******************************************************************

inline MEDDLY::node_handle
MEDDLY::node_headers::lastUsedHandle() const
{
  return a_last;
}

// ******************************************************************

inline bool
MEDDLY::node_headers::isActive(node_handle p) const
{
  return !isDeleted(p);
}

// ******************************************************************

/*
inline bool
MEDDLY::node_headers::isZombie(node_handle p) const
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  return (0==address[p].offset) && (0!=address[p].level);
#else
  MEDDLY_DCASSERT(addresses);
  MEDDLY_DCASSERT(levels);
  return (0==addresses->get(size_t(p))) && (0!=levels->get(size_t(p)));
#endif
}
*/

// ******************************************************************

inline bool
MEDDLY::node_headers::isDeleted(node_handle p) const
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  return (0==address[p].level);
#else
  MEDDLY_DCASSERT(levels);
  return (0==levels->get(size_t(p)));
#endif
}

// ******************************************************************

/*
inline bool
MEDDLY::node_headers::isDeactivated(node_handle p) const
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  return (0==address[p].level);
#else
  MEDDLY_DCASSERT(levels);
  return (0==levels->get(size_t(p)));
#endif
}
*/

// ******************************************************************

inline void
MEDDLY::node_headers::deactivate(node_handle p)
{
  MEDDLY_DCASSERT(isActive(p));
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  address[p].level = 0;
#else
  MEDDLY_DCASSERT(levels);
  levels->set(size_t(p), 0);
#endif
}


// ******************************************************************

inline MEDDLY::node_address 
MEDDLY::node_headers::getNodeAddress(node_handle p) const
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  return address[p].offset;
#else
  MEDDLY_DCASSERT(addresses);
  return addresses->get(size_t(p));
#endif
}

// ******************************************************************

inline void MEDDLY::node_headers::setNodeAddress(node_handle p, node_address a)
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  address[p].offset = a;
#else
  MEDDLY_DCASSERT(addresses);
  addresses->set(size_t(p), a);
#endif
}

// ******************************************************************

inline void MEDDLY::node_headers::moveNodeAddress(node_handle p, 
  node_address old_addr, node_address new_addr)
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  MEDDLY_DCASSERT(old_addr == address[p].offset);
  address[p].offset = new_addr;
#else
  MEDDLY_DCASSERT(addresses);
  MEDDLY_DCASSERT(old_addr == addresses->get(size_t(p)));
  addresses->set(size_t(p), new_addr);
#endif
}

// ******************************************************************

inline int
MEDDLY::node_headers::getNodeLevel(node_handle p) const
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  return address[p].level;
#else
  MEDDLY_DCASSERT(levels);
  return levels->get(size_t(p));
#endif
}

// ******************************************************************

inline void
MEDDLY::node_headers::setNodeLevel(node_handle p, int k)
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  address[p].level = k;
#else
  MEDDLY_DCASSERT(levels);
  levels->set(size_t(p), k);
#endif
}

// ******************************************************************

/*
inline bool
MEDDLY::node_headers::trackingCacheCounts() const
{
#ifdef OLD_NODE_HEADERS
  return usesCacheCounts;
#else
  return cache_counts;
#endif
}
*/

// ******************************************************************

inline unsigned long
MEDDLY::node_headers::getNodeCacheCount(node_handle p) const
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
#ifdef OLD_NODE_HEADERS
  // MEDDLY_DCASSERT(usesCacheCounts); // or do we just return 0?  TBD
  MEDDLY_DCASSERT(address);
  return address[p].cache_count;
#else
  MEDDLY_DCASSERT(cache_counts);
  return cache_counts->get(size_t(p));
#endif
}

// ******************************************************************

inline void
MEDDLY::node_headers::cacheNode(node_handle p)
{
  if (p<1) return;    // terminal node

  MEDDLY_DCASSERT(isActive(p));
#ifdef OLD_NODE_HEADERS

  MEDDLY_DCASSERT(address);
  address[p].cache_count++;
#ifdef TRACK_CACHECOUNT
  fprintf(stdout, "\t+Node %d is in %d caches\n", p, address[p].cache_count);
  fflush(stdout);
#endif

#else

  MEDDLY_DCASSERT(cache_counts);
  cache_counts->increment(size_t(p));
#ifdef TRACK_CACHECOUNT
  fprintf(stdout, "\t+Node %d is in %u caches\n", p, cache_counts->get(size_t(p)));
  fflush(stdout);
#endif

#endif
}

// ******************************************************************

inline void
MEDDLY::node_headers::uncacheNode(MEDDLY::node_handle p)
{
  if (p<1) return;
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);

#ifdef OLD_NODE_HEADERS

  MEDDLY_DCASSERT(address);
  MEDDLY_DCASSERT(address[p].cache_count > 0);
  address[p].cache_count--;

#ifdef TRACK_CACHECOUNT
  fprintf(stdout, "\t-Node %d is in %d caches\n", p, address[p].cache_count);
  fflush(stdout);
#endif

  if (address[p].cache_count) return;

  //
  // Still here?  Might need to clean up
  //
  
  if (isDeleted(p)) {
      // we were already disconnected; must be using pessimistic
#ifdef TRACK_UNREACHABLE_NODES
      parent.stats.unreachable_nodes--;
#endif
      recycleNodeHandle(p);
  } else {
      // We're still active
      // See if we're now completely disconnected
      // and if so, tell parent to recycle node storage
      if (0==address[p].incoming_count) {
#ifdef TRACK_UNREACHABLE_NODES
        parent.stats.unreachable_nodes--;
#endif

        parent.deleteNode(p);
        recycleNodeHandle(p);
      }
  }

#else

  MEDDLY_DCASSERT(addresses);
  MEDDLY_DCASSERT(cache_counts);

  if (cache_counts->isPositiveAfterDecrement(size_t(p))) {
#ifdef TRACK_CACHECOUNT
    fprintf(stdout, "\t-Node %d is in %lu caches\n", p, cache_counts->get(size_t(p)));
    fflush(stdout);
#endif
    return;
  }

#ifdef TRACK_CACHECOUNT
  fprintf(stdout, "\t-Node %d is not in any caches\n", p);
  fflush(stdout);
#endif

  //
  // Still here?  Cache count now zero; might need to clean up
  //

  if (isDeleted(p)) {
    //
    // Must be using pessimistic
    //
#ifdef TRACK_UNREACHABLE_NODES
    parent.stats.unreachable_nodes--;
#endif
    recycleNodeHandle(p);
  } else {
    // we were active.
    // See if we're now completely disconnected
    // and if so, tell parent to recycle node storage
    if (
          (incoming_counts && (0==incoming_counts->get(size_t(p)))) 
          ||
          (is_reachable && (0==is_reachable->get(size_t(p))))
    ) {
#ifdef TRACK_UNREACHABLE_NODES
      parent.stats.unreachable_nodes--;
#endif
      parent.deleteNode(p);
      recycleNodeHandle(p);
    }
  }

#endif
}

// ******************************************************************

inline void
MEDDLY::node_headers::setInCacheBit(node_handle p)
{
  if (p<1) return;    // terminal node
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);

#ifdef OLD_NODE_HEADERS

  MEDDLY_DCASSERT(address);
  address[p].cache_count = 1;

#else

  MEDDLY_DCASSERT(is_in_cache);
  is_in_cache->set(size_t(p), 1);

#endif
}

// ******************************************************************

inline void
MEDDLY::node_headers::clearAllInCacheBits()
{
#ifdef OLD_NODE_HEADERS

  MEDDLY_DCASSERT(address || (0==a_last));

  for (node_handle p=1; p<=a_last; p++) {
    address[p].cache_count = 0;
  }

#else

  MEDDLY_DCASSERT(is_in_cache);
  is_in_cache->clearAll();

#endif
}

// ******************************************************************

/*
inline bool
MEDDLY::node_headers::trackingIncomingCounts() const
{
#ifdef OLD_NODE_HEADERS
  return usesIncomingCounts;
#else
  return incoming_counts;
#endif
}
*/

// ******************************************************************

inline unsigned long 
MEDDLY::node_headers::getIncomingCount(node_handle p) const
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
#ifdef OLD_NODE_HEADERS
  // MEDDLY_DCASSERT(usesIncomingCounts); // or do we just return 0?  TBD
  MEDDLY_DCASSERT(address);
  return address[p].incoming_count;
#else
  MEDDLY_DCASSERT(incoming_counts);
  return incoming_counts->get(size_t(p));
#endif
}

// ******************************************************************

inline MEDDLY::node_handle
MEDDLY::node_headers::linkNode(node_handle p)
{
  if (p<1) return p;    // terminal node

  MEDDLY_DCASSERT(isActive(p));

#ifdef OLD_NODE_HEADERS

  // MEDDLY_DCASSERT(usesIncomingCounts); // or do we just return?  TBD

  MEDDLY_DCASSERT(address);
  MEDDLY_DCASSERT(address[p].offset);

  if (0==address[p].incoming_count) {
    // Reclaim an unreachable
    parent.stats.reclaimed_nodes++;
#ifdef TRACK_UNREACHABLE_NODES
    parent.stats.unreachable_nodes--;
#endif
  }

  address[p].incoming_count++;

#ifdef TRACK_DELETIONS
  fprintf(stdout, "\t+Node %d count now %ld\n", p, address[p].incoming_count);
  fflush(stdout);
#endif

#else // OLD_NODE_HEADERS

  MEDDLY_DCASSERT(incoming_counts);
  if (incoming_counts->isZeroBeforeIncrement(size_t(p))) {
    // Reclaim an unreachable 
    parent.stats.reclaimed_nodes++;
#ifdef TRACK_UNREACHABLE_NODES
    parent.stats.unreachable_nodes--;
#endif
  }

#ifdef TRACK_DELETIONS
  fprintf(stdout, "\t+Node %d count now %ld\n", p, incoming_counts->get(size_t(p)));
  fflush(stdout);
#endif

#endif // OLD_NODE_HEADERS

  return p;
}

// ******************************************************************

inline void
MEDDLY::node_headers::unlinkNode(node_handle p)
{
  if (p<1) return;    // terminal node

  MEDDLY_DCASSERT(isActive(p));

#ifdef OLD_NODE_HEADERS

  // MEDDLY_DCASSERT(usesIncomingCounts); // or do we just return?  TBD
  MEDDLY_DCASSERT(address);
  MEDDLY_DCASSERT(address[p].offset);
  MEDDLY_DCASSERT(address[p].incoming_count>0);

  address[p].incoming_count--;

#ifdef TRACK_DELETIONS
  fprintf(stdout, "\t+Node %d count now %ld\n", p, address[p].incoming_count);
  fflush(stdout);
#endif

  if (address[p].incoming_count) return;

  //
  // Node p just became unreachable.
  // Various cleanups below.
  //

  //
  // If we're not in any caches, delete
  //
  if (0==address[p].cache_count) {
        parent.deleteNode(p);
        recycleNodeHandle(p);
        return;
  }

  //
  // Still in some cache somewhere.
  //
  
  if (pessimistic) {
    parent.deleteNode(p);
  } else {
#ifdef TRACK_UNREACHABLE_NODES
    parent.stats.unreachable_nodes++;
#endif
  }


#else // OLD_NODE_HEADERS

  MEDDLY_DCASSERT(incoming_counts);
  MEDDLY_DCASSERT(addresses);

  if (incoming_counts->isPositiveAfterDecrement(size_t(p))) {
#ifdef TRACK_DELETIONS
    fprintf(stdout, "\t+Node %d count now %ld\n", p, incoming_counts->get(size_t(p)));
    fflush(stdout);
#endif
    return;
  }

#ifdef TRACK_DELETIONS
  fprintf(stdout, "\t+Node %d count now zero\n", p);
  fflush(stdout);
#endif

  //
  // Node is unreachable.  See if we need to do some cleanup
  //

  if (
        (cache_counts && (0==cache_counts->get(size_t(p)))) 
        ||
        (is_in_cache && (0==is_in_cache->get(size_t(p))))
  ) {
    //
    // Unreachable and not in any caches.  Delete and recycle handle.
    //
    parent.deleteNode(p);
    recycleNodeHandle(p);
    return;
  }
  
  if (pessimistic) {
    //
    // Delete; keep handle until caches are cleared.
    //
    parent.deleteNode(p);
  } else {
    //
    // Optimistic.  Keep unreachables around.
    //
#ifdef TRACK_UNREACHABLE_NODES
    parent.stats.unreachable_nodes++;
#endif
  }

#endif // OLD_NODE_HEADERS
}

// ******************************************************************

inline void
MEDDLY::node_headers::setReachableBit(node_handle p)
{
  if (p<1) return;    // terminal node
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);

  MEDDLY_DCASSERT(isActive(p));

#ifdef OLD_NODE_HEADERS

  MEDDLY_DCASSERT(address);
  address[p].incoming_count = 1;

#else

  MEDDLY_DCASSERT(is_reachable);
  is_reachable->set(size_t(p), 1);

#endif
}

// ******************************************************************

inline bool
MEDDLY::node_headers::hasReachableBit(node_handle p) const
{
  if (p<1) return 1;    // terminal node
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);

#ifdef OLD_NODE_HEADERS

  MEDDLY_DCASSERT(address);
  return address[p].incoming_count;

#else

  MEDDLY_DCASSERT(is_reachable);
  return is_reachable->get(size_t(p));

#endif
}

// ******************************************************************

inline void
MEDDLY::node_headers::clearAllReachableBits()
{
#ifdef OLD_NODE_HEADERS

  MEDDLY_DCASSERT(address || (0==a_last));
  for (node_handle p=1; p<=a_last; p++) {
    address[p].incoming_count = 0;
  }

#else

  MEDDLY_DCASSERT(is_reachable);
  is_reachable->clearAll();

#endif
}

// ******************************************************************

inline int 
MEDDLY::node_headers::getNodeImplicitFlag(node_handle p) const
{
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  return address[p].is_implicit;
#else
  return (implicit_bits) ? implicit_bits->get(size_t(p)) : false;
#endif
}

// ******************************************************************

inline void 
MEDDLY::node_headers::setNodeImplicitFlag(node_handle p, bool flag)
{
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  address[p].is_implicit = flag;
#else
  if (implicit_bits) {
    implicit_bits->set(size_t(p), flag);
  } else {
    if (!flag) return;
    implicit_bits = new bitvector(*this);
    implicit_bits->expand(a_size);
    implicit_bits->set(size_t(p), flag);
  }
#endif
}
                         
// ******************************************************************                       

#ifdef OLD_NODE_HEADERS
inline MEDDLY::node_handle
MEDDLY::node_headers::getNextOf(node_handle p) const
#else
inline size_t
MEDDLY::node_headers::getNextOf(size_t p) const
#endif
{
  MEDDLY_DCASSERT(p>0);
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(p<=a_size);
  MEDDLY_DCASSERT(address);
  MEDDLY_DCASSERT(0==address[p].level);
  return address[p].offset;
#else
  MEDDLY_DCASSERT(addresses);
  MEDDLY_DCASSERT(levels);
  MEDDLY_DCASSERT(0==levels->get(size_t(p)));
  return addresses->get(size_t(p));
#endif
}

// ******************************************************************

inline void
#ifdef OLD_NODE_HEADERS
MEDDLY::node_headers::setNextOf(node_handle p, node_handle n)
#else
MEDDLY::node_headers::setNextOf(size_t p, size_t n)
#endif
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  MEDDLY_DCASSERT(n>=0);
#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  MEDDLY_DCASSERT(0==address[p].level);
  address[p].offset = node_address(n);
#else
  MEDDLY_DCASSERT(addresses);
  MEDDLY_DCASSERT(levels);
  MEDDLY_DCASSERT(0==levels->get(size_t(p)));
  addresses->set(size_t(p), node_address(n));
#endif
}

// ******************************************************************
// *                                                                *
// *              inlined memory_manager_style methods              *
// *                                                                *
// ******************************************************************

inline const char* MEDDLY::memory_manager_style::getName() const
{
  return name;
}

// ******************************************************************
// *                                                                *
// *                 inlined memory_manager methods                 *
// *                                                                *
// ******************************************************************

inline const char* MEDDLY::memory_manager::getStyleName() const
{
  return style_name;
}

inline void* MEDDLY::memory_manager::getChunkAddress(node_address h) const
{
  MEDDLY_DCASSERT(isValidHandle(h));

  return chunk_multiplier 
    ?  chunk_base + chunk_multiplier * h
    :  slowChunkAddress(h);
}

inline void MEDDLY::memory_manager::incMemUsed(size_t b)
{
  my_mem.incMemUsed(b);
}

inline void MEDDLY::memory_manager::decMemUsed(size_t b)
{
  my_mem.decMemUsed(b);
}

inline void MEDDLY::memory_manager::incMemAlloc(size_t b)
{
  my_mem.incMemAlloc(b);
}

inline void MEDDLY::memory_manager::decMemAlloc(size_t b)
{
  my_mem.decMemAlloc(b);
}

inline void MEDDLY::memory_manager::zeroMemUsed()
{
  my_mem.zeroMemUsed();
}

inline void MEDDLY::memory_manager::zeroMemAlloc()
{
  my_mem.zeroMemAlloc();
}

inline void MEDDLY::memory_manager::setChunkBase(void* p)
{
  chunk_base = (char*) p;
}

inline void MEDDLY::memory_manager::setChunkMultiplier(unsigned int m)
{
  chunk_multiplier = m;
}


// ******************************************************************
// *                                                                *
// *               inlined node_storage_style methods               *
// *                                                                *
// ******************************************************************

inline const char* MEDDLY::node_storage_style::getName() const
{
  return name;
}

// ******************************************************************
// *                                                                *
// *                  inlined node_storage methods                  *
// *                                                                *
// ******************************************************************


inline const char*
MEDDLY::node_storage::getStyleName() const
{
  return style_name;
}

inline const MEDDLY::expert_forest*
MEDDLY::node_storage::getParent() const
{
  MEDDLY_DCASSERT(parent);
  return parent;
}

inline MEDDLY::expert_forest*
MEDDLY::node_storage::getParent()
{
  MEDDLY_DCASSERT(parent);
  return parent;
}

inline void
MEDDLY::node_storage::moveNodeOffset(MEDDLY::node_handle node, node_address old_addr,
    node_address new_addr)
{
  MEDDLY_DCASSERT(parent);
  parent->moveNodeOffset(node, old_addr, new_addr);
}

// ******************************************************************
// *                                                                *
// *                 inlined  expert_forest methods                 *
// *                                                                *
// ******************************************************************

inline MEDDLY::node_handle
MEDDLY::expert_forest::bool_Tencoder::value2handle(bool v)
{
  return v ? -1 : 0;
}

inline bool
MEDDLY::expert_forest::bool_Tencoder::handle2value(MEDDLY::node_handle h)
{
  if (-1 == h)
    return true;
  if (0 == h)
    return false;
  throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::int_Tencoder::value2handle(int v)
{
  MEDDLY_DCASSERT(4 == sizeof(MEDDLY::node_handle));
  if (v < -1073741824 || v > 1073741823) {
    // Can't fit in 31 bits (signed)
    throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
  }
  if (v)
    v |= 0x80000000; // sets the sign bit
  return v;
}

inline int
MEDDLY::expert_forest::int_Tencoder::handle2value(MEDDLY::node_handle h)
{
  // << 1 kills the sign bit
  // >> 1 puts us back, and extends the (new) sign bit
  return (h << 1) >> 1;
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::float_Tencoder::value2handle(float v)
{
  MEDDLY_DCASSERT(4 == sizeof(MEDDLY::node_handle));
  MEDDLY_DCASSERT(sizeof(float) <= sizeof(MEDDLY::node_handle));
  if (0.0 == v)
    return 0;
  intfloat x;
  x.real = v;
  // strip lsb in fraction, and add sign bit
  return node_handle( (unsigned(x.integer) >> 1) | 0x80000000 );
}

inline float
MEDDLY::expert_forest::float_Tencoder::handle2value(MEDDLY::node_handle h)
{
  MEDDLY_DCASSERT(4 == sizeof(MEDDLY::node_handle));
  MEDDLY_DCASSERT(sizeof(float) <= sizeof(MEDDLY::node_handle));
  if (0 == h)
    return 0.0;
  intfloat x;
  x.integer = (h << 1); // remove sign bit
  return x.real;
}

template<typename T>
inline size_t
MEDDLY::expert_forest::EVencoder<T>::edgeBytes()
{
  return sizeof(T);
}
template<typename T>
inline void
MEDDLY::expert_forest::EVencoder<T>::writeValue(void* ptr, T val)
{
  memcpy(ptr, &val, sizeof(T));
}
template<typename T>
inline void
MEDDLY::expert_forest::EVencoder<T>::readValue(const void* ptr, T &val)
{
  memcpy(&val, ptr, sizeof(T));
}

template<typename T>
inline void MEDDLY::expert_forest::EVencoder<T>::show(output &s, const void* ptr)
{
  T val;
  readValue(ptr, val);
  s << val;
}

namespace MEDDLY {

template<>
inline void expert_forest::EVencoder<int>::write(output &s, const void* ptr)
{
  int val;
  readValue(ptr, val);
  s << val;
}

template<>
inline void expert_forest::EVencoder<long>::write(output &s, const void* ptr)
{
  long val;
  readValue(ptr, val);
  s << val;
}

template<>
inline void expert_forest::EVencoder<float>::write(output &s, const void* ptr)
{
  float val;
  readValue(ptr, val);
  s.put(val, 8, 8, 'e');
}

template<>
inline void expert_forest::EVencoder<double>::write(output &s, const void* ptr)
{
  double val;
  readValue(ptr, val);
  s.put(val, 8, 8, 'e');
}

template<>
inline void expert_forest::EVencoder<int>::read(input &s, void* ptr)
{
  writeValue(ptr, int(s.get_integer()));
}

template<>
inline void expert_forest::EVencoder<long>::read(input &s, void* ptr)
{
  writeValue(ptr, long(s.get_integer()));
}

template<>
inline void expert_forest::EVencoder<float>::read(input &s, void* ptr)
{
  writeValue(ptr, float(s.get_real()));
}

template<>
inline void expert_forest::EVencoder<double>::read(input &s, void* ptr)
{
  writeValue(ptr, double(s.get_real()));
}

}

template<typename T>
  inline MEDDLY::node_handle
  MEDDLY::expert_forest::handleForValue(T v) const
  {
    switch (getRangeType()) {
      case BOOLEAN:
        return bool_Tencoder::value2handle(v);
      case INTEGER:
        return int_Tencoder::value2handle(v);
      case REAL:
        return float_Tencoder::value2handle(v);
      default:
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }
  }

template<typename T>
  inline void
  MEDDLY::expert_forest::getValueFromHandle(MEDDLY::node_handle n, T& v) const
  {
    MEDDLY_DCASSERT(isTerminalNode(n));
    switch (getRangeType()) {
      case BOOLEAN:
        v = bool_Tencoder::handle2value(n);
        return;
      case INTEGER:
        v = int_Tencoder::handle2value(n);
        return;
      case REAL:
        v = float_Tencoder::handle2value(n);
        return;
      default:
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }
  }

inline bool
MEDDLY::expert_forest::getBooleanFromHandle(MEDDLY::node_handle n) const
{
  MEDDLY_DCASSERT(isTerminalNode(n));
  switch (getRangeType()) {
    case BOOLEAN:
      return bool_Tencoder::handle2value(n);
    case INTEGER:
      return int_Tencoder::handle2value(n);
    case REAL:
      return float_Tencoder::handle2value(n);
    default:
      throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  }
}

inline int
MEDDLY::expert_forest::getIntegerFromHandle(MEDDLY::node_handle n) const
{
  MEDDLY_DCASSERT(isTerminalNode(n));
  switch (getRangeType()) {
    case BOOLEAN:
      return bool_Tencoder::handle2value(n);
    case INTEGER:
      return int_Tencoder::handle2value(n);
    case REAL:
      return float_Tencoder::handle2value(n);
    default:
      throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  }
}

inline float
MEDDLY::expert_forest::getRealFromHandle(MEDDLY::node_handle n) const
{
  MEDDLY_DCASSERT(isTerminalNode(n));
  switch (getRangeType()) {
    case BOOLEAN:
      return bool_Tencoder::handle2value(n);
    case INTEGER:
      return int_Tencoder::handle2value(n);
    case REAL:
      return float_Tencoder::handle2value(n);
    default:
      throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  }
}

inline MEDDLY::memstats&
MEDDLY::expert_forest::changeMemStats()
{
  return mstats;
}

inline unsigned char
MEDDLY::expert_forest::edgeBytes() const
{
  return edge_bytes;
}

inline bool
MEDDLY::expert_forest::areEdgeValuesHashed() const
{
  return hash_edge_values;
}

inline unsigned char
MEDDLY::expert_forest::unhashedHeaderBytes() const
{
  return unhashed_bytes;
}

inline unsigned char
MEDDLY::expert_forest::hashedHeaderBytes() const
{
  return hashed_bytes;
}

inline const MEDDLY::expert_domain*
MEDDLY::expert_forest::getExpertDomain() const
{
  return (MEDDLY::expert_domain*) getDomain();
}

inline MEDDLY::expert_domain*
MEDDLY::expert_forest::useExpertDomain()
{
  return (MEDDLY::expert_domain*) useDomain();
}

// --------------------------------------------------
// Node address information
// --------------------------------------------------


inline MEDDLY::node_address
MEDDLY::expert_forest::getNodeAddress(node_handle p) const
{
  return nodeHeaders.getNodeAddress(p);
}

inline void
MEDDLY::expert_forest::setNodeAddress(node_handle p, node_address a)
{
  nodeHeaders.setNodeAddress(p, a);
}

// --------------------------------------------------
// Node level information
// --------------------------------------------------

inline int
MEDDLY::expert_forest::getNodeLevel(node_handle p) const
{
  if (isTerminalNode(p)) return 0;
  return nodeHeaders.getNodeLevel(p);
}

inline bool
MEDDLY::expert_forest::isPrimedNode(node_handle p) const
{
  return getNodeLevel(p) < 0;
}

inline bool
MEDDLY::expert_forest::isUnprimedNode(node_handle p) const
{
  return getNodeLevel(p) > 0;
}

inline int
MEDDLY::expert_forest::getNumVariables() const
{
  return getDomain()->getNumVariables();
}

inline int
MEDDLY::expert_forest::getMinLevelIndex() const
{
  return isForRelations() ? -getNumVariables() : 0;
}

inline bool
MEDDLY::expert_forest::isValidLevel(int k) const
{
  return (k >= getMinLevelIndex()) && (k <= getNumVariables());
}

inline bool
MEDDLY::expert_forest::isExtensibleLevel(int k) const
{
  MEDDLY_DCASSERT(isValidLevel(k));
  return getDomain()->getVar(k < 0? -k: k)->isExtensible();
}

inline int
MEDDLY::expert_forest::getLevelSize(int lh) const
{
  MEDDLY_DCASSERT(isValidLevel(lh));
  int var=getVarByLevel(lh);
  if (var < 0) {
    return getDomain()->getVariableBound(-var, true);
  }
  else {
    return getDomain()->getVariableBound(var, false);
  }
}

inline int
MEDDLY::expert_forest::getVariableSize(int var) const
{
  return getDomain()->getVariableBound(var, false);
}

inline void
MEDDLY::expert_forest::setNodeLevel(node_handle p, int level)
{
  nodeHeaders.setNodeLevel(p, level);
}

// --------------------------------------------------
// Managing incoming edge counts
// --------------------------------------------------

/*
inline bool
MEDDLY::expert_forest::trackingInCounts() const
{
  return nodeHeaders.trackingIncomingCounts();
}
*/

inline unsigned long
MEDDLY::expert_forest::getNodeInCount(MEDDLY::node_handle p) const
{
  MEDDLY_DCASSERT(deflt.useReferenceCounts);
  return nodeHeaders.getIncomingCount(p);
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::linkNode(MEDDLY::node_handle p)
{
  if (deflt.useReferenceCounts) {
    return nodeHeaders.linkNode(p);
  } else {
    return p;
  }
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::linkNode(const MEDDLY::dd_edge &p)
{
  MEDDLY_DCASSERT(p.getForest() == this);
  if (deflt.useReferenceCounts) {
    return nodeHeaders.linkNode(p.getNode());
  } else {
    return p.getNode();
  }
}

inline void
MEDDLY::expert_forest::unlinkNode(MEDDLY::node_handle p)
{
  if (deflt.useReferenceCounts) {
    nodeHeaders.unlinkNode(p);
  } 
}

inline void
MEDDLY::expert_forest::markNode(node_handle p)
{
  if (deflt.useReferenceCounts) return;
  if (p<1) return;
  if (nodeHeaders.hasReachableBit(p)) return;

#ifdef DEBUG_MARK_SWEEP
  printf("Marking node %d\n", p);
#endif

  nodeHeaders.setReachableBit(p);

  MEDDLY_DCASSERT(nodeMan);
  nodeMan->markDownPointers( nodeHeaders.getNodeAddress(p) );
}

inline bool
MEDDLY::expert_forest::hasReachableBit(node_handle p) const
{
  MEDDLY_DCASSERT(!deflt.useReferenceCounts);
  return nodeHeaders.hasReachableBit(p);
}

// --------------------------------------------------
// Managing cache counts
// --------------------------------------------------

/*
inline bool
MEDDLY::expert_forest::trackingCacheCounts() const
{
  return nodeHeaders.trackingCacheCounts();
}
*/

inline void
MEDDLY::expert_forest::cacheNode(MEDDLY::node_handle p)
{
  if (deflt.useReferenceCounts) {
    nodeHeaders.cacheNode(p);
    return;
  } 
}

inline void
MEDDLY::expert_forest::uncacheNode(MEDDLY::node_handle p)
{
  if (deflt.useReferenceCounts) {
    nodeHeaders.uncacheNode(p);
  }
}

inline void
MEDDLY::expert_forest::setCacheBit(MEDDLY::node_handle p)
{
  if (!deflt.useReferenceCounts) {
    nodeHeaders.setInCacheBit(p);
  }
}

inline void
MEDDLY::expert_forest::clearAllCacheBits()
{
  if (!deflt.useReferenceCounts) {
#ifdef DEBUG_MARK_SWEEP
    printf("Clearing cache bits for forest %u\n", FID());
#endif
    nodeHeaders.clearAllInCacheBits();
  }
}

inline void
MEDDLY::expert_forest::sweepAllCacheBits()
{
  if (!deflt.useReferenceCounts) {
#ifdef DEBUG_MARK_SWEEP
    printf("Sweeping cache bits for forest %u\n", FID());
#endif
    nodeHeaders.sweepAllInCacheBits();
  }
}

// --------------------------------------------------
// Node status
// --------------------------------------------------

inline bool
MEDDLY::expert_forest::isActiveNode(node_handle p) const
{
  return nodeHeaders.isActive(p);
}

/*
inline bool
MEDDLY::expert_forest::isZombieNode(node_handle p) const
{
  return nodeHeaders.isZombie(p);
}
*/

inline bool
MEDDLY::expert_forest::isDeletedNode(node_handle p) const
{
  return nodeHeaders.isDeleted(p);
}

inline bool
MEDDLY::expert_forest::isTerminalNode(MEDDLY::node_handle p)
{
  return (p < 1);
}


inline bool
MEDDLY::expert_forest::isValidNonterminalIndex(MEDDLY::node_handle node) const
{
  const node_handle a_last = nodeHeaders.lastUsedHandle();
  return (node > 0) && (node <= a_last);
}

inline bool
MEDDLY::expert_forest::isValidNodeIndex(MEDDLY::node_handle node) const
{
  const node_handle a_last = nodeHeaders.lastUsedHandle();
  return node <= a_last;
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::getLastNode() const
{
  const node_handle a_last = nodeHeaders.lastUsedHandle();
  return a_last;
}

// --------------------------------------------------
// Extensible Node Information:
// --------------------------------------------------
inline bool
MEDDLY::expert_forest::isExtensible(node_handle p) const
{
  return nodeMan->isExtensible(getNodeAddress(p));
}

//
// Unorganized from here
//

inline int
MEDDLY::expert_forest::getIndexSetCardinality(MEDDLY::node_handle node) const
{
  MEDDLY_DCASSERT(isIndexSet());
  if (isTerminalNode(node)) return (node != 0) ? 1 : 0;
  // yes iff the unhashed extra header is non-zero.
  const int* uhh = (const int*) nodeMan->getUnhashedHeaderOf(getNodeAddress(node));
  MEDDLY_DCASSERT(*uhh > 0);
  return *uhh;
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::getNext(MEDDLY::node_handle p) const
{
  return nodeMan->getNextOf(getNodeAddress(p));
}

inline void
MEDDLY::expert_forest::setNext(MEDDLY::node_handle p, MEDDLY::node_handle n)
{
  nodeMan->setNextOf(getNodeAddress(p), n);
}
inline bool 
MEDDLY::expert_forest::isImplicit(node_handle p) const
{
  return nodeHeaders.getNodeImplicitFlag(p);
}

inline unsigned
MEDDLY::expert_forest::hash(MEDDLY::node_handle p) const
{
  return hashNode(p);
}

inline MEDDLY::forest::node_status
MEDDLY::expert_forest::getNodeStatus(MEDDLY::node_handle node) const
{
  if (isMarkedForDeletion()) {
    return MEDDLY::forest::DEAD;
  }
  if (isTerminalNode(node)) {
    return terminalNodesStatus;
  }
  if (isDeletedNode(node)) {
    return MEDDLY::forest::DEAD;
  }
  // Active node.

  // If we're using reference counts,
  // and the incoming count is zero,
  // then we must be using optimistic
  // and the node is stale but recoverable.

  // If we're NOT using reference counts,
  // since we're not a deleted node,
  // assume we are still active.

  if (deflt.useReferenceCounts) {
    if (getNodeInCount(node) == 0) {
      return MEDDLY::forest::RECOVERABLE;
    } 
  }

  return MEDDLY::forest::ACTIVE;
}

inline unsigned
MEDDLY::expert_forest::hashNode(MEDDLY::node_handle p) const
{
  return nodeMan->hashNode(getNodeLevel(p), getNodeAddress(p));
}

inline int
MEDDLY::expert_forest::getSingletonIndex(MEDDLY::node_handle p, MEDDLY::node_handle &down) const
{
  return nodeMan->getSingletonIndex(getNodeAddress(p), down);
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::getSingletonDown(MEDDLY::node_handle node, int index) const
{
  MEDDLY::node_handle down;
  if (getSingletonIndex(node, down) == index) return down;
  return 0;
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::getDownPtr(MEDDLY::node_handle p, int index) const
{
  return nodeMan->getDownPtr(getNodeAddress(p), index);
}

inline void
MEDDLY::expert_forest::getDownPtr(MEDDLY::node_handle p, int index, int& ev,
    MEDDLY::node_handle& dn) const
{
  nodeMan->getDownPtr(getNodeAddress(p), index, ev, dn);
}

inline void
MEDDLY::expert_forest::getDownPtr(MEDDLY::node_handle p, int index, long& ev,
    MEDDLY::node_handle& dn) const
{
  nodeMan->getDownPtr(getNodeAddress(p), index, ev, dn);
}

inline void
MEDDLY::expert_forest::getDownPtr(MEDDLY::node_handle p, int index, float& ev,
    MEDDLY::node_handle& dn) const
{
  nodeMan->getDownPtr(getNodeAddress(p), index, ev, dn);
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::getTransparentNode() const
{
  return transparent;
}

inline void 
MEDDLY::expert_forest::fillUnpacked(MEDDLY::unpacked_node &un, MEDDLY::node_handle node, unpacked_node::storage_style st2) 
const
{
  const int level = getNodeLevel(node);
  MEDDLY_DCASSERT(0 != level);
  un.bind_to_forest(this, level, unsigned(getLevelSize(level)), true); 
  MEDDLY_DCASSERT(getNodeAddress(node));
  nodeMan->fillUnpacked(un, getNodeAddress(node), st2);
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::createReducedNode(int in, MEDDLY::unpacked_node *un)
{
  MEDDLY_DCASSERT(un);
  MEDDLY_DCASSERT(un->isBuildNode());
  un->computeHash();
  MEDDLY::node_handle q = createReducedHelper(in, *un);
#ifdef TRACK_DELETIONS
  printf("Created node %d\n", q);
#endif
  unpacked_node::recycle(un);
  return q;
}

inline unsigned
MEDDLY::expert_forest::getImplicitTableCount() const
{
  unsigned q = getImplTableCount();
  return q;
}

inline MEDDLY::relation_node* 
MEDDLY::expert_forest::buildImplicitNode(node_handle rnh)
{
  return buildImplNode(rnh);
}

inline MEDDLY::node_handle 
MEDDLY::expert_forest::getImplicitTerminalNode() const
{
  return getImplTerminalNode();
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::createRelationNode(MEDDLY::relation_node *un)
{
  MEDDLY_DCASSERT(un);
  MEDDLY::node_handle q = createImplicitNode(*un);
#ifdef TRACK_DELETIONS
  printf("Created node %d\n", q);
#endif
  return q;
}


template<class T>
inline void
MEDDLY::expert_forest::createReducedNode(int in, MEDDLY::unpacked_node *un, T& ev,
      MEDDLY::node_handle& node)
{
  MEDDLY_DCASSERT(un);
  MEDDLY_DCASSERT(un->isBuildNode());
  normalize(*un, ev);
  MEDDLY_DCASSERT(ev >= 0);
  un->computeHash();
  node = createReducedHelper(in, *un);
#ifdef TRACK_DELETIONS
  printf("Created node %d\n", node);
#endif
  unpacked_node::recycle(un);
}

inline bool
MEDDLY::expert_forest::areDuplicates(MEDDLY::node_handle node, const MEDDLY::unpacked_node &nr) const
{
  MEDDLY_DCASSERT(node > 0);
  if (nodeHeaders.getNodeLevel(node) != nr.getLevel()) {
    return false;
  }
  return nodeMan->areDuplicates(nodeHeaders.getNodeAddress(node), nr);
}

inline void
MEDDLY::expert_forest::setEdgeSize(unsigned char ebytes, bool hashed)
{
  MEDDLY_DCASSERT(0 == edge_bytes);
  edge_bytes = ebytes;
  hash_edge_values = hashed;
}

inline void
MEDDLY::expert_forest::setUnhashedSize(unsigned char ubytes)
{
  MEDDLY_DCASSERT(0 == unhashed_bytes);
  unhashed_bytes = ubytes;
}

inline void
MEDDLY::expert_forest::setHashedSize(unsigned char hbytes)
{
  MEDDLY_DCASSERT(0 == hashed_bytes);
  hashed_bytes = hbytes;
}

/*
inline bool
MEDDLY::expert_forest::isTimeToGc() const
{
  return isPessimistic() ? (stats.zombie_nodes > deflt.zombieTrigger)
      : (stats.orphan_nodes > deflt.orphanTrigger);
}
*/

inline void
MEDDLY::expert_forest::moveNodeOffset(MEDDLY::node_handle node, node_address old_addr,
    node_address new_addr)
{
  nodeHeaders.moveNodeAddress(node, old_addr, new_addr);
}

inline void
MEDDLY::expert_forest::getVariableOrder(int* level2var) const
{
  // Assume sufficient space has been allocated for order
  level2var[0] = 0;
  for (int i = 1; i < getNumVariables() + 1; i++) {
    level2var[i] = var_order->getVarByLevel(i);
  }
}

inline std::shared_ptr<const MEDDLY::variable_order>
MEDDLY::expert_forest::variableOrder() const
{
  return var_order;
}

// ******************************************************************
// *                                                                *
// *                     inlined opname methods                     *
// *                                                                *
// ******************************************************************

inline int
MEDDLY::opname::getIndex() const
{
  return index;
}

inline const char*
MEDDLY::opname::getName() const
{
  return name;
}

// ******************************************************************
// *                                                                *
// *               inlined specialized_opname methods               *
// *                                                                *
// ******************************************************************

inline void
MEDDLY::specialized_opname::arguments::setAutoDestroy(bool destroy)
{
  destroyWhenDone = destroy;
}

inline bool
MEDDLY::specialized_opname::arguments::autoDestroy() const
{
  return destroyWhenDone;
}

// ******************************************************************
// *                                                                *
// *                inlined numerical_opname methods                *
// *                                                                *
// ******************************************************************

inline MEDDLY::specialized_operation*
MEDDLY::numerical_opname::buildOperation(const dd_edge &x_ind,
    const dd_edge &A, const dd_edge &y_ind) const
{
  numerical_args na(x_ind, A, y_ind);
  na.setAutoDestroy(false); // na will be destroyed when we return
  return buildOperation(&na);
}


// ******************************************************************
// *                                                                *
// *                inlined satpregen_opname methods                *
// *                                                                *
// ******************************************************************


inline bool
MEDDLY::satpregen_opname::pregen_relation::isFinalized() const
{
  return 0 == next;
}

inline MEDDLY::dd_edge*
MEDDLY::satpregen_opname::pregen_relation::arrayForLevel(int k) const
{
  MEDDLY_DCASSERT(isFinalized());
  MEDDLY_CHECK_RANGE(1, k, K + 1);
  if (level_index) {
    // "by events"
    if (level_index[k - 1] > level_index[k]) {
      return events + level_index[k];
    }
    else {
      // empty list
      return 0;
    }
  }
  else {
    // "by levels"
    return events+k;
  }
}

inline unsigned
MEDDLY::satpregen_opname::pregen_relation::lengthForLevel(int k) const
{
  MEDDLY_DCASSERT(isFinalized());
  MEDDLY_CHECK_RANGE(1, k, K + 1);
  if (level_index) {
    // "by events"
    return level_index[k - 1] - level_index[k];
  }
  else {
    // "by levels"
    return events[k].getNode() ? 1 : 0;
  }
}

inline MEDDLY::forest*
MEDDLY::satpregen_opname::pregen_relation::getInForest() const
{
  return insetF;
}

inline MEDDLY::forest*
MEDDLY::satpregen_opname::pregen_relation::getRelForest() const
{
  return mxdF;
}

inline MEDDLY::forest*
MEDDLY::satpregen_opname::pregen_relation::getOutForest() const
{
  return outsetF;
}


// ******************************************************************
// *                                                                *
// *                 inlined  satotf_opname methods                 *
// *                                                                *
// ******************************************************************


inline MEDDLY::expert_forest*
MEDDLY::satotf_opname::subevent::getForest() {
  return f;
} 

inline int
MEDDLY::satotf_opname::subevent::getNumVars() const {
  return num_vars;
}

inline const int*
MEDDLY::satotf_opname::subevent::getVars() const {
  return vars;
}

inline int
MEDDLY::satotf_opname::subevent::getTop() const {
  return top;
}

inline bool
MEDDLY::satotf_opname::subevent::isFiring() const {
  return is_firing;
}

inline bool
MEDDLY::satotf_opname::subevent::isEnabling() const {
  return !is_firing;
}

inline const MEDDLY::dd_edge&
MEDDLY::satotf_opname::subevent::getRoot() const {
  return root;
}

inline bool
MEDDLY::satotf_opname::subevent::usesExtensibleVariables() const {
  return uses_extensible_variables;
}

// ****************************************************************************

inline MEDDLY::expert_forest*
MEDDLY::satotf_opname::otf_relation::getInForest() const
{
  return insetF;
}

inline MEDDLY::expert_forest*
MEDDLY::satotf_opname::otf_relation::getRelForest() const
{
  return mxdF;
}

inline MEDDLY::expert_forest*
MEDDLY::satotf_opname::otf_relation::getOutForest() const
{
  return outsetF;
}

inline bool
MEDDLY::satotf_opname::otf_relation::isConfirmed(int level, int i) const
{
  if (level < num_levels &&  i >= 0) {
    return (i < insetF->getLevelSize(level) && confirmed[level][i]);
  }
  throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
}

inline int
MEDDLY::satotf_opname::otf_relation::getNumOfEvents(int level) const
{
  MEDDLY_CHECK_RANGE(1, level, num_levels);
  return num_events_by_top_level[level];
}

inline const MEDDLY::dd_edge&
MEDDLY::satotf_opname::otf_relation::getEvent(int level, int i) const
{
  MEDDLY_CHECK_RANGE(0, i, getNumOfEvents(level));
  return events_by_top_level[level][i]->getRoot();
}


inline bool
MEDDLY::satotf_opname::otf_relation::rebuildEvent(int level, int i)
{
  MEDDLY_CHECK_RANGE(0, i, getNumOfEvents(level));
  return events_by_top_level[level][i]->rebuild();
}

inline const bool*
MEDDLY::satotf_opname::otf_relation::getLocalStates(int level)
{
  MEDDLY_CHECK_RANGE(0, level, num_levels);
  return confirmed[level];
}

inline int
MEDDLY::satotf_opname::otf_relation::getNumConfirmed(int level) const
{
  MEDDLY_CHECK_RANGE(0, level, num_levels);
  return num_confirmed[level];
}

// ******************************************************************
// *                                                                *
// *                 inlined  relation_node methods                *
// *                                                                *
// ******************************************************************


inline unsigned long
MEDDLY::relation_node::getSignature() const
{
  return signature;
}

inline int
MEDDLY::relation_node::getLevel() const
{
  return level;
}

inline node_handle
MEDDLY::relation_node::getDown() const
{
  return down;
}

inline void
MEDDLY::relation_node::setDown(node_handle d) 
{
  down = d;
}

inline node_handle
MEDDLY::relation_node::getID() const
{
  return ID;
}

inline void
MEDDLY::relation_node::setID(node_handle n_ID)
{
  ID=n_ID;
}

inline long
MEDDLY::relation_node::getFire() const
{
  return fire;
}

inline void
MEDDLY::relation_node::setFire(long fire_val)
{
  fire = fire_val;
}

inline long
MEDDLY::relation_node::getEnable() const
{
  return enable;
}

inline void
MEDDLY::relation_node::setEnable(long enable_val)
{
  enable = enable_val;
}


inline long
MEDDLY::relation_node::getPieceSize() const
{
  return piece_size;
}

inline void
MEDDLY::relation_node::setPieceSize(long pS)
{
  piece_size=pS;
}

inline long*
MEDDLY::relation_node::getTokenUpdate() const
{
  return token_update;
}

inline
void
MEDDLY::relation_node::setTokenUpdate(long* n_token_update)
{
  token_update = n_token_update;
}

// ******************************************************************
// *                                                                *
// *                 inlined  satimpl_opname methods                *
// *                                                                *
// ******************************************************************


/*inline MEDDLY::relation_node*
MEDDLY::satimpl_opname::implicit_relation::nodeExists(node_handle n)
{
  std::unordered_map<node_handle, relation_node*>::iterator finder = impl_unique.find(n);
  if(finder!=impl_unique.end())
    return finder->second;
  else
    return NULL;
}

inline bool
MEDDLY::satimpl_opname::implicit_relation::isReserved(node_handle n)
{
  return (n==1);
}*/

//************************************************************************

inline MEDDLY::expert_forest*
MEDDLY::satimpl_opname::implicit_relation::getInForest() const
{
  return insetF;
}

inline MEDDLY::expert_forest*
MEDDLY::satimpl_opname::implicit_relation::getOutForest() const
{
  return outsetF;
}

inline MEDDLY::expert_forest*
MEDDLY::satimpl_opname::implicit_relation::getMixRelForest() const
{
  return mixRelF;
}

// ***********************************************************************

inline long
MEDDLY::satimpl_opname::implicit_relation::getTotalEvent(int level)
{
  int total_event = 0;
  for(int i=1;i<=level;i++)
    total_event +=  lengthForLevel(i);
  
  return total_event;
}

inline long
MEDDLY::satimpl_opname::implicit_relation::lengthForLevel(int level) const
{
  return event_added[level];
}

inline node_handle*
MEDDLY::satimpl_opname::implicit_relation::arrayForLevel(int level) const
{
  return event_list[level];
}

// ****************************************************************************

inline MEDDLY::expert_forest*
MEDDLY::satimpl_opname::implicit_relation::getRelForest() const
{
  return mixRelF;
}


inline long
MEDDLY::satimpl_opname::implicit_relation::getConfirmedStates(int level) const
{
  return confirm_states[level];
}

inline void
MEDDLY::satimpl_opname::implicit_relation::setConfirmedStates(int level,int i)
{
  resizeConfirmedArray(level,i);
  MEDDLY_DCASSERT(confirmed_array_size[level]>i);
  if(!isConfirmedState(level,i))
    {
      confirmed[level][i]=true;
      confirm_states[level]++;
    }
}

inline bool
MEDDLY::satimpl_opname::implicit_relation::isConfirmedState(int level,int i)
{

  return (i < insetF->getLevelSize(level) && confirmed[level][i]);
}

// ****************************************************************************


inline MEDDLY::expert_forest*
MEDDLY::satimpl_opname::subevent::getForest() {
  return f;
} 

inline int
MEDDLY::satimpl_opname::subevent::getNumVars() const {
  return num_vars;
}

inline const int*
MEDDLY::satimpl_opname::subevent::getVars() const {
  return vars;
}

inline int
MEDDLY::satimpl_opname::subevent::getTop() const {
  return top;
}

inline bool
MEDDLY::satimpl_opname::subevent::isFiring() const {
  return is_firing;
}

inline bool
MEDDLY::satimpl_opname::subevent::isEnabling() const {
  return !is_firing;
}

inline const MEDDLY::dd_edge&
MEDDLY::satimpl_opname::subevent::getRoot() const {
  return root;
}

inline bool
MEDDLY::satimpl_opname::subevent::usesExtensibleVariables() const {
  return uses_extensible_variables;
}

// ******************************************************************
// *                                                                *
// *                   inlined ct_object  methods                   *
// *                                                                *
// ******************************************************************


// ******************************************************************
// *                                                                *
// *                 inlined  compute_table methods                 *
// *                                                                *
// ******************************************************************

inline void
MEDDLY::compute_table::entry_key::setup(const compute_table::entry_type* et, unsigned repeats)
{
  MEDDLY_DCASSERT(et);
  etype = et;
  num_repeats = repeats;
  MEDDLY_DCASSERT( 0==repeats || et->isRepeating() );
  total_slots = et->getKeySize(repeats);
  if (total_slots > data_alloc) {
    data_alloc = (1+(data_alloc / 8)) * 8;   // allocate in chunks of size 8
    data = (entry_item*) realloc(data, data_alloc*sizeof(entry_item));
    if (0==data) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(data, 0, total_slots * sizeof(entry_item));
  currslot = 0;
#ifdef DEVELOPMENT_CODE
  has_hash = false;
#endif
}

inline const MEDDLY::compute_table::entry_type*
MEDDLY::compute_table::entry_key::getET() const
{
  return etype;
}

inline void MEDDLY::compute_table::entry_key::writeN(node_handle nh) 
{
  MEDDLY_CHECK_RANGE(0, currslot, total_slots);
  MEDDLY_DCASSERT(compute_table::NODE == theSlotType());
  data[currslot++].N = nh;
}

inline void MEDDLY::compute_table::entry_key::writeI(int i) 
{
  MEDDLY_CHECK_RANGE(0, currslot, total_slots);
  MEDDLY_DCASSERT(compute_table::INTEGER == theSlotType());
  data[currslot++].I = i;
}

inline void MEDDLY::compute_table::entry_key::writeL(long i) 
{
  MEDDLY_CHECK_RANGE(0, currslot, total_slots);
  MEDDLY_DCASSERT(compute_table::LONG == theSlotType());
  data[currslot++].L = i;
}

inline void MEDDLY::compute_table::entry_key::writeF(float f) 
{
  MEDDLY_CHECK_RANGE(0, currslot, total_slots);
  MEDDLY_DCASSERT(compute_table::FLOAT == theSlotType());
  data[currslot++].F = f;
}

inline const MEDDLY::compute_table::entry_item* 
MEDDLY::compute_table::entry_key::rawData() const 
{ 
  return data;
}

inline unsigned MEDDLY::compute_table::entry_key::dataLength() const
{ 
  return total_slots;
}

inline unsigned MEDDLY::compute_table::entry_key::numRepeats() const
{ 
  return num_repeats;
}

inline const void*
MEDDLY::compute_table::entry_key::readTempData() const
{
  return temp_data;
}

inline unsigned
MEDDLY::compute_table::entry_key::numTempBytes() const
{
  return temp_bytes;
}

inline void*
MEDDLY::compute_table::entry_key::allocTempData(unsigned bytes)
{
  temp_bytes = bytes;
  if (bytes > temp_alloc) {
    temp_alloc = (1+(temp_bytes/64)) * 64;    // allocate in chunks of 64 bytes
    temp_data = realloc(temp_data, temp_alloc);
    if (0==temp_data) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  return temp_data;
}

inline void
MEDDLY::compute_table::entry_key::cacheNodes() const
{
  for (unsigned i=0; i<total_slots; i++) {
    expert_forest* f = etype->getKeyForest(i);
    if (f) {
      f->cacheNode(data[i].N);
    }
  }
}

inline unsigned MEDDLY::compute_table::entry_key::getHash() const
{
  MEDDLY_DCASSERT(has_hash);
  return hash_value;
}

inline void MEDDLY::compute_table::entry_key::setHash(unsigned h)
{
  hash_value = h;
#ifdef DEVELOPMENT_CODE
  has_hash = true;
#endif
}

inline MEDDLY::compute_table::typeID MEDDLY::compute_table::entry_key::theSlotType() const
{
  //
  // Adjust currslot for OP entry, and number of repeats entry
  //
  // return etype->getKeyType(currslot - (etype->isRepeating() ? 2 : 1) );
  return etype->getKeyType(currslot);
}

// ******************************************************************

inline MEDDLY::node_handle MEDDLY::compute_table::entry_result::readN() 
{
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(compute_table::NODE == etype->getResultType(currslot));
  return data[currslot++].N;
}

inline int MEDDLY::compute_table::entry_result::readI()
{
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(compute_table::INTEGER == etype->getResultType(currslot));
  return data[currslot++].I;
}

inline float MEDDLY::compute_table::entry_result::readF()
{
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(compute_table::FLOAT == etype->getResultType(currslot));
  return data[currslot++].F;
}

inline long MEDDLY::compute_table::entry_result::readL()
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(compute_table::LONG == etype->getResultType(currslot));
  return data[currslot++].L;
}

inline double MEDDLY::compute_table::entry_result::readD()
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(compute_table::DOUBLE == etype->getResultType(currslot));
  return data[currslot++].D;
}

inline MEDDLY::ct_object* MEDDLY::compute_table::entry_result::readG()
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(compute_table::GENERIC == etype->getResultType(currslot));
  return data[currslot++].G;
}


inline void MEDDLY::compute_table::entry_result::reset() 
{
  currslot = 0;
}

inline void MEDDLY::compute_table::entry_result::writeN(node_handle nh)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(compute_table::NODE == etype->getResultType(currslot));
  build[currslot++].N = nh;
}

inline void MEDDLY::compute_table::entry_result::writeI(int i)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(compute_table::INTEGER == etype->getResultType(currslot));
  build[currslot++].I = i;
}

inline void MEDDLY::compute_table::entry_result::writeF(float f)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(compute_table::FLOAT == etype->getResultType(currslot));
  build[currslot++].F = f;
}

inline void MEDDLY::compute_table::entry_result::writeL(long L)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(compute_table::LONG == etype->getResultType(currslot));
  build[currslot++].L = L;
}

inline void MEDDLY::compute_table::entry_result::writeD(double D)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(compute_table::DOUBLE == etype->getResultType(currslot));
  build[currslot++].D = D;
}

inline void MEDDLY::compute_table::entry_result::writeG(ct_object* G)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(compute_table::GENERIC == etype->getResultType(currslot));
  build[currslot++].G = G;
}

inline void
MEDDLY::compute_table::entry_result::setValid()
{
  is_valid = true;
  data = build;
}

inline void
MEDDLY::compute_table::entry_result::setValid(const entry_item* d)
{
  is_valid = true;
  data = d;
}

inline void
MEDDLY::compute_table::entry_result::setInvalid()
{
  is_valid = false;
}

inline
MEDDLY::compute_table::entry_result::operator bool() const
{
  return is_valid;
}

inline void
MEDDLY::compute_table::entry_result::cacheNodes() const
{
  for (unsigned i=0; i<etype->getResultSize(); i++) {
    expert_forest* f = etype->getResultForest(i);
    if (f) {
      f->cacheNode(build[i].N);
    }
  }
}

inline const MEDDLY::compute_table::entry_item* 
MEDDLY::compute_table::entry_result
::rawData() const
{
  return build;
}

inline unsigned MEDDLY::compute_table::entry_result
::dataLength() const
{
  return etype->getResultSize();
}


// ******************************************************************

inline unsigned MEDDLY::compute_table::entry_type::getID() const
{
  return etID;
}

inline void MEDDLY::compute_table::entry_type::mightUpdateResults()
{
  updatable_result = true;
}

inline bool MEDDLY::compute_table::entry_type::isResultUpdatable() const
{
  return updatable_result;
}

inline const char* MEDDLY::compute_table::entry_type
::getName() const
{
  return name;
}

inline bool MEDDLY::compute_table::entry_type::isRepeating() const
{
  return len_kr_type;
}

inline unsigned MEDDLY::compute_table::entry_type
::getKeySize(unsigned reps) const
{
  return len_ks_type + (reps * len_kr_type);
}

inline unsigned MEDDLY::compute_table::entry_type
::getKeyBytes(unsigned reps) const
{
  return ks_bytes + (reps * kr_bytes);
}

inline void MEDDLY::compute_table::entry_type
::getKeyType(unsigned i, typeID &t, expert_forest* &f) const
{
  if (i<len_ks_type) {
    MEDDLY_DCASSERT(ks_type);
    MEDDLY_DCASSERT(ks_forest);
    t = ks_type[i];
    f = ks_forest[i];
    return;
  }
  MEDDLY_DCASSERT(len_kr_type);
  i -= len_ks_type;
  i %= len_kr_type;
  MEDDLY_DCASSERT(kr_type);
  MEDDLY_DCASSERT(kr_forest);
  t = kr_type[i];
  f = kr_forest[i];
}

inline MEDDLY::compute_table::typeID MEDDLY::compute_table::entry_type
::getKeyType(unsigned i) const
{
  if (i<len_ks_type) {
    MEDDLY_DCASSERT(ks_type);
    return ks_type[i];
  }
  MEDDLY_DCASSERT(len_kr_type);
  i -= len_ks_type;
  i %= len_kr_type;
  MEDDLY_DCASSERT(kr_type);
  return kr_type[i];
}

inline MEDDLY::expert_forest* MEDDLY::compute_table::entry_type
::getKeyForest(unsigned i) const
{
  MEDDLY_DCASSERT(ks_forest);
  if (i<len_ks_type) {
    return ks_forest[i];
  }
  MEDDLY_DCASSERT(len_kr_type);
  i -= len_ks_type;
  i %= len_kr_type;
  return kr_forest[i];
}

inline unsigned MEDDLY::compute_table::entry_type
::getResultSize() const
{
  return len_r_type;
}

inline unsigned MEDDLY::compute_table::entry_type
::getResultBytes() const
{
  return r_bytes;
}

inline void MEDDLY::compute_table::entry_type
::getResultType(unsigned i, typeID &t, expert_forest* &f) const
{
  MEDDLY_CHECK_RANGE(0, i, len_r_type);
  MEDDLY_DCASSERT(r_type);
  MEDDLY_DCASSERT(r_forest);
  t = r_type[i];
  f = r_forest[i];
}

inline MEDDLY::compute_table::typeID MEDDLY::compute_table::entry_type
::getResultType(unsigned i) const
{
  MEDDLY_CHECK_RANGE(0, i, len_r_type);
  MEDDLY_DCASSERT(r_type);
  return r_type[i];
}

inline MEDDLY::expert_forest* MEDDLY::compute_table::entry_type
::getResultForest(unsigned i) const
{
  MEDDLY_CHECK_RANGE(0, i, len_r_type);
  MEDDLY_DCASSERT(r_forest);
  return r_forest[i];
}

inline void MEDDLY::compute_table::entry_type
::markForDeletion()
{
  is_marked_for_deletion = true;
}

inline void MEDDLY::compute_table::entry_type
::unmarkForDeletion()
{
  is_marked_for_deletion = false;
}

inline bool MEDDLY::compute_table::entry_type
::isMarkedForDeletion() const
{
  return is_marked_for_deletion;
}

// ******************************************************************

inline bool
MEDDLY::compute_table::isOperationTable() const
{
  return global_et;
}

// convenience methods, for grabbing edge values
inline void
MEDDLY::compute_table::readEV(const MEDDLY::node_handle* p, int &ev)
{
  ev = p[0];
}
inline void
MEDDLY::compute_table::readEV(const MEDDLY::node_handle* p, long &ev)
{
  long* l = (long*) p;
  ev = l[0];
}
inline void
MEDDLY::compute_table::readEV(const MEDDLY::node_handle* p, float &ev)
{
  float* f = (float*) p;
  ev = f[0];
}

inline MEDDLY::compute_table::entry_key*
MEDDLY::compute_table::useEntryKey(const entry_type* et, unsigned repeats)
{
  if (0==et) return 0;
  MEDDLY_DCASSERT( (0==repeats) || et->isRepeating() );

  entry_key* k;
  if (free_keys) {
    k = free_keys;
    free_keys = free_keys->next;
  } else {
    k = new entry_key();
  }
  k->setup(et, repeats);
  return k;
}

inline void
MEDDLY::compute_table::recycle(entry_key* k)
{
  if (k) {
    k->next = free_keys;
    free_keys = k;
  }
}

inline const MEDDLY::compute_table::stats&
MEDDLY::compute_table::getStats()
{
  return perf;
}

inline const MEDDLY::compute_table::entry_type*
MEDDLY::compute_table::getEntryType(operation* op, unsigned slot)
{
  MEDDLY_DCASSERT(op);
  MEDDLY_CHECK_RANGE(0, slot, op->getNumETids());
  unsigned etid = op->getFirstETid() + slot;
  MEDDLY_CHECK_RANGE(0, etid, entryInfoSize);
  return entryInfo[etid];
}

inline const MEDDLY::compute_table::entry_type*
MEDDLY::compute_table::getEntryType(unsigned etid)
{
  MEDDLY_CHECK_RANGE(0, etid, entryInfoSize);
  return entryInfo[etid];
}

inline void
MEDDLY::compute_table::setHash(entry_key *k, unsigned h)
{
  MEDDLY_DCASSERT(k);
  k->setHash(h);
}

// ******************************************************************
// *                                                                *
// *                   inlined  operation methods                   *
// *                                                                *
// ******************************************************************


inline void
MEDDLY::operation::registerInForest(MEDDLY::forest* f)
{
  if (f)
    f->registerOperation(this);
}

inline void
MEDDLY::operation::unregisterInForest(MEDDLY::forest* f)
{
  if (f)
    f->unregisterOperation(this);
}

inline void
MEDDLY::operation::registerEntryType(unsigned slot, compute_table::entry_type* et)
{
  MEDDLY_CHECK_RANGE(0, slot, num_etids);
  MEDDLY_DCASSERT(etype);
  MEDDLY_DCASSERT(0==etype[slot]);
  etype[slot] = et;
  compute_table::registerEntryType(first_etid + slot, et);
}

inline bool
MEDDLY::operation::isMarkedForDeletion() const
{
  return is_marked_for_deletion;
}

inline void
MEDDLY::operation::setNext(operation* n)
{
  next = n;
}

inline MEDDLY::operation*
MEDDLY::operation::getNext()
{
  return next;
}

inline bool
MEDDLY::operation::usesMonolithicComputeTable()
{
  return Monolithic_CT;
}

inline unsigned
MEDDLY::operation::getIndex() const
{
  return oplist_index;
}

inline MEDDLY::operation*
MEDDLY::operation::getOpWithIndex(unsigned i)
{
  MEDDLY_CHECK_RANGE(0, i, list_size);
  return op_list[i];
}

inline unsigned
MEDDLY::operation::getOpListSize()
{
  return list_size;
}

inline void
MEDDLY::operation::setFirstETid(unsigned slot)
{
  first_etid = slot;
}

inline unsigned
MEDDLY::operation::getFirstETid() const
{
  return first_etid;
}

inline unsigned
MEDDLY::operation::getNumETids() const
{
  return num_etids;
}

inline const char*
MEDDLY::operation::getName() const
{
  return theOpName->getName();
}

inline const MEDDLY::opname*
MEDDLY::operation::getOpName() const
{
  return theOpName;
}

// ******************************************************************
// *                                                                *
// *                inlined  unary_operation methods                *
// *                                                                *
// ******************************************************************

inline bool
MEDDLY::unary_operation::matches(const MEDDLY::expert_forest* arg,
    const MEDDLY::expert_forest* res) const
{
  return (arg == argF && res == resF);
}

inline bool
MEDDLY::unary_operation::matches(const MEDDLY::expert_forest* arg, opnd_type res) const
{
  return (arg == argF && resultType == res);
}

inline bool
MEDDLY::unary_operation::checkForestCompatibility() const
{
  if (resultType == FOREST) {
    auto o1 = argF->variableOrder();
    auto o2 = resF->variableOrder();
    return o1->is_compatible_with(*o2);
  }
  else {
    return true;
  }
}

inline void
MEDDLY::unary_operation::compute(const dd_edge &arg, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(arg, res);
}

// ******************************************************************
// *                                                                *
// *                inlined binary_operation methods                *
// *                                                                *
// ******************************************************************


inline bool
MEDDLY::binary_operation::matches(const MEDDLY::expert_forest* arg1,
    const MEDDLY::expert_forest* arg2, const MEDDLY::expert_forest* res) const
{
  return (arg1 == arg1F && arg2 == arg2F && res == resF);
}

inline void
MEDDLY::binary_operation::operationCommutes()
{
  can_commute = (arg1F == arg2F);
}

inline bool
MEDDLY::binary_operation::checkForestCompatibility() const
{
  auto o1 = arg1F->variableOrder();
  auto o2 = arg2F->variableOrder();
  auto o3 = resF->variableOrder();
  return o1->is_compatible_with(*o2) && o1->is_compatible_with(*o3);
}

inline void
MEDDLY::binary_operation::compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(ar1, ar2, res);
}

// ******************************************************************
// *                                                                *
// *                       inlined  functions                       *
// *                                                                *
// ******************************************************************


inline MEDDLY::unary_operation*
MEDDLY::getOperation(const unary_opname* code, const dd_edge& arg,
    const dd_edge& res)
{
  return getOperation(code, (MEDDLY::expert_forest*) arg.getForest(),
      (MEDDLY::expert_forest*) res.getForest());
}

inline MEDDLY::unary_operation*
MEDDLY::getOperation(const unary_opname* code, const dd_edge& arg,
    opnd_type res)
{
  return getOperation(code, (MEDDLY::expert_forest*) arg.getForest(), res);
}

inline MEDDLY::binary_operation*
MEDDLY::getOperation(const binary_opname* code, const dd_edge& arg1,
    const dd_edge& arg2, const dd_edge& res)
{
  printf("\n HERER");
  return getOperation(code, (MEDDLY::expert_forest*) arg1.getForest(),
      (MEDDLY::expert_forest*) arg2.getForest(), (MEDDLY::expert_forest*) res.getForest());
}

#endif




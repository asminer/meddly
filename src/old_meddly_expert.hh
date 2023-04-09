
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

#include "defines.h"
#include "io.h"
#include "memstats.h"

#include "forest.h"

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
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1, k, K + 1);
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
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1, k, K + 1);
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
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1, level, num_levels);
  return num_events_by_top_level[level];
}

inline const MEDDLY::dd_edge&
MEDDLY::satotf_opname::otf_relation::getEvent(int level, int i) const
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, getNumOfEvents(level));
  return events_by_top_level[level][i]->getRoot();
}


inline bool
MEDDLY::satotf_opname::otf_relation::rebuildEvent(int level, int i)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, getNumOfEvents(level));
  return events_by_top_level[level][i]->rebuild();
}

inline const bool*
MEDDLY::satotf_opname::otf_relation::getLocalStates(int level)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, level, num_levels);
  return confirmed[level];
}

inline int
MEDDLY::satotf_opname::otf_relation::getNumConfirmed(int level) const
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, level, num_levels);
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

inline MEDDLY::expert_forest*
MEDDLY::relation_node::getForest() {
  return f;
}

inline int
MEDDLY::relation_node::getLevel() const
{
  return level;
}

inline rel_node_handle
MEDDLY::relation_node::getDown() const
{
  return down;
}

inline void
MEDDLY::relation_node::setDown(rel_node_handle d)
{
  down = d;
}


inline rel_node_handle
MEDDLY::relation_node::getID() const
{
  return ID;
}

inline void
MEDDLY::relation_node::setID(rel_node_handle n_ID)
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

inline void
MEDDLY::relation_node::setInhibit(long inh)
{
  inhibit = inh;
}

inline long
MEDDLY::relation_node::getInhibit() const
{
  return inhibit;
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


inline MEDDLY::relation_node*
MEDDLY::satimpl_opname::implicit_relation::nodeExists(rel_node_handle n)
{
  std::unordered_map<rel_node_handle, relation_node*>::iterator finder = impl_unique.find(n);
  if(finder!=impl_unique.end())
    return finder->second;
  else
    return NULL;
}

inline bool
MEDDLY::satimpl_opname::implicit_relation::isReserved(rel_node_handle n)
{
  return (n==1);
}

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

inline rel_node_handle*
MEDDLY::satimpl_opname::implicit_relation::arrayForLevel(int level) const
{
  return event_list[level];
}

// ****************************************************************************

inline MEDDLY::expert_forest*
MEDDLY::satimpl_opname::implicit_relation::getRelForest() const
{
  return mxdF;
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

// ******************************************************************
// *                                                                *
// *                 inlined  sathyb_opname methods                 *
// *                                                                *
// ******************************************************************

inline MEDDLY::expert_forest*
MEDDLY::sathyb_opname::hybrid_relation::getInForest() const
{
  return insetF;
}

inline MEDDLY::expert_forest*
MEDDLY::sathyb_opname::hybrid_relation::getOutForest() const
{
  return outsetF;
}

inline MEDDLY::expert_forest*
MEDDLY::sathyb_opname::hybrid_relation::getHybridForest() const
{
  return hybRelF;
}

// ***********************************************************************

inline long
MEDDLY::sathyb_opname::hybrid_relation::getTotalEvent(int level)
{
  int total_event = 0;
  for(int i=1;i<=level;i++)
    total_event +=  lengthForLevel(i);

  return total_event;
}

inline long
MEDDLY::sathyb_opname::hybrid_relation::lengthForLevel(int level) const
{
  return num_events_by_top_level[level];
}

inline MEDDLY::sathyb_opname::event**
MEDDLY::sathyb_opname::hybrid_relation::arrayForLevel(int level) const
{
  return events_by_top_level[level];
}

// ****************************************************************************

inline int
MEDDLY::sathyb_opname::hybrid_relation::getConfirmedStates(int level) const
{
  return confirm_states[level];
}


inline bool
MEDDLY::sathyb_opname::hybrid_relation::isConfirmedState(int level,int i)
{
  return (i < insetF->getLevelSize(level) && (i < confirmed_array_size[level]) && confirmed[level][i]);
}

// ******************************************************************


inline MEDDLY::expert_forest*
MEDDLY::sathyb_opname::subevent::getForest() {
  return f;
}

inline int
MEDDLY::sathyb_opname::subevent::getNumVars() const {
  return num_vars;
}

inline const int*
MEDDLY::sathyb_opname::subevent::getVars() const {
  return vars;
}

inline int
MEDDLY::sathyb_opname::subevent::getTop() const {
  return top;
}

inline bool
MEDDLY::sathyb_opname::subevent::isFiring() const {
  return is_firing;
}

inline bool
MEDDLY::sathyb_opname::subevent::isEnabling() const {
  return !is_firing;
}

inline bool
MEDDLY::sathyb_opname::subevent::isImplicit() const {
  return num_vars == 1;
}

inline const MEDDLY::dd_edge&
MEDDLY::sathyb_opname::subevent::getRoot() const {
  return root;
}

inline const MEDDLY::node_handle
MEDDLY::sathyb_opname::subevent::getRootHandle() const {
  return root_handle;
}

inline void
MEDDLY::sathyb_opname::subevent::setRootHandle( node_handle ID ) {
  root_handle = ID;
}

inline MEDDLY::node_handle
MEDDLY::sathyb_opname::subevent::getDown() const {
  return down;
}

inline void
MEDDLY::sathyb_opname::subevent::setDown( node_handle d_ID ) {
   down = d_ID;
}

inline int
MEDDLY::sathyb_opname::subevent::getEnable() const {
  return enable;
}

inline int
MEDDLY::sathyb_opname::subevent::getFire() const {
  return fire;
}

inline bool
MEDDLY::sathyb_opname::subevent::usesExtensibleVariables() const {
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
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(compute_table::NODE == theSlotType());
  data[currslot++].N = nh;
}

inline void MEDDLY::compute_table::entry_key::writeI(int i)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(compute_table::INTEGER == theSlotType());
  data[currslot++].I = i;
}

inline void MEDDLY::compute_table::entry_key::writeL(long i)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(compute_table::LONG == theSlotType());
  data[currslot++].L = i;
}

inline void MEDDLY::compute_table::entry_key::writeF(float f)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
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
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, len_r_type);
  MEDDLY_DCASSERT(r_type);
  MEDDLY_DCASSERT(r_forest);
  t = r_type[i];
  f = r_forest[i];
}

inline MEDDLY::compute_table::typeID MEDDLY::compute_table::entry_type
::getResultType(unsigned i) const
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, len_r_type);
  MEDDLY_DCASSERT(r_type);
  return r_type[i];
}

inline MEDDLY::expert_forest* MEDDLY::compute_table::entry_type
::getResultForest(unsigned i) const
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, len_r_type);
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
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, slot, op->getNumETids());
  unsigned etid = op->getFirstETid() + slot;
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, etid, entryInfoSize);
  return entryInfo[etid];
}

inline const MEDDLY::compute_table::entry_type*
MEDDLY::compute_table::getEntryType(unsigned etid)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, etid, entryInfoSize);
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
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, slot, num_etids);
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
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, list_size);
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
  if (resultType == opnd_type::FOREST) {
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
  computeDDEdge(arg, res, true);
}

inline void
MEDDLY::unary_operation::computeTemp(const dd_edge &arg, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(arg, res, false);
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
  computeDDEdge(ar1, ar2, res, true);
}

inline void
MEDDLY::binary_operation::computeTemp(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(ar1, ar2, res, false);
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
  return getOperation(code, (MEDDLY::expert_forest*) arg1.getForest(),
      (MEDDLY::expert_forest*) arg2.getForest(), (MEDDLY::expert_forest*) res.getForest());
}

#endif





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

#include "defines.h"
#include "hash_stream.h"

#define MERGE_AND_SPLIT_HOLES
// #define DEBUG_COMPACTION
// #define DEBUG_SLOW
// #define MEMORY_TRACE
// #define DEEP_MEMORY_TRACE



#ifndef TEST_ENCODING // HACK - turn off lots of code, for testing

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_reader  methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_reader::node_reader()
{
  extra_hashed = 0;
  ext_alloc = 0;
  ext_size = 0;
  down = 0;
  index = 0;
  edge = 0;
  alloc = 0;
  ealloc = 0;
  size = 0;
  nnzs = 0;
  level = 0;
#ifdef DEVELOPMENT_CODE
  has_hash = false;
#endif
}

MEDDLY::node_reader::~node_reader()
{
  clear();
}

void MEDDLY::node_reader::clear()
{
  free(extra_hashed);
  free(down);
  free(index);
  free(edge);
  down = 0;
  index = 0;
  edge = 0;
  alloc = 0;
  ealloc = 0;
  size = 0;
  nnzs = 0;
  level = 0;
}

void MEDDLY::node_reader
::resize(int k, int ns, char eb, bool full)
{
  level = k;
  is_full = full;
  size = ns;
  edge_bytes = eb;
  if (size > alloc) {
    int nalloc = ((ns/8)+1)*8;
    MEDDLY_DCASSERT(nalloc > ns);
    MEDDLY_DCASSERT(nalloc>0);
    MEDDLY_DCASSERT(nalloc>alloc);
    down = (long*) realloc(down, nalloc*sizeof(long));
    if (0==down) throw error(error::INSUFFICIENT_MEMORY);
    index = (int*) realloc(index, nalloc*sizeof(int));
    if (0==index) throw error(error::INSUFFICIENT_MEMORY);
    alloc = nalloc;
  }
  if (edge_bytes * size > ealloc) {
    int nalloc = ((edge_bytes * size)/8+1)*8;
    MEDDLY_DCASSERT(nalloc>0);
    MEDDLY_DCASSERT(nalloc>ealloc);
    edge = realloc(edge, nalloc);
    if (0==edge) throw error(error::INSUFFICIENT_MEMORY);
    ealloc = nalloc;
  }
}

void MEDDLY::node_reader
::resize_header(int extra_slots)
{
  ext_size = extra_slots;
  if (ext_size > ext_alloc) {
    ext_alloc = ((ext_size/8)+1)*8;
    MEDDLY_DCASSERT(ext_alloc > ext_size);
    MEDDLY_DCASSERT(ext_alloc>0);
    extra_hashed = (int*) realloc(extra_hashed, ext_alloc * sizeof(int));
    if (0==extra_hashed) throw error(error::INSUFFICIENT_MEMORY);
  }
}

void MEDDLY::node_reader::computeHash(bool hashEdgeValues)
{
#ifdef DEVELOPMENT_CODE
  MEDDLY_DCASSERT(!has_hash);
#endif
  
  hash_stream s;
  s.start(level);

  for (int e=0; e<ext_size; e++) {
    MEDDLY_DCASSERT(extra_hashed);
    s.push(extra_hashed[e]);
  }
  
  if (isSparse()) {
    if (hashEdgeValues) {
      for (int z=0; z<nnzs; z++) {
        MEDDLY_DCASSERT(d(z));
        s.push(i(z), d(z), ei(z));
      }
    } else {
      for (int z=0; z<nnzs; z++) {
        MEDDLY_DCASSERT(d(z));
        s.push(i(z), d(z));
      }
    }
  } else {
    if (hashEdgeValues) {
      for (int n=0; n<size; n++) {
        if (0==d(n)) continue;
        s.push(n, d(n), ei(n));
      }
    } else {
      for (int n=0; n<size; n++) {
        if (0==d(n)) continue;
        s.push(n, d(n));
      }
    }
  }
  h = s.finish();
#ifdef DEVELOPMENT_CODE
  has_hash = true;
#endif
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_builder methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_builder::node_builder()
{
  parent = 0;
  extra_hashed = 0;
  extra_unhashed = 0;
  down = 0;
  indexes = 0;
  edge = 0;
}

MEDDLY::node_builder::~node_builder()
{
  free(extra_hashed);
  free(extra_unhashed);
  free(down);
  free(indexes);
  free(edge);
}

void MEDDLY::node_builder::init(int k, const expert_forest* p)
{
  MEDDLY_DCASSERT(p);
  parent = p;
  level = k;
  hhsize = parent->hashedHeaderSize();
  if (hhsize) {
    extra_hashed = (int*) malloc(uhsize * sizeof(int));
    if (0==extra_hashed)
      throw error(error::INSUFFICIENT_MEMORY);
  }
  uhsize = parent->unhashedHeaderSize();
  if (uhsize) {
    extra_unhashed = (int*) malloc(hhsize * sizeof(int));
    if (0==extra_unhashed)
      throw error(error::INSUFFICIENT_MEMORY);
  }
  size = 0;
  alloc = 0;
  lock = false;
  hashEdgeValues = parent->areEdgeValuesHashed();
  edge_bytes = parent->edgeSize() * sizeof(int);
}


void MEDDLY::node_builder::computeHash()
{
#ifdef DEVELOPMENT_CODE
  MEDDLY_DCASSERT(!has_hash);
#endif
  
  hash_stream s;
  s.start(level);

  for (int e=0; e<hhsize; e++) {
    MEDDLY_DCASSERT(extra_hashed);
    s.push(extra_hashed[e]);
  }
 /* 
  if (is_sparse) {
    if (hashEdgeValues) {
      for (int z=0; z<size; z++) {
        MEDDLY_DCASSERT(down[z]);
        int* ep = (int*) ( (char*)edge + z * edge_bytes );
        s.push(indexes[z], down[z], *ep);
      }
    } else {
      for (int z=0; z<size; z++) {
        MEDDLY_DCASSERT(down[z]);
        s.push(indexes[z], down[z]);
      }
    }
  } else {
    if (hashEdgeValues) {
      for (int i=0; i<size; i++) {
        if (0==down[i]) continue;
        int* ep = (int*) ( (char*)edge + i * edge_bytes );
        s.push(i, down[i], *ep);
      }
    } else {
      for (int i=0; i<size; i++) {
        if (0==down[i]) continue;
        s.push(i, down[i]);
      }
    }
  }
  */
  if (isSparse()) {
    if (hashEdgeValues) {
      for (int z=0; z<size; z++) {
        MEDDLY_DCASSERT(d(z));
        s.push(i(z), d(z), ei(z));
      }
    } else {
      for (int z=0; z<size; z++) {
        MEDDLY_DCASSERT(d(z));
        s.push(i(z), d(z));
      }
    }
  } else {
    if (hashEdgeValues) {
      for (int n=0; n<size; n++) {
        if (0==d(n)) continue;
        s.push(n, d(n), ei(n));
      }
    } else {
      for (int n=0; n<size; n++) {
        if (0==d(n)) continue;
        s.push(n, d(n));
      }
    }
  }
  h = s.finish();
#ifdef DEVELOPMENT_CODE
  has_hash = true;
#endif
}


void MEDDLY::node_builder::enlarge()
{
  if (size <= alloc) return;
  alloc = ((size / 8)+1) * 8;
  MEDDLY_DCASSERT(alloc > size);
  down = (long*) realloc(down, alloc * sizeof(long));
  if (0==down) throw error(error::INSUFFICIENT_MEMORY);
  indexes = (int*) realloc(indexes, alloc * sizeof(int));
  if (0==indexes) throw error(error::INSUFFICIENT_MEMORY);
  int es = parent->edgeSize();
  if (es>0) {
    edge = realloc(edge, alloc * es * sizeof(int));
    if (0==edge) throw error(error::INSUFFICIENT_MEMORY);
  }
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_storage methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************


/*

  Hole management.
  ==============================================================
  There are two kinds of holes depending on their location in the "hole grid":
  Index Holes and Non-index Holes.

  The hole grid structure:
  ------------------------
  (holes_bottom)
  holes_of_size_0 (index) -- holes_of_size_0 (non_index) -- (non_index) -- NULL
  |
  holes_of_size_1 -- ""
  |
  :
  :
  (holes_top)

  TBD:
  Note that the grid only stores holes up to the largest hole
  requested.  Larger holes are stored in the "large holes list".

  Index holes are represented as follows:
  ---------------------------------------
  [0] -size (number of slots in hole)     
  [1] up
  [2] down 
  [3] next pointer (nodes of same size)
  [4..size-2] Unused
  :
  :
  [size-1] -size

  Non-index holes are represented as follows:
  [0] -size (number of slots in hole)     
  [1] flag (<0, indicates non-index node)
  [2] prev pointer (nodes of same size)
  [3] next pointer (nodes of same size)
  [4..size-2] Unused
  :
  :
  [size-1] -size
*/

MEDDLY::node_storage::node_storage()
{
  parent = 0;
}

MEDDLY::node_storage::~node_storage()
{
  if (parent) parent->changeStats().decMemAlloc(size*sizeof(int));
  free(data);
}

void MEDDLY::node_storage::init(expert_forest* p)
{
  MEDDLY_DCASSERT(0==parent);
  parent = p;
  data = 0;
  data_down = 0;
  size = 0;
  last = 0;
  max_request = 0;
  large_holes = 0;
  holes_top = 0;
  holes_bottom = 0;
  hole_slots = 0;
  edgeSize = parent->edgeSize();
  unhashedHeader = parent->unhashedHeaderSize();
  hashedHeader = parent->hashedHeaderSize();
}

void MEDDLY::node_storage::unlinkDown(int addr)
{
  scanner ns;
  initScanner(ns, addr);
  for (; ns; ++ns) {
    int dp = ns.d();
    if (dp) parent->unlinkNode(dp);
  }
}

void MEDDLY::node_storage::compact(bool shrink)
{
  //
  // Should we even bother?
  //
  if (0==data || 0==hole_slots) return;

  if (hole_slots <= parent->getPolicies().compact_min)  return;

  if (hole_slots <  parent->getPolicies().compact_max) {

    // If percentage of holes is below trigger, then don't compact

    if (100 * hole_slots < last * parent->getPolicies().compact_frac) return;

  }

#ifdef DEBUG_SLOW
  fprintf(stderr, "Compacting forest level\n");
#endif
#ifdef MEMORY_TRACE
  printf("Compacting\n");
#endif
#ifdef DEBUG_COMPACTION
  printf("Before compaction:\n");
  dumpInternal(stdout);
  printf("\n");
#endif

  //
  // Scan the whole array of data, copying over itself and skipping holes.
  // 
  int *node_ptr = data + 1;  // since we leave [0] empty
  int *end_ptr = data + last + 1;
  int *curr_ptr = node_ptr;

  while (node_ptr < end_ptr) {
    if (*node_ptr < 0) {
      //
      // This is a hole, skip it
      // 
      MEDDLY_DCASSERT(*node_ptr == node_ptr[-(*node_ptr)-1]);
      node_ptr += -(*node_ptr);
      continue;
    } 
    //
    // A real node, move it
    //
    MEDDLY_DCASSERT(!parent->isPessimistic() || *node_ptr != 0);
    
    int old_off = node_ptr - data;
    int new_off = curr_ptr - data;

    // copy the node, except for the padding.
    int datalen = slotsForNode(sizeOf(node_ptr - data)) - 1;
    if (node_ptr != curr_ptr) {
      memmove(curr_ptr, node_ptr, datalen * sizeof(int));
    }
    node_ptr += datalen;
    curr_ptr += datalen;
    //
    // Skip any padding
    //
    if (*node_ptr < 0) {
      node_ptr -= *node_ptr;  
    }
    //
    // Copy trailer, the node number
    //
    *curr_ptr = *node_ptr;
    parent->moveNodeOffset(*curr_ptr, old_off, new_off);
    curr_ptr++;
    node_ptr++;

  } // while
  MEDDLY_DCASSERT(node_ptr == end_ptr);

  last = (curr_ptr - 1 - data);

  // set up hole pointers and such
  holes_top = holes_bottom = 0;
  hole_slots = 0;
  large_holes = 0;

  parent->changeStats().num_compactions++;

  if (shrink && size > min_size && last < size/2) {
    int new_size = size/2;
    while (new_size > min_size && new_size > last * 3) { new_size /= 2; }
    resize(new_size);
  }

#ifdef DEBUG_COMPACTION
  printf("After compaction:\n");
  dumpInternal(stdout);
  printf("\n");
#endif
}

void MEDDLY::node_storage::dumpInternal(FILE *s) const
{
  if (0==data) return; // nothing to display
  
#ifdef NODE_STORAGE_PER_LEVEL
  fprintf(s, "Level %ld: ", parent->getLevelNumber(this));
#endif
  fprintf(s, "Last slot used: %d\n", last);
  fprintf(s, "large_holes: %d\n", large_holes);
  fprintf(s, "Grid: top = %d bottom = %d\n", holes_top, holes_bottom);

  fprintf(s, "Data array by record: \n");
  int awidth = digits(parent->getLastNode());
  int a;
  for (a=1; a<=last; ) {
    fflush(s);
    a = dumpNode(s, a);
  } // for a
  fprintf(s, "%*d : free slots\n", awidth, a);
  fflush(s);
  MEDDLY_DCASSERT(a == (last)+1);
}

void MEDDLY::node_storage::dumpInternal(FILE *s, int a) const
{
  MEDDLY_DCASSERT(data);
  
#ifdef NODE_STORAGE_PER_LEVEL
  fprintf(s, "Level %ld: ", parent->getLevelNumber(this));
#endif
  fprintf(s, "Last slot used: %d\n", last);
  fprintf(s, "large_holes: %d\n", large_holes);
  fprintf(s, "Grid: top = %d bottom = %d\n", holes_top, holes_bottom);

  dumpNode(s, a);
}

void MEDDLY::node_storage
::addToChainCounts(std::map<int, int> &chainLengths) const
{
  for (int curr = holes_bottom; curr; curr = holeUp(curr)) {
    int currHoleOffset = curr;
    int count = 0;
    // traverse this chain to get its length
    for (count = 0; currHoleOffset; count++) {
      currHoleOffset = holeNext(currHoleOffset);
    }
    int currHoleSize = -data[curr];
    chainLengths[currHoleSize] += count;
  }
  // add the large hole list
  int count = 0;
  for (int curr = large_holes; curr; curr = holeNext(curr)) {
    count++;
  }
  chainLengths[-1] += count;
}

void MEDDLY::node_storage::fillReader(int addr, node_reader &nr) const
{
  scanner ns;
  initScanner(ns, addr);

  // Copy hashed header

  if (hashedHeader) {
    nr.resize_header(hashedHeader);
    memcpy(nr.extra_hashed, HH(addr), hashedHeader * sizeof(int));
  }

  // Copy everything else

  if (ns.isFull()) {
    if (nr.isFull()) {
      for (; ns; ++ns) {
        nr.down[ns.count()] = ns.d();
      }
      for (int i=ns.fullSize(); i<nr.size; ++i) {
        nr.down[i] = 0;
      }
      if (nr.edge_bytes) {
        ns.restart();
        memcpy(nr.edge, ns.eptr(), ns.sizeOf() * nr.edge_bytes);
        void* evext = (char*)nr.edge + (ns.sizeOf() * nr.edge_bytes);
        memset(evext, 0, (nr.size-ns.sizeOf()) * nr.edge_bytes);
      }
      return;
    } 
    int& z = nr.nnzs;
    z = 0;
    if (nr.edge_bytes) {
      void* nev = nr.edge;
      for (; ns; ++ns) {
        if (0==ns.d()) continue;
        nr.down[z] = ns.d();
        nr.index[z] = ns.count();
        memcpy(nev, ns.eptr(), nr.edge_bytes);
        nev = (char*)nev + nr.edge_bytes;
        z++;
      } // for i
    } else {
      for (; ns; ++ns) if (ns.d()) {
        nr.down[z] = ns.d();
        nr.index[z] = ns.count();
        z++;
      } // for i
    } // if ev
    return; 
  }
  // 
  // ns must be sparse

  if (nr.isFull()) {
    if (nr.edge_bytes) {
      memset(nr.down, 0, nr.size * sizeof(int));
      memset(nr.edge, 0, nr.size * nr.edge_bytes);
      for (; ns; ++ns) {
        nr.down[ns.i()] = ns.d();
        int off = ns.i() * nr.edge_bytes;
        memcpy((char*)nr.edge + off, ns.eptr(), nr.edge_bytes);
      }
    } else {
      int i=0; 
      for (; ns; ++ns) {
        for (; i<ns.i(); i++) nr.down[i] = 0;
        nr.down[i] = ns.d();
        i++;
      }
      for (; i<nr.size; i++) nr.down[i] = 0;
    } // if ev
    return;
  }

  nr.nnzs = ns.sparseSize();
  for (; ns; ++ns) {
    nr.down[ns.count()] = ns.d();
    nr.index[ns.count()] = ns.i();
  }
  if (nr.edge_bytes) {
    ns.restart();
    memcpy(nr.edge, ns.eptr(), ns.sizeOf() * nr.edge_bytes);
  }
}

int MEDDLY::node_storage
::makeNode(int p, const node_builder &nb, forest::policies::node_storage opt)
{
  int nnz, truncsize;
  switch (opt) {
    case forest::policies::FULL_STORAGE:
          if (nb.isSparse()) {
            truncsize = -1;
            for (int z=0; z<nb.getNNZs(); z++) {
              truncsize = MAX(truncsize, nb.i(z));
            }
            if (truncsize<0) return 0;
            return makeFullNode(p, truncsize+1, nb);
          } else {
            for (int i=nb.getSize()-1; i>=0; i--) {
              if (nb.d(i)) return makeFullNode(p, i+1, nb);
            }
          }
          return 0;

    case forest::policies::SPARSE_STORAGE:
          if (nb.isSparse()) {
            return makeSparseNode(p, nb.getNNZs(), nb);
          } else {
            nnz = 0;
            for (int i=0; i<nb.getSize(); i++) {
              if (nb.d(i)) nnz++;
            }
            return makeSparseNode(p, nnz, nb);
          }

    case forest::policies::FULL_OR_SPARSE_STORAGE:
          truncsize = -1;
          if (nb.isSparse()) {
            nnz = nb.getNNZs();
            for (int z=0; z<nnz; z++) {
              truncsize = MAX(truncsize, nb.i(z));
            }
            if (truncsize<0) return 0;
            truncsize++;
          } else {
            nnz = 0;
            for (int i=0; i<nb.getSize(); i++) if (nb.d(i)) {
              nnz++;
              truncsize = i;
            }
            truncsize++;
          }
          if (slotsForNode(-nnz) < slotsForNode(truncsize)) { 
            return makeSparseNode(p, nnz, nb);
          } else {
            return makeFullNode(p, truncsize, nb);
          }

    default:
        throw error(error::MISCELLANEOUS);
  }
}

unsigned MEDDLY::node_storage::hashNode(const node_header& node) const
{
  hash_stream s;
  s.start(node.level);

  for (int e=0; e<hashedHeader; e++) {
    s.push(hashedHeaderOf(node.offset)[e]);
  }

  scanner ns;
  initScanner(ns, node.offset);

  if (ns.isSparse()) {
    if (parent->areEdgeValuesHashed()) {
      for (; ns; ++ns) {
        int ev;
        memcpy(&ev, ns.eptr(), sizeof(int));
        s.push(ns.i(), ns.d(), ev);
      }
    } else {
      for (; ns; ++ns) {
        MEDDLY_DCASSERT(ns.d());
        s.push(ns.i(), ns.d());
      }
    }
  } else {
    if (parent->areEdgeValuesHashed()) {
      for (; ns; ++ns) {
        if (0==ns.d()) continue;
        int ev;
        memcpy(&ev, ns.eptr(), sizeof(int));
        s.push(ns.count(), ns.d(), ev);
      }
    } else {
      for (; ns; ++ns) {
        if (0==ns.d()) continue;
        s.push(ns.count(), ns.d());
      }
    }
  }
  return s.finish();
}

int MEDDLY::node_storage::getSingletonIndex(int addr, long &down) const
{
  scanner ns;
  initScanner(ns, addr);

  if (ns.isFull()) {
    // full node
    for (; ns; ++ns) {
      if (0==ns.d()) continue;
      if (ns.count()+1 != ns.fullSize()) return -1;
      down = ns.d();
      return ns.count();
    }
    return -1;
  }
  // sparse node --- easy
  down = ns.d();
  int i = ns.i();
  ++ns;
  if (ns) return -1;
  return i;
}

long MEDDLY::node_storage
::getDownPtr(int addr, int index) const
{
  scanner ns;
  initScanner(ns, addr);
  return ns.moveToIndex(index) ? ns.d() : 0;
}

void MEDDLY::node_storage
::getDownPtr(int addr, int index, int& ev, long& dn) const
{
  scanner ns;
  initScanner(ns, addr);
  if (ns.moveToIndex(index)) {
    dn = ns.d();
    ev = ns.ei();
  } else {
    ev = 0;
    dn = 0;
  }
}

void MEDDLY::node_storage
::getDownPtr(int addr, int index, float& ev, long& dn) const
{
  scanner ns;
  initScanner(ns, addr);
  if (ns.moveToIndex(index)) {
    dn = ns.d();
    ev = ns.ef();
  } else {
    ev = 0;
    dn = 0;
  }
}

void MEDDLY::node_storage
::showNode(FILE* s, int addr, bool verbose) const
{
  scanner R;
  initScanner(R, addr);
  if (R.isFull()) {
    // Full node
    if (verbose) fprintf(s, " size: %d", R.fullSize());
    fprintf(s, " down: [");
    for (; R; ++R) {
      if (R.count()) fprintf(s, "|"); 
      if (R.hasEdges()) {
        fprintf(s, "<");
        parent->showEdgeValue(s, R.eptr(), 0);
        fprintf(s, ", ");
      } 
      if (parent->isTerminalNode(R.d())) {
        parent->showTerminal(s, R.d());
      } else {
        fprintf(s, "%d", R.d());
      }
      if (R.hasEdges()) fprintf(s, ">");
    } // for i
    fprintf(s, "]");
  } else {
    // Sparse node
    if (verbose) fprintf(s, " nnz : %d", R.sparseSize());
    fprintf(s, " down: (");
    for (; R; ++R) {
      if (R.count()) fprintf(s, ", ");
      fprintf(s, "%d:", R.i());
      if (R.hasEdges()) {
        fprintf(s, "<");
        parent->showEdgeValue(s, R.eptr(), 0);
        fprintf(s, ", ");
      } 
      if (parent->isTerminalNode(R.d())) {
        parent->showTerminal(s, R.d());
      } else {
        fprintf(s, "%d", R.d());
      }
      if (R.hasEdges()) fprintf(s, ">");
    } // for z
    fprintf(s, ")");
  }

  // show extra header stuff
  if (unhashedHeader) {
    parent->showUnhashedHeader(s, unhashedHeaderOf(addr));
  }
  if (hashedHeader) {
    parent->showHashedHeader(s, hashedHeaderOf(addr));
  }
}

//
//
// Private
//
//

int MEDDLY::node_storage::dumpNode(FILE *s, int a) const
{
  MEDDLY_DCASSERT(data);
  if (0==a) return 0;
  int awidth = digits(parent->getLastNode());
  fprintf(s, "%*d : [%d", awidth, a, data[a]);
  for (int i=1; i<3; i++) {
    fprintf(s, "|%d", data[a+i]);
  }
  if (data[a]<0) { 
    // hole
    fprintf(s, "| ... ");
    a -= data[a];  
  } else {
    // proper node
    int reqd = slotsForNode(sizeOf(a));
    int nElements = activeNodeActualSlots(a);
    for (int i=3; i<reqd-1; i++) {
      fprintf(s, "|%d", data[a+i]);
    }
    if (reqd < nElements) {
      fprintf(s, "| ... %d unused ... ", -data[a+reqd-1]);
    }
    a += activeNodeActualSlots(a);
  }
  fprintf(s, "|%d]\n", data[a-1]);
  return a;
}


int MEDDLY::node_storage
::makeFullNode(int p, int size, const node_builder &nb)
{
  int addr = allocNode(size, p, false);
  MEDDLY_DCASSERT(1==getCountOf(addr));
  int* down = FD(addr);
  if (edgeSize) {
      MEDDLY_DCASSERT(nb.hasEdges());
      char* edge = (char*) FE(addr);
      int edge_bytes = edgeSize * sizeof(int);
      if (nb.isSparse()) {
        memset(down, 0, size * sizeof(int));
        memset(edge, 0, size * edge_bytes);
        for (int z=0; z<nb.getNNZs(); z++) {
          int i = nb.i(z);
          MEDDLY_CHECK_RANGE(0, i, size);
          down[i] = nb.d(z);
          memcpy(edge + i * edge_bytes, nb.eptr(z), edge_bytes);
        }
      } else {
        for (int i=0; i<size; i++) down[i] = nb.d(i);
        // kinda hacky
        memcpy(edge, nb.eptr(0), size * edge_bytes);
      }
  } else {
      MEDDLY_DCASSERT(!nb.hasEdges());
      if (nb.isSparse()) {
        memset(down, 0, size * sizeof(int));
        for (int z=0; z<nb.getNNZs(); z++) {
          MEDDLY_CHECK_RANGE(0, nb.i(z), size);
          down[nb.i(z)] = nb.d(z);
        }
      } else {
        for (int i=0; i<size; i++) down[i] = nb.d(i);
      }
  }
  copyExtraHeader(addr, nb);
  return addr;
}


int MEDDLY::node_storage
::makeSparseNode(int p, int size, const node_builder &nb)
{
  int addr = allocNode(-size, p, false);
  MEDDLY_DCASSERT(1==getCountOf(addr));
  int* index = SI(addr);
  int* down  = SD(addr);
  if (nb.hasEdges()) {
      MEDDLY_DCASSERT(nb.hasEdges());
      char* edge = (char*) SE(addr);
      int edge_bytes = edgeSize * sizeof(int);
      if (nb.isSparse()) {
        for (int z=0; z<size; z++) {
          down[z] = nb.d(z);
          index[z] = nb.i(z);
        }
        // kinda hacky
        memcpy(edge, nb.eptr(0), size * edge_bytes);
      } else {
        int z = 0;
        for (int i=0; i<nb.getSize(); i++) if (nb.d(i)) {
          MEDDLY_CHECK_RANGE(0, z, size);
          down[z] = nb.d(i);
          index[z] = i;
          memcpy(edge + z * edge_bytes, nb.eptr(i), edge_bytes);
          z++;
        }
      }
  } else {
      MEDDLY_DCASSERT(!nb.hasEdges());
      if (nb.isSparse()) {
        for (int z=0; z<size; z++) {
          down[z] = nb.d(z);
          index[z] = nb.i(z);
        }
      } else {
        int z = 0;
        for (int i=0; i<nb.getSize(); i++) if (nb.d(i)) {
          MEDDLY_CHECK_RANGE(0, z, size);
          down[z] = nb.d(i);
          index[z] = i;
          z++;
        }
      }
  }
  copyExtraHeader(addr, nb);
  return addr;
}


void MEDDLY::node_storage::copyExtraHeader(int addr, const node_builder &nb)
{
  // copy extra header info, if any
  if (unhashedHeader) {
    int* uh = UH(addr);
    for (int i=0; i<unhashedHeader; i++) {
      uh[i] = nb.uh(i);
    }
  }
  if (hashedHeader) {
    int* hh = HH(addr);
    for (int i=0; i<hashedHeader; i++) {
      hh[i] = nb.hh(i);
    }
  }
}


//
//  NOTE: we set the first slot to contain
//  the NEGATIVE of the number of slots in the hole.
//  For recycled holes, usually that's a no-op :^)
//
int MEDDLY::node_storage::getHole(int slots)
{
  MEDDLY_DCASSERT(parent);

  if (slots > max_request) {
    max_request = slots;
    // Traverse the large hole list, and re-insert
    // all the holes, in case some are now smaller
    // than max_request.
    int curr = large_holes;
    large_holes = 0;
    for (; curr; ) {
      int next = holeNext(curr);
      gridInsert(curr);
      curr = next;
    }
  }

  // First, try for a hole exactly of this size
  // by traversing the index nodes in the hole grid
  int chain = 0;
  int curr = holes_bottom;
  while (curr) {
    if (slots == -(data[curr])) break;
    if (slots < -(data[curr])) {
      // no exact match possible
      curr = 0;
      break;
    }
    // move up the hole grid
    curr = holeUp(curr);
    chain++;
  }

  if (curr) {
    // perfect fit
    hole_slots -= slots;
    // try to not remove the "index" node
    int next = holeNext(curr);
    if (next) {
      midRemove(next);
#ifdef MEMORY_TRACE
      printf("Removed Non-Index hole %d\n", next);
#ifdef DEEP_MEMORY_TRACE
      dumpInternal(stdout);
#else 
      dumpInternal(stdout, next);
#endif
#endif
      return next;
    }
    indexRemove(curr);
#ifdef MEMORY_TRACE
    printf("Removed Index hole %d\n", curr);
#ifdef DEEP_MEMORY_TRACE
    dumpInternal(stdout);
#else 
    dumpInternal(stdout, curr);
#endif
#endif
    return curr;
  }

#ifdef MERGE_AND_SPLIT_HOLES
  // No hole with exact size.
  // Take the first large hole.
  if (large_holes) {
    curr = large_holes;
    large_holes = holeNext(curr);
    if (large_holes) holePrev(large_holes) = 0;

    // Sanity check: the hole is large enough
    MEDDLY_DCASSERT(slots < -data[curr]);
    
    // Is there space for a leftover hole?
    const int min_node_size = slotsForNode(0);
    if (slots + min_node_size >= -data[curr]) {
      // leftover space is too small to be useful,
      // just send the whole thing.
#ifdef MEMORY_TRACE
      printf("Removed entire large hole %d\n", curr);
#ifdef DEEP_MEMORY_TRACE
      dumpInternal(stdout);
#else 
      dumpInternal(stdout, curr);
#endif
#endif
      hole_slots += data[curr]; // SUBTRACTS the number of slots
      return curr;
    }

    // This is a large hole, save the leftovers
    int newhole = curr + slots;
    int newsize = -(data[curr]) - slots;
    data[newhole] = -newsize;
    data[newhole + newsize - 1] = -newsize;
    gridInsert(newhole);

#ifdef MEMORY_TRACE
    printf("Removed part of hole %d\n", curr);
#ifdef DEEP_MEMORY_TRACE
    dumpInternal(stdout);
#else 
    dumpInternal(stdout, curr);
#endif
#endif
    data[curr] = -slots;
    hole_slots -= slots;
    return curr;
  }
#endif

  // 
  // Still here?  We couldn't recycle a node.
  // 

  //
  // First -- try to compact if we need to expand
  //
  if (parent->getPolicies().compactBeforeExpand) {
    if (last + slots >= size) compact(false);
  }

  //
  // Do we need to expand?
  //
  if (last + slots >= size) {
    // new size is 50% more than previous (37.5% * 4 = 1.5 => 50% growth)
    int new_size = MAX( size, last + slots ) * 1.5;

    resize(new_size);
  }

  // 
  // Grab node from the end
  //
  int h = last + 1;
  last += slots;
  data[h] = -slots;
  return h;
}

void MEDDLY::node_storage::makeHole(int addr, int slots)
{
#ifdef MEMORY_TRACE
  printf("Calling makeHole(%d, %d)\n", addr, slots);
#endif
  MEDDLY_DCASSERT(parent);

  parent->changeStats().decMemUsed(slots * sizeof(int));

  hole_slots += slots;
  data[addr] = data[addr+slots-1] = -slots;

  if (!parent->getPolicies().recycleNodeStorageHoles) return;

  // Check for a hole to the left
#ifdef MERGE_AND_SPLIT_HOLES
  if (data[addr-1] < 0) {
    // Merge!
#ifdef MEMORY_TRACE
    printf("Left merging\n");
#endif
    int lefthole = addr + data[addr-1];
    MEDDLY_DCASSERT(data[lefthole] == data[addr-1]);
    if (non_index_hole == holeUp(lefthole)) midRemove(lefthole);
    else indexRemove(lefthole);
    slots += (-data[lefthole]);
    addr = lefthole;
    data[addr] = data[addr+slots-1] = -slots;
  }
#endif

  // if addr is the last hole, absorb into free part of array
  MEDDLY_DCASSERT(addr + slots - 1 <= last);
  if (addr+slots-1 == last) {
    last -= slots;
    hole_slots -= slots;
    if (size > min_size && (last + 1) < size/2) {
      int new_size = size/2;
      while (new_size > (last + 1) * 2) new_size /= 2;
      if (new_size < min_size) new_size = min_size;
      resize(new_size);
    }
#ifdef MEMORY_TRACE
    printf("Made Last Hole %d, last %d\n", addr, last);
#ifdef DEEP_MEMORY_TRACE
    dumpInternal(stdout);
#else 
    dumpInternal(stdout, 0);
#endif
#endif
    return;
  }

#ifdef MERGE_AND_SPLIT_HOLES
  // Check for a hole to the right
  if (data[addr+slots]<0) {
    // Merge!
#ifdef MEMORY_TRACE
    printf("Right merging\n");
#endif
    int righthole = addr+slots;
    if (non_index_hole == holeUp(righthole)) midRemove(righthole);
    else indexRemove(righthole);
    slots += (-data[righthole]);
    data[addr] = data[addr+slots-1] = -slots;
  }
#endif

  // Add hole to grid
  gridInsert(addr); 

#ifdef MEMORY_TRACE
  printf("Made Hole %d\n", addr);
#ifdef DEEP_MEMORY_TRACE
  dumpInternal(stdout);
#else
  dumpInternal(stdout, addr);
#endif
#endif
}

void MEDDLY::node_storage::gridInsert(int p_offset)
{
#ifdef MEMORY_TRACE
  printf("gridInsert(%d)\n", p_offset);
#endif

  // sanity check to make sure that the first and last slots in this hole
  // have the same value, i.e. -(# of slots in the hole)
  MEDDLY_DCASSERT(data[p_offset] == data[p_offset - data[p_offset] - 1]);

  // Check if we belong in the grid, or the large hole list
  if (-data[p_offset] > max_request) {
#ifdef MEMORY_TRACE
    printf("\tAdding to large_holes: %d\n", large_holes);
#endif
    // add to the large hole list
    holeUp(p_offset) = non_index_hole;
    holeNext(p_offset) = large_holes;
    holePrev(p_offset) = 0;
    if (large_holes) holePrev(large_holes) = p_offset;
    large_holes = p_offset;
    return;
  }

  // special case: empty
  if (0 == holes_bottom) {
#ifdef MEMORY_TRACE
    printf("\tAdding to empty grid\n");
#endif
    // index hole
    holeUp(p_offset) = 0;
    holeDown(p_offset) = 0;
    holeNext(p_offset) = 0;
    holes_top = holes_bottom = p_offset;
    return;
  }
  // special case: at top
  if (data[p_offset] < data[holes_top]) {
#ifdef MEMORY_TRACE
    printf("\tAdding new chain at top\n");
#endif
    // index hole
    holeUp(p_offset) = 0;
    holeNext(p_offset) = 0;
    holeDown(p_offset) = holes_top;
    holeUp(holes_top) = p_offset;
    holes_top = p_offset;
    return;
  }
  int above = holes_bottom;
  int below = 0;
  while (data[p_offset] < data[above]) {
    below = above;
    above = holeUp(below);
    MEDDLY_DCASSERT(holeDown(above) == below);
    MEDDLY_DCASSERT(above);  
  }
  if (data[p_offset] == data[above]) {
#ifdef MEMORY_TRACE
    printf("\tAdding to chain\n");
#endif
    // Found, add this to chain
    // making a non-index hole
    int right = holeNext(above);
    holeUp(p_offset) = non_index_hole;
    holePrev(p_offset) = above;
    holeNext(p_offset) = right;
    if (right) holePrev(right) = p_offset;
    holeNext(above) = p_offset;
    return; 
  }
#ifdef MEMORY_TRACE
  printf("\tAdding new chain\n");
#endif
  // we should have above < p_offset < below  (remember, -sizes)
  // create an index hole since there were no holes of this size
  holeUp(p_offset) = above;
  holeDown(p_offset) = below;
  holeNext(p_offset) = 0;
  holeDown(above) = p_offset;
  if (below) {
    holeUp(below) = p_offset;
  } else {
    MEDDLY_DCASSERT(above == holes_bottom);
    holes_bottom = p_offset;
  }
}

void MEDDLY::node_storage::midRemove(int p_offset)
{
#ifdef MEMORY_TRACE
  printf("midRemove(%d)\n", p_offset);
#endif

  // Remove a "middle" node, either in the grid
  // or in the large hole list.
  //
  MEDDLY_DCASSERT(isHoleNonIndex(p_offset));

  int left = holePrev(p_offset); 
  int right = holeNext(p_offset);

#ifdef MEMORY_TRACE
  printf("\tIN left: %d  right: %d  large_holes: %d\n", 
    left, right, large_holes
  );
#endif

  if (left) {
    holeNext(left) = right;
  } else {
    // MUST be head of the large hole list
    MEDDLY_DCASSERT(large_holes == p_offset);
    large_holes = right;
  }
  if (right) holePrev(right) = left;

#ifdef MEMORY_TRACE
  printf("\tOUT large_holes: %d\n", large_holes);
#endif

  // Sanity checks
#ifdef DEVELOPMENT_CODE
  if (large_holes) {
    MEDDLY_CHECK_RANGE(1, large_holes, last+1);
    MEDDLY_DCASSERT(data[large_holes] < 0);
  }
#endif
}

void MEDDLY::node_storage::indexRemove(int p_offset)
{
#ifdef MEMORY_TRACE
  printf("indexRemove(%d)\n", p_offset);
#endif

  MEDDLY_DCASSERT(!isHoleNonIndex(p_offset));
  int above = holeUp(p_offset); 
  int below = holeDown(p_offset); 
  int right = holeNext(p_offset); 

  if (right >= 1) {
    // there are nodes to the right!
    MEDDLY_DCASSERT(holeUp(right) < 0);
    holeUp(right) = above;
    holeDown(right) = below;
    // update the pointers of the holes (index) above and below it
    if (above) {
      holeDown(above) = right;
    } else {
      holes_top = right;
    }

    if (below) {
      holeUp(below) = right;
    } else {
      holes_bottom = right;
    }
    
  } else {
    // there are no non-index nodes
    MEDDLY_DCASSERT(right < 1);

    // this was the last node of its size
    // update the pointers of the holes (index) above and below it
    if (above) {
      holeDown(above) = below;
    } else {
      holes_top = below;
    }

    if (below) {
      holeUp(below) = above;
    } else {
      holes_bottom = above;
    }
  }
  // Sanity checks
#ifdef DEVELOPMENT_CODE
  if (large_holes) {
    MEDDLY_CHECK_RANGE(1, large_holes, last+1);
    MEDDLY_DCASSERT(data[large_holes] < 0);
  }
#endif
}

void MEDDLY::node_storage::resize(int new_size)
{
  int* new_data = (int *) realloc(data, new_size * sizeof(int));
  if (0 == new_data) throw error(error::INSUFFICIENT_MEMORY);
  if (0== data) new_data[0] = 0;
  if (new_size > size) {
    parent->changeStats().incMemAlloc((new_size - size) * sizeof(int));
  } else {
    parent->changeStats().decMemAlloc((size - new_size) * sizeof(int));
  }
#ifdef MEMORY_TRACE
    printf("Resized data[]. Old size: %d, New size: %d, Last: %d.\n", 
      size, new_size, last
    );
#endif
  size = new_size;
  if (data != new_data) {
    // update pointers
    data = new_data;
    data_down = data + commonHeaderLength + unhashedHeader + hashedHeader;
  }
}

int MEDDLY::node_storage::allocNode(int sz, int tail, bool clear)
{
  MEDDLY_DCASSERT(parent);
  int slots = slotsForNode(sz);
  MEDDLY_DCASSERT(slots >= commonExtra + unhashedHeader + hashedHeader);

  int off = getHole(slots);
  int got = -data[off];
  parent->changeStats().incMemUsed(got * sizeof(int));
  MEDDLY_DCASSERT(got >= slots);
  if (clear) memset(data+off, 0, slots*sizeof(int));
  setCountOf(off, 1);                     // #incoming
  setNextOf(off, temp_node_value);        // mark as a temp node
  setSizeOf(off, sz);                     // size
  data[off+slots-1] = slots - got;        // negative padding
  data[off+got-1] = tail;                 // tail entry
#ifdef MEMORY_TRACE
  printf("Allocated new node, asked %d, got %d, position %d (size %d)\n", slots, got, off, sz);
#ifdef DEEP_MEMORY_TRACE
  dumpInternal(stdout);
#else
  dumpInternal(stdout, off);
#endif
#endif
  return off;
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                     node_compacted helpers                     *
// *                                                                *
// *                                                                *
// ******************************************************************
#endif  // #ifndef TEST_ENCODING

template <int bytes, class INT>
inline void valToData(INT a, unsigned char* b)
{
  //
  // The compiler is smart enough to optimize out these if's
  //
  b[0] = a & 0xff;
  if (bytes > 1)  { a >>= 8;  b[1] = a & 0xff;  }
  if (bytes > 2)  { a >>= 8;  b[2] = a & 0xff;  }
  if (bytes > 3)  { a >>= 8;  b[3] = a & 0xff;  }
  if (bytes > 4)  { a >>= 8;  b[4] = a & 0xff;  }
  if (bytes > 5)  { a >>= 8;  b[5] = a & 0xff;  }
  if (bytes > 6)  { a >>= 8;  b[6] = a & 0xff;  }
  if (bytes > 7)  { a >>= 8;  b[7] = a & 0xff;  }
}

template <int bytes, class INT>
inline void dataToVal(const unsigned char* b, INT &a)
{
  //
  // The compiler is smart enough to optimize out these if's
  //
  if (bytes > 7)  { a |= b[7];  a <<= 8;  }
  if (bytes > 6)  { a |= b[6];  a <<= 8;  }
  if (bytes > 5)  { a |= b[5];  a <<= 8;  }
  if (bytes > 4)  { a |= b[4];  a <<= 8;  }
  if (bytes > 3)  { a |= b[3];  a <<= 8;  }
  if (bytes > 2)  { a |= b[2];  a <<= 8;  }
  if (bytes > 1)  { a |= b[1];  a <<= 8;  }
  a |= b[0];
}

template <int bytes>
inline void longToData(long L, unsigned char* d)
{
  valToData<bytes>(L, d);
}

template <int bytes>
inline void dataToLong(const unsigned char* d, long& L)
{
  // deal with negatives properly
  if (d[bytes-1] & 0x80) {
    L = (~0L) << 8;
  } else {
    L = 0;
  }
  dataToVal<bytes>(d, L);
}

template <int bytes>
inline void intToData(int L, unsigned char* d)
{
  valToData<bytes>(L, d);
}

template <int bytes>
inline void dataToInt(const unsigned char* d, int& L)
{
  // deal with negatives properly
  if (d[bytes-1] & 0x80) {
    L = (~0) << 8;
  } else {
    L = 0;
  }
  dataToVal<bytes>(d, L);
}

template <int bytes>
inline void downToData(long P, unsigned char* d)
{
  // positive P: as usual.
  if (P >= 0) {
    valToData<bytes>(P, d);
    return;
  }

  // negative P: this is a terminal pointer.
  //              next msb set - terminal value is negative.
  //
  //  No conversion necessary because msb propogates when we shift
  static const unsigned long nmsb = (0x40L) << ((sizeof(long)-1)*8);
  if (P & nmsb) {
    valToData<bytes>(P, d);
    return;
  }

  // negative P: this is a terminal pointer.
  //              next msb clr - terminal value is positive.
  //              
  //  The thing to do here is deal with the msb manually:
  //  clear msb, encode, set msb.
  static const unsigned long msboff = ~ ((0x80L) << ((sizeof(long)-1)*8));
  valToData<bytes>(P & msboff, d);
  d[bytes-1] |= 0x80;
}

template <int bytes>
inline void dataToDown(const unsigned char* d, long& P)
{
  // Is this a terminal value?
  if (d[bytes-1] & 0x80) {
    // YES.
    // Is this a negative terminal value?
    if (d[bytes-1] & 0x40) {
      // YES.
      // Easy case: same as ordinary negatives.
      P = (~0L) << 8;
      dataToVal<bytes>(d, P);
      return;
    }
    // NO.
    // Positive terminal value.
    P = 0;
    dataToVal<bytes>(d, P);
    if (bytes != 8) {
      // Move MSB
      static const unsigned long bmsboff = ~ ((0x80L) << ((bytes-1)*8));
      static const unsigned long msbon = (0x80L) << ((sizeof(long)-1)*8);
      P = (P & bmsboff) | msbon;
    }
    return;
  }

  // non-terminal value: as usual
  P = 0;
  dataToVal<bytes>(d, P);
}

template <class INT>
inline void bytesRequiredForVal(INT a, int& bytes)
{
  if (sizeof(INT) == bytes) return;
  long la = (a<0) ? -(a<<1) : a<<1; 
  static unsigned long msbyte = (0xffL) << (sizeof(INT)-1)*8;
  unsigned long mask = msbyte;
  for (int b=sizeof(INT); b>bytes; b--) {
    if (la & mask) {
      bytes = b;
      return;
    }
    mask = (mask >> 8) & ~msbyte;
  }
}

inline void bytesRequiredForDown(long a, int& bytes)
{
  if (a<0) {
    // terminal value
    a <<= 1;
    if (a<0) {
      a <<= 1;
      a = -a;
    } else {
      a <<= 1;
    }
  } else {
    a <<= 1;
  }
  static unsigned long msbyte = (0xffL) << (sizeof(long)-1)*8;
  unsigned long mask = msbyte;
  for (int b=sizeof(long); b>bytes; b--) {
    if (a & mask) {
      bytes = b;
      return;
    }
    mask = (mask >> 8) & ~msbyte;
  }
}

#ifndef TEST_ENCODING // turn off code

// ******************************************************************
// *                                                                *
// *                                                                *
// *                     node_compacted methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************

/*
    Data is stored in a char* array, so we have to copy stuff off
    (since things might not be aligned properly).
    Nodes are stored as follows (numbers are offsets into the "address").
    Hole storage details follow.

      long: >=0:  The incoming count; this is an active node.
            < 0:  This is a hole; number of bytes in the hole.

      long: >=0:  Next pointer in the unique table.
            < 0:  This node is not in the unique table yet;
                  some special value is here to indicate status.
       
      uint(32b):  Size.  Number of downward pointers stored.

      uchar(8b):  Compression info.
                  High 4 bits:  dpbytes, the number of bytes used 
                                to store downward pointers.
                                Range of possible values is 0-15, 
                                but currently, only the values 1-8 are used.

                  Low 4 bits:   ixbytes, number of bytes used to store indexes.
                                Range of possible values is 0-15,
                                but currently, only the values 0-4 are used.
                                (Max node size is 2^31-1).
                                A count of 0 means the node is stored in full;
                                otherwise, the node is stored sparsely.

      uhbytes:    Unhashed header info; might be 0 bytes.
      hhbytes:    Hashed header info; might be 0 bytes.

      dpbytes:    Down pointer [0]
        :
        :
      dpbytes:    Down pointer [size-1]


      ixbytes:    Index value [0]
        :
        :
      ixbytes:    Index value [size-1]


      evbytes:    Edge value [0]
        :
        :
      evbytes:    Edge value [size-1]

      uchar(8b):  Padding length (can be 0).
                  Number of unused extra bytes, at most 255.

                  (extra bytes here, if any)

      long: > 0:  The very last slot is a non-negative long, giving
                  the forest node number.
                  (A negative value here indicates a hole.)


      The absolute smallest number of bytes for a node,
      assuming a zero size is allowed, is therefore:
        2 longs + 1 int + 1 char + 1 char + 1 long
      

      Holes are stored using a grid.
      Currently, holes require more storage than the smallest possible
      node.  If we switch to lists, we can save one long per hole,
      and then they will be very close.

      Holes are stored as follows.

        long: < 0:  -number of bytes in the hole

        long: >=0:  This is an index hole,
                    and this is an upward pointer in the grid.
              < 0:  This is a non-index hole.  Value is ignored.

        long:       down or previous pointer.
        long:       next pointer.
          :
          :         (unused bytes)
          :
        long: < 0:  -number of bytes in the hole



*/

MEDDLY::node_compacted::node_compacted()
{
  parent = 0;
}

MEDDLY::node_compacted::~node_compacted()
{
  if (parent) parent->changeStats().decMemAlloc(size);
  free(data);
}

void MEDDLY::node_compacted::init(expert_forest* p)
{
  MEDDLY_DCASSERT(0==parent);
  parent = p;
  data = 0;
  data_down = 0;
  size = 0;
  last = 0;
  max_request = 0;
  large_holes = 0;
  holes_top = 0;
  holes_bottom = 0;
  hole_slots = 0;
  evbytes = parent->edgeBytes();
  uhbytes = parent->unhashedHeaderBytes();
  hhbytes = parent->hashedHeaderBytes();
}

void MEDDLY::node_compacted::compact(bool shrink)
{
  //
  // Should we even bother?
  //
  if (0==data || 0==hole_slots) return;

  if (hole_slots <= parent->getPolicies().compact_min)  return;

  if (hole_slots <  parent->getPolicies().compact_max) {

    // If percentage of holes is below trigger, then don't compact

    if (100 * hole_slots < last * parent->getPolicies().compact_frac) return;

  }

#ifdef DEBUG_SLOW
  fprintf(stderr, "Compacting forest level\n");
#endif
#ifdef MEMORY_TRACE
  printf("Compacting\n");
#endif
#ifdef DEBUG_COMPACTION
  printf("Before compaction:\n");
  dumpInternal(stdout);
  printf("\n");
#endif

  //
  // Scan the whole array of data, copying over itself and skipping holes.
  // 
  long node = 1;  // since we leave [0] empty
  long curr = 1;

  while (node < last) {
    long L;
    // read first slot; determine if we are a hole or not.
    dataToLong<sizeof(long)>(data+node, L);
    if (L<0) {
      //
      // This is a hole, skip it
      // 
      MEDDLY_DCASSERT(
        0==memcmp(data+node, data+node-L-sizeof(long), sizeof(long))
      );
      node -= L;
      continue;
    }

    //
    // A real node, move it
    //
    long old_off = node;
    long new_off = curr;
    MEDDLY_DCASSERT(!parent->isPessimistic() || L>0);

    int size;
    dataToInt<sizeof(int)>(data_size + node, size);

    long nodebytes = size;
    nodebytes *= dpbytes(node) + ixbytes(node) + evbytes;
    nodebytes += uhbytes + hhbytes + commonHeaderBytes;
    
    // move everything up to the padding value
    if (node != curr) {
      memmove(data + curr, data + node, nodebytes);
    }
    node += nodebytes;
    curr += nodebytes;

    // skip padding
    data[curr] = 0;
    node += data[node];
    node++;

    // copy trailer
    memmove(data + curr, data + node, commonTailBytes);
    node += commonTailBytes;
    curr += commonTailBytes;
    
    // update node header
    dataToLong<sizeof(long)>(data + curr - sizeof(long), L);
    parent->moveNodeOffset(L, old_off, new_off);
  } // while
  MEDDLY_DCASSERT(node == 1+last);
  last = curr-1;

  // set up hole pointers and such
  holes_top = holes_bottom = 0;
  hole_slots = 0;
  large_holes = 0;

  parent->changeStats().num_compactions++;

  if (shrink && size > min_size && last < size/2) {
    long new_size = size/2;
    while (new_size > min_size && new_size > last * 3) { new_size /= 2; }
    resize(new_size);
  }

#ifdef DEBUG_COMPACTION
  printf("After compaction:\n");
  dumpInternal(stdout);
  printf("\n");
#endif
}


void MEDDLY::node_compacted::dumpInternal(FILE *s) const
{
  if (0==data) return; // nothing to display
  
  fprintf(s, "Last slot used: %ld\n", last);
  fprintf(s, "large_holes: %ld\n", large_holes);
  fprintf(s, "Grid: top = %ld bottom = %ld\n", holes_top, holes_bottom);

  fprintf(s, "Data array by record: \n");
  int awidth = digits(parent->getLastNode());
  long a;
  for (a=1; a<=last; ) {
    fflush(s);
    a = dumpNode(s, a);
  } // for a
  fprintf(s, "%*ld : free slots\n", awidth, a);
  fflush(s);
  MEDDLY_DCASSERT(a == (last)+1);
}

void MEDDLY::node_compacted::dumpInternal(FILE *s, long a) const
{
  MEDDLY_DCASSERT(data);
  
  fprintf(s, "Last slot used: %ld\n", last);
  fprintf(s, "large_holes: %ld\n", large_holes);
  fprintf(s, "Grid: top = %ld bottom = %ld\n", holes_top, holes_bottom);

  dumpNode(s, a);
}


//
//
// Private
//
//

long MEDDLY::node_compacted::dumpNode(FILE* s, long a) const
{
  MEDDLY_DCASSERT(data);
  if (0==a) return 0;
  int awidth = digits(parent->getLastNode());
  fprintf(s, "%*ld : [", awidth, a);

  const unsigned char* d = data + a;
  const unsigned char* stop;

  // first slot: definitely a long
  long L;
  dataToLong<sizeof(long)>(d, L);
  fprintf(s, "%ld", L);
  bool is_hole = (L<0);

  if (is_hole) {
    //
    // Hole.
    //
    stop = d - L;
    // get next 3 longs
    for (int i=3; i; i--) {
      d += sizeof(long);
      dataToLong<sizeof(long)>(d, L);
      fprintf(s, "|%ld", L);
    }
    fprintf(s, "| ... ");
  } else {
    //
    // Node.
    //
    
    // get next: long
    d += sizeof(long);
    dataToLong<sizeof(long)>(d, L);
    fprintf(s, "|n=%ld", L);
    d += sizeof(long);

    // get size: int
    int size;
    dataToInt<sizeof(int)>(d, size);
    fprintf(s, "|s=%d", size);
    d += sizeof(int);

    // get compression info
    int dpb = dpbytesOf(d[0]);
    MEDDLY_DCASSERT(dpbytes(a) == dpb);
    int ixb = ixbytesOf(d[0]);
    MEDDLY_DCASSERT(ixbytes(a) == ixb);
    d++;
    fprintf(s, "|dpb=%d, ixb=%d", dpb, ixb);

    // show unhashed header
    if (uhbytes) {
      fprintf(s, "|uh");
      for (int i=0; i<uhbytes; i++) {
        fprintf(s, " %02x", d[0]);
        d++;
      }
    }

    // show hashed header
    if (hhbytes) {
      fprintf(s, "|hh");
      for (int i=0; i<hhbytes; i++) {
        fprintf(s, " %02x", d[0]);
        d++;
      }
    }

    // show next chunk: downward pointers
    for (int i=0; i<size; i++) {
      switch (dpb) {
        case 1: dataToDown<1>(d, L);  break;
        case 2: dataToDown<2>(d, L);  break;
        case 3: dataToDown<3>(d, L);  break;
        case 4: dataToDown<4>(d, L);  break;
        case 5: dataToDown<5>(d, L);  break;
        case 6: dataToDown<6>(d, L);  break;
        case 7: dataToDown<7>(d, L);  break;
        case 8: dataToDown<8>(d, L);  break;
        default:
          throw error(error::MISCELLANEOUS);
      }
      fprintf(s, "|d %ld", L);
      d += dpb;
    }

    // show next chunk: indexes
    if (ixb) for (int i=0; i<size; i++) {
      int I;
      switch (ixb) {
        case 1: dataToInt<1>(d, I); break;
        case 2: dataToInt<2>(d, I); break;
        case 3: dataToInt<3>(d, I); break;
        case 4: dataToInt<4>(d, I); break;
        default:
          throw error(error::MISCELLANEOUS);
      }
      fprintf(s, "|i %d", I);
      d += ixb;
    }

    // show next chunk: edge values (raw)
    if (evbytes) for (int i=0; i<size; i++) {
      fprintf(s, "|e");
      for (int b=0; b<evbytes; b++) {
        fprintf(s, " %02x", d[0]);
        d++;
      }
    }

    // get the padding
    unsigned int padding = d[0];
    d++;
    fprintf(s, "| ... %d unused ... ", padding);
    stop = d + padding + sizeof(long);
  }
  // show the last, long slot
  dataToLong<sizeof(long)>(stop-sizeof(long), L);
  fprintf(s, "|%ld]\n", L);
  return stop - (data+a);
}


void MEDDLY::node_compacted::resize(long new_size)
{
  unsigned char* new_data = (unsigned char *) realloc(data, new_size);
  if (0 == new_data) throw error(error::INSUFFICIENT_MEMORY);
  if (0== data) new_data[0] = 0;
  if (new_size > size) {
    parent->changeStats().incMemAlloc(new_size - size);
  } else {
    parent->changeStats().decMemAlloc(size - new_size);
  }
#ifdef MEMORY_TRACE
    printf("Resized data[]. Old size: %ld, New size: %ld, Last: %ld.\n", 
      size, new_size, last
    );
#endif
  size = new_size;
  if (data != new_data) {
    // update pointers
    data = new_data;
    data_size = data + 2*sizeof(long);
    data_info = data_size + sizeof(int);
    data_down = data_info + 1 + uhbytes + hhbytes;
  }
}

#endif  // ifndef TEST_ENCODING


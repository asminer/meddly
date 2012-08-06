
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
// #define VALIDATE_INCOUNTS

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_reader  methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_reader::node_reader()
{
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
    down = (int*) realloc(down, nalloc*sizeof(int));
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
  down = (int*) realloc(down, alloc * sizeof(int));
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
  // use clear() to free memory
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
  zombie_nodes = 0;
  edgeSize = parent->edgeSize();
  unhashedHeader = parent->unhashedHeaderSize();
  hashedHeader = parent->hashedHeaderSize();
}

void MEDDLY::node_storage::unlinkDown(int addr)
{
#ifdef VALIDATE_INCOUNTS
  int* downptr;
#else
  const int* downptr;
#endif
  int size;
  if (isSparse(addr)) {
    downptr = sparseDownOf(addr);
    size = sparseSizeOf(addr); 
  } else {
    downptr = fullDownOf(addr);
    size = fullSizeOf(addr);
  }
#ifdef VALIDATE_INCOUNTS
  for (int i=0; i<size; i++) {
    int temp = downptr[i];
    downptr[i] = 0;
    parent->unlinkNode(temp);
  }
#else
  for (int i=0; i<size; i++) {
    parent->unlinkNode(downptr[i]);
  }
#endif
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

  if (0==a) return;
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
  if (isFull(addr)) {
    int i;
    int stop = fullSizeOf(addr);
    const int* dn = fullDownOf(addr);
    if (nr.isFull()) {
      memcpy(nr.down, dn, stop * sizeof(int));
      int* nrdext = nr.down + stop;
      memset(nrdext, 0, (nr.size-stop) * sizeof(int));
      if (nr.edge_bytes) {
        const void* ev = fullEdgeOf(addr);
        memcpy(nr.edge, ev, stop * nr.edge_bytes);
        void* evext = (char*)nr.edge + (stop * nr.edge_bytes);
        memset(evext, 0, (nr.size-stop) * nr.edge_bytes);
      }
    } else {
      int& z = nr.nnzs;
      z = 0;
      if (nr.edge_bytes) {
        void* nev = nr.edge;
        for (i=0; i<stop; i++) {
          if (0==dn[i]) continue;
          nr.down[z] = dn[i];
          nr.index[z] = i;
          const void* ev = (const char*)fullEdgeOf(addr) + i * nr.edge_bytes;
          memcpy(nev, ev, nr.edge_bytes);
          nev = (char*)nev + nr.edge_bytes;
          z++;
        } // for i
      } else {
        for (i=0; i<stop; i++) if (dn[i]) {
          nr.down[z] = dn[i];
          nr.index[z] = i;
          z++;
        } // for i
      } // if ev
    }
  } else {
    int i = 0;
    int nnz = sparseSizeOf(addr);
    const int* dn = sparseDownOf(addr);
    const int* ix = sparseIndexesOf(addr);
    if (nr.isFull()) {
      if (nr.edge_bytes) {
        const void* ev = sparseEdgeOf(addr);
        memset(nr.down, 0, nr.size * sizeof(int));
        memset(nr.edge, 0, nr.size * nr.edge_bytes);
        for (int z=0; z<nnz; z++) {
          nr.down[ix[z]] = dn[z];
          int off = ix[z] * nr.edge_bytes;
          memcpy((char*)nr.edge + off, (char*)ev + off, nr.edge_bytes);
          i++;
        }
      } else {
        for (int z=0; z<nnz; z++) {
          for (; i<ix[z]; i++) nr.down[i] = 0;
          nr.down[i] = dn[z];
          i++;
        }
        for (; i<nr.size; i++) nr.down[i] = 0;
      } // if ev
    } else {
      nr.nnzs = nnz;
      memcpy(nr.down, dn, nnz * sizeof(int));
      memcpy(nr.index, ix, nnz * sizeof(int));
      if (nr.edge_bytes) {
        memcpy(nr.edge, sparseEdgeOf(addr), 
          nnz * nr.edge_bytes);
      }
    }
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
          if (nb.isSparse()) {
            truncsize = -1;
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
//
// Private
//
//

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



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

#include "old_scheme.h"


#define MERGE_AND_SPLIT_HOLES
// #define DEBUG_COMPACTION
// #define DEBUG_SLOW
// #define MEMORY_TRACE
// #define DEEP_MEMORY_TRACE



// ******************************************************************
// *                                                                *
// *                                                                *
// *                    old_node_storage methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::old_node_storage::old_node_storage() : node_storage()
{
}

MEDDLY::old_node_storage::~old_node_storage()
{
  decMemAlloc(size*sizeof(long));
  free(data);
}

MEDDLY::node_storage* MEDDLY::old_node_storage
::createForForest(expert_forest* f) const
{
  old_node_storage* nns = new old_node_storage;
  nns->initForForest(f);
  nns->data = 0;
  nns->size = 0;
  nns->last = 0;
  nns->max_request = 0;
  nns->large_holes = 0;
  nns->holes_top = 0;
  nns->holes_bottom = 0;
  nns->hole_slots = 0;
  nns->edgeSize = f->edgeSize();
  nns->unhashedHeader = f->unhashedHeaderSize();
  nns->hashedHeader = f->hashedHeaderSize();
  return nns;
}

void MEDDLY::old_node_storage::collectGarbage(bool shrink)
{
  //
  // Should we even bother?
  //
  if (0==data || 0==hole_slots) return;

  if (hole_slots <= getParent()->getPolicies().compact_min)  return;

  if (hole_slots <  getParent()->getPolicies().compact_max) {

    // If percentage of holes is below trigger, then don't compact

    if (100 * hole_slots < last * getParent()->getPolicies().compact_frac) 
      return;

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
  long *node_ptr = (data + 1);  // since we leave [0] empty
  long *end_ptr =  (data + last + 1);
  long *curr_ptr = node_ptr;

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
    MEDDLY_DCASSERT(!getParent()->isPessimistic() || *node_ptr != 0);
    
    long old_off = node_ptr - data;
    long new_off = curr_ptr - data;

    // copy the node, except for the tail
    int datalen = longSlotsForNode(sizeOf(old_off)) - longTailSize;
    if (node_ptr != curr_ptr) {
      memmove(curr_ptr, node_ptr, datalen * sizeof(long));
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
    moveNodeOffset(*curr_ptr, old_off, new_off);
    curr_ptr++;
    node_ptr++;

  } // while
  MEDDLY_DCASSERT(node_ptr == end_ptr);

  last = (curr_ptr - 1 - data);

  // set up hole pointers and such
  holes_top = holes_bottom = 0;
  hole_slots = 0;
  large_holes = 0;

  incCompactions(); // parent->changeStats().num_compactions++;

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

void MEDDLY::old_node_storage
::reportMemoryUsage(FILE* s, const char* pad, int verb) const
{
  if (verb<=6) return;

  fprintf(s, "%s", pad);
#ifdef NODE_STORAGE_PER_LEVEL
  fprintf(s, "Level %ld ", parent->getLevelNumber(this));
#endif
  fprintf(s, "Memory stats for node storage\n");

  if (verb>7) {
    long holemem = hole_slots * sizeof(long);
    fprintf(s, "%s  Hole Memory Usage:\t%ld\n", pad, holemem);
  }

  // Compute chain length histogram
  std::map<int, int> chainLengths;
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

  // Display the histogram

  fprintf(s, "%s  Hole Chains (size, count):\n", pad);
  for (std::map<int, int>::iterator iter = chainLengths.begin();
      iter != chainLengths.end(); ++iter)
    {
      if (iter->first<0)
        fprintf(s, "%s\tlarge: %d\n", pad, iter->second);
      else
        fprintf(s, "%s\t%5d: %d\n", pad, iter->first, iter->second);
    }
  fprintf(s, "%s  End of Hole Chains\n", pad);
}

void MEDDLY::old_node_storage::showNode(FILE* s, long addr, bool verb) const
{
  if (sizeOf(addr) < 0) {
    // Sparse node
    int nnz = -sizeOf(addr);
    if (verb) fprintf(s, " nnz : %d", nnz);
    fprintf(s, " down: (");
    for (int z=0; z<nnz; z++) {
      if (z) fprintf(s, ", ");
      fprintf(s, "%d:", SI(addr)[z]);
      if (edgeSize) {
        fprintf(s, "<");
        getParent()->showEdgeValue(s, SEP(addr, z), 0);
        fprintf(s, ", ");
      } 
      int d = SD(addr)[z];
      if (getParent()->isTerminalNode(d)) {
        getParent()->showTerminal(s, d);
      } else {
        fprintf(s, "%d", d);
      }
      if (edgeSize) fprintf(s, ">");
    } // for z
    fprintf(s, ")");
  } else {
    // Full node
    int size = sizeOf(addr);
    if (verb) fprintf(s, " size: %d", size);
    fprintf(s, " down: [");
    for (int i=0; i<size; i++) {
      if (i) fprintf(s, "|"); 
      if (edgeSize) {
        fprintf(s, "<");
        getParent()->showEdgeValue(s, FEP(addr, i), 0);
        fprintf(s, ", ");
      } 
      int d = FD(addr)[i];
      if (getParent()->isTerminalNode(d)) {
        getParent()->showTerminal(s, d);
      } else {
        fprintf(s, "%d", d);
      }
      if (edgeSize) fprintf(s, ">");
    } // for i
    fprintf(s, "]");
  }

  // show extra header stuff
  if (unhashedHeader) {
    getParent()->showUnhashedHeader(s, UH(addr));
  }
  if (hashedHeader) {
    getParent()->showHashedHeader(s, HH(addr));
  }
}

long MEDDLY::old_node_storage
::makeNode(long p, const node_builder &nb, node_storage_flags opt)
{
  int nnz, truncsize;

  //
  // Easy case - sparse nodes disabled
  //
  if (0==forest::policies::ALLOW_SPARSE_STORAGE & opt) {

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
  }

  //
  // Easy case - full nodes disabled
  //
  if (0==forest::policies::ALLOW_FULL_STORAGE & opt) {

          if (nb.isSparse()) {
            return makeSparseNode(p, nb.getNNZs(), nb);
          } else {
            nnz = 0;
            for (int i=0; i<nb.getSize(); i++) {
              if (nb.d(i)) nnz++;
            }
            return makeSparseNode(p, nnz, nb);
          }
  }

  //
  // Full and sparse are allowed, determine which is more compact
  //
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
  if (longSlotsForNode(-nnz) < longSlotsForNode(truncsize)) { 
    return makeSparseNode(p, nnz, nb);
  } else {
    return makeFullNode(p, truncsize, nb);
  }

}

void MEDDLY::old_node_storage::unlinkDownAndRecycle(long addr)
{
  //
  // Unlink down pointers
  //
  int* down;
  int size = sizeOf(addr);
  if (size < 0) {
    size = -size;
    down = SD(addr);
  } else {
    down = FD(addr);
  }
  for (int i=0; i<size; i++) {
    getParent()->unlinkNode(down[i]);
  }

  //
  // Recycle
  //
  makeHole(addr, activeNodeActualLongSlots(addr));
}


void MEDDLY::old_node_storage::fillReader(long addr, node_reader &nr) const
{
  // Copy hashed header

  if (hashedHeader) {
    resize_header(nr, hashedHeader);
    memcpy(extra_hashed(nr), HH(addr), hashedHeader * sizeof(int));
  }

  // Copy everything else

  int size = sizeOf(addr);
  if (size < 0) {
    //
    // Node is sparse
    //
    int nnz = -size;
    int* down = SD(addr);
    int* index = SI(addr);
    if (nr.isFull()) {
      if (nr.edgeBytes()) {
        memset(down_of(nr), 0, nr.getSize() * sizeof(long));
        memset(edge_of(nr), 0, nr.getSize() * nr.edgeBytes());
        for (int z=0; z<nnz; z++) {
          int i = index[z];
          down_of(nr)[i] = down[z];
          int off = i * nr.edgeBytes();
          memcpy(edge_of(nr) + off, SEP(addr, z), nr.edgeBytes());
        } // for z
      } else {
        memset(down_of(nr), 0, nr.getSize() * sizeof(long));
        for (int z=0; z<nnz; z++) {
          int i = index[z];
          down_of(nr)[i] = down[z];
        }
      } // if ev
      return;
    }

    // nr is sparse
    nnzs_of(nr) = nnz;
    for (int z=0; z<nnz; z++) {
      down_of(nr)[z] = down[z];
      index_of(nr)[z] = index[z];
    } // for z
    if (nr.edgeBytes()) {
      memcpy(edge_of(nr), SE(addr), nnz * nr.edgeBytes());
    }
    return;
  } 
  //
  // Node is full
  //
  int* down = FD(addr);
  if (nr.isFull()) {
    int i;
    for (i=0; i<size; i++) {
      down_of(nr)[i] = down[i];
    } 
    for (; i<nr.getSize(); i++) {
      down_of(nr)[i] = 0;
    }
    if (edgeSize) {
      long bytes = size * nr.edgeBytes();
      memcpy(edge_of(nr), FE(addr), bytes);
      void* evext = edge_of(nr) + bytes;
      memset(evext, 0, (nr.getSize()-size) * nr.edgeBytes());
    }
    return;
  }
  // nr is sparse
  int& z = nnzs_of(nr);
  z = 0;
  if (nr.edgeBytes()) {
    char* nev = edge_of(nr);
    for (int i=0; i<size; i++) if (down[i]) {
      down_of(nr)[z] = down[i];
      index_of(nr)[z] = i;
      memcpy(nev, FEP(addr, i), nr.edgeBytes());
      nev = nev + nr.edgeBytes();
      z++;
    } // for i
  } else {
    for (int i=0; i<size; i++) if (down[i]) {
      down_of(nr)[z] = down[i];
      index_of(nr)[z] = i;
      z++;
    } // for i
  } // if ev
  return; 
}


unsigned MEDDLY::old_node_storage::hashNode(const node_header& node) const
{
  hash_stream s;
  s.start(node.level);

  // Do the hashed header part, if any
  if (hashedHeader) {
    int* hhptr = HH(node.offset);
    for (int e=0; e<hashedHeader; e++) {
      s.push(hhptr[e]);
    }
  }

  //
  // Hash the node itself
  
  int size = sizeOf(node.offset);
  if (size < 0) {
    // Node is sparse
    int nnz = -size;
    int* down = SD(node.offset);
    int* index = SI(node.offset);
    if (getParent()->areEdgeValuesHashed()) {
      for (int z=0; z<nnz; z++) {
        s.push(index[z], down[z], ((int*)SEP(node.offset, z))[0]);
      } // for z
    } else {
      for (int z=0; z<nnz; z++) {
        s.push(index[z], down[z]);
      } // for z
    }
  } else {
    // Node is full
    int* down = FD(node.offset);
    if (getParent()->areEdgeValuesHashed()) {
      for (int i=0; i<size; i++) if (down[i]) {
        s.push(i, down[i], ((int*)FEP(node.offset, i))[0]);
      } // for z
    } else {
      for (int i=0; i<size; i++) if (down[i]) {
        s.push(i, down[i]);
      } // for z
    }
  }

  return s.finish();
}


int MEDDLY::old_node_storage::getSingletonIndex(long addr, long &down) const
{
  int size = sizeOf(addr);
  if (size<0) {
    // sparse node --- easy
    if (size != -1) return -1;
    down = SD(addr)[0];
    return SI(addr)[0];
  } 
  
  // full node
  int* dnptr = FD(addr);
  for (int i=0; i<size; i++) {
    if (0==dnptr[i]) continue;
    if (i+1 != size) return -1;
    down = dnptr[i];
    return i;
  }
  return -1;
}


long MEDDLY::old_node_storage
::getDownPtr(long addr, int index) const
{
  int size = sizeOf(addr);
  if (size<0) {
    int z = findSparseIndex(addr, index);
    if (z<0) return 0;
    return SD(addr)[z];
  } else {
    if (index < size) {
      return FD(addr)[index];
    } else {
      return 0;
    }
  }
}

void MEDDLY::old_node_storage
::getDownPtr(long addr, int index, int& ev, long& dn) const
{
  if (sizeOf(addr)<0) {
    int z = findSparseIndex(addr, index);
    if (z<0) {
      dn = 0;
      ev = 0;
    } else {
      dn = SD(addr)[z];
      ev = ((int*)SEP(addr, z))[0];
    }
  } else {
    if (index < size) {
      dn = FD(addr)[index];
      ev = ((int*)FEP(addr, index))[0];
    } else {
      dn = 0;
      ev = 0;
    }
  }
}

void MEDDLY::old_node_storage
::getDownPtr(long addr, int index, float& ev, long& dn) const
{
  if (sizeOf(addr)<0) {
    int z = findSparseIndex(addr, index);
    if (z<0) {
      dn = 0;
      ev = 0;
    } else {
      dn = SD(addr)[z];
      ev = ((float*)SEP(addr, z))[0];
    }
  } else {
    if (index < size) {
      dn = FD(addr)[index];
      ev = ((float*)FEP(addr, index))[0];
    } else {
      dn = 0;
      ev = 0;
    }
  }
}


int MEDDLY::old_node_storage::getUnhashedHeaderOf(long addr, int ind) const
{
  return UH(addr)[ind];
}

int MEDDLY::old_node_storage::getHashedHeaderOf(long addr, int ind) const
{
  return HH(addr)[ind];
}




void MEDDLY::old_node_storage::dumpInternalInfo(FILE* s) const
{
  if (0==data) return; // nothing to display
#ifdef NODE_STORAGE_PER_LEVEL
  fprintf(s, "Level %ld: ", getParent()->getLevelNumber(this));
#endif
  fprintf(s, "Last slot used: %d\n", last);
  fprintf(s, "large_holes: %d\n", large_holes);
  fprintf(s, "Grid: top = %d bottom = %d\n", holes_top, holes_bottom);
}


long MEDDLY::old_node_storage::dumpInternalNode(FILE *s, long a) const
{
  MEDDLY_DCASSERT(data);
  if (a<=0) return 0;
  int awidth = digits(getParent()->getLastNode());
  if (a > last) {
    fprintf(s, "%*ld : free slots\n", awidth, a);
    return 0;
  }
  fprintf(s, "%*ld : [%ld|", awidth, a, data[a]);
  if (data[a]<0) {
    // hole
    const int* hole = holeOf(a);
    for (int i=hole_up_index; i<=hole_next_index; i++) {
      fprintf(s, "%d", hole[i]);
      if ((i+1) % slots_per_long) fputc(',', s); 
      else fputc('|', s);
    }
    fprintf(s, "| ... ");
    a -= data[a];  
  } else {
    // not a hole
    const int* node = chunkOf(a);
    fprintf(s, "%ld|", data[a+1]);  // next
    int slots = activeNodeActualLongSlots(a);
    int nLongs = slots - longTailSize;
    int i=size_index;
    for (int L=size_index / slots_per_long; L<nLongs; L++) {
      fprintf(s, "%d,%d|", node[i], node[i+1]);
      i += 2;
    }
    MEDDLY_DCASSERT(0==slots % slots_per_long);
    a += slots;
  }
  fprintf(s, "%ld]\n", data[a-1]);
  return a;
}








//
//
// Private
//
//

long MEDDLY::old_node_storage
::makeFullNode(long p, int size, const node_builder &nb)
{
  long addr = allocNode(size, p, false);
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

long MEDDLY::old_node_storage
::makeSparseNode(long p, int size, const node_builder &nb)
{
  long addr = allocNode(-size, p, false);
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


void MEDDLY::old_node_storage
::copyExtraHeader(long addr, const node_builder &nb)
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
int MEDDLY::old_node_storage::getHole(int slots)
// TBD
{
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
    const int min_node_size = longSlotsForNode(0);
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
  if (getParent()->getPolicies().compactBeforeExpand) {
    if (last + slots >= size) collectGarbage(false);
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

void MEDDLY::old_node_storage::makeHole(int addr, int slots)
// TBD
{
#ifdef MEMORY_TRACE
  printf("Calling makeHole(%d, %d)\n", addr, slots);
#endif

  decMemUsed(slots * sizeof(long));

  hole_slots += slots;
  data[addr] = data[addr+slots-1] = -slots;

  if (!getParent()->getPolicies().recycleNodeStorageHoles) return;

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

void MEDDLY::old_node_storage::gridInsert(int p_offset)
// TBD
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

void MEDDLY::old_node_storage::midRemove(int p_offset)
// TBD
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

void MEDDLY::old_node_storage::indexRemove(int p_offset)
// TBD
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

void MEDDLY::old_node_storage::resize(int new_size)
{
  long* new_data = (long*) realloc(data, new_size * sizeof(long));
  if (0 == new_data) throw error(error::INSUFFICIENT_MEMORY);
  if (0== data) new_data[0] = 0;
  if (new_size > size) {
    incMemAlloc((new_size - size) * sizeof(long));
  } else {
    decMemAlloc((size - new_size) * sizeof(long));
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
    updateCountArray(data + count_index);
    updateNextArray(data + next_index);
  }
}

long MEDDLY::old_node_storage::allocNode(int sz, long tail, bool clear)
{
  int slots = longSlotsForNode(sz);
  MEDDLY_DCASSERT(2*slots >= intExtra + unhashedHeader + hashedHeader);

  int off = getHole(slots);
  int got = -data[off];
  incMemUsed(got * sizeof(long));
  MEDDLY_DCASSERT(got >= slots);
  if (clear) memset(data+off, 0, slots*sizeof(long));
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
// *                   front-end global variables                   *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
  // node storage mechanism used for versions < 0.10 of the library
  old_node_storage THE_OLD_NODE_STORAGE;

  const node_storage* OLD_STORAGE = &THE_OLD_NODE_STORAGE;
};


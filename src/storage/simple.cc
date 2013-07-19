
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

#include "simple.h"

#include "hm_grid.h"
#include "hm_array.h"
#include "hm_heap.h"
#include "hm_none.h"

#define DEBUG_ENCODING
// #define DEBUG_COMPACTION
// #define DEBUG_SLOW
// #define MEMORY_TRACE
// #define DEEP_MEMORY_TRACE

inline char slotsForBytes(int bytes) 
{
  int sl = bytes / sizeof(MEDDLY::node_handle);
  if (bytes % sizeof(MEDDLY::node_handle)) sl++;
  return sl;
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                     simple_storage methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_storage::simple_storage(holeman* hm) : node_storage()
{
  holeManager = hm;
  data = 0;
}

MEDDLY::simple_storage::~simple_storage()
{
  delete holeManager;
}

void MEDDLY::simple_storage::collectGarbage(bool shrink)
{
  //
  // Should we even bother?
  //
  node_handle wasted = holeManager->holeSlots();
  if (0==data || 0==wasted) return;
  if (wasted <= getParent()->getPolicies().compact_min) {
    return;
  }
  if (wasted <  getParent()->getPolicies().compact_max) {

    // If percentage of wasted slots is below trigger, then don't compact
    if (100 * wasted < 
        holeManager->lastSlot() * getParent()->getPolicies().compact_frac) 
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
  node_handle* node_ptr = (data + 1);  // since we leave [0] empty
  node_handle* end_ptr =  (data + holeManager->lastSlot() + 1);
  node_handle* curr_ptr = node_ptr;

  while (node_ptr < end_ptr) {
    if (*node_ptr < 0) {
      //
      // This is a hole, skip it
      // 
      node_ptr = data + holeManager->chunkAfterHole(node_ptr - data);
      continue;
    } 
    //
    // A real node, move it
    //
    MEDDLY_DCASSERT(!getParent()->isPessimistic() || *node_ptr != 0);
    
    long old_off = node_ptr - data;
    long new_off = curr_ptr - data;

    // copy the node, except for the tail
    int datalen = slotsForNode(sizeOf(old_off)) - 1;
    if (node_ptr != curr_ptr) {
      memmove(curr_ptr, node_ptr, datalen * sizeof(node_handle));
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

  holeManager->clearHolesAndShrink( (curr_ptr - 1 - data), shrink );
  MEDDLY_DCASSERT(0==holeManager->holeSlots());

  incCompactions(); 

#ifdef DEBUG_COMPACTION
  printf("After compaction:\n");
  dumpInternal(stdout);
  printf("\n");
#endif
}

void MEDDLY::simple_storage
::reportStats(FILE* s, const char* pad, unsigned flags) const
{
  static unsigned STORAGE = 
    expert_forest::STORAGE_STATS | expert_forest::STORAGE_DETAILED;

  if (flags & STORAGE) {
    fprintf(s, "%s", pad);
#ifdef NODE_STORAGE_PER_LEVEL
    fprintf(s, "Level %ld ", parent->getLevelNumber(this));
#endif
    fprintf(s, "Stats for %s\n", getStorageName());

    // anything for us?
  }

  holeManager->reportStats(s, pad, flags);

#ifdef DEVELOPMENT_CODE
  verifyStats();
#endif
}

void MEDDLY::simple_storage::showNode(FILE* s, node_address addr, bool verb) const
{
  if (sizeOf(addr) < 0) {
    // Sparse node
    int nnz = -sizeOf(addr);
    if (verb) fprintf(s, " nnz : %d", nnz);
    fprintf(s, " down: (");
    for (int z=0; z<nnz; z++) {
      if (z) fprintf(s, ", ");
      fprintf(s, "%ld:", long(SI(addr)[z]));
      if (edgeSlots) {
        fprintf(s, "<");
        getParent()->showEdgeValue(s, SEP(addr, z), 0);
        fprintf(s, ", ");
      } 
      node_handle d = SD(addr)[z];
      if (getParent()->isTerminalNode(d)) {
        getParent()->showTerminal(s, d);
      } else {
        fprintf(s, "%ld", long(d));
      }
      if (edgeSlots) fprintf(s, ">");
    } // for z
    fprintf(s, ")");
  } else {
    // Full node
    int size = sizeOf(addr);
    if (verb) fprintf(s, " size: %d", size);
    fprintf(s, " down: [");
    for (int i=0; i<size; i++) {
      if (i) fprintf(s, "|"); 
      if (edgeSlots) {
        fprintf(s, "<");
        getParent()->showEdgeValue(s, FEP(addr, i), 0);
        fprintf(s, ", ");
      } 
      node_handle d = FD(addr)[i];
      if (getParent()->isTerminalNode(d)) {
        getParent()->showTerminal(s, d);
      } else {
        fprintf(s, "%ld", long(d));
      }
      if (edgeSlots) fprintf(s, ">");
    } // for i
    fprintf(s, "]");
  }

  // show extra header stuff
  if (unhashedSlots) {
    getParent()->showUnhashedHeader(s, UH(addr));
  }
  if (hashedSlots) {
    getParent()->showHashedHeader(s, HH(addr));
  }
}

MEDDLY::node_address MEDDLY::simple_storage
::makeNode(node_handle p, const node_builder &nb, node_storage_flags opt)
{
#ifdef DEBUG_ENCODING
  printf("simple_storage making node\n        temp:  ");
  nb.show(stdout, true);
#endif
  int nnz, truncsize;

  //
  // Easy case - sparse nodes disabled
  //
  if (0==(forest::policies::ALLOW_SPARSE_STORAGE & opt)) {

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
  if (0==(forest::policies::ALLOW_FULL_STORAGE & opt)) {

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
  if (slotsForNode(-nnz) < slotsForNode(truncsize)) { 
    return makeSparseNode(p, nnz, nb);
  } else {
    return makeFullNode(p, truncsize, nb);
  }

}

void MEDDLY::simple_storage::unlinkDownAndRecycle(node_address addr)
{
  //
  // Unlink down pointers
  //
  const node_handle* down;
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
  holeManager->recycleChunk(addr, activeNodeActualSlots(addr));
}

bool MEDDLY::simple_storage
::areDuplicates(node_address addr, const node_builder &nb) const
{
  return areDupsTempl(addr, nb);
}

bool MEDDLY::simple_storage
::areDuplicates(node_address addr, const node_reader &nr) const
{
  return areDupsTempl(addr, nr);
}

void MEDDLY::simple_storage::fillReader(node_address addr, node_reader &nr) const
{
#ifdef DEBUG_ENCODING
  printf("simple_storage filling reader\n    internal: ");
  dumpInternalNode(stdout, addr, 0x03);
  printf("        node: ");
  showNode(stdout, addr, true);
#endif
  /*
  // Copy hashed header

  if (hashedSlots) {
    resize_header(nr, hashedHeader);
    memcpy(extra_hashed(nr), HH(addr), hashedHeader * sizeof(int));
  }
  */

  // Copy everything else

  int size = sizeOf(addr);
  if (size < 0) {
    //
    // Node is sparse
    //
    int nnz = -size;
    node_handle* down = SD(addr);
    node_handle* index = SI(addr);
    if (nr.isFull()) {
      if (nr.hasEdges()) {
        memset(down_of(nr), 0, nr.getSize() * sizeof(node_handle));
        memset(edge_of(nr), 0, nr.getSize() * nr.edgeBytes());
        for (int z=0; z<nnz; z++) {
          int i = index[z];
          down_of(nr)[i] = down[z];
          int off = i * nr.edgeBytes();
          memcpy(edge_of(nr) + off, SEP(addr, z), nr.edgeBytes());
        } // for z
      } else {
        memset(down_of(nr), 0, nr.getSize() * sizeof(node_handle));
        for (int z=0; z<nnz; z++) {
          int i = index[z];
          down_of(nr)[i] = down[z];
        }
      } // if ev
#ifdef DEBUG_ENCODING
      printf("\n        temp:  ");
      nr.show(stdout, getParent(), true);
      printf("\n");
#endif
      return;
    }

    // nr is sparse
    nnzs_of(nr) = nnz;
    for (int z=0; z<nnz; z++) {
      down_of(nr)[z] = down[z];
      index_of(nr)[z] = index[z];
    } // for z
    if (nr.hasEdges()) {
      memcpy(edge_of(nr), SE(addr), nnz * nr.edgeBytes());
    }
#ifdef DEBUG_ENCODING
    printf("\n        temp:  ");
    nr.show(stdout, getParent(), true);
    printf("\n");
#endif
    return;
  } 
  //
  // Node is full
  //
  node_handle* down = FD(addr);
  if (nr.isFull()) {
    int i;
    for (i=0; i<size; i++) {
      down_of(nr)[i] = down[i];
    } 
    for (; i<nr.getSize(); i++) {
      down_of(nr)[i] = 0;
    }
    if (edgeSlots) {
      long bytes = size * nr.edgeBytes();
      memcpy(edge_of(nr), FE(addr), bytes);
      void* evext = edge_of(nr) + bytes;
      memset(evext, 0, (nr.getSize()-size) * nr.edgeBytes());
    }
#ifdef DEBUG_ENCODING
    printf("\n        temp:  ");
    nr.show(stdout, getParent(), true);
    printf("\n");
#endif
    return;
  }
  // nr is sparse
  int& z = nnzs_of(nr);
  z = 0;
  if (nr.hasEdges()) {
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
#ifdef DEBUG_ENCODING
  printf("\n        temp:  ");
  nr.show(stdout, getParent(), true);
  printf("\n");
#endif
  return; 
}


unsigned MEDDLY::simple_storage::hashNode(const node_header& node) const
{
  hash_stream s;
  s.start(node.level);

  // Do the hashed header part, if any
  /*
  if (hashedHeader) {
    int* hhptr = HH(node.offset);
    for (int e=0; e<hashedHeader; e++) {
      s.push(hhptr[e]);
    }
  }
  */

  //
  // Hash the node itself
  
  int size = sizeOf(node.offset);
  if (size < 0) {
    // Node is sparse
    int nnz = -size;
    node_handle* down = SD(node.offset);
    node_handle* index = SI(node.offset);
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
    node_handle* down = FD(node.offset);
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


int MEDDLY::simple_storage
::getSingletonIndex(node_address addr, node_handle &down) const
{
  int size = sizeOf(addr);
  if (size<0) {
    // sparse node --- easy
    if (size != -1) return -1;
    down = SD(addr)[0];
    return SI(addr)[0];
  } 
  
  // full node
  const node_handle* dnptr = FD(addr);
  for (int i=0; i<size; i++) {
    if (0==dnptr[i]) continue;
    if (i+1 != size) return -1;
    down = dnptr[i];
    return i;
  }
  return -1;
}


MEDDLY::node_handle 
MEDDLY::simple_storage
::getDownPtr(node_address addr, int index) const
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

void MEDDLY::simple_storage
::getDownPtr(node_address addr, int index, int& ev, node_handle& dn) const
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
    if (index < sizeOf(addr)) {
      dn = FD(addr)[index];
      ev = ((int*)FEP(addr, index))[0];
    } else {
      dn = 0;
      ev = 0;
    }
  }
}

void MEDDLY::simple_storage
::getDownPtr(node_address addr, int index, float& ev, node_handle& dn) const
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
    if (index < sizeOf(addr)) {
      dn = FD(addr)[index];
      ev = ((float*)FEP(addr, index))[0];
    } else {
      dn = 0;
      ev = 0;
    }
  }
}


const void* MEDDLY::simple_storage
::getUnhashedHeaderOf(node_address addr) const
{
  return UH(addr);
}

const void* MEDDLY::simple_storage
::getHashedHeaderOf(node_address addr) const
{
  return HH(addr);
}


//
//
// Protected
//
//

void MEDDLY::simple_storage::localInitForForest(const expert_forest* f)
{
  edgeSlots = slotsForBytes(f->edgeBytes());
  unhashedSlots = slotsForBytes(f->unhashedHeaderBytes());
  hashedSlots = slotsForBytes(f->hashedHeaderBytes());
  MEDDLY_DCASSERT(holeManager);
  holeManager->setParent(this);
}

void MEDDLY::simple_storage::updateData(node_handle* d)
{
  data = d;
  updateCountArray(data + count_index);
  updateNextArray(data + next_index);
}

int MEDDLY::simple_storage::smallestNode() const
{
  return slotsForNode(0);
}

void MEDDLY::simple_storage::dumpInternalInfo(FILE* s) const
{
#ifdef NODE_STORAGE_PER_LEVEL
  fprintf(s, "Level %ld: ", getParent()->getLevelNumber(this));
#endif
  holeManager->dumpInternalInfo(s);
}

void MEDDLY::simple_storage::dumpInternalTail(FILE* s) const
{
  holeManager->dumpInternalTail(s);
}



MEDDLY::node_address 
MEDDLY::simple_storage
::dumpInternalNode(FILE *s, node_address a, unsigned flags) const
{
  if (a<=0) return 0;
  int awidth = digits(getParent()->getLastNode());
  if (a > holeManager->lastSlot()) {
    fprintf(s, "%*ld : free slots\n", awidth, long(a));
    return 0;
  }
  MEDDLY_DCASSERT(data);
  if (data[a]<0) { 
    // hole
    if (flags & 0x02) {
      fprintf(s, "%*ld : ", awidth, a);
      holeManager->dumpHole(s, a);
    }
    a = holeManager->chunkAfterHole(a);
  } else {
    // proper node
    if (flags & 0x01) {
      fprintf(s, "%*ld : ", awidth, a);
      int nElements = activeNodeActualSlots(a);
      fprintf(s, "[%ld", long(data[a]));
      for (int i=1; i<nElements; i++) {
        fprintf(s, "|%ld", long(data[a+i]));
      }
      fprintf(s, "]\n");
    }
    a += activeNodeActualSlots(a);
  }
  return a;
}








//
//
// Private
//
//

MEDDLY::node_handle MEDDLY::simple_storage
::makeFullNode(node_handle p, int size, const node_builder &nb)
{
#if 0
  node_address addr = allocNode(size, p, false);
  MEDDLY_DCASSERT(1==getCountOf(addr));
  node_handle* down = FD(addr);
  if (edgeSlots) {
      MEDDLY_DCASSERT(nb.hasEdges());
      char* edge = (char*) FE(addr);
      int edge_bytes = edgeSlots * sizeof(node_handle);
      if (nb.isSparse()) {
        memset(down, 0, size * sizeof(node_handle));
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
        memset(down, 0, size * sizeof(node_handle));
        for (int z=0; z<nb.getNNZs(); z++) {
          MEDDLY_CHECK_RANGE(0, nb.i(z), size);
          down[nb.i(z)] = nb.d(z);
        }
      } else {
        for (int i=0; i<size; i++) down[i] = nb.d(i);
      }
  }
#else
  node_address addr = allocNode(size, p, true);
  MEDDLY_DCASSERT(1==getCountOf(addr));
  node_handle* down = FD(addr);
  if (edgeSlots) {
      MEDDLY_DCASSERT(nb.hasEdges());
      char* edge = (char*) FE(addr);
      int edge_bytes = edgeSlots * sizeof(node_handle);
      if (nb.isSparse()) {
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
        for (int z=0; z<nb.getNNZs(); z++) {
          MEDDLY_CHECK_RANGE(0, nb.i(z), size);
          down[nb.i(z)] = nb.d(z);
        }
      } else {
        for (int i=0; i<size; i++) down[i] = nb.d(i);
      }
  }
#endif
  copyExtraHeader(addr, nb);
#ifdef DEBUG_ENCODING
  printf("\n        made: ");
  showNode(stdout, addr, true);
  printf("\n    internal: ");
  dumpInternalNode(stdout, addr, 0x03);
#endif
  return addr;
}

MEDDLY::node_handle MEDDLY::simple_storage
::makeSparseNode(node_handle p, int size, const node_builder &nb)
{
  node_address addr = allocNode(-size, p, false);
  MEDDLY_DCASSERT(1==getCountOf(addr));
  node_handle* index = SI(addr);
  node_handle* down  = SD(addr);
  if (nb.hasEdges()) {
      MEDDLY_DCASSERT(nb.hasEdges());
      char* edge = (char*) SE(addr);
      int edge_bytes = edgeSlots * sizeof(node_handle);
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
#ifdef DEBUG_ENCODING
  printf("\n        made: ");
  showNode(stdout, addr, true);
  printf("\n    internal: ");
  dumpInternalNode(stdout, addr, 0x03);
#endif
  return addr;
}


void MEDDLY::simple_storage
::copyExtraHeader(node_address addr, const node_builder &nb)
{
  // copy extra header info, if any
  if (unhashedSlots) {
    nb.getUH(UH(addr));
  }
  if (hashedSlots) {
    nb.getHH(HH(addr));
  }
}


MEDDLY::node_handle 
MEDDLY::simple_storage::allocNode(int sz, node_handle tail, bool clear)
{
  int slots = slotsForNode(sz);
  MEDDLY_DCASSERT(slots >= extraSlots + unhashedSlots + hashedSlots);

  node_handle off = holeManager->requestChunk(slots);
  node_handle got = -data[off];
  incMemUsed(got * sizeof(node_handle));
  MEDDLY_DCASSERT(got >= slots);
  if (clear) memset(data+off, 0, slots*sizeof(node_handle));
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


#ifdef DEVELOPMENT_CODE
void MEDDLY::simple_storage::verifyStats() const
{
  int holes = 0;
  node_address hole_count = 0;
  for (node_address a=1; a<=holeManager->lastSlot(); ) {
    if (data[a]<0) {
      // hole
      node_address anext = holeManager->chunkAfterHole(a);
      hole_count += (anext - a);
      holes++;
      a = anext;
      continue;
    }
    // not hole
    a += activeNodeActualSlots(a);
  }
  // done scan, compare stats

  if ((hole_count == holeManager->holeSlots()) &&
     (holes == holeManager->numHoles())
     )   return;

  printf("Counted holes: %d stat: %d\n", 
    holes, holeManager->numHoles()
  );
  
  printf("Counted hole slots: %ld stat: %ld\n", 
    hole_count, holeManager->holeSlots()
  );

  dumpInternal(stdout, 0x03);
  
  MEDDLY_DCASSERT(0);
}
#endif

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      simple_grid  methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_grid::simple_grid(holeman* hm) : simple_storage(hm)
{
}

MEDDLY::simple_grid::~simple_grid()
{
}

MEDDLY::node_storage* MEDDLY::simple_grid
::createForForest(expert_forest* f) const
{
  simple_storage* nns = new simple_grid(new hm_grid);
  nns->initForForest(f);
  return nns;
}

const char* MEDDLY::simple_grid::getStorageName() const
{
  return "simple node storage with grid for holes";
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                      simple_array methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_array::simple_array(holeman* hm) : simple_storage(hm)
{
}

MEDDLY::simple_array::~simple_array()
{
}

MEDDLY::node_storage* MEDDLY::simple_array
::createForForest(expert_forest* f) const
{
  simple_storage* nns = new simple_array(new hm_array);
  nns->initForForest(f);
  return nns;
}

const char* MEDDLY::simple_array::getStorageName() const
{
  return "simple node storage with array of lists for holes";
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      simple_heap  methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_heap::simple_heap(holeman* hm) : simple_storage(hm)
{
}

MEDDLY::simple_heap::~simple_heap()
{
}

MEDDLY::node_storage* MEDDLY::simple_heap
::createForForest(expert_forest* f) const
{
  simple_storage* nns = new simple_heap(new hm_heap);
  nns->initForForest(f);
  return nns;
}

const char* MEDDLY::simple_heap::getStorageName() const
{
  return "simple node storage with heaps for holes";
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                      simple_none  methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_none::simple_none(holeman* hm) : simple_storage(hm)
{
}

MEDDLY::simple_none::~simple_none()
{
}

MEDDLY::node_storage* MEDDLY::simple_none
::createForForest(expert_forest* f) const
{
  simple_storage* nns = new simple_none(new hm_none);
  nns->initForForest(f);
  return nns;
}

const char* MEDDLY::simple_none::getStorageName() const
{
  return "simple node storage with no hole management";
}


// ******************************************************************
// *                                                                *
// *                   front-end global variables                   *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
  simple_grid THE_SIMPLE_GRID(0);
  simple_array THE_SIMPLE_ARRAY(0);
  simple_heap THE_SIMPLE_HEAP(0);
  simple_none THE_SIMPLE_NONE(0);

  const node_storage* SIMPLE_GRID = &THE_SIMPLE_GRID;
  const node_storage* SIMPLE_ARRAY = &THE_SIMPLE_ARRAY;
  const node_storage* SIMPLE_HEAP = &THE_SIMPLE_HEAP;
  const node_storage* SIMPLE_NONE = &THE_SIMPLE_NONE;
};


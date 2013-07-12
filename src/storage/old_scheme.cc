
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

#include "hm_grid.h"

#define MERGE_AND_SPLIT_HOLES
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

MEDDLY::node_handle MEDDLY::old_node_storage::verify_hole_slots;

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
  decMemAlloc(size*sizeof(node_handle));
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
  nns->fragment_slots = 0;
  nns->edgeSlots = slotsForBytes(f->edgeBytes());
  nns->unhashedSlots = slotsForBytes(f->unhashedHeaderBytes());
  nns->hashedSlots = slotsForBytes(f->hashedHeaderBytes());
  return nns;
}

// original version
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
  node_handle* node_ptr = (data + 1);  // since we leave [0] empty
  node_handle* end_ptr =  (data + last + 1);
  node_handle* curr_ptr = node_ptr;

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

  last = (curr_ptr - 1 - data);

  // set up hole pointers and such
  holes_top = holes_bottom = 0;
  hole_slots = 0;
  fragment_slots = 0;
  large_holes = 0;

  incCompactions(); // parent->changeStats().num_compactions++;

  if (shrink && size > min_size && last < size/2) {
    node_handle new_size = size/2;
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
::reportStats(FILE* s, const char* pad, unsigned flags) const
{
  static unsigned STORAGE =
    expert_forest::STORAGE_STATS | expert_forest::STORAGE_DETAILED;
  static unsigned HOLE_MANAGER =
    expert_forest::HOLE_MANAGER_STATS | expert_forest::HOLE_MANAGER_DETAILED;


  if (flags & expert_forest::STORAGE_STATS) {
    fprintf(s, "%s", pad);
#ifdef NODE_STORAGE_PER_LEVEL
    fprintf(s, "Level %ld ", parent->getLevelNumber(this));
#endif
    fprintf(s, "Stats for \"classic\" node storage:\n");

    // anything?
  }

  if (! (flags & HOLE_MANAGER)) return;

  fprintf(s, "%sStats for \"classic\" hole management:\n", pad);

  if (flags & expert_forest::HOLE_MANAGER_STATS) {
    unsigned long holemem = hole_slots * sizeof(node_handle);
    fprintf(s, "%s    ", pad);
    fprintmem(s, holemem, flags & expert_forest::HUMAN_READABLE_MEMORY);
    fprintf(s, " wasted in holes\n");
    unsigned long fragmem = fragment_slots * sizeof(node_handle);
    fprintf(s, "%s    ", pad);
    fprintmem(s, fragmem, flags & expert_forest::HUMAN_READABLE_MEMORY);
    fprintf(s, " wasted in fragments\n");
  }

  if (! (flags & expert_forest::HOLE_MANAGER_DETAILED)) return;

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

  fprintf(s, "%s    Hole Chains (size, count):\n", pad);
  for (std::map<int, int>::iterator iter = chainLengths.begin();
      iter != chainLengths.end(); ++iter)
    {
      if (iter->first<0)
        fprintf(s, "%s\tlarge: %d\n", pad, iter->second);
      else
        fprintf(s, "%s\t%5d: %d\n", pad, iter->first, iter->second);
    }
  fprintf(s, "%s    End of Hole Chains\n", pad);
}

void MEDDLY::old_node_storage::showNode(FILE* s, node_address addr, bool verb) const
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

MEDDLY::node_address MEDDLY::old_node_storage
::makeNode(node_handle p, const node_builder &nb, node_storage_flags opt)
{
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

void MEDDLY::old_node_storage::unlinkDownAndRecycle(node_address addr)
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
  node_handle padding;
  makeHole(addr, activeNodeActualSlots(addr, padding));
  fragment_slots += padding;
}

bool MEDDLY::old_node_storage
::areDuplicates(node_address addr, const node_builder &nb) const
{
  return areDupsTempl(addr, nb);
}

bool MEDDLY::old_node_storage
::areDuplicates(node_address addr, const node_reader &nr) const
{
  return areDupsTempl(addr, nr);
}

void MEDDLY::old_node_storage::fillReader(node_address addr, node_reader &nr) const
{
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
  return; 
}


unsigned MEDDLY::old_node_storage::hashNode(const node_header& node) const
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


int MEDDLY::old_node_storage
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
MEDDLY::old_node_storage
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

void MEDDLY::old_node_storage
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

void MEDDLY::old_node_storage
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


const void* MEDDLY::old_node_storage
::getUnhashedHeaderOf(node_address addr) const
{
  return UH(addr);
}

const void* MEDDLY::old_node_storage
::getHashedHeaderOf(node_address addr) const
{
  return HH(addr);
}



void MEDDLY::old_node_storage::updateData(node_handle* d)
{
}

int MEDDLY::old_node_storage::smallestNode() const
{
  return slotsForNode(0);
}

void MEDDLY::old_node_storage::dumpInternalInfo(FILE* s) const
{
  if (0==data) return; // nothing to display
#ifdef NODE_STORAGE_PER_LEVEL
  fprintf(s, "Level %ld: ", getParent()->getLevelNumber(this));
#endif
  fprintf(s, "Last slot used: %ld\n", long(last));
  fprintf(s, "large_holes: %ld\n", long(large_holes));
  fprintf(s, "Grid: top = %ld bottom = %ld\n", long(holes_top), long(holes_bottom));
  verify_hole_slots = 0;
}


MEDDLY::node_address 
MEDDLY::old_node_storage
::dumpInternalNode(FILE *s, node_address a, unsigned flags) const
{
  if (a<=0) return 0;
  int awidth = digits(getParent()->getLastNode());
  if (a > last) {
    fprintf(s, "%*ld : free slots\n", awidth, long(a));
    return 0;
  }
  MEDDLY_DCASSERT(data);
  bool print = (data[a] < 0) ? (flags & 0x02) : (flags & 0x01);
  if (print) fprintf(s, "%*ld : ", awidth, a);
  if (data[a]<0) { 
    // hole
    if (print) {
      fprintf(s, "[%ld", long(data[a]));
      for (int i=1; i<3; i++) {
        fprintf(s, "|%ld", long(data[a+i]));
      }
      fprintf(s, "| ... ");
    }
    verify_hole_slots += -data[a];
    a -= data[a];  
    if (print) fprintf(s, "%ld]\n", long(data[a-1]));
  } else {
    // proper node
    if (print) {
      fprintf(s, "%*ld : ", awidth, a);
      int nElements = activeNodeActualSlots(a);
      fprintf(s, "[%ld", long(data[a]));
      for (int i=1; i<nElements-1; i++) {
        fprintf(s, "|%ld", long(data[a+i]));
      }
    }
    a += activeNodeActualSlots(a);
    if (print) fprintf(s, "%ld]\n", long(data[a-1]));
  }
  return a;
}

void MEDDLY::old_node_storage::dumpInternalTail(FILE* s) const
{
  if (verify_hole_slots != hole_slots) {
    fprintf(s, "Counted hole slots: %ld\n", long(verify_hole_slots));
    MEDDLY_DCASSERT(false);
  }
}







//
//
// Private
//
//

MEDDLY::node_handle MEDDLY::old_node_storage
::makeFullNode(node_handle p, int size, const node_builder &nb)
{
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
  copyExtraHeader(addr, nb);
  return addr;
}

MEDDLY::node_handle MEDDLY::old_node_storage
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
  return addr;
}


void MEDDLY::old_node_storage
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

//
//  NOTE: we set the first slot to contain
//  the NEGATIVE of the number of slots in the hole.
//  For recycled holes, usually that's a no-op :^)
//
MEDDLY::node_handle MEDDLY::old_node_storage::getHole(int slots)
{
  if (slots > max_request) {
    max_request = slots;
    // Traverse the large hole list, and re-insert
    // all the holes, in case some are now smaller
    // than max_request.
    node_handle curr = large_holes;
    large_holes = 0;
    for (; curr; ) {
      node_handle next = holeNext(curr);
      gridInsert(curr);
      curr = next;
    }
  }

  // First, try for a hole exactly of this size
  // by traversing the index nodes in the hole grid
  node_handle chain = 0;
  node_handle curr = holes_bottom;
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
    node_handle next = holeNext(curr);
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
    node_handle newhole = curr + slots;
    node_handle newsize = -(data[curr]) - slots;
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
    node_handle new_size = MAX( size, last + slots ) * 1.5;

    resize(new_size);
  }

  // 
  // Grab node from the end
  //
  node_handle h = last + 1;
  last += slots;
  data[h] = -slots;
  return h;
}

void MEDDLY::old_node_storage::makeHole(node_handle addr, int slots)
{
#ifdef MEMORY_TRACE
  printf("Calling makeHole(%d, %d)\n", addr, slots);
#endif

  decMemUsed(slots * sizeof(node_handle));

  hole_slots += slots;
  data[addr] = data[addr+slots-1] = -slots;

  // Check for a hole to the left
#ifdef MERGE_AND_SPLIT_HOLES
  if (data[addr-1] < 0) {
    // Merge!
#ifdef MEMORY_TRACE
    printf("Left merging\n");
#endif
    node_handle lefthole = addr + data[addr-1];
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
      node_handle new_size = size/2;
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
    node_handle righthole = addr+slots;
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

void MEDDLY::old_node_storage::gridInsert(node_handle p_offset)
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
  node_handle above = holes_bottom;
  node_handle below = 0;
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
    node_handle right = holeNext(above);
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

void MEDDLY::old_node_storage::midRemove(node_handle p_offset)
{
#ifdef MEMORY_TRACE
  printf("midRemove(%d)\n", p_offset);
#endif

  // Remove a "middle" node, either in the grid
  // or in the large hole list.
  //
  MEDDLY_DCASSERT(isHoleNonIndex(p_offset));

  node_handle left = holePrev(p_offset); 
  node_handle right = holeNext(p_offset);

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

void MEDDLY::old_node_storage::indexRemove(node_handle p_offset)
{
#ifdef MEMORY_TRACE
  printf("indexRemove(%d)\n", p_offset);
#endif

  MEDDLY_DCASSERT(!isHoleNonIndex(p_offset));
  node_handle above = holeUp(p_offset); 
  node_handle below = holeDown(p_offset); 
  node_handle right = holeNext(p_offset); 

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

void MEDDLY::old_node_storage::resize(node_handle new_size)
{
  node_handle* new_data = (node_handle*) realloc(data, new_size * sizeof(node_handle));
  if (0 == new_data) throw error(error::INSUFFICIENT_MEMORY);
  if (0== data) new_data[0] = 0;
  if (new_size > size) {
    incMemAlloc((new_size - size) * sizeof(node_handle));
  } else {
    decMemAlloc((size - new_size) * sizeof(node_handle));
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


MEDDLY::node_handle 
MEDDLY::old_node_storage::allocNode(int sz, node_handle tail, bool clear)
{
  int slots = slotsForNode(sz);
  MEDDLY_DCASSERT(slots >= extraSlots + unhashedSlots + hashedSlots);

  node_handle off = getHole(slots);
  node_handle got = -data[off];
  incMemUsed(got * sizeof(node_handle));
  MEDDLY_DCASSERT(got >= slots);
  if (clear) memset(data+off, 0, slots*sizeof(node_handle));
  setCountOf(off, 1);                     // #incoming
  setNextOf(off, temp_node_value);        // mark as a temp node
  setSizeOf(off, sz);                     // size
  data[off+slots-1] = slots - got;        // negative padding
  fragment_slots += got - slots;          
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

  const node_storage* CLASSIC_STORAGE = &THE_OLD_NODE_STORAGE;
};


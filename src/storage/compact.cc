
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

#include "compact.h"
#include "bytepack.h"

#include "hm_grid.h"

// #define DEBUG_ENCODING
// #define DEBUG_COMPACTION
// #define DEBUG_SLOW
// #define MEMORY_TRACE
// #define DEEP_MEMORY_TRACE


inline void fprintRaw(FILE* s, const char* what, unsigned char* x, int bytes)
{
  if (bytes<1) return;
  fputs(what, s);
  for (int i=0; i<bytes; i++) {
    fprintf(s, " %x", x[i]);
  }
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                    compact_storage  methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::compact_storage::compact_storage(holeman* hm) : node_storage()
{
  holeManager = hm;
  data = 0;
  memchunk = 0;
}

// ******************************************************************

MEDDLY::compact_storage::~compact_storage()
{
  delete holeManager;
}

// ******************************************************************
void MEDDLY::compact_storage::collectGarbage(bool shrink)
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
    int pbytes, ibytes;
    getStyleOf(old_off, pbytes, ibytes);
    int datalen = slotsForNode(sizeOf(old_off), pbytes, ibytes) - 1;
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

// ******************************************************************
void MEDDLY::compact_storage
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
  // verifyStats();
#endif
}

// ******************************************************************
void MEDDLY::compact_storage
::showNode(FILE* s, node_address addr, bool verb) const
{
  int pbytes, ibytes;
  getStyleOf(addr, pbytes, ibytes);
  if (verb) {
    fprintf(s, " pb : %d", pbytes);
    fprintf(s, " ib : %d", ibytes);
  }
  if (sizeOf(addr) < 0) {
    //
    // Sparse node
    //
    int nnz = -sizeOf(addr);
    if (verb)  fprintf(s, " nnz : %d", nnz);
    fprintf(s, " down: (");
    unsigned char* rawd = sparseDown(addr);
    unsigned char* rawi = sparseIndex(addr);
    unsigned char* rawe = sparseEdge(addr);
    for (int i=0; i<nnz; i++) {
      if (i) fprintf(s, ", ");

      node_handle down; 
      int index;
      dataToDown(rawd, pbytes, down);
      dataToUnsigned(rawi, ibytes, index);

      fprintf(s, "%d:", index);
      if (edgeBytes) {
        fprintf(s, "<");
        getParent()->showEdgeValue(s, rawe, 0);
        fprintf(s, ", ");
      }
      if (getParent()->isTerminalNode(down)) {
        getParent()->showTerminal(s, down);
      } else {
        fprintf(s, "%ld", long(down));
      }
      if (edgeBytes) {
        fprintf(s, ">");
        rawe += edgeBytes;
      }
      rawd += pbytes;
      rawi += ibytes;
    } // for i
    fprintf(s, ")");
  } else {
    //
    // Full node
    //
    int size = sizeOf(addr);
    if (verb) fprintf(s, " size: %d", size);
    fprintf(s, " down: [");
    unsigned char* rawd = fullDown(addr);
    unsigned char* rawe = fullEdge(addr);
    for (int i=0; i<size; i++) {
      if (i) fprintf(s, "|"); 
      node_handle down;
      int index;
      dataToDown(rawd, pbytes, down);

      if (edgeBytes) {
        fprintf(s, "<");
        getParent()->showEdgeValue(s, rawe, 0);
        fprintf(s, ", ");
      } 
      if (getParent()->isTerminalNode(down)) {
        getParent()->showTerminal(s, down);
      } else {
        fprintf(s, "%ld", long(down));
      }
      if (edgeBytes) {
        fprintf(s, ">");
        rawe += edgeBytes;
      }
      rawd += pbytes;
    } // for i
    fprintf(s, "]");
  }


  // show extra header stuff
  if (unhashedBytes) {
    getParent()->showUnhashedHeader(s, UH(addr));
  }
  if (hashedBytes) {
    getParent()->showHashedHeader(s, HH(addr));
  }
}

// ******************************************************************
MEDDLY::node_address 
MEDDLY::compact_storage
::makeNode(node_handle p, const node_builder &nb, node_storage_flags opt)
{
#ifdef DEBUG_ENCODING
  printf("compact_storage making node\n        temp:  ");
  nb.show(stdout, true);
#endif
  int nnz;
  int imax = -1;
  int pbytes = 1;

  //
  // Scan node, get node stats
  //

  if (nb.isSparse()) {
    nnz = nb.getNNZs();
    for (int z=0; z<nnz; z++) {
      MEDDLY_DCASSERT(nb.d(z));
      pbytes = MAX(pbytes, bytesRequiredForDown(nb.d(z)));
      MEDDLY_DCASSERT(nb.i(z) > imax);
      imax = nb.i(z);
    } 
  } else {
    nnz = 0;
    for (int i=0; i<nb.getSize(); i++) {
      if (nb.d(i)) {
        pbytes = MAX(pbytes, bytesRequiredForDown(nb.d(i)));
        imax = i;
        nnz++;
      }
    }
  }

  node_address addr;
  //
  // Easy case - sparse nodes disabled
  //
  if (0==(forest::policies::ALLOW_SPARSE_STORAGE & opt)) {

      addr = makeFullNode(p, imax+1, pbytes, nb);
#ifdef DEBUG_ENCODING
      printf("\n        made: ");
      showNode(stdout, addr, true);
      printf("\n    internal: ");
      dumpInternalNode(stdout, addr, 0x03);
      MEDDLY_DCASSERT(areDuplicates(addr, nb));
#endif
      return addr;
  }

  //
  // Determine byte requirements for indexes
  //
  MEDDLY_DCASSERT(imax>=0);
  int ibytes = bytesRequired4(imax);

  //
  // Easy case - full nodes disabled
  //
  if (0==(forest::policies::ALLOW_FULL_STORAGE & opt)) {

      addr = makeSparseNode(p, nnz, pbytes, ibytes, nb);
#ifdef DEBUG_ENCODING
      printf("\n        made: ");
      showNode(stdout, addr, true);
      printf("\n    internal: ");
      dumpInternalNode(stdout, addr, 0x03);
      MEDDLY_DCASSERT(areDuplicates(addr, nb));
#endif
      return addr;
  }

  //
  // Full and sparse are allowed, determine which is more compact
  //
  if (slotsForNode(-nnz, pbytes, ibytes) < slotsForNode(imax+1, pbytes, ibytes))
  {
      addr = makeSparseNode(p, nnz, pbytes, ibytes, nb);
  } else {
      addr = makeFullNode(p, imax+1, pbytes, nb);
  }
#ifdef DEBUG_ENCODING
  printf("\n        made: ");
  showNode(stdout, addr, true);
  printf("\n    internal: ");
  dumpInternalNode(stdout, addr, 0x03);
  MEDDLY_DCASSERT(areDuplicates(addr, nb));
#endif
  return addr;
}

// ******************************************************************

void MEDDLY::compact_storage::unlinkDownAndRecycle(node_address addr)
{
  //
  // Unlink down pointers
  //
  int size = sizeOf(addr);
  int pbytes = pointerBytesOf(addr);
  unsigned char* down;
  if (size < 0) {
    size = -size;
    down = sparseDown(addr);
  } else {
    down = fullDown(addr);
  }
  unlinkDown(pbytes, down, size);

  //
  // Recycle
  //
  holeManager->recycleChunk(addr, activeNodeActualSlots(addr));
}

// ******************************************************************

bool MEDDLY::compact_storage
::areDuplicates(node_address addr, const node_builder &nb) const
{
  return areDupsTempl(addr, nb);
}

// ******************************************************************

bool MEDDLY::compact_storage::
areDuplicates(node_address addr, const node_reader &nr) const
{
  return areDupsTempl(addr, nr);
}

// ******************************************************************

void MEDDLY::compact_storage
::fillReader(node_address addr, node_reader &nr) const
{
#ifdef DEBUG_ENCODING
  printf("compact_storage filling reader\n    internal: ");
  dumpInternalNode(stdout, addr, 0x03);
  printf("        node: ");
  showNode(stdout, addr, true);
#endif
  // TBD: Copy hashed header

  int size = sizeOf(addr);
  if (size < 0) {
    //
    // Node is sparse
    //
    int pbytes;
    int ibytes;
    getStyleOf(addr, pbytes, ibytes);
    if (nr.isFull()) {
      readFullFromSparse(pbytes, ibytes, addr, -size, nr);
    } else {
      readSparseFromSparse(pbytes, ibytes, addr, -size, nr);
    }
  } else {
    //
    // Node is full
    //
    if (nr.isFull()) {
      readFullFromFull(pointerBytesOf(addr), addr, size, nr);
    } else {
      readSparseFromFull(pointerBytesOf(addr), addr, size, nr);
    }
  }
#ifdef DEBUG_ENCODING
  printf("\n        temp:  ");
  nr.show(stdout, getParent(), true);
  printf("\n");
  MEDDLY_DCASSERT(areDuplicates(addr, nr));
#endif
}

// ******************************************************************

unsigned MEDDLY::compact_storage::hashNode(const node_header& node) const
{
  hash_stream s;
  s.start(node.level);

  // Do the hashed header part, if any
  // TBD

  //
  // Hash the node itself
  int size = sizeOf(node.offset);
  if (size < 0) {
    hashSparse(s, node.offset, -size);
  } else {
    hashFull(s, node.offset, size);
  }

  return s.finish();
}

// ******************************************************************

int MEDDLY::compact_storage::
getSingletonIndex(node_address addr, node_handle &down) const
{
  int size = sizeOf(addr);
  if (size < 0) {
    int pbytes, ibytes;
    getStyleOf(addr, pbytes, ibytes);
    dataToDown(sparseDown(addr), pbytes, down);
    int index;
    dataToUnsigned(sparseIndex(addr), ibytes, index);
    return index;
  } else {
    return getSingletonFull(pointerBytesOf(addr), addr, size, down);
  }
}

// ******************************************************************

MEDDLY::node_handle 
MEDDLY::compact_storage::getDownPtr(node_address addr, int index) const
{
  int size = sizeOf(addr);
  if (size<0) {
    int pbytes, ibytes;
    getStyleOf(addr, pbytes, ibytes);
    int z = findSparseIndex(ibytes, addr, -size, index);
    if (z<0) return 0;
    node_handle down;
    dataToDown(sparseDown(addr) + z*pbytes, pbytes, down);
    return down;
  } else {
    if (index < size) {
      int pbytes = pointerBytesOf(addr);
      node_handle down;
      dataToDown(fullDown(addr) + index*pbytes, pbytes, down);
      return down;
    } else {
      return 0;
    }
  }
}

// ******************************************************************

void MEDDLY::compact_storage
::getDownPtr(node_address addr, int index, int& ev, node_handle& dn) const
{
  int size = sizeOf(addr);
  if (size<0) {
    int pbytes, ibytes;
    getStyleOf(addr, pbytes, ibytes);
    int z = findSparseIndex(ibytes, addr, -size, index);
    if (z<0) {
      dn = 0;
      ev = 0;
    } else {
      dataToDown(sparseDown(addr) + z*pbytes, pbytes, dn);
      ev = ((int*)(sparseEdge(addr) + z*edgeBytes))[0];
    }
  } else {
    if (index < size) {
      int pbytes = pointerBytesOf(addr);
      dataToDown(fullDown(addr) + index*pbytes, pbytes, dn);
      ev = ((int*)(sparseEdge(addr) + index*edgeBytes))[0];
    } else {
      dn = 0;
      ev = 0;
    }
  }
}

// ******************************************************************

void MEDDLY::compact_storage
::getDownPtr(node_address addr, int index, float& ev, node_handle& dn) const
{
  int size = sizeOf(addr);
  if (size<0) {
    int pbytes, ibytes;
    getStyleOf(addr, pbytes, ibytes);
    int z = findSparseIndex(ibytes, addr, -size, index);
    if (z<0) {
      dn = 0;
      ev = 0;
    } else {
      dataToDown(sparseDown(addr) + z*pbytes, pbytes, dn);
      ev = ((float*)(sparseEdge(addr) + z*edgeBytes))[0];
    }
  } else {
    if (index < size) {
      int pbytes = pointerBytesOf(addr);
      dataToDown(fullDown(addr) + index*pbytes, pbytes, dn);
      ev = ((float*)(sparseEdge(addr) + index*edgeBytes))[0];
    } else {
      dn = 0;
      ev = 0;
    }
  }
}

// ******************************************************************

const void* MEDDLY::compact_storage
::getUnhashedHeaderOf(node_address addr) const
{
  return UH(addr);
}

// ******************************************************************

const void* MEDDLY::compact_storage
::getHashedHeaderOf(node_address addr) const
{
  return HH(addr);
}

//
//
// Protected
//
//

void MEDDLY::compact_storage::localInitForForest(const expert_forest* f)
{
  edgeBytes = f->edgeBytes();
  unhashedBytes = f->unhashedHeaderBytes();
  hashedBytes = f->hashedHeaderBytes();
  MEDDLY_DCASSERT(holeManager);
  holeManager->setParent(this);
}

// ******************************************************************

void MEDDLY::compact_storage::updateData(node_handle* d)
{
  data = d;
  updateCountArray(data + count_index);
  updateNextArray(data + next_index);
  memchunk = data + mem_index;
}

// ******************************************************************

int MEDDLY::compact_storage::smallestNode() const
{
  return slotsForNode(0, 0, 0);
}

// ******************************************************************

void MEDDLY::compact_storage::dumpInternalInfo(FILE* s) const
{
#ifdef NODE_STORAGE_PER_LEVEL
  fprintf(s, "Level %ld: ", getParent()->getLevelNumber(this));
#endif
  holeManager->dumpInternalInfo(s);
}

// ******************************************************************

void MEDDLY::compact_storage::dumpInternalTail(FILE* s) const
{
  holeManager->dumpInternalTail(s);
}

// ******************************************************************

MEDDLY::node_address 
MEDDLY::compact_storage
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
      int nElements = activeNodeActualSlots(a) - 1;
      fprintf(s, "[in %ld, next %ld", long(data[a]), long(data[a+1]));
      int sz = sizeOf(a);
      fprintf(s, ", size %d", sz);
      unsigned char x = rawStyleOf(a);
      fprintf(s, ", style %c%c%c:%c%c:%c%c%c",
        (x & 0x80) ? '1' : '0',
        (x & 0x40) ? '1' : '0',
        (x & 0x20) ? '1' : '0',
        (x & 0x10) ? '1' : '0',
        (x & 0x08) ? '1' : '0',
        (x & 0x04) ? '1' : '0',
        (x & 0x02) ? '1' : '0',
        (x & 0x01) ? '1' : '0'
      );
      fprintRaw(s, ", uh", UH(a), unhashedBytes);
      fprintRaw(s, ", hh", HH(a), hashedBytes);

      unsigned char* end;
      if (sz<0) {
        // sparse
        fprintRaw(s, ", down", sparseDown(a), -sz * pointerBytesOf(a));
        fprintRaw(s, ", index", sparseIndex(a), -sz * indexBytesOf(a));
        fprintRaw(s, ", edge", sparseEdge(a), -sz * edgeBytes);
        end = sparseEdge(a) - sz * edgeBytes;
      } else {
        // full
        fprintRaw(s, ", down", fullDown(a), sz * pointerBytesOf(a));
        fprintRaw(s, ", edge", fullEdge(a), sz * edgeBytes);
        end = fullEdge(a) + sz * edgeBytes;
      }

      unsigned char* tail = (unsigned char*) (data + a + nElements);
      fprintRaw(s, ", pad", end, tail-end);
      fprintf(s, ", tail %ld]\n", long(data[a+nElements]));
    }
    a += activeNodeActualSlots(a);
  }
  return a;
}

// ******************************************************************
//
//
// Private
//
//
// ******************************************************************

MEDDLY::node_address
MEDDLY::compact_storage
::makeFullNode(node_handle p, int size, int pbytes, const node_builder &nb)
{
  int slots = slotsForNode(size, pbytes, 0);
  node_address addr = allocNode(slots, p, true);
  MEDDLY_DCASSERT(1==getCountOf(addr));

  setSizeOf(addr, size);
  setStyleOf(addr, pbytes, 1);
  copyExtraHeader(addr, nb);

  if (nb.isSparse()) {
    return copySparseIntoFull(pbytes, nb, size, addr);
  } else {
    return copyFullIntoFull(pbytes, nb, size, addr);
  }
}

// ******************************************************************

MEDDLY::node_address
MEDDLY::compact_storage::makeSparseNode(node_handle p, int size, 
  int pbytes, int ibytes, const node_builder &nb)
{
  int slots = slotsForNode(-size, pbytes, ibytes);
  node_address addr = allocNode(slots, p, true);
  MEDDLY_DCASSERT(1==getCountOf(addr));

  setSizeOf(addr, -size);
  setStyleOf(addr, pbytes, ibytes);
  copyExtraHeader(addr, nb);

  if (nb.isSparse()) {
    return copySparseIntoSparse(pbytes, ibytes, nb, size, addr);
  } else {
    return copyFullIntoSparse(pbytes, ibytes, nb, size, addr);
  }
}

// ******************************************************************

MEDDLY::node_address
MEDDLY::compact_storage::allocNode(int slots, node_handle tail, bool clear)
{
  node_handle off = holeManager->requestChunk(slots);
  node_handle got = -data[off];
  incMemUsed(got * sizeof(node_handle));
  MEDDLY_DCASSERT(got >= slots);
  if (clear) memset(data+off, 0, slots*sizeof(node_handle));
  setCountOf(off, 1);                     // #incoming
  setNextOf(off, -1);                     // mark as a temp node
  data[off+slots-1] = slots - got;        // negative padding
  data[off+got-1] = tail;                 // tail entry
#ifdef MEMORY_TRACE
  printf("Allocated new node, asked %d, got %d, position %d\n", slots, got, off);
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
// *                      compact_grid methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::compact_grid::compact_grid(holeman* hm) : compact_storage(hm)
{
}

MEDDLY::compact_grid::~compact_grid()
{
}

MEDDLY::node_storage* MEDDLY::compact_grid
::createForForest(expert_forest* f) const
{
  compact_storage* nns = new compact_grid(new hm_grid);
  nns->initForForest(f);
  return nns;
}

const char* MEDDLY::compact_grid::getStorageName() const
{
  return "compact node storage with grid for holes";
}


// ******************************************************************
// *                                                                *
// *                   front-end global variables                   *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
  compact_grid THE_COMPACT_GRID(0);

  const node_storage* COMPACT_GRID = &THE_COMPACT_GRID;
};




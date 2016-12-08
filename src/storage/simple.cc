
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

// #define DEBUG_ENCODING
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

inline int bytesForSlots(int slots)
{
  return slots * sizeof(MEDDLY::node_handle);
}

namespace MEDDLY {
  class simple_separated;

  // for now...

  class simple_storage;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                     simple_separated class                     *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Original node storage mechanism in a forest.
    Memory management is completely separated from the class,
    which makes this implementation a little different
    from the original.

    Details of node storage are left to the derived forests.
    However, every active node is stored in the following format.

      common   {  slot[0] : incoming count (MSB cleared), >= 0.
      header --{  slot[1] : next pointer in unique table or special value.
               {  slot[2] : size.  >=0 for full storage, <0 for sparse.

      unhashed    {       : slots used for any extra information
      header    --{       : as needed on a forest-by-forest basis.
      (optional)  {       : Info does NOT affect node uniqueness.

      hashed      {       : slots used for any extra information
      header    --{       : as needed on a forest-by-forest basis.
      (optional)  {       : Info DOES affect node uniqueness.

                  {       : Downward pointers.
                  {       : If full storage, there are size pointers
      down -------{       : and entry i gives downward pointer i.
                  {       : If sparse storage, there are -size pointers
                  {       : and entry i gives a pointer but the index
                  {       : corresponds to index[i], below.
          
                  {       : Index entries.
                  {       : Unused for full storage.
      index ------{       : If sparse storage, entry i gives the
      (sparse)    {       : index for outgoing edge i, and there are
                  {       : -size entries.

                  {       : Edge values.
      edge        {       : If full storage, there are size * edgeSize
      values -----{       : slots; otherwise there are -size * edgeSize
                  {       : slots.  Derived forests are responsible
                  {       : for packing information into these slots.

                { -padlen : Any node padding to allow for future
                {         : expansion, or for memory management purposes
                {         : (e.g., memory hold is larger than requested).
      padding --{         : padlen is number of padded slots.  If the
                {         : first entry after the node proper is negative,
                {         : then it specifies the number of long padding
                {         : slots; otherwise, there is no padding.

      tail    --{ slot[L] : The forest node number,
                            guaranteed to be non-negative
                            (MSB cleared).


*/
class MEDDLY::simple_separated : public node_storage {
  public:
    simple_separated(const char* n, expert_forest* f, const memory_manager_style* mst);
    virtual ~simple_separated();

  // required interface
  public:
    virtual void collectGarbage(bool shrink);
    virtual void reportStats(output &s, const char* pad, unsigned flags) const;

    virtual node_address makeNode(node_handle p, const unpacked_node &nb, 
        node_storage_flags opt);

    virtual void unlinkDownAndRecycle(node_address addr);

    virtual bool areDuplicates(node_address addr, const unpacked_node &nr) const;
    virtual void fillUnpacked(unpacked_node &nr, node_address addr, unpacked_node::storage_style) const;

    virtual unsigned hashNode(int level, node_address addr) const;
    virtual int getSingletonIndex(node_address addr, node_handle &down) const;
    virtual node_handle getDownPtr(node_address addr, int index) const;
    virtual void getDownPtr(node_address addr, int ind, int& ev, node_handle& dn) const;
    virtual void getDownPtr(node_address addr, int ind, float& ev, node_handle& dn) const;
    virtual const void* getUnhashedHeaderOf(node_address addr) const;
    virtual const void* getHashedHeaderOf(node_address addr) const;

  protected:
    virtual void updateData(node_handle* d);
    virtual int smallestNode() const;
    virtual void dumpInternalInfo(output &) const;
    virtual node_address firstNodeAddress() const;
    virtual node_address 
    dumpInternalNode(output &, node_address addr, unsigned flags) const;
    virtual void dumpInternalTail(output &) const;


  // --------------------------------------------------------
  // |  Misc. helpers.
  private:
      /// How many int slots would be required for a node with given size.
      ///   @param  sz  negative for sparse storage, otherwise full.
      inline int slotsForNode(int sz) const {
          int node_slots = (sz<0) ? (2+slots_per_edge) * (-sz) : (1+slots_per_edge) * sz;
          return header_slots + unhashed_slots + hashed_slots + node_slots;
      }

      /** Create a new node, stored as truncated full.
          Space is allocated for the node, and data is copied.
            @param  p     Node handle number.
            @param  size  Number of downward pointers.
            @param  nb    Node data is copied from here.
            @return       The "address" of the new node.
      */
      node_handle makeFullNode(node_handle p, int size, const unpacked_node &nb);

      /** Create a new node, stored sparsely.
          Space is allocated for the node, and data is copied.
            @param  p     Node handle number.
            @param  size  Number of nonzero downward pointers.
            @param  nb    Node data is copied from here.
            @return       The "address" of the new node.
      */
      node_handle makeSparseNode(node_handle p, int size, 
        const unpacked_node &nb);


      //
      // binary search for an index
      //
      inline static int findSparseIndex(int i, const node_handle* index, int N) {
        int low = 0;
        int high = N;
        while (low < high) {
          int z = (low+high)/2;
          MEDDLY_CHECK_RANGE(0, z, N);
          if (index[z] == i) return z;
          if (index[z] < i) low = z + 1;
          else              high = z;
        }
        return -1;
      }

  private:
    memory_manager* MM;

    //
    // Header indexes that are fixed
    //
    static const int count_slot = 0;
    static const int next_slot = 1;
    static const int size_slot = 2;
    static const int header_slots = size_slot+1;

    //
    // Header info that varies by forest
    //
    char unhashed_start;
    char unhashed_slots;
    char hashed_start;
    char hashed_slots;
    char down_start;
    char slots_per_edge;
};


// ******************************************************************
// *                                                                *
// *                                                                *
// *                    simple_separated methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::simple_separated
::simple_separated(const char* n, expert_forest* f, const memory_manager_style* mst)
: node_storage(n, f)
{
  MM = mst->initManager(sizeof(node_handle), slotsForNode(0));

  unhashed_start = header_slots;
  unhashed_slots = slotsForBytes(f->unhashedHeaderBytes());
  hashed_start = unhashed_start + unhashed_slots;
  hashed_slots = slotsForBytes(f->hashedHeaderBytes());
  down_start = hashed_start + hashed_slots;
  slots_per_edge = slotsForBytes(f->edgeBytes());
}

MEDDLY::simple_separated::~simple_separated()
{
  // TBD - special steps to recycle all nodes?
  delete MM;
}

void MEDDLY::simple_separated::collectGarbage(bool shrink)
{
  // TBD
}

void MEDDLY::simple_separated::reportStats(output &s, const char* pad, 
  unsigned flags) const
{
  static unsigned STORAGE = 
    expert_forest::STORAGE_STATS | expert_forest::STORAGE_DETAILED;

  if (flags & STORAGE) {
    s << pad << "Stats for " << getStyleName() << "\n";

    // anything for us?

    MM->reportStats(s, pad, flags & expert_forest::STORAGE_DETAILED);
  }


  /*
#ifdef DEVELOPMENT_CODE
  verifyStats();
#endif
  */
}

MEDDLY::node_address MEDDLY::simple_separated
::makeNode(node_handle p, const unpacked_node &nb, node_storage_flags opt)
{
#ifdef DEBUG_ENCODING
  printf("simple_separated making node\n        temp:  ");
  nb.show(stdout, true);
#endif
  const node_handle tv = getParent()->getTransparentNode();

  //
  // Easy case - sparse nodes disabled
  //
  if (0==(forest::policies::ALLOW_SPARSE_STORAGE & opt)) {
    if (nb.isSparse()) {
      int truncsize = -1;
      for (int z=0; z<nb.getNNZs(); z++) {
        truncsize = MAX(truncsize, nb.i(z));
      }
      return (truncsize<0) ? tv : makeFullNode(p, truncsize+1, nb);
    } else {
      for (int i=nb.getSize()-1; i>=0; i--) {
        if (nb.d(i)!=tv) {
          return makeFullNode(p, i+1, nb);
        }
      }
      return tv;
    }
  }

  //
  // Easy case - full nodes disabled
  //
  if (0==(forest::policies::ALLOW_FULL_STORAGE & opt)) {
    if (nb.isSparse()) {
      return makeSparseNode(p, nb.getNNZs(), nb);
    } else {
      int nnz = 0;
      for (int i=0; i<nb.getSize(); i++) {
        if (nb.d(i)!=tv) {
        	nnz++;
        }
      }
      return makeSparseNode(p, nnz, nb);
    }
  }

  //
  // Full and sparse are allowed, determine which is more compact
  //
  int truncsize = -1;
  int nnz;
  if (nb.isSparse()) {
    nnz = nb.getNNZs();
    if (0==nnz) {
      return tv;
    }
    for (int z=0; z<nnz; z++) {
      truncsize = MAX(truncsize, nb.i(z));
    }
//    if (truncsize<0) {
//      return tv;
//    }
    truncsize++;
  } else {
    nnz = 0;
    for (int i=0; i<nb.getSize(); i++) {
      if (nb.d(i)!=tv) {
        nnz++;
        truncsize = i;
      }
    }
    if (0==nnz) {
      return tv;
    }
    truncsize++;
  }

  return (slotsForNode(-nnz) < slotsForNode(truncsize)) ? makeSparseNode(p, nnz, nb) : makeFullNode(p, truncsize, nb);
}



void MEDDLY::simple_separated::unlinkDownAndRecycle(node_address addr)
{
  MEDDLY_DCASSERT(MM);
  node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  
  const int size = chunk[size_slot];

  //
  // Unlink down pointers
  //
  const int num_down = size < 0 ? -size : size;
  node_handle* down = chunk + down_start;
  for (int i=0; i<num_down; i++) {
    getParent()->unlinkNode(down[i]);
  }

  //
  // Determine number of slots in this node
  //
  size_t actual_slots = slotsForNode(size);
  if (chunk[actual_slots] < 0) {
    // padding
    actual_slots += (-chunk[actual_slots]);
  }

  //
  // Recycle
  //
  MM->recycleChunk(addr, actual_slots);
}



bool MEDDLY::simple_separated
::areDuplicates(node_address addr, const unpacked_node &n) const
{
  MEDDLY_DCASSERT(MM);
  const node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  MEDDLY_DCASSERT( n.hasEdges() == (slots_per_edge>0) );
  
  const node_handle size = chunk[size_slot];

  //
  // Compare any extra header information
  //

  if (unhashed_slots) {
    if (memcmp(chunk + unhashed_start, n.UHptr(), n.UHbytes())) return false;
  }
  if (hashed_slots) {
    if (memcmp(chunk + hashed_start, n.HHptr(), n.HHbytes())) return false;
  }

  
  //
  // Compare edges
  //

  const node_handle tv = getParent()->getTransparentNode();
  const node_handle* down = chunk + down_start;
  if (size < 0) {
    //
    // Node is stored sparsely
    //
    const int nnz = -size;
    const node_handle* index = down + nnz;
    const node_handle* edge = slots_per_edge ? (index + nnz) : 0;

    if (n.isFull()) {
      //
      // n is full
      //
      int i = 0;
      for (int z=0; z<nnz; z++) {
        if (index[z] >= n.getSize()) return false;  // too large
        // loop - skipped edges must be transparent
        for (; i<index[z]; i++) { 
          if (slots_per_edge) {
            if (!getParent()->isTransparentEdge(n.d(i), n.eptr(i))) return false;
          } else {
            if (n.d(i)!=tv) return false;
          }
        } // for
        // compare this[z] with n[i]
        if (down[z] != n.d(i)) return false;
        if (edge) if (!getParent()->areEdgeValuesEqual(
          edge + z*slots_per_edge, n.eptr(index[z])
        )) return false;
        i++;
      }
      return true;
    } else {
      //
      // n is sparse
      //
      if (n.getNNZs() != nnz) return false;
      // check that down matches
      for (int z=0; z<nnz; z++) {
        if (index[z] != n.i(z)) return false;
        if (down[z] != n.d(z))  return false;
      }
      // check that edges match
      if (n.hasEdges()) {
        for (int z=0; z<nnz; z++) {
          if (!getParent()->areEdgeValuesEqual(edge + z*slots_per_edge, n.eptr(z))) {
            return false;
          }
        }
      }
      // must be equal
      return true;
    }
  } else {
    //
    // Node is stored truncated full
    //
    const node_handle* edge = slots_per_edge ? (down + size) : 0;

    if (n.isFull()) {
      //
      // n is full
      //
      if (size > n.getSize()) return false;
      // check down
      int i;
      for (i=0; i<size; i++) {
        if (down[i] != n.d(i)) return false;
      }
      for (; i<n.getSize(); i++) {
          if (slots_per_edge) {
            if (!getParent()->isTransparentEdge(n.d(i), n.eptr(i))) return false;
          } else {
            if (n.d(i)!=tv) return false;
          }
      }
      // check edges
      if (n.hasEdges()) {
        for (int i=0; i<size; i++) if (down[i]) {
          if (!getParent()->areEdgeValuesEqual(edge + i*slots_per_edge, n.eptr(i))) {
            return false;
          }
        }
      }
      // must be equal
      return true;
    } else {
      //
      // n is sparse
      //
      int i = 0;
      // check down
      for (int z=0; z<n.getNNZs(); z++) {
        if (n.i(z) >= size) return false;
        // loop - skipped edges must be transparent
        for (; i<n.i(z); i++) {
          if (slots_per_edge) {
            if (!getParent()->isTransparentEdge(down[i], edge+i*slots_per_edge)) return false;
          } else {
            if (down[i]!=tv) return false;
          }
        }
        if (n.d(z) != down[i]) return false;
        if (n.hasEdges()) {
          if (!getParent()->areEdgeValuesEqual(edge + n.i(z) * slots_per_edge, n.eptr(z))) {
            return false;
          }
        }
        i++;
      }
      if (i<size) return false; // there WILL be a non-zero down
      // must be equal
      return true;
    }
  }
}


void MEDDLY::simple_separated
::fillUnpacked(unpacked_node &nr, node_address addr, unpacked_node::storage_style st2) const
{
#ifdef DEBUG_ENCODING
  printf("simple_separeted filling reader\n    internal: ");
  dumpInternalNode(stdout, addr, 0x03);
  printf("        node: ");
  showNode(stdout, addr, true);
#endif

  MEDDLY_DCASSERT(MM);
  const node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  MEDDLY_DCASSERT( n.hasEdges() == (slots_per_edge>0) );
  
  //
  // Copy any extra header information
  //

  if (unhashed_slots) {
    memcpy(nr.UHdata(), chunk + unhashed_start, nr.UHbytes());
  }
  if (hashed_slots) {
    memcpy(nr.HHdata(), chunk + hashed_start, nr.HHbytes());
  }
  
  //
  // Copy everything else
  //

  const node_handle size = chunk[size_slot];
  const node_handle* down = chunk + down_start;

  /*
      Set the unpacked node storage style based on settings
  */
  switch (st2) {
      case unpacked_node::FULL_NODE:     nr.bind_as_full(true);    break;
      case unpacked_node::SPARSE_NODE:   nr.bind_as_full(false);   break;
      case unpacked_node::AS_STORED:     nr.bind_as_full(size>=0); break;

      default:            assert(0);
  };


  if (size < 0) {
    //
    // Node is sparse
    //
    int nnz = -size;
    const node_handle* index = down + nnz;
    const node_handle* edge = slots_per_edge ? (index + nnz) : 0;

    if (nr.isFull()) {
      //
      // Copying into a full node
      //
      for (int i=0; i<nr.getSize(); i++) {
        if (nr.hasEdges()) {
          getParent()->getTransparentEdge(nr.d_ref(i), nr.eptr_write(i));
        } else {
          nr.d_ref(i) = getParent()->getTransparentNode();
        }
      }
      for (int z=0; z<nnz; z++) {
        int i = index[z];
        nr.d_ref(i) = down[z];
        if (nr.hasEdges()) {
          memcpy(nr.eptr_write(i), edge + z*slots_per_edge, nr.edgeBytes());
        }
      }
#ifdef DEBUG_ENCODING
      printf("\n        temp:  ");
      nr.show(stdout, getParent(), true);
      printf("\n");
#endif
      return;
    }

    //
    // Copying into a sparse node
    //

    for (int z=0; z<nnz; z++) {
      nr.d_ref(z) = down[z];
      nr.i_ref(z) = index[z];
      if (nr.hasEdges()) {
        memcpy(nr.eptr_write(z), edge + z*slots_per_edge, nr.edgeBytes());
      }
    } // for z
    nr.shrinkSparse(nnz);
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
  const node_handle* edge = slots_per_edge ? (down + size) : 0;

  if (nr.isFull()) {
    //
    // Copying into a full node
    //
    int i;
    for (i=0; i<size; i++) {
      nr.d_ref(i) = down[i];
      if (nr.hasEdges()) {
        memcpy(nr.eptr_write(i), edge + i*slots_per_edge, nr.edgeBytes());
      }
    } 
    if (unpacked_node::AS_STORED == st2) {
      nr.shrinkFull(size);
    } else {
      for (; i<nr.getSize(); i++) {
        getParent()->getTransparentEdge(nr.d_ref(i), nr.eptr_write(i));
      }
    }
#ifdef DEBUG_ENCODING
    printf("\n        temp:  ");
    nr.show(stdout, getParent(), true);
    printf("\n");
#endif
    return;
  }

  //
  // Copying into a sparse node
  //

  int z = 0;
  for (int i=0; i<size; i++) if (down[i]) {
    nr.d_ref(z) = down[i];
    nr.i_ref(z) = i;
    if (nr.hasEdges()) {
      memcpy(nr.eptr_write(z), edge + i*slots_per_edge, nr.edgeBytes());
    }
    z++;
  } // for i
  nr.shrinkSparse(z);
#ifdef DEBUG_ENCODING
  printf("\n        temp:  ");
  nr.show(stdout, getParent(), true);
  printf("\n");
#endif
}


unsigned MEDDLY::simple_separated::hashNode(int level, node_address addr) const
{
  MEDDLY_DCASSERT(MM);
  const node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  hash_stream s;
  s.start(0);

  //
  // Hash any header portion
  //

  if (hashed_slots) {
    s.push(chunk + hashed_start, getParent()->hashedHeaderBytes());
  }

  //
  // Hash the node itself
  //
  
  const int size = chunk[size_slot];
  const node_handle* down = chunk + down_start;

  if (size < 0) {
    //
    // Node is sparse
    //
    int nnz = -size;
    const node_handle* index = down + nnz;
    const node_handle* edge = slots_per_edge ? (index + nnz) : 0;
    if (getParent()->areEdgeValuesHashed()) {
      const int edge_bytes = bytesForSlots(slots_per_edge);
      for (int z=0; z<nnz; z++) {
        s.push(index[z], down[z]);
        s.push(edge + z * slots_per_edge, edge_bytes);
      } // for z
    } else {
      for (int z=0; z<nnz; z++) {
        s.push(index[z], down[z]);
      } // for z
    }
  } else {
    //
    // Node is full
    //
    const node_handle* edge = slots_per_edge ? (down + size) : 0;
    if (getParent()->areEdgeValuesHashed()) {
      const int edge_bytes = bytesForSlots(slots_per_edge);
      for (int i=0; i<size; i++) {
        if (!getParent()->isTransparentEdge(down[i], edge + i*slots_per_edge)) {
          s.push(i, down[i]);
          s.push(edge + i * slots_per_edge, edge_bytes);
    	  }
      } // for z
    } else {
      const node_handle tv=getParent()->getTransparentNode();
      for (int i=0; i<size; i++) {
    	  if (down[i]!=tv) {
          s.push(i, down[i]);
    	  }
      } // for z
    }
  }

  return s.finish();
}


int MEDDLY::simple_separated
::getSingletonIndex(node_address addr, node_handle &down) const
{
  MEDDLY_DCASSERT(MM);
  const node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  const node_handle size = chunk[size_slot];

  if (size<0) {
    //
    // sparse node --- easy
    //
    if (size != -1) return -1;
    down = chunk[0+down_start];         // 0+ stops a compiler warning
    return chunk[down_start - size];    // -size is number of nonzeroes
  } 
  
  //
  // full node
  //
  const node_handle tv=getParent()->getTransparentNode();
  const node_handle* dnptr = chunk + down_start;
  for (int i=0; i<size; i++) {
    if (tv==dnptr[i]) continue;
    if (i+1 != size) return -1;
    down = dnptr[i];
    return i;
  }
  return -1;
}


MEDDLY::node_handle 
MEDDLY::simple_separated
::getDownPtr(node_address addr, int index) const
{
  if (index<0) throw error(error::INVALID_VARIABLE);

  MEDDLY_DCASSERT(MM);
  const node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  const node_handle size = chunk[size_slot];
  const node_handle* down = chunk + down_start;

  if (size<0) {
    int z = findSparseIndex(index, down - size, -size);
    if (z<0) return getParent()->getTransparentNode();
    return down[z];
  } else {
    if (index < size) {
      return down[index];
    } else {
      return getParent()->getTransparentNode();
    }
  }
}


void MEDDLY::simple_separated
::getDownPtr(node_address addr, int index, int& ev, node_handle& dn) const
{
  if (index<0) throw error(error::INVALID_VARIABLE);

  MEDDLY_DCASSERT(MM);
  const node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(slots_per_edge>0);

  const node_handle size = chunk[size_slot];
  const node_handle* down = chunk + down_start;

  if (size<0) {
    int z = findSparseIndex(index, down - size, -size);
    if (z<0) {
      getParent()->getTransparentEdge(dn, &ev);
    } else {
      const node_handle* edge = down - 2*size;
      dn = down[z];
      ev = ((int*) (edge + z*slots_per_edge)) [0];
    }
  } else {
    if (index < size) {
      const node_handle* edge = slots_per_edge ? (down + size) : 0;
      dn = down[index];
      ev = ((int*) (edge + index*slots_per_edge)) [0];
    } else {
      getParent()->getTransparentEdge(dn, &ev);
    }
  }
}


void MEDDLY::simple_separated
::getDownPtr(node_address addr, int index, float& ev, node_handle& dn) const
{
  if (index<0) throw error(error::INVALID_VARIABLE);

  MEDDLY_DCASSERT(MM);
  const node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(slots_per_edge>0);

  const node_handle size = chunk[size_slot];
  const node_handle* down = chunk + down_start;

  if (size<0) {
    int z = findSparseIndex(index, down - size, -size);
    if (z<0) {
      getParent()->getTransparentEdge(dn, &ev);
    } else {
      const node_handle* edge = down - 2*size;
      dn = down[z];
      ev = ((float*) (edge + z*slots_per_edge)) [0];
    }
  } else {
    if (index < size) {
      const node_handle* edge = slots_per_edge ? (down + size) : 0;
      dn = down[index];
      ev = ((float*) (edge + index*slots_per_edge)) [0];
    } else {
      getParent()->getTransparentEdge(dn, &ev);
    }
  }
}



const void* MEDDLY::simple_separated
::getUnhashedHeaderOf(node_address addr) const
{
  MEDDLY_DCASSERT(MM);
  const node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(unhashed_slots);

  return chunk + unhashed_start;
}


const void* MEDDLY::simple_separated
::getHashedHeaderOf(node_address addr) const
{
  MEDDLY_DCASSERT(MM);
  const node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(hashed_slots);

  return chunk + hashed_start;
}



void MEDDLY::simple_separated::updateData(node_handle* d)
{
  MEDDLY_DCASSERT(0);

  //
  // Required for old holeman class; eventually discard?
  //
}

int MEDDLY::simple_separated::smallestNode() const
{
  return slotsForNode(0);
}

void MEDDLY::simple_separated::dumpInternalInfo(output &s) const
{
  s << "simple_separated::dumpInternalInfo\n";
  // MM->dumpInternalInfo(s);
}

void MEDDLY::simple_separated::dumpInternalTail(output &s) const
{
  s << "simple_separated::dumpInternalTail\n";
  // MM->dumpInternalTail(s);
}


MEDDLY::node_address 
MEDDLY::simple_separated::firstNodeAddress() const
{
  return MM->getFirstAddress();
}

MEDDLY::node_address 
MEDDLY::simple_separated
::dumpInternalNode(output &s, node_address a, unsigned flags) const
{
  if (a<=0) return 0;

  MEDDLY_DCASSERT(MM);

  int awidth = digits(getParent()->getLastNode());

  if (!MM->isAddressInUse(a)) {
    //
    // Hole, let memory manager deal with it
    //
    if (flags & 0x02) {
      s.put(long(a), awidth);
      s << " : ";
      MM->dumpInternalUnused(s, a);
    }
    return MM->getNextAddress(a);
  }

  //
  // proper node
  //
  node_handle* end = (node_handle*) MM->getChunkAddress(a);
  const node_handle* chunk = end;
  MEDDLY_DCASSERT(chunk);
  const bool show_node = flags & 0x01;

  if (show_node) {
    s.put(long(a), awidth);
    s << " : [";
  }
  const node_handle size = chunk[size_slot];

  //
  // common header
  //

  if (show_node) {
    s.put(long(chunk[0]));
    for (int i=1; i<header_slots; i++) {
      s.put(", ");
      s.put(long(chunk[i]));
    }
    s.put(';');
  }

  //
  // unhashed header if any
  //

  if (show_node && unhashed_slots) {
    s.put(" uh ");
    s.put_hex(long(chunk[0+unhashed_start]));
    for (int i=1; i<unhashed_slots; i++) {
      s.put(", ");
      s.put_hex(long(chunk[unhashed_start+i]));
    }
    s.put(';');
  }

  //
  // hashed header if any
  //
  
  if (show_node && hashed_slots) {
    s.put(" hh ");
    s.put_hex(long(chunk[0+hashed_start]));
    for (int i=1; i<hashed_slots; i++) {
      s.put(", ");
      s.put_hex(long(chunk[hashed_start+i]));
    }
    s.put(';');
  }

  //
  // Down pointers
  //
  end += down_start;
  const node_handle* down = end;
  const int dnlen = (size < 0) ? -size : size;

  if (show_node) {
    s.put(" down ");
    s.put(long(down[0]));
    for (int i=1; i<dnlen; i++) {
      s.put(", ");
      s.put(long(down[i]));
    }
    s.put(';');
  }
  end += dnlen;

  //
  // Indexes, if any
  //
  if (size<0) {
    if (show_node) {
      s.put(" index ");
      const node_handle* index = down - size;
      s.put(long(index[0]));
      for (int i=1; i<dnlen; i++) {
        s.put(", ");
        s.put(long(index[i]));
      }
      s.put(';');
    }
    end += dnlen;
  }

  //
  // Edge values, if any
  //
  if (slots_per_edge) {
    if (show_node) {
      s.put(" ev ");
      for (int i=0; i<dnlen; i++) {
      if (i) s.put(", ");
        for (int j=0; j<slots_per_edge; j++) {
          if (j) s.put('|');
          s.put_hex(*end);
          end++;
        } // for each slot
      } // for each edge value
      s.put(';');
    } else {
      end += slots_per_edge * dnlen;
    }
  }

  //
  // Is there any padding?
  //
  if (*end < 0) {
    if (show_node)  s << " padded with " << -(*end) << " slots;";
    end -= *end;
  }

  //
  // Last slot
  //
  if (show_node) {
    s.put(" tail ");
    s << *end << "]\n";
  }

  a += (end - chunk);
  return a;
}





//
// Helpers
//

MEDDLY::node_handle MEDDLY::simple_separated
::makeFullNode(node_handle p, int size, const unpacked_node &nb)
{
  //
  // Determine amount of memory we need, and request it
  //

  int slots_req = slotsForNode(size);
  MEDDLY_DCASSERT(slots_req > 0);
  size_t slots_given = slots_req;
  MEDDLY_DCASSERT(MM);
  unsigned long addr = MM->requestChunk(slots_given);
  if (0==addr) {
    throw error(error::INSUFFICIENT_MEMORY);
  }

  node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  //
  // Set size, incoming count, etc.
  //
  chunk[count_slot] = 1;
  chunk[size_slot] = size;

  //
  // Copy any extra header information
  //

  if (unhashed_slots) {
    memcpy(chunk + unhashed_start, nb.UHptr(), nb.UHbytes());
  }
  if (hashed_slots) {
    memcpy(chunk + hashed_start, nb.HHptr(), nb.HHbytes());
  }

  //
  // Copy downward pointers and edge values (if any)
  //

  node_handle* down = chunk + down_start;
  if (slots_per_edge) {
      //
      // There's edge values
      //
      MEDDLY_DCASSERT(nb.hasEdges());
      char* edge = (char*) (down + size); 
      int edge_bytes = bytesForSlots(slots_per_edge);
      if (nb.isSparse()) {
        for (int i=0; i<size; i++) {
          getParent()->getTransparentEdge(down[i], edge + i * edge_bytes);
        }
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
      //
      // No edge values
      //
      MEDDLY_DCASSERT(!nb.hasEdges());
      if (nb.isSparse()) {
        const node_handle tv = getParent()->getTransparentNode();
        for (int i=0; i<size; i++) {
          down[i] = tv;
        }
        for (int z=0; z<nb.getNNZs(); z++) {
          MEDDLY_CHECK_RANGE(0, nb.i(z), size);
          down[nb.i(z)] = nb.d(z);
        }
      } else {
        for (int i=0; i<size; i++) down[i] = nb.d(i);
      }
  }

  //
  // Deal with any padding and the tail
  //

  node_handle* tail = down + size + slots_per_edge * size;
  int delta = slots_given - slots_req;
  MEDDLY_DCASSERT(delta>=0);
  tail[0] = -delta;
  tail[delta] = p;  // if delta = 0, we just overwrote but that's fine

#ifdef DEBUG_ENCODING
  printf("\n        made: ");
  showNode(stdout, addr, true);
  printf("\n    internal: ");
  dumpInternalNode(stdout, addr, 0x03);
#endif
  return addr;
}

MEDDLY::node_handle MEDDLY::simple_separated
::makeSparseNode(node_handle p, int size, const unpacked_node &nb)
{
  //
  // Determine amount of memory we need, and request it
  //

  int slots_req = slotsForNode(-size);
  MEDDLY_DCASSERT(slots_req > 0);
  size_t slots_given = slots_req;
  unsigned long addr = MM->requestChunk(slots_given);
  if (0==addr) {
    throw error(error::INSUFFICIENT_MEMORY);
  }

  node_handle* chunk = (node_handle*) MM->getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  //
  // Set size, incoming count, etc.
  //
  chunk[count_slot] = 1;
  chunk[size_slot] = -size;

  //
  // Copy any extra header information
  //

  if (unhashed_slots) {
    memcpy(chunk + unhashed_start, nb.UHptr(), nb.UHbytes());
  }
  if (hashed_slots) {
    memcpy(chunk + hashed_start, nb.HHptr(), nb.HHbytes());
  }

  //
  // Copy downward pointers, indexes, and edge values (if any)
  //

  node_handle* down = chunk + down_start;
  node_handle* index = down + size;

  if (slots_per_edge) {
      //
      // There's edge values
      //
      MEDDLY_DCASSERT(nb.hasEdges());
      char* edge = (char*) (index + size); 
      int edge_bytes = bytesForSlots(slots_per_edge);
      if (nb.isSparse()) {
        for (int z=0; z<size; z++) {
          down[z] = nb.d(z);
          index[z] = nb.i(z);
        }
        // kinda hacky
        memcpy(edge, nb.eptr(0), size * edge_bytes);
      } else {
        int z = 0;
        for (int i=0; i<nb.getSize(); i++) {
          if (getParent()->isTransparentEdge(nb.d(i), nb.eptr(i))) continue;
          MEDDLY_CHECK_RANGE(0, z, size);
          down[z] = nb.d(i);
          index[z] = i;
          memcpy(edge + z * edge_bytes, nb.eptr(i), edge_bytes);
          z++;
        }
      }
  } else {
      //
      // No edge values
      //
      MEDDLY_DCASSERT(!nb.hasEdges());
      if (nb.isSparse()) {
        for (int z=0; z<size; z++) {
          down[z] = nb.d(z);
          index[z] = nb.i(z);
        }
      } else {
        int z = 0;
        const node_handle tv = getParent()->getTransparentNode();
        for (int i=0; i<nb.getSize(); i++) {
          if (nb.d(i)!=tv) {
            MEDDLY_CHECK_RANGE(0, z, size);
            down[z] = nb.d(i);
            index[z] = i;
            z++;
          }
        }
      }
  }

  //
  // Deal with any padding and the tail
  //

  node_handle* tail = down + size + slots_per_edge * size;
  int delta = slots_given - slots_req;
  MEDDLY_DCASSERT(delta>=0);
  tail[0] = -delta;
  tail[delta] = p;  // if delta = 0, we just overwrote but that's fine

#ifdef DEBUG_ENCODING
  printf("\n        made: ");
  showNode(stdout, addr, true);
  printf("\n    internal: ");
  dumpInternalNode(stdout, addr, 0x03);
#endif
  return addr;
}



// ******************************************************************
// *                                                                *
// *                                                                *
// *                 simple_separated_style methods                 *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_separated_style::simple_separated_style(const char* n)
 : node_storage_style(n)
{
}

MEDDLY::simple_separated_style::~simple_separated_style()
{
}

MEDDLY::node_storage* MEDDLY::simple_separated_style
::createForForest(expert_forest* f, const memory_manager_style* mst) const
{
//  return 0;
  return new simple_separated(getName(), f, mst);
}


//
//
//
//
//
// OLD node storage classes below
//
//
//
//
//
//
// TBD - very long term - remove from here
// 


// ******************************************************************
// *                                                                *
// *                                                                *
// *                      simple_storage class                      *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Original node storage mechanism in a forest.
    The hole manager is separated from the class,
    which makes this implementation a little different
    from the original.

    Limits: offsets must be the same as node_handles, 
    because the data structure uses those for hole data.

    Details of node storage are left to the derived forests.
    However, every active node is stored in the following format.

      common   {  slot[0] : incoming count, >= 0.
      header --{  slot[1] : next pointer in unique table or special value.
               {  slot[2] : size.  >=0 for full storage, <0 for sparse.

      unhashed    {       : slots used for any extra information
      header    --{       : as needed on a forest-by-forest basis.
      (optional)  {       : Info does NOT affect node uniqueness.

      hashed      {       : slots used for any extra information
      header    --{       : as needed on a forest-by-forest basis.
      (optional)  {       : Info DOES affect node uniqueness.

                  {       : Downward pointers.
                  {       : If full storage, there are size pointers
      down -------{       : and entry i gives downward pointer i.
                  {       : If sparse storage, there are -size pointers
                  {       : and entry i gives a pointer but the index
                  {       : corresponds to index[i], below.
          
                  {       : Index entries.
                  {       : Unused for full storage.
      index ------{       : If sparse storage, entry i gives the
      (sparse)    {       : index for outgoing edge i, and there are
                  {       : -size entries.

                  {       : Edge values.
      edge        {       : If full storage, there are size * edgeSize
      values -----{       : slots; otherwise there are -size * edgeSize
                  {       : slots.  Derived forests are responsible
                  {       : for packing information into these slots.

                { -padlen : Any node padding to allow for future
                {         : expansion, or for memory management purposes
                {         : (e.g., memory hold is larger than requested).
      padding --{         : padlen is number of padded slots.  If the
                {         : first entry after the node proper is negative,
                {         : then it specifies the number of long padding
                {         : slots; otherwise, there is no padding.

      tail    --{ slot[L] : The forest node number,
                            guaranteed to be non-negative.


    When nodes are deleted, the memory slots are marked as a "hole",
    using the following format.

          slot[0] : -numslots, the number of slots in the hole
            .
            .
            .
          slot[L] : -numslots, with L = numslots-1
        
    The first few slots of the hole are used for a hole management
    data structure, described below.
    Note that a hole is guaranteed to be at least 5 slots long
    (assuming a node of size 0, with no extra header info, is impossible).

    Hole management details are handled by another class.

*/
class MEDDLY::simple_storage : public node_storage {
  // required interface
  public:
    simple_storage(const char* n, expert_forest* f, holeman* hm);
    virtual ~simple_storage();

    virtual void collectGarbage(bool shrink);
    virtual void reportStats(output &s, const char* pad, unsigned flags) const;

    virtual node_address makeNode(node_handle p, const unpacked_node &nb, 
        node_storage_flags opt);

    virtual void unlinkDownAndRecycle(node_address addr);

    virtual bool areDuplicates(node_address addr, const unpacked_node &nr) const;
    virtual void fillUnpacked(unpacked_node &nr, node_address addr, unpacked_node::storage_style) const;
    virtual unsigned hashNode(int level, node_address addr) const;
    virtual int getSingletonIndex(node_address addr, node_handle &down) const;
    virtual node_handle getDownPtr(node_address addr, int index) const;
    virtual void getDownPtr(node_address addr, int ind, int& ev, node_handle& dn) const;
    virtual void getDownPtr(node_address addr, int ind, float& ev, node_handle& dn) const;
    virtual const void* getUnhashedHeaderOf(node_address addr) const;
    virtual const void* getHashedHeaderOf(node_address addr) const;

  protected:
    virtual void updateData(node_handle* d);
    virtual int smallestNode() const;
    virtual void dumpInternalInfo(output &) const;
    virtual node_address firstNodeAddress() const;
    virtual node_address 
    dumpInternalNode(output &, node_address addr, unsigned flags) const;
    virtual void dumpInternalTail(output &) const;

  /*
  private:
      // For debugging/display purposes
      const char* storageName;
      */

  private:
      static const long temp_node_value = -5;

      // header indexes (relative to chunk start)
      static const int count_index = 0;
      static const int next_index = 1;    
      static const int size_index = 2;
      static const int headerSlots = size_index+1;

      // Counts for extra slots
      static const int tailSlots = 1;
      static const int extraSlots = headerSlots + tailSlots;

  private:
      /// copy of the data array
      node_handle* data;
      holeman*  holeManager;

  // header sizes; vary by forest.
  private:
      /// Number of slots required for each outgoing edge's value (can be 0).
      char edgeSlots;
      /// Number of slots for extra unhashed data (typically 0).
      char unhashedSlots;
      /// Number of slots for extra hashed data (typically 0).
      char hashedSlots;


  // --------------------------------------------------------
  // |  helpers for inlines.
  private:
      inline node_handle* chunkOf(node_handle addr) const {
        MEDDLY_DCASSERT(data);
        MEDDLY_DCASSERT(holeManager);
        MEDDLY_CHECK_RANGE(1, addr, holeManager->lastSlot()+1);
        MEDDLY_DCASSERT(data[addr] >= 0);  // it's not a hole
        return data + addr;
      }
      /*
      inline node_handle* holeOf(node_handle addr) const {
        MEDDLY_DCASSERT(data);
        MEDDLY_DCASSERT(holeManager);
        MEDDLY_CHECK_RANGE(1, addr, holeManager->lastSlot()+1);
        MEDDLY_DCASSERT(data[addr] < 0);  // it's a hole
        return data + addr;
      }
      */
      inline node_handle& rawSizeOf(node_handle addr) const {
        return chunkOf(addr)[size_index];
      }
      inline node_handle  sizeOf(node_handle addr) const { 
        return rawSizeOf(addr); 
      }
      inline void setSizeOf(node_handle addr, node_handle sz) { 
        rawSizeOf(addr) = sz; 
      }

      inline node_handle* UH(node_handle addr) const {
          return chunkOf(addr) + headerSlots;
      }
      inline node_handle* HH(node_handle addr) const {
          return UH(addr) + unhashedSlots;
      }
      // full down, as a pointer
      inline node_handle* FD(node_handle addr) const {
          MEDDLY_DCASSERT(rawSizeOf(addr)>0);
          return HH(addr) + hashedSlots;
      }
      // full edges, as a pointer
      inline node_handle* FE(node_handle addr) const {
          return FD(addr) + rawSizeOf(addr);
      }
      // a particular full edge, as a pointer
      inline void* FEP(node_handle addr, int p) const {
        return FE(addr) + p * edgeSlots;
      }
      // sparse down, as a pointer
      inline node_handle* SD(node_handle addr) const {
          MEDDLY_DCASSERT(rawSizeOf(addr)<0);
          return HH(addr) + hashedSlots;
      }
      // sparse indexes, as a pointer
      inline node_handle* SI(node_handle addr) const {
          return SD(addr) - rawSizeOf(addr);  // sparse size is negative
      }
      // sparse edges, as a pointer
      inline node_handle* SE(node_handle addr) const {
          return SD(addr) - 2*rawSizeOf(addr);  // sparse size is negative
      }
      // a particular sparse edge, as a pointer
      inline void* SEP(node_handle addr, int p) const {
        return SE(addr) + p * edgeSlots;
      }
      // binary search for an index
      inline int findSparseIndex(node_handle addr, int i) const {
        int low = 0;
        int nnz = -rawSizeOf(addr);
        MEDDLY_DCASSERT(nnz>=0);
        int high = nnz;
        node_handle* index = SI(addr);
        while (low < high) {
          int z = (low+high)/2;
          MEDDLY_CHECK_RANGE(0, z, nnz);
          if (index[z] == i) return z;
          if (index[z] < i) low = z + 1;
          else              high = z;
        }
        return -1;
      }


  // --------------------------------------------------------
  // |  inlines.
  private:
      /// How many int slots would be required for a node with given size.
      ///   @param  sz  negative for sparse storage, otherwise full.
      inline int slotsForNode(int sz) const {
          int nodeSlots = (sz<0) ? (2+edgeSlots) * (-sz) : (1+edgeSlots) * sz;
          return extraSlots + unhashedSlots + hashedSlots + nodeSlots;
      }

      /// Find actual number of slots used for this active node.
      inline int activeNodeActualSlots(node_handle addr) const {
          int end = addr + slotsForNode(sizeOf(addr))-1;
          // account for any padding
          if (data[end] < 0) {
            end -= data[end];
          }
          return end - addr + 1;
      }



  // --------------------------------------------------------
  // |  Misc. helpers.
  private:

      /** Create a new node, stored as truncated full.
          Space is allocated for the node, and data is copied.
            @param  p     Node handle number.
            @param  size  Number of downward pointers.
            @param  nb    Node data is copied from here.
            @return       The "address" of the new node.
      */
      node_handle makeFullNode(node_handle p, int size, const unpacked_node &nb);

      /** Create a new node, stored sparsely.
          Space is allocated for the node, and data is copied.
            @param  p     Node handle number.
            @param  size  Number of nonzero downward pointers.
            @param  nb    Node data is copied from here.
            @return       The "address" of the new node.
      */
      node_handle makeSparseNode(node_handle p, int size, 
        const unpacked_node &nb);

      void copyExtraHeader(node_address addr, const unpacked_node &nb);

      /** Allocate enough slots to store a node with given size.
          Also, stores the node size in the node.
            @param  sz      negative for sparse storage, otherwise full.
            @param  tail    Node id
            @param  clear   Should the node be zeroed.
            @return         Offset in the data array.
      */
      node_handle allocNode(int sz, node_handle tail, bool clear);


#ifdef DEVELOPMENT_CODE
      void verifyStats() const;
#endif
}; 


// ******************************************************************
// *                                                                *
// *                                                                *
// *                     simple_storage methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_storage::simple_storage(const char* n, expert_forest* f, holeman* hm)
 : node_storage(n, f)
{
  holeManager = hm;
  // storageName = sN;
  data = 0;

  edgeSlots = slotsForBytes(f->edgeBytes());
  unhashedSlots = slotsForBytes(f->unhashedHeaderBytes());
  hashedSlots = slotsForBytes(f->hashedHeaderBytes());
  MEDDLY_DCASSERT(holeManager);
  holeManager->setParent(this);
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
      memmove(curr_ptr, node_ptr, bytesForSlots(datalen));
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
::reportStats(output &s, const char* pad, unsigned flags) const
{
  static unsigned STORAGE = 
    expert_forest::STORAGE_STATS | expert_forest::STORAGE_DETAILED;

  if (flags & STORAGE) {
    s << pad << "Stats for " << getStyleName() << "\n";

    // anything for us?
  }

  holeManager->reportStats(s, pad, flags);

#ifdef DEVELOPMENT_CODE
  verifyStats();
#endif
}



MEDDLY::node_address MEDDLY::simple_storage
::makeNode(node_handle p, const unpacked_node &nb, node_storage_flags opt)
{
#ifdef DEBUG_ENCODING
  printf("simple_storage making node\n        temp:  ");
  nb.show(stdout, true);
#endif
  node_handle tv = getParent()->getTransparentNode();

  //
  // Easy case - sparse nodes disabled
  //
  if (0==(forest::policies::ALLOW_SPARSE_STORAGE & opt)) {
    if (nb.isSparse()) {
      int truncsize = -1;
      for (int z=0; z<nb.getNNZs(); z++) {
        truncsize = MAX(truncsize, nb.i(z));
      }
      return (truncsize<0) ? tv : makeFullNode(p, truncsize+1, nb);
    } else {
      for (int i=nb.getSize()-1; i>=0; i--) {
        if (nb.d(i)!=tv) {
          return makeFullNode(p, i+1, nb);
        }
      }
      return tv;
    }
  }

  //
  // Easy case - full nodes disabled
  //
  if (0==(forest::policies::ALLOW_FULL_STORAGE & opt)) {
    if (nb.isSparse()) {
      return makeSparseNode(p, nb.getNNZs(), nb);
    } else {
      int nnz = 0;
      for (int i=0; i<nb.getSize(); i++) {
        if (nb.d(i)!=tv) {
        	nnz++;
        }
      }
      return makeSparseNode(p, nnz, nb);
    }
  }

  //
  // Full and sparse are allowed, determine which is more compact
  //
  int truncsize = -1;
  int nnz;
  if (nb.isSparse()) {
    nnz = nb.getNNZs();
    if(nnz==0) {
      return tv;
    }
    for (int z=0; z<nnz; z++) {
      truncsize = MAX(truncsize, nb.i(z));
    }
//    if (truncsize<0) {
//      return tv;
//    }
    truncsize++;
  } else {
    nnz = 0;
    for (int i=0; i<nb.getSize(); i++) {
      if (nb.d(i)!=tv) {
        nnz++;
        truncsize = i;
      }
    }
    if(nnz==0) {
      return tv;
    }
    truncsize++;
  }

  return (slotsForNode(-nnz) < slotsForNode(truncsize)) ? makeSparseNode(p, nnz, nb) : makeFullNode(p, truncsize, nb);
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
::areDuplicates(node_address addr, const unpacked_node &n) const
{
  if (n.HHbytes()) {
    if (memcmp(HH(addr), n.HHptr(), n.HHbytes())) return false;
  }

  if (n.UHbytes()) {
    if (memcmp(UH(addr), n.UHptr(), n.UHbytes())) return false;
  }

  node_handle tv=getParent()->getTransparentNode();
  int size = sizeOf(addr);
  if (size<0) {
    //
    // Node is sparse
    //
    int nnz = -size;
    const node_handle* down = SD(addr);
    const node_handle* index = SI(addr);
    if (n.isFull()) {
      // check that down matches
      int i = 0;
      for (int z=0; z<nnz; z++) {
        if (index[z] >= n.getSize()) return false;
        for (; i<index[z]; i++) {
          if (n.d(i)!=tv) {
            return false;
          }
        }
        if (n.d(i) != down[z]) return false;
        i++;
      }
      for (; i<n.getSize(); i++) {
        if (n.d(i)!=tv) {
          return false;
        }
      }
      // check that edges match
      if (n.hasEdges()) {
        for (int z=0; z<nnz; z++) {
          if (!getParent()->areEdgeValuesEqual( 
                  SEP(addr, z),   n.eptr(index[z])
              )) return false;
        } // for z
      }
      // must be equal
      return true;
    }
    else {
      // n is sparse
      if (n.getNNZs() != nnz) return false;
      // check that down matches
      for (int z=0; z<nnz; z++) {
        if (index[z] != n.i(z)) return false;
        if (down[z] != n.d(z))  return false;
      }
      // check that edges match
      if (n.hasEdges()) {
        for (int z=0; z<nnz; z++) {
          if (!getParent()->areEdgeValuesEqual(SEP(addr, z), n.eptr(z))) {
            return false;
          }
        }
      }
      // must be equal
      return true;
    }
  }
  else {
    //
    // Node is truncated full
    //
    const node_handle* down = FD(addr);
    if (n.isFull()) {
      if (size > n.getSize()) return false;
      // check down
      int i;
      for (i=0; i<size; i++) {
        if (down[i] != n.d(i)) return false;
      }
      for (; i<n.getSize(); i++) {
        if (n.d(i)!=tv) {
          return false;
        }
      }
      // check edges
      if (n.hasEdges()) {
        for (int i=0; i<size; i++) if (down[i]) {
          if (!getParent()->areEdgeValuesEqual(FEP(addr, i), n.eptr(i))) {
            return false;
          }
        }
      }
      // must be equal
      return true;
    }
    else {
      // n is sparse
      int i = 0;
      // check down
      for (int z=0; z<n.getNNZs(); z++) {
        if (n.i(z) >= size) return false;
        for (; i<n.i(z); i++) {
          if (down[i]!=tv) {
            return false;
          }
        }
        if (n.d(z) != down[i]) return false;
        i++;
      }
      if (i<size) return false; // there WILL be a non-zero down
      // check edges
      if (n.hasEdges()) {
        for (int z=0; z<n.getNNZs(); z++) {
          if (!getParent()->areEdgeValuesEqual(FEP(addr, n.i(z)),  n.eptr(z))) {
            return false;
          }
        } // for z
      }
      // must be equal
      return true;
    }
  }
}

void MEDDLY::simple_storage
::fillUnpacked(unpacked_node &nr, node_address addr, unpacked_node::storage_style st2) const
{
#ifdef DEBUG_ENCODING
  printf("simple_storage filling reader\n    internal: ");
  dumpInternalNode(stdout, addr, 0x03);
  printf("        node: ");
  showNode(stdout, addr, true);
#endif

  const node_handle tv = getParent()->getTransparentNode();
  
  // Copy extra header

  if (hashedSlots) {
    memcpy(nr.HHdata(), HH(addr), getParent()->hashedHeaderBytes());
  }

  if (unhashedSlots) {
    memcpy(nr.UHdata(), UH(addr), getParent()->unhashedHeaderBytes());
  }

  // Copy everything else

  int size = sizeOf(addr);

  /*
      Set the unpacked node storage style based on settings
  */
  switch (st2) {
      case unpacked_node::FULL_NODE:     nr.bind_as_full(true);    break;
      case unpacked_node::SPARSE_NODE:   nr.bind_as_full(false);   break;
      case unpacked_node::AS_STORED:     nr.bind_as_full(size>=0); break;

      default:            assert(0);
  };

  if (size < 0) {
    //
    // Node is sparse
    //
    int nnz = -size;
    const node_handle* down = SD(addr);
    const node_handle* index = SI(addr);

    if (nr.isFull()) {
      for (int i=0; i<nr.getSize(); i++) {
        nr.d_ref(i) = tv;
        if (nr.hasEdges()) {
          memset(nr.eptr_write(i), 0, nr.edgeBytes());
        }
      }
      for (int z=0; z<nnz; z++) {
        int i = index[z];
        nr.d_ref(i) = down[z];
        if (nr.hasEdges()) {
          memcpy(nr.eptr_write(i), SEP(addr, z), nr.edgeBytes());
        }
      }
#ifdef DEBUG_ENCODING
      printf("\n        temp:  ");
      nr.show(stdout, getParent(), true);
      printf("\n");
#endif
      return;
    }

    // nr is sparse

    for (int z=0; z<nnz; z++) {
      nr.d_ref(z) = down[z];
      nr.i_ref(z) = index[z];
      if (nr.hasEdges()) {
        memcpy(nr.eptr_write(z), SEP(addr, z), nr.edgeBytes());
      }
    } // for z
    nr.shrinkSparse(nnz);
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
  const node_handle* down = FD(addr);
  if (nr.isFull()) {
    int i;
    for (i=0; i<size; i++) {
      nr.d_ref(i) = down[i];
      if (nr.hasEdges()) {
        memcpy(nr.eptr_write(i), FEP(addr, i), nr.edgeBytes());
      }
    } 
    if (unpacked_node::AS_STORED == st2) {
      nr.shrinkFull(size);
    } else {
      for (; i<nr.getSize(); i++) {
        nr.d_ref(i) = tv;
        if (nr.hasEdges()) {
          memset(nr.eptr_write(i), 0, nr.edgeBytes());
        }
      }
    }
#ifdef DEBUG_ENCODING
    printf("\n        temp:  ");
    nr.show(stdout, getParent(), true);
    printf("\n");
#endif
    return;
  }

  // nr is sparse
  int z = 0;
  for (int i=0; i<size; i++) if (down[i]) {
    nr.d_ref(z) = down[i];
    nr.i_ref(z) = i;
    if (nr.hasEdges()) {
      memcpy(nr.eptr_write(z), FEP(addr, i), nr.edgeBytes());
    }
    z++;
  } // for i
  nr.shrinkSparse(z);
#ifdef DEBUG_ENCODING
  printf("\n        temp:  ");
  nr.show(stdout, getParent(), true);
  printf("\n");
#endif
}


unsigned MEDDLY::simple_storage::hashNode(int level, node_address addr) const
{
  hash_stream s;
  s.start(0);

  // Do the hashed header part, if any

  if (hashedSlots) {
    s.push(HH(addr), getParent()->hashedHeaderBytes());
  }

  //
  // Hash the node itself
  
  int size = sizeOf(addr);
  if (size < 0) {
    // Node is sparse
    int nnz = -size;
    node_handle* down = SD(addr);
    node_handle* index = SI(addr);
    if (getParent()->areEdgeValuesHashed()) {
      for (int z=0; z<nnz; z++) {
        s.push(index[z], down[z], ((int*)SEP(addr, z))[0]);
      } // for z
    } else {
      for (int z=0; z<nnz; z++) {
        s.push(index[z], down[z]);
      } // for z
    }
  } else {
    // Node is full
    node_handle* down = FD(addr);
    node_handle tv=getParent()->getTransparentNode();
    if (getParent()->areEdgeValuesHashed()) {
      for (int i=0; i<size; i++) {
    	if (down[i]!=tv) {
          s.push(i, down[i], ((int*)FEP(addr, i))[0]);
    	}
      } // for z
    } else {
      for (int i=0; i<size; i++) {
    	if (down[i]!=tv) {
          s.push(i, down[i]);
    	}
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
  if (index<0) throw error(error::INVALID_VARIABLE);
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
  if (index<0) throw error(error::INVALID_VARIABLE);
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
  if (index<0) throw error(error::INVALID_VARIABLE);
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

void MEDDLY::simple_storage::dumpInternalInfo(output &s) const
{
  holeManager->dumpInternalInfo(s);
}

void MEDDLY::simple_storage::dumpInternalTail(output &s) const
{
  holeManager->dumpInternalTail(s);
}


MEDDLY::node_address 
MEDDLY::simple_storage::firstNodeAddress() const
{
  return 1;
}

MEDDLY::node_address 
MEDDLY::simple_storage
::dumpInternalNode(output &s, node_address a, unsigned flags) const
{
  if (a<=0) return 0;
  int awidth = digits(getParent()->getLastNode());
  if (a > holeManager->lastSlot()) {
    s.put(long(a), awidth);
    s << " : free slots\n";
    return 0;
  }
  MEDDLY_DCASSERT(data);
  if (data[a]<0) { 
    // hole
    if (flags & 0x02) {
      s.put(long(a), awidth);
      s << " : ";
      holeManager->dumpHole(s, a);
    }
    a = holeManager->chunkAfterHole(a);
  } else {
    // proper node
    if (flags & 0x01) {
      s.put(long(a), awidth);
      s << " : ";
      int nElements = activeNodeActualSlots(a);
      s << "[" << long(data[a]);
      for (int i=1; i<nElements; i++) {
        s.put('|');
        s.put(long(data[a+i]));
      }
      s << "]\n";
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
::makeFullNode(node_handle p, int size, const unpacked_node &nb)
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
::makeSparseNode(node_handle p, int size, const unpacked_node &nb)
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
        node_handle tv = getParent()->getTransparentNode();
        for (int i=0; i<nb.getSize(); i++) {
          if (nb.d(i)!=tv) {
            MEDDLY_CHECK_RANGE(0, z, size);
            down[z] = nb.d(i);
            index[z] = i;
            z++;
          }
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
::copyExtraHeader(node_address addr, const unpacked_node &nb)
{
  // copy extra header info, if any
  if (unhashedSlots) {
    memcpy(UH(addr), nb.UHptr(), nb.UHbytes());
  }
  if (hashedSlots) {
    memcpy(HH(addr), nb.HHptr(), nb.HHbytes());
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
  if (clear) {
//	memset(data+off, 0, slots*sizeof(node_handle));
	node_handle tv = getParent()->getTransparentNode();
	for(int i=0; i<slots; i++){
		data[off+i]=tv;
	}
  }
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

  FILE_output myout(stdout);
  dumpInternal(myout, 0x03);
  
  MEDDLY_DCASSERT(0);
}
#endif

// ******************************************************************
// *                                                                *
// *                                                                *
// *                   simple_grid_style  methods                   *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_grid_style::simple_grid_style(const char* n)
 : node_storage_style(n)
{
}

MEDDLY::simple_grid_style::~simple_grid_style()
{
}

MEDDLY::node_storage* MEDDLY::simple_grid_style
::createForForest(expert_forest* f, const memory_manager_style*) const
{
  return new simple_storage("simple_grid", f, new hm_grid);
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                   simple_array_style methods                   *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_array_style::simple_array_style(const char* n)
 : node_storage_style(n)
{
}

MEDDLY::simple_array_style::~simple_array_style()
{
}

MEDDLY::node_storage* MEDDLY::simple_array_style
::createForForest(expert_forest* f, const memory_manager_style*) const
{
  return new simple_storage("simple_array", f, new hm_array);
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                   simple_heap_style  methods                   *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_heap_style::simple_heap_style(const char* n)
 : node_storage_style(n)
{
}

MEDDLY::simple_heap_style::~simple_heap_style()
{
}

MEDDLY::node_storage* MEDDLY::simple_heap_style
::createForForest(expert_forest* f, const memory_manager_style*) const
{
  return new simple_storage("simple_heap", f, new hm_heap);
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                   simple_none_style  methods                   *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::simple_none_style::simple_none_style(const char* n)
 : node_storage_style(n)
{
}

MEDDLY::simple_none_style::~simple_none_style()
{
}

MEDDLY::node_storage* MEDDLY::simple_none_style
::createForForest(expert_forest* f, const memory_manager_style*) const
{
  return new simple_storage("simple_none", f, new hm_none);
}


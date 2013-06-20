
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

#ifndef OLD_SCHEME_H
#define OLD_SCHEME_H

#include "../defines.h"
#include "../hash_stream.h"

namespace MEDDLY {
  // node storage mechanism used for versions < 0.10 of the library
  class old_node_storage;
};


/** Original node storage mechanism in a forest.

    Limits: offsets must be ints, because the data structure
    for holes uses ints.

    Details of node storage are left to the derived forests.
    However, every active node is stored in the following format.

    Below, slots are described as "ints"; note that some fields
    are stored as longs, which could require 2 slots.
    The long fields must be aligned properly, so we may need to
    pad the node with an extra, empty slot to take care of this.
      

      common   {  slot[0] : long, incoming count, >= 0.
      header --{  slot[1] : long, next pointer in unique table or special value.
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

                { -padlen : (long) Any node padding to allow for future
                {         : expansion, or for memory management purposes
                {         : (e.g., memory hold is larger than requested).
      padding --{         : padlen is number of padded slots.  If the
                {         : first entry after the node proper is negative,
                {         : then it specifies the number of long padding
                {         : slots; otherwise, there is no padding.

      tail    --{ slot[L] : long, the forest node number,
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


    Hole management.
    ==============================================================
    There are two kinds of holes depending on their location in the grid:
    Index Holes and Non-index Holes.
    Rows in the grid correspond to holes of the same size.
    The left-most column of the grid is a (vertical) list,
    and these are the index holes.
    Nodes in the middle are not connected vertically,
    and these are the non-index holes.

    The hole grid structure:
    ------------------------
    (holes_bottom)
    holes_of_size_0 (index) -- (non_index) -- (non_index) -- NULL
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
    [0] long, -size (number of slots in hole)     
    [1] up
    [2] down 
    [3] next pointer (nodes of same size)
    [4..size-2] Unused
    :
    :
    [size-1] long, -size

    Non-index holes are represented as follows:
    [0] long, -size (number of slots in hole)     
    [1] flag (<0, indicates non-index node)
    [2] prev pointer (nodes of same size)
    [3] next pointer (nodes of same size)
    [4..size-2] Unused
    :
    :
    [size-1] long, -size

*/
class MEDDLY::old_node_storage : public node_storage {
  // required interface
  public:
    old_node_storage();
    virtual ~old_node_storage();

    virtual node_storage* createForForest(expert_forest* f) const;
    virtual void collectGarbage(bool shrink);
    virtual void reportMemoryUsage(FILE* s, const char* pad, int vL) const;

    virtual void showNode(FILE* s, long addr, bool verb) const;

    virtual long makeNode(long p, const node_builder &nb, 
        node_storage_flags opt);

    virtual void unlinkDownAndRecycle(long addr);

    virtual bool areDuplicates(long addr, const node_builder &nb) const;
    virtual bool areDuplicates(long addr, const node_reader &nr) const;
    virtual void fillReader(long addr, node_reader &nr) const;
    virtual unsigned hashNode(const node_header& p) const;
    virtual int getSingletonIndex(long addr, long &down) const;
    virtual long getDownPtr(long addr, int index) const;
    virtual void getDownPtr(long addr, int ind, int& ev, long& dn) const;
    virtual void getDownPtr(long addr, int ind, float& ev, long& dn) const;
    virtual int getUnhashedHeaderOf(long addr, int ind) const;
    virtual int getHashedHeaderOf(long addr, int ind) const;

  protected:
    virtual void dumpInternalInfo(FILE*) const;
    virtual long dumpInternalNode(FILE*, long addr) const;


  private:
      static const int slots_per_long = sizeof(long) / sizeof(int);

      static const int min_size = 1024;

      /// Special values
      static const int non_index_hole = -2;
      static const long temp_node_value = -5;

      // long header indexes (relative to chunk start, as a long*)
      static const int count_index = 0;
      static const int next_index = 1;    

      // int header indexes (relative to chunk start, as an int*)
      static const int size_index = 2 * slots_per_long;
      static const int intHeaderSize = size_index+1;

      // Counts for extra slots
      static const int longTailSize = 1;
      static const int intTailSize = longTailSize * slots_per_long;
      static const int intExtra = intHeaderSize + intTailSize;

      // hole indexes
      static const int hole_up_index = slots_per_long;
      static const int hole_down_index = hole_up_index+1;
      static const int hole_prev_index = hole_down_index;
      static const int hole_next_index = hole_prev_index+1;

      /// data array
      long* data;
      /// Size of data array.
      int size;
      /// Last used data slot.  Also total number of longs "allocated"
      int last;

  // Holes grid info
  private:
      /// Largest hole ever requested
      int max_request;
      /// List of large holes
      int large_holes;
      /// Pointer to top of holes grid
      int holes_top;
      /// Pointer to bottom of holes grid
      int holes_bottom;
      /// Total ints in holes
      int hole_slots;

  // header sizes; vary by forest.
  private:
      /// Size of each outgoing edge's value (can be 0).
      char edgeSize;
      /// Size of extra unhashed data (typically 0).
      char unhashedHeader;
      /// Size of extra hashed data (typically 0).
      char hashedHeader;


  // --------------------------------------------------------
  // |  helpers for inlines.
  private:
      inline int* chunkOf(int addr) const {
        MEDDLY_DCASSERT(data);
        MEDDLY_CHECK_RANGE(1, addr, last+1);
        MEDDLY_DCASSERT(data[addr] >= 0);  // it's not a hole
        return (int*)(data + addr);
      }
      inline int* holeOf(int addr) const {
        MEDDLY_DCASSERT(data);
        MEDDLY_CHECK_RANGE(1, addr, last+1);
        MEDDLY_DCASSERT(data[addr] < 0);  // it's a hole
        return (int*)(data + addr);
      }
      inline int& rawSizeOf(int addr) const {
          return chunkOf(addr)[size_index];
      }
      inline int  sizeOf(int addr) const        { return rawSizeOf(addr); }
      inline void setSizeOf(int addr, int sz)   { rawSizeOf(addr) = sz; }
      inline int* UH(int addr) const {
          return chunkOf(addr) + intHeaderSize;
      }
      inline int* HH(int addr) const {
          return chunkOf(addr) + intHeaderSize + unhashedHeader;
      }
      // full down, as a pointer
      inline int* FD(int addr) const {
          MEDDLY_DCASSERT(rawSizeOf(addr)>0);
          return HH(addr) + hashedHeader;
      }
      // full edges, as a pointer
      inline int* FE(int addr) const {
          return FD(addr) + rawSizeOf(addr);
      }
      // a particular full edge, as a pointer
      inline void* FEP(int addr, int p) const {
        return FE(addr) + p * edgeSize;
      }
      // sparse down, as a pointer
      inline int* SD(int addr) const {
          MEDDLY_DCASSERT(rawSizeOf(addr)<0);
          return HH(addr) + hashedHeader;
      }
      // sparse indexes, as a pointer
      inline int* SI(int addr) const {
          return SD(addr) - rawSizeOf(addr);  // sparse size is negative
      }
      // sparse edges, as a pointer
      inline int* SE(int addr) const {
          return SD(addr) - 2*rawSizeOf(addr);  // sparse size is negative
      }
      // a particular sparse edge, as a pointer
      inline void* SEP(int addr, int p) const {
        return SE(addr) + p * edgeSize;
      }
      // binary search for an index
      inline int findSparseIndex(int addr, int i) const {
        int low = 0;
        int nnz = -rawSizeOf(addr);
        MEDDLY_DCASSERT(nnz>=0);
        int high = nnz;
        int* index = SI(addr);
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
      inline int intSlotsForNode(int sz) const {
          int edges = (sz<0) ? (2+edgeSize) * (-sz) : (1+edgeSize) * sz;
          return intExtra + unhashedHeader + hashedHeader + edges;
      }
      /// How many long slots would be required for a node with given size.
      ///   @param  sz  negative for sparse storage, otherwise full.
      inline int longSlotsForNode(int sz) const {
        int is = intSlotsForNode(sz);
        is += is % slots_per_long;
        return is / slots_per_long;
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
      long makeFullNode(long p, int size, const node_builder &nb);

      /** Create a new node, stored sparsely.
          Space is allocated for the node, and data is copied.
            @param  p     Node handle number.
            @param  size  Number of nonzero downward pointers.
            @param  nb    Node data is copied from here.
            @return       The "address" of the new node.
      */
      long makeSparseNode(long p, int size, const node_builder &nb);

      void copyExtraHeader(long addr, const node_builder &nb);

  // --------------------------------------------------------
  // |  Hole management helpers.
  private:
      inline int& h_up(int off) const {
        return holeOf(off)[hole_up_index];
      }
      inline int& h_down(int off) const {
        return holeOf(off)[hole_down_index];
      }
      inline int& h_prev(int off) const {
        return holeOf(off)[hole_prev_index];
      }
      inline int& h_next(int off) const {
        return holeOf(off)[hole_next_index];
      }

      inline int& holeUp(int off)       { return h_up(off); }
      inline int  holeUp(int off) const { return h_up(off); }

      inline int& holeDown(int off)       { return h_down(off); }
      inline int  holeDown(int off) const { return h_down(off); }

      inline int& holePrev(int off)       { return h_prev(off); }
      inline int  holePrev(int off) const { return h_prev(off); }

      inline int& holeNext(int off)       { return h_next(off); }
      inline int  holeNext(int off) const { return h_next(off); }

      inline bool isHoleNonIndex(int p_offset) const {
          return (non_index_hole == h_up(p_offset));
      }

      /// Find actual number of slots used for this active node.
      inline long activeNodeActualLongSlots(long off) const {
          long end = off + longSlotsForNode(sizeOf(off))-1;
          // account for any padding
          if (data[end] < 0) {
            end -= data[end];
          }
          return end - off + 1;
      }

      // returns offset to the hole found in level
      int getHole(int slots);

      // makes a hole of size == slots, at the specified offset
      void makeHole(int p_offset, int slots);

      // add a hole to the hole grid
      void gridInsert(int p_offset);

      // remove a non-index hole from the hole grid
      void midRemove(int p_offset);

      // remove an index hole from the hole grid
      void indexRemove(int p_offset);

      // resize the data array.
      void resize(int new_slots);

      /** Allocate enough slots to store a node with given size.
          Also, stores the node size in the node.
            @param  sz      negative for sparse storage, otherwise full.
            @param  tail    Node id
            @param  clear   Should the node be zeroed.
            @return         Offset in the data array.
      */
      long allocNode(int sz, long tail, bool clear);

  // --------------------------------------------------------
  // |  Node comparison as a template
  private:
    template <class nodetype>
    inline bool areDupsTempl(long addr, const nodetype &n) const {
      int size = sizeOf(addr);
      if (size<0) {
        //
        // Node is sparse
        //
        int nnz = -size;
        int* down = SD(addr);
        int* index = SI(addr);
        if (n.isFull()) {
          // check that down matches
          int i = 0;
          for (int z=0; z<nnz; z++) {
            if (index[z] >= n.getSize()) return false;
            for (; i<index[z]; i++) if (n.d(i)) return false;
            if (n.d(i) != down[z]) return false;
            i++;
          }
          for (; i<n.getSize(); i++) if (n.d(i)) return false;
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
              if (!getParent()->areEdgeValuesEqual( 
                      SEP(addr, z),   n.eptr(z)
                  )) return false;
          }
        } 
        // must be equal
        return true;
      }
      //
      // Node is truncated full
      //
      int* down = FD(addr);
      if (n.isFull()) {
        if (size > n.getSize()) return false;
        // check down
        int i;
        for (i=0; i<size; i++) {
          if (down[i] != n.d(i)) return false;
        }
        for ( ; i<n.getSize(); i++) {
          if (n.d(i)) return false;
        }
        // check edges
        if (n.hasEdges()) {
          for (int i=0; i<size; i++) if (down[i]) {
              if (!getParent()->areEdgeValuesEqual( 
                      FEP(addr, i),   n.eptr(i)
                  )) return false;
          }
        }
        // must be equal
        return true;
      }
      // n is sparse
      int i = 0;
      // check down
      for (int z=0; z<n.getNNZs(); z++) {
        if (n.i(z) >= size) return false;
        for (; i<n.i(z); i++) if (down[i]) return false;
        if (n.d(z) != down[i]) return false;
        i++;
      }
      if (i<size) return false; // there WILL be a non-zero down
      // check edges
      if (n.hasEdges()) {
        for (int z=0; z<n.getNNZs(); z++) {
          if (!getParent()->areEdgeValuesEqual(
                FEP(addr, n.i(z)),  n.eptr(z)
              )) return false;
        } // for z
      }
      // must be equal
      return true;
    }
}; 

#endif


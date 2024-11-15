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

#include "best.h"
#include "../defines.h"
#include "../unpacked_node.h"
#include "../memory.h"
#include "../forest.h"
#include "../operators.h"
#include "../node_marker.h"
#include <set>
#include <map>
#include <cassert>

#define MAX_PATTERN_LEN 10

#include "../hash_stream.h"

// #define DEBUG_ENCODING
// #define DEBUG_DECODING
// #define VERIFY_ENCODING
// #define VERIFY_DECODING
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
  class best_storage;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                     best_storage class                     *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Original node storage mechanism in a forest.
 Memory management is completely separated from the class,
 which makes this implementation a little different
 from the original.

 Details of node storage are left to the derived forests.
 However, every active node is stored in the following format.



 common   {  slot[0] : next pointer >=0 in unique table.
 header --{  slot[1] : Bit 0 (LSB) : 1 for extensible, 0 for not.
 Bit 1       : 1 for sparse, 0 for full.
 Bit MSB..2  : size

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
class MEDDLY::best_storage : public node_storage {
public:
  best_storage(const char* n, forest* f, const memory_manager_style* mst,
          memstats &ms);
  virtual ~best_storage();

  // required interface
public:
  virtual void collectGarbage(bool shrink);
  virtual void reportStats(output &s, const char* pad, reporting_flags flags)
      const;

  virtual node_address makeNode(node_handle p, const unpacked_node &nb,
                                node_storage_flags opt);

  virtual void unlinkDownAndRecycle(node_address addr);
  virtual void addDownToQueue(node_marker &m, node_address addr) const;

  virtual bool areDuplicates(node_address addr, const unpacked_node &nr) const;
  virtual void fillUnpacked(unpacked_node &nr, node_address addr, node_storage_flags) const;

  virtual unsigned hashNode(int level, node_address addr) const;
  virtual bool isSingletonNode(node_address addr, unsigned &index,
          node_handle &down) const;

#ifdef ALLOW_DEPRECATED_SINGLETON
  virtual int getSingletonIndex(node_address addr, node_handle &down) const;
#endif
  virtual node_handle getDownPtr(node_address addr, int index) const;
  virtual void getDownPtr(node_address addr, int ind, int& ev, node_handle& dn) const;
  virtual void getDownPtr(node_address addr, int ind, long& ev, node_handle& dn) const;
  virtual void getDownPtr(node_address addr, int ind, float& ev, node_handle& dn) const;
  virtual int getExtensibleIndex(node_address addr) const;
  virtual const void* getUnhashedHeaderOf(node_address addr) const;
  virtual const void* getHashedHeaderOf(node_address addr) const;

  virtual node_handle getNextOf(node_address addr) const;
  virtual void setNextOf(node_address addr, node_handle n);

protected:
  virtual void updateData(node_handle* d);
  virtual node_address firstNodeAddress() const;
  virtual void dumpInternalInfo(output &) const;
  virtual node_address
  dumpInternalNode(output &, node_address addr, unsigned flags) const;
  virtual void dumpInternalTail(output &) const;


  // --------------------------------------------------------
  // |  Misc. helpers.
private:
  /// How many int slots would be required for a node with given size.
  ///   @param  sz      Number of downward pointers.
  ///   @param  sparse  True for sparse storage, otherwise full.
  inline int slotsForNode(int sz, bool sparse) const {
    int node_slots = sparse ? (2+slots_per_edge) * sz : (1+slots_per_edge) * sz;
    return extra_slots + unhashed_slots + hashed_slots + node_slots;
  }
  inline int slotsForNode(int sz) const {
    int node_slots = sz;
    return extra_slots + unhashed_slots + hashed_slots + node_slots;
  }
  inline node_handle* getChunkAddress(node_address addr) const {
    MEDDLY_DCASSERT(MM);
    return (node_handle*) MM->getChunkAddress(addr);
  }
  inline static unsigned int getRawSize(const node_handle* chunk) {
    return chunk[size_or_id_slot];
  }
  inline static unsigned int getRawSize(int size, bool sparse, bool extensible, bool pattern) {
    return (((unsigned int)(size) << 3) | (sparse? 4U: 0) | (extensible? 2U: 0) | (pattern? 1U: 0));
  }
  inline static unsigned int getSize(unsigned int raw_size) {  // Returns size or ID
    return ((unsigned int)raw_size) >> 3;
  }
  inline static bool isSparse(unsigned int raw_size) {
    return ((raw_size & 4U) != 0);
  }
  inline static bool isExtensible(unsigned int raw_size) {
    return ((raw_size & 2U) != 0);
  }
  inline static bool isPattern(unsigned int raw_size) {
    return ((raw_size & 1U) != 0);
  }
  virtual bool isExtensible(node_address addr) const {
    MEDDLY_DCASSERT(MM);
    const node_handle* chunk = getChunkAddress(addr);
    MEDDLY_DCASSERT(chunk);
    const unsigned int raw_size = getRawSize(chunk);
    return isExtensible(raw_size);
  }

  //
  // binary search for an index
  //
  inline static int findSparseIndex(int i, const node_handle* index, int N) {
    int low = 0;
    int high = N;
    while (low < high) {
      int z = (low+high)/2;
      MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, z, N);
      if (index[z] == i) return z;
      if (index[z] < i) low = z + 1;
      else              high = z;
    }
    return -1;
  }
  // binary search for an index
  inline int findSparseIndex(int i, const node_handle* chunk) const {
    const unsigned int raw_size = getRawSize(chunk);
    MEDDLY_DCASSERT(isSparse(raw_size));
    int nnz = getSize(raw_size);
    const node_handle* index = chunk + down_start + nnz;
    if (isExtensible(raw_size) && i >= index[nnz-1]) return nnz-1;
    return findSparseIndex(i, index, nnz);
  }

  /** Create a new node, stored as truncated full.
   Space is allocated for the node, and data is copied.
   @param  p     Node handle number.
   @param  size  Number of downward pointers.
   @param  nb    Node data is copied from here.
   @return       The "address" of the new node.
   */
  node_address makeFullNode(node_handle p, int size, const unpacked_node &nb);

  /** Create a new node, stored sparsely.
   Space is allocated for the node, and data is copied.
   @param  p     Node handle number.
   @param  size  Number of nonzero downward pointers.
   @param  nb    Node data is copied from here.
   @return       The "address" of the new node.
   */
  node_address makeSparseNode(node_handle p, int size,
                              const unpacked_node &nb);

  /** Generate pattern from pattern-identifier.
   @param  identifier     Unique Identifier of some pattern.
   @return                The truncated full pattern.
   */
  std::string  generatePatternFromIndex(node_handle identifier) const;

  /** Generate pattern from node.
   @param  nb    Node data is copied from here.
   @return       The "pattern" of the node of length MAX_PATTERN_LEN.
   */
  std::string  generatePatternFromNode(const unpacked_node &nb) const;

  /** Generate pattern-identifier from node
   Pattern is generated for the node, and pattern-identifier is calculated.
   @param  nb    Node data is copied from here.
   @return       The pattern-identifier of the node.
   */
  node_handle generateIndexFromNode(const unpacked_node &nb);

  /** Reverse a pattern to ensure all full-truncated patterns have same pattern-identifier.
   The formula for identifier calculation generates same numbers for patterns of type "t(*)X"
   The needful is to have same identfier for "X(t*)",
   where X is a pattern formed by downward pointers of node and 't' is transaparent-value
   Hence, the need for reversal.
   @param  pattern     Raw-Pattern as obtained from a node.
   @return             Pattern satisfying "unique-identifier-full-truncated" condition
   */
  std::string specialReverse(std::string pattern) const;

  /** Create a new node, stored as pattern.
   Space is allocated for the node, and data is copied.
   @param  p     Node handle number.
   @param  size  Number of unique downward pointers.
   @param  nb    Node data is copied from here.
   @return       The "address" of the new node.
   */
  node_address makePatternNode(node_handle p, int size, const unpacked_node &nb, const std::map<int,node_handle> nhmap);

private:
  memory_manager* MM;

  //
  // Header indexes that are fixed
  //
  // static const int count_slot = 0;
  static const int next_slot = 0;
  static const int size_or_id_slot = 1;
  static const int header_slots = size_or_id_slot+1;

  // Slots at the end
  static const int tail_slots = 1;
  static const int extra_slots = header_slots + tail_slots;

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
// *                    best_storage methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::best_storage
::best_storage(const char* n, forest* f, const memory_manager_style* mst,
        memstats &ms)
: node_storage(n, f)
{
  MM = mst->initManager(sizeof(node_handle), slotsForNode(0, false), ms);

  unhashed_start = header_slots;
  unhashed_slots = slotsForBytes(f->unhashedHeaderBytes());
  hashed_start = unhashed_start + unhashed_slots;
  hashed_slots = slotsForBytes(f->hashedHeaderBytes());
  down_start =  hashed_start + hashed_slots;
    switch (f->getEdgeType()) {
        case edge_type::VOID:
            slots_per_edge = 0;
            break;

        case edge_type::INT:
            slots_per_edge = slotsForBytes(sizeof(int));
            break;

        case edge_type::LONG:
            slots_per_edge = slotsForBytes(sizeof(long));
            break;

        case edge_type::FLOAT:
            slots_per_edge = slotsForBytes(sizeof(float));
            break;

        case edge_type::DOUBLE:
            slots_per_edge = slotsForBytes(sizeof(double));
            break;

        default:
           MEDDLY_DCASSERT(false);
    }
}

MEDDLY::best_storage::~best_storage()
{
  // TBD - special steps to recycle all nodes?
  delete MM;
}

void MEDDLY::best_storage::collectGarbage(bool shrink)
{
  // TBD
}

void MEDDLY::best_storage::reportStats(output &s, const char* pad,
                                       reporting_flags flags) const
{
  static unsigned STORAGE = STORAGE_STATS | STORAGE_DETAILED;

  if (flags & STORAGE) {
    // s << pad << "Stats for " << getStyleName() << "\n";

  }

  static unsigned HOLEMAN = HOLE_MANAGER_STATS | HOLE_MANAGER_DETAILED;

  if (flags & HOLEMAN) {
    MM->reportStats(s, pad,
                    flags & HUMAN_READABLE_MEMORY,
                    flags & HOLE_MANAGER_DETAILED);
  }


  /*
   #ifdef DEVELOPMENT_CODE
   verifyStats();
   #endif
   */
}

MEDDLY::node_address MEDDLY::best_storage
::makeNode(node_handle p, const unpacked_node &nb, node_storage_flags opt)
{
#ifdef DEBUG_ENCODING
  printf("best_storage making node\n        temp:  ");
  FILE_output out(stdout);
  nb.show(out, true);
#endif
  //
  // Determine truncsize and nnzs
  //
  unsigned nnzs = 0;
  unsigned truncsize = 0;
  if (nb.isSparse()) {
    //
    // nb is sparse, easy to get nnzs
    //
    nnzs = nb.getSize();
    MEDDLY_DCASSERT(nb.isSorted());
    truncsize = nb.index(nnzs-1)+1;
  } else {
    //
    // nb is full, need to scan it and count
    //
    if (slots_per_edge) {
      //
      // EV, check for transparent edges
      //
      for (unsigned i=0; i<nb.getSize(); i++) {
        if (!parent->isTransparentEdge(nb.down(i), nb.edgeval(i))) {
          nnzs++;
          truncsize = i+1;
        }
      }
    } else {
      //
      // MT, check for pointers to transparent nodes
      //
      const node_handle tv = parent->getTransparentNode();
      for (unsigned i=0; i<nb.getSize(); i++) {
        if (nb.down(i) != tv) {
          nnzs++;
          truncsize = i+1;
        }
      }
    }
  }
  //
  // Determine unique non-zero pointers : uniqnnzs
  //
  int uniqnnzs = 0;
  std::map<int,node_handle> nhmap;
  std::set<node_handle> nhset;

  if( ( !nb.isExtensible() && !nb.hasEdges() ) && ( (nb.isFull() && (nb.getSize()<MAX_PATTERN_LEN)) || (nb.isSparse() && (nb.index(nb.getSize()-1)+1<MAX_PATTERN_LEN)) ) ){
    if (nb.isSparse()) {
      //
      // nb is sparse
      //
      MEDDLY_DCASSERT(nb.index(nb.getSize()-1)+1<MAX_PATTERN_LEN);
      int nhseq = 1;
      for (unsigned i=0; i<nb.getSize(); i++)
        {
        int old_size = nhset.size();
        nhset.insert(nb.down(i));
        if(nhset.size()>old_size) // if set size increases, then this node_handle is unique
          {
          nhmap.insert(std::pair<int,node_handle>(nhseq,nb.down(i))); // set contains only unique pointers
          nhseq+=1;
          }
        }
      uniqnnzs = nhmap.size();
    } else {
      //
      // nb is full and MT
      //
      MEDDLY_DCASSERT(nb.getSize()<MAX_PATTERN_LEN);

      if (!slots_per_edge) {
        const node_handle tv = parent->getTransparentNode();

        int nhseq = 1;
        for (int i=0; i<nb.getSize(); i++) {
          if (nb.down(i) != tv) {
            int old_size = nhset.size();
            nhset.insert(nb.down(i));
            if(nhset.size()>old_size) // if set size increases, then this node_handle is unique
              {
              nhmap.insert(std::pair<int,node_handle>(nhseq,nb.down(i))); // insert only the unique non-tv pointers
              nhseq+=1;
              }
          }
        }
        uniqnnzs = nhmap.size();
      }
    }
  }
  std::cout<<"\n nb.hasEdges"<<nb.hasEdges()<<" nb.isExtensible()"<<nb.isExtensible()<<" uniqnnzs"<<uniqnnzs;

  node_address addr = 0;

  if (FULL_ONLY == opt) {
    //
    // Sparse nodes disabled?
    //
    if (uniqnnzs == 0)
      //
      // Pattern not possible? Just build a full one
      //
      addr = makeFullNode(p, truncsize, nb);
    else if( slotsForNode(truncsize, false) < slotsForNode(uniqnnzs) )
      //
      // Full node is more compact
      //
      addr = makeFullNode(p, truncsize, nb);
    else
      //
      // Pattern is more compact
      //
      addr = makePatternNode(p, uniqnnzs, nb, nhmap);
  } else if (SPARSE_ONLY == opt) {
    //
    // Full nodes disabled?
    //
    if (uniqnnzs == 0)
      //
      // Pattern not possible? Just build a sparse one
      //
      addr = makeSparseNode(p, nnzs, nb);
    else if( slotsForNode(nnzs, true) < slotsForNode(uniqnnzs) )
      //
      // Sparse node is more compact
      //
      addr = makeSparseNode(p, nnzs, nb);
    else
      //
      // Pattern is more compact
      //
      addr = makePatternNode(p, uniqnnzs, nb, nhmap);
  } else {
    //
    // Full, sparse and pattern are allowed, determine which is more compact
    //
    if (slotsForNode(truncsize, false) < slotsForNode(nnzs, true)) {
      //
      // Full nodes are more compact
      //
      if (uniqnnzs == 0)
        //
        // Pattern not possible? Just build a full one
        //
        addr = makeFullNode(p, truncsize, nb);
      else if( slotsForNode(truncsize, false) < slotsForNode(uniqnnzs) )
        //
        // Full node is more compact
        //
        addr = makeFullNode(p, truncsize, nb);
      else
        //
        // Pattern is more compact
        //
        addr = makePatternNode(p, uniqnnzs, nb, nhmap);
    } else {
      //
      // Sparse nodes are more compact
      //
      if (uniqnnzs == 0)
        //
        // Pattern not possible? Just build a sparse one
        //
        addr = makeSparseNode(p, nnzs, nb);
      else if( slotsForNode(nnzs, true) < slotsForNode(uniqnnzs) )
        //
        // Sparse node is more compact
        //
        addr = makeSparseNode(p, nnzs, nb);
      else
        //
        // Pattern is more compact
        //
        addr = makePatternNode(p, uniqnnzs, nb, nhmap);
    }
  }


#ifdef VERIFY_ENCODING
  if (!areDuplicates(addr, nb)) {
    FILE_output out(stdout);
    printf("Error encoding unpacked node: ");
    nb.show(out, true);
    printf("packed: ");
    dumpInternalNode(out, addr, 0x03);
    printf("\n");
    fflush(stdout);
  }
#endif
  MEDDLY_DCASSERT(areDuplicates(addr, nb));
  return addr;
}



void MEDDLY::best_storage::unlinkDownAndRecycle(node_address addr)
{
#ifdef MEMORY_TRACE
  printf("recycling node at address %ld\n", addr);
  FILE_output out(stdout);
  dumpInternalNode(out, addr, 0x03);
#endif
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  const unsigned int raw_size = getRawSize(chunk);
  const unsigned int size = getSize(raw_size);
  const unsigned int is_sparse = isSparse(raw_size);
  const unsigned int is_pattern = isPattern(raw_size);

  //
  // Unlink down pointers
  //
  const node_handle* down = chunk + down_start;
  int uniqnnzs = 0;
  if(is_pattern)
    {
    std::string pattern_from_index = generatePatternFromIndex(size);

    int trunc_pattern_size = 0;
    std::set<char> uniqnh;
    for(int i=0;i<MAX_PATTERN_LEN;i++)
      {
      if(pattern_from_index[i]!='t')
        {
        trunc_pattern_size = i+1;
        uniqnh.insert(pattern_from_index[i]);
        }
      uniqnnzs = uniqnh.size();
      }

    for (unsigned int i=0; i<trunc_pattern_size; i++) {
      if(pattern_from_index[i]!='t')
        parent->unlinkNode(down[pattern_from_index[i]-'A']);
    }
    } else {
      for (unsigned int i=0; i<size; i++) {
        parent->unlinkNode(down[i]);
      }
    }

  // Can this move after unlinking?
  MEDDLY_DCASSERT(MM->getChunkAddress(addr) == chunk);

  //
  // Determine number of slots in this node
  //

  size_t actual_slots = is_pattern ? slotsForNode(uniqnnzs) : slotsForNode(size, is_sparse);
  if (chunk[actual_slots-1] < 0) {
    // padding
    actual_slots += (-chunk[actual_slots-1]);
  }

  //
  // Recycle
  //
  MM->recycleChunk(addr, actual_slots);
}



void MEDDLY::best_storage::addDownToQueue(node_marker &nm, node_address addr) const
{
#ifdef DEBUG_MARK_SWEEP
  printf("marking children at address %ld\n", addr);
  FILE_output out(stdout);
  dumpInternalNode(out, addr, 0x03);
#endif
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  const unsigned int raw_size = getRawSize(chunk);
  const unsigned int size = getSize(raw_size);
  // const unsigned int is_sparse = isSparse(raw_size);
  const unsigned int is_pattern = isPattern(raw_size);

  //
  // Mark down pointers
  //
  const node_handle* down = chunk + down_start;
  // int uniqnnzs = 0;
  if(is_pattern)
    {
    std::string pattern_from_index = generatePatternFromIndex(size);

    int trunc_pattern_size = 0;
    std::set<char> uniqnh;
    for(int i=0;i<MAX_PATTERN_LEN;i++)
      {
      if(pattern_from_index[i]!='t')
        {
        trunc_pattern_size = i+1;
        uniqnh.insert(pattern_from_index[i]);
        }
      // uniqnnzs = uniqnh.size();
      }

    for (unsigned int i=0; i<trunc_pattern_size; i++)
      {
      if(pattern_from_index[i]!='t')
        nm.addToQueue(down[pattern_from_index[i]-'A']);
      }
    } else {
      for (unsigned int i=0; i<size; i++) {
        nm.addToQueue(down[i]);
      }
    }
}


bool MEDDLY::best_storage
::areDuplicates(node_address addr, const unpacked_node &n) const
{
  MEDDLY_DCASSERT(n.isTrim());
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  MEDDLY_DCASSERT( n.hasEdges() == (slots_per_edge>0) );

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


  const unsigned int raw_size = getRawSize(chunk);
  const unsigned int size = getSize(raw_size);

  const bool is_extensible = isExtensible(raw_size);
  const bool is_pattern = isPattern(raw_size);
  const bool is_sparse = isSparse(raw_size);

  if (n.isExtensible()) {
    if (!is_extensible) return false;
  } else {
    if (is_extensible) return false;
  }

  const node_handle tv = parent->getTransparentNode();
  const node_handle* down = chunk + down_start;


  if(is_pattern)
    {
    //retrieve the actual pattern formed by a node and represented by an index_no.
    std::string pattern_from_node =  specialReverse(generatePatternFromNode(n));
    std::string pattern_from_index = generatePatternFromIndex(size);

    int trunc_pattern_size = 0;
    for(int i=0;i<MAX_PATTERN_LEN;i++)
      {
      if(pattern_from_index[i]!='t')
        trunc_pattern_size = i+1;
      }
    const node_handle tv = parent->getTransparentNode();
    if(n.isFull())
      {
      //
      // n is full
      //
      if (trunc_pattern_size > unsigned(n.getSize())) return false;
      // check down
      unsigned int i;
      for (i=0; i<trunc_pattern_size; i++) {
        if ( ( (n.down(i) != tv) && (down[pattern_from_index[i]-'A'] != n.down(i) ) )
            || ( (n.down(i) == tv) && (pattern_from_index[i]!='t') )
            ) return false;
      }

      for (; i<unsigned(n.getSize()); i++) {
        if (n.down(i)!=tv) return false;
      }
      //must be equal
      return true;
      }else {
        //
        // n is sparse
        //
        int i = 0;
        // check down
        for (unsigned z=0; z<n.getSize(); z++) {
          if (unsigned(n.index(z)) >= trunc_pattern_size) return false;
          // loop - skipped edges must be transparent
          for (; i<n.index(z); i++) {
            if (pattern_from_index[i]!='t') return false;
          }
          if (n.down(z) != down[pattern_from_index[i]-'A']) return false;
          i++;
        } // for z
        if (unsigned(i)<trunc_pattern_size) return false; // there WILL be a non-zero down

        //must be equal
        return true;
      }
    }else if (is_sparse) {
      //
      // Node is stored sparsely
      //
      const unsigned int nnz = size;
      const node_handle* index = down + nnz;
      const node_handle* edge = slots_per_edge ? (index + nnz) : 0;

      if (n.isFull()) {
        //
        // n is full
        //
        node_handle i = 0;
        for (unsigned int z=0; z<nnz; z++) {
          if (index[z] >= n.getSize()) return false;  // too large
                                                      // loop - skipped edges must be transparent
          for (; i<index[z]; i++) {
            if (n.down(i)!=tv) return false;
          } // for
            // compare this[z] with n[i]
          if (down[z] != n.down(i)) return false;
          i++;
        } // for z
          //
          // Anything beyond must be transparent
        for (; i<n.getSize(); i++) {
          if (n.down(i)!=tv) {
            return false;
          }
        }
        //
        // now check edges
        //
        if (edge) {
          for (unsigned int z=0; z<nnz; z++) {
            if (! n.edgeval(index[z]).equals(edge + z*slots_per_edge)) {
              return false;
            }
          } // for z
        }
        // must be equal
        return true;
      } else {
        //
        // n is sparse
        //
        if (n.getSize() != nnz) return false;
        // check that down matches
        for (unsigned int z=0; z<nnz; z++) {
          if (index[z] != n.index(z)) return false;
          if (down[z] != n.down(z))  return false;
        }
        // check that edges match
        if (n.hasEdges()) {
          for (unsigned int z=0; z<nnz; z++) {
            if (! n.edgeval(z).equals(edge + z*slots_per_edge)) {
              return false;
            }
          }
        }
        // must be equal
        return true;
      }
    } else{
      //
      // Node is stored truncated full
      //
      const node_handle* edge = slots_per_edge ? (down + size) : 0;

      if (n.isFull()) {
        //
        // n is full
        //
        if (size > unsigned(n.getSize())) return false;
        // check down
        unsigned int i;
        for (i=0; i<size; i++) {
          if (down[i] != n.down(i)) return false;
        }
        for (; i<unsigned(n.getSize()); i++) {
          if (n.down(i)!=tv) return false;
        }
        // check edges
        if (n.hasEdges()) {
          for (unsigned int i=0; i<size; i++) if (down[i]) {
            if (! n.edgeval(i).equals(edge + i*slots_per_edge)) {
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
        unsigned i = 0;
        // check down
        for (unsigned z=0; z<n.getSize(); z++) {
          if (unsigned(n.index(z)) >= size) return false;
          // loop - skipped edges must be transparent
          for (; i<n.index(z); i++) {
            if (down[i]!=tv) return false;
          }
          if (n.down(z) != down[i]) return false;
          i++;
        } // for z
        if (unsigned(i)<size) return false; // there WILL be a non-zero down
                                            //
                                            // check edges
                                            //
        if (n.hasEdges()) {
          for (unsigned z=0; z<n.getSize(); z++) {
            if (! n.edgeval(z).equals(edge + n.index(z) * slots_per_edge)) {
              return false;
            }
          } // for z
        }
        // must be equal
        return true;
      }
    }

}


void MEDDLY::best_storage
::fillUnpacked(unpacked_node &nr, node_address addr, node_storage_flags st2) const
{
#ifdef DEBUG_DECODING
  FILE_output out(stdout);
  printf("best_storage filling reader\n    internal: ");
  dumpInternalNode(out, addr, 0x03);
#endif

  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  MEDDLY_DCASSERT(nr.hasEdges() == (slots_per_edge>0) );

  //
  // Copy any extra header information
  //

  if (unhashed_slots) {
    nr.setUHdata(chunk + unhashed_start);
    // memcpy(nr.UHdata(), chunk + unhashed_start, nr.UHbytes());
  }
  if (hashed_slots) {
    nr.setHHdata(chunk + hashed_start);
    // memcpy(nr.HHdata(), chunk + hashed_start, nr.HHbytes());
  }

  //
  // Copy everything else
  //

  const unsigned int raw_size = getRawSize(chunk);
  const unsigned int size = getSize(raw_size);
  const bool is_sparse = isSparse(raw_size);
  const bool is_extensible = isExtensible(raw_size);
  const bool is_pattern = isPattern(raw_size);
  const node_handle* down = chunk + down_start;

  /*
   Set the unpacked node storage style based on settings
   */
  switch (st2) {
    case FULL_ONLY:
        nr.setFull();
        break;

    case SPARSE_ONLY:
        nr.setSparse();
        break;

    case FULL_OR_SPARSE:
        if (is_pattern || !is_sparse) {
            nr.setFull();
        } else {
            nr.setSparse();
        }
        break;

    default:            assert(0);
  };

#ifdef ALLOW_EXTENSIBLE
  // if (is_extensible) nr.markAsExtensible(); else nr.markAsNotExtensible();
  if (is_extensible && parent->isExtensibleLevel(nr.getLevel()))
    nr.markAsExtensible();
  else
    nr.markAsNotExtensible();
#endif

  // Make sure that when an extensible node is unpacked that the trailing edges
  // are filled correctly (regardless of the storage scheme of the unpacked node)

  const node_handle tv = parent->getTransparentNode();

  if(is_pattern){

    MEDDLY_DCASSERT(parent->isExtensibleLevel(nr.getLevel()) == false);

    //pattern node is stored with size of unique pointers + the unique pointers + identifier for the pattern thus formed

    if(nr.isFull()) {
      //
      //Copying into full Node
      //

      //Obtain the pattern formed by child pointers
      std::string node_pattern = generatePatternFromIndex(size);

      //Obtain size for nr
      int trunc_pattern_size = 0;
      for(int i=0;i<MAX_PATTERN_LEN;i++)
        {
        if(node_pattern[i]!='t')
          trunc_pattern_size=i+1;
        }

      //populate nr.down

      nr.clear(0, nr.getSize());

      for(int i=0;i<trunc_pattern_size;i++)
        {
        if(node_pattern[i] != 't') //If not transparent
          nr.setFull(i, down[node_pattern[i] - 'A']);
          // nr.d_ref(i) = down[node_pattern[i] - 'A'];
        }

    } else {

      //
      //Copying into sparse node
      //

      //Obtain the pattern formed by child pointers
      std::string node_pattern = generatePatternFromIndex(size);

      //Obtain size for nr
      int trunc_pattern_size = 0;
      for(int i=0;i<MAX_PATTERN_LEN;i++)
        {
        if(node_pattern[i]!='t')
          {
          trunc_pattern_size = i+1;
          }
        }

      //populate nr.down() and nr.index()
      int z = 0;
      for(int i=0;i<trunc_pattern_size;i++)
        {
        if(node_pattern[i]!='t') //If transparent
          {
              nr.setSparse(z, i, down[node_pattern[i] - 'A']);
          // nr.i_ref(z) = i;
          // nr.d_ref(z) = down[node_pattern[i] - 'A'];
              z++;
          }
        }
      if (nr.isExtensible() == false) {nr.shrink(z);}

    }


  } else if (is_sparse) {
    //
    // Node is sparse
    //
    int nnz = size;
    const node_handle* index = down + nnz;
    const node_handle* edge = slots_per_edge ? (index + nnz) : 0;
    const int ext_i = index[nnz-1];
    const int ext_d = is_extensible? down[nnz-1]: tv;
    const void* ext_ptr = is_extensible? (edge + (nnz-1)*slots_per_edge): 0;


    if (nr.isFull()) {
      //
      // Copying into a full node
      //

      // Write the sparse edges
      for (int z=0; z<nnz; z++) {
        if (nr.hasEdges()) {
          nr.setFull(index[z], edge + z*slots_per_edge, down[z]);
        } else {
          nr.setFull(index[z], down[z]);
        }
      }
      // Write the extensible edge to the rest (trailing locations)
      for (int i=ext_i+1; i<nr.getSize(); i++) {
        if (nr.hasEdges()) {
          if (ext_ptr) {
            nr.setFull(i, ext_ptr, ext_d);
          } else {
            nr.setFull(i, parent->getTransparentEdge(), ext_d);
          }
        } else {
          nr.setFull(i, ext_d);
        }

        /*
        nr.d_ref(i) = ext_d;
        if (nr.hasEdges()) {
          if (ext_ptr)
            nr.setEdgeRaw(i, ext_ptr);
          else
            nr.set_edgeval(i, parent->getTransparentEdge());
        }
        */
      }
    } else {
      //
      // Copying into a sparse node
      //

      for (int z=0; z<nnz; z++) {
        if (nr.hasEdges()) {
            nr.setSparse(z, index[z], edge + z*slots_per_edge, down[z]);
        } else {
            nr.setSparse(z, index[z], down[z]);
        }

        /*
        nr.d_ref(z) = down[z];
        nr.i_ref(z) = index[z];
        if (nr.hasEdges()) {
          nr.setEdgeRaw(z, edge + z*slots_per_edge);
        }
        */
      } // for z
      if (nr.isExtensible() == is_extensible) {
        nr.shrink(nnz);
      } else {
        MEDDLY_DCASSERT(is_extensible);
        MEDDLY_DCASSERT(!nr.isExtensible());
        int i = ext_i+1;
        for (unsigned z=nnz; z<nr.getSize(); z++, i++) {
            if (nr.hasEdges()) {
                if (ext_ptr) {
                    nr.setSparse(z, i, ext_ptr, ext_d);
                } else {
                    nr.setSparse(z, i, parent->getTransparentEdge(), ext_d);
                }
            } else {
                nr.setSparse(z, i, ext_d);
            }
            /*
          nr.i_ref(z) = i;
          nr.d_ref(z) = ext_d;
          if (nr.hasEdges()) {
            if (ext_ptr)
              nr.setEdgeRaw(z, ext_ptr);
            else
              nr.set_edgeval(z, parent->getTransparentEdge());
          }
          */
        }
      }
    }
  } else {
    //
    // Node is full
    //
    const node_handle* edge = slots_per_edge ? (down + size) : 0;

    if (nr.isFull()) {
      //
      // Copying into a full node
      //
      unsigned i;
      for (i=0; i<size; i++) {
        if (nr.hasEdges()) {
            nr.setFull(i, edge + i*slots_per_edge, down[i]);
        } else {
            nr.setFull(i, down[i]);
        }
          /*
        nr.d_ref(i) = down[i];
        if (nr.hasEdges()) {
          nr.setEdgeRaw(i, edge + i*slots_per_edge);
        }
        */
      }
      if (FULL_OR_SPARSE == st2 && is_extensible == nr.isExtensible()) {
        nr.shrink(size);
      } else {
        const int ext_d = is_extensible? down[size-1]: tv;
        const void* ext_ptr = is_extensible? (edge + (size-1)*slots_per_edge): 0;
        for (; i<unsigned(nr.getSize()); i++) {
            if (nr.hasEdges()) {
                if (ext_ptr) {
                    nr.setFull(i, ext_ptr, ext_d);
                } else {
                    nr.setFull(i, parent->getTransparentEdge(), ext_d);
                }
            } else {
                nr.setFull(i, ext_d);
            }
            /*
          nr.d_ref(i) = ext_d;
          if (nr.hasEdges()) {
            if (ext_ptr)
              nr.setEdgeRaw(i, ext_ptr);
            else
              nr.set_edgeval(i, parent->getTransparentEdge());
          }
          */
        }
      }
    } else {
      //
      // Copying into a sparse node
      //

      unsigned z = 0;
      for (unsigned int i=0; i<size; i++) if (down[i]) {
        if (nr.hasEdges()) {
            nr.setSparse(z, i, edge + i*slots_per_edge, down[i]);
        } else {
            nr.setSparse(z, i, down[i]);
        }

          /*
        nr.d_ref(z) = down[i];
        nr.i_ref(z) = i;
        if (nr.hasEdges()) {
          nr.setEdgeRaw(z, edge + i*slots_per_edge);
        }
        */

        z++;
      } // for i
      if (nr.isExtensible() == is_extensible) nr.shrink(z);
      else {
        MEDDLY_DCASSERT(is_extensible);
        MEDDLY_DCASSERT(!nr.isExtensible());
        const int ext_i = size-1;
        const int ext_d = down[ext_i];
        const void* ext_ptr = (edge + (ext_i)*slots_per_edge);
        for (int i = ext_i + 1 ; z<nr.getSize(); z++, i++) {
            if (nr.hasEdges()) {
                if (ext_ptr) {
                    nr.setSparse(z, i, ext_ptr, ext_d);
                } else {
                    nr.setSparse(z, i, parent->getTransparentEdge(), ext_d);
                }
            } else {
                nr.setSparse(z, i, ext_d);
            }
            /*
          nr.i_ref(z) = i;
          nr.d_ref(z) = ext_d;
          if (nr.hasEdges()) {
            if (ext_ptr)
              nr.setEdgeRaw(z, ext_ptr);
            else
              nr.set_edgeval(z, parent->getTransparentEdge());
          }
          */
        }
      }
    }
  }

#ifdef DEBUG_DECODING
  printf("\n        temp:  ");
  nr.show(out, true);
  printf("\n");
  fflush(stdout);
#endif

#ifdef DEVELOPMENT_CODE
  if (!is_extensible || parent->isExtensibleLevel(nr.getLevel())) {
#ifdef VERIFY_DECODING
    if (!areDuplicates(addr, nr)) {
      FILE_output out(stdout);
      printf("Error decoding packed node: ");
      dumpInternalNode(out, addr, 0x03);
      printf("Unpacked to: ");
      nr.show(out, true);
      printf("\n");
      fflush(stdout);
    }
#else
    MEDDLY_DCASSERT(!nr.isTrim() || areDuplicates(addr, nr));
#endif
  }
#endif
}


unsigned MEDDLY::best_storage::hashNode(int level, node_address addr) const
{
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  hash_stream s;
  s.start(0);

  //
  // Hash any header portion
  //

  if (hashed_slots) {
    s.push(chunk + hashed_start, parent->hashedHeaderBytes());
  }

  //
  // Hash the node itself
  //

  const unsigned int raw_size = getRawSize(chunk);
  const unsigned int size = getSize(raw_size);
  const bool is_sparse = isSparse(raw_size);
  const bool is_pattern = isPattern(raw_size);
  const node_handle* down = chunk + down_start;


  if(is_pattern)
    {

    std::string pattern_from_index = generatePatternFromIndex(size);

    int trunc_pattern_size = 0;
    for(int i=0;i<MAX_PATTERN_LEN;i++)
      {
      if(pattern_from_index[i]!='t')
        trunc_pattern_size = i+1;
      }

    for(int i=0;i<trunc_pattern_size;i++)
      {
      if(pattern_from_index[i]!='t')
        {
        s.push(i, down[pattern_from_index[i] - 'A']);
        }
      }
    }else if (is_sparse) {
      //
      // Node is sparse
      //
      int nnz = size;
      const node_handle* index = down + nnz;
      const node_handle* edge = slots_per_edge ? (index + nnz) : 0;
      if (parent->areEdgeValuesHashed()) {
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
      const node_handle tv=parent->getTransparentNode();
      const node_handle* edge = slots_per_edge ? (down + size) : 0;
      if (parent->areEdgeValuesHashed()) {
        const int edge_bytes = bytesForSlots(slots_per_edge);
        for (unsigned i=0; i<size; i++) {
          if (down[i]!=tv) {
            s.push(i, down[i]);
            s.push(edge + i * slots_per_edge, edge_bytes);
          }
        } // for z
      } else {
        for (unsigned i=0; i<size; i++) {
          if (down[i]!=tv) {
            s.push(i, down[i]);
          }
        } // for z
      }
    }

  return s.finish();
}


int MEDDLY::best_storage
::getSingletonIndex(node_address addr, node_handle &down) const
{
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  const unsigned int raw_size = getRawSize(chunk);
  if (isExtensible(raw_size)) return -1;

  const unsigned int size = getSize(raw_size);
  const bool is_sparse = isSparse(raw_size);
  const bool is_pattern = isPattern(raw_size);

  if(is_pattern){
    std::string pattern = generatePatternFromIndex(size);
    int count_opaque=0;
    int singleton_idx=-1;
    for(int i=0;i<pattern.length();i++)
      {
      if(pattern[i]=='A')
        {
        if(count_opaque==0)
          {
          count_opaque++;
          singleton_idx = i;
          }
        else
          return -1;
        }else if(pattern[i]!='t')
          return -1;
      }

    down = chunk[down_start+0];
    return singleton_idx;
  }

  if (is_sparse) {
    //
    // sparse node --- easy
    //
    if (size != 1) return -1;
    down = chunk[0 + down_start];       // 0+ stops a compiler warning
    return chunk[down_start + size];    // size is number of nonzeroes
  }

  //
  // full node
  //
  const node_handle tv=parent->getTransparentNode();
  const node_handle* dnptr = chunk + down_start;
  for (unsigned i=0; i<size; i++) {
    if (tv==dnptr[i]) continue;
    if (i+1 != size) return -1;
    down = dnptr[i];
    return i;
  }
  return -1;
}


/// Extensible Index is the index of the last edge
int
MEDDLY::best_storage
::getExtensibleIndex(node_address addr) const
{
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  const unsigned int raw_size = getRawSize(chunk);
  if (!isExtensible(raw_size)) return -1;
  if (isPattern(raw_size)) return -1;

  int sz = getSize(raw_size);
  const node_handle* down = chunk + down_start;
  if (isSparse(raw_size)) {
    MEDDLY_DCASSERT(down[sz-1] != parent->getTransparentNode());
    const node_handle* index = down + sz;
    return index[sz-1];
  } else {
    MEDDLY_DCASSERT(down[sz-1] != parent->getTransparentNode());
    return sz-1;
  }
}


MEDDLY::node_handle
MEDDLY::best_storage
::getDownPtr(node_address addr, int i) const
{
  if (i<0) throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);

  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  const unsigned int raw_size = getRawSize(chunk);
  const unsigned size = getSize(raw_size);
  const node_handle* down = chunk + down_start;

  int z = i;

  if (isPattern(raw_size)){
    std::string pattern = generatePatternFromIndex(size);
    return
    (pattern[i]=='t') ?
    parent->getTransparentNode()
    : down[pattern[i]-'A'];
  }

  if (isSparse(raw_size)) {
    const node_handle* index = down + size;
    z =
    (isExtensible(raw_size) && i >= index[size-1])
    ? (size - 1)
    : findSparseIndex(i, index, size);
  } else {
    if (unsigned(i) >= size) {
      z = isExtensible(raw_size) ? size - 1 : -1;
    }
  }
  return
  (z < 0)
  ? parent->getTransparentNode()
  : down[z];
}


void
MEDDLY::best_storage
::getDownPtr(node_address addr, int i, float& ev, node_handle& dn) const
{
  if (i<0) throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);

  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(slots_per_edge>0);

  const unsigned int raw_size = getRawSize(chunk);
  const unsigned size = getSize(raw_size);
  const bool is_sparse = isSparse(raw_size);
  const node_handle* down = chunk + down_start;

  int z = i;

  if (is_sparse) {
    const node_handle* index = down + size;
    z =
    (isExtensible(raw_size) && i >= index[size-1])
    ? (size - 1)
    : findSparseIndex(i, index, size);
  } else {
    if (unsigned(i) >= size) {
      z = isExtensible(raw_size) ? size - 1 : -1;
    }
  }

  if (z < 0) {
    dn = 0;
    ev = 0;
  } else {
    const node_handle* edge = down + (is_sparse? 2*size: size);
    dn = down[z];
    ev = ((float*) (edge + z*slots_per_edge)) [0];
  }
}

void
MEDDLY::best_storage
::getDownPtr(node_address addr, int i, int& ev, node_handle& dn) const
{
  if (i<0) throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);

  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(slots_per_edge>0);

  const unsigned int raw_size = getRawSize(chunk);
  const unsigned size = getSize(raw_size);
  const bool is_sparse = isSparse(raw_size);
  const node_handle* down = chunk + down_start;

  int z = i;

  if (is_sparse) {
    const node_handle* index = down + size;
    z =
    (isExtensible(raw_size) && i >= index[size-1])
    ? (size - 1)
    : findSparseIndex(i, index, size);
  } else {
    if (unsigned(i) >= size) {
      z = isExtensible(raw_size) ? size - 1 : -1;
    }
  }

  if (z < 0) {
    dn = 0;
    ev = 0;
  } else {
    const node_handle* edge = down + (is_sparse? 2*size: size);
    dn = down[z];
    ev = ((int*) (edge + z*slots_per_edge)) [0];
  }
}

void MEDDLY::best_storage
::getDownPtr(node_address addr, int i, long& ev, node_handle& dn) const
{
  if (i<0) throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);

  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(slots_per_edge>0);

  const unsigned int raw_size = getRawSize(chunk);
  const unsigned size = getSize(raw_size);
  const bool is_sparse = isSparse(raw_size);
  const node_handle* down = chunk + down_start;

  int z = i;

  if (is_sparse) {
    const node_handle* index = down + size;
    z =
    (isExtensible(raw_size) && i >= index[size-1])
    ? (size - 1)
    : findSparseIndex(i, index, size);
  } else {
    if (unsigned(i) >= size) {
      z = isExtensible(raw_size) ? size - 1 : -1;
    }
  }

  if (z < 0) {
    dn = 0;
    ev = 0;
  } else {
    const node_handle* edge = down + (is_sparse? 2*size: size);
    dn = down[z];
    ev = ((long*) (edge + z*slots_per_edge)) [0];
  }
}


const void* MEDDLY::best_storage
::getUnhashedHeaderOf(node_address addr) const
{
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(unhashed_slots);

  return chunk + unhashed_start;
}


const void* MEDDLY::best_storage
::getHashedHeaderOf(node_address addr) const
{
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(hashed_slots);

  return chunk + hashed_start;
}


MEDDLY::node_handle MEDDLY::best_storage
::getNextOf(node_address addr) const
{
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  return chunk[next_slot];
}


void MEDDLY::best_storage
::setNextOf(node_address addr, node_handle n)
{
  node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(n>=0);
  chunk[next_slot] = n;
}


void MEDDLY::best_storage::updateData(node_handle* d)
{
  MEDDLY_DCASSERT(0);

  //
  // Required for old holeman class; eventually discard?
  //
}


void MEDDLY::best_storage::dumpInternalInfo(output &s) const
{
  s << "best_storage::dumpInternalInfo\n";
  // MM->dumpInternalInfo(s);
}


void MEDDLY::best_storage::dumpInternalTail(output &s) const
{
  s << "best_storage::dumpInternalTail\n";
  // MM->dumpInternalTail(s);
}


MEDDLY::node_address
MEDDLY::best_storage::firstNodeAddress() const
{
  return MM->getFirstAddress();
}


MEDDLY::node_address
MEDDLY::best_storage
::dumpInternalNode(output &s, node_address a, unsigned flags) const
{
  if (a<=0) return 0;

  MEDDLY_DCASSERT(MM);

  int awidth = digits(parent->getLastNode());

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
  node_handle* end = getChunkAddress(a);
  const node_handle* chunk = end;
  MEDDLY_DCASSERT(chunk);
  const bool show_node = flags & 0x01;

  if (show_node) {
    s.put(long(a), awidth);
    s << " : [";
  }
  const unsigned int raw_size = chunk[size_or_id_slot];
  const unsigned int size = getSize(raw_size);

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
  int dlen;

  if(isPattern(raw_size)){
    std::string pattern = generatePatternFromIndex(size);
    std::set<char> uniqnh;
    for (int i=0; i<MAX_PATTERN_LEN; i++) {
      if(pattern[i]!='t')
        uniqnh.insert(pattern[i]);
    }
    dlen = uniqnh.size();
  }

  const unsigned int dnlen = isPattern(raw_size)?dlen:size;


  if (show_node) {
    s.put(" down ");
    s.put(long(down[0]));
    for (unsigned int i=1; i<dnlen; i++) {
      s.put(", ");
      s.put(long(down[i]));
    }
    if (isExtensible(raw_size)) s.put(" (ext)");
    s.put(';');
  }
  end += dnlen;

  //
  // Indexes, if any
  //
  if (isSparse(raw_size)) {
    if (show_node) {
      s.put(" index ");
      const node_handle* index = down + size;
      s.put(long(index[0]));
      for (unsigned int i=1; i<dnlen; i++) {
        s.put(", ");
        s.put(long(index[i]));
      }
      if (isExtensible(raw_size)) s.put(" (ext)");
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
      for (unsigned int i=0; i<dnlen; i++) {
        if (i) s.put(", ");

        edge_value ev;
        ev.set(parent->getEdgeType(), end);
        ev.write(s);
        s << ")";
      } // for each edge value
      if (isExtensible(raw_size)) s.put(" (ext)");
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
    if (*end < -1024) {
      // Sanity check
      s.put('\n');
      MEDDLY_DCASSERT(0);
    }
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


MEDDLY::node_address MEDDLY::best_storage
::makeFullNode(node_handle p, int size, const unpacked_node &nb)
{

#ifdef DEBUG_ENCODING
  printf("\nBuilding full node with handle %ld\n", long(p));
#endif
  //
  // Determine amount of memory we need, and request it
  //

  size_t slots_req = slotsForNode(size, false);
  MEDDLY_DCASSERT(slots_req > 0);
  size_t slots_given = slots_req;
  node_address addr = MM->requestChunk(slots_given);
  if (0==addr) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  MEDDLY_DCASSERT(slots_given >= slots_req);

  node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  //
  // Set size
  //
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, (size_t)0, (size_t)size_or_id_slot, slots_given);
  chunk[size_or_id_slot] = getRawSize(size, false, nb.isExtensible(), false);

  //
  // Copy any extra header information
  //

  if (unhashed_slots) {
    nb.getUHdata(chunk + unhashed_start);
    // memcpy(chunk + unhashed_start, nb.UHptr(), nb.UHbytes());
  }
  if (hashed_slots) {
    nb.getHHdata(chunk + unhashed_start);
    // memcpy(chunk + hashed_start, nb.HHptr(), nb.HHbytes());
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
      const node_handle tv = parent->getTransparentNode();
      const edge_value &te = parent->getTransparentEdge();
      for (int i=0; i<size; i++) {
        down[i] = tv;
        te.get(parent->getEdgeType(), edge + i * edge_bytes);
      }
      for (int z=0; z<nb.getSize(); z++) {
        int i = nb.index(z);
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, size);
        down[i] = nb.down(z);
        nb.edgeval(z).get(parent->getEdgeType(), edge + i * edge_bytes);
      }
    } else {
      for (int i=0; i<size; i++) {
          down[i] = nb.down(i);
          nb.edgeval(i).get(parent->getEdgeType(), edge + i * edge_bytes);
      }
    }
  } else {
    //
    // No edge values
    //
    MEDDLY_DCASSERT(!nb.hasEdges());
    if (nb.isSparse()) {
      const node_handle tv = parent->getTransparentNode();
      for (int i=0; i<size; i++) {
        down[i] = tv;
      }
      for (int z=0; z<nb.getSize(); z++) {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, (int)nb.index(z), size);
        down[nb.index(z)] = nb.down(z);
      }
    } else {
      for (int i=0; i<size; i++) down[i] = nb.down(i);
    }
  }

  //
  // Deal with any padding and the tail
  //

  // int tail = down_start + size + slots_per_edge * size;
  long delta = slots_given - slots_req;
  MEDDLY_DCASSERT(delta>=0);
  MEDDLY_DCASSERT(delta<1024);    // Sanity check
  MEDDLY_DCASSERT(size_t(down_start + size + slots_per_edge * size + 1) == slots_req);
  chunk[slots_req-1] = -delta;  // Where we expect the node to end
  chunk[slots_given-1] = p;     // Where the node actually ends
                                // Note: if slots_req == slots_given, then the second statement
                                // overwrites chunk[slots_req-1], but that's exactly what we want.


  /*
   node_handle* tail = down + size + slots_per_edge * size;
   int delta = slots_given - slots_req;
   MEDDLY_DCASSERT(delta>=0);
   MEDDLY_DCASSERT( (tail-chunk)+1 == slots_req );
   tail[0] = -delta;
   tail[delta] = p;  // if delta = 0, we just overwrote but that's fine
   */

#ifdef DEBUG_ENCODING
  printf("\n    internal: ");
  FILE_output out(stdout);
  dumpInternalNode(out, addr, 0x03);
  fflush(stdout);
#endif
  return addr;
}


MEDDLY::node_address MEDDLY::best_storage
::makeSparseNode(node_handle p, int size, const unpacked_node &nb)
{
#ifdef DEBUG_ENCODING
  printf("\nBuilding sparse node with handle %ld\n", long(p));
#endif
  //
  // Determine amount of memory we need, and request it
  //

  size_t slots_req = slotsForNode(size, true);
  MEDDLY_DCASSERT(slots_req > 0);
  size_t slots_given = slots_req;
  node_address addr = MM->requestChunk(slots_given);
  if (0==addr) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  MEDDLY_DCASSERT(slots_given >= slots_req);

  node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);

  //
  // Set size
  //
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, (size_t)0, (size_t)size_or_id_slot, slots_given);
  chunk[size_or_id_slot] = getRawSize(size, true, nb.isExtensible(), false);

  //
  // Copy any extra header information
  //
  if (unhashed_slots) {
    nb.getUHdata(chunk + unhashed_start);
    // memcpy(chunk + unhashed_start, nb.UHptr(), nb.UHbytes());
  }
  if (hashed_slots) {
    nb.getHHdata(chunk + unhashed_start);
    // memcpy(chunk + hashed_start, nb.HHptr(), nb.HHbytes());
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
        down[z] = nb.down(z);
        index[z] = nb.index(z);
        nb.edgeval(z).get(parent->getEdgeType(), edge + z * edge_bytes);
      }
    } else {
      int z = 0;
      for (int i=0; i<nb.getSize(); i++) {
        if (parent->isTransparentEdge(nb.down(i), nb.edgeval(i))) continue;
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, z, size);
        down[z] = nb.down(i);
        index[z] = i;
        nb.edgeval(i).get(parent->getEdgeType(), edge + z * edge_bytes);
        z++;
      }
      MEDDLY_DCASSERT(size == z);
    }
  } else {
    //
    // No edge values
    //
    MEDDLY_DCASSERT(!nb.hasEdges());
    if (nb.isSparse()) {
      for (int z=0; z<size; z++) {
        down[z] = nb.down(z);
        index[z] = nb.index(z);
      }
    } else {
      int z = 0;
      const node_handle tv = parent->getTransparentNode();
      for (int i=0; i<nb.getSize(); i++) {
        if (nb.down(i)!=tv) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, z, size);
          down[z] = nb.down(i);
          index[z] = i;
          z++;
        }
      }
      MEDDLY_DCASSERT(size == z);
    }
  }

#ifdef DEVELOPMENT_CODE
  // check if the sparse node is sorted
  for (int z=1; z<size; z++) {
    MEDDLY_DCASSERT(index[z-1] < index[z]);
  }
#endif

  //
  // Deal with any padding and the tail
  //
  long delta = slots_given - slots_req;
  MEDDLY_DCASSERT(delta>=0);
  MEDDLY_DCASSERT(delta<1024);    // Sanity check
  MEDDLY_DCASSERT(size_t(down_start + 2*size + slots_per_edge * size + 1) == slots_req);
  chunk[slots_req-1] = -delta;  // Where we expect the node to end
  chunk[slots_given-1] = p;     // Where the node actually ends
                                // Note: if slots_req == slots_given, then the second statement
                                // overwrites chunk[slots_req-1], but that's exactly what we want.

  /*
   node_handle* tail = down + 2*size + slots_per_edge * size;
   MEDDLY_DCASSERT(delta>=0);
   tail[0] = -delta;
   tail[delta] = p;  // if delta = 0, we just overwrote but that's fine
   */

#ifdef DEBUG_ENCODING
  /*
   printf("\n        made: ");
   showNode(stdout, addr, true);
   */
  printf("\n    internal: ");
  FILE_output out(stdout);
  dumpInternalNode(out, addr, 0x03);
  fflush(stdout);
#endif
  return addr;
}


//
// Pattern Helpers
//

MEDDLY::node_address MEDDLY::best_storage
::makePatternNode(node_handle p, int size, const unpacked_node &nb, std::map<int,node_handle> nhmap)
{
#ifdef DEBUG_ENCODING
  printf("\nBuilding pattern node with handle %ld\n", long(p));
#endif
  //
  // Determine amount of memory we need, and request it
  //

  size_t slots_req = slotsForNode(size);
  MEDDLY_DCASSERT(slots_req > 0);
  size_t slots_given = slots_req;
  node_address addr = MM->requestChunk(slots_given);
  if (0==addr) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  MEDDLY_DCASSERT(slots_given >= slots_req);
  node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(!nb.isExtensible());

  //
  // Copy any extra header information
  //

  if (unhashed_slots) {
    nb.getUHdata(chunk + unhashed_start);
    // memcpy(chunk + unhashed_start, nb.UHptr(), nb.UHbytes());
  }
  if (hashed_slots) {
    nb.getHHdata(chunk + unhashed_start);
    // memcpy(chunk + hashed_start, nb.HHptr(), nb.HHbytes());
  }


  //
  // Set identifier
  //
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, (size_t)0, (size_t)size_or_id_slot, slots_given);
  int pattern_id = generateIndexFromNode(nb);
  chunk[size_or_id_slot] = getRawSize(pattern_id, false, false, true);
  //
  // Copy downward pointers, indexes, and edge values (if any)
  //

  node_handle* down = chunk + down_start;  // what is the point of copying these? Since node contains slots,


  if (slots_per_edge) {
    //
    // There's edge values
    //
    /** Not available for best_storage yet **/
  } else {
    //
    // No edge values
    //
    MEDDLY_DCASSERT(!nb.hasEdges());
    int i = 0;
    for (std::map<int,node_handle>::iterator it = nhmap.begin(); it!=nhmap.end(); it++) {
      // iterate through nhmap to retrieve unique pointers
      down[i] = it->second;
      i+=1;
    }

    MEDDLY_DCASSERT(i==size);

  }


  //
  // Deal with any padding and the tail
  //
  long delta = slots_given - slots_req;
  MEDDLY_DCASSERT(delta>=0);
  MEDDLY_DCASSERT(delta<1024);    // Sanity check
  MEDDLY_DCASSERT(size_t(down_start + size + slots_per_edge * size + 1) == slots_req);
  chunk[slots_req-1] = -delta;  // Where we expect the node to end
  chunk[slots_given-1] = p;     // Where the node actually ends
                                // Note: if slots_req == slots_given, then the second statement
                                // overwrites chunk[slots_req-1], but that's exactly what we want.

  /*
   node_handle* tail = down + 2*size + slots_per_edge * size;
   MEDDLY_DCASSERT(delta>=0);
   tail[0] = -delta;
   tail[delta] = p;  // if delta = 0, we just overwrote but that's fine
   */

#ifdef DEBUG_ENCODING
  /*
   printf("\n        made: ");
   showNode(stdout, addr, true);
   */
  printf("\n    internal: ");
  FILE_output out(stdout);
  dumpInternalNode(out, addr, 0x03);
  fflush(stdout);
#endif
  return addr;
}


std::string MEDDLY::best_storage
::specialReverse(std::string pattern) const
{
  std::string rev;
  int max = 0;
  std::map<char,char> rev_rep;
  char current='A';

  // Normal reversal
  for(int i = 0;i<pattern.length();i++)
    {
    rev+= pattern[pattern.length() - 1 -i];
    if((pattern[i]!='t')&&(pattern[i]>max))
      max = pattern[i];
    }


  for(int i = 0;i<pattern.length();i++)
    {
    if((rev[i]!='t')&& !(rev_rep.find(rev[i])!=rev_rep.end())) // If not transparent and not mapped already
      {
      rev_rep[rev[i]]=current;
      current++;

      }
    }

  //Replace opaque pointers accordingly
  for(int i = 0;i<pattern.length();i++)
    {
    if((rev[i]!='t'))
      {
      std::map<char,char>::iterator it = rev_rep.find(rev[i]);
      rev[i] = it->second;
      }
    }

  return rev;
}


std::string MEDDLY::best_storage
::generatePatternFromNode(const unpacked_node &nb) const
{

  std::string node_pattern;
  char seqid='A';
  std::map<node_handle,char> node_rep;
  const node_handle tv = parent->getTransparentNode();
  //build the pattern
  if(nb.isSparse())
    {
    if(slots_per_edge)
      {
      // not defined
      }else {
        //obtain the last truncated index
        node_handle full_truncated_size = nb.index(nb.getSize()-1) + 1;
        //
        // All entries before last_index are either given in down or are transparent
        //
        for(int i=0,z=0;i<full_truncated_size;i++)
          {
          if((z<nb.getSize()) && (nb.index(z) == i)) //the index of sparse node is same the current index being read
            {
            //assign a new character to the node_handle nb.down(z)
            //increment z
            std::map<node_handle,char>::iterator it = node_rep.find(nb.down(z));
            if(it !=node_rep.end()) // Already exists in the map
              {
              node_pattern +=it->second;
              }
            else // Add to the map
              {
              node_rep.insert (std::pair<node_handle,char>(nb.down(z),seqid));
              node_pattern +=seqid;
              seqid +=1;
              }
            z++;
            } else {
              // the current index is tv
              // assign 't' to this child
              node_pattern +='t';
            }
          }
      }
    } else {
      if(slots_per_edge)
        {
        // not defined
        } else {
          //full MT node

          //obtain the last truncated index
          node_handle full_truncated_size = 0;
          for (int i=0; i<nb.getSize(); i++) {
            if (nb.down(i) != tv) {
              full_truncated_size = i+1;
            }
          }

          //
          // All entries before last_index are either given in down or are transparent
          //
          for(int i=0;i<full_truncated_size;i++)
            {
            if(nb.down(i) !=tv) //the down pointer of full node not transparent
              {
              //assign a new character to the node_handle nb.down(z)
              //increment z
              std::map<node_handle,char>::iterator it = node_rep.find(nb.down(i));
              if(it !=node_rep.end()) // Already exists in the map
                {
                node_pattern +=it->second;
                }
              else // Add to the map
                {
                node_rep.insert (std::pair<node_handle,char>(nb.down(i),seqid));
                node_pattern +=seqid;
                seqid +=1;
                }
              } else {
                // the current index is tv
                // assign 't' to this child
                node_pattern +='t';
              }
            }
        }
    }

  std::string rev_pattern = specialReverse(node_pattern);

  return rev_pattern;
}


MEDDLY::node_handle MEDDLY::best_storage
::generateIndexFromNode(const unpacked_node &nb)
{


  std::string numPatt = generatePatternFromNode(nb);
  //All the patterns received are in full-truncated format
  int len = numPatt.length();

  int **arr  = (int**)malloc((len+1) * sizeof(int*));
  for( int i = 1; i <= len; i++)
    {
    arr[i] = (int*)malloc(len * sizeof(int));
    arr[i][0] = i+1;
    }


  /* arr[i][j] gives the A-weightage in (j+1)th position from right
   1 is the value of A-weightage in 1st position (rightmost)
   Every i has a fit-size of i+1
   -> If leftmost has lettervalue < i, then use arr[i][j-1] for next right
   -> If leftmost has lettervalue = i, then use arr[i+1][j-1] for next right
   */

  for( int j = 1; j <= len; j++)
    {
      for(int i = 1; i <= len; i++)
        {
        if(i + j <=len)
          {
          arr[i][j] = i*arr[i][j-1] + arr[i+1][j-1];
          }
        }
    }

  int idx = 0;
  int nxtG = 1;
  for( int i = len; i >= 1; i--)
    {
    int val = numPatt[len-i]=='t'?0:(int)(numPatt[len-i]) - 64;

    assert((val>=0)&&(val<=26)); // To ensure the pattern is entered in correct format

    idx += i-2>=0?val*arr[nxtG][i-2]:val*1;

    if(val == nxtG)
      nxtG++;
    }

  return idx;
}


//
//From the index build the pattern
//Using the unique pointers stored in down
//build the sequence of children

std::string MEDDLY::best_storage
::generatePatternFromIndex(node_handle index_of_pattern) const
{
  int len = MAX_PATTERN_LEN;
  int **arr  = (int**)malloc((len+1) * sizeof(int*));
  for( int i = 1; i <= len; i++)
    {
    arr[i] = (int*)malloc(len * sizeof(int));
    arr[i][0] = i+1;
    }


  /* arr[i][j] gives the A-weightage in (j+1)th position from right */
  /* 1 is the value of A-weightage in 1st position (rightmost) */
  /* Every i has a fit-size of i+1
   -> If leftmost has lettervalue < i, then use arr[i][j-1] for next right
   -> If leftmost has lettervalue = i, then use arr[i+1][j-1] for next right
   */
  for( int j = 1; j <= len; j++)
    {
      for(int i = 1; i <= len; i++)
        {
        if(i + j <=len)
          {
          arr[i][j] = i*arr[i][j-1] + arr[i+1][j-1];
          }
        }
    }

  assert(index_of_pattern <= arr[1][len-1]); // To ensure index is within the limits.

  std::string patt="";
  int nxtG = arr[1][len-2];
  int k=1;
  for( int i = 1; i <= len;i++)
    {
    int val = index_of_pattern / nxtG;

    if(val>=k)
      val = k;

    index_of_pattern = index_of_pattern - val*nxtG;

    char ch = val==0?'t':val + 64;
    patt += ch ;//std::to_string(ch);

    if(val == k)
      {
      k=k+1;
      nxtG = len-i-2>=0?arr[k][len-i-2]:1;
      }

    else
      nxtG = len-i-2>=0?arr[k][len-i-2]:1;
    }

  std::string rev_patt = specialReverse(patt);


  return rev_patt;

}




// ******************************************************************
// *                                                                *
// *                                                                *
// *                 best_storage_style methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::best_storage_style::best_storage_style(const char* n)
: node_storage_style(n)
{
}

MEDDLY::best_storage_style::~best_storage_style()
{
}

MEDDLY::node_storage* MEDDLY::best_storage_style
::createForForest(forest* f, const memory_manager_style* mst, memstats &ms) const
{
  //  return 0;
  return new best_storage(getName(), f, mst, ms);
}





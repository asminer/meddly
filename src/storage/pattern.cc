
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

#include "pattern.h"
#include <set>
#include <map>

#include "../hash_stream.h"
#define MAX_PATTERN_LEN 10

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
  class pattern_storage;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                     pattern_storage class                     *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Pattern node storage mechanism in a forest.
 Memory management is completely separated from the class,
 which makes this implementation a little different
 from the original.
 
 Details of node storage are left to the derived forests.
 However, every active node is stored in the following format.
 
 
 
 common   {  slot[0] : next pointer >=0 in unique table.
 header --{  slot[1] : Bit MSB..0  : size
 
 unhashed    {       : slots used for any extra information
 header    --{       : as needed on a forest-by-forest basis.
 (optional)  {       : Info does NOT affect node uniqueness.
 
 hashed      {       : slots used for any extra information
 header    --{       : as needed on a forest-by-forest basis.
 (optional)  {       : Info DOES affect node uniqueness.
 
 {       : Unique downward pointers.
 {       : If full storage, there are size pointers
 down -------{       : and entry i gives downward pointer i.
 {       : If sparse storage, there are -size pointers
 {       : and entry i gives a pointer but the index
 {       : corresponds to index[i], below.
 
 
 {       : A identifier for the pattern formed by the pointers.
 {       : Full or truncated-full storage, both have same identifier for trailing true values.
 identifier--{       : and entry i gives downward pointer i.
 {       : If sparse storage, there are -size pointers
 {       : and entry i gives a pointer but the index
 {       : .
 
 
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
class MEDDLY::pattern_storage : public node_storage {
public:
  pattern_storage(const char* n, expert_forest* f, const memory_manager_style* mst);
  virtual ~pattern_storage();
  
  // required interface
public:
  virtual void collectGarbage(bool shrink);
  virtual void reportStats(output &s, const char* pad, unsigned flags) const;
  
  virtual node_address makeNode(node_handle p, const unpacked_node &nb, 
                                node_storage_flags opt);
  
  virtual void unlinkDownAndRecycle(node_address addr);
  virtual void markDownPointers(node_address addr);
  
  virtual bool areDuplicates(node_address addr, const unpacked_node &nr) const;
  virtual void fillUnpacked(unpacked_node &nr, node_address addr, unpacked_node::storage_style) const;
  
  virtual unsigned hashNode(int level, node_address addr) const;
  virtual int getSingletonIndex(node_address addr, node_handle &down) const;
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
  ///   @param  sz      Number of unique downward pointers.
  inline int slotsForNode(int sz) const {
    int node_slots = sz + 1; // 1 extra to store the identifier.
    return extra_slots + unhashed_slots + hashed_slots + node_slots;
  }
  
  
  inline node_handle* getChunkAddress(node_address addr) const {
    MEDDLY_DCASSERT(MM);
    return (node_handle*) MM->getChunkAddress(addr);
  }
  
  
  virtual bool isExtensible(node_address addr) const {
    //TBD
    return false;
  }
  
private:
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
  static const int header_slots = next_slot + 1; // header starts after next_slot
  
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
  char identifier_start;
  char down_start;
  char slots_per_edge;
};


// ******************************************************************
// *                                                                *
// *                                                                *
// *                    pattern_storage methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::pattern_storage
::pattern_storage(const char* n, expert_forest* f, const memory_manager_style* mst)
: node_storage(n, f)
{
  
  MM = mst->initManager(sizeof(node_handle), slotsForNode(0), f->changeMemStats());
  
  unhashed_start = header_slots;
  unhashed_slots = slotsForBytes(f->unhashedHeaderBytes());
  hashed_start = unhashed_start + unhashed_slots;
  hashed_slots = slotsForBytes(f->hashedHeaderBytes());
  identifier_start = hashed_start + hashed_slots;
  down_start = identifier_start + 1;
  slots_per_edge = slotsForBytes(f->edgeBytes()); // 0 if not a part of edge-valued forest
}

MEDDLY::pattern_storage::~pattern_storage()
{
  // TBD - special steps to recycle all nodes?
  delete MM;
}

void MEDDLY::pattern_storage::collectGarbage(bool shrink)
{
  // TBD
}

void MEDDLY::pattern_storage::reportStats(output &s, const char* pad, 
                                          unsigned flags) const
{
  
  static unsigned STORAGE = 
  expert_forest::STORAGE_STATS | expert_forest::STORAGE_DETAILED;
  
  if (flags & STORAGE) {
    // s << pad << "Stats for " << getStyleName() << "\n";
    
  }
  
  static unsigned HOLEMAN =
  expert_forest::HOLE_MANAGER_STATS | expert_forest::HOLE_MANAGER_DETAILED;
  
  if (flags & HOLEMAN) {
    MM->reportStats(s, pad, 
                    flags & expert_forest::HUMAN_READABLE_MEMORY,
                    flags & expert_forest::HOLE_MANAGER_DETAILED);
  }
  
  
  /*
   #ifdef DEVELOPMENT_CODE
   verifyStats();
   #endif
   */
}

MEDDLY::node_address MEDDLY::pattern_storage
::makeNode(node_handle p, const unpacked_node &nb, node_storage_flags opt)
{
  
#ifdef DEBUG_ENCODING
  printf("pattern_storage making node\n        temp:  ");
  FILE_output out(stdout);
  nb.show(out, true);
#endif
  //
  // Determine unique non-zero pointers : uniqnnzs
  //
  int uniqnnzs = 0;
  std::map<int,node_handle> nhmap;
  std::set<node_handle> nhset;
  
  if (nb.isSparse()) {
    //
    // nb is sparse
    //
    MEDDLY_DCASSERT(nb.getNNZs()<MAX_PATTERN_LEN);
    int nhseq = 1;
    for (int i=0; i<nb.getNNZs(); i++) 
      {
      int old_size = nhset.size();
      nhset.insert(nb.d(i));
      if(nhset.size()>old_size) // if set size increases, then this node_handle is unique
        {
        nhmap.insert(std::pair<int,node_handle>(nhseq,nb.d(i))); // set contains only unique pointers
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
      const node_handle tv = getParent()->getTransparentNode();
      
      int nhseq = 1;
      for (int i=0; i<nb.getSize(); i++) {
        if (nb.d(i) != tv) {
          int old_size = nhset.size();
          nhset.insert(nb.d(i));
          if(nhset.size()>old_size) // if set size increases, then this node_handle is unique
            {
            nhmap.insert(std::pair<int,node_handle>(nhseq,nb.d(i))); // insert only the unique non-tv pointers
            nhseq+=1;
            }
        }
      }
      uniqnnzs = nhmap.size();
    }
  }
  
  node_address addr = 0;
  addr = makePatternNode(p, uniqnnzs, nb, nhmap);
  
  
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



void MEDDLY::pattern_storage::unlinkDownAndRecycle(node_address addr)
{
  
#ifdef MEMORY_TRACE
  printf("recycling node at address %ld\n", addr);
  FILE_output out(stdout);
  dumpInternalNode(out, addr, 0x03);
#endif
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  
  
  //
  // Unlink down pointers
  //
  const node_handle* index = chunk + identifier_start;
  const node_handle* down = chunk + down_start;
  
  std::string pattern_from_index = generatePatternFromIndex(index[0]);
  
  
  int trunc_pattern_size = 0;
  std::set<char> uniqnh;
  for(int i=0;i<MAX_PATTERN_LEN;i++)
    {
    if(pattern_from_index[i]!='t')
      {
      trunc_pattern_size = i+1;
      uniqnh.insert(pattern_from_index[i]);
      }
    }
  const unsigned int uniqnnzs = uniqnh.size();
  
  for (int i=0; i<trunc_pattern_size; i++) {
    if(pattern_from_index[i]!='t')
      getParent()->unlinkNode(down[pattern_from_index[i]-'A']);
  }
  //
  // Unlink the index for pattern
  //
  //getParent()->unlinkNode(index[0]);
  
  // Can this move after unlinking?
  MEDDLY_DCASSERT(MM->getChunkAddress(addr) == chunk);
  
  //
  // Determine number of slots in this node
  //
  size_t actual_slots = slotsForNode(uniqnnzs);
  if (chunk[actual_slots-1] < 0) {
    // padding
    actual_slots += (-chunk[actual_slots-1]);
  }
  
  //
  // Recycle
  //
  MM->recycleChunk(addr, actual_slots);
}


void MEDDLY::pattern_storage::markDownPointers(node_address addr)
{
#ifdef DEBUG_MARK_SWEEP
  printf("marking children at address %ld\n", addr);
  FILE_output out(stdout);
  dumpInternalNode(out, addr, 0x03);
#endif
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  
  //
  // Mark down pointers
  //
  const node_handle* index = chunk + identifier_start;
  const node_handle* down = chunk + down_start;
  
  std::string pattern_from_index = generatePatternFromIndex(index[0]);
  
  
  int trunc_pattern_size = 0;
  std::set<char> uniqnh;
  for(int i=0;i<MAX_PATTERN_LEN;i++)
    {
    if(pattern_from_index[i]!='t')
      {
      trunc_pattern_size = i+1;
      uniqnh.insert(pattern_from_index[i]);
      }
    }
  const unsigned int uniqnnzs = uniqnh.size();
  
  for (int i=0; i<trunc_pattern_size; i++) {
    if(pattern_from_index[i]!='t')
      getParent()->markNode(down[pattern_from_index[i]-'A']);
  }
}


bool MEDDLY::pattern_storage
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
  
  //const unsigned int raw_size = getRawSize(chunk);
  
  if (n.isExtensible()) return false;
  
  const node_handle* index = chunk + identifier_start;
  const node_handle* down = chunk + down_start;
  
  
  //retrieve the actual pattern formed by a node and reppresented by an index_no.
  std::string pattern_from_node =  specialReverse(generatePatternFromNode(n));
  std::string pattern_from_index = generatePatternFromIndex(index[0]);
  
  int trunc_pattern_size = 0;
  for(int i=0;i<MAX_PATTERN_LEN;i++)
    {
    if(pattern_from_index[i]!='t')
      trunc_pattern_size = i+1;
    }
  const node_handle tv = getParent()->getTransparentNode();
  if(n.isFull())
    {
    //
    // n is full
    //
    if (trunc_pattern_size > unsigned(n.getSize())) return false;
    // check down
    unsigned int i;
    for (i=0; i<trunc_pattern_size; i++) {
      if ( ( (n.d(i) != tv) && (down[pattern_from_index[i]-'A'] != n.d(i) ) )
          || ( (n.d(i) == tv) && (pattern_from_index[i]!='t') )
          ) return false;
    }
    for (; i<unsigned(n.getSize()); i++) {
      if (n.d(i)!=tv) return false;
    }
    //must be equal
    return true;
    }else {
      //
      // n is sparse
      //
      int i = 0;
      // check down
      for (int z=0; z<n.getNNZs(); z++) {
        if (unsigned(n.i(z)) >= trunc_pattern_size) return false;
        // loop - skipped edges must be transparent
        for (; i<n.i(z); i++) {
          if (pattern_from_index[i]!='t') return false;
        }
        if (n.d(z) != down[pattern_from_index[i]-'A']) return false;
        i++;
      } // for z
      if (unsigned(i)<trunc_pattern_size) return false; // there WILL be a non-zero down
      
      //must be equal
      return true;
    }
}


void MEDDLY::pattern_storage
::fillUnpacked(unpacked_node &nr, node_address addr, unpacked_node::storage_style st2) const
{
  
#ifdef DEBUG_DECODING
  FILE_output out(stdout);
  printf("pattern_storage filling reader\n    internal: ");
  dumpInternalNode(out, addr, 0x03);
#endif
  
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  
  MEDDLY_DCASSERT(nr.hasEdges() == (slots_per_edge>0) );
  
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
  
  const node_handle* index = chunk + identifier_start; // Verify this
  const node_handle* down = chunk + down_start;
  
  /*
   Set the unpacked node storage style based on settings
   */
  
  switch (st2) {
    case unpacked_node::FULL_NODE:     nr.bind_as_full(true);    break;
    case unpacked_node::SPARSE_NODE:   nr.bind_as_full(false);   break;
    case unpacked_node::AS_STORED:     nr.bind_as_full(true);    break; // Since is_sparse makes no sense here, lets go with full      
    default:            assert(0);
  };
  
  
  // Make sure that when an extensible node is unpacked that the trailing edges
  // are filled correctly (regardless of the storage scheme of the unpacked node)
  
  const node_handle tv = getParent()->getTransparentNode();
  
  MEDDLY_DCASSERT(getParent()->isExtensibleLevel(nr.getLevel()) == false);
  nr.markAsNotExtensible();
  
  //pattern node is stored with size of unique pointers + the unique pointers + identifier for the pattern thus formed
  
  if(nr.isFull()) {
    //
    //Copying into full Node
    //
    
    //Obtain the pattern formed by child pointers
    std::string node_pattern = generatePatternFromIndex(index[0]);
    
    //Obtain size for nr
    int trunc_pattern_size = 0; 
    for(int i=0;i<MAX_PATTERN_LEN;i++)
      {
      if(node_pattern[i]!='t')
        trunc_pattern_size=i+1;
      }
    
    //populate nr.down
    
    for(int i=0;i<nr.getSize();i++)
      {
      nr.d_ref(i) = tv;
      }
    
    for(int i=0;i<trunc_pattern_size;i++)
      {
      if(node_pattern[i] != 't') //If not transparent
        nr.d_ref(i) = down[node_pattern[i] - 'A'];
      }
    
  } else {
    
    //
    //Copying into sparse node
    //
    
    //Obtain the pattern formed by child pointers
    std::string node_pattern = generatePatternFromIndex(index[0]);
    
    //Obtain size for nr
    int trunc_pattern_size = 0;
    for(int i=0;i<MAX_PATTERN_LEN;i++)
      {
      if(node_pattern[i]!='t')
        {
        trunc_pattern_size = i+1;
        }
      }
    
    //populate nr.d() and nr.i()
    int z = 0;
    for(int i=0;i<trunc_pattern_size;i++)
      {
      if(node_pattern[i]!='t') //If transparent
        {
        nr.i_ref(z) = i;
        nr.d_ref(z) = down[node_pattern[i] - 'A'];
        z++;
        }
      }
    if (nr.isExtensible() == false) { nr.shrinkSparse(z);}
    
    
  }
  
#ifdef DEBUG_DECODING
  printf("\n  temp:  ");
  nr.show(out, true);
  printf("\n");
  fflush(stdout);
#endif
  
#ifdef DEVELOPMENT_CODE
  if (getParent()->isExtensibleLevel(nr.getLevel())) {
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


unsigned MEDDLY::pattern_storage::hashNode(int level, node_address addr) const
{
  
  const node_handle* chunk = getChunkAddress(addr);
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
  
  
  const node_handle* index = chunk + identifier_start;
  const node_handle* down = chunk + down_start;
  
  std::string pattern_from_index = generatePatternFromIndex(index[0]);
  
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
  
  
  return s.finish();
}


int MEDDLY::pattern_storage
::getSingletonIndex(node_address addr, node_handle &down) const
{
  
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  
  const node_handle* index = chunk + identifier_start;
  
  std::string pattern = generatePatternFromIndex(index[0]);
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


/// Extensible Index is the index of the last edge
int
MEDDLY::pattern_storage
::getExtensibleIndex(node_address addr) const
{
  //T.B.D , currently does not handle extensible
  return -1;
}


MEDDLY::node_handle 
MEDDLY::pattern_storage
::getDownPtr(node_address addr, int i) const
{
  
  if (i<0) throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);
  
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  
  const node_handle* index = chunk + identifier_start;
  const node_handle* down = chunk + down_start;
  
  std::string pattern = generatePatternFromIndex(index[0]);
  return 
  (pattern[i]=='t') ? 
  getParent()->getTransparentNode()
  : down[pattern[i]-'A'];
  
}


void
MEDDLY::pattern_storage
::getDownPtr(node_address addr, int i, float& ev, node_handle& dn) const
{
  
  //T.B.D
}

void
MEDDLY::pattern_storage
::getDownPtr(node_address addr, int i, int& ev, node_handle& dn) const
{
  
  //T.B.D
}

void MEDDLY::pattern_storage
::getDownPtr(node_address addr, int i, long& ev, node_handle& dn) const
{
  
  //T.B.D
}


const void* MEDDLY::pattern_storage
::getUnhashedHeaderOf(node_address addr) const
{
  
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(unhashed_slots);
  
  return chunk + unhashed_start;
}


const void* MEDDLY::pattern_storage
::getHashedHeaderOf(node_address addr) const
{
  
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(hashed_slots);
  
  return chunk + hashed_start;
}


MEDDLY::node_handle MEDDLY::pattern_storage
::getNextOf(node_address addr) const
{
  
  const node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  
  return chunk[next_slot];
}


void MEDDLY::pattern_storage
::setNextOf(node_address addr, node_handle n)
{
  node_handle* chunk = getChunkAddress(addr);
  MEDDLY_DCASSERT(chunk);
  MEDDLY_DCASSERT(n>=0);
  chunk[next_slot] = n;
}


void MEDDLY::pattern_storage::updateData(node_handle* d)
{
  
  MEDDLY_DCASSERT(0);
  
  //
  // Required for old holeman class; eventually discard?
  //
}


void MEDDLY::pattern_storage::dumpInternalInfo(output &s) const
{
  
  s << "pattern_storage::dumpInternalInfo\n";
  // MM->dumpInternalInfo(s);
}


void MEDDLY::pattern_storage::dumpInternalTail(output &s) const
{
  
  s << "pattern_storage::dumpInternalTail\n";
  // MM->dumpInternalTail(s);
}


MEDDLY::node_address 
MEDDLY::pattern_storage::firstNodeAddress() const
{
  
  return MM->getFirstAddress();
}


MEDDLY::node_address 
MEDDLY::pattern_storage
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
  node_handle* end = getChunkAddress(a);
  const node_handle* chunk = end;
  MEDDLY_DCASSERT(chunk);
  const bool show_node = flags & 0x01;
  
  if (show_node) {
    s.put(long(a), awidth);
    s << " : [";
  }
  
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
  // Pattern Identifier
  //
  end += identifier_start;
  const node_handle* index = end;
  
  std::string pattern = generatePatternFromIndex(index[0]);
  std::set<char> uniqnh;
  for (int i=0; i<MAX_PATTERN_LEN; i++) { 
    if(pattern[i]!='t')
      uniqnh.insert(pattern[i]);
  }
  const unsigned int dnlen = uniqnh.size();
  
  if (show_node) {
    s.put(" index ");
    s.put(long(index[0]));
    s.put(';');
  }
  
  //
  // Down pointers
  //
  
  end+=1;
  const node_handle* down = end;
  
  
  if (show_node) {
    s.put(" down ");
    s.put(long(down[0]));
    for (int i=1; i<dnlen; i++) {
      s.put(", ");
      s.put(long(down[i]));
    }
    //if (isExtensible(raw_size)) s.put(" (ext)");
    s.put(';');
  }
  end += dnlen;
  
  
  
  
  //
  // Edge values, if any
  //
  if (slots_per_edge) {
    //T.B.D
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

MEDDLY::node_address MEDDLY::pattern_storage
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
  node_handle* index = chunk + identifier_start;
  
  //
  // Set identifier
  //
  //MEDDLY_CHECK_RANGE(0, identifier_slot, slots_given);
  index[0] = generateIndexFromNode(nb);
  
  node_handle* down = chunk + down_start;  // what is the point of copying these? Since node contains slots,
  
  
  if (slots_per_edge) { 
    //
    // There's edge values
    //
    /** Not available for pattern_storage yet **/
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
  MEDDLY_DCASSERT(size_t(identifier_start + size + 1 + slots_per_edge * size + 1) == slots_req);
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


std::string MEDDLY::pattern_storage
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


std::string MEDDLY::pattern_storage
::generatePatternFromNode(const unpacked_node &nb) const
{
  
  std::string node_pattern;
  char seqid='A';
  std::map<node_handle,char> node_rep; 
  const node_handle tv = getParent()->getTransparentNode();
  //build the pattern
  if(nb.isSparse())
    {
    if(slots_per_edge)
      {
      // not defined
      }else {
        //obtain the last truncated index
        node_handle full_truncated_size = nb.i(nb.getNNZs()-1) + 1;
        //
        // All entries before last_index are either given in down or are transparent
        //
        for(int i=0,z=0;i<full_truncated_size;i++)
          {
          if((z<nb.getNNZs()) && (nb.i(z) == i)) //the index of sparse node is same the current index being read
            {
            //assign a new character to the node_handle nb.d(z)
            //increment z
            std::map<node_handle,char>::iterator it = node_rep.find(nb.d(z));
            if(it !=node_rep.end()) // Already exists in the map
              {
              node_pattern +=it->second;
              }
            else // Add to the map
              {
              node_rep.insert (std::pair<node_handle,char>(nb.d(z),seqid));
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
            if (nb.d(i) != tv) {
              full_truncated_size = i+1;
            }
          }
          
          //
          // All entries before last_index are either given in down or are transparent
          //
          for(int i=0;i<full_truncated_size;i++)
            {
            if(nb.d(i) !=tv) //the down pointer of full node not transparent
              {
              //assign a new character to the node_handle nb.d(z)
              //increment z
              std::map<node_handle,char>::iterator it = node_rep.find(nb.d(i));
              if(it !=node_rep.end()) // Already exists in the map
                {
                node_pattern +=it->second;
                }
              else // Add to the map
                {
                node_rep.insert (std::pair<node_handle,char>(nb.d(i),seqid));
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


MEDDLY::node_handle MEDDLY::pattern_storage
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

std::string MEDDLY::pattern_storage
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
// *                 pattern_storage_style methods                 *
// *                                                                *
// *                                                                *
// ******************************************************************


MEDDLY::pattern_storage_style::pattern_storage_style(const char* n)
: node_storage_style(n)
{
}

MEDDLY::pattern_storage_style::~pattern_storage_style()
{
}

MEDDLY::node_storage* MEDDLY::pattern_storage_style
::createForForest(expert_forest* f, const memory_manager_style* mst) const
{
  //  return 0;
  return new pattern_storage(getName(), f, mst);
}




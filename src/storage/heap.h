
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


#ifndef HEAP_H
#define HEAP_H

namespace MEDDLY {

// ******************************************************************

/** Things we want shared by all heaps
*/
class heapCommon {
  protected:
    static long pathPtr[64];
    static bool pathRight[64];
    static int pathlen;
};

long heapCommon::pathPtr[64];
bool heapCommon::pathRight[64];
int heapCommon::pathlen;


// ******************************************************************

/** Generic heap data structure as a template.

    Your basic heap data structure from undergrad days,
    EXCEPT we do NOT store the heap in an array, but instead
    allow the heap nodes to be spread throughout memory.
    (Typical application: chunks of memory that we need to
    keep track of for various reasons in a "priority queue".)

    Details of the node storage are abstracted away using
    the heapManager (the template parameter).  Nodes are
    referred to using "handles", which are non-zero.
    Each heapManager must provide the following functions:

        unsigned long ID(long handle) const
            Get the ID number of the given node.
            This specifies the position in the heap.
            The root node has ID 1.  In general, 
            a node with ID i has left child with ID 2i, 
            right child with ID 2i+1, and parent with ID i/2.


        long key(long handle) const
            Get the key for the given node.
            This is what we base the heap property on.


      
        void reset(long handle, unsigned long ID, long left, long right)
            Set the ID and children nodes for the given node.

        long& left(long handle)
            Reference to the left child of the given node
    
        long& right(long handle)
            Reference to the right child of the given node
*/

template <class HM>
class myheap : heapCommon {
  public:
    myheap(HM& heapMan);

    /// Return the number of elements in the heap
    inline unsigned long getSize() const { return size; }

    /// Dump the heap contents to the given stream, for debugging. 
    void dumpHeap(FILE* s) const;

    /// Add an item to the heap
    void addToHeap(long handle);

    /// Remove an item from the heap, by its ID number
    long removeByID(unsigned long ID);


  protected:
    void findPathToID(unsigned long ID) const;
    long downHeap(unsigned long ID, long left, long right, long replace);

  private:
    HM& heapMan;
    unsigned long size;
    long root;
};


// ------------------------------------------------------------

template <class HM>
myheap<HM>::myheap(HM& hm) : heapMan(hm)
{
  size = 0;
  root = 0;
}

// ------------------------------------------------------------

template <class HM>
void myheap<HM>::dumpHeap(FILE* s) const
{
  int bit = 1;
  unsigned long two2b = 0x1;
  for (unsigned long I=1; I<=size; I++) {
    if (I>=two2b) {
      fprintf(s, "\nLevel %2d: ", bit);
      bit++;
      two2b <<= 1;
    }
    findPathToID(I);
    assert(I == heapMan.ID(pathPtr[pathlen]));
    fprintf(s, "%ld ", heapMan.key(pathPtr[pathlen]));
  }
  fprintf(s, "\n");
}

// ------------------------------------------------------------

template <class HM>
void myheap<HM>::addToHeap(long ptr)
{
#ifdef DEBUG_ADD
  fprintf(stderr, "addToHeap(%ld)\n", ptr);
#endif

  unsigned long newID = ++size;
  findPathToID(newID);
  if (0==pathlen) {
    root = ptr;
    heapMan.reset(ptr, newID, 0, 0);
    return;
  }

#ifdef DEBUG_ADD
  fprintf(stderr, "addToHeap Path: ");
  for (int i=0; i<pathlen; i++) {
    fprintf(stderr, "%ld (%c) ", pathPtr[i], pathRight[i]?'r':'l');
  }
  fprintf(stderr, "\n");
#endif

  // determine position in path to insert node
  int ins;
  long kh = heapMan.key(ptr);
  unsigned long ID = 0x1;
  for (ins=0; ins<pathlen; ins++) {
    if (kh < heapMan.key(pathPtr[ins])) break;
    ID *= 2;
    if (pathRight[ins]) ID++;
  }

#ifdef DEBUG_ADD
  fprintf(stderr, "addToHeap path position is %d ID %lu\n", ins, ID);
#endif

  // adjust incoming pointer
  if (ins) {
    if (pathRight[ins-1])   heapMan.right(pathPtr[ins-1]) = ptr;
    else                    heapMan.left(pathPtr[ins-1]) = ptr;
  } else {
    root = ptr;
  }

  // Adjust nodes along the path, downward
  for (; ins<pathlen; ins++) {
    long child = pathPtr[ins];
    if (pathRight[ins]) {
      heapMan.reset(ptr, ID, heapMan.left(child), child);
      ID *= 2;
      ID++;
    } else {
      heapMan.reset(ptr, ID, child, heapMan.right(child));
      ID *= 2;
    }
    ptr = child;
  }
  heapMan.reset(ptr, ID, 0, 0);
  MEDDLY_DCASSERT(ID == newID);
}

// ------------------------------------------------------------

template <class HM>
long myheap<HM>::removeByID(unsigned long ID)
{
#ifdef DEBUG_REMOVE
  fprintf(stderr, "removeByID(%lu)\n", ID);
#endif
  if (size<1) return -1;

  // remove the rightmost leaf
  findPathToID(size);
  size--;
  long replace = pathPtr[pathlen];
  if (pathlen) {
    long parent = pathPtr[pathlen-1];
    if (pathRight[pathlen-1]) {
      assert(replace == heapMan.right(parent));
      heapMan.right(parent) = 0;
    } else {
      assert(replace == heapMan.left(parent));
      heapMan.left(parent) = 0;
    }
  } else {
#ifdef DEBUG_REMOVE
    fprintf(stderr, "removeByID: removing root %ld\n", replace);
    assert(replace == root);
#endif
    root = 0;
    return replace;
  }
#ifdef DEBUG_REMOVE
  fprintf(stderr, "removeByID: got replace %ld\n", replace);
#endif

  // find spot we really want to remove
  findPathToID(ID);
  long remove = pathPtr[pathlen];

#ifdef DEBUG_REMOVE
  fprintf(stderr, "removeByID Path: ");
  for (int i=0; i<=pathlen; i++) {
    fprintf(stderr, "%ld (%c) ", pathPtr[i], pathRight[i]?'r':'l');
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "removeByID: removing %ld\n", remove);
#endif

  long chl = heapMan.left(remove);
  long chr = heapMan.right(remove);
  if (pathlen) {
    if (pathRight[pathlen-1]) heapMan.right(pathPtr[pathlen-1]) = downHeap(ID, chl, chr, replace);
    else                      heapMan.left(pathPtr[pathlen-1]) = downHeap(ID, chl, chr, replace);
  } else {
    root = downHeap(ID, chl, chr, replace);
  }

  return remove;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

template <class HM>
void myheap<HM>::findPathToID(unsigned long ID) const
{
  unsigned long two2b = 0x1;
  int bit = 0;
  for (; two2b <= ID; two2b <<=1) { bit++; }
  two2b >>=1; // go up one level
  long n = root; 
  pathlen = 0;
  for (;;) {
    two2b >>=1;
    pathPtr[pathlen] = n;
    pathRight[pathlen] = ID & two2b;
    if (0==--bit) return;
    n = (pathRight[pathlen]) ? heapMan.right(n) : heapMan.left(n);
    pathlen++;
  }
}

// ----------------------------------------------------------------------

template <class HM>
long myheap<HM>::downHeap(unsigned long ID, long left, long right, long replace)
{
  assert(replace);
  long kp = heapMan.key(replace);
  long kl = left ? heapMan.key(left) : kp+1;    // sane value for null pointer
  long kr = right ? heapMan.key(right) : kp+1;  // sane value for null pointer

  if (kl < kr) {
    if (kp < kl) {
        // kp is smallest
        heapMan.reset(replace, ID, left, right);
#ifdef DEBUG_REMOVE
        fprintf(stderr, "downHeap(%lu, %ld, %ld, %ld) = %ld\n", ID, left, right, replace, replace);
#endif
        return replace;
    } else {
        // kl is smallest
        long newleft = downHeap(2*ID, heapMan.left(left), heapMan.right(left), replace);
        heapMan.reset(left, ID, newleft, right);
#ifdef DEBUG_REMOVE
        fprintf(stderr, "downHeap(%lu, %ld, %ld, %ld) = %ld\n", ID, left, right, replace, left);
#endif
        return left;
      }
  } else {
      if (kp < kr) {
        // kp is smallest
        heapMan.reset(replace, ID, left, right);
#ifdef DEBUG_REMOVE
        fprintf(stderr, "downHeap(%lu, %ld, %ld, %ld) = %ld\n", ID, left, right, replace, replace);
#endif
        return replace;
      } else {
        // kr is smallest
        long newright = downHeap(2*ID+1, heapMan.left(right), heapMan.right(right), replace);
        heapMan.reset(right, ID, left, newright);
#ifdef DEBUG_REMOVE
        fprintf(stderr, "downHeap(%lu, %ld, %ld, %ld) = %ld\n", ID, left, right, replace, right);
#endif
        return right;
      }
  }
}



// ******************************************************************

} // namespace MEDDLY

#endif

// $Id: hm_grid.cc 436 2013-07-01 21:43:08Z asminer $

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



#include "hm_heap.h"


#define MERGE_AND_SPLIT_HOLES
// #define DEBUG_COMPACTION
// #define DEBUG_SLOW
// #define MEMORY_TRACE
// #define DEEP_MEMORY_TRACE

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        hm_heap  methods                        *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::hm_heap::hm_heap(node_storage* p) : holeman(6, p)
{
  /*
  max_request = 0;
  large_holes = 0;
  holes_top = 0;
  holes_bottom = 0;
  */
}

// ******************************************************************

MEDDLY::hm_heap::~hm_heap()
{
}

// ******************************************************************
// *                                                                *
// *                        private  helpers                        *
// *                                                                *
// ******************************************************************

void MEDDLY::hm_heap::insertHole(node_handle p_offset)
{
  // TBD
}

// ******************************************************************

void MEDDLY::hm_heap::removeHole(node_handle p_offset)
{
  // TBD
}

// ******************************************************************

void MEDDLY::hm_heap
::addToHeap(node_handle& root, node_handle& size, node_handle p)
{
  MEDDLY_DCASSERT(p);
  // size remains negative
  size--;
  node_handle IDp = -size;
  node_handle pp = takePathToID(root, IDp);
  if (!pp) {
    // empty heap
    MEDDLY_DCASSERT(1==IDp);
    makeHeapNode(p, IDp, 0, 0, 0);
    root = p;
    return;
  }
  MEDDLY_DCASSERT(ID(pp) == IDp);
  makeHeapNode(p, IDp, pp, 0, 0);
  if (IDp % 2) {
    Right(pp) = p;
  } else {
    Left(pp) = p;
  }
  // p is in position; now perform an "upheap" operation
  node_handle* node = holeOf(p);
  while (pp) {
    if (pp < p) return; // heap condition is satisfied; stop

    // "swap" p with pp
    node_handle* pptr = holeOf(pp);

    // sanity checks
    MEDDLY_DCASSERT(pptr[ID_index] > 0);
    MEDDLY_DCASSERT(node[ID_index] > 0);
    MEDDLY_DCASSERT(node[ID_index] / 2 == pptr[ID_index]);
    MEDDLY_DCASSERT(node[parent_index] == pp);
    // determine which child
    if (node[ID_index] % 2) {
      // node is a right child
      MEDDLY_DCASSERT(pptr[right_index] == p);
      pptr[right_index] = node[right_index];
      node[right_index] = pp;
      SWAP(pptr[left_index], node[left_index]);
    } else {
      // node is a left child
      MEDDLY_DCASSERT(pptr[left_index] == p);
      pptr[left_index] = node[left_index];
      node[left_index] = pp;
      SWAP(pptr[right_index], node[right_index]);
    }
    // update the rest of the node
    pp = node[parent_index] = pptr[parent_index];
    pptr[parent_index] = p;
    pptr[ID_index] = node[ID_index];
    node[ID_index] /= 2;
  }
  // we made it all the way to the root
  root = p;
}

// ******************************************************************

void MEDDLY::hm_heap
::removeFromHeap(node_handle& root, node_handle& size, node_handle p)
{
  MEDDLY_DCASSERT(p);
  MEDDLY_DCASSERT(p == takePathToID(root, ID(p)));

  // use rightmost leaf node as the replacement node
  node_handle replace = takePathToID(root, -size);
  size++;

  // disconnect rightmost leaf node
  node_handle rmp = Parent(replace);
  if (rmp) {
    if (replace == Right(rmp)) {
      Right(rmp) = 0;
    } else {
      MEDDLY_DCASSERT(replace = Left(rmp));
      Left(rmp) = 0;
    }
  }

  // Lucky special case:
  // if we're deleting the rightmost leaf node, we're done!
  if (replace == p) {
    if (0==size) root = 0;  // empty heap now
    return;
  }
  
  // repair the heap downward:
  // where p was, put either left child, right child, or replacement
  node_handle* pptr = holeOf(p);
  node_handle newp = downHeap(pptr[parent_index], pptr[ID_index], 
    pptr[left_index], pptr[right_index], replace);

  // link this to where it goes
  if (0==pptr[parent_index]) {
    // p is the root; this is a common case
    root = newp;
  } else if (Right(pptr[parent_index]) == p) {
    // p was the right child
    Right(pptr[parent_index]) = newp;
  } else {
    // p was the left child
    MEDDLY_DCASSERT(Left(pptr[parent_index]) == p);
    Left(pptr[parent_index]) = newp;
  }
}

// ******************************************************************

MEDDLY::node_handle MEDDLY::hm_heap::downHeap(node_handle parent,
  node_handle ID, node_handle left, node_handle right, node_handle replace)
{
  MEDDLY_DCASSERT(replace);
  node_handle kl = left ? left : replace+1;
  node_handle kr = right ? right : replace+1;

  if (kl < kr) {
    if (replace < kl) {
        // replace is smallest
        makeHeapNode(replace, ID, parent, left, right);
        return replace;
    } else {
        // kl is smallest
        node_handle* lptr = holeOf(left);
        node_handle newleft = 
          downHeap(left, 2*ID, lptr[left_index], lptr[right_index], replace);
        fillHeapNode(lptr, ID, parent, newleft, right);
        return left;
      }
  } else {
      if (replace < kr) {
        // replace is smallest
        makeHeapNode(replace, ID, parent, left, right);
        return replace;
      } else {
        // kr is smallest
        node_handle* rptr = holeOf(right);
        node_handle newright = 
          downHeap(right, 2*ID+1, rptr[left_index], rptr[right_index], replace);
        fillHeapNode(rptr, ID, parent, left, newright);
        return right;
      }
  }
}




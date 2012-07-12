
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

/*! \file temp_dd_edge.cc

    Implementation for class temp_dd_edge described in meddly_expert.h
*/


#include "defines.h"

#ifdef USE_EXPERIMENTAL_TEMPEDGES

MEDDLY::temp_dd_edge::temp_dd_edge()
: levelHandle(0), forestHandle(0), iValue(0), rValue(0), size(0),
  downpointers(0), iEdgeValues(0), rEdgeValues(0)
{ }


MEDDLY::temp_dd_edge::~temp_dd_edge()
{
  if (downpointers) {
    /*
    for (int i = 0; i < size; ++i) {
      if (downpointers[i] != 0) delete downpointers[i];
    }
    */
    free(downpointers);
  }
  if (iEdgeValues) free(iEdgeValues);
  if (rEdgeValues) free(rEdgeValues);
}


bool MEDDLY::temp_dd_edge::convertToDDEdge(dd_edge& result) const
{
  // Pre-order traversal -- children first then parent.
  // i.e. depth-first.
  int node = 0;
  if (forestHandle->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    if (reduce(node)) {
      result.set(node, 0, forestHandle->getNodeLevel(node));
      return true;
    }
#if 0
  } else if (forestHandle->getEdgeLabeling() == forest::EVPLUS) {
    int ev = 0;
    if (reduceAndNormalize(node, ev)) {
      result.set(node, ev, forestHandle->getNodeLevel(node));
      return true;
    }
  } else {
    float ev = 0;
    if (reduceAndNormalize(node, ev)) {
      result.set(node, ev, forestHandle->getNodeLevel(node));
      return true;
    }
#endif
  }
  return false;
}


bool MEDDLY::temp_dd_edge::reduce(int& result) const
{
  // MDD
  // All temp nodes are full nodes.
  // Differentiate between MDD and MXD?

  // Compute table
  std::map<temp_dd_edge*, int> ct;

  // Set special case nodes

  int zero = 0;
  switch (forestHandle->getRangeType()) {
    //case forest::range_type::BOOLEAN:
    case forest::BOOLEAN:
      MEDDLY_DCASSERT(0 == forestHandle->getTerminalNode(false));
      zero = 0;
      break;

    //case forest::range_type::INTEGER:
    case forest::INTEGER:
      zero = forestHandle->getTerminalNode(0);
      break;

    //case forest::range_type::REAL:
    case forest::REAL:
      zero = forestHandle->getTerminalNode(0.0f);
      break;

    default:
      assert(false);
  }

  bool retVal = reduce(ct, zero, result);
  forestHandle->unlinkNode(zero);
  return retVal;
}


bool MEDDLY::temp_dd_edge::reduce(std::map<temp_dd_edge*, int>& ct, int zero,
    int& result) const
{
  MEDDLY_DCASSERT(forestHandle->getEdgeLabeling() == forest::MULTI_TERMINAL);

  // Special case: terminal node
  if (levelHandle == 0) {
    switch (forestHandle->getRangeType()) {
      //case forest::range_type::BOOLEAN:
      case forest::BOOLEAN:
        MEDDLY_DCASSERT(-1 == forestHandle->getTerminalNode(true));
        result = -1;
        break;

        //case forest::range_type::INTEGER:
      case forest::INTEGER:
        result = forestHandle->getTerminalNode((int)iValue);
        break;

        //case forest::range_type::REAL:
      case forest::REAL:
        result = forestHandle->getTerminalNode((float)rValue);
        break;

      default:
        assert(false);
    }
    return true;
  }

  MEDDLY_DCASSERT(levelHandle != 0 && size > 0);

  // Look for e in ct
  // If found return cached node
  // If not found build result and store in ct

  temp_dd_edge* e = const_cast<temp_dd_edge*>(&(*this));
  std::map<temp_dd_edge*, int>::iterator iter = ct.find(e);
  if (ct.end() != iter) {
    // found
    result = iter->second;
    forestHandle->linkNode(result);
    return true;
  }

  // Not found in compute table
  // Build the node
  result = forestHandle->createTempNode(levelHandle, size);
  int i = 0;
  for ( ; i < size; ++i)
  {
    if (0 == downpointers[i]) {
      forestHandle->setDownPtrWoUnlink(result, i, zero);
    }
    else {
      int temp = 0;
      if(!downpointers[i]->reduce(ct, zero, temp)) {
        // Not Well-Formed
        break;
      }
      forestHandle->setDownPtrWoUnlink(result, i, temp);
      forestHandle->unlinkNode(temp);
    }
  }

  if (i < size) {
    // Not Well-Formed
    // Delete temporary node: set all downpointers to zero and reduce.
    int j = 0;
    while (j < i) {
      forestHandle->setDownPtr(result, j++, zero);
    }
    while (j < size) {
      forestHandle->setDownPtrWoUnlink(result, j++, zero);
    }
  }

  result = forestHandle->reduceNode(result);
  ct[e] = result;

  // Return true for Well-Formed
  return i == size;
}


void MEDDLY::temp_dd_edge::add(const int* vlist, const int* vplist)
{
  // Add element to this tree.
  // Nodes in the tree are not-reduced.

  /*
     Find the index in vlist or vplist for this temp_dd_edge's level.

     If necessary, initialize downpointers[] and other arrays.
     The arrays should be at least of size index+1.

     If downpointers[index] == 0,
     Initialize downpointers[index] with an empty temp_dd_edge.

     Call downpointers[index]->add(vlist, vplist).
   */

  if (levelHandle == 0) return;

  MEDDLY_DCASSERT(forestHandle->isForRelations() || levelHandle > 0);

  int index = levelHandle > 0? vlist[levelHandle]: vplist[-levelHandle];

  if (index >= forestHandle->getLevelSize(levelHandle))
    throw error(error::INVALID_ASSIGNMENT);

  MEDDLY_DCASSERT((size > 0 && downpointers != 0) ||
      (size == 0 && downpointers == 0));

  if (size <= index) {
    // size could be 0.
    int newSize = MIN( forestHandle->getLevelSize(levelHandle),
        MAX( size*2, index+1 ) );
    downpointers = (temp_dd_edge**) realloc(downpointers,
        sizeof(temp_dd_edge*) * newSize);
    if (downpointers == 0) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);
    memset(downpointers + size, 0, sizeof(temp_dd_edge*) * (newSize - size));
    size = newSize;
  }

  if (downpointers[index] == 0) {
    downpointers[index] = new temp_dd_edge();
    downpointers[index]->forestHandle = forestHandle;
    if (forestHandle->isForRelations()) {
      downpointers[index]->levelHandle = 
        levelHandle < 0? (-levelHandle)-1 : -levelHandle;
    } else {
      downpointers[index]->levelHandle = levelHandle-1;
    }
  }

  downpointers[index]->add(vlist, vplist);
}

#endif

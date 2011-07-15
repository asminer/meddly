
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "union.h"
#include "apply_base.h"

namespace MEDDLY {
  class union_mdd;
  class union_mxd;

  class union_opname;
};

// ******************************************************************
// *                                                                *
// *                        union_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::union_mdd : public generic_binary_mdd {
  public:
    union_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::union_mdd::union_mdd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mdd(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::union_mdd::checkTerminals(int a, int b, int& c)
{
  if (a == -1 || b == -1) {
    c = -1;
    return true;
  }
  if (a == 0) {
    if (b==0) {
      c = 0;
      return true;
    }
    if (arg2F == resF) {
      c = resF->linkNode(b);
      return true;
    } 
    return false;
  }
  if (b == 0) {
    if (arg1F == resF) {
      c = resF->linkNode(a);
      return true;
    } 
    return false;
  }
  if (a == b) {
    if (arg1F == arg2F && arg1F == resF) {
      c = resF->linkNode(b);
      return true;
    } 
    return false;
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                        union_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::union_mxd : public generic_binary_mxd {
  public:
    union_mxd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
    virtual int computeIdent(int a, int b);
    virtual int computeIdentExpandA(int a, int b);
    virtual int computeIdentExpandOneLevel(int a, int b);
};

MEDDLY::union_mxd::union_mxd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mxd(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::union_mxd::checkTerminals(int a, int b, int& c)
{
  if (a == -1 && b == -1) {
    c = -1;
    return true;
  }
  if (0 == a) {
    if (0 == b) {
      c = 0;
      return true;
    }
    if (arg2F == resF) {
      c = resF->linkNode(b);
      return true;
    } else {
      return false;
    }
  }
  if (a == b) 
    if (arg1F == arg2F && arg1F == resF) {
      c = resF->linkNode(b);
      return true;
    } else {
      return false;
    }
  if (b == 0) {
    if (arg1F == resF) {
      c = resF->linkNode(a);
      return true;
    } else {
      return false;
    }
  }
  return false;
}

int MEDDLY::union_mxd::computeIdent(int a, int b)
{
  // Terminal conditions for recursion
  if (-1 == a || -1 == b) {
    return -1;
  }
  if (0 == a) {
    if (0 == b) return 0;
    if (arg2F == resF) return resF->linkNode(b);
  }
  if (a == b) {
    if (arg1F == arg2F && arg1F == resF) return resF->linkNode(b);
  }
  if (b == 0) {
    if (arg1F == resF) return resF->linkNode(a);
  }

  // Search compute table
  int result = 0;
#ifdef IGNORE_MAPPED_HEIGHT
  int aHeight = arg1F->getMappedNodeHeight(a);
  int bHeight = arg2F->getMappedNodeHeight(b);
  if (aHeight >= IGNORE_MAPPED_HEIGHT && bHeight >= IGNORE_MAPPED_HEIGHT) {
    if ((a > b
          ? owner->cc->find(owner, a, b, result)
          : owner->cc->find(owner, b, a, result))) {
      resF->linkNode(result);
      return result;
    }
  }
#else
  if (findResult(a, b, result)) {
    return result;
  }
  int aHeight = arg1F->getMappedNodeHeight(a);
  int bHeight = arg2F->getMappedNodeHeight(b);
#endif

  result = 0;

  if (aHeight > bHeight) {
    // result[i][j] = a[i][j]
    // except result[i][i] = a[i][i] + b
    result = computeIdentExpandA(a, b);
  } else if (aHeight < bHeight) {
    result = computeIdentExpandA(b, a);
  } else {
    result = computeIdentExpandOneLevel(b, a);
  }

  result = resF->reduceNode(result);

  // Save result to compute table
#ifdef IGNORE_MAPPED_HEIGHT
  if (aHeight < IGNORE_MAPPED_HEIGHT || bHeight < IGNORE_MAPPED_HEIGHT) {
    return result;
  }
#endif

#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d)=%d\n", getName(), a, b, result);
#endif

  saveResult(a, b, result);
  return result;
}

int MEDDLY::union_mxd::computeIdentExpandA(int a, int b)
{
  int resultLevel = arg1F->getNodeLevel(a);
  int resultSize = arg1F->getLevelSize(resultLevel);

  MEDDLY_DCASSERT(arg1F == resF);

  // Copy a into result.
  int result = arg1F->makeACopy(a, resultSize);
  int* resultDptrs = 0;
  assert(resF->getDownPtrs(result, resultDptrs));

  bool noChange = true;

  // Add b to result[i][i]
  int* rDptrs = resultDptrs;
  for (int i = 0; i < resultSize; ++i, ++rDptrs) {
    if (*rDptrs == 0) {
      int pNode = resF->createTempNode(-resultLevel, i + 1, true);
      resF->setDownPtrWoUnlink(pNode, i, b);
      *rDptrs = resF->reduceNode(pNode);
      noChange = false;
    } else {
      int mxdII = resF->getDownPtr(*rDptrs, i);
      int temp = 0;
      if (mxdII == 0) {
        temp = b; arg2F->linkNode(b);
      } else {
        temp = computeIdent(mxdII, b);
      }
      if (temp != mxdII) {
        int pNode = resF->makeACopy(*rDptrs, i + 1);
        resF->setDownPtr(pNode, i, temp);
        resF->unlinkNode(*rDptrs);
        *rDptrs = resF->reduceNode(pNode);
        noChange = false;
      }
      resF->unlinkNode(temp);
    }
  }

  if (noChange) {
    // result is the same as node a.
    // Don't call reduce; discard result and return a.
    arg1F->linkNode(a);
    resF->unlinkNode(result);
    result = a;
  } else {
    result = resF->reduceNode(result);
  }

  return result;
}

int MEDDLY::union_mxd::computeIdentExpandOneLevel(int a, int b)
{
  MEDDLY_DCASSERT(!arg1F->isTerminalNode(a));
  MEDDLY_DCASSERT(!arg2F->isTerminalNode(b));
  MEDDLY_DCASSERT(arg2F->getNodeLevel(b) == arg1F->getNodeLevel(a));

  MEDDLY_DCASSERT(arg1F == resF);

  int resultLevel = arg1F->getNodeLevel(a);
  int resultSize = arg1F->getLevelSize(resultLevel);

  // Set-up function pointer for call to the next level.
  // If a and b are at the unprime level, call computeIdentExpandOneLevel().
  // Otherwise, call computeIdent().
  int (union_mxd::*function)(int, int) =
    resultLevel > 0
    ? &union_mxd::computeIdentExpandOneLevel
    : &union_mxd::computeIdent;

  // Copy a into result.
  int result = arg1F->makeACopy(a, resultSize);
  int* resultDptrs = 0;
  assert(resF->getDownPtrs(result, resultDptrs));

  bool noChange = true;

  // Add b to result.
  int* rDptrs = resultDptrs;
  const int* bDptrs = 0;
  assert(arg2F->getDownPtrs(b, bDptrs));
  if (arg2F->isFullNode(b)) {
    int bSize = arg2F->getFullNodeSize(b);
    for (const int* bEnd = bDptrs + bSize; bDptrs != bEnd; ++rDptrs, ++bDptrs)
    {
      // Terminal conditions.
      if (*rDptrs == *bDptrs || 0 == *bDptrs) continue;
      if (0 == *rDptrs) {
        *rDptrs = *bDptrs;
        resF->linkNode(*rDptrs);
        noChange = false;
        continue;
      }
      // Expand *rDptrs and *bDptrs.
      int pNode = (this->*function)(*rDptrs, *bDptrs);
      if (*rDptrs != pNode) {
        resF->unlinkNode(*rDptrs);
        *rDptrs = pNode;
        noChange = false;
      } else {
        resF->unlinkNode(pNode);
      }
    }
  }
  else {
    MEDDLY_DCASSERT(arg2F->isSparseNode(b));
    int nDptrs = arg2F->getSparseNodeSize(b);
    const int* bIndexes = 0;
    assert(arg2F->getSparseNodeIndexes(b, bIndexes));
    for (const int* bEnd = bDptrs + nDptrs; bDptrs != bEnd; ++bDptrs)
    {
      // Terminal conditions
      rDptrs = resultDptrs + *bIndexes++;
      MEDDLY_DCASSERT(*bDptrs != 0);
      if (*rDptrs == *bDptrs) continue;
      if (0 == *rDptrs) {
        *rDptrs = *bDptrs;
        resF->linkNode(*rDptrs);
        noChange = false;
        continue;
      }
      // Expand *rDptrs and *bDptrs.
      int pNode = (this->*function)(*rDptrs, *bDptrs);
      if (*rDptrs != pNode) {
        resF->unlinkNode(*rDptrs);
        *rDptrs = pNode;
        noChange = false;
      } else {
        resF->unlinkNode(pNode);
      }
    }

  }

  if (noChange) {
    // result is the same as node a.
    // Don't call reduce; discard result and return a.
    arg1F->linkNode(a);
    resF->unlinkNode(result);
    result = a;
  } else {
    result = resF->reduceNode(result);
  }

  return result;
}




// ******************************************************************
// *                                                                *
// *                       union_opname class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::union_opname : public binary_opname {
  public:
    union_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::union_opname::union_opname()
 : binary_opname("Union")
{
}

MEDDLY::binary_operation* 
MEDDLY::union_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH);

  if (
    (a1->isForRelations() != r->isForRelations()) ||
    (a2->isForRelations() != r->isForRelations()) ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a2->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (r->getRangeType() != forest::BOOLEAN)
  )
    throw error(error::TYPE_MISMATCH);

  if (r->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    if (r->isForRelations())
      return new union_mxd(this, a1, a2, r);
    else
      return new union_mdd(this, a1, a2, r);
  }

  throw error(error::NOT_IMPLEMENTED);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeUnion(const settings &s)
{
  return new union_opname;
}



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
#include "prepostimage.h"

// #define TRACE_ALL_OPS

namespace MEDDLY {
  class image_op;

  class preimage_mdd;
  class postimage_mdd;

  class preimage_opname;
  class postimage_opname;
};

// ******************************************************************
// *                                                                *
// *                         image_op class                         *
// *                                                                *
// ******************************************************************

/// Abstract base class for all pre/post image operations.
class MEDDLY::image_op : public binary_operation {
  public:
    image_op(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual bool isStaleEntry(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;

    inline bool findResult(int a, int b, int &c) {
      CTsrch.key(0) = a;
      CTsrch.key(1) = b;
      const int* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[2]);
      return true;
    }
    inline int saveResult(int a, int b, int c) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      entry.key(0) = arg1F->cacheNode(a); 
      entry.key(1) = arg2F->cacheNode(b);
      entry.result(0) = resF->cacheNode(c);
      CT->addEntry();
      return c;
    }
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual int compute(int a, int b);
  protected:
    binary_operation* unionOp;
    virtual int compute_rec(int a, int b) = 0;
};

MEDDLY::image_op::image_op(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc, 2, 1, a1, a2, res)
{
  unionOp = 0;
}

bool MEDDLY::image_op::isStaleEntry(const int* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}

void MEDDLY::image_op::discardEntry(const int* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void
MEDDLY::image_op::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d, %d): %d]", getName(), data[0], data[1], data[2]);
}

void MEDDLY::image_op
::compute(const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  int cnode = compute(a.getNode(), b.getNode());
  c.set(cnode, 0, resF->getNodeLevel(cnode));
}

int MEDDLY::image_op::compute(int a, int b)
{
  if (resF->getRangeType() == forest::BOOLEAN) {
    unionOp = getOperation(UNION, resF, resF, resF);
  } else {
    unionOp = getOperation(MAXIMUM, resF, resF, resF);
  }
  MEDDLY_DCASSERT(unionOp);
  return compute_rec(a, b);
}

// ******************************************************************
// *                                                                *
// *                       preimage_mdd class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::preimage_mdd : public image_op {
  public:
    preimage_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual int compute_rec(int a, int b);
};

MEDDLY::preimage_mdd::preimage_mdd(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : image_op(oc, a1, a2, res)
{
}

int MEDDLY::preimage_mdd::compute_rec(int mdd, int mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (arg2F->isTerminalNode(mxd)) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->getTerminalNode(true);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  int result = 0;
  if (findResult(mdd, mxd, result)) {
    return result;
  }

  // check if mxd and mdd are at the same level
  int mddLevel = arg1F->getNodeLevel(mdd);
  int mxdLevel = arg2F->getNodeLevel(mxd);
  int rLevel = MAX(ABS(mxdLevel), mddLevel);
  int rSize = resF->getLevelSize(rLevel);
  expert_forest::nodeBuilder& nb = resF->useNodeBuilder(rLevel, rSize);

  // Initialize mdd reader
  node_reader* A = (mddLevel < rLevel)
    ? arg1F->initRedundantReader(rLevel, mdd, true)
    : arg1F->initNodeReader(mdd, true);

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (int i=0; i<rSize; i++) {
      nb.d(i) = compute_rec(A->d(i), mxd);
    }
  } else {
    // 
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (int i=0; i<rSize; i++) nb.d(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    node_reader* Ru = (mxdLevel < 0)
      ? arg2F->initRedundantReader(rLevel, mxd, false)
      : arg2F->initNodeReader(mxd, false);

    node_reader* Rp = arg2F->useNodeReader();

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->d(iz)))) {
        arg2F->initIdentityReader(*Rp, rLevel, i, Ru->d(iz), false);
      } else {
        arg2F->initNodeReader(*Rp, Ru->d(iz), false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
        if (0==A->d(j))   continue; 
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        int newstates = compute_rec(A->d(j), Rp->d(jz));
        if (0==newstates) continue;
        if (0==nb.d(i)) {
          nb.d(i) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        int oldi = nb.d(i);
        nb.d(i) = unionOp->compute(newstates, oldi);
        resF->unlinkNode(oldi);
        resF->unlinkNode(newstates);
      } // for j
  
    } // for i

    arg2F->recycle(Rp);
    arg2F->recycle(Ru);
  } // else

  // cleanup mdd reader
  arg1F->recycle(A);

  result = resF->createReducedNode(-1, nb);
#ifdef TRACE_ALL_OPS
  printf("computed preimage(%d, %d) = %d\n", mdd, mxd, result);
#endif
  return saveResult(mdd, mxd, result); 
}


// ******************************************************************
// *                                                                *
// *                      postimage_mdd  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::postimage_mdd : public image_op {
  public:
    postimage_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual int compute_rec(int a, int b);
};

MEDDLY::postimage_mdd::postimage_mdd(const binary_opname* oc, 
  expert_forest* a1, expert_forest* a2, expert_forest* res)
: image_op(oc, a1, a2, res)
{
}

int MEDDLY::postimage_mdd::compute_rec(int mdd, int mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (arg2F->isTerminalNode(mxd)) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->getTerminalNode(true);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  int result = 0;
  if (findResult(mdd, mxd, result)) {
    return result;
  }

  // check if mxd and mdd are at the same level
  int mddLevel = arg1F->getNodeLevel(mdd);
  int mxdLevel = arg2F->getNodeLevel(mxd);
  int rLevel = MAX(ABS(mxdLevel), mddLevel);
  int rSize = resF->getLevelSize(rLevel);
  expert_forest::nodeBuilder& nb = resF->useNodeBuilder(rLevel, rSize);

  // Initialize mdd reader
  node_reader* A = (mddLevel < rLevel)
    ? arg1F->initRedundantReader(rLevel, mdd, true)
    : arg1F->initNodeReader(mdd, true);

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (int i=0; i<rSize; i++) {
      nb.d(i) = compute_rec(A->d(i), mxd);
    }
  } else {
    // 
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (int i=0; i<rSize; i++) nb.d(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    node_reader* Ru = (mxdLevel < 0)
      ? arg2F->initRedundantReader(rLevel, mxd, false)
      : arg2F->initNodeReader(mxd, false);

    node_reader* Rp = arg2F->useNodeReader();

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      if (0==A->d(i))   continue; 
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->d(iz)))) {
        arg2F->initIdentityReader(*Rp, rLevel, i, Ru->d(iz), false);
      } else {
        arg2F->initNodeReader(*Rp, Ru->d(iz), false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        int newstates = compute_rec(A->d(i), Rp->d(jz));
        if (0==newstates) continue;
        if (0==nb.d(j)) {
          nb.d(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        int oldj = nb.d(j);
        nb.d(j) = unionOp->compute(newstates, oldj);
        resF->unlinkNode(oldj);
        resF->unlinkNode(newstates);
      } // for j
  
    } // for i

    arg2F->recycle(Rp);
    arg2F->recycle(Ru);
  } // else

  // cleanup mdd reader
  arg1F->recycle(A);

  result = resF->createReducedNode(-1, nb);
#ifdef TRACE_ALL_OPS
  printf("computed new postimage(%d, %d) = %d\n", mdd, mxd, result);
#endif
  return saveResult(mdd, mxd, result); 
}


// ******************************************************************
// *                                                                *
// *                     preimage_opname  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::preimage_opname : public binary_opname {
  public:
    preimage_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::preimage_opname::preimage_opname()
 : binary_opname("Pre-image")
{
}

MEDDLY::binary_operation* 
MEDDLY::preimage_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a2->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH);

  return new preimage_mdd(this, a1, a2, r);
}


// ******************************************************************
// *                                                                *
// *                     postimage_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::postimage_opname : public binary_opname {
  public:
    postimage_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::postimage_opname::postimage_opname()
 : binary_opname("Post-image")
{
}

MEDDLY::binary_operation* 
MEDDLY::postimage_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a2->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH);

  return new postimage_mdd(this, a1, a2, r);
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializePreImage(const settings &s)
{
  return new preimage_opname;
}

MEDDLY::binary_opname* MEDDLY::initializePostImage(const settings &s)
{
  return new postimage_opname;
}


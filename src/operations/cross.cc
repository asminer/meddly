
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
#include "cross.h"
#include "../compute_cache.h"

// ******************************************************************
// *                                                                *
// *                         cross_op class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

/** Abstract base class for cross operations.
*/
class cross_op : public operation {
  public:
    cross_op()                              { }
    virtual ~cross_op()                     { }

    virtual const char* getName() const     { return "Cross"; }
    virtual int getKeyLength() const        { return 3; }
    virtual int getKeyLengthInBytes() const { return 3*sizeof(int); }

    virtual int getAnsLength() const        { return 1; }
    virtual int getAnsLengthInBytes() const { return sizeof(int); }  

    virtual int getCacheEntryLength() const         { return 4; }
    virtual int getCacheEntryLengthInBytes() const  { return 4*sizeof(int); }

    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);

    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    inline bool findResult(op_info* owner, int k, int a, int b, int& c) {
      static int key[3];
      key[0] = k;
      key[1] = a;
      key[2] = b;
      const int* cacheEntry = owner->cc->find(owner, key);
      if (0==cacheEntry) return false;
      c = cacheEntry[3];
      owner->p[2].getForest()->linkNode(c);
      return true;
    }

    inline void saveResult(op_info* owner, int k, int a, int b, int c) {
      owner->p[0].getForest()->cacheNode(a);
      owner->p[1].getForest()->cacheNode(b);
      owner->p[2].getForest()->cacheNode(c);
      static int entry[4];
      entry[0] = k;
      entry[1] = a;
      entry[2] = b;
      entry[3] = c;
      owner->cc->add(owner, entry);
    }

  protected:
    inline void 
    type_check(const op_info* o, forest::range_type t, forest::edge_labeling e) 
    {
        if (o == 0)
          throw error(error::UNKNOWN_OPERATION);
        if (o->op == 0 || o->p == 0 || o->cc == 0)
          throw error(error::TYPE_MISMATCH);
        if (o->nParams != 3)
          throw error(error::WRONG_NUMBER);
        if (o->p[0].isForestOf(false, t, e) &&
            o->p[1].isForestOf(false, t, e) &&
            o->p[2].isForestOf(true,  t, e)) 
        {
          return;
        } else {
          throw error(error::TYPE_MISMATCH);
        }
    }
};

} // namespace MEDDLY

bool MEDDLY::cross_op::
isEntryStale(const op_info* owner, const int* data)
{
  DCASSERT(owner->nParams == 3);
  // data[0] is the level number
  return owner->p[0].getForest()->isStale(data[1]) ||
         owner->p[1].getForest()->isStale(data[2]) ||
         owner->p[2].getForest()->isStale(data[3]);
}

void MEDDLY::cross_op::
discardEntry(op_info* owner, const int* data)
{
  DCASSERT(owner->nParams == 3);
  // data[0] is the level number
  owner->p[0].getForest()->uncacheNode(data[1]);
  owner->p[1].getForest()->uncacheNode(data[2]);
  owner->p[2].getForest()->uncacheNode(data[3]);
}

void
MEDDLY::cross_op::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  DCASSERT(owner->nParams == 3);
  fprintf(strm, "[%s(%d, %d, %d): %d]",
      owner->op->getName(), data[0], data[1], data[2], data[3]
  );
}


// ******************************************************************
// *                                                                *
// *                        bool_cross class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

/// Boolean cross product operation.
class bool_cross : public cross_op {
  public:
    bool_cross()            { }
    virtual ~bool_cross()   { }

    static bool_cross* getInstance();

    virtual void typeCheck(const op_info* owner) {
      return type_check(owner, forest::BOOLEAN, forest::MULTI_TERMINAL);
    }

    virtual void compute(op_info* owner, const dd_edge &a,
      const dd_edge &b, dd_edge &c);

    int compute_unprimed(op_info* owner, int ht, int a, int b);
    int compute_primed(op_info* owner, int ht, int a, int b);
};

} // namespace MEDDLY

MEDDLY::bool_cross*
MEDDLY::bool_cross::
getInstance()
{
  static bool_cross instance;
  return &instance;
}

void
MEDDLY::bool_cross::
compute(op_info* owner, const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  if (0==owner) throw error(error::TYPE_MISMATCH);
  int L = owner->p[0].getForest()->getDomain()->getNumVariables();
  int cnode = compute_unprimed(owner, L, a.getNode(), b.getNode());
  c.set(cnode, 0, owner->p[2].getForest()->getNodeLevel(cnode));
}

int
MEDDLY::bool_cross::
compute_unprimed(op_info* owner, int lh, int a, int b)
{
  if (0==lh) {
    return a;
  }
  if (0==a || 0==b) return 0;

  // convert height to level handle
  expert_forest* fc = owner->p[2].getForest();

  // check compute table
  int c;
  if (findResult(owner, lh, a, b, c)) return c;

  // build new result node
  c = fc->createTempNodeMaxSize(lh, true);

  // recurse
  expert_forest* fa = owner->p[0].getForest();
  if (fa->getNodeHeight(a) < lh) {
    // skipped level
    int d = compute_primed(owner, lh, a, b);
    for (int i=fc->getLevelSize(lh)-1; i>=0; i--) {
      fc->setDownPtr(c, i, d);
    } 
  } else {
    // not skipped level
    if (fa->isFullNode(a)) {
      // Full storage
      for (int i=fa->getFullNodeSize(a)-1; i>=0; i--) {
        int ai = fa->getFullNodeDownPtr(a, i);
        if (0==ai) continue;
        fc->setDownPtr(c, i, compute_primed(owner, lh, ai, b));
      }
    } else {
      // Sparse storage
      for (int z=fa->getSparseNodeSize(a)-1; z>=0; z--) {
        int i = fa->getSparseNodeIndex(a, z);
        int ai = fa->getSparseNodeDownPtr(a, z);
        fc->setDownPtr(c, i, compute_primed(owner, lh, ai, b));
      }
    }
  }

  // reduce
  c = fc->reduceNode(c);

  // add to compute table
  saveResult(owner, lh, a, b, c);

  return c;
}

int
MEDDLY::bool_cross::
compute_primed(op_info* owner, int ht, int a, int b)
{
  if (0==a || 0==b) return 0;

  // convert height to level handle
  expert_forest* fc = owner->p[2].getForest();
  int lh = -ht;

  // check compute table
  int c;
  if (findResult(owner, lh, a, b, c)) return c;

  // build new result node
  c = fc->createTempNodeMaxSize(lh, true);

  // recurse
  ht--;
  expert_forest* fb = owner->p[0].getForest();
  if (fb->getNodeHeight(b) <= ht) {
    // skipped level
    int d = compute_unprimed(owner, ht, a, b);
    for (int i=fc->getLevelSize(lh)-1; i>=0; i--) {
      fc->setDownPtr(c, i, d);
    } 
  } else {
    // not skipped level
    if (fb->isFullNode(b)) {
      // Full storage
      for (int i=fb->getFullNodeSize(b)-1; i>=0; i--) {
        int bi = fb->getFullNodeDownPtr(b, i);
        if (0==bi) continue;
        fc->setDownPtr(c, i, compute_unprimed(owner, ht, a, bi));
      }
    } else {
      // Sparse storage
      for (int z=fb->getSparseNodeSize(b)-1; z>=0; z--) {
        int i = fb->getSparseNodeIndex(b, z);
        int bi = fb->getSparseNodeDownPtr(b, z);
        fc->setDownPtr(c, i, compute_unprimed(owner, ht, a, bi));
      }
    }
  }

  // reduce
  c = fc->reduceNode(c);

  // add to compute table
  saveResult(owner, lh, a, b, c);

  return c;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::operation* MEDDLY::getCrossOperation(const op_param &ot)
{
  if (!ot.isForest()) return 0;
  if (!ot.isBoolForest()) return 0;

  return bool_cross::getInstance();
}



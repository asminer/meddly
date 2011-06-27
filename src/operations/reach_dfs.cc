
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
#include "reach_dfs.h"
#include "../compute_table.h"

#include <vector>
#include <map>
#include <set>

#define ALT_SATURATE_HELPER

namespace MEDDLY {
  class common_dfs_mt;
  class forwd_dfs_mt;
  class bckwd_dfs_mt;

  class forwd_dfs_opname;
  class bckwd_dfs_opname;
};


// ******************************************************************
// *                                                                *
// *                      common_dfs_mt  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::common_dfs_mt : public binary_operation {
  public:
    common_dfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int* entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual int compute(int a, int b) = 0;
};

MEDDLY::common_dfs_mt::common_dfs_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : binary_operation(oc, a1, a2, res)
{
  key_length = 2;
  ans_length = 1;
}

bool MEDDLY::common_dfs_mt::isEntryStale(const int* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}

void MEDDLY::common_dfs_mt::discardEntry(const int* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void MEDDLY::common_dfs_mt::showEntry(FILE* strm, const int* data) const
{
  fprintf(strm, "[%s(%d, %d): %d]", getName(), data[0], data[1], data[2]);
}

void MEDDLY::common_dfs_mt
::compute(const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  int cnode = compute(a.getNode(), b.getNode());
  c.set(cnode, 0, resF->getNodeLevel(cnode));
}


// ******************************************************************
// *                                                                *
// *                       forwd_dfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_dfs_mt : public common_dfs_mt {
  public:
    forwd_dfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);
  protected:
    virtual int compute(int a, int b);

    virtual void initialize();
    virtual void clear();
    virtual void clearSplitMxdComputeTableEntries();

    virtual void splitMxd(int mxd);
    virtual int saturate(int mdd);
#ifndef ALT_SATURATE_HELPER
    virtual void saturateHelper(int mddLevel, std::vector<int>& mdd);
    virtual int recFire(int mdd, int mxd);
#else
    void saturateHelper(int mdd);
    int recFire(int mdd, int mxd);

    void recFireExpandMdd(int mdd, int mxd, int& result);
    void recFireExpandMxd(int mdd, int mxd, int& result);
    void recFireExpand(int mdd, int mxd, int& result);

#endif

    bool getMxdAsVec(int mxd, int**& vec, int& size);

    // Helper for getMxdAsVec().
    // Use getMxd(0, 0, true) to clear the static arrays in getMxd().
    int** getMatrix(unsigned nodeLevel, unsigned size, bool clear);


    inline int getMddUnion(int a, int b) {
      DCASSERT(ddf->isTerminalNode(a) || ddf->isReducedNode(a));
      DCASSERT(ddf->isTerminalNode(b) || ddf->isReducedNode(b));
      DCASSERT(mddUnion);
      return mddUnion->compute(a, b);
    };
    inline int getMxdIntersection(int a, int b) {
      DCASSERT(xdf->isTerminalNode(a) || xdf->isReducedNode(a));
      DCASSERT(xdf->isTerminalNode(b) || xdf->isReducedNode(b));
      DCASSERT(mxdIntersection);
      return mxdIntersection->compute(a, b);
    }
    inline int getMxdDifference(int a, int b) {
      DCASSERT(xdf->isTerminalNode(a) || xdf->isReducedNode(a));
      DCASSERT(xdf->isTerminalNode(b) || xdf->isReducedNode(b));
      DCASSERT(mxdDifference);
      return mxdDifference->compute(a, b);
    }

    inline bool findSaturateResult(int a, int& b) {
      std::map<int, int>::iterator iter = saturateCT.find(a);
      if (iter == saturateCT.end()) return false;
      if (resF->isStale(iter->second)) {
        arg1F->uncacheNode(iter->first);
        resF->uncacheNode(iter->second);
        saturateCT.erase(iter);
        return false;
      }
      b = iter->second;
      return true;
    }
    inline void saveSaturateResult(int a, int b) {
      DCASSERT(! findSaturateResult(a, b));
      saturateCT[a] = b;
      arg1F->cacheNode(a);
      resF->cacheNode(b);
    }
    inline void clearSaturateCT() {
      std::map<int, int>::iterator iter = saturateCT.begin();
      std::map<int, int>::iterator end = saturateCT.end();
      while (iter != end) {
        arg1F->uncacheNode(iter->first);
        resF->uncacheNode(iter->second);
        saturateCT.erase(iter++);
      }
    }

  private:

    expert_domain* ed;            // domain

    // Next-state function is split and stored here (see Saturation algorithm).
    std::vector<int> splits;

    // scratch.size () == number of variable handles in the domain
    // scratch[level_handle].size() == level_bound(level_handle)
    int**             scratch0;
    int**             scratch1;
    int**             scratch2;
    int               scratchSize;

    binary_operation* mddUnion;
    binary_operation* mxdIntersection;
    binary_operation* mxdDifference;

    std::map<int, int>  saturateCT;
};

// ******************************************************************
// *                                                                *
// *                     forwd_dfs_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_dfs_opname : public binary_opname {
  public:
    forwd_dfs_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::forwd_dfs_opname::forwd_dfs_opname()
 : binary_opname("ReachableDFS")
{
}

MEDDLY::binary_operation* 
MEDDLY::forwd_dfs_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
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

  return new forwd_dfs_mt(this, a1, a2, r);
}

// ******************************************************************
// *                                                                *
// *                     bckwd_dfs_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_dfs_opname : public binary_opname {
  public:
    bckwd_dfs_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::bckwd_dfs_opname::bckwd_dfs_opname()
 : binary_opname("ReverseReachableDFS")
{
}

MEDDLY::binary_operation* 
MEDDLY::bckwd_dfs_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH);

  if (a1 != r)
    throw error(error::FOREST_MISMATCH);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    (a1->getRangeType() != a2->getRangeType()) ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) 
  )
    throw error(error::TYPE_MISMATCH);

  throw error(error::NOT_IMPLEMENTED);
  // return new bckwd_dfs_mt(this, a1, a2, r);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeForwardDFS(const settings &s)
{
  // TBD
  return 0;
}

MEDDLY::binary_opname* MEDDLY::initializeBackwardDFS(const settings &s)
{
  // TBD
  return 0;
}


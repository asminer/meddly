
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



#ifndef REACHABILITY_H
#define REACHABILITY_H

#include "prepostimage.h"
#include "../compute_cache.h"
#include <vector>
#include <map>

#define ALT_SATURATE_HELPER


// **********************************************************************
//
//                        Helper functions
//
// **********************************************************************

unsigned getNodeSize(expert_forest* f, int node);
void clearVector(std::vector<int>& v, unsigned sz);




// **********************************************************************
//
//            Forward Reachability: Traditional Algorithm
//
// **********************************************************************

class mdd_reachability_bfs : public mdd_mxd_image_operation {
    binary_operation* unionOp;
  public:
    static mdd_reachability_bfs* getInstance();
    int compute(op_info* owner, int a, int b);
    virtual const char* getName() const { return "Mdd-Mxd Reachability BFS"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mdd_reachability_bfs();
    mdd_reachability_bfs(const mdd_reachability_bfs& copy);
    mdd_reachability_bfs& operator=(const mdd_reachability_bfs& copy);
    virtual ~mdd_reachability_bfs();
};




// **********************************************************************
//
//            Forward Reachability: Saturation
//
// **********************************************************************

class mdd_reachability_dfs : public mdd_mxd_image_operation {
  public:
    static mdd_reachability_dfs* getInstance();
    int compute(op_info* owner, int a, int b);
    virtual const char* getName() const { return "Mdd-Mxd Reachability DFS"; }
    virtual bool isCommutative() const { return false; }
    virtual bool findResult(op_info* owner, int a, int b, int& c);
    virtual void saveResult(op_info* owner, int a, int b, int c);

  public:
    mdd_reachability_dfs();
    mdd_reachability_dfs(const mdd_reachability_dfs& copy);
    mdd_reachability_dfs& operator=(const mdd_reachability_dfs& copy);
    virtual ~mdd_reachability_dfs();

    virtual void initialize(op_info* owner);
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

    virtual int getMddUnion(int a, int b);
    virtual int getMxdIntersection(int a, int b);
    virtual int getMxdDifference(int a, int b);

    op_info*       owner;         // pointer to dfs reachability operation
    expert_forest* ddf;           // MDD forest
    expert_forest* xdf;           // MXD forest
    expert_domain* ed;            // domain
    expert_compute_manager* ecm;  // compute manager

    // Next-state function is split and stored here (see Saturation algorithm).
    std::vector<int> splits;

    // scratch.size () == number of variable handles in the domain
    // scratch[level_handle].size() == level_bound(level_handle)
    int**             scratch0;
    int**             scratch1;
    int**             scratch2;
    int               scratchSize;

    op_info*          mddUnionOp;
    op_info*          mxdIntersectionOp;
    op_info*          mxdDifferenceOp;

    mdd_union*        mddUnion;
    mxd_intersection* mxdIntersection;
    mxd_difference*   mxdDifference;

    std::map<int, int>  saturateCT;
    bool findSaturateResult(int a, int& b);
    void saveSaturateResult(int a, int b);
    void clearSaturateCT();

    bool getMxdAsVec(int mxd, int**& vec, int& size);

    // Helper for getMxdAsVec().
    // Use getMxd(0, 0, true) to clear the static arrays in getMxd().
    int** getMatrix(unsigned nodeLevel, unsigned size, bool clear);
};




// **********************************************************************
//
//            Backward Reachability: Traditional Algorithm
//
// **********************************************************************

class mdd_backward_reachability_bfs : public mdd_reachability_bfs {
  public:
    static mdd_backward_reachability_bfs* getInstance();
    int compute(op_info* owner, int a, int b);
    virtual const char* getName() const {
      return "Mdd-Mxd Reverse Reachability BFS"; }
};




// **********************************************************************
//
//            Backward Reachability: Saturation
//
// **********************************************************************

// TODO: Make this faster by applying the techniques used with
// mdd_reachability_dfs (based on per-operation compute table and getMatrix).
class mdd_backward_reachability_dfs : public mdd_reachability_dfs {
  public:
    static mdd_backward_reachability_dfs* getInstance();
    virtual const char* getName() const {
      return "Mdd-Mxd Reverse Reachability DFS";
    }

    virtual int saturate(int mdd);
    virtual void saturateHelper(int mddLevel, std::vector<int>& mdd) {
      reverseSaturateHelper(mddLevel, mdd);
    }
    virtual void reverseSaturateHelper(int mddLevel, std::vector<int>& mdd);
    virtual int reverseRecFire(int mdd, int mxd);
};











// **********************************************************************
//
//                            Inlined Methods
//
// **********************************************************************


inline
unsigned getNodeSize(expert_forest* f, int node)
{
  return unsigned(
      f->isFullNode(node)
      ? f->getFullNodeSize(node)
      : 1 + f->getSparseNodeIndex(node, f->getSparseNodeSize(node) - 1));
}

inline
void clearVector(std::vector<int>& v, unsigned sz)
{
  if (sz < v.size()) {
    fill(v.begin(), v.begin() + sz, 0);
  } else {
    fill(v.begin(), v.end(), 0);
  }
}

inline 
int mdd_reachability_dfs::getMddUnion(int a, int b)
{
  DCASSERT(ddf->isTerminalNode(a) || ddf->isReducedNode(a));
  DCASSERT(ddf->isTerminalNode(b) || ddf->isReducedNode(b));
  return mddUnion->compute(mddUnionOp, a, b);
}

inline
int mdd_reachability_dfs::getMxdIntersection(int a, int b)
{
  DCASSERT(xdf->isTerminalNode(a) || xdf->isReducedNode(a));
  DCASSERT(xdf->isTerminalNode(b) || xdf->isReducedNode(b));
  return mxdIntersection->compute(mxdIntersectionOp, a, b);
}

inline
int mdd_reachability_dfs::getMxdDifference(int a, int b)
{
  DCASSERT(xdf->isTerminalNode(a) || xdf->isReducedNode(a));
  DCASSERT(xdf->isTerminalNode(b) || xdf->isReducedNode(b));
  return mxdDifference->compute(mxdDifferenceOp, a, b);
}


inline
void mdd_reachability_dfs::saveSaturateResult(int a, int b)
{
  DCASSERT(! findSaturateResult(a, b));
  saturateCT[a] = b;
  ddf->cacheNode(a);
  ddf->cacheNode(b);
}

inline
bool mdd_reachability_dfs::findSaturateResult(int a, int& b)
{
  std::map<int, int>::iterator iter = saturateCT.find(a);
  if (iter == saturateCT.end()) return false;
  if (ddf->isStale(iter->second)) {
    ddf->uncacheNode(iter->first);
    ddf->uncacheNode(iter->second);
    saturateCT.erase(iter);
    return false;
  }
  b = iter->second;
  return true;
}

inline
void mdd_reachability_dfs::clearSaturateCT()
{
  std::map<int, int>::iterator iter = saturateCT.begin();
  std::map<int, int>::iterator end = saturateCT.end();

  while (iter != end) {
    ddf->uncacheNode(iter->first);
    ddf->uncacheNode(iter->second);
    saturateCT.erase(iter++);
  }
}

#endif


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


//#define ALT_SATURATE_HELPER

// ---------------------- Forward Reachability -------------------

// Traditional reachability analysis
class mdd_reachability_bfs : public mdd_mxd_image_operation {
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


// Reachability via "saturation" algorithm
class mdd_reachability_dfs : public mdd_mxd_image_operation {
  public:
    static mdd_reachability_dfs* getInstance();
    int compute(op_info* owner, int a, int b);
    virtual const char* getName() const { return "Mdd-Mxd Reachability DFS"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mdd_reachability_dfs();
    mdd_reachability_dfs(const mdd_reachability_dfs& copy);
    mdd_reachability_dfs& operator=(const mdd_reachability_dfs& copy);
    virtual ~mdd_reachability_dfs();

    virtual void initialize(op_info* owner);
    virtual void clear();

    virtual void splitMxd(int mxd);
    virtual int saturate(int mdd);
    virtual int recFire(int mdd, int mxd);
#ifdef ALT_SATURATE_HELPER
    virtual void saturateHelper(int mdd);
    virtual void saturateHelperUnPrimeFull(int mdd, int mxd);
    virtual void saturateHelperUnPrimeSparse(int mdd, int mxd);
    virtual void saturateHelperPrimeFull(int mdd, int i, int mxdI,
      std::vector<bool>& next);
    virtual void saturateHelperPrimeSparse(int mdd, int i, int mxdI,
      std::vector<bool>& next);

    virtual void recFireExpandMdd(int mdd, int mxd, int result);
    virtual void recFireExpandMxd(int mdd, int mxd, int result);
    virtual void recFireFF(int mdd, int mxd, int result);
    virtual void recFireFS(int mdd, int mxd, int result);
    virtual void recFireSF(int mdd, int mxd, int result);
    virtual void recFireSS(int mdd, int mxd, int result);
    virtual void recFirePrime(int mdd, int mxd, int result);
#else
    virtual void saturateHelper(int mddLevel, std::vector<int>& mdd);
#endif

    virtual int getMddUnion(int a, int b);
    virtual int getMxdIntersection(int a, int b);
    virtual int getMxdDifference(int a, int b);

    virtual unsigned getNodeSize(expert_forest* f, int node) const;
    virtual void clearVector(std::vector<int>& v, unsigned sz) const;

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
};

inline
unsigned mdd_reachability_dfs::getNodeSize(expert_forest* f, int node) const
{
  return unsigned(
      f->isFullNode(node)
      ? f->getFullNodeSize(node)
      : 1 + f->getSparseNodeIndex(node, f->getSparseNodeSize(node) - 1));
}

inline
void mdd_reachability_dfs::clearVector(std::vector<int>& v, unsigned sz) const
{
  if (sz < v.size()) {
    fill(v.begin(), v.begin() + sz, 0);
  } else {
    fill(v.begin(), v.end(), 0);
  }
}

// ---------------------- Backward Reachability -------------------

// Traditional backward reachability analysis
class mdd_reversereach_bfs : public mdd_mxd_image_operation {
  public:
    static mdd_reversereach_bfs* getInstance();
    int compute(op_info* owner, int a, int b);
    virtual const char* getName() const { return "Mdd-Mxd Reverse Reachability BFS"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mdd_reversereach_bfs();
    mdd_reversereach_bfs(const mdd_reversereach_bfs& copy);
    mdd_reversereach_bfs& operator=(const mdd_reversereach_bfs& copy);
    virtual ~mdd_reversereach_bfs();
};



#endif

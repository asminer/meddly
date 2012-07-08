
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

#include <vector>
#include <map>
#include <set>

// #define TRACE_RECFIRE
// #define DEBUG_DFS

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

    virtual bool isStaleEntry(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int* entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual int compute(int a, int b) = 0;

  protected:
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
};

MEDDLY::common_dfs_mt::common_dfs_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc, 2, 1, a1, a2, res)
{
}

bool MEDDLY::common_dfs_mt::isStaleEntry(const int* data)
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

    void saturateHelper(expert_forest::nodeBuilder& mdd);
    int recFire(int mdd, int mxd);

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
      MEDDLY_DCASSERT(! findSaturateResult(a, b));
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

  protected:
    // Next-state function is split and stored here (see Saturation algorithm).
    std::vector<int> splits;

    binary_operation* mddUnion;
    binary_operation* mxdIntersection;
    binary_operation* mxdDifference;

    std::map<int, int>  saturateCT;
};

MEDDLY::forwd_dfs_mt::forwd_dfs_mt(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : common_dfs_mt(opcode, arg1, arg2, res)
{
}

int MEDDLY::forwd_dfs_mt::compute(int mdd, int mxd)
{
  // Initialize class members and helper operations
  initialize();

  // Depth-first reachability analysis (Saturation)

#ifdef DEBUG_DFS
  printf("Consolidated Next-State Function:\n");
  arg2F->showNodeGraph(stdout, mxd);
  printf("\n");

  printf("Initial State:\n");
  arg1F->showNodeGraph(stdout, mdd);
  printf("\n");
#endif

  // Split the next-state function: each level has its own next-state function
  // The nsf is stored into the vector splits
  splitMxd(mxd);
  
  // Clear the splitMxd() related compute table entries
  clearSplitMxdComputeTableEntries();

#ifdef DEBUG_DFS
  printf("Split Next-State Function:\n");
  for (int i = splits.size() - 1; i >= 0; i--)
  {
    printf("Level %d, Node %d\n", i, splits[i]);
    arg2F->showNodeGraph(stdout, splits[i]);
    printf("\n");
  }

  fflush(stdout);
#endif

#ifdef DEBUG_SPLITS
  for (unsigned i = 0; i < splits.size(); i++) {
    fprintf(stdout, "# of nodes in the next-state function: %1.6e\n",
        double(arg2F->getNodeCount(splits[i])));
  }

  for (unsigned i = 0; i < splits.size(); i++) {
    fprintf(stdout, "# of edges in the next-state function: %1.6e\n",
        double(arg2F->getEdgeCount(splits[i], false)));
  }
#endif

  // Saturate the node
  int result = saturate(mdd);

  // clear pointers to dd nodes, class members and operation pointers
  clear();

  return result;
}

void MEDDLY::forwd_dfs_mt::initialize()
{
  // set up mdd operations
  mddUnion = getOperation(UNION, resF, resF, resF);
  assert(mddUnion != 0);

  mxdIntersection = getOperation(INTERSECTION, arg2F, arg2F, arg2F);
  assert(mxdIntersection != 0);

  mxdDifference = getOperation(DIFFERENCE, arg2F, arg2F, arg2F);
  assert(mxdDifference != 0);

  // Initialize the splits vector
  splits.resize(1+resF->getNumVariables(), 0);
}

void MEDDLY::forwd_dfs_mt::clear()
{
  clearSplitMxdComputeTableEntries();

  // Clear compute table for saturate()
  // (not the same as the one used by recFire)
  clearSaturateCT();

  // Clear compute table for recFire()
  removeAllComputeTableEntries();

  // Clear mdd union compute table entries
  if (mddUnion) mddUnion->removeAllComputeTableEntries();

  // clear pointer to dd nodes
  for (unsigned i = 0u; i < splits.size(); i++) arg2F->unlinkNode(splits[i]);

  // clear class members and pointers to operations
  splits.clear();
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;
}

void MEDDLY::forwd_dfs_mt::clearSplitMxdComputeTableEntries()
{
  if (mxdIntersection)
    mxdIntersection->removeAllComputeTableEntries();
  if (mxdDifference)
    mxdDifference->removeAllComputeTableEntries();
}

// Partition the nsf based on "top level"
void MEDDLY::forwd_dfs_mt::splitMxd(int mxd)
{
  MEDDLY_DCASSERT(arg2F);

  // we'll be unlinking later, so...
  arg2F->linkNode(mxd);

  // Build from top down
  for (int level = arg2F->getNumVariables(); level; level--) {

    if (0==mxd) {
      // common and easy special case
      splits[level] = 0;
      continue;
    }

    int mxdLevel = arg2F->getNodeLevel(mxd);
    MEDDLY_DCASSERT(ABS(mxdLevel <= level));

    // Initialize row reader
    expert_forest::nodeReader* Mu = isLevelAbove(level, mxdLevel)
      ? arg2F->initRedundantReader(level, mxd)
      : arg2F->initNodeReader(mxd);

    bool first = true;
    int maxDiag;

    // Read "rows"
    for (int i=0; i<Mu->getSize(); i++) {
      // Initialize column reader
      int mxdPLevel = arg2F->getNodeLevel((*Mu)[i]);
      expert_forest::nodeReader* Mp = isLevelAbove(-level, mxdPLevel)
        ? arg2F->initIdentityReader(-level, i, (*Mu)[i])
        : arg2F->initNodeReader((*Mu)[i]);

      // Intersect along the diagonal
      if (first) {
        maxDiag = arg2F->linkNode((*Mp)[i]);
        first = false;
      } else {
        int nmd = mxdIntersection->compute(maxDiag, (*Mp)[i]);
        arg2F->unlinkNode(maxDiag);
        maxDiag = nmd;
      }

      // cleanup
      arg2F->recycle(Mp);
    } // for i

    // maxDiag is what we can split from here
    splits[level] = mxdDifference->compute(mxd, maxDiag);
    arg2F->unlinkNode(mxd);
    mxd = maxDiag;

    // Cleanup
    arg2F->recycle(Mu);
  } // for level
}


int MEDDLY::forwd_dfs_mt::saturate(int mdd)
{
#ifdef DEBUG_DFS
  printf("mdd: %d\n", mdd);
#endif

  // how does saturateHelper get called?
  // bottom-up i.e. call helper on children before calling helper for parent

  MEDDLY_DCASSERT(arg1F->isReducedNode(mdd));

  // terminal condition for recursion
  if (arg1F->isTerminalNode(mdd)) return mdd;

  // search compute table
  int n = 0;
  if (findSaturateResult(mdd, n)) {
    resF->linkNode(n);
    return n;
  }

  int k = arg1F->getNodeLevel(mdd);      // level
  int sz = arg1F->getLevelSize(k);       // size

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d\n", mdd, k, sz);
#endif

  expert_forest::nodeBuilder& nb = resF->useNodeBuilder(k, sz);
  expert_forest::nodeReader* mddDptrs = arg1F->initNodeReader(mdd);
  for (int i=0; i<sz; i++) {
    nb.d(i) = (*mddDptrs)[i] ? saturate((*mddDptrs)[i]) : 0;
  }
  arg1F->recycle(mddDptrs);
  saturateHelper(nb);
  n = resF->createReducedNode(-1, nb);

  // save in compute table
  saveSaturateResult(mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  return n;
}

void MEDDLY::forwd_dfs_mt::saturateHelper(expert_forest::nodeBuilder& nb)
{
  int mxd = splits[nb.getLevel()];
  if (mxd == 0) return;

  int mxdLevel = arg2F->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd row reader, note we might skip the unprimed level
  expert_forest::nodeReader* Ru = (mxdLevel<0)
    ? arg2F->initRedundantReader(nb.getLevel(), mxd)
    : arg2F->initNodeReader(mxd);

  // indexes to explore (TBD: write a fast class for this)
  std::set<int> enabledIndexes;
  std::set<int>::iterator iter;
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.d(i)) enabledIndexes.insert(i);
  }

  // explore indexes
  while (!enabledIndexes.empty()) {
    iter = enabledIndexes.begin();
    int i = *iter;
    enabledIndexes.erase(iter);

    MEDDLY_DCASSERT(nb.d(i));
    if (0==(*Ru)[i]) continue;  // row i is empty

    // grab column (TBD: build these ahead of time?)
    int dlevel = arg2F->getNodeLevel((*Ru)[i]);

    expert_forest::nodeReader* Rp = (dlevel == -nb.getLevel())
      ? arg2F->initNodeReader((*Ru)[i])
      : arg2F->initIdentityReader(-nb.getLevel(), i, (*Ru)[i]);

    for (int j=0; j<nb.getSize(); j++) {
      if (0==(*Rp)[j]) continue;
      if (-1==nb.d(j)) continue;  // nothing can be added to this set

      int rec = recFire(nb.d(i), (*Rp)[j]);

      if (rec == 0) continue;
      if (rec == nb.d(j)) { 
        resF->unlinkNode(rec); 
        continue; 
      }

      bool updated = true;

      if (0 == nb.d(j)) {
        nb.d(j) = rec;
      }
      else if (rec == -1) {
        resF->unlinkNode(nb.d(j));
        nb.d(j) = -1;
      }
      else {
        int acc = mddUnion->compute(nb.d(j), rec);
        resF->unlinkNode(rec);
        if (acc != nb.d(j)) {
          resF->unlinkNode(nb.d(j));
          nb.d(j) = acc;
        } else {
          resF->unlinkNode(acc);
          updated = false;
        }
      }

      if (updated) {
        if (j == i) {
          // Restart inner for-loop.
          j = -1;
        } else {
          enabledIndexes.insert(j);
        }
      }

    } // for j

    // cleanup
    arg2F->recycle(Rp);
  } // while there are indexes to explore

  // cleanup
  arg2F->recycle(Ru);
}

int MEDDLY::forwd_dfs_mt::recFire(int mdd, int mxd)
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

#ifdef TRACE_RECFIRE
  printf("computing recFire(%d, %d)\n", mdd, mxd);
  printf("  node %3d ", mdd);
  arg1F->showNode(stdout, mdd, 1);
  printf("\n  node %3d ", mxd);
  arg2F->showNode(stdout, mxd, 1);
  printf("\n");
#endif

  // check if mxd and mdd are at the same level
  int mddLevel = arg1F->getNodeLevel(mdd);
  int mxdLevel = arg2F->getNodeLevel(mxd);
  int rLevel = MAX(ABS(mxdLevel), mddLevel);
  int rSize = resF->getLevelSize(rLevel);
  expert_forest::nodeBuilder& nb = resF->useNodeBuilder(rLevel, rSize);

  // Initialize mdd reader
  expert_forest::nodeReader* A = (mddLevel < rLevel)
    ? arg1F->initRedundantReader(rLevel, mdd)
    : arg1F->initNodeReader(mdd);

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i=0; i<rSize; i++) {
      nb.d(i) = recFire((*A)[i], mxd);
    }

  } else {
    // 
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (int i=0; i<rSize; i++) nb.d(i) = 0;

    // Initialize mxd reader, note we might skip the unprimed level
    expert_forest::nodeReader* Ru = (mxdLevel < 0)
      ? arg2F->initRedundantReader(rLevel, mxd)
      : arg2F->initNodeReader(mxd);

    // loop over mxd "rows"
    for (int i=0; i<rSize; i++) {
      if (0==(*A)[i])   continue; 
      if (0==(*Ru)[i])  continue; 
      expert_forest::nodeReader* Rp;
      Rp = (isLevelAbove(-rLevel, arg2F->getNodeLevel((*Ru)[i])))
        ? arg2F->initIdentityReader(rLevel, i, (*Ru)[i])
        : arg2F->initNodeReader((*Ru)[i]);

      // loop over mxd "columns"
      for (int j=0; j<rSize; j++) {
        if (0==(*Rp)[j])  continue;
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        int newstates = recFire((*A)[i], (*Rp)[j]);
        if (0==newstates) continue;
        if (0==nb.d(j)) {
          nb.d(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        int oldj = nb.d(j);
        nb.d(j) = mddUnion->compute(newstates, oldj);
        resF->unlinkNode(oldj);
        resF->unlinkNode(newstates);
      } // for j
  
      arg2F->recycle(Rp);
    } // for i

    arg2F->recycle(Ru);
  } // else

  // cleanup mdd reader
  arg1F->recycle(A);

  saturateHelper(nb);
  result = resF->createReducedNode(-1, nb);
#ifdef TRACE_ALL_OPS
  printf("computed recfire(%d, %d) = %d\n", mdd, mxd, result);
#endif
  return saveResult(mdd, mxd, result); 
}




// ******************************************************************
// *                                                                *
// *                       bckwd_dfs_mt class                       *
// *                                                                *
// ******************************************************************

#if 0

class MEDDLY::bckwd_dfs_mt : public forwd_dfs_mt {
  public:
    bckwd_dfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual ~bckwd_dfs_mt() { }
  protected:
    virtual int saturate(int mdd);
    virtual void saturateHelper(int mddLevel, std::vector<int>& mdd) {
      reverseSaturateHelper(mddLevel, mdd);
    }
    virtual void reverseSaturateHelper(int mddLevel, std::vector<int>& mdd);
    virtual int reverseRecFire(int mdd, int mxd);
};

MEDDLY::bckwd_dfs_mt::bckwd_dfs_mt(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : forwd_dfs_mt(opcode, arg1, arg2, res)
{
}

int MEDDLY::bckwd_dfs_mt::saturate(int mdd)
{
#ifdef DEBUG_DFS
  printf("mdd: %d\n", mdd);
#endif

  // how does saturateHelper get called?
  // bottom-up i.e. call helper on children before calling helper for parent

  MEDDLY_DCASSERT(resF->isReducedNode(mdd));

  // terminal condition for recursion
  if (resF->isTerminalNode(mdd)) return mdd;

  // search compute table
  int n = 0;
  if (findSaturateResult(mdd, n)) {
    resF->linkNode(n);
    return n;
  }

  int k = resF->getNodeLevel(mdd);      // level
  int sz = resF->getLevelSize(k);       // size

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d\n", mdd, k, sz);
#endif

  std::vector<int> nDptrs(sz, 0);
  assert(resF->getDownPtrs(mdd, nDptrs));
  for (std::vector<int>::iterator iter = nDptrs.begin();
      iter != nDptrs.end(); ++iter) {
    if (*iter) *iter = saturate(*iter);
  }
  saturateHelper(k, nDptrs);
  n = resF->reduceNode(resF->createTempNode(k, nDptrs));

  // save in compute table
  saveSaturateResult(mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  return n;
}


void MEDDLY::bckwd_dfs_mt::reverseSaturateHelper(int mddLevel,
    std::vector<int>& mdd)
{
  MEDDLY_DCASSERT(unsigned(resF->getLevelSize(mddLevel)) == mdd.size());

  int mxd = splits[mddLevel];
  if (arg2F->isTerminalNode(mxd)) return;

  std::vector<int> mxdDptrs;
  if (!arg2F->getDownPtrs(mxd, mxdDptrs)) return;

  // Get hold of the arrays for this level
  int levelSize = resF->getLevelSize(mddLevel);
  int* curr = scratch1[mddLevel];
  int* next = scratch2[mddLevel];
  assert(unsigned(levelSize) == mdd.size());
  for (int i = 0; i < levelSize; i++) next[i] = 1;

  bool repeat = true;
  while (repeat)
  {
    int* temp = curr;
    curr = next;
    next = temp;
    memset(next, 0, levelSize * sizeof(int));

    // For each mxd[i][j] != 0,
    //    If mdd[j] != 0 AND curr[j] == true
    //      mdd[i] += reverseRecFire(mdd[j], mxd[i][j])
    //      If mdd[i] is updated,
    //        next[i] = true

    // for each mxd[i] != 0
    for (unsigned i = 0u; i < mxdDptrs.size(); i++)
    {
      if (mxdDptrs[i] == 0) continue;
      MEDDLY_DCASSERT(!arg2F->isTerminalNode(mxdDptrs[i]));

      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(arg2F->getDownPtrs(mxdDptrs[i], mxdIDptrs));

      // For each mxd[i][j] != 0
      // mdd[i] += reverseRecFire(mdd[j], mxd[i][j])

      if (arg2F->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = arg2F->getFullNodeSize(mxdDptrs[i]);
        for (unsigned j = 0; j < mxdISize; j++)
        {
          if (mxdIDptrs[j] == 0 || mdd[j] == 0 || curr[j] == 0) continue;
          int f = reverseRecFire(mdd[j], mxdIDptrs[j]);
          if (f == 0) continue;
          int u = getMddUnion(mdd[i], f);
          resF->unlinkNode(f);
          if (u != mdd[i]) {
            // update mdd[i] and mark for next iteration
            resF->unlinkNode(mdd[i]);
            mdd[i] = u;
            next[i] = 1;
          }
          else {
            resF->unlinkNode(u);
          }
        }
      }
      else {
        // Sparse node
        unsigned mxdISize = arg2F->getSparseNodeSize(mxdDptrs[i]);
        const int* mxdIIptrs = 0;
        assert(arg2F->getSparseNodeIndexes(mxdDptrs[i], mxdIIptrs));
        for (int k = 0; k < mxdISize; k++)
        {
          MEDDLY_DCASSERT(0 != mxdIDptrs[k]);
          unsigned j = mxdIIptrs[k];
          if (mdd[j] == 0 || curr[j] == 0) continue;
          int f = reverseRecFire(mdd[j], mxdIDptrs[k]);
          if (f == 0) continue;
          int u = getMddUnion(mdd[i], f);
          resF->unlinkNode(f);
          if (u != mdd[i]) {
            // Update mdd[i] and mark for next iteration
            resF->unlinkNode(mdd[i]);
            mdd[i] = u;
            next[i] = 1;
          }
          else {
            resF->unlinkNode(u);
          }
        }
      }
    }

    // Check if the loop should repeat.
    repeat = false;
    int* nextEnd = next + levelSize;
    for (int* iter = next; iter != nextEnd; )
    {
      if (*iter++ != 0) {
        repeat = true;
        break;
      }
    }
  }
}


int MEDDLY::bckwd_dfs_mt::reverseRecFire(int mdd, int mxd)
{
  MEDDLY_DCASSERT(resF->isReducedNode(mdd));
  MEDDLY_DCASSERT(arg2F->isReducedNode(mxd));

  if (mxd == -1) {
    resF->linkNode(mdd);
    return mdd;
  }
  if (mxd == 0 || mdd == 0) return 0;

  int result = 0;
  if (findResult(mdd, mxd, result)) {
    return result;
  }

  int mxdHeight = arg2F->getNodeHeight(mxd);
  int mddHeight = resF->getNodeHeight(mdd);
  int nodeHeight = MAX(mxdHeight, mddHeight);
  int newSize = resF->getLevelSize(nodeHeight);
  std::vector<int> node(newSize, 0);

  if (mxdHeight < mddHeight) {
    std::vector<int> mddDptrs;
    resF->getDownPtrs(mdd, mddDptrs);
    for (unsigned i = 0; i < mddDptrs.size(); i++)
    {
      if (mddDptrs[i] != 0) node[i] = reverseRecFire(mddDptrs[i], mxd);
    }
  } else if (mxdHeight > mddHeight) {
    std::vector<int> mxdIDptrs;
    std::vector<int> mxdDptrs;
    arg2F->getDownPtrs(mxd, mxdDptrs);
    for (unsigned i = 0; i < mxdDptrs.size(); i++)
    {
      if (mxdDptrs[i] == 0) continue;

      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(arg2F->getDownPtrs(mxdDptrs[i], mxdIDptrs));

      if (arg2F->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = arg2F->getFullNodeSize(mxdDptrs[i]);
        for (unsigned j = 0; j < mxdISize; j++)
        {
          if (mxdIDptrs[j] != 0) continue;
          int f = reverseRecFire(mdd, mxdIDptrs[j]);
          if (f == 0) continue;
          int u = getMddUnion(node[i], f);
          resF->unlinkNode(f);
          resF->unlinkNode(node[i]);
          node[i] = u;
        }
      }
      else {
        // Sparse node
        unsigned mxdISize = arg2F->getSparseNodeSize(mxdDptrs[i]);
        for (unsigned k = 0; k < mxdISize; k++)
        {
          MEDDLY_DCASSERT(0 != mxdIDptrs[k]);
          int f = reverseRecFire(mdd, mxdIDptrs[k]);
          if (f == 0) continue;
          int u = getMddUnion(node[i], f);
          resF->unlinkNode(f);
          resF->unlinkNode(node[i]);
          node[i] = u;
        }
      }
    }
  } else {
    MEDDLY_DCASSERT(mxdHeight == mddHeight);
    std::vector<int> mddDptrs;
    std::vector<int> mxdDptrs;
    std::vector<int> mxdIDptrs;
    resF->getDownPtrs(mdd, mddDptrs);
    arg2F->getDownPtrs(mxd, mxdDptrs);

    for (unsigned i = 0; i < mxdDptrs.size(); i++)
    {
      if (mxdDptrs[i] == 0) continue;

      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(arg2F->getDownPtrs(mxdDptrs[i], mxdIDptrs));
      if (arg2F->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = arg2F->getFullNodeSize(mxdDptrs[i]);
        unsigned min = mddDptrs.size() < mxdISize? mddDptrs.size(): mxdISize;
        for (unsigned j = 0; j < min; j++)
        {
          if (mxdIDptrs[j] == 0 || mddDptrs[j] == 0) continue;
          int f = reverseRecFire(mddDptrs[j], mxdIDptrs[j]);
          if (f == 0) continue;
          int u = getMddUnion(node[i], f);
          resF->unlinkNode(f);
          resF->unlinkNode(node[i]);
          node[i] = u;
        }
      }
      else {
        // Sparse node
        unsigned mxdISize = arg2F->getSparseNodeSize(mxdDptrs[i]);
        const int* mxdIIptrs = 0;
        assert(arg2F->getSparseNodeIndexes(mxdDptrs[i], mxdIIptrs));
        for (int k = 0; k < mxdISize; k++)
        {
          MEDDLY_DCASSERT(0 != mxdIDptrs[k]);
          unsigned j = mxdIIptrs[k];
          if (j >= mddDptrs.size()) break;
          int f = reverseRecFire(mddDptrs[j], mxdIDptrs[k]);
          if (f == 0) continue;
          int u = getMddUnion(node[i], f);
          resF->unlinkNode(f);
          resF->unlinkNode(node[i]);
          node[i] = u;
        }
      }
    }
  }

  unsigned i = 0u;
  for ( ; i < node.size() && node[i] == 0; i++);
  if (i != node.size()) reverseSaturateHelper(nodeHeight, node);

  int n = resF->createTempNode(nodeHeight, node);
  result = resF->reduceNode(n);

  saveResult(mdd, mxd, result);
  return result;
}
#endif

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

  // return new bckwd_dfs_mt(this, a1, a2, r);
  throw error(error::NOT_IMPLEMENTED);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeForwardDFS(const settings &s)
{
  return new forwd_dfs_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeBackwardDFS(const settings &s)
{
  return new bckwd_dfs_opname;
}


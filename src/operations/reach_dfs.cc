
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

// #define DEBUG_RECFIRE

#define NEW_REDUCTIONS

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

#ifdef NEW_REDUCTIONS
    
    void saturateHelper(expert_forest::nodeBuilder& mdd);
    int recFire(int mdd, int mxd);

    expert_forest::nodeBuilder& recFireExpandMdd(int mdd, int mxd);
    expert_forest::nodeBuilder& recFireExpandMxd(int mdd, int mxd);
    expert_forest::nodeBuilder& recFireExpand(int mdd, int mxd);

#else

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
#endif

    bool getMxdAsVec(int mxd, int**& vec, int& size);

    // Helper for getMxdAsVec().
    // Use getMxd(0, 0, true) to clear the static arrays in getMxd().
    int** getMatrix(unsigned nodeLevel, unsigned size, bool clear);


    inline int getMddUnion(int a, int b) {
      MEDDLY_DCASSERT(resF->isTerminalNode(a) || resF->isReducedNode(a));
      MEDDLY_DCASSERT(resF->isTerminalNode(b) || resF->isReducedNode(b));
      MEDDLY_DCASSERT(mddUnion);
      return mddUnion->compute(a, b);
    };
    inline int getMxdIntersection(int a, int b) {
      MEDDLY_DCASSERT(arg2F->isTerminalNode(a) || arg2F->isReducedNode(a));
      MEDDLY_DCASSERT(arg2F->isTerminalNode(b) || arg2F->isReducedNode(b));
      MEDDLY_DCASSERT(mxdIntersection);
      return mxdIntersection->compute(a, b);
    }
    inline int getMxdDifference(int a, int b) {
      MEDDLY_DCASSERT(arg2F->isTerminalNode(a) || arg2F->isReducedNode(a));
      MEDDLY_DCASSERT(arg2F->isTerminalNode(b) || arg2F->isReducedNode(b));
      MEDDLY_DCASSERT(mxdDifference);
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

MEDDLY::forwd_dfs_mt::forwd_dfs_mt(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : common_dfs_mt(opcode, arg1, arg2, res), 
    scratch0(0), scratch1(0), scratch2(0), scratchSize(0)
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
  // set up aliases
  ed = smart_cast<expert_domain*>(resF->useDomain());
  assert(ed != 0);

  // set up mdd operation: union
  mddUnion = getOperation(UNION, resF, resF, resF);
  assert(mddUnion != 0);

  // set up mxd operations: intersection and difference
  mxdIntersection = getOperation(INTERSECTION, arg2F, arg2F, arg2F);
  assert(mxdIntersection != 0);

  mxdDifference = getOperation(DIFFERENCE, arg2F, arg2F, arg2F);
  assert(mxdDifference != 0);

  // Intialize the scratch 2-D vector (i.e. make it the correct size)
  int nLevels = ed->getNumVariables() + 1;

  // Initialize the splits vector
  splits.resize(nLevels, 0);

  // Initialize the scratch arrays to be used during saturation.
  assert(scratchSize == 0);
  assert(scratch0 == 0);
  assert(scratch1 == 0);
  assert(scratch2 == 0);
  scratchSize = nLevels;
  scratch0 = (int**) malloc(scratchSize * sizeof(int*));
  scratch1 = (int**) malloc(scratchSize * sizeof(int*));
  scratch2 = (int**) malloc(scratchSize * sizeof(int*));
  memset(scratch0, 0, scratchSize * sizeof(int*));
  memset(scratch1, 0, scratchSize * sizeof(int*));
  memset(scratch2, 0, scratchSize * sizeof(int*));
  for (unsigned height = ed->getNumVariables(); height > 0; --height)
  {
    int sz = ed->getVariableBound(height);
    assert(height < scratchSize);
    scratch0[height] = (int*) malloc(sz * sizeof(int));
    scratch1[height] = (int*) malloc(sz * sizeof(int));
    scratch2[height] = (int*) malloc(sz * sizeof(int));
    memset(scratch0[height], 0, sz * sizeof(int));
    memset(scratch1[height], 0, sz * sizeof(int));
    memset(scratch2[height], 0, sz * sizeof(int));
  }
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
  if (scratchSize > 0) {
    for (int i = 0; i < scratchSize; i++)
    {
      if (scratch0[i]) free(scratch0[i]);
      if (scratch1[i]) free(scratch1[i]);
      if (scratch2[i]) free(scratch2[i]);
    }
    if (scratch0) { free(scratch0); scratch0 = 0; }
    if (scratch1) { free(scratch1); scratch1 = 0; }
    if (scratch2) { free(scratch2); scratch2 = 0; }
    scratchSize = 0;
  }
  ed = 0;
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;

  // Clear the memory allocated within getMatrix()'s static vectors.
  getMatrix(0, 0, true);
}

void MEDDLY::forwd_dfs_mt::clearSplitMxdComputeTableEntries()
{
  if (mxdIntersection)
    mxdIntersection->removeAllComputeTableEntries();
  if (mxdDifference)
    mxdDifference->removeAllComputeTableEntries();
}

// split is used to split a mxd for the saturation algorithm
void MEDDLY::forwd_dfs_mt::splitMxd(int mxd)
{
#if 0
  arg2F->linkNode(mxd);
  splits[arg2F->getNodeLevel(mxd)] = mxd;
#else
  MEDDLY_DCASSERT(arg2F != 0);

  int falseNode = arg2F->getTerminalNode(false);
  int trueNode = arg2F->getTerminalNode(true);

  // find intersection for all mxd[i][i]
  // -- if mxd is smaller than max level size, then some mxd[i] is zero,
  //    therefore intersection is 0.
  //    -- if node is sparse, then there is some mxd[i] = 0,
  //       therefore intersection is 0.
  //    -- if node is full, if size < max level size, intersection is 0.
  // 
  // if intersection == 0, add mxd to level_mxd[level], return.
  // otherwise, 
  // -- create mxdSize nodes at primed level with a copy of corresponding
  //    mxd[i].
  // -- for each new_mxd[i][i],
  //    -- mxd[i][i] = mxd[i][i] - intersection 
  // -- set new_mxd[i] after reducing the primed level nodes
  //    note that new_mxd will never be 0 since mxd is not an identity node
  //
  // add new_mxd to level_mxd[level]
  // 
  // repeat the above for intersection
  //
  int level = 0;
  int intersection = falseNode;
  int mxdSize = 0;
  int mxdI = falseNode;

  arg2F->linkNode(mxd);

  while (!arg2F->isTerminalNode(mxd)) {
    level = arg2F->getNodeLevel(mxd);
    MEDDLY_DCASSERT(level > 0); // we only deal with unprimed levels

    // Find intersection for all mxd[i][i]
    // Note: only do this if mxd is a full node; when it is sparse, some
    // mxd[i] is 0 therefore the intersection will always be 0 (falseNode).
    intersection = falseNode;
    if (arg2F->isFullNode(mxd)) {
      mxdSize = arg2F->getFullNodeSize(mxd);
      if (mxdSize == arg2F->getLevelSize(level)) {
        // for all i, mxd[i] != 0
        intersection = trueNode;
        bool first = true;
        for (int i = 0; i < mxdSize && intersection != falseNode; ++i)
        {
          mxdI = arg2F->getFullNodeDownPtr(mxd, i);

          // If any mxd[i] is a terminal (according to Identity Reduced rules)
          // it must be node 0, and mxd[i][i] is also 0. Therefore,
          // the intersection is 0. So check for this condition, and break
          // out of the loop it true.

          // if mxdI is a terminal node it must be a 0 (falseNode)
          MEDDLY_DCASSERT((arg2F->isTerminalNode(mxdI) && mxdI == falseNode) ||
              !arg2F->isTerminalNode(mxdI));

          int mxdII = falseNode;

          if (!arg2F->isTerminalNode(mxdI)) {
            if (arg2F->isFullNode(mxdI)) {
              if (arg2F->getFullNodeSize(mxdI) > i)
                mxdII = arg2F->getFullNodeDownPtr(mxdI, i);
            } else {
              MEDDLY_DCASSERT(arg2F->isSparseNode(mxdI));
              // search for ith index
              int found = -1;
              int mxdINnz = arg2F->getSparseNodeSize(mxdI);

              if (mxdINnz > 8) {
                // binary search
                int start = 0;
                int stop = mxdINnz - 1;

                while (start < stop) {
                  int mid = (start + stop) / 2;
                  int midIndex = arg2F->getSparseNodeIndex(mxdI, mid);
                  if (midIndex < i) {
                    start = (mid == start)? mid + 1: mid;
                  } else {
                    stop = mid;
                  }
                }

                assert(start == stop);
                if (arg2F->getSparseNodeIndex(mxdI, start) == i) {
                  found = start;
                }
              }
              else {
                // linear search
                for (int j = 0; j < mxdINnz; ++j)
                {
                  if (arg2F->getSparseNodeIndex(mxdI, j) == i) {
                    found = j;
                    break;
                  }
                }
              }

              if (found != -1)
                mxdII = arg2F->getSparseNodeDownPtr(mxdI, found);
            }
          }

          if (!first) {
            int temp = getMxdIntersection(intersection, mxdII);
            arg2F->unlinkNode(intersection);
            intersection = temp;
          } else {
            first = false;
            arg2F->linkNode(mxdII);
            arg2F->unlinkNode(intersection);
            intersection = mxdII;
          }

#ifdef DEBUG_DFS
          printf("intersection: %d level: %d\n",
              intersection, arg2F->getNodeLevel(intersection));
#endif
        }
      }
    }

    MEDDLY_DCASSERT(splits[level] == falseNode);

    MEDDLY_DCASSERT(intersection == falseNode ||
        arg2F->getNodeLevel(mxd) > arg2F->getNodeLevel(intersection));

    if (intersection != falseNode) {
      splits[level] = getMxdDifference(mxd, intersection);
        // mxdDifference->compute(mxdDifferenceOp, mxd, intersection);
    } else {
      splits[level] = mxd;
      arg2F->linkNode(mxd);
    }

    // intersection becomes the mxd for the next iteration
    arg2F->unlinkNode(mxd);
    mxd = intersection;
  }

  MEDDLY_DCASSERT(arg2F->isTerminalNode(mxd));
  arg2F->unlinkNode(mxd);
#endif
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

#ifdef NEW_REDUCTIONS
  expert_forest::nodeBuilder& nb = resF->useNodeBuilder(k, sz);
  expert_forest::nodeReader* mddDptrs = arg1F->initNodeReader(mdd);
  for (int i=0; i<sz; i++) {
    nb.d(i) = (*mddDptrs)[i] ? saturate((*mddDptrs)[i]) : 0;
  }
  arg1F->recycle(mddDptrs);
  saturateHelper(nb);
  n = resF->createReducedNode(-1, nb);

#else // OLD STUFF

#ifndef ALT_SATURATE_HELPER
  std::vector<int> nDptrs(sz, 0);
  assert(arg1F->getDownPtrs(mdd, nDptrs));
  for (std::vector<int>::iterator iter = nDptrs.begin();
      iter != nDptrs.end(); ++iter) {
    if (*iter) *iter = saturate(*iter);
  }
  saturateHelper(k, nDptrs);
  n = resF->reduceNode(resF->createTempNode(k, nDptrs));
#else
  n = resF->createTempNode(k, sz, true);
  int * nDptrs = 0;
  assert(resF->getDownPtrs(n, nDptrs));
  const int* mddDptrs = 0;
  assert(arg1F->getDownPtrs(mdd, mddDptrs));

  if (arg1F->isFullNode(mdd)) {
    int mddSize = arg1F->getFullNodeSize(mdd);
    for (int i = 0; i < mddSize; ++i) {
      if (mddDptrs[i]) nDptrs[i] = saturate(mddDptrs[i]);
    }
  }
  else {
    int mddNDptrs = arg1F->getSparseNodeSize(mdd);
    const int* mddIndexes = 0;
    assert(arg1F->getSparseNodeIndexes(mdd, mddIndexes));
    for (int i = 0; i < mddNDptrs; ++i) {
      MEDDLY_DCASSERT(mddDptrs[i]);
      nDptrs[mddIndexes[i]] = saturate(mddDptrs[i]);
    }
  }

  // call saturateHelper for n
#ifdef DEBUG_DFS
  printf("Calling saturate: level %d\n", k);
#endif

  saturateHelper(n);

  // reduce and return
  n = resF->reduceNode(n);
#endif
#endif // END OF OLD STUFF

  // save in compute table
  saveSaturateResult(mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  return n;
}

#ifdef NEW_REDUCTIONS

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
        int acc = getMddUnion(nb.d(j), rec);
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
  MEDDLY_DCASSERT(resF->isReducedNode(mdd));
  MEDDLY_DCASSERT(arg2F->isReducedNode(mxd));

  if (mxd == -1) {
    resF->linkNode(mdd);
    return mdd;
  }
  if (mxd == 0 || mdd == 0) return 0;

#ifdef DEBUG_RECFIRE
  printf("calling recFire(%d, %d)\n", mdd, mxd);
#endif

  MEDDLY_DCASSERT(mdd == -1 || resF->getNodeLevel(mdd) > 0);
  MEDDLY_DCASSERT(arg2F->getNodeLevel(mxd) > 0);

  int result = 0;
  if (findResult(mdd, mxd, result)) {
    return result;
  }

  int mxdLevel = arg2F->getNodeLevel(mxd);
  int mddLevel = resF->getNodeLevel(mdd);
  expert_forest::nodeBuilder& nb = 
    (mxdLevel < mddLevel) ? recFireExpandMdd(mdd, mxd) :
    (mddLevel < mxdLevel) ? recFireExpandMxd(mdd, mxd) :
    recFireExpand(mdd, mxd);

  // TBD: check for all zeroes?
  saturateHelper(nb);
  result = resF->createReducedNode(-1, nb);

  // Save result and return it.
  saveResult(mdd, mxd, result);
  return result;
}


MEDDLY::expert_forest::nodeBuilder& 
MEDDLY::forwd_dfs_mt::recFireExpandMdd(int mdd, int mxd)
{
  // for i = 0 to levelSize-1
  //   result[i] = mdd[i] == 0? 0: recFire(mdd[i], mxd);
  MEDDLY_DCASSERT(resF->isReducedNode(mdd));
  MEDDLY_DCASSERT(!resF->isTerminalNode(mdd));

  int rLevel = resF->getNodeLevel(mdd);
  int rSize = resF->getLevelSize(rLevel);
  expert_forest::nodeBuilder& 
    result = resF->useNodeBuilder(rLevel, rSize);

  expert_forest::nodeReader* m = resF->initNodeReader(mdd);
  for (int i=0; i<m->getSize(); i++) {
    result.d(i) = ((*m)[i]) ? recFire((*m)[i], mxd) : 0;
  }
  resF->recycle(m);
  return result;
}


MEDDLY::expert_forest::nodeBuilder& 
MEDDLY::forwd_dfs_mt::recFireExpandMxd(int mdd, int mxd)
{
  // for i = 0 to levelSize-1
  //    if (mdd[i] != 0 && mxd[i] != 0)
  //      for j = 0 to levelSize-1
  //        result[j] += recFire(mdd[i], mxd[i][j])
  //
  // In this case, mdd[i] == mdd. Therefore,
  //
  // for i = 0 to levelSize-1
  //    if (mxd[i] != 0)
  //      for j = 0 to levelSize-1
  //        result[j] += recFire(mdd, mxd[i][j])
  //

  MEDDLY_DCASSERT(arg2F->isReducedNode(mxd) && !arg2F->isTerminalNode(mxd));
  MEDDLY_DCASSERT(arg2F->getNodeLevel(mxd) > 0);

  int rLevel = arg2F->getNodeLevel(mxd);
  int rSize = resF->getLevelSize(rLevel);

  expert_forest::nodeBuilder& 
    result = resF->useNodeBuilder(rLevel, rSize);
  for (int i=0; i<rSize; i++) result.d(i) = 0;

  int** mxdVec = 0;
  int mxdSize = 0;
  assert(getMxdAsVec(mxd, mxdVec, mxdSize));
  assert(mxdSize == rSize);

  for (int i = 0; i < mxdSize; ++i) {
    const int* mxdI = mxdVec[i];
    if (mxdI == 0) continue;

    for (int j = 0; j < mxdSize; ++j) {
      if (0 == mxdI[j])       continue;
      if (-1 == result.d(j))  continue;
      int rec = recFire(mdd, mxdI[j]);
      if (0 == rec)           continue;
      if (rec == result.d(j)) { 
        resF->unlinkNode(rec); 
        continue;
      }
      if (0 == result.d(j)) { 
        result.d(j) = rec; 
        continue;
      }
      if (-1 == rec) { 
        resF->unlinkNode(result.d(j)); 
        result.d(j) = -1; 
        continue;
      }
      int acc = getMddUnion(result.d(j), rec);
      resF->unlinkNode(rec);
      resF->unlinkNode(result.d(j));
      result.d(j) = acc;
    } // for j
  } // for i

  return result;
}


MEDDLY::expert_forest::nodeBuilder &
MEDDLY::forwd_dfs_mt::recFireExpand(int mdd, int mxd)
{
  // for i = 0 to levelSize-1
  //    if (mdd[i] != 0 && mxd[i] != 0)
  //      for j = 0 to levelSize-1
  //        result[j] += recFire(mdd[i], mxd[i][j])

  MEDDLY_DCASSERT(resF->isReducedNode(mdd) && !resF->isTerminalNode(mdd));
  MEDDLY_DCASSERT(arg2F->isReducedNode(mxd) && !arg2F->isTerminalNode(mxd));
  MEDDLY_DCASSERT(resF->getNodeLevel(mdd) == arg2F->getNodeLevel(mxd));
  MEDDLY_DCASSERT(resF->getNodeLevel(mdd) > 0);

  int rLevel = resF->getNodeLevel(mdd);
  int rSize = resF->getLevelSize(rLevel);
  expert_forest::nodeBuilder& 
    result = resF->useNodeBuilder(rLevel, rSize);
  for (int i=0; i<rSize; i++) result.d(i) = 0;

  expert_forest::nodeReader* m = resF->initNodeReader(mdd);

  int** mxdVec = 0;
  int mxdSize = 0;
  assert(getMxdAsVec(mxd, mxdVec, mxdSize));
  assert(mxdSize == rSize);

  // recFireExpandMxdPrime(mdd, mxd, result[])
  // At the primed level
  // for i = 0 to mxdSize-1
  //    if (mxd[i]) result[i] = UNION(result[i], recFire(mdd, mxd[i]));

  for (int i=0; i<rSize; i++) {
    if (0 == (*m)[i]) continue;
    const int* mxdI = mxdVec[i];
    if (0 == mxdI) continue;

    for (int j = 0; j < mxdSize; ++j) {
      if (0 == mxdI[j])       continue;
      if (-1 == result.d(j))  continue;
      int rec = recFire((*m)[i], mxdI[j]);
      if (0 == rec)           continue;
      if (rec == result.d(j)) { 
        resF->unlinkNode(rec); 
        continue;
      }
      if (0 == result.d(j)) { 
        result.d(j) = rec; 
        continue;
      }
      if (-1 == rec) { 
        resF->unlinkNode(result.d(j)); 
        result.d(j) = -1; 
        continue;
      }
      int acc = getMddUnion(result.d(j), rec);
      resF->unlinkNode(rec);
      resF->unlinkNode(result.d(j));
      result.d(j) = acc;
    } // for j

  } // for i


  resF->recycle(m);
  return result;
}


#else // OLD CODE...

#ifndef ALT_SATURATE_HELPER

// ----------------------------------------------------
//
//    Vector-based implementation -- works but slow
//
// ----------------------------------------------------

void MEDDLY::forwd_dfs_mt::saturateHelper(int mddLevel, std::vector<int>& mdd)
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
  for (int i = 0; i < levelSize; i) next[i++] = 1;

  bool repeat = true;
  while (repeat)
  {
    int* temp = curr;
    curr = next;
    next = temp;
    memset(next, 0, levelSize * sizeof(int));

    // for each mxd[i1 != 0
    for (unsigned i = 0u; i < mxdDptrs.size(); i++)
    {
      if (mxdDptrs[i] == 0 || mdd[i] == 0 || curr[i] == 0) continue;
      MEDDLY_DCASSERT(!arg2F->isTerminalNode(mxdDptrs[i]));

      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(arg2F->getDownPtrs(mxdDptrs[i], mxdIDptrs));

      if (arg2F->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = arg2F->getFullNodeSize(mxdDptrs[i]);
        // For each mxd[i][j] != 0
        for (unsigned j = 0; j < mxdISize; j++)
        {
          if (mxdIDptrs[j] != 0) {
            int f = recFire(mdd[i], mxdIDptrs[j]);
            if (f != 0) {
              int u = getMddUnion(mdd[j], f);
              resF->unlinkNode(f);
              if (u != mdd[j]) {
                // update mdd[j] and mark for next iteration
                resF->unlinkNode(mdd[j]);
                mdd[j] = u;
                next[j] = 1;
              }
              else {
                resF->unlinkNode(u);
              }
            }
          }
        }
      }
      else {
        // Sparse node
        unsigned mxdISize = arg2F->getSparseNodeSize(mxdDptrs[i]);
        const int* mxdIIptrs = 0;
        assert(arg2F->getSparseNodeIndexes(mxdDptrs[i], mxdIIptrs));
        // For each mxd[i][j] != 0
        for (int k = 0; k < mxdISize; k++)
        {
          MEDDLY_DCASSERT(0 != mxdIDptrs[k]);
          int f = recFire(mdd[i], mxdIDptrs[k]);
          if (f != 0) {
            unsigned j = mxdIIptrs[k];
            int u = getMddUnion(mdd[j], f);
            resF->unlinkNode(f);
            if (u != mdd[j]) {
              // Update mdd[j] and mark for next iteration
              resF->unlinkNode(mdd[j]);
              mdd[j] = u;
              next[j] = 1;
            }
            else {
              resF->unlinkNode(u);
            }
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


int MEDDLY::forwd_dfs_mt::recFire(int mdd, int mxd)
{
  MEDDLY_DCASSERT(resF->isReducedNode(mdd));
  MEDDLY_DCASSERT(arg2F->isReducedNode(mxd));

  if (mxd == -1) {
    resF->linkNode(mdd);
    return mdd;
  }
  if (mxd == 0 || mdd == 0) return 0;

  int result = 0;
  if (findResult(owner, mdd, mxd, result)) {
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
      if (mddDptrs[i] != 0) node[i] = recFire(mddDptrs[i], mxd);
    }
  } else if (mxdHeight > mddHeight) {
    std::vector<int> mxdIDptrs;
    std::vector<int> mxdDptrs;
    arg2F->getDownPtrs(mxd, mxdDptrs);
    for (unsigned i = 0; i < mxdDptrs.size(); i++)
    {
      if (mxdDptrs[i] == 0) continue;

#ifdef USE_GET_VECTOR_DOWN_POINTERS
      unsigned mxdISize = getNodeSize(arg2F, mxdDptrs[i]);
      clearVector(mxdIDptrs, mxdISize);
      if (!arg2F->getDownPtrs(mxdDptrs[i], mxdIDptrs)) continue;
      assert(mxdISize <= mxdIDptrs.size());

      for (unsigned j = 0; j < mxdISize; j++)
      {
        if (mxdIDptrs[j] == 0) continue;
        int f = recFire(mdd, mxdIDptrs[j]);
        if (f == 0) continue;
        int u = getMddUnion(node[j], f);
        resF->unlinkNode(f);
        resF->unlinkNode(node[j]);
        node[j] = u;
      }
#else
      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(arg2F->getDownPtrs(mxdDptrs[i], mxdIDptrs));

      if (arg2F->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = arg2F->getFullNodeSize(mxdDptrs[i]);
        for (unsigned j = 0; j < mxdISize; j++)
        {
          if (mxdIDptrs[j] != 0) {
            int f = recFire(mdd, mxdIDptrs[j]);
            if (f != 0) {
              int u = getMddUnion(node[j], f);
              resF->unlinkNode(f);
              resF->unlinkNode(node[j]);
              node[j] = u;
            }
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
          int f = recFire(mdd, mxdIDptrs[k]);
          if (f != 0) {
            unsigned j = mxdIIptrs[k];
            int u = getMddUnion(node[j], f);
            resF->unlinkNode(f);
            resF->unlinkNode(node[j]);
            node[j] = u;
          }
        }
      }
#endif
    }
  } else {
    MEDDLY_DCASSERT(mxdHeight == mddHeight);
    std::vector<int> mddDptrs;
    std::vector<int> mxdDptrs;
    std::vector<int> mxdIDptrs;
    resF->getDownPtrs(mdd, mddDptrs);
    arg2F->getDownPtrs(mxd, mxdDptrs);
    unsigned min = MIN(mddDptrs.size(), mxdDptrs.size());
    for (unsigned i = 0; i < min; i++)
    {
      if (mxdDptrs[i] == 0 || mddDptrs[i] == 0) continue;

#ifdef USE_GET_VECTOR_DOWN_POINTERS
      unsigned mxdISize = getNodeSize(arg2F, mxdDptrs[i]);
      clearVector(mxdIDptrs, mxdISize);
      if (!arg2F->getDownPtrs(mxdDptrs[i], mxdIDptrs)) continue;
      assert(mxdISize <= mxdIDptrs.size());

      for (unsigned j = 0; j < mxdISize; j++)
      {
        if (mxdIDptrs[j] == 0) continue;
        int f = recFire(mddDptrs[i], mxdIDptrs[j]);
        if (f == 0) continue;
        int u = getMddUnion(node[j], f);
        resF->unlinkNode(f);
        resF->unlinkNode(node[j]);
        node[j] = u;
      }
#else
      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(arg2F->getDownPtrs(mxdDptrs[i], mxdIDptrs));
      if (arg2F->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = arg2F->getFullNodeSize(mxdDptrs[i]);
        for (unsigned j = 0; j < mxdISize; j++)
        {
          if (mxdIDptrs[j] != 0) {
            int f = recFire(mddDptrs[i], mxdIDptrs[j]);
            if (f != 0) {
              int u = getMddUnion(node[j], f);
              resF->unlinkNode(f);
              resF->unlinkNode(node[j]);
              node[j] = u;
            }
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
          int f = recFire(mddDptrs[i], mxdIDptrs[k]);
          if (f != 0) {
            unsigned j = mxdIIptrs[k];
            int u = getMddUnion(node[j], f);
            resF->unlinkNode(f);
            resF->unlinkNode(node[j]);
            node[j] = u;
          }
        }
      }
#endif
    }
  }

  unsigned i = 0u;
  for ( ; i < node.size() && node[i] == 0; i++);
  if (i != node.size()) saturateHelper(nodeHeight, node);

  int n = resF->createTempNode(nodeHeight, node);
  result = resF->reduceNode(n);

  saveResult(owner, mdd, mxd, result);
  return result;
}


#else

// ----------------------------------------------------
//
//            getMxdAsVec() + sets for curr[]
//
// ----------------------------------------------------

void MEDDLY::forwd_dfs_mt::saturateHelper(int mdd)
{
  MEDDLY_DCASSERT(!resF->isReducedNode(mdd));
  int mddLevel = resF->getNodeLevel(mdd);
  int levelSize = resF->getLevelSize(mddLevel);
  MEDDLY_DCASSERT(levelSize == resF->getFullNodeSize(mdd));

  int mxd = splits[mddLevel];
  if (mxd == 0) return;

  MEDDLY_DCASSERT(arg2F->getNodeLevel(mxd) == mddLevel);

  int** mxdVec = 0;
  int mxdSize = 0;
  assert(getMxdAsVec(mxd, mxdVec, mxdSize));

  int* mddDptrs = 0;
  assert(resF->getDownPtrs(mdd, mddDptrs));
  int mddSize = resF->getFullNodeSize(mdd);

  MEDDLY_DCASSERT(mddSize == mxdSize);

  // For each mdd[i] != 0 && mxd[i][j] != 0,
  //    mdd[j] += recfire(mdd[i], mxd[i][j])

  std::set<int> enabledIndexes;
  std::set<int>::iterator iter;

  for (int i = 0; i < mxdSize; enabledIndexes.insert(i++));

  while (!enabledIndexes.empty()) {
    iter = enabledIndexes.begin();
    int i = *iter;
    enabledIndexes.erase(iter);

    const int* mddI = mddDptrs + i;
    const int* mxdI = mxdVec[i];
    if (*mddI == 0 || mxdI == 0) continue;

    for (int j = 0; j < mxdSize; ++j) {

      if (mxdI[j] == 0) continue;
      if (mddDptrs[j] == -1) continue;

      int rec = recFire(*mddI, mxdI[j]);

      if (rec == 0) continue;
      if (rec == mddDptrs[j]) { resF->unlinkNode(rec); continue; }

      bool updated = true;

      if (0 == mddDptrs[j]) {
        mddDptrs[j] = rec;
      }
      else if (rec == -1) {
        resF->unlinkNode(mddDptrs[j]);
        mddDptrs[j] = -1;
      }
      else {
        int acc = getMddUnion(mddDptrs[j], rec);
        resF->unlinkNode(rec);
        if (acc != mddDptrs[j]) {
          resF->unlinkNode(mddDptrs[j]);
          mddDptrs[j] = acc;
        } else {
          resF->unlinkNode(acc);
          updated = false;
        }
      }

      if (updated) {
        if (j == i) {
          // Re-run inner for-loop.
          j = -1;
        } else {
          enabledIndexes.insert(j);
        }
      }

    } // For each i in enabledIndexes

  } // For each mdd[i] and mxd[i]

}


void MEDDLY::forwd_dfs_mt::recFireExpandMdd(int mdd, int mxd, int& result)
{
  // for i = 0 to levelSize-1
  //   result[i] = mdd[i] == 0? 0: recFire(mdd[i], mxd);
  MEDDLY_DCASSERT(resF->isReducedNode(mdd));
  MEDDLY_DCASSERT(!resF->isTerminalNode(mdd));
  MEDDLY_DCASSERT(result == 0);

  int rLevel = resF->getNodeLevel(mdd);
  int rSize = resF->getLevelSize(rLevel);
  result = resF->createTempNode(rLevel, rSize, true);

  int* rDptrs = 0;
  assert(resF->getDownPtrs(result, rDptrs));

  const int* mddDptrs = 0;
  assert(resF->getDownPtrs(mdd, mddDptrs));

  if (resF->isFullNode(mdd)) {
    int mddSize = resF->getFullNodeSize(mdd);
    const int* mddStop = mddDptrs + mddSize;
    for ( ; mddDptrs != mddStop; ++rDptrs, ++mddDptrs) {
      if (*mddDptrs) *rDptrs = recFire(*mddDptrs, mxd);
    }
  }
  else {
    MEDDLY_DCASSERT(resF->isSparseNode(mdd));
    int mddNDptrs = resF->getSparseNodeSize(mdd);
    const int* mddStop = mddDptrs + mddNDptrs; 
    const int* mddIndexes = 0;
    assert(resF->getSparseNodeIndexes(mdd, mddIndexes));
    for ( ; mddDptrs != mddStop; ++mddIndexes, ++mddDptrs) {
      MEDDLY_DCASSERT(*mddDptrs != 0);
      rDptrs[*mddIndexes] = recFire(*mddDptrs, mxd);
    }
  }
}


void MEDDLY::forwd_dfs_mt::recFireExpandMxd(int mdd, int mxd, int& result)
{
  // for i = 0 to levelSize-1
  //    if (mdd[i] != 0 && mxd[i] != 0)
  //      for j = 0 to levelSize-1
  //        result[j] += recFire(mdd[i], mxd[i][j])
  //
  // In this case, mdd[i] == mdd. Therefore,
  //
  // for i = 0 to levelSize-1
  //    if (mxd[i] != 0)
  //      for j = 0 to levelSize-1
  //        result[j] += recFire(mdd, mxd[i][j])
  //

  MEDDLY_DCASSERT(arg2F->isReducedNode(mxd) && !arg2F->isTerminalNode(mxd));
  MEDDLY_DCASSERT(arg2F->getNodeLevel(mxd) > 0);
  MEDDLY_DCASSERT(result == 0);

  int rLevel = arg2F->getNodeLevel(mxd);
  int rSize = resF->getLevelSize(rLevel);
  result = resF->createTempNode(rLevel, rSize, true);

  int* rDptrs = 0;
  assert(resF->getDownPtrs(result, rDptrs));

  int** mxdVec = 0;
  int mxdSize = 0;
  assert(getMxdAsVec(mxd, mxdVec, mxdSize));
  assert(mxdSize == rSize);

  for (int i = 0; i < mxdSize; ++i) {
    const int* mxdI = mxdVec[i];
    if (mxdI == 0) continue;

    for (int j = 0; j < mxdSize; ++j) {
      if (mxdI[j] == 0) {}
      else if (rDptrs[j] == -1) {}
      else {
        int rec = recFire(mdd, mxdI[j]);
        if (rec == 0) {}
        else if (rec == rDptrs[j]) { resF->unlinkNode(rec); }
        else if (0 == rDptrs[j]) { rDptrs[j] = rec; }
        else if (rec == -1) { resF->unlinkNode(rDptrs[j]); rDptrs[j] = -1; }
        else {
          int acc = getMddUnion(rDptrs[j], rec);
          resF->unlinkNode(rec);
          resF->unlinkNode(rDptrs[j]);
          rDptrs[j] = acc;
        }
      }
    }
  }

}


void MEDDLY::forwd_dfs_mt::recFireExpand(int mdd, int mxd, int& result)
{
  // for i = 0 to levelSize-1
  //    if (mdd[i] != 0 && mxd[i] != 0)
  //      for j = 0 to levelSize-1
  //        result[j] += recFire(mdd[i], mxd[i][j])

  MEDDLY_DCASSERT(resF->isReducedNode(mdd) && !resF->isTerminalNode(mdd));
  MEDDLY_DCASSERT(arg2F->isReducedNode(mxd) && !arg2F->isTerminalNode(mxd));
  MEDDLY_DCASSERT(resF->getNodeLevel(mdd) == arg2F->getNodeLevel(mxd));
  MEDDLY_DCASSERT(resF->getNodeLevel(mdd) > 0);
  MEDDLY_DCASSERT(result == 0);

  int rLevel = resF->getNodeLevel(mdd);
  int rSize = resF->getLevelSize(rLevel);
  result = resF->createTempNode(rLevel, rSize, true);

  int* rDptrs = 0;
  const int* mddDptrs = 0;

  assert(resF->getDownPtrs(result, rDptrs));
  assert(resF->getDownPtrs(mdd, mddDptrs));

  int** mxdVec = 0;
  int mxdSize = 0;
  assert(getMxdAsVec(mxd, mxdVec, mxdSize));
  assert(mxdSize == rSize);

  // recFireExpandMxdPrime(mdd, mxd, result[])
  // At the primed level
  // for i = 0 to mxdSize-1
  //    if (mxd[i]) result[i] = UNION(result[i], recFire(mdd, mxd[i]));

  if (resF->isFullNode(mdd)) {
    int mddSize = resF->getFullNodeSize(mdd);
    for (int i = 0; i < mddSize; ++i) {
      const int mddI = mddDptrs[i];
      const int* mxdI = mxdVec[i];
      if (mddI == 0 || mxdI == 0) continue;

      for (int j = 0; j < mxdSize; ++j) {
        if (mxdI[j] == 0) {}
        else if (rDptrs[j] == -1) {}
        else {
          int rec = recFire(mddI, mxdI[j]);
          if (rec == 0) {}
          else if (rec == rDptrs[j]) { resF->unlinkNode(rec); }
          else if (0 == rDptrs[j]) { rDptrs[j] = rec; }
          else if (rec == -1) { resF->unlinkNode(rDptrs[j]); rDptrs[j] = -1; }
          else {
            int acc = getMddUnion(rDptrs[j], rec);
            resF->unlinkNode(rec);
            resF->unlinkNode(rDptrs[j]);
            rDptrs[j] = acc;
          }
        }
      }
    }
  } else {
    MEDDLY_DCASSERT(resF->isSparseNode(mdd));
    int mddNDptrs = resF->getSparseNodeSize(mdd);
    const int* mddIndexes = 0;
    assert(resF->getSparseNodeIndexes(mdd, mddIndexes));

    for (int i = 0; i < mddNDptrs; ++i) {
      const int mddI = mddDptrs[i];
      const int* mxdI = mxdVec[mddIndexes[i]];
      if (mddI == 0 || mxdI == 0) continue;

      for (int j = 0; j < mxdSize; ++j) {
        if (mxdI[j] == 0) {}
        else if (rDptrs[j] == -1) {}
        else {
          int rec = recFire(mddI, mxdI[j]);
          if (rec == 0) {}
          else if (rec == rDptrs[j]) { resF->unlinkNode(rec); }
          else if (0 == rDptrs[j]) { rDptrs[j] = rec; }
          else if (rec == -1) { resF->unlinkNode(rDptrs[j]); rDptrs[j] = -1; }
          else {
            int acc = getMddUnion(rDptrs[j], rec);
            resF->unlinkNode(rec);
            resF->unlinkNode(rDptrs[j]);
            rDptrs[j] = acc;
          }
        }
      }
    }
  }

}


int MEDDLY::forwd_dfs_mt::recFire(int mdd, int mxd)
{
  MEDDLY_DCASSERT(resF->isReducedNode(mdd));
  MEDDLY_DCASSERT(arg2F->isReducedNode(mxd));

  if (mxd == -1) {
    resF->linkNode(mdd);
    return mdd;
  }
  if (mxd == 0 || mdd == 0) return 0;

  MEDDLY_DCASSERT(mdd == -1 || resF->getNodeLevel(mdd) > 0);
  MEDDLY_DCASSERT(arg2F->getNodeLevel(mxd) > 0);

  int result = 0;
  if (findResult(mdd, mxd, result)) {
    return result;
  }

  int mxdLevel = arg2F->getNodeLevel(mxd);
  int mddLevel = resF->getNodeLevel(mdd);
  result = 0;
  if (mxdLevel < mddLevel) {
    recFireExpandMdd(mdd, mxd, result);
  } else if (mddLevel < mxdLevel) {
    recFireExpandMxd(mdd, mxd, result);
  } else {
    MEDDLY_DCASSERT(mxdLevel == mddLevel);
    recFireExpand(mdd, mxd, result);
  }

  // If any result[i] != 0, call saturateHelper()
  int* rDptrs = 0;
  assert(resF->getDownPtrs(result, rDptrs));
  int* rStop = rDptrs + resF->getFullNodeSize(result);
  for ( ; rDptrs != rStop && *rDptrs == 0; ++rDptrs);

  if (rDptrs == rStop) {
    // All result[i] == 0. This nodes will reduce to 0.
    // Instead of calling reduceNode(), unlink it and set result to 0.
    resF->unlinkNode(result);
    result = 0;
  } else {
    // Some result[i] != 0.
    // Saturate the node, and then reduce it.
    saturateHelper(result);
    result = resF->reduceNode(result);
  }

  // Save result and return it.
  saveResult(mdd, mxd, result);
  return result;
}

#endif
#endif // ...OLD CODE


// Returns an mxd as array of int[]
// If mxd[i] == 0, vec[i] = 0
// Otherwise, (vec[i])[j] = mxd[i][j]
// size is sizeof(int[])
bool MEDDLY::forwd_dfs_mt::getMxdAsVec(int mxd, int**& vec, int& size)
{
  assert(vec == 0);
  if (arg2F->isTerminalNode(mxd)) return false;

  MEDDLY_DCASSERT(arg2F->getNodeLevel(mxd) > 0);

  // Build the arrays
  size = arg2F->getLevelSize(arg2F->getNodeLevel(mxd));
  vec = getMatrix(arg2F->getNodeLevel(mxd), size, false);

  // Build the mxd
  const int* mxdDptrs = 0;
  assert(arg2F->getDownPtrs(mxd, mxdDptrs));

  if (arg2F->isFullNode(mxd)) {

    int mxdSize = arg2F->getFullNodeSize(mxd);
    memset(vec + mxdSize, 0, (size - mxdSize) * sizeof(int*));

    for (int i = 0; i < mxdSize; ++i) {
      int mxdI = mxdDptrs[i];

      if (mxdI == 0) {
        vec[i] = 0;
        continue;
      }

      // Build mxd[i]
      int* currVec = vec[i];
      const int* mxdIDptrs = 0;
      assert(arg2F->getDownPtrs(mxdI, mxdIDptrs));

      if (arg2F->isFullNode(mxdI)) {
        int mxdISize = arg2F->getFullNodeSize(mxdI);
        memcpy(currVec, mxdIDptrs, mxdISize * sizeof(int));
      } else {
        MEDDLY_DCASSERT(arg2F->isSparseNode(mxdI));
        int mxdINDptrs = arg2F->getSparseNodeSize(mxdI);
        const int* mxdIIndexes = 0;
        assert(arg2F->getSparseNodeIndexes(mxdI, mxdIIndexes));
        for (int j = 0 ; j < mxdINDptrs; ++j) {
          currVec[mxdIIndexes[j]] = mxdIDptrs[j];
        }

      } // End Build mxd[i]
    } // End for each mxd[i]

  } else {

    MEDDLY_DCASSERT(arg2F->isSparseNode(mxd));
    int mxdNDptrs = arg2F->getSparseNodeSize(mxd);
    const int* mxdIndexes = 0;
    assert(arg2F->getSparseNodeIndexes(mxd, mxdIndexes));

    int mxdSize = mxdIndexes[mxdNDptrs - 1] + 1;
    memset(vec + mxdSize, 0, (size - mxdSize) * sizeof(int*));

    int expectedIndex = 0;
    for (int i = 0; i < mxdNDptrs; ++i) {
      MEDDLY_DCASSERT(mxdDptrs[i] != 0);

      if (mxdIndexes[i] > expectedIndex) {
        memset(vec + expectedIndex, 0,
            (mxdIndexes[i] - expectedIndex) * sizeof(int*));
        expectedIndex = mxdIndexes[i] + 1;
      }

      // Build mxd[i]
      int* currVec = vec[mxdIndexes[i]];
      int mxdI = mxdDptrs[i];
      const int* mxdIDptrs = 0;
      assert(arg2F->getDownPtrs(mxdI, mxdIDptrs));

      if (arg2F->isFullNode(mxdI)) {
        int mxdISize = arg2F->getFullNodeSize(mxdI);
        memcpy(currVec, mxdIDptrs, mxdISize * sizeof(int));
      } else {
        MEDDLY_DCASSERT(arg2F->isSparseNode(mxdI));
        int mxdINDptrs = arg2F->getSparseNodeSize(mxdI);
        const int* mxdIIndexes = 0;
        assert(arg2F->getSparseNodeIndexes(mxdI, mxdIIndexes));
        for (int j = 0 ; j < mxdINDptrs; ++j) {
          currVec[mxdIIndexes[j]] = mxdIDptrs[j];
        }

      } // End Build mxd[i]
    } // End for each mxd[i]

  }

  return true;
}


int** MEDDLY::forwd_dfs_mt::getMatrix(unsigned nodeLevel,
    unsigned size, bool clear)
{
  // Enlarge static arrays.
  static std::vector< int** > temp;
  static std::vector< int > sizes;
  static std::vector< int** > matrices;

  if (clear) {
    for (unsigned i = 0; i < temp.size(); i++) {
      if (temp[i]) {
        for (unsigned j = 0; j < sizes[i]; j++) {
          free(temp[i][j]);
        }
        free(temp[i]);
        temp[i] = 0;
        sizes[i] = 0;
        free(matrices[i]);
        matrices[i] = 0;
      }
    }
    return 0;
  }

  // Resize temp and sizes
  if (temp.size() <= nodeLevel) {
    temp.resize(nodeLevel + 1, 0);
    sizes.resize(nodeLevel + 1, 0);
    matrices.resize(nodeLevel + 1, 0);
  }

  // Resize the matrix at temp[nodeLevel]
  if (sizes[nodeLevel] < size) {

    // Update temp
    temp[nodeLevel] =
      (int**) realloc(temp[nodeLevel], size * sizeof(int*));
    assert(temp[nodeLevel] != 0);
    memset(temp[nodeLevel] + sizes[nodeLevel], 0,
        (size - sizes[nodeLevel]) * sizeof(int*));

    for (int i = 0; i < size; i++) {
      temp[nodeLevel][i] =
        (int*) realloc(temp[nodeLevel][i], size * sizeof(int));
      assert(temp[nodeLevel][i] != 0);
    }

    // Update sizes
    sizes[nodeLevel] = size;

    // Update matrices
    matrices[nodeLevel] = 
      (int**) realloc(matrices[nodeLevel], size * sizeof(int*));
    assert(matrices[nodeLevel] != 0);
  }

  // Clear the matrix at temp[nodeLevel]
  for (int i = 0; i < size; i++) {
    memset(temp[nodeLevel][i], 0, size * sizeof(int));
  }
  memcpy(matrices[nodeLevel], temp[nodeLevel], size * sizeof(int*));

  return matrices[nodeLevel];
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


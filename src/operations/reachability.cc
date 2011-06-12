
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


#include "reachability.h"

#include <set>

//#define IGNORE_TERMS
//#define IGNORE_INCOUNT 2
//#define DEBUG_DFS

//#define USE_GET_VECTOR_DOWN_POINTERS


// **********************************************************************
//
//            Forward Reachability: Traditional Algorithm
//
// **********************************************************************


mdd_reachability_bfs* mdd_reachability_bfs::getInstance()
{
  static mdd_reachability_bfs instance;
  return &instance;
}


mdd_reachability_bfs::mdd_reachability_bfs()
{ }


mdd_reachability_bfs::~mdd_reachability_bfs() {}


int mdd_reachability_bfs::compute(op_info* owner, int mdd, int mxd)
{
  // set up aliases
  DCASSERT(owner->nParams == 3 && owner->p[0] == owner->p[2]);
  const int nOperands = 3;
  op_param plist[nOperands] = {owner->p[0], owner->p[0], owner->p[0]};
  expert_compute_manager* ecm = 
    smart_cast<expert_compute_manager*>(MEDDLY::getComputeManager());
  assert(ecm != 0);
  op_info* unionOp =
    ecm->getOpInfo(compute_manager::UNION, plist, nOperands);
  assert(unionOp != 0);
  plist[1] = owner->p[1];
  op_info* postImageOp =
    ecm->getOpInfo(compute_manager::POST_IMAGE, plist, nOperands);
  assert(postImageOp != 0);
  expert_forest* mddNm = getExpertForest(owner, 0);

  // Traditional (breadth-first) reachability analysis

#if 1
  expert_forest* mxdNm = getExpertForest(owner, 1);
  dd_edge nsf(mxdNm);
  mxdNm->linkNode(mxd);
  nsf.set(mxd, 0, mxdNm->getNodeLevel(mxd));
  dd_edge reachableStates(mddNm);
  mddNm->linkNode(mdd);
  reachableStates.set(mdd, 0, mddNm->getNodeLevel(mdd));
  dd_edge prevReachableStates(mddNm);

  while(prevReachableStates != reachableStates)
  {
    prevReachableStates = reachableStates;
    // printf("\nPost-Image (mdd:%d, mxd:%d): ",
    //    reachableStates.getNode(), nsf.getNode());
    dd_edge postImage(mddNm);
    ecm->apply(postImageOp, reachableStates, nsf, postImage);
    // printf("%d\n", postImage.getNode());
    // postImage.show(stdout, 2);
    // printf("\nUnion (mdd:%d, mdd:%d): ",
    //    reachableStates.getNode(), postImage.getNode());
    ecm->apply(unionOp, reachableStates, postImage, reachableStates);
    // printf("%d\n", reachableStates.getNode());
  }

  int result = reachableStates.getNode();
  mddNm->linkNode(result);
#else

  mdd_union* unionOpPtr = smart_cast<mdd_union*>(unionOp->op);
  DCASSERT(unionOpPtr != 0);
  mdd_post_image* postImageOpPtr =
    smart_cast<mdd_post_image*>(postImageOp->op);
  DCASSERT(postImageOpPtr != 0);

  int nsf = mxd;
  int reachableStates = mdd;
  int prevReachableStates = 0;
  int postImage = mdd;

  mddNm->linkNode(reachableStates);
  mddNm->linkNode(postImage);

  do {
    int prevPostImage = postImage;
    postImage = postImageOpPtr->compute(postImageOp, postImage, nsf);
    mddNm->unlinkNode(prevPostImage);

    prevReachableStates = reachableStates;
    reachableStates = unionOpPtr->compute(unionOp, reachableStates, postImage);
    mddNm->unlinkNode(prevReachableStates);
  } while (reachableStates != prevReachableStates);

  mddNm->unlinkNode(postImage);

  int result = reachableStates;
  // no need for linkNode(reachableStates) because that has already been
  // called once.

#endif

  return result;
}





// **********************************************************************
//
//            Forward Reachability: Saturation
//
// **********************************************************************

mdd_reachability_dfs* mdd_reachability_dfs::getInstance()
{
  static mdd_reachability_dfs instance;
  return &instance;
}


mdd_reachability_dfs::mdd_reachability_dfs()
: scratch0(0), scratch1(0), scratch2(0), scratchSize(0)
{ }


mdd_reachability_dfs::~mdd_reachability_dfs()
{
  if (scratchSize > 0) {
    for (int i = 0; i < scratchSize; i++)
    {
      if (scratch0[i]) free(scratch0[i]);
      if (scratch1[i]) free(scratch1[i]);
      if (scratch2[i]) free(scratch2[i]);
    }
    if (scratch0) free(scratch0);
    if (scratch1) free(scratch1);
    if (scratch2) free(scratch2);
  }
}


int mdd_reachability_dfs::compute(op_info* owner, int mdd, int mxd)
{
  DCASSERT(owner->nParams == 3 && owner->p[0] == owner->p[2]);

  // Initialize class members and helper operations
  initialize(owner);

  // Depth-first reachability analysis (Saturation)

#ifdef DEBUG_DFS
  printf("Consolidated Next-State Function:\n");
  xdf->showNodeGraph(stdout, mxd);
  printf("\n");

  printf("Initial State:\n");
  ddf->showNodeGraph(stdout, mdd);
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
    xdf->showNodeGraph(stdout, splits[i]);
    printf("\n");
  }

  fflush(stdout);
#endif

#ifdef DEBUG_SPLITS
  for (unsigned i = 0; i < splits.size(); i++) {
    fprintf(stdout, "# of nodes in the next-state function: %1.6e\n",
        double(xdf->getNodeCount(splits[i])));
  }

  for (unsigned i = 0; i < splits.size(); i++) {
    fprintf(stdout, "# of edges in the next-state function: %1.6e\n",
        double(xdf->getEdgeCount(splits[i], false)));
  }
#endif

  // Saturate the node
  int result = saturate(mdd);

  // clear pointers to dd nodes, class members and operation pointers
  clear();

  return result;
}


void mdd_reachability_dfs::initialize(op_info* o)
{
  // set up aliases
  owner = o;
  ecm = smart_cast<expert_compute_manager*>(MEDDLY::getComputeManager());
  assert(ecm != 0);
  ddf = getExpertForest(owner, 0);
  assert(ddf != 0);
  xdf = getExpertForest(owner, 1);
  assert(xdf != 0);
  DCASSERT(ddf->getDomain() == xdf->getDomain());
  ed = smart_cast<expert_domain*>(ddf->useDomain());
  assert(ed != 0);

  // set up mdd operation: union
  const int nOperands = 3;
  op_param plist[nOperands] = {owner->p[0], owner->p[0], owner->p[0]};

  mddUnionOp = ecm->getOpInfo(compute_manager::UNION, plist, nOperands);
  assert(mddUnionOp != 0);
  mddUnion = smart_cast<mdd_union*>(mddUnionOp->op);
  assert(mddUnion != 0);

  // set up mxd operations: intersection and difference
  plist[0] = owner->p[1]; plist[1] = owner->p[1]; plist[2] = owner->p[1];

  mxdIntersectionOp =
    ecm->getOpInfo(compute_manager::INTERSECTION, plist, nOperands);
  assert(mxdIntersectionOp != 0);
  mxdIntersection = smart_cast<mxd_intersection*>(mxdIntersectionOp->op);
  assert(mxdIntersection != 0);

  mxdDifferenceOp =
    ecm->getOpInfo(compute_manager::DIFFERENCE, plist, nOperands);
  assert(mxdDifferenceOp != 0);
  mxdDifference = smart_cast<mxd_difference*>(mxdDifferenceOp->op);
  assert(mxdDifference != 0);

  // Intialize the scratch 2-D vector (i.e. make it the correct size)
  int nLevels = ed->getTopVariable() + 1;

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
    int lh = ed->getVariableWithHeight(height);
    int sz = ed->getVariableBound(lh);
    assert(lh < scratchSize);
    scratch0[lh] = (int*) malloc(sz * sizeof(int));
    scratch1[lh] = (int*) malloc(sz * sizeof(int));
    scratch2[lh] = (int*) malloc(sz * sizeof(int));
    memset(scratch0[lh], 0, sz * sizeof(int));
    memset(scratch1[lh], 0, sz * sizeof(int));
    memset(scratch2[lh], 0, sz * sizeof(int));
  }
}


void mdd_reachability_dfs::clearSplitMxdComputeTableEntries()
{
  if (mxdIntersectionOp)
    mxdIntersectionOp->cc->removeEntries(mxdIntersectionOp);
  if (mxdDifferenceOp)
    mxdDifferenceOp->cc->removeEntries(mxdDifferenceOp);
}


void mdd_reachability_dfs::clear()
{
  clearSplitMxdComputeTableEntries();

  // Clear compute table for saturate()
  // (not the same as the one used by recFire)
  clearSaturateCT();

  // Clear compute table for recFire()
  owner->cc->removeEntries(owner);

  // Clear mdd union compute table entries
  if (mddUnion) mddUnionOp->cc->removeEntries(mddUnionOp);

  // clear pointer to dd nodes
  for (unsigned i = 0u; i < splits.size(); i++) xdf->unlinkNode(splits[i]);

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
  owner = 0;
  ecm = 0;
  ddf = 0;
  xdf = 0;
  ed = 0;
  mddUnionOp = 0;
  mxdIntersectionOp = 0;
  mxdDifferenceOp = 0;
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;

  // Clear the memory allocated within getMatrix()'s static vectors.
  getMatrix(0, 0, true);
}


// split is used to split a mxd for the saturation algorithm
void mdd_reachability_dfs::splitMxd(int mxd)
{
#if 0
  xdf->linkNode(mxd);
  splits[xdf->getNodeLevel(mxd)] = mxd;
#else
  DCASSERT(xdf != 0);

  int falseNode = xdf->getTerminalNode(false);
  int trueNode = xdf->getTerminalNode(true);

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

  xdf->linkNode(mxd);

  while (!xdf->isTerminalNode(mxd)) {
    level = xdf->getNodeLevel(mxd);
    DCASSERT(level > 0); // we only deal with unprimed levels

    // Find intersection for all mxd[i][i]
    // Note: only do this if mxd is a full node; when it is sparse, some
    // mxd[i] is 0 therefore the intersection will always be 0 (falseNode).
    intersection = falseNode;
    if (xdf->isFullNode(mxd)) {
      mxdSize = xdf->getFullNodeSize(mxd);
      if (mxdSize == xdf->getLevelSize(level)) {
        // for all i, mxd[i] != 0
        intersection = trueNode;
        bool first = true;
        for (int i = 0; i < mxdSize && intersection != falseNode; ++i)
        {
          mxdI = xdf->getFullNodeDownPtr(mxd, i);

          // If any mxd[i] is a terminal (according to Identity Reduced rules)
          // it must be node 0, and mxd[i][i] is also 0. Therefore,
          // the intersection is 0. So check for this condition, and break
          // out of the loop it true.

          // if mxdI is a terminal node it must be a 0 (falseNode)
          DCASSERT((xdf->isTerminalNode(mxdI) && mxdI == falseNode) ||
              !xdf->isTerminalNode(mxdI));

          int mxdII = falseNode;

          if (!xdf->isTerminalNode(mxdI)) {
            if (xdf->isFullNode(mxdI)) {
              if (xdf->getFullNodeSize(mxdI) > i)
                mxdII = xdf->getFullNodeDownPtr(mxdI, i);
            } else {
              DCASSERT(xdf->isSparseNode(mxdI));
              // search for ith index
              int found = -1;
              int mxdINnz = xdf->getSparseNodeSize(mxdI);

              if (mxdINnz > 8) {
                // binary search
                int start = 0;
                int stop = mxdINnz - 1;

                while (start < stop) {
                  int mid = (start + stop) / 2;
                  int midIndex = xdf->getSparseNodeIndex(mxdI, mid);
                  if (midIndex < i) {
                    start = (mid == start)? mid + 1: mid;
                  } else {
                    stop = mid;
                  }
                }

                assert(start == stop);
                if (xdf->getSparseNodeIndex(mxdI, start) == i) {
                  found = start;
                }
              }
              else {
                // linear search
                for (int j = 0; j < mxdINnz; ++j)
                {
                  if (xdf->getSparseNodeIndex(mxdI, j) == i) {
                    found = j;
                    break;
                  }
                }
              }

              if (found != -1)
                mxdII = xdf->getSparseNodeDownPtr(mxdI, found);
            }
          }

          if (!first) {
            int temp = getMxdIntersection(intersection, mxdII);
            //mxdIntersection->compute(mxdIntersectionOp, intersection, mxdII);
            xdf->unlinkNode(intersection);
            intersection = temp;
          } else {
            first = false;
            xdf->linkNode(mxdII);
            xdf->unlinkNode(intersection);
            intersection = mxdII;
          }

#ifdef DEBUG_DFS
          printf("intersection: %d level: %d\n",
              intersection, xdf->getNodeLevel(intersection));
#endif
        }
      }
    }

    DCASSERT(splits[level] == falseNode);

    DCASSERT(intersection == falseNode ||
        xdf->getNodeLevel(mxd) > xdf->getNodeLevel(intersection));

    if (intersection != falseNode) {
      splits[level] = getMxdDifference(mxd, intersection);
        // mxdDifference->compute(mxdDifferenceOp, mxd, intersection);
    } else {
      splits[level] = mxd;
      xdf->linkNode(mxd);
    }

    // intersection becomes the mxd for the next iteration
    xdf->unlinkNode(mxd);
    mxd = intersection;
  }

  DCASSERT(xdf->isTerminalNode(mxd));
  xdf->unlinkNode(mxd);
#endif
}


int mdd_reachability_dfs::saturate(int mdd)
{
#ifdef DEBUG_DFS
  printf("mdd: %d\n", mdd);
#endif

  // how does saturateHelper get called?
  // bottom-up i.e. call helper on children before calling helper for parent

  DCASSERT(ddf->isReducedNode(mdd));

  // terminal condition for recursion
  if (ddf->isTerminalNode(mdd)) return mdd;

  // search compute table
  int n = 0;
  if (findSaturateResult(mdd, n)) {
    ddf->linkNode(n);
    return n;
  }

  int k = ddf->getNodeLevel(mdd);      // level
  int sz = ddf->getLevelSize(k);       // size

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d\n", mdd, k, sz);
#endif

#ifndef ALT_SATURATE_HELPER
  std::vector<int> nDptrs(sz, 0);
  assert(ddf->getDownPtrs(mdd, nDptrs));
  for (std::vector<int>::iterator iter = nDptrs.begin();
      iter != nDptrs.end(); ++iter) {
    if (*iter) *iter = saturate(*iter);
  }
  saturateHelper(k, nDptrs);
  n = ddf->reduceNode(ddf->createTempNode(k, nDptrs));
#else
  n = ddf->createTempNode(k, sz, true);
  int * nDptrs = 0;
  assert(ddf->getDownPtrs(n, nDptrs));
  const int* mddDptrs = 0;
  assert(ddf->getDownPtrs(mdd, mddDptrs));

  if (ddf->isFullNode(mdd)) {
    int mddSize = ddf->getFullNodeSize(mdd);
    for (int i = 0; i < mddSize; ++i) {
      if (mddDptrs[i]) nDptrs[i] = saturate(mddDptrs[i]);
    }
  }
  else {
    int mddNDptrs = ddf->getSparseNodeSize(mdd);
    const int* mddIndexes = 0;
    assert(ddf->getSparseNodeIndexes(mdd, mddIndexes));
    for (int i = 0; i < mddNDptrs; ++i) {
      DCASSERT(mddDptrs[i]);
      nDptrs[mddIndexes[i]] = saturate(mddDptrs[i]);
    }
  }

  // call saturateHelper for n
#ifdef DEBUG_DFS
  printf("Calling saturate: level %d\n", k);
#endif

  saturateHelper(n);

  // reduce and return
  n = ddf->reduceNode(n);
#endif

  // save in compute table
  saveSaturateResult(mdd, n);

#ifdef DEBUG_DFS
  ddf->showNodeGraph(stdout, n);
#endif

  return n;
}


#ifndef ALT_SATURATE_HELPER

// ----------------------------------------------------
//
//    Vector-based implementation -- works but slow
//
// ----------------------------------------------------

void mdd_reachability_dfs::saturateHelper(int mddLevel, std::vector<int>& mdd)
{
  DCASSERT(unsigned(ddf->getLevelSize(mddLevel)) == mdd.size());

  int mxd = splits[mddLevel];
  if (xdf->isTerminalNode(mxd)) return;

  std::vector<int> mxdDptrs;
  if (!xdf->getDownPtrs(mxd, mxdDptrs)) return;

  // Get hold of the arrays for this level
  int levelSize = ddf->getLevelSize(mddLevel);
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
      DCASSERT(!xdf->isTerminalNode(mxdDptrs[i]));

      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs));

      if (xdf->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = xdf->getFullNodeSize(mxdDptrs[i]);
        // For each mxd[i][j] != 0
        for (unsigned j = 0; j < mxdISize; j++)
        {
          if (mxdIDptrs[j] != 0) {
            int f = recFire(mdd[i], mxdIDptrs[j]);
            if (f != 0) {
              int u = getMddUnion(mdd[j], f);
              ddf->unlinkNode(f);
              if (u != mdd[j]) {
                // update mdd[j] and mark for next iteration
                ddf->unlinkNode(mdd[j]);
                mdd[j] = u;
                next[j] = 1;
              }
              else {
                ddf->unlinkNode(u);
              }
            }
          }
        }
      }
      else {
        // Sparse node
        unsigned mxdISize = xdf->getSparseNodeSize(mxdDptrs[i]);
        const int* mxdIIptrs = 0;
        assert(xdf->getSparseNodeIndexes(mxdDptrs[i], mxdIIptrs));
        // For each mxd[i][j] != 0
        for (int k = 0; k < mxdISize; k++)
        {
          DCASSERT(0 != mxdIDptrs[k]);
          int f = recFire(mdd[i], mxdIDptrs[k]);
          if (f != 0) {
            unsigned j = mxdIIptrs[k];
            int u = getMddUnion(mdd[j], f);
            ddf->unlinkNode(f);
            if (u != mdd[j]) {
              // Update mdd[j] and mark for next iteration
              ddf->unlinkNode(mdd[j]);
              mdd[j] = u;
              next[j] = 1;
            }
            else {
              ddf->unlinkNode(u);
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


int mdd_reachability_dfs::recFire(int mdd, int mxd)
{
  DCASSERT(ddf->isReducedNode(mdd));
  DCASSERT(xdf->isReducedNode(mxd));

  if (mxd == -1) {
    ddf->linkNode(mdd);
    return mdd;
  }
  if (mxd == 0 || mdd == 0) return 0;

  int result = 0;
  if (findResult(owner, mdd, mxd, result)) {
    return result;
  }

  int mxdHeight = xdf->getNodeHeight(mxd);
  int mddHeight = ddf->getNodeHeight(mdd);
  int nodeHeight = MAX(mxdHeight, mddHeight);
  int nodeLevel = ed->getVariableWithHeight(nodeHeight);
  int newSize = ddf->getLevelSize(nodeLevel);
  std::vector<int> node(newSize, 0);

  if (mxdHeight < mddHeight) {
    std::vector<int> mddDptrs;
    ddf->getDownPtrs(mdd, mddDptrs);
    for (unsigned i = 0; i < mddDptrs.size(); i++)
    {
      if (mddDptrs[i] != 0) node[i] = recFire(mddDptrs[i], mxd);
    }
  } else if (mxdHeight > mddHeight) {
    std::vector<int> mxdIDptrs;
    std::vector<int> mxdDptrs;
    xdf->getDownPtrs(mxd, mxdDptrs);
    for (unsigned i = 0; i < mxdDptrs.size(); i++)
    {
      if (mxdDptrs[i] == 0) continue;

#ifdef USE_GET_VECTOR_DOWN_POINTERS
      unsigned mxdISize = getNodeSize(xdf, mxdDptrs[i]);
      clearVector(mxdIDptrs, mxdISize);
      if (!xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs)) continue;
      assert(mxdISize <= mxdIDptrs.size());

      for (unsigned j = 0; j < mxdISize; j++)
      {
        if (mxdIDptrs[j] == 0) continue;
        int f = recFire(mdd, mxdIDptrs[j]);
        if (f == 0) continue;
        int u = getMddUnion(node[j], f);
        ddf->unlinkNode(f);
        ddf->unlinkNode(node[j]);
        node[j] = u;
      }
#else
      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs));

      if (xdf->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = xdf->getFullNodeSize(mxdDptrs[i]);
        for (unsigned j = 0; j < mxdISize; j++)
        {
          if (mxdIDptrs[j] != 0) {
            int f = recFire(mdd, mxdIDptrs[j]);
            if (f != 0) {
              int u = getMddUnion(node[j], f);
              ddf->unlinkNode(f);
              ddf->unlinkNode(node[j]);
              node[j] = u;
            }
          }
        }
      }
      else {
        // Sparse node
        unsigned mxdISize = xdf->getSparseNodeSize(mxdDptrs[i]);
        const int* mxdIIptrs = 0;
        assert(xdf->getSparseNodeIndexes(mxdDptrs[i], mxdIIptrs));
        for (int k = 0; k < mxdISize; k++)
        {
          DCASSERT(0 != mxdIDptrs[k]);
          int f = recFire(mdd, mxdIDptrs[k]);
          if (f != 0) {
            unsigned j = mxdIIptrs[k];
            int u = getMddUnion(node[j], f);
            ddf->unlinkNode(f);
            ddf->unlinkNode(node[j]);
            node[j] = u;
          }
        }
      }
#endif
    }
  } else {
    DCASSERT(mxdHeight == mddHeight);
    std::vector<int> mddDptrs;
    std::vector<int> mxdDptrs;
    std::vector<int> mxdIDptrs;
    ddf->getDownPtrs(mdd, mddDptrs);
    xdf->getDownPtrs(mxd, mxdDptrs);
    unsigned min = MIN(mddDptrs.size(), mxdDptrs.size());
    for (unsigned i = 0; i < min; i++)
    {
      if (mxdDptrs[i] == 0 || mddDptrs[i] == 0) continue;

#ifdef USE_GET_VECTOR_DOWN_POINTERS
      unsigned mxdISize = getNodeSize(xdf, mxdDptrs[i]);
      clearVector(mxdIDptrs, mxdISize);
      if (!xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs)) continue;
      assert(mxdISize <= mxdIDptrs.size());

      for (unsigned j = 0; j < mxdISize; j++)
      {
        if (mxdIDptrs[j] == 0) continue;
        int f = recFire(mddDptrs[i], mxdIDptrs[j]);
        if (f == 0) continue;
        int u = getMddUnion(node[j], f);
        ddf->unlinkNode(f);
        ddf->unlinkNode(node[j]);
        node[j] = u;
      }
#else
      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs));
      if (xdf->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = xdf->getFullNodeSize(mxdDptrs[i]);
        for (unsigned j = 0; j < mxdISize; j++)
        {
          if (mxdIDptrs[j] != 0) {
            int f = recFire(mddDptrs[i], mxdIDptrs[j]);
            if (f != 0) {
              int u = getMddUnion(node[j], f);
              ddf->unlinkNode(f);
              ddf->unlinkNode(node[j]);
              node[j] = u;
            }
          }
        }
      }
      else {
        // Sparse node
        unsigned mxdISize = xdf->getSparseNodeSize(mxdDptrs[i]);
        const int* mxdIIptrs = 0;
        assert(xdf->getSparseNodeIndexes(mxdDptrs[i], mxdIIptrs));
        for (int k = 0; k < mxdISize; k++)
        {
          DCASSERT(0 != mxdIDptrs[k]);
          int f = recFire(mddDptrs[i], mxdIDptrs[k]);
          if (f != 0) {
            unsigned j = mxdIIptrs[k];
            int u = getMddUnion(node[j], f);
            ddf->unlinkNode(f);
            ddf->unlinkNode(node[j]);
            node[j] = u;
          }
        }
      }
#endif
    }
  }

  unsigned i = 0u;
  for ( ; i < node.size() && node[i] == 0; i++);
  if (i != node.size()) saturateHelper(nodeLevel, node);

  int n = ddf->createTempNode(nodeLevel, node);
  result = ddf->reduceNode(n);

  saveResult(owner, mdd, mxd, result);
  return result;
}


#else

// ----------------------------------------------------
//
//            getMxdAsVec() + sets for curr[]
//
// ----------------------------------------------------

void mdd_reachability_dfs::saturateHelper(int mdd)
{
  DCASSERT(!ddf->isReducedNode(mdd));
  int mddLevel = ddf->getNodeLevel(mdd);
  int levelSize = ddf->getLevelSize(mddLevel);
  DCASSERT(levelSize == ddf->getFullNodeSize(mdd));

  int mxd = splits[mddLevel];
  if (mxd == 0) return;

  DCASSERT(xdf->getNodeLevel(mxd) == mddLevel);

  int** mxdVec = 0;
  int mxdSize = 0;
  assert(getMxdAsVec(mxd, mxdVec, mxdSize));

  int* mddDptrs = 0;
  assert(ddf->getDownPtrs(mdd, mddDptrs));
  int mddSize = ddf->getFullNodeSize(mdd);

  DCASSERT(mddSize == mxdSize);

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
      if (rec == mddDptrs[j]) { ddf->unlinkNode(rec); continue; }

      bool updated = true;

      if (0 == mddDptrs[j]) {
        mddDptrs[j] = rec;
      }
      else if (rec == -1) {
        ddf->unlinkNode(mddDptrs[j]);
        mddDptrs[j] = -1;
      }
      else {
        int acc = getMddUnion(mddDptrs[j], rec);
        ddf->unlinkNode(rec);
        if (acc != mddDptrs[j]) {
          ddf->unlinkNode(mddDptrs[j]);
          mddDptrs[j] = acc;
        } else {
          ddf->unlinkNode(acc);
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


void mdd_reachability_dfs::recFireExpandMdd(int mdd, int mxd, int& result)
{
  // for i = 0 to levelSize-1
  //   result[i] = mdd[i] == 0? 0: recFire(mdd[i], mxd);
  DCASSERT(ddf->isReducedNode(mdd));
  DCASSERT(!ddf->isTerminalNode(mdd));
  DCASSERT(result == 0);

  int rLevel = ddf->getNodeLevel(mdd);
  int rSize = ddf->getLevelSize(rLevel);
  result = ddf->createTempNode(rLevel, rSize, true);

  int* rDptrs = 0;
  assert(ddf->getDownPtrs(result, rDptrs));

  const int* mddDptrs = 0;
  assert(ddf->getDownPtrs(mdd, mddDptrs));

  if (ddf->isFullNode(mdd)) {
    int mddSize = ddf->getFullNodeSize(mdd);
    const int* mddStop = mddDptrs + mddSize;
    for ( ; mddDptrs != mddStop; ++rDptrs, ++mddDptrs) {
      if (*mddDptrs) *rDptrs = recFire(*mddDptrs, mxd);
    }
  }
  else {
    DCASSERT(ddf->isSparseNode(mdd));
    int mddNDptrs = ddf->getSparseNodeSize(mdd);
    const int* mddStop = mddDptrs + mddNDptrs; 
    const int* mddIndexes = 0;
    assert(ddf->getSparseNodeIndexes(mdd, mddIndexes));
    for ( ; mddDptrs != mddStop; ++mddIndexes, ++mddDptrs) {
      DCASSERT(*mddDptrs != 0);
      rDptrs[*mddIndexes] = recFire(*mddDptrs, mxd);
    }
  }
}


void mdd_reachability_dfs::recFireExpandMxd(int mdd, int mxd, int& result)
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

  DCASSERT(xdf->isReducedNode(mxd) && !xdf->isTerminalNode(mxd));
  DCASSERT(xdf->getNodeLevel(mxd) > 0);
  DCASSERT(result == 0);

  int rLevel = xdf->getNodeLevel(mxd);
  int rSize = ddf->getLevelSize(rLevel);
  result = ddf->createTempNode(rLevel, rSize, true);

  int* rDptrs = 0;
  assert(ddf->getDownPtrs(result, rDptrs));

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
        else if (rec == rDptrs[j]) { ddf->unlinkNode(rec); }
        else if (0 == rDptrs[j]) { rDptrs[j] = rec; }
        else if (rec == -1) { ddf->unlinkNode(rDptrs[j]); rDptrs[j] = -1; }
        else {
          int acc = getMddUnion(rDptrs[j], rec);
          ddf->unlinkNode(rec);
          ddf->unlinkNode(rDptrs[j]);
          rDptrs[j] = acc;
        }
      }
    }
  }

}


void mdd_reachability_dfs::recFireExpand(int mdd, int mxd, int& result)
{
  // for i = 0 to levelSize-1
  //    if (mdd[i] != 0 && mxd[i] != 0)
  //      for j = 0 to levelSize-1
  //        result[j] += recFire(mdd[i], mxd[i][j])

  DCASSERT(ddf->isReducedNode(mdd) && !ddf->isTerminalNode(mdd));
  DCASSERT(xdf->isReducedNode(mxd) && !xdf->isTerminalNode(mxd));
  DCASSERT(ddf->getNodeLevel(mdd) == xdf->getNodeLevel(mxd));
  DCASSERT(ddf->getNodeLevel(mdd) > 0);
  DCASSERT(result == 0);

  int rLevel = ddf->getNodeLevel(mdd);
  int rSize = ddf->getLevelSize(rLevel);
  result = ddf->createTempNode(rLevel, rSize, true);

  int* rDptrs = 0;
  const int* mddDptrs = 0;

  assert(ddf->getDownPtrs(result, rDptrs));
  assert(ddf->getDownPtrs(mdd, mddDptrs));

  int** mxdVec = 0;
  int mxdSize = 0;
  assert(getMxdAsVec(mxd, mxdVec, mxdSize));
  assert(mxdSize == rSize);

  // recFireExpandMxdPrime(mdd, mxd, result[])
  // At the primed level
  // for i = 0 to mxdSize-1
  //    if (mxd[i]) result[i] = UNION(result[i], recFire(mdd, mxd[i]));

  if (ddf->isFullNode(mdd)) {
    int mddSize = ddf->getFullNodeSize(mdd);
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
          else if (rec == rDptrs[j]) { ddf->unlinkNode(rec); }
          else if (0 == rDptrs[j]) { rDptrs[j] = rec; }
          else if (rec == -1) { ddf->unlinkNode(rDptrs[j]); rDptrs[j] = -1; }
          else {
            int acc = getMddUnion(rDptrs[j], rec);
            ddf->unlinkNode(rec);
            ddf->unlinkNode(rDptrs[j]);
            rDptrs[j] = acc;
          }
        }
      }
    }
  } else {
    DCASSERT(ddf->isSparseNode(mdd));
    int mddNDptrs = ddf->getSparseNodeSize(mdd);
    const int* mddIndexes = 0;
    assert(ddf->getSparseNodeIndexes(mdd, mddIndexes));

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
          else if (rec == rDptrs[j]) { ddf->unlinkNode(rec); }
          else if (0 == rDptrs[j]) { rDptrs[j] = rec; }
          else if (rec == -1) { ddf->unlinkNode(rDptrs[j]); rDptrs[j] = -1; }
          else {
            int acc = getMddUnion(rDptrs[j], rec);
            ddf->unlinkNode(rec);
            ddf->unlinkNode(rDptrs[j]);
            rDptrs[j] = acc;
          }
        }
      }
    }
  }

}


int mdd_reachability_dfs::recFire(int mdd, int mxd)
{
  DCASSERT(ddf->isReducedNode(mdd));
  DCASSERT(xdf->isReducedNode(mxd));

  if (mxd == -1) {
    ddf->linkNode(mdd);
    return mdd;
  }
  if (mxd == 0 || mdd == 0) return 0;

  DCASSERT(mdd == -1 || ddf->getNodeLevel(mdd) > 0);
  DCASSERT(xdf->getNodeLevel(mxd) > 0);

  int result = 0;
  if (findResult(owner, mdd, mxd, result)) {
    return result;
  }

  int mxdLevel = xdf->getNodeLevel(mxd);
  int mddLevel = ddf->getNodeLevel(mdd);
  result = 0;
  if (mxdLevel < mddLevel) {
    recFireExpandMdd(mdd, mxd, result);
  } else if (mddLevel < mxdLevel) {
    recFireExpandMxd(mdd, mxd, result);
  } else {
    DCASSERT(mxdLevel == mddLevel);
    recFireExpand(mdd, mxd, result);
  }

  // If any result[i] != 0, call saturateHelper()
  int* rDptrs = 0;
  assert(ddf->getDownPtrs(result, rDptrs));
  int* rStop = rDptrs + ddf->getFullNodeSize(result);
  for ( ; rDptrs != rStop && *rDptrs == 0; ++rDptrs);

  if (rDptrs == rStop) {
    // All result[i] == 0. This nodes will reduce to 0.
    // Instead of calling reduceNode(), unlink it and set result to 0.
    ddf->unlinkNode(result);
    result = 0;
  } else {
    // Some result[i] != 0.
    // Saturate the node, and then reduce it.
    saturateHelper(result);
    result = ddf->reduceNode(result);
  }

  // Save result and return it.
  saveResult(owner, mdd, mxd, result);
  return result;
}


// Returns an mxd as array of int[]
// If mxd[i] == 0, vec[i] = 0
// Otherwise, (vec[i])[j] = mxd[i][j]
// size is sizeof(int[])
bool mdd_reachability_dfs::getMxdAsVec(int mxd, int**& vec, int& size)
{
  assert(vec == 0);
  if (xdf->isTerminalNode(mxd)) return false;

  DCASSERT(xdf->getNodeLevel(mxd) > 0);

  // Build the arrays
  size = xdf->getLevelSize(xdf->getNodeLevel(mxd));
  vec = getMatrix(xdf->getNodeLevel(mxd), size, false);

  // Build the mxd
  const int* mxdDptrs = 0;
  assert(xdf->getDownPtrs(mxd, mxdDptrs));

  if (xdf->isFullNode(mxd)) {

    int mxdSize = xdf->getFullNodeSize(mxd);
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
      assert(xdf->getDownPtrs(mxdI, mxdIDptrs));

      if (xdf->isFullNode(mxdI)) {
        int mxdISize = xdf->getFullNodeSize(mxdI);
        memcpy(currVec, mxdIDptrs, mxdISize * sizeof(int));
      } else {
        DCASSERT(xdf->isSparseNode(mxdI));
        int mxdINDptrs = xdf->getSparseNodeSize(mxdI);
        const int* mxdIIndexes = 0;
        assert(xdf->getSparseNodeIndexes(mxdI, mxdIIndexes));
        for (int j = 0 ; j < mxdINDptrs; ++j) {
          currVec[mxdIIndexes[j]] = mxdIDptrs[j];
        }

      } // End Build mxd[i]
    } // End for each mxd[i]

  } else {

    DCASSERT(xdf->isSparseNode(mxd));
    int mxdNDptrs = xdf->getSparseNodeSize(mxd);
    const int* mxdIndexes = 0;
    assert(xdf->getSparseNodeIndexes(mxd, mxdIndexes));

    int mxdSize = mxdIndexes[mxdNDptrs - 1] + 1;
    memset(vec + mxdSize, 0, (size - mxdSize) * sizeof(int*));

    int expectedIndex = 0;
    for (int i = 0; i < mxdNDptrs; ++i) {
      DCASSERT(mxdDptrs[i] != 0);

      if (mxdIndexes[i] > expectedIndex) {
        memset(vec + expectedIndex, 0,
            (mxdIndexes[i] - expectedIndex) * sizeof(int*));
        expectedIndex = mxdIndexes[i] + 1;
      }

      // Build mxd[i]
      int* currVec = vec[mxdIndexes[i]];
      int mxdI = mxdDptrs[i];
      const int* mxdIDptrs = 0;
      assert(xdf->getDownPtrs(mxdI, mxdIDptrs));

      if (xdf->isFullNode(mxdI)) {
        int mxdISize = xdf->getFullNodeSize(mxdI);
        memcpy(currVec, mxdIDptrs, mxdISize * sizeof(int));
      } else {
        DCASSERT(xdf->isSparseNode(mxdI));
        int mxdINDptrs = xdf->getSparseNodeSize(mxdI);
        const int* mxdIIndexes = 0;
        assert(xdf->getSparseNodeIndexes(mxdI, mxdIIndexes));
        for (int j = 0 ; j < mxdINDptrs; ++j) {
          currVec[mxdIIndexes[j]] = mxdIDptrs[j];
        }

      } // End Build mxd[i]
    } // End for each mxd[i]

  }

  return true;
}

#endif


int** mdd_reachability_dfs::getMatrix(unsigned nodeLevel,
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


bool
mdd_reachability_dfs::
findResult(op_info* owner, int a, int b, int& c)
{
#ifndef USE_BINARY_COMPUTE_CACHE
  static int key[2];

  // create cache entry
  key[0] = a; key[1] = b;

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  c = cacheEntry[2];
  getExpertForest(owner, 2)->linkNode(c);
  return true;

#else

  if (!owner->cc->find(owner, a, b, c)) return false;
  getExpertForest(owner, 2)->linkNode(c);
  return true;

#endif
}


void
mdd_reachability_dfs::
saveResult(op_info* owner, int a, int b, int c)
{
  // a and c belong to mdd, b belongs to mxd.

#ifndef USE_BINARY_COMPUTE_CACHE
  static int cacheEntry[3];

  // create cache entry
  cacheEntry[0] = a; cacheEntry[1] = b; cacheEntry[2] = c;

  getExpertForest(owner, 0)->cacheNode(a);
  getExpertForest(owner, 1)->cacheNode(b);
  getExpertForest(owner, 2)->cacheNode(c);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));

#else
  owner->cc->add(owner, a, b, c);
#endif
}




// **********************************************************************
//
//            Backward Reachability: Traditional Algorithm
//
// **********************************************************************

mdd_backward_reachability_bfs* mdd_backward_reachability_bfs::getInstance()
{
  static mdd_backward_reachability_bfs instance;
  return &instance;
}


int mdd_backward_reachability_bfs::compute(op_info* owner, int mdd, int mxd)
{
  // set up aliases
  DCASSERT(owner->nParams == 3 && owner->p[0] == owner->p[2]);
  const int nOperands = 3;
  op_param plist[nOperands] = {owner->p[0], owner->p[0], owner->p[0]};
  expert_compute_manager* ecm = 
    smart_cast<expert_compute_manager*>(MEDDLY::getComputeManager());
  assert(ecm != 0);
  op_info* unionOp =
    ecm->getOpInfo(compute_manager::UNION, plist, nOperands);
  assert(unionOp != 0);
  plist[1] = owner->p[1];
  op_info* preImageOp =
    ecm->getOpInfo(compute_manager::PRE_IMAGE, plist, nOperands);
  assert(preImageOp != 0);
  expert_forest* mddNm = getExpertForest(owner, 0);

  // Traditional (breadth-first) reachability analysis

#if 1
  expert_forest* mxdNm = getExpertForest(owner, 1);
  dd_edge nsf(mxdNm);
  mxdNm->linkNode(mxd);
  nsf.set(mxd, 0, mxdNm->getNodeLevel(mxd));
  dd_edge reachableStates(mddNm);
  mddNm->linkNode(mdd);
  reachableStates.set(mdd, 0, mddNm->getNodeLevel(mdd));
  dd_edge prevReachableStates(mddNm);

  while(prevReachableStates != reachableStates)
  {
    prevReachableStates = reachableStates;
#ifdef DEBUG_REV_SAT_BFS
    printf("\nPre-Image (mdd:%d, mxd:%d): ",
        reachableStates.getNode(), nsf.getNode());
#endif
    dd_edge preImage(mddNm);
    ecm->apply(preImageOp, reachableStates, nsf, preImage);
#ifdef DEBUG_REV_SAT_BFS
    printf("%d\n", preImage.getNode());
    preImage.show(stdout, 2);
    printf("\nUnion (mdd:%d, mdd:%d): ",
        reachableStates.getNode(), preImage.getNode());
#endif
    ecm->apply(unionOp, reachableStates, preImage, reachableStates);
#ifdef DEBUG_REV_SAT_BFS
    printf("%d\n", reachableStates.getNode());
#endif
  }

  int result = reachableStates.getNode();
  mddNm->linkNode(result);
#else

  mdd_union* unionOpPtr = smart_cast<mdd_union*>(unionOp->op);
  DCASSERT(unionOpPtr != 0);
  mdd_pre_image* preImageOpPtr =
    smart_cast<mdd_pre_image*>(preImageOp->op);
  DCASSERT(preImageOpPtr != 0);

  int nsf = mxd;
  int reachableStates = mdd;
  int prevReachableStates = 0;
  int preImage = mdd;

  mddNm->linkNode(reachableStates);
  mddNm->linkNode(preImage);

  do {
    int prevPreImage = preImage;
    preImage = preImageOpPtr->compute(preImageOp, preImage, nsf);
    mddNm->unlinkNode(prevPreImage);

    prevReachableStates = reachableStates;
    reachableStates = unionOpPtr->compute(unionOp, reachableStates, preImage);
    mddNm->unlinkNode(prevReachableStates);
  } while (reachableStates != prevReachableStates);

  mddNm->unlinkNode(preImage);

  int result = reachableStates;
  // no need for linkNode(reachableStates) because that has already been
  // called once.

#endif

  return result;
}





// **********************************************************************
//
//            Backward Reachability: Saturation
//
// **********************************************************************


mdd_backward_reachability_dfs* mdd_backward_reachability_dfs::getInstance()
{
  static mdd_backward_reachability_dfs instance;
  return &instance;
}


int mdd_backward_reachability_dfs::saturate(int mdd)
{
#ifdef DEBUG_DFS
  printf("mdd: %d\n", mdd);
#endif

  // how does saturateHelper get called?
  // bottom-up i.e. call helper on children before calling helper for parent

  DCASSERT(ddf->isReducedNode(mdd));

  // terminal condition for recursion
  if (ddf->isTerminalNode(mdd)) return mdd;

  // search compute table
  int n = 0;
  if (findSaturateResult(mdd, n)) {
    ddf->linkNode(n);
    return n;
  }

  int k = ddf->getNodeLevel(mdd);      // level
  int sz = ddf->getLevelSize(k);       // size

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d\n", mdd, k, sz);
#endif

  std::vector<int> nDptrs(sz, 0);
  assert(ddf->getDownPtrs(mdd, nDptrs));
  for (std::vector<int>::iterator iter = nDptrs.begin();
      iter != nDptrs.end(); ++iter) {
    if (*iter) *iter = saturate(*iter);
  }
  saturateHelper(k, nDptrs);
  n = ddf->reduceNode(ddf->createTempNode(k, nDptrs));

  // save in compute table
  saveSaturateResult(mdd, n);

#ifdef DEBUG_DFS
  ddf->showNodeGraph(stdout, n);
#endif

  return n;
}


void mdd_backward_reachability_dfs::reverseSaturateHelper(int mddLevel,
    std::vector<int>& mdd)
{
  DCASSERT(unsigned(ddf->getLevelSize(mddLevel)) == mdd.size());

  int mxd = splits[mddLevel];
  if (xdf->isTerminalNode(mxd)) return;

  std::vector<int> mxdDptrs;
  if (!xdf->getDownPtrs(mxd, mxdDptrs)) return;

  // Get hold of the arrays for this level
  int levelSize = ddf->getLevelSize(mddLevel);
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
      DCASSERT(!xdf->isTerminalNode(mxdDptrs[i]));

      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs));

      // For each mxd[i][j] != 0
      // mdd[i] += reverseRecFire(mdd[j], mxd[i][j])

      if (xdf->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = xdf->getFullNodeSize(mxdDptrs[i]);
        for (unsigned j = 0; j < mxdISize; j++)
        {
          if (mxdIDptrs[j] == 0 || mdd[j] == 0 || curr[j] == 0) continue;
          int f = reverseRecFire(mdd[j], mxdIDptrs[j]);
          if (f == 0) continue;
          int u = getMddUnion(mdd[i], f);
          ddf->unlinkNode(f);
          if (u != mdd[i]) {
            // update mdd[i] and mark for next iteration
            ddf->unlinkNode(mdd[i]);
            mdd[i] = u;
            next[i] = 1;
          }
          else {
            ddf->unlinkNode(u);
          }
        }
      }
      else {
        // Sparse node
        unsigned mxdISize = xdf->getSparseNodeSize(mxdDptrs[i]);
        const int* mxdIIptrs = 0;
        assert(xdf->getSparseNodeIndexes(mxdDptrs[i], mxdIIptrs));
        for (int k = 0; k < mxdISize; k++)
        {
          DCASSERT(0 != mxdIDptrs[k]);
          unsigned j = mxdIIptrs[k];
          if (mdd[j] == 0 || curr[j] == 0) continue;
          int f = reverseRecFire(mdd[j], mxdIDptrs[k]);
          if (f == 0) continue;
          int u = getMddUnion(mdd[i], f);
          ddf->unlinkNode(f);
          if (u != mdd[i]) {
            // Update mdd[i] and mark for next iteration
            ddf->unlinkNode(mdd[i]);
            mdd[i] = u;
            next[i] = 1;
          }
          else {
            ddf->unlinkNode(u);
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


int mdd_backward_reachability_dfs::reverseRecFire(int mdd, int mxd)
{
  DCASSERT(ddf->isReducedNode(mdd));
  DCASSERT(xdf->isReducedNode(mxd));

  if (mxd == -1) {
    ddf->linkNode(mdd);
    return mdd;
  }
  if (mxd == 0 || mdd == 0) return 0;

  int result = 0;
  if (findResult(owner, mdd, mxd, result)) {
    return result;
  }

  int mxdHeight = xdf->getNodeHeight(mxd);
  int mddHeight = ddf->getNodeHeight(mdd);
  int nodeHeight = MAX(mxdHeight, mddHeight);
  int nodeLevel = ed->getVariableWithHeight(nodeHeight);
  int newSize = ddf->getLevelSize(nodeLevel);
  std::vector<int> node(newSize, 0);

  if (mxdHeight < mddHeight) {
    std::vector<int> mddDptrs;
    ddf->getDownPtrs(mdd, mddDptrs);
    for (unsigned i = 0; i < mddDptrs.size(); i++)
    {
      if (mddDptrs[i] != 0) node[i] = reverseRecFire(mddDptrs[i], mxd);
    }
  } else if (mxdHeight > mddHeight) {
    std::vector<int> mxdIDptrs;
    std::vector<int> mxdDptrs;
    xdf->getDownPtrs(mxd, mxdDptrs);
    for (unsigned i = 0; i < mxdDptrs.size(); i++)
    {
      if (mxdDptrs[i] == 0) continue;

      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs));

      if (xdf->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = xdf->getFullNodeSize(mxdDptrs[i]);
        for (unsigned j = 0; j < mxdISize; j++)
        {
          if (mxdIDptrs[j] != 0) continue;
          int f = reverseRecFire(mdd, mxdIDptrs[j]);
          if (f == 0) continue;
          int u = getMddUnion(node[i], f);
          ddf->unlinkNode(f);
          ddf->unlinkNode(node[i]);
          node[i] = u;
        }
      }
      else {
        // Sparse node
        unsigned mxdISize = xdf->getSparseNodeSize(mxdDptrs[i]);
        for (unsigned k = 0; k < mxdISize; k++)
        {
          DCASSERT(0 != mxdIDptrs[k]);
          int f = reverseRecFire(mdd, mxdIDptrs[k]);
          if (f == 0) continue;
          int u = getMddUnion(node[i], f);
          ddf->unlinkNode(f);
          ddf->unlinkNode(node[i]);
          node[i] = u;
        }
      }
    }
  } else {
    DCASSERT(mxdHeight == mddHeight);
    std::vector<int> mddDptrs;
    std::vector<int> mxdDptrs;
    std::vector<int> mxdIDptrs;
    ddf->getDownPtrs(mdd, mddDptrs);
    xdf->getDownPtrs(mxd, mxdDptrs);

    for (unsigned i = 0; i < mxdDptrs.size(); i++)
    {
      if (mxdDptrs[i] == 0) continue;

      // Breaking it up into two cases (a) mxd[i] is full, or (b) sparse.
      const int* mxdIDptrs = 0;
      assert(xdf->getDownPtrs(mxdDptrs[i], mxdIDptrs));
      if (xdf->isFullNode(mxdDptrs[i])) {
        // Full node
        unsigned mxdISize = xdf->getFullNodeSize(mxdDptrs[i]);
        unsigned min = mddDptrs.size() < mxdISize? mddDptrs.size(): mxdISize;
        for (unsigned j = 0; j < min; j++)
        {
          if (mxdIDptrs[j] == 0 || mddDptrs[j] == 0) continue;
          int f = reverseRecFire(mddDptrs[j], mxdIDptrs[j]);
          if (f == 0) continue;
          int u = getMddUnion(node[i], f);
          ddf->unlinkNode(f);
          ddf->unlinkNode(node[i]);
          node[i] = u;
        }
      }
      else {
        // Sparse node
        unsigned mxdISize = xdf->getSparseNodeSize(mxdDptrs[i]);
        const int* mxdIIptrs = 0;
        assert(xdf->getSparseNodeIndexes(mxdDptrs[i], mxdIIptrs));
        for (int k = 0; k < mxdISize; k++)
        {
          DCASSERT(0 != mxdIDptrs[k]);
          unsigned j = mxdIIptrs[k];
          if (j >= mddDptrs.size()) break;
          int f = reverseRecFire(mddDptrs[j], mxdIDptrs[k]);
          if (f == 0) continue;
          int u = getMddUnion(node[i], f);
          ddf->unlinkNode(f);
          ddf->unlinkNode(node[i]);
          node[i] = u;
        }
      }
    }
  }

  unsigned i = 0u;
  for ( ; i < node.size() && node[i] == 0; i++);
  if (i != node.size()) reverseSaturateHelper(nodeLevel, node);

  int n = ddf->createTempNode(nodeLevel, node);
  result = ddf->reduceNode(n);

  saveResult(owner, mdd, mxd, result);
  return result;
}


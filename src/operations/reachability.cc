
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

//#define IGNORE_TERMS 0
//#define IGNORE_INCOUNT 2
//#define DEBUG_DFS

//#define USE_GET_VECTOR_DOWN_POINTERS

inline expert_forest* getExpertForest(op_info* op, int index) {
  return op->p[index].getForest();
  // return smart_cast<expert_forest*>(op->f[index]);
}

inline const expert_forest* getExpertForest(const op_info* op, int index) {
  return op->p[index].readForest();
  // return smart_cast<const expert_forest*>(op->f[index]);
}

// ---------------------- MDD Traditional Reachability -------------------


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
    smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager());
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


// ---------------------- MDD Saturation-based Reachability -------------------


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

  for (unsigned i = 0; i < splits.size(); i++) {
    fprintf(stdout, "# of nodes in the next-state function: %1.6e\n",
        double(xdf->getNodeCount(splits[i])));
  }

  for (unsigned i = 0; i < splits.size(); i++) {
    fprintf(stdout, "# of edges in the next-state function: %1.6e\n",
        double(xdf->getEdgeCount(splits[i], false)));
  }

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
  ecm = smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager());
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


void mdd_reachability_dfs::clear()
{
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
}


// split is used to split a mxd for the saturation algorithm
void mdd_reachability_dfs::splitMxd(int mxd)
{
#if 1
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
#else
  xdf->linkNode(mxd);
  splits[xdf->getNodeLevel(mxd)] = mxd;
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

  int k = ddf->getNodeLevel(mdd);      // level
  int sz = ddf->getLevelSize(k);       // size

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d\n", mdd, k, sz);
#endif

  std::vector<int> node(sz, 0);
  std::vector<int> mddDptrs;
  ddf->getDownPtrs(mdd, mddDptrs);

  std::vector<int>::iterator nodeIter = node.begin();
  std::vector<int>::iterator mddIter = mddDptrs.begin();
  for ( ; mddIter != mddDptrs.end(); ++nodeIter, ++mddIter)
  {
    if (*mddIter != 0) *nodeIter = saturate(*mddIter);
  }
  
  // call saturateHelper for n
#ifdef DEBUG_DFS
  printf("Calling saturate: level %d\n", k);
#endif
  saturateHelper(k, node);

  // reduce and return
  int n = ddf->createTempNode(k, node);
  n = ddf->reduceNode(n);

#ifdef DEBUG_DFS
  ddf->showNodeGraph(stdout, n);
#endif

  return n;
}


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
  for (int i = 0; i < levelSize; i++) next[i] = 1;

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


int mdd_reachability_dfs::getMddUnion(int a, int b)
{
  return mddUnion->compute(mddUnionOp, a, b);
}

int mdd_reachability_dfs::getMxdIntersection(int a, int b)
{
  return mxdIntersection->compute(mxdIntersectionOp, a, b);
}

int mdd_reachability_dfs::getMxdDifference(int a, int b)
{
  return mxdDifference->compute(mxdDifferenceOp, a, b);
}


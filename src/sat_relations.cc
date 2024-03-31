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

#include "sat_relations.h"
#include "ops_builtin.h"

// ******************************************************************
// *                                                                *
// *                    pregen_relation  methods                    *
// *                                                                *
// ******************************************************************


MEDDLY::pregen_relation::pregen_relation(forest* mxd, unsigned nevents)
{
    if (!mxd) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    mxdF = mxd;
    if  (
            !mxdF->isForRelations() ||
            (mxdF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
        )
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    K = mxd->getMaxLevelIndex();

    num_events = nevents;
    if (num_events) {
        events = new dd_edge[num_events];
        next = new unsigned[num_events];
        for (unsigned e=0; e<num_events; e++) {
            events[e].attach(mxd);
        }
    } else {
        events = nullptr;
        next = nullptr;
    }
    last_event = 0;

    level_index = new unsigned[K+1];
    for (unsigned k=0; k<=K; k++) level_index[k] = 0;   // null pointer
}

MEDDLY::pregen_relation::pregen_relation(forest* mxd)
{
    if (!mxd) throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    mxdF = mxd;
    K = mxd->getMaxLevelIndex();

    events = new dd_edge[K+1];
    for (unsigned k=0; k<=K; k++) {
        events[k].attach(mxd);
    }

    next = nullptr;
    level_index = nullptr;

    num_events = 0;
    last_event = 0;
}


MEDDLY::pregen_relation::~pregen_relation()
{
    delete[] events;
    delete[] next;
    delete[] level_index;
}


void MEDDLY::pregen_relation::addToRelation(const dd_edge &r)
{
  MEDDLY_DCASSERT(mxdF);

  if (r.getForest() != mxdF)  throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

  int k = r.getLevel();
  if (0==k) return;
  if (k<0) k = -k;

  if (0==level_index) {
    // relation is "by levels"

    apply(UNION, events[k], r, events[k]);

  } else {
    // relation is "by events"

    if (isFinalized())            throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    if (last_event >= num_events) throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);

    events[last_event] = r;
    next[last_event] = level_index[k];
    level_index[k] = ++last_event;
  }
}


void MEDDLY::pregen_relation::splitMxd(splittingOption split)
{
  if (split == None) return;
  if (split == MonolithicSplit) {
    unionLevels();
    split = SplitOnly;
  }

  // For each level k, starting from the top level
  //    MXD(k) is the matrix diagram at level k.
  //    Calculate ID(k), the intersection of the diagonals of MXD(k).
  //    Subtract ID(k) from MXD(k).
  //    Add ID(k) to level(ID(k)).

#ifdef DEBUG_FINALIZE_SPLIT
  printf("Splitting events in finalize()\n");
  printf("events array: [");
  for (int i=0; i<=K; i++) {
    if (i) printf(", ");
    printf("%d", events[i]);
  }
  printf("]\n");
#endif

  // Initialize operations
  binary_operation* mxdUnion = UNION(mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdUnion);

  binary_operation* mxdIntersection = INTERSECTION(mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdIntersection);

  binary_operation* mxdDifference = DIFFERENCE(mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdDifference);

  dd_edge maxDiag(mxdF);

  for (int k = int(K); k > 1; k--) {
    if (0 == events[k].getNode()) continue;

    MEDDLY_DCASSERT(ABS(events[k].getLevel() <= k));

    // Initialize unpacked nodes
    unpacked_node* Mu = (isLevelAbove(k, events[k].getLevel()))
      ?   unpacked_node::newRedundant(mxdF, k, events[k].getNode(), FULL_ONLY)
      :   mxdF->newUnpacked(events[k].getNode(), FULL_ONLY)
    ;

    unpacked_node* Mp = unpacked_node::New(mxdF);

    // Read "rows"
    for (unsigned i = 0; i < Mu->getSize(); i++) {
      // Initialize column reader
      if (isLevelAbove(-k, mxdF->getNodeLevel(Mu->down(i)))) {
        Mp->initIdentity(mxdF, -k, i, Mu->down(i), FULL_ONLY);
      } else {
        mxdF->unpackNode(Mp, Mu->down(i), FULL_ONLY);
      }

      // Intersect along the diagonal
      if (0==i) {
        maxDiag.set( mxdF->linkNode(Mp->down(i)) );
      } else {
        dd_edge mpd(mxdF);
        mpd.set( mxdF->linkNode(Mp->down(i)) );
        mxdIntersection->computeTemp(maxDiag, mpd, maxDiag);
      }
    } // for i

    // Cleanup
    unpacked_node::Recycle(Mp);
    unpacked_node::Recycle(Mu);

    if (0 == maxDiag.getNode()) {
#ifdef DEBUG_FINALIZE_SPLIT
      printf("splitMxd: event %d, maxDiag %d\n", events[k], maxDiag);
#endif
      continue;
    }

    // Subtract maxDiag from events[k]
    // Do this only for SplitOnly. Other cases are handled later.

    if (split == SplitOnly) {
      mxdDifference->computeTemp(events[k], maxDiag, events[k]);
#ifdef DEBUG_FINALIZE_SPLIT
      printf("SplitOnly: event %d = event %d - maxDiag %d\n",
          events[k], tmp, maxDiag);
#endif
    }

    // Add maxDiag to events[level(maxDiag)]
    int maxDiagLevel = ABS(maxDiag.getLevel());

    mxdUnion->computeTemp(maxDiag, events[maxDiagLevel], events[maxDiagLevel]);

    // Subtract events[maxDiagLevel] from events[k].
    // Do this only for SplitSubtract. SplitSubtractAll is handled later.
    if (split == SplitSubtract) {
      mxdDifference->computeTemp(events[k], events[maxDiagLevel], events[k]);

#ifdef DEBUG_FINALIZE_SPLIT
      printf("SplitSubtract: event %d = event %d - event[maxDiagLevel] %d\n",
          events[k], tmp, maxDiag);
#endif
    }

  } // for k

  if (split == SplitSubtractAll) {
    // Subtract event[i] from all event[j], where j > i.
    for (int i = 1; i < K; i++) {
      if (0==events[i].getNode()) continue;

      for (int j = i + 1; j <= K; j++) {
        if (0==events[j].getNode()) continue;

        mxdDifference->computeTemp(events[j], events[i], events[j]);

#ifdef DEBUG_FINALIZE_SPLIT
        printf("SplitSubtractAll: event %d = event %d - event %d\n",
              events[j], tmp, events[i]);
#endif
      } // for j
    } // for i
  }

#ifdef DEBUG_FINALIZE_SPLIT
  printf("After splitting events in finalize()\n");
  printf("events array: [");
  for (int i=0; i<=K; i++) {
    if (i) printf(", ");
    printf("%d", events[i]);
  }
  printf("]\n");
#endif
}

// HERE!


void MEDDLY::pregen_relation::unionLevels()
{
  if (K < 1) return;

  binary_operation* mxdUnion = UNION(mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdUnion);

  dd_edge u(mxdF);
  for (unsigned k=1; k<=K; k++) {
    apply(UNION, u, events[k], u);
    events[k].set(0);
  }
  events[u.getLevel()] = u;
}


void MEDDLY::pregen_relation::finalize(splittingOption split)
{
  if (0==level_index) {
    // by levels
    switch (split) {
      case SplitOnly: printf("Split: SplitOnly\n"); break;
      case SplitSubtract: printf("Split: SplitSubtract\n"); break;
      case SplitSubtractAll: printf("Split: SplitSubtractAll\n"); break;
      case MonolithicSplit: printf("Split: MonolithicSplit\n"); break;
      default: printf("Split: None\n");
    }
    splitMxd(split);
    if (split != None && split != MonolithicSplit) {
#ifdef DEBUG_FINALIZE_SPLIT
      // Union the elements, and then re-run.
      // Result must be the same as before.
      dd_edge* old_events = new dd_edge[K+1];
      for(unsigned k = 0; k <= K; k++) {
        old_events[k] = events[k];
      }
      splitMxd(MonolithicSplit);
      binary_operation* mxdDifference = DIFFERENCE(mxdF, mxdF, mxdF);
      MEDDLY_DCASSERT(mxdDifference);
      for(unsigned k = 0; k <= K; k++) {
        if (old_events[k] != events[k]) {
          node_handle diff1 = mxdDifference->computeTemp(old_events[k], events[k]);
          node_handle diff2 = mxdDifference->computeTemp(events[k], old_events[k]);
          printf("error at level %d, n:k %d:%d, %d:%d\n",
              k,
              diff1, mxdF->getNodeLevel(diff1),
              diff2, mxdF->getNodeLevel(diff2)
              );
          mxdF->unlinkNode(diff1);
          mxdF->unlinkNode(diff2);
        }
      }

      delete [] old_events;
#endif
    }
    return;
  }

  //
  // Still here?  Must be by events.
  //

#ifdef DEBUG_FINALIZE
  printf("Finalizing pregen relation\n");
  printf("%u events total\n", last_event);
  printf("events array: [");
  for (unsigned i=0; i<last_event; i++) {
    if (i) printf(", ");
    printf("%ld", long(events[i].getNode()));
  }
  printf("]\n");
  printf("next array: [");
  for (unsigned i=0; i<last_event; i++) {
    if (i) printf(", ");
    printf("%d", next[i]);
  }
  printf("]\n");
  printf("level_index array: [%d", level_index[1]);
  for (int i=2; i<=K; i++) {
    printf(", %d", level_index[i]);
  }
  printf("]\n");
#endif

  //
  // Convert from array of linked lists to contiguous array.
  //
  dd_edge* new_events = new dd_edge[last_event];
  unsigned P = 0;
  for (unsigned k=K; k; k--) {
    unsigned L = level_index[k];
    level_index[k] = P;
    while (L) {
      // L+1 is index of an element at level k
      L--;
      new_events[P++] = events[L];
      L = next[L];
    } // while L
  }
  level_index[0] = P;
  delete[] events;
  events = new_events;

  // done with next pointers
  delete[] next;
  next = 0;

#ifdef DEBUG_FINALIZE
  printf("\nAfter finalization\n");
  printf("events array: [");
  for (unsigned i=0; i<last_event; i++) {
    if (i) printf(", ");
    printf("%ld", long(events[i].getNode()));
  }
  printf("]\n");
  printf("level_index array: [%d", level_index[1]);
  for (int i=2; i<=K; i++) {
    printf(", %d", level_index[i]);
  }
  printf("]\n");
#endif
}


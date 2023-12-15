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

#include "../defines.h"
#include "sat_pregen.h"
#include <typeinfo> // for "bad_cast" exception

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../oper_special.h"
#include "../opname_satur.h"
#include "../ops_builtin.h"

#define DEBUG_FINALIZE
// #define DEBUG_FINALIZE_SPLIT
// #define DEBUG_EVENT_MASK

namespace MEDDLY {
  class saturation_by_events_opname;
  class saturation_by_events_op;

  class common_dfs_by_events_mt;
  class forwd_dfs_by_events_mt;
  class bckwd_dfs_by_events_mt;

  class fb_saturation_opname;

} // Namespace MEDDLY



// ******************************************************************
// *                                                                *
// *                    satpregen_opname methods                    *
// *                                                                *
// ******************************************************************


MEDDLY::satpregen_opname::satpregen_opname(const char* n)
 : specialized_opname(n)
{
}

MEDDLY::satpregen_opname::~satpregen_opname()
{
}

void
MEDDLY::satpregen_opname::pregen_relation
::setForests(forest* inf, forest* mxd, forest* outf)
{
  insetF = inf;
  outsetF = outf;
  mxdF = smart_cast <MEDDLY::expert_forest*>(mxd);
  if (0==insetF || 0==outsetF || 0==mxdF) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);

  // Check for same domain
  if (
    (insetF->getDomain() != mxdF->getDomain()) ||
    (outsetF->getDomain() != mxdF->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  // for now, anyway, inset and outset must be same forest
  if (insetF != outsetF)
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

  // Check forest types
  if (
    insetF->isForRelations()    ||
    !mxdF->isForRelations()     ||
    outsetF->isForRelations()   ||
    (insetF->getRangeType() != mxdF->getRangeType())        ||
    (outsetF->getRangeType() != mxdF->getRangeType())       ||
    (insetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)   ||
    (outsetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)  ||
    (mxdF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)     ||
    (outsetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  // Forests are good; set number of variables
  K = mxdF->getMaxLevelIndex();
}

MEDDLY::satpregen_opname::pregen_relation
::pregen_relation(forest* inf, forest* mxd, forest* outf, unsigned nevents)
{
  setForests(inf, mxd, outf);

  num_events = nevents;
  if (num_events) {
    events = new dd_edge[num_events];
    next = new unsigned[num_events];
    for (unsigned e=0; e<num_events; e++) {
      events[e].attach(mxdF);
    }
  } else {
    events = 0;
    next = 0;
  }
  last_event = 0;

  level_index = new unsigned[K+1];
  for (unsigned k=0; k<=K; k++) level_index[k] = 0;   // null pointer
}

MEDDLY::satpregen_opname::pregen_relation
::pregen_relation(forest* inf, forest* mxd, forest* outf)
{
  setForests(inf, mxd, outf);

  events = new dd_edge[K+1];
  for (unsigned k=0; k<=K; k++) {
      events[k].attach(mxdF);
  }

  next = 0;
  level_index = 0;

  num_events = 0;
  last_event = 0;
}


MEDDLY::satpregen_opname::pregen_relation
::~pregen_relation()
{
  delete[] events;
  delete[] next;
  delete[] level_index;
}

void
MEDDLY::satpregen_opname::pregen_relation
::addToRelation(const dd_edge &r)
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


void
MEDDLY::satpregen_opname::pregen_relation
::splitMxd(splittingOption split)
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
  binary_operation* mxdUnion = getOperation(UNION, mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdUnion);

  binary_operation* mxdIntersection =
    getOperation(INTERSECTION, mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdIntersection);

  binary_operation* mxdDifference = getOperation(DIFFERENCE, mxdF, mxdF, mxdF);
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


void
MEDDLY::satpregen_opname::pregen_relation
::unionLevels()
{
  if (K < 1) return;

  binary_operation* mxdUnion = getOperation(UNION, mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdUnion);

  dd_edge u(mxdF);
  for (unsigned k=1; k<=K; k++) {
    apply(UNION, u, events[k], u);
    events[k].set(0);
  }
  events[u.getLevel()] = u;
}


void
MEDDLY::satpregen_opname::pregen_relation
::finalize(splittingOption split)
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
      binary_operation* mxdDifference =
        getOperation(DIFFERENCE, mxdF, mxdF, mxdF);
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

// ******************************************************************
// *                                                                *
// *               saturation_by_events_opname  class               *
// *                                                                *
// ******************************************************************

/** Simple class to keep compute table happy.
*/
class MEDDLY::saturation_by_events_opname : public unary_opname {
  static saturation_by_events_opname* instance;
  public:
    saturation_by_events_opname();

    static saturation_by_events_opname* getInstance();

};

MEDDLY::saturation_by_events_opname* MEDDLY::saturation_by_events_opname::instance = 0;

MEDDLY::saturation_by_events_opname::saturation_by_events_opname()
 : unary_opname("Saturate_by_events")
{
}

MEDDLY::saturation_by_events_opname* MEDDLY::saturation_by_events_opname::getInstance()
{
  if (0==instance) instance = new saturation_by_events_opname;
  return instance;
}

// ******************************************************************
// *                                                                *
// *                      saturation_by_events_op  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_by_events_op : public unary_operation {
    common_dfs_by_events_mt* parent;
  public:
    saturation_by_events_op(common_dfs_by_events_mt* p,
      expert_forest* argF, expert_forest* resF);
    virtual ~saturation_by_events_op();

    void saturate(const dd_edge& in, dd_edge& out);
    node_handle saturate(node_handle mdd, int level);

  protected:
    inline ct_entry_key*
    findSaturateResult(node_handle a, int level, node_handle& b) {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      if (argF->isFullyReduced()) CTsrch->writeI(level);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      b = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveSaturateResult(ct_entry_key* Key,
      node_handle a, node_handle b)
    {
      CTresult[0].reset();
      CTresult[0].writeN(b);
      CT0->addEntry(Key, CTresult[0]);
      return b;
    }
};


// ******************************************************************
// *                                                                *
// *            common_dfs_by_events_mt  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::common_dfs_by_events_mt : public specialized_operation {
  public:
    common_dfs_by_events_mt(satpregen_opname* opcode,
      satpregen_opname::pregen_relation* rel);
    virtual ~common_dfs_by_events_mt();

    virtual void compute(const dd_edge& a, dd_edge &c);
    virtual void saturateHelper(unpacked_node& mdd) = 0;

  protected:
    inline ct_entry_key*
    findResult(node_handle a, node_handle b, node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CTsrch->writeN(b);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(ct_entry_key* Key,
      node_handle a, node_handle b, node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
      return c;
    }

  protected:
    binary_operation* mddUnion;
    binary_operation* mxdIntersection;
    binary_operation* mxdDifference;

    satpregen_opname::pregen_relation* rel;

    expert_forest* arg1F;
    expert_forest* arg2F;
    expert_forest* resF;

  protected:
    class indexq {
        static const int NULPTR = -1;
        static const int NOTINQ = -2;
        int* data;
        unsigned size;
        int head;
        int tail;
      public:
        // used by parent for recycling
        indexq* next;
      public:
        indexq();
        ~indexq();
        void resize(unsigned sz);
        inline bool isEmpty() const {
          return NULPTR == head;
        }
        inline void add(unsigned i) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i, size);
          if (NOTINQ != data[i]) return;
          if (NULPTR == head) {
            // empty list
            head = int(i);
          } else {
            // not empty list
              MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, (unsigned) tail, size);
            data[tail] = int(i);
          }
          tail = int(i);
          data[i] = NULPTR;
        }
        inline unsigned remove() {
          MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, (unsigned)head, size);
          unsigned ans = unsigned(head);
          head = data[head];
          data[ans] = NOTINQ;
          MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, ans, size);
          return ans;
        }
    };

  protected:
    class charbuf {
      public:
        char* data;
        unsigned size;
        charbuf* next;
      public:
        charbuf();
        ~charbuf();
        void resize(unsigned sz);
    };

  private:
    indexq* freeqs;
    charbuf* freebufs;

  protected:
    inline indexq* useIndexQueue(unsigned sz) {
      indexq* ans;
      if (freeqs) {
        ans = freeqs;
        freeqs = freeqs->next;
      } else {
        ans = new indexq();
      }
      MEDDLY_DCASSERT(ans);
      ans->resize(sz);
      ans->next = 0;
      return ans;
    }
    inline void recycle(indexq* a) {
      MEDDLY_DCASSERT(a);
      MEDDLY_DCASSERT(a->isEmpty());
      a->next = freeqs;
      freeqs = a;
    }

    inline charbuf* useCharBuf(unsigned sz) {
      charbuf* ans;
      if (freebufs) {
        ans = freebufs;
        freebufs = freebufs->next;
      } else {
        ans = new charbuf();
      }
      MEDDLY_DCASSERT(ans);
      ans->resize(sz);
      ans->next = 0;
      return ans;
    }
    inline void recycle(charbuf* a) {
      MEDDLY_DCASSERT(a);
      a->next = freebufs;
      freebufs = a;
    }

    inline virtual bool checkForestCompatibility() const {
      return true;
    }
};

// ******************************************************************
// *                                                                *
// *           saturation_by_events_op  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_by_events_op
::saturation_by_events_op(common_dfs_by_events_mt* p,
  expert_forest* argF, expert_forest* resF)
  : unary_operation(saturation_by_events_opname::getInstance(), 1, argF, resF)
{
  parent = p;

  const char* name = saturation_by_events_opname::getInstance()->getName();
  ct_entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new ct_entry_type(name, "NI:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(3, resF);
  } else {
    et = new ct_entry_type(name, "N:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(2, resF);
  }
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::saturation_by_events_op::~saturation_by_events_op()
{
  removeAllComputeTableEntries();
}

void MEDDLY::saturation_by_events_op::saturate(const dd_edge& in, dd_edge& out)
{
  out.set( saturate(in.getNode(), argF->getMaxLevelIndex()) );
}

MEDDLY::node_handle
MEDDLY::saturation_by_events_op::saturate(node_handle mdd, int k)
{
#ifdef DEBUG_DFS
  printf("mdd: %d, k: %d\n", mdd, k);
#endif

  // terminal condition for recursion
  if (argF->isTerminalNode(mdd)) return mdd;

  // search compute table
  node_handle n = 0;
  ct_entry_key* Key = findSaturateResult(mdd, k, n);
  if (0==Key) return n;

  unsigned sz = unsigned(argF->getLevelSize(k));    // size
  int mdd_level = argF->getNodeLevel(mdd);          // mdd level

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d, mdd_level: %d\n",
      mdd, k, sz, mdd_level);
#endif

  unpacked_node* nb = unpacked_node::newFull(resF, k, sz);
  // Initialize mdd reader
  unpacked_node *mddDptrs = unpacked_node::New(argF);
  if (mdd_level < k) {
    mddDptrs->initRedundant(argF, k, mdd, FULL_ONLY);
  } else {
    argF->unpackNode(mddDptrs, mdd, FULL_ONLY);
  }

  // Do computation
  for (unsigned i=0; i<sz; i++) {
    nb->d_ref(i) = mddDptrs->down(i) ? saturate(mddDptrs->down(i), k-1) : 0;
  }

  // Cleanup
  unpacked_node::Recycle(mddDptrs);

  parent->saturateHelper(*nb);
  n = resF->createReducedNode(-1, nb);

  // save in compute table
  saveSaturateResult(Key, mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  return n;
}


// ******************************************************************
// *                                                                *
// *           common_dfs_by_events_mt  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::common_dfs_by_events_mt::common_dfs_by_events_mt(
  satpregen_opname* opcode,
  satpregen_opname::pregen_relation* relation)
: specialized_operation(opcode, 1)
{
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;
  freeqs = 0;
  freebufs = 0;
  rel = relation;
  arg1F = static_cast<expert_forest*>(rel->getInForest());
  arg2F = static_cast<expert_forest*>(rel->getRelForest());
  resF = static_cast<expert_forest*>(rel->getOutForest());

  registerInForest(arg1F);
  registerInForest(arg2F);
  registerInForest(resF);

  ct_entry_type* et = new ct_entry_type(opcode->getName(), "NN:N");
  et->setForestForSlot(0, arg1F);
  et->setForestForSlot(1, arg2F);
  et->setForestForSlot(3, resF);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::common_dfs_by_events_mt::~common_dfs_by_events_mt()
{
  if (rel->autoDestroy()) delete rel;
  unregisterInForest(arg1F);
  unregisterInForest(arg2F);
  unregisterInForest(resF);
}

void MEDDLY::common_dfs_by_events_mt
::compute(const dd_edge &a, dd_edge &c)
{
  // Initialize operations
  mddUnion = getOperation(UNION, resF, resF, resF);
  MEDDLY_DCASSERT(mddUnion);

  mxdIntersection = getOperation(INTERSECTION, arg2F, arg2F, arg2F);
  MEDDLY_DCASSERT(mxdIntersection);

  mxdDifference = getOperation(DIFFERENCE, arg2F, arg2F, arg2F);
  MEDDLY_DCASSERT(mxdDifference);

#ifdef DEBUG_INITIAL
  printf("Calling saturate for states:\n");
  a.showGraph(stdout);
#endif
#ifdef DEBUG_NSF
  printf("Calling saturate for NSF:\n");
  // b.showGraph(stdout);
#endif

  // Execute saturation operation
  if (!rel->isFinalized()) {
    printf("Transition relation has not been finalized.\n");
    printf("Finalizing using default options... ");
    rel->finalize();
    printf("done.\n");
  }
  saturation_by_events_op* so = new saturation_by_events_op(this, arg1F, resF);
  so->saturate(a, c);

  // Cleanup
  while (freeqs) {
    indexq* t = freeqs;
    freeqs = t->next;
    delete t;
  }
  while (freebufs) {
    charbuf* t = freebufs;
    freebufs = t->next;
    delete t;
  }
  delete so;
}

// ******************************************************************
// *       common_dfs_by_events_mt::indexq  methods                 *
// ******************************************************************

MEDDLY::common_dfs_by_events_mt::indexq::indexq()
{
  data = 0;
  size = 0;
  head = NULPTR;
}

MEDDLY::common_dfs_by_events_mt::indexq::~indexq()
{
  free(data);
}

void MEDDLY::common_dfs_by_events_mt::indexq::resize(unsigned sz)
{
  if (sz <= size) return;
  data = (int*) realloc(data, sz * sizeof(int));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);

  for (; size < sz; size++) data[size] = NOTINQ;
}

// ******************************************************************
// *       common_dfs_by_events_mt::charbuf methods                 *
// ******************************************************************

MEDDLY::common_dfs_by_events_mt::charbuf::charbuf()
{
  data = 0;
  size = 0;
}

MEDDLY::common_dfs_by_events_mt::charbuf::~charbuf()
{
  free(data);
}

void MEDDLY::common_dfs_by_events_mt::charbuf::resize(unsigned sz)
{
  if (sz <= size) return;
  data = (char*) realloc(data, sz * sizeof(char));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *             forwd_dfs_by_events_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_dfs_by_events_mt : public common_dfs_by_events_mt {
  public:
    forwd_dfs_by_events_mt(satpregen_opname* opcode,
    satpregen_opname::pregen_relation* rel);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::forwd_dfs_by_events_mt::forwd_dfs_by_events_mt(
  satpregen_opname* opcode,
  satpregen_opname::pregen_relation* rel)
  : common_dfs_by_events_mt(opcode, rel)
{
}


void MEDDLY::forwd_dfs_by_events_mt::saturateHelper(unpacked_node& nb)
{
  unsigned nEventsAtThisLevel = rel->lengthForLevel(nb.getLevel());
  if (0 == nEventsAtThisLevel) return;

  // Initialize mxd readers, note we might skip the unprimed level
  dd_edge* events = rel->arrayForLevel(nb.getLevel());
  unpacked_node** Ru = new unpacked_node*[nEventsAtThisLevel];
  for (unsigned ei = 0; ei < nEventsAtThisLevel; ei++) {
    Ru[ei] = unpacked_node::New(arg2F);
    int eventLevel = events[ei].getLevel();
    MEDDLY_DCASSERT(ABS(eventLevel) == nb.getLevel());
    if (eventLevel<0) {
      Ru[ei]->initRedundant(arg2F, nb.getLevel(), events[ei].getNode(), FULL_ONLY);
    } else {
      arg2F->unpackNode(Ru[ei], events[ei].getNode(), FULL_ONLY);
    }
  }
  unpacked_node* Rp = unpacked_node::New(arg2F);

  dd_edge nbdj(resF), newst(resF);

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (unsigned i = 0; i < nb.getSize(); i++) {
    if (nb.d(i)) queue->add(i);
  }

  // explore indexes
  while (!queue->isEmpty()) {
    unsigned i = queue->remove();

    MEDDLY_DCASSERT(nb.d(i));

    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
      if (0 == Ru[ei]->down(i)) continue;  // row i of the event ei is empty

      // grab column (TBD: build these ahead of time?)
      int dlevel = arg2F->getNodeLevel(Ru[ei]->down(i));

      if (dlevel == -nb.getLevel()) {
        arg2F->unpackNode(Rp, Ru[ei]->down(i), SPARSE_ONLY);
      } else {
        Rp->initIdentity(arg2F, -nb.getLevel(), i, Ru[ei]->down(i), SPARSE_ONLY);
      }

      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        unsigned j = Rp->index(jz);
        if (-1==nb.d(j)) continue;  // nothing can be added to this set

        node_handle rec = recFire(nb.d(i), Rp->down(jz));

        if (rec == 0) continue;
        if (rec == nb.d(j)) {
          resF->unlinkNode(rec);
          continue;
        }

        bool updated = true;

        if (0 == nb.d(j)) {
          nb.d_ref(j) = rec;
        }
        else if (rec == -1) {
          resF->unlinkNode(nb.d(j));
          nb.d_ref(j) = -1;
        }
        else {
          nbdj.set(nb.d(j));  // clobber
          newst.set(rec);     // clobber
          mddUnion->computeTemp(nbdj, newst, nbdj);
          updated = (nbdj.getNode() != nb.d(j));
          nb.setFull(j, nbdj);
        }

        if (updated) queue->add(j);
      } // for j
    } // for all events, ei
  } // while there are indexes to explore

  // cleanup
  unpacked_node::Recycle(Rp);
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) unpacked_node::Recycle(Ru[ei]);
  delete[] Ru;
  recycle(queue);
}


// Same as post-image, except we saturate before reducing.
MEDDLY::node_handle MEDDLY::forwd_dfs_by_events_mt::recFire(
  node_handle mdd, node_handle mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (arg2F->isTerminalNode(mxd)) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->handleForValue(1);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  node_handle result = 0;
  ct_entry_key* Key = findResult(mdd, mxd, result);
  if (0==Key) return result;

#ifdef TRACE_RECFIRE
  printf("computing recFire(%d, %d)\n", mdd, mxd);
  printf("  node %3d ", mdd);
  arg1F->showNode(stdout, mdd, 1);
  printf("\n  node %3d ", mxd);
  arg2F->showNode(stdout, mxd, 1);
  printf("\n");
#endif

  // check if mxd and mdd are at the same level
  const int mddLevel = arg1F->getNodeLevel(mdd);
  const int mxdLevel = arg2F->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);

  dd_edge nbdj(resF), newst(resF);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(arg1F);
  if (mddLevel < rLevel) {
    A->initRedundant(arg1F, rLevel, mdd, FULL_ONLY);
  } else {
    arg1F->unpackNode(A, mdd, FULL_ONLY);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (unsigned i=0; i<rSize; i++) {
      nb->d_ref(i) = recFire(A->down(i), mxd);
    }

  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New(arg2F);
    unpacked_node *Rp = unpacked_node::New(arg2F);
    if (mxdLevel < 0) {
      Ru->initRedundant(arg2F, rLevel, mxd, SPARSE_ONLY);
    } else {
      arg2F->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      unsigned i = Ru->index(iz);
      if (0==A->down(i))   continue;
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(arg2F, rLevel, i, Ru->down(iz), SPARSE_ONLY);
      } else {
        arg2F->unpackNode(Rp, Ru->down(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        unsigned j = Rp->index(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = recFire(A->down(i), Rp->down(jz));
        if (0==newstates) continue;
        if (0==nb->down(j)) {
          nb->d_ref(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        nbdj.set(nb->down(j));
        newst.set(newstates);
        mddUnion->computeTemp(nbdj, newst, nbdj);
        nb->setFull(j, nbdj);
      } // for j

    } // for i

    unpacked_node::Recycle(Rp);
    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);

  saturateHelper(*nb);
  result = resF->createReducedNode(-1, nb);
#ifdef TRACE_ALL_OPS
  printf("computed recfire(%d, %d) = %d\n", mdd, mxd, result);
#endif
#ifdef TRACE_RECFIRE
  printf("computed recfire(%d, %d) = %d\n", mdd, mxd, result);
  printf("  node %3d ", result);
  resF->showNode(stdout, result, 1);
  printf("\n");
#endif
  return saveResult(Key, mdd, mxd, result);
}




// ******************************************************************
// *                                                                *
// *             bckwd_dfs_by_events_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_dfs_by_events_mt : public common_dfs_by_events_mt {
  public:
    bckwd_dfs_by_events_mt(satpregen_opname* opcode,
    satpregen_opname::pregen_relation* rel);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::bckwd_dfs_by_events_mt::bckwd_dfs_by_events_mt(
  satpregen_opname* opcode,
  satpregen_opname::pregen_relation* rel)
  : common_dfs_by_events_mt(opcode, rel)
{
}

void MEDDLY::bckwd_dfs_by_events_mt::saturateHelper(unpacked_node& nb)
{
  unsigned nEventsAtThisLevel = rel->lengthForLevel(nb.getLevel());
  if (0 == nEventsAtThisLevel) return;

  // Initialize mxd readers, note we might skip the unprimed level
  dd_edge* events = rel->arrayForLevel(nb.getLevel());
  unpacked_node** Ru = new unpacked_node*[nEventsAtThisLevel];
  for (unsigned ei = 0; ei < nEventsAtThisLevel; ei++) {
    Ru[ei] = unpacked_node::New(arg2F);
    int eventLevel = events[ei].getLevel();
    MEDDLY_DCASSERT(ABS(eventLevel) == nb.getLevel());
    if (eventLevel<0) {
      Ru[ei]->initRedundant(arg2F, nb.getLevel(), events[ei].getNode(), FULL_ONLY);
    } else {
      arg2F->unpackNode(Ru[ei], events[ei].getNode(), FULL_ONLY);
    }
  }
  unpacked_node* Rp = unpacked_node::New(arg2F);

  dd_edge nbdi(resF), newst(resF);

  // indexes to explore
  charbuf* expl = useCharBuf(nb.getSize());
  for (unsigned i = 0; i < nb.getSize(); i++) expl->data[i] = 2;
  bool repeat = true;

  // explore
  while (repeat) {
    // "advance" the explore list
    for (unsigned i=0; i<nb.getSize(); i++) if (expl->data[i]) expl->data[i]--;
    repeat = false;

    // explore all events
    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
      // explore all rows
      for (unsigned iz=0; iz<Ru[ei]->getSize(); iz++) {
        unsigned i = Ru[ei]->index(iz);
        // grab column (TBD: build these ahead of time?)
        int dlevel = arg2F->getNodeLevel(Ru[ei]->down(iz));

        if (dlevel == -nb.getLevel()) {
          arg2F->unpackNode(Rp, Ru[ei]->down(iz), SPARSE_ONLY);
        } else {
          Rp->initIdentity(arg2F, -nb.getLevel(), i, Ru[ei]->down(iz), SPARSE_ONLY);
        }

        for (unsigned jz=0; jz<Rp->getSize(); jz++) {
          unsigned j = Rp->index(jz);
          if (0==expl->data[j]) continue;
          if (0==nb.d(j))       continue;
          // We have an i->j edge to explore
          node_handle rec = recFire(nb.d(j), Rp->down(jz));

          if (0==rec) continue;
          if (rec == nb.d(i)) {
            resF->unlinkNode(rec);
            continue;
          }

          bool updated = true;

          if (0 == nb.d(i)) {
            nb.d_ref(i) = rec;
          }
          else if (-1 == rec) {
            resF->unlinkNode(nb.d(i));
            nb.d_ref(i) = -1;
          }
          else {
            nbdi.set(nb.d(i));
            newst.set(rec);
            mddUnion->computeTemp(nbdi, newst, nbdi);
            updated = (nbdi.getNode() != nb.d(i));
            nb.setFull(i, nbdi);
          }
          if (updated) {
            expl->data[i] = 2;
            repeat = true;
          }
        } // for j
      } // for i
    } // for each event
  } // while repeat

  // cleanup
  unpacked_node::Recycle(Rp);
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) unpacked_node::Recycle(Ru[ei]);
  delete[] Ru;
  recycle(expl);
}

MEDDLY::node_handle MEDDLY::bckwd_dfs_by_events_mt::recFire(node_handle mdd,
  node_handle mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (arg2F->isTerminalNode(mxd)) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->handleForValue(1);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  node_handle result = 0;
  ct_entry_key* Key = findResult(mdd, mxd, result);
  if (0==Key) return result;

  // check if mxd and mdd are at the same level
  const int mddLevel = arg1F->getNodeLevel(mdd);
  const int mxdLevel = arg2F->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);

  dd_edge nbdi(resF), newst(resF);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(arg1F);
  if (mddLevel < rLevel) {
    A->initRedundant(arg1F, rLevel, mdd, FULL_ONLY);
  } else {
    arg1F->unpackNode(A, mdd, FULL_ONLY);
  }


  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (unsigned i=0; i<rSize; i++) {
      nb->d_ref(i) = recFire(A->down(i), mxd);
    }
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New(arg2F);
    unpacked_node *Rp = unpacked_node::New(arg2F);
    if (mxdLevel < 0) {
      Ru->initRedundant(arg2F, rLevel, mxd, SPARSE_ONLY);
    } else {
      arg2F->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      unsigned i = Ru->index(iz);
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(arg2F, rLevel, i, Ru->down(iz), SPARSE_ONLY);
      } else {
        arg2F->unpackNode(Rp, Ru->down(iz), SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        unsigned j = Rp->index(jz);
        if (0==A->down(j))   continue;
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = recFire(A->down(j), Rp->down(jz));
        if (0==newstates) continue;
        if (0==nb->down(i)) {
          nb->d_ref(i) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        nbdi.set(nb->down(i));
        newst.set(newstates);
        mddUnion->computeTemp(nbdi, newst, nbdi);
        nb->setFull(i, nbdi);
      } // for j

    } // for i

    unpacked_node::Recycle(Rp);
    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);

  saturateHelper(*nb);
  result = resF->createReducedNode(-1, nb);
#ifdef TRACE_ALL_OPS
  printf("computed recFire(%d, %d) = %d\n", mdd, mxd, result);
#endif
  return saveResult(Key, mdd, mxd, result);
}


// ******************************************************************
// *                                                                *
// *                   fb_saturation_opname class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::fb_saturation_opname : public satpregen_opname {
    bool forward;
  public:
    fb_saturation_opname(bool fwd);
    virtual specialized_operation* buildOperation(arguments* a);
};

MEDDLY::fb_saturation_opname::fb_saturation_opname(bool fwd)
 : satpregen_opname(fwd ? "SaturationFwd" : "SaturationBack")
{
  forward = fwd;
}

MEDDLY::specialized_operation*
MEDDLY::fb_saturation_opname::buildOperation(arguments* a)
{
  pregen_relation* rel = dynamic_cast<pregen_relation*>(a);
  if (0==rel) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);

  //
  // No sanity checks needed here; we did them already when constructing a.
  //

  MEDDLY::specialized_operation* op = 0;
  if (forward)
    op = new forwd_dfs_by_events_mt(this, rel);
  else
    op = new bckwd_dfs_by_events_mt(this, rel);

  // Do we need to delete rel here?
  // No, if needed, do this in the destructor for op.

  return op;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::satpregen_opname* MEDDLY::initSaturationForward()
{
  return new fb_saturation_opname(true);
}

MEDDLY::satpregen_opname* MEDDLY::initSaturationBackward()
{
  return new fb_saturation_opname(false);
}



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

namespace MEDDLY {
  class saturation_by_events_opname;
  class saturation_by_events_op;

  class common_dfs_by_events_mt;
  class forwd_dfs_by_events_mt;
  class bckwd_dfs_by_events_mt;

  class fb_saturation_opname;
};


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

    static const saturation_by_events_opname* getInstance();
 
};

MEDDLY::saturation_by_events_opname* MEDDLY::saturation_by_events_opname::instance = 0;

MEDDLY::saturation_by_events_opname::saturation_by_events_opname()
 : unary_opname("Saturate_by_events")
{
}

const MEDDLY::saturation_by_events_opname* MEDDLY::saturation_by_events_opname::getInstance()
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

    node_handle saturate(node_handle mdd);
    node_handle saturate(node_handle mdd, int level);

  protected:
    inline compute_table::entry_key* 
    findSaturateResult(node_handle a, int level, node_handle& b) {
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      if (argF->isFullyReduced()) CTsrch->writeI(level);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      b = resF->linkNode(CTresult[0].readN()); 
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveSaturateResult(compute_table::entry_key* Key,
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
    common_dfs_by_events_mt(const satpregen_opname* opcode,
      satpregen_opname::pregen_relation* rel);
    virtual ~common_dfs_by_events_mt();

    virtual void compute(const dd_edge& a, dd_edge &c);
    virtual void saturateHelper(unpacked_node& mdd) = 0;

  protected:
    inline compute_table::entry_key* 
    findResult(node_handle a, node_handle b, node_handle &c) 
    {
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CTsrch->writeN(b);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(compute_table::entry_key* Key,
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
        int size;
        int head;
        int tail;
      public:
        // used by parent for recycling
        indexq* next;
      public:
        indexq();
        ~indexq();
        void resize(int sz);
        inline bool isEmpty() const {
          return NULPTR == head;
        }
        inline void add(int i) {
          MEDDLY_CHECK_RANGE(0, i, size);
          if (NOTINQ != data[i]) return;
          if (NULPTR == head) {
            // empty list
            head = i;
          } else {
            // not empty list
            MEDDLY_CHECK_RANGE(0, tail, size);
            data[tail] = i;
          }
          tail = i;
          data[i] = NULPTR;
        }
        inline int remove() {
          MEDDLY_CHECK_RANGE(0, head, size);
          int ans = head;
          head = data[head];
          data[ans] = NOTINQ;
          MEDDLY_CHECK_RANGE(0, ans, size);
          return ans;
        }
    };

  protected:
    class charbuf {
      public:
        char* data;
        int size;
        charbuf* next;
      public:
        charbuf();
        ~charbuf();
        void resize(int sz);
    };

  private:
    indexq* freeqs;
    charbuf* freebufs;

  protected:
    inline indexq* useIndexQueue(int sz) {
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

    inline charbuf* useCharBuf(int sz) {
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
  compute_table::entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new compute_table::entry_type(name, "NI:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(3, resF);
  } else {
    et = new compute_table::entry_type(name, "N:N");
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

MEDDLY::node_handle MEDDLY::saturation_by_events_op::saturate(node_handle mdd)
{
  return saturate(mdd, argF->getNumVariables());
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
  compute_table::entry_key* Key = findSaturateResult(mdd, k, n);
  if (0==Key) return n;

  int sz = argF->getLevelSize(k);               // size
  int mdd_level = argF->getNodeLevel(mdd);      // mdd level

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d, mdd_level: %d\n",
      mdd, k, sz, mdd_level);
#endif

  unpacked_node* nb = unpacked_node::newFull(resF, k, sz);
  // Initialize mdd reader
  unpacked_node *mddDptrs = unpacked_node::useUnpackedNode();
  if (mdd_level < k) {
    mddDptrs->initRedundant(argF, k, mdd, true);
  } else {
    mddDptrs->initFromNode(argF, mdd, true);
  }

  // Do computation
  for (int i=0; i<sz; i++) {
    nb->d_ref(i) = mddDptrs->d(i) ? saturate(mddDptrs->d(i), k-1) : 0;
  }

  // Cleanup
  unpacked_node::recycle(mddDptrs);

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
  const satpregen_opname* opcode,
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

  compute_table::entry_type* et = new compute_table::entry_type(opcode->getName(), "NN:N");
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
  a.show(stdout, 2);
#endif
#ifdef DEBUG_NSF
  printf("Calling saturate for NSF:\n");
  // b.show(stdout, 2);
#endif

  // Execute saturation operation
  if (!rel->isFinalized()) {
    printf("Transition relation has not been finalized.\n");
    printf("Finalizing using default options... ");
    rel->finalize();
    printf("done.\n");
  }
  saturation_by_events_op* so = new saturation_by_events_op(this, arg1F, resF);
  node_handle cnode = so->saturate(a.getNode());
  c.set(cnode);

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

void MEDDLY::common_dfs_by_events_mt::indexq::resize(int sz)
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

void MEDDLY::common_dfs_by_events_mt::charbuf::resize(int sz)
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
    forwd_dfs_by_events_mt(const satpregen_opname* opcode,
    satpregen_opname::pregen_relation* rel);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::forwd_dfs_by_events_mt::forwd_dfs_by_events_mt(
  const satpregen_opname* opcode,
  satpregen_opname::pregen_relation* rel)
  : common_dfs_by_events_mt(opcode, rel)
{
}


void MEDDLY::forwd_dfs_by_events_mt::saturateHelper(unpacked_node& nb)
{
  int nEventsAtThisLevel = rel->lengthForLevel(nb.getLevel());
  if (0 == nEventsAtThisLevel) return;

  // Initialize mxd readers, note we might skip the unprimed level
  node_handle* events = rel->arrayForLevel(nb.getLevel());
  unpacked_node** Ru = new unpacked_node*[nEventsAtThisLevel];
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    Ru[ei] = unpacked_node::useUnpackedNode();
    int eventLevel = arg2F->getNodeLevel(events[ei]);
    MEDDLY_DCASSERT(ABS(eventLevel) == nb.getLevel());
    if (eventLevel<0) {
      Ru[ei]->initRedundant(arg2F, nb.getLevel(), events[ei], true);
    } else {
      Ru[ei]->initFromNode(arg2F, events[ei], true);
    }
  }
  unpacked_node* Rp = unpacked_node::useUnpackedNode();

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.d(i)) queue->add(i);
  }

  // explore indexes
  while (!queue->isEmpty()) {
    int i = queue->remove();

    MEDDLY_DCASSERT(nb.d(i));

    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
      if (0 == Ru[ei]->d(i)) continue;  // row i of the event ei is empty

      // grab column (TBD: build these ahead of time?)
      int dlevel = arg2F->getNodeLevel(Ru[ei]->d(i));

      if (dlevel == -nb.getLevel()) {
        Rp->initFromNode(arg2F, Ru[ei]->d(i), false);
      } else {
        Rp->initIdentity(arg2F, -nb.getLevel(), i, Ru[ei]->d(i), false);
      }

      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
        if (-1==nb.d(j)) continue;  // nothing can be added to this set

        node_handle rec = recFire(nb.d(i), Rp->d(jz));

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
          node_handle acc = mddUnion->compute(nb.d(j), rec);
          resF->unlinkNode(rec);
          if (acc != nb.d(j)) {
            resF->unlinkNode(nb.d(j));
            nb.d_ref(j) = acc;
          } else {
            resF->unlinkNode(acc);
            updated = false;
          }
        }

        if (updated) queue->add(j);
      } // for j
    } // for all events, ei
  } // while there are indexes to explore

  // cleanup
  unpacked_node::recycle(Rp);
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) unpacked_node::recycle(Ru[ei]);
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
  compute_table::entry_key* Key = findResult(mdd, mxd, result);
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
  const int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::useUnpackedNode();
  if (mddLevel < rLevel) {
    A->initRedundant(arg1F, rLevel, mdd, true);
  } else {
    A->initFromNode(arg1F, mdd, true);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i=0; i<rSize; i++) {
      nb->d_ref(i) = recFire(A->d(i), mxd);
    }

  } else {
    // 
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (int i=0; i<rSize; i++) nb->d_ref(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::useUnpackedNode();
    unpacked_node *Rp = unpacked_node::useUnpackedNode();
    if (mxdLevel < 0) {
      Ru->initRedundant(arg2F, rLevel, mxd, false);
    } else {
      Ru->initFromNode(arg2F, mxd, false);
    }

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      if (0==A->d(i))   continue; 
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(arg2F, rLevel, i, Ru->d(iz), false);
      } else {
        Rp->initFromNode(arg2F, Ru->d(iz), false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = recFire(A->d(i), Rp->d(jz));
        if (0==newstates) continue;
        if (0==nb->d(j)) {
          nb->d_ref(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        int oldj = nb->d(j);
        nb->d_ref(j) = mddUnion->compute(newstates, oldj);
        resF->unlinkNode(oldj);
        resF->unlinkNode(newstates);
      } // for j
  
    } // for i

    unpacked_node::recycle(Rp);
    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);

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
    bckwd_dfs_by_events_mt(const satpregen_opname* opcode,
    satpregen_opname::pregen_relation* rel);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::bckwd_dfs_by_events_mt::bckwd_dfs_by_events_mt(
  const satpregen_opname* opcode,
  satpregen_opname::pregen_relation* rel)
  : common_dfs_by_events_mt(opcode, rel)
{
}

void MEDDLY::bckwd_dfs_by_events_mt::saturateHelper(unpacked_node& nb)
{
  int nEventsAtThisLevel = rel->lengthForLevel(nb.getLevel());
  if (0 == nEventsAtThisLevel) return;

  // Initialize mxd readers, note we might skip the unprimed level
  node_handle* events = rel->arrayForLevel(nb.getLevel());
  unpacked_node** Ru = new unpacked_node*[nEventsAtThisLevel];
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    Ru[ei] = unpacked_node::useUnpackedNode();
    int eventLevel = arg2F->getNodeLevel(events[ei]);
    MEDDLY_DCASSERT(ABS(eventLevel) == nb.getLevel());
    if (eventLevel<0) {
      Ru[ei]->initRedundant(arg2F, nb.getLevel(), events[ei], true);
    } else {
      Ru[ei]->initFromNode(arg2F, events[ei], true);
    }
  }
  unpacked_node* Rp = unpacked_node::useUnpackedNode();

  // indexes to explore
  charbuf* expl = useCharBuf(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) expl->data[i] = 2;
  bool repeat = true;

  // explore 
  while (repeat) {
    // "advance" the explore list
    for (int i=0; i<nb.getSize(); i++) if (expl->data[i]) expl->data[i]--;
    repeat = false;

    // explore all events
    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
      // explore all rows
      for (int iz=0; iz<Ru[ei]->getNNZs(); iz++) {
        int i = Ru[ei]->i(iz);
        // grab column (TBD: build these ahead of time?)
        int dlevel = arg2F->getNodeLevel(Ru[ei]->d(iz));

        if (dlevel == -nb.getLevel()) {
          Rp->initFromNode(arg2F, Ru[ei]->d(iz), false); 
        } else {
          Rp->initIdentity(arg2F, -nb.getLevel(), i, Ru[ei]->d(iz), false);
        }

        for (int jz=0; jz<Rp->getNNZs(); jz++) {
          int j = Rp->i(jz);
          if (0==expl->data[j]) continue;
          if (0==nb.d(j))       continue;
          // We have an i->j edge to explore
          node_handle rec = recFire(nb.d(j), Rp->d(jz));

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
            node_handle acc = mddUnion->compute(nb.d(i), rec);
            resF->unlinkNode(rec);
            if (acc != nb.d(i)) {
              resF->unlinkNode(nb.d(i));
              nb.d_ref(i) = acc;
            } else {
              resF->unlinkNode(acc);
              updated = false;
            }
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
  unpacked_node::recycle(Rp);
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) unpacked_node::recycle(Ru[ei]);
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
  compute_table::entry_key* Key = findResult(mdd, mxd, result);
  if (0==Key) return result;

  // check if mxd and mdd are at the same level
  const int mddLevel = arg1F->getNodeLevel(mdd);
  const int mxdLevel = arg2F->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::useUnpackedNode();
  if (mddLevel < rLevel) {
    A->initRedundant(arg1F, rLevel, mdd, true);
  } else {
    A->initFromNode(arg1F, mdd, true);
  }


  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (int i=0; i<rSize; i++) {
      nb->d_ref(i) = recFire(A->d(i), mxd);
    }
  } else {
    // 
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (int i=0; i<rSize; i++) nb->d_ref(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::useUnpackedNode();
    unpacked_node *Rp = unpacked_node::useUnpackedNode();
    if (mxdLevel < 0) {
      Ru->initRedundant(arg2F, rLevel, mxd, false);
    } else {
      Ru->initFromNode(arg2F, mxd, false);
    }

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(arg2F, rLevel, i, Ru->d(iz), false);
      } else {
        Rp->initFromNode(arg2F, Ru->d(iz), false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
        if (0==A->d(j))   continue; 
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = recFire(A->d(j), Rp->d(jz));
        if (0==newstates) continue;
        if (0==nb->d(i)) {
          nb->d_ref(i) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        const node_handle oldi = nb->d(i);
        nb->d_ref(i) = mddUnion->compute(newstates, oldi);
        resF->unlinkNode(oldi);
        resF->unlinkNode(newstates);
      } // for j
  
    } // for i

    unpacked_node::recycle(Rp);
    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);

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
    virtual specialized_operation* buildOperation(arguments* a) const;
};

MEDDLY::fb_saturation_opname::fb_saturation_opname(bool fwd)
 : satpregen_opname(fwd ? "SaturationFwd" : "SaturationBack")
{
  forward = fwd;
}

MEDDLY::specialized_operation*
MEDDLY::fb_saturation_opname::buildOperation(arguments* a) const
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


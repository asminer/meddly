
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
#include "sat_otf.h"
#include <typeinfo> // for "bad_cast" exception

namespace MEDDLY {
  class otfsat_by_events_opname;
  class otfsat_by_events_op;

  class common_otf_dfs_by_events_mt;
  class forwd_otf_dfs_by_events_mt;
  class bckwd_otf_dfs_by_events_mt;

  class fb_otf_saturation_opname;
};


// ******************************************************************
// *                                                                *
// *                    otfsat_by_events_opname  class              *
// *                                                                *
// ******************************************************************

/** Simple class to keep compute table happy.
*/
class MEDDLY::otfsat_by_events_opname : public unary_opname {
  static otfsat_by_events_opname* instance;
  public:
    otfsat_by_events_opname();

    static const otfsat_by_events_opname* getInstance();
 
};

MEDDLY::otfsat_by_events_opname* MEDDLY::otfsat_by_events_opname::instance = 0;

MEDDLY::otfsat_by_events_opname::otfsat_by_events_opname()
 : unary_opname("Otf_Saturate_by_events")
{
}

const MEDDLY::otfsat_by_events_opname* MEDDLY::otfsat_by_events_opname::getInstance()
{
  if (0==instance) instance = new otfsat_by_events_opname;
  return instance;
}

// ******************************************************************
// *                                                                *
// *                      otfsat_by_events_op  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::otfsat_by_events_op : public unary_operation {
    common_otf_dfs_by_events_mt* parent;
  public:
    otfsat_by_events_op(common_otf_dfs_by_events_mt* p,
      expert_forest* argF, expert_forest* resF);
    virtual ~otfsat_by_events_op();

    node_handle saturate(node_handle mdd);
    node_handle saturate(node_handle mdd, int level);

#ifdef OLD_OP_CT
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
#endif
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData, bool key_only) const;
#endif

  protected:
    inline compute_table::entry_key* 
    findSaturateResult(node_handle a, int level, node_handle& b) {
#ifdef OLD_OP_CT
      compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      if (argF->isFullyReduced()) CTsrch->writeI(level);
#ifdef OLD_OP_CT
      compute_table::entry_result& cacheFind = CT0->find(CTsrch);
      if (!cacheFind) return CTsrch;
      b = resF->linkNode(cacheFind.readN()); 
#else
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      b = resF->linkNode(CTresult[0].readN()); 
#endif
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveSaturateResult(compute_table::entry_key* Key,
      node_handle a, node_handle b) 
    {
#ifdef OLD_OP_CT
      argF->cacheNode(a);
      resF->cacheNode(b);
      static compute_table::entry_result result(1);
      result.reset();
      result.writeN(b);
      CT0->addEntry(Key, result);
#else
      CTresult[0].reset();
      CTresult[0].writeN(b);
      CT0->addEntry(Key, CTresult[0]);
#endif
      return b;
    }
};


// ******************************************************************
// *                                                                *
// *            common_otf_dfs_by_events_mt  class                  *
// *                                                                *
// ******************************************************************

class MEDDLY::common_otf_dfs_by_events_mt : public specialized_operation {
  public:
    common_otf_dfs_by_events_mt(const satotf_opname* opcode,
      satotf_opname::otf_relation* rel);
    virtual ~common_otf_dfs_by_events_mt();

#ifdef OLD_OP_CT
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle*);
#endif
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData, bool key_only) const;
#endif
    virtual void compute(const dd_edge& a, dd_edge &c);
    virtual void saturateHelper(unpacked_node& mdd) = 0;

  protected:
    inline compute_table::entry_key* 
    findResult(node_handle a, node_handle b, node_handle &c) 
    {
#ifdef OLD_OP_CT
      compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CTsrch->writeN(b);
#ifdef OLD_OP_CT
      compute_table::entry_result& cacheFind = CT0->find(CTsrch);
      if (!cacheFind) return CTsrch;
      c = resF->linkNode(cacheFind.readN());
#else
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
#endif
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(compute_table::entry_key* Key,
      node_handle a, node_handle b, node_handle c) 
    {
#ifdef OLD_OP_CT
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
      resF->cacheNode(c);
      static compute_table::entry_result result(1);
      result.reset();
      result.writeN(c);
      CT0->addEntry(Key, result);
#else
      CTresult[0].reset();
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
#endif
      return c;
    }

  protected:
    binary_operation* mddUnion;
    binary_operation* mxdIntersection;
    binary_operation* mxdDifference;

    satotf_opname::otf_relation* rel;

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
// *           otfsat_by_events_op  methods                     *
// *                                                                *
// ******************************************************************

#ifdef OLD_OP_CT
MEDDLY::otfsat_by_events_op
::otfsat_by_events_op(common_otf_dfs_by_events_mt* p,
  expert_forest* argF, expert_forest* resF)
  : unary_operation(otfsat_by_events_opname::getInstance(),
    ((argF != 0 && argF->isFullyReduced())? 2: 1), 1, argF, resF)
{
  parent = p;
}
#else
MEDDLY::otfsat_by_events_op
::otfsat_by_events_op(common_otf_dfs_by_events_mt* p,
  expert_forest* argF, expert_forest* resF)
  : unary_operation(otfsat_by_events_opname::getInstance(), 1, argF, resF)
{
  parent = p;

  const char* name = otfsat_by_events_opname::getInstance()->getName();
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
#endif

MEDDLY::otfsat_by_events_op::~otfsat_by_events_op()
{
  removeAllComputeTableEntries();
}

MEDDLY::node_handle MEDDLY::otfsat_by_events_op::saturate(node_handle mdd)
{
  // Saturate
  return saturate(mdd, argF->getNumVariables());
}

MEDDLY::node_handle
MEDDLY::otfsat_by_events_op::saturate(node_handle mdd, int k)
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

  const int sz = argF->getLevelSize(k);               // size
  const int mdd_level = argF->getNodeLevel(mdd);      // mdd level

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

#ifdef OLD_OP_CT

#ifndef USE_NODE_STATUS
bool MEDDLY::otfsat_by_events_op::isStaleEntry(const node_handle* data)
{
  return (argF->isFullyReduced()
  ? (argF->isStale(data[0]) || resF->isStale(data[2]))
  : (argF->isStale(data[0]) || resF->isStale(data[1])));
}
#else
MEDDLY::forest::node_status
MEDDLY::otfsat_by_events_op::getStatusOfEntry(const node_handle* data)
{
  MEDDLY::forest::node_status a = argF->getNodeStatus(data[0]);
  MEDDLY::forest::node_status c =
    resF->getNodeStatus(data[ (argF->isFullyReduced()? 2: 1) ]);

  if (a == MEDDLY::forest::DEAD ||
      c == MEDDLY::forest::DEAD)
    return MEDDLY::forest::DEAD;
  else if (a == MEDDLY::forest::RECOVERABLE ||
      c == MEDDLY::forest::RECOVERABLE)
    return MEDDLY::forest::RECOVERABLE;
  else
    return MEDDLY::forest::ACTIVE;
}
#endif

void MEDDLY::otfsat_by_events_op::discardEntry(const node_handle* data)
{
  if (argF->isFullyReduced()) {
    argF->uncacheNode(data[0]);
    resF->uncacheNode(data[2]);
  } else {
    argF->uncacheNode(data[0]);
    resF->uncacheNode(data[1]);
  }
}

void MEDDLY::otfsat_by_events_op::showEntry(output &strm,
  const node_handle* data, bool key_only) const
{
  if (argF->isFullyReduced()) {
    strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) << "): ";
    if (key_only) strm << "?]";
    else          strm << long(data[2]) << "]";
  } else {
    strm << "[" << getName() << "(" << long(data[0]) << "): ";
    if (key_only) strm << "?]";
    else          strm << long(data[1]) << "]";
  }
}

#endif // OLD_OP_CT

// ******************************************************************
// *                                                                *
// *           common_otf_dfs_by_events_mt  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::common_otf_dfs_by_events_mt::common_otf_dfs_by_events_mt(
  const satotf_opname* opcode,
  satotf_opname::otf_relation* relation)
#ifdef OLD_OP_CT
: specialized_operation(opcode, 2, 1)
#else
: specialized_operation(opcode, 1)
#endif
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
#ifdef OLD_OP_CT
  setAnswerForest(resF);
#else
  compute_table::entry_type* et = new compute_table::entry_type(opcode->getName(), "NN:N");
  et->setForestForSlot(0, arg1F);
  et->setForestForSlot(1, arg2F);
  et->setForestForSlot(3, resF);
  registerEntryType(0, et);
  buildCTs();
#endif
}

MEDDLY::common_otf_dfs_by_events_mt::~common_otf_dfs_by_events_mt()
{
  if (rel->autoDestroy()) delete rel;
  unregisterInForest(arg1F);
  unregisterInForest(arg2F);
  unregisterInForest(resF);
}

#ifdef OLD_OP_CT

#ifndef USE_NODE_STATUS
bool MEDDLY::common_otf_dfs_by_events_mt::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}
#else
MEDDLY::forest::node_status
MEDDLY::common_otf_dfs_by_events_mt::getStatusOfEntry(const node_handle* data)
{
  MEDDLY::forest::node_status a = arg1F->getNodeStatus(data[0]);
  MEDDLY::forest::node_status b = arg2F->getNodeStatus(data[1]);
  MEDDLY::forest::node_status c = resF->getNodeStatus(data[2]);

  if (a == MEDDLY::forest::DEAD ||
      b == MEDDLY::forest::DEAD ||
      c == MEDDLY::forest::DEAD)
    return MEDDLY::forest::DEAD;
  else if (a == MEDDLY::forest::RECOVERABLE ||
      b == MEDDLY::forest::RECOVERABLE ||
      c == MEDDLY::forest::RECOVERABLE)
    return MEDDLY::forest::RECOVERABLE;
  else
    return MEDDLY::forest::ACTIVE;
}
#endif

void MEDDLY::common_otf_dfs_by_events_mt::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void MEDDLY::common_otf_dfs_by_events_mt::showEntry(output &strm,
  const node_handle* data, bool key_only) const
{
  strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) << "): ";
  if (key_only) strm << "?]";
  else          strm << long(data[2]) << "]";
}

#endif // OLD_OP_CT

void MEDDLY::common_otf_dfs_by_events_mt
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
  otfsat_by_events_op* so = new otfsat_by_events_op(this, arg1F, resF);
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
// *       common_otf_dfs_by_events_mt::indexq  methods                 *
// ******************************************************************

MEDDLY::common_otf_dfs_by_events_mt::indexq::indexq()
{
  data = 0;
  size = 0;
  head = NULPTR;
}

MEDDLY::common_otf_dfs_by_events_mt::indexq::~indexq()
{
  free(data);
}

void MEDDLY::common_otf_dfs_by_events_mt::indexq::resize(int sz)
{
  if (sz <= size) return;
  data = (int*) realloc(data, sz * sizeof(int));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);

  for (; size < sz; size++) data[size] = NOTINQ;
}

// ******************************************************************
// *       common_otf_dfs_by_events_mt::charbuf methods             *
// ******************************************************************

MEDDLY::common_otf_dfs_by_events_mt::charbuf::charbuf()
{
  data = 0;
  size = 0;
}

MEDDLY::common_otf_dfs_by_events_mt::charbuf::~charbuf()
{
  free(data);
}

void MEDDLY::common_otf_dfs_by_events_mt::charbuf::resize(int sz)
{
  if (sz <= size) return;
  data = (char*) realloc(data, sz * sizeof(char));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *             forwd_otf_dfs_by_events_mt class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_otf_dfs_by_events_mt : public common_otf_dfs_by_events_mt {
  public:
    forwd_otf_dfs_by_events_mt(const satotf_opname* opcode,
    satotf_opname::otf_relation* rel);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
    void recFireHelper(const int, const int, const node_handle, const node_handle,
        unpacked_node*, unpacked_node*);
};

MEDDLY::forwd_otf_dfs_by_events_mt::forwd_otf_dfs_by_events_mt(
  const satotf_opname* opcode,
  satotf_opname::otf_relation* rel)
  : common_otf_dfs_by_events_mt(opcode, rel)
{
}


void MEDDLY::forwd_otf_dfs_by_events_mt::saturateHelper(unpacked_node& nb)
{
  const int level = nb.getLevel();
  const int nEventsAtThisLevel = rel->getNumOfEvents(level);
  if (0 == nEventsAtThisLevel) return;

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node** Ru = new unpacked_node*[nEventsAtThisLevel];
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    rel->rebuildEvent(level, ei);
    node_handle mxd = rel->getEvent(level, ei);
    if (0==mxd) {
      Ru[ei] = 0;
    } else {
      Ru[ei] = unpacked_node::useUnpackedNode();
      const int eventLevel = arg2F->getNodeLevel(mxd);
      if (ABS(eventLevel) < level || eventLevel < 0) {
        // Takes care of two situations:
        // - skipped unprimed level (due to Fully Reduced)
        // - skipped unprimed and primed levels (due to Fully Identity Reduced)
        Ru[ei]->initRedundant(arg2F, level, mxd, true);
      } else {
        Ru[ei]->initFromNode(arg2F, mxd, true);
      }
    }
  }
  unpacked_node* Rp = unpacked_node::useUnpackedNode();

  //      Node reader auto expands when passed by reference
  //      Node builder can be expanded via a call to unpacked_node::resize()
  //      Queue can be expanded via a call to indexq::resize()

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.d(i)) queue->add(i);
  }

  // explore indexes
  while (!queue->isEmpty()) {
    const int i = queue->remove();

    MEDDLY_DCASSERT(nb.d(i));

    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
      // If event i needs rebuilding,
      //    Rebuild it, and
      //    Update the event reader
      if (rel->rebuildEvent(level, ei)) {
        node_handle mxd = rel->getEvent(level, ei);
        if (0==mxd) {
          if (Ru[ei]) {
            unpacked_node::recycle(Ru[ei]);
            Ru[ei] = 0;
          }
        } else {
          if (0==Ru[ei]) {
            Ru[ei] = unpacked_node::useUnpackedNode();
          }
          const int eventLevel = arg2F->getNodeLevel(mxd);
          if (ABS(eventLevel) < level || eventLevel < 0) {
            // Takes care of two situations:
            // - skipped unprimed level (due to Fully Reduced)
            // - skipped unprimed and primed levels (due to Fully Identity Reduced)
            Ru[ei]->initRedundant(arg2F, level, mxd, true);
          } else {
            Ru[ei]->initFromNode(arg2F, mxd, true);
          }
        }
      }
      // check if row i of the event ei is empty
      if (0 == Ru[ei]) continue;
      MEDDLY_DCASSERT(!Ru[ei]->isExtensible());
      node_handle ei_i = (i < Ru[ei]->getSize())
                        ? Ru[ei]->d(i)
                        : (Ru[ei]->isExtensible() ? Ru[ei]->ext_d() : 0);
      if (0 == ei_i) continue;

      // grab column (TBD: build these ahead of time?)
      const int dlevel = arg2F->getNodeLevel(ei_i);

      if (dlevel == -level) {
        Rp->initFromNode(arg2F, ei_i, false);
      } else {
        Rp->initIdentity(arg2F, -level, i, ei_i, false);
      }

      MEDDLY_DCASSERT(!Rp->isExtensible());

      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        const int j = Rp->i(jz);
        if (j < nb.getSize() && -1==nb.d(j)) continue;  // nothing can be added to this set

        node_handle newstates = recFire(nb.d(i), Rp->d(jz));
        if (newstates == 0) continue;

        // Confirm local state
        rel->confirm(level, j);

        if (j >= nb.getSize()) {
          int oldSize = nb.getSize();
          // resize the node builder, and clear out the new entries
          nb.resize(j+1);
          while(oldSize < nb.getSize()) { nb.d_ref(oldSize++) = 0; }
          // resize the queue, and clear out the new entries
          queue->resize(nb.getSize());
        }

        if (0 == nb.d(j)) {
          nb.d_ref(j) = newstates;
          queue->add(j);
        } else {
          // there's new states and existing states; union them.
          const node_handle oldj = nb.d_ref(j);
          nb.d_ref(j) = mddUnion->compute(nb.d(j), newstates);
          if (oldj != nb.d(j)) queue->add(j);
          resF->unlinkNode(oldj);
          resF->unlinkNode(newstates); 
        }

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
MEDDLY::node_handle MEDDLY::forwd_otf_dfs_by_events_mt::recFire(
  node_handle mdd, node_handle mxd)
{
  //      Node builder expansion
  //          - done
  //      MXD doesnt expand but unconfirmed > confirmed
  //          - so mdd may need to expand to accomodate rec_fire result
  //          - No need: primed and unprimed variables are of the same size
  //      Check if createReduce and handle a unpacked_node of size
  //          larger than the variable size
  //          - No
  //      Can we build a unpacked_node with the size of the primed variable
  //          and save some trouble?
  //          - No need: primed and unprimed variables are of the same size

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

  mxd = arg2F->linkNode(mxd);

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

#if 0
    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);
      if (0==A->d(i)) continue; 
      const node_handle pnode = Ru->d(iz);
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(pnode))) {
        Rp->initIdentity(arg2F, rLevel, i, pnode, false);
      } else {
        Rp->initFromNode(arg2F, pnode, false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        const int j = Rp->i(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = recFire(A->d(i), Rp->d(jz));
        if (0==newstates) continue;

        if (rel->confirm(rLevel, j)) {
          // Newly confirmed local state
          // Node builder must expand to accomodate j
          if (j >= nb->getSize()) {
            int oldSize = nb->getSize();
            // resize the node builder, and clear out the new entries
            nb->resize(j+1);
            while(oldSize < nb->getSize()) { nb->d_ref(oldSize++) = 0; }
          }
        }

        if (0==nb->d(j)) {
          nb->d_ref(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        const int oldj = nb->d(j);
        nb->d_ref(j) = mddUnion->compute(newstates, oldj);
        resF->unlinkNode(oldj);
        resF->unlinkNode(newstates);
      } // for j
    } // for i
#else
    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);
      if (0==A->d(i)) continue; 
      recFireHelper(i, rLevel, Ru->d(iz), A->d(i), Rp, nb);
    }
    // loop over the extensible portion of mxd (if any)
    MEDDLY_DCASSERT(!Ru->isExtensible());
    if (Ru->isExtensible()) {
      const node_handle pnode = Ru->ext_d();
      for (int i = Ru->ext_i()+1; i < A->getSize(); i++) {
        if (0 == A->d(i)) continue;
        recFireHelper(i, rLevel, pnode, A->d(i), Rp, nb);
      }
    }
#endif

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


void MEDDLY::forwd_otf_dfs_by_events_mt::recFireHelper(
  const int i,
  const int rLevel,
  const node_handle Ru_i,
  const node_handle A_i,
  unpacked_node *Rp,
  unpacked_node* nb)
{
  if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru_i))) {
    Rp->initIdentity(arg2F, rLevel, i, Ru_i, false);
  } else {
    Rp->initFromNode(arg2F, Ru_i, false);
  }

  MEDDLY_DCASSERT(!Rp->isExtensible());

  // loop over mxd "columns"
  for (int jz=0; jz<Rp->getNNZs(); jz++) {
    const int j = Rp->i(jz);
    // ok, there is an i->j "edge".
    // determine new states to be added (recursively)
    // and add them
    node_handle newstates = recFire(A_i, Rp->d(jz));
    if (0==newstates) continue;

    // Confirm local state
    rel->confirm(rLevel, j);

    if (j >= nb->getSize()) {
      int oldSize = nb->getSize();
      // resize the node builder, and clear out the new entries
      nb->resize(j+1);
      while(oldSize < nb->getSize()) { nb->d_ref(oldSize++) = 0; }
    }

    if (0==nb->d(j)) {
      nb->d_ref(j) = newstates;
    } else {
      // there's new states and existing states; union them.
      const node_handle oldj = nb->d(j);
      nb->d_ref(j) = mddUnion->compute(oldj, newstates);
      resF->unlinkNode(oldj);
      resF->unlinkNode(newstates);
    }
  } // for j
}



// ******************************************************************
// *                                                                *
// *                   fb_otf_saturation_opname class               *
// *                                                                *
// ******************************************************************

class MEDDLY::fb_otf_saturation_opname : public satotf_opname {
    bool forward;
  public:
    fb_otf_saturation_opname(bool fwd);
    virtual specialized_operation* buildOperation(arguments* a) const;
};

MEDDLY::fb_otf_saturation_opname::fb_otf_saturation_opname(bool fwd)
 : satotf_opname(fwd ? "OtfSaturationFwd" : "OtfSaturationBack")
{
  forward = fwd;
}

MEDDLY::specialized_operation*
MEDDLY::fb_otf_saturation_opname::buildOperation(arguments* a) const
{
  otf_relation* rel = dynamic_cast<otf_relation*>(a);
  if (0==rel) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);

  //
  // No sanity checks needed here; we did them already when constructing a.
  //

  MEDDLY::specialized_operation* op = 0;
  if (forward)
    op = new forwd_otf_dfs_by_events_mt(this, rel);
  else {
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    // op = new bckwd_otf_dfs_by_events_mt(this, rel);
  }

  // Do we need to delete rel here?
  // No, if needed, do this in the destructor for op.

  return op;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::satotf_opname* MEDDLY::initOtfSaturationForward()
{
  return new fb_otf_saturation_opname(true);
}

MEDDLY::satotf_opname* MEDDLY::initOtfSaturationBackward()
{
  return new fb_otf_saturation_opname(false);
}


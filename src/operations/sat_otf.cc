
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

    virtual bool isStaleEntry(const node_handle* entryData);
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData) const;

  protected:
    inline compute_table::search_key* 
    findSaturateResult(node_handle a, int level, node_handle& b) {
      compute_table::search_key* CTsrch = useCTkey();
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->reset();
      CTsrch->writeNH(a);
      if (argF->isFullyReduced()) CTsrch->write(level);
      compute_table::search_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      b = resF->linkNode(cacheFind.readNH()); 
      doneCTkey(CTsrch);
      return 0;
    }
    inline node_handle saveSaturateResult(compute_table::search_key* Key,
      node_handle a, node_handle b) 
    {
      argF->cacheNode(a);
      compute_table::entry_builder &entry = CT->startNewEntry(Key);
      entry.writeResultNH(resF->cacheNode(b));
      CT->addEntry();
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

    virtual bool isStaleEntry(const node_handle* entryData);
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData) const;
    virtual void compute(const dd_edge& a, dd_edge &c);
    virtual void saturateHelper(node_builder& mdd) = 0;

  protected:
    inline compute_table::search_key* 
    findResult(node_handle a, node_handle b, node_handle &c) 
    {
      compute_table::search_key* CTsrch = useCTkey();
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->reset();
      CTsrch->writeNH(a);
      CTsrch->writeNH(b);
      compute_table::search_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      c = resF->linkNode(cacheFind.readNH());
      doneCTkey(CTsrch);
      return 0;
    }
    inline node_handle saveResult(compute_table::search_key* Key,
      node_handle a, node_handle b, node_handle c) 
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
      compute_table::entry_builder &entry = CT->startNewEntry(Key);
      entry.writeResultNH(resF->cacheNode(c));
      CT->addEntry();
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
};

// ******************************************************************
// *                                                                *
// *           otfsat_by_events_op  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::otfsat_by_events_op
::otfsat_by_events_op(common_otf_dfs_by_events_mt* p,
  expert_forest* argF, expert_forest* resF)
  : unary_operation(otfsat_by_events_opname::getInstance(),
    ((argF != 0 && argF->isFullyReduced())? 2: 1), 1, argF, resF)
{
  parent = p;

}

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
  compute_table::search_key* Key = findSaturateResult(mdd, k, n);
  if (0==Key) return n;

  int sz = argF->getLevelSize(k);               // size
  int mdd_level = argF->getNodeLevel(mdd);      // mdd level

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d, mdd_level: %d\n",
      mdd, k, sz, mdd_level);
#endif

  node_builder& nb = resF->useNodeBuilder(k, sz);
  node_reader* mddDptrs =
    (mdd_level < k)
    ? argF->initRedundantReader(k, mdd, true)
    : argF->initNodeReader(mdd, true);
  for (int i=0; i<sz; i++) {
    nb.d(i) = mddDptrs->d(i) ? saturate(mddDptrs->d(i), k-1) : 0;
  }
  node_reader::recycle(mddDptrs);
  parent->saturateHelper(nb);
  n = resF->createReducedNode(-1, nb);

  // save in compute table
  saveSaturateResult(Key, mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  return n;
}

bool MEDDLY::otfsat_by_events_op::isStaleEntry(const node_handle* data)
{
  return (argF->isFullyReduced()
  ? (argF->isStale(data[0]) || resF->isStale(data[2]))
  : (argF->isStale(data[0]) || resF->isStale(data[1])));
}

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
  const node_handle* data) const
{
  if (argF->isFullyReduced()) {
    strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) << "): " << long(data[2]) << "]";
  } else {
    strm << "[" << getName() << "(" << long(data[0]) << "): " << long(data[1]) << "]";
  }
}


// ******************************************************************
// *                                                                *
// *           common_otf_dfs_by_events_mt  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::common_otf_dfs_by_events_mt::common_otf_dfs_by_events_mt(
  const satotf_opname* opcode,
  satotf_opname::otf_relation* relation)
: specialized_operation(opcode, 2, 1)
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
  setAnswerForest(resF);
}

MEDDLY::common_otf_dfs_by_events_mt::~common_otf_dfs_by_events_mt()
{
  if (rel->autoDestroy()) delete rel;
  unregisterInForest(arg1F);
  unregisterInForest(arg2F);
  unregisterInForest(resF);
}

bool MEDDLY::common_otf_dfs_by_events_mt::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}

void MEDDLY::common_otf_dfs_by_events_mt::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void MEDDLY::common_otf_dfs_by_events_mt::showEntry(output &strm,
  const node_handle* data) const
{
  strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) << "): " << long(data[2]) << "]";
}

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
#if 0
  node_handle cnode = so->saturate(a.getNode());
  c.set(cnode);
#else
  node_handle prev = 0;
  node_handle curr = a.getNode();
  while (prev != curr) {
    prev = curr;
    curr = so->saturate(prev);
  }
  c.set(curr);
#endif

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
    throw error(error::INSUFFICIENT_MEMORY);

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
    throw error(error::INSUFFICIENT_MEMORY);
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
    virtual void saturateHelper(node_builder& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::forwd_otf_dfs_by_events_mt::forwd_otf_dfs_by_events_mt(
  const satotf_opname* opcode,
  satotf_opname::otf_relation* rel)
  : common_otf_dfs_by_events_mt(opcode, rel)
{
}


void MEDDLY::forwd_otf_dfs_by_events_mt::saturateHelper(node_builder& nb)
{
  int level = nb.getLevel();
  int nEventsAtThisLevel = rel->getNumOfEvents(level);
  if (0 == nEventsAtThisLevel) return;

  // Initialize mxd readers, note we might skip the unprimed level
  node_reader** Ru = new node_reader*[nEventsAtThisLevel];
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    rel->rebuildEvent(level, ei);
    node_handle mxd = rel->getEvent(level, ei);
    if (0==mxd) {
      printf("ZERO EVENT FOUND!\n");
      Ru[ei] = 0;
    } else {
      printf("NON-ZERO EVENT FOUND!\n");
      int eventLevel = arg2F->getNodeLevel(mxd);
      MEDDLY_DCASSERT(ABS(eventLevel) == level);
      Ru[ei] = (eventLevel<0)
        ? arg2F->initRedundantReader(level, mxd, true)
        : arg2F->initNodeReader(mxd, true);
    }
  }
  node_reader* Rp = node_reader::useReader();

  //      Node reader auto expands when passed by reference
  //      Node builder can be expanded via a call to node_builder::resize()
  //      Queue can be expanded via a call to indexq::resize()

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
      // If event i needs rebuilding,
      //    Rebuild it, and
      //    Update the event reader
      if (rel->rebuildEvent(level, ei)) {
        node_handle mxd = rel->getEvent(level, ei);
        if (0==mxd) {
          printf("ZERO EVENT FOUND!\n");
          if (Ru[ei]) {
            node_reader::recycle(Ru[ei]);
            Ru[ei] = 0;
          }
        } else {
          printf("NON-ZERO EVENT FOUND!\n");
          int eventLevel = arg2F->getNodeLevel(mxd);
          MEDDLY_DCASSERT(ABS(eventLevel) == level);
          if (Ru[ei]){
            if (eventLevel<0)
              arg2F->initRedundantReader(*Ru[ei], level, mxd, true);
            else
              arg2F->initNodeReader(*Ru[ei], mxd, true);
          } else {
            Ru[ei] = (eventLevel<0)
              ? arg2F->initRedundantReader(level, mxd, true)
              : arg2F->initNodeReader(mxd, true);
          }
        }
      }
      // check if row i of the event ei is empty
      if (0 == Ru[ei] || i >= Ru[ei]->getSize() || 0 == Ru[ei]->d(i)) continue;

      // grab column (TBD: build these ahead of time?)
      int dlevel = arg2F->getNodeLevel(Ru[ei]->d(i));

      if (dlevel == -level) {
        arg2F->initNodeReader(*Rp, Ru[ei]->d(i), false);
      } else {
        arg2F->initIdentityReader(*Rp, -level, i, Ru[ei]->d(i), false);
      }

      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
        if (j < nb.getSize() && -1==nb.d(j)) continue;  // nothing can be added to this set

        node_handle rec = recFire(nb.d(i), Rp->d(jz));

        if (rec == 0) continue;

        if (rel->confirm(level, j)) {
          // Newly confirmed local state
          // Node builder must expand to accomodate j
          if (j >= nb.getSize()) {
            int oldSize = nb.getSize();
            // resize the node builder, and clear out the new entries
            nb.resize(j+1);
            while(oldSize < nb.getSize()) { nb.d(oldSize++) = 0; }
            // resize the queue, and clear out the new entries
            queue->resize(nb.getSize());
          }
        }

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
          node_handle acc = mddUnion->compute(nb.d(j), rec);
          resF->unlinkNode(rec);
          if (acc != nb.d(j)) {
            resF->unlinkNode(nb.d(j));
            nb.d(j) = acc;
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
  node_reader::recycle(Rp);
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) node_reader::recycle(Ru[ei]);
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
  //      Check if createReduce and handle a node_builder of size
  //          larger than the variable size
  //          - No
  //      Can we build a node_builder with the size of the primed variable
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
  compute_table::search_key* Key = findResult(mdd, mxd, result);
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
  int mddLevel = arg1F->getNodeLevel(mdd);
  int mxdLevel = arg2F->getNodeLevel(mxd);
  int rLevel = MAX(ABS(mxdLevel), mddLevel);
  int rSize = resF->getLevelSize(rLevel);
  node_builder& nb = resF->useNodeBuilder(rLevel, rSize);

  // Initialize mdd reader
  node_reader* A = (mddLevel < rLevel)
    ? arg1F->initRedundantReader(rLevel, mdd, true)
    : arg1F->initNodeReader(mdd, true);

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i=0; i<rSize; i++) {
      nb.d(i) = recFire(A->d(i), mxd);
    }

  } else {
    // 
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (int i=0; i<rSize; i++) nb.d(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    node_reader* Ru = (mxdLevel < 0)
      ? arg2F->initRedundantReader(rLevel, mxd, false)
      : arg2F->initNodeReader(mxd, false);
    node_reader* Rp = node_reader::useReader();

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      if (0==A->d(i))   continue; 
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->d(iz)))) {
        arg2F->initIdentityReader(*Rp, rLevel, i, Ru->d(iz), false);
      } else {
        arg2F->initNodeReader(*Rp, Ru->d(iz), false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = recFire(A->d(i), Rp->d(jz));
        if (0==newstates) continue;

        if (rel->confirm(rLevel, j)) {
          // Newly confirmed local state
          // Node builder must expand to accomodate j
          if (j >= nb.getSize()) {
            int oldSize = nb.getSize();
            // resize the node builder, and clear out the new entries
            nb.resize(j+1);
            while(oldSize < nb.getSize()) { nb.d(oldSize++) = 0; }
          }
        }

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
  
    } // for i

    node_reader::recycle(Rp);
    node_reader::recycle(Ru);
  } // else

  // cleanup mdd reader
  node_reader::recycle(A);

  saturateHelper(nb);
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



#if 0

// ******************************************************************
// *                                                                *
// *             bckwd_otf_dfs_by_events_mt class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_otf_dfs_by_events_mt : public common_otf_dfs_by_events_mt {
  public:
    bckwd_otf_dfs_by_events_mt(const satotf_opname* opcode,
    satotf_opname::otf_relation* rel);
  protected:
    virtual void saturateHelper(node_builder& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::bckwd_otf_dfs_by_events_mt::bckwd_otf_dfs_by_events_mt(
  const satotf_opname* opcode,
  satotf_opname::otf_relation* rel)
  : common_otf_dfs_by_events_mt(opcode, rel)
{
}

void MEDDLY::bckwd_otf_dfs_by_events_mt::saturateHelper(node_builder& nb)
{
  int nEventsAtThisLevel = rel->getNumOfEvents(nb.getLevel());
  if (0 == nEventsAtThisLevel) return;

  // Initialize mxd readers, note we might skip the unprimed level
  node_handle* events = rel->get(nb.getLevel());
  node_reader** Ru = new node_reader*[nEventsAtThisLevel];
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    int eventLevel = arg2F->getNodeLevel(events[ei]);
    MEDDLY_DCASSERT(ABS(eventLevel) == nb.getLevel());
    Ru[ei] = (eventLevel<0)
      ? arg2F->initRedundantReader(nb.getLevel(), events[ei], false)
      : arg2F->initNodeReader(events[ei], false);
  }
  node_reader* Rp = node_reader::useReader();

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
          arg2F->initNodeReader(*Rp, Ru[ei]->d(iz), false); 
        } else {
          arg2F->initIdentityReader(*Rp, -nb.getLevel(), i,
            Ru[ei]->d(iz), false);
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
            nb.d(i) = rec;
          }
          else if (-1 == rec) {
            resF->unlinkNode(nb.d(i));
            nb.d(i) = -1;
          } 
          else {
            node_handle acc = mddUnion->compute(nb.d(i), rec);
            resF->unlinkNode(rec);
            if (acc != nb.d(i)) {
              resF->unlinkNode(nb.d(i));
              nb.d(i) = acc;
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
  node_reader::recycle(Rp);
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) node_reader::recycle(Ru[ei]);
  delete[] Ru;
  recycle(expl);
}

MEDDLY::node_handle MEDDLY::bckwd_otf_dfs_by_events_mt::recFire(node_handle mdd,
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
  compute_table::search_key* Key = findResult(mdd, mxd, result);
  if (0==Key) return result;

  // check if mxd and mdd are at the same level
  int mddLevel = arg1F->getNodeLevel(mdd);
  int mxdLevel = arg2F->getNodeLevel(mxd);
  int rLevel = MAX(ABS(mxdLevel), mddLevel);
  int rSize = resF->getLevelSize(rLevel);
  node_builder& nb = resF->useNodeBuilder(rLevel, rSize);

  // Initialize mdd reader
  node_reader* A = (mddLevel < rLevel)
    ? arg1F->initRedundantReader(rLevel, mdd, true)
    : arg1F->initNodeReader(mdd, true);

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (int i=0; i<rSize; i++) {
      nb.d(i) = recFire(A->d(i), mxd);
    }
  } else {
    // 
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // clear out result (important!)
    for (int i=0; i<rSize; i++) nb.d(i) = 0;

    // Initialize mxd readers, note we might skip the unprimed level
    node_reader* Ru = (mxdLevel < 0)
      ? arg2F->initRedundantReader(rLevel, mxd, false)
      : arg2F->initNodeReader(mxd, false);
    node_reader* Rp = node_reader::useReader();

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->d(iz)))) {
        arg2F->initIdentityReader(*Rp, rLevel, i, Ru->d(iz), false);
      } else {
        arg2F->initNodeReader(*Rp, Ru->d(iz), false);
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
        if (0==nb.d(i)) {
          nb.d(i) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        node_handle oldi = nb.d(i);
        nb.d(i) = mddUnion->compute(newstates, oldi);
        resF->unlinkNode(oldi);
        resF->unlinkNode(newstates);
      } // for j
  
    } // for i

    node_reader::recycle(Rp);
    node_reader::recycle(Ru);
  } // else

  // cleanup mdd reader
  node_reader::recycle(A);

  saturateHelper(nb);
  result = resF->createReducedNode(-1, nb);
#ifdef TRACE_ALL_OPS
  printf("computed recFire(%d, %d) = %d\n", mdd, mxd, result);
#endif
  return saveResult(Key, mdd, mxd, result); 
}

#endif
















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
  if (0==rel) throw error(error::INVALID_ARGUMENT);

  //
  // No sanity checks needed here; we did them already when constructing a.
  //

  MEDDLY::specialized_operation* op = 0;
  if (forward)
    op = new forwd_otf_dfs_by_events_mt(this, rel);
  else {
    throw MEDDLY::error::NOT_IMPLEMENTED;
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

MEDDLY::satotf_opname* MEDDLY::initOtfSaturationForward(const settings &s)
{
  return new fb_otf_saturation_opname(true);
}

MEDDLY::satotf_opname* MEDDLY::initOtfSaturationBackward(const settings &s)
{
  return new fb_otf_saturation_opname(false);
}


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
#include <set>

#include "../operators.h"
#include "../minterms.h"
#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../oper_satur.h"
#include "../sat_relations.h"
#include "../ops_builtin.h"
#include "../forest_levels.h"

namespace MEDDLY {
    class otfsat_by_events_op;

    class common_otf_dfs_by_events_mt;
    class forwd_otf_dfs_by_events_mt;
    class bckwd_otf_dfs_by_events_mt;
};



// ******************************************************************
// *                                                                *
// *                      otfsat_by_events_op  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::otfsat_by_events_op : public saturation_operation {
    common_otf_dfs_by_events_mt* parent;
  public:
    otfsat_by_events_op(common_otf_dfs_by_events_mt* p,
        const char* name, forest* argF, forest* resF);
    virtual ~otfsat_by_events_op();

    virtual void compute(const dd_edge& in, dd_edge& out);
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
// *            common_otf_dfs_by_events_mt  class                  *
// *                                                                *
// ******************************************************************

class MEDDLY::common_otf_dfs_by_events_mt : public saturation_operation {
  public:
    common_otf_dfs_by_events_mt(const char* name, forest* argF,
            otf_relation* rel, forest* resF);
    virtual ~common_otf_dfs_by_events_mt();

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

    otf_relation* rel;
    forest* relF;

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
        inline void add(int i) {
            MEDDLY_CHECK_RANGE(0u, (unsigned)i, size);
          if (NOTINQ != data[i]) return;
          if (NULPTR == head) {
            // empty list
            head = i;
          } else {
            // not empty list
              MEDDLY_CHECK_RANGE(0u, (unsigned)tail, size);
            data[tail] = i;
          }
          tail = i;
          data[i] = NULPTR;
        }
        inline int remove() {
            MEDDLY_CHECK_RANGE(0u, (unsigned)head, size);
          int ans = head;
          head = data[head];
          data[ans] = NOTINQ;
          MEDDLY_CHECK_RANGE(0u, (unsigned)ans, size);
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
// *           otfsat_by_events_op  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::otfsat_by_events_op
::otfsat_by_events_op(common_otf_dfs_by_events_mt* p,
    const char* name, forest* argF, forest* resF)
  : saturation_operation(name, 1, argF, resF)
{
  parent = p;

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

MEDDLY::otfsat_by_events_op::~otfsat_by_events_op()
{
  // removeAllComputeTableEntries();
}

void MEDDLY::otfsat_by_events_op::compute(const dd_edge& in, dd_edge& out)
{
  // Saturate
  out.set( saturate(in.getNode(), argF->getMaxLevelIndex()) );
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
  ct_entry_key* Key = findSaturateResult(mdd, k, n);
  if (0==Key) return n;

  const unsigned sz = unsigned(argF->getLevelSize(k));    // size
  const int mdd_level = argF->getNodeLevel(mdd);          // mdd level

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d, mdd_level: %d\n",
      mdd, k, sz, mdd_level);
#endif

  unpacked_node* nb = unpacked_node::newWritable(resF, k, sz, FULL_ONLY);
  // Initialize mdd reader
  unpacked_node *mddDptrs = unpacked_node::New(argF, FULL_ONLY);
  if (mdd_level < k) {
    mddDptrs->initRedundant(k, mdd);
  } else {
    mddDptrs->initFromNode(mdd);
  }

  // Do computation
  for (unsigned i=0; i<sz; i++) {
    if (mddDptrs->down(i)) {
      nb->setFull(i, saturate(mddDptrs->down(i), k-1));
    }
    // nb->d_ref(i) = mddDptrs->down(i) ? saturate(mddDptrs->down(i), k-1) : 0;
  }

  // Cleanup
  unpacked_node::Recycle(mddDptrs);

  parent->saturateHelper(*nb);
  edge_value ev;
  resF->createReducedNode(nb, ev, n);
  MEDDLY_DCASSERT(ev.isVoid());

  // save in compute table
  saveSaturateResult(Key, mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  return n;
}


// ******************************************************************
// *                                                                *
// *           common_otf_dfs_by_events_mt  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::common_otf_dfs_by_events_mt::common_otf_dfs_by_events_mt(
  const char* name, forest* argF, otf_relation* relation, forest* resF)
: saturation_operation(name, 1, argF, resF)
{
    if (!argF || !relation || !resF) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }

    rel = relation;
    relF = relation->getRelForest();

    // for now, anyway, inset and outset must be same forest
    if (argF != resF) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }

    // Check for same domain
    if  (
            (argF->getDomain() != relF->getDomain()) ||
            (resF->getDomain() != relF->getDomain())
        )
    {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    mddUnion = nullptr;
    mxdIntersection = nullptr;
    mxdDifference = nullptr;

    freeqs = 0;
    freebufs = 0;

    registerInForest(relF);

    ct_entry_type* et = new ct_entry_type(name, "NN:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(1, relF);
    et->setForestForSlot(3, resF);
    registerEntryType(0, et);
    buildCTs();
}

MEDDLY::common_otf_dfs_by_events_mt::~common_otf_dfs_by_events_mt()
{
    unregisterInForest(relF);
    delete rel;
}

void MEDDLY::common_otf_dfs_by_events_mt
::compute(const dd_edge &a, dd_edge &c)
{
  // Initialize operations
  mddUnion = build(UNION, resF, resF, resF);
  MEDDLY_DCASSERT(mddUnion);

  mxdIntersection = build(INTERSECTION, relF, relF, relF);
  MEDDLY_DCASSERT(mxdIntersection);

  mxdDifference = build(DIFFERENCE, relF, relF, relF);
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
  otfsat_by_events_op* so = new otfsat_by_events_op(this, "Otf_Saturate_by_events", argF, resF);
  so->compute(a, c);

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

void MEDDLY::common_otf_dfs_by_events_mt::indexq::resize(unsigned sz)
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

void MEDDLY::common_otf_dfs_by_events_mt::charbuf::resize(unsigned sz)
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
    forwd_otf_dfs_by_events_mt(const char* name, forest* argF,
            otf_relation* rel, forest* resF);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
    void recFireHelper(const unsigned, const int, const node_handle, const node_handle,
        unpacked_node*, unpacked_node*);
};

MEDDLY::forwd_otf_dfs_by_events_mt::forwd_otf_dfs_by_events_mt(
    const char* name, forest* argF, otf_relation* rel, forest* resF)
  : common_otf_dfs_by_events_mt(name, argF, rel, resF)
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
    const dd_edge& mxd = rel->getEvent(level, ei);
    if (0==mxd.getNode()) {
      Ru[ei] = 0;
    } else {
      Ru[ei] = unpacked_node::New(relF, FULL_ONLY);
      const int eventLevel = mxd.getLevel();
      if (ABS(eventLevel) < level || eventLevel < 0) {
        // Takes care of two situations:
        // - skipped unprimed level (due to Fully Reduced)
        // - skipped unprimed and primed levels (due to Fully Identity Reduced)
        Ru[ei]->initRedundant(level, mxd.getNode());
      } else {
        Ru[ei]->initFromNode(mxd.getNode());
      }
    }
  }
  unpacked_node* Rp = unpacked_node::New(relF, SPARSE_ONLY);

  dd_edge nbdj(resF), newst(resF);

  //      Node reader auto expands when passed by reference
  //      Node builder can be expanded via a call to unpacked_node::resize()
  //      Queue can be expanded via a call to indexq::resize()

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (unsigned i = 0; i < nb.getSize(); i++) {
    if (nb.down(i)) queue->add(int(i));
  }

  // explore indexes
  while (!queue->isEmpty()) {
    const unsigned i = unsigned(queue->remove());

    MEDDLY_DCASSERT(nb.down(i));

    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
      // If event i needs rebuilding,
      //    Rebuild it, and
      //    Update the event reader
      if (rel->rebuildEvent(level, ei)) {
        const dd_edge& mxd = rel->getEvent(level, ei);
        if (0==mxd.getNode()) {
          if (Ru[ei]) {
            unpacked_node::Recycle(Ru[ei]);
            Ru[ei] = 0;
          }
        } else {
          if (0==Ru[ei]) {
            Ru[ei] = unpacked_node::New(relF, FULL_ONLY);
          }
          const int eventLevel = mxd.getLevel();
          if (ABS(eventLevel) < level || eventLevel < 0) {
            // Takes care of two situations:
            // - skipped unprimed level (due to Fully Reduced)
            // - skipped unprimed and primed levels (due to Fully Identity Reduced)
            Ru[ei]->initRedundant(level, mxd.getNode());
          } else {
            Ru[ei]->initFromNode(mxd.getNode());
          }
        }
      }
      // check if row i of the event ei is empty
      if (0 == Ru[ei]) continue;
      node_handle ei_i = (i < Ru[ei]->getSize())
                        ? Ru[ei]->down(i)
                        : 0;
      if (0 == ei_i) continue;

      // grab column (TBD: build these ahead of time?)
      const int dlevel = relF->getNodeLevel(ei_i);

      if (dlevel == -level) {
        Rp->initFromNode(ei_i);
      } else {
        Rp->initIdentity(-level, i, ei_i);
      }

      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        const unsigned j = Rp->index(jz);
        if (j < nb.getSize() && -1==nb.down(j)) continue;  // nothing can be added to this set

        node_handle newstates = recFire(nb.down(i), Rp->down(jz));
        if (newstates == 0) continue;

        // Confirm local state
        rel->confirm(level, int(j));

        if (j >= nb.getSize()) {
          const unsigned oldSize = nb.getSize();
          // resize the node builder, and clear out the new entries
          nb.resize(j+1);
          nb.clear(oldSize, nb.getSize());
          // while(oldSize < nb.getSize()) { nb.d_ref(oldSize++) = 0; }
          // resize the queue, and clear out the new entries
          queue->resize(nb.getSize());
        }

        bool updated = true;

        if (0 == nb.down(j)) {
          nb.setFull(j, newstates);
          // nb.d_ref(j) = newstates;
        } else {
          nbdj.set(nb.down(j));      // clobber
          newst.set(newstates);   // clobber
          mddUnion->computeTemp(nbdj, newst, nbdj);
          updated = (nbdj.getNode() != nb.down(j));
          nb.setFull(j, nbdj);
        }

        if (updated) queue->add((int)j);
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
  if (relF->isTerminalNode(mxd)) {
    if (argF->isTerminalNode(mdd)) {
      return resF->handleForValue(1);
    }
    // mxd is identity
    if (argF == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  node_handle result = 0;
  ct_entry_key* Key = findResult(mdd, mxd, result);
  if (0==Key) return result;

#ifdef TRACE_RECFIRE
  printf("computing recFire(%d, %d)\n", mdd, mxd);
  printf("  node %3d ", mdd);
  argF->showNode(stdout, mdd, 1);
  printf("\n  node %3d ", mxd);
  relF->showNode(stdout, mxd, 1);
  printf("\n");
#endif

  mxd = relF->linkNode(mxd);

  // check if mxd and mdd are at the same level
  const int mddLevel = argF->getNodeLevel(mdd);
  const int mxdLevel = relF->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* nb = unpacked_node::newWritable(resF, rLevel, rSize, FULL_ONLY);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(argF, FULL_ONLY);
  if (mddLevel < rLevel) {
    A->initRedundant(rLevel, mdd);
  } else {
    A->initFromNode(mdd);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (unsigned i=0; i<rSize; i++) {
      nb->setFull(i, recFire(A->down(i), mxd));
      // nb->d_ref(i) = recFire(A->down(i), mxd);
    }

  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New(relF, SPARSE_ONLY);
    unpacked_node *Rp = unpacked_node::New(relF, SPARSE_ONLY);
    if (mxdLevel < 0) {
      Ru->initRedundant(rLevel, mxd);
    } else {
      Ru->initFromNode(mxd);
    }

#if 0
    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getSize(); iz++) {
      const int i = Ru->index(iz);
      if (0==A->down(i)) continue;
      const node_handle pnode = Ru->down(iz);
      if (isLevelAbove(-rLevel, relF->getNodeLevel(pnode))) {
        Rp->initIdentity(relF, rLevel, i, pnode, SPARSE_ONLY);
      } else {
        Rp->initFromNode(pnode);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getSize(); jz++) {
        const int j = Rp->index(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = recFire(A->down(i), Rp->down(jz));
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

        if (0==nb->down(j)) {
          nb->d_ref(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        const int oldj = nb->down(j);
        nb->d_ref(j) = mddUnion->computeTemp(newstates, oldj);
        resF->unlinkNode(oldj);
        resF->unlinkNode(newstates);
      } // for j
    } // for i
#else
    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      const unsigned i = Ru->index(iz);
      if (0==A->down(i)) continue;
      recFireHelper(i, rLevel, Ru->down(iz), A->down(i), Rp, nb);
    }
    // loop over the extensible portion of mxd (if any)
#endif

    unpacked_node::Recycle(Rp);
    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);

  saturateHelper(*nb);
  edge_value ev;
  resF->createReducedNode(nb, ev, result);
  MEDDLY_DCASSERT(ev.isVoid());
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
  const unsigned i,
  const int rLevel,
  const node_handle Ru_i,
  const node_handle A_i,
  unpacked_node *Rp,
  unpacked_node* nb)
{
  if (isLevelAbove(-rLevel, relF->getNodeLevel(Ru_i))) {
    Rp->initIdentity(rLevel, i, Ru_i);
  } else {
    Rp->initFromNode(Ru_i);
  }

  dd_edge nbdj(resF), newst(resF);

  // loop over mxd "columns"
  for (unsigned jz=0; jz<Rp->getSize(); jz++) {
    const unsigned j = Rp->index(jz);
    // ok, there is an i->j "edge".
    // determine new states to be added (recursively)
    // and add them
    node_handle newstates = recFire(A_i, Rp->down(jz));
    if (0==newstates) continue;

    // Confirm local state
    rel->confirm(rLevel, int(j));

    if (j >= nb->getSize()) {
      const unsigned oldSize = nb->getSize();
      // resize the node builder, and clear out the new entries
      nb->resize(j+1);
      nb->clear(oldSize, nb->getSize());
      // while(oldSize < nb->getSize()) { nb->d_ref(oldSize++) = 0; }
    }

    if (0==nb->down(j)) {
      nb->setFull(j, newstates);
      // nb->d_ref(j) = newstates;
    } else {
      // there's new states and existing states; union them.
      nbdj.set(nb->down(j));
      newst.set(newstates);
      mddUnion->computeTemp(nbdj, newst, nbdj);
      nb->setFull(j, nbdj);
    }
  } // for j
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_operation*
MEDDLY::SATURATION_OTF_FORWARD(forest* inF, otf_relation* nsf, forest* outF)
{
    if (!nsf) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    if (nsf->getOutForest() != outF) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    return new forwd_otf_dfs_by_events_mt("OtfSaturationFwd", inF, nsf, outF);
}

MEDDLY::saturation_operation*
MEDDLY::SATURATION_OTF_BACKWARD(forest* inF, otf_relation* nsf, forest* outf)
{
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

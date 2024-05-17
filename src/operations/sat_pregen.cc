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

#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../oper_satur.h"
#include "../ops_builtin.h"
#include "../sat_relations.h"

#define DEBUG_FINALIZE
// #define DEBUG_FINALIZE_SPLIT
// #define DEBUG_EVENT_MASK

namespace MEDDLY {
  class saturation_by_events_op;

  class common_dfs_by_events_mt;
  class forwd_dfs_by_events_mt;
  class bckwd_dfs_by_events_mt;
} // Namespace MEDDLY



// ******************************************************************
// *                                                                *
// *                 saturation_by_events_op  class                 *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_by_events_op : public saturation_operation {
    common_dfs_by_events_mt* parent;
  public:
    saturation_by_events_op(common_dfs_by_events_mt* p,
        forest* argF, forest* resF);
    virtual ~saturation_by_events_op();

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
// *            common_dfs_by_events_mt  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::common_dfs_by_events_mt : public saturation_operation {
  public:
    common_dfs_by_events_mt(const char* name, forest* inf,
        pregen_relation* rel, forest* outf);
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

    pregen_relation* rel;
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
    forest* argF, forest* resF)
  : saturation_operation("Pregen-sat", 1, argF, resF)
{
  parent = p;

  ct_entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new ct_entry_type("Pregen-sat", "NI:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(3, resF);
  } else {
    et = new ct_entry_type("Pregen-sat", "N:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(2, resF);
  }
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::saturation_by_events_op::~saturation_by_events_op()
{
  // removeAllComputeTableEntries();
}

void MEDDLY::saturation_by_events_op::compute(const dd_edge& in, dd_edge& out)
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
    if (mddDptrs->down(i)) {
      nb->setFull(i, saturate(mddDptrs->down(i), k-1));
    }
    // nb->d_ref(i) = mddDptrs->down(i) ? saturate(mddDptrs->down(i), k-1) : 0;
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
  const char* name, forest* inF, pregen_relation* relation, forest* outF)
: saturation_operation(name, 1, inF, outF)
{
    if (!inF || !relation || !outF) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }

    rel = relation;
    relF = relation->getRelForest();

    // for now, anyway, inset and outset must be same forest
    if (inF != outF) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }

    // Check for same domain
    if  (
            (inF->getDomain() != relF->getDomain()) ||
            (outF->getDomain() != relF->getDomain())
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

MEDDLY::common_dfs_by_events_mt::~common_dfs_by_events_mt()
{
    unregisterInForest(relF);
    delete rel;
}

void MEDDLY::common_dfs_by_events_mt
::compute(const dd_edge &a, dd_edge &c)
{
  // Initialize operations
  mddUnion = UNION(resF, resF, resF);
  MEDDLY_DCASSERT(mddUnion);

  mxdIntersection = INTERSECTION(relF, relF, relF);
  MEDDLY_DCASSERT(mxdIntersection);

  mxdDifference = DIFFERENCE(relF, relF, relF);
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
  saturation_by_events_op* so = new saturation_by_events_op(this, argF, resF);
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
    forwd_dfs_by_events_mt(forest* inf, pregen_relation* relation,
            forest* outf);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::forwd_dfs_by_events_mt::forwd_dfs_by_events_mt(forest* inf,
    pregen_relation* rel, forest* outf)
  : common_dfs_by_events_mt("SaturationFwd", inf, rel, outf)
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
    Ru[ei] = unpacked_node::New(relF);
    int eventLevel = events[ei].getLevel();
    MEDDLY_DCASSERT(ABS(eventLevel) == nb.getLevel());
    if (eventLevel<0) {
      Ru[ei]->initRedundant(relF, nb.getLevel(), events[ei].getNode(), FULL_ONLY);
    } else {
      relF->unpackNode(Ru[ei], events[ei].getNode(), FULL_ONLY);
    }
  }
  unpacked_node* Rp = unpacked_node::New(relF);

  dd_edge nbdj(resF), newst(resF);

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (unsigned i = 0; i < nb.getSize(); i++) {
    if (nb.down(i)) queue->add(i);
  }

  // explore indexes
  while (!queue->isEmpty()) {
    unsigned i = queue->remove();

    MEDDLY_DCASSERT(nb.down(i));

    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
      if (0 == Ru[ei]->down(i)) continue;  // row i of the event ei is empty

      // grab column (TBD: build these ahead of time?)
      int dlevel = relF->getNodeLevel(Ru[ei]->down(i));

      if (dlevel == -nb.getLevel()) {
        relF->unpackNode(Rp, Ru[ei]->down(i), SPARSE_ONLY);
      } else {
        Rp->initIdentity(relF, -nb.getLevel(), i, Ru[ei]->down(i), SPARSE_ONLY);
      }

      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        unsigned j = Rp->index(jz);
        if (-1==nb.down(j)) continue;  // nothing can be added to this set

        node_handle rec = recFire(nb.down(i), Rp->down(jz));

        if (rec == 0) continue;
        if (rec == nb.down(j)) {
          resF->unlinkNode(rec);
          continue;
        }

        bool updated = true;

        if (0 == nb.down(j)) {
          nb.setFull(j, rec);
          // nb.d_ref(j) = rec;
        }
        else if (rec == -1) {
          resF->unlinkNode(nb.down(j));
          nb.setFull(j, rec);
          // nb.d_ref(j) = -1;
        }
        else {
          nbdj.set(nb.down(j));  // clobber
          newst.set(rec);     // clobber
          mddUnion->computeTemp(nbdj, newst, nbdj);
          updated = (nbdj.getNode() != nb.down(j));
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
  arg1F->showNode(stdout, mdd, 1);
  printf("\n  node %3d ", mxd);
  relF->showNode(stdout, mxd, 1);
  printf("\n");
#endif

  // check if mxd and mdd are at the same level
  const int mddLevel = argF->getNodeLevel(mdd);
  const int mxdLevel = relF->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);

  dd_edge nbdj(resF), newst(resF);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(argF);
  if (mddLevel < rLevel) {
    A->initRedundant(argF, rLevel, mdd, FULL_ONLY);
  } else {
    argF->unpackNode(A, mdd, FULL_ONLY);
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
    unpacked_node *Ru = unpacked_node::New(relF);
    unpacked_node *Rp = unpacked_node::New(relF);
    if (mxdLevel < 0) {
      Ru->initRedundant(relF, rLevel, mxd, SPARSE_ONLY);
    } else {
      relF->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      unsigned i = Ru->index(iz);
      if (0==A->down(i))   continue;
      if (isLevelAbove(-rLevel, relF->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(relF, rLevel, i, Ru->down(iz), SPARSE_ONLY);
      } else {
        relF->unpackNode(Rp, Ru->down(iz), SPARSE_ONLY);
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
          nb->setFull(j, newstates);
          // nb->d_ref(j) = newstates;
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
    bckwd_dfs_by_events_mt(forest* inf, pregen_relation* relation,
            forest* outf);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::bckwd_dfs_by_events_mt::bckwd_dfs_by_events_mt(forest* inf,
    pregen_relation* rel, forest* outf)
  : common_dfs_by_events_mt("SaturationBack", inf, rel, outf)
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
    Ru[ei] = unpacked_node::New(relF);
    int eventLevel = events[ei].getLevel();
    MEDDLY_DCASSERT(ABS(eventLevel) == nb.getLevel());
    if (eventLevel<0) {
      Ru[ei]->initRedundant(relF, nb.getLevel(), events[ei].getNode(), FULL_ONLY);
    } else {
      relF->unpackNode(Ru[ei], events[ei].getNode(), FULL_ONLY);
    }
  }
  unpacked_node* Rp = unpacked_node::New(relF);

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
        int dlevel = relF->getNodeLevel(Ru[ei]->down(iz));

        if (dlevel == -nb.getLevel()) {
          relF->unpackNode(Rp, Ru[ei]->down(iz), SPARSE_ONLY);
        } else {
          Rp->initIdentity(relF, -nb.getLevel(), i, Ru[ei]->down(iz), SPARSE_ONLY);
        }

        for (unsigned jz=0; jz<Rp->getSize(); jz++) {
          unsigned j = Rp->index(jz);
          if (0==expl->data[j]) continue;
          if (0==nb.down(j))       continue;
          // We have an i->j edge to explore
          node_handle rec = recFire(nb.down(j), Rp->down(jz));

          if (0==rec) continue;
          if (rec == nb.down(i)) {
            resF->unlinkNode(rec);
            continue;
          }

          bool updated = true;

          if (0 == nb.down(i)) {
            nb.setFull(i, rec);
            // nb.d_ref(i) = rec;
          }
          else if (-1 == rec) {
            resF->unlinkNode(nb.down(i));
            nb.setFull(i, rec);
            // nb.d_ref(i) = -1;
          }
          else {
            nbdi.set(nb.down(i));
            newst.set(rec);
            mddUnion->computeTemp(nbdi, newst, nbdi);
            updated = (nbdi.getNode() != nb.down(i));
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

  // check if mxd and mdd are at the same level
  const int mddLevel = argF->getNodeLevel(mdd);
  const int mxdLevel = relF->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);

  dd_edge nbdi(resF), newst(resF);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(argF);
  if (mddLevel < rLevel) {
    A->initRedundant(argF, rLevel, mdd, FULL_ONLY);
  } else {
    argF->unpackNode(A, mdd, FULL_ONLY);
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
    unpacked_node *Ru = unpacked_node::New(relF);
    unpacked_node *Rp = unpacked_node::New(relF);
    if (mxdLevel < 0) {
      Ru->initRedundant(relF, rLevel, mxd, SPARSE_ONLY);
    } else {
      relF->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      unsigned i = Ru->index(iz);
      if (isLevelAbove(-rLevel, relF->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(relF, rLevel, i, Ru->down(iz), SPARSE_ONLY);
      } else {
        relF->unpackNode(Rp, Ru->down(iz), SPARSE_ONLY);
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
          nb->setFull(i, newstates);
          // nb->d_ref(i) = newstates;
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
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_operation* MEDDLY::SATURATION_FORWARD(forest* inF,
        pregen_relation* rel, forest* outF)
{
    return new forwd_dfs_by_events_mt(inF, rel, outF);
}

MEDDLY::saturation_operation* MEDDLY::SATURATION_BACKWARD(forest* inF,
        pregen_relation* rel, forest* outF)
{
    return new bckwd_dfs_by_events_mt(inF, rel, outF);
}

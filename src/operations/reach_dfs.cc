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
#include "reach_dfs.h"

#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"

// #define TRACE_RECFIRE
// #define DEBUG_DFS
// #define DEBUG_INITIAL
// #define DEBUG_NSF
// #define DEBUG_SPLIT

namespace MEDDLY {
    class saturation_op;
    class saturation_evplus_op;

    class common_dfs;

    class common_dfs_mt;

    class forwd_dfs_mt;
    class bckwd_dfs_mt;

    class common_dfs_evplus;

    class forwd_dfs_evplus;

    binary_list FWD_DFS_cache;
    binary_list REV_DFS_cache;
};


// ******************************************************************
// *                                                                *
// *                      saturation_op  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_op : public unary_operation {
    common_dfs_mt* parent;
  public:
    saturation_op(common_dfs_mt* p, unary_list &c, forest* argF, forest* resF);
    virtual ~saturation_op();

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
// *                   saturation_evplus_op class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_evplus_op : public unary_operation {
    common_dfs_evplus* parent;
  public:
    saturation_evplus_op(common_dfs_evplus* p, unary_list &c, forest* argF, forest* resF);
    virtual ~saturation_evplus_op();

    void saturate(const dd_edge& in, dd_edge& out);
    // void saturate(long ev, node_handle evmdd, long& resEv, node_handle& resEvmdd);
    void saturate(long ev, node_handle evmdd, int level, long& resEv, node_handle& resEvmdd);

  protected:
    inline ct_entry_key*
    findSaturateResult(long aev, node_handle a, int level, long& bev, node_handle& b) {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      if (argF->isFullyReduced()) CTsrch->writeI(level);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      bev = CTresult[0].readL();
      bev += aev;
      b = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveSaturateResult(ct_entry_key* Key,
      long aev, node_handle a, long bev, node_handle b)
    {
      CTresult[0].reset();
      CTresult[0].writeL(bev - aev);
      CTresult[0].writeN(b);
      CT0->addEntry(Key, CTresult[0]);
      return b;
    }
};

// ******************************************************************
// *                                                                *
// *                        common_dfs class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::common_dfs : public binary_operation {
  public:
    common_dfs(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res);

  protected:
    dd_edge* splits;
    binary_operation* mddUnion;
    binary_operation* mxdIntersection;
    binary_operation* mxdDifference;

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
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, unsigned(head), size);
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

    void splitMxd(node_handle mxd);

    void cleanup();

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
};

// ******************************************************************
// *                                                                *
// *                      common_dfs_mt  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::common_dfs_mt : public common_dfs {
  public:
    common_dfs_mt(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res);

    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);
    virtual void saturateHelper(unpacked_node &mdd) = 0;

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

};

// ******************************************************************
// *                                                                *
// *                    common_dfs_evplus  class                    *
// *                                                                *
// ******************************************************************

class MEDDLY::common_dfs_evplus : public common_dfs {
  public:
    common_dfs_evplus(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res);

    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);
    virtual void saturateHelper(unpacked_node &mdd) = 0;

  protected:
    inline ct_entry_key*
    findResult(long aev, node_handle a, node_handle b, long& cev, node_handle& c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CTsrch->writeN(b);

      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      cev = CTresult[0].readL();
      c = resF->linkNode(CTresult[0].readN());
      if (c != 0) {
        cev += aev;
      }
      CT0->recycle(CTsrch);
      return 0;
    }
    inline void saveResult(ct_entry_key* Key,
      long aev, node_handle a, node_handle b, long cev, node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeL(c == 0 ? 0L : cev - aev);
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
    }
};

// ******************************************************************
// *                                                                *
// *                     saturation_op  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_op
::saturation_op(common_dfs_mt* p, unary_list &c, forest* argF, forest* resF)
  : unary_operation(c, 1, argF, resF)
{
  parent = p;

  ct_entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new ct_entry_type(c.getName(), "NI:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(3, resF);
  } else {
    et = new ct_entry_type(c.getName(), "N:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(2, resF);
  }
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::saturation_op::~saturation_op()
{
}

void MEDDLY::saturation_op::saturate(const dd_edge& in, dd_edge& out)
{
  out.set( saturate(in.getNode(), argF->getMaxLevelIndex()) );
}

MEDDLY::node_handle MEDDLY::saturation_op::saturate(node_handle mdd, int k)
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

  unpacked_node* C = unpacked_node::newWritable(resF, k, sz, FULL_ONLY);
  // Initialize mdd reader
  unpacked_node *mddDptrs = unpacked_node::New(argF, FULL_ONLY);
  if (mdd_level < k) {
    mddDptrs->initRedundant(k, mdd);
  } else {
    mddDptrs->initFromNode(mdd);
  }

  // Do computation
  for (unsigned i=0; i<sz; i++) {
    C->setFull(i,
      mddDptrs->down(i) ? saturate(mddDptrs->down(i), k-1) : 0
    );
  }

  // Cleanup
  unpacked_node::Recycle(mddDptrs);

  parent->saturateHelper(*C);
  edge_value ev;
  resF->createReducedNode(C, ev, n);
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
// *                 saturation_evplus_op  methods                  *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_evplus_op
::saturation_evplus_op(common_dfs_evplus* p, unary_list &c, forest* argF, forest* resF)
  : unary_operation(c, 1, argF, resF)
{
  parent = p;

  ct_entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new ct_entry_type(c.getName(), "NI:LN");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(4, resF);
  } else {
    et = new ct_entry_type(c.getName(), "N:LN");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(3, resF);
  }
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::saturation_evplus_op::~saturation_evplus_op()
{
}

void MEDDLY::saturation_evplus_op::saturate(const dd_edge& in, dd_edge& out)
{
  long aev = Inf<long>();
  in.getEdgeValue(aev);
  long cev = Inf<long>();
  node_handle cnode = 0;
  saturate(aev, in.getNode(), argF->getMaxLevelIndex(), cev, cnode);
  out.set(cev, cnode);
}

void MEDDLY::saturation_evplus_op::saturate(long ev, node_handle evmdd, int k, long& resEv, node_handle& resEvmdd)
{
#ifdef DEBUG_DFS
  printf("evmdd: %d, ev: %ld, k: %d\n", evmdd, ev, k);
#endif

  // terminal condition for recursion
  if (argF->isTerminalNode(evmdd)) {
    resEv = ev;
    resEvmdd = evmdd;
    return;
  }

  // search compute table
  ct_entry_key* Key = findSaturateResult(ev, evmdd, k, resEv, resEvmdd);
  if (0==Key) {
    return;
  }

  const unsigned sz = unsigned(argF->getLevelSize(k));    // size
  const int evmdd_level = argF->getNodeLevel(evmdd);      // evmdd level

#ifdef DEBUG_DFS
  printf("evmdd: %d, level: %d, size: %d, evmdd_level: %d\n",
      evmdd, k, sz, evmdd_level);
#endif

  unpacked_node* C = unpacked_node::newWritable(resF, k, sz, FULL_ONLY);
  // Initialize evmdd reader
  unpacked_node *evmddDptrs = unpacked_node::New(argF, FULL_ONLY);
  if (evmdd_level < k) {
    evmddDptrs->initRedundant(k, evmdd);
  } else {
    evmddDptrs->initFromNode(evmdd);
  }

  // Do computation
  for (unsigned i=0; i<sz; i++) {
    long cev = Inf<long>();
    node_handle c = 0;
    if (evmddDptrs->down(i) != 0) {
      saturate(evmddDptrs->edgeval(i).getLong(), evmddDptrs->down(i), k-1, cev, c);
    }
    C->setFull(i, edge_value(cev), c);
    // C->setEdge(i, cev);
    // C->d_ref(i) = c;
  }

  // Cleanup
  unpacked_node::Recycle(evmddDptrs);

  parent->saturateHelper(*C);
  edge_value Cev;
  resF->createReducedNode(C, Cev, resEvmdd);
  resEv = Cev.getLong() + ev;

  // save in compute table
  saveSaturateResult(Key, ev, evmdd, resEv, resEvmdd);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, resEvmdd);
#endif
}


// ******************************************************************
// *                                                                *
// *                       common_dfs methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::common_dfs::common_dfs(binary_list& oc, forest* a1,
  forest* a2, forest* res)
: binary_operation(oc, 1, a1, a2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkRelations(__FILE__, __LINE__, SET, RELATION, SET);
    if (a1 != res) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    splits = 0;
    mddUnion = 0;
    mxdIntersection = 0;
    mxdDifference = 0;
    freeqs = 0;
    freebufs = 0;
}


// Partition the nsf based on "top level"
void MEDDLY::common_dfs::splitMxd(node_handle mxd_nh)
{
  MEDDLY_DCASSERT(arg2F);
  MEDDLY_DCASSERT(0==splits);
  splits = new dd_edge[arg2F->getNumVariables()+1];

  dd_edge maxDiag(arg2F);
  dd_edge mpdi(arg2F);
  dd_edge mxd(arg2F);
  mxd.set( arg2F->linkNode(mxd_nh) );

  // Build from top down
  for (int level = arg2F->getMaxLevelIndex(); level; level--) {
    splits[level].attach(arg2F);

    if (0==mxd.getNode()) {
      // common and easy special case
      continue;
    }

    int mxdLevel = mxd.getLevel();
    MEDDLY_DCASSERT(ABS(mxdLevel) <= level);

    // Initialize readers
    unpacked_node *Mu = unpacked_node::New(arg2F, FULL_ONLY);
    unpacked_node *Mp = unpacked_node::New(arg2F, FULL_ONLY);
    if (isLevelAbove(level, mxdLevel)) {
      Mu->initRedundant(level, mxd.getNode());
    } else {
      Mu->initFromNode(mxd.getNode());
    }

    // Read "rows"
    for (unsigned i=0; i<Mu->getSize(); i++) {
      // Initialize column reader
      int mxdPLevel = arg2F->getNodeLevel(Mu->down(i));
      if (isLevelAbove(-level, mxdPLevel)) {
        Mp->initIdentity(-level, i, Mu->down(i));
      } else {
        Mp->initFromNode(Mu->down(i));
      }

      // Intersect along the diagonal
      if (0==i) {
        maxDiag.set( arg2F->linkNode(Mp->down(i)) );
      } else {
        mpdi.set( arg2F->linkNode(Mp->down(i)) );
        mxdIntersection->computeTemp(maxDiag, mpdi, maxDiag);
      }

      // cleanup
    } // for i

    // maxDiag is what we can split from here
    mxdDifference->computeTemp(mxd, maxDiag, splits[level]);
    mxd = maxDiag;

    // Cleanup
    unpacked_node::Recycle(Mp);
    unpacked_node::Recycle(Mu);
  } // for level

#ifdef DEBUG_SPLIT
  printf("After splitting monolithic event in msat\n");
  printf("splits array: [");
  for (unsigned k=0; k <= arg2F->getNumVariables(); k++) {
    if (k) printf(", ");
    printf("%d", splits[k]);
  }
  printf("]\n");
#endif
}


void MEDDLY::common_dfs::cleanup()
{
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
  // for (int i = arg2F->getNumVariables(); i; i--) arg2F->unlinkNode(splits[i]);
  // Shouldn't need to do the above b/c we're using dd_edges now
  delete[] splits;
  splits = 0;
}

// ******************************************************************
// *                   common_dfs::indexq methods                   *
// ******************************************************************

MEDDLY::common_dfs::indexq::indexq()
{
  data = 0;
  size = 0;
  head = NULPTR;
}

MEDDLY::common_dfs::indexq::~indexq()
{
  free(data);
}

void MEDDLY::common_dfs::indexq::resize(unsigned sz)
{
  if (sz <= size) return;
  data = (int*) realloc(data, sz * sizeof(int));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);

  for (; size < sz; size++) data[size] = NOTINQ;
}

// ******************************************************************
// *                  common_dfs::charbuf  methods                  *
// ******************************************************************

MEDDLY::common_dfs::charbuf::charbuf()
{
  data = 0;
  size = 0;
}

MEDDLY::common_dfs::charbuf::~charbuf()
{
  free(data);
}

void MEDDLY::common_dfs::charbuf::resize(unsigned sz)
{
  if (sz <= size) return;
  data = (char*) realloc(data, sz * sizeof(char));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *                     common_dfs_mt  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::common_dfs_mt::common_dfs_mt(binary_list& oc, forest* a1,
  forest* a2, forest* res)
: common_dfs(oc, a1, a2, res)
{
  ct_entry_type* et = new ct_entry_type(oc.getName(), "NN:N");
  et->setForestForSlot(0, a1);
  et->setForestForSlot(1, a2);
  et->setForestForSlot(3, res);
  registerEntryType(0, et);
  buildCTs();
}

void MEDDLY::common_dfs_mt
::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag)
{
  // Initialize operations
  mddUnion = UNION(resF, resF, resF);
  MEDDLY_DCASSERT(mddUnion);

  mxdIntersection = INTERSECTION(arg2F, arg2F, arg2F);
  MEDDLY_DCASSERT(mxdIntersection);

  mxdDifference = DIFFERENCE(arg2F, arg2F, arg2F);
  MEDDLY_DCASSERT(mxdDifference);

#ifdef DEBUG_INITIAL
  printf("Calling saturate for states:\n");
  a.show(stdout, 2);
#endif
#ifdef DEBUG_NSF
  printf("Calling saturate for NSF:\n");
  b.show(stdout, 2);
#endif

  // Partition NSF by levels
  splitMxd(b.getNode());

  // Execute saturation operation
  unary_list dummy("Saturation");
  saturation_op *so = new saturation_op(this, dummy, arg1F, resF);
  so->saturate(a, c);

  // Cleanup
  cleanup();
  // so->removeAllComputeTableEntries();
  delete so;
}


// ******************************************************************
// *                                                                *
// *                       forwd_dfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_dfs_mt : public common_dfs_mt {
  public:
    forwd_dfs_mt(forest* arg1, forest* arg2, forest* res);
  protected:
    virtual void saturateHelper(unpacked_node &mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::forwd_dfs_mt::forwd_dfs_mt(forest* arg1, forest* arg2, forest* res)
  : common_dfs_mt(FWD_DFS_cache, arg1, arg2, res)
{
}

void MEDDLY::forwd_dfs_mt::saturateHelper(unpacked_node &nb)
{
  node_handle mxd = splits[nb.getLevel()].getNode();
  if (mxd == 0) return;

  const int mxdLevel = arg2F->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node *Ru = unpacked_node::New(arg2F, FULL_ONLY);
  unpacked_node *Rp = unpacked_node::New(arg2F, SPARSE_ONLY);
  if (mxdLevel < 0) {
    Ru->initRedundant(nb.getLevel(), mxd);
  } else {
    Ru->initFromNode(mxd);
  }

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (unsigned i = 0; i < nb.getSize(); i++) {
    if (nb.down(i)) queue->add(i);
  }

  dd_edge nbdj(resF), temp(resF);

  // explore indexes
  while (!queue->isEmpty()) {
    const unsigned i = queue->remove();

    MEDDLY_DCASSERT(nb.down(i));
    if (0==Ru->down(i)) continue;  // row i is empty

    // grab column (TBD: build these ahead of time?)
    const int dlevel = arg2F->getNodeLevel(Ru->down(i));

    if (dlevel == -nb.getLevel()) {
      Rp->initFromNode(Ru->down(i));
    } else {
      Rp->initIdentity(-nb.getLevel(), i, Ru->down(i));
    }

    for (int jz=0; jz<Rp->getSize(); jz++) {
      MEDDLY_DCASSERT(jz>=0);
      const unsigned j = Rp->index(unsigned(jz));
      if (-1==nb.down(j)) continue;  // nothing can be added to this set

      node_handle rec = recFire(nb.down(i), Rp->down(unsigned(jz)));

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
      }
      else {
        temp.set(rec);  // clobbers rec; that's what we want
        nbdj.set( nb.down(j) );
        mddUnion->computeTemp(nbdj, temp, nbdj);
        updated = (nbdj.getNode() != nb.down(j));
        nb.setFull(j, nbdj);
      }

      if (updated) {
        if (j == i) {
          // Restart inner for-loop.
          jz = -1;
        } else {
          queue->add(j);
        }
      }

    } // for j

  } // while there are indexes to explore

  // cleanup
  unpacked_node::Recycle(Rp);
  unpacked_node::Recycle(Ru);
  recycle(queue);
}

// Same as post-image, except we saturate before reducing.
MEDDLY::node_handle MEDDLY::forwd_dfs_mt::recFire(node_handle mdd, node_handle mxd)
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
  unpacked_node* nb = unpacked_node::newWritable(resF, rLevel, rSize, FULL_ONLY);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(arg1F, FULL_ONLY);
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
    unpacked_node *Ru = unpacked_node::New(arg2F, SPARSE_ONLY);
    unpacked_node *Rp = unpacked_node::New(arg2F, SPARSE_ONLY);
    if (mxdLevel < 0) {
      Ru->initRedundant(rLevel, mxd);
    } else {
      Ru->initFromNode(mxd);
    }

    dd_edge nbdj(resF), newst(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      const unsigned i = Ru->index(iz);
      if (0==A->down(i))   continue;
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(rLevel, i, Ru->down(iz));
      } else {
        Rp->initFromNode(Ru->down(iz));
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
        newst.set(newstates); // clobber when done
        nbdj.set(nb->down(j));   // also clobber when done
        mddUnion->computeTemp(newst, nbdj, nbdj);
        nb->setFull(j, nbdj);
      } // for j

    } // for i

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




// ******************************************************************
// *                                                                *
// *                       bckwd_dfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_dfs_mt : public common_dfs_mt {
  public:
    bckwd_dfs_mt(forest* arg1, forest* arg2, forest* res);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::bckwd_dfs_mt::bckwd_dfs_mt(forest* arg1, forest* arg2, forest* res)
  : common_dfs_mt(REV_DFS_cache, arg1, arg2, res)
{
}

void MEDDLY::bckwd_dfs_mt::saturateHelper(unpacked_node& nb)
{
  node_handle mxd = splits[nb.getLevel()].getNode();
  if (mxd == 0) return;

  const int mxdLevel = arg2F->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node *Ru = unpacked_node::New(arg2F, SPARSE_ONLY);
  unpacked_node *Rp = unpacked_node::New(arg2F, SPARSE_ONLY);
  if (mxdLevel < 0) {
    Ru->initRedundant(nb.getLevel(), mxd);
  } else {
    Ru->initFromNode(mxd);
  }

  // indexes to explore
  charbuf* expl = useCharBuf(nb.getSize());
  for (unsigned i = 0; i < nb.getSize(); i++) expl->data[i] = 2;
  bool repeat = true;

  dd_edge nbdi(resF), temp(resF);

  // explore
  while (repeat) {
    // "advance" the explore list
    for (unsigned i=0; i<nb.getSize(); i++) if (expl->data[i]) expl->data[i]--;
    repeat = false;

    // explore all rows
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      const unsigned i = Ru->index(iz);
      // grab column (TBD: build these ahead of time?)
      const int dlevel = arg2F->getNodeLevel(Ru->down(iz));

      if (dlevel == -nb.getLevel()) {
        Rp->initFromNode(Ru->down(iz));
      } else {
        Rp->initIdentity(-nb.getLevel(), i, Ru->down(iz));
      }

      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        const unsigned j = Rp->index(jz);
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
          temp.set(rec);
          mddUnion->computeTemp(nbdi, temp, nbdi);
          updated = nbdi.getNode() != nb.down(i);
          nb.setFull(i, nbdi);
        }
        if (updated) {
          expl->data[i] = 2;
          repeat = true;
        }
      } // for j
    } // for i
  } // while repeat
  // cleanup
  unpacked_node::Recycle(Rp);
  unpacked_node::Recycle(Ru);
  recycle(expl);
}

MEDDLY::node_handle MEDDLY::bckwd_dfs_mt::recFire(node_handle mdd, node_handle mxd)
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
  unpacked_node* nb = unpacked_node::newWritable(resF, rLevel, rSize, FULL_ONLY);

  dd_edge nbdi(resF), temp(resF);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(arg1F, FULL_ONLY);
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
    unpacked_node *Ru = unpacked_node::New(arg2F, SPARSE_ONLY);
    unpacked_node *Rp = unpacked_node::New(arg2F, SPARSE_ONLY);
    if (mxdLevel < 0) {
      Ru->initRedundant(rLevel, mxd);
    } else {
      Ru->initFromNode(mxd);
    }

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      const unsigned i = Ru->index(iz);
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(rLevel, i, Ru->down(iz));
      } else {
        Rp->initFromNode(Ru->down(iz));
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        const unsigned j = Rp->index(jz);
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
        temp.set(newstates);
        mddUnion->computeTemp(temp, nbdi, nbdi);
        nb->setFull(i, nbdi);
      } // for j

    } // for i

    unpacked_node::Recycle(Rp);
    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);

  saturateHelper(*nb);
  edge_value ev;;
  resF->createReducedNode(nb, ev, result);
  MEDDLY_DCASSERT(ev.isVoid());
#ifdef TRACE_ALL_OPS
  printf("computed recFire(%d, %d) = %d\n", mdd, mxd, result);
#endif
  return saveResult(Key, mdd, mxd, result);
}

// ******************************************************************
// *                                                                *
// *                     common_dfs_evplus class                    *
// *                                                                *
// ******************************************************************

MEDDLY::common_dfs_evplus::common_dfs_evplus(binary_list& oc, forest* a1,
  forest* a2, forest* res)
: common_dfs(oc, a1, a2, res)
{
  ct_entry_type* et = new ct_entry_type(oc.getName(), "NN:LN");
  et->setForestForSlot(0, a1);
  et->setForestForSlot(1, a2);
  et->setForestForSlot(4, res);
  registerEntryType(0, et);
  buildCTs();
}


void MEDDLY::common_dfs_evplus
::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag)
{
  // Initialize operations
  mddUnion = UNION(resF, resF, resF);
  MEDDLY_DCASSERT(mddUnion);

  mxdIntersection = INTERSECTION(arg2F, arg2F, arg2F);
  MEDDLY_DCASSERT(mxdIntersection);

  mxdDifference = DIFFERENCE(arg2F, arg2F, arg2F);
  MEDDLY_DCASSERT(mxdDifference);

#ifdef DEBUG_INITIAL
  printf("Calling saturate for states:\n");
  a.show(stdout, 2);
#endif
#ifdef DEBUG_NSF
  printf("Calling saturate for NSF:\n");
  b.show(stdout, 2);
#endif

  // Partition NSF by levels
  splitMxd(b.getNode());

  // Execute saturation operation
  unary_list dummy("Saturation");
  saturation_evplus_op *so = new saturation_evplus_op(this, dummy, arg1F, resF);
  so->saturate(a, c);

  // Cleanup
  cleanup();
  // so->removeAllComputeTableEntries();
  delete so;
}


// ******************************************************************
// *                                                                *
// *                     forwd_dfs_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_dfs_evplus : public common_dfs_evplus {
  public:
  forwd_dfs_evplus(forest* arg1, forest* arg2, forest* res);
  protected:
    virtual void saturateHelper(unpacked_node &mdd);
    void recFire(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
};

MEDDLY::forwd_dfs_evplus::forwd_dfs_evplus(forest* arg1, forest* arg2,
        forest* res) : common_dfs_evplus(FWD_DFS_cache, arg1, arg2, res)
{
}

void MEDDLY::forwd_dfs_evplus::saturateHelper(unpacked_node &nb)
{
  node_handle mxd = splits[nb.getLevel()].getNode();
  if (mxd == 0) return;

  const int mxdLevel = arg2F->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node *Ru = unpacked_node::New(arg2F, FULL_ONLY);
  unpacked_node *Rp = unpacked_node::New(arg2F, SPARSE_ONLY);
  if (mxdLevel < 0) {
    Ru->initRedundant(nb.getLevel(), mxd);
  } else {
    Ru->initFromNode(mxd);
  }

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (unsigned i = 0; i < nb.getSize(); i++) {
    if (nb.down(i) != 0) {
      queue->add(i);
    }
  }

  dd_edge nbdj(resF), temp(resF);

  // explore indexes
  while (!queue->isEmpty()) {
    const unsigned i = queue->remove();

    MEDDLY_DCASSERT(nb.down(i));
    if (0==Ru->down(i)) continue;  // row i is empty

    // grab column (TBD: build these ahead of time?)
    const int dlevel = arg2F->getNodeLevel(Ru->down(i));

    if (dlevel == -nb.getLevel()) {
      Rp->initFromNode(Ru->down(i));
    } else {
      Rp->initIdentity(-nb.getLevel(), i, Ru->down(i));
    }

    for (int jz=0; jz<Rp->getSize(); jz++) {
      MEDDLY_DCASSERT(jz >= 0);
      const unsigned j = Rp->index(unsigned(jz));

      long recev = Inf<long>();
      node_handle rec = 0;
      recFire(nb.edgeval(i).getLong(), nb.down(i), Rp->down(unsigned(jz)), recev, rec);

      if (rec == 0) continue;

      // Increase the distance
      recev++;

      if (rec == nb.down(j)) {
        if (recev < nb.edgeval(j).getLong()) {
          nb.edgeval(j) = recev;
        }
        resF->unlinkNode(rec);
        continue;
      }

      bool updated = true;

      if (0 == nb.down(j)) {
        nb.setFull(j, edge_value(recev), rec);
        // nb.setEdge(j, recev);
        // nb.d_ref(j) = rec;
      }
//      else if (rec == -1) {
//        resF->unlinkNode(nb.down(j));
//        nb.d_ref(j) = -1;
//      }
      else {
        nbdj.set(nb.edgeval(j).getLong(), nb.down(j));
        temp.set(recev, rec);  // clobbers rec; that's what we want
        mddUnion->computeTemp(nbdj, temp, nbdj);
        updated = (nbdj.getNode() != nb.down(j));
        nb.setFull(j, nbdj);
      }

      if (updated) {
        if (j == i) {
          // Restart inner for-loop.
          jz = -1;
        } else {
          queue->add(j);
        }
      }
    } // for j
  } // while there are indexes to explore

  // cleanup
  unpacked_node::Recycle(Rp);
  unpacked_node::Recycle(Ru);
  recycle(queue);
}

// Same as post-image, except we saturate before reducing.
void MEDDLY::forwd_dfs_evplus::recFire(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  // termination conditions
  if (mxd == 0 || evmdd == 0) {
    resEv = 0;
    resEvmdd = 0;
    return;
  }
  if (arg2F->isTerminalNode(mxd)) {
    // XXX: Can be merged???
    if (arg1F->isTerminalNode(evmdd)) {
      resEv = ev;
      resEvmdd = evmdd;
      return;
    }
    // mxd is identity
    if (arg1F == resF) {
      resEv = ev;
      resEvmdd = resF->linkNode(evmdd);
      return;
    }
  }

  // check the cache
  ct_entry_key* Key = findResult(ev, evmdd, mxd, resEv, resEvmdd);
  if (0==Key) return;

#ifdef TRACE_RECFIRE
  printf("computing recFire(<%d, %d>, %d)\n", ev, evmdd, mxd);
  printf("  node <%d, %d> ", ev, evmdd);
  arg1F->showNode(stdout, evmdd, 1);
  printf("\n  node %3d ", mxd);
  arg2F->showNode(stdout, mxd, 1);
  printf("\n");
#endif

  // check if mxd and evmdd are at the same level
  const int evmddLevel = arg1F->getNodeLevel(evmdd);
  const int mxdLevel = arg2F->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), evmddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* nb = unpacked_node::newWritable(resF, rLevel, rSize, FULL_ONLY);

  // Initialize evmdd reader
  unpacked_node *A = unpacked_node::New(arg1F, FULL_ONLY);
  if (evmddLevel < rLevel) {
    A->initRedundant(rLevel, evmdd);
  } else {
    A->initFromNode(evmdd);
  }

  if (evmddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (unsigned i=0; i<rSize; i++) {
      long nev = Inf<long>();
      node_handle n = 0;
      recFire(A->edgeval(i).getLong() + ev, A->down(i), mxd, nev, n);
      nb->setFull(i, edge_value(nev), n);
      // nb->setEdge(i, nev);
      // nb->d_ref(i) = n;
    }

  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= evmddLevel);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New(arg2F, SPARSE_ONLY);
    unpacked_node *Rp = unpacked_node::New(arg2F, SPARSE_ONLY);
    if (mxdLevel < 0) {
      Ru->initRedundant(rLevel, mxd);
    } else {
      Ru->initFromNode(mxd);
    }

    dd_edge nbdj(resF), newst(resF);

    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getSize(); iz++) {
      const unsigned i = Ru->index(iz);
      if (0==A->down(i))   continue;
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->down(iz)))) {
        Rp->initIdentity(rLevel, i, Ru->down(iz));
      } else {
        Rp->initFromNode(Ru->down(iz));
      }

      // loop over mxd "columns"
      for (unsigned jz=0; jz<Rp->getSize(); jz++) {
        unsigned j = Rp->index(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        long nev = Inf<long>();
        node_handle n = 0;
        recFire(A->edgeval(i).getLong() + ev, A->down(i), Rp->down(jz), nev, n);

        if (0==n) continue;
        if (0==nb->down(j)) {
          nb->setFull(j, edge_value(nev), n);
          // nb->setEdge(j, nev);
          // nb->d_ref(j) = n;
          continue;
        }

        // there's new states and existing states; union them.
        newst.set(nev, n); // clobber when done
        nbdj.set(nb->edgeval(j).getLong(), nb->down(j));   // also clobber when done
        mddUnion->computeTemp(newst, nbdj, nbdj);
        nb->setFull(j, nbdj);
      } // for j

    } // for i

    unpacked_node::Recycle(Rp);
    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);

  saturateHelper(*nb);
  edge_value rev;
  resF->createReducedNode(nb, rev, resEvmdd);
  resEv = rev.getLong();
#ifdef TRACE_ALL_OPS
  printf("computed recfire(<%d, %d>, %d) = <%d, %d>\n", ev, evmdd, mxd, resEv, resEvmdd);
#endif
#ifdef TRACE_RECFIRE
  printf("computed recfire(<%d, %d>, %d) = <%d, %d>\n", ev, evmdd, mxd, resEv, resEvmdd);
  printf("  node <%d, %d> ", resEv, resEvmdd);
  resF->showNode(stdout, resEvmdd, 1);
  printf("\n");
#endif

  saveResult(Key, ev, evmdd, mxd, resEv, resEvmdd);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::REACHABLE_STATES_DFS(forest* a, forest* b,
        forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  FWD_DFS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if ( (b->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL) ) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        return FWD_DFS_cache.add(new forwd_dfs_mt(a, b, c));
    }
    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {
        return FWD_DFS_cache.add(new forwd_dfs_evplus(a, b, c));
    }
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::REACHABLE_STATES_DFS_init()
{
    FWD_DFS_cache.reset("ReachableDFS");
}

void MEDDLY::REACHABLE_STATES_DFS_done()
{
    MEDDLY_DCASSERT(FWD_DFS_cache.isEmpty());
}

MEDDLY::binary_operation* MEDDLY::REVERSE_REACHABLE_DFS(forest* a, forest* b,
        forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  REV_DFS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if ( (b->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL) ) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        return REV_DFS_cache.add(new bckwd_dfs_mt(a, b, c));
    }
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::REVERSE_REACHABLE_DFS_init()
{
    REV_DFS_cache.reset("ReverseReachableDFS");
}

void MEDDLY::REVERSE_REACHABLE_DFS_done()
{
    MEDDLY_DCASSERT(REV_DFS_cache.isEmpty());
}


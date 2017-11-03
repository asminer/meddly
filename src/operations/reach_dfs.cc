
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


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "reach_dfs.h"
#include <deque>
#include <algorithm>

// #define TRACE_RECFIRE
// #define DEBUG_DFS
// #define DEBUG_INITIAL
// #define DEBUG_NSF
// #define DEBUG_SPLIT

//#define COUNT_UNION_LENGTH
//#define BATCHED_UNLINK        // unlink intermediate results as a batch
#define BATCHED_UNIONS
//#define USE_NODE_COUNT

#ifdef BATCHED_UNLINK
std::deque<MEDDLY::node_handle> to_be_unlinked;
#endif

#ifdef BATCHED_UNIONS

#define BATCHED_UNIONS_FOR_SATURATE   // enables batch unions in saturate()
#define BATCHED_UNIONS_FOR_REC_FIRE   // enables batch unions in rec_fire()

#define DQ  // acc starts with 0, and d(i) is placed in queue
            // otherwise, acc starts with d(i)

#ifdef USE_NODE_COUNT

// Node count heuristic for batch unions
#define COMPARISON_OPERATION NODE_COUNT
#define RESULT_TYPE long
#define MEDDLY_RESULT_TYPE INTEGER
#define DESCENDING_ORDER false

#else

// Cardinality heuristic for batch unions
#define COMPARISON_OPERATION CARDINALITY
#define RESULT_TYPE double
#define MEDDLY_RESULT_TYPE REAL
#define DESCENDING_ORDER true

#endif  // USE_NODE_COUNT

#endif  // BATCHED_UNIONS

namespace MEDDLY {
  class saturation_opname;
  class saturation_op;
  class saturation_evplus_op;

  class common_dfs_mt;

  class forwd_dfs_mt;
  class bckwd_dfs_mt;

  class common_dfs_evplus;

  class forwd_dfs_evplus;

  class forwd_dfs_opname;
  class bckwd_dfs_opname;
};

#ifdef BATCHED_UNIONS
class sorter {
  MEDDLY::unary_operation *op;
  int level;
  bool descending;

  public:

  sorter(MEDDLY::unary_operation *op)
    : op(op), level(0), descending(false) {}

  void setLevel(int k) { level = k; }
  void setDescending(bool d) { descending = d; }

  // returns true if a must be placed before b in the sorted order
  bool operator()(MEDDLY::node_handle a, MEDDLY::node_handle b) const {
    if (descending) {
      // descending order
      // returns true if a > b
      if (a == b || b == -1 || a == 0) return false;
      if (b == 0 || a == -1) return true;
      RESULT_TYPE card_a = 0, card_b = 0;
      op->compute(level, a, card_a);
      op->compute(level, b, card_b);
      return card_a > card_b;
    }
    // ascending order
    // returns true if a < b
    if (a == b || a == -1 || b == 0) return false;
    if (a == 0 || b == -1) return true;
    RESULT_TYPE card_a = 0, card_b = 0;
    op->compute(level, a, card_a);
    op->compute(level, b, card_b);
    return card_a < card_b;
  }
};
#endif

// ******************************************************************
// *                                                                *
// *                    saturation_opname  class                    *
// *                                                                *
// ******************************************************************

/** Simple class to keep compute table happy.
*/
class MEDDLY::saturation_opname : public unary_opname {
  static saturation_opname* instance;
  public:
    saturation_opname();

    static const saturation_opname* getInstance();
 
};

MEDDLY::saturation_opname* MEDDLY::saturation_opname::instance = 0;

MEDDLY::saturation_opname::saturation_opname()
 : unary_opname("Saturate")
{
}

const MEDDLY::saturation_opname* MEDDLY::saturation_opname::getInstance()
{
  if (0==instance) instance = new saturation_opname;
  return instance;
}

// ******************************************************************
// *                                                                *
// *                      saturation_op  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_op : public unary_operation {
    common_dfs_mt* parent;
  public:
    saturation_op(common_dfs_mt* p, expert_forest* argF, expert_forest* resF);
    virtual ~saturation_op();

    node_handle saturate(node_handle mdd);
    node_handle saturate(node_handle mdd, int level);

#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
#endif
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
// *                      saturation_evplus_op  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_evplus_op : public unary_operation {
    common_dfs_evplus* parent;
  public:
    saturation_evplus_op(common_dfs_evplus* p, expert_forest* argF, expert_forest* resF);
    virtual ~saturation_evplus_op();

    void saturate(long ev, node_handle evmdd, long& resEv, node_handle& resEvmdd);
    void saturate(long ev, node_handle evmdd, int level, long& resEv, node_handle& resEvmdd);

    virtual bool isStaleEntry(const node_handle* entryData);
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData) const;

  protected:
    inline compute_table::search_key*
    findSaturateResult(long aev, node_handle a, int level, long& bev, node_handle& b) {
      compute_table::search_key* CTsrch = useCTkey();
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->reset();
      CTsrch->writeNH(a);
      if (argF->isFullyReduced()) CTsrch->write(level);
      compute_table::search_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      cacheFind.read(bev);
      bev += aev;
      b = resF->linkNode(cacheFind.readNH());
      doneCTkey(CTsrch);
      return 0;
    }
    inline node_handle saveSaturateResult(compute_table::search_key* Key,
      long aev, node_handle a, long bev, node_handle b)
    {
      argF->cacheNode(a);
      compute_table::entry_builder &entry = CT->startNewEntry(Key);
      entry.writeResult(bev - aev);
      entry.writeResultNH(resF->cacheNode(b));
      CT->addEntry();
      return b;
    }
};

// ******************************************************************
// *                                                                *
// *                      common_dfs_mt  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::common_dfs_mt : public binary_operation {
  public:
    common_dfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle*);
#endif
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData) const;
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual void saturateHelper(unpacked_node &mdd) = 0;

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
    void splitMxd(node_handle mxd);

  protected:
    node_handle* splits;
    binary_operation* mddUnion;
    binary_operation* mxdIntersection;
    binary_operation* mxdDifference;

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
// *                    common_dfs_evplus  class                    *
// *                                                                *
// ******************************************************************

class MEDDLY::common_dfs_evplus : public binary_operation {
  public:
    common_dfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual bool isStaleEntry(const node_handle* entryData);
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData) const;
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual void saturateHelper(unpacked_node &mdd) = 0;

  protected:
    inline compute_table::search_key*
    findResult(long aev, node_handle a, node_handle b, long& cev, node_handle& c)
    {
      compute_table::search_key* CTsrch = useCTkey();
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->reset();
      CTsrch->writeNH(a);
      CTsrch->writeNH(b);

      compute_table::search_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      cacheFind.read(cev);
      c = resF->linkNode(cacheFind.readNH());
      if (c != 0) {
        cev += aev;
      }
      doneCTkey(CTsrch);
      return 0;
    }
    inline void saveResult(compute_table::search_key* Key,
      long aev, node_handle a, node_handle b, long cev, node_handle c)
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
      compute_table::entry_builder &entry = CT->startNewEntry(Key);
      // Write long
      entry.writeResult(c == 0 ? 0L : cev - aev);
      entry.writeResultNH(resF->cacheNode(c));
      CT->addEntry();
    }
    void splitMxd(node_handle mxd);

  protected:
    node_handle* splits;
    binary_operation* evplusUnionMin;
    binary_operation* mxdIntersection;
    binary_operation* mxdDifference;

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
// *                     saturation_op  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_op
::saturation_op(common_dfs_mt* p, expert_forest* argF, expert_forest* resF)
  : unary_operation(saturation_opname::getInstance(),
  ((argF != 0 && argF->isFullyReduced())? 2: 1),
  1, argF, resF)
{
  parent = p;
}

MEDDLY::saturation_op::~saturation_op()
{
}

MEDDLY::node_handle MEDDLY::saturation_op::saturate(node_handle mdd)
{
  return saturate(mdd, argF->getNumVariables());
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
  compute_table::search_key* Key = findSaturateResult(mdd, k, n);
  if (0==Key) return n;

  const int sz = argF->getLevelSize(k);               // size
  const int mdd_level = argF->getNodeLevel(mdd);      // mdd level

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d, mdd_level: %d\n",
      mdd, k, sz, mdd_level);
#endif

  unpacked_node* C = unpacked_node::newFull(resF, k, sz);
  // Initialize mdd reader
  unpacked_node *mddDptrs = unpacked_node::useUnpackedNode();
  if (mdd_level < k) {
    mddDptrs->initRedundant(argF, k, mdd, true);
  } else {
    mddDptrs->initFromNode(argF, mdd, true);
  }

  // Do computation
  for (int i=0; i<sz; i++) {
    C->d_ref(i) = mddDptrs->d(i) ? saturate(mddDptrs->d(i), k-1) : 0;
  }

  // Cleanup
  unpacked_node::recycle(mddDptrs);

  parent->saturateHelper(*C);
  n = resF->createReducedNode(-1, C);

  // save in compute table
  saveSaturateResult(Key, mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  return n;
}

#ifndef USE_NODE_STATUS
bool MEDDLY::saturation_op::isStaleEntry(const node_handle* data)
{
  return (argF->isFullyReduced()
  ? (argF->isStale(data[0]) || resF->isStale(data[2]))
  : (argF->isStale(data[0]) || resF->isStale(data[1])));
}
#else
MEDDLY::forest::node_status
MEDDLY::saturation_op::getStatusOfEntry(const node_handle* data)
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

void MEDDLY::saturation_op::discardEntry(const node_handle* data)
{
  if (argF->isFullyReduced()) {
    argF->uncacheNode(data[0]);
    resF->uncacheNode(data[2]);
  } else {
    argF->uncacheNode(data[0]);
    resF->uncacheNode(data[1]);
  }
}

void MEDDLY::saturation_op::showEntry(output &strm, const node_handle* data) const
{
  if (argF->isFullyReduced()) {
    strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) << "): " << long(data[2]) << "]";
  } else {
    strm << "[" << getName() << "(" << long(data[0]) << "): " << long(data[1]) << "]";
  }
}


// ******************************************************************
// *                                                                *
// *                 saturation_evplus_op  methods                  *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_evplus_op
::saturation_evplus_op(common_dfs_evplus* p, expert_forest* argF, expert_forest* resF)
  : unary_operation(saturation_opname::getInstance(),
  (argF != 0 && argF->isFullyReduced())
      ? (sizeof(node_handle) + sizeof(int)) / sizeof(node_handle)
      : (sizeof(node_handle) / sizeof(node_handle)),
  (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle),
  argF, resF)
{
  parent = p;
}

MEDDLY::saturation_evplus_op::~saturation_evplus_op()
{
}

void MEDDLY::saturation_evplus_op::saturate(long ev, node_handle evmdd, long& resEv, node_handle& resEvmdd)
{
  saturate(ev, evmdd, argF->getNumVariables(), resEv, resEvmdd);
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
  compute_table::search_key* Key = findSaturateResult(ev, evmdd, k, resEv, resEvmdd);
  if (0==Key) {
    return;
  }

  const int sz = argF->getLevelSize(k);               // size
  const int evmdd_level = argF->getNodeLevel(evmdd);      // evmdd level

#ifdef DEBUG_DFS
  printf("evmdd: %d, level: %d, size: %d, evmdd_level: %d\n",
      evmdd, k, sz, evmdd_level);
#endif

  unpacked_node* C = unpacked_node::newFull(resF, k, sz);
  // Initialize evmdd reader
  unpacked_node *evmddDptrs = unpacked_node::useUnpackedNode();
  if (evmdd_level < k) {
    evmddDptrs->initRedundant(argF, k, evmdd, true);
  } else {
    evmddDptrs->initFromNode(argF, evmdd, true);
  }

  // Do computation
  for (int i=0; i<sz; i++) {
    long cev = Inf<long>();
    node_handle c = 0;
    if (evmddDptrs->d(i) != 0) {
      saturate(evmddDptrs->ei(i), evmddDptrs->d(i), k-1, cev, c);
    }
    C->setEdge(i, cev);
    C->d_ref(i) = c;
  }

  // Cleanup
  unpacked_node::recycle(evmddDptrs);

  parent->saturateHelper(*C);
  resF->createReducedNode(-1, C, resEv, resEvmdd);
  resEv += ev;

  // save in compute table
  saveSaturateResult(Key, ev, evmdd, resEv, resEvmdd);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, resEvmdd);
#endif
}

bool MEDDLY::saturation_evplus_op::isStaleEntry(const node_handle* data)
{
  return (argF->isFullyReduced()
  ? (argF->isStale(data[0]) || resF->isStale(data[(sizeof(node_handle) + sizeof(int) + sizeof(long)) / sizeof(node_handle)]))
  : (argF->isStale(data[0]) || resF->isStale(data[(sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)])));
}

void MEDDLY::saturation_evplus_op::discardEntry(const node_handle* data)
{
  if (argF->isFullyReduced()) {
    argF->uncacheNode(data[0]);
    resF->uncacheNode(data[(sizeof(node_handle) + sizeof(int) + sizeof(long)) / sizeof(node_handle)]);
  } else {
    argF->uncacheNode(data[0]);
    resF->uncacheNode(data[(sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)]);
  }
}

void MEDDLY::saturation_evplus_op::showEntry(output &strm, const node_handle* data) const
{
  if (argF->isFullyReduced()) {
    strm << "[" << getName()
      << "(" << long(data[0])
      << ", " << long(data[sizeof(node_handle) / sizeof(node_handle)])
      << "): " << long(data[(sizeof(node_handle) + sizeof(int) + sizeof(long)) / sizeof(node_handle)])
      << "]";
  } else {
    strm << "[" << getName()
      << "(" << long(data[0])
      << "): " << long(data[(sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)])
      << "]";
  }
}


// ******************************************************************
// *                                                                *
// *                     common_dfs_mt  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::common_dfs_mt::common_dfs_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc, 2, 1, a1, a2, res)
{
  splits = 0;
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;
  freeqs = 0;
  freebufs = 0;
}

#ifndef USE_NODE_STATUS
bool MEDDLY::common_dfs_mt::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}
#else
MEDDLY::forest::node_status
MEDDLY::common_dfs_mt::getStatusOfEntry(const node_handle* data)
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

void MEDDLY::common_dfs_mt::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void MEDDLY::common_dfs_mt::showEntry(output &strm, const node_handle* data) const
{
  strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) << "): " << long(data[2]) << "]";
}

#ifdef BATCHED_UNIONS
std::vector<std::vector<std::deque<MEDDLY::node_handle>>> g_union_nodes;
#endif
#ifdef COUNT_UNION_LENGTH
// print histogram
std::vector<long> union_length_counts;
std::vector<std::vector<long>> g_counts;
#endif

void MEDDLY::common_dfs_mt
::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c)
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
  b.show(stdout, 2);
#endif

  // Partition NSF by levels
  splitMxd(b.getNode());

#ifdef BATCHED_UNIONS
  // Initialize sets for each level
  g_union_nodes.clear();
  g_union_nodes.push_back(std::vector<std::deque<node_handle>>());
  for (int i = 1; i <= resF->getNumVariables(); i++) {
    std::vector<std::deque<node_handle>> x(resF->getVariableSize(i));
    g_union_nodes.push_back(x);
  }
#endif

#ifdef COUNT_UNION_LENGTH
  union_length_counts.clear();
  union_length_counts.resize(10, 0);
  g_counts.clear();
  g_counts.push_back(std::vector<long>());
  for (int i = 1; i <= resF->getNumVariables(); i++) {
    std::vector<long> x(resF->getVariableSize(i), 0);
    g_counts.push_back(x);
  }
  assert(g_counts.size() == resF->getNumVariables()+1);
#endif

  // Execute saturation operation
  saturation_op *so = new saturation_op(this, arg1F, resF);
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
  so->removeAllComputeTableEntries();
  for (int i = arg2F->getNumVariables(); i; i--) arg2F->unlinkNode(splits[i]);
  delete[] splits;
  splits = 0;
  delete so;

#ifdef COUNT_UNION_LENGTH
  // print histogram
  printf("\nLength :\tFrequency\tFreq*Length\tPercentage of Saturation Unions\n");
  double total_unions = 0;
  for (int i = 0; i < union_length_counts.size(); i++) {
    total_unions += i*union_length_counts[i];
  }
  if (total_unions < 1) total_unions = 1;
  for (int i = 0; i < union_length_counts.size(); i++) {
    printf("\t%d :\t%*ld\t%*ld\t%*.*lf\n", i, 8, union_length_counts[i],
        8, union_length_counts[i]*i,
        8, 1, union_length_counts[i]*i*100.0/total_unions);
  }

  union_length_counts.clear();
  g_counts.clear();
#endif

#ifdef BATCHED_UNIONS
  g_union_nodes.clear();
#endif

#ifdef BATCHED_UNLINK
  for (auto n : to_be_unlinked) resF->unlinkNode(n);
  to_be_unlinked.clear();
#endif
}

// Partition the nsf based on "top level"
void MEDDLY::common_dfs_mt::splitMxd(node_handle mxd)
{
  MEDDLY_DCASSERT(arg2F);
  MEDDLY_DCASSERT(0==splits);
  splits = new node_handle[arg2F->getNumVariables()+1];
  splits[0] = 0;

  // we'll be unlinking later, so...
  arg2F->linkNode(mxd);

  // Build from top down
  for (int level = arg2F->getNumVariables(); level; level--) {

    if (0==mxd) {
      // common and easy special case
      splits[level] = 0;
      continue;
    }

    int mxdLevel = arg2F->getNodeLevel(mxd);
    MEDDLY_DCASSERT(ABS(mxdLevel) <= level);

    // Initialize readers
    unpacked_node *Mu = unpacked_node::useUnpackedNode();
    unpacked_node *Mp = unpacked_node::useUnpackedNode();
    if (isLevelAbove(level, mxdLevel)) {
      Mu->initRedundant(arg2F, level, mxd, true);
    } else {
      Mu->initFromNode(arg2F, mxd, true);
    }

    bool first = true;
    node_handle maxDiag;

    // Read "rows"
    for (int i=0; i<Mu->getSize(); i++) {
      // Initialize column reader
      int mxdPLevel = arg2F->getNodeLevel(Mu->d(i));
      if (isLevelAbove(-level, mxdPLevel)) {
        Mp->initIdentity(arg2F, -level, i, Mu->d(i), true);
      } else {
        Mp->initFromNode(arg2F, Mu->d(i), true);
      }

      // Intersect along the diagonal
      if (first) {
        maxDiag = arg2F->linkNode(Mp->d(i));
        first = false;
      } else {
        node_handle nmd = mxdIntersection->compute(maxDiag, Mp->d(i));
        arg2F->unlinkNode(maxDiag);
        maxDiag = nmd;
      }

      // cleanup
    } // for i

    // maxDiag is what we can split from here
    splits[level] = mxdDifference->compute(mxd, maxDiag);
    arg2F->unlinkNode(mxd);
    mxd = maxDiag;

    // Cleanup
    unpacked_node::recycle(Mp);
    unpacked_node::recycle(Mu);
  } // for level

#ifdef DEBUG_SPLIT
  printf("After splitting monolithic event in msat\n");
  printf("splits array: [");
  for (int k=0; k <= arg2F->getNumVariables(); k++) {
    if (k) printf(", ");
    printf("%d", splits[k]);
  }
  printf("]\n");
#endif
}

// ******************************************************************
// *                 common_dfs_mt::indexq  methods                 *
// ******************************************************************

MEDDLY::common_dfs_mt::indexq::indexq()
{
  data = 0;
  size = 0;
  head = NULPTR;
}

MEDDLY::common_dfs_mt::indexq::~indexq()
{
  free(data);
}

void MEDDLY::common_dfs_mt::indexq::resize(int sz)
{
  if (sz <= size) return;
  data = (int*) realloc(data, sz * sizeof(int));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);

  for (; size < sz; size++) data[size] = NOTINQ;
}

// ******************************************************************
// *                 common_dfs_mt::charbuf methods                 *
// ******************************************************************

MEDDLY::common_dfs_mt::charbuf::charbuf()
{
  data = 0;
  size = 0;
}

MEDDLY::common_dfs_mt::charbuf::~charbuf()
{
  free(data);
}

void MEDDLY::common_dfs_mt::charbuf::resize(int sz)
{
  if (sz <= size) return;
  data = (char*) realloc(data, sz * sizeof(char));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *                       forwd_dfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_dfs_mt : public common_dfs_mt {
  public:
    forwd_dfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);
    ~forwd_dfs_mt();

  protected:
    virtual void saturateHelper(unpacked_node &mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);

#ifdef BATCHED_UNIONS
    // Sort a queue of nodes by node count
    void sort(bool descending_order, int level, std::deque<node_handle>& q);
    unary_operation* mdd_op;
    sorter *mdd_sorter;
#endif
};

MEDDLY::forwd_dfs_mt::forwd_dfs_mt(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : common_dfs_mt(opcode, arg1, arg2, res)
{
#ifdef BATCHED_UNIONS
  mdd_op = getOperation(COMPARISON_OPERATION, resF, MEDDLY_RESULT_TYPE);
  MEDDLY_DCASSERT(mdd_op);
  mdd_sorter = new sorter(mdd_op);
#endif
}

MEDDLY::forwd_dfs_mt::~forwd_dfs_mt()
{
#ifdef BATCHED_UNIONS
  delete mdd_sorter;
#endif
}

#ifdef BATCHED_UNIONS

// Sort a queue of nodes by cardinality
void MEDDLY::forwd_dfs_mt::sort(bool descending_order,
    int level, std::deque<node_handle>& q)
{
  mdd_sorter->setLevel(level);
  mdd_sorter->setDescending(descending_order);
  std::sort(q.begin(), q.end(), *mdd_sorter);
}

#endif

#ifdef BATCHED_UNIONS_FOR_SATURATE

// SaturateHelper with batched unions
void MEDDLY::forwd_dfs_mt::saturateHelper(unpacked_node &nb)
{
  node_handle mxd = splits[nb.getLevel()];
  if (mxd == 0) return;

  const int mxdLevel = arg2F->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node *Ru = unpacked_node::useUnpackedNode();
  unpacked_node *Rp = unpacked_node::useUnpackedNode();
  if (mxdLevel < 0) {
    Ru->initRedundant(arg2F, nb.getLevel(), mxd, true);
  } else {
    Ru->initFromNode(arg2F, mxd, true);
  }

  // batches of unions per node index
  std::vector<std::deque<node_handle>>& union_nodes = g_union_nodes[nb.getLevel()];
  if (union_nodes.size() < nb.getSize()) union_nodes.resize(nb.getSize());
  for (int i = 0; i < union_nodes.size(); i++) assert(union_nodes[i].empty());

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.d(i)) {
      queue->add(i);
      // move d(i) to the associated queue, and reset d(i).
      // otherwise, batched-union algorithm will not start.
      union_nodes[i].push_back(nb.d(i));
      nb.d_ref(i) = 0;
    }
  }

  // explore indexes
  while (!queue->isEmpty()) {
    const int i = queue->remove();

    // accumulate nb.d(i) and all nodes in union_nodes[i]
    node_handle acc = nb.d(i);
    resF->linkNode(acc);
    if (!union_nodes[i].empty()) {
      // Sort queue by COMPARISON_OPERATION
#ifdef DQ
      union_nodes[i].push_back(acc);
      acc = 0;
#endif
      sort(DESCENDING_ORDER, nb.getLevel(), union_nodes[i]);
      // Add the nodes to acc
      for (auto mdd : union_nodes[i]) {
        // printf("%d\n", mdd);
        node_handle temp = mddUnion->compute(acc, mdd);
        resF->unlinkNode(acc);
#ifdef BATCHED_UNLINK
        to_be_unlinked.push_back(mdd);
#else
        resF->unlinkNode(mdd);
#endif
        acc = temp;
      }
      union_nodes[i].clear();
    }

    if (acc == nb.d(i)) {
      resF->unlinkNode(acc);
      continue; // nothing new to explore in acc
    }

    // update d(i) with the result of accumulation.
    resF->unlinkNode(nb.d(i));
    nb.d_ref(i) = acc;

    MEDDLY_DCASSERT(nb.d(i));
    if (0==Ru->d(i)) continue;  // row i is empty

    // grab column (TBD: build these ahead of time?)
    const int dlevel = arg2F->getNodeLevel(Ru->d(i));

    if (dlevel == -nb.getLevel()) {
      Rp->initFromNode(arg2F, Ru->d(i), false);
    } else {
      Rp->initIdentity(arg2F, -nb.getLevel(), i, Ru->d(i), false);
    }

    for (int jz=0; jz<Rp->getNNZs(); jz++) {
      const int j = Rp->i(jz);
      if (-1==nb.d(j)) continue;  // nothing can be added to this set

      node_handle rec = recFire(nb.d(i), Rp->d(jz));

      if (0 == rec) continue;
      if (rec == nb.d(j)) { 
        resF->unlinkNode(rec); 
        continue; 
      }

#if 0
      if (0 == nb.d(j)) {
        nb.d_ref(j) = rec;
      }
      else if (-1 == rec) {
        resF->unlinkNode(nb.d(j));
        nb.d_ref(j) = -1;
        // nothing can be added to nb.d(j), therefore, empty its union_queue
        for (auto n : union_nodes[j]) resF->unlinkNode(n);
        union_nodes[j].clear();
      }
      else {
        union_nodes[j].push_back(rec);
      }
      queue->add(j);
#else
      queue->add(j);
      union_nodes[j].push_back(rec);
#endif
    } // for j

  } // while there are indexes to explore

  // cleanup
  unpacked_node::recycle(Rp);
  unpacked_node::recycle(Ru);
  recycle(queue);
}

#else

// SaturateHelper with immediate unions
void MEDDLY::forwd_dfs_mt::saturateHelper(unpacked_node &nb)
{
  node_handle mxd = splits[nb.getLevel()];
  if (mxd == 0) return;

  const int mxdLevel = arg2F->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node *Ru = unpacked_node::useUnpackedNode();
  unpacked_node *Rp = unpacked_node::useUnpackedNode();
  if (mxdLevel < 0) {
    Ru->initRedundant(arg2F, nb.getLevel(), mxd, true);
  } else {
    Ru->initFromNode(arg2F, mxd, true);
  }

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.d(i)) queue->add(i);
  }

#ifdef COUNT_UNION_LENGTH
  std::vector<long>& counts = g_counts[nb.getLevel()];
  if (counts.size() < nb.getSize()) counts.resize(nb.getSize());
  std::fill(counts.begin(), counts.end(), 0);
#endif

  // explore indexes
  while (!queue->isEmpty()) {
    const int i = queue->remove();

#ifdef COUNT_UNION_LENGTH
    if (0 != counts[i]) {
      union_length_counts[(counts[i] > 9? 9: counts[i])]++;
      counts[i] = 0;
    }
#endif

    MEDDLY_DCASSERT(nb.d(i));
    if (0==Ru->d(i)) continue;  // row i is empty

    // grab column (TBD: build these ahead of time?)
    const int dlevel = arg2F->getNodeLevel(Ru->d(i));

    if (dlevel == -nb.getLevel()) {
      Rp->initFromNode(arg2F, Ru->d(i), false);
    } else {
      Rp->initIdentity(arg2F, -nb.getLevel(), i, Ru->d(i), false);
    }

    for (int jz=0; jz<Rp->getNNZs(); jz++) {
      const int j = Rp->i(jz);
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
#ifdef BATCHED_UNLINK
        to_be_unlinked.push_back(rec);
#else
        resF->unlinkNode(rec);
#endif
        if (acc != nb.d(j)) {
          resF->unlinkNode(nb.d(j));
          nb.d_ref(j) = acc;
        } else {
          resF->unlinkNode(acc);
          updated = false;
        }
      }

      if (updated) {
        if (j == i) {
          // Restart inner for-loop.
          jz = -1;
        } else {
          queue->add(j);
#ifdef COUNT_UNION_LENGTH
          counts[j]++;
#endif
        }
      }

    } // for j

  } // while there are indexes to explore

  // cleanup
  unpacked_node::recycle(Rp);
  unpacked_node::recycle(Ru);
  recycle(queue);
}

#endif  // End of SaturateHelper

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

#ifdef BATCHED_UNIONS_FOR_REC_FIRE
    // batches of unions per node index
    std::vector<std::deque<node_handle>>& union_nodes = g_union_nodes[nb->getLevel()];
    if (union_nodes.size() < nb->getSize()) union_nodes.resize(nb->getSize());
    // for (int i = 0; i < union_nodes.size(); i++) assert(union_nodes[i].empty());
#endif

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

#ifdef COUNT_UNION_LENGTH
    std::vector<long>& counts = g_counts[nb->getLevel()];
    if (counts.size() < rSize) counts.resize(rSize);
    std::fill(counts.begin(), counts.end(), 0);
#endif

    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);

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
#ifdef BATCHED_UNIONS_FOR_REC_FIRE
        // batched union
        /*
        if (0 == nb->d(j) && union_nodes[j].empty()) {
          nb->d_ref(j) = newstates;
          continue;
        }
        */
        union_nodes[j].push_back(newstates);
#else
        // no batches
        if (0==nb->d(j)) {
          nb->d_ref(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        const int oldj = nb->d(j);
        nb->d_ref(j) = mddUnion->compute(newstates, oldj);
        resF->unlinkNode(oldj);
#ifdef BATCHED_UNLINK
        to_be_unlinked.push_back(newstates);
#else
        resF->unlinkNode(newstates);
#endif
#endif

#ifdef COUNT_UNION_LENGTH
        counts[j]++;
#endif
      } // for j
  
    } // for i

#ifdef COUNT_UNION_LENGTH
    for (int i = 0; i < counts.size(); i++) {
      if (0 != counts[i]) {
        union_length_counts[(counts[i] > 9? 9: counts[i])]++;
      }
    }
#endif

#ifdef BATCHED_UNIONS_FOR_REC_FIRE
    for (int i = 0; i < union_nodes.size(); i++) {
      assert(0 == nb->d(i));
      node_handle acc = 0;
      if (!union_nodes[i].empty()) {
        // Sort queue by COMPARISON_OPERATION
        sort(DESCENDING_ORDER, nb->getLevel(), union_nodes[i]);
        // Add the nodes to acc
        for (auto mdd : union_nodes[i]) {
          // printf("%d\n", mdd);
          node_handle temp = mddUnion->compute(acc, mdd);
          resF->unlinkNode(acc);
#ifdef BATCHED_UNLINK
          to_be_unlinked.push_back(mdd);
#else
          resF->unlinkNode(mdd);
#endif
          acc = temp;
        }
        union_nodes[i].clear();
      }
      nb->d_ref(i) = acc;
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




// ******************************************************************
// *                                                                *
// *                       bckwd_dfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_dfs_mt : public common_dfs_mt {
  public:
    bckwd_dfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
};

MEDDLY::bckwd_dfs_mt::bckwd_dfs_mt(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : common_dfs_mt(opcode, arg1, arg2, res)
{
}

void MEDDLY::bckwd_dfs_mt::saturateHelper(unpacked_node& nb)
{
  node_handle mxd = splits[nb.getLevel()];
  if (mxd == 0) return;

  const int mxdLevel = arg2F->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node *Ru = unpacked_node::useUnpackedNode();
  unpacked_node *Rp = unpacked_node::useUnpackedNode();
  if (mxdLevel < 0) {
    Ru->initRedundant(arg2F, nb.getLevel(), mxd, false);
  } else {
    Ru->initFromNode(arg2F, mxd, false);
  }

  // indexes to explore
  charbuf* expl = useCharBuf(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) expl->data[i] = 2;
  bool repeat = true;

  // explore 
  while (repeat) {
    // "advance" the explore list
    for (int i=0; i<nb.getSize(); i++) if (expl->data[i]) expl->data[i]--;
    repeat = false;

    // explore all rows
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);
      // grab column (TBD: build these ahead of time?)
      const int dlevel = arg2F->getNodeLevel(Ru->d(iz));

      if (dlevel == -nb.getLevel()) {
        Rp->initFromNode(arg2F, Ru->d(iz), false);
      } else {
        Rp->initIdentity(arg2F, -nb.getLevel(), i, Ru->d(iz), false);
      }

      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        const int j = Rp->i(jz);
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
  } // while repeat
  // cleanup
  unpacked_node::recycle(Rp);
  unpacked_node::recycle(Ru);
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
  compute_table::search_key* Key = findResult(mdd, mxd, result);
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
      const int i = Ru->i(iz);
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru->d(iz)))) {
        Rp->initIdentity(arg2F, rLevel, i, Ru->d(iz), false);
      } else {
        Rp->initFromNode(arg2F, Ru->d(iz), false);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        const int j = Rp->i(jz);
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
        node_handle oldi = nb->d(i);
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
// *                     common_dfs_evplus class                    *
// *                                                                *
// ******************************************************************

MEDDLY::common_dfs_evplus::common_dfs_evplus(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc,
    (sizeof(node_handle) + sizeof(node_handle)) / sizeof(node_handle),
    (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle),
    a1, a2, res)
{
  splits = 0;
  evplusUnionMin = 0;
  mxdIntersection = 0;
  mxdDifference = 0;
  freeqs = 0;
  freebufs = 0;
}

bool MEDDLY::common_dfs_evplus::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[sizeof(node_handle) / sizeof(node_handle)]) ||
         resF->isStale(data[(2 * sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)]);
}

void MEDDLY::common_dfs_evplus::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[sizeof(node_handle) / sizeof(node_handle)]);
  resF->uncacheNode(data[(2 * sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)]);
}

void MEDDLY::common_dfs_evplus::showEntry(output &strm, const node_handle* data) const
{
  strm << "[" << getName()
    << "(" << long(data[0])
    << ", " << long(data[sizeof(node_handle) / sizeof(node_handle)])
    << "): " << long(data[(2 * sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)])
    << "]";
}

void MEDDLY::common_dfs_evplus
::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  // Initialize operations
  evplusUnionMin = getOperation(UNION, resF, resF, resF);
  MEDDLY_DCASSERT(evplusUnionMin);

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
  b.show(stdout, 2);
#endif

  // Partition NSF by levels
  splitMxd(b.getNode());

  // Execute saturation operation
  saturation_evplus_op *so = new saturation_evplus_op(this, arg1F, resF);
  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long cev = Inf<long>();
  node_handle cnode = 0;
  so->saturate(aev, a.getNode(), cev, cnode);
  c.set(cnode, cev);

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
  so->removeAllComputeTableEntries();
  for (int i = arg2F->getNumVariables(); i; i--) arg2F->unlinkNode(splits[i]);
  delete[] splits;
  splits = 0;
  delete so;
}

// Partition the nsf based on "top level"
void MEDDLY::common_dfs_evplus::splitMxd(node_handle mxd)
{
  MEDDLY_DCASSERT(arg2F);
  MEDDLY_DCASSERT(0==splits);
  splits = new node_handle[arg2F->getNumVariables()+1];
  splits[0] = 0;

  // we'll be unlinking later, so...
  arg2F->linkNode(mxd);

  // Build from top down
  for (int level = arg2F->getNumVariables(); level; level--) {

    if (0==mxd) {
      // common and easy special case
      splits[level] = 0;
      continue;
    }

    int mxdLevel = arg2F->getNodeLevel(mxd);
    MEDDLY_DCASSERT(ABS(mxdLevel) <= level);

    // Initialize readers
    unpacked_node *Mu = unpacked_node::useUnpackedNode();
    unpacked_node *Mp = unpacked_node::useUnpackedNode();
    if (isLevelAbove(level, mxdLevel)) {
      Mu->initRedundant(arg2F, level, mxd, true);
    } else {
      Mu->initFromNode(arg2F, mxd, true);
    }

    bool first = true;
    node_handle maxDiag;

    // Read "rows"
    for (int i=0; i<Mu->getSize(); i++) {
      // Initialize column reader
      int mxdPLevel = arg2F->getNodeLevel(Mu->d(i));
      if (isLevelAbove(-level, mxdPLevel)) {
        Mp->initIdentity(arg2F, -level, i, Mu->d(i), true);
      } else {
        Mp->initFromNode(arg2F, Mu->d(i), true);
      }

      // Intersect along the diagonal
      if (first) {
        maxDiag = arg2F->linkNode(Mp->d(i));
        first = false;
      } else {
        node_handle nmd = mxdIntersection->compute(maxDiag, Mp->d(i));
        arg2F->unlinkNode(maxDiag);
        maxDiag = nmd;
      }

      // cleanup
    } // for i

    // maxDiag is what we can split from here
    splits[level] = mxdDifference->compute(mxd, maxDiag);
    arg2F->unlinkNode(mxd);
    mxd = maxDiag;

    // Cleanup
    unpacked_node::recycle(Mp);
    unpacked_node::recycle(Mu);
  } // for level

#ifdef DEBUG_SPLIT
  printf("After splitting monolithic event in msat\n");
  printf("splits array: [");
  for (int k=0; k <= arg2F->getNumVariables(); k++) {
    if (k) printf(", ");
    printf("%d", splits[k]);
  }
  printf("]\n");
#endif
}

// ******************************************************************
// *                 common_dfs_evplus::indexq  methods                 *
// ******************************************************************

MEDDLY::common_dfs_evplus::indexq::indexq()
{
  data = 0;
  size = 0;
  head = NULPTR;
}

MEDDLY::common_dfs_evplus::indexq::~indexq()
{
  free(data);
}

void MEDDLY::common_dfs_evplus::indexq::resize(int sz)
{
  if (sz <= size) return;
  data = (int*) realloc(data, sz * sizeof(int));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY);

  for (; size < sz; size++) data[size] = NOTINQ;
}

// ******************************************************************
// *                 common_dfs_evplus::charbuf methods                 *
// ******************************************************************

MEDDLY::common_dfs_evplus::charbuf::charbuf()
{
  data = 0;
  size = 0;
}

MEDDLY::common_dfs_evplus::charbuf::~charbuf()
{
  free(data);
}

void MEDDLY::common_dfs_evplus::charbuf::resize(int sz)
{
  if (sz <= size) return;
  data = (char*) realloc(data, sz * sizeof(char));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY);
}

// ******************************************************************
// *                                                                *
// *                     forwd_dfs_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_dfs_evplus : public common_dfs_evplus {
  public:
  forwd_dfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);
  protected:
    virtual void saturateHelper(unpacked_node &mdd);
    void recFire(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
};

MEDDLY::forwd_dfs_evplus::forwd_dfs_evplus(const binary_opname* opcode,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : common_dfs_evplus(opcode, arg1, arg2, res)
{
}

void MEDDLY::forwd_dfs_evplus::saturateHelper(unpacked_node &nb)
{
  node_handle mxd = splits[nb.getLevel()];
  if (mxd == 0) return;

  const int mxdLevel = arg2F->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node *Ru = unpacked_node::useUnpackedNode();
  unpacked_node *Rp = unpacked_node::useUnpackedNode();
  if (mxdLevel < 0) {
    Ru->initRedundant(arg2F, nb.getLevel(), mxd, true);
  } else {
    Ru->initFromNode(arg2F, mxd, true);
  }

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.d(i) != 0) {
      queue->add(i);
    }
  }

  // explore indexes
  while (!queue->isEmpty()) {
    const int i = queue->remove();

    MEDDLY_DCASSERT(nb.d(i));
    if (0==Ru->d(i)) continue;  // row i is empty

    // grab column (TBD: build these ahead of time?)
    const int dlevel = arg2F->getNodeLevel(Ru->d(i));

    if (dlevel == -nb.getLevel()) {
      Rp->initFromNode(arg2F, Ru->d(i), false);
    } else {
      Rp->initIdentity(arg2F, -nb.getLevel(), i, Ru->d(i), false);
    }

    for (int jz=0; jz<Rp->getNNZs(); jz++) {
      const int j = Rp->i(jz);

      long recev = Inf<long>();
      node_handle rec = 0;
      recFire(nb.ei(i), nb.d(i), Rp->d(jz), recev, rec);

      if (rec == 0) continue;

      // Increase the distance
      recev++;

      if (rec == nb.d(j)) {
        if (recev < nb.ei(j)) {
          nb.setEdge(j, recev);
        }
        resF->unlinkNode(rec);
        continue;
      }

      bool updated = true;

      if (0 == nb.d(j)) {
        nb.setEdge(j, recev);
        nb.d_ref(j) = rec;
      }
//      else if (rec == -1) {
//        resF->unlinkNode(nb.d(j));
//        nb.d_ref(j) = -1;
//      }
      else {
        long accev = Inf<long>();
        node_handle acc = 0;
        evplusUnionMin->compute(nb.ei(j), nb.d(j), recev, rec, accev, acc);
        resF->unlinkNode(rec);
        if (acc != nb.d(j)) {
          resF->unlinkNode(nb.d(j));
          nb.setEdge(j, accev);
          nb.d_ref(j) = acc;
        } else {
          MEDDLY_DCASSERT(accev == nb.ei(j));
          resF->unlinkNode(acc);
          updated = false;
        }
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
  unpacked_node::recycle(Rp);
  unpacked_node::recycle(Ru);
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
  compute_table::search_key* Key = findResult(ev, evmdd, mxd, resEv, resEvmdd);
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
  const int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize evmdd reader
  unpacked_node *A = unpacked_node::useUnpackedNode();
  if (evmddLevel < rLevel) {
    A->initRedundant(arg1F, rLevel, evmdd, true);
  } else {
    A->initFromNode(arg1F, evmdd, true);
  }

  if (evmddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i=0; i<rSize; i++) {
      long nev = Inf<long>();
      node_handle n = 0;
      recFire(A->ei(i) + ev, A->d(i), mxd, nev, n);
      nb->setEdge(i, nev);
      nb->d_ref(i) = n;
    }

  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= evmddLevel);

    // clear out result (important!)
    for (int i=0; i<rSize; i++) {
      nb->setEdge(i, 0L);
      nb->d_ref(i) = 0;
    }

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
      const int i = Ru->i(iz);
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
        long nev = Inf<long>();
        node_handle n = 0;
        recFire(A->ei(i) + ev, A->d(i), Rp->d(jz), nev, n);

        if (0==n) continue;
        if (0==nb->d(j)) {
          nb->setEdge(j, nev);
          nb->d_ref(j) = n;
          continue;
        }

        // there's new states and existing states; union them.
        long oldev = nb->ei(j);
        node_handle oldj = nb->d(j);
        long newev = Inf<long>();
        node_handle newstates = 0;
        evplusUnionMin->compute(nev, n, oldev, oldj, newev, newstates);
        nb->setEdge(j, newev);
        nb->d_ref(j) = newstates;

        resF->unlinkNode(oldj);
        resF->unlinkNode(n);
      } // for j

    } // for i

    unpacked_node::recycle(Rp);
    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);

  saturateHelper(*nb);
  resF->createReducedNode(-1, nb, resEv, resEvmdd);
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
// *                     forwd_dfs_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_dfs_opname : public binary_opname {
  public:
    forwd_dfs_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::forwd_dfs_opname::forwd_dfs_opname()
 : binary_opname("ReachableDFS")
{
}

MEDDLY::binary_operation* 
MEDDLY::forwd_dfs_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (a1->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    return new forwd_dfs_mt(this, a1, a2, r);
  }
  else if (a1->getEdgeLabeling() == forest::EVPLUS) {
    return new forwd_dfs_evplus(this, a1, a2, r);
  }
  else {
    throw error(error::TYPE_MISMATCH);
  }
}



// ******************************************************************
// *                                                                *
// *                     bckwd_dfs_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_dfs_opname : public binary_opname {
  public:
    bckwd_dfs_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::bckwd_dfs_opname::bckwd_dfs_opname()
 : binary_opname("ReverseReachableDFS")
{
}

MEDDLY::binary_operation* 
MEDDLY::bckwd_dfs_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (a1 != r)
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) 
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  return new bckwd_dfs_mt(this, a1, a2, r);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeForwardDFS()
{
  return new forwd_dfs_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeBackwardDFS()
{
  return new bckwd_dfs_opname;
}


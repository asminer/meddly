
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


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "reach_dfs.h"

// #define TRACE_RECFIRE
// #define DEBUG_DFS
// #define DEBUG_INITIAL
// #define DEBUG_NSF
// #define DEBUG_SPLIT

namespace MEDDLY {
  class saturation_opname;
  class saturation_op;

  class common_dfs_mt;

  class forwd_dfs_mt;
  class bckwd_dfs_mt;

  class forwd_dfs_opname;
  class bckwd_dfs_opname;
};


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

    long saturate(long mdd);

    virtual bool isStaleEntry(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int* entryData) const;

  protected:
    inline bool findSaturateResult(long a, long& b) {
      CTsrch.key(0) = a;
      const int* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      b = resF->linkNode(cacheFind[2]);
      return true;
    }
    inline long saveSaturateResult(long a, long b) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      entry.key(0) = argF->cacheNode(a);
      entry.result(0) = resF->cacheNode(b);
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

    virtual bool isStaleEntry(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int* entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual void saturateHelper(node_builder& mdd) = 0;

  protected:
    inline bool findResult(long a, long b, long &c) {
      CTsrch.key(0) = a;
      CTsrch.key(1) = b;
      const int* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[2]);
      return true;
    }
    inline long saveResult(long a, long b, long c) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      entry.key(0) = arg1F->cacheNode(a); 
      entry.key(1) = arg2F->cacheNode(b);
      entry.result(0) = resF->cacheNode(c);
      CT->addEntry();
      return c;
    }
    void splitMxd(long mxd);

  protected:
    long* splits;
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
// *                     saturation_op  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_op
::saturation_op(common_dfs_mt* p, expert_forest* argF, expert_forest* resF)
  : unary_operation(saturation_opname::getInstance(), 1, 1, argF, resF)
{
  parent = p;
}

MEDDLY::saturation_op::~saturation_op()
{
}

long MEDDLY::saturation_op::saturate(long mdd)
{
#ifdef DEBUG_DFS
  printf("mdd: %d\n", mdd);
#endif

  // MEDDLY_DCASSERT(argF->isReducedNode(mdd));

  // terminal condition for recursion
  if (argF->isTerminalNode(mdd)) return mdd;

  // search compute table
  long n = 0;
  if (findSaturateResult(mdd, n)) {
    resF->linkNode(n);
    return n;
  }

  int k = argF->getNodeLevel(mdd);      // level
  int sz = argF->getLevelSize(k);       // size

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d\n", mdd, k, sz);
#endif

  node_builder& nb = resF->useNodeBuilder(k, sz);
  node_reader* mddDptrs = argF->initNodeReader(mdd, true);
  for (int i=0; i<sz; i++) {
    nb.d(i) = mddDptrs->d(i) ? saturate(mddDptrs->d(i)) : 0;
  }
  node_reader::recycle(mddDptrs);
  parent->saturateHelper(nb);
  n = resF->createReducedNode(-1, nb);

  // save in compute table
  saveSaturateResult(mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  return n;
}

bool MEDDLY::saturation_op::isStaleEntry(const int* data)
{
  return argF->isStale(data[0]) ||
         resF->isStale(data[1]);
}

void MEDDLY::saturation_op::discardEntry(const int* data)
{
  argF->uncacheNode(data[0]);
  resF->uncacheNode(data[1]);
}

void MEDDLY::saturation_op::showEntry(FILE* strm, const int* data) const
{
  fprintf(strm, "[%s(%d): %d]", getName(), data[0], data[1]);
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
  binary_operation* mddUnion = 0;
  binary_operation* mxdIntersection = 0;
  binary_operation* mxdDifference = 0;
  freeqs = 0;
  freebufs = 0;
}

bool MEDDLY::common_dfs_mt::isStaleEntry(const int* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}

void MEDDLY::common_dfs_mt::discardEntry(const int* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void MEDDLY::common_dfs_mt::showEntry(FILE* strm, const int* data) const
{
  fprintf(strm, "[%s(%d, %d): %d]", getName(), data[0], data[1], data[2]);
}

void MEDDLY::common_dfs_mt
::compute(const dd_edge &a, const dd_edge &b, dd_edge &c)
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

  // Execute saturation operation
  saturation_op *so = new saturation_op(this, arg1F, resF);
  long cnode = so->saturate(a.getNode());
  c.set(cnode, 0, resF->getNodeLevel(cnode));

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
void MEDDLY::common_dfs_mt::splitMxd(long mxd)
{
  MEDDLY_DCASSERT(arg2F);
  MEDDLY_DCASSERT(0==splits);
  splits = new long[arg2F->getNumVariables()+1];
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
    MEDDLY_DCASSERT(ABS(mxdLevel <= level));

    // Initialize readers
    node_reader* Mu = isLevelAbove(level, mxdLevel)
      ? arg2F->initRedundantReader(level, mxd, true)
      : arg2F->initNodeReader(mxd, true);
    node_reader* Mp = node_reader::useReader();

    bool first = true;
    long maxDiag;

    // Read "rows"
    for (int i=0; i<Mu->getSize(); i++) {
      // Initialize column reader
      int mxdPLevel = arg2F->getNodeLevel(Mu->d(i));
      if (isLevelAbove(-level, mxdPLevel)) {
        arg2F->initIdentityReader(*Mp, -level, i, Mu->d(i), true);
      } else {
        arg2F->initNodeReader(*Mp, Mu->d(i), true);
      }

      // Intersect along the diagonal
      if (first) {
        maxDiag = arg2F->linkNode(Mp->d(i));
        first = false;
      } else {
        long nmd = mxdIntersection->compute(maxDiag, Mp->d(i));
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
    node_reader::recycle(Mp);
    node_reader::recycle(Mu);
  } // for level

#ifdef DEBUG_SPLIT
  for (int k=arg2F->getNumVariables(); k; k--) {
    if (splits[k]) {
      printf("------------------------------------------------------------\n");
      printf("Level %d nsf: %d\n", k, splits[k]);
      printf("------------------------------------------------------------\n");
      arg2F->showNodeGraph(stdout, splits[k]);
    }
  }
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
    throw error(error::INSUFFICIENT_MEMORY);

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
    throw error(error::INSUFFICIENT_MEMORY);
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
  protected:
    virtual void saturateHelper(node_builder& mdd);
    long recFire(long mdd, long mxd);
};

MEDDLY::forwd_dfs_mt::forwd_dfs_mt(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : common_dfs_mt(opcode, arg1, arg2, res)
{
}

void MEDDLY::forwd_dfs_mt::saturateHelper(node_builder& nb)
{
  long mxd = splits[nb.getLevel()];
  if (mxd == 0) return;

  int mxdLevel = arg2F->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  node_reader* Ru = (mxdLevel<0)
    ? arg2F->initRedundantReader(nb.getLevel(), mxd, true)
    : arg2F->initNodeReader(mxd, true);
  node_reader* Rp = node_reader::useReader();

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.d(i)) queue->add(i);
  }

  // explore indexes
  while (!queue->isEmpty()) {
    int i = queue->remove();

    MEDDLY_DCASSERT(nb.d(i));
    if (0==Ru->d(i)) continue;  // row i is empty

    // grab column (TBD: build these ahead of time?)
    int dlevel = arg2F->getNodeLevel(Ru->d(i));

    if (dlevel == -nb.getLevel()) {
      arg2F->initNodeReader(*Rp, Ru->d(i), false);
    } else {
      arg2F->initIdentityReader(*Rp, -nb.getLevel(), i, Ru->d(i), false);
    }

    for (int jz=0; jz<Rp->getNNZs(); jz++) {
      int j = Rp->i(jz);
      if (-1==nb.d(j)) continue;  // nothing can be added to this set

      long rec = recFire(nb.d(i), Rp->d(jz));

      if (rec == 0) continue;
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
        long acc = mddUnion->compute(nb.d(j), rec);
        resF->unlinkNode(rec);
        if (acc != nb.d(j)) {
          resF->unlinkNode(nb.d(j));
          nb.d(j) = acc;
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
        }
      }

    } // for j

  } // while there are indexes to explore

  // cleanup
  node_reader::recycle(Rp);
  node_reader::recycle(Ru);
  recycle(queue);
}

// Same as post-image, except we saturate before reducing.
long MEDDLY::forwd_dfs_mt::recFire(long mdd, long mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (arg2F->isTerminalNode(mxd)) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->getTerminalNode(true);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  long result = 0;
  if (findResult(mdd, mxd, result)) {
    return result;
  }

#ifdef TRACE_RECFIRE
  printf("computing recFire(%d, %d)\n", mdd, mxd);
  printf("  node %3d ", mdd);
  arg1F->showNode(stdout, mdd, 1);
  printf("\n  node %3d ", mxd);
  arg2F->showNode(stdout, mxd, 1);
  printf("\n");
#endif

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
        long newstates = recFire(A->d(i), Rp->d(jz));
        if (0==newstates) continue;
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
  return saveResult(mdd, mxd, result); 
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
    virtual void saturateHelper(node_builder& mdd);
    long recFire(long mdd, long mxd);
};

MEDDLY::bckwd_dfs_mt::bckwd_dfs_mt(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : common_dfs_mt(opcode, arg1, arg2, res)
{
}

void MEDDLY::bckwd_dfs_mt::saturateHelper(node_builder& nb)
{
  long mxd = splits[nb.getLevel()];
  if (mxd == 0) return;

  int mxdLevel = arg2F->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  node_reader* Ru = (mxdLevel<0)
    ? arg2F->initRedundantReader(nb.getLevel(), mxd, false)
    : arg2F->initNodeReader(mxd, false);
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

    // explore all rows
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      int i = Ru->i(iz);
      // grab column (TBD: build these ahead of time?)
      int dlevel = arg2F->getNodeLevel(Ru->d(iz));

      if (dlevel == -nb.getLevel()) {
        arg2F->initNodeReader(*Rp, Ru->d(iz), false); 
      } else {
        arg2F->initIdentityReader(*Rp, -nb.getLevel(), i, Ru->d(iz), false);
      }

      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);
        if (0==expl->data[j]) continue;
        if (0==nb.d(j))       continue;
        // We have an i->j edge to explore
        long rec = recFire(nb.d(j), Rp->d(jz));

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
          long acc = mddUnion->compute(nb.d(i), rec);
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
  } // while repeat
  // cleanup
  node_reader::recycle(Rp);
  node_reader::recycle(Ru);
  recycle(expl);
}

long MEDDLY::bckwd_dfs_mt::recFire(long mdd, long mxd)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (arg2F->isTerminalNode(mxd)) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->getTerminalNode(true);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  long result = 0;
  if (findResult(mdd, mxd, result)) {
    return result;
  }

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
        long newstates = recFire(A->d(j), Rp->d(jz));
        if (0==newstates) continue;
        if (0==nb.d(i)) {
          nb.d(i) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        long oldi = nb.d(i);
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
  return saveResult(mdd, mxd, result); 
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
    throw error(error::DOMAIN_MISMATCH);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a2->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (r->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH);

  return new forwd_dfs_mt(this, a1, a2, r);
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
    throw error(error::DOMAIN_MISMATCH);

  if (a1 != r)
    throw error(error::FOREST_MISMATCH);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    (a1->getRangeType() != a2->getRangeType()) ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) 
  )
    throw error(error::TYPE_MISMATCH);

  return new bckwd_dfs_mt(this, a1, a2, r);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeForwardDFS(const settings &s)
{
  return new forwd_dfs_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeBackwardDFS(const settings &s)
{
  return new bckwd_dfs_opname;
}


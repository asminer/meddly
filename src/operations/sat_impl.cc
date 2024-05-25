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
#include "relation_node.h"
#include "sat_impl.h"
#include <typeinfo> // for "bad_cast" exception
#include <set>
#include <map>

#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../oper_satur.h"
#include "../sat_relations.h"
#include "../ops_builtin.h"


namespace MEDDLY {
    class saturation_impl_by_events_op;

    class common_impl_dfs_by_events_mt;
    class forwd_impl_dfs_by_events_mt;
};

// #define DEBUG_INITIAL
// #define DEBUG_IS_REACHABLE



// ******************************************************************
// *                                                                *
// *             saturation_impl_by_events_op  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_impl_by_events_op : public saturation_operation {
  common_impl_dfs_by_events_mt* parent;
public:
  saturation_impl_by_events_op(common_impl_dfs_by_events_mt* p,
          const char* name, forest* argF, forest* resF);
  virtual ~saturation_impl_by_events_op();

  virtual void compute(const dd_edge& a, dd_edge& r);
  node_handle saturate(node_handle mdd, int level);

  // for reachable state in constraint detection
  virtual bool isReachable(const dd_edge& mdd, const dd_edge& constr);
  bool isReachable(
    node_handle mdd,
    int k,
    node_handle constraint,
    node_handle& saturation_result);

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
  inline void recycleCTKey(ct_entry_key* CTsrch) {
    CT0->recycle(CTsrch);
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
// *            common_impl_dfs_by_events_mt  class                 *
// *                                                                *
// ******************************************************************

class MEDDLY::common_impl_dfs_by_events_mt : public saturation_operation {
public:
  common_impl_dfs_by_events_mt(const char* name, implicit_relation* rel);
  virtual ~common_impl_dfs_by_events_mt();

  virtual void compute(const dd_edge& a, dd_edge &c);
  virtual bool isReachable(const dd_edge& a, const dd_edge& constraint);
  virtual void saturateHelper(unpacked_node& mdd) = 0;
  // for detecting reachable state in constraint
  virtual bool saturateHelper(unpacked_node& mdd, node_handle constraint) = 0;

protected:
  inline ct_entry_key*
  findResult(node_handle a, rel_node_handle b, node_handle &c)
  {
    ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeN(a);
    CTsrch->writeL(b);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    c = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
  }
  inline void recycleCTKey(ct_entry_key* CTsrch) {
    CT0->recycle(CTsrch);
  }
  inline node_handle saveResult(ct_entry_key* Key,
                                node_handle a, rel_node_handle b, node_handle c)
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

  implicit_relation* rel;

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
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, size);
      if (NOTINQ != data[i]) return;
      if (NULPTR == head) {
        // empty list
        head = i;
      } else {
        // not empty list
          MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, tail, size);
        data[tail] = i;
      }
      tail = i;
      data[i] = NULPTR;
    }
    inline int remove() {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, head, size);
      int ans = head;
      head = data[head];
      data[ans] = NOTINQ;
      MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, ans, size);
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
// *               forwd_impl_dfs_by_events_mt class                *
// *                                                                *
// ******************************************************************


class MEDDLY::forwd_impl_dfs_by_events_mt : public common_impl_dfs_by_events_mt {
public:
  forwd_impl_dfs_by_events_mt(const char* name, implicit_relation* rel);
protected:
  virtual void saturateHelper(unpacked_node& mdd);
  node_handle recFire(node_handle mdd, rel_node_handle mxd);
  MEDDLY::node_handle recFireSet(node_handle mdd, std::vector<rel_node_handle> mxd);

  // for reachable state in constraint detection
  bool saturateHelper(
      unpacked_node& nb,
      node_handle constraint);
  bool recFire(
       node_handle mdd,
       rel_node_handle mxd,
       node_handle constraint,
       node_handle& result);
 };

MEDDLY::forwd_impl_dfs_by_events_mt::forwd_impl_dfs_by_events_mt(
        const char* name, implicit_relation* rel)
: common_impl_dfs_by_events_mt(name, rel)
{
}


void MEDDLY::forwd_impl_dfs_by_events_mt::saturateHelper(unpacked_node& nb)
{

  int nEventsAtThisLevel = rel->lengthForLevel(nb.getLevel());

  if (0 == nEventsAtThisLevel) return;

  // Initialize mxd readers, note we might skip the unprimed level
  const int level = nb.getLevel();
  node_handle* events = rel->arrayForLevel(level);
  relation_node** Ru = new relation_node*[nEventsAtThisLevel];
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    Ru[ei] = rel->nodeExists(events[ei]);
#ifdef DEVELOPMENT_CODE
    int eventLevel = Ru[ei]->getLevel();
    MEDDLY_DCASSERT(ABS(eventLevel) == level);
#endif
  }

  dd_edge nbdj(resF), newst(resF);

  domain* dm = resF->getDomain();

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.down(i)) {
      queue->add(i);
    }
  }

  // explore indexes
  while (!queue->isEmpty()) {
    int i = queue->remove();

    MEDDLY_DCASSERT(nb.down(i));

    bool is_union = rel->isUnionPossible(nb.getLevel(),i,Ru);

    if(!is_union)
      {
      for (int ei = 0; ei < nEventsAtThisLevel; ei++) {

        int j = Ru[ei]->nextOf(i);
        if(j==-1) continue;
        if (j < nb.getSize() && -1==nb.down(j)) continue; // nothing can be added to this set

        node_handle rec = recFire(nb.down(i), Ru[ei]->getDown());

        if (rec == 0) continue;

        //confirm local state
        rel->confirm(level,j);

        if(j>=nb.getSize())
          {
          int new_var_bound = resF->isExtensibleLevel(nb.getLevel())? -(j+1): (j+1);
          dm->enlargeVariableBound(nb.getLevel(), false, new_var_bound);
          const unsigned oldSize = nb.getSize();
          nb.resize(j+1);
          nb.clear(oldSize, nb.getSize());
          // while(oldSize < nb.getSize()) { nb.d_ref(oldSize++) = 0; }
          queue->resize(nb.getSize());
          }

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

      } // for all events, ei

    } // No union possible
    else {
    std::unordered_map<long,std::vector<rel_node_handle>> list_of_j = rel->getListOfNexts(nb.getLevel(),i, Ru);

    for (std::unordered_map<long,std::vector<rel_node_handle>>::iterator jt=list_of_j.begin(); jt!=list_of_j.end(); ++jt) {
     //For each j get the list of different events
      int next = jt->first;

      int j = next;
      if(j==-1) continue;
      if (j < nb.getSize() && -1==nb.down(j)) continue; // nothing can be added to this set

      node_handle rec = recFireSet(nb.down(i), jt->second);
      if (rec == 0) continue;

      //confirm local state
      rel->confirm(level,j);

      if(j>=nb.getSize())
        {
        int new_var_bound = resF->isExtensibleLevel(nb.getLevel())? -(j+1): (j+1);
        dm->enlargeVariableBound(nb.getLevel(), false, new_var_bound);
        const unsigned oldSize = nb.getSize();
        nb.resize(j+1);
        nb.clear(oldSize, nb.getSize());
        // while(oldSize < nb.getSize()) { nb.d_ref(oldSize++) = 0; }
        queue->resize(nb.getSize());
        }

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

    }// for all j's

   } //if union possible

 }// more indexes to explore

  delete[] Ru;
  recycle(queue);
}

MEDDLY::node_handle MEDDLY::forwd_impl_dfs_by_events_mt::recFireSet(
                                                                    MEDDLY::node_handle mdd,
                                                                    std::vector<rel_node_handle> vector_mxd)
{
  std::vector<node_handle> array_rec(vector_mxd.size());

  dd_edge ans(resF), union_rec(resF);

  for(int rn = 0; rn < vector_mxd.size(); rn ++){
    ans.set( recFire(mdd,vector_mxd[rn]) );
    mddUnion->computeTemp(union_rec, ans, union_rec);
  }

  return resF->linkNode(union_rec.getNode());
}


// Same as post-image, except we saturate before reducing.
MEDDLY::node_handle MEDDLY::forwd_impl_dfs_by_events_mt::recFire(
    MEDDLY::node_handle mdd, rel_node_handle mxd)
{

  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;

  if (mxd==1) {
    if (argF->isTerminalNode(mdd)) {
      return resF->handleForValue(1);
    }
    // mxd is identity
    if (argF == resF)
      return resF->linkNode(mdd);
  }

  relation_node* relNode = rel->nodeExists(mxd); // The relation node

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
  const int mddLevel = argF->getNodeLevel(mdd);
  const int mxdLevel = relNode->getLevel();
  const int rLevel = MAX(mxdLevel, mddLevel);
   int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);
  domain* dm = resF->getDomain();

  dd_edge nbdj(resF), newst(resF);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(argF);
  if (mddLevel < rLevel) {
    A->initRedundant(argF, rLevel, mdd, FULL_ONLY);
  } else {
    argF->unpackNode(A, mdd, FULL_ONLY);
  }

  //Re-Think
  if (mddLevel > mxdLevel) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i=0; i<rSize; i++) {
      nb->setFull(i, recFire(A->down(i), mxd));
      // nb->d_ref(i) = recFire(A->down(i), mxd);
    }

  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(mxdLevel >= mddLevel);

    // Initialize mxd readers, note we might skip the unprimed level

    // loop over mxd "rows"
        for (int iz=0; iz<rSize; iz++) {
          int i = iz; // relation_node enabling condition
          if (0==A->down(i))   continue;

          // loop over mxd "columns"
          int j = relNode->nextOf(i);
          if(j==-1) continue;

          node_handle newstates = recFire(A->down(i), relNode->getDown());
          if (0==newstates) continue;

          //confirm local state
          if(!rel->isConfirmedState(rLevel,j)) // if not confirmed before
            {
            rel->confirm(rLevel,j); // confirm and enlarge
            if (j >= nb->getSize()) {
              int new_var_bound = resF->isExtensibleLevel(nb->getLevel())? -(j+1): (j+1);
              dm->enlargeVariableBound(nb->getLevel(), false, new_var_bound);
              const unsigned oldSize = nb->getSize();
              nb->resize(j+1);
              // while(oldSize < nb->getSize()) { nb->d_ref(oldSize++) = 0; }
              nb->clear(oldSize, nb->getSize());
            }
            }
          // ok, there is an i->j "edge".
          // determine new states to be added (recursively)
          // and add them

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

        } // for i


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

  //saveResult(Key, mdd, mxd, result);

  //node_handle resultTest = 0;
  // findResult(mdd,mxd,resultTest);



  return saveResult(Key, mdd, mxd, result);
}

// ******************************************************************
// *                                                                *
// *             common_impl_dfs_by_events_mt  methods              *
// *                                                                *
// ******************************************************************

MEDDLY::common_impl_dfs_by_events_mt::common_impl_dfs_by_events_mt(
    const char* name, implicit_relation* R)
    : saturation_operation(name, 1, R->getInForest(), R->getOutForest())
{
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;
  freeqs = 0;
  freebufs = 0;
  rel = R;

  ct_entry_type* et = new ct_entry_type(name, "NL:N");
  et->setForestForSlot(0, argF);
  et->setForestForSlot(3, resF);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::common_impl_dfs_by_events_mt::~common_impl_dfs_by_events_mt()
{
  delete rel;
}

bool MEDDLY::common_impl_dfs_by_events_mt
::isReachable(const dd_edge &a, const dd_edge& constraint)
{
  // Initialize operations
  if (0 == mddUnion) {
    mddUnion = UNION(resF, resF, resF);
    MEDDLY_DCASSERT(mddUnion);
  }

#ifdef DEBUG_INITIAL
  printf("Calling isReachable for states:\n");
  ostream_output s(std::cout);
  a.showGraph(s);
  std::cout.flush();
#endif

  // Execute saturation operation
  unary_list dummy("Saturation_by_events");
  saturation_impl_by_events_op* so = new saturation_impl_by_events_op(
          this, "Impl_sat", argF, resF);
  bool is_reachable = so->isReachable(a, constraint);

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
  return is_reachable;
}

void MEDDLY::common_impl_dfs_by_events_mt
::compute(const dd_edge &a, dd_edge &c)
{
  // Initialize operations
  mddUnion = UNION(resF, resF, resF);
  MEDDLY_DCASSERT(mddUnion);

  /*mxdIntersection = INTERSECTION(arg2F, arg2F, arg2F);
   MEDDLY_DCASSERT(mxdIntersection);

   mxdDifference = DIFFERENCE(arg2F, arg2F, arg2F);
   MEDDLY_DCASSERT(mxdDifference);*/

#ifdef DEBUG_INITIAL
  printf("Calling saturate for states:\n");
  ostream_output s(std::cout);
  a.showGraph(s);
  std::cout.flush();
#endif
#ifdef DEBUG_NSF
  printf("Calling saturate for NSF:\n");
  // b.showGraph(stdout);
#endif

  // Execute saturation operation
  /* if (!rel->isFinalized()) {
   printf("Transition relation has not been finalized.\n");
   printf("Finalizing using default options... ");
   rel->finalize();
   printf("done.\n");
   }*/

  saturation_impl_by_events_op* so = new saturation_impl_by_events_op(this,
          "Impl_sat", argF, resF);
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
// *       common_impl_dfs_by_events_mt::indexq  methods                 *
// ******************************************************************

MEDDLY::common_impl_dfs_by_events_mt::indexq::indexq()
{
  data = 0;
  size = 0;
  head = NULPTR;
}

MEDDLY::common_impl_dfs_by_events_mt::indexq::~indexq()
{
  free(data);
}

void MEDDLY::common_impl_dfs_by_events_mt::indexq::resize(int sz)
{
  if (sz <= size) return;
  data = (int*) realloc(data, sz * sizeof(int));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);

  for (; size < sz; size++) data[size] = NOTINQ;
}

// ******************************************************************
// *       common_impl_dfs_by_events_mt::charbuf methods                 *
// ******************************************************************

MEDDLY::common_impl_dfs_by_events_mt::charbuf::charbuf()
{
  data = 0;
  size = 0;
}

MEDDLY::common_impl_dfs_by_events_mt::charbuf::~charbuf()
{
  free(data);
}

void MEDDLY::common_impl_dfs_by_events_mt::charbuf::resize(int sz)
{
  if (sz <= size) return;
  data = (char*) realloc(data, sz * sizeof(char));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *               saturation_impl_by_events_op  methods            *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_impl_by_events_op
::saturation_impl_by_events_op(common_impl_dfs_by_events_mt* p,
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

MEDDLY::saturation_impl_by_events_op::~saturation_impl_by_events_op()
{
  // removeAllComputeTableEntries();
}


void MEDDLY::saturation_impl_by_events_op
::compute(const dd_edge& a, dd_edge& r)
{
  // Saturate
  const node_handle mdd = a.getNode();
#ifdef DEBUG_INITIAL
  printf("Calling saturate for states:\n");
  ostream_output s(std::cout);
  argF->showNodeGraph(s, &mdd, 1);
  std::cout.flush();
#endif
  r.set( saturate(mdd, argF->getMaxLevelIndex()) );
}

MEDDLY::node_handle
MEDDLY::saturation_impl_by_events_op::saturate(node_handle mdd, int k)
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

  const int sz = argF->getLevelSize(k);               // size
  const int mdd_level = argF->getNodeLevel(mdd);      // mdd level

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
  for (int i=0; i<sz; i++) {
    nb->setFull(i, mddDptrs->down(i) ? saturate(mddDptrs->down(i), k-1) : 0);
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


// ------------------
// Deadlock detection
// ------------------

struct node_pair {
  MEDDLY::node_handle first;
  MEDDLY::node_handle second;
};

bool operator<(const node_pair& l, const node_pair& r) {
  return (l.first < r.first || (l.first == r.first && l.second < r.second));
}

std::map<node_pair, bool> intersection_cache;
int time_since_gc = 0;

bool isIntersectionEmpty(
    MEDDLY::forest* mddF,
    MEDDLY::node_handle node_A,
    MEDDLY::node_handle node_B)
{
    using namespace MEDDLY;

  if (forest::isTerminalNode(node_A)) {
    // if (node_A == 0) return true; else return (node_B == 0);
    return (node_A == 0) || (node_B == 0);
  } else if (forest::isTerminalNode(node_B)) {
    // if (node_B == 0) return true; else return false;
    return (node_B == 0);
  }

  if (node_A > node_B) {
    // lexicographically ordered to improve cache hits
    // intersection is commutative
    node_handle temp = node_A;
    node_A = node_B;
    node_B = temp;
  }

  // search cache
  node_pair key = {node_A, node_B};
  auto entry_result = intersection_cache.find(key);
  if (entry_result != intersection_cache.end()) {
    // found cached entry
    return entry_result->second;
  }

  // unpack nodes
  unpacked_node* unp_A =
    mddF->getNodeLevel(node_A) >= mddF->getNodeLevel(node_B)
    ? mddF->newUnpacked(node_A, FULL_ONLY)
    : unpacked_node::newRedundant(mddF, mddF->getNodeLevel(node_B), node_A, FULL_ONLY);
  unpacked_node* unp_B =
    mddF->getNodeLevel(node_B) >= mddF->getNodeLevel(node_A)
    ? mddF->newUnpacked(node_B, FULL_ONLY)
    : unpacked_node::newRedundant(mddF, mddF->getNodeLevel(node_A), node_B, FULL_ONLY);

  // compute result
  bool result = true;
  int min_size = unp_A->getSize() < unp_B->getSize()? unp_A->getSize(): unp_B->getSize();
  for (int i = 0; i < min_size && result; i++) {
    if (!isIntersectionEmpty(mddF, unp_A->down(i), unp_B->down(i))) result = false;
  }

  // recycle unpacked nodes
  unpacked_node::Recycle(unp_A);
  unpacked_node::Recycle(unp_B);

  // cache result
  mddF->cacheNode(node_A);
  mddF->cacheNode(node_B);
  intersection_cache.emplace(std::make_pair(key, result));
  MEDDLY_DCASSERT(result == intersection_cache[key]);

  // garbage collection on intersection cache
  if (++time_since_gc > 1000000) {
    std::map<node_pair, bool> temp_intersection_cache;
    for (auto& i : intersection_cache) {
      const node_pair& p = i.first;
      if (!mddF->isActiveNode(p.first) || !mddF->isActiveNode(p.second)) {
        mddF->uncacheNode(p.first);
        mddF->uncacheNode(p.second);
      } else {
        temp_intersection_cache.emplace_hint(
            temp_intersection_cache.end(),
            std::make_pair(p, i.second));
      }
    }
    intersection_cache.swap(temp_intersection_cache);
    time_since_gc = 0;
  }

  return result;
}

bool MEDDLY::saturation_impl_by_events_op::isReachable(const dd_edge &a,
        const dd_edge &c)
{
  // Saturate and check is any element in constraint is reachable
  const node_handle mdd = a.getNode();
  const node_handle constraint = c.getNode();
  node_handle saturation_result = 0;
  bool result =
    isIntersectionEmpty(argF, mdd, constraint)
    ? isReachable(mdd, argF->getMaxLevelIndex(), constraint, saturation_result)
      // nothing reachable in initial state, continue exploring
    : true;
      // found reachable state in initial state

  // clear cache
  for (auto& i : intersection_cache) {
    argF->unlinkNode(i.first.first);
    argF->unlinkNode(i.first.second);
  }
  intersection_cache.clear();

#ifdef DEBUG_IS_REACHABLE
  if (false == result) {
    printf("isReachable return false, reachable states:\n");
    ostream_output s(std::cout);
    argF->showNodeGraph(s, &saturation_result, 1);
    std::cout.flush();
  }
#endif

  if (saturation_result) argF->unlinkNode(saturation_result);

  return result;
}

// Modified version of saturate() for detecting reachability.
// Recursively calls isReachable (similar to saturate() calling saturate()),
// and finally saturates the new node before returning result of saturation.
// Note that if isReachable returns false, then the information returned via
// saturation_result is the saturated node,
// and if isReachable returns true, then there are no guarantees for the information
// returned via saturated_result.
bool
MEDDLY::saturation_impl_by_events_op::isReachable(
  node_handle mdd,
  int k,
  node_handle constraint,         // set of states we are looking to reach
  node_handle& saturation_result)
{
#ifdef DEBUG_DFS
  printf("mdd: %d, k: %d, constraint: %d\n", mdd, k, constraint);
#endif

  // terminal conditions for recursion
  if (argF->isTerminalNode(mdd)) {
    saturation_result = mdd;
    return !isIntersectionEmpty(argF, saturation_result, constraint);
  }

  // If no states in constraint are reachable then return false
  if (argF->isTerminalNode(constraint)) {
    // if constraint is 0, then no state in constraint can be discovered, return false.
    // if constraint is 1, then there is at least one reachable state in constraint, return true.
    MEDDLY_DCASSERT(!argF->isTerminalNode(mdd));
    if (0 == constraint) {
      saturation_result = saturate(mdd, k);
      return false;
    }
    return true;
  }

  // search compute table
  node_handle n = 0;
  ct_entry_key* Key = findSaturateResult(mdd, k, n);
  if (0==Key) {
    saturation_result = n;
    return !isIntersectionEmpty(argF, saturation_result, constraint);
  }

  const int mdd_level = argF->getNodeLevel(mdd);
  const int constraint_level = argF->getNodeLevel(constraint);

  // Initialize mdd reader
  unpacked_node *mddDptrs = unpacked_node::New(argF);
  if (mdd_level < k) {
    mddDptrs->initRedundant(argF, k, mdd, FULL_ONLY);
  } else {
    argF->unpackNode(mddDptrs, mdd, FULL_ONLY);
  }
  MEDDLY_DCASSERT(!mddDptrs->isExtensible());

  // Initialize constraint reader
  unpacked_node* consDptrs = unpacked_node::New(argF);
  if (constraint_level < k) {
    consDptrs->initRedundant(argF, k, constraint, FULL_ONLY);
  } else {
    argF->unpackNode(consDptrs, constraint, FULL_ONLY);
  }

  // Initialize writer for result node, nb.
  const int sz = mddDptrs->getSize();
  unpacked_node* nb = unpacked_node::newFull(resF, k, sz);

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d, mdd_level: %d\n",
      mdd, k, sz, mdd_level);
#endif

#ifdef ALLOW_EXTENSIBLE
  const node_handle ext_d = consDptrs->isExtensible() ? consDptrs->ext_d() : 0;
#else
  const node_handle ext_d = 0;
#endif

  // Do computation
  for (int i=0; i<sz; i++) {
    node_handle temp = 0;
    node_handle cons_i = (i < consDptrs->getSize() ? consDptrs->down(i) : ext_d);
    if (isReachable(mddDptrs->down(i), k-1, cons_i, temp)) {
      // found reachable state in constraint, cleanup and return true
      for (int j = 0; j < i; j++) {
        if (nb->down(j)) {
            argF->unlinkNode(nb->down(j));
            nb->setFull(j, 0);
            // nb->d_ref(j) = 0;
        }
      }
      unpacked_node::Recycle(nb);
      unpacked_node::Recycle(mddDptrs);
      unpacked_node::Recycle(consDptrs);
      recycleCTKey(Key);
      return true;
    } else {
      nb->setFull(i, temp);
      // nb->d_ref(i) = temp;
    }
  }

  // Cleanup
  unpacked_node::Recycle(mddDptrs);
  unpacked_node::Recycle(consDptrs);

  // Reduce nb and save in compute table
  if (parent->saturateHelper(*nb, constraint)) {
    // found a reachable state in constraint
    for (int j = 0; j < sz; j++) {
      if (nb->down(j)) {
          argF->unlinkNode(nb->down(j));
          nb->setFull(j, 0);
         // nb->d_ref(j) = 0;
      }
    }
    unpacked_node::Recycle(nb);
    recycleCTKey(Key);
    return true;
  }

  n = resF->createReducedNode(-1, nb);
  saveSaturateResult(Key, mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  saturation_result = n;
  return !isIntersectionEmpty(argF, saturation_result, constraint);
}


bool MEDDLY::forwd_impl_dfs_by_events_mt::saturateHelper(
    unpacked_node& nb,
    node_handle constraint)
{

  if (0 == constraint) {
    // no reachable state in constraint possible
    saturateHelper(nb);
    return false;
  }

  int nEventsAtThisLevel = rel->lengthForLevel(nb.getLevel());
  if (0 == nEventsAtThisLevel) return false;

  const int constraint_level = argF->getNodeLevel(constraint);
  const int level = nb.getLevel();

  // Initialize constraint reader
  unpacked_node* consDptrs = unpacked_node::New(argF);
  if (constraint_level < level) {
    consDptrs->initRedundant(argF, level, constraint, FULL_ONLY);
  } else {
    argF->unpackNode(consDptrs, constraint, FULL_ONLY);
  }
#ifdef ALLOW_EXTENSIBLE
  const node_handle cons_ext_d = consDptrs->isExtensible() ? consDptrs->ext_d() : 0;
#else
  const node_handle cons_ext_d = 0;
#endif

  // Initialize mxd readers, note we might skip the unprimed level
  node_handle* events = rel->arrayForLevel(level);
  relation_node** Ru = new relation_node*[nEventsAtThisLevel];
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    Ru[ei] = rel->nodeExists(events[ei]);
#ifdef DEVELOPMENT_CODE
    int eventLevel = Ru[ei]->getLevel();
    MEDDLY_DCASSERT(ABS(eventLevel) == level);
#endif
  }

  dd_edge nbdj(resF), newst(resF);

  domain* dm = resF->getDomain();

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.down(i)) {
      queue->add(i);
    }
  }

  // explore indexes
  while (!queue->isEmpty()) {
    int i = queue->remove();

    MEDDLY_DCASSERT(nb.down(i));

    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {

      int j = Ru[ei]->nextOf(i);
      if(j==-1) continue;
      if (j < nb.getSize() && -1==nb.down(j)) continue; // nothing can be added to this set

      const node_handle cons_j = (j < consDptrs->getSize() ? consDptrs->down(j) : cons_ext_d);
      node_handle rec = 0;
      if (recFire(nb.down(i), Ru[ei]->getDown(), cons_j, rec)) {
        // found reachable state in constraint
        unpacked_node::Recycle(consDptrs);
        delete[] Ru;
        while (!queue->isEmpty()) queue->remove();
        recycle(queue);
        return true;
      }

      if (rec == 0) continue;

      //confirm local state
      rel->confirm(level,j);

      if(j>=nb.getSize())
      {
        int new_var_bound = resF->isExtensibleLevel(nb.getLevel())? -(j+1): (j+1);
        dm->enlargeVariableBound(nb.getLevel(), false, new_var_bound);
        const unsigned oldSize = nb.getSize();
        nb.resize(j+1);
        nb.clear(oldSize, nb.getSize());
        // while(oldSize < nb.getSize()) { nb.d_ref(oldSize++) = 0; }
        queue->resize(nb.getSize());
      }

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

    } // for all events, ei

  }// more indexes to explore

  delete[] Ru;
  recycle(queue);
  return false;
}


// Same as post-image, except we saturate before reducing.
bool MEDDLY::forwd_impl_dfs_by_events_mt::recFire(
     MEDDLY::node_handle mdd,
     rel_node_handle mxd,
     MEDDLY::node_handle constraint,
     MEDDLY::node_handle& result)
{

  if (0 == constraint) {
    // no reachable state in constraint possible
    result = recFire(mdd, mxd);
    return false;
  }

  // termination conditions
  if (mxd == 0 || mdd == 0) {
    result = 0;
    return false;
  }

  if (mxd==1) {
    if (argF->isTerminalNode(mdd) || argF == resF) {
      result =
        argF->isTerminalNode(mdd)
        ? resF->handleForValue(1)
        : resF->linkNode(mdd);
      return !isIntersectionEmpty(argF, result, constraint);
    }
  }

  relation_node* relNode = rel->nodeExists(mxd); // The relation node

  // check the cache
  result = 0;
  ct_entry_key* Key = findResult(mdd, mxd, result);
  if (0==Key) {
    return !isIntersectionEmpty(argF, result, constraint);
  }

#ifdef TRACE_RECFIRE
  printf("computing recFire(%d, %d)\n", mdd, mxd);
  printf("  node %3d ", mdd);
  argF->showNode(stdout, mdd, 1);
  printf("\n  node %3d ", mxd);
  arg2F->showNode(stdout, mxd, 1);
  printf("\n");
#endif

  // check if mxd and mdd are at the same level
  const int mddLevel = argF->getNodeLevel(mdd);
  const int mxdLevel = relNode->getLevel();
  const int rLevel = MAX(mxdLevel, mddLevel);
  const int constraint_level = argF->getNodeLevel(constraint);
  int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);
  domain* dm = resF->getDomain();

  dd_edge nbdj(resF), newst(resF);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(argF);
  if (mddLevel < rLevel) {
    A->initRedundant(argF, rLevel, mdd, FULL_ONLY);
  } else {
    argF->unpackNode(A, mdd, FULL_ONLY);
  }

  // Initialize constraint reader
  unpacked_node* consDptrs = unpacked_node::New(argF);
  if (constraint_level < rLevel) {
    consDptrs->initRedundant(argF, rLevel, constraint, FULL_ONLY);
  } else {
    argF->unpackNode(consDptrs, constraint, FULL_ONLY);
  }
#ifdef ALLOW_EXTENSIBLE
  const node_handle cons_ext_d = consDptrs->isExtensible() ? consDptrs->ext_d() : 0;
#else
  const node_handle cons_ext_d = 0;
#endif

  //Re-Think
  if (mddLevel > mxdLevel) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i=0; i<rSize; i++) {
      const node_handle cons_i = (i < consDptrs->getSize() ? consDptrs->down(i) : cons_ext_d);
      node_handle temp = 0;
      if (recFire(A->down(i), mxd, cons_i, temp)) {
        // found reachable state in constraint: abort, cleanup and return true
        for (int j = 0; j < i; j++) {
          if (nb->down(j)) {
              argF->unlinkNode(nb->down(j));
              nb->setFull(j, 0);
              //nb->d_ref(j) = 0;
          }
        }
        unpacked_node::Recycle(nb);
        unpacked_node::Recycle(A);
        unpacked_node::Recycle(consDptrs);
        return true;
      }
      nb->setFull(i, temp);
      // nb->d_ref(i) = temp;
    }

  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(mxdLevel >= mddLevel);

    // Initialize mxd readers, note we might skip the unprimed level

    // loop over mxd "rows"
    for (int iz=0; iz<rSize; iz++) {
      int i = iz; // relation_node enabling condition
      if (0==A->down(i))   continue;

      // loop over mxd "columns"
      int j = relNode->nextOf(i);
      if(j==-1) continue;

      node_handle newstates = 0;
      const node_handle cons_j = (j < consDptrs->getSize() ? consDptrs->down(j) : cons_ext_d);
      if (recFire(A->down(i), relNode->getDown(), cons_j, newstates)) {
        // found reachable state in constraint: abort, cleanup and return true
        for (int k = 0; k < nb->getSize(); k++) {
          if (nb->down(k)) {
              argF->unlinkNode(nb->down(k));
              nb->setFull(k, 0);
              // nb->d_ref(k) = 0;
          }
        }
        unpacked_node::Recycle(nb);
        unpacked_node::Recycle(A);
        unpacked_node::Recycle(consDptrs);
        return true;
      }

      if (0==newstates) continue;

      //confirm local state
      if(!rel->isConfirmedState(rLevel,j)) // if not confirmed before
      {
        rel->confirm(rLevel,j); // confirm and enlarge
        if (j >= nb->getSize()) {
          int new_var_bound = resF->isExtensibleLevel(nb->getLevel())? -(j+1): (j+1);
          dm->enlargeVariableBound(nb->getLevel(), false, new_var_bound);
          int oldSize = nb->getSize();
          nb->resize(j+1);
          nb->clear(oldSize, nb->getSize());
          // while(oldSize < nb->getSize()) { nb->d_ref(oldSize++) = 0; }
        }
      }
      // ok, there is an i->j "edge".
      // determine new states to be added (recursively)
      // and add them

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
    } // for i


  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);
  unpacked_node::Recycle(consDptrs);

  if (saturateHelper(*nb, constraint)) {
    // found reachable state in constraint
    const int sz = nb->getSize();
    for (int j = 0; j < sz; j++) {
      if (nb->down(j)) {
          argF->unlinkNode(nb->down(j));
          nb->setFull(j, 0);
          // nb->d_ref(j) = 0;
      }
    }
    unpacked_node::Recycle(nb);
    recycleCTKey(Key);
    return true;
  }

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

  saveResult(Key, mdd, mxd, result);
  return !isIntersectionEmpty(argF, result, constraint);
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************


MEDDLY::saturation_operation* MEDDLY::SATURATION_IMPL_FORWARD(forest* inF,
        implicit_relation* nsf, forest* outF)
{
    if (!nsf)
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    if (inF != nsf->getInForest())
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    if (outF != nsf->getOutForest())
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

    return new forwd_impl_dfs_by_events_mt("ImplSat", nsf);
}



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
#include "sat_impl.h"
#include <typeinfo> // for "bad_cast" exception
#include <set>
#include <map>

namespace MEDDLY {
  class saturation_impl_by_events_opname;
  class saturation_impl_by_events_op;
  
  class common_impl_dfs_by_events_mt;
  class forwd_impl_dfs_by_events_mt;
};

// #define DEBUG_INITIAL
// #define DEBUG_IS_REACHABLE

// ******************************************************************
// *                                                                *
// *               saturation_impl_by_events_opname  class               *
// *                                                                *
// ******************************************************************

/** Simple class to keep compute table happy.
 */
class MEDDLY::saturation_impl_by_events_opname : public unary_opname {
  static saturation_impl_by_events_opname* instance;
public:
  saturation_impl_by_events_opname();
  
  static const saturation_impl_by_events_opname* getInstance();
  
};

MEDDLY::saturation_impl_by_events_opname* MEDDLY::saturation_impl_by_events_opname::instance = 0;

MEDDLY::saturation_impl_by_events_opname::saturation_impl_by_events_opname()
: unary_opname("Saturate_by_events")
{
}

const MEDDLY::saturation_impl_by_events_opname* MEDDLY::saturation_impl_by_events_opname::getInstance()
{
  if (0==instance) instance = new saturation_impl_by_events_opname;
  return instance;
}

// ******************************************************************
// *                                                                *
// *             saturation_impl_by_events_op  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_impl_by_events_op : public unary_operation {
  common_impl_dfs_by_events_mt* parent;
public:
  saturation_impl_by_events_op(common_impl_dfs_by_events_mt* p,
                               expert_forest* argF, expert_forest* resF);
  virtual ~saturation_impl_by_events_op();
  
  node_handle saturate(node_handle mdd);
  node_handle saturate(node_handle mdd, int level);

  // for reachable state in constraint detection
  bool isReachable(
    node_handle mdd,
    node_handle constraint);
  bool isReachable(
    node_handle mdd,
    int k,
    node_handle constraint,
    node_handle& saturation_result);

  
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
  inline void recycleCTKey(compute_table::search_key* CTsrch) {
    doneCTkey(CTsrch);
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
// *            common_impl_dfs_by_events_mt  class                 *
// *                                                                *
// ******************************************************************

class MEDDLY::common_impl_dfs_by_events_mt : public specialized_operation {
public:
  common_impl_dfs_by_events_mt(const satimpl_opname* opcode,
                               satimpl_opname::implicit_relation* rel);
  virtual ~common_impl_dfs_by_events_mt();
  
#ifndef USE_NODE_STATUS
  virtual bool isStaleEntry(const node_handle* entryData);
#else
  virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle*);
#endif
  virtual void discardEntry(const node_handle* entryData);
  virtual void showEntry(output &strm, const node_handle* entryData) const;
  virtual void compute(const dd_edge& a, dd_edge &c);
  virtual bool isReachable(const dd_edge& a, const dd_edge& constraint);
  virtual void saturateHelper(unpacked_node& mdd) = 0;
  // for detecting reachable state in constraint
  virtual bool saturateHelper(unpacked_node& mdd, node_handle constraint) = 0;
  
protected:
  inline compute_table::search_key*
  findResult(node_handle a, rel_node_handle b, node_handle &c)
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
  inline void recycleCTKey(compute_table::search_key* CTsrch) {
    doneCTkey(CTsrch);
  }
  inline node_handle saveResult(compute_table::search_key* Key,
                                node_handle a, rel_node_handle b, node_handle c)
  {
  arg1F->cacheNode(a);
  compute_table::entry_builder &entry = CT->startNewEntry(Key);
  entry.writeResultNH(resF->cacheNode(c));
  CT->addEntry();
  return c;
  }
  
protected:
  binary_operation* mddUnion;
  binary_operation* mxdIntersection;
  binary_operation* mxdDifference;
  
  satimpl_opname::implicit_relation* rel;
  
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
// *               forwd_impl_dfs_by_events_mt class                *
// *                                                                *
// ******************************************************************


class MEDDLY::forwd_impl_dfs_by_events_mt : public common_impl_dfs_by_events_mt {
public:
  forwd_impl_dfs_by_events_mt(const satimpl_opname* opcode,
                              satimpl_opname::implicit_relation* rel);
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
                                                                 const satimpl_opname* opcode,
                                                                 satimpl_opname::implicit_relation* rel)
: common_impl_dfs_by_events_mt(opcode, rel)
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
  
  expert_domain* dm = static_cast<expert_domain*>(resF->useDomain());
  
  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.d(i)) {
      queue->add(i);
    }
  }
  
  // explore indexes
  while (!queue->isEmpty()) {
    int i = queue->remove();
    
    MEDDLY_DCASSERT(nb.d(i));
    
    bool is_union = rel->isUnionPossible(nb.getLevel(),i,Ru);
    
    if(!is_union)
      {
      for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
        
        int j = Ru[ei]->nextOf(i);
        if(j==-1) continue;
        if (j < nb.getSize() && -1==nb.d(j)) continue; // nothing can be added to this set
        
        node_handle rec = recFire(nb.d(i), Ru[ei]->getDown());
        
        if (rec == 0) continue;
        
        //confirm local state
        rel->setConfirmedStates(level,j);
        
        if(j>=nb.getSize())
          {
          int new_var_bound = resF->isExtensibleLevel(nb.getLevel())? -(j+1): (j+1);
          dm->enlargeVariableBound(nb.getLevel(), false, new_var_bound);
          int oldSize = nb.getSize();
          nb.resize(j+1);
          while(oldSize < nb.getSize()) { nb.d_ref(oldSize++) = 0; }
          queue->resize(nb.getSize());
          }
        
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
        
      } // for all events, ei
      
    } // No union possible
    else { 
    std::unordered_map<long,std::vector<rel_node_handle>> list_of_j = rel->getListOfNexts(nb.getLevel(),i, Ru);
     
    for (std::unordered_map<long,std::vector<rel_node_handle>>::iterator jt=list_of_j.begin(); jt!=list_of_j.end(); ++jt) {
     //For each j get the list of different events
      int next = jt->first;
      
      int j = next;
      if(j==-1) continue;
      if (j < nb.getSize() && -1==nb.d(j)) continue; // nothing can be added to this set
        
      node_handle rec = recFireSet(nb.d(i), jt->second);
      if (rec == 0) continue;
        
      //confirm local state
      rel->setConfirmedStates(level,j);
        
      if(j>=nb.getSize())
        {
        int new_var_bound = resF->isExtensibleLevel(nb.getLevel())? -(j+1): (j+1);
        dm->enlargeVariableBound(nb.getLevel(), false, new_var_bound);
        int oldSize = nb.getSize();
        nb.resize(j+1);
        while(oldSize < nb.getSize()) { nb.d_ref(oldSize++) = 0; }
        queue->resize(nb.getSize());
        }
      
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
  node_handle union_rec = 0;
  
  for(int rn = 0; rn < vector_mxd.size(); rn ++){
    int ans=recFire(mdd,vector_mxd[rn]);
    union_rec = mddUnion->compute(union_rec,ans);
  }
  
  return union_rec;
}


// Same as post-image, except we saturate before reducing.
MEDDLY::node_handle MEDDLY::forwd_impl_dfs_by_events_mt::recFire(
                                                                 MEDDLY::node_handle mdd, rel_node_handle mxd)
{
  
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  
  if (mxd==1) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->handleForValue(1);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }
  
  relation_node* relNode = rel->nodeExists(mxd); // The relation node
  
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
  const int mxdLevel = relNode->getLevel();
  const int rLevel = MAX(mxdLevel, mddLevel);
   int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);
  expert_domain* dm = static_cast<expert_domain*>(resF->useDomain());
  
  // Initialize mdd reader
  unpacked_node *A = unpacked_node::useUnpackedNode();
  if (mddLevel < rLevel) {
    A->initRedundant(arg1F, rLevel, mdd, true);
  } else {
    A->initFromNode(arg1F, mdd, true);
  }
  
  //Re-Think
  if (mddLevel > mxdLevel) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    
    for (int i=0; i<rSize; i++) {
      nb->d_ref(i) = recFire(A->d(i), mxd);
    }
    
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(mxdLevel >= mddLevel);
    // clear out result (important!)
    for (int i=0; i<rSize; i++) {
      nb->d_ref(i) = 0;
    }
    
    // Initialize mxd readers, note we might skip the unprimed level
    
    // loop over mxd "rows"
        for (int iz=0; iz<rSize; iz++) {
          int i = iz; // relation_node enabling condition
          if (0==A->d(i))   continue;
          
          // loop over mxd "columns"
          int j = relNode->nextOf(i);
          if(j==-1) continue;
          
          node_handle newstates = recFire(A->d(i), relNode->getDown());
          if (0==newstates) continue;
          
          //confirm local state
          if(!rel->isConfirmedState(rLevel,j)) // if not confirmed before
            {
            rel->setConfirmedStates(rLevel,j); // confirm and enlarge
            if (j >= nb->getSize()) {
              int new_var_bound = resF->isExtensibleLevel(nb->getLevel())? -(j+1): (j+1);
              dm->enlargeVariableBound(nb->getLevel(), false, new_var_bound);
              int oldSize = nb->getSize();
              nb->resize(j+1);
              while(oldSize < nb->getSize()) { nb->d_ref(oldSize++) = 0; }
            }
            }
          // ok, there is an i->j "edge".
          // determine new states to be added (recursively)
          // and add them
          
          if (0==nb->d(j)) {
            nb->d_ref(j) = newstates;
            continue;
          }
          // there's new states and existing states; union them.
          int oldj = nb->d(j);
          nb->d_ref(j) = mddUnion->compute(newstates, oldj);
          resF->unlinkNode(oldj);
          resF->unlinkNode(newstates);
          
          
        } // for i

    
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
                                                                   const satimpl_opname* opcode,
                                                                   satimpl_opname::implicit_relation* relation)
: specialized_operation(opcode, 2, 1)
{
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;
  freeqs = 0;
  freebufs = 0;
  rel = relation;
  arg1F = static_cast<expert_forest*>(rel->getInForest());
  //arg2F = static_cast<expert_forest*>(rel->getInForest());
  resF = static_cast<expert_forest*>(rel->getOutForest());
  
  registerInForest(arg1F);
  //registerInForest(arg2F);
  registerInForest(resF);
  setAnswerForest(resF);
}

MEDDLY::common_impl_dfs_by_events_mt::~common_impl_dfs_by_events_mt()
{
  if (rel->autoDestroy()) delete rel;
  unregisterInForest(arg1F);
  //unregisterInForest(arg2F);
  unregisterInForest(resF);
}

#ifndef USE_NODE_STATUS
bool MEDDLY::common_impl_dfs_by_events_mt::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[0]) ||
  //arg2F->isStale(data[1]) ||
  resF->isStale(data[2]);
}
#else
MEDDLY::forest::node_status
MEDDLY::common_impl_dfs_by_events_mt::getStatusOfEntry(const node_handle* data)
{
  MEDDLY::forest::node_status a = arg1F->getNodeStatus(data[0]);
  MEDDLY::forest::node_status c = resF->getNodeStatus(data[2]);

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

void MEDDLY::common_impl_dfs_by_events_mt::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[0]);
  //arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void MEDDLY::common_impl_dfs_by_events_mt::showEntry(output &strm,
                                                     const node_handle* data) const
{
  strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) << "): " << long(data[2]) << "]";
}

bool MEDDLY::common_impl_dfs_by_events_mt
::isReachable(const dd_edge &a, const dd_edge& constraint)
{
  // Initialize operations
  if (0 == mddUnion) {
    mddUnion = getOperation(UNION, resF, resF, resF);
    MEDDLY_DCASSERT(mddUnion);
  }
  
#ifdef DEBUG_INITIAL
  printf("Calling isReachable for states:\n");
  ostream_output s(std::cout);
  a.show(s, 2);
  std::cout.flush();
#endif
  
  // Execute saturation operation
  saturation_impl_by_events_op* so = new saturation_impl_by_events_op(this, arg1F, resF);
  bool is_reachable = so->isReachable(a.getNode(), constraint.getNode());

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
  mddUnion = getOperation(UNION, resF, resF, resF);
  MEDDLY_DCASSERT(mddUnion);
  
  /*mxdIntersection = getOperation(INTERSECTION, arg2F, arg2F, arg2F);
   MEDDLY_DCASSERT(mxdIntersection);
   
   mxdDifference = getOperation(DIFFERENCE, arg2F, arg2F, arg2F);
   MEDDLY_DCASSERT(mxdDifference);*/
  
#ifdef DEBUG_INITIAL
  printf("Calling saturate for states:\n");
  ostream_output s(std::cout);
  a.show(s, 2);
  std::cout.flush();
#endif
#ifdef DEBUG_NSF
  printf("Calling saturate for NSF:\n");
  // b.show(stdout, 2);
#endif
  
  // Execute saturation operation
  /* if (!rel->isFinalized()) {
   printf("Transition relation has not been finalized.\n");
   printf("Finalizing using default options... ");
   rel->finalize();
   printf("done.\n");
   }*/
  
  saturation_impl_by_events_op* so = new saturation_impl_by_events_op(this, arg1F, resF);
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
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::satimpl_opname* MEDDLY::initImplSaturationForward()
{
  return new satimpl_opname("SaturationFwd");
}


MEDDLY::specialized_operation*
MEDDLY::satimpl_opname::buildOperation(arguments* a) const
{
  
  implicit_relation* rel = dynamic_cast<implicit_relation*>(a);
  if (0==rel) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  
  MEDDLY::specialized_operation* op = 0;
  op = new forwd_impl_dfs_by_events_mt(this, rel);
  
  return op;
}

// ******************************************************************
// *                                                                *
// *               saturation_impl_by_events_op  methods            *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_impl_by_events_op
::saturation_impl_by_events_op(common_impl_dfs_by_events_mt* p,
                               expert_forest* argF, expert_forest* resF)
: unary_operation(saturation_impl_by_events_opname::getInstance(),
                  ((argF != 0 && argF->isFullyReduced())? 2: 1), 1, argF, resF)
{
  parent = p;
  
}

MEDDLY::saturation_impl_by_events_op::~saturation_impl_by_events_op()
{
  removeAllComputeTableEntries();
}


MEDDLY::node_handle MEDDLY::saturation_impl_by_events_op::saturate(MEDDLY::node_handle mdd)
{
  // Saturate
#ifdef DEBUG_INITIAL
  printf("Calling saturate for states:\n");
  ostream_output s(std::cout);
  argF->showNodeGraph(s, &mdd, 1);
  std::cout.flush();
#endif
  return saturate(mdd, argF->getNumVariables());
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
  compute_table::search_key* Key = findSaturateResult(mdd, k, n);
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

#ifndef USE_NODE_STATUS
bool MEDDLY::saturation_impl_by_events_op::isStaleEntry(const node_handle* data)
{
  return (argF->isFullyReduced()
          ? (argF->isStale(data[0]) || resF->isStale(data[2]))
          : (argF->isStale(data[0]) || resF->isStale(data[1])));
}
#else
MEDDLY::forest::node_status
MEDDLY::saturation_impl_by_events_op::getStatusOfEntry(const node_handle* data)
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

void MEDDLY::saturation_impl_by_events_op::discardEntry(const node_handle* data)
{
  if (argF->isFullyReduced()) {
    argF->uncacheNode(data[0]);
    resF->uncacheNode(data[2]);
  } else {
    argF->uncacheNode(data[0]);
    resF->uncacheNode(data[1]);
  }
}

void MEDDLY::saturation_impl_by_events_op::showEntry(output &strm,
                                                     const node_handle* data) const
{
  if (argF->isFullyReduced()) {
    strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) << "): " << long(data[2]) << "]";
  } else {
    strm << "[" << getName() << "(" << long(data[0]) << "): " << long(data[1]) << "]";
  }
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
    MEDDLY::expert_forest* mddF,
    MEDDLY::node_handle node_A,
    MEDDLY::node_handle node_B)
{
  if (MEDDLY::expert_forest::isTerminalNode(node_A)) {
    // if (node_A == 0) return true; else return (node_B == 0);
    return (node_A == 0) || (node_B == 0);
  } else if (MEDDLY::expert_forest::isTerminalNode(node_B)) {
    // if (node_B == 0) return true; else return false;
    return (node_B == 0);
  }

  if (node_A > node_B) {
    // lexicographically ordered to improve cache hits
    // intersection is commutative
    MEDDLY::node_handle temp = node_A;
    node_A = node_B;
    node_B = temp;
  }

  // search cache
  node_pair key = {node_A, node_B};
  auto search_result = intersection_cache.find(key);
  if (search_result != intersection_cache.end()) {
    // found cached entry
    return search_result->second;
  }

  // unpack nodes
  MEDDLY::unpacked_node* unp_A =
    mddF->getNodeLevel(node_A) >= mddF->getNodeLevel(node_B)
    ? MEDDLY::unpacked_node::newFromNode(mddF, node_A, true)
    : MEDDLY::unpacked_node::newRedundant(mddF, mddF->getNodeLevel(node_B), node_A, true);
  MEDDLY::unpacked_node* unp_B =
    mddF->getNodeLevel(node_B) >= mddF->getNodeLevel(node_A)
    ? MEDDLY::unpacked_node::newFromNode(mddF, node_B, true)
    : MEDDLY::unpacked_node::newRedundant(mddF, mddF->getNodeLevel(node_A), node_B, true);

  // compute result
  bool result = true;
  int min_size = unp_A->getSize() < unp_B->getSize()? unp_A->getSize(): unp_B->getSize();
  for (int i = 0; i < min_size && result; i++) {
    if (!isIntersectionEmpty(mddF, unp_A->d(i), unp_B->d(i))) result = false;
  }

  // recycle unpacked nodes
  MEDDLY::unpacked_node::recycle(unp_A);
  MEDDLY::unpacked_node::recycle(unp_B);

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

bool MEDDLY::saturation_impl_by_events_op::isReachable(MEDDLY::node_handle mdd, MEDDLY::node_handle constraint)
{
  // Saturate and check is any element in constraint is reachable
  MEDDLY::node_handle saturation_result = 0;
  bool result = 
    isIntersectionEmpty(argF, mdd, constraint)
    ? isReachable(mdd, argF->getNumVariables(), constraint, saturation_result)
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
  compute_table::search_key* Key = findSaturateResult(mdd, k, n);
  if (0==Key) {
    saturation_result = n;
    return !isIntersectionEmpty(argF, saturation_result, constraint);
  }

  const int mdd_level = argF->getNodeLevel(mdd);
  const int constraint_level = argF->getNodeLevel(constraint);

  // Initialize mdd reader
  unpacked_node *mddDptrs = unpacked_node::useUnpackedNode();
  if (mdd_level < k) {
    mddDptrs->initRedundant(argF, k, mdd, true);
  } else {
    mddDptrs->initFromNode(argF, mdd, true);
  }
  MEDDLY_DCASSERT(!mddDptrs->isExtensible());

  // Initialize constraint reader
  unpacked_node* consDptrs = unpacked_node::useUnpackedNode();
  if (constraint_level < k) {
    consDptrs->initRedundant(argF, k, constraint, true);
  } else {
    consDptrs->initFromNode(argF, constraint, true);
  }

  // Initialize writer for result node, nb.
  const int sz = mddDptrs->getSize();
  unpacked_node* nb = unpacked_node::newFull(resF, k, sz);

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d, mdd_level: %d\n",
      mdd, k, sz, mdd_level);
#endif

  const node_handle ext_d = consDptrs->isExtensible() ? consDptrs->ext_d() : 0;

  // Do computation
  for (int i=0; i<sz; i++) {
    node_handle temp = 0;
    node_handle cons_i = (i < consDptrs->getSize() ? consDptrs->d(i) : ext_d);
    if (isReachable(mddDptrs->d(i), k-1, cons_i, temp)) {
      // found reachable state in constraint, cleanup and return true
      for (int j = 0; j < i; j++) {
        if (nb->d(j)) { argF->unlinkNode(nb->d(j)); nb->d_ref(j) = 0; }
      }
      unpacked_node::recycle(nb);
      unpacked_node::recycle(mddDptrs);
      unpacked_node::recycle(consDptrs);
      recycleCTKey(Key);
      return true;
    } else {
      nb->d_ref(i) = temp;
    }
  }

  // Cleanup
  unpacked_node::recycle(mddDptrs);
  unpacked_node::recycle(consDptrs);

  // Reduce nb and save in compute table
  if (parent->saturateHelper(*nb, constraint)) {
    // found a reachable state in constraint
    for (int j = 0; j < sz; j++) {
      if (nb->d(j)) { argF->unlinkNode(nb->d(j)); nb->d_ref(j) = 0; }
    }
    unpacked_node::recycle(nb);
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

  const int constraint_level = arg1F->getNodeLevel(constraint);
  const int level = nb.getLevel();

  // Initialize constraint reader
  unpacked_node* consDptrs = unpacked_node::useUnpackedNode();
  if (constraint_level < level) {
    consDptrs->initRedundant(arg1F, level, constraint, true);
  } else {
    consDptrs->initFromNode(arg1F, constraint, true);
  }
  const node_handle cons_ext_d = consDptrs->isExtensible() ? consDptrs->ext_d() : 0;

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

  expert_domain* dm = static_cast<expert_domain*>(resF->useDomain());

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.d(i)) {
      queue->add(i);
    }
  }

  // explore indexes
  while (!queue->isEmpty()) {
    int i = queue->remove();

    MEDDLY_DCASSERT(nb.d(i));

    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {

      int j = Ru[ei]->nextOf(i);
      if(j==-1) continue;
      if (j < nb.getSize() && -1==nb.d(j)) continue; // nothing can be added to this set

      const node_handle cons_j = (j < consDptrs->getSize() ? consDptrs->d(j) : cons_ext_d);
      node_handle rec = 0;
      if (recFire(nb.d(i), Ru[ei]->getDown(), cons_j, rec)) {
        // found reachable state in constraint
        unpacked_node::recycle(consDptrs);
        delete[] Ru;
        while (!queue->isEmpty()) queue->remove();
        recycle(queue);
        return true;
      }

      if (rec == 0) continue;

      //confirm local state
      rel->setConfirmedStates(level,j);

      if(j>=nb.getSize())
      {
        int new_var_bound = resF->isExtensibleLevel(nb.getLevel())? -(j+1): (j+1);
        dm->enlargeVariableBound(nb.getLevel(), false, new_var_bound);
        int oldSize = nb.getSize();
        nb.resize(j+1);
        while(oldSize < nb.getSize()) { nb.d_ref(oldSize++) = 0; }
        queue->resize(nb.getSize());
      }

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
    if (arg1F->isTerminalNode(mdd) || arg1F == resF) {
      result = 
        arg1F->isTerminalNode(mdd)
        ? resF->handleForValue(1)
        : resF->linkNode(mdd);
      return !isIntersectionEmpty(arg1F, result, constraint);
    }
  }

  relation_node* relNode = rel->nodeExists(mxd); // The relation node

  // check the cache
  result = 0;
  compute_table::search_key* Key = findResult(mdd, mxd, result);
  if (0==Key) {
    return !isIntersectionEmpty(arg1F, result, constraint);
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
  const int mddLevel = arg1F->getNodeLevel(mdd);
  const int mxdLevel = relNode->getLevel();
  const int rLevel = MAX(mxdLevel, mddLevel);
  const int constraint_level = arg1F->getNodeLevel(constraint);
  int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);
  expert_domain* dm = static_cast<expert_domain*>(resF->useDomain());

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::useUnpackedNode();
  if (mddLevel < rLevel) {
    A->initRedundant(arg1F, rLevel, mdd, true);
  } else {
    A->initFromNode(arg1F, mdd, true);
  }

  // Initialize constraint reader
  unpacked_node* consDptrs = unpacked_node::useUnpackedNode();
  if (constraint_level < rLevel) {
    consDptrs->initRedundant(arg1F, rLevel, constraint, true);
  } else {
    consDptrs->initFromNode(arg1F, constraint, true);
  }
  const node_handle cons_ext_d = consDptrs->isExtensible() ? consDptrs->ext_d() : 0;

  //Re-Think
  if (mddLevel > mxdLevel) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i=0; i<rSize; i++) {
      const node_handle cons_i = (i < consDptrs->getSize() ? consDptrs->d(i) : cons_ext_d);
      node_handle temp = 0;
      if (recFire(A->d(i), mxd, cons_i, temp)) {
        // found reachable state in constraint: abort, cleanup and return true
        for (int j = 0; j < i; j++) {
          if (nb->d(j)) { arg1F->unlinkNode(nb->d(j)); nb->d_ref(j) = 0; }
        }
        unpacked_node::recycle(nb);
        unpacked_node::recycle(A);
        unpacked_node::recycle(consDptrs);
        return true;
      }
      nb->d_ref(i) = temp;
    }

  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(mxdLevel >= mddLevel);
    // clear out result (important!)
    for (int i=0; i<rSize; i++) {
      nb->d_ref(i) = 0;
    }

    // Initialize mxd readers, note we might skip the unprimed level

    // loop over mxd "rows"
    for (int iz=0; iz<rSize; iz++) {
      int i = iz; // relation_node enabling condition
      if (0==A->d(i))   continue;

      // loop over mxd "columns"
      int j = relNode->nextOf(i);
      if(j==-1) continue;

      node_handle newstates = 0;
      const node_handle cons_j = (j < consDptrs->getSize() ? consDptrs->d(j) : cons_ext_d);
      if (recFire(A->d(i), relNode->getDown(), cons_j, newstates)) {
        // found reachable state in constraint: abort, cleanup and return true
        for (int k = 0; k < nb->getSize(); k++) {
          if (nb->d(k)) { arg1F->unlinkNode(nb->d(k)); nb->d_ref(k) = 0; }
        }
        unpacked_node::recycle(nb);
        unpacked_node::recycle(A);
        unpacked_node::recycle(consDptrs);
        return true;
      }

      if (0==newstates) continue;

      //confirm local state
      if(!rel->isConfirmedState(rLevel,j)) // if not confirmed before
      {
        rel->setConfirmedStates(rLevel,j); // confirm and enlarge
        if (j >= nb->getSize()) {
          int new_var_bound = resF->isExtensibleLevel(nb->getLevel())? -(j+1): (j+1);
          dm->enlargeVariableBound(nb->getLevel(), false, new_var_bound);
          int oldSize = nb->getSize();
          nb->resize(j+1);
          while(oldSize < nb->getSize()) { nb->d_ref(oldSize++) = 0; }
        }
      }
      // ok, there is an i->j "edge".
      // determine new states to be added (recursively)
      // and add them

      if (0==nb->d(j)) {
        nb->d_ref(j) = newstates;
        continue;
      }
      // there's new states and existing states; union them.
      int oldj = nb->d(j);
      nb->d_ref(j) = mddUnion->compute(newstates, oldj);
      resF->unlinkNode(oldj);
      resF->unlinkNode(newstates);


    } // for i


  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);
  unpacked_node::recycle(consDptrs);

  if (saturateHelper(*nb, constraint)) {
    // found reachable state in constraint
    const int sz = nb->getSize();
    for (int j = 0; j < sz; j++) {
      if (nb->d(j)) { arg1F->unlinkNode(nb->d(j)); nb->d_ref(j) = 0; }
    }
    unpacked_node::recycle(nb);
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
  return !isIntersectionEmpty(arg1F, result, constraint);
}


bool
MEDDLY::satimpl_opname::implicit_relation::isReachable(const dd_edge& initial_states, const dd_edge& constraint)
{
  // build implicit saturation operation operation
  specialized_operation* satop = SATURATION_IMPL_FORWARD->buildOperation(this);
  MEDDLY_DCASSERT(satop);
  forwd_impl_dfs_by_events_mt* op = dynamic_cast<forwd_impl_dfs_by_events_mt*>(satop);
  MEDDLY_DCASSERT(op);
#ifdef DEBUG_IS_REACHABLE
  std::cout << "[In " << __func__ << "] constraint dd_edge:\n";
  ostream_output s(std::cout);
  constraint.show(s, 2);
  std::cout.flush();
#endif
  return op->isReachable(initial_states, constraint);
}


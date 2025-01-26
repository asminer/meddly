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
//#include "event_ordering.h"
#include "sat_hyb.h"
#include <typeinfo> // for "bad_cast" exception
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "../operators.h"
#include "../minterms.h"
#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../oper_satur.h"
#include "../ops_builtin.h"
#include "../sat_relations.h"

   #define OUT_OF_BOUNDS -1
   #define NOT_KNOWN -2
   #define TERMINAL_NODE 1
  // #define CONJUNCT_SUBEVENTS
  // #define ENABLE_CONSTRAINED_SAT
  // #define NONREDUNDANT
   #define MIX_RELATION



namespace MEDDLY {
  class saturation_hyb_by_events_opname;
  class saturation_hyb_by_events_op;

  class common_hyb_dfs_by_events_mt;
  class forwd_hyb_dfs_by_events_mt;

};

// #define DEBUG_INITIAL
// #define DEBUG_IS_REACHABLE


// ******************************************************************
// *                                                                *
// *             saturation_hyb_by_events_op  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_hyb_by_events_op : public saturation_operation {
  common_hyb_dfs_by_events_mt* parent;
public:
  saturation_hyb_by_events_op(common_hyb_dfs_by_events_mt* p,
          const char* name, forest* argF, forest* resF);
  virtual ~saturation_hyb_by_events_op();

  virtual void compute(const dd_edge &arg, dd_edge &res);
  node_handle saturate(node_handle mdd);
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
// *            common_hyb_dfs_by_events_mt  class                 *
// *                                                                *
// ******************************************************************

class MEDDLY::common_hyb_dfs_by_events_mt : public saturation_operation {
public:
  common_hyb_dfs_by_events_mt(const char* name, forest* argF,
          hybrid_relation* rel, forest* resF);
  virtual ~common_hyb_dfs_by_events_mt();

  virtual void compute(const dd_edge& a, dd_edge &c);
  virtual void saturateHelper(unpacked_node& mdd) = 0;


protected:
  inline ct_entry_key*
  findResult(node_handle a, node_handle* b, int num_se, node_handle &c)
  {
    ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], num_se);
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeN(a);
    for(int i = 0; i<num_se; i++)
      CTsrch->writeN(b[i]);
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
                                node_handle a, node_handle* b, int num_se, node_handle c)
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

  hybrid_relation* rel;

  inline forest* relF() const { return rel->getHybridForest(); }

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
// *               forwd_hyb_dfs_by_events_mt class                *
// *                                                                *
// ******************************************************************


class MEDDLY::forwd_hyb_dfs_by_events_mt : public common_hyb_dfs_by_events_mt {
public:
  forwd_hyb_dfs_by_events_mt(const char* name, forest* argF,
          hybrid_relation* rel, forest* resF);
protected:
  virtual void saturateHelper(unpacked_node& nb);
 MEDDLY::node_handle recFire(MEDDLY::node_handle mdd, MEDDLY::node_handle* seHandles, int num_se);
  void recFireHelper(const unsigned, const int, const MEDDLY::node_handle, const MEDDLY::node_handle,
        unpacked_node*, unpacked_node*, MEDDLY::node_handle*, bool, int, int, int);

  //MEDDLY::node_handle recFireSet(node_handle mdd, std::set<node_handle> mxd);
  std::pair<MEDDLY::node_handle,int> getHighestNodeHandles(MEDDLY::node_handle* seHandles, int num_se);


 };


MEDDLY::forwd_hyb_dfs_by_events_mt::forwd_hyb_dfs_by_events_mt(
        const char* name, forest* argF, hybrid_relation* rel, forest* resF)
: common_hyb_dfs_by_events_mt(name, argF, rel, resF)
{
}


std::pair<MEDDLY::node_handle,int>
MEDDLY::forwd_hyb_dfs_by_events_mt::getHighestNodeHandles(MEDDLY::node_handle* seHandles, int num_se) {

    if(num_se == 1)
      return  std::make_pair(seHandles[0],0);

    int max_level = 0;
    int highest_nh = 0;
    int which_se = -1;

    for(int it = 0; it < num_se; it++)
    {
	    int p_lvl = relF()->getNodeLevel(seHandles[it]);
      if(p_lvl<0) p_lvl = -p_lvl;
    	if(p_lvl>=max_level)
        {
          max_level = p_lvl;
          highest_nh = seHandles[it];
          which_se = it;
        }
    }

      return std::make_pair(highest_nh,which_se);  // Not implemented the general case.
}

void MEDDLY::forwd_hyb_dfs_by_events_mt::saturateHelper(unpacked_node& nb)
{
  int nEventsAtThisLevel = rel->lengthForLevel(nb.getLevel());
  if (0 == nEventsAtThisLevel) return;

  // Initialize mxd readers, note we might skip the unprimed level
  const int level = nb.getLevel();
  hybrid_event** event_handles = rel->arrayForLevel(level);

  dd_edge nbdj(resF), newst(resF);

  domain* dm = resF->getDomain();

  int* event_Ru_Rn_index = (int*)malloc(nEventsAtThisLevel*sizeof(int));
  unpacked_node** Ru = new unpacked_node*[nEventsAtThisLevel];
  relation_node** Rn = new relation_node*[nEventsAtThisLevel];
  unpacked_node* Rp = unpacked_node::New(relF(), SPARSE_ONLY);
  int i_Ru = 0;
  int i_Rn = 0;

  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    event_handles[ei]->rebuild();
    //std::set<node_handle> seHandles = event_handles[ei]->getComponentAt(level);
    node_handle se_nh = event_handles[ei]->getTopComponent();
    if(relF()->isImplicit(se_nh)) {
      Rn[i_Rn] = relF()->buildImplicitNode(se_nh);
      event_Ru_Rn_index[ei] = i_Rn;
      i_Rn++;
    } else {
      Ru[i_Ru] = unpacked_node::New(relF(), SPARSE_ONLY);
      if(relF()->getNodeLevel(se_nh) == level)
          relF()->unpackNode(Ru[i_Ru], se_nh, FULL_ONLY);  // node is present at unprime-level
      else
         Ru[i_Ru]->initRedundant(relF(), level, se_nh, SPARSE_ONLY);
      event_Ru_Rn_index[ei] = i_Ru;
      i_Ru++;
    }
  }

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

      bool is_rebuilt = false;
      if(event_handles[ei]->getNumOfSubevents()>0)
       is_rebuilt = event_handles[ei]->rebuild();

      //std::set<node_handle> seHandles = event_handles[ei]->getComponentAt(level); // Obtain only one handle as per the design chosen
      // MEDDLY_DCASSERT(seHandles.size() == 1);
      int does_rn_exist = event_handles[ei]->getNumOfRelnodes()>0?1:0;
      int which_se = 0;
      node_handle se_nh = event_handles[ei]->getTopComponent();//*(seHandles.begin());
      // Based on whether it is relation-node or mxd node
      // Find the next of "i"
      // Find the next down handle
      // Prepare handles to be taken care of during recFire


      bool imFlag = relF()->isImplicit(se_nh);
      int jC;
      if(imFlag)
      {
      //  int j = Rn[event_Ru_Rn_index[ei]]->nextOf(i);
        jC = 1;

      } else {
        for(int w = 0; w<event_handles[ei]->getNumOfSubevents(); w++) {
        hybrid_subevent* my_se = event_handles[ei]->getSubevents()[w];
        if(se_nh == my_se->getRootHandle())
          {
           which_se = w; break;
          }
        }

        if(is_rebuilt)
        {
        if(relF()->getNodeLevel(se_nh) == level)
          relF()->unpackNode(Ru[event_Ru_Rn_index[ei]], se_nh, FULL_ONLY);  // node is present at unprime-level
        else
          Ru[event_Ru_Rn_index[ei]]->initRedundant(relF(), level, se_nh, SPARSE_ONLY);  // node was at prime-level, so build redudant node at unprime-level
        }

        node_handle ei_i_p = Ru[event_Ru_Rn_index[ei]]->down(i);
        if( 0 == ei_i_p) continue;

        const int dlevel = relF()->getNodeLevel(ei_i_p);

        if (dlevel == -level)
          relF()->unpackNode(Rp, ei_i_p, SPARSE_ONLY);
        else
          Rp->initIdentity(relF(), -level, i, ei_i_p, SPARSE_ONLY);

        jC = Rp->getSize();

      }
      int num_se = event_handles[ei]->getNumOfSubevents()+(event_handles[ei]->getNumOfRelnodes()>0?1:0);
      node_handle* seHandlesForRecursion = (node_handle*)malloc(num_se*sizeof(node_handle));
      memcpy(seHandlesForRecursion,event_handles[ei]->getAllComponents(),num_se*sizeof(node_handle));


      for(int jz = 0; jz < jC; jz ++) {
        int j;
        if(imFlag) {
         j =  Rn[event_Ru_Rn_index[ei]]->nextOf(i);
         if(j==-1) continue; // Not enabled
         if (j < nb.getSize() && -1==nb.down(j)) continue;
         seHandlesForRecursion[0]=Rn[event_Ru_Rn_index[ei]]->getDown();
        } else {
          j = Rp->index(jz);
          if(j==-1) continue; // Not enabled
          if (j < nb.getSize() && -1==nb.down(j)) continue;
          seHandlesForRecursion[which_se+does_rn_exist] = Rp->down(jz);
        }

        node_handle rec = recFire(nb.down(i), seHandlesForRecursion, num_se);

        if (rec == 0) continue;

        //confirm local state
        if(!rel->isConfirmedState(level,j))
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
          // nb.d_ref(j) = rec;
        }
        else {
          nbdj.set(nb.down(j));  // clobber
          newst.set(rec);     // clobber
          mddUnion->computeTemp(nbdj, newst, nbdj);
          updated = (nbdj.getNode() != nb.down(j));
          nb.setFull(j, nbdj);
        }
        if (updated) queue->add(j);

      } //for all j's

  } // for all events, ei

 }// more indexes to explore

   recycle(queue);

}

#if 0
MEDDLY::node_handle MEDDLY::forwd_hyb_dfs_by_events_mt::recFireSet(MEDDLY::node_handle mdd,
                                                                    std::set<node_handle> vector_mxd)
{
  dd_edge ans(resF), union_rec(resF);

  for(auto rn_it = vector_mxd.begin(); rn_it != vector_mxd.end(); rn_it ++){
    std::set<node_handle>  one_item;
    one_item.insert(*rn_it);
    node_handle r_ans = recFire(mdd,one_item);
    ans.set(r_ans);
    mddUnion->compute(union_rec, ans, union_rec);
    one_item.clear();
  }

  return union_rec.getNode();
}
#endif


// Same as post-image, except we saturate before reducing.
MEDDLY::node_handle MEDDLY::forwd_hyb_dfs_by_events_mt::recFire(
        MEDDLY::node_handle mdd, node_handle* seHandles, int num_se)
{
  // termination conditions
  std::pair<MEDDLY::node_handle, int> mxd_whichse_pair = getHighestNodeHandles(seHandles, num_se);
  MEDDLY::node_handle mxd = mxd_whichse_pair.first;
  int which_se  = mxd_whichse_pair.second;

  if (mxd == 0 || mdd == 0) return 0;

 // printf("\n recFire: mxd<%d>: %d",arg2F->getNodeLevel(mxd), mxd);
  if (mxd==-1) {
    if (argF->isTerminalNode(mdd)) {
      return resF->handleForValue(1);
    }
    // mxd is identity
    if (argF == resF)
      {
       return resF->linkNode(mdd);
      }
  }


  // check the cache
  MEDDLY::node_handle result = 0;
  ct_entry_key* Key = findResult(mdd, seHandles, num_se, result);
  if (0==Key) return result;

  #ifdef TRACE_RECFIRE
  printf("computing recFire(%d, %d)\n", mdd, mxd);
  printf("  node %3d ", mdd);
  argF->showNode(stdout, mdd, 1);
  printf("\n  node %3d ", mxd);
  relF()->showNode(stdout, mxd, 1);
  printf("\n");
  #endif

  // check if mxd and mdd are at the same level
  const int mddLevel = argF->getNodeLevel(mdd);
  const int mxdLevel = relF()->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);
  // domain* dm = resF->getDomain();

  dd_edge nbdj(resF), newst(resF);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New(argF, FULL_ONLY);
  if (mddLevel < rLevel) {
    A->initRedundant(argF, rLevel, mdd, FULL_ONLY);
  } else {
    argF->unpackNode(A, mdd, FULL_ONLY);
  }

  //Re-Think
  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (int i=0; i<rSize; i++) {
      nb->setFull(i, recFire(A->down(i), seHandles, num_se));
      // nb->d_ref(i) = recFire(A->down(i), seHandles, num_se);
    }

  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // Initialize mxd readers, note we might skip the unprimed level

    bool imFlag = relF()->isImplicit(mxd);
    int row_size = rSize;
    relation_node* relNode;
    unpacked_node *Ru = unpacked_node::New(relF(), SPARSE_ONLY);
    unpacked_node *Rp = unpacked_node::New(relF(), SPARSE_ONLY);

    if(imFlag) {
      relNode = relF()->buildImplicitNode(mxd);
      row_size = rSize; // Get the update array size;
    }
    else {
       if (mxdLevel < 0) {
        Ru->initRedundant(relF(), rLevel, mxd, SPARSE_ONLY);
      } else {
        relF()->unpackNode(Ru, mxd, SPARSE_ONLY);
      }
      row_size = Ru->getSize();
    }


    // for intersection to happen : For each i, check across the set_mxd what are the common j's.
    // Retain the identity behavior
    // Then use those j's only.
    node_handle* seHandlesLower = (node_handle*)malloc(sizeof(MEDDLY::node_handle)*num_se);
    memcpy(seHandlesLower, seHandles, sizeof(MEDDLY::node_handle)*num_se);
  //  seHandlesLower[0] = 0;

    //for (int iz=0; iz<row_size; iz++) { // for each i -> get common j's downhandle

      // Initialize mxd readers, note we might skip the unprimed level
      if(imFlag) {
        for (int iz=0; iz<row_size; iz++) {
        const unsigned i = iz;
        int j = relNode->nextOf(i);
        if(j==-1) continue;
        recFireHelper(i, rLevel, relNode->getDown(), A->down(i), Rp, nb, seHandlesLower, imFlag, j, num_se, which_se);
       }
       } else  {
        for (int iz=0; iz<row_size; iz++) {
        const unsigned i = Ru->index(iz);
        if (0==A->down(i)) continue;
        recFireHelper(i, rLevel, Ru->down(iz), A->down(i), Rp, nb, seHandlesLower, imFlag, -1, num_se, which_se);
        }
        unpacked_node::Recycle(Rp);
        unpacked_node::Recycle(Ru);
    }

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

  return saveResult(Key, mdd, seHandles, num_se, result);
}

void MEDDLY::forwd_hyb_dfs_by_events_mt::recFireHelper(
  const unsigned i,
  const int rLevel,
  const node_handle Ru_i,
  const node_handle A_i,
  unpacked_node *Rp,
  unpacked_node* nb,
  node_handle* recHandles,
  bool imFlag,
  int j,
  int num_se, int which_se)
{


  if(!imFlag) {
  if(-rLevel == relF()->getNodeLevel(Ru_i))
     relF()->unpackNode(Rp, Ru_i, SPARSE_ONLY);
  else
     return;
  }

  dd_edge nbdj(resF), newst(resF);

  unsigned jc = imFlag ? 1 : Rp->getSize();

  for (unsigned jz=0; jz<jc; jz++) {
    if((jz == 0) && imFlag)
      {
        if (j == -1) continue;
        recHandles[which_se] = Ru_i;
      }
    else {
        j = Rp->index(jz);
        if (j == -1) continue;
        recHandles[which_se] = Rp->down(jz);
      }

      MEDDLY::node_handle newstates = recFire(A_i, recHandles, num_se);
      if (0==newstates) continue;

      //confirm local state
      if(!rel->isConfirmedState(rLevel,j)) // if not confirmed before
      {
        rel->confirm(rLevel,j); // confirm and enlarge
        if (j >= nb->getSize()) {
          int new_var_bound = resF->isExtensibleLevel(nb->getLevel())? -(j+1) : (j+1);
          domain* dm = resF->getDomain();
          dm->enlargeVariableBound(nb->getLevel(), false, new_var_bound);
          const unsigned oldSize = nb->getSize();
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

  }


}


// ******************************************************************
// *                                                                *
// *             common_hyb_dfs_by_events_mt  methods              *
// *                                                                *
// ******************************************************************

MEDDLY::common_hyb_dfs_by_events_mt::common_hyb_dfs_by_events_mt(
    const char* name, forest* argF, hybrid_relation* relation, forest* resF)
    : saturation_operation(name, 1, argF, resF)
{
    if (!argF || !relation || !resF) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }

    rel = relation;

    if (relation->getInForest() != argF) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    if (relation->getOutForest() != resF) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }

    // for now, anyway, inset and outset must be same forest
    if (argF != resF) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }

    // Check for same domain
    if  (
            (argF->getDomain() != relF()->getDomain()) ||
            (resF->getDomain() != relF()->getDomain())
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

  registerInForest(relF());

  ct_entry_type* et = new ct_entry_type(name, "N.N:N");
  et->setForestForSlot(0, argF);
  et->setForestForSlot(2, relF());
  et->setForestForSlot(4, resF);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::common_hyb_dfs_by_events_mt::~common_hyb_dfs_by_events_mt()
{
  unregisterInForest(relF());
  delete rel;
}



void MEDDLY::common_hyb_dfs_by_events_mt
::compute(const dd_edge &a, dd_edge &c)
{
  // Initialize operations
  mddUnion = UNION(resF, resF, resF);
  MEDDLY_DCASSERT(mddUnion);


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

  saturation_hyb_by_events_op* so = new saturation_hyb_by_events_op(this,
          "Hybrid_saturation", argF, resF);

  MEDDLY::node_handle cnode = so->saturate(a.getNode());
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
// *       common_hyb_dfs_by_events_mt::indexq  methods                 *
// ******************************************************************

MEDDLY::common_hyb_dfs_by_events_mt::indexq::indexq()
{
  data = 0;
  size = 0;
  head = NULPTR;
}

MEDDLY::common_hyb_dfs_by_events_mt::indexq::~indexq()
{
  free(data);
}

void MEDDLY::common_hyb_dfs_by_events_mt::indexq::resize(int sz)
{
  if (sz <= size) return;
  data = (int*) realloc(data, sz * sizeof(int));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);

  for (; size < sz; size++) data[size] = NOTINQ;
}

// ******************************************************************
// *       common_hyb_dfs_by_events_mt::charbuf methods                 *
// ******************************************************************

MEDDLY::common_hyb_dfs_by_events_mt::charbuf::charbuf()
{
  data = 0;
  size = 0;
}

MEDDLY::common_hyb_dfs_by_events_mt::charbuf::~charbuf()
{
  free(data);
}

void MEDDLY::common_hyb_dfs_by_events_mt::charbuf::resize(int sz)
{
  if (sz <= size) return;
  data = (char*) realloc(data, sz * sizeof(char));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *               saturation_hyb_by_events_op  methods            *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_hyb_by_events_op
::saturation_hyb_by_events_op(common_hyb_dfs_by_events_mt* p,
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

MEDDLY::saturation_hyb_by_events_op::~saturation_hyb_by_events_op()
{
  // removeAllComputeTableEntries();
}


void MEDDLY::saturation_hyb_by_events_op::compute(const dd_edge &arg,
        dd_edge &res)
{
    const node_handle mdd = arg.getNode();
  // Saturate

#ifdef DEBUG_INITIAL
    printf("Calling saturate for states:\n");
    ostream_output s(std::cout);
    argF->showNodeGraph(s, &mdd, 1);
    std::cout.flush();
#endif
    res.set(saturate(mdd, argF->getMaxLevelIndex()));
}

MEDDLY::node_handle
MEDDLY::saturation_hyb_by_events_op::saturate(node_handle mdd)
{
  // Saturate

#ifdef DEBUG_INITIAL
    printf("Calling saturate for states:\n");
    ostream_output s(std::cout);
    argF->showNodeGraph(s, &mdd, 1);
    std::cout.flush();
#endif
    return saturate(mdd, argF->getMaxLevelIndex());
}

MEDDLY::node_handle
MEDDLY::saturation_hyb_by_events_op::saturate(node_handle mdd, int k)
{
#ifdef DEBUG_DFS
  printf("mdd: %d, k: %d\n", mdd, k);
#endif


  // terminal condition for recursion
  if (argF->isTerminalNode(mdd)) return mdd;

  // search compute table
  MEDDLY::node_handle n = 0;
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
  unpacked_node *mddDptrs = unpacked_node::New(argF, FULL_ONLY);
  if (mdd_level < k) {
    mddDptrs->initRedundant(argF, k, mdd, FULL_ONLY);
  } else {
    argF->unpackNode(mddDptrs, mdd, FULL_ONLY);
  }

  // Do computation
  for (int i=0; i<sz; i++) {
    nb->setFull(i, mddDptrs->down(i) ? saturate(mddDptrs->down(i), k-1) : 0);
//     nb->d_ref(i) = mddDptrs->down(i) ? saturate(mddDptrs->down(i), k-1) : 0;
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
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

/*

MEDDLY::sathyb_opname* MEDDLY::initHybSaturationForward()
{
  return new sathyb_opname("SaturationFwd");
}


MEDDLY::specialized_operation*
MEDDLY::sathyb_opname::buildOperation(arguments* a)
{

  hybrid_relation* rel = dynamic_cast<hybrid_relation*>(a);
  if (0==rel) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);

  MEDDLY::specialized_operation* op = 0;
  op = new forwd_hyb_dfs_by_events_mt(this, rel);

  return op;
}
*/

MEDDLY::saturation_operation*
MEDDLY::SATURATION_HYB_FORWARD(forest* iF, hybrid_relation* nsf, forest* oF)
{
    return new forwd_hyb_dfs_by_events_mt("HybridSat", iF, nsf, oF);
}



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
#include "old_meddly.h"
#include "old_meddly.hh"
#include "old_meddly_expert.h"
#include "old_meddly_expert.hh"
#include "relation_node.h"
#include "sat_impl.h"
#include <typeinfo> // for "bad_cast" exception
#include <set>
#include <map>

#include "ct_entry_result.h"
#include "opname_satur.h"

   #define OUT_OF_BOUNDS -1
   #define NOT_KNOWN -2
   #define TERMINAL_NODE 1


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
// *                     satimpl_opname  methods                    *
// *                                                                *
// ******************************************************************



MEDDLY::satimpl_opname::satimpl_opname(const char* n)
: specialized_opname(n)
{
}

MEDDLY::satimpl_opname::~satimpl_opname()
{
}

// ******************************************************************

MEDDLY::satimpl_opname::implicit_relation::implicit_relation(forest* inmdd, forest* relmxd,
                                                             forest* outmdd)
: insetF(static_cast<expert_forest*>(inmdd)), outsetF(static_cast<expert_forest*>(outmdd)), mixRelF(static_cast<expert_forest*>(relmxd))
{
  mixRelF = static_cast<expert_forest*>(relmxd);

  if (0==insetF || 0==outsetF || 0==mixRelF ) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);

  // Check for same domain
  if (insetF->getDomain() != outsetF->getDomain())
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  // for now, anyway, inset and outset must be same forest
  if (insetF != outsetF)
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

  // Check forest types
  if (
      insetF->isForRelations()    ||
      outsetF->isForRelations()   ||
      (insetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)   ||
      (outsetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
      )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  // Forests are good; set number of variables
  num_levels = insetF->getDomain()->getNumVariables();



  //Allocate event_list
  event_list = (rel_node_handle**)malloc(unsigned(num_levels+1)*sizeof(rel_node_handle*));
  event_list_alloc = (long*)malloc(unsigned(num_levels+1)*sizeof(long));
  event_added = (long*)malloc(unsigned(num_levels+1)*sizeof(long));


  confirm_states = (long*)malloc(unsigned(num_levels+1)*sizeof(long));
  confirmed_array_size = (long*)malloc(unsigned(num_levels+1)*sizeof(long));
  confirmed = new bool*[num_levels+1];

  confirmed[0]=0;
  for(int i = 1;i<=num_levels;i++)
    {
    event_list[i] = (rel_node_handle*)malloc(8*sizeof(rel_node_handle));
    confirmed[i] = (bool*)malloc(insetF->getVariableSize(i)*sizeof(bool));
    event_list_alloc[i] = 8;
    event_added[i] = 0;
    confirm_states[i] = 0;

    confirmed_array_size[i]=insetF->getVariableSize(i);
    for(int j = 0;j<insetF->getVariableSize(i);j++)
      confirmed[i][j]=false;
    }



  //create the terminal node
  relation_node *Terminal = new MEDDLY::relation_node(0,mixRelF,0,-1,0,0,-1);
  //mixRelF->createRelationNode(Terminal);
  Terminal->setID(TERMINAL_NODE);
  std::pair<rel_node_handle, relation_node*> TerminalNode(TERMINAL_NODE,Terminal);
  impl_unique.insert(TerminalNode);
  last_in_node_array = TERMINAL_NODE;

}


void
MEDDLY::satimpl_opname::implicit_relation::resizeEventArray(int level)
{
  event_added[level] += 1;
  if (event_added[level] > event_list_alloc[level]) {
    int nalloc = ((event_added[level]/8)+1)*8;
    MEDDLY_DCASSERT(nalloc > 0);
    MEDDLY_DCASSERT(nalloc > event_added[level]);
    MEDDLY_DCASSERT(nalloc > event_list_alloc[level]);
    event_list[level] = (rel_node_handle*) realloc(event_list[level], nalloc*sizeof(rel_node_handle));
    if (0==event_list[level]) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    event_list_alloc[level] = nalloc;
  }
}

void
MEDDLY::satimpl_opname::implicit_relation::resizeConfirmedArray(int level,int index)
{
  int nalloc = index+1;
 if(nalloc>confirmed_array_size[level])
    {

       MEDDLY_DCASSERT(nalloc > 0);
       MEDDLY_DCASSERT(confirmed_array_size[level] >= 0);
       if(confirmed_array_size[level]==0)
         {
           confirmed[level] = (bool*)malloc(nalloc*sizeof(bool));
           if (0==confirmed[level]) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
         }
        else
          {
            confirmed[level] = (bool*)realloc(confirmed[level], nalloc*sizeof(bool));
            if (0==confirmed[level]) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
          }

        for(int i = confirmed_array_size[level];i<nalloc;i++)
          confirmed[level][i]=false;

         confirmed_array_size[level]=nalloc;
    }

}

void findConfirmedStatesImpl(MEDDLY::satimpl_opname::implicit_relation* rel,
                             bool** confirmed, long* confirm_states,
                             MEDDLY::node_handle mdd, int level,
                             std::set<MEDDLY::node_handle>& visited) {
  if (level == 0) return;
  if (visited.find(mdd) != visited.end()) return;

  MEDDLY::expert_forest* insetF = rel->getInForest();
  int mdd_level = insetF->getNodeLevel(mdd);
  if (MEDDLY::isLevelAbove(level, mdd_level)) {
    // skipped level; confirm all local states at this level
    // go to the next level
    int level_size = insetF->getLevelSize(level);
    for (int i = 0; i < level_size; i++) {
      if (!confirmed[level][i]) {
        rel->setConfirmedStates(level, i);
      }
    }
    findConfirmedStatesImpl(rel, confirmed, confirm_states, mdd, level-1, visited);
  } else {
    if (MEDDLY::isLevelAbove(mdd_level, level)) {
      throw MEDDLY::error(MEDDLY::error::INVALID_VARIABLE, __FILE__, __LINE__);
    }
    // mdd_level == level
    visited.insert(mdd);
    MEDDLY::unpacked_node *nr = insetF->newUnpacked(mdd, MEDDLY::SPARSE_ONLY);
    for (int i = 0; i < nr->getNNZs(); i++) {
      if (!confirmed[level][nr->i(i)]) {
        rel->setConfirmedStates(level, nr->i(i));
      }
      findConfirmedStatesImpl(rel, confirmed, confirm_states, nr->d(i), level-1, visited);
    }
    MEDDLY::unpacked_node::recycle(nr);
  }
}

void MEDDLY::satimpl_opname::implicit_relation::setConfirmedStates(const dd_edge& set)
{
  // Perform a depth-first traversal of set:
  //    At each level, mark all enabled states as confirmed.

  // Enlarge the confirmed arrays if needed
  for (int i = 1 ; i<=num_levels; i++)
    {
      int levelSize = getInForest()->getLevelSize(i);
      resizeConfirmedArray(i, levelSize);
    }

    std::set<node_handle> visited;
    findConfirmedStatesImpl(const_cast<implicit_relation*>(this),
                      confirmed, confirm_states, set.getNode(), num_levels, visited);

}



MEDDLY::satimpl_opname::implicit_relation::~implicit_relation()
{
  last_in_node_array = 0;
  impl_unique.clear();

  for(int i = 0; i <=num_levels; i++) {delete[] event_list[i]; delete[] confirmed[i];}
  delete[] event_list;
  delete[] event_added;
  delete[] event_list_alloc;
  delete[] confirmed;
  delete[] confirm_states;
  delete[] confirmed_array_size;
}


    MEDDLY::rel_node_handle
MEDDLY::satimpl_opname::implicit_relation::isUniqueNode(relation_node* n)
{
  bool is_unique_node = true;
  std::unordered_map<rel_node_handle, relation_node*>::iterator it = impl_unique.begin();
  while(it != impl_unique.end())
    {
    is_unique_node = !((it->second)->equals(n));
    if(is_unique_node==false)
      return (it->second)->getID();
    ++it;
    }
  return 0;
}

    MEDDLY::rel_node_handle
MEDDLY::satimpl_opname::implicit_relation::registerNode(bool is_event_top, relation_node* n)
{

  rel_node_handle nLevel = n->getLevel();

#ifdef DEVELOPMENT_CODE
  rel_node_handle downHandle = n->getDown();
  relation_node* downNode = nodeExists(downHandle);
  rel_node_handle downLevel = downNode->getLevel();
  MEDDLY_DCASSERT( ( ( downNode!=NULL ) && ( nLevel > downLevel ) )
                    ||
                   ( downLevel == 0 ) );
#endif

  rel_node_handle n_ID = isUniqueNode(n);

  if(n_ID==0) // Add new node
   {
    n_ID  = last_in_node_array + 1;
    std::pair<rel_node_handle, relation_node*> add_node(n_ID,n);
    impl_unique.insert(add_node);
    if(impl_unique.find(n_ID) != impl_unique.end())
    {
      last_in_node_array = n_ID;
      n->setID(n_ID);
    }
    //mixRelF->createRelationNode(n);
  }
  else //Delete the node
    {
     delete n;
    }

  if(is_event_top)
    {
    resizeEventArray(nLevel);
    event_list[nLevel][event_added[nLevel] - 1] = n_ID;
    }

  return n_ID;
}

void
MEDDLY::satimpl_opname::implicit_relation::show()
{
  rel_node_handle** event_list_copy = (rel_node_handle**)malloc((num_levels+1)*sizeof(rel_node_handle*));
  if (0==event_list_copy) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  long total_events = 0;
  for(int i = 1;i<=num_levels;i++) total_events +=event_added[i];
  for(int i = 1;i<=num_levels;i++)
    {
     event_list_copy[i] = (rel_node_handle*)malloc(total_events*sizeof(rel_node_handle));
     if (0==event_list_copy[i]) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }

  for(int i = num_levels;i>=1;i--)
    for(int j=0;j<total_events;j++)
      event_list_copy[i][j]=0;


  int eid = 0;
  for(int i = num_levels;i>=1;i--)
    {
     int k = 0;
     std::cout<<"\n [";
     for(int j=0;j<total_events;j++)
      {

        if((event_list_copy[i][eid]==0)&&(k<event_added[i]))
          {
            event_list_copy[i][eid] = event_list[i][k];
            relation_node* hold_it = nodeExists(event_list[i][k]);
            relation_node* hold_down = nodeExists(hold_it->getDown());
            event_list_copy[hold_down->getLevel()][eid] = hold_down->getID();
          k++;eid++;
          }


      int dig_ctr = event_list_copy[i][j]>1000?4:(event_list_copy[i][j]>100?3:(event_list_copy[i][j]>10?2:1));

      int spc_bef =(6 - dig_ctr)/2;
      int spc_aft = 6 - dig_ctr - spc_bef;

      for(int s=0;s<spc_bef;s++) std::cout<<" ";
      if(event_list_copy[i][j] != 0) std::cout<<event_list_copy[i][j];
      else std::cout<<"_";
      for(int s=0;s<spc_aft;s++) std::cout<<" ";
      if(j!=total_events-1)
          std::cout<<"|";
      }
     std::cout<<"]";
    }

  for(int i = 0;i<num_levels+1;i++) delete event_list_copy[i];
  delete[] event_list_copy;

}

void MEDDLY::satimpl_opname::implicit_relation::bindExtensibleVariables() {
  //
  // Find the bounds for each extensbile variable
  //
  expert_domain* ed = static_cast<expert_domain*>(outsetF->useDomain());

  for (int k = 1; k <= num_levels; k++) {
    int bound = 0;
    int n_confirmed = 0;

    for (int j = 0; j < confirmed_array_size[k]; j++) {
      if (confirmed[k][j]) { bound = j+1; n_confirmed++; }
    }

    MEDDLY_DCASSERT(bound > 0);
    MEDDLY_DCASSERT(n_confirmed == confirm_states[k]);
    ed->enlargeVariableBound(k, false, bound);
  }
}

MEDDLY::node_handle
MEDDLY::satimpl_opname::implicit_relation::buildMxdForest()
{

  //Get number of Variables and Events
  int nVars = outsetF->getDomain()->getNumVariables();
  int nEvents = getTotalEvent(nVars);


  rel_node_handle* event_tops = (rel_node_handle*)malloc((nEvents)*sizeof(rel_node_handle));
  int e = 0;

  for(int i = 1 ;i<=nVars;i++)
    {
    int num_events_at_this_level = lengthForLevel(i);
    for(int j = 0;j<num_events_at_this_level;j++)
      event_tops[e++]=arrayForLevel(i)[j];
    }

  domain *d = outsetF->useDomain();

  forest* mxd = d->createForest(true,range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
  dd_edge nsf(mxd);

  dd_edge* monolithic_nsf = new dd_edge(mxd);
  for(int i=0;i<nEvents;i++)
    {
    (*monolithic_nsf) += buildEventMxd(event_tops[i],mxd);
    }

  node_handle monolithic_nsf_handle = monolithic_nsf->getNode();
  mxdF = (expert_forest*)mxd;

  /*for(int i = 0; i<nEvents;i++)
   {
   dd_edge nsf_ev(mxd);
   nsf_ev = buildEventMxd(event_tops[i],mxd);
   apply(UNION, nsf, nsf_ev, nsf);
   }*/

  return monolithic_nsf_handle;
}


MEDDLY::dd_edge
MEDDLY::satimpl_opname::implicit_relation::buildEventMxd(rel_node_handle eventTop, forest *mxd)
{
  //mxd is built on a domain obtained from result of saturation
  int nVars = outsetF->getDomain()->getNumVariables();
  //int* sizes = new int[nVars];
  relation_node* Rnode = nodeExists(eventTop);
  rel_node_handle* rnh_array = (rel_node_handle*)malloc((nVars+1)*sizeof(rel_node_handle));
  // int top_level = Rnode->getLevel();

  expert_forest* ef = (expert_forest*) mxd;

  //Get relation node handles
  for (int i=nVars; i>=1; i--)
    {
      if(Rnode->getLevel()==i)// if variable i is a part of this event
        {
          rnh_array[i] = Rnode->getID(); // keep track of node_handles that are part of this event
          Rnode = nodeExists(Rnode->getDown()); // move to next variable in the event
        }
      else // if not, then
        {
        rnh_array[i] = -1; // node handle of the variable i in the event
        continue;
        }
    }

  node_handle below = ef->handleForValue(true); // Terminal true node

  for (int i=1; (i<=nVars)&&(below!=0); i++)
    {
        if(rnh_array[i]!=-1)
          {
            const int maxVar = ef->getDomain()->getVariableBound(i);
            Rnode = nodeExists(rnh_array[i]);
            //Create a new unprimed node for variable i
            MEDDLY_DCASSERT(outsetF->getVariableSize(i)>=Rnode->getPieceSize());
            unpacked_node* UP_var = unpacked_node::newFull(ef, i, Rnode->getPieceSize());

            for (int j=0; j<Rnode->getPieceSize(); j++) {

              // long new_j = confirmed[i][Rnode->getTokenUpdate()[j]]? Rnode->getTokenUpdate()[j] : -2;
              long new_j = Rnode->getTokenUpdate()[j];

              if(new_j>=0 && new_j<maxVar) // do not exceed variable bounds
                {
                   //Create primed node for each valid index of the unprimed node
                  unpacked_node* P_var = unpacked_node::newSparse(ef, -i, 1);
                  P_var->i_ref(0) = new_j;
                  P_var->d_ref(0) = ef->linkNode(below); // primed node for new_j index points to terminal or unprime node
                  UP_var->d_ref(j) = ef->createReducedNode(j, P_var);
                }
              else
                UP_var->d_ref(j) = ef->handleForValue(false); // unprimed node for j index points to false
              }

              ef->unlinkNode(below);
              below = ef->createReducedNode(-1, UP_var);
          }
    }

  dd_edge nsf(mxd);
  nsf.set(below);

  return nsf;
}

// ******************************************************************


std::unordered_map<long,std::vector<MEDDLY::rel_node_handle>>
MEDDLY::satimpl_opname::implicit_relation::getListOfNexts(int level, long i, relation_node **R)
{
  std::unordered_map<long,std::vector<rel_node_handle>> jList;
  // atleast as many j's as many events
  for(int k=0;k<lengthForLevel(level);k++)
    {
    long key = R[k]->nextOf(i);
    jList[key].reserve(lengthForLevel(level));
    int rnh_dwn = R[k]->getDown();
    jList[key].push_back(rnh_dwn);
    }

  return jList;
}

bool
MEDDLY::satimpl_opname::implicit_relation::isUnionPossible(int level, long i, relation_node **R)
{
  if(lengthForLevel(level)==1)
     return false;

  std::vector<int> jset(lengthForLevel(level), 0);

   int last_j = 0;
   for(int k=0;k<lengthForLevel(level);k++)
    {
    long key = R[k]->nextOf(i);
    int flag = 0;
    for(int m=0;m<last_j;m++)
      if(jset[m]==key)
        {
          flag=1;
          break;
        }

      if(flag==0)
        {
          jset[k]=key;
          last_j++;
        }
    }
  if(lengthForLevel(level)==last_j)
   return false;
  else
    return true;
}

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

class MEDDLY::common_impl_dfs_by_events_mt : public specialized_operation {
public:
  common_impl_dfs_by_events_mt(const satimpl_opname* opcode,
                               satimpl_opname::implicit_relation* rel);
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

  dd_edge nbdj(resF), newst(resF);

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
          nbdj.set(nb.d(j));  // clobber
          newst.set(rec);     // clobber
          mddUnion->computeTemp(nbdj, newst, nbdj);
          updated = (nbdj.getNode() != nb.d(j));
          nb.set_d(j, nbdj);
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
        nbdj.set(nb.d(j));  // clobber
        newst.set(rec);     // clobber
        mddUnion->computeTemp(nbdj, newst, nbdj);
        updated = (nbdj.getNode() != nb.d(j));
        nb.set_d(j, nbdj);
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
  const int mxdLevel = relNode->getLevel();
  const int rLevel = MAX(mxdLevel, mddLevel);
   int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);
  expert_domain* dm = static_cast<expert_domain*>(resF->useDomain());

  dd_edge nbdj(resF), newst(resF);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New();
  if (mddLevel < rLevel) {
    A->initRedundant(arg1F, rLevel, mdd, true);
  } else {
    arg1F->unpackNode(A, mdd, FULL_ONLY);
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
          nbdj.set(nb->d(j));
          newst.set(newstates);
          mddUnion->computeTemp(nbdj, newst, nbdj);
          nb->set_d(j, nbdj);

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
: specialized_operation(opcode, 1)
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
  ct_entry_type* et = new ct_entry_type(opcode->getName(), "NL:N");
  et->setForestForSlot(0, arg1F);
  et->setForestForSlot(3, resF);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::common_impl_dfs_by_events_mt::~common_impl_dfs_by_events_mt()
{
  if (rel->autoDestroy()) delete rel;
  unregisterInForest(arg1F);
  //unregisterInForest(arg2F);
  unregisterInForest(resF);
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
: unary_operation(saturation_impl_by_events_opname::getInstance(), 1, argF, resF)
{
  parent = p;

  const char* name = saturation_impl_by_events_opname::getInstance()->getName();
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
  unpacked_node *mddDptrs = unpacked_node::New();
  if (mdd_level < k) {
    mddDptrs->initRedundant(argF, k, mdd, true);
  } else {
    argF->unpackNode(mddDptrs, mdd, FULL_ONLY);
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
  auto entry_result = intersection_cache.find(key);
  if (entry_result != intersection_cache.end()) {
    // found cached entry
    return entry_result->second;
  }

  // unpack nodes
  MEDDLY::unpacked_node* unp_A =
    mddF->getNodeLevel(node_A) >= mddF->getNodeLevel(node_B)
    ? mddF->newUnpacked(node_A, MEDDLY::FULL_ONLY)
    : MEDDLY::unpacked_node::newRedundant(mddF, mddF->getNodeLevel(node_B), node_A, true);
  MEDDLY::unpacked_node* unp_B =
    mddF->getNodeLevel(node_B) >= mddF->getNodeLevel(node_A)
    ? mddF->newUnpacked(node_B, MEDDLY::FULL_ONLY)
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
  ct_entry_key* Key = findSaturateResult(mdd, k, n);
  if (0==Key) {
    saturation_result = n;
    return !isIntersectionEmpty(argF, saturation_result, constraint);
  }

  const int mdd_level = argF->getNodeLevel(mdd);
  const int constraint_level = argF->getNodeLevel(constraint);

  // Initialize mdd reader
  unpacked_node *mddDptrs = unpacked_node::New();
  if (mdd_level < k) {
    mddDptrs->initRedundant(argF, k, mdd, true);
  } else {
    argF->unpackNode(mddDptrs, mdd, FULL_ONLY);
  }
  MEDDLY_DCASSERT(!mddDptrs->isExtensible());

  // Initialize constraint reader
  unpacked_node* consDptrs = unpacked_node::New();
  if (constraint_level < k) {
    consDptrs->initRedundant(argF, k, constraint, true);
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
  unpacked_node* consDptrs = unpacked_node::New();
  if (constraint_level < level) {
    consDptrs->initRedundant(arg1F, level, constraint, true);
  } else {
    arg1F->unpackNode(consDptrs, constraint, FULL_ONLY);
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

  dd_edge nbdj(resF), newst(resF);

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
        nbdj.set(nb.d(j));  // clobber
        newst.set(rec);     // clobber
        mddUnion->computeTemp(nbdj, newst, nbdj);
        updated = (nbdj.getNode() != nb.d(j));
        nb.set_d(j, nbdj);
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
  ct_entry_key* Key = findResult(mdd, mxd, result);
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

  dd_edge nbdj(resF), newst(resF);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New();
  if (mddLevel < rLevel) {
    A->initRedundant(arg1F, rLevel, mdd, true);
  } else {
    arg1F->unpackNode(A, mdd, FULL_ONLY);
  }

  // Initialize constraint reader
  unpacked_node* consDptrs = unpacked_node::New();
  if (constraint_level < rLevel) {
    consDptrs->initRedundant(arg1F, rLevel, constraint, true);
  } else {
    arg1F->unpackNode(consDptrs, constraint, FULL_ONLY);
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
      nbdj.set(nb->d(j));
      newst.set(newstates);
      mddUnion->computeTemp(nbdj, newst, nbdj);
      nb->set_d(j, nbdj);
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


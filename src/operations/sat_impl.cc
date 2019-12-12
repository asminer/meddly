
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
#include <vector>

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

MEDDLY::relation_node::relation_node(unsigned long sign, forest* f, int lvl, node_handle d, long enable_val, long fire_val) :
f(static_cast<expert_forest*>(f))
{
  signature  = sign;
  level = lvl;
  down = d;
  enable = enable_val;
  fire = fire_val;
  piece_size = 0;
  token_update = NULL;
}

MEDDLY::relation_node::~relation_node()
{
}


long MEDDLY::relation_node::nextOf(long i)
{
  //to be defined for the example you use & comment this definition
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
  //if(i >= enable) return i + fire;
  //else return -1;
}

bool
MEDDLY::relation_node::equals(const relation_node* n) const
{
  if((level == n->getLevel()) && (down == n->getDown()) 
     && (fire == n->getFire()) && (enable == n->getEnable()))
    return true;
  else
    return false;
}

void
MEDDLY::relation_node::expandTokenUpdate(long i)
{
  if(getPieceSize()==0)
  {
    token_update = (long*)malloc(1*sizeof(long));
    piece_size = 1;
    token_update[0]=NOT_KNOWN;
  }
  if(i>0)
  {
    token_update = (long*)realloc(token_update,(i+1)*sizeof(long));
    for(int j = piece_size;j<=i;j++)
      token_update[j]=NOT_KNOWN;
    piece_size = i+1;
  }
}

void
MEDDLY::relation_node::setTokenUpdateAtIndex(long i,long val)
{
  MEDDLY_DCASSERT(i<getPieceSize());
  token_update[i] = val;
}
// ******************************************************************

MEDDLY::satimpl_opname::implicit_relation::implicit_relation(forest* inmdd, forest* relmxd, forest* outmdd, event** E, int ne) 
: insetF(static_cast<expert_forest*>(inmdd)), outsetF(static_cast<expert_forest*>(outmdd))
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
      !mixRelF->isForRelations()  ||
      (insetF->getEdgeLabeling() != forest::MULTI_TERMINAL)   ||
      (outsetF->getEdgeLabeling() != forest::MULTI_TERMINAL)  ||
      (mixRelF->getEdgeLabeling() != forest::MULTI_TERMINAL)  
      )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
  
  // Forests are good; set number of variables
  num_levels = insetF->getDomain()->getNumVariables() + 1;
  
  //Allocate event_list
  event_list = (node_handle**)malloc(unsigned(num_levels)*sizeof(node_handle*));
  event_list_alloc = (long*)malloc(unsigned(num_levels)*sizeof(long));
  event_added = (long*)malloc(unsigned(num_levels)*sizeof(long));

  //Confirmed local states
  confirm_states = (long*)malloc(unsigned(num_levels)*sizeof(long));
  confirmed_array_size = (long*)malloc(unsigned(num_levels)*sizeof(long));
  confirmed = new bool*[num_levels];
  confirmed[0]=0;
  for(int i = 1;i<num_levels;i++)
    {
    event_list[i] = (node_handle*)malloc(8*sizeof(node_handle));
    confirmed[i] = (bool*)malloc(insetF->getVariableSize(i)*sizeof(bool));
    event_list_alloc[i] = 8;
    event_added[i] = 0; 
    confirm_states[i] = 0;
    
    confirmed_array_size[i]=insetF->getVariableSize(i);
    for(int j = 0;j<insetF->getVariableSize(i);j++)
      confirmed[i][j]=false;
    }
  
  // Build the events-per-level data structure
  // (0) Initialize
  num_events_by_top_level = new int[num_levels];
  num_events_by_level = new int[num_levels];
  memset(num_events_by_top_level, 0, sizeof(int)*unsigned(num_levels));
  memset(num_events_by_level, 0, sizeof(int)*unsigned(num_levels));
  
  // (1) Determine the number of events per level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nVars = E[i]->getNumVars();
    const int* vars = E[i]->getVars();
    for (int j = 0; j < nVars; j++) {
      num_events_by_level[vars[j]]++;
    }
    num_events_by_top_level[E[i]->getTop()]++;
  }
  
  // (2) Allocate events[i]
  events_by_top_level = new event**[num_levels];
  events_by_level = new event**[num_levels];
  for (int i = 0; i < num_levels; i++) {
    events_by_top_level[i] = num_events_by_top_level[i] > 0
    ? new event*[num_events_by_top_level[i]]: 0;
    events_by_level[i] = num_events_by_level[i] > 0
    ? new event*[num_events_by_level[i]]: 0;
    num_events_by_top_level[i] = 0; // reset this; to be used by the next loop
    num_events_by_level[i] = 0; // reset this; to be used by the next loop
  }
  
  // (3) Store events by level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nVars = E[i]->getNumVars();
    const int* vars = E[i]->getVars();
    for (int j = 0; j < nVars; j++) {
      int level = vars[j];
      events_by_level[level][num_events_by_level[level]++] = E[i];
    }
    int level = E[i]->getTop();
    events_by_top_level[level][num_events_by_top_level[level]++] = E[i];
    E[i]->markForRebuilding();
  }
  
  // Build the subevents-per-level data structure
  // (0) Initialize
  num_subevents_by_level = new int[num_levels];
  memset(num_subevents_by_level, 0, sizeof(int)*unsigned(num_levels));
  
  // (1) Determine the number of subevents per level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nse = E[i]->getNumOfSubevents();
    subevent** se = E[i]->getSubevents();
    for (int j = 0; j < nse; j++) {
      int nVars = se[j]->getNumVars();
      const int* vars = se[j]->getVars();
      for (int k = 0; k < nVars; k++) {
        num_subevents_by_level[vars[k]]++;
      }
    }
  }
  
  // (2) Allocate subevents[i]
  subevents_by_level = new subevent**[num_levels];
  for (int i = 0; i < num_levels; i++) {
    subevents_by_level[i] = num_subevents_by_level[i] > 0
    ? new subevent*[num_subevents_by_level[i]]: 0;
    num_subevents_by_level[i] = 0; // reset this; to be used by the next loop
  }
  
  // (3) Store subevents by level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nse = E[i]->getNumOfSubevents();
    subevent** se = E[i]->getSubevents();
    for (int j = 0; j < nse; j++) {
      int nVars = se[j]->getNumVars();
      const int* vars = se[j]->getVars();
      for (int k = 0; k < nVars; k++) {
        int level = vars[k];
        subevents_by_level[level][num_subevents_by_level[level]++] = se[j];
      }
    }
  }
  
  
  // Build the relnodes-per-level data structure
  // (0) Initialize
  num_relnodes_by_level = new int[num_levels];
  memset(num_relnodes_by_level, 0, sizeof(int)*unsigned(num_levels));
  
  // (1) Determine the number of relnodes per level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nrn = E[i]->getNumOfRelnodes();
    relation_node** rn = E[i]->getRelNodes();
    for (int j = 0; j < nrn; j++) {
      int var = rn[j]->getLevel();
      num_relnodes_by_level[var]++;
      }
    }
  
  // (2) Allocate relnodes[i]
  relnodes_by_level = new MEDDLY::relation_node**[num_levels];
  for (int i = 0; i < num_levels; i++) {
    relnodes_by_level[i] = num_relnodes_by_level[i] > 0
    ? new relation_node*[num_relnodes_by_level[i]]: 0;
    num_relnodes_by_level[i] = 0; // reset this; to be used by the next loop
  }
  
  // (3) Store relnodes by level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nrn = E[i]->getNumOfRelnodes();
    relation_node** rn = E[i]->getRelNodes();
    for (int j = 0; j < nrn; j++) {
      int level = rn[j]->getLevel();
      relnodes_by_level[level][num_relnodes_by_level[level]++] = rn[j];
    }
  }
}

// *************************************Old Stuff*******************************************
void
MEDDLY::satimpl_opname::implicit_relation::resizeEventArray(int level)
{
  event_added[level] += 1;
  if (event_added[level] > event_list_alloc[level]) {
    int nalloc = ((event_added[level]/8)+1)*8;
    MEDDLY_DCASSERT(nalloc > 0);
    MEDDLY_DCASSERT(nalloc > event_added[level]);
    MEDDLY_DCASSERT(nalloc > event_list_alloc[level]);
    event_list[level] = (node_handle*) realloc(event_list[level], nalloc*sizeof(node_handle));
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
  
   printf("\nInside findConfirmedStatesImpl:");
  
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
    MEDDLY::unpacked_node *nr = MEDDLY::unpacked_node::newFromNode(insetF, mdd, false);
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
  for (int i = 1 ; i<num_levels; i++) 
    {
      int levelSize = getInForest()->getLevelSize(i);
      resizeConfirmedArray(i, levelSize);
    }
  
    std::set<node_handle> visited;
    printf("\nfindConfirmedStatesImpl:");
    findConfirmedStatesImpl(const_cast<implicit_relation*>(this),
                      confirmed, confirm_states, set.getNode(), num_levels - 1, visited);
  
}

void MEDDLY::satimpl_opname::implicit_relation::setConfirmedStates(int level, int index)
{
  // For each subevent that affects this level:
  //    (1) call subevent::confirm()
  //    (2) for each level k affected by the subevent,
  //        (a) enlarge variable bound of k to variable bound of k'
  
  printf("\n Confirming index %d at level %d", index, level);
  
  resizeConfirmedArray(level, index+1);
  
  MEDDLY_DCASSERT(confirmed_array_size[level] > index);
  //if (isConfirmedState(level, index)) return false; 
  
  // Get subevents affected by this level, and rebuild them.
  int nSubevents = num_subevents_by_level[level];
  for (int i = 0; i < nSubevents; i++) {
    subevents_by_level[level][i]->confirm(const_cast<implicit_relation&>(*this),
                                          level, index);
   }
  
  const int nEvents = num_events_by_level[level];
  for (int i = 0; i < nEvents; i++) {
    // events_by_level[level][i]->enlargeVariables();
    events_by_level[level][i]->markForRebuilding();
  }
  
  
  
  confirmed[level][index] = true;
  confirm_states[level]++;
  //return true;
}



MEDDLY::satimpl_opname::implicit_relation::~implicit_relation()
{
  //last_in_node_array = 0;
  //impl_unique.clear();
  
  for(int i = 0; i <num_levels; i++) {delete[] event_list[i]; delete[] confirmed[i];}
  delete[] event_list;
  delete[] event_added;
  delete[] event_list_alloc;
  delete[] confirmed;
  delete[] confirm_states;
  delete[] confirmed_array_size;
}


/*node_handle
MEDDLY::satimpl_opname::implicit_relation::isUniqueNode(relation_node* n)
{
  bool is_unique_node = true;
  std::unordered_map<node_handle, relation_node*>::iterator it = impl_unique.begin();
  while(it != impl_unique.end())
    {
    is_unique_node = !((it->second)->equals(n));
    if(is_unique_node==false)
      return (it->second)->getID();
    ++it;
    }
  return 0;
}

node_handle
MEDDLY::satimpl_opname::implicit_relation::registerNode(bool is_event_top, relation_node* n)
{
  
  node_handle nLevel = n->getLevel();

#ifdef DEVELOPMENT_CODE
  node_handle downHandle = n->getDown();
  relation_node* downNode = nodeExists(downHandle);
  node_handle downLevel = downNode->getLevel();
  MEDDLY_DCASSERT( ( ( downNode!=NULL ) && ( nLevel > downLevel ) ) 
                    || 
                   ( downLevel == 0 ) );
#endif

  node_handle n_ID = isUniqueNode(n);
  
  if(n_ID==0) // Add new node
   {
    n_ID  = last_in_node_array + 1;
    std::pair<node_handle, relation_node*> add_node(n_ID,n);
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
*/

void MEDDLY::satimpl_opname::implicit_relation::addTopEvent(node_handle rnh, int level)
{
    resizeEventArray(level);
    event_list[level][event_added[level] - 1] = rnh;
}

void
MEDDLY::satimpl_opname::implicit_relation::show()
{
  node_handle** event_list_copy = (node_handle**)malloc((num_levels)*sizeof(node_handle*));
  if (0==event_list_copy) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  long total_events = 0;
  for(int i = 1;i<num_levels;i++) total_events +=event_added[i];
  for(int i = 1;i<num_levels;i++)
    {
     event_list_copy[i] = (node_handle*)malloc(total_events*sizeof(node_handle));
     if (0==event_list_copy[i]) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }

  for(int i = num_levels-1;i>=1;i--)
    for(int j=0;j<total_events;j++)
      event_list_copy[i][j]=0;
  
  
  int eid = 0;
  for(int i = num_levels-1;i>=1;i--)
    {
     int k = 0;
     std::cout<<"\n [";
     for(int j=0;j<total_events;j++)
      {
      
        if((event_list_copy[i][eid]==0)&&(k<event_added[i]))
          {
            event_list_copy[i][eid] = event_list[i][k];
            relation_node* hold_it = getRelForests()->buildImplicitNode(event_list[i][k]);
            relation_node* hold_down = getRelForests()->buildImplicitNode(hold_it->getDown());
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
  
  for(int i = 0;i<num_levels;i++) delete event_list_copy[i];
  delete[] event_list_copy;
  
}

void MEDDLY::satimpl_opname::implicit_relation::bindExtensibleVariables() {
  //
  // Find the bounds for each extensbile variable
  //
  expert_domain* ed = static_cast<expert_domain*>(outsetF->useDomain());
  
  for (int k = 1; k < num_levels; k++) {
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

/*MEDDLY::node_handle
MEDDLY::satimpl_opname::implicit_relation::buildMxdForest()
{
  
  //Get number of Variables and Events
  int nVars = outsetF->getDomain()->getNumVariables();
  int nEvents = getTotalEvent(nVars);
  
  
  node_handle* event_tops = (node_handle*)malloc((nEvents)*sizeof(node_handle));
  int e = 0;
  
  for(int i = 1 ;i<=nVars;i++)
    {
    int num_events_at_this_level = lengthForLevel(i);
    for(int j = 0;j<num_events_at_this_level;j++)
      event_tops[e++]=arrayForLevel(i)[j];
    }
  
  domain *d = outsetF->useDomain();
  
  forest* mxd = d->createForest(true,forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge nsf(mxd);
  
  dd_edge* monolithic_nsf = new dd_edge(mxd);
  for(int i=0;i<nEvents;i++)
    {
    (*monolithic_nsf) += buildEventMxd(event_tops[i],mxd);
    }
  
  node_handle monolithic_nsf_handle = monolithic_nsf->getNode();
  mxdF = (expert_forest*)mxd;
  
  //for(int i = 0; i<nEvents;i++)
  // {
  // dd_edge nsf_ev(mxd);
  // nsf_ev = buildEventMxd(event_tops[i],mxd);
  // apply(UNION, nsf, nsf_ev, nsf);
  // }
  
  return monolithic_nsf_handle;
}


MEDDLY::dd_edge
MEDDLY::satimpl_opname::implicit_relation::buildEventMxd(node_handle eventTop, forest *mxd)
{
  //mxd is built on a domain obtained from result of saturation
  int nVars = outsetF->getDomain()->getNumVariables();
  //int* sizes = new int[nVars];
  relation_node* Rnode = getRelForest()->buildImplicitNode(eventTop);
  node_handle* rnh_array = (node_handle*)malloc((nVars+1)*sizeof(node_handle));
  // int top_level = Rnode->getLevel();
  
  // domain* d = outsetF->useDomain();
  expert_forest* ef = (expert_forest*) mxd;
  
  //Get relation node handles
  for (int i=nVars; i>=1; i--)
    {
    
      if(Rnode->getLevel()==i)// if variable i is a part of this event
        {
          rnh_array[i] = Rnode->getID(); // keep track of node_handles that are part of this event
          Rnode = getRelForest()->buildImplicitNode(Rnode->getDown()); // move to next variable in the event
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
            Rnode = getRelForest()->buildImplicitNode(rnh_array[i]);
            //Create a new unprimed node for variable i
            MEDDLY_DCASSERT(outsetF->getVariableSize(i)>=Rnode->getPieceSize());
            unpacked_node* UP_var = unpacked_node::newFull(ef, i, Rnode->getPieceSize());
          
            for (int j=0; j<Rnode->getPieceSize(); j++) {
          
              long new_j = confirmed[i][Rnode->getTokenUpdate()[j]]? Rnode->getTokenUpdate()[j] : -2;
              
              if(new_j>=0) 
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
}*/


// *******************************************************************************


std::unordered_map<long,std::vector<node_handle>>
MEDDLY::satimpl_opname::implicit_relation::getListOfNexts(int level, long i, relation_node **R)
{
  std::unordered_map<long,std::vector<node_handle>> jList;
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
  
   int* jset = (int*)malloc(lengthForLevel(level)*sizeof(int));
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

// **************************************************************************************


MEDDLY::satimpl_opname::subevent::subevent(forest* f, int* v, int nv, bool firing, long eval, long fval)
: vars(0), num_vars(nv), root(dd_edge(f)), top(0),
f(static_cast<expert_forest*>(f)), is_firing(firing), enable(eval), fire(fval)
{
  MEDDLY_DCASSERT(f != 0);
  MEDDLY_DCASSERT(v != 0);
  MEDDLY_DCASSERT(nv > 0);
  
  vars = new int[num_vars];
  memcpy(vars, v, unsigned(num_vars) * sizeof(int));
  
  // find top
  top = vars[0];
  for (int i = 1; i < num_vars; i++) {
    if (isLevelAbove(top, vars[i])) top = vars[i];
  }
  
  uses_extensible_variables = false;
  for (int i = 0; i < num_vars; i++) {
    if (this->f->isExtensibleLevel(vars[i])) {
      uses_extensible_variables = true;
      break;
    }
  }
  
  unpminterms = pminterms = 0;
  num_minterms = size_minterms = 0;
  process_minterm_pos = -1;
  processed_minterm_pos = -1;
}

MEDDLY::satimpl_opname::subevent::~subevent()
{
  if (vars) delete [] vars;
  for (int i=0; i<num_minterms; i++) {
    delete[] unpminterms[i];
    delete[] pminterms[i];
  }
  free(unpminterms);
  free(pminterms);
}

void MEDDLY::satimpl_opname::subevent::clearMinterms()
{
  for (int i=0; i<num_minterms; i++) {
    delete[] unpminterms[i];
    delete[] pminterms[i];
  }
  free(unpminterms);
  free(pminterms);
  unpminterms = pminterms = 0;
  num_minterms = 0;
}


void MEDDLY::satimpl_opname::subevent::confirm(implicit_relation& rel, int v, int i) {
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}


bool MEDDLY::satimpl_opname::subevent::addMinterm(const int* from, const int* to)
{
  
  /*
   ostream_output out(std::cout);
   out << "Adding minterm: [";
   for (int i = f->getNumVariables(); i >= 0; i--) {
   out << from[i] << " -> " << to[i] << " , ";
   }
   out << "]\n";
   */
  
  if (num_minterms >= size_minterms) {
    int old_size = size_minterms;
    size_minterms = (0==size_minterms)? 8: MIN(2*size_minterms, 256 + size_minterms);
    unpminterms = (int**) realloc(unpminterms, unsigned(size_minterms) * sizeof(int**));
    pminterms = (int**) realloc(pminterms, unsigned(size_minterms) * sizeof(int**));
    if (0==unpminterms || 0==pminterms) return false; // malloc or realloc failed
    for (int i=old_size; i<size_minterms; i++) {
      unpminterms[i] = 0;
      pminterms[i] = 0;
    }
  }
  if (0==unpminterms[num_minterms]) {
    unpminterms[num_minterms] = new int[f->getNumVariables() + 1];
    MEDDLY_DCASSERT(0==pminterms[num_minterms]);
    pminterms[num_minterms] = new int[f->getNumVariables() + 1];
  }
  // out << "Added minterm: [";
  for (int i = f->getNumVariables(); i >= 0; i--) {
    unpminterms[num_minterms][i] = from[i];
    pminterms[num_minterms][i] = to[i];
    // out << unpminterms[num_minterms][i] << " -> " << pminterms[num_minterms][i] << " , ";
  }
  // out << "]\n";
  expert_domain* d = static_cast<expert_domain*>(f->useDomain());
  for (int i = num_vars - 1; i >= 0; i--) {
    int level = vars[i];
    // expand "to" since the set of unconfirmed local states is always larger
    if (to[level] > 0 && to[level] >= f->getLevelSize(-level)) {
      if (f->isExtensibleLevel(level))
        {
          d->enlargeVariableBound(level, false, -(1+to[level]));
          }
      else
        {
         d->enlargeVariableBound(level, false, 1+to[level]);
        }
    }
  }
  num_minterms++;
  process_minterm_pos +=1;
  return true;
}

void MEDDLY::satimpl_opname::subevent::buildRoot() {
  //printf("\n num_minterms in this se = %d, to be done = %d ",num_minterms,process_minterm_pos-processed_minterm_pos );
  if (0 == num_minterms) return;
  if (processed_minterm_pos == process_minterm_pos) return;
  //if (1 == num_vars) return ;
  
   /*ostream_output out(std::cout);
   out << "\nBuilding subevent from " << num_minterms << " minterms\n";
   out << "New minterm is added at position " << process_minterm_pos << " minterms\n";
   for (int i = 0; i < num_minterms; i++) {
   out << "minterm[" << i << "]: [ ";
   for (int j = f->getNumVariables(); j >= 0; j--) {
   out << unpminterms[i][j] << " -> " << pminterms[i][j] << ", ";
   }
   out << " ]\n";
   }*/
   
  std::vector<std::vector<int>> terms;
  std::vector<int> pterms;
  std::vector<int> unpterms;
  
 /* for(int i=1;i<=f->getNumVariables();i++)
    {
      unpterms.push_back(unpminterms[process_minterm_pos][i]);
      pterms.push_back(pminterms[process_minterm_pos][i]);
    }
  terms.push_back(unpterms);
  terms.push_back(pterms);
  */
  
  node_handle rnh = 0;
  while(processed_minterm_pos!=process_minterm_pos)
    {
      for(int i=1;i<=f->getNumVariables();i++)
        {
        unpterms.push_back(unpminterms[processed_minterm_pos+1][i]);
        pterms.push_back(pminterms[processed_minterm_pos+1][i]);
        }
      terms.push_back(unpterms);
      terms.push_back(pterms);
    
      rnh = f->unionOneMinterm(root.getNode(), terms);
      if(rnh!=root.getNode()) root.set(rnh);
      processed_minterm_pos +=1;
    
      unpterms.clear();
      pterms.clear();
      terms.clear();
    
    }
  
  #if 0
  if (usesExtensibleVariables()) {
    dd_edge sum(root);
    f->createEdge(unpminterms, pminterms, num_minterms, sum);
    num_minterms = 0;
    root += sum;
  } else {
    f->createEdge(unpminterms, pminterms, num_minterms, root);
  }
  #endif
  
  
  root_handle = root.getNode();
  printf("\n My root handle is %d",root_handle);
  // out << "\nEquivalent event: " << root.getNode() << "\n";
  // out << "Result: ";
  // root.show(out, 2);
}


void MEDDLY::satimpl_opname::subevent::showInfo(output& out) const {
  int num_levels = f->getDomain()->getNumVariables();
  for (int i = 0; i < num_minterms; i++) {
    out << "minterm[" << i << "]: ";
    for (int lvl = num_levels; lvl > 0; lvl--) {
      out << unpminterms[i][lvl] << " -> " << pminterms[i][lvl] << ", ";
    }
    out << "]\n";
  }
  root.show(out, 2);
}

long MEDDLY::satimpl_opname::subevent::mintermMemoryUsage() const {
  long n_minterms = 0L;
  for (int i = 0; i < size_minterms; i++) {
    if (unpminterms[i] != 0) n_minterms++;
  }
  return long(n_minterms * 2) * long(f->getDomain()->getNumVariables()) * long(sizeof(int));
}

// ============================================================

#if 1
MEDDLY::satimpl_opname::event::event(subevent** p, int np, relation_node** r, int nr)
{
  if (np == 0 && nr == 0)
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
    
  f = np == 0 ? r[0]->getForest() : p[0]->getForest();
  
  for (int i=1; i<np; i++) {
    if (p[i]->getForest() != f) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }
  
  for (int i=1; i<nr; i++) {
    if (r[i]->getForest() != f) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }
  
  /*if (p == 0 || np <= 0 || p[0]->getForest() == 0)
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  f = p[0]->getForest();
  for (int i=1; i<np; i++) {
    if (p[i]->getForest() != f) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }*/
  
  top = 0;
  
  num_subevents = np;
  if(np>0)
    {
      subevents = new subevent*[np];
      std::vector<int> sorted_by_top_se;
      for (int i = 0; i < np; i++) sorted_by_top_se.push_back(p[i]->getTop());
      std::sort(sorted_by_top_se.begin(),sorted_by_top_se.end());
      for (int i = 0; i < np; i++) 
        {
          int pos = -1;
          for (int j = 0; j < np; j++) 
            if(sorted_by_top_se[j] == p[i]->getTop())
              pos = j;
          subevents[pos] = p[i];
        }
      top = sorted_by_top_se[np-1];
    }
  if(np==0)
    subevents = NULL;
  
  num_relnodes = nr;
  if(nr>0)
    {
      relnodes = new relation_node*[nr];
      std::vector<int> sorted_by_top_rn;
      for (int i = 0; i < nr; i++) sorted_by_top_rn.push_back(r[i]->getLevel());
      std::sort(sorted_by_top_rn.begin(),sorted_by_top_rn.end());
      for (int i = 0; i < nr; i++) 
        {
        int pos = -1;
        for (int j = 0; j < nr; j++) 
          if(sorted_by_top_rn[j] == r[i]->getLevel())
            pos = j;
        relnodes[pos] = r[i];
        }
    top = sorted_by_top_rn[nr-1] > top ?  sorted_by_top_rn[nr-1]:top;
    }
  if(nr==0)
    relnodes = NULL;
  
  
  printf("\n top level = %d", top);
  /*top = np>0? p[0]->getTop() : r[0]->getLevel(); 
  for (int i = 1; i < np; i++) {
    if (top < p[i]->getTop()) top = p[i]->getTop();
  }  
  for (int i = 0; i < nr; i++) {
    if (top < r[i]->getLevel()) top = r[i]->getLevel();
  }*/
  
  
  num_components = num_subevents + num_relnodes;
  
  // Find the variable that effect this event from the list of subevents.
  // TODO:
  // Not efficient. p[i] is a sorted list of integers.
  // Should be able to insert in O(n) time
  // where n is the sum(p[i]->getNumVars).
#if 0
  bool all_firing_subevents = true;
#endif
#ifdef DEVELOPMENT_CODE
  bool all_enabling_subevents = true;
#endif
  std::set<int> rVars;
  std::set<int> sVars;
  std::set<int> firingVars;
  
  //each relation node depends only on a single variable
    if(nr>0)
      for (int i = 0; i < nr; i++)
        rVars.insert(r[i]->getLevel()); 
  
  if(np>0)
    {
      for (int i = 0; i < np; i++) {
        const int* subeventVars = p[i]->getVars();
        sVars.insert(subeventVars, subeventVars+p[i]->getNumVars());
        if (p[i]->isFiring()) {
    #ifdef DEVELOPMENT_CODE
          all_enabling_subevents = false;
    #endif
          firingVars.insert(subeventVars, subeventVars+p[i]->getNumVars());
        } else {
    #if 0
          all_firing_subevents = false;
    #endif
        }
      }
    }
  
  MEDDLY_DCASSERT(all_enabling_subevents || !firingVars.empty());
  
#if 0
  is_disabled = (all_enabling_subevents || all_firing_subevents);
#else
  is_disabled = false;
#endif
  num_vars = sVars.size() + rVars.size();
  vars = new int[num_vars];
  int* curr = &vars[0];
  std::set<int>::iterator rit=rVars.begin();
  std::set<int>::iterator sit=sVars.begin();
  while(rit!=rVars.end() && sit!=sVars.end()) {
    if(*rit < *sit)
      *curr++ = *rit++;
    else  
      *curr++ = *sit++;
  }
  while(rit!=rVars.end())
    *curr++ = *rit++;
  while(sit!=sVars.end())
    *curr++ = *sit++;
  
  num_firing_vars = firingVars.size();
  firing_vars = new int[num_firing_vars];
  curr = &firing_vars[0];
  for (std::set<int>::iterator it=firingVars.begin(); it!=firingVars.end(); ) {
    *curr++ = *it++;
  }
  
  
  num_rel_vars = rVars.size();
  relNode_vars = new int[num_rel_vars];
  curr = &relNode_vars[0];
  for (std::set<int>::iterator it=rVars.begin(); it!=rVars.end(); ) {
    *curr++ = *it++;
  }
  
  // Create the implicit nodes already
  // Make sure the nodes are sorted by level bottom to top
  node_handle nh_next = -1;
  for (int i = 0; i < nr; i++)
    {
      relnodes[i]->setDown(nh_next);
      nh_next = f->createRelationNode(relnodes[i]);
      relnodes[i]->setID(nh_next);
    }
  
  rebuild();
  root = dd_edge(f);
  event_mask = dd_edge(f);
  event_mask_from_minterm = 0;
  event_mask_to_minterm = 0;
  needs_rebuilding = is_disabled? false: true;
  
  printf("\n This event is done building");
  printf("\n");
}
#endif

/****************/
 #if 0
MEDDLY::satimpl_opname::event::event(subevent** p, int np)
{
  if (np == 0)
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  
  f = p[0]->getForest();
  
  for (int i=1; i<np; i++) {
    if (p[i]->getForest() != f) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }
  
  /*if (p == 0 || np <= 0 || p[0]->getForest() == 0)
   throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
   f = p[0]->getForest();
   for (int i=1; i<np; i++) {
   if (p[i]->getForest() != f) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
   }*/
  
  top = 0;
  num_components = np;
  
  
  std::vector<subevent*> dep_se;
  std::vector<subevent*> indep_se;
  
  for (int i=0; i<np; i++) {
    if(p[i]->isImplicit()) indep_se.push_back(p[i]);
    else dep_se.push_back(p[i]);
    }
  
  num_relnodes = indep_se.size();
  num_subevents = dep_se.size();
  
  
  // get the top of event
  if(num_subevents>0)
    {
    subevents = new subevent*[num_subevents];
    std::vector<int> sorted_by_top_se;
    for (int i = 0; i < num_subevents; i++) sorted_by_top_se.push_back(dep_se[i]->getTop());
    std::sort(sorted_by_top_se.begin(),sorted_by_top_se.end());
    for (int i = 0; i < num_subevents; i++) 
      {
      int pos = -1;
      for (int j = 0; j < num_subevents; j++) 
        if(sorted_by_top_se[j] == dep_se[i]->getTop())
          pos = j;
      subevents[pos] = dep_se[i];
      }
    top = sorted_by_top_se[num_subevents-1];
    }
  if(num_subevents==0)
    subevents = NULL;
  
  if(num_relnodes>0)
    {
    relnodes = new relation_node*[num_relnodes];
    std::vector<int> sorted_by_top_rn;
    for (int i = 0; i < num_relnodes; i++) sorted_by_top_rn.push_back(indep_se[i]->getTop());
    std::sort(sorted_by_top_rn.begin(),sorted_by_top_rn.end());
    for (int i = 0; i < num_relnodes; i++) 
      {
      int pos = -1;
      for (int j = 0; j < num_relnodes; j++) 
        if(sorted_by_top_rn[j] == indep_se[i]->getTop())
          pos = j;
      relnodes[pos] = new relation_node(14401,indep_se[i]->getForest(),indep_se[i]->getTop(),indep_se[i]->getDown(), indep_se[i]->getEnable(),indep_se[i]->getFire());
      }
    top = sorted_by_top_rn[num_relnodes-1] > top ?  sorted_by_top_rn[num_relnodes-1]:top;
    }
  if(num_relnodes==0)
    relnodes = NULL;
  
  
  /*top = np>0? p[0]->getTop() : r[0]->getLevel(); 
   for (int i = 1; i < np; i++) {
   if (top < p[i]->getTop()) top = p[i]->getTop();
   }  
   for (int i = 0; i < nr; i++) {
   if (top < r[i]->getLevel()) top = r[i]->getLevel();
   }*/
  
  
  
  // Find the variables that effect this event from the list of subevents.
  // TODO:
  // Not efficient. p[i] is a sorted list of integers.
  // Should be able to insert in O(n) time
  // where n is the sum(p[i]->getNumVars).
#if 0
  bool all_firing_subevents = true;
#endif
#ifdef DEVELOPMENT_CODE
  bool all_enabling_subevents = true;
#endif
  std::set<int> rVars;
  std::set<int> sVars;
  std::set<int> firingVars;
  
  //each relation node depends only on a single variable
  if(num_relnodes>0)
    for (int i = 0; i < num_relnodes; i++)
      rVars.insert(indep_se[i]->getTop()); 
  
  if(num_subevents>0)
    {
    for (int i = 0; i < num_subevents; i++) {
      const int* subeventVars = dep_se[i]->getVars();
      sVars.insert(subeventVars, subeventVars+dep_se[i]->getNumVars());
      if (dep_se[i]->isFiring()) {
#ifdef DEVELOPMENT_CODE
        all_enabling_subevents = false;
#endif
        firingVars.insert(subeventVars, subeventVars+dep_se[i]->getNumVars());
      } else {
#if 0
        all_firing_subevents = false;
#endif
      }
    }
    }
  
  MEDDLY_DCASSERT(all_enabling_subevents || !firingVars.empty());
  
#if 0
  is_disabled = (all_enabling_subevents || all_firing_subevents);
#else
  is_disabled = false;
#endif
  num_vars = sVars.size() + rVars.size();
  vars = new int[num_vars];
  int* curr = &vars[0];
  std::set<int>::iterator rit=rVars.begin();
  std::set<int>::iterator sit=sVars.begin();
  while(rit!=rVars.end() && sit!=sVars.end()) {
    if(*rit < *sit)
      *curr++ = *rit++;
    else  
      *curr++ = *sit++;
  }
  while(rit!=rVars.end())
    *curr++ = *rit++;
  while(sit!=sVars.end())
    *curr++ = *sit++;
  
  num_firing_vars = firingVars.size();
  firing_vars = new int[num_firing_vars];
  curr = &firing_vars[0];
  for (std::set<int>::iterator it=firingVars.begin(); it!=firingVars.end(); ) {
    *curr++ = *it++;
  }
  
  
  num_rel_vars = rVars.size();
  relNode_vars = new int[num_rel_vars];
  curr = &relNode_vars[0];
  for (std::set<int>::iterator it=rVars.begin(); it!=rVars.end(); ) {
    *curr++ = *it++;
  }
  
  // Create the implicit nodes already
  // Make sure the nodes are sorted by level bottom to top
  node_handle nh_next = -1;
  for (int i = 0; i < num_relnodes; i++)
    {
    relnodes[i]->setDown(nh_next);
    nh_next = f->createRelationNode(relnodes[i]);
    relnodes[i]->setID(nh_next);
    printf("\n Bult implicit node %d",nh_next);
    }
  
  rebuild();
  root = dd_edge(f);
  event_mask = dd_edge(f);
  event_mask_from_minterm = 0;
  event_mask_to_minterm = 0;
  needs_rebuilding = is_disabled? false: true;
}
#endif
/****************/


MEDDLY::satimpl_opname::event::~event()
{
  for (int i=0; i<num_subevents; i++) delete subevents[i];
  for (int i=0; i<num_relnodes; i++) delete relnodes[i];
  delete[] subevents;
  delete[] relnodes;
  delete[] vars;
  delete[] relNode_vars;
  delete[] firing_vars;
  delete[] event_mask_from_minterm;
  delete[] event_mask_to_minterm;
}


void MEDDLY::satimpl_opname::event::buildEventMask()
{
  MEDDLY_DCASSERT(num_subevents > 0);
  MEDDLY_DCASSERT(f);
  
  if (0 == event_mask_from_minterm) {
    const size_t minterm_size = size_t(f->getNumVariables()+1);
    event_mask_from_minterm = new int[minterm_size];
    event_mask_to_minterm = new int[minterm_size];
    
    for (unsigned i = 0; i < minterm_size; i++) {
      event_mask_from_minterm[i] = MEDDLY::DONT_CARE;
      event_mask_to_minterm[i] = MEDDLY::DONT_CHANGE;
    }
    
    for (int i = 0; i < num_firing_vars; i++) {
      event_mask_to_minterm[firing_vars[i]] = MEDDLY::DONT_CARE;
    }
  }
  
   f->createEdge(&event_mask_from_minterm, &event_mask_to_minterm, 1, event_mask);
  // See how to build the event mask.
  /*
   node_handle next = value2handle(0);
   bool impl_flag = false;
   for( int i = 1; i <= f->getNumVariables(); i++)
    {
        // if this variable can be implicit relation, the event mask has it.  
          for( int j = 0; j < num_rel_vars; j++)
           if(relnodes[j]->getLevel() == i)
                {
                    next = f->createRelationNode(relnodes[j]);
                    impl_flag = true;
                    break;
                }
    
            if(!impl_flag) {
              
            node_handle nextpr;
            if (DONT_CARE == event_mask_to_minterm[i]) {
            if (f->isFullyReduced()) {
              // DO NOTHING
              nextpr = next;
            } else {
              // build redundant node
              unpacked_node* nb = 0;
              if (f->isExtensibleLevel(i)) {
                nb = unpacked_node::newFull(f, -i, 1);
                nb->d_ref(0) = next;
                // link count should be unchanged
                nb->markAsExtensible();
              } else {
                int sz = f->getLevelSize(-i);
                nb = unpacked_node::newFull(f, -i, sz);
                for (int v=0; v<sz; v++) {
                  nb->d_ref(v) = f->linkNode(next);
                }
                f->unlinkNode(next);
              }
              nextpr = f->createReducedNode(-1, nb);
            }
          } else if (DONT_CHANGE == event_mask_to_minterm[i]) {
            //
            // Identity node
            //
            if(DONT_CARE == event_mask_from_minterm[i]){
              if (f->isIdentityReduced()) continue;
              next = makeIdentityEdgeForDontCareDontChange(i, next);
              continue;
            }
            
            MEDDLY_DCASSERT(event_mask_from_minterm[i]>=0);
            
            if (f->isIdentityReduced()) {
              // DO NOTHING
              nextpr = next;
            }
            else if(f->isQuasiReduced() && f->getTransparentNode()!=value2handle(0)){
              unsigned sz = unsigned(f->getLevelSize(-i));
              unpacked_node* nbp = unpacked_node::newFull(f, -i, sz);
              node_handle zero=makeOpaqueZeroNodeAtLevel(i-1);
              for(unsigned v=0; v<sz; v++){
                nbp->d_ref(v)=(v==event_mask_from_minterm[i] ? f->linkNode(next) : f->linkNode(zero));
              }
              f->unlinkNode(zero);
              
              nextpr = f->createReducedNode(event_mask_from_minterm[i], nbp);
            }
            else {
              unpacked_node* nbp = unpacked_node::newSparse(f, -i, 1);
              nbp->i_ref(0) = event_mask_from_minterm[i];
              nbp->d_ref(0) = next;
              // link count should be unchanged
              
              nextpr = f->createReducedNode(event_mask_from_minterm[i], nbp);
            }
          }
          }
    
          impl_flag = false;
    }
  */  
  
  
#ifdef DEBUG_EVENT_MASK
  printf("event_mask: %d\n" , event_mask.getNode());
  ostream_output out(std::cout);
  event_mask.show(out, 2);
#endif
}


bool MEDDLY::satimpl_opname::event::rebuild()
{
  //MEDDLY_DCASSERT(num_subevents > 0);
  /*if (is_disabled) return false;*/
  if(num_subevents > 0)
    {
      if (!needs_rebuilding) return false;
      needs_rebuilding = false;
    }
  
  // An event is a conjunction of sub-events (or sub-functions).
  printf("\n I am in event %d", this);
  for (int i = 0; i < num_subevents; i++)
    {
      printf("\n Rebuild is about to happen for sub-event %d",i);
      subevents[i]->buildRoot();
    }
     
  printf("\n See the chain of rel_nodes:%d nodes",num_relnodes);
  if(num_relnodes>0)
    {
      int kk = 0;
      relation_node *temp = NULL;
      node_handle nnh = relnodes[num_relnodes-1]->getID();
      while(nnh!=-1)
        {
          temp = f->buildImplicitNode(nnh);
          printf("%d<%d>->",temp->getID(),temp->getLevel());
          nnh = temp->getDown();
        }
     }
  
  
  //replace the existing subevent's previous handles
  level_component.clear();
   for (int i = 0; i < num_subevents; i++) {
     int level_of_subevent = subevents[i]->getTop();
     
     if(level_component.find(level_of_subevent)!=level_component.end())
         level_component.at(level_of_subevent).insert(subevents[i]->getRootHandle());
     else 
       {
         std::set<node_handle> handles_at_this_top;
         handles_at_this_top.insert(subevents[i]->getRootHandle());
         level_component.insert(std::make_pair(level_of_subevent,handles_at_this_top));
       }
    }
  
   // Modify it to only hold the topmost relNode!
   for (int i = 0; i < num_relnodes; i++)
     {
        int level_of_relNode = relnodes[i]->getLevel();
       
       if(level_component.find(level_of_relNode)!=level_component.end())
         level_component.at(level_of_relNode).insert(relnodes[i]->getID());
       else 
         {
          std::set<node_handle> handles_at_this_top;
          handles_at_this_top.insert(relnodes[i]->getID());
          level_component.insert(std::make_pair(level_of_relNode,handles_at_this_top));
         }
     } 
  
  
  
  //buildEventMask();
  
  /*dd_edge e(event_mask);
  for (int i = 0; i < num_subevents; i++) {
    e *= subevents[i]->getRoot();
  }*/
  
  /*
   if (e.getNode() == 0) {
   ostream_output out(std::cout);
   f->useDomain()->showInfo(out);
   out << "subevent: " << event_mask.getNode() << "\n";
   event_mask.show(out, 2);
   for (int i = 0; i < num_subevents; i++) {
   out << "subevent: " << subevents[i]->getRoot().getNode() << "\n";
   subevents[i]->getRoot().show(out, 2);
   }
   }*/
   ostream_output out(std::cout);
   for (int i = 0; i < num_subevents; i++) {
   out << "subevent: " << subevents[i]->getRoot().getNode() << "\n";
   if(subevents[i]->getRoot().getNode()) subevents[i]->getRoot().show(out, 2);
   }
   //e.show(out, 2);
   
  
  if (level_component[top] == root_handle) return false;
  root_handle = level_component[top];
  return true;
}

int MEDDLY::satimpl_opname::event::downLevel(int level) const{ 
  for(int i = num_vars; i>0; i--)
    {
      if(vars[i] == level)
        return vars[i-1];
    }
  return 0;
}

void MEDDLY::satimpl_opname::event::enlargeVariables()
{
  expert_domain* ed = static_cast<expert_forest*>(f)->useExpertDomain();
  for (int i = 0; i < num_vars; i++) {
    int unprimed = ABS(vars[i]);
    int primed = -unprimed;
    int unprimedSize = f->getLevelSize(unprimed);
    int primedSize = f->getLevelSize(primed);
    if (unprimedSize < primedSize) {
      expert_variable* vh = ed->getExpertVar(unprimed);
      if (vh->isExtensible())
        vh->enlargeBound(false, -primedSize);
      else
        vh->enlargeBound(false, primedSize);
    }
    MEDDLY_DCASSERT(f->getLevelSize(unprimed) == f->getLevelSize(primed));
  }
}


void MEDDLY::satimpl_opname::event::showInfo(output& out) const {
  for (int i = 0; i < num_subevents; i++) {
    out << "subevent " << i << "\n";
    subevents[i]->showInfo(out);
  }
}

long MEDDLY::satimpl_opname::event::mintermMemoryUsage() const {
  long usage = 0;
  for (int i = 0; i < num_subevents; i++) {
    usage += subevents[i]->mintermMemoryUsage();
  }
  return usage;
}


// ******************************************************************


// ******************************************************************
// *                                                                *
// *               saturation_impl_by_events_opname  class          *
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
  inline compute_table::entry_key*
  findSaturateResult(node_handle a, int level, node_handle& b) {
    compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeN(a);
    if (argF->isFullyReduced()) CTsrch->writeI(level);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    b = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
  }
  inline void recycleCTKey(compute_table::entry_key* CTsrch) {
    CT0->recycle(CTsrch);
  }
  inline node_handle saveSaturateResult(compute_table::entry_key* Key,
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
// *            common_impl_dfs_by_events_mt class                 *
// *                                                                *
// ******************************************************************

class MEDDLY::common_impl_dfs_by_events_mt : public specialized_operation {
public:
  common_impl_dfs_by_events_mt(const satimpl_opname* opcode,
                               satimpl_opname::implicit_relation* rel);
  virtual ~common_impl_dfs_by_events_mt();
  
  virtual void compute(const dd_edge& a, dd_edge &c);
  /*virtual bool isReachable(const dd_edge& a, const dd_edge& constraint);*/
  virtual void saturateHelper(unpacked_node& mdd) = 0;
  // for detecting reachable state in constraint
  /*virtual bool saturateHelper(unpacked_node& mdd, node_handle constraint) = 0;*/
  
protected:
  inline compute_table::entry_key*
  findResult(node_handle a, std::vector<node_handle> b, node_handle &c)
  {
    // needs to be searched in that CT where the forest matches!
    compute_table::entry_key* CTsrch = 0;
    CTsrch = CT0->useEntryKey(etype[0], b.size());
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeN(a);
    for(auto it = b.begin(); it!=b.end(); it++ )
       CTsrch->writeN(*it);
    //CTsrch->writeL(b);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    c = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
  }
  
  inline void recycleCTKey(compute_table::entry_key* CTsrch) {
    CT0->recycle(CTsrch);
  }
  inline node_handle saveResult(compute_table::entry_key* Key,
                                node_handle a, std::vector<node_handle> b, node_handle c)
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
  node_handle recFire(node_handle mdd, std::map<int,std::set<node_handle>> seHandles); //mxd, satimpl_opname::event* e);
  //MEDDLY::node_handle recFireSet(node_handle mdd, std::vector<node_handle> mxd);

  // for reachable state in constraint detection
  /*bool saturateHelper(
      unpacked_node& nb,
      node_handle constraint);
  bool recFire(
       node_handle mdd,
       node_handle mxd,
       node_handle constraint,
       node_handle& result);*/
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
  
  // Obtain handles to events for this level
  const int level = nb.getLevel();
  satimpl_opname::event** event_handles = rel->arrayForLevel(level); 
  
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
    
    // for each event
    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
      
      // rebuild the marking-dependent subevents, if needed
      // change the if-condition
      event_handles[ei]->rebuild();
      
      printf("\n COmpleted rebuild");
      // build a set of pairs, where each element contains: {level,node_handle} 
      // an updated set is passed to next recursion based on the level processed in current recursion
      std::map<int, std::set<node_handle>> seHandles = event_handles[ei]->getComponents();
      const int* eVars = event_handles[ei]->getVars();
      int num_eVars = event_handles[ei]->getNumVars();
      
      int jc = 0;
      std::vector<int> jvector;
      std::vector<node_handle> dvector;
      std::map<int, std::vector<node_handle>> j_d_map;
      
      // intersectio of the subevents
       printf("\n Intersection of subevents at level %d",level);
      for(auto nse_it = seHandles[level].begin(); nse_it != seHandles[level].end(); nse_it++) {
         printf("\n Handle:%d",*nse_it);
        
        if(*nse_it==0) continue;
        bool imFlag = arg2F->isImplicit(*nse_it);
        printf("\n Am i implicit?:%d",imFlag);
        
        if(imFlag)
          {
          relation_node* Ru = arg2F->buildImplicitNode(*nse_it);
          int j = Ru->nextOf(i);
          
          if(j_d_map.find(j)!=j_d_map.end())
            j_d_map.at(j).push_back(Ru->getDown());
          else 
            {
              std::vector<node_handle> handles_down_j;
              handles_down_j.push_back(Ru->getDown());
              j_d_map.insert(std::make_pair(j,handles_down_j));
            }
          
          } else {
            printf("\n I am about to build the mxd Node for %d on level %d", *nse_it,arg2F->getNodeLevel(*nse_it));
            unpacked_node* Ru = unpacked_node::useUnpackedNode();
            Ru->initFromNode(arg2F, *nse_it, true);
             printf("\n I did so");
            node_handle ei_p = Ru->d(i);
            if (0 == ei_p) continue;
            unpacked_node* Rp = unpacked_node::useUnpackedNode();
            Rp->initFromNode(arg2F, ei_p, false);
            for(int jz = 0; jz < Rp->getNNZs(); jz++)
              {
                jvector.push_back(Rp->i(jz));
                dvector.push_back(Rp->d(jz));
              if(j_d_map.find(Rp->i(jz))!=j_d_map.end())
                j_d_map.at(Rp->i(jz)).push_back(Rp->d(jz));
              else 
                {
                  std::vector<node_handle> handles_down_j;
                  handles_down_j.push_back(Rp->d(jz));
                  j_d_map.insert(std::make_pair(Rp->i(jz),handles_down_j));
                }
              }
            printf("\n I built the mxd Node");
          }
        
      }
      
      std::map<int, std::set<node_handle>> j_dset_map;
      
      printf("\n Intresetion vector size  %d", seHandles[level].size());
      for(auto jit = j_d_map.begin();jit!=j_d_map.end();jit++)
        {
          if(jit->second.size()==seHandles[level].size())
            { 
              std::set<node_handle> set_of_handles_down_j(jit->second.begin(), jit->second.end());
              j_dset_map.insert(std::make_pair(jit->first,set_of_handles_down_j));
              printf("\n jd_map: [%d]:",jit->first);
              for(auto kit = j_dset_map[jit->first].begin();kit!= j_dset_map[jit->first].end();kit++)
              printf("\n %d,",*kit);
            }
         }
      
      
      for(auto jit = j_dset_map.begin(); jit!=j_dset_map.end(); jit++)
        { 
          int j = jit->first;
          if(j==-1) continue;
          if (j < nb.getSize() && -1==nb.d(j)) continue; // nothing can be added to this set
          
          std::map<int, std::set<node_handle>> seHandlesForRecursion = seHandles;
          for(auto jdit = jit->second.begin(); jdit != jit->second.end(); jdit++)
            {
              int level_of_downj = arg2F->getNodeLevel(*jdit);
              if( seHandlesForRecursion.find(level_of_downj) != seHandlesForRecursion.end())
                  seHandlesForRecursion.at(level_of_downj).insert(*jdit);
              else
                {
                  std::set<node_handle> handles_at_level;
                  handles_at_level.insert(*jdit);
                  seHandlesForRecursion.insert(std::make_pair(level_of_downj, handles_at_level));
                }
            }
          seHandlesForRecursion.erase(level); 
         
          node_handle rec = recFire(nb.d(i), seHandlesForRecursion);
          
          
          if (rec == 0) continue;
          
          //confirm local state
          printf("\n Confirming j=%d from i=%d",j,i);
          if(!rel->isConfirmedState(level,j))
            rel->setConfirmedStates(level, j);
          
          if(j>=nb.getSize())
            {
            int new_var_bound = resF->isExtensibleLevel(nb.getLevel())? -(j+1): (j+1);
            dm->enlargeVariableBound(nb.getLevel(), false, new_var_bound);
            int oldSize = nb.getSize();
            nb.resize(j+1);
            while(oldSize < nb.getSize()) { nb.d_ref(oldSize++) = 0; }
            queue->resize(nb.getSize());
             printf("\n Resizing  queue to %d", nb.getSize());
            }  
          if (rec == nb.d(j)) {
            resF->unlinkNode(rec);
            continue;
          }
          
          bool updated = true;
          
          if (0 == nb.d(j)) {
            nb.d_ref(j) = rec;
            printf("\n Assigning %dth child of node %d as %d", j,&nb,rec);
          }
          else if (rec == -1) {
            resF->unlinkNode(nb.d(j));
            nb.d_ref(j) = -1;
            printf("\n Assigning %dth child of node %d as %d", j,&nb,nb.d(j));
          }
          else {
            nbdj.set(nb.d(j));  // clobber
            newst.set(rec);     // clobber
            mddUnion->computeTemp(nbdj, newst, nbdj);
            updated = (nbdj.getNode() != nb.d(j));
            nb.set_d(j, nbdj);
            }
          
          if (updated) {printf("\n Added index %d for level %d", j, level);queue->add(j);}
          
        } //for all j's
      printf("\n No more j's to look at for this ei %d", ei);
    } // for all events, ei
  }// more indexes to explore
  
  //delete[] Ru;
  recycle(queue);
}

// Same as post-image, except we saturate before reducing.
MEDDLY::node_handle MEDDLY::forwd_impl_dfs_by_events_mt::recFire(
                                                                 MEDDLY::node_handle mdd, std::map<int, std::set<node_handle>> seHandles)
{
  // current mxd handle is at highest level from this set
  std::map<int,std::set<node_handle>>::reverse_iterator rit = seHandles.rbegin();
  std::set<node_handle> set_mxd = rit->second;
  printf("\n {}{}{}{}{} setmxd :");
  for(auto iit = set_mxd.begin(); iit != set_mxd.end();iit++ )
    printf("%d,",*iit);
    
  node_handle mxd = 0;
  if(set_mxd.size() == 1)
    {
      auto sit = set_mxd.begin();
      mxd = *sit;
  
     // termination conditions
     if (mxd == 0 || mdd == 0) return 0;
  
      if (mxd == -1) {
        if (arg1F->isTerminalNode(mdd))
        {
          return resF->handleForValue(1);
        }
        // mxd is identity
        if (arg1F == resF)
         return resF->linkNode(mdd);
      }
    }
  
  // keep track of compute_table entries
  std::vector<node_handle> mxd_vector;
  
  printf("\n ComputeTableEntry : <mdd:%d,mxd(",mdd);
  for(auto it = seHandles.begin(); it!= seHandles.end(); it++)
    for(auto hit = it->second.begin(); hit!= it->second.end(); hit++)
        mxd_vector.push_back(*hit);
    
  printf(")>");
 
  // check the cache
  FILE_output meddlyout(stdout);
  node_handle result = 0;
  compute_table::entry_key* Key = findResult(mdd, mxd_vector, result); // need to see if i might need to store forest into compute_table or maintain separate compute tables for each event forest
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
  const int mxdLevel = rit->first; //arg2F->getNodeLevel(mxd);
  const int rLevel = MAX(mxdLevel, mddLevel);
  int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);
  expert_domain* dm = static_cast<expert_domain*>(resF->useDomain());
  
  dd_edge nbdj(resF), newst(resF);
  
  
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
    printf("\n Skipped levels in the MXD");
    for (int i=0; i<rSize; i++) {
       printf("\n recFire(%d)",A->d(i));
      nb->d_ref(i) = recFire(A->d(i), seHandles);
    }
    
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(mxdLevel >= mddLevel);
    
    // clear out result (important!)    
    for (unsigned i=0; i<rSize; i++) nb->d_ref(i) = 0;
    
    std::map< std::pair<int, int>, std::vector<node_handle> > from_to_handle;
    
    for(auto sit = set_mxd.begin(); sit!=set_mxd.end(); sit++) // set of mxds at level k
    {
      mxd = *sit;
      int jc = 0;
      std::vector<int> jvector;
      std::vector<node_handle> dvector;
      bool imFlag = arg2F->isImplicit(mxd);
      int row_size = rSize;
      
      if(imFlag) {
        relation_node* relNode = arg2F->buildImplicitNode(mxd);
        row_size = rSize; // Get the update array size;
      } 
      else {
        unpacked_node *Ru = unpacked_node::useUnpackedNode();
        if (mxdLevel < 0) {
          Ru->initRedundant(arg2F, rLevel, mxd, false);
        } else {
          Ru->initFromNode(arg2F, mxd, false);
        }
        row_size = Ru->getNNZs();
      }
      
    
      // for intersection to happen : For each i, check across the set_mxd what are the common j's. Then use those j's only. 
      
      for (int iz=0; iz<row_size; iz++) { // for each i -> get common j's downhandle
        jvector.clear();
        dvector.clear();
        int i = iz;
        // Initialize mxd readers, note we might skip the unprimed level
        if(imFlag) {
          
          relation_node* relNode = arg2F->buildImplicitNode(mxd);
          jc = 1;
          int j = relNode->nextOf(i);
          jvector.push_back(relNode->nextOf(i));
          dvector.push_back(relNode->getDown());
          
          std::pair <int,int> from_to(i,j);
          
          if( from_to_handle.find(from_to) != from_to_handle.end())
            from_to_handle.at(from_to).push_back(relNode->getDown());
          else
            {
            std::vector<node_handle> handles_at_change;
            handles_at_change.push_back(relNode->getDown());
            from_to_handle.insert(std::make_pair(from_to, handles_at_change));
            }
          
        } else  { 
          
          unpacked_node *Ru = unpacked_node::useUnpackedNode();
          unpacked_node *Rp = unpacked_node::useUnpackedNode();
          if (mxdLevel < 0) {
            Ru->initRedundant(arg2F, rLevel, mxd, false);
          } else {
            Ru->initFromNode(arg2F, mxd, false);
          }
          
          i = Ru->i(iz);
          
          if (0==A->d(i)) continue; 
          const node_handle pnode = Ru->d(iz);
          if (isLevelAbove(-rLevel, arg2F->getNodeLevel(pnode))) {
            Rp->initIdentity(arg2F, rLevel, i, pnode, false);
          } else {
            Rp->initFromNode(arg2F, pnode, false);
          }
          jc = Rp->getNNZs();
          // loop over mxd "columns"
          for(int jz = 0; jz < jc; jz++)
            {
              jvector.push_back(Rp->i(jz));
              dvector.push_back(Rp->d(jz));
              std::pair <int,int> from_to(i,Rp->i(jz));
            
              if( from_to_handle.find(from_to) != from_to_handle.end())
                from_to_handle.at(from_to).push_back(Rp->d(jz));
              else
              {
                std::vector<node_handle> handles_at_change;
                handles_at_change.push_back(Rp->d(jz));
                from_to_handle.insert(std::make_pair(from_to, handles_at_change));
              }
            
            }
         }
      }
    }
    
    std::map< std::pair<int, int>, std::set<node_handle> > from_to_uniqhandle; 
                                                 
   for(auto pit = from_to_handle.begin();pit!=from_to_handle.end();pit++)
   {
   if(pit->second.size()==set_mxd.size())
     { 
       std::set<node_handle> set_of_handles_down_ij(pit->second.begin(), pit->second.end());
       from_to_uniqhandle.insert(std::make_pair(pit->first,set_of_handles_down_ij));
     }
   }                                              
                                                 
    // loop over from-to pairs
                                                 
    for (auto pit  = from_to_uniqhandle.begin(); pit!= from_to_uniqhandle.end(); pit++) {
        
        int i = pit->first.first;
        int j = pit->first.second;
        if(j==-1) continue;
        
        
        
        std::map<int, std::set<node_handle>> seHandlesForRecursion = seHandles;
        //seHandlesForRecursion.insert(std::make_pair(arg2F->getNodeLevel(dvector[jz]), dvector[jz]));
        for(auto pdit = pit->second.begin(); pdit != pit->second.end(); pdit++)
          {
          int level_of_downij = arg2F->getNodeLevel(*pdit);
          
          if( seHandlesForRecursion.find(level_of_downij) != seHandlesForRecursion.end())
            seHandlesForRecursion.at(level_of_downij).insert(*pdit);
          else
            {
            std::set<node_handle> handles_at_ij;
            handles_at_ij.insert(*pdit);
            seHandlesForRecursion.insert(std::make_pair(level_of_downij, handles_at_ij));
            }
          }
        seHandlesForRecursion.erase(mxdLevel); 
        
        node_handle newstates = recFire(A->d(i), seHandlesForRecursion);
        if (0==newstates) continue;
        
        //confirm local state
        if(!rel->isConfirmedState(rLevel,j)) // if not confirmed before
          {
          rel->setConfirmedStates(rLevel,j); // confirm and enlarge
          if (j >= nb->getSize()) {
            int new_var_bound = resF->isExtensibleLevel(nb->getLevel())? -(j+1) : (j+1);
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
        
      } // for i-j pairs
   
    
  } // else
  
  // cleanup mdd reader
  unpacked_node::recycle(A);
  
  printf("\nrecFire Calling saturateHelper on node at level %d", nb->getLevel());
  saturateHelper(*nb);
  printf("\n Called saturateHelper on node at level %d", nb->getLevel());
  result = resF->createReducedNode(-1, nb);
  printf("\n Called saturateHelper & got handle  %d at level %d", result, nb->getLevel());
  
#ifdef TRACE_ALL_OPS
  printf("computed recfire(%d, %d) = %d\n", mdd, mxd, result);
#endif
#ifdef TRACE_RECFIRE
  printf("computed recfire(%d, %d) = %d\n", mdd, mxd, result);
  printf("  node %3d ", result);
  resF->showNode(stdout, result, 1);
  printf("\n");
#endif
  
  //saveResult(Key, mdd, mxd_vector, result);
  //node_handle resultTest = 0;
  //findResult(mdd,mxd,resultTest);
  return saveResult(Key, mdd, mxd_vector,result);
}


















// ****************************************************
//
//                  Old code
//
// ****************************************************

#if 0
void MEDDLY::forwd_impl_dfs_by_events_mt::saturateHelper(unpacked_node& nb)
{
  int nEventsAtThisLevel = rel->lengthForLevel(nb.getLevel());
  printf("\n I am in satHelper for level %d",nb.getLevel());
  if (0 == nEventsAtThisLevel) return;
  
  // Initialize mxd readers, note we might skip the unprimed level
  const int level = nb.getLevel();
  satimpl_opname::event** event_handles = rel->arrayForLevel(level);
  
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
    
    //bool is_union = rel->isUnionPossible(nb.getLevel(),i,Ru);
    //if(!is_union)
      {
      for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
        
        if(event_handles[ei]->getNumOfSubevents()>0) event_handles[ei]->rebuild(); // rebuild only if any change has happened!
        
        // obtain handles of relnodes/subevents of this event 
        // into a vector of length = number of levels
        // initialized for handle 0
        std::vector<int> eVector;
        const int* eVars = event_handles[ei]->getVars();
        int num_eVars = event_handles[ei]->getNumVars();
        
        for(int k=0,l=0;k<arg2F->getNumVariables()+1; k++)
          {
            if(k == 0)
              eVector.push_back(0);
            else if( (l < num_eVars) && (k == eVars[l]) )
              {
                // check if the top handle is just extensible 
                // if so, then put the down handle which is not extensible and is actually useful
                node_handle se_nh = event_handles[ei]->getComponentAt(eVars[l]);
                printf("\n get se_nh = %d from level %d\n",se_nh,eVars[l]);
                 printf("\n &&&&&&&&&&&& Adding minterm to forest %d &&&&&&&&&&&&", arg2F);
                  eVector.push_back(se_nh);
                /*if( arg2F->isImplicit(se_nh) || !arg2F->isExtensible(se_nh) || arg2F->getExtensibleIndex(se_nh)>0 ) eVector.push_back(se_nh);
                else  {
                    printf("\n Need to do this node handle %d is at level %d", se_nh, arg2F->getNodeLevel(se_nh));
                    node_handle notext_se_nh = arg2F->getDownPtr(se_nh, 0);
                    while(arg2F->isExtensible(notext_se_nh))
                        notext_se_nh = arg2F->getDownPtr(notext_se_nh, 0);
                    eVector.push_back(notext_se_nh);
                 }*/
                l++;
              }
            else  
               eVector.push_back(0);
          }
        
       
        int jc = 0;
        std::vector<int> jvector;
        std::vector<node_handle> dvector;
        bool imFlag = arg2F->isImplicit(eVector[level]);
        if(imFlag)
          {
            relation_node* Ru = arg2F->buildImplicitNode(eVector[level]);
            jc = 1;
            printf("\n I get here about to call nextOf???\n");
            jvector.push_back(Ru->nextOf(i)); //user must have nextOf defined as per mapping. so virtual.
            printf("\nCould call nextOf???\n");
            dvector.push_back(Ru->getDown());
          } else {
            unpacked_node* Ru = unpacked_node::useUnpackedNode();
            Ru->initFromNode(arg2F, eVector[level], true);
            node_handle ei_p = Ru->d(i);
            printf("\n I get here new handle from Ru->d(%d) = %d\n",i,ei_p);
            if (0 == ei_p) continue;
            unpacked_node* Rp = unpacked_node::useUnpackedNode();
            Rp->initFromNode(arg2F, ei_p, false);
            jc = Rp->getNNZs();
            for(int jz = 0; jz < jc; jz++)
              {
                jvector.push_back(Rp->i(jz));
                dvector.push_back(Rp->d(jz));
              }
          }
       printf("\n I get here jcount = %d???\n",jc);
       for(int jz = 0; jz < jc; jz++)
       { 
        int j = jvector[jz];
        printf("\n i=%d  j = %d\n",i, j);
        if(j==-1) continue;
        if (j < nb.getSize() && -1==nb.d(j)) continue; // nothing can be added to this set
        
        int event_down_level = event_handles[ei]->downLevel(level);
        printf("\n Current level = %d, event_down_level = %d where handle is %d",level,event_down_level,dvector[jz]);
        int subevent_down_level = arg2F->getNodeLevel(dvector[jz]); // we may get empty extensible nodes
         printf("\n subevent_down_level level = %d",subevent_down_level);
         
        if( (subevent_down_level == event_down_level) && 
            (subevent_down_level == 0 || !arg2F->isImplicit(dvector[jz]) ) &&
            //(arg2F->isExtensible(dvector[jz]) && (arg2F->getExtensibleIndex(dvector[jz]) != 0)) )
            ( (eVector[event_down_level] != 0) || arg2F->isImplicit(eVector[event_down_level]) ) )
          {
            // Should not be replacing an implicit node
            MEDDLY_DCASSERT(!arg2F->isImplicit(eVector[event_down_level]) || (eVector[event_down_level] == 0) );
            eVector[event_down_level] = dvector[jz] ;
          }
        
        // put the handle corresponding to next level at eVector[0]
         eVector[0] = event_down_level;
         
        // eVector[0] = eVector[event_down_level];
        
         printf("\n Going to ");
        node_handle rec = recFire(nb.d(i), eVector, event_handles[ei]);
         printf("\nrec = %d",rec);
        if (rec == 0) continue;
        
        //confirm local state
        if(!rel->isConfirmedState(level,j))
          rel->setConfirmedStates(level, j);
        
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
          
      } //for all j's
        
      } // for all events, ei
      } // No union possible
    /*else { 
    std::unordered_map<long,std::vector<node_handle>> list_of_j = rel->getListOfNexts(nb.getLevel(),i, Ru);
     
    for (std::unordered_map<long,std::vector<node_handle>>::iterator jt=list_of_j.begin(); jt!=list_of_j.end(); ++jt) {
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
    
   } //if union possible*/
  
 }// more indexes to explore
  
  //delete[] Ru;
  recycle(queue);
}

MEDDLY::node_handle MEDDLY::forwd_impl_dfs_by_events_mt::recFireSet(
                                                                    MEDDLY::node_handle mdd, 
                                                                    std::vector<node_handle> vector_mxd)
{
  std::vector<node_handle> array_rec(vector_mxd.size());

  dd_edge ans(resF), union_rec(resF);
  
  /*for(int rn = 0; rn < vector_mxd.size(); rn ++){
    ans.set( recFire(mdd,vector_mxd[rn]) );
    mddUnion->computeTemp(union_rec, ans, union_rec);
  }*/
  
  return union_rec.getNode();
}


// Same as post-image, except we saturate before reducing.
MEDDLY::node_handle MEDDLY::forwd_impl_dfs_by_events_mt::recFire(
                                                                 MEDDLY::node_handle mdd, std::vector<int> eVector,
                                                                 satimpl_opname::event* ei
                                                                 )
{
  
  int mxd = eVector[0]==0?-1:eVector[eVector[0]];
  std::vector<int> mxd_vector;
  //int mxd = eVector[0];
  
  // Special meaning
  // eVector[0] holds the level that current recFire call must operate on for the event
  // terminal conditions can be found at eVector[eVector[0]]
  
  
  
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  
  if (mxd == -1) {
    printf("\n THIS MIGHT GO");
    if (arg1F->isTerminalNode(mdd)) {
      { printf("\n THIS WILL GO");
      return resF->handleForValue(1);
      }
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }
  
  printf("\n mddLevel is %d mxdLevel is %d\n",arg1F->getNodeLevel(mdd),arg2F->getNodeLevel(mxd) );
  
  if(mxd!=-1)
  for(int i = eVector[0]; i < eVector.size(); i++)
    {
    //if(!arg2F->isImplicit(eVector[i]))
      {
      mxd_vector.push_back(eVector[i]);
      }
    }
  else mxd_vector.push_back(-1);
  
  //expert_forest* mixedRelF = rel->getRelForest(); 
  //bool imFlag = mixedRelF->isImplicit(mxd);
  
  
  // check the cache
  FILE_output meddlyout(stdout);
  node_handle result = 0;
  compute_table::entry_key* Key = findResult(mdd, mxd_vector, result);
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
  const int rLevel = MAX(mxdLevel, mddLevel);
  int rSize = resF->getLevelSize(rLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);
  expert_domain* dm = static_cast<expert_domain*>(resF->useDomain());

  dd_edge nbdj(resF), newst(resF);
  
  
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
    printf("\n mddLvel is higher\n");
    for (int i=0; i<rSize; i++) {
      nb->d_ref(i) = recFire(A->d(i), eVector, ei);
    }
    
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(mxdLevel >= mddLevel);
    
    // clear out result (important!)    
    for (unsigned i=0; i<rSize; i++) nb->d_ref(i) = 0;
    
    
    int jc = 0;
    std::vector<int> jvector;
    std::vector<node_handle> dvector;
    bool imFlag = arg2F->isImplicit(mxd);
    int row_size = rSize;
    
    if(imFlag) {
        relation_node* relNode = rel->getRelForest()->buildImplicitNode(mxd);
        row_size = rSize;
      } 
    else {
        unpacked_node *Ru = unpacked_node::useUnpackedNode();
        if (mxdLevel < 0) {
          Ru->initRedundant(arg2F, rLevel, mxd, false);
        } else {
          Ru->initFromNode(arg2F, mxd, false);
        }
        row_size = Ru->getNNZs();
     }
  
    
    for (int iz=0; iz<row_size; iz++) {
        jvector.clear();
        dvector.clear();
        int i = iz;
        // Initialize mxd readers, note we might skip the unprimed level
        if(imFlag) {
          
            relation_node* relNode = rel->getRelForest()->buildImplicitNode(mxd);
            jc = 1;
            jvector.push_back(relNode->nextOf(i));
            dvector.push_back(relNode->getDown());
            } else  { 
            
              unpacked_node *Ru = unpacked_node::useUnpackedNode();
              unpacked_node *Rp = unpacked_node::useUnpackedNode();
              if (mxdLevel < 0) {
                Ru->initRedundant(arg2F, rLevel, mxd, false);
              } else {
                Ru->initFromNode(arg2F, mxd, false);
              }
              
              i = Ru->i(iz);
          
              if (0==A->d(i)) continue; 
              const node_handle pnode = Ru->d(iz);
              if (isLevelAbove(-rLevel, arg2F->getNodeLevel(pnode))) {
                Rp->initIdentity(arg2F, rLevel, i, pnode, false);
              } else {
                Rp->initFromNode(arg2F, pnode, false);
              }
              jc = Rp->getNNZs();
              // loop over mxd "columns"
              for(int jz = 0; jz < jc; jz++)
              {
                jvector.push_back(Rp->i(jz));
                dvector.push_back(Rp->d(jz));
              }
            
          }
       
       // loop over mxd "columns"
       for (int jz=0; jz<jc; jz++) {
         
          int j = jvector[jz];
          if(j==-1) continue;
          
          
          int subevent_down_level = arg2F->getNodeLevel(dvector[jz]);
          if( (subevent_down_level!=0) && ( (eVector[subevent_down_level] != 0) || arg2F->isImplicit(eVector[subevent_down_level])) )
            {
              eVector[subevent_down_level] = dvector[jz] ;
              //if(rLevel == subevent_down_level + 1) eVector[0] = dvector[jz];
            }
         
           int dec = mxdLevel - 1;
           // Go down eVector and find the first non-zero entry.
           while( (dec>=1) && (!eVector[dec]) ) {
             dec--;}
           if(dec!=0) 
              eVector[0] = dec;
           else 
             eVector[0] = 0;
          
         
          node_handle newstates = recFire(A->d(i), eVector, ei);
          if (0==newstates) continue;
          
          //confirm local state
          if(!rel->isConfirmedState(rLevel,j)) // if not confirmed before
            {
            rel->setConfirmedStates(rLevel,j); // confirm and enlarge
            if (j >= nb->getSize()) {
              int new_var_bound = resF->isExtensibleLevel(nb->getLevel())? -(j+1) : (j+1);
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
          
        } // for j

      } // for i
    
  } // else
  
  // cleanup mdd reader
  unpacked_node::recycle(A);
  
  
  saturateHelper(*nb);
  result = resF->createReducedNode(-1, nb);
  printf("\n The node is created with handle %d", result);
  
  #ifdef TRACE_ALL_OPS
  printf("computed recfire(%d, %d) = %d\n", mdd, mxd, result);
  #endif
#ifdef TRACE_RECFIRE
  printf("computed recfire(%d, %d) = %d\n", mdd, mxd, result);
  printf("  node %3d ", result);
  resF->showNode(stdout, result, 1);
  printf("\n");
#endif
  
  //saveResult(Key, mdd, mxd_vector, result);
  //node_handle resultTest = 0;
  //findResult(mdd,mxd,resultTest);
  return saveResult(Key, mdd, mxd_vector, result);
}

#endif

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
  arg2F = static_cast<expert_forest*>(rel->getMixRelForest());
  resF = static_cast<expert_forest*>(rel->getOutForest());
  
  registerInForest(arg1F);
  registerInForest(arg2F);
  registerInForest(resF);
  compute_table::entry_type* et = new compute_table::entry_type(opcode->getName(), "N.N:N");
  et->setForestForSlot(0, arg1F);
  et->setForestForSlot(2, arg2F);
  et->setForestForSlot(4, resF);
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

/*bool MEDDLY::common_impl_dfs_by_events_mt
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
}*/

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
  compute_table::entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new compute_table::entry_type(name, "NI:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(3, resF);
  } else {
    et = new compute_table::entry_type(name, "N:N");
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
  compute_table::entry_key* Key = findSaturateResult(mdd, k, n);
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
  printf("\n Calling saturateHelper on node at level %d", nb->getLevel());
  parent->saturateHelper(*nb);
  printf("\n Abt Completion Called saturateHelper on node at level %d", nb->getLevel());
  
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

/*struct node_pair {
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
}*/

// Modified version of saturate() for detecting reachability.
// Recursively calls isReachable (similar to saturate() calling saturate()),
// and finally saturates the new node before returning result of saturation.
// Note that if isReachable returns false, then the information returned via
// saturation_result is the saturated node,
// and if isReachable returns true, then there are no guarantees for the information
// returned via saturated_result.
/*bool
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
  compute_table::entry_key* Key = findSaturateResult(mdd, k, n);
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
     node_handle mxd,
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
  compute_table::entry_key* Key = findResult(mdd, mxd, result);
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
*/

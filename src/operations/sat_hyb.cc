
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

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../oper_special.h"
#include "../opname_satur.h"
#include "../ops_builtin.h"

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
// *                     sathyb_opname  methods                    *
// *                                                                *
// ******************************************************************


MEDDLY::sathyb_opname::sathyb_opname(const char* n)
: specialized_opname(n)
{
}

MEDDLY::sathyb_opname::~sathyb_opname()
{
}

/*long MEDDLY::relation_node::delta()
{
  //to be defined for the example you use & comment this definition
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}*/

// ******************************************************************

MEDDLY::sathyb_opname::hybrid_relation::hybrid_relation(forest* inmdd, forest* relmxd,
                                                             forest* outmdd, event** E, int ne)
: insetF(static_cast<expert_forest*>(inmdd)), outsetF(static_cast<expert_forest*>(outmdd)), hybRelF(static_cast<expert_forest*>(relmxd))
{
  hybRelF = static_cast<expert_forest*>(relmxd);

  if (0==insetF || 0==outsetF || 0==hybRelF ) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);

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
      !hybRelF->isForRelations()  ||
      (insetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)   ||
      (outsetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)  ||
      (hybRelF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
      )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  // Forests are good; set number of variables
  num_levels = insetF->getDomain()->getNumVariables() + 1;

  //Allocate event_list
  event_list = (node_handle**)malloc(unsigned(num_levels)*sizeof(node_handle*));
  event_list_alloc = (long*)malloc(unsigned(num_levels)*sizeof(long));
  event_added = (long*)malloc(unsigned(num_levels)*sizeof(long));

  //Confirmed local states
  confirm_states = (int*)malloc(unsigned(num_levels)*sizeof(int));
  confirmed_array_size = (int*)malloc(unsigned(num_levels)*sizeof(int));
  confirmed = new bool*[num_levels];
  confirmed[0]=0;
  for(int i = 1;i<num_levels;i++)
    {
    event_list[i] = (node_handle*)malloc(8*sizeof(node_handle));
    int plevel_size = hybRelF->getLevelSize(-i);
    int unplevel_size = hybRelF->getLevelSize(i);
    int level_size = plevel_size>unplevel_size ? plevel_size:unplevel_size;
    confirmed[i] = (bool*)malloc(unsigned(level_size)*sizeof(bool));
    event_list_alloc[i] = 8;
    event_added[i] = 0;
    confirm_states[i] = 0;

    confirmed_array_size[i]=level_size;
    for(int j = 0;j<level_size;j++)
      confirmed[i][j]=false;
    }
  // Build the events-per-level data structure
  // (0) Initialize
  num_events_by_top_level = new int[num_levels]; // number of events that have this level as its top-level participant
  num_events_by_level = new int[num_levels]; // number of events that have this level as its participant
  memset(num_events_by_top_level, 0, sizeof(int)*unsigned(num_levels));
  memset(num_events_by_level, 0, sizeof(int)*unsigned(num_levels));
  // (1) Determine the number of events per level
  for (int i = 0; i < ne; i++) {
    //if (E[i]->isDisabled()) continue;
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
    //if (E[i]->isDisabled()) continue;
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
    //if (E[i]->isDisabled()) continue;
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
    //if (E[i]->isDisabled()) continue;
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
    //if (E[i]->isDisabled()) continue;
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
   //if (E[i]->isDisabled()) continue;
    int nrn = E[i]->getNumOfRelnodes();
    relation_node** rn = E[i]->getRelNodes();
    for (int j = 0; j < nrn; j++) {
      int level = rn[j]->getLevel();
      relnodes_by_level[level][num_relnodes_by_level[level]++] = rn[j];
    }
  }
}


void
MEDDLY::sathyb_opname::hybrid_relation::resizeEventArray(int level)
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
MEDDLY::sathyb_opname::hybrid_relation::resizeConfirmedArray(int level, int index)
{
  #if 0
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
        {
          printf("\n Setting level %d index %d as false", level, i);
          confirmed[level][i]=false;
        }

         confirmed_array_size[level]=nalloc;
    }
    #else
    if (index <= confirmed_array_size[level]) return;
    index = MAX( index, confirmed_array_size[level]*2);
    confirmed[level] = (bool*) realloc(confirmed[level], unsigned(index) * sizeof(bool));
    if (confirmed[level] == 0) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    for (int i = confirmed_array_size[level]; i < index; i++) confirmed[level][i] = false;
    confirmed_array_size[level] = index;
  #endif

}

void findConfirmedStatesImpl(MEDDLY::sathyb_opname::hybrid_relation* rel,
                             bool** confirmed, int* confirm_states,
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

void MEDDLY::sathyb_opname::hybrid_relation::setConfirmedStates(const dd_edge& set)
{
  // Perform a depth-first traversal of set:
  //    At each level, mark all enabled states as confirmed.

  // Enlarge the confirmed arrays if needed
  for (int i = 1 ; i<num_levels; i++)
    {
      int levelSize = hybRelF->getLevelSize(-i);
      resizeConfirmedArray(i, levelSize);
    }

    std::set<node_handle> visited;
    findConfirmedStatesImpl(const_cast<hybrid_relation*>(this),
                      confirmed, confirm_states, set.getNode(), num_levels-1, visited);

}

void
MEDDLY::sathyb_opname::hybrid_relation::setConfirmedStates(int level, int index)
{
  // For each subevent that affects this level:
  //    (1) call subevent::confirm()
  //    (2) for each level k affected by the subevent,
  //        (a) enlarge variable bound of k to variable bound of k'


  resizeConfirmedArray(level, index+1);

  MEDDLY_DCASSERT(confirmed_array_size[level] > index);
  if (isConfirmedState(level, index)) return;


  // Get subevents affected by this level, and rebuild them.
  int nSubevents = num_subevents_by_level[level];

  if (nSubevents>0) {

    for (int i = 0; i < nSubevents; i++) {
      subevents_by_level[level][i]->confirm(const_cast<hybrid_relation&>(*this),
                                            level, index);
     }


    const int nEvents = num_events_by_level[level];
    for (int i = 0; i < nEvents; i++) {
      if(events_by_level[level][i]->getNumOfSubevents()>0)
        events_by_level[level][i]->markForRebuilding();

    }

  }

  confirmed[level][index] = true;
  confirm_states[level]++;
  //return true;
}



MEDDLY::sathyb_opname::hybrid_relation::~hybrid_relation()
{
  /*last_in_node_array = 0;
  impl_unique.clear();

  for(int i = 0; i <=num_levels; i++) {delete[] event_list[i]; delete[] confirmed[i];}
  delete[] event_list;
  delete[] event_added;
  delete[] event_list_alloc;
  delete[] confirmed;
  delete[] confirm_states;
  delete[] confirmed_array_size;*/
}

void
MEDDLY::sathyb_opname::hybrid_relation::show()
{
  /*node_handle** event_list_copy = (node_handle**)malloc((num_levels)*sizeof(node_handle*));
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
  delete[] event_list_copy;*/

}

void MEDDLY::sathyb_opname::hybrid_relation::bindExtensibleVariables() {
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

#if 0
MEDDLY::node_handle
MEDDLY::sathyb_opname::hybrid_relation::buildMxdForest()
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

  forest* mxd = d->createForest(true,forest::BOOLEAN, edge_labeling::MULTI_TERMINAL);
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
MEDDLY::sathyb_opname::hybrid_relation::buildEventMxd(node_handle eventTop, forest *mxd)
{
  //mxd is built on a domain obtained from result of saturation
  int nVars = outsetF->getDomain()->getNumVariables();
  //int* sizes = new int[nVars];
  relation_node* Rnode = nodeExists(eventTop);
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
            Rnode = nodeExists(rnh_array[i]);
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
}
#endif

// ******************************************************************


std::unordered_map<long,std::vector<MEDDLY::node_handle>>
MEDDLY::sathyb_opname::hybrid_relation::getListOfNexts(int level, long i, relation_node **R)
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
MEDDLY::sathyb_opname::hybrid_relation::isUnionPossible(int level, long i, relation_node **R)
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

MEDDLY::sathyb_opname::subevent::subevent(forest* f, int* v, int nv, bool firing)
: vars(0), num_vars(nv), root(dd_edge(f)), top(0),
f(static_cast<expert_forest*>(f)), is_firing(firing)
{
  MEDDLY_DCASSERT(f != 0);
  MEDDLY_DCASSERT(v != 0);
  MEDDLY_DCASSERT(nv > 0);

  vars = new int[num_vars];
  memcpy(vars, v, unsigned(num_vars) * sizeof(int));

  // find top
  top = vars[0];
  for (int i = 1; i < num_vars; i++) {
    if (isLevelAbove(vars[i], top)) top = vars[i];
  }

  uses_extensible_variables = false;
  for (int i = 0; i < num_vars; i++) {
    if (this->f->isExtensibleLevel(vars[i])) {
      uses_extensible_variables = true;
      break;
    }
  }

  down = -1;
  unpminterms = pminterms = 0;
  num_minterms = size_minterms = 0;
  process_minterm_pos = -1;
  processed_minterm_pos = -1;
}

MEDDLY::sathyb_opname::subevent::~subevent()
{
  if (vars) delete [] vars;
  for (int i=0; i<num_minterms; i++) {
    delete[] unpminterms[i];
    delete[] pminterms[i];
  }
  free(unpminterms);
  free(pminterms);
}

void MEDDLY::sathyb_opname::subevent::clearMinterms()
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


void MEDDLY::sathyb_opname::subevent::confirm(hybrid_relation& rel, int v, int i) {
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}


bool MEDDLY::sathyb_opname::subevent::addMinterm(const int* from, const int* to)
{


   ostream_output out(std::cout);
   /*out << "Adding MEDDLY minterm: [";
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
 //  out << "]\n";
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

void MEDDLY::sathyb_opname::subevent::buildRoot() {
//  printf("\n num_minterms in this se = %d, to be done = %d \n",num_minterms,process_minterm_pos-processed_minterm_pos );
  if (0 == num_minterms) return;
  if (1 == num_vars) return ;


   ostream_output out(std::cout);
   /*out << "\nBuilding subevent from " << num_minterms << " minterms\n";
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

// Arrange minterms before union-ing

 #if 0
 bool semi_union = false;

 // adding first minterm
 if(root_handle == 0)
 {
      dd_edge sum(root);
      f->createEdge(unpminterms, pminterms, 1, sum);
      num_minterms -= 1;
      root += sum;
      root_handle = root.getNode();
      semi_union = true;
 } else {
  //already exist
  //lets do CT-less union
  while(num_minterms != 0)
    {
      for(int i=1;i<=f->getNumVariables();i++)
        {
        //unpterms.push_back(unpminterms[processed_minterm_pos+1][i]);
        //pterms.push_back(pminterms[processed_minterm_pos+1][i]);

        unpterms.push_back(unpminterms[num_minterms-1][i]);
        pterms.push_back(pminterms[num_minterms-1][i]);

        }
      terms.push_back(unpterms);
      terms.push_back(pterms);

      rnh = f->unionOneMinterm(root_handle, terms);
      out << "\nEquivalent event: " << rnh << "\n";
      root.set(rnh);
      root_handle = root.getNode();
      //processed_minterm_pos +=1;
      //process_minterm_pos = -1;
      num_minterms --;
      unpterms.clear();
      pterms.clear();
      terms.clear();
    }
 }

 if(semi_union && num_minterms>0)
 {
   while(num_minterms != 0)
    {
      for(int i=1;i<=f->getNumVariables();i++)
        {
        //unpterms.push_back(unpminterms[processed_minterm_pos+1][i]);
        //pterms.push_back(pminterms[processed_minterm_pos+1][i]);

        unpterms.push_back(unpminterms[num_minterms][i]);
        pterms.push_back(pminterms[num_minterms][i]);

        }
      terms.push_back(unpterms);
      terms.push_back(pterms);

      rnh = f->unionOneMinterm(root_handle, terms);
      out << "\nEquivalent event: " << rnh << "\n";
      root.set(rnh);
      root_handle = root.getNode();
      //processed_minterm_pos +=1;
      //process_minterm_pos = -1;
      num_minterms --;
      unpterms.clear();
      pterms.clear();
      terms.clear();
    }
 } else {
    if (usesExtensibleVariables()) {
      dd_edge sum(root);
      f->createEdge(unpminterms, pminterms, num_minterms, sum);
      num_minterms = 0;
      root += sum;
    } else {
      f->createEdge(unpminterms, pminterms, num_minterms, root);

      //dd_edge sum(root);
      //f->createEdge(unpminterms, pminterms, num_minterms, sum);

      //num_minterms = 0;
      //root += sum;

    }
    //processed_minterm_pos = process_minterm_pos;

  }
  #endif


  // Older version: Create mxd, union, destroy mxd
  #if 1
  if (usesExtensibleVariables()) {
      dd_edge sum(root);
      f->createEdge(unpminterms, pminterms, num_minterms, sum);
      num_minterms = 0;
      root += sum;
    } else {
       f->createEdge(unpminterms, pminterms, num_minterms, root);
    }
  #endif


  // Union minterm one-by-one w/o bulding mxd
  #if 0
  for( int w = 0;w <num_minterms; w++){
    int* pminterms1 = new int[f->getNumVariables() + 1];
    int* unpminterms1 = new int[f->getNumVariables() + 1];

    for( int kk=0; kk<=f->getNumVariables(); kk++ )
    {
      pminterms1[kk] = pminterms[w][kk];
      unpminterms1[kk] = unpminterms[w][kk];
    }

    rnh = f->unionOneMinterm(root_handle, unpminterms1, pminterms1, f->getNumVariables());
    root.set(rnh);
    //if(w == num_minterms-1) root.show(out, 2);
    root_handle = root.getNode();
  }
  num_minterms = 0;
  #endif

  root_handle = root.getNode();
  //out << "\nEquivalent event: " << root.getNode() << "\n";
  //out << "Result: ";
  //root.show(out, 2);
}


void MEDDLY::sathyb_opname::subevent::showInfo(output& out) const {
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

long MEDDLY::sathyb_opname::subevent::mintermMemoryUsage() const {
  long n_minterms = 0L;
  for (int i = 0; i < size_minterms; i++) {
    if (unpminterms[i] != 0) n_minterms++;
  }
  return long(n_minterms * 2) * long(f->getDomain()->getNumVariables()) * long(sizeof(int));
}

// ============================================================


MEDDLY::sathyb_opname::event::event(subevent** p, int np, relation_node** r, int nr)
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
      int pos = 0;
      for (int i = 0; i < np; i++)
          {
            for (int k = 0; k < np; k++){
              if ((p[k]->getDown()!=0) && (sorted_by_top_se[i] == p[k]->getTop()))
              {
                subevents[pos] = p[k];
                p[k]->setDown(0);
                pos++;
                break;
              }
            }
          }
      top = subevents[np-1]->getTop(); //sorted_by_top_se[np-1];
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

  // Get all the variables that are part of sub-events
  // Assumption: A variable which is in a sub-event cannot be part of relation_node;
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
  // set number of variables of this event
  num_vars = sVars.size() + rVars.size();

  // set the variables of this event
  vars = new int[num_vars];
  int* curr = &vars[0];

  if((np>0) && (nr>0)) {
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
  } else if(np>0) {
    std::set<int>::iterator sit=sVars.begin();
    while(sit!=sVars.end())
        *curr++ = *sit++;
  } else {
    std::set<int>::iterator rit=rVars.begin();
    while(rit!=rVars.end())
        *curr++ = *rit++;
  }

   // set number of variables that undergo firing of this event
  num_firing_vars = firingVars.size();

  // set the variables that undergo firing of this event
  // all relation nodes participate in firing; Nodetypes = (inhibitor+fire) / (enable+fire) / (fire)
  //
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

  // Create the implicit nodes as a chain
  // Make sure the nodes are sorted by level bottom to top
  node_handle nh_next = -1;
  for (int i = 0; i < nr; i++)
    {
      relnodes[i]->setDown(nh_next); // assign down handle of this node
      nh_next = f->createRelationNode(relnodes[i]); // register the node in the forest & get node_handle
      relnodes[i]->setID(nh_next); // assign node handle to this node
    }


  first_time_build = true;
  all_components = (node_handle*)malloc(sizeof(node_handle)*(num_subevents+(num_relnodes>0?1:0)));
  rebuild();
  root = dd_edge(f);
  event_mask = dd_edge(f);
  event_mask_from_minterm = 0;
  event_mask_to_minterm = 0;
  needs_rebuilding = is_disabled ? false: true;

  printf("\n Event cumulative = %d",f->getImplicitTableCount());

}

MEDDLY::sathyb_opname::event::~event()
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


void MEDDLY::sathyb_opname::event::buildEventMask()
{
  MEDDLY_DCASSERT(num_subevents > 0);
  MEDDLY_DCASSERT(f);

  if (0 == event_mask_from_minterm) {
    const size_t minterm_size = size_t(f->getNumVariables()+1);
    event_mask_from_minterm = new int[minterm_size];
    event_mask_to_minterm = new int[minterm_size];

    // all non-affected variables get DONT_CARE
    // MXD is stored with FullyReduced rule
    // Allows on-th-fly intersection w/ relation_nodes since non-affected variable levels will be skipped
    for (unsigned i = 0; i < minterm_size; i++) {
      event_mask_from_minterm[i] = MEDDLY::DONT_CARE;
      event_mask_to_minterm[i] = MEDDLY::DONT_CARE;
    }

    for (int i = 0; i < num_firing_vars; i++) {
      event_mask_to_minterm[firing_vars[i]] = MEDDLY::DONT_CARE;//MEDDLY::DONT_CHANGE;
    }
  }

   f->createEdge(&event_mask_from_minterm, &event_mask_to_minterm, 1, event_mask);

//#ifdef DEBUG_EVENT_MASK
  ostream_output out(std::cout);
  event_mask.show(out, 2);
//#endif
}


bool MEDDLY::sathyb_opname::event::rebuild()
{
  /*printf("\n Before rebuilding, subevent handles :");
  for(auto m = all_components.begin(); m!=all_components.end();m++)
    printf("->%d", *m);*/

    //MEDDLY_DCASSERT(num_subevents > 0);
  if  (!first_time_build && !num_subevents) return false;
  /*if (is_disabled) return false;*/
  if((first_time_build)&&(num_relnodes>0))
    {
      int level_of_relNode = relnodes[num_relnodes-1]->getLevel(); //get level of highest relation node

      std::set<node_handle> handles_at_this_top = level_component[level_of_relNode];
      handles_at_this_top.insert(relnodes[num_relnodes-1]->getID());
      component_se_type.insert(std::make_pair(relnodes[num_relnodes-1]->getID(),true));
      level_component[level_of_relNode] = handles_at_this_top;
      all_components[0] = relnodes[num_relnodes-1]->getID();
    }

  if(num_subevents > 0)
    {
      //needs_rebuilding is set from within saturation
      if (!needs_rebuilding) return false;
      needs_rebuilding = false;
    }


    // Should not store both subevent handle & event handle in level_component map else how to distinguish!
    #ifdef CONJUNCT_SUBEVENTS
      // conjunct all sub-relations
      int max_of_se_level = 0;
      for (int i = 0; i < num_subevents; i++)
      {
        int level_of_subevent = abs(subevents[i]->getTop());
        if(level_of_subevent>max_of_se_level) max_of_se_level = level_of_subevent;
      }

      // Remove old event handle
      level_component[max_of_se_level].erase(partial_root.getNode());
      component_se_type.erase(partial_root.getNode());
      partial_root = dd_edge(event_mask);
      for (int i = 0; i < num_subevents; i++)
      {
        subevents[i]->buildRoot();
        partial_root *= subevents[i]->getRoot();
      }
      // Modify it to only hold the topmost nodehandle!
      level_component[max_of_se_level].insert(partial_root.getNode());
      component_se_type.insert(std::make_pair(partial_root.getNode(),true));

      #if 0
      ostream_output out(std::cout);
      f->useDomain()->showInfo(out);
      for (int i = 0; i < num_subevents; i++) {
        out << "subevent: " << subevents[i]->getRoot().getNode() << "\n";
        subevents[i]->getRoot().show(out, 2);
      }
      out << "event: " << partial_root.getNode() << "\n";
      partial_root.show(out, 2);
      #endif

    #else
      int idx_offset = 0;
      if(num_relnodes>0) idx_offset = 1;
       // Retain the sub-events separately
      for (int i = 0; i < num_subevents; i++)
      {
        int level_of_subevent = abs(subevents[i]->getTop());

        // erase existing handle of subevent
        // This is possible because one level is only part of one sub-event
        if(level_component.find(level_of_subevent)!=level_component.end())
          {
            level_component.at(level_of_subevent).erase(subevents[i]->getRootHandle());
            component_se_type.erase(subevents[i]->getRootHandle());
            all_components[i+idx_offset] = 0;
          }

         subevents[i]->buildRoot();

        // insert new rebuilt handle of subevent
        std::set<node_handle> handles_at_this_top = level_component[level_of_subevent];
        handles_at_this_top.insert(subevents[i]->getRootHandle());
        component_se_type.insert(std::make_pair(subevents[i]->getRootHandle(),subevents[i]->isFiring()));
        level_component[level_of_subevent] = handles_at_this_top;
        all_components[i+idx_offset] = subevents[i]->getRootHandle();
        }

        #if 0
        ostream_output out(std::cout);
        for (int i = 0; i < num_subevents; i++) {
        out << "subevent: " << subevents[i]->getRoot().getNode() << "\n";
        if(subevents[i]->getRoot().getNode()) subevents[i]->getRoot().show(out, 2);
        }
        #endif
    #endif

  // event stores a list of root_handle. Why? because of existense of multiple subevents at same top-level (for general case)
 /* printf("\n After rebuilding:, subevent handles :");
  for(auto m = all_components.begin(); m!=all_components.end();m++)
    printf("->%d", *m);*/

  if (level_component[top] == root_handle) return false;
  root_handle = level_component[top];
  first_time_build = false;
  return true;
}

int MEDDLY::sathyb_opname::event::downLevel(int level) const{
  for(int i = num_vars; i>0; i--)
    {
      if(vars[i] == level)
        return vars[i-1];
    }
  return 0;
}

void MEDDLY::sathyb_opname::event::enlargeVariables()
{
  expert_domain* ed = static_cast<expert_forest*>(f)->useExpertDomain();
  for (int i = 0; i < num_vars; i++) {
    int unprimed = ABS(vars[i]);
    int primed = -unprimed;
    int unprimedSize = f->getLevelSize(unprimed);
    int primedSize = f->getLevelSize(primed);
    if (unprimedSize < primedSize) {
      variable* vh = ed->getExpertVar(unprimed);
      if (vh->isExtensible())
        vh->enlargeBound(false, -primedSize);
      else
        vh->enlargeBound(false, primedSize);
    }
    MEDDLY_DCASSERT(f->getLevelSize(unprimed) == f->getLevelSize(primed));
  }
}


void MEDDLY::sathyb_opname::event::showInfo(output& out) const {
  for (int i = 0; i < num_subevents; i++) {
    out << "subevent " << i << "\n";
    subevents[i]->showInfo(out);
  }
}

long MEDDLY::sathyb_opname::event::mintermMemoryUsage() const {
  long usage = 0;
  for (int i = 0; i < num_subevents; i++) {
    usage += subevents[i]->mintermMemoryUsage();
  }
  return usage;
}


// ******************************************************************
// *                                                                *
// *               saturation_hyb_by_events_opname  class               *
// *                                                                *
// ******************************************************************

/** Simple class to keep compute table happy.
 */
class MEDDLY::saturation_hyb_by_events_opname : public unary_opname {
  static saturation_hyb_by_events_opname* instance;
public:
  saturation_hyb_by_events_opname();

  static saturation_hyb_by_events_opname* getInstance();

};

MEDDLY::saturation_hyb_by_events_opname* MEDDLY::saturation_hyb_by_events_opname::instance = 0;

MEDDLY::saturation_hyb_by_events_opname::saturation_hyb_by_events_opname()
: unary_opname("Saturate_by_events")
{
}

MEDDLY::saturation_hyb_by_events_opname* MEDDLY::saturation_hyb_by_events_opname::getInstance()
{
  if (0==instance) instance = new saturation_hyb_by_events_opname;
  return instance;
}

// ******************************************************************
// *                                                                *
// *             saturation_hyb_by_events_op  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_hyb_by_events_op : public unary_operation {
  common_hyb_dfs_by_events_mt* parent;
public:
  saturation_hyb_by_events_op(common_hyb_dfs_by_events_mt* p,
                               expert_forest* argF, expert_forest* resF);
  virtual ~saturation_hyb_by_events_op();

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

class MEDDLY::common_hyb_dfs_by_events_mt : public specialized_operation {
public:
  common_hyb_dfs_by_events_mt(sathyb_opname* opcode,
                               sathyb_opname::hybrid_relation* rel);
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

  sathyb_opname::hybrid_relation* rel;

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
// *               forwd_hyb_dfs_by_events_mt class                *
// *                                                                *
// ******************************************************************


class MEDDLY::forwd_hyb_dfs_by_events_mt : public common_hyb_dfs_by_events_mt {
public:
  forwd_hyb_dfs_by_events_mt(sathyb_opname* opcode,
                              sathyb_opname::hybrid_relation* rel);
protected:
  virtual void saturateHelper(unpacked_node& nb);
 MEDDLY::node_handle recFire(MEDDLY::node_handle mdd, MEDDLY::node_handle* seHandles, int num_se);
  void recFireHelper(const unsigned, const int, const MEDDLY::node_handle, const MEDDLY::node_handle,
        unpacked_node*, unpacked_node*, MEDDLY::node_handle*, bool, int, int, int);

  //MEDDLY::node_handle recFireSet(node_handle mdd, std::set<node_handle> mxd);
  std::pair<MEDDLY::node_handle,int> getHighestNodeHandles(MEDDLY::node_handle* seHandles, int num_se);


 };


MEDDLY::forwd_hyb_dfs_by_events_mt::forwd_hyb_dfs_by_events_mt(
                                                                 sathyb_opname* opcode,
                                                                 sathyb_opname::hybrid_relation* rel)
: common_hyb_dfs_by_events_mt(opcode, rel)
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
	    int p_lvl = arg2F->getNodeLevel(seHandles[it]);
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
  sathyb_opname::event** event_handles = rel->arrayForLevel(level);

  dd_edge nbdj(resF), newst(resF);

  expert_domain* dm = static_cast<expert_domain*>(resF->useDomain());

  int* event_Ru_Rn_index = (int*)malloc(nEventsAtThisLevel*sizeof(int));
  unpacked_node** Ru = new unpacked_node*[nEventsAtThisLevel];
  relation_node** Rn = new relation_node*[nEventsAtThisLevel];
  unpacked_node* Rp = unpacked_node::New();
  int i_Ru = 0;
  int i_Rn = 0;

  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    event_handles[ei]->rebuild();
    //std::set<node_handle> seHandles = event_handles[ei]->getComponentAt(level);
    node_handle se_nh = event_handles[ei]->getTopComponent();
    if(arg2F->isImplicit(se_nh)) {
      Rn[i_Rn] = arg2F->buildImplicitNode(se_nh);
      event_Ru_Rn_index[ei] = i_Rn;
      i_Rn++;
    } else {
      Ru[i_Ru] = unpacked_node::New();
      if(arg2F->getNodeLevel(se_nh) == level)
          arg2F->unpackNode(Ru[i_Ru], se_nh, FULL_ONLY);  // node is present at unprime-level
      else
         Ru[i_Ru]->initRedundant(arg2F, level, se_nh, false);
      event_Ru_Rn_index[ei] = i_Ru;
      i_Ru++;
    }
  }

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


      bool imFlag = arg2F->isImplicit(se_nh);
      int jC;
      if(imFlag)
      {
        int j = Rn[event_Ru_Rn_index[ei]]->nextOf(i);
        jC = 1;

      } else {
        for(int w = 0; w<event_handles[ei]->getNumOfSubevents(); w++) {
        sathyb_opname::subevent* my_se = event_handles[ei]->getSubevents()[w];
        if(se_nh == my_se->getRootHandle())
          {
           which_se = w; break;
          }
        }

        if(is_rebuilt)
        {
        if(arg2F->getNodeLevel(se_nh) == level)
          arg2F->unpackNode(Ru[event_Ru_Rn_index[ei]], se_nh, FULL_ONLY);  // node is present at unprime-level
        else
          Ru[event_Ru_Rn_index[ei]]->initRedundant(arg2F, level, se_nh, false);  // node was at prime-level, so build redudant node at unprime-level
        }

        node_handle ei_i_p = Ru[event_Ru_Rn_index[ei]]->d(i);
        if( 0 == ei_i_p) continue;

        const int dlevel = arg2F->getNodeLevel(ei_i_p);

        if (dlevel == -level)
          arg2F->unpackNode(Rp, ei_i_p, SPARSE_ONLY);
        else
          Rp->initIdentity(arg2F, -level, i, ei_i_p, false);

        jC = Rp->getNNZs();

      }
      int num_se = event_handles[ei]->getNumOfSubevents()+(event_handles[ei]->getNumOfRelnodes()>0?1:0);
      node_handle* seHandlesForRecursion = (node_handle*)malloc(num_se*sizeof(node_handle));
      memcpy(seHandlesForRecursion,event_handles[ei]->getAllComponents(),num_se*sizeof(node_handle));


      for(int jz = 0; jz < jC; jz ++) {
        int j;
        if(imFlag) {
         j =  Rn[event_Ru_Rn_index[ei]]->nextOf(i);
         if(j==-1) continue; // Not enabled
         if (j < nb.getSize() && -1==nb.d(j)) continue;
         seHandlesForRecursion[0]=Rn[event_Ru_Rn_index[ei]]->getDown();
        } else {
          j = Rp->i(jz);
          if(j==-1) continue; // Not enabled
          if (j < nb.getSize() && -1==nb.d(j)) continue;
          seHandlesForRecursion[which_se+does_rn_exist] = Rp->d(jz);
        }

        node_handle rec = recFire(nb.d(i), seHandlesForRecursion, num_se);

        if (rec == 0) continue;

        //confirm local state
        if(!rel->isConfirmedState(level,j))
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
    if (arg1F->isTerminalNode(mdd)) {
      return resF->handleForValue(1);
    }
    // mxd is identity
    if (arg1F == resF)
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
  arg1F->showNode(stdout, mdd, 1);
  printf("\n  node %3d ", mxd);
  arg2F->showNode(stdout, mxd, 1);
  printf("\n");
  #endif

  // check if mxd and mdd are at the same level
  const int mddLevel = arg1F->getNodeLevel(mdd);
  const int mxdLevel = arg2F->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
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
  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.
    for (int i=0; i<rSize; i++) {
      nb->d_ref(i) = recFire(A->d(i), seHandles, num_se);
    }

  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(mxdLevel >= mddLevel);

    // Initialize mxd readers, note we might skip the unprimed level

    bool imFlag = arg2F->isImplicit(mxd);
    int row_size = rSize;
    relation_node* relNode;
    unpacked_node *Ru = unpacked_node::New();
    unpacked_node *Rp = unpacked_node::New();

    if(imFlag) {
      relNode = arg2F->buildImplicitNode(mxd);
      row_size = rSize; // Get the update array size;
    }
    else {
       if (mxdLevel < 0) {
        Ru->initRedundant(arg2F, rLevel, mxd, false);
      } else {
        arg2F->unpackNode(Ru, mxd, SPARSE_ONLY);
      }
      row_size = Ru->getNNZs();
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
        recFireHelper(i, rLevel, relNode->getDown(), A->d(i), Rp, nb, seHandlesLower, imFlag, j, num_se, which_se);
       }
       } else  {
        for (int iz=0; iz<row_size; iz++) {
        const unsigned i = Ru->i(iz);
        if (0==A->d(i)) continue;
        recFireHelper(i, rLevel, Ru->d(iz), A->d(i), Rp, nb, seHandlesLower, imFlag, -1, num_se, which_se);
        }
        unpacked_node::recycle(Rp);
        unpacked_node::recycle(Ru);
    }

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
  if(-rLevel == arg2F->getNodeLevel(Ru_i))
     arg2F->unpackNode(Rp, Ru_i, SPARSE_ONLY);
  else
     return;
  }

  dd_edge nbdj(resF), newst(resF);

  unsigned jc = imFlag ? 1 : Rp->getNNZs();

  for (unsigned jz=0; jz<jc; jz++) {
    if((jz == 0) && imFlag)
      {
        if (j == -1) continue;
        recHandles[which_se] = Ru_i;
      }
    else {
        j = Rp->i(jz);
        if (j == -1) continue;
        recHandles[which_se] = Rp->d(jz);
      }

      MEDDLY::node_handle newstates = recFire(A_i, recHandles, num_se);
      if (0==newstates) continue;

      //confirm local state
      if(!rel->isConfirmedState(rLevel,j)) // if not confirmed before
      {
        rel->setConfirmedStates(rLevel,j); // confirm and enlarge
        if (j >= nb->getSize()) {
          int new_var_bound = resF->isExtensibleLevel(nb->getLevel())? -(j+1) : (j+1);
          expert_domain* dm = static_cast<expert_domain*>(resF->useDomain());
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

  }


}


// ******************************************************************
// *                                                                *
// *             common_hyb_dfs_by_events_mt  methods              *
// *                                                                *
// ******************************************************************

MEDDLY::common_hyb_dfs_by_events_mt::common_hyb_dfs_by_events_mt(
                                                                   sathyb_opname* opcode,
                                                                   sathyb_opname::hybrid_relation* relation)
: specialized_operation(opcode, 1)
{
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;
  freeqs = 0;
  freebufs = 0;
  rel = relation;
  arg1F = static_cast<expert_forest*>(rel->getInForest());
  arg2F = static_cast<expert_forest*>(rel->getHybridForest());
  resF = static_cast<expert_forest*>(rel->getOutForest());

  registerInForest(arg1F);
  registerInForest(arg2F);
  registerInForest(resF);
  ct_entry_type* et = new ct_entry_type(opcode->getName(), "N.N:N");
  et->setForestForSlot(0, arg1F);
  et->setForestForSlot(2, arg2F);
  et->setForestForSlot(4, resF);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::common_hyb_dfs_by_events_mt::~common_hyb_dfs_by_events_mt()
{
  if (rel->autoDestroy()) delete rel;
  unregisterInForest(arg1F);
  //unregisterInForest(arg2F);
  unregisterInForest(resF);
}



void MEDDLY::common_hyb_dfs_by_events_mt
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

  saturation_hyb_by_events_op* so = new saturation_hyb_by_events_op(this, arg1F, resF);

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
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

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

// ******************************************************************
// *                                                                *
// *               saturation_hyb_by_events_op  methods            *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_hyb_by_events_op
::saturation_hyb_by_events_op(common_hyb_dfs_by_events_mt* p,
                               expert_forest* argF, expert_forest* resF)
: unary_operation(saturation_hyb_by_events_opname::getInstance(), 1, argF, resF)
{
  parent = p;

  const char* name = saturation_hyb_by_events_opname::getInstance()->getName();
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
  removeAllComputeTableEntries();
}


MEDDLY::node_handle MEDDLY::saturation_hyb_by_events_op::saturate(MEDDLY::node_handle mdd)
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

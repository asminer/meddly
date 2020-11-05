
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
#include "sat_md_impl.h"
#include <typeinfo> // for "bad_cast" exception
#include <set>
#include <map>

   #define OUT_OF_BOUNDS -1
   #define NOT_KNOWN -2
   #define TERMINAL_NODE 1


namespace MEDDLY {
  class saturation_md_impl_by_events_opname;
  class saturation_md_impl_by_events_op;
  
  class common_md_impl_dfs_by_events_mt;
  class forwd_md_impl_dfs_by_events_mt;
};

// #define DEBUG_INITIAL
// #define DEBUG_IS_REACHABLE

// ******************************************************************
// *                                                                *
// *                     satmdimpl_opname  methods                  *
// *                                                                *
// ******************************************************************



MEDDLY::satmdimpl_opname::satmdimpl_opname(const char* n)
: specialized_opname(n)
{
}

MEDDLY::satmdimpl_opname::~satmdimpl_opname()
{
}


MEDDLY::gen_relation_node::gen_relation_node(std::string sign1, forest* f, int lvl, node_handle d, long* const_wgts, int* d_vars, int noof_dep_vars) : f(static_cast<expert_forest*>(f))
{
  signature  = sign1;
  level = lvl;
  down = d;
  if(const_wgts != NULL)
     {
       constant_wgts = (long*)malloc(sizeof(long)*3);
       constant_wgts[0] = const_wgts[0];
       constant_wgts[1] = const_wgts[1];
       constant_wgts[2] = const_wgts[2];
     }
  else
    constant_wgts = NULL;
  
  noof_deps = noof_dep_vars;
  dep_vars = (int*)malloc(sizeof(int)*noof_dep_vars);
  for(int i = 0; i < noof_deps; i++)
  {
    dep_vars[i] = d_vars[i];
  }
}

MEDDLY::gen_relation_node::~gen_relation_node()
{
}

std::set<std::vector<long>> MEDDLY::gen_relation_node::buildFreeValues(MEDDLY::expert_forest* resF)
{
  //to be defined for the example you use & comment this definition
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

long MEDDLY::gen_relation_node::delta(std::vector<long> v_in, std::vector<long> v_free, long i)
{
  //to be defined for the example you use & comment this definition
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

std::vector<long> MEDDLY::gen_relation_node::omega(std::vector<long> v_in, std::vector<long> v_free, long i)
{
  //to be defined for the example you use & comment this definition
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

bool
MEDDLY::gen_relation_node::equals(const gen_relation_node* n) const
{
  
  //bool const_flag = (constant_wgts[0] == n->getEnable()) && (constant_wgts[1] == n->getInhibit()) && (constant_wgts[2] == n->getFire());
  
  
  bool dep_flag = (noof_deps == n->getCountDep());
  
  if(dep_flag>0)
  {
    for(int i = 0; i < noof_deps; i++ )
    {
      dep_flag = dep_flag && (dep_vars[i] == n->getDepVariables()[i]);
      
      if(!dep_flag) break;
    }
  }
  
  if( (signature == n->getSignature()) && (level == n->getLevel()) && (down == n->getDown()) && dep_flag )
  {
    return true;
  }
  else
  {
    return false;
  }
}

// ******************************************************************

MEDDLY::satmdimpl_opname::md_implicit_relation::md_implicit_relation(forest* inmdd, forest* relmxd,
                                                             forest* outmdd)
: insetF(static_cast<expert_forest*>(inmdd)), outsetF(static_cast<expert_forest*>(outmdd)), relF(static_cast<expert_forest*>(relmxd))
{
  relF = static_cast<expert_forest*>(relmxd);
  
  if (0==insetF || 0==outsetF || 0==relF ) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  
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
      (insetF->getEdgeLabeling() != forest::MULTI_TERMINAL)   ||
      (outsetF->getEdgeLabeling() != forest::MULTI_TERMINAL)
      )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
  
  // Forests are good; set number of variables
  num_levels = insetF->getDomain()->getNumVariables();
  total_events = 0;
    
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
  /*gen_relation_node *Terminal = new MEDDLY::gen_relation_node(0, relF, 0, -1, NULL, NULL, 0);
  Terminal->setID(TERMINAL_NODE);
  last_in_node_array = TERMINAL_NODE;
  */
}


void
MEDDLY::satmdimpl_opname::md_implicit_relation::resizeEventArray(int level)
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
MEDDLY::satmdimpl_opname::md_implicit_relation::resizeConfirmedArray(int level, int index)
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

void findConfirmedStatesImpl(MEDDLY::satmdimpl_opname::md_implicit_relation* rel,
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

void MEDDLY::satmdimpl_opname::md_implicit_relation::setConfirmedStates(const dd_edge& set)
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
    findConfirmedStatesImpl(const_cast<md_implicit_relation*>(this),
                      confirmed, confirm_states, set.getNode(), num_levels, visited);
}



MEDDLY::satmdimpl_opname::md_implicit_relation::~md_implicit_relation()
{
  last_in_node_array = 0;
  
  for(int i = 0; i <=num_levels; i++) {delete[] event_list[i]; delete[] confirmed[i];}
  delete[] event_list;
  delete[] event_added;
  delete[] event_list_alloc;
  delete[] confirmed;
  delete[] confirm_states;
  delete[] confirmed_array_size;
}

int setOutputVariablesLocally(int* input_vars, int noof_inp, std::set<int> free_vars,
                        std::map<std::pair<int,int>, std::set<int>> &EL_set, int e_no, int lvl, int** output_vars) {
  
  std::set<int> op_set;
  
  for(auto i = free_vars.begin(); i != free_vars.end(); i++)
  {
    if(*i < lvl) op_set.insert(*i);
    else
    {
      auto it = EL_set[{e_no,*i}].find(lvl);
      if(it != EL_set[{e_no,*i}].end()) // curr_level is used, erase it
      {
        EL_set[{e_no,*i}].erase(it);
        if(!EL_set[{e_no,*i}].empty()) op_set.insert(*i);
      }
    }
  }
  
  for(int i = 0; i < noof_inp; i++)
  {
    if(input_vars[i] < lvl) op_set.insert(input_vars[i]);
    else if(input_vars[i] == lvl)
    {
      if(!EL_set[{e_no,input_vars[i]}].empty()) op_set.insert(input_vars[i]);
      
    }else
    {
      auto it = EL_set[{e_no,input_vars[i]}].find(lvl); // Am I using it and/or relaying it.
      if(it != EL_set[{e_no,input_vars[i]}].end()) // curr_level is used, erase it
      {
        EL_set[{e_no,input_vars[i]}].erase(it);
        if(!EL_set[{e_no,input_vars[i]}].empty()) op_set.insert(input_vars[i]);
      }
    }
  }
  

  if(!EL_set[{e_no,lvl}].empty()) op_set.insert(lvl);
  
  
  *output_vars = (int*)malloc(op_set.size()*sizeof(int));
  int w = 0;
  for(auto it = op_set.begin(); it!=op_set.end(); it++, w++)
  { (*output_vars)[w] = *it;
  }
    
  return op_set.size();
}

void
MEDDLY::satmdimpl_opname::md_implicit_relation::registerEventNodes(gen_relation_node** nb, int sz)
{
  //traverse the nodes and get an overall idea
  std::set<int> empty_set;
  empty_set.clear();
  rel_node_handle nh_next = -1;
  for(int i = 0;i < sz; i++)
  {
    nb[i]->setDown(nh_next); // assign down handle of this node
    nh_next = relF->createGenRelationNode(nb[i]); // register the node in the forest & get node_handle
    nb[i]->setID(nh_next); // assign node handle to this node
  }
  
  // for current event ordered at "total_events"
  // Assign each level to the set of levels where it has usage(empty for now)
  for(int i = 1; i<=num_levels; i++)
    event_level_affectlist.insert({std::make_pair(total_events,i), empty_set});
  
  for(int i = 0;i < sz; i++)
  {
    for(int d = 0; d < nb[i]->getCountDep(); d++)
    {
      event_level_affectlist[{total_events,nb[i]->getDepVariables()[d]}].insert(nb[i]->getLevel());
    }
  }
  //decide all params of the nodes
  for(int i = sz-1; i >= 0; i--)
  {
    int noof_dep = nb[i]->getCountDep();
    int curr_level = nb[i]->getLevel();
    
    if(i==sz-1) // top-node
    {
      int noof_inp = 0;
      int* input_vars = (int*)malloc(noof_inp*sizeof(int));
     
      std::set<int> free_set;
      for(int j = 0;j < noof_dep; j++)
      {
        int dep_level = nb[i]->getDepVariables()[j];
        if(dep_level < curr_level)
          free_set.insert(dep_level);
      }
      
      int* output_vars;
      int noof_op = setOutputVariablesLocally(input_vars, noof_inp, free_set,
                                       event_level_affectlist, total_events, curr_level, &output_vars);
     
      nb[i]->setInputVariables(input_vars, noof_inp);
      nb[i]->setOutputVariables(output_vars, noof_op);
      nb[i]->setFreeVariables(free_set);
    
    } else {
      
        int noof_inp = nb[i+1]->getCountOutput(); // Input of this is output of above
        int* input_vars = nb[i+1]->getOutputVariables();
        
        std::set<int> free_set;
        // Of all deps some are in input and some are free and one is self. Inputs are known. Free are the remaning ones.
        int j = 0; // dep
        int k = 0; // input
        
        
        while((j<noof_dep)&&(k<noof_inp)){
          if(nb[i]->getDepVariables()[j] == input_vars[k]){ // this dep var is specified in input & hence not free.
            k++;
            j++;
          } else if ((nb[i]->getDepVariables()[j] < input_vars[k]) && (nb[i]->getDepVariables()[j]!=nb[i]->getLevel())) { // this dep variable is a free.var
            free_set.insert(nb[i]->getDepVariables()[j]);
            j++;
          } else { // this input variable is not a dep.var & hence not free
            k++;
          }
        }
      
        if(j<noof_dep) // still more deps to read. & no more inp : these all free
        {
          while(j<noof_dep){
            if(nb[i]->getDepVariables()[j]!=nb[i]->getLevel()){
              free_set.insert(nb[i]->getDepVariables()[j]);
              j++;
            }
          }
        }
        
        int* output_vars;
        int noof_op = setOutputVariablesLocally(input_vars, noof_inp, free_set,
                                         event_level_affectlist, total_events, curr_level, &output_vars);
       
        nb[i]->setInputVariables(input_vars, noof_inp);
        nb[i]->setOutputVariables(output_vars, noof_op);
        nb[i]->setFreeVariables(free_set);
    }
    
  }
  
  
  int nLevel = nb[sz-1]->getLevel();
  resizeEventArray(nLevel);
  event_list[nLevel][event_added[nLevel] - 1] = nb[sz-1]->getID();
  total_events +=1;
  
}



// ******************************************************************
// *                                                                *
// *               saturation_impl_by_events_opname  class               *
// *                                                                *
// ******************************************************************

/** Simple class to keep compute table happy.
 */
class MEDDLY::saturation_md_impl_by_events_opname : public unary_opname {
  static saturation_md_impl_by_events_opname* instance;
public:
  saturation_md_impl_by_events_opname();
  
  static const saturation_md_impl_by_events_opname* getInstance();
  
};

MEDDLY::saturation_md_impl_by_events_opname* MEDDLY::saturation_md_impl_by_events_opname::instance = 0;

MEDDLY::saturation_md_impl_by_events_opname::saturation_md_impl_by_events_opname()
: unary_opname("Saturate_by_events")
{
}

const MEDDLY::saturation_md_impl_by_events_opname* MEDDLY::saturation_md_impl_by_events_opname::getInstance()
{
  if (0==instance) instance = new saturation_md_impl_by_events_opname;
  return instance;
}

// ******************************************************************
// *                                                                *
// *             saturation_impl_by_events_op  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::saturation_md_impl_by_events_op : public unary_operation {
  common_md_impl_dfs_by_events_mt* parent;
public:
  saturation_md_impl_by_events_op(common_md_impl_dfs_by_events_mt* p,
                               expert_forest* argF, expert_forest* resF);
  virtual ~saturation_md_impl_by_events_op();
  
  node_handle saturate(node_handle mdd);
  node_handle saturate(node_handle mdd, int level);


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
// *            common_md_impl_dfs_by_events_mt  class                 *
// *                                                                *
// ******************************************************************

class MEDDLY::common_md_impl_dfs_by_events_mt : public specialized_operation {
public:
  common_md_impl_dfs_by_events_mt(const satmdimpl_opname* opcode,
                               satmdimpl_opname::md_implicit_relation* rel);
  virtual ~common_md_impl_dfs_by_events_mt();
  
  virtual void compute(const dd_edge& a, dd_edge &c);
  virtual void saturateHelper(unpacked_node& mdd) = 0;
  
protected:
  inline compute_table::entry_key*
  findResult(node_handle a, rel_node_handle b, std::vector<long> input_values, node_handle &c)
  {
    compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], input_values.size());
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeN(a);
    CTsrch->writeN(b);
    for(int i = 0; i<input_values.size(); i++)
    {
      CTsrch->writeL(input_values[i]);
    }
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
                                node_handle a, rel_node_handle b, std::vector<long> input_values, node_handle c)
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
  
  satmdimpl_opname::md_implicit_relation* rel;
  
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
// *               forwd_md_impl_dfs_by_events_mt class                *
// *                                                                *
// ******************************************************************


class MEDDLY::forwd_md_impl_dfs_by_events_mt : public common_md_impl_dfs_by_events_mt {
public:
  forwd_md_impl_dfs_by_events_mt(const satmdimpl_opname* opcode,
                              satmdimpl_opname::md_implicit_relation* rel);
protected:
  virtual void saturateHelper(unpacked_node& mdd);
  node_handle recFire(node_handle mdd, rel_node_handle mxd, std::vector<long> inputs);
  
 };

MEDDLY::forwd_md_impl_dfs_by_events_mt::forwd_md_impl_dfs_by_events_mt(
                                                                 const MEDDLY::satmdimpl_opname* opcode,
                                                                 MEDDLY::satmdimpl_opname::md_implicit_relation* rel)
: common_md_impl_dfs_by_events_mt(opcode, rel)
{
}

long extractOwnValue(MEDDLY::gen_relation_node* relNode, std::vector<long> input_values)
{
  int node_level = relNode->getLevel();
  int no_of_inp = relNode->getCountInput();
  int* ip_vars = relNode->getInputVariables();
  for(int i = 0; i< no_of_inp; i++)
  {
    if(node_level==ip_vars[i])
      return input_values[i];
    else if(node_level>ip_vars[i])
      return -1;
  }
  
  return -1;
}


void MEDDLY::forwd_md_impl_dfs_by_events_mt::saturateHelper(unpacked_node& nb)
{
  int nEventsAtThisLevel = rel->lengthForLevel(nb.getLevel());
  if (0 == nEventsAtThisLevel) return;
  
  // Initialize mxd readers, note we might skip the unprimed level
  const int level = nb.getLevel();
  node_handle* events = rel->arrayForLevel(level);
  gen_relation_node** Ru = new gen_relation_node*[nEventsAtThisLevel];
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    Ru[ei] = arg2F->buildGenImplicitNode(events[ei]);
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
  
  //i' = i + a1 - a2
  
  // explore indexes
  while (!queue->isEmpty()) {
    int i = queue->remove();
    MEDDLY_DCASSERT(nb.d(i));
    
      for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
        // Should this be calculated several times? or Can i precalculate and re-use it?
        // Do I use the full cartesian_product or just the news ones?
        //std::set<std::vector<long>> cart_prod = calculateCartesianProd(resF, Ru[ei]->getFreeVariables(), Ru[ei]->getCountFree(), rel);
        std::set<std::vector<long>> cart_prod = Ru[ei]->buildFreeValues(resF);

        for(auto guess_it = cart_prod.begin(); guess_it != cart_prod.end(); guess_it++){
          std::vector<long> empty_vector;
          
          int j = Ru[ei]->delta(empty_vector, *guess_it, i);
          if(j==-1) continue;
          if (j < nb.getSize() && -1==nb.d(j)) continue; // nothing can be added to this set
          
          node_handle rec = recFire(nb.d(i), Ru[ei]->getDown(), Ru[ei]->omega(empty_vector, *guess_it, i));
          
          std::vector<long> whatis =Ru[ei]->omega(empty_vector, *guess_it, i);
          
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
        } // free variables' guesses
      } // for all events, ei
    
  
 }// more indexes to explore
  
  delete[] Ru;
  recycle(queue);
}


// Same as post-image, except we saturate before reducing.
MEDDLY::node_handle MEDDLY::forwd_md_impl_dfs_by_events_mt::recFire(
                                                                 MEDDLY::node_handle mdd, rel_node_handle mxd, std::vector<long> input_values)
{
  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  
  if (mxd==-1) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->handleForValue(1);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }
  
  gen_relation_node* relNode = arg2F->buildGenImplicitNode(mxd); // The relation node
  
  // check the cache
  node_handle result = 0;
  compute_table::entry_key* Key = findResult(mdd, mxd, input_values, result);  // Place the inputs into this CT.
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
      nb->d_ref(i) = recFire(A->d(i), mxd, input_values);
    }
    
  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(mxdLevel >= mddLevel);
    
    // does input_values have my own value??
    int own_value = extractOwnValue(relNode, input_values);
    
    // Initialize mxd readers, note we might skip the unprimed level
    
    // loop over mxd "rows" if the value is not specified via input
        for (int iz=0; iz<rSize; iz++) {
          if((own_value != -1) && (own_value!=iz)) continue;
          int i = iz; // relation_node enabling condition
          if (0==A->d(i))   continue;
          
          // loop over mxd "columns"
          // loop over free values
          //std::set<std::vector<long>> cart_prod = calculateCartesianProd(resF, relNode->getFreeVariables(), relNode->getCountFree(), rel);
          std::set<std::vector<long>> cart_prod = relNode->buildFreeValues(resF);
          for(auto guess_it = cart_prod.begin(); guess_it != cart_prod.end(); guess_it++){
            
          int j = relNode->delta(input_values, *guess_it, i);
          if(j==-1) continue;
          
          node_handle newstates = recFire(A->d(i), relNode->getDown(), relNode->omega(input_values, *guess_it, i));
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
          
          } // for free-variables
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
  
  
  
  return saveResult(Key, mdd, mxd, input_values, result);
}

// ******************************************************************
// *                                                                *
// *             common_impl_dfs_by_events_mt  methods              *
// *                                                                *
// ******************************************************************

MEDDLY::common_md_impl_dfs_by_events_mt::common_md_impl_dfs_by_events_mt(
                                                                   const satmdimpl_opname* opcode,
                                                                   satmdimpl_opname::md_implicit_relation* relation)
: specialized_operation(opcode, 1)
{
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;
  freeqs = 0;
  freebufs = 0;
  rel = relation;
  arg1F = static_cast<expert_forest*>(rel->getInForest());
  arg2F = static_cast<expert_forest*>(rel->getRelForest());
  resF = static_cast<expert_forest*>(rel->getOutForest());
  
  registerInForest(arg1F);
  //registerInForest(arg2F);
  registerInForest(resF);
  compute_table::entry_type* et = new compute_table::entry_type(opcode->getName(), "NN.L:N");
  et->setForestForSlot(0, arg1F);
  et->setForestForSlot(1, arg2F);
  et->setForestForSlot(5, resF);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::common_md_impl_dfs_by_events_mt::~common_md_impl_dfs_by_events_mt()
{
  if (rel->autoDestroy()) delete rel;
  unregisterInForest(arg1F);
  unregisterInForest(arg2F);
  unregisterInForest(resF);
}

void MEDDLY::common_md_impl_dfs_by_events_mt
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
  
  saturation_md_impl_by_events_op* so = new saturation_md_impl_by_events_op(this, arg1F, resF);
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
// *       common_md_impl_dfs_by_events_mt::indexq  methods                 *
// ******************************************************************

MEDDLY::common_md_impl_dfs_by_events_mt::indexq::indexq()
{
  data = 0;
  size = 0;
  head = NULPTR;
}

MEDDLY::common_md_impl_dfs_by_events_mt::indexq::~indexq()
{
  free(data);
}

void MEDDLY::common_md_impl_dfs_by_events_mt::indexq::resize(int sz)
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

MEDDLY::common_md_impl_dfs_by_events_mt::charbuf::charbuf()
{
  data = 0;
  size = 0;
}

MEDDLY::common_md_impl_dfs_by_events_mt::charbuf::~charbuf()
{
  free(data);
}

void MEDDLY::common_md_impl_dfs_by_events_mt::charbuf::resize(int sz)
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

MEDDLY::satmdimpl_opname* MEDDLY::initMDImplSaturationForward()
{
  return new satmdimpl_opname("SaturationFwd");
}


MEDDLY::specialized_operation*
MEDDLY::satmdimpl_opname::buildOperation(arguments* a) const
{
  
  md_implicit_relation* rel = dynamic_cast<md_implicit_relation*>(a);
  if (0==rel) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  
  MEDDLY::specialized_operation* op = 0;
  op = new forwd_md_impl_dfs_by_events_mt(this, rel);
  
  return op;
}

// ******************************************************************
// *                                                                *
// *               saturation_impl_by_events_op  methods            *
// *                                                                *
// ******************************************************************

MEDDLY::saturation_md_impl_by_events_op
::saturation_md_impl_by_events_op(common_md_impl_dfs_by_events_mt* p,
                               expert_forest* argF, expert_forest* resF)
: unary_operation(saturation_md_impl_by_events_opname::getInstance(), 1, argF, resF)
{
  parent = p;

  const char* name = saturation_md_impl_by_events_opname::getInstance()->getName();
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

MEDDLY::saturation_md_impl_by_events_op::~saturation_md_impl_by_events_op()
{
  removeAllComputeTableEntries();
}


MEDDLY::node_handle MEDDLY::saturation_md_impl_by_events_op::saturate(MEDDLY::node_handle mdd)
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
MEDDLY::saturation_md_impl_by_events_op::saturate(node_handle mdd, int k)
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
  parent->saturateHelper(*nb);
  n = resF->createReducedNode(-1, nb);
  
  // save in compute table
  saveSaturateResult(Key, mdd, n);
  
  #ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
  #endif
  
  return n;
}


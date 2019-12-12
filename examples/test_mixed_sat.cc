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


/*
 State space generation using On-The-Fly Saturation.
 
 Model: Swimming Pool Petri Net
 
 Places: A, B, C, D, F, G, H, I
 Transitions: GetK, GetB, RelK, GetK2, RBag, RKey, Enter;
 
 Initial state: G = 20*N, H = 15*N, I = 10*N
*/

#if 1
#include "../src/meddly.h"
#include "../src/meddly_expert.h"

using namespace MEDDLY;

class pn;
class enabling_subevent;
class firing_subevent;
class derived_relation_node;
int usage(const char*);

class pn {
public:
  pn(int nTokens);
  ~pn();
  dd_edge& getReachableStateSet();
  void clearReachableStateSet();
  
private:
  void initializeMeddly();
  void buildDomain();
  void buildMdd();
  void buildMxd();
  void buildTransition();
  void buildInitialState();
  void buildOtfSaturationOp();
  void doOtfSaturation();
  
  bool addValue(int vh, int value, int& index);
  bool indexOf(int vh, int value, int &index);
  bool valueOf(int vh, int index, int &value);
  
  satimpl_opname::implicit_relation* impl_rel;
  
  satimpl_opname::event* GetK;
  satimpl_opname::event* GetB;
  satimpl_opname::event* RelK;
  satimpl_opname::event* GetK2;
  satimpl_opname::event* RBag;
  satimpl_opname::event* RKey;
  satimpl_opname::event* Enter;
  
  satimpl_opname::subevent* GetK_1;
  satimpl_opname::subevent* GetK_9;
  relation_node* GetK_2;
  
  satimpl_opname::subevent* GetB_2;
  satimpl_opname::subevent* GetB_8;
  relation_node* GetB_3;
  
  satimpl_opname::subevent* RelK_3;
  relation_node* RelK_4;
  relation_node* RelK_9;
  
  satimpl_opname::subevent* GetK2_4;
  satimpl_opname::subevent* GetK2_9;
  relation_node* GetK2_5;
  
  satimpl_opname::subevent* RBag_5;
  relation_node* RBag_6;
  relation_node* RBag_8;
  
  satimpl_opname::subevent* RKey_6;
  relation_node* RKey_7;
  relation_node* RKey_9;
  
  satimpl_opname::subevent* Enter_7;
  relation_node* Enter_1;
  
  
  
  specialized_operation* impl_sat_op;
  
  int n_tokens;
  domain* dom;
  expert_forest* mdd;
  expert_forest* mxd;
  dd_edge* initial_state;
  dd_edge* reachable_states;
  
  const int n_vars = 9;
  const int A_level = 1;
  const int B_level = 2;
  const int C_level = 3;
  const int D_level = 4;
  const int E_level = 5;
  const int F_level = 6;
  const int G_level = 7;
  const int H_level = 8;
  const int I_level = 9;
  
  
  // index to value
  std::vector<std::vector<int>> index_to_value;
  std::vector<std::vector<int>> value_to_index;
  
  friend enabling_subevent;
  friend firing_subevent;
  friend derived_relation_node;
  //friend derived_event;
};

/*class derived_event: public satimpl_opname::event {
 public:
 derived_event(pn* model, satimpl_opname::subevent** se, int n_subevent_vars);
 ~derived_event();
 virtual int getValue(int level, int i);
 virtual int addValue(int level, int i);
 protected:
 pn* model;
 };*/

class enabling_subevent : public satimpl_opname::subevent {
public:
  enabling_subevent(pn* model, forest* mxd,
                    int* subevent_vars, int n_subevent_vars,
                    int* event_vars, int n_event_vars);
  ~enabling_subevent();
  virtual void confirm(satimpl_opname::implicit_relation &rel, int v, int index);
  void initializeMinterm();
protected:
  pn* model;
  forest* mxd;
  int subevent_var;
  int* event_vars;
  int n_event_vars;
  int* unp_minterm;
  int* p_minterm;
  int minterm_size;
};

class firing_subevent : public satimpl_opname::subevent {
public:
  firing_subevent(pn* model, forest* mxd,
                  int* subevent_vars, int n_subevent_vars,
                  int* event_vars, int n_event_vars,
                  int eoffset, int foffset);
  ~firing_subevent();
  virtual void confirm(satimpl_opname::implicit_relation &rel, int v, int index);
  void initializeMinterm();
protected:
  pn* model;
  forest* mxd;
  int subevent_var;
  int* event_vars;
  int n_event_vars;
  int eoffset;
  int foffset;  
  int* unp_minterm;
  int* p_minterm;
  int minterm_size;
};


class derived_relation_node : public relation_node{
public:
  derived_relation_node(pn* model, forest* mxd,
                        int rn_var,
                        int en_offset, int f_offset);
  ~derived_relation_node();
  //virtual void confirm(satimpl_opname::implicit_relation &rel, int v, int index);
  virtual long nextOf(long i);
protected:
  pn* model;
  forest* f;
  int rn_var;
  int en_offset;
  int f_offset;
  
};

int main(int argc, char* argv[]) {
  int nTokens = 0;  // Number of tokens
  
  for (int i = 1; i < argc; i++) {
    const char* cmd = argv[i];
    if (strncmp(cmd, "-n", 2) == 0) {
      nTokens = strtol(cmd+2, NULL, 10);
      if (nTokens < 1) {
        return 1+usage(argv[0]);
      }
      continue;
    }
    return 1+usage(argv[0]);
  }
  
  if (nTokens < 1) return 1+usage(argv[0]);
  std::cout << "Tokens: " << nTokens << std::endl;
  
  pn model(nTokens);
  dd_edge rss = model.getReachableStateSet();
  ostream_output s(std::cout);
  rss.show(s, 2);
  double rss_card = rss.getCardinality();
  std::cout << "\n\nNumber of states: " << rss_card << std::endl;
}


int usage(const char* who)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) name = ptr+1;
  }
  
  printf("\nUsage: %s [options]\n\n", name);
  
  printf("\t-n<#tokens>: set number of tokens\n\n");
  printf("\n");
  return 0;
}

// -------------- class pn implementation ---------------

pn::pn(int n_tokens) {
  // must:
  //    initialize meddly
  //    create a domain with vars A and B
  //    create mdd forest
  //    create mxd forest
  //    build T_ab
  
  this->n_tokens = n_tokens;
  
  impl_rel = 0;
  GetK_1 = 0;
  GetK_9 = 0;
  GetK_2 = 0;
  
  GetB_2 = 0;
  GetB_8 = 0;
  GetB_3 = 0;
  
  RelK_3 = 0;
  RelK_4 = 0;
  RelK_9 = 0;
  
  GetK2_4 = 0;
  GetK2_9 = 0;
  GetK2_5 = 0;
  
  RBag_5 = 0;
  RBag_6 = 0;
  RBag_8 = 0;
  
  RKey_6 = 0;
  RKey_7 = 0;
  RKey_9 = 0;
  
  Enter_7 = 0;
  Enter_1 = 0;
  
  /*************/  
  
  GetK = 0;
  GetB = 0;
  RelK = 0;
  GetK2 = 0;
  RBag = 0;
  RKey = 0;
  Enter = 0;
  impl_sat_op = 0;
  
  dom = 0;
  mdd = 0;
  mxd = 0;
  initial_state = 0;
  reachable_states = 0;
  
  // initialize meddly
  initializeMeddly();
  
  // create a domain with vars A and B
  buildDomain();
  
  // create mdd forest
  buildMdd();
  
  // create mxd forest
  buildMxd();
  
  // build an initial state with A = n_tokens, B = 0
  buildInitialState();
  
  // build T_ab
  buildTransition();
  
  // clear the set of reachable states
  clearReachableStateSet();
}

pn::~pn() {
  // must:
  // delete the otf relation
  // delete the events
  // delete the domain
  
  // delete the otf relation
  // delete otf_rel;
  
  // delete the events
  delete GetK;
  delete GetB;
  delete RelK;
  delete GetK2;
  delete RBag;
  delete RKey;
  delete Enter;
  
  
  // delete the domain
  MEDDLY::cleanup();
}

void pn::initializeMeddly() {
  MEDDLY::initialize();
}

void pn::buildDomain() {
  int bounds[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  dom = createDomainBottomUp((int*)&bounds[0], n_vars);
  assert(dom);
}

void pn::buildMdd() {
  forest::policies p(false);
  p.setQuasiReduced();
  mdd = static_cast<expert_forest*>(dom->createForest(false, forest::BOOLEAN, forest::MULTI_TERMINAL, p));
  assert(mdd);
}

void pn::buildMxd() {
  forest::policies p(true);
  p.setFullStorage();
  mxd = static_cast<expert_forest*>(dom->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL, p));
  assert(mxd);
}

void pn::buildTransition() {
  // Only one transition in the system: T_ab: A--, B--, C++
  // build A-- subevent
  // build C++ subevent
  // build B-- relation node
  // build array of variables this transition depends on
  
  //if (A_ABC_firing) delete A_ABC_firing;
  //if (A_ABC_enabling) delete A_ABC_enabling;
  //if (C_ABC_firing) delete C_ABC_firing;
  //if (T_ABC) delete T_ABC;
  
  int  dep_GetK_1[] = {A_level};
  int  dep_GetK_9[] = {I_level};
  int  dep_GetK_2[] = {B_level};
  int dep_GetK[] = {A_level, B_level, I_level};
  
  int  dep_GetB_2[] = {B_level};
  int  dep_GetB_8[] = {H_level};
  int  dep_GetB_3[] = {C_level};
  int dep_GetB[] = {B_level, C_level, H_level};
  
  int  dep_RelK_3[] = {C_level};
  int  dep_RelK_4[] = {D_level};
  int  dep_RelK_9[] = {I_level};
  int dep_RelK[] = {C_level, D_level, I_level};  
  
  int  dep_GetK2_4[] = {D_level};
  int  dep_GetK2_9[] = {I_level};
  int  dep_GetK2_5[] = {E_level};
  int dep_GetK2[] = {D_level, E_level, I_level};
  
  int  dep_RBag_5[] = {E_level};
  int  dep_RBag_6[] = {F_level};
  int  dep_RBag_8[] = {H_level};
  int dep_RBag[] = {E_level, F_level, H_level};
  
  int  dep_RKey_6[] = {F_level};
  int  dep_RKey_7[] = {G_level};
  int  dep_RKey_9[] = {I_level};
  int dep_RKey[] = {F_level, G_level, I_level};
  
  int  dep_Enter_7[] = {G_level};
  int  dep_Enter_1[] = {A_level};
  int dep_Enter[] = {A_level, G_level};  
  
  GetK_1 = new firing_subevent(this, mxd, (int*)(&dep_GetK_1[0]), 1, (int*)(&dep_GetK[0]), 3, 1, -1);
  GetK_9 = new firing_subevent(this, mxd, (int*)(&dep_GetK_9[0]), 1, (int*)(&dep_GetK[0]), 3, 1, -1);
  GetK_2 = new derived_relation_node(this, mxd, dep_GetK_2[0], 0, 1);
  satimpl_opname::subevent* subevents_GetK[] = {GetK_1, GetK_9};
  relation_node* relNodes_GetK[] = {GetK_2};
  GetK = new satimpl_opname::event((satimpl_opname::subevent**)(&subevents_GetK[0]), 2, (relation_node**)(&relNodes_GetK[0]), 1);
  
  
  GetB_2 = new firing_subevent(this, mxd, (int*)(&dep_GetB_2[0]), 1, (int*)(&dep_GetB[0]), 3, 1, -1);
  GetB_8 = new firing_subevent(this, mxd, (int*)(&dep_GetB_8[0]), 1, (int*)(&dep_GetB[0]), 3, 1, -1); 
  GetB_3 = new derived_relation_node(this, mxd, dep_GetB_3[0], 0, 1);
  satimpl_opname::subevent* subevents_GetB[] = {GetB_2,GetB_8};
  relation_node* relNodes_GetB[] = {GetB_3};
  GetB = new satimpl_opname::event((satimpl_opname::subevent**)(&subevents_GetB[0]), 2, (relation_node**)(&relNodes_GetB[0]), 1);
  
  
  RelK_3 = new firing_subevent(this, mxd, (int*)(&dep_RelK_3[0]), 1, (int*)(&dep_RelK[0]), 3, 1, -1); 
  RelK_4 = new derived_relation_node(this, mxd, dep_RelK_4[0], 0, 1);
  RelK_9 = new derived_relation_node(this, mxd, dep_RelK_9[0], 0, 1);
  satimpl_opname::subevent* subevents_RelK[] = {RelK_3};
  relation_node* relNodes_RelK[] = {RelK_4, RelK_9};
  RelK = new satimpl_opname::event((satimpl_opname::subevent**)(&subevents_RelK[0]), 1, (relation_node**)(&relNodes_RelK[0]), 2);
  
  
  GetK2_4 = new firing_subevent(this, mxd, (int*)(&dep_GetK2_4[0]), 1, (int*)(&dep_GetK2[0]), 3, 1, -1); 
  GetK2_9 = new firing_subevent(this, mxd, (int*)(&dep_GetK2_9[0]), 1, (int*)(&dep_GetK2[0]), 3, 1, -1);
  GetK2_5 = new derived_relation_node(this, mxd, dep_GetK2_5[0], 0, 1);
  satimpl_opname::subevent* subevents_GetK2[] = {GetK2_4,GetK2_9};
  relation_node* relNodes_GetK2[] = {GetK2_5};
  GetK2 = new satimpl_opname::event((satimpl_opname::subevent**)(&subevents_GetK2[0]), 2, (relation_node**)(&relNodes_GetK2[0]), 1);
  
  
  RBag_5 = new firing_subevent(this, mxd, (int*)(&dep_RBag_5[0]), 1, (int*)(&dep_RBag[0]), 3, 1, -1);
  RBag_6 = new derived_relation_node(this, mxd, dep_RBag_6[0], 0, 1);
  RBag_8 = new derived_relation_node(this, mxd, dep_RBag_8[0], 0, 1);
  satimpl_opname::subevent* subevents_RBag[] = {RBag_5};
  relation_node* relNodes_RBag[] = {RBag_6, RBag_8};
  RBag = new satimpl_opname::event((satimpl_opname::subevent**)(&subevents_RBag[0]), 1, (relation_node**)(&relNodes_RBag[0]), 2);
  
  
  RKey_6 = new firing_subevent(this, mxd, (int*)(&dep_RKey_6[0]), 1, (int*)(&dep_RKey[0]), 3, 1, -1);
  RKey_7 = new derived_relation_node(this, mxd, dep_RKey_7[0], 0, 1);
  RKey_9 = new derived_relation_node(this, mxd, dep_RKey_9[0], 0, 1);
  satimpl_opname::subevent* subevents_RKey[] = {RKey_6};
  relation_node* relNodes_RKey[] = {RKey_7, RKey_9};
  RKey = new satimpl_opname::event((satimpl_opname::subevent**)(&subevents_RKey[0]), 1, (relation_node**)(&relNodes_RKey[0]), 2);
  
  
  Enter_7 = new firing_subevent(this, mxd, (int*)(&dep_Enter_7[0]), 1, (int*)(&dep_Enter[0]), 2, 1, -1);
  Enter_1 = new derived_relation_node(this, mxd, dep_Enter_1[0], 0, 1);
  satimpl_opname::subevent* subevents_Enter[] = {Enter_7};
  relation_node* relNodes_Enter[] = {Enter_1};
  Enter = new satimpl_opname::event((satimpl_opname::subevent**)(&subevents_Enter[0]), 1, (relation_node**)(&relNodes_Enter[0]), 1);
  
  
  
}

void pn::buildInitialState() {
  initial_state = new dd_edge(mdd);
  int init_state[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int* init[] = {init_state};
  mdd->createEdge((int**)(&init[0]), 1, *initial_state);
  assert(initial_state);
  ostream_output s(std::cout);
  initial_state->show(s, 2);
  int index = -1;
  assert(addValue(A_level, 0, index));
  assert(addValue(B_level, 0, index));
  assert(addValue(C_level, 0, index));
  assert(addValue(D_level, 0, index));
  assert(addValue(E_level, 0, index));
  assert(addValue(F_level, 0, index));
  assert(addValue(G_level, 20*n_tokens, index));
  assert(addValue(H_level, 15*n_tokens, index));
  assert(addValue(I_level, 10*n_tokens, index));
}

void pn::buildOtfSaturationOp() {
  if (impl_sat_op == 0) {
    if (impl_rel == 0) {
      // build otf saturation relation
      satimpl_opname::event* events[] = {GetK, GetB, RelK, GetK2, RBag, RKey, Enter};
      forest* mxdd = mxd;
      impl_rel = new satimpl_opname::implicit_relation(mdd, mxd, mdd, &events[0], 7);
    }
    assert(impl_rel);
    impl_sat_op = SATURATION_IMPL_FORWARD->buildOperation(impl_rel);
    assert(impl_sat_op);
  }
}

void pn::doOtfSaturation() {
  if (impl_sat_op == 0) buildOtfSaturationOp();
  assert(impl_sat_op);
  
  if (reachable_states) delete reachable_states;
  reachable_states = new dd_edge(mdd);
  assert(reachable_states);
  
  impl_rel->setConfirmedStates(*initial_state);
  impl_sat_op->compute(*initial_state, *reachable_states);
  assert(reachable_states);
}

dd_edge& pn::getReachableStateSet() {
  if (reachable_states == 0) doOtfSaturation();
  assert(reachable_states);
  return *reachable_states;
}

void pn::clearReachableStateSet() {
  if (reachable_states) delete reachable_states;
  reachable_states = 0;
}

bool pn::indexOf(int vh, int value, int &index) {
  
  if (value_to_index.size() <= unsigned(vh)) return false;
  if (value_to_index[vh].size() <= unsigned(value)) return false;
  index = value_to_index[vh][value];
  return true;
}

bool pn::valueOf(int vh, int index, int &value) {
  if (index_to_value.size() <= unsigned(vh)) return false;
  if (index_to_value[vh].size() <= unsigned(index)) return false;
  value = index_to_value[vh][index];
  return true;
}

bool pn::addValue(int vh, int value, int &index) {
  if (indexOf(vh, value, index) && index >= 0) return false;
  
  if (value_to_index.size() <= unsigned(vh)) value_to_index.resize(vh+1);
  if (index_to_value.size() <= unsigned(vh)) index_to_value.resize(vh+1);
  if (value_to_index[vh].size() <= unsigned(value)) {
    for (int j = (value + 1 - value_to_index[vh].size()); j > 0; j--)
      value_to_index[vh].push_back(-1);
  }
  
  // get next available index, and map value to that index
  index = index_to_value[vh].size();
  value_to_index[vh][value] = index;
  index_to_value[vh].resize(index+1);
  index_to_value[vh][index] = value;
  
  return true;
}

// ----------------- end of class pn implementation -----------------

// -------------- class enabling_subevent implementation ------------

enabling_subevent::enabling_subevent(
                                     pn* model, forest* mxd,
                                     int* se_vars, int n_se_vars,  // subevent vars
                                     int* e_vars, int n_e_vars)    // event vars
: satimpl_opname::subevent(mxd, se_vars, n_se_vars, false, 0, 0)
{
  this->model = model;
  this->mxd = mxd;
  assert(n_se_vars == 1);
  subevent_var = se_vars[0];
  n_event_vars = n_e_vars;
  event_vars = new int[n_event_vars];
  for (int i = 0; i < n_event_vars; i++) event_vars[i] = e_vars[i];
  minterm_size = 0;
  unp_minterm = 0;
  p_minterm = 0;
}


enabling_subevent::~enabling_subevent() {
  delete [] event_vars;
  delete [] unp_minterm;
  delete [] p_minterm;
}


void enabling_subevent::initializeMinterm() {
  assert(minterm_size == 0);
  minterm_size = mxd->getDomain()->getNumVariables() + 1;
  unp_minterm = new int[minterm_size];
  p_minterm = new int[minterm_size];
  for (int i = 1; i < minterm_size; i++) {
    unp_minterm[i] = DONT_CARE;
    p_minterm[i] = DONT_CHANGE;
  }
  for (int i = 0; i < n_event_vars; i++) {
    p_minterm[event_vars[i]] = DONT_CARE;
  }
}

void enabling_subevent::confirm(satimpl_opname::implicit_relation &rel, int v, int index) {
  std::cout << "Confirming (enabling): (" << v << ", " << index << ")\n";
  // add minterm:
  // if (value(v[index]) > 0)
  //    add minterm with, minterm[v] = index, minterm[v'] = dont-care
  MEDDLY_DCASSERT(v == subevent_var);
  if (minterm_size == 0) initializeMinterm();
  int value = 0;
  if (!model->valueOf(v, index, value)) {
    // unknown index!
    assert(false);
  }
  if (value > 0) {
    // add enabling minterm
    unp_minterm[v] = index;
    p_minterm[v] = DONT_CARE;
    addMinterm(unp_minterm, p_minterm);
  }
}


// ---------- end of class enabling_subevent implementation ---------

// -------------- class firing_subevent implementation --------------


firing_subevent::firing_subevent(
                                 pn* model, forest* mxd,
                                 int* se_vars, int n_se_vars,  // subevent vars
                                 int* e_vars, int n_e_vars,    // event vars
                                 int eoffset, int foffset)
: satimpl_opname::subevent(mxd, se_vars, n_se_vars, true, eoffset, foffset)
{
  this->model = model;
  this->mxd = mxd;
  assert(n_se_vars == 1);
  subevent_var = se_vars[0];
  n_event_vars = n_e_vars;
  event_vars = new int[n_event_vars];
  for (int i = 0; i < n_event_vars; i++) event_vars[i] = e_vars[i];
  this->eoffset = eoffset;
  this->foffset = foffset;
  minterm_size = 0;
  unp_minterm = 0;
  p_minterm = 0;
}


firing_subevent::~firing_subevent() {
  delete [] event_vars;
  delete [] unp_minterm;
  delete [] p_minterm;
}


void firing_subevent::initializeMinterm() {
  assert(minterm_size == 0);
  minterm_size = mxd->getDomain()->getNumVariables() + 1;
  unp_minterm = new int[minterm_size];
  p_minterm = new int[minterm_size];
  unp_minterm[0] = 0;p_minterm[0] = 0;
  for (int i = 1; i < minterm_size; i++) {
    unp_minterm[i] = DONT_CARE;
    p_minterm[i] = DONT_CARE;
  }
  for (int i = 0; i < n_event_vars; i++) {
    p_minterm[event_vars[i]] = DONT_CHANGE;
  }
}

void firing_subevent::confirm(satimpl_opname::implicit_relation &rel, int v, int index) {
  std::cout << "Confirming (firing): (" << v << ", " << index << ")\n";
  // add minterm:
  // if it is a known index
  //    find value, next_value, and next_value_index
  //    add minterm with, minterm[v] = index, minterm[v'] = next_value_index
  MEDDLY_DCASSERT(v == subevent_var);
  if (minterm_size == 0) initializeMinterm();
  int value = 0;
  if (!model->valueOf(v, index, value)) {
    // unknown index!
    assert(false);
  }
  
  if (value < eoffset) 
    {
    printf("\n Cannot add this minterm because token(%d) is not atleast >= enable_offset(%d)", value, eoffset);
    return;
    }
  int next_value = value + foffset;
  
  
  int next_value_index = -1;
/*#ifdef DEVELOPMENT_CODE
  assert(model->addValue(v, next_value, next_value_index));
#else*/
  model->addValue(v, next_value, next_value_index);
//#endif
  // add enabling minterm
  unp_minterm[v] = index;
  p_minterm[v] = next_value_index;
  printf("\n AT level %d : minterm from=%d, to=%d",v,index,next_value_index); 
  addMinterm(unp_minterm, p_minterm);
  std::cout << "\t\tk: " << v
  << ", index: " << index
  << ", value: " << value
  << ", next_value: " << next_value
  << ", next_value_index: " << next_value_index
  << "\n";
}

// ----------- end of class firing_subevent implementation ----------

/*
 derived_event::derived_event(pn* m, satimpl_opname::subevent** se, int nse) : satimpl_opname::event((satimpl_opname::subevent**) (se) ,nse) {
 model = m;
 }
 derived_event::~derived_event() {}
 
 int derived_event::getValue(int v, int index) {
 int value = 0;
 model->valueOf(v, index, value);
 return value;
 }
 
 int derived_event::addValue(int v, int value) {
 int index = -1;
 model->addValue(v, value, index);
 return index;
 }*/


// ------------ class derived_relation_node implmentation -----------

derived_relation_node::derived_relation_node(
                                             pn* model, forest* f,
                                             int rn_var,  // relation node var
                                             int en_offset, int f_offset)
: relation_node(141010, f, rn_var, -1, en_offset, f_offset)
//relation_node(unsigned long signature, forest* f, int level, node_handle down, long e_val, long f_val);
{
  this->model = model;
  this->f = f;
  this->rn_var = rn_var;
  this->en_offset = en_offset;
  this->f_offset = f_offset;
}


derived_relation_node::~derived_relation_node() {
}

long derived_relation_node::nextOf(long i) {
  
  int val_at_i = -1;
  model->valueOf(rn_var, i, val_at_i);
  assert(val_at_i!=-1);
  
  /*if(val_at_i>=getPieceSize()) //Array needs to be allocated
   expandTokenUpdate(val_at_i);*/
  
  long val_at_j = val_at_i >= en_offset? val_at_i + f_offset : -1;
  printf("\n To Token: %d",val_at_j);
  
  if(val_at_j<0) return -1;
  
  int j;
  if(!model->addValue(rn_var, val_at_j, j))
    model->indexOf(rn_var, val_at_j, j);
  printf(", Index: %d",j);
  
  /* if(getTokenUpdate()[val_at_i]==-2) //Array needs to be updated
   {  
   int val = val_at_j<0? -1: val_at_j;
   setTokenUpdateAtIndex(i,val_at_j);
   }
   
   return getTokenUpdate()[val_at_i];*/
  
  return j;
  
  /*
   if (val<0)
   return -1;
   else
   return val;
   */
  /*
   if(i>=getPieceSize()) //Array needs to be allocated
   expandTokenUpdate(i);
   
   if(getTokenUpdate()[i]==NOT_KNOWN) //Array needs to be updated
   {  
   long result = i+nxtList[getID()];
   long val = result>=0?result:OUT_OF_BOUNDS;
   setTokenUpdateAtIndex(i,val);
   }
   
   return getTokenUpdate()[i];
   */
}

// --------- end of class derived_relation_node implmentation ---------
#endif


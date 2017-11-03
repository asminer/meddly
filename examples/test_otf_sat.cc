
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
 
  Model: A simple petri net
 
  Places: A, B
  Transitions: T_ab
 
  T_ab: A--, B++
  Initial state: A = N, B = 0
*/

#include "meddly.h"
#include "meddly_expert.h"

using namespace MEDDLY;

class pn;
class enabling_subevent;
class firing_subevent;
int usage(const char*);

class pn {
  public:
    pn(int nTokens);
    ~pn();
    dd_edge& getReachableStateSet();
    void clearRechableStateSet();

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

    satotf_opname::otf_relation* otf_rel;
    satotf_opname::subevent* A_firing;
    satotf_opname::subevent* B_firing;
    satotf_opname::subevent* A_enabling;
    satotf_opname::event* T_AB;
    specialized_operation* otf_sat_op;

    int n_tokens;
    domain* dom;
    expert_forest* mdd;
    expert_forest* mxd;
    dd_edge* initial_state;
    dd_edge* reachable_states;

    const int n_vars = 3;
    const int A_level = 3;
    const int B_level = 2;
    const int C_level = 1;

    // index to value
    std::vector<std::vector<int>> index_to_value;
    std::vector<std::vector<int>> value_to_index;

    friend enabling_subevent;
    friend firing_subevent;
};

class enabling_subevent : public satotf_opname::subevent {
  public:
    enabling_subevent(pn* model, forest* mxd,
        int* subevent_vars, int n_subevent_vars,
        int* event_vars, int n_event_vars);
    ~enabling_subevent();
    virtual void confirm(satotf_opname::otf_relation &rel, int v, int index);
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

class firing_subevent : public satotf_opname::subevent {
  public:
    firing_subevent(pn* model, forest* mxd,
        int* subevent_vars, int n_subevent_vars,
        int* event_vars, int n_event_vars,
        int offset);
    ~firing_subevent();
    virtual void confirm(satotf_opname::otf_relation &rel, int v, int index);
    void initializeMinterm();
  protected:
    pn* model;
    forest* mxd;
    int subevent_var;
    int* event_vars;
    int n_event_vars;
    int offset;
    int* unp_minterm;
    int* p_minterm;
    int minterm_size;
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

  otf_rel = 0;
  A_firing = 0;
  B_firing = 0;
  A_enabling = 0;
  T_AB = 0;
  otf_sat_op = 0;

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
  clearRechableStateSet();
}

pn::~pn() {
  // must:
  // delete the otf relation
  // delete the events
  // delete the domain

  // delete the otf relation
  // delete otf_rel;

  // delete the events
  delete T_AB;

  // delete the domain
  MEDDLY::cleanup();
}

void pn::initializeMeddly() {
  MEDDLY::initialize();
}

void pn::buildDomain() {
  int bounds[] = {-1, -1, -1};
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
  mxd = static_cast<expert_forest*>(dom->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL));
  assert(mxd);
}

void pn::buildTransition() {
  // Only one transition in the system: T_ab: A--, B++
  // build A-- subevent
  // build B++ subevent
  // build array of variables this transition depends on

  if (A_firing) delete A_firing;
  if (A_enabling) delete A_enabling;
  if (B_firing) delete B_firing;
  if (T_AB) delete T_AB;

  int dep_vars_A[] = {A_level};
  int dep_vars_B[] = {B_level};
  int dep_vars_TAB[] = {A_level, B_level};

  A_firing = new firing_subevent(this, mxd, (int*)(&dep_vars_A[0]), 1, (int*)(&dep_vars_TAB[0]), 2, -1);
  B_firing = new firing_subevent(this, mxd, (int*)(&dep_vars_B[0]), 1, (int*)(&dep_vars_TAB[0]), 2, 1);
  A_enabling = new enabling_subevent(
      this,                       // pn model
      mxd,                        // forest
      (int*)(&dep_vars_A[0]),     // dependant variable
      1,                          // # dependant variables
      (int*)(&dep_vars_TAB[0]),   // event dependant variables
      2                           // # event dependant variables
      );

  satotf_opname::subevent* subevents[] = {A_firing, B_firing, A_enabling};
  T_AB = new satotf_opname::event((satotf_opname::subevent**)(&subevents[0]), 3);
  assert(T_AB);
}

void pn::buildInitialState() {
  initial_state = new dd_edge(mdd);
  int init_state[] = {0, 0, 0, 0};
  int* init[] = {init_state};
  mdd->createEdge((int**)(&init[0]), 1, *initial_state);
  assert(initial_state);
  int index = -1;
  assert(addValue(A_level, n_tokens, index));
  assert(addValue(B_level, 0, index));
  assert(addValue(C_level, 0, index));
}

void pn::buildOtfSaturationOp() {
  if (otf_sat_op == 0) {
    if (otf_rel == 0) {
      // build otf saturation relation
      satotf_opname::event* events[] = {T_AB};
      otf_rel = new satotf_opname::otf_relation(mdd, mxd, mdd, &events[0], 1);
    }
    assert(otf_rel);
    otf_sat_op = SATURATION_OTF_FORWARD->buildOperation(otf_rel);
    assert(otf_sat_op);
  }
}

void pn::doOtfSaturation() {
  if (otf_sat_op == 0) buildOtfSaturationOp();
  assert(otf_sat_op);

  if (reachable_states) delete reachable_states;
  reachable_states = new dd_edge(mdd);
  assert(reachable_states);

  otf_rel->confirm(*initial_state);

  otf_sat_op->compute(*initial_state, *reachable_states);
  assert(reachable_states);
}

dd_edge& pn::getReachableStateSet() {
  if (reachable_states == 0) doOtfSaturation();
  assert(reachable_states);
  return *reachable_states;
}

void pn::clearRechableStateSet() {
  if (reachable_states) delete reachable_states;
  reachable_states = 0;
}

bool pn::indexOf(int vh, int value, int &index) {
  if (value_to_index.size() <= vh) return false;

  if (value_to_index[vh].size() <= value) return false;
  index = value_to_index[vh][value];
  return true;
}

bool pn::valueOf(int vh, int index, int &value) {
  if (index_to_value.size() <= vh) return false;
  if (index_to_value[vh].size() <= index) return false;
  value = index_to_value[vh][index];
  return true;
}

bool pn::addValue(int vh, int value, int &index) {
  if (indexOf(vh, value, index) && index >= 0) return false;

  if (value_to_index.size() <= vh) value_to_index.resize(vh+1);
  if (index_to_value.size() <= vh) index_to_value.resize(vh+1);
  if (value_to_index[vh].size() <= value) {
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
: satotf_opname::subevent(mxd, se_vars, n_se_vars, false)
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

void enabling_subevent::confirm(satotf_opname::otf_relation &rel, int v, int index) {
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
    int offset)
: satotf_opname::subevent(mxd, se_vars, n_se_vars, true)
{
  this->model = model;
  this->mxd = mxd;
  assert(n_se_vars == 1);
  subevent_var = se_vars[0];
  n_event_vars = n_e_vars;
  event_vars = new int[n_event_vars];
  for (int i = 0; i < n_event_vars; i++) event_vars[i] = e_vars[i];
  this->offset = offset;
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
  for (int i = 1; i < minterm_size; i++) {
    unp_minterm[i] = DONT_CARE;
    p_minterm[i] = DONT_CHANGE;
  }
  for (int i = 0; i < n_event_vars; i++) {
    p_minterm[event_vars[i]] = DONT_CARE;
  }
}

void firing_subevent::confirm(satotf_opname::otf_relation &rel, int v, int index) {
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

  int next_value = value + offset;
  if (next_value < 0) return;

  int next_value_index = -1;
#ifdef DEVELOPMENT_CODE
  assert(model->addValue(v, next_value, next_value_index));
#else
  model->addValue(v, next_value, next_value_index);
#endif
  // add enabling minterm
  unp_minterm[v] = index;
  p_minterm[v] = next_value_index;
  addMinterm(unp_minterm, p_minterm);
  std::cout << "\t\tk: " << v
  << ", index: " << index
  << ", value: " << value
  << ", next_value: " << next_value
  << ", next_value_index: " << next_value_index
  << "\n";
}

// ----------- end of class firing_subevent implementation ----------


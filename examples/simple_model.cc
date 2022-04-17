
/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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

#include "../src/meddly.h"
#include "../src/meddly_expert.h"
#include "simple_model.h"

#define VERBOSE
#define NOT_KNOWN     -2
#define OUT_OF_BOUNDS -1



// #define DEBUG_GENERATE
// #define DEBUG_EVENTS
// #define DEBUG_EVENTS_LONG

// #define SHOW_EVENT_HANDLES

// #define SAME_FOREST_OPERATIONS

inline int MAX(int a, int b) {
  return (a>b) ? a : b;
}

/*
    Common core for building the next state function.

      @param  events    Array of dimension \a nEvents.
                        Each entry is a string where
                        str[j] tells how this event affects variable j:
                        '+' means the event increases the variable by one
                            (unless this would violate the bound)
                        '-' means the event decreases the variable by one
                            (unless the value is already zero)
                        '.' means no change.
      @param  nEvents   Size of events array.
      @param  mxd       The forest to use; should be a boolean one for relations.

      @param  pnsf      If not null, store the relation for each event here.

      @param  mono      If not null, accumulate the monolithic relation here.

      @param  verb      Verbosity level.
*/

void buildNextStateFunction(const char* const* events, int nEvents,
  MEDDLY::forest* mxd, 
  MEDDLY::satpregen_opname::pregen_relation* pnsf,
  MEDDLY::dd_edge* mono, int verb)
{
  using namespace MEDDLY;
  if (verb) fprintf(stderr, "Building next-state function\n");

  // set up auxiliary mtmxd forest and edges
  domain* d = mxd->useDomain();
  forest* mtmxd = d->createForest(
    true, forest::INTEGER, forest::MULTI_TERMINAL
  );
  int nVars = d->getNumVariables();
  int maxBound = d->getVariableBound(1, false);
  for (int i=2; i<=nVars; i++) {
    maxBound = MAX(maxBound, d->getVariableBound(i, false));
  }
  maxBound++;
  long* temp = new long[maxBound];
  int* minterm = new int[nVars+1];
  int* mtprime = new int[nVars+1];
  dd_edge** varP  = new dd_edge*[nVars+1];
  varP[0] = 0;
  dd_edge** inc   = new dd_edge*[nVars+1];
  inc[0] = 0;
  dd_edge** dec   = new dd_edge*[nVars+1];
  dec[0] = 0;

  //  Create edge for each variable xi'
  for (int i=1; i<=nVars; i++) {
    varP[i] = new dd_edge(mtmxd);
    mtmxd->createEdgeForVar(i, true, varP[i][0]);
  }

  // Create edge for each function xi+1
  for (int i=0; i<maxBound; i++) temp[i] = i+1;
  for (int i=1; i<=nVars; i++) {
    inc[i] = new dd_edge(mtmxd);
    mtmxd->createEdgeForVar(i, false, temp, inc[i][0]);
  }

  // Create edge for each function xi-1
  for (int i=0; i<maxBound; i++) temp[i] = i-1;
  for (int i=1; i<=nVars; i++) {
    dec[i] = new dd_edge(mtmxd);
    mtmxd->createEdgeForVar(i, false, temp, dec[i][0]);
  }

  //
  // Initialize accumulators
  //
  if (mono) mxd->createEdge(false, *mono);

  for (int e=0; e<nEvents; e++) {
    const char* ev = events[e];
    if (2==verb) fprintf(stderr, "%5d", e);
    if (verb>2) fprintf(stderr, "Event %5d", e);

    dd_edge nsf_ev(mxd);
    dd_edge term(mxd);

    //
    // build mask for this event
    //
    for (int i=1; i<=nVars; i++) {
      if ('.' == ev[i]) {
        minterm[i] = DONT_CARE;
        mtprime[i] = DONT_CHANGE;
      } else {
        minterm[i] = DONT_CARE;
        mtprime[i] = DONT_CARE;
      }
    }
    mxd->createEdge(&minterm, &mtprime, 1, nsf_ev);
#ifdef DEBUG_EVENTS
    printf("Initial nsf for event %d\n", e);
    nsf_ev.show(stdout, 2);
#endif
    if (verb>2) fprintf(stderr, " : ");
    
    //
    // 'and' with the "do care" levels
    //
    for (int i=1; i<=nVars; i++) {
#ifdef SAME_FOREST_OPERATIONS
      dd_edge docare(mtmxd);
#endif

      if ('.' == ev[i]) {
        if (verb>3) fputc('.', stderr);
        continue;
      } else {
        if (verb>2) fputc(ev[i], stderr);
      }
      switch (ev[i]) {
        case '+':   
#ifdef SAME_FOREST_OPERATIONS
                    apply(EQUAL, varP[i][0], inc[i][0], docare);
#else
                    apply(EQUAL, varP[i][0], inc[i][0], term);
#endif
                    break;

        case '-':   
#ifdef SAME_FOREST_OPERATIONS
                    apply(EQUAL, varP[i][0], dec[i][0], docare);
#else
                    apply(EQUAL, varP[i][0], dec[i][0], term);
#endif
                    break;

        default:    throw 1;
      } // switch
#ifdef SAME_FOREST_OPERATIONS
      apply(COPY, docare, term);
#endif
#ifdef DEBUG_EVENTS
      printf("Term for event %d, level %d\n", e, i);
      term.show(stdout, 2);
#endif
      nsf_ev *= term;
    } // for i

#ifdef DEBUG_EVENTS
    printf("Complete nsf for event %d:\n", e);
    nsf_ev.show(stdout, 2);
#endif
    if (verb>2) fputc(' ', stderr);
#ifdef SHOW_EVENT_HANDLES
    fprintf(stderr, "%d", nsf_ev.getNode());
#endif

    if (pnsf) {
      //
      //  Add event to relation
      //
      pnsf->addToRelation(nsf_ev);
    }

    if (mono) {
      //
      //  union with overall
      //
      *mono += nsf_ev;
#ifdef DEBUG_EVENTS_LONG
      printf("Complete after adding event %d:\n", e);
      mono->show(stdout, 2);
#endif
    }

    if (verb>2) fputc('\n', stderr);
  } // for e
  if (verb==2) fputc('\n', stderr);

  // cleanup
  delete[] mtprime;
  delete[] minterm;
  delete[] temp;
  for (int i=1; i<=nVars; i++) {
    delete varP[i];
    delete inc[i];
    delete dec[i];
  }
  delete[] varP;
  delete[] inc;
  delete[] dec;
  destroyForest(mtmxd);

#ifdef DEBUG_EVENTS
  printf("Complete NSF:\n");
  mono->show(stdout, 2); 
#endif
}


void buildNextStateFunction(const char* const* events, int nEvents,
  MEDDLY::forest* mxd, MEDDLY::dd_edge &nsf, int verb)
{
  buildNextStateFunction(events, nEvents, mxd, 0, &nsf, verb);
}

void buildNextStateFunction(const char* const* events, int nEvents,
  MEDDLY::satpregen_opname::pregen_relation* pnsf, int verb)
{
  if (0==pnsf) return;
  buildNextStateFunction(events, nEvents, pnsf->getRelForest(), pnsf, 0, verb);
  pnsf->finalize();

#ifdef SHOW_EVENT_HANDLES
  using namespace MEDDLY;
  // check what we got
  for (int k=16; k; k--) {
    int len = pnsf->lengthForLevel(k);
    if (0==len) continue;
    printf("Events at level %d:\n\t", k);
    node_handle* List = pnsf->arrayForLevel(k);
    for (int i=0; i<len; i++)
      printf("%d ", List[i]);
    printf("\n");
  }
#endif
}




//
//  Explicit RS construction
//

bool fireEvent(const char* event, const int* current, int* next, int nVars)
{
  for (int i=nVars; i; i--) {
    if ('.' == event[i]) {
      next[i] = current[i];
      continue;
    }
    if ('-' == event[i]) {
      next[i] = current[i] - 1;
      if (next[i] < 0) return false;
      continue;
    }
    if ('+' == event[i]) {
      next[i] = current[i] + 1;
      // TBD ... check for overflow
      continue;
    }
    throw 1;  // bad event string
  }
  return true;
}

void explicitReachset(const char* const* events, int nEvents, 
  MEDDLY::forest* f, MEDDLY::dd_edge &expl, MEDDLY::dd_edge &RS, int batchsize)
{
  int b;
  int nVars = f->getDomain()->getNumVariables();
  if (batchsize < 1) batchsize = 256;

  // initialize batch memory
  int** minterms = new int*[batchsize];
  for (b=0; b<batchsize; b++) {
    minterms[b] = new int[1+nVars];
  }
  
  // unexplored states
  MEDDLY::dd_edge unexplored(f);
  // batch of states
  MEDDLY::dd_edge batch(f);
  b = 0; 
  // exploration loop.
  MEDDLY::enumerator I(expl);
  for (;;) {
    unexplored.clear();
    I.start(expl);
    if (!I) break;    // nothing left to explore, bail out
    // explore everything in expl
    for (; I; ++I) {
      const int* curr = I.getAssignments();
#ifdef DEBUG_GENERATE
      printf("Exploring state: (%d", curr[1]);
      for (int n=2; n<=nVars; n++) printf(", %d", curr[n]);
      printf(")\n");
#endif
      // what's enabled?
      for (int e=0; e<nEvents; e++) {
        if (!fireEvent(events[e], curr, minterms[b], nVars)) continue;
#ifdef DEBUG_GENERATE
        printf("  -- (event %d) --> (%d", e, minterms[b][1]);
        for (int n=2; n<=nVars; n++) printf(", %d", minterms[b][n]);
        printf(")\n");
#endif
        bool seen;
        f->evaluate(RS, minterms[b], seen);
        if (seen) continue;     // already known in RS
        f->evaluate(unexplored, minterms[b], seen);
        if (seen) continue;     // already in unexplored list
        b++;
        if (b>=batchsize) {
          // Buffer is full; flush it
          f->createEdge(minterms, b, batch);
          unexplored += batch;
          RS += batch;
          b = 0;
        }
      }
    }
    // expl is empty.
    // Flush the buffer
    if (b) {
      f->createEdge(minterms, b, batch);
      unexplored += batch;
      b = 0;
    }
    RS += unexplored;
    expl = unexplored;
  }

  // cleanup batch memory
  for (b=0; b<batchsize; b++) {
    delete[] minterms[b];
  }
  delete[] minterms;
}


/*************************************************************/
int* nxtList;
class derRelNode : public MEDDLY::relation_node
{
  int en;
  int fr;
  int h;
public:
  derRelNode(MEDDLY::forest* mxdF, int level, int dwn, int e, int f, int inh): relation_node(141010, mxdF, level, dwn, e, f, inh)
  {
    en = e;
    fr = f;
    h = inh;
  }
  
  long nextOf(long i) override
  {
  
    if(i>=getPieceSize()) { //Array needs to be allocated
       expandTokenUpdate(i);
    }
  
    if(getTokenUpdate()[i]==NOT_KNOWN) //Array needs to be updated
    {	
    	long result = ((i>=en) && (h==-1?true:i<h))? i+fr : OUT_OF_BOUNDS;
    	long val = result;
     	setTokenUpdateAtIndex(i,val);
    }
   return getTokenUpdate()[i];
  }
};


void buildImplicitRelation(const int* const* events, int nEvents,int nPlaces, int bounds, MEDDLY::forest* mddF, MEDDLY::forest* mxdF, MEDDLY::satimpl_opname::implicit_relation* T)
{
  
  unsigned node_count = 0;
  int* tops_of_events = (int*)malloc(size_t(nEvents)*sizeof(int));
  
  for(int e = 0;e < nEvents; e++)
    {
    bool done = false;
    for( int p = 1; p <= nPlaces; p++)
      {
      if(events[e][p]!='0') node_count +=1;
      if((events[e][nPlaces-p+1]!=0)&&(!done))
        {
        tops_of_events[e] = nPlaces-p+1;
        done = true;
        }
      }
    }
  
  
  derRelNode** rNode = (derRelNode**)malloc(node_count*sizeof(derRelNode*));
  
  // Add/Subtract Tokens
  nxtList = (int*)malloc((node_count+2)*sizeof(int));
  nxtList[0] = 0;
  nxtList[1] = 0;
  
  int rctr = 0;
  for(int e = 0;e < nEvents; e++)
    {
    unsigned long sign = 0;
    int previous_node_handle = 1;
    for( int p = 1; p <= nPlaces; p++)
      {
       sign = events[e][p]>=0?(sign*10)+events[e][p]:(sign*100)+events[e][p];
        if(events[e][p]!=0)
        {
          rNode[rctr] = new derRelNode(mxdF,p,previous_node_handle,events[e][p]<0?-events[e][p]:0,events[e][p],-1);
          previous_node_handle = T->registerNode((tops_of_events[e]==p),rNode[rctr]);
          
          rctr++;
          nxtList[previous_node_handle] = events[e][p];
        }
      }
    }
}


MEDDLY::sathyb_opname::event** buildHybridRelation(const int* const* events, int nEvents,int nPlaces, int bounds, MEDDLY::forest* mddF, MEDDLY::forest* mxdF) {

  unsigned node_count = 0;
  int* tops_of_events = (int*)malloc(size_t(nEvents)*sizeof(int));
  
  for(int e = 0;e < nEvents; e++)
    {
    bool done = false;
    for( int p = 1; p <= nPlaces; p++)
      {
      if(events[e][p]!='0') node_count +=1;
      if((events[e][nPlaces-p+1]!=0)&&(!done))
        {
        tops_of_events[e] = nPlaces-p+1;
        done = true;
        }
      }
    }
  
  derRelNode** rNode = (derRelNode**)malloc(node_count*sizeof(derRelNode*));
  MEDDLY::sathyb_opname::event** T = (MEDDLY::sathyb_opname::event**)malloc(nEvents*sizeof(MEDDLY::sathyb_opname::event*));
  
  // Add/Subtract Tokens
  nxtList = (int*)malloc((node_count+2)*sizeof(int));
  nxtList[0] = 0;
  nxtList[1] = 0;
  
  int rctr = 0;
  for(int e = 0;e < nEvents; e++)
    {
    unsigned long sign = 0;
    int previous_node_handle = 1;
    int e_rn = 0;
    for( int p = 1; p <= nPlaces; p++)
      {
       sign = events[e][p]>=0?(sign*10)+events[e][p]:(sign*100)+events[e][p];
        if(events[e][p]!=0)
        {
          e_rn++;
          rNode[rctr] = new derRelNode(mxdF,p,-1,events[e][p]<0?-events[e][p]:0,events[e][p],-1);
          rctr++;
          nxtList[previous_node_handle] = events[e][p];
        }
      }
       T[e] = new MEDDLY::sathyb_opname::event(NULL, 0, (MEDDLY::relation_node**)(&rNode[rctr-e_rn]), e_rn);
    }

    return T;

}

/*************************************************************/

#if 1
class param_rel_node : public MEDDLY::gen_relation_node
{
  std::vector<int> dep_variable;
  int noofdep;
  std::vector<int> coeffs;
  long* constant_wgts;
  std::vector<int> vars_in;
  std::vector<int> vars_free;
  std::vector<int> vars_out;
  long self_coeff;
  
  
public:
  param_rel_node(unsigned long signature, MEDDLY::forest* f, int level, MEDDLY::node_handle down, 
    long* const_wgts, int* dep_vars, int noof_dep_vars, long* coefficients, long own_coeff) : gen_relation_node(signature, f, level, down)
  {
    noofdep = noof_dep_vars;
    if(const_wgts != NULL)
    {
      constant_wgts = (long*)malloc(sizeof(long)*3);
      constant_wgts[0] = const_wgts[0];
      constant_wgts[1] = const_wgts[1];
      constant_wgts[2] = const_wgts[2];
    }
    for(int i = 0; i<noofdep; i++)
    {
      dep_variable.push_back(dep_vars[i]); // what variables;
      coeffs.push_back(coefficients[i]);  // 2*p1 + (-4*p2) - (3*p3)
    }
    self_coeff = own_coeff;
    
  }
  
  void setInputsFreeOut(int* ip_var, int isz, int* free_var, int fsz, int* out_var, int osz) {
    printf("\n I am here inside setInputsFreeOut \n");
    printf("\n I am here inside setInputsFreeOut \n");
    for(int i = 0; i < isz; i++)
    {
      vars_in.push_back(ip_var[i]);
    }
    
    for(int i = 0; i < fsz; i++)
    {
      vars_free.push_back(free_var[i]);
    }
    
    for(int i = 0; i < osz; i++)
    {
      vars_out.push_back(out_var[i]);
    }
  }
  
  
  std::set<std::vector<long>> cartesianProductTwo(std::set<std::vector<long>> set_a, int var, MEDDLY::expert_forest* resF) {
    std::set<std::vector<long>> result;
    for (auto it = set_a.begin(); it!= set_a.end(); it++)
    {
      MEDDLY::expert_domain* dom = static_cast<MEDDLY::expert_domain*>(resF->useDomain());
      int set_b = dom->getVariableBound(var);
      for(int i = 0; i <set_b; i++){
        //*it = vector = element of a set
          std::vector<long> new_item = *it;
          new_item.push_back(i);
          result.insert(new_item);
      }
    }
    return result;
  }
  
  
  long* buildFreeValues(MEDDLY::forest* resF, long* inp_values, long i, int& jsize) override {
    std::set<std::vector<long>> resultProd;
     std::vector<long> empty_vector;
    for(int i = 0; i<vars_free.size(); i++)
     {
       int var = vars_free[i];
       MEDDLY::expert_domain* dom = static_cast<MEDDLY::expert_domain*>(resF->useDomain());
       int dm = dom->getVariableBound(var);
       if(i == 0){
         for(int j = 0; j < dm ; j++)
         {
             std::vector<long> new_item;
             new_item.push_back(j);
             resultProd.insert(new_item);
         }
       } else {
         resultProd = cartesianProductTwo(resultProd, var, resF);
       }
     }
     if(resultProd.empty()) resultProd.insert(empty_vector);
     return resultProd;
  }
  
  long delta(long* v_in, long* v_free, long i) override {
    int expr = 0;
    
    int d_k = 0;
    int i_k = 0;
    int f_k = 0;
    
    if(constant_wgts!=NULL)
    {
      expr += constant_wgts[2] + constant_wgts[0];
    }
    /*
    printf("\n here to get next of var%d:%d", this->getLevel(),i);
    if(dep_variable.size()>0)printf("\n How many deps of this %d level? %d %d", this->getLevel(),dep_variable.size(),dep_variable[0]);
    if(value_free.size()>0) printf("\n here to get next of freear?:%d (myvar:%d)", value_free[0], vars_free[0]);
    if(value_in.size()>0) printf("\n here to get next of inpvar?:%d (myvar:%d)", value_in[0], vars_in[0]);
    if(value_in.size()>1) printf("\n here to get next of inpvar?:%d (myvar:%d)", value_in[1], vars_in[1]);
    if(vars_out.size()>1) printf("\n here to get outpvar is %d",vars_out[0]);
    
    for(int x=0;x<dep_variable.size();x++)
    {
      printf("\n I depend on variable %d", dep_variable[x]);
    }
     */
    
    while((d_k<dep_variable.size()) && (i_k<vars_in.size()) && (f_k<vars_free.size()))
    {
      if(dep_variable[d_k] < vars_free[f_k]) // dep must exist in vars_in
      {
        while(dep_variable[d_k] != vars_in[i_k])
        {
          i_k++;
        }
        expr += coeffs[d_k]*value_in[i_k];
        d_k++;
        i_k++;
      }else if(dep_variable[d_k] == vars_free[f_k])
      {
        expr += coeffs[d_k]*value_free[f_k];
        d_k++;
        f_k++;
      }
    }
    
    if(d_k<dep_variable.size())
    {
      if(f_k<vars_free.size()) //remaining dependents must be free because those aren't in input
      {
        while(d_k<dep_variable.size()){
        //printf("\n dep_variable[d_k=%d] = %d & free_vars[f_k=%d] = %d", d_k,dep_variable[d_k], f_k, vars_free[f_k]);
        expr += coeffs[d_k]*value_free[f_k];
        d_k++;
        f_k++;
        }
      } else if (i_k<vars_in.size()) // some of the inputs cover the remaining dependents
      {
        while(d_k<dep_variable.size()){
          if(dep_variable[d_k] == vars_in[i_k])
          {
           expr += coeffs[d_k]*value_in[i_k];
           d_k++;
           i_k++;
          }else i_k++;
        }
      }
    }
    
    int j = i*(1+self_coeff) + expr;
     
     return j>=0?j:-1;
  }
  
  long* omega(long* v_in, long i, long j) override {
  
     std::vector<long> value_out;
    
     int o_k = 0;
     int i_k = 0;
     int f_k = 0;
    
     while((o_k<vars_out.size()) && (i_k<vars_in.size()) && (f_k<vars_free.size()))
     {
       if(vars_out[o_k] < vars_free[f_k]) // op must exist in vars_in
       {
         while(vars_out[o_k] != vars_in[i_k])
         {
           i_k++;
         }
         value_out.push_back(value_in[i_k]);
         o_k++;
         i_k++;
       }else if(vars_out[o_k] == vars_free[f_k])
       {
         value_out.push_back(value_free[f_k]);
         o_k++;
         f_k++;
       }else if(vars_out[o_k] == this->getLevel())
       {
         value_out.push_back(i);
         o_k++;
       }
     }
     
     if(o_k<vars_out.size())
     {
       if(f_k<vars_free.size()) //remaining dependents must be free because those aren't in input
       {
         while(o_k<vars_out.size()){
         if(vars_out[o_k] == vars_free[f_k])
         {
           value_out.push_back(value_free[f_k]);
           o_k++;
           f_k++;
         }else if(vars_out[o_k] == this->getLevel())
         {
           value_out.push_back(i);
           o_k++;
         }
         }
       } else if (i_k<vars_in.size()) // some of the inputs cover the remaining dependents
       {
         while(o_k<vars_out.size()){
           if(vars_out[o_k] == vars_in[i_k])
           {
            value_out.push_back(value_in[i_k]);
            o_k++;
            i_k++;
           }else if(vars_out[o_k] == this->getLevel())
           {
             value_out.push_back(i);
             o_k++;
           }else i_k++;
         }
       }
     }
    
    return value_out;
  }
};

/*
 unsigned long signature, MEDDLY::forest* f, int level, node_handle down, long* constant_wgts, long* dep_vars, int noof_dep_vars, long* coefficients
 */
//std::vector< std::map<int, std::map<int, int> > > events
// event = std::map<int, std::map<int, int>> arcs
// arc = place : map of place with coeff

void buildGenericImplicitRelation(std::vector<std::map<int,std::map<int, int>>> events, int nEvents, MEDDLY::forest* mxdF, MEDDLY::satmdimpl_opname::md_implicit_relation* T)
{
  
  for(int e = 0; e<nEvents; e++)
    {
      std::map<int, std::map<int, int>> this_event = events[e];
      int rctr = 0;
      param_rel_node** rNode = (param_rel_node**)malloc(this_event.size()*sizeof(param_rel_node*));
      for(auto arc_it = this_event.begin(); arc_it != this_event.end(); arc_it++, rctr++ )
      {
        int aff_level = arc_it->first;
        long* constant_wgts = (long*)malloc(3*sizeof(long));
        int* dep_vars = (int*)malloc(arc_it->second.size()*sizeof(int));
        long* coefficients = (long*)malloc(arc_it->second.size()*sizeof(long));
        int dep_var_ct = 0; // keep track of dep_vars
        long own_coeff = 0;
        std::string node_sign = "p"+std::to_string(aff_level);
        for(auto wgt_it = arc_it->second.begin(); wgt_it != arc_it->second.end(); wgt_it++)
        {
          if(wgt_it->first == 0) // variable = 0 meaning const wgts
          {
            if(wgt_it->second > 0)
            {
              constant_wgts[2] = (long)wgt_it->second;
              constant_wgts[0] = 0;
              node_sign +=std::to_string(constant_wgts[2]);
            }
            else
            {
              constant_wgts[0] = (long)wgt_it->second;
              constant_wgts[2] = 0;
              node_sign += "-"+std::to_string(constant_wgts[0]);
            }
          }else{
            constant_wgts = NULL;
            if(aff_level != wgt_it->first)
            {
              dep_vars[dep_var_ct] = wgt_it->first;
              coefficients[dep_var_ct] = wgt_it->second;
              node_sign += std::to_string(coefficients[dep_var_ct])+"* p"+std::to_string(dep_vars[dep_var_ct]);
              dep_var_ct++;
            } else {
              own_coeff = wgt_it->second;
              node_sign += std::to_string(own_coeff)+"* p"+std::to_string(aff_level);
            }
          }
        } // effect of this_event on aff_level is all known.
        // build the node for it.
        rNode[rctr] = new param_rel_node(node_sign, mxdF, aff_level, -1, constant_wgts, dep_vars, dep_var_ct, coefficients, own_coeff);
      }
      // Use all the nodes of this event to register the event
      T->registerEventNodes((MEDDLY::gen_relation_node**)rNode, this_event.size());
      for(int k=0; k <this_event.size(); k++)
      {
        rNode[k]->setInputsFreeOut(rNode[k]->getInputVariables(), rNode[k]->getCountInput(), rNode[k]->getFreeVariables(), rNode[k]->getCountFree(), rNode[k]->getOutputVariables(), rNode[k]->getCountOutput());
      }
      
    }
}
#endif
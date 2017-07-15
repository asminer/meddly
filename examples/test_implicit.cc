// $Id$

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

#include <cstdlib>
#include <string.h>
#include <fstream>

#define _MEDDLY_WITHOUT_IOSTREAM_

#include "meddly.h"
#include "meddly_expert.h"
#include "simple_model.h"
#include "timer.h"
#include "loggers.h"

// #define DUMP_NSF
// #define DUMP_REACHABLE


char** sample_model;
int p1_position;
int PLACES = 3;
int TRANS = 1;
int N = -1;


using namespace MEDDLY;

FILE_output meddlyout(stdout);

void buildModel(const char* order)
{
  const char* sample_model1[] = {
    "X-++"  // T1
  };
  
  p1_position = 1;
  for (int i=0; i<PLACES; i++)
    {
    if (order[i]=='1')
      p1_position = i+1;
    }
  
  sample_model = (char**) malloc(TRANS * sizeof(char*));
  
  for(int i=0;i<TRANS;i++)
    {
    sample_model[i] = (char*) malloc((PLACES+2) * sizeof(char));
    for(int j=1;j<(PLACES+1);j++)
      {
      sample_model[i][j]=sample_model1[i][order[j-1]-'0'];
      std::cout<< order[j-1]-'0' << "\n";
      std::cout<< sample_model1[i][order[j-1]-'0'] << "\n";
      }
    sample_model[i][PLACES+1] = '\0';
    std::cout<<sample_model[i];
    }
  
}


int usage(const char* who)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) name = ptr+1;
  }
  printf("\nUsage: %s [options] nnnn \n\n", name);
  printf("\tnnnn: number of initial tokens\n");
  printf("\t-impl: use implicit nodes\n");
  printf("\t-bfs: use traditional iterations\n");
  printf("\t-dfs: use saturation\n");
  printf("\t-esat: use saturation by events\n");
  printf("\t-ksat: use saturation by levels\n");
  printf("\t-msat: use monolithic saturation (default)\n\n");
  printf("\t-exp: use explicit (very slow)\n\n");
  printf("\t--batch b: specify explicit batch size\n\n");
  printf("\t -l lfile: Write logging information to specified file\n\n");
  printf("\t -O order: Write the order of variables\n\n");
  
  return 1;
}

void printStats(const char* who, const forest* f)
{
  printf("%s stats:\n", who);
  const expert_forest* ef = (expert_forest*) f;
  ef->reportStats(meddlyout, "\t",
                  expert_forest::HUMAN_READABLE_MEMORY  |
                  expert_forest::BASIC_STATS | expert_forest::EXTRA_STATS |
                  expert_forest::STORAGE_STATS | expert_forest::HOLE_MANAGER_STATS
                  );
}


int main(int argc, const char** argv)
{
  char method = 'i';
  int batchsize = 256;
  const char* lfile = 0;
  
  for (int i=1; i<argc; i++) {
    
    /*if (strcmp("-O", argv[i])==0) {
     if(i+1 < argc)
     {
     buildModel(argv[i+1]);
     i++;
     }
     continue;
     }
     if (strcmp("-bfs", argv[i])==0) {
     method = 'b';
     continue;
     }
     if (strcmp("-dfs", argv[i])==0) {
     method = 'm';
     continue;
     }
     if (strcmp("-esat", argv[i])==0) {
     method = 's';
     continue;
     }
     if (strcmp("-ksat", argv[i])==0) {
     method = 'k';
     continue;
     }
     if (strcmp("-msat", argv[i])==0) {
     method = 'm';
     continue;
     }
     if (strcmp("-exp", argv[i])==0) {
     method = 'e';
     continue;
     }
     if (strcmp("-l", argv[i])==0) {
     lfile = argv[i+1];
     i++;
     continue;
     }
     if (strcmp("--batch", argv[i])==0) {
     i++;
     if (argv[i]) batchsize = atoi(argv[i]);
     continue;
     }*/
    N = atoi(argv[i]);
  }
  
  
  if (N<0) return usage(argv[0]);
  
  domain* d = 0;
  try {
    
    MEDDLY::initialize();
    
    timer start;
    
    printf("+----------------------------------------------------+\n");
    printf("|   Initializing sample_model  model with %-4d  |\n", N);
    printf("+----------------------------------------------------+\n");
    fflush(stdout);
    
    // Initialize domain
    int* sizes = new int[PLACES];
    for (int i=PLACES-1; i>=0; i--) sizes[i] = N+1;
    d = createDomainBottomUp(sizes, PLACES);
    delete[] sizes;
    forest::policies pr(true);
    pr.setPessimistic();
    forest::policies p(false);
    p.setPessimistic();
    
    // associate loggers
    /*std::ofstream log;
     forest::logger* LOG = 0;
     if (lfile) {
     log.open(lfile, std::ofstream::out);
     if (!log) {
     printf("Couldn't open %s for writing, no logging\n", lfile);
     } else {
     LOG = new simple_logger(log);
     LOG->recordNodeCounts();
     LOG->recordTimeStamps();
     char comment[80];
     snprintf(comment, 80, "Automatically generated by sample_model construction (N=%d)", N);
     LOG->addComment(comment);
     mdd->setLogger(LOG, "MDD");
     mxd->setLogger(LOG, "MxD");
     }
     }*/
    
    if('i' == method)
      {
      
      /*Single Event
       satimpl_opname::relation_node* T11 = new satimpl_opname::relation_node(1,true,1,3,3);
       satimpl_opname::relation_node* T12 = new satimpl_opname::relation_node(1,true,1,2,2);
       satimpl_opname::relation_node* T13 = new satimpl_opname::relation_node(1,false,1,1,1);*/
      
      
      /*2 Events*/
      satimpl_opname::relation_node* T11 = new satimpl_opname::relation_node(10,1,1);
      satimpl_opname::relation_node* T12 = new satimpl_opname::relation_node(11,2,2);
      satimpl_opname::relation_node* T22 = new satimpl_opname::relation_node(12,2,1);
      
      
      //satimpl_opname::relation_node* T23 = new satimpl_opname::relation_node(1,true,1,3,4);
      
      //CREATE FORESTS
      forest* inmdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL,p);
      forest* outmdd = inmdd;
      
      //ADD INITIAL STATE
      int* elementList[1];
      elementList[0] = new int[PLACES + 1];
      elementList[0][1] = N; elementList[0][2] = 0; elementList[0][3] = 0;
      dd_edge first(inmdd);
      dd_edge reachable(outmdd);
      inmdd->createEdge(elementList, 1, first);
      outmdd->createEdge(elementList, 1, reachable);
      
      
      //CREATE RELATION
      satimpl_opname::implicit_relation* T = new satimpl_opname::implicit_relation(inmdd,outmdd);
      
      /*T->registerNode(false,T13);
       T->registerNode(false,T12);
       T->registerNode(true,T11);*/
      T->registerNode(false,T11);
      T->registerNode(true,T12);
      T->registerNode(true,T22);
      //T->registerNode(true,T23);
      
      specialized_operation* sat = 0;
      
      printf("Building reachability set using saturation implicit relation");
      if (0==SATURATION_IMPL_FORWARD) {
        throw error(error::UNKNOWN_OPERATION);
      }
      sat = SATURATION_IMPL_FORWARD->buildOperation(T);
      
      if (0==sat) {
        throw error(error::INVALID_OPERATION);
      }
      sat->compute(first, reachable);
      
      start.note_time();
      printf("Done\n");
      printf("Reachability set construction took %.4e seconds\n",
             start.get_last_interval() / 1000000.0);
      fflush(stdout);
      
      //#ifdef DUMP_REACHABLE
      printf("Reachable states:\n");
      reachable.show(meddlyout, 2);
      //#endif
      
      printStats("MDD", outmdd);
      fflush(stdout);
      
      double c;
      apply(CARDINALITY, reachable, c);
      operation::showAllComputeTables(meddlyout, 3);
      
      printf("Approx. %g reachable states\n", c);
      
      
      }
    
    /*
     
     // Build initial state
     if (LOG) LOG->newPhase(mdd, "Building initial state");
     int* initial = new int[PLACES+1];
     for (int i=PLACES; i; i--) initial[i] = 0;
     initial[p1_position] = N;
     dd_edge init_state(mdd);
     mdd->createEdge(&initial, 1, init_state);
     delete[] initial;
     
     //
     // Build next-state function
     if (LOG) LOG->newPhase(mxd, "Building next-state function");
     dd_edge nsf(mxd);
     satpregen_opname::pregen_relation* ensf = 0;
     specialized_operation* sat = 0;
     
     if ('s' == method) {
     ensf = new satpregen_opname::pregen_relation(mdd, mxd, mdd, 4);
     }
     if ('k' == method) {
     ensf = new satpregen_opname::pregen_relation(mdd, mxd, mdd);
     }
     
     if ('e' != method) {
     
     if (ensf) {
     start.note_time();
     buildNextStateFunction(sample_model, TRANS, ensf, 4);
     start.note_time();
     } else {
     start.note_time();
     buildNextStateFunction(sample_model, TRANS, mxd, nsf, 4);
     start.note_time();
     #ifdef DUMP_NSF
     printf("Next-state function:\n");
     nsf.show(meddlyout, 2);
     #endif
     }
     printf("Next-state function construction took %.4e seconds\n",
     start.get_last_interval() / 1000000.0);
     printStats("MxD", mxd);
     }
     
     if (LOG) LOG->newPhase(mdd, "Building reachability set");
     dd_edge reachable(mdd);
     start.note_time();
     switch (method) {
     case 'b':
     printf("Building reachability set using traditional algorithm\n");
     fflush(stdout);
     apply(REACHABLE_STATES_BFS, init_state, nsf, reachable);
     break;
     
     case 'm':
     printf("Building reachability set using saturation, monolithic relation\n");
     fflush(stdout);
     apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
     break;
     
     case 'e':
     printf("Building reachability set using explicit search\n");
     printf("Using batch size: %d\n", batchsize);
     fflush(stdout);
     explicitReachset(sample_model, 4, mdd, init_state, reachable, batchsize);
     break;
     
     case 'k':
     case 's':
     printf("Building reachability set using saturation, relation");
     if ('k'==method)  printf(" by levels\n");
     else              printf(" by events\n");
     fflush(stdout);
     if (0==SATURATION_FORWARD) {
     throw error(error::UNKNOWN_OPERATION);
     }
     sat = SATURATION_FORWARD->buildOperation(ensf);
     if (0==sat) {
     throw error(error::INVALID_OPERATION);
     }
     sat->compute(init_state, reachable);
     break;
     
     default:
     printf("Error - unknown method\n");
     exit(2);
     };
     start.note_time();
     printf("Done\n");
     printf("Reachability set construction took %.4e seconds\n",
     start.get_last_interval() / 1000000.0);
     fflush(stdout);
     
     #ifdef DUMP_REACHABLE
     printf("Reachable states:\n");
     reachable.show(meddlyout, 2);
     #endif
     
     printStats("MDD", outmdd);
     fflush(stdout);
     
     double c;
     apply(CARDINALITY, outmdd, c);
     operation::showAllComputeTables(meddlyout, 3);
     
     printf("Approx. %g reachable states\n", c);
     */
    /*
     // cleanup
     if (LOG) {
     LOG->newPhase(mdd, "Cleanup");
     LOG->newPhase(mxd, "Cleanup");
     MEDDLY::destroyDomain(d);
     delete LOG;
     }
     MEDDLY::cleanup();*/
    return 0;
  }
  catch (MEDDLY::error e) {
    printf("Caught MEDDLY error: %s\n", e.getName());
    return 1;
  }
}


long MEDDLY::satimpl_opname::relation_node::nextOf(long i)
{
  
  //Assume variable domain is [0,2]
  
  
  if(this->ID == 2) //t1: Level1
    {
    if(1<=i) return i-1;
    else return i;
    }
  else if(this->ID == 3) //t1: Level2
    {
    if(1+i<N+1) return i+1;
    else return i;
    }
  else if (this->ID == 4) //t2: Level2
    {
    if(1<=i) return i-1;
    else return i;
    }
}



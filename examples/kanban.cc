
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
#include <gmp.h>

#define _MEDDLY_WITHOUT_IOSTREAM_

#include "../src/meddly.h"
#include "simple_model.h"
#include "../timing/timer.h"
#include "../src/log_simple.h"

// #define DUMP_NSF
// #define DUMP_REACHABLE

const char* kanban[] = {
  "X-+..............",  // Tin1 TA
  "X.-+.............",  // Tr1 TB
  "X.+-.............",  // Tb1 TC
  "X.-.+............",  // Tg1 TD
  "X.....-+.........",  // Tr2 TE
  "X.....+-.........",  // Tb2 TF
  "X.....-.+........",  // Tg2 TG
  "X+..--+..-+......",  // Ts1_23 TH
  "X.........-+.....",  // Tr3 TI
  "X.........+-.....",  // Tb3 TJ
  "X.........-.+....",  // Tg3 TK
  "X....+..-+..--+..",  // Ts23_4 TL
  "X.............-+.",  // Tr4 TM
  "X.............+-.",  // Tb4 TN
  "X............+..-",  // Tout4 TO
  "X.............-.+"   // Tg4 TP
};

using namespace MEDDLY;

FILE_output meddlyout(stdout);

int usage(const char* who)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) name = ptr+1;
  }
  printf("\nUsage: %s [options] nnnn \n\n", name);
  printf("\tnnnn: number of parts\n");
  printf("\t-bfs: use traditional iterations\n");
  printf("\t-dfs: use saturation\n");
  printf("\t-esat: use saturation by events\n");
  printf("\t-ksat: use saturation by levels\n");
  printf("\t-msat: use monolithic saturation (default)\n\n");
  printf("\t-approx: approximate the number of states (default)\n");
  printf("\t-exact:  determine the exact number of states (requires gmp)\n\n");
  printf("\t-exp: use explicit (very slow)\n\n");
  printf("\t--batch b: specify explicit batch size\n\n");
  printf("\t -l lfile: Write logging information to specified file\n\n");
  printf("\t-pdf: write the MDD representing the reachable states to Kanban.pdf\n\n");
  return 1;
}

void printStats(const char* who, const forest* f)
{
  printf("%s stats:\n", who);
  f->reportStats(meddlyout, "\t",
    HUMAN_READABLE_MEMORY  |
    BASIC_STATS | EXTRA_STATS |
    STORAGE_STATS | STORAGE_DETAILED |
    HOLE_MANAGER_STATS | HOLE_MANAGER_DETAILED
  );
}


int main(int argc, const char** argv)
{
  int N = -1;
  char method = 'm';
  unsigned batchsize = 256;
  const char* lfile = 0;
  bool build_pdf = false;
  bool approx_count = true;

  for (int i=1; i<argc; i++) {
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
    if (strcmp("-approx", argv[i])==0) {
      approx_count = true;
      continue;
    }
    if (strcmp("-exact", argv[i])==0) {
      approx_count = false;
      continue;
    }
    if (strcmp("-l", argv[i])==0) {
      lfile = argv[i+1];
      i++;
      continue;
    }
    if (strcmp("--batch", argv[i])==0) {
      i++;
      if (argv[i]) {
          int bs = atoi(argv[i]);
          if (bs>0) batchsize = unsigned(bs);
      }
      continue;
    }
    if (strcmp("-pdf", argv[i])==0) {
      build_pdf = true;
      continue;
    }
    N = atoi(argv[i]);
  }

  if (N<0) return usage(argv[0]);

  domain* d = 0;
  try {

    MEDDLY::initialize();

    timer start;

    printf("+-------------------------------------------+\n");
    printf("|   Initializing Kanban model for N = %-4d  |\n", N);
    printf("+-------------------------------------------+\n");
    fflush(stdout);

    // Initialize domain
    int* sizes = new int[16];
    for (int i=15; i>=0; i--) sizes[i] = N+1;
    d = domain::createBottomUp(sizes, 16);
    delete[] sizes;

    // Initialize forests
    forest* mdd = forest::create(d, 0, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
    forest* mxd = forest::create(d, 1, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);

    // associate loggers
    FILE_output log;
    logger* LOG = 0;
    if (lfile) {
      FILE* fout = fopen(lfile, "w");
      if (!fout) {
        printf("Couldn't open %s for writing, no logging\n", lfile);
      } else {
        log.setFILE(fout);
        LOG = new simple_logger(log);
        LOG->recordNodeCounts();
        LOG->recordTimeStamps();
        char comment[80];
        snprintf(comment, 80, "Automatically generated by kanban (N=%d)", N);
        LOG->addComment(comment);
        mdd->setLogger(LOG, "MDD");
        mxd->setLogger(LOG, "MxD");
      }
    }

    // Build initial state
    if (LOG) LOG->newPhase(mdd, "Building initial state");
    int* initial = new int[17];
    for (int i=16; i; i--) initial[i] = 0;
    initial[1] = initial[5] = initial[9] = initial[13] = N;
    dd_edge init_state(mdd);
    mdd->createEdge(&initial, 1, init_state);
    delete[] initial;

    //
    // Build next-state function
    if (LOG) LOG->newPhase(mxd, "Building next-state function");
    dd_edge nsf(mxd);
    pregen_relation* ensf = 0;
    saturation_operation* sat = 0;

    if ('s' == method) {
      ensf = new pregen_relation(mxd, 16);
    }
    if ('k' == method) {
      ensf = new pregen_relation(mxd);
    }

    if ('e' != method) {

      if (ensf) {
        start.note_time();
        buildNextStateFunction(kanban, 16, ensf, 4);
        start.note_time();
      } else {
        start.note_time();
        buildNextStateFunction(kanban, 16, mxd, nsf, 4);
        start.note_time();
#ifdef DUMP_NSF
        printf("Next-state function:\n");
        nsf.show(meddlyout, 2);
#endif
      }
      printf("Next-state function construction took %.4e seconds\n",
          start.get_last_seconds());
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
        explicitReachset(kanban, 16, mdd, init_state, reachable, batchsize);
        break;

      case 'k':
      case 's':
        printf("Building reachability set using saturation, relation");
        if ('k'==method)  printf(" by levels\n");
        else              printf(" by events\n");
        fflush(stdout);
        sat = SATURATION_FORWARD(mdd, ensf, mdd);
        if (0==sat) {
          throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
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
        start.get_last_seconds() );
    fflush(stdout);

#ifdef DUMP_REACHABLE
    printf("Reachable states:\n");
    reachable.show(meddlyout, 2);
#endif

    printStats("MDD", mdd);
    fflush(stdout);

#ifdef HAVE_LIBGMP
    if (approx_count) {
#endif
      double c;
      apply(CARDINALITY, reachable, c);
      operation::showAllComputeTables(meddlyout, 3);
      printf("Approx. %g reachable states\n", c);
#ifdef HAVE_LIBGMP
    } else {
      mpz_t c;
      mpz_init(c);
      apply(CARDINALITY, reachable, c);
      operation::showAllComputeTables(meddlyout, 3);
      printf("Exactly ");
      mpz_out_str(stdout, 10, c);
      printf(" states\n");
    }
#endif

    if (build_pdf) {
      reachable.setLabel("reachable");
      dot_maker mdd_dot(mdd, "kanban");
      mdd_dot.addRootEdge(reachable);
      mdd_dot.doneGraph();
      mdd_dot.runDot("pdf");
      if ('m' == method) {
        nsf.setLabel("next-state");
        dot_maker mxd_dot(mxd, "kanban-nsf");
        mxd_dot.addRootEdge(nsf);
        mxd_dot.doneGraph();
        mxd_dot.runDot("pdf");
      }
    }

    // cleanup
    if (LOG) {
      LOG->newPhase(mdd, "Cleanup");
      LOG->newPhase(mxd, "Cleanup");
      domain::destroy(d);
      delete LOG;
    }
    MEDDLY::cleanup();
    return 0;
  }
  catch (MEDDLY::error e) {
    printf("Caught MEDDLY error: %s\n", e.getName());
    return 1;
  }
}


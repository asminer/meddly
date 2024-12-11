
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
#include <unistd.h>

#include "../src/meddly.h"
#include "simple_model.h"

#define WRITE_MXD
#define WRITE_MDD
#define WRITE_EVMDD
// #define DEBUG_FILE
#define PROGRESS

const char* kanban[] = {
  "X-+..............",  // Tin1
  "X.-+.............",  // Tr1
  "X.+-.............",  // Tb1
  "X.-.+............",  // Tg1
  "X.....-+.........",  // Tr2
  "X.....+-.........",  // Tb2
  "X.....-.+........",  // Tg2
  "X+..--+..-+......",  // Ts1_23
  "X.........-+.....",  // Tr3
  "X.........+-.....",  // Tb3
  "X.........-.+....",  // Tg3
  "X....+..-+..--+..",  // Ts23_4
  "X.............-+.",  // Tr4
  "X.............+-.",  // Tb4
  "X............+..-",  // Tout4
  "X.............-.+"   // Tg4
};

long expected[] = {
  1, 160, 4600, 58400, 454475, 2546432, 11261376,
  41644800, 133865325, 384392800, 1005927208
};

using namespace MEDDLY;


// Build the domain
inline domain* buildKanbanDomain(int N)
{
    int sizes[16];
    for (int i=15; i>=0; i--) sizes[i] = N+1;
    return domain::createBottomUp(sizes, 16);
}

// Build the initial state
inline void buildInitial(int N, forest* mdd, dd_edge &init_state)
{
    minterm initial(mdd);
    initial.setAllVars(0);
    initial.setVar(1, N);
    initial.setVar(5, N);
    initial.setVar(9, N);
    initial.setVar(13, N);
    initial.buildFunction(init_state);
}

/*
    Generate transition relation, reachability set for given N,
    and write them to a file.
*/
long writeReachset(FILE* s, int N)
{
    // Build domain
    domain* d = buildKanbanDomain(N);

    // Build initial state
    forest* mdd = forest::create(d, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL);
    dd_edge init_state(mdd);
    buildInitial(N, mdd, init_state);

#ifdef PROGRESS
    fputc('i', stdout);
    fflush(stdout);
#endif

    // Build next-state function
    forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL);
    dd_edge nsf(mxd);
    buildNextStateFunction(kanban, 16, mxd, nsf);

#ifdef PROGRESS
    fputc('n', stdout);
    fflush(stdout);
#endif

    // Build reachable states
    dd_edge reachable(mdd);
    apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);

#ifdef PROGRESS
    fputc('r', stdout);
    fflush(stdout);
#endif

    // Build index set for reachable states
    forest* evmdd = forest::create(d, SET, range_type::INTEGER,
                        edge_labeling::INDEX_SET);
    dd_edge reach_index(evmdd);
    apply(CONVERT_TO_INDEX_SET, reachable, reach_index);

#ifdef PROGRESS
    fputc('x', stdout);
    fflush(stdout);
#endif

    long c;
    apply(CARDINALITY, reachable, c);

#ifdef PROGRESS
    fputc('c', stdout);
    fflush(stdout);
#endif

    FILE_output mys(s);
#ifdef WRITE_MXD
    mdd_writer wmxd(mys, mxd);
    wmxd.writeRootEdge(nsf);
    wmxd.finish();
#endif

#ifdef WRITE_MDD
    // Write initial state + reachable states
    mdd_writer wmdd(mys, mdd);
    wmdd.writeRootEdge(init_state);
    wmdd.writeRootEdge(reachable);
    wmdd.finish();
#endif

#ifdef WRITE_EVMDD
    mdd_writer wevmdd(mys, evmdd);
    wevmdd.writeRootEdge(reach_index);
    wevmdd.finish();
#endif

    domain::destroy(d);

    return c;
}


/*
    Generate transition relation, reachability set for given N,
    then read from a file and check for equality.
*/
bool generateAndRead(FILE* s, int N)
{
  // Build domain
  domain* d = buildKanbanDomain(N);

  // Build initial state
  forest* mdd = forest::create(d, SET, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
  dd_edge init_state(mdd);
  buildInitial(N, mdd, init_state);

  // Build next-state function
  forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
  dd_edge nsf(mxd);
  buildNextStateFunction(kanban, 16, mxd, nsf);

  // Build reachable states
  dd_edge reachable(mdd);
  apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);

  // Build index set for reachable states
  forest* evmdd = forest::create(d, SET, range_type::INTEGER, edge_labeling::INDEX_SET);
  dd_edge reach_index(evmdd);
  apply(CONVERT_TO_INDEX_SET, reachable, reach_index);

  // Now, read from the file and verify
  FILE_input mys(s);

  dd_edge list[2];
#ifdef WRITE_MXD
  mdd_reader R_mxd(mys, mxd);
  dd_edge r_nsf;
  R_mxd.readRootEdge(r_nsf);
  if (r_nsf != nsf) {
    printf("Failed to generate and read MXD\n");
    return false;
  }
#endif

#ifdef WRITE_MDD
  mdd_reader R_mdd(mys, mdd);
  dd_edge r_init, r_reach;
  R_mdd.readRootEdge(r_init);
  R_mdd.readRootEdge(r_reach);
  if (r_init != init_state) {
    printf("Failed to generate and read initial state\n");
    return false;
  }
  if (r_reach != reachable) {
    printf("Failed to generate and read reachable states\n");
    return false;
  }
#endif

#ifdef WRITE_EVMDD
  mdd_reader R_evmdd(mys, evmdd);
  dd_edge r_index;
  R_evmdd.readRootEdge(r_index);
  if (r_index != reach_index) {
    printf("Failed to generate and read reachable state indexes\n");
    FILE_output myout(stdout);
    myout << "reach_index: " << reach_index << "\n";
    myout << "from file  : " << r_index << "\n";
#ifdef DEBUG_FILE
    reach_index.showGraph(myout);
    r_index.showGraph(myout);
#endif
    return false;
  }
#endif

  domain::destroy(d);
  return true;
}


/*
    Read the transition relation, reachability set from a file,
    then generate them for a given N and check for equality.
*/
bool readAndGenerate(FILE* s, int N)
{
  FILE_input mys(s);

  // Build domain
  domain* d = buildKanbanDomain(N);

  // Build (empty) forests
  forest* mdd = forest::create(d, SET, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
  forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
  forest* evmdd = forest::create(d, SET, range_type::INTEGER, edge_labeling::INDEX_SET);

  dd_edge mxdfile;
  dd_edge mddfile0, mddfile1;
  dd_edge indexfile;
  // Read from the file
#ifdef WRITE_MXD
  mdd_reader R_mxd(mys, mxd);
  R_mxd.readRootEdge(mxdfile);
#endif
#ifdef WRITE_MDD
  mdd_reader R_mdd(mys, mdd);
  R_mdd.readRootEdge(mddfile0);
  R_mdd.readRootEdge(mddfile1);
#endif
#ifdef WRITE_EVMDD
  mdd_reader R_evmdd(mys, evmdd);
  R_evmdd.readRootEdge(indexfile);
#endif

  // Build initial state
  dd_edge init_state(mdd);
  buildInitial(N, mdd, init_state);

  // Build next-state function
  dd_edge nsf(mxd);
  buildNextStateFunction(kanban, 16, mxd, nsf);

  // Build reachable states
  dd_edge reachable(mdd);
  apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);

  // Build index set for reachable states
  dd_edge reach_index(evmdd);
  apply(CONVERT_TO_INDEX_SET, reachable, reach_index);

  // Verify
#ifdef WRITE_MXD
  if (mxdfile != nsf) {
    printf("Failed to read and generate MXD\n");
    return false;
  }
#endif
#ifdef WRITE_MDD
  if (mddfile0 != init_state) {
    printf("Failed to generate and read initial state\n");
    return false;
  }
  if (mddfile1 != reachable) {
    printf("Failed to generate and read reachable states\n");
    return false;
  }
#endif
#ifdef WRITE_EVMDD
  if (indexfile != reach_index) {
    printf("Failed to generate and read reachable state indexes\n");
    return false;
  }
#endif

  domain::destroy(d);
  return true;
}



/*
    Main program, of course
*/
int main()
{
#ifdef DEBUG_FILE
  const int N = 2;
#else
  const int N = 11;
#endif
  MEDDLY::initialize();

  char filename[20];
  strcpy(filename, "kan_io.data.XXXXXX");
  mkstemp(filename);
#ifdef DEBUG_FILE
  printf("Creating file %s\n", filename);
#endif

  FILE* s = fopen(filename, "w");
  if (0==s) {
    printf("Couldn't open file %s for writing\n", filename);
    return 2;
  }

  try {
    printf("Saving dds to file...\n");
    for (int n=1; n<N; n++) {
      printf("N=%2d:  ", n);
      fflush(stdout);
      long c = writeReachset(s, n);
      printf("%12ld states\n", c);
      if (c !=expected[n]) {
        printf("Wrong number of states!\n");
        throw 1;
      }
    }
    fclose(s);

    // read back the file into already built forest

    s = fopen(filename, "r");
    if (0==s) {
      printf("Couldn't open file %s for reading\n", filename);
      throw 2;
    }

    printf("Generate and read...\n");
    for (int n=1; n<N; n++) {
      printf("N=%2d:  ", n);
      fflush(stdout);
      if (!generateAndRead(s, n)) throw 3;
      printf("%19s\n", "verified");
    }
    fclose(s);

    // read back the file and then build up forest

    s = fopen(filename, "r");
    if (0==s) {
      printf("Couldn't open file %s for reading\n", filename);
      throw 2;
    }

    printf("Read and generate...\n");
    for (int n=1; n<N; n++) {
      printf("N=%2d:  ", n);
      fflush(stdout);
      if (!readAndGenerate(s, n)) throw 4;
      printf("%19s\n", "verified");
    }
    fclose(s);


    // Cleanup
#ifndef DEBUG_FILE
    remove(filename);
#endif
    MEDDLY::cleanup();
    printf("Done\n");
    return 0;
  }
  catch (int e) {
    // cleanup
#ifndef DEBUG_FILE
    remove(filename);
#endif
    MEDDLY::cleanup();
    printf("\nError %d\n", e);
    return e;
  }
  catch (MEDDLY::error e) {
    // cleanup
#ifndef DEBUG_FILE
    remove(filename);
#endif
    MEDDLY::cleanup();
    printf("\nError: %s, thrown at %s line %u\n", e.getName(),
            e.getFile(), e.getLine());
    return 1;
  }
  catch (...) {
    // cleanup
#ifndef DEBUG_FILE
    remove(filename);
#endif
    MEDDLY::cleanup();
    printf("\nFailed (caught exception)\n");
    return 1;
  }
}


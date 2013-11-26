
// $Id$

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
    This program reads in a grid representing a generation in
    Conway's Game Of Life, and determines possible previous generations.
*/

#include <iostream>
#include <string.h>

#include "meddly.h"
#include "meddly_expert.h"

// #define DEBUG_INPUT
// #define DEBUG_NUM_TRUE
// #define DEBUG_CONSTRAINTS

struct lifegrid {
    int width;
    int height;
    bool** grid;
  public:
    lifegrid();
    ~lifegrid();
    inline int getLevelOf(int r, int c) const {
      if (r<0 || r>=height) return 0;
      if (c<0 || c>=width) return 0;
      return r*width + c + 1;
    }

  private:
    void operator=(const lifegrid &) { }
    lifegrid(const lifegrid &) { }
};

lifegrid::lifegrid()
{
  width = 0;
  height = 0;
  grid = 0;
}

lifegrid::~lifegrid()
{
  if (grid) {
    for (int i=0; i<height; i++) delete[] grid[i];
  }
  delete[] grid;
}

const char* stripPath(const char* name)
{
  const char* ptr;
  for (ptr=name; *ptr; ptr++) {
    if ('/' == *ptr) name = ptr+1;
  }
  return name;
}

inline bool isAlive(char c)
{
  if (' ' == c) return 0;
  if ('.' == c) return 0;
  if ('`' == c) return 0;
  if (',' == c) return 0;
  if ('\'' == c) return 0;
  if ('\r' == c) return 0;
  if ('\n' == c) return 0;
  return 1;
}

int helpScreen(const char* who)
{
  using namespace std;
  cerr << "\n    Usage: " << stripPath(who) << " [switches]\n";
  cerr << "\n    Reads a grid, from standard input, representing a generation\n";
  cerr << "    of Conway's Game of Life.  Then...\n";
  cerr << "\n";
  cerr << "        -h       This help screen\n";
  cerr << "\n";
  cerr << "        -a       Show ALL possible previous generations\n";
  cerr << "        -m       Show all minimal previous generations\n";
  cerr << "        -s       Show some minimal previous generation (default)\n";
  cerr << "\n";
  return 0;
}

void parseInput(lifegrid &G)
{
  using namespace std;

  //
  // Parse the input
  //

#ifdef DEBUG_INPUT
  cerr << "Parsing grid\n";
#endif

  const int maxline = 65536;
  char line[maxline];
  bool proper_grid = true;
  int width = -1;
  int maxrows = 256;
  bool** grid = new bool*[maxrows];

  int height;
  for (height=0; ; height++) {
    //
    // Grab the current line
    //
    cin.getline(line, maxline);
    if (cin.eof()) break;
    int this_width = strlen(line);
    if (width < 0) {
      width = this_width;
    } else if (proper_grid) {
      if (width != this_width) {
        cerr << "Error: lines are not all the same width\n";
        exit(1);
      }
    }
    //
    // Build the life row for this line
    //
    bool* row = new bool[this_width];
    for (int i=0; i<this_width; i++) {
      row[i] = isAlive(line[i]);
    }
    //
    // Add to the grid
    //
    if (height+1 >= maxrows) {
      // expand
      maxrows += 256;
      bool** newgrid = new bool*[maxrows];
      for (int i=0; i<height; i++) {
        newgrid[i] = grid[i];
      }
      delete[] grid;
      grid = newgrid;
    }
    grid[height] = row;

  } // for
  //
  // Fill the struct
  //
  G.width = width;
  G.height = height;
  G.grid = grid;
  //
  // Zero out the rest of the grid
  //
  for (; height<maxrows; height++) {
    grid[height] = 0;
  }
}

using namespace MEDDLY;

void numTrueGeneral(const int* levels, int* bottom, int n, int deflt, dd_edge &e)
{
#ifdef DEBUG_NUM_TRUE
  using namespace std;
  cerr << "inside numTrueGeneral, levels array: [";
  for (int i=0; ; i++) {
    if (i) cerr << ", ";
    cerr << levels[i];
    if (0==levels[i]) break;
  }
  cerr << "]\n";
#endif
  //
  // Get length of levels array
  //
  int num_levels = 0;
  for (; levels[num_levels]; num_levels++);

  //
  // Build temp storage for "current level"
  //
  int* current = new int[n];
  expert_forest* ef = (expert_forest*) e.getForest();

  //
  // Loop through levels, from bottom up, and build
  // the DD as follows.
  // For each level, we have a node for "i set vars",
  // for various possible i.
  // The node at this level for i, has a pointer to i
  // and to i+1 below.
  //
  for (int ki = num_levels-1; ki>=0; ki--) {
#ifdef DEBUG_NUM_TRUE
    cerr << "building for level " << levels[ki] << "\n";
    cerr << "  below nodes: [" << bottom[0];
    for (int i=1; i<n; i++) {
      cerr << ", " << bottom[i];
    }
    cerr << "]\n";
#endif
    for (int i=0; i<n; i++) { 
      node_builder& nb = ef->useNodeBuilder(levels[ki], 2);
      nb.d(0) = ef->linkNode(bottom[i]);
      nb.d(1) = ef->linkNode( (i+1<n) ? bottom[i+1] : deflt );
      current[i] = ef->createReducedNode(-1, nb);
    } // for i
    //
    // now, recycle the level below us
    //
    for (int i=0; i<n; i++) {
      ef->unlinkNode(bottom[i]);
      bottom[i] = current[i];
    }
  }

#ifdef DEBUG_NUM_TRUE
  cerr << "  final nodes: [ (" << bottom[0] << ") ";
  for (int i=1; i<n; i++) {
    cerr << ", " << bottom[i];
  }
  cerr << "]\n";
#endif

  //
  // bottom[0] is our answer, set it up
  //
  e.set(bottom[0], 0);

  //
  // Recycle all other nodes
  //
  for (int i=1; i<n; i++) ef->unlinkNode(bottom[i]);

  // 
  // Cleanup
  //
  delete[] current;

#ifdef DEBUG_NUM_TRUE
  cerr << "leaving numTrueGeneral, answer is " << bottom[0] << ":\n";
  cerr.flush();
  e.show(stderr, 2);
#endif
}

void numTrueEquals(const int* levels, int n, dd_edge &e)
{
  expert_forest* ef = (expert_forest*) e.getForest();
  int* bottom = new int[n+1];
  for (int i=0; i<n; i++) bottom[i] = 0;
  bottom[n] = ef->getTerminalNode(1);
  numTrueGeneral(levels, bottom, n+1, 0, e);
  delete[] bottom;
}

void numTrueGreaterThan(const int* levels, int n, dd_edge &e)
{
  expert_forest* ef = (expert_forest*) e.getForest();
  int* bottom = new int[n+1];
  for (int i=0; i<n+1; i++) bottom[i] = 0;
  numTrueGeneral(levels, bottom, n+1, ef->getTerminalNode(1), e);
  delete[] bottom;
}

void numTrueLessThan(const int* levels, int n, dd_edge &e)
{
  expert_forest* ef = (expert_forest*) e.getForest();
  int* bottom = new int[n];
  for (int i=0; i<n; i++) bottom[i] = ef->getTerminalNode(1);
  numTrueGeneral(levels, bottom, n, 0, e);
  delete[] bottom;
}

void buildNeighborList(const lifegrid &G, int r, int c, int* levels)
{
  int p = 0;
  levels[p] = G.getLevelOf(r+1, c+1);
  if (levels[p]) p++;
  levels[p] = G.getLevelOf(r+1, c);
  if (levels[p]) p++;
  levels[p] = G.getLevelOf(r+1, c-1);
  if (levels[p]) p++;
  levels[p] = G.getLevelOf(r, c+1);
  if (levels[p]) p++;
  levels[p] = G.getLevelOf(r, c-1);
  if (levels[p]) p++;
  levels[p] = G.getLevelOf(r-1, c+1);
  if (levels[p]) p++;
  levels[p] = G.getLevelOf(r-1, c);
  if (levels[p]) p++;
  levels[p] = G.getLevelOf(r-1, c-1);
  if (levels[p]) p++;
  levels[p] = 0;
}

void buildConstraints(const lifegrid &G, int row, int col, dd_edge &e)
{
  if (!G.getLevelOf(row, col)) return;

  static int neighbors[9];
  buildNeighborList(G, row, col, neighbors);

#ifdef DEBUG_CONSTRAINTS
  using namespace std;
  cerr << "Building constraint for (" << row << ", " << col << "), ";
  if (G.grid[row][col]) cerr << "alive\n"; else cerr << "dead\n";
#endif

  forest* f = e.getForest();

  if (G.grid[row][col]) {
    //
    // Alive.
    //
    // We're alive in this generation iff, in the previous generation:
    //
    //    (we're alive AND we have 2 alive neighbors) OR
    //    (3 alive neighbors)
    //

    dd_edge alive(f);
    f->createEdgeForVar(G.getLevelOf(row, col), false, alive);
    dd_edge two_N(f);
    numTrueEquals(neighbors, 2, two_N);
    dd_edge three_N(f);
    numTrueEquals(neighbors, 3, three_N);

    apply(INTERSECTION, alive, two_N, e);
    apply(UNION, e, three_N, e);

  } else {
    //
    // Dead.
    //
    // We're dead in this generation iff, in the previous generation:
    //
    //  (less than 2 alive neighbors) OR
    //  (more than 3 alive neighbors) OR
    //  (we're dead AND we have 2 alive neighbors)
    // 

    bool inv_terms[2] = {true, false};
    dd_edge dead(f);
    f->createEdgeForVar(G.getLevelOf(row, col), false, inv_terms, dead);
    dd_edge lt_2(f);
    numTrueLessThan(neighbors, 2, lt_2);
    dd_edge gt_3(f);
    numTrueGreaterThan(neighbors, 3, gt_3);
    dd_edge two_N(f);
    numTrueEquals(neighbors, 2, two_N);

    apply(INTERSECTION, dead, two_N, e);
    apply(UNION, e, gt_3, e);
    apply(UNION, e, lt_2, e);

  }

#ifdef DEBUG_CONSTRAINTS
  cerr << "finished constraint for (" << row << ", " << col << "), ";
  if (G.grid[row][col]) cerr << "alive\n"; else cerr << "dead\n";
  cerr.flush();
  e.show(stderr, 2);
#endif
}


void AndAlongRowsThenCols(dd_edge *C, const lifegrid &G)
{
  std::cerr << "Combining constraints for each row\n";
  for (int r=0; r<G.height; r++) {
    for (int c=1; c<G.width; c++) {
      int i = G.getLevelOf(r, c-1);
      int j = G.getLevelOf(r, c);
      std::cerr << ".";
      std::cerr.flush();
      apply(INTERSECTION, C[i], C[j], C[j]);
      C[i].clear();
    } // for c
    std::cerr << "\n";
  } // for r

  std::cerr << "Combining overall row constraints\n";
  // TBD
}

void AndAlongColsThenRows(dd_edge *C, const lifegrid &G)
{
  std::cerr << "Combining constraints for each column\n";
  for (int c=0; c<G.width; c++) {
    for (int r=1; r<G.height; r++) {
      int i = G.getLevelOf(r-1, c);
      int j = G.getLevelOf(r, c);
      std::cerr << ".";
      std::cerr.flush();
      apply(INTERSECTION, C[i], C[j], C[j]);
      C[i].clear();
    } // for c
    std::cerr << "\n";
  } // for r

  std::cerr << "Combining overall column constraints\n\t";
  int r = G.height-1;
  for (int c=1; c<G.width; c++) {
    int i = G.getLevelOf(r, c-1);
    int j = G.getLevelOf(r, c);
    if (c>1) std::cerr << ", ";
    std::cerr << c;
    std::cerr.flush();
    apply(INTERSECTION, C[i], C[j], C[j]);
    C[i].clear();
  }
  std::cerr << "\nDone!\n";
}

void showGeneration(const lifegrid &G, const int* minterm, const char* prefix)
{
  using namespace std;
  for (int j=0; j<G.height; j++) {
    if (prefix) cout << prefix;
    for (int i=0; i<G.width; i++) {
      cout << (minterm[G.getLevelOf(j, i)] ? 'A' : '.');
    } // for i
    cout << "\n";
  } // for j
}

int sum(const int* minterm, int N)
{
  int a = 0;
  for (; N; N--) {
    a += minterm[N];
  }
  return a;
}

int main(int argc, const char** argv)
{
  enum {
    all_previous, all_minimal, some_minimal
  } display = some_minimal;
  using namespace std;
  initialize();
  cerr << "Using " << getLibraryInfo(0) << "\n";

  for (int i=1; i<argc; i++) {
    if (argv[i][0] != '-' || argv[i][2] != 0) {
      return 1+helpScreen(argv[0]);
    }
    switch (argv[i][1]) {
      case 'h': return helpScreen(argv[0]);
      case 'a':
          display = all_previous;
          break;
      case 'm':
          display = all_minimal;
          break;
      case 's':
          display = some_minimal;
          break;

      default : return 1+helpScreen(argv[0]);
    } // switch
  } // for i

  //
  // Grab the grid
  //
  lifegrid G;
  parseInput(G);
  cerr << "Got " << G.width << "x" << G.height << " grid\n";
  if (0==G.width || 0==G.height) return 0;

#ifdef DEBUG_INPUT
  for (int j=0; j<G.height; j++) {
    for (int i=0; i<G.width; i++) {
      cerr << (G.grid[j][i] ? 'A' : ' ');
    } // for i
    cerr << "\n";
  } // for j
#endif

  //
  // Build a forest
  //
  int N = G.height * G.width;
  int* scratch = new int[N+1];
  for (int i=0; i<=N; i++) scratch[i] = 2;
  domain* d = createDomainBottomUp(scratch, N);
  assert(d);
  forest* f = 
    d->createForest(false, forest::BOOLEAN, forest::MULTI_TERMINAL);
  assert(f);

  //
  // Build constraints for each cell
  //
  cerr << "Building constraints for each cell   ";
  cerr.flush();
  dd_edge* C = new dd_edge[N+1];  
  dd_edge answer(f);
  for (int i=0; i<=N; i++) C[i] = answer;    // set forest
  for (int r=0; r<G.height; r++) {
    for (int c=0; c<G.width; c++) {
      int i = G.getLevelOf(r, c);
      buildConstraints(G, r, c, C[i]);
    } // for c
  } // for r
  cerr << "    Done\n";

  //
  // Combine constraints
  //
  AndAlongColsThenRows(C, G);
  answer = C[N];
  C[N].clear();
  delete[] C;

  //
  //
  // Outputs
  //
  //

  long c;
  apply(CARDINALITY, answer, c);
  cerr << "There are " << c << " possible previous generations\n";
  if (0==c) {
    // Bail!
    cleanup();
    return 0;
  }

  //
  // Compute stats and maybe show all (selected by a switch)
  // 

  int min = N;
  int max = 0;
  int num_max = 0;
  int num_min = 0;
  enumerator first(answer);
  long n = 0;
  for (; first; ++first) {
      const int* minterm = first.getAssignments();
      int alive = sum(minterm, N);
      if (alive==max) num_max++;
      if (alive>max) {
        max = alive;
        num_max = 1;
      }
      if (alive==min) num_min++;
      if (alive<min) {
        min = alive;
        num_min = 1;
      }
      n++;
      if (all_previous == display) {
        cout << "Previous #" << n << "\n";
        showGeneration(G, minterm, 0);      
      }
  }
  cerr << "Max alive cells: " << max << " (" << num_max << " ways)\n";
  cerr << "Min alive cells: " << min << " (" << num_min << " ways)\n";

  if (all_minimal == display || some_minimal == display) {
    first.start(answer);
    n=0;
    for (; first; ++first) {
      const int* minterm = first.getAssignments();
      int alive = sum(minterm, N);
      if (alive > min) continue;
      if (all_minimal == display) {
        n++;
        cout << "Minimal previous #" << n << "\n";
      } else {
        cerr << "One minimal previous generation:\n";
      }
      showGeneration(G, minterm, 0);      
      if (some_minimal == display) break;
    }
  }

  // Done!
  cleanup();
  return 0;
}



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
    Builds the set of solutions to the queen cover problem for 
    user-specified board size NxN.

    In other words, finds all possible ways to put queens onto
    an NxN chessboard so that every square either contains a queen,
    or is attached by a queen.

    State: 
      For each square, does it hold a queen.

    Constraints:
      conjunction over all squares:
        is this square covered

*/

#include <cstdlib>
#include <fstream>

#include "meddly.h"
#include "meddly_expert.h"
#include "loggers.h"
#include "timer.h"

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <queue>

// #define USE_EMPTY_CHECK
// #define PUSH_ROW_CONSTRAINTS
#define USE_SMALLER_SPANS

using namespace MEDDLY;

int N;
FILE* outfile;
const char* lfile;

dd_edge** qic;
dd_edge** qidp;
dd_edge** qidm;

bool use_batches;
bool inner_batch;
bool outer_batch;
bool monolithic_batch;

typedef enum {
  node_count,
  edge_count,
  cardinality,
  height,
  n_ary_intersection,
  spans,
  accumulate,
  folding
} batch_heuristic;
batch_heuristic heuristic;

typedef struct {
  node_handle a;
  node_handle b;
} node_pair;

class sorter {
  public:
    sorter(forest* f, batch_heuristic h);
    void setLevel(int k) { level = k; }
    void setDescending(bool d) { smaller_is_better = d; }

    template<typename T>
      bool lessThan(MEDDLY::node_handle a, MEDDLY::node_handle b) const {
        T val_a;
        T val_b;
        op->compute(level, a, val_a);
        op->compute(level, b, val_b);
        return (val_a < val_b);
      }

    // returns true if a must be placed before b in the sorted order
    bool operator()(MEDDLY::node_handle a, MEDDLY::node_handle b) const {
      // returns true if a < b
      // note: priority queue places the largest item on the top
      // so a < b means do b first, than do a.
      if (a == b) return true;
      if (a == -1 || b == 0) return smaller_is_better;
      if (a == 0 || b == -1) return !smaller_is_better;
      if (op_type == INTEGER)
        return lessThan<long int>(a, b) != smaller_is_better;
      else 
        return lessThan<double>(a, b) != smaller_is_better;
    }

  protected:
    MEDDLY::expert_forest* f;
    batch_heuristic heuristic;
    MEDDLY::opnd_type op_type;
    MEDDLY::unary_operation *op;
    int level;
    bool smaller_is_better;
};

class span_checker {
  public:
    span_checker(forest* f);
    void gc();
    unsigned getSpan(node_handle node);
    std::map<node_handle, unsigned> getNode2SpanMap(std::set<node_handle>& v);
    node_pair getBestPair(std::set<node_handle>& s);

  protected:
    MEDDLY::expert_forest* f;
    std::map<node_handle, unsigned> span_cache;
    std::map<node_pair, unsigned> pair_span_cache;
};

class emptiness_checker {
  public:
    emptiness_checker(forest* f) {
      this->f = dynamic_cast<expert_forest*>(f);
      mddMultiply = getOperation(MULTIPLY, f, f, f);
      assert(mddMultiply);
      n_constraints = 0;
      pings = hits = 0;
    }

    ~emptiness_checker() {
      for (auto i : intersect_cache) {
        for (auto j : i.first) {
          f->unlinkNode(j);
        }
        f->unlinkNode(i.second);
      }
    }

    bool isIntersectionEmpty(std::vector<node_handle>& v) {
#ifdef USE_EMPTY_CHECK
      bool is_empty = false;

      std::set<node_handle> s;
      for (unsigned i = 0; i < v.size() && !is_empty; ) {
        s.insert(v[i++]);
        s.insert(v[i++]);
        is_empty = isEmpty(s);
      }
      return is_empty;
#else
      return false;
#endif
    }

    node_handle getIntersection(std::vector<node_handle>& v) {
      std::set<node_handle> s;
      s.insert(v.begin(), v.end());
      return intersect(s);
    }

    unsigned getPings() const { return pings; }

    unsigned getHits() const { return hits; }

  protected:
    void saveResult(std::set<node_handle>& s, bool r) {
      cache[s] = r;
    }

    bool findResult(std::set<node_handle>& s, bool& r) {
      pings++;
      std::map<std::set<node_handle>, bool>::iterator it = cache.find(s);
      if (it != cache.end()) {
        r = it->second;
        hits++;
        return true;
      }
      return false;
    }

    void clearResult(std::set<node_handle>& s) {
      cache.erase(s);
    }

    void saveIntersectResult(std::set<node_handle>& s, node_handle r) {
      for (auto i : s) { f->linkNode(i); }
      intersect_cache[s] = f->linkNode(r);
#if 0
      printf("Saved: intersection(");
      for (auto i : s) { printf(" %d", i); }
      printf(" ) = %d\n", r);
#endif
    }

    bool findIntersectResult(std::set<node_handle>& s, node_handle& r) {
      pings++;
      std::map<std::set<node_handle>, node_handle>::iterator it = intersect_cache.find(s);
      if (it != intersect_cache.end()) {
        r = f->linkNode(it->second);
        hits++;
        return true;
      }
      return false;
    }

    void clearIntersectResult(std::set<node_handle>& s) {
      std::map<std::set<node_handle>, node_handle>::iterator it = intersect_cache.find(s);
      if (it != intersect_cache.end()) {
        f->unlinkNode(it->second);
        for (auto i : s) {
          f->unlinkNode(i);
        }
        intersect_cache.erase(it);
      }
    }

    bool isEmpty(std::set<node_handle>& s) {
      if (s.empty()) return true;
      if (s.size() == 1) {
        return 0 == *(s.begin());
      }
      int k = 0;
      for (auto n : s) {
        if (n == 0) return true;
        if (k < f->getNodeLevel(n)) k = f->getNodeLevel(n);
      }

      bool is_empty = true;
      if (findResult(s, is_empty)) return is_empty;

      // for each i in 0..n-1, where n==levelSize(k)
      //    build a set s.t. s_next[j] = s[j][i]
      //    if (isEmpty(s_next)) return true; else continue;

      // a s_next per i
      std::vector<std::set<node_handle>> s_next(f->getLevelSize(k));
      bool is_zero[f->getLevelSize(k)] = {false};
      for (auto n : s) {
        if (f->getNodeLevel(n) == k) {
          unpacked_node* node = unpacked_node::newFromNode(f, n, true);
          for (int i = 0; i < node->getSize(); i++) {
            if (is_zero[i]) continue;
            if (node->d(i)) {
              s_next[i].insert(node->d(i));
            } else {
              is_zero[i] = true;
              s_next[i].clear();
            }
          }
          unpacked_node::recycle(node);
        } else {
          MEDDLY_DCASSERT(f->isFullyReduced());
          unsigned index = 0;
          for (auto i : s_next) if (!is_zero[index++]) i.insert(n);
        }
      }

      for (unsigned index = 0; index < s_next.size() && is_empty; index++) {
        if (!is_zero[index]) is_empty = isEmpty(s_next[index]);
#if 0
        printf("k: %d, index: %d, [ ", k, index);
        for (auto j : s_next[index]) {
          printf("%d ", j);
        }
        printf("]: %s\n", (is_empty? "empty": "not-empty"));
#endif
      }

      saveResult(s, is_empty);
      return is_empty;
    }

    node_handle intersect(std::set<node_handle>& s) {
#if 0
      printf("At [ ");
      for (auto j : s) {
        printf("%d:%d ", j, f->getNodeLevel(j));
      }
      printf("]\n");
#endif

      if (s.empty()) return 0;
      if (s.size() == 1) return f->linkNode(*s.begin());

      int k = 0;
      for (auto n : s) {
        if (n == 0) return 0;
        if (k < f->getNodeLevel(n)) k = f->getNodeLevel(n);
      }

      assert(k != 0);

      node_handle result = 0;
#if 1
      if (findIntersectResult(s, result)) {
        // printf("Found result: %d\n", result);
        return result;
      }
#endif

      // for each i in 0..n-1, where n==levelSize(k)
      //    build a set s.t. s_next[j] = s[j][i]
      //    if (isEmpty(s_next)) return true; else continue;

      // a s_next per i
      std::vector<std::set<node_handle>> s_next(f->getLevelSize(k));
      bool is_zero[f->getLevelSize(k)] = {false};
      for (auto n : s) {
        if (f->getNodeLevel(n) == k) {
          unpacked_node* node = unpacked_node::newFromNode(f, n, true);
          for (int i = 0; i < node->getSize(); i++) {
            if (is_zero[i]) continue;
            if (node->d(i)) {
              s_next[i].insert(f->linkNode(node->d(i)));
            } else {
              is_zero[i] = true;
              for (auto z : s_next[i]) {
                f->unlinkNode(z);
              }
              s_next[i].clear();
            }
          }
          unpacked_node::recycle(node);
        } else {
          MEDDLY_DCASSERT(f->isFullyReduced());
          for (int i = 0; i < f->getLevelSize(k); i++) {
            if (is_zero[i]) continue;
            s_next[i].insert(f->linkNode(n));
          }
        }
      }

#if 0
      for (unsigned index = 0; index < s_next.size(); index++) {
        printf("k: %d, index: %d, [ ", k, index);
        for (auto j : s_next[index]) {
          printf("%d ", j);
        }
        printf("]\n");
      }
#endif

      unpacked_node* result_node = unpacked_node::newFull(f, k, s_next.size());
      for (unsigned index = 0; index < s_next.size(); index++) {
        result_node->d_ref(index) = is_zero[index]? 0: intersect(s_next[index]);
      }

      for (auto i : s_next) {
        for (auto j : i) {
          if (j) f->unlinkNode(j);
        }
      }

      result = f->createReducedNode(-1, result_node);
#if 1
      saveIntersectResult(s, result);
#endif
      return result;
    }

  protected:
    expert_forest* f;
    binary_operation* mddMultiply;
    unsigned n_constraints;
    std::map<std::set<node_handle>, bool> cache;
    std::map<std::set<node_handle>, node_handle> intersect_cache;
    unsigned pings, hits;

};

forest* buildQueenForest(forest::policies &p)
{
  printf("Initializing domain and forest\n");
  const char* ndp = "unknown node deletion";
  switch (p.deletion) {
    case forest::policies::NEVER_DELETE:
        ndp = "never delete";
        break;

    case forest::policies::OPTIMISTIC_DELETION:
        ndp = "optimistic node deletion";
        break;

    case forest::policies::PESSIMISTIC_DELETION:
        ndp = "pessimistic node deletion";
        break;
  }
  printf("\tUsing %s policy\n", ndp);


  int* vars = new int[N*N];
  for (int i=0; i<N*N; i++) {
    vars[i] = 2;
  }
  domain* d = createDomainBottomUp(vars, N*N);
  assert(d);
#if 1
  forest* f = 
    d->createForest(false, forest::INTEGER, forest::MULTI_TERMINAL, p);
#else
  p.setQuasiReduced();
  forest* f = 
    d->createForest(false, forest::INTEGER, forest::MULTI_TERMINAL, p);
#endif

  assert(f);

  delete[] vars;
  return f;
}

inline int ijmap(int i, int j)
{
  return i*N+j+1;
}

void hasQueen(forest* f, int i, int j, dd_edge &e)
{
  f->createEdgeForVar(ijmap(i, j), false, e);
}

void queenInRow(forest* f, int r, dd_edge &e)
{
  f->createEdge(int(0), e);
  if (r<0 || r>=N) return;
  for (int j=0; j<N; j++) {
    dd_edge tmp(f);
    hasQueen(f, r, j, tmp);
    apply(MAXIMUM, e, tmp, e);
  }
}

void queenInCol(forest* f, int c, dd_edge &e)
{
  if (c>=0 && c<N && qic[c]) {
    e = *qic[c];
    return;
  }
  f->createEdge(int(0), e);
  if (c<0 || c>=N) return;
  for (int i=0; i<N; i++) {
    dd_edge tmp(f);
    hasQueen(f, i, c, tmp);
    apply(MAXIMUM, e, tmp, e);
  }
  qic[c] = new dd_edge(e);
}

void queenInDiagP(forest* f, int d, dd_edge &e)
{
  if (d>=0 && d<2*N-1 && qidp[d]) {
    e = *qidp[d];
    return;
  }
  f->createEdge(int(0), e);
  if (d<0 || d>=2*N-1) return;
  for (int i=0; i<N; i++) for (int j=0; j<N; j++) if (i+j==d) {
    dd_edge tmp(f);
    hasQueen(f, i, j, tmp);
    apply(MAXIMUM, e, tmp, e);
  }
  qidp[d] = new dd_edge(e);
}

void queenInDiagM(forest* f, int d, dd_edge &e)
{
  if (d>-N && d<N && qidm[d+N-1]) {
    e = *qidm[d+N-1];
    return;
  }
  f->createEdge(int(0), e);
  if (d<=-N || d>=N) return;
  for (int i=0; i<N; i++) for (int j=0; j<N; j++) if (i-j==d) {
    dd_edge tmp(f);
    hasQueen(f, i, j, tmp);
    apply(MAXIMUM, e, tmp, e);
  }
  qidm[d+N-1] = new dd_edge(e);
}

node_handle batchMultiply(
    expert_forest* f,
    batch_heuristic heuristic,
    std::vector<node_handle>& v)
{
  if (v.empty()) return 0;
  if (v.size() == 1) {
    node_handle result = v[0];
    v.clear();
    return result;
  }

  if (heuristic == n_ary_intersection) {
    emptiness_checker ec(f);
    node_handle result = ec.getIntersection(v);

    printf("Emptiness Checker Stats:\n");
    printf("\t\t Pings: %u\n", ec.getPings());
    printf("\t\t Hits: %u\n", ec.getHits());
    printf("\t\t Ratio: %lf\n", (double)ec.getHits()/(double)ec.getPings());

    return result;
  }
  
  binary_operation* mddMultiply = getOperation(MULTIPLY, f, f, f);
  assert(mddMultiply);

  if (heuristic == accumulate) {
    node_handle result = f->linkNode(v[0]);
    for (auto i : v) {
      if (result) {
        node_handle temp = result;
        result = mddMultiply->compute(result, i);
        f->unlinkNode(temp);
      }
      f->unlinkNode(i);
    }
    return result;
  }

  if (heuristic == folding) {
    assert(v.size() >= 2);
    while (v.size() > 1) {
      std::set<node_handle> s;
      node_handle result = v[0];
      printf("\n");
      for (unsigned i = 0; (i+1) < v.size() && result; ) {
        printf("\rFolding %lu nodes: [", v.size());
        for (unsigned j = 0; j < i; j++) printf("=");
        fflush(stdout);

        result = mddMultiply->compute(v[i], v[i+1]);
        s.insert(result);
        i += 2;

        printf("==>");
        for (unsigned j = i+1; j < v.size(); j++) printf(" ");
        printf("]");
        fflush(stdout);
      }
      printf("\n");

      // unlink nodes in v
      for (auto i : v) f->unlinkNode(i);

      // return 0 if result is 0
      if (0 == result) {
        // unlink nodes in s
        for (auto i : s) f->unlinkNode(i);
        v.clear();
        return 0;
      }

      // otherwise copy s into v
      unsigned j = 0;
      for (auto i : s) v[j++] = i;
      v.resize(j);
    }
    assert(v.size() == 1);
    return v[0];
  }

  sorter mdd_sorter(f, heuristic);
  mdd_sorter.setLevel(f->getNumVariables());

  if (heuristic == spans) {
    //
    // (1) Find span for each node in v
    // (2) Find span for each pair of nodes
    // (3) Multiply the pair of nodes with the smallest combined span
    // (4) Remove the pair of nodes from v
    // (5) Add the result of mutiplication to v
    // (6) Repeat (1)-(5) till v has only one node
    //

    span_checker sc(f);
    std::set<node_handle> nodes;
    for (auto i : v) nodes.insert(i);
    for (auto i : nodes) f->linkNode(i);

    while (nodes.size() > 1) {
      node_pair pair = sc.getBestPair(nodes);
      node_handle product = mddMultiply->compute(pair.a, pair.b);
      sc.gc();
      unsigned product_span = sc.getSpan(product);
#if 0
      printf("Product of (%d:%d:%u) and (%d:%d:%u) : (%d:%d:%u)\n",
      pair.a, f->getNodeLevel(pair.a), sc.getSpan(pair.a),
      pair.b, f->getNodeLevel(pair.b), sc.getSpan(pair.b),
      product, f->getNodeLevel(product), sc.getSpan(product)
      );
#endif
      nodes.erase(pair.a);
      nodes.erase(pair.b);
      nodes.insert(product);
      f->unlinkNode(pair.a);
      f->unlinkNode(pair.b);
      if (0 == product_span) {
        for (auto i : nodes) f->unlinkNode(i);
        break;
      }
    }

    node_handle result = nodes.empty()? 0: *(nodes.begin());

    return result;

  } else {

    std::priority_queue<node_handle, std::vector<node_handle>, sorter> pq(mdd_sorter);
    printf("\n[");
    for (auto n : v) { 
      printf(" (%d)", n);
      pq.push(n);
    }
    printf(" ]\n");
    v.clear();

    printf("\n");
    while (pq.size() > 1 && pq.top() != 0) {
      node_handle first = pq.top(); pq.pop();
      node_handle second = pq.top(); pq.pop();
      printf("\t%lu (%d) (%d): ", 2+pq.size(), first, second);
      pq.push(mddMultiply->compute(first, second));
      printf("(%d)\n", pq.top());
      f->unlinkNode(first);
      f->unlinkNode(second);
    }

    assert(!pq.empty());
    node_handle result = pq.top(); pq.pop();
    while (!pq.empty()) { f->unlinkNode(pq.top()); pq.pop(); }

    return result;
  }
}

bool processArgs(int argc, const char** argv, forest::policies &p)
{
  lfile = 0;
  p.setPessimistic();
  N = -1;
  int i;
  use_batches = false;
  inner_batch = false;
  outer_batch = false;
  monolithic_batch = false;
  heuristic = node_count;
  for (i=1; i<argc; i++) {
    if (strcmp("-opt", argv[i])==0) {
      p.setOptimistic();
      continue;
    }
    if (strcmp("-pess", argv[i])==0) {
      p.setPessimistic();
      continue;
    }
    if (strcmp("-acc", argv[i])==0) {
      use_batches = false;
      continue;
    }
    if (strcmp("-acc_batch", argv[i])==0) {
      use_batches = true;
      heuristic = accumulate;
      continue;
    }
    if (strcmp("-folding", argv[i])==0) {
      use_batches = true;
      heuristic = folding;
      continue;
    }
    if (strcmp("-node_count", argv[i])==0) {
      use_batches = true;
      heuristic = node_count;
      continue;
    }
    if (strcmp("-cardinality", argv[i])==0) {
      use_batches = true;
      heuristic = cardinality;
      continue;
    }
    if (strcmp("-spans", argv[i])==0) {
      use_batches = true;
      heuristic = spans;
      continue;
    }
    if (strcmp("-n_ary", argv[i])==0) {
      use_batches = true;
      heuristic = n_ary_intersection;
      continue;
    }
    if (strcmp("-inner", argv[i])==0) {
      use_batches = true;
      inner_batch = true;
      continue;
    }
    if (strcmp("-outer", argv[i])==0) {
      use_batches = true;
      outer_batch = true;
      continue;
    }
    if (strcmp("-monolithic", argv[i])==0) {
      use_batches = true;
      monolithic_batch = true;
      continue;
    }
    if (strcmp("-l", argv[i])==0) {
      lfile = argv[i+1];
      i++;
      continue;
    }
    N = atoi(argv[i]);
    break;
  }
  if (N<1) return false;
  if (++i < argc) {
    outfile = fopen(argv[2], "w");
    if (0==outfile) {
      printf("Couldn't open %s for writing, no solutions will be written\n", argv[2]);
    }
  } else {
    outfile = 0;
  }
  if (use_batches && !inner_batch && !monolithic_batch) outer_batch = true;
  return true;
}

int usage(const char* who)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) name = ptr+1;
  }
  printf("Usage: %s <-opt|-pess> <-acc|-node_count|-cardinality> <-l lfile> N <outfile>\n\n",
      name);
  printf("\t           N:  board dimension\n");
  printf("\t        -opt:  Optimistic node deletion\n");
  printf("\t       -pess:  Pessimistic node deletion (default)\n");
  printf("\t        -acc:  Accumulate constraints in order (default)\n");
  printf("\t -node_count:  Accumulate constraints in batches, with node_count as heuristic\n");
  printf("\t-cardinality:  Accumulate constraints in batches, with cardinality as heuristic\n");
  printf("\t      -spans:  Accumulate constraints in batches, with spans as heuristic\n");
  printf("\t      -inner:  Apply heuristic when combining column constraints (1 batch/row)\n");
  printf("\t      -outer:  Apply heuristic when combining row constraints\n");
  printf("\t -monolithic:  Apply heuristic when combining all constraints using a single batch\n");
  printf("\t    -l lfile:  Write logging information to specified file\n");
  printf("\t   <outfile>:  if specified, we write all solutions to this file\n\n");
  return 1;
}

int main(int argc, const char** argv)
{
  initialize();
  forest::policies p(false);
  if (!processArgs(argc, argv, p)) return usage(argv[0]);
  printf("Using %s\n", getLibraryInfo(0));

  timer stopwatch;
  stopwatch.note_time();

  printf("\nDetermining queen covers for %dx%d chessboard.\n", N, N);

  expert_forest* f = dynamic_cast<expert_forest*>(buildQueenForest(p));
  
  std::ofstream log;
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
      snprintf(comment, 80, "Automatically generated by queen_cover (N=%d)", N);
      LOG->addComment(comment);
      f->setLogger(LOG, "forest name");
      LOG->newPhase(f, "Initializing");
    }
  }

  dd_edge num_queens(f);
  f->createEdge(int(0), num_queens);
  for (int i=0; i<N; i++) for (int j=0; j<N; j++) {
    dd_edge tmp(f);
    hasQueen(f, i, j, tmp);
    apply(PLUS, num_queens, tmp, num_queens);
  } // for i,j

  // "By-hand" compute tables for queen in col, diag
  qic = new dd_edge*[N];
  for (int i=0; i<N; i++) qic[i] = 0;
  qidp = new dd_edge*[2*N-1];
  qidm = new dd_edge*[2*N-1];
  for (int i=0; i<2*N-1; i++) {
    qidp[i] = qidm[i] = 0;
  }

  dd_edge solutions(f);

  long time_for_multiply = 0L;

  int q;
  for (q=1; q<=N; q++) {

    timer timer_for_q;
    timer_for_q.note_time();
    long time_for_q = 0L;

    char buffer[80];
    snprintf(buffer, 80, "Trying to cover with %d queens", q);
    printf("\n%s\n", buffer);
    if (LOG) LOG->newPhase(f, buffer);

    // Build all possible placements of q queens:
    dd_edge qqueens(f);
    f->createEdge(q, qqueens);
    apply(EQUAL, qqueens, num_queens, qqueens);
    solutions = qqueens;

    if (!use_batches) {
      printf("\tDefault accumulation of constraints\n");
      printf("\tAdding constraints by row\n\t\t");
      for (int i=0; i<N && solutions.getNode() != 0; i++) {
        printf("%d", N-i);
        fflush(stdout);
        dd_edge rowcov = qqueens;
        dd_edge qir(f);
        queenInRow(f, i, qir);
        for (int j=0; j<N; j++) {
          dd_edge col(f);
          queenInCol(f, j, col);
          dd_edge dgp(f);
          queenInDiagP(f, i+j, dgp);
          dd_edge dgm(f);
          queenInDiagM(f, i-j, dgm);
          // "OR" together the possible attackers for this square
          apply(MAXIMUM, col, dgp, col);
          apply(MAXIMUM, col, dgm, col);
          apply(MAXIMUM, col, qir, col);
          // "AND" with this row
          timer_for_q.note_time();
          apply(MULTIPLY, rowcov, col, rowcov);
          timer_for_q.note_time();
          time_for_q += timer_for_q.get_last_interval();
        } // for j
        printf(", ");
        fflush(stdout);
        // "AND" directly with set of covers
        timer_for_q.note_time();
        apply(MULTIPLY, solutions, rowcov, solutions);
        timer_for_q.note_time();
        time_for_q += timer_for_q.get_last_interval();
      } // for i
    } else {
      printf("\tAccumulation of constraints using heuristics\n");
      std::vector<MEDDLY::node_handle> constraints;
      if (outer_batch || monolithic_batch) {
        constraints.push_back(f->linkNode(solutions.getNode()));
      }
      printf("\tAdding constraints by row\n\t\t");
      for (int i=0; i<N; i++) {
        printf("%d", N-i);
        fflush(stdout);
        std::vector<MEDDLY::node_handle> row_constraints;
        dd_edge row_cov = qqueens;
        dd_edge qir(f);
        queenInRow(f, i, qir);
        for (int j=0; j<N; j++) {
          dd_edge col(f);
          queenInCol(f, j, col);
          dd_edge dgp(f);
          queenInDiagP(f, i+j, dgp);
          dd_edge dgm(f);
          queenInDiagM(f, i-j, dgm);
          // "OR" together the possible attackers for this square
          apply(MAXIMUM, col, dgp, col);
          apply(MAXIMUM, col, dgm, col);
          apply(MAXIMUM, col, qir, col);
          if (inner_batch || monolithic_batch) {
            row_constraints.push_back(f->linkNode(col.getNode()));
          } else {
            timer_for_q.note_time();
            apply(MULTIPLY, row_cov, col, row_cov);
            timer_for_q.note_time();
            time_for_q += timer_for_q.get_last_interval();
          }
        } // for j
        printf(", ");
        fflush(stdout);
        if (monolithic_batch) {
          for (auto n : row_constraints) constraints.push_back(n);
        } else {
          if (inner_batch) {
            row_constraints.push_back(f->linkNode(qqueens.getNode()));
            timer_for_q.note_time();
            row_cov.set(batchMultiply(f, heuristic, row_constraints));
            timer_for_q.note_time();
            time_for_q += timer_for_q.get_last_interval();
          }
          if (outer_batch) {
            constraints.push_back(f->linkNode(row_cov.getNode()));
          } else {
            // "AND" directly with set of covers
            timer_for_q.note_time();
            apply(MULTIPLY, solutions, row_cov, solutions);
            timer_for_q.note_time();
            time_for_q += timer_for_q.get_last_interval();
          }
        }
      } // for i

      if (monolithic_batch || outer_batch) {
        timer_for_q.note_time();
        printf("Computing solution without emptiness check.\n");
        solutions.set(batchMultiply(f, heuristic, constraints));
        timer_for_q.note_time();
        time_for_q += timer_for_q.get_last_interval();
      }
    }


    printf("\n");

    time_for_multiply += time_for_q;

    printf("\nq: %d, %lg seconds CPU time elapsed\n",
      q,
      time_for_q / 1000000.0);

    // Any solutions?
    if (0==solutions.getNode()) {
      printf("\tNo solutions\n");
      continue;
    } else {
      printf("\tSuccess\n");
      break;
    }
  } // for q

  // Cleanup
  for (int i=0; i<N; i++) delete qic[i];
  for (int i=0; i<2*N-1; i++) {
    delete qidp[i];
    delete qidm[i];
  }
  delete[] qic;
  delete[] qidp;
  delete[] qidm;

  printf("\n%lg seconds CPU time elapsed\n", 
    stopwatch.get_last_interval() / 1000000.0);
  printf("Forest stats:\n");
  FILE_output myout(stdout);
  expert_forest* ef = (expert_forest*)f;
  ef->reportStats(myout, "\t", 
    expert_forest::HUMAN_READABLE_MEMORY  |
    expert_forest::BASIC_STATS | expert_forest::EXTRA_STATS |
    expert_forest::STORAGE_STATS | expert_forest::HOLE_MANAGER_STATS | 
    expert_forest::HOLE_MANAGER_DETAILED
  );
  operation::showAllComputeTables(myout, 3);

  long c;
  apply(CARDINALITY, solutions, c);
  printf("\nFor a %dx%d chessboard, ", N, N);
  printf("there are %ld covers with %d queens\n\n", c, q);

  if (outfile) {
    printf("Writing solutions to file %s\n", argv[2]);
  
    fprintf(outfile, "%d # Board dimension\n\n", N);
    // show the solutions
    enumerator iter(solutions);
    long counter;
    for (counter = 1; iter; ++iter, ++counter) {
      fprintf(outfile, "solution %5ld:  ", counter);
      const int* minterm = iter.getAssignments();
      for (int i=0; i<N; i++) for (int j=0; j<N; j++) {
        if (minterm[ijmap(i,j)]) {
          fprintf(outfile, "(%2d, %2d) ", i+1, j+1);
        }
      }
      fprintf(outfile, "\n");
    } // for iter
    fprintf(outfile, "\n");
  }

  if (LOG) LOG->newPhase(f, "Cleanup");
  cleanup();
  delete LOG;

  printf("\n%lg seconds CPU time elapsed for multiply\n", 
    time_for_multiply / 1000000.0);

  return 0;
}

// ---------------------------------------------------------------------------
//
// class: sorter
//
// ---------------------------------------------------------------------------

sorter::sorter(forest* f, batch_heuristic h)
: heuristic(h), level(0)
{
  this->f = dynamic_cast<expert_forest*>(f);
  switch (heuristic) {
    case node_count:
      op_type = INTEGER;
      op = getOperation(NODE_COUNT, f, op_type);
      smaller_is_better = true;
      break;
    case cardinality:
      op_type = REAL;
      op = getOperation(CARDINALITY, f, op_type);
      smaller_is_better = true;
      break;
    default:
      break;
  }
}

// ---------------------------------------------------------------------------
//
// class: span_checker
//
// ---------------------------------------------------------------------------

span_checker::span_checker(forest* f)
{
  this->f = dynamic_cast<expert_forest*>(f);
}

void span_checker::gc() {
  std::map<node_handle, unsigned> new_span_cache;
  for (auto i : span_cache) {
    if (!f->isStale(i.first)) {
      new_span_cache.emplace_hint(new_span_cache.end(), i.first, i.second);
    }
  }
  span_cache.swap(new_span_cache);

  pair_span_cache.clear();
}


unsigned span_checker::getSpan(node_handle node) {
  // terminal cases
  int k = f->getNodeLevel(node);
  if (k == 0) return 0;

  // check in cache
  std::map<node_handle, unsigned>::iterator it = span_cache.find(node);
  if (it != span_cache.end()) return it->second;

  // compute span
  unsigned span = 0;
  unpacked_node* n = unpacked_node::newFromNode(f, node, false);
  for (int i = 0; i < n->getNNZs(); i++) {
    assert(n->d(i));
    int child_k = f->getNodeLevel(n->d(i));
    assert(k > child_k);
    unsigned child_span = getSpan(n->d(i)) + unsigned(k - child_k);
    if (child_span > span) span = child_span;
  }

  // save and return result
  unpacked_node::recycle(n);
  span_cache[node] = span;
  return span;
}

std::map<node_handle, unsigned> span_checker::getNode2SpanMap(std::set<node_handle>& s) {
  std::map<node_handle, unsigned> node2span;
  for (auto i : s) {
    node2span[i] = getSpan(i);
  }
  return node2span;
}

node_pair span_checker::getBestPair(std::set<node_handle>& s) {
  // Find the span for each node
  // Find the span for each pair of nodes
  // Find the best span and return it
  std::map<node_handle, unsigned> nodes2spans = getNode2SpanMap(s);
#if 0
  printf("Spans:\n");
  for (auto i : nodes2spans) {
    printf("\t%d: %d : %u\n", i.first, f->getNodeLevel(i.first), i.second);
  }
#endif
  std::map<node_handle, unsigned>::iterator iter1 = nodes2spans.begin();

  node_pair best_pair;
#ifdef USE_SMALLER_SPANS
  unsigned best_span = unsigned(f->getNumVariables());
#else
  unsigned best_span = 0;
#endif
  bool first_time = true;

  for ( ; iter1 != nodes2spans.end(); iter1++) {
    int level_a = f->getNodeLevel(iter1->first);
    unsigned span_a = iter1->second;
    for (std::map<node_handle, unsigned>::iterator iter2 = iter1;
        ++iter2 != nodes2spans.end(); ) {
      int level_b = f->getNodeLevel(iter2->first);
      unsigned span_b = iter2->second;
      unsigned span_a_b = 
        (level_a > level_b)
        ? MAX( span_a, span_b + unsigned(level_a - level_b) )
        : MAX( span_b, span_a + unsigned(level_b - level_a) )
        ;
#ifdef USE_SMALLER_SPANS
      if (span_a_b < best_span || first_time) {
#else
      if (span_a_b > best_span || first_time) {
#endif
        best_pair.a = iter1->first;
        best_pair.b = iter2->first;
        best_span = span_a_b;
        first_time = false;
      }
    }
  }

#if 0
  printf("Best pair (node:k:span): ( %d:%d:%u , %d:%d:%u ), span: %u\n",
  best_pair.a, f->getNodeLevel(best_pair.a), nodes2spans[best_pair.a],
  best_pair.b, f->getNodeLevel(best_pair.b), nodes2spans[best_pair.b],
  best_span);
#endif
  return best_pair;
}


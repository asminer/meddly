
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

#include <cstdlib>
#include <iostream>
#include <cassert>
#include "../src/meddly.h"

using namespace MEDDLY;

domain* initDomain() {
  MEDDLY::initialize();
  int bounds[] = {-1, -1, -1};
  const int n_vars = 3;
  domain* d = createDomainBottomUp(bounds, n_vars);
  assert(d);
  return d;
}

forest* createForest(domain *d, bool relation) {
  forest* mxd = d->createForest(relation, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
  assert(mxd);
  return mxd;
}

void createEdge(forest* mdd, int var, const std::vector<int>& indices, dd_edge& result) {
  // one minterm per index
  const int nVars = mdd->useDomain()->getNumVariables();
  int** minterms = new int*[indices.size()];
  for (unsigned i = 0; i < indices.size(); i++) {
    minterms[i] = new int[nVars+1];
    for (int j = 0; j <= nVars; j++) {
      minterms[i][j] = DONT_CARE;
    }
    minterms[i][var] = indices[i];
  }
  mdd->createEdge(minterms, indices.size(), result);
  for (unsigned i = 0; i < indices.size(); i++) delete [] minterms[i];
  delete [] minterms;
}

void createEdge(forest* mdd, int var, int value, dd_edge& result) {
  // one minterm per index
  const int nVars = mdd->useDomain()->getNumVariables();
  int** minterms = new int*[value];
  for (int i = 0; i < value; i++) {
    minterms[i] = new int[nVars+1];
    for (int j = 0; j <= nVars; j++) {
      minterms[i][j] = DONT_CARE;
    }
    minterms[i][var] = i;
  }
  mdd->createEdge(minterms, value, result);
  for (int i = 0; i < value; i++) delete [] minterms[i];
  delete [] minterms;
}

// disabling expressions for transitions:
// {
// (k_a:v_a, k_b:v_b, ...), // t1: k_a < v_a | k_b < v_b | ...
// (k_c:v_c, k_d:v_d, ...), // t2: k_c < v_c | k_d < v_d | ...
// ...
// }
//
// potential deadlock states = conjunction of transition disabling expressions

struct int_pair {
  int level;
  int value;
  int_pair(int k, int v) : level(k), value(v) {}
};

void buidPotentialDeadlockStates(forest* mdd,
    std::vector<std::vector<int_pair>>& disabling_expressions,
    dd_edge& result) {
  MEDDLY_DCASSERT(mdd && !mdd->isForRelations());

  if (disabling_expressions.size() == 0) return;

  std::vector<dd_edge> disabling_ddedges;

  for (auto i : disabling_expressions) {
    dd_edge i_union(mdd);
    for (auto j : i) {
      dd_edge expr(mdd);
      createEdge(mdd, j.level, j.value, expr);
      i_union += expr;
    }
    disabling_ddedges.push_back(i_union);
  }

  ostream_output s(std::cout);
  MEDDLY_DCASSERT(disabling_ddedges.size() > 0);
  result = disabling_ddedges[0];
  for (auto i : disabling_ddedges) {
    std::cout << "\nDisabling dd_edge:\n";
    i.showGraph(s);
    result *= i;
  }
}


void testStateSetForest(domain *d) {
  forest* mdd = createForest(d, false);

  std::vector<int> indices;
  for (int i = 0; i < 3; i++) indices.push_back(i);
  dd_edge level_1_le_3(mdd);
  dd_edge level_2_le_3(mdd);
  dd_edge level_3_le_3(mdd);
  createEdge(mdd, 1, indices, level_1_le_3);
  createEdge(mdd, 2, indices, level_2_le_3);
  createEdge(mdd, 3, indices, level_3_le_3);

  dd_edge result = level_1_le_3;
  result += level_2_le_3;
  result += level_3_le_3;

  ostream_output s(std::cout);
  level_1_le_3.showGraph(s);
  level_2_le_3.showGraph(s);
  level_3_le_3.showGraph(s);
  result.showGraph(s);

  // building potential deadlock states
  // disabling expressions per transition:
  // t1 : k1 < 3 or k2 < 4
  // t2 : k1 < 2 or k3 < 2
  // t3 : k2 < 1
  std::vector<std::vector<int_pair>> disabling_expressions;
  std::vector<int_pair> t1;
  t1.emplace_back(1, 3);
  t1.emplace_back(2, 4);
  std::vector<int_pair> t2;
  t2.emplace_back(1, 2);
  t2.emplace_back(3, 2);
  std::vector<int_pair> t3;
  t3.emplace_back(2, 1);
  disabling_expressions.push_back(t1);
  disabling_expressions.push_back(t2);
  disabling_expressions.push_back(t3);

  buidPotentialDeadlockStates(mdd, disabling_expressions, result);
  std::cout << "\nPotential deadlock states dd_edge:\n";
  result.showGraph(s);
}


void testRelationForest(domain *d) {
  const int n_vars = d->getNumVariables();
  forest* mxd = createForest(d, true);

  int* unp_minterm_1 = new int[n_vars+1];
  int* unp_minterm_2 = new int[n_vars+1];
  int* p_minterm_1 = new int[n_vars+1];
  int* p_minterm_2 = new int[n_vars+1];
  int* unp_minterms[] = {unp_minterm_1, unp_minterm_2};
  int* p_minterms[] = {p_minterm_1, p_minterm_2};

  for (int i = 0; i <= n_vars; i++) {
    unp_minterm_1[i] = DONT_CARE;
    p_minterm_1[i] = DONT_CARE;
    unp_minterm_2[i] = DONT_CARE;
    p_minterm_2[i] = DONT_CARE;
  }
  unp_minterm_1[2] = 0;
  p_minterm_1[2] = 1;
  unp_minterm_2[2] = 1;
  p_minterm_2[2] = 2;

  dd_edge with_dont_change(mxd);
  dd_edge only_dont_cares(mxd);

  mxd->createEdge((int**)(&unp_minterm_1), (int**)(&p_minterm_1), 1, with_dont_change);
  mxd->createEdge((int**)(&unp_minterm_2), (int**)(&p_minterm_2), 1, only_dont_cares);

  dd_edge result = with_dont_change;
  result += only_dont_cares;

  ostream_output s(std::cout);
  with_dont_change.showGraph(s);
  only_dont_cares.showGraph(s);
  result.showGraph(s);

  mxd->createEdge((int**)(&unp_minterms[0]), (int**)(&p_minterms[0]), 2, result);
  result.showGraph(s);

  delete [] unp_minterm_1;
  delete [] unp_minterm_2;
  delete [] p_minterm_1;
  delete [] p_minterm_2;
}


int main(int argc, char* argv[]) {
  domain* d = initDomain();
  testRelationForest(d);
  testStateSetForest(d);
  MEDDLY::cleanup();
}

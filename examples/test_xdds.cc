
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
#include "meddly.h"

using namespace MEDDLY;

int main(int argc, char* argv[]) {
  MEDDLY::initialize();
  int bounds[] = {-1, -1, -1};
  const int n_vars = 3;
  domain* d = createDomainBottomUp(bounds, n_vars);
  assert(d);
  forest* mxd = d->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL);
  assert(mxd);

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
  with_dont_change.show(s, 2);
  only_dont_cares.show(s, 2);
  result.show(s, 2);

  mxd->createEdge((int**)(&unp_minterms[0]), (int**)(&p_minterms[0]), 2, result);
  result.show(s, 2);

  delete [] unp_minterm_1;
  delete [] unp_minterm_2;
  delete [] p_minterm_1;
  delete [] p_minterm_2;

  MEDDLY::cleanup();
}

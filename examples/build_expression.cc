
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



/**
 * build_expression.cc
 *
 * Demonstrates use of
 *   createEdgeForVar(int, bool, dd_edge&), and
 *   createEdgeForVar(int, bool, int*, dd_edge&)
 *   createEdgeForVar(int, bool, float*, dd_edge&)
 */

#include <iostream>

#include "../include/mddlib.h"
#include "../src/operation_ext.h"
#include "../src/timer.h"

#define RELATION 0
#define USE_REALS 1

#if USE_REALS
  typedef float element_type;
#else
  typedef int element_type;
#endif


// Builds the expression y1 + 2*y2 + 3*y3
// Assumes states has atleast 3 levels, each of size 2.
// Assumed variable order: y3 y2 y1 TERMINALS.
// Uses createEdgeForVar(int, bool, dd_edge&)
dd_edge buildExpression(forest* states);

// Similar to buildExpression, but uses a array to define the terminals.
// Uses createEdgeForVar(int, bool, int*/float*, dd_edge&)
dd_edge buildExpressionWithTerms(forest* states);

int main(int argc, char *argv[])
{
  const int nVariables = 4;
  int variableBound = 2;

  // initialize the variable bounds array to provide to the domain
  int bounds[nVariables];
  for (int i = 0; i < nVariables; ++i)
    bounds[i] = variableBound;

  // Create a domain
  domain *d = MDDLIB_createDomain();
  assert(d != 0);
  assert(domain::SUCCESS == d->createVariablesBottomUp(bounds, nVariables));

  // Create an MXD forest in this domain (to store states)

#if RELATION
  bool relation = true;
#else
  bool relation = false;
#endif

#if USE_REALS
  forest::range_type range = forest::REAL;
#else
  forest::range_type range = forest::INTEGER;
#endif

  forest* states =
    d->createForest(relation, range, forest::MULTI_TERMINAL);
  assert(states != 0);

  dd_edge expr = buildExpression(states);
  expr.show(stdout, 2);

  expr = buildExpressionWithTerms(states);
  expr.show(stdout, 2);

  // Cleanup; in this case simply delete the domain
  delete d;

  return 0;
}

dd_edge buildExpression(forest* states)
{
  // Building expression y1 + 2*y2 + 3*y3

  // ---- Building y1 ----

  // y1
  dd_edge y1(states);
  states->createEdgeForVar(1, false, y1);

  // ---- Building 2*y2 ----

  // y2
  dd_edge y2(states);
  states->createEdgeForVar(2, false, y2);

  // constant 2
  dd_edge cons2(states);
#if USE_REALS
  states->createEdge(2.0f, cons2);
#else
  states->createEdge(2, cons2);
#endif

  // 2y2
  y2 *= cons2;

  // ---- Building 3*y3 ----

  // y3
  dd_edge y3(states);
  states->createEdgeForVar(3, false, y3);

  // constant 3
  dd_edge cons3(states);
#if USE_REALS
  states->createEdge(3.0f, cons3);
#else
  states->createEdge(3, cons3);
#endif
  
  // 3y3
  y3 *= cons3;

  // ---- Building y1 + 2*y2 + 3*y3 ----

  y1 = y1 + y2 + y3;
  
  return y1;
}

dd_edge buildExpressionWithTerms(forest* states)
{
  // Building expression y1 + 2*y2 + 3*y3

#if USE_REALS
  float terms[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
#else
  int terms[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
#endif

  // ---- Building y1 ----

  // y1
  dd_edge y1(states);
  assert(forest::SUCCESS == states->createEdgeForVar(1, false, terms, y1));

  // ---- Building 2*y2 ----

  // y2
  dd_edge y2(states);
  assert(forest::SUCCESS == states->createEdgeForVar(2, false, terms, y2));

  // constant 2
  dd_edge cons2(states);
#if USE_REALS
  states->createEdge(2.0f, cons2);
#else
  states->createEdge(2, cons2);
#endif

  // 2y2
  y2 *= cons2;

  // ---- Building 3*y3 ----

  // y3
  dd_edge y3(states);
  assert(forest::SUCCESS == states->createEdgeForVar(3, false, terms, y3));

  // constant 3
  dd_edge cons3(states);
#if USE_REALS
  states->createEdge(3.0f, cons3);
#else
  states->createEdge(3, cons3);
#endif
  
  // 3y3
  y3 *= cons3;

  // ---- Building y1 + 2*y2 + 3*y3 ----

  y1 = y1 + y2 + y3;

  return y1;
}



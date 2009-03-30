
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



// TODO: Testing.
// TODO: Finish the advanced functions in expert_domain.


#include "../include/mddexpert.h"
#include "../src/mdds_ext.h"



// ----------------- domain ------------------------

domain::domain() {}
domain::~domain() {}

const char* domain::getErrorCodeName(domain::error e)
{
  switch (e) {
    case SUCCESS:
        return "Operation returned successfully";
    case NOT_IMPLEMENTED:
        return "Operation not implemented";
    case INSUFFICIENT_MEMORY:
        return "Operation failed -- lack of memory";
    case DOMAIN_NOT_EMPTY:
        return "Operation failed -- needs empty domain";
    case INVALID_VARIABLE:
        return "Operation failed -- invalid variable handle";
    case INVALID_BOUND:
        return "Operation falied -- variable bound out of range";
    default:
        return "Unknown error code";
  }
}

domain* MDDLIB_createDomain()
{
  return new expert_domain();
}


// ----------------- expert_domain ------------------------

expert_domain::expert_domain()
: nForests(0), forests(0), szForests(0),
  nVars(0), levelBounds(0), allocatedLevels(0), topLevel(0),
  levelsToHeightsMap(0), heightsToLevelsMap(0)
{}


expert_domain::~expert_domain()
{
  // delete registered forests
  for (int i = 0; i < szForests; ++i) {
    if (forests[i] != 0) {
      delete forests[i];
      assert(forests[i] == 0);
    }
  }
  
  // free arrays
  free(levelBounds);
  free(levelsToHeightsMap);
  free(heightsToLevelsMap);
}


domain::error expert_domain::createVariablesBottomUp(const int* bounds, int N)
{
  // domain must be empty -- no variables defined so far
  DCASSERT(nForests == 0);
  DCASSERT(nVars == 0);
  if (nForests != 0 || nVars != 0)
    return domain::DOMAIN_NOT_EMPTY;

  // N is the number of variables (0: TERMINALS, 1: bottom , N: top)
  allocatedLevels = 2;
  while (allocatedLevels < (N + 1))
    allocatedLevels *= 2;
  nVars = N;
  topLevel = N;

  levelBounds = (int *) malloc(allocatedLevels * sizeof(int));
  if (levelBounds == 0) {
    return domain::INSUFFICIENT_MEMORY;
  }

  levelsToHeightsMap = (int *) malloc(allocatedLevels * sizeof(int));
  if (levelsToHeightsMap == 0) {
    free(levelBounds);
    return domain::INSUFFICIENT_MEMORY;
  }
  
  heightsToLevelsMap = (int *) malloc(allocatedLevels * sizeof(int));
  if (heightsToLevelsMap == 0) {
    free(levelBounds);
    free(levelsToHeightsMap);
    return domain::INSUFFICIENT_MEMORY;
  }

  levelBounds[0] = 0;
  levelsToHeightsMap[0] = 0;
  heightsToLevelsMap[0] = 0;
  for (int i = 1; i < N + 1; ++i) {
    levelBounds[i] = bounds[i-1];
    levelsToHeightsMap[i] = i;
    heightsToLevelsMap[i] = i;
  }
  // identify invalid levels with -1
  for (int i = N + 1; i < allocatedLevels; ++i) {
    levelBounds[i] = -1;
    levelsToHeightsMap[i] = -1;
    heightsToLevelsMap[i] = -1;
  }

  return domain::SUCCESS;
}


domain::error expert_domain::createVariablesTopDown(const int* bounds, int N)
{
  // use createVariablesBottomUp, and then fix the levels to heights mapping
  domain::error err = createVariablesBottomUp(bounds, N);
  if (err != SUCCESS)
    return err;

  // modify levelsToHeightsMap and heightsToLevelsMap such that
  // height(level1) = N, height(level2) = N - 1, ..., height(levelN) = 1
  for (int i = 1; i < N + 1; ++i) {
    levelsToHeightsMap[i] = N + 1 - i;
    heightsToLevelsMap[N + 1 - i] = i;
  }
  // you could do that above with N/2 swaps but that is less readable
  // and this method is not a significant factor in library performance.

  return domain::SUCCESS;
}


void expert_domain::showInfo(FILE* strm)
{
  // list variables handles, their bounds and heights.
  fprintf(strm, "Domain info:\n");
  fprintf(strm, "  #variables: %d\n", nVars);
  fprintf(strm, "  Variables listed in height-order (ascending):\n");
  fprintf(strm, "    height\\tthandle\t\tbound\n");
  for (int i = 1; i < nVars + 1; ++i) {
    fprintf(strm, "    %d\t\t%d\t\t%d\n",
            i, heightsToLevelsMap[i], levelBounds[heightsToLevelsMap[i]]);
  }
  // call showNodes for each of the forests in this domain.
  for (int i = 0; i < szForests; i++) {
    if (forests[i] != 0)
      forests[i]->showInfo(strm);
  }
}


forest* expert_domain::createForest(bool rel, forest::range_type t,
    forest::edge_labeling e)
{
  expert_forest* f = 0;

  if (rel) {
    if(e == forest::MULTI_TERMINAL) {
      if (t == forest::BOOLEAN) {
        f = new mxd_node_manager(this);
      } else if (t == forest::INTEGER || t == forest::REAL) {
        f = new mtmxd_node_manager(this, t);
      }
    }
  } else {
    if (e == forest::MULTI_TERMINAL) {
      if (t == forest::BOOLEAN) {
        f = new mdd_node_manager(this);
      } else if (t == forest::INTEGER || t == forest::REAL) {
        f = new mtmdd_node_manager(this, t);
      }
    } else if (e == forest::EVPLUS) {
      if (t == forest::INTEGER || t == forest::REAL) {
        f = new evmdd_node_manager(this, t);
      }
    }
  }

  if (f != 0) {
    // Add it to the list of forests
    assert (nForests <= szForests);
    if (szForests == 0) {
      // initialize forests[]
      szForests = 4;
      forests = (expert_forest **) malloc(szForests * sizeof(expert_forest *));
      assert(forests != 0);
      memset(forests, 0, szForests * sizeof(expert_forest *));
      assert(nForests == 0);
      forests[nForests] = f;
    } else if (nForests == szForests) {
      // expand forests[]
      expert_forest** temp = (expert_forest **) realloc(forests,
          szForests * 2 * sizeof (expert_forest *));
      assert(temp != 0);
      forests = temp;
      memset(forests + szForests, 0, szForests * sizeof(expert_forest*));
      szForests *= 2;
      forests[nForests] = f;
    } else {
      // find the first hole in forests[] -- this array should be small,
      // modifications will be rare and therefore O(n^2) complexity is
      // still not significant.
      // TODO: modify scheme to the one that dd_edges and forests share.
      int i = 0;
      for (i = 0; i < szForests && forests[i] != 0; ++i)
        ;
      assert (i < szForests);
      forests[i] = f;
    }
    nForests++;
  }
  
  return f;
}


void expert_domain::unlinkForest(expert_forest* f)
{
  // find forest
  int i = 0;
  for (i = 0; i < szForests; ++i) {
    if (forests[i] == f) {
      // unlink forest
      forests[i] = 0;
      return;
    }
  }
  assert(false);
}


// TODO: finish the rest

domain::error expert_domain::createVariable(int below, int &vh)
{
  return domain::NOT_IMPLEMENTED;
}


domain::error expert_domain::destroyVariable(int vh)
{
  return domain::NOT_IMPLEMENTED;
}


domain::error expert_domain::swapOrderOfVariables(int vh1, int vh2)
{
  return domain::NOT_IMPLEMENTED;
}

int expert_domain::findVariableBound(int vh) const
{
  // TODO: not implemented
  printf("expert_domain::findVariableBound() not implemented;");
  printf(" use getVariableBound().\n");
  return getVariableBound(vh);
}

domain::error expert_domain::enlargeVariableBound(int vh, int b)
{
  // delete registered forests
  for (int i = 0; i < szForests; ++i)
  {
    if (forests[i] != 0) {
      if (forests[i]->getReductionRule() != forest::QUASI_REDUCED)
        return domain::NOT_IMPLEMENTED;
    }
  }

  if (getVariableBound(vh) == -1) return domain::NOT_IMPLEMENTED;

  if (levelBounds[vh] < b) levelBounds[vh] = b;
  return domain::SUCCESS;
}

domain::error expert_domain::shrinkVariableBound(int vh, int b, bool force)
{
  return domain::NOT_IMPLEMENTED;
}



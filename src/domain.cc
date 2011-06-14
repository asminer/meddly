
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


#include "defines.h"
#if 0
#include "forests/mdds_ext.h"
#else
#include "forests/mtmdd.h"
#include "forests/mtmxd.h"
#include "forests/evmdd.h"
#endif

// ----------------------------------------------------------------------
// varaiable
// ----------------------------------------------------------------------

MEDDLY::variable::variable(int b, char* n)
{
  un_bound = b;
  pr_bound = b;
  name = n;
}

MEDDLY::variable::~variable()
{
  delete[] name;
}

// ----------------------------------------------------------------------
// expert_varaiable
// ----------------------------------------------------------------------

MEDDLY::expert_variable::expert_variable(int b, char* n)
 : variable(b, n)
{
  domlist = 0;
  dl_alloc = 0;
  dl_used = 0;
}

MEDDLY::expert_variable::~expert_variable()
{
  free(domlist);
}

void MEDDLY::expert_variable::addToList(domain* d)
{
  if (dl_used >= dl_alloc) {
    int ns = dl_alloc+8;
    domain** dl = (domain**) realloc(domlist, ns * sizeof(void*));
    if (0==dl) throw error(error::INSUFFICIENT_MEMORY);
    dl_alloc = ns;
    domlist = dl;
  }
  domlist[dl_used] = d;
  dl_used++;
}

void MEDDLY::expert_variable::removeFromList(const domain* d)
{
  int find;
  for (find=0; find<dl_used; find++) {
    if (d == domlist[find]) break;
  }
  if (find >= dl_used) return;  // not found; should we throw something?
  domlist[find] = domlist[dl_used-1];
  dl_used--;
}

void MEDDLY::expert_variable::enlargeBound(bool prime, int b)
{
  throw error(error::NOT_IMPLEMENTED);
}

void MEDDLY::expert_variable::shrinkBound(int b, bool force)
{
  throw error(error::NOT_IMPLEMENTED);
}

// ----------------------------------------------------------------------
// domain
// ----------------------------------------------------------------------


MEDDLY::domain::domain(variable** v, int N) 
{
  vars = v;
  nVars = N;
  for (int i=0; i<N; i++) {
    ((expert_variable*)vars[i])->addToList(this);
  }
}

MEDDLY::domain::~domain() 
{
  for (int i=0; i<nVars; i++) {
    ((expert_variable*)vars[i])->removeFromList(this);
  }
  free(vars);
}

// ----------------------------------------------------------------------
// expert_domain
// ----------------------------------------------------------------------

MEDDLY::expert_domain::expert_domain(variable** x, int n)
: domain(x, n)
{
  nForests = 0;
  forests = 0;
  szForests = 0;
}


MEDDLY::expert_domain::~expert_domain()
{
  // delete registered forests
  for (int i = 0; i < szForests; ++i) {
    if (forests[i] != 0) {
      delete forests[i];
      DCASSERT(forests[i] == 0);
    }
  }
}


void MEDDLY::expert_domain::createVariablesBottomUp(const int* bounds, int N)
{
  // domain must be empty -- no variables defined so far
  if (nForests != 0 || nVars != 0)
    throw error(error::DOMAIN_NOT_EMPTY);

  vars = (variable**) malloc((1+N) * sizeof(void*));
  if (0==vars) throw error(error::INSUFFICIENT_MEMORY);
  nVars = N;

  vars[0] = 0;
  for (int i=1; i<=N; i++) {
    vars[i] = MEDDLY::createVariable(bounds[i-1], 0);
  }
}


void MEDDLY::expert_domain::createVariablesTopDown(const int* bounds, int N)
{
  // domain must be empty -- no variables defined so far
  if (nForests != 0 || nVars != 0)
    throw error(error::DOMAIN_NOT_EMPTY);

  vars = (variable**) malloc((1+N) * sizeof(void*));
  if (0==vars) throw error(error::INSUFFICIENT_MEMORY);
  nVars = N;

  vars[0] = 0;
  for (int i=N; i; i--) {
    vars[N-i+1] = MEDDLY::createVariable(bounds[i], 0);
  }
}


void MEDDLY::expert_domain::showInfo(FILE* strm)
{
  // list variables handles, their bounds and heights.
  fprintf(strm, "Domain info:\n");
  fprintf(strm, "  #variables: %d\n", nVars);
  fprintf(strm, "  Variables listed in height-order (ascending):\n");
  fprintf(strm, "    height\t\tname\t\tbound\t\tprime-bound\n");
  for (int i = 1; i < nVars + 1; ++i) {
    const char* name = vars[i]->getName();
    if (0==name) name = "null";
    fprintf(strm, "    %d\t\t%s\t\t%d\t\t%d\n",
            i, name, vars[i]->getBound(0), vars[i]->getBound(1));
  }

  // call showNodes for each of the forests in this domain.
  for (int i = 0; i < szForests; i++) {
    if (forests[i] != 0)
      forests[i]->showInfo(strm, 2);
  }
}


MEDDLY::forest* MEDDLY::expert_domain::createForest(bool rel, forest::range_type t,
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
      f = new evplusmdd_node_manager(this);
    } else if (e == forest::EVTIMES) {
      f = new evtimesmdd_node_manager(this);
    }
  }

  if (f != 0) {
    // Add it to the list of forests
    DCASSERT (nForests <= szForests);
#if 0
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
#else
    // Previous strategy of re-using forest handles can unexpected errors
    // (user node handles, cached entries, etc.)
    // Therefore, the new scheme does not recycling forest handles.
    // Considering that will be relatively few forest handles created,
    // this should not be an issue.

    if (nForests == szForests) {
      // expand forests[]
      int newSize = szForests == 0? 4: szForests * 2;
      expert_forest** temp = (expert_forest **) realloc(forests,
          newSize * sizeof (expert_forest *));
      if (0 == temp) outOfMemory();
      forests = temp;
      memset(forests + szForests, 0,
          (newSize - szForests) * sizeof(expert_forest*));
      szForests = newSize;
    }
    forests[nForests++] = f;
#endif
  }
  
  return f;
}


void MEDDLY::expert_domain::unlinkForest(expert_forest* f)
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




// TODO: not implemented
void MEDDLY::expert_domain::swapOrderOfVariables(int vh1, int vh2)
{
  throw error(error::NOT_IMPLEMENTED);
}

// TODO: not implemented
int MEDDLY::expert_domain::findVariableBound(int vh) const
{
  printf("expert_domain::findVariableBound() not implemented;");
  printf(" use getVariableBound().\n");
  return getVariableBound(vh, false);
}

/*

void MEDDLY::expert_domain::enlargeVariableBound(int vh, bool prime, int b)
{
  // if !prime, expand both prime and unprime
  // else expand only prime

  if (getVariableBound(vh, false) == -1) throw error(error::NOT_IMPLEMENTED);

  if (prime) {
    if (pLevelBounds[vh] < b) pLevelBounds[vh] = b;
  } else {
    if (levelBounds[vh] < b) levelBounds[vh] = b;
    if (pLevelBounds[vh] < b) pLevelBounds[vh] = b;
  }
}

// TODO: not implemented
void MEDDLY::expert_domain::shrinkVariableBound(int vh, int b, bool force)
{
  throw error(error::NOT_IMPLEMENTED);
}

*/

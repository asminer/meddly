
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
#include "forests/mtmddbool.h"
#include "forests/mtmddint.h"
#include "forests/mtmddreal.h"

#include "forests/mtmxdbool.h"
#include "forests/mtmxdint.h"
#include "forests/mtmxdreal.h"

#include "forests/evmdd.h"
#endif

// #define DEBUG_CLEANUP

namespace MEDDLY {
  extern settings meddlySettings;
}

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
  // if that was the last domain...
  if (0==dl_used) delete this;  
}

void MEDDLY::expert_variable::enlargeBound(bool prime, int b)
{
  if (prime) {
    if (pr_bound < b) pr_bound = b;
  } else {
    if (un_bound < b) un_bound = b;
    if (pr_bound < b) pr_bound = b;
  }
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
  for (int i=1; i<N; i++) {
    ((expert_variable*)vars[i])->addToList(this);
  }
  forests = 0;
  szForests = 0;

  //
  // Add myself to the master list
  //
  if (-1 == free_list) {
    expandDomList();
  }
  MEDDLY_DCASSERT(free_list != -1);
  my_index = free_list;
  free_list = dom_free[free_list];
  dom_list[my_index] = this;
  dom_free[my_index] = -1;

#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating domain #%d\n", my_index);
#endif
}

MEDDLY::domain::~domain() 
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Deleting domain #%d\n", my_index);
#endif

  for (int i=1; i<nVars; i++) {
    ((expert_variable*)vars[i])->removeFromList(this);
  }
  free(vars);

  for (int i=0; i<szForests; i++) {
    delete forests[i];
    MEDDLY_DCASSERT(0==forests[i]);
  }
  free(forests);

  //
  // Remove myself from the master list
  //
  dom_list[my_index] = 0;
  dom_free[my_index] = free_list;
  free_list = my_index;
}

void MEDDLY::domain::expandDomList()
{
  int ndls = dom_list_size + 16;
  domain** tmp_dl = (domain**) realloc(dom_list, ndls * sizeof(domain*));
  if (0==tmp_dl) throw error(error::INSUFFICIENT_MEMORY);
  dom_list = tmp_dl;
  int* tmp_df = (int*) realloc(dom_free, ndls * sizeof(int));
  if (0==tmp_df) throw error(error::INSUFFICIENT_MEMORY);
  dom_free = tmp_df;
  for (int i=dom_list_size; i<ndls; i++) {
    dom_list[i] = 0;
    dom_free[i] = i+1;
  }
  dom_free[ndls-1] = -1;
  free_list = dom_list_size;
  dom_list_size = ndls;
}

void MEDDLY::domain::deleteDomList()
{
  for (int i=0; i<dom_list_size; i++) {
    delete dom_list[i];
  }
  free(dom_list);
  free(dom_free);
  dom_list_size = 0;
  dom_list = 0;
  dom_free = 0;
  free_list = -1;
}

MEDDLY::forest* MEDDLY::domain::createForest(bool rel, forest::range_type t, 
    forest::edge_labeling e, const forest::policies &p)
{
  int slot = findEmptyForestSlot();

  expert_forest* f = 0;

  switch (e) {
    case forest::MULTI_TERMINAL:
        switch (t) {
            case forest::BOOLEAN:
                if (rel)  f = new mt_mxd_bool(slot, this, p);
                else      f = new mt_mdd_bool(slot, this, p);
                break;

            case forest::INTEGER:
                if (rel)  f = new mt_mxd_int(slot, this, p);
                else      f = new mt_mdd_int(slot, this, p);
                break;

            case forest::REAL:
                if (rel)  f = new mt_mxd_real(slot, this, p);
                else      f = new mt_mdd_real(slot, this, p);
                break;

            default:
                throw error(error::TYPE_MISMATCH);
        }; // range type switch
        break;

    case forest::EVPLUS:
      f = new evp_mdd_int(slot, this, p);
      break;

    case forest::EVTIMES:
      f = new evt_mdd_real(slot, this, p);
      break;

    default:
      throw error(error::TYPE_MISMATCH);
  } // edge label switch

  MEDDLY_DCASSERT(f);
  forests[slot] = f;
  return f;
}

MEDDLY::forest* 
MEDDLY::domain
::createForest(bool rel, forest::range_type t, forest::edge_labeling e)
{
  return createForest(rel, t, e, 
    rel ? meddlySettings.mxdDefaults : meddlySettings.mddDefaults);
}

void MEDDLY::domain::showInfo(FILE* strm)
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

void MEDDLY::domain::unlinkForest(forest* f, int slot)
{
  if (forests[slot] != f)
    throw error(error::MISCELLANEOUS);
  forests[slot] = 0;
}

int MEDDLY::domain::findEmptyForestSlot()
{
  for (int slot=0; slot<szForests; slot++) {
    if (0==forests[slot]) return slot;
  }
  // need to expand
  int newSize;
  if (szForests) {
    if (szForests > 16) newSize = szForests + 16;
    else                newSize = szForests * 2;
  } else {
    newSize = 4;
  }
  forest** temp = (forest **) realloc(
    forests, newSize * sizeof (expert_forest *)
  );
  if (0 == temp) throw error(error::INSUFFICIENT_MEMORY);
  forests = temp;
  memset(forests + szForests, 0,
      (newSize - szForests) * sizeof(expert_forest*));
  int slot = szForests;
  szForests = newSize;
  return slot;
}

void MEDDLY::domain::markForDeletion()
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Marking domain #%d for deletion\n", my_index);
#endif
  for (int slot=0; slot<szForests; slot++) 
    if (forests[slot]) forests[slot]->markForDeletion();
}



// ----------------------------------------------------------------------
// expert_domain
// ----------------------------------------------------------------------

MEDDLY::expert_domain::expert_domain(variable** x, int n)
: domain(x, n)
{
}


MEDDLY::expert_domain::~expert_domain()
{
}


void MEDDLY::expert_domain::createVariablesBottomUp(const int* bounds, int N)
{
  // domain must be empty -- no variables defined so far
  if (hasForests() || nVars != 0)
    throw error(error::DOMAIN_NOT_EMPTY);

  vars = (variable**) malloc((1+N) * sizeof(void*));
  if (0==vars) throw error(error::INSUFFICIENT_MEMORY);
  nVars = N;

  vars[0] = 0;
  for (int i=1; i<=N; i++) {
    vars[i] = MEDDLY::createVariable(bounds[i-1], 0);
    ((expert_variable*)vars[i])->addToList(this);
  }
}


void MEDDLY::expert_domain::createVariablesTopDown(const int* bounds, int N)
{
  // domain must be empty -- no variables defined so far
  if (hasForests() || nVars != 0)
    throw error(error::DOMAIN_NOT_EMPTY);

  vars = (variable**) malloc((1+N) * sizeof(void*));
  if (0==vars) throw error(error::INSUFFICIENT_MEMORY);
  nVars = N;

  vars[0] = 0;
  for (int i=N; i; i--) {
    vars[N-i+1] = MEDDLY::createVariable(bounds[i], 0);
    ((expert_variable*)vars[i])->addToList(this);
  }
}

void MEDDLY::expert_domain::insertVariableAboveLevel(int lev, variable* v)
{
  throw error(error::NOT_IMPLEMENTED);
}

void MEDDLY::expert_domain::removeVariableAtLevel(int lev)
{
  throw error(error::NOT_IMPLEMENTED);
}

int MEDDLY::expert_domain::findLevelOfVariable(const variable *v) const
{
  // TBD: more efficient implementation based on binary search?
  int i;
  for (i=nVars; i; i--) {
    if (vars[i] == v) break;
  }
  return i;
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


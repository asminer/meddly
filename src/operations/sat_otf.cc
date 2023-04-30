
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

#include "../defines.h"
#include "sat_otf.h"
#include <typeinfo> // for "bad_cast" exception
#include <set>

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../oper_special.h"
#include "../opname_satur.h"

namespace MEDDLY {
  class otfsat_by_events_opname;
  class otfsat_by_events_op;

  class common_otf_dfs_by_events_mt;
  class forwd_otf_dfs_by_events_mt;
  class bckwd_otf_dfs_by_events_mt;

  class fb_otf_saturation_opname;
};


// ******************************************************************
// *                                                                *
// *                     satotf_opname  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::satotf_opname::satotf_opname(const char* n)
 : specialized_opname(n)
{
}

MEDDLY::satotf_opname::~satotf_opname()
{
}

// ============================================================

MEDDLY::satotf_opname::subevent::subevent(forest* f, int* v, int nv, bool firing)
: vars(0), num_vars(nv), root(dd_edge(f)), top(0),
  f(static_cast<expert_forest*>(f)), is_firing(firing)
{
  MEDDLY_DCASSERT(f != 0);
  MEDDLY_DCASSERT(v != 0);
  MEDDLY_DCASSERT(nv > 0);

  vars = new int[num_vars];
  memcpy(vars, v, unsigned(num_vars) * sizeof(int));

  // find top
  top = vars[0];
  for (int i = 1; i < num_vars; i++) {
    if (isLevelAbove(vars[i],top)) top = vars[i];
  }

  uses_extensible_variables = false;
  for (int i = 0; i < num_vars; i++) {
    if (this->f->isExtensibleLevel(vars[i])) {
      uses_extensible_variables = true;
      break;
    }
  }

  unpminterms = pminterms = 0;
  num_minterms = size_minterms = 0;
}

MEDDLY::satotf_opname::subevent::~subevent()
{
  if (vars) delete [] vars;
  for (int i=0; i<num_minterms; i++) {
    delete[] unpminterms[i];
    delete[] pminterms[i];
  }
  free(unpminterms);
  free(pminterms);
}

void MEDDLY::satotf_opname::subevent::clearMinterms()
{
  for (int i=0; i<num_minterms; i++) {
    delete[] unpminterms[i];
    delete[] pminterms[i];
  }
  free(unpminterms);
  free(pminterms);
  unpminterms = pminterms = 0;
  num_minterms = 0;
}


void MEDDLY::satotf_opname::subevent::confirm(otf_relation& rel, int v, int i) {
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}


bool MEDDLY::satotf_opname::subevent::addMinterm(const int* from, const int* to)
{
  /*
  ostream_output out(std::cout);
  out << "Adding minterm: [";
  for (int i = f->getNumVariables(); i >= 0; i--) {
    out << from[i] << " -> " << to[i] << " , ";
  }
  out << "]\n";
  */

  if (num_minterms >= size_minterms) {
    int old_size = size_minterms;
    size_minterms = (0==size_minterms)? 8: MIN(2*size_minterms, 256 + size_minterms);
    unpminterms = (int**) realloc(unpminterms, unsigned(size_minterms) * sizeof(int**));
    pminterms = (int**) realloc(pminterms, unsigned(size_minterms) * sizeof(int**));
    if (0==unpminterms || 0==pminterms) return false; // malloc or realloc failed
    for (int i=old_size; i<size_minterms; i++) {
      unpminterms[i] = 0;
      pminterms[i] = 0;
    }
  }
  if (0==unpminterms[num_minterms]) {
    unpminterms[num_minterms] = new int[f->getNumVariables() + 1];
    MEDDLY_DCASSERT(0==pminterms[num_minterms]);
    pminterms[num_minterms] = new int[f->getNumVariables() + 1];
  }
  // out << "Added minterm: [";
  for (int i = f->getNumVariables(); i >= 0; i--) {
    unpminterms[num_minterms][i] = from[i];
    pminterms[num_minterms][i] = to[i];
    // out << unpminterms[num_minterms][i] << " -> " << pminterms[num_minterms][i] << " , ";
  }
  // out << "]\n";
  expert_domain* d = static_cast<expert_domain*>(f->useDomain());
  for (int i = num_vars - 1; i >= 0; i--) {
    int level = vars[i];
    // expand "to" since the set of unconfirmed local states is always larger
    if (to[level] > 0 && to[level] >= f->getLevelSize(-level)) {
      if (f->isExtensibleLevel(level))
        d->enlargeVariableBound(level, false, -(1+to[level]));
      else
        d->enlargeVariableBound(level, false, 1+to[level]);
    }
  }
  num_minterms++;
  return true;
}

void MEDDLY::satotf_opname::subevent::buildRoot() {
  if (0 == num_minterms) return;
  /*
  ostream_output out(std::cout);
  out << "Building subevent from " << num_minterms << " minterms\n";
  for (int i = 0; i < num_minterms; i++) {
    out << "minterm[" << i << "]: [ ";
    for (int j = f->getNumVariables(); j >= 0; j--) {
      out << unpminterms[i][j] << " -> " << pminterms[i][j] << ", ";
    }
    out << " ]\n";
  }
  */
  if (usesExtensibleVariables()) {
    dd_edge sum(root);
    f->createEdge(unpminterms, pminterms, num_minterms, sum);
    num_minterms = 0;
    root += sum;
  } else {
    f->createEdge(unpminterms, pminterms, num_minterms, root);
  }
  // out << "Equivalent event: " << root.getNode() << "\n";
  // out << "Result: ";
  // root.show(out, 2);
}


void MEDDLY::satotf_opname::subevent::showInfo(output& out) const {
  int num_levels = f->getDomain()->getNumVariables();
  for (int i = 0; i < num_minterms; i++) {
    out << "minterm[" << i << "]: ";
    for (int lvl = num_levels; lvl > 0; lvl--) {
      out << unpminterms[i][lvl] << " -> " << pminterms[i][lvl] << ", ";
    }
    out << "]\n";
  }
  root.show(out, 2);
}

long MEDDLY::satotf_opname::subevent::mintermMemoryUsage() const {
  long n_minterms = 0L;
  for (int i = 0; i < size_minterms; i++) {
    if (unpminterms[i] != 0) n_minterms++;
  }
  return long(n_minterms * 2) * long(f->getDomain()->getNumVariables()) * long(sizeof(int));
}

// ============================================================

MEDDLY::satotf_opname::event::event(subevent** p, int np)
{
  if (p == 0 || np <= 0 || p[0]->getForest() == 0)
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  f = p[0]->getForest();
  for (int i=1; i<np; i++) {
    if (p[i]->getForest() != f) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }

  num_subevents = np;
  subevents = new subevent*[np];
  for (int i = 0; i < np; i++) subevents[i] = p[i];

  top = p[0]->getTop();
  for (int i = 1; i < np; i++) {
    if (top < p[i]->getTop()) top = p[i]->getTop();
  }

  // Find the variable that effect this event from the list of subevents.
  // TODO:
  // Not efficient. p[i] is a sorted list of integers.
  // Should be able to insert in O(n) time
  // where n is the sum(p[i]->getNumVars).
#if 0
  bool all_firing_subevents = true;
#endif
#ifdef DEVELOPMENT_CODE
  bool all_enabling_subevents = true;
#endif
  std::set<int> sVars;
  std::set<int> firingVars;

  for (int i = 0; i < np; i++) {
    const int* subeventVars = p[i]->getVars();
    sVars.insert(subeventVars, subeventVars+p[i]->getNumVars());
    if (p[i]->isFiring()) {
#ifdef DEVELOPMENT_CODE
      all_enabling_subevents = false;
#endif
      firingVars.insert(subeventVars, subeventVars+p[i]->getNumVars());
    } else {
#if 0
      all_firing_subevents = false;
#endif
    }
  }

  MEDDLY_DCASSERT(all_enabling_subevents || !firingVars.empty());

#if 0
  is_disabled = (all_enabling_subevents || all_firing_subevents);
#else
  is_disabled = false;
#endif

  num_vars = sVars.size();
  vars = new int[num_vars];
  int* curr = &vars[0];
  for (std::set<int>::iterator it=sVars.begin(); it!=sVars.end(); ) {
    *curr++ = *it++;
  }

  num_firing_vars = firingVars.size();
  firing_vars = new int[num_firing_vars];
  curr = &firing_vars[0];
  for (std::set<int>::iterator it=firingVars.begin(); it!=firingVars.end(); ) {
    *curr++ = *it++;
  }

  root = dd_edge(f);
  event_mask = dd_edge(f);
  event_mask_from_minterm = 0;
  event_mask_to_minterm = 0;
  needs_rebuilding = is_disabled? false: true;
}

MEDDLY::satotf_opname::event::~event()
{
  for (int i=0; i<num_subevents; i++) delete subevents[i];
  delete[] subevents;
  delete[] vars;
  delete[] firing_vars;
  delete[] event_mask_from_minterm;
  delete[] event_mask_to_minterm;
}

void MEDDLY::satotf_opname::event::buildEventMask()
{
  MEDDLY_DCASSERT(num_subevents > 0);
  MEDDLY_DCASSERT(f);

  if (0 == event_mask_from_minterm) {
    const size_t minterm_size = size_t(f->getNumVariables()+1);
    event_mask_from_minterm = new int[minterm_size];
    event_mask_to_minterm = new int[minterm_size];

    for (unsigned i = 0; i < minterm_size; i++) {
      event_mask_from_minterm[i] = MEDDLY::DONT_CARE;
      event_mask_to_minterm[i] = MEDDLY::DONT_CHANGE;
    }

    for (int i = 0; i < num_firing_vars; i++) {
      event_mask_to_minterm[firing_vars[i]] = MEDDLY::DONT_CARE;
    }
  }

  f->createEdge(&event_mask_from_minterm, &event_mask_to_minterm, 1, event_mask);
#ifdef DEBUG_EVENT_MASK
  printf("event_mask: %d\n" , event_mask.getNode());
  ostream_output out(std::cout);
  event_mask.show(out, 2);
#endif
}


bool MEDDLY::satotf_opname::event::rebuild()
{
  MEDDLY_DCASSERT(num_subevents > 0);
  if (is_disabled) return false;
  if (!needs_rebuilding) return false;
  needs_rebuilding = false;

  // An event is a conjunction of sub-events (or sub-functions).
  for (int i = 0; i < num_subevents; i++) {
    subevents[i]->buildRoot();
  }
  buildEventMask();

  dd_edge e(event_mask);
  for (int i = 0; i < num_subevents; i++) {
    e *= subevents[i]->getRoot();
  }

  /*
  if (e.getNode() == 0) {
    ostream_output out(std::cout);
    f->useDomain()->showInfo(out);
    out << "subevent: " << event_mask.getNode() << "\n";
    event_mask.show(out, 2);
    for (int i = 0; i < num_subevents; i++) {
    out << "subevent: " << subevents[i]->getRoot().getNode() << "\n";
    subevents[i]->getRoot().show(out, 2);
    }
  }
  ostream_output out(std::cout);
  for (int i = 0; i < num_subevents; i++) {
  out << "subevent: " << subevents[i]->getRoot().getNode() << "\n";
  subevents[i]->getRoot().show(out, 2);
  }
  e.show(out, 2);
  */
  if (e == root) return false;
  root = e;
  return true;
}

void MEDDLY::satotf_opname::event::enlargeVariables()
{
  expert_domain* ed = static_cast<expert_forest*>(f)->useExpertDomain();
  for (int i = 0; i < num_vars; i++) {
    int unprimed = ABS(vars[i]);
    int primed = -unprimed;
    int unprimedSize = f->getLevelSize(unprimed);
    int primedSize = f->getLevelSize(primed);
    if (unprimedSize < primedSize) {
      variable* vh = ed->getExpertVar(unprimed);
      if (vh->isExtensible())
        vh->enlargeBound(false, -primedSize);
      else
        vh->enlargeBound(false, primedSize);
    }
    MEDDLY_DCASSERT(f->getLevelSize(unprimed) == f->getLevelSize(primed));
  }
}


void MEDDLY::satotf_opname::event::showInfo(output& out) const {
  for (int i = 0; i < num_subevents; i++) {
    out << "subevent " << i << "\n";
    subevents[i]->showInfo(out);
  }
}

long MEDDLY::satotf_opname::event::mintermMemoryUsage() const {
  long usage = 0;
  for (int i = 0; i < num_subevents; i++) {
    usage += subevents[i]->mintermMemoryUsage();
  }
  return usage;
}

// ============================================================

MEDDLY::satotf_opname::otf_relation::otf_relation(forest* inmdd,
  forest* mxd, forest* outmdd, event** E, int ne)
: insetF(static_cast<expert_forest*>(inmdd)),
  mxdF(static_cast<expert_forest*>(mxd)),
  outsetF(static_cast<expert_forest*>(outmdd))
{
  if (0==insetF || 0==outsetF || 0==mxdF) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);

  // Check for same domain
  if (
    (insetF->getDomain() != mxdF->getDomain()) ||
    (outsetF->getDomain() != mxdF->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  // for now, anyway, inset and outset must be same forest
  if (insetF != outsetF)
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

  // Check forest types
  if (
    insetF->isForRelations()    ||
    !mxdF->isForRelations()     ||
    outsetF->isForRelations()   ||
    (insetF->getRangeType() != mxdF->getRangeType())        ||
    (outsetF->getRangeType() != mxdF->getRangeType())       ||
    (insetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)   ||
    (mxdF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)     ||
    (outsetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  // Forests are good; set number of variables
  num_levels = mxdF->getDomain()->getNumVariables() + 1;

  // Build the events-per-level data structure
  // (0) Initialize
  num_events_by_top_level = new int[num_levels];
  num_events_by_level = new int[num_levels];
  memset(num_events_by_top_level, 0, sizeof(int)*unsigned(num_levels));
  memset(num_events_by_level, 0, sizeof(int)*unsigned(num_levels));

  // (1) Determine the number of events per level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nVars = E[i]->getNumVars();
    const int* vars = E[i]->getVars();
    for (int j = 0; j < nVars; j++) {
      num_events_by_level[vars[j]]++;
    }
    num_events_by_top_level[E[i]->getTop()]++;
  }

  // (2) Allocate events[i]
  events_by_top_level = new event**[num_levels];
  events_by_level = new event**[num_levels];
  for (int i = 0; i < num_levels; i++) {
    events_by_top_level[i] = num_events_by_top_level[i] > 0
      ? new event*[num_events_by_top_level[i]]: 0;
    events_by_level[i] = num_events_by_level[i] > 0
      ? new event*[num_events_by_level[i]]: 0;
    num_events_by_top_level[i] = 0; // reset this; to be used by the next loop
    num_events_by_level[i] = 0; // reset this; to be used by the next loop
  }

  // (3) Store events by level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nVars = E[i]->getNumVars();
    const int* vars = E[i]->getVars();
    for (int j = 0; j < nVars; j++) {
      int level = vars[j];
      events_by_level[level][num_events_by_level[level]++] = E[i];
    }
    int level = E[i]->getTop();
    events_by_top_level[level][num_events_by_top_level[level]++] = E[i];
    E[i]->markForRebuilding();
  }

  // Build the subevents-per-level data structure
  // (0) Initialize
  num_subevents_by_level = new int[num_levels];
  memset(num_subevents_by_level, 0, sizeof(int)*unsigned(num_levels));

  // (1) Determine the number of subevents per level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nse = E[i]->getNumOfSubevents();
    subevent** se = E[i]->getSubevents();
    for (int j = 0; j < nse; j++) {
      int nVars = se[j]->getNumVars();
      const int* vars = se[j]->getVars();
      for (int k = 0; k < nVars; k++) {
        num_subevents_by_level[vars[k]]++;
      }
    }
  }

  // (2) Allocate subevents[i]
  subevents_by_level = new subevent**[num_levels];
  for (int i = 0; i < num_levels; i++) {
    subevents_by_level[i] = num_subevents_by_level[i] > 0
      ? new subevent*[num_subevents_by_level[i]]: 0;
    num_subevents_by_level[i] = 0; // reset this; to be used by the next loop
  }

  // (3) Store subevents by level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nse = E[i]->getNumOfSubevents();
    subevent** se = E[i]->getSubevents();
    for (int j = 0; j < nse; j++) {
      int nVars = se[j]->getNumVars();
      const int* vars = se[j]->getVars();
      for (int k = 0; k < nVars; k++) {
        int level = vars[k];
        subevents_by_level[level][num_subevents_by_level[level]++] = se[j];
      }
    }
  }

  // confirmed local states
  confirmed = new bool*[num_levels];
  size_confirmed = new int[num_levels];
  num_confirmed = new int[num_levels];
  confirmed[0] = 0;
  for (int i = 1; i < num_levels; i++) {
    int level_size = mxdF->getLevelSize(-i);
    confirmed[i] = (bool*) malloc(unsigned(level_size) * sizeof(bool));
    for (int j = 0; j < level_size; j++) {
      confirmed[i][j] = false;
    }
    size_confirmed[i] = level_size;
    num_confirmed[i] = 0;
  }

  // TBD - debug here
}

MEDDLY::satotf_opname::otf_relation::~otf_relation()
{
  // ostream_output out(std::cout);
  // showInfo(out);
  for (int i = 0; i < num_levels; i++) {
    delete[] subevents_by_level[i];
    delete[] events_by_level[i];
    delete[] events_by_top_level[i];
    free(confirmed[i]);
  }
  delete[] subevents_by_level;
  delete[] events_by_level;
  delete[] events_by_top_level;
  delete[] num_subevents_by_level;
  delete[] num_events_by_level;
  delete[] num_events_by_top_level;
  delete[] confirmed;
  delete[] size_confirmed;
  delete[] num_confirmed;
}

void MEDDLY::satotf_opname::otf_relation::clearMinterms()
{
  for (int level = 1; level < num_levels; level++) {
    // Get subevents affected by this level, and rebuild them.
    int nSubevents = num_subevents_by_level[level];
    for (int i = 0; i < nSubevents; i++) {
      subevents_by_level[level][i]->clearMinterms();
    }
  }
}


bool MEDDLY::satotf_opname::otf_relation::confirm(int level, int index)
{
  // For each subevent that affects this level:
  //    (1) call subevent::confirm()
  //    (2) for each level k affected by the subevent,
  //        (a) enlarge variable bound of k to variable bound of k'

  enlargeConfirmedArrays(level, index+1);

  MEDDLY_DCASSERT(size_confirmed[level] > index);
  if (isConfirmed(level, index)) return false;

  // Get subevents affected by this level, and rebuild them.
  int nSubevents = num_subevents_by_level[level];
  for (int i = 0; i < nSubevents; i++) {
    subevents_by_level[level][i]->confirm(const_cast<otf_relation&>(*this),
        level, index);
  }

  // Get events affected by this level, and mark them stale.
  const int nEvents = num_events_by_level[level];
  for (int i = 0; i < nEvents; i++) {
    // events_by_level[level][i]->enlargeVariables();
    events_by_level[level][i]->markForRebuilding();
  }

  confirmed[level][index] = true;
  num_confirmed[level]++;
  return true;
}


void findConfirmedStates(MEDDLY::satotf_opname::otf_relation* rel,
    bool** confirmed, int* num_confirmed,
    MEDDLY::node_handle mdd, int level,
    std::set<MEDDLY::node_handle>& visited) {
  if (level == 0) return;
  if (visited.find(mdd) != visited.end()) return;

  MEDDLY::expert_forest* insetF = rel->getInForest();
  int mdd_level = insetF->getNodeLevel(mdd);
  if (MEDDLY::isLevelAbove(level, mdd_level)) {
    // skipped level; confirm all local states at this level
    // go to the next level
    int level_size = insetF->getLevelSize(level);
    for (int i = 0; i < level_size; i++) {
      if (!confirmed[level][i]) {
        rel->confirm(level, i);
      }
    }
    findConfirmedStates(rel, confirmed, num_confirmed, mdd, level-1, visited);
  } else {
    if (MEDDLY::isLevelAbove(mdd_level, level)) {
      throw MEDDLY::error(MEDDLY::error::INVALID_VARIABLE, __FILE__, __LINE__);
    }
    // mdd_level == level
    visited.insert(mdd);
    MEDDLY::unpacked_node *nr = insetF->newUnpacked(mdd, MEDDLY::SPARSE_ONLY);
    for (unsigned i = 0; i < nr->getNNZs(); i++) {
      if (!confirmed[level][nr->i(i)]) {
        rel->confirm(level, int(nr->i(i)));
      }
      findConfirmedStates(rel, confirmed, num_confirmed, nr->d(i), level-1, visited);
    }
    MEDDLY::unpacked_node::recycle(nr);
  }
}

void MEDDLY::satotf_opname::otf_relation::confirm(const dd_edge& set)
{
  // Perform a depth-first traversal of set:
  //    At each level, mark all enabled states as confirmed.

  // Enlarge the confirmed arrays if needed
  for (int i = 1; i < num_levels; i++) {
      enlargeConfirmedArrays(i, mxdF->getLevelSize(-i));
  }
  std::set<node_handle> visited;
  findConfirmedStates(const_cast<otf_relation*>(this),
      confirmed, num_confirmed, set.getNode(), num_levels-1, visited);

  // ostream_output out(std::cout);
  // showInfo(out);
}



void MEDDLY::satotf_opname::otf_relation::enlargeConfirmedArrays(int level, int sz)
{
#if 0
  int curr = size_confirmed[level];
  int needed = mxdF->getLevelSize(-level);
  int new_size = curr*2 < needed? needed: curr*2;
  confirmed[level] = (bool*) realloc(confirmed[level], new_size * sizeof(bool));
  if (confirmed[level] == 0) throw error::INSUFFICIENT_MEMORY;
  for ( ; curr < new_size; ) confirmed[level][curr++] = false;
  size_confirmed[level] = curr;
#else
  if (sz <= size_confirmed[level]) return;
  sz = MAX( sz , size_confirmed[level]*2 );
  confirmed[level] = (bool*) realloc(confirmed[level], unsigned(sz) * sizeof(bool));
  if (confirmed[level] == 0) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  for (int i = size_confirmed[level]; i < sz; i++) confirmed[level][i] = false;
  size_confirmed[level] = sz;
#endif
}


void MEDDLY::satotf_opname::otf_relation::showInfo(output &strm) const
{
  for (int level = 1; level < num_levels; level++) {
    for (int ei = 0; ei < getNumOfEvents(level); ei++) {
      events_by_top_level[level][ei]->rebuild();
    }
  }
  strm << "On-the-fly relation info:\n";
  for (int level = 1; level < num_levels; level++) {
    for (int ei = 0; ei < getNumOfEvents(level); ei++) {
      strm << "Level: " << level << ", Event:[" << ei << "]: needs rebuilding? "
        << (events_by_top_level[level][ei]->needsRebuilding()? "true": "false")
        << "\n";
      strm << "Depends on levels: [";
      for (int j = 0; j < events_by_top_level[level][ei]->getNumVars(); j++) {
        strm << " " << events_by_top_level[level][ei]->getVars()[j] << " ";
      }
      strm << "]\n";
      events_by_top_level[level][ei]->showInfo(strm);
      events_by_top_level[level][ei]->getRoot().show(strm, 2);
    }
  }
}

long MEDDLY::satotf_opname::otf_relation::mintermMemoryUsage() const {
  long usage = 0;
  for (int level = 1; level < num_levels; level++) {
    for (int ei = 0; ei < getNumOfEvents(level); ei++) {
      usage += events_by_top_level[level][ei]->mintermMemoryUsage();
    }
  }
  return usage;
}

void MEDDLY::satotf_opname::otf_relation::bindExtensibleVariables() {
  //
  // Find the bounds for each extensbile variable
  //
  expert_domain* ed = mxdF->useExpertDomain();
  for (int k = 1; k < num_levels; k++) {
    int bound = 0;
    int n_confirmed = 0;

    for (int j = 0; j < size_confirmed[k]; j++) {
      if (confirmed[k][j]) { bound = j+1; n_confirmed++; }
    }

    MEDDLY_DCASSERT(bound > 0);
    MEDDLY_DCASSERT(n_confirmed == num_confirmed[k]);
    ed->enlargeVariableBound(k, false, bound);
  }
}


double MEDDLY::satotf_opname::otf_relation::getArcCount(
  const dd_edge& mask,
  bool count_duplicates)
{
  MEDDLY_DCASSERT(outsetF->isQuasiReduced());
  MEDDLY_DCASSERT(mxdF->isIdentityReduced());

  double arc_count = 0;
  dd_edge mxd_mask(mxdF);

  // Build confirmed mask
  dd_edge confirmed_local_states(outsetF);
  confirmed_local_states.set(bool_Tencoder::value2handle(true));
  for (int k = 1; k < num_levels; k++) {
    node_handle current_node = confirmed_local_states.getNode();
    int current_level = outsetF->getNodeLevel(current_node);
    int next_level = outsetF->upLevel(outsetF->upLevel(current_level));
    MEDDLY_DCASSERT(next_level >= 0);
    unpacked_node* node =
      unpacked_node::newFull(outsetF, next_level, unsigned(outsetF->getLevelSize(next_level)));
    for (unsigned i = 0; i < node->getSize(); i++) {
      node->d_ref(i) =
        confirmed[next_level][i]
        ? outsetF->linkNode(current_node)
        : 0;
    }
    node_handle next_node = outsetF->createReducedNode(-1, node);
    confirmed_local_states.set(next_node);
  }

  dd_edge confirmed_local_states_mask = mask * confirmed_local_states;
  MEDDLY::apply(MEDDLY::CROSS, confirmed_local_states_mask, confirmed_local_states_mask, mxd_mask);

  if (count_duplicates) {
    for (int k = 1; k < num_levels; k++) {
      for (int ei = 0; ei < getNumOfEvents(k); ei++) {
        // start with (num_level-1) to correctly count edges in skipped levels
        dd_edge rg_ei = events_by_top_level[k][ei]->getRoot();
        rg_ei *= mxd_mask;
        arc_count += rg_ei.getCardinality();
      }
    }
  } else {
    // build monolithic
    dd_edge monolithic_nsf(mxdF);
    for (int k = 1; k < num_levels; k++) {
      dd_edge nsf_i(mxdF);
      for (int ei = 0; ei < getNumOfEvents(k); ei++) {
        nsf_i += events_by_top_level[k][ei]->getRoot();
      }
      monolithic_nsf += nsf_i;
    }
    monolithic_nsf *= mxd_mask;
    arc_count = monolithic_nsf.getCardinality();
  }
  return arc_count;
}


void MEDDLY::satotf_opname::otf_relation::getBoundedMonolithicNSF(dd_edge &root) const
{
  fprintf(stderr, "otf_relation::getBoundedMonolithicNSF not implemented yet\n");
  throw 42;

  /*

  TBD - convert to dd_edge based plus intersection

  */

  /*

  //
  // Build the union all events
  //
  dd_edge* monolithic_nsf = new dd_edge(mxdF);
  for (int k = 1; k < num_levels; k++) {
    for (int ei = 0; ei < getNumOfEvents(k); ei++) {
      (*monolithic_nsf) += events_by_top_level[k][ei]->getRoot();
    }
  }

  //
  // Find the bounds for each extensbile variable
  //
  int i = 1;
  for ( ; i < num_levels && !mxdF->isExtensibleLevel(i); i++);
  const bool found_extensible_variable = (i < num_levels);

  MEDDLY::node_handle bounded_monolithic_nsf = 0;

  if (!found_extensible_variable) {
    bounded_monolithic_nsf = mxdF->linkNode(monolithic_nsf->getNode());
  } else {
    int bounds[num_levels];
    bounds[0] = 0;
    for (int i = 1; i < num_levels; i++) {
      int j = 0;
      for (j = size_confirmed[i]-1; j >= 0 && !confirmed[i][j]; j--);
      bounds[i] = j+1;
      MEDDLY_DCASSERT(bounds[i] > 0);
    }

    //
    // Convert monolithic_nsf to a bounded event based on confirmed local states
    // - Recursively build a new MDD while pruning all extensible nodes
    //   to the size of the bounded variables.
    //
    std::unordered_map<MEDDLY::node_handle, MEDDLY::node_handle> cache;
    bounded_monolithic_nsf =
      getBoundedMxd(monolithic_nsf->getNode(), bounds, num_levels, cache);

    // clear cache
    for (auto i : cache) mxdF->unlinkNode(i.second);

    // set bounded variable sizes
    expert_domain* ed = mxdF->useExpertDomain();
    for (int i = 1; i < num_levels; i++) {
      ed->enlargeVariableBound(i, false, bounds[i]);
    }
  }
  delete monolithic_nsf;

  return bounded_monolithic_nsf;
  */
}

/*

MEDDLY::node_handle MEDDLY::satotf_opname::otf_relation::getBoundedMxd(
    MEDDLY::node_handle mxd,
    const int* bounds,
    int num_levels,
    std::unordered_map<MEDDLY::node_handle, MEDDLY::node_handle>& cache
    )
{
  if (mxdF->isTerminalNode(mxd)) return mxd;
  if (!mxdF->isExtensible(mxd)) return mxdF->linkNode(mxd);

  // mxd is an extensible node:
  //    use bounds[] to make build a non-extensible version and return it

  std::unordered_map<MEDDLY::node_handle, MEDDLY::node_handle>::iterator iter = cache.find(mxd);
  if (iter != cache.end()) {
    return mxdF->linkNode(iter->second);
  }

  int mxd_level = mxdF->getNodeLevel(mxd);
  int result_size = bounds[ABS(mxd_level)];
  unpacked_node *mxd_node = mxdF->newUnpacked(mxd, FULL_ONLY);
  unpacked_node *bounded_node = unpacked_node::newFull(mxdF, mxd_level, result_size);
  MEDDLY::node_handle ext_d = mxd_node->ext_d();

  int i = 0;
  for ( ; i < mxd_node->getSize(); i++) bounded_node->d_ref(i) = mxdF->linkNode(mxd_node->d(i));
  for ( ; i < result_size; i++) bounded_node->d_ref(i) = mxdF->linkNode(ext_d);
  bounded_node->markAsNotExtensible();

  MEDDLY::node_handle result = mxdF->createReducedNode(-1, bounded_node);
  cache[mxd] = mxdF->linkNode(result);
  unpacked_node::recycle(mxd_node);

  return result;
}

*/

// ******************************************************************
// *                                                                *
// *                    otfsat_by_events_opname  class              *
// *                                                                *
// ******************************************************************

/** Simple class to keep compute table happy.
*/
class MEDDLY::otfsat_by_events_opname : public unary_opname {
  static otfsat_by_events_opname* instance;
  public:
    otfsat_by_events_opname();

    static const otfsat_by_events_opname* getInstance();

};

MEDDLY::otfsat_by_events_opname* MEDDLY::otfsat_by_events_opname::instance = 0;

MEDDLY::otfsat_by_events_opname::otfsat_by_events_opname()
 : unary_opname("Otf_Saturate_by_events")
{
}

const MEDDLY::otfsat_by_events_opname* MEDDLY::otfsat_by_events_opname::getInstance()
{
  if (0==instance) instance = new otfsat_by_events_opname;
  return instance;
}

// ******************************************************************
// *                                                                *
// *                      otfsat_by_events_op  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::otfsat_by_events_op : public unary_operation {
    common_otf_dfs_by_events_mt* parent;
  public:
    otfsat_by_events_op(common_otf_dfs_by_events_mt* p,
      expert_forest* argF, expert_forest* resF);
    virtual ~otfsat_by_events_op();

    void saturate(const dd_edge& in, dd_edge& out);
    // node_handle saturate(node_handle mdd);
    node_handle saturate(node_handle mdd, int level);

  protected:
    inline ct_entry_key*
    findSaturateResult(node_handle a, int level, node_handle& b) {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      if (argF->isFullyReduced()) CTsrch->writeI(level);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      b = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveSaturateResult(ct_entry_key* Key,
      node_handle a, node_handle b)
    {
      CTresult[0].reset();
      CTresult[0].writeN(b);
      CT0->addEntry(Key, CTresult[0]);
      return b;
    }
};


// ******************************************************************
// *                                                                *
// *            common_otf_dfs_by_events_mt  class                  *
// *                                                                *
// ******************************************************************

class MEDDLY::common_otf_dfs_by_events_mt : public specialized_operation {
  public:
    common_otf_dfs_by_events_mt(const satotf_opname* opcode,
      satotf_opname::otf_relation* rel);
    virtual ~common_otf_dfs_by_events_mt();

    virtual void compute(const dd_edge& a, dd_edge &c);
    virtual void saturateHelper(unpacked_node& mdd) = 0;

  protected:
    inline ct_entry_key*
    findResult(node_handle a, node_handle b, node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CTsrch->writeN(b);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(ct_entry_key* Key,
      node_handle a, node_handle b, node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
      return c;
    }

  protected:
    binary_operation* mddUnion;
    binary_operation* mxdIntersection;
    binary_operation* mxdDifference;

    satotf_opname::otf_relation* rel;

    expert_forest* arg1F;
    expert_forest* arg2F;
    expert_forest* resF;

  protected:
    class indexq {
        static const int NULPTR = -1;
        static const int NOTINQ = -2;
        int* data;
        unsigned size;
        int head;
        int tail;
      public:
        // used by parent for recycling
        indexq* next;
      public:
        indexq();
        ~indexq();
        void resize(unsigned sz);
        inline bool isEmpty() const {
          return NULPTR == head;
        }
        inline void add(int i) {
          MEDDLY_CHECK_RANGE(0, i, size);
          if (NOTINQ != data[i]) return;
          if (NULPTR == head) {
            // empty list
            head = i;
          } else {
            // not empty list
            MEDDLY_CHECK_RANGE(0, tail, size);
            data[tail] = i;
          }
          tail = i;
          data[i] = NULPTR;
        }
        inline int remove() {
          MEDDLY_CHECK_RANGE(0, head, size);
          int ans = head;
          head = data[head];
          data[ans] = NOTINQ;
          MEDDLY_CHECK_RANGE(0, ans, size);
          return ans;
        }
    };

  protected:
    class charbuf {
      public:
        char* data;
        unsigned size;
        charbuf* next;
      public:
        charbuf();
        ~charbuf();
        void resize(unsigned sz);
    };

  private:
    indexq* freeqs;
    charbuf* freebufs;

  protected:
    inline indexq* useIndexQueue(unsigned sz) {
      indexq* ans;
      if (freeqs) {
        ans = freeqs;
        freeqs = freeqs->next;
      } else {
        ans = new indexq();
      }
      MEDDLY_DCASSERT(ans);
      ans->resize(sz);
      ans->next = 0;
      return ans;
    }
    inline void recycle(indexq* a) {
      MEDDLY_DCASSERT(a);
      MEDDLY_DCASSERT(a->isEmpty());
      a->next = freeqs;
      freeqs = a;
    }

    inline charbuf* useCharBuf(unsigned sz) {
      charbuf* ans;
      if (freebufs) {
        ans = freebufs;
        freebufs = freebufs->next;
      } else {
        ans = new charbuf();
      }
      MEDDLY_DCASSERT(ans);
      ans->resize(sz);
      ans->next = 0;
      return ans;
    }
    inline void recycle(charbuf* a) {
      MEDDLY_DCASSERT(a);
      a->next = freebufs;
      freebufs = a;
    }

    inline virtual bool checkForestCompatibility() const {
      return true;
    }
};

// ******************************************************************
// *                                                                *
// *           otfsat_by_events_op  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::otfsat_by_events_op
::otfsat_by_events_op(common_otf_dfs_by_events_mt* p,
  expert_forest* argF, expert_forest* resF)
  : unary_operation(otfsat_by_events_opname::getInstance(), 1, argF, resF)
{
  parent = p;

  const char* name = otfsat_by_events_opname::getInstance()->getName();
  ct_entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new ct_entry_type(name, "NI:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(3, resF);
  } else {
    et = new ct_entry_type(name, "N:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(2, resF);
  }
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::otfsat_by_events_op::~otfsat_by_events_op()
{
  removeAllComputeTableEntries();
}

void MEDDLY::otfsat_by_events_op::saturate(const dd_edge& in, dd_edge& out)
{
  // Saturate
  out.set( saturate(in.getNode(), argF->getNumVariables()) );
}

/*
MEDDLY::node_handle MEDDLY::otfsat_by_events_op::saturate(node_handle mdd)
{
  // Saturate
  return saturate(mdd, argF->getNumVariables());
}
*/

MEDDLY::node_handle
MEDDLY::otfsat_by_events_op::saturate(node_handle mdd, int k)
{
#ifdef DEBUG_DFS
  printf("mdd: %d, k: %d\n", mdd, k);
#endif

  // terminal condition for recursion
  if (argF->isTerminalNode(mdd)) return mdd;

  // search compute table
  node_handle n = 0;
  ct_entry_key* Key = findSaturateResult(mdd, k, n);
  if (0==Key) return n;

  const unsigned sz = unsigned(argF->getLevelSize(k));    // size
  const int mdd_level = argF->getNodeLevel(mdd);          // mdd level

#ifdef DEBUG_DFS
  printf("mdd: %d, level: %d, size: %d, mdd_level: %d\n",
      mdd, k, sz, mdd_level);
#endif

  unpacked_node* nb = unpacked_node::newFull(resF, k, sz);
  // Initialize mdd reader
  unpacked_node *mddDptrs = unpacked_node::New();
  if (mdd_level < k) {
    mddDptrs->initRedundant(argF, k, mdd, true);
  } else {
    argF->unpackNode(mddDptrs, mdd, FULL_ONLY);
  }

  // Do computation
  for (unsigned i=0; i<sz; i++) {
    nb->d_ref(i) = mddDptrs->d(i) ? saturate(mddDptrs->d(i), k-1) : 0;
  }

  // Cleanup
  unpacked_node::recycle(mddDptrs);

  parent->saturateHelper(*nb);
  n = resF->createReducedNode(-1, nb);

  // save in compute table
  saveSaturateResult(Key, mdd, n);

#ifdef DEBUG_DFS
  resF->showNodeGraph(stdout, n);
#endif

  return n;
}


// ******************************************************************
// *                                                                *
// *           common_otf_dfs_by_events_mt  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::common_otf_dfs_by_events_mt::common_otf_dfs_by_events_mt(
  const satotf_opname* opcode,
  satotf_opname::otf_relation* relation)
: specialized_operation(opcode, 1)
{
  mddUnion = 0;
  mxdIntersection = 0;
  mxdDifference = 0;
  freeqs = 0;
  freebufs = 0;
  rel = relation;
  arg1F = static_cast<expert_forest*>(rel->getInForest());
  arg2F = static_cast<expert_forest*>(rel->getRelForest());
  resF = static_cast<expert_forest*>(rel->getOutForest());

  registerInForest(arg1F);
  registerInForest(arg2F);
  registerInForest(resF);

  ct_entry_type* et = new ct_entry_type(opcode->getName(), "NN:N");
  et->setForestForSlot(0, arg1F);
  et->setForestForSlot(1, arg2F);
  et->setForestForSlot(3, resF);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::common_otf_dfs_by_events_mt::~common_otf_dfs_by_events_mt()
{
  if (rel->autoDestroy()) delete rel;
  unregisterInForest(arg1F);
  unregisterInForest(arg2F);
  unregisterInForest(resF);
}

void MEDDLY::common_otf_dfs_by_events_mt
::compute(const dd_edge &a, dd_edge &c)
{
  // Initialize operations
  mddUnion = getOperation(UNION, resF, resF, resF);
  MEDDLY_DCASSERT(mddUnion);

  mxdIntersection = getOperation(INTERSECTION, arg2F, arg2F, arg2F);
  MEDDLY_DCASSERT(mxdIntersection);

  mxdDifference = getOperation(DIFFERENCE, arg2F, arg2F, arg2F);
  MEDDLY_DCASSERT(mxdDifference);

#ifdef DEBUG_INITIAL
  printf("Calling saturate for states:\n");
  a.show(stdout, 2);
#endif
#ifdef DEBUG_NSF
  printf("Calling saturate for NSF:\n");
  // b.show(stdout, 2);
#endif

  // Execute saturation operation
  otfsat_by_events_op* so = new otfsat_by_events_op(this, arg1F, resF);
  so->saturate(a, c);

  // Cleanup
  while (freeqs) {
    indexq* t = freeqs;
    freeqs = t->next;
    delete t;
  }
  while (freebufs) {
    charbuf* t = freebufs;
    freebufs = t->next;
    delete t;
  }
  delete so;
}

// ******************************************************************
// *       common_otf_dfs_by_events_mt::indexq  methods                 *
// ******************************************************************

MEDDLY::common_otf_dfs_by_events_mt::indexq::indexq()
{
  data = 0;
  size = 0;
  head = NULPTR;
}

MEDDLY::common_otf_dfs_by_events_mt::indexq::~indexq()
{
  free(data);
}

void MEDDLY::common_otf_dfs_by_events_mt::indexq::resize(unsigned sz)
{
  if (sz <= size) return;
  data = (int*) realloc(data, sz * sizeof(int));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);

  for (; size < sz; size++) data[size] = NOTINQ;
}

// ******************************************************************
// *       common_otf_dfs_by_events_mt::charbuf methods             *
// ******************************************************************

MEDDLY::common_otf_dfs_by_events_mt::charbuf::charbuf()
{
  data = 0;
  size = 0;
}

MEDDLY::common_otf_dfs_by_events_mt::charbuf::~charbuf()
{
  free(data);
}

void MEDDLY::common_otf_dfs_by_events_mt::charbuf::resize(unsigned sz)
{
  if (sz <= size) return;
  data = (char*) realloc(data, sz * sizeof(char));
  if (0==data)
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *             forwd_otf_dfs_by_events_mt class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_otf_dfs_by_events_mt : public common_otf_dfs_by_events_mt {
  public:
    forwd_otf_dfs_by_events_mt(const satotf_opname* opcode,
    satotf_opname::otf_relation* rel);
  protected:
    virtual void saturateHelper(unpacked_node& mdd);
    node_handle recFire(node_handle mdd, node_handle mxd);
    void recFireHelper(const unsigned, const int, const node_handle, const node_handle,
        unpacked_node*, unpacked_node*);
};

MEDDLY::forwd_otf_dfs_by_events_mt::forwd_otf_dfs_by_events_mt(
  const satotf_opname* opcode,
  satotf_opname::otf_relation* rel)
  : common_otf_dfs_by_events_mt(opcode, rel)
{
}


void MEDDLY::forwd_otf_dfs_by_events_mt::saturateHelper(unpacked_node& nb)
{
  const int level = nb.getLevel();
  const int nEventsAtThisLevel = rel->getNumOfEvents(level);
  if (0 == nEventsAtThisLevel) return;

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node** Ru = new unpacked_node*[nEventsAtThisLevel];
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
    rel->rebuildEvent(level, ei);
    const dd_edge& mxd = rel->getEvent(level, ei);
    if (0==mxd.getNode()) {
      Ru[ei] = 0;
    } else {
      Ru[ei] = unpacked_node::New();
      const int eventLevel = mxd.getLevel();
      if (ABS(eventLevel) < level || eventLevel < 0) {
        // Takes care of two situations:
        // - skipped unprimed level (due to Fully Reduced)
        // - skipped unprimed and primed levels (due to Fully Identity Reduced)
        Ru[ei]->initRedundant(arg2F, level, mxd.getNode(), true);
      } else {
        arg2F->unpackNode(Ru[ei], mxd.getNode(), FULL_ONLY);
      }
    }
  }
  unpacked_node* Rp = unpacked_node::New();

  dd_edge nbdj(resF), newst(resF);

  //      Node reader auto expands when passed by reference
  //      Node builder can be expanded via a call to unpacked_node::resize()
  //      Queue can be expanded via a call to indexq::resize()

  // indexes to explore
  indexq* queue = useIndexQueue(nb.getSize());
  for (unsigned i = 0; i < nb.getSize(); i++) {
    if (nb.d(i)) queue->add(int(i));
  }

  // explore indexes
  while (!queue->isEmpty()) {
    const unsigned i = unsigned(queue->remove());

    MEDDLY_DCASSERT(nb.d(i));

    for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
      // If event i needs rebuilding,
      //    Rebuild it, and
      //    Update the event reader
      if (rel->rebuildEvent(level, ei)) {
        const dd_edge& mxd = rel->getEvent(level, ei);
        if (0==mxd.getNode()) {
          if (Ru[ei]) {
            unpacked_node::recycle(Ru[ei]);
            Ru[ei] = 0;
          }
        } else {
          if (0==Ru[ei]) {
            Ru[ei] = unpacked_node::New();
          }
          const int eventLevel = mxd.getLevel();
          if (ABS(eventLevel) < level || eventLevel < 0) {
            // Takes care of two situations:
            // - skipped unprimed level (due to Fully Reduced)
            // - skipped unprimed and primed levels (due to Fully Identity Reduced)
            Ru[ei]->initRedundant(arg2F, level, mxd.getNode(), true);
          } else {
            arg2F->unpackNode(Ru[ei], mxd.getNode(), FULL_ONLY);
          }
        }
      }
      // check if row i of the event ei is empty
      if (0 == Ru[ei]) continue;
      MEDDLY_DCASSERT(!Ru[ei]->isExtensible());
      node_handle ei_i = (i < Ru[ei]->getSize())
                        ? Ru[ei]->d(i)
                        : (Ru[ei]->isExtensible() ? Ru[ei]->ext_d() : 0);
      if (0 == ei_i) continue;

      // grab column (TBD: build these ahead of time?)
      const int dlevel = arg2F->getNodeLevel(ei_i);

      if (dlevel == -level) {
        arg2F->unpackNode(Rp, ei_i, SPARSE_ONLY);
      } else {
        Rp->initIdentity(arg2F, -level, i, ei_i, false);
      }

      MEDDLY_DCASSERT(!Rp->isExtensible());

      for (unsigned jz=0; jz<Rp->getNNZs(); jz++) {
        const unsigned j = Rp->i(jz);
        if (j < nb.getSize() && -1==nb.d(j)) continue;  // nothing can be added to this set

        node_handle newstates = recFire(nb.d(i), Rp->d(jz));
        if (newstates == 0) continue;

        // Confirm local state
        rel->confirm(level, int(j));

        if (j >= nb.getSize()) {
          unsigned oldSize = nb.getSize();
          // resize the node builder, and clear out the new entries
          nb.resize(j+1);
          while(oldSize < nb.getSize()) { nb.d_ref(oldSize++) = 0; }
          // resize the queue, and clear out the new entries
          queue->resize(nb.getSize());
        }

        bool updated = true;

        if (0 == nb.d(j)) {
          nb.d_ref(j) = newstates;
        } else {
          nbdj.set(nb.d(j));      // clobber
          newst.set(newstates);   // clobber
          mddUnion->computeTemp(nbdj, newst, nbdj);
          updated = (nbdj.getNode() != nb.d(j));
          nb.set_d(j, nbdj);
        }

        if (updated) queue->add((int)j);
      } // for j
    } // for all events, ei
  } // while there are indexes to explore

  // cleanup
  unpacked_node::recycle(Rp);
  for (int ei = 0; ei < nEventsAtThisLevel; ei++) unpacked_node::recycle(Ru[ei]);
  delete[] Ru;
  recycle(queue);
}


// Same as post-image, except we saturate before reducing.
MEDDLY::node_handle MEDDLY::forwd_otf_dfs_by_events_mt::recFire(
  node_handle mdd, node_handle mxd)
{
  //      Node builder expansion
  //          - done
  //      MXD doesnt expand but unconfirmed > confirmed
  //          - so mdd may need to expand to accomodate rec_fire result
  //          - No need: primed and unprimed variables are of the same size
  //      Check if createReduce and handle a unpacked_node of size
  //          larger than the variable size
  //          - No
  //      Can we build a unpacked_node with the size of the primed variable
  //          and save some trouble?
  //          - No need: primed and unprimed variables are of the same size

  // termination conditions
  if (mxd == 0 || mdd == 0) return 0;
  if (arg2F->isTerminalNode(mxd)) {
    if (arg1F->isTerminalNode(mdd)) {
      return resF->handleForValue(1);
    }
    // mxd is identity
    if (arg1F == resF)
      return resF->linkNode(mdd);
  }

  // check the cache
  node_handle result = 0;
  ct_entry_key* Key = findResult(mdd, mxd, result);
  if (0==Key) return result;

#ifdef TRACE_RECFIRE
  printf("computing recFire(%d, %d)\n", mdd, mxd);
  printf("  node %3d ", mdd);
  arg1F->showNode(stdout, mdd, 1);
  printf("\n  node %3d ", mxd);
  arg2F->showNode(stdout, mxd, 1);
  printf("\n");
#endif

  mxd = arg2F->linkNode(mxd);

  // check if mxd and mdd are at the same level
  const int mddLevel = arg1F->getNodeLevel(mdd);
  const int mxdLevel = arg2F->getNodeLevel(mxd);
  const int rLevel = MAX(ABS(mxdLevel), mddLevel);
  const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
  unpacked_node* nb = unpacked_node::newFull(resF, rLevel, rSize);

  // Initialize mdd reader
  unpacked_node *A = unpacked_node::New();
  if (mddLevel < rLevel) {
    A->initRedundant(arg1F, rLevel, mdd, true);
  } else {
    arg1F->unpackNode(A, mdd, FULL_ONLY);
  }

  if (mddLevel > ABS(mxdLevel)) {
    //
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (unsigned i=0; i<rSize; i++) {
      nb->d_ref(i) = recFire(A->d(i), mxd);
    }

  } else {
    //
    // Need to process this level in the MXD.
    MEDDLY_DCASSERT(ABS(mxdLevel) >= mddLevel);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node *Ru = unpacked_node::New();
    unpacked_node *Rp = unpacked_node::New();
    if (mxdLevel < 0) {
      Ru->initRedundant(arg2F, rLevel, mxd, false);
    } else {
      arg2F->unpackNode(Ru, mxd, SPARSE_ONLY);
    }

#if 0
    // loop over mxd "rows"
    for (int iz=0; iz<Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);
      if (0==A->d(i)) continue;
      const node_handle pnode = Ru->d(iz);
      if (isLevelAbove(-rLevel, arg2F->getNodeLevel(pnode))) {
        Rp->initIdentity(arg2F, rLevel, i, pnode, false);
      } else {
        arg2F->unpackNode(Rp, pnode, SPARSE_ONLY);
      }

      // loop over mxd "columns"
      for (int jz=0; jz<Rp->getNNZs(); jz++) {
        const int j = Rp->i(jz);
        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle newstates = recFire(A->d(i), Rp->d(jz));
        if (0==newstates) continue;

        if (rel->confirm(rLevel, j)) {
          // Newly confirmed local state
          // Node builder must expand to accomodate j
          if (j >= nb->getSize()) {
            int oldSize = nb->getSize();
            // resize the node builder, and clear out the new entries
            nb->resize(j+1);
            while(oldSize < nb->getSize()) { nb->d_ref(oldSize++) = 0; }
          }
        }

        if (0==nb->d(j)) {
          nb->d_ref(j) = newstates;
          continue;
        }
        // there's new states and existing states; union them.
        const int oldj = nb->d(j);
        nb->d_ref(j) = mddUnion->computeTemp(newstates, oldj);
        resF->unlinkNode(oldj);
        resF->unlinkNode(newstates);
      } // for j
    } // for i
#else
    // loop over mxd "rows"
    for (unsigned iz=0; iz<Ru->getNNZs(); iz++) {
      const unsigned i = Ru->i(iz);
      if (0==A->d(i)) continue;
      recFireHelper(i, rLevel, Ru->d(iz), A->d(i), Rp, nb);
    }
    // loop over the extensible portion of mxd (if any)
    MEDDLY_DCASSERT(!Ru->isExtensible());
    if (Ru->isExtensible()) {
      const node_handle pnode = Ru->ext_d();
      for (unsigned i = Ru->ext_i()+1; i < A->getSize(); i++) {
        if (0 == A->d(i)) continue;
        recFireHelper(i, rLevel, pnode, A->d(i), Rp, nb);
      }
    }
#endif

    unpacked_node::recycle(Rp);
    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);

  saturateHelper(*nb);
  result = resF->createReducedNode(-1, nb);
#ifdef TRACE_ALL_OPS
  printf("computed recfire(%d, %d) = %d\n", mdd, mxd, result);
#endif
#ifdef TRACE_RECFIRE
  printf("computed recfire(%d, %d) = %d\n", mdd, mxd, result);
  printf("  node %3d ", result);
  resF->showNode(stdout, result, 1);
  printf("\n");
#endif
  return saveResult(Key, mdd, mxd, result);
}


void MEDDLY::forwd_otf_dfs_by_events_mt::recFireHelper(
  const unsigned i,
  const int rLevel,
  const node_handle Ru_i,
  const node_handle A_i,
  unpacked_node *Rp,
  unpacked_node* nb)
{
  if (isLevelAbove(-rLevel, arg2F->getNodeLevel(Ru_i))) {
    Rp->initIdentity(arg2F, rLevel, i, Ru_i, false);
  } else {
    arg2F->unpackNode(Rp, Ru_i, SPARSE_ONLY);
  }

  MEDDLY_DCASSERT(!Rp->isExtensible());

  dd_edge nbdj(resF), newst(resF);

  // loop over mxd "columns"
  for (unsigned jz=0; jz<Rp->getNNZs(); jz++) {
    const unsigned j = Rp->i(jz);
    // ok, there is an i->j "edge".
    // determine new states to be added (recursively)
    // and add them
    node_handle newstates = recFire(A_i, Rp->d(jz));
    if (0==newstates) continue;

    // Confirm local state
    rel->confirm(rLevel, int(j));

    if (j >= nb->getSize()) {
      unsigned oldSize = nb->getSize();
      // resize the node builder, and clear out the new entries
      nb->resize(j+1);
      while(oldSize < nb->getSize()) { nb->d_ref(oldSize++) = 0; }
    }

    if (0==nb->d(j)) {
      nb->d_ref(j) = newstates;
    } else {
      // there's new states and existing states; union them.
      nbdj.set(nb->d(j));
      newst.set(newstates);
      mddUnion->computeTemp(nbdj, newst, nbdj);
      nb->set_d(j, nbdj);
    }
  } // for j
}



// ******************************************************************
// *                                                                *
// *                   fb_otf_saturation_opname class               *
// *                                                                *
// ******************************************************************

class MEDDLY::fb_otf_saturation_opname : public satotf_opname {
    bool forward;
  public:
    fb_otf_saturation_opname(bool fwd);
    virtual specialized_operation* buildOperation(arguments* a) const;
};

MEDDLY::fb_otf_saturation_opname::fb_otf_saturation_opname(bool fwd)
 : satotf_opname(fwd ? "OtfSaturationFwd" : "OtfSaturationBack")
{
  forward = fwd;
}

MEDDLY::specialized_operation*
MEDDLY::fb_otf_saturation_opname::buildOperation(arguments* a) const
{
  otf_relation* rel = dynamic_cast<otf_relation*>(a);
  if (0==rel) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);

  //
  // No sanity checks needed here; we did them already when constructing a.
  //

  MEDDLY::specialized_operation* op = 0;
  if (forward)
    op = new forwd_otf_dfs_by_events_mt(this, rel);
  else {
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    // op = new bckwd_otf_dfs_by_events_mt(this, rel);
  }

  // Do we need to delete rel here?
  // No, if needed, do this in the destructor for op.

  return op;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::satotf_opname* MEDDLY::initOtfSaturationForward()
{
  return new fb_otf_saturation_opname(true);
}

MEDDLY::satotf_opname* MEDDLY::initOtfSaturationBackward()
{
  return new fb_otf_saturation_opname(false);
}


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

#include "sat_relations.h"
#include "ops_builtin.h"
#include "operators.h"
#include "minterms.h"
#include "relation_node.h"

#define OUT_OF_BOUNDS -1
#define NOT_KNOWN -2
#define TERMINAL_NODE 1


// ******************************************************************
// *                                                                *
// *                    pregen_relation  methods                    *
// *                                                                *
// ******************************************************************


MEDDLY::pregen_relation::pregen_relation(forest* mxd, unsigned nevents)
{
    if (!mxd) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    mxdF = mxd;
    if  (
            !mxdF->isForRelations() ||
            (mxdF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
        )
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    K = mxd->getMaxLevelIndex();

    num_events = nevents;
    if (num_events) {
        events = new dd_edge[num_events];
        next = new unsigned[num_events];
        for (unsigned e=0; e<num_events; e++) {
            events[e].attach(mxd);
        }
    } else {
        events = nullptr;
        next = nullptr;
    }
    last_event = 0;

    level_index = new unsigned[K+1];
    for (unsigned k=0; k<=K; k++) level_index[k] = 0;   // null pointer
}

MEDDLY::pregen_relation::pregen_relation(forest* mxd)
{
    if (!mxd) throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    mxdF = mxd;
    K = mxd->getMaxLevelIndex();

    events = new dd_edge[K+1];
    for (unsigned k=0; k<=K; k++) {
        events[k].attach(mxd);
    }

    next = nullptr;
    level_index = nullptr;

    num_events = 0;
    last_event = 0;
}


MEDDLY::pregen_relation::~pregen_relation()
{
    delete[] events;
    delete[] next;
    delete[] level_index;
}


void MEDDLY::pregen_relation::addToRelation(const dd_edge &r)
{
  MEDDLY_DCASSERT(mxdF);

  if (r.getForest() != mxdF)  throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

  int k = r.getLevel();
  if (0==k) return;
  if (k<0) k = -k;

  if (0==level_index) {
    // relation is "by levels"

    apply(UNION, events[k], r, events[k]);

  } else {
    // relation is "by events"

    if (isFinalized())            throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    if (last_event >= num_events) throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);

    events[last_event] = r;
    next[last_event] = level_index[k];
    level_index[k] = ++last_event;
  }
}


void MEDDLY::pregen_relation::splitMxd(splittingOption split)
{
  if (split == None) return;
  if (split == MonolithicSplit) {
    unionLevels();
    split = SplitOnly;
  }

  // For each level k, starting from the top level
  //    MXD(k) is the matrix diagram at level k.
  //    Calculate ID(k), the intersection of the diagonals of MXD(k).
  //    Subtract ID(k) from MXD(k).
  //    Add ID(k) to level(ID(k)).

#ifdef DEBUG_FINALIZE_SPLIT
  printf("Splitting events in finalize()\n");
  printf("events array: [");
  for (int i=0; i<=K; i++) {
    if (i) printf(", ");
    printf("%d", events[i]);
  }
  printf("]\n");
#endif

  // Initialize operations
  binary_operation* mxdUnion = UNION(mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdUnion);

  binary_operation* mxdIntersection = INTERSECTION(mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdIntersection);

  binary_operation* mxdDifference = DIFFERENCE(mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdDifference);

  dd_edge maxDiag(mxdF);

  for (int k = int(K); k > 1; k--) {
    if (0 == events[k].getNode()) continue;

    MEDDLY_DCASSERT(ABS(events[k].getLevel() <= k));

    // Initialize unpacked nodes
    unpacked_node* Mu = (isLevelAbove(k, events[k].getLevel()))
      ?   unpacked_node::newRedundant(mxdF, k, events[k].getNode(), FULL_ONLY)
      :   mxdF->newUnpacked(events[k].getNode(), FULL_ONLY)
    ;

    unpacked_node* Mp = unpacked_node::New(mxdF);

    // Read "rows"
    for (unsigned i = 0; i < Mu->getSize(); i++) {
      // Initialize column reader
      if (isLevelAbove(-k, mxdF->getNodeLevel(Mu->down(i)))) {
        Mp->initIdentity(mxdF, -k, i, Mu->down(i), FULL_ONLY);
      } else {
        mxdF->unpackNode(Mp, Mu->down(i), FULL_ONLY);
      }

      // Intersect along the diagonal
      if (0==i) {
        maxDiag.set( mxdF->linkNode(Mp->down(i)) );
      } else {
        dd_edge mpd(mxdF);
        mpd.set( mxdF->linkNode(Mp->down(i)) );
        mxdIntersection->computeTemp(maxDiag, mpd, maxDiag);
      }
    } // for i

    // Cleanup
    unpacked_node::Recycle(Mp);
    unpacked_node::Recycle(Mu);

    if (0 == maxDiag.getNode()) {
#ifdef DEBUG_FINALIZE_SPLIT
      printf("splitMxd: event %d, maxDiag %d\n", events[k], maxDiag);
#endif
      continue;
    }

    // Subtract maxDiag from events[k]
    // Do this only for SplitOnly. Other cases are handled later.

    if (split == SplitOnly) {
      mxdDifference->computeTemp(events[k], maxDiag, events[k]);
#ifdef DEBUG_FINALIZE_SPLIT
      printf("SplitOnly: event %d = event %d - maxDiag %d\n",
          events[k], tmp, maxDiag);
#endif
    }

    // Add maxDiag to events[level(maxDiag)]
    int maxDiagLevel = ABS(maxDiag.getLevel());

    mxdUnion->computeTemp(maxDiag, events[maxDiagLevel], events[maxDiagLevel]);

    // Subtract events[maxDiagLevel] from events[k].
    // Do this only for SplitSubtract. SplitSubtractAll is handled later.
    if (split == SplitSubtract) {
      mxdDifference->computeTemp(events[k], events[maxDiagLevel], events[k]);

#ifdef DEBUG_FINALIZE_SPLIT
      printf("SplitSubtract: event %d = event %d - event[maxDiagLevel] %d\n",
          events[k], tmp, maxDiag);
#endif
    }

  } // for k

  if (split == SplitSubtractAll) {
    // Subtract event[i] from all event[j], where j > i.
    for (int i = 1; i < K; i++) {
      if (0==events[i].getNode()) continue;

      for (int j = i + 1; j <= K; j++) {
        if (0==events[j].getNode()) continue;

        mxdDifference->computeTemp(events[j], events[i], events[j]);

#ifdef DEBUG_FINALIZE_SPLIT
        printf("SplitSubtractAll: event %d = event %d - event %d\n",
              events[j], tmp, events[i]);
#endif
      } // for j
    } // for i
  }

#ifdef DEBUG_FINALIZE_SPLIT
  printf("After splitting events in finalize()\n");
  printf("events array: [");
  for (int i=0; i<=K; i++) {
    if (i) printf(", ");
    printf("%d", events[i]);
  }
  printf("]\n");
#endif
}

// HERE!


void MEDDLY::pregen_relation::unionLevels()
{
  if (K < 1) return;

  binary_operation* mxdUnion = UNION(mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdUnion);

  dd_edge u(mxdF);
  for (unsigned k=1; k<=K; k++) {
    apply(UNION, u, events[k], u);
    events[k].set(0);
  }
  events[u.getLevel()] = u;
}


void MEDDLY::pregen_relation::finalize(splittingOption split)
{
  if (0==level_index) {
    // by levels
    switch (split) {
      case SplitOnly: printf("Split: SplitOnly\n"); break;
      case SplitSubtract: printf("Split: SplitSubtract\n"); break;
      case SplitSubtractAll: printf("Split: SplitSubtractAll\n"); break;
      case MonolithicSplit: printf("Split: MonolithicSplit\n"); break;
      default: printf("Split: None\n");
    }
    splitMxd(split);
    if (split != None && split != MonolithicSplit) {
#ifdef DEBUG_FINALIZE_SPLIT
      // Union the elements, and then re-run.
      // Result must be the same as before.
      dd_edge* old_events = new dd_edge[K+1];
      for(unsigned k = 0; k <= K; k++) {
        old_events[k] = events[k];
      }
      splitMxd(MonolithicSplit);
      binary_operation* mxdDifference = DIFFERENCE(mxdF, mxdF, mxdF);
      MEDDLY_DCASSERT(mxdDifference);
      for(unsigned k = 0; k <= K; k++) {
        if (old_events[k] != events[k]) {
          node_handle diff1 = mxdDifference->computeTemp(old_events[k], events[k]);
          node_handle diff2 = mxdDifference->computeTemp(events[k], old_events[k]);
          printf("error at level %d, n:k %d:%d, %d:%d\n",
              k,
              diff1, mxdF->getNodeLevel(diff1),
              diff2, mxdF->getNodeLevel(diff2)
              );
          mxdF->unlinkNode(diff1);
          mxdF->unlinkNode(diff2);
        }
      }

      delete [] old_events;
#endif
    }
    return;
  }

  //
  // Still here?  Must be by events.
  //

#ifdef DEBUG_FINALIZE
  printf("Finalizing pregen relation\n");
  printf("%u events total\n", last_event);
  printf("events array: [");
  for (unsigned i=0; i<last_event; i++) {
    if (i) printf(", ");
    printf("%ld", long(events[i].getNode()));
  }
  printf("]\n");
  printf("next array: [");
  for (unsigned i=0; i<last_event; i++) {
    if (i) printf(", ");
    printf("%d", next[i]);
  }
  printf("]\n");
  printf("level_index array: [%d", level_index[1]);
  for (int i=2; i<=K; i++) {
    printf(", %d", level_index[i]);
  }
  printf("]\n");
#endif

  //
  // Convert from array of linked lists to contiguous array.
  //
  dd_edge* new_events = new dd_edge[last_event];
  unsigned P = 0;
  for (unsigned k=K; k; k--) {
    unsigned L = level_index[k];
    level_index[k] = P;
    while (L) {
      // L+1 is index of an element at level k
      L--;
      new_events[P++] = events[L];
      L = next[L];
    } // while L
  }
  level_index[0] = P;
  delete[] events;
  events = new_events;

  // done with next pointers
  delete[] next;
  next = 0;

#ifdef DEBUG_FINALIZE
  printf("\nAfter finalization\n");
  printf("events array: [");
  for (unsigned i=0; i<last_event; i++) {
    if (i) printf(", ");
    printf("%ld", long(events[i].getNode()));
  }
  printf("]\n");
  printf("level_index array: [%d", level_index[1]);
  for (int i=2; i<=K; i++) {
    printf(", %d", level_index[i]);
  }
  printf("]\n");
#endif
}

// ******************************************************************
// *                                                                *
// *                      otf_subevent methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::otf_subevent::otf_subevent(forest* f, int* v, int nv, bool firing)
: vars(0), num_vars(nv), root(dd_edge(f)), top(0),
  f(f), is_firing(firing)
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

MEDDLY::otf_subevent::~otf_subevent()
{
  if (vars) delete [] vars;
  for (int i=0; i<num_minterms; i++) {
    delete[] unpminterms[i];
    delete[] pminterms[i];
  }
  free(unpminterms);
  free(pminterms);
}

void MEDDLY::otf_subevent::clearMinterms()
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


void MEDDLY::otf_subevent::confirm(otf_relation& rel, int v, int i) {
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}


bool MEDDLY::otf_subevent::addMinterm(const int* from, const int* to)
{
  /*
  ostream_output out(std::cout);
  out << "Adding minterm: [";
  for (int i = f->getMaxLevelIndex(); i >= 0; i--) {
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
  for (int i = f->getMaxLevelIndex(); i >= 0; i--) {
    unpminterms[num_minterms][i] = from[i];
    pminterms[num_minterms][i] = to[i];
    // out << unpminterms[num_minterms][i] << " -> " << pminterms[num_minterms][i] << " , ";
  }
  // out << "]\n";
  domain* d = f->getDomain();
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

void MEDDLY::otf_subevent::buildRoot() {
  if (0 == num_minterms) return;
  /*
  ostream_output out(std::cout);
  out << "Building subevent from " << num_minterms << " minterms\n";
  for (int i = 0; i < num_minterms; i++) {
    out << "minterm[" << i << "]: [ ";
    for (int j = f->getMaxLevelIndex(); j >= 0; j--) {
      out << unpminterms[i][j] << " -> " << pminterms[i][j] << ", ";
    }
    out << " ]\n";
  }
  */
  if (usesExtensibleVariables()) {
    dd_edge sum(root);
    f->createEdge(unpminterms, pminterms, num_minterms, sum);
    num_minterms = 0;
    //
    // root += sum;
    binary_operation* opPlus = UNION(root.getForest(), sum.getForest(), root.getForest());
    MEDDLY_DCASSERT(opPlus);
    opPlus->computeTemp(root, sum, root);
  } else {
    f->createEdge(unpminterms, pminterms, num_minterms, root);
  }
  // out << "Equivalent event: " << root.getNode() << "\n";
  // out << "Result: ";
  // root.showGraph(out);
}


void MEDDLY::otf_subevent::showInfo(output& out) const {
  int num_levels = f->getMaxLevelIndex();
  for (int i = 0; i < num_minterms; i++) {
    out << "minterm[" << i << "]: ";
    for (int lvl = num_levels; lvl > 0; lvl--) {
      out << unpminterms[i][lvl] << " -> " << pminterms[i][lvl] << ", ";
    }
    out << "]\n";
  }
  root.showGraph(out);
}

long MEDDLY::otf_subevent::mintermMemoryUsage() const {
  long n_minterms = 0L;
  for (int i = 0; i < size_minterms; i++) {
    if (unpminterms[i] != 0) n_minterms++;
  }
  return long(n_minterms * 2) * long(f->getMaxLevelIndex()) * long(sizeof(int));
}

// ******************************************************************
// *                                                                *
// *                       otf_event  methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::otf_event::otf_event(otf_subevent** p, int np)
{
  if (p == 0 || np <= 0 || p[0]->getForest() == 0)
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  f = p[0]->getForest();
  for (int i=1; i<np; i++) {
    if (p[i]->getForest() != f) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }

  num_subevents = np;
  subevents = new otf_subevent*[np];
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

MEDDLY::otf_event::~otf_event()
{
  for (int i=0; i<num_subevents; i++) delete subevents[i];
  delete[] subevents;
  delete[] vars;
  delete[] firing_vars;
  delete[] event_mask_from_minterm;
  delete[] event_mask_to_minterm;
}

void MEDDLY::otf_event::buildEventMask()
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
  event_mask.showGraph(out);
#endif
}


bool MEDDLY::otf_event::rebuild()
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
    binary_operation* opAnd = INTERSECTION(
        e.getForest(), subevents[i]->getRoot().getForest(), e.getForest()
    );
    MEDDLY_DCASSERT(opAnd);
    opAnd->compute(e, subevents[i]->getRoot(), e);
    //
    // e *= subevents[i]->getRoot();
  }

  /*
  if (e.getNode() == 0) {
    ostream_output out(std::cout);
    f->useDomain()->showInfo(out);
    out << "subevent: " << event_mask.getNode() << "\n";
    event_mask.showGraph(out);
    for (int i = 0; i < num_subevents; i++) {
    out << "subevent: " << subevents[i]->getRoot().getNode() << "\n";
    subevents[i]->getRoot().showGraph(out);
    }
  }
  ostream_output out(std::cout);
  for (int i = 0; i < num_subevents; i++) {
  out << "subevent: " << subevents[i]->getRoot().getNode() << "\n";
  subevents[i]->getRoot().showGraph(out);
  }
  e.showGraph(out);
  */
  if (e == root) return false;
  root = e;
  return true;
}

void MEDDLY::otf_event::enlargeVariables()
{
  domain* ed = f->getDomain();
  for (int i = 0; i < num_vars; i++) {
    int unprimed = ABS(vars[i]);
    int primed = -unprimed;
    int unprimedSize = f->getLevelSize(unprimed);
    int primedSize = f->getLevelSize(primed);
    if (unprimedSize < primedSize) {
      variable* vh = ed->getVar(unprimed);
      if (vh->isExtensible())
        vh->enlargeBound(false, -primedSize);
      else
        vh->enlargeBound(false, primedSize);
    }
    MEDDLY_DCASSERT(f->getLevelSize(unprimed) == f->getLevelSize(primed));
  }
}


void MEDDLY::otf_event::showInfo(output& out) const {
  for (int i = 0; i < num_subevents; i++) {
    out << "subevent " << i << "\n";
    subevents[i]->showInfo(out);
  }
}

long MEDDLY::otf_event::mintermMemoryUsage() const {
  long usage = 0;
  for (int i = 0; i < num_subevents; i++) {
    usage += subevents[i]->mintermMemoryUsage();
  }
  return usage;
}

// ******************************************************************
// *                                                                *
// *                      otf_relation methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::otf_relation::otf_relation(forest* mxd, forest* mdd,
    otf_event** E, int ne) : mxdF(mxd), resF(mdd)
{
    if (!mxdF) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    if (!resF) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);

    if (resF->getDomain() != mxdF->getDomain()) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

/*
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
    */

  // Check forest types
  if (
    !mxdF->isForRelations()     ||
    resF->isForRelations()   ||
    (mxdF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL) ||
    (resF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
  ) {
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
  }

  // Forests are good; set number of variables
  num_levels = mxdF->getMaxLevelIndex() + 1;

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
  events_by_top_level = new otf_event**[num_levels];
  events_by_level = new otf_event**[num_levels];
  for (int i = 0; i < num_levels; i++) {
    events_by_top_level[i] = num_events_by_top_level[i] > 0
      ? new otf_event*[num_events_by_top_level[i]]: nullptr;
    events_by_level[i] = num_events_by_level[i] > 0
      ? new otf_event*[num_events_by_level[i]]: nullptr;
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
    otf_subevent** se = E[i]->getSubevents();
    for (int j = 0; j < nse; j++) {
      int nVars = se[j]->getNumVars();
      const int* vars = se[j]->getVars();
      for (int k = 0; k < nVars; k++) {
        num_subevents_by_level[vars[k]]++;
      }
    }
  }

  // (2) Allocate subevents[i]
  subevents_by_level = new otf_subevent**[num_levels];
  for (int i = 0; i < num_levels; i++) {
    subevents_by_level[i] = num_subevents_by_level[i] > 0
      ? new otf_subevent*[num_subevents_by_level[i]]: nullptr;
    num_subevents_by_level[i] = 0; // reset this; to be used by the next loop
  }

  // (3) Store subevents by level
  for (int i = 0; i < ne; i++) {
    if (E[i]->isDisabled()) continue;
    int nse = E[i]->getNumOfSubevents();
    otf_subevent** se = E[i]->getSubevents();
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

MEDDLY::otf_relation::~otf_relation()
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

void MEDDLY::otf_relation::clearMinterms()
{
  for (int level = 1; level < num_levels; level++) {
    // Get subevents affected by this level, and rebuild them.
    int nSubevents = num_subevents_by_level[level];
    for (int i = 0; i < nSubevents; i++) {
      subevents_by_level[level][i]->clearMinterms();
    }
  }
}


bool MEDDLY::otf_relation::confirm(int level, int index)
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

void MEDDLY::otf_relation::confirm(const dd_edge& set)
{
  // Perform a depth-first traversal of set:
  //    At each level, mark all enabled states as confirmed.

  // Enlarge the confirmed arrays if needed
  for (int i = 1; i < num_levels; i++) {
      enlargeConfirmedArrays(i, mxdF->getLevelSize(-i));
  }
  std::set<node_handle> visited;
  findConfirmedStates(confirmed, num_confirmed, set.getNode(),
          num_levels-1, visited);

  // ostream_output out(std::cout);
  // showInfo(out);
}

void MEDDLY::otf_relation::findConfirmedStates(bool** confirmed,
    int* num_confirmed, node_handle mdd, int level,
    std::set<MEDDLY::node_handle>& visited)
{
  if (level == 0) return;
  if (visited.find(mdd) != visited.end()) return;

  int mdd_level = resF->getNodeLevel(mdd);
  if (MEDDLY::isLevelAbove(level, mdd_level)) {
    // skipped level; confirm all local states at this level
    // go to the next level
    int level_size = resF->getLevelSize(level);
    for (int i = 0; i < level_size; i++) {
      if (!confirmed[level][i]) {
        confirm(level, i);
      }
    }
    findConfirmedStates(confirmed, num_confirmed, mdd, level-1, visited);
  } else {
    if (MEDDLY::isLevelAbove(mdd_level, level)) {
      throw MEDDLY::error(MEDDLY::error::INVALID_VARIABLE, __FILE__, __LINE__);
    }
    // mdd_level == level
    visited.insert(mdd);
    MEDDLY::unpacked_node *nr = resF->newUnpacked(mdd, MEDDLY::SPARSE_ONLY);
    for (unsigned i = 0; i < nr->getSize(); i++) {
      if (!confirmed[level][nr->index(i)]) {
        confirm(level, int(nr->index(i)));
      }
      findConfirmedStates(confirmed, num_confirmed, nr->down(i), level-1, visited);
    }
    MEDDLY::unpacked_node::Recycle(nr);
  }
}



void MEDDLY::otf_relation::enlargeConfirmedArrays(int level, int sz)
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


void MEDDLY::otf_relation::showInfo(output &strm) const
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
      events_by_top_level[level][ei]->getRoot().showGraph(strm);
    }
  }
}

long MEDDLY::otf_relation::mintermMemoryUsage() const {
  long usage = 0;
  for (int level = 1; level < num_levels; level++) {
    for (int ei = 0; ei < getNumOfEvents(level); ei++) {
      usage += events_by_top_level[level][ei]->mintermMemoryUsage();
    }
  }
  return usage;
}

void MEDDLY::otf_relation::bindExtensibleVariables() {
  //
  // Find the bounds for each extensbile variable
  //
  domain* ed = mxdF->getDomain();
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


double MEDDLY::otf_relation::getArcCount(
  const dd_edge& mask,
  bool count_duplicates)
{
  MEDDLY_DCASSERT(resF->isQuasiReduced());
  MEDDLY_DCASSERT(mxdF->isIdentityReduced());

  double arc_count = 0;
  dd_edge mxd_mask(mxdF);

  // Build confirmed mask
  dd_edge confirmed_local_states(resF);
  terminal one_terminal(true);
  confirmed_local_states.set(one_terminal.getHandle());
  for (int k = 1; k < num_levels; k++) {
    node_handle current_node = confirmed_local_states.getNode();
    int current_level = resF->getNodeLevel(current_node);
    int next_level = resF->upLevel(resF->upLevel(current_level));
    MEDDLY_DCASSERT(next_level >= 0);
    unpacked_node* node = unpacked_node::newFull(resF, next_level,
            unsigned(resF->getLevelSize(next_level))
    );
    for (unsigned i = 0; i < node->getSize(); i++) {
      if (confirmed[next_level][i]) {
          node->setFull(i, resF->linkNode(current_node));
      }
      /*
      node->d_ref(i) =
        confirmed[next_level][i]
        ? resF->linkNode(current_node)
        : 0;
      */
    }
    node_handle next_node = resF->createReducedNode(-1, node);
    confirmed_local_states.set(next_node);
  }

  dd_edge confirmed_local_states_mask(mask);
  apply(INTERSECTION, confirmed_local_states_mask,
        confirmed_local_states, confirmed_local_states_mask);

  apply(CROSS, confirmed_local_states_mask, confirmed_local_states_mask, mxd_mask);

  if (count_duplicates) {
    for (int k = 1; k < num_levels; k++) {
      for (int ei = 0; ei < getNumOfEvents(k); ei++) {
        // start with (num_level-1) to correctly count edges in skipped levels
        dd_edge rg_ei = events_by_top_level[k][ei]->getRoot();
        apply(INTERSECTION, rg_ei, mxd_mask, rg_ei);
        double c;
        apply(CARDINALITY, rg_ei, c);
        arc_count += c;
      }
    }
  } else {
    // build monolithic
    dd_edge monolithic_nsf(mxdF);
    for (int k = 1; k < num_levels; k++) {
      dd_edge nsf_i(mxdF);
      for (int ei = 0; ei < getNumOfEvents(k); ei++) {
        apply(UNION, nsf_i, events_by_top_level[k][ei]->getRoot(), nsf_i);
      }
      apply(UNION, monolithic_nsf, nsf_i, monolithic_nsf);
    }
    apply(INTERSECTION, monolithic_nsf, mxd_mask, monolithic_nsf);
    apply(CARDINALITY, monolithic_nsf, arc_count);
  }
  return arc_count;
}


void MEDDLY::otf_relation::getBoundedMonolithicNSF(dd_edge &root) const
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

MEDDLY::node_handle MEDDLY::otf_relation::getBoundedMxd(
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
  for ( ; i < mxd_node->getSize(); i++) bounded_node->d_ref(i) = mxdF->linkNode(mxd_node->down(i));
  for ( ; i < result_size; i++) bounded_node->d_ref(i) = mxdF->linkNode(ext_d);
  bounded_node->markAsNotExtensible();

  MEDDLY::node_handle result = mxdF->createReducedNode(-1, bounded_node);
  cache[mxd] = mxdF->linkNode(result);
  unpacked_node::Recycle(mxd_node);

  return result;
}

*/


// ******************************************************************
// *                                                                *
// *                   implicit_relation  methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::implicit_relation::implicit_relation(forest* inmdd, forest* relmxd,
                                                             forest* outmdd)
: insetF(inmdd), outsetF(outmdd), mixRelF(relmxd)
{
  mixRelF = relmxd;

  if (0==insetF || 0==outsetF || 0==mixRelF ) throw error(error::MISCELLANEOUS, __FILE__, __LINE__);

  // Check for same domain
  if (insetF->getDomain() != outsetF->getDomain())
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  // for now, anyway, inset and outset must be same forest
  if (insetF != outsetF)
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

  // Check forest types
  if (
      insetF->isForRelations()    ||
      outsetF->isForRelations()   ||
      (insetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)   ||
      (outsetF->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
      )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  // Forests are good; set number of variables
  num_levels = insetF->getMaxLevelIndex();



  //Allocate event_list
  event_list = (rel_node_handle**)malloc(unsigned(num_levels+1)*sizeof(rel_node_handle*));
  event_list_alloc = (long*)malloc(unsigned(num_levels+1)*sizeof(long));
  event_added = (long*)malloc(unsigned(num_levels+1)*sizeof(long));


  confirm_states = (long*)malloc(unsigned(num_levels+1)*sizeof(long));
  confirmed_array_size = (long*)malloc(unsigned(num_levels+1)*sizeof(long));
  confirmed = new bool*[num_levels+1];

  confirmed[0]=0;
  for(int i = 1;i<=num_levels;i++)
    {
    event_list[i] = (rel_node_handle*)malloc(8*sizeof(rel_node_handle));
    confirmed[i] = (bool*)malloc(insetF->getVariableSize(i)*sizeof(bool));
    event_list_alloc[i] = 8;
    event_added[i] = 0;
    confirm_states[i] = 0;

    confirmed_array_size[i]=insetF->getVariableSize(i);
    for(int j = 0;j<insetF->getVariableSize(i);j++)
      confirmed[i][j]=false;
    }



  //create the terminal node
  relation_node *Terminal = new MEDDLY::relation_node(0,mixRelF,0,-1,0,0,-1);
  //mixRelF->createRelationNode(Terminal);
  Terminal->setID(TERMINAL_NODE);
  std::pair<rel_node_handle, relation_node*> TerminalNode(TERMINAL_NODE,Terminal);
  impl_unique.insert(TerminalNode);
  last_in_node_array = TERMINAL_NODE;

}


long
MEDDLY::implicit_relation::getTotalEvent(int level)
{
  int total_event = 0;
  for(int i=1;i<=level;i++)
    total_event +=  lengthForLevel(i);

  return total_event;
}


void
MEDDLY::implicit_relation::resizeEventArray(int level)
{
  event_added[level] += 1;
  if (event_added[level] > event_list_alloc[level]) {
    int nalloc = ((event_added[level]/8)+1)*8;
    MEDDLY_DCASSERT(nalloc > 0);
    MEDDLY_DCASSERT(nalloc > event_added[level]);
    MEDDLY_DCASSERT(nalloc > event_list_alloc[level]);
    event_list[level] = (rel_node_handle*) realloc(event_list[level], nalloc*sizeof(rel_node_handle));
    if (0==event_list[level]) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    event_list_alloc[level] = nalloc;
  }
}

void
MEDDLY::implicit_relation::resizeConfirmedArray(int level,int index)
{
  int nalloc = index+1;
 if(nalloc>confirmed_array_size[level])
    {

       MEDDLY_DCASSERT(nalloc > 0);
       MEDDLY_DCASSERT(confirmed_array_size[level] >= 0);
       if(confirmed_array_size[level]==0)
         {
           confirmed[level] = (bool*)malloc(nalloc*sizeof(bool));
           if (0==confirmed[level]) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
         }
        else
          {
            confirmed[level] = (bool*)realloc(confirmed[level], nalloc*sizeof(bool));
            if (0==confirmed[level]) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
          }

        for(int i = confirmed_array_size[level];i<nalloc;i++)
          confirmed[level][i]=false;

         confirmed_array_size[level]=nalloc;
    }

}

void findConfirmedStatesImpl(MEDDLY::implicit_relation* rel,
                             bool** confirmed, long* confirm_states,
                             MEDDLY::node_handle mdd, int level,
                             std::set<MEDDLY::node_handle>& visited) {
  if (level == 0) return;
  if (visited.find(mdd) != visited.end()) return;

  MEDDLY::forest* insetF = rel->getInForest();
  int mdd_level = insetF->getNodeLevel(mdd);
  if (MEDDLY::isLevelAbove(level, mdd_level)) {
    // skipped level; confirm all local states at this level
    // go to the next level
    int level_size = insetF->getLevelSize(level);
    for (int i = 0; i < level_size; i++) {
      if (!confirmed[level][i]) {
        rel->setConfirmedStates(level, i);
      }
    }
    findConfirmedStatesImpl(rel, confirmed, confirm_states, mdd, level-1, visited);
  } else {
    if (MEDDLY::isLevelAbove(mdd_level, level)) {
      throw MEDDLY::error(MEDDLY::error::INVALID_VARIABLE, __FILE__, __LINE__);
    }
    // mdd_level == level
    visited.insert(mdd);
    MEDDLY::unpacked_node *nr = insetF->newUnpacked(mdd, MEDDLY::SPARSE_ONLY);
    for (int i = 0; i < nr->getSize(); i++) {
      if (!confirmed[level][nr->index(i)]) {
        rel->setConfirmedStates(level, nr->index(i));
      }
      findConfirmedStatesImpl(rel, confirmed, confirm_states, nr->down(i), level-1, visited);
    }
    MEDDLY::unpacked_node::Recycle(nr);
  }
}

void MEDDLY::implicit_relation::setConfirmedStates(const dd_edge& set)
{
  // Perform a depth-first traversal of set:
  //    At each level, mark all enabled states as confirmed.

  // Enlarge the confirmed arrays if needed
  for (int i = 1 ; i<=num_levels; i++)
    {
      int levelSize = getInForest()->getLevelSize(i);
      resizeConfirmedArray(i, levelSize);
    }

    std::set<node_handle> visited;
    findConfirmedStatesImpl(const_cast<implicit_relation*>(this),
                      confirmed, confirm_states, set.getNode(), num_levels, visited);

}



MEDDLY::implicit_relation::~implicit_relation()
{
  last_in_node_array = 0;
  impl_unique.clear();

  for(int i = 0; i <=num_levels; i++) {delete[] event_list[i]; delete[] confirmed[i];}
  delete[] event_list;
  delete[] event_added;
  delete[] event_list_alloc;
  delete[] confirmed;
  delete[] confirm_states;
  delete[] confirmed_array_size;
}


    MEDDLY::rel_node_handle
MEDDLY::implicit_relation::isUniqueNode(relation_node* n)
{
  bool is_unique_node = true;
  std::unordered_map<rel_node_handle, relation_node*>::iterator it = impl_unique.begin();
  while(it != impl_unique.end())
    {
    is_unique_node = !((it->second)->equals(n));
    if(is_unique_node==false)
      return (it->second)->getID();
    ++it;
    }
  return 0;
}



    MEDDLY::rel_node_handle
MEDDLY::implicit_relation::registerNode(bool is_event_top, relation_node* n)
{

  rel_node_handle nLevel = n->getLevel();

#ifdef DEVELOPMENT_CODE
  rel_node_handle downHandle = n->getDown();
  relation_node* downNode = nodeExists(downHandle);
  rel_node_handle downLevel = downNode->getLevel();
  MEDDLY_DCASSERT( ( ( downNode!=NULL ) && ( nLevel > downLevel ) )
                    ||
                   ( downLevel == 0 ) );
#endif

  rel_node_handle n_ID = isUniqueNode(n);

  if(n_ID==0) // Add new node
   {
    n_ID  = last_in_node_array + 1;
    std::pair<rel_node_handle, relation_node*> add_node(n_ID,n);
    impl_unique.insert(add_node);
    if(impl_unique.find(n_ID) != impl_unique.end())
    {
      last_in_node_array = n_ID;
      n->setID(n_ID);
    }
    //mixRelF->createRelationNode(n);
  }
  else //Delete the node
    {
     delete n;
    }

  if(is_event_top)
    {
    resizeEventArray(nLevel);
    event_list[nLevel][event_added[nLevel] - 1] = n_ID;
    }

  return n_ID;
}

void
MEDDLY::implicit_relation::show()
{
  rel_node_handle** event_list_copy = (rel_node_handle**)malloc((num_levels+1)*sizeof(rel_node_handle*));
  if (0==event_list_copy) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  long total_events = 0;
  for(int i = 1;i<=num_levels;i++) total_events +=event_added[i];
  for(int i = 1;i<=num_levels;i++)
    {
     event_list_copy[i] = (rel_node_handle*)malloc(total_events*sizeof(rel_node_handle));
     if (0==event_list_copy[i]) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }

  for(int i = num_levels;i>=1;i--)
    for(int j=0;j<total_events;j++)
      event_list_copy[i][j]=0;


  int eid = 0;
  for(int i = num_levels;i>=1;i--)
    {
     int k = 0;
     std::cout<<"\n [";
     for(int j=0;j<total_events;j++)
      {

        if((event_list_copy[i][eid]==0)&&(k<event_added[i]))
          {
            event_list_copy[i][eid] = event_list[i][k];
            relation_node* hold_it = nodeExists(event_list[i][k]);
            relation_node* hold_down = nodeExists(hold_it->getDown());
            event_list_copy[hold_down->getLevel()][eid] = hold_down->getID();
          k++;eid++;
          }


      int dig_ctr = event_list_copy[i][j]>1000?4:(event_list_copy[i][j]>100?3:(event_list_copy[i][j]>10?2:1));

      int spc_bef =(6 - dig_ctr)/2;
      int spc_aft = 6 - dig_ctr - spc_bef;

      for(int s=0;s<spc_bef;s++) std::cout<<" ";
      if(event_list_copy[i][j] != 0) std::cout<<event_list_copy[i][j];
      else std::cout<<"_";
      for(int s=0;s<spc_aft;s++) std::cout<<" ";
      if(j!=total_events-1)
          std::cout<<"|";
      }
     std::cout<<"]";
    }

  for(int i = 0;i<num_levels+1;i++) delete event_list_copy[i];
  delete[] event_list_copy;

}

void MEDDLY::implicit_relation::bindExtensibleVariables() {
  //
  // Find the bounds for each extensbile variable
  //
  domain* ed = outsetF->getDomain();

  for (int k = 1; k <= num_levels; k++) {
    int bound = 0;
    int n_confirmed = 0;

    for (int j = 0; j < confirmed_array_size[k]; j++) {
      if (confirmed[k][j]) { bound = j+1; n_confirmed++; }
    }

    MEDDLY_DCASSERT(bound > 0);
    MEDDLY_DCASSERT(n_confirmed == confirm_states[k]);
    ed->enlargeVariableBound(k, false, bound);
  }
}

MEDDLY::node_handle
MEDDLY::implicit_relation::buildMxdForest()
{

  //Get number of Variables and Events
  int nVars = outsetF->getMaxLevelIndex();
  int nEvents = getTotalEvent(nVars);


  rel_node_handle* event_tops = (rel_node_handle*)malloc((nEvents)*sizeof(rel_node_handle));
  int e = 0;

  for(int i = 1 ;i<=nVars;i++)
    {
    int num_events_at_this_level = lengthForLevel(i);
    for(int j = 0;j<num_events_at_this_level;j++)
      event_tops[e++]=arrayForLevel(i)[j];
    }

  domain *d = outsetF->getDomain();

  forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
                                edge_labeling::MULTI_TERMINAL);
  dd_edge nsf(mxd);

  dd_edge* monolithic_nsf = new dd_edge(mxd);

    binary_operation* opUnion = UNION(mxd, mxd, mxd);
    for(int i=0;i<nEvents;i++) {
        opUnion->computeTemp(
            *monolithic_nsf,  buildEventMxd(event_tops[i],mxd),
            *monolithic_nsf
        );
        // OLD
        // (*monolithic_nsf) += buildEventMxd(event_tops[i],mxd);
    }

  node_handle monolithic_nsf_handle = monolithic_nsf->getNode();
  mxdF = (forest*)mxd;

  /*for(int i = 0; i<nEvents;i++)
   {
   dd_edge nsf_ev(mxd);
   nsf_ev = buildEventMxd(event_tops[i],mxd);
   apply(UNION, nsf, nsf_ev, nsf);
   }*/

  return monolithic_nsf_handle;
}


MEDDLY::dd_edge
MEDDLY::implicit_relation::buildEventMxd(rel_node_handle eventTop, forest *mxd)
{
  //mxd is built on a domain obtained from result of saturation
  int nVars = outsetF->getMaxLevelIndex();
  //int* sizes = new int[nVars];
  relation_node* Rnode = nodeExists(eventTop);
  rel_node_handle* rnh_array = (rel_node_handle*)malloc((nVars+1)*sizeof(rel_node_handle));
  // int top_level = Rnode->getLevel();

  forest* ef = (forest*) mxd;

  //Get relation node handles
  for (int i=nVars; i>=1; i--)
    {
      if(Rnode->getLevel()==i)// if variable i is a part of this event
        {
          rnh_array[i] = Rnode->getID(); // keep track of node_handles that are part of this event
          Rnode = nodeExists(Rnode->getDown()); // move to next variable in the event
        }
      else // if not, then
        {
        rnh_array[i] = -1; // node handle of the variable i in the event
        continue;
        }
    }

  node_handle below = ef->handleForValue(true); // Terminal true node

  for (int i=1; (i<=nVars)&&(below!=0); i++)
    {
        if(rnh_array[i]!=-1)
          {
            const int maxVar = ef->getDomain()->getVariableBound(i);
            Rnode = nodeExists(rnh_array[i]);
            //Create a new unprimed node for variable i
            MEDDLY_DCASSERT(outsetF->getVariableSize(i)>=Rnode->getPieceSize());
            unpacked_node* UP_var = unpacked_node::newFull(ef, i, Rnode->getPieceSize());

            for (int j=0; j<Rnode->getPieceSize(); j++) {

              // long new_j = confirmed[i][Rnode->getTokenUpdate()[j]]? Rnode->getTokenUpdate()[j] : -2;
              long new_j = Rnode->getTokenUpdate()[j];

              if(new_j>=0 && new_j<maxVar) // do not exceed variable bounds
                {
                   //Create primed node for each valid index of the unprimed node
                  unpacked_node* P_var = unpacked_node::newSparse(ef, -i, 1);
                  P_var->setSparse(0, new_j, ef->linkNode(below));

                  // P_var->i_ref(0) = new_j;
                  // P_var->d_ref(0) = ef->linkNode(below); // primed node for new_j index points to terminal or unprime node
                  UP_var->setFull(j, ef->createReducedNode(j, P_var));
                  // UP_var->d_ref(j) = ef->createReducedNode(j, P_var);
                }
              else
                UP_var->setFull(j, ef->handleForValue(false));
                // UP_var->d_ref(j) = ef->handleForValue(false); // unprimed node for j index points to false
              }

              ef->unlinkNode(below);
              below = ef->createReducedNode(-1, UP_var);
          }
    }

  dd_edge nsf(mxd);
  nsf.set(below);

  return nsf;
}

// ******************************************************************


std::unordered_map<long,std::vector<MEDDLY::rel_node_handle>>
MEDDLY::implicit_relation::getListOfNexts(int level, long i, relation_node **R)
{
  std::unordered_map<long,std::vector<rel_node_handle>> jList;
  // atleast as many j's as many events
  for(int k=0;k<lengthForLevel(level);k++)
    {
    long key = R[k]->nextOf(i);
    jList[key].reserve(lengthForLevel(level));
    int rnh_dwn = R[k]->getDown();
    jList[key].push_back(rnh_dwn);
    }

  return jList;
}

bool
MEDDLY::implicit_relation::isUnionPossible(int level, long i, relation_node **R)
{
  if(lengthForLevel(level)==1)
     return false;

  std::vector<int> jset(lengthForLevel(level), 0);

   int last_j = 0;
   for(int k=0;k<lengthForLevel(level);k++)
    {
    long key = R[k]->nextOf(i);
    int flag = 0;
    for(int m=0;m<last_j;m++)
      if(jset[m]==key)
        {
          flag=1;
          break;
        }

      if(flag==0)
        {
          jset[k]=key;
          last_j++;
        }
    }
  if(lengthForLevel(level)==last_j)
   return false;
  else
    return true;
}


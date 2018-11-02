
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
    Implementation of operation framework.

    Actual operations are in separate files, in operations/ directory.
*/

#include "defines.h"
#include <set>
// #include "compute_table.h"

   #define OUT_OF_BOUNDS -1
   #define NOT_KNOWN -2
   #define TERMINAL_NODE 1

// #define DEBUG_CLEANUP
// #define DEBUG_FINALIZE
// #define DEBUG_FINALIZE_SPLIT
// #define DEBUG_EVENT_MASK

namespace MEDDLY {

  /*
      Convert from array of linked lists to contiguous array with indexes.
  
      Idea:
      traverse the array in order.
      everything after this point is still "linked lists".
      everything before this point is contiguous, but
        it is possible that one of the lists points here.
        if so, the next pointer holds the "forwarding address".
  */
  template <typename TYPE>
  void defragLists(TYPE* data, int* next, int* Lists, int NL)
  {
    int P = 0;  // "this point"
    for (int i=NL; i; i--) {
      int L = Lists[i];
      Lists[i] = P;
      while (L>=0) {
        // if L<P then there must be a forwarding address; follow it
        while (L<P) L = next[L];
        // ok, we're at the right slot now
        int nxt = next[L];
        if (L != P) {
          // element L belongs in slot P, swap them...
          SWAP(data[L], data[P]);
          next[L] = next[P];
          // ...and set up the forwarding address
          next[P] = L;
        }
        P++;
        L = nxt;
      }
    } // for k
    Lists[0] = P;
  }
}


// ******************************************************************
// *                         opname methods                         *
// ******************************************************************

int MEDDLY::opname::next_index;

MEDDLY::opname::opname(const char* n)
{
  name = n;
  index = next_index;
  next_index++;
}

MEDDLY::opname::~opname()
{
  // library must be closing
}

// ******************************************************************
// *                      unary_opname methods                      *
// ******************************************************************

MEDDLY::unary_opname::unary_opname(const char* n) : opname(n)
{
}

MEDDLY::unary_opname::~unary_opname()
{
}

MEDDLY::unary_operation* 
MEDDLY::unary_opname::buildOperation(expert_forest* ar, expert_forest* rs) const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);  
}

MEDDLY::unary_operation* 
MEDDLY::unary_opname::buildOperation(expert_forest* ar, opnd_type res) const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}


// ******************************************************************
// *                     binary_opname  methods                     *
// ******************************************************************

MEDDLY::binary_opname::binary_opname(const char* n) : opname(n)
{
}

MEDDLY::binary_opname::~binary_opname()
{
}

// ******************************************************************
// *                   specialized_opname methods                   *
// ******************************************************************

MEDDLY::specialized_opname::arguments::arguments()
{
  setAutoDestroy(true);
}

MEDDLY::specialized_opname::arguments::~arguments()
{
}


MEDDLY::specialized_opname::specialized_opname(const char* n) : opname(n)
{
}

MEDDLY::specialized_opname::~specialized_opname()
{
}

// ******************************************************************
// *                    numerical_opname methods                    *
// ******************************************************************

MEDDLY::numerical_opname::numerical_args
::numerical_args(const dd_edge &xi, const dd_edge &a, const dd_edge &yi)
 : x_ind(xi), A(a), y_ind(yi)
{
}

MEDDLY::numerical_opname::numerical_args::~numerical_args()
{
}


MEDDLY::numerical_opname::numerical_opname(const char* n)
 : specialized_opname(n)
{
}

MEDDLY::numerical_opname::~numerical_opname()
{
}

// ******************************************************************
// *                                                                *
// *                    satpregen_opname methods                    *
// *                                                                *
// ******************************************************************


MEDDLY::satpregen_opname::satpregen_opname(const char* n)
 : specialized_opname(n)
{
}

MEDDLY::satpregen_opname::~satpregen_opname()
{
}

void
MEDDLY::satpregen_opname::pregen_relation
::setForests(forest* inf, forest* mxd, forest* outf)
{
  insetF = inf;
  outsetF = outf;
  mxdF = smart_cast <MEDDLY::expert_forest*>(mxd);
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
    (insetF->getEdgeLabeling() != forest::MULTI_TERMINAL)   ||
    (outsetF->getEdgeLabeling() != forest::MULTI_TERMINAL)  ||
    (mxdF->getEdgeLabeling() != forest::MULTI_TERMINAL)     ||
    (outsetF->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  // Forests are good; set number of variables
  K = mxdF->getDomain()->getNumVariables();
}

MEDDLY::satpregen_opname::pregen_relation
::pregen_relation(forest* inf, forest* mxd, forest* outf, int nevents)
{
  setForests(inf, mxd, outf);

  num_events = nevents;
  if (num_events) {
    events = new node_handle[num_events];
    next = new int[num_events];
  } else {
    events = 0;
    next = 0;
  }
  last_event = -1;

  level_index = new int[K+1];
  for (int k=0; k<=K; k++) level_index[k] = -1;   // null pointer
}

MEDDLY::satpregen_opname::pregen_relation
::pregen_relation(forest* inf, forest* mxd, forest* outf)
{
  setForests(inf, mxd, outf);

  events = new node_handle[K+1];
  for (int k=0; k<=K; k++) events[k] = 0;

  next = 0;
  level_index = 0;

  num_events = -1;
  last_event = -1;
}


MEDDLY::satpregen_opname::pregen_relation
::~pregen_relation()
{
  if (events && mxdF) {
    for (int k=0; k<=last_event; k++) if(events[k]) mxdF->unlinkNode(events[k]);
  }
  delete[] events;
  delete[] next;
  delete[] level_index;
}

void
MEDDLY::satpregen_opname::pregen_relation
::addToRelation(const dd_edge &r)
{
  MEDDLY_DCASSERT(mxdF);

  if (r.getForest() != mxdF)  throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

  int k = r.getLevel(); 
  if (0==k) return;
  if (k<0) k = -k;   

  if (0==level_index) {
    // relation is "by levels"

    if (0==events[k]) {
      events[k] = mxdF->linkNode(r.getNode());
    } else {
      // already have something at this level; perform the union
      dd_edge tmp(mxdF);
      tmp.set(events[k]); // does not increase incoming count
      tmp += r;
      events[k] = mxdF->linkNode(tmp.getNode());
      // tmp should be destroyed here
    }

  } else {
    // relation is "by events"

    if (isFinalized())              throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    if (last_event+1 >= num_events) throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);

    last_event++;

    events[last_event] = mxdF->linkNode(r.getNode());
    next[last_event] = level_index[k];
    level_index[k] = last_event;
  }
}


void
MEDDLY::satpregen_opname::pregen_relation
::splitMxd(splittingOption split)
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
  binary_operation* mxdUnion = getOperation(UNION, mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdUnion);

  binary_operation* mxdIntersection =
    getOperation(INTERSECTION, mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdIntersection);

  binary_operation* mxdDifference = getOperation(DIFFERENCE, mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdDifference);

  for (int k = K; k > 1; k--) {
    if (0 == events[k]) continue;

    MEDDLY_DCASSERT(ABS(mxdF->getNodeLevel(events[k])) <= k);

    // Initialize unpacked nodes
    unpacked_node* Mu = (isLevelAbove(k, mxdF->getNodeLevel(events[k]))) 
      ?   unpacked_node::newRedundant(mxdF, k, events[k], true)
      :   unpacked_node::newFromNode(mxdF, events[k], true)
    ;

    unpacked_node* Mp = unpacked_node::useUnpackedNode();

    bool first = true;
    node_handle maxDiag = 0;

    // Read "rows"
    for (int i = 0; i < Mu->getSize(); i++) {
      // Initialize column reader
      if (isLevelAbove(-k, mxdF->getNodeLevel(Mu->d(i)))) {
        Mp->initIdentity(mxdF, -k, i, Mu->d(i), true);
      } else {
        Mp->initFromNode(mxdF, Mu->d(i), true);
      }

      // Intersect along the diagonal
      if (first) {
        maxDiag = mxdF->linkNode(Mp->d(i));
        first = false;
      } else {
        node_handle nmd = mxdIntersection->compute(maxDiag, Mp->d(i));
        mxdF->unlinkNode(maxDiag);
        maxDiag = nmd;
      }
    } // for i

    // Cleanup
    unpacked_node::recycle(Mp);
    unpacked_node::recycle(Mu);

    if (0 == maxDiag) {
#ifdef DEBUG_FINALIZE_SPLIT
      printf("splitMxd: event %d, maxDiag %d\n", events[k], maxDiag);
#endif
      continue;
    }

    // Subtract maxDiag from events[k]
    // Do this only for SplitOnly. Other cases are handled later.
    if (split == SplitOnly) {
      int tmp = events[k];
      events[k] = mxdDifference->compute(events[k], maxDiag);
      mxdF->unlinkNode(tmp);
#ifdef DEBUG_FINALIZE_SPLIT
      printf("SplitOnly: event %d = event %d - maxDiag %d\n",
          events[k], tmp, maxDiag);
#endif
    } 

    // Add maxDiag to events[level(maxDiag)]
    int maxDiagLevel = ABS(mxdF->getNodeLevel(maxDiag));
    int tmp = events[maxDiagLevel];
    events[maxDiagLevel] = mxdUnion->compute(maxDiag, events[maxDiagLevel]);
    mxdF->unlinkNode(tmp);
    mxdF->unlinkNode(maxDiag);

    // Subtract events[maxDiagLevel] from events[k].
    // Do this only for SplitSubtract. SplitSubtractAll is handled later.
    if (split == SplitSubtract) {
      int tmp = events[k];
      events[k] = mxdDifference->compute(events[k], events[maxDiagLevel]);
      mxdF->unlinkNode(tmp);
#ifdef DEBUG_FINALIZE_SPLIT
      printf("SplitSubtract: event %d = event %d - event[maxDiagLevel] %d\n",
          events[k], tmp, maxDiag);
#endif
    }

  } // for k

  if (split == SplitSubtractAll) {
    // Subtract event[i] from all event[j], where j > i.
    for (int i = 1; i < K; i++) {
      for (int j = i + 1; j <= K; j++) {
        if (events[i] && events[j]) {
          int tmp = events[j];
          events[j] = mxdDifference->compute(events[j], events[i]);
          mxdF->unlinkNode(tmp);
#ifdef DEBUG_FINALIZE_SPLIT
          printf("SplitSubtractAll: event %d = event %d - event %d\n",
              events[j], tmp, events[i]);
#endif
        }
      }
    }
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


void
MEDDLY::satpregen_opname::pregen_relation
::unionLevels()
{
  if (K < 1) return;

  binary_operation* mxdUnion = getOperation(UNION, mxdF, mxdF, mxdF);
  MEDDLY_DCASSERT(mxdUnion);

  node_handle u = events[1];
  events[1] = 0;

  for (int k = 2; k <= K; k++) {
    if (0 == events[k]) continue;
    node_handle temp = mxdUnion->compute(u, events[k]);
    mxdF->unlinkNode(events[k]);
    events[k] = 0;
    mxdF->unlinkNode(u);
    u = temp;
  }

  events[mxdF->getNodeLevel(u)] = u;
}


void
MEDDLY::satpregen_opname::pregen_relation
::finalize(splittingOption split)
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
      node_handle* old_events = new node_handle[K+1];
      for(int k = 0; k <= K; k++) {
        old_events[k] = events[k];
        if (old_events[k]) mxdF->linkNode(old_events[k]);
      }
      splitMxd(MonolithicSplit);
      binary_operation* mxdDifference =
        getOperation(DIFFERENCE, mxdF, mxdF, mxdF);
      MEDDLY_DCASSERT(mxdDifference);
      for(int k = 0; k <= K; k++) {
        if (old_events[k] != events[k]) {
          node_handle diff1 = mxdDifference->compute(old_events[k], events[k]);
          node_handle diff2 = mxdDifference->compute(events[k], old_events[k]);
          printf("error at level %d, n:k %d:%d, %d:%d\n",
              k,
              diff1, mxdF->getNodeLevel(diff1),
              diff2, mxdF->getNodeLevel(diff2)
              );
          mxdF->unlinkNode(diff1);
          mxdF->unlinkNode(diff2);
        }
        if (old_events[k]) mxdF->unlinkNode(old_events[k]);
      }

      delete [] old_events;
#endif
    }
    return; 
  }

#ifdef DEBUG_FINALIZE
  printf("Finalizing pregen relation\n");
  printf("%d events total\n", last_event+1);
  printf("events array: [");
  for (int i=0; i<=last_event; i++) {
    if (i) printf(", ");
    printf("%d", events[i]);
  }
  printf("]\n");
  printf("next array: [");
  for (int i=0; i<=last_event; i++) {
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

  // convert from array of linked lists to contiguous array.
  defragLists(events, next, level_index, K);

  // done with next pointers
  delete[] next;
  next = 0;

#ifdef DEBUG_FINALIZE
  printf("\nAfter finalization\n");
  printf("events array: [");
  for (int i=0; i<=last_event; i++) {
    if (i) printf(", ");
    printf("%d", events[i]);
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
  memcpy(vars, v, num_vars * sizeof(int));

  // find top
  top = vars[0];
  for (int i = 1; i < num_vars; i++) {
    if (isLevelAbove(top, vars[i])) top = vars[i];
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
    unpminterms = (int**) realloc(unpminterms, size_minterms * sizeof(int**));
    pminterms = (int**) realloc(pminterms, size_minterms * sizeof(int**));
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
  return long(n_minterms * 2) * long(f->getDomain()->getNumVariables()) * sizeof(int);
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
    const size_t minterm_size = f->getNumVariables()+1;
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
      expert_variable* vh = ed->getExpertVar(unprimed);
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
    (insetF->getEdgeLabeling() != forest::MULTI_TERMINAL)   ||
    (mxdF->getEdgeLabeling() != forest::MULTI_TERMINAL)     ||
    (outsetF->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  // Forests are good; set number of variables
  num_levels = mxdF->getDomain()->getNumVariables() + 1;

  // Build the events-per-level data structure
  // (0) Initialize
  num_events_by_top_level = new int[num_levels];
  num_events_by_level = new int[num_levels];
  memset(num_events_by_top_level, 0, sizeof(int)*num_levels);
  memset(num_events_by_level, 0, sizeof(int)*num_levels);

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
  memset(num_subevents_by_level, 0, sizeof(int)*num_levels);

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
    confirmed[i] = (bool*) malloc(level_size * sizeof(bool));
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
    MEDDLY::unpacked_node *nr = MEDDLY::unpacked_node::newFromNode(insetF, mdd, false);
    for (int i = 0; i < nr->getNNZs(); i++) {
      if (!confirmed[level][nr->i(i)]) {
        rel->confirm(level, nr->i(i));
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
  confirmed[level] = (bool*) realloc(confirmed[level], sz * sizeof(bool));
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
  confirmed_local_states.set(MEDDLY::expert_forest::bool_Tencoder::value2handle(true));
  for (int i = 1; i < num_levels; i++) {
    node_handle current_node = confirmed_local_states.getNode();
    int current_level = outsetF->getNodeLevel(current_node);
    int next_level = outsetF->upLevel(outsetF->upLevel(current_level));
    MEDDLY_DCASSERT(next_level >= 0);
    unpacked_node* node =
      unpacked_node::newFull(outsetF, next_level, outsetF->getLevelSize(next_level));
    for (int i = 0; i < node->getSize(); i++) {
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
    for (int i = 1; i < num_levels; i++) {
      for (int ei = 0; ei < getNumOfEvents(i); ei++) {
        // start with (num_level-1) to correctly count edges in skipped levels
        dd_edge rg_ei = events_by_top_level[i][ei]->getRoot();
        rg_ei *= mxd_mask;
        arc_count += rg_ei.getCardinality();
      }
    }
  } else {
    // build monolithic 
    dd_edge monolithic_nsf(mxdF);
    for (int i = 1; i < num_levels; i++) {
      dd_edge nsf_i(mxdF);
      for (int ei = 0; ei < getNumOfEvents(i); ei++) {
        nsf_i += events_by_top_level[i][ei]->getRoot();
      }
      monolithic_nsf += nsf_i;
    }
    monolithic_nsf *= mxd_mask;
    arc_count = monolithic_nsf.getCardinality();
  }
  return arc_count;
}


MEDDLY::node_handle MEDDLY::satotf_opname::otf_relation::getBoundedMonolithicNSF() {
  //
  // Build the union all events
  //
  dd_edge* monolithic_nsf = new dd_edge(mxdF);
  for (int i = 1; i < num_levels; i++) {
    for (int ei = 0; ei < getNumOfEvents(i); ei++) {
      (*monolithic_nsf) += events_by_top_level[i][ei]->getRoot();
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
}

MEDDLY::node_handle MEDDLY::satotf_opname::otf_relation::getBoundedMxd(
    MEDDLY::node_handle mxd,
    const int* bounds,
    int num_levels,
    std::unordered_map<MEDDLY::node_handle, MEDDLY::node_handle>& cache
    ) {
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
  unpacked_node *mxd_node = unpacked_node::newFromNode(mxdF, mxd, true);
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

// ******************************************************************
// *                                                                *
// *                 minimum_witness_opname methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_opname::constrained_opname(const char* n)
  : specialized_opname(n)
{
}

// ******************************************************************
// *                                                                *
// *                     satimpl_opname  methods                    *
// *                                                                *
// ******************************************************************



MEDDLY::satimpl_opname::satimpl_opname(const char* n)
: specialized_opname(n)
{
}

MEDDLY::satimpl_opname::~satimpl_opname()
{
}


MEDDLY::relation_node::relation_node(unsigned long sign, int lvl, rel_node_handle d)
{
  signature  = sign;
  level = lvl;
  down = d;
  piece_size = 0;
  token_update = NULL;
}

MEDDLY::relation_node::~relation_node()
{
}


long MEDDLY::relation_node::nextOf(long i)
{
  //to be defined for the example you use & comment this definition
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

bool
MEDDLY::relation_node::equals(const relation_node* n) const
{
  if((signature == n->getSignature()) && (level == n->getLevel()) && (down == n->getDown()))
    return true;
  else
    return false;
}

void
MEDDLY::relation_node::expandTokenUpdate(long i)
{
  if(getPieceSize()==0)
  {
    token_update = (long*)malloc(1*sizeof(long));
    piece_size = 1;
    token_update[0]=NOT_KNOWN;
  }
  if(i>0)
  {
    token_update = (long*)realloc(token_update,(i+1)*sizeof(long));
    for(int j = piece_size;j<=i;j++)
      token_update[j]=NOT_KNOWN;
    piece_size = i+1;
  }
}

void
MEDDLY::relation_node::setTokenUpdateAtIndex(long i,long val)
{
  MEDDLY_DCASSERT(i<getPieceSize());
  token_update[i] = val;
}
// ******************************************************************

MEDDLY::satimpl_opname::implicit_relation::implicit_relation(forest* inmdd, forest* relmxd,
                                                             forest* outmdd)
: insetF(static_cast<expert_forest*>(inmdd)), outsetF(static_cast<expert_forest*>(outmdd)), mixRelF(static_cast<expert_forest*>(relmxd))
{
  
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
      (insetF->getEdgeLabeling() != forest::MULTI_TERMINAL)   ||
      (outsetF->getEdgeLabeling() != forest::MULTI_TERMINAL)
      )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
  
  // Forests are good; set number of variables
  num_levels = insetF->getDomain()->getNumVariables();
  
  
  
  //Allocate event_list
  event_list = (rel_node_handle**)malloc((num_levels+1)*sizeof(rel_node_handle*));
  event_list_alloc = (long*)malloc((num_levels+1)*sizeof(long));
  event_added = (long*)malloc((num_levels+1)*sizeof(long));

  
  confirm_states = (long*)malloc((num_levels+1)*sizeof(long));
  confirmed_array_size = (long*)malloc((num_levels+1)*sizeof(long));
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
  relation_node *Terminal = new relation_node(0,0,TERMINAL_NODE);
  //mixRelF->createRelationNode(Terminal);
  Terminal->setID(TERMINAL_NODE);
  std::pair<rel_node_handle, relation_node*> TerminalNode(TERMINAL_NODE,Terminal);
  impl_unique.insert(TerminalNode);
  last_in_node_array = TERMINAL_NODE;
  
}


void
MEDDLY::satimpl_opname::implicit_relation::resizeEventArray(int level)
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
MEDDLY::satimpl_opname::implicit_relation::resizeConfirmedArray(int level,int index)
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

void findConfirmedStatesImpl(MEDDLY::satimpl_opname::implicit_relation* rel,
                             bool** confirmed, long* confirm_states,
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
    MEDDLY::unpacked_node *nr = MEDDLY::unpacked_node::newFromNode(insetF, mdd, false);
    for (int i = 0; i < nr->getNNZs(); i++) {
      if (!confirmed[level][nr->i(i)]) {
        rel->setConfirmedStates(level, nr->i(i));
      }
      findConfirmedStatesImpl(rel, confirmed, confirm_states, nr->d(i), level-1, visited);
    }
    MEDDLY::unpacked_node::recycle(nr);
  }
}

void MEDDLY::satimpl_opname::implicit_relation::setConfirmedStates(const dd_edge& set)
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



MEDDLY::satimpl_opname::implicit_relation::~implicit_relation()
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


rel_node_handle
MEDDLY::satimpl_opname::implicit_relation::isUniqueNode(relation_node* n)
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

rel_node_handle
MEDDLY::satimpl_opname::implicit_relation::registerNode(bool is_event_top, relation_node* n)
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
    mixRelF->createRelationNode(n);
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
MEDDLY::satimpl_opname::implicit_relation::show()
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

void MEDDLY::satimpl_opname::implicit_relation::bindExtensibleVariables() {
  //
  // Find the bounds for each extensbile variable
  //
  expert_domain* ed = static_cast<expert_domain*>(outsetF->useDomain());
  
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
MEDDLY::satimpl_opname::implicit_relation::buildMxdForest()
{
  
  //Get number of Variables and Events
  int nVars = outsetF->getDomain()->getNumVariables();
  int nEvents = getTotalEvent(nVars);
  
  
  rel_node_handle* event_tops = (rel_node_handle*)malloc((nEvents)*sizeof(rel_node_handle));
  int e = 0;
  
  for(int i = 1 ;i<=nVars;i++)
    {
    int num_events_at_this_level = lengthForLevel(i);
    for(int j = 0;j<num_events_at_this_level;j++)
      event_tops[e++]=arrayForLevel(i)[j];
    }
  
  domain *d = outsetF->useDomain();
  
  forest* mxd = d->createForest(true,forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge nsf(mxd);
  
  dd_edge* monolithic_nsf = new dd_edge(mxd);
  for(int i=0;i<nEvents;i++)
    {
    (*monolithic_nsf) += buildEventMxd(event_tops[i],mxd);
    }
  
  node_handle monolithic_nsf_handle = monolithic_nsf->getNode();
  mxdF = (expert_forest*)mxd;
  
  /*for(int i = 0; i<nEvents;i++)
   {
   dd_edge nsf_ev(mxd);
   nsf_ev = buildEventMxd(event_tops[i],mxd);
   apply(UNION, nsf, nsf_ev, nsf);
   }*/
  
  return monolithic_nsf_handle;
}


MEDDLY::dd_edge
MEDDLY::satimpl_opname::implicit_relation::buildEventMxd(rel_node_handle eventTop, forest *mxd)
{
  //mxd is built on a domain obtained from result of saturation
  int nVars = outsetF->getDomain()->getNumVariables();
  //int* sizes = new int[nVars];
  relation_node* Rnode = nodeExists(eventTop);
  rel_node_handle* rnh_array = (rel_node_handle*)malloc((nVars+1)*sizeof(rel_node_handle));
  // int top_level = Rnode->getLevel();
  
  // domain* d = outsetF->useDomain();
  expert_forest* ef = (expert_forest*) mxd;
  
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
            Rnode = nodeExists(rnh_array[i]);
            //Create a new unprimed node for variable i
            MEDDLY_DCASSERT(outsetF->getVariableSize(i)>=Rnode->getPieceSize());
            unpacked_node* UP_var = unpacked_node::newFull(ef, i, Rnode->getPieceSize());
          
            for (int j=0; j<Rnode->getPieceSize(); j++) {
          
              long new_j = confirmed[i][Rnode->getTokenUpdate()[j]]? Rnode->getTokenUpdate()[j] : -2;
              
              if(new_j>=0) 
                {
                   //Create primed node for each valid index of the unprimed node
                  unpacked_node* P_var = unpacked_node::newSparse(ef, -i, 1);
                  P_var->i_ref(0) = new_j;
                  P_var->d_ref(0) = ef->linkNode(below); // primed node for new_j index points to terminal or unprime node
                  UP_var->d_ref(j) = ef->createReducedNode(j, P_var);
                }
              else
                UP_var->d_ref(j) = ef->handleForValue(false); // unprimed node for j index points to false
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


std::unordered_map<long,std::vector<rel_node_handle>>
MEDDLY::satimpl_opname::implicit_relation::getListOfNexts(int level, long i, relation_node **R)
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
MEDDLY::satimpl_opname::implicit_relation::isUnionPossible(int level, long i, relation_node **R)
{
  if(lengthForLevel(level)==1)
     return false;
  
   int* jset = (int*)malloc(lengthForLevel(level)*sizeof(int));
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

// ******************************************************************
// *                       operation  methods                       *
// ******************************************************************

#ifdef OLD_OP_CT

MEDDLY::operation::operation(const opname* n, int kl, int al)
{
#ifdef DEBUG_CLEANUP
  fprintf(stdout, "Creating operation %p\n", this);
  fflush(stdout);
#endif
  MEDDLY_DCASSERT(kl>=0);
  MEDDLY_DCASSERT(al>=0);
  theOpName = n;
  key_length = kl;
  ans_length = al;
  is_marked_for_deletion = false;
  next = 0;
  discardStaleHits = true;

  // 
  // assign an index to this operation
  //
  if (free_list>=0) {
    oplist_index = free_list;
    free_list = op_holes[free_list];
  } else {
    if (list_size >= list_alloc) {
      int nla = list_alloc + 256;
      op_list = (operation**) realloc(op_list, nla * sizeof(void*));
      op_holes = (int*) realloc(op_holes, nla * sizeof(int));
      if (0==op_list || 0==op_holes) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      for (int i=list_size; i<list_alloc; i++) {
        op_list[i] = 0;
        op_holes[i] = -1;
      }
      list_alloc = nla;
    }
    oplist_index = list_size;
    list_size++;
  }
  op_list[oplist_index] = this;

  if (key_length) { // this op uses the CT
    // 
    // Initialize CT 
    //
    if (Monolithic_CT) {
      CT0 = Monolithic_CT;
    } else {
      CT0 = ct_initializer::createForOp(this);
    }
  } else {
    MEDDLY_DCASSERT(0==ans_length);
    CT0 = 0;
  }
}

#else   // OLD_OP_CT

MEDDLY::operation::operation(const opname* n, unsigned et_slots)
{
#ifdef DEBUG_CLEANUP
  fprintf(stdout, "Creating operation %p\n", this);
  fflush(stdout);
#endif
  theOpName = n;
  num_etids = et_slots;

  is_marked_for_deletion = false;
  next = 0;

  // 
  // assign an index to this operation
  //
  if (free_list>=0) {
    oplist_index = free_list;
    free_list = op_holes[free_list];
  } else {
    if (list_size >= list_alloc) {
      int nla = list_alloc + 256;
      op_list = (operation**) realloc(op_list, nla * sizeof(void*));
      op_holes = (int*) realloc(op_holes, nla * sizeof(int));
      if (0==op_list || 0==op_holes) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      for (int i=list_size; i<list_alloc; i++) {
        op_list[i] = 0;
        op_holes[i] = -1;
      }
      list_alloc = nla;
    }
    oplist_index = list_size;
    list_size++;
  }
  op_list[oplist_index] = this;

  //
  // Delay CT initialization!
  // The derived class hasn't set up the entry types yet!
  //
  CT = 0;

  //
  // Set up slots to save our entry_types.
  //
  if (et_slots) {
    etype = new compute_table::entry_type* [et_slots];
    for (unsigned i=0; i<et_slots; i++) {
      etype[i] = 0;
    }
  } else {
    etype = 0;
  }

  //
  // Allocate CTresults
  //
  if (et_slots) {
    CTresult = new compute_table::entry_result [et_slots];
  } else {
    CTresult = 0;
  }

  //
  // Allocate our slots
  //
  compute_table::registerOp(this, et_slots);
}

void MEDDLY::operation::buildCTs()
{
  if (0==num_etids) return;

  CT = new compute_table* [num_etids];

  if (Monolithic_CT) {
    for (unsigned i=0; i<num_etids; i++) {
      CT[i] = Monolithic_CT;
    }
  } else {
    for (unsigned i=0; i<num_etids; i++) {
      CT[i] = ct_initializer::createForOp(this, i);
    }
  }

  //
  // Initialize CTresults
  //
  for (unsigned i=0; i<num_etids; i++) {
    CTresult[i].initialize(etype[i]);
  }

  //
  // Most operations use only one slot
  //
  CT0 = CT[0];  
}

#endif  // OLD_OP_CT


MEDDLY::operation::~operation()
{
#ifdef DEBUG_CLEANUP
  fprintf(stdout, "Deleting operation %p %s\n", this, getName());
  fflush(stdout);
#endif

#ifdef OLD_OP_CT
  if (CT0 && (CT0!=Monolithic_CT)) delete CT0;
#else
  if (CT) {
    for (unsigned i=0; i<num_etids; i++) {
      if (CT[i] != Monolithic_CT)
        delete CT[i];
    }
    delete[] CT;
  }
  // Don't delete the entries in etype, they're owned by compute_table.
  delete[] etype;
  delete[] CTresult;
  compute_table::unregisterOp(this, num_etids);
#endif

  if (oplist_index >= 0) {
    MEDDLY_DCASSERT(op_list[oplist_index] == this);
    op_list[oplist_index] = 0;
    op_holes[oplist_index] = free_list;
    free_list = oplist_index;
  }
#ifdef DEBUG_CLEANUP
  fprintf(stdout, "Deleted operation %p %s\n", this, getName());
  fflush(stdout);
#endif
}

void MEDDLY::operation::removeStalesFromMonolithic()
{
  if (Monolithic_CT) Monolithic_CT->removeStales();
}

void MEDDLY::operation::removeAllFromMonolithic()
{
  if (Monolithic_CT) Monolithic_CT->removeAll();
}


void MEDDLY::operation::markForDeletion()
{
#ifdef DEBUG_CLEANUP
  fprintf(stdout, "Marking operation %p %s for deletion\n", this, getName());
  fflush(stdout);
#endif
  if (is_marked_for_deletion) return;
  is_marked_for_deletion = true;
#ifdef OLD_OP_CT
  if (CT0 && CT0->isOperationTable()) CT0->removeStales();
#else
  for (unsigned i=0; i<num_etids; i++) {
    etype[i]->markForDeletion();
  }
  if (CT) {
    for (unsigned i=0; i<num_etids; i++) {
      if (CT[i] && CT[i]->isOperationTable()) CT[i]->removeStales();
    }
  }
#endif
}

void MEDDLY::operation::destroyAllOps() 
{
  for (int i=0; i<list_size; i++) delete op_list[i];
  free(op_list);
  free(op_holes);
  op_list = 0;
  op_holes = 0;
  list_size = 0;
  list_alloc = 0;
  free_list = -1;
}

void MEDDLY::operation::removeStaleComputeTableEntries()
{
#ifdef OLD_OP_CT
  if (CT0) CT0->removeStales();
#else
  bool has_monolithic = false;
  if (CT) {
    for (unsigned i=0; i<num_etids; i++) {
      if (0==CT[i]) continue;
      if (CT[i]->isOperationTable()) {
        CT[i]->removeStales();
      } else {
        has_monolithic = true;
      }
    }
  }
  if (has_monolithic) {
    Monolithic_CT->removeStales();
  }
#endif
}

void MEDDLY::operation::removeAllComputeTableEntries()
{
#ifdef DEBUG_CLEANUP
  fprintf(stdout, "Removing entries for operation %p %s\n", this, getName());
  fflush(stdout);
#endif
  if (is_marked_for_deletion) return;
  is_marked_for_deletion = true;
#ifndef OLD_OP_CT
  for (unsigned i=0; i<num_etids; i++) {
    etype[i]->markForDeletion();
  }
#endif
  removeStaleComputeTableEntries();
#ifndef OLD_OP_CT
  for (unsigned i=0; i<num_etids; i++) {
    etype[i]->unmarkForDeletion();
  }
#endif
  is_marked_for_deletion = false;
#ifdef DEBUG_CLEANUP
  fprintf(stdout, "Removed entries for operation %p %s\n", this, getName());
  fflush(stdout);
#endif
}

void MEDDLY::operation::showMonolithicComputeTable(output &s, int verbLevel)
{
  if (Monolithic_CT) Monolithic_CT->show(s, verbLevel);
}

void MEDDLY::operation::showAllComputeTables(output &s, int verbLevel)
{
  if (Monolithic_CT) {
    Monolithic_CT->show(s, verbLevel);
    return;
  }
  for (int i=0; i<list_size; i++) 
    if (op_list[i]) {
      op_list[i]->showComputeTable(s, verbLevel);
    }
}

void MEDDLY::operation::showComputeTable(output &s, int verbLevel) const
{
#ifdef OLD_OP_CT
  if (CT0) CT0->show(s, verbLevel);
#else
  bool has_monolithic = false;
  if (CT) {
    for (unsigned i=0; i<num_etids; i++) {
      if (0==CT[i]) continue;
      if (CT[i]->isOperationTable()) {
        CT[i]->show(s, verbLevel);
      } else {
        has_monolithic = true;
      }
    }
  }
  if (has_monolithic) {
    Monolithic_CT->show(s, verbLevel);
  }
#endif
}

// ******************************************************************
// *                    unary_operation  methods                    *
// ******************************************************************

#ifdef OLD_OP_CT

MEDDLY::unary_operation::unary_operation(const unary_opname* code, int kl, 
  int al, expert_forest* arg, expert_forest* res) : operation(code, kl, al)
{
  argF = arg;
  resultType = FOREST;
  resF = res;

  registerInForest(argF);
  registerInForest(resF);

  setAnswerForest(resF);
}

MEDDLY::unary_operation::unary_operation(const unary_opname* code, int kl, 
  int al, expert_forest* arg, opnd_type res) : operation(code, kl, al)
{
  argF = arg;
  resultType = res;
  resF = 0;

  registerInForest(argF);

  setAnswerForest(0);
}

#else

MEDDLY::unary_operation::unary_operation(const unary_opname* code, 
  unsigned et_slots, expert_forest* arg, expert_forest* res)
: operation(code, et_slots)
{
  argF = arg;
  resultType = FOREST;
  resF = res;

  registerInForest(argF);
  registerInForest(resF);
}

MEDDLY::unary_operation::unary_operation(const unary_opname* code, 
  unsigned et_slots, expert_forest* arg, opnd_type res)
: operation(code, et_slots)
{
  argF = arg;
  resultType = res;
  resF = 0;

  registerInForest(argF);
}

#endif

MEDDLY::unary_operation::~unary_operation()
{
  unregisterInForest(argF);
  unregisterInForest(resF);
}

void MEDDLY::unary_operation::computeDDEdge(const dd_edge &arg, dd_edge &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, long &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, double &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, ct_object &c)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

// ******************************************************************
// *                    binary_operation methods                    *
// ******************************************************************

#ifdef OLD_OP_CT
MEDDLY::binary_operation::binary_operation(const binary_opname* op, int kl, 
  int al, expert_forest* arg1, expert_forest* arg2, expert_forest* res)
: operation(op, kl, al)
{
  arg1F = arg1;
  arg2F = arg2;
  resF = res;

  registerInForest(arg1F);
  registerInForest(arg2F);
  registerInForest(resF);

  setAnswerForest(resF);
  can_commute = false;
}

#else

MEDDLY::binary_operation::binary_operation(const binary_opname* op, 
  unsigned et_slots, expert_forest* arg1, expert_forest* arg2, expert_forest* res)
: operation(op, et_slots)
{
  arg1F = arg1;
  arg2F = arg2;
  resF = res;

  registerInForest(arg1F);
  registerInForest(arg2F);
  registerInForest(resF);

  can_commute = false;
}

#endif

MEDDLY::binary_operation::~binary_operation()
{
  unregisterInForest(arg1F);
  unregisterInForest(arg2F);
  unregisterInForest(resF);
}

MEDDLY::node_handle 
MEDDLY::binary_operation::compute(node_handle a, node_handle b)
{
  throw error(error::WRONG_NUMBER, __FILE__, __LINE__);
}

MEDDLY::node_handle 
MEDDLY::binary_operation::compute(int k, node_handle a, node_handle b)
{
  throw error(error::WRONG_NUMBER, __FILE__, __LINE__);
}

void MEDDLY::binary_operation::compute(int av, node_handle ap,
  int bv, node_handle bp, int &cv, node_handle &cp)
{
  throw error(error::WRONG_NUMBER, __FILE__, __LINE__);
}

void MEDDLY::binary_operation::compute(long av, node_handle ap,
  long bv, node_handle bp, long &cv, node_handle &cp)
{
  throw error(error::WRONG_NUMBER);
}

void MEDDLY::binary_operation::compute(long av, node_handle ap,
  node_handle bp, long &cv, node_handle &cp)
{
  throw error(error::WRONG_NUMBER);
}

void MEDDLY::binary_operation::compute(float av, node_handle ap,
  float bv, node_handle bp, float &cv, node_handle &cp)
{
  throw error(error::WRONG_NUMBER, __FILE__, __LINE__);
}

// ******************************************************************
// *                 specialized_operation  methods                 *
// ******************************************************************

#ifdef OLD_OP_CT
MEDDLY::
specialized_operation::
specialized_operation(const specialized_opname* op, int kl, int al) 
 : operation(op, kl, al)
{
}

#else
MEDDLY::
specialized_operation::
specialized_operation(const specialized_opname* op, unsigned et_slots)
 : operation(op, et_slots)
{
}

#endif

MEDDLY::specialized_operation::~specialized_operation()
{
}

void MEDDLY::specialized_operation::compute(const dd_edge &arg, dd_edge &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::specialized_operation::compute(const dd_edge &arg, bool &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::specialized_operation::compute(const dd_edge &ar1, 
  const dd_edge &ar2, dd_edge &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::specialized_operation::compute(double* y, const double* x)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::specialized_operation::compute(const dd_edge &ar1,
  const dd_edge &ar2, const dd_edge &ar3, dd_edge &res)
{
  throw error(error::TYPE_MISMATCH);
}

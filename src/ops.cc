
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
#include "old_meddly.h"
#include "old_meddly.hh"
#include "old_meddly_expert.h"
#include "old_meddly_expert.hh"
// #include "compute_table.h"

// #define DEBUG_CLEANUP

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
// *                 minimum_witness_opname methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_opname::constrained_opname(const char* n)
  : specialized_opname(n)
{
}

// ******************************************************************
// *                       operation  methods                       *
// ******************************************************************

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
  if (free_list) {
    oplist_index = free_list;
    free_list = op_holes[free_list];
  } else {
    oplist_index = ++list_size;
    if (list_size >= list_alloc) {
      unsigned nla = list_alloc + 256;
      op_list = (operation**) realloc(op_list, nla * sizeof(void*));
      op_holes = (unsigned*) realloc(op_holes, nla * sizeof(unsigned));
      if (0==op_list || 0==op_holes) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      list_alloc = nla;
      for (unsigned i=list_size; i<list_alloc; i++) {
        op_list[i] = 0;
        op_holes[i] = 0;
      }
    }
    if (0==list_size) {
      // Never use slot 0
      list_size++;
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


MEDDLY::operation::~operation()
{
#ifdef DEBUG_CLEANUP
  fprintf(stdout, "Deleting operation %p %s\n", this, getName());
  fflush(stdout);
#endif

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

  if (oplist_index) {
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
  for (unsigned i=0; i<num_etids; i++) {
    etype[i]->markForDeletion();
  }
  if (CT) {
    for (unsigned i=0; i<num_etids; i++) {
      if (CT[i] && CT[i]->isOperationTable()) CT[i]->removeStales();
    }
  }
}

void MEDDLY::operation::destroyAllOps()
{
  for (unsigned i=0; i<list_size; i++) delete op_list[i];
  free(op_list);
  free(op_holes);
  op_list = 0;
  op_holes = 0;
  list_size = 0;
  list_alloc = 0;
  free_list = 0;
}

void MEDDLY::operation::removeStaleComputeTableEntries()
{
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
}

void MEDDLY::operation::removeAllComputeTableEntries()
{
#ifdef DEBUG_CLEANUP
  fprintf(stdout, "Removing entries for operation %p %s\n", this, getName());
  fflush(stdout);
#endif
  if (is_marked_for_deletion) return;
  is_marked_for_deletion = true;
  for (unsigned i=0; i<num_etids; i++) {
    etype[i]->markForDeletion();
  }
  removeStaleComputeTableEntries();
  for (unsigned i=0; i<num_etids; i++) {
    etype[i]->unmarkForDeletion();
  }
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
  for (unsigned i=0; i<list_size; i++)
    if (op_list[i]) {
      op_list[i]->showComputeTable(s, verbLevel);
    }
}

void MEDDLY::operation::countAllNodeEntries(const expert_forest* f, size_t* counts)
{
  if (Monolithic_CT) {
    Monolithic_CT->countNodeEntries(f, counts);
  }
  for (unsigned i=0; i<list_size; i++)
    if (op_list[i]) {
      op_list[i]->countCTEntries(f, counts);
    }
}

void MEDDLY::operation::showComputeTable(output &s, int verbLevel) const
{
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
}

void MEDDLY::operation::countCTEntries(const expert_forest* f, size_t* counts) const
{
  if (CT) {
    for (unsigned i=0; i<num_etids; i++) {
      if (0==CT[i]) continue;
      if (CT[i]->isOperationTable()) {
        CT[i]->countNodeEntries(f, counts);
      }
    }
  }
}

// ******************************************************************
// *                    unary_operation  methods                    *
// ******************************************************************

MEDDLY::unary_operation::unary_operation(const unary_opname* code,
  unsigned et_slots, expert_forest* arg, expert_forest* res)
: operation(code, et_slots)
{
  argF = arg;
  resultType = opnd_type::FOREST;
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

MEDDLY::unary_operation::~unary_operation()
{
  unregisterInForest(argF);
  unregisterInForest(resF);
}

void MEDDLY::unary_operation::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
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

MEDDLY::binary_operation::~binary_operation()
{
  unregisterInForest(arg1F);
  unregisterInForest(arg2F);
  unregisterInForest(resF);
}

#ifdef KEEP_LL_COMPUTES

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

#endif

// ******************************************************************
// *                 specialized_operation  methods                 *
// ******************************************************************

MEDDLY::
specialized_operation::
specialized_operation(const specialized_opname* op, unsigned et_slots)
 : operation(op, et_slots)
{
}

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

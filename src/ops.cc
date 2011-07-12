
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

/*
    Implementation of operation framework.

    Actual operations are in separate files, in operations/ directory.
*/

#include "defines.h"
#include "compute_table.h"

namespace MEDDLY {
  extern settings meddlySettings;
}

// ******************************************************************
// *                       ct_object  methods                       *
// ******************************************************************

MEDDLY::ct_object::ct_object()
{
}

MEDDLY::ct_object::~ct_object()
{
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
  throw error(error::TYPE_MISMATCH);  
}

MEDDLY::unary_operation* 
MEDDLY::unary_opname::buildOperation(expert_forest* ar, opnd_type res) const
{
  throw error(error::TYPE_MISMATCH);
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
// *                    numerical_opname methods                    *
// ******************************************************************

MEDDLY::numerical_opname::numerical_opname(const char* n) : opname(n)
{
}

MEDDLY::numerical_opname::~numerical_opname()
{
}

// ******************************************************************
// *                       operation  methods                       *
// ******************************************************************

MEDDLY::operation::operation(const opname* n, int kl, int al)
{
  theOpName = n;
  key_length = kl;
  ans_length = al;
  is_marked_for_deletion = false;
  next = 0;
  discardStaleHits = true;

  if (0==key_length) {
    oplist_index = -1;
    CT = 0;
  } else {
    // 
    // assign an index to this operation, for use in CT
    //
    if (free_list>=0) {
      oplist_index = free_list;
      free_list = op_holes[free_list];
    } else {
      if (list_size >= list_alloc) {
        int nla = list_alloc + 256;
        op_list = (operation**) realloc(op_list, nla * sizeof(void*));
        op_holes = (int*) realloc(op_holes, nla * sizeof(int));
        if (0==op_list || 0==op_holes) throw error(error::INSUFFICIENT_MEMORY);
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
    // Initialize CT 
    //
    if (Monolithic_CT)
      CT = Monolithic_CT;
    else {
      CT = createOperationTable(meddlySettings.computeTable, this); 
    }

    //
    // Initialize CT search structure
    //
    CT->initializeSearchKey(CTsrch, this);

  }
}

MEDDLY::operation::~operation()
{
  if (CT && CT->isOperationTable()) delete CT;
  delete next;
  if (oplist_index >= 0) {
    DCASSERT(op_list[oplist_index] == this);
    op_list[oplist_index] = 0;
    op_holes[oplist_index] = free_list;
    free_list = oplist_index;
  }
}

void MEDDLY::operation::removeStalesFromMonolithic()
{
  if (Monolithic_CT) Monolithic_CT->removeStales();
}

void MEDDLY::operation::markForDeletion()
{
  is_marked_for_deletion = true;
  if (CT && CT->isOperationTable()) CT->removeStales();
}

void MEDDLY::operation::destroyOpList()
{
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
  if (CT) CT->removeStales();
}

void MEDDLY::operation::removeAllComputeTableEntries()
{
  if (is_marked_for_deletion) return;
  is_marked_for_deletion = true;
  if (CT) CT->removeStales();
  is_marked_for_deletion = false;
}

void MEDDLY::operation::showMonolithicComputeTable(FILE* s, int verbLevel)
{
  if (Monolithic_CT) Monolithic_CT->show(s, verbLevel);
}

void MEDDLY::operation::showAllComputeTables(FILE* s, int verbLevel)
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

void MEDDLY::operation::showComputeTable(FILE* s, int verbLevel) const
{
  if (CT) CT->show(s, verbLevel);
}

// ******************************************************************
// *                    unary_operation  methods                    *
// ******************************************************************

MEDDLY::unary_operation::unary_operation(const unary_opname* code, int kl, 
  int al, expert_forest* arg, expert_forest* res) : operation(code, kl, al)
{
  argF = arg;
  resultType = FOREST;
  resF = res;

  argFslot = argF->registerOperation(this);
  resFslot = resF->registerOperation(this);

  setAnswerForest(resF);
}

MEDDLY::unary_operation::unary_operation(const unary_opname* code, int kl, 
  int al, expert_forest* arg, opnd_type res) : operation(code, kl, al)
{
  argF = arg;
  resultType = res;
  resF = 0;

  argFslot = argF->registerOperation(this);
  resFslot = -1;

  setAnswerForest(0);
}

MEDDLY::unary_operation::~unary_operation()
{
  argF->unregisterOperation(this, argFslot);
  if (resF) resF->unregisterOperation(this, resFslot);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, dd_edge &res)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, long &res)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, double &res)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, ct_object &c)
{
  throw error(error::TYPE_MISMATCH);
}

// ******************************************************************
// *                    binary_operation methods                    *
// ******************************************************************

MEDDLY::binary_operation::binary_operation(const binary_opname* op, int kl, 
  int al, expert_forest* arg1, expert_forest* arg2, expert_forest* res)
: operation(op, kl, al)
{
  arg1F = arg1;
  arg2F = arg2;
  resF = res;

  arg1Fslot = arg1F->registerOperation(this);
  arg2Fslot = arg2F->registerOperation(this);
  resFslot = res->registerOperation(this);

  setAnswerForest(resF);
}

MEDDLY::binary_operation::~binary_operation()
{
  arg1F->unregisterOperation(this, arg1Fslot);
  arg2F->unregisterOperation(this, arg2Fslot);
  resF->unregisterOperation(this, resFslot);
}

int MEDDLY::binary_operation::compute(int a, int b)
{
  throw error(error::WRONG_NUMBER);
}

int MEDDLY::binary_operation::compute(int k, int a, int b)
{
  throw error(error::WRONG_NUMBER);
}

// ******************************************************************
// *                  numerical_operation  methods                  *
// ******************************************************************

MEDDLY::numerical_operation::numerical_operation(const numerical_opname* op) 
 : operation(op, 0, 0)
{
}

MEDDLY::numerical_operation::~numerical_operation()
{
}

void MEDDLY::numerical_operation::compute(double* y, const double* x)
{
  throw error(error::TYPE_MISMATCH);
}

// ******************************************************************
// *                     op_initializer methods                     *
// ******************************************************************

MEDDLY::op_initializer::op_initializer(op_initializer* bef)
{
  before = bef;
}

MEDDLY::op_initializer::~op_initializer()
{
  delete before;
}

void MEDDLY::op_initializer::initChain(const settings &s)
{
  if (before) before->initChain(s);
  init(s);
}

void MEDDLY::op_initializer::cleanupChain()
{
  cleanup();
  if (before) before->cleanupChain();
}


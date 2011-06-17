
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
    Implementation of classes for unary operations.
*/

#include "defines.h"


// ******************************************************************
// *                       operation  methods                       *
// ******************************************************************

MEDDLY::operation::operation()
{
  key_length = 0;
  ans_length = 0;
}

MEDDLY::operation::~operation()
{
}

// ******************************************************************
// *                    unary_operation  methods                    *
// ******************************************************************

MEDDLY::unary_operation::unary_operation(const unary_opcode* code, 
        expert_forest* arg, expert_forest* res)
{
  next = 0;
  opcode = code;
  argF = arg;
  resultType = FOREST;
  resF = res;
}

MEDDLY::unary_operation::unary_operation(const unary_opcode* code,
        expert_forest* arg, opnd_type res)
{
  next = 0;
  opcode = code;
  argF = arg;
  resultType = res;
  resF = 0;
}

MEDDLY::unary_operation::~unary_operation()
{
  delete next;
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
// *                     unary_builder  methods                     *
// ******************************************************************

MEDDLY::unary_builder::unary_builder(unary_opcode* oc)
{
  next = oc->builders;
  oc->builders = next;
  opcode = oc;
}

MEDDLY::unary_builder::~unary_builder()
{
  // library must be closing
  delete next;
}

bool 
MEDDLY::unary_builder::canBuild(const forest* arg, const forest* res) const
{
  return false;
}

bool MEDDLY::unary_builder::canBuild(const forest* arg, opnd_type res) const
{
  return false;
}

MEDDLY::unary_operation* 
MEDDLY::unary_builder::build(const forest* ar, const forest* rs) const
{
  throw error(error::UNKNOWN_OPERATION);  
}

MEDDLY::unary_operation* 
MEDDLY::unary_builder::build(const forest* ar, opnd_type res) const
{
  throw error(error::UNKNOWN_OPERATION);
}

// ******************************************************************
// *                      unary_opcode methods                      *
// ******************************************************************

int MEDDLY::unary_opcode::next_index;
MEDDLY::unary_opcode* MEDDLY::unary_opcode::list;

MEDDLY::unary_opcode::unary_opcode(const char* n)
{
  name = n;
  index = next_index;
  next_index++;
  builders = 0;
  next = list;
  list = this;
}

MEDDLY::unary_opcode::~unary_opcode()
{
  // library must be closing
  delete builders;
  delete next;
}

MEDDLY::unary_operation* 
MEDDLY::unary_opcode::buildOperation(const forest* arg, const forest* res) const
{
  unary_builder* curr;
  for (curr=builders; curr; curr=curr->getNext()) {
    if (curr->canBuild(arg, res)) {
      return curr->build(arg, res);
    }
  } // for
  throw error(error::TYPE_MISMATCH);
}

MEDDLY::unary_operation* 
MEDDLY::unary_opcode::buildOperation(const forest* arg, opnd_type res) const
{
  unary_builder* curr;
  for (curr=builders; curr; curr=curr->getNext()) {
    if (curr->canBuild(arg, res)) {
      return curr->build(arg, res);
    }
  } // for
  throw error(error::TYPE_MISMATCH);
}



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
#include "compute_table.h"

// ******************************************************************
// *                     unary_builder  methods                     *
// ******************************************************************

MEDDLY::unary_builder::unary_builder(unary_opcode* oc)
{
  next = oc->builders;
  oc->builders = this;
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


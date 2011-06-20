
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
    Implementation of the "global" functions given in
    meddly.h and meddly_expert.h.

*/

#include "defines.h"
#include "revision.h"
#include "compute_table.h"

// unary operations
#include "operations/cardinality.h"

namespace MEDDLY {
  // "global" variables

  settings meddlySettings;

  statistics meddlyStats;

  expert_compute_manager* ECM = 0;

  // unary operation "codes"

  const unary_opname* COPY = 0;
  const unary_opname* CARDINALITY = 0;
  const unary_opname* COMPLEMENT = 0;
  const unary_opname* MAX_RANGE = 0;
  const unary_opname* MIN_RANGE = 0;
  const unary_opname* CONVERT_TO_INDEX_SET = 0;

  // cache of operations
  operation** op_cache = 0;
  // size of cache
  int op_cache_size = 0;

  // operation settings
  bool& operation::useMonolithicCT(meddlySettings.useMonolithicComputeTable);

  // Monolithic compute table, if used
  compute_table* operation::Monolithic_CT = 0;

};

// helpers

void initStats(MEDDLY::statistics &s)
{
}

//----------------------------------------------------------------------
// front end - unary operations
//----------------------------------------------------------------------

MEDDLY::unary_operation* MEDDLY::getOperation(const unary_opname* code, 
  const expert_forest* arg, const expert_forest* res)
{
  if (code->getIndex()<0 || code->getIndex()>=op_cache_size)
    throw error(error::INVALID_OPERATION);

  operation* curr;
  operation* prev = 0;
  for (curr=op_cache[code->getIndex()]; curr; curr=curr->getNext()) {
    if (((unary_operation*)curr)->matches(arg, res)) {
      // move to front of list...
      if (prev) {
        prev->setNext(curr->getNext());
        curr->setNext(op_cache[code->getIndex()]);
        op_cache[code->getIndex()] = curr;
      }
      // ... and return
      return (unary_operation*) curr;
    }
    prev = curr;
  } // for

  // none present, build a new one...
  curr = code->buildOperation(arg, res);
  // ...move it to the front...
  curr->setNext(op_cache[code->getIndex()]);
  op_cache[code->getIndex()] = curr;
  // ...and return
  return (unary_operation*) curr;
}

MEDDLY::unary_operation* MEDDLY::getOperation(const unary_opname* code, 
  const expert_forest* arg, opnd_type res)
{
  if (code->getIndex()<0 || code->getIndex()>=op_cache_size)
    throw error(error::INVALID_OPERATION);

  operation* curr;
  operation* prev = 0;
  for (curr=op_cache[code->getIndex()]; curr; curr=curr->getNext()) {
    if (((unary_operation*)curr)->matches(arg, res)) {
      // move to front of list...
      if (prev) {
        prev->setNext(curr->getNext());
        curr->setNext(op_cache[code->getIndex()]);
        op_cache[code->getIndex()] = curr;
      }
      // ... and return
      return (unary_operation*) curr;
    }
    prev = curr;
  } // for

  // none present, build a new one...
  curr = code->buildOperation(arg, res);
  // ...move it to the front...
  curr->setNext(op_cache[code->getIndex()]);
  op_cache[code->getIndex()] = curr;
  // ...and return
  return (unary_operation*) curr;
}

void MEDDLY::removeOperationFromCache(operation* op)
{
  if (0==op) return;
  const opname* code = op->getOpName();

  operation* curr;
  operation* prev = 0;
  for (curr=op_cache[code->getIndex()]; curr; curr=curr->getNext()) {
    if (curr == op) break;
    prev = curr;
  } // for
  if (0==curr) return;  // not found
  // remove curr
  if (prev) {
    prev->setNext(curr->getNext());
  } else {
    op_cache[code->getIndex()] = curr->getNext();
  }
  curr->setNext(0);
}

void MEDDLY::apply(const unary_opname* code, const dd_edge &a, dd_edge &c)
{
  expert_forest* aF = (expert_forest*) a.getForest();
  expert_forest* cF = (expert_forest*) c.getForest();
  unary_operation* op = getOperation(code, aF, cF);
  op->compute(a, c);
}

void MEDDLY::apply(const unary_opname* code, const dd_edge &a, long &c)
{
  expert_forest* aF = (expert_forest*) a.getForest();
  unary_operation* op = getOperation(code, aF, INTEGER);
  op->compute(a, c);
}

void MEDDLY::apply(const unary_opname* code, const dd_edge &a, double &c)
{
  expert_forest* aF = (expert_forest*) a.getForest();
  unary_operation* op = getOperation(code, aF, REAL);
  op->compute(a, c);
}

void MEDDLY::apply(const unary_opname* code, const dd_edge &a, opnd_type cr,
  ct_object &c)
{
  expert_forest* aF = (expert_forest*) a.getForest();
  unary_operation* op = getOperation(code, aF, cr);
  op->compute(a, c);
}

//----------------------------------------------------------------------
// front end - create and destroy objects
//----------------------------------------------------------------------

MEDDLY::variable* MEDDLY::createVariable(int bound, char* name)
{
  if (0==ECM) throw error(error::UNINITIALIZED);
  return new expert_variable(bound, name);
}

MEDDLY::domain* MEDDLY::createDomain(variable** vars, int N)
{
  if (0==ECM) throw error(error::UNINITIALIZED);
  return new expert_domain(vars, N);
}

MEDDLY::domain* MEDDLY::createDomainBottomUp(const int* bounds, int N)
{
  if (0==ECM) throw error(error::UNINITIALIZED);
  domain* d = new expert_domain(0, 0);
  d->createVariablesBottomUp(bounds, N);
  return d;
}

void MEDDLY::destroyDomain(MEDDLY::domain* &d)
{
  if (0==d) return;
  if (0==ECM) throw error(error::UNINITIALIZED);
  expert_domain* ed = (expert_domain*) d;
  ed->markForDeletion();
  operation::removeStalesFromMonolithic();
  delete ed;
  d = 0;
}

void MEDDLY::destroyForest(MEDDLY::forest* &f)
{
  if (0==f) return;
  if (0==ECM) throw error(error::UNINITIALIZED);
  expert_forest* ef = (expert_forest*) f;
  ef->markForDeletion();
  operation::removeStalesFromMonolithic();
  delete ef;
  f = 0;
}

void MEDDLY::destroyOperation(MEDDLY::operation* &op)
{
  if (0==op) return;
  if (0==ECM) throw error(error::UNINITIALIZED);
  removeOperationFromCache(op);
  op->markForDeletion();
  operation::removeStalesFromMonolithic();
  delete op;
  op = 0;
}

MEDDLY::compute_manager* MEDDLY::getComputeManager()
{
  if (0==ECM) throw error(error::UNINITIALIZED);
  return ECM;
}

//----------------------------------------------------------------------
// front end - initialize and cleanup
//----------------------------------------------------------------------

void MEDDLY::initialize(settings s)
{
  if (ECM) throw error(error::ALREADY_INITIALIZED);
  meddlySettings = s;
  initStats(meddlyStats);

  // set up monolithic compute table, if needed
  if (meddlySettings.useMonolithicComputeTable) {
    compute_table::settings s;
    operation::Monolithic_CT = createMonolithicTable(s);
  }

  opname::next_index = 0;
  opname::list = 0;

  // Initialize opnames here

  CARDINALITY = initializeCardinality(s);
  // ...

  // set up operation cache
  op_cache_size = opname::next_index;
  op_cache = new operation*[op_cache_size];
  for (int i=0; i<op_cache_size; i++) {
    op_cache[i] = 0;
  }

  ECM = new expert_compute_manager(meddlySettings);
}

void MEDDLY::cleanup()
{
  if (0==ECM) throw error(error::UNINITIALIZED);
  delete ECM;
  ECM = 0;

  // TBD: cleanup unary stuff
}

//----------------------------------------------------------------------
// front end - library info
//----------------------------------------------------------------------

const MEDDLY::settings& MEDDLY::getLibrarySettings()
{
  return meddlySettings;
}

const MEDDLY::statistics& MEDDLY::getLibraryStats()
{
  return meddlyStats;
}

const char* MEDDLY::getLibraryInfo(int what)
{
  static char* title = 0;
  switch (what) {
    case 0:
      if (!title) {
        title = new char[80];
        if (REVISION_NUMBER) {
          snprintf(title, 80, 
            "%s version %s.%d (32-bit and 64-bit compatible)", 
            PACKAGE_NAME, VERSION, REVISION_NUMBER
          );
        } else {
          snprintf(title, 80, 
            "%s version %s (32-bit and 64-bit compatible)", 
            PACKAGE_NAME, VERSION
          );
        }
      }
      return title;

    case 1:
      return "Copyright (C) 2009, Iowa State University Research Foundation, Inc.";

    case 2:
      return "Released under the GNU Lesser General Public License, version 3";
 
    case 3:
      return "http://meddly.sourceforge.net/";

    case 4:
      return "Data Structures and operations available:\n\
(1) MDDs: Union, Intersection, Difference.\n\
(2) Matrix Diagrams (MXDs): Union, Intersection, Difference.\n\
(3) Multi-Terminal MDDs (MTMDDs) with integer or real terminals:\n\
    Arithmetic: Plus, Minus, Multiply, Divide, Min, Max.\n\
    Logical: <, <=, >, >=, ==, !=.\n\
    Conversion to and from MDDs.\n\
(4) Multi-Terminal MXDs (MTMXDs) with integer or real terminals:\n\
    Arithmetic: Plus, Minus, Multiply, Divide, Min, Max.\n\
    Logical: <, <=, >, >=, ==, !=.\n\
    Conversion to and from MXDs.\n\
";
  }
  return 0;
}



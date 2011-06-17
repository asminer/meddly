
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


namespace MEDDLY {
  // "global" variables

  settings meddlySettings;

  statistics meddlyStats;

  expert_compute_manager* ECM = 0;

  // unary operation "codes"

  const unary_opcode* COPY = 0;
  const unary_opcode* CARDINALITY = 0;
  const unary_opcode* COMPLEMENT = 0;
  const unary_opcode* MAX_RANGE = 0;
  const unary_opcode* MIN_RANGE = 0;
  const unary_opcode* CONVERT_TO_INDEX_SET = 0;

  // cache of unary operations
  unary_operation** unary_cache = 0;
  // size of unary cache
  int unary_cache_size = 0;


  // Monolithic compute table, if used
  compute_table* Monolithic_CT = 0;

};

// helpers

void initStats(MEDDLY::statistics &s)
{
}

//----------------------------------------------------------------------
// front end - unary operations
//----------------------------------------------------------------------

MEDDLY::unary_operation* MEDDLY::getOperation(const unary_opcode* code, 
  const forest* arg, const forest* res)
{
  if (code->getIndex()<0 || code->getIndex()>=unary_cache_size)
    throw error(error::INVALID_OPERATION);

  unary_operation* curr;
  unary_operation* prev = 0;
  for (curr=unary_cache[code->getIndex()]; curr; curr=curr->getNext()) {
    if (curr->matches((expert_forest*)arg, (expert_forest*)res)) {
      // move to front of list...
      if (prev) {
        prev->setNext(curr->getNext());
        curr->setNext(unary_cache[code->getIndex()]);
        unary_cache[code->getIndex()] = curr;
      }
      // ... and return
      return curr;
    }
    prev = curr;
  } // for

  // none present, build a new one...
  curr = code->buildOperation(arg, res);
  // ...move it to the front...
  curr->setNext(unary_cache[code->getIndex()]);
  unary_cache[code->getIndex()] = curr;
  // ...and return
  return curr;
}

MEDDLY::unary_operation* MEDDLY::getOperation(const unary_opcode* code, 
  const forest* arg, opnd_type res)
{
  if (code->getIndex()<0 || code->getIndex()>=unary_cache_size)
    throw error(error::INVALID_OPERATION);

  unary_operation* curr;
  unary_operation* prev = 0;
  for (curr=unary_cache[code->getIndex()]; curr; curr=curr->getNext()) {
    if (curr->matches((expert_forest*)arg, res)) {
      // move to front of list...
      if (prev) {
        prev->setNext(curr->getNext());
        curr->setNext(unary_cache[code->getIndex()]);
        unary_cache[code->getIndex()] = curr;
      }
      // ... and return
      return curr;
    }
    prev = curr;
  } // for

  // none present, build a new one...
  curr = code->buildOperation(arg, res);
  // ...move it to the front...
  curr->setNext(unary_cache[code->getIndex()]);
  unary_cache[code->getIndex()] = curr;
  // ...and return
  return curr;
}

void MEDDLY::apply(const unary_opcode* code, const dd_edge &a, dd_edge &c)
{
  unary_operation* op = getOperation(code, a.getForest(), c.getForest());
  op->compute(a, c);
}

void MEDDLY::apply(const unary_opcode* code, const dd_edge &a, long &c)
{
  unary_operation* op = getOperation(code, a.getForest(), INTEGER);
  op->compute(a, c);
}

void MEDDLY::apply(const unary_opcode* code, const dd_edge &a, double &c)
{
  unary_operation* op = getOperation(code, a.getForest(), REAL);
  op->compute(a, c);
}

void MEDDLY::apply(const unary_opcode* code, const dd_edge &a, opnd_type cr,
  ct_object &c)
{
  unary_operation* op = getOperation(code, a.getForest(), cr);
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
  // delete registered forests
  for (int i = 0; i < ed->nForests; ++i) {
    delete ed->forests[i];
  }
  delete d;
  d = 0;
}

void MEDDLY::destroyForest(MEDDLY::forest* &f)
{
  if (0==f) return;
  if (0==ECM) throw error(error::UNINITIALIZED);
  expert_domain* ed = (expert_domain*)f->useDomain();
  ed->unlinkForest(f);
  delete f;
  f = 0;
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

  unary_opcode::next_index = 0;
  unary_opcode::list = 0;

  // initialize unary ops here

  // set up operation cache
  unary_cache_size = unary_opcode::next_index;
  unary_cache = new unary_operation*[unary_cache_size];
  for (int i=0; i<unary_cache_size; i++) {
    unary_cache[i] = 0;
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



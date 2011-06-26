
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



#include "defines.h"
#include "operations/operation_ext.h"
// #include "operations/reachability.h"
#include "operations/cardinality.h"
#include "operations/maxmin_range.h"
#include "operations/vect_matr.h"
#include "operations/cross.h"
#include "compute_cache.h"

// ------------------ compute_manager --------------------------
MEDDLY::compute_manager::compute_manager() {}


MEDDLY::compute_manager::~compute_manager() {}




// ------------------------ expert_compute_manager --------------------

MEDDLY::expert_compute_manager::expert_compute_manager(const settings &s)
: useCTchaining(s.doComputeTablesUseChaining),
  maxCTsize(s.maxComputeTableSize)
{
  // initialize compute cache
  cc = new compute_cache();
  cc->setPolicy(useCTchaining, maxCTsize);
  // initialize builtinOpEntries and customOpEntries
  builtinOpEntries = new std::map<builtin_op_key, op_info>();
  customOpEntries = new std::map<custom_op_key, op_info>();
}


MEDDLY::expert_compute_manager::~expert_compute_manager()
{
  delete customOpEntries;
  delete builtinOpEntries;
  delete cc;
}


const char* MEDDLY::expert_compute_manager::getOperationName(
    compute_manager::op_code op) const
{
  return "Unknown operation";
}

template <class TYPE>
inline void
unary_apply(op_info* owner, const dd_edge &a, TYPE &b)
{
  // type check
  owner->op->typeCheck(owner);
  owner->op->compute(owner, a, b);
}

void MEDDLY::expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, dd_edge &b)
{
  return unary_apply(owner, a, b);
}

void MEDDLY::expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, long &b)
{
  return unary_apply(owner, a, b);
}

void MEDDLY::expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, double &b)
{
  return unary_apply(owner, a, b);
}

void MEDDLY::expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, ct_object &b)
{
#ifdef HAVE_LIBGMP
  unary_apply(owner, a, b);
#else
  throw error(error::UNKNOWN_OPERATION);
#endif
}


void MEDDLY::expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  // type check
  owner->op->typeCheck(owner);
  owner->op->compute(owner, a, b, c);
}


template <class TYPE>
inline void 
unary_apply(expert_compute_manager* CM, compute_manager::op_code op, 
            const dd_edge &a, TYPE &b)
{
  static op_param plist[2];
  plist[0].set(a);
  plist[1].set(b);
  op_info* opInfo = CM->getOpInfo(op, plist, 2);
  if (0==opInfo) throw error(error::UNKNOWN_OPERATION);
  CM->apply(opInfo, a, b);
}

void MEDDLY::expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, dd_edge &b)
{
  unary_apply(this, op, a, b);
}

void MEDDLY::expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, long &b)
{
  unary_apply(this, op, a, b);
}

void MEDDLY::expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, double &b)
{
  unary_apply(this, op, a, b);
}

void MEDDLY::expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, ct_object &b)
{
#ifdef HAVE_LIBGMP
  unary_apply(this, op, a, b);
#else
  throw error(error::UNKNOWN_OPERATION);
#endif
}


void MEDDLY::expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, const dd_edge &b,
    dd_edge &c)
{
  static const int nForests = 3;
  static op_param plist[nForests];
  plist[0].set(a);
  plist[1].set(b);
  plist[2].set(c);
  op_info* opInfo = getOpInfo(op, plist, nForests);
  if (0==opInfo) throw error(error::UNKNOWN_OPERATION);
  apply(opInfo, a, b, c);
}


void MEDDLY::expert_compute_manager::showComputeTable(FILE* strm) const
{
  if (cc != 0) if (cc->getNumEntries() > 0) { cc->show(strm, true); }

  // Other compute tables may exist (!= cc). Clear those as well.
  std::map<builtin_op_key, op_info>::iterator curr = builtinOpEntries->begin();
  std::map<builtin_op_key, op_info>::iterator end = builtinOpEntries->end();
  for ( ; curr != end; ++curr) {
    compute_cache* cache = (curr->second).cc;
    if (cache != cc) cache->show(strm, true);
  }
}


long MEDDLY::expert_compute_manager::getNumCacheEntries() const
{
  long sum = (cc == 0)? 0: cc->getNumEntries();

  // Other compute tables may exist (!= cc). Clear those as well.
  std::map<builtin_op_key, op_info>::const_iterator curr, end;
  curr = builtinOpEntries->begin();
  end = builtinOpEntries->end();
  for ( ; curr != end; ++curr) {
    compute_cache* cache = (curr->second).cc;
    if (cache != cc) sum += cache->getNumEntries();
  }

  return sum;
}


void MEDDLY::expert_compute_manager::removeStales(expert_forest* f)
{
  // remove stales from the main computa table (cc).
  if (cc) cc->removeStales();

  // Other compute tables may exist (!= cc). Remove stales from those also.
  std::map<builtin_op_key, op_info>::iterator curr, end;
  curr = builtinOpEntries->begin();
  end = builtinOpEntries->end();
  for ( ; curr != end; ++curr) {
    op_info& opInfo = curr->second;
    if (opInfo.cc == cc) continue;
    bool clearCC = false;
    for (int i = 0; i < opInfo.nParams; ++i) {
      if (! opInfo.p[i].isForest()) continue;
      if (opInfo.p[i].getForest() == f) { clearCC = true; break; }
    }
    if (clearCC) {
      // fprintf(stderr, "Removing stale entries for %s\n", opInfo.op->getName());
      opInfo.cc->removeStales();
    }
  }
}


void MEDDLY::expert_compute_manager::removeStales(op_info* owner)
{
  if (owner == 0) {
    if (cc) cc->removeStales();
  } else if (owner->cc == 0 || owner->cc == cc) {
    // Entries for this owner are stored in the main compute table.
    if (cc) cc->removeStales(owner);
  } else {
    DCASSERT(owner->cc);
    // Entries for this owner are stored in custom compute table.
    fprintf(stderr, "Clearing cache table for %s\n", owner->op->getName());
    owner->cc->removeStales();
  }
}


void expert_compute_manager::removeEntries(op_info* owner)
{
  if (owner == 0) {
    if (cc) cc->removeEntries();
  } else if (owner->cc == 0 || owner->cc == cc) {
    // Entries for this owner are stored in the main compute table.
    if (cc) cc->removeEntries(owner);
  } else {
    DCASSERT(owner->cc);
    // Entries for this owner are stored in custom compute table.
    fprintf(stderr, "Clearing cache table for %s\n", owner->op->getName());
    owner->cc->removeEntries();
  }
}


void expert_compute_manager::clearComputeTable()
{
  if (cc) cc->clear();

  // Other compute tables may exist (!= cc). Clear those as well.
  std::map<builtin_op_key, op_info>::iterator curr = builtinOpEntries->begin();
  std::map<builtin_op_key, op_info>::iterator end = builtinOpEntries->end();
  for ( ; curr != end; ++curr) {
    compute_cache* cache = (curr->second).cc;
    if (cache != cc) cache->clear();
  }
}


void MEDDLY::expert_compute_manager::addBuiltinOp(const builtin_op_key& key,
  const old_operation* op, const op_param* plist, int n)
{
  compute_cache* cache = cc;

#ifdef USE_BINARY_COMPUTE_CACHE
  if (op->getCacheEntryLength() == 3 && n == 3) {
    bool allForests = true;
    for (int i = 0; i < n; ++i) {
      if (! plist[i].isForest()) { allForests = false; break; }
    }
    if (allForests) {
      // fprintf(stderr, "Using binary compute cache for %s\n", op->getName());
      cache = new binary_compute_cache(op, plist, n);
      cache->setPolicy(useCTchaining, maxCTsize);
    }
  }
#endif

  op_info entry(const_cast<old_operation*>(op), const_cast<op_param*>(plist),
      n, cache);

  (*builtinOpEntries)[key] = entry;
#ifdef DEVELOPMENT_CODE
  std::map<builtin_op_key, op_info>::iterator itr = builtinOpEntries->find(key);
  if (itr == builtinOpEntries->end()) {
    // print the entries
    itr = builtinOpEntries->begin();
    while (itr != builtinOpEntries->end()) {
      printf("{");
      itr->first.print(stdout);
      itr->second.print(stdout);
      printf("}\n");
      ++itr;
    }
    assert(false);
  }
  assert(itr->second == entry);
#endif
}


op_info* MEDDLY::expert_compute_manager::getOpInfo(compute_manager::op_code op,
    const op_param* const plist, int N)
{
  return 0;
}


op_info* MEDDLY::expert_compute_manager::getOpInfo(const old_operation* op,
    const op_param* const plist, int N)
{
  // search in custom op entries
  custom_op_key key(op, plist, N);
  std::map<custom_op_key, op_info>::iterator itr = customOpEntries->find(key);
  if (itr == customOpEntries->end()) {
    // add new entry
    op_info entry(const_cast<old_operation*>(op), 
        const_cast<op_param*>(plist), N, cc);
    (*customOpEntries)[key] = entry;
    itr = customOpEntries->find(key);
#ifdef DEVELOPMENT_CODE
    assert (itr != customOpEntries->end());
    assert(itr->second == entry);
#endif
  }
  return &(itr->second);
}

void 
MEDDLY::expert_compute_manager::vectorMatrixMultiply(double* y, const dd_edge &y_ind,
                      const double* x, const dd_edge &x_ind, const dd_edge &A)
{
  const expert_forest* const fy = (expert_forest*) y_ind.getForest();
  const expert_forest* const fx = (expert_forest*) x_ind.getForest();
  const expert_forest* const fA = (expert_forest*) A.getForest();

  if (
           (fy->getRangeType() != forest::INTEGER) 
        || (fy->isForRelations())
        || (fx->getRangeType() != forest::INTEGER)
        || (fx->isForRelations())
        || (fA->getRangeType() != forest::REAL)
        || (!fA->isForRelations())
      ) 
  {
    throw error(error::TYPE_MISMATCH);
  }

  // A can't be fully reduced.
  if (forest::FULLY_REDUCED == fA->getReductionRule()) {
    throw error(error::TYPE_MISMATCH);
  }

  // For now, fy and fx must be EV+MDDs.
  if (     (fy->getEdgeLabeling() != forest::EVPLUS) 
        || (fx->getEdgeLabeling() != forest::EVPLUS) )
  {
    throw error(error::NOT_IMPLEMENTED);
  }

  //everyone must use the same domain
  if (      (fx->getDomain() != fy->getDomain()) 
        ||  (fx->getDomain() != fA->getDomain())  )
  {
    throw error(error::TYPE_MISMATCH);
  }

  static op_param plist[5];
  plist[0].set(y);
  plist[1].set(y_ind);
  plist[2].set(x);
  plist[3].set(x_ind);
  plist[4].set(A);

  switch (fA->getEdgeLabeling()) {
    case forest::MULTI_TERMINAL:
      return vectorMatrixMult_evplus_mt(
        plist, fy->getDomain()->getNumVariables(), y, y_ind.getNode(),
        x, x_ind.getNode(), A.getNode()
      );

    case forest::EVTIMES:
      return vectorMatrixMult_evplus_evtimes(
        plist, fy->getDomain()->getNumVariables(), y, y_ind.getNode(),
        x, x_ind.getNode(), A.getNode()
      );

    default:
      throw error(error::TYPE_MISMATCH);
  };


}


void 
MEDDLY::expert_compute_manager::matrixVectorMultiply(double* y, const dd_edge &y_ind,
                      const dd_edge &A, const double* x, const dd_edge &x_ind)
{
  const expert_forest* const fy = (expert_forest*) y_ind.getForest();
  const expert_forest* const fA = (expert_forest*) A.getForest();
  const expert_forest* const fx = (expert_forest*) x_ind.getForest();

  if (
           (fy->getRangeType() != forest::INTEGER) 
        || (fy->isForRelations())
        || (fx->getRangeType() != forest::INTEGER)
        || (fx->isForRelations())
        || (fA->getRangeType() != forest::REAL)
        || (!fA->isForRelations())
      ) 
  {
    throw error(error::TYPE_MISMATCH);
  }

  // A can't be fully reduced.
  if (forest::FULLY_REDUCED == fA->getReductionRule()) {
    throw error(error::TYPE_MISMATCH);
  }

  // For now, fy and fx must be EV+MDDs.
  if (     (fy->getEdgeLabeling() != forest::EVPLUS) 
        || (fx->getEdgeLabeling() != forest::EVPLUS) )
  {
    throw error(error::NOT_IMPLEMENTED);
  }

  //everyone must use the same domain
  if (      (fx->getDomain() != fy->getDomain()) 
        ||  (fx->getDomain() != fA->getDomain())  )
  {
    throw error(error::TYPE_MISMATCH);
  }

  static op_param plist[5];
  plist[0].set(y);
  plist[1].set(y_ind);
  plist[2].set(A);
  plist[3].set(x);
  plist[4].set(x_ind);

  switch (fA->getEdgeLabeling()) {
    case forest::MULTI_TERMINAL:
      return matrixVectorMult_evplus_mt(
        plist, fy->getDomain()->getNumVariables(), y, y_ind.getNode(),
        A.getNode(), x, x_ind.getNode()
      );

    case forest::EVTIMES:
      return matrixVectorMult_evplus_evtimes(
        plist, fy->getDomain()->getNumVariables(), y, y_ind.getNode(),
        A.getNode(), x, x_ind.getNode()
      );

    default:
      throw error(error::TYPE_MISMATCH);
  };


}




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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_LIBGMP
#include <gmp.h>
#endif
#include "../defines.h"
#include "cardinality.h"
#include "mpz_object.h"

// #define DEBUG_CARD

namespace MEDDLY {
  class card_int;
  class card_mdd_int;
  class card_mxd_int;

  class card_real;
  class card_mdd_real;
  class card_mxd_real;

#ifdef HAVE_LIBGMP
  class card_mpz;
  class card_mdd_mpz;
  class card_mxd_mpz;
#endif

  class card_opname;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        Card  operations                        *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                         card_int class                         *
// *                                                                *
// ******************************************************************

//  Abstract base class: cardinality that returns an integer
class MEDDLY::card_int : public unary_operation {
public:
  card_int(const unary_opname* oc, expert_forest* arg);

  // common
  virtual bool isStaleEntry(const node_handle* entryData);
  virtual void discardEntry(const node_handle* entryData);
  virtual void showEntry(FILE* strm, const node_handle* entryData) const;

protected:
  static inline void overflow_acc(long &a, long x) {
    a += x;
    if (a < x) throw error(error::OVERFLOW);
  }
  static inline long overflow_mult(long a, long x) {
    a *= x;
    if (a < x) throw error(error::OVERFLOW);
    return a;
  }
};

MEDDLY::card_int::card_int(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, 1, sizeof(long) / sizeof(int), arg, INTEGER)
{
}

bool MEDDLY::card_int::isStaleEntry(const node_handle* data)
{
  return argF->isStale(data[0]);
}

void MEDDLY::card_int::discardEntry(const node_handle* data)
{
  argF->uncacheNode(data[0]);
}

void MEDDLY::card_int::showEntry(FILE* strm, const node_handle* data) const
{
  long answer;
  memcpy(&answer, data+1, sizeof(long));
  fprintf(strm, "[%s(%d): %ld(L)]", getName(), data[0], answer);
}

// ******************************************************************
// *                                                                *
// *                       card_mdd_int class                       *
// *                                                                *
// ******************************************************************

//  Cardinality on MDDs, returning integer
class MEDDLY::card_mdd_int : public card_int {
public:
  card_mdd_int(const unary_opname* oc, expert_forest* arg)
    : card_int(oc, arg) { }
  virtual void compute(const dd_edge &arg, long &res) {
    res = compute(argF->getDomain()->getNumVariables(), arg.getNode());
  }
  long compute(int k, node_handle a);
};

long MEDDLY::card_mdd_int::compute(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0;
  if (0==k) return 1;

  // Quickly deal with skipped levels
  if (argF->getNodeLevel(a) < k) {
    return overflow_mult(compute(k-1, a), argF->getLevelSize(k));
  }
  
  // Check compute table
  CTsrch.key(0) = a; 
  const node_handle* cacheEntry = CT->find(CTsrch);
  if (cacheEntry) {
    // ugly but portable
    long answer;
    memcpy(&answer, cacheEntry+1, sizeof(long));
    return answer;
  }

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, false);

  // Recurse
  long card = 0;
  int kdn = k-1;
  for (int z=0; z<A->getNNZs(); z++) {
    overflow_acc(card, compute(kdn, A->d(z)));
  }
  
  // Cleanup
  node_reader::recycle(A);

  // Add entry to compute table
  compute_table::temp_entry &entry = CT->startNewEntry(this);
  entry.key(0) = argF->cacheNode(a);
  entry.copyResult(0, &card, sizeof(long));
  CT->addEntry();

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %ld(L)\n", a, card);
#endif
  return card;
}


// ******************************************************************
// *                                                                *
// *                       card_mxd_int class                       *
// *                                                                *
// ******************************************************************

//  Cardinality on MxDs, returning integer
class MEDDLY::card_mxd_int : public card_int {
public:
  card_mxd_int(const unary_opname* oc, expert_forest* arg)
    : card_int(oc, arg) { }
  virtual void compute(const dd_edge &arg, long &res) {
    res = compute(argF->getDomain()->getNumVariables(), arg.getNode());
  }
  long compute(int k, node_handle a);
};

long MEDDLY::card_mxd_int::compute(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0;
  if (0==k) return 1;
  
  // Quickly deal with skipped levels
  if (isLevelAbove(k, argF->getNodeLevel(a))) {
    if (k<0 && argF->isIdentityReduced()) {
      // identity node
      return compute(-k-1, a);
    }
    // redundant node
    return overflow_mult(compute(-k, a), argF->getLevelSize(k));
  }
  
  // Check compute table
  CTsrch.key(0) = a; 
  const node_handle* cacheEntry = CT->find(CTsrch);
  if (cacheEntry) {
    // ugly but portable
    long answer;
    memcpy(&answer, cacheEntry+1, sizeof(long));
    return answer;
  }

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, false);

  // Recurse
  long card = 0;
  int kdn = (k<0) ? -k-1 : -k;
  for (int z=0; z<A->getNNZs(); z++) {
    overflow_acc(card, compute(kdn, A->d(z)));
  }
  
  // Cleanup
  node_reader::recycle(A);

  // Add entry to compute table
  compute_table::temp_entry &entry = CT->startNewEntry(this);
  entry.key(0) = argF->cacheNode(a);
  entry.copyResult(0, &card, sizeof(long));
  CT->addEntry();

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %ld(L)\n", a, card);
#endif
  return card;
}


// ******************************************************************
// *                                                                *
// *                        card_real  class                        *
// *                                                                *
// ******************************************************************

//  Abstract base class: cardinality that returns a real
class MEDDLY::card_real : public unary_operation {
public:
  card_real(const unary_opname* oc, expert_forest* arg);

  // common
  virtual bool isStaleEntry(const node_handle* entryData);
  virtual void discardEntry(const node_handle* entryData);
  virtual void showEntry(FILE* strm, const node_handle* entryData) const;
};

MEDDLY::card_real::card_real(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, 1, sizeof(double) / sizeof(int), arg, REAL)
{
}

bool MEDDLY::card_real::isStaleEntry(const node_handle* data)
{
  return argF->isStale(data[0]);
}

void MEDDLY::card_real::discardEntry(const node_handle* data)
{
  argF->uncacheNode(data[0]);
}

void MEDDLY::card_real::showEntry(FILE* strm, const node_handle* data) const
{
  double answer;
  memcpy(&answer, data+1, sizeof(double));
  fprintf(strm, "[%s(%d): %le]", getName(), data[0], answer);
}


// ******************************************************************
// *                                                                *
// *                      card_mdd_real  class                      *
// *                                                                *
// ******************************************************************

//  Cardinality on MDDs, returning real
class MEDDLY::card_mdd_real : public card_real {
public:
  card_mdd_real(const unary_opname* oc, expert_forest* arg)
    : card_real(oc, arg) { }
  virtual void compute(const dd_edge &arg, double &res) {
    res = compute(argF->getDomain()->getNumVariables(), arg.getNode());
  }
  double compute(int ht, node_handle a);
};

double MEDDLY::card_mdd_real::compute(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0.0;
  if (0==k) return 1.0;

  // Quickly deal with skipped levels
  if (argF->getNodeLevel(a) < k) {
    return compute(k-1, a) * argF->getLevelSize(k);
  }
  
  // Check compute table
  CTsrch.key(0) = a; 
  const node_handle* cacheEntry = CT->find(CTsrch);
  if (cacheEntry) {
    // ugly but portable
    double answer;
    memcpy(&answer, cacheEntry+1, sizeof(double));
    return answer;
  }

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, false);

  // Recurse
  double card = 0;
  int kdn = k-1;
  for (int z=0; z<A->getNNZs(); z++) {
    card += compute(kdn, A->d(z));
  }
  
  // Cleanup
  node_reader::recycle(A);

  // Add entry to compute table
  compute_table::temp_entry &entry = CT->startNewEntry(this);
  entry.key(0) = argF->cacheNode(a);
  entry.copyResult(0, &card, sizeof(double));
  CT->addEntry();

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %le(L)\n", a, card);
#endif
  return card;
}



// ******************************************************************
// *                                                                *
// *                      card_mxd_real  class                      *
// *                                                                *
// ******************************************************************

//  Cardinality on MxDs, returning real
class MEDDLY::card_mxd_real : public card_real {
public:
  card_mxd_real(const unary_opname* oc, expert_forest* arg)
    : card_real(oc, arg) { }
  virtual void compute(const dd_edge &arg, double &res) {
    res = compute(argF->getDomain()->getNumVariables(), arg.getNode());
  }
  double compute(int k, node_handle a);
};

double MEDDLY::card_mxd_real::compute(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0.0;
  if (0==k) return 1.0;

  // Quickly deal with skipped levels
  if (isLevelAbove(k, argF->getNodeLevel(a))) {
    if (k<0 && argF->isIdentityReduced()) {
      // identity node
      return compute(-k-1, a);
    }
    // redundant node
    return compute(-k, a) * argF->getLevelSize(k);
  }
  
  // Check compute table
  CTsrch.key(0) = a; 
  const node_handle* cacheEntry = CT->find(CTsrch);
  if (cacheEntry) {
    // ugly but portable
    double answer;
    memcpy(&answer, cacheEntry+1, sizeof(double));
    return answer;
  }

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, false);

  // Recurse
  double card = 0;
  int kdn = (k<0) ? -k-1 : -k;
  for (int z=0; z<A->getNNZs(); z++) {
    card += compute(kdn, A->d(z));
  }
  
  // Cleanup
  node_reader::recycle(A);

  // Add entry to compute table
  compute_table::temp_entry &entry = CT->startNewEntry(this);
  entry.key(0) = argF->cacheNode(a);
  entry.copyResult(0, &card, sizeof(double));
  CT->addEntry();

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %le\n", a, card);
#endif
  return card;
}




// ******************************************************************
// *                                                                *
// *                         card_mpz class                         *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP

//  Abstract base class: cardinality that returns large (mpz) integers.
class MEDDLY::card_mpz : public unary_operation {
public:
  card_mpz(const unary_opname* oc, expert_forest* arg);

  // common
  virtual bool isStaleEntry(const node_handle* entryData);
  virtual void discardEntry(const node_handle* entryData);
  virtual void showEntry(FILE* strm, const node_handle* entryData) const;
};

MEDDLY::card_mpz::card_mpz(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, 1, sizeof(ct_object*) / sizeof(int), arg, HUGEINT)
{
}

bool MEDDLY::card_mpz::isStaleEntry(const node_handle* data)
{
  return argF->isStale(data[0]);
}

void MEDDLY::card_mpz::discardEntry(const node_handle* data)
{
  argF->uncacheNode(data[0]);
  mpz_object* answer;
  memcpy(&answer, data+1, sizeof(mpz_object*));
  delete answer;
}

void MEDDLY::card_mpz::showEntry(FILE* strm, const node_handle* entryData) const
{
  mpz_object* answer;
  memcpy(&answer, entryData+1, sizeof(mpz_object*));
  fprintf(strm, "[%s(%d): ", getName(), entryData[0]);
  answer->show(strm);
  fprintf(strm, "]");
}

#endif

// ******************************************************************
// *                                                                *
// *                       card_mdd_mpz class                       *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP

/// Cardinality of MDDs, returning large (mpz) integers.
class MEDDLY::card_mdd_mpz : public card_mpz {
public:
  card_mdd_mpz(const unary_opname* oc, expert_forest* arg)
    : card_mpz(oc, arg) { }
  virtual void compute(const dd_edge& a, ct_object &res) {
    mpz_object& mcard = dynamic_cast <mpz_object &> (res);
    compute(argF->getDomain()->getNumVariables(), a.getNode(), mcard);
  }
  void compute(int k, node_handle a, mpz_object &b);
};

void MEDDLY::card_mdd_mpz::compute(int k, node_handle a, mpz_object &card)
{
  // Terminal cases
  if (0==a) {
    card.setValue(0);
    return;
  }
  if (0==k) {
    card.setValue(1);
    return;
  }

  // Quickly deal with skipped levels
  if (argF->getNodeLevel(a) < k) {
    // skipped level
    compute(k-1, a, card);
    card.multiply(argF->getLevelSize(k));
    return;
  }
  
  // Check compute table
  CTsrch.key(0) = a;
  const node_handle* cacheEntry = CT->find(CTsrch);
  if (cacheEntry) {
    mpz_object* answer;
    memcpy(&answer, cacheEntry+1, sizeof(mpz_object*));
    answer->copyInto(card);
    return;
  }

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, false);

  // Recurse
  mpz_object tmp;
  tmp.setValue(0);
  card.setValue(0);
  int kdn = k-1;
  for (int z=0; z<A->getNNZs(); z++) {
    compute(kdn, A->d(z), tmp);
    card.add(tmp);
  }
  
  // Cleanup
  node_reader::recycle(A);

  // Add entry to compute table
  compute_table::temp_entry &entry = CT->startNewEntry(this);
  entry.key(0) = argF->cacheNode(a);
  mpz_object* answer = new mpz_object(card);
  entry.copyResult(0, &answer, sizeof(mpz_object*));
  CT->addEntry();

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is ", a);
  card.show(stderr);
  fprintf(stderr, "\n");
#endif
}

#endif

// ******************************************************************
// *                                                                *
// *                       card_mxd_mpz class                       *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP

/// Cardinality of MxDs, returning large (mpz) integers.
class MEDDLY::card_mxd_mpz : public card_mpz {
public:
  card_mxd_mpz(const unary_opname* oc, expert_forest* arg)
    : card_mpz(oc, arg) { }
  virtual void compute(const dd_edge& a, ct_object &res) {
    mpz_object& mcard = dynamic_cast <mpz_object &> (res);
    compute(argF->getDomain()->getNumVariables(), a.getNode(), mcard);
  }
  void compute(int k, node_handle a, mpz_object &b);
};

void MEDDLY::card_mxd_mpz::compute(int k, node_handle a, mpz_object &card)
{
  // Terminal cases
  if (0==a) {
    card.setValue(0);
    return;
  }
  if (0==k) {
    card.setValue(1);
    return;
  }

  // Quickly deal with skipped levels
  if (isLevelAbove(k, argF->getNodeLevel(a))) {
    if (k<0 && argF->isIdentityReduced()) {
      // identity node
      compute(-k-1, a, card);
      return;
    }
    // redundant node
    compute(-k, a, card);
    card.multiply(argF->getLevelSize(k));
    return;
  }
  
  // Check compute table
  CTsrch.key(0) = a;
  const node_handle* cacheEntry = CT->find(CTsrch);
  if (cacheEntry) {
    mpz_object* answer;
    memcpy(&answer, cacheEntry+1, sizeof(mpz_object*));
    answer->copyInto(card);
    return;
  }

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, false);

  // Recurse
  mpz_object tmp;
  tmp.setValue(0);
  card.setValue(0);
  int kdn = (k<0) ? -k-1 : -k;
  for (int z=0; z<A->getNNZs(); z++) {
    compute(kdn, A->d(z), tmp);
    card.add(tmp);
  }
  
  // Cleanup
  node_reader::recycle(A);

  // Add entry to compute table
  compute_table::temp_entry &entry = CT->startNewEntry(this);
  entry.key(0) = argF->cacheNode(a);
  mpz_object* answer = new mpz_object(card);
  entry.copyResult(0, &answer, sizeof(mpz_object*));
  CT->addEntry();

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is ", a);
  card.show(stderr);
  fprintf(stderr, "\n");
#endif
}

#endif



// ******************************************************************
// *                                                                *
// *                       card_opname  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::card_opname : public unary_opname {
  public:
    card_opname();
    virtual unary_operation* 
      buildOperation(expert_forest* ar, opnd_type res) const;
};

MEDDLY::card_opname::card_opname()
 : unary_opname("Card")
{
}

MEDDLY::unary_operation* 
MEDDLY::card_opname::buildOperation(expert_forest* arg, opnd_type res) const
{
  if (0==arg) return 0;
  switch (res) {
    case INTEGER:
      if (arg->isForRelations())
        return new card_mxd_int(this, arg);
      else                        
        return new card_mdd_int(this, arg);

    case REAL:
      if (arg->isForRelations())
        return new card_mxd_real(this, arg);
      else                        
        return new card_mdd_real(this, arg);

#ifdef HAVE_LIBGMP
    case HUGEINT:
      if (arg->isForRelations())
        return new card_mxd_mpz(this, arg);
      else                        
        return new card_mdd_mpz(this, arg);
#endif

    default:
      throw error(error::TYPE_MISMATCH);
  }
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeCardinality(const settings &s)
{
  return new card_opname;
}


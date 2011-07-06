
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
#include "../compute_table.h"
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
  virtual bool isEntryStale(const int* entryData);
  virtual void discardEntry(const int* entryData);
  virtual void showEntry(FILE* strm, const int *entryData) const;

protected:
  static inline void overflow_acc(long &a, long x) {
    a += x;
    if (a < x) throw error(error::OVERFLOW);
  }

  compute_table::search_key CTsrch;
};

MEDDLY::card_int::card_int(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, true, arg, INTEGER)
{
  key_length = 1;
  ans_length = sizeof(long) / sizeof(int); 
  CT->initializeSearchKey(CTsrch, this);
}

bool MEDDLY::card_int::isEntryStale(const int* data)
{
  return argF->isStale(data[0]);
}

void MEDDLY::card_int::discardEntry(const int* data)
{
  argF->uncacheNode(data[0]);
}

void MEDDLY::card_int::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d): %ld(L)]",
      getName(), data[0], ((const long*)(data+1))[0]
  );
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
  long compute(int ht, int a);
};

long MEDDLY::card_mdd_int::compute(int ht, int a)
{
  if (0==a) return 0;
  if (argF->getNodeHeight(a) < ht) {
    // skipped level
    long card = compute(ht-1, a);
    long ans = card * argF->getLevelSize(ht);
    if (ans < card) 
      throw error(error::OVERFLOW);
    return ans;
  }
  
  // Terminal case
  if (argF->isTerminalNode(a)) return 1;

  // Check compute table
  CTsrch.key(0) = a; 
  const int* cacheEntry = CT->find(CTsrch);
  if (cacheEntry) {
    // ugly but portable
    long answer;
    memcpy(&answer, cacheEntry+1, sizeof(long));
    return answer;
  }

  // recurse
  long card = 0;
  if (argF->isFullNode(a)) {
    int asize = argF->getFullNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      overflow_acc(card, compute(ht - 1, argF->getFullNodeDownPtr(a, i)));
    } // for i
  } else {
    int asize = argF->getSparseNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      overflow_acc(card, compute(ht - 1, argF->getSparseNodeDownPtr(a, i)));
    }
  }

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
    res = compute(argF->getDomain()->getNumVariables(), false, arg.getNode());
  }
  long compute(int ht, bool primed, int a);
};

long MEDDLY::card_mxd_int::compute(int ht, bool primed, int a)
{
  if (0==a) return 0;
  if (argF->getNodeHeight(a) < ht) {
    // skipped level
    long card;
    if (argF->getReductionRule() == forest::IDENTITY_REDUCED) {
      assert(!primed);
      card = compute(ht-1, false, a);
    } else {
      card = primed 
                ? compute(ht-1, false, a)
                : compute(ht, true, a);
    }
    long ans = card * argF->getLevelSize(ht);
    if (ans < card) 
      throw error(error::OVERFLOW);
    return ans;
  }

  // Terminal case
  if (argF->isTerminalNode(a)) return 1;

  // Check compute table
  const int* cacheEntry = CT->find(this, &a);
  if (cacheEntry) {
    return ((const long*)(cacheEntry+1))[0];
  }

// mdd version
  
  // recurse
  int htm1 = primed ? ht-1 : ht;
  bool nextpr = !primed;
  long card = 0;
  if (argF->isFullNode(a)) {
    int asize = argF->getFullNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      overflow_acc(
        card, compute(htm1, nextpr, argF->getFullNodeDownPtr(a, i))
      );
    } // for i
  } else {
    int asize = argF->getSparseNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      overflow_acc(
        card, compute(htm1, nextpr, argF->getSparseNodeDownPtr(a, i))
      );
    }
  }

  // Add entry to compute table
  static int ansEntry[1+sizeof(long)/sizeof(int)];
  ansEntry[0] = argF->cacheNode(a);
  ((long*)(ansEntry+1))[0] = card;

  CT->add(this, const_cast<const int*>(ansEntry));
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
  virtual bool isEntryStale(const int* entryData);
  virtual void discardEntry(const int* entryData);
  virtual void showEntry(FILE* strm, const int *entryData) const;
};

MEDDLY::card_real::card_real(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, true, arg, REAL)
{
  key_length = 1;
  ans_length = sizeof(double) / sizeof(int); 
}

bool MEDDLY::card_real::isEntryStale(const int* data)
{
  return argF->isStale(data[0]);
}

void MEDDLY::card_real::discardEntry(const int* data)
{
  argF->uncacheNode(data[0]);
}

void MEDDLY::card_real::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d): %le]",
      getName(), data[0], ((const double*)(data+1))[0]
  );
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
  double compute(int ht, int a);
};

double MEDDLY::card_mdd_real::compute(int ht, int a)
{
  if (0==a) return 0.0;
  if (argF->getNodeHeight(a) < ht) {
    // skipped level
    return argF->getLevelSize(ht) * compute(ht - 1, a);
  }
  
  // Terminal case
  if (argF->isTerminalNode(a)) return 1.0;

  // Check compute table
  const int* cacheEntry = CT->find(this, &a);
  if (cacheEntry) {
    return ((const double*)(cacheEntry+1))[0];
  }

  // recurse
  double card = 0.0;
  if (argF->isFullNode(a)) {
    int asize = argF->getFullNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      card += compute(ht - 1, argF->getFullNodeDownPtr(a, i));
    } // for i
  } else {
    int asize = argF->getSparseNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      card += compute(ht - 1, argF->getSparseNodeDownPtr(a, i));
    }
  }

  // Add entry to compute table
  static int ansEntry[1+sizeof(long)/sizeof(int)];
  ansEntry[0] = argF->cacheNode(a);
  ((double*)(ansEntry+1))[0] = card;

  CT->add(this, const_cast<const int*>(ansEntry));
#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %le\n", a, card);
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
    res = compute(argF->getDomain()->getNumVariables(), false, arg.getNode());
  }
  double compute(int ht, bool primed, int a);
};

double MEDDLY::card_mxd_real::compute(int ht, bool primed, int a)
{
  if (0==a) return 0.0;
  if (argF->getNodeHeight(a) < ht) {
    // skipped level
    double card;
    if (argF->getReductionRule() == forest::IDENTITY_REDUCED) {
      assert(!primed);
      card = compute(ht-1, false, a);
    } else {
      card = primed 
                ? compute(ht-1, false, a)
                : compute(ht, true, a);
    }
    return card * argF->getLevelSize(ht);
  }
  
  // Terminal case
  if (argF->isTerminalNode(a)) return 1.0;

  // Check compute table
  const int* cacheEntry = CT->find(this, &a);
  if (cacheEntry) {
    return ((const double*)(cacheEntry+1))[0];
  }

  // recurse
  int htm1 = primed ? ht-1 : ht;
  bool nextprimed = !primed;
  double card = 0.0;
  if (argF->isFullNode(a)) {
    int asize = argF->getFullNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      card += compute(htm1, nextprimed, argF->getFullNodeDownPtr(a, i));
    } // for i
  } else {
    int asize = argF->getSparseNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      card += compute(htm1, nextprimed, argF->getSparseNodeDownPtr(a, i));
    }
  }

  // Add entry to compute table
  static int ansEntry[1+sizeof(long)/sizeof(int)];
  ansEntry[0] = argF->cacheNode(a);
  ((double*)(ansEntry+1))[0] = card;

  CT->add(this, const_cast<const int*>(ansEntry));
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
  virtual bool isEntryStale(const int* entryData);
  virtual void discardEntry(const int* entryData);
  virtual void showEntry(FILE* strm, const int *entryData) const;
};

MEDDLY::card_mpz::card_mpz(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, true, arg, HUGEINT)
{
  key_length = 1;
  ans_length = sizeof(ct_object*) / sizeof(int); 
}

bool MEDDLY::card_mpz::isEntryStale(const int* data)
{
  return argF->isStale(data[0]);
}

void MEDDLY::card_mpz::discardEntry(const int* data)
{
  argF->uncacheNode(data[0]);
  delete( ((ct_object**)(data+1))[0] );
}

void MEDDLY::card_mpz::showEntry(FILE* strm, const int *entryData) const
{
  fprintf(strm, "[%s(%d): ", getName(), entryData[0]);
  ((const mpz_object**)(entryData+1))[0]->show(strm);
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
  void compute(int ht, int a, mpz_object &b);
};

void MEDDLY::card_mdd_mpz::compute(int ht, int a, mpz_object &card)
{
  if (0==a) {
    card.setValue(0);
    return;
  }
  if (argF->getNodeHeight(a) < ht) {
    // skipped level
    compute(ht-1, a, card);
    card.multiply(argF->getLevelSize(ht));
    return;
  }
  
  // Terminal case
  if (argF->isTerminalNode(a)) {
    card.setValue(1);
    return;
  }

  // Check compute table
  const int* cacheEntry = CT->find(this, &a);
  if (cacheEntry) {
    ((const mpz_object**)(cacheEntry+1))[0]->copyInto(card);
    return;
  }

  // recurse
  mpz_object tmp;
  tmp.setValue(0);
  card.setValue(0);
  if (argF->isFullNode(a)) {
    int asize = argF->getFullNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      compute(ht - 1, argF->getFullNodeDownPtr(a, i), tmp);
      card.add(tmp);
    } // for i
  } else {
    int asize = argF->getSparseNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      compute(ht - 1, argF->getSparseNodeDownPtr(a, i), tmp);
      card.add(tmp);
    }
  }

  // Add entry to compute table
  static int ansEntry[1+sizeof(mpz_t)/sizeof(int)];
  ansEntry[0] = argF->cacheNode(a);
  ((mpz_object**)(ansEntry+1))[0] = new mpz_object(card);

  CT->add(this, const_cast<const int*>(ansEntry));
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
    compute(argF->getDomain()->getNumVariables(), false, a.getNode(), mcard);
  }
  void compute(int ht, bool primed, int a, mpz_object &b);
};

void MEDDLY::card_mxd_mpz::compute(int ht, bool primed, int a, mpz_object &card)
{
  if (0==a) {
    card.setValue(0);
    return;
  }
  if (argF->getNodeHeight(a) < ht) {
    // skipped level
    if (argF->getReductionRule() == forest::IDENTITY_REDUCED) {
      assert(!primed);
      compute(ht-1, false, a, card);
    } else {
      if (primed) compute(ht-1, false, a, card);
      else        compute(ht, true, a, card);
    }
    card.multiply(argF->getLevelSize(ht));
    return;
  }
  
  // Terminal case
  if (argF->isTerminalNode(a)) {
    card.setValue(1);
    return;
  }

  // Check compute table
  const int* cacheEntry = CT->find(this, &a);
  if (cacheEntry) {
    ((const mpz_object**)(cacheEntry+1))[0]->copyInto(card);
    return;
  }

  // recurse
  int htm1 = primed ? ht-1 : ht;
  bool nextprimed = !primed;
  mpz_object tmp;
  tmp.setValue(0);
  card.setValue(0);
  if (argF->isFullNode(a)) {
    int asize = argF->getFullNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      compute(htm1, nextprimed, argF->getFullNodeDownPtr(a, i), tmp);
      card.add(tmp);
    } // for i
  } else {
    int asize = argF->getSparseNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      compute(htm1, nextprimed, argF->getSparseNodeDownPtr(a, i), tmp);
      card.add(tmp);
    }
  }

  // Add entry to compute table
  static int ansEntry[1+sizeof(mpz_t)/sizeof(int)];
  ansEntry[0] = argF->cacheNode(a);
  ((mpz_object**)(ansEntry+1))[0] = new mpz_object(card);

  CT->add(this, const_cast<const int*>(ansEntry));
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



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

// for now
#include "../compute_cache.h"

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
};

MEDDLY::card_int::card_int(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, arg, INTEGER)
{
  key_length = 1;
  ans_length = sizeof(long) / sizeof(int); 
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
  const int* cacheEntry = CT->find(this, &a);
  if (cacheEntry) {
    return ((const long*)(cacheEntry+1))[0];
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
  static int ansEntry[1+sizeof(long)/sizeof(int)];
  argF->cacheNode(a);
  ansEntry[0] = a;
  ((long*)(ansEntry+1))[0] = card;

  CT->add(this, const_cast<const int*>(ansEntry));
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
  argF->cacheNode(a);
  ansEntry[0] = a;
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
 : unary_operation(oc, arg, REAL)
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
  argF->cacheNode(a);
  ansEntry[0] = a;
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
  argF->cacheNode(a);
  ansEntry[0] = a;
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
 : unary_operation(oc, arg, HUGEINT)
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
  argF->cacheNode(a);
  ansEntry[0] = a;
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
  argF->cacheNode(a);
  ansEntry[0] = a;
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
// *                       Card_opname  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::card_opname : public unary_opname {
  public:
    card_opname();
    virtual unary_operation* 
      buildOperation(const forest* ar, opnd_type res) const;
};

MEDDLY::card_opname::card_opname()
 : unary_opname("Card")
{
}

MEDDLY::unary_operation* 
MEDDLY::card_opname::buildOperation(const forest* arg, opnd_type res) const
{
  if (0==arg) return 0;
  switch (res) {
    case INTEGER:
      if (arg->isForRelations())
        return new card_mxd_int(this, (expert_forest*)arg);
      else                        
        return new card_mdd_int(this, (expert_forest*)arg);

    case REAL:
      if (arg->isForRelations())
        return new card_mxd_real(this, (expert_forest*)arg);
      else                        
        return new card_mdd_real(this, (expert_forest*)arg);

#ifdef HAVE_LIBGMP
    case HUGEINT:
      if (arg->isForRelations())
        return new card_mxd_mpz(this, (expert_forest*)arg);
      else                        
        return new card_mdd_mpz(this, (expert_forest*)arg);
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

const MEDDLY::unary_opname* MEDDLY::initializeCardinality(const settings &s)
{
  return new card_opname;
}

// ******************************************************************
// ******************************************************************
// ******************************************************************
// ***                         OLD  STUFF                         ***
// ******************************************************************
// ******************************************************************
// ******************************************************************

namespace MEDDLY {

// ******************************************************************
// *                                                                *
// *                      cardinality_op class                      *
// *                                                                *
// ******************************************************************

/** Abstract base class for cardinality operations.
*/
class cardinality_op : public old_operation {
  public:
    cardinality_op()                        { }
    virtual ~cardinality_op()               { }

    virtual const char* getName() const     { return "Cardinality"; }
    virtual int getKeyLength() const        { return 1; }
    virtual int getKeyLengthInBytes() const { return sizeof(int); }

    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);

  protected:
    inline void 
    type_check(const op_info* owner, op_param::type p1type) 
    {
        if (owner == 0)
          throw error(error::UNKNOWN_OPERATION);
        if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
          throw error(error::TYPE_MISMATCH);
        if (owner->nParams != 2)
          throw error(error::WRONG_NUMBER);
        if (owner->p[1].getType() != p1type)
          throw error(error::TYPE_MISMATCH);
    }
};

bool cardinality_op::
isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[0] <--> DD node to take cardinality of
  // data[1] <--> result value (not a DD node)
  DCASSERT(owner->nParams == 2);
  return owner->p[0].getForest()->isStale(data[0]);
}

void cardinality_op::
discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[0] <--> DD node to take cardinality of
  // data[1] <--> result value (not a DD node)
  DCASSERT(owner->nParams == 2);
  owner->p[0].getForest()->uncacheNode(data[0]);
}

// ******************************************************************
// *                                                                *
// *                       mdd_int_card class                       *
// *                                                                *
// ******************************************************************

/// Cardinality of MDDs, returning integers.
class mdd_int_card : public cardinality_op {
  public:
    mdd_int_card()          { }
    virtual ~mdd_int_card() { }

    virtual int getAnsLength() const  { 
      return sizeof(long) / sizeof(int); 
    }
    virtual int getCacheEntryLength() const { 
      return 1 + sizeof(long) / sizeof(int);
    }

    virtual int getAnsLengthInBytes() const { 
      return sizeof(long); 
    }
    virtual int getCacheEntryLengthInBytes() const { 
      return sizeof(int) + sizeof(long);
    }

    static mdd_int_card* getInstance();

    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    virtual void typeCheck(const op_info* owner) {
      if (owner->p[0].getForest()->isForRelations())
        throw error(error::TYPE_MISMATCH);
      type_check(owner, op_param::INTEGER);
    }

    virtual void compute(op_info* owner, const dd_edge& a,
      long& b);

    long compute(op_info* owner, int ht, int a);
};

mdd_int_card* 
mdd_int_card::
getInstance()
{
  static mdd_int_card instance;
  return &instance;
}

void
mdd_int_card::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  DCASSERT(owner->nParams == 2);
  fprintf(strm, "[%s(%d): %ld(L)]",
      owner->op->getName(), data[0], ((const long*)(data+1))[0]
  );
}

void 
mdd_int_card::
compute(op_info* owner, const dd_edge& a, long& b)
{
  if (0==owner) 
    throw error(error::TYPE_MISMATCH);
  expert_forest* f = owner->p[0].getForest();
  b = compute(owner, f->getDomain()->getNumVariables(), a.getNode());
  return;
}

long mdd_int_card::compute(op_info* owner, int ht, int a)
{
  if (0==a) return 0;
  expert_forest* f = owner->p[0].getForest();
  if (f->getNodeHeight(a) < ht) {
    // skipped level
    const expert_domain* d = smart_cast <const expert_domain*> (f->getDomain());
    long card = compute(owner, ht-1, a);
    if (card<0) return -1;
    long ans = card * f->getLevelSize(ht);
    if (ans < card) return -1; // overflow
    return ans;
  }
  
  // Terminal case
  if (f->isTerminalNode(a)) return 1;

  // Check compute table
  const int* cacheEntry = owner->cc->find(owner, &a);
  if (cacheEntry) {
    return ((const long*)(cacheEntry+1))[0];
  }

  // recurse
  bool a_is_full = f->isFullNode(a);
  int asize = a_is_full ? f->getFullNodeSize(a) : f->getSparseNodeSize(a);
  long card = 0;
  for (int i = 0; i < asize; ++i) {
    long d = a_is_full 
              ? compute(owner, ht - 1, f->getFullNodeDownPtr(a, i))
              : compute(owner, ht - 1, f->getSparseNodeDownPtr(a, i));

    // tricky, because we need to keep any overflows at -1
    d += card;
    if (d < card) { // must be overflow somewhere
      card = -1;
      break;
    } else {
      card = d;
    }
  } // for i

  // Add entry to compute table
  static int ansEntry[1+sizeof(long)/sizeof(int)];
  owner->p[0].getForest()->cacheNode(a);
  ansEntry[0] = a;
  ((long*)(ansEntry+1))[0] = card;

  owner->cc->add(owner, const_cast<const int*>(ansEntry));
#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %ld(L)\n", a, card);
#endif
  return card;
}

// ******************************************************************
// *                                                                *
// *                      mdd_real_card  class                      *
// *                                                                *
// ******************************************************************

/// Cardinality of MDDs, returning reals.
class mdd_real_card : public cardinality_op {
  public:
    mdd_real_card()          { }
    virtual ~mdd_real_card() { }

    virtual int getAnsLength() const  { 
      return sizeof(double) / sizeof(int); 
    }
    virtual int getCacheEntryLength() const { 
      return 1 + sizeof(double) / sizeof(int);
    }

    virtual int getAnsLengthInBytes() const { 
      return sizeof(double); 
    }
    virtual int getCacheEntryLengthInBytes() const { 
      return sizeof(int) + sizeof(double);
    }

    static mdd_real_card* getInstance();

    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    virtual void typeCheck(const op_info* owner) {
      if (owner->p[0].getForest()->isForRelations())
        throw error(error::TYPE_MISMATCH);
      type_check(owner, op_param::REAL);
    }

    virtual void compute(op_info* owner, const dd_edge& a,
      double& b);

    double compute(op_info* owner, int ht, int a);
};

mdd_real_card* 
mdd_real_card::
getInstance()
{
  static mdd_real_card instance;
  return &instance;
}

void
mdd_real_card::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  DCASSERT(owner->nParams == 2);
  fprintf(strm, "[%s(%d): %le]",
      owner->op->getName(), data[0], ((const double*)(data+1))[0]
  );
}

void 
mdd_real_card::
compute(op_info* owner, const dd_edge& a, double& b)
{
  if (0==owner) 
    throw error(error::TYPE_MISMATCH);
  expert_forest* f = owner->p[0].getForest();
  b = compute(owner, f->getDomain()->getNumVariables(), a.getNode());
}

double mdd_real_card::compute(op_info* owner, int ht, int a)
{
  if (0==a) return 0.0;
  expert_forest* f = owner->p[0].getForest();
  if (f->getNodeHeight(a) < ht) {
    // skipped level
    const expert_domain* d = smart_cast <const expert_domain*> (f->getDomain());
    return f->getLevelSize(ht) * compute(owner, ht - 1, a);
  }
  
  // Terminal case
  if (f->isTerminalNode(a)) return 1.0;

  // Check compute table
  const int* cacheEntry = owner->cc->find(owner, &a);
  if (cacheEntry) {
    return ((const double*)(cacheEntry+1))[0];
  }

  // recurse
  double card = 0.0;
  if (f->isFullNode(a)) {
    int asize = f->getFullNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      card += compute(owner, ht - 1, f->getFullNodeDownPtr(a, i));
    } // for i
  } else {
    int asize = f->getSparseNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      card += compute(owner, ht - 1, f->getSparseNodeDownPtr(a, i));
    }
  }
  
  // Add entry to compute table
  static int ansEntry[1+sizeof(double)/sizeof(int)];
  owner->p[0].getForest()->cacheNode(a);
  ansEntry[0] = a;
  ((double*)(ansEntry+1))[0] = card;

  owner->cc->add(owner, const_cast<const int*>(ansEntry));
#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %le\n", a, card);
#endif
  return card;
}

// ******************************************************************
// *                                                                *
// *                       mdd_mpz_card class                       *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP

/// Cardinality of MDDs, returning large (mpz) integers.
class mdd_mpz_card : public cardinality_op {
  public:
    mdd_mpz_card()          { }
    virtual ~mdd_mpz_card() { }

    virtual int getAnsLength() const  { 
      return sizeof(ct_object*) / sizeof(int); 
    }
    virtual int getCacheEntryLength() const { 
      return 1 + sizeof(ct_object*) / sizeof(int);
    }

    virtual int getAnsLengthInBytes() const { 
      return sizeof(ct_object*); 
    }
    virtual int getCacheEntryLengthInBytes() const { 
      return sizeof(int) + sizeof(ct_object*);
    }

    static mdd_mpz_card* getInstance();

    virtual void discardEntry(op_info* owner, const int* entryData);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    virtual void typeCheck(const op_info* owner) {
      if (owner->p[0].getForest()->isForRelations())
        throw error(error::TYPE_MISMATCH);
      type_check(owner, op_param::HUGEINT);
    }

    virtual void compute(op_info* owner, const dd_edge& a,
      ct_object &b);

    void compute(op_info* owner, int ht, int a, mpz_object &b);
};

mdd_mpz_card* 
mdd_mpz_card::
getInstance()
{
  static mdd_mpz_card instance;
  return &instance;
}

void mdd_mpz_card::
discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[0] <--> DD node to take cardinality of
  // data[1] <--> result value (not a DD node)
  DCASSERT(owner->nParams == 2);
  owner->p[0].getForest()->uncacheNode(data[0]);
  delete( ((ct_object**)(data+1))[0] );
}

void
mdd_mpz_card::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  DCASSERT(owner->nParams == 2);
  fprintf(strm, "[%s(%d): ", owner->op->getName(), data[0]);
  ((const mpz_object**)(data+1))[0]->show(strm);
  fprintf(strm, "]");
}

void 
mdd_mpz_card::
compute(op_info* owner, const dd_edge& a, ct_object &card)
{
  if (0==owner) 
    throw error(error::TYPE_MISMATCH);
  expert_forest* f = owner->p[0].getForest();
  mpz_object& mcard = dynamic_cast <mpz_object &> (card);
  compute(owner, f->getDomain()->getNumVariables(), a.getNode(), mcard);
}

void mdd_mpz_card::compute(op_info* owner, int ht, int a, mpz_object &card)
{
  if (0==a) {
    card.setValue(0);
    return;
  }
  expert_forest* f = owner->p[0].getForest();
  if (f->getNodeHeight(a) < ht) {
    // skipped level
    const expert_domain* d = smart_cast <const expert_domain*> (f->getDomain());
    compute(owner, ht-1, a, card);
    card.multiply(f->getLevelSize(ht));
    return;
  }
  
  // Terminal case
  if (f->isTerminalNode(a)) {
    card.setValue(1);
    return;
  }

  // Check compute table
  const int* cacheEntry = owner->cc->find(owner, &a);
  if (cacheEntry) {
    ((const mpz_object**)(cacheEntry+1))[0]->copyInto(card);
    return;
  }

  // recurse
  mpz_object tmp;
  tmp.setValue(0);
  card.setValue(0);
  if (f->isFullNode(a)) {
    int asize = f->getFullNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      compute(owner, ht - 1, f->getFullNodeDownPtr(a, i), tmp);
      card.add(tmp);
    } // for i
  } else {
    int asize = f->getSparseNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      compute(owner, ht - 1, f->getSparseNodeDownPtr(a, i), tmp);
      card.add(tmp);
    }
  }

  // Add entry to compute table
  static int ansEntry[1+sizeof(mpz_t)/sizeof(int)];
  owner->p[0].getForest()->cacheNode(a);
  ansEntry[0] = a;
  ((mpz_object**)(ansEntry+1))[0] = new mpz_object(card);

  owner->cc->add(owner, const_cast<const int*>(ansEntry));
#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is ", a);
  card.show(stderr);
  fprintf(stderr, "\n");
#endif
}

#endif

// ******************************************************************
// *                                                                *
// *                       mxd_int_card class                       *
// *                                                                *
// ******************************************************************

/// Cardinality of MxDs, returning integers.
class mxd_int_card : public cardinality_op {
  public:
    mxd_int_card()          { }
    virtual ~mxd_int_card() { }

    virtual int getAnsLength() const  { 
      return sizeof(long) / sizeof(int); 
    }
    virtual int getCacheEntryLength() const { 
      return 1 + sizeof(long) / sizeof(int);
    }

    virtual int getAnsLengthInBytes() const { 
      return sizeof(long); 
    }
    virtual int getCacheEntryLengthInBytes() const { 
      return sizeof(int) + sizeof(long);
    }

    static mxd_int_card* getInstance();

    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    virtual void typeCheck(const op_info* owner) {
      if (!owner->p[0].getForest()->isForRelations())
        throw error(error::TYPE_MISMATCH);
      type_check(owner, op_param::INTEGER);
    }

    virtual void compute(op_info* owner, const dd_edge& a,
      long& b);

    long compute(op_info* owner, int ht, bool primed, int a);
};

mxd_int_card* 
mxd_int_card::
getInstance()
{
  static mxd_int_card instance;
  return &instance;
}

void
mxd_int_card::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  DCASSERT(owner->nParams == 2);
  fprintf(strm, "[%s(%d): %ld(L)]",
      owner->op->getName(), data[0], ((const long*)(data+1))[0]
  );
}

void 
mxd_int_card::
compute(op_info* owner, const dd_edge& a, long& b)
{
  if (0==owner) 
    throw error(error::TYPE_MISMATCH);
  expert_forest* f = owner->p[0].getForest();
  b = compute(owner, f->getDomain()->getNumVariables(), false, a.getNode());
}

long mxd_int_card::compute(op_info* owner, int ht, bool primed, int a)
{
  if (0==a) return 0;
  expert_forest* f = owner->p[0].getForest();
  if (f->getNodeHeight(a) < ht) {
    // skipped level
    long card;
    if (f->getReductionRule() == forest::IDENTITY_REDUCED) {
      assert(!primed);
      card = compute(owner, ht-1, false, a);
    } else {
      card = primed 
                ? compute(owner, ht-1, false, a)
                : compute(owner, ht, true, a);
    }
    if (card<0) return -1;
    const expert_domain* d = smart_cast <const expert_domain*> (f->getDomain());
    long ans = card * f->getLevelSize(ht);
    if (ans < card) return -1; // overflow
    return ans;
  }
  
  // Terminal case
  if (f->isTerminalNode(a)) return 1;

  // Check compute table
  const int* cacheEntry = owner->cc->find(owner, &a);
  if (cacheEntry) {
    return ((const long*)(cacheEntry+1))[0];
  }

  // recurse
  bool a_is_full = f->isFullNode(a);
  int asize = a_is_full ? f->getFullNodeSize(a) : f->getSparseNodeSize(a);
  long card = 0;
  int htm1 = primed ? ht-1 : ht;
  bool nextprimed = !primed;
  for (int i = 0; i < asize; ++i) {
    long d = a_is_full 
              ? compute(owner, htm1, nextprimed, f->getFullNodeDownPtr(a, i))
              : compute(owner, htm1, nextprimed, f->getSparseNodeDownPtr(a, i));

    // tricky, because we need to keep any overflows at -1
    d += card;
    if (d < card) { // must be overflow somewhere
      card = -1;
      break;
    } else {
      card = d;
    }
  } // for i

  // Add entry to compute table
  static int ansEntry[1+sizeof(long)/sizeof(int)];
  owner->p[0].getForest()->cacheNode(a);
  ansEntry[0] = a;
  ((long*)(ansEntry+1))[0] = card;

  owner->cc->add(owner, const_cast<const int*>(ansEntry));
#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %ld(L)\n", a, card);
#endif
  return card;
}

// ******************************************************************
// *                                                                *
// *                      mxd_real_card  class                      *
// *                                                                *
// ******************************************************************

/// Cardinality of MDDs, returning reals.
class mxd_real_card : public cardinality_op {
  public:
    mxd_real_card()          { }
    virtual ~mxd_real_card() { }

    virtual int getAnsLength() const  { 
      return sizeof(double) / sizeof(int); 
    }
    virtual int getCacheEntryLength() const { 
      return 1 + sizeof(double) / sizeof(int);
    }

    virtual int getAnsLengthInBytes() const { 
      return sizeof(double); 
    }
    virtual int getCacheEntryLengthInBytes() const { 
      return sizeof(int) + sizeof(double);
    }

    static mxd_real_card* getInstance();

    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    virtual void typeCheck(const op_info* owner) {
      if (!owner->p[0].getForest()->isForRelations())
        throw error(error::TYPE_MISMATCH);
      type_check(owner, op_param::REAL);
    }

    virtual void compute(op_info* owner, const dd_edge& a,
      double& b);

    double compute(op_info* owner, int ht, bool primed, int a);
};

mxd_real_card* 
mxd_real_card::
getInstance()
{
  static mxd_real_card instance;
  return &instance;
}

void
mxd_real_card::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  DCASSERT(owner->nParams == 2);
  fprintf(strm, "[%s(%d): %le]",
      owner->op->getName(), data[0], ((const double*)(data+1))[0]
  );
}

void 
mxd_real_card::
compute(op_info* owner, const dd_edge& a, double& b)
{
  if (0==owner) 
    throw error(error::TYPE_MISMATCH);
  expert_forest* f = owner->p[0].getForest();
  b = compute(owner, f->getDomain()->getNumVariables(), false, a.getNode());
}

double mxd_real_card::compute(op_info* owner, int ht, bool primed, int a)
{
  if (0==a) return 0.0;
  expert_forest* f = owner->p[0].getForest();
  if (f->getNodeHeight(a) < ht) {
    // skipped level
    double card;
    if (f->getReductionRule() == forest::IDENTITY_REDUCED) {
      assert(!primed);
      card = compute(owner, ht-1, false, a);
    } else {
      card = primed 
                ? compute(owner, ht-1, false, a)
                : compute(owner, ht, true, a);
    }
    const expert_domain* d = smart_cast <const expert_domain*> (f->getDomain());
    return card * f->getLevelSize(ht);
  }
  
  // Terminal case
  if (f->isTerminalNode(a)) return 1.0;

  // Check compute table
  const int* cacheEntry = owner->cc->find(owner, &a);
  if (cacheEntry) {
    return ((const double*)(cacheEntry+1))[0];
  }

  // recurse
  double card = 0.0;
  int htm1 = primed ? ht-1 : ht;
  bool nextprimed = !primed;
  if (f->isFullNode(a)) {
    int asize = f->getFullNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      card += compute(owner, htm1, nextprimed, f->getFullNodeDownPtr(a, i));
    } // for i
  } else {
    int asize = f->getSparseNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      card += compute(owner, htm1, nextprimed, f->getSparseNodeDownPtr(a, i));
    }
  }
  
  // Add entry to compute table
  static int ansEntry[1+sizeof(double)/sizeof(int)];
  owner->p[0].getForest()->cacheNode(a);
  ansEntry[0] = a;
  ((double*)(ansEntry+1))[0] = card;

  owner->cc->add(owner, const_cast<const int*>(ansEntry));
#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %le\n", a, card);
#endif
  return card;
}

// ******************************************************************
// *                                                                *
// *                       mxd_mpz_card class                       *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP

/// Cardinality of MDDs, returning large (mpz) integers.
class mxd_mpz_card : public cardinality_op {
  public:
    mxd_mpz_card()          { }
    virtual ~mxd_mpz_card() { }

    virtual int getAnsLength() const  { 
      return sizeof(ct_object*) / sizeof(int); 
    }
    virtual int getCacheEntryLength() const { 
      return 1 + sizeof(ct_object*) / sizeof(int);
    }

    virtual int getAnsLengthInBytes() const { 
      return sizeof(ct_object*); 
    }
    virtual int getCacheEntryLengthInBytes() const { 
      return sizeof(int) + sizeof(ct_object*);
    }

    static mxd_mpz_card* getInstance();

    virtual void discardEntry(op_info* owner, const int* entryData);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    virtual void typeCheck(const op_info* owner) {
      if (!owner->p[0].getForest()->isForRelations())
        throw error(error::TYPE_MISMATCH);
      type_check(owner, op_param::HUGEINT);
    }

    virtual void compute(op_info* owner, const dd_edge& a,
      ct_object &b);

    void compute(op_info* owner, int ht, bool primed, int a, mpz_object &b);
};

mxd_mpz_card* 
mxd_mpz_card::
getInstance()
{
  static mxd_mpz_card instance;
  return &instance;
}

void mxd_mpz_card::
discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[0] <--> DD node to take cardinality of
  // data[1] <--> result value (not a DD node)
  DCASSERT(owner->nParams == 2);
  owner->p[0].getForest()->uncacheNode(data[0]);
  delete( ((ct_object**)(data+1))[0] );
}

void
mxd_mpz_card::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  DCASSERT(owner->nParams == 2);
  fprintf(strm, "[%s(%d): ", owner->op->getName(), data[0]);
  ((const mpz_object**)(data+1))[0]->show(strm);
  fprintf(strm, "]");
}

void 
mxd_mpz_card::
compute(op_info* owner, const dd_edge& a, ct_object &card)
{
  if (0==owner) 
    throw error(error::TYPE_MISMATCH);
  expert_forest* f = owner->p[0].getForest();
  mpz_object& mcard = dynamic_cast <mpz_object &> (card);
  compute(owner, f->getDomain()->getNumVariables(), false, a.getNode(), mcard);
}

void mxd_mpz_card::compute(op_info* owner, int ht, bool primed, int a, mpz_object &card)
{
  if (0==a) {
    card.setValue(0);
    return;
  }
  expert_forest* f = owner->p[0].getForest();
  if (f->getNodeHeight(a) < ht) {
    // skipped level
    if (f->getReductionRule() == forest::IDENTITY_REDUCED) {
      assert(!primed);
      compute(owner, ht-1, false, a, card);
    } else {
      if (primed) compute(owner, ht-1, false, a, card);
      else        compute(owner, ht, true, a, card);
    }
    const expert_domain* d = smart_cast <const expert_domain*> (f->getDomain());
    card.multiply(f->getLevelSize(ht));
    return;
  }
  
  // Terminal case
  if (f->isTerminalNode(a)) {
    card.setValue(1);
    return;
  }

  // Check compute table
  const int* cacheEntry = owner->cc->find(owner, &a);
  if (cacheEntry) {
    ((const mpz_object**)(cacheEntry+1))[0]->copyInto(card);
    return;
  }

  // recurse
  mpz_object tmp;
  tmp.setValue(0);
  card.setValue(0);
  int htm1 = primed ? ht-1 : ht;
  bool nextprimed = !primed;
  if (f->isFullNode(a)) {
    int asize = f->getFullNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      compute(owner, htm1, nextprimed, f->getFullNodeDownPtr(a, i), tmp);
      card.add(tmp);
    } // for i
  } else {
    int asize = f->getSparseNodeSize(a);
    for (int i = 0; i < asize; ++i) {
      compute(owner, htm1, nextprimed, f->getSparseNodeDownPtr(a, i), tmp);
      card.add(tmp);
    }
  }

  // Add entry to compute table
  static int ansEntry[1+sizeof(mpz_t)/sizeof(int)];
  owner->p[0].getForest()->cacheNode(a);
  ansEntry[0] = a;
  ((mpz_object**)(ansEntry+1))[0] = new mpz_object(card);

  owner->cc->add(owner, const_cast<const int*>(ansEntry));
#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is ", a);
  card.show(stderr);
  fprintf(stderr, "\n");
#endif
}

#endif

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

old_operation* getCardinalityOperation(const op_param &ft, const op_param &rt)
{
  if (!ft.isForest()) return 0;
  switch (rt.getType()) {
    case op_param::INTEGER:
        if (ft.readForest()->isForRelations()) {
          return mxd_int_card::getInstance();
        } else {
          return mdd_int_card::getInstance();
        }
        return 0;

    case op_param::REAL:
        if (ft.readForest()->isForRelations()) {
          return mxd_real_card::getInstance();
        } else {
          return mdd_real_card::getInstance();
        }
        return 0;

#ifdef HAVE_LIBGMP

    case op_param::HUGEINT:
        if (ft.readForest()->isForRelations()) {
          return mxd_mpz_card::getInstance();
        } else {
          return mdd_mpz_card::getInstance();
        }
        return 0;

#endif

    default:
        return 0;
  }
}

} // namespace MEDDLY


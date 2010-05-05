
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
#include "../compute_cache.h"
#include "mpz_object.h"

// #define DEBUG_CARD

// ******************************************************************
// *                                                                *
// *                      cardinality_op class                      *
// *                                                                *
// ******************************************************************

/** Abstract base class for cardinality operations.
*/
class cardinality_op : public operation {
  public:
    cardinality_op()                        { }
    virtual ~cardinality_op()               { }

    virtual const char* getName() const     { return "Cardinality"; }
    virtual int getKeyLength() const        { return 1; }
    virtual int getKeyLengthInBytes() const { return sizeof(int); }

    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);

  protected:
    inline compute_manager::error 
    type_check(const op_info* owner, bool mdd, op_param::type p1type) 
    {
        if (owner == 0)
          return compute_manager::UNKNOWN_OPERATION;
        if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
          return compute_manager::TYPE_MISMATCH;
        if (owner->nParams != 2)
          return compute_manager::WRONG_NUMBER;
        if (owner->p[0].isMdd()!=mdd || owner->p[1].getType() != p1type)
          return compute_manager::TYPE_MISMATCH;
        return compute_manager::SUCCESS;
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

    virtual compute_manager::error typeCheck(const op_info* owner) {
      return type_check(owner, true, op_param::INTEGER);
    }

    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
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

compute_manager::error 
mdd_int_card::
compute(op_info* owner, const dd_edge& a, long& b)
{
  if (0==owner) return compute_manager::TYPE_MISMATCH;
  expert_forest* f = owner->p[0].getForest();
  b = compute(owner, f->getDomain()->getNumVariables(), a.getNode());
  return compute_manager::SUCCESS;
}

long mdd_int_card::compute(op_info* owner, int ht, int a)
{
  if (0==a) return 0;
  expert_forest* f = owner->p[0].getForest();
  if (f->getNodeHeight(a) < ht) {
    // skipped level
    const expert_domain* d = smart_cast <const expert_domain*> (f->getDomain());
    int k = d->getVariableWithHeight(ht);
    long card = compute(owner, ht-1, a);
    if (card<0) return -1;
    long ans = card * f->getLevelSize(k);
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

    virtual compute_manager::error typeCheck(const op_info* owner) {
      return type_check(owner, true, op_param::REAL);
    }

    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
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

compute_manager::error 
mdd_real_card::
compute(op_info* owner, const dd_edge& a, double& b)
{
  if (0==owner) return compute_manager::TYPE_MISMATCH;
  expert_forest* f = owner->p[0].getForest();
  b = compute(owner, f->getDomain()->getNumVariables(), a.getNode());
  return compute_manager::SUCCESS;
}

double mdd_real_card::compute(op_info* owner, int ht, int a)
{
  if (0==a) return 0.0;
  expert_forest* f = owner->p[0].getForest();
  if (f->getNodeHeight(a) < ht) {
    // skipped level
    const expert_domain* d = smart_cast <const expert_domain*> (f->getDomain());
    int k = d->getVariableWithHeight(ht);
    return f->getLevelSize(k) * compute(owner, ht - 1, a);
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

    virtual compute_manager::error typeCheck(const op_info* owner) {
      return type_check(owner, true, op_param::HUGEINT);
    }

    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
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

compute_manager::error 
mdd_mpz_card::
compute(op_info* owner, const dd_edge& a, ct_object &card)
{
  if (0==owner) return compute_manager::TYPE_MISMATCH;
  expert_forest* f = owner->p[0].getForest();
  mpz_object& mcard = dynamic_cast <mpz_object &> (card);
  compute(owner, f->getDomain()->getNumVariables(), a.getNode(), mcard);
  return compute_manager::SUCCESS;
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
    int k = d->getVariableWithHeight(ht);
    compute(owner, ht-1, a, card);
    card.multiply(f->getLevelSize(k));
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

    virtual compute_manager::error typeCheck(const op_info* owner) {
      return type_check(owner, false, op_param::INTEGER);
    }

    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
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

compute_manager::error 
mxd_int_card::
compute(op_info* owner, const dd_edge& a, long& b)
{
  if (0==owner) return compute_manager::TYPE_MISMATCH;
  expert_forest* f = owner->p[0].getForest();
  b = compute(owner, f->getDomain()->getNumVariables(), false, a.getNode());
  return compute_manager::SUCCESS;
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
    int k = d->getVariableWithHeight(ht);
    long ans = card * f->getLevelSize(k);
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

    virtual compute_manager::error typeCheck(const op_info* owner) {
      return type_check(owner, false, op_param::REAL);
    }

    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
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

compute_manager::error 
mxd_real_card::
compute(op_info* owner, const dd_edge& a, double& b)
{
  if (0==owner) return compute_manager::TYPE_MISMATCH;
  expert_forest* f = owner->p[0].getForest();
  b = compute(owner, f->getDomain()->getNumVariables(), false, a.getNode());
  return compute_manager::SUCCESS;
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
    int k = d->getVariableWithHeight(ht);
    return card * f->getLevelSize(k);
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

    virtual compute_manager::error typeCheck(const op_info* owner) {
      return type_check(owner, false, op_param::HUGEINT);
    }

    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
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

compute_manager::error 
mxd_mpz_card::
compute(op_info* owner, const dd_edge& a, ct_object &card)
{
  if (0==owner) return compute_manager::TYPE_MISMATCH;
  expert_forest* f = owner->p[0].getForest();
  mpz_object& mcard = dynamic_cast <mpz_object &> (card);
  compute(owner, f->getDomain()->getNumVariables(), false, a.getNode(), mcard);
  return compute_manager::SUCCESS;
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
    int k = d->getVariableWithHeight(ht);
    card.multiply(f->getLevelSize(k));
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

operation* getCardinalityOperation(const op_param &ft, const op_param &rt)
{
  if (!ft.isForest()) return 0;
  switch (rt.getType()) {
    case op_param::INTEGER:
        if (ft.isMdd()) {
          if (ft.isMT())      return mdd_int_card::getInstance();
          if (ft.isEVPlus())  return 0;
          if (ft.isEVTimes()) return 0;
        } else {
          if (ft.isMT())      return mxd_int_card::getInstance();
          if (ft.isEVPlus())  return 0;
          if (ft.isEVTimes()) return 0;
        }
        return 0;

    case op_param::REAL:
        if (ft.isMdd()) {
          if (ft.isMT())      return mdd_real_card::getInstance();
          if (ft.isEVPlus())  return 0;
          if (ft.isEVTimes()) return 0;
        } else {
          if (ft.isMT())      return mxd_real_card::getInstance();
          if (ft.isEVPlus())  return 0;
          if (ft.isEVTimes()) return 0;
        }
        return 0;

#ifdef HAVE_LIBGMP

    case op_param::HUGEINT:
        if (ft.isMdd()) {
          if (ft.isMT())      return mdd_mpz_card::getInstance();
          if (ft.isEVPlus())  return 0;
          if (ft.isEVTimes()) return 0;
        } else {
          if (ft.isMT())      return mxd_mpz_card::getInstance();
          if (ft.isEVPlus())  return 0;
          if (ft.isEVTimes()) return 0;
        }
        return 0;

#endif

    default:
        return 0;
  }
}


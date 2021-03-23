
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
  // long* belowcount;
  // int lastNode;
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

protected:
  static inline void overflow_acc(long &a, long x) {
    a += x;
    if (a < x) throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
  }
  static inline long overflow_mult(long a, long x) {
    a *= x;
    if (a < x) throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
    return a;
  }
};

MEDDLY::card_int::card_int(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, 1, arg, INTEGER)
{
  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:L");
  et->setForestForSlot(0, arg);
  registerEntryType(0, et);
  buildCTs();
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
    res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode());
  }
  long compute_r(int k, node_handle a);
};

long MEDDLY::card_mdd_int::compute_r(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0;
  if (0==k) return 1;

  // Quickly deal with skipped levels
  if (argF->getNodeLevel(a) < k) {
    return overflow_mult(compute_r(k-1, a), argF->getLevelSize(k));
  }

  // Check compute table
  compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    CT0->recycle(CTsrch);
    return CTresult[0].readL();
  }

  // Initialize node reader
  unpacked_node* A = unpacked_node::newFromNode(argF, a, false);

  // Recurse
  long card = 0;
  int kdn = k-1;
  for (unsigned z=0; z<A->getNNZs(); z++) {
    overflow_acc(card, compute_r(kdn, A->d(z)));
  }

  // Cleanup
  unpacked_node::recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeL(card);
  CT0->addEntry(CTsrch, CTresult[0]);

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
    res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode());
  }
  long compute_r(int k, node_handle a);
};

long MEDDLY::card_mxd_int::compute_r(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0;
  if (0==k) return 1;

  // Quickly deal with skipped levels
  if (isLevelAbove(k, argF->getNodeLevel(a))) {
    if (k<0 && argF->isIdentityReduced()) {
      // identity node
      return compute_r(argF->downLevel(k), a);
    }
    // redundant node
    return overflow_mult(compute_r(argF->downLevel(k), a), argF->getLevelSize(k));
  }
  
  // Check compute table
  compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a); 
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    CT0->recycle(CTsrch);
    return CTresult[0].readL();
  }

  // Initialize node reader
  unpacked_node* A = unpacked_node::newFromNode(argF, a, false);

  // Recurse
  long card = 0;
  int kdn = argF->downLevel(k);
  for (unsigned z=0; z<A->getNNZs(); z++) {
    overflow_acc(card, compute_r(kdn, A->d(z)));
  }
  
  // Cleanup
  unpacked_node::recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeL(card);
  CT0->addEntry(CTsrch, CTresult[0]);

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
protected:
    inline compute_table::entry_key*
    findResult(node_handle a, double &b)
    {
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      b = CTresult[0].readI();
      CT0->recycle(CTsrch);
      return 0;
    }
};

MEDDLY::card_real::card_real(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, 1, arg, REAL)
{
  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:D");
  et->setForestForSlot(0, arg);
  registerEntryType(0, et);
  buildCTs();
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
      // if(lastNode>0)
      // {
      // belowcount=new long[lastNode];//argF->getCurrentNumNodes()];
      // for(int i=0;i<lastNode;i++)
      // belowcount[i]=0;
      // resultmap.clear();
      // }

    res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode());
    // FILE_output meddlyout(stdout);
    //CT0->show(meddlyout,2);
    // printf("lastNode in CARD %d\n",lastNode );
    // delete belowcount;

//     for(int i=1;i<int(res)+1;i++){
//     compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
//     MEDDLY_DCASSERT(CTsrch);
//     CTsrch->writeN(i);
// //     printf("getCurrentNumNodes: %ld\n",argF->getCurrentNumNodes() );
// // printf("getKeyType: %d\n", etype[0]->getKeyType(0));
//     CT0->find(CTsrch, CTresult[0]);
//     if (CTresult[0]) {
//
//         // int bres=
//          printf("Card result for %d is %d \n",i ,(int)CTresult[0].readD());
//
//      belowcount[i-1]=(int)CTresult[0].readD();
//     }
// }
    // belowcount=new int[(int)res];



    // if(lastNode>0){
    // std::map<int, long>::iterator itr;
    // for (itr = resultmap.begin(); itr != resultmap.end(); ++itr) {
    //     belowcount[itr->first]=itr->second;
    // }
    // resultmap.clear();
    // }


    // delete belowcount;
  }
  double compute_r(int ht, node_handle a);
// private:
//     std::map<int, long> resultmap;
};

double MEDDLY::card_mdd_real::compute_r(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0.0;
  if (0==k) return 1.0;

  // Quickly deal with skipped levels
  if (argF->getNodeLevel(a) < k) {
    return compute_r(k-1, a) * argF->getLevelSize(k);
  }
  
  // Check compute table
  compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a); 
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    CT0->recycle(CTsrch);
    return CTresult[0].readD();
  }

  // Initialize node reader
  unpacked_node* A = unpacked_node::newFromNode(argF, a, false);

  // Recurse
  double card = 0;
  int kdn = k-1;
  for (unsigned z=0; z<A->getNNZs(); z++) {
    card += compute_r(kdn, A->d(z));
  }

  // Cleanup
  unpacked_node::recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeD(card);
  // resultmap[a]=card;
  CT0->addEntry(CTsrch, CTresult[0]);

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
    res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode());
  }
  double compute_r(int k, node_handle a);
};

double MEDDLY::card_mxd_real::compute_r(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0.0;
  if (0==k) return 1.0;

  // Quickly deal with skipped levels
  if (isLevelAbove(k, argF->getNodeLevel(a))) {
    if (k<0 && argF->isIdentityReduced()) {
      // identity node
      return compute_r(argF->downLevel(k), a);
    }
    // redundant node
    return compute_r(argF->downLevel(k), a) * argF->getLevelSize(k);
  }
  
  // Check compute table
  compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a); 
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    CT0->recycle(CTsrch);
    return CTresult[0].readD();
  }

  // Initialize node reader
  unpacked_node* A = unpacked_node::newFromNode(argF, a, false);

  // Recurse
  double card = 0;
  int kdn = argF->downLevel(k);
  for (unsigned z=0; z<A->getNNZs(); z++) {
    card += compute_r(kdn, A->d(z));
  }
  
  // Cleanup
  unpacked_node::recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeD(card);
  CT0->addEntry(CTsrch, CTresult[0]);


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
};

MEDDLY::card_mpz::card_mpz(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, 1, arg, HUGEINT)
{
  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:G");
  et->setForestForSlot(0, arg);
  registerEntryType(0, et);
  buildCTs();
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
    compute_r(argF->getDomain()->getNumVariables(), a.getNode(), mcard);
  }
  void compute_r(int k, node_handle a, mpz_object &b);
};

void MEDDLY::card_mdd_mpz::compute_r(int k, node_handle a, mpz_object &card)
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
    compute_r(k-1, a, card);
    card.multiply(argF->getLevelSize(k));
    return;
  }
  
  // Check compute table
  compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    ct_object* G = CTresult[0].readG();
    mpz_object* answer = smart_cast <mpz_object*> (G);
    MEDDLY_DCASSERT(answer);
    answer->copyInto(card);
    CT0->recycle(CTsrch);
    return;
  }

  // Initialize node reader
  unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
  MEDDLY_DCASSERT(!A->isExtensible());

  // Recurse
  mpz_object tmp;
  tmp.setValue(0);
  card.setValue(0);
  int kdn = k-1;
  for (unsigned z=0; z<A->getNNZs(); z++) {
    compute_r(kdn, A->d(z), tmp);
    card.add(tmp);
  }
  
  // Cleanup
  unpacked_node::recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeG(new mpz_object(card));
  CT0->addEntry(CTsrch, CTresult[0]);

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
    compute_r(argF->getDomain()->getNumVariables(), a.getNode(), mcard);
  }
  void compute_r(int k, node_handle a, mpz_object &b);
};

void MEDDLY::card_mxd_mpz::compute_r(int k, node_handle a, mpz_object &card)
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
      compute_r(argF->downLevel(k), a, card);
      return;
    }
    // redundant node
    compute_r(argF->downLevel(k), a, card);
    card.multiply(argF->getLevelSize(k));
    return;
  }
  
  // Check compute table
  compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    ct_object* G = CTresult[0].readG();
    mpz_object* answer = smart_cast <mpz_object*> (G);
    MEDDLY_DCASSERT(answer);
    answer->copyInto(card);
    CT0->recycle(CTsrch);
    return;
  }

  // Initialize node reader
  unpacked_node* A = unpacked_node::newFromNode(argF, a, false);

  // Recurse
  mpz_object tmp;
  tmp.setValue(0);
  card.setValue(0);
  int kdn = argF->downLevel(k);
  for (unsigned z=0; z<A->getNNZs(); z++) {
    compute_r(kdn, A->d(z), tmp);
    card.add(tmp);
  }
  
  // Cleanup
  unpacked_node::recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeG(new mpz_object(card));
  CT0->addEntry(CTsrch, CTresult[0]);

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
      throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
  }
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeCardinality()
{
  return new card_opname;
}


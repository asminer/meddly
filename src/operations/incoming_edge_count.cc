
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
#include "incoming_edge_count.h"
#endif
#ifdef HAVE_LIBGMP
#include <gmp.h>
#endif
#include "../defines.h"
#include <set>
//#include "mpz_object.h"

// #define DEBUG_IEC

namespace MEDDLY {
  class iec_int;
  class iec_mdd_int;
 // class card_mxd_int;
  int* incomingedgecount;
  class iec_real;
  class iec_mdd_real;
  /*  class card_mxd_real;

#ifdef HAVE_LIBGMP
  class card_mpz;
  class card_mdd_mpz;
  class card_mxd_mpz;
#endif
 */
  class iec_opname;

};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        IEC  operations                        *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                         iec_int class                         *
// *                                                                *
// ******************************************************************

//  Abstract base class: cardinality that returns an integer
class MEDDLY::iec_int : public unary_operation {
public:
  iec_int(const unary_opname* oc, expert_forest* arg);
  std::set<node_handle> visitedNode;
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
MEDDLY::iec_int::iec_int(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, 1, arg, INTEGER)
{
  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:L");
  et->setForestForSlot(0, arg);
  registerEntryType(0, et);
  buildCTs();
  visitedNode.clear();
}

// ******************************************************************
// *                                                                *
// *                       iec_mdd_int class                       *
// *                                                                *
// ******************************************************************

//  incoming edge count on MDDs, returning integer
class MEDDLY::iec_mdd_int : public iec_int {
public:
  iec_mdd_int(const unary_opname* oc, expert_forest* arg)
    : iec_int(oc, arg) { }
  virtual void compute(const dd_edge &arg, long &res) {
    res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode());
  }
  long compute_r(int k, node_handle a);
protected:
  inline compute_table::entry_key*
  findResult(node_handle a, long &b)
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
  inline long saveResult(compute_table::entry_key* Key,
    node_handle a, long &b)
  {
    CTresult[0].reset();
    CTresult[0].writeI(b);
    CT0->addEntry(Key, CTresult[0]);
    return b;
  }
};

long MEDDLY::iec_mdd_int::compute_r(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0;
  if (0==k) return 0;
  // Quickly deal with skipped levels
  if (argF->getNodeLevel(a) < k) {
    return overflow_mult(compute_r(k-1, a), argF->getLevelSize(k));
  }

  // Check compute table
  long iec = 0;
  if(visitedNode.find(a)==visitedNode.end()){

    // Initialize node reader
    unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
    int kdn = k-1;
    for (unsigned z=0; z<A->getNNZs(); z++) {
      if(kdn>0){
        long iec;
        compute_table::entry_key* CTsrch = findResult(A->d(z), iec);
        if (0==CTsrch) iec++;
        else iec=1;
        saveResult(CTsrch, A->d(z), iec);
        overflow_acc(iec, compute_r(kdn, A->d(z)));
      }
    }
    visitedNode.insert(a);

  // compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  // MEDDLY_DCASSERT(CTsrch);
  // CTsrch->writeN(a);
  // CT0->find(CTsrch, CTresult[0]);
  // if (CTresult[0]) {
  //   CT0->recycle(CTsrch);
  //   iec=CTresult[0].readL();
  //   overflow_acc(iec,1);
  //   return iec;
  // }
  // else{
	//   iec=1;
  // }



  // Recurse

  // int kdn = k-1;
  // for (unsigned z=0; z<A->getNNZs(); z++) {
	// if(kdn>0)
  //   overflow_acc(iec, compute_r(kdn, A->d(z)));
  // }

  // Cleanup
  unpacked_node::recycle(A);

  // Add entry to compute table
  // CTresult[0].reset();
  // CTresult[0].writeL(iec);
  // CT0->addEntry(CTsrch, CTresult[0]);

#ifdef DEBUG_CARD
  fprintf(stderr, "iec of node %d is %ld(L)\n", a, iec);
#endif
  return iec;
  }
}


// // ******************************************************************
// // *                                                                *
// // *                       card_mxd_int class                       *
// // *                                                                *
// // ******************************************************************
//
// //  Cardinality on MxDs, returning integer
// class MEDDLY::card_mxd_int : public card_int {
// public:
//   card_mxd_int(const unary_opname* oc, expert_forest* arg)
//     : card_int(oc, arg) { }
//   virtual void compute(const dd_edge &arg, long &res) {
//     res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode());
//   }
//   long compute_r(int k, node_handle a);
// };
//
// long MEDDLY::card_mxd_int::compute_r(int k, node_handle a)
// {
//   // Terminal cases
//   if (0==a) return 0;
//   if (0==k) return 1;
//
//   // Quickly deal with skipped levels
//   if (isLevelAbove(k, argF->getNodeLevel(a))) {
//     if (k<0 && argF->isIdentityReduced()) {
//       // identity node
//       return compute_r(argF->downLevel(k), a);
//     }
//     // redundant node
//     return overflow_mult(compute_r(argF->downLevel(k), a), argF->getLevelSize(k));
//   }
//
//   // Check compute table
//   compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
//   MEDDLY_DCASSERT(CTsrch);
//   CTsrch->writeN(a);
//   CT0->find(CTsrch, CTresult[0]);
//   if (CTresult[0]) {
//     CT0->recycle(CTsrch);
//     return CTresult[0].readL();
//   }
//
//   // Initialize node reader
//   unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
//
//   // Recurse
//   long card = 0;
//   int kdn = argF->downLevel(k);
//   for (unsigned z=0; z<A->getNNZs(); z++) {
//     overflow_acc(card, compute_r(kdn, A->d(z)));
//   }
//
//   // Cleanup
//   unpacked_node::recycle(A);
//
//   // Add entry to compute table
//   CTresult[0].reset();
//   CTresult[0].writeL(card);
//   CT0->addEntry(CTsrch, CTresult[0]);
//
// #ifdef DEBUG_CARD
//   fprintf(stderr, "Cardinality of node %d is %ld(L)\n", a, card);
// #endif
//   return card;
// }
//
//
 // ******************************************************************
 // *                                                                *
 // *                        iec_real  class                        *
 // *                                                                *
 // ******************************************************************

 //  Abstract base class: cardinality that returns a real
 class MEDDLY::iec_real : public unary_operation {
 public:
   iec_real(const unary_opname* oc, expert_forest* arg);
protected:
    bool* visitedNode;
    //int* iec;
 };

 MEDDLY::iec_real::iec_real(const unary_opname* oc, expert_forest* arg)
  : unary_operation(oc, 1, arg, REAL)
 {
   compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:D");
   et->setForestForSlot(0, arg);
   registerEntryType(0, et);
   buildCTs();
 }

 // ******************************************************************
 // *                                                                *
 // *                      iec_mdd_real  class                      *
 // *                                                                *
 // ******************************************************************

 //  Cardinality on MDDs, returning real
 class MEDDLY::iec_mdd_real : public iec_real {
 public:
   iec_mdd_real(const unary_opname* oc, expert_forest* arg)
     : iec_real(oc, arg) { }
   virtual void compute(const dd_edge &arg, double &res) {
	  incomingedgecount=new int[(int)arg.getCardinality()+1];
     visitedNode=new bool[(int)arg.getCardinality()+1];
     res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode());
#ifdef DEBUG_IEC
     for(int i=0;i<arg.getCardinality()+1;i++)
     {
    	 printf("i %d %d \n",i, incomingedgecount[i]);
     }
#endif
   }
   double compute_r(int ht, node_handle a);
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
  inline long saveResult(compute_table::entry_key* Key,
    node_handle a, double &b)
  {
    CTresult[0].reset();
    CTresult[0].writeI(b);
    CT0->addEntry(Key, CTresult[0]);
    return b;
  }
 };

 double MEDDLY::iec_mdd_real::compute_r(int k, node_handle a)
 {
   // Terminal cases
   if (0==a) return 0.0;
   if (0==k) return 0.0;

   // Quickly deal with skipped levels
   if (argF->getNodeLevel(a) < k) {
     return compute_r(k-1, a) * argF->getLevelSize(k);
   }
   if(!visitedNode[a]){
		unpacked_node* A = unpacked_node::newFromNode(argF, a, false);

		// Recurse

		int kdn = k - 1;
		for (unsigned z = 0; z < A->getNNZs(); z++) {
			if(kdn>0){
				incomingedgecount[int(A->d(z))]++;
				compute_r(kdn, A->d(z));
			}
		}
		visitedNode[a]=true;
   }



//   // Check compute table
//   compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
//   MEDDLY_DCASSERT(CTsrch);
//   CTsrch->writeN(a);
//   CT0->find(CTsrch, CTresult[0]);
//   if (CTresult[0]) {
//     CT0->recycle(CTsrch);
//     return CTresult[0].readD();
//   }
//
//   // Initialize node reader
//   unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
//
//   // Recurse
//   double card = 0;
//   int kdn = k-1;
//   for (unsigned z=0; z<A->getNNZs(); z++) {
//     card += compute_r(kdn, A->d(z));
//   }
//
//   // Cleanup
//   unpacked_node::recycle(A);
//
//   // Add entry to compute table
//   CTresult[0].reset();
//   CTresult[0].writeD(card);
//   CT0->addEntry(CTsrch, CTresult[0]);

 #ifdef DEBUG_IEC
   fprintf(stderr, "IEC of node %d is %le(L)\n", a, iec[a]);
 #endif
   return incomingedgecount[a];
 }



// // ******************************************************************
// // *                                                                *
// // *                      card_mxd_real  class                      *
// // *                                                                *
// // ******************************************************************
//
// //  Cardinality on MxDs, returning real
// class MEDDLY::card_mxd_real : public card_real {
// public:
//   card_mxd_real(const unary_opname* oc, expert_forest* arg)
//     : card_real(oc, arg) { }
//   virtual void compute(const dd_edge &arg, double &res) {
//     res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode());
//   }
//   double compute_r(int k, node_handle a);
// };
//
// double MEDDLY::card_mxd_real::compute_r(int k, node_handle a)
// {
//   // Terminal cases
//   if (0==a) return 0.0;
//   if (0==k) return 1.0;
//
//   // Quickly deal with skipped levels
//   if (isLevelAbove(k, argF->getNodeLevel(a))) {
//     if (k<0 && argF->isIdentityReduced()) {
//       // identity node
//       return compute_r(argF->downLevel(k), a);
//     }
//     // redundant node
//     return compute_r(argF->downLevel(k), a) * argF->getLevelSize(k);
//   }
//
//   // Check compute table
//   compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
//   MEDDLY_DCASSERT(CTsrch);
//   CTsrch->writeN(a);
//   CT0->find(CTsrch, CTresult[0]);
//   if (CTresult[0]) {
//     CT0->recycle(CTsrch);
//     return CTresult[0].readD();
//   }
//
//   // Initialize node reader
//   unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
//
//   // Recurse
//   double card = 0;
//   int kdn = argF->downLevel(k);
//   for (unsigned z=0; z<A->getNNZs(); z++) {
//     card += compute_r(kdn, A->d(z));
//   }
//
//   // Cleanup
//   unpacked_node::recycle(A);
//
//   // Add entry to compute table
//   CTresult[0].reset();
//   CTresult[0].writeD(card);
//   CT0->addEntry(CTsrch, CTresult[0]);
//
//
// #ifdef DEBUG_CARD
//   fprintf(stderr, "Cardinality of node %d is %le\n", a, card);
// #endif
//   return card;
// }
//
//
//
//
// // ******************************************************************
// // *                                                                *
// // *                         card_mpz class                         *
// // *                                                                *
// // ******************************************************************
//
// #ifdef HAVE_LIBGMP
//
// //  Abstract base class: cardinality that returns large (mpz) integers.
// class MEDDLY::card_mpz : public unary_operation {
// public:
//   card_mpz(const unary_opname* oc, expert_forest* arg);
// };
//
// MEDDLY::card_mpz::card_mpz(const unary_opname* oc, expert_forest* arg)
//  : unary_operation(oc, 1, arg, HUGEINT)
// {
//   compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:G");
//   et->setForestForSlot(0, arg);
//   registerEntryType(0, et);
//   buildCTs();
// }
//
// #endif
//
// // ******************************************************************
// // *                                                                *
// // *                       card_mdd_mpz class                       *
// // *                                                                *
// // ******************************************************************
//
// #ifdef HAVE_LIBGMP
//
// /// Cardinality of MDDs, returning large (mpz) integers.
// class MEDDLY::card_mdd_mpz : public card_mpz {
// public:
//   card_mdd_mpz(const unary_opname* oc, expert_forest* arg)
//     : card_mpz(oc, arg) { }
//   virtual void compute(const dd_edge& a, ct_object &res) {
//     mpz_object& mcard = dynamic_cast <mpz_object &> (res);
//     compute_r(argF->getDomain()->getNumVariables(), a.getNode(), mcard);
//   }
//   void compute_r(int k, node_handle a, mpz_object &b);
// };
//
// void MEDDLY::card_mdd_mpz::compute_r(int k, node_handle a, mpz_object &card)
// {
//   // Terminal cases
//   if (0==a) {
//     card.setValue(0);
//     return;
//   }
//   if (0==k) {
//     card.setValue(1);
//     return;
//   }
//
//   // Quickly deal with skipped levels
//   if (argF->getNodeLevel(a) < k) {
//     // skipped level
//     compute_r(k-1, a, card);
//     card.multiply(argF->getLevelSize(k));
//     return;
//   }
//
//   // Check compute table
//   compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
//   MEDDLY_DCASSERT(CTsrch);
//   CTsrch->writeN(a);
//   CT0->find(CTsrch, CTresult[0]);
//   if (CTresult[0]) {
//     ct_object* G = CTresult[0].readG();
//     mpz_object* answer = smart_cast <mpz_object*> (G);
//     MEDDLY_DCASSERT(answer);
//     answer->copyInto(card);
//     CT0->recycle(CTsrch);
//     return;
//   }
//
//   // Initialize node reader
//   unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
//   MEDDLY_DCASSERT(!A->isExtensible());
//
//   // Recurse
//   mpz_object tmp;
//   tmp.setValue(0);
//   card.setValue(0);
//   int kdn = k-1;
//   for (unsigned z=0; z<A->getNNZs(); z++) {
//     compute_r(kdn, A->d(z), tmp);
//     card.add(tmp);
//   }
//
//   // Cleanup
//   unpacked_node::recycle(A);
//
//   // Add entry to compute table
//   CTresult[0].reset();
//   CTresult[0].writeG(new mpz_object(card));
//   CT0->addEntry(CTsrch, CTresult[0]);
//
// #ifdef DEBUG_CARD
//   fprintf(stderr, "Cardinality of node %d is ", a);
//   card.show(stderr);
//   fprintf(stderr, "\n");
// #endif
// }
//
// #endif
//
// // ******************************************************************
// // *                                                                *
// // *                       card_mxd_mpz class                       *
// // *                                                                *
// // ******************************************************************
//
// #ifdef HAVE_LIBGMP
//
// /// Cardinality of MxDs, returning large (mpz) integers.
// class MEDDLY::card_mxd_mpz : public card_mpz {
// public:
//   card_mxd_mpz(const unary_opname* oc, expert_forest* arg)
//     : card_mpz(oc, arg) { }
//   virtual void compute(const dd_edge& a, ct_object &res) {
//     mpz_object& mcard = dynamic_cast <mpz_object &> (res);
//     compute_r(argF->getDomain()->getNumVariables(), a.getNode(), mcard);
//   }
//   void compute_r(int k, node_handle a, mpz_object &b);
// };
//
// void MEDDLY::card_mxd_mpz::compute_r(int k, node_handle a, mpz_object &card)
// {
//   // Terminal cases
//   if (0==a) {
//     card.setValue(0);
//     return;
//   }
//   if (0==k) {
//     card.setValue(1);
//     return;
//   }
//
//   // Quickly deal with skipped levels
//   if (isLevelAbove(k, argF->getNodeLevel(a))) {
//     if (k<0 && argF->isIdentityReduced()) {
//       // identity node
//       compute_r(argF->downLevel(k), a, card);
//       return;
//     }
//     // redundant node
//     compute_r(argF->downLevel(k), a, card);
//     card.multiply(argF->getLevelSize(k));
//     return;
//   }
//
//   // Check compute table
//   compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
//   MEDDLY_DCASSERT(CTsrch);
//   CTsrch->writeN(a);
//   CT0->find(CTsrch, CTresult[0]);
//   if (CTresult[0]) {
//     ct_object* G = CTresult[0].readG();
//     mpz_object* answer = smart_cast <mpz_object*> (G);
//     MEDDLY_DCASSERT(answer);
//     answer->copyInto(card);
//     CT0->recycle(CTsrch);
//     return;
//   }
//
//   // Initialize node reader
//   unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
//
//   // Recurse
//   mpz_object tmp;
//   tmp.setValue(0);
//   card.setValue(0);
//   int kdn = argF->downLevel(k);
//   for (unsigned z=0; z<A->getNNZs(); z++) {
//     compute_r(kdn, A->d(z), tmp);
//     card.add(tmp);
//   }
//
//   // Cleanup
//   unpacked_node::recycle(A);
//
//   // Add entry to compute table
//   CTresult[0].reset();
//   CTresult[0].writeG(new mpz_object(card));
//   CT0->addEntry(CTsrch, CTresult[0]);
//
// #ifdef DEBUG_CARD
//   fprintf(stderr, "Cardinality of node %d is ", a);
//   card.show(stderr);
//   fprintf(stderr, "\n");
// #endif
// }
//
// #endif



// ******************************************************************
// *                                                                *
// *                       card_opname  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::iec_opname : public unary_opname {
  public:
    iec_opname();
    virtual unary_operation*
      buildOperation(expert_forest* ar, opnd_type res) const;
};

MEDDLY::iec_opname::iec_opname()
 : unary_opname("IEC")
{
}

MEDDLY::unary_operation*
MEDDLY::iec_opname::buildOperation(expert_forest* arg, opnd_type res) const
{
  if (0==arg) return 0;
  switch (res) {

    case INTEGER:

      if (arg->isForRelations())
    	  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        //return new card_mxd_int(this, arg);
        //return new iec_mdd_int(this, arg);

      else
        return new iec_mdd_int(this, arg);

     case REAL:
       if (arg->isForRelations());
         //return new iec_mxd_real(this, arg);
       else
         return new iec_mdd_real(this, arg);

// #ifdef HAVE_LIBGMP
//     case HUGEINT:
//       if (arg->isForRelations())
//         return new card_mxd_mpz(this, arg);
//       else
//         return new card_mdd_mpz(this, arg);
// #endif

    default:
      throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
  }
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeIncomingEdgeCount()
{
  return new iec_opname;
}
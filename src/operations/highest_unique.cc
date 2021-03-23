
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
#include "highest_unique.h"
#endif
#ifdef HAVE_LIBGMP
#include <gmp.h>
#endif
#include "../defines.h"
#include <set>
#include <map>
//#include "mpz_object.h"

// #define DEBUG_AC

namespace MEDDLY {
  class hu_int;
  class hu_mdd_int;
 // class card_mxd_int;
 std::map< int, std::set<int>> highestunique;
  //int* abovecount;
  class hu_real;
  class hu_mdd_real;
  /*  class card_mxd_real;

#ifdef HAVE_LIBGMP
  class card_mpz;
  class card_mdd_mpz;
  class card_mxd_mpz;
#endif
 */
  class hu_opname;

};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        AC  operations                          *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                         ac_int class                         *
// *                                                                *
// ******************************************************************

//  Abstract base class: cardinality that returns an integer
class MEDDLY::hu_int : public unary_operation {
public:
  hu_int(const unary_opname* oc, expert_forest* arg);
  //std::set<node_handle> visitedNode;
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
MEDDLY::hu_int::hu_int(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, 1, arg, INTEGER)
{
  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:L");
  et->setForestForSlot(0, arg);
  registerEntryType(0, et);
  buildCTs();
  // visitedNode.clear();
}

// ******************************************************************
// *                                                                *
// *                       ac_mdd_int class                       *
// *                                                                *
// ******************************************************************

//  incoming edge count on MDDs, returning integer
class MEDDLY::hu_mdd_int : public hu_int {
public:
  hu_mdd_int(const unary_opname* oc, expert_forest* arg)
    : hu_int(oc, arg) { }
  virtual void compute(const dd_edge &arg, long &res) {

    /*  std::set<int> a={10};
      highestunique.insert ( std::pair<int,std::set<int>>(1,a) );
      std::set<int> rset=highestunique[1];
      printf("HAHA\n" );
      for (auto it=rset.begin(); it != rset.end(); ++it)
        printf("%d\n", *it);*/

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

long MEDDLY::hu_mdd_int::compute_r(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0;
  if (0==k) return 0;
  // Quickly deal with skipped levels
  if (argF->getNodeLevel(a) < k) {
    return overflow_mult(compute_r(k-1, a), argF->getLevelSize(k));
  }
  return 0;

  // Check compute table
  // long ac = 0;
  // if(visitedNode.find(a)==visitedNode.end()){
  //
  //   // Initialize node reader
  //   unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
  //   int kdn = k-1;
  //   for (unsigned z=0; z<A->getNNZs(); z++) {
  //     if(kdn>0){
  //       long ac;
  //       compute_table::entry_key* CTsrch = findResult(A->d(z), ac);
  //       if (0==CTsrch) ac++;
  //       else ac=1;
  //       saveResult(CTsrch, A->d(z), ac);
  //       overflow_acc(ac, compute_r(kdn, A->d(z)));
  //     }
  //   }
  //   visitedNode.insert(a);

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
  // unpacked_node::recycle(A);

  // Add entry to compute table
  // CTresult[0].reset();
  // CTresult[0].writeL(iec);
  // CT0->addEntry(CTsrch, CTresult[0]);

// #ifdef DEBUG_CARD
//   fprintf(stderr, "highest unique of node %d is %ld(L)\n", a, ac);
// #endif
//   return ac;
  // }
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
 // *                        ac_real  class                        *
 // *                                                                *
 // ******************************************************************

 //  Abstract base class: cardinality that returns a real
 class MEDDLY::hu_real : public unary_operation {
 public:
   hu_real(const unary_opname* oc, expert_forest* arg);
//protected:
    //int* iec;
 };

 MEDDLY::hu_real::hu_real(const unary_opname* oc, expert_forest* arg)
  : unary_operation(oc, 1, arg, REAL)
 {
   compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:D");
   et->setForestForSlot(0, arg);
   registerEntryType(0, et);
   buildCTs();
 }

 // ******************************************************************
 // *                                                                *
 // *                      ac_mdd_real  class                      *
 // *                                                                *
 // ******************************************************************

 //  Cardinality on MDDs, returning real
 class MEDDLY::hu_mdd_real : public hu_real {
 public:
   hu_mdd_real(const unary_opname* oc, expert_forest* arg)
     : hu_real(oc, arg) { }
   virtual void compute(const dd_edge &arg, double &res) {
       // for(int i=0;i<arg.getCardinality()+1;i++)
       // {
      	//  printf("beforeAC %d %d \n",i, incomingedgecount[i]);
       // }
// 	  abovecount=new int[(int)arg.getCardinality()+1];
//      iec=new int[(int)arg.getCardinality()+1];
//      for(int i=0;i<arg.getCardinality()+1;i++)
//      {
//     	 iec[i]=incomingedgecount[i];
//      }


int card = lastNode;//argF->getCurrentNumNodes();// instead of arg.getCardinality()
if(card<1)
{
    throw error(error::INVALID_SEQUENCE, __FILE__, __LINE__);
}

// std::set<int> a={10};
// highestunique.insert ( std::pair<int,std::set<int>>(1,a) );
// std::set<int> rset=highestunique[1];
// printf("HAHA\n" );
// updateInsert(1,10);
// for (auto it=rset.begin(); it != rset.end(); ++it)
//   printf("%d\n", *it);
// rset.insert(20);
// highestunique.at(1)=rset;
// rset=highestunique[1];
// printf("NEW\n" );
// updateInsert(1,20);
// for (auto it=rset.begin(); it != rset.end(); ++it)
//   printf("%d\n", *it);
// int arrsize=sizeof(incomingedgecount)/sizeof(incomingedgecount[0]);
// printf("ARRSIZE %d XXXX\n",arrsize );
// int card = arrsize-1;
    explored=new bool[card];
    incomingedgecountHU=new int[card];
    lastPosition=new int[card];
    startInterval= new int[card];
    stopInterval=new int[card];
    if(explored==0 || incomingedgecountHU==0 || lastPosition==0||startInterval==0|| stopInterval==0)
    {printf("ERROR IN making arrays\n" ); char c=getchar();}
    for(int i=0;i<card;i++)
    {
        explored[i]=false;
        incomingedgecountHU[i]=0;
        lastPosition[i]=0;
        startInterval[i]=0;
        stopInterval[i]=0;
    }
    for(int i=0;i<card;i++){
    incomingedgecountHU[i]=incomingedgecount[i];
    }
    pset = new std::set<int>[card];
     // res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode(),pset );
     res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode() );

     // printf("COMPUTED %d\n",card );
     // printf("***HU***\n" );
     // for(int i=0;i<=(int)card;i++){
   	 //  printf("HU[%d]= \t", i);
     //    std::set<int> rset=highestunique[i];
     //    for (auto it=rset.begin(); it != rset.end(); ++it){
     //       printf(", %d", *it);
     //
     //   }
     //   printf("[%d, \t %d ,\t %d ]\t \n",startInterval[i],stopInterval[i],lastPosition[i] );
     //     printf("\n" );
     //  }
     //  printf("**HU****\n" );
     for(int i=0;i<card;i++){
     pset[i].clear();
     }
     delete[] pset;//=NULL;
     delete[] incomingedgecountHU;
     delete[] explored;
     delete[] startInterval;
     delete[] stopInterval;
     delete[] lastPosition;
     //  std::set<int> rset=highestunique[2];
     //  for (auto it=rset.begin(); it != rset.end(); ++it){
     //     printf(", %d", *it);
     //     printf("\n" );
     // }
#ifdef DEBUG_AC
     for(int i=0;i<card;i++)
     {
    	 printf("HU %d %d \n",i, highestunique[i]);
     }
#endif
   }
   // double compute_r(int ht, node_handle a,std::set<int>);
   double compute_r(int ht, node_handle a);

   // std::set<int> CheckInterval(std::set<int> pset,node_handle a );
   std::set<int> CheckInterval(node_handle a );

 protected:
     bool* explored;
     // std::set<int> pset;
     std::set<int>* pset;

     int* incomingedgecountHU;
     int* lastPosition;
     int* startInterval;
     int* stopInterval;
     int c=1;
  inline void updateInsert(node_handle a, node_handle b){
      //printf("COMING TO updateInsert\n");
      if ( highestunique.find(a) == highestunique.end() ) {
        // not found
        std::set<int> rset={b};
        highestunique.insert ( std::pair<int,std::set<int>>(a, rset) );
        // for (auto it=rset.begin(); it != rset.end(); ++it)
        //   printf("%d\n", *it);
      } else {
        // found
      std::set<int> rset=highestunique[a];
      rset.insert(b);
      highestunique.at(a)=rset;
      // for (auto it=rset.begin(); it != rset.end(); ++it)
      //   printf("%d\n", *it);
  }

  }
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
 // std::set<int> MEDDLY::hu_mdd_real::CheckInterval(std::set<int> pset,node_handle a ){
 std::set<int> MEDDLY::hu_mdd_real::CheckInterval(node_handle a ){

std::set<int> shouldBeRemoved;
// printf("COMING CheckInterval %d \n",a );
for (auto l: pset[a] ){
    // printf("L %d \t %d \t %d \t %d\n",l,  startInterval[l],stopInterval[l],lastPosition[l]);
    // printf("P %d \t %d \t %d \t %d \n",a,  startInterval[a],stopInterval[a],lastPosition[a]);

    if((startInterval[l]>=startInterval[a])&&(stopInterval[l]<=stopInterval[a])
    &&(((lastPosition[l]!=0)&&(lastPosition[l]>=startInterval[a])&&(lastPosition[l]<=stopInterval[a]))||(lastPosition[l]==0))){
        // printf("TRUE\n" );
        updateInsert(a,l);
        shouldBeRemoved.insert(l);
    }

}
// printf("END CheckInterval %d\n",a );

return shouldBeRemoved;


 }

 // double MEDDLY::hu_mdd_real::compute_r(int k, node_handle a,std::set<int>pset)
 double MEDDLY::hu_mdd_real::compute_r(int k, node_handle a)

 {
     // printf("COMING TO hu_mdd_real \n" );
     if(explored[a]==1){
         // printf("EXP %d is true\n",a-1 );
         // printf("else  %d\n",a-1);

         lastPosition[a]=c;
         c++;
        //  for (auto it=pset[a-1].begin(); it != pset[a-1].end(); ++it)
        // printf("else Pset is %d\n",*it );
        //  printf("Done else  %d\n",a-1);
     }else{

         // printf("C is %d\n",c );
         startInterval[a]=c;
         // printf("startInterval[ %d ] =%d\n",a-1,c );
         c++;
         int kdn = k - 1;
         unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
         for (unsigned z = 0; z < A->getNNZs(); z++) {
             if (kdn>0)
             {
                incomingedgecountHU[A->d(z)]--;
                // std::set<int> childpset;
                // compute_r(k-1,A->d(z),childpset);
                compute_r(k-1,A->d(z));

                if(incomingedgecountHU[A->d(z)]==0)
                {
                     // printf("ADDED %d to pset of %d \n",int(A->d(z))-1, a-1 );
                    // pset.insert(int(A->d(z))-1);
                    // pset.insert(childpset.begin(), childpset.end());
                    pset[a].insert(A->d(z));
                    pset[a].insert(pset[A->d(z)].begin(), pset[A->d(z)].end());
                   //  for (auto it=pset[a-1].begin(); it != pset[a-1].end(); ++it)
                   // printf("Pset is**** %d\n",*it );

                }

             }
         }
         stopInterval[a]=c;
         // printf("stopInterval[ %d ] =%d\n",a-1,c );

         c++;
         explored[a]=true;
         // printf("explored[ %d ] is true\n",a-1 );
         unpacked_node::recycle(A);
     }
     if(incomingedgecountHU[a]==0)
     {
         // printf("incomingedgecountHU %d  is ZERO\n", a-1);
         std::set<int> result;
        //  for (auto it=pset[a-1].begin(); it != pset[a-1].end(); ++it)
        // printf("Pset is %d\n",*it );
    // std::set<int> toRemove= CheckInterval(pset,a-1);
    std::set<int> toRemove= CheckInterval(a);

    for (auto r: toRemove ){
        pset[a].erase(r);
    }
          }
     return 0;
 }
   // Terminal cases
   // if (0==a) return 0.0;
   // if (0==k) return 0.0;

   // // Quickly deal with skipped levels
   // if (argF->getNodeLevel(a) < k) {
   //   return compute_r(k-1, a) * argF->getLevelSize(k);
   // }
   //
	// 	unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
   //
	// 	// Recurse
   //
	// 	int kdn = k - 1;
	// 	for (unsigned z = 0; z < A->getNNZs(); z++) {
	// 		if(kdn>0){
	// 			abovecount[int(A->d(z))]+=abovecount[a];
	// 			iec[int(A->d(z))]--;
	// 			if (iec[int(A->d(z))]==0){
	// 			compute_r(kdn, A->d(z));
	// 			}
	// 		}
	// 	}
   //
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

 // #ifdef DEBUG_AC
 //   fprintf(stderr, "HU of node %d is %le(L)\n", a, abovecount[a]);
 // #endif
 //   return abovecount[a];
 // }



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

class MEDDLY::hu_opname : public unary_opname {
  public:
    hu_opname();
    virtual unary_operation*
      buildOperation(expert_forest* ar, opnd_type res) const;
};

MEDDLY::hu_opname::hu_opname()
 : unary_opname("HU")
{
}

MEDDLY::unary_operation*
MEDDLY::hu_opname::buildOperation(expert_forest* arg, opnd_type res) const
{
  if (0==arg) return 0;
  switch (res) {

    case INTEGER:

      if (arg->isForRelations())
    	  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        //return new card_mxd_int(this, arg);
        //return new iec_mdd_int(this, arg);

      else
        return new hu_mdd_int(this, arg);

     case REAL:
       if (arg->isForRelations());
         //return new iec_mxd_real(this, arg);
       else
         return new hu_mdd_real(this, arg);

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

MEDDLY::unary_opname* MEDDLY::initializeHighestUnique()
{
  return new hu_opname;
}

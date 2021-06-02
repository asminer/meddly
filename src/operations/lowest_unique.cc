
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
#include "lowest_unique.h"
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
  class lu_int;
  class lu_mdd_int;
 // class card_mxd_int;
 std::map< int, std::set<int>> lowestunique;
  //int* abovecount;
  class lu_real;
  class lu_mdd_real;
  /*  class card_mxd_real;

#ifdef HAVE_LIBGMP
  class card_mpz;
  class card_mdd_mpz;
  class card_mxd_mpz;
#endif
 */
  class lu_opname;

};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        LU  operations                          *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                         lu_int class                         *
// *                                                                *
// ******************************************************************

//  Abstract base class: cardinality that returns an integer
class MEDDLY::lu_int : public unary_operation {
public:
  lu_int(const unary_opname* oc, expert_forest* arg);
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
MEDDLY::lu_int::lu_int(const unary_opname* oc, expert_forest* arg)
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
class MEDDLY::lu_mdd_int : public lu_int {
public:
  lu_mdd_int(const unary_opname* oc, expert_forest* arg)
    : lu_int(oc, arg) { }
  virtual void compute(const dd_edge &arg, long &res) {

    /*  std::set<int> a={10};
      lowestunique.insert ( std::pair<int,std::set<int>>(1,a) );
      std::set<int> rset=lowestunique[1];
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

long MEDDLY::lu_mdd_int::compute_r(int k, node_handle a)
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
 class MEDDLY::lu_real : public unary_operation {
 public:
   lu_real(const unary_opname* oc, expert_forest* arg);
protected:
     std::map<std::pair<int, int>, int> cache;
 };

 MEDDLY::lu_real::lu_real(const unary_opname* oc, expert_forest* arg)
  : unary_operation(oc, 1, arg, REAL)
 {
   // compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:D");
   // et->setForestForSlot(0, arg);
   // registerEntryType(0, et);
   // buildCTs();
   compute_table::entry_type* et2 = new compute_table::entry_type(oc->getName(), "NN:N");
   et2->setForestForSlot(0, arg);
   et2->setForestForSlot(1, arg);
   et2->setForestForSlot(3, resF);
   registerEntryType(0, et2);
   buildCTs();
 }

 // ******************************************************************
 // *                                                                *
 // *                      ac_mdd_real  class                      *
 // *                                                                *
 // ******************************************************************

 //  Cardinality on MDDs, returning real
 class MEDDLY::lu_mdd_real : public lu_real {
 public:
   lu_mdd_real(const unary_opname* oc, expert_forest* arg)
     : lu_real(oc, arg) { }
   virtual void compute(const dd_edge &arg, double &res) {


int card = lastNode;
if(card<1)
{
    throw error(error::INVALID_SEQUENCE, __FILE__, __LINE__);
}

    explored=new bool[card];
    LUreverse=new int[card];
    if(explored==0||LUreverse==0)
    {printf("ERROR IN making arrays\n" ); getchar();}
    for(int i=0;i<card;i++)
    {
        explored[i]=false;
        LUreverse[i]=-1;
    }
    initialize(argF->getDomain()->getNumVariables(), arg.getNode() );
    for(int i=0;i<card;i++)
    {
        explored[i]=false;
    }
    // for(int i=0;i<=(int)card;i++){
    //  printf("LUR[%d]= %d\t", i,LUreverse[i]);
    // }
    // printf("Initialized Complete\n" );
    // expert_forest* ef = static_cast<expert_forest*>(arg);
    // expert_forest* ef=(expert_forest*) arg.getForest();
    // node_handle n=arg.getNode();

    // node_handle* list = ef->markNodesInSubgraph(&n, 1, false);
    // delete ef;
    // getchar();
    // for(int i=0;i<card;i++){
    // incomingedgecountHU[i]=incomingedgecount[i];
    // }
    // pset = new std::set<int>[card];
     // res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode(),pset );
     cache.clear();
    compute_r(argF->getDomain()->getNumVariables(), arg.getNode() );
    // for(int i=0;i<=(int)card;i++){
    //  printf("LUR[%d]= %d\t", i,LUreverse[i]);
    // }
    // printf("\nComputation Complete\n" );
    for(int i=0;i<card;i++)
    {
        if(LUreverse[i]!=-1)
        {
            updateInsert(LUreverse[i],i);
        }
    }
    // for(int i=0;i<=(int)card;i++){
    //  printf("LUR[%d]= %d\t", i,LUreverse[i]);
    // }
    // getchar();
    // printf("Computation Complete after update\n" );

     // printf("COMPUTED %d\n",card );
   //   printf("***LU***\n" );
   //   for(int i=0;i<(int)card;i++){
   // 	  printf("LU[%d]= \t", i);
   //      std::set<int> rset=lowestunique[i];
   //      for (auto it=rset.begin(); it != rset.end(); ++it){
   //         printf(", %d", *it);
   //
   //     }
   //     printf("\n" );
   // }
     //   printf("[%d, \t %d ,\t %d ]\t \n",startInterval[i],stopInterval[i],lastPosition[i] );
     //     printf("\n" );
     //  }
     //  printf("**HU****\n" );
     // for(int i=0;i<card;i++){
     // pset[i].clear();
     // }
     cache.clear();
     delete[] explored;
     delete[] LUreverse;

#ifdef DEBUG_LU
     for(int i=0;i<card;i++)
     {
    	 printf("LU %d %d \n",i, lowestunique[i]);
     }
#endif
   }
   void compute_r(int ht, node_handle a);
   void initialize(int ht, node_handle a);

   // std::set<int> CheckInterval(std::set<int> pset,node_handle a );
   // std::set<int> CheckInterval(node_handle a );

 protected:
     bool* explored;
     int* LUreverse;
  inline void updateInsert(node_handle a, node_handle b){
      //printf("COMING TO updateInsert\n");
      if ( lowestunique.find(a) == lowestunique.end() ) {
        // not found
        std::set<int> rset={b};
        lowestunique.insert ( std::pair<int,std::set<int>>(a, rset) );
        // for (auto it=rset.begin(); it != rset.end(); ++it)
        //   printf("%d\n", *it);
      } else {
        // found
      std::set<int> rset=lowestunique[a];
      rset.insert(b);
      lowestunique.at(a)=rset;
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
  inline compute_table::entry_key*
  findResult(node_handle a, node_handle b, node_handle &c)
  {
    compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    MEDDLY_DCASSERT(CTsrch);

      CTsrch->writeN(a);
      CTsrch->writeN(b);

    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    c = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
  }

  inline void saveResult(compute_table::entry_key* K,
    node_handle a, node_handle b, node_handle c)
  {
    CTresult[0].reset();
    CTresult[0].writeN(c);
    CT0->addEntry(K, CTresult[0]);
  }

  inline int ChildsLeastCommonUnique(node_handle firstChild,node_handle secondChild)
  {
      // printf("ChildsLeastCommonUnique\n");

      int firstChildLVL=argF->getNodeLevel(firstChild);
      int secondChildLVL=argF->getNodeLevel(secondChild);
      // printf("firstChildLVL %d secondChildLVL %d\n", firstChildLVL,secondChildLVL);
      // printf("LU firstChildLVL %d LU secondChildLVL %d\n",LUreverse[ firstChild],LUreverse[secondChild]);
      if(firstChildLVL==0 || secondChildLVL==0) return 0;
      if(LUreverse[firstChild]==0 || LUreverse[secondChild]==0) return 0;


     /* if( LUreverse[ firstChild]==-1)
      compute_r(firstChildLVL,firstChild);
      if( LUreverse[ secondChild]==-1)
      compute_r(secondChildLVL,secondChild);
      if(LUreverse[firstChild]==0 || LUreverse[secondChild]==0) return 0;

      printf("firstChildLVL %d secondChildLVL %d\n", firstChildLVL,secondChildLVL);
      printf("LU firstChildLVL %d LU secondChildLVL %d\n",LUreverse[ firstChild],LUreverse[secondChild]);

      if(firstChildLVL==1||secondChildLVL==1) return 0;*/
      if(firstChildLVL<secondChildLVL){
          return ChildsLeastCommonUnique(firstChild,LUreverse[secondChild]);
      }
      else if(firstChildLVL>secondChildLVL){
          return ChildsLeastCommonUnique(LUreverse[firstChild],secondChild);
      }
      else if(LUreverse[firstChild]==LUreverse[secondChild]){
          return LUreverse[firstChild];
      }
      else{
          return ChildsLeastCommonUnique(LUreverse[firstChild],LUreverse[secondChild]);
      }
  }
  inline int LeastCommonUnique(int lvl ,node_handle a)
  {
       int firstChild=0;
       int secondChild=0;
       // compute_table::entry_key* Key=NULL;
       int kdn = lvl - 1;
       // if(kdn==0)
       //  return 0;
       unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
       for(unsigned z = 0; z < A->getNNZs(); z++){
           // if(kdn>0){
               if(firstChild==0){
                 firstChild=A->d(z);
                 if(LUreverse[firstChild]==0) return 0;
                 // compute_r(kdn,firstChild);
                 // printf("Set firstChild %d\n", firstChild);
               }
               else{
                   secondChild=A->d(z);
                   // printf("LU secondChild\n",LUreverse[secondChild] );
                   if(LUreverse[secondChild]==0) {
                       return 0;
                   }

                   // printf("Set secondChild %d\n", secondChild);
                    // node_handle result = 0;
                    // printf("result is set\n" );
                   //  Key = findResult(firstChild, secondChild, result);
                   //  printf("Key set\n" );
                   // if (0==Key){
                   //      unpacked_node::recycle(A);
                   //      return result;
                   //  }
                    // printf("Could not find result\n" );
                   // int oldfirstChild=firstChild;

                   //  if (LUreverse[secondChild]==-1)
                   // compute_r(kdn,secondChild);
                   // else
                   // getchar();
                   int prevFirstChild=firstChild;
                 std::map<std::pair<int,int>,int>::iterator res = cache.find(std::make_pair(prevFirstChild,secondChild));

                  if(res != cache.end())
                  {

                     // printf("FOUND %d %d %d\n",prevFirstParent, secondParent,X);
                     firstChild=res->second;

                  }
                  else{
                      // printf("Call PairLowestSeparatorAbove %d %d\n", prevFirstParent,secondParent);
                      firstChild=ChildsLeastCommonUnique(prevFirstChild,secondChild);
                 cache[std::make_pair(prevFirstChild,secondChild)]=firstChild;
                 // printf("SAVE %d %d %d\n",prevFirstParent, secondParent,firstParent);
                   }


                   // getchar();
                   // printf("Result is %d\n",firstChild );
                   // saveResult(Key, oldfirstChild, secondChild, firstChild);
                   // printf("Saved\n" );
                   // saveResult(Key, secondChild, oldfirstChild, firstChild);
                   // printf("Saved\n" );

                   if(firstChild==0) {
                        unpacked_node::recycle(A);
                       return 0;
                   }
               }
       }
        unpacked_node::recycle(A);
       return firstChild;
  }
 };
 void MEDDLY::lu_mdd_real::initialize(int k, node_handle a)

 {
     if(explored[a]==0){
         int firstChild=-1;
         bool singleChild=true;
         bool childWithLUE=false;
         unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
         int kdn = k - 1;
         if(kdn==0)
         LUreverse[a]=0;
         else if(kdn>0){
             for(unsigned z = 0; z < A->getNNZs(); z++){

                     if (firstChild==-1){
                         firstChild=A->d(z);
                     }else if(A->d(z)!=firstChild&&singleChild){
                         singleChild=false;
                     }
                     initialize(kdn,A->d(z));
                     if(LUreverse[A->d(z)]==0&&!singleChild){
                         childWithLUE=true;
                         singleChild=false;
                     }
             }
         //}
         if(singleChild/*&&kdn>0*/){
             // updateInsert(firstChild,a);
             LUreverse[a]=firstChild;
             // LUisSet[a]=true;
         }
         else if(childWithLUE/*||kdn==0*/){
             // updateInsert(1,(int)a);
             LUreverse[(int)a]=0;
             // LUisSet[a]=true;
         }
     }
         explored[a]=true;
         unpacked_node::recycle(A);
     }
 }


 void MEDDLY::lu_mdd_real::compute_r(int k, node_handle a)
 {
     // node_handle* list = this->markNodesInSubgraph(&a, 1, false);
     // getchar();
     if(!explored[a]/*&&/*a>0&& *//*LUreverse[a]==-1*/){
         // printf("a is %d and lvl is %d\n",a,k );
         unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
         int kdn = k - 1;
          if(kdn>0){
              for(unsigned z = 0; z < A->getNNZs(); z++){
                 // if(LUreverse[A->d(z)]==-1)
                 // if(kdn>0)
                 {
                     // printf("Child is %d and LU is %d\n",A->d(z),LUreverse[A->d(z)] );
                     compute_r(kdn,A->d(z));
                 }
              }

              if(LUreverse[a]==-1)
              {
                  LUreverse[a]=LeastCommonUnique(k,a);
                  // printf(" LUreverse[%d]=%d;\n",a,LUreverse[a] );
              }
          }
          else{
              LUreverse[a]=0;
          }

         unpacked_node::recycle(A);
     }
     explored[a]=true;
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

class MEDDLY::lu_opname : public unary_opname {
  public:
    lu_opname();
    virtual unary_operation*
      buildOperation(expert_forest* ar, opnd_type res) const;
};

MEDDLY::lu_opname::lu_opname()
 : unary_opname("HU")
{
}

MEDDLY::unary_operation*
MEDDLY::lu_opname::buildOperation(expert_forest* arg, opnd_type res) const
{
  if (0==arg) return 0;
  switch (res) {

    case INTEGER:

      if (arg->isForRelations())
    	  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        //return new card_mxd_int(this, arg);
        //return new iec_mdd_int(this, arg);

      else
        return new lu_mdd_int(this, arg);

     case REAL:
       if (arg->isForRelations());
         //return new iec_mxd_real(this, arg);
       else
         return new lu_mdd_real(this, arg);

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

MEDDLY::unary_opname* MEDDLY::initializeLowestUnique()
{
  return new lu_opname;
}

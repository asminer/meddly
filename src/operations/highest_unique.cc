
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
// #define HRec
 #define HTrav
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
protected:
    std::map<std::pair<int, int>, int> cache;

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
// printf("HighestUnique start\n" );
int card = lastNode;//argF->getCurrentNumNodes();// instead of arg.getCardinality()
if(card<1)
{
    throw error(error::INVALID_SEQUENCE, __FILE__, __LINE__);
}

    explored=new bool[card];
    incomingedgecountHU=new int[card];
    #ifdef HRec
    LSA=new int[card];
    firstParent= new int[card];
    setToRoot= new bool[card];
    singleParent= new bool[card];
    #endif
    #ifdef HTrav
    lastPosition=new int[card];
    startInterval= new int[card];
    stopInterval=new int[card];
    blockedBy= new int[card];
    #endif
    #ifdef HRec
    if(explored==0)
    {printf("ERROR IN making arrays\n" ); char c=getchar();}
    #endif
    #ifdef HTrav
    if(explored==0 || incomingedgecountHU==0 || lastPosition==0||startInterval==0|| stopInterval==0)
    {printf("ERROR IN making arrays\n" ); char c=getchar();}
    #endif
    for(int i=0;i<card;i++)
    {
         explored[i]=false;
        incomingedgecountHU[i]=0;
        #ifdef HRec
        firstParent[i]=-1;
        setToRoot[i]=false;
        singleParent[i]=true;
        LSA[i]=-1;
        #endif
        #ifdef HTrav
         lastPosition[i]=0;
         startInterval[i]=0;
         stopInterval[i]=0;
         blockedBy[i]=0;
         #endif
    }
    for(int i=0;i<card;i++){
    incomingedgecountHU[i]=incomingedgecount[i];
    }
    #ifdef HTrav
    pset = new std::set<int>[card];
    for(int i=0;i<card;i++){
    pset[i].clear();
    }
 // res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode(),pset );

        res = compute_rt(argF->getDomain()->getNumVariables(), arg.getNode() );
        // getchar();
        // getchar();
 // printf("COMPUTED %d\n",card );
 // printf("***HU***\n" );
 // for(int i=0;i<=(int)card;i++){
 //  printf("HU[%d]= \t", i);
 //    std::set<int> rset=highestunique[i];
 //    for (auto it=rset.begin(); it != rset.end(); ++it){
 //       printf(", %d", *it);
 //
 //   }
 //
 //   // printf("[%d, \t %d ,\t %d ]\t \n",startInterval[i],stopInterval[i],lastPosition[i] );
 //     printf("\n" );
 //  }
 //  getchar();
 //  getchar();
 //  printf("**HU****\n" );

    #endif
    #ifdef HRec
    parents = new std::set<int>[card];
    initialize(argF->getDomain()->getNumVariables(), arg.getNode(),arg.getNode());
    // printf("Initialize completed\n" );
    for(int i=0;i<card;i++){
    incomingedgecountHU[i]=incomingedgecount[i];
    }
    cache.clear();
    compute_r(argF->getDomain()->getNumVariables(), arg.getNode(),arg.getNode() );
    for(int i=0;i<card;i++){
    parents[i].clear();
    }
    for(int i=1;i<card;i++)
    {
        // if(LSA[i]==0)
        // {
        //     printf("LSA is 0 for %d\n",i );
        //     getchar();
        // }
        if(LSA[i]!=-1)
        {
            updateInsert(LSA[i],i);
        }
        else{

            if(i!=arg.getNode()&&incomingedgecount[i]>0){
            printf("ERR LSA %d is -1\n", i );
            getchar();
            }
        }
    }

    // printf("***HU-REC***\n" );
    // for(int i=0;i<=(int)card;i++){
    //  printf("HU[%d]= \t", i);
    //    std::set<int> rset=highestunique[i];
    //    for (auto it=rset.begin(); it != rset.end(); ++it){
    //       printf(", %d", *it);
    //
    //   }
    //
    //   // printf("[%d, \t %d ,\t %d ]\t \n",startInterval[i],stopInterval[i],lastPosition[i] );
    //     printf("\n" );
    //  }
    //  getchar();
    //  getchar();
    cache.clear();

    #endif

     delete[] incomingedgecountHU;
     #ifdef HRec
     delete[] parents;//=NULL;
     delete[] LSA;
     delete[] setToRoot;
     delete[] singleParent;
     delete[]firstParent;
     #endif
     #ifdef HTrav
     delete[] startInterval;
     delete[] stopInterval;
     delete[] lastPosition;
     delete[] blockedBy;
     #endif

#ifdef DEBUG_AC
     for(int i=0;i<card;i++)
     {
    	 printf("HU %d %d \n",i, highestunique[i]);
     }
#endif
// printf("HighestUnique End\n" );
   }
   /*

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


   // printf("ROOT is %d\n",arg.getNode() );
    // getchar();
     // res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode(),pset );
    // for(int i=0;i<card;i++){
    //     if(incomingedgecountHU[i]<0)
    //     printf("ERR incomingedgecountHU\n" );
    // }
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

    //  std::set<int> rset=highestunique[2];
    //  for (auto it=rset.begin(); it != rset.end(); ++it){
    //     printf(", %d", *it);
    //     printf("\n" );
    // }
    */
   // double compute_r(int ht, node_handle a,std::set<int>);
   void compute_r(int ht, node_handle a,node_handle root);
   double compute_rt(int ht, node_handle a);
   void initialize(int ht, node_handle a, node_handle root);

    // double compute_rn(int ht, node_handle a);
   // std::set<int> CheckInterval(std::set<int>* pset,node_handle a );
    std::set<int> CheckInterval(node_handle a );

 protected:
      bool* explored;



     int* incomingedgecountHU;
     #ifdef HRec
      std::set<int>* parents;
     int* LSA;
     bool* setToRoot;
     bool* singleParent;
     int* firstParent;
     #endif
     #ifdef HTrav
      std::set<int>* pset;
     int* lastPosition;
     int* startInterval;
     int* stopInterval;
     int* blockedBy;
     #endif
     int c=1;
  inline int PairLowestSeparatorAbove(node_handle firstParent,node_handle secondParent,node_handle root){
      #ifdef HRec
      int firstParentLVL=argF->getNodeLevel(firstParent);
      int secondParentLVL=argF->getNodeLevel(secondParent);
      int rootLVL=argF->getNodeLevel(root);
      // printf("firstChildLVL %d secondChildLVL %d\n", firstChildLVL,secondChildLVL);
      // printf("LU firstChildLVL %d LU secondChildLVL %d\n",LUreverse[ firstChild],LUreverse[secondChild]);
      if(firstParent==root||secondParent==root) return root;
      if(firstParentLVL==rootLVL-1 || secondParentLVL==rootLVL-1) return root;
      if(LSA[firstParent]==root || LSA[secondParent]==root) return root;
      if(LSA[firstParent]==-1 ||LSA[secondParent]==-1) {
          printf("ERR PairLowestSeparatorAbove\n" );
          printf("firstParent %d LSA %d\n",firstParent,LSA[firstParent] );
          printf("secondParent %d LSA %d\n",secondParent,LSA[secondParent] );
          getchar();
      }
      if(firstParentLVL<secondParentLVL){
          return PairLowestSeparatorAbove(LSA[firstParent],secondParent,root);
      }
      else if(firstParentLVL>secondParentLVL){
          return PairLowestSeparatorAbove(firstParent,LSA[secondParent],root);
      }
      else if(LSA[firstParent]==LSA[secondParent]){
          return LSA[firstParent];
      }
      else{
          return PairLowestSeparatorAbove(LSA[firstParent],LSA[secondParent],root);
      }
      #endif
  }
  inline int LowestSeparatorAboveSet(int lvl,std::set<int> parents,node_handle root){
      #ifdef HRec
      int firstParent=0;
      int secondParent=0;
      int kdn = lvl - 1;

      // unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
      // for(unsigned z = 0; z < A->getNNZs(); z++){
        for(auto a : parents){
            // printf("a is %d\n", a);
            // getchar();
              if(firstParent==0){
                firstParent=a;
                // printf("firstParent %d\n",a );
                if(LSA[firstParent]==root) return root;
                }
              else{
                  secondParent=a;
                  // printf("LU secondChild\n",LUreverse[secondChild] );
                  if(LSA[secondParent]==root) return root;
                  int prevFirstParent=firstParent;
                  std::map<std::pair<int,int>,int>::iterator res = cache.find(std::make_pair(prevFirstParent,secondParent));

                   if(res != cache.end())
                   {

                      // printf("FOUND %d %d %d\n",prevFirstParent, secondParent,X);
                      firstParent=res->second;

                   }
                   else{
                       // printf("Call PairLowestSeparatorAbove %d %d\n", prevFirstParent,secondParent);
                  firstParent=PairLowestSeparatorAbove(prevFirstParent,secondParent,root);
                  cache[std::make_pair(prevFirstParent,secondParent)]=firstParent;
                  // printf("SAVE %d %d %d\n",prevFirstParent, secondParent,firstParent);
                    }
                  // getchar();
                  // printf("Result is %d\n",firstChild );
                  // saveResult(Key, oldfirstChild, secondChild, firstChild);
                  // printf("Saved\n" );
                  // saveResult(Key, secondChild, oldfirstChild, firstChild);
                  // printf("Saved\n" );

                  if(firstParent==root) {
                       // unpacked_node::recycle(A);
                      return root;
                  }
              }
      }
       // unpacked_node::recycle(A);
      return firstParent;
     #endif
  }

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
 // std::set<int> MEDDLY::hu_mdd_real::CheckInterval(std::set<int>* pset,node_handle a ){
  std::set<int> MEDDLY::hu_mdd_real::CheckInterval(node_handle a ){
#ifdef HTrav
std::set<int> shouldBeRemoved;
// printf("COMING CheckInterval %d \n",a );
// if(a==1535){
//     for (auto l: pset[a] ){
//         printf("1535 pset %d\n",l );
//     }
//     getchar();
//  }
for (auto l: pset[a] ){
    // if(a==1535){
    // printf("L %d \t ST %d \t SP %d \t LP %d BB %d, INCBB %d\n",l,  startInterval[l],stopInterval[l],lastPosition[l],blockedBy[l],incomingedgecountHU[blockedBy[l]]);
    // printf("P %d \t ST %d \t SP %d \t LP %d BB %d\n",a,  startInterval[a],stopInterval[a],lastPosition[a],blockedBy[a]);
    // }
    if((blockedBy[l]==0)||(blockedBy[l]!=0 && incomingedgecountHU[blockedBy[l]]==0)||(argF->getNodeLevel(l)+1==argF->getNodeLevel(a)))
    if((startInterval[l]>=startInterval[a])&&(stopInterval[l]<=stopInterval[a])
    &&(((lastPosition[l]!=0)&&(lastPosition[l]>=startInterval[a])&&(lastPosition[l]<=stopInterval[a]))||(lastPosition[l]==0))){
        // printf("TRUE\n" );
        // if(a==1535)
        // {
        //     printf("added l %d to a %d\n",l,a );
        //     printf("L IEC %d\t %d \t %d \t %d \t %d\n",l,  incomingedgecountHU[l],startInterval[l],stopInterval[l],lastPosition[l]);
        //     printf("P IEC %d\t %d \t %d \t %d \t %d \n",a, incomingedgecountHU[a], startInterval[a],stopInterval[a],lastPosition[a]);
        //
        // }

        updateInsert(a,l);
        shouldBeRemoved.insert(l);
    }

}
// printf("END CheckInterval %d\n",a );

return shouldBeRemoved;

#endif
 }

  double MEDDLY::hu_mdd_real::compute_rt(int k, node_handle a){
      #ifdef HTrav
       // printf("COMING TO hu_mdd_real traverse \n" );
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
             //FOR Testing
             // if(a==1534)
             // printf("1534 children %d, lvl %d\n",A->d(z),kdn );
             // if(A->d(z)==1181)
             // printf("Parent 1181 %d lvl %d\n",a, k);
              //FOR Testing
            incomingedgecountHU[A->d(z)]--;
            // std::set<int> childpset;
            // compute_r(k-1,A->d(z),childpset);
            compute_rt(k-1,A->d(z));
            if(blockedBy[a]!=0){
                //FOR Testing
                // if(A->d(z)==1528)
                // printf("blockedBy %d\n",a );
                //FOR Testing
                if((blockedBy[A->d(z)]==0)||(blockedBy[A->d(z)]!=0&& argF->getNodeLevel(blockedBy[a])>argF->getNodeLevel(blockedBy[A->d(z)]))){
                    blockedBy[A->d(z)]=blockedBy[a];
                    //FOR Testing
                    if(A->d(z)==1528)
                    printf("blockedBy set to %d\n",blockedBy[a] );
                    //FOR Testing
                }
            }
            if(incomingedgecount[a]!=0){
                //FOR Testing
                // if(A->d(z)==1528)
                // printf("incomingedgecount[a] %d %d\n",a, incomingedgecount[a] );
                //FOR Testing
                if((blockedBy[A->d(z)]==0)||(blockedBy[A->d(z)]!=0&& argF->getNodeLevel(a)>argF->getNodeLevel(blockedBy[A->d(z)]))){
                    blockedBy[A->d(z)]=a;
                    //FOR Testing
                    // if(A->d(z)==1528)
                    // printf("blockedBy set to %d\n",a );
                    //FOR Testing
                }
            }
            if(incomingedgecountHU[A->d(z)]==0)
            {
                //FOR Testing
                // if(A->d(z)==1528)
                // printf("ADDED %d to pset of %d \n",int(A->d(z)), a );
                //FOR Testing
                // pset.insert(int(A->d(z))-1);
                // pset.insert(childpset.begin(), childpset.end());
                pset[a].insert(A->d(z));
                pset[a].insert(pset[A->d(z)].begin(), pset[A->d(z)].end());
                //FOR Testing
               //  if(A->d(z)==1528)
               //
               //  for (auto it=pset[a].begin(); it != pset[a].end(); ++it)
               // printf("Pset %d is**** %d\n",a,*it );
               //FOR Testing

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
 // }
           }
      return 0;
#endif
  }
 void MEDDLY::hu_mdd_real::initialize(int k, node_handle a, node_handle root)

 {
     #ifdef HRec
    int kdn = k - 1;
    unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
    for (unsigned z = 0; z < A->getNNZs(); z++) {
        if (kdn>0){
            if((a==root)||(LSA[a]==root)){
                setToRoot[A->d(z)]=true;
            }
            if(firstParent[A->d(z)]==-1){
                firstParent[A->d(z)]=a;
            }else{
                if(firstParent[A->d(z)]!=a){
                    singleParent[A->d(z)]=false;
                }
            }
            incomingedgecountHU[A->d(z)]--;
            if(incomingedgecountHU[A->d(z)]==0){
                if(singleParent[A->d(z)]){
                    LSA[A->d(z)]=firstParent[A->d(z)];
                }else{
                    if(setToRoot[A->d(z)]){
                        LSA[A->d(z)]=root;
                    }
                }
                initialize(kdn,A->d(z),root);
            }
        }
    }
    unpacked_node::recycle(A);
    #endif
 }

 void MEDDLY::hu_mdd_real::compute_r(int k, node_handle a,node_handle root)
 {
     #ifdef HRec
    int kdn = k - 1;

    unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
    for (unsigned z = 0; z < A->getNNZs(); z++) {
        if (kdn>0){
            incomingedgecountHU[A->d(z)]--;
            if(LSA[A->d(z)]==-1){
                if(LSA[a]==root||a==root){
                    LSA[A->d(z)]=root;
                    parents[A->d(z)].clear();
                }else{
                parents[A->d(z)].insert(a);
                if(LSA[a]==-1 && a!=root){printf("ERR Parent is %d LSA is -1\n",a );getchar();}
                }
            }
            if(incomingedgecountHU[A->d(z)]==0){
                if(LSA[A->d(z)]==-1){
                    // printf("Call LowestSeparatorAboveSet for %d\n",A->d(z) );
                    for(auto p:parents[A->d(z)]){
                        if(LSA[p]==-1)
                        { printf("LSA %d is -1\n", p);
                        getchar();}
                    }
                    LSA[A->d(z)]=LowestSeparatorAboveSet(kdn,parents[A->d(z)],root);
                }
                compute_r(kdn,A->d(z),root);
                parents[A->d(z)].clear();
            }
        }
    }
    unpacked_node::recycle(A);
    #endif
 }

     // printf("COMING TO hu_mdd_real \n" );
//      if(explored[a]==1){
//          // printf("EXP %d is true\n",a-1 );
//          // printf("else  %d\n",a-1);
//
//          lastPosition[a]=c;
//          c++;
//         //  for (auto it=pset[a-1].begin(); it != pset[a-1].end(); ++it)
//         // printf("else Pset is %d\n",*it );
//         //  printf("Done else  %d\n",a-1);
//      }else{
//
//          // printf("C is %d\n",c );
//          startInterval[a]=c;
//          // printf("startInterval[ %d ] =%d\n",a-1,c );
//          c++;
//          int kdn = k - 1;
//          unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
//          for (unsigned z = 0; z < A->getNNZs(); z++) {
//              if (kdn>0)
//              {
//
//
//                 incomingedgecountHU[A->d(z)]--;
//
//                 // std::set<int> childpset;
//                 // compute_r(k-1,A->d(z),childpset);
//                 compute_r(kdn,A->d(z));
//                 bool singleParent=false;
//                 if(incomingedgecount[A->d(z)]==1){
//                     if(A->d(z)==352)
//                     {
//                         printf("single parent a is %d\n",a );
//                         getchar();
//                     }
//                            updateInsert(a,A->d(z));
//                            singleParent=true;
//                            pset[a].insert(pset[A->d(z)].begin(), pset[A->d(z)].end());
//
//                        }
//                 else if(incomingedgecountHU[A->d(z)]==0)
//                 {
//                     if(A->d(z)==353)
//                     {
//                         printf("a is %d\n",a );
//                         getchar();
//                     }
//                      // printf("ADDED %d to pset of %d \n",int(A->d(z))-1, a-1 );
//                     // pset.insert(int(A->d(z))-1);
//                     // pset.insert(childpset.begin(), childpset.end());
//                     pset[a].insert(A->d(z));
//                     pset[a].insert(pset[A->d(z)].begin(), pset[A->d(z)].end());
//                    //  for (auto it=pset[a-1].begin(); it != pset[a-1].end(); ++it)
//                    // printf("Pset is**** %d\n",*it );
//
//                 }
//
//              }
//          }
//          stopInterval[a]=c;
//          // printf("stopInterval[ %d ] =%d\n",a-1,c );
//
//          c++;
//          explored[a]=true;
//          // printf("explored[ %d ] is true\n",a-1 );
//          unpacked_node::recycle(A);
//      }
//      if(incomingedgecountHU[a]==0)
//      {
//          unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
//          bool dooperation=true;
//          for (unsigned z = 0; z < A->getNNZs(); z++) {
//              if(incomingedgecountHU[A->d(z)]!=0)
//              dooperation=false;
//          }
//          if(a==14210)
//          printf("dooperation 14210 %d\n",dooperation );
//              unpacked_node::recycle(A);
//
//          // printf("incomingedgecountHU %d  is ZERO\n", a-1);
//          // std::set<int> result;
//         //  for (auto it=pset[a-1].begin(); it != pset[a-1].end(); ++it)
//         // printf("Pset is %d\n",*it );
//     // std::set<int> toRemove= CheckInterval(pset,a-1);
//     if(dooperation){
//         if(a==14210)
//         printf("dooperation 14210\n" );
//     std::set<int> toRemove= CheckInterval(a);
//
//     for (auto r: toRemove ){
//         pset[a].erase(r);
//     }
// }
//           }
 //     return 0;
 // }
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

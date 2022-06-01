
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
#include <list>
#include <map>
//#include "mpz_object.h"

// #define DEBUG_AC
// #define LRec
// #define LTrav //Wrong never should be turned on
#define Ldom

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

    #ifdef LRec
    res = compute_r(argF->getDomain()->getNumVariables(), arg.getNode());
    #endif
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
    incomingedgecountHU=new int[card];
    RincomingedgecountHU=new int[card];

    #ifdef LRec
    LUreverse=new int[card];
    if(explored==0||LUreverse==0)
    {printf("ERROR IN making arrays\n" ); getchar();}
    #endif
    #ifdef LTrav
    lastPosition=new int[card];
    startInterval= new int[card];
    stopInterval=new int[card];
    blockedBy= new int[card];
    if(explored==0 ||RincomingedgecountHU==0 ||incomingedgecountHU==0 || lastPosition==0||startInterval==0|| stopInterval==0)
    {printf("ERROR IN making arrays\n" ); char c=getchar();}
    #endif

    for(int i=0;i<card;i++)
    {
        explored[i]=false;
        incomingedgecountHU[i]=0;
        RincomingedgecountHU[i]=0;
        #ifdef LRec
        LUreverse[i]=-1;
        #endif
        #ifdef LTrav
        lastPosition[i]=0;
        startInterval[i]=0;
        stopInterval[i]=0;
        blockedBy[i]=0;
        #endif
    }
    for(int i=0;i<card;i++){
    incomingedgecountHU[i]=incomingedgecount[i];
    }
    #ifdef LTrav
    pset = new std::set<int>[card];
    RMDD= new std::list<int>[card];
    for(int i=0;i<card;i++){
    pset[i].clear();
    RMDD[i].clear();
    }

    compute_RMDD(argF->getDomain()->getNumVariables(),arg.getNode());
    //Testing
      //   for (int i = 0; i < card; i++){
      //       printf("N %d NL %d\n",i,argF->getNodeLevel(i) );
      //       std::list<int>::iterator iter = RMDD[i].begin();
      //       while (iter != RMDD[i].end())
      //       {
      //           printf("%d NL:%d, ",*iter,argF->getNodeLevel(*iter) );
      //         // std::cout << '\t' << *iter << '\n';
      //         iter++;
      //       }
      //       printf("\n" );
      //   }
      // getchar();
    //Tesing end
    for(int i=0;i<card;i++){
        explored[i]=false;
    }
    res = compute_rt(argF->getDomain()->getNumVariables(), 1/*arg.getNode()*/ );
    for(int i=0;i<card;i++){
    pset[i].clear();
    RMDD[i].clear();
    }
    #endif
    #ifdef LRec
    initialize(argF->getDomain()->getNumVariables(), arg.getNode() );
    for(int i=0;i<card;i++)
    {
        explored[i]=false;
    }
    #endif
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
     #ifdef LRec
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
    #endif
    #ifdef Ldom
g=new std::vector<int>[card];// dominator
rg=new std::vector<int> [card]; // revese graph
bucket=new std::vector<int>[card]; // bucket stores the set of semidom
sdom=new int[card]; // semi dominator
par=new int[card]; // parent information
dom=new int[card]; //dominator
dsu=new int[card];
label=new int[card];
arr=new int[card];
rev=new int[card];
for(int i=0; i<card;i++){
  arr[i]=0;
}
compute_rgraph(argF->getDomain()->getNumVariables(),arg.getNode());
// getchar();
for(int i=0; i<card;i++){
  arr[i]=0;
}
T=0;
compute_rdom(1);
int n=T;
for(int i=n;i>=1;i--){
  for(int j=0; j<rg[i].size();j++){
    sdom[i]=std::min(sdom[i],sdom[Find(rg[i][j])]);
  }
  if(i>1) bucket[sdom[i]].push_back(i);
  for(int j=0;j<bucket[i].size();j++)
{
  int w = bucket[i][j];
  int v = Find(w);
  if(sdom[v]==sdom[w])dom[w]=sdom[w];
  else dom[w] = v;
}
if(i>1)Union(par[i],i);
}
for(int i=2;i<=n;i++)
 {
if(dom[i]!=sdom[i])
  dom[i]=dom[dom[i]];
  // if(rev[i]==arg.getNode()&& rev[dom[i]]==arg.getNode());
  // else
updateInsert(rev[dom[i]],rev[i]);
}
delete[] rg; // revese graph
delete[] bucket; // bucket stores the set of semidom
delete[] sdom; // semi dominator
delete[] par; // parent information
delete[] dom; //dominator
delete[] dsu;
delete[] label;
delete[] arr;
delete[] rev;

#endif


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
   // getchar();
     //   printf("[%d, \t %d ,\t %d ]\t \n",startInterval[i],stopInterval[i],lastPosition[i] );
     //     printf("\n" );
     //  }
     //  printf("**HU****\n" );
     // for(int i=0;i<card;i++){
     // pset[i].clear();
     // }
      delete[] incomingedgecountHU;
      delete[] RincomingedgecountHU;
     #ifdef LRec
     cache.clear();
     delete[] LUreverse;
     #endif
     delete[] explored;
     #ifdef LTrav
     delete[] pset;
     delete[] RMDD;
     delete[] startInterval;
     delete[] stopInterval;
     delete[] lastPosition;
     delete[] blockedBy;
     #endif

#ifdef DEBUG_LU
     for(int i=0;i<card;i++)
     {
    	 printf("LU %d %d \n",i, lowestunique[i]);
     }
#endif
   }
   void compute_r(int ht, node_handle a);
   void initialize(int ht, node_handle a);
   double compute_rt(int ht, node_handle a);
    std::set<int> CheckInterval(node_handle a );
    void compute_RMDD(int ht,node_handle a);
    void compute_rdom(node_handle a);
void compute_rgraph(int ht,node_handle a);
   // std::set<int> CheckInterval(std::set<int> pset,node_handle a );
   // std::set<int> CheckInterval(node_handle a );

 protected:
     bool* explored;
     int* RincomingedgecountHU;
     int* incomingedgecountHU;

     #ifdef LTrav

     std::set<int>* pset;
     std::list<int>* RMDD;
    int* lastPosition;
    int* startInterval;
    int* stopInterval;
    int* blockedBy;
    int c=1;

     #endif
     #ifdef LRec
     int* LUreverse;
     #endif
     #ifdef Ldom
int T=0;
// std::vector<int>* LSAtree;// dominator
std::vector<int>* g;
std::vector<int>* rg; // revese graph
std::vector<int>* bucket; // bucket stores the set of semidom
int* sdom; // semi dominator
int* par; // parent information
int* dom; //dominator
int* dsu;
int* label;
int* arr;
int* rev;
#endif
  inline void updateInsert(node_handle a, node_handle b){
      //printf("COMING TO updateInsert\n");
      // #ifdef LTrav
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
// #endif
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
  inline int Find(int u,int x=0)
 {
     #ifdef Ldom
   if(u==dsu[u])return x?-1:u;
   int v = Find(dsu[u],x+1);
   if(v<0)return u;
   if(sdom[label[dsu[u]]] < sdom[label[u]])
       label[u] = label[dsu[u]];
   dsu[u] = v;
   return x?v:label[u];
   #endif
 }
 inline void Union(int u,int v) //Add an edge u-->v
 {
       #ifdef Ldom
   dsu[v]=u; 	//yup,its correct :)
   #endif
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
      #ifdef LRec
      // printf("ChildsLeastCommonUnique\n");

      int firstChildLVL=argF->getNodeLevel(firstChild);
      int secondChildLVL=argF->getNodeLevel(secondChild);
      // printf("firstChildLVL %d secondChildLVL %d\n", firstChildLVL,secondChildLVL);
      // printf("LU firstChildLVL %d LU secondChildLVL %d\n",LUreverse[ firstChild],LUreverse[secondChild]);
      // if(firstChild==secondChild) return firstChild;
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
      #endif
  }
  inline int LeastCommonUnique(int lvl ,node_handle a)
  {
      #ifdef LRec
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
       #endif
  }
 };
 void MEDDLY::lu_mdd_real::compute_RMDD( int lvl,node_handle a){
     #ifdef LTrav
     explored[a] = true;
     if(lvl>0){
        int kdn=lvl-1;

        unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
        for(unsigned z = 0; z < A->getNNZs(); z++){
            if(A->d(z)>0){
                RMDD[A->d(z)].push_back(a);
                RincomingedgecountHU[a]++;
                if(!explored[A->d(z)]){
                    compute_RMDD(kdn,A->d(z));
                }
            }
        }
        unpacked_node::recycle(A);
    }
     #endif
 }

 std::set<int> MEDDLY::lu_mdd_real::CheckInterval(node_handle a ){
#ifdef LTrav
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
   if((blockedBy[l]==0)||(blockedBy[l]!=0 && RincomingedgecountHU[blockedBy[l]]==0)||(argF->getNodeLevel(l)+1==argF->getNodeLevel(a)))
   if((startInterval[l]>=startInterval[a])&&(stopInterval[l]<=stopInterval[a])
   &&(((lastPosition[l]!=0)&&(lastPosition[l]>=startInterval[a])&&(lastPosition[l]<=stopInterval[a]))||(lastPosition[l]==0))){
      //NOT needed RMDD
   // if((blockedBy[l]==0)||(blockedBy[l]!=0 && incomingedgecountHU[blockedBy[l]]==0)||(argF->getNodeLevel(l)+1==argF->getNodeLevel(a)))
   // if((startInterval[l]>=startInterval[a])&&(stopInterval[l]<=stopInterval[a])
   // &&(((lastPosition[l]!=0)&&(lastPosition[l]>=startInterval[a])&&(lastPosition[l]<=stopInterval[a]))||(lastPosition[l]==0))){
//NOT needed RMDD
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

double MEDDLY::lu_mdd_real::compute_rt(int k, node_handle a){
    #ifdef LTrav
     // printf("COMING TO lu_mdd_real traverse \n" );
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
   //NOT needed RMDD
   // unpacked_node* A = unpacked_node::newFromNode(argF, a, false);

   // for (unsigned z = 0; z < A->getNNZs(); z++) {
   //NOT needed RMDD
   for(auto const &z: RMDD[a]){
       if (kdn>0)
       {
           //FOR Testing
           // if(a==1534)
           // printf("1534 children %d, lvl %d\n",A->d(z),kdn );
           // if(A->d(z)==1181)
           // printf("Parent 1181 %d lvl %d\n",a, k);
            //FOR Testing

            //NOT needed RMDD
          // incomingedgecountHU[A->d(z)]--;
          //NOT needed RMDD
          RincomingedgecountHU[z]--;
          // std::set<int> childpset;
          // compute_r(k-1,A->d(z),childpset);

          //NOT needed RMDD
          // compute_rt(k-1,A->d(z));
          //NOT needed RMDD
           compute_rt(k-1,z);
          if(blockedBy[a]!=0){
              //FOR Testing
              // if(A->d(z)==1528)
              // printf("blockedBy %d\n",a );
              //FOR Testing
              if((blockedBy[z]==0)||(blockedBy[z]!=0&& argF->getNodeLevel(blockedBy[a])>argF->getNodeLevel(blockedBy[z]))){
                  blockedBy[z]=blockedBy[a];
                   //NOT needed RMDD
              // if((blockedBy[A->d(z)]==0)||(blockedBy[A->d(z)]!=0&& argF->getNodeLevel(blockedBy[a])>argF->getNodeLevel(blockedBy[A->d(z)]))){
              //     blockedBy[A->d(z)]=blockedBy[a];
               //NOT needed RMDD
                  //FOR Testing
                  // if(A->d(z)==1528)
                  // printf("blockedBy set to %d\n",blockedBy[a] );
                  //FOR Testing
              }
          }
          if(incomingedgecount[a]!=0){
              //FOR Testing
              // if(A->d(z)==1528)
              // printf("incomingedgecount[a] %d %d\n",a, incomingedgecount[a] );
              //FOR Testing
              if((blockedBy[z]==0)||(blockedBy[z]!=0&& argF->getNodeLevel(a)>argF->getNodeLevel(blockedBy[z]))){
                  blockedBy[z]=a;
                  //NOT needed RMDD
              // if((blockedBy[A->d(z)]==0)||(blockedBy[A->d(z)]!=0&& argF->getNodeLevel(a)>argF->getNodeLevel(blockedBy[A->d(z)]))){
              //     blockedBy[A->d(z)]=a;
              //NOT needed RMDD
                  //FOR Testing
                  // if(A->d(z)==1528)
                  // printf("blockedBy set to %d\n",a );
                  //FOR Testing
              }
          }
          if(RincomingedgecountHU[z]==0)
          //NOT needed RMDD
          // if(incomingedgecountHU[A->d(z)]==0)
          //NOT needed RMDD
          {
              //FOR Testing
              // if(A->d(z)==1528)
              // printf("ADDED %d to pset of %d \n",int(A->d(z)), a );
              //FOR Testing
              // pset.insert(int(A->d(z))-1);
              // pset.insert(childpset.begin(), childpset.end());
              //NOT needed RMDD
              // pset[a].insert(A->d(z));
              // pset[a].insert(pset[A->d(z)].begin(), pset[A->d(z)].end());
              //NOT needed RMDD
              pset[a].insert(z);
              pset[a].insert(pset[z].begin(), pset[z].end());
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
//NOT needed RMDD
// unpacked_node::recycle(A);
//NOT needed RMDD
}
 //NOT needed RMDD
// if(incomingedgecountHU[a]==0)
 //NOT needed RMDD
 if(RincomingedgecountHU[a]==0)
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
 void MEDDLY::lu_mdd_real::initialize(int k, node_handle a)

 {
     #ifdef LRec
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
     #endif
 }


 void MEDDLY::lu_mdd_real::compute_r(int k, node_handle a)
 {
     #ifdef LRec
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
     #endif
 }

 void MEDDLY::lu_mdd_real::compute_rgraph(int k, node_handle a){
   #ifdef Ldom
   if (k>0){
   arr[a]=T;

   unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
   int kdn = k - 1;
     for (unsigned z = 0; z < A->getNNZs(); z++) {
         if(A->d(z)>0){

     // for(int i=0;i<g[a].size();i++){
       int w=A->d(z);
       if(!arr[w]){
           // if(a==224){
           //     printf("*Node %d child %d\n",a,A->d(z) );
           //     // getchar();
           // }
         compute_rgraph(kdn,w);
         // par[arr[w]]=arr[a];
       }
       g[w].push_back(a);
   }
     }
   unpacked_node::recycle(A);
}

   #endif
 }

 void MEDDLY::lu_mdd_real::compute_rdom( node_handle a){
   #ifdef Ldom
   T++; arr[a]=T; rev[T]=a;
   label[T]=T; sdom[T]=T; dsu[T]=T;
   // unpacked_node* A = unpacked_node::newFromNode(argF, a, false);
     // for (unsigned z = 0; z < A->getNNZs(); z++) {
     for(int i=0;i<g[a].size();i++){
       int w=g[a][i];
       if(!arr[w]){
         compute_rdom(w);
         par[arr[w]]=arr[a];
       }
       rg[arr[w]].push_back(arr[a]);
     }
     // }
   // unpacked_node::recycle(A);
   #endif
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

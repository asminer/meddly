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
 #include "../defines.h"
 #include "covr.h"
 #include <typeinfo> // for "bad_cast" exception
 #include <set>
 #include <list>

 #include "../minterms.h"
 #include "../ct_entry_result.h"
 #include "../compute_table.h"
 #include "../oper_unary.h"
 #include "../oper_binary.h"
 #include "../oper_special.h"
 #include "../opname_satur.h"
 #include "../ops_builtin.h"
 #include "../oper_unary.h"

 namespace MEDDLY {
 // class cov_by_events_opname;
 // class cov_by_events_op;
 //
 // class common_cov_by_events_mt;
 // class cov_by_events_mt;
 //
 // class cov_opname;
 class compareij;
 class covrExtractCover_opname;
 class covrExtractCoverMT;

 class covrExtractFrom_opname;
 class covrExtractFromMT;
 };

 // ******************************************************************
 // *                                                                *
 // *                        covrExtractFromMT class                 *
 // *                                                                *
 // ******************************************************************
 class MEDDLY::covrExtractFromMT:  public unary_operation {
 public:
   covrExtractFromMT(unary_opname* oc, expert_forest* arg, expert_forest* res);

   virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);
   binary_operation* mddUnion;

 protected:
   virtual void compute_r( node_handle a,int k, node_handle& b);

   inline ct_entry_key*
   findResult(node_handle a,node_handle &b)
   {
     ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
     MEDDLY_DCASSERT(CTsrch);
     CTsrch->writeN(a);
     CT0->find(CTsrch, CTresult[0]);
     if (!CTresult[0]) return CTsrch;
     b = resF->linkNode(CTresult[0].readN());
     CT0->recycle(CTsrch);
     return 0;
   }
   inline node_handle saveResult(ct_entry_key* Key,
   node_handle a,node_handle b)
   {
     CTresult[0].reset();
     CTresult[0].writeN(b);
     CT0->addEntry(Key, CTresult[0]);
     return b;
   }
};
MEDDLY::covrExtractFromMT::covrExtractFromMT(unary_opname* oc, expert_forest* arg, expert_forest* res)
  : unary_operation(oc, 1, arg, res)
{
  MEDDLY_DCASSERT(!resF->isForRelations());

  ct_entry_type* et = new ct_entry_type(oc->getName(), "N:N");
  et->setForestForSlot(0, arg);
  et->setForestForSlot(2, res);
  registerEntryType(0, et);
  buildCTs();
  dd_edge er(resF);
  mddUnion=getOperation(UNION, er,er,er);
  MEDDLY_DCASSERT(mddUnion);
}
void MEDDLY::covrExtractFromMT::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{

  node_handle bnode = 0;
  compute_r( arg.getNode(),resF->getNumVariables(), bnode);
  // printf("bnode is %d\n",bnode );
  res.set(bnode);
  // ostream_output meddlyout(std::cout);
  // const int mddLevel = resF->getNodeLevel(bnode);
  // printf("mddLevel %d\n",mddLevel );
  // unpacked_node* A=resF->newUnpacked(bnode,FULL_ONLY);
  // A->show(meddlyout,true);
  // printf("BNODE^^^\n" );
}

void MEDDLY::covrExtractFromMT::compute_r( node_handle a,int k, node_handle& res)
{
    ostream_output meddlyout(std::cout);
    // printf("covrExtractFromMT : a is %d\n",a );
    if ( (!resF->isQuasiReduced()||k==0)&& argF->isTerminalNode(a)){
        if(a==0)res=0;
        else {
            long val;
            argF->getValueFromHandle(a, val);
            res=resF->handleForValue(val);
        }
        return;
    }
    if(ABS(argF->getNodeLevel(a))<k){
        int size=resF->getLevelSize(k);
        unpacked_node* T =unpacked_node::newFull(resF,k,size);
        node_handle newstates=0;
        compute_r( a,k-1, newstates);
        for(int i=0;i<size; i++){
            T->d_ref(i)=newstates;
        }
        res=resF->createReducedNode(-1,T);
        return;
    }
    // check the cache
    ct_entry_key* Key = findResult(a,res);
    if (0==Key) {
      return;
    }
    const int aLevel = argF->getNodeLevel(a);
    const int rLevel=ABS(aLevel);
    const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
    unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);
    for (unsigned j = 0; j < rSize; j++) {
      C->d_ref(j) = 0;
    }
    unpacked_node* A =  aLevel<0
      ? unpacked_node::newRedundant(argF, rLevel,a, true)
      : argF->newUnpacked(a, FULL_ONLY);

    for (unsigned i = 0; i < rSize; i++) {
        // printf("i %d, a %d,  A->d(i) %d \n",i,a,A->d(i) );
        int pLevel = argF->getNodeLevel(A->d(i));
        // printf("rLevel %d pLevel %d\n",rLevel,pLevel );
        unpacked_node* B = isLevelAbove(-rLevel, pLevel)
        ? unpacked_node::newIdentity(argF, -rLevel, i, A->d(i), true)
        : argF->newUnpacked(A->d(i), FULL_ONLY);
        dd_edge newstatesE(resF), djp(resF);

        for (unsigned j = 0; j < rSize; j++) {
            if (0 == B->d(j)) {
              continue;
            }
            // printf("j %d\n",j );
            // printf("B->d(j) %d \n", B->d(j) );
            node_handle newstates = 0;
            compute_r( B->d(j),k-1, newstates);
            if (0==newstates) {
              continue;
            }
            // printf("newstates  %d\n",newstates );
            // getchar();
            newstatesE.set(newstates);
            djp.set(C->d(i));
            mddUnion->computeTemp(newstatesE, djp, djp);
            C->set_d(i, djp);

            // C->d_ref(i) = newstates;
            // C->show(meddlyout,true);
            // printf("\nADDED %d %d\n",i,newstates );
        }
        unpacked_node::recycle(B);
        // printf("recycling B\n" );
    }
    // cleanup mdd reader
    unpacked_node::recycle(A);
    // printf("recycling A\n" );
    // printf("SHOW C\n" );
    // C->show(meddlyout,true);
    // getchar();
    res=resF->createReducedNode(-1, C);
    // printf("SET result%d\n", res );

    // #ifdef TRACE_ALL_OPS
    //   printf("computed new tcXrel(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmxd, mxd, resEv, resEvmdd);
    // #endif
    saveResult(Key, a,res);
    printf("Save result\n" );
}

// ******************************************************************
// *                                                                *
// *                       covrExtractCoverMT class                 *
// *                                                                *
// ******************************************************************
enum RV
{
 G=-1,
 E= -2,
 L= -3,
 N= -4

};
class MEDDLY::covrExtractCoverMT:  public unary_operation {
public:
  covrExtractCoverMT(unary_opname* oc, expert_forest* arg, expert_forest* res);

  virtual void computeDDEdge(const dd_edge &arg, dd_edge &res,markcmp* cij, bool userFlag);
  binary_operation* mddUnion;
   void compute_r( node_handle a,node_handle b,RV above,markcmp* _cij, node_handle& c);

protected:
  virtual void compute_r( node_handle a,int k, node_handle& b){};
  markcmp* cij;
  inline ct_entry_key*
  findResult(node_handle a,node_handle b,int c,node_handle &d)
  {
    ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeN(a);
    CTsrch->writeN(b);
    CTsrch->writeI(c);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    d = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
  }
  inline node_handle saveResult(ct_entry_key* Key,
  node_handle a,node_handle b,int c,node_handle d)
  {
    CTresult[0].reset();
    CTresult[0].writeN(d);
    CT0->addEntry(Key, CTresult[0]);
    return b;
  }
};
MEDDLY::covrExtractCoverMT::covrExtractCoverMT(unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, 1, arg, res)
{
 MEDDLY_DCASSERT(!resF->isForRelations());

 ct_entry_type* et = new ct_entry_type(oc->getName(), "NNI:N");
 et->setForestForSlot(0, arg);
 et->setForestForSlot(1,arg);
 et->setForestForSlot(4, res);
 registerEntryType(0, et);
 buildCTs();
 dd_edge er(resF);
 mddUnion=getOperation(UNION, er,er,er);
 MEDDLY_DCASSERT(mddUnion);
}
void MEDDLY::covrExtractCoverMT::computeDDEdge(const dd_edge &arg, dd_edge &res,markcmp* _cij, bool userFlag)
{

    node_handle bnode = 0;
    RV above=E;
    compute_r( arg.getNode(),arg.getNode(),above,_cij,  bnode);
    printf("bnode is %d\n",bnode );
    cij=_cij;
    res.set(bnode);
    // // ostream_output meddlyout(std::cout);
    ostream_output meddlyout(std::cout);
    res.showGraph(meddlyout);
    printf("RESULT COVR EXTRACTCOVERED\n" );
    // // argF->showNodeGraph(meddlyout, arg.getNode());
    // // arg.show(stdout, 2);
    // printf("NODE!\n" );
    // printf("arg.getNode() %d\n",arg.getNode() );
    getchar();

}

void MEDDLY::covrExtractCoverMT::compute_r( node_handle a,node_handle b,RV above,markcmp* cmp, node_handle& res)
{
   ostream_output meddlyout(std::cout);
   printf("covrExtractCoverMT a %d b %d\n",a,b );

   if(a==0||b==0){
       res=0;
       return;
   }
   if(argF->isTerminalNode(a)&&argF->isTerminalNode(b)){
       printf("terminal %d %d\n",a,b );
       if(a==0||b==0) res=0;
       else if(above==L){
            long aval;
            resF->getValueFromHandle(a, aval);
            long bval;
            resF->getValueFromHandle(b, bval);
            res=resF->handleForValue(aval*bval);
            printf("res%d\n",aval*bval );
       }else{
           res=0;
       }
       return;
   }
   //check the cache
   int iabove=above;
   ct_entry_key* Key=findResult(a,b,above,res);
   if(0==Key){
       return;
   }

   const int aLevel=resF->getNodeLevel(a);
   const int bLevel=resF->getNodeLevel(b);
   const int rLevel=MAX(ABS(aLevel),ABS(bLevel));
   if (aLevel<0||bLevel<0){
       printf("ERROR a %d b %d aLevel %d bLevel %d\n",a,b,aLevel,bLevel );
   }
   const unsigned rSize= unsigned(resF->getLevelSize(rLevel));
   unpacked_node* C=unpacked_node::newFull(resF,rLevel,rSize);
   for (unsigned j = 0; j < rSize; j++) {
     C->d_ref(j) = 0;
   }
   unpacked_node* A=resF->newUnpacked(a,FULL_ONLY);
   // A->show(meddlyout,true);
   // printf("\n a^^^\n" );
   unpacked_node* B = resF->newUnpacked(b, FULL_ONLY);
   // B->show(meddlyout,true);
   // printf("\n b^^^\n" );

   //////////////////////////////
   for(unsigned i=0; i<rSize; i++){
       // printf("i %d\n",i );
       if(rLevel>aLevel){
           printf("INSIDE if\n" );
           dd_edge newstatesE(resF), djp(resF);
           for(unsigned j=0; j<rSize; j++){
               RV nabove=above;
               int compareij=cmp->compare(i,j,rLevel);
               if(i==j);
               else if(compareij>0){
                   if(above==L||above==E){nabove=L;}
                   else nabove=N;
               }else {nabove=N;}
               node_handle newstates=0;
               // printf("rec call if\n" );
               compute_r(a,B->d(j),nabove, cmp, newstates);
               if(newstates==0){continue;}
               newstatesE.set(newstates);
               djp.set(C->d(i));
               mddUnion->computeTemp(newstatesE, djp, djp);
               C->set_d(i, djp);
           }
       }
       else{
           dd_edge newstatesE(resF), djp(resF);

           for(unsigned j=0; j<rSize; j++){
               printf("ELSE a %d b %d i%d j %d rLevel %d rsize %d\n",a,b,i,j,rLevel,rSize );
               if(0==B->d(j)){
                   continue;
               }
               RV nabove=above;
               int compareij=cmp->compare(i,j,rLevel);
               if(i==j);
               else if(compareij>0){
                   if(above==L||above==E){nabove=L;}
                   else nabove=N;
               }else {nabove=N;}
               node_handle newstates=0;
               // printf("i %d, j %d\n",i,j );
               // printf("rec call else %d, %d \n",A->d(i),B->d(j) );
               compute_r(A->d(i),B->d(j),nabove, cmp, newstates);
               if(newstates==0){
                   printf("ZERO\n" );
                   continue;
               }
               printf(" NOT ZERO\n" );
               newstatesE.set(newstates);
               djp.set(C->d(i));
               mddUnion->computeTemp(newstatesE, djp, djp);
               C->set_d(i, djp);
           }
       }
   }


   unpacked_node::recycle(A);
   unpacked_node::recycle(B);
   res=resF->createReducedNode(-1,C);
   saveResult(Key, a,b,above,res);



   // printf("covrExtractFromMT : a is %d\n",a );
   // if ( (!resF->isQuasiReduced()||k==0)&& argF->isTerminalNode(a)){
   //     if(a==0)res=0;
   //     else {
   //         long val;
   //         argF->getValueFromHandle(a, val);
   //         res=resF->handleForValue(val);
   //     }
   //     return;
   // }
   // if(ABS(argF->getNodeLevel(a))<k){
   //     int size=resF->getLevelSize(k);
   //     unpacked_node* T =unpacked_node::newFull(resF,k,size);
   //     node_handle newstates=0;
   //     compute_r( a,b,above,k-1, newstates);
   //     for(int i=0;i<size; i++){
   //         T->d_ref(i)=newstates;
   //     }
   //     res=resF->createReducedNode(-1,T);
   //     return;
   // }
   // // check the cache
   // ct_entry_key* Key = findResult(a,res);
   // if (0==Key) {
   //   return;
   // }
   // const int aLevel = argF->getNodeLevel(a);
   // const int rLevel=ABS(aLevel);
   // const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
   // unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);
   // for (unsigned j = 0; j < rSize; j++) {
   //   C->d_ref(j) = 0;
   // }
   // unpacked_node* A =  aLevel<0
   //   ? unpacked_node::newRedundant(argF, rLevel,a, true)
   //   : argF->newUnpacked(a, FULL_ONLY);
   //
   // for (unsigned i = 0; i < rSize; i++) {
   //     // printf("i %d, a %d,  A->d(i) %d \n",i,a,A->d(i) );
   //     int pLevel = argF->getNodeLevel(A->d(i));
   //     // printf("rLevel %d pLevel %d\n",rLevel,pLevel );
   //     unpacked_node* B = isLevelAbove(-rLevel, pLevel)
   //     ? unpacked_node::newIdentity(argF, -rLevel, i, A->d(i), true)
   //     : argF->newUnpacked(A->d(i), FULL_ONLY);
   //     dd_edge newstatesE(resF), djp(resF);
   //
   //     for (unsigned j = 0; j < rSize; j++) {
   //         if (0 == B->d(j)) {
   //           continue;
   //         }
   //         // printf("j %d\n",j );
   //         // printf("B->d(j) %d \n", B->d(j) );
   //         node_handle newstates = 0;
   //         compute_r( B->d(j),b,above,k-1, newstates);
   //         if (0==newstates) {
   //           continue;
   //         }
   //         // printf("newstates  %d\n",newstates );
   //         // getchar();
   //         newstatesE.set(newstates);
   //         djp.set(C->d(i));
   //         mddUnion->computeTemp(newstatesE, djp, djp);
   //         C->set_d(i, djp);
   //
   //         // C->d_ref(i) = newstates;
   //         // C->show(meddlyout,true);
   //         // printf("\nADDED %d %d\n",i,newstates );
   //     }
   //     unpacked_node::recycle(B);
   //     // printf("recycling B\n" );
   // }
   // // cleanup mdd reader
   // unpacked_node::recycle(A);
   // // printf("recycling A\n" );
   // // printf("SHOW C\n" );
   // // C->show(meddlyout,true);
   // // getchar();
   // res=resF->createReducedNode(-1, C);
   // // printf("SET result%d\n", res );
   //
   // // #ifdef TRACE_ALL_OPS
   // //   printf("computed new tcXrel(<%ld, %d>, %d) = <%ld, %d>\n", ev, evmxd, mxd, resEv, resEvmdd);
   // // #endif
   // saveResult(Key, a,res);
   // printf("Save result\n" );
}

 // ******************************************************************
 // *                                                                *
 // *                   covrExtractFrom_opname class                 *
 // *                                                                *
 // ******************************************************************

 class MEDDLY::covrExtractFrom_opname : public unary_opname {

 public:
 covrExtractFrom_opname();
 virtual unary_operation* buildOperation(expert_forest* ar, expert_forest* res);
 };

 MEDDLY::covrExtractFrom_opname::covrExtractFrom_opname()
         : unary_opname("Extract From")
 {
 }

 MEDDLY::unary_operation*
 MEDDLY::covrExtractFrom_opname::buildOperation(expert_forest* arg, expert_forest* res)
 {
         if (arg->isForRelations() )
         {
             printf("covrExtractFrom_opname CORRECT!\n" );
            return new covrExtractFromMT(this,arg,res);
         }

 };



 // ******************************************************************
 // *                                                                *
 // *                  covrExtractCover_opname class                 *
 // *                                                                *
 // ******************************************************************

 class MEDDLY::covrExtractCover_opname : public unary_opname {

 public:
 covrExtractCover_opname();
 virtual unary_operation* buildOperation(expert_forest* ar, expert_forest* res);
 };

 MEDDLY::covrExtractCover_opname::covrExtractCover_opname()
         : unary_opname("Extract Cover")
 {
 }

 MEDDLY::unary_operation*
 MEDDLY::covrExtractCover_opname::buildOperation(expert_forest* arg, expert_forest* res)
 {
    if (arg->isForRelations() )
    {
        printf("covrExtractCover_opname CORRECT!\n" );
        return new covrExtractCoverMT(this,arg,res);
    }

 };


 // ******************************************************************
 // *                                                                *
 // *                           Front  end                           *
 // *                                                                *
 // ******************************************************************

 MEDDLY::unary_opname* MEDDLY::initializeExtractFrom()
 {
         printf("This is COVR EXTRACTFROM!\n" );
         return new covrExtractFrom_opname();
 }

 MEDDLY::unary_opname* MEDDLY::initializeExtractCovered(){
        printf("This is Extract Covered!\n" );
        return new covrExtractCover_opname();
 }

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

 class covrCoveredFrom_opname;
 class covrCoveredFromMT;

 class covrCoveredTo_opname;
 class covrCoveredToMT;
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
    // printf("Save result\n" );
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
    bool debug=false;
    //////
    ostream_output meddlyout(std::cout);
    // unpacked_node* A=argF->newUnpacked(arg.getNode(),FULL_ONLY);
    // A->show(meddlyout,true);
    // printf("\n a^^^\n" );
    // getchar();
    //////
    compute_r( arg.getNode(),arg.getNode(),above,_cij,  bnode);
    // printf("bnode is %d\n",bnode );
    cij=_cij;
    res.set(bnode);
    // // ostream_output meddlyout(std::cout);
    // ostream_output meddlyout(std::cout);
    if(debug){
    res.showGraph(meddlyout);
    printf("RESULT COVR EXTRACTCOVERED\n" );
    }
    // // argF->showNodeGraph(meddlyout, arg.getNode());
    // // arg.show(stdout, 2);
    // printf("NODE!\n" );
    // printf("arg.getNode() %d\n",arg.getNode() );
    // getchar();

}

void MEDDLY::covrExtractCoverMT::compute_r( node_handle a,node_handle b,RV above,markcmp* cmp, node_handle& res)
{
    bool debug=false;
    if(debug){
   ostream_output meddlyout(std::cout);
   printf("covrExtractCoverMT a %d b %d\n",a,b );
    }
   if(a==0||b==0){
       res=0;
       return;
   }
   if(argF->isTerminalNode(a)&&argF->isTerminalNode(b)){
       if(debug){
       printf("terminal %d %d\n",a,b );
        }
       if(a==0||b==0) res=0;
       else if(above==L){
            long aval;
            argF->getValueFromHandle(a, aval);
            long bval;
            argF->getValueFromHandle(b, bval);
            res=resF->handleForValue(aval*bval);
            if(debug){
            printf("res%d\n",aval*bval );
            }
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

   const int aLevel=argF->getNodeLevel(a);
   const int bLevel=argF->getNodeLevel(b);
   const int rLevel=MAX(ABS(aLevel),ABS(bLevel));
   if (aLevel<0||bLevel<0){
       printf("ERROR a %d b %d aLevel %d bLevel %d\n",a,b,aLevel,bLevel );
   }
   const unsigned rSize= unsigned(resF->getLevelSize(rLevel));
   unpacked_node* C=unpacked_node::newFull(resF,rLevel,rSize);
   for (unsigned j = 0; j < rSize; j++) {
     C->d_ref(j) = 0;
   }
   unpacked_node* A=argF->newUnpacked(a,FULL_ONLY);
   // A->show(meddlyout,true);
   // printf("\n a^^^\n" );
   unpacked_node* B = argF->newUnpacked(b, FULL_ONLY);
   // B->show(meddlyout,true);
   // printf("\n b^^^\n" );

   //////////////////////////////
   for(unsigned i=0; i<rSize; i++){
       // printf("i %d\n",i );
       if(rLevel>aLevel){
           if(debug){
           printf("INSIDE if\n" );
            }
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
               if(debug){
               printf("ELSE a %d b %d i%d j %d rLevel %d rsize %d\n",a,b,i,j,rLevel,rSize );
                }
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
                   if(debug){
                   printf("ZERO\n" );
                    }
                   continue;
               }
               if(debug){
               printf(" NOT ZERO\n" );
                }
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
}

// ******************************************************************
// *                                                                *
// *                       covrCoveredFromMT class                  *
// *                                                                *
// ******************************************************************
class MEDDLY::covrCoveredFromMT:  public binary_operation {
public:
  covrCoveredFromMT(binary_opname* oc, expert_forest* arg1,expert_forest* arg2, expert_forest* res);

  virtual void computeDDEdge(const dd_edge &arg, const dd_edge & arg2, dd_edge &res, bool userFlag);
  binary_operation* mddUnion;
   void compute_r( node_handle a,node_handle b, node_handle& c);

protected:
  inline ct_entry_key*
  findResult(node_handle a,node_handle b,node_handle &d)
  {
    ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeN(a);
    CTsrch->writeN(b);
    // CTsrch->writeI(c);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    d = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
  }
  inline node_handle saveResult(ct_entry_key* Key,
  node_handle a,node_handle b,node_handle d)
  {
    CTresult[0].reset();
    CTresult[0].writeN(d);
    CT0->addEntry(Key, CTresult[0]);
    return b;
  }
};
MEDDLY::covrCoveredFromMT::covrCoveredFromMT(binary_opname* oc, expert_forest* arg,expert_forest* arg2, expert_forest* res)
 : binary_operation(oc, 1, arg,arg2, res)
{
 MEDDLY_DCASSERT(resF->isForRelations());

 ct_entry_type* et = new ct_entry_type(oc->getName(), "NN:N");
 et->setForestForSlot(0, arg);
 et->setForestForSlot(1,arg);
 et->setForestForSlot(3, res);
 registerEntryType(0, et);
 buildCTs();
 dd_edge er(resF);
 mddUnion=getOperation(UNION, er,er,er);
 MEDDLY_DCASSERT(mddUnion);
}
void MEDDLY::covrCoveredFromMT::computeDDEdge(const dd_edge &arg, const dd_edge & arg2, dd_edge &res, bool userFlag)
{

    // printf("covrCoveredFromMT computeDDEdge\n" );
    node_handle bnode = 0;
    // RV above=E;
    // //////
    ostream_output meddlyout(std::cout);
    //
    compute_r( arg.getNode(),arg2.getNode(),  bnode);

    // printf("bnode is %d\n",bnode );
    // cij=_cij;
    res.set(bnode);
    bool debug=false;
    if(debug){
    arg.showGraph(meddlyout);
    printf("^^^arg1^^^\n" );
    arg2.showGraph(meddlyout);
    printf("arg2^^^^\n" );

    res.showGraph(meddlyout);
    printf("res\n" );
    printf("RESULT COVR coveredFrom\n" );
    //
    getchar();
    }

}

void MEDDLY::covrCoveredFromMT::compute_r( node_handle a,node_handle e, node_handle& res)
{   ostream_output meddlyout(std::cout);

    bool debug=false;
    if(debug){
   printf("covrCoveredFromMT a %d e %d\n",a,e );
    }

   if(a==0||e==0){
       res=0;
       return;
   }
   if(arg1F->isTerminalNode(a)&&arg2F->isTerminalNode(e)){
       if(debug){
       printf("terminal %d %d\n",a,e );
        }
       if(a==0||e==0) res=0;
       else{
            long aval;
            arg1F->getValueFromHandle(a, aval);
            long eval;
            arg2F->getValueFromHandle(e, eval);
            res=resF->handleForValue(aval*eval);
            if(debug){
            printf("res%ld\n",aval*eval );
            }
       }
       return;
   }
   //check the cache
   ct_entry_key* Key=findResult(a,e,res);
   if(0==Key){
       return;
   }

   const int aLevel=arg1F->getNodeLevel(a);
   const int eLevel=arg2F->getNodeLevel(e);
   const int rLevel=MAX(ABS(aLevel),ABS(eLevel));
   if(debug){
   printf("alevel %d, elevel %d\n",aLevel,eLevel );
    }

   const unsigned rSize= unsigned(resF->getLevelSize(rLevel));
    const unsigned ESize= unsigned(arg2F->getLevelSize(eLevel));
   unpacked_node* C=unpacked_node::newFull(resF,rLevel,rSize);

   unpacked_node* A = isLevelAbove(rLevel, aLevel)
    ? unpacked_node::newRedundant(arg1F, rLevel, a, true)
    : arg1F->newUnpacked(a, FULL_ONLY);
    unpacked_node* E = arg2F->newUnpacked(e, FULL_ONLY);
    if(debug){
    A->show(meddlyout,true);
    printf("\n a^^^\n" );
    }
    for (unsigned i = 0; i < rSize; i++) {
        // if(E->d(i)==0) continue;
        if(debug){
        printf("i %d a%d e %d\n",i,a, e );
        }
        int pLevel = arg1F->getNodeLevel(A->d(i));
        if(debug){
        printf("pLevel %d\n",pLevel );
        }
        unpacked_node* B = isLevelAbove(-rLevel, pLevel)
        ? unpacked_node::newIdentity(arg1F, -rLevel, i,  A->d(i), true)
        : arg1F->newUnpacked(A->d(i), FULL_ONLY);

        unpacked_node* D = unpacked_node::newFull(resF, -rLevel, rSize);
        if (rLevel > ABS(aLevel)) {
            if(debug){
            printf("IF\n" );
            }
            for(unsigned j = 0; j < rSize; j++){
                node_handle nres=0;
                // compute_r(B->d(j),e,nres);
                compute_r(a,E->d(j),nres);
                D->d_ref(j)=nres;
            }
        }else{
            if(debug){
            printf("ELSE\n" );
            }
            MEDDLY_DCASSERT(ABS(eLevel) >= ABS(pLevel));
            for (unsigned j = 0; j < rSize; j++) {
              D->d_ref(j) = 0;
            }
            // unpacked_node* E = arg2F->newUnpacked(e, FULL_ONLY);
            if(debug){
            E->show(meddlyout,true);
            printf("\n e^^^\n" );
            }
            dd_edge newstatesE(resF), djp(resF);

            for(unsigned j = 0; j < rSize; j++){
                if(B->d(j)==0) continue;
                node_handle newstates=0;
                compute_r(B->d(j),E->d(i),newstates);
                if (0==newstates) {
                  continue;
                }
                if (0 == D->d(j)) {
                  D->d_ref(j) = newstates;
                  continue;
                }
                // there's new states and existing states; union them.
                newstatesE.set(newstates);
                djp.set(D->d(j));
                mddUnion->computeTemp(newstatesE, djp, djp);
                D->set_d(j, djp);

            }
            // unpacked_node::recycle(E);

        }

        node_handle cnode = 0;
        cnode=resF->createReducedNode(int(i), D);
        C->d_ref(i) = cnode;

        unpacked_node::recycle(B);

    }
    unpacked_node::recycle(E);
    unpacked_node::recycle(A);


    res=resF->createReducedNode(-1,C);
    if(a!=-1&& e!=-1)
    saveResult(Key, a,e,res);
}


// ******************************************************************
// *                                                                *
// *                       covrCoveredToMT class                  *
// *                                                                *
// ******************************************************************
class MEDDLY::covrCoveredToMT:  public binary_operation {
public:
  covrCoveredToMT(binary_opname* oc, expert_forest* arg1,expert_forest* arg2, expert_forest* res);

  virtual void computeDDEdge(const dd_edge &arg, const dd_edge & arg2, dd_edge &res, bool userFlag);
  binary_operation* mddUnion;
   void compute_r( node_handle a,node_handle b, node_handle& c);

protected:
  inline ct_entry_key*
  findResult(node_handle a,node_handle b,node_handle &d)
  {
    ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeN(a);
    CTsrch->writeN(b);
    // CTsrch->writeI(c);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    d = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
  }
  inline node_handle saveResult(ct_entry_key* Key,
  node_handle a,node_handle b,node_handle d)
  {
    CTresult[0].reset();
    CTresult[0].writeN(d);
    CT0->addEntry(Key, CTresult[0]);
    return b;
  }
};
MEDDLY::covrCoveredToMT::covrCoveredToMT(binary_opname* oc, expert_forest* arg,expert_forest* arg2, expert_forest* res)
 : binary_operation(oc, 1, arg,arg2, res)
{
 MEDDLY_DCASSERT(resF->isForRelations());

 ct_entry_type* et = new ct_entry_type(oc->getName(), "NN:N");
 et->setForestForSlot(0, arg);
 et->setForestForSlot(1,arg);
 et->setForestForSlot(3, res);
 registerEntryType(0, et);
 buildCTs();
 dd_edge er(resF);
 mddUnion=getOperation(UNION, er,er,er);
 MEDDLY_DCASSERT(mddUnion);
}
void MEDDLY::covrCoveredToMT::computeDDEdge(const dd_edge &arg, const dd_edge & arg2, dd_edge &res, bool userFlag)
{

    // printf("covrCoveredToMT computeDDEdge\n" );
    node_handle bnode = 0;
    // RV above=E;
    // //////
    ostream_output meddlyout(std::cout);
    //
    compute_r( arg.getNode(),arg2.getNode(),  bnode);

    // printf("bnode is %d\n",bnode );
    // cij=_cij;
    res.set(bnode);
    bool debug=false;
    if(debug){
    arg.showGraph(meddlyout);
    printf("^^^arg1^^^\n" );
    arg2.showGraph(meddlyout);
    printf("arg2^^^^\n" );

    res.showGraph(meddlyout);
    printf("res\n" );
    printf("RESULT COVR coveredTo\n" );
    //
    getchar();
    }

}

void MEDDLY::covrCoveredToMT::compute_r( node_handle a,node_handle e, node_handle& res)
{
    bool debug=false;
    ostream_output meddlyout(std::cout);

    if(debug){
   printf("covrCoveredToMT a %d e %d\n",a,e );
    }
   if(a==0||e==0){
       res=0;
       return;
   }
   if(arg1F->isTerminalNode(a)&&arg2F->isTerminalNode(e)){
       if(debug){
       printf("terminal %d %d\n",a,e );
        }
       if(a==0||e==0) res=0;
       else{
            long aval;
            arg1F->getValueFromHandle(a, aval);
            long eval;
            arg2F->getValueFromHandle(e, eval);
            res=resF->handleForValue(aval*eval);
            if(debug){
            printf("res%ld\n",aval*eval );
            }
       }
       return;
   }
   //check the cache
   ct_entry_key* Key=findResult(a,e,res);
   if(0==Key){
       return;
   }

   const int aLevel=arg1F->getNodeLevel(a);
   const int eLevel=arg2F->getNodeLevel(e);
   const int rLevel=MAX(ABS(aLevel),ABS(eLevel));
   if(debug){
   printf("alevel %d, elevel %d\n",aLevel,eLevel );
    }

   const unsigned rSize= unsigned(resF->getLevelSize(rLevel));
    const unsigned ESize= unsigned(arg2F->getLevelSize(eLevel));
   unpacked_node* C=unpacked_node::newFull(resF,rLevel,rSize);

   unpacked_node* A = isLevelAbove(rLevel, aLevel)
    ? unpacked_node::newRedundant(arg1F, rLevel, a, true)
    : arg1F->newUnpacked(a, FULL_ONLY);
    unpacked_node* E = arg2F->newUnpacked(e, FULL_ONLY);
    if(debug){
    A->show(meddlyout,true);
    printf("\n a^^^\n" );
    }
    for (unsigned i = 0; i < rSize; i++) {
        // if(E->d(i)==0) continue;
        if(debug){
        printf("i %d a%d e %d\n",i,a, e );
        }
        int pLevel = arg1F->getNodeLevel(A->d(i));
        if(debug){
        printf("pLevel %d\n",pLevel );
        }
        unpacked_node* B = isLevelAbove(-rLevel, pLevel)
        ? unpacked_node::newIdentity(arg1F, -rLevel, i,  A->d(i), true)
        : arg1F->newUnpacked(A->d(i), FULL_ONLY);

        unpacked_node* D = unpacked_node::newFull(resF, -rLevel, rSize);
        if (rLevel > ABS(aLevel)) {
            if(debug){
            printf("IF\n" );
            }
            for(unsigned j = 0; j < rSize; j++){
                node_handle nres=0;
                // compute_r(B->d(j),e,nres);
                compute_r(a,E->d(j),nres);
                D->d_ref(j)=nres;
            }
        }else{
            if(debug){
            printf("ELSE\n" );
            }
            MEDDLY_DCASSERT(ABS(eLevel) >= ABS(pLevel));
            for (unsigned j = 0; j < rSize; j++) {
              D->d_ref(j) = 0;
            }
            // unpacked_node* E = arg2F->newUnpacked(e, FULL_ONLY);
            if(debug){
            E->show(meddlyout,true);
            printf("\n e^^^\n" );
            }
            dd_edge newstatesE(resF), djp(resF);

            for(unsigned j = 0; j < rSize; j++){
                if(B->d(j)==0) continue;
                node_handle newstates=0;
                compute_r(B->d(j),E->d(j),newstates);
                if (0==newstates) {
                  continue;
                }
                if (0 == D->d(j)) {
                  D->d_ref(j) = newstates;
                  continue;
                }
                // there's new states and existing states; union them.
                newstatesE.set(newstates);
                djp.set(D->d(j));
                mddUnion->computeTemp(newstatesE, djp, djp);
                D->set_d(j, djp);

            }
            // unpacked_node::recycle(E);

        }

        node_handle cnode = 0;
        cnode=resF->createReducedNode(int(i), D);
        C->d_ref(i) = cnode;

        unpacked_node::recycle(B);

    }
    unpacked_node::recycle(E);
    unpacked_node::recycle(A);


    res=resF->createReducedNode(-1,C);
    if(a!=-1&& e!=-1)
    saveResult(Key, a,e,res);
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
             // printf("covrExtractFrom_opname CORRECT!\n" );
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
    if (!arg->isForRelations() )
    {
        // printf("covrExtractCover_opname CORRECT!\n" );
        return new covrExtractCoverMT(this,arg,res);
    }

 };


 // ******************************************************************
 // *                                                                *
 // *                  covrCoveredFrom_opname class                  *
 // *                                                                *
 // ******************************************************************

 class MEDDLY::covrCoveredFrom_opname : public binary_opname {

 public:
 covrCoveredFrom_opname();
 virtual binary_operation* buildOperation(expert_forest* arg1,expert_forest* arg2, expert_forest* res);
 };

 MEDDLY::covrCoveredFrom_opname::covrCoveredFrom_opname()
         : binary_opname("Covered from")
 {
 }

 MEDDLY::binary_operation*
 MEDDLY::covrCoveredFrom_opname::buildOperation(expert_forest* arg1,expert_forest* arg2, expert_forest* res)
 {
    if (arg1->isForRelations()&& !arg2->isForRelations() && res->isForRelations())
    {
        // printf("covrCoveredFrom_opnameCORRECT!\n" );
        return new covrCoveredFromMT(this,arg1,arg2,res);
    }

 };

 // ******************************************************************
 // *                                                                *
 // *                  covrCoveredTo_opname class                  *
 // *                                                                *
 // ******************************************************************

 class MEDDLY::covrCoveredTo_opname : public binary_opname {

 public:
 covrCoveredTo_opname();
 virtual binary_operation* buildOperation(expert_forest* arg1,expert_forest* arg2, expert_forest* res);
 };

 MEDDLY::covrCoveredTo_opname::covrCoveredTo_opname()
         : binary_opname("Covered to")
 {
 }

 MEDDLY::binary_operation*
 MEDDLY::covrCoveredTo_opname::buildOperation(expert_forest* arg1,expert_forest* arg2, expert_forest* res)
 {
    if (arg1->isForRelations()&& !arg2->isForRelations() && res->isForRelations())
    {
        // printf("covrCoveredTo_opname CORRECT!\n" );
        // return nullptr;
        return new covrCoveredToMT(this,arg1,arg2,res);
    }

 };

 // ******************************************************************
 // *                                                                *
 // *                           Front  end                           *
 // *                                                                *
 // ******************************************************************

 MEDDLY::unary_opname* MEDDLY::initializeExtractFrom()
 {
         // printf("This is COVR EXTRACTFROM!\n" );
         return new covrExtractFrom_opname();
 }

 MEDDLY::unary_opname* MEDDLY::initializeExtractCovered(){
        // printf("This is Extract Covered!\n" );
        return new covrExtractCover_opname();
 }

 MEDDLY::binary_opname* MEDDLY::initializeCoveredFrom(){
        // printf("This is covered from\n" );
        // return nullptr;
        return new covrCoveredFrom_opname();
 }

 MEDDLY::binary_opname* MEDDLY::initializeCoveredTo(){
        // printf("This is covered to\n" );
        return new covrCoveredTo_opname();
 }

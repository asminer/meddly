
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


#include "mtmddbool.h"
#include "../unique_table.h"

MEDDLY::mt_mdd_bool::mt_mdd_bool(unsigned dsl, domain *d, const policies &p, int* level_reduction_rule, bool tv)
: mtmdd_forest(dsl, d, BOOLEAN, p, level_reduction_rule)
{
  initializeForest();

  transparent=bool_Tencoder::value2handle(tv);
}

MEDDLY::mt_mdd_bool::~mt_mdd_bool()
{ }

void MEDDLY::mt_mdd_bool::createEdge(bool term, dd_edge& e)
{
  createEdgeTempl<bool_Tencoder, bool>(term, e);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mdd_bool::createEdge(const int* const* vlist, int N, dd_edge &e)
{
  binary_operation* unionOp = getOperation(UNION, this, this, this);
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);

  int num_vars=getNumVariables();

  // Create vlist following the mapping between variable and level
  int** ordered_vlist=static_cast<int**>(malloc(N*sizeof(int*)+(num_vars+1)*N*sizeof(int)));
  if(ordered_vlist==0){
	  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }

  ordered_vlist[0]=reinterpret_cast<int*>(&ordered_vlist[N]);
  for(int i=1; i<N; i++) {
	  ordered_vlist[i]=(ordered_vlist[i-1]+num_vars+1);
  }
  for(int i=0; i<=num_vars; i++) {
	  int level=getLevelByVar(i);
	  for(int j=0; j<N; j++) {
		  ordered_vlist[j][level]=vlist[j][i];
	  }
  }

  mtmdd_edgemaker<bool_Tencoder, bool>
  EM(this, ordered_vlist, 0, order, N, getDomain()->getNumVariables(), unionOp);

  e.set(EM.createEdge());

  free(ordered_vlist);

#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mdd_bool::
createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a)
{
  createEdgeForVarTempl<bool_Tencoder, bool>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mdd_bool
::evaluate(const dd_edge &f, const int* vlist, bool &term) const
{
  term = bool_Tencoder::handle2value(evaluateRaw(f, vlist));
}

void MEDDLY::mt_mdd_bool::showTerminal(output &s, node_handle tnode) const
{
  bool_Tencoder::show(s, tnode);
}

void MEDDLY::mt_mdd_bool::writeTerminal(output &s, node_handle tnode) const
{
  bool_Tencoder::write(s, tnode);
}

MEDDLY::node_handle MEDDLY::mt_mdd_bool::readTerminal(input &s)
{
  return bool_Tencoder::read(s);
}

const char* MEDDLY::mt_mdd_bool::codeChars() const
{
  return "dd_tvb";
}

void MEDDLY::mt_mdd_bool::underApproximate(dd_edge &e, int Threashold)
{
    printf("underApproximate begin \n" );
    FILE_output meddlyout(stdout);

    showNode(meddlyout, e.getNode(), SHOW_DETAILS);
    // initialize the belowcount
    //
    int num_vars=getNumVariables();
    printf("num_vars: %d XXX \n",num_vars );
    int maxid=0;

    for(int ll=1; ll<(int)num_vars+1;ll++)
    {
        int sizeunique=unique->getNumEntries(ll);
        node_handle* uniqueNodes= new node_handle[sizeunique];
        unique->getItems(getVarByLevel(ll),uniqueNodes,sizeunique);
        for(int k=0;k<sizeunique;k++){
             // printf("%d\n",uniqueNodes[k]);
             if (uniqueNodes[k]>maxid)
             maxid=uniqueNodes[k];
        }


    }
    // printf("XXXXUNX\n" );
    belowcount=0;
    incomingedgecount=0;
    abovecount=0;
    highestunique.clear();
    uniquecount=0;
    // initialize arrays.
    // belowcount=new int[maxid];
    // incomingedgecount=new int[maxid];
    // abovecount=new int[maxid];
    // uniquecount= new int[maxid];

    printf("maxid %d\n",maxid );
    lastNode=maxid;
    // int count = sizeof(belowcount) / sizeof(*belowcount);
    // int count = *(&belowcount + 1) - belowcount;
    // printf("XXXXX SIZE = %d\n", MEDDLY::lastNode);
    operation* op = operation::getOpWithIndex(CARDINALITY->getIndex());
    removeOperationFromCache(op);
     op = operation::getOpWithIndex(IEC->getIndex());
     removeOperationFromCache(op);
     op = operation::getOpWithIndex(AC->getIndex());
     removeOperationFromCache(op);
     op = operation::getOpWithIndex(HU->getIndex());
     removeOperationFromCache(op);
     op = operation::getOpWithIndex(UC->getIndex());
     removeOperationFromCache(op);
    //op->removeAllComputeTableEntries();

    // printf("getCurrentNumNodes() %d\n", );
    int cC=unique->getNumEntries();//this->getCurrentNumNodes();
    while((int)cC>Threashold){
        printf("INSIDE CC>Threashold\n" );
        int diff=cC-Threashold;
        double cCard;
        apply(CARDINALITY,e,cCard);
        printf("CARDDDD %f \n", cCard );
        for(int i=0;i<(int)cC;i++){
           printf("BELOW %d \t %d\n", i, belowcount[i]);
        }
    for(int i=0;i<(int)cC;i++){
       printf("BELOW %d \t %d\n", i, belowcount[i]);
    }
    double cI;
    apply(IEC, e, cI);

    for(int i=0;i<(int)cC;i++){
 	  printf("IEC [%d] \t %d\n", i, incomingedgecount[i]);
    }
   // printf("Exactly %f\n", c1);
   double cA;
   apply(AC, e, cA);
   for(int i=0;i<(int)cC;i++){
     printf("AC [%d] \t %d\n", i, abovecount[i]);
   }
   double cH;
   apply(HU, e, cH);
   double cU;
   apply(UC, e, cU);

   for(int i=0;i<(int)cC;i++){
     printf("UC [ %d]= \t %d\n", i, uniquecount[i]);
   }
   int minIndex=0;
   double mindensity=DBL_MAX;
   bool findBasedOnDiff=false;
   for(int i=0;i<(int)cC;i++){
       if(uniquecount[i]==diff)
       {
          minIndex=i;
          printf("FOUND %d\n",i );
          findBasedOnDiff=true;
          // break;
       }
   }

  if(!findBasedOnDiff){
   double* density=new double[(int)cC];

   for(int i=0;i<(int)cC;i++){
       if(!((abovecount[i]==0)||(belowcount[i]==0)||(uniquecount[i]==0)))
     density[i]=(abovecount[i]*belowcount[i])/((double)uniquecount[i]);
     else
     density[i]=0;
     printf("%d %d %d %f\n", abovecount[i], belowcount[i],uniquecount[i], density[i]);
     if((density[i]!=0)&&(density[i]<mindensity))
     {
         minIndex=i;
         mindensity=density[i];
     }
   }


   printf("MinIndex: %d, mindensity: %f\n", minIndex+1,mindensity );
}
   printf("MinIndex is %d and getNodeLevel %d\n",minIndex+1,getNodeLevel(minIndex+1) );
   std::map<int,int> map;
   unsigned lvl=getNodeLevel(minIndex+1);
   int sizeunique=unique->getNumEntries(lvl);
   node_handle* uniqueNodes= new node_handle[sizeunique];
   // unique->getItems(uniqueNodes,sizeunique);
   printf("getSize() Unique : %d \n",sizeunique);
   FILE_output meddlyout(stdout);

   unique->show(meddlyout);
   unique->getItems(getVarByLevel(lvl),uniqueNodes,sizeunique);

  // unique->getItemsBylevel(uniqueNodes,1);

   for(int k=0;k<sizeunique;k++)
   {
       printf("%d\n",uniqueNodes[k]);
       map[uniqueNodes[k]]=uniqueNodes[k];
   }
   printf("MinIndex %d\n",(minIndex+1) );
   map[(minIndex+1)]=0;

   for(auto it = map.cbegin(); it != map.cend(); ++it)
{
     printf("map[%d]= %d\n",it->first ,it->second);
}

// expert_forest* ef = (expert_forest*) e.getForest();
// ef->deleteNode((minIndex+1));
// unlinkNode((minIndex+1));
// this->deleteNode((minIndex+1));
// uncacheNode((minIndex+1));
op = operation::getOpWithIndex(CARDINALITY->getIndex());
removeAllComputeTableEntries();
removeOperationFromCache(op);
showComputeTable(meddlyout,2);
printf("CT^^^\n" );
unsigned h = hashNode((minIndex+1));
int g=unique->remove(h,(minIndex+1));
// nodeMan->unlinkDownAndRecycle(getNodeAddress((minIndex+1)));
RemoveDuplicate2(lvl,map,e);
printf("%d\n",getNodeStatus((minIndex+1)));
cC=this->getCurrentNumNodes();
e.show(meddlyout,2);
printf("CC %d Threashold %d\n",cC,Threashold );
printf("unique_count %d\n",unique->getNumEntries() );
break;
}
// cleanup();
printf("this->getCurrentNumNodes() is %ld \n",this->getCurrentNumNodes() );

printf("underApproximate End\n" );
}
void MEDDLY::mt_mdd_bool::HeuristicUnderApproximate(dd_edge &e, int Threashold)
{
    printf("HunderApproximate\n" );

}
void MEDDLY::mt_mdd_bool::RemoveDuplicate2(int lvl, std::map<int,int> map,dd_edge &e){
    expert_forest* ef = (expert_forest*) e.getForest();
    bool changeInLevel=false;
    int num_vars=getNumVariables();
    FILE_output meddlyout(stdout);
    std::map<int,int> nmap;
    for(int l=0;l<=num_vars;l++){
    int sizeunique=unique->getNumEntries(l);
    node_handle* uniqueNodes= new node_handle[sizeunique];
    unique->getItems(getVarByLevel(l),uniqueNodes,sizeunique);
    for(int i=0;i<sizeunique;i++)
    showNode(meddlyout,uniqueNodes[i], SHOW_DETAILS);
}
    printf("DONE!\n" );
    // showNode(meddlyout,ef->getNode(), SHOW_DETAILS);
    printf("XXSTARTXX\n" );
    for(int l=lvl;l<num_vars;l++){
        changeInLevel=false;
        int sizeunique=unique->getNumEntries(lvl+1);
        printf("Number of entries  %d\n",sizeunique );
        node_handle* uniqueNodes= new node_handle[sizeunique];
        unique->getItems(getVarByLevel(lvl+1),uniqueNodes,sizeunique);
        int b=getDomain()->getVariableBound(lvl+1, false);
        for(int k=0;k<sizeunique;k++){
            unsigned h = hashNode(uniqueNodes[k]);
            unique->remove(h,uniqueNodes[k]);
            // printf("unique node is removed is %d\n",g );
            unpacked_node* un = unpacked_node::newFromNode(ef, uniqueNodes[k], unpacked_node::AS_STORED);
            // newRedundant(argF, level, a, false)
            showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);
            printf("\n");
            for(int ik=0;ik<=b; ik++){
                int dpt=getDownPtr(uniqueNodes[k],ik);
                if(dpt!=map[dpt]){
                    printf("Changed %d th ptr which was %d to %d\n", ik, dpt,map[dpt] );
                    un->d_ref(ik) = ef->linkNode(map[dpt]);
                }
            }
            node_handle q=unique->find(*un, getVarByLevel(un->getLevel()));
            printf("p is %d\n",q );
            if(q!=0){
                printf("CAME FIND\n" );
                changeInLevel=true;
                nmap[uniqueNodes[k]]=q;
            }
            else
            modifyReducedNodeInPlaceWithoutRecyclingDownPtr(un,uniqueNodes[k]);
            // ef->unlinkNode((node_handle)uniqueNodes[k]);
            // uniqueNodes[k]=ef->createReducedNode(-1, un);//result;
            showNode(meddlyout,uniqueNodes[k], SHOW_DETAILS);

            printf("Inside end RD\n" );
            // modifyReducedNodeInPlace(un, uniqueNodes[k]);
            // unpacked_node::recycle(un);

        }

    }
    // showNode(meddlyout,e.getNode(), SHOW_DETAILS);
    // printf("XXENDXX\n" );
    // e.show(meddlyout,2);
    for(int l=0;l<=num_vars;l++){
    int sizeunique=unique->getNumEntries(l);
    node_handle* uniqueNodes= new node_handle[sizeunique];
    unique->getItems(getVarByLevel(l),uniqueNodes,sizeunique);
    for(int i=0;i<sizeunique;i++)
    showNode(meddlyout,uniqueNodes[i], SHOW_DETAILS);
    printf("DONE!\n" );
}

}

void MEDDLY::mt_mdd_bool::RemoveDuplicate(int lvl, std::map<int,int> map){
    int num_vars=getNumVariables();
    FILE_output meddlyout(stdout);

    printf("num_vars %d\n",num_vars );
    bool changeInLevel=false;
    std::map<int,int> nmap;
    for(int l=lvl;l<=num_vars;l++){
        changeInLevel=false;
        int sizeunique=unique->getNumEntries(lvl+1);
        printf(" sizeuniquesizeunique inside RD %d\n",sizeunique);

        node_handle* uniqueNodes= new node_handle[sizeunique];
        unique->getItems(getVarByLevel(lvl+1),uniqueNodes,sizeunique);
        for(int k=0;k<sizeunique;k++)
        {
            int size=unique->getNumEntries(lvl+1);
            printf(" uniquesize inside RD %d\n",size);
            unsigned h = hashNode(uniqueNodes[k]);
            int g=unique->remove(h,uniqueNodes[k]);
            printf("g is %d\n",g );
             size=unique->getNumEntries(lvl+1);
            printf(" uniquesize inside RD %d\n",size);
            // unique->add(h,uniqueNodes[k]);
            // size=unique->getNumEntries(lvl+1);
           // printf(" uniquesize inside RD %d\n",size);
           int b=getDomain()->getVariableBound(lvl+1, false);
           printf("getVariableBound %d\n",b);
           unpacked_node* un = unpacked_node::newFromNode(this, uniqueNodes[k], unpacked_node::AS_STORED);
           showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);

           for(int ik=0;ik<=b; ik++){
               // printf("DPTR %d\n",getDownPtr(uniqueNodes[k],ik) );
               int dpt=getDownPtr(uniqueNodes[k],ik);
                if(dpt!=map[dpt])
                {

                    // dd_edge mapdpt(this);
                    // mapdpt.set(map[dpt]);
                    // un->set_d(ik,mapdpt);
                    // unpacked_node* nkp = unpacked_node::newFromNode(this, -k, unpacked_node::AS_STORED);
                    // nkp->i_ref(0) = pa;
                    // nkp->d_ref(0) = bottom;

                //     if(map[dpt]==0)
                //     {
                //         node_handle bottom = this->handleForValue(0);
                //
                //     un->d_ref(ik) = this->linkNode(bottom);//this->createReducedNode(ik, nkp);
                //
                // }
                //     else

                ////////////////////////////////////////////////////////////////////////
                    un->d_ref(ik) = linkNode(map[dpt]);//this->createReducedNode(ik, nkp);
                    printf("Changed %d th ptr which was %d to %d\n", ik, dpt,map[dpt] );
                ///////////////////////////////////////////////////////////////////

                    // un->show(meddlyout, true);
                    // printf("DONE SHOW IN LOOP\n" );
                    // node_handle node = createReducedNode(-1, un);
                    // // MEDDLY_DCASSERT(getNodeInCount(node) == 1 && getNodeLevel(node) == level + 1);
                    // uniqueNodes[k]=node;
                    // // swapNodes(uniqueNodes[ik], node);
                    // showNode(meddlyout, node, SHOW_DETAILS);
                    // showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);
                    //
                    // unlinkNode(node);
                    //
                    // showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);
                    //  printf("uniqueNodes DONE SHOW IN LOOP\n" );
                    // unpacked_node* AX = unpacked_node::newFromNode(this, uniqueNodes[k], false);

                }
                // un->show(meddlyout, true);
                // printf("DONE SHOW IN LOOP\n" );
           }
           // un->show(meddlyout, true);
           // printf("DONE SHOW IN LOOP\n" );


           // node_handle node = createReducedNode(-1, un);
           // // uniqueNodes[k]=node;
           // swapNodes(uniqueNodes[k], node);
           // showNode(meddlyout, node, SHOW_DETAILS);
           // showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);
           //
           // unlinkNode(node);


           uniqueNodes[k]=createReducedNode(-1, un);
           modifyReducedNodeInPlace(un, uniqueNodes[k]);

           // showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);
           //  printf("uniqueNodes DONE SHOW IN LOOP\n" );

           node_handle q=unique->find(*un, getVarByLevel(un->getLevel()));
           printf("p is %d\n",q );
           if(q!=0){
               printf("CAME FIND\n" );
               changeInLevel=true;
               nmap[uniqueNodes[k]]=q;
           }else{
               // node_handle node = createReducedNode(-1, un);
               // un->show(meddlyout, true);
               // printf("DONE SHOW IN LOOP\n" );
              // node_handle nx = createReducedNode(-1, un);
              printf("HERE ZEROOOO!!\n");
              showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);
               printf("HERE ZEROOOO!!\n");
            unsigned h = hashNode(uniqueNodes[k]);
               unique->add(h,uniqueNodes[k]);
               nmap[uniqueNodes[k]]=uniqueNodes[k];
           }


           // showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);
           printf("DONE!\n" );
           size=unique->getNumEntries(lvl+1);
            printf(" unique inside RD %d\n",size);
        }
        if(!changeInLevel)
        break;
        map=nmap;
    }
}

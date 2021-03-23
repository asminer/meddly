
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
#include <random>

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


void MEDDLY::mt_mdd_bool::underApproximate(dd_edge &e, long Threashold, long maxThreshold,float desiredPercentage,int option)
{
    if(option==1) return;
     clock_t start, end;
     start = clock();
    // printf("underApproximate begin \n" );
    FILE_output meddlyout(stdout);
    int num_vars=getNumVariables();
    // printf("num_vars %d\n",num_vars );
    int cC=e.getNodeCount();//unique->getNumEntries();//this->getCurrentNumNodes();
    int initialRootCardinality=e.getCardinality();
    double initialRootDensity=(initialRootCardinality/cC)*desiredPercentage;

    if(cC<maxThreshold) return;

    #if 0
    std::random_device rd;
std::default_random_engine generator;
generator.seed( rd() );
    std::normal_distribution<> d{(double)Threashold,((double)Threashold*0.05)};
    Threashold=std::round(d(generator));
    printf("Threashold %d\n",Threashold );

    #endif
    // Threashold=std::rand()%((int)(0.1*Threashold) + 1) +(int)(0.9*Threashold) ;
    /////////////*******
    // printf("NodeCount is %d\n",cC );
    // if(option==2){
    //     Threashold=initialRootDensity*desiredPercentage;
    // }
    while((cC>Threashold&&(option!=2))||(option==2)){
    //     printf("num_vars %d\n",num_vars );
    //     printf("cC%d\n",cC );
    //     printf( "Start cc>Threashold :");
    // printf("XXXXUNX\n" );
     maxid=e.getLastHandle();
    // printf("lastUsedHandle %d\n",maxid );
    //
    // printf("maxid %d\n",maxid );
    lastNode=maxid+1;
        // int diff=cC-Threashold;
        // printf("diff %d\n",diff );
        double cCard;
        apply(BC,e,cCard);
        // printf("CARDDDD %f \n", cCard );

    double cI;
    apply(IEC, e, cI);
    // printf("IEC Calculated\n" );

   double cA;
   apply(AC, e, cA);
   // printf("AC Calculated\n" );
   double cH;
   apply(HU, e, cH);
   // printf("HU Calculated\n" );
   double cU;
   apply(UC, e, cU);
   // printf("UC Calculated\n" );

   // for(int i=0;i<(int)maxid;i++){
   //     if((incomingedgecount[i]>0||i==e.getNode())&&(belowcount[i]<(long)1))
   //     {
   //         printf("ERROR i %d AC %ld BC %ld UC %d iEC %d \n", i,abovecount[i], belowcount[i],uniquecount[i]);
   //         char c=getchar();
   //      }
   // }
   int* levelcount=new int[num_vars+1];
   for(int i=0;i<=num_vars; i++){
       levelcount[i]=0;
   }
   for(int i=0;i<=(int)maxid;i++){
       if((abovecount[i]>0)||(i==e.getNode()))
        levelcount[getNodeLevel(i)]++;
   }

   int minIndex=0;
   if((option==0)||(option==2)){
   double mindensity=DBL_MAX;

      double density;
      mpz_object densitympz;
      densitympz.setValue(LONG_MAX);

      mpz_object option2comparingmpz;
      double intpart;
      double fractpart = std::modf (initialRootDensity , &intpart);
      option2comparingmpz.setValue(intpart);
      option2comparingmpz.setReminder(fractpart);

   for(int i=0;i<=(int)maxid;i++){
       if((uniquecount[i]>0)&&(levelcount[getNodeLevel(i)]>1))//&&(abovecount[i]>0)&&(belowcount[i]>0)&&(uniquecount[i]>0)&&(incomingedgecount[i]>0||i==e.getNode()))
       {
           mpz_object mulmpz;
           mulmpz.setValue(abovecount[i]);
           mulmpz.multiply(belowcount[i]);
           mulmpz.division(uniquecount[i]);
           if((mulmpz.compare(mulmpz, densitympz)<=0)&&((option!=2)||(option==2 &&mulmpz.compare(mulmpz, option2comparingmpz)<=0)))
           {
               minIndex=i;
               mulmpz.copyInto(densitympz);
               densitympz.setReminder(mulmpz.rdvalue);
           }
       }
   }

   // densitympz.showwithreminder(meddlyout);
   if(minIndex==0)
   {

    printf("ZERO %d\n",num_vars );
    // char c=getchar();
    // int k=0;
    //      for(int i=0;i<(int)maxid;i++){
    //
    //          printf("%d %d AC %ld BC %ld UC %d \n",k, i,abovecount[i], belowcount[i],uniquecount[i]);
    //          k++;
    //      }
    //      c=getchar();
    //      for(int i=0;i<=num_vars; i++){
    //          printf("level %d is %d\n", i,levelcount[i]);
    //      }
    //      c=getchar();
    //
    //      return;
   }
   delete[]levelcount;
   if(option==2 &&minIndex==0 )
    {
        end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        printf("Time taken %f \n",time_taken );
        return;
    }
    }
    if(option==1){
        node_handle n=e.getNode();
     node_handle* list = markNodesInSubgraph(&n, 1, false);
     long i;
     // for (i=0; list[i]; i++) { }
     // printf("i %d\n",i );
     for (i=0; list[i]; i++) {
         // printf("%d\n",list[i] );
    // if(uniquecount[list[i]]<1){printf("ERRROR %d\n",list[i] ); getchar();}
     if((uniquecount[list[i]]<=(cC-Threashold))&&(levelcount[getNodeLevel(list[i])]>1))
     {
         minIndex=list[i];
         // printf("Selected minIndex is %d and level is %d and UC is %d CC is %d \n",minIndex, getNodeLevel(minIndex), uniquecount[minIndex], cC );
         break;
     }
     // if(i==0){printf("%d %d\n",list[i],n );}
    }
    if(minIndex==0)
    printf("ERRROR \n" );
    // printf("Done\n" );
    // getchar();
     // printf("IIIIII %d\n",i );
     free(list);
    }
    // if(option==2){
    //
    // }
    unsigned lvl=getNodeLevel(minIndex);
   // printf("\nFINAL\n" );
    // printf("density MinIndex: %d, UC %d, BC %d AC %d\n", minIndex,uniquecount[minIndex] , belowcount[minIndex],abovecount[minIndex]);


   // printf("MinIndex is %d and getNodeLevel %d and IEC is %d\n",minIndex,lvl ,incomingedgecount[(minIndex)]);
   // if(lvl==1){
       // getchar();
       // for(int i=0;i<(int)maxid;i++){
       //     if(abovecount[i]>0&&getNodeLevel(i)==lvl)
       //     printf("%d AC %d BC %d UC %d lvl %d \n", i,abovecount[i], belowcount[i],uniquecount[i], getNodeLevel(i));
       // }
       // getchar();
   // }
   // c = getchar( );
   std::map<int,int> map;

   int sizeunique=unique->getNumEntries(lvl);
   node_handle* uniqueNodes= new node_handle[sizeunique];
   unique->getItems(getVarByLevel(lvl),uniqueNodes,sizeunique);

  // unique->getItemsBylevel(uniqueNodes,1);
   for(int k=0;k<sizeunique;k++)
   {
      if(incomingedgecount[uniqueNodes[k]]>0|| uniqueNodes[k]==e.getNode())
       map[uniqueNodes[k]]=uniqueNodes[k];
   }
   // printf("MinIndex %d\n",(minIndex+1) );
   map[(minIndex)]=0;
   delete uniqueNodes;
   // printf("MAke MAP! \n" );

// printf("getNodeStatus %d\n",getNodeStatus((minIndex)));
// int deleted=uniquecount[minIndex];
// printf("CALL RemoveDuplicate2 \n" );
RemoveDuplicate2(lvl,map,e);
cC=e.getNodeCount();
// printf("End CALL RemoveDuplicate2 \n" );

// cC--;
if(belowcount!=0)delete[] belowcount; else {printf("ERRR belowcount\n"); getchar();}
if(incomingedgecount!=0)delete [] incomingedgecount;else {printf("ERRR incomingedgecount\n"); getchar();}
if(abovecount!=0)delete[] abovecount;else {printf("ERRR abovecount\n"); getchar();}
highestunique.clear();
if(uniquecount!=0)delete[] uniquecount;else {printf("ERRR uniquecount\n"); getchar();}

// printf("DELETED \n" );


map.clear();
// e.show(meddlyout,0);
// printf("NodeCount is %d\n",cC );
}
    /////////////*******
// printf("NddeCountAfter %d\n",e.getNodeCount() );
end = clock();
double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
printf("Time taken %f \n",time_taken );
 // printf("underApproximate End\n" );
}
void MEDDLY::mt_mdd_bool::uniqueNodesforp(node_handle a/*,bool* visitedNode*/,std::set<int> &result ){
result.insert(a);
// if(visitedNode[a]){
//     return uniquecount[a];
// }
std::set<int> rset=highestunique[a];
for (auto it=rset.begin(); it != rset.end(); ++it){
   // printf("Res %d\n",*it );
   uniqueNodesforp((*it),result);
   // uniquecount[a]+=compute_r((*it));

}
}

void MEDDLY::mt_mdd_bool::HeuristicUnderApproximate(dd_edge &e, long Threashold, long maxThreshold,float desiredPercentage, int option)
{
    if(option==1) return;
    printf("HunderApproximate\n" );

         clock_t start, end;
         start = clock();
        FILE_output meddlyout(stdout);
        int num_vars=getNumVariables();
        int cC=e.getNodeCount();//unique->getNumEntries();//this->getCurrentNumNodes();
        if(cC<maxThreshold) {
            // printf("LESS THAN Threashold\n" );
            return;
        }
        int initialRootCardinality=e.getCardinality();
        double initialRootDensity=((double)initialRootCardinality/cC)*desiredPercentage;
        #if 0
        std::random_device rd;
    std::default_random_engine generator;
    generator.seed( rd() );
        std::normal_distribution<> d{(double)Threashold,((double)Threashold*0.05)};
        Threashold=std::round(d(generator));
        printf("Threashold %d\n",Threashold );

        #endif
        // while(cC>Threashold){
         maxid=e.getLastHandle();
        lastNode=maxid+1;

            double cCard;
            apply(BC,e,cCard);

        double cI;
        apply(IEC, e, cI);

       double cA;
       apply(AC, e, cA);
       double cH;
       apply(HU, e, cH);
       double cU;
       apply(UC, e, cU);

       int* levelcount=new int[num_vars+1];
       int* copylevelcount=new int[num_vars+1];
       for(int i=0;i<=num_vars; i++){
           levelcount[i]=0;
           copylevelcount[i]=0;
       }
       for(int i=0;i<=(int)maxid;i++){
           if((abovecount[i]>0)||(i==e.getNode()))
            {
                levelcount[getNodeLevel(i)]++;
                copylevelcount[getNodeLevel(i)]++;
            }
       }

       int minIndex=0;
       std::set<int> m;
       std::set<int> removedNode;
       std::set<int> neverdelete;
       if(option==0||(option==2)||(option==3)){
       // double mindensity=DBL_MAX;
       mpz_object* arrdensity=new mpz_object[maxid+1];
       for(int i=0;i<=maxid;i++)
       arrdensity[i].setValue(LONG_MAX);

          // double density;
          // mpz_object densitympz;
          // densitympz.setValue(LONG_MAX);

       for(int i=0;i<=(int)maxid;i++){
           if((uniquecount[i]>0)&&(levelcount[getNodeLevel(i)]>1))//&&(abovecount[i]>0)&&(belowcount[i]>0)&&(uniquecount[i]>0)&&(incomingedgecount[i]>0||i==e.getNode()))
           {
               // mpz_object mulmpz;
               // mulmpz.setValue(abovecount[i]);
               // mulmpz.multiply(belowcount[i]);
               // mulmpz.division(uniquecount[i]);
               arrdensity[i].setValue(abovecount[i]);
               arrdensity[i].multiply(belowcount[i]);
               arrdensity[i].division(uniquecount[i]);
               // if(mulmpz.compare(mulmpz, densitympz)<=0)
               // {
               //     minIndex=i;
               //     mulmpz.copyInto(densitympz);
               //     densitympz.setReminder(mulmpz.rdvalue);
               // }
           }
       }

       int ck=cC;
       int root=e.getNode();
       mpz_object option2comparingmpz;
       double intpart;
       double fractpart = std::modf (initialRootDensity , &intpart);
       option2comparingmpz.setValue(intpart);
       option2comparingmpz.setReminder(fractpart);
       bool option3densitycheck=false;
       removedNode.clear();
       // printf("**IN WHILE LOOP option %d\n", option );
        while(((cC-removedNode.size()>Threashold)&&((option==0)||(option==3)))||(option==2)){
            // printf("**IN WHILE LOOP\n" );
            mpz_object densitympz;
            densitympz.setValue(LONG_MAX);
            minIndex=0;
            for(int i=0;i<=(int)maxid;i++){

                // arrdensity[i].showwithreminder(meddlyout);
                // printf("\n");
                // option2comparingmpz.showwithreminder(meddlyout);
                // printf("\n");
                // printf("______\n" );
                if(((incomingedgecount[i]>0)&&(arrdensity[i].compare(arrdensity[i], densitympz)<=0)&&(copylevelcount[getNodeLevel(i)]>1)&&((option==0)||(option==3 && option3densitycheck))||
                ((incomingedgecount[i]>0)&&(arrdensity[i].compare(arrdensity[i], option2comparingmpz)<=0)&&(copylevelcount[getNodeLevel(i)]>1)&&((option==2)||(option==3 && !option3densitycheck)))))
                //&&(!(removedNode.find(i))))
                {
                    // printf("option3densitycheck %d size %d\n",option3densitycheck,removedNode.size());
                     printf("SELECTED\n" );
                    bool is_in = removedNode.find(i) != removedNode.end();
                    bool neverDelete_isin = neverdelete.find(i) != neverdelete.end();
                    if(!is_in && !neverDelete_isin){
                    minIndex=i;
                    arrdensity[i].copyInto(densitympz);
                    densitympz.setReminder(arrdensity[i].rdvalue);
                    }
                }
            }
            if(minIndex!=0){

            // bool* visitedNode= new bool[maxid+1];
            // for(int i=0;i<=maxid;i++){
            //     visitedNode[i]=false;
            // }
              std::set<int> resultp;
            uniqueNodesforp(minIndex/*,visitedNode*/,resultp);
            if(resultp.size()!=uniquecount[minIndex])
            {printf("ERRR in uniqueNodesforp\n");getchar();}
            bool shouldadd=true;
            for(auto i= resultp.begin();i!=resultp.end();++i){
                if(copylevelcount[getNodeLevel(*i)]<=1)
                {
                    shouldadd=false;
                    printf("should not add\n" );
                    //getchar();
                    neverdelete.insert(minIndex);
                    break;
                }
            }
            if(shouldadd){

                m.insert(minIndex);
            // printf("%d %d %d\n",minIndex,uniquecount[minIndex], resultp.size() );
            // printf("ADDED %d\n", option3densitycheck );
            // printf("ADDED\n" );
            // getchar();
            removedNode.insert(resultp.begin(), resultp.end());
            for(int i=0;i<=num_vars; i++)
            copylevelcount[i]=levelcount[i];
            for(auto i= removedNode.begin();i!=removedNode.end();++i){
                copylevelcount[getNodeLevel(*i)]--;
            }
            // for(auto i = removedNode.begin(); i != removedNode.end(); ++i)
            // printf("%d, ",(*i) );
            // cC=ck-removedNode.size();
            // printf("\n%d %d %d\n",ck,removedNode.size(), cC );
            }
            else{
                neverdelete.insert(minIndex);
            }
            // getchar();
            }else{
                // printf("CAME ELSE\n" );
                // getchar();
                if(option==2 )
                 {
                     printf("Deleted node by density %d\n",removedNode.size() );
                     // end = clock();
                     // double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
                     // printf("Time taken %f \n",time_taken );
                     break;
                 }
                if(option==3 && !option3densitycheck){
                    option3densitycheck=true;
                    printf("Deleted node by density %d\n",removedNode.size() );
                    if(removedNode.size()==0)
                    { getchar();

                    }
                    continue;
                }else{
                // printf("Cannot Delete more\n" );
                // getchar();
                // for(int i=0;i<num_vars;i++)
                // printf("%d %d\n",i, levelcount[i] );
                // getchar();
                break;
                // printf("Error\n" );
                }
            }
            // for(int i=0;i<=num_vars; i++)
            // copylevelcount[i]=levelcount[i];
        }
        printf("Expected %d\n", cC-removedNode.size());
        // printf("rn %d\n",removedNode.size() );

        // printf("cC%d rn %d\n",cC, removedNode.size() );
        // for(auto i = removedNode.begin(); i != removedNode.end(); ++i)
        // printf("%d, ",(*i) );
        // getchar();

       // densitympz.showwithreminder(meddlyout);
       // if(minIndex==0)
       // {
       //  printf("ZERO %d\n",num_vars );
       //  char c=getchar();
       //  int k=0;
       //       for(int i=0;i<(int)maxid;i++){
       //
       //           printf("%d %d AC %ld BC %ld UC %d \n",k, i,abovecount[i], belowcount[i],uniquecount[i]);
       //           k++;
       //       }
       //       c=getchar();
       //       for(int i=0;i<=num_vars; i++){
       //           printf("level %d is %d\n", i,levelcount[i]);
       //       }
       //       c=getchar();
       //
       //       return;
       // }
       delete[]levelcount;
       delete[]copylevelcount;
       delete[] arrdensity;

        }
        if(option==1){
        //     node_handle n=e.getNode();
        //  node_handle* list = markNodesInSubgraph(&n, 1, false);
        //  long i;
        //  for (i=0; list[i]; i++) {
        //  if((uniquecount[list[i]]<=(cC-Threashold))&&(levelcount[getNodeLevel(list[i])]>1))
        //  {
        //      minIndex=list[i];
        //      break;
        //  }
        // }
        // if(minIndex==0)
        // printf("ERRROR \n" );
        //
        //  free(list);
        }


       std::map<int,int> map;
       std::set<int> levels;
       int minlvl=INT_MAX;
       node_handle n=e.getNode();
       node_handle* list = markNodesInSubgraph(&n, 1, false);
       for (long i=0; list[i]; i++) {
            map[list[i]]=list[i];
        }
       delete[] list;

       // std::set<int> tempres;
       // for(auto i = m.begin(); i != m.end(); ++i){
       //     std::set<int>resultp;
       //     uniqueNodesforp((*i)/*,visitedNode*/,resultp);
       //     tempres.insert(resultp.begin(), resultp.end());
       // }
       // printf("%d=====%d\n",tempres.size(), removedNode.size() );
       for(auto i = m.begin(); i != m.end(); ++i){
          map[(*i)]=0;
          int nodelvl=getNodeLevel((*i));
          levels.insert(nodelvl);
          if(nodelvl<minlvl)
            minlvl=nodelvl;
        }

       // for(int k=0;k<sizeunique;k++)
       // {
       //    if(incomingedgecount[uniqueNodes[k]]>0|| uniqueNodes[k]==e.getNode())
       //     map[uniqueNodes[k]]=uniqueNodes[k];
       // }


       // unsigned lvl=getNodeLevel(minIndex);
       // int sizeunique=unique->getNumEntries(lvl);
       // node_handle* uniqueNodes= new node_handle[sizeunique];
       // unique->getItems(getVarByLevel(lvl),uniqueNodes,sizeunique);
       //
       // for(int k=0;k<sizeunique;k++)
       // {
       //    if(incomingedgecount[uniqueNodes[k]]>0|| uniqueNodes[k]==e.getNode())
       //     map[uniqueNodes[k]]=uniqueNodes[k];
       // }
       // // printf("MinIndex %d\n",(minIndex+1) );
       // map[(minIndex)]=0;
       // delete uniqueNodes;


       // printf("MAke MAP! \n" );

    // printf("getNodeStatus %d\n",getNodeStatus((minIndex)));
    // int deleted=uniquecount[minIndex];
    // printf("CALL RemoveDuplicate2 \n" );
     // for(auto i = levels.begin(); i != levels.end(); ++i){
     //     printf("%d , \n",(*i) );
     // }
     // getchar();
    ////////////////********************

    // for(auto i = m.begin(); i != m.end(); ++i){
    //     printf("%d ,",(*i) );
    // }
    // printf("\n" );
    RemoveDuplicateSet(minlvl,levels,map,e);
    cC=e.getNodeCount();
    // printf("End CALL RemoveDuplicate2 \n" );

    // cC--;
    if(belowcount!=0)delete[] belowcount; else {printf("ERRR belowcount\n"); getchar();}
    if(incomingedgecount!=0)delete [] incomingedgecount;else {printf("ERRR incomingedgecount\n"); getchar();}
    if(abovecount!=0)delete[] abovecount;else {printf("ERRR abovecount\n"); getchar();}
    highestunique.clear();
    if(uniquecount!=0)delete[] uniquecount;else {printf("ERRR uniquecount\n"); getchar();}

    // printf("DELETED \n" );


    map.clear();
    // e.show(meddlyout,0);
    // printf("HNodeCount is %d\n",cC );
    // getchar();

    // }
        /////////////*******
    // printf("NddeCountAfter %d\n",e.getNodeCount() );
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    printf("Time taken %f \n",time_taken );
    printf("HunderApproximate End\n" );
}
 int MEDDLY::mt_mdd_bool::RemoveDuplicateSet(int lvl, std::set<int>levels,std::map<int,int>map,dd_edge &e){
     int root=e.getNode();
      // printf("size is %d\n",e.getNodeCount() );
     bool changeInLevel=false;
     int num_vars=getNumVariables();
     FILE_output meddlyout(stdout);
     // std::map<int,int> nmap;
     // printf("RemoveDuplicateSet\n" );
     int result=0;
     node_handle* uniqueNodes;
     for(int l=lvl;l<num_vars;l++){
         // printf("l is %d\n",l );
         // getchar();
         bool is_in = levels.find(l) != levels.end();
         // printf("is_in is %d\n",is_in );
         // getchar();
         if(is_in){
             levels.erase(l);
             // for(auto i = levels.begin(); i != levels.end(); ++i){
             //     printf("%d , \n",(*i) );
             // }
         }
         changeInLevel=false;

         int sizeunique=unique->getNumEntries(l+1);
         uniqueNodes= new node_handle[sizeunique];
         unique->getItems(getVarByLevel(l+1),uniqueNodes,sizeunique);
         int b=getDomain()->getVariableBound(l+1, false);
         for(int k=0;k<sizeunique;k++){
         if(((incomingedgecount[uniqueNodes[k]]>0)||(uniqueNodes[k]==root))&&(uniqueNodes[k]<=maxid)&&(map[uniqueNodes[k]]!=0)){
             unpacked_node* un =  unpacked_node::newRedundant(this, l+1, uniqueNodes[k], true);
             bool nodeChanged=false;
             for(int ik=0;ik<=b; ik++){
                 int dpt=getDownPtr(uniqueNodes[k],ik);
                 if(dpt>=1 && dpt!=map[dpt]){
                     un->d_ref(ik) = this->linkNode(map[dpt]);
                     nodeChanged=true;
                     // printf("\nChanged %d th ptr which was %d to %d and lvl is %d \n", ik, dpt,map[dpt], getNodeLevel(dpt) );
                     // getchar();
                 }
                 else{
                     un->d_ref(ik)=this->linkNode(dpt);
                 }
             }
             if(nodeChanged){
                 // showNode(meddlyout,uniqueNodes[k], SHOW_DETAILS);
                 changeInLevel=true;
                 node_handle q=unique->find(*un, getVarByLevel(un->getLevel()));
                 if(q!=0){
                     result++;
                     // printf("Found\n" );
                     map[uniqueNodes[k]]=q;
                     unpacked_node::recycle(un);
                    }else{
                        if(l+1<num_vars){
                            node_handle temporaryNode=this->createReducedNode(-1,un);
                            this->linkNode(temporaryNode);
                              // showNode(meddlyout, temporaryNode, SHOW_DETAILS);
                             map[uniqueNodes[k]]=temporaryNode;
                              // printf("DONE MODIFY IN PLACE %d to %d\n",uniqueNodes[k],temporaryNode );
                              // getchar();

                         }else{
                             // printf("HERE WANT TO set root\n");
                             node_handle node=this->createReducedNode(-1,un);
                             node = makeNodeAtTop(node);
                             this->linkNode(node);
                             e.set(node);
                             // showNode(meddlyout, node, SHOW_DETAILS);
                             // printf("Root set!\n");
                             // printf("end size is %d\n",e.getNodeCount() );
                             //nmap.clear();
                             map.clear();
                             // changeInLevel=false;
                             return result;
                         }
                     }
                 }else{
                 // map[uniqueNodes[k]]=uniqueNodes[k];
                 unpacked_node::recycle(un);
             }
         }
         else{
             ;//nothing
         }
        }
        delete[] uniqueNodes;
        if(changeInLevel==false){
            if(levels.size()==0)
            {
                printf("NOT CORRECT\n" );
                getchar();
                break;
            }
        }
        else{
            ;
            // if(levels.size()>0)
            // l=*(levels.begin())-1;
            // printf("new l is\n", );
        }
     }
return 0;
printf("Done RemoveDuplicateSet\n" );

 }
int MEDDLY::mt_mdd_bool::RemoveDuplicate2(int lvl, std::map<int,int> map,dd_edge &e){
    int root=e.getNode();
    // printf("ROOT is \n",root );
    int result=0;
    bool changeInLevel=false;
    int num_vars=getNumVariables();
    FILE_output meddlyout(stdout);
    std::map<int,int> nmap;
    // printf("num_vars %d\n",num_vars );
    node_handle* uniqueNodes;
    // for(int i=0;i<(int)maxid;i++){
    //      printf(" i %d AC %ld BC %ld UC %d iEC %d \n", i,abovecount[i], belowcount[i],uniquecount[i]);
    //     // if((incomingedgecount[i]>0||i==e.getNode())&&(belowcount[i]<(long)1))
    //     // {
    //     //     printf("ERROR i %d AC %ld BC %ld UC %d iEC %d \n", i,abovecount[i], belowcount[i],uniquecount[i]);
    //     //     char c=getchar();
    //     //  }
    // }
    for(int l=lvl;l<num_vars;l++){
        // nmap.clear();
        // printf("lvl %d\n",l );
        changeInLevel=false;
        int sizeunique=unique->getNumEntries(l+1);
         // printf("Number of entries  %d\n",sizeunique );
        uniqueNodes= new node_handle[sizeunique];
        unique->getItems(getVarByLevel(l+1),uniqueNodes,sizeunique);
        // printf("get unique items %d\n",sizeunique );
        int b=getDomain()->getVariableBound(l+1, false);
        // printf("domain is %d\n",b );
        for(int k=0;k<sizeunique;k++)
        {
            // nmap[uniqueNodes[k]]=uniqueNodes[k];
        if(((incomingedgecount[uniqueNodes[k]]>0)||(uniqueNodes[k]==root))&&(uniqueNodes[k]<=maxid))//abovecount[uniqueNodes[k]]>0)//
        {
            // printf("Node %d %d \n",uniqueNodes[k], maxid );
            nmap[uniqueNodes[k]]=uniqueNodes[k];
            unpacked_node* un =  unpacked_node::newRedundant(this, l+1, uniqueNodes[k], true);//unpacked_node::newFull(this, l+1, b);
            bool nodeChanged=false;
            for(int ik=0;ik<=b; ik++){
                int dpt=getDownPtr(uniqueNodes[k],ik);
                if(dpt>=1 && dpt!=map[dpt]){
                    // printf("\nChanged %d th ptr which was %d to %d and lvl is %d \n", ik, dpt,map[dpt], getNodeLevel(dpt) );
                    un->d_ref(ik) = this->linkNode(map[dpt]);
                    nodeChanged=true;
                }
                else{
                    un->d_ref(ik)=this->linkNode(dpt);
                }
            }
            if(nodeChanged){
                // showNode(meddlyout,uniqueNodes[k], SHOW_DETAILS);
                changeInLevel=true;
            node_handle q=unique->find(*un, getVarByLevel(un->getLevel()));
             // printf("p is %d\n",q );
             if(q!=0){
                result++;
                nmap[uniqueNodes[k]]=q;
                unpacked_node::recycle(un);
                }
             else{
                 if(l+1<num_vars){
                node_handle temporaryNode=this->createReducedNode(-1,un);
                this->linkNode(temporaryNode);
                // showNode(meddlyout, temporaryNode, SHOW_DETAILS);
                nmap[uniqueNodes[k]]=temporaryNode;
                // printf("DONE MODIFY IN PLACE %d to %d\n",uniqueNodes[k],temporaryNode );

                }
                else
                {
                    // printf("HERE WANT TO set root\n");
                    node_handle node=this->createReducedNode(-1,un);
                    node = makeNodeAtTop(node);
                    this->linkNode(node);
                    e.set(node);
                    // printf("Root set!\n");
                    nmap.clear();
                    map.clear();
                    // changeInLevel=false;
                    return result;
                }
                }

            }else{
                nmap[uniqueNodes[k]]=uniqueNodes[k];
                unpacked_node::recycle(un);
            }
        }else{;
            // printf("Node %d has iec %d\n",uniqueNodes[k],incomingedgecount[uniqueNodes[k]] );
        }

    }
    delete[] uniqueNodes;
    map.clear();

    // printf("CHANGE IN LVL is %d\n",changeInLevel );
    if(changeInLevel==false){
        // printf("BREAK\n" );
        break;
    }
    else{
        std::map<int, int>::iterator itr;
        map.clear();
        for (itr = nmap.begin(); itr != nmap.end(); ++itr) {
        map[ itr->first ]= itr->second ;
        }
    }
    nmap.clear();
    }
    // printf("DONE! before check\n" );
    // printf("ROOT IS %d\n",e.getNode() );
    //     printf("DONE!\n" );


}
/*
// showNode(meddlyout,uniqueNodes[k], SHOW_DETAILS);
// printf("\n" );
// c=getchar();
//map=nmap;
// printf("MAp %d %d \n",itr->first,itr->second );
// ef->unlinkNode(uniqueNodes[k]);
//unpacked_node::recycle(un);
// char c=getchar();
// e.setForest(this);
// this->registerEdge(e);
// e.set(node);
// unlinkNode(root);
// this->unlinkNode(prevroot);
printf("CCCC %d\n",e.getNodeCount() );

*/
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
           // printf("p is %d\n",q );
           if(q!=0){
               // printf("CAME FIND\n" );
               changeInLevel=true;
               nmap[uniqueNodes[k]]=q;
           }else{
               // node_handle node = createReducedNode(-1, un);
               // un->show(meddlyout, true);
               // printf("DONE SHOW IN LOOP\n" );
              // node_handle nx = createReducedNode(-1, un);
              // printf("HERE ZEROOOO!!\n");
              // showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);
               // printf("HERE ZEROOOO!!\n");
            unsigned h = hashNode(uniqueNodes[k]);
               unique->add(h,uniqueNodes[k]);
               nmap[uniqueNodes[k]]=uniqueNodes[k];
           }


           // showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);
           // printf("DONE!\n" );
           size=unique->getNumEntries(lvl+1);
            // printf(" unique inside RD %d\n",size);
        }
        if(!changeInLevel)
        break;
        map=nmap;
    }
}

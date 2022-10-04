
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
#include <algorithm>
#include <climits>
#ifdef HAVE_LIBGMP
#include <gmp.h>
#endif

#define PRINTON

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


    srand(time(0));
printf("UAStart %d\n",option );
if(option==1) return;
srand(time(0));
clock_t start, end;
start = clock();
// printf("underApproximate begin \n" );
FILE_output meddlyout(stdout);
int num_vars=getNumVariables();
// printf("num_vars %d\n",num_vars );
int cC=e.getNodeCount();                      //unique->getNumEntries();//this->getCurrentNumNodes();
// printf("cC %d\n",cC );
if(cC<maxThreshold) return;
int initialRootCardinality=e.getCardinality();
double initialRootDensity=(initialRootCardinality/cC)*desiredPercentage;



// #if 0
// std::random_device rd;
// std::default_random_engine generator;
// generator.seed( rd() );
// std::normal_distribution<> d{(double)Threashold,((double)Threashold*0.05)};
// Threashold=std::round(d(generator));
// printf("Threashold %d\n",Threashold );
//
// #endif
// Threashold=std::rand()%((int)(0.1*Threashold) + 1) +(int)(0.9*Threashold) ;
/////////////*******
// printf("NodeCount is %d\n",cC );
// if(option==2){
//     Threashold=initialRootDensity*desiredPercentage;
// }
int num_deletedNode=0;
while((cC>Threashold&&(option!=2))||(option==2)) {
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
#ifdef HAVE_LIBGMP
        mpz_object mpzd;
        mpz_object mpz0;
        mpzd.setValue(0);
        mpzd.setReminder(0);
        mpz0.setValue(0);
        mpz0.setReminder(0);
        ACBC=new mpz_object[lastNode];
#else
        ACBC=new long double[lastNode];
#endif
        // printf("ACBC created\n" );
        node_handle nl=e.getNode();

        node_handle* list = markNodesInSubgraph(&nl, 1, false);
        std::list<int>* uniquelist=markNodesInSubgraphByLvl(&nl,1);
        for(long i=0; i<lastNode; i++) {
#ifdef HAVE_LIBGMP
                ACBC[i].setValue(0);
                ACBC[i].setReminder(0);
#else
                ACBC[i]=0;
#endif
        }
        for (int i=0; list[i]; i++) {
                // for(int i=0;i<lastNode;i++){
                // if(abovecount[i]>1)
#ifdef HAVE_LIBGMP
                ACBC[list[i]].setValue(abovecount[list[i]]);
                ACBC[list[i]].multiply(belowcount[list[i]]);
#else
                ACBC[list[i]]=(long double) abovecount[list[i]]*belowcount[list[i]];

                if (ACBC[list[i]]/abovecount[list[i]]!=belowcount[list[i]])
                { printf("NOT CORRECT%llu,%lu, %lu\n",ACBC[i],abovecount[list[i]],belowcount[list[i]] );
                  exit(0);
                  // getchar();
                }
#endif
                // else
                // ACBC[i]=0;
        }
        delete[] abovecount;
        // printf("abovecount deleted\n" );
        delete[] belowcount;
        // printf("belowcount deleted\n" );
        double cH;
        apply(HU, e, cH);
        // printf("HU Calculated\n" );
        double cU;
        apply(UC, e, cU);

        double cL;
        apply(LU, e, cL);

        double cUA;
        apply(UAC, e, cUA);
        // double caU;
        // apply(UAC, e, caU);
        // printf("UC Calculated\n" );

        // for(int i=0;i<(int)maxid;i++){
        //     if((incomingedgecount[i]>0||i==e.getNode())&&(belowcount[i]<(long)1))
        //     {
        //         printf("ERROR i %d AC %ld BC %ld UC %d iEC %d \n", i,abovecount[i], belowcount[i],uniquecount[i]);
        //         char c=getchar();
        //      }
        // }
        int* levelcount=new int[num_vars+1];
        for(int i=0; i<=num_vars; i++) {
                levelcount[i]=uniquelist[i].size();         //0;
        }
        // for(int i=0;i<=(int)maxid;i++){
        //     if((abovecount[i]>0)||(i==e.getNode()))
        //      levelcount[getNodeLevel(i)]++;
        // }

        int minIndex=0;
        if((option==0)||(option==2)) {
#ifdef HAVE_LIBGMP
mpz_object mindensity;
mindensity.setValue(INT_MAX);
mindensity.setReminder(0);
mpz_object density;
density.setValue(INT_MAX);
density.setReminder(0);
#else
                double mindensity=DBL_MAX;

                double density=DBL_MAX;
#endif
                // mpz_object densitympz;
                std::list<int> resultSet;
                // mpz_object resultsetValue;
                // densitympz.setValue(LONG_MAX);
                // resultsetValue.setValue(LONG_MAX);
                // mpz_object option2comparingmpz;
                // double intpart;
                // double fractpart = std::modf (initialRootDensity , &intpart);
                // option2comparingmpz.setValue(intpart);
                // option2comparingmpz.setReminder(fractpart);
                // doubleDensityClass* doubleDensityArray=new doubleDensityClass[lastNode];
                //
                // for (i=0; list[i]; i++) {
                //     doubleDensityArray[list[i]].density=(long double)ACBC[list[i]]/(double)(uniquecount[list[i]]+uniqueAbovecount[list[i]]);//((double)(abovecount[list[i]])/(double)(uniquecount[list[i]]+uniqueAbovecount[list[i]]))*(double)belowcount[list[i]];
                //             doubleDensityArray[list[i]].index=list[i];
                //             doubleDensityArray[list[i]].isIn=false;
                //             doubleDensityArray[list[i]].neverIn=false;
                // }
                // std::sort(doubleDensityArray, doubleDensityArray+lastNode);

                for(int i=0; i<=(int)maxid; i++) {
#ifdef HAVE_LIBGMP
              if((mpzd.compare(ACBC[i],mpz0)==1)&&(levelcount[getNodeLevel(i)]>1))
#else
                        if((ACBC[i]>0)&&(levelcount[getNodeLevel(i)]>1))     //&&(abovecount[i]>0)&&(belowcount[i]>0)&&(uniquecount[i]>0)&&(incomingedgecount[i]>0||i==e.getNode()))
#endif
                        {
#ifdef HAVE_LIBGMP

                                mpz_object densityi;
                                densityi=ACBC[i];
                                densityi.setReminder(0);
                                densityi.division((uniquecount[i]+uniqueAbovecount[i]));
                                // mpz_object mulmpz;
                                // mulmpz.setValue(abovecount[i]);
                                // mulmpz.multiply(belowcount[i]);
                                // mulmpz.division(uniquecount[i]+uniqueAbovecount[i]);
                                if((mpzd.compare(densityi,density)<=0 /*mulmpz.compare(mulmpz, densitympz)<=0*/)&&((option!=2)||(option==2 /*&&mulmpz.compare(mulmpz, option2comparingmpz)<=0*/)))

#else
                                double densityi=ACBC[i]/(uniquecount[i]+uniqueAbovecount[i]);
                                // mpz_object mulmpz;
                                // mulmpz.setValue(abovecount[i]);
                                // mulmpz.multiply(belowcount[i]);
                                // mulmpz.division(uniquecount[i]+uniqueAbovecount[i]);
                                if((densityi<=density /*mulmpz.compare(mulmpz, densitympz)<=0*/)&&((option!=2)||(option==2 /*&&mulmpz.compare(mulmpz, option2comparingmpz)<=0*/)))
#endif
                                {
                                        if(resultSet.size()>0) {
                                                if(option==0) {
#ifdef HAVE_LIBGMP
                if(mpzd.compare(density,densityi)==0 /*resultsetValue.compare(resultsetValue,mulmpz)==0*/)

#else
                                                        if(density==densityi /*resultsetValue.compare(resultsetValue,mulmpz)==0*/)
#endif
                                                        {
                                                                resultSet.push_back(i);
                                                                // mulmpz.copyInto(resultsetValue);
                                                                // resultsetValue.setReminder(mulmpz.rdvalue);
                                                                // density=densityi
                                                        }

                                                        /* else{
                                                             resultSet.clear();
                                                             resultSet.insert(i);
                                                              mulmpz.copyInto(densitympz);
                                                              densitympz.setReminder(mulmpz.rdvalue);
                                                             mulmpz.copyInto(resultsetValue);
                                                             resultsetValue.setReminder(mulmpz.rdvalue);
                                                           }
                                                           }else if(option==2){
                                                           if(resultsetValue.compare(resultsetValue,mulmpz)==0){
                                                             resultSet.insert(i);
                                                             mulmpz.copyInto(resultsetValue);
                                                             resultsetValue.setReminder(mulmpz.rdvalue);
                                                           }else{
                                                             resultSet.clear();
                                                             mulmpz.copyInto(densitympz);
                                                             densitympz.setReminder(mulmpz.rdvalue);
                                                             resultSet.insert(i);
                                                             mulmpz.copyInto(resultsetValue);
                                                             resultsetValue.setReminder(mulmpz.rdvalue);
                                                           }*/
                                                }
                                        }else{
                                                resultSet.push_back(i);
                                                density=densityi;
                                                // mulmpz.copyInto(resultsetValue);
                                                // mulmpz.copyInto(densitympz);
                                                // densitympz.setReminder(mulmpz.rdvalue);
                                                // resultsetValue.setReminder(mulmpz.rdvalue);
                                        }
                                        // minIndex=i;
                                        // mulmpz.copyInto(densitympz);
                                        // densitympz.setReminder(mulmpz.rdvalue);
                                }
                        }
                }
                if(resultSet.size()>0) {
                        // printf("Density %lf\n", density);
                        std::list<int>::iterator iter = resultSet.begin();
                        int dgen=rand() % (resultSet.size());
                        std::advance(iter, dgen);
                        minIndex=(*(iter));
                        resultSet.clear();
                }
                // printf("minIndex %d\n",minIndex );
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
        if(option==1) {
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
        for(int k=0; k<sizeunique; k++)
        {
                if(incomingedgecount[uniqueNodes[k]]>0|| uniqueNodes[k]==e.getNode())
                        map[uniqueNodes[k]]=uniqueNodes[k];
        }
        // printf("MinIndex %d\n",(minIndex+1) );
        map[(minIndex)]=0;
        std::set<int> removedNodeA;
        std::set<int> removedNodeB;
        std::list<int> resultp;
        uniqueNodesforp(minIndex,resultp);
        std::list<int> resultap;
        uniqueAboveNodesforp(minIndex /*,visitedNode*/,resultap);
        if(resultap.size()!=uniqueAbovecount[minIndex])
        {printf("ERRR in uniqueAboveNodesforp %d %d\n",resultap.size(),uniqueAbovecount[minIndex]); getchar();}
        if(resultp.size()!=uniquecount[minIndex])
        {printf("ERRR in uniqueNodesforp\n"); getchar();}
        removedNodeB.insert(resultp.begin(), resultp.end());
        // removedNodeB.erase(minIndex);
        removedNodeA.insert(resultap.begin(), resultap.end());
        delete uniqueNodes;
        // printf("MAke MAP! \n" );

        // printf("getNodeStatus %d\n",getNodeStatus((minIndex)));
        // int deleted=uniquecount[minIndex];
        // printf("CALL RemoveDuplicate2 \n" );
        std::set<int> s(removedNodeA);
        // printf("MinIndex is %d\n",minIndex );
        // removedNodeB.erase(minIndex);
        // printf("******\n" );
                  #ifdef DBG_MTMDD
        s.insert(removedNodeB.begin(), removedNodeB.end());
        for (auto it = removedNodeA.begin(); it !=
             removedNodeA.end(); ++it)
                // if(incomingedgecount[(*it)]==1)
                printf("%d map %d iec %d lvl %d out of %d\n",(*it),map[(*it)],incomingedgecount[(*it)],getNodeLevel((*it)),getNumVariables() );
        printf("BEfore**RNA***\n" );
        for (auto it = removedNodeB.begin(); it !=
             removedNodeB.end(); ++it)
                // if(incomingedgecount[(*it)]==1)
                printf("%d map %d iec %d lvl %d out of %d\n",(*it),map[(*it)],incomingedgecount[(*it)],getNodeLevel((*it)),getNumVariables() );
        printf("BEfore**RNB***\n" );

        printf("BEfore RNA %d RNB %d Union %d, minlvl %d\n",removedNodeA.size(),removedNodeB.size(),s.size(),getNodeLevel(minIndex) );
                  #endif
        // printf("Cardinality before RD %ld\n", e.getCardinality());
        RemoveDuplicate2(lvl,map,e,removedNodeA,removedNodeB);
        num_deletedNode++;
        // printf("Cardinality after RD %ld\n", e.getCardinality());

                  #ifdef DBG_MTMDD
        printf("After RNA %d RNB %d\n",removedNodeA.size(),removedNodeB.size() );
        for (auto it = removedNodeB.begin(); it !=
             removedNodeB.end(); ++it)
                // if(incomingedgecount[(*it)]==1)
                printf("%d map %d iec %d lvl %d out of %d\n",(*it),map[(*it)],incomingedgecount[(*it)],getNodeLevel((*it)),getNumVariables() );
        printf("**RNB***\n" );
        // getchar();
                  #endif
        cC=e.getNodeCount();
        // printf("End CALL RemoveDuplicate2 \n" );

        // cC--;
        // if(belowcount!=0)delete[] belowcount; else {printf("ERRR belowcount\n"); getchar();}
        if(incomingedgecount!=0) delete [] incomingedgecount; else {printf("ERRR incomingedgecount\n"); getchar();}
        // if(abovecount!=0)delete[] abovecount;else {printf("ERRR abovecount\n"); getchar();}
        highestunique.clear();
        if(uniquecount!=0) delete[] uniquecount; else {printf("ERRR uniquecount\n"); getchar();}
        lowestunique.clear();
        if(uniqueAbovecount!=0) delete[] uniqueAbovecount; else{printf("ERR uniqueAbovecount\n" ); getchar();}
        // printf("DELETED \n" );
        delete[] ACBC;


        map.clear();
        // e.show(meddlyout,0);
        // printf("NodeCount is %d\n",cC );
}
/////////////*******
// printf("NddeCountAfter %d\n",e.getNodeCount() );
printf("deleted_node %d\n",num_deletedNode );
end = clock();
double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
printf("Time taken %f \n",time_taken );
// printf("underApproximate End\n" );
}



void MEDDLY::mt_mdd_bool::uniqueNodesforp(node_handle a,std::list<int> &result ){
result.push_back(a);// result.insert(a);
// printf("a in uniqueNodesforp is %d\n",a );

// std::set<int> rset=highestunique[a];
// for (auto it=rset.begin(); it != rset.end(); ++it){
for (std::list<int>::iterator it=highestuniquelist[a].begin(); it != highestuniquelist[a].end(); ++it){
   // printf("Res %d\n",*it );
   uniqueNodesforp((*it),result);
}
}
void MEDDLY::mt_mdd_bool::uniqueAboveNodesforp(node_handle a,std::list<int> &result ){

// std::set<int> rset=lowestunique[a];
// for (auto it=rset.begin(); it != rset.end(); ++it){
for (std::list<int>::iterator it=lowestuniquelist[a].begin(); it != lowestuniquelist[a].end(); ++it){

   result.push_back((*it));// result.insert((*it));
   uniqueAboveNodesforp((*it),result);

}
}
void MEDDLY::mt_mdd_bool::getNC(int k,node_handle a,bool* visitedNode,std::set<int> &result){
result.insert(a);
// printf("getNC inserted %d\n",a );
// printf("VN %d\n",visitedNode[a] );
if(visitedNode[a]==false){
     int kdn = k-1;

 int b=getDomain()->getVariableBound(k, false);
 // printf("b is %d\n",b );

     for(int ik=0;ik<=b; ik++){

         int dpt=getDownPtr(a,ik);
 // printf("child is %d kdn is %d\n",dpt,kdn );
         // if(dpt==1868)
         // {
         //     printf("1868 Reached from %d\n",a );
         //     getchar();
         // }
         // if(dpt==34348)
         // {
         //     printf("34348 Reached from %d\n",a );
         //     getchar();
         // }
         // if(dpt==34366)
         // {
         //     printf("34348 Reached from %d\n",a );
         //     getchar();
         // }
         // if(kdn>0&&dpt>1){
         //     // printf("Child\n",dpt );
         //
         //     getNC(kdn,dpt,visitedNode,result);
         // }
     }
     visitedNode[a]=true;

}
}
void MEDDLY::mt_mdd_bool::merge(doubleDensityClass* array, int const left, int const mid, int const right)
{
	auto const subArrayOne = mid - left + 1;
	auto const subArrayTwo = right - mid;

	// Create temp arrays
	auto *leftArray = new doubleDensityClass[subArrayOne],
		*rightArray = new doubleDensityClass[subArrayTwo];

	// Copy data to temp arrays leftArray[] and rightArray[]
	for (auto i = 0; i < subArrayOne; i++)
		{
            leftArray[i] = array[left + i];
            // if(array[left + i].density>0)
            //     printf("%lf %lf\n", array[left + i].density, leftArray[i].density);
        }
	for (auto j = 0; j < subArrayTwo; j++)
		rightArray[j] = array[mid + 1 + j];

	auto indexOfSubArrayOne = 0, // Initial index of first sub-array
		indexOfSubArrayTwo = 0; // Initial index of second sub-array
	int indexOfMergedArray = left; // Initial index of merged array

	// Merge the temp arrays back into array[left..right]
	while (indexOfSubArrayOne < subArrayOne && indexOfSubArrayTwo < subArrayTwo) {
		if (leftArray[indexOfSubArrayOne].density<=rightArray[indexOfSubArrayTwo].density) {
			array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
			indexOfSubArrayOne++;
		}
		else {
			array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
			indexOfSubArrayTwo++;
		}
		indexOfMergedArray++;
	}
	// Copy the remaining elements of
	// left[], if there are any
	while (indexOfSubArrayOne < subArrayOne) {
		array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
		indexOfSubArrayOne++;
		indexOfMergedArray++;
	}
	// Copy the remaining elements of
	// right[], if there are any
	while (indexOfSubArrayTwo < subArrayTwo) {
		array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
		indexOfSubArrayTwo++;
		indexOfMergedArray++;
	}
}

// begin is for left index and end is
// right index of the sub-array
// of arr to be sorted */
void MEDDLY::mt_mdd_bool::mergeSort(doubleDensityClass* array, int const begin, int const end)
{
	if (begin >= end)
		return; // Returns recursively

	auto mid = begin + (end - begin) / 2;
	mergeSort(array, begin, mid);
	mergeSort(array, mid + 1, end);
	merge(array, begin, mid, end);
}

void MEDDLY::mt_mdd_bool::HeuristicUnderApproximate(dd_edge &e, long Threashold, long maxThreshold,float desiredPercentage, int option, int deletedApproach, float rootStatePercentage)
{
        unsigned nodeCount=e.getNodeCount();
        if(option!=0|| nodeCount<maxThreshold) return;
        srand(time(0));
        clock_t start, end;
        start = clock();
        FILE_output meddlyout(stdout);
        int num_vars=getNumVariables();
#ifdef PRINTON
        printf("HunderApproximate\n" );
#endif
        while((deletedApproach==1)||(nodeCount>Threashold)) {
                maxid=e.getLastHandle();
                lastNode=maxid+1;
#ifdef PRINTON
                clock_t startc= clock();
#endif
                double nodeCountard;
                apply(BC,e,nodeCountard);
#ifdef PRINTON
                printf("Below count time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
                startc = clock();
#endif
                double cI;
                apply(IEC, e, cI);
#ifdef PRINTON
                printf("IEC time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
                startc = clock();
#endif
                double cA;
                apply(AC, e, cA);
#ifdef PRINTON
                printf("AC time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
                startc = clock();
#endif
#ifdef HAVE_LIBGMP
                mpz_object mpzd;
                mpzd.setValue(0);
                mpzd.setReminder(0);
                mpz_object mpz0;
                mpz0.setValue(0);
                mpz0.setReminder(0);
                mpz_object mpz1;
                mpz1.setValue(1);
                mpz1.setReminder(0);
                ACBC=new mpz_object[lastNode];
#else
                ACBC=new long double[lastNode];
#endif
                node_handle root=e.getNode();
                node_handle* list = markNodesInSubgraph(&root, 1, false);
                std::list<int>* uniquelist=markNodesInSubgraphByLvl(&root,1);
                for(long i=0; i<lastNode; i++){
#ifdef HAVE_LIBGMP
                ACBC[i].setValue(0);
                ACBC[i].setReminder(0);
#else
                ACBC[i]=0.0;
#endif
                }

#ifdef HAVE_LIBGMP
                mpz_object sumOfstateforselectedNode;
                sumOfstateforselectedNode.setValue(0);
#else
                long double sumOfstateforselectedNode=0.0;
#endif
                long long int rootNumberofState=belowcount[root];
                for (long i=0; list[i]; i++) {
#ifdef HAVE_LIBGMP
                        ACBC[list[i]].setValue(abovecount[list[i]]);
                        ACBC[list[i]].multiply(belowcount[list[i]]);
#else
                        ACBC[list[i]]=(long double) abovecount[list[i]]*belowcount[list[i]];
                        if (ACBC[list[i]]/abovecount[list[i]]!=belowcount[list[i]]) {
#ifdef PRINTON
                                printf("NOT CORRECT ACBC %Lf, %llu, %llu\n",ACBC[i],abovecount[list[i]],belowcount[list[i]]);
#endif
                                throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
                                exit(0);
                        }
#endif
                }
                delete[] abovecount;
                delete[] belowcount;
                double cH;
                apply(HU, e, cH);
#ifdef PRINTON
                printf("HU time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
                startc = clock();
#endif
                double cU;
                apply(UC, e, cU);
#ifdef PRINTON
                printf("UC time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
                startc = clock();
#endif
                double cL;
                apply(LU, e, cL);
#ifdef PRINTON
                printf("LU time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
                startc = clock();
#endif
                double caU;
                apply(UAC, e, caU);
#ifdef PRINTON
                printf("UAC time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
#endif
                int* levelcount=new int[num_vars+1];
                for(int i=0; i<=num_vars; i++) {
                        levelcount[i]=uniquelist[i].size();
                }
                int minIndex=0;
                std::list<int> selectedNodeforDeletion;
                std::list<int> removedNode;
#ifdef PRINTON
                startc = clock();
#endif

#ifdef HAVE_LIBGMP
                mpzDensityClass* doubleDensityArray=new mpzDensityClass[lastNode];
#else
                doubleDensityClass* doubleDensityArray=new doubleDensityClass[lastNode];
#endif
                bool* isInArr= new bool[lastNode];
                bool* neverInArr= new bool[lastNode];
                for(long i=0; i<lastNode; i++)
                {
#ifdef HAVE_LIBGMP
                        doubleDensityArray[i].density.setValue(-1);
                        doubleDensityArray[i].density.setReminder(0);
                        doubleDensityArray[i].index=-1;
#else
                        doubleDensityArray[i].density=-1;
                        doubleDensityArray[i].index=-1;
#endif
                        isInArr[i]=false;
                        neverInArr[i]=false;
                }
                for (long i=0; list[i]; i++) {
                        if(list[i]>lastNode) {
#ifdef PRINTON
                                printf("ERRROR LASTNODE %d %d\n", list[i],lastNode);
#endif
                        }
                        if(list[i]!=root) {
#ifdef HAVE_LIBGMP
                                doubleDensityArray[list[i]].density=ACBC[list[i]];
                                doubleDensityArray[list[i]].density.division((long)(uniquecount[list[i]]+uniqueAbovecount[list[i]]));
#else
                                doubleDensityArray[list[i]].density=(long double)ACBC[list[i]]/(double)(uniquecount[list[i]]+uniqueAbovecount[list[i]]); //((double)(abovecount[list[i]])/(double)(uniquecount[list[i]]+uniqueAbovecount[list[i]]))*(double)belowcount[list[i]];
#endif
                                doubleDensityArray[list[i]].index=list[i];
                                isInArr[list[i]]=false;
                                neverInArr[list[i]]=false;
                        }
                }
#ifdef PRINTON
                printf("Density time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
                startc = clock();
#endif
                std::random_device rd;
                std::mt19937 g(rd());
                std::shuffle(doubleDensityArray, doubleDensityArray+lastNode, g);
#ifdef HAVE_LIBGMP
std::sort(doubleDensityArray, doubleDensityArray+lastNode);

// std::sort(doubleDensityArray, doubleDensityArray+lastNode, [](const doubleDensityArray& x, const doubleDensityArray& y) {
//            mpz_object dmpz;
//            if(dmpz.compare(x,y)==-1)
//            {
//            return true;
//             }else{
//            return false;
//             }
//     });
#else
                std::sort(doubleDensityArray, doubleDensityArray+lastNode);
#endif
#ifdef PRINTON
                printf("Sorting time doubleDensityArray %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
                startc = clock();
#endif
                removedNode.clear();
                selectedNodeforDeletion.clear();
                int shouldberemoved=0;
#ifdef HAVE_LIBGMP
                // int x=0;
                // for(int i=0; i<lastNode;i++){
                //   if(x==100) break;
                //   mpz0.setValue(0);
                //   mpz0.setReminder(0);
                //   if(mpzd.compare(doubleDensityArray[i].density,mpz0)==1){
                //     doubleDensityArray[i].density.showwithreminder(meddlyout);
                //     printf("\n" );
                //     x++;
                //   }
                // }
                // getchar();
                sumOfstateforselectedNode.setValue(0);
#else
                sumOfstateforselectedNode=0.0;
#endif
                int dsaIndexxk=0;
                while(nodeCount-shouldberemoved>Threashold) {
                        minIndex=0;
                        int numberDeleted=0;
                        for ( int dsaIndexx=dsaIndexxk; dsaIndexx<lastNode; dsaIndexx++) {
                                int i= doubleDensityArray[dsaIndexx].index;
                                if(i>lastNode) {
#ifdef PRINTON
                                        printf("ERRROR LASTNODE %d %d\n", i,lastNode);
#endif
                                }
#ifdef HAVE_LIBGMP
                                mpzd.setValue(0);
                                mpzd.setReminder(0);
                                if (i>0&& mpzd.compare(doubleDensityArray[dsaIndexx].density,mpz0)==1)
#else
                                if (i>0&& doubleDensityArray[dsaIndexx].density>0.0)
#endif
                                {
                                        if(i!=root)
#ifdef HAVE_LIBGMP
                                                if((mpzd.compare(ACBC[i],mpz0)==1)&&(levelcount[getNodeLevel(i)]>1))
#else
                                                if((ACBC[i]>0)&&(levelcount[getNodeLevel(i)]>1))
#endif
                                                {
                                                        bool is_in =isInArr[i];
                                                        bool neverDelete_isin =neverInArr[i];
                                                        if(!is_in && !neverDelete_isin &&(i>1)) {
                                                                minIndex=i;
                                                                numberDeleted++;
                                                                dsaIndexxk=dsaIndexx+1;
                                                                break;
                                                        }else{
                                                                numberDeleted++;
                                                        }
                                                }
                                }
                        }
                        if(minIndex!=0) {
                                std::list<int> resultp;
                                uniqueNodesforp(minIndex,resultp);
                                std::list<int> resultap;
                                uniqueAboveNodesforp(minIndex,resultap);
                                if(resultap.size()!=uniqueAbovecount[minIndex])
                                {
#ifdef PRINTON
                                        printf("ERRR in uniqueAboveNodesforp %ld %d\n",resultap.size(),uniqueAbovecount[minIndex]);
#endif
                                        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
                                        exit(0);
                                        // getchar();
                                }
                                if(resultp.size()!=uniquecount[minIndex])
                                {
#ifdef PRINTON
                                        printf("ERRR in uniqueNodesforp\n"); //getchar();
#endif
                                        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
                                        exit(0);
                                }
                                if(levelcount[getNodeLevel(minIndex)]>1) { //should add
                                        // startf= clock();
#ifdef HAVE_LIBGMP
                                        sumOfstateforselectedNode.add(ACBC[minIndex]);
                                        mpz_object rootNumberofStatePercentage;
                                        // mpz_object mpz1;
                                        // mpz1.setValue(1);
                                        // mpz1.setReminder(0);
                                        // mpz_object mpzd;
                                        // mpzd.setValue(0);
                                        // mpzd.setReminder(0);
                                        // mpzd=mpz1;
                                        // mpzd.division(rootStatePercentage);
                                        long drootStatePercentage=(double)1.0/rootStatePercentage;
                                        // printf("%f 1/rootPerc %ld\n",rootStatePercentage, drootStatePercentage);
                                        // getchar();

                                        rootNumberofStatePercentage.setValue(rootNumberofState);
                                        rootNumberofStatePercentage.setReminder(0);
                                        rootNumberofStatePercentage.division(drootStatePercentage);
                                          // printf("r, %d \n",rootNumberofState );
                                          // rootNumberofStatePercentage.showwithreminder(meddlyout);
                                          // getchar();
                                        // =(long double)((long double)rootNumberofState/ ((long double)1.0/rootStatePercentage));

#else
                                        sumOfstateforselectedNode+=(long double)ACBC[minIndex];
                                        long double rootNumberofStatePercentage=(long double)((long double)rootNumberofState/ ((long double)1.0/rootStatePercentage));


                                        if (rootNumberofStatePercentage<0)
                                        {
#ifdef PRINTON
                                                printf("NOT CORRECT %Lf, %lld, %f\n",rootNumberofStatePercentage,rootNumberofState,rootStatePercentage );
#endif
                                                throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
                                                exit(0);
                                        }
#endif
#ifdef HAVE_LIBGMP
                                        // mpz_object mpzd;
                                        mpzd.setValue(0);
                                        mpzd.setReminder(0);
                                        if((deletedApproach==1)&&(selectedNodeforDeletion.size()>0)&&(mpzd.compare(sumOfstateforselectedNode,rootNumberofStatePercentage)==1)) {
#else
                                        if((deletedApproach==1)&&(selectedNodeforDeletion.size()>0)&&(sumOfstateforselectedNode>rootNumberofStatePercentage)) {
#endif
                                                break;
                                        }
                                        else{
                                                selectedNodeforDeletion.push_back(minIndex);
                                        }
                                        for(auto i= resultp.begin(); i!=resultp.end(); ++i) {
                                                if(!isInArr[(*i)]) {
                                                        isInArr[(*i)]=true;
                                                        levelcount[getNodeLevel(*i)]--;
                                                        shouldberemoved++;
                                                }
                                        }
                                        for(auto i= resultap.begin(); i!=resultap.end(); ++i) {
                                                if(!isInArr[(*i)]) {
                                                        isInArr[(*i)]=true;
                                                        levelcount[getNodeLevel(*i)]--;
                                                        shouldberemoved++;
                                                }
                                        }
                                }
                                else{neverInArr[minIndex]=true;}
                        }else{
                                break;
                        }
                }
#ifdef PRINTON
                printf("Selecting nodes for deletion time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
#endif
                if(selectedNodeforDeletion.size()==0) {
#ifdef PRINTON
                        printf("No node can get selected.\n" );
#endif
                        return;
                }
#ifdef PRINTON
                printf("selectedNodeforDeletion size %d\n",selectedNodeforDeletion.size() );
                printf("Removed size %d\n",shouldberemoved);
                printf("Expected %d\n", nodeCount-shouldberemoved);
#ifdef HAVE_LIBGMP
                printf("sumOfstate for selectedNodex ");
                sumOfstateforselectedNode.showwithreminder(meddlyout);
                printf("\n");
#else
                printf("sumOfstate for selectedNodex %Lf\n",sumOfstateforselectedNode );
#endif
#endif
                delete[]levelcount;
                delete[] doubleDensityArray;
                delete[] isInArr;
                delete[] neverInArr;
                std::map<int,int> map;
                bool* levelarray=new bool[num_vars+1];
                for(int i=0; i<=num_vars; i++)
                        levelarray[i]=false;
                int minlvl=INT_MAX;
                for(long i=0; i<=lastNode; i++)
                        map[i]=0;
                for (long i=0; list[i]; i++) {
                        map[list[i]]=list[i];
                }
                for(auto i = selectedNodeforDeletion.begin(); i != selectedNodeforDeletion.end(); ++i) {
                        map[(*i)]=0;
                        int nodelvl=getNodeLevel((*i));
                        levelarray[nodelvl]=true;
                        if(nodelvl<minlvl)
                                minlvl=nodelvl;
                }
#ifdef PRINTON
                printf("Cardinality before RD %llu\n", rootNumberofState);
                startc=clock();
#endif
                int merged=RemoveDuplicateSet(minlvl,levelarray,uniquelist, map,e );
#ifdef PRINTON
                printf("RemoveDuplicateSet time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
#endif
                nodeCount=e.getNodeCount();
                // printf("Node count is %d\n",nodeCount);
                if((deletedApproach==1)&&(nodeCount<=Threashold)) {
                        deletedApproach=0;
                }
                if(incomingedgecount!=0) delete [] incomingedgecount; else {
#ifdef PRINTON
                        printf("ERRR incomingedgecount\n");
#endif
                }
                highestunique.clear();
                if(highestuniquelist!=0) delete[] highestuniquelist; else {
#ifdef PRINTON
                        printf("ERRR highestuniquelist\n");
#endif
                }
                if(uniquecount!=0) delete[] uniquecount; else {
#ifdef PRINTON
                        printf("ERRR uniquecount\n");
#endif
                }
                lowestunique.clear();
                if(lowestuniquelist!=0) delete[] lowestuniquelist; else {
#ifdef PRINTON
                        printf("ERRR lowestuniquelist\n");
#endif
                }
                if(uniqueAbovecount!=0) delete[] uniqueAbovecount; else{
#ifdef PRINTON
                        printf("ERR uniqueAbovecount\n" );
#endif
                }
                delete[] ACBC;
                delete[] levelarray;
                delete[] list;
                map.clear();
        }
#ifdef PRINTON
        end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        printf("Time taken %f \n",time_taken );
        printf("HunderApproximate EndZ\n" );
#endif

}

////////////////////////////////REMOVE
      // void MEDDLY::mt_mdd_bool::HeuristicUnderApproximateDeleted(dd_edge &e, long Threashold, long maxThreshold,float desiredPercentage, int option, int deletedApproach, float rootStatePercentage)
      // {
      //     if(option==1) return;
      //     srand(time(0));
      //          clock_t start, end;
      //          start = clock();
      //         FILE_output meddlyout(stdout);
      //         int num_vars=getNumVariables();
      //         int cC=e.getNodeCount();//unique->getNumEntries();//this->getCurrentNumNodes();
      //         if(cC<maxThreshold) {
      //             // printf("LESS THAN Threashold\n" );
      //             return;
      //         }
      //          printf("HunderApproximate\n" );
      //
      //         int initialRootCardinality=e.getCardinality();
      //         double initialRootDensity=((double)initialRootCardinality/cC)*desiredPercentage;
      //         #if 0
      //         std::random_device rd;
      //     std::default_random_engine generator;
      //     generator.seed( rd() );
      //         std::normal_distribution<> d{(double)Threashold,((double)Threashold*0.05)};
      //         Threashold=std::round(d(generator));
      //         printf("Threashold %d\n",Threashold );
      //
      // #endif
      //         // while(cC>Threashold){
      //         while((deletedApproach==1)||((e.getNodeCount()>Threashold)&&((option==0)||(option==3)))||(option==2)){
      //             // printf("START calculation \n" );
      //             cC=e.getNodeCount();
      //          maxid=e.getLastHandle();
      //         lastNode=maxid+1;
      //         clock_t startc= clock();
      //             double cCard;
      //             apply(BC,e,cCard);
      //         printf("Below count time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      //         startc = clock();
      //         double cI;
      //         apply(IEC, e, cI);
      //         printf("IEC time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      //         startc = clock();
      //        double cA;
      //        apply(AC, e, cA);
      //        printf("AC time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      //        startc = clock();
      //        ACBC=new long double[lastNode];
      //        printf("ACBC created\n" );
      //        node_handle nl=e.getNode();
      //
      //        node_handle* list = markNodesInSubgraph(&nl, 1, false);
      //        std::list<int>* uniquelist=markNodesInSubgraphByLvl(&nl,1);
      //        for(long i=0;i<lastNode;i++)
      //        ACBC[i]=0.0;
      //        // for (int i=0; list[i]; i++) {
      //        //     if(abovecount[i]<1)
      //        //     {printf("Abovecount not correct %d %d %llu\n",i, incomingedgecount[i],abovecount[i] );
      //        //     getchar();
      //        // }
      //        // }
      //        long double sumOfstateforselectedNode=0.0;
      //        long long int rootNumberofState=belowcount[e.getNode()];
      //        for (int i=0; list[i]; i++) {
      //        // for(int i=0;i<lastNode;i++){
      //            // if(abovecount[i]>1)
      //            ACBC[list[i]]=(long double) abovecount[list[i]]*belowcount[list[i]];
      //            if (ACBC[list[i]]/abovecount[list[i]]!=belowcount[list[i]])
      //            { printf("NOT CORRECT%llu,%lu, %lu\n",ACBC[i],abovecount[list[i]],belowcount[list[i]] );
      //            exit(0);
      //            // getchar();
      //             }
      //            // else
      //            // ACBC[i]=0;
      //        }
      //        // printf("ACBC calculated\n" );
      //        delete[] abovecount;
      //         printf("abovecount deleted\n" );
      //        delete[] belowcount;
      //        printf("belowcount deleted\n" );
      //        // delete[] incomingedgecount;
      //        printf("incomingedgecount deleted\n" );
      //        // if(incomingedgecount!=0)delete [] incomingedgecount;else {printf("ERRR incomingedgecount\n"); getchar();}
      //
      //
      //        double cH;
      //        apply(HU, e, cH);
      //        printf("HU time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      //        startc = clock();
      //        double cU;
      //        apply(UC, e, cU);
      //        printf("UC time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      //        startc = clock();
      //        double cL;
      //        apply(LU, e, cL);
      //        printf("LU time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      //        startc = clock();
      //        double caU;
      //        apply(UAC, e, caU);
      //        printf("UAC time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      //
      //
      //
      //
      //        // printf("rootNumberofState %ld \n", rootNumberofState);
      //
      //
      //        int* levelcount=new int[num_vars+1];
      //        // int* copylevelcount=new int[num_vars+1];
      //        for(int i=0;i<=num_vars; i++){
      //            levelcount[i]=uniquelist[i].size();
      //            // copylevelcount[i]=0;
      //        }
      //
      //        // for(int i=0;i<=(int)maxid;i++){
      //        //     if((abovecount/*ACBC*/[i]>0)||(i==e.getNode()))
      //        //      {
      //        //          levelcount[getNodeLevel(i)]++;
      //        //          // copylevelcount[getNodeLevel(i)]++;
      //        //      }
      //        // }
      //
      //        int minIndex=0;
      //        std::list<int> m;
      //        std::list<int> removedNode;
      //        // std::set<int> removedNodeA;
      //        // std::set<int> removedNodeB;
      //        // std::set<int> neverdelete;
      //
      //
      //        if(option==0||(option==2)||(option==3)){
      //        // double mindensity=DBL_MAX;
      //        // mpz_object* arrdensity=new mpz_object[maxid+1];
      //        // std::vector<mpz_object> vectordensity;
      //
      //        // for(int i=0;i<=maxid;i++){
      //        //     arrdensity[i].setValue(LONG_MAX);
      //        //     arrdensity[i].setIndex(i);
      //        //     arrdensity[i].setReminder(0);
      //        //  }
      //
      //           // double density;
      //           // mpz_object densitympz;
      //           // densitympz.setValue(LONG_MAX);
      // ///////////////////////////////////START RETURN////
      //        // for(int i=0;i<=(int)maxid;i++){
      //        //     if((uniquecount[i]>0)&&(levelcount[getNodeLevel(i)]>1))//&&(abovecount[i]>0)&&(belowcount[i]>0)&&(uniquecount[i]>0)&&(incomingedgecount[i]>0||i==e.getNode()))
      //        //     {
      //        //         // mpz_object mulmpz;
      //        //         // mulmpz.setValue(abovecount[i]);
      //        //         // mulmpz.multiply(belowcount[i]);
      //        //         // mulmpz.division(uniquecount[i]);
      //        //
      //        //
      //        //         arrdensity[i].setValue(abovecount[i]);
      //        //         arrdensity[i].multiply(belowcount[i]);
      //        //         arrdensity[i].division(uniquecount[i]+uniqueAbovecount[i]);
      //        //         arrdensity[i].setIndex(i);
      //        //
      //        //     }
      //        // }
      // ///////////////////////////////END RETURN ///////////////////////////////
      //        // printf("vectordensity size %d\n",vectordensity.size() );
      //        // vectordensity.clear();
      //        startc = clock();
      //        long i;
      //        node_handle root=e.getNode();
      //        // struct densityStruct {
      //        //     mpz_object density;
      //        //     int index;
      //        //     bool removed;
      //        //     bool neverShouldRemove;
      //        //  };
      //         // class doubleDensityStruct {
      //         //     double density;
      //         //     int index;
      //         //     bool isIn;
      //         //     bool neverIn;
      //         //  };
      //          doubleDensityClass* doubleDensityArray=new doubleDensityClass[lastNode];
      //          bool* isInArr= new bool[lastNode];
      //          bool* neverInArr= new bool[lastNode];
      // 		 for(long i=0;i<lastNode;i++)
      //          {
      // 			 doubleDensityArray[i].density=-1;
      // 			 doubleDensityArray[i].index=-1;
      // 			 // doubleDensityArray[i].isIn=false;
      // 			 // doubleDensityArray[i].neverIn=false;
      //              isInArr[i]=false;
      //              neverInArr[i]=false;
      // 		 }
      //         // densityStruct* densityStructArray=new densityStruct[lastNode];
      //         // mpz_object* densityArray=new mpz_object[lastNode];
      //         // double* doubleDensityArray= new double[lastNode];
      //         // bool* isInDensityArray=new bool[lastNode];
      //         // bool* neverInDensityArray=new bool[lastNode];
      //
      //         for (i=0; list[i]; i++) {
      //             if(list[i]>lastNode){
      //                 printf("ERRROR LASTNODE %d %d\n", list[i],lastNode);
      //             }
      //             if(list[i]!=root){
      //                 doubleDensityArray[list[i]].density=(long double)ACBC[list[i]]/(double)(uniquecount[list[i]]+uniqueAbovecount[list[i]]);//((double)(abovecount[list[i]])/(double)(uniquecount[list[i]]+uniqueAbovecount[list[i]]))*(double)belowcount[list[i]];
      //                 doubleDensityArray[list[i]].index=list[i];
      //                 // doubleDensityArray[list[i]].isIn=false;
      //                 // doubleDensityArray[list[i]].neverIn=false;
      //                 isInArr[list[i]]=false;
      //                 neverInArr[list[i]]=false;
      //                 // mpz_object vdensity;
      //                 // densityArray[list[i]].setValue(abovecount[list[i]]);
      //                 // densityArray[list[i]].multiply(belowcount[list[i]]);
      //                 // densityArray[list[i]].division(uniquecount[list[i]]+uniqueAbovecount[list[i]]);
      //                 // densityArray[list[i]].setIndex(list[i]);
      //                 // densityArray[list[i]].isIn=false;
      //                 // densityArray[list[i]].neverIn=false;
      //                 // densityArray[list[i]].showwithreminder(meddlyout);
      //             }
      //         }
      //         // getchar();
      //
      //         // int dsaIndex=0;
      //         // for (i=0; list[i]; i++) {
      //         //     if(list[i]!=root){
      //         //         mpz_object vdensity;
      //         //         vdensity.setValue(abovecount[list[i]]);
      //         //         vdensity.multiply(belowcount[list[i]]);
      //         //         vdensity.division(uniquecount[list[i]]+uniqueAbovecount[list[i]]);
      //         //         vdensity.setIndex(list[i]);
      //         //         densityStructArray[list[i]].density=vdensity;
      //         //         densityStructArray[list[i]].index=list[i];
      //         //         densityStructArray[list[i]].removed=false;
      //         //         densityStructArray[list[i]].neverShouldRemove=false;
      //         //     }
      //         // }
      //
      // /////////////////////////////////////Should Delete
      //        // for (i=0; list[i]; i++) {
      //        //     if(list[i]!=root){
      //        //         mpz_object vdensity;
      //        //         vdensity.setValue(abovecount[list[i]]);
      //        //         vdensity.multiply(belowcount[list[i]]);
      //        //         vdensity.division(uniquecount[list[i]]+uniqueAbovecount[list[i]]);
      //        //         vdensity.setIndex(list[i]);
      //        //         vectordensity.push_back(vdensity);
      //        //     }
      //        // }
      //        /////////////////////////////////////Should Delete End
      //
      //        printf("Density time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      //        // printf("vectordensity size %d\n",vectordensity.size() );
      //        // getchar();
      //        // printf("vectordensity Stored\n" );
      //        // for (auto i: vectordensity){
      //        //     i.showwithreminder(meddlyout);
      //        //              printf("\n");
      //        // }
      //        // auto rng = std::default_random_engine {};
      //        // std::shuffle(vectordensity.begin(), vectordensity.end(), rng);
      //        // printf("lastNode %ld,size %ld \n",lastNode, sizeof(densityStructArray) / sizeof( densityStructArray[0]) );
      //    //     for ( int i=0;i<lastNode;i++){
      //    //
      //    //     densityArray[i].showwithreminder(meddlyout);
      //    //              printf("\n");
      //    // }
      //    // getchar();
      //    //     getchar();
      //    /////////////////////////////////GMP/////////////////
      //     //    startc = clock();
      //     //    // mergeSort(densityArray, 0, lastNode-1);
      //     //    std::sort(densityArray, densityArray+lastNode, [](const mpz_object& x, const mpz_object& y) {
      //     //        mpz_object dmpz;
      //     //        if(dmpz.compare(x,y)<0)
      //     //        {
      //     //        return true;
      //     //         }else{
      //     //        return false;
      //     //         }
      //     // });
      //     // printf("Sorting time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      // //////////////////////////////////////////////////////////
      // // int y=0;
      // // for ( int i=0;i<lastNode;i++){
      // //
      // //
      // // if((y<100)&(incomingedgecount[doubleDensityArray[i].index]>0)){
      // // printf("%lf, %d\n",doubleDensityArray[i].density,doubleDensityArray[i].index);
      // //          y++;
      // //      }
      // //  }
      // //  getchar();
      //     startc = clock();
      //     // mergeSort(doubleDensityArray, 0, lastNode-1);
      //     // getchar();
      // 	std::random_device rd;
      //    std::mt19937 g(rd());
      //
      //    std::shuffle(doubleDensityArray, doubleDensityArray+lastNode, g);
      //     std::sort(doubleDensityArray, doubleDensityArray+lastNode);
      //
      //
      //
      //     // , [](const doubleDensityClass& x, const doubleDensityClass& y) {
      //     //
      //     //        // if(x.density<y.density)
      //     //        // {
      //     //        //      if(x.density>0 && y.density>0){
      //     //        //     printf("%lf < %lf \n",x.density,y.density );
      //     //        //     getchar();
      //     //        // }
      //     //        return true;
      //     //         }else{
      //     //        return false;
      //     //         }
      //     // });
      //  printf("Sorting time doubleDensityArray %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      //  // int x=0;
      //  //
      //  // // for (auto i: doubleDensityArray){
      //  // doubleDensityClass d;
      //  //     for ( int i=0;i<lastNode;i++){
      //  //
      //  //     // i.clearBuffer();
      //  //     // if(densityArray[i].index>0 )
      //  //     if((incomingedgecount[doubleDensityArray[i].index]>0)){
      //  //         d=doubleDensityArray[i];
      //  //         if(doubleDensityArray[i].density==0.0){
      //  //             printf("%ld, %ld, %d\n",abovecount[doubleDensityArray[i].index],belowcount[doubleDensityArray[i].index],uniquecount[doubleDensityArray[i].index]+uniqueAbovecount[doubleDensityArray[i].index] );
      //  //             getchar();
      //  //         }
      //  //     printf("%lf, %d\n",doubleDensityArray[i].density,doubleDensityArray[i].index);
      //  //              // printf("\n");
      //  //              x++;
      //  //          }
      //  //      }
      //               // getchar();
      //  // }
      //  // getchar();
      //
      //     //    std::sort(densityStructArray, densityStructArray+lastNode, [](const densityStruct& x, const densityStruct& y) {
      //     //        mpz_object dmpz;
      //     //        if(dmpz.compare(x.density,y.density)<0)
      //     //        {
      //     //        return true;
      //     //         }else{
      //     //        return false;
      //     //         }
      //     // });
      //     /////////////////////////////////////Should Delete
      //
      //     //    std::sort(vectordensity.begin(), vectordensity.end(), [](const mpz_object& x, const mpz_object& y) {
      //     //        mpz_object dmpz;
      //     //        if(dmpz.compare(x,y)<0)
      //     //        {
      //     //        return true;
      //     //         }else{
      //     //        return false;
      //     //         }
      //     // });
      //     /////////////////////////////////////Should Delete End
      //
      //     // for (auto i: densityStructArray){
      //     // int x=0;
      //     //     for ( int i=0;i<lastNode;i++){
      //     //
      //     //     // i.clearBuffer();
      //     //     if(densityArray[i].index>0 )
      //     //     if(x<100){
      //     //     densityArray[i].showwithreminder(meddlyout);
      //     //              printf("\n");
      //     //              x++;
      //     //          }
      //     //              // getchar();
      //     // }
      //     // getchar();
      //
      //     // vectordensity.clear();
      //     // getchar();
      //     // for (auto i: vectordensity){
      //     //     i.showwithreminder(meddlyout);
      //     //              printf("\n");
      //     //              // getchar();
      //     // }
      //     // getchar();
      //     // printf("Sorted\n" );
      //     // for(int i=0;i<=(int)maxid;i++){
      //     //         arrdensity[i].showwithreminder(meddlyout);
      //     //         printf("\n");
      //     // }
      //     // getchar();
      //        // int ck=cC;
      //        // int root=e.getNode();
      //
      //        ////Not needed
      //        // mpz_object option2comparingmpz;
      //        // double intpart;
      //        // double fractpart = std::modf (initialRootDensity , &intpart);
      //        // option2comparingmpz.setValue(intpart);
      //        // option2comparingmpz.setReminder(fractpart);
      //        ////Not needed
      //        bool option3densitycheck=false;
      //        removedNode.clear();
      //         ////Not needed
      //        // removedNodeA.clear();
      //        // removedNodeB.clear();
      //         ////Not needed
      //        m.clear();
      //         // printf("**IN WHILE LOOP option %d %d\n", option, cC-removedNode.size()>Threashold );
      //         // for (auto i: vectordensity){
      //         //     printf("%d ,",i.index );
      //         // }
      //         // printf("\n" );
      //         // getchar();
      // ////////////////////////////////////////////////////////////////////////////////////////////////////
      // startc = clock();
      // int shouldberemoved=0;
      // sumOfstateforselectedNode=0.0;
      // for(int j=0;j<lastNode;j++){
      //     // if(densityStructArray[j].removed)
      //     // if(doubleDensityArray[j].isIn)
      // 	// doubleDensityArray[j].isIn=false;
      // 	// doubleDensityArray[j].neverIn=false;
      //     isInArr[j]=false;
      //     neverInArr[j]=false;
      // 	// if(densityArray[j].isIn)
      //     // shouldberemoved++;
      // }
      // int dsaIndexxk=0;
      // cC=e.getNodeCount();
      // long double selectednodeDensity=0.0;
      // while(((cC-shouldberemoved/*removedNode.size()*/>Threashold)&&((option==0)||(option==3)))||(option==2)){
      //     minIndex=0;
      //     clock_t startf= clock();
      //     int numberDeleted=0;
      //     selectednodeDensity=0.0;
      //     for ( int dsaIndexx=dsaIndexxk;dsaIndexx<lastNode;dsaIndexx++){
      //         int i= doubleDensityArray[dsaIndexx].index;//densityArray[dsaIndexx].index;
      //         if(i>lastNode){
      //             printf("ERRROR LASTNODE %d %d\n", i,lastNode);
      //         }
      //         if (i>0&& doubleDensityArray[dsaIndexx].density>0.0){
      //         if(i!=root)
      //         if(((ACBC[i]>0)&&(levelcount[getNodeLevel(i)]>1)&&((option==0)||(option==3 && option3densitycheck))||
      //         ((ACBC[i]>0)&&(levelcount[getNodeLevel(i)]>1)&&((option==2)||(option==3 && !option3densitycheck)))))
      //         {
      //             bool is_in =isInArr[i];//doubleDensityArray[dsaIndexx].isIn;//densityArray[dsaIndexx].isIn; //densityStructArray[dsaIndex].removed;//removedNode.find(i) != removedNode.end();
      //             bool neverDelete_isin =neverInArr[i];//doubleDensityArray[dsaIndexx].neverIn;//densityArray[dsaIndexx].neverIn; //densityStructArray[dsaIndex].neverShouldRemove; //neverdelete.find(i) != neverdelete.end();
      //             if(!is_in && !neverDelete_isin &&(i>1)){
      //             minIndex=i;
      //             selectednodeDensity=doubleDensityArray[dsaIndexx].density;
      //             numberDeleted++;
      //             dsaIndexxk=dsaIndexx+1;
      //             break;
      //         }else{
      //             numberDeleted++;
      //         }
      //          }
      //     }
      //     }
      //
      //     if(minIndex!=0){
      //
      //       std::list<int> resultp;
      //     uniqueNodesforp(minIndex/*,visitedNode*/,resultp);
      //
      //     std::list<int> resultap;
      //   uniqueAboveNodesforp(minIndex/*,visitedNode*/,resultap);
      //
      //   if(resultap.size()!=uniqueAbovecount[minIndex])
      //   {printf("ERRR in uniqueAboveNodesforp %d %d\n",resultap.size(),uniqueAbovecount[minIndex]);getchar();}
      //     if(resultp.size()!=uniquecount[minIndex])
      //     {printf("ERRR in uniqueNodesforp\n");getchar();}
      //     bool shouldadd=true;
      //
      //     if(shouldadd)
      //     if(levelcount[getNodeLevel(minIndex)]<=1){
      //         shouldadd=false;
      //     }
      //
      //     if(shouldadd){
      //         printf("Density  %Lf \n", selectednodeDensity);
      //         startf= clock();
      //         sumOfstateforselectedNode+=ACBC[minIndex];//(abovecount[minIndex]*belowcount[minIndex]);
      //         long double rootNumberofStatePercentage=(long double)(rootNumberofState/ (1.0/rootStatePercentage));
      //         if (rootNumberofStatePercentage<0)
      //         { printf("NOT CORRECT %Lf, %Lf, %f\n",rootNumberofStatePercentage,rootNumberofState,rootStatePercentage );
      //         exit(0);
      //          }
      //         if((deletedApproach==1)&&(m.size()>0)&&(sumOfstateforselectedNode>(long double)rootNumberofStatePercentage/*(rootNumberofState*rootStatePercentage)*/)){
      //             break;
      //         }
      //         else{
      //         m.push_back(minIndex);
      //         }
      //
      //         for(auto i= resultp.begin();i!=resultp.end();++i){
      //             // if(!doubleDensityArray[(*i)].isIn){
      //             if(!isInArr[(*i)]){
      //             isInArr[(*i)]=true;
      //             // doubleDensityArray[(*i)].isIn=true;//densityArray[(*i)].isIn=true;
      //             levelcount[getNodeLevel(*i)]--;
      //             shouldberemoved++;
      //             }
      //         }
      //         for(auto i= resultap.begin();i!=resultap.end();++i){
      //             // if(!doubleDensityArray[(*i)].isIn){
      //             if(!isInArr[(*i)]){
      //             isInArr[(*i)]=true;
      //             // doubleDensityArray[(*i)].isIn=true;//densityArray[(*i)].isIn=true;
      //             levelcount[getNodeLevel(*i)]--;
      //             shouldberemoved++;
      //             }
      //         }
      //
      //     }
      //     else{
      //         neverInArr[minIndex]=true;
      //         // doubleDensityArray[minIndex].neverIn=true;// densityArray[minIndex].neverIn=true;
      //     }
      //     }else{
      // 		break;
      //
      //     }
      //
      //    //  shouldberemoved=0;
      //    // for(int j=0;j<lastNode;j++){
      //    //     if(doubleDensityArray[j].isIn)// if(densityArray[j].isIn)
      //    //     shouldberemoved++;
      //    // }
      // }
      //
      // printf("Selecting nodes for deletion time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      // if(m.size()==0){
      //     printf("No node can get selected.\n" );
      //     return;
      // }
      //
      // ///////////////////////////////////////////////////BEGIN RETURN IF NEEDED////////////////////////////////////////////////
      //         // while(((cC-removedNode.size()>Threashold)&&((option==0)||(option==3)))||(option==2)){
      //         //     // printf("**IN WHILE LOOP\n" );
      //         //     mpz_object densitympz;
      //         //     densitympz.setValue(LONG_MAX);
      //         //     std::set<int> resultSet;
      //         //     mpz_object resultsetValue;
      //         //     resultsetValue.setValue(LONG_MAX);
      //         //     minIndex=0;
      //         //
      //         //     for(int i=0;i<=(int)maxid;i++){
      //         //
      //         //         // arrdensity[i].showwithreminder(meddlyout);
      //         //         // printf("\n");
      //         //         // option2comparingmpz.showwithreminder(meddlyout);
      //         //         // printf("\n");
      //         //         // printf("______\n" );
      //         //         if(i!=root)
      //         //         if(((incomingedgecount[i]>0)&&(arrdensity[i].compare(arrdensity[i], densitympz)<=0)&&(copylevelcount[getNodeLevel(i)]>1)&&((option==0)||(option==3 && option3densitycheck))||
      //         //         ((incomingedgecount[i]>0)&&(arrdensity[i].compare(arrdensity[i], option2comparingmpz)<=0)&&(copylevelcount[getNodeLevel(i)]>1)&&((option==2)||(option==3 && !option3densitycheck)))))
      //         //         //&&(!(removedNode.find(i))))
      //         //         {
      //         //             // printf("option3densitycheck %d size %d\n",option3densitycheck,removedNode.size());
      //         //              // printf("SELECTED\n" );
      //         //             bool is_in = removedNode.find(i) != removedNode.end();
      //         //             bool neverDelete_isin = neverdelete.find(i) != neverdelete.end();
      //         //             if(!is_in && !neverDelete_isin){
      //         //                 if(resultSet.size()>0){
      //         //                     if((option==0)||(option==3 && option3densitycheck)){
      //         //                         if(resultsetValue.compare(resultsetValue,arrdensity[i])==0){
      //         //                             resultSet.insert(i);
      //         //                             arrdensity[i].copyInto(resultsetValue);
      //         //                             resultsetValue.setReminder(arrdensity[i].rdvalue);
      //         //                         }else{
      //         //                             resultSet.clear();
      //         //                             resultSet.insert(i);
      //         //                              arrdensity[i].copyInto(densitympz);
      //         //                              densitympz.setReminder(arrdensity[i].rdvalue);
      //         //                             arrdensity[i].copyInto(resultsetValue);
      //         //                             resultsetValue.setReminder(arrdensity[i].rdvalue);
      //         //                         }
      //         //                     }else if((option==2)||(option==3 && !option3densitycheck)){
      //         //                         if(resultsetValue.compare(resultsetValue,arrdensity[i])==0){
      //         //                             resultSet.insert(i);
      //         //                             arrdensity[i].copyInto(resultsetValue);
      //         //                             resultsetValue.setReminder(arrdensity[i].rdvalue);
      //         //                         }else{
      //         //                             resultSet.clear();
      //         //                             arrdensity[i].copyInto(densitympz);
      //         //                             densitympz.setReminder(arrdensity[i].rdvalue);
      //         //                             resultSet.insert(i);
      //         //                             arrdensity[i].copyInto(resultsetValue);
      //         //                             resultsetValue.setReminder(arrdensity[i].rdvalue);
      //         //                         }
      //         //                     }
      //         //                 }else{
      //         //                    resultSet.insert(i);
      //         //                    arrdensity[i].copyInto(resultsetValue);
      //         //                    arrdensity[i].copyInto(densitympz);
      //         //                    densitympz.setReminder(arrdensity[i].rdvalue);
      //         //                    resultsetValue.setReminder(arrdensity[i].rdvalue);
      //         //                 }
      //         //
      //         //             // minIndex=i;
      //         //             // arrdensity[i].copyInto(densitympz);
      //         //             // densitympz.setReminder(arrdensity[i].rdvalue);
      //         //             }
      //         //         }
      //         //     }
      //         //     if(resultSet.size()>0){
      //         //        std::set<int>::iterator iter = resultSet.begin();
      //         //        int dgen=rand() % (resultSet.size());
      //         //       std::advance(iter, dgen);
      //         //        minIndex=(*(iter));
      //         //        // arrdensity[minIndex].showwithreminder(meddlyout);
      //         //        // getchar();
      //         //       resultSet.clear();
      //         //     }
      //         //     if(minIndex!=0){
      //         //
      //         //     // bool* visitedNode= new bool[maxid+1];
      //         //     // for(int i=0;i<=maxid;i++){
      //         //     //     visitedNode[i]=false;
      //         //     // }
      //         //       std::set<int> resultp;
      //         //     uniqueNodesforp(minIndex/*,visitedNode*/,resultp);
      //         //     std::set<int> resultap;
      //         //   uniqueAboveNodesforp(minIndex/*,visitedNode*/,resultap);
      //         //   // std::set<int> bresultp;
      //         //   // for(auto p:resultap){
      //         //   //     std::set<int> belowresultap;
      //         //   //     uniqueNodesforp(p/*,visitedNode*/,belowresultap);
      //         //   //     bresultp.insert(belowresultap.begin(),belowresultap.end());
      //         //   // }
      //         //   if(resultap.size()!=uniqueAbovecount[minIndex])
      //         //   {printf("ERRR in uniqueAboveNodesforp %d %d\n",resultap.size(),uniqueAbovecount[minIndex]);getchar();}
      //         //     if(resultp.size()!=uniquecount[minIndex])
      //         //     {printf("ERRR in uniqueNodesforp\n");getchar();}
      //         //     bool shouldadd=true;
      //         //     // for(auto i= bresultp.begin();i!=bresultp.end();++i){
      //         //     //     if(copylevelcount[getNodeLevel(*i)]<=1)
      //         //     //     {
      //         //     //         shouldadd=false;
      //         //     //         printf("should not add bresultp\n" );
      //         //     //         //getchar();
      //         //     //         neverdelete.insert(minIndex);
      //         //     //         break;
      //         //     //     }
      //         //     // }
      //         //
      //         //     if(shouldadd)
      //         //     for(auto i= resultp.begin();i!=resultp.end();++i){
      //         //         if(copylevelcount[getNodeLevel(*i)]<=1)
      //         //         {
      //         //             shouldadd=false;
      //         //             printf("should not add resultp\n" );
      //         //             //getchar();
      //         //             neverdelete.insert(minIndex);
      //         //             break;
      //         //         }
      //         //     }
      //         //     if(shouldadd)
      //         //     for(auto i= resultap.begin();i!=resultap.end();++i){
      //         //         if(copylevelcount[getNodeLevel(*i)]<=1)
      //         //         {
      //         //             shouldadd=false;
      //         //             printf("should not add resultap\n" );
      //         //             //getchar();
      //         //             neverdelete.insert(minIndex);
      //         //             break;
      //         //         }
      //         //     }
      //         //
      //         //     if(shouldadd){
      //         //         // printf("ADDED?\n" );
      //         //         sumOfstateforselectedNode+=(abovecount[minIndex]*belowcount[minIndex]);
      //         //         if((deletedApproach==1)&&(m.size()>0)&&(sumOfstateforselectedNode>(long int)(rootStatePercentage*rootNumberofState))){
      //         //             break;
      //         //         }
      //         //         else{
      //         //         m.insert(minIndex);
      //         //         }
      //         //
      //         //     removedNode.insert(resultp.begin(), resultp.end());
      //         //     removedNode.insert(resultap.begin(), resultap.end());
      //         //     removedNodeB.insert(resultp.begin(), resultp.end());
      //         //     // removedNodeB.insert(minIndex);
      //         //     removedNodeA.insert(resultap.begin(), resultap.end());
      //         //     for(int i=0;i<=num_vars; i++)
      //         //     copylevelcount[i]=levelcount[i];
      //         //     for(auto i= removedNode.begin();i!=removedNode.end();++i){
      //         //         copylevelcount[getNodeLevel(*i)]--;
      //         //     }
      //         //     // for(auto i = removedNode.begin(); i != removedNode.end(); ++i)
      //         //     // printf("%d, ",(*i) );
      //         //     // cC=ck-removedNode.size();
      //         //     // printf("\n%d %d %d\n",ck,removedNode.size(), cC );
      //         //     }
      //         //     else{
      //         //         neverdelete.insert(minIndex);
      //         //     }
      //         //     // getchar();
      //         //     }else{
      //         //         // printf("CAME ELSE\n" );
      //         //         // getchar();
      //         //         if(option==2 )
      //         //          {
      //         //              printf("Deleted node by density %d\n",removedNode.size() );
      //         //              // end = clock();
      //         //              // if(removedNode.size()==0){
      //         //              //     for(int i=0;i<(int)maxid;i++){
      //         //              //         if((incomingedgecount[i]>0||i==e.getNode())&&(uniqueAbovecount[i]>0))
      //         //              //         {
      //         //              //             printf("ERROR i %d AC %ld BC %ld UC %d UAC %d \n", i,abovecount[i], belowcount[i],uniquecount[i],uniqueAbovecount[i]);
      //         //              //             // char c=getchar();
      //         //              //          }
      //         //              //     }
      //         //              //     printf("NC is %ld UAC[0] is %ld\n",e.getNodeCount(), uniqueAbovecount[0]);
      //         //              //     getchar();
      //         //              // }
      //         //              // double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      //         //              // printf("Time taken %f \n",time_taken );
      //         //              break;
      //         //          }
      //         //         if(option==3 && !option3densitycheck){
      //         //             option3densitycheck=true;
      //         //             printf("Deleted node by density %d\n",removedNode.size() );
      //         //             if(removedNode.size()==0)
      //         //             { getchar();
      //         //
      //         //             }
      //         //             continue;
      //         //         }else{
      //         //         // printf("Cannot Delete more\n" );
      //         //         // getchar();
      //         //         // for(int i=0;i<num_vars;i++)
      //         //         // printf("%d %d\n",i, levelcount[i] );
      //         //         // getchar();
      //         //         break;
      //         //         // printf("Error\n" );
      //         //         }
      //         //     }
      //         //     // for(int i=0;i<=num_vars; i++)
      //         //     // copylevelcount[i]=levelcount[i];
      //         // }
      //         ////////////////////////////////////////////ENDRETURN IF NEEDED/////
      //         //  shouldberemoved=0;
      //         // for(int j=0;j<lastNode;j++){
      //         //     // if(densityStructArray[j].removed)
      //         //     if(doubleDensityArray[j].isIn)// if(densityArray[j].isIn)
      //         //     shouldberemoved++;
      //         // }
      //         printf("m size %d\n",m.size() );
      //         printf("Removed size %d\n",shouldberemoved);
      //         printf("Expected %d\n", cC-shouldberemoved);
      //         printf("sumOfstate for selectedNode %lu\n",sumOfstateforselectedNode );
      //         // printf("Removed size %d\n",removedNode.size() );
      //         // printf("Expected %d\n", cC-removedNode.size());
      //         // printf("rn %d\n",removedNode.size() );
      //
      //         // printf("cC%d rn %d\n",cC, removedNode.size() );
      //         // for(auto i = removedNode.begin(); i != removedNode.end(); ++i)
      //         // printf("%d, ",(*i) );
      //         // getchar();
      //
      //        // densitympz.showwithreminder(meddlyout);
      //        // if(minIndex==0)
      //        // {
      //        //  printf("ZERO %d\n",num_vars );
      //        //  char c=getchar();
      //        //  int k=0;
      //        //       for(int i=0;i<(int)maxid;i++){
      //        //
      //        //           printf("%d %d AC %ld BC %ld UC %d \n",k, i,abovecount[i], belowcount[i],uniquecount[i]);
      //        //           k++;
      //        //       }
      //        //       c=getchar();
      //        //       for(int i=0;i<=num_vars; i++){
      //        //           printf("level %d is %d\n", i,levelcount[i]);
      //        //       }
      //        //       c=getchar();
      //        //
      //        //       return;
      //        // }
      //        delete[]levelcount;
      //        // delete[]copylevelcount;
      //        // delete[] arrdensity;
      //        // delete[] densityStructArray;
      //        // delete[] densityArray;
      //        delete[] doubleDensityArray;
      //        delete[] isInArr;
      //        delete[] neverInArr;
      //         // vectordensity.clear();
      //         }
      //         if(option==1){
      //         //     node_handle n=e.getNode();
      //         //  node_handle* list = markNodesInSubgraph(&n, 1, false);
      //         //  long i;
      //         //  for (i=0; list[i]; i++) {
      //         //  if((uniquecount[list[i]]<=(cC-Threashold))&&(levelcount[getNodeLevel(list[i])]>1))
      //         //  {
      //         //      minIndex=list[i];
      //         //      break;
      //         //  }
      //         // }
      //         // if(minIndex==0)
      //         // printf("ERRROR \n" );
      //         //
      //         //  free(list);
      //     }
      //
      //
      //        std::map<int,int> map;
      //        // std::set<int> levels;
      //        int num_vars=getNumVariables();
      //        bool* levelarray=new bool[num_vars+1];
      //        for(int i=0;i<=num_vars;i++)
      //        levelarray[i]=false;
      //        int minlvl=INT_MAX;
      //        // node_handle n=e.getNode();
      //        // node_handle* list = markNodesInSubgraph(&n, 1, false);
      //        for(int i=0;i<=lastNode;i++)
      //        map[i]=0;
      //        for (long i=0; list[i]; i++) {
      //             map[list[i]]=list[i];
      //         }
      //        // delete[] list;
      //        int minnode=0;
      //        for(auto i = m.begin(); i != m.end(); ++i){
      //           map[(*i)]=0;
      //           int nodelvl=getNodeLevel((*i));
      //           levelarray[nodelvl]=true;
      //           // levels.insert(nodelvl);
      //           if(nodelvl<minlvl)
      //             minlvl=nodelvl;
      //         }
      //         for(auto i = m.begin(); i != m.end(); ++i){
      //             int nodelvl=getNodeLevel((*i));
      //             if(nodelvl==minlvl)
      //             minnode+=uniquecount[(*i)];
      //         }
      //        // for(int k=0;k<sizeunique;k++)
      //        // {
      //        //    if(incomingedgecount[uniqueNodes[k]]>0|| uniqueNodes[k]==e.getNode())
      //        //     map[uniqueNodes[k]]=uniqueNodes[k];
      //        // }
      //
      //
      //        // unsigned lvl=getNodeLevel(minIndex);
      //        // int sizeunique=unique->getNumEntries(lvl);
      //        // node_handle* uniqueNodes= new node_handle[sizeunique];
      //        // unique->getItems(getVarByLevel(lvl),uniqueNodes,sizeunique);
      //        //
      //        // for(int k=0;k<sizeunique;k++)
      //        // {
      //        //    if(incomingedgecount[uniqueNodes[k]]>0|| uniqueNodes[k]==e.getNode())
      //        //     map[uniqueNodes[k]]=uniqueNodes[k];
      //        // }
      //        // // printf("MinIndex %d\n",(minIndex+1) );
      //        // map[(minIndex)]=0;
      //        // delete uniqueNodes;
      //
      //
      //        // printf("MAke MAP! \n" );
      //
      //     // printf("getNodeStatus %d\n",getNodeStatus((minIndex)));
      //     // int deleted=uniquecount[minIndex];
      //     // printf("CALL RemoveDuplicate2 \n" );
      //      // for(auto i = levels.begin(); i != levels.end(); ++i){
      //      //     printf("%d , \n",(*i) );
      //      // }
      //      // getchar();
      //     ////////////////********************
      //
      //     // for(auto i = m.begin(); i != m.end(); ++i){
      //     //     printf("%d ,",(*i) );
      //     // }
      //     // printf("\n" );
      //     // getchar();
      //     // printf("*****\n" );
      //     // for (auto it = removedNode.begin(); it !=
      //     //                          removedNode.end(); ++it)
      //     // printf("%d\n",(*it) );
      //     // printf("*****\n" );
      //     // getchar();
      //     ////NotNeeded
      //     // std::set<int> s(removedNodeA);
      //    // s.insert(removedNodeB.begin(), removedNodeB.end());
      //    ///NotNeeded
      //    // for (auto it = removedNodeB.begin(); it !=
      //    //                          removedNodeB.end(); ++it)
      //    // // if(incomingedgecount[(*it)]==1)
      //    // printf("%d map %d iec %d lvl %d out of %d\n",(*it),map[(*it)],incomingedgecount[(*it)],getNodeLevel((*it)),getNumVariables() );
      //    // printf("**RNB***\n" );
      //    // getchar();
      //    #ifdef DBG_MTMDD
      //     printf("BEfore RNA %d RNB %d Union %d, minlvl %d\n",removedNodeA.size(),removedNodeB.size(),s.size(),minlvl );
      //     #endif
      //     printf("Cardinality before RD %lu\n", rootNumberofState/*belowcount[e.getNode()]*/);
      //     // std::list<int>* uniquelist=markNodesInSubgraphByLvl(&nl,1);
      //
      //     startc=clock();
      //     int merged=RemoveDuplicateSet(minlvl,levelarray,uniquelist,/*levels,*/map,e/*,removedNodeA,removedNodeB*/);
      //     printf("RemoveDuplicateSet time %f\n", double(clock() - startc) / double(CLOCKS_PER_SEC));
      //
      //
      //     // getchar();
      //     // if(belowcount!=0)delete[] belowcount; else {printf("ERRR belowcount\n"); getchar();}
      //     // maxid=e.getLastHandle();
      //     // lastNode=maxid+1;
      //     // apply(BC,e,cCard);
      //     //
      //     // printf("Cardinality after RD %ld\n", belowcount[e.getNode()]);
      //     // printf("root is %d\n",e.getNode() );
      //     printf("Node count is %d\n",e.getNodeCount() );
      //     if((deletedApproach==1)&&(e.getNodeCount()<=Threashold)){
      //     deletedApproach=0;
      //     }
      //     // std::set<int> result;
      //     // maxid=e.getLastHandle();
      //     // bool* visitedNodex= new bool[maxid];
      //     // for(int i=0;i<=(int)maxid;i++){
      //     //     visitedNodex[i]=false;
      //     // }
      //     // node_handle n1=e.getNode();
      //     // printf("Node LVL is %d & maxid is %d\n",getNodeLevel(n1),maxid );
      //     // getNC(getNodeLevel(n1),n1,visitedNodex,result);
      //     // delete[] visitedNodex;
      //     // printf("NCNCNC %d\n",result.size() );
      //     // getchar();
      //     #ifdef DBG_MTMDD
      //     printf("merged %d\n",merged );
      //     printf("After RNA %d RNB %d\n",removedNodeA.size(),removedNodeB.size() );
      //     // for (auto it = removedNodeB.begin(); it !=
      //     //                          removedNodeB.end(); ++it)
      //     // printf("%d map %d iec %d lvl %d out of %d\n",(*it),map[(*it)],incomingedgecount[(*it)],getNodeLevel((*it)),getNumVariables() );
      //     printf("**RNB***\n" );
      //     #endif
      //     // std::set<int> notRemoved;
      //     // int num_vars=getNumVariables();
      //     // for(int l=1;l<num_vars;l++){
      //     // int sizeunique=unique->getNumEntries(l);
      //     // node_handle* uniqueNodes= new node_handle[sizeunique];
      //     // unique->getItems(getVarByLevel(l),uniqueNodes,sizeunique);
      //     //  for(int k=0;k<sizeunique;k++){
      //     //      const bool is_in = removedNodeB.find(uniqueNodes[k]) != removedNodeB.end();
      //     //      if(is_in){
      //     //          notRemoved.insert(uniqueNodes[k]);
      //     //      }
      //     //  }
      //     // delete uniqueNodes;
      //     // }
      //     #ifdef DBG_MTMDD
      //     /*node_handle*/ n=e.getNode();
      //     /*node_handle**/ list = markNodesInSubgraph(&n, 1, false);
      //     int ix=0;
      //     for (long i=0; list[i]; i++) {
      //         // printf("%d\n",list[i] );
      //          const bool is_in = removedNodeB.find(list[i]) != removedNodeB.end();
      //          if(is_in)
      //         {   ix++;
      //         notRemoved.insert(list[i]);
      //         }
      //          // map[list[i]]=list[i];
      //      }
      //      delete[] list;
      //      printf("ix is %d\n",ix );
      //     for (auto it = notRemoved.begin(); it !=
      //                              notRemoved.end(); ++it)
      //     // if(incomingedgecount[(*it)]==1)
      //     {showNode(meddlyout,(*it), SHOW_DETAILS);
      //     // printf("%d***%d\n",(*it),incomingedgecount[(*it)] );
      //     // printf("%d map %d iec %d lvl %d out of %d\n",(*it),map[(*it)],incomingedgecount[(*it)],getNodeLevel((*it)),getNumVariables() );
      //     }
      //     printf("**notRemoved***\n" );
      //     #endif
      //     // if(ix!=0)
      //     // getchar();
      //     // double cI;
      //     // if(incomingedgecount!=0)delete [] incomingedgecount;else {printf("ERRR incomingedgecount\n"); getchar();}
      //     //
      //     // apply(IEC, e, cI);
      //     // getchar();
      //
      //     // cC=e.getNodeCount();
      //     // getchar();
      //
      //
      //     // for(int i=0;i<=lastNode;i++)
      //     // map[i]=0;
      //
      //
      //     // getchar();
      //     // printf("End CALL RemoveDuplicate2 \n" );
      //
      //     // cC--;
      //     // if(belowcount!=0)delete[] belowcount; else {printf("ERRR belowcount\n"); getchar();}
      //     if(incomingedgecount!=0)delete [] incomingedgecount;else {printf("ERRR incomingedgecount\n"); getchar();}
      //     // if(abovecount!=0)delete[] abovecount;else {printf("ERRR abovecount\n"); getchar();}
      //     highestunique.clear();
      //     if(highestuniquelist!=0)delete[] highestuniquelist;else {printf("ERRR highestuniquelist\n"); getchar();}
      //     if(uniquecount!=0)delete[] uniquecount;else {printf("ERRR uniquecount\n"); getchar();}
      //     lowestunique.clear();
      //     if(lowestuniquelist!=0)delete[] lowestuniquelist;else {printf("ERRR lowestuniquelist\n"); getchar();}
      //     if(uniqueAbovecount!=0) delete[] uniqueAbovecount; else{printf("ERR uniqueAbovecount\n" );getchar();}
      //     // printf("DELETED \n" );
      //     delete[] ACBC;
      //     delete[] levelarray;
      //     delete[] list;
      //     map.clear();
      //     }
      //     // e.show(meddlyout,0);
      //     #ifdef DBG_MTMDD
      //     printf("HNodeCount is %d\n",cC );
      //     getchar();
      //     #endif
      //     // }
      //         /////////////*******
      //     // printf("NddeCountAfter %d\n",e.getNodeCount() );
      //     end = clock();
      //     double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      //     printf("Time taken %f \n",time_taken );
      //     printf("HunderApproximate End\n" );
      // }

      // //////////////////////////////////////
 int MEDDLY::mt_mdd_bool::RemoveDuplicateSet(int lvl, bool* levelarray,std::list<int>* uniquelist/*, std::set<int>levels*/,std::map<int,int> map,dd_edge &e/*,std::set<int>&RNA,std::set<int>&RNB*/){
     int root=e.getNode();
     // printf("RemoveDuplicateSet\n" );
      // printf("size is %d\n",e.getNodeCount() );
     bool changeInLevel=false;
     int num_vars=getNumVariables();
     FILE_output meddlyout(stdout);
     // std::map<int,int> nmap;
     // printf("RemoveDuplicateSet\n" );
     int merged=0;
     int nodeChangedNumber=0;
     // node_handle* uniqueNodes;
     for(int l=lvl;l<num_vars;l++){
         #ifdef DBG_MTMDD
         printf("l is %d\n",l );
         // printf("RNA %d RNB %d\n",RNA.size(),RNB.size() );
 #endif
         // getchar();
         bool is_in =levelarray[l];// levels.find(l) != levels.end();
         // printf("is_in is %d\n",is_in );
         // getchar();
         if(is_in){
             // levels.erase(l);
             // for(auto i = levels.begin(); i != levels.end(); ++i){
             //     printf("%d , \n",(*i) );
             // }
         }
         changeInLevel=false;

         // int sizeunique=unique->getNumEntries(l+1);
         // uniqueNodes= new node_handle[sizeunique];
         // clock_t startc=clock();
         // unique->getItems(getVarByLevel(l+1),uniqueNodes,sizeunique);
         // printf("Time taken unique->getItems %f \n",double(clock() - startc) / double(CLOCKS_PER_SEC) );
         // clock_t startc=clock();
         // uniquelist[l+1];
         // printf("Time taken uniquelist[l] %f \n",double(clock() - startc) / double(CLOCKS_PER_SEC) );
         // printf("Size %d, Size %d\n", sizeunique,uniquelist[l+1].size());
         // getchar();
         int uniquelistsize=uniquelist[l+1].size();
         int b=getDomain()->getVariableBound(l+1, false);
        // clock_t startd=clock();
         // for(int k=0;k<sizeunique;k++){
         // if(((incomingedgecount[uniqueNodes[k]]>0)||(uniqueNodes[k]==root))&&(uniqueNodes[k]<=maxid)&&(map[uniqueNodes[k]]!=0)){
         for(const auto& unique : uniquelist[l+1]){
             // if(((ACBC/*incomingedgecount*/[unique]>0)||(unique==root))&&(unique<=maxid)&&(map[unique]!=0)){
             if(map[unique]!=0){

             unpacked_node* un =  unpacked_node::newFull(this, l+1,b);
             bool nodeChanged=false;
             // clock_t starte=clock();
             for(int ik=0;ik<b; ik++){
                 // int dpt=getDownPtr(uniqueNodes[k],ik);
                 int dpt=getDownPtr(unique,ik);

                 // if(dpt==3463){
                 //     printf("Parent 3463 is %d\n",uniqueNodes[k] );
                 //     getchar();
                 // }
                 if(dpt>=1 && dpt!=map[dpt]){
                     // if(map[dpt]!=0)
                     // showNode(meddlyout,uniqueNodes[k], SHOW_DETAILS);
                     #ifdef DBG_MTMDD
                     printf("C uniqueNodes %d\n",dpt );
             #endif
                     ///////////////Remove RNA RNB
                     // RNA.erase(dpt);
                     // RNB.erase(dpt);
                     ///////////////Remove RNA RNB

                     // printf("%d\n",uniqueNodes[k] );
                     un->d_ref(ik) = this->linkNode(map[dpt]);
                     nodeChanged=true;

                 }
                 else{
                     un->d_ref(ik)=this->linkNode(dpt);
                 }
             }
             // printf("Time taken check node all pointers %f \n",double(clock() - starte) / double(CLOCKS_PER_SEC) );

             if(nodeChanged){
                 // clock_t startk=clock();
                 nodeChangedNumber++;
                 ///////////////Remove RNA RNB

                 // RNA.erase(uniqueNodes[k]);
                 // RNB.erase(uniqueNodes[k]);
                 ///////////////Remove RNA RNB

                 // showNode(meddlyout,uniqueNodes[k], SHOW_DETAILS);
                 // printf("\n" );
#ifdef DBG_MTMDD
                 printf("D uniqueNodes %d",uniqueNodes[k] );
#endif
                 changeInLevel=true;
                 node_handle q=0;//=unique->find(*un, getVarByLevel(un->getLevel()));
                 if(q!=0){
#ifdef DBG_MTMDD
                     printf("merged\n" );
#endif
                     merged++;
                     // printf("Found\n" );
                     // printf("Dup %d\n",uniqueNodes[k] );
                     // map[uniqueNodes[k]]=q;
                     map[unique]=q;
                     unpacked_node::recycle(un);
                    }else{

                        if(l+1<num_vars){
                             // clock_t startr=clock();
                            node_handle temporaryNode=this->createReducedNode(-1,un);
                            // printf("Time taken createReducedNode %f \n",double(clock() - startr) / double(CLOCKS_PER_SEC) );

#ifdef DBG_MTMDD
                            showNode(meddlyout, uniqueNodes[k], SHOW_DETAILS);
                            showNode(meddlyout, temporaryNode, SHOW_DETAILS);

                             printf("newNode%d\n",temporaryNode );
#endif
                             ///////////////Remove RNA RNB

                             // RNA.erase(temporaryNode);
                             // RNB.erase(temporaryNode);
                             ///////////////Remove RNA RNB

                            // printf("%d\n",uniqueNodes[k] );
                            this->linkNode(temporaryNode);
                              // showNode(meddlyout, temporaryNode, SHOW_DETAILS);
                             // map[uniqueNodes[k]]=temporaryNode;
                              map[unique]=temporaryNode;
                              // printf("DONE MODIFY IN PLACE %d to %d\n",uniqueNodes[k],temporaryNode );
                              // getchar();


                         }else{
#ifdef DBG_MTMDD
                             printf("HERE WANT TO set root\n");
#endif
                             node_handle node=this->createReducedNode(-1,un);
                             // printf("%d\n",uniqueNodes[k] );
                             node = makeNodeAtTop(node);
                             this->linkNode(node);

                             e.set(node);
#ifdef DBG_MTMDD
                             showNode(meddlyout, node, SHOW_DETAILS);
                             printf("Root set! %d\n",node);
                             getchar();
#endif
                             // printf("end size is %d\n",e.getNodeCount() );
                             //nmap.clear();
                             map.clear();
                             // changeInLevel=false;
#ifdef DBG_MTMDD
                             printf("NodeChangedNumber %d\n",nodeChangedNumber );
#endif
                             return merged;
                         }
                     }
                     // printf("Time taken node changed %f \n",double(clock() - startk) / double(CLOCKS_PER_SEC) );

                 }else{
                 // map[uniqueNodes[k]]=uniqueNodes[k];
                 unpacked_node::recycle(un);
             }
         }
         else{
             // if((map[uniqueNodes[k]]==0)&&(incomingedgecount[uniqueNodes[k]]>0))
#ifdef HAVE_LIBGMP
             mpz_object mpzd;
             mpz_object mpz0;
             mpz0.setValue(0);
             mpz0.setReminder(0);
             if((map[unique]==0)&&(mpzd.compare(ACBC[unique],mpz0)==1))
#else
             if((map[unique]==0)&&(ACBC[unique]>0))
#endif
             // printf("MAP %d is 0\n",uniqueNodes[k] );
             {
                 // showNode(meddlyout,uniqueNodes[k], SHOW_DETAILS);
                 // printf("\n" );
                 #ifdef DBG_MTMDD
                 printf("0 uniqueNodes %d\n",uniqueNodes[k] );
         #endif
                 ///////////////Remove RNA RNB

                 // RNA.erase(uniqueNodes[k]);
                 // RNB.erase(uniqueNodes[k]);
                 ///////////////Remove RNA RNB

             }
             ;//nothing
         }
        }
        // printf("Time taken forloop unique %f \n",double(clock() - startd) / double(CLOCKS_PER_SEC) );

        // delete[] uniqueNodes;
        #ifdef DBG_MTMDD
        printf("changeInLevel %d\n",changeInLevel );
#endif
        if(changeInLevel==false){
            bool correct=false;
            for(int k=0;k<=num_vars; k++)
                if (levelarray[k]==true)
                    correct=true;
            // if(levels.size()==0)
            if(correct==false)
            {
                printf("NOT CORRECT\n" );
                exit(0);
                // getchar();
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
     printf("Done RemoveDuplicateSet\n" );

return merged;
#ifdef DBG_MTMDD
printf("Done RemoveDuplicateSet\n" );
#endif
 }
int MEDDLY::mt_mdd_bool::RemoveDuplicate2(int lvl, std::map<int,int> map,dd_edge &e,std::set<int>&RNA,std::set<int>&RNB){
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
            unpacked_node* un =  unpacked_node::newFull(this, l+1, b);//unpacked_node::newFull(this, l+1, b);
            bool nodeChanged=false;
            for(int ik=0;ik<b; ik++){
                int dpt=getDownPtr(uniqueNodes[k],ik);
                if(dpt>=1 && dpt!=map[dpt]){
                    RNA.erase(dpt);
                    RNB.erase(dpt);
                    // printf("\nChanged %d th ptr which was %d to %d and lvl is %d \n", ik, dpt,map[dpt], getNodeLevel(dpt) );
                    un->d_ref(ik) = this->linkNode(map[dpt]);
                    nodeChanged=true;
                }
                else{
                    un->d_ref(ik)=this->linkNode(dpt);
                }
            }
            if(nodeChanged){
                RNA.erase(uniqueNodes[k]);
                RNB.erase(uniqueNodes[k]);
                // showNode(meddlyout,uniqueNodes[k], SHOW_DETAILS);
                changeInLevel=true;
            node_handle q=0;//=unique->find(*un, getVarByLevel(un->getLevel()));
             // printf("p %d is %d\n",uniqueNodes[k],q );
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


// bool MEDDLY::mt_mdd_bool::compareLong(mpz_object a, mpz_object b)
// {
//     return true;
//     // return (i1.start < i2.start);
// }
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

           for(int ik=0;ik<b; ik++){
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
                    // printf("Changed %d th ptr which was %d to %d\n", ik, dpt,map[dpt] );
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

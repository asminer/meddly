
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
#include "config.h"
#endif
#include "../defines.h"
#include "reach_bfs.h"
#include <list>
#include <algorithm>
// #define DEBUG_BFS
// #define VERBOSE_BFS
// #define BFSDFS
#define Normal
// #define CHKBFS
// #define maxThreshold 100000
// #define minThreshold 100000
// #define iterationcount 50
// #define desiredPercentage 0.5
// #define DUP

namespace MEDDLY {
    long maxThreshold;
    long minThreshold;
    float desiredPercentage;
    unsigned int optionsForUA;
    float rootStatePercentage;
    int deletedApproach;
    long timeForUA;
    long timeT;
  class common_bfs;
  class common_bfs_ua;
  class common_bfs_hua;
  class common_bfs_event;
  // class common_bfs_mt;
  class forwd_bfs_mt;
  class bckwd_bfs_mt;
  class forwd_bfs_ua_mt;
  class forwd_bfs_hua_mt;
  class forwd_bfs_mt_event;


  // class common_bfs_evplus;
  class forwd_bfs_evplus;
  class bckwd_bfs_evplus;
  class forwd_bfs_ua_evplus;
  class forwd_bfs_hua_evplus;

  class forwd_bfs_opname;
  class bckwd_bfs_opname;
  class forwd_bfs_ua_opname;
  class forwd_bfs_hua_opname;
  class forwd_bfs_opname_event;
};

// ******************************************************************
// *                                                                *
// *                        common_bfs class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::common_bfs : public binary_operation {
  public:
    common_bfs(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);

  protected:
    inline void setUnionOp(binary_operation* uop)
    {
      MEDDLY_DCASSERT(uop);
      MEDDLY_DCASSERT(0==unionOp);
      unionOp = uop;
    }

    inline void setImageOp(binary_operation* iop)
    {
      MEDDLY_DCASSERT(iop);
      MEDDLY_DCASSERT(0==imageOp);
      imageOp = iop;
    }

  private:
    binary_operation* unionOp;
    binary_operation* imageOp;

};


MEDDLY::common_bfs::common_bfs(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc, 0, a1, a2, res)
{
  unionOp = 0;
  imageOp = 0;
}

void MEDDLY::common_bfs::computeDDEdge(const dd_edge &init, const dd_edge &R, dd_edge &reachableStates, bool userFlag)
{
    clock_t start, end;
    start = clock();
  MEDDLY_DCASSERT(unionOp);
  MEDDLY_DCASSERT(imageOp);

  reachableStates = init;
  dd_edge prevReachable(resF);
  dd_edge front(resF);
  FILE_output meddlyout(stdout);
  long peakreachable=0;
  long lastreachable=0;
#ifdef DEBUG_BFS
  FILE_output debug(stderr);
  debug << "Relation: ";
  R.show(debug, 2);
  debug << "Initial states: ";
  init.show(debug, 2);
  long iters = 0;
#endif
#ifdef VERBOSE_BFS
  long iters = 0;
  FILE_OUTPUT verbose(stderr);
#endif
  while (prevReachable != reachableStates) {
      // printf("XXXXXIt Done0000\n" );
      // R.show(meddlyout,0);
#ifdef VERBOSE_BFS
    iters++;
    verbose << "Iteration " << iters << ":\n";
#endif
    // printf("XXXXXIt Done0\n" );
    prevReachable = reachableStates;
    // printf("XXXXXIt Done 1\n" );

     // reachableStates.show(meddlyout, 0);

    imageOp->computeDDEdge(reachableStates, R, front, userFlag);
    // printf("XXXXXIt Done2\n" );

#ifdef VERBOSE_BFS
    verbose << "\timage done ";
    front.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    iters++;
    debug << "Iteration " << iters << "\npseudo-frontier: ";
    front.show(debug, 2);
#endif
    unionOp->computeDDEdge(reachableStates, front, reachableStates, userFlag);
#ifdef VERBOSE_BFS
    verbose << "\tunion done ";
    reachableStates.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    debug << "Reachable so far: ";
    reachableStates.show(debug, 2);
#endif
printf("cardafter %f\n",reachableStates.getCardinality() );
lastreachable=reachableStates.getNodeCount();
#ifdef CHKBFS
if(lastreachable>maxThreshold)
return;
#endif
if(lastreachable>peakreachable)
{peakreachable=lastreachable;}
printf("XXXX %ld\t %ld\n",lastreachable, peakreachable );

 // reachableStates.getForest()->underApproximate(reachableStates,1100);
//1000
// reachableStates.show(meddlyout,0);
 // printf("XXXX\n" );
 end = clock();
 double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
 if(time_taken>timeT){
 printf("TimeOut\n" );
 return;
}
  }




}

// ******************************************************************
// *                                                                *
// *                        common_bfs_ua class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::common_bfs_ua : public binary_operation {
  public:
    common_bfs_ua(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);

  protected:
    inline void setUnionOp(binary_operation* uop)
    {
      MEDDLY_DCASSERT(uop);
      MEDDLY_DCASSERT(0==unionOp);
      unionOp = uop;
    }

    inline void setImageOp(binary_operation* iop)
    {
      MEDDLY_DCASSERT(iop);
      MEDDLY_DCASSERT(0==imageOp);
      imageOp = iop;
    }

  private:
    binary_operation* unionOp;
    binary_operation* imageOp;

};


MEDDLY::common_bfs_ua::common_bfs_ua(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc, 0, a1, a2, res)
{
  unionOp = 0;
  imageOp = 0;
}

void MEDDLY::common_bfs_ua::computeDDEdge(const dd_edge &init, const dd_edge &R, dd_edge &reachableStates, bool userFlag)
{
    #ifdef BFSDFS
    MEDDLY_DCASSERT(unionOp);
    MEDDLY_DCASSERT(imageOp);
    reachableStates = init;
    binary_operation* opMinus = getOperation(DIFFERENCE, reachableStates, reachableStates, reachableStates);
    MEDDLY_DCASSERT(opMinus);

    dd_edge from(resF);
    from=init;
    dd_edge to(resF);
    dd_edge newstate(resF);
    dd_edge reachablenew(resF);
    // dd_edge front(resF);
    FILE_output meddlyout(stdout);
    // dd_edge* arrddedge= new dd_edge[10];
    // int k=0;
    // int size=0;
    #ifdef DEBUG_BFS
    FILE_output debug(stderr);
    debug << "Relation: ";
    R.show(debug, 2);
    debug << "Initial states: ";
    init.show(debug, 2);
    long iters = 0;
    #endif
    #ifdef VERBOSE_BFS
    long iters = 0;
    FILE_OUTPUT verbose(stderr);
    #endif
    while (true) {
        printf("XXXXXIt Done0000\n" );
        // R.show(meddlyout,0);
    #ifdef VERBOSE_BFS
      iters++;
      verbose << "Iteration " << iters << ":\n";
    #endif
      // printf("XXXXXIt Done0\n" );
      // prevReachable = reachableStates;
      // printf("XXXXXIt Done 1\n" );

       // reachableStates.show(meddlyout, 0);

      imageOp->computeDDEdge(from, R, to, userFlag);
      // if(front.getNodeCount()>size)
      // size=front.getNodeCount();
      // printf("XXXX %d %d\n",front.getNodeCount(),size );
      // to.getForest()->underApproximate(to,1500,15000);

      printf("XXXXXIt Done2\n" );

    #ifdef VERBOSE_BFS
      verbose << "\timage done ";
      front.show(verbose, 0);
      verbose << "\n";
    #endif
    #ifdef DEBUG_BFS
      iters++;
      debug << "Iteration " << iters << "\npseudo-frontier: ";
      front.show(debug, 2);
    #endif
    opMinus->computeTemp(to,reachableStates,newstate);
    printf("XXXXXIt Done3\n" );

    if(to==reachableStates){//(newstate.getNodeCount()==0){
        printf("XXXXXIt Done3-1\n" );

        unionOp->computeDDEdge(reachableStates, newstate, reachablenew, userFlag);

        imageOp->computeDDEdge(reachablenew, R, to, userFlag);
        if(to==reachableStates){ return;}
        else{
            printf("XXXXXIt Done3-2\n" );
            opMinus->computeTemp(to,reachableStates,newstate);
            printf("XXXXXIt Done3-3\n" );
        }
    }
    printf("XXXXXIt Done4\n" );

    from=newstate;
    printf("XXXXXIt Done5\n" );

    from.getForest()->underApproximate(from,5000,18000,0,0);

      unionOp->computeDDEdge(reachableStates, from, reachableStates, userFlag);
      // printf("XXXXXIt Done3\n" );
    #ifdef VERBOSE_BFS
      verbose << "\tunion done ";
      reachableStates.show(verbose, 0);
      verbose << "\n";
    #endif
    #ifdef DEBUG_BFS
      debug << "Reachable so far: ";
      reachableStates.show(debug, 2);
    #endif
    // expert_forest* ef=(expert_forest*)reachableStates.getForest();
    // printf("XXXX %d\n",reachableStates.getNodeCount() );
    //2389790
    //getchar();
    // reachableStates.getForest()->underApproximate(reachableStates,1000,1500);
    // for(int i=0;i<10;i++){
    //   if(reachableStates==arrddedge[i])
    //   {
    //       printf("duplicate\n" );
    //       getchar();
    //       getchar();
    //   }
    // }
    // arrddedge[k%10]=reachableStates;
    // k++;
    //1063000
    // reachableStates.getForest()->underApproximate(reachableStates, 500,500);
    //1000
    // reachableStates.show(meddlyout,0);
    // printf("XXXX\n" );
    // getchar();
    }

    // delete [] arrddedge;
#endif

#ifdef Normal
clock_t start, end;
start = clock();
  MEDDLY_DCASSERT(unionOp);
  MEDDLY_DCASSERT(imageOp);
  // binary_operation* opMinus = getOperation(DIFFERENCE, reachableStates, reachableStates, reachableStates);
  reachableStates = init;
  dd_edge prevReachable(resF);
  dd_edge front(resF);
  FILE_output meddlyout(stdout);
    #ifdef DUP
  dd_edge* arrddedge= new dd_edge[10];
  int k=0;
  #endif
  #ifdef iterationcount
  int uacall=0;
  #endif
  int step=0;
  int size=0;
#ifdef DEBUG_BFS
  FILE_output debug(stderr);
  debug << "Relation: ";
  R.show(debug, 2);
  debug << "Initial states: ";
  init.show(debug, 2);
  long iters = 0;
#endif
#ifdef VERBOSE_BFS
  long iters = 0;
  FILE_OUTPUT verbose(stderr);
#endif
  while (prevReachable != reachableStates) {
      // printf("XXXXXIt Done0000\n" );
      // R.show(meddlyout,0);
#ifdef VERBOSE_BFS
    iters++;
    verbose << "Iteration " << iters << ":\n";
#endif
    printf("XXXXXIt Done0\n" );
    prevReachable = reachableStates;
    // printf("XXXXXIt Done 1\n" );

     // reachableStates.show(meddlyout, 0);

    imageOp->computeDDEdge(reachableStates, R, front, userFlag);

    step++;
    // if(front.getNodeCount()>size)
    // size=front.getNodeCount();
    // printf("XXXX %d %d %d\n",step, front.getNodeCount(),size );
    // front.getForest()->underApproximate(front,20,200,0);
    // printf("XXXX %d\n",step);
    // printf("XXXXXIt Done2\n" );

#ifdef VERBOSE_BFS
    verbose << "\timage done ";
    front.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    iters++;
    debug << "Iteration " << iters << "\npseudo-frontier: ";
    front.show(debug, 2);
#endif
    unionOp->computeDDEdge(reachableStates, front, reachableStates, userFlag);
    printf("before %d\n",reachableStates.getNodeCount() );
    printf("cardbefore %f\n",reachableStates.getCardinality() );
    #ifdef iterationcount
    if(reachableStates.getNodeCount()>maxThreshold){
    if(uacall<iterationcount)
    uacall++;
    else
    return;
    }
    #endif
    // reachableStates.getForest()->underApproximate(reachableStates,minThreshold,maxThreshold,0,0);
    reachableStates.getForest()->underApproximate(reachableStates,minThreshold,maxThreshold,desiredPercentage,optionsForUA);

    // reachableStates.getForest()->HeuristicUnderApproximate(reachableStates,minThreshold,maxThreshold,0);
    printf("after %d\n",reachableStates.getNodeCount() );
    printf("cardafter %f\n",reachableStates.getCardinality() );
    printf("XXXX %d\n",step);
    unionOp->computeDDEdge(reachableStates, init, reachableStates, userFlag);

    // printf("XXXXXIt Done3\n" );
#ifdef VERBOSE_BFS
    verbose << "\tunion done ";
    reachableStates.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    debug << "Reachable so far: ";
    reachableStates.show(debug, 2);
#endif
// expert_forest* ef=(expert_forest*)reachableStates.getForest();
// printf("XXXX %d\n",reachableStates.getNodeCount() );
//2389790
//getchar();
// reachableStates.getForest()->underApproximate(reachableStates,1000,1500);
#ifdef DUP
for(int i=0;i<10;i++){
    if(reachableStates==arrddedge[i])
    {
        printf("duplicate\n" );
        getchar();
        getchar();
    }
}
arrddedge[k%10]=reachableStates;
k++;
#endif
//1063000
// reachableStates.getForest()->underApproximate(reachableStates, 500,500);
//1000
// reachableStates.show(meddlyout,0);
 // printf("XXXX\n" );
// getchar();
end = clock();
double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
if(time_taken>timeForUA){
printf("TimeOut\n" );
return;
}
  }
  #ifdef DUP
delete [] arrddedge;
#endif
#endif
}



// ******************************************************************
// *                                                                *
// *                        common_bfs_hua class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::common_bfs_hua : public binary_operation {
  public:
    common_bfs_hua(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);

  protected:
    inline void setUnionOp(binary_operation* uop)
    {
      MEDDLY_DCASSERT(uop);
      MEDDLY_DCASSERT(0==unionOp);
      unionOp = uop;
    }

    inline void setImageOp(binary_operation* iop)
    {
      MEDDLY_DCASSERT(iop);
      MEDDLY_DCASSERT(0==imageOp);
      imageOp = iop;
    }

  private:
    binary_operation* unionOp;
    binary_operation* imageOp;

};


MEDDLY::common_bfs_hua::common_bfs_hua(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc, 0, a1, a2, res)
{
  unionOp = 0;
  imageOp = 0;
}

void MEDDLY::common_bfs_hua::computeDDEdge(const dd_edge &init, const dd_edge &R, dd_edge &reachableStates, bool userFlag)
{
    #ifdef BFSDFS
    MEDDLY_DCASSERT(unionOp);
    MEDDLY_DCASSERT(imageOp);
    reachableStates = init;
    binary_operation* opMinus = getOperation(DIFFERENCE, reachableStates, reachableStates, reachableStates);
    MEDDLY_DCASSERT(opMinus);

    dd_edge from(resF);
    from=init;
    dd_edge to(resF);
    dd_edge newstate(resF);
    dd_edge reachablenew(resF);
    // dd_edge front(resF);
    FILE_output meddlyout(stdout);
    // dd_edge* arrddedge= new dd_edge[10];
    // int k=0;
    // int size=0;
    #ifdef DEBUG_BFS
    FILE_output debug(stderr);
    debug << "Relation: ";
    R.show(debug, 2);
    debug << "Initial states: ";
    init.show(debug, 2);
    long iters = 0;
    #endif
    #ifdef VERBOSE_BFS
    long iters = 0;
    FILE_OUTPUT verbose(stderr);
    #endif
    while (true) {
        printf("XXXXXIt Done0000\n" );
        // R.show(meddlyout,0);
    #ifdef VERBOSE_BFS
      iters++;
      verbose << "Iteration " << iters << ":\n";
    #endif
      // printf("XXXXXIt Done0\n" );
      // prevReachable = reachableStates;
      // printf("XXXXXIt Done 1\n" );

       // reachableStates.show(meddlyout, 0);

      imageOp->computeDDEdge(from, R, to, userFlag);
      // if(front.getNodeCount()>size)
      // size=front.getNodeCount();
      // printf("XXXX %d %d\n",front.getNodeCount(),size );
      // to.getForest()->underApproximate(to,1500,15000);

      printf("XXXXXIt Done2\n" );

    #ifdef VERBOSE_BFS
      verbose << "\timage done ";
      front.show(verbose, 0);
      verbose << "\n";
    #endif
    #ifdef DEBUG_BFS
      iters++;
      debug << "Iteration " << iters << "\npseudo-frontier: ";
      front.show(debug, 2);
    #endif
    opMinus->computeTemp(to,reachableStates,newstate);
    printf("XXXXXIt Done3\n" );

    if(to==reachableStates){//(newstate.getNodeCount()==0){
        printf("XXXXXIt Done3-1\n" );

        unionOp->computeDDEdge(reachableStates, newstate, reachablenew, userFlag);

        imageOp->computeDDEdge(reachablenew, R, to, userFlag);
        if(to==reachableStates){ printf("Done to=reachableStates\n" );return;}
        else{
            printf("XXXXXIt Done3-2\n" );
            opMinus->computeTemp(to,reachableStates,newstate);
            printf("XXXXXIt Done3-3\n" );
        }
    }
    printf("XXXXXIt Done4\n" );

    from=newstate;
    printf("XXXXXIt Done5\n" );

    from.getForest()->underApproximate(from,5000,18000);

      unionOp->computeDDEdge(reachableStates, from, reachableStates, userFlag);
      // printf("XXXXXIt Done3\n" );
    #ifdef VERBOSE_BFS
      verbose << "\tunion done ";
      reachableStates.show(verbose, 0);
      verbose << "\n";
    #endif
    #ifdef DEBUG_BFS
      debug << "Reachable so far: ";
      reachableStates.show(debug, 2);
    #endif
    // expert_forest* ef=(expert_forest*)reachableStates.getForest();
    // printf("XXXX %d\n",reachableStates.getNodeCount() );
    //2389790
    //getchar();
    // reachableStates.getForest()->underApproximate(reachableStates,1000,1500);
    // for(int i=0;i<10;i++){
    //   if(reachableStates==arrddedge[i])
    //   {
    //       printf("duplicate\n" );
    //       getchar();
    //       getchar();
    //   }
    // }
    // arrddedge[k%10]=reachableStates;
    // k++;
    //1063000
    // reachableStates.getForest()->underApproximate(reachableStates, 500,500);
    //1000
    // reachableStates.show(meddlyout,0);
    // printf("XXXX\n" );
    // getchar();
    }

    // delete [] arrddedge;
#endif

#ifdef Normal
clock_t start, end;
start = clock();
  MEDDLY_DCASSERT(unionOp);
  MEDDLY_DCASSERT(imageOp);
  // binary_operation* opMinus = getOperation(DIFFERENCE, reachableStates, reachableStates, reachableStates);
  reachableStates = init;
  dd_edge prevReachable(resF);
  dd_edge front(resF);
  FILE_output meddlyout(stdout);
    #ifdef DUP
  dd_edge* arrddedge= new dd_edge[10];
  int k=0;
  #endif
  #ifdef iterationcount
  int uacall=0;
  #endif
  int step=0;
  int size=0;
#ifdef DEBUG_BFS
  FILE_output debug(stderr);
  debug << "Relation: ";
  R.show(debug, 2);
  debug << "Initial states: ";
  init.show(debug, 2);
  long iters = 0;
#endif
#ifdef VERBOSE_BFS
  long iters = 0;
  FILE_OUTPUT verbose(stderr);
#endif
std::list<std::pair<long int ,int>> PairStateNode;
std::pair<long int ,int> lastBFSStateNode;
int BFSbetween=0;
//std::list<long int> stateCount;
//std::list<int> nodeCount;
bool oneMoreBFS=false;
bool BFSExec=false;
  while (prevReachable != reachableStates) {
      BFSExec=false;
      // printf("XXXXXIt Done0000\n" );
      // R.show(meddlyout,0);
#ifdef VERBOSE_BFS
    iters++;
    verbose << "Iteration " << iters << ":\n";
#endif
    printf("XXXXXIt Done0\n" );
    prevReachable = reachableStates;
    // printf("XXXXXIt Done 1\n" );

     // reachableStates.show(meddlyout, 0);

    imageOp->computeDDEdge(reachableStates, R, front, userFlag);

    step++;
    // if(front.getNodeCount()>size)
    // size=front.getNodeCount();
    // printf("XXXX %d %d %d\n",step, front.getNodeCount(),size );
    // front.getForest()->underApproximate(front,20,200,0);
    // printf("XXXX %d\n",step);
    // printf("XXXXXIt Done2\n" );

#ifdef VERBOSE_BFS
    verbose << "\timage done ";
    front.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    iters++;
    debug << "Iteration " << iters << "\npseudo-frontier: ";
    front.show(debug, 2);
#endif
    unionOp->computeDDEdge(reachableStates, front, reachableStates, userFlag);
    BFSExec=true;
    int currentNC=reachableStates.getNodeCount();
    printf("before %d\n",reachableStates.getNodeCount() );
    double currentRS=reachableStates.getCardinality();
    printf("cardbefore %f\n",currentRS );

    #ifdef iterationcount
    if(reachableStates.getNodeCount()>maxThreshold){
    if(uacall<iterationcount)
    uacall++;
    else
    return;
    }
    #endif

    // reachableStates.getForest()->underApproximate(reachableStates,minThreshold,maxThreshold,0);
    /*if(reachableStates.getNodeCount()>maxThreshold && updatePercentage){
        if(lastRS>0 && pprercentage>0){
            if(lastRS<currentRS){;}
            else if(lastRS==currentRS){
                if(lastNC<currentNC){
                    desiredPercentage-=pprercentage;
                }
                if(lastNC>currentNC){
                    desiredPercentage+=pprercentage;
                }
                if(lastNC==currentNC){
                     printf("this approach might not works!\n" );
                     getchar();
                }

            }else if(lastRS>currentRS){
                desiredPercentage-=pprercentage;
            }

        }
        pprercentage=currentRS/lastRS;
        lastRS=currentRS;
        lastNC=currentNC;

    }*/
    // reachableStates.getForest()->HeuristicUnderApproximate(reachableStates,minThreshold,maxThreshold,desiredPercentage,optionsForUA,deletedApproach, rootStatePercentage);
    long int getcard=reachableStates.getCardinality();
      int getnode=reachableStates.getNodeCount();
      if(oneMoreBFS){
          if(lastBFSStateNode.first<getcard){
              printf("More state!\n" );
              if(lastBFSStateNode.second<getnode){
                  printf("More node!" );
                  return;
                  // getchar();
              }else{
                  printf("LESS node!->good" );
                  // getchar();
              }
          }
          else{
              printf("LESS state!\n" );
              // getchar();
          }
      }
      else{
      auto it=std::find(PairStateNode.begin(),PairStateNode.end(),std::make_pair(getcard,getnode));
      if(it!=PairStateNode.end()){
          printf("Might be stuck!\n" );
          oneMoreBFS=true;
          lastBFSStateNode=std::make_pair(getcard,getnode);
          // getchar();
      }
      }
      /*auto it = std::find(stateCount.begin(), stateCount.end(), getcard);
      if(it != stateCount.end()){
          // int dist=std::distance(stateCount.begin(), it);
          auto itn = std::find(nodeCount.begin(), nodeCount.end(),getnode);
          if(itn != nodeCount.end()){
          // if(nodeCount[dist]==getno

              std::srand(time(0)); // use current time as seed for random generator
              float random_variable = std::rand()%100;
              float random_variable2 = std::rand()%100;
              printf("Might be stuck! %f %f\n",random_variable,random_variable2 );
              getchar();
              if(random_variable2>50){
                  desiredPercentage-=0.1;
              }
              else{
                  if(random_variable>50)
                  minThreshold=(int)(0.90*minThreshold);
                  else{
                      minThreshold=(int)(0.95*minThreshold);
                  }
              }
          }
      }
      */
      if(getnode >maxThreshold && oneMoreBFS==false ){
          PairStateNode.push_back(std::make_pair(getcard,getnode));
          // stateCount.push_back(getcard);
          // nodeCount.push_back(getnode);
      }

      // reachableStates.getForest()->underApproximate(reachableStates,minThreshold,maxThreshold,0);
      if(getnode>maxThreshold){
         BFSExec=false;
        }
      if(oneMoreBFS==false)
      {
          long originalminThreshold=minThreshold;
          // if(BFSbetween<=1){
          //     minThreshold=(int)(0.90*minThreshold);
          //     desiredPercentage-=0.1;
          //     if(desiredPercentage<0.0) desiredPercentage=0.1;
          //     maxThreshold=maxThreshold+ (int)(0.01*maxThreshold);
          //     printf("%d %d, %f\n", minThreshold,maxThreshold,desiredPercentage);
          // }
          // if(getnode<maxThreshold){
          //     BFSbetween++;
          // }
          reachableStates.getForest()->HeuristicUnderApproximate(reachableStates,minThreshold,maxThreshold,desiredPercentage,optionsForUA,deletedApproach, rootStatePercentage);
          minThreshold=originalminThreshold;

          // if(getnode<maxThreshold){
          //     BFSbetween++;
          // }else
          BFSbetween=0;
      }else{
          BFSbetween++;
      }

      /*else{
          oneMoreBFS=false;
      }*/
    // Option 0 is minThreshold and maxThreshold
    // Option 2 is maxThreshold and remove until all nodes with density < rootdensity*desiredPercentage is deleted
    // Option 3 is maxThreshold and remove until
    //  1-all nodes with density <rootdensity*desiredPercentage is deleted, and
    //  2-number of nodes<maxThreshold
    printf("after %d\n",reachableStates.getNodeCount() );
    // getchar();
    printf("cardafter %f\n",reachableStates.getCardinality() );
    printf("XXXX %d\n",step);
    unionOp->computeDDEdge(reachableStates, init, reachableStates, userFlag);
    if(reachableStates==prevReachable){
        printf("Equal %d\n", BFSExec);
    }
    // printf("XXXXXIt Done3\n" );
#ifdef VERBOSE_BFS
    verbose << "\tunion done ";
    reachableStates.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    debug << "Reachable so far: ";
    reachableStates.show(debug, 2);
#endif
// expert_forest* ef=(expert_forest*)reachableStates.getForest();
// printf("XXXX %d\n",reachableStates.getNodeCount() );
//2389790
//getchar();
// reachableStates.getForest()->underApproximate(reachableStates,1000,1500);
#ifdef DUP
for(int i=0;i<10;i++){
    if(reachableStates==arrddedge[i])
    {
        printf("duplicate\n" );
        getchar();
        getchar();
    }
}
arrddedge[k%10]=reachableStates;
k++;
#endif
//1063000
// reachableStates.getForest()->underApproximate(reachableStates, 500,500);
//1000
// reachableStates.show(meddlyout,0);
 // printf("XXXX\n" );
// getchar();
end = clock();
double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
if(time_taken>timeForUA){
printf("TimeOut\n" );
return;
}
  }
  #ifdef DUP
delete [] arrddedge;
#endif
#endif
}


// ******************************************************************
// *                                                                *
// *                       forwd_bfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_mt : public common_bfs {
  public:
    forwd_bfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);
};

MEDDLY::forwd_bfs_mt::forwd_bfs_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : common_bfs(oc, a1, a2, res)
{
  if (res->getRangeType() == forest::BOOLEAN) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    setUnionOp( getOperation(MAXIMUM, res, res, res) );
  }
  setImageOp( getOperation(POST_IMAGE, a1, a2, res) );
  // FILE_output meddlyout(stdout);
  //
  // res->showInfo(meddlyout,2);
}

// ******************************************************************
// *                                                                *
// *                       forwd_bfs_ua_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_ua_mt : public common_bfs_ua {
  public:
    forwd_bfs_ua_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);
};

MEDDLY::forwd_bfs_ua_mt::forwd_bfs_ua_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : common_bfs_ua(oc, a1, a2, res)
{
  if (res->getRangeType() == forest::BOOLEAN) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    setUnionOp( getOperation(MAXIMUM, res, res, res) );
  }
  setImageOp( getOperation(POST_IMAGE, a1, a2, res) );
  // FILE_output meddlyout(stdout);
  //
  // res->showInfo(meddlyout,2);
}

// ******************************************************************
// *                                                                *
// *                       forwd_bfs_ua_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_hua_mt : public common_bfs_hua {
  public:
    forwd_bfs_hua_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);
};

MEDDLY::forwd_bfs_hua_mt::forwd_bfs_hua_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : common_bfs_hua(oc, a1, a2, res)
{
  if (res->getRangeType() == forest::BOOLEAN) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    setUnionOp( getOperation(MAXIMUM, res, res, res) );
  }
  setImageOp( getOperation(POST_IMAGE, a1, a2, res) );
  // FILE_output meddlyout(stdout);
  //
  // res->showInfo(meddlyout,2);
}

// ******************************************************************
// *                                                                *
// *                       bckwd_bfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_bfs_mt : public common_bfs {
  public:
    bckwd_bfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

};

MEDDLY::bckwd_bfs_mt::bckwd_bfs_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : common_bfs(oc, a1, a2, res)
{
  if (res->getRangeType() == forest::BOOLEAN) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    setUnionOp( getOperation(MAXIMUM, res, res, res) );
  }
  setImageOp( getOperation(PRE_IMAGE, a1, a2, res) );
}


// ******************************************************************
// *                                                                *
// *                    common_bfs_evplus  class                    *
// *                                                                *
// ******************************************************************

/*

class MEDDLY::common_bfs_evplus : public binary_operation {
  public:
  common_bfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual void compute(long aev, node_handle a, node_handle b, long& resEv, node_handle& resEvmdd) = 0;
  protected:
    binary_operation* unionMinOp;
    binary_operation* imageOp;

    inline void iterate(long ev, node_handle init, node_handle mxd, long& resEv, node_handle& resEvmdd) {
      resEv = ev;
      resEvmdd = arg1F->linkNode(init);

      node_handle prevReachable = 0;
#ifdef DEBUG_BFS
      fprintf(stderr, "Relation: %d\n", mxd);
      arg2F->showNodeGraph(stderr, mxd);
      fprintf(stderr, "Initial states: <%ld, %d>\n", ev, init);
      arg1F->showNodeGraph(stderr, init);
      long iters = 0;
#endif
#ifdef VERBOSE_BFS
      long iters = 0;
#endif
      while (prevReachable != resEvmdd) {
#ifdef VERBOSE_BFS
        iters++;
        fprintf(stimageOpderr, "Iteration %d:\n", iters);
#endif
        resF->unlinkNode(prevReachable);
        prevReachable = resEvmdd;
        long front_ev = Inf<long>();
        node_handle front = 0;
        imageOp->computeTemp(resEv, resEvmdd, mxd, front_ev, front);
#ifdef VERBOSE_BFS
        fprintf(stderr, "\timage done <%ld, %d>\n", front_ev, front);
#endif
#ifdef DEBUG_BFS
        iters++;
        fprintf(stderr, "Iteration %d\npseudo-frontier: <%ld, %d>\n", iters, front_ev, front);
        arg1F->showNodeGraph(stderr, front);
#endif
        unionMinOp->computeTemp(resEv, resEvmdd, front_ev, front, resEv, resEvmdd);
#ifdef VERBOSE_BFS
        fprintf(stderr, "\tunion done <%ld, %d>\n", resEv, resEvmdd);
#endif
#ifdef DEBUG_BFS
        fprintf(stderr, "Reachable so far: <%ld, %d>\n", resEv, resEvmdd);
        arg1F->showNodeGraph(stderr, resEvmdd);
#endif
        resF->unlinkNode(front);
      }
      resF->unlinkNode(prevReachable);
    }
};


MEDDLY::common_bfs_evplus::common_bfs_evplus(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc, 0, a1, a2, res)
{
  unionMinOp = 0;
  imageOp = 0;
}


void MEDDLY::common_bfs_evplus::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long cev = Inf<long>();
  node_handle cnode = 0;
  compute(aev, a.getNode(), b.getNode(), cev, cnode);
  c.set(cnode, cev);
}

*/

// ******************************************************************
// *                                                                *
// *                     forwd_bfs_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_evplus : public common_bfs {
  public:
  forwd_bfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

//     virtual void compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
};

MEDDLY::forwd_bfs_evplus::forwd_bfs_evplus(const binary_opname* oc, expert_forest* a1,
  // expert_forest* a2, expert_forest* res) : common_bfs_evplus(oc, a1, a2, res)
  expert_forest* a2, expert_forest* res) : common_bfs(oc, a1, a2, res)
{
  if (res->getRangeType() == forest::INTEGER) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    throw error(error::INVALID_OPERATION);
  }
  setImageOp( getOperation(POST_IMAGE, a1, a2, res) );
}

/*
void MEDDLY::forwd_bfs_evplus::compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  if (resF->getRangeType() == forest::INTEGER) {
    unionMinOp = getOperation(UNION, resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION);
  }
  imageOp = getOperation(POST_IMAGE, arg1F, arg2F, resF);

  iterate(ev, evmdd, mxd, resEv, resEvmdd);
}
*/


// ******************************************************************
// *                                                                *
// *                     forwd_bfs_ua_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_ua_evplus : public common_bfs_ua {
  public:
  forwd_bfs_ua_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

//     virtual void compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
};

MEDDLY::forwd_bfs_ua_evplus::forwd_bfs_ua_evplus(const binary_opname* oc, expert_forest* a1,
  // expert_forest* a2, expert_forest* res) : common_bfs_evplus(oc, a1, a2, res)
  expert_forest* a2, expert_forest* res) : common_bfs_ua(oc, a1, a2, res)
{
  if (res->getRangeType() == forest::INTEGER) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    throw error(error::INVALID_OPERATION);
  }
  setImageOp( getOperation(POST_IMAGE, a1, a2, res) );
}


// ******************************************************************
// *                                                                *
// *                     forwd_bfs_ua_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_hua_evplus : public common_bfs_hua {
  public:
  forwd_bfs_hua_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

//     virtual void compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
};

MEDDLY::forwd_bfs_hua_evplus::forwd_bfs_hua_evplus(const binary_opname* oc, expert_forest* a1,
  // expert_forest* a2, expert_forest* res) : common_bfs_evplus(oc, a1, a2, res)
  expert_forest* a2, expert_forest* res) : common_bfs_hua(oc, a1, a2, res)
{
  if (res->getRangeType() == forest::INTEGER) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    throw error(error::INVALID_OPERATION);
  }
  setImageOp( getOperation(POST_IMAGE, a1, a2, res) );
}
// ******************************************************************
// *                                                                *
// *                     bckwd_bfs_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_bfs_evplus : public common_bfs {
  public:
    bckwd_bfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    // virtual void compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
};

MEDDLY::bckwd_bfs_evplus::bckwd_bfs_evplus(const binary_opname* oc, expert_forest* a1,
  // expert_forest* a2, expert_forest* res) : common_bfs_evplus(oc, a1, a2, res)
  expert_forest* a2, expert_forest* res) : common_bfs(oc, a1, a2, res)
{
  if (res->getRangeType() == forest::INTEGER) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    throw error(error::INVALID_OPERATION);
  }
  setImageOp( getOperation(PRE_IMAGE, a1, a2, res) );
}

/*
void MEDDLY::bckwd_bfs_evplus::compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  if (resF->getRangeType() == forest::INTEGER) {
    unionMinOp = getOperation(UNION, resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION);
  }
  imageOp = getOperation(PRE_IMAGE, arg1F, arg2F, resF);

  iterate(ev, evmdd, mxd, resEv, resEvmdd);
}
*/


// ******************************************************************
// *                                                                *
// *                     forwd_bfs_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_opname : public binary_opname {
  public:
    forwd_bfs_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::forwd_bfs_opname::forwd_bfs_opname()
 : binary_opname("ReachableBFS")
{
}

MEDDLY::binary_operation* 
MEDDLY::forwd_bfs_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (a1->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    return new forwd_bfs_mt(this, a1, a2, r);
  }
  else if (a1->getEdgeLabeling() == forest::EVPLUS) {
    return new forwd_bfs_evplus(this, a1, a2, r);
  }
  else {
    throw error(error::TYPE_MISMATCH);
  }
}

// ******************************************************************
// *                                                                *
// *                     forwd_bfs_ua_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_ua_opname : public binary_opname {
  public:
    forwd_bfs_ua_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::forwd_bfs_ua_opname::forwd_bfs_ua_opname()
 : binary_opname("ReachableBFSUA")
{
}

MEDDLY::binary_operation*
MEDDLY::forwd_bfs_ua_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (a1->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    return new forwd_bfs_ua_mt(this, a1, a2, r);
  }
  else if (a1->getEdgeLabeling() == forest::EVPLUS) {
    return new forwd_bfs_ua_evplus(this, a1, a2, r);
  }
  else {
    throw error(error::TYPE_MISMATCH);
  }
}


// ******************************************************************
// *                                                                *
// *                     forwd_bfs_hua_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_hua_opname : public binary_opname {
  public:
    forwd_bfs_hua_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::forwd_bfs_hua_opname::forwd_bfs_hua_opname()
 : binary_opname("ReachableBFSHUA")
{
}

MEDDLY::binary_operation*
MEDDLY::forwd_bfs_hua_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (a1->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    return new forwd_bfs_hua_mt(this, a1, a2, r);
  }
  else if (a1->getEdgeLabeling() == forest::EVPLUS) {
    return new forwd_bfs_hua_evplus(this, a1, a2, r);
  }
  else {
    throw error(error::TYPE_MISMATCH);
  }
}


// ******************************************************************
// *                                                                *
// *                     bckwd_bfs_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_bfs_opname : public binary_opname {
  public:
    bckwd_bfs_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::bckwd_bfs_opname::bckwd_bfs_opname()
 : binary_opname("ReverseReachableBFS")
{
}

MEDDLY::binary_operation* 
MEDDLY::bckwd_bfs_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (a1 != r)
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) 
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (a1->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    return new bckwd_bfs_mt(this, a1, a2, r);
  }
  else if (a1->getEdgeLabeling() == forest::EVPLUS) {
    return new bckwd_bfs_evplus(this, a1, a2, r);
  }
  else {
    throw error(error::TYPE_MISMATCH);
  }
}


// ******************************************************************
// *                                                                *
// *                        common_bfs_event class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::common_bfs_event : public binary_operation_event {
  public:
    common_bfs_event(const binary_opname_event* opcode, expert_forest* arg1,
      dd_edge* arg2, int arg3, expert_forest* res);

    virtual void computeDDEdgeEvent(const dd_edge& a, const dd_edge* b, int c,dd_edge &d, bool userFlag);

  protected:
    inline void setUnionOp(binary_operation* uop)
    {
      MEDDLY_DCASSERT(uop);
      MEDDLY_DCASSERT(0==unionOp);
      unionOp = uop;
    }

    inline void setImageOp(binary_operation* iop)
    {
      MEDDLY_DCASSERT(iop);
      MEDDLY_DCASSERT(0==imageOp);
      imageOp = iop;
    }

  private:
    binary_operation* unionOp;
    binary_operation* imageOp;

};


MEDDLY::common_bfs_event::common_bfs_event(const binary_opname_event* oc, expert_forest* a1,
  dd_edge* a2,int a3, expert_forest* res)
: binary_operation_event(oc, 0, a1, a2,a3, res)
{
  unionOp = 0;
  imageOp = 0;
}

void MEDDLY::common_bfs_event::computeDDEdgeEvent(const dd_edge &init, const dd_edge *R,int n, dd_edge &reachableStates, bool userFlag)
{
    printf("Need to be write down!\n" );
    clock_t start, end;
    start = clock();
  MEDDLY_DCASSERT(unionOp);
  MEDDLY_DCASSERT(imageOp);

  reachableStates = init;
  dd_edge prevReachable(resF);
  dd_edge front(resF);
  FILE_output meddlyout(stdout);
  long peakreachable=0;
  long lastreachable=0;
#ifdef DEBUG_BFS
  FILE_output debug(stderr);
  debug << "Relation: ";
  R.show(debug, 2);
  debug << "Initial states: ";
  init.show(debug, 2);
  long iters = 0;
#endif
#ifdef VERBOSE_BFS
  long iters = 0;
  FILE_OUTPUT verbose(stderr);
#endif
  while (prevReachable != reachableStates) {
      // printf("XXXXXIt Done0000\n" );
      // R.show(meddlyout,0);
#ifdef VERBOSE_BFS
    iters++;
    verbose << "Iteration " << iters << ":\n";
#endif
    // printf("XXXXXIt Done0\n" );
    prevReachable = reachableStates;
    // printf("XXXXXIt Done 1\n" );

     // reachableStates.show(meddlyout, 0);
     for(int e=0; e<n; e++){
    imageOp->computeDDEdge(reachableStates, R[e], front, userFlag);
    // printf("XXXXXIt Done2\n" );

#ifdef VERBOSE_BFS
    verbose << "\timage done ";
    front.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    iters++;
    debug << "Iteration " << iters << "\npseudo-frontier: ";
    front.show(debug, 2);
#endif
    unionOp->computeDDEdge(reachableStates, front, reachableStates, userFlag);
    }
#ifdef VERBOSE_BFS
    verbose << "\tunion done ";
    reachableStates.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    debug << "Reachable so far: ";
    reachableStates.show(debug, 2);
#endif
printf("cardafter %f\n",reachableStates.getCardinality() );
lastreachable=reachableStates.getNodeCount();
#ifdef CHKBFS
if(lastreachable>maxThreshold)
return;
#endif
if(lastreachable>peakreachable)
{peakreachable=lastreachable;}
printf("XXXX %ld\t %ld\n",lastreachable, peakreachable );

 // reachableStates.getForest()->underApproximate(reachableStates,1100);
//1000
// reachableStates.show(meddlyout,0);
 // printf("XXXX\n" );
 end = clock();
 double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
 if(time_taken>timeT){
 printf("TimeOut\n" );
 return;
}
  }




}


// ******************************************************************
// *                                                                *
// *                       forwd_bfs_mt_event class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_mt_event : public common_bfs_event {
  public:
    forwd_bfs_mt_event(const binary_opname_event* opcode, expert_forest* arg1,
      dd_edge* arg2, int arg3, expert_forest* res);
};

MEDDLY::forwd_bfs_mt_event::forwd_bfs_mt_event(const binary_opname_event* oc, expert_forest* a1,
 dd_edge* a2, int a3, expert_forest* res) : common_bfs_event(oc, a1, a2,a3, res)
{
    printf("start setting image and union operation\n" );

  if (res->getRangeType() == forest::BOOLEAN) {
    setUnionOp( getOperation(UNION, res, res, res) );
    printf("Done setting union operation\n" );
  } else {
    setUnionOp( getOperation(MAXIMUM, res, res, res) );
  }
  printf("start setting image operation\n" );

  setImageOp( getOperation(POST_IMAGE, a1, a2[0], res) );
  printf("Done setting image and union operation\n" );
  // FILE_output meddlyout(stdout);
  //
  // res->showInfo(meddlyout,2);
}



// ******************************************************************
// *                                                                *
// *                     forwd_bfs_opname_event class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_opname_event : public binary_opname_event {
  public:
    forwd_bfs_opname_event();
    virtual binary_operation_event* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const{return 0;};
      virtual binary_operation_event* buildOperation(expert_forest* a1,
        dd_edge* a2, int a3, expert_forest* r)const;
};

MEDDLY::forwd_bfs_opname_event::forwd_bfs_opname_event()
 : binary_opname_event("ReachableBFSEvent")
{
}

MEDDLY::binary_operation_event*
MEDDLY::forwd_bfs_opname_event::buildOperation(expert_forest* a1, dd_edge* a2,int a3,
  expert_forest* r) const
{
    printf("XXXX came here buildOperation\n" );
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain())// ||
    // (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    // !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) //||
    // (a2->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    {
        printf("throw mismatch\n" );
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

  if (a1->getEdgeLabeling() == forest::MULTI_TERMINAL) {
      printf("Should call something\n");
      // return 0;
    return new forwd_bfs_mt_event(this, a1, a2,a3, r);
  }
  else if (a1->getEdgeLabeling() == forest::EVPLUS) {
      printf("Should call something2\n");

    // return new forwd_bfs_hua_evplus(this, a1, a2, r);
  }
  else {
    throw error(error::TYPE_MISMATCH);
  }
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeForwardBFS()
{
  return new forwd_bfs_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeBackwardBFS()
{
  return new bckwd_bfs_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeForwardBFSUA()
{
  return new forwd_bfs_ua_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeForwardBFSHUA()
{
  return new forwd_bfs_hua_opname;
}

MEDDLY::binary_opname_event* MEDDLY::initializeAllBFSGen()
{
    printf("Need to get completed!!! MEDDLY::initializeAllBFSGen\n" );
    return new forwd_bfs_opname_event;
}

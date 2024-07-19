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
 #include <climits>

int TimeLimit=120000;//INT_MAX;
#include "../defines.h"
#include "cov.h"
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
namespace MEDDLY {
class cov_by_events_opname;
class cov_by_events_op;

class common_cov_by_events_mt;
class cov_by_events_mt;

class cov_opname;
class compareij;
};

// ******************************************************************
// *                                                                *
// *                    cov_by_events_opname  class              *
// *                                                                *
// ******************************************************************

/** Simple class to keep compute table happy.
 */
class MEDDLY::cov_by_events_opname : public unary_opname {
static cov_by_events_opname* instance;
public:
cov_by_events_opname();

static cov_by_events_opname* getInstance();

};

MEDDLY::cov_by_events_opname* MEDDLY::cov_by_events_opname::instance = 0;

MEDDLY::cov_by_events_opname::cov_by_events_opname()
        : unary_opname("Cov_by_events")
{
}

MEDDLY::cov_by_events_opname* MEDDLY::cov_by_events_opname::getInstance()
{
        if (0==instance) instance = new cov_by_events_opname;
        return instance;
}

// ******************************************************************
// *                                                                *
// *                      cov_by_events_op  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::cov_by_events_op : public unary_operation {
common_cov_by_events_mt* parent;
public:
cov_by_events_op(common_cov_by_events_mt* p,
                 expert_forest* argF, expert_forest* resF);
virtual ~cov_by_events_op();

void RSBFSfly(const dd_edge& in, dd_edge& out, satotf_opname::otf_relation* rel,  markcmp* cij);
void Coverability(const dd_edge& in, dd_edge& out, satotf_opname::otf_relation* rel,  markcmp* cij);

protected:
inline ct_entry_key*
findSaturateResult(node_handle a, int level, node_handle& b) {
        ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->writeN(a);
        if (argF->isFullyReduced()) CTsrch->writeI(level);
        CT0->find(CTsrch, CTresult[0]);
        if (!CTresult[0]) return CTsrch;
        b = resF->linkNode(CTresult[0].readN());
        CT0->recycle(CTsrch);
        return 0;
}
inline node_handle saveSaturateResult(ct_entry_key* Key,
                                      node_handle a, node_handle b)
{
        CTresult[0].reset();
        CTresult[0].writeN(b);
        CT0->addEntry(Key, CTresult[0]);
        return b;
}
};

// ******************************************************************
// *                                                                *
// *            common_cov_by_events_mt  class                  *
// *                                                                *
// ******************************************************************

class MEDDLY::common_cov_by_events_mt : public specialized_operation {
public:
common_cov_by_events_mt(satotf_opname* opcode,
                        satotf_opname::otf_relation* rel,markcmp* cij);
virtual ~common_cov_by_events_mt();

virtual void compute(const dd_edge& a, dd_edge &c /*,compareij *cij*/);
virtual void computeMXD(const dd_edge& a, dd_edge &c);

virtual void setComparefunction(compareij* cij);

protected:
inline ct_entry_key*
findResult(node_handle a, node_handle b, node_handle &c)
{
        ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->writeN(a);
        CTsrch->writeN(b);
        CT0->find(CTsrch, CTresult[0]);
        if (!CTresult[0]) return CTsrch;
        c = resF->linkNode(CTresult[0].readN());
        CT0->recycle(CTsrch);
        return 0;
}
inline node_handle saveResult(ct_entry_key* Key,
                              node_handle a, node_handle b, node_handle c)
{
        CTresult[0].reset();
        CTresult[0].writeN(c);
        CT0->addEntry(Key, CTresult[0]);
        return c;
}

public:
binary_operation* mddUnion;
binary_operation* mddImage;
binary_operation* mrcc;
binary_operation* covtc;
unary_operation* covrExtractFrom;
unary_operation* covrExtractCover;
binary_operation* covrCoveredFrom;
binary_operation* covrCoveredTo;
markcmp* cij;
protected:
binary_operation* mxdIntersection;
binary_operation* mxdDifference;

satotf_opname::otf_relation* rel;

expert_forest* arg1F;
expert_forest* arg2F;
expert_forest* resF;

protected:
class indexq {
static const int NULPTR = -1;
static const int NOTINQ = -2;
int* data;
unsigned size;
int head;
int tail;
public:
// used by parent for recycling
indexq* next;
public:
indexq();
~indexq();
void resize(unsigned sz);
inline bool isEmpty() const {
        return NULPTR == head;
}
inline void add(int i) {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, (unsigned)i, size);
        if (NOTINQ != data[i]) return;
        if (NULPTR == head) {
                // empty list
                head = i;
        } else {
                // not empty list
                MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, (unsigned)tail, size);
                data[tail] = i;
        }
        tail = i;
        data[i] = NULPTR;
}
inline int remove() {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, (unsigned)head, size);
        int ans = head;
        head = data[head];
        data[ans] = NOTINQ;
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, (unsigned)ans, size);
        return ans;
}
};

protected:
class charbuf {
public:
char* data;
unsigned size;
charbuf* next;
public:
charbuf();
~charbuf();
void resize(unsigned sz);
};

private:
indexq* freeqs;
charbuf* freebufs;

protected:
inline indexq* useIndexQueue(unsigned sz) {
        indexq* ans;
        if (freeqs) {
                ans = freeqs;
                freeqs = freeqs->next;
        } else {
                ans = new indexq();
        }
        MEDDLY_DCASSERT(ans);
        ans->resize(sz);
        ans->next = 0;
        return ans;
}
inline void recycle(indexq* a) {
        MEDDLY_DCASSERT(a);
        MEDDLY_DCASSERT(a->isEmpty());
        a->next = freeqs;
        freeqs = a;
}

inline charbuf* useCharBuf(unsigned sz) {
        charbuf* ans;
        if (freebufs) {
                ans = freebufs;
                freebufs = freebufs->next;
        } else {
                ans = new charbuf();
        }
        MEDDLY_DCASSERT(ans);
        ans->resize(sz);
        ans->next = 0;
        return ans;
}
inline void recycle(charbuf* a) {
        MEDDLY_DCASSERT(a);
        a->next = freebufs;
        freebufs = a;
}

inline virtual bool checkForestCompatibility() const {
        return true;
}
};

// ******************************************************************
// *                                                                *
// *           cov_by_events_op  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::cov_by_events_op
::cov_by_events_op(common_cov_by_events_mt* p,
                   expert_forest* argF, expert_forest* resF)
        : unary_operation(cov_by_events_opname::getInstance(), 1, argF, resF)
{
        parent = p;

        const char* name = cov_by_events_opname::getInstance()->getName();
        ct_entry_type* et;

        if (argF->isFullyReduced()) {
                // CT entry includes level info
                et = new ct_entry_type(name, "NI:N");
                et->setForestForSlot(0, argF);
                et->setForestForSlot(3, resF);
        } else {
                et = new ct_entry_type(name, "N:N");
                et->setForestForSlot(0, argF);
                et->setForestForSlot(2, resF);
        }
        registerEntryType(0, et);
        buildCTs();
}

MEDDLY::cov_by_events_op::~cov_by_events_op()
{
        removeAllComputeTableEntries();
}


void MEDDLY::cov_by_events_op::RSBFSfly(const dd_edge& init, dd_edge& reachableStates, satotf_opname::otf_relation* rel,markcmp* cij)
{
        reachableStates=init;
        dd_edge prevReachable(resF);
        dd_edge front(resF);
        while (prevReachable != reachableStates) {
                prevReachable = reachableStates;
                for(int level=1; level<argF->getNumVariables()+1; level++) {
                        for (int ei = 0; ei < rel->getNumOfEvents(level); ei++) {
                                rel->rebuildEvent(level, ei);
                        }
                }
                for(int level=1; level<argF->getNumVariables()+1; level++) {
                        int nEventsAtThisLevel = rel->getNumOfEvents(level);
                        for (int ei = 0; ei < nEventsAtThisLevel; ei++) {
                                rel->rebuildEvent(level, ei);
                                const dd_edge& mxd = rel->getEvent(level, ei);
                                if (0==mxd.getNode()) {}
                                else{
                                        std::list<int>* shouldConfirm=new std::list<int>[argF->getNumVariables()+1];
                                        parent->mddImage->computeDDEdgeSC(reachableStates, mxd, front, true,shouldConfirm);
                                        for(int listidx=0; listidx<argF->getNumVariables()+1; listidx++) {
                                                for (auto const &scidx: shouldConfirm[listidx]) {
                                                        rel->confirm(listidx, int(scidx));
                                                        cij->compare(int(scidx),int(scidx),listidx);
                                                        printf("compare %d\n",listidx );

                                                        getchar();
                                                }
                                        }

                                        parent->mddUnion->computeDDEdge(reachableStates, front, reachableStates, true);
                                        delete[]shouldConfirm;
                                }
                        }
                }
        }
}


#include <chrono>
using namespace std::chrono;

void MEDDLY::cov_by_events_op::Coverability(const dd_edge& init, dd_edge& reachableStates, satotf_opname::otf_relation* rel,markcmp* cij)
{
        reachableStates=init;
        dd_edge prevReachable(resF);
        dd_edge front(resF);
        int i=0;
        bool first=false;
        ostream_output meddlyout(std::cout);
        bool debug=false;

        // bool test=true;
        if(debug){
        init.showGraph(meddlyout);
       printf("init^^^");
        }
        // while (prevReachable != reachableStates) {
        auto start = high_resolution_clock::now();
        do{

                first=false;
                i++;
                prevReachable = reachableStates;
                for(int level=1; level<argF->getNumVariables()+1; level++) {
                        for (int ei = 0; ei < rel->getNumOfEvents(level); ei++) {
                                rel->rebuildEvent(level, ei);
                                // const dd_edge& mxd = rel->getEvent(level, ei);
                                // mxd.showGraph(meddlyout);
                                // printf("MXD %d ei!\n",ei );
                                // getchar();
                        }
                }
                // dd_edge urel(argF);
                // urel=rel->getMonolithicNSF();
                std::list<int>* shouldConfirms=new std::list<int>[argF->getNumVariables()+1];
                dd_edge nfront(resF);
                dd_edge efront(resF);
                dd_edge lfront(resF);
                dd_edge tlfront(resF);
                dd_edge updatecovered(resF);
                int gfront=0;
                nfront=init;
                efront=init;
                lfront=init;
                tlfront=init;
                updatecovered=reachableStates;

                // reachableStates.showGraph(meddlyout);
                // printf("reachableStates^^BEFORE!\n");
                // getchar();
                // urel.showGraph(meddlyout);
                // printf("transition relation^^\n" );
                // getchar();
                for(int level=1; level<argF->getNumVariables()+1; level++,gfront=0) {
                        for (int ei = 0; ei < rel->getNumOfEvents(level); ei++) {
                                const dd_edge& urel = rel->getEvent(level, ei);
                                printf("transition relation LVL %d %d^^\n",level,ei );

                                if(debug){
                                 urel.showGraph(meddlyout);
                                printf("transition relation LVL %d %d^^\n",level,ei );
                                 getchar();
                                 reachableStates.showGraph(meddlyout);
                                printf("reachableStates^^" );
                                getchar();
                                }
                                // printf("Before calling mrrc\n" );
                                parent->mrcc->computeDDEdgeSC(/*reachableStates*/prevReachable, urel, nfront,efront, lfront,tlfront, gfront, true,shouldConfirms,parent->cij);
                                // prevReachable=reachableStates;
                                // printf("mrrc computed^^^ " );
                                long cardetest=0;
                                // apply(CARDINALITY, lfront, cardetest);
                                // if(cardetest==0){
                                //     printf("%d\n",lfront.getNode() );
                                //     getchar();
                                // }
                                // printf("%d\n",lfront.getNode() );
                                // printf("CARD lfront %ld\n",cardetest );
                                if(gfront) first=true;
                                if(lfront.getNode()>0 /*gfront||cardetest>0*/)
                                {

                                 printf("OMEGA FOUND\n" );

                                resF->reportStats(meddlyout, "\t",
                                  expert_forest::HUMAN_READABLE_MEMORY |
                                  expert_forest::BASIC_STATS
                                );
                                auto stop = high_resolution_clock::now();
                                auto duration = duration_cast<microseconds>(stop - start);
                                printf("duration %ld microseconds \n", duration.count());
                                // getchar();
                                    }
                                        // printf("gfront %d\n",gfront );
                                // dd_edge etest(resF);
                                // apply(EQUANT,ffront, etest);
                                if(debug){
                                apply(CARDINALITY, nfront, cardetest);
                                printf("CARD nfront %ld\n",cardetest );
                                apply(CARDINALITY, lfront, cardetest);
                                printf("CARD lfront %ld\n",cardetest );
                                apply(CARDINALITY, tlfront, cardetest);
                                printf("CARD tlfront %ld\n",cardetest );
                                apply(CARDINALITY, efront, cardetest);
                                printf("CARD efront %ld\n",cardetest );
                                }
                                // getchar();
                                // etest.showGraph(meddlyout);
                                // printf("ETEST^^\n" );
                                // efront.showGraph(meddlyout);
                                // printf("EFRONT^^\n" );
                                // getchar();
                                parent->mddUnion->computeDDEdge(reachableStates, nfront, reachableStates, true);

                                parent->mddUnion->computeDDEdge(reachableStates, efront, reachableStates, true);

                                // parent->mddUnion->computeDDEdgeOmega(reachableStates, ffront, reachableStates, true, parent->cij);

                                // reachableStates.showGraph(meddlyout);
                                // printf("reachableStates^^\n");
                                // getchar();

                                // parent->mddUnion->computeDDEdge(reachableStates, lfront, reachableStates, true);
                                // reachableStates.showGraph(meddlyout);
                                // printf("reachableStates^^\n");
                                // long cardetest=0;

                                if(debug){
                                apply(CARDINALITY, reachableStates, cardetest);
                                printf("CARD reachset %ld\n",cardetest );
                                }
                                // getchar();
                                for(int listidx=0; listidx<argF->getNumVariables()+1; listidx++) {
                                        for (auto const &scidx: shouldConfirms[listidx]) {
                                                rel->confirm(listidx, int(scidx));
                                                // printf("listidx %d, scidx %d\n",listidx, int(scidx) );

                                        }
                                }
                                // delete[]shouldConfirms;
                                // if(test){
                                dd_edge luniontl(resF);
                                luniontl=lfront;
                                parent->mddUnion->computeDDEdge(lfront, tlfront, luniontl, true);

                                if(debug){
                                lfront.showGraph(meddlyout);
                                printf("l front^^^" );
                                tlfront.showGraph(meddlyout);
                                printf("tl front ^^^" );
                                luniontl.showGraph(meddlyout);
                               printf("l union tl ^^^" );
                               getchar();
                                }

                                parent->covtc->computeDDEdge(reachableStates/*prevReachable*/, luniontl, updatecovered,true);
                                if(debug){
                                updatecovered.showGraph(meddlyout);
                                printf("updatecovered^^^ " );
                                }
                               ///newly added to union the result with reachablestates
                               parent->mddUnion->computeDDEdge(reachableStates, updatecovered, reachableStates, true);
                                // }
                               // getchar();
                               dd_edge efrom(resF);

                               parent->covrExtractFrom->computeDDEdge(reachableStates,efrom,true);
                               // long pcs=0;
                               //  apply(CARDINALITY, efrom, pcs);
                               // printf("PCS %ld\n",pcs );
                               // printf("number of nodes %ld\n",efrom.getNodeCount() );
                               // printf("number of edge %ld\n",efrom.getEdgeCount() );

                               // if(test){
                               if(debug){
                               printf("covrExtractFrom^^^\n" );
                               efrom.showGraph(meddlyout);
                                }
                               // dd_edge eefrom(argF);
                               // eefrom=efrom;
                               if(debug){
                               efrom.showGraph(meddlyout);
                               printf("eefrom^^^\n" );
                                }

                               dd_edge ecovered(resF);
                               parent->covrExtractCover->computeDDEdge(efrom,ecovered,parent->cij,true);
                               // getchar();
                               efrom-=ecovered;
                               long pcs=0;
                                apply(CARDINALITY, efrom, pcs);
                               printf("PCS %ld\n",pcs );
                               printf("number of nodes %ld\n",efrom.getNodeCount() );
                               printf("number of edge %ld\n",efrom.getEdgeCount() );

                               // // printf("mddsize %d\n",efrom.getNodeCount()*sizeof(efrom) );
                               // printf("mddsize %ld\n",pcs*argF->getNumVariables()* sizeof(int));

                               if(debug){
                               efrom.showGraph(meddlyout);
                               printf("Efrom minus ecovered^^\n" );
                                }

                               // dd_edge coveredFrom(resF);
                               // coveredFrom=init;
                               // parent->covrCoveredFrom->computeDDEdge(reachableStates,efrom,coveredFrom,true);
                               // parent->covrCoveredTo->computeDDEdge(reachableStates,efrom,coveredFrom,true);

                               parent->covrCoveredFrom->computeDDEdge(reachableStates,efrom,reachableStates,true);
                               parent->covrCoveredTo->computeDDEdge(reachableStates,efrom,reachableStates,true);
                                // }
                               // getchar();
                               // if(gfront){
                               //     // resF->showInfo(meddlyout,2);
                               //     resF->reportStats(meddlyout, "\t",
                               //       expert_forest::HUMAN_READABLE_MEMORY |
                               //       expert_forest::BASIC_STATS
                               //     );
                               //     auto stop = high_resolution_clock::now();
                               //     auto duration = duration_cast<microseconds>(stop - start);
                               //     printf("duration %ld microseconds \n", duration.count());
                               //
                               // }
                               cardetest=0;
                               apply(CARDINALITY, reachableStates, cardetest);
                               printf("CARD CR %ld\n",cardetest );
                               printf("number of nodes CR %ld\n",reachableStates.getNodeCount() );
                               printf("number of edge CR %ld\n",reachableStates.getEdgeCount() );
                        }
                }
                // reachableStates.showGraph(meddlyout);
                // printf("reachableStates^^\n");
                // long cardier;
                // apply(CARDINALITY, reachableStates, cardier);
                // printf("i %d CARD CR %ld\n",i,cardier );
                // getchar();

                // reachableStates.showGraph(meddlyout);
                // getchar();
                delete[]shouldConfirms;
                // printf("PR %d, F %d\n",prevReachable == reachableStates, first==true );
                // if(prevReachable == reachableStates&& first==true) {
                //         first=false;
                // }else if(prevReachable == reachableStates&& first==false) {
                //         first=true;
                // }
                auto stop = high_resolution_clock::now();
                auto duration = duration_cast<minutes>(stop - start);
                if(duration.count()>TimeLimit) {
                        printf("TimeOut\n" );
                        return;
                }
                // printf("prevReachable %d reachableStates\n",prevReachable == reachableStates );
        }while (prevReachable != reachableStates || first==true);
        printf("DONEDONE %d\n",i );

}

// ******************************************************************
// *                                                                *
// *           common_cov_by_events_mt  methods                     *
// *                                                                *
// ******************************************************************


MEDDLY::common_cov_by_events_mt::common_cov_by_events_mt(
        satotf_opname* opcode,
        satotf_opname::otf_relation* relation,markcmp* _cij)
        : specialized_operation(opcode, 1)
{
        mddUnion = 0;
        mddImage = 0;
        mrcc=0;
        covtc=0;
        covrExtractFrom=0;
        covrExtractCover=0;
        covrCoveredFrom=0;
        covrCoveredTo=0;
        mxdIntersection = 0;
        mxdDifference = 0;
        freeqs = 0;
        freebufs = 0;
        rel = relation;
        cij=_cij;
        arg1F = static_cast<expert_forest*>(rel->getInForest());
        arg2F = static_cast<expert_forest*>(rel->getRelForest());
        resF = static_cast<expert_forest*>(rel->getOutForest());
        registerInForest(arg1F);
        registerInForest(arg2F);
        registerInForest(resF);

        ct_entry_type* et = new ct_entry_type(opcode->getName(), "NN:N");
        et->setForestForSlot(0, arg1F);
        et->setForestForSlot(1, arg2F);
        et->setForestForSlot(3, resF);
        registerEntryType(0, et);
        buildCTs();
}

MEDDLY::common_cov_by_events_mt::~common_cov_by_events_mt()
{
        if (rel->autoDestroy()) delete rel;
        unregisterInForest(arg1F);
        unregisterInForest(arg2F);
        unregisterInForest(resF);
}

void MEDDLY::common_cov_by_events_mt
::setComparefunction(compareij* cij){

}

void MEDDLY::common_cov_by_events_mt
::computeMXD(const dd_edge &a, dd_edge &c)
{
        // Initialize operations
        mddUnion = getOperation(UNION, arg2F, arg2F, arg2F);
        MEDDLY_DCASSERT(mddUnion);

        mxdIntersection = getOperation(INTERSECTION, arg2F, arg2F, arg2F);
        MEDDLY_DCASSERT(mxdIntersection);

        mxdDifference = getOperation(DIFFERENCE, arg2F, arg2F, arg2F);
        MEDDLY_DCASSERT(mxdDifference);

        mddImage=getOperation(POST_IMAGE, arg1F, arg2F, resF);
        MEDDLY_DCASSERT(mddImage);

        dd_edge d(resF);
        dd_edge e(resF);
        dd_edge f(resF);
        mrcc=getOperationW(MRC_POST_IMAGE,a,arg2F,c,d,e,f);
        MEDDLY_DCASSERT(mrcc);
        covtc=getOperation(COV_TC, a,arg2F,c);
        MEDDLY_DCASSERT(covtc);
        dd_edge efrom(resF);
        covrExtractFrom=getOperation(EXTRACTFROM,a,efrom);
        MEDDLY_DCASSERT(covrExtractFrom);
        dd_edge afrom(arg1F);
        covrExtractCover=getOperation(EXTRACTCOVERED,afrom,efrom);
        MEDDLY_DCASSERT(covrExtractCover);
        // dd_edge coveredFrom(resF);
        covrCoveredFrom=getOperation(COVEREDFROM,a,afrom, c);
        MEDDLY_DCASSERT(covrCoveredFrom);
        // dd_edge coveredTo(resF);
        covrCoveredTo=getOperation(COVEREDTO,a,afrom, c);
        MEDDLY_DCASSERT(covrCoveredTo);
        ostream_output meddlyout(std::cout);
#ifdef DEBUG_INITIAL
        printf("Calling saturate for states:\n");
        a.showGraph(stdout);
#endif
#ifdef DEBUG_NSF
        printf("Calling saturate for NSF:\n");
        // b.showGraph(stdout);
#endif

        // Execute saturation operation
        cov_by_events_op* so = new cov_by_events_op(this, arg1F, resF);
        so->Coverability(a,c,rel,cij);

        // Cleanup
        while (freeqs) {
                indexq* t = freeqs;
                freeqs = t->next;
                delete t;
        }
        while (freebufs) {
                charbuf* t = freebufs;
                freebufs = t->next;
                delete t;
        }
        delete so;
}


void MEDDLY::common_cov_by_events_mt
::compute(const dd_edge &a, dd_edge &c)
{
        printf("MEDDLY::common_cov_by_events_mt\n" );
        // Initialize operations
        mddUnion = getOperation(UNION, resF, resF, resF);
        MEDDLY_DCASSERT(mddUnion);

        mxdIntersection = getOperation(INTERSECTION, arg2F, arg2F, arg2F);
        MEDDLY_DCASSERT(mxdIntersection);

        mxdDifference = getOperation(DIFFERENCE, arg2F, arg2F, arg2F);
        MEDDLY_DCASSERT(mxdDifference);

        mddImage=getOperation(POST_IMAGE, arg1F, arg2F, resF);
        MEDDLY_DCASSERT(mddImage);
        ostream_output meddlyout(std::cout);
#ifdef DEBUG_INITIAL
        printf("Calling saturate for states:\n");
        a.showGraph(stdout);
#endif
#ifdef DEBUG_NSF
        printf("Calling saturate for NSF:\n");
        // b.showGraph(stdout);
#endif

        // Execute saturation operation
        cov_by_events_op* so = new cov_by_events_op(this, arg1F, resF);
        so->RSBFSfly(a,c,rel,cij);

        // Cleanup
        while (freeqs) {
                indexq* t = freeqs;
                freeqs = t->next;
                delete t;
        }
        while (freebufs) {
                charbuf* t = freebufs;
                freebufs = t->next;
                delete t;
        }
        delete so;
}

// ******************************************************************
// *       common_cov_by_events_mt::indexq  methods                 *
// ******************************************************************

MEDDLY::common_cov_by_events_mt::indexq::indexq()
{
        data = 0;
        size = 0;
        head = NULPTR;
}

MEDDLY::common_cov_by_events_mt::indexq::~indexq()
{
        free(data);
}

void MEDDLY::common_cov_by_events_mt::indexq::resize(unsigned sz)
{
        if (sz <= size) return;
        data = (int*) realloc(data, sz * sizeof(int));
        if (0==data)
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);

        for (; size < sz; size++) data[size] = NOTINQ;
}

// ******************************************************************
// *       common_cov_by_events_mt::charbuf methods             *
// ******************************************************************

MEDDLY::common_cov_by_events_mt::charbuf::charbuf()
{
        data = 0;
        size = 0;
}

MEDDLY::common_cov_by_events_mt::charbuf::~charbuf()
{
        free(data);
}

void MEDDLY::common_cov_by_events_mt::charbuf::resize(unsigned sz)
{
        if (sz <= size) return;
        data = (char*) realloc(data, sz * sizeof(char));
        if (0==data)
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *             cov_by_events_mt class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::cov_by_events_mt : public common_cov_by_events_mt {
public:
cov_by_events_mt(satotf_opname* opcode,
                 satotf_opname::otf_relation* rel, markcmp* cij);
protected:
};

MEDDLY::cov_by_events_mt::cov_by_events_mt(
        satotf_opname* opcode,
        satotf_opname::otf_relation* rel,markcmp* cij)
        : common_cov_by_events_mt(opcode, rel,cij)
{

}

// ******************************************************************
// *                                                                *
// *                       cov_opname class                         *
// *                                                                *
// ******************************************************************

class MEDDLY::cov_opname : public satotf_opname {
bool forward;
markcmp* mcmp;
public:
cov_opname(bool fwd);
virtual specialized_operation* buildOperation(arguments* a);
virtual specialized_operation* buildOperationX(arguments* a,markcmp* a2);
};

MEDDLY::cov_opname::cov_opname(bool fwd)
        : satotf_opname(fwd ? "cov" : "cov")
{
        forward = fwd;
}

MEDDLY::specialized_operation*
MEDDLY::cov_opname::buildOperation(arguments* a)
{
        otf_relation* rel = dynamic_cast<otf_relation*>(a);
        if (0==rel) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);

        //
        // No sanity checks needed here; we did them already when constructing a.
        //

        MEDDLY::specialized_operation* op = 0;
        if (forward)
                op = new cov_by_events_mt(this, rel,0);
        else {
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }

        // Do we need to delete rel here?
        // No, if needed, do this in the destructor for op.

        return op;
}


MEDDLY::specialized_operation*
MEDDLY::cov_opname::buildOperationX(arguments* a,markcmp* a2)
{
        otf_relation* rel = dynamic_cast<otf_relation*>(a);
        if (0==rel) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
        //
        // No sanity checks needed here; we did them already when constructing a.
        //


        MEDDLY::specialized_operation* op = 0;
        if (forward)
                op = new cov_by_events_mt(this, rel,a2);
        else {
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }

        // Do we need to delete rel here?
        // No, if needed, do this in the destructor for op.

        return op;
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::satotf_opname* MEDDLY::initCov()
{
        printf("This is COV!\n" );
        return new cov_opname(true);
}

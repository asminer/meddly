
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



/*
  State space generation for Dining Philosphers (N=2).

  The model has 2 philosophers and 2 forks.

  Each philosopher can be in state {I, W, L, R, E} where

  I:    idle philosopher
  WB:   philosopher is waiting for both forks
  HL:   philosopher has left fork
  HR:   philosopher has right fork
  E:    philosopher is eating

  Each fork can be in state {A, NA} where

  A:    fork is available
  NA:   fork is not available

  Philosphers can move from one state to another as:

  I -> WB

  The synchronization between philosopher 1 and the forks:

  WB1 ->  HR1
  A1  ->  NA1

  WB1 ->  HL1
  A2  ->  NA2

  HR1 ->  E1
  A2  ->  NA2

  HL1 ->  E1
  A1  ->  NA1

  E1  ->  I1
  NA1 ->  A1
  NA2 ->  A2

  The synchronization between philosopher 2 and the forks:

  WB2 ->  HR2
  A2  ->  NA2

  WB2 ->  HL2
  A1  ->  NA1

  HR2 ->  E2
  A1  ->  NA1

  HL2 ->  E2
  A2  ->  NA2

  E2  ->  I2
  NA1 ->  A1
  NA2 ->  A2

  Initially, all forks are in state "A" and all philosophers are in state "I".
  How many reachable states?
  Exceptions?
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <fstream>

#include "../config.h"
#if HAVE_LIBGMP
#include <gmp.h>
#endif

#include "../src/meddly.h"
#include "../timing/timer.h"
#include "../src/log_simple.h"

// #define OLD_ITERATORS

using namespace MEDDLY;

// #define NAME_VARIABLES
// #define SHOW_MXD
// #define SHOW_MDD

// Specify the variable order, top down:
enum varorder {
    fpfp,       // fork, phil, fork, phil, ...
    pfpf,       // phil, fork, phil, fork, ...
    ppff,       // phil, ..., phil, fork, ..., fork
    ffpp        // fork, ..., fork, phil, ..., phil
};

struct switches {
        bool mark_sweep;
        bool pessimistic;
        bool exact;
        char method;
        // bool chaining;
        bool printReachableStates;
        varorder vord;
        bool build_pdf;

    public:
        switches() {
            mark_sweep = false;
            pessimistic = false;
            exact = false;
            method = 'b';
            // chaining = true;
            printReachableStates = false;
            vord = fpfp;
            build_pdf = false;
        }
};

// **********************************************************************
// Class to deal with different variable orders
// **********************************************************************

class model2var {
    public:
        model2var(varorder vo, int phils);

        inline int philVar(int p) const {
            switch (vord) {

                case ffpp:  return p+1;

                case ppff:  return nPhils+p+1;

                case pfpf:  return 2*p +2;  // 0th philosopher at 2nd variable

                case fpfp:
                default:    return 2*p +1;  // 0th philosopher at 1st variable
            }
        }

        inline int forkVar(int p) const {
            switch (vord) {

                case ffpp:  return nPhils+p+1;

                case ppff:  return p+1;

                case pfpf:  return 2*p +1;  // 0th fork at 1st variable

                case fpfp:
                default:    return 2*p +2;  // 0th fork at 2nd variable
            }

        }

        inline int leftVar(int p) const {
            return forkVar( (p + nPhils - 1) % nPhils );
        }
        inline int rightVar(int p) const {
            return forkVar(p);
        }

        inline int numPhils() const {
            return nPhils;
        }

        inline int numLevels() const {
            return 2*nPhils;
        }

    private:
        varorder vord;
        int nPhils;
};


model2var::model2var(varorder vo, int phils) : vord(vo)
{
    nPhils = phils;
}

// **********************************************************************
// Class for building the transition relation
// **********************************************************************

class philsModel {
    public:
        philsModel(int nPhils, const model2var& m2v, forest& mxd);
        ~philsModel();

        // event builders
        inline void setPhilosopher(int phil) {
            ph = M2V.philVar(phil);
            rf = M2V.rightVar(phil);
            lf = M2V.leftVar(phil);
        };

        void Idle2WaitBoth(dd_edge &e);
        void WaitBoth2HaveRight(dd_edge &e);
        void WaitBoth2HaveLeft(dd_edge &e);
        void HaveRight2Eat(dd_edge &e);
        void HaveLeft2Eat(dd_edge &e);
        void Eat2Idle(dd_edge &e);

        // build everything for a given phil
        void eventsForPhil(int phil, dd_edge &e);

    private:
        minterm mint;
        int nPhils;
        // int sz;
        forest& mxd;

        int ph;
        int rf;
        int lf;

        const model2var& M2V;
};


philsModel::philsModel(int nP, const model2var& m2v, forest& _mxd)
: mint(_mxd.getDomain(), RELATION), mxd(_mxd), M2V(m2v)
{
    nPhils = nP;

    for (unsigned i=mint.getNumVars(); i; --i)
    {
        mint.setVars(i, DONT_CARE, DONT_CHANGE);
    }

    // we always set the minterm back when we're done :^)
}

philsModel::~philsModel()
{
}

void philsModel::Idle2WaitBoth(dd_edge &e)
{
    /* I(ph) -> WB(ph) */
    mint.setVars(ph, 0, 1);
    mint.buildFunction(false, e);
    mint.setVars(ph, DONT_CARE, DONT_CHANGE);
}

void philsModel::WaitBoth2HaveRight(dd_edge &e)
{
    /* WB(ph) -> HR(ph), A(rf) -> NA(rf) */
    mint.setVars(ph, 1, 3);
    mint.setVars(rf, 0, 1);
    mint.buildFunction(false, e);
    mint.setVars(ph, DONT_CARE, DONT_CHANGE);
    mint.setVars(rf, DONT_CARE, DONT_CHANGE);
}

void philsModel::WaitBoth2HaveLeft(dd_edge &e)
{
    /* WB(ph) -> HR(ph), A(lf) -> NA(lf) */
    mint.setVars(lf, 0, 1);
    mint.setVars(ph, 1, 2);
    mint.buildFunction(false, e);
    mint.setVars(lf, DONT_CARE, DONT_CHANGE);
    mint.setVars(ph, DONT_CARE, DONT_CHANGE);
}

void philsModel::HaveRight2Eat(dd_edge &e)
{
    /* HR(ph) -> E(ph), A(lf) -> NA(lf) */
    mint.setVars(lf, 0, 1);
    mint.setVars(ph, 3, 4);
    mint.buildFunction(false, e);
    mint.setVars(lf, DONT_CARE, DONT_CHANGE);
    mint.setVars(ph, DONT_CARE, DONT_CHANGE);
}

void philsModel::HaveLeft2Eat(dd_edge &e)
{
    /* HL(ph) -> E(ph), A(rf) -> NA(rf) */
    mint.setVars(ph, 2, 4);
    mint.setVars(rf, 0, 1);
    mint.buildFunction(false, e);
    mint.setVars(ph, DONT_CARE, DONT_CHANGE);
    mint.setVars(rf, DONT_CARE, DONT_CHANGE);
}

void philsModel::Eat2Idle(dd_edge &e)
{
    /* E(ph) -> I(ph), NA(rf) -> A(rf), NA(lf) -> A(lf) */
    mint.setVars(lf, 1, 0);
    mint.setVars(ph, 4, 0);
    mint.setVars(rf, 1, 0);
    mint.buildFunction(false, e);
    mint.setVars(lf, DONT_CARE, DONT_CHANGE);
    mint.setVars(ph, DONT_CARE, DONT_CHANGE);
    mint.setVars(rf, DONT_CARE, DONT_CHANGE);
}

void philsModel::eventsForPhil(int phil, dd_edge &e)
{
    setPhilosopher(phil);
    dd_edge temp(&mxd);
    Idle2WaitBoth(e);
    WaitBoth2HaveRight(temp);
    e += temp;
    WaitBoth2HaveLeft(temp);
    e += temp;
    HaveRight2Eat(temp);
    e += temp;
    HaveLeft2Eat(temp);
    e += temp;
    Eat2Idle(temp);
    e += temp;
}

// **********************************************************************
// Helper functions
// **********************************************************************

void printStats(const char* who, const forest* f)
{
    printf("%s stats:\n", who);
    FILE_output mout(stdout);
    f->reportStats(mout, "\t",
            HUMAN_READABLE_MEMORY  |
            BASIC_STATS | EXTRA_STATS |
            STORAGE_STATS | STORAGE_DETAILED |
            HOLE_MANAGER_STATS | HOLE_MANAGER_DETAILED
    );
}

variable** initializeVariables(const model2var &M2V)
{
    // set bounds for each variable
    // forks: 2 states, philosophers: 5 states
    variable** vars = new variable* [1+M2V.numLevels()];
    for (int i=0; i<=M2V.numLevels(); i++) {
        vars[i] = nullptr;
    }

#ifdef NAME_VARIABLES
    char buffer[32];
#endif

    for (int i=0; i<M2V.numPhils(); i++) {
        /*
           Create ith fork
        */
        int v = M2V.forkVar(i);

        if (vars[v]) {
            fprintf(stderr, "Error: forkVar(%d) = %d index in use\n", i, v);
            exit(1);
        }

        std::string name;
#ifdef NAME_VARIABLES
        buffer[0] = 0;
        snprintf(buffer, 32, "fork%d", i);
        name = buffer;
#endif
        vars[v] = createVariable(2, name);

        /*
           Create ith philosopher
        */
        v = M2V.philVar(i);

        if (vars[v]) {
            fprintf(stderr, "Error: philVar(%d) = %d index in use\n", i, v);
            exit(1);
        }

#ifdef NAME_VARIABLES
        buffer[0] = 0;
        snprintf(buffer, 32, "phil%d", i);
        name = buffer;
#endif
        vars[v] = createVariable(5, name);

    } // for i

    return vars;
}


void buildInitialState(dd_edge &init)
{
    minterm mt(init.getForest());
    mt.setAllVars(0);
    mt.buildFunction(false, init);
}



// Test Index Set
void testIndexSet(const dd_edge& mdd, dd_edge& indexSet)
{
    apply(CONVERT_TO_INDEX_SET, mdd, indexSet);
    FILE_output mout(stdout);
    double card = 0;
    apply(CARDINALITY, indexSet, card);
    mout << "Index set has cardinality " << card << "\n";
    indexSet.showGraph(mout);
}

domain* runWithOptions(int nPhilosophers, const switches &sw, logger* LOG)
{
    timer start;
    FILE_output meddlyout(stdout);

    // Number of levels in domain (excluding terminals)
    // int nLevels = nPhilosophers * 2;


    const char* order_description = 0;
    switch (sw.vord) {
        case pfpf:  order_description = "phil, fork, phil, fork, ...";
                    break;
        case fpfp:  order_description = "fork, phil, fork, phil, ...";
                    break;
        case ppff:  order_description = "phil, ..., phil, fork, ..., fork";
                    break;
        case ffpp:  order_description = "fork, ..., fork, phil, ..., phil";
                    break;
        default:    order_description = "unknown order?";
    }
    printf("Using variable order (from top down): %s\n", order_description);
    model2var M2V(sw.vord, nPhilosophers);

    // Set up arrays bounds based on nPhilosophers
    variable** vars = initializeVariables(M2V);

    printf("Initiailzing forests\n");

    // Create a domain and set up the state variables.
    domain *d = domain::create(vars, M2V.numLevels());
    assert(d);

    // Set up MDD options
    policies pmdd(false), pmxd(true);

    if (sw.mark_sweep) {
        pmdd.useReferenceCounts = false;
        pmxd.useReferenceCounts = false;
        printf("Using mark and sweep\n");
    } else {
        pmdd.useReferenceCounts = true;
        pmxd.useReferenceCounts = true;
        if (sw.pessimistic) {
            printf("Using pessimistic node deletion\n");
            pmdd.setPessimistic();
            pmxd.setPessimistic();
        } else {
            printf("Using optimistic node deletion\n");
            pmdd.setOptimistic();
            pmxd.setOptimistic();
        }
    }

    // Create an MDD forest in this domain (to store states)
    forest* mdd = forest::create(d, SET, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL, pmdd);
    assert(mdd);
    mdd->setLogger(LOG, "MDD");

    // Create a MXD forest in domain (to store transition diagrams)
    forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL, pmxd);
    assert(mxd);
    mxd->setLogger(LOG, "MxD");

    // Set up initial state array based on nPhilosophers
    if (LOG) LOG->newPhase(mdd, "Building initial state");
    dd_edge initialStates(mdd);
    buildInitialState(initialStates);

    if (LOG) LOG->newPhase(mxd, "Building next-state function");
    printf("Building next-state function for %d dining philosophers\n",
            nPhilosophers);
    fflush(stdout);
    start.note_time();

    // Create a matrix diagram to represent the next-state function
    // Next-State function is computed by performing a union of next-state
    // functions that each represent how a philosopher's state can change.
    philsModel model(nPhilosophers, M2V, *mxd);
    dd_edge nsf(mxd);
    pregen_relation* ensf = nullptr;
    saturation_operation* sat = nullptr;

    if ('s' == sw.method) {
        ensf = new pregen_relation(mxd, 6*nPhilosophers);
    }
    if ('k' == sw.method) {
        ensf = new pregen_relation(mxd);
    }

    if (ensf) {
        dd_edge temp(mxd);
        for (int i = 0; i < nPhilosophers; i++) {
            model.setPhilosopher(i);
            model.Idle2WaitBoth(temp);
            ensf->addToRelation(temp);
            model.WaitBoth2HaveRight(temp);
            ensf->addToRelation(temp);
            model.WaitBoth2HaveLeft(temp);
            ensf->addToRelation(temp);
            model.HaveRight2Eat(temp);
            ensf->addToRelation(temp);
            model.HaveLeft2Eat(temp);
            ensf->addToRelation(temp);
            model.Eat2Idle(temp);
            ensf->addToRelation(temp);
        }
        ensf->finalize();
    } else {
        dd_edge phil(mxd);
        for (int i = 0; i < nPhilosophers; i++) {
            model.eventsForPhil(i, phil);
            nsf += phil;
        }
    }
    start.note_time();
    printf("Next-state function construction took %.4e seconds\n",
            start.get_last_seconds());
    if (!ensf) {
        printf("Next-state function MxD has\n\t%lu nodes\n\t%lu edges\n",
                nsf.getNodeCount(), nsf.getEdgeCount());
    }

    //
    // Show stats for nsf construction
    //
    printStats("MxD", mxd);

#ifdef SHOW_MXD
    printf("Next-State Function:\n");
    nsf.showGraph(meddlyout);
#endif


    //
    // Build reachable states
    //
    if (LOG) LOG->newPhase(mdd, "Building reachability set");
    dd_edge reachableStates(initialStates);
    start.note_time();

    switch (sw.method) {
        case 'b':
            printf("Building reachability set using traditional algorithm\n");
            fflush(stdout);
            apply(REACHABLE_STATES_BFS, initialStates, nsf, reachableStates);
            break;

        case 'm':
            printf("Building reachability set using saturation, monolithic relation\n");
            fflush(stdout);
            apply(REACHABLE_STATES_DFS, initialStates, nsf, reachableStates);
            break;

        case 'k':
        case 's':
            printf("Building reachability set using saturation, relation");
            if ('k'==sw.method) printf(" by levels\n");
            else                printf(" by events\n");
            fflush(stdout);
            sat = SATURATION_FORWARD(mdd, ensf, mdd);
            if (0==sat) {
                throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
            }
            sat->compute(initialStates, reachableStates);
            break;

        default:
            printf("Error - unknown method\n");
            exit(2);
    };
    start.note_time();
    printf("Done\n");

    printf("Reachability set construction took %.4e seconds\n",
            start.get_last_seconds() );
    fflush(stdout);
    printf("#Nodes: %lu\n", reachableStates.getNodeCount());
    printf("#Edges: %lu\n", reachableStates.getEdgeCount());

#ifdef SHOW_MDD
    printf("Reachability set:\n");
    reachableStates.showGraph(meddlyout);
#endif

    // Show stats for rs construction
    printStats("MDD", mdd);

    compute_table::showAll(meddlyout, 3);

    double c = 0;
    apply(CARDINALITY, reachableStates, c);
    printf("Approximately %e reachable states\n", c);
    fflush(stdout);

#if HAVE_LIBGMP
    if (sw.exact) {
        mpz_t nrs;
        mpz_init(nrs);
        apply(CARDINALITY, reachableStates, nrs);
        printf("Exactly ");
        mpz_out_str(0, 10, nrs);
        printf(" reachable states\n");
        fflush(stdout);
        mpz_clear(nrs);
    }
#endif

    if (sw.printReachableStates) {
#ifdef OLD_ITERATORS
        // Create a EV+MDD forest in this domain (to store index set)
        forest* evplusmdd =
            forest::create(d, false, range_type::INTEGER, edge_labeling::INDEX_SET);
        assert(evplusmdd != NULL);

        // Test Convert MDD to Index Set EV+MDD
        dd_edge indexSet(evplusmdd);
        testIndexSet(reachableStates, indexSet);
        int* element = (int *) malloc((nLevels + 1) * sizeof(int));

        double cardinality;
        apply(CARDINALITY, indexSet, cardinality);
        for (int index = 0; index < int(cardinality); index++)
        {
            evplusmdd->getElement(indexSet, index, element);
            printf("Element at index %d: [ ", index);
            for (int i = nLevels; i > 0; i--)
            {
                printf("%d ", element[i]);
            }
            printf("]\n");
        }
#else
        FILE_output out(stdout);
        long c=0;
        for (auto I = reachableStates.begin(); I; ++I) {
            out << "State at index " << c << ": ";
            (*I).show(out);
            out << "\n";
            ++c;
        }
        out << c << " reachable states\n";
#endif
    }


    if (sw.build_pdf) {
        reachableStates.setLabel("reachable");
        dot_maker dm(mdd, "out");
        dm.addRootEdge(reachableStates);
        dm.doneGraph();
        dm.runDot("pdf");
    }

    /*
       Next will be cleanup.
       Log that.
    */
    if (LOG) {
        LOG->newPhase(mdd, "Cleanup");
        LOG->newPhase(mxd, "Cleanup");
    }

    printf("Destroying dd_edges\n");
    fflush(stdout);
    return d;
}



int usage(const char* who)
{
    /* Strip leading directory, if any: */
    const char* name = who;
    for (const char* ptr=who; *ptr; ptr++) {
        if ('/' == *ptr) name = ptr+1;
    }

    printf("\nUsage: %s [options]\n\n", name);

    printf("\t-n<#phils>: set number of philosophers\n");
    printf("\n");

    printf("\t-exact:     display the exact number of states\n");
    printf("\t-cs<cache>: set cache size (0 for library default)\n");
    printf("\n");

    printf("\t-ms:        use mark and sweep node deletion\n");
    printf("\t-opt:       use optimistic node deletion (default)\n");
    printf("\t-pess:      use pessimistic node deletion (lower mem usage)\n");
    printf("\n");

    printf("\t-bfs:       use traditional iterations (default)\n\n");
    printf("\t-dfs:       use fastest saturation (currently, -msat)\n");
    printf("\t-esat:      use saturation by events\n");
    printf("\t-ksat:      use saturation by levels\n");
    printf("\t-msat:      use monolithic saturation\n");
    printf("\n");

    printf("\t-l lfile:   Write logging information to specified file\n");
    printf("\t-pdf:       Write MDD for reachable states to out.pdf\n");
    printf("\n");

    printf("\t-ofpfp:     (default) Variable order: fork, phil, fork, phil, ...\n");
    printf("\t-opfpf:     Variable order: phil, fork, phil, fork, ...\n");
    printf("\t-oppff:     Variable order: phil, ..., phil, fork, ..., fork\n");
    printf("\t-offpp:     Variable order: fork, ..., fork, phil, ..., phil\n");
    printf("\n");

    printf("\t-print:     Print the reachable states\n");
    printf("\n");
    return 0;
}


int main(int argc, char *argv[])
{
    int nPhilosophers = 0; // number of philosophers
    switches sw;
    int cacheSize = 0;
    const char* lfile = 0;

    for (int i=1; i<argc; i++) {
        const char* cmd = argv[i];
        if (strcmp(cmd, "-opt") == 0) {
            sw.mark_sweep = false;
            sw.pessimistic = false;
            continue;
        }
        if (strcmp(cmd, "-pess") == 0) {
            sw.mark_sweep = false;
            sw.pessimistic = true;
            continue;
        }
        if (strcmp(cmd, "-ms") == 0) {
            sw.mark_sweep = true;
            continue;
        }
        if (strncmp(cmd, "-cs", 3) == 0) {
            cacheSize = strtol(&cmd[3], NULL, 10);
            if (cacheSize < 1) {
                return 1+usage(argv[0]);
            }
            continue;
        }
        if (strcmp(cmd, "-exact") == 0) {
            sw.exact = true;
            continue;
        }
        if (strcmp(cmd, "-print") == 0) {
            sw.printReachableStates = true;
            continue;
        }
        if (strncmp(cmd, "-n", 2) == 0) {
            nPhilosophers = strtol(cmd+2, NULL, 10);
            if (nPhilosophers < 1) {
                return 1+usage(argv[0]);
            }
            continue;
        }

        if (strcmp(cmd, "-bfs") == 0) {
            sw.method = 'b';
            continue;
        }

        if (strcmp(cmd, "-dfs") == 0) {
            sw.method = 'm';
            continue;
        }

        if (strcmp(cmd, "-msat") == 0) {
            sw.method = 'm';
            continue;
        }
        if (strcmp(cmd, "-esat") == 0) {
            sw.method = 's';
            continue;
        }
        if (strcmp(cmd, "-ksat") == 0) {
            sw.method = 'k';
            continue;
        }
        if (strcmp("-l", argv[i])==0) {
            lfile = argv[i+1];
            i++;
            continue;
        }
        if (strcmp("-pdf", argv[i])==0) {
            sw.build_pdf = true;
            continue;
        }

        if (strcmp(cmd, "-ofpfp") == 0) {
            sw.vord = fpfp;
            continue;
        }
        if (strcmp(cmd, "-opfpf") == 0) {
            sw.vord = pfpf;
            continue;
        }
        if (strcmp(cmd, "-oppff") == 0) {
            sw.vord = ppff;
            continue;
        }
        if (strcmp(cmd, "-offpp") == 0) {
            sw.vord = ffpp;
            continue;
        }

        return 1+usage(argv[0]);
    }

    while (nPhilosophers < 2) {
        printf("Enter the number of philosophers (at least 2): ");
        scanf("%d", &nPhilosophers);
    }

    // Initialize MEDDLY

    initializer_list* L = defaultInitializerList(0);
    if (cacheSize > 0) {
        ct_initializer::setMaxSize(cacheSize);
    }
    MEDDLY::initialize(L);

    //
    // Set up logger, if any
    //

    std::ofstream lstrm;
    ostream_output log(lstrm);
    logger* LOG = 0;
    if (lfile) {
        lstrm.open(lfile, std::ofstream::out);
        if (!lstrm) {
            printf("Couldn't open %s for writing, no logging\n", lfile);
        } else {
            LOG = new simple_logger(log);
            LOG->recordNodeCounts();
            LOG->recordTimeStamps();
            char comment[80];
            snprintf(comment, 80, "Automatically generated by dining_phils (N=%d)",
                    nPhilosophers);
            LOG->addComment(comment);
        }
    }

    try {
        domain* d = runWithOptions(nPhilosophers, sw, LOG);
        printf("Cleaning up\n");
        fflush(stdout);
        if (LOG) {
            domain::destroy(d);
            delete LOG;
        }
        MEDDLY::cleanup();
        printf("\n\nDONE\n");
        return 0;
    }
    catch (error e) {
        printf("Caught MEDDLY error: %s\n", e.getName());
        return 1;
    }
}


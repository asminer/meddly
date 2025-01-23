
#include "../src/meddly.h"
#include "simple_model.h"

#include <assert.h>

#define VERBOSE

using namespace MEDDLY;

const char* mdl[] = { "X-+." };

int main()
{
    try {
        MEDDLY::initialize();

        // Initialize domain
        const int varsizes[] = { 3, 3, 3};
        domain* d = domain::createBottomUp(varsizes, 3);
        assert(d);

        // Initialize forests
        forest* mdd = forest::create(d, SET, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);
        forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);

        // build initial state
        minterm init_mt(mdd);
        dd_edge init_state(mdd);
        init_mt.setVar(1, 2);
        init_mt.setVar(2, 0);
        init_mt.setVar(3, 0);
        init_mt.buildFunction(false, init_state);

        // build next-state function
        dd_edge nsf(mxd);
        buildNextStateFunction(mdl, 1, mxd, nsf);

        minterm mask_mt(mxd);
        dd_edge mask(mxd);
        mask_mt.setVars(3, 0, 0);
        mask_mt.setVars(2, DONT_CARE, DONT_CARE);
        mask_mt.setVars(1, DONT_CARE, DONT_CARE);
        mask_mt.buildFunction(false, mask);
        nsf *= mask;

        // build rs using traditional & saturation
        dd_edge reachable1(mdd);
        dd_edge reachable2(mdd);
        apply(REACHABLE_STATES_DFS, init_state, nsf, reachable1);
        apply(REACHABLE_STATES_BFS, init_state, nsf, reachable2);

        // Display everything
#ifdef VERBOSE
        FILE_output out(stdout);
        printf("Initial state:\n");
        init_state.showGraph(out);

        printf("NSF:\n");
        nsf.showGraph(out);

        printf("DFS states:\n");
        reachable1.showGraph(out);

        printf("BFS states:\n");
        reachable2.showGraph(out);
#endif

        int retval;
        if (reachable1 == reachable2) {
            retval = 0;
        } else {
            retval = 1;
        }

#ifdef VERBOSE
        if (retval) {
            printf("\nReachable states DO NOT match\n\n");
        } else {
            printf("\nReachable states match\n\n");
        }
#endif

        // cleanup
        MEDDLY::cleanup();
        return retval;
    }
    catch (MEDDLY::error e) {
        fprintf(stderr, "Caught error %s, thrown in %s line %u\n",
                e.getName(), e.getFile(), e.getLine());
        return 1;
    }
}


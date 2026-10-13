#include <iostream>
#include "../src/meddly.h"

using namespace MEDDLY;

int main() {
    MEDDLY::initialize();

    // 1. Domain with 2 variables {0, 1}
    int bounds[2] = {2, 2};
    domain* d = domain::createBottomUp(bounds, 2);

    // 2. Identity-reduced relation forest (unweighted MULTI_TERMINAL)
    policies p(true);
    p.setIdentityReduced();
    forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
                                 edge_labeling::MULTI_TERMINAL, p);

    // 3. Create a "by levels" pregen_relation
    pregen_relation* rel = new pregen_relation(mxd);

    // 4. Add a transition: variable 2 transitions 0 -> 1, variable 1 unchanged
    minterm mt(mxd);
    mt.setVars(2, 0, 1);
    mt.setVars(1, DONT_CARE, DONT_CHANGE);

    dd_edge e(mxd);
    mt.buildFunction(false, e);
    rel->addToRelation(e);

    // 5. Finalize triggers splitMxd(SplitSubtract)
    std::cout << "Calling rel->finalize()..." << std::endl;
    rel->finalize(); // <-- SIGSEGV at sat_relations.cc:232

    std::cout << "Completed successfully!" << std::endl;

    delete rel;
    MEDDLY::cleanup();
    return 0;
}


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

#include <iostream>
#include "meddly/meddly.h"

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


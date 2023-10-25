/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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

#include "../src/meddly.h"
#include <iostream>

using namespace MEDDLY;

void checkCount(unsigned long expected, const dd_edge &e)
{
    using namespace std;

    cout << "Checking ref count (node " << e.getNode() << ", expecting " << expected << "): ";

    const expert_forest* f = dynamic_cast <expert_forest*>(e.getForest());
    if (!f) throw "Bad forest";

    const unsigned long cnt = f->getNodeInCount(e.getNode());
    cout << cnt << endl;

    if (cnt != expected) throw "Bad count";
}

int main()
{
    try {
        MEDDLY::initialize();

        ostream_output s(std::cout);

        const int bounds[] = { 5, 5, 5 };

        domain* d = createDomainBottomUp(bounds, 3);
        policies p;
        p.useDefaults(false);
        p.useReferenceCounts = true;
        p.setPessimistic();

        forest* mdd = d->createForest(false, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL, p);

        dd_edge x1(mdd), x2(mdd), x3(mdd), x123(mdd), tmp(mdd);

        mdd->createEdgeForVar(1, false, x1);
        mdd->createEdgeForVar(2, false, x2);
        mdd->createEdgeForVar(3, false, x3);

        checkCount(1, x1);
        checkCount(1, x2);
        checkCount(1, x3);

        x123 = x1;

        checkCount(2, x1);
        checkCount(2, x123);

        x123 *= x2;

        x123.showGraph(s);

        checkCount(5, x1);

        x123 *= x3;

        checkCount(1, x123);

        tmp = x1;

        checkCount(6, x1);

        tmp *= x123;

        checkCount(2, x123);

        checkCount(5, x1);

        mdd->createEdge(false, tmp);

        checkCount(1, x123);

        mdd->createEdge(false, x123);

        checkCount(1, x1);

        MEDDLY::cleanup();
        return 0;
    }
    catch (MEDDLY::error e) {
        std::cerr << "Caught meddly error '" << e.getName() << "'\n";
        std::cerr << "    thrown in " << e.getFile() << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        std::cerr << "Caught our own error: " << e << "\n";
        return 2;
    }
    std::cerr << "Some other error?\n";
    return 3;
}

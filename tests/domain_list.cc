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

void flip(MEDDLY::domain* &d)
{
    const int bounds[] = { 5, 5, 5, 5, 5 };
    if (d) {
        domain::destroy(d);
    } else {
        d = domain::createBottomUp(bounds, 5);
    }
}

void test_domain_list()
{
    domain* D[10];
    for (unsigned i=0; i<10; i++) {
        D[i] = nullptr;
    }

    // test the list
    const char* pi = "3.1415926535 8979323846 2643383279 5028841971 6939937510 5820974944 5923078164 0628620899 8628034825";

    for (unsigned i=0; pi[i]; i++) {
        if ('.' == pi[i]) continue;
        if (' ' == pi[i]) continue;

        flip(D[pi[i]-'0']);
    }

    std::cout << "Non-null domains:\n";
    for (unsigned i=0; i<10; i++) {
        if (D[i])
            std::cout << "  #" << i << "\n";
    }

    domain::testMarkAllDomains(true);

    for (unsigned i=0; i<10; i++) {
        if (D[i]) {
            if (! D[i]->isMarkedForDeletion()) throw "should be marked";
        }
    }

    domain::testMarkAllDomains(false);

    for (unsigned i=0; i<10; i++) {
        if (D[i]) {
            if (D[i]->isMarkedForDeletion()) throw "should not be marked";
        }
    }

}

int main()
{
    using namespace std;
    try {
        MEDDLY::initialize();
        test_domain_list();
        MEDDLY::cleanup();
        cout << "Passed\n";
        return 0;
    }
    catch (MEDDLY::error e) {
        cerr << "Caught meddly error '" << e.getName() << "'\n";
        cerr << "    thrown in " << e.getFile() << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        cerr << "Caught our own error: " << e << "\n";
        return 2;
    }
    cerr << "Some other error?\n";
    return 3;
}


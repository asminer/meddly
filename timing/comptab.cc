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

#include <iostream>
#include <iomanip>

#include "../src/meddly.h"
#include "timer.h"
#include "reporting.h"

using namespace MEDDLY;

const unsigned LEVELS = 128;
const unsigned DOTS = 16;

#ifdef DEVELOPMENT_CODE
const unsigned CREATIONS = 32*1024;
const unsigned RCREATE   = 8*1024;
#else
const unsigned CREATIONS = 1024*1024;
const unsigned RCREATE   = 256*1024;
#endif

// Request and recycle lots of compute table keys
void reqKeys(const ct_entry_type *CTE)
{
    static ct_entry_key* keys[LEVELS];

    using namespace std;
    cout << "\nTiming test: CT key creation / recycling\n\n";
    cout << "Requesting " << setw(10) << DOTS * CREATIONS * LEVELS << " keys ";

    timer T;

    for (unsigned i=0; i<DOTS; i++) {
        cout << ".";
        cout.flush();
        for (unsigned j=0; j<CREATIONS; j++) {
            unsigned k = 0;
            while (k<LEVELS) {
                keys[k++] = compute_table::useEntryKey(CTE, 0);
            }
            while (k) {
                compute_table::recycle(keys[--k]);
            }
        } // for j
    } // for i

    T.note_time();
    show_sec(std::cout, T, 3, 3);

    if (startReport(T, __FILE__, "req key")) {
        report  << "Requested " << DOTS * CREATIONS * LEVELS << " keys"
                << std::endl;
    }
}

// Request, build, and recycle lots of compute table keys
void makeKeys(const ct_entry_type *CTE)
{
    static ct_entry_key* keys[LEVELS];

    using namespace std;
    cout << "\nTiming test: CT key building\n\n";
    cout << "Building   " << setw(10) << DOTS * CREATIONS * LEVELS << " keys ";

    timer T;

    for (unsigned i=0; i<DOTS; i++) {
        cout << ".";
        cout.flush();
        for (unsigned j=0; j<CREATIONS; j++) {
            unsigned k = 0;
            while (k<LEVELS) {
                keys[k] = compute_table::useEntryKey(CTE, 0);
                keys[k]->writeI(int(k));
                keys[k]->writeN(i);
                keys[k]->writeN(j);
                k++;
            }
            while (k) {
                compute_table::recycle(keys[--k]);
            }
        } // for j
    } // for i

    T.note_time();
    show_sec(std::cout, T, 3, 3);

    if (startReport(T, __FILE__, "make key")) {
        report  << "Built " << DOTS * CREATIONS * LEVELS << " keys"
                << std::endl;
    }
}

// Request, build, and recycle lots of 'repeating' compute table keys
void makeRKeys(const ct_entry_type *CTE)
{
    static ct_entry_key* keys[LEVELS];

    using namespace std;
    cout << "\nTiming test: Repeating CT key building\n\n";
    cout << "Building   " << setw(10) << DOTS * RCREATE * LEVELS << " keys ";

    timer T;

    for (unsigned i=0; i<DOTS; i++) {
        cout << ".";
        cout.flush();
        for (unsigned j=0; j<RCREATE; j++) {
            unsigned k = 0;
            while (k<LEVELS) {
                const unsigned repeats = k % 16;
                keys[k] = compute_table::useEntryKey(CTE, repeats);
                keys[k]->writeI(int(k));
                keys[k]->writeN(i);
                for (unsigned r=0; r<repeats; r++) {
                    keys[k]->writeI(int(r));
                    keys[k]->writeN(j);
                }
                k++;
            }
            while (k) {
                compute_table::recycle(keys[--k]);
            }
        } // for j
    } // for i

    T.note_time();
    show_sec(std::cout, T, 3, 3);

    if (startReport(T, __FILE__, "make rkey")) {
        report  << "Built " << DOTS * RCREATE * LEVELS << " repeating keys"
                << std::endl;
    }
}

int main(int argc, const char** argv)
{
    try {
        setReport(argc, argv);
        MEDDLY::initialize();
        forest* F = nullptr;

        // ct_entry_type cte("test_entry", "INN:N");
        ct_entry_type* cte = new ct_entry_type("test_entry");
        cte->setFixed('I', F, F);
        cte->setResult(F);

        // ct_entry_type ctr("test_repeating", "IN.IN:N");
        ct_entry_type* ctr = new ct_entry_type("test_repeating");
        ctr->setFixed('I', F);
        ctr->setRepeat('I', F);
        ctr->setResult(F);

        reqKeys(cte);
        makeKeys(cte);
        makeRKeys(ctr);

        MEDDLY::cleanup();

        std::cout << "\n";
        return 0;
    }
    catch (MEDDLY::error e) {
        std::cerr   << "\nCaught meddly error '" << e.getName()
                    << "'\n    thrown in " << e.getFile()
                    << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        std::cerr << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    std::cerr << "\nSome other error?\n";
    return 3;

}


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
const unsigned CREATIONS = 16*1024;
const unsigned RCREATE   = 8*1024;
#else
const unsigned CREATIONS = 512*1024;
const unsigned RCREATE   =  64*1024;
#endif

// #define USE_NEW_CT_INTERFACE

int global_i;
int global_j;

void _reqKeys(const ct_entry_type* CTE, unsigned L)
{
    if (0==L) return;

#ifdef USE_NEW_CT_INTERFACE
    ct_vector key(CTE->getKeySize());
#else
    ct_entry_key* key = compute_table::useEntryKey(CTE, 0);
#endif

    _reqKeys(CTE, --L);

#ifndef USE_NEW_CT_INTERFACE
    compute_table::recycle(key);
#endif
}


// Request and recycle lots of compute table keys
void reqKeys(const ct_entry_type *CTE)
{
    using namespace std;
    cout << "\nTiming test: CT key creation / recycling\n\n";
    cout << "Requesting " << setw(10) << DOTS * CREATIONS * LEVELS << " keys ";

    timer T;

    for (global_i=0; global_i<DOTS; global_i++) {
        cout << ".";
        cout.flush();
        for (global_j=0; global_j<CREATIONS; global_j++) {
            _reqKeys(CTE, LEVELS);
        } // for j
    } // for i

    T.note_time();
    show_sec(std::cout, T, 3, 3);

    if (startReport(T, __FILE__, "req key")) {
        report  << "Requested " << DOTS * CREATIONS * LEVELS << " keys"
                << std::endl;
    }
}

void _makeKeys(const ct_entry_type* CTE, unsigned L)
{
    if (0==L) return;
#ifdef USE_NEW_CT_INTERFACE
    ct_vector key(CTE->getKeySize());
    key[0].setI( int(L) );
    key[1].setN(global_i);
    key[2].setN(global_j);
#else
    ct_entry_key* key = compute_table::useEntryKey(CTE, 0);
    key->writeI( int(L) );
    key->writeN(global_i);
    key->writeN(global_j);
#endif


    _makeKeys(CTE, --L);

#ifndef USE_NEW_CT_INTERFACE
    compute_table::recycle(key);
#endif
}

// Request, build, and recycle lots of compute table keys
void makeKeys(const ct_entry_type *CTE)
{
    using namespace std;
    cout << "\nTiming test: CT key building\n\n";
    cout << "Building   " << setw(10) << DOTS * CREATIONS * LEVELS << " keys ";

    timer T;

    for (global_i=0; global_i<DOTS; global_i++) {
        cout << ".";
        cout.flush();
        for (global_j=0; global_j<CREATIONS; global_j++) {
            _makeKeys(CTE, LEVELS);
        } // for j
    } // for i

    T.note_time();
    show_sec(std::cout, T, 3, 3);

    if (startReport(T, __FILE__, "make key")) {
        report  << "Built " << DOTS * CREATIONS * LEVELS << " keys"
                << std::endl;
    }
}

void _makeRKeys(const ct_entry_type* CTE, unsigned L)
{
    if (0==L) return;
    const unsigned repeats = L % 16;
#ifdef USE_NEW_CT_INTERFACE
    ct_vector key(CTE->getKeySize(repeats));
    unsigned slot = 0;
    key[slot++].setI(int(L));
    key[slot++].setN(global_i);
    for (unsigned r=0; r<repeats; r++) {
        key[slot++].setI(int(r));
        key[slot++].setN(global_j);
    }
#else
    ct_entry_key* key = compute_table::useEntryKey(CTE, repeats);
    key->writeI(int(L));
    key->writeN(global_i);
    for (unsigned r=0; r<repeats; r++) {
        key->writeI(int(r));
        key->writeN(global_j);
    }

#endif

    _makeRKeys(CTE, --L);

#ifndef USE_NEW_CT_INTERFACE
    compute_table::recycle(key);
#endif
}

// Request, build, and recycle lots of 'repeating' compute table keys
void makeRKeys(const ct_entry_type *CTE)
{
    using namespace std;
    cout << "\nTiming test: Repeating CT key building\n\n";
    cout << "Building   " << setw(10) << DOTS * RCREATE * LEVELS << " keys ";

    timer T;

    for (global_i=0; global_i<DOTS; global_i++) {
        cout << ".";
        cout.flush();
        for (global_j=0; global_j<RCREATE; global_j++) {
            _makeRKeys(CTE, LEVELS);
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


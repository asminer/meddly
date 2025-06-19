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

#include "defines.h"
#include "initializer.h"
#include "error.h"
#include "domain.h"
#include "relforest.h"
#include "oper.h"
#include "unpacked_node.h"
#include "ct_initializer.h"
#include "ct_vector.h"
#include "ops_builtin.h"
#include "memory_managers/init_managers.h"
#include "storage/init_storage.h"
#include "init_forests.h"
#include "revision.h"

// #define DEBUG_INITLIST

// ******************************************************************
// *                                                                *
// *                    initializer_list methods                    *
// *                                                                *
// ******************************************************************

MEDDLY::initializer_list* MEDDLY::initializer_list::meddlyInitializers;
bool MEDDLY::initializer_list::isRunning = false;

MEDDLY::initializer_list::initializer_list(initializer_list* prev)
{
    previous = prev;
}

MEDDLY::initializer_list::~initializer_list()
{
    // DON'T delete previous
}

void MEDDLY::initializer_list::initializeLibrary(initializer_list* L)
{
    if (libraryIsRunning()) {
        throw error(error::ALREADY_INITIALIZED, __FILE__, __LINE__);
    }

    // Hard-coded statics
    memstats::initGlobalStats();
    operation::initializeStatics();
    unpacked_node::initStatics();
    domain::initDomList();
    relforest::initStatics();
    ct_vector::initStatics();
    ct_entry_type::initStatics();

    // Reverse the list
    initializer_list* reverse = nullptr;
    while (L) {
        initializer_list* next = L->previous;
        L->previous = reverse;
        reverse = L;
        L = next;
    }
    // Run through intializers
    for (initializer_list* curr = reverse; curr; curr = curr->previous) {
#ifdef DEBUG_INITLIST
        printf("setup() on object %p\n", curr);
#endif
        curr->setup();
    }
    // Reverse it back
    while (reverse) {
        initializer_list* next = reverse->previous;
        reverse->previous = L;
        L = reverse;
        reverse = next;
    }
    meddlyInitializers = L;
    isRunning = true;
}

void MEDDLY::initializer_list::cleanupLibrary()
{
    if (!libraryIsRunning()) {
        throw error(error::UNINITIALIZED, __FILE__, __LINE__);
    }

#ifdef STATS_ON_DESTROY
    if (operation::Monolithic_CT) {
        fprintf(stderr, "Compute table (before destroy):\n");
        operation::Monolithic_CT->show(stderr, false);
    }

#endif

    domain::markDomList();

    operation::destroyAllOps();

    domain::deleteDomList();
    relforest::freeStatics();

    //
    // Run through initializers
    //
    while (meddlyInitializers) {
        initializer_list* prev = meddlyInitializers->previous;
#ifdef DEBUG_INITLIST
        printf("cleanup() on object %p\n", meddlyInitializers);
#endif
        meddlyInitializers->cleanup();
        delete meddlyInitializers;
        meddlyInitializers = prev;
    }
    // ^ should we traverse this list before clobbering forests, etc?

    unpacked_node::doneStatics();
    ct_vector::doneStatics();
    ct_entry_type::doneStatics();

    isRunning = false;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::initializer_list* MEDDLY::defaultInitializerList(initializer_list* prev)
{
    prev = new memman_initializer(prev);
    prev = new ct_initializer(prev);
    prev = new storage_initializer(prev);
    prev = makeBuiltinInitializer(prev);
    prev = new forest_initializer(prev);

    return prev;
}

void MEDDLY::initialize()
{
    initializer_list::initializeLibrary( defaultInitializerList(0) );
}

void MEDDLY::initialize(initializer_list* L)
{
    initializer_list::initializeLibrary(L);
}

void MEDDLY::cleanup()
{
    initializer_list::cleanupLibrary();
}

const char* MEDDLY::getLibraryInfo(int what)
{
  static char* title = 0;
  switch (what) {
    case 0:
      if (!title) {
        title = new char[80];
        snprintf(title, 80,
#ifdef DEVELOPMENT_CODE
          "%s version %s.dev",
#else
          "%s version %s",
#endif
            PACKAGE_NAME, VERSION
        );
      }
      return title;

    case 1:
      return "Copyright (C) 2009, Iowa State University Research Foundation, Inc.";

    case 2:
      return "Released under the GNU Lesser General Public License, version 3";

    case 3:
      return PACKAGE_URL;

    case 4:
      return "Data Structures and operations available:\n\
(1) MDDs: Union, Intersection, Difference.\n\
(2) Matrix Diagrams (MXDs): Union, Intersection, Difference.\n\
(3) Multi-Terminal MDDs (MTMDDs) with integer or real terminals:\n\
    Arithmetic: Plus, Minus, Multiply, Divide, Min, Max.\n\
    Logical: <, <=, >, >=, ==, !=.\n\
    Conversion to and from MDDs.\n\
(4) Multi-Terminal MXDs (MTMXDs) with integer or real terminals:\n\
    Arithmetic: Plus, Minus, Multiply, Divide, Min, Max.\n\
    Logical: <, <=, >, >=, ==, !=.\n\
    Conversion to and from MXDs.\n\
";

    case 5:
      return MEDDLY_DATE;
  }
  return 0;
}


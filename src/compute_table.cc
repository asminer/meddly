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
#include "error.h"

#include "ct_initializer.h"
#include "compute_table.h"
#include "ct_entry_type.h"
#include "ct_entry_key.h"

// **********************************************************************
// *                                                                    *
// *                    compute_table_style  methods                    *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table_style::compute_table_style(bool um)
{
    uses_mono = um;
}

MEDDLY::compute_table_style::~compute_table_style()
{
}

MEDDLY::compute_table*
MEDDLY::compute_table_style::create(const ct_settings &s) const
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}


MEDDLY::compute_table*
MEDDLY::compute_table_style::create(const ct_settings &s,
      unsigned etid) const
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

// **********************************************************************
// *                                                                    *
// *                       compute_table  statics                       *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table* MEDDLY::compute_table::Monolithic_CT;
MEDDLY::ct_entry_key* MEDDLY::compute_table::free_keys;

MEDDLY::compute_table* MEDDLY::compute_table::front;
MEDDLY::compute_table* MEDDLY::compute_table::back;

// **********************************************************************
// *                                                                    *
// *                       compute_table  methods                       *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table::compute_table(const ct_settings &s, unsigned etid)
    : global_etid(etid)
{
    maxSize = s.maxSize;
    if (0==maxSize) {
        throw error(error::INVALID_ASSIGNMENT, __FILE__, __LINE__);
    }

    switch (s.staleRemoval) {
        case staleRemovalOption::Aggressive:
            checkStalesOnFind = true;
            checkStalesOnResize = true;
            break;
        case staleRemovalOption::Moderate:
            checkStalesOnFind = false;
            checkStalesOnResize = true;
            break;
        case staleRemovalOption::Lazy:
            checkStalesOnFind = false;
            checkStalesOnResize = false;
            break;
    }

    perf.numEntries = 0;
    perf.hits = 0;
    perf.pings = 0;
    perf.numLargeSearches = 0;
    perf.maxSearchLength = 0;
    for (unsigned i=0; i<perf.searchHistogramSize; i++) {
        perf.searchHistogram[i] = 0;
    }
    perf.resizeScans = 0;

    //
    // Add to list of CTs
    //
    next = nullptr;
    if (back) {
        back->next = this;
        prev = back;
    } else {
        front = this;
        prev = nullptr;
    }
    back = this;

#ifdef DEBUG_CLEANUP
    std::cout << "new compute table " << this << "\n";
#endif
}

MEDDLY::compute_table::~compute_table()
{
    //
    // Remove from list of CTs
    //
    if (next) {
        next->prev = prev;
    } else {
        MEDDLY_DCASSERT(back == this);
        back = prev;
    }
    if (prev) {
        prev->next = next;
    } else {
        MEDDLY_DCASSERT(front == this);
        front = next;
    }

#ifdef DEBUG_CLEANUP
    std::cout << "done compute table " << this << "\n";
#endif
}

#ifdef ALLOW_DEPRECATED_0_17_6
MEDDLY::ct_entry_key*
MEDDLY::compute_table::useEntryKey(const ct_entry_type* et, unsigned repeats)
{
            if (!et) return nullptr;
            MEDDLY_DCASSERT( (0==repeats) || et->isRepeating() );

            ct_entry_key* k;
            if (free_keys) {
                k = free_keys;
                free_keys = free_keys->next;
            } else {
                k = new ct_entry_key();
            }
            k->setup(et, repeats);
            return k;
}

void MEDDLY::compute_table::recycle(ct_entry_key* k)
{
            if (k) {
                k->next = free_keys;
                free_keys = k;
            }
}
#endif

void MEDDLY::compute_table::clearForestCTBits(std::vector <bool> &skipF) const
{
    if (global_etid) {
        const ct_entry_type* et = ct_entry_type::getEntryType(global_etid);
        if (et) et->clearForestCTBits(skipF);
    } else {
        ct_entry_type::clearAllForestCTBits(skipF);
    }
}

void MEDDLY::compute_table::sweepForestCTBits(std::vector <bool> &whichF) const
{
    if (global_etid) {
        const ct_entry_type* et = ct_entry_type::getEntryType(global_etid);
        if (et) et->sweepForestCTBits(whichF);
    } else {
        ct_entry_type::sweepAllForestCTBits(whichF);
    }
}

void MEDDLY::compute_table::updateEntry(ct_entry_key*, const ct_entry_result &)
{
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}


void MEDDLY::compute_table::countAllNodeEntries(const forest* f,
                std::vector <unsigned long> &counts)
{
    for (compute_table* ct = front; ct; ct=ct->next) {
        ct->countNodeEntries(f, counts);
    }
}

void MEDDLY::compute_table::showAll(output &s, int verbLevel)
{
    for (compute_table* ct = front; ct; ct=ct->next) {
        ct->show(s, verbLevel);
    }
}

void MEDDLY::compute_table::initStatics(const compute_table_style* ct_factory,
        const ct_settings &the_settings)
{
#ifdef ALLOW_DEPRECATED_0_17_6
    free_keys = nullptr;
#endif

    front = nullptr;
    back = nullptr;
    Monolithic_CT = nullptr;

    if (ct_factory && ct_factory->usesMonolithic()) {
        Monolithic_CT = ct_factory->create(the_settings);
#ifdef DEBUG_CLEANUP
        std::cout << "Created monolithic CT " << Monolithic_CT << "\n";
#endif
    }
}

void MEDDLY::compute_table::doneStatics()
{
#ifdef ALLOW_DEPRECATED_0_17_6
    while (free_keys) {
        ct_entry_key* n = free_keys->next;
        delete free_keys;
        free_keys = n;
    }
#endif

    if (Monolithic_CT) {
#ifdef DEBUG_CLEANUP
        std::cout << "Destroying monolithic CT " << Monolithic_CT << "\n";
#endif
//        Monolithic_CT->removeAll();
//
//        ^ This shouldn't be necessary, since we're tearing down the
//        entire library, which means we're deleting all forests.
//
        delete Monolithic_CT;
        Monolithic_CT = nullptr;
    }
}


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

#include "ct_initializer.h"
#include "compute_table.h"

// #define DEBUG_ENTRY_TYPE
// #define DEBUG_ENTRY_REGISTRY

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
// *                       compute_table  methods                       *
// *                                                                    *
// **********************************************************************


// std::vector<MEDDLY::ct_entry_type*> MEDDLY::compute_table::entryInfo;
MEDDLY::ct_entry_key* MEDDLY::compute_table::free_keys;

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

}

MEDDLY::compute_table::~compute_table()
{
}

void MEDDLY::compute_table::initStatics()
{
    free_keys = 0;
    //
    // Initialize entryInfo list
    //
    // entryInfo.resize(0);
    // entryInfo.resize(1);
    // entryInfo[0] = nullptr; // not sure if we need to reserve etid 0
}

void MEDDLY::compute_table::doneStatics()
{
    while (free_keys) {
        ct_entry_key* n = free_keys->next;
        delete free_keys;
        free_keys = n;
    }
    // delete the items?  TBD
    // entryInfo.resize(0);
}

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

/*
void MEDDLY::compute_table::unregisterOp(operation* op, unsigned num_ids)
{
    if (0==op) return;
    if (0==num_ids) return;
    unsigned stopID = op->getFirstETid()+num_ids;
    for (unsigned i=op->getFirstETid(); i<stopID; i++) {
        delete entryInfo.at(i);
        entryInfo[i] = nullptr;
    }
    //
    // absorb trailing nulls (except slot 0) in entryInfo
    //
    while (entryInfo.size() > 1) {
        if (entryInfo.back()) break;
        entryInfo.pop_back();
    }
}
*/

void MEDDLY::compute_table::updateEntry(ct_entry_key*, const ct_entry_result &)
{
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}


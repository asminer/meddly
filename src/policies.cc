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

#include "policies.h"

// ******************************************************************
// *                                                                *
// *                        policies methods                        *
// *                                                                *
// ******************************************************************

MEDDLY::policies::policies()
{
    nodemm = nullptr;   //
    nodestor = nullptr; // should cause an exception later
}

MEDDLY::policies::policies(set_or_rel rel)
{
    useDefaults(rel);
}

void MEDDLY::policies::useDefaults(set_or_rel rel)
{
    reduction = rel
        ? reduction_rule::IDENTITY_REDUCED
        : reduction_rule::FULLY_REDUCED;

    storage_flags = FULL_OR_SPARSE;

    deletion = node_deletion::OPTIMISTIC;

    // compact_min = 100;
    // compact_max = 1000000;
    // compact_frac = 40;
    // zombieTrigger = 1000000;
    // orphanTrigger = 500000;
    // compactAfterGC = false;
    // compactBeforeExpand = true;

    useReferenceCounts = true;

    // nodemm = ORIGINAL_GRID;
    nodemm = ARRAY_PLUS_GRID;
    // nodemm = MALLOC_MANAGER;
    // nodemm = HEAP_MANAGER;

    nodestor = SIMPLE_STORAGE;
    //nodestor = PATTERN_STORAGE;
    //nodestor = BEST_STORAGE;

    reorder = reordering_type::SINK_DOWN;
    swap = variable_swap_type::VAR;
}



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

#include "reach_index.h"

#define DEBUG_BASE

// **********************************************************************
// *                                                                    *
// *                     sat_index_explorer methods                     *
// *                                                                    *
// **********************************************************************

MEDDLY::sat_index_explorer::sat_index_explorer(rel_node &_RN, bool _forwd)
    : RN(_RN), forwd(_forwd)
{
    if (forwd) {
        // determine current size
        // and call expandRows
    } else {
        // determine current size
        // build everything but its the transpose
    }
}

MEDDLY::sat_index_explorer::~sat_index_explorer()
{
}

void MEDDLY::sat_index_explorer::exploreRow(unsigned i)
{

    // TBD

    finishRow(i);
}

void MEDDLY::sat_index_explorer::finishRow(unsigned)
{
    // Do nothing
}

void MEDDLY::sat_index_explorer::expandRows(unsigned newsz)
{
    // TBD
}

void MEDDLY::sat_index_explorer::finishExpandRows(unsigned, unsigned)
{
    // Do nothing
}

// **********************************************************************
// *                                                                    *
// *                             Front  end                             *
// *                                                                    *
// **********************************************************************

MEDDLY::sat_index_explorer* MEDDLY::makeSatIndexExplorer(char which)
{
    // TBD
    return nullptr;
}


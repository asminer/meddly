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

#include "oper_satur.h"
#include "error.h"
#include "forest.h"
#include "dd_edge.h"

// ******************************************************************
// *                  saturation_operation methods                  *
// ******************************************************************

MEDDLY::saturation_operation::saturation_operation(const char* name,
    unsigned et_slots, forest* arg, forest* res)
    : operation(name, et_slots)
{
    argF = arg;
    resF = res;

    registerInForest(argF);
    registerInForest(resF);
}

MEDDLY::saturation_operation::~saturation_operation()
{
    unregisterInForest(argF);
    unregisterInForest(resF);
}

bool
MEDDLY::saturation_operation::checkForestCompatibility() const
{
    auto o1 = argF->variableOrder();
    auto o2 = resF->variableOrder();
    return o1->is_compatible_with(*o2);
}

bool
MEDDLY::saturation_operation::isReachable(const dd_edge& a, const dd_edge& c)
{
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}


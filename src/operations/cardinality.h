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

#ifndef MEDDLY_CARDINALITY_H
#define MEDDLY_CARDINALITY_H

#include "../oper.h"

namespace MEDDLY {
    class unary_operation;

    /// Builder for 'cardinality' operations.
    unary_operation* CARDINALITY(forest* arg, opnd_type res);
    void CARDINALITY_init();
    void CARDINALITY_done();
};

#endif


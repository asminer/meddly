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

#ifndef MEDDLY_ARITH_MAX_H
#define MEDDLY_ARITH_MAX_H

#include "../oper.h"

/*
 * Code for element-wise arithmetic and like operations.
 */

namespace MEDDLY {
    class forest;
    class binary_operation;

    /// The 'maximum' operation builder.
    binary_operation* MAXIMUM(forest* a, forest* b, forest* c);
    void MAXIMUM_init();
    void MAXIMUM_done();
}

#endif



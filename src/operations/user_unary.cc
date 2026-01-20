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

#include "../defines.h"
#include "../ct_vector.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../forest_levels.h"
#include "user_unary.h"

MEDDLY::user_unary_operation::user_unary_operation(forest* arg, forest* res,
        user_defined_unary _F) : unary_operation(arg, res), F(_F)
#ifdef TRACE_USER
            , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__);

    // TBD!
    ct = nullptr;
}

MEDDLY::user_unary_operation::~user_unary_operation()
{
    ct->markForDestroy();
}


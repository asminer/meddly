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

#include "ct_entry_result.h"

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

MEDDLY::ct_entry_result::ct_entry_result()
{
    build = 0;
    data = 0;
    etype = 0;
}

MEDDLY::ct_entry_result::~ct_entry_result()
{
    delete[] build;
}

void MEDDLY::ct_entry_result::initialize(const ct_entry_type* et)
{
    ASSERT(__FILE__, __LINE__, et);
    etype = et;
    const unsigned slots = etype->getResultSize();
    ASSERT(__FILE__, __LINE__, 0==build);
    build = new ct_entry_item[slots];
}

#endif

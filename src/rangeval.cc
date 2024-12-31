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

#include "rangeval.h"

MEDDLY::rangeval::rangeval(bool v, range_type rt)
{
    the_type = rt;
    s_value = range_special::NORMAL;

    if (isReal()) {
        d_value = v;
    } else {
        l_value = v;
    }
}

MEDDLY::rangeval::rangeval(long v, range_type rt)
{
    the_type = rt;
    s_value = range_special::NORMAL;

    if (isReal()) {
        d_value = v;
    } else {
        l_value = v;
    }
}

MEDDLY::rangeval::rangeval(double v, range_type rt)
{
    the_type = rt;
    s_value = range_special::NORMAL;

    if (isReal()) {
        d_value = v;
    } else {
        l_value = int(v);
    }
}

MEDDLY::rangeval::rangeval(range_special v, range_type rt)
{
    the_type = rt;
    s_value = v;
    if (range_special::NORMAL == v) {
        // set value to 0 as default
        if (isReal()) {
            d_value = 0.0;
        } else {
            l_value = 0;
        }
    }
}


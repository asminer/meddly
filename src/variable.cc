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
#include "variable.h"
#include "error.h"

// ******************************************************************
// *                        variable methods                        *
// ******************************************************************

MEDDLY::variable::variable(int b, char* n)
{
    name = n;
    is_extensible = false;

    if (b < 0) {
        is_extensible = true;
        b = -b;
    }
    un_bound = b;
    pr_bound = b;
}

MEDDLY::variable::~variable()
{
#ifdef DEBUG_CLEANUP
    printf("destroying variable %s\n", name);
#endif
    delete[] name;
}

void MEDDLY::variable::setName(char *n)
{
    delete[] name;
    name = n;
}

void MEDDLY::variable::removeFromList(const domain* d)
{
    const unsigned size = domlist.size();
    if (size<1) return;
    const unsigned last = size-1;

    unsigned find;
    for (find=0; find<size; find++) {
        if (d == domlist[find]) break;
    }

    if (find != last) {
        SWAP(domlist[find], domlist[last]);
    }
    domlist.pop_back();
}

void MEDDLY::variable::enlargeBound(bool prime, int b)
{
    if (b < 1) {
        is_extensible = true;
        b = -b;
    } else if (is_extensible && b > 0) {
        // changing from extensible to non-extensible
        is_extensible = false;
    }

    if (prime) {
        // set prime bound
        if (b > pr_bound) pr_bound = b;
    } else {
        // set prime and unprime bound
        if (b > un_bound) un_bound = b;
        if (b > pr_bound) pr_bound = b;
    }
}

void MEDDLY::variable::shrinkBound(int b, bool force)
{
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}


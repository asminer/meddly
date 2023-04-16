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

    domlist = 0;
    dl_alloc = 0;
    dl_used = 0;
}

MEDDLY::variable::~variable()
{
#ifdef DEBUG_CLEANUP
    printf("destroying variable %s\n", name);
#endif
    delete[] name;
    free(domlist);
}

void MEDDLY::variable::setName(char *n)
{
    delete[] name;
    name = n;
}

void MEDDLY::variable::addToList(domain* d)
{
    if (dl_used >= dl_alloc) {
        int ns = dl_alloc+8;
        domain** dl = (domain**) realloc(domlist, ns * sizeof(void*));
        if (0==dl) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        dl_alloc = ns;
        domlist = dl;
    }
    domlist[dl_used] = d;
    dl_used++;
}

void MEDDLY::variable::removeFromList(const domain* d)
{
    int find;
    for (find=0; find<dl_used; find++) {
        if (d == domlist[find]) break;
    }
    if (find >= dl_used) return;  // not found; should we throw something?
    domlist[find] = domlist[dl_used-1];
    dl_used--;
    // if that was the last domain...
    if (0==dl_used) delete this;
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


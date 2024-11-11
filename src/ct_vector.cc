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

#include "ct_vector.h"
#include "ct_generics.h"

// ******************************************************************
// *                                                                *
// *                         ct_item methods                        *
// *                                                                *
// ******************************************************************

MEDDLY::ct_item::ct_item()
{
    type = ct_typeID::ERROR;
    next = nullptr;
}

void MEDDLY::ct_item::show(output &s) const
{
    switch (type) {
        case ct_typeID::NODE:
            s.put(the_node);
            s.put(" N");
            return;

        case ct_typeID::INTEGER:
            s.put(the_int);
            s.put(" I");
            return;

        case ct_typeID::LONG:
            s.put(the_long);
            s.put(" L");
            return;

        case ct_typeID::FLOAT:
            s.put(the_float);
            s.put(" F");
            return;

        case ct_typeID::DOUBLE:
            s.put(the_double);
            s.put(" D");
            return;

        case ct_typeID::GENERIC:
            if (the_generic) {
                the_generic->show(s);
            } else {
                s.put("null G");
            }
            return;

        default:
            MEDDLY_DCASSERT(false);
    }

}

// ******************************************************************
// *                                                                *
// *                        ct_vector methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::ct_vector::ct_vector(unsigned sz) : _size(sz)
{
    data = useArray(_size);
#ifdef DEVELOPMENT_CODE
    hasHashVal = false;
#endif

    my_entry = 0;
}

MEDDLY::ct_vector::~ct_vector()
{
    recycleArray(data, _size);
}

void MEDDLY::ct_vector::show(output &s) const
{
    s.put("[");
    for (unsigned i=0; i<_size; i++) {
        if (i) s.put(", ");
        data[i].show(s);
    }
    s.put("]");
}

// ******************************************************************

MEDDLY::ct_item* MEDDLY::ct_vector::lists[16];

void MEDDLY::ct_vector::initStatics()
{
    for (unsigned u=0; u<16; u++) {
        lists[u] = nullptr;
    }
}

void MEDDLY::ct_vector::doneStatics()
{
    for (unsigned u=0; u<16; u++) {
        while (lists[u]) {
            ct_item* nxt = lists[u][0].getNext();
            delete[] lists[u];
            lists[u] = nxt;
        }
    }
}

MEDDLY::ct_item* MEDDLY::ct_vector::useArray(unsigned sz)
{
    if (0==sz) return nullptr;
    if (sz<=16) {
        ct_item* curr = lists[sz-1];
        if (curr) {
            lists[sz-1] = curr->getNext();
            return curr;
        }
    }
    return new ct_item[sz];
}

void MEDDLY::ct_vector::recycleArray(ct_item* v, unsigned sz)
{
    if (!v) return;
    MEDDLY_DCASSERT(sz);
    if (sz<=16) {
        v->setNext(lists[sz-1]);
        lists[sz-1] = v;
    } else {
        delete[] v;
    }
}


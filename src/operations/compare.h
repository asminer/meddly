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

#ifndef MEDDLY_COMPARE_H
#define MEDDLY_COMPARE_H

#include "../oper.h"

namespace MEDDLY {
    class binary_operation;

    //
    // Builder for 'equal' operations.
    //
    binary_operation* EQUAL(forest* a, forest* b, forest* c);
    void EQUAL_init();
    void EQUAL_done();

    //
    // The 'not equal' operation builder.
    //
    binary_operation* NOT_EQUAL(forest* a, forest* b, forest* c);
    void NOT_EQUAL_init();
    void NOT_EQUAL_done();

    //
    // The 'greater than' operation builder.
    //
    binary_operation* GREATER_THAN(forest* a, forest* b, forest* c);
    void GREATER_THAN_init();
    void GREATER_THAN_done();

    //
    // The 'greater or equal' operation builder.
    //
    binary_operation* GREATER_THAN_EQUAL(forest* a, forest* b, forest* c);
    void GREATER_THAN_EQUAL_init();
    void GREATER_THAN_EQUAL_done();

    //
    // The 'less than' operation builder.
    //
    binary_operation* LESS_THAN(forest* a, forest* b, forest* c);
    void LESS_THAN_init();
    void LESS_THAN_done();

    //
    // The 'less or equal' operation builder.
    //
    binary_operation* LESS_THAN_EQUAL(forest* a, forest* b, forest* c);
    void LESS_THAN_EQUAL_init();
    void LESS_THAN_EQUAL_done();
};

#endif


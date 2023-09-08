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

#ifndef MEDDLY_OPERATORS_H
#define MEDDLY_OPERATORS_H

#include "dd_edge.h"

/*
 *  Overloading operators for dd_edges.
 *
 *  All of these assume that the edges are "top-level" nodes
 *  (matters only for quasi-reduced).
 *
 *  For binary operators +,-,*,/,&&,||
 *      both operand edges must be attached to the same forest
 *      the result edge will also belong to that forest
 *
 *  For update operators +=, etc
 *      the operands do not need to belong to the same forest.
 *
 *  These simply call apply on the appropriate built-in operation;
 *  see ops_builtin.h
 */

namespace MEDDLY {

    dd_edge operator+(const dd_edge &e1, const dd_edge &e2);
    dd_edge operator-(const dd_edge &e1, const dd_edge &e2);
    dd_edge operator*(const dd_edge &e1, const dd_edge &e2);
    dd_edge operator/(const dd_edge &e1, const dd_edge &e2);

    dd_edge operator&(const dd_edge &e1, const dd_edge &e2);
    dd_edge operator|(const dd_edge &e1, const dd_edge &e2);
    dd_edge operator!(const dd_edge &e);

    dd_edge operator+=(dd_edge &e1, const dd_edge &e2);
    dd_edge operator-=(dd_edge &e1, const dd_edge &e2);
    dd_edge operator*=(dd_edge &e1, const dd_edge &e2);
    dd_edge operator/=(dd_edge &e1, const dd_edge &e2);

    dd_edge operator&=(dd_edge &e1, const dd_edge &e2);
    dd_edge operator|=(dd_edge &e1, const dd_edge &e2);

};

#endif // include guard

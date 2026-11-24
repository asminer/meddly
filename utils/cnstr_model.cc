/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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

#include "utils/cnstr_model.h"

// ******************************************************************

inline void show_comparison(MEDDLY::output &out, cnstr_model::comparison cmp)
{
    switch (cmp) {
        case cnstr_model::comparison::EQ:   return out.put("==");
        case cnstr_model::comparison::NE:   return out.put("!=");
        case cnstr_model::comparison::LT:   return out.put("< ");
        case cnstr_model::comparison::LE:   return out.put("<=");
        case cnstr_model::comparison::GT:   return out.put("> ");
        case cnstr_model::comparison::GE:   return out.put(">=");
    }
}

// ******************************************************************
// *                                                                *
// *                cnstr_model::constraint  methods                *
// *                                                                *
// ******************************************************************

cnstr_model::constraint::constraint(unsigned L) : depends(L)
{
    for (unsigned i=0; i<L; i++) {
        depends[i] = false;
    }
}

// ******************************************************************

cnstr_model::constraint::~constraint()
{
}

// ******************************************************************

void cnstr_model::constraint::show_deps(MEDDLY::output &out) const
{
    for (unsigned i=0; i<depends.size(); i++) {
        if (depends[i]) out.put('X');
        else            out.put('.');
    }
}

// ******************************************************************

bool cnstr_model::constraint::deps_less_than(const constraint* C) const
{
    if (!C) return true;
    MEDDLY_DCASSERT(depends.size() == C->depends.size());
    for (unsigned i=depends.size(); i; ) {
        --i;
        if (depends[i]) {
            if (C->depends[i]) continue;
            return false;
        }
        if (C->depends[i]) {
            return true;
        }
    }
    // equal
    return false;
}

// ******************************************************************
// *                                                                *
// *             cnstr_model::compare_two_vars  methods             *
// *                                                                *
// ******************************************************************

cnstr_model::compare_two_vars::compare_two_vars(unsigned L, unsigned _v1,
        comparison c, unsigned _v2, int d) : constraint(L)
{
    v1 = _v1;
    v2 = _v2;
    delta = d;
    cmp = c;

    add_dependency(v1);
    add_dependency(v2);
}

// ******************************************************************

void cnstr_model::compare_two_vars::build_edge(MEDDLY::dd_edge &r)
{
    using namespace MEDDLY;
    forest* mdd = r.getForest();

    // TBD
}

// ******************************************************************

void cnstr_model::compare_two_vars::show(MEDDLY::output &out) const
{
    out << "x_" << v1;
    show_comparison(out, cmp);
    out << "x_" << v2;
    if (delta<0) {
        out << delta;
    }
    if (delta>0) {
        out << "+" << delta;
    }
}


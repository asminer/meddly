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

#ifndef CNSTR_MODEL_H
#define CNSTR_MODEL_H

#include "meddly/meddly.h"

/**
    Class for a "model" to help coordinate the various
    constraint problem examples.
 */
class cnstr_model {
    public:
        enum class comparison {
            EQ,
            NE,
            LT,
            LE,
            GT,
            GE
        };
    public:
        class constraint {
            protected:
                std::vector <bool> depends;
            public:
                /**
                    Base class constructor.
                        @param  L   Total number of variables;
                                    needed to initialize depends array.
                 */
                constraint(unsigned L);
                virtual ~constraint();

                /**
                    Build the function for this constraint,
                    that returns true iff the constraint holds.
                */
                virtual void build_edge(MEDDLY::dd_edge &r) = 0;

                /**
                    Show the constraint.
                */
                virtual void show(MEDDLY::output &out) const = 0;

                /**
                    Show the dependencies, as a vector.
                */
                void show_deps(MEDDLY::output &out) const;

                /**
                    Return true iff this < C,
                    using dependencies to determine the order.
                 */
                bool deps_less_than(const constraint* C) const;

                inline void add_dependency(unsigned k) {
                    MEDDLY_DCASSERT(k>0);
                    depends[k-1] = true;
                }

        }; // class constraint
    public:
        class compare_two_vars : public constraint {
                unsigned v1;
                unsigned v2;
                int delta;
                comparison cmp;
            public:
                /**
                    The constraint is
                        v1 <c> v2 + d
                    where <c> is a comparison operation.
                 */
                compare_two_vars(unsigned L, unsigned v1,
                        comparison c, unsigned v2, int d);

                inline bool eval(unsigned val1, unsigned val2) const {
                    int i1 = int(val1);
                    int i2 = int(val2);

                    switch (cmp) {
                        case comparison::EQ:    return i1 == (i2 + delta);
                        case comparison::NE:    return i1 != (i2 + delta);
                        case comparison::LT:    return i1 <  (i2 + delta);
                        case comparison::LE:    return i1 <= (i2 + delta);
                        case comparison::GT:    return i1 >  (i2 + delta);
                        case comparison::GE:    return i1 >= (i2 + delta);
                        default:                throw "bad comparison";
                    }
                }

                virtual void build_edge(MEDDLY::dd_edge &r);
                virtual void show(MEDDLY::output &out) const;
        };
};

#endif


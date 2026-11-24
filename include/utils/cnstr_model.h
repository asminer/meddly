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
    public:
        cnstr_model();
        virtual ~cnstr_model();

        /**
            Process command line argument i,
            including any parameters that it takes.
                @param  i       Argument to process. Will be updated.
                @param  argv    Entire argument list.

                @return true,   if we processed the argument
                        false,  if we don't know how to process the argument.
        */
        virtual bool process_arg(unsigned &i, const char** argv);

        /**
            Show all command line switches that we know about.
                @param  out     Stream to write to
                @param  arg0    Executable name to show usage,
                                or nullptr to not show usage.
        */
        virtual void show_usage(std::ostream &out, const char* arg0 = nullptr);

        /**
            Build the model.
            This will be called after we have processed all
            command line arguments with process_arg.

            This should, among other things,
                build the domain and call set_domain(),
                build the constraints and call add_constraint() on each.

            @param  out     Stream to write to
            @param  verb    Verbosity level (0 for quiet).

            @return true on success, false on error
                    (assumed to be bad parameters, so
                    we will then display usage info).
        */
        virtual bool build_model(MEDDLY::output &out) = 0;

        // TBD, more here?

        /**
            Tear down the model.
            Allows us to then cal build_model() again after
            adjusting any model instance parameters.
            Children classes may override this, but should call
            the parent's method also.
        */
        virtual void destroy_model();

    protected:
        /// Options helper: section name
        void show_section(std::ostream &out, const char* secname) const;

        /// Options helper: synopsis
        void show_synopsis(std::ostream &out, const char* arg0,
                const char* usage) const;

        /// Options helper: show program description
        void show_description(std::ostream &out, const char* sec, const char* para) const;

        /// Options helper: show a switch and description
        /// (description lines are automatically split).
        /// switches should be 10 chars or fewer.
        void show_switch(std::ostream &out, const char* sw, const char* desc)
            const;

        inline void set_domain(MEDDLY::domain* d) {
            MEDDLY_DCASSERT(d);
            MEDDLY_DCASSERT(!D);
            D = d;
        }

        inline MEDDLY::domain* get_domain() const {
            return D;
        }

        inline void add_constraint(constraint* c) {
            if (c) constraints.push_back(c);
        }

        /**
            Get the theoretical number of solutions.
                @param  card    Store the result here
                @return true    IFF we can determine this.
            The default behavior is to simply return false.
        */
        virtual bool known_cardinality(double &card) const;

#ifdef HAVE_LIBGMP
        virtual bool known_cardinality(mpz_t card) const;
#endif

        void show_cardinality(MEDDLY::output &out, const MEDDLY::dd_edge &e,
                bool with_timing=true);

        inline static const char* basename(const char* path) {
            const char* b = path;
            for (; *path; path++) {
                if ('/' == *path) b = path+1;
            }
            return b;
        }
        inline static const char* copy_next_word(const char* line,
                char* w, unsigned maxw)
        {
            maxw--;
            unsigned i = 0;
            bool leading = true;
            for (;; ++line) {
                if (*line == 0) {
                    w[i] = 0;
                    break;
                }
                if (' ' == *line) {
                    if (leading) continue;
                    w[i] = 0;
                    break;
                }
                leading = false;
                w[i] = *line;
                i++;
                if (i >= maxw) {
                    w[i] = 0;
                    break;
                }
            }
            return w[0] ? line : nullptr;
        }

        static void printStats(MEDDLY::output &out, const MEDDLY::forest* F);

    protected:  // set defaults for command-line switches
        inline char get_solve_method() const {
            return howsolve;
        }
        inline void set_solve_method(char hs) {
            howsolve = hs;
        }
        inline char get_garbage() const {
            return garbage;
        }
        inline void set_garbage(char g) {
            garbage = g;
        }

        inline void default_approx() {
            approx = true;
        }
        inline void default_exact() {
            approx = false;
        }
        inline bool approx_cardinality() const {
            return approx;
        }


    private:
        std::vector <constraint*> constraints;    // model constraints

        MEDDLY::domain* D;      // domain for the forests
        MEDDLY::forest* Fset;   // boolean forest for constraints + solutions

    private:    // various command-line switches

        // how to solve constraints
        char howsolve;
            // 'a': accumulate in depend order
            // 'f': fold (accumulate in pairs)
            // 'r': accumulate in reverse depend order
            // '<': smallest size (#nodes) first
            // '>': largest size (#nodes) first

        // garbage collection
        char garbage;
            // 'm': mark and sweep
            // 'o': optimistic ref counts
            // 'p': pessimistic ref counts

        // number of solutions to display
        char displaysols;
            // '0': none
            // '1': pick one and show it
            // 'a': show all solutions

        // approximate cardinalities
        bool approx;

        // verbosity level
        int verb;

        // output file: the solutions (or null for none)
        const char* file_set;

};

#endif


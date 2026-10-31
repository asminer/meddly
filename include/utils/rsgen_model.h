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

#ifndef RSGEN_MODEL_H
#define RSGEN_MODEL_H

#include "meddly/meddly.h"

/**
    Class for a "model" to help coordinate the various
    reachability set / distance function generation examples.
 */
class rsgen_model {
    public:
        struct level_delta {
                unsigned level;
                int delta;
            public:
                void set(unsigned lvl, int d) {
                    MEDDLY_DCASSERT(lvl > 0);
                    level = lvl;
                    delta = d;
                }
        };
    public:
        class event {
                const char* name;
            public:
                event(const char* name=nullptr);
                virtual ~event();

                inline const char* get_name() const { return name; }

                /**
                    Explicitly fire an event.
                        @param  curr    Current state
                        @param  next    Next state

                        @return true,   if the event can be fired in the
                                        current state
                                false,  otherwise
                */
                virtual bool fire(const MEDDLY::minterm &curr,
                        MEDDLY::minterm &) = 0;

                /**
                    Build the next state function for this event.
                        @param  mtmxd   Forest for the relation
                        @param  R       Edge to store the result
                        @return true,   iff we could build the relation.
                */
                virtual bool build_nsf(MEDDLY::forest* mtmxd,
                    MEDDLY::dd_edge &R) = 0;

                // something for implicit relations here?
        };
    public:
        class simple_event : public event {
                std::vector <level_delta> updates;
            public:
                simple_event(const char* name,
                        const std::vector <level_delta> &_upd);

                virtual bool fire(const MEDDLY::minterm &curr,
                        MEDDLY::minterm &next);

                virtual bool build_nsf(MEDDLY::forest* mtmxd,
                        MEDDLY::dd_edge &R);
        };
    public:
        rsgen_model();
        virtual ~rsgen_model();

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
                build the events and call add_event() on each.

            @return true on success, false on error
                    (assumed to be bad parameters, so
                    we will then display usage info).
        */
        virtual bool build_model() = 0;

        /// Build the initial state, as a minterm.
        virtual void build_initial(MEDDLY::minterm &init) = 0;

        /**
            Generate.
            Will generate the state space or distance function,
            using the appropriate methods, as dictated by the
            command line switches.
        */
        void generate();

        /// Build a monolithic relation.
        void build_monolithic(MEDDLY::dd_edge &R);

    protected:
        inline void set_domain(MEDDLY::domain* d) {
            MEDDLY_DCASSERT(d);
            MEDDLY_DCASSERT(!D);
            D = d;
        }

        inline void add_event(event* e) {
            events.push_back(e);
        }

        inline static const char* basename(const char* path) {
            const char* b = path;
            for (; *path; path++) {
                if ('/' == *path) b = path+1;
            }
            return b;
        }

        static void printStats(const char* who, const MEDDLY::forest* F);

    private:
        std::vector <event*> events;    // model events

        MEDDLY::domain* D;      // domain for the forests

        MEDDLY::forest* Fset;   // forest for reachset or distances
        MEDDLY::forest* Fmxd;   // forest for relation, or null
        MEDDLY::forest* Fmtmxd; // auxiliary forest for relation, or null

        // generation method
        char method;

        // approximate cardinalities
        bool approx;

        // verbosity level
        int verb;
};

#endif

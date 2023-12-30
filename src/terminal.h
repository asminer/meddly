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

#ifndef MEDDLY_TERMINAL_H
#define MEDDLY_TERMINAL_H

#include "defines.h"
#include "error.h"
// #include <cstring>  // for memcpy

namespace MEDDLY {
    class terminal;

    /// Type of terminal value.
    enum class terminal_type {
        /// Null type
        NOTHING,
        /// Nothing; for EV
        OMEGA,
        /// For MT boolean
        BOOLEAN,
        /// For MT integer
        INTEGER,
        /// For MT real
        REAL
    };

    class input;
    class output;
};


/**
    Unified terminal object.
*/
class MEDDLY::terminal {

    public:
        /// Don't initialize
        terminal() {
            mytype = terminal_type::NOTHING;
        }

        /// Initialize from a handle.
        terminal(terminal_type t, node_handle h) {
            setFromHandle(t, h);
        }

        /// Initialize from a boolean label
        terminal(bool v) {
            setBoolean(v);
        }
        /// Initialize from an integer label
        terminal(int v) {
            setInteger(v);
        }
        /// Initialize from an integer label
        terminal(long v) {
            setInteger(v);
        }
        /// Initialize from a real label
        terminal(float v) {
            setReal(v);
        }

        /// Initialize from a boolean label, but convert types if needed
        terminal(bool v, terminal_type t) {
            setFromValue(t, v);
        }
        /// Initialize from an integer label, but convert types if needed
        terminal(int v, terminal_type t) {
            setFromValue(t, v);
        }
        /// Initialize from an integer label, but convert types if needed
        terminal(long v, terminal_type t) {
            setFromValue(t, v);
        }
        /// Initialize from a real label, but convert types if needed
        terminal(float v, terminal_type t) {
            setFromValue(t, v);
        }

        //
        // Getters for the type
        //

        inline bool isOmega() const {
            return (terminal_type::OMEGA == mytype);
        }
        inline bool isBoolean() const {
            return (terminal_type::BOOLEAN == mytype);
        }
        inline bool isInteger() const {
            return (terminal_type::INTEGER == mytype);
        }
        inline bool isReal() const {
            return (terminal_type::REAL == mytype);
        }
        inline bool hasType(terminal_type t) const {
            return (t == mytype);
        }

        //
        // Getters for the label
        //

        inline node_handle getOmega() const {
            MEDDLY_DCASSERT(isOmega());
            return t_omega;
        }
        inline bool getBoolean() const {
            MEDDLY_DCASSERT(isBoolean());
            return t_boolean;
        }
        inline int getInteger() const {
            MEDDLY_DCASSERT(isInteger());
            return t_integer;
        }
        inline float getReal() const {
            MEDDLY_DCASSERT(isReal());
            return t_real;
        }

        //
        // Get the label into any type
        // Allows for type conversions
        //
        template <typename T>
        inline void getValue(T &v) const {
            switch (mytype) {
                case terminal_type::OMEGA:
                        v = T(t_omega);
                        return;

                case terminal_type::BOOLEAN:
                        v = T(t_boolean);
                        return;

                case terminal_type::INTEGER:
                        v = T(t_integer);
                        return;

                case terminal_type::REAL:
                        v = T(t_real);
                        return;

                default:
                        MEDDLY_DCASSERT(false);
            }

        }

        //
        // Get the node handle for this terminal
        //
        inline node_handle getHandle() const {
            switch (mytype) {
                case terminal_type::OMEGA:
                        return t_omega;

                case terminal_type::BOOLEAN:
                        return t_boolean ? -1 : 0;

                case terminal_type::INTEGER:
                        if (t_integer) {
                            return t_integer | -2147483648; // set sign bit
                        } else {
                            return 0;
                        }

                case terminal_type::REAL:
                        if (t_real) {
                            union {
                                node_handle h;
                                float f;
                            } x;
                            x.f = t_real;
                            // const unsigned x = * (unsigned*) (&t_real);
                            // strip the lsb in fraction, and add sign bit
                            return (x.h>>1) | -2147483648;
                            /*
                            unsigned x;
                            memcpy(&x, &t_real, sizeof(unsigned));
                            return node_handle(x>>1) | -2147483648;
                            */
                        } else {
                            return 0;
                        }

                default:
                        MEDDLY_DCASSERT(false);
                        return 0;
            }
        }

        //
        // Setters for the label
        //

        inline void setOmega(node_handle h=-1) {
            mytype = terminal_type::OMEGA;
            t_omega = h;
            if (h > 0) {
                throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
            }
        }
        inline void setBoolean(bool v) {
            mytype = terminal_type::BOOLEAN;
            t_boolean = v;
        }
        inline void setInteger(int v) {
            mytype = terminal_type::INTEGER;
            if (v < -1073741824 || v > 1073741823) {
                // Can't fit in 31 bits (signed)
                throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
            }
            t_integer = v;
        }
        inline void setReal(float v) {
            mytype = terminal_type::REAL;
            // TBD: should we mask off the lsb here?
            t_real = v;
        }

        //
        // Set the terminal from any value
        // Allows for type conversions
        //
        template <typename T>
        inline void setFromValue(terminal_type t, T v) {
            switch (mytype=t) {
                case terminal_type::OMEGA:
                        t_omega = v;
                        return;

                case terminal_type::BOOLEAN:
                        t_boolean = v;
                        return;

                case terminal_type::INTEGER:
                        t_integer = int(v);
                        return;

                case terminal_type::REAL:
                        t_real = float(v);
                        return;

                default:
                        MEDDLY_DCASSERT(false);
            }

        }

        //
        // Set the terminal from a node handle
        //
        inline void setFromHandle(terminal_type t, node_handle h) {
            switch (mytype=t) {
                case terminal_type::OMEGA:
                        t_omega = h;
                        return;

                case terminal_type::BOOLEAN:
                        if (h<-1 || h>0) {
                            throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
                        }
                        t_boolean = (h!=0);
                        return;

                case terminal_type::INTEGER:
                        // << 1 kills the sign bit
                        // >> 1 puts us back, and extends the (new) sign bit
                        // Works correctly also for 0
                        t_integer = (h << 1) >> 1;
                        return;

                case terminal_type::REAL:
                        union {
                            node_handle h;
                            float f;
                        } x;
                        // Strip sign bit
                        x.h = (h << 1);
                        t_real = x.f;
                        return;
                        /*
                        h <<= 1;
                        memcpy(&t_real, &h, sizeof(float));
                        */

                default:
                        MEDDLY_DCASSERT(false);
            }
        }

        //
        // Comparisons, so we can have sets of terminals
        //
        inline bool operator<(const terminal &t) const {
            switch (mytype) {
                case terminal_type::OMEGA:
                        return t_omega < t.getOmega();

                case terminal_type::BOOLEAN:
                        return t_boolean < t.getBoolean();

                case terminal_type::INTEGER:
                        return t_integer < t.getInteger();

                case terminal_type::REAL:
                        return t_real < t.getReal();

                default:
                        MEDDLY_DCASSERT(false);
            }
            return false;
        }

        //
        // File I/O
        //

        /**
            Read a terminal node from an input file (stream).
            This will set the type appropriately.
        */
        void read(input &s);

        /**
            Write a terminal node to an output file (stream).
            Formatted so that we can read the value back
            using read().
        */
        void write(output &s) const;

    private:
        union {
            node_handle     t_omega;        // for special values
            bool            t_boolean;
            int             t_integer;
            float           t_real;
        };

        terminal_type mytype;
};


#endif // include guard

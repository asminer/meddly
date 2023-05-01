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

#ifndef MEDDLY_OPNAME_H
#define MEDDLY_OPNAME_H

#include "defines.h"

namespace MEDDLY {
    class opname;
    class unary_opname;
    class binary_opname;
    class specialized_opname;

    class unary_operation;
    class binary_operation;
    class specialized_operation;

    class initializer_list;
    class expert_forest;

    void initialize(initializer_list *);
    void cleanup();

    class opname_init;      // hidden in opname.cc

    // ******************************************************************
    // *                     Named unary operations                     *
    // ******************************************************************

    /// Copy a function across forests, with the same domain.
    const unary_opname* COPY();

    /// Unary operation.  Return the number of variable assignments
    /// so that the function evaluates to non-zero.
    const unary_opname* CARDINALITY();

    /// For BOOLEAN forests, flip the return values.
    const unary_opname* COMPLEMENT();

    /// Find the largest value returned by the function.
    const unary_opname* MAX_RANGE();

    /// Find the smallest value returned by the function.
    const unary_opname* MIN_RANGE();

    /// Convert MDD to EV+MDD index set.  A special case of COPY, really.
    const unary_opname* CONVERT_TO_INDEX_SET();

    /// Extract cycles (EV+MDD) from transitive closure (EV+MxD)
    const unary_opname* CYCLE();

    /// Randomly select one state from a set of states
    const unary_opname* SELECT();
};

// ******************************************************************
// *                                                                *
// *                          opname class                          *
// *                                                                *
// ******************************************************************

/**
    Class (hierarchy) for names of operations.
    These are used as factories, to build specific operations
    (tied to forests) to be applied to decision diagrams
*/
class MEDDLY::opname {
    public:
        inline static const unary_opname* COPY() {
            return _COPY;
        }
        inline static const unary_opname* CARDINALITY() {
            return _CARD;
        }
        inline static const unary_opname* COMPLEMENT() {
            return _COMPL;
        }
        inline static const unary_opname* MAX_RANGE() {
            return _MAXRANGE;
        }
        inline static const unary_opname* MIN_RANGE() {
            return _MINRANGE;
        }
        inline static const unary_opname* CONVERT_TO_INDEX_SET() {
            return _MDD2INDEX;
        }
        inline static const unary_opname* CYCLE() {
            return _CYCLE;
        }
        inline static const unary_opname* SELECT() {
            return _SELECT;
        }

    public:
        opname(const char* n);
        virtual ~opname();

        inline size_t getIndex() const { return index; }
        inline const char* getName() const { return name; }

        inline static size_t nextIndex() { return next_index; }

        static initializer_list* makeInitializer(initializer_list* prev);

    private:
        const char* name;
        size_t index;
        static size_t next_index;

    private:
        static unary_opname* _COPY;
        static unary_opname* _CARD;
        static unary_opname* _COMPL;
        static unary_opname* _MAXRANGE;
        static unary_opname* _MINRANGE;
        static unary_opname* _MDD2INDEX;
        static unary_opname* _CYCLE;
        static unary_opname* _SELECT;

    friend class MEDDLY::opname_init;

        friend void MEDDLY::initialize(initializer_list *);
        friend void MEDDLY::cleanup();
};

// ******************************************************************
// *                                                                *
// *                       unary_opname class                       *
// *                                                                *
// ******************************************************************

/// Unary operation names.
class MEDDLY::unary_opname : public opname {
    public:
        unary_opname(const char* n);
        virtual ~unary_opname();

        virtual unary_operation*
            buildOperation(expert_forest* arg, expert_forest* res) const;

        virtual unary_operation*
            buildOperation(expert_forest* arg, opnd_type res) const;
};

// ******************************************************************
// *                                                                *
// *                      binary_opname  class                      *
// *                                                                *
// ******************************************************************

/// Binary operation names.
class MEDDLY::binary_opname : public opname {
    public:
        binary_opname(const char* n);
        virtual ~binary_opname();

        virtual binary_operation* buildOperation(expert_forest* arg1,
            expert_forest* arg2, expert_forest* res) const = 0;
};


// ******************************************************************
// *                                                                *
// *                    specialized_opname class                    *
// *                                                                *
// ******************************************************************

/// Specialized operation names.
class MEDDLY::specialized_opname : public opname {
    public:
        /**
            Abstract base class for arguments to buildOperation().
            Derived operation names must provide derived classes for arguments.
        */
        class arguments {
            public:
                arguments();
                virtual ~arguments();

                /**
                    Specify if arguments should be destroyed or not.
                    If yes (the default), the operation will destroy
                    the arguments once they are no longer needed.
                */
                inline void setAutoDestroy(bool destroy) {
                    destroyWhenDone = destroy;
                }
                inline bool autoDestroy() const {
                    return destroyWhenDone;
                }

            private:
                bool destroyWhenDone;
        };

    public:
        specialized_opname(const char* n);
        virtual ~specialized_opname();

        /** Note - unlike the more general binary and unary ops,
            a specialized operation might be crafted to the specific
            arguments (passed as an abstract class).
            Examples are:
                -   operations that will be called several times with
                    several of the arguments unchanged, and we want
                    to do some preprocessing on the unchanged arguments.

                -   operations with very bizarre and/or user-defined
                    parameters (like on-the-fly saturation).

            @param  a   Arguments.  Will be destroyed when we are finished,
                        if autoDestroy() is set for the arguments.
        */
        virtual specialized_operation* buildOperation(arguments* a) const = 0;
};

#endif // #include guard

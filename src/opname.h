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

    class operation;
    class unary_operation;
    class binary_operation;
    class specialized_operation;

    class ct_object;
    class dd_edge;

    class initializer_list;

    void initialize(initializer_list *);
    void cleanup();

    class opname_init;      // hidden in opname.cc

    /// Argument and result types for apply operations.
    enum class opnd_type {
        FOREST      = 0,
        BOOLEAN     = 1,
        INTEGER     = 2,
        REAL        = 3,
        HUGEINT     = 4,
        FLOATVECT   = 5,
        DOUBLEVECT  = 6
    };

    // ******************************************************************
    // *                          Unary  apply                          *
    // ******************************************************************

    /** Apply a unary operator.
        The operand and the result are not necessarily in the same forest,
        but they must belong to forests that share the same domain.
        This is useful, for instance, for copying a function to a new forest.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result, where \a c = \a op \a a.
    */
    void apply(unary_opname* (*op)(), const dd_edge &a, dd_edge &c);


    /** Apply a unary operator.
        For operators whose result is an integer.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result, where \a c = \a op \a a.
    */
    void apply(unary_opname* (*op)(), const dd_edge &a, long &c);

    /** Apply a unary operator.
        For operators whose result is a real.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result, where \a c = \a op \a a.
    */
    void apply(unary_opname* (*op)(), const dd_edge &a, double &c);

    /** Apply a unary operator.
        For operators whose result is a anything else.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result, where \a c = \a op \a a.
    */
    void apply(unary_opname* (*op)(), const dd_edge &a, opnd_type cr,
        ct_object &c);

#ifdef __GMP_H__
    /** Apply a unary operator.
        For operators whose result is an arbitrary-precision integer
        (as supplied by the GNU MP library).
            @param  op  Operator handle.
            @param  a   Operand.
            @param  c   Input: an initialized MP integer.
                        Output parameter: the result, where \a c = \a op \a a.
    */
    void apply(unary_opname* (*op)(), const dd_edge &a, mpz_t &c);
#endif

    // ******************************************************************
    // *                          Binary apply                          *
    // ******************************************************************

    /** Apply a binary operator.
        \a a, \a b and \a c are not required to be in the same forest,
        but they must have the same domain. The result will be in the
        same forest as \a result. The operator decides the type of forest
        for each \a dd_edge.
        Useful, for example, for constructing comparisons
        where the resulting type is "boolean" but the operators are not,
        e.g., c = f EQUALS g.
            @param  op    Operator handle.
            @param  a     First operand.
            @param  b     Second operand.
            @param  c     Output parameter: the result,
                          where \a c = \a a \a op \a b.
    */
    void apply(binary_opname* (*op)(), const dd_edge &a, const dd_edge &b,
        dd_edge &c);


    // ******************************************************************
    // *                     Named unary operations                     *
    // ******************************************************************

    /// Copy a function across forests, with the same domain.
    unary_opname* COPY();

    /// Unary operation.  Return the number of variable assignments
    /// so that the function evaluates to non-zero.
    unary_opname* CARDINALITY();

    /// For BOOLEAN forests, flip the return values.
    unary_opname* COMPLEMENT();

    /// Find the largest value returned by the function.
    unary_opname* MAX_RANGE();

    /// Find the smallest value returned by the function.
    unary_opname* MIN_RANGE();

    /// Convert MDD to EV+MDD index set.  A special case of COPY, really.
    unary_opname* CONVERT_TO_INDEX_SET();

    /// Extract cycles (EV+MDD) from transitive closure (EV+MxD)
    unary_opname* CYCLE();

    /// Randomly select one state from a set of states
    unary_opname* SELECT();

    // ******************************************************************
    // *                    Named  binary operations                    *
    // ******************************************************************

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
        inline static unary_opname* COPY() {
            return _COPY;
        }
        inline static unary_opname* CARDINALITY() {
            return _CARD;
        }
        inline static unary_opname* COMPLEMENT() {
            return _COMPL;
        }
        inline static unary_opname* MAX_RANGE() {
            return _MAXRANGE;
        }
        inline static unary_opname* MIN_RANGE() {
            return _MINRANGE;
        }
        inline static unary_opname* CONVERT_TO_INDEX_SET() {
            return _MDD2INDEX;
        }
        inline static unary_opname* CYCLE() {
            return _CYCLE;
        }
        inline static unary_opname* SELECT() {
            return _SELECT;
        }

    public:
        opname(const char* n);
        virtual ~opname();

//        inline size_t getIndex() const { return index; }
        inline const char* getName() const { return name; }

//        inline static size_t nextIndex() { return next_index; }

        static initializer_list* makeInitializer(initializer_list* prev);

    protected:
        void removeFromCache(operation* op);

    protected:
        operation* cache;

    private:
        const char* name;
//        size_t index;   // TBD - needed?
//        static size_t next_index;

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

        unary_operation* getOperation(const dd_edge &arg, const dd_edge &res);
        unary_operation* getOperation(const dd_edge &arg, opnd_type res);

        void removeOperationFromCache(unary_operation* op);

    protected:
        virtual unary_operation*
            buildOperation(const dd_edge &arg, const dd_edge &res) const;

        virtual unary_operation*
            buildOperation(const dd_edge &arg, opnd_type res) const;

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

        binary_operation* getOperation(const dd_edge &arg1,
                const dd_edge &arg2, const dd_edge &res);

        void removeOperationFromCache(binary_operation* op);

    protected:
        virtual binary_operation* buildOperation(const dd_edge &arg1,
            const dd_edge &arg2, const dd_edge &res) const = 0;
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

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
    // *                    Wrapper for GMP integers                    *
    // ******************************************************************

#ifdef __GMP_H__
    ct_object& get_mpz_wrapper();
    void unwrap(const ct_object &, mpz_t &value);
#endif

    // ******************************************************************
    // *                          Unary  apply                          *
    // ******************************************************************

    /** Apply a unary operator.
        The operand and the result are not necessarily in the same forest,
        but they must belong to forests that share the same domain.
        This is useful, for instance, for copying a function to a new forest.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result,
                            where \a c = \a op \a a.
    */
    void apply(unary_opname* (*op)(), const dd_edge &a, dd_edge &c);


    /** Apply a unary operator.
        For operators whose result is an integer.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result,
                            where \a c = \a op \a a.
    */
    void apply(unary_opname* (*op)(), const dd_edge &a, long &c);

    /** Apply a unary operator.
        For operators whose result is a real.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result,
                            where \a c = \a op \a a.
    */
    void apply(unary_opname* (*op)(), const dd_edge &a, double &c);

    /** Apply a unary operator.
        For operators whose result is a anything else.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result,
                            where \a c = \a op \a a.
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
                        Output: the result, where \a c = \a op \a a.
    */
    inline void apply(unary_opname* (*op)(), const dd_edge &a, mpz_t &c) {
        ct_object& x = get_mpz_wrapper();
        apply(op, a, opnd_type::HUGEINT, x);
        unwrap(x, c);
    }
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

    /// Unary operation.  Return the number of variable assignments
    /// so that the function evaluates to non-zero.
    unary_opname* CARDINALITY();

    /// For BOOLEAN forests, flip the return values.
    unary_opname* COMPLEMENT();

    /// Convert MDD to EV+MDD index set.  A special case of COPY, really.
    unary_opname* CONVERT_TO_INDEX_SET();

    /// Copy a function across forests, with the same domain.
    unary_opname* COPY();

    /// Extract cycles (EV+MDD) from transitive closure (EV+MxD)
    unary_opname* CYCLE();

    /// Find the largest value returned by the function.
    unary_opname* MAX_RANGE();

    /// Find the smallest value returned by the function.
    unary_opname* MIN_RANGE();

    /// Randomly select one state from a set of states
    unary_opname* SELECT();

    // ******************************************************************
    // *                    Named  binary operations                    *
    // ******************************************************************

    /// Set union operation for forests with range_type of BOOLEAN
    binary_opname* UNION();

    /// Set intersection operation for forests with range_type of BOOLEAN
    binary_opname* INTERSECTION();

    /// Set difference operation for forests with range_type of BOOLEAN
    binary_opname* DIFFERENCE();

    /// Combine two functions into a single one, where the operands are MDDs
    /// and the result is an MXD.  Specifically, for MDD operands f and g,
    /// produces MXD h where
    ///     h(xn, x'n, ..., x1, x'1) = f(xn, ..., x1) * g(x'n, ..., x'1).
    /// Works for BOOLEAN forests.
    binary_opname* CROSS();

    /// The minimum of two functions, with range_type INTEGER or REAL
    binary_opname* MINIMUM();

    /// The maximum of two functions, with range_type INTEGER or REAL
    binary_opname* MAXIMUM();

    /// Add two functions, with range type INTEGER and REAL
    binary_opname* PLUS();

    /// Subtract two functions, with range type INTEGER and REAL
    binary_opname* MINUS();

    /// Multiply two functions, with range type INTEGER and REAL
    binary_opname* MULTIPLY();

    /// Divide two functions, with range type INTEGER and REAL
    binary_opname* DIVIDE();

    /// Take the remainder of two functions, with range type INTEGER
    binary_opname* MODULO();

    /// Compare for equality, two functions with range type INTEGER or REAL
    binary_opname* EQUAL();

    /// Compare for inequality, two functions with range type INTEGER or REAL
    binary_opname* NOT_EQUAL();

    /// Compare for <, two functions with range type INTEGER or REAL
    binary_opname* LESS_THAN();

    /// Compare for <=, two functions with range type INTEGER or REAL
    binary_opname* LESS_THAN_EQUAL();

    /// Compare for >, two functions with range type INTEGER or REAL
    binary_opname* GREATER_THAN();

    /// Compare for >=, two functions with range type INTEGER or REAL
    binary_opname* GREATER_THAN_EQUAL();

    /** Plus operation used to compute transitive closure and further
        minimum witness. The first operand must be an EV+MxD and the second
        operand must be an EV+MDD. The result is an EV+MxD.
    */
    binary_opname* PRE_PLUS();
    binary_opname* POST_PLUS();

    /** Image operations on a set-of-states with a transition function.
        The first operand must be the set-of-states and the second operand
        must be the transition function. The result is a set-of-states that
        must be stored in the same forest as the first operand.

        Applies to:
        PRE_IMAGE, POST_IMAGE,
        TC_POST_IMAGE,
        REACHABLE_STATES_DFS, REACHABLE_STATES_BFS,
        REVERSE_REACHABLE_DFS, REVERSE_REACHABLE_BFS.
    */
    binary_opname* PRE_IMAGE();
    binary_opname* POST_IMAGE();
    binary_opname* TC_POST_IMAGE();
    binary_opname* REACHABLE_STATES_DFS();
    binary_opname* REACHABLE_STATES_BFS();
    binary_opname* REVERSE_REACHABLE_DFS();
    binary_opname* REVERSE_REACHABLE_BFS();

    /** Vector matrix multiply, where the first argument is vector (MDD),
        the second argument is a matrix (MXD),
        and the result is a vector (MDD).
    */
    binary_opname* VM_MULTIPLY();

    /** Matrix vector multiply, where the first argument is a matrix (MXD),
        the second argument is a vector (MDD),
        and the result is a vector (MDD).
    */
    binary_opname* MV_MULTIPLY();

    /** Matrix multiplication, where the first argument is a matrix (MXD),
        the second argument is a matrix (MXD),
        and the result is a matrix (MXD),
        such that, C[m][n] += A[m][i] * B[i][n], for all m, n and i.
    */
    binary_opname* MM_MULTIPLY();
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
        friend class MEDDLY::opname_init;
    public:
        opname(const char* n);
        virtual ~opname();

        inline const char* getName() const { return name; }

        static initializer_list* makeInitializer(initializer_list* prev);

    protected:
        void removeFromCache(operation* op);

    protected:
        operation* cache;

    private:
        const char* name;

    private:
        static unary_opname* _CARD;
        static unary_opname* _COMPL;
        static unary_opname* _MDD2INDEX;
        static unary_opname* _COPY;
        static unary_opname* _CYCLE;
        static unary_opname* _MAXRANGE;
        static unary_opname* _MINRANGE;
        static unary_opname* _SELECT;
    private:
        static binary_opname* _UNION;
        static binary_opname* _INTERSECT;
        static binary_opname* _DIFFERENCE;
        static binary_opname* _CROSS;
        static binary_opname* _MIN;
        static binary_opname* _MAX;
        static binary_opname* _PLUS;
        static binary_opname* _MINUS;
        static binary_opname* _MULTIPLY;
        static binary_opname* _DIVIDE;
        static binary_opname* _MODULO;
        static binary_opname* _EQ;
        static binary_opname* _NE;
        static binary_opname* _LT;
        static binary_opname* _LE;
        static binary_opname* _GT;
        static binary_opname* _GE;
        static binary_opname* _PRE_PLUS;
        static binary_opname* _POST_PLUS;
        static binary_opname* _PRE_IMAGE;
        static binary_opname* _POST_IMAGE;
        static binary_opname* _TC_POST_IMAGE;
        static binary_opname* _FORWARD_DFS;
        static binary_opname* _FORWARD_BFS;
        static binary_opname* _BACKWARD_DFS;
        static binary_opname* _BACKWARD_BFS;
        static binary_opname* _VM_MULTIPLY;
        static binary_opname* _MV_MULTIPLY;
        static binary_opname* _MM_MULTIPLY;
    public:
        inline static unary_opname* CARDINALITY() {
            return _CARD;
        }
        inline static unary_opname* COMPLEMENT() {
            return _COMPL;
        }
        inline static unary_opname* CONVERT_TO_INDEX_SET() {
            return _MDD2INDEX;
        }
        inline static unary_opname* COPY() {
            return _COPY;
        }
        inline static unary_opname* CYCLE() {
            return _CYCLE;
        }
        inline static unary_opname* MAX_RANGE() {
            return _MAXRANGE;
        }
        inline static unary_opname* MIN_RANGE() {
            return _MINRANGE;
        }
        inline static unary_opname* SELECT() {
            return _SELECT;
        }
    public:
        inline static binary_opname* UNION() {
            return _UNION;
        }
        inline static binary_opname* INTERSECTION() {
            return _INTERSECT;
        }
        inline static binary_opname* DIFFERENCE() {
            return _DIFFERENCE;
        }
        inline static binary_opname* CROSS() {
            return _CROSS;
        }
        inline static binary_opname* MINIMUM() {
            return _MIN;
        }
        inline static binary_opname* MAXIMUM() {
            return _MAX;
        }
        inline static binary_opname* PLUS() {
            return _PLUS;
        }
        inline static binary_opname* MINUS() {
            return _MINUS;
        }
        inline static binary_opname* MULTIPLY() {
            return _MULTIPLY;
        }
        inline static binary_opname* DIVIDE() {
            return _DIVIDE;
        }
        inline static binary_opname* MODULO() {
            return _MODULO;
        }
        inline static binary_opname* EQUAL() {
            return _EQ;
        }
        inline static binary_opname* NOT_EQUAL() {
            return _NE;
        }
        inline static binary_opname* LESS_THAN() {
            return _LT;
        }
        inline static binary_opname* LESS_THAN_EQUAL() {
            return _LE;
        }
        inline static binary_opname* GREATER_THAN() {
            return _GT;
        }
        inline static binary_opname* GREATER_THAN_EQUAL() {
            return _GE;
        }
        inline static binary_opname* PRE_PLUS() {
            return _PRE_PLUS;
        }
        inline static binary_opname* POST_PLUS() {
            return _POST_PLUS;
        }
        inline static binary_opname* PRE_IMAGE() {
            return _PRE_IMAGE;
        }
        inline static binary_opname* POST_IMAGE() {
            return _POST_IMAGE;
        }
        inline static binary_opname* TC_POST_IMAGE() {
            return _TC_POST_IMAGE;
        }
        inline static binary_opname* REACHABLE_STATES_DFS() {
            return _FORWARD_DFS;
        }
        inline static binary_opname* REACHABLE_STATES_BFS() {
            return _FORWARD_BFS;
        }
        inline static binary_opname* REVERSE_REACHABLE_DFS() {
            return _BACKWARD_DFS;
        }
        inline static binary_opname* REVERSE_REACHABLE_BFS() {
            return _BACKWARD_BFS;
        }
        inline static binary_opname* VM_MULTIPLY() {
            return _VM_MULTIPLY;
        }
        inline static binary_opname* MV_MULTIPLY() {
            return _MV_MULTIPLY;
        }
        inline static binary_opname* MM_MULTIPLY() {
            return _MM_MULTIPLY;
        }
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

// ******************************************************************
// *                                                                *
// *                      Old interface  stuff                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    unary_operation* getOperation(unary_opname* (*op)(),
            const dd_edge &a, const dd_edge &c)
    {
        unary_opname* uop = op();
        MEDDLY_DCASSERT(uop);
        return uop->getOperation(a, c);
    }

    unary_operation* getOperation(unary_opname* (*op)(),
            const dd_edge &arg, opnd_type res)
    {
        unary_opname* uop = op();
        MEDDLY_DCASSERT(uop);
        return uop->getOperation(arg, res);
    }

    binary_operation* getOperation(binary_opname* (*op)(),
            const dd_edge &a1, const dd_edge &a2, const dd_edge &r)
    {
        binary_opname* bop = op();
        MEDDLY_DCASSERT(bop);
        return bop->getOperation(a1, a2, r);
    }
};

#endif // #include guard

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

    typedef unary_opname* (*unary_handle)();
    typedef binary_opname* (*binary_handle)();

    class operation;
    class unary_operation;
    class binary_operation;
    class specialized_operation;

    class ct_object;
    class dd_edge;
    class expert_forest;

    class initializer_list;

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
    void apply(unary_handle op, const dd_edge &a, dd_edge &c);


    /** Apply a unary operator.
        For operators whose result is an integer.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result,
                            where \a c = \a op \a a.
    */
    void apply(unary_handle op, const dd_edge &a, long &c);

    /** Apply a unary operator.
        For operators whose result is a real.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result,
                            where \a c = \a op \a a.
    */
    void apply(unary_handle op, const dd_edge &a, double &c);

    /** Apply a unary operator.
        For operators whose result is a anything else.
            @param  op    Operator handle.
            @param  a     Operand.
            @param  c     Output parameter: the result,
                            where \a c = \a op \a a.
    */
    void apply(unary_handle op, const dd_edge &a, opnd_type cr,
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
    inline void apply(unary_handle op, const dd_edge &a, mpz_t &c) {
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
    void apply(binary_handle op, const dd_edge &a, const dd_edge &b,
        dd_edge &c);

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
        opname(const char* n);
        virtual ~opname();

        inline const char* getName() const { return name; }

    protected:
        void removeFromCache(operation* op);

    protected:
        operation* cache;

    private:
        const char* name;
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
            buildOperation(const dd_edge &arg, const dd_edge &res);

        virtual unary_operation*
            buildOperation(const dd_edge &arg, opnd_type res);

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
        // OLD interface
        virtual binary_operation* buildOperation(expert_forest* arg1,
            expert_forest* arg2, expert_forest* res) = 0;
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
        virtual specialized_operation* buildOperation(arguments* a) = 0;
};

// ******************************************************************
// *                                                                *
// *                      Old interface  stuff                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    inline unary_operation* getOperation(unary_handle op,
            const dd_edge &a, const dd_edge &c)
    {
        unary_opname* uop = op();
        MEDDLY_DCASSERT(uop);
        return uop->getOperation(a, c);
    }

    inline unary_operation* getOperation(unary_handle op,
            const dd_edge &arg, opnd_type res)
    {
        unary_opname* uop = op();
        MEDDLY_DCASSERT(uop);
        return uop->getOperation(arg, res);
    }

    inline binary_operation* getOperation(binary_handle op,
            const dd_edge &a1, const dd_edge &a2, const dd_edge &r)
    {
        binary_opname* bop = op();
        MEDDLY_DCASSERT(bop);
        return bop->getOperation(a1, a2, r);
    }
};

#endif // #include guard

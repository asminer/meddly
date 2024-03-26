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

#ifndef MEDDLY_OPER_UNARY_H
#define MEDDLY_OPER_UNARY_H

#include "dd_edge.h"
#include "oper.h"

namespace MEDDLY {
    class ct_object;
    class unary_operation;
    class unary_list;
    class forest;
};


// ******************************************************************
// *                                                                *
// *                     unary_operation  class                     *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a unary operation in a specific forest.
    Specific operations will be derived from this class.
*/
class MEDDLY::unary_operation : public operation {
    public:
        unary_operation(unary_list& owner, unsigned et_slots,
            forest* arg, forest* res);

        unary_operation(unary_list& owner, unsigned et_slots,
            forest* arg, opnd_type res);

    protected:
        virtual ~unary_operation();

    public:
        /**
            Checks forest comatability and then calls computeDDEdge().
        */
        void compute(const dd_edge &arg, dd_edge &res);
        void computeTemp(const dd_edge &arg, dd_edge &res);

        virtual void compute(const dd_edge &arg, long &res);
        virtual void compute(const dd_edge &arg, double &res);
        virtual void compute(const dd_edge &arg, ct_object &c);
        virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);

    protected:
        virtual bool checkForestCompatibility() const;


    protected:
        forest* argF;
        forest* resF;
        opnd_type resultType;

    private:
        unary_list& parent;
        unary_operation* next;

        friend void destroyOperation(unary_operation* &op);
        friend class unary_list;
};

// ******************************************************************
// *                                                                *
// *                        unary_list class                        *
// *                                                                *
// ******************************************************************

/**
    List of unary operations of the same type; used when
    building unary operations for specific forests.
*/
class MEDDLY::unary_list {
        const char* name;
        unary_operation* front;
    public:
        unary_list(const char* n = nullptr);

        inline void setName(const char* n) { name = n; }
        inline const char* getName() const { return name; }

        inline unary_operation* addOperation(unary_operation* uop) {
            if (uop) {
                uop->next = front;
                front = uop;
            }
            return uop;
        }

        inline void removeOperation(unary_operation* uop)
        {
            if (front == uop) {
                front = front->next;
                return;
            }
            searchRemove(uop);
        }

        inline unary_operation* findOperation(const forest* argF, const forest* resF)
        {
            if (!front) return nullptr;
            if ((front->argF == argF) && (front->resF == resF)) return front;
            return mtfUnary(argF, resF);
        }

        inline unary_operation* findOperation(const forest* argF, opnd_type resType)
        {
            if (!front) return nullptr;
            if ((front->argF == argF) && (front->resultType == resType)) return front;
            return mtfUnary(argF, resType);
        }

    private:
        unary_operation* mtfUnary(const forest* argF, const forest* resF);
        unary_operation* mtfUnary(const forest* argF, opnd_type resType);
        void searchRemove(unary_operation* uop);
};

// ******************************************************************
// *                                                                *
// *                          Unary  apply                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

    typedef unary_operation* (*unary_builtin1)(forest* arg, forest* res);
    typedef unary_operation* (*unary_builtin2)(forest* arg, opnd_type res);

    // ******************************************************************
    // *                          Unary  apply                          *
    // ******************************************************************

    /** Apply a unary operator.
        The operand and the result are not necessarily in the same forest,
        but they must belong to forests that share the same domain.
        This is useful, for instance, for copying a function to a new forest.
            @param  bu  Built-in unary operator (function)
            @param  a   Operand.
            @param  c   Output parameter: the result,
                        where \a c = \a op \a a.
    */
    inline void apply(unary_builtin1 bu, const dd_edge &a, dd_edge &c)
    {
        unary_operation* uop = bu(a.getForest(), c.getForest());
        uop->compute(a, c);
    }

    /** Apply a unary operator.
        For operators whose result is an integer.
            @param  bu  Built-in unary operator (function)
            @param  a   Operand.
            @param  c   Output parameter: the result,
                          where \a c = \a op \a a.
    */
    inline void apply(unary_builtin2 bu, const dd_edge &a, long &c)
    {
        unary_operation* uop = bu(a.getForest(), opnd_type::INTEGER);
        uop->compute(a, c);
    }

    /** Apply a unary operator.
        For operators whose result is a real.
            @param  bu  Built-in unary operator (function)
            @param  a   Operand.
            @param  c   Output parameter: the result,
                          where \a c = \a op \a a.
    */
    inline void apply(unary_builtin2 bu, const dd_edge &a, double &c)
    {
        unary_operation* uop = bu(a.getForest(), opnd_type::REAL);
        uop->compute(a, c);
    }

    /** Apply a unary operator.
        For operators whose result is a anything else.
            @param  bu  Built-in unary operator (function)
            @param  a   Operand.
            @param  rt  Type of result
            @param  c   Output parameter: the result,
                          where \a c = \a op \a a.
    */
    inline void apply(unary_builtin2 bu, const dd_edge &a, opnd_type rt,
            ct_object &c)
    {
        unary_operation* uop = bu(a.getForest(), rt);
        uop->compute(a, c);
    }

#ifdef __GMP_H__
    /** Apply a unary operator.
        For operators whose result is an arbitrary-precision integer
        (as supplied by the GNU MP library).
            @param  bu  Built-in unary operator (function)
            @param  a   Operand.
            @param  c   Input: an initialized MP integer.
                        Output: the result, where \a c = \a op \a a.
    */
    inline void apply(unary_builtin2 bu, const dd_edge &a, mpz_t &c)
    {
        ct_object& x = get_mpz_wrapper();
        unary_operation* uop = bu(a.getForest(), opnd_type::HUGEINT);
        uop->compute(a, x);
        unwrap(x, c);
    }
#endif

#ifdef ALLOW_DEPRECATED_0_17_5
    inline unary_operation* getOperation(unary_builtin1 bu,
            const dd_edge &a, const dd_edge &c)
    {
        return bu(a.getForest(), c.getForest());
    }

    inline unary_operation* getOperation(unary_builtin2 bu,
            const dd_edge &arg, opnd_type res)
    {
        return bu(arg.getForest(), res);
    }
#endif

};



#endif // #include guard

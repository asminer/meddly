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

#include "oper.h"
#include "forest.h"

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

    protected:
        // Fairly standard checks; call these in the operator's constructor.
        //
        /// Make sure all three domains are the same.
        inline void checkDomains(const char* file, unsigned line) const
        {
            if  ( (argF->getDomain() != resF->getDomain()) )
            {
                throw error(error::DOMAIN_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments set/relation status matches
        inline void checkAllRelations(const char* file, unsigned line,
                set_or_rel a) const
        {
            if  (
                    (argF->isForRelations() != a)  ||
                    (resF->isForRelations() != a)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure all arguments match the edge labeling rule
        inline void checkAllLabelings(const char* file, unsigned line,
                edge_labeling a) const
        {
            if  (
                    (argF->getEdgeLabeling() != a)  ||
                    (resF->getEdgeLabeling() != a)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments edge labeling rules match
        inline void checkLabelings(const char* file, unsigned line,
                edge_labeling a, edge_labeling r) const
        {
            if  (
                    (argF->getEdgeLabeling() != a)  ||
                    (resF->getEdgeLabeling() != r)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments match the range type
        inline void checkAllRanges(const char* file, unsigned line,
                range_type rt) const
        {
            if  (
                    (argF->getRangeType() != rt)  ||
                    (resF->getRangeType() != rt)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments match the range type
        inline void checkRanges(const char* file, unsigned line,
                range_type a, range_type r) const
        {
            if  (
                    (argF->getRangeType() != a)  ||
                    (resF->getRangeType() != r)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }

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

        // friend void destroyOperation(unary_operation* &op);
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

        void reset(const char* n);

        inline const char* getName() const { return name; }
        inline bool isEmpty() const { return (!front); }

        inline unary_operation* add(unary_operation* uop) {
            if (uop) {
                uop->next = front;
                front = uop;
            }
            return uop;
        }

        inline void remove(unary_operation* uop)
        {
            if (front == uop) {
                front = front->next;
                return;
            }
            searchRemove(uop);
        }

        inline unary_operation* find(const forest* argF, const forest* resF)
        {
            if (!front) return nullptr;
            if ((front->argF == argF) && (front->resF == resF)) return front;
            return mtfUnary(argF, resF);
        }

        inline unary_operation* find(const forest* argF, opnd_type resType)
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

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
#include "oper_item.h"
#include "forest.h"

#define OLD_UNARY

namespace MEDDLY {
    class oper_item;
    class unary_operation;
    class unary_factory;
    class forest;

    typedef unary_factory&  (*unary_builtin0)();

#ifdef OLD_UNARY
    typedef unary_operation* (*unary_builtin1)(forest* arg, forest* res);
    typedef unary_operation* (*unary_builtin2)(forest* arg, opnd_type res);
    class unary_list;
    class ct_object;
#endif
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
#ifdef ALLOW_DEPRECATED_0_17_6
        /// Old constructor, returns DD
        unary_operation(unary_list& owner, unsigned et_slots,
            forest* arg, forest* res);

        /// Old constructor, returns a number
        unary_operation(unary_list& owner, unsigned et_slots,
            forest* arg, opnd_type res);
#endif
        /** New constructor for a unary op that returns a DD.
            Compute table entry information is managed 'by hand'
            in the derived class.
                @param  arg     Forest containing the argument.

                @param  res     Forest containing the result.
                                Must be compatible with the argument forest.
        */
        unary_operation(forest* arg, forest* res);

        /** New constructor for a unary op that returns a DD.
            Compute table entry information is managed 'by hand'
            in the derived class.
                @param  arg     Forest containing the argument.

                @param  res     Type of the result; should NOT be FOREST.
        */
        unary_operation(forest* arg, opnd_type res);

    protected:
        virtual ~unary_operation();

    protected:
        // Fairly standard checks; call these in the operator's constructor.
        //
        /// Make sure all three domains are the same.
        inline void checkDomains(const char* file, unsigned line) const
        {
            MEDDLY_DCASSERT(argF);
            MEDDLY_DCASSERT(resF);
            if  ( (argF->getDomain() != resF->getDomain()) )
            {
                throw error(error::DOMAIN_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments set/relation status match each other
        inline void checkAllRelations(const char* file, unsigned line) const
        {
            MEDDLY_DCASSERT(argF);
            if  (argF->isForRelations() != resF->isForRelations())
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments set/relation status matches a
        inline void checkAllRelations(const char* file, unsigned line,
                set_or_rel a) const
        {
            MEDDLY_DCASSERT(argF);
            if  (
                    (argF->isForRelations() != a)  ||
                    (resF && resF->isForRelations() != a)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure all arguments match the edge labeling rule
        inline void checkAllLabelings(const char* file, unsigned line,
                edge_labeling a) const
        {
            MEDDLY_DCASSERT(argF);
            if  (
                    (argF->getEdgeLabeling() != a)  ||
                    (resF && resF->getEdgeLabeling() != a)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments edge labeling rules match
        inline void checkLabelings(const char* file, unsigned line,
                edge_labeling a, edge_labeling r) const
        {
            MEDDLY_DCASSERT(argF);
            MEDDLY_DCASSERT(resF);
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
            MEDDLY_DCASSERT(argF);
            if  (
                    (argF->getRangeType() != rt)  ||
                    (resF && resF->getRangeType() != rt)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments match the range type
        inline void checkRanges(const char* file, unsigned line,
                range_type a, range_type r) const
        {
            MEDDLY_DCASSERT(argF);
            MEDDLY_DCASSERT(resF);
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

#ifdef ALLOW_DEPRECATED_0_17_6
        void computeTemp(const dd_edge &arg, dd_edge &res);
        virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);
#endif

        /**
            New virtual compute method for DDs.

                @param  L       Recursion level.
                                If all forests are quasi-reduced,
                                then this is the top level of the
                                operand and result.
                                Ignored for some reduction rules (e.g.,
                                fully reduced) but important for others
                                (e.g., quasi reduced).

                @param  in      Incoming edge index.
                                Important only for identity-reduced
                                relations when L is positive.
                                Use ~0 if there is no edge index.

                @param  av      Edge value for operand
                @param  ap      Node for operand, must be below L.

                @param  cv      Edge value of result
                @param  cp      Node for result, will be below L.
         */
        virtual void compute(
                int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp);

        /**
            New virtual compute method for all other result types.
                @param  L       Recursion level.
                @param  in      Incoming edge index.
                @param  av      Edge value for operand
                @param  ap      Node for operand
                @param  res     Result of the operation is stored here.
         */
        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap, oper_item &res);


        inline void compute(const dd_edge &arg, long &res) {
            oper_item tmp(res);
            compute(argF->getMaxLevelIndex(), ~0,
                    arg.getEdgeValue(), arg.getNode(), tmp);
            res = tmp.getInteger();
        }

        inline void compute(const dd_edge &arg, double &res) {
            oper_item tmp(res);
            compute(argF->getMaxLevelIndex(), ~0,
                    arg.getEdgeValue(), arg.getNode(), tmp);
            res = tmp.getReal();
        }

#ifdef HAVE_LIBGMP
        inline void compute(const dd_edge &arg, mpz_ptr v) {
            oper_item tmp(v);
            compute(argF->getMaxLevelIndex(), ~0,
                    arg.getEdgeValue(), arg.getNode(), tmp);
        }
#endif

        inline void compute(const dd_edge &arg, oper_item &res) {
            compute(argF->getMaxLevelIndex(), ~0,
                    arg.getEdgeValue(), arg.getNode(), res);
        }

    protected:
        inline bool checkForestCompatibility() const
        {
            if (resultType == opnd_type::FOREST) {
                auto o1 = argF->variableOrder();
                auto o2 = resF->variableOrder();
                return o1->is_compatible_with(*o2);
            } else {
                return true;
            }
        }

    protected:
        forest* argF;
        forest* resF;
        opnd_type resultType;

    private:
        /// Who built us
        unary_factory* factory;
        unary_list* parent;
        unary_operation* next;
#ifdef ALLOW_DEPRECATED_0_17_6
        bool new_style;
#endif

        friend class unary_factory;
        friend class unary_list;
};

// ******************************************************************
// *                                                                *
// *                      unary_factory  class                      *
// *                                                                *
// ******************************************************************

/**
    Mechanism to create specific unary operations.
    Derive a subclass from this, and override methods setup(),
        one of the build() methods, and optionally cleanup().
 */
class MEDDLY::unary_factory {
    public:
        // EMPTY constructor, always
        unary_factory() { }

        // EMPTY destructor, always
        ~unary_factory() { }

        /**
            Must be overridden.
            Essentially, the constructor, but called explicitly
            during library initialization.
        */
        virtual void setup() = 0;

        /**
            Build a specific operation instance.
            For this version, the result is another DD.
            The default behavior returns NULL.
        */
        virtual unary_operation* build(forest* arg, forest* res);

        /**
            Build a specific operation instance.
            For this version, the result is a value.
            The default behavior returns NULL.
        */
        virtual unary_operation* build(forest* arg, opnd_type res);

        /// Default: just calls _cleanup()
        virtual void cleanup();

        inline const char* getName() const { return _name; }
        inline const char* getDocs() const { return _doc; }

        /// Should ONLY be called in unary_operation destructor
        inline void remove(unary_operation* uop)
        {
            if (front == uop) {
                front = front->next;
                return;
            }
            searchRemove(uop);
        }


        /** Apply a unary operator.
            The operand and the result are not necessarily in the same forest,
            but they must belong to forests that share the same domain.
            This is useful, for instance, for copying a function
            from one forest to another.
                @param  a   Operand.
                @param  c   Output parameter: the result,
                            where \a c = \a op \a a.
        */
        inline void apply(const dd_edge &a, dd_edge &c)
        {
            unary_operation* uop = build(a.getForest(), c.getForest());
            if (!uop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            uop->compute(a, c);
        }

        /** Apply a unary operator.
            For operators whose result is an integer.
                @param  a   Operand.
                @param  c   Output parameter: the result,
                            where \a c = \a op \a a.
        */
        inline void apply(const dd_edge &a, long &c)
        {
            unary_operation* uop = build(a.getForest(), opnd_type::INTEGER);
            if (!uop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            uop->compute(a, c);
        }

        /** Apply a unary operator.
            For operators whose result is a real.
                @param  a   Operand.
                @param  c   Output parameter: the result,
            where \a c = \a op \a a.
        */
        inline void apply(const dd_edge &a, double &c)
        {
            unary_operation* uop = build(a.getForest(), opnd_type::REAL);
            if (!uop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            uop->compute(a, c);
        }

        /** Apply a unary operator.
            For operators whose result is a real.
                @param  a   Operand.
                @param  c   Output parameter: the result,
                            where \a c = \a op \a a.
        */
        inline void apply(const dd_edge &a, oper_item &c)
        {
            unary_operation* uop = build(a.getForest(), c.getType());
            if (!uop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            uop->compute(a, c);
        }

#ifdef HAVE_LIBGMP
        /** Apply a unary operator.
            For operators whose result is an arbitrary-precision integer
            (as supplied by the GNU MP library).
                @param  a   Operand.
                @param  c   Input: an initialized MP integer.
                            Output: the result, where \a c = \a op \a a.
        */
        inline void apply(const dd_edge &a, mpz_ptr c)
        {
            unary_operation* uop = build(a.getForest(), opnd_type::HUGEINT);
            if (!uop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            uop->compute(a, c);
        }
#endif


    protected:
        /** Initialize base class.
            Should be called by method setup().
            Called during library initialization.
                @param  name    Name of operator.
                @param  doc     Documentation.
        */
        void _setup(const char* name, const char* doc);

        /** Clean up base class.
            Called during library cleanup.
        */
        void _cleanup();

        /// Add an operation to the list.
        inline unary_operation* add(unary_operation* uop) {
            if (uop) {
                MEDDLY_DCASSERT(nullptr == uop->factory);
                uop->factory = this;
                uop->setName(_name);
                uop->next = front;
                front = uop;
            }
            return uop;
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

    private:
        const char* _name;
        const char* _doc;
        unary_operation* front;

};

#ifdef OLD_UNARY

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
                if (uop->parent) {
                    // REMOVE EVENTUALLY
                    MEDDLY_DCASSERT(uop->parent == this);
                } else {
                    uop->parent = this;
                    uop->setName(name);
                }
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

#endif

// ******************************************************************
// *                                                                *
// *                          Unary  apply                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

    // ******************************************************************
    // *                          Unary  apply                          *
    // ******************************************************************
    template <class RES_T>
    inline void apply(unary_builtin0 bu, const dd_edge &a, RES_T &c)
    {
        bu().apply(a, c);
    }

#ifdef HAVE_LIBGMP
    inline void apply(unary_builtin0 bu, const dd_edge &a, mpz_ptr c)
    {
        bu().apply(a, c);
    }
#endif

    inline unary_operation* build(unary_builtin0 bu, forest* arg, forest* res)
    {
        return bu().build(arg, res);
    }

    inline unary_operation* build(unary_builtin0 bu, forest* arg, opnd_type res)
    {
        return bu().build(arg, res);
    }

#ifdef OLD_UNARY

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
        if (!uop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
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
        if (!uop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
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
        if (!uop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        uop->compute(a, c);
    }

    /** Apply a unary operator.
        For operators whose result is a real.
            @param  bu  Built-in unary operator (function)
            @param  a   Operand.
            @param  c   Output parameter: the result,
                          where \a c = \a op \a a.
    */
    inline void apply(unary_builtin2 bu, const dd_edge &a, oper_item &c)
    {
        unary_operation* uop = bu(a.getForest(), c.getType());
        if (!uop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        uop->compute(a, c);
    }

#ifdef HAVE_LIBGMP
    /** Apply a unary operator.
        For operators whose result is an arbitrary-precision integer
        (as supplied by the GNU MP library).
            @param  bu  Built-in unary operator (function)
            @param  a   Operand.
            @param  c   Input: an initialized MP integer.
                        Output: the result, where \a c = \a op \a a.
    */
    inline void apply(unary_builtin2 bu, const dd_edge &a, mpz_ptr c)
    {
        unary_operation* uop = bu(a.getForest(), opnd_type::HUGEINT);
        if (!uop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        uop->compute(a, c);
    }
#endif

#endif

};



#endif // #include guard

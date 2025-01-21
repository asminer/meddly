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

#ifndef MEDDLY_OPER_BINARY_H
#define MEDDLY_OPER_BINARY_H

#include "oper.h"
#include "forest.h"

namespace MEDDLY {
    class dd_edge;
    class binary_operation;
    class binary_list;
};

// ******************************************************************
// *                                                                *
// *                     binary_operation class                     *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a binary operation in a specific forest.
    Specific operations will be derived from this class.
*/
class MEDDLY::binary_operation : public operation {
    public:
#ifdef ALLOW_DEPRECATED_0_17_6
        /**
            OLD constructor.
         */
        binary_operation(binary_list& owner, unsigned et_slots,
            forest* arg1, forest* arg2, forest* res);
#endif

        /**
            NEW constructor.
            Compute table entry information is managed 'by hand'
            in the derived class.
                @param  arg1    Forest containing the first argument
                                of the binary operation.

                @param  arg2    Forest containing the second argument
                                of the binary operation. May equal or
                                not arg1, and if unequal they must
                                be "compatible" (same domain, same
                                variable order).

                @param  res     Forest that holds the result.
                                Must be compatible with the argument
                                forests.
        */
        binary_operation(forest* arg1, forest* arg2, forest* res);

    protected:
        virtual ~binary_operation();

    protected:
        // Fairly standard checks; call these in the operator's constructor.
        //
        /// Make sure all three domains are the same.
        inline void checkDomains(const char* file, unsigned line) const
        {
            if  (
                    (arg1F->getDomain() != resF->getDomain()) ||
                    (arg2F->getDomain() != resF->getDomain())
                )
            {
                throw error(error::DOMAIN_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments set/relation status match each other
        inline void checkAllRelations(const char* file, unsigned line) const
        {
            if  (
                    (arg1F->isForRelations() != resF->isForRelations())  ||
                    (arg2F->isForRelations() != resF->isForRelations())
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments set/relation status matches
        inline void checkAllRelations(const char* file, unsigned line,
                set_or_rel a) const
        {
            if  (
                    (arg1F->isForRelations() != a)  ||
                    (arg2F->isForRelations() != a)  ||
                    (resF->isForRelations() != a)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments set/relation status matches
        inline void checkRelations(const char* file, unsigned line,
                set_or_rel a1, set_or_rel a2, set_or_rel r) const
        {
            if  (
                    (arg1F->isForRelations() != a1)  ||
                    (arg2F->isForRelations() != a2)  ||
                    (resF->isForRelations() != r)
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
                    (arg1F->getEdgeLabeling() != a)  ||
                    (arg2F->getEdgeLabeling() != a)  ||
                    (resF->getEdgeLabeling() != a)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments edge labeling rules match
        inline void checkLabelings(const char* file, unsigned line,
                edge_labeling a1, edge_labeling a2, edge_labeling r) const
        {
            if  (
                    (arg1F->getEdgeLabeling() != a1)  ||
                    (arg2F->getEdgeLabeling() != a2)  ||
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
                    (arg1F->getRangeType() != rt)  ||
                    (arg2F->getRangeType() != rt)  ||
                    (resF->getRangeType() != rt)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments match the edge type
        inline void checkAllEdgeTypes(const char* file, unsigned line,
                edge_type et) const
        {
            if  (
                    (arg1F->getEdgeType() != et)  ||
                    (arg2F->getEdgeType() != et)  ||
                    (resF->getEdgeType() != et)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }

    public:
        /**
            Checks forest comatability and then calls computeDDEdge() (OLD)
            or the new virtual compute method (NEW).
        */
        void compute(const dd_edge &ar1, const dd_edge &ar2,
                dd_edge &res);

#ifdef ALLOW_DEPRECATED_0_17_6
        void computeTemp(const dd_edge &ar1, const dd_edge &ar2,
                dd_edge &res);

        virtual void computeDDEdge(const dd_edge &ar1, const dd_edge &ar2,
                dd_edge &res, bool userFlag);
#endif

        /**
            New virtual compute method.

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

                @param  av      Edge value for operand 1
                @param  ap      Node for operand 1, must be below L.

                @param  bv      Edge value for operand 2
                @param  bp      Node for operand 2, must be below L.

                @param  cv      Edge value of result
                @param  cp      Node for result, will be below L.
         */
        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp)
#ifndef ALLOW_DEPRECATED_0_17_6
                = 0
#endif
                ;

#ifdef ALLOW_DEPRECATED_0_17_6
        /**
            During the transition, callers may need to know
            which version of compute() to call.
         */
        inline bool useNewCompute() const { return new_style; }
#endif

    protected:
#ifdef ALLOW_DEPRECATED_0_17_6
        inline bool canCommute() const {
            return can_commute;
        }
        /// Call this in the constructor of any operation that commutes.
        inline void operationCommutes() {
            can_commute = (arg1F == arg2F);
        }
#endif

        // Check if the variables orders of relevant forests are compatible
        inline bool checkForestCompatibility() const
        {
            auto o1 = arg1F->variableOrder();
            auto o2 = arg2F->variableOrder();
            auto o3 = resF->variableOrder();
            return o1->is_compatible_with(*o2) && o1->is_compatible_with(*o3);
        }

    protected:
        forest* arg1F;
        forest* arg2F;
        forest* resF;

    private:
        binary_list* parent;
        binary_operation* next;
#ifdef ALLOW_DEPRECATED_0_17_6
        bool can_commute;
        bool new_style;
#endif

        friend class binary_list;
};

// ******************************************************************
// *                                                                *
// *                       binary_list  class                       *
// *                                                                *
// ******************************************************************

/**
    List of binary operations of the same type; used when
    building binary operations for specific forests.
*/
class MEDDLY::binary_list {
        const char* name;
        binary_operation* front;
    public:
        binary_list(const char* n = nullptr);

        void reset(const char* n);

        inline const char* getName() const { return name; }
        inline bool isEmpty() const { return (!front); }

        inline binary_operation* add(binary_operation* bop) {
            if (bop) {
                if (bop->parent) {
                    // REMOVE EVENTUALLY
                    MEDDLY_DCASSERT(bop->parent == this);
                } else {
                    bop->parent = this;
                    bop->setName(name);
                }
                bop->next = front;
                front = bop;
            }
            return bop;
        }

        inline void remove(binary_operation* bop)
        {
            if (front == bop) {
                front = front->next;
                return;
            }
            searchRemove(bop);
        }

        inline binary_operation* find(const forest* arg1F,
                const forest* arg2F, const forest* resF)
        {
            if (!front) return nullptr;
            if ((front->arg1F == arg1F) && (front->arg2F == arg2F) && (front->resF == resF)) return front;
            return mtfBinary(arg1F, arg2F, resF);
        }

    private:
        binary_operation* mtfBinary(const forest* arg1F, const forest* arg2F, const forest* resF);
        void searchRemove(binary_operation* uop);
};


// ******************************************************************
// *                                                                *
// *                          Binary apply                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

    typedef binary_operation* (*binary_builtin)(forest* arg1,
            forest* arg2, forest* res);

    /** Apply a binary operator.
        \a a, \a b and \a c are not required to be in the same forest,
        but they must have the same domain. The result will be in the
        same forest as \a result. The operator decides the type of forest
        for each \a dd_edge.
        Useful, for example, for constructing comparisons
        where the resulting type is "boolean" but the operators are not,
        e.g., c = f EQUALS g.
            @param  bb    Built-in binary operation.
            @param  a     First operand.
            @param  b     Second operand.
            @param  c     Output parameter: the result,
                          where \a c = \a a \a op \a b.
    */
    inline void apply(binary_builtin bb, const dd_edge &a, const dd_edge &b,
        dd_edge &c)
    {
        binary_operation* bop = bb(a.getForest(), b.getForest(), c.getForest());
        bop->compute(a, b, c);
    }

#ifdef ALLOW_DEPRECATED_0_17_5
    inline binary_operation* getOperation(binary_builtin bb,
            const dd_edge &a, const dd_edge &b, const dd_edge &c)
    {
        return bb(a.getForest(), b.getForest(), c.getForest());
    }
#endif

};

#endif // #include guard

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

#ifndef MEDDLY_OPER_TERNARY_H
#define MEDDLY_OPER_TERNARY_H

#include "oper.h"
#include "forest.h"

namespace MEDDLY {
    class dd_edge;
    class ternary_operation;
    class ternary_list;
};

// ******************************************************************
// *                                                                *
// *                    ternary_operation  class                    *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a ternary operation in a specific forest.
    Examples include constrained saturation (for now, anyway)
    and IF-THEN-ELSE (tbd).
    Specific operations will be derived from this class.
*/
class MEDDLY::ternary_operation : public operation {
    public:
        ternary_operation(ternary_list& owner, unsigned et_slots,
            forest* arg1, forest* arg2, forest* arg3, forest* res);

    protected:
        virtual ~ternary_operation();

    protected:
        // Fairly standard checks; call these in the operator's constructor.
        //
        /// Make sure all domains are the same.
        inline void checkDomains(const char* file, unsigned line) const
        {
            if  (
                    (arg1F->getDomain() != resF->getDomain()) ||
                    (arg2F->getDomain() != resF->getDomain()) ||
                    (arg3F->getDomain() != resF->getDomain())
                )
            {
                throw error(error::DOMAIN_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments set/relation status matches
        inline void checkAllRelations(const char* file, unsigned line,
                set_or_rel a) const
        {
            if  (
                    (arg1F->isForRelations() != a)  ||
                    (arg2F->isForRelations() != a)  ||
                    (arg3F->isForRelations() != a)  ||
                    (resF->isForRelations() != a)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments set/relation status matches
        inline void checkRelations(const char* file, unsigned line,
                set_or_rel a1, set_or_rel a2, set_or_rel a3, set_or_rel r) const
        {
            if  (
                    (arg1F->isForRelations() != a1)  ||
                    (arg2F->isForRelations() != a2)  ||
                    (arg3F->isForRelations() != a3)  ||
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
                    (arg3F->getEdgeLabeling() != a)  ||
                    (resF->getEdgeLabeling() != a)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments edge labeling rules match
        inline void checkLabelings(const char* file, unsigned line,
                edge_labeling a1, edge_labeling a2, edge_labeling a3,
                edge_labeling r) const
        {
            if  (
                    (arg1F->getEdgeLabeling() != a1)  ||
                    (arg2F->getEdgeLabeling() != a2)  ||
                    (arg3F->getEdgeLabeling() != a3)  ||
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
                    (arg3F->getRangeType() != rt)  ||
                    (resF->getRangeType() != rt)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }

    public:
        /**
            Checks forest comatability and then calls computeDDEdge().
        */
        void compute(const dd_edge &ar1, const dd_edge &ar2,
                const dd_edge &ar3, dd_edge &res);
        void computeTemp(const dd_edge &ar1, const dd_edge &ar2,
                const dd_edge &ar3, dd_edge &res);

        virtual void computeDDEdge(const dd_edge &ar1, const dd_edge &ar2,
                const dd_edge &ar3, dd_edge &res, bool userFlag) = 0;

    protected:
        // Check if the variables orders of relevant forests are compatible
        virtual bool checkForestCompatibility() const;

    protected:
        forest* arg1F;
        forest* arg2F;
        forest* arg3F;
        forest* resF;

    private:
        ternary_list& parent;
        ternary_operation* next;

        friend class ternary_list;
};

// ******************************************************************
// *                                                                *
// *                       ternary_list class                       *
// *                                                                *
// ******************************************************************

/**
    List of ternary operations of the same type; used when
    building ternary operations for specific forests.
*/
class MEDDLY::ternary_list {
        const char* name;
        ternary_operation* front;
    public:
        ternary_list(const char* n = nullptr);

        void reset(const char* n);

        inline const char* getName() const { return name; }
        inline bool isEmpty() const { return (!front); }

        inline ternary_operation* add(ternary_operation* top) {
            if (top) {
                top->next = front;
                front = top;
            }
            return top;
        }

        inline void remove(ternary_operation* top)
        {
            if (front == top) {
                front = front->next;
                return;
            }
            searchRemove(top);
        }

        inline ternary_operation* find(const forest* arg1F,
                const forest* arg2F, const forest* arg3F, const forest* resF)
        {
            if (!front) return nullptr;
            if ((front->arg1F == arg1F) &&
                (front->arg2F == arg2F) &&
                (front->arg3F == arg3F) &&
                (front->resF == resF)) return front;
            return mtfTernary(arg1F, arg2F, arg3F, resF);
        }

    private:
        ternary_operation* mtfTernary(const forest* arg1F, const forest* arg2F,
            const forest* arg3F, const forest* resF);
        void searchRemove(ternary_operation* top);
};


// ******************************************************************
// *                                                                *
// *                         Ternary  apply                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

    typedef ternary_operation* (*ternary_builtin)(forest* arg1,
            forest* arg2, forest* arg3, forest* res);

    /** Apply a ternary operator.
        \a a, \a b, \a c, and \a d are not required to be in the same forest,
        but they must have the same domain. The operator decides the type of
        forest for each \a dd_edge.
            @param  bt    Built-in ternary operation.
            @param  a     First operand.
            @param  b     Second operand.
            @param  c     Third operand.
            @param  d     Output parameter: the result,
                          where \a c = \a a \a op \a b.
    */
    inline void apply(ternary_builtin bt, const dd_edge &a, const dd_edge &b,
        const dd_edge &c, dd_edge &d)
    {
        ternary_operation* top = bt(a.getForest(), b.getForest(),
            c.getForest(), d.getForest());
        top->compute(a, b, c, d);
    }

};

#endif // #include guard

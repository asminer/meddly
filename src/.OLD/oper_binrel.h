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

#ifndef MEDDLY_OPER_BINREL_H
#define MEDDLY_OPER_BINREL_H

#include "oper.h"
#include "forest.h"
#include "relforest.h"

namespace MEDDLY {
    class dd_edge;
    class binrel_operation;
    class binrel_list;
};

// ******************************************************************
// *                                                                *
// *                     binrel_operation class                     *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a binary operation between an MDD and a
    generic relation (from relforest).
    Specific operations will be derived from this class.
*/
class MEDDLY::binrel_operation : public operation {
    public:
        /**
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
        binrel_operation(forest* arg1, relforest* arg2, forest* res);

    protected:
        virtual ~binrel_operation();

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
            if  (arg1F->isForRelations() != resF->isForRelations())
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
                    (arg1F->getEdgeLabeling() != a)  ||
                    (resF->getEdgeLabeling() != a)
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
        void compute(const dd_edge &ar1, const relforest::edge &ar2,
                dd_edge &res);


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
                edge_value &cv, node_handle &cp) = 0;

    protected:
        // Check if the variables orders of relevant forests are compatible
        inline bool checkForestCompatibility() const
        {
            auto o1 = arg1F->variableOrder();
            auto o3 = resF->variableOrder();
            return o1->is_compatible_with(*o3);
        }

    protected:
        forest* arg1F;
        relforest* arg2F;
        forest* resF;

    private:
        binrel_list* parent;
        binrel_operation* next;

        friend class binrel_list;
};

// ******************************************************************
// *                                                                *
// *                       binrel_list  class                       *
// *                                                                *
// ******************************************************************

/**
    List of binrel operations of the same type; used when
    building binrel operations for specific forests.
*/
class MEDDLY::binrel_list {
        const char* name;
        binrel_operation* front;
    public:
        binrel_list(const char* n = nullptr);

        void reset(const char* n);

        inline const char* getName() const { return name; }
        inline bool isEmpty() const { return (!front); }

        inline binrel_operation* add(binrel_operation* bop) {
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

        inline void remove(binrel_operation* bop)
        {
            if (front == bop) {
                front = front->next;
                return;
            }
            searchRemove(bop);
        }

        inline binrel_operation* find(const forest* arg1F,
                const relforest* arg2F, const forest* resF)
        {
            if (!front) return nullptr;
            if ((front->arg1F == arg1F) &&
                (front->arg2F == arg2F) &&
                (front->resF == resF))
            {
                return front;
            }
            return mtfBinary(arg1F, arg2F, resF);
        }

    private:
        binrel_operation* mtfBinary(const forest* arg1F,
                const relforest* arg2F, const forest* resF);
        void searchRemove(binrel_operation* uop);
};


// ******************************************************************
// *                                                                *
// *                          Binary apply                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

    typedef binrel_operation* (*binrel_builtin)(forest* arg1,
            relforest* arg2, forest* res);

    /** Apply a binrel operator.
        \a a and \a c are not required to be in the same forest,
        but they must have the same domain, and have the same domain
        as relation \a b. The operator decides the type of forest
        for each \a dd_edge.
            @param  bb    Built-in binrel operation.
            @param  a     First operand.
            @param  b     Second operand.
            @param  c     Output parameter: the result,
                          where \a c = \a a \a op \a b.
    */
    inline void apply(binrel_builtin bb, const dd_edge &a,
        const relforest::edge &b, dd_edge &c)
    {
        binrel_operation* bop = bb(a.getForest(), b.parent, c.getForest());
        if (!bop) throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        bop->compute(a, b, c);
    }

};

#endif // #include guard

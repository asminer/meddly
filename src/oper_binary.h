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

namespace MEDDLY {
    class dd_edge;
    class binary_operation;

    // ******************************************************************
    // *                      Operation management                      *
    // ******************************************************************

#if 0
    /** Find, or build if necessary, a binary operation.
            @param  code    Operation we want
            @param  arg1    Argument 1 forest
            @param  arg2    Argument 2 forest
            @param  res     Result forest
            @return         The matching operation, if it already exists;
                            a new operation, otherwise.
    */
    binary_operation* getOperation(const binary_opname* code,
        expert_forest* arg1, expert_forest* arg2, expert_forest* res);

    /** Find, or build if necessary, a binary operation.
            @param  code    Operation we want
            @param  arg1    Argument 1 forest taken from this dd_edge
            @param  arg2    Argument 2 forest taken from this dd_edge
            @param  res     Result forest taken from this dd_edge
            @return         The matching operation, if it already exists;
                            a new operation, otherwise.
    */
    binary_operation* getOperation(const binary_opname* code,
        const dd_edge& arg1, const dd_edge& arg2, const dd_edge& res);

#endif

    /** Safely destroy the given binary operation.
        It should be unnecessary to call this directly.
    */
    void destroyOperation(binary_operation* &op);

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
        binary_operation(const binary_opname* code, unsigned et_slots,
            expert_forest* arg1, expert_forest* arg2, expert_forest* res);

    protected:
        virtual ~binary_operation();

    public:
        bool matches(const dd_edge &arg1, const dd_edge &arg2,
                const dd_edge &res) const;

        // high-level front-end

        /**
            Checks forest comatability and then calls computeDDEdge().
        */
        void compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);
        void computeTemp(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);

        virtual void computeDDEdge(const dd_edge &ar1, const dd_edge &ar2,
                dd_edge &res, bool userFlag) = 0;

    protected:
        inline void operationCommutes() {
            can_commute = (arg1F == arg2F);
        }

        // Check if the variables orders of relevant forests are compatible
        virtual bool checkForestCompatibility() const;

    protected:
        bool can_commute;
        expert_forest* arg1F;
        expert_forest* arg2F;
        expert_forest* resF;
        opnd_type resultType;
};


#endif // #include guard

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
#include <list>
#include "markcmp.h"
namespace MEDDLY {
    class dd_edge;
    class binary_operation;
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
        friend void destroyOperation(binary_operation* &op);
    public:
        binary_operation(binary_opname* code, unsigned et_slots,
            expert_forest* arg1, expert_forest* arg2, expert_forest* res);
        binary_operation(binary_opname* code, unsigned et_slots,
            expert_forest* arg1, expert_forest* arg2, expert_forest* res,expert_forest* res2,expert_forest* res3);

    protected:
        virtual ~binary_operation();

    public:
        bool matches(const dd_edge &arg1, const dd_edge &arg2,
                const dd_edge &res) const;
        bool matches(const dd_edge &arg1, const dd_edge &arg2,
                const dd_edge &res,const dd_edge &res2,const dd_edge &res3) const;

        // high-level front-end

        /**
            Checks forest comatability and then calls computeDDEdge().
        */
        void compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);
        void computeTemp(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);

        virtual void computeDDEdge(const dd_edge &ar1, const dd_edge &ar2,
                dd_edge &res, bool userFlag) = 0;
        virtual void computeDDEdgeSC(const dd_edge &ar1, const dd_edge &ar2,
                dd_edge &res, bool userFlag,std::list<int>* shouldConfirm){};
        virtual void computeDDEdgeSC(const dd_edge &ar1, const dd_edge &ar2,
                dd_edge &res, dd_edge &res2, bool userFlag,std::list<int>* shouldConfirm){};
        virtual void computeDDEdgeSC(const dd_edge &ar1, const dd_edge &ar2,
                dd_edge &res, dd_edge &res2, dd_edge &res3, int &res4, bool userFlag,std::list<int>* shouldConfirm, markcmp* cij){};
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
        expert_forest* resF2;
        expert_forest* resF3;
        opnd_type resultType;
};


#endif // #include guard

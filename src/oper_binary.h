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
#include "dd_edge.h"

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
            forest* arg1, forest* arg2, forest* res);

    protected:
        virtual ~binary_operation();

    public:
        inline bool matches(const dd_edge &arg1, const dd_edge &arg2,
                const dd_edge &res) const
        {
            return
                arg1.isAttachedTo(arg1F) &&
                arg2.isAttachedTo(arg2F) &&
                res.isAttachedTo(resF);
        }

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
        forest* arg1F;
        forest* arg2F;
        forest* resF;
        opnd_type resultType;
};


#endif // #include guard

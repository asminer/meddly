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
        binary_operation(binary_list& owner, unsigned et_slots,
            forest* arg1, forest* arg2, forest* res);

    protected:
        virtual ~binary_operation();

    public:
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

    private:
        binary_list& parent;
        binary_operation* next;

        friend class binary_list;
        friend void destroyOperation(binary_operation* &op);
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
        binary_list(const char* n);

        inline const char* getName() const { return name; }

        inline binary_operation* addOperation(binary_operation* bop) {
            if (bop) {
                bop->next = front;
                front = bop;
            }
            return bop;
        }

        inline void removeOperation(binary_operation* bop)
        {
            if (front == bop) {
                front = front->next;
                return;
            }
            searchRemove(bop);
        }

        inline binary_operation* findOperation(const forest* arg1F,
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


#endif // #include guard

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

namespace MEDDLY {
    class dd_edge;
    class ct_object;
    class unary_operation;
    class unary_list;
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
        // bool matches(const dd_edge &arg, const dd_edge &res) const;
        // bool matches(const dd_edge &arg, opnd_type res) const;

        // high-level front-ends

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
        unary_list(const char* n);

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

#endif // #include guard

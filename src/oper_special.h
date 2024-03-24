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

#ifndef MEDDLY_OPER_SPECIAL_H
#define MEDDLY_OPER_SPECIAL_H

#include "oper.h"

namespace MEDDLY {
    class dd_edge;
    class specialized_operation;

    /// Safely destroy the given specialized operation.
    void destroyOperation(specialized_operation* &op);
};

// ******************************************************************
// *                                                                *
// *                  specialized_operation  class                  *
// *                                                                *
// ******************************************************************

/** Mechanism to apply specialized operations.
*/
class MEDDLY::specialized_operation : public operation {
    public:
        specialized_operation(const char* name, unsigned et_slots);
    protected:
        virtual ~specialized_operation();
    public:
        /** For unary (like) operations.
            Note that there could be other "built in" operands.
            Default behavior is to throw an exception.
        */
        virtual void compute(const dd_edge &arg, dd_edge &res);

        /** For unary (like) operations with boolean results.
            Note that there could be other "built in" operands.
            Default behavior is to throw an exception.
        */
        virtual void compute(const dd_edge &arg, bool &res);

        /** For binary (like) operations.
            Note that there could be other "built in" operands.
            Default behavior is to throw an exception.
        */
        virtual void compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);

        /** For numerical operations.
            compute y += some function of x, depending on the operation.
            Default behavior is to throw an exception.
        */
        virtual void compute(double* y, const double* x);

        /** For tenary (like) operations.
            Note that there could be other "built in" operands.
            Default behavior is to throw an exception.
        */
        virtual void compute(const dd_edge &ar1, const dd_edge &ar2, const dd_edge &ar3, dd_edge &res);

    private:

        friend void destroyOperation(specialized_operation* &op);
};

#endif // #include guard

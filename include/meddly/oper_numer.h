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

#ifndef MEDDLY_OPER_NUMER_H
#define MEDDLY_OPER_NUMER_H

#include "oper.h"

namespace MEDDLY {
    class dd_edge;
    class numerical_operation;
};

// ******************************************************************
// *                                                                *
// *                   numerical_operation  class                   *
// *                                                                *
// ******************************************************************

/**
    Mechanism to apply numerical operations.
*/
class MEDDLY::numerical_operation : public operation {
    public:
        numerical_operation(const char* name);
    protected:
        virtual ~numerical_operation();
    public:
        /** For numerical operations.
            compute y += some function of x, depending on the operation.
        */
        virtual void compute(double* y, const double* x) = 0;

};

#endif // #include guard

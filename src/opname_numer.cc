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

#include "defines.h"
#include "opname_numer.h"
#include "initializer.h"

//
// For initializing numerical operations
//
#include "operations/vect_matr.h"


// ******************************************************************
// *                                                                *
// *                  numerical_opname_init  class                  *
// *                                                                *
// ******************************************************************

class MEDDLY::numerical_opname_init : public initializer_list {
    public:
        numerical_opname_init(initializer_list* p);
    protected:
        virtual void setup();
        virtual void cleanup();
    private:
        inline void mydelete(numerical_opname* &p) {
            delete p;
            p = nullptr;
        }
};

MEDDLY::numerical_opname_init::numerical_opname_init(initializer_list* p)
    : initializer_list(p)
{
    numerical_opname::_EXPLVECT_MATR_MULT = nullptr;
    numerical_opname::_MATR_EXPLVECT_MULT = nullptr;
}

void MEDDLY::numerical_opname_init::setup()
{
    numerical_opname::_EXPLVECT_MATR_MULT   =   initExplVectorMatrixMult()  ;
    numerical_opname::_MATR_EXPLVECT_MULT   =   initMatrixExplVectorMult()  ;
}

void MEDDLY::numerical_opname_init::cleanup()
{
    mydelete(numerical_opname::_EXPLVECT_MATR_MULT);
    mydelete(numerical_opname::_MATR_EXPLVECT_MULT);
}

// ******************************************************************
// *                                                                *
// *                    numerical_opname methods                    *
// *                                                                *
// ******************************************************************

MEDDLY::numerical_opname* MEDDLY::numerical_opname::_EXPLVECT_MATR_MULT;
MEDDLY::numerical_opname* MEDDLY::numerical_opname::_MATR_EXPLVECT_MULT;

MEDDLY::numerical_opname::numerical_args
::numerical_args(const dd_edge &xi, const dd_edge &a, const dd_edge &yi)
 : x_ind(xi), A(a), y_ind(yi)
{
}

MEDDLY::numerical_opname::numerical_args::~numerical_args()
{
}


MEDDLY::numerical_opname::numerical_opname(const char* n)
 : specialized_opname(n)
{
}

MEDDLY::numerical_opname::~numerical_opname()
{
}

MEDDLY::initializer_list*
MEDDLY::numerical_opname::makeInitializer(initializer_list* prev)
{
    return new numerical_opname_init(prev);
}


// ******************************************************************
// *                                                                *
// *                front end:  numerical operations                *
// *                                                                *
// ******************************************************************

MEDDLY::numerical_opname* MEDDLY::EXPLVECT_MATR_MULT()
{
    return numerical_opname::EXPLVECT_MATR_MULT();
}

MEDDLY::numerical_opname* MEDDLY::MATR_EXPLVECT_MULT()
{
    return numerical_opname::MATR_EXPLVECT_MULT();
}


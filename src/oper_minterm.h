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

#ifndef MEDDLY_OPER_MINTERM_H
#define MEDDLY_OPER_MINTERM_H

#ifdef ALLOW_MINTERM_OPS

#include "oper.h"
#include "forest.h"

namespace MEDDLY {
    class dd_edge;
    class minterm_coll;
    class minterm_operation;
};

// ******************************************************************
// *                                                                *
// *                    minterm_operation  class                    *
// *                                                                *
// ******************************************************************

/** Mechanism to apply an operation between a forest and a
    collection of minterms. Specific operations will be derived
    from this class.
*/
class MEDDLY::minterm_operation : public operation {
    public:
        /**
            Compute table entry information is managed 'by hand'
            in the derived class.
                @param  arg1    Forest containing the first argument
                                of the minterm operation.

                @param  arg2    Minterm collection for the second argument.

                @param  res     Forest that holds the result.
                                Must be compatible with the argument
                                forests.
        */
        minterm_operation(forest* arg1, minterm_coll &arg2, forest* res);

        virtual ~minterm_operation();

    public:
        /**
            Checks forest comatability and then calls the
            virtual compute method.
        */
        void compute(const dd_edge &ar1, dd_edge &res);

        /**
            Virtual compute method.

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

                @param  low     Low element in minterm collection.
                @param  high    1+high element in minterm collection.

                @param  cv      Edge value of result
                @param  cp      Node for result, will be below L.
         */
        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                unsigned low, unsigned high,
                edge_value &cv, node_handle &cp) = 0;

    protected:
        forest* arg1F;
        minterm_coll &arg2;
        forest* resF;
};

// ******************************************************************
// *                                                                *
// *                         Minterm  apply                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

    typedef minterm_operation* (*minterm_builtin)(forest* arg1,
            minterm_coll &arg2, forest* res);

    /** Apply a minterm operator.
            @param  OP    Built-in minterm operation.
            @param  a     First operand.
            @param  b     Second operand: minterm collection.
            @param  c     Output parameter: the result,
                          where \a c = \a a \a op \a b.
    */
    inline void apply(minterm_builtin OP, const dd_edge &a, minterm_coll &b,
        dd_edge &c)
    {
        minterm_operation* mop = OP(a.getForest(), b, c.getForest());
        mop->compute(a, c);
        delete mop;
    }
};

#endif

#endif // #include guard

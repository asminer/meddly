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

#include "../defines.h"
#include "arith_distmin.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

#include "arith_templ.h"

namespace MEDDLY {
    class DISTMIN_factory;
};


// ******************************************************************
// *                                                                *
// *                        mt_distmin class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct mt_distmin {
        inline static const char* name() {
            return "min";
        }
        inline static bool commutes() {
            return true;
        }
        inline static bool stopOnEqualArgs() {
            return true;
        }
        inline static void makeEqualResult(int L, unsigned in,
                const forest* fa, node_handle a,
                forest* fc, edge_value &cv, node_handle &c,
                unary_operation* copier)
        {
            // min(a, a) = a; just copy it
            //
            MEDDLY_DCASSERT(copier);
            edge_value zero;
            copier->compute(L, in, zero, a, cv, c);
        }
        inline static bool simplifiesToFirstArg(int L,
                const forest* fa, node_handle &a,
                const forest* fb, node_handle b)
        {
            return false;
        }
        inline static bool simplifiesToSecondArg(int L,
                const forest* fa, node_handle a,
                const forest* fb, node_handle &b)
        {
            return false;
        }

        inline static void apply(const forest* fa, node_handle a,
                const forest* fb, node_handle b,
                const forest* fc, node_handle &c)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            if (av < 0) {
                if (bv < 0) {
                    c = fc->handleForValue(-1);
                } else {
                    c = b;
                }
            } else {
                if (bv < 0) {
                    c = a;
                } else {
                    c = fc->handleForValue( MIN(av, bv) );
                }
            }
        }

        inline static void apply(const edge_value &a, const edge_value &b,
                                    edge_value &c)
        {
            MEDDLY_DCASSERT(a.isVoid());
            MEDDLY_DCASSERT(b.isVoid());
            c.set();
        }
    };
};


// ******************************************************************
// *                                                                *
// *                     DISTMIN_factory  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::DISTMIN_factory : public binary_factory {
    public:
        virtual void setup();
        virtual binary_operation* build_new(forest* a, forest* b, forest* c);
};

// ******************************************************************

void MEDDLY::DISTMIN_factory::setup()
{
    _setup(__FILE__, "DIST_MIN", "Minimum value, except negatives are treated as infinity. Forest ranges must be integer or real. Forests must be multi-terminal, and over the same domain.");
}

MEDDLY::binary_operation*
MEDDLY::DISTMIN_factory::build_new(forest* a, forest* b, forest* c)
{
    if (c->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            return new arith_compat<EdgeOp_none, mt_distmin<float> > (a,b,c);
        } else {
            return new arith_compat<EdgeOp_none, mt_distmin<long> > (a,b,c);
        }
    }

    return nullptr;
}

// ******************************************************************
// *                                                                *
// *                       DIST_MIN front end                       *
// *                                                                *
// ******************************************************************

MEDDLY::binary_factory& MEDDLY::DIST_MIN()
{
    static DISTMIN_factory F;
    return F;
}


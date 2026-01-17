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
#include "arith_mod.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

#include "arith_templ.h"

namespace MEDDLY {
    class MOD_factory;
};


// ******************************************************************
// *                                                                *
// *                          mt_mod class                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct mt_mod {
        inline static const char* name() {
            return "mod";
        }
        inline static bool commutes() {
            return false;
        }
        inline static bool stopOnEqualArgs() {
            return true;
        }
        inline static void makeEqualResult(int L, unsigned in,
                const forest* fa, node_handle a,
                forest* fc, edge_value &cv, node_handle &c,
                unary_operation* copier)
        {
            //
            // NOTE:
            //   this (perhaps wrongly) assumes that the function
            //   encoded by a has no zero values;
            //   otherwise we're returning 0 for 0%0
            //
            cv.set();
            c = fc->makeRedundantsTo( fc->handleForValue(RANGE(0)), 0, L );
        }
        inline static bool simplifiesToFirstArg(int L,
                const forest* fa, node_handle &a,
                const forest* fb, node_handle b)
        {
            return (0==a);
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
            if (0 == bv) throw error(error::DIVIDE_BY_ZERO, __FILE__, __LINE__);
            c = fc->handleForValue( av % bv );
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
// *                        evplus_mod class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EDGETYPE>
    struct evplus_mod {
        inline static const char* name() {
            return "mod";
        }
        inline static bool commutes() {
            return false;
        }
        inline static bool stopOnEqualArgs() {
            return true;
        }
        inline static void makeEqualResult(int L, unsigned in,
                const edge_value &av, const node_handle &a,
                forest* fc, edge_value &cv, node_handle &c,
                unary_operation* copyop)
        {
            //
            // NOTE:
            //   this (perhaps wrongly) assumes that the function
            //   encoded by a has no zero values;
            //   otherwise we're returning 0 for 0%0
            //
            cv = EDGETYPE(0);
            c = fc->makeRedundantsTo( OMEGA_NORMAL, 0, L );
        }
        inline static bool simplifiesToFirstArg(int L,
                const forest* f1, edge_value &av, node_handle &an,
                const forest* f2, const edge_value &bv, node_handle bn)
        {
            if (OMEGA_NORMAL == an) {
                if (0 == EDGETYPE(av)) return true;
            }
            if (OMEGA_INFINITY == bn) {
                return true;
            }
            return false;
        }
        inline static bool simplifiesToSecondArg(int L,
                const forest* f1, edge_value &av, node_handle &an,
                const forest* f2, const edge_value &bv, node_handle bn)
        {
            return false;
        }

        inline static void apply(const edge_value &av, node_handle an,
                const edge_value &bv, node_handle bn,
                edge_value &cv, node_handle &cn)
        {
            //
            // Special case: modding by infinity
            //
            if (OMEGA_INFINITY == bn)
            {
                if (OMEGA_INFINITY == an) {
                    throw error(error::INFINITY_DIV_INFINITY, __FILE__, __LINE__);
                }
                cn = OMEGA_NORMAL;
                cv = av;
                return;
            }

            //
            // Special case: modding by zero
            //
            MEDDLY_DCASSERT(OMEGA_NORMAL == bn);
            if (0 == EDGETYPE(bv)) {
                throw error(error::DIVIDE_BY_ZERO, __FILE__, __LINE__);
            }

            //
            // Still going? We're modding by something finite.
            //

            if (OMEGA_INFINITY == an)
            {
                cn = OMEGA_INFINITY;
                cv = EDGETYPE(0);
                return;
            }

            MEDDLY_DCASSERT(OMEGA_NORMAL == an);

            cn = OMEGA_NORMAL;
            cv.set(EDGETYPE(av) % EDGETYPE(bv));
        }

    };
};

// ******************************************************************
// *                                                                *
// *                       MOD_factory  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::MOD_factory : public binary_factory {
    public:
        virtual void setup();
        virtual binary_operation* build_new(forest* a, forest* b, forest* c);
};

// ******************************************************************

void MEDDLY::MOD_factory::setup()
{
    _setup(__FILE__, "MODULO", "Integer modulo. Behavior is undefined for division by zero. Forest ranges must be integer. Forests should be all MT, EV+, or EV*, over the same domain.");
}

MEDDLY::binary_operation*
MEDDLY::MOD_factory::build_new(forest* a, forest* b, forest* c)
{
    if (c->isMultiTerminal()) {
        if (a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL )
        {
            return nullptr;
        }
        return new arith_compat<EdgeOp_none, mt_mod<long> > (a,b,c);
    }

    if (c->isEVPlus()) {
        switch (c->getEdgeType()) {
            case edge_type::INT:
                return new arith_pushdn<EdgeOp_plus<int>, evplus_mod<int> >
                        (a,b,c);

            case edge_type::LONG:
                return new arith_pushdn<EdgeOp_plus<long>, evplus_mod<long> >
                        (a,b,c);

            default:
                return nullptr;
        }
    }

    return nullptr;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_factory& MEDDLY::MODULO()
{
    static MOD_factory F;
    return F;
}


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
#include "arith_div.h"

#include "../ops_builtin.h" // for COPY
#include "../oper_binary.h"
#include "../oper_unary.h"
#include "../ct_vector.h"
#include "../forest_levels.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

#include "arith_templ.h"

namespace MEDDLY {
    class DIV_factory;
};


// ******************************************************************
// *                                                                *
// *                          mt_div class                          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct mt_div {
        inline static const char* name() {
            return "div";
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
            //   otherwise we're returning 1 for 0/0
            //
            cv.set();
            c = fc->makeRedundantsTo( fc->handleForValue(RANGE(1)), 0, L );
        }
        inline static bool simplifiesToFirstArg(int L,
                const forest* fa, node_handle &a,
                const forest* fb, node_handle b)
        {
            if (0 == a) return true;
            if (fb->isIdentityReduced()) return false;
            terminal one( RANGE(1) );
            return (one.getHandle() == b);
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
            c = fc->handleForValue( av / bv );
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
// *                        evplus_div class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EDGETYPE>
    struct evplus_div {
        inline static const char* name() {
            return "div";
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
            //   otherwise we're returning 1 for 0/0
            //
            cv = EDGETYPE(1);
            c = fc->makeRedundantsTo( OMEGA_NORMAL, 0, L );
        }
        inline static bool simplifiesToFirstArg(int L,
                const forest* f1, edge_value &av, node_handle &an,
                const forest* f2, const edge_value &bv, node_handle bn)
        {
            if (OMEGA_NORMAL == an) {
                if (0 == EDGETYPE(av)) return true;
            }
            if (OMEGA_NORMAL == bn) {
                if (1 == EDGETYPE(bv)) return (!f2->isIdentityReduced());
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
            // Special case: dividing by infinity
            //
            if (OMEGA_INFINITY == bn)
            {
                if (OMEGA_INFINITY == an) {
                    throw error(error::INFINITY_DIV_INFINITY, __FILE__, __LINE__);
                }
                cn = OMEGA_NORMAL;
                cv = EDGETYPE(0);
                return;
            }

            //
            // Special case: dividing by zero
            //
            MEDDLY_DCASSERT(OMEGA_NORMAL == bn);
            if (0 == EDGETYPE(bv)) {
                throw error(error::DIVIDE_BY_ZERO, __FILE__, __LINE__);
            }

            //
            // Still going? We're dividing by something finite.
            //

            if (OMEGA_INFINITY == an)
            {
                cn = OMEGA_INFINITY;
                cv = EDGETYPE(0);
                return;
            }

            MEDDLY_DCASSERT(OMEGA_NORMAL == an);
            cn = OMEGA_NORMAL;
            cv.set(EDGETYPE(av) / EDGETYPE(bv));
        }

    };
};

// ******************************************************************
// *                                                                *
// *                        evstar_div class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EDGETYPE>
    struct evstar_div {
        inline static const char* name() {
            return "div";
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
            //   otherwise we're returning 1 for 0/0
            //
            cv = EDGETYPE(1);
            c = fc->makeRedundantsTo( OMEGA_NORMAL, 0, L );
        }
        inline static bool simplifiesToFirstArg(int L,
                const forest* fa, node_handle &a,
                const forest* fb, node_handle b)
        {
            return (OMEGA_ZERO == a);
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
            if (OMEGA_ZERO == b) {
                throw error(error::DIVIDE_BY_ZERO, __FILE__, __LINE__);
            }
            if (OMEGA_ZERO == a)
            {
                c = OMEGA_ZERO;
            }
            else
            {
                c = OMEGA_NORMAL;
            }
        }

        inline static void apply(const edge_value &a, const edge_value &b,
                                    edge_value &c)
        {
            if (0 == EDGETYPE(b)) {
                throw error(error::DIVIDE_BY_ZERO, __FILE__, __LINE__);
            }
            c.set(EDGETYPE(a) / EDGETYPE(b));
        }

    };
};

// ******************************************************************
// *                                                                *
// *                       DIV_factory  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::DIV_factory : public binary_factory {
    public:
        virtual void setup();
        virtual binary_operation* build_new(forest* a, forest* b, forest* c);
};

// ******************************************************************

void MEDDLY::DIV_factory::setup()
{
    _setup(__FILE__, "DIVIDE", "Division. Behavior is undefined for division by zero. Forest ranges must be integer or real. Forests should be all MT, EV+, or EV*, over the same domain.");
}

MEDDLY::binary_operation*
MEDDLY::DIV_factory::build_new(forest* a, forest* b, forest* c)
{
    if (c->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            return new arith_compat<EdgeOp_none, mt_div<float> > (a,b,c);
        } else {
            return new arith_compat<EdgeOp_none, mt_div<long> > (a,b,c);
        }
    }

    if (c->isEVPlus()) {
        switch (c->getEdgeType()) {
            case edge_type::INT:
                return new arith_pushdn<EdgeOp_plus<int>, evplus_div<int> >
                        (a,b,c);

            case edge_type::LONG:
                return new arith_pushdn<EdgeOp_plus<long>, evplus_div<long> >
                        (a,b,c);

            default:
                return nullptr;
        }
    }

    if (c->isEVTimes()) {
        switch (c->getEdgeType()) {
            case edge_type::FLOAT:
                return new arith_compat<EdgeOp_times<float>, evstar_div<float> >
                        (a,b,c);

            case edge_type::DOUBLE:
                return new arith_compat<EdgeOp_times<double>, evstar_div<double> >
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

MEDDLY::binary_factory& MEDDLY::DIVIDE()
{
    static DIV_factory F;
    return F;
}


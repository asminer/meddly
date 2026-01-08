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
#include "arith_mult.h"

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
    class MULT_factory;
};


// ******************************************************************
// *                                                                *
// *                         mt_mult  class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct mt_mult {
        inline static const char* name() {
            return "mult";
        }
        inline static bool commutes() {
            return true;
        }
        inline static bool stopOnEqualArgs() {
            return false;
        }
        inline static void makeEqualResult(int L, unsigned in,
                const forest* fa, node_handle a,
                forest* fc, edge_value &cv, node_handle &c,
                unary_operation* copier)
        {
            FAIL(__FILE__, __LINE__);
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
            if (0 == b) return true;
            if (fa->isIdentityReduced()) return false;
            terminal one( RANGE(1) );
            return (one.getHandle() == a);
        }

        inline static void apply(const forest* fa, node_handle a,
                const forest* fb, node_handle b,
                const forest* fc, node_handle &c)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            c = fc->handleForValue( av * bv );
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
// *                       evplus_mult  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EDGETYPE>
    struct evplus_mult {
        inline static const char* name() {
            return "mult";
        }
        inline static bool commutes() {
            return true;
        }
        inline static bool stopOnEqualArgs() {
            return false;
        }
        inline static void makeEqualResult(int L, unsigned in,
                const edge_value &av, const node_handle &a,
                forest* fc, edge_value &cv, node_handle &c,
                unary_operation* copyop)
        {
            FAIL(__FILE__, __LINE__);
        }
        inline static bool simplifiesToFirstArg(int L,
                const forest* f1, edge_value &av, node_handle &an,
                const forest* f2, const edge_value &bv, node_handle bn)
        {
            if (OMEGA_INFINITY == an) return (!f1->isIdentityReduced());
            if (OMEGA_NORMAL == an) {
                EDGETYPE aev;
                av.get(aev);
                if (0 == aev) return true;
            }
            if (OMEGA_NORMAL == bn) {
                EDGETYPE bev;
                bv.get(bev);
                if (1 == bev) return (!f2->isIdentityReduced());
            }
            return false;
        }
        inline static bool simplifiesToSecondArg(int L,
                const forest* f1, edge_value &av, node_handle &an,
                const forest* f2, const edge_value &bv, node_handle bn)
        {
            if (OMEGA_INFINITY == bn) return (!f2->isIdentityReduced());
            if (OMEGA_NORMAL == bn) {
                EDGETYPE bev;
                bv.get(bev);
                if (0 == bev) return true;
            }
            if (OMEGA_NORMAL == an) {
                EDGETYPE aev;
                av.get(aev);
                if (1 == aev) return (!f1->isIdentityReduced());
            }
            return false;
        }

        inline static void apply(const edge_value &av, node_handle an,
                const edge_value &bv, node_handle bn,
                edge_value &cv, node_handle &cn)
        {
            if ((OMEGA_INFINITY == an) || (OMEGA_INFINITY == bn))
            {
                cn = OMEGA_INFINITY;
                cv = EDGETYPE(0);
                return;
            }
            MEDDLY_DCASSERT(OMEGA_NORMAL == an);
            MEDDLY_DCASSERT(OMEGA_NORMAL == bn);
            cn = OMEGA_NORMAL;
            EDGETYPE aev, bev;
            av.get(aev);
            bv.get(bev);
            cv.set(aev * bev);
        }

    };
};

// ******************************************************************
// *                                                                *
// *                       evstar_mult  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EDGETYPE>
    struct evstar_mult {
        inline static const char* name() {
            return "mult";
        }
        inline static bool commutes() {
            return true;
        }
        inline static bool stopOnEqualArgs() {
            return false;
        }
        inline static void makeEqualResult(int L, unsigned in,
                const forest* fa, node_handle a,
                forest* fc, edge_value &cv, node_handle &c,
                unary_operation* copier)
        {
            FAIL(__FILE__, __LINE__);
        }
        inline static bool simplifiesToFirstArg(int L,
                const forest* fa, node_handle &a,
                const forest* fb, node_handle b)
        {
            if (OMEGA_ZERO == a) return true;
            if (fb->isIdentityReduced()) return false;
            return (OMEGA_NORMAL == b);
        }
        inline static bool simplifiesToSecondArg(int L,
                const forest* fa, node_handle a,
                const forest* fb, node_handle &b)
        {
            if (OMEGA_ZERO == b) return true;
            if (fa->isIdentityReduced()) return false;
            return (OMEGA_NORMAL == a);
        }

        inline static void apply(const forest* fa, node_handle a,
                const forest* fb, node_handle b,
                const forest* fc, node_handle &c)
        {
            if ((OMEGA_ZERO == a) || (OMEGA_ZERO == b))
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
            EDGETYPE av, bv;
            a.get(av);
            b.get(bv);
            c.set(av * bv);
        }

    };
};

// ******************************************************************
// *                                                                *
// *                       MULTIPLY front end                       *
// *                                                                *
// ******************************************************************
/*
MEDDLY::binary_operation* MEDDLY::MULTIPLY(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  MULT_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (c->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            bop = new arith_compat<EdgeOp_none, mt_mult<float> > (a,b,c);
        } else {
            bop = new arith_compat<EdgeOp_none, mt_mult<long> > (a,b,c);
        }
        return MULT_cache.add( bop );
    }

    if (c->isEVPlus()) {
        switch (c->getEdgeType()) {
            case edge_type::INT:
                bop = new arith_pushdn<EdgeOp_plus<int>, evplus_mult<int> >
                        (a,b,c);
                break;

            case edge_type::LONG:
                bop = new arith_pushdn<EdgeOp_plus<long>, evplus_mult<long> >
                        (a,b,c);
                break;

            default:
                return nullptr;
        }
        return MULT_cache.add( bop );
    }

    if (c->isEVTimes()) {
        switch (c->getEdgeType()) {
            case edge_type::FLOAT:
                bop = new arith_compat<EdgeOp_times<float>, evstar_mult<float> >
                        (a,b,c);
                break;

            case edge_type::DOUBLE:
                bop = new arith_compat<EdgeOp_times<double>, evstar_mult<double> >
                        (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
        return MULT_cache.add( bop );
    }

    return nullptr;
}

void MEDDLY::MULTIPLY_init()
{
    MULT_cache.reset("Multiply");
}

void MEDDLY::MULTIPLY_done()
{
    MEDDLY_DCASSERT(MULT_cache.isEmpty());
}
*/
// ******************************************************************
// *                                                                *
// *                       MULT_factory class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::MULT_factory : public binary_factory {
    public:
        virtual void setup();
        virtual binary_operation* build_new(forest* a, forest* b, forest* c);
};

// ******************************************************************

void MEDDLY::MULT_factory::setup()
{
    _setup(__FILE__, "MULTIPLY", "Multiplication. Forest ranges must be integer or real. Forests should be all MT, EV+, or EV*, over the same domain.");
}

MEDDLY::binary_operation*
MEDDLY::MULT_factory::build_new(forest* a, forest* b, forest* c)
{
    if (c->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            return new arith_compat<EdgeOp_none, mt_mult<float> > (a,b,c);
        } else {
            return new arith_compat<EdgeOp_none, mt_mult<long> > (a,b,c);
        }
    }

    if (c->isEVPlus()) {
        switch (c->getEdgeType()) {
            case edge_type::INT:
                return new arith_pushdn<EdgeOp_plus<int>, evplus_mult<int> >
                        (a,b,c);

            case edge_type::LONG:
                return new arith_pushdn<EdgeOp_plus<long>, evplus_mult<long> >
                        (a,b,c);

            default:
                return nullptr;
        }
    }

    if (c->isEVTimes()) {
        switch (c->getEdgeType()) {
            case edge_type::FLOAT:
                return new arith_compat<EdgeOp_times<float>, evstar_mult<float> >
                        (a,b,c);

            case edge_type::DOUBLE:
                return new arith_compat<EdgeOp_times<double>, evstar_mult<double> >
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

MEDDLY::binary_factory& MEDDLY::MULTIPLY()
{
    static MULT_factory F;
    return F;
}


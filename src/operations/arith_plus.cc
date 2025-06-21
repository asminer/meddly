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
#include "arith_plus.h"

#include "../ops_builtin.h" // for COPY
#include "../oper_item.h"
#include "../oper_binary.h"
#include "../oper_unary.h"
#include "../ct_vector.h"
#include "../forest_levels.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

#include "arith_templ.h"

//
// Operation instance cache
//
namespace MEDDLY {
    binary_list PLUS_cache;
};


// ******************************************************************
// *                                                                *
// *                         mt_plus  class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct mt_plus {
        inline static const char* name() {
            return "plus";
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
            return (0 == b);
        }
        inline static bool simplifiesToSecondArg(int L,
                const forest* fa, node_handle a,
                const forest* fb, node_handle &b)
        {
            return (0 == a);
        }

        inline static void apply(const forest* fa, node_handle a,
                const forest* fb, node_handle b,
                const forest* fc, node_handle &c)
        {
            RANGE av, bv;
            fa->getValueFromHandle(a, av);
            fb->getValueFromHandle(b, bv);
            c = fc->handleForValue( av + bv );
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
// *                       evplus_plus  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EDGETYPE>
    struct evplus_plus {
        inline static const char* name() {
            return "plus";
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
            if (fa->isIdentityReduced() || fb->isIdentityReduced()) {
                return false;
            }
            return (OMEGA_INFINITY == a) || (OMEGA_NORMAL == b);
        }
        inline static bool simplifiesToSecondArg(int L,
                const forest* fa, node_handle a,
                const forest* fb, node_handle &b)
        {
            if (fa->isIdentityReduced() || fb->isIdentityReduced()) {
                return false;
            }
            return (OMEGA_INFINITY == b) || (OMEGA_NORMAL == a);
        }

        inline static void apply(const forest* fa, node_handle a,
                const forest* fb, node_handle b,
                const forest* fc, node_handle &c)
        {
            if ((OMEGA_INFINITY == a) || (OMEGA_INFINITY == b))
            {
                c = OMEGA_INFINITY;
            }
            else
            {
                MEDDLY_DCASSERT(OMEGA_NORMAL == a);
                MEDDLY_DCASSERT(OMEGA_NORMAL == b);
                c = OMEGA_NORMAL;
            }
        }

        inline static void apply(const edge_value &a, const edge_value &b,
                                    edge_value &c)
        {
            EDGETYPE av, bv;
            a.get(av);
            b.get(bv);
            c.set(av+bv);
        }
    };
};

// ******************************************************************
// *                                                                *
// *                       evstar_plus  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EDGETYPE>
    struct evstar_plus {
        inline static const char* name()
        {
            return "plus";
        }
        inline static bool commutes()
        {
            return true;
        }
        inline static bool stopOnEqualArgs()
        {
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
                const forest*, const edge_value &av, node_handle &a,
                const forest*, const edge_value &bv, node_handle b)
        {
            return (OMEGA_ZERO == b);
        }
        inline static bool simplifiesToSecondArg(int L,
                const forest*, const edge_value &av, node_handle a,
                const forest*, const edge_value &bv, node_handle &b)
        {
            return (OMEGA_ZERO == a);
        }

        inline static void factor(
                edge_value &c, node_handle d,
                edge_value &e, node_handle f,
                edge_value &a)
        {
            if (OMEGA_ZERO == d) {
                //
                // left operand is zero;
                // factor using the right operand.
                //
                if (OMEGA_ZERO == f) {
                    a = edge_value(EDGETYPE(1));
                } else {
                    a = e;
                    e = EDGETYPE(1);
                }
            } else {
                //
                // Factor on left operand
                //
                EDGETYPE av;
                c.get(av);
                c = EDGETYPE(1);
                a = av;
                e.divide(av);
            }
        }

        inline static bool alwaysFactorsToIdentity()
        {
            return true;
        }

        inline static void apply(const edge_value &c, node_handle d,
                          const edge_value &e, node_handle f,
                          edge_value &a, node_handle &b)
        {
            EDGETYPE av, cv, ev;
            c.get(cv);
            e.get(ev);
            av = cv + ev;
            a = av;
            b = av ? OMEGA_NORMAL : OMEGA_ZERO;
        }

    };
};

// ******************************************************************
// *                                                                *
// *                         PLUS front end                         *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::PLUS(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  PLUS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (c->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            bop = new arith_compat<EdgeOp_none, mt_plus<float> > (a,b,c);
        } else {
            bop = new arith_compat<EdgeOp_none, mt_plus<long> > (a,b,c);
        }
        return PLUS_cache.add( bop );
    }

    if (c->isEVPlus()) {
        switch (c->getEdgeType()) {
            case edge_type::INT:
                bop = new arith_compat<EdgeOp_plus<int>, evplus_plus<int> >
                        (a,b,c);
                break;

            case edge_type::LONG:
                bop = new arith_compat<EdgeOp_plus<long>, evplus_plus<long> >
                        (a,b,c);
                break;

            default:
                return nullptr;
        }
        return PLUS_cache.add( bop );
    }

    if (c->isEVTimes()) {
        switch (c->getEdgeType()) {
            case edge_type::FLOAT:
                bop = new arith_factor<EdgeOp_times<float>, evstar_plus<float> >
                        (a,b,c);
                break;

            case edge_type::DOUBLE:
                bop = new arith_factor<EdgeOp_times<double>, evstar_plus<double> >
                        (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
        return PLUS_cache.add( bop );
    }

    return nullptr;
}

void MEDDLY::PLUS_init()
{
    PLUS_cache.reset("Plus");
}

void MEDDLY::PLUS_done()
{
    MEDDLY_DCASSERT(PLUS_cache.isEmpty());
}


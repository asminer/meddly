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
#include "arith_minus.h"

#include "../ops_builtin.h" // for COPY
#include "../oper_item.h"
#include "../oper_binary.h"
#include "../oper_unary.h"
#include "../ct_vector.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

#include "arith_templ.h"

//
// Operation instance cache
//
namespace MEDDLY {
    binary_list MINUS_cache;
};


// ******************************************************************
// *                                                                *
// *                         mt_minus class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RANGE>
    struct mt_minus {
        inline static const char* name() {
            return "minus";
        }
        inline static bool commutes() {
            return false;
        }
        inline static bool stopOnEqualArgs() {
            return true;
        }
        inline static void makeEqualResult(const forest* F, node_handle &a) {
            a = F->handleForValue(RANGE(0));
        }
        inline static bool simplifiesToFirstArg(
                const forest* fa, node_handle &a,
                const forest* fb, node_handle b)
        {
            return (0 == b);
        }
        inline static bool simplifiesToSecondArg(
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
            c = fc->handleForValue( av - bv );
        }
    };
};

// ******************************************************************
// *                                                                *
// *                       evplus_minus class                       *
// *                                                                *
// ******************************************************************

#if 0

namespace MEDDLY {
    template <class EDGETYPE>
    struct evplus_minus {
        inline static const char* name()
        {
            return "minus";
        }
        inline static bool commutes()
        {
            return false;
        }
        inline static bool stopOnEqualArgs()
        {
            return true;
        }
        inline static void makeEqualResult(edge_value &av, node_handle &a)
        {
            // minus(a, a) = a; do nothing to a
        }
        inline static bool simplifiesToFirstArg(
                const edge_value &av, node_handle &a,
                const edge_value &bv, node_handle b)
        {
            return (OMEGA_INFINITY == a);
        }
        inline static bool simplifiesToSecondArg(
                const edge_value &av, node_handle a,
                const edge_value &bv, node_handle &b)
        {
            return (OMEGA_INFINITY == b);
        }

        inline static void factor(
                edge_value &c, node_handle d,
                edge_value &e, node_handle f,
                edge_value &a)
        {
            if (OMEGA_INFINITY == d) {
                //
                // left operand is infinity;
                // factor using the right operand.
                //
                if (OMEGA_INFINITY == f) {
                    a = edge_value(EDGETYPE(0));
                } else {
                    a = e;
                    e = edge_value(EDGETYPE(0));
                }
            } else {
                //
                // Factor on left operand
                //
                a = c;
                c = edge_value(EDGETYPE(0));
                EDGETYPE av;
                a.get(av);
                e.subtract(av);
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
            if ((OMEGA_INFINITY == d) || (OMEGA_INFINITY == f))
            {
                a = edge_value(EDGETYPE(0));
                b = OMEGA_INFINITY;
            } else {
                //
                // both finite
                //
                b = OMEGA_NORMAL;
                EDGETYPE cv, ev;
                c.get(cv);
                e.get(ev);
                a = edge_value(MAX(cv, ev));
            }
        }

    };
};

#endif

// ******************************************************************
// *                                                                *
// *                       evstar_minus class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EDGETYPE>
    struct evstar_minus {
        inline static const char* name()
        {
            return "minus";
        }
        inline static bool commutes()
        {
            return false;
        }
        inline static bool stopOnEqualArgs()
        {
            return true;
        }
        inline static void makeEqualResult(edge_value &av, node_handle &a)
        {
            av = EDGETYPE(0);
            a = OMEGA_ZERO;
        }
        inline static bool simplifiesToFirstArg(
                const edge_value &av, node_handle &a,
                const edge_value &bv, node_handle b)
        {
            return (OMEGA_ZERO == b);
        }
        inline static bool simplifiesToSecondArg(
                const edge_value &av, node_handle a,
                const edge_value &bv, node_handle &b)
        {
            return false;
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
            av = cv - ev;
            a = av;
            b = av ? OMEGA_NORMAL : OMEGA_ZERO;
        }

    };
};


// ******************************************************************
// *                                                                *
// *                        MINUS  front end                        *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::MINUS(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  MINUS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (c->isMultiTerminal()) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            bop = new arith_compat<EdgeOp_none, mt_minus<float> > (a,b,c);
        } else {
            bop = new arith_compat<EdgeOp_none, mt_minus<long> > (a,b,c);
        }
        return MINUS_cache.add( bop );
    }

    /*

    if (c->isEVPlus()) {
        switch (c->getEdgeType()) {
            case edge_type::INT:
                bop = new arith_factor<EdgeOp_plus<int>, evplus_minus<int> >
                        (a,b,c);
                break;

            case edge_type::LONG:
                bop = new arith_factor<EdgeOp_plus<long>, evplus_minus<long> >
                        (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
        return MINUS_cache.add( bop );
    }

    */

    if (c->isEVTimes()) {
        switch (c->getEdgeType()) {
            case edge_type::FLOAT:
                bop = new arith_factor<EdgeOp_times<float>, evstar_minus<float> >
                        (a,b,c);
                break;

            case edge_type::DOUBLE:
                bop = new arith_factor<EdgeOp_times<double>, evstar_minus<double> >
                        (a,b,c);
                break;

            default:
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
        return MINUS_cache.add( bop );
    }

    return nullptr;
}

void MEDDLY::MINUS_init()
{
    MINUS_cache.reset("Minus");
}

void MEDDLY::MINUS_done()
{
    MEDDLY_DCASSERT(MINUS_cache.isEmpty());
}


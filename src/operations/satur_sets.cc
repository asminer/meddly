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
#include "satur_sets.h"
#include "satur_index.h"

// #define RECFIRE_THEN_SAT

// #define TRACE
// #define DEBUG_SPLIT
// #define DEBUG_SPLIT_FULL
// #define TRACE_RECFIRE

// #define COUNT_CALLS

#include "satur_sets_v1.h"
#include "satur_sets_v2.h"


// ******************************************************************
// *                                                                *
// *                  reachset_satur_factory class                  *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <bool FWD, int VER>
    class reachset_satur_factory : public binary_factory {
        public:
            virtual void setup();

            virtual binary_operation*
                build_new(forest* a, forest* b, forest* c);

            template <class EOP, class ATYPE>
            inline static binary_operation* new_sat_set_mt(bool fwd,
                    forest* a, forest* b, forest* c)
            {
                switch (VER) {
                    case 1:
                        return new saturation1_set_mtrel<EOP, ATYPE>
                            (FWD, a, b, c);

                    case 2:
                        return new saturation2_set_mtrel<EOP, ATYPE>
                            (FWD, a, b, c);
                    default:
                        MEDDLY_DCASSERT(false);
                        return nullptr;
                }
            }
    };
};

// ******************************************************************

template <bool FWD, int VER>
void MEDDLY::reachset_satur_factory <FWD, VER>::setup()
{
    switch (VER) {
        case 1:
            if (FWD) {
                _setup(__FILE__, "REACHABLE_SATUR(true, 1)", "Build forward reachability set using saturation version 1 (the algorithm of Ciardo, Lüttgen, and Siminiceanu, 2001). The first argument is the set of initial states, and the second argument is the transition relation.");
            } else {
                _setup(__FILE__, "REACHABLE_SATUR(false, 1)", "Build backward reachability set using saturation version 1 (the algorithm of Ciardo, Lüttgen, and Siminiceanu, 2001). The first argument is the set of initial states, and the second argument is the transition relation.");
            }
            return;

        case 2:
            if (FWD) {
                _setup(__FILE__, "REACHABLE_SATUR(true, 2)", "Build forward reachability set using saturation version 2 (the algorithm of Molnár and Majzik, 2019). The first argument is the set of initial states, and the second argument is the transition relation.");
            } else {
                _setup(__FILE__, "REACHABLE_SATUR(false, 2)", "Build backward reachability set using saturation version 2 (the algorithm of Molnár and Majzik, 2019). The first argument is the set of initial states, and the second argument is the transition relation.");
            }
            return;

        default:
            MEDDLY_DCASSERT(false);
            if (FWD) {
                _setup(__FILE__, "REACHABLE_SATUR(true, ?)", "Unsupported version");
            } else {
                _setup(__FILE__, "REACHABLE_SATUR(false, ?)", "Unsupported version");
            }
    }
}

template <bool FWD, int VER>
MEDDLY::binary_operation*
MEDDLY::reachset_satur_factory <FWD, VER>::build_new(forest* a, forest* b, forest* c)
{
    MEDDLY_DCASSERT(1==VER);

    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {

        switch (c->getRangeType()) {
            case range_type::BOOLEAN:
                return new_sat_set_mt<EdgeOp_none, mt_prepost>
                    (FWD, a, b, c);

            case range_type::INTEGER:
                if (c->isFullyReduced())  {
                    return new_sat_set_mt<EdgeOp_none, mt_distance>
                        (FWD, a, b, c);
                }

            default:
                return nullptr;
        }
    }

    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {

        switch (a->getEdgeType()) {
            case edge_type::INT:
                return new_sat_set_mt<EdgeOp_plus<int>, ev_prepost<int> >
                    (FWD, a, b, c);

            case edge_type::LONG:
                return new_sat_set_mt<EdgeOp_plus<long>, ev_prepost<long> >
                    (FWD, a, b, c);

            default:
                return nullptr;
        };
    }
    return nullptr;
}




// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_factory& MEDDLY::REACHABLE_SATUR(bool fwd, int version)
{
    //
    // Version 1
    //
    static reachset_satur_factory<true, 1>  forwd1;
    static reachset_satur_factory<false, 1> bckwd1;

    //
    // Version 2
    //

    static reachset_satur_factory<true, 2>  forwd2;
    static reachset_satur_factory<false, 2> bckwd2;

    switch (version) {
        case 1:
            if (fwd) return forwd1;
            else     return bckwd1;

        case 2:
            if (fwd) return forwd2;
            else     return bckwd2;

        default:
            return BOGUS_BINARY();
    }
}


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

#ifndef PREPOST_COMMON_H
#define PREPOST_COMMON_H

//
// Structs for use with preimage, postimage, and saturation operations.
//

// ******************************************************************
// *                                                                *
// *                       mt_prepost  struct                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    struct mt_prepost {
        inline static const char* name(bool forwd)
        {
            return forwd ? "post-image" : "pre-image";
        }

        /// Get the accumulate operation for result nodes.
        inline static binary_operation* accumulateOp(forest* resF)
        {
            return build(UNION, resF, resF, resF);
        }

        /// Does the edge correspond to an unreachable state?
        inline static bool isUnreachable(const edge_value &ev, node_handle p)
        {
            MEDDLY_DCASSERT(ev.isVoid());
            return 0==p;
        }

        /// Does the edge correspond to everything is reachable?
        inline static bool areAllReachable(const edge_value &ev, node_handle p)
        {
            MEDDLY_DCASSERT(ev.isVoid());
            return -1==p;
        }

        /// Set edge to be unreachable states
        inline static void setUnreachable(edge_value &v, node_handle &p)
        {
            v.set();
            p = 0;
        }

        /// Set all edges to unreachable
        inline static void setAllUnreachable(unpacked_node *U)
        {
            U->clear(0, U->getSize());
        }

        /// Apply the operation when b is a terminal node.
        static void apply(forest* fa, const edge_value &av, node_handle a,
                          const forest* fb, node_handle b,
                          forest* fc, edge_value &cv, node_handle &c)
        {
            MEDDLY_DCASSERT(b<0);
            unary_operation* copy = build(COPY, fa, fc);
            MEDDLY_DCASSERT(copy);
            copy->compute(fa->getNodeLevel(a), ~0, av, a, cv, c);
        }

    };
};

// ******************************************************************
// *                                                                *
// *                       mt_distance struct                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    struct mt_distance {
        inline static const char* name(bool forwd)
        {
            return forwd ? "post-image" : "pre-image";
        }

        /// Get the accumulate operation for result nodes.
        inline static binary_operation* accumulateOp(forest* resF)
        {
            return build(DIST_MIN, resF, resF, resF);
        }

        /// Does the edge correspond to an unreachable state?
        inline static bool isUnreachable(const edge_value &ev, node_handle p)
        {
            MEDDLY_DCASSERT(ev.isVoid());
            if (p >= 0) return false;
            terminal t(terminal_type::INTEGER, p);
            return t.getInteger() < 0;
        }

        /// Does the edge correspond to everything is reachable?
        inline static bool areAllReachable(const edge_value &ev, node_handle p)
        {
            MEDDLY_DCASSERT(ev.isVoid());
            return 0==p;
        }

        /// Set edge to be unreachable states
        inline static void setUnreachable(edge_value &v, node_handle &p)
        {
            v.set();
            terminal t = -1;
            p = t.getIntegerHandle();
        }

        /// Set all edges to unreachable
        inline static void setAllUnreachable(unpacked_node *U)
        {
            terminal t = -1;
            node_handle neg1 = t.getIntegerHandle();
            for (unsigned i=0; i<U->getSize(); i++) {
                U->setFull(i, neg1);
            }
        }

        /// Apply the operation when b is a terminal node.
        static void apply(forest* fa, const edge_value &av, node_handle a,
                          const forest* fb, node_handle b,
                          forest* fc, edge_value &cv, node_handle &c)
        {
            MEDDLY_DCASSERT(b<0);
            unary_operation* increment = build(DIST_INC, fa, fc);
            MEDDLY_DCASSERT(increment);
            increment->compute(fa->getNodeLevel(a), ~0, av, a, cv, c);
        }

    };
};


// ******************************************************************
// *                                                                *
// *                       ev_prepost  struct                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <typename INT>
    struct ev_prepost {
        inline static const char* name(bool forwd)
        {
            return forwd ? "post-image" : "pre-image";
        }

        /// Get the accumulate operation for result nodes.
        inline static binary_operation* accumulateOp(forest* resF)
        {
            return build(MINIMUM, resF, resF, resF);
        }

        /// Does the edge correspond to an unreachable state?
        inline static bool isUnreachable(const edge_value &ev, node_handle p)
        {
            return OMEGA_INFINITY == p;
        }

        /// Does the edge correspond to everything is reachable?
        inline static bool areAllReachable(const edge_value &ev, node_handle p)
        {
            if (OMEGA_NORMAL != p) return false;
            return INT(ev) == 0;
        }

        /// Set edge to be unreachable states
        inline static void setUnreachable(edge_value &v, node_handle &p)
        {
            p = OMEGA_INFINITY;
            v = INT(0);
        }

        /// Set all edges to unreachable
        inline static void setAllUnreachable(unpacked_node *U)
        {
            U->clear(0, U->getSize());
        }

        /// Apply the operation when b is a terminal node.
        static void apply(forest* fa, const edge_value &av, node_handle a,
                          const forest* fb, node_handle b,
                          forest* fc, edge_value &cv, node_handle &c)
        {
            MEDDLY_DCASSERT(b<0);
            unary_operation* copy = build(COPY, fa, fc);
            MEDDLY_DCASSERT(copy);
            copy->compute(fa->getNodeLevel(a), ~0, av, a, cv, c);

            if (OMEGA_INFINITY != c) {
                cv.add(INT(1));
            }
        }

    };
};

#endif

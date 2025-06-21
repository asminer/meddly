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

#ifndef MEDDLY_FOREST_EDGERULES_H
#define MEDDLY_FOREST_EDGERULES_H

#include "defines.h"

// ******************************************************************
// *                                                                *
// *                       EdgeOp_none  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    /**
        Small class for MT behavior on edges.
    */
    class EdgeOp_none {
        public:
            /// Do we have edge values
            static inline bool hasEdgeValues()
            {
                return false;
            }
            /// Edge value type; for CT entries
            static inline char edgeValueTypeLetter()
            {
                return 0;
            }
            /// Clear an edge value.
            static inline void clear(edge_value &v)
            {
                v.set();
            }
            /// Return a + b
            static inline edge_value applyOp(const edge_value &a,
                    const edge_value &b)
            {
                return edge_value();
            }
            /// Accumulate a += b
            static inline void accumulateOp(edge_value &a,
                    const edge_value &b)
            {
            }
            /// Accumulate a += -b
            static inline void accumulateInverse(edge_value &a,
                    const edge_value &b)
            {
            }
            static inline terminal buildTerm(
                            const edge_value &val, node_handle p)
            {
                FAIL(__FILE__, __LINE__, "shouldn't be called");
                // Keep compiler happy:
                return terminal(p, terminal_type::OMEGA);
            }
            static inline bool isZeroFunction(const edge_value &val,
                    node_handle p)
            {
                return (0==p);
            }

    };
};


// ******************************************************************
// *                                                                *
// *                       EdgeOp_plus  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    /**
        Small template class for EV+ behavior on edges.
    */
    template <class TYPE>
    class EdgeOp_plus {
        public:
            /// Do we have edge values
            static inline bool hasEdgeValues()
            {
                return true;
            }
            /// Edge value type; for CT entries
            static char edgeValueTypeLetter();
            /// Clear an edge value.
            static inline void clear(edge_value &v)
            {
                v.set(TYPE(0));
            }
            /// Accumulate two edge values. Returns a + b
            static inline edge_value applyOp(const edge_value &a,
                    const edge_value &b)
            {
                TYPE av, bv;
                a.get(av);
                b.get(bv);
                return edge_value(av+bv);
            }
            /// Accumulate one edge value into another: a += b
            static inline void accumulateOp(edge_value &a,
                    const edge_value &b)
            {
                TYPE bv;
                b.get(bv);
                a.add(bv);
            }
            /// Accumulate an inverse: a += -b
            static inline void accumulateInverse(edge_value &a,
                    const edge_value &b)
            {
                TYPE bv;
                b.get(bv);
                a.subtract(bv);
            }
            /// Build terminal value for edge <val, p>
            static inline terminal buildTerm(
                            const edge_value &val, node_handle p)
            {
                if (OMEGA_NORMAL == p) {
                    TYPE v;
                    val.get(v);
                    return terminal(v);
                }
                return terminal(p, terminal_type::OMEGA);
            }
            static inline bool isZeroFunction(const edge_value &val,
                    node_handle p)
            {
                if (p != OMEGA_NORMAL) return false;
                TYPE bv;
                val.get(bv);
                return (0==bv);
            }
    };
    template <>
    inline char EdgeOp_plus<int>::edgeValueTypeLetter() {
        return 'I';
    }
    template <>
    inline char EdgeOp_plus<long>::edgeValueTypeLetter() {
        return 'L';
    }
    template <>
    inline char EdgeOp_plus<float>::edgeValueTypeLetter() {
        return 'F';
    }
    template <>
    inline char EdgeOp_plus<double>::edgeValueTypeLetter() {
        return 'D';
    }
};

// ******************************************************************
// *                                                                *
// *                       EdgeOp_times class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    /**
      Small template class for EV+ behavior on edges.
    */
    template <class TYPE>
    class EdgeOp_times {
        public:
            /// Do we have edge values
            static inline bool hasEdgeValues()
            {
                return true;
            }
            /// Edge value type; for CT entries
            static char edgeValueTypeLetter();
            /// Clear an edge value.
            static inline void clear(edge_value &v)
            {
                v.set(TYPE(1));
            }
            /// Accumulate two edge values. Returns a * b
            static inline edge_value applyOp(const edge_value &a,
                    const edge_value &b)
            {
                TYPE av, bv;
                a.get(av);
                b.get(bv);
                return edge_value(av*bv);
            }
            /// Accumulate one edge value into another: a *= b
            static inline void accumulateOp(edge_value &a,
                    const edge_value &b)
            {
                TYPE bv;
                b.get(bv);
                a.multiply(bv);
            }
            /// Accumulate an inverse: a *= 1/b
            static inline void accumulateInverse(edge_value &a,
                    const edge_value &b)
            {
                TYPE bv;
                b.get(bv);
                a.divide(bv);
            }
            /// Build terminal value for edge <val, p>
            static inline terminal buildTerm(
                            const edge_value &val, node_handle p)
            {
                if (OMEGA_NORMAL == p) {
                    TYPE v;
                    val.get(v);
                    return terminal(v);
                }
                return terminal(p, terminal_type::OMEGA);
            }
            static inline bool isZeroFunction(const edge_value &val,
                    node_handle p)
            {
                return (OMEGA_ZERO == p);
            }
    };
    template <>
    inline char EdgeOp_times<float>::edgeValueTypeLetter() {
        return 'F';
    }
    template <>
    inline char EdgeOp_times<double>::edgeValueTypeLetter() {
        return 'D';
    }
};

#endif // include guard

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
#include "maxmin_range.h"

#include "../oper_item.h"
#include "../oper_unary.h"
#include "../ct_vector.h"
#include "../ct_generics.h"
#include "../forest_levels.h"


namespace MEDDLY {
    class intrange;
    class intmin;
    class intmax;

    class realrange;
    class realmin;
    class realmax;

    class MAXRANGE_factory;
    class MINRANGE_factory;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        Range operations                        *
// *                                                                *
// *                                                                *
// ******************************************************************

/*
    Required methods for RTYPE classes (should be inlined):

        /// Get the operand type
        static opnd_type getOpndType();

        /// Return the compute table letter for the type
        static char getCTletter();

        /// Sets result from a compute table item
        static void set(oper_item &result, const ct_item &cached);

        /// Sets a compute table item from the result
        static void set(ct_item &cached, const oper_item &result);

        /// Initialize an item from a terminal
        static void initItem(oper_item &result, node_handle term);

        /// Update a result from another result
        static void updateItem(oper_item &result, const oper_item &val);

 */

// ******************************************************************
// *                                                                *
// *                       range_templ  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RTYPE>
    class range_templ : public unary_operation {
        public:
            range_templ(forest* arg);
            virtual ~range_templ();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap, oper_item &result);

        protected:
            void _compute(node_handle A, oper_item &result);

        private:
            ct_entry_type* ct;
#ifdef TRACE
            ostream_output out;
#endif
    };
};

// ******************************************************************

template <class RTYPE>
MEDDLY::range_templ<RTYPE>::range_templ(forest* arg)
    : unary_operation(arg, RTYPE::getOpndType())
#ifdef TRACE
      , out(std::cout)
#endif
{
    ct = new ct_entry_type("range");
    ct->setFixed(arg);
    ct->setResult(RTYPE::getCTletter());
    ct->doneBuilding();
}

template <class RTYPE>
MEDDLY::range_templ<RTYPE>::~range_templ()
{
    ct->markForDestroy();
}

template <class RTYPE>
void MEDDLY::range_templ<RTYPE>::compute(int L, unsigned in,
        const edge_value &av, node_handle ap, oper_item &result)
{
    MEDDLY_DCASSERT(result.hasType(RTYPE::getOpndType()));
#ifdef TRACE
    out.indentation(0);
#endif
    _compute(ap, result);
}

template <class RTYPE>
void MEDDLY::range_templ<RTYPE>::_compute(node_handle A, oper_item &r)
{
    //
    // Terminal case
    //
    if (argF->isTerminalNode(A)) {
        RTYPE::initItem(r, A);
        return;
    }

    //
    // Check compute table
    //
    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        RTYPE::set(r, res[0]);
        return;
    }

    //
    // Do computation
    //
    unpacked_node* Au = unpacked_node::newFromNode(argF, A, SPARSE_ONLY);
    _compute(Au->down(0), r);
    oper_item tmp(RTYPE::getOpndType());
    for (unsigned i=1; i<Au->getSize(); i++) {
        _compute(Au->down(i), tmp);
        RTYPE::updateItem(r, tmp);
    }

    //
    // Cleanup
    //
    unpacked_node::Recycle(Au);

    //
    // Save result in CT
    //
    RTYPE::set(res[0], r);
    ct->addCT(key, res);
}

// ******************************************************************
// *                                                                *
// *                     integer helper classes                     *
// *                                                                *
// ******************************************************************

class MEDDLY::intrange {
    public:
        static inline opnd_type getOpndType() {
            return opnd_type::INTEGER;
        }

        static inline char getCTletter() {
            return 'L';
        }

        static inline void set(oper_item &result, const ct_item &cached) {
            result.integer() = cached.getL();
        }

        static inline void set(ct_item &cached, const oper_item &result) {
            cached.setL(result.getInteger());
        }

        static inline void initItem(oper_item &result, node_handle term) {
            terminal t;
            t.setFromHandle(terminal_type::INTEGER, term);
            result.integer() = t.getInteger();
        }
};

class MEDDLY::intmin : public MEDDLY::intrange {
    public:
        static inline void updateItem(oper_item &result, const oper_item &val)
        {
            result.integer() = MIN(result.integer(), val.getInteger());
        }
};

class MEDDLY::intmax : public MEDDLY::intrange {
    public:
        static inline void updateItem(oper_item &result, const oper_item &val)
        {
            result.integer() = MAX(result.integer(), val.getInteger());
        }
};

// ******************************************************************
// *                                                                *
// *                      real helper  classes                      *
// *                                                                *
// ******************************************************************

class MEDDLY::realrange {
    public:
        static inline opnd_type getOpndType() {
            return opnd_type::REAL;
        }

        static inline char getCTletter() {
            return 'D';
        }

        static inline void set(oper_item &result, const ct_item &cached) {
            result.real() = cached.getD();
        }

        static inline void set(ct_item &cached, const oper_item &result) {
            cached.setD(result.getReal());
        }

        static inline void initItem(oper_item &result, node_handle term) {
            terminal t;
            t.setFromHandle(terminal_type::REAL, term);
            result.real() = t.getReal();
        }
};

class MEDDLY::realmin : public MEDDLY::realrange {
    public:
        static inline void updateItem(oper_item &result, const oper_item &val)
        {
            result.real() = MIN(result.real(), val.getReal());
        }
};

class MEDDLY::realmax : public MEDDLY::realrange {
    public:
        static inline void updateItem(oper_item &result, const oper_item &val)
        {
            result.real() = MAX(result.real(), val.getReal());
        }
};

// ******************************************************************
// *                                                                *
// *                     MAXRANGE_factory class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::MAXRANGE_factory : public unary_factory {
    public:
        virtual void setup();
        virtual unary_operation* build(forest* arg, opnd_type res);
};

// ******************************************************************

void MEDDLY::MAXRANGE_factory::setup()
{
    _setup(__FILE__, "MAX_RANGE", "Find the maximum value returned by a function. The result type should match the input forest range type.");
}

MEDDLY::unary_operation*
MEDDLY::MAXRANGE_factory::build(forest* arg, opnd_type res)
{
    if (!arg) return nullptr;
    unary_operation* uop = cache_find(arg, res);
    if (uop) {
        return uop;
    }

    if (arg->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
        return nullptr;

    switch (res) {
        case opnd_type::INTEGER:
            if (range_type::INTEGER != arg->getRangeType())
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            return cache_add(new range_templ<intmax>(arg));

        case opnd_type::REAL:
            if (range_type::REAL != arg->getRangeType())
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            return cache_add(new range_templ<realmax>(arg));

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    } // switch

    return nullptr;
}

// ******************************************************************
// *                                                                *
// *                     MINRANGE_factory class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::MINRANGE_factory : public unary_factory {
    public:
        virtual void setup();
        virtual unary_operation* build(forest* arg, opnd_type res);
};

// ******************************************************************

void MEDDLY::MINRANGE_factory::setup()
{
    _setup(__FILE__, "MIN_RANGE", "Find the minimum value returned by a function. The result type should match the input forest range type.");
}

MEDDLY::unary_operation*
MEDDLY::MINRANGE_factory::build(forest* arg, opnd_type res)
{
    if (!arg) return nullptr;
    unary_operation* uop = cache_find(arg, res);
    if (uop) {
        return uop;
    }

    if (arg->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
        return nullptr;

    switch (res) {
        case opnd_type::INTEGER:
            if (range_type::INTEGER != arg->getRangeType())
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            return cache_add(new range_templ<intmin>(arg));

        case opnd_type::REAL:
            if (range_type::REAL != arg->getRangeType())
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            return cache_add(new range_templ<realmin>(arg));

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    } // switch

    return nullptr;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_factory& MEDDLY::MAX_RANGE()
{
    static MAXRANGE_factory F;
    return F;
}

MEDDLY::unary_factory& MEDDLY::MIN_RANGE()
{
    static MINRANGE_factory F;
    return F;
}


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
#include "cardinality.h"
#ifdef HAVE_LIBGMP
#include <gmp.h>
#endif

#include "../oper_item.h"
#include "../oper_unary.h"
#include "../ct_vector.h"
#include "../ct_generics.h"
#include "../forest_levels.h"

namespace MEDDLY {
    class intcard;
    class realcard;
    class mpzcard;

    class CARD_factory;
};

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        Card  operations                        *
// *                                                                *
// *                                                                *
// ******************************************************************

/*
    Required methods for RTYPE classes (should be inlined):

        /// Get the operand type
        static opnd_type getOpndType();

        /// Return the compute table letter for the type
        static char getCTletter();

        /// Sets result = val
        static void set(oper_item &result, long val);

        /// Sets result from a compute table item
        static void set(oper_item &result, const ct_item &cached);

        /// Sets a compute table item from the result
        static void set(ct_item &cached, const oper_item &result);

        /// Sets result *= scalar
        static void scaleBy(oper_item &result, long scalar);

        /// Sets result += val
        static void addTo(oper_item &result, const oper_item &val);

        /// Initialize a temporary value
        static void initTemp(oper_item &temp);

        /// Done with a temporary value
        static void doneTemp(oper_item &temp);

        /// For debugging
        static void show(output &out, const oper_item &val);

 */

// ******************************************************************
// *                                                                *
// *                        card_templ class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RTYPE>
    class card_templ : public unary_operation {
        public:
            card_templ(forest* arg);
            virtual ~card_templ();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap, oper_item &result);

        protected:
            void _compute(int L, node_handle A, oper_item &result);

        private:
            ct_entry_type* ct;
#ifdef TRACE
            ostream_output out;
#endif
    };
};

// ******************************************************************

template <class RTYPE>
MEDDLY::card_templ<RTYPE>::card_templ(forest* arg)
    : unary_operation(arg, RTYPE::getOpndType())
#ifdef TRACE
      , out(std::cout)
#endif
{
    ct = new ct_entry_type("cardinality");
    ct->setFixed(arg);
    ct->setResult(RTYPE::getCTletter());
    ct->doneBuilding();
}

template <class RTYPE>
MEDDLY::card_templ<RTYPE>::~card_templ()
{
    ct->markForDestroy();
}

template <class RTYPE>
void MEDDLY::card_templ<RTYPE>::compute(int L, unsigned in,
        const edge_value &av, node_handle ap, oper_item &result)
{
    MEDDLY_DCASSERT(result.hasType(RTYPE::getOpndType()));
#ifdef TRACE
    out.indentation(0);
#endif
    _compute(L, ap, result);
}

template <class RTYPE>
void MEDDLY::card_templ<RTYPE>::_compute(int L, node_handle A,
        oper_item &result)
{
    //
    // Terminal cases
    //
    if (0==A) {
        RTYPE::set(result, 0);
        return;
    }
    if (0==L) {
        RTYPE::set(result, 1);
        return;
    }

    const int Ldn = argF->isForRelations()
        ?   MXD_levels::downLevel(L)
        :   MDD_levels::downLevel(L);

    //
    // Quickly deal with skipped levels
    //
    if (argF->getNodeLevel(A) != L) {
        _compute(Ldn, A, result);
        if (L>0 || !argF->isIdentityReduced()) {
            RTYPE::scaleBy(result, argF->getLevelSize(L));
        }
        return;
    }

#ifdef TRACE
    out << "card_templ::_compute(" << L << ", " << A << ")\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        RTYPE::set(result, res[0]);
#ifdef TRACE
        out << "  CT hit: ";
        RTYPE::show(out, result);
        out << '\n';
#endif
        return;
    }

    //
    // Do computation
    //
    unpacked_node* Au = unpacked_node::newFromNode(argF, A, SPARSE_ONLY);

#ifdef TRACE
    out << "    node " << A << ": ";
    Au->show(out, false);
    out.indent_more();
    out.put('\n');
#endif
    RTYPE::set(result, 0);
    oper_item temp;
    RTYPE::initTemp(temp);
    for (unsigned z=0; z<Au->getSize(); z++) {
        _compute(Ldn, Au->down(z), temp);
        RTYPE::addTo(result, temp);
    }
    RTYPE::doneTemp(temp);
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "card_templ::_compute(" << L << ", " << A << ") = ";
    RTYPE::show(out, result);
    out << '\n';
#endif

    //
    // Cleanup
    //
    unpacked_node::Recycle(Au);

    //
    // Save result in CT
    //
    RTYPE::set(res[0], result);
    ct->addCT(key, res);
}


// ******************************************************************
// *                                                                *
// *                   integer cardinality  class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::intcard {
    public:
        static inline opnd_type getOpndType() {
            return opnd_type::INTEGER;
        }

        static inline char getCTletter() {
            return 'L';
        }

        static inline void set(oper_item &result, long val) {
            result.integer() = val;
        }

        static inline void set(oper_item &result, const ct_item &cached) {
            result.integer() = cached.getL();
        }

        static inline void set(ct_item &cached, const oper_item &result) {
            cached.setL(result.getInteger());
        }

        static inline void scaleBy(oper_item &result, long scalar) {
            result.integer() *= scalar;
        }

        static inline void addTo(oper_item &result, const oper_item &val) {
            result.integer() += val.getInteger();
        }

        static inline void initTemp(oper_item& temp) {
            temp.init(0L);
        }

        static inline void doneTemp(oper_item& temp) {
            // nothing to do
        }

        static inline void show(output &out, const oper_item &val) {
            out.put(val.getInteger());
        }
};

// ******************************************************************
// *                                                                *
// *                     real cardinality class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::realcard {
    public:
        static inline opnd_type getOpndType() {
            return opnd_type::REAL;
        }

        static inline char getCTletter() {
            return 'D';
        }

        static inline void set(oper_item &result, long val) {
            result.real() = val;
        }

        static inline void set(oper_item &result, const ct_item &cached) {
            result.real() = cached.getD();
        }

        static inline void set(ct_item &cached, const oper_item &result) {
            cached.setD(result.getReal());
        }

        static inline void scaleBy(oper_item &result, long scalar) {
            result.real() *= scalar;
        }

        static inline void addTo(oper_item &result, const oper_item &val) {
            result.real() += val.getReal();
        }

        static inline void initTemp(oper_item& temp) {
            temp.init(0.0);
        }

        static inline void doneTemp(oper_item& temp) {
            // nothing to do
        }

        static inline void show(output &out, const oper_item &val) {
            out.put(val.getReal());
        }
};


// ******************************************************************
// *                                                                *
// *                     mpz  cardinality class                     *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP

class MEDDLY::mpzcard {
    public:
        static inline opnd_type getOpndType() {
            return opnd_type::HUGEINT;
        }

        static inline char getCTletter() {
            return 'G';
        }

        static inline void set(oper_item &result, long val) {
            mpz_set_si(result.hugeint(), val);
        }

        static inline void set(oper_item &result, const ct_item &cached) {
            cached_mpz* G = smart_cast <cached_mpz*> (cached.getG());
            MEDDLY_DCASSERT(G);
            mpz_set(result.hugeint(), G->mpz());
        }

        static inline void set(ct_item &cached, const oper_item &result) {
            cached_mpz* G = new cached_mpz(result.getHugeint());
            cached.setG(G);
        }

        static inline void scaleBy(oper_item &result, long scalar) {
            mpz_mul_si(result.hugeint(), result.hugeint(), scalar);
        }

        static inline void addTo(oper_item &result, const oper_item &val) {
            mpz_add(result.hugeint(), result.hugeint(), val.getHugeint());
        }

        static inline void initTemp(oper_item& temp) {
            mpz_ptr t = new mpz_t;
            mpz_init(t);
            temp.init(t);
        }

        static inline void doneTemp(oper_item& temp) {
            mpz_ptr t = temp.hugeint();
            mpz_clear(t);
            delete t;
        }

        static inline void show(output &out, const oper_item &val) {
            unsigned digits = mpz_sizeinbase(val.getHugeint(), 10) + 2;
            char* str = new char[digits];
            mpz_get_str(str, 10, val.getHugeint());
            out.put(str);
            delete[] str;
        }
};

#endif

// ******************************************************************
// *                                                                *
// *                       CARD_factory class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::CARD_factory : public unary_factory {
    public:
        virtual void setup();
        virtual unary_operation* build(forest* arg, opnd_type res);
};

// ******************************************************************

void MEDDLY::CARD_factory::setup()
{
    _setup("CARDINALITY", "Set cardinality. Determines the number of variable assignments causing the function to evaluate to non-zero (or non-infinity, for EV+MDDs). The result is allowed to be type long, double, or mpz_t (requires GMP library).");
}

MEDDLY::unary_operation* MEDDLY::CARD_factory::build(forest* arg, opnd_type res)
{
    if (!arg) return nullptr;
    unary_operation* uop =  cache_find(arg, res);
    if (uop) {
        return uop;
    }

    switch (res) {
        case opnd_type::INTEGER:
            return cache_add(new card_templ<intcard>(arg));

        case opnd_type::REAL:
            return cache_add(new card_templ<realcard>(arg));

#ifdef HAVE_LIBGMP
        case opnd_type::HUGEINT:
            return cache_add(new card_templ<mpzcard>(arg));
#endif

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_factory& MEDDLY::CARDINALITY()
{
    static CARD_factory F;
    return F;
}


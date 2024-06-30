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

#include "mpz_object.h"
#include "../ct_entry_key.h"
#include "../ct_entry_result.h"

// #include "../compute_table.h"
#include "../oper_item.h"
#include "../oper_unary.h"
#include "../ct_vector.h"
#include "../ct_generics.h"

// #define DEBUG_CARD

namespace MEDDLY {
    class card_int;
    class card_mdd_int;
    class card_mxd_int;

    class card_real;
    class card_mdd_real;
    class card_mxd_real;

#ifdef HAVE_LIBGMP
    class card_mpz;
    class card_mdd_mpz;
    class card_mxd_mpz;
#endif
    // ^ kill all these

    class intcard;
    class realcard;
    class mpzcard;

    unary_list CARD_cache;
};

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

// #define OLD_CARD


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
        static void set(oper_item &result, int val);

        /// Sets result from a compute table item
        static void set(oper_item &result, const ct_item &cached);

        /// Sets a compute table item from the result
        static void set(ct_item &cached, const oper_item &result);

        /// Sets result *= scalar
        static void scaleBy(oper_item &result, int scalar);

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
// *                    card_mdd  template class                    *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class RTYPE>
    class card_mdd : public unary_operation {
        public:
            card_mdd(forest* arg);
            virtual ~card_mdd();

            virtual void compute(int L, const edge_value &av, node_handle ap,
                    oper_item &result);

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
MEDDLY::card_mdd<RTYPE>::card_mdd(forest* arg)
    : unary_operation(arg, RTYPE::getOpndType())
#ifdef TRACE
      , out(std::cout)
#endif
{
    ct = new ct_entry_type("mdd_cardinality");
    ct->setFixed(arg);
    ct->setResult(RTYPE::getCTletter());
    ct->doneBuilding();
}

template <class RTYPE>
MEDDLY::card_mdd<RTYPE>::~card_mdd()
{
    ct->markForDestroy();
}

template <class RTYPE>
void MEDDLY::card_mdd<RTYPE>::compute(int L, const edge_value &av,
        node_handle ap, oper_item &result)
{
    MEDDLY_DCASSERT(result.hasType(RTYPE::getOpndType()));
#ifdef TRACE
    out.indentation(0);
#endif
    _compute(L, ap, result);
}

template <class RTYPE>
void MEDDLY::card_mdd<RTYPE>::_compute(int L, node_handle A, oper_item &result)
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

    //
    // Quickly deal with skipped levels
    //
    if (argF->getNodeLevel(A) < L) {
        _compute(L-1, A, result);
        RTYPE::scaleBy(result, argF->getLevelSize(L));
        return;
    }

#ifdef TRACE
    out << "card_mdd::_compute(" << L << ", " << A << ")\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(1);
    ct_vector res(1);
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
    unpacked_node* Au = argF->newUnpacked(A, SPARSE_ONLY);

#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    RTYPE::set(result, 0);
    oper_item temp;
    RTYPE::initTemp(temp);
    for (unsigned z=0; z<Au->getSize(); z++) {
        _compute(L-1, Au->down(z), temp);
        RTYPE::addTo(result, temp);
    }
    RTYPE::doneTemp(temp);
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "card_mdd::_compute(" << L << ", " << A << ") = ";
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
// *                    card_mxd  template class                    *
// *                                                                *
// ******************************************************************

// TBD

// ******************************************************************


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

        static inline void set(oper_item &result, int val) {
            result.integer() = val;
        }

        static inline void set(oper_item &result, const ct_item &cached) {
            result.integer() = cached.getL();
        }

        static inline void set(ct_item &cached, const oper_item &result) {
            cached.setL(result.getInteger());
        }

        static inline void scaleBy(oper_item &result, int scalar) {
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

        static inline void set(oper_item &result, int val) {
            result.real() = val;
        }

        static inline void set(oper_item &result, const ct_item &cached) {
            result.real() = cached.getD();
        }

        static inline void set(ct_item &cached, const oper_item &result) {
            cached.setD(result.getReal());
        }

        static inline void scaleBy(oper_item &result, int scalar) {
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

        static inline void set(oper_item &result, int val) {
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

        static inline void scaleBy(oper_item &result, int scalar) {
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
            // TBD
            //
            // out.put(val.getReal());
            out.put("hugeint");
        }
};

#endif

// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// *                                                                *
// *                       OLD IMPLEMENTATION                       *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                         card_int class                         *
// *                                                                *
// ******************************************************************

//  Abstract base class: cardinality that returns an integer
class MEDDLY::card_int : public unary_operation {
    public:
        card_int(forest* arg);

    protected:
        static inline void overflow_acc(long &a, long x) {
            a += x;
            if (a < x) throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
        }
        static inline long overflow_mult(long a, long x) {
            a *= x;
            if (a < x) throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
            return a;
        }
};

MEDDLY::card_int::card_int(forest* arg)
    : unary_operation(CARD_cache, 1, arg, opnd_type::INTEGER)
{
    ct_entry_type* et = new ct_entry_type(CARD_cache.getName(), "N:L");
    et->setForestForSlot(0, arg);
    registerEntryType(0, et);
    buildCTs();
}

// ******************************************************************
// *                                                                *
// *                       card_mdd_int class                       *
// *                                                                *
// ******************************************************************

//  Cardinality on MDDs, returning integer
class MEDDLY::card_mdd_int : public card_int {
    public:
        card_mdd_int(forest* arg) : card_int(arg) { }
        virtual void compute(const dd_edge &arg, long &res) {
            res = compute_r(argF->getMaxLevelIndex(), arg.getNode());
        }
        long compute_r(int k, node_handle a);
};

long MEDDLY::card_mdd_int::compute_r(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0;
  if (0==k) return 1;

  // Quickly deal with skipped levels
  if (argF->getNodeLevel(a) < k) {
    return overflow_mult(compute_r(k-1, a), argF->getLevelSize(k));
  }

  // Check compute table
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    CT0->recycle(CTsrch);
    return CTresult[0].readL();
  }

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // Recurse
  long card = 0;
  int kdn = k-1;
  for (unsigned z=0; z<A->getSize(); z++) {
    overflow_acc(card, compute_r(kdn, A->down(z)));
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeL(card);
  CT0->addEntry(CTsrch, CTresult[0]);

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %ld(L)\n", a, card);
#endif
  return card;
}


// ******************************************************************
// *                                                                *
// *                       card_mxd_int class                       *
// *                                                                *
// ******************************************************************

//  Cardinality on MxDs, returning integer
class MEDDLY::card_mxd_int : public card_int {
    public:
        card_mxd_int(forest* arg) : card_int(arg) { }
        virtual void compute(const dd_edge &arg, long &res) {
            res = compute_r(argF->getMaxLevelIndex(), arg.getNode());
        }
        long compute_r(int k, node_handle a);
};

long MEDDLY::card_mxd_int::compute_r(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0;
  if (0==k) return 1;

  // Quickly deal with skipped levels
  if (isLevelAbove(k, argF->getNodeLevel(a))) {
    if (k<0 && argF->isIdentityReduced()) {
      // identity node
      return compute_r(argF->downLevel(k), a);
    }
    // redundant node
    return overflow_mult(compute_r(argF->downLevel(k), a), argF->getLevelSize(k));
  }

  // Check compute table
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    CT0->recycle(CTsrch);
    return CTresult[0].readL();
  }

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // Recurse
  long card = 0;
  int kdn = argF->downLevel(k);
  for (unsigned z=0; z<A->getSize(); z++) {
    overflow_acc(card, compute_r(kdn, A->down(z)));
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeL(card);
  CT0->addEntry(CTsrch, CTresult[0]);

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %ld(L)\n", a, card);
#endif
  return card;
}


// ******************************************************************
// *                                                                *
// *                        card_real  class                        *
// *                                                                *
// ******************************************************************

//  Abstract base class: cardinality that returns a real
class MEDDLY::card_real : public unary_operation {
    public:
        card_real(forest* arg);
};

MEDDLY::card_real::card_real(forest* arg)
    : unary_operation(CARD_cache, 1, arg, opnd_type::REAL)
{
    ct_entry_type* et = new ct_entry_type(CARD_cache.getName(), "N:D");
    et->setForestForSlot(0, arg);
    registerEntryType(0, et);
    buildCTs();
}

// ******************************************************************
// *                                                                *
// *                      card_mdd_real  class                      *
// *                                                                *
// ******************************************************************

//  Cardinality on MDDs, returning real
class MEDDLY::card_mdd_real : public card_real {
    public:
        card_mdd_real(forest* arg) : card_real(arg) { }
        virtual void compute(const dd_edge &arg, double &res) {
            res = compute_r(argF->getMaxLevelIndex(), arg.getNode());
        }
        double compute_r(int ht, node_handle a);
};

double MEDDLY::card_mdd_real::compute_r(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0.0;
  if (0==k) return 1.0;

  // Quickly deal with skipped levels
  if (argF->getNodeLevel(a) < k) {
    return compute_r(k-1, a) * argF->getLevelSize(k);
  }

  // Check compute table
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    CT0->recycle(CTsrch);
    return CTresult[0].readD();
  }

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // Recurse
  double card = 0;
  int kdn = k-1;
  for (unsigned z=0; z<A->getSize(); z++) {
    card += compute_r(kdn, A->down(z));
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeD(card);
  CT0->addEntry(CTsrch, CTresult[0]);

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %le(L)\n", a, card);
#endif
  return card;
}



// ******************************************************************
// *                                                                *
// *                      card_mxd_real  class                      *
// *                                                                *
// ******************************************************************

//  Cardinality on MxDs, returning real
class MEDDLY::card_mxd_real : public card_real {
    public:
        card_mxd_real(forest* arg) : card_real(arg) { }
        virtual void compute(const dd_edge &arg, double &res) {
            res = compute_r(argF->getMaxLevelIndex(), arg.getNode());
        }
        double compute_r(int k, node_handle a);
};

double MEDDLY::card_mxd_real::compute_r(int k, node_handle a)
{
  // Terminal cases
  if (0==a) return 0.0;
  if (0==k) return 1.0;

  // Quickly deal with skipped levels
  if (isLevelAbove(k, argF->getNodeLevel(a))) {
    if (k<0 && argF->isIdentityReduced()) {
      // identity node
      return compute_r(argF->downLevel(k), a);
    }
    // redundant node
    return compute_r(argF->downLevel(k), a) * argF->getLevelSize(k);
  }

  // Check compute table
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    CT0->recycle(CTsrch);
    return CTresult[0].readD();
  }

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // Recurse
  double card = 0;
  int kdn = argF->downLevel(k);
  for (unsigned z=0; z<A->getSize(); z++) {
    card += compute_r(kdn, A->down(z));
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeD(card);
  CT0->addEntry(CTsrch, CTresult[0]);


#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is %le\n", a, card);
#endif
  return card;
}




// ******************************************************************
// *                                                                *
// *                         card_mpz class                         *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP

//  Abstract base class: cardinality that returns large (mpz) integers.
class MEDDLY::card_mpz : public unary_operation {
    public:
        card_mpz(forest* arg);
};

MEDDLY::card_mpz::card_mpz(forest* arg)
    : unary_operation(CARD_cache, 1, arg, opnd_type::HUGEINT)
{
    ct_entry_type* et = new ct_entry_type(CARD_cache.getName(), "N:G");
    et->setForestForSlot(0, arg);
    registerEntryType(0, et);
    buildCTs();
}

#endif

// ******************************************************************
// *                                                                *
// *                       card_mdd_mpz class                       *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP

/// Cardinality of MDDs, returning large (mpz) integers.
class MEDDLY::card_mdd_mpz : public card_mpz {
    public:
        card_mdd_mpz(forest* arg) : card_mpz(arg) { }
        virtual void compute(const dd_edge& a, ct_object &res) {
            mpz_object& mcard = dynamic_cast <mpz_object &> (res);
            compute_r(argF->getMaxLevelIndex(), a.getNode(), mcard);
        }
        void compute_r(int k, node_handle a, mpz_object &b);
};

void MEDDLY::card_mdd_mpz::compute_r(int k, node_handle a, mpz_object &card)
{
  // Terminal cases
  if (0==a) {
    card.setValue(0);
    return;
  }
  if (0==k) {
    card.setValue(1);
    return;
  }

  // Quickly deal with skipped levels
  if (argF->getNodeLevel(a) < k) {
    // skipped level
    compute_r(k-1, a, card);
    card.multiply(argF->getLevelSize(k));
    return;
  }

  // Check compute table
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    ct_object* G = CTresult[0].readG();
    mpz_object* answer = smart_cast <mpz_object*> (G);
    MEDDLY_DCASSERT(answer);
    answer->copyInto(card);
    CT0->recycle(CTsrch);
    return;
  }

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);
  MEDDLY_DCASSERT(!A->isExtensible());

  // Recurse
  mpz_object tmp;
  tmp.setValue(0);
  card.setValue(0);
  int kdn = k-1;
  for (unsigned z=0; z<A->getSize(); z++) {
    compute_r(kdn, A->down(z), tmp);
    card.add(tmp);
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeG(new mpz_object(card));
  CT0->addEntry(CTsrch, CTresult[0]);

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is ", a);
  card.show(stderr);
  fprintf(stderr, "\n");
#endif
}

#endif

// ******************************************************************
// *                                                                *
// *                       card_mxd_mpz class                       *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP

/// Cardinality of MxDs, returning large (mpz) integers.
class MEDDLY::card_mxd_mpz : public card_mpz {
    public:
        card_mxd_mpz(forest* arg) : card_mpz(arg) { }
        virtual void compute(const dd_edge& a, ct_object &res) {
            mpz_object& mcard = dynamic_cast <mpz_object &> (res);
            compute_r(argF->getMaxLevelIndex(), a.getNode(), mcard);
        }
        void compute_r(int k, node_handle a, mpz_object &b);
};

void MEDDLY::card_mxd_mpz::compute_r(int k, node_handle a, mpz_object &card)
{
  // Terminal cases
  if (0==a) {
    card.setValue(0);
    return;
  }
  if (0==k) {
    card.setValue(1);
    return;
  }

  // Quickly deal with skipped levels
  if (isLevelAbove(k, argF->getNodeLevel(a))) {
    if (k<0 && argF->isIdentityReduced()) {
      // identity node
      compute_r(argF->downLevel(k), a, card);
      return;
    }
    // redundant node
    compute_r(argF->downLevel(k), a, card);
    card.multiply(argF->getLevelSize(k));
    return;
  }

  // Check compute table
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    ct_object* G = CTresult[0].readG();
    mpz_object* answer = smart_cast <mpz_object*> (G);
    MEDDLY_DCASSERT(answer);
    answer->copyInto(card);
    CT0->recycle(CTsrch);
    return;
  }

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // Recurse
  mpz_object tmp;
  tmp.setValue(0);
  card.setValue(0);
  int kdn = argF->downLevel(k);
  for (unsigned z=0; z<A->getSize(); z++) {
    compute_r(kdn, A->down(z), tmp);
    card.add(tmp);
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Add entry to compute table
  CTresult[0].reset();
  CTresult[0].writeG(new mpz_object(card));
  CT0->addEntry(CTsrch, CTresult[0]);

#ifdef DEBUG_CARD
  fprintf(stderr, "Cardinality of node %d is ", a);
  card.show(stderr);
  fprintf(stderr, "\n");
#endif
}

#endif


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_operation* MEDDLY::CARDINALITY(forest* arg, opnd_type res)
{
    if (!arg) return nullptr;
    unary_operation* uop =  CARD_cache.find(arg, res);
    if (uop) {
        return uop;
    }

    switch (res) {
        case opnd_type::INTEGER:
            if (arg->isForRelations())
                return CARD_cache.add(new card_mxd_int(arg));
            else
#ifdef OLD_CARD
                return CARD_cache.add(new card_mdd_int(arg));
#else
                return CARD_cache.add(new card_mdd<intcard>(arg));
#endif

        case opnd_type::REAL:
            if (arg->isForRelations())
                return CARD_cache.add(new card_mxd_real(arg));
            else
#ifdef OLD_CARD
                return CARD_cache.add(new card_mdd_real(arg));
#else
                return CARD_cache.add(new card_mdd<realcard>(arg));
#endif

#ifdef HAVE_LIBGMP
        case opnd_type::HUGEINT:
            if (arg->isForRelations())
                return CARD_cache.add(new card_mxd_mpz(arg));
            else
#ifdef OLD_CARD
                return CARD_cache.add(new card_mdd_mpz(arg));
#else
                return CARD_cache.add(new card_mdd<mpzcard>(arg));
#endif
#endif

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
}

void MEDDLY::CARDINALITY_init()
{
    CARD_cache.reset("Card");
}

void MEDDLY::CARDINALITY_done()
{
    MEDDLY_DCASSERT( CARD_cache.isEmpty() );
}

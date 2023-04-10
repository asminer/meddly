
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


/*! \file meddly.hh

    Implementation details for interface in meddly.h.
*/

#ifndef MEDDLY_HH
#define MEDDLY_HH

#include "defines.h"
#include "error.h"

//---------------------- Inlines ---------------------------------------------

// MEDDLY::

inline MEDDLY::domain* MEDDLY::createDomain() {
  return createDomain((variable**) 0, 0);
}

#ifdef __GMP_H__
inline void MEDDLY::apply(const unary_opname* op, const dd_edge &a, mpz_t &c) {
  ct_object& x = get_mpz_wrapper();
  apply(op, a, opnd_type::HUGEINT, x);
  unwrap(x, c);
}
#endif




// ******************************************************************
// *                                                                *
// *                                                                *
// *                        enumerator class                        *
// *                                                                *
// *                                                                *
// ******************************************************************


// enumerator::iterator::
inline bool MEDDLY::enumerator::iterator::build_error() const {
  return 0==F;
}

inline int MEDDLY::enumerator::iterator::levelChanged() const {
  return level_change;
}

inline const int* MEDDLY::enumerator::iterator::getAssignments() const {
  return index;
}

// enumerator::
inline MEDDLY::enumerator::operator bool() const {
  return is_valid;
}

inline void MEDDLY::enumerator::operator++() {
  is_valid &= I->next();
}

inline const int* MEDDLY::enumerator::getAssignments() const {
  if (I && is_valid) return I->getAssignments(); else return 0;
}

inline const int* MEDDLY::enumerator::getPrimedAssignments() const {
  if (I && is_valid) return I->getPrimedAssignments(); else return 0;
}

inline void MEDDLY::enumerator::getValue(int &v) const {
  if (I && is_valid) I->getValue(v);
}

inline void MEDDLY::enumerator::getValue(long &v) const {
  if (I && is_valid) I->getValue(v);
}

inline void MEDDLY::enumerator::getValue(float &v) const {
  if (I && is_valid) I->getValue(v);
}

inline int MEDDLY::enumerator::levelChanged() const {
  if (I) return I->levelChanged();
  return 0;
}

inline MEDDLY::enumerator::type MEDDLY::enumerator::getType() const {
  return T;
}

#endif

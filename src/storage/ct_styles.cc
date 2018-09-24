
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


#include "ct_styles.h"


// #define DEBUG_SLOW

// #define DEBUG_CT
// #define DEBUG_CT_SEARCHES
// #define DEBUG_CT_SLOTS

// #define DEBUG_REHASH
// #define DEBUG_TABLE2LIST
// #define DEBUG_LIST2TABLE
// #define DEBUG_CTALLOC

// #define DEBUG_ISDEAD
// #define DEBUG_ISSTALE
// #define DEBUG_REMOVESTALES
// #define DEBUG_REMOVESTALES_DETAILS


#define INTEGRATED_MEMMAN

#include "../hash_stream.h"

#ifdef OLD_OP_CT
#include "ct_old.h"
#else
#include "ct_typebased.h"
#include "ct_none.h"
#endif


// **********************************************************************
// *                                                                    *
// *                  monolithic_chained_style methods                  *
// *                                                                    *
// **********************************************************************

MEDDLY::monolithic_chained_style::monolithic_chained_style() 
{ 
}

MEDDLY::compute_table* 
MEDDLY::monolithic_chained_style::create(const ct_initializer::settings &s) const 
{
#ifdef OLD_OP_CT
  return new ct_old<true, true>(s, 0, 0);
#else
  switch (s.compression) {
    case ct_initializer::None:
                                    return new ct_none<true, true>(s, 0, 0);
    case ct_initializer::TypeBased:
                                    return new ct_typebased<true, true>(s, 0, 0);
    default:
                                    return 0;
  }
#endif
}

bool MEDDLY::monolithic_chained_style::usesMonolithic() const 
{
  return true;
}

// **********************************************************************
// *                                                                    *
// *                 monolithic_unchained_style methods                 *
// *                                                                    *
// **********************************************************************


MEDDLY::monolithic_unchained_style::monolithic_unchained_style() 
{ 
}

MEDDLY::compute_table* 
MEDDLY::monolithic_unchained_style::create(const ct_initializer::settings &s) const 
{
#ifdef OLD_OP_CT
  return new ct_old<true, false>(s, 0, 0);
#else
  switch (s.compression) {
    case ct_initializer::None:
                                    return new ct_none<true, false>(s, 0, 0);
    case ct_initializer::TypeBased:
                                    return new ct_typebased<true, false>(s, 0, 0);
    default:
                                    return 0;
  }
#endif
}

bool MEDDLY::monolithic_unchained_style::usesMonolithic() const 
{
  return true;
}

// **********************************************************************
// *                                                                    *
// *                  operation_chained_style  methods                  *
// *                                                                    *
// **********************************************************************

MEDDLY::operation_chained_style::operation_chained_style() 
{ 
}

MEDDLY::compute_table* 
MEDDLY::operation_chained_style::create(const ct_initializer::settings &s, operation* op, unsigned slot) const 
{
#ifdef OLD_OP_CT
  return new ct_old<false, true>(s, op, slot);
#else
  switch (s.compression) {
    case ct_initializer::None:
                                    return new ct_none<false, true>(s, 0, 0);
    case ct_initializer::TypeBased:
                                    return new ct_typebased<false, true>(s, 0, 0);
    default:
                                    return 0;
  }
#endif
}

bool MEDDLY::operation_chained_style::usesMonolithic() const 
{
  return false;
}


// **********************************************************************
// *                                                                    *
// *                 operation_unchained_style  methods                 *
// *                                                                    *
// **********************************************************************


MEDDLY::operation_unchained_style::operation_unchained_style() 
{ 
}

MEDDLY::compute_table* 
MEDDLY::operation_unchained_style::create(const ct_initializer::settings &s, operation* op, unsigned slot) const 
{
#ifdef OLD_OP_CT
  return new ct_old<false, false>(s, op, slot);
#else
  switch (s.compression) {
    case ct_initializer::None:
                                    return new ct_none<false, false>(s, 0, 0);
    case ct_initializer::TypeBased:
                                    return new ct_typebased<false, false>(s, 0, 0);
    default:
                                    return 0;
  }
#endif
}

bool MEDDLY::operation_unchained_style::usesMonolithic() const 
{
  return false;
}



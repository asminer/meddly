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

#include "ct_entry_type.h"

// **********************************************************************
// *                                                                    *
// *                         ct_object  methods                         *
// *                                                                    *
// **********************************************************************

MEDDLY::ct_object::ct_object()
{
}

MEDDLY::ct_object::~ct_object()
{
}

// **********************************************************************
// *                                                                    *
// *                          Helper functions                          *
// *                                                                    *
// **********************************************************************

inline unsigned bytes4typeID(MEDDLY::compute_table::typeID t)
{
  switch (t) {
    case MEDDLY::compute_table::NODE      : return sizeof(MEDDLY::node_handle);
    case MEDDLY::compute_table::INTEGER   : return sizeof(int);
    case MEDDLY::compute_table::LONG      : return sizeof(long);
    case MEDDLY::compute_table::FLOAT     : return sizeof(float);
    case MEDDLY::compute_table::DOUBLE    : return sizeof(double);
    case MEDDLY::compute_table::GENERIC   : return sizeof(MEDDLY::ct_object*);
    default:    return 0;
  }
}

inline char typeID2char(MEDDLY::compute_table::typeID t)
{
  switch (t) {
    case MEDDLY::compute_table::NODE      : return 'N';
    case MEDDLY::compute_table::INTEGER   : return 'I';
    case MEDDLY::compute_table::LONG      : return 'L';
    case MEDDLY::compute_table::FLOAT     : return 'F';
    case MEDDLY::compute_table::DOUBLE    : return 'D';
    case MEDDLY::compute_table::GENERIC   : return 'G';
    default:    return '?';
  }
}

inline MEDDLY::compute_table::typeID char2typeID(char c)
{
  switch (c) {
    case 'N':   return MEDDLY::compute_table::NODE;
    case 'I':   return MEDDLY::compute_table::INTEGER;
    case 'L':   return MEDDLY::compute_table::LONG;
    case 'F':   return MEDDLY::compute_table::FLOAT;
    case 'D':   return MEDDLY::compute_table::DOUBLE;
    case 'G':   return MEDDLY::compute_table::GENERIC;
    default:    return MEDDLY::compute_table::ERROR;
  }
}

// **********************************************************************
// *                                                                    *
// *                       ct_entry_type  methods                       *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table::entry_type::entry_type(const char* _name, const char* pattern)
{
  name = _name;
  is_marked_for_deletion = false;

  updatable_result = false;

  bool saw_dot = false;
  bool saw_colon = false;

  unsigned dot_slot = 0;
  unsigned colon_slot = 0;

  //
  // First scan: find '.' and ':'
  //

  unsigned length;
  for (length=0; pattern[length]; length++) {
    if ('.' == pattern[length]) {
      if (saw_dot || saw_colon) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
      saw_dot = true;
      dot_slot = length;
      continue;
    }
    if (':' == pattern[length]) {
      if (saw_colon) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
      saw_colon = true;
      colon_slot = length;
      continue;
    }
  }

  //
  // Determine lengths for each portion
  //

  if (!saw_colon) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);

  if (saw_dot) {
    // "012.456:89"
    len_ks_type = dot_slot;
    MEDDLY_DCASSERT(colon_slot > dot_slot);
    len_kr_type = (colon_slot - dot_slot) - 1;
  } else {
    // "01234:67"
    len_ks_type = colon_slot;
    len_kr_type = 0;
  }
  len_r_type = (length - colon_slot) - 1;

  if (
    (saw_dot && 0==len_kr_type)             // "foo.:bar" is bad
    ||
    (len_ks_type + len_kr_type == 0)        // no key?
    ||
    (len_r_type == 0)                       // no result?
  ) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);

  //
  // Build starting portion of key
  //
  ks_bytes = 0;
  if (len_ks_type) {
    ks_type = new typeID[len_ks_type];
    ks_forest = new expert_forest*[len_ks_type];
    for (unsigned i=0; i<len_ks_type; i++) {
      ks_type[i] = char2typeID(pattern[i]);
      ks_bytes += bytes4typeID(ks_type[i]);
      ks_forest[i] = 0;
    }
  } else {
    // This is possible if the pattern begins with .
    ks_type = 0;
    ks_forest = 0;
  }

  //
  // Build repeating portion of key
  //
  kr_bytes = 0;
  if (len_kr_type) {
    kr_type = new typeID[len_kr_type];
    kr_forest = new expert_forest*[len_kr_type];
    for (unsigned i=0; i<len_kr_type; i++) {
      kr_type[i] = char2typeID(pattern[i + dot_slot + 1]);
      kr_bytes += bytes4typeID(kr_type[i]);
      kr_forest[i] = 0;
    }
  } else {
    kr_type = 0;
    kr_forest = 0;
  }

  //
  // Build result
  //
  MEDDLY_DCASSERT(len_r_type);
  r_bytes = 0;
  r_type = new typeID[len_r_type];
  r_forest = new expert_forest*[len_r_type];
  for (unsigned i=0; i<len_r_type; i++) {
    r_type[i] = char2typeID(pattern[i + colon_slot + 1]);
    r_bytes += bytes4typeID(r_type[i]);
    r_forest[i] = 0;
  }

#ifdef DEBUG_ENTRY_TYPE
  printf("Built entry type %s with pattern '%s'\n", name, pattern);
  printf("Key start: \"");
  for (unsigned i=0; i<len_ks_type; i++) {
    fputc(typeID2char(ks_type[i]), stdout);
  }
  printf("\"  (%u bytes)\n", ks_bytes);
  printf("Key repeat: \"");
  for (unsigned i=0; i<len_kr_type; i++) {
    fputc(typeID2char(kr_type[i]), stdout);
  }
  printf("\"  (%u bytes)\n", kr_bytes);
  printf("Result: \"");
  for (unsigned i=0; i<len_r_type; i++) {
    fputc(typeID2char(r_type[i]), stdout);
  }
  printf("\"  (%u bytes)\n", r_bytes);
#endif
}

MEDDLY::compute_table::entry_type::~entry_type()
{
  delete[] ks_type;
  delete[] ks_forest;
  delete[] kr_type;
  delete[] kr_forest;
  delete[] r_type;
  delete[] r_forest;
}

void MEDDLY::compute_table::entry_type::setForestForSlot(unsigned i, expert_forest* f)
{
  if (i<len_ks_type) {
    if (NODE != ks_type[i]) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
    ks_forest[i] = f;
    return;
  }
  i -= len_ks_type;

  if (len_kr_type) {
    // adjust for the .
    i--;
    if (i<len_kr_type) {
      if (NODE != kr_type[i]) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
      kr_forest[i] = f;
      return;
    }
  }

  i -= len_kr_type;
  // adjust for :
  i--;

  if (i < len_r_type) {
    if (NODE != r_type[i]) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
    r_forest[i] = f;
    return;
  }

  // i is too large
  throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
}

void MEDDLY::compute_table::entry_type::clearForestCTBits(bool* skipF, unsigned N) const
{
  unsigned i;
  for (i=0; i<len_ks_type; i++) {
    expert_forest* f = ks_forest[i];
    if (0==f) continue;
    MEDDLY_DCASSERT(f->FID() < N);
    if (skipF[ f->FID() ]) continue;
    f->clearAllCacheBits();
    skipF[ f->FID() ] = 1;
  }
  for (i=0; i<len_kr_type; i++) {
    expert_forest* f = kr_forest[i];
    if (0==f) continue;
    MEDDLY_DCASSERT(f->FID() < N);
    if (skipF[ f->FID() ]) continue;
    f->clearAllCacheBits();
    skipF[ f->FID() ] = 1;
  }
  for (i=0; i<len_r_type; i++) {
    expert_forest* f = r_forest[i];
    if (0==f) continue;
    MEDDLY_DCASSERT(f->FID() < N);
    if (skipF[ f->FID() ]) continue;
    f->clearAllCacheBits();
    skipF[ f->FID() ] = 1;
  }
}

void MEDDLY::compute_table::entry_type::sweepForestCTBits(bool* whichF, unsigned N) const
{
  unsigned i;
  for (i=0; i<len_ks_type; i++) {
    expert_forest* f = ks_forest[i];
    if (0==f) continue;
    MEDDLY_DCASSERT(f->FID() < N);
    if (whichF[ f->FID() ]) {
      f->sweepAllCacheBits();
      whichF[ f->FID() ] = 0;
    }
  }
  for (i=0; i<len_kr_type; i++) {
    expert_forest* f = kr_forest[i];
    if (0==f) continue;
    MEDDLY_DCASSERT(f->FID() < N);
    if (whichF[ f->FID() ]) {
      f->sweepAllCacheBits();
      whichF[ f->FID() ] = 0;
    }
  }
  for (i=0; i<len_r_type; i++) {
    expert_forest* f = r_forest[i];
    if (0==f) continue;
    MEDDLY_DCASSERT(f->FID() < N);
    if (whichF[ f->FID() ]) {
      f->sweepAllCacheBits();
      whichF[ f->FID() ] = 0;
    }
  }
}


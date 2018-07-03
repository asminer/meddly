
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


#include "defines.h"

#include "storage/ct_classic.h"

#define DEBUG_ENTRY_TYPE

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

// ******************************************************************
// *                                                                *
// *                     ct_initializer  methods                    *
// *                                                                *
// ******************************************************************

MEDDLY::ct_initializer::settings MEDDLY::ct_initializer::the_settings;
const MEDDLY::compute_table_style* MEDDLY::ct_initializer::ct_factory;
MEDDLY::compute_table_style* MEDDLY::ct_initializer::builtin_ct_factory;

MEDDLY::ct_initializer::ct_initializer(initializer_list* prev) : initializer_list(prev)
{
  ct_factory = 0;
  builtin_ct_factory = 0;

  setBuiltinStyle(MonolithicUnchainedHash);
  setMaxSize(16777216);
  setStaleRemoval(Moderate);

  //
  // Set to null for now.
  // Set to the proper manager in setup()
  //
  setMemoryManager(0);  
}

MEDDLY::ct_initializer::~ct_initializer()
{
  delete builtin_ct_factory;
  builtin_ct_factory = 0;
}

void MEDDLY::ct_initializer::setup()
{
  if (0==ct_factory) throw error(error::INVALID_ASSIGNMENT, __FILE__, __LINE__);

  MEDDLY_DCASSERT(FREELISTS);
  setMemoryManager(FREELISTS);

  if (ct_factory->usesMonolithic()) {
    operation::Monolithic_CT = ct_factory->create(the_settings);
  }

  compute_table::initialize();
}

void MEDDLY::ct_initializer::cleanup()
{
  delete operation::Monolithic_CT;
  operation::Monolithic_CT = 0;

  compute_table::destroy();
}

void MEDDLY::ct_initializer::setStaleRemoval(staleRemovalOption sro)
{
  the_settings.staleRemoval = sro;
}

void MEDDLY::ct_initializer::setMaxSize(unsigned ms)
{
  the_settings.maxSize = ms;
}

void MEDDLY::ct_initializer::setBuiltinStyle(builtinCTstyle cts)
{
  delete builtin_ct_factory;
  builtin_ct_factory = 0;
  switch (cts) {
    case MonolithicUnchainedHash:
          builtin_ct_factory = new monolithic_unchained_style;
          break;

    case MonolithicChainedHash:
          builtin_ct_factory = new monolithic_chained_style;
          break;

    case OperationUnchainedHash:
          builtin_ct_factory = new operation_unchained_style;
          break;

    case OperationChainedHash:
          builtin_ct_factory = new operation_chained_style;
          break;
  }

  ct_factory = builtin_ct_factory;
}

void MEDDLY::ct_initializer::setUserStyle(const compute_table_style* cts)
{
  delete builtin_ct_factory;
  builtin_ct_factory = 0;
  ct_factory = cts;
}

void MEDDLY::ct_initializer::setMemoryManager(const memory_manager_style* mms)
{
  the_settings.MMS = mms;
}

#ifdef OLD_OP_CT
MEDDLY::compute_table* MEDDLY::ct_initializer::createForOp(operation* op)
{
  if (ct_factory) {
    return ct_factory->create(the_settings, op, 0);
  } else {
    return 0;
  }
}
#else
MEDDLY::compute_table* MEDDLY::ct_initializer::createForOp(operation* op, unsigned slot)
{
  if (ct_factory) {
    return ct_factory->create(the_settings, op, slot);
  } else {
    return 0;
  }
}
#endif

// **********************************************************************
// *                                                                    *
// *                    compute_table_style  methods                    *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table_style::compute_table_style()
{
}

MEDDLY::compute_table_style::~compute_table_style()
{
}

MEDDLY::compute_table* 
MEDDLY::compute_table_style::create(const ct_initializer::settings &s)
      const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}


MEDDLY::compute_table* 
MEDDLY::compute_table_style::create(const ct_initializer::settings &s, 
      operation* op, unsigned slot) const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

// **********************************************************************
// *                                                                    *
// *                       compute_table  methods                       *
// *                                                                    *
// **********************************************************************

#ifndef OLD_OP_CT
MEDDLY::compute_table::entry_type** MEDDLY::compute_table::entryInfo;
unsigned MEDDLY::compute_table::entryInfoAlloc;
unsigned MEDDLY::compute_table::entryInfoSize;
#endif
MEDDLY::compute_table::entry_key* MEDDLY::compute_table::free_keys;

MEDDLY::compute_table::compute_table(const ct_initializer::settings &s)
{
  maxSize = s.maxSize;
  if (0==maxSize)
    throw error(error::INVALID_ASSIGNMENT, __FILE__, __LINE__);

  switch (s.staleRemoval) {
    case ct_initializer::Aggressive:
            checkStalesOnFind = true;
            checkStalesOnResize = true;
            break;
    case ct_initializer::Moderate:
            checkStalesOnFind = false;
            checkStalesOnResize = true;
            break;
    case ct_initializer::Lazy:
            checkStalesOnFind = false;
            checkStalesOnResize = false;
            break;
  }

  perf.numEntries = 0;
  perf.hits = 0;
  perf.pings = 0;
  perf.numLargeSearches = 0;
  perf.maxSearchLength = 0;
  for (int i=0; i<perf.searchHistogramSize; i++)
    perf.searchHistogram[i] = 0;
}

MEDDLY::compute_table::~compute_table()
{
}

void MEDDLY::compute_table::initialize()
{
  free_keys = 0;
#ifndef OLD_OP_CT
  //
  // Initialize entryInfo list
  //
  entryInfo = 0;
  entryInfoAlloc = 0;
  entryInfoSize = 0;
  // zero that array?
#endif
}

void MEDDLY::compute_table::destroy()
{
  while (free_keys) {
    entry_key* n = free_keys->next;
    delete free_keys;
    free_keys = n;
  }
#ifndef OLD_OP_CT
  // delete the items?  TBD
  delete[] entryInfo;
#endif
}

#ifndef OLD_OP_CT

void MEDDLY::compute_table::registerOp(operation* op, unsigned num_ids)
{
  if (0==op) return;
  if (0==num_ids) return;

  if (1==num_ids) {
    //
    // Most common case, and easiest to find a hole.
    //
    for (unsigned i=0; i<entryInfoSize; i++) {
      if (0==entryInfo[i]) {
        op->setFirstETid(i);
        return;
      }
    } // for i
  } // if 1==num_ids

  //
  // No holes, or we want more than one slot so we didn't bother looking for holes.
  // Grab slots off the end
  //
  op->setFirstETid(entryInfoSize);
  entryInfoSize += num_ids;

  if (entryInfoSize > entryInfoAlloc) {
    //
    // Need to enlarge
    //
    entry_type** net = new entry_type* [entryInfoAlloc+256];
    unsigned i;
    for (i=0; i<entryInfoAlloc; i++) {
      net[i] = entryInfo[i];
    }
    entryInfoAlloc += 256;
    for (; i<entryInfoAlloc; i++) {
      net[i] = 0;
    }
    delete[] entryInfo;
    entryInfo = net;
  }
}

void MEDDLY::compute_table::registerEntryType(unsigned etid, entry_type* et)
{
  MEDDLY_CHECK_RANGE(0, etid, entryInfoSize);
  MEDDLY_DCASSERT(0==entryInfo[etid]);
  entryInfo[etid] = et;
}

void MEDDLY::compute_table::unregisterOp(operation* op, unsigned num_ids)
{
  if (0==op) return;
  if (0==num_ids) return;
  MEDDLY_CHECK_RANGE(0, op->getFirstETid(), entryInfoSize);
  unsigned stopID = op->getFirstETid()+num_ids;
  for (unsigned i=op->getFirstETid(); i<stopID; i++) {
    delete entryInfo[i];
    entryInfo[i] = 0;
  }
  //
  // decrement entryInfoSize if array ends in zeroes
  //
  while (entryInfoSize && (0==entryInfo[entryInfoSize-1])) entryInfoSize--;
}

#endif

// **********************************************************************

MEDDLY::compute_table::entry_key::entry_key()
{
  op = 0;
  data_alloc = 8;
  data = (int*) malloc(data_alloc * sizeof(int));
  // malloc: because realloc later
}

MEDDLY::compute_table::entry_key::~entry_key()
{
  free(data);
}

// **********************************************************************

MEDDLY::compute_table::entry_result::entry_result()
{
  data = 0;
  build = 0;
  ansLength = 0;
}

MEDDLY::compute_table::entry_result::entry_result(unsigned slots)
{
  build = new node_handle[slots];
  data = build;
  ansLength = slots;
}

MEDDLY::compute_table::entry_result::~entry_result()
{
  delete[] build;
}

// **********************************************************************

MEDDLY::compute_table::entry_type::entry_type(const char* _name, const char* pattern)
{
  name = _name;

  bool saw_dot = false;
  bool saw_colon = false;

  unsigned dot_slot = 0;
  unsigned colon_slot = 0;

  //
  // First scan: find '.' and ':'
  //

  unsigned length;
  for (length=0; name[length]; length++) {
    if ('.' == name[length]) {
      if (saw_dot || saw_colon) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
      saw_dot = true;
      dot_slot = length;
      continue;
    }
    if (':' == name[length]) {
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
  if (len_ks_type) {
    ks_type = new char[len_ks_type];
    ks_forest = new forest*[len_ks_type];
    for (unsigned i=0; i<len_ks_type; i++) {
      ks_type[i] = name[i];
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
  if (len_kr_type) {
    kr_type = new char[len_kr_type];
    kr_forest = new forest*[len_kr_type];
    for (unsigned i=0; i<len_kr_type; i++) {
      kr_type[i] = name[i + dot_slot + 1];
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
  r_type = new char[len_r_type];
  r_forest = new forest*[len_r_type];
  for (unsigned i=0; i<len_r_type; i++) {
    r_type[i] = name[i + colon_slot + 1];
    r_forest[i] = 0;
  }

#ifdef DEBUG_ENTRY_TYPE
  printf("Built entry type %s with pattern %s:\n", name, pattern);
  printf("Key start: \"");
  for (unsigned i=0; i<len_ks_type; i++) {
    fputc(ks_type[i], stdout);
  }
  printf("\"\n");
  printf("Key repeat: \"");
  for (unsigned i=0; i<len_kr_type; i++) {
    fputc(kr_type[i], stdout);
  }
  printf("\"\n");
  printf("Result: \"");
  for (unsigned i=0; i<len_r_type; i++) {
    fputc(r_type[i], stdout);
  }
  printf("\"\n");
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

void MEDDLY::compute_table::entry_type::setForestForSlot(unsigned i, forest* f)
{
  if (i<len_ks_type) {
    if ('N' != ks_type[i]) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
    ks_forest[i] = f;
    return;
  }
  i -= len_ks_type;

  if (len_kr_type) {
    // adjust for the .
    i--;
    if (i<len_kr_type) {
      if ('N' != kr_type[i]) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
      kr_forest[i] = f;
      return;
    }
  }
  
  i -= len_kr_type;
  // adjust for :
  i--;

  if (i < len_r_type) {
    if ('N' != r_type[i]) throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
    r_forest[i] = f;
    return;
  }

  // i is too large
  throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
}


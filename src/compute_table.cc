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
#include "old_meddly.h"
#include "old_meddly.hh"
#include "old_meddly_expert.h"
#include "old_meddly_expert.hh"

#include "storage/ct_styles.h"

// #define DEBUG_ENTRY_TYPE
// #define DEBUG_ENTRY_REGISTRY

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
//  setCompression(None);
  setCompression(TypeBased);

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

void MEDDLY::ct_initializer::setCompression(compressionOption co)
{
  the_settings.compression = co;
}

void MEDDLY::ct_initializer::setMemoryManager(const memory_manager_style* mms)
{
  the_settings.MMS = mms;
}

MEDDLY::compute_table* MEDDLY::ct_initializer::createForOp(operation* op, unsigned slot)
{
  if (ct_factory) {
    return ct_factory->create(the_settings, op, slot);
  } else {
    return 0;
  }
}

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


MEDDLY::ct_entry_type** MEDDLY::compute_table::entryInfo;
unsigned MEDDLY::compute_table::entryInfoAlloc;
unsigned MEDDLY::compute_table::entryInfoSize;
MEDDLY::ct_entry_key* MEDDLY::compute_table::free_keys;

MEDDLY::compute_table::compute_table(const ct_initializer::settings &s,
  operation* op, unsigned slot)
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
  for (unsigned i=0; i<perf.searchHistogramSize; i++) {
    perf.searchHistogram[i] = 0;
  }
  perf.resizeScans = 0;

  //
  // Global operation vs monolithic
  //
  if (op) {
    global_et = getEntryType(op, slot);
    MEDDLY_DCASSERT(global_et);
  } else {
    global_et = 0;
  }
}

MEDDLY::compute_table::~compute_table()
{
}

void MEDDLY::compute_table::initialize()
{
  free_keys = 0;
  //
  // Initialize entryInfo list
  //
  entryInfo = 0;
  entryInfoAlloc = 0;
  entryInfoSize = 0;
  // zero that array?
}

void MEDDLY::compute_table::destroy()
{
  while (free_keys) {
    ct_entry_key* n = free_keys->next;
    delete free_keys;
    free_keys = n;
  }
  // delete the items?  TBD
  delete[] entryInfo;
}

void MEDDLY::compute_table::clearForestCTBits(bool* skipF, unsigned N) const
{
  if (global_et) {
    // Operation cache
    global_et->clearForestCTBits(skipF, N);
    return;
  }
  //
  // Monolithic cache.
  //
  for (unsigned i=0; i<entryInfoSize; i++) {
    if (entryInfo[i]) {
      entryInfo[i]->clearForestCTBits(skipF, N);
    }
  }
}

void MEDDLY::compute_table::sweepForestCTBits(bool* whichF, unsigned N) const
{
  if (global_et) {
    // Operation cache
    global_et->sweepForestCTBits(whichF, N);
    return;
  }
  //
  // Monolithic cache.
  //
  for (unsigned i=0; i<entryInfoSize; i++) {
    if (entryInfo[i]) {
      entryInfo[i]->sweepForestCTBits(whichF, N);
    }
  }
}

void MEDDLY::compute_table::registerOp(operation* op, unsigned num_ids)
{
  if (0==op) return;
  if (0==num_ids) return;

#ifdef DEBUG_ENTRY_REGISTRY
  printf("Requesting %u entry slots for operation %s\n", num_ids, op->getName());
#endif

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
    ct_entry_type** net = new ct_entry_type* [entryInfoAlloc+256];
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

void MEDDLY::compute_table::registerEntryType(unsigned etid, ct_entry_type* et)
{
    MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, etid, entryInfoSize);
  MEDDLY_DCASSERT(0==entryInfo[etid]);
  entryInfo[etid] = et;
  et->etID = etid;
#ifdef DEBUG_ENTRY_REGISTRY
  printf("Registering entry %s in slot %u\n", et->getName(), etid);
  printf("Updated Table:\n");
  for (unsigned i=0; i<entryInfoSize; i++) {
    printf("    %2u: %s\n", i, entryInfo[i] ? entryInfo[i]->getName() : " ");
  }
#endif
}

void MEDDLY::compute_table::unregisterOp(operation* op, unsigned num_ids)
{
  if (0==op) return;
  if (0==num_ids) return;
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, op->getFirstETid(), entryInfoSize);
  unsigned stopID = op->getFirstETid()+num_ids;
  for (unsigned i=op->getFirstETid(); i<stopID; i++) {
    delete entryInfo[i];
    entryInfo[i] = 0;
  }
  //
  // decrement entryInfoSize if array ends in zeroes
  //
  while (entryInfoSize && (0==entryInfo[entryInfoSize-1])) entryInfoSize--;
#ifdef DEBUG_ENTRY_REGISTRY
  printf("Unregistering %u slots for operation %s\n", num_ids, op->getName());
  printf("Updated Table:\n");
  for (unsigned i=0; i<entryInfoSize; i++) {
    printf("    %2u: %s\n", i, entryInfo[i] ? entryInfo[i]->getName() : " ");
  }
#endif
}

// **********************************************************************


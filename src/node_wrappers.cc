
// $Id$

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



// TODO: Testing

#include "defines.h"
#include "hash_stream.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_reader  methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_reader::node_reader()
{
  extra_hashed = 0;
  ext_alloc = 0;
  ext_size = 0;
  down = 0;
  index = 0;
  edge = 0;
  alloc = 0;
  ealloc = 0;
  size = 0;
  nnzs = 0;
  level = 0;
#ifdef DEVELOPMENT_CODE
  has_hash = false;
#endif
}

MEDDLY::node_reader::~node_reader()
{
  clear();
}

void MEDDLY::node_reader::clear()
{
  free(extra_hashed);
  free(down);
  free(index);
  free(edge);
  down = 0;
  index = 0;
  edge = 0;
  alloc = 0;
  ealloc = 0;
  size = 0;
  nnzs = 0;
  level = 0;
}

void MEDDLY::node_reader
::resize(int k, int ns, char eb, bool full)
{
  level = k;
  is_full = full;
  size = ns;
  edge_bytes = eb;
  if (size > alloc) {
    int nalloc = ((ns/8)+1)*8;
    MEDDLY_DCASSERT(nalloc > ns);
    MEDDLY_DCASSERT(nalloc>0);
    MEDDLY_DCASSERT(nalloc>alloc);
    down = (long*) realloc(down, nalloc*sizeof(long));
    if (0==down) throw error(error::INSUFFICIENT_MEMORY);
    index = (int*) realloc(index, nalloc*sizeof(int));
    if (0==index) throw error(error::INSUFFICIENT_MEMORY);
    alloc = nalloc;
  }
  if (edge_bytes * size > ealloc) {
    int nalloc = ((edge_bytes * size)/8+1)*8;
    MEDDLY_DCASSERT(nalloc>0);
    MEDDLY_DCASSERT(nalloc>ealloc);
    edge = realloc(edge, nalloc);
    if (0==edge) throw error(error::INSUFFICIENT_MEMORY);
    ealloc = nalloc;
  }
}

void MEDDLY::node_reader
::resize_header(int extra_slots)
{
  ext_size = extra_slots;
  if (ext_size > ext_alloc) {
    ext_alloc = ((ext_size/8)+1)*8;
    MEDDLY_DCASSERT(ext_alloc > ext_size);
    MEDDLY_DCASSERT(ext_alloc>0);
    extra_hashed = (int*) realloc(extra_hashed, ext_alloc * sizeof(int));
    if (0==extra_hashed) throw error(error::INSUFFICIENT_MEMORY);
  }
}

void MEDDLY::node_reader::computeHash(bool hashEdgeValues)
{
#ifdef DEVELOPMENT_CODE
  MEDDLY_DCASSERT(!has_hash);
#endif
  
  hash_stream s;
  s.start(level);

  for (int e=0; e<ext_size; e++) {
    MEDDLY_DCASSERT(extra_hashed);
    s.push(extra_hashed[e]);
  }
  
  if (isSparse()) {
    if (hashEdgeValues) {
      for (int z=0; z<nnzs; z++) {
        MEDDLY_DCASSERT(d(z));
        s.push(i(z), d(z), ei(z));
      }
    } else {
      for (int z=0; z<nnzs; z++) {
        MEDDLY_DCASSERT(d(z));
        s.push(i(z), d(z));
      }
    }
  } else {
    if (hashEdgeValues) {
      for (int n=0; n<size; n++) {
        if (0==d(n)) continue;
        s.push(n, d(n), ei(n));
      }
    } else {
      for (int n=0; n<size; n++) {
        if (0==d(n)) continue;
        s.push(n, d(n));
      }
    }
  }
  h = s.finish();
#ifdef DEVELOPMENT_CODE
  has_hash = true;
#endif
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_builder methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_builder::node_builder()
{
  parent = 0;
  extra_hashed = 0;
  extra_unhashed = 0;
  down = 0;
  indexes = 0;
  edge = 0;
}

MEDDLY::node_builder::~node_builder()
{
  free(extra_hashed);
  free(extra_unhashed);
  free(down);
  free(indexes);
  free(edge);
}

void MEDDLY::node_builder::init(int k, const expert_forest* p)
{
  MEDDLY_DCASSERT(p);
  parent = p;
  level = k;
  hhsize = parent->hashedHeaderSize();
  if (hhsize) {
    extra_hashed = (int*) malloc(uhsize * sizeof(int));
    if (0==extra_hashed)
      throw error(error::INSUFFICIENT_MEMORY);
  }
  uhsize = parent->unhashedHeaderSize();
  if (uhsize) {
    extra_unhashed = (int*) malloc(hhsize * sizeof(int));
    if (0==extra_unhashed)
      throw error(error::INSUFFICIENT_MEMORY);
  }
  size = 0;
  alloc = 0;
  lock = false;
  hashEdgeValues = parent->areEdgeValuesHashed();
  edge_bytes = parent->edgeSize() * sizeof(int);
}


void MEDDLY::node_builder::computeHash()
{
#ifdef DEVELOPMENT_CODE
  MEDDLY_DCASSERT(!has_hash);
#endif
  
  hash_stream s;
  s.start(level);

  for (int e=0; e<hhsize; e++) {
    MEDDLY_DCASSERT(extra_hashed);
    s.push(extra_hashed[e]);
  }
  if (isSparse()) {
    if (hashEdgeValues) {
      for (int z=0; z<size; z++) {
        MEDDLY_DCASSERT(d(z));
        s.push(i(z), d(z), ei(z));
      }
    } else {
      for (int z=0; z<size; z++) {
        MEDDLY_DCASSERT(d(z));
        s.push(i(z), d(z));
      }
    }
  } else {
    if (hashEdgeValues) {
      for (int n=0; n<size; n++) {
        if (0==d(n)) continue;
        s.push(n, d(n), ei(n));
      }
    } else {
      for (int n=0; n<size; n++) {
        if (0==d(n)) continue;
        s.push(n, d(n));
      }
    }
  }
  h = s.finish();
#ifdef DEVELOPMENT_CODE
  has_hash = true;
#endif
}


void MEDDLY::node_builder::enlarge()
{
  if (size <= alloc) return;
  alloc = ((size / 8)+1) * 8;
  MEDDLY_DCASSERT(alloc > size);
  down = (long*) realloc(down, alloc * sizeof(long));
  if (0==down) throw error(error::INSUFFICIENT_MEMORY);
  indexes = (int*) realloc(indexes, alloc * sizeof(int));
  if (0==indexes) throw error(error::INSUFFICIENT_MEMORY);
  int es = parent->edgeSize();
  if (es>0) {
    edge = realloc(edge, alloc * es * sizeof(int));
    if (0==edge) throw error(error::INSUFFICIENT_MEMORY);
  }
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_storage methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_storage::node_storage()
{
  parent = 0;
  stats = 0;
  counts = 0;
  nexts = 0;
}

MEDDLY::node_storage::~node_storage()
{
  // nothing, derived classes must handle everything
}

void MEDDLY::node_storage::initForForest(expert_forest* f)
{
  MEDDLY_DCASSERT(0==parent);
  parent = f;
  stats = &parent->changeStats();
}


void MEDDLY::node_storage::dumpInternal(FILE* s) const
{
  dumpInternalInfo(s);
  fprintf(s, "Data array by record:\n");
  for (long a=1; a > 0; ) {
    fflush(s);
    a = dumpInternalNode(s, a);
  } // for a
  fflush(s);
}


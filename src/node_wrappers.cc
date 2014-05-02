
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
::show(FILE* s, const expert_forest* parent, bool verb) const
{
  int stop;
  if (isSparse()) {
    if (verb) fprintf(s, "nnzs: %d ", size);
    fprintf(s, "down: (");
    stop = nnzs;
  } else {
    if (verb) fprintf(s, "size: %d ", size);
    fprintf(s, "down: [");
    stop = size;
  }

  for (int z=0; z<stop; z++) {
    if (isSparse()) {
      if (z) fprintf(s, ", ");
      fprintf(s, "%d:", i(z));
    } else {
      if (z) fprintf(s, "|");
    }
    if (parent->edgeBytes()) {
      fprintf(s, "<");
      parent->showEdgeValue(s, eptr(z));
      fprintf(s, ", ");
    }
    if (parent->isTerminalNode(d(z))) {
      parent->showTerminal(s, d(z));
    } else {
      fprintf(s, "%ld", long(d(z)));
    }
    if (parent->edgeBytes()) fprintf(s, ">");
  }

  if (isSparse()) {
    fputc(')', s);
  } else {
    fputc(']', s);
  }
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
    down = (node_handle*) realloc(down, nalloc*sizeof(node_handle));
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
::resize_header(int extra_bytes)
{
  ext_size = extra_bytes;
  if (ext_size > ext_alloc) {
    ext_alloc = ((ext_size/8)+1)*8;
    MEDDLY_DCASSERT(ext_alloc > ext_size);
    MEDDLY_DCASSERT(ext_alloc>0);
    extra_hashed =  realloc(extra_hashed, ext_alloc);
    if (0==extra_hashed) throw error(error::INSUFFICIENT_MEMORY);
  }
}

void MEDDLY::node_reader::computeHash(bool hashEdgeValues, node_handle tv)
{
#ifdef DEVELOPMENT_CODE
  MEDDLY_DCASSERT(!has_hash);
#endif
  
  hash_stream s;
  s.start(level);

  if (ext_size) {
    s.push(extra_hashed, ext_size);
  }
  
  if (isSparse()) {
    if (hashEdgeValues) {
      for (int z=0; z<nnzs; z++) {
        MEDDLY_DCASSERT(d(z)!=tv);
        s.push(i(z), d(z), ei(z));
      }
    } else {
      for (int z=0; z<nnzs; z++) {
        MEDDLY_DCASSERT(d(z)!=tv);
        s.push(i(z), d(z));
      }
    }
  } else {
    if (hashEdgeValues) {
      for (int n=0; n<size; n++) {
        if (d(n)!=tv) {
          s.push(n, d(n), ei(n));
        }
      }
    } else {
      for (int n=0; n<size; n++) {
        if (d(n)!=tv) {
          s.push(n, d(n));
        }
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
  hhbytes = 0;
  extra_unhashed = 0;
  down = 0;
  indexes = 0;
  edge = 0;
#ifdef DEVELOPMENT_CODE
  has_hash = false;
#endif
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
  hhbytes = parent->hashedHeaderBytes();
  if (hhbytes) {
    extra_hashed = malloc(hhbytes);
    if (0==extra_hashed)
      throw error(error::INSUFFICIENT_MEMORY);
  }
  if (parent->unhashedHeaderBytes()) {
    extra_unhashed = malloc(parent->unhashedHeaderBytes());
    if (0==extra_unhashed)
      throw error(error::INSUFFICIENT_MEMORY);
  }
  size = 0;
  alloc = 0;
  lock = false;
  edge_bytes = parent->edgeBytes();
}

void MEDDLY::node_builder::show(FILE* s, bool verb) const
{
  if (isSparse()) {
    if (verb) fprintf(s, "nnzs: %d ", size);
    fprintf(s, "down: (");
  } else {
    if (verb) fprintf(s, "size: %d ", size);
    fprintf(s, "down: [");
  }

  for (int z=0; z<size; z++) {
    if (isSparse()) {
      if (z) fprintf(s, ", ");
      fprintf(s, "%d:", i(z));
    } else {
      if (z) fprintf(s, "|");
    }
    if (parent->edgeBytes()) {
      fprintf(s, "<");
      parent->showEdgeValue(s, eptr(z));
      fprintf(s, ", ");
    }
    if (parent->isTerminalNode(d(z))) {
      parent->showTerminal(s, d(z));
    } else {
      fprintf(s, "%ld", long(d(z)));
    }
    if (parent->edgeBytes()) fprintf(s, ">");
  }

  if (isSparse()) {
    fputc(')', s);
  } else {
    fputc(']', s);
  }
}

void MEDDLY::node_builder::computeHash()
{
#ifdef DEVELOPMENT_CODE
  MEDDLY_DCASSERT(!has_hash);
#endif
  
  hash_stream s;
  s.start(level);

  if (hhbytes) {
    s.push(extra_hashed, hhbytes);
  }

  if (isSparse()) {
    if (parent->areEdgeValuesHashed()) {
      for (int z=0; z<size; z++) {
        MEDDLY_DCASSERT(d(z)!=parent->getTransparentNode());
        s.push(i(z), d(z), ei(z));
      }
    } else {
      for (int z=0; z<size; z++) {
        MEDDLY_DCASSERT(d(z)!=parent->getTransparentNode());
        s.push(i(z), d(z));
      }
    }
  } else {
	node_handle tv=parent->getTransparentNode();
    if (parent->areEdgeValuesHashed()) {
      for (int n=0; n<size; n++) {
        if (d(n)!=tv) {
          s.push(n, d(n), ei(n));
        }
      }
    } else {
      for (int n=0; n<size; n++) {
        if (d(n)!=tv) {
          s.push(n, d(n));
        }
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
  down = (node_handle*) realloc(down, alloc * sizeof(node_handle));
  if (0==down) throw error(error::INSUFFICIENT_MEMORY);
  indexes = (int*) realloc(indexes, alloc * sizeof(int));
  if (0==indexes) throw error(error::INSUFFICIENT_MEMORY);
  if (parent->edgeBytes()>0) {
    edge = realloc(edge, alloc * parent->edgeBytes());
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
  localInitForForest(f);
}

void MEDDLY::node_storage
::writeNode(FILE* s, node_address, const node_handle*) const
{
  throw error(error::NOT_IMPLEMENTED);
}

void MEDDLY::node_storage::dumpInternal(FILE* s, unsigned flags) const
{
  dumpInternalInfo(s);
  fprintf(s, "Data array by record:\n");
  for (node_address a=1; a > 0; ) {
    fflush(s);
    a = dumpInternalNode(s, a, flags);
  } // for a
  dumpInternalTail(s);
  fflush(s);
}

void MEDDLY::node_storage::localInitForForest(const expert_forest* f)
{
  // default - do nothing
}

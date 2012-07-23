
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

/*
void MEDDLY::node_reader::computeHash()
{
  MEDDLY_DCASSERT(!has_hash);
  
  hash_stream s;
  s.start(level);

  if (is_full) {
    if (parent->areEdgeValuesHashed(level)) {
      for (int i=0; i<size; i++) {
        if (0==down[i]) continue;
        int* ep = (int*) ( (char*)edge + i * edge_bytes );
        s.push(i, down[i], *ep);
      }
    } else {
      for (int i=0; i<size; i++) {
        if (0==down[i]) continue;
        s.push(i, down[i]);
      }
    }
  } else {
    if (parent->areEdgeValuesHashed(level)) {
      for (int z=0; z<nnzs; z++) {
        MEDDLY_DCASSERT(down[z]);
        int* ep = (int*) ( (char*)edge + z * edge_bytes );
        s.push(index[z], down[z], *ep);
      }
    } else {
      for (int z=0; z<nnzs; z++) {
        MEDDLY_DCASSERT(down[z]);
        s.push(index[z], down[z]);
      }
    }
  }
  h = s.finish();
  has_hash = true;
}
*/

/*
void MEDDLY::node_reader::dump(FILE* s) const
{
  if (is_full) {
    fprintf(s, "[%d", down[0]);
    for (int i=1; i<size; i++)
      fprintf(s, ", %d", down[i]);
    fprintf(s, "]");
  } else {
    fprintf(s, "(%d:%d", index[0], down[0]);
    for (int z=1; z<nnzs; z++) 
      fprintf(s, ", %d:%d", index[z], down[z]);
    fprintf(s, ")");
  }
}
*/

void MEDDLY::node_reader
// ::resize(const expert_forest* p, int k, int ns, bool full)
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
    down = (int*) realloc(down, nalloc*sizeof(int));
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
  hhsize = parent->hashedHeaderSize(k);
  if (hhsize) {
    extra_hashed = (int*) malloc(uhsize * sizeof(int));
    if (0==extra_hashed)
      throw error(error::INSUFFICIENT_MEMORY);
  }
  uhsize = parent->unhashedHeaderSize(k);
  if (uhsize) {
    extra_unhashed = (int*) malloc(hhsize * sizeof(int));
    if (0==extra_unhashed)
      throw error(error::INSUFFICIENT_MEMORY);
  }
  size = 0;
  alloc = 0;
  lock = false;
  hashEdgeValues = parent->areEdgeValuesHashed(k);
  edge_bytes = parent->edgeSize(k) * sizeof(int);
}

/*
bool MEDDLY::node_builder::equals(int p) const
{
  const expert_forest::nodeData &node = parent->getNode(p);
  if (node.level != level) return false;

  MEDDLY_DCASSERT(0==edge);  // edge value stuff not implemented yet

  if (parent->hashedHeader) {
    const int* hh = parent->hashedHeaderOf(p);
    if (memcmp(extra_hashed, hh, parent->hashedHeader * sizeof(int)))
      return false;
  }

  // p is full:
  if (parent->isFull(node.offset)) {
    int fs = parent->fullSizeOf(node.offset);
    const int* pd = parent->fullDownOf(node.offset);
    if (is_sparse) {
      int i = 0;
      for (int z=0; z<size; z++) {
        if (indexes[z] >= fs) return false;
        for (; i<indexes[z]; i++) if (pd[i]) return false;
        if (down[z] != pd[i]) return false;
        i++;
      } // for z
      for (; i<fs; i++) if (pd[i]) return false;
      return true;
    }
    // we're also full
    if (fs > size) return false;
    if (memcmp(down, pd, fs * sizeof(int))) return false;
    for (int i=fs; i<size; i++) if (down[i]) return false;
    return true;
  }

  // p is sparse:
  int nnz = parent->sparseSizeOf(node.offset);
  const int* pd = parent->sparseDownOf(node.offset);
  const int* pi = parent->sparseIndexesOf(node.offset);

  if (is_sparse) {
    if (nnz != size) return false;
    if (memcmp(down, pd, nnz * sizeof(int))) return false;
    if (memcmp(indexes, pi, nnz * sizeof(int))) return false;
    return true;
  }

  // We must be full
  int i = 0;
  for (int z=0; z<nnz; z++) {
    for (; i<pi[z]; i++) if (down[i]) return false;
    if (down[i] != pd[z]) return false;
    i++;
  }
  for (; i<size; i++) if (down[i]) return false;
  return true;
}
*/

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
  
  if (is_sparse) {
    if (hashEdgeValues) {
      for (int z=0; z<size; z++) {
        MEDDLY_DCASSERT(down[z]);
        int* ep = (int*) ( (char*)edge + z * edge_bytes );
        s.push(indexes[z], down[z], *ep);
      }
    } else {
      for (int z=0; z<size; z++) {
        MEDDLY_DCASSERT(down[z]);
        s.push(indexes[z], down[z]);
      }
    }
  } else {
    if (hashEdgeValues) {
      for (int i=0; i<size; i++) {
        if (0==down[i]) continue;
        int* ep = (int*) ( (char*)edge + i * edge_bytes );
        s.push(i, down[i], *ep);
      }
    } else {
      for (int i=0; i<size; i++) {
        if (0==down[i]) continue;
        s.push(i, down[i]);
      }
    }
  }
  h = s.finish();
#ifdef DEVELOPMENT_CODE
  has_hash = true;
#endif
}


void MEDDLY::node_builder::copyIntoFull(int* d, int N) const
{
  MEDDLY_DCASSERT(0==edge_bytes);
  if (is_sparse) {
    memset(d, 0, N * sizeof(int));
    for (int z=0; z<size; z++) {
      if (indexes[z] < N) d[indexes[z]] = down[z];
    }
  } else {
    memcpy(d, down, N * sizeof(int));
  }
}

void MEDDLY::node_builder::copyIntoFull(int* d, void* e, int N) const
{
  MEDDLY_DCASSERT(edge_bytes);
  if (is_sparse) {
    memset(d, 0, N * sizeof(int));
    memset(e, 0, N * edge_bytes);
    for (int z=0; z<size; z++) {
      if (indexes[z] < N) {
        int i = indexes[z];
        d[i] = down[z];
        memcpy(((char*)e) + i * edge_bytes, ((char*)edge) + z * edge_bytes, edge_bytes);
      }
    }
  } else {
    memcpy(d, down, N * sizeof(int));
    memcpy(e, edge, N * edge_bytes);
  }
}

void MEDDLY::node_builder
::copyIntoSparse(int* d, int* ix, int Z) const
{
  MEDDLY_DCASSERT(0==edge_bytes);
  if (is_sparse) {
    memcpy(d, down, Z * sizeof(int));
    memcpy(ix, indexes, Z * sizeof(int));
  } else {
    int z=0;
    MEDDLY_DCASSERT(Z>0);
    for (int i=0; i<size; i++)
      if (down[i]) {
        d[z] = down[i];
        ix[z] = i;
        z++;
        if (z>=Z) return;
      }
  }
}

void MEDDLY::node_builder
::copyIntoSparse(int* d, int* ix, void* e, int Z) const
{
  MEDDLY_DCASSERT(edge_bytes);
  if (is_sparse) {
    memcpy(d, down, Z * sizeof(int));
    memcpy(ix, indexes, Z * sizeof(int));
    memcpy(e, edge, Z * edge_bytes);
  } else {
    int z=0;
    MEDDLY_DCASSERT(Z>0);
    for (int i=0; i<size; i++)
      if (down[i]) {
        d[z] = down[i];
        ix[z] = i;
        memcpy(((char*)e) + z * edge_bytes, ((char*)edge) + i * edge_bytes, edge_bytes);
        z++;
        if (z>=Z) return;
      }
  }
}

void MEDDLY::node_builder::enlarge()
{
  if (size <= alloc) return;
  alloc = ((size / 8)+1) * 8;
  MEDDLY_DCASSERT(alloc > size);
  down = (int*) realloc(down, alloc * sizeof(int));
  if (0==down) throw error(error::INSUFFICIENT_MEMORY);
  indexes = (int*) realloc(indexes, alloc * sizeof(int));
  if (0==indexes) throw error(error::INSUFFICIENT_MEMORY);
  int es = parent->edgeSize(level);
  if (es>0) {
    edge = realloc(edge, alloc * es * sizeof(int));
    if (0==edge) throw error(error::INSUFFICIENT_MEMORY);
  }
}


// ******************************************************************
// *                                                                *
// *                      node_finder  methods                      *
// *                                                                *
// ******************************************************************

/*
MEDDLY::node_finder::node_finder(const expert_forest* p, int n)
{
  parent = p;
  node = n;
  parent->initNodeReader(thisnode, node, false);
  h = parent->hashNode(node);
}

MEDDLY::node_finder::~node_finder()
{
}

bool MEDDLY::node_finder::equals(int p)
{
  parent->initNodeReader(compare, p, false);

  if (thisnode.getLevel() != compare.getLevel()) return false;
  if (thisnode.getNNZs() != compare.getNNZs()) return false;

  for (int z=0; z<thisnode.getNNZs(); z++) {
    if (thisnode.d(z) != compare.d(z)) return false;
    if (thisnode.i(z) != compare.i(z)) return false;
  }
  // TBD : edge values

  return true;
}

*/

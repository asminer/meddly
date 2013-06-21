
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
#include "unique_table.h"
#include "hash_stream.h"
//#include <set>
//#include <queue>
//#include <vector>

// #define DEBUG_CLEANUP

// #define DEBUG_ADDRESS_RESIZE
// #define DEBUG_GC

// #define TRACK_DELETIONS

// #define DEBUG_CREATE_REDUCED

// Thoroughly check reference counts.
// Very slow.  Use only for debugging.
// #define VALIDATE_INCOUNTS
// #define VALIDATE_INCOUNTS_ON_DELETE

// #define GC_OFF
// #define DEBUG_GC

// #define REPORT_ON_DESTROY

// #define LAZY_NODE_HASHING

const int a_min_size = 1024;

// ******************************************************************
// *                                                                *
// *                    forest::statset  methods                    *
// *                                                                *
// ******************************************************************

MEDDLY::forest::statset::statset()
{
  reclaimed_nodes = 0;
  num_compactions = 0;
  garbage_collections = 0;
  zombie_nodes = 0;
  orphan_nodes = 0;
  active_nodes = 0;
  peak_active = 0;
  memory_used = 0;
  memory_alloc = 0;
  peak_memory_used = 0;
  peak_memory_alloc = 0;
  memory_UT = 0;
  peak_memory_UT = 0;
  max_UT_chain = 0;
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                         forest methods                         *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::forest
::forest(int ds, domain* _d, bool rel, range_type t, edge_labeling ev, 
  const policies &p) : deflt(p)
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating forest #%d in domain #%d\n", ds, _d->ID());
#endif
  d_slot = ds;
  is_marked_for_deletion = false;
  d = _d;
  isRelation = rel;
  rangeType = t;
  edgeLabel = ev;
  // check policies
  if (!isRelation) {
    if (policies::IDENTITY_REDUCED == deflt.reduction)
      throw error(error::INVALID_POLICY);
  } 
  //
  // Initialize array of operations
  //
  opCount = 0;
  szOpCount = 0;
  //
  // Initialize list of registered dd_edges
  //
  firstHole = -1; // firstHole < 0 indicates no holes.
  firstFree = 0;
  sz = 256;

  // Create an array to store pointers to dd_edges.
  edge = (edge_data *) malloc(sz * sizeof(edge_data));
  for (unsigned i = 0; i < sz; ++i) {
    edge[i].nextHole = -1;
    edge[i].edge = 0;
  }
}

MEDDLY::forest::~forest()
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Deleting forest #%d in domain #%d\n", d_slot, d->ID());
#endif
  // operations are deleted elsewhere...
  free(opCount);
  // Make SURE our edges are orphaned
  for (unsigned i = 0; i < firstFree; ++i) {
    if (edge[i].edge) edge[i].edge->orphan();
  }
  free(edge);
  // NOTE: since the user is provided with the dd_edges instances (as opposed
  // to a pointer), the user program will automatically call the
  // destructor for each dd_edge when the corresponding variable goes out of
  // scope. Therefore there is no need to destruct dd_edges from here.
  d->unlinkForest(this, d_slot);
}

void MEDDLY::forest::markForDeletion()
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Marking forest #%d for deletion in domain #%d\n", d_slot, d->ID());
#endif
  is_marked_for_deletion = true;
  // deal with operations associated with this forest
  for (int i=0; i<szOpCount; i++) 
    if (opCount[i]) {
      operation* op = operation::getOpWithIndex(i);
      op->markForDeletion();
    }
  unregisterDDEdges();
}

void MEDDLY::forest::createEdgeForVar(int vh, bool vp, bool* terms, dd_edge& a)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdgeForVar(int vh, bool vp, int* terms, dd_edge& a)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdgeForVar(int vh, bool vp, float* terms, dd_edge& a)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdge(const int* const* vlist, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const int* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const float* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::createEdge(const int* const* vl, const int* const* vpl, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const int* const* vplist,
      const int* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const int* const* vplist,
      const float* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdge(bool val, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdge(int val, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdge(float val, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::evaluate(const dd_edge &f, const int* vl, bool &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::evaluate(const dd_edge &f, const int* vl, int &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::evaluate(const dd_edge &f, const int* vl, float &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::evaluate(const dd_edge& f, const int* vl, const int* vpl, bool &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::evaluate(const dd_edge& f, const int* vl, const int* vpl, int &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::evaluate(const dd_edge& f, const int* vl, const int* vpl, float &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::getElement(const dd_edge& a, int index, int* e)
{
  throw error(error::INVALID_OPERATION);
}

void MEDDLY::forest::removeStaleComputeTableEntries()
{
  if (operation::useMonolithicComputeTable()) {
    operation::removeStalesFromMonolithic();
  } else {
    for (int i=0; i<szOpCount; i++) 
      if (opCount[i]) {
        operation* op = operation::getOpWithIndex(i);
        op->removeStaleComputeTableEntries();
      }
  }
}

void MEDDLY::forest::removeAllComputeTableEntries()
{
  if (is_marked_for_deletion) return;
  if (operation::useMonolithicComputeTable()) {
    is_marked_for_deletion = true;
    operation::removeStalesFromMonolithic();
    is_marked_for_deletion = false;
  } else {
    for (int i=0; i<szOpCount; i++) 
      if (opCount[i]) {
        operation* op = operation::getOpWithIndex(i);
        op->removeAllComputeTableEntries();
      }
  }
}

void MEDDLY::forest::showComputeTable(FILE* s, int verbLevel) const
{
  if (operation::useMonolithicComputeTable()) {
    operation::showMonolithicComputeTable(s, verbLevel);
  } else {
    for (int i=0; i<szOpCount; i++) 
      if (opCount[i]) {
        operation* op = operation::getOpWithIndex(i);
        op->showComputeTable(s, verbLevel);
      }
  }
}

void MEDDLY::forest::registerOperation(const operation* op)
{
  MEDDLY_DCASSERT(op->getIndex() >= 0);
  if (op->getIndex() >= szOpCount) {
    // need to expand
    int newSize = ((op->getIndex() / 16) +1 )*16; // expand in chunks of 16
    int* tmp = (int*) realloc(opCount, newSize * sizeof(int));
    if (0==tmp) throw error(error::INSUFFICIENT_MEMORY);
    for ( ; szOpCount < newSize; szOpCount++) {
      tmp[szOpCount] = 0;
    }
    opCount = tmp;
  }
  opCount[op->getIndex()] ++;
}

void MEDDLY::forest::unregisterOperation(const operation* op)
{
  MEDDLY_DCASSERT(op->getIndex() >= 0);
  MEDDLY_DCASSERT(szOpCount > op->getIndex());
  MEDDLY_DCASSERT(opCount[op->getIndex()]>0);
  opCount[op->getIndex()] --;
}

void MEDDLY::forest::registerEdge(dd_edge& e) 
{
  // add to collection of edges for this forest.
  // change e.index to help find this edge at a later time.
  if (firstHole >= 0) {
    // hole available; fill it up
    int index = firstHole;
    firstHole = edge[firstHole].nextHole;
    edge[index].edge = &e;
    edge[index].nextHole = -1;
    e.setIndex(index);
  } else {
    // no holes available, add to end of array
    if (firstFree >= sz) {
      // expand edge[]
      int new_sz = sz * 2;
      edge_data* new_edge =
          (edge_data*) realloc(edge, new_sz * sizeof(edge_data));
      if (0 == new_edge) throw error(error::INSUFFICIENT_MEMORY);
      edge = new_edge;
      for (int i = sz; i < new_sz; ++i)
      {
        edge[i].nextHole = -1;
        edge[i].edge = 0;
      }
      sz = new_sz;
    }
    MEDDLY_DCASSERT(firstFree < sz);
    edge[firstFree].nextHole = -1;
    edge[firstFree].edge = &e;
    e.setIndex(firstFree);
    ++firstFree;
  }
}


void MEDDLY::forest::unregisterEdge(dd_edge& e) 
{
  // remove this edge from the collection of edges for this forest.
  // change e.index to -1.
  MEDDLY_DCASSERT(e.getIndex() >= 0);
  int index = e.getIndex();
  MEDDLY_DCASSERT(edge[index].edge == &e);
  edge[index].edge = 0;
  edge[index].nextHole = firstHole;
  firstHole = index;
  e.setIndex(-1);
}


void MEDDLY::forest::unregisterDDEdges() 
{
  // Go through the list of valid edges (value > 0), and set
  // the e.index to -1 (indicating unregistered edge).

  // ignore the NULLs; release the rest
  for (unsigned i = 0; i < firstFree; ++i) {
    if (edge[i].edge != 0) {
      MEDDLY_DCASSERT(edge[i].nextHole == -1);
      edge[i].edge->set(0, 0);
      edge[i].edge->setIndex(-1);
    }
  }

  // firstHole < 0 indicates no holes.
  for (unsigned i = 0; i < firstFree; ++i) {
    edge[i].nextHole = -1;
    edge[i].edge = 0;
  }
  firstHole = -1;
  firstFree = 0;
}

// ******************************************************************

MEDDLY::forest::edge_visitor::edge_visitor()
{
}

MEDDLY::forest::edge_visitor::~edge_visitor()
{
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                     expert_forest  methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************



MEDDLY::expert_forest::expert_forest(int ds, domain *d, bool rel, range_type t,
  edge_labeling ev, const policies &p)
: forest(ds, d, rel, t, ev, p)
{
  //
  // Inltialize address array
  //
  a_size = a_min_size;
  address = (node_header *) malloc(a_size * sizeof(node_header));
  if (0 == address) throw error(error::INSUFFICIENT_MEMORY);
  stats.incMemAlloc(a_size * sizeof(node_header));
  memset(address, 0, a_size * sizeof(node_header));
  a_last = a_unused = a_next_shrink = 0;
  // peak_nodes 0;
  

  //
  // Initialize level array
  //
  int N = getNumVariables();
#ifdef NODE_STORAGE_PER_LEVEL
  if (rel) {
    raw_levels = new node_storage*[2*N+1];
    levels = raw_levels + N;
    stats.incMemAlloc(2*N+1 * sizeof(void*));
  } else {
    raw_levels = new node_storage*[N+1];
    levels = raw_levels;
    stats.incMemAlloc(N+1 * sizeof(void*));
  }
#endif

  //
  // Initialize builders array
  //
  if (rel) {
    raw_builders = new node_builder[2*N+1];
    builders = raw_builders + N;
  } else {
    raw_builders = new node_builder[N+1];
    builders = raw_builders;
  }

  //
  // Initialize misc. protected data
  //
  terminalNodesAreStale = false;

  //
  // Initialize misc. private data
  //
  unique = new unique_table(this);
  performing_gc = false;
  in_validate = 0;
  in_val_size = 0;

  //
  // Initialize node characteristics to defaults
  //
  edge_bytes = 0;
  hash_edge_values = false;
  unhashed_bytes = 0;
  hashed_bytes = 0;
}


MEDDLY::expert_forest::~expert_forest() 
{
#ifdef REPORT_ON_DESTROY
  printf("Destroying forest.  Stats:\n");
  reportMemoryUsage(stdout, "\t", 9);
#endif
  // Address array
  free(address);

#ifdef NODE_STORAGE_PER_LEVEL
  // Level array
  for (int k=getMinLevelIndex(); k<=getNumVariables(); k++) {
    delete levels[k];
  }
  delete[] raw_levels;
#else
  delete nodeMan;
#endif

  // builders array
  delete[] raw_builders;

  // unique table
  delete unique;

  // Misc. private data
  free(in_validate);
}

void MEDDLY::expert_forest::initializeForest()
{
  if (!deflt.nodestor) {
    throw error(error::MISCELLANEOUS);
  }

  //
  // Initialize node storage
  //
#ifdef NODE_STORAGE_PER_LEVEL
  for (int k=getMinLevelIndex(); k<=getNumVariables(); k++) {
    levels[k] = k ? deflt.nodestor->createForForest(this) : 0;
  }
#else
  nodeMan = deflt.nodestor->createForForest(this);
#endif

  //
  // Initialize builders array
  //
  for (int k=getMinLevelIndex(); k<=getNumVariables(); k++) {
    builders[k].init(k, this);
  }
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                   public debugging   methods                   '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

void MEDDLY::expert_forest::dump(FILE *s) const
{
  int nwidth = digits(a_last);
  for (long p=0; p<=a_last; p++) {
    fprintf(s, "%*ld\t", nwidth, p);
    showNode(s, p, 1);
    fprintf(s, "\n");
    fflush(s);
  }
}

void MEDDLY::expert_forest::dumpInternal(FILE *s) const
{
  fprintf(s, "Internal forest storage\n");
  fprintf(s, "First unused node index: %ld\n", long(a_unused));
  int awidth = digits(a_last);
  fprintf(s, " Node# :  ");
  for (node_handle p=1; p<=a_last; p++) {
    if (p) fprintf(s, " ");
    fprintf(s, "%*ld", awidth, long(p));
  }
  fprintf(s, "\nLevel  : [");
  for (node_handle p=1; p<=a_last; p++) {
    if (p) fprintf(s, "|");
    fprintf(s, "%*d", awidth, address[p].level);
  }
  fprintf(s, "]\n");
  fprintf(s, "\nOffset : [");
  for (node_handle p=1; p<=a_last; p++) {
    if (p) fprintf(s, "|");
    fprintf(s, "%*ld", awidth, address[p].offset);
  }
  fprintf(s, "]\n");
  fprintf(s, "\nCache  : [");
  for (node_handle p=1; p<=a_last; p++) {
    if (p) fprintf(s, "|");
    fprintf(s, "%*d", awidth, address[p].cache_count);
  }
  fprintf(s, "]\n\n");

#ifdef NODE_STORAGE_PER_LEVEL
  for (int i=1; i<=getNumVariables(); i++) {
    levels[i]->dumpInternal(s);
    if (isForRelations()) levels[-i]->dumpInternal(s);
  }
#else
  nodeMan->dumpInternal(s);
#endif
 
  unique->show(s);
  fflush(s);
}

void MEDDLY::expert_forest::dumpUniqueTable(FILE *s) const
{
  unique->show(s);
}

void MEDDLY::expert_forest::validateIncounts(bool exact)
{
#ifdef VALIDATE_INCOUNTS
  static int idnum = 0;
  idnum++;

  // Inspect every active node's down pointers to determine
  // the incoming count for every active node.
  
  node_handle sz = getLastNode() + 1;
  if (sz > in_val_size) {
    in_validate = (node_handle*) 
                  realloc(in_validate, a_size * sizeof(node_handle));
    in_val_size = a_size;
  }
  MEDDLY_DCASSERT(sz <= in_val_size);
  memset(in_validate, 0, sizeof(node_handle) * sz);
  node_reader P;
  for (node_handle i = 1; i < sz; ++i) {
    MEDDLY_DCASSERT(!isTerminalNode(i));
    if (!isActiveNode(i)) continue;
    initNodeReader(P, i, false);

    // add to reference counts
    for (int z=0; z<P.getNNZs(); z++) {
      if (isTerminalNode(P.d(z))) continue;
      MEDDLY_CHECK_RANGE(0, P.d(z), sz);
      in_validate[P.d(z)]++;
    }
  } // for i

  // Add counts for registered dd_edges
  nodecounter foo(this, in_validate);
  visitRegisteredEdges(foo);
  
  // Validate the incoming count stored with each active node using the
  // in_count array computed above
  for (node_handle i = 1; i < sz; ++i) {
    MEDDLY_DCASSERT(!isTerminalNode(i));
    if (!isActiveNode(i)) continue;
    bool fail = exact
      ?  in_validate[i] != readInCount(i)
      :  in_validate[i] >  readInCount(i);
    if (fail) {
      printf("Validation #%d failed\n", idnum);
      printf("For node %d\n    we got %d\n    node says %d\n",
        i, in_validate[i], readInCount(i));
      dump(stdout);
      throw error(error::MISCELLANEOUS);
    }
    // Note - might not be exactly equal
    // because there could be dd_edges that refer to nodes
    // and we didn't count them.
  }

#ifdef TRACK_DELETIONS
  printf("Incounts validated #%d\n", idnum);
#endif
#endif
}


// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                      public handy methods                      '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

MEDDLY::node_handle* 
MEDDLY::expert_forest::markNodesInSubgraph(node_handle root, bool sort) const
{
  if (isTerminalNode(root)) return 0;

  // initialize lists
  bool* inList = new bool[a_last];
  for (int i=0; i<a_last; i++) inList[i] = false;
  inList--;
  int mlen = 0;
  int msize = 1024;
  node_handle* marked = (node_handle*) malloc(msize * sizeof(node_handle));

  // Initialize search
  marked[mlen] = root;
  mlen++;
  inList[root] = true;

  node_reader *M = node_reader::useReader();

  // Breadth-first search
  for (int mexpl=0; mexpl<mlen; mexpl++) {
    // explore node marked[mexpl]
    initNodeReader(*M, marked[mexpl], false);
    for (int i=0; i<M->getNNZs(); i++) {
      if (isTerminalNode(M->d(i))) continue;
      MEDDLY_CHECK_RANGE(0, M->d(i)-1, a_last);
      if (inList[M->d(i)]) continue;
      // add dn to list
      if (mlen+1 >= msize) { 
          // expand.  Note we're leaving an extra slot
          // at the end, for the terminal 0.
          msize += 1024;
          marked = (node_handle*) realloc(marked, msize*sizeof(node_handle));
          if (0==marked) throw error(error::INSUFFICIENT_MEMORY);
      }
      inList[M->d(i)] = true;
      marked[mlen] = M->d(i);
      mlen++;
    } // for i
  } // for mexpl
  node_reader::recycle(M);

  // sort
  if (sort && mlen>0) {
    mlen = 0;
    for (int i=1; i<=a_last; i++) if (inList[i]) {
      marked[mlen] = i;
      mlen++;
    }
  }

  // cleanup
  inList++;
  delete[] inList;
  if (0==mlen) {
    free(marked);
    return 0;
  }
  // add 0 to the list
  marked[mlen] = 0;
  return marked;
}

long MEDDLY::expert_forest::getNodeCount(node_handle p) const
{
  node_handle* list = markNodesInSubgraph(p, true);
  if (0==list) return 0;
  long i;
  for (i=0; list[i]; i++) { }
  free(list);
  return i;
}

long MEDDLY::expert_forest::getEdgeCount(node_handle p, bool countZeroes) const
{
  node_handle* list = markNodesInSubgraph(p, true);
  if (0==list) return 0;
  long ec=0;
  node_reader *M = node_reader::useReader();
  for (long i=0; list[i]; i++) {
    initNodeReader(*M, list[i], countZeroes);
    ec += countZeroes ? M->getSize() : M->getNNZs();
  }
  node_reader::recycle(M);
  free(list);
  return ec;
}

void MEDDLY::expert_forest::showNode(FILE* s, node_handle p, int verbose) const
{
  if (isTerminalNode(p)) {
    fprintf(s, "(terminal)");
    return;
  }
  if (isDeletedNode(p)) {
    fprintf(s, "DELETED");
    return;
  }
  if (isZombieNode(p)) {
    fprintf(s, "Zombie cc: %d", -address[p].cache_count);
    return;
  }
  const node_header& node = getNode(p);
#ifdef NODE_STORAGE_PER_LEVEL
  const node_storage* nodeMan = levels[node.level];
#endif
  if (verbose) {
    // node: was already written.
    const variable* v = getDomain()->getVar(ABS(node.level));
    if (v->getName()) {
      fprintf(s, " level: %s", v->getName());
    } else {
      fprintf(s, " level: %d", ABS(node.level));
    }
    if (getNodeLevel(p) < 0)
      fprintf(s, "'");
    else
      fprintf(s, " ");
    fprintf(s, " in: %ld", long(nodeMan->getCountOf(node.offset)));
    fprintf(s, " cc: %d", node.cache_count);
  } else {
    fprintf(s, "node: %ld", long(p));
  }
  nodeMan->showNode(s, node.offset, verbose);
}

void MEDDLY::expert_forest::showNodeGraph(FILE *s, node_handle p) const
{
  node_handle* list = markNodesInSubgraph(p, true);
  if (0==list) return;

  // Print by levels
  for (int k = getNumVariables(); k; )
  {
    bool printed = false;
    for (long i=0; list[i]; i++) {
      if (getNodeLevel(list[i]) != k) continue;

      if (!printed) {
        const variable* v = getDomain()->getVar(ABS(k));
        char primed = (k>0) ? ' ' : '\'';
        if (v->getName()) {
          fprintf(s, "Level: %s%c\n", v->getName(), primed);
        } else {
          fprintf(s, "Level: %d%c\n", ABS(k), primed);
        }
        printed = true;
      }

      fprintf(s, "  ");
      showNode(s, list[i]);
      fprintf(s, "\n");
    }
    
    // next level
    k *= -1;
    if (k>0) k--;
  } // for k

  free(list);
}

void MEDDLY::expert_forest
::reportMemoryUsage(FILE * s, const char* pad, int verb) const
{
  if (verb>0) {
    fprintf(s, "%sPeak Nodes:             %ld\n", pad, getPeakNumNodes());
    fprintf(s, "%sActive Nodes:           %ld\n", pad, getCurrentNumNodes());
  }
#if 0
  unsigned count = 0;
  for (int i = 1; i <= getLastNode(); ++i) if (isActiveNode(i)) ++count;
  fprintf(s, "%cActive Nodes (manual):\t\t%d\n", pad, count);
  fprintf(s, "%c%cZombie Nodes:\t\t%d\n", pad, filler,
      getZombieNodeCount());
  fprintf(s, "%c%cTemp Nodes:\t\t%d\n", pad, filler, getTempNodeCount());
  fprintf(s, "%c%cOrphan Nodes:\t\t%d\n", pad, filler,
      getOrphanNodeCount());
#endif
  if (verb>2) {
    fprintf(s, "%sReclaimed Nodes:        %ld\n", pad, stats.reclaimed_nodes);
  }
  fprintf(s, "%sMem Used:               %ld\n", pad, getCurrentMemoryUsed());
  fprintf(s, "%sPeak Mem Used:          %ld\n", pad, getPeakMemoryUsed());
  if (verb>1) {
    fprintf(s, "%sMem Allocated:          %ld\n", pad,
      getCurrentMemoryAllocated());
    fprintf(s, "%sPeak Mem Allocated:     %ld\n",
      pad, getPeakMemoryAllocated());
  }
  if (verb>3) {
    fprintf(s, "%sUnique Tbl Mem Used:    %u\n", pad, 
      unique->getMemUsed());
  }
  if (verb>5) {
    fprintf(s, "%sCompactions:            %ld\n", pad, stats.num_compactions);
    fprintf(s, "%sGarbage Collections:    %ld\n", pad, stats.garbage_collections);
  }
#ifdef NODE_STORAGE_PER_LEVEL
  for (int i=getMinLevelIndex(); i<=getNumVariables(); i++) {
    if (levels[i]) levels[i]->reportMemoryUsage(s, pad, verb);
  }
#else
  nodeMan->reportMemoryUsage(s, pad, verb);
#endif
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                   methods for  reading nodes                   '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

void MEDDLY::expert_forest
::initRedundantReader(node_reader &nr, int k, node_handle p, bool full) const
{
  MEDDLY_DCASSERT(0==edgeBytes());
  int nsize = getLevelSize(k);
  nr.resize(k, nsize, 0, full);
  for (int i=0; i<nsize; i++) {
    nr.down[i] = p;
  }
  if (!full) {
    for (int i=0; i<nsize; i++) nr.index[i] = i;
    nr.nnzs = nsize;
  }
}

void MEDDLY::expert_forest
::initRedundantReader(node_reader &nr, int k, int ev, node_handle np, bool full) const
{
  MEDDLY_DCASSERT(sizeof(int)==edgeBytes());
  int nsize = getLevelSize(k);
  nr.resize(k, nsize, sizeof(int), full);
  for (int i=0; i<nsize; i++) {
    nr.down[i] = np;
    ((int*)nr.edge)[i] = ev;
  }
  if (!full) {
    for (int i=0; i<nsize; i++) nr.index[i] = i;
    nr.nnzs = nsize;
  }
}

void MEDDLY::expert_forest
::initIdentityReader(node_reader &nr, int k, int i, node_handle node, bool full) const
{
  MEDDLY_DCASSERT(0==edgeBytes());
  int nsize = getLevelSize(k);
  if (full) {
    nr.resize(k, nsize, 0, full);
    memset(nr.down, 0, nsize * sizeof(node_handle));
    nr.down[i] = node;
  } else {
    nr.resize(k, 1, 0, full);
    nr.nnzs = 1;
    nr.down[0] = node;
    nr.index[0] = i;
  }
}


void MEDDLY::expert_forest
::initIdentityReader(node_reader &nr, int k, int i, int ev, node_handle node, 
  bool full) const
{
  MEDDLY_DCASSERT(sizeof(int)==edgeBytes());
  int nsize = getLevelSize(k);
  if (full) {
    nr.resize(k, nsize, sizeof(int), full);
    memset(nr.down, 0, nsize * sizeof(node_handle));
    memset(nr.edge, 0, nsize * sizeof(int));
    nr.down[i] = node;
    ((int*)nr.edge)[i] = ev;
  } else {
    nr.resize(k, 1, sizeof(int), full);
    nr.nnzs = 1;
    nr.down[0] = node;
    ((int*)nr.edge)[0] = ev;
    nr.index[0] = i;
  }
}



// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '   public virtual methods provided here,  no need to override   '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


void MEDDLY::expert_forest::garbageCollect()
{
  if (performing_gc) return;
  performing_gc = true;
  stats.garbage_collections++;

#ifdef DEBUG_GC
  printf("Garbage collection in progress... \n");
  fflush(stdout);
#endif

  if (isPessimistic()) {
#ifdef DEBUG_GC
    printf("Zombie nodes: %ld\n", stats.zombie_nodes);
#endif
    // remove the stale nodes entries from caches
    removeStaleComputeTableEntries();
#ifdef DEBUG_GC
    printf("Zombie nodes: %ld\n", stats.zombie_nodes);
#endif
#ifdef DEVELOPMENT_CODE
    if (stats.zombie_nodes != 0) {
      showInfo(stderr, 1);
      showComputeTable(stderr, 6);
    }
    MEDDLY_DCASSERT(stats.zombie_nodes == 0);
#endif
  } else {
    // remove the stale nodes entries from caches
    removeStaleComputeTableEntries();
  }

  if (deflt.compactAfterGC) {
#ifdef DEBUG_GC
    printf("Compacting levels...\n");
    fflush(stdout);
#endif

    compactMemory(); 

#ifdef DEBUG_GC
    printf("  done compacting.\n");
    fflush(stdout);
#endif
  }

  performing_gc = false;
}

void MEDDLY::expert_forest::compactMemory()
{
#ifdef NODE_STORAGE_PER_LEVEL
  for (int i=getMinLevelIndex(); i<=getNumVariables(); i++) {
    levels[i]->collectGarbage(true);
  }
#else 
  nodeMan->collectGarbage(true);
#endif
}

void MEDDLY::expert_forest::showInfo(FILE* s, int verb)
{
  // Show forest with appropriate level of detail
  if (1==verb)  dump(s);
  else          dumpInternal(s); 
  fprintf(s, "DD stats:\n");
  reportMemoryUsage(s, "    ", verb);
  fprintf(s, "Unique table stats:\n");
  fprintf(s, "\t%-24s%u\n", "Current size:", unique->getSize());
  fprintf(s, "\t%-24s%u\n", "Current entries:", unique->getNumEntries());
}


// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '    virtual methods to be overridden by some derived classes    '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

bool MEDDLY::expert_forest
::areEdgeValuesEqual(const void* eva, const void* evb) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::normalize(node_builder &nb, int& ev) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::normalize(node_builder &nb, float& ev) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::showTerminal(FILE* s, node_handle tnode) const
{
  throw error(error::NOT_IMPLEMENTED);
}

void MEDDLY::expert_forest::showEdgeValue(FILE* s, const void* edge, int i) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::showHashedHeader(FILE* s, const void* hh) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::showUnhashedHeader(FILE* s, const void* uh) const
{
  throw error(error::TYPE_MISMATCH);
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                        private  methods                        '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

void MEDDLY::expert_forest::handleNewOrphanNode(node_handle p)
{
  MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));
  MEDDLY_DCASSERT(isActiveNode(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(readInCount(p) == 0);

  // insted of orphan_nodes++ here; do it only when the orphan is not going
  // to get deleted or converted into a zombie

  // Two possible scenarios:
  // (1) a reduced node, or
  // (2) a temporary node ready to be deleted.
  // MEDDLY_DCASSERT(isReducedNode(p) || getCacheCount(p) == 0);

  if (getCacheCount(p) == 0) {
    // delete node
    // this should take care of the temporary nodes also
#ifdef TRACK_DELETIONS
    printf("Deleting node %d\n\t", p);
    showNode(stdout, p);
    printf("\n");
    fflush(stdout);
#endif
    deleteNode(p);
  }
  else if (isPessimistic()) {
    // zombify node
#ifdef TRACK_DELETIONS
    printf("Zombifying node %d\n", p);
    fflush(stdout);
#endif
    zombifyNode(p);
  }
  else {
    stats.orphan_nodes++;
  }

}

void MEDDLY::expert_forest::deleteOrphanNode(node_handle p) 
{
  MEDDLY_DCASSERT(!isPessimistic());
  MEDDLY_DCASSERT(getCacheCount(p) == 0 && readInCount(p) == 0);
#ifdef TRACK_DELETIONS
  printf("Deleting node %d\n\t", p);
  showNode(stdout, p);
  printf("\n");
  fflush(stdout);
#endif
  stats.orphan_nodes--;
  deleteNode(p);
}

void MEDDLY::expert_forest::deleteNode(node_handle p)
{
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(readInCount(p) == 0);
  MEDDLY_DCASSERT(isActiveNode(p));

#ifdef VALIADE_INCOUNTS_ON_DELETE
  validateIncounts(false);
#endif

  unsigned h = hashNode(p);
#ifdef DEVELOPMENT_CODE
  node_reader key; 
  initNodeReader(key, p, false);
  key.computeHash(areEdgeValuesHashed());
  if (unique->find(key) != p) {
    fprintf(stderr, "Error in deleteNode\nFind: %ld\np: %ld\n",
      long(unique->find(key)), long(p));
    dumpInternal(stdout);
    MEDDLY_DCASSERT(false);
  }
  node_handle x = unique->remove(h, p);
  MEDDLY_DCASSERT(p == x);
#else
  unique->remove(h, p);
#endif

#ifdef TRACK_DELETIONS
  printf("%s: p = %d, unique->remove(p) = %d\n", __func__, p, x);
  fflush(stdout);
#endif

  MEDDLY_DCASSERT(address[p].cache_count == 0);

  node_address addr = getNode(p).offset;

#ifdef NODE_STORAGE_PER_LEVEL
  node_storage* nodeMan = levels[getNode(p).level];
#endif

  // unlink children and recycle node memory
  nodeMan->unlinkDownAndRecycle(addr);

  // recycle the index
  freeActiveNode(p);

  // if (nodeMan.compactLevel) nodeMan.compact(false);

#ifdef VALIDATE_INCOUNTS_ON_DELETE
  validateIncounts(false);
#endif

}

void MEDDLY::expert_forest::zombifyNode(node_handle p)
{
  MEDDLY_DCASSERT(isActiveNode(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(getCacheCount(p) > 0);  // otherwise this node should be deleted
  MEDDLY_DCASSERT(readInCount(p) == 0);
  MEDDLY_DCASSERT(address[p].cache_count > 0);

  stats.zombie_nodes++;
  stats.decActive(1);

  // mark node as zombie
  address[p].cache_count = -address[p].cache_count;

  unsigned h = hashNode(p);
#ifdef DEVELOPMENT_CODE 
  node_reader key; 
  initNodeReader(key, p, false);
  key.computeHash(areEdgeValuesHashed());
  if (unique->find(key) != p) {
    fprintf(stderr, "Fail: can't find reduced node %ld; got %ld\n", 
      long(p), long(unique->find(key)));
    dumpInternal(stderr);
    MEDDLY_DCASSERT(false);
  }
  node_handle x = unique->remove(h, p);
  MEDDLY_DCASSERT(x==p);
#else
  unique->remove(h, p);
#endif

  node_address addr = getNode(p).offset;

#ifdef NODE_STORAGE_PER_LEVEL
  node_storage &nodeMan = levels[getNode(p).level];
#endif

  // unlink children and recycle node memory
  nodeMan->unlinkDownAndRecycle(addr);
}

MEDDLY::node_handle MEDDLY::expert_forest
::createReducedHelper(int in, const node_builder &nb)
{
#ifdef DEVELOPMENT_CODE
  validateDownPointers(nb);
#endif

  // get sparse, truncated full sizes and check
  // for redundant / identity reductions.
  int nnz;
  if (nb.isSparse()) {
    // Reductions for sparse nodes
    nnz = nb.getNNZs();
#ifdef DEVELOPMENT_CODE
    for (int z=0; z<nnz; z++) {
      MEDDLY_DCASSERT(nb.d(z));
    } // for z
#endif

    // Check for identity nodes
    if (1==nnz && in==nb.i(0)) {
      if (isIdentityEdge(nb, 0)) {
        return nb.d(0);
      }
    }

    // Check for redundant nodes
    if (nnz == getLevelSize(nb.getLevel())) {
      if (isRedundant(nb)) {
        // unlink downward pointers, except the one we're returning.
        for (int i = 1; i<nb.getNNZs(); i++)  unlinkNode(nb.d(i));
        return nb.d(0);
      }
    }

  } else {
    // Reductions for full nodes
    MEDDLY_DCASSERT(nb.getSize() == getLevelSize(nb.getLevel()));
    nnz = 0;
    for (int i=0; i<nb.getSize(); i++) {
      if (nb.d(i)) nnz++;
    } // for i

    // Check for identity nodes
    if (1==nnz) {
      if (isIdentityEdge(nb, in)) {
        return nb.d(in);
      }
    }

    // Check for redundant nodes
    if (nnz == nb.getSize()) {
      if (isRedundant(nb)) {
        // unlink downward pointers, except the one we're returning.
        for (int i = 1; i<nb.getSize(); i++)  unlinkNode(nb.d(i));
        return nb.d(0);
      }
    }
  }

  // Is this a zero node?
  if (0==nnz) {
    // no need to unlink
    return 0;
  }

  // check for duplicates in unique table
  node_handle q = unique->find(nb);
  if (q) {
    // unlink all downward pointers
    for (int i = 0; i<nb.rawSize(); i++)  unlinkNode(nb.d(i));
    return linkNode(q);
  }

  // 
  // Not eliminated by reduction rule.
  // Not a duplicate.
  //
  // We need to create a new node for this.

  // NOW is the best time to run the garbage collector, if necessary.
#ifndef GC_OFF
  if (isTimeToGc()) garbageCollect();
#endif

  // Grab a new node
  node_handle p = getFreeNodeHandle();
  address[p].level = nb.getLevel();
  MEDDLY_DCASSERT(0 == address[p].cache_count);
#ifdef NODE_STORAGE_PER_LEVEL
  node_storage* nodeMan = levels[nb.getLevel()];
#endif

  // All of the work is in nodeMan now :^)
  address[p].offset = nodeMan->makeNode(p, nb, getNodeStorage());

#ifdef SAVE_HASHES
  address[p].hash = nb.hash();
#endif
  
  // add to UT 
  unique->add(nb.hash(), p);

#ifdef DEVELOPMENT_CODE
  node_reader key;
  initNodeReader(key, p, false);
  key.computeHash(areEdgeValuesHashed());
  MEDDLY_DCASSERT(key.hash() == nb.hash());
  MEDDLY_DCASSERT(unique->find(key) == p);
#endif
#ifdef DEBUG_CREATE_REDUCED
  printf("Created node %d\n", p);
  dump(stdout);
#endif
  return p;
}

void MEDDLY::expert_forest::validateDownPointers(const node_builder &nb) const
{
  int nextLevel;
  switch (getReductionRule()) {
    case policies::IDENTITY_REDUCED:
    case policies::FULLY_REDUCED:
      if (nb.isSparse()) {
        for (int z=0; z<nb.getNNZs(); z++) {
          MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), getNodeLevel(nb.d(z))));
        } 
      } else {
        for (int i=0; i<nb.getSize(); i++) {
          MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), getNodeLevel(nb.d(i))));
        }
      }
      break;

    case policies::QUASI_REDUCED:
      if (isForRelations()) 
        nextLevel = (nb.getLevel()<0) ? -(nb.getLevel()+1) : -nb.getLevel();
      else
        nextLevel = nb.getLevel()-1;
      if (nb.isSparse()) {
        for (int z=0; z<nb.getNNZs(); z++) {
          MEDDLY_DCASSERT(getNodeLevel(nb.d(z)) == nextLevel);
        } 
      } else {
        for (int i=0; i<nb.getSize(); i++) {
          if (0==nb.d(i)) continue;
          MEDDLY_DCASSERT(getNodeLevel(nb.d(i)) == nextLevel);
        }
      }
      break;

    default:
      throw error(error::NOT_IMPLEMENTED);
  }

}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '          Methods for managing  available node handles          '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

void MEDDLY::expert_forest::expandHandleList()
{
  // increase size by 50%
  int delta = a_size / 2;
  MEDDLY_DCASSERT(delta>=0);
  address = (node_header*) realloc(address, (a_size+delta) * sizeof(node_header));
  if (0==address) {
    throw error(error::INSUFFICIENT_MEMORY);
  }
  stats.incMemAlloc(delta * sizeof(node_header));
  memset(address + a_size, 0, delta * sizeof(node_header));
  a_size += delta;
  a_next_shrink = a_size / 2;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged address array, new size %d\n", a_size);
#endif
}

void MEDDLY::expert_forest::shrinkHandleList()
{
  // Determine new size
  int new_size = a_min_size;
  while (a_last >= new_size) new_size += new_size/2;
  int delta = a_size - new_size;
  if (0==delta) {
    a_next_shrink = 0;
    return;
  }
  // shrink the array
  MEDDLY_DCASSERT(delta>=0);
  MEDDLY_DCASSERT(a_size-delta>=a_min_size);
  address = (node_header*) realloc(address, new_size * sizeof(node_header));
  if (0==address) {
    throw error(error::INSUFFICIENT_MEMORY);
  }
  stats.decMemAlloc(delta * sizeof(node_header));
  a_size -= delta;
  a_next_shrink = a_size / 2;
  MEDDLY_DCASSERT(a_last < a_size);
  // rebuild the free list
  a_unused = 0;
  for (int i=a_last; i; i--) {
    if (address[i].isDeleted()) {
      address[i].setNextDeleted(a_unused);
      a_unused = i;
    }
  }
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Shrank address array, new size %d\n", a_size);
#endif
}

// ******************************************************************
// *                                                                *
// *               expert_forest::nodecounter methods               *
// *                                                                *
// ******************************************************************

MEDDLY::expert_forest::nodecounter::nodecounter(expert_forest *p, int* c)
 : edge_visitor()
{
  parent = p;
  counts = c;
}

MEDDLY::expert_forest::nodecounter::~nodecounter()
{
  // DO NOT delete counts.
}

void MEDDLY::expert_forest::nodecounter::visit(dd_edge &e)
{
  int n = e.getNode();
  if (parent->isTerminalNode(n)) return;
  MEDDLY_DCASSERT(n>0);
  MEDDLY_DCASSERT(n<=parent->getLastNode());
  counts[n]++;
}


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

#include <fstream>
#include <sstream>
#include "defines.h"
#include "initializer.h"
#include "unique_table.h"
#include "relation_node.h"
#include "impl_unique_table.h"
#include "hash_stream.h"
// #include "storage/bytepack.h"
#include "reordering/reordering_factory.h"

#include "oper.h"
#include "operators.h"
#include "node_marker.h"

// for timestamps.
// to do - check during configuration that these are present,
// and act accordingly here

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

//
// For constructing forests.
//
#include "forests/mtmddbool.h"
#include "forests/mtmddint.h"
#include "forests/mtmddreal.h"

#include "forests/mtmxdbool.h"
#include "forests/mtmxdint.h"
#include "forests/mtmxdreal.h"

#include "forests/evmdd_pluslong.h"
#include "forests/evmdd_timesreal.h"

#include "forests/evmxd_pluslong.h"
#include "forests/evmxd_timesreal.h"


#ifdef ALLOW_DEPRECATED_0_17_3
#include "io_dot.h"
#endif

//#include <set>
//#include <queue>
//#include <vector>

// #define DEBUG_CLEANUP
// #define DEBUG_ADDRESS_RESIZE
// #define DEBUG_CREATE_REDUCED
// #define DEBUG_GC
// #define DEBUG_WRITE
// #define DEBUG_READ

// #define TRACK_DELETIONS

// Thoroughly check reference counts.
// Very slow.  Use only for debugging.
// #define VALIDATE_INCOUNTS
// #define VALIDATE_INCOUNTS_ON_DELETE

// #define SHOW_VALIDATE_CACHECOUNTS

// #define GC_OFF

// #define REPORT_ON_DESTROY
// #define DUMP_ON_FOREST_DESTROY

// ******************************************************************
// *                                                                *
// *                                                                *
// *                          forest class                          *
// *                                                                *
// *                                                                *
// ******************************************************************

//
// Static members
//

MEDDLY::forest** MEDDLY::forest::all_forests;
unsigned MEDDLY::forest::max_forests;
unsigned MEDDLY::forest::gfid;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Static "constructors"
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MEDDLY::forest*
MEDDLY::forest::create(domain* d, set_or_rel sr, range_type t,
            edge_labeling el, const policies &p,
            int* level_reduction_rule, int tv)
{
    switch (el) {
        case edge_labeling::MULTI_TERMINAL:
            switch (t) {
                case range_type::BOOLEAN:
                    if (sr)
                        return new mt_mxd_bool(d, p, level_reduction_rule, tv);
                    else
                        return new mt_mdd_bool(d, p, level_reduction_rule, tv);

                case range_type::INTEGER:
                    if (sr)
                        return new mt_mxd_int(d, p, level_reduction_rule, tv);
                    else
                        return new mt_mdd_int(d, p, level_reduction_rule, tv);

                case range_type::REAL:
                    if (sr)
                        return new mt_mxd_real(d, p, level_reduction_rule, (float)tv);
                    else
                        return new mt_mdd_real(d, p, level_reduction_rule, (float)tv);

                default:
                    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            }; // range type switch
            // should never get here
            return nullptr;


        case edge_labeling::EVPLUS:
            if (range_type::INTEGER != t) {
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            }
            if (sr)
                return new evmxd_pluslong(d, p, level_reduction_rule);
            else
                return new evmdd_pluslong(d, p, level_reduction_rule);


        case edge_labeling::INDEX_SET:
            if (range_type::INTEGER != t || sr) {
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            }
            return new evmdd_index_set_long(d, p, level_reduction_rule);


        case edge_labeling::EVTIMES:
            if (range_type::REAL != t || !sr) {
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            }
            return new evmxd_timesreal(d, p);


        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    } // edge label switch

    // should never get here
    return nullptr;
}


//
// Forest registry methods
//

void MEDDLY::forest::initStatics()
{
    all_forests = new forest* [256];
    max_forests = 256;
    for (unsigned i=0; i<max_forests; i++) {
        all_forests[i] = nullptr;
    }
    gfid = 0;
}

void MEDDLY::forest::freeStatics()
{
    delete[] all_forests;
    max_forests = 0;
}

void MEDDLY::forest::registerForest(forest* f)
{
    // Assign global ID to f
    gfid++;
    f->fid = gfid;

    // Get slot for f
    if (gfid >= max_forests) {
        // enlarge array
        unsigned newmax = max_forests * 2;
        if (newmax > 1000000000) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
        forest** newall = new forest* [newmax];
        unsigned i;
        for (i=0; i<max_forests; i++) {
            newall[i] = all_forests[i];
        }
        for ( ; i<newmax; i++) {
            newall[i] = nullptr;
        }
        delete[] all_forests;
        all_forests = newall;
        max_forests = newmax;
    }
    all_forests[gfid] = f;

    // Register in the domain
    f->d->registerForest(f);
}

void MEDDLY::forest::unregisterForest(forest* f)
{
    // Remove from forest slot
    if (f->fid <= max_forests) {
        all_forests[f->fid] = nullptr;
    }

    // Unregister in the domain
    f->d->unregisterForest(f);
}

//
// Root edge registry methods
//

void MEDDLY::forest::registerEdge(dd_edge& e)
{
    unsigned index;
    if (roots_hole) {
        // Pull from the free list
        //
        index = roots_hole;
        roots_hole = roots[roots_hole].nextHole;
    } else {
        // Pull from the array.
        // Check if we need to expand
        if (roots_next >= roots_size) {
            unsigned new_size = roots_size * 2;
            if (new_size > 1000000000) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            edge_data* new_roots = new edge_data[new_size];
            unsigned i;
            for (i=0; i<roots_size; i++) {
                new_roots[i] = roots[i];
            }
            for (; i<new_size; i++) {
                new_roots[i].nextHole = 0;
                new_roots[i].edge = nullptr;
            }
            delete[] roots;
            roots = new_roots;
            roots_size = new_size;
        }
        index = roots_next++;
    }

    roots[index].edge = &e;
    roots[index].nextHole = 0;
    e.index = index;
    e.parentFID = fid;
}


void MEDDLY::forest::unregisterEdge(dd_edge& e)
{
    // remove a root edge.
    MEDDLY_DCASSERT(e.parentFID == fid);
    MEDDLY_DCASSERT(e.index > 0);
    MEDDLY_DCASSERT(roots[e.index].edge == &e);

    if (e.index+1 == roots_next) {
        //
        // Instead of adding to the free list,
        // absorb this "hole" at the end of the array
        //
        MEDDLY_DCASSERT(roots_next);
        roots_next--;
        MEDDLY_DCASSERT(roots_next);
        roots[e.index].nextHole = 0;
    } else {
        //
        // Add to the front of the free list
        //
        roots[e.index].nextHole = roots_hole;
        roots_hole = e.index;
    }

    roots[e.index].edge = nullptr;
    e.index = 0;
    e.parentFID = 0;
}


void MEDDLY::forest::unregisterDDEdges()
{
    // Unregister ALL root edges
    // (e.g., because we're destroying the forest)

    MEDDLY_DCASSERT(roots);
    MEDDLY_DCASSERT(0==roots[0].edge);

    // ignore the NULLs; release the rest
    for (unsigned i = 1; i < roots_next; ++i) {
        if (roots[i].edge) {
            MEDDLY_DCASSERT(0==roots[i].nextHole);
            roots[i].edge->set(0);
            roots[i].edge->index = 0;
            roots[i].edge->parentFID = 0;
        }
    }

    roots_hole = 0;
    roots_next = 1;
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *              OLD,  still unorganized forest stuff              *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                    forest::statset  methods                    *
// *                                                                *
// ******************************************************************

MEDDLY::forest::statset::statset()
{
  reachable_scans = 0;
  reclaimed_nodes = 0;
  num_compactions = 0;
  garbage_collections = 0;
#ifdef TRACK_UNREACHABLE_NODES
  unreachable_nodes = 0;
#endif
  active_nodes = 0;
  peak_active = 0;
  /*
  memory_used = 0;
  memory_alloc = 0;
  peak_memory_used = 0;
  peak_memory_alloc = 0;
  */
  /*
  memory_UT = 0;
  peak_memory_UT = 0;
  max_UT_chain = 0;
  */
}

// ******************************************************************
// *                                                                *
// *                     forest::logger methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::forest::logger::logger()
{
  nfix = true;
  /*
    All other settings to false;
    the default logger does nothing!
  */
  node_counts = false;
  time_stamps = false;

  /* Do this now, regardless */
#ifdef HAVE_SYS_TIME_H
  static struct timeval curr_time;
  static struct timezone tz;
  gettimeofday(&curr_time, &tz);
  startsec = curr_time.tv_sec;
  startusec = curr_time.tv_usec;
#else
  startsec = 0;
  startusec = 0;
#endif
}

MEDDLY::forest::logger::~logger()
{
}

void MEDDLY::forest::logger::currentTime(long &sec, long &usec)
{
#ifdef HAVE_SYS_TIME_H
  static struct timeval curr_time;
  static struct timezone tz;
  gettimeofday(&curr_time, &tz);
  sec = curr_time.tv_sec - startsec;
  usec = curr_time.tv_usec - startusec;
  if (usec<0) {
    sec--;
    usec += 1000000;
  }
#else
  sec = 0;
  usec = 0;
#endif
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                         forest methods                         *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::policies MEDDLY::forest::mddDefaults;
MEDDLY::policies MEDDLY::forest::mxdDefaults;


MEDDLY::forest
::forest(domain* _d, bool rel, range_type t, edge_labeling ev,
  const policies &p,int* lrr) : deflt(p)
{
    d = _d;
    isRelation = rel;
    rangeType = t;
    edgeLabel = ev;

    registerForest(this);

#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating forest #%u in domain #%d\n", ds, _d->ID());
#endif


  is_marked_for_deletion = false;


    if(lrr==NULL){

        if(isUserDefinedReduced()){
		throw error(error::INVALID_POLICY, __FILE__, __LINE__);
        }

	else{

        if(isRelation)
        {
            int span = 2*(d->getNumVariables()) + 1;
            lrr=(int*)malloc(sizeof(int)*span);

            if(isQuasiReduced())
                for(int i=0;i<span;i++)
                {lrr[i]=i==0?-4:(i%2==0?-2:-1);}

            else if (isFullyReduced())
                for(int i=0;i<span;i++)
                {lrr[i]=i==0?-4:(i%2==0?-1:-1);}

            else if (isIdentityReduced())
                for(int i=0;i<span;i++)
                {lrr[i]=i==0?-4:(i%2==0?-3:-1);}
        }
        else
        {
            lrr=(int*)malloc(sizeof(int)*(d->getNumVariables()+1));
            lrr[0]=-4;
            if(isQuasiReduced())
                for(int i=1;i<=d->getNumVariables();i++)
                    lrr[i]=-2;

            else
                if (isFullyReduced())
                for(int i=1;i<=d->getNumVariables();i++)
                    lrr[i]=-1;

        }

	}

	level_reduction_rule=lrr;
    }
	else if(isUserDefinedReduced()){

    	level_reduction_rule=lrr;
    }
	else
		throw error(error::INVALID_POLICY, __FILE__, __LINE__);




  // check policies
  if (!isRelation) {
    if (reduction_rule::IDENTITY_REDUCED == deflt.reduction)
      throw error(error::INVALID_POLICY, __FILE__, __LINE__);

    for(int i=1;i<=d->getNumVariables();i++)
        if(level_reduction_rule[i]==-3)              //isIdentityReduced()
         throw error(error::INVALID_POLICY, __FILE__, __LINE__);
  }


    //
    // Initialize array of operations
    //
    opCount = 0;
    szOpCount = 0;


    //
    // Initialize the root edges
    //
    roots_hole = 0;  // firstHole == 0 indicates no holes.
    roots_next = 1;  // never use slot 0
    roots_size = 1024;

    // Create an array to store pointers to dd_edges.
    roots = new edge_data [roots_size];
    for (unsigned i = 0; i < roots_size; ++i) {
        roots[i].nextHole = 0;
        roots[i].edge = nullptr;
    }

    //
    // Empty logger
    //

    theLogger = 0;
}

MEDDLY::forest::~forest()
{
  // operations are deleted elsewhere...
  free(opCount);

  for (unsigned i = 0; i < roots_size; ++i) {
      if (roots[i].edge) {
          roots[i].edge->parentFID = 0;
          roots[i].edge->node = 0;
          roots[i].edge->edgeval.set();
      }
  }
  delete[] roots;

  // NOTE: since the user is provided with the dd_edges instances (as opposed
  // to a pointer), the user program will automatically call the
  // destructor for each dd_edge when the corresponding variable goes out of
  // scope. Therefore there is no need to destruct dd_edges from here.

    unregisterForest(this);
}

void MEDDLY::forest::markForDeletion()
{
  if (is_marked_for_deletion) return;
  is_marked_for_deletion = true;
  // deal with operations associated with this forest
  for (unsigned i=0; i<szOpCount; i++)
    if (opCount[i]) {
      operation* op = operation::getOpWithIndex(i);
      op->markForDeletion();
    }
  unregisterDDEdges();
}

void MEDDLY::forest::createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest::createEdgeForVar(int vh, bool vp, const long* terms, dd_edge& a)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest::createEdgeForVar(int vh, bool vp, const float* terms, dd_edge& a)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

MEDDLY::node_handle MEDDLY::forest::unionOneMinterm(node_handle a,  int* from,  int* to, int level)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}


void MEDDLY::forest::createEdge(const int* const* vlist, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const long* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const float* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest
::createEdge(const int* const* vl, const int* const* vpl, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const int* const* vplist, const long* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest::createEdge(bool val, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest::createEdge(long val, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest::createEdge(float val, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest::evaluate(const dd_edge &f, const int* vl, bool &t) const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest::evaluate(const dd_edge &f, const int* vl, long &t) const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest::evaluate(const dd_edge &f, const int* vl, float &t) const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest
::evaluate(const dd_edge& f, const int* vl, const int* vpl, bool &t) const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest
::evaluate(const dd_edge& f, const int* vl, const int* vpl, long &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::evaluate(const dd_edge& f, const int* vl, const int* vpl, float &t) const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest::getElement(const dd_edge& a, int index, int* e)
{
  throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
}

void MEDDLY::forest::getElement(const dd_edge& a, long index, int* e)
{
  throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
}

//



void MEDDLY::forest::removeStaleComputeTableEntries()
{
  if (operation::usesMonolithicComputeTable()) {
    operation::removeStalesFromMonolithic();
  } else {
    for (unsigned i=0; i<szOpCount; i++)
      if (opCount[i]) {
        operation* op = operation::getOpWithIndex(i);
        op->removeStaleComputeTableEntries();
      }
  }
}

void MEDDLY::forest::removeAllComputeTableEntries()
{
  if (is_marked_for_deletion) return;
  if (operation::usesMonolithicComputeTable()) {
    is_marked_for_deletion = true;
    operation::removeStalesFromMonolithic();
    is_marked_for_deletion = false;
  } else {
    for (unsigned i=0; i<szOpCount; i++)
      if (opCount[i]) {
        operation* op = operation::getOpWithIndex(i);
        op->removeAllComputeTableEntries();
      }
  }
}

void MEDDLY::forest::showComputeTable(output &s, int verbLevel) const
{
  if (operation::usesMonolithicComputeTable()) {
    operation::showMonolithicComputeTable(s, verbLevel);
  } else {
    for (unsigned i=0; i<szOpCount; i++)
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
    unsigned newSize = ((op->getIndex() / 16) +1 )*16; // expand in chunks of 16
    unsigned* tmp = (unsigned*) realloc(opCount, newSize * sizeof(unsigned));
    if (0==tmp) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
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



// ******************************************************************

#ifdef EDGE_VISITORS

MEDDLY::forest::edge_visitor::edge_visitor()
{
}

MEDDLY::forest::edge_visitor::~edge_visitor()
{
}

#endif

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      expert_forest  stuff                      *
// *                                                                *
// *                                                                *
// ******************************************************************

#ifdef EDGE_VISITORS

// ******************************************************************
// *                                                                *
// *                expert_forest::nodecounter class                *
// *                                                                *
// ******************************************************************

class MEDDLY::expert_forest::nodecounter: public edge_visitor {
    expert_forest* parent;
    int* counts;
  public:
    nodecounter(expert_forest*p, int* c);
    virtual ~nodecounter();
    virtual void visit(dd_edge &e);
};


// ******************************************************************
// *               expert_forest::nodecounter methods               *
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

// ******************************************************************
// *                                                                *
// *                expert_forest::nodemarker  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::expert_forest::nodemarker: public edge_visitor {
    expert_forest* parent;
  public:
    nodemarker(expert_forest *p);
    virtual ~nodemarker();
    virtual void visit(dd_edge &e);
};


// ******************************************************************
// *               expert_forest::nodemarker  methods               *
// ******************************************************************

MEDDLY::expert_forest::nodemarker::nodemarker(expert_forest *p)
 : edge_visitor()
{
  parent = p;
}

MEDDLY::expert_forest::nodemarker::~nodemarker()
{
  // nothing to do
}

void MEDDLY::expert_forest::nodemarker::visit(dd_edge &e)
{
  if (e.getForest() != parent) return;
#ifdef DEBUG_MARK_SWEEP
  printf("Traversing root node %ld\n", e.getNode());
#endif
  parent->markNode(e.getNode());
}

#endif


// ******************************************************************
// *                                                                *
// *                                                                *
// *                     expert_forest  methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************

const unsigned MEDDLY::expert_forest::HUMAN_READABLE_MEMORY   = 0x0001;
const unsigned MEDDLY::expert_forest::BASIC_STATS             = 0x0002;
const unsigned MEDDLY::expert_forest::EXTRA_STATS             = 0x0004;
const unsigned MEDDLY::expert_forest::FOREST_STATS            = 0x0008;
const unsigned MEDDLY::expert_forest::STORAGE_STATS           = 0x0010;
const unsigned MEDDLY::expert_forest::STORAGE_DETAILED        = 0x0020;
const unsigned MEDDLY::expert_forest::UNIQUE_TABLE_STATS      = 0x0040;
const unsigned MEDDLY::expert_forest::UNIQUE_TABLE_DETAILED   = 0x0080;
const unsigned MEDDLY::expert_forest::HOLE_MANAGER_STATS      = 0x0100;
const unsigned MEDDLY::expert_forest::HOLE_MANAGER_DETAILED   = 0x0200;

//
// Display flags
//

const unsigned int MEDDLY::expert_forest::SHOW_DELETED      = 0x10;
const unsigned int MEDDLY::expert_forest::SHOW_UNREACHABLE  = 0x08;
const unsigned int MEDDLY::expert_forest::SHOW_DETAILS      = 0x04;
const unsigned int MEDDLY::expert_forest::SHOW_INDEX        = 0x02;
const unsigned int MEDDLY::expert_forest::SHOW_TERMINALS    = 0x01;


MEDDLY::expert_forest::expert_forest(domain *d, bool rel, range_type t,
  edge_labeling ev, const policies &p, int* level_reduction_rule)
: forest(d, rel, t, ev, p, level_reduction_rule), nodeHeaders(*this)
{
  nodeHeaders.setPessimistic(isPessimistic());

  // Initialize variable order
  var_order = d->makeDefaultVariableOrder();

  //
  // Initialize misc. protected data
  //
  terminalNodesStatus = MEDDLY::forest::ACTIVE;

  //
  // Initialize misc. private data
  //
  unique = new unique_table(this);
  implUT = new impl_unique_table(this);
  performing_gc = false;
  in_validate = 0;
  in_val_size = 0;
  delete_depth = 0;

  //
  // Initialize node characteristics to defaults
  //
  unhashed_bytes = 0;
  hashed_bytes = 0;
}


MEDDLY::expert_forest::~expert_forest()
{
#ifdef REPORT_ON_DESTROY
  printf("Destroying forest.  Stats:\n");
  reportMemoryUsage(stdout, "\t", 9);
#endif

  delete nodeMan;

  // unique table
  delete unique;
  delete implUT;

  // Misc. private data
  free(in_validate);
}

void MEDDLY::expert_forest::initializeForest()
{
  if (!deflt.nodestor) {
    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  }

  //
  // Initialize node storage
  //
  nodeMan = deflt.nodestor->createForForest(this, deflt.nodemm);

}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                  public mark & sweep  methods                  '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

void MEDDLY::expert_forest::markAllRoots()
{
    if (deflt.useReferenceCounts) return;

    stats.reachable_scans++;

#ifdef DEBUG_MARK_SWEEP
    printf("Determining which nodes are reachable in forest %u\n", FID());
#endif

    nodeHeaders.clearAllReachableBits();

    for (unsigned i=1; i<roots_next; i++) {
        dd_edge* r = roots[i].edge;
        if (r) markNode(r->getNode());
    }

    // TBD: this goes away when unpacked nodes use dd_edges :)
    unpacked_node::markBuildListChildren(this);
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                   public debugging   methods                   '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

void MEDDLY::expert_forest::dump(output &s, unsigned int flags) const
{
  for (long p=0; p<=nodeHeaders.lastUsedHandle(); p++) {
    if (showNode(s, p, flags | SHOW_INDEX)) {
      s.put('\n');
      s.flush();
    }
  }
}

void MEDDLY::expert_forest::dumpInternal(output &s) const
{
  s << "Internal forest storage\n";

  nodeHeaders.dumpInternal(s);

  nodeMan->dumpInternal(s, 0x03);

  // unique->show(s);
  s.flush();
}

void MEDDLY::expert_forest::dumpUniqueTable(output &s) const
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
      MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, P.d(z), sz);
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
      ?  in_validate[i] != getNodeInCount(i)
      :  in_validate[i] >  getNodeInCount(i);
    if (fail) {
      printf("Validation #%d failed\n", idnum);
      long l_i = i;
      long l_v = in_validate[i];
      long l_c = getNodeInCount(i);
      printf("For node %ld\n\tcount: %ld\n\tnode:  %ld\n", l_i, l_v, l_c);
      dump(stdout, SHOW_DETAILS);
      MEDDLY_DCASSERT(0);
      throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
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


void MEDDLY::expert_forest::validateCacheCounts() const
{
  if (!deflt.useReferenceCounts) return;

#ifdef DEVELOPMENT_CODE
#ifdef SHOW_VALIDATE_CACHECOUNTS
  printf("Validating cache counts for %ld handles\n", getLastNode());
#endif
  const node_handle N = getLastNode()+1;
  size_t* counts = new size_t[N];

#ifdef SHOW_VALIDATE_CACHECOUNTS
  printf("  Counting...\n");
  fflush(stdout);
#endif
  for (node_handle i=0; i<N; i++) counts[i] = 0;
  operation::countAllNodeEntries(this, counts);

#ifdef SHOW_VALIDATE_CACHECOUNTS
  printf("  Validating...\n");
#endif
  for (node_handle i=1; i<N; i++) {
    if (nodeHeaders.getNodeCacheCount(i) == counts[i]) continue;
    printf("\tCount mismatch node %ld\n", long(i));
    printf("\t  We counted %lu\n", counts[i]);
    printf("\t  node  says %lu\n", nodeHeaders.getNodeCacheCount(i));
  }
  node_handle maxi = 1;
  for (node_handle i=2; i<N; i++) {
    if (counts[i] > counts[maxi]) {
      maxi = i;
    }
  }
#ifdef SHOW_VALIDATE_CACHECOUNTS
  if (maxi < N) printf("  Largest count: %lu for node %ld at level %d\n",
    counts[maxi], maxi, nodeHeaders.getNodeLevel(maxi)
  );
#endif

  delete[] counts;
#endif
}


void MEDDLY::expert_forest::countNodesByLevel(long* active) const
{
  int L = getNumVariables();
  int l;
  if (isForRelations()) {
    l = -L;
  } else {
    l = 0;
  }

  for (; l<=L; l++) active[l] = 0;

  for (long p=1; p<=nodeHeaders.lastUsedHandle(); p++) {
    if (nodeHeaders.isDeleted(p)) continue;
    active[nodeHeaders.getNodeLevel(p)]++;
  }
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                      public handy methods                      '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

MEDDLY::node_marker*
MEDDLY::expert_forest::makeNodeMarker()
{
    node_marker* nm = new node_marker(false, nodeHeaders, nodeMan, this);
    nm->expand(nodeHeaders.lastUsedHandle()+1);
    nm->unmarkAll();
    return nm;
}

MEDDLY::node_handle*
MEDDLY::expert_forest
::markNodesInSubgraph(const node_handle* root, int N, bool sort) const
{
  MEDDLY_DCASSERT(root);

  const node_handle a_last = nodeHeaders.lastUsedHandle();
  // initialize lists
  bool* inList = new bool[a_last];
  for (int i=0; i<a_last; i++) inList[i] = false;
  inList--;

  int mlen = 0;
  int msize = 0;
  node_handle* marked = 0;

  // Initialize search
  for (int i=0; i<N; i++) {
    if (isTerminalNode(root[i])) continue;
    if (inList[root[i]]) continue;

    // add dn to list
    if (mlen+1 >= msize) {
      // expand.  Note we're leaving an extra slot
      // at the end, for the terminal 0.
      msize += 1024;
      node_handle* new_marked = (node_handle*)
        realloc(marked, msize*sizeof(node_handle));
      if (0==new_marked) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      marked = new_marked;
    }

    marked[mlen] = root[i];
    mlen++;
    inList[root[i]] = true;
  }

  unpacked_node *M = unpacked_node::New();

  // Breadth-first search
  for (int mexpl=0; mexpl<mlen; mexpl++) {
    // explore node marked[mexpl]
    unpackNode(M, marked[mexpl], SPARSE_ONLY);
    for (unsigned i=0; i<M->getNNZs(); i++) {
      if (isTerminalNode(M->d(i))) continue;
      MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, M->d(i)-1, a_last);
      if (inList[M->d(i)]) continue;
      // add dn to list
      if (mlen+1 >= msize) {
          // expand.  Note we're leaving an extra slot
          // at the end, for the terminal 0.
          msize += 1024;
          node_handle* new_marked = (node_handle*)
            realloc(marked, msize*sizeof(node_handle));
          if (0==new_marked) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
          marked = new_marked;
      }
      inList[M->d(i)] = true;
      marked[mlen] = M->d(i);
      mlen++;
    } // for i
  } // for mexpl

  unpacked_node::recycle(M);

  // sort
  if (sort && mlen>0) {
    mlen = 0;
    for (int i=1; i<=a_last; i++) {
      if (inList[i]) {
        marked[mlen] = i;
        mlen++;
      }
    }
  }

  // cleanup
  inList++;
  delete[] inList;

  if (0 == mlen) {
    if (marked) free(marked);
    return 0;
  }

  // add 0 to the list
  marked[mlen] = 0;
  return marked;
}

unsigned long MEDDLY::expert_forest::getNodeCount(node_handle p) const
{
    node_handle* list = markNodesInSubgraph(&p, 1, true);
    if (0==list) return 0;
    unsigned long i;
    for (i=0; list[i]; i++) { }
    free(list);
    return i;
}

unsigned long MEDDLY::expert_forest::getNodeCount(const node_handle* roots,
        int N) const
{
    node_handle* list = markNodesInSubgraph(roots, N, false);
    if (0==list) return 0;
    unsigned long i;
    for (i=0; list[i]; i++) { }
    free(list);
    return i;
}

unsigned long MEDDLY::expert_forest::getEdgeCount(node_handle p,
        bool countZeroes) const
{
    node_handle* list = markNodesInSubgraph(&p, 1, true);
    if (0==list) return 0;
    unsigned long ec=0;
    unpacked_node *M = unpacked_node::New();
    for (unsigned long i=0; list[i]; i++) {
        unpackNode(M, list[i], countZeroes ? FULL_ONLY : SPARSE_ONLY);
        ec += countZeroes ? M->getSize() : M->getNNZs();
    }
    unpacked_node::recycle(M);
    free(list);
    return ec;
}

bool MEDDLY::expert_forest
::showNode(output &s, node_handle p, unsigned int flags) const
{
  /*
    Deal with cases where nothing will be displayed.
  */
  bool isReachable =
    isTerminalNode(p) ||
    (
        deflt.useReferenceCounts ? (getNodeInCount(p)) : (hasReachableBit(p))
    );

  if (isTerminalNode(p)) {
    if (!(flags & SHOW_TERMINALS))  return false;
  } else
  if (isDeletedNode(p)) {
    if (!(flags & SHOW_DELETED))    return false;
  } else
  if (!isReachable) {
    if (!(flags & SHOW_UNREACHABLE))     return false;
  }

  /*
    Show the node index, if selected.
  */
  if (flags & SHOW_INDEX) {
    int nwidth = digits(nodeHeaders.lastUsedHandle());
    s.put(long(p), nwidth);
    s.put('\t');
  }

  /*
    Deal with special cases
  */
  if (isTerminalNode(p)) {
    s << "(terminal)";
    return true;
  }
  if (isDeletedNode(p)) {
    s << "DELETED";
    return true;
  }
  if (!isReachable) {
    s << "Unreachable ";
  }

  /*
    Ordinary node
  */
  if (flags & SHOW_DETAILS) {
    // node: was already written.
    nodeHeaders.showHeader(s, p);
  } else {
    s << "node: " << long(p);
  }

  s.put(' ');
  unpacked_node* un = newUnpacked(p, FULL_OR_SPARSE);
  un->show(s, flags & SHOW_DETAILS);
  unpacked_node::recycle(un);

  return true;
}


#ifdef ALLOW_DEPRECATED_0_17_3

void MEDDLY::expert_forest
::showNodeGraph(output &s, const node_handle* p, int n)
{
    node_marker* M = makeNodeMarker();
    MEDDLY_DCASSERT(M);

    for (int i=0; i<n; i++) {
        M->mark(p[i]);
    }

    M->showByLevels(s);

    delete M;
}


void MEDDLY::expert_forest
::writeNodeGraphPicture(const char* filename, const char *ext,
    const node_handle* p, const char* const* labels, int n)
{
    if (filename == NULL || ext == NULL || p == NULL) return;
    if (!isMultiTerminal()) {
        fprintf(stderr,
            "%s: Error. Only implemented for Multi-Terminal MDDs\n",
            __func__);
        return;
    }

    dot_maker DM(this, filename);

    for (unsigned i=0; i<n; i++) {
        dd_edge E(this);

        if (labels) E.setLabel(labels[i]);
        E.set_and_link(p[i]);

        DM.addRootEdge(E);
    }
    DM.doneGraph();
    DM.runDot(ext);
}
#endif

void MEDDLY::expert_forest
::reportStats(output &s, const char* pad, unsigned flags) const
{
  if (flags & BASIC_STATS) {
    bool human = flags & HUMAN_READABLE_MEMORY;
    s << pad << getCurrentNumNodes() << " current nodes\n";
    s << pad << getPeakNumNodes() << " peak nodes\n" << pad;
    s.put_mem(getCurrentMemoryUsed(), human);
    s << " current memory used\n" << pad;
    s.put_mem(getPeakMemoryUsed(), human);
    s << " peak memory used\n" << pad;
    s.put_mem(getCurrentMemoryAllocated(), human);
    s << " current memory allocated\n" << pad;
    s.put_mem(getPeakMemoryAllocated(), human);
    s << " peak memory allocated\n";
  }
  if (flags & EXTRA_STATS) {
    s << pad << stats.reachable_scans << " scans for reachable nodes\n";
    s << pad << stats.reclaimed_nodes << " reclaimed nodes\n";
    s << pad << stats.num_compactions << " compactions\n";
    s << pad << stats.garbage_collections << " garbage collections\n";
  }
  // forest specific
  reportForestStats(s, pad);
  // header storage
  nodeHeaders.reportStats(s, pad, flags);
  // node storage
  nodeMan->reportStats(s, pad, flags);
  // unique table
  unique->reportStats(s, pad, flags);
}


// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '   public virtual methods provided here,  no need to override   '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

void MEDDLY::expert_forest
::writeEdges(output &s, const dd_edge* E, int n) const
{
    node_handle* eRaw = new node_handle[n];
    for (int i=0; i<n; i++) {
        eRaw[i] = E[i].getNode();
    }
    node_handle* output2index = markNodesInSubgraph(eRaw, n, false);
    delete[] eRaw;

    // Make a dummy array for fringe case where
    // there are no non-zero nodes in the subgraph.
    if (0 == output2index) {
        output2index = (node_handle*) malloc(1*sizeof(node_handle));
        output2index[0] = 0;
    }

    // move a pointer to the end of the list, and
    // find the largest node index we're writing
    int maxnode = 0;
    int last;
    for (last = 0; output2index[last]; last++) {
        maxnode = MAX(maxnode, output2index[last]);
    };
    int num_nodes = last;
    last--;

    // arrange nodes to output, by levels
    for (int k=getNumVariables(); k; ) {
        int i = 0;
        while (i < last) {
            // move last to the left, until we have something not at level k
            for (; i < last; last--) {
                if (getNodeLevel(output2index[last]) != k) break;
            }
            // move i to the right, until we have something at level k
            for (; i < last; i++) {
                if (getNodeLevel(output2index[i]) == k) break;
            }
            if (i < last) {
                SWAP(output2index[i], output2index[last]);
            }
        }

        // next level
        k *= -1;
        if (k>0) k--;
    } // loop over levels

    // build the inverse mapping
    node_handle* index2output = new node_handle[maxnode+1];
    for (int i=0; i<maxnode; i++) index2output[i] = 0;
    for (int i=0; output2index[i]; i++) {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1, output2index[i], maxnode+1);
        index2output[output2index[i]] = i+1;
    }

#ifdef DEBUG_WRITE
    printf("Writing edges:\n");
    showNodeGraph(stdout, eRaw, n);
    printf("Got list of nodes:\n");
    for (int i=0; output2index[i]; i++) {
        if (i) printf(", ");
        printf("%ld", long(output2index[i]));
    }
    printf("\n");
    printf("Got inverse list:\n");
    for (int i=0; i<=maxnode; i++) {
        if (i) printf(", ");
        printf("%ld", long(index2output[i]));
    }
    printf("\n");
#endif

    // Write the nodes
    const char* block = codeChars();
    s << block << " " << num_nodes << "\n";

    unpacked_node* un = unpacked_node::New();
    for (int i=0; output2index[i]; i++) {
        s << getNodeLevel(output2index[i]) << " ";
        unpackNode(un, output2index[i], FULL_OR_SPARSE);
        un->write(s, index2output);
    }
    unpacked_node::recycle(un);

    // reverse the block
    for (int i=strlen(block); i; ) {
        i--;
        s.put(block[i]);
    }
    s.put('\n');

    // Write the actual edge pointers
    s << "ptrs " << n << "\n";
    for (int i=0; i<n; i++) {
        s.put('\t');
        E[i].write(s, index2output);
        s.put('\n');
    }
    s << "srtp\n";


    // Cleanup
    delete[] index2output;
    free(output2index);
}

void MEDDLY::expert_forest::readEdges(input &s, dd_edge* E, int n)
{
  try {
    const char* block = codeChars();

    s.stripWS();
    s.consumeKeyword(block);
    s.stripWS();
    int num_nodes = s.get_integer();
#ifdef DEBUG_READ_DD
    std::cerr << "Reading " << num_nodes << " nodes in " << block << " forest\n";
#endif
#ifdef DEBUG_READ
    printf("Reading %d nodes in forest %s\n", num_nodes, codeChars());
#endif

    // start a mapping
    node_handle* map = new node_handle[num_nodes+1];
    for (int i=0; i<=num_nodes; i++) map[i] = 0;

    for (int node_index=1; node_index<=num_nodes; node_index++) {
      // Read this node

      // read the level number
      s.stripWS();
      int k = s.get_integer();
      if (!isValidLevel(k)) {
        throw error(error::INVALID_LEVEL, __FILE__, __LINE__);
      }

#ifdef DEBUG_READ
      printf("Reading %d ", k);
#endif

      // read the node size (sparse/full)
      s.stripWS();
      int rawsize = s.get_integer();
      int n;
      unpacked_node* nb = (rawsize < 0)
        ? unpacked_node::newSparse(this, k, n=-rawsize)
        : unpacked_node::newFull(this, k, n=rawsize);

#ifdef DEBUG_READ
      printf("%d ", rawsize);
#endif

      //
      // read the node
      //
      nb->read(s, map);


      // ok, done reading; time to reduce it
      map[node_index] = createReducedNode(-1, nb);

#ifdef DEBUG_READ
      printf("\nNode index %d reduced to ", node_index);
      FILE_output myout(stdout);
      showNode(myout, map[node_index], SHOW_DETAILS | SHOW_INDEX);
      printf("\n");
#endif

  } // for node_index

    // reverse the block
    static char buffer[40];
    int blocklen = strlen(block);
    MEDDLY_DCASSERT(blocklen < 40);
    for (int i=0; i<blocklen; i++) {
      buffer[i] = block[blocklen-i-1];
    }
    buffer[blocklen] = 0;
#ifdef DEBUG_READ
    printf("Done reading, expecting %s keyword\n", buffer);
#endif

    // match the reversed block
    s.stripWS();
    s.consumeKeyword(buffer);
#ifdef DEBUG_READ
    printf("Got %s\n", buffer);
#endif
#ifdef DEBUG_READ_DD
    std::cerr << "Finished " << buffer << " forest\n";
#endif

    // Read the pointers
    s.stripWS();
    s.consumeKeyword("ptrs");
#ifdef DEBUG_READ
    printf("Got ptrs\n");
#endif

    s.stripWS();
    int num_ptrs = s.get_integer();
#ifdef DEBUG_READ
    printf("Reading %d pointers\n", num_ptrs);
#endif
#ifdef DEBUG_READ_DD
    std::cerr << "Reading " << num_ptrs << " pointers\n";
#endif
    MEDDLY_DCASSERT(num_ptrs <= n);
    if (num_ptrs > n) {
#ifdef DEBUG_READ
      printf("Error at %s:%d, E[] is of size %d, needs to be at least %d\n",
        __FILE__, __LINE__, n, num_ptrs);
      fflush(stdout);
#endif
      throw error(error::INVALID_ASSIGNMENT, __FILE__, __LINE__);
    }
    for (int i=0; i<num_ptrs; i++) {
        E[i].attach(this);
        E[i].read(s, map);
    }

    s.stripWS();
    s.consumeKeyword("srtp");

#ifdef DEBUG_READ_DD
    std::cerr << "Done reading pointers\n";
#endif

    // unlink map pointers
    for (int i=0; i<=num_nodes; i++) unlinkNode(map[i]);
    delete[] map;

#ifdef DEVELOPMENT_CODE
    validateIncounts(true);
#endif
  } // try
  catch (error& e) {
#ifdef DEBUG_READ
    printf("Read failed (error: %s)\n", e.getName());
    printf("Failed, next few characters of file:\n");
    for (int i=0; i<10; i++) {
      int c = s.get_char();
      if (EOF == c) {
        printf("EOF");
        break;
      }
      printf("%c (%d) ", c, c);
    }
    fputc('\n', stdout);
#endif
    throw e;
  }
}

/*

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
      FILE_output mystderr(stderr);
      showInfo(mystderr, 1);
      showComputeTable(mystderr, 6);
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
  nodeMan->collectGarbage(true);
}

*/

void MEDDLY::expert_forest::showInfo(output &s, int verb)
{
  // Show forest with appropriate level of detail
  if (1==verb)  dump(s, SHOW_DETAILS);
  else          dumpInternal(s);
  s << "DD stats:\n";
  reportStats(s, "    ", ~0);
  // reportMemoryUsage(s, "    ", verb);
  s << "Unique table stats:\n\t";
  s.put("Current size:", -24);
  s << long(unique->getSize()) << "\n\t";
  s.put("Current entries:", -24);
  s << long(unique->getNumEntries()) << "\n";
}


// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '    virtual methods to be overridden by some derived classes    '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

MEDDLY::enumerator::iterator* MEDDLY::expert_forest::makeFixedRowIter() const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

MEDDLY::enumerator::iterator*
MEDDLY::expert_forest::makeFixedColumnIter() const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

const char* MEDDLY::expert_forest::codeChars() const
{
  return "unknown dd";
}

void MEDDLY::expert_forest::normalize(unpacked_node &nb, int& ev) const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::expert_forest::normalize(unpacked_node &nb, long& ev) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::normalize(unpacked_node &nb, float& ev) const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::expert_forest::reportForestStats(output &s, const char* pad) const
{
  // default - do nothing
}


void MEDDLY::expert_forest::showHeaderInfo(output &s, const unpacked_node &) const
{
}

void MEDDLY::expert_forest::writeHeaderInfo(output &s, const unpacked_node &) const
{
}

void MEDDLY::expert_forest::readHeaderInfo(input &s, unpacked_node &) const
{
}


void MEDDLY::expert_forest::reorderVariables(const int* level2var)
{
  removeAllComputeTableEntries();

  // Create a temporary variable order
  // Support in-place update and avoid interfering other forests
  var_order = std::make_shared<variable_order>(*var_order);

  auto reordering = reordering_factory::create(getPolicies().reorder);
  reordering->reorderVariables(this, level2var);

  var_order = useDomain()->makeVariableOrder(*var_order);
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                                                                '
// '                        private  methods                        '
// '                                                                '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


void MEDDLY::expert_forest::deleteNode(node_handle p)
{
#ifdef TRACK_DELETIONS
  for (int i=0; i<delete_depth; i++) printf(" ");
  printf("Forest %u deleting node ", FID());
  FILE_output s(stdout);
  showNode(s, p, SHOW_INDEX | SHOW_DETAILS | SHOW_UNREACHABLE);
  printf("\n");
  fflush(stdout);
#endif
#ifdef VALIDATE_INCOUNTS_ON_DELETE
  delete_depth++;
#endif

  MEDDLY_DCASSERT(isValidNonterminalIndex(p));
  MEDDLY_DCASSERT(isActiveNode(p));
  if (deflt.useReferenceCounts) {
    MEDDLY_DCASSERT(getNodeInCount(p) == 0);
  } else {
    MEDDLY_DCASSERT(!hasReachableBit(p));
  }

  unsigned h = hashNode(p);
#ifdef DEVELOPMENT_CODE
  if (!isExtensible(p) || isExtensibleLevel(getNodeLevel(p))) {
    unpacked_node* key = newUnpacked(p, SPARSE_ONLY);
    key->computeHash();
    if (unique->find(*key, getVarByLevel(key->getLevel())) != p) {
      fprintf(stderr, "Error in deleteNode\nFind: %ld\np: %ld\n",
          static_cast<long>(unique->find(*key, getVarByLevel(key->getLevel()))), static_cast<long>(p));
      FILE_output myout(stdout);
      dumpInternal(myout);
      MEDDLY_DCASSERT(false);
    }
    node_handle x = unique->remove(h, p);
    MEDDLY_DCASSERT(p == x);
    unpacked_node::recycle(key);
  }
#else
  unique->remove(h, p);
#endif

  stats.decActive(1);

#ifdef TRACK_DELETIONS
  // start at one, because we have incremented the depth
  // for (int i=1; i<delete_depth; i++) printf(" ");
  // printf("%s: p = %d, unique->remove(p) = %d\n", __func__, p, x);
  // fflush(stdout);
#endif

  // unlink children and recycle node memory
  nodeMan->unlinkDownAndRecycle(getNodeAddress(p));
  setNodeAddress(p, 0);
  nodeHeaders.deactivate(p);

  // if (nodeMan.compactLevel) nodeMan.compact(false);

#ifdef VALIDATE_INCOUNTS_ON_DELETE
  delete_depth--;
  if (0==delete_depth) {
    validateIncounts(false);
  }
#endif

}



MEDDLY::node_handle MEDDLY::expert_forest
::createReducedHelper(int in, unpacked_node &nb)
{
    nb.computeHash();

#ifdef DEVELOPMENT_CODE
  validateDownPointers(nb);
#endif

  // check is the node is written in order,
  // if not rearrange it in ascending order of indices.
  if (nb.isSparse()) nb.sort();

  if (nb.isExtensible()) return createReducedExtensibleNodeHelper(in, nb);

  // get sparse, truncated full sizes and check
  // for redundant / identity reductions.
  int nnz;
  if (nb.isSparse()) {
    // Reductions for sparse nodes
    nnz = nb.getNNZs();
#ifdef DEVELOPMENT_CODE
    for (int z=0; z<nnz; z++) {
      MEDDLY_DCASSERT(nb.d(z)!=getTransparentNode());
    } // for z
#endif

    // Check for identity nodes
    if (1==nnz && in==nb.i(0)) {
      if (isIdentityEdge(nb, 0)) {
#ifdef DEBUG_CREATE_REDUCED
        printf("Identity node ");
        FILE_output s(stdout);
        showNode(s, nb.d(0), SHOW_DETAILS | SHOW_INDEX);
        printf("\n");
#endif
        return nb.d(0);
      }
    }

    // Check for redundant nodes
    if (nnz == getLevelSize(nb.getLevel()) && !isExtensibleLevel(nb.getLevel())) {
      if (isRedundant(nb)) {
        // unlink downward pointers, except the one we're returning.
        for (int i = 1; i<nnz; i++)  unlinkNode(nb.d(i));
#ifdef DEBUG_CREATE_REDUCED
        printf("Redundant node ");
        FILE_output s(stdout);
        showNode(s, nb.d(0), SHOW_DETAILS | SHOW_INDEX);
        printf("\n");
#endif
        return nb.d(0);
      }
    }

  } else {
    // Reductions for full nodes
    MEDDLY_DCASSERT(nb.getSize() <= getLevelSize(nb.getLevel()));
    nnz = 0;
    for (unsigned i=0; i<nb.getSize(); i++) {
      if (nb.d(i)!=getTransparentNode()) nnz++;
    } // for i

    // Check for identity nodes
    if (1==nnz) {
      if (in < nb.getSize() && isIdentityEdge(nb, in)) {
#ifdef DEBUG_CREATE_REDUCED
        printf("Identity node ");
        FILE_output s(stdout);
        showNode(s, nb.d(0), SHOW_DETAILS | SHOW_INDEX);
        printf("\n");
#endif
        return nb.d(in);
      }
    }

    // Check for redundant nodes
    if (nnz == getLevelSize(nb.getLevel()) && !isExtensibleLevel(nb.getLevel())) {
      if (isRedundant(nb)) {
        // unlink downward pointers, except the one we're returning.
        for (unsigned i = 1; i<nb.getSize(); i++)  unlinkNode(nb.d(i));
#ifdef DEBUG_CREATE_REDUCED
        printf("Redundant node ");
        FILE_output s(stdout);
        showNode(s, nb.d(0), SHOW_DETAILS | SHOW_INDEX);
        printf("\n");
#endif
        return nb.d(0);
      }
    }
  }

  // Is this a transparent node?
  if (0==nnz) {
    // no need to unlink
    return getTransparentNode();
  }

  // check for duplicates in unique table
  node_handle q = unique->find(nb, getVarByLevel(nb.getLevel()));
  if (q) {
    // unlink all downward pointers
    int rawsize = nb.isSparse() ? nb.getNNZs() : nb.getSize();
    for (int i = 0; i<rawsize; i++)  unlinkNode(nb.d(i));
    return linkNode(q);
  }

  //
  // Not eliminated by reduction rule.
  // Not a duplicate.
  //
  // We need to create a new node for this.

  // NOW is the best time to run the garbage collector, if necessary.
#ifndef GC_OFF
  // if (isTimeToGc()) garbageCollect();
#endif

  // Grab a new node
  node_handle p = nodeHeaders.getFreeNodeHandle();
  nodeHeaders.setNodeLevel(p, nb.getLevel());
  if (deflt.useReferenceCounts) {
    MEDDLY_DCASSERT(0 == nodeHeaders.getIncomingCount(p));
    MEDDLY_DCASSERT(0 == nodeHeaders.getNodeCacheCount(p));
  } else {
    nodeHeaders.setReachableBit(p);
    nodeHeaders.setInCacheBit(p);
  }

  stats.incActive(1);
  if (theLogger && theLogger->recordingNodeCounts()) {
    theLogger->addToActiveNodeCount(this, nb.getLevel(), 1);
  }

  // All of the work is in nodeMan now :^)
  nodeHeaders.setNodeAddress(p, nodeMan->makeNode(p, nb, getNodeStorage()));
  linkNode(p);

  // add to UT
  unique->add(nb.hash(), p);

#ifdef DEVELOPMENT_CODE
  unpacked_node* key = newUnpacked(p, SPARSE_ONLY);
  key->computeHash();
  MEDDLY_DCASSERT(key->hash() == nb.hash());
  node_handle f = unique->find(*key, getVarByLevel(key->getLevel()));
  MEDDLY_DCASSERT(f == p);
  unpacked_node::recycle(key);
#endif
#ifdef DEBUG_CREATE_REDUCED
  printf("Created node ");
  FILE_output s(stdout);
  showNode(s, p, SHOW_DETAILS | SHOW_INDEX);
  printf("\n");
#endif

  return p;
}

unsigned MEDDLY::expert_forest
::getImplTableCount() const
{
    return implUT->getNumEntries();
}

MEDDLY::relation_node* MEDDLY::expert_forest
::buildImplNode(node_handle rnh)
{
  return implUT->getNode(rnh);
}

MEDDLY::node_handle MEDDLY::expert_forest
::getImplTerminalNode() const
{
  return implUT->getLastHandle();
}


MEDDLY::node_handle MEDDLY::expert_forest
::createImplicitNode(MEDDLY::relation_node &nb)
{
  // check for duplicates in unique table

  node_handle q = implUT->isDuplicate(&nb);
  if (q) return q;

  // Not eliminated by reduction rule.
  // Not a duplicate.
  //
  // We need to create a new node for this.

  // NOW is the best time to run the garbage collector, if necessary.
#ifndef GC_OFF
  // if (isTimeToGc()) garbageCollect();
#endif

  // Grab a new node
  node_handle p = nodeHeaders.getFreeNodeHandle();
  nodeHeaders.setNodeImplicitFlag(p, true);
  nodeHeaders.setNodeLevel(p, nb.getLevel());
  MEDDLY_DCASSERT(0 == nodeHeaders.getNodeCacheCount(p));
  MEDDLY_DCASSERT(0 == nodeHeaders.getIncomingCount(p));

  stats.incActive(1);
  if (theLogger && theLogger->recordingNodeCounts()) {
    theLogger->addToActiveNodeCount(this, nb.getLevel(), 1);
  }

  #if 0
  // All of the work is in satimpl_opname::implicit_relation now :^)
  nodeHeaders.setNodeAddress(p, nb.getID());
  linkNode(p);

  // add to UT
  unique->add(nb.getSignature(), p);
  #endif

  // All of the work is in implUT now :^)
  // If it is not duplicate, universal handle will become handle to node in implUT.
  nb.setID(p);
  nodeHeaders.setNodeAddress(p, implUT->add(p, &nb));
  linkNode(p);


#ifdef DEBUG_CREATE_REDUCED
  printf("Created node ");
  FILE_output s(stdout);
  showNode(s, p, SHOW_DETAILS | SHOW_INDEX);
  printf("\n");
#endif

  return p;
}



MEDDLY::node_handle MEDDLY::expert_forest
::createReducedExtensibleNodeHelper(int in, unpacked_node &nb)
{
#ifdef DEVELOPMENT_CODE
  validateDownPointers(nb);
#endif
  MEDDLY_DCASSERT(nb.isExtensible());
  MEDDLY_DCASSERT(nb.isTrim());

  // NOTE: Identity reduction not possible for nodes marked as extensible.
  //       Fully-Identity reduction is still possible when
  //       prime-level nodes are non-extensible, and get Identity reduced.

  // get sparse, truncated full sizes and check
  // for redundant / identity reductions.
  int nnz = 0;
  const int rawsize = nb.isSparse()? nb.getNNZs(): nb.getSize();
  for (int i=0; i<rawsize; i++) {
    if (nb.d(i)!=getTransparentNode()) nnz++;
  } // for i

  // Is this a transparent node?
  if (0==nnz) {
    // no need to unlink
    return getTransparentNode();
  }

  // Check for redundant nodes
  if (isRedundant(nb)) {
    MEDDLY_DCASSERT(nnz == 1 && nb.ext_i() == 0);
#ifdef DEBUG_CREATE_REDUCED
    printf("Redundant node ");
    FILE_output s(stdout);
    showNode(s, nb.ext_d(), SHOW_DETAILS | SHOW_INDEX);
    printf("\n");
#endif
    return nb.ext_d();
  }

  // check for duplicates in unique table
  node_handle q = unique->find(nb, getVarByLevel(nb.getLevel()));
  if (q) {
    // unlink all downward pointers
    for (int i = 0; i<rawsize; i++)  unlinkNode(nb.d(i));
    return linkNode(q);
  }

  //
  // Not eliminated by reduction rule.
  // Not a duplicate.
  //
  // We need to create a new node for this.

  // NOW is the best time to run the garbage collector, if necessary.
#ifndef GC_OFF
  // if (isTimeToGc()) garbageCollect();
#endif

  // Expand level size
  const int nb_ext_i = nb.ext_i();
  if (nb_ext_i >= getLevelSize(nb.getLevel())) {
    getDomain()->enlargeVariableBound(nb.getLevel(), false, -(nb_ext_i+1));
  }

  // Grab a new node
  node_handle p = nodeHeaders.getFreeNodeHandle();
  nodeHeaders.setNodeLevel(p, nb.getLevel());
  MEDDLY_DCASSERT(0 == nodeHeaders.getNodeCacheCount(p));
  MEDDLY_DCASSERT(0 == nodeHeaders.getIncomingCount(p));

  stats.incActive(1);
  if (theLogger && theLogger->recordingNodeCounts()) {
    theLogger->addToActiveNodeCount(this, nb.getLevel(), 1);
  }

  // All of the work is in nodeMan now :^)
  nodeHeaders.setNodeAddress(p, nodeMan->makeNode(p, nb, getNodeStorage()));
  // TODO: need to link?
  linkNode(p);

  // add to UT
  unique->add(nb.hash(), p);

#ifdef DEVELOPMENT_CODE
  unpacked_node* key = newUnpacked(p, SPARSE_ONLY);
  key->computeHash();
  MEDDLY_DCASSERT(key->hash() == nb.hash());
  node_handle f = unique->find(*key, getVarByLevel(key->getLevel()));
  MEDDLY_DCASSERT(f == p);
  unpacked_node::recycle(key);
#endif
#ifdef DEBUG_CREATE_REDUCED
  printf("Created node ");
  FILE_output s(stdout);
  showNode(s, p, SHOW_DETAILS | SHOW_INDEX);
  printf("\n");
#endif

  return p;
}

void MEDDLY::expert_forest::swapNodes(node_handle p, node_handle q)
{
  unique->remove(hashNode(p), p);
  unique->remove(hashNode(q), q);

  nodeHeaders.swapNodes(p, q, false);

  unique->add(hashNode(p), p);
  unique->add(hashNode(q), q);
}

MEDDLY::node_handle MEDDLY::expert_forest::modifyReducedNodeInPlace(unpacked_node* un, node_handle p)
{
  unique->remove(hashNode(p), p);
  nodeMan->unlinkDownAndRecycle(nodeHeaders.getNodeAddress(p));

  un->computeHash();

  nodeHeaders.setNodeLevel(p, un->getLevel());
  node_address addr = nodeMan->makeNode(p, *un, getNodeStorage());
  nodeHeaders.setNodeAddress(p, addr);
  // incoming count, cache count remains unchanged

  unique->add(un->hash(), p);

#ifdef DEVELOPMENT_CODE
  unpacked_node* key = newUnpacked(p, SPARSE_ONLY);
  key->computeHash();
  MEDDLY_DCASSERT(key->hash() == un->hash());
  node_handle f = unique->find(*key, getVarByLevel(key->getLevel()));
  MEDDLY_DCASSERT(f == p);
  unpacked_node::recycle(key);
#endif
#ifdef DEBUG_CREATE_REDUCED
  printf("Created node ");
  FILE_output s(stdout);
  showNode(s, p, SHOW_DETAILS | SHOW_INDEX);
  printf("\n");
#endif

  unpacked_node::recycle(un);
  return p;
}

void MEDDLY::expert_forest::validateDownPointers(const unpacked_node &nb) const
{
  switch (getReductionRule()) {
    case reduction_rule::IDENTITY_REDUCED:
    case reduction_rule::FULLY_REDUCED:
      if (nb.isSparse()) {
        for (unsigned z=0; z<nb.getNNZs(); z++) {
          if (isTerminalNode(nb.d(z))) continue;
          MEDDLY_DCASSERT(!isDeletedNode(nb.d(z)));
          if (isLevelAbove(nb.getLevel(), getNodeLevel(nb.d(z)))) continue;
          FILE_output s(stdout);
          s << "Down pointer violation in created node at level " << nb.getLevel() << ":\n";
          nb.show(s, true);
          s << "\nPointer " << nb.d(z) << " at level " << getNodeLevel(nb.d(z)) << ":\n";
          showNode(s, nb.d(z), SHOW_DETAILS);
          MEDDLY_DCASSERT(0);
        }
      } else {
        for (unsigned i=0; i<nb.getSize(); i++) {
          if (isTerminalNode(nb.d(i))) continue;
          MEDDLY_DCASSERT(!isDeletedNode(nb.d(i)));
          if (isLevelAbove(nb.getLevel(), getNodeLevel(nb.d(i)))) continue;
          FILE_output s(stdout);
          s << "Down pointer violation in created node at level " << nb.getLevel() << ":\n";
          nb.show(s, true);
          s << "\nPointer " << nb.d(i) << " at level " << getNodeLevel(nb.d(i)) << ":\n";
          showNode(s, nb.d(i), SHOW_DETAILS);
          MEDDLY_DCASSERT(0);
        }
      }
      break;

    case reduction_rule::QUASI_REDUCED:
#ifdef DEVELOPMENT_CODE
      int nextLevel;
      if (isForRelations())
        nextLevel = (nb.getLevel()<0) ? -(nb.getLevel()+1) : -nb.getLevel();
      else
        nextLevel = nb.getLevel()-1;
#endif
      if (nb.isSparse()) {
        for (unsigned z=0; z<nb.getNNZs(); z++) {
          MEDDLY_DCASSERT(getNodeLevel(nb.d(z)) == nextLevel);
        }
      } else {
        node_handle tv=getTransparentNode();
        for (unsigned i=0; i<nb.getSize(); i++) {
          if (nb.d(i)==tv) continue;
          MEDDLY_DCASSERT(getNodeLevel(nb.d(i)) == nextLevel);
        }
      }
      break;

    default:
      throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
  }

}

void MEDDLY::expert_forest::recycle(unpacked_node* n)
{
    unpacked_node::recycle(n);
}

//
// Stuff that used to be inlined but now can't
//

void
MEDDLY::expert_forest::unpackNode(MEDDLY::unpacked_node* un,
    MEDDLY::node_handle node, node_storage_flags st2) const
{
    MEDDLY_DCASSERT(un);
    const int level = getNodeLevel(node);
    MEDDLY_DCASSERT(0 != level);
    un->bind_to_forest(this, level, unsigned(getLevelSize(level)), true);
    MEDDLY_DCASSERT(getNodeAddress(node));
    nodeMan->fillUnpacked(*un, getNodeAddress(node), st2);
}

//----------------------------------------------------------------------
// front end - create and destroy objects
//----------------------------------------------------------------------


void MEDDLY::destroyForest(MEDDLY::forest* &f)
{
  if (0==f) return;
  if (!initializer_list::libraryIsRunning()) {
      throw error(error::UNINITIALIZED, __FILE__, __LINE__);
  }
  f->markForDeletion();
  operation::purgeAllMarked();
  delete f;
  f = 0;
}



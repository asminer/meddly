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
#include "logger.h"
#include "compute_table.h"
#include "ct_entry_type.h"  // for invalidateAllWithForest

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

std::vector <MEDDLY::forest*> MEDDLY::forest::all_forests;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Static "constructors" and "destructors"
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


void MEDDLY::forest::destroy(forest* &f)
{
    if (!f) return;
    if (!initializer_list::libraryIsRunning()) {
        throw error(error::UNINITIALIZED, __FILE__, __LINE__);
    }
    // f->markForDeletion();
    // operation::purgeAllMarked();
    delete f;
    f = nullptr;
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Node unpacking methods
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::forest::unpackNode(MEDDLY::unpacked_node* un,
    node_handle node, node_storage_flags st2) const
{
    MEDDLY_DCASSERT(un);
    MEDDLY_DCASSERT(un->isAttachedTo(this));
    un->setLevel( getNodeLevel(node) );
    MEDDLY_DCASSERT(un->getLevel() != 0);
    un->resize( getLevelSize( un->getLevel() ) );
    MEDDLY_DCASSERT(getNodeAddress(node));
    nodeMan->fillUnpacked(*un, getNodeAddress(node), st2);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Methods to add/remove nodes
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MEDDLY::relation_node* MEDDLY::forest
::buildImplicitNode(node_handle rnh)
{
    return implUT->getNode(rnh);
}

void MEDDLY::forest::createReducedNode(unpacked_node *un, edge_value &ev,
                node_handle &node, int in)
{
    MEDDLY_DCASSERT(un);
#ifdef DEBUG_CREATE_REDUCED
    ostream_output out(std::cout);
    out << "Reducing unpacked node ";
    un->show(out, true);
    out << " in forest " << FID() << "\n";
#endif
    //
    //
    // Normalize the node and count nonzeroes.
    //
    //
    unsigned minplusone = 0;
    unsigned nnz = 0;
    switch (edgeLabel) {
        case edge_labeling::MULTI_TERMINAL:
                ev.set();
                for (unsigned i=0; i<un->getSize(); i++) {
                    if (un->down(i)) ++nnz;
                }
                break;

        case edge_labeling::EVPLUS:
        case edge_labeling::INDEX_SET:
                MEDDLY_DCASSERT(isRangeType(range_type::INTEGER));
                ev.set(0L);
                for (unsigned i=0; i<un->getSize(); i++) {
                    if (0 == un->down(i)) continue;
                    ++nnz;
                    if (minplusone) {
                        if (un->edgeval(i).getLong() < ev.getLong()) {
                            ev = un->edgeval(i);
                        }
                    } else {
                        minplusone = i+1;
                        ev = un->edgeval(i);
                    }
                }
                if (ev.getLong()) {
                    //
                    // non-zero adjustment
                    //
                    for (unsigned i=0; i<un->getSize(); i++) {
                        if (0 == un->down(i)) continue;
                        un->subtractFromEdge(i, ev.getLong());
                    }
                }
                break;

        case edge_labeling::EVTIMES:
                MEDDLY_DCASSERT(isRangeType(range_type::REAL));
                ev.set(0.0f);
                minplusone = 0;
                for (unsigned i=0; i<un->getSize(); i++) {
                    if (0 == un->down(i)) continue;
                    if (un->edgeval(i).getFloat()) {
                        ++nnz;
                        if (!minplusone) {
                            ev = un->edgeval(i);
                            minplusone = i+1;
                        }
                    }
                }
                if (ev.getFloat()) {
                    for (unsigned i=0; i<un->getSize(); i++) {
                        un->divideEdge(i, ev.getFloat());
                    }
                }
                break;

        default:
                MEDDLY_DCASSERT(false);
    }

    //
    // Is this a transparent node?
    //
    if (0==nnz) {
#ifdef DEBUG_CREATE_REDUCED
        out << "    ===> Transparent node\n";
#endif
        // nothing to unlink
        node = getTransparentNode();
        unpacked_node::Recycle(un);
        return;
    }

    //
    // check is the node is written in order,
    // if not rearrange it in ascending order of indices.
    //
    if (un->isSparse()) {
        un->sort();
    }
#ifdef DEVELOPMENT_CODE
    validateDownPointers(*un);
#endif

#ifdef ALLOW_EXTENSIBLE
//    if (un->isExtensible()) return createReducedExtensibleNodeHelper(in, *un);
#endif

    //
    // Eliminate identity patterns
    //
    if (1==nnz && isIdentityReduced() && un->getLevel() < 0) {

        //
        // Check identity pattern
        //
        if (un->isSparse()) {
            if (in == long(un->index(0))) {
#ifdef DEBUG_CREATE_REDUCED
                out << "    ===> Identity node\n";
#endif
                node = un->down(0);
                unpacked_node::Recycle(un);
                return;
            }
        } else {
            if (in>=0 && in<long(un->getSize()) && un->down(in))
            {
#ifdef DEBUG_CREATE_REDUCED
                out << "    ===> Identity node\n";
#endif
                node = un->down(in);
                unpacked_node::Recycle(un);
                return;
            }
        }

    } // check identity patterns


    //
    // Eliminate redundant patterns
    //
    if ( isFullyReduced() || (isIdentityReduced() && un->getLevel() > 0))
    {
        //
        // Redundant node check
        //
        bool redundant = (nnz == getLevelSize(un->getLevel()));
        if (redundant) {
            node = un->down(0);
            for (unsigned i=1; i<nnz; i++) {
                if (un->down(i) == node) continue;
                redundant = false;
                break;
            }
        }
        if (redundant && un->hasEdges()) {
            for (unsigned i=1; i<nnz; i++) {
                if (un->edgeval(i) == un->edgeval(0)) continue;
                redundant = false;
                break;
            }
        }
        if (redundant) {
#ifdef DEBUG_CREATE_REDUCED
            out << "    ===> Redundant node\n";
#endif
            unlinkAllDown(*un, 1);  // unlink all children but one
            unpacked_node::Recycle(un);
            return;
        }
    } // isFullyReduced()

    //
    // Check for a duplicate node in the unique table
    //
    un->computeHash();
    node = unique->find(*un, getVarByLevel(un->getLevel()));
    if (node) {
#ifdef DEBUG_CREATE_REDUCED
        out << "    ===> Duplicate of node " << node << "\n";
#endif
        // unlink all downward pointers
        unlinkAllDown(*un);
        unpacked_node::Recycle(un);
        linkNode(node);
        return;
    }

    //
    // This node is new. Grab a node handle we can use.
    //
    node = nodeHeaders.getFreeNodeHandle();
    nodeHeaders.setNodeLevel(node, un->getLevel());
    if (reachable) {
        reachable->setMarked(node);
        nodeHeaders.setInCacheBit(node);
    } else {
        MEDDLY_DCASSERT(0 == nodeHeaders.getIncomingCount(node));
        MEDDLY_DCASSERT(0 == nodeHeaders.getNodeCacheCount(node));
        linkNode(node);
    }

    stats.incActive(1);
    if (theLogger && theLogger->recordingNodeCounts()) {
        theLogger->addToActiveNodeCount(this, un->getLevel(), 1);
    }

    //
    // Copy the temp node into long-term node storage
    //
    nodeHeaders.setNodeAddress(node,
        nodeMan->makeNode(node, *un, getPolicies().storage_flags)
    );

    //
    // add to the unique table
    //
    unique->add(un->hash(), node);

#ifdef DEBUG_CREATE_REDUCED
    out << "    ===> New node " << node << "\n\t";
    showNode(out, node, SHOW_DETAILS | SHOW_INDEX);
    out << '\n';
#endif

    //
    // Sanity check: can we find the node we just created
    //
#ifdef DEVELOPMENT_CODE
    unpacked_node* key = newUnpacked(node, SPARSE_ONLY);
    key->computeHash();
    MEDDLY_DCASSERT(key->hash() == un->hash());
    node_handle f = unique->find(*key, getVarByLevel(key->getLevel()));
    MEDDLY_DCASSERT(f == node);
    unpacked_node::Recycle(key);
#endif

    //
    // Cleanup
    //
    unpacked_node::Recycle(un);
}

void MEDDLY::forest::deleteNode(node_handle p)
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

    CHECK_RANGE(__FILE__, __LINE__, 1, p, 1+nodeHeaders.lastUsedHandle());
    MEDDLY_DCASSERT(isActiveNode(p));
    if (reachable) {
        MEDDLY_DCASSERT(!reachable->isMarked(p));
    } else {
        MEDDLY_DCASSERT(getNodeInCount(p) == 0);
    }

    unsigned h = hashNode(p);
#ifdef DEVELOPMENT_CODE
    if (!isExtensible(p) || isExtensibleLevel(getNodeLevel(p))) {
        unpacked_node* key = newUnpacked(p, SPARSE_ONLY);
        key->computeHash();
        if (unique->find(*key, getVarByLevel(key->getLevel())) != p) {
            fprintf(stderr, "Error in deleteNode\nFind: %ld\np: %ld\n",
            static_cast<long>(unique->find(*key,
                    getVarByLevel(key->getLevel()))), static_cast<long>(p));
            FILE_output myout(stdout);
            dumpInternal(myout);
            MEDDLY_DCASSERT(false);
        }
        node_handle x = unique->remove(h, p);
        MEDDLY_DCASSERT(p == x);
        unpacked_node::Recycle(key);
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


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Node packing helper methods
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#ifdef ALLOW_EXTENSIBLE
MEDDLY::node_handle MEDDLY::forest
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
  for (unsigned i=0; i<nb.getSize(); i++) {
    if (nb.down(i)!=getTransparentNode()) nnz++;
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
    unlinkAllDown(nb);
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
  nodeHeaders.setNodeAddress(p, nodeMan->makeNode(p, nb, getPolicies().storage_flags));
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
  unpacked_node::Recycle(key);
#endif
#ifdef DEBUG_CREATE_REDUCED
  printf("Created node ");
  FILE_output s(stdout);
  showNode(s, p, SHOW_DETAILS | SHOW_INDEX);
  printf("\n");
#endif

  return p;
}
#endif

MEDDLY::node_handle MEDDLY::forest
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



MEDDLY::node_handle
MEDDLY::forest::_makeRedundantsTo(node_handle p, int K, int L)
{
    MEDDLY_DCASSERT(L);
    unpacked_node* U;
    if (isMultiTerminal()) {
        //
        // Multi-terminal
        //
        while (K != L) {
            if (isForRelations()) {
                K = upLevel(K);
            } else {
                K++;
            }
            U = unpacked_node::newRedundant(this, K, p, FULL_ONLY);
            linkAllDown(*U, 1);
            edge_value ev;
            createReducedNode(U, ev, p);
            MEDDLY_DCASSERT(ev.isVoid());
        }
    } else {
        //
        // Edge-valued
        //

        MEDDLY_DCASSERT(false);
    }


    return p;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Moving nodes around
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::forest::swapNodes(node_handle p, node_handle q)
{
    unique->remove(hashNode(p), p);
    unique->remove(hashNode(q), q);

    nodeHeaders.swapNodes(p, q, false);

    unique->add(hashNode(p), p);
    unique->add(hashNode(q), q);
}

MEDDLY::node_handle MEDDLY::forest
    ::modifyReducedNodeInPlace(unpacked_node* un, node_handle p)
{
    unique->remove(hashNode(p), p);
    nodeMan->unlinkDownAndRecycle(nodeHeaders.getNodeAddress(p));

    un->computeHash();

    nodeHeaders.setNodeLevel(p, un->getLevel());
    node_address addr = nodeMan->makeNode(p, *un, getPolicies().storage_flags);
    nodeHeaders.setNodeAddress(p, addr);
    // incoming count, cache count remains unchanged

    unique->add(un->hash(), p);

#ifdef DEVELOPMENT_CODE
    unpacked_node* key = newUnpacked(p, SPARSE_ONLY);
    key->computeHash();
    MEDDLY_DCASSERT(key->hash() == un->hash());
    node_handle f = unique->find(*key, getVarByLevel(key->getLevel()));
    MEDDLY_DCASSERT(f == p);
    unpacked_node::Recycle(key);
#endif
#ifdef DEBUG_CREATE_REDUCED
    printf("Created node ");
    FILE_output s(stdout);
    showNode(s, p, SHOW_DETAILS | SHOW_INDEX);
    printf("\n");
#endif

    unpacked_node::Recycle(un);
    return p;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Node manager initialization
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::forest::initializeStorage()
{
  //
  // Initialize node storage
  //

  if (!deflt.nodestor) {
    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  }
  nodeMan = deflt.nodestor->createForForest(this, deflt.nodemm, mstats);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// I/O methods
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::forest::showHeaderInfo(output &s, const unpacked_node &) const
{
}

void MEDDLY::forest::writeHeaderInfo(output &s, const unpacked_node &) const
{
}

void MEDDLY::forest::readHeaderInfo(input &s, unpacked_node &) const
{
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Root edge registry methods
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::forest::registerEdge(dd_edge& e)
{
    e.next = roots;
    if (roots) {
        roots->prev = &e;
    }
    e.prev = nullptr;
    roots = &e;
    e.parentFID = fid;
}


void MEDDLY::forest::unregisterEdge(dd_edge& e)
{
    // remove a root edge.
    MEDDLY_DCASSERT(e.parentFID == fid);

    if (e.prev) {
        e.prev->next = e.next;
    } else {
        MEDDLY_DCASSERT(&e == roots);
        roots = e.next;
    }
    if (e.next) {
        e.next->prev = e.prev;
    }
    e.parentFID = 0;
}

unsigned MEDDLY::forest::countRegisteredEdges() const
{
    unsigned count = 0;
    for (const dd_edge* r = roots; r; r=r->next) {
        ++count;
    }
    return count;
}

void MEDDLY::forest::markAllRoots()
{
    if (!reachable) return;

    stats.reachable_scans++;

#ifdef DEBUG_MARK_SWEEP
    printf("Determining which nodes are reachable in forest %u\n", FID());
#endif

    reachable->unmarkAll();

    for (const dd_edge* r = roots; r; r=r->next) {
        reachable->mark(r->getNode());
    }

    unpacked_node::MarkWritable(*reachable);
}


void MEDDLY::forest::unregisterDDEdges()
{
    // Unregister ALL root edges
    // (e.g., because we're destroying the forest)

    while (roots) {
        dd_edge* n = roots->next;

        // Clear out the edge.
        roots->next = nullptr;
        roots->prev = nullptr;
        roots->parentFID = 0;
        roots->node = 0;

        roots = n;
    }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Operation registry methods
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::forest::removeAllComputeTableEntries()
{
    if (is_marked_for_deletion) return;
    //
    // Remove all compute table entries.
    //
    is_marked_for_deletion = true;
    if (compute_table::Monolithic()) {
        compute_table::Monolithic()->removeStales();
    }
    is_marked_for_deletion = false;
    ct_entry_type::removeAllCTEntriesWithForest(this);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Forest registry methods
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::forest::initStatics()
{
    all_forests.clear();
    all_forests.push_back(nullptr);
    // reserve index 0 for 'no forest'
}

void MEDDLY::forest::freeStatics()
{
    all_forests.clear();
}

void MEDDLY::forest::registerForest(forest* f)
{
#ifdef DEBUG_CLEANUP
    std::cout << "Registering forest " << f << ", #" << all_forests.size() << "\n";
#endif
    // Add f to the forest registry
    f->fid = all_forests.size();
    all_forests.push_back(f);

    // Register in the domain
    f->d->registerForest(f);

    // Initialize unpacked node entries for f
    unpacked_node::initForest(f);
}

void MEDDLY::forest::unregisterForest(forest* f)
{
#ifdef DEBUG_CLEANUP
    std::cout << "Unregistering forest " << f << ", #" << f->fid << "\n";
#endif
    // Remove from forest slot
    if (f->fid < all_forests.size()) {
#ifdef DEVELOPMENT_CODE
        all_forests.at(f->fid) = nullptr;
#else
        all_forests[f->fid] = nullptr;
#endif
    }

    // Unregister in the domain
    f->d->unregisterForest(f);

    // Clear out unpacked node entries for f
    unpacked_node::doneForest(f);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Methods for debugging / logging
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool MEDDLY::forest
::showNode(output &s, node_handle p, display_flags flags) const
{
  /*
    Deal with cases where nothing will be displayed.
  */
  bool isReachable =
    isTerminalNode(p) ||
    (
        reachable
            ? reachable->isMarked(p)
            : getNodeInCount(p)
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
  unpacked_node::Recycle(un);

  return true;
}



void MEDDLY::forest
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



void MEDDLY::forest::dump(output &s, display_flags flags) const
{
  for (long p=0; p<=nodeHeaders.lastUsedHandle(); p++) {
    if (showNode(s, p, flags | SHOW_INDEX)) {
      s.put('\n');
      s.flush();
    }
  }
}

void MEDDLY::forest::dumpInternal(output &s) const
{
  s << "Internal forest storage\n";

  nodeHeaders.dumpInternal(s);

  nodeMan->dumpInternal(s, 0x03);

  // unique->show(s);
  s.flush();
}

void MEDDLY::forest::dumpUniqueTable(output &s) const
{
  unique->show(s);
}

void MEDDLY::forest::validateIncounts(bool exact)
{
#ifndef VALIDATE_INCOUNTS
    return;
#endif

    static int idnum = 0;
    idnum++;

    // Inspect every active node's down pointers to determine
    // the incoming count for every active node.

    node_handle sz = getLastNode() + 1;
    std::vector <unsigned> in_validate(sz);

    unpacked_node P(this);
    for (node_handle i = 1; i < sz; ++i) {
        if (!isActiveNode(i)) continue;
        unpackNode(&P, i, SPARSE_ONLY);

        // add to reference counts
        for (unsigned z=0; z<P.getSize(); z++) {
            if (isTerminalNode(P.down(z))) continue;
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, P.down(z), sz);
            in_validate[P.down(z)]++;
        }
    } // for i

    // Add counts for registered dd_edges
    for (const dd_edge* r = roots; r; r=r->next) {
        if (isTerminalNode(r->getNode())) continue;
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, r->getNode(), sz);
        in_validate[r->getNode()]++;
    }

    //
    // TBD: add counts for down pointers in unpacked_nodes
    // in the build list for this forest.
    //

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
            FILE_output fout(stdout);
            dump(fout, SHOW_DETAILS);
            fout.flush();
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
}


void MEDDLY::forest::validateCacheCounts() const
{
  if (!deflt.useReferenceCounts) return;

#ifdef DEVELOPMENT_CODE
#ifdef SHOW_VALIDATE_CACHECOUNTS
  printf("Validating cache counts for %ld handles\n", getLastNode());
#endif
  const node_handle N = getLastNode()+1;

#ifdef SHOW_VALIDATE_CACHECOUNTS
  printf("  Counting...\n");
  fflush(stdout);
#endif

  std::vector <unsigned long> counts(N, 0);
  compute_table::countAllNodeEntries(this, counts);

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

#endif
}


void MEDDLY::forest::countNodesByLevel(long* active) const
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


void MEDDLY::forest::validateDownPointers(const unpacked_node &nb) const
{
    switch (getReductionRule()) {
        case reduction_rule::IDENTITY_REDUCED:
        case reduction_rule::FULLY_REDUCED:
            for (unsigned i=0; i<nb.getSize(); i++) {
                if (isTerminalNode(nb.down(i))) continue;
                MEDDLY_DCASSERT(!isDeletedNode(nb.down(i)));
                if (isLevelAbove(nb.getLevel(), getNodeLevel(nb.down(i))))
                    continue;

                FILE_output s(stdout);
                s << "Down pointer violation in created node at level "
                  << nb.getLevel() << ":\n";
                nb.show(s, true);
                s << "\nPointer " << nb.down(i) << " at level "
                  << getNodeLevel(nb.down(i)) << ":\n";
                showNode(s, nb.down(i), SHOW_DETAILS);
                MEDDLY_DCASSERT(0);
            }
            break;

        case reduction_rule::QUASI_REDUCED:
#ifdef DEVELOPMENT_CODE
            int nextLevel;
            if (isForRelations()) {
                nextLevel = downLevel(nb.getLevel());
            } else {
                nextLevel = nb.getLevel()-1;
            }
#endif
            for (unsigned i=0; i<nb.getSize(); i++) {
                if (nb.down(i)==getTransparentNode()) continue;
                MEDDLY_DCASSERT(getNodeLevel(nb.down(i)) == nextLevel);
            }
            break;

        default:
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    };  // switch
}

void MEDDLY::forest::showInfo(output &s, int verb)
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

void MEDDLY::forest::setLogger(logger* L, const char* name)
{
    theLogger = L;
    if (theLogger) theLogger->logForestInfo(this, name);
}

void MEDDLY::forest::reportForestStats(output &s, const char* pad) const
{
  // default - do nothing
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Methods for variable reordering
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::forest::reorderVariables(const int* level2var)
{
    removeAllComputeTableEntries();


  // Create a temporary variable order
  // Support in-place update and avoid interfering other forests
  var_order = std::make_shared<variable_order>(*var_order);

  auto reordering = reordering_factory::create(getPolicies().reorder);
  reordering->reorderVariables(this, level2var);

  var_order = getDomain()->makeVariableOrder(*var_order);
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
// *                                                                *
// *                         forest methods                         *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::policies MEDDLY::forest::mddDefaults;
MEDDLY::policies MEDDLY::forest::mxdDefaults;


MEDDLY::forest
::forest(domain* _d, bool rel, range_type t, edge_labeling ev,
  const policies &p, int* lrr) : nodeHeaders(*this, mstats, stats), deflt(p)
{
    // Set up domain
    d = _d;

    // Initialize variable order
    var_order = d->makeDefaultVariableOrder();

    isRelation = rel;
    rangeType = t;
    edgeLabel = ev;

    registerForest(this);

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
                for(unsigned i=1;i<=d->getNumVariables();i++)
                    lrr[i]=-2;

            else
                if (isFullyReduced())
                for(unsigned i=1;i<=d->getNumVariables();i++)
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

    for(unsigned i=1;i<=d->getNumVariables();i++)
        if(level_reduction_rule[i]==-3)              //isIdentityReduced()
         throw error(error::INVALID_POLICY, __FILE__, __LINE__);
  }

#ifdef ALLOW_DEPRECATED_0_17_4
    //
    // Initialize misc. protected data
    //
    terminalNodesStatus = MEDDLY::forest::ACTIVE;
#endif

    //
    // Initialize the root edges
    //
    roots = nullptr;

    //
    // Empty logger
    //
    theLogger = 0;

    //
    // Unique table(s)
    //
    unique = new unique_table(this);
    implUT = new impl_unique_table(this);

    //
    // Finish initializing nodeHeaders
    //
    nodeHeaders.initialize();

    if (deflt.useReferenceCounts) {
        reachable = nullptr;
    } else {
        reachable = new node_marker(this, &nodeHeaders);
        nodeHeaders.linkReachable(reachable->linkBits());
    }

    //
    // nodeMan is initialized in initializeStorage()
    //
    nodeMan = nullptr;

    //
    // Initialize node characteristics to defaults
    //
    unhashed_bytes = 0;
    hashed_bytes = 0;

    //
    // For debugging:
    //
    delete_depth = 0;
}

MEDDLY::forest::~forest()
{
#ifdef DEBUG_CLEANUP
    std::cout << "Deleting forest " << this << ", #" << FID() << "\n";
#endif
#ifdef REPORT_ON_DESTROY
    printf("Destroying forest.  Stats:\n");
    reportMemoryUsage(stdout, "\t", 9);
#endif

    unregisterDDEdges();

    // unique table
    delete unique;
    delete implUT;

    // node storage
    delete nodeMan;

    unregisterForest(this);
    ct_entry_type::invalidateAllWithForest(this);
    operation::destroyAllWithForest(this);
}

void MEDDLY::forest::markForDeletion()
{
    if (is_marked_for_deletion) return;
    is_marked_for_deletion = true;
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

MEDDLY::enumerator::iterator* MEDDLY::forest::makeFixedRowIter() const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

MEDDLY::enumerator::iterator* MEDDLY::forest::makeFixedColumnIter() const
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}


//



// ===================================================================
//
// Deprecated as of version 0.17.4
//
// ===================================================================

#ifdef ALLOW_DEPRECATED_0_17_4
MEDDLY::forest::node_status
MEDDLY::forest::getNodeStatus(MEDDLY::node_handle node) const
{
    if (isMarkedForDeletion()) {
        return MEDDLY::forest::DEAD;
    }
    if (isTerminalNode(node)) {
        return terminalNodesStatus;
    }
    if (isDeletedNode(node)) {
        return MEDDLY::forest::DEAD;
    }
    // Active node.

    // If we're using reference counts,
    // and the incoming count is zero,
    // then we must be using optimistic
    // and the node is stale but recoverable.

    // If we're NOT using reference counts,
    // since we're not a deleted node,
    // assume we are still active.

    if (deflt.useReferenceCounts) {
        if (getNodeInCount(node) == 0) {
            return MEDDLY::forest::RECOVERABLE;
        }
    }

    return MEDDLY::forest::ACTIVE;
}


// ******************************************************************
// *                     expert_forest  methods                     *
// ******************************************************************

MEDDLY::expert_forest::expert_forest(domain *d, bool rel, range_type t,
  edge_labeling ev, const policies &p, int* level_reduction_rule)
: forest(d, rel, t, ev, p, level_reduction_rule)
{
}

MEDDLY::expert_forest::~expert_forest()
{
}

#endif


// ===================================================================
//
// Deprecated as of version 0.17.3
//
// ===================================================================

#ifdef ALLOW_DEPRECATED_0_17_3

MEDDLY::node_handle*
MEDDLY::forest
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

  unpacked_node::Recycle(M);

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


unsigned long MEDDLY::forest::getNodeCount(node_handle p) const
{
    node_marker nm(this);
    nm.mark(p);
    return nm.countMarked();
}

unsigned long MEDDLY::forest::getNodeCount(const node_handle* roots, int N) const
{
    node_marker nm(this);
    for (int i=0; i<N; i++) {
        nm.mark(roots[i]);
    }
    return nm.countMarked();
}

unsigned long MEDDLY::forest::getEdgeCount(node_handle p, bool countZeroes) const
{
    node_marker nm(this);
    nm.mark(p);
    return countZeroes ? nm.countEdges() : nm.countNonzeroEdges();
}


void MEDDLY::forest
::showNodeGraph(output &s, const node_handle* p, int n) const
{
    node_marker M(this);
    for (int i=0; i<n; i++) {
        M.mark(p[i]);
    }
    M.showByLevels(s);
}


void MEDDLY::forest
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

void MEDDLY::forest
::writeEdges(output &s, const dd_edge* E, unsigned n) const
{
    //
    // Mark all reachable nodes
    //
    node_marker M(this);
    for (unsigned i=0; i<n; i++) {
        M.mark(E[i].getNode());
    }

    //
    // Build list of nodes to write, in bottom-up order.
    // When finished, output2handle[i] gives the handle
    // of the ith node to write.
    //
    vector <node_handle> output2handle;
    for (int k=1; k <= getNumVariables(); k++) {
        if (isForRelations()) {
            M.getNodesAtLevel(-k, output2handle);
        }

        M.getNodesAtLevel(k, output2handle);
    } // for k
    const unsigned num_nodes = output2handle.size();

    //
    // Build the inverse list.
    // When finished, handle2output[h] gives 0 if the node with handle h
    // will not be written; otherwise it gives the order it appears in
    // the output file.
    //
    node_handle maxnode = 0;
    for (unsigned i=0; i<output2handle.size(); i++) {
        UPDATEMAX(maxnode, output2handle[i]);
    }
    vector <unsigned> handle2output (1+maxnode);
    for (unsigned i=0; i<output2handle.size(); i++) {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1, output2handle[i], maxnode+1);
        handle2output[output2handle[i]] = i+1;
    }

#ifdef DEBUG_WRITE
    ostream_output out(std::cout);
    out << "Got list of nodes:\n";
    for (int i=0; output2handle[i]; i++) {
        if (i) out << ", ";
        out << output2handle[i];
    }
    out << "\n";
    out << "Got inverse list:\n";
    for (int i=0; i<=maxnode; i++) {
        if (i) out << ", ";
        out << handle2output[i];
    }
    out << "\n";
#endif

    //
    // Write the nodes
    //
    const char* block = codeChars();
    s << block << " " << num_nodes << "\n";

    unpacked_node* un = unpacked_node::New();
    for (unsigned i=0; i<output2handle.size(); i++) {
        s << getNodeLevel(output2handle[i]) << " ";
        unpackNode(un, output2handle[i], FULL_OR_SPARSE);
        un->write(s, handle2output);
    }
    unpacked_node::Recycle(un);

    // reverse the block
    for (int i=strlen(block); i; ) {
        i--;
        s.put(block[i]);
    }
    s.put('\n');

    //
    // Write the actual edge pointers
    //
    s << "ptrs " << n << "\n";
    for (unsigned i=0; i<n; i++) {
        s.put('\t');
        E[i].write(s, handle2output);
        s.put('\n');
    }
    s << "srtp\n";
}

void MEDDLY::forest::readEdges(input &s, dd_edge* E, unsigned n)
{
    try {
        const char* block = codeChars();

        s.stripWS();
        s.consumeKeyword(block);
        s.stripWS();
        unsigned num_nodes = unsigned(s.get_integer());
#ifdef DEBUG_READ_DD
        std::cerr << "Reading " << num_nodes << " nodes in " << block << " forest\n";
#endif
#ifdef DEBUG_READ
        std::cerr << "Reading " << num_nodes << " nodes in forest " << block << "\n";

#endif

        //
        // Build translation from file node# to forest node handles
        //
        std::vector <node_handle> map(1+num_nodes);

        for (unsigned node_index=1; node_index<=num_nodes; node_index++) {
            // Read this node

            //
            // read the level number
            //
            s.stripWS();
            int k = s.get_integer();
            if (!isValidLevel(k)) {
                throw error(error::INVALID_LEVEL, __FILE__, __LINE__);
            }

#ifdef DEBUG_READ
            std::cerr << "Reading: level " << k;
#endif

            //
            // read the node size (sparse/full)
            //
            s.stripWS();
            int rawsize = s.get_integer();
            int n;
            unpacked_node* nb = (rawsize < 0)
                ? unpacked_node::newSparse(this, k, n=-rawsize)
                : unpacked_node::newFull(this, k, n=rawsize);

#ifdef DEBUG_READ
            std::cerr << " rawsize " << rawsize << '\n';
#endif

            //
            // read the node
            //
            nb->read(s, map);

            //
            // Reduce the node, and update the translation
            //
            map[node_index] = createReducedNode(-1, nb);

#ifdef DEBUG_READ
            std::cerr << "File node " << node_index << " reduced to ";
            ostream_output myout(std::cerr);
            showNode(myout, map[node_index], SHOW_DETAILS | SHOW_INDEX);
            std::cerr << std::endl;
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
        std::cerr << "Done reading, expecting " << buffer << " keyword\n";
#endif

        //
        // match the reversed block
        //
        s.stripWS();
        s.consumeKeyword(buffer);
#ifdef DEBUG_READ
        std::cerr << "  got " << buffer << "\n";
#endif
#ifdef DEBUG_READ_DD
        std::cerr << "Finished " << buffer << " forest\n";
#endif

        //
        // Read the pointers
        //
        s.stripWS();
        s.consumeKeyword("ptrs");
#ifdef DEBUG_READ
        std::cerr << "Got ptrs\n";
#endif

        s.stripWS();
        unsigned num_ptrs = unsigned(s.get_integer());
#ifdef DEBUG_READ
        std::cerr << "Reading " << num_ptrs << " pointers\n";
#endif
#ifdef DEBUG_READ_DD
        std::cerr << "Reading " << num_ptrs << " pointers\n";
#endif
        MEDDLY_DCASSERT(num_ptrs <= n);
        if (num_ptrs > n) {
#ifdef DEBUG_READ
            std::cerr   << "Error at " << __FILE__ << ":" << __LINE__
                        << ", E[] is of size " << n
                        << ", needs to be at least " << num_ptrs << std::endl;
#endif
            throw error(error::INVALID_ASSIGNMENT, __FILE__, __LINE__);
        }
        for (unsigned i=0; i<num_ptrs; i++) {
            E[i].attach(this);
            E[i].read(s, map);
        }

        s.stripWS();
        s.consumeKeyword("srtp");

#ifdef DEBUG_READ_DD
        std::cerr << "Done reading pointers\n";
#endif

        //
        // unlink map pointers
        //
        for (unsigned i=1; i<=num_nodes; i++) {
            unlinkNode(map[i]);
        }

#ifdef DEVELOPMENT_CODE
        validateIncounts(true);
#endif
    } // try
    catch (error& e) {
#ifdef DEBUG_READ
        std::cerr << "Read failed (error: " << e.getName() << ")\n";
        std::cerr << "Next few characters of file:\n";
        for (unsigned i=0; i<20; i++) {
            int c = s.get_char();
            if (EOF == c) {
                std::cerr << "EOF";
                break;
            }
            std::cerr << char(c) << " (" << c << ") ";
        }
        std::cerr << std::endl;
#endif
        throw e;
    }
}


const char* MEDDLY::forest::codeChars() const
{
  return "unknown dd";
}

#endif


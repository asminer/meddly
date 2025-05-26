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
#include <cmath>  // for round()

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


// #define DEBUG_ADDRESS_RESIZE
// #define DEBUG_CREATE_REDUCED
// #define DEBUG_GC
// #define DEBUG_WRITE
// #define DEBUG_READ

// #define TRACK_DELETIONS


// Thoroughly check reference counts.
// Very slow.  Use only for debugging.
// #define VALIDATE_INCOUNTS_ON_DELETE

// Check reference counts after operations.
// Very slow. Use only for debugging.
// #define ACTUALLY_VALIDATE_INCOUNTS

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
            // int* level_reduction_rule,
            int tv)
{
    switch (el) {
        case edge_labeling::MULTI_TERMINAL:
            switch (t) {
                case range_type::BOOLEAN:
                    if (sr)
                        return new mt_mxd_bool(d, p, tv);
                    else
                        return new mt_mdd_bool(d, p, tv);

                case range_type::INTEGER:
                    if (sr)
                        return new mt_mxd_int(d, p, tv);
                    else
                        return new mt_mdd_int(d, p, tv);

                case range_type::REAL:
                    if (sr)
                        return new mt_mxd_real(d, p, (float)tv);
                    else
                        return new mt_mdd_real(d, p, (float)tv);

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
                return new evmxd_pluslong(d, p);
            else
                return new evmdd_pluslong(d, p);


        case edge_labeling::INDEX_SET:
            if (range_type::INTEGER != t || sr) {
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            }
            return new evmdd_index_set_long(d, p);


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

#ifdef ALLOW_DEPRECATED_0_17_8

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

#endif

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Methods to add/remove nodes
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MEDDLY::relation_node* MEDDLY::forest
::buildImplicitNode(node_handle rnh)
{
    return implUT->getNode(rnh);
}

namespace MEDDLY {

    //
    // Template helper: normalizing EV+
    //
    template <typename EDGETYPE>
    inline void normalize_evplus(unpacked_node &un, edge_value &ev,
            unsigned &nnz)
        {
            nnz = 0;
            ev.set(EDGETYPE(0));
            EDGETYPE minval = 0;
            //
            // First scan: find smallest value
            //
            for (unsigned i=0; i<un.getSize(); i++) {
                if (0 == un.down(i)) {
                    un.edgeval(i) = ev;
                    continue;
                }
                EDGETYPE uni;
                un.edgeval(i).get(uni);
                if (nnz) {
                    if (uni < minval)
                    {
                        minval = uni;
                    }
                } else {
                    minval = uni;
                }
                ++nnz;
            }
            //
            // Second scan: adjust; unless adjustment is zero
            //
            if (minval) {
                for (unsigned i=0; i<un.getSize(); i++) {
                    if (un.down(i)) {
                        un.edgeval(i).subtract(minval);
                    }
                }
            }
            ev.set(minval);
        }

    //
    // Template helper: normalizing EV*
    //
    template <typename EDGETYPE>
    inline void normalize_evstar(unpacked_node &un, edge_value &ev,
            unsigned &nnz)
        {
            nnz = 0;
            ev.set(EDGETYPE(0));
            EDGETYPE firstval = 0;
            //
            // First scan: find the first non-zero edge
            //
            for (unsigned i=0; i<un.getSize(); i++) {
                if (0 == un.down(i)) {
                    un.edgeval(i) = ev;
                    continue;
                }
                EDGETYPE uni;
                un.edgeval(i).get(uni);
                MEDDLY_DCASSERT(uni);
                //
                // Actual nonzero edge
                //
                if (!nnz) {
                    firstval = uni;
                }
                ++nnz;
            }
            //
            // Second scan: adjust; unless adjustment is zero or one
            //
            if (firstval && (firstval != 1)) {
                for (unsigned i=0; i<un.getSize(); i++) {
                    un.edgeval(i).divide(firstval);
                }
            }
            ev.set(firstval);
        }

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
    unsigned nnz = 0;
    switch (edgeLabel) {
        case edge_labeling::MULTI_TERMINAL:
                ev.set();
                if (termprec) {
                    // Round any terminal values based on termprec
                    // AND count number of nonzeroes.
                    for (unsigned i=0; i<un->getSize(); i++) {
                        if (un->down(i) == 0) {
                            // Terminal 0
                            continue;
                        }
                        ++nnz;
                        if (un->down(i) > 0) {
                            // Nonterminal
                            continue;
                        }
                        // Terminal that needs rounding
                        terminal T(terminal_type::REAL, un->down(i));

                        T.adjustReal(
                                round(T.getReal() / termprec) * termprec
                        );

                        un->down(i) = T.getHandle();
                    }
                } else {
                    // Just count nonzeroes
                    for (unsigned i=0; i<un->getSize(); i++) {
                        if (un->down(i)) ++nnz;
                    }
                }
                break;

        case edge_labeling::EVPLUS:
        case edge_labeling::INDEX_SET:
                MEDDLY_DCASSERT(isRangeType(range_type::INTEGER));

                normalize_evplus<long>(*un, ev, nnz);

                /*
                // TBD: allow other edge types other than long?
                MEDDLY_DCASSERT(edge_type::LONG == the_edge_type);
                ev.set(0L);
                for (unsigned i=0; i<un->getSize(); i++) {
                    if (0 == un->down(i)) {
                        un->setEdgeval(i, 0L);
                        continue;
                    }
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
                        if (un->down(i)) {
                            un->subtractFromEdge(i, ev.getLong());
                        }
                    }
                }
                */
                break;

        case edge_labeling::EVTIMES:
                MEDDLY_DCASSERT(isRangeType(range_type::REAL));

                normalize_evstar<float>(*un, ev, nnz);

                /*
                // TBD: allow other edge types other than float?
                MEDDLY_DCASSERT(edge_type::FLOAT == the_edge_type);
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
                */
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
    unpacked_node* key = unpacked_node::newFromNode(this, node, SPARSE_ONLY);
    key->computeHash();
    if (key->hash() != un->hash())
    {
        FILE_output s(stderr);
        s << "Hash mismatch\n";
        s << "Original node: ";
        un->show(s, true);
        s << "\n";
#ifdef DEBUG_UNPACKED_HASH
        un->debugHash(s);
#endif
        s << "Sparse version: ";
        key->show(s, true);
        s << "\n";
#ifdef DEBUG_UNPACKED_HASH
        key->debugHash(s);
#endif
    }
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
        unpacked_node* key = unpacked_node::newFromNode(this, p, SPARSE_ONLY);
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
        validateIncounts(false, __FILE__, __LINE__, "delete");
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
    unpacked_node* U = nullptr;

    //
    // Multi-terminal
    //

    // Special case for identity reduced:
    //   if p is at a primed level,
    //   and K is the same level as p,
    //   then we need to take care if p is a singleton node.
    //
    bool check_singleton = isIdentityReduced() && (K<0) && (getNodeLevel(p) == K);

    while (K != L) {
        if (isForRelations()) {
            K = MXD_levels::upLevel(K);
        } else {
            K = MDD_levels::upLevel(K);
        }

        if (check_singleton) {
            unsigned sind;
            node_handle sdwn;
            if (isSingletonNode(p, sind, sdwn)) {
                U = unpacked_node::newWritable(this, K, FULL_ONLY);

                // unsigned size = getLevelSize(K);
                // U = unpacked_node::newFull(this, K, size);
                for (unsigned i=0; i<U->getSize(); i++) {
                    U->setFull(i, noop_edge, linkNode( (i==sind) ? sdwn : p ));
                }
                unlinkNode(p);

                edge_value ev;
                createReducedNode(U, ev, p);
                MEDDLY_DCASSERT(ev == noop_edge);
                check_singleton = false;
                continue;
            }
            check_singleton = false;
        }
        if (!isIdentityReduced() || K<0)
        {
            /*
               Build a redundant node.
               Necessary for quasi reduced always,
               and for identity reduced at primed levels
               (unprimed is fully reduced).
               */
            U = unpacked_node::newRedundant(this, K, noop_edge, p, FULL_ONLY);
            linkAllDown(*U, 1);
            edge_value ev;
            createReducedNode(U, ev, p);
            MEDDLY_DCASSERT(ev == noop_edge);
        }
    }

    return p;
}


MEDDLY::node_handle
MEDDLY::forest::_makeIdentitiesTo(node_handle p, int K, int L, int in)
{
    MEDDLY_DCASSERT(!isIdentityReduced());
    MEDDLY_DCASSERT(isForRelations());

    MEDDLY_DCASSERT(L!=0);
    unpacked_node* Uun;
    unpacked_node* Upr;
    edge_value ev;
    //
    // Multi-terminal
    //

    if (K<0) {
        //
        // Add a redundant layer
        //
        if ( !isFullyReduced() ) {
            Uun = unpacked_node::newRedundant(this, -K, noop_edge, p, FULL_ONLY);
            linkAllDown(*Uun, 1);
            createReducedNode(Uun, ev, p);
            MEDDLY_DCASSERT(ev == noop_edge);
        }
        K = MXD_levels::upLevel(K);
    }

    MEDDLY_DCASSERT(K>=0);
    const int Lstop = (L<0) ? MXD_levels::downLevel(L) : L;

    //
    // Proceed in unprimed, primed pairs
    //
    for (K++; K<=Lstop; K++) {
        Uun = unpacked_node::newWritable(this, K, FULL_ONLY);
        // Uun = unpacked_node::newFull(this, K, getLevelSize(K));

        // build primed level nodes
        for (unsigned i=0; i<Uun->getSize(); i++) {
            Upr = unpacked_node::newWritable(this, -K, 1, SPARSE_ONLY);
            // Upr = unpacked_node::newSparse(this, -K, 1);
            Upr->setSparse(0, i, noop_edge, linkNode(p));
            node_handle h;
            createReducedNode(Upr, ev, h);
            MEDDLY_DCASSERT(ev == noop_edge);
            Uun->setFull(i, noop_edge, h);
        }

        unlinkNode(p);
        createReducedNode(Uun, ev, p);
        MEDDLY_DCASSERT(ev == noop_edge);
    } // for k

    //
    // Add top identity node, if L is negative
    //
    if (L<0) {
        MEDDLY_DCASSERT(-K == L);
        Upr = unpacked_node::newWritable(this, L, 1, SPARSE_ONLY);
        MEDDLY_DCASSERT(in>=0);
        Upr->setSparse(0, unsigned(in), noop_edge, p);
        createReducedNode(Upr, ev, p);
        MEDDLY_DCASSERT(ev == noop_edge);
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
    unpacked_node* key = unpacked_node::newFromNode(this, p, SPARSE_ONLY);
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
//
// value to edge
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::forest::getEdgeForValue(rangeval T, edge_value &v, node_handle &p)
    const
{
    if (!T.hasType(rangeType)) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
    switch (edgeLabel) {

        case edge_labeling::MULTI_TERMINAL:
                if (!T.isNormal()) {
                    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
                }

                v.set();
                switch (rangeType) {
                    case range_type::BOOLEAN:
                    {
                        terminal t = terminal(bool(T));
                        p = t.getHandle();
                        return;
                    }

                    case range_type::INTEGER:
                    {
                        terminal t = terminal(long(T));
                        p = t.getHandle();
                        return;
                    }

                    case range_type::REAL:
                    {
                        terminal t = terminal(double(T));
                        p = t.getHandle();
                        return;
                    }

                    default:
                        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

                } // switch
                // shouldn't get here
                MEDDLY_DCASSERT(false);
                return;

        case edge_labeling::INDEX_SET:
        case edge_labeling::EVPLUS:
                if (T.isPlusInfinity()) {
                    v.setTempl(the_edge_type, 0);
                    p = OMEGA_INFINITY;
                    return;
                }
                if (!T.isNormal()) {
                    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
                }
                p = OMEGA_NORMAL;

                switch (the_edge_type) {
                    case edge_type::INT:
                        v.set(int(T));
                        return;

                    case edge_type::LONG:
                        v.set(long(T));
                        return;

                    default:
                        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

                } // switch
                // shouldn't get here
                MEDDLY_DCASSERT(false);
                return;


        case edge_labeling::EVTIMES:
                if (!T.isNormal()) {
                    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
                }

                switch (the_edge_type) {
                    case edge_type::FLOAT:
                        v.set(float(T));
                        p = (0.0 == v.getFloat()) ? OMEGA_ZERO : OMEGA_NORMAL;
                        return;

                    case edge_type::DOUBLE:
                        v.set(double(T));
                        p = (0.0 == v.getDouble()) ? OMEGA_ZERO : OMEGA_NORMAL;
                        return;

                    default:
                        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

                } // switch
                // shouldn't get here
                MEDDLY_DCASSERT(false);
                return;

        default:
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// edge to value
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::forest::getValueForEdge(const edge_value &v, node_handle p,
        rangeval &T) const
{
    if (!isTerminalNode(p)) {
        throw error(error::INVALID_LEVEL, __FILE__, __LINE__);
    }
    switch (edgeLabel) {

        case edge_labeling::MULTI_TERMINAL:
            {
                MEDDLY_DCASSERT(v.isVoid());
                terminal t(the_terminal_type, p);
                switch (the_terminal_type) {
                    case terminal_type::BOOLEAN:
                    {
                        T = rangeval(t.getBoolean());
                        return;
                    }

                    case terminal_type::INTEGER:
                    {
                        T = rangeval(t.getInteger());
                        return;
                    }

                    case terminal_type::REAL:
                    {
                        T = rangeval(t.getReal());
                        return;
                    }

                    default:
                        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

                } // switch
                  // shouldn't get here
                MEDDLY_DCASSERT(false);
                return;
            }

        case edge_labeling::INDEX_SET:
        case edge_labeling::EVPLUS:
                if (OMEGA_INFINITY == p) {
                    T = rangeval(range_special::PLUS_INFINITY,
                                    range_type::INTEGER);
                    return;
                }
                MEDDLY_DCASSERT(OMEGA_NORMAL == p);
                switch (the_edge_type) {
                    case edge_type::INT:
                        T = rangeval(v.getInt());
                        return;

                    case edge_type::LONG:
                        T = rangeval(v.getLong());
                        return;

                    default:
                        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

                } // switch
                // shouldn't get here
                MEDDLY_DCASSERT(false);
                return;


        case edge_labeling::EVTIMES:
                if (OMEGA_ZERO == p) {
                    T = rangeval(0.0);
                    return;
                }
                MEDDLY_DCASSERT(OMEGA_NORMAL == p);

                switch (the_edge_type) {
                    case edge_type::FLOAT:
                        T = rangeval(v.getFloat());
                        return;

                    case edge_type::DOUBLE:
                        T = rangeval(v.getDouble());
                        return;

                    default:
                        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

                } // switch
                // shouldn't get here
                MEDDLY_DCASSERT(false);
                return;

        default:
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    }
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

void MEDDLY::forest::showEdge(output &s, const edge_value &ev, node_handle d)
    const
{
    if (isMultiTerminal()) {
            if (d>0) {
                s.put('#');
                s.put(d);
            } else {
                terminal t(the_terminal_type, d);
                t.show(s);
            }
            return;
    }

    //
    // Edge valued
    //

    if (d == 0) {
        if (isEVTimes()) {
            s.put("<0, zero>");
        } else {
            s.put("<oo, w>");
        }
    } else {
        s.put('<');
        ev.show(s);
        s.put(", ");
        if (d < 0) {
            s.put('w');
        } else {
            s.put('#');
            s.put(d);
        }
        s.put('>');
    }
}


void MEDDLY::forest::showHeaderInfo(output &s, const unpacked_node &u) const
{
    if (isIndexSet()) {
        s.put(" card: ");
        s.put(static_cast<const long*>(u.UHptr())[0]);
    }
}

void MEDDLY::forest::writeHeaderInfo(output &s, const unpacked_node &u) const
{
    if (isIndexSet()) {
        s.put('\t');
        s.put(static_cast<const long*>(u.UHptr())[0]);
        s.put('\n');
    }
}

void MEDDLY::forest::readHeaderInfo(input &s, unpacked_node &u) const
{
    if (isIndexSet()) {
        long card = s.get_integer();
        u.setUHdata(&card);
#ifdef DEBUG_READ_DD
        std::cerr << "    got cardinality " << card << "\n";
#endif
    }
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
    // unreduced_node::initForest(f);
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
    // unreduced_node::doneForest(f);
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
  unpacked_node* un = unpacked_node::newFromNode(this, p, FULL_OR_SPARSE);
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

void MEDDLY::forest::validateIncounts(bool exact, const char* FN, unsigned LN,
        const char* opname) const
{
#ifndef ACTUALLY_VALIDATE_INCOUNTS
    return;
#endif

    static int idnum = 0;
    idnum++;

    // Inspect every active node's down pointers to determine
    // the incoming count for every active node.

    node_handle sz = getLastNode() + 1;
    std::vector <unsigned> in_validate(sz);

    unpacked_node* P = unpacked_node::New(this, SPARSE_ONLY);
    for (node_handle i = 1; i < sz; ++i) {
        if (!isActiveNode(i)) continue;
        P->initFromNode(i);

        // add to reference counts
        for (unsigned z=0; z<P->getSize(); z++) {
            if (isTerminalNode(P->down(z))) continue;
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, P->down(z), sz);
            in_validate[P->down(z)]++;
        }
    } // for i
    unpacked_node::Recycle(P);

    // Add counts for registered dd_edges
    for (const dd_edge* r = roots; r; r=r->next) {
        if (isTerminalNode(r->getNode())) continue;
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, r->getNode(), sz);
        in_validate[r->getNode()]++;
    }

    //
    // Add counts for down pointers in unpacked_nodes
    // in the build list for this forest.
    //
    unpacked_node::AddToIncomingCounts(this, in_validate);

    // Validate the incoming count stored with each active node using the
    // in_count array computed above
    for (node_handle i = 1; i < sz; ++i) {
        MEDDLY_DCASSERT(!isTerminalNode(i));
        if (!isActiveNode(i)) continue;
        bool fail = exact
            ?  in_validate[i] != getNodeInCount(i)
            :  in_validate[i] >  getNodeInCount(i);

        if (fail) {
            FILE_output fout(stdout);
            fout << "Validation #" << idnum << " failed for\n";
            fout << "\tnode " << i << "\n";
            fout << "\tnode's count " << getNodeInCount(i) << "\n";
            fout << "\tactual count " << in_validate[i] << "\n";
            dump(fout, SHOW_DETAILS);
            fout << "Requested from " << FN << " line " << LN;
            if (opname) fout << " operation " << opname;
            fout << '\n';
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
            // Check for identity rule violations
            for (unsigned i=0; i<nb.getSize(); i++) {
                unsigned out = nb.isSparse() ? nb.index(i) : i;
                unsigned index;
                node_handle down;
                const int dnlvl = getNodeLevel(nb.down(i));
                if (dnlvl>0) continue;
                if (!isSingletonNode(nb.down(i), index, down)) continue;
                if ((dnlvl != -nb.getLevel()) || (index == out))
                {
                    FILE_output s(stdout);
                    s   << "Identity violation in node created at level "
                        << nb.getLevel() << ":\n";
                    nb.show(s, true);
                    s   << "\nPointer " << out << " to " << nb.down(i)
                        << " at level " << dnlvl << ":\n";
                    showNode(s, nb.down(i), SHOW_DETAILS);
                    s   << "\n";
                    if (dnlvl != -nb.getLevel()) {
                        s << "Jumps too far, from " << nb.getLevel()
                          << " to " << dnlvl << "\n";
                    }
                    if (index == i) {
                        s << "Illegal edge into singleton\n";
                    }
                    s.flush();
                    MEDDLY_DCASSERT(0);
                }
            }
            // NO break here; we need to fall through to below

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
                nextLevel = MXD_levels::downLevel(nb.getLevel());
            } else {
                nextLevel = MDD_levels::downLevel(nb.getLevel());
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


MEDDLY::forest::forest(domain* _d, bool rel, range_type t, edge_labeling ev,
  const policies &p) : nodeHeaders(*this, mstats, stats), deflt(p)
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

    //
    // Other defaults
    //
    setTerminalPrecision(0);
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

    // free(level_reduction_rule);

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

#ifdef ALLOW_DEPRECATED_0_17_9

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

#endif

MEDDLY::node_handle MEDDLY::forest::unionOneMinterm(node_handle a,  int* from,  int* to, int level)
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::forest::createEdgeForVar(int vh, bool pr,
            const rangeval* terms, dd_edge& result)
{
    //
    // Sanity checks
    //
    if (vh < 0 || vh > getNumVariables())
        throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);
    if (!result.isAttachedTo(this))
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    if (!isForRelations() && pr)
        throw error(error::INVALID_ASSIGNMENT, __FILE__, __LINE__);

    const int level = getLevelByVar(vh);


    //
    // Get info for node we're building
    //
    const int k = pr ? -level : level;
    const int km1 = isForRelations()
        ? MXD_levels::downLevel(k)
        : MDD_levels::downLevel(k)
    ;

    //
    // Make the node
    //
    edge_value ev;
    node_handle node;
    unpacked_node* nb = unpacked_node::newWritable(this, k, FULL_ONLY);

    if (terms) {
        for (unsigned i=0; i<nb->getSize(); i++) {
            if ( ! terms[i].hasType(rangeType) )
            {
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            }
            getEdgeForValue(terms[i], ev, node);
            nb->setFull(i, ev, makeRedundantsTo(node, 0, km1) );
        }
    } else {
        for (unsigned i=0; i<nb->getSize(); i++) {
            switch (rangeType) {
                case range_type::INTEGER:
                    getEdgeForValue(long(i), ev, node);
                    break;
                case range_type::REAL:
                    getEdgeForValue(double(i), ev, node);
                    break;
                default:
                    getEdgeForValue(bool(i), ev, node);
                    break;
            }
            nb->setFull(i, ev, makeRedundantsTo(node, 0, km1) );
        }
    }

    //
    // Reduce, add redundants if necessary
    //
    createReducedNode(nb, ev, node);
    node = makeRedundantsTo(node, k, getNumVariables());
    result.set(ev, node);
}


// ===================================================================
//
// Deprecated as of version 0.17.7
//
// ===================================================================

#ifdef ALLOW_DEPRECATED_0_17_7

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

void MEDDLY::forest::createEdge(double val, dd_edge &e)
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
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
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

#endif


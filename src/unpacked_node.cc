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

//
// TBD: add zero-out range [low, high)
// TBD: on node recycle, set down + edge to transparent on [low, high)
//

// TODO: Testing

#include <map>
#include "defines.h"
#include "unpacked_node.h"
#include "forest.h"
#include "hash_stream.h"
#include "terminal.h"
#include "node_marker.h"
#include "operators.h"

// #define DEBUG_FORLISTS

// ******************************************************************
// *                     unpacked_lists  struct                     *
// ******************************************************************

struct MEDDLY::unpacked_lists {
    unpacked_node* free_full_only;
    unpacked_node* free_full_or_sparse;
    unpacked_node* free_sparse_only;
    unpacked_node* building;

    inline void init() {
        free_full_only = nullptr;
        free_full_or_sparse = nullptr;
        free_sparse_only = nullptr;
        building = nullptr;
    }

    inline unpacked_node* popFree(node_storage_flags ns) {
        unpacked_node* x;
        switch (ns) {
            case FULL_ONLY:
                x = free_full_only;
                if (x) free_full_only = x->next;
                return x;

            case SPARSE_ONLY:
                x = free_sparse_only;
                if (x) free_sparse_only = x->next;
                return x;

            case FULL_OR_SPARSE:
                x = free_full_or_sparse;
                if (x) free_full_or_sparse = x->next;
                return x;

            default:
                MEDDLY_DCASSERT(false);
        }
        return nullptr;
    }
    inline void pushFree(unpacked_node* x, node_storage_flags ns) {
        MEDDLY_DCASSERT(x);
        switch (ns) {
            case FULL_ONLY:
                x->next = free_full_only;
                free_full_only = x;
                return;

            case SPARSE_ONLY:
                x->next = free_sparse_only;
                free_sparse_only = x;
                return;

            case FULL_OR_SPARSE:
                x->next = free_full_or_sparse;
                free_full_or_sparse = x;
                return;

            default:
                MEDDLY_DCASSERT(false);
        }
    }
};


// ******************************************************************
// *                     unpacked_node  statics                     *
// ******************************************************************

MEDDLY::unpacked_lists* MEDDLY::unpacked_node::ForLists;
unsigned MEDDLY::unpacked_node::ForListsAlloc;

// ******************************************************************
// *                                                                *
// *                                                                *
// *                     unpacked_node  methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::unpacked_node::unpacked_node(const forest* f, node_storage_flags fs)
    : parent(f),
        pFID(f->FID()),
        extra_unhashed_size(f->unhashedHeaderBytes()),
        extra_hashed_size(f->hashedHeaderBytes()),
        nodestor(fs),
        the_edge_type(f->getEdgeType())
{
    MEDDLY_DCASSERT(f);

    next = nullptr;
    prev = nullptr;

    modparent = nullptr;

    _down = nullptr;
    _index = nullptr;
    _edge = nullptr;

    mark_extra = 0;

    alloc = 0;
    size = 0;

    level = 0;

#ifdef DEVELOPMENT_CODE
    can_be_recycled = false;
    has_hash = false;
#endif

    //
    // Allocate extra headers if needed
    //

    if (extra_unhashed_size) {
        extra_unhashed = malloc(extra_unhashed_size);
        if (!extra_unhashed) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
    }
    if (extra_hashed_size) {
        extra_hashed = malloc(extra_hashed_size);
        if (!extra_hashed) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
    }

    orig_was_identity = false;
}

#ifdef ALLOW_DEPRECATED_0_17_8

MEDDLY::unpacked_node::unpacked_node(const forest* f)
    : nodestor(FULL_OR_SPARSE)
{
    next = nullptr;
    prev = nullptr;

    parent = nullptr;
    modparent = nullptr;
    pFID = 0;

    _down = nullptr;
    _index = nullptr;
    _edge = nullptr;

    mark_extra = 0;

    alloc = 0;
    size = 0;

    extra_unhashed = nullptr;
    extra_unhashed_size = 0;

    extra_hashed = nullptr;
    extra_hashed_size = 0;

    level = 0;

#ifdef ALLOW_EXTENSIBLE
    can_be_extensible = false;
    is_extensible = false;
#endif

#ifdef DEVELOPMENT_CODE
    can_be_recycled = false;
    has_hash = false;
#endif

    if (f) attach(f);
}

#endif

MEDDLY::unpacked_node::~unpacked_node()
{
    free(_down);
    free(_index);
    free(_edge);

    free(extra_unhashed);
    free(extra_hashed);
}

#ifdef ALLOW_DEPRECATED_0_17_8

void MEDDLY::unpacked_node::attach(const forest* f)
{
    MEDDLY_DCASSERT(!parent);
    MEDDLY_DCASSERT(f);

    parent = f;
    modparent = nullptr;
    pFID = f->FID();
    the_edge_type = parent->getEdgeType();

    //
    // Allocate extra headers if needed
    //

    extra_unhashed_size =   f->unhashedHeaderBytes();
    extra_hashed_size =     f->hashedHeaderBytes();

    if (extra_unhashed_size) {
        extra_unhashed = malloc(extra_unhashed_size);
        if (!extra_unhashed) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
    }
    if (extra_hashed_size) {
        extra_hashed = malloc(extra_hashed_size);
        if (!extra_hashed) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
    }

    orig_was_identity = false;
}

#endif

void MEDDLY::unpacked_node::initFromNode(node_handle node)
{
    level = parent->getNodeLevel(node);
    MEDDLY_DCASSERT(level != 0);
    resize(parent->getLevelSize(level));
    MEDDLY_DCASSERT(parent->getNodeAddress(node));
    MEDDLY_DCASSERT(parent->getNodeManager());
    parent->getNodeManager()->fillUnpacked(*this, parent->getNodeAddress(node), nodestor);
}


void MEDDLY::unpacked_node::initRedundant(int k, node_handle node)
{
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(parent->isTerminalNode(node) || !parent->isDeletedNode(node));
#ifdef ALLOW_EXTENSIBLE
    is_extensible = parent->isExtensibleLevel(k);
    resize( is_extensible ? 1 : unsigned(parent->getLevelSize(k)) );
#else
    resize( unsigned(parent->getLevelSize(k)) );
#endif
    level = k;

    if (is_full) {
        for (unsigned i=0; i<getSize(); i++) {
            setFull(i, node);
        }
    } else {
        for (unsigned i=0; i<getSize(); i++) {
            setSparse(i, i, node);
        }
    }

    orig_was_identity = false;
}

void MEDDLY::unpacked_node::initRedundant(int k, const edge_value &ev,
        node_handle node)
{
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(parent->isTerminalNode(node) || !parent->isDeletedNode(node));
#ifdef ALLOW_EXTENSIBLE
    is_extensible = parent->isExtensibleLevel(k);
    resize( is_extensible ? 1 : unsigned(parent->getLevelSize(k)) );
#else
    resize( unsigned(parent->getLevelSize(k)) );
#endif
    level = k;

    if (is_full) {
        for (unsigned i=0; i<getSize(); i++) {
            setFull(i, ev, node);
        }
    } else {
        for (unsigned i=0; i<getSize(); i++) {
            setSparse(i, i, ev, node);
        }
    }

    orig_was_identity = false;
}



void MEDDLY::unpacked_node::initIdentity(int k, unsigned i, node_handle node)
{
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(parent->isTerminalNode(node) || !parent->isDeletedNode(node));
    level = k;

    if (is_full) {
        resize(parent->getLevelSize(k));
        clear(0, getSize());
        setFull(i, node);
    } else {
        resize(1);
        setSparse(0, i, node);
    }

    orig_was_identity = true;
}

void MEDDLY::unpacked_node::initIdentity(int k, unsigned i,
        const edge_value &ev, node_handle node)
{
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(parent->isTerminalNode(node) || !parent->isDeletedNode(node));
    level = k;

    if (is_full) {
        resize(parent->getLevelSize(k));
        clear(0, getSize());
        setFull(i, ev, node);
    } else {
        resize(1);
        setSparse(0, i, ev, node);
    }

    orig_was_identity = true;
}

/*
    Build blank writable nodes
*/

MEDDLY::unpacked_node* MEDDLY::unpacked_node::newWritable(forest* f, int lvl,
        node_storage_flags fs)
{
    MEDDLY_DCASSERT(f);
    unsigned tsz = f->getLevelSize(lvl);

    //
    // Same as the other newWritable() below here
    //
    unpacked_node* U = New(f, fs);
    MEDDLY_DCASSERT(U);
    U->level = lvl;
    U->resize(tsz);
    U->is_full = (fs != SPARSE_ONLY);
    if (U->is_full) {
        U->clear(0, tsz);
    }
    U->allowWrites(f);
    return U;
}

MEDDLY::unpacked_node* MEDDLY::unpacked_node::newWritable(forest* f, int lvl,
        unsigned tsz, node_storage_flags fs)
{
    unpacked_node* U = New(f, fs);
    MEDDLY_DCASSERT(U);
    U->level = lvl;
    U->resize(tsz);
    U->is_full = (fs != SPARSE_ONLY);
    if (U->is_full) {
        U->clear(0, tsz);
    }
    U->allowWrites(f);
    return U;
}

/*
    Display methods
*/

void MEDDLY::unpacked_node::show(output &s, bool details) const
{
    if (details) {
        s << "level: " << level;
        s << (isSparse() ? " nnzs: " : " size: ") << long(getSize());
#ifdef ALLOW_EXTENSIBLE
        if (isExtensible()) s.put('*');
#endif
        s.put(' ');
    }
    s << "down: " << (isSparse() ? '(' : '[');

    for (unsigned z=0; z<getSize(); z++) {
        if (isSparse()) {
            if (z) s << ", ";
            s << (unsigned long) index(z) << ":";
        } else {
            if (z) s.put('|');
        }
        if (hasEdges()) {
            parent->showEdge(s, edgeval(z), down(z));
        } else {
            edge_value none;
            parent->showEdge(s, none, down(z));
        }
    }

#ifdef ALLOW_EXTENSIBLE
    if (isExtensible()) s.put('*');
#endif

    s.put( isSparse() ? ')' : ']' );

    // show extra header stuff
    parent->showHeaderInfo(s, *this);
}


void MEDDLY::unpacked_node::write(output &s, const std::vector <unsigned> &map)
    const
{
    MEDDLY_DCASSERT(parent);

    if (isSparse()) {
        s.put(-long(getSize()));
    } else {
        s.put(long(getSize()));
    }

    //
    // write indexes (sparse only)
    //
    if (isSparse()) {
        s.put('\n');
        s.put('\t');
        for (unsigned z=0; z<getSize(); z++) {
            s.put(' ');
            s.put((unsigned long) index(z));
        }
    }

    //
    // write down pointers
    //
    s.put('\n');
    s.put('\t');
    for (unsigned z=0; z<getSize(); z++) {
        s.put(' ');
        const long d = down(z);
        if (d <= 0) {
            // terminal
            terminal t(parent->getTerminalType(), d);
            t.write(s);
        } else {
            // non-terminal
            s.put("n ");
            s.put(map[(unsigned long) d]);
        }
    }

    //
    // write edge values, if any
    //
    if (hasEdges()) {
        s.put('\n');
        s.put('\t');
        for (unsigned z=0; z<getSize(); z++) {
            s.put(' ');
            edgeval(z).write(s);
        }
    }
    s.put('\n');


    // write extra header stuff, should no-op if there isn't any
    parent->writeHeaderInfo(s, *this);
}


void MEDDLY::unpacked_node::read(input &s, const std::vector<node_handle> &map)
{
#ifdef DEBUG_READ_DD
    std::cerr << "  in unpacked_node::read\n";
#endif
    MEDDLY_DCASSERT(modparent);
    //
    // We already read the size
    //

    //
    // read indexes (sparse only)
    //
    if (isSparse()) {
#ifdef DEBUG_READ_DD
        std::cerr << "    reading indexes\n";
#endif
        for (unsigned z=0; z<getSize(); z++) {
            s.stripWS();
            long ndx = s.get_integer();
            if (ndx < 0) {
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            _index[z] = unsigned(ndx);
        }
    }

    //
    // read down pointers
    //
#ifdef DEBUG_READ_DD
    std::cerr << "    reading " << getSize() << " down pointers\n";
#endif
    for (unsigned z=0; z<getSize(); z++) {
        s.stripWS();
        int c = s.get_char();
        if ('n' == c) {
            // non-terminal
            node_handle d = s.get_integer();
#ifdef DEBUG_READ_DD
            std::cerr << "\tgot nonterm " << d << "\n";
#endif
            if (d < 0) {
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            _down[z] = modparent->linkNode(map[(unsigned long) d]);
        } else {
            // terminal
            s.unget(c);
            terminal t;
            t.read(s);
            _down[z] = t.getHandle();
        }
    }

    //
    // read edge values, if any
    //
    if (hasEdges()) {
#ifdef DEBUG_READ_DD
        std::cerr << "    reading " << stop << " edge values\n";
#endif
        for (unsigned z=0; z<getSize(); z++) {
            s.stripWS();
            _edge[z].read(s);
        }
    }

    //
    // read extra header stuff, if any
    //
    parent->readHeaderInfo(s, *this);

#ifdef DEBUG_READ_DD
    std::cerr << "  done unpacked_node::read\n";
#endif
}

//
// Helpers for checking in a node
//

void MEDDLY::unpacked_node::computeHash()
{
    MEDDLY_DCASSERT(!has_hash);
    trim();

    hash_stream s;
    s.start(0);

    if (extra_hashed_size) {
        s.push(extra_hashed, extra_hashed_size);
    }

    if (isSparse()) {
        if (parent->areEdgeValuesHashed()) {
            for (unsigned z=0; z<getSize(); z++) {
                MEDDLY_DCASSERT(!parent->isTransparentEdge(edgeval(z), down(z)));
                s.push(index(z), unsigned(down(z)));
                edgeval(z).hash(s);
            }
        } else {
            for (unsigned z=0; z<getSize(); z++) {
                MEDDLY_DCASSERT(down(z)!=parent->getTransparentNode());
                s.push(index(z), unsigned(down(z)));
            }
        }
    } else {
        if (parent->areEdgeValuesHashed()) {
            for (unsigned n=0; n<getSize(); n++) {
                if (!parent->isTransparentEdge(edgeval(n), down(n))) {
                    s.push(n, unsigned(down(n)));
                    edgeval(n).hash(s);
                }
            }
        } else {
            for (unsigned n=0; n<getSize(); n++) {
                if (down(n)!=parent->getTransparentNode()) {
                    s.push(n, unsigned(down(n)));
                }
            }
        }
    }

    the_hash = s.finish();
#ifdef DEVELOPMENT_CODE
    has_hash = true;
#endif
}

#ifdef DEBUG_UNPACKED_HASH
void MEDDLY::unpacked_node::debugHash(output &debug) const
{
    debug << "Hash computation:\n";
    hash_stream s;
    s.start(0);

    if (extra_hashed_size) {
        debug << "    push extra hashed info\n";
        s.push(extra_hashed, extra_hashed_size);
    }

    if (isSparse()) {
        if (parent->areEdgeValuesHashed()) {
            for (unsigned z=0; z<getSize(); z++) {
                MEDDLY_DCASSERT(!parent->isTransparentEdge(edgeval(z), down(z)));
                debug << "    push " << index(z) << ", "
                      << down(z) << ", ";
                edgeval(z).show(debug);
                debug << "\n";

                s.push(index(z), unsigned(down(z)));
                edgeval(z).hash(s);
            }
        } else {
            for (unsigned z=0; z<getSize(); z++) {
                MEDDLY_DCASSERT(down(z)!=parent->getTransparentNode());
                debug << "    push " << index(z) << ", " << down(z) << "\n";

                s.push(index(z), unsigned(down(z)));
            }
        }
    } else {
        if (parent->areEdgeValuesHashed()) {
            for (unsigned n=0; n<getSize(); n++) {
                if (!parent->isTransparentEdge(edgeval(n), down(n))) {

                    debug << "    push " << n << ", "
                        << down(n) << ", ";
                    edgeval(n).show(debug);
                    debug << "\n";

                    s.push(n, unsigned(down(n)));
                    edgeval(n).hash(s);
                }
            }
        } else {
            for (unsigned n=0; n<getSize(); n++) {
                if (down(n)!=parent->getTransparentNode()) {
                    debug << "    push " << index(n) << ", " << down(n) << "\n";
                    s.push(n, unsigned(down(n)));
                }
            }
        }
    }

    debug << "stream finish: " << s.finish() << "\n";
    debug << "stored hash  : " << the_hash << "\n";
}
#endif

// remove all edges starting at the given index
void MEDDLY::unpacked_node::trim()
{
#ifdef ALLOW_EXTENSIBLE
    if (isTrim()) return;

    //
    // If extensible edge is transparent,
    // mark the node as not-extensible and return
    //
    if (parent->getTransparentNode() == _down[getSize()-1]) {
        markAsNotExtensible();
        // TBD: do we need to shrink the size by 1?
        // --size;
        return;
    }

    MEDDLY_DCASSERT(modparent);

    if (isSparse()) {
        unsigned z = getSize()-1;
        while (z && (_index[z-1]+1 == _index[z]) && (_down[z-1] == _down[z])) {
            modparent->unlinkNode(_down[z]);
            z--;
        }
        shrink(z+1);
    } else {
        unsigned z = getSize()-1;
        while (z && (_down[z-1] == _down[z])) {
            modparent->unlinkNode(_down[z]);
            z--;
        }
        shrink(z+1);
    }

    MEDDLY_DCASSERT(isExtensible() && isTrim());
#endif
}


// check is the node is written in order,
// if not rearrange it in ascending order of indices.
void MEDDLY::unpacked_node::sort()
{
    if (!isSparse()) return;
    MEDDLY_DCASSERT(_index);

    //
    // First, scan indexes to see if we're already sorted,
    // and obtain the largest index.
    //
    unsigned maxind = 0;
    bool sorted = true;
    for (unsigned i=0; i<getSize(); i++) {
        unsigned ip1 = 1+index(i);
        if (ip1 > maxind) {
            maxind = ip1;
        } else {
            sorted = false;
        }
    }
    if (sorted) return;

    //
    // Build an array that gives, for each index,
    // its current position +1 in the list.
    //
    unsigned* position = new unsigned[maxind];
    for (unsigned i=0; i<maxind; i++) {
        position[i] = 0;
    }
    for (unsigned i=0; i<getSize(); i++) {
        if (position[index(i)]) {
            // two of the same indexes, that's bad
            throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
        }
        position[index(i)] = i+1;
    }

    //
    // Do the sort.
    // Scan the position array, skip zeroes (indexes not present).
    // Then, swap the current position of the element with its
    // desired position.
    //
    unsigned zd = 0;
    for (unsigned i=0; i<maxind; i++) {
        if (position[i]) {
            const unsigned zn = position[i]-1;
            if (zn != zd) {
                SWAP(position[index(zd)], position[index(zn)]);
                SWAP(_edge[zd], _edge[zn]);
                SWAP(_down[zd], _down[zn]);
                SWAP(_index[zd], _index[zn]);
            }
            ++zd;
        }
    }

    //
    // Cleanup
    //
    delete[] position;
}

bool MEDDLY::unpacked_node::isSorted() const
{
    if (!isSparse()) {
        return true;
    }
    for (unsigned z = 1; z < getSize(); z++) {
        if (index(z-1)>= index(z)) return false;
    }
    return true;
}

void MEDDLY::unpacked_node::clear(unsigned low, unsigned high)
{
    CHECK_RANGE(__FILE__, __LINE__, 0u, low, alloc);
    CHECK_RANGE(__FILE__, __LINE__, 0u, high, alloc+1);
    MEDDLY_DCASSERT(_down);

    if (hasEdges()) {
        MEDDLY_DCASSERT(_edge);
        for (unsigned i=low; i<high; i++) {
            parent->getTransparentEdge(_edge[i], _down[i]);
        }
    } else {
        for (unsigned i=low; i<high; i++) {
            _down[i] = parent->getTransparentNode();
        }
    }
#ifdef DEVELOPMENT_CODE
    has_hash = false;
#endif
}

//
// Static methods for free lists
//

MEDDLY::unpacked_node* MEDDLY::unpacked_node::New(const forest* f,
        node_storage_flags ns)
{
    MEDDLY_DCASSERT(f);
    const unsigned FID = f->FID();
    CHECK_RANGE(__FILE__, __LINE__, 1u, FID, ForListsAlloc);
    MEDDLY_DCASSERT(ForLists);

    unpacked_node* n = ForLists[FID].popFree(ns);
    if (n) {
        //
        // We have a recycled node.
        //
        MEDDLY_DCASSERT(n->isAttachedTo(f));
        n->prev = nullptr;
        n->next = nullptr;
    } else {
        n = new unpacked_node(f, ns);
    }
    MEDDLY_DCASSERT(n->isAttachedTo(f));
    n->orig_was_identity = false;
#ifdef DEVELOPMENT_CODE
    n->has_hash = false;
    n->can_be_recycled = true;
#endif
    n->is_full = (FULL_ONLY == ns);
    return n;
}

void MEDDLY::unpacked_node::AddToBuildList(unpacked_node* n)
{
    if (!n) return;

    CHECK_RANGE(__FILE__, __LINE__, 1u, n->pFID, ForListsAlloc);
    MEDDLY_DCASSERT(ForLists);

    if (ForLists[n->pFID].building) {
        ForLists[n->pFID].building->prev = n;
    }
    n->next = ForLists[n->pFID].building;
    n->prev = nullptr;
    ForLists[n->pFID].building = n;

#ifdef DEBUG_FORLISTS
    std::cerr << "Added " << n << " to build list # " << n->pFID
              << "\nList is now: -> ";
    showDoubly(ForLists[n->pFID].building);
#endif
}

void MEDDLY::unpacked_node::MarkWritable(node_marker &M)
{
    MEDDLY_DCASSERT(M.getParent());
    const unsigned FID = M.getParent()->FID();

    for (const unpacked_node* curr = ForLists[FID].building; curr; curr=curr->next)
    {
#ifdef DEBUG_MARK_SWEEP
        std::cerr   << "Traversing unpacked node at level "
                    << curr->getLevel() << "\n\t";
        stream_output s(std::cerr);
        curr->show(s, true);
        std::cerr << '\n';
#endif
        for (unsigned i=0; i<curr->getSize(); i++) {
            M.mark(curr->down(i));
        }
        M.mark(curr->mark_extra);
    } // for curr
}

void MEDDLY::unpacked_node::AddToIncomingCounts(const forest* F,
        std::vector <unsigned> &incount)
{
    MEDDLY_DCASSERT(F);
    for (const unpacked_node* curr = ForLists[F->FID()].building;
            curr; curr=curr->next)
    {
#ifdef DEBUG_MARK_SWEEP
        std::cerr   << "Traversing unpacked node at level "
                    << curr->getLevel() << "\n\t";
        stream_output s(std::cerr);
        curr->show(s, true);
        std::cerr << '\n';
#endif
        for (unsigned i=0; i<curr->getSize(); i++) {
            if (curr->down(i) > 0) {
                ++incount[curr->down(i)];
            }
        }
        if (curr->mark_extra > 0) {
            ++incount[curr->mark_extra];
        }
    } // for curr
}

void MEDDLY::unpacked_node::Recycle(unpacked_node* r)
{
    if (!r) return;

    if (!ForLists) {
        //
        // Lists have all been destroyed; this must be a late recycle.
        // Just delete it.
        //
        delete r;
        return;
    }

#ifdef DEVELOPMENT_CODE
    MEDDLY_DCASSERT(r->can_be_recycled);
#endif
    CHECK_RANGE(__FILE__, __LINE__, 1u, r->pFID, ForListsAlloc);
    MEDDLY_DCASSERT(ForLists);

    if (r->modparent) {
        // remove from building list (it's doubly linked)
        unpacked_node* n = r->next;
        if (r->prev) {
            unpacked_node* p = r->prev;
            p->next = n;
            if (n) n->prev = p;
        } else {
            // r is the front of the list
            MEDDLY_DCASSERT(ForLists[r->pFID].building == r);
            ForLists[r->pFID].building = n;
            if (n) n->prev = nullptr;
        }
        r->prev = nullptr;

        r->modparent = nullptr;

#ifdef DEBUG_FORLISTS
        std::cerr << "Recycle removed " << r << " from buildlist #" << r->pFID
                  << "\nList is now: -> ";
        showDoubly(ForLists[r->pFID].building);
#endif
    }

    //
    // Push onto recycled list
    //
    ForLists[r->pFID].pushFree(r, r->nodestor);
}

void MEDDLY::unpacked_node::initForest(const forest* f)
{
    if (!f) return;
    const unsigned FID = f->FID();
    if (FID >= ForListsAlloc) {
        unsigned newalloc = (FID/16 + 1) * 16;

        ForLists = (unpacked_lists*)
                   realloc(ForLists, newalloc * sizeof (unpacked_lists));

        if (!ForLists) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
        ForListsAlloc = newalloc;
    }

    ForLists[FID].init();
}

void MEDDLY::unpacked_node::doneForest(const forest* f)
{
    if (!f) return;
    const unsigned FID = f->FID();
    CHECK_RANGE(__FILE__, __LINE__, 1u, FID, ForListsAlloc);

    deleteList(ForLists[FID].building);
    deleteList(ForLists[FID].free_full_only);
    deleteList(ForLists[FID].free_full_or_sparse);
    deleteList(ForLists[FID].free_sparse_only);
}

void MEDDLY::unpacked_node::initStatics()
{
    ForLists = nullptr;
    ForListsAlloc = 0;
}

void MEDDLY::unpacked_node::doneStatics()
{
    free(ForLists);
    ForLists = nullptr;
    ForListsAlloc = 0;
}

//
// Private helpers
//

void MEDDLY::unpacked_node::expand(unsigned ns)
{
    if (ns > alloc) {
        unsigned nalloc = ((ns/16)+1)*16;
        MEDDLY_DCASSERT(nalloc > ns);
        MEDDLY_DCASSERT(nalloc>0);
        MEDDLY_DCASSERT(nalloc>alloc);
        _down = (node_handle*) realloc(_down, nalloc*sizeof(node_handle));
        if (!_down) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
        if (FULL_ONLY == nodestor) {
            MEDDLY_DCASSERT(nullptr == _index);
        } else {
            _index = (unsigned*) realloc(_index, nalloc*sizeof(unsigned));
            if (!_index) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
        }
        if (edge_type::VOID == the_edge_type) {
            MEDDLY_DCASSERT(nullptr == _edge);
        } else {
            _edge = (edge_value*) realloc(_edge, nalloc*sizeof(edge_value));
            if (!_edge) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
        }
        alloc = nalloc;
    }
    size = ns;
#ifdef DEVELOPMENT_CODE
    has_hash = false;
#endif
}

void MEDDLY::unpacked_node::showSingly(const unpacked_node* list)
{
    while (list) {
        std::cerr << " -> " << list;
        list = list->next;
    }
    std::cerr << " -|\n";
}

void MEDDLY::unpacked_node::showDoubly(const unpacked_node* list)
{
    while (list) {
        std::cerr << " -> (prev " << list->prev << ") " << list;
        list = list->next;
    }
    std::cerr << " -|\n";
}

#if 0

// ******************************************************************
// *                                                                *
// *                     unreduced_lists struct                     *
// *                                                                *
// ******************************************************************

struct MEDDLY::unreduced_lists {
    unreduced_node* recycled;
    unreduced_node* building;
};

// ******************************************************************
// *                                                                *
// *                 unreduced_node  static members                 *
// *                                                                *
// ******************************************************************

MEDDLY::unreduced_lists*    MEDDLY::unreduced_node::ForLists;
unsigned                    MEDDLY::unreduced_node::ForListsAlloc;

char*                       MEDDLY::unreduced_node::free_headers[16];
MEDDLY::node_handle*        MEDDLY::unreduced_node::free_down[16];
unsigned*                   MEDDLY::unreduced_node::free_index[16];
MEDDLY::edge_value*         MEDDLY::unreduced_node::free_edge[16];


// ******************************************************************
// *                                                                *
// *                                                                *
// *                     unreduced_node methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::unreduced_node::unreduced_node()
{
    parent = nullptr;
    build_list_FID = 0;

    next = nullptr;
    prev = nullptr;

    _header = nullptr;
    _down = nullptr;
    _index = nullptr;
    _edge = nullptr;

    mark_extra = 0;
    size = 0;
    alloc = 0;

    _header_slot = 0;
    unhashed_header_bytes = 0;
    hashed_header_bytes = 0;

#ifdef DEVELOPMENT_CODE
    has_hash = false;
#endif
}

MEDDLY::unreduced_node::~unreduced_node()
{
    freeNode();
}

void MEDDLY::unreduced_node::initFromNode(const forest* f, node_handle node,
                node_storage_flags fs)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));


    // TBD
    orig_was_identity = false;
}

void MEDDLY::unreduced_node::initRedundant(const forest *f, int k,
        const edge_value &ev, node_handle node, node_storage_flags fs)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    MEDDLY_DCASSERT(ev.hasType(f->getEdgeType()));

    allocNode(f,  f->getLevelSize(k), fs);
    level = k;

    if (ev.isVoid()) {
        if (isFull()) {
            for (unsigned i=0; i<getSize(); i++) {
                setFull(i, node);
            }
        } else {
            for (unsigned i=0; i<getSize(); i++) {
                setSparse(i, i, node);
            }
        }
    } else {
        if (isFull()) {
            for (unsigned i=0; i<getSize(); i++) {
                setFull(i, ev, node);
            }
        } else {
            for (unsigned i=0; i<getSize(); i++) {
                setSparse(i, i, ev, node);
            }
        }
    }

    orig_was_identity = false;
}

void MEDDLY::unreduced_node::initIdentity(const forest *f, int k, unsigned i,
        const edge_value &ev, node_handle node, node_storage_flags fs)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    MEDDLY_DCASSERT(ev.hasType(f->getEdgeType()));

    level = k;
    if (FULL_ONLY == fs) {
        allocNode(f, f->getLevelSize(k), fs);
        clear(0, getSize());
        setFull(i, ev, node);
    } else {
        allocNode(f, 1, fs);
        setSparse(0, i, ev, node);
    }

    orig_was_identity = true;
}

void MEDDLY::unreduced_node::initEmpty(forest* f, int k, unsigned size,
        node_storage_flags fs)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);

    level = k;
    allocNode(f, size, fs);
    if (FULL_ONLY == fs) {
        clear(0, getSize());
    }
    mark_extra = 0;
    AddToBuildList(f, this);

    orig_was_identity = false;
}

void MEDDLY::unreduced_node::clear(unsigned low, unsigned high)
{
    CHECK_RANGE(__FILE__, __LINE__, 0u, low, alloc);
    CHECK_RANGE(__FILE__, __LINE__, 0u, high, 1+alloc);
    if (hasEdges()) {
        for (unsigned i=low; i<high; i++) {
            parent->getTransparentEdge(edgeval(i), down(i));
        }
    } else {
        for (unsigned i=low; i<high; i++) {
            down(i) = parent->getTransparentNode();
        }
    }
#ifdef DEVELOPMENT_CODE
    has_hash = false;
#endif
}

void MEDDLY::unreduced_node::freeNode()
{
    if (parent) {
        //
        // Recycle or free _header
        //
        if (_header) {
           pushFreeHeader(_header, _header_slot);
           _header = nullptr;
        }
        //
        // Recycle or free _down
        //
        const unsigned _down_slot = size2slot(alloc);

        if (_down) {
            pushFreeDown(_down, _down_slot);
            _down = nullptr;
        }
        //
        // Recycle or free _index
        //
        if (_index) {
            pushFreeIndex(_index, _down_slot);
            _index = nullptr;
        }
        //
        // Recycle or free _edge
        //
        if (_edge) {
            pushFreeEdge(_edge, _down_slot);
            _edge = nullptr;
        }
        parent = nullptr;

        if (build_list_FID) {
            RemoveFromBuildList(this);
        }

        alloc = 0;
    }
    MEDDLY_DCASSERT(!_header);
    MEDDLY_DCASSERT(!_down);
    MEDDLY_DCASSERT(!_index);
    MEDDLY_DCASSERT(!_edge);
    MEDDLY_DCASSERT(!build_list_FID);
    MEDDLY_DCASSERT(!next);
    MEDDLY_DCASSERT(!prev);
}

void MEDDLY::unreduced_node::allocNode(const forest* f, unsigned _size,
        node_storage_flags fs)
{
    MEDDLY_DCASSERT(f);

    if (parent != f) {
        if (parent) {
            freeNode();
        }
        MEDDLY_DCASSERT(!parent);
        MEDDLY_DCASSERT(!_header);

        parent = f;
        //
        // Allocate header data if needed
        //
        CHECK_RANGE(__FILE__, __LINE__, 0u, f->unhashedHeaderBytes(), 256u);
        CHECK_RANGE(__FILE__, __LINE__, 0u, f->hashedHeaderBytes(), 256u);

        unhashed_header_bytes = f->unhashedHeaderBytes();
        hashed_header_bytes   = f->hashedHeaderBytes();

        const unsigned hbytes = unhashed_header_bytes + hashed_header_bytes;
        _header_slot = size2slot(hbytes);
        if (hbytes) {
            _header = popFreeHeader(_header_slot);
            if (!_header) {
                const unsigned halloc = slot2size(_header_slot);
                _header = new char[halloc];
            }
        }
    } else {
        if (build_list_FID) {
            RemoveFromBuildList(this);
        }
    }

    expand(_size, (fs != FULL_ONLY), f->getEdgeType() != edge_type::VOID);
    size = _size;
}

void MEDDLY::unreduced_node::expand(unsigned ns, bool make_index,
        bool make_edge)
{
    if (ns > 123456789) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    const unsigned oldslot = size2slot(alloc);
    unsigned newalloc;
    const unsigned newslot = size2slot(ns);
    if (newslot < 15) {
        newalloc = slot2size(newslot);
    } else {
        newalloc = size2slot(15);
        while (newalloc < ns) {
            newalloc += newalloc / 2;
        }
    }
    MEDDLY_DCASSERT(newalloc >= ns);

    _down = (node_handle*) realloc(_down, newalloc*sizeof(node_handle));
    if (!_down) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }

    if (make_index) {
        _index = (unsigned*) realloc(_index, newalloc*sizeof(unsigned));
        if (!_index) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
    } else {
        if (_index) {
            pushFreeIndex(_index, oldslot);
            _index = nullptr;
        }
    }

    if (make_edge) {
        _edge = (edge_value*) realloc(_edge, newalloc*sizeof(edge_value));
        if (!_edge) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
    } else {
        if (_edge) {
            pushFreeEdge(_edge, oldslot);
            _edge = nullptr;
        }
    }
    alloc = newalloc;
#ifdef DEVELOPMENT_CODE
    has_hash = false;
#endif
}

// ******************************************************************
// *                                                                *
// *                 static  unreduced_node methods                 *
// *                                                                *
// ******************************************************************

void MEDDLY::unreduced_node::AddToBuildList(forest* F, unreduced_node* n)
{
    MEDDLY_DCASSERT(F);
    MEDDLY_DCASSERT(n);
    MEDDLY_DCASSERT(F->FID());
    n->build_list_FID = F->FID();

    CHECK_RANGE(__FILE__, __LINE__, 1u, F->FID(), ForListsAlloc);
    MEDDLY_DCASSERT(ForLists);

    if (ForLists[F->FID()].building) {
        ForLists[F->FID()].building->prev = n;
    }
    n->next = ForLists[F->FID()].building;
    n->prev = nullptr;
    ForLists[F->FID()].building = n;

#ifdef DEBUG_FORLISTS
    std::cerr << "Added " << n << " to build list # " << F->FID()
              << "\nList is now: -> ";
    showDoubly(ForLists[F->FID()].building);
#endif
}

void MEDDLY::unreduced_node::RemoveFromBuildList(unreduced_node* n)
{
    MEDDLY_DCASSERT(n);
    if (!ForLists) {
        //
        // Lists have all been destroyed; this must be a late recycle.
        //

        n->next = nullptr;
        n->prev = nullptr;
        n->build_list_FID = 0;
        return;
    }

    //
    // Remove n from its (doubly-linked) building list
    //

    CHECK_RANGE(__FILE__, __LINE__, 1u, n->build_list_FID, ForListsAlloc);
    MEDDLY_DCASSERT(ForLists);

    unreduced_node* next = n->next;
    unreduced_node* prev = n->prev;
    if (prev) {
        prev->next = next;
    } else {
        // we're at the front of the list
        MEDDLY_DCASSERT(ForLists[n->build_list_FID].building == n);
        ForLists[n->build_list_FID].building = next;
    }
    if (next) {
        next->prev = prev;
    }

#ifdef DEBUG_FORLISTS
    std::cerr << "Removed " << n << " from buildlist #" << n->build_list_FID
              << "\nList is now: -> ";
    showDoubly(ForLists[n->build_list_FID].building);
#endif

    n->next = nullptr;
    n->prev = nullptr;
    n->build_list_FID = 0;
}

void MEDDLY::unreduced_node::initForest(const forest* f)
{
    if (!f) return;
    const unsigned FID = f->FID();
    if (FID >= ForListsAlloc) {
        unsigned newalloc = (FID/16 + 1) * 16;

        ForLists = (unreduced_lists*)
                   realloc(ForLists, newalloc * sizeof (unreduced_lists));

        if (!ForLists) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
        ForListsAlloc = newalloc;
    }

    ForLists[FID].recycled = nullptr;
    ForLists[FID].building = nullptr;

#ifdef DEBUG_FORLISTS
    std::cerr << "initForest #" << FID << " lists are " <<
        ForLists[FID].building << " and " << ForLists[FID].recycled << "\n";
#endif
}

void MEDDLY::unreduced_node::doneForest(const forest* f)
{
    if (!f) return;
    const unsigned FID = f->FID();
    CHECK_RANGE(__FILE__, __LINE__, 1u, FID, ForListsAlloc);

    // Delete the building list
    while (ForLists[FID].building) {
        unreduced_node* p = ForLists[FID].building;
        ForLists[FID].building = p->next;
        delete p;
    }
    // Delete the recycled list
    while (ForLists[FID].recycled) {
        unreduced_node* p = ForLists[FID].recycled;
        ForLists[FID].recycled = p->next;
        delete p;
    }
#ifdef DEBUG_FORLISTS
    std::cerr << "doneForest #" << FID << " lists are " <<
        ForLists[FID].building << " and " << ForLists[FID].recycled << "\n";
#endif
}

void MEDDLY::unreduced_node::showDoubly(const unreduced_node* list)
{
    while (list) {
        std::cerr << " -> (prev " << list->prev << ") " << list;
        list = list->next;
    }
    std::cerr << " -|\n";
}


void MEDDLY::unreduced_node::initStatics()
{
    ForLists = nullptr;
    ForListsAlloc = 0;

    for (unsigned i=0; i<16; i++) {
        free_headers[i] = nullptr;
        free_down[i]    = nullptr;
        free_index[i]   = nullptr;
        free_edge[i]    = nullptr;
    }
}

void MEDDLY::unreduced_node::doneStatics()
{
    free(ForLists);
    ForLists = nullptr;
    ForListsAlloc = 0;

    for (unsigned i=0; i<16; i++) {
        for (;;) {
            char* hdr = popFreeHeader(i);
            if (!hdr) break;
            delete[] hdr;   // won't be resized
        }
        for (;;) {
            node_handle* dn = popFreeDown(i);
            if (!dn) break;
            free(dn);
        }
        for (;;) {
            unsigned* ix = popFreeIndex(i);
            if (!ix) break;
            free(ix);
        }
        for (;;) {
            edge_value* ev = popFreeEdge(i);
            if (!ev) break;
            free(ev);
        }
    }
}

#endif // turn off unreduced_node

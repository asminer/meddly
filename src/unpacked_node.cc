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
    unpacked_node* recycled;
    unpacked_node* building;
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

MEDDLY::unpacked_node::unpacked_node(const forest* f)
{
    next = nullptr;
    prev = nullptr;

    parent = nullptr;
    modparent = nullptr;
    pFID = 0;

#ifndef USE_VECTOR
    _down = nullptr;
    _index = nullptr;
    _edge = nullptr;

    alloc = 0;
    size = 0;
#endif

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

MEDDLY::unpacked_node::~unpacked_node()
{
#ifndef USE_VECTOR
    delete[] _down;
    delete[] _index;
    delete[] _edge;
#endif

    free(extra_unhashed);
    free(extra_hashed);
}

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
}


void MEDDLY::unpacked_node::initRedundant(const forest *f, int k,
    node_handle node, node_storage_flags fs)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(isAttachedTo(f));
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
#ifdef ALLOW_EXTENSIBLE
    is_extensible = f->isExtensibleLevel(k);
    resize( is_extensible ? 1 : unsigned(f->getLevelSize(k)) );
#else
    resize( unsigned(f->getLevelSize(k)) );
#endif
    level = k;

    if (SPARSE_ONLY == fs) {
        setSparse();
        for (unsigned i=0; i<getSize(); i++) {
            setSparse(i, i, node);
        }
    } else {
        setFull();
        for (unsigned i=0; i<getSize(); i++) {
            setFull(i, node);
        }
        is_full = true;
    }
}

void MEDDLY::unpacked_node::initRedundant(const forest *f, int k,
    const edge_value &ev, node_handle node, node_storage_flags fs)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(isAttachedTo(f));
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
#ifdef ALLOW_EXTENSIBLE
    is_extensible = f->isExtensibleLevel(k);
    resize( is_extensible ? 1 : unsigned(f->getLevelSize(k)) );
#else
    resize( unsigned(f->getLevelSize(k)) );
#endif
    level = k;

    if (SPARSE_ONLY == fs) {
        setSparse();
        for (unsigned i=0; i<getSize(); i++) {
            setSparse(i, i, ev, node);
        }
    } else {
        setFull();
        for (unsigned i=0; i<getSize(); i++) {
            setFull(i, ev, node);
        }
    }
}



void MEDDLY::unpacked_node::initIdentity(const forest *f, int k,
  unsigned i, node_handle node, node_storage_flags fs)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(isAttachedTo(f));
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    level = k;

    if (FULL_ONLY == fs) {
        setFull();
        resize(f->getLevelSize(k));
        clear(0, getSize());

        setFull(i, node);
    } else {
        setSparse();
        resize(1);

        setSparse(0, i, node);
    }
}

void MEDDLY::unpacked_node::initIdentity(const forest *f, int k, unsigned i,
        const edge_value &ev, node_handle node, node_storage_flags fs)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(isAttachedTo(f));
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    level = k;

    if (FULL_ONLY == fs) {
        setFull();
        resize(f->getLevelSize(k));
        clear(0, getSize());

        setFull(i, ev, node);
    } else {
        setSparse();
        resize(1);

        setSparse(0, i, ev, node);
    }
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
            s << (unsigned long) _index[z] << ":";
        } else {
            if (z) s.put('|');
        }
        parent->showEdge(s, _edge[z], _down[z]);
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
            s.put((unsigned long) _index[z]);
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
            _edge[z].write(s);
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
                MEDDLY_DCASSERT(!parent->isTransparentEdge(_down[z], _edge[z]));
                s.push(_index[z], unsigned(_down[z]));
                _edge[z].hash(s);
            }
        } else {
            for (unsigned z=0; z<getSize(); z++) {
                MEDDLY_DCASSERT(_down[z]!=parent->getTransparentNode());
                s.push(_index[z], unsigned(_down[z]));
            }
        }
    } else {
        if (parent->areEdgeValuesHashed()) {
            for (unsigned n=0; n<getSize(); n++) {
                if (!parent->isTransparentEdge(_down[n], _edge[n])) {
                    s.push(n, unsigned(_down[n]));
                    _edge[n].hash(s);
                }
            }
        } else {
            for (unsigned n=0; n<getSize(); n++) {
                if (_down[n]!=parent->getTransparentNode()) {
                    s.push(n, unsigned(_down[n]));
                }
            }
        }
    }

    the_hash = s.finish();
#ifdef DEVELOPMENT_CODE
    has_hash = true;
#endif
}


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
#ifndef USE_VECTOR
    MEDDLY_DCASSERT(_index);
#endif

    //
    // First, scan indexes to see if we're already sorted,
    // and obtain the largest index.
    //
    unsigned maxind = 0;
    bool sorted = true;
    for (unsigned i=0; i<getSize(); i++) {
        unsigned ip1 = 1+_index[i];
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
        if (position[_index[i]]) {
            // two of the same indexes, that's bad
            throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
        }
        position[_index[i]] = i+1;
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
                SWAP(_edge[zd], _edge[zn]);
                SWAP(_down[zd], _down[zn]);
                SWAP(position[_index[zd]], position[_index[zn]]);
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
        if (_index[z-1] >= _index[z]) return false;
    }
    return true;
}

void MEDDLY::unpacked_node::clear(unsigned low, unsigned high)
{
#ifndef USE_VECTOR
    CHECK_RANGE(__FILE__, __LINE__, 0u, low, alloc);
    CHECK_RANGE(__FILE__, __LINE__, 0u, high, alloc+1);
    MEDDLY_DCASSERT(_down);
#endif

    if (hasEdges()) {
#ifndef USE_VECTOR
        MEDDLY_DCASSERT(_edge);
#endif
        for (unsigned i=low; i<high; i++) {
            parent->getTransparentEdge(_down[i], _edge[i]);
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

MEDDLY::unpacked_node* MEDDLY::unpacked_node::New(const forest* f)
{
    MEDDLY_DCASSERT(f);
    const unsigned FID = f->FID();
    CHECK_RANGE(__FILE__, __LINE__, 1u, FID, ForListsAlloc);
    MEDDLY_DCASSERT(ForLists);

    if (ForLists[FID].recycled) {
        // pull from free list
        unpacked_node* n = ForLists[FID].recycled;
        ForLists[FID].recycled = n->next;
        n->prev = nullptr;
        n->next = nullptr;
#ifdef DEVELOPMENT_CODE
        n->has_hash = false;
#endif
        MEDDLY_DCASSERT(n->isAttachedTo(f));

#ifdef DEBUG_FORLISTS
        std::cerr << "New removed " << n << " from recycled list #" << FID
                  << "\nList is now: -> ";
        showSingly(ForLists[FID].recycled);
#endif

        return n;
    } else {
        // Our free list is empty
        unpacked_node* n = new unpacked_node(f);
#ifdef DEVELOPMENT_CODE
        n->can_be_recycled = true;
#endif
        MEDDLY_DCASSERT(n->isAttachedTo(f));
        return n;
    }
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
            M.mark(curr->_down[i]);
        }
    } // for curr
}

void MEDDLY::unpacked_node::Recycle(unpacked_node* r)
{
    if (!r) return;

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
    // Push onto recycled list (it's singly linked)
    //
    r->next = ForLists[r->pFID].recycled;
    ForLists[r->pFID].recycled = r;

#ifdef DEBUG_FORLISTS
    std::cerr << "Recycle added " << r << " to freelist #" << r->pFID
              << "\nList is now: -> ";
    showSingly(ForLists[r->pFID].recycled);
#endif
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

    ForLists[FID].recycled = nullptr;
    ForLists[FID].building = nullptr;

#ifdef DEBUG_FORLISTS
    std::cerr << "initForest #" << FID << " lists are " <<
        ForLists[FID].building << " and " << ForLists[FID].recycled << "\n";
#endif
}

void MEDDLY::unpacked_node::doneForest(const forest* f)
{
    if (!f) return;
    const unsigned FID = f->FID();
    CHECK_RANGE(__FILE__, __LINE__, 1u, FID, ForListsAlloc);

    deleteList(ForLists[FID].building);
    deleteList(ForLists[FID].recycled);
#ifdef DEBUG_FORLISTS
    std::cerr << "doneForest #" << FID << " lists are " <<
        ForLists[FID].building << " and " << ForLists[FID].recycled << "\n";
#endif
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

#ifndef USE_VECTOR
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
        _index = (unsigned*) realloc(_index, nalloc*sizeof(unsigned));
        if (!_index) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
        _edge = (edge_value*) realloc(_edge, nalloc*sizeof(edge_value));
        if (!_edge) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
        alloc = nalloc;
    }
    size = ns;
#ifdef DEVELOPMENT_CODE
    has_hash = false;
#endif
}
#endif

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


// ===================================================================
//
// Deprecated as of version 0.17.4
//
// ===================================================================

#ifdef ALLOW_DEPRECATED_0_17_4

/*
  Initializers

  Extensible nodes
        + Every node at level k, where level k represents an extensible
          variable, is represented by an extensible node.
        + Whether a node is extensible or not is determined by querying
          the corresponding level's property.
        + The last downpointer in an extensible node is considered to
          repeat for all indices till +infinity.
*/


void MEDDLY::unpacked_node::initRedundant(const forest *f, int k,
  float ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(isAttachedTo(f));
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    is_extensible = f->isExtensibleLevel(k);
    resize( is_extensible ? 1 : unsigned(f->getLevelSize(k)) );
    level = k;

    for (unsigned i=0; i<getSize(); i++) {
        _down[i] = node;
        _edge[i].set(ev);
    }
    if (!full) {
        for (unsigned i=0; i<getSize(); i++) _index[i] = i;
    }
    is_full = full;
}

void MEDDLY::unpacked_node::initRedundant(const forest *f, int k,
  int ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(isAttachedTo(f));
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    is_extensible = f->isExtensibleLevel(k);
    resize( is_extensible ? 1 : unsigned(f->getLevelSize(k)) );
    level = k;

    for (unsigned i=0; i<getSize(); i++) {
        _down[i] = node;
        _edge[i].set(ev);
    }
    if (!full) {
        for (unsigned i=0; i<getSize(); i++) _index[i] = i;
    }
    is_full = full;
}

void MEDDLY::unpacked_node::initRedundant(const forest *f, int k,
  long ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(isAttachedTo(f));
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    is_extensible = f->isExtensibleLevel(k);
    resize( is_extensible ? 1 : unsigned(f->getLevelSize(k)) );
    level = k;

    for (unsigned i=0; i<getSize(); i++) {
        _down[i] = node;
        _edge[i].set(ev);
    }
    if (!full) {
        for (unsigned i=0; i<getSize(); i++) _index[i] = i;
    }
    is_full = full;
}


void MEDDLY::unpacked_node::initIdentity(const forest *f, int k,
  unsigned i, int ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(isAttachedTo(f));
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    level = k;

    if (full) {
        setFull();
        resize(f->getLevelSize(k));
        clear(0, getSize());

        _down[i] = node;
        _edge[i].set(ev);
    } else {
        setSparse();
        resize(1);

        _down[0] = node;
        _index[0] = i;
        _edge[0].set(ev);
    }
}

void MEDDLY::unpacked_node::initIdentity(const forest *f, int k,
  unsigned i, long ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(isAttachedTo(f));
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    level = k;

    if (full) {
        setFull();
        resize(f->getLevelSize(k));
        clear(0, getSize());

        _down[i] = node;
        _edge[i].set(ev);
    } else {
        setSparse();
        resize(1);

        _down[0] = node;
        _index[0] = i;
        _edge[0].set(ev);
    }
}

void MEDDLY::unpacked_node::initIdentity(const forest *f, int k,
  unsigned i, float ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(isAttachedTo(f));
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    level = k;

    if (full) {
        setFull();
        resize(f->getLevelSize(k));
        clear(0, getSize());

        _down[i] = node;
        _edge[i].set(ev);
    } else {
        setSparse();
        resize(1);

        _down[0] = node;
        _index[0] = i;
        _edge[0].set(ev);
    }
}

#endif

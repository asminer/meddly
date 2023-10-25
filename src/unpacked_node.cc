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

#include <map>
#include "defines.h"
#include "unpacked_node.h"
#include "forest.h"
#include "hash_stream.h"
#include "terminal.h"

// ******************************************************************
// *                     unpacked_node  statics                     *
// ******************************************************************

MEDDLY::unpacked_node* MEDDLY::unpacked_node::freeList;
MEDDLY::unpacked_node* MEDDLY::unpacked_node::buildList;


// ******************************************************************
// *                                                                *
// *                                                                *
// *                     unpacked_node  methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::unpacked_node::unpacked_node()
{
    parent = nullptr;
    eparent = nullptr;
    modparent = nullptr;
    extra_unhashed = 0;
    ext_uh_alloc = 0;
    ext_uh_size = 0;
    extra_hashed = 0;
    ext_h_alloc = 0;
    ext_h_size = 0;
    down = nullptr;
    index = nullptr;
    edge = nullptr;
    is_extensible = false;
    can_be_extensible = false;
    alloc = 0;
    // ealloc = 0;
    size = 0;
    nnzs = 0;
    level = 0;
    is_in_build_list = false;
#ifdef DEVELOPMENT_CODE
    has_hash = false;
#endif
}

MEDDLY::unpacked_node::~unpacked_node()
{
    clear();
}

void MEDDLY::unpacked_node::clear()
{
    free(extra_unhashed);
    free(extra_hashed);
    free(down);
    free(index);
    free(edge);
    down = nullptr;
    index = nullptr;
    edge = nullptr;
    is_extensible = false;
    alloc = 0;
    // ealloc = 0;
    size = 0;
    nnzs = 0;
    level = 0;
}

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

void MEDDLY::unpacked_node::initRedundant(const expert_forest *f, int k,
  node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    unsigned nsize = f->isExtensibleLevel(k) ? 1 : unsigned(f->getLevelSize(k));
    bind_to_forest(f, k, nsize, full);
    for (unsigned i=0; i<nsize; i++) {
        down[i] = node;
    }
    if (!full) {
        for (unsigned i=0; i<nsize; i++) index[i] = i;
        nnzs = nsize;
    }
    is_extensible = f->isExtensibleLevel(k);
}

void MEDDLY::unpacked_node::initRedundant(const expert_forest *f, int k,
  float ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    unsigned nsize = unsigned(f->getLevelSize(k));
    bind_to_forest(f, k, nsize, full);
    for (unsigned i=0; i<nsize; i++) {
        down[i] = node;
        edge[i].set(ev);
    }
    if (!full) {
        for (unsigned i=0; i<nsize; i++) index[i] = i;
        nnzs = nsize;
    }
}

void MEDDLY::unpacked_node::initRedundant(const expert_forest *f, int k,
  int ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    unsigned nsize = unsigned(f->getLevelSize(k));
    bind_to_forest(f, k, nsize, full);
    for (unsigned i=0; i<nsize; i++) {
        down[i] = node;
        edge[i].set(ev);
    }
    if (!full) {
        for (unsigned i=0; i<nsize; i++) index[i] = i;
        nnzs = nsize;
    }
}

void MEDDLY::unpacked_node::initRedundant(const expert_forest *f, int k,
  long ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    unsigned nsize = unsigned(f->getLevelSize(k));
    bind_to_forest(f, k, nsize, full);
    for (unsigned i=0; i<nsize; i++) {
        down[i] = node;
        edge[i].set(ev);
    }
    if (!full) {
        for (unsigned i=0; i<nsize; i++) index[i] = i;
        nnzs = nsize;
    }
}

void MEDDLY::unpacked_node::initIdentity(const expert_forest *f, int k,
  unsigned i, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    unsigned nsize = unsigned(f->getLevelSize(k));
    if (full) {
        bind_to_forest(f, k, nsize, full);
        clearFullEdges();
        down[i] = node;
    } else {
        bind_to_forest(f, k, 1, full);
        nnzs = 1;
        down[0] = node;
        index[0] = i;
    }
}

void MEDDLY::unpacked_node::initIdentity(const expert_forest *f, int k,
  unsigned i, int ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    unsigned nsize = unsigned(f->getLevelSize(k));
    if (full) {
        bind_to_forest(f, k, nsize, full);
        clearFullEdges();
        down[i] = node;
        edge[i].set(ev);
    } else {
        bind_to_forest(f, k, 1, full);
        nnzs = 1;
        down[0] = node;
        edge[0].set(ev);
        index[0] = i;
    }
}

void MEDDLY::unpacked_node::initIdentity(const expert_forest *f, int k,
  unsigned i, long ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    unsigned nsize = unsigned(f->getLevelSize(k));
    if (full) {
        bind_to_forest(f, k, nsize, full);
        clearFullEdges();
        down[i] = node;
        edge[i].set(ev);
    } else {
        bind_to_forest(f, k, 1, full);
        nnzs = 1;
        down[0] = node;
        edge[0].set(ev);
        index[0] = i;
    }
}

void MEDDLY::unpacked_node::initIdentity(const expert_forest *f, int k,
  unsigned i, float ev, node_handle node, bool full)
{
    MEDDLY_DCASSERT(f);
    MEDDLY_DCASSERT(k);
    MEDDLY_DCASSERT(f->isTerminalNode(node) || !f->isDeletedNode(node));
    unsigned nsize = unsigned(f->getLevelSize(k));
    if (full) {
        bind_to_forest(f, k, nsize, full);
        clearFullEdges();
        down[i] = node;
        edge[i].set(ev);
    } else {
        bind_to_forest(f, k, 1, full);
        nnzs = 1;
        down[0] = node;
        edge[0].set(ev);
        index[0] = i;
    }
}

/*
  Usage
*/

void MEDDLY::unpacked_node::show(output &s, bool details) const
{
    unsigned stop;
    if (isSparse()) {
        if (details) s << "nnzs: " << (unsigned long)nnzs << (isExtensible()? "*": "") << " ";
        s << "down: (";
        stop = nnzs;
    } else {
        if (details) s << "size: " << (unsigned long)size << (isExtensible()? "*": "") << " ";
        s << "down: [";
        stop = size;
    }

    for (unsigned z=0; z<stop; z++) {
        if (isSparse()) {
            if (z) s << ", ";
            s << (unsigned long) i(z) << ":";
        } else {
            if (z) s.put('|');
        }
        eparent->showEdge(s, edge[z], down[z]);
    }

    if (isExtensible()) s.put('*');

    if (isSparse()) {
        s.put(')');
    } else {
        s.put(']');
    }

    // show extra header stuff
    eparent->showHeaderInfo(s, *this);
}

void MEDDLY::unpacked_node::write(output &s, const node_handle* map) const
{
    unsigned stop;
    if (isSparse()) {
        s.put(-long(nnzs));
        stop = nnzs;
    } else {
        s.put(long(size));
        stop = size;
    }

    //
    // write indexes (sparse only)
    //
    if (isSparse()) {
        s.put('\n');
        s.put('\t');
        for (unsigned z=0; z<nnzs; z++) {
            s.put(' ');
            s.put((unsigned long) i(z));
        }
    }

    //
    // write down pointers
    //
    s.put('\n');
    s.put('\t');
    for (unsigned z=0; z<stop; z++) {
        s.put(' ');
        if (d(z) <= 0) {
            // terminal
            terminal t;
            t.setFromHandle(the_terminal_type, d(z));
            t.write(s);
        } else {
            // non-terminal
            s.put("n ");
            s.put(long( map ? map[d(z)] : d(z) ));
        }
    }

    //
    // write edge values, if any
    //
    if (edge_type::VOID != the_edge_type) {
        s.put('\n');
        s.put('\t');
        for (unsigned z=0; z<stop; z++) {
            s.put(' ');
            edge[z].write(s);
        }
    }
    s.put('\n');


    // write extra header stuff, should no-op if there isn't any
    eparent->writeHeaderInfo(s, *this);
}

void MEDDLY::unpacked_node::read(input &s, const node_handle* map)
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
    unsigned stop;
    if (isSparse()) {
#ifdef DEBUG_READ_DD
        std::cerr << "    reading indexes\n";
#endif
        for (unsigned z=0; z<nnzs; z++) {
            s.stripWS();
            long ndx = s.get_integer();
            if (ndx < 0) {
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            i_ref(z) = unsigned(ndx);
        }
        stop = nnzs;
    } else {
        stop = size;
    }

    //
    // read down pointers
    //
#ifdef DEBUG_READ_DD
    std::cerr << "    reading " << stop << " down pointers\n";
#endif
    for (unsigned z=0; z<stop; z++) {
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
            d_ref(z) = modparent->linkNode(map ? map[d] : d);
        } else {
            // terminal
            s.unget(c);
            terminal t;
            t.read(s);
            d_ref(z) = t.getHandle();
        }
    }

    //
    // read edge values, if any
    //
    if (edge_type::VOID != the_edge_type) {
#ifdef DEBUG_READ_DD
        std::cerr << "    reading " << stop << " edge values\n";
#endif
        for (unsigned z=0; z<stop; z++) {
            s.stripWS();
            edge[z].read(s);
        }
    }

    //
    // read extra header stuff, if any
    //
    eparent->readHeaderInfo(s, *this);

#ifdef DEBUG_READ_DD
    std::cerr << "  done unpacked_node::read\n";
#endif
}

void MEDDLY::unpacked_node
::resize(unsigned ns)
{
    size = ns;
    nnzs = ns;
    if (size > alloc) {
        unsigned nalloc = ((ns/8)+1)*8;
        MEDDLY_DCASSERT(nalloc > ns);
        MEDDLY_DCASSERT(nalloc>0);
        MEDDLY_DCASSERT(nalloc>alloc);
        down = (node_handle*) realloc(down, nalloc*sizeof(node_handle));
        if (!down) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        index = (unsigned*) realloc(index, nalloc*sizeof(unsigned));
        if (!index) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        alloc = nalloc;
        edge = (edge_value*) realloc(edge, nalloc*sizeof(edge_value));
        if (!edge) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
}

void MEDDLY::unpacked_node::bind_to_forest(const expert_forest* f,
    int k, unsigned ns, bool full)
{
    parent = f;
    eparent = f;
    modparent = 0;
    level = k;
    is_full = full;

    can_be_extensible = f->isExtensibleLevel(k);
    markAsNotExtensible();
    the_edge_type = f->getEdgeType();
    the_terminal_type = f->getTerminalType();
    resize(ns);

    // Allocate headers
    ext_h_size = f->hashedHeaderBytes();
    if (ext_h_size > ext_h_alloc) {
        ext_h_alloc = ((ext_h_size/8)+1)*8;
        MEDDLY_DCASSERT(ext_h_alloc > ext_h_size);
        MEDDLY_DCASSERT(ext_h_alloc>0);
        extra_hashed =  realloc(extra_hashed, ext_h_alloc);
        if (!extra_hashed) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }

    ext_uh_size = f->unhashedHeaderBytes();
    if (ext_uh_size > ext_uh_alloc) {
        ext_uh_alloc = ((ext_uh_size/8)+1)*8;
        MEDDLY_DCASSERT(ext_uh_alloc > ext_uh_size);
        MEDDLY_DCASSERT(ext_uh_alloc>0);
        extra_unhashed =  realloc(extra_unhashed, ext_uh_alloc);
        if (!extra_unhashed) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
}



void MEDDLY::unpacked_node::removeFromBuildList(unpacked_node* b)
{
    MEDDLY_DCASSERT(b);
    MEDDLY_DCASSERT(b->is_in_build_list);
    MEDDLY_DCASSERT(buildList);
    if (b == buildList) {
#ifdef DEBUG_BUILDLIST
        printf("Removing unpacked node (level %d) from front of build list\n", b->getLevel());
#endif
        // should always happen if we're recursively building
        buildList = b->next;
        b->next = 0;
        return;
    }
#ifdef DEBUG_BUILDLIST
    printf("Removing unpacked node (level %d) from middle of build list\n", b->getLevel());
#endif
    unpacked_node* prev = buildList;
    for (unpacked_node* curr = buildList->next; curr; curr=curr->next) {
        if (b == curr) {
            prev->next = b->next;
            b->next = 0;
            return;
        }
        prev = curr;
    }
    MEDDLY_DCASSERT(0);
}

void MEDDLY::unpacked_node::markBuildListChildren(expert_forest* F)
{
    for (unpacked_node* curr = buildList; curr; curr=curr->next) {
        if (curr->parent != F) {
#ifdef DEBUG_MARK_SWEEP
            printf("Skipping unpacked node for other forest\n");
#endif
            continue;
        }
#ifdef DEBUG_MARK_SWEEP
        printf("Traversing unpacked node at level %d\n\t", curr->getLevel());
        FILE_output s(stdout);
        curr->show(s, true);
        printf("\n");
#endif
        if (curr->isSparse()) {
            // Sparse
            for (unsigned i=0; i<curr->getNNZs(); i++) {
                F->markNode(curr->d(i));
            }
        } else {
            // Full
            for (unsigned i=0; i<curr->getSize(); i++) {
                F->markNode(curr->d(i));
            }
        }
        if (curr->isExtensible()) {
            F->markNode(curr->ext_d());
        }
    } // for curr
}

void MEDDLY::unpacked_node::computeHash()
{
    MEDDLY_DCASSERT(!has_hash);
    trim();

    hash_stream s;
    s.start(0);

    if (ext_h_size) {
        s.push(extra_hashed, ext_h_size);
    }

    if (isSparse()) {
        if (eparent->areEdgeValuesHashed()) {
            for (unsigned z=0; z<nnzs; z++) {
                MEDDLY_DCASSERT(!eparent->isTransparentEdge(d(z), edge[z]));
                s.push(i(z), unsigned(d(z)));
                edge[z].hash(s);
            }
        } else {
            for (unsigned z=0; z<nnzs; z++) {
                MEDDLY_DCASSERT(d(z)!=eparent->getTransparentNode());
                s.push(i(z), unsigned(d(z)));
            }
        }
    } else {
        if (eparent->areEdgeValuesHashed()) {
            for (unsigned n=0; n<size; n++) {
                if (!eparent->isTransparentEdge(d(n), edge[n])) {
                    s.push(n, unsigned(d(n)));
                    edge[n].hash(s);
                }
            }
        } else {
            for (unsigned n=0; n<size; n++) {
                if (d(n)!=eparent->getTransparentNode()) {
                    s.push(n, unsigned(d(n)));
                }
            }
        }
    }

    h = s.finish();
#ifdef DEVELOPMENT_CODE
    has_hash = true;
#endif
}


// check is the node is written in order,
// if not rearrange it in ascending order of indices.
void MEDDLY::unpacked_node::sort()
{
    if (!isSparse()) return;

    unsigned k = 1;
    for (k = 1; k < getNNZs() && i(k-1) < i(k) ; k++);
    if (k == getNNZs()) return; // already sorted

    // sort from (k-1) to (nnz-1)
    --k;
    std::map<unsigned, unsigned> sorter;
    for (unsigned m = k; m < getNNZs(); m++) {
        sorter[i(m)] = m;
    }

    // allocate new arrays for index, node handles and edge-values
    node_handle* old_down = down;
    unsigned* old_index = index;
    edge_value* old_edge = edge;
    unsigned old_nnzs = nnzs;

    down = 0;
    index = 0;
    edge = 0;
    size = 0;
    nnzs = 0;
    alloc = 0;
    resize(old_nnzs);

    // copy into new arrays
    memcpy(down, old_down, sizeof(node_handle) * k);
    memcpy(index, old_index, sizeof(unsigned) * k);
    memcpy(edge, old_edge, sizeof(edge_value) * k);

    for (auto s_iter = sorter.begin(); s_iter != sorter.end(); s_iter++, k++) {
        unsigned old_location = s_iter->second;
        index[k] = old_index[old_location];
        down[k] = old_down[old_location];
        edge[k] = old_edge[old_location];
    }

    free(old_down);
    free(old_index);
    free(old_edge);
}

// remove all edges starting at the given index
void MEDDLY::unpacked_node::trim()
{
    if (!isExtensible()) return;
    if (isTrim()) return;

    // If extensible edge is transparent, mark the node as not-extensible and return
    if (d((isSparse()? getNNZs() : getSize()) - 1) == eparent->getTransparentNode()) {
        markAsNotExtensible();
        return;
    }

    MEDDLY_DCASSERT(isExtensible() && !isTrim());
    MEDDLY_DCASSERT(modparent);

    if (isSparse()) {
        unsigned z = getNNZs()-1;
        while (z > 0 && (i(z-1)+1) == i(z) && d(z-1) == d(z)) {
            modparent->unlinkNode(d(z));
            z--;
        }
        if (z != (getNNZs() - 1)) {
            // node is smaller than before, shrink it to the correct size.
            shrinkSparse(z+1);
        }
    } else {
        unsigned z = getSize()-1;
        while (z > 0 && d(z-1) == d(z)) {
            modparent->unlinkNode(d(z));
            z--;
        }
        if (z != (getSize()-1)) {
            shrinkFull(z+1);
        }
    }

    MEDDLY_DCASSERT(isExtensible() && isTrim());
}

// checks if the node is has no trailing redundant edges
bool MEDDLY::unpacked_node::isTrim() const
{
    if (!isExtensible()) return true;

    if (isSparse()) {
        unsigned nnz = getNNZs();
        return (nnz < 2 || i(nnz-1) != (i(nnz-2)+1) || d(nnz-1) != d(nnz-2));
    } else {
        unsigned sz = getSize();
        return (sz < 2 || d(sz-1) != d(sz-2));
    }
}

// checks if the node indices are in ascending order
bool MEDDLY::unpacked_node::isSorted() const
{
    if (!isSparse()) return true;

    for (unsigned z = 1; z < getNNZs(); z++) {
        if (i(z-1) >= i(z)) return false;
    }

    return true;
}

void MEDDLY::unpacked_node::clearEdges(unsigned stop)
{
    MEDDLY_DCASSERT(down);
    MEDDLY_DCASSERT(edge);
    if (edge_type::VOID == the_edge_type) {
        for (unsigned i=0; i<stop; i++) {
            down[i] = eparent->getTransparentNode();
        }
    } else {
        for (unsigned i=0; i<stop; i++) {
            eparent->getTransparentEdge(down[i], edge[i]);
        }
    }
}

void
MEDDLY::unpacked_node::freeRecycled()
{
    while (freeList) {
        MEDDLY::unpacked_node* n = freeList->next;
        delete freeList;
        freeList = n;
    }
}

void MEDDLY::unpacked_node::initStatics()
{
    freeList = nullptr;
    buildList = nullptr;
}


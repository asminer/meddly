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

#include "node_marker.h"
#include "operators.h"

// #define DEBUG_MARK

MEDDLY::node_marker::node_marker(const forest* F, array_watcher* w)
    : marked(w)
{
    S_top = nullptr;
    S_free = nullptr;

    For = F;
    if (For) {
        marked.expand(unsigned(For->getLastNode() + 2));
    }
    marked.clearAll();
}


MEDDLY::node_marker::~node_marker()
{
    while (S_top) {
        mystack* n = S_top->next;
        delete S_top;
        S_top = n;
    }
    while (S_free) {
        mystack* n = S_free->next;
        delete S_free;
        S_free = n;
    }
}

//
// Public methods
//

size_t MEDDLY::node_marker::countEdges() const
{
    MEDDLY_DCASSERT(For);

    size_t ec = 0;
    unpacked_node* M = unpacked_node::New(For);
    size_t i=0;
    while( (i=marked.firstOne(i+1)) < marked.getSize() )
    {
        For->unpackNode(M, i, FULL_ONLY);
        ec += M->getSize();
    }
    unpacked_node::Recycle(M);
    return ec;
}

size_t MEDDLY::node_marker::countNonzeroEdges() const
{
    MEDDLY_DCASSERT(For);

    size_t ec = 0;
    unpacked_node* M = unpacked_node::New(For);
    size_t i=0;
    while( (i=marked.firstOne(i+1)) < marked.getSize() )
    {
        For->unpackNode(M, i, SPARSE_ONLY);
        ec += M->getSize();
    }
    unpacked_node::Recycle(M);
    return ec;
}

void MEDDLY::node_marker::showByLevels(output &s) const
{
    MEDDLY_DCASSERT(For);

    unpacked_node* M = unpacked_node::New(For);

    const int lwid = (int) digits(For->getNumVariables());
    const int nwid = (int) digits(getSize());

    for (int k=int(For->getNumVariables()); k; k = MXD_levels::downLevel(k)) {

        bool level_printed = false;

        size_t i=0;
        while( (i=marked.firstOne(i+1)) < marked.getSize() )
        {
            if (For->getNodeLevel(i) != k) continue;
            //
            // Node i is at level k.
            //

            if (!level_printed) {
                //
                // Show level name
                //
                s << "Level: ";
                s.put(k, lwid);
                s << " Var: ";
                const variable* v = For->getDomain()->getVar(
                    unsigned(For->getVarByLevel(ABS(k)))
                );
                char primed = (k>0) ? ' ' : '\'';
                if (v->hasName()) {
                    s << v->getName() << primed << '\n';
                } else {
                    s << For->getVarByLevel(ABS(k)) << primed << '\n';
                }
                level_printed = true;
            }

            //
            // Show the node
            //
            s << "    node:";
            s.put(i, nwid);
            s << " incount: ";
            s.put(For->getNodeInCount(i));
            s.put(' ');
            For->unpackNode(M, i, FULL_OR_SPARSE);
            M->show(s, true);
            s.put('\n');

        } // for i

    } // for k
    unpacked_node::Recycle(M);
}

void MEDDLY::node_marker::getNodesAtLevel(int k, std::vector <node_handle> &v)
    const
{
    MEDDLY_DCASSERT(For);

    size_t i=0;
    while( (i=marked.firstOne(i+1)) < marked.getSize() )
    {
        if (For->getNodeLevel(i) != k) continue;
        v.push_back(i);
    } // for i
}

void MEDDLY::node_marker::getTerminals(std::set <node_handle> &v) const
{
    MEDDLY_DCASSERT(For);

    unpacked_node* M = unpacked_node::New(For);

    size_t i=0;
    while( (i=marked.firstOne(i+1)) < marked.getSize() )
    {
        For->unpackNode(M, i, SPARSE_ONLY);
        for (unsigned j=0; j<M->getSize(); j++) {
            if (M->down(j)>0) continue;
            // M->down(j) is a terminal node handle
            v.insert(M->down(j));
        }
    }

    unpacked_node::Recycle(M);
}

//
// Private methods
//

void MEDDLY::node_marker::_mark(node_handle p)
{
    MEDDLY_DCASSERT(For);

    CHECK_RANGE(__FILE__, __LINE__, 1, p, (node_handle) marked.getSize());
    const node_storage* nodeMan = For->getNodeManager();
    MEDDLY_DCASSERT(nodeMan);

    MEDDLY_DCASSERT(!S_top);

    push();
    addToQueue(p);

    while (S_top) {
#ifdef DEBUG_MARK
        debug(S_top);
#endif
        if (S_top->queue.empty()) {
#ifdef DEBUG_MARK
            std::cerr << "\tempty; popping\n";
#endif
            pop();
            continue;
        }

        p = S_top->queue.back();
        S_top->queue.pop_back();
#ifdef DEBUG_MARK
        std::cerr << "\texploring " << p << "\n";
#endif
        CHECK_RANGE(__FILE__, __LINE__, 1L, long(p), long(marked.getSize()));

        //
        // If S_top is empty, we can re-use it; otherwise
        // we need to push an empty queue onto the stack.
        //
        if (!S_top->queue.empty()) {
            push();
        }
        nodeMan->addDownToQueue( *this, For->getNodeAddress(p) );
    }

}

void MEDDLY::node_marker::debug(const mystack *s)
{
    std::cerr << "Exploring " << s << ": {";
    for (unsigned i=0; i<s->queue.size(); i++) {
        if (i) std::cerr << ", ";
        std::cerr << s->queue[i];
    }
    std::cerr << "}\n";
}

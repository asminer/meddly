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

#include "satur_index.h"
#include "../variable.h"
#include "../domain.h"
#include "../forest.h"
#include "../rel_node.h"
#include "../io.h"

// #define TRACE

// **********************************************************************

namespace MEDDLY {

    // helper objects

    class index_fifo;

    // various explorers

    class explore_fifo;
}

// **********************************************************************
// *                                                                    *
// *                     sat_index_explorer methods                     *
// *                                                                    *
// **********************************************************************

MEDDLY::sat_index_explorer::sat_index_explorer(forest* F, int lvl, bool fwd)
    : For(F), level(lvl), forwd(fwd)
{
    MEDDLY_DCASSERT(For);
    var = For->getDomain()->getVar(level);
    RN = nullptr;
    U = unpacked_node::New(For, SPARSE_ONLY);
}

MEDDLY::sat_index_explorer::~sat_index_explorer()
{
    if (RN) {
        For->doneRelNode(RN);
        RN = nullptr;
    }
    unpacked_node::Recycle(U);
}

void MEDDLY::sat_index_explorer::restart(node_handle n)
{
#ifdef TRACE
    std::cout << "saturating level " << level
              << " (mxd " << n << ")\n";
#endif
    //
    // Clear old
    //
    for (unsigned i=0; i<diagonals.size(); i++) {
        diagonals[i] = 0;
    }
    for (unsigned i=0; i<rows.size(); i++) {
        if (!rows[i].explored) continue;
        rows[i].elements.clear();
        rows[i].explored = false;
    }
    if (RN) {
        For->doneRelNode(RN);
    }

    //
    // Build new
    //
    RN = For->buildRelNode(n);

    // update size
    expandRows(var->getBound(!forwd));

    //
    // For forward exploration, build rows on the fly as needed.
    //
    if (forwd) return;

    //
    // For backward exploration, build everything at the beginning.
    // We explore the relation node elements [i,j]
    // and store the transpose [j,i] internally.
    //
    row_element elem;
    for (unsigned i=0; i<var->getBound(false); i++) {
        if (!RN->outgoing(i, *U)) continue;
        elem.index = i;

        for (unsigned z=0; z<U->getSize(); z++) {
            const unsigned j = U->index(z);
            if (j==i) {
                // Set diagonal
                diagonals[i] = U->down(z);
            } else {
                // append (i, down) to row j
                elem.down = U->down(z);
                rows[j].elements.push_back(elem);
            }
        } // for z
    } // for row i
    for (unsigned i=0; i<rows.size(); i++) {
        rows[i].explored = true;
    }

    finishAllRows();
}

void MEDDLY::sat_index_explorer::show(output &s) const
{
    s.put("explorer at level ");
    s.put(level);
    s.put("\n");
    for (unsigned i=0; i<rows.size(); i++) {
        s.put("Row ");
        s.put(i);
        s.put(":");
        if (!rows[i].explored) {
            s.put(" unexplored\n");
            continue;
        }
        s.put("\n");
        if (diagonals[i]) {
            s.put("    diag: ");
            s.put(diagonals[i]);
            s.put("\n");
        }
        for (unsigned z=0; z<rows[i].elements.size(); z++) {
            s.put("    ");
            s.put(rows[i].elements[z].index);
            s.put(": ");
            s.put(rows[i].elements[z].down);
            s.put("\n");
        }
    }
}

void MEDDLY::sat_index_explorer::clear()
{
    // Do nothing
}

void MEDDLY::sat_index_explorer::exploreRow(unsigned i)
{
    MEDDLY_DCASSERT(forwd);

    unsigned max_index = i;
    if (RN->outgoing(i, *U)) {
        for (unsigned z=0; z<U->getSize(); z++) {
            const unsigned j = U->index(z);
            if (j==i) {
                // Set diagonal
                diagonals[i] = U->down(z);
            } else {
                // append (j, down) to row i
                row_element elem;
                elem.index = j;
                elem.down = U->down(z);
                rows[i].elements.push_back(elem);
                max_index = MAX(max_index, j);
            }
        } // for z
    }
    rows[i].explored = true;
    if (max_index > rows.size()) {
        expandRows(1+max_index);
    }

    finishRow(i);
}

void MEDDLY::sat_index_explorer::finishRow(unsigned)
{
    // Do nothing
}

void MEDDLY::sat_index_explorer::finishAllRows()
{
    // Do nothing
}


void MEDDLY::sat_index_explorer::expandRows(unsigned newsz)
{
    if (newsz <= rows.size()) return;

    const unsigned oldsize = rows.size();

    rows.resize(newsz);
    diagonals.resize(newsz);

    finishExpandRows(oldsize, newsz);
}

void MEDDLY::sat_index_explorer::finishExpandRows(unsigned, unsigned)
{
    // Do nothing
}

void MEDDLY::sat_index_explorer::finishUpdate(unsigned)
{
    // Do nothing
}

// **********************************************************************
// *                                                                    *
// *                         index_fifo   class                         *
// *                                                                    *
// **********************************************************************

#define NULPTR -1
#define NOTINQ -2

class MEDDLY::index_fifo {
        /**
            A linked list of nodes.
            Element i is NOTINQ if it is not in the queue.
            Otherwise, element i is the index of the next item
            in the queue, or -1 if it is the end of the list.
         */
        std::vector <int> data;

        /// Index (in data) of the front of the queue.
        int head;
        /// Index (in data) of the end of the queue.
        int tail;
    public:
        index_fifo() {
            head = tail = NULPTR;
        }

        inline void resize(unsigned sz) {
            data.resize(sz, NOTINQ);
        }
        inline void clear() {
            data.assign(data.size(), NOTINQ);
        }
        inline bool isEmpty() const {
            return NULPTR == head;
        }
        inline int front() const {
            return (head >= 0) ? head : -1;
        }
        inline void add(int i) {
            MEDDLY_CHECK_RANGE(0, i, int(data.size()));
            if (NOTINQ != data[i]) {
                // already in queue
                return;
            }
            if (NULPTR == head) {
                // queue is empty
                head = i;
            } else {
                // queue is not empty
                MEDDLY_CHECK_RANGE(0, tail, int(data.size()));
                data[tail] = i;
            }
            tail = i;
            data[i] = NULPTR;
        }
        inline int remove() {
            MEDDLY_CHECK_RANGE(0, head, int(data.size()));
            int ans = head;
            head = data[head];
            data[ans] = NOTINQ;
            return ans;
        }
};

// **********************************************************************
// *                                                                    *
// *                        explore_fifo   class                        *
// *                                                                    *
// **********************************************************************

class MEDDLY::explore_fifo : public sat_index_explorer {
    public:
        explore_fifo(forest* _F, int _level, bool _forwd);
        virtual ~explore_fifo();
        virtual bool nextEdge(unsigned &i, unsigned &j, node_handle &down);

    protected:
        virtual void clear();
        virtual void finishExpandRows(unsigned oldsz, unsigned newsz);
        virtual void finishUpdate(unsigned i);

    private:
        index_fifo queue;
        int rptr;
};

MEDDLY::explore_fifo::explore_fifo(forest* _F, int _level, bool _forwd)
    : sat_index_explorer(_F, _level, _forwd)
{
    rptr = -1;
}

MEDDLY::explore_fifo::~explore_fifo()
{
}

bool MEDDLY::explore_fifo::nextEdge(unsigned &i, unsigned &j, node_handle &dn)
{
    for (;;) {
        if (queue.isEmpty()) return false;

        int ii = queue.front();
        MEDDLY_DCASSERT(ii>=0);
        i = unsigned(ii);

        if (rptr < 0) {
            rptr = 0;
            if (diagonals[i]) {
                j = i;
                dn = diagonals[i];
                break;
            }
        }

        if (rptr >= rows[i].elements.size()) {
            queue.remove();
            rptr = -1;
            continue;
        }

        j = rows[i].elements[rptr].index;
        dn = rows[i].elements[rptr].down;
        rptr++;
        break;
    }

#ifdef TRACE
    std::cout << "level " << level << ": " << i << " -> " << j
              << " (down " << dn << ")\n";
#endif
    return true;
}

void MEDDLY::explore_fifo::clear()
{
    queue.clear();
}

void MEDDLY::explore_fifo::finishExpandRows(unsigned oldsz, unsigned newsz)
{
    queue.resize(newsz);
}

void MEDDLY::explore_fifo::finishUpdate(unsigned i)
{
    queue.add(int(i));
}


// **********************************************************************
// *                                                                    *
// *                             Front  end                             *
// *                                                                    *
// **********************************************************************

MEDDLY::sat_index_explorer* MEDDLY::makeSatIndexExplorer(char which,
        forest* F, int level, bool forwd)
{
    switch (which) {
        // default:

        case 'q':
        default:

            return new explore_fifo(F, level, forwd);

    }

    // Fail safe
    return nullptr;
}


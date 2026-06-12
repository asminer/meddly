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

#ifndef MEDDLY_SATUR_GRAPH_H
#define MEDDLY_SATUR_GRAPH_H

#include "../defines.h"
#include <vector>

namespace MEDDLY {
    class output;
    class variable;
    class forest;
    class rel_node;
    class unpacked_node;

    class satur_graph {
        public:
            struct sparse_element {
                int index;
                node_handle down;

                sparse_element(int ndx, node_handle dn) {
                    index = ndx;
                    down = dn;
                }
            };

        public:
            satur_graph();
            ~satur_graph();

            void attach(forest* F, int level, bool forwd);
            inline void restart(node_handle n) {
                if (node != n) _restart(n);
            }

            inline forest* getForest() const {
                return For;
            }
            inline unsigned numRows() const {
                return rowptr.size();
            }

            inline bool isRowUnexplored(unsigned i) const {
#ifdef DEVELOPMENT_CODE
                return rowptr.at(i) < 0;
#else
                return rowptr[i] < 0;
#endif
            }
            inline int getRowStart(unsigned i) const {
#ifdef DEVELOPMENT_CODE
                MEDDLY_DCASSERT(rowptr.at(i) >= 0);
#endif
                return rowptr[i];
            }
            inline bool getElement(int p, unsigned &j, node_handle &dn) const {
                MEDDLY_DCASSERT(p>=0);
#ifdef DEVELOPMENT_CODE
                if (elements.at(p).index < 0) return false;
                j = unsigned(elements.at(p).index);
                dn = elements.at(p).down;
#else
                if (elements[p].index < 0) return false;
                j = unsigned(elements[p].index);
                dn = elements[p].down;
#endif
                return true;
            }
            inline node_handle getDiagonal(unsigned i) {
#ifdef DEVELOPMENT_CODE
                return diagonals.at(i) < 0;
#else
                return diagonals[i] < 0;
#endif
            }
            inline void ensureRowExplored(unsigned i) {
                if (isRowUnexplored(i)) {
                    exploreRow(i);
                }
            }

            void show(output &s) const;

        protected:
            void _restart(node_handle n);
            void exploreRow(unsigned i);
            void buildTranspose();

            inline void expandRows(unsigned sz) {
                if (sz > rowptr.size()) {
                    diagonals.resize(sz, 0);
                    rowptr.resize(sz, -1);
                }
            }

        private:
            // Graph data
            std::vector <int> rowptr;
            std::vector <sparse_element> elements;
            std::vector <node_handle> diagonals;

            // MxD data
            forest* For;
            variable* var;
            rel_node* RN;
            unpacked_node* U;

            node_handle node;
            int level;
            bool forwd;
    };

};

#endif

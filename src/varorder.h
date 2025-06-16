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

#ifndef MEDDLY_VARORDER_H
#define MEDDLY_VARORDER_H

#include <vector>

namespace MEDDLY {
    class variable;
    class variable_order;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      variable order class                      *
// *                                                                *
// *                                                                *
// ******************************************************************

class MEDDLY::variable_order {
    protected:
        std::vector<int> level2var;
        std::vector<int> var2level;

    public:
        variable_order(const int* order, int size);
        variable_order(const variable_order& order);

        inline int getVarByLevel(int level) const {
            ASSERT(__FILE__, __LINE__, level>=0);
            return level2var[std::size_t(level)];
        }
        inline int getLevelByVar(int var) const {
            ASSERT(__FILE__, __LINE__, var>=0);
            return var2level[std::size_t(var)];
        }

        // Exchange two variables
        // The two variables don't have to be adjacent
        void exchange(int var1, int var2);

        bool is_compatible_with(const int* order) const;
        bool is_compatible_with(const variable_order& order) const;
};

#endif

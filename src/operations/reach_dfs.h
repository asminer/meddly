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

#ifndef MEDDLY_REACH_DFS_H
#define MEDDLY_REACH_DFS_H

namespace MEDDLY {
    class forest;
    class binary_operation;
    class binary_list;

    /// The 'dfs forward reachable' operation builder.
    binary_operation* REACHABLE_STATES_DFS(forest* a, forest* b, forest* c);
    void REACHABLE_STATES_DFS_init();
    void REACHABLE_STATES_DFS_done();

    /// The 'dfs backward reachable' operation builder.
    binary_operation* REVERSE_REACHABLE_DFS(forest* a, forest* b, forest* c);
    void REVERSE_REACHABLE_DFS_init();
    void REVERSE_REACHABLE_DFS_done();
}

#endif


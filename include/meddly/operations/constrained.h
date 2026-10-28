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

#ifndef MEDDLY_SAT_CONSTRAINED_H
#define MEDDLY_SAT_CONSTRAINED_H

namespace MEDDLY {
    class forest;
    class ternary_operation;

    ternary_operation* CONSTRAINED_BACKWARD_BFS(forest* consF, forest* inF,
            forest* relF, forest* outF);

    void CONSTRAINED_BACKWARD_BFS_init();
    void CONSTRAINED_BACKWARD_BFS_done();

    ternary_operation* CONSTRAINED_FORWARD_DFS(forest* consF, forest* inF,
            forest* relF, forest* outF);

    void CONSTRAINED_FORWARD_DFS_init();
    void CONSTRAINED_FORWARD_DFS_done();

    ternary_operation* CONSTRAINED_BACKWARD_DFS(forest* consF, forest* inF,
            forest* relF, forest* outF);

    void CONSTRAINED_BACKWARD_DFS_init();
    void CONSTRAINED_BACKWARD_DFS_done();

}

#endif

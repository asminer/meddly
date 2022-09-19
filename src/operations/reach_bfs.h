
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

#ifndef REACH_BFS_H
#define REACH_BFS_H

namespace MEDDLY {
  class binary_opname;

  class binary_opname_event;
  /// Set up a binary_opname for the "reachable bfs" operation.
  binary_opname* initializeForwardBFS();

  /// Set up a binary_opname for the "reverse reachable bfs" operation.
  binary_opname* initializeBackwardBFS();

  /// Set up a binary_opname for the  "reachable bfs with under approximation" opertaion.
  binary_opname* initializeForwardBFSUA();

  /// Set up a binary_opname for the  "reachable bfs with heuristic under approximation" opertaion.
  binary_opname* initializeForwardBFSHUA();

  binary_opname_event* initializeAllBFSGen();

  binary_opname_event* initializeAllBFSGenHUA();
}

#endif

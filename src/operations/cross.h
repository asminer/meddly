
// $Id$

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

#ifndef CROSS_H
#define CROSS_H

namespace MEDDLY {

  class old_operation;
  class op_param;

  /** Minimalist front-end interface.
      Return the appropriate operation for 
      the given forest and desired return type.
        @param  ot  Operand type
        @return The appropriate operation, or 0 on error.
  */
  old_operation* getCrossOperation(const op_param &ot);

}

#endif



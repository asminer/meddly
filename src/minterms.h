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

#ifndef MEDDLY_MINTERMS_H
#define MEDDLY_MINTERMS_H

//
// Idea: define minterm objects, for both full and sparse storage,
// and for set vs relation.
//
// For now: defines that should eventually be part of minterm objects
//

//
// A minterm should be a std::vector of ints, plus a rangeval.
//
// Should we define collections of minterms here?
//

namespace MEDDLY {
  /** Special value for minterms: don't care what this variable does.
      I.e., do the same thing for all possible assignments for a variable.
  */
  const int DONT_CARE  = -1;
  /** Special value for primed minterms: don't change this variable.
      Forces the primed value to equal the unprimed value for a variable.
      Undefined for unprimed variables.
  */
  const int DONT_CHANGE = -2;

};  // namespace MEDDLY

#endif

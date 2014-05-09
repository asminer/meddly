
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



/*! \file compute_table.h

    Compute table interface.

    This interface is for "expert" interface users who wish to implement
    user-defined operations.

    A compute table is used to cache the results of operations on MDDs.
    An expert-user wishing to cache results of user-defined operations may
    need to use this interface.

*/

#ifndef COMPUTE_TABLE_H
#define COMPUTE_TABLE_H

#if 0
#include "defines.h"

namespace MEDDLY {

  /** Build a new, monolithic table.
      Monolithic means that the table stores entries for several
      (ideally, all) operations.
  */
  compute_table* createMonolithicTable(
    const settings::computeTableSettings &s
  ); 

  /** Build a new table for a single operation.
  */
  compute_table* createOperationTable(
    const settings::computeTableSettings &s, operation* op
  );
}
#endif

#endif


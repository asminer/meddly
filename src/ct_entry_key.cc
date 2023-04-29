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

#include "ct_entry_key.h"


MEDDLY::ct_entry_key::ct_entry_key()
{
  data_alloc = 8;
  etype = 0;
  data = (ct_entry_item*) malloc(data_alloc * sizeof(ct_entry_item));
  temp_data = 0;
  temp_bytes = 0;
  temp_alloc = 0;
  // malloc: because realloc later
}

MEDDLY::ct_entry_key::~ct_entry_key()
{
  free(data);
  free(temp_data);
}

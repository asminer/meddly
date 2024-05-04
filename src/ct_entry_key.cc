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

// **********************************************************************
// *                                                                    *
// *                        ct_entry_key methods                        *
// *                                                                    *
// **********************************************************************


MEDDLY::ct_entry_key::ct_entry_key()
{
    data_alloc = 8;
    etype = nullptr;
    data = (ct_entry_item*) malloc(data_alloc * sizeof(ct_entry_item));
    temp_data = nullptr;
    temp_bytes = 0;
    temp_alloc = 0;
    // malloc: because realloc later

    my_entry = 0;
    entry_slots = 0;
    result_shift = 0;
}

MEDDLY::ct_entry_key::~ct_entry_key()
{
    free(data);
    free(temp_data);
}

void MEDDLY::ct_entry_key::setup(const ct_entry_type* et, unsigned repeats)
{
    MEDDLY_DCASSERT(et);
    etype = et;
    num_repeats = repeats;
    MEDDLY_DCASSERT( 0==repeats || et->isRepeating() );
    total_slots = et->getKeySize(repeats);
    if (total_slots > data_alloc) {
        // Allocate in chunks of size 8
        data_alloc = (1+(data_alloc / 8)) * 8;
        data = (ct_entry_item*) realloc(data, data_alloc*sizeof(ct_entry_item));
        if (!data) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
    }
    memset(data, 0, total_slots * sizeof(ct_entry_item));
    currslot = 0;
#ifdef DEVELOPMENT_CODE
    has_hash = false;
#endif
}

void* MEDDLY::ct_entry_key::allocTempData(unsigned bytes)
{
    temp_bytes = bytes;
    if (bytes > temp_alloc) {
        // Allocate in chunks of 64 bytes
        temp_alloc = (1+(temp_bytes/64)) * 64;
        temp_data = realloc(temp_data, temp_alloc);
        if (!temp_data) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
    }
    return temp_data;
}


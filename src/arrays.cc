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

#include "io.h"
#include "error.h"

#include "arrays.h"

// ******************************************************************
// *                                                                *
// *                     array_watcher  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::array_watcher::array_watcher()
{
}

MEDDLY::array_watcher::~array_watcher()
{
}

void MEDDLY::array_watcher::expandElementSize(unsigned, unsigned)
{
}

void MEDDLY::array_watcher::shrinkElementSize(unsigned, unsigned)
{
}

// ******************************************************************
// *                                                                *
// *                      level_array  methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::level_array::level_array(int max_level, array_watcher* w)
{
    watch = w;

    data8 = nullptr;
    data16 = nullptr;
    data32 = nullptr;
    size = 0;

    if (max_level < 128) {
        bytes = 1;
        MEDDLY_DCASSERT(1 == sizeof(char));
    } else if (max_level < 32768) {
        bytes = 2;
        MEDDLY_DCASSERT(2 == sizeof(short));
    } else {
        bytes = 4;
        MEDDLY_DCASSERT(4 == sizeof(int));
    }

    if (watch) watch->expandElementSize(0, bytes*8);
}

MEDDLY::level_array::~level_array()
{
    free(data8);
    free(data16);
    free(data32);
    if (watch) watch->shrinkElementSize(bytes*8, 0);
}

void MEDDLY::level_array::expand(size_t ns)
{
    if (ns <= size) return;

    char* d8 = nullptr;
    short* d16 = nullptr;
    int* d32 = nullptr;

    switch (bytes) {
        case 1:   // char array
            d8 = (char*) realloc(data8, ns * bytes);
            if (!d8) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            memset(d8 + size, 0, (ns-size) * bytes );
            data8 = d8;
            break;

        case 2:   // short array
            d16 = (short*) realloc(data16, ns * bytes);
            if (!d16) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            memset(d16 + size, 0, (ns-size) * bytes );
            data16 = d16;
            break;

        case 4:   // int array
            d32 = (int*) realloc(data32, ns * bytes);
            if (!d32) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            memset(d32 + size, 0, (ns-size) * bytes );
            data32 = d32;
            break;

        default:
            MEDDLY_DCASSERT(0);
    }

    size = ns;
#ifdef DEBUG_LEVEL_RESIZE
    std::cerr << "Enlarged level array, new size " << size << std::endl;
#endif
}

void MEDDLY::level_array::shrink(size_t ns)
{
    if (ns >= size) return;

    char* d8 = nullptr;
    short* d16 = nullptr;
    int* d32 = nullptr;

    switch (bytes) {
        case 1:   // char array
            d8 = (char*) realloc(data8, ns * bytes);
            if (!d8) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            data8 = d8;
            break;

        case 2:   // short array
            d16 = (short*) realloc(data16, ns * bytes);
            if (!d16) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            data16 = d16;
            break;

        case 4:   // int array
            d32 = (int*) realloc(data32, ns * bytes);
            if (!d32) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            data32 = d32;
            break;

        default:
            MEDDLY_DCASSERT(0);
    }

    size = ns;
#ifdef DEBUG_LEVEL_RESIZE
    std::cerr << "Reduced level array, new size " << size << std::endl;
#endif
}

void MEDDLY::level_array::show(output &s, size_t first, size_t last,
        int width) const
{
    s.put('[');
    s.put(long(get(first)), width);
    for (size_t i=first+1; i<=last; i++) {
        s.put('|');
        s.put(long(get(i)), width);
    }
    s.put(']');
}



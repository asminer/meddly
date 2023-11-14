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

// #define DEBUG_COUNTER_RESIZE
// #define DEBUG_ADDRESS_RESIZE


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

// ******************************************************************
// *                                                                *
// *                     counter_array  methods                     *
// *                                                                *
// ******************************************************************


MEDDLY::counter_array::counter_array(array_watcher *w)
{
    watch = w;

    data8 = nullptr;
    data16 = nullptr;
    data32 = nullptr;

    size = 0;
    counts_09bit = 0;
    counts_17bit = 0;
    bytes = sizeof(unsigned char);
    if (watch) watch->expandElementSize(0, bytes*8);
}

MEDDLY::counter_array::~counter_array()
{
    free(data8);
    free(data16);
    free(data32);
    if (watch) watch->shrinkElementSize(bytes*8, 0);
}

void MEDDLY::counter_array::expand(size_t ns)
{
    if (ns <= size) return;

    unsigned char* d8 = nullptr;
    unsigned short* d16 = nullptr;
    unsigned int* d32 = nullptr;

    switch (bytes) {
        case 1:   // char array
            MEDDLY_DCASSERT(0==data16);
            MEDDLY_DCASSERT(0==data32);
            d8 = (unsigned char*) realloc(data8, ns * bytes);
            if (!d8) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            memset(d8 + size, 0, (ns-size) * bytes );
            data8 = d8;
            break;

        case 2:   // short array
            if (counts_09bit) {
                d16 = (unsigned short*) realloc(data16, ns * bytes);
                if (!d16) {
                    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
                }
                memset(d16 + size, 0, (ns-size) * bytes );
                data16 = d16;
                MEDDLY_DCASSERT(0==data8);
                MEDDLY_DCASSERT(0==data32);
            } else {
                shrink16to8(ns);
                memset(data8 + size, 0, (ns-size) * bytes );
            }
            break;

        case 4:   // int array
            if (counts_17bit) {
                d32 = (unsigned int*) realloc(data32, ns * bytes);
                if (!d32) {
                    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
                }
                memset(d32 + size, 0, (ns-size) * bytes );
                data32 = d32;
                MEDDLY_DCASSERT(0==data8);
                MEDDLY_DCASSERT(0==data16);
            } else if (counts_09bit) {
                shrink32to16(ns);
                memset(data16 + size, 0, (ns-size) * bytes );
            } else {
                shrink32to8(ns);
                memset(data8 + size, 0, (ns-size) * bytes );
            }
            break;

        default:
            MEDDLY_DCASSERT(0);
    }; // switch bytes

    size = ns;
#ifdef DEBUG_COUNTER_RESIZE
    std::cerr << "Enlarged counter array, new size " << size << std::endl;
#endif
}

void MEDDLY::counter_array::shrink(size_t ns)
{
    if (ns >= size) return;

    unsigned char* d8 = nullptr;
    unsigned short* d16 = nullptr;
    unsigned int* d32 = nullptr;

    switch (bytes) {
        case 1:   // char array
            MEDDLY_DCASSERT(0==counts_09bit);
            MEDDLY_DCASSERT(0==counts_17bit);
            d8 = (unsigned char*) realloc(data8, ns * bytes);
            if (!d8) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            data8 = d8;
            break;

        case 2:   // short array
            MEDDLY_DCASSERT(0==counts_17bit);
            if (counts_09bit) {
                d16 = (unsigned short*) realloc(data16, ns * bytes);
                if (!d16) {
                    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
                }
                data16 = d16;
            } else {
                shrink16to8(ns);
            }
            break;

        case 4:   // int array
            if (counts_17bit) {
                d32 = (unsigned int*) realloc(data32, ns * bytes);
                if (!d32) {
                    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
                }
                data32 = d32;
            } else if (counts_09bit) {
                shrink32to16(ns);
            } else {
                shrink32to8(ns);
            }
            break;

        default:
            MEDDLY_DCASSERT(0);
    }; // switch bytes

    size = ns;
#ifdef DEBUG_COUNTER_RESIZE
    std::cerr << "Reduced counter array, new size " << size << std::endl;
#endif
}

void MEDDLY::counter_array::show(output &s, size_t first, size_t last,
        int width) const
{
    s.put('[');
    s.put((unsigned long) get(first), width);
    for (size_t i=first+1; i<=last; i++) {
        s.put('|');
        s.put((unsigned long) get(i), width);
    }
    s.put(']');
}

void MEDDLY::counter_array::expand8to16(size_t j)
{
#ifdef DEBUG_COUNTER_RESIZE
    std::cerr << "Widening counter array from 8 to 16 bits (due to element "
              << j << ")\n";
#endif
    MEDDLY_DCASSERT(1==bytes);
    MEDDLY_DCASSERT(!data16);
    MEDDLY_DCASSERT(!data32);
    MEDDLY_DCASSERT(!counts_09bit);
    MEDDLY_DCASSERT(!counts_17bit);

    counts_09bit = 1;
    bytes = sizeof(unsigned short);
    data16 = (unsigned short*) malloc(size * bytes);
    if (!data16) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    for (size_t i=0; i<size; i++) {
        data16[i] = data8[i];
    }
    data16[j] = 256;
    free(data8);
    data8 = nullptr;
    if (watch) watch->expandElementSize(8, 16);
}

void MEDDLY::counter_array::expand16to32(size_t j)
{
#ifdef DEBUG_COUNTER_RESIZE
    std::cerr << "Widening counter array from 16 to 32 bits (due to element "
              << j << ")\n";
#endif
    MEDDLY_DCASSERT(2==bytes);
    MEDDLY_DCASSERT(!data8);
    MEDDLY_DCASSERT(!data32);
    MEDDLY_DCASSERT(!counts_17bit);

    counts_17bit = 1;
    bytes = sizeof(unsigned int);
    data32 = (unsigned int*) malloc(size * bytes);
    if (!data32) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    for (size_t i=0; i<size; i++) {
        data32[i] = data16[i];
    }
    data32[j] = 65536;
    free(data16);
    data16 = nullptr;
    if (watch) watch->expandElementSize(16, 32);
}

void MEDDLY::counter_array::shrink16to8(size_t ns)
{
#ifdef DEBUG_COUNTER_RESIZE
    std::cerr << "Narrowing counter array from 16 to 8 bits\n";
#endif
    MEDDLY_DCASSERT(2==bytes);
    MEDDLY_DCASSERT(data16);
    MEDDLY_DCASSERT(!data8);
    MEDDLY_DCASSERT(!data32);
    MEDDLY_DCASSERT(!counts_09bit);
    MEDDLY_DCASSERT(!counts_17bit);

    bytes = sizeof(unsigned char);
    data8 = (unsigned char*) malloc(ns * bytes);
    if (!data8) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    size_t stop = (ns < size) ? ns : size;
    for (size_t i=0; i<stop; i++) {
        data8[i] = data16[i];
    }
    free(data16);
    data16 = nullptr;
    if (watch) watch->shrinkElementSize(16, 8);
}

void MEDDLY::counter_array::shrink32to16(size_t ns)
{
#ifdef DEBUG_COUNTER_RESIZE
    std::cerr << "Narrowing counter array from 16 to 8 bits\n";
#endif
    MEDDLY_DCASSERT(4==bytes);
    MEDDLY_DCASSERT(data32);
    MEDDLY_DCASSERT(!data8);
    MEDDLY_DCASSERT(!data16);
    MEDDLY_DCASSERT(!counts_17bit);

    bytes = sizeof(unsigned short);
    data16 = (unsigned short*) malloc(ns * bytes);
    if (!data16) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    size_t stop = (ns < size) ? ns : size;
    for (size_t i=0; i<stop; i++) {
        data16[i] = data32[i];
    }
    free(data32);
    data32 = nullptr;
    if (watch) watch->shrinkElementSize(32, 16);
}

void MEDDLY::counter_array::shrink32to8(size_t ns)
{
#ifdef DEBUG_COUNTER_RESIZE
    std::cerr << "Narrowing counter array from 16 to 8 bits\n";
#endif
    MEDDLY_DCASSERT(4==bytes);
    MEDDLY_DCASSERT(data32);
    MEDDLY_DCASSERT(!data8);
    MEDDLY_DCASSERT(!data16);
    MEDDLY_DCASSERT(!counts_09bit);
    MEDDLY_DCASSERT(!counts_17bit);

    bytes = sizeof(unsigned char);
    data8 = (unsigned char*) malloc(ns * bytes);
    if (!data8) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    size_t stop = (ns < size) ? ns : size;
    for (size_t i=0; i<stop; i++) {
        data8[i] = data32[i];
    }
    free(data32);
    data32 = nullptr;
    if (watch) watch->shrinkElementSize(32, 8);
}


// ******************************************************************
// *                                                                *
// *                     address_array  methods                     *
// *                                                                *
// ******************************************************************


MEDDLY::address_array::address_array(array_watcher *w)
{
    watch = w;
    data32 = nullptr;
    data64 = nullptr;
    num_large_elements = 0;
    bytes = 4;
    size = 0;
    if (watch) watch->expandElementSize(0, bytes*8);
}

MEDDLY::address_array::~address_array()
{
    free(data32);
    free(data64);
    if (watch) watch->shrinkElementSize(bytes*8, 0);
}

void MEDDLY::address_array::expand(size_t ns)
{
    if (ns <= size) return;

    if (4==bytes) {
        MEDDLY_DCASSERT(!data64);
        unsigned int* d = (unsigned int*) realloc(data32, ns * bytes);
        if (!d) {
            throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
        memset(d + size, 0, (ns-size) * bytes);
        data32 = d;
    } else {
        MEDDLY_DCASSERT(!data32);
        MEDDLY_DCASSERT(8==bytes);
        if (num_large_elements) {
            // Stay at 64 bits
            unsigned long* d = (unsigned long*) realloc(data64, ns * bytes);
            if (!d) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            memset(d + size, 0, (ns-size) * sizeof(unsigned long) );
            data64 = d;
        } else {
            shrink64to32(ns);
        }
    }

    size = ns;
#ifdef DEBUG_ADDRESS_RESIZE
    std::cerr << "Enlarged address array, new size " << size << std::endl;
#endif
}

void MEDDLY::address_array::shrink(size_t ns)
{
    if (ns >= size) return;

    if (4==bytes) {
        MEDDLY_DCASSERT(data32);
        MEDDLY_DCASSERT(!data64);
        unsigned int* d = (unsigned int*) realloc(data32, ns * bytes);
        if (!d) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
        }
        data32 = d;
    } else {
        MEDDLY_DCASSERT(8==bytes);
        MEDDLY_DCASSERT(0==data32);
        if (num_large_elements) {
            // Stay at 64 bits
            unsigned long* d = (unsigned long*)
                realloc(data64, ns * sizeof(unsigned long));
            if (!d) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            data64 = d;
        } else {
            shrink64to32(ns);
        }
    }

    size = ns;
#ifdef DEBUG_ADDRESS_RESIZE
    std::cerr << "Reduced address array, new size " << size << std::endl;
#endif
}

void MEDDLY::address_array::show(output &s, size_t first, size_t last,
                int width) const
{
    s.put('[');
    s.put(get(first), width);
    for (size_t i=first+1; i<=last; i++) {
        s.put('|');
        s.put(get(i), width);
    }
    s.put(']');
}

void MEDDLY::address_array::expand32to64()
{
#ifdef DEBUG_ADDRESS_RESIZE
    std::cerr << "Widening counter array from 32 to 64 bits\n";
#endif
    MEDDLY_DCASSERT(4==bytes);
    MEDDLY_DCASSERT(!data64);
    MEDDLY_DCASSERT(!num_large_elements);

    bytes = 8;
    MEDDLY_DCASSERT(sizeof(unsigned long) == bytes);
    data64 = (unsigned long*) malloc(size * bytes);
    if (!data64) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    for (size_t i=0; i<size; i++) {
        data64[i] = data32[i];
    }
    free(data32);
    data32 = nullptr;
    num_large_elements = 1;
    if (watch) watch->expandElementSize(32, 64);
}

void MEDDLY::address_array::shrink64to32(size_t ns)
{
#ifdef DEBUG_ADDRESS_RESIZE
    std::cerr << "Narrowing counter array from 64 to 32 bits\n";
#endif
    MEDDLY_DCASSERT(8==bytes);
    MEDDLY_DCASSERT(data64);
    MEDDLY_DCASSERT(!data32);
    MEDDLY_DCASSERT(!num_large_elements);

    bytes = 4;
    MEDDLY_DCASSERT(sizeof(unsigned int) == bytes);
    data32 = (unsigned int*) malloc(ns * bytes);
    if (!data32) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    size_t stop = (ns < size) ? ns : size;
    for (size_t i=0; i<stop; i++) {
        data32[i] = data64[i];
    }
    free(data64);
    data64 = nullptr;
    if (watch) watch->shrinkElementSize(64, 32);
}


// ******************************************************************
// *                                                                *
// *                       bitvector  methods                       *
// *                                                                *
// ******************************************************************


MEDDLY::bitvector::bitvector(array_watcher *w)
{
    watch = w;
    data = nullptr;
    size = 0;
    if (watch) watch->expandElementSize(0, sizeof(bool)*8);
}

MEDDLY::bitvector::~bitvector()
{
    free(data);
    if (watch) watch->shrinkElementSize(sizeof(bool)*8, 0);
}

void MEDDLY::bitvector::expand(size_t ns)
{
    if (ns <= size) return;

    bool* d = (bool*) realloc(data, ns * sizeof(bool));
    if (!d) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    memset(d + size, 0, (ns-size) * sizeof(bool) );
    size = ns;
    data = d;
#ifdef DEBUG_BITVECTOR_RESIZE
    std::cerr << "Enlarged bitvector, new size " << size << std::endl;
#endif
}

void MEDDLY::bitvector::shrink(size_t ns)
{
    if (ns >= size) return;
    bool* d = (bool*) realloc(data, ns * sizeof(bool));
    if (!d) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    size = ns;
    data = d;
#ifdef DEBUG_BITVECTOR_RESIZE
    std::cerr << "Reduced bitvector, new size " << size << std::endl;
#endif
}


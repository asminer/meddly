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

#ifndef MEDDLY_ARRAYS_H
#define MEDDLY_ARRAYS_H

#include "defines.h"
#include <cstring>  // for memset

namespace MEDDLY {
    class output;

    class array_watcher;

    class level_array;
    class counter_array;
};


// ******************************************************************
// *                                                                *
// *                      array_watcher  class                      *
// *                                                                *
// ******************************************************************

/*
    Used to receive notifications when an array changes size.
    The default behavior here is to do nothing; derive a class
    from this one to change that.
*/
class MEDDLY::array_watcher {
    public:
        array_watcher();
        virtual ~array_watcher();

        /**
            Called when array elements get larger.
                @param  oldbits     Previous number of bits per element
                @param  newbits     New number of bits per element
        */
        virtual void expandElementSize(unsigned oldbits, unsigned newbits);

        /**
            Called when array elements get smaller.
                @param  oldbits     Previous number of bits per element
                @param  newbits     New number of bits per element
        */
        virtual void shrinkElementSize(unsigned oldbits, unsigned newbits);
};


// ******************************************************************
// *                                                                *
// *                       level_array  class                       *
// *                                                                *
// ******************************************************************

/**
    An array of levels.
    The element size is chosen dynamically based on the number of levels,
    but this is done only once during object construction.
*/
class MEDDLY::level_array {
        array_watcher* watch;
        char* data8;
        short* data16;
        int* data32;
        size_t size;
        unsigned bytes;
    public:
        level_array(int max_level, array_watcher* w=nullptr);
        ~level_array();

        void expand(size_t ns);
        void shrink(size_t ns);

        void show(output &s, size_t first, size_t last, int width) const;

        inline size_t entry_bits() const {
            return bytes * 8;
        }

        inline int get(size_t i) const {
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), i, size);
            if (data8) {
                MEDDLY_DCASSERT(!data16);
                MEDDLY_DCASSERT(!data32);
                return data8[i];
            }
            if (data16) {
                MEDDLY_DCASSERT(!data32);
                return data16[i];
            }
            MEDDLY_DCASSERT(data32);
            return data32[i];
        }

        inline void set(size_t i, int v) {
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), i, size);
            if (data8) {
                CHECK_RANGE(__FILE__, __LINE__, -128, v, 128);
                MEDDLY_DCASSERT(!data16);
                MEDDLY_DCASSERT(!data32);
                data8[i] = v;
                return;
            }
            if (data16) {
                CHECK_RANGE(__FILE__, __LINE__, -32768, v, 32768);
                MEDDLY_DCASSERT(!data32);
                data16[i] = v;
                return;
            }
            MEDDLY_DCASSERT(data32);
            data32[i] = v;
        }

        inline void swap(size_t i, size_t j) {
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), i, size);
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), j, size);
            if (data8) {
                MEDDLY_DCASSERT(!data16);
                MEDDLY_DCASSERT(!data32);
                SWAP(data8[i], data8[j]);
                return;
            }
            if (data16) {
                MEDDLY_DCASSERT(!data32);
                SWAP(data16[i], data16[j]);
                return;
            }
            MEDDLY_DCASSERT(data32);
            SWAP(data32[i], data32[j]);
        }

};


// ******************************************************************
// *                                                                *
// *                      counter_array  class                      *
// *                                                                *
// ******************************************************************

/**
    An array of counters.
    The element size is expanded as needed, based on the maximum count.
*/
class MEDDLY::counter_array {
        array_watcher* watch;
        unsigned char* data8;
        unsigned short* data16;
        unsigned int* data32;
        size_t size;
        size_t counts_09bit;  // #counts requiring at least 9 bits
        size_t counts_17bit;  // #counts requiring at least 17 bits
        unsigned bytes;
    public:
        counter_array(array_watcher* w = nullptr);
        ~counter_array();

        void expand(size_t ns);
        void shrink(size_t ns);

        void show(output &s, size_t first, size_t last, int width) const;

        inline size_t entry_bits() const { return bytes * 8; }

        inline unsigned int get(size_t i) const {
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), i, size);
            if (data8) {
                MEDDLY_DCASSERT(!data16);
                MEDDLY_DCASSERT(!data32);
                return data8[i];
            }
            if (data16) {
                MEDDLY_DCASSERT(!data32);
                return data16[i];
            }
            MEDDLY_DCASSERT(data32);
            return data32[i];
        }

        inline void swap(size_t i, size_t j) {
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), i, size);
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), j, size);
            if (data8) {
                MEDDLY_DCASSERT(!data16);
                MEDDLY_DCASSERT(!data32);
                SWAP(data8[i], data8[j]);
                return;
            }
            if (data16) {
                MEDDLY_DCASSERT(!data32);
                SWAP(data16[i], data16[j]);
                return;
            }
            MEDDLY_DCASSERT(data32);
            SWAP(data32[i], data32[j]);
        }

        inline void increment(size_t i) {
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), i, size);
            if (data8) {
                MEDDLY_DCASSERT(!data16);
                MEDDLY_DCASSERT(!data32);
                if (0 == ++data8[i]) expand8to16(i);
                return;
            }
            if (data16) {
                MEDDLY_DCASSERT(!data32);
                ++data16[i];
                if (256 == data16[i]) ++counts_09bit;
                if (0 == data16[i]) expand16to32(i);
                return;
            }
            MEDDLY_DCASSERT(data32);
            data32[i]++;
            if (256 == data32[i]) ++counts_09bit;
            if (65536 == data32[i]) ++counts_17bit;
        }

        inline void decrement(size_t i) {
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), i, size);
            if (data8) {
                MEDDLY_DCASSERT(!data16);
                MEDDLY_DCASSERT(!data32);
                MEDDLY_DCASSERT(data8[i]);
                --data8[i];
                return;
            }
            if (data16) {
                MEDDLY_DCASSERT(!data32);
                MEDDLY_DCASSERT(data16[i]);
                if (256 == data16[i]) {
                    MEDDLY_DCASSERT(counts_09bit);
                    --counts_09bit;
                }
                --data16[i];
                return;
            }
            MEDDLY_DCASSERT(data32);
            MEDDLY_DCASSERT(data32[i]);
            if (256 == data32[i]) {
                MEDDLY_DCASSERT(counts_09bit);
                --counts_09bit;
            }
            if (65536 == data32[i]) {
                MEDDLY_DCASSERT(counts_17bit);
                --counts_17bit;
            }
            --data32[i];
        }

        inline bool isZeroBeforeIncrement(size_t i) {
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), i, size);
            if (data8) {
                MEDDLY_DCASSERT(!data16);
                MEDDLY_DCASSERT(!data32);
                if (0==data8[i]) {
                    data8[i] = 1;
                    return true;
                }
                if (0 == ++data8[i]) expand8to16(i);
                return false;
            }
            if (data16) {
                MEDDLY_DCASSERT(!data32);
                if (0==data16[i]) {
                    data16[i] = 1;
                    return true;
                }
                ++data16[i];
                if (256 == data16[i]) ++counts_09bit;
                if (0 == data16[i]) expand16to32(i);
                return false;
            }
            MEDDLY_DCASSERT(data32);
            if (0==data32[i]) {
                data32[i] = 1;
                return true;
            }
            data32[i]++;
            if (256 == data32[i]) ++counts_09bit;
            if (65536 == data32[i]) ++counts_17bit;
            return false;
        }

        inline bool isPositiveAfterDecrement(size_t i) {
            CHECK_RANGE(__FILE__, __LINE__, size_t(0), i, size);
            if (data8) {
                MEDDLY_DCASSERT(!data16);
                MEDDLY_DCASSERT(!data32);
                MEDDLY_DCASSERT(data8[i]);
                return --data8[i] > 0;
            }
            if (data16) {
                MEDDLY_DCASSERT(!data32);
                MEDDLY_DCASSERT(data16[i]);
                if (256 == data16[i]) {
                    MEDDLY_DCASSERT(counts_09bit);
                    --counts_09bit;
                }
                return --data16[i] > 0;
            }
            MEDDLY_DCASSERT(data32);
            MEDDLY_DCASSERT(data32[i]);
            if (256 == data32[i]) {
                MEDDLY_DCASSERT(counts_09bit);
                --counts_09bit;
            }
            if (65536 == data32[i]) {
                MEDDLY_DCASSERT(counts_17bit);
                --counts_17bit;
            }
            return --data32[i] > 0;
        }

    private:

        // Expand from 8-bit to 16-bit entries because of element i
        void expand8to16(size_t i);
        // Expand from 16-bit to 32-bit entries because of element i
        void expand16to32(size_t i);

        void shrink16to8(size_t ns);
        void shrink32to16(size_t ns);
        void shrink32to8(size_t ns);
};


#endif

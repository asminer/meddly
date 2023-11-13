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
};



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


/*
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
            return size_t(bytes) * 8;
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



#endif

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

#ifndef MEDDLY_MEMSTATS_H
#define MEDDLY_MEMSTATS_H

#include "defines.h"

#include <cstddef>  // for size_t

namespace MEDDLY {
    class trackpeak;
    class memstats;
};

/**
    Class to track peak size
*/
class MEDDLY::trackpeak {
    public:
        inline void operator += (size_t s) {
            current += s;
            UPDATEMAX(peak, current);
        }
        inline void operator -= (size_t s) {
            ASSERT(__FILE__, __LINE__, current >= s);
            current -= s;
        }
        inline void operator -= (const trackpeak &s) {
            ASSERT(__FILE__, __LINE__, current >= s.current);
            current -= s.current;
        }
        inline void zero() {
            current = 0;
        }
        inline void init() {
            current = 0;
            peak = 0;
        }

        inline size_t getCurrent()  const { return current; }
        inline size_t getPeak()     const { return peak; }
    private:
        size_t current;
        size_t peak;
};

/**
    Class for memory statistics.
*/
class MEDDLY::memstats {
    public:
        memstats();

        inline void incMemUsed(size_t b) {
            memused += b;
            global_memused += b;
        }
        inline void decMemUsed(size_t b) {
            memused -= b;
            global_memused -= b;
        }

        inline void incMemAlloc(size_t b) {
            memalloc += b;
            global_memalloc += b;
        }

        inline void decMemAlloc(size_t b) {
            memalloc -= b;
            global_memalloc -= b;
        }

        inline void zeroMemUsed() {
            global_memused -= memused;
            memused.zero();
        }
        inline void zeroMemAlloc() {
            global_memalloc -= memalloc;
            memalloc.zero();
        }

        inline size_t getMemUsed() const {
            return memused.getCurrent();
        }
        inline size_t getMemAlloc() const {
            return memalloc.getCurrent();
        }
        inline size_t getPeakMemUsed() const {
            return memused.getPeak();
        }
        inline size_t getPeakMemAlloc() const {
            return memalloc.getPeak();
        }

        static inline size_t getGlobalMemUsed() {
            return global_memused.getCurrent();
        }
        static inline size_t getGlobalMemAlloc() {
            return global_memalloc.getCurrent();
        }
        static inline size_t getGlobalPeakMemUsed() {
            return global_memused.getPeak();
        }
        static inline size_t getGlobalPeakMemAlloc() {
            return global_memalloc.getPeak();
        }

        static void initGlobalStats();

        static void resetGlobalMemory();

    private:
        trackpeak memused;
        trackpeak memalloc;

        static trackpeak global_memused;
        static trackpeak global_memalloc;
};


#endif

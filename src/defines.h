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

#ifndef MEDDLY_DEFINES_H
#define MEDDLY_DEFINES_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// #define ALLOW_DEPRECATED_0_17_2  "11/11/2023"
// #define ALLOW_DEPRECATED_0_17_3  "12/03/2023"
// #define ALLOW_DEPRECATED_0_17_4  "12/19/2023"
#define ALLOW_DEPRECATED_0_17_5

// #define DEBUG_CLEANUP
// #define DEBUG_SLOW
// #define TRACE_ALL_OPS
// #define DEBUG_READ_DD

#include <limits>

// Handy Constants

namespace MEDDLY {
  // const int INF = std::numeric_limits<int>::max();
  // const float NAN = std::numeric_limits<float>::quiet_NaN();
  template<typename T>
  inline T Inf()    { return std::numeric_limits<T>::max(); }
  inline float Nan()  { return std::numeric_limits<float>::quiet_NaN(); }
  inline bool isNan(float t) { return t != t; }
  inline bool isNan(int t) { return false; }

  // Handy Macros

    /// Standard MAX "macro".
    template <class T> inline T MAX(T X,T Y) { return ((X>Y)?X:Y); }
    /// Standard MIN "macro".
    template <class T> inline T MIN(T X,T Y) { return ((X<Y)?X:Y); }
    /// Standard ABS "macro".
    template <class T> inline T ABS(T X) { return ((X<0)?(-X):(X)); }
    /// SWAP "macro".
    template <class T> inline void SWAP(T &x, T &y) { T tmp=x; x=y; y=tmp; }
    /// POSITIVE "macro".
    template <class T> inline bool POSITIVE(T X) { return (X>0) ? true : false; }

    /// Update a maximum
    template <class T> inline void UPDATEMAX(T &X, T Y) { if (Y > X) X = Y; }

    // Number of digits
    template <class T>
    inline unsigned digits(T a) {
        unsigned d;
        for (d=1; a; d++) { a /= 10; }
        return d;
    }

  /// Get the "top level" of an operation.  Works for primed & unprimed.
  inline int topLevel(int k1, int k2) {
    if (ABS(k1) > ABS(k2)) return k1;
    if (ABS(k2) > ABS(k1)) return k2;
    return MAX(k1, k2);
  }

  /// Determine if level k1 is above k2.  Works for primed & unprimed.
  inline bool isLevelAbove(int k1, int k2) {
    if (ABS(k1) > ABS(k2)) return true;
    if (ABS(k2) > ABS(k1)) return false;
    return k1 > k2;
  }

}

/*

   There are now two modes of code generation:
   "DEVELOPMENT_CODE" and "RELEASE_CODE".

   If "DEVELOPMENT_CODE" is defined (usually done in the makefile) then
   debugging macros and assertions will be turned on.  Otherwise we assume
   that we have "RELEASE_CODE" and they are turned off.

   Macros useful for debugging "development code" that are turned off
   for release code (for speed):

   MEDDLY_DCASSERT()
   MEDDLY_CHECK_RANGE(low, x, high+1)

*/

#ifdef DEBUG_PRINTS_ON
#define DEBUG_HASH_H
#define DEBUG_HASH_EXPAND_H
#define DEBUG_MDD_HASH_H
#define DEBUG_MDD_HASH_EXPAND_H
#define TRACE_REDUCE
#define MEMORY_TRACE
#endif

// Safe typecasting for development code;  fast casting otherwise

#ifdef DEVELOPMENT_CODE
#define smart_cast	dynamic_cast
#else
#define smart_cast	static_cast
#endif


// Flags for development version only. Significant reduction in performance.
#ifdef DEVELOPMENT_CODE
#define RANGE_CHECK_ON
#define DCASSERTS_ON
#endif

// #define TRACK_DELETIONS
// #define TRACK_CACHECOUNT
// #define TRACK_UNREACHABLE_NODES


// Use this for assertions that will fail only when your
// code is wrong.  Handy for debugging.
#ifdef DCASSERTS_ON
#include <cassert>
#define MEDDLY_DCASSERT(X) assert(X)
#else
#define MEDDLY_DCASSERT(X)
#endif

// Use this for range checking assertions that should succeed.
#ifdef RANGE_CHECK_ON
#include <cassert>
#include <iostream>
#endif
namespace MEDDLY {
    template <class INT>
    inline void CHECK_RANGE(const char* fn, unsigned ln,
            INT min, INT value, INT max)
    {
#ifdef RANGE_CHECK_ON
        if (value < min || value >= max ) {
            std::cerr << "Check range at " << fn << " line " << ln;
            std::cerr << " failed:\n    min: " << min;
            std::cerr << "\n    val: " << value;
            std::cerr << "\n    max: " << max << '\n';
            assert(false);
        }
#endif
    }
}


//
// Typedefs and constants
//

namespace MEDDLY {

    /** Handles for nodes.
        This should be either int or long, and effectively limits
        the number of possible nodes per forest.
        As an int, we get 2^32-1 possible nodes per forest,
        which should be enough for most applications.
        As a long on a 64-bit machine, we get 2^64-1 possible nodes
        per forest, at the expense of nearly doubling the memory used.
        This also specifies the incoming count range for each node.
    */
    typedef int  node_handle;

    /** Handles for relation nodes.
        TBD: can we just use node_handle everywhere?
     */
    typedef int  rel_node_handle;

    /** Node addresses.
        This is used for internal storage of a node,
        and should probably not be changed.
        The typedef is given simply to clarify the code
        (hopefully :^)
    */
    typedef unsigned long node_address;

};


#endif // #include guard



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



/** @name defines.h
    @type File
    @args \ 

  The base of all files.  So if you change this, everything gets to recompile.

  This file is for good global defines, such as ASSERT and TRACE and crud.

  Since this file is only intended for global definitions, there is no
  associated defines.c or defines.cc file.
 */

//@{

#ifndef DEFINES_H
#define DEFINES_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// #define DEBUG_SLOW
// #define TRACE_ALL_OPS

// Things that everyone will need
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cstdio>
#include <limits>

// Meddly
#include "meddly.h"
#include "meddly_expert.h"

// Macro to handle extern "C"
#ifdef __cplusplus
#  define BEGIN_C_DECLS extern "C" {
#  define END_C_DECLS   }
#else /* !__cplusplus */
#  define BEGIN_C_DECLS  
#  define END_C_DECLS  
#endif /* __cplusplus */

// Handy Constants

namespace MEDDLY {
  const int INF = std::numeric_limits<int>::max();
  const float NAN = std::numeric_limits<float>::quiet_NaN();
  inline bool isNan(float t) { return t != t; }
  inline bool isNan(int t) { return t != t; }

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

  // maxUlps would be a small +ve integer.
  inline bool isAlmostEqual(int A, int B) {
    static const int maxUlps = 16;
  
    if (isNan(A) && isNan(B)) return true;
  
    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.
    MEDDLY_CHECK_RANGE(1, maxUlps, 4 * 1024 * 1024);
  
    // Make aInt lexicographically ordered as a twos-complement int
    // int aInt = *(int*)&A;
    // if (aInt < 0) aInt = 0x80000000 - aInt;
  
    // Make bInt lexicographically ordered as a twos-complement int
    // int bInt = *(int*)&B;
    // if (bInt < 0) bInt = 0x80000000 - bInt;
  
    // return ABS(aInt - bInt) <= maxUlps;
  
    // return A < 0
    //         ? B < 0 ? ABS(B - A) <= maxUlps: ABS(mask - A - B) <= maxUlps
    //         : B < 0 ? ABS(A - mask + B) <= maxUlps: ABS(A - B) <= maxUlps;
#if 0
    static const int mask = 0x80000000;
    return ABS(A < 0? B < 0? B - A: mask - A - B: B < 0? A - mask + B: A - B)
            <= maxUlps;
#else
    if (A < 0) A = 0x80000000 - A;
    if (B < 0) B = 0x80000000 - B;
    return ABS(A - B) <= maxUlps;
#endif
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



#endif

//@}


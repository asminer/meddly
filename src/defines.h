
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

  // Number of digits
  inline int digits(long a) {
    int d;
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



#endif

//@}


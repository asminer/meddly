
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

// Things that everyone will need
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cstdio>

// Macro to handle extern "C"
#ifdef __cplusplus
#  define BEGIN_C_DECLS extern "C" {
#  define END_C_DECLS   }
#else /* !__cplusplus */
#  define BEGIN_C_DECLS  
#  define END_C_DECLS  
#endif /* __cplusplus */

// Handy Constants

const int INF = INT_MAX;

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
template <class T> inline bool POSITIVE(T X) { return (X>0) ? 1 : 0; }

/// SIGN "macro".
template <class T> inline int SIGN(T X) { return (X<0) ? -1 : POSITIVE(X); }

/*

   There are now two modes of code generation:
   "DEVELOPMENT_CODE" and "RELEASE_CODE".

   If "DEVELOPMENT_CODE" is defined (usually done in the makefile) then 
   debugging macros and assertions will be turned on.  Otherwise we assume
   that we have "RELEASE_CODE" and they are turned off.

   Macros useful for debugging "development code" that are turned off 
   for release code (for speed):

   DCASSERT()
   CHECK_RANGE(low, x, high+1)

*/

#ifdef DEVELOPMENT_CODE
#if 0
#define DEBUG_HASH_H
#define DEBUG_HASH_EXPAND_H
#define DEBUG_MDD_HASH_H
#define DEBUG_MDD_HASH_EXPAND_H
#define DEBUG_MDD_H
#define DEBUG_MDD_SET
#define TRACK_DELETIONS
#define TRACE_REDUCE
#define MEMORY_TRACE
#endif
#endif

#ifdef DEVELOPMENT_CODE
#define MEM_TRACE_ON
#define RANGE_CHECK_ON
#define DCASSERTS_ON
#endif

// Use this for assertions that you always check for
#define ASSERT(X) assert(X)

// Use this for assertions that will fail only when your
// code is wrong.  Handy for debugging.
#ifdef DCASSERTS_ON 
#define DCASSERT(X) assert(X)
#else
#define DCASSERT(X)
#endif

// Also useful for debugging.
#ifdef RANGE_CHECK_ON
#define CHECK_RANGE(MIN, VALUE, MAX)  {ASSERT(VALUE<MAX);ASSERT(VALUE>=MIN);}
#else
#define CHECK_RANGE(MIN, VALUE, MAX)
#endif

// Safe typecasting for development code;  fast casting otherwise

#ifdef DEVELOPMENT_CODE
#define smart_cast	dynamic_cast
#else
#define smart_cast	static_cast
#endif

#endif

//@}


# $Id$

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

## AC_PREREQ(2.60)
AC_PREREQ(2.59)
AC_INIT([meddly], [0.2], [junaid@cs.iastate.edu])
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([include/meddly.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_HEADER_TIME
## AC_TYPE_UINT32_T

# Checks for library functions.
AC_FUNC_FORK
AC_HEADER_STDC
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_TYPE_SIGNAL
AC_FUNC_STRTOD
AC_CHECK_FUNCS([gettimeofday memmove memset strtol])

AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile])
AC_OUTPUT

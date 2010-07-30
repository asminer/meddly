## $Id$

## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

check_PROGRAMS = \
  nqueens check_xA

TESTS = \
  nqueens check_xA

AM_CPPFLAGS = -I../src -I../src/operations
AM_CXXFLAGS = -Wall


nqueens_SOURCES = nqueens.cc
nqueens_LDADD = ../src/libmeddly.la

check_xA_SOURCES = check_xA.cc
check_xA_LDADD = ../src/libmeddly.la


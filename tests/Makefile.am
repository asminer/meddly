## $Id$

## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

check_PROGRAMS = \
  nqueens check_xA chk_cross kanban

TESTS = \
  nqueens check_xA chk_cross kanban

AM_CPPFLAGS = -I../src -I../src/operations
AM_CXXFLAGS = -Wall


nqueens_SOURCES = nqueens.cc
nqueens_LDADD = ../src/libmeddly.la

check_xA_SOURCES = check_xA.cc
check_xA_LDADD = ../src/libmeddly.la

chk_cross_SOURCES = chk_cross.cc
chk_cross_LDADD = ../src/libmeddly.la

kanban_SOURCES = kanban.cc simple_model.h simple_model.cc
kanban_LDADD = ../src/libmeddly.la

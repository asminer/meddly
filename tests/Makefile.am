## $Id$

## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

check_PROGRAMS = \
  sat_test nqueens check_xA chk_cross kanban

TESTS = \
  sat_test nqueens check_xA chk_cross kanban

AM_CPPFLAGS = -I../src -I../src/operations
AM_CXXFLAGS = -Wall

sat_test_SOURCES = sat_test.cc
sat_test_LDADD = ../src/libmeddly.la

nqueens_SOURCES = nqueens.cc
nqueens_LDADD = ../src/libmeddly.la

check_xA_SOURCES = check_xA.cc
check_xA_LDADD = ../src/libmeddly.la

chk_cross_SOURCES = chk_cross.cc
chk_cross_LDADD = ../src/libmeddly.la

kanban_SOURCES = kanban.cc simple_model.h simple_model.cc
kanban_LDADD = ../src/libmeddly.la

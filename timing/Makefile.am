
## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

noinst_PROGRAMS = \
    ddedge \
    mark \
    reduce

AM_CXXFLAGS = -Wall

ddedge_SOURCES = timer.h ddedge.cc reporting.h reporting.cc
ddedge_LDADD = ../src/libmeddly.la

mark_SOURCES = timer.h mark.cc reporting.h reporting.cc simple_model.h simple_model.cc
mark_LDADD = ../src/libmeddly.la

reduce_SOURCES = timer.h reduce.cc reporting.h reporting.cc
reduce_LDADD = ../src/libmeddly.la


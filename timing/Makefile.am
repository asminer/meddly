
## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

noinst_PROGRAMS = \
    mark \
    reduce

AM_CXXFLAGS = -Wall

mark_SOURCES = timer.h mark.cc simple_model.h simple_model.cc
mark_LDADD = ../src/libmeddly.la

reduce_SOURCES = timer.h reduce.cc
reduce_LDADD = ../src/libmeddly.la


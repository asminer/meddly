## $Id$

## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

SUBDIRS = forests operations storage .

lib_LTLIBRARIES = libmeddly.la

_revision:
	./BuildRevNo

revision.h:
	./BuildRevNo

BUILT_SOURCES = _revision revision.h 
CLEANFILES = revision.h

libmeddly_la_SOURCES = \
  defines.h \
  meddly.h \
  meddly_expert.h \
  meddly.cc \
  node_wrappers.cc \
  ops.cc \
  dd_edge.cc \
  enumerator.cc \
  domain.cc \
  forest.cc \
  compute_table.h \
  compute_table.cc \
  unique_table.h \
  unique_table.cc \
  hash_stream.h \
  timer.h 

libmeddly_la_LIBADD = \
  forests/libforests.la \
  operations/libops.la \
  storage/libstorage.la

## compile apps

AM_CXXFLAGS = -Wall

include_HEADERS = meddly.h meddly_expert.h 

EXTRA_DIST = BuildRevNo revision.h

## $Id$

## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

SUBDIRS = forests operations storage reordering .

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
  meddly.hh \
  meddly_expert.h \
  meddly_expert.hh \
  meddly.cc \
  io.cc \
  node_wrappers.cc \
  ops.cc \
  dd_edge.cc \
  enumerator.cc \
  domain.cc \
  forest.cc \
  compute_table.cc \
  unique_table.h \
  unique_table.cc \
  hash_stream.h \
  loggers.h \
  loggers.cc \
  timer.h \
  global_rebuilder.cc

libmeddly_la_LIBADD = \
  forests/libforests.la \
  operations/libops.la \
  storage/libstorage.la \
  reordering/libreordering.la

## compile apps

AM_CXXFLAGS = -Wall

include_HEADERS = meddly.h meddly_expert.h meddly.hh meddly_expert.hh

EXTRA_DIST = BuildRevNo revision.h

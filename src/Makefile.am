## $Id$

## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

SUBDIRS = forests operations .

lib_LTLIBRARIES = libmeddly.la

libmeddly_la_SOURCES = \
  revision.h \
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
  unique_table.cc
	defines.h \
  hash_stream.h \
	timer.h 

libmeddly_la_LIBADD = \
	forests/libforests.la \
	operations/libops.la

## compile apps

AM_CXXFLAGS = -Wall

include_HEADERS = meddly.h meddly_expert.h 


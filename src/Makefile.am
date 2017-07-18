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
	dd_edge.cc \
	temp_dd_edge.cc \
	domain.cc \
	forest.cc \
	compute_manager.cc \
	compute_cache.cc \
	compute_cache.h \
	defines.h \
	hash.h \
	fixed_size_hash.h \
	chained_hash.h \
	mdd_hash.h \
	timer.h 

libmeddly_la_LIBADD = \
	forests/libforests.la \
	operations/libops.la

## compile apps

AM_CXXFLAGS = -Wall

include_HEADERS = meddly.h meddly_expert.h 


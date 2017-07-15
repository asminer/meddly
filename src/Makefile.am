## $Id$

## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

##
## This version: don't use recursive make
##

## SUBDIRS = forests operations storage .

lib_LTLIBRARIES = libmeddly.la

_revision:
	./BuildRevNo

revision.h:
	./BuildRevNo

BUILT_SOURCES = _revision revision.h 
CLEANFILES = revision.h

libmeddly_la_SOURCES = \
  meddly.cc \
  io.cc \
  node_wrappers.cc \
  ops.cc \
  dd_edge.cc \
  enumerator.cc \
  domain.cc \
  forest.cc \
  compute_table.cc \
  unique_table.cc \
  loggers.cc \
  global_rebuilder.cc \
  \
  forests/mt.cc \
  forests/mtmdd.cc \
  forests/mtmddbool.cc \
  forests/mtmddint.cc \
  forests/mtmddreal.cc \
  forests/mtmxd.cc \
  forests/mtmxdbool.cc \
  forests/mtmxdint.cc \
  forests/mtmxdreal.cc \
  forests/ev.cc \
  forests/evmdd.cc \
  forests/evmdd_plusint.cc \
  forests/evmxd.cc \
  forests/evmxd_timesreal.cc \
  forests/init_forests.cc \
  \
  operations/apply_base.cc   \
  operations/mpz_object.cc   \
  operations/copy.cc         \
  operations/cardinality.cc  \
  operations/complement.cc   \
  operations/maxmin_range.cc \
  operations/mdd2index.cc    \
  operations/union.cc        \
  operations/intersection.cc \
  operations/difference.cc   \
  operations/cross.cc        \
  operations/maxmin.cc       \
  operations/plus.cc         \
  operations/minus.cc        \
  operations/multiply.cc     \
  operations/divide.cc       \
  operations/modulo.cc       \
  operations/comp_eq.cc      \
  operations/comp_ne.cc      \
  operations/comp_lt.cc      \
  operations/comp_le.cc      \
  operations/comp_gt.cc      \
  operations/comp_ge.cc      \
  operations/prepostimage.cc \
  operations/reach_bfs.cc    \
  operations/reach_dfs.cc    \
  operations/sat_pregen.cc   \
  operations/sat_otf.cc      \
  operations/sat_impl.cc   \
  operations/vect_matr.cc    \
  operations/mm_mult.cc      \
  operations/init_builtin.cc \
  \
  storage/holeman.cc \
  storage/hm_array.cc \
  storage/hm_grid.cc \
  storage/hm_heap.cc \
  storage/hm_none.cc \
  storage/old_scheme.cc \
  storage/simple.cc \
  storage/compact.cc \
  storage/ct_classic.cc \
  storage/init_storage.cc

  
## compile apps

include_HEADERS = meddly.h meddly_expert.h meddly.hh meddly_expert.hh

EXTRA_DIST = BuildRevNo revision.h

AM_CXXFLAGS = -Wall

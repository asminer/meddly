## $Id$

## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

## bin_PROGRAMS = \
noinst_PROGRAMS = \
testmddlib \
graphs_to_mdds \
test_union_mdd \
test_mtmdd \
test_mtmxd \
test_evmdd \
build_expression \
kan1_batch \
dining_phils \
rubiks_cube \
nqueens \
queen_cover

## CXXFLAGS = 
CXXFLAGS += -Wall
## CXXFLAGS += -g
# CXXFLAGS += -ggdb
## CXXFLAGS += -pg
# CXXFLAGS += -DDEVELOPMENT_CODE
## CXXFLAGS += -DDIFF_SAT
## CXXFLAGS += -O2
## CXXFLAGS += -O3
## CXXFLAGS += -m64

## testmddlib_SOURCES = `$(CXX) testmddlib.cc -MM`
testmddlib_SOURCES = testmddlib.cc
testmddlib_LDADD = ../src/libmeddly.la

graphs_to_mdds_SOURCES = graphs_to_mdds.cc
graphs_to_mdds_LDADD = ../src/libmeddly.la

test_union_mdd_SOURCES = test_union_mdd.cc
test_union_mdd_LDADD = ../src/libmeddly.la

test_mtmdd_SOURCES = test_mtmdd.cc
test_mtmdd_LDADD = ../src/libmeddly.la

test_mtmxd_SOURCES = test_mtmxd.cc
test_mtmxd_LDADD = ../src/libmeddly.la

build_expression_SOURCES = build_expression.cc
build_expression_LDADD = ../src/libmeddly.la

kan1_batch_SOURCES = kan1_batch.cc
kan1_batch_LDADD = ../src/libmeddly.la

dining_phils_SOURCES = dining_phils.cc
dining_phils_LDADD = ../src/libmeddly.la

rubiks_cube_SOURCES = rubiks_cube.cc
rubiks_cube_LDADD = ../src/libmeddly.la

test_evmdd_SOURCES = test_evmdd.cc
test_evmdd_LDADD = ../src/libmeddly.la

nqueens_SOURCES = nqueens.cc
nqueens_LDADD = ../src/libmeddly.la

queen_cover_SOURCES = queen_cover.cc
queen_cover_LDADD = ../src/libmeddly.la


## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

noinst_PROGRAMS = \
    enctest \
    \
    bug_00 \
    bug_01 \
    bug_02 \
    check_xA \
    chk_copy \
    chk_cross \
    chk_ct \
    chk_ct_old \
    chk_iters \
    chk_minterms \
    chk_reorder \
    ddedge \
    domain_list \
    kan_batch \
    kan_index \
    kan_io \
    kan_show \
    kanban \
    nodestor \
    nqueens \
    ops_arith \
    ops_comp \
    ops_set \
    refcount \
    sat_test \
    randtest

TESTS = \
    bug_00 \
    bug_01 \
    bug_02 \
    check_xA \
    chk_copy \
    chk_cross \
    chk_ct \
    chk_ct_old \
    chk_iters \
    chk_minterms \
    chk_reorder \
    ddedge \
    domain_list \
    kan_batch \
    kan_index \
    kan_io \
    kan_show \
    kanban \
    nodestor \
    nqueens \
    ops_arith \
    ops_comp \
    ops_set \
    refcount \
    sat_test

AM_CXXFLAGS = -Wall

#
# non-test programs
#

enctest_SOURCES = enctest.cc
enctest_LDADD =

#
# test programs
#

bug_00_SOURCES = bug_00.cc simple_model.h simple_model.cc
bug_00_LDADD = ../src/libmeddly.la

bug_01_SOURCES = bug_01.cc
bug_01_LDADD = ../src/libmeddly.la

bug_02_SOURCES = bug_02.cc
bug_02_LDADD = ../src/libmeddly.la

check_xA_SOURCES = check_xA.cc
check_xA_LDADD = ../src/libmeddly.la

chk_copy_SOURCES = chk_copy.cc
chk_copy_LDADD = ../src/libmeddly.la

chk_cross_SOURCES = chk_cross.cc
chk_cross_LDADD = ../src/libmeddly.la

chk_ct_SOURCES = chk_ct.cc
chk_ct_LDADD = ../src/libmeddly.la

chk_ct_old_SOURCES = chk_ct_old.cc
chk_ct_old_LDADD = ../src/libmeddly.la

chk_iters_SOURCES = chk_iters.cc randomize.cc
chk_iters_LDADD = ../src/libmeddly.la

chk_minterms_SOURCES = chk_minterms.cc randomize.cc
chk_minterms_LDADD = ../src/libmeddly.la

chk_reorder_SOURCES = chk_reorder.cc
chk_reorder_LDADD = ../src/libmeddly.la

ddedge_SOURCES = ddedge.cc
ddedge_LDADD = ../src/libmeddly.la

domain_list_SOURCES = domain_list.cc
domain_list_LDADD = ../src/libmeddly.la

kan_batch_SOURCES = kan_batch.cc simple_model.h simple_model.cc \
                  kan_rs1.h kan_rs2.h kan_rs3.h
kan_batch_LDADD = ../src/libmeddly.la

kan_index_SOURCES = kan_index.cc simple_model.h simple_model.cc
kan_index_LDADD = ../src/libmeddly.la

kan_io_SOURCES = kan_io.cc simple_model.h simple_model.cc
kan_io_LDADD = ../src/libmeddly.la

kan_show_SOURCES = kan_show.cc simple_model.h simple_model.cc \
                  kan_rs1.h kan_rs2.h kan_rs3.h
kan_show_LDADD = ../src/libmeddly.la

kanban_SOURCES = kanban.cc simple_model.h simple_model.cc
kanban_LDADD = ../src/libmeddly.la

nodestor_SOURCES = nodestor.cc
nodestor_LDADD = ../src/libmeddly.la

nqueens_SOURCES = nqueens.cc
nqueens_LDADD = ../src/libmeddly.la

ops_arith_SOURCES = ops_arith.cc randomize.cc
ops_arith_LDADD = ../src/libmeddly.la

ops_comp_SOURCES = ops_comp.cc randomize.cc
ops_comp_LDADD = ../src/libmeddly.la

ops_set_SOURCES = ops_set.cc randomize.cc
ops_set_LDADD = ../src/libmeddly.la

refcount_SOURCES = refcount.cc
refcount_LDADD = ../src/libmeddly.la

sat_test_SOURCES = sat_test.cc
sat_test_LDADD = ../src/libmeddly.la

randtest_SOURCES = randomize.cc randtest.cc
randtest_LDADD = ../src/libmeddly.la

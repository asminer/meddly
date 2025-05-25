
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
    chk_mtcoll \
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
    ops_comp \
    ops_maxmin \
    ops_plusminus \
    ops_mult \
    ops_divmod \
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
    \
    chk_minterms.--set.--fully._sh \
    chk_minterms.--set.--quasi._sh \
    chk_minterms.--rel.--fully._sh \
    chk_minterms.--rel.--quasi._sh \
    chk_minterms.--rel.--ident._sh \
    \
    chk_mtcoll.--set.--fully._sh \
    chk_mtcoll.--set.--quasi._sh \
    chk_mtcoll.--rel.--fully._sh \
    chk_mtcoll.--rel.--quasi._sh \
    chk_mtcoll.--rel.--ident._sh \
    \
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
    \
    ops_comp.--set.--MT.--int._sh \
    ops_comp.--set.--MT.--float._sh \
    ops_comp.--set.--EVp.--long._sh \
    ops_comp.--rel.--MT.--int._sh \
    ops_comp.--rel.--MT.--float._sh \
    ops_comp.--rel.--EVp.--long._sh \
    ops_comp.--rel.--EVs.--float._sh \
    \
    ops_maxmin.--set.--MT.--int._sh \
    ops_maxmin.--set.--MT.--float._sh \
    ops_maxmin.--set.--EVp.--long._sh \
    ops_maxmin.--rel.--MT.--int._sh \
    ops_maxmin.--rel.--MT.--float._sh \
    ops_maxmin.--rel.--EVp.--long._sh \
    ops_maxmin.--rel.--EVs.--float._sh \
    \
    ops_plusminus.--set.--MT.--int._sh \
    ops_plusminus.--set.--MT.--float._sh \
    ops_plusminus.--set.--EVp.--long._sh \
    ops_plusminus.--rel.--MT.--int._sh \
    ops_plusminus.--rel.--MT.--float._sh \
    ops_plusminus.--rel.--EVp.--long._sh \
    ops_plusminus.--rel.--EVs.--float._sh \
    \
    ops_mult.--set.--MT.--int._sh \
    ops_mult.--set.--MT.--float._sh \
    ops_mult.--set.--EVp.--long._sh \
    ops_mult.--rel.--MT.--int._sh \
    ops_mult.--rel.--MT.--float._sh \
    ops_mult.--rel.--EVp.--long._sh \
    ops_mult.--rel.--EVs.--float._sh \
    \
    ops_divmod.--set.--MT.--int._sh \
    ops_divmod.--set.--MT.--float._sh \
    ops_divmod.--set.--EVp.--long._sh \
    ops_divmod.--rel.--MT.--int._sh \
    ops_divmod.--rel.--MT.--float._sh \
    ops_divmod.--rel.--EVp.--long._sh \
    ops_divmod.--rel.--EVs.--float._sh \
    \
    ops_set.--set._sh \
    ops_set.--rel._sh \
    \
    refcount \
    sat_test

AM_CXXFLAGS = -Wall

EXTRA_DIST = buildWrapper.sh

#
# Automatically build a wrapper script,
# for tests with command line arguments.
#
###########################################

chk_minterms.--set.--fully._sh: buildWrapper.sh
	./buildWrapper.sh $@

chk_minterms.--set.--quasi._sh: buildWrapper.sh
	./buildWrapper.sh $@

chk_minterms.--rel.--fully._sh: buildWrapper.sh
	./buildWrapper.sh $@

chk_minterms.--rel.--quasi._sh: buildWrapper.sh
	./buildWrapper.sh $@

chk_minterms.--rel.--ident._sh: buildWrapper.sh
	./buildWrapper.sh $@

###########################################

chk_mtcoll.--set.--fully._sh: buildWrapper.sh
	./buildWrapper.sh $@

chk_mtcoll.--set.--quasi._sh: buildWrapper.sh
	./buildWrapper.sh $@

chk_mtcoll.--rel.--fully._sh: buildWrapper.sh
	./buildWrapper.sh $@

chk_mtcoll.--rel.--quasi._sh: buildWrapper.sh
	./buildWrapper.sh $@

chk_mtcoll.--rel.--ident._sh: buildWrapper.sh
	./buildWrapper.sh $@

###########################################

ops_comp.--set.--MT.--int._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_comp.--set.--MT.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_comp.--set.--EVp.--long._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_comp.--rel.--MT.--int._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_comp.--rel.--MT.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_comp.--rel.--EVp.--long._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_comp.--rel.--EVs.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

###########################################

ops_maxmin.--set.--MT.--int._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_maxmin.--set.--MT.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_maxmin.--set.--EVp.--long._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_maxmin.--rel.--MT.--int._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_maxmin.--rel.--MT.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_maxmin.--rel.--EVp.--long._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_maxmin.--rel.--EVs.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

###########################################

ops_plusminus.--set.--MT.--int._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_plusminus.--set.--MT.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_plusminus.--set.--EVp.--long._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_plusminus.--rel.--MT.--int._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_plusminus.--rel.--MT.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_plusminus.--rel.--EVp.--long._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_plusminus.--rel.--EVs.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

###########################################

ops_mult.--set.--MT.--int._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_mult.--set.--MT.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_mult.--set.--EVp.--long._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_mult.--rel.--MT.--int._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_mult.--rel.--MT.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_mult.--rel.--EVp.--long._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_mult.--rel.--EVs.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

###########################################

ops_divmod.--set.--MT.--int._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_divmod.--set.--MT.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_divmod.--set.--EVp.--long._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_divmod.--rel.--MT.--int._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_divmod.--rel.--MT.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_divmod.--rel.--EVp.--long._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_divmod.--rel.--EVs.--float._sh: buildWrapper.sh
	./buildWrapper.sh $@

###########################################

ops_set.--set._sh: buildWrapper.sh
	./buildWrapper.sh $@

ops_set.--rel._sh: buildWrapper.sh
	./buildWrapper.sh $@

###########################################

clean-local:
	rm -f *._sh

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

chk_cross_SOURCES = chk_cross.cc randomize.cc
chk_cross_LDADD = ../src/libmeddly.la

chk_ct_SOURCES = chk_ct.cc
chk_ct_LDADD = ../src/libmeddly.la

chk_ct_old_SOURCES = chk_ct_old.cc
chk_ct_old_LDADD = ../src/libmeddly.la

chk_iters_SOURCES = chk_iters.cc randomize.cc
chk_iters_LDADD = ../src/libmeddly.la

chk_minterms_SOURCES = chk_minterms.cc randomize.h randomize.cc
chk_minterms_LDADD = ../src/libmeddly.la

chk_mtcoll_SOURCES = chk_mtcoll.cc randomize.h randomize.cc
chk_mtcoll_LDADD = ../src/libmeddly.la

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

ops_comp_SOURCES = ops_comp.cc randomize.h randomize.cc
ops_comp_LDADD = ../src/libmeddly.la

ops_maxmin_SOURCES = ops_maxmin.cc randomize.h randomize.cc
ops_maxmin_LDADD = ../src/libmeddly.la

ops_plusminus_SOURCES = ops_plusminus.cc randomize.h randomize.cc
ops_plusminus_LDADD = ../src/libmeddly.la

ops_mult_SOURCES = ops_mult.cc randomize.h randomize.cc
ops_mult_LDADD = ../src/libmeddly.la

ops_divmod_SOURCES = ops_divmod.cc randomize.h randomize.cc
ops_divmod_LDADD = ../src/libmeddly.la

ops_set_SOURCES = ops_set.cc randomize.h randomize.cc
ops_set_LDADD = ../src/libmeddly.la

refcount_SOURCES = refcount.cc
refcount_LDADD = ../src/libmeddly.la

sat_test_SOURCES = sat_test.cc
sat_test_LDADD = ../src/libmeddly.la

randtest_SOURCES = randomize.h randomize.cc randtest.cc
randtest_LDADD = ../src/libmeddly.la


## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

##
## This version: don't use recursive make
##

## SUBDIRS = forests operations storage .

lib_LTLIBRARIES = libmeddly.la

libmeddly_la_SOURCES = \
    compute_table.cc \
    ct_entry_key.cc \
    ct_entry_result.cc \
    ct_entry_type.cc \
    ct_initializer.cc \
    dd_edge.cc \
    domain.cc \
    encoders.cc \
    enumerator.cc \
    error.cc \
    forest.cc \
    global_rebuilder.cc \
    impl_unique_table.cc \
    initializer.cc \
    io.cc \
    loggers.cc \
    memory.cc \
    memstats.cc \
    node_headers.cc \
    node_storage.cc \
    oper.cc \
    oper_binary.cc \
    oper_special.cc \
    oper_unary.cc \
    opname.cc \
    opname_numer.cc \
    opname_satur.cc \
    ops_builtin.cc \
    policies.cc \
    relation_node.cc \
    unique_table.cc \
    unpacked_node.cc \
    variable.cc \
    varorder.cc \
    \
    forests/ev.cc \
    forests/evmdd.cc \
    forests/evmdd_pluslong.cc \
    forests/evmdd_timesreal.cc \
    forests/evmxd.cc \
    forests/evmxd_pluslong.cc \
    forests/evmxd_timesreal.cc \
    forests/init_forests.cc \
    forests/mt.cc \
    forests/mtmdd.cc \
    forests/mtmddbool.cc \
    forests/mtmddint.cc \
    forests/mtmddreal.cc \
    forests/mtmxd.cc \
    forests/mtmxdbool.cc \
    forests/mtmxdint.cc \
    forests/mtmxdreal.cc \
    \
    memory_managers/array_grid.cc  \
    memory_managers/freelists.cc  \
    memory_managers/heap_manager.cc  \
    memory_managers/init_managers.cc  \
    memory_managers/malloc_style.cc  \
    memory_managers/orig_grid.cc  \
    \
    operations/apply_base.cc  \
    operations/cardinality.cc  \
    operations/comp_eq.cc      \
    operations/comp_ge.cc      \
    operations/comp_gt.cc      \
    operations/comp_le.cc      \
    operations/comp_lt.cc      \
    operations/comp_ne.cc      \
    operations/complement.cc   \
    operations/constrained.cc   \
    operations/copy.cc         \
    operations/cross.cc        \
    operations/cycle.cc        \
    operations/difference.cc   \
    operations/divide.cc       \
    operations/intersection.cc \
    operations/maxmin.cc       \
    operations/maxmin_range.cc \
    operations/mdd2index.cc    \
    operations/minus.cc        \
    operations/mm_mult.cc      \
    operations/modulo.cc       \
    operations/mpz_object.cc   \
    operations/multiply.cc     \
    operations/plus.cc         \
    operations/prepostimage.cc \
    operations/prepostplus.cc  \
    operations/reach_bfs.cc    \
    operations/reach_dfs.cc    \
    operations/sat_hyb.cc    \
    operations/sat_impl.cc    \
    operations/sat_otf.cc      \
    operations/sat_pregen.cc   \
    operations/sccgraph.cc     \
    operations/select.cc       \
    operations/transitive_closure.cc \
    operations/union.cc        \
    operations/vect_matr.cc    \
    \
    storage/best.cc \
    storage/ct_styles.cc \
    storage/init_storage.cc \
    storage/pattern.cc \
    storage/simple.cc


## compile apps

include_HEADERS = meddly.h

pkginclude_HEADERS = \
    compute_table.h \
    ct_entry_key.h \
    ct_entry_result.h \
    ct_entry_type.h \
    ct_initializer.h \
    dd_edge.h \
    defines.h \
    domain.h \
    encoders.h \
    enumerator.h \
    error.h \
    forest.h \
    global_rebuilder.h \
    hash_stream.h \
    heap.h \
    initializer.h \
    io.h \
    memory.h \
    memstats.h \
    minterms.h \
    node_headers.h \
    node_storage.h \
    oper.h  \
    oper_binary.h  \
    oper_special.h  \
    oper_unary.h  \
    opname.h  \
    opname_numer.h  \
    opname_satur.h  \
    policies.h \
    relation_node.h \
    revision.h \
    timer.h \
    unique_table.h \
    unpacked_node.h \
    variable.h \
    varorder.h

AM_CXXFLAGS = -Wall

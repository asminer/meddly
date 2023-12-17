
## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

##
## This version: don't use recursive make
##

## SUBDIRS = forests operations storage .

lib_LTLIBRARIES = libmeddly.la

##
## This list should include any source file that should be part of
## the library's tarball
##

libmeddly_la_SOURCES = \
    arrays.h                arrays.cc \
    compute_table.h         compute_table.cc \
    ct_entry_key.h          ct_entry_key.cc \
    ct_entry_result.h       ct_entry_result.cc \
    ct_entry_type.h         ct_entry_type.cc \
    ct_initializer.h        ct_initializer.cc \
    dd_edge.h               dd_edge.cc \
    defines.h               \
    domain.h                domain.cc \
    edge_value.h            edge_value.cc \
    enumerator.h            enumerator.cc \
    error.h                 error.cc \
    forest.h                forest.cc \
    global_rebuilder.h      global_rebuilder.cc \
    hash_stream.h           \
    heap.h                  \
    impl_unique_table.h     impl_unique_table.cc \
    initializer.h           initializer.cc \
    io.h                    io.cc \
    io_dot.h                io_dot.cc \
    io_mdds.h               io_mdds.cc \
    loggers.h               loggers.cc \
    memory.h                memory.cc \
    memstats.h              memstats.cc \
    minterms.h              \
    node_headers.h          node_headers.cc \
    node_marker.h           node_marker.cc \
    node_storage.h          node_storage.cc \
    oper.h                  oper.cc \
    oper_binary.h           oper_binary.cc \
    oper_special.h          oper_special.cc \
    oper_unary.h            oper_unary.cc \
    operators.h             operators.cc \
    opname.h                opname.cc \
    opname_numer.h          opname_numer.cc \
    opname_satur.h          opname_satur.cc \
    ops_builtin.h           ops_builtin.cc \
    policies.h              policies.cc \
    relation_node.h         relation_node.cc \
    statset.h               statset.cc \
    terminal.h              terminal.cc \
    unique_table.h          unique_table.cc \
    unpacked_node.h         unpacked_node.cc \
    variable.h              variable.cc \
    varorder.h              varorder.cc \
    \
    forests/ev.h                    forests/ev.cc \
    forests/evmdd.h                 forests/evmdd.cc \
    forests/evmdd_pluslong.h        forests/evmdd_pluslong.cc \
    forests/evmdd_timesreal.h       forests/evmdd_timesreal.cc \
    forests/evmxd.h                 forests/evmxd.cc \
    forests/evmxd_pluslong.h        forests/evmxd_pluslong.cc \
    forests/evmxd_timesreal.h       forests/evmxd_timesreal.cc \
    forests/init_forests.h          forests/init_forests.cc \
    forests/mt.h                    forests/mt.cc \
    forests/mtmdd.h                 forests/mtmdd.cc \
    forests/mtmddbool.h             forests/mtmddbool.cc \
    forests/mtmddint.h              forests/mtmddint.cc \
    forests/mtmddreal.h             forests/mtmddreal.cc \
    forests/mtmxd.h                 forests/mtmxd.cc \
    forests/mtmxdbool.h             forests/mtmxdbool.cc \
    forests/mtmxdint.h              forests/mtmxdint.cc \
    forests/mtmxdreal.h             forests/mtmxdreal.cc \
    \
    memory_managers/array_grid.h        memory_managers/array_grid.cc  \
    memory_managers/freelists.h         memory_managers/freelists.cc  \
    memory_managers/heap_manager.h      memory_managers/heap_manager.cc  \
    memory_managers/hole_base.h         \
    memory_managers/init_managers.h     memory_managers/init_managers.cc  \
    memory_managers/malloc_style.h      memory_managers/malloc_style.cc  \
    memory_managers/orig_grid.h         memory_managers/orig_grid.cc  \
    \
    operations/apply_base.h             operations/apply_base.cc  \
    operations/cardinality.h            operations/cardinality.cc  \
    operations/comp_eq.h                operations/comp_eq.cc      \
    operations/comp_ge.h                operations/comp_ge.cc      \
    operations/comp_gt.h                operations/comp_gt.cc      \
    operations/comp_le.h                operations/comp_le.cc      \
    operations/comp_lt.h                operations/comp_lt.cc      \
    operations/comp_ne.h                operations/comp_ne.cc      \
    operations/complement.h             operations/complement.cc   \
    operations/constrained.h            operations/constrained.cc   \
    operations/copy.h                   operations/copy.cc         \
    operations/cross.h                  operations/cross.cc        \
    operations/cycle.h                  operations/cycle.cc        \
    operations/difference.h             operations/difference.cc   \
    operations/divide.h                 operations/divide.cc       \
    operations/intersection.h           operations/intersection.cc \
    operations/maxmin.h                 operations/maxmin.cc       \
    operations/maxmin_range.h           operations/maxmin_range.cc \
    operations/mdd2index.h              operations/mdd2index.cc    \
    operations/minus.h                  operations/minus.cc        \
    operations/mm_mult.h                operations/mm_mult.cc      \
    operations/modulo.h                 operations/modulo.cc       \
    operations/mpz_object.h             operations/mpz_object.cc   \
    operations/multiply.h               operations/multiply.cc     \
    operations/plus.h                   operations/plus.cc         \
    operations/prepostimage.h           operations/prepostimage.cc \
    operations/prepostplus.h            operations/prepostplus.cc  \
    operations/reach_bfs.h              operations/reach_bfs.cc    \
    operations/reach_dfs.h              operations/reach_dfs.cc    \
    operations/sat_hyb.h                operations/sat_hyb.cc    \
    operations/sat_impl.h               operations/sat_impl.cc    \
    operations/sat_otf.h                operations/sat_otf.cc      \
    operations/sat_pregen.h             operations/sat_pregen.cc   \
    operations/sccgraph.h               operations/sccgraph.cc     \
    operations/select.h                 operations/select.cc       \
    operations/transitive_closure.h     operations/transitive_closure.cc \
    operations/union.h                  operations/union.cc        \
    operations/vect_matr.h              operations/vect_matr.cc    \
    \
    reordering/bring_up_reordering.h            \
    reordering/highest_inversion_reordering.h   \
    reordering/larc_reordering.h                \
    reordering/lowest_cost_reordering.h         \
    reordering/lowest_inversion_reordering.h    \
    reordering/lowest_memory_reordering.h       \
    reordering/random_reordering.h              \
    reordering/reordering_base.h                \
    reordering/reordering_factory.h             \
    reordering/sink_down_reordering.h           \
    \
    storage/best.h              storage/best.cc \
    storage/bytepack.h          \
    storage/ct_none.h           \
    storage/ct_styles.h         storage/ct_styles.cc \
    storage/ct_typebased.h      \
    storage/init_storage.h      storage/init_storage.cc \
    storage/pattern.h           storage/pattern.cc \
    storage/simple.h            storage/simple.cc


## compile apps

# include_HEADERS = meddly.h

##
## List of headers to be copied into (prefix)/include/meddly
## during installation
##

pkginclude_HEADERS = \
    compute_table.h \
    ct_entry_key.h \
    ct_entry_result.h \
    ct_entry_type.h \
    ct_initializer.h \
    dd_edge.h \
    defines.h \
    domain.h \
    enumerator.h \
    error.h \
    forest.h \
    global_rebuilder.h \
    hash_stream.h \
    heap.h \
    initializer.h \
    io.h \
    meddly.h \
    memory.h \
    memstats.h \
    minterms.h \
    node_headers.h \
    node_storage.h \
    oper.h  \
    oper_binary.h  \
    oper_special.h  \
    oper_unary.h  \
    operators.h \
    opname.h  \
    opname_numer.h  \
    opname_satur.h  \
    ops_builtin.h \
    policies.h \
    relation_node.h \
    revision.h \
    unique_table.h \
    unpacked_node.h \
    variable.h \
    varorder.h

install-data-local:
	rm -f $(pkgincludedir)/*

AM_CXXFLAGS = -Wall

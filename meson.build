project(
    'Meddly',
    'cpp',
    version: run_command('head', '-n', '1', 'include/meddly/version.h', check: true).stdout().strip().split()[2].strip('"'),
    meson_version: '>= 0.49',
    default_options: ['warning_level=2', 'cpp_std=c++11']
)

if get_option('buildtype') == 'debug'
    add_project_arguments('-DDEVELOPMENT_CODE', language: 'cpp')
endif

# add_project_arguments(my_args, language: 'cpp')

deplist = [ ]

#
# Check for GMP
#

gmp_dep = dependency('gmp', required: false)
if gmp_dep.found()
    add_project_arguments('-DHAVE_LIBGMP', language: 'cpp')
else
    message('GMP library not found, disabling')
endif


#
# Check for timer support
#

mycc = meson.get_compiler('cpp')
if mycc.has_function('getrusage', prefix: '#include <sys/resource.h>')
    add_project_arguments('-DHAVE_GETRUSAGE', language: 'cpp')
endif
if mycc.has_function('gettimeofday', prefix: '#include <sys/time.h>')
    add_project_arguments('-DHAVE_SYS_TIME_H', language: 'cpp')
endif

#
# Set up source files
#

meddly_src = files(
    'src/arrays.cc',
    'src/compute_table.cc',
    'src/ct_entry_key.cc',
    'src/ct_entry_result.cc',
    'src/ct_entry_type.cc',
    'src/ct_generics.cc',
    'src/ct_initializer.cc',
    'src/ct_vector.cc',
    'src/dd_edge.cc',
    'src/domain.cc',
    'src/edge_value.cc',
    'src/error.cc',
    'src/forest.cc',
    'src/global_rebuilder.cc',
    'src/impl_unique_table.cc',
    'src/init_forests.cc',
    'src/initializer.cc',
    'src/io.cc',
    'src/io_dot.cc',
    'src/io_mdds.cc',
    'src/logger.cc',
    'src/log_json.cc',
    'src/log_simple.cc',
    'src/memory.cc',
    'src/memstats.cc',
    'src/minterms.cc',
    'src/node_headers.cc',
    'src/node_marker.cc',
    'src/node_storage.cc',
    'src/oper.cc',
    'src/oper_binary.cc',
    'src/oper_item.cc',
    'src/oper_numer.cc',
    'src/oper_satur.cc',
    'src/oper_ternary.cc',
    'src/oper_unary.cc',
    'src/operators.cc',
    'src/ops_builtin.cc',
    'src/policies.cc',
    'src/rangeval.cc',
    'src/rel_node.cc',
    'src/relation_node.cc',
    'src/sat_relations.cc',
    'src/statset.cc',
    'src/terminal.cc',
    'src/unique_table.cc',
    'src/unpacked_node.cc',
    'src/variable.cc',
    'src/varorder.cc',

    'src/forests/ev.cc',
    'src/forests/evmdd.cc',
    'src/forests/evmdd_pluslong.cc',
    'src/forests/evmdd_timesreal.cc',
    'src/forests/evmxd.cc',
    'src/forests/evmxd_pluslong.cc',
    'src/forests/evmxd_timesreal.cc',
    'src/forests/mt.cc',
    'src/forests/mtmdd.cc',
    'src/forests/mtmddbool.cc',
    'src/forests/mtmddint.cc',
    'src/forests/mtmddreal.cc',
    'src/forests/mtmxd.cc',
    'src/forests/mtmxdbool.cc',
    'src/forests/mtmxdint.cc',
    'src/forests/mtmxdreal.cc',

    'src/memory_managers/array_grid.cc',
    'src/memory_managers/freelists.cc',
    'src/memory_managers/heap_manager.cc',
    'src/memory_managers/init_managers.cc',
    'src/memory_managers/malloc_style.cc',
    'src/memory_managers/orig_grid.cc',

    'src/operations/apply_base.cc',
    'src/operations/arith_max.cc',
    'src/operations/arith_min.cc',
    'src/operations/arith_distmin.cc',
    'src/operations/arith_plus.cc',
    'src/operations/arith_minus.cc',
    'src/operations/arith_mult.cc',
    'src/operations/arith_div.cc',
    'src/operations/arith_mod.cc',
    'src/operations/cardinality.cc',
    'src/operations/compare.cc',
    'src/operations/complement.cc',
    'src/operations/constrained.cc',
    'src/operations/copy.cc',
    'src/operations/cross.cc',
    'src/operations/cycle.cc',
    'src/operations/difference.cc',
    'src/operations/dist_inc.cc',
    'src/operations/intersection.cc',
    'src/operations/maxmin_range.cc',
    'src/operations/mdd2index.cc',
    'src/operations/mm_mult.cc',
    'src/operations/prepost_rels.cc',
    'src/operations/prepost_sets.cc',
    'src/operations/prepostplus.cc',
    'src/operations/reach_dfs.cc',
    'src/operations/reach_trad.cc',
    'src/operations/satur_index.cc',
    'src/operations/satur_graph.cc',
    'src/operations/satur_queue.cc',
    'src/operations/satur_sets.cc',
    'src/operations/sat_hyb.cc',
    'src/operations/sat_impl.cc',
    'src/operations/sat_otf.cc',
    'src/operations/sat_pregen.cc',
    'src/operations/sccgraph.cc',
    'src/operations/select.cc',
    'src/operations/transitive_closure.cc',
    'src/operations/union.cc',
    'src/operations/user_unary.cc',
    'src/operations/vect_matr.cc',

    'src/storage/ct_styles.cc',
    'src/storage/init_storage.cc',
    'src/storage/simple.cc'
)

#
# Set up header files
#

meddly_inc = include_directories('include')

headers_meddly = files(
    'include/meddly/arrays.h',
    'include/meddly/compute_table.h',
    'include/meddly/ct_entry_key.h',
    'include/meddly/ct_entry_result.h',
    'include/meddly/ct_entry_type.h',
    'include/meddly/ct_initializer.h',
    'include/meddly/ct_vector.h',
    'include/meddly/dd_edge.h',
    'include/meddly/defines.h',
    'include/meddly/domain.h',
    'include/meddly/edge_value.h',
    'include/meddly/error.h',
    'include/meddly/forest.h',
    'include/meddly/forest_edgerules.h',
    'include/meddly/forest_levels.h',
    'include/meddly/global_rebuilder.h',
    'include/meddly/hash_stream.h',
    'include/meddly/heap.h',
    'include/meddly/initializer.h',
    'include/meddly/io.h',
    'include/meddly/io_dot.h',
    'include/meddly/io_mdds.h',
    'include/meddly/meddly.h',
    'include/meddly/memory.h',
    'include/meddly/memstats.h',
    'include/meddly/minterms.h',
    'include/meddly/node_headers.h',
    'include/meddly/node_marker.h',
    'include/meddly/node_storage.h',
    'include/meddly/oper.h',
    'include/meddly/oper_binary.h',
    'include/meddly/oper_item.h',
    'include/meddly/oper_numer.h',
    'include/meddly/oper_satur.h',
    'include/meddly/oper_ternary.h',
    'include/meddly/oper_unary.h',
    'include/meddly/operators.h',
    'include/meddly/ops_builtin.h',
    'include/meddly/policies.h',
    'include/meddly/rangeval.h',
    'include/meddly/rel_node.h',
    'include/meddly/relation_node.h',
    'include/meddly/sat_relations.h',
    'include/meddly/statset.h',
    'include/meddly/terminal.h',
    'include/meddly/unique_table.h',
    'include/meddly/unpacked_node.h',
    'include/meddly/variable.h',
    'include/meddly/varorder.h'
)

install_headers(headers_meddly, subdir: 'meddly')

#
# Set up library
#

meddly_lib = library('meddly',
    meddly_src,
    dependencies: gmp_dep,
    include_directories: meddly_inc,
    install: true
)

meddly_dep = declare_dependency(
    include_directories: meddly_inc,
    link_with: meddly_lib
)

#
# Build in subdirectories
#

subdir('examples')
subdir('timing')
subdir('tests')



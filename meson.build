project(
    'Meddly',
    'cpp',
    version: run_command('cat', 'include/meddly/version.h', check: true).stdout().strip().split()[2].strip('"'),
    meson_version: '>= 0.49',
    default_options: ['warning_level=2']
)

my_args = [ ]

if get_option('buildtype') == 'debug'
    my_args += [ '-DDEVELOPMENT_CODE' ]
endif

add_project_arguments(my_args, language: 'cpp')

deplist = [ ]

#
# Figure out dependencies and write to a config file
#

conf_data = configuration_data()
gmp_dep = dependency('gmp', required: false)

if gmp_dep.found()
    conf_data.set('HAVE_GMP', 1)
    deplist += [ gmp_dep ]
else
    conf_data.set('HAVE_GMP', 0)
endif

configure_file(
    output: 'config.h',
    configuration: conf_data
)

#
# Set up files
#

meddly_inc = include_directories('include')

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
# Set up library
#

meddly_lib = library('meddly',
    meddly_src,
    include_directories: meddly_inc,
    install: true
)

meddly_dep = declare_dependency(
    include_directories: meddly_inc,
    link_with: meddly_lib
)

# install_headers(headers, subdir: 'meddly')

#
# TBD: recurse into subdirectories
#



######################################################################

simple_model_lib = static_library('simple_model',
    sources: 'simple_model.cc',
    dependencies: [ gmp_dep, meddly_dep ],
    install: false
)

simple_model_dep = declare_dependency(
    dependencies: meddly_dep,
    include_directories: meddly_inc,
    link_with: simple_model_lib
)

######################################################################

zstream_lib = static_library('zstream',
    sources: 'zstream.cc',
    dependencies: meddly_dep,
    install: false
)

zstream_dep = declare_dependency(
    dependencies: meddly_dep,
    include_directories: meddly_inc,
    link_with: zstream_lib
)


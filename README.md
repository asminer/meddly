
# About MEDDLY

MEDDLY stands for Multi-terminal and Edge-valued Decision Diagram LibrarY.
It is an open-source software library that supports several types of
decision diagrams, including binary decision diagrams (BDDs).
See the
<a href="https://asminer.github.io/meddly/" target="_blank">
MEDDLY webpage
</a>
for more information.

# Building the library

You will need a C++ compiler and meson.

1.  Configure the build directory.
    Generally you need to run
    ```meson setup (builddir)```
    with the following suggestions for additional arguments.

        + `--buildtype=debug` for a debug version of the library

        + `--buildtype=release` for a release version of the library

        + `--prefix=(alternate location)` if you don't want to
            install Meddly in the usual system location.

    For developers, it is recommended to set up a development build
    directory and a release build directory. Build directory names
    should begin with `build` so they are ignored by git.
    You should only need to do this once, per build directory.


2.  Build the library.
    From within the appropriate build directory, run
    ```meson compile```


3.  Run regression tests (optional).
    From within the appropriate build directory, run
    ```meson test```


4.  Installing.
    From within the appropriate build directory, run
    ```meson install```


# Building the documentation

You will need doxygen to build the documentation.

1.  ```cd doxygen```

2.  ```./GenerateDocs.sh```


# Subdirectories

| Subdirectory  | Purpose |
|:------------- |:------- |
| `developers/` | Documentation and scripts useful to Meddly developers |
| `docs/`       | Github pages sources |
| `doxygen/`    | Generated doxygen docs from source code |
| `examples/`   | Stand-alone applications that use the library |
| `include/`    | Include files |
| `src/`        | Soource files |
| `tests/`      | Test applications used by "meson test" |
| `timing/`     | Some timing benchmarks |
| `www/`        | Sourceforge web page sources |


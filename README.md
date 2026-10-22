
# About MEDDLY

MEDDLY stands for Multi-terminal and Edge-valued Decision Diagram LibrarY.
It is an open-source software library that supports several types of
decision diagrams, including binary decision diagrams (BDDs).
See the [MEDDLY webpage](https://asminer.github.io/meddly/)
for more information.

# Building the library

You will need g++, autoconf, automake and libtool.

1.  Configure.
    Run script `developers/Config_generic.sh`
    for common situations;
    to see available switches for this script,
    run `developers/Config_generic.sh --`.

    OR, you can do this step by hand, by first running
    `./autogen.sh`
    to create the configure script, and then by running
    `./configure`.
    There are several configuration options available,
    discussed below.

2.  Build the library with `make`.
    Library will be in src/. Applications in src/ and examples/. This may be
    sufficient for most users.

3.  Run regression tests with `make check`.

4.  Install with `make install`.
    Library in lib/, include files in include/.  This can be changed
    in step (1) by selecting an alternate location, for example using
        `--prefix=/usr/local`
    will install libraries in /usr/local/lib and headers in /usr/local/include.
    If you reconfigure, make sure you run "make clean" before running "make".

Other convenient make targets:

+   `make uninstall`:
    Removes the files installed in bin/ and lib/.

+   `make clean`:
    Removes the results of compilation from src/ and examples/.

# Building the documentation

You will need doxygen to build the documentation.

1.  `cd doxygen`

2.  If you see doxygen/html/index.html and doxygen/refman.pdf, your docs have
    already been built and you can ignore the rest of the steps. Otherwise,
    run `make`.
    This will create the folders doxygen/html and doxygen/latex.
        + doxygen/html: HTML documentation. See index.html.
        + doxygen/latex: Tex-based documentation. See refman.tex.
          If you have latex installed:
            a.  `cd latex`
            b.  `make`
                This will create refman.dvi.
            c.  `make pdf`
                This will create refman.pdf (using ps2pdf).

3.  Technical documentation can be generated the same way, in doxygen-devel/

4.  Running `make clean` will remove the documentation.


# Configuration Options

Several environment variables may be set before running configure.

+   To build an optimized version of the library:
    `./configure CXXFLAGS="-O3"`
    or
    `developers/Config_generic.sh`

+   To build a debuggable version of the library:
    `./configure CXXFLAGS="-ggdb -DDEVELOPMENT_CODE"`
    or
    `developers/Config_generic.sh --debug`

+   To set the location of the gmp library, by hand
    (only necessary if ./configure fails to find it):
    `./configure CPPFLAGS="-I/path/to/gmp/header" LDFLAGS="-L/path/to/gmp/lib"`

+   To disable gmp support in Meddly:
    `./configure --without-gmp`
    or
    `developers/Config_generic.sh --without-gmp`


# Subdirectories

| Subdirectory  | Purpose |
|:------------- |:------- |
| `developers/` | Documentation and scripts useful to Meddly developers |
| `docs/`       | Online documentation for github pages |
| `doxygen/`    | Generated doxygen docs from source code |
| `examples/`   | Stand-alone applications that use the library |
| `include/`    | Include files installed here, unless prefix is changed |
| `lib/`        | Compiled library installed here, unless prefix is changed |
| `m4/`         | Directory of configuration helper scripts |
| `manual/`     | User manual source files |
| `src/`        | Library source code |
| `tests/`      | Test applications used by "make check" |
| `timing/`     | Some timing benchmarks |
| `visual/`     | Visualization utilities that work with the library. |
| `www/`        |  Web page sources |


---
title: Building
shorttitle: Building
---

## Packages required

* gcc and g++, or other C/C++ compiler with support for C++11
* meson
* doxygen, if you want to build documentation

## Packages strongly recommended

* [GMP](https://gmplib.org), for very large integer support.


## Build Steps

### 1. Configure the build directory

MEDDLY uses `meson` as its configure and build system.
For configuration, run
```bash
meson build build-dir (options)
```
with the following suggestions for additional options.
+ Use option `--buildtype=debug` for a debug version of the library
+ Use option `--buildtype=release` for a release version of the library
+ Use option `--prefix=(alternate location)` if you don't want to
            install Meddly in the usual system location.
You can run this command more than once, with different build
directories for different versions (for example, a debugging version
and a release version).

### 2.  Build

Switch to your build directory and run
```bash
meson compile
```


### 3.  Test

To run regression tests, from the build directory run
```bash
meson test
```


### 4.  Install

From within the build directory, run
```bash
meson install
```

### 5.  Build documentation (optional)

```bash
cd doxygen
./GenerateDocs.sh
```
This will generate both html and LaTeX documentation.
You can then view the html documentation by opening ```html/index.html```,
or you can run LaTeX on ```latex/refman.tex```.


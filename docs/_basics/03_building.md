---
title: Building
shorttitle: Building
---

## Packages required 

* gcc and g++, or other C/C++ compiler with support for C++11
* automake
* autoconfig
* libtool
* doxygen, if you want to build documentation

## Packages strongly recommended 

* [GMP](https://gmplib.org), for very large integer support.


## Build Steps

### 0.  Create configure script

Run the script
```bash
./autogen.sh
```
for the first build,
and whenever the source code has changes
in its file or directory layout.


### 1.  Configure

Run
```bash
./configure
```
to setup the build files.
This can be used to set various compile options;
see **Configuration Options** below.


### 2.  Run make

```bash
make
```

This builds everything.
The library will be in ```src/```,
and example applications will be in ```examples/```.


### 3.  Test

```bash
make check
```

This builds and runs several regression tests.


### 4.  Install

```bash
make install
```

This copies the compiled library into ```lib/```
and include files into ```include/```.
These locations can be changed in the ```./configure```
step by specifying an alternate location;
see **Configuration Options** below.


### 5.  Build documentation (optional)

For end-user documentation, use directory ```doxygen/```.
For complete documentation, use directory ```doxygen-devel/```.
Switch to the documentation directory and run make:
```bash
cd doxygen
make
```
This will generate both html and LaTeX documentation.
You can then view the html documentation by opening ```html/index.html```,
or you can run LaTeX on ```latex/refman.tex```.


## Other make targets

```bash
make uninstall
```
Removes the files installed in ```lib/```.


```bash
make clean
```
Removes the results of compilation.




## Configuration Options

Several environment variables may be set before running configure.
Also, configure accepts parameters to alter the library configuration.
These are listed below.

To build an optimized version of the library:
```bash
./configure CXXFLAGS="-O3"
```

To build a version of the library with support for
Extensible Decision Diagrams (XDDs):
```bash
./configure CPPFLAGS="-DUSE_XDDS"
```

To build a debuggable version of the library:
```bash
./configure CPPFLAGS="-DDEVELOPMENT_CODE" CXXFLAGS="-ggdb"
```
Note that turning on development code may produce
a significantly slower library.


To set the location of the gmp library, by hand 
(only necessary if configure fails to find it):
```bash
./configure CPPFLAGS="-I/path/to/gmp/header" LDFLAGS="-L/path/to/gmp/lib"
```

To disable gmp support in Meddly:
```bash
./configure --without-gmp
```

To change the install locations for the compiled library and include files
(for example, to ```/path/for/lib``` and ```/path/for/include```):
```bash
./configure --prefix=/path/for 
```

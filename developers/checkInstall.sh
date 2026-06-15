#!/bin/bash

#
# Check for gmp at the given locations
#   Arg1: path for include file
#   Arg2: path for compiled library
#
try_gmp()
{
    if [ -f $1/gmp.h ]; then
        for gmplib in libgmp.a libgmp.so libgmp.dll.a; do
            if [ -f $2/$gmplib ]; then
                printf "GMP Library found\n"
                printf "    gmp.h found at $1\n"
                printf "    %s found at $2\n" "$gmplib"

                gmpincl="-I$1"
                gmplink="-L$2"
                return 0
            fi
        done
    fi
    return 1
}

find_gmp()
{
    try_gmp $GMP_INCLUDE $GMP_LIBRARY && return 0
    try_gmp /opt/local/include /opt/local/lib && return 0
    try_gmp /opt/homebrew/include /opt/homebrew/lib && return 0
    try_gmp /usr/local/include /usr/local/lib && return 0
    try_gmp /usr/local/include /usr/local/lib64 && return 0
    try_gmp /usr/include /usr/lib && return 0
    try_gmp /usr/include /usr/lib64 && return 0
    try_gmp /usr/include/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu && return 0
    try_gmp /mingw/include /mingw/lib && return 0
    printf "Could not find gmp Library.\n"
    printf "You can specify the location using environment variables\n"
    printf "GMP_INCLUDE (the path to gmp.h) and GMP_LIBRARY.\n\n"
    printf "Or you can run this script with argument --without-gmp\n\n"
    exit 1
}

#
# main script
#

if [ ! -f ../configure.ac ]
then
    printf "\nRun this in the developers directory.\n\n"
    exit 1
fi
find_gmp

cd ..
rm -f include/meddly/* lib/*
if !  make install
then
    printf "\nmake install failed in parent directory.\n\n"
    cd developers
    exit 1
fi
cd developers
rm -f nqueens
printf "Rebuilding nqueens\n"
if g++ -I../include/meddly $gmpincl nqueens.cc -L../lib -lmeddly $gmplink -o nqueens
then
    printf "Build succeeded\n"
else
    printf "Build failed?\n"
    exit 1
fi

printf "Running, just to test library\n\n"
env LD_LIBRARY_PATH=../lib ./nqueens

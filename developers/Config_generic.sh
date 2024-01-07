#!/bin/bash

#
# See if the argument is of the form "--prefix=..."
#
try_prefix()
{
    local tp=`grep '^--prefix=.*' <<< "$1"`
    if [ "x" != "x$tp" ]; then
        prefix="$tp"
        return 0
    else
        return 1
    fi
}

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

                gmpargs="CPPFLAGS=-I$1 LDFLAGS=-L$2"
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

usage()
{
    printf "\n"
    printf "Usage: %s  [switches]\n\n" "$1"
    printf "Attempt to locate the GMP library, set compiler flags\n"
    printf "accordingly, and run the configure script.\n\n"
    printf "Switches:\n"
    printf "\t--debug           Turn on debugging code\n"
    printf "\t--strict          Extra compiler warnings\n"
    printf "\t--without-gmp     Turn off gmp support\n"
    printf "\t--prefix=...      Pass prefix option through to configure\n"
    printf "\n"
    exit 1
}

if [ ! -f "configure" ]; then
    echo "Run this in the root directory"
    exit 1
fi

gmpargs=""
cxxdebug="-O3"
cxxwarn="-Wall -Wno-sign-conversion -Wno-sign-compare -Wno-shadow"
prefix=""

for args; do
    if [ "--debug" == "$args" ]; then
        cxxdebug="-ggdb -DDEVELOPMENT_CODE"
        continue
    fi
    if [ "--strict" == "$args" ]; then
        cxxwarn="-Wall -Wsign-conversion -Wsign-compare -Wshadow"
        continue
    fi
    if [ "--without_gmp" == "$args" ]; then
        gmpargs="--without_gmp"
        continue
    fi
    if try_prefix "$args"; then
        continue
    fi
    usage $0
done

if [ "x" == "x$gmpargs" ]; then
    find_gmp
fi

./configure CXXFLAGS="$cxxwarn $cxxdebug" $gmpargs $prefix ||  exit 1

printf "\n"
printf "Configuration complete\n"
printf "    CXXFLAGS: %s\n" "$cxxdebug $cxxwarn"
printf "    GMP     : %s\n" "$gmpargs"
printf "    %s\n" "$prefix"
printf "\n"

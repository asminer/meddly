#!/bin/bash

# Look for gmp
# Args: directories to search
#   will echo the first matching directory.
#   If none, then error out.
#
find_gmp_path()
{
    for args; do
        echo "Looking for gmp in $args" > /dev/stderr
        if [ ! -d $args ]; then
            continue
        fi
        echo "    directory exists" > /dev/stderr
        if [ ! -d $args/include ]; then
            continue
        fi
        echo "    include/ exists" > /dev/stderr
        if [ ! -d $args/lib ]; then
            continue
        fi
        echo "    lib/ exists" > /dev/stderr
        if [ ! -f $args/include/gmp.h ]; then
            continue
        fi
        echo "    include/gmp.h exists" > /dev/stderr
        if [ ! -f $args/lib/libgmp.a ]; then
            continue
        fi
        echo "    lib/libgmp.a exists" > /dev/stderr
        echo $args
        return 0
    done
}

GMP=`find_gmp_path /opt/local /opt/homebrew`
if [ ! "$GMP" ]; then
    echo
    echo "Could not find gmp"
    echo
    exit 1
fi

./configure CXXFLAGS="-ggdb -DDEVELOPMENT_CODE" LDFLAGS=-L/$GMP/lib CPPFLAGS=-I/$GMP/include

echo
echo Configured using $GMP for gmp, with development tests and debug info on
echo

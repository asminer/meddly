#!/bin/bash
#
#

#
# Arg1: git url
# Arg2: dirname (without Meddly_version.number)
# Arg3: version
#
clone_old()
{
    declare -r url="$1"
    declare -r dirname="$2"
    declare -r version="$3"
    declare -r targdir="$dirname/Meddly_$version"

    if [ "$dry" ]; then
        echo "Clone (old) into $targdir"
        return 0
    fi
    if [ -d "$targdir" ]; then
        echo "$targdir exists, skipping version $version"
        return 0
    fi

    git clone $url -b releases/v$version $targdir
}

#
# Arg1: git url
# Arg2: dirname (without Meddly_version.number)
# Arg3: version
#
clone_new()
{
    declare -r url="$1"
    declare -r dirname="$2"
    declare -r version="$3"
    declare -r targdir="$dirname/Meddly_$version"

    if [ "$dry" ]; then
        echo "Clone (new) into $targdir"
        return 0
    fi
    if [ -d "$targdir" ]; then
        echo "$targdir exists, skipping version $version"
        return 0
    fi

    git clone $url -b v$version $targdir
}

#
# Arg1: key to search for
# Arg2..Argn: list to search
#
contains()
{
    local key="$1"
    shift
    for args; do
        if [ "$key" == "$args" ]; then
            return 0
        fi
    done
    return 1
}

usage()
{
    echo
    echo "Usage: $1 [switches] dirname"
    echo
    echo "Will clone a copy of the library for each release, each into directory"
    echo "    dirname/Meddly_version.number"
    echo
    echo "Valid switches:"
    echo
    echo "    -d"
    echo "        Dry run. Just show what versions would be copied."
    echo
    echo "    -v version"
    echo "        Can be repeated. Specify versions to copy. For example,"
    echo "        to copy versions 0.02 and 0.03 into your home directory:"
    echo "            # $1 -v 0.02 -v 0.03 ~"
    echo
    echo "    -x version"
    echo "        Can be repeated. Specify versions to exclude."
    echo
    exit 0
}

#
# Process command line arguments
#

dry=""
vlist=""
xlist=""
dirname=""

scriptname="$0"
workdir=`pwd`

while [ $# -gt 0 ]; do
    if [ "x$1" == "x-d" ]; then
        dry="y"
        shift
        continue
    fi
    if [ "x$1" == "x-v" ]; then
        vlist="$vlist $2"
        shift
        shift
        continue
    fi
    if [ "x$1" == "x-x" ]; then
        xlist="$xlist $2"
        shift
        shift
        continue
    fi
    if [ "$dirname" ]; then
        usage $scriptname
    fi
    dirname="$1"
    shift
done

if [ ! "$dirname" ]; then
    usage $scriptname
fi

#
# Get repository info:
#   url to fetch from
#   list of old version numbers (branches origin/releases/vN.NN -> N.NN)
#   list of new version numbers (tags     VN.NN.NN -> N.NN.NN)
#

declare -r gitfetch=`git remote -v | grep fetch | grep origin | awk '{print $2}'`
declare -r allold=`git branch -a | grep "origin/releases/" | awk -F/ '{print $NF}' | sed 's/^v//' | xargs `
declare -r allnew=`git tag | grep 'v[0-9.]*' | sed 's/^v//' | xargs`

#
# Determine old and new lists, based on specified versions.
# If no versions specified, use the entire list.
#

if [ "$vlist" ]; then
    oldlist=""
    newlist=""
    for v in $allold; do
        if contains $v $vlist; then
            oldlist="$oldlist $v"
        fi
    done
    for v in $allnew; do
        if contains $v $vlist; then
            newlist="$newlist $v"
        fi
    done
else
    oldlist="$allold"
    newlist="$allnew"
fi

#
# Remove any excluded versions.
#

tmplist="$oldlist"
oldlist=""
for v in $tmplist; do
    if contains $v $xlist; then
        continue
    fi
    oldlist="$oldlist $v"
done
tmplist="$newlist"
newlist=""
for v in $tmplist; do
    if contains $v $xlist; then
        continue
    fi
    newlist="$newlist $v"
done

#
# Actually clone
#

for v in $oldlist; do
    clone_old $gitfetch $dirname $v
done
for v in $newlist; do
    clone_new $gitfetch $dirname $v
done


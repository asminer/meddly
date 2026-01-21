#!/bin/bash

if [ ! -f ../configure.ac ]
then
    printf "\nRun this in the developers directory.\n\n"
    exit 1
fi

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
if g++ -I../include/meddly nqueens.cc -L../lib -lmeddly -o nqueens
then
    printf "Build succeeded\n"
else
    printf "Build failed?\n"
    exit 1
fi

printf "Running, just to test library\n\n"
env LD_LIBRARY_PATH=../lib ./nqueens

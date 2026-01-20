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
if g++ -I../include/meddly -L../lib nqueens.cc -lmeddly -o nqueens
then
    printf "Build succeeded\n"
else
    printf "Build failed?\n"
fi

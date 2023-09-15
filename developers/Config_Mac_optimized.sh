#!/bin/bash

echo Configuring with development tests off, compiler optimizations on
./configure CXXFLAGS="-O3" LDFLAGS=-L/opt/local/lib CPPFLAGS=-I/opt/local/include

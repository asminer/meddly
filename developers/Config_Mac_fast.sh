#!/bin/bash

echo Configuring with development tests, compiler optimizations
./configure CXXFLAGS="-O3 -DDEVELOPMENT_CODE" LDFLAGS=-L/opt/local/lib CPPFLAGS=-I/opt/local/include

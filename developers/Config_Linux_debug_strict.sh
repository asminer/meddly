#!/bin/bash

echo Configuring with development tests, debug info
./configure CXXFLAGS="-Wsign-conversion -Wshadow -ggdb -DDEVELOPMENT_CODE"

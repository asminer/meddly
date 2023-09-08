#!/bin/bash

echo Configuring with development tests, debug info
./configure CXXFLAGS="-ggdb -DDEVELOPMENT_CODE"

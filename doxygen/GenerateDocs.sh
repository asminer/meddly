#!/bin/bash

VER=$(awk -F\" '{print $2}' ../include/meddly/version.h | head -n 1)
sed "s/MEDDLYVERSION/$VER/" Doxyfile.skel > Doxyfile
doxygen Doxyfile

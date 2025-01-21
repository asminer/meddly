#!/bin/bash

if [ "$1" ];
then
    cmdline=`sed 's/\./ /g' <<< "$1"`
    ./$cmdline
fi

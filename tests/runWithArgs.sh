#!/bin/bash
#

dn=`dirname $0`
script="$1"
cmdline=`sed 's/\./ /g' <<< "$script"`
cmdline="$dn/$cmdline"
printf 'Running %s/%s\n' "$dn" "$cmdline"
$dn/$cmdline


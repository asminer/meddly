#!/bin/bash
#


#
# Process arguments
#

SCRIPTDIR="`pwd`/`dirname $0`"

if [ "$#" -lt 1 ]; then
    echo > /dev/stderr
    echo "Usage: $0 file.html dir dir dir ..." > /dev/stderr
    echo > /dev/stderr
    echo "Run benchmarks in each directory, and collect the results into" > /dev/stderr
    echo "an html table written to standard output." > /dev/stderr
    echo > /dev/stderr
    exit 1
fi

OUTFILE="$1"
shift
echo "Writing to $OUTFILE"

for darg; do
    $SCRIPTDIR/RunBench.sh -d $darg -h $OUTFILE
done


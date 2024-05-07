#!/bin/bash
#
# Script to run tests in directory timing/
#

#
# List of executables to run
#

RUNLIST="comptab ddedge fids mark opids reduce unpacked"
PID="$$"
RUNDATE=`date +"%Y %h %d %H:%M:%S"`

#
# Check for a valid directory
#
checkDir() {
    if [ ! -d $1 ]; then
        echo "Can't find `basename $1`/ directory; run this script in"
        echo "the root directory or within developers/."
        exit 1
    fi
}

#
# Cleanup
#
cleanup() {
    rm -f report.$PID.*.txt
}

#
# Terminate cleanly on ctrl-C
#
ctrl_c() {
    printf "\n\nCaught CTRL-C, cleaning up\n\n"
    cleanup
    exit 1
}

trap ctrl_c INT

#
# Process command-line switches
#

BPATH=""
MPATH=""
while [ $# -gt 0 ]; do
    if [ "x$1" == "x-b" ]; then
        BPATH=$2
        shift
        shift
        continue
    fi

    if [ "x$1" == "x-d" ]; then
        MPATH=$2
        shift
        shift
        continue
    fi

    printf "\nUsage: $0 [options]\n\nOptions:\n";
    printf "\t-b <path>   Save benchmark info to the given directory.\n"
    printf "\t-d <path>   Directory of version to test; otherwise . or ..\n"
    printf "\n"
    exit 1
done

#
# Figure out where things are
#

if [ "$MPATH" ]; then
    EXDIR="$MPATH/examples"
    TDIR="$MPATH/timing"
elif [ -d "examples" ]; then
    EXDIR=examples
    TDIR=timing
else
    EXDIR=../examples
    TDIR=../timing
fi

if [ "$BPATH" ]; then
    if [ ! -d $BPATH ]; then
        echo "Didn't find benchmark directory \"$BPATH\"."
        exit 1
    fi

    HTMLFILE="$BPATH/timing.html"
    counter=0
    if [ -f "$BPATH/timing.count" ]; then
        read counter < "$BPATH/timing.count"
    fi
    counter="$[ counter + 1 ]"
    TEXTFILE=`printf "%s/timing.%03d.txt" "$BPATH" "$counter"
else
    HTMLFILE=""
    TEXTFILE=""
fi

checkDir $EXDIR
checkDir $TDIR

#
# Grab revision information
#
declare -r VERSION=$($EXDIR/libinfo | awk '{print $3}')
declare -r REVDATE=$($EXDIR/libinfo 5)
if [ ! "$MPATH" ]; then
    declare -r BRANCH=$(git status | head -n 1 | awk '{print $3}')
fi


#
# Make sure executables can be found
#
for r in $RUNLIST; do
    if [ ! -x $TDIR/$r ]; then
        echo "Couldn't find executable $r,"
        echo "Are you running this in the main directory?"
        echo
        exit 1
    fi
done

#
# Run everything
#

for r in $RUNLIST; do
    echo "======================================================================"
    echo "Running $r"
    echo "======================================================================"
    $TDIR/$r -r report.$PID.$r.txt
done

#
# Tabulate everything
#

myPrint() {
    echo "$1 times"
    echo "======================================================================"
    while read time prog line; do
        details=`awk -F$ '{print $2}' <<< $line`
        printf "    %7.3f ........ %s\n" "$time" "$details"
    done < report.$PID.$1.txt
    echo
}

printSummary() {
    while read time prog line; do
        sumry=`awk -F$ '{print $1}' <<< $line`
        echo "    <th>$sumry</th>"
    done < report.$PID.$1.txt
}

textInfo() {
    printf "\n\nTiming tests for $VERSION (branch: $BRANCH):\n"
    printf "\tReleased $REVDATE\n"
    printf "\tRun on $RUNDATE\n\n"

    for r in $RUNLIST; do
        myPrint $r
    done
}

appendHTML() {
    if [ ! "$1" ]; then
        return
    fi

    if [ ! -f "$1" ]; then
        echo "Initializing $1"
        cat > $1 <<EOF
<!DOCTYPE html>
<html>
<head>
<style>
table {
    border-collapse:collapse;
}
th {
    padding-top:2px;
    padding-bottom:2px;
    padding-left:5px;
    padding-right:5px;
}
td {
    text-align:right;
    padding-top:2px;
    padding-bottom:2px;
    padding-left:5px;
    padding-right:5px;
}
</style>
</head>

<body>
<p>
Table of MEDDLY timing tests,
primarily for tracking performance changes over time.
</p>
<table border=1>
<tr>
    <th rowspan=2>Test date</th>
    <th rowspan=2>Release date</th>
    <th rowspan=2>Version</th>
EOF
        for r in $RUNLIST; do
            lines=`wc -l < report.$PID.$r.txt`
            printf "    <th colspan=%d>$r</th>\n" "$lines" >> $1
        done
        echo "</tr>" >> $1
        echo "<tr>" >> $1
        for r in $RUNLIST; do
            printSummary $r >> $1
        done
        echo "</tr>" >> $1
        echo "</table>" >> $1
        echo "</body>" >> $1
        echo "</html>" >> $1
    fi

    #
    # Strip off the bottom </table> </body> </html> lines
    #

    mv -f $1 $1.old
    grep -v "</table>" $1.old | grep -v "</body>" | grep -v "</html>" > $1

    #
    # Add new lines
    #
    echo "<tr>" >> $1
    echo "    <td>$RUNDATE</td>" >> $1
    echo "    <td>$REVDATE</td>" >> $1
    echo "    <td>$VERSION $BRANCH</td>" >> $1
    for r in $RUNLIST; do
        while read time other; do
            printf "    <td> %3.3f </td>\n" "$time" >> $1
        done < report.$PID.$r.txt
    done
    echo "</tr>" >> $1

    echo "</table>" >> $1
    echo "</body>" >> $1
    echo "</html>" >> $1

    rm -f $1.old
}

textInfo | tee $TEXTFILE

appendHTML $HTMLFILE
if [ "$counter" ]; then
    echo "$counter" > $BPATH/timing.count
fi

cleanup


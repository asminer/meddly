#!/bin/bash
#

declare -r TIMEOUT=300
declare -r TIMECHK=2

declare -r OUTTXT="out.$$.txt"
declare -r TIMETXT="time.$$.txt"

#
# Arrays for stats.
#   (1) Name of the stat
#   (2) command to build (#s for source dir, #e for examples dir)
#   (3) printf formatting
#   (4) Result of the command
#
statName=(
    'Version'
    'Released'
    '# Files'
    '# Lines'
    '# Chars'
)
statCmd=(
    "#e/libinfo | awk '{print \$3}'"
    "#e/libinfo 5"
    "ls #s/*.h #s/*.hh #s/*.cc #s/*/*.h #s/*/*.cc | wc -l"
    "wc -l #s/*.h #s/*.hh #s/*.cc #s/*/*.h #s/*/*.cc | tail -n 1 | awk '{print \$1}'"
    "wc -c #s/*.h #s/*.hh #s/*.cc #s/*/*.h #s/*/*.cc | tail -n 1 | awk '{print \$1}'"
)
statPrint=(
    "%s"
    "%s"
    "%'d"
    "%'d"
    "%'d"
)
statResult=("")

#
# Arrays for benchmarks.
#   (1) Name
#   (2) Command (#e for example dir)
#   (3) timing result
#   (4) peak memory result
#
benchName=(
    "nqueens 14"
    "qcover 12"
    #
    "phils 800 (trad)"
    "kanban 75 (trad)"
    "slot 20 (trad)"
    #
    "phils 10000 (sat)"
    "kanban 200 (sat)"
    "slot 100 (sat)"
    #
    "kanban 6 (exp)"
    "slot 7 (exp)"
)
benchCmd=(
    "#e/nqueens 14"
    "#e/queen_cover 12"
    #
    "#e/dining_phils -n800"
    "#e/kanban 75 -bfs"
    "#e/slot 20 -bfs"
    #
    "#e/dining_phils -n10000 -dfs"
    "#e/kanban 200 -dfs"
    "#e/slot 100 -dfs"
    #
    "#e/kanban 6 -exp"
    "#e/slot 7 -exp"
)
benchTime=("")
benchPeak=("")

#
# Show all results
#
showResults()
{
    local i=""

    if [ "$branch" ]; then
        printf "\nSource code stats ($branch branch):\n"
    else
        printf "\nSource code stats:\n"
    fi
    for i in ${!statName[@]}; do
        printf "%16s: ${statPrint[$i]}\n" "${statName[$i]}" "${statResult[$i]}"
    done

    printf "\nBenchmark results:\n"
    local totstring="0"
    for i in ${!benchName[@]}; do
        if [ "${benchTime[$i]}" ]; then
            printf "%20s:  %7.2f sec  " "${benchName[$i]}" "${benchTime[$i]}"
            totstring="$totstring + ${benchTime[$i]}"
            if [ "TO" == "${benchPeak[$i]}" ]; then
                printf "(timeout)\n"
            else
                printf "peak: %s\n" "${benchPeak[$i]}"
            fi
        fi
    done
    printf "\n%20s:  %7.2f sec\n\n" "Total" `bc <<< "$totstring"`
}

#
# Append results to html table
#
showHtml()
{
    printf "<tr>\n"
    for i in ${!statName[@]}; do
        printf "    <td>${statPrint[$i]}</td>\n" "${statResult[$i]}"
    done
    for i in ${!benchName[@]}; do
        if [ "TO" == "${benchPeak[$i]}" ]; then
            printf "    <td>timeout</td>\n"
        else
            printf "    <td>%7.2f sec</td>\n" "${benchTime[$i]}"
        fi
    done

    printf "</tr>\n"
}


#
# Run with a timeout
#
# Arg1: command to run (possibly with #e)
#
runWithTimeout()
{
    local cmd=`sed "s|#e|$EXDIR|g" <<< "$1"`
    local barecmd=`sed "s|#e/||g" <<< "$1" | awk '{print $1}'`

    echo
    echo "======================================================================"
    echo "Running $cmd"
    echo "======================================================================"
    echo

    time ( $cmd 2>&1 | tee $OUTTXT ) 2> $TIMETXT &
    PID=$!

    SECONDS=0
    while true; do
        sleep $TIMECHK
        if ps $PID > /dev/null; then
            #
            # still running
            if [ $SECONDS -ge $TIMEOUT ]; then
                echo
                echo "  Timeout exceeded ($TIMEOUT seconds), terminating"
                echo
                local kidpid=`ps | grep "$barecmd" | awk '{print $1}'`
                kill -9 $kidpid

                return 1
            fi
        else
            # completed

            return 0
        fi
    done
}


#
# Clean up temp files
#
cleanup()
{
    touch $OUTTXT $TIMETXT
    rm $OUTTXT $TIMETXT
}

#
# Catch ctrl-c
#
bailout()
{
    echo
    echo "Caught CTRL-C, terminating script"
    echo
    cleanup
    exit 1
}

trap bailout INT

#
# Process switches (these are new!)
#

MPATH=""
HTMLFILE=""
TEXTFILE=""

while [ $# -gt 0 ]; do
    if [ "x$1" == "x-d" ]; then
        MPATH=$2
        shift
        shift
        continue
    fi

    if [ "x$1" == "x-h" ]; then
        HTMLFILE=$2
        shift
        shift
        continue
    fi

    if [ "x$1" == "x-t" ]; then
        TEXTFILE=$2
        shift
        shift
        continue
    fi

    printf "\nUsage: $0 [options]\n\nOptions:\n";
    printf "\t-d <path>   Directory of version to test; otherwise . or ..\n"
    printf "\t-h <file>   Append output, in html format, to file.\n"
    printf "\t-t <file>   Write  output, in text format, to file.\n"
    printf "\n"
    exit 1
done

#
# Find examples and source directories
#

if [ "$MPATH" ]; then
    EXDIR="$MPATH/examples"
    SRCDIR="$MPATH/src"
elif [ -d "examples" ]; then
    EXDIR="examples"
    SRCDIR="src"
else
    EXDIR="../examples"
    SRCDIR="../src"
fi

if [ ! -d $EXDIR ]; then
    echo "Didn't find examples directory (tried \"$EXDIR\")."
    echo "Run this script in the root directory or within developers,"
    echo "or use the -d switch."
    exit 1
fi
if [ ! -d $SRCDIR ]; then
    echo "Didn't find source directory (tried \"$EXDIR\")."
    echo "Run this script in the root directory or within developers,"
    echo "or use the -d switch."
    exit 1
fi

if [ ! -f $EXDIR/kanban ]; then
  echo "Didn't find executable $EXDIR/kanban"
  exit 1
fi


if [ "$TEXTFILE" ]; then
    echo Writing text summary to file: $TEXTFILE
fi
if [ "$HTMLFILE" ]; then
    echo Writing html summary to file: $HTMLFILE
fi

#
# Determine stats
#

for i in ${!statName[@]}; do
    cmd=`sed "s|#e|$EXDIR|g" <<< "${statCmd[$i]}" | sed "s|#s|$SRCDIR|g"`
    statResult[$i]=`echo "$cmd" | bash 2> /dev/null`
done
if [ ! "$MPATH" ]; then
  branch=$(git status | head -n 1 | awk '{print $3}')
fi

#
# Run benchmarks
#

for i in ${!benchName[@]}; do
    showResults

    cleanup
    if runWithTimeout "${benchCmd[$i]}"; then
        benchTime[$i]=`awk '/user/{print $2}' $TIMETXT | tr 'ms' '  ' | awk '{print $1*60+$2}'`
        benchPeak[$i]=`awk '/peak memory allocated/{print $1" "$2";"}' $OUTTXT | xargs`
    else
        benchTime[$i]=$[ $TIMEOUT ]
        benchPeak[$i]="TO"
    fi
done
cleanup

#
# Show results
#

showResults | tee $TEXTFILE

if [ "$HTMLFILE" ]; then
    showHtml >> $HTMLFILE
fi


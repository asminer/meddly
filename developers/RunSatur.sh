#!/bin/bash

METHOD1="msat"
METHOD2="sat1"

first_word()
{
    echo $1
}

#
# Args: directories to try
#
find_kanban()
{
    for arg; do
        if [ -x "$arg/kanban" ]
        then
            echo "$arg"
            return
        fi
    done
}

#
# Arg1: model
# Arg2: trad time
# Arg3: front time
# Arg4: trad final
# Arg5: front final
# Arg6: trad peak
# Arg7: front peak
table_row()
{
    printf "%-11s || %7.3f | %7.3f || %6d | %6d || %8d | %8d\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7"
}

table_header()
{
    printf -- "            ||  Saturation Time  ||   Final Nodes   ||     Peak   Nodes   \n"
    printf -- " Model      || %7s | %7s || %6s | %6s || %8s | %8s\n" "$METHOD1" "$METHOD2" "$METHOD1" "$METHOD2" "$METHOD1" "$METHOD2"
    printf -- "------------++---------+---------++--------+--------++----------+---------\n"
}

show_table()
{
    echo

    table_header

    for m in 'phils 10000' 'kanban 200' 'slot 100' 'tiles 3x4'
    do
        base=`first_word $m`

        file1=/tmp/$base.$METHOD1.txt
        file2=/tmp/$base.$METHOD2.txt

        # echo "$m  /tmp/$base.trad.txt"

        if [ "tiles" == "$base" ]
        then
            time1=`awk '/reachable states construction took/{print $5}' $file1`
            time2=`awk '/reachable states construction took/{print $5}' $file2`
        else
            time1=`awk '/Reachability set construction took/{print $5}' $file1`
            time2=`awk '/Reachability set construction took/{print $5}' $file2`
        fi

        if [ "phils" == "$base" ]
        then
            final1=`awk '/#Nodes:/{print $2}' $file1`
            final2=`awk '/#Nodes:/{print $2}' $file2`
        fi
        if [ "kanban" == "$base" ]
        then
            final1=`awk '/Reachability set uses/{print $4}' $file1`
            final2=`awk '/Reachability set uses/{print $4}' $file2`
        fi
        if [ "slot" == "$base" ]
        then
            final1=`awk '/Reachability set uses/{print $4}' $file1`
            final2=`awk '/Reachability set uses/{print $4}' $file2`
        fi
        if [ "tiles" == "$base" ]
        then
            final1=`awk '/MDD for reachable states/{print $6}' $file1`
            final2=`awk '/MDD for reachable states/{print $6}' $file2`
        fi

        peak1=`awk '/MDD stats:/{mdd=1} /peak nodes/{if (mdd) print $1}' $file1`
        peak2=`awk '/MDD stats:/{mdd=1} /peak nodes/{if (mdd) print $1}' $file2`


        # echo "    trad $ttime sec $tfinal final nodes $tpeak peak"
        # echo "    front $ftime sec $ffinal final nodes $fpeak peak"

        table_row "$m" "$time1" "$time2" "$final1" "$final2" "$peak1" "$peak2"
    done
    echo
}

path=`find_kanban examples build-release/examples ../build-release/examples`
if [ ! "$path" ]
then
    echo
    echo "Couldn't find examples directory with executables."
    echo "Run this script in the build directory, or elsewhere"
    echo "if the build directory is named build-release/"
    echo
    exit 1
fi

for prog in dining_phils kanban slot tiles
do
    if [ ! -x $path/$prog ]
    then
        echo "Couldn't find $prog executable"
        exit 1
    fi
done

$path/dining_phils -$METHOD1 -n10000 | tee /tmp/phils.$METHOD1.txt
$path/dining_phils -$METHOD2 -n10000 | tee /tmp/phils.$METHOD2.txt

$path/kanban -$METHOD1 200 | tee /tmp/kanban.$METHOD1.txt
$path/kanban -$METHOD2 200 | tee /tmp/kanban.$METHOD2.txt

$path/slot -$METHOD1 100 | tee /tmp/slot.$METHOD1.txt
$path/slot -$METHOD2 100 | tee /tmp/slot.$METHOD2.txt

$path/tiles --pos --$METHOD1 3 4 | tee /tmp/tiles.$METHOD1.txt
$path/tiles --pos --$METHOD2 3 4 | tee /tmp/tiles.$METHOD2.txt

show_table | tee .satur.txt

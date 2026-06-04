#!/bin/bash

first_word()
{
    echo $1
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
    printf "%-10s || %10s | %10s || %6s | %6s || %10s | %10s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7"
}

table_header()
{
    printf -- "           || RSS Generation Time (s) ||   Final Nodes   ||       Peak   Nodes     \n"
    printf -- " Model     || Trad.      | Frontier   || Trad,  | Front. || Trad.      | Frontier  \n"
    printf -- "-----------++------------+------------++--------+--------++------------+-----------\n"
}

path="examples"
if [ ! -d examples ]
then
    path="../examples"
    if [ ! -d ../examples ]
    then
        echo
        echo "Couldn't find examples directory."
        echo "Run this script in the root or developers directory."
        echo
        exit 1
    fi
fi

for prog in dining_phils kanban slot tiles
do
    if [ ! -x $path/$prog ]
    then
        echo "Couldn't find $prog executable"
        exit 1
    fi
done

$path/dining_phils -trad -n800 | tee /tmp/phils.trad.txt
$path/dining_phils -front -n800 | tee /tmp/phils.front.txt

$path/kanban -trad 75 | tee /tmp/kanban.trad.txt
$path/kanban -front 75 | tee /tmp/kanban.front.txt

$path/slot -trad 20 | tee /tmp/slot.trad.txt
$path/slot -front 20 | tee /tmp/slot.front.txt

$path/tiles --pos --trad 3 3 | tee /tmp/tiles.trad.txt
$path/tiles --pos --front 3 3 | tee /tmp/tiles.front.txt


echo

table_header

for m in 'phils 800' 'kanban 75' 'slot 20' 'tiles 3x3'
do
    base=`first_word $m`

    tfile=/tmp/$base.trad.txt
    ffile=/tmp/$base.front.txt

    # echo "$m  /tmp/$base.trad.txt"

    if [ "tiles" == "$base" ]
    then
        ttime=`awk '/reachable states construction took/{print $5}' $tfile`
        ftime=`awk '/reachable states construction took/{print $5}' $ffile`
    else
        ttime=`awk '/Reachability set construction took/{print $5}' $tfile`
        ftime=`awk '/Reachability set construction took/{print $5}' $ffile`
    fi

    if [ "phils" == "$base" ]
    then
        tfinal=`awk '/#Nodes:/{print $2}' $tfile`
        ffinal=`awk '/#Nodes:/{print $2}' $ffile`
    fi
    if [ "kanban" == "$base" ]
    then
        tfinal=`awk '/Reachability set uses/{print $4}' $tfile`
        ffinal=`awk '/Reachability set uses/{print $4}' $ffile`
    fi
    if [ "slot" == "$base" ]
    then
        tfinal=`awk '/Reachability set uses/{print $4}' $tfile`
        ffinal=`awk '/Reachability set uses/{print $4}' $ffile`
    fi
    if [ "tiles" == "$base" ]
    then
        tfinal=`awk '/MDD for reachable states/{print $6}' $tfile`
        ffinal=`awk '/MDD for reachable states/{print $6}' $ffile`
    fi

    tpeak=`awk '/MDD stats:/{mdd=1} /peak nodes/{if (mdd) print $1}' $tfile`
    fpeak=`awk '/MDD stats:/{mdd=1} /peak nodes/{if (mdd) print $1}' $ffile`


    # echo "    trad $ttime sec $tfinal final nodes $tpeak peak"
    # echo "    front $ftime sec $ffinal final nodes $fpeak peak"

    table_row "$m" "$ttime" "$ftime" "$tfinal" "$ffinal" "$tpeak" "$fpeak"
done
echo

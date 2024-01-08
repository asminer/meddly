#!/bin/bash
#

write_header()
{
    cat <<EOF
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
Gigantic table of MEDDLY stats,
primarily for comparing library versions.
Tested on
EOF
    date +"%d %B %Y."
    cat <<EOF
</p>
<table border=1>
<tr>
    <th colspan=5>Source Information</th>
    <th colspan=2>Constraint benchmarks</th>
    <th colspan=3>Traditional reachability</th>
    <th colspan=3>Saturation reachability</th>
    <th colspan=2>Explicit reachability</th>
</tr>
<tr>
    <th>Version</th>
    <th>Release Date</th>
    <th>files</th>
    <th>lines</th>
    <th>chars</th>
    <th>nqueens 14</th>
    <th>qcover 12</th>
    <th>phils 800</th>
    <th>kanban 75</th>
    <th>slot 20</th>
    <th>phils 10k</th>
    <th>kanban 200</th>
    <th>slot 100</th>
    <th>kanban 6</th>
    <th>slot 7</th>
</tr>
<!-- Data here -->
EOF
}

write_footer()
{
    cat <<EOF
<!-- Done data -->
</table>
</body>
</html>
EOF
}

#
# Process arguments
#

SCRIPTDIR=`pwd`

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

write_header > $OUTFILE
for darg; do
    $SCRIPTDIR/RunBench.sh -d $darg -h $OUTFILE
done
write_footer >> $OUTFILE


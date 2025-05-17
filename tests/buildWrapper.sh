#!/bin/bash
#

for args; do
    script="$args"
    noext=`basename -s ._sh "$script"`
    if [ "$script" != "$noext" ]
    then
        cmdline=`sed 's/\./ /g' <<< "$noext"`
        echo '#!/bin/bash' > $script
        echo 'echo "Running $cmdline $@"' >> $script
        echo 'DIR=`dirname $0`' >> $script
        echo '$DIR/$cmdline $@' >> $script
        chmod +x $script
    fi
done

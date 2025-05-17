#!/bin/bash
#

for args; do
    script="$args"
    noext=`basename -s ._sh "$script"`
    if [ "$script" != "$noext" ]
    then
        cmdline=`sed 's/\./ /g' <<< "$noext"`
        echo '#!/bin/bash' > $script
        echo 'DIR=`dirname $0`' >> $script
        printf 'echo "Running $DIR/%s $@"\n' "$cmdline" >> $script
        printf '$DIR/%s $@\n' "$cmdline" >> $script
        chmod +x $script
    fi
done

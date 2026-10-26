#!/bin/bash
#
#

DRYRUN=""   # is this a dry run.
RELTYPE=""  # release type, m or M.

DEFFILE="include/meddly/defines.h"
REVFILE="include/meddly/version.h"
RELDIR="docs/_releases"
UNRELNOTES="$RELDIR/unreleased.md"
NAVFILE="docs/_data/navigation.yml"

OLDVER=""   # version number currently
THISVER=""  # version number of this release
NEXTVER=""  # next version number after release

#
# Display usage and exit
#
usage()
{
    printf "\nUsage: $1 [-d] -M | -m \n\n"
    printf "Creates a release.\n"
    printf "Updates the version number and release date in sources.\n"
    printf "Creates a tag for the release, named v[version].\n\n"
    printf "Before running, be sure that all modified files have been committed\n"
    printf "to the repository.\n\n"
    printf "Switches:\n"
    printf "    -d: dry run. Just show what would be done.\n\n"

    printf "    -m: minor release. If configure.ac specifies version x.y.z, then\n"
    printf "        the current release will be x.y.z, and configure.ac will be\n"
    printf "        changed to version x.y.z+1.\n\n"

    printf "    -M: major release. If configure.ac specifies version x.y.z, then\n"
    printf "        the current release will be x.y+1.0, and configure.ac will be\n"
    printf "        changed to version x.y+1.1.\n\n"
    exit 1
}

#
# Determine current and next version numbers
#
determine_versions()
{
    OLDVER=`awk -F\" '{print $2}' ../$REVFILE | head -n 1`

    if [ "M" == "$RELTYPE" ]
    then
        THISVER=`awk -F. '{print $1"."($2+1)".0"}' <<< "$OLDVER"`
    else
        THISVER="$OLDVER"
    fi
    NEXTVER=`awk -F. '{print $1"."$2"."($3+1)}' <<< "$THISVER"`

    echo "Current version: $OLDVER"
    echo "Release version: $THISVER"
    echo "Next    version: $NEXTVER"
}


#
# Arg1: new version
# Arg2: empty string, or "unreleased"
#
update_revfile()
{
    local tab
    local outfile
    local rdate

    if [ "$DRYRUN" ]
    then
        tab="    "
        outfile="/dev/stdout"
        echo "Updated $REVFILE would be:"
    else
        tab=""
        outfile="../$REVFILE"
        echo "Updating $REVFILE"
    fi

    if [ "$2" ]
    then
        rdate="($2)"
    else
        rdate=`date +"%Y %B %d"`
    fi

    printf "%s#define VERSION \"%s\"\n%sconst char* MEDDLY_DATE = \"%s\"\n"  "$tab" "$1"  "$tab" "$rdate" > $outfile

    if [ "$DRYRUN" ]
    then
        echo "would now do: git add ../$REVFILE"
    else
        git add ../$REVFILE
    fi
}

#
# Update version in ALLOW_DEPRECATED
#
# Arg1: original version
# Arg2: updated version
# Arg3: file to update
#
update_deprecated()
{
    if [ ! -f $3 ]
    then
        return
    fi

    local srch="ALLOW_DEPRECATED_$1"
    local repl="ALLOW_DEPRECATED_$2"

    sed "s/$srch/$repl/" $3 > $3.new

    if diff -q $3 $3.new > /dev/null
    then
        if [ "$DRYRUN" ]
        then
            rm $3.new
        fi
    else
        if [ "$DRYRUN" ]
        then
            echo "Updating $3:"
            diff $3 $3.new
            rm $3.new
        else
            mv -f $3.new $3
            git add $3
        fi
    fi
}


#
# No arguments
#
update_deffile()
{
    if [ ! -f "../$DEFFILE" ]; then
        echo "Couldn't find $DEFFILE, could not update!"
        return
    fi
    local ddate=`date +"%m/%d/%Y"`
    local oldversid=`tr '.' '_' <<< "$OLDVER"`
    local thisversid=`tr '.' '_' <<< "$THISVER"`

    sed "1,/=======/s|$oldversid.*|$thisversid  \"$ddate\"|" ../$DEFFILE > ../$DEFFILE.new

    if [ "$DRYRUN" ]
    then
        echo "The followung updates would be made to $DEFFILE:"
        diff ../$DEFFILE ../$DEFFILE.new
        rm ../$DEFFILE.new
    else
        mv -f ../$DEFFILE.new ../$DEFFILE
        git add ../$DEFFILE
    fi

    if [ "$oldversid" != "$thisversid" ]
    then
        find ../include ../src | while read fname
        do
            update_deprecated $oldversid $thisversid $fname
        done
    fi

}



#
# Parse command line arguments
#

for arg
do
    if [ "-m" == "$arg" ]
    then
        RELTYPE="m"
        continue
    fi
    if [ "-M" == "$arg" ]
    then
        RELTYPE="M"
        continue
    fi
    if [ "-d" == "$arg" ]
    then
        DRYRUN="d"
        continue
    fi
    echo "Unknown argument or switch: $arg"
    usage $0
done

#
# Check command line arguments
#

if [ ! $RELTYPE ]
then
    echo "Missing required -m or -M."
    usage $0
fi
if [ "m" == "$RELTYPE" ]
then
    echo "Minor release"
else
    echo "Major release"
fi
if [ "$DRYRUN" ]
then
    echo "DRY RUN - no changes will be made"
fi
echo

#
# Sanity checks
#

if [ ! -f Release.sh ]
then
    echo "Run this in the developers directory."
    exit 1
fi
if [ ! -f ../$UNRELNOTES ]
then
    echo "Missing release docs; check file $UNRELNOTES"
    exit 1
fi
if diff -q blank_unreleased.md ../$UNRELNOTES > /dev/null
then
    echo "Update release documentation $UNRELNOTES"
    exit 1
fi
remotes=`git remote`
for r in origin sourceforge
do
    if ! grep $r > /dev/null <<< $remotes
    then
        echo "Missing expected remote: $r"
        exit 1
    fi
done

#
# Start "making changes"
#

determine_versions

if [ -f ../$RELDIR/$THISVER.md ]; then
    printf "Can't update documentation; file $RELDIR/$THISVER.md exists\n"
    exit 1
fi

update_revfile $THISVER
update_deffile

#
# TBD HERE
#

exit 1


#
#
confirm_or_quit()
{
    while [ 0 ]
    do
        printf "Is this OK? Type yes to continue, no to quit\n"
        read resp
        if [ "$resp" == "no" ]
        then
            exit 1
        fi
        if [ "$resp" == "yes" ]
        then
            return
        fi
    done
}

#
# Arg1: version
#
update_defines()
{
    printf "    updating $DEFFILE\n"
    if [ ! -f "../$DEFFILE" ]; then
        printf "\t$DEFFILE not found\n";
        return
    fi
    local ddate=`date +"%m/%d/%Y"`
    local versid=`tr '.' '_' <<< "$1"`
    mv -f ../$DEFFILE ../$DEFFILE.old
    sed "1,/=======/s|$versid.*|$versid  \"$ddate\"|" ../$DEFFILE.old > ../$DEFFILE
    git add ../$DEFFILE
}

#
# Arg1: version
#
update_docs()
{
    local ddate=`date +"%Y-%m-%d"`
    echo "---"
    echo "title: New in Version $1"
    echo "shorttitle: Version $1"
    echo "number: $1"
    echo "date: $ddate"
    echo "---"

    sed '/^---$/,/^---$/d' ../$UNRELNOTES
}

#
# Arg1: version
#
update_nav()
{
    # Copy first part
    #
    sed -n '1,/url: "\/releases\/unreleased"/p' ../$NAVFILE

    #
    # Add new entry
    echo "      - title: \"Version $1\""
    echo "        url: \"/releases/$1\""

    # Copy remaining part
    #
    sed '1,/url: "\/releases\/unreleased"/d' ../$NAVFILE
}

#
# Arg1: version
#
update_config()
{
    sed "s|^\(AC_INIT(MEDDLY\), [0-9.]*, \(.*\)|\1, $1, \2|" ../configure.ac
}


#
# Build version numbers
#

version=`awk -F, '/AC_INIT/{print $2}' ../configure.ac | tr -d '[],'`
oldversion="$version"
if [ "$major" ]; then
    version=`awk -F. '{for (i=1; i<NF-1; i++) printf($i"."); print $(NF-1)+1".0"}' <<< $version`
fi
nextver=`awk -F. '{for (i=1; i<NF; i++) printf($i"."); print $NF+1}' <<< $version`
printf "Creating release for version $version\n"
printf "Next version will be $nextver\n\n"

confirm_or_quit

printf "Checking release...\n"

#
# Release checks and whatnot
#

printf "Preparing relase\n"

cp ../configure.ac ../configure.ac.old
if [ "$version" != "$oldversion" ]; then
    printf "    updating configure.ac\n"
    update_config $version > ../configure.ac.new
    mv -f ../configure.ac.new ../configure.ac
fi

update_revision "$version"
update_defines "$version"

printf "    creating $RELDIR/$version.md\n"
update_docs "$version"  > ../$RELDIR/$version.md
cp -f blank_unreleased.md ../$UNRELNOTES

printf "    updating $NAVFILE\n"
update_nav $version > newnav.yml
mv -f newnav.yml ../$NAVFILE


#
# Save changes
#

git add ../$RELDIR/$version.md ../$UNRELNOTES ../$NAVFILE ../configure.ac
git commit -m "Updated for $version release"
git push
git push sourceforge

#
# Create and publish a tag for $version
#

git tag -a "v$version" -m "Tag for release $version"
git push origin v$version
git push sourceforge v$version


#
# Build distribution
#

printf "Building tarball for upload to sourceforge\n"
make -C .. dist

#
# Automatically bump the version number
#

printf "Updating version to $nextver\n"
update_config "$nextver" > ../configure.ac.new
mv -f ../configure.ac.new ../configure.ac

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
# Indent some text
#
indent()
{
    sed 's/^/    |    /'
}

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
    echo
}


#
# Arg1: new version
# Arg2: date prefix
#
update_revfile()
{
    local rdate=`date +"%Y %B %d"`

    if [ "$2" ]
    then
        rdate="$2 $rdate"
    fi

    printf "#define VERSION \"%s\"\nconst char* MEDDLY_DATE = \"%s\"\n"  "$1" "$rdate" > ../$REVFILE.new

    if [ "$DRYRUN" ]
    then
        echo "Updated $REVFILE:"
        indent < ../$REVFILE.new
        rm ../$REVFILE.new
    else
        mv -f ../$REVFILE.new ../$REVFILE
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
            diff $3 $3.new | indent
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
        diff ../$DEFFILE ../$DEFFILE.new | indent
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
# No arguments
#
addto_deffile()
{
    if [ ! -f "../$DEFFILE" ]; then
        echo "Couldn't find $DEFFILE, could not update!"
        return
    fi
    local nextvid=`tr '.' '_' <<< "$NEXTVER"`

    sed "1,/=======/s|==>|#define ALLOW_DEPRECATED_$nextvid\n// ==>|" ../$DEFFILE > ../$DEFFILE.new

    if [ "$DRYRUN" ]
    then
        echo "We would add this to $DEFFILE:"
        diff ../$DEFFILE ../$DEFFILE.new | indent
        rm ../$DEFFILE.new
    else
        mv -f ../$DEFFILE.new ../$DEFFILE
        git add ../$DEFFILE
    fi
}

#
# No arguments
#
update_docs()
{
    local ddate=`date +"%Y-%m-%d"`
    echo "---"
    echo "title: New in Version $THISVER"
    echo "shorttitle: Version $THISVER"
    echo "number: $THISVER"
    echo "date: $ddate"
    echo "---"

    sed '/^---$/,/^---$/d' ../$UNRELNOTES
}

#
# No arguments
#
update_nav()
{
    # Copy first part
    #
    sed -n '1,/url: "\/releases\/unreleased"/p' ../$NAVFILE

    #
    # Add new entry
    echo "      - title: \"Version $THISVER\""
    echo "        url: \"/releases/$THISVER\""

    # Copy remaining part
    #
    sed '1,/url: "\/releases\/unreleased"/d' ../$NAVFILE
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
# update web docs
#
newdoc="$RELDIR/$THISVER.md"
if [ "$DRYRUN" ]
then
    echo "Created $newdoc would be:"
    update_docs | indent
else
    echo "Creating $newdoc"
    update_docs > ../$newdoc
    cp -f blank_unreleased.md ../$UNRELNOTES
    git add ../$newdoc ../$UNRELNOTES
fi

#
# update web navigation
#
update_nav > newnav.yml
if [ "$DRYRUN" ]
then
    echo "Changes to $NAVFILE:"
    diff newnav.yml ../$NAVFILE | indent
    rm -f newnav.yml
else
    echo "Updating $NAVFILE"
    mv -f newnav.yml ../$NAVFILE
    git add ../$NAVFILE
fi

#
# Commit git changes, tag new release
#
if [ "$DRYRUN" ]
then
    echo "(git commit)"
    echo "(git push)"
    echo "(git tag -a v$THISVER)"
    echo "(git push)"
else
    echo "Committing changes"
    git commit -m "Updated for $THISVER release"
    git push
    git push sourceforge

    echo "Tagging release"
    git tag -a "v$THISVER" -m "Tag for release $THISVER"
    git push origin v$THISVER
    git push sourceforge v$THISVER
fi

#
# Build the tarball for this release
#

tarbase="meddly-$THISVER"
tarfile="$tarbase.tar.gz"
echo "Building $tarfile"
tar czf ../$tarfile -C .. --exclude .OLD -s "|^|$tarbase/|" doxygen/ examples/ include/ src/ tests/ timing/ AUTHORS COPYING COPYING.LESSER meson.build README.md

#
# Automatically bump the version number
#

update_revfile $NEXTVER after
addto_deffile


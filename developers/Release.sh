#!/bin/bash
#
#

DEFFILE="src/defines.h"
REVFILE="src/revision.h"
RELDIR="docs/_releases"
UNRELNOTES="$RELDIR/unreleased.md"
NAVFILE="docs/_data/navigation.yml"

#
# Arg1: version
#
update_revision()
{
    printf "    updating $REVFILE\n"
    local rdate=`date +"%Y %B %d"`
    printf "const char* MEDDLY_DATE = \"%s\";\nconst char* MEDDLY_VERS = \"%s\";\n" "$rdate" "$1" > ../$REVFILE
    git add ../$REVFILE
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

usage()
{
    printf "\nUsage: $1 -M | -m \n\n"
    printf "Creates a release and update the version number in configure.ac.\n"
    printf "Automatically updates the revision date in the source,\n"
    printf "and creates a tag for the release, named v[version].\n\n"
    printf "Before running, be sure that all modified files have been committed\n"
    printf "to the repository.\n\n"
    printf "Switches:\n"
    printf "    -m: minor release. If configure.ac specifies version x.y.z, then\n"
    printf "        the current release will be x.y.z, and configure.ac will be\n"
    printf "        changed to version x.y.z+1.\n\n"
    printf "    -M: major release. If configure.ac specifies version x.y.z, then\n"
    printf "        the current release will be x.y+1.0, and configure.ac will be\n"
    printf "        changed to version x.y+1.1.\n\n"
    exit 1
}

#
# Checks
#

if [ $# -ne 1 ]; then
    usage $0
fi

major="x"
if [ "x-m" == "x$1" ]; then
    major=""
fi
if [ "x-M" == "x$1" ]; then
    major="y"
fi
if [ "x" == "$major" ]; then
    usage $0
fi

if [ ! -f ../configure.ac ]; then
    printf "\nRun this in the developers directory.\n\n"
    exit 1
fi

if [ ! -f ../$UNRELNOTES ]; then
    printf "Missing release docs; check file $UNRELNOTES\n\n"
    exit 1
fi

if diff -q blank_unreleased.md ../$UNRELNOTES; then
    printf "Update release documentation $UNRELNOTES\n\n"
    exit 1
fi

remotes=`git remote`
for r in origin sourceforge; do
    if ! grep $r > /dev/null <<< $remotes; then
        printf "Missing expected remote: $r\n\n"
        exit 1
    fi
done

#
# Determine version number
#

version=`awk '/AC_INIT/{print $2}' ../configure.ac | tr -d '[],'`
oldversion="$version"
if [ "$major" ]; then
    version=`awk -F. '{for (i=1; i<NF-1; i++) printf($i"."); print $(NF-1)+1".0"}' <<< $version`
fi

if [ -f ../$RELDIR/$version.md ]; then
    printf "Can't update documentation; file $RELDIR/$version.md exists\n"
    exit 1
fi


printf "Creating release for version $version\n"

cp ../configure.ac ../configure.ac.old
if [ "$version" != "$oldversion" ]; then
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

nextver=`awk -F. '{for (i=1; i<NF; i++) printf($i"."); print $NF+1}' <<< $version`

printf "\nNext version will be $nextver\n"

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

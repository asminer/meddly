#!/bin/bash
#
#

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

version=`awk '/AC_INIT/{print $2}' ../configure.ac | tr -d '[],'`
if [ "$major" ]; then
    version=`awk -F. '{for (i=1; i<NF-1; i++) printf($i"."); print $(NF-1)+1".0"}' <<< $version`
fi


printf "Creating release for version $version\n"

rdate=`date +"%Y %B %d"`
printf "const char* MEDDLY_DATE = \"%s\";\n" "$rdate" > ../src/revision.h
printf "const char* MEDDLY_VERS = \"%s\";\n" "$version" >> ../src/revision.h

nextver=`awk -F. '{for (i=1; i<NF; i++) printf($i"."); print $NF+1}' <<< $version`

echo "New version would be $nextver"

#
# Save changes
#

git add ../src/revision.h
git commit -m "Updated release date for $version"
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

make dist

#
# Automatically bump the version number
#

printf "Updating version to $nextver\n"
sed "/AC_INIT/s/$version/$nextver/" ../configure.ac > ../configure.new
mv -f ../configure.ac ../configure.ac.old
mv ../configure.new ../configure.ac


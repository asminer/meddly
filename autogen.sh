# $Id$

# aclocal -I config
aclocal
autoheader
libtoolize --force
automake --foreign --add-missing --copy --force-missing
autoconf
# ./configure --prefix=$PWD

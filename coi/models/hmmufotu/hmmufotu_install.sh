#!/usr/bin/env bash
set -e -x

cd $TMPDIR

module load boost/1.79.0

EIGENVERSION=3.3.5
rm -rf eigen-$EIGENVERSION
wget -O-\
    "https://gitlab.com/libeigen/eigen/-/archive/$EIGENVERSION/eigen-$EIGENVERSION.tar.gz" |
tar -xzf-

VERSION=1.5.1
rm -rf HmmUFOtu-$VERSION
wget -O- https://github.com/Grice-Lab/HmmUFOtu/archive/refs/tags/v$VERSION.tar.gz |
  tar -xzf -

cd HmmUFOtu-$VERSION
./configure --prefix=/projappl/project_2005718/\
            --with-boost-libdir=$BOOST_INSTALL_ROOT/lib\
            CPPFLAGS="-I$TMPDIR/eigen-$EIGENVERSION"\
            CFLAGS="-O3 -march=native -mtune=native"\
            CXXFLAGS="-O3 -march=native -mtune=native"

make
make install

/projappl/project_2005718/bin/hmmufotu -h

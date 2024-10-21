#!/usr/bin/env bash
set -x

VERSION=4.0.4

cd $TMPDIR/$USER
rm -rf veryfasttree-$VERSION
wget -O-\
 https://github.com/citiususc/veryfasttree/archive/refs/tags/v${VERSION}.tar.gz |
tar -xzf -

cd veryfasttree-4.0.4
cmake -DCMAKE_INSTALL_PREFIX:PATH=/projappl/project_2005718\
      -DUSE_NATIVE=ON\
      -DUSE_AVX512=ON\
      .
make
make install

/projappl/project_2005718/bin/VeryFastTree -h

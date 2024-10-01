#!/usr/bin/env bash

# need to download eigen first

wget -O- https://github.com/Grice-Lab/HmmUFOtu/archive/refs/tags/v1.5.1.tar.gz |
  tar -xzf -

cd HmmUFOtu-1.5.1
./configure --prefix=/projappl/project_2005718/
export CPPFLAGS="$CPPFLAGS -I/usr/include/eigen3"
export CFLAGS="$CFLAGS -O3 -march=native -mtune=native"
export CXXFLAGS="$CFLAGS -O3 -march=native -mtune=native"

make
make install
cd ..

#!/usr/bin/env bash
export CFLAGS="-O3 -march=native -mtune=native"
export CXXFLAGS="-O3 -march=native -mtune=native"

wget -O- https://github.com/pierrebarbera/epa-ng/archive/refs/tags/v0.3.8.tar.gz |
tar -xzf -

cd epa-ng-0.3.8
make

cd ..

wget -O- https://github.com/lczech/gappa/archive/refs/tags/v0.8.5.tar.gz |
tar -xzf -

cd gappa-0.8.5
make

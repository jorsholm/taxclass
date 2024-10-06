#!/usr/bin/env bash
set -e -x

cd $TMPDIR
INSTALL_PREFIX="/projappl/project_2005718"
export CFLAGS="-O3 -march=native -mtune=native"
export CXXFLAGS="-O3 -march=native -mtune=native"
PATH="$INSTALL_PREFIX/bin:$PATH"

EPANGVERSION=0.3.8
rm -rf epa-ng-$EPANGVERSION

wget -O- https://github.com/pierrebarbera/epa-ng/archive/refs/tags/v0.3.8.tar.gz |
tar -xzf -

cd epa-ng-0.3.8
make
cp bin/epa-ng $INSTALL_PREFIX/bin
cd ..

epa-ng -v

GAPPAVERSION=0.8.5
wget -O- https://github.com/lczech/gappa/archive/refs/tags/v$GAPPAVERSION.tar.gz |
tar -xzf -

cd gappa-$GAPPAVERSION
make
cp bin/gappa $INSTALL_PREFIX/bin
cd ..

gappa -v

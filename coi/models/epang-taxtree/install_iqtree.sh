#!/usr/bin/env bash
set -e -x

# load module for compiler and boost libraries
module load gcc/13.2.0 boost/1.85.0

cd $TMPDIR
INSTALL_PREFIX="/projappl/project_2005718"
export CFLAGS="-O3 -march=native -mtune=native"
export CXXFLAGS="-O3 -march=native -mtune=native"
PATH="$INSTALL_PREFIX/bin:$PATH"

# download Eigen3 headers
EIGENVERSION=3.4.0
wget -O- https://gitlab.com/libeigen/eigen/-/archive/${EIGENVERSION}/eigen-${EIGENVERSION}.tar.gz |
tar -xzf -

# download IQTREE2 from github
# (source tarball does not include required submodules)
IQTREEVERSION=2.4.0
git clone --recursive https://github.com/iqtree/iqtree2.git
cd iqtree2

# set specified release version
git checkout v${IQTREEVERSION}

# compile
mkdir build
cd build
cmake -DEIGEN3_INCLUDE_DIR=../../eigen-3.4.0/Eigen  ..
make

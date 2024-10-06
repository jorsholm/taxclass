#!/usr/bin/env bash
set -x

VERSION=4.0.4

module load gcc/13.2.0 cuda/12.6.0
cd $TMPDIR/$USER

rm -rf veryfasttree-$VERSION
wget -O-\
 https://github.com/citiususc/veryfasttree/archive/refs/tags/v${VERSION}.tar.gz |
tar -xzf -

cd veryfasttree-$VERSION
cmake -DCMAKE_CUDA_ARCHITECTURES=70\
      -DCMAKE_EXECUTABLE_SUFFIX=_cuda\
      -DCMAKE_INSTALL_PREFIX:PATH=/projappl/project_2005718\
      -DUSE_NATIVE=ON\
      -DUSE_AVX512=ON\
      -DUSE_CUDA=ON\
      .
make
make install

/projappl/project_2005718/bin/VeryFastTree_cuda -h

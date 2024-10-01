#!/usr/bin/env bash

wget https://github.com/citiususc/veryfasttree/archive/refs/tags/v4.0.4.tar.gz
tar -xzf v4.0.4.tar.gz

cd veryfasttree-4.0.4
cmake -DUSE_NATIVE=ON -DUSE_AVX512=ON -DUSE_CUDA=ON
make

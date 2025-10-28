#!/usr/bin/env bash

# Launch all CPU tests for epang-taxtree
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=epang-taxtree_$CPUS.out \
         --error=epang-taxtree_$CPUS.out \
         run_epang_taxtree.sh
done


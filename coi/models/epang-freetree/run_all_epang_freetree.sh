#!/usr/bin/env bash

# Launch all CPU tests for epang-freetree
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=epang-freetree_$CPUS.out \
         --error=epang-freetree_$CPUS.out \
         run_epang_freetree.sh
done


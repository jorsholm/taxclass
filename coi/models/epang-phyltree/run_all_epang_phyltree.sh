#!/usr/bin/env bash

# Launch all CPU tests for epang-phyltree
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=epang-phyltree_$CPUS.out \
         --error=epang-phyltree_$CPUS.out \
         run_epang_phyltree.sh
done


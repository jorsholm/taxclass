#!/usr/bin/env bash

# Launch all CPU tests for protax
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=protax_$CPUS.out \
         --error=protax_$CPUS.out \
         run_protax.sh
done


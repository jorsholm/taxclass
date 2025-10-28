#!/usr/bin/env bash

# Launch all CPU tests for protax-a
for CPUS in 1
do
  sbatch --cpus-per-task=$CPUS \
         --output=protax-a_$CPUS.out \
         --error=protax-a_$CPUS.out \
         run_protax.sh
done


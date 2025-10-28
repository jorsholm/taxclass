#!/usr/bin/env bash

# Launch all CPU tests for sintax
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=sintax_$CPUS.out \
         --error=sintax_$CPUS.out \
         run_sintax.sh
done


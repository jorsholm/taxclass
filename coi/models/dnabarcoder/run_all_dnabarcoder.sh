#!/usr/bin/env bash

# Launch all CPU tests for dnabarcoder
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=dnabarcoder_$CPUS.out \
         --error=dnabarcoder_$CPUS.out \
         run_dnabarcoder.sh
done


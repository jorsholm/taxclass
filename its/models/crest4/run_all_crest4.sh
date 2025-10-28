#!/usr/bin/env bash

# Launch all CPU tests for crest4
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=crest4_$CPUS.out \
         --error=crest4_$CPUS.out \
         run_crest4.sh
done


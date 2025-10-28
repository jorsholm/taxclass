#!/usr/bin/env bash

# Launch all CPU tests for blast
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=blast_$CPUS.out \
         --error=blast_$CPUS.out \
         run_blast.sh
done


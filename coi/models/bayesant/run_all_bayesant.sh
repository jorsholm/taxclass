#!/usr/bin/env bash

# Launch all CPU tests for bayesant
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS\
         --output=bayesant_$CPUS.out\
         --error=bayesand_$CPUS.out\
         run_bayesant.sh
done


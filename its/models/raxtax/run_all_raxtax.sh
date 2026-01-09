#!/usr/bin/env bash

# Launch all CPU tests for sintax
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=raxtax_$CPUS.out \
         --error=raxtax_$CPUS.out \
         run_raxtax.sh
done


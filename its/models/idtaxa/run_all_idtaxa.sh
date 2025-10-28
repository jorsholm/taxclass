#!/usr/bin/env bash

# Launch all CPU tests for idtaxa
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=idtaxa_$CPUS.out \
         --error=idtaxa_$CPUS.out \
         run_idtaxa.sh
done


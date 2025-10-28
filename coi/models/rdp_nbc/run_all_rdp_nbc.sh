#!/usr/bin/env bash

# Launch all CPU tests for rdp_nbc
for CPUS in 1
do
  sbatch --cpus-per-task=$CPUS \
         --output=rdp_nbc_$CPUS.out \
         --error=rdp_nbc_$CPUS.out \
         run_rdp_nbc.sh
done


#!/usr/bin/env bash

# Launch all CPU tests for mycoai_cnn
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=mycoai_cnn_$CPUS.out \
         --error=mycoai_cnn_$CPUS.out \
         run_mycoai_cnn.sh
done

# Launch GPU test for mycoai_cnn
sbatch run_mycoai_cnn_gpu.sh


#!/usr/bin/env bash

# Launch all CPU tests for mycoai_bert
for CPUS in 1 4 16 40
do
  sbatch --cpus-per-task=$CPUS \
         --output=mycoai_bert_$CPUS.out \
         --error=mycoai_bert_$CPUS.out \
         run_mycoai_bert.sh
done

# Launch GPU test for mycoai_bert
sbatch run_mycoai_bert_gpu.sh

#!/usr/bin/env bash
#SBATCH --job-name=bayesant
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --array=16
#SBATCH --output=bayesant_%a.out
#SBATCH --error=bayesant_%a.out
#SBATCH --mail-type=ALL

# installation commands (in R)
# install.packages("renv")
# renv::init(bare=TRUE)
# renv::install("alessandrozito/BayesANT@38d0ac4")
# renv::snapshot()

# put project bin directory (including R) on the path
export PATH="/projappl/project_2005718/bin:$PATH"

# set environmental variables to tell various parallel computation libraries
# how many CPU cores we have
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# GNU time is not on the path by default on compute nodes
export PATH="/appl/opt/time/1.9/bin:$PATH"
TIME="$(which time) --verbose"

# define file names
MODEL=bayesant
DATA=../../data
RESULTS=../../results/$MODEL

# make sure the results directory exists
mkdir -p $RESULTS

for alphabet in nt
do
  # train the model
  export TRAIN_FILE=$DATA/train_${alphabet}_label.fasta
  export MODEL_FILE=train_${alphabet}_$SLURM_ARRAY_TASK_ID.rds

  $TIME Rscript fit_BayesANT.R

  # classify the test sequences
  for TEST in test testshort
  do
    export TEST_FILE=$DATA/${TEST}_${alphabet}.fasta
    export RESULT_FILE=${RESULTS}/${MODEL}_${TEST}_${alphabet}_$SLURM_ARRAY_TASK_ID.tsv

    $TIME Rscript predict_BayesANT.R
  done
done

echo Job $SLURM_JOB_ID completed

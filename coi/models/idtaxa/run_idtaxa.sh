#!/usr/bin/env bash
#SBATCH --job-name=idtaxa
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1
#SBATCH --output=idtaxa_%a.out
#SBATCH --error=idtaxa_%a.out
#SBATCH --mail-type=ALL

# installation commands (in R)
# install.packages("renv")
# renv::init()
# renv::install("bioc::DECIPHER")

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
MODEL=idtaxa
DATA=../../data
RESULTS=../../results/$MODEL

# make sure the results directory exists
mkdir -p $RESULTS

for alphabet in nt aa
do
  # train the model
  TRAIN_FILE=$DATA/train_$alphabet.fasta
  TAX_FILE=$DATA/train_tax.tsv
  MODEL_FILE=train_$alphabet_$SLURM_ARRAY_TASK_ID.rds

  $TIME Rscript learn_taxa.R

  # classify the test sequences
  for TEST in test testshort
  do
    TEST_FILE=$DATA/${TEST}_nt.fasta
    RESULT_FILE=${RESULTS}/${MODEL}_${TEST}_${alphabet}_$SLURM_ARRAY_TASK_ID.tsv

    $TIME Rscript id_taxa.R
  done
done

echo Job $SLURM_JOB_ID completed
